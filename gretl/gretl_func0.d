/* Port of the shared head of gretl_func.c (11,869 lines in the C): user-defined functions and function packages.

   The file is ported in five modules, as geneval.c was: this one holds the shared head (the private structs, flags, the file-level state
   and the macros of the first 290 lines of the preprocessed source), gretl_func1.d .. gretl_func4.d hold the functions in their original order
   (cut points at function boundaries). There is no `static` across modules in D, so the functions of the slices are all public (names are
   unique since they were one C file) and the slices import what they need from each other by name, cyclically. gretl.gretl_func is a facade
   that publicly imports the four slices, so that importers of "gretl_func" (the C header gretl_func.h, gretl_funch.d for the types) find
   everything.

   Naming: the typedefs fn_param, fn_arg, obsinfo, linegen, ufunc, fncall and fnpkg are the structs FnParam, FnArg, ObsInfo, LineGen, Ufunc,
   Fncall and Fnpkg; fn_line is Stmt (cmd_privateh.d); the anonymous enums and enum ln_flags / FCFlags are FpFlags, LnFlags and FCFlags with lower case
   members (FP_CONST -> FpFlags.const_, LINE_IGNORE -> LnFlags.ignore, FC_RECURSING -> FCFlags.recursing). The function-like macros are functions.
   The file-level statics are public __gshared variables here.

   Preprocessing (as in the other ports): the debug switches (CALL_DEBUG, EXEC_DEBUG, ARGS_DEBUG, PKG_DEBUG, FN_DEBUG, COMP_DEBUG, REC_DEBUG,
   FNPARSE_DEBUG, VALDEBUG) are off and their output is not ported; COMPILE_RECURSIVE is 0; ONE_WORKSPACE is 1; HAVE_MPI and USE_RLIB (the
   embedded R library check) are not defined, so those branches are not in the port. */

module gretl.gretl_func0;

import gretl.cmd_privateh : Stmt;
import gretl.genmainh : Generator;
import gretl.glib : guint8;
import gretl.glist : GList;
import gretl.gretl_bundleh;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_funch;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.uservarh : UserVar;

enum INT_USE_XLIST = -999;
enum INT_USE_MYLIST = -777;

enum FpFlags {
    const_    = 1 << 0, /* explicitly marked as "const" */
    optional  = 1 << 1, /* marked as optional (null is OK) */
    auto_     = 1 << 2  /* marked as automatic (has a default value) */
}

/* structure representing a parameter of a user-defined function */

struct FnParam {
    char * name;     /* the name of the parameter */
    GretlType type;  /* its type (C: char) */
    char * descrip;  /* its description */
    char ** labels;  /* value labels, if applicable */
    int nlabels;     /* number of value labels */
    char flags;      /* additional information (e.g. "const" flag) */
    double deflt;    /* default value */
    double min;      /* minimum value (scalar parameters only) */
    double max;      /* maximum value (scalar parameters only) */
    double step;     /* step increment (scalars only) */
}

bool param_is_const(const(FnParam) * p) { return (p.flags & FpFlags.const_) != 0; }
bool param_is_optional(const(FnParam) * p) { return (p.flags & FpFlags.optional) != 0; }
bool param_is_auto(const(FnParam) * p) { return (p.flags & FpFlags.auto_) != 0; }

void set_param_optional(FnParam * p) { p.flags |= FpFlags.optional; }
void set_param_auto(FnParam * p) { p.flags |= (FpFlags.optional | FpFlags.auto_); }

enum LnFlags {
    ignore = 1 << 0,
    nocomp = 1 << 1  /* "don't compile genr" */
}

bool line_has_loop(const(Stmt) * l) { return l.ci == GretlCmdIndex.loop && l.ptr != null; }
bool ignore_line(const(Stmt) * l) { return (l.flags & LnFlags.ignore) != 0; }
bool line_no_comp(const(Stmt) * l) { return (l.flags & LnFlags.nocomp) != 0; }

enum UNSET_VALUE = -1.0e200;

bool default_unset(const(FnParam) * p) { return p.deflt == UNSET_VALUE; }

/* structure representing sample information at start of a function call */

struct ObsInfo {
    int structure;       /* time-series, etc. */
    int pd;              /* data frequency */
    int t1, t2;          /* starting and ending observations */
    int added;           /* number of observations added within function */
    char[OBSLEN] stobs_; /* string representation of starting obs (see Ufunc: property stobs) */
    @property char * stobs() { return stobs_.ptr; }
    @property const(char) * stobs() const { return stobs_.ptr; }
    int panel_pd;        /* panel time frequency, if applicable */
    double panel_sd0;    /* panel time starting point */
}

/* structure representing a call to a user-defined function */

enum FCFlags {
    recursing = 1 << 0, /* call to f() has a call to f() in its ancestry */
    prev_msgs = 1 << 1, /* record of "messages" setting before call */
    prev_echo = 1 << 2, /* record of "echo" setting before call */
    preserve  = 1 << 3  /* function call should be preserved */
}

/* Note: the FC_PRESERVE flag is set when we want to allow for repeated
   calls to a given function. In the non-recursive case we just reuse
   the stored fncall struct. In the recursive case we need to ensure
   there's a distinct call struct for each depth of function execution,
   and the function then acquires a GList of fncalls.
*/

struct LineGen {
    int idx;
    Generator * genr;
}

struct Fncall {
    Ufunc * fun;     /* the function called */
    int argc;        /* argument count */
    int orig_v;      /* number of series defined on entry */
    FnArg * args;    /* argument array */
    int * ptrvars;   /* list of pointer arguments */
    int * listvars;  /* list of series included in a list argument */
    GList * lists;   /* list of names of list arguments */
    char * retname;  /* name of return value (or dummy string) */
    GretlType rtype; /* return type (when not fixed in advance) */
    ObsInfo obs;     /* sample info */
    int flags;       /* indicators for recursive call, etc (FCFlags) */
    Stmt * line;     /* currently executing line */
    LineGen * lgen;  /* array of per-line genrs */
    int n_lgen;      /* number of linegens */
}

/* structure representing a user-defined function */

struct Ufunc {
    char[FN_NAMELEN] name_; /* identifier (C name: name; the property name gives the pointer the C code uses) */
    @property char * name() { return name_.ptr; }
    @property const(char) * name() const { return name_.ptr; }
    Fnpkg * pkg;           /* pointer to parent package, or NULL */
    Fncall * call;         /* pointer to current call, or NULL */
    GList * calls;         /* for use in recursion */
    int pkg_role;          /* printer, plotter, etc. */
    int flags;             /* private, plugin, etc. (UfunAttrs) */
    int line_idx;          /* current line index (compiling) */
    int n_lines;           /* number of lines of code */
    Stmt * lines;          /* array of lines of code */
    int n_params;          /* number of parameters */
    FnParam * params;      /* parameter info array */
    GretlType rettype;     /* return type (if any) */
    char * docstr;         /* short (optional) documentation string */
}

/* structure representing a function package */

struct Fnpkg {
    char[FN_NAMELEN] name_; /* package name (C name: name; see Ufunc) */
    @property char * name() { return name_.ptr; }
    @property const(char) * name() const { return name_.ptr; }
    char * fname;          /* filename */
    char * author;         /* author's name */
    char * email;          /* author's email address */
    char * version_;       /* package version string (C name: version) */
    char * date;           /* last revision date */
    char * descrip;        /* package description */
    char * help;           /* package help text */
    char * gui_help;       /* GUI-specific help (optional) */
    char * Rdeps;          /* R dependencies (if any) */
    char * sample;         /* sample caller script */
    char * help_fname;     /* filename: package help text */
    char * gui_help_fname; /* filename: GUI-specific help text */
    char * sample_fname;   /* filename: sample caller script */
    char * tags;           /* tag string(s) */
    char * label;          /* for use in GUI menus */
    char * mpath;          /* menu path in GUI */
    int minver;            /* minimum required gretl version */
    guint8 uses_subdir;    /* lives in subdirectory (0/1) */
    guint8 prechecked;     /* already checked for data requirement */
    guint8 data_access;    /* wants access to full data range */
    DataReq dreq;          /* data requirement */
    int modelreq;          /* required model type, if applicable */
    Ufunc ** pub;          /* pointers to public interfaces */
    Ufunc ** priv;         /* pointers to private functions */
    int n_pub;             /* number of public functions */
    int n_priv;            /* number of private functions */
    guint8 overrides;      /* number of overrides of built-in functions */
    char ** datafiles;     /* names of packaged data files */
    char ** depends;       /* names of dependencies */
    char * provider;       /* name of "provider" package, if applicable */
    char * conflict;       /* name of comflicting package, if any */
    int n_files;           /* number of data files */
    int n_depends;         /* number of dependencies */
    void * editor;         /* for GUI use */
}

/* acceptable types for parameters of user-defined functions */

bool ok_function_arg_type(GretlType t)
{
    return t == GretlType.bool_ ||
           t == GretlType.int_ ||
           t == GretlType.obs ||
           t == GretlType.double_ ||
           t == GretlType.series ||
           t == GretlType.list ||
           t == GretlType.matrix ||
           t == GretlType.string_ ||
           t == GretlType.bundle ||
           t == GretlType.scalar_ref ||
           t == GretlType.series_ref ||
           t == GretlType.matrix_ref ||
           t == GretlType.bundle_ref ||
           t == GretlType.string_ref ||
           t == GretlType.strings ||
           t == GretlType.matrices ||
           t == GretlType.bundles ||
           t == GretlType.lists ||
           t == GretlType.arrays ||
           t == GretlType.strings_ref ||
           t == GretlType.matrices_ref ||
           t == GretlType.bundles_ref ||
           t == GretlType.lists_ref ||
           t == GretlType.arrays_ref ||
           t == GretlType.numeric;
}

/* structure representing an argument to a user-defined function */

struct FnArg {
    GretlType type;       /* argument type (C: char) */
    char shifted;         /* level was shifted for execution */
    char * upname;        /* name of supplied arg at caller level */
    UserVar * uvar;       /* reference to "parent", if any */
    union Val {
        int idnum;        /* named series arg (series ID) */
        double x;         /* scalar arg */
        double * px;      /* anonymous series arg */
        GretlMatrixRaw * m; /* matrix arg */
        char * str;       /* string arg */
        int * list;       /* list arg */
        GretlBundle * b;  /* anonymous bundle pointer */
        GretlArray * a;   /* array argument */
    }
    Val val;
}

__gshared int n_ufuns;          /* number of user-defined functions in memory */
__gshared Ufunc ** ufuns;       /* array of pointers to user-defined functions */
__gshared Ufunc * current_fdef; /* pointer to function currently being defined */
__gshared GList * callstack;    /* stack of function calls */
__gshared int n_pkgs;           /* number of loaded function packages */
__gshared Fnpkg ** pkgs;        /* array of pointers to loaded packages */
__gshared Fnpkg * current_pkg;  /* pointer to package currently being edited */

/* record of state, and communication of state with outside world */

__gshared int compiling;        /* boolean: are we compiling a function currently? */
__gshared int fn_executing;     /* depth of function call stack */
__gshared int compiling_python;

bool function_is_private(const(Ufunc) * f) { return (f.flags & UfunAttrs.private_) != 0; }
bool function_is_noprint(const(Ufunc) * f) { return (f.flags & UfunAttrs.noprint) != 0; }
bool function_is_menu_only(const(Ufunc) * f) { return (f.flags & UfunAttrs.menu_only) != 0; }
bool function_must_assign(const(Ufunc) * f) { return (f.flags & UfunAttrs.assign) != 0; }
bool function_is_recursive(const(Ufunc) * f) { return (f.flags & UfunAttrs.recurses) != 0; }

void set_call_recursing(Fncall * c) { c.flags |= FCFlags.recursing; }
bool is_recursing(const(Fncall) * c) { return (c.flags & FCFlags.recursing) != 0; }

struct FlagAndKey {
    int flag;
    const(char) * key;
}

__gshared FlagAndKey[12] pkg_lookups = [
    FlagAndKey(UfunRole.bundle_print,  BUNDLE_PRINT.ptr),
    FlagAndKey(UfunRole.bundle_plot,   BUNDLE_PLOT.ptr),
    FlagAndKey(UfunRole.bundle_test,   BUNDLE_TEST.ptr),
    FlagAndKey(UfunRole.bundle_fcast,  BUNDLE_FCAST.ptr),
    FlagAndKey(UfunRole.bundle_extra,  BUNDLE_EXTRA.ptr),
    FlagAndKey(UfunRole.gui_main,      GUI_MAIN.ptr),
    FlagAndKey(UfunRole.gui_precheck,  GUI_PRECHECK.ptr),
    FlagAndKey(UfunRole.plot_precheck, PLOT_PRECHECK.ptr),
    FlagAndKey(UfunRole.list_maker,    LIST_MAKER.ptr),
    FlagAndKey(UfunRole.r_setup,       R_SETUP.ptr),
    FlagAndKey(UfunRole.ui_maker,      UI_MAKER.ptr),
    FlagAndKey(-1,                     null)
];

bool pkg_aux_role(int r)
{
    return r == UfunRole.bundle_print ||
           r == UfunRole.bundle_plot ||
           r == UfunRole.bundle_test ||
           r == UfunRole.bundle_fcast ||
           r == UfunRole.bundle_extra ||
           r == UfunRole.r_setup ||
           r == UfunRole.ui_maker;
}
