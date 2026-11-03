/* Port of gretl_func.c, part 1 of 4: see gretl_func0.d. */

module gretl.gretl_func1;

import gretl.genlex : delete_function_override;

import gretl.glist : g_list_free_full;

import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gretl_func_ext : get_command_index, gretl_cmd_free, gretl_cmd_init, gretl_loop_destroy;
import gretl.xmlh : xmlNodeListGetString, xmlStrcmp;

import gretl.dataseth : VarFlags;
import gretl.gretl_pathsh : PkgType;
import gretl.gretl_prnh : PrnType;

import gretl.cmd_privateh : Cmd;

import gretl.cmd_privateh : ExecState, Stmt;
import gretl.garray : GArray;
import gretl.genmainh : Generator;
import gretl.gretl_string_table : SeriesTable;
import gretl.objstackh : GretlObjType;
import gretl.xmlh : xmlDocPtr, xmlNodePtr;

import gretl.cmd_privateh : ExecState;
import gretl.dataset : series_ensure_level_zero, series_set_flag, series_unset_flag;
import gretl.genmain : destroy_genr;
import gretl.glib : gchar, gint8;
import gretl.glist : GList, g_list_append, g_list_first, g_list_free, g_list_last, g_list_remove;
import gretl.gmem : g_free;
import gretl.gretl_array : gretl_array_get_type;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_func2 : function_package_free;
import gretl.gretl_func3 : get_function_package_by_filename, get_function_package_by_name, ufunc_get_structure;
import gretl.gretl_func4 : adjust_indent, normalize_hansl;
import gretl.gretl_paths : gretl_function_package_get_path;
import gretl.gretl_prn : gretl_print_destroy, gretl_print_get_buffer, gretl_print_new, gretl_print_steal_buffer, pprintf, pputc, pputs;
import gretl.gretl_string_table : gretl_insert_builtin_string;
import gretl.gretl_typemap : gretl_is_array_type, gretl_type_from_string, gretl_type_get_name, gretl_type_get_ref_type;
import gretl.gretl_xml : gretl_xml_child_get_string, gretl_xml_child_get_strings_array, gretl_xml_get_prop_as_bool, gretl_xml_get_prop_as_double, gretl_xml_get_prop_as_int, gretl_xml_get_prop_as_string, gretl_xml_get_string, gretl_xml_put_string, gretl_xml_put_strings_array_quoted, gretl_xml_put_tagged_string;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.interacth : CMD_COMMENT;
import gretl.libset : gretl_debugging_on;
import gretl.missingh : NADBL, na;
import gretl.options : destroy_option_params_at_level;
import gretl.printout : bufgets, bufgets_finalize, bufgets_init, bufspace;
import gretl.strutils : gretl_strdup, string_is_blank, strings_array_free, strings_array_new, tailstrip;
import gretl.uservar : set_previous_depth, switch_uservar_hash;

import gretl.cmd_privateh : Stmt;
import gretl.glib : gboolean;
import gretl.uservarh : UserVar;
import gretl.xmlh : xmlDocPtr, xmlNodePtr;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc, qsort, atoi, atof;
import core.stdc.stdio : fprintf, printf, sprintf, snprintf, stderr, stdout, FILE, fopen, fclose, fgets, fputs;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr, strstr, strncmp, strrchr;

import core.stdc.stdarg : va_arg, va_end, va_list, va_start;
import core.stdc.ctype : isalpha, isspace, isdigit;
import core.stdc.errno : errno;
import core.stdc.stdio : fputc, rewind, sscanf, remove;
import core.stdc.stdlib : strtod, strtol;
import core.stdc.string : memmove, strcspn, strncat, strspn;
version (Posix) import core.sys.posix.sys.stat : stat_t; else import gretl.winport : stat_t;
import gretl.genlex : delete_function_override;
import gretl.gretl_func0;
import gretl.xmlh : xmlChar;
import gretl.gretl_funch;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;

private alias XUC = const(xmlChar) *;

int pkg_key_get_role(const(char) * key)
{
    int i;

    if (key != null && *key != '\0') {
        for (i=0; pkg_lookups[i].flag > 0; i++) {
            if (!strcmp(key, pkg_lookups[i].key)) {
                return pkg_lookups[i].flag;
            }
        }
    }

    return UfunRole.none;
}

extern (C) const(char) * package_role_get_key(int flag)
{
    int i;

    for (i=0; pkg_lookups[i].flag > 0; i++) {
        if (flag == pkg_lookups[i].flag) {
            return pkg_lookups[i].key;
        }
    }

    return null;
}

void set_function_private(Ufunc * u, gboolean s)
{
    if (s) {
        u.flags |= UfunAttrs.private_;
    } else {
        u.flags &= ~UfunAttrs.private_;
    }
}

extern (C) int gretl_compiling_function()
{
    return compiling;
}

extern (C) Ufunc * get_docstring_target()
{
    return current_fdef;
}

extern (C) int gretl_compiling_python(const(char) * line)
{
    if (compiling_python) {
        char[4] s1;
        char[8] s2;

        if (sscanf(line, "%3s %7s", s1.ptr, s2.ptr) == 2 &&
            !strcmp(s1.ptr, "end") && !strcmp(s2.ptr, "foreign")) {
            compiling_python = 0;
        }
    }

    return compiling_python;
}

void set_compiling_on()
{
    compiling = 1;
}

void set_compiling_off()
{
    compiling = compiling_python = 0;
}

extern (C) int gretl_function_depth()
{
    return fn_executing;
}

void adjust_array_arg_type(FnArg * arg)
{
    GretlType t = gretl_array_get_type(arg.val.a);

    if (arg.type == GretlType.array_ref) {
        arg.type = gretl_type_get_ref_type(t);
    } else {
        arg.type = t;
    }
}

void nullify_upname(FnArg * arg)
{
    if (arg != null && arg.upname != null) {
        free(arg.upname);
        arg.upname = null;
    }
}

int fn_arg_set_data(FnArg * arg, const(char) * name, UserVar * uvar, GretlType type, void * p)
{
    int err = 0;

    arg.type = type;
    arg.shifted = 0;
    arg.uvar = uvar;
    arg.upname = (name != null)? gretl_strdup(name) : null;

    if (type == GretlType.none) {
        arg.val.x = 0;
    } else if (type == GretlType.double_ ||
               type == GretlType.scalar_ref) {
        arg.val.x = *cast(double *) p;
    } else if (type == GretlType.int_ ||
               type == GretlType.obs) {
        arg.val.x = *cast(int *) p;
    } else if (type == GretlType.series) {
        arg.val.px = cast(double *) p;
    } else if (type == GretlType.matrix ||
               type == GretlType.matrix_ref) {
        arg.val.m = cast(GretlMatrixRaw *) p;
    } else if (type == GretlType.string_ ||
               type == GretlType.string_ref) {
        arg.val.str = cast(char *) p;
    } else if (type == GretlType.list) {
        arg.val.list = cast(int *) p;
    } else if (type == GretlType.series_ref ||
               type == GretlType.useries) {
        arg.val.idnum = *cast(int *) p;
    } else if (type == GretlType.bundle ||
               type == GretlType.bundle_ref) {
        arg.val.b = cast(GretlBundle *) p;
    } else if (type == GretlType.array ||
               type == GretlType.array_ref) {
        arg.val.a = cast(GretlArray *) p;
        adjust_array_arg_type(arg);
    } else {
        err = GretlError.types;
    }

    return err;
}

int fncall_add_args_array(Fncall * fc)
{
    int i;
    int np = fc.fun.n_params;
    int err = 0;

    fc.args = cast(FnArg *) malloc(np * (*fc.args).sizeof);

    if (fc.args == null) {
        err = GretlError.alloc;
    } else {
        for (i=0; i<np; i++) {
            fc.args[i].type = GretlType.none;
            fc.args[i].shifted = 0;
            fc.args[i].upname = null;
            fc.args[i].uvar = null;
        }
    }

    return err;
}

void fncall_clear_args_array(Fncall * fc)
{
    int i;
    int np = fc.fun.n_params;
    FnArg * arg;

    for (i=0; i<np; i++) {
        arg = &fc.args[i];
        arg.type = GretlType.none;
        arg.shifted = 0;
        nullify_upname(arg);
        arg.uvar = null;
        arg.val.px = null;
    }

    fc.argc = 0;
}

void fncall_destroy_args_array(Fncall * fc)
{
    if (fc.args != null) {
	int i;

	for (i=0; i<fc.fun.n_params; i++) {
	    if (fc.args[i].upname != null) {
		free(fc.args[i].upname);
	    }
	}
	free(fc.args);
    }
}

void maybe_set_param_const(FnParam * fp)
{
    fp.flags |= FpFlags.const_;
}

/**
 * push_function_arg:
 * @fc: pointer to function call.
 * @name: name of variable (or null for anonymous).
 * @uvar: reference to UserVar or null.
 * @type: type of argument to add.
 * @value: pointer to value to add.
 *
 * Writes a new argument of the specified type and value into the
 * argument array of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

extern (C) int push_function_arg(Fncall * fc, const(char) * name, void * uvar, GretlType type, void * value)
{
    int err = 0;

    if (fc == null || fc.fun == null) {
        err = GretlError.data;
    } else if (fc.argc >= fc.fun.n_params) {
        fprintf(stderr, "function %s has %d parameters but argc = %d\n",
                fc.fun.name, fc.fun.n_params, fc.argc);
        err = GretlError.data;
    } else if (fc.args == null) {
        err = fncall_add_args_array(fc);
    }

    if (!err) {
        err = fn_arg_set_data(&fc.args[fc.argc], name, cast(UserVar *) uvar, type, value);
        fc.argc += 1;
    }

    return err;
}

/**
 * push_anon_function_arg:
 * @fc: pointer to function call.
 * @type: type of argument to add.
 * @value: pointer to value to add.
 *
 * Writes a new argument of the specified type and value into the
 * argument array of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

extern (C) int push_anon_function_arg(Fncall * fc, GretlType type, void * value)
{
    return push_function_arg(fc, null, null, type, value);
}

/**
 * set_anon_function_arg:
 * @fc: pointer to function call.
 * @i: index of argument.
 * @type: type of argument.
 * @value: pointer to value to set.
 *
 * Sets type and value of (anonymous) argument @i in @fc. The index @i
 * must be >= 0 and less than the argc member of @fc.
 *
 * Returns: 0 on success, non-zero on failure.
 */

extern (C) int set_anon_function_arg(Fncall * fc, int i, GretlType type, void * value)
{
    if (i >= 0 && i < fc.argc) {
        fn_arg_set_data(&fc.args[i], null, null, type, value);
        return 0;
    } else {
        return GretlError.data;
    }
}

/**
 * push_function_args:
 * @fc: pointer to function call.
 *
 * Writes multiple entries into the argument array of @fc.
 * Each argument must be given in the form {type, value},
 * The list of entries must be terminated with -1.0
 *
 * Returns: 0 on success, non-zero on failure.
 */

extern (C) int push_function_args(Fncall * fc, ...)
{
    va_list ap;
    GretlType argtype;
    void * value;
    int err = 0;

    va_start(ap, fc);
    while (err == 0) {
        argtype = cast(GretlType) va_arg!(int)(ap);
        if (cast(int) argtype < 0) {
            /* reached the end of the args */
            break;
        }
        value = va_arg!(void *)(ap);
        err = push_function_arg(fc, null, null, argtype, value);
    }
    va_end(ap);

    return err;
}

/* fncall_new: if @preserve is non-zero, we don't destroy the Fncall
   struct automatically on completion of function execution -- in
   which case the caller of this function is responsible for calling
   fncall_destroy() at a suitable point in the proceedings.
*/

extern (C) Fncall * fncall_new(Ufunc * fun, int preserve)
{
    Fncall * call = cast(Fncall *) malloc(Fncall.sizeof);

    if (call != null) {
        call.fun = fun;
        call.ptrvars = null;
        call.listvars = null;
        call.lists = null;
        call.retname = null;
        call.rtype = fun.rettype;
        call.argc = 0;
        call.args = null;
        call.flags = preserve ? FCFlags.preserve : 0;
        call.line = null;
        call.lgen = null;
        call.n_lgen = 0;
    }

    return call;
}

int call_is_in_use(Fncall * fc)
{
    if (callstack != null) {
        GList *tmp = g_list_last(callstack);

        while (tmp != null) {
            if (fc == tmp.data) {
                return 1;
            }
            tmp = tmp.prev;
        }
    }

    return 0;
}

void clear_fncall_data(Fncall * fc)
{
    if (fc.args != null) {
        fncall_clear_args_array(fc);
    }
    if (fc.ptrvars != null) {
        free(fc.ptrvars);
        fc.ptrvars = null;
    }
    if (fc.listvars != null) {
        free(fc.listvars);
        fc.listvars = null;
    }
    if (fc.lists != null) {
        g_list_free(fc.lists);
        fc.lists = null;
    }
    free(fc.retname);
    fc.retname = null;
    fc.line = null;
}

/* This function is called (only) from geneval.c, in the function
   eval_ufunc(). The primary role of the Fncall struct is to serve as
   a vehicle for conveying the arguments to the user function in
   question, but it also carries state information about the call,
   including whether or not it is subject to recursion.
*/

extern (C) Fncall * user_func_get_fncall(Ufunc * fun)
{
    GList *tmp = g_list_first(fun.calls);
    Fncall * fc = null;

    while (tmp != null) {
        if (!call_is_in_use(cast(Fncall *) tmp.data)) {
            fc = cast(Fncall *) tmp.data;
            break;
        }
        tmp = tmp.next;
    }

    if (fc != null) {
        clear_fncall_data(fc);
        fun.call = fc;
    } else {
        fun.call = fncall_new(fun, 1);
        fun.calls = g_list_append(fun.calls, fun.call);
    }


    return fun.call;
}

extern (C) void fncall_destroy(void * data)
{
    if (data != null) {
        Fncall * call = cast(Fncall *) data;

        fncall_destroy_args_array(call);
        free(call.ptrvars);
        free(call.listvars);
        free(call.retname);
        g_list_free(call.lists);
        if (call.n_lgen > 0) {
            int i;

            for (i=0; i<call.n_lgen; i++) {
                destroy_genr(call.lgen[i].genr);
            }
            free(call.lgen);
        }
        free(call);
    }
}

void maybe_destroy_fncall(Fncall ** pcall)
{
    Fncall * call = *pcall;

    if (call != null && !(call.flags & FCFlags.preserve)) {
        if (call.fun != null) {
            if (call.fun.calls != null) {
                call.fun.calls = g_list_remove(call.fun.calls, call);
            }
            if (call.fun.call == call) {
                call.fun.call = null;
            }
        }
        fncall_destroy(call);
        *pcall = null;
    }
}

/* For use with overloaded functions whose return type is
   'numeric'. The specific return type is not known
   until execution is completed; we can supply it here.
*/

extern (C) GretlType fncall_get_return_type(Fncall * call)
{
    if (call != null) {
        return call.rtype;
    }

    return GretlType.none;
}

/* Portmanteau function to get a caller struct for a function named
   @funcname from a function package named @pkgname. We first check if
   the specified package is already in memory; if not we try to find
   it in the local filesystem, and load it into memory if successful.

   If/once the package is in fact loaded we look up the specified
   function; and if that's successful we allocate a caller struct and
   return it.

   The @pkgpath argument can be given as null, but if the path to the
   package is known to the caller and is provided via this argument
   this will speed up the look-up in case the package is not already
   loaded.
*/

extern (C) Fncall * get_pkg_function_call(const(char) * funcname, const(char) * pkgname, const(char) * pkgpath)
{
    Fncall * fc = null;
    Ufunc * uf = null;
    Fnpkg * pkg;

    /* is the package already loaded? */
    pkg = get_function_package_by_name(pkgname);

    if (pkg == null) {
        /* no, so look it up */
        int err = 0;

        if (pkgpath != null) {
            /* path was supplied by caller */
            pkg = get_function_package_by_filename(pkgpath, &err);
        } else {
            /* we need to search */
            char * mypath;

            mypath = gretl_function_package_get_path(pkgname, PkgType.all);
            if (mypath != null) {
                pkg = get_function_package_by_filename(mypath, &err);
                free(mypath);
            }
        }
    }

    if (pkg != null) {
        uf = get_function_from_package(funcname, pkg);
    }

    if (uf == null) {
        gretl_errmsg_sprintf(_("Couldn't find function %s").ptr, funcname);
    } else {
        fc = fncall_new(uf, 0); /* FIXME second arg? */
    }

    return fc;
}

Fnpkg * function_package_alloc(const(char) * fname)
{
    Fnpkg * pkg = cast(Fnpkg *) malloc(Fnpkg.sizeof);

    if (pkg == null) {
        return null;
    }

    pkg.fname = gretl_strdup(fname);
    if (pkg.fname == null) {
        free(pkg);
        return null;
    }


    pkg.name[0] = '\0';
    pkg.author = null;
    pkg.email = null;
    pkg.version_ = null;
    pkg.date = null;
    pkg.descrip = null;
    pkg.help = null;
    pkg.gui_help = null;
    pkg.Rdeps = null;
    pkg.sample = null;
    pkg.help_fname = null;
    pkg.gui_help_fname = null;
    pkg.sample_fname = null;
    pkg.tags = null;
    pkg.label = null;
    pkg.mpath = null;
    pkg.dreq = DataReq.needs_data;
    pkg.modelreq = 0;
    pkg.minver = 0;
    pkg.uses_subdir = 0;
    pkg.prechecked = 0;
    pkg.data_access = 0;

    pkg.pub = pkg.priv = null;
    pkg.n_pub = pkg.n_priv = 0;
    pkg.overrides = 0;
    pkg.datafiles = null;
    pkg.n_files = 0;
    pkg.depends = null;
    pkg.n_depends = 0;
    pkg.provider = null;
    pkg.conflict = null;
    pkg.editor = null;

    return pkg;
}

/* For function call @call, set the 'listarg' status of any named
   series provided to the called function via list arguments.  This is
   required for correct handling of namespaces.
*/

void set_listargs_from_call(Fncall * call, DatasetRaw * dset)
{
    int i;
    int vi;

    if (dset == null) {
        return;
    }

    for (i=1; i<dset.v; i++) {
        series_unset_flag(dset, i, VarFlags.listarg);
    }

    if (call != null && call.listvars != null) {
        for (i=1; i<=call.listvars[0]; i++) {
            vi = call.listvars[i];
            series_set_flag(dset, vi, VarFlags.listarg);
        }
    }
}

void get_prior_call_for_function(Ufunc * fun)
{
    GList *tmp = g_list_last(callstack);
    Fncall * call;

    tmp = tmp.prev;
    while (tmp != null) {
        call = cast(Fncall *) tmp.data;
        if (call.fun == fun) {
            /* reattach */
            fun.call = call;
            break;
        }
        tmp = tmp.prev;
    }
}

void set_executing_off(Fncall ** pcall, DatasetRaw * dset, Prn * prn)
{
    int dbg = gretl_debugging_on();
    Fncall * thiscall = *pcall;
    Fncall * popcall = null;
    Ufunc * fun = thiscall.fun;


    destroy_option_params_at_level(fn_executing);
    set_previous_depth(fn_executing);
    fn_executing--;

    if (is_recursing(thiscall)) {
        get_prior_call_for_function(fun);
    }

    callstack = g_list_remove(callstack, thiscall);


    if (dbg) {
        pputs(prn, "*** ");
        bufspace(gretl_function_depth(), prn);
        pprintf(prn, "exiting function %s, ", fun.name);
    }

    maybe_destroy_fncall(pcall);

    if (fn_executing > 0) {
        GList *tmp = g_list_last(callstack);

        popcall = cast(Fncall *) tmp.data;
    } else {
        g_list_free(callstack);
        callstack = null;
        gretl_insert_builtin_string("pkgdir", null);
    }

    if (dset != null) {
        set_listargs_from_call(popcall, dset);
    }

    if (popcall == null) {
        /* returning to main */
        switch_uservar_hash(0);
        if (dset != null) {
            series_ensure_level_zero(dset);
        }
    }

    if (dbg) {
        if (popcall != null) {
            pprintf(prn, "returning to %s\n", popcall.fun.name);
        } else {
            pputs(prn, "returning to main\n");
        }
    }
}

/**
 * n_user_functions:
 *
 * Returns: the number of hansl functions currently loaded in memory.
 */

extern (C) int n_user_functions()
{
    return n_ufuns;
}

/**
 * n_free_functions:
 *
 * Returns: the number of functions loaded in memory
 * that are not currently attached to any function package,
 * and are therefore available for packaging.
 */

extern (C) int n_free_functions()
{
    int i;
    int n = 0;

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i].pkg == null) {
            n++;
        }
    }

    return n;
}

/**
 * get_user_function_by_index:
 * @idx: index number.
 *
 * Returns: pointer to the user-function that currently
 * occupies (0-based) slot @idx in the array of loaded
 * functions, or %null if @idx is out of bounds.
 */

extern (C) Ufunc * get_user_function_by_index(int idx)
{
    return (idx < 0 || idx >= n_ufuns)? null : ufuns[idx];
}

/**
 * fn_n_params:
 * @fun: pointer to user-function.
 *
 * Returns: the number of parameters associated with @fun.
 */

extern (C) int fn_n_params(Ufunc * fun)
{
    return fun.n_params;
}

/**
 * fn_param_type:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the type of parameter @i of function
 * @fun.
 */

extern (C) int fn_param_type(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? 0 :
        fun.params[i].type;
}

/**
 * fn_param_name:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the name of parameter @i of function
 * @fun.
 */

extern (C) const(char) * fn_param_name(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? null :
        fun.params[i].name;
}

/**
 * fn_param_descrip:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the description of parameter @i of function
 * @fun (if any), otherwise %null.
 */

extern (C) const(char) * fn_param_descrip(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? null :
        fun.params[i].descrip;
}

/**
 * fn_param_value_labels:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 * @n: location to receive number of labels.
 *
 * Returns: the value-labels associated with parameter @i
 * of function @fun (if any), otherwise %null.
 */

extern (C) const(char) ** fn_param_value_labels(Ufunc * fun, int i, int * n)
{
    if (i >= 0 && i < fun.n_params) {
        *n = fun.params[i].nlabels;
        return cast(const(char) **) fun.params[i].labels;
    } else {
        *n = 0;
        return null;
    }
}

/**
 * fn_param_has_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if the (scalar) parameter @i of function @fun
 * (if any) has a default value set, otherwise 0.0
 */

extern (C) int fn_param_has_default(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return 0;
    } else {
        FnParam * fp = &fun.params[i];

        return !default_unset(fp);
    }
}

/**
 * fn_param_default:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the default value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

extern (C) double fn_param_default(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return NADBL;
    } else {
        FnParam * fp = &fun.params[i];

        return default_unset(fp)? NADBL : fp.deflt;
    }
}

/**
 * fn_param_minval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the minimum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

extern (C) double fn_param_minval(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? NADBL :
        fun.params[i].min;
}

/**
 * fn_param_maxval:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the maximum value of (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

extern (C) double fn_param_maxval(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? NADBL :
        fun.params[i].max;
}

/**
 * fn_param_step:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: the step value for (scalar) parameter @i of
 * function @fun (if any), otherwise #NADBL.
 */

extern (C) double fn_param_step(Ufunc * fun, int i)
{
    return (i < 0 || i >= fun.n_params)? NADBL :
        fun.params[i].step;
}

/* An "optional" parameter ("null" is OK and omission is
   also OK) must be a pointery thing (not a scalar).
*/

int arg_may_be_optional(GretlType t)
{
    return gretl_ref_type(t) || gretl_is_array_type(t) ||
        t == GretlType.series ||
        t == GretlType.matrix ||
        t == GretlType.bundle ||
        t == GretlType.list ||
        t == GretlType.string_;
}

/**
 * fn_param_automatic:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is marked "auto",
 * meaning that it is optional and has a default value,
 * otherwise 0.0
 */

extern (C) int fn_param_automatic(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return 0;
    }

    return arg_may_be_optional(fun.params[i].type) &&
        (fun.params[i].flags & FpFlags.auto_);
}

/**
 * fn_param_optional:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is optional,
 * otherwise 0.0
 */

extern (C) int fn_param_optional(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return 0;
    }

    return arg_may_be_optional(fun.params[i].type) &&
        (fun.params[i].flags & FpFlags.optional);
}

/**
 * fn_param_uses_xlist:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is
 * designed to select an integer based on a gretl
 * model's list of regressors, otherwise 0.0
 */

extern (C) int fn_param_uses_xlist(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return 0;
    }

    return (fun.params[i].type == GretlType.int_ &&
            fun.params[i].deflt == INT_USE_XLIST);
}

/**
 * fn_param_uses_mylist:
 * @fun: pointer to user-function.
 * @i: 0-based parameter index.
 *
 * Returns: 1 if parameter @i of function @fun is
 * designed to select an integer based on a custom
 * list, constructed by the function, otherwise 0.0
 */

extern (C) int fn_param_uses_mylist(Ufunc * fun, int i)
{
    if (i < 0 || i >= fun.n_params) {
        return 0;
    }

    return (fun.params[i].type == GretlType.int_ &&
            fun.params[i].deflt == INT_USE_MYLIST);
}

/**
 * user_func_get_return_type:
 * @fun: pointer to user-function.
 *
 * Returns: the return type of function @fun.
 */

extern (C) int user_func_get_return_type(Ufunc * fun)
{
    if (fun == null) {
        return GretlType.none;
    } else {
        return fun.rettype;
    }
}

/**
 * user_func_is_noprint:
 * @fun: pointer to user-function.
 *
 * Returns: 1 if the function is not designed to print anything.
 */

extern (C) int user_func_is_noprint(Ufunc * fun)
{
    if (fun == null) {
        return 0;
    } else {
        return function_is_noprint(fun);
    }
}

/**
 * user_func_is_menu_only:
 * @fun: pointer to user-function.
 *
 * Returns: 1 if the function is not designed to be called
 * other than via a GUI menu.
 */

extern (C) int user_func_is_menu_only(Ufunc * fun)
{
    if (fun == null) {
        return 0;
    } else {
        return function_is_menu_only(fun);
    }
}

/**
 * user_func_must_assign:
 * @fun: pointer to user-function.
 *
 * Returns: 1 if the function return value must be assigned in
 * GUI usage.
 */

extern (C) int user_func_must_assign(Ufunc * fun)
{
    if (fun == null) {
        return 0;
    } else {
        return function_must_assign(fun);
    }
}

/**
 * user_func_add_docstr:
 * @fun: pointer to user-function.
 * @str: the string to append.
 *
 * Adds @str as a line to @fun's docstr.

 * Returns: 0 on success, non-zero on failure.
 */

extern (C) int user_func_add_docstr(Ufunc * fun, const(char) * str)
{
    if (fun == null) {
        return GretlError.data;
    } else if (fun.docstr == null) {
	fun.docstr = g_strdup(str);
    } else {
	gchar *tmp = g_strdup_printf("%s\n%s", fun.docstr, str);

	g_free(fun.docstr);
	fun.docstr = tmp;
    }

    return 0;
}

/**
 * user_func_get_docstr:
 * @fun: pointer to user-function.
 *
 * Returns: the doc string set on to @fun, or null if
 * nothing has been set.
 */

extern (C) const(char) * user_func_get_docstr(Ufunc * fun)
{
    if (fun != null) {
        return fun.docstr;
    } else {
	return null;
    }
}

/**
 * user_function_name_by_index:
 * @i: the position of a user-function in the array of
 * loaded functions.
 *
 * Returns: the name of the function, or %null if
 * @i is out of bounds.
 */

extern (C) const(char) * user_function_name_by_index(int i)
{
    if (i >= 0 && i < n_ufuns) {
        return ufuns[i].name;
    } else {
        return null;
    }
}

/**
 * user_function_index_by_name:
 * @name: function name.
 * @pkg: reference function package, or null.
 *
 * Looks up the 0-based index of the named function
 * in the current array of user-functions. If @pkg is
 * non-null the search for @name is confined to functions
 * that belong to @pkg; otherwise the index of the first
 * match is returned.
 *
 * Returns: 0-based index or -1 on failure.
 */

extern (C) int user_function_index_by_name(const(char) * name, Fnpkg * pkg)
{
    int i;

    for (i=0; i<n_ufuns; i++) {
        if ((pkg == null || ufuns[i].pkg == pkg) &&
            !strcmp(name, ufuns[i].name)) {
            return i;
        }
    }

    return -1;
}

__gshared int fname_idx;

extern (C) void function_names_init()
{
    fname_idx = 0;
}

/* Apparatus used in the GUI selector for composing a new function
   package, or editing an existing one.  In the first case we want a
   list of names of currently unpackaged functions (and the @pkg
   argument should be null); in the second the list should include
   both unpackaged functions and those belonging to the package
   in question (specified via the @pkg arg).

   The caller should first invoke function_names_init(), then keep
   calling next_available_function_name() until it returns %null. The
   pointer argument @idxp provides a means to grab the "index number"
   (position in the current functions array) corresponding to the
   returned function name.
*/

extern (C) const(char) * next_available_function_name(Fnpkg * pkg, int * idxp)
{
    const(char) * ret = null;
    Ufunc * fun;

    if (n_ufuns == 0) {
        fname_idx = 0;
        return null;
    }

    while (fname_idx < n_ufuns) {
        fun = ufuns[fname_idx++];
        if (fun.pkg == null || fun.pkg == pkg) {
            ret = fun.name;
            *idxp = fname_idx - 1;
            break;
        }
    }

    return ret;
}

/* end GUI function-name apparatus */

Fncall * current_function_call()
{
    if (callstack != null) {
        GList *tmp = g_list_last(callstack);

        return cast(Fncall *) tmp.data;
    } else {
        return null;
    }
}

Ufunc * currently_called_function()
{
    Fncall *call = current_function_call();

    return (call != null)? call.fun : null;
}

extern (C) void current_function_info(const(char) ** funcname, const(char) ** pkgname)
{
    Ufunc *u = currently_called_function();

    if (u != null) {
        if (funcname != null) {
            *funcname = u.name;
        }
        if (pkgname != null && u.pkg != null) {
            *pkgname = u.pkg.name;
        }
    }
}

extern (C) int function_is_executing(const(char) * funcname)
{
    int ret = 0;

    if (callstack != null) {
        GList *tmp = g_list_last(callstack);
        Fncall * call;

        while (tmp != null) {
            call = cast(Fncall *) tmp.data;
            if (!strcmp(call.fun.name, funcname)) {
                ret = 1;
                break;
            }
            tmp = tmp.prev;
        }
    }

    return ret;
}

extern (C) Fnpkg * get_active_function_package(Gretlopt opt)
{
    Ufunc *u = currently_called_function();

    if (u != null && u.pkg != null) {
        if (opt == Gretlopt.OPT_NONE) {
            return u.pkg;
        } else if ((opt & Gretlopt.OPT_O) && u.pkg.overrides) {
            /* in this case we're only interested if the
               package overrides any built-in functions
            */
            return u.pkg;
        }
    }

    return null;
}

extern (C) Fnpkg * gretl_function_get_package(Ufunc * fun)
{
    return fun == null ? null : fun.pkg;
}


/* see if a function is currently employed in the call stack */

int function_in_use(Ufunc * fun)
{
    GList * tmp = callstack;
    Fncall * call;

    while (tmp != null) {
        call = cast(Fncall *) tmp.data;
        if (fun == call.fun) {
            return 1;
        }
        tmp = tmp.next;
    }

    return 0;
}

extern (C) int gretl_function_recursing()
{
    if (callstack != null) {
        GList *tmp = g_list_last(callstack);
        Fncall * call = cast(Fncall *) tmp.data;

        return is_recursing(call);
    } else {
        return 0;
    }
}


__gshared Fnpkg * bundle_pkg;

extern (C) void set_bundle_pkg(const(char) * pkgname)
{
    if (pkgname == null) {
	bundle_pkg = null;
    } else {
	bundle_pkg = get_function_package_by_name(pkgname);
    }
}

/**
 * get_user_function_by_name:
 * @name: name to test.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given name and it is accessible in
 * context (i.e. it's not private to a package other than
 * the one that's currently active, if any), otherwise
 * %null.
 */

extern (C) Ufunc * get_user_function_by_name(const(char) * name)
{
    Fnpkg * pkg = null;
    Ufunc * fun = null;
    int i;

    if (n_ufuns == 0) {
        return null;
    }

    if (bundle_pkg != null) {
	pkg = bundle_pkg;
    } else {
	pkg = current_pkg;
    }

    if (pkg == null) {
        fun = currently_called_function();
        if (fun != null) {
            pkg = fun.pkg;
            fun = null;
        }
    }


    if (pkg != null) {
        /* There's an active function package: try first
           for functions that belong to the package.
        */
        for (i=0; i<pkg.n_pub; i++) {
            /* public members */
            if (!strcmp(name, pkg.pub[i].name)) {
                fun = pkg.pub[i];
                break;
            }
        }
        if (fun == null) {
            /* private members */
            for (i=0; i<pkg.n_priv; i++) {
                if (!strcmp(name, pkg.priv[i].name)) {
                    fun = pkg.priv[i];
                    break;
                }
            }
        }
        if (fun == null && pkg.provider != null) {
            /* functions shared by provider */
            Fnpkg *prv = get_function_package_by_name(pkg.provider);

            if (prv != null) {
                for (i=0; i<prv.n_priv; i++) {
                    if (!strcmp(name, prv.priv[i].name)) {
                        fun = prv.priv[i];
                        break;
                    }
                }
            }
        }
    }

    if (fun == null) {
        /* Match any function (WAS any non-private function) */
        for (i=0; i<n_ufuns; i++) {
            if (!function_is_private(ufuns[i]) &&
                !strcmp(name, ufuns[i].name)) {
                fun = ufuns[i];
                break;
            }
        }
    }


    return fun;
}

/**
 * is_user_function:
 * @name: name to test.
 *
 * Returns: 1 if @name identifies a hansl function accessible
 * in the current context, otherwise 0.0
 */

extern (C) int is_user_function(const(char) * name)
{
    return get_user_function_by_name(name) != null;
}

/**
 * get_function_from_package:
 * @funname: name of function to retrieve.
 * @pkg: function package.
 *
 * Returns: pointer to a user-function, if there exists a
 * function of the given @funname that is associated with
 * function package @pkg, otherwise null.  This is used
 * in the gretl function package editor.
 */

extern (C) Ufunc * get_function_from_package(const(char) * funname, Fnpkg * pkg)
{
    if (pkg != null) {
        int i;

        for (i=0; i<pkg.n_pub; i++) {
            if (!strcmp(funname, pkg.pub[i].name)) {
                return pkg.pub[i];
            }
        }
        for (i=0; i<pkg.n_priv; i++) {
            if (!strcmp(funname, pkg.priv[i].name)) {
                return pkg.priv[i];
            }
        }
    }

    return null;
}

/* allocate and initialize a new array of @n parameters */

FnParam * allocate_params(int n)
{
    FnParam * params;
    int i;

    params = cast(FnParam *) malloc(n * (*params).sizeof);
    if (params == null) {
        return null;
    }

    for (i=0; i<n; i++) {
        params[i].name = null;
        params[i].type = GretlType.none;
        params[i].descrip = null;
        params[i].labels = null;
        params[i].nlabels = 0;
        params[i].flags = 0;
        params[i].deflt = UNSET_VALUE;
        params[i].min = NADBL;
        params[i].max = NADBL;
        params[i].step = NADBL;
    }

    return params;
}

Ufunc * ufunc_new()
{
    Ufunc * fun = cast(Ufunc *) malloc(Ufunc.sizeof);

    if (fun == null) {
        return null;
    }

    fun.name[0] = '\0';
    fun.pkg = null;
    fun.call = null;
    fun.calls = null;
    fun.flags = 0;
    fun.pkg_role = 0;

    fun.n_lines = 0;
    fun.line_idx = 1;
    fun.lines = null;

    fun.n_params = 0;
    fun.params = null;

    fun.rettype = GretlType.none;
    fun.docstr = null;

    return fun;
}

void free_lines_array(Stmt * lines, int n)
{
    Stmt * line;
    int i;

    if (lines == null) {
        return;
    }


    for (i=0; i<n; i++) {
        line = &lines[i];
        free(line.s);
        if (line_has_loop(line)) {
            gretl_loop_destroy(line.ptr);
        }
    }

    free(lines);
}

void free_params_array(FnParam * params, int n)
{
    int i;

    if (params == null) return;

    for (i=0; i<n; i++) {
        free(params[i].name);
        free(params[i].descrip);
        strings_array_free(params[i].labels, params[i].nlabels);
    }
    free(params);
}

void destroy_ufunc_calls(Ufunc * fun)
{
    if (fun.calls != null) {
        g_list_free_full(fun.calls, &fncall_destroy);
        fun.calls = null;
        fun.call = null;
    } else if (fun.call != null) {
        fncall_destroy(fun.call);
        fun.call = null;
    }
}

void clear_ufunc_data(Ufunc * fun)
{
    free_lines_array(fun.lines, fun.n_lines);
    free_params_array(fun.params, fun.n_params);

    fun.lines = null;
    fun.params = null;

    fun.n_lines = 0;
    fun.line_idx = 1;
    fun.n_params = 0;

    destroy_ufunc_calls(fun);

    fun.rettype = GretlType.none;
    g_free(fun.docstr);
    fun.docstr = null;
}

void ufunc_free(Ufunc * fun)
{
    /* free_lines_array() handles any attached loops */
    free_lines_array(fun.lines, fun.n_lines);
    free_params_array(fun.params, fun.n_params);
    /* destroy_ufunc_calls() handles attached genrs */
    destroy_ufunc_calls(fun);
    g_free(fun.docstr);
    free(fun);
}

int add_allocated_ufunc(Ufunc * fun)
{
    int nf = n_ufuns;
    Ufunc **myfuns;

    myfuns = cast(Ufunc* *) realloc(ufuns, (nf + 1) * (*myfuns).sizeof);

    if (myfuns == null) {
        return GretlError.alloc;
    }

    ufuns = myfuns;
    ufuns[nf] = fun;
    n_ufuns++;


    return 0;
}

Ufunc * add_ufunc(const(char) * fname)
{
    Ufunc *fun = ufunc_new();


    if (fun != null) {
        strncat(fun.name, fname, FN_NAMELEN - 1);
        if (add_allocated_ufunc(fun)) {
            ufunc_free(fun);
            fun = null;
        }
    }

    return fun;
}

int no_scalar_default(FnParam * fp)
{
    int ret = 0;

    if (default_unset(fp)) {
        ret = 1;
    } else if (fp.type != GretlType.double_ && na(fp.deflt)) {
        ret = 1;
    }

    return ret;
}

/* handling of XML function packages */

enum FuncsTask {
    info,
    load,
    code,
    sample,
    help,
    query
}

const(char) * arg_type_xml_string(int t)
{
    if (t == GretlType.scalar_ref) {
        return "scalarref";
    } else if (t == GretlType.series_ref) {
        return "seriesref";
    } else if (t == GretlType.matrix_ref) {
        return "matrixref";
    } else if (t == GretlType.bundle_ref) {
        return "bundleref";
    } else if (t == GretlType.string_ref) {
        return "stringref";
    } else if (t == GretlType.strings_ref) {
        return "stringsref";
    } else if (t == GretlType.matrices_ref) {
        return "matricesref";
    } else if (t == GretlType.bundles_ref) {
        return "bundlesref";
    } else if (t == GretlType.lists_ref) {
        return "listsref"; /* not actually allowed */
    } else {
        return gretl_type_get_name(cast(GretlType) t);
    }
}

GretlType return_type_from_string(const(char) * s, int * err)
{
    GretlType t;

    if (!strcmp(s, "void")) {
        /* not OK as arg type, but OK as return */
        t = GretlType.void_;
    } else {
        t = gretl_type_from_string(s);
    }

    if (!ok_function_return_type(t)) {
        if (*s == '\0') {
            gretl_errmsg_sprintf(_("Missing function return type").ptr);
        } else if (t == GretlType.none) {
            gretl_errmsg_sprintf(_("Expected a function return type, found '%s'").ptr,
                                 s);
        } else {
            gretl_errmsg_sprintf(_("%s: invalid return type for function").ptr,
                                 s);
        }
        *err = GretlError.types;
    }

    return t;
}

GretlType param_field_to_type(const(char) * s, const(char) * funname, int * err)
{
    GretlType t = gretl_type_from_string(s);

    if (!ok_function_arg_type(t)) {
        gretl_errmsg_sprintf(_("function %s: invalid parameter type '%s'").ptr,
                             funname, s);
        *err = GretlError.invarg;
    }

    return t;
}

gboolean special_int_default(FnParam * param, xmlNodePtr np)
{
    gboolean ret = 0;

    if (param.type == GretlType.int_) {
        char * s = null;

        gretl_xml_get_prop_as_string(np, "default", &s);
        if (s != null) {
            if (strstr(s, "$mylist")) {
                param.deflt = INT_USE_MYLIST;
                ret = 1;
            } else if (strstr(s, "$xlist")) {
                param.deflt = INT_USE_XLIST;
                ret = 1;
            }
            free(s);
        }
    }

    return ret;
}

/* read the parameter info for a function from XML file */

int func_read_params(xmlNodePtr node, xmlDocPtr doc, Ufunc * fun)
{
    xmlNodePtr cur;
    char * field;
    int n;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n) || n < 0) {
        fprintf(stderr, "Couldn't read param count\n");
        return GretlError.data;
    }

    if (n == 0) {
        return 0;
    }

    fun.params = allocate_params(n);
    if (fun.params == null) {
        return GretlError.alloc;
    }

    fun.n_params = n;

    gretl_push_c_numeric_locale();

    cur = node.children;
    n = 0;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "param".ptr)) {
            FnParam * param = &fun.params[n++];

            if (gretl_xml_get_prop_as_string(cur, "name", &field)) {
		if (strlen(field) >= FN_NAMELEN) {
		    gretl_errmsg_set(_("Identifier exceeds the maximum of 31 characters").ptr);
		    free(field);
		    err = GretlError.data;
		    break;
		}
                param.name = field;
            } else {
                err = GretlError.data;
                break;
            }
            if (gretl_xml_get_prop_as_string(cur, "type", &field)) {
                param.type = param_field_to_type(field, fun.name, &err);
                free(field);
                if (special_int_default(param, cur)) {
                    ; /* handled */
                } else if (gretl_scalar_type(param.type)) {
                    double x;

                    if (gretl_xml_get_prop_as_double(cur, "default", &x)) {
                        param.deflt = x;
                    } else {
                        param.deflt = UNSET_VALUE;
                    }
                    if (param.type != GretlType.bool_) {
                        gretl_xml_get_prop_as_double(cur, "min", &param.min);
                        gretl_xml_get_prop_as_double(cur, "max", &param.max);
                        gretl_xml_get_prop_as_double(cur, "step", &param.step);
                    }
                }
                if (gretl_xml_get_prop_as_bool(cur, "auto")) {
                    set_param_auto(param);
                } else if (gretl_xml_get_prop_as_bool(cur, "optional")) {
                    set_param_optional(param);
                }
                if (gretl_xml_get_prop_as_bool(cur, "const")) {
                    maybe_set_param_const(param);
                }
            } else {
                err = GretlError.data;
                break;
            }
            gretl_xml_child_get_string(cur, doc, "description",
                                       &param.descrip);
            gretl_xml_child_get_strings_array(cur, doc, "labels",
                                              &param.labels,
                                              &param.nlabels);
        }
        cur = cur.next;
    }

    gretl_pop_c_numeric_locale();

    if (!err && n != fun.n_params) {
        err = GretlError.data;
    }

    return err;
}

/* When we come here @fun may already have an array of lines allocated
   (if we know its size in advance). In that case @n_alloced will be
   non-null and will contain the number of lines allocated.  Otherwise
   the array of lines will be null initially, and @n_alloced will be
   null.

   In both cases fun.n_lines will be zero on the first call to this
   function, and it must be incremented each time a line is saved.
*/

int push_function_line(Ufunc * fun, char * s, int ci, int donate, int * n_alloced)
{
    Stmt * ln = null;
    int n = fun.n_lines + 1;
    int err = 0;

    if (n_alloced == null || n > *n_alloced) {
        /* create or enlarge the array of lines */
        Stmt * lines = cast(Stmt *) realloc(fun.lines, n * Stmt.sizeof);

        if (lines == null) {
            return GretlError.alloc;
        }
        fun.lines = lines;
        fun.n_lines = n;
        if (n_alloced != null) {
            *n_alloced = n;
        }
    } else {
        /* we already have sufficient storage */
        fun.n_lines += 1;
    }

    ln = &fun.lines[n-1];
    ln.idx = cast(short) fun.line_idx;

    if (donate) {
        /* take ownership of @s */
        ln.s = s;
    } else {
        /* make a copy of @s */
        ln.s = gretl_strdup(s);
        if (ln.s == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        ln.ptr = null;
        ln.ci = cast(short) ci;
        ln.next = 0;
        /* be on the safe side wrt string substitution */
        ln.flags = strchr(ln.s, '$') ? LnFlags.nocomp : 0;
        if (ci == CMD_COMMENT) {
            ln.flags |= LnFlags.ignore;
        }
        fun.line_idx += 1;
    }

    return err;
}

/* Get a quick reading on the number of non-blank lines in @buf, so we
   can pre-allocate all the lines needed to store a function
   definition. This may be an overestimate if some "blank" lines
   actually contain spaces.
*/

int buf_get_n_lines(const(char) * buf)
{
    const(char) * p = buf;
    const(char) * prev = null;
    int n = 1;

    while ((p = strchr(p, '\n')) != null) {
        if (prev != null && p - prev > 1) {
            n++;
        }
        prev = p++;
    }

    return n;
}

/* does a line that starts with "set" actually change
   the program state? heuristic
*/

int changes_state(const(char) * line)
{
    if (string_is_blank(line)) {
	return 0;
    } else {
	char[16] setarg;

	sscanf(line, "%15s", setarg.ptr);
	return strcmp(setarg.ptr, "stopwatch");
    }
}

/* Read the lines of hansl code from the <code> element of the XML
   representation of a function.
*/

int func_read_code(xmlNodePtr node, xmlDocPtr doc, Ufunc * fun)
{
    ExecState state;

    memset(&state, 0, ExecState.sizeof);
    Cmd cmd;

    memset(&cmd, 0, Cmd.sizeof);
    char[MAXLINE] line;
    char * buf;
    char * s;
    gint8 uses_set = 0;
    gint8 has_flow = 0;
    int save_comments = 1; /* added 2023-03-16 */
    int n_lines = 0;
    int err = 0;

    buf = cast(char *) xmlNodeListGetString(doc, node.children, 1);
    if (buf == null) {
        return 1;
    }

    n_lines = buf_get_n_lines(buf);
    bufgets_init(buf);
    gretl_cmd_init(&cmd);
    state.cmd = &cmd;

    /* pre-allocate storage */
    fun.lines = cast(Stmt *) calloc(n_lines, (*fun.lines).sizeof);

    while (bufgets(line.ptr, line.sizeof, buf)) {
        s = line.ptr;
        while (isspace(*s)) s++;
        state.line = s;
        err = get_command_index(&state, GretlCmdIndex.func, save_comments);
        if (err) {
            break;
        }
        if (cmd.ci == GretlCmdIndex.set && changes_state(s + 3)) {
            uses_set = 1;
        } else if (cmd.ci == GretlCmdIndex.if_ || cmd.ci == GretlCmdIndex.loop) {
            has_flow = 1;
        }
        tailstrip(s);
        if (string_is_blank(s)) {
            fun.line_idx += 1;
        } else {
            err = push_function_line(fun, s, cmd.ci, 0, &n_lines);
        }
    }

    if (fun.n_lines < n_lines) {
        /* shrink to avoid a memory leak later */
        fun.lines = cast(Stmt *) realloc(fun.lines, n_lines * (*fun.lines).sizeof);
    }

    if (uses_set) {
        fun.flags |= UfunAttrs.uses_set;
    }
    if (has_flow) {
        fun.flags |= UfunAttrs.has_flow;
    }

    bufgets_finalize(buf);
    free(buf);

    if (!err && has_flow) {
        err = ufunc_get_structure(fun);
    }

    gretl_cmd_free(&cmd);

    return err;
}

void print_opt_flags(FnParam * param, Prn * prn)
{
    if (param_is_auto(param)) {
        pputs(prn, "[auto]");
    } else if (param_is_optional(param)) {
        pputs(prn, "[null]");
    }
}

void print_param_description(FnParam * param, Prn * prn)
{
    if (param.descrip != null && *param.descrip != '\0') {
        pprintf(prn, " \"%s\"", param.descrip);
    }
}

void print_param_labels(FnParam * param, Prn * prn)
{
    int i;

    pputs(prn, " {");
    for (i=0; i<param.nlabels; i++) {
        pprintf(prn, "\"%s\"", param.labels[i]);
        if (i < param.nlabels - 1) {
            pputs(prn, ", ");
        }
    }
    pputc(prn, '}');
}

void print_min_max_deflt(FnParam * param, Prn * prn)
{
    if (na(param.min) && na(param.max) && default_unset(param)) {
        return; /* no-op */
    } else if (na(param.min) && na(param.max)) {
        /* got a default value only? */
        if (!default_unset(param)) {
            if (na(param.deflt)) {
                pputs(prn, "[NA]");
            } else if (param.deflt == INT_USE_XLIST) {
                pputs(prn, "[$xlist]");
            } else if (param.deflt == INT_USE_MYLIST) {
                pputs(prn, "[$mylist]");
            } else {
                pprintf(prn, "[%g]", param.deflt);
            }
        }
        return;
    }

    pputc(prn, '[');

    /* minimum */
    if (!na(param.min)) pprintf(prn, "%g", param.min);
    pputc(prn, ':');

    /* maximum */
    if (!na(param.max)) pprintf(prn, "%g", param.max);
    pputc(prn, ':');

    /* default */
    if (!default_unset(param)) {
        if (na(param.deflt)) {
            pputs(prn, "NA");
        } else {
            pprintf(prn, "%g", param.deflt);
        }
    }

    if (!na(param.step)) {
        /* step */
        pputc(prn, ':');
        pprintf(prn, "%g", param.step);
    }

    pputc(prn, ']');
}

/* free @fun and also remove it from the list of loaded
   functions */

void ufunc_unload(Ufunc * fun)
{
    int i;
    int j;
    int found = 0;

    if (n_ufuns == 0 || fun == null) {
        return;
    }


    /* remove this function from the array of loaded functions */

    for (i=0; i<n_ufuns; i++) {
        if (ufuns[i] == fun) {
            for (j=i; j<n_ufuns-1; j++) {
                ufuns[j] = ufuns[j+1];
            }
            found = 1;
            break;
        }
    }

    if (fun.pkg != null && fun.pkg.overrides) {
        delete_function_override(fun.name, fun.pkg.name);
    }
    ufunc_free(fun);

    if (found) {
        n_ufuns--;
    }
}

/* remove a package from the current listing and free it;
   if @full is non-zero, also unload any member functions
*/

void real_function_package_unload(Fnpkg * pkg, int full)
{
    int i;
    int j;
    int found = 0;


    /* unload children? */
    if (full) {
        for (i=0; i<pkg.n_priv; i++) {
            ufunc_unload(pkg.priv[i]);
        }
        for (i=0; i<pkg.n_pub; i++) {
            ufunc_unload(pkg.pub[i]);
        }
    }

    /* remove this package from the array of loaded packages */
    for (i=0; i<n_pkgs; i++) {
        if (pkgs[i] == pkg) {
            for (j=i; j<n_pkgs-1; j++) {
                pkgs[j] = pkgs[j+1];
            }
            found = 1;
            break;
        }
    }

    /* free the package itself */
    function_package_free(pkg);

    if (found) {
        n_pkgs--;
    }
}

/* Append a pointer to @fun to the array of child-pointers in @pkg: we
   do this when reading the definition of a packaged function from an
   XML file.  Note that this action does not add the functions to the
   array of loaded functions -- that's done separately, if we're
   loading the package 'for real'.
*/

int attach_ufunc_to_package(Ufunc * fun, Fnpkg * pkg)
{
    Ufunc **ufs;
    int n;
    int err = 0;

    if (function_is_private(fun)) {
        n = pkg.n_priv;
        ufs = cast(Ufunc* *) realloc(pkg.priv, (n + 1) * (*ufs).sizeof);
        if (ufs == null) {
            err = GretlError.alloc;
        } else {
            pkg.priv = ufs;
            pkg.priv[n] = fun;
            pkg.n_priv += 1;
        }
    } else {
        n = pkg.n_pub;
        ufs = cast(Ufunc* *) realloc(pkg.pub, (n + 1) * (*ufs).sizeof);
        if (ufs == null) {
            err = GretlError.alloc;
        } else {
            pkg.pub = ufs;
            pkg.pub[n] = fun;
            pkg.n_pub += 1;
        }
    }


    return err;
}

void maybe_set_menu_only(Ufunc * fun, Fnpkg * pkg)
{
    if (pkg.mpath != null && strstr(pkg.mpath, "MODELWIN")) {
        if (fun.pkg_role == UfunRole.gui_main) {
            fun.flags |= UfunAttrs.menu_only;
        }
    }
}

/* Read a single user-function definition from XML file: if the
   function is a child of a package, the @pkg argument will be
   non-null.
*/

int read_ufunc_from_xml(xmlNodePtr node, xmlDocPtr doc, Fnpkg * pkg)
{
    Ufunc *fun = ufunc_new();
    xmlNodePtr cur;
    char * tmp;
    int err = 0;

    if (fun == null) {
        return GretlError.alloc;
    }

    if (!gretl_xml_get_prop_as_string(node, "name", &tmp)) {
        ufunc_free(fun);
        return GretlError.data;
    }

    strncat(fun.name, tmp, FN_NAMELEN - 1);
    free(tmp);

    if (pkg != null) {
        fun.pkg = pkg;
    }

    if (gretl_xml_get_prop_as_string(node, "type", &tmp)) {
        fun.rettype = return_type_from_string(tmp, &err);
        free(tmp);
    } else {
        fun.rettype = GretlType.void_;
    }

    if (err) {
        ufunc_free(fun);
        return err;
    }

    if (gretl_xml_get_prop_as_bool(node, "private")) {
        fun.flags |= UfunAttrs.private_;
    }
    if (gretl_xml_get_prop_as_bool(node, "no-print")) {
        fun.flags |= UfunAttrs.noprint;
    }
    if (gretl_xml_get_prop_as_bool(node, "menu-only")) {
        fun.flags |= UfunAttrs.menu_only;
    }
    if (gretl_xml_get_prop_as_bool(node, "must-assign")) {
        fun.flags |= UfunAttrs.assign;
    }
    if (gretl_xml_get_prop_as_string(node, "pkg-role", &tmp)) {
        fun.pkg_role = pkg_key_get_role(tmp);
        free(tmp);
    }
    if (!(fun.flags & UfunAttrs.menu_only) && pkg != null) {
        maybe_set_menu_only(fun, pkg);
    }


    if (pkg == null && function_is_private(fun)) {
        fprintf(stderr, "unpackaged function: can't be private\n");
        ufunc_free(fun);
        return GretlError.data;
    }

    cur = node.children;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "help".ptr)) {
            /* backward compatibility: help used to be attached to
               a package's public interface */
            if (pkg.help != null) {
                free(pkg.help);
            }
            pkg.help = gretl_xml_get_string(cur, doc);
        } else if (!xmlStrcmp(cur.name, cast(XUC) "params".ptr)) {
            err = func_read_params(cur, doc, fun);
            if (err) {
                fprintf(stderr, "%s: error parsing function parameters\n",
                        fun.name);
            }
        } else if (!xmlStrcmp(cur.name, cast(XUC) "return".ptr)) {
            gretl_errmsg_set(_("Old-style function definitions no longer supported").ptr);
            err = GretlError.data;
        } else if (!xmlStrcmp(cur.name, cast(XUC) "code".ptr)) {
            err = func_read_code(cur, doc, fun);
        }
        cur = cur.next;
    }

    if (!err) {
        if (pkg != null) {
            /* function belongs to a package */
            err = attach_ufunc_to_package(fun, pkg);
        } else {
            /* reading a standalone function from session file */
            err = add_allocated_ufunc(fun);
        }
    }

    if (err) {
        ufunc_free(fun);
    }


    return err;
}

void print_param_limit(FnParam * param, int i, Prn * prn)
{
    if (param.type == GretlType.int_) {
        if (i == 0) {
            pprintf(prn, " min=\"%.0f\"", param.min);
        } else if (i == 1) {
            pprintf(prn, " max=\"%.0f\"", param.max);
        } else {
            pprintf(prn, " default=\"%.0f\"", param.deflt);
        }
    } else {
        if (i == 0) {
            pprintf(prn, " min=\"%.15g\"", param.min);
        } else if (i == 1) {
            pprintf(prn, " max=\"%.15g\"", param.max);
        } else {
            pprintf(prn, " default=\"%.15g\"", param.deflt);
        }
    }
}

bool parm_has_children(FnParam * p)
{
    return p.descrip != null || p.nlabels > 0;
}

/* write out a single user-defined function as XML, according to
   gretlfunc.dtd */

int write_function_xml(Ufunc * fun, Prn * prn, int mpi)
{

    GretlType rtype = fun.rettype;
    int this_indent = 0;
    int next_indent = 0;
    int i;

    if (rtype == GretlType.none) {
        rtype = GretlType.void_;
    }

    pprintf(prn, "<gretl-function name=\"%s\" type=\"%s\"",
            fun.name, gretl_type_get_name(rtype));

    if (function_is_private(fun)) {
        pputs(prn, " private=\"1\"");
    }
    if (function_is_noprint(fun)) {
        pputs(prn, " no-print=\"1\"");
    }
    if (function_is_menu_only(fun)) {
        pputs(prn, " menu-only=\"1\"");
    }
    if (function_must_assign(fun)) {
        pputs(prn, " must-assign=\"1\"");
    }

    if (fun.pkg_role) {
        pprintf(prn, " pkg-role=\"%s\"", package_role_get_key(fun.pkg_role));
    }

    pputs(prn, ">\n");

    if (fun.n_params > 0) {
        FnParam * param;

        gretl_push_c_numeric_locale();
        pprintf(prn, " <params count=\"%d\">\n", fun.n_params);
        for (i=0; i<fun.n_params; i++) {
            param = &fun.params[i];
            pprintf(prn, "  <param name=\"%s\" type=\"%s\"",
                    param.name, arg_type_xml_string(param.type));
            if (!na(param.min)) {
                print_param_limit(param, 0, prn);
            }
            if (!na(param.max)) {
                print_param_limit(param, 1, prn);
            }
            if (!default_unset(param)) {
                if (na(param.deflt)) {
                    pputs(prn, " default=\"NA\"");
                } else if (param.deflt == INT_USE_MYLIST) {
                    pputs(prn, " default=\"$mylist\"");
                } else if (param.deflt == INT_USE_XLIST) {
                    pputs(prn, " default=\"$xlist\"");
                } else {
                    print_param_limit(param, 2, prn);
                }
            }
            if (!na(param.step)) {
                pprintf(prn, " step=\"%g\"", param.step);
            }
            if (param_is_auto(param)) {
                pputs(prn, " auto=\"true\"");
            } else if (param_is_optional(param)) {
                pputs(prn, " optional=\"true\"");
            }
            if (param_is_const(param)) {
                pputs(prn, " const=\"true\"");
            }
            if (parm_has_children(param)) {
                pputs(prn, ">\n"); /* terminate opening tag */
                if (param.descrip != null) {
                    gretl_xml_put_tagged_string("description",
                                                param.descrip,
                                                prn);
                }
                if (param.nlabels > 0) {
                    gretl_xml_put_strings_array_quoted("labels",
                                                       cast(const(char) **) param.labels,
                                                       param.nlabels, prn);
                }
                pputs(prn, "  </param>\n");
            } else {
                pputs(prn, "/>\n"); /* terminate opening tag */
            }
        }
        pputs(prn, " </params>\n");
        gretl_pop_c_numeric_locale();
    }

    pputs(prn, "<code>");

    if (mpi) {
        /* minimal output of function code */
        for (i=0; i<fun.n_lines; i++) {
            gretl_xml_put_string(fun.lines[i].s, prn);
            pputc(prn, '\n');
        }
    } else {
        /* pay some attention to readability */
        for (i=0; i<fun.n_lines; i++) {
            if (i > 0 && fun.lines[i].idx - fun.lines[i-1].idx > 1) {
                /* reinstate single blank lines */
                pputc(prn, '\n');
            }
            /* basic indentation adjustment */
            adjust_indent(fun.lines[i].s, &this_indent, &next_indent);
            bufspace(2 * this_indent, prn);
            gretl_xml_put_string(fun.lines[i].s, prn);
            pputc(prn, '\n');
        }
    }

    pputs(prn, "</code>\n");
    pputs(prn, "</gretl-function>\n");

    return 0;
}

/* script-style output */

extern (C) void print_function_signature(Ufunc * fun, Prn * prn)
{
    const(char) * s;
    int i;
    int pos = 0;

    if (fun.rettype == GretlType.none) {
        pos += pprintf(prn, "function void %s ", fun.name);
    } else {
        const(char) * typestr = gretl_type_get_name(fun.rettype);

        pos += pprintf(prn, "function %s %s (", typestr, fun.name);
    }

    gretl_push_c_numeric_locale();

    if (fun.n_params == 0) {
        pputs(prn, "void)");
    }

    for (i=0; i<fun.n_params; i++) {
        FnParam * fp = &fun.params[i];

        if (fp.flags & FpFlags.const_) {
            pputs(prn, "const ");
        }
        s = gretl_type_get_name(fp.type);
        if (s[strlen(s) - 1] == '*') {
            pprintf(prn, "%s%s", s, fp.name);
        } else {
            pprintf(prn, "%s %s", s, fp.name);
        }
        if (fp.type == GretlType.bool_) {
            if (!default_unset(fp) && !na(fp.deflt)) {
                pprintf(prn, "[%g]", fp.deflt); /* FIXME? */
            }
        } else if (gretl_scalar_type(fp.type)) {
            print_min_max_deflt(fp, prn);
        } else if (arg_may_be_optional(fp.type)) {
            print_opt_flags(fp, prn);
        }
        print_param_description(fp, prn);
        if (fp.nlabels > 0) {
            print_param_labels(fp, prn);
        }
        if (i == fun.n_params - 1) {
            pputc(prn, ')');
        } else {
            pputs(prn, ",\n");
            bufspace(pos, prn);
        }
    }

    pputc(prn, '\n');
    gretl_pop_c_numeric_locale();
}

extern (C) char * gretl_function_get_signature(Ufunc * u)
{
    Prn * prn = null;
    char * buf = null;

    prn = gretl_print_new(PrnType.buffer, null);
    print_function_signature(u, prn);
    buf = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);

    return buf;
}

/* The following is an attempt at generic indentation based
   on occurrences of '{' and '}', which should probably
   work for R code.
*/

void foreign_adjust_indent(const(char) * s, int * this_indent, int * next_indent)
{
    int n = cast(int) strlen(s);

    *this_indent = *next_indent;

    if (strncmp(s, "} else", 6) == 0) {
        *next_indent = *this_indent;
        *this_indent = *this_indent - 1;
    } else if (s[n-1] == '{') {
        *next_indent = *this_indent + 1;
    } else if (strcmp(s, "}") == 0 ||
               strncmp(s, "} else", 6) == 0) {
        *this_indent = *next_indent = *this_indent - 1;
    }
}

void maybe_toggle_foreign(const(char) * s, int * in_foreign)
{
    if (!*in_foreign) {
        if (strncmp(s, "foreign ", 8) == 0) {
            *in_foreign = 1;
        }
    } else if (strncmp(s, "end foreign", 11) == 0) {
        *in_foreign = 0;
    }
}

/**
 * gretl_function_print_code:
 * @u: pointer to user-function.
 * @tabwidth: number of spaces per "tab" (logical indent).
 * @prn: printing struct.
 *
 * Prints out function @fun to @prn, script-style. Note
 * that this function may be called for several functions
 * in turn, with a single @prn as target.
 *
 * Returns: 0 on success, non-zero if @fun is %null.
 */

extern (C) int gretl_function_print_code(Ufunc * u, int tabwidth, Prn * prn)
{
    Prn * ptmp;
    const(char) * buf;
    int this_indent = 0;
    int next_indent = 0;
    int in_foreign = 0;
    int i;
    int j;

    if (u == null) {
        return GretlError.data;
    }

    if (tabwidth == 0) {
        tabwidth = 2;
    }

    ptmp = gretl_print_new(PrnType.buffer, null);
    print_function_signature(u, ptmp);

    for (i=0; i<u.n_lines; i++) {
        maybe_toggle_foreign(u.lines[i].s, &in_foreign);
        if (in_foreign) {
            foreign_adjust_indent(u.lines[i].s, &this_indent, &next_indent);
        } else {
            adjust_indent(u.lines[i].s, &this_indent, &next_indent);
        }
        for (j=0; j<=this_indent; j++) {
            bufspace(tabwidth, ptmp);
        }
        pputs(ptmp, u.lines[i].s);
        if (i < u.n_lines - 1) {
            if (u.lines[i+1].idx > u.lines[i].idx + 1) {
                pputc(ptmp, '\n');
            }
        }
        pputc(ptmp, '\n');
    }

    pputs(ptmp, "end function\n");

    /* now refine the indentation */
    buf = gretl_print_get_buffer(ptmp);
    normalize_hansl(buf, tabwidth, prn);
    gretl_print_destroy(ptmp);

    return 0;
}

/* called from kalman.c for time-varying matrices */

extern (C) char ** gretl_function_retrieve_code(Ufunc * u, int * nlines)
{
    char **S = null;
    Stmt * line;
    int i;
    int j = 0;

    for (i=0; i<u.n_lines; i++) {
        line = &(u.lines[i]);
        if (!ignore_line(line)) {
            j++;
        }
    }

    if (j > 0) {
        S = strings_array_new(j);
    }

    if (S != null) {
        *nlines = j;
        j = 0;
        for (i=0; i<u.n_lines; i++) {
            line = &(u.lines[i]);
            if (!ignore_line(line)) {
                S[j++] = line.s;
            }
        }
    }

    return S;
}
