/* Port of lib/src/cmd_private.h: the command-execution structures (CMD, ExecState, OpenOp, stmt) and the CmdFlags
   enumeration. The prototypes of the C header belong to cmd_private.c, which is not ported yet. Differences from the
   C: CMD is `Cmd`, ExecState keeps its name (interacth.d publicly imports it from here), stmt is `Stmt`,
   CmdFlags loses its CMD_ prefix (`catch` becomes `catch_`), the macros cmd_ccmt, cmd_subst and cmd_nosub are
   functions, and LOOPSET and cmd_token are opaque (defined in monte_carlo.c and tokenizer.c). */

module gretl.cmd_privateh;


import core.stdc.stdio : FILENAME_MAX;

import gretl.gretl_prn : Prn;
import gretl.gretl_restricth : GretlRestriction;
import gretl.interacth : ExecFlags, MAXSAVENAME;
import gretl.libgretlh;
import gretl.objstackh : GretlObjType;
import gretl.systemh : EquationSystem;
import gretl.varh : GretlVar;

enum COMPILE_FEVAL = 1;

enum CmdFlags {
    ccmt   = 1 << 0, /* line is in a C-style multi-line comment */
    subst  = 1 << 1, /* string substitution has been done on command */
    prog   = 1 << 2, /* command is in context of progressive loop */
    catch_ = 1 << 3, /* error from command should be "caught" */
    nosub  = 1 << 4, /* no @-substitution called for (pre-checked) */
    endfun = 1 << 5  /* line terminates a function definition */
}

struct CmdToken;   /* opaque: defined in tokenizer.c */
struct LoopSet;    /* opaque: defined in monte_carlo.c */

struct Cmd {
    int ci;          /* current command index */
    int err;         /* error code */
    int context;     /* for block commands, index of current context */
    int ciflags;     /* see CIFlags in tokenizer.c */
    Gretlopt opt;    /* option(s) for command */
    CmdFlags flags;  /* status flags for command invocation */
    GretlType gtype; /* specified type for "genr" */
    int order;       /* lag order, where appropriate */
    int auxint;      /* auxiliary int (e.g. VECM rank) */
    int cstart;      /* token index of start of 'real' command */
    int ntoks;       /* number of tokens actually used */
    int nt_alloced;  /* number of tokens allocated */
    CmdToken * toks; /* tokens */
    char * vstart;   /* portion of line where varargs or expr start */
    char * param;    /* basic parameter string */
    char * parm2;    /* second parameter string */
    int * list;      /* list of series and/or control integers */
    int * auxlist;   /* needed for "gappy" lag lists */
    char[MAXSAVENAME] savename; /* for object-saving mechanism */
}

bool cmd_ccmt(const(Cmd) * c)
{
    return (c.flags & CmdFlags.ccmt) != 0;
}

bool cmd_subst(const(Cmd) * c)
{
    return (c.flags & CmdFlags.subst) != 0;
}

bool cmd_nosub(const(Cmd) * c)
{
    return (c.flags & CmdFlags.nosub) != 0;
}

alias ExecCallback = extern (C) int function(ExecState *, void *, GretlObjType type);

struct ExecState {
    ExecFlags flags;
    Cmd * cmd;
    Prn * prn;
    char * line;
    char[MAXLEN] runfile;
    Model * model;          /* "workspace" model */
    Model * pmod;           /* set if new model is estimated */
    LoopSet * loop;
    EquationSystem * sys;
    GretlRestriction * rset;
    GretlVar * var;
    void * prev_model;
    GretlObjType prev_type;
    int prev_model_count;
    char * submask;         /* record of incoming sub-sample for functions */
    ubyte padded;           /* record of incoming panel padding, if any */
    ubyte in_comment;
    ubyte free_line;
    ExecCallback callback;
}

struct OpenOp {
    char[FILENAME_MAX] fname;
    int quiet;
    int http;
    int dbdata;
    int ftype;
}

/* structure representing a line of a function or loop */
struct Stmt {
    char * s;      /* text of command line */
    short ci;      /* command index */
    short next;    /* line index to skip to after loop or conditional */
    short idx;     /* 1-based line index allowing for blanks (functions only) */
    short flags;   /* state flags (different for functions vs loops) */
    void * ptr;    /* attachment for LOOPSET or GENERATOR */
}
