/* Port of interact.h */

module gretl.interacth;


/* ExecState is defined in cmd_private.h (cmd_privateh.d); it is re-exported from here for the modules that
   imported the former opaque declaration */
public import gretl.cmd_privateh : ExecState;

import gretl.gretl_commandsh;

enum MAXSAVENAME = 32;
enum CMD_NULL = -1;
enum CMD_COMMENT = -2;
enum CMD_MASKED = -3;

enum ProgramOptions {
    batch   = 1 << 0,
    help    = 1 << 1,
    version_ = 1 << 2,
    runit   = 1 << 3,
    dbopen  = 1 << 4,
    webdb   = 1 << 5,
    fnpkg   = 1 << 6,
    dump    = 1 << 7,
    debug_  = 1 << 8,
    quiet   = 1 << 9,
    english = 1 << 10,
    basque  = 1 << 11,
    makepkg = 1 << 12,
    instpkg = 1 << 13,
    tool    = 1 << 14,
    no_plot = 1 << 15
}

enum ExecFlags {
    console_exec  = 1 << 0,
    script_exec   = 1 << 1,
    include_exec  = 1 << 2,
    function_exec = 1 << 3,
    debug_exec    = 1 << 4,
    callback_exec = 1 << 5,
    init_exec     = 1 << 6
}

bool HIDDEN_COMMAND(GretlCmdIndex c) {
    return (c == GretlCmdIndex.funcerr || c == GretlCmdIndex.funcret);
}
