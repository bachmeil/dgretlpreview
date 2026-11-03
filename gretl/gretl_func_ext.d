/* Declarations of the functions that gretl_func.c (gretl_func0 .. gretl_func4.d) calls and that are not ported yet: the command-line / execution
   state (interact.c, tokenizer), the flow control and loop subsystems (flow_control.c, monte_carlo.c), the markdown converter (gretl_mdconv.c),
   the foreign-language and add-on helpers (gretl_foreign.c, addons_utils.c), the function-statement analysis. The signatures are those of the C headers where those are in the tree (interact.h) and are otherwise taken from the call
   sites in gretl_func.c. They are extern (C), so they bind to the aborting stubs of testing/stubs_unported.d until the files are ported. */

module gretl.gretl_func_ext;

import gretl.cmd_privateh : Cmd, ExecCallback, ExecState, Stmt;
import gretl.gretl_prn : Prn;
import gretl.libgretlh : DatasetRaw, Model;
import gretl.glist : GList;

extern (C) {
    /* interact.c */
    int get_command_index(ExecState * s, int cmode, int preserve);
    int gretl_cmd_init(Cmd * cmd);
    void gretl_cmd_free(Cmd * cmd);
    int parse_command_line(ExecState * s, DatasetRaw * dset, void * ptr);
    int is_plotting_command(Cmd * cmd);
    int maybe_exec_line(ExecState * s, DatasetRaw * dset, void * ptr);
    int flow_control(ExecState * s, DatasetRaw * dset, void * pgenr);
    void gretl_exec_state_init(ExecState * s, int flags, char * line, Cmd * cmd, Model * model, Prn * prn);
    void gretl_exec_state_clear(ExecState * s);
    ExecCallback get_gui_callback();

    /* flow_control.c */
    int gretl_if_state_check(int indent0);
    int gretl_if_state_false();
    int gretl_if_state_record();
    void gretl_if_state_reset(int indent0);

    /* monte_carlo.c */
    void gretl_loop_destroy(void * loop);
    int gretl_loop_exec(ExecState * s, DatasetRaw * dset);
    void loop_reset_uvars(void * loop);
    int get_loop_renaming();

    /* gretl_mdconv.c */
    void md_to_gretl(const(char) * buf, Prn * prn);
    int simple_markdown_detect(const(char) * buf);

    /* gretl_foreign.c, addons_utils.c, function statements (function_lookup and the function overrides are in genlex.d) */
    int check_R_depends(const(char) * pkgname, const(char) * Rdeps, Prn * prn);
    int is_gretl_addon(const(char) * pkgname);
    int statements_get_structure(Stmt * lines, int n_lines, const(char) * fname, int * recurses);
}
