/* Minimal working versions of libgretl functions that the ported code needs
   at run time but that have not been ported yet. Everything else that is
   still unported is covered by the aborting stubs in stubs_unported.d.

   These are test scaffolding only: they implement just enough behaviour for
   the smoke tests, and should be deleted (and the real module imported) as
   the corresponding C files are ported. */

module stubs;

import core.stdc.stdio : sprintf;
import core.stdc.stdlib : calloc, malloc;
import core.stdc.string : strcmp;

import core.stdc.string : memset, strncmp;
import gretl.cmd_privateh : Cmd, CmdFlags, ExecState, Stmt;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_prn : Prn, pprintf, pputc;
import gretl.libgretlh;

extern (C):




int printmodel(Model * pmod, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    pprintf(prn, "Model (%d coefficients, %d observations)\n", pmod.ncoeff, pmod.nobs);
    foreach (i; 0 .. pmod.ncoeff) {
        pprintf(prn, "  coeff[%d] = %g (se %g)\n", i, pmod.coeff[i], pmod.sderr[i]);
    }
    return 0;
}

/* function statements (gretl_func_ext.d declares it): a function body that does not recurse */

int statements_get_structure(Stmt * lines, int n_lines, const(char) * fname, int * recurses)
{
    *recurses = 0;
    return 0;
}

/* sub-sampling (subsample.c): the tests never sub-sample the dataset */

void * fetch_full_dataset() { return null; }
int dataset_is_subsampled(const(DatasetRaw) * d) { return 0; }
int dataset_is_resampled(const(DatasetRaw) * d) { return 0; }
int complex_subsampled() { return 0; }
void sync_dataset_shared_members(DatasetRaw * d) {}
int gretl_lists_revise(const(int) * dlist, int dmin) { return 0; }   /* renumbers saved lists after a deletion */

/* matrix_extra.c: vector from a C array, optionally squared elementwise */

import gretl.gretl_matrix : gretl_matrix_alloc;
import gretl.gretl_matrixh;

GretlMatrixRaw * gretl_vector_from_array(const(double) * x, int n, GretlMatrixMod mod)
{
    GretlMatrixRaw * v = gretl_matrix_alloc(n, 1);

    foreach (i; 0 .. n) {
        v.val[i] = (mod == GretlMatrixMod.square) ? x[i] * x[i] : x[i];
    }
    return v;
}

/* variable deletion bookkeeping (dataset deletion checks saved lists etc.) */
int check_variable_deletion_list(const(int) * list, const(DatasetRaw) * dset) { return 0; }

/* no full-dataset bookkeeping in the tests */
void maybe_free_full_dataset(DatasetRaw * dset) { }

/* no function-name tables in the tests: nothing is reserved */
int * get_list_by_name(const(char) * name) { return null; }

/* CSV export settings: defaults (C locale, comma delimiter, "NA") */
const(char) * get_csv_na_write_string() { return "NA"; }
char get_data_export_decpoint() { return '.'; }
char get_data_export_delimiter() { return ','; }
int printing_to_standard_stream(Prn * prn) { return 0; }

/* csv import settings: defaults */
const(char) * get_csv_na_read_string() { return "default"; }
int try_http(const(char) * s, char * fname, int * err) { return 0; }


/* start-up helpers that gretl_set_paths() calls: no-ops in the tests */
void set_gretl_tex_preamble() { }
void set_graph_color_from_string(int i, const(char) * s) { }
void set_plotstyle(const(char) * s) { }
int gretl_is_updated(const(char) * s) { return 0; }
int is_gretl_addon(const(char) * name) { return 0; }
char * gretl_addon_get_path(const(char) * name) { return null; }
void update_addons_index(const(char) ** pkgs) { }

/* printout.c / graphing.c helpers used by describe.c's printing code */
const(char) * plotname(const(DatasetRaw) * dset, int v, int full)
{
    return dset.varname[v];
}

char * gretl_fix_exponent(char * s) { return s; }
int get_gretl_digits() { return 6; }

/* getsymb() has the full table of symbol names in genlex.c: the test version names nothing */


/* printout.c is not ported: simplified stand-in (no exponent fix-up) */
/* nls.c is not ported: no NLS estimation is in progress in the tests */
int doing_nls() { return 0; }





/* objstack.c is not ported: no "last model" bookkeeping in the tests */
import gretl.objstackh : GretlObjType;
void set_as_last_model(void * ptr, GretlObjType type) { }

/* genlex.c is not ported: no constants are read-only in the tests */
int scalar_is_read_only_index(const(char) * name) { return 0; }

/* threading set-up (gretl_openmp/threads code is not ported): single-threaded, nothing to configure */
void num_threads_init(int blas_type) { }

/* clean-up of the saved-object stack and the user-function tables (objstack/function code is not ported): there
   is nothing registered to release */
void gretl_saved_objects_cleanup() { }

/* objstack.c and monte_carlo.c are not ported: the genr lexer and evaluator ask whether a name is a saved object (none is) and whether a loop
   is renaming variables (none is); this used to be repeated in the tests that run genr (the ARMA initialization does) */
extern (C) void * gretl_get_object_by_name(const(char) * name) { return null; }
extern (C) int get_loop_renaming() { return 0; }

/* glibc's mallinfo(), which some tests use to look for leaks, for Windows: only uordblks (bytes in use) is filled in,
   from the busy blocks of the process heap, which is where the C runtime's malloc gets its memory. */
version (Windows) {
    extern (C) struct MallInfoWin {
        int arena, ordblks, smblks, hblks, hblkhd, usmblks, fsmblks, uordblks, fordblks, keepcost;
    }

    extern (C) MallInfoWin mallinfo()
    {
        import core.sys.windows.windows;

        MallInfoWin r;
        PROCESS_HEAP_ENTRY e;
        HANDLE h = GetProcessHeap();
        size_t used = 0;

        HeapLock(h);
        while (HeapWalk(h, &e)) {
            if (e.wFlags & PROCESS_HEAP_ENTRY_BUSY) {
                used += e.cbData;
            }
        }
        HeapUnlock(h);
        r.uordblks = cast(int) used;
        return r;
    }
}

/* definition-time stand-ins for the unported command interpreter (used by gretl_func) */
int gretl_cmd_init(Cmd * cmd)
{
    memset(cmd, 0, Cmd.sizeof);
    return 0;
}

void gretl_cmd_free(Cmd * cmd) {}


/* classify the line by its first word: enough for the structure checks of gretl_function_append_line */
int get_command_index(ExecState * s, int cmode, int preserve)
{
    const(char) * p = s.line;

    s.cmd.flags = cast(CmdFlags) 0;
    while (*p == ' ' || *p == '\t') p++;
    if (*p == '#' || *p == '\0') {
        s.cmd.ci = 0;
        return 0;
    }
    if (!strncmp(p, "end function".ptr, 12)) {
        s.cmd.ci = GretlCmdIndex.end;
        s.cmd.flags = CmdFlags.endfun;
    } else if (!strncmp(p, "return".ptr, 6)) {
        s.cmd.ci = GretlCmdIndex.funcret;
    } else if (!strncmp(p, "if ".ptr, 3)) {
        s.cmd.ci = GretlCmdIndex.if_;
    } else if (!strncmp(p, "else".ptr, 4)) {
        s.cmd.ci = GretlCmdIndex.else_;
    } else if (!strncmp(p, "endif".ptr, 5)) {
        s.cmd.ci = GretlCmdIndex.endif;
    } else if (!strncmp(p, "quit".ptr, 4)) {
        s.cmd.ci = GretlCmdIndex.quit;
    } else {
        s.cmd.ci = GretlCmdIndex.genr;
    }
    return 0;
}
