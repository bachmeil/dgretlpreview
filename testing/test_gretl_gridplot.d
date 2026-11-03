/* Smoke tests for gretl_gridplot.d: the grid-layout arithmetic (set_mp_grid), the graphic-file-name test, the
   filter that drops "set term" lines from a sub-plot, the matrix index helper, option checking and the
   collecting state.

   Not exercised: collecting plots into a strings array, finalizing a multiplot script and building one from an
   array (they need the array/user-variable machinery and gnuplot), and the file-name search.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : FILE, _IOLBF, fclose, fgets, printf, rewind, setvbuf, stdout, tmpfile;
import core.stdc.string : strlen;

import gretl.gretl_gridplot;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

bool close(double a, double b, double tol = 1e-8)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

string slurp(FILE * fp)
{
    char[256] line;
    string out_;

    rewind(fp);
    while (fgets(line.ptr, line.length, fp) !is null) {
        out_ ~= line[0 .. strlen(line.ptr)].idup;
    }
    return out_;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    check(gretl_gridplot_collecting() == 0, "nothing is being collected at the start");
    check(gretl_gridplot_add_plot(cast(char *) "plot x".ptr) != 0, "adding a plot outside a gpbuild block is an error");

    /* automatic grid for 5 plots: 3 rows (ceil sqrt 5), 2 columns */
    mp_rows = 0;
    mp_cols = 0;
    check(set_mp_grid(5) == 0 && mp_rows == 3 && mp_cols == 2, "fully automatic grid for 5 plots is 3 x 2");
    mp_rows = 2;
    mp_cols = 0;
    check(set_mp_grid(5) == 0 && mp_cols == 3, "with two rows requested, five plots need three columns");
    mp_rows = 0;
    mp_cols = 4;
    check(set_mp_grid(9) == 0 && mp_rows == 3, "with four columns requested, nine plots need three rows");
    mp_rows = 7;
    mp_cols = 0;
    check(set_mp_grid(5) != 0, "more rows than plots is an error");
    mp_rows = 0;
    mp_cols = 0;

    /* graphic file names */
    check(is_graphic_filename("fig.png") == 1 && is_graphic_filename("a/b.pdf") == 1 && is_graphic_filename("data.csv") == 0,
          "is_graphic_filename");

    /* sub-plot filtering */
    FILE * fp = tmpfile();
    filter_subplot("set term png\nplot x\nset term pdf\nplot y\n", fp);
    check(slurp(fp) == "plot x\nplot y\n", "filter_subplot drops the set term lines");
    fclose(fp);

    /* matrix index helper: entries are 1-based plot numbers */
    GretlMatrixRaw * m = gretl_matrix_alloc(2, 2);

    m.val[0] = 1;
    m.val[1] = 3;
    m.val[2] = 2;
    m.val[3] = 4;
    check(get_subplot_index(m, 0, 1) == 1 && get_subplot_index(m, 1, 0) == 2, "get_subplot_index returns a zero-based index");
    gretl_matrix_free(m);

    /* defaults and option checking */
    mp_fontsize = 99;
    set_multiplot_defaults();
    check(mp_fontsize == 10 && mp_width == 800 && mp_height == 600, "set_multiplot_defaults");
    check(check_gridplot_options(Gretlopt.OPT_NONE) == 0, "no options is fine");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
