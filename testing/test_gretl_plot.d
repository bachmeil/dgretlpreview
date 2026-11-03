/* Smoke tests for gretl_plot.d: the "plot" block state machine: starting a plot with a series as data source,
   appending literal lines, the error returns for an unknown source, an unknown field and a duplicated start,
   and plot_command_active.

   Not exercised: gretl_plot_finalize (it builds and runs a gnuplot plot), "option" lines (they are checked
   against the full option table), printf/strings fields and matrix or list data sources (user variables).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp;

import gretl.dataset : create_new_dataset;
import gretl.gretl_plot;
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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    DatasetRaw * d = create_new_dataset(2, 5, 0);

    d.varname[1][0 .. 1] = "x";
    d.varname[1][1] = 0;

    check(plot_command_active() == 0, "no plot block is active at the start");
    check(gretl_plot_append_line("literal set title 'a'", d) != 0, "appending without a started plot is an error");

    check(gretl_plot_start("nosuchthing", d) != 0 && plot_command_active() == 0, "an unknown data source is an error and starts nothing");

    check(gretl_plot_start("x", d) == 0 && plot_command_active() == 1, "starting a plot with the series x");
    check(plot.datatype == GretlType.series && plot.datasource !is null && strcmp(plot.datasource, "x") == 0,
          "the data source is recorded as a series");

    check(gretl_plot_append_line("literal set title 'My plot'", d) == 0 && plot.nlines == 1, "a literal line is collected");
    check(plot.lines[0] !is null && strcmp(plot.lines[0], "set title 'My plot'") == 0, "the literal text is kept without the field name");
    check(gretl_plot_append_line("literal set xlabel 'X'", d) == 0 && plot.nlines == 2, "a second literal line");

    check(gretl_plot_start("x", d) != 0 && plot_command_active() == 0, "starting while a plot is in progress is an error and clears it");

    check(gretl_plot_start("x", d) == 0, "a new plot can be started after the error");
    check(gretl_plot_append_line("bogus field", d) != 0 && plot_command_active() == 0, "an unknown field is an error and clears the plot");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
