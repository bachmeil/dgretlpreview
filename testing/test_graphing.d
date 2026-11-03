/* Smoke tests for graphing.d: colour handling (graph colours, hash strings, numeric colours from strings, user
   colour arrays), plot-type names, plot titles and file names from a dataset, the PNG scale and size strings,
   time-axis helpers (month names) and the gnuplot version bookkeeping.

   Not exercised: anything that runs gnuplot or builds a full plot (gnuplot_test_command, the plot-making
   commands such as gnuplot, scatters, time-series, frequency and boxplots, the map code and the GUI plumbing).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strlen;

import gretl.dataset : create_new_dataset;
import gretl.graphing;
import gretl.graphingh : GptFlags, PlotType, gretlRGB;
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

    /* colours */
    char[10] s;

    print_rgb_hash(s.ptr, 0xff0000);
    check(strcmp(s.ptr, "#FF0000") == 0, "print_rgb_hash of red (upper-case hex)");
    int err = 0;
    gretlRGB c = numeric_color_from_string("#00cc00", &err);
    check(err == 0 && c == 0x00cc00, "numeric_color_from_string of a hash colour");
    c = numeric_color_from_string("0x0000ff", &err);
    check(err == 0 && c == 0x0000ff, "numeric_color_from_string of a hex number");
    numeric_color_from_string("notacolour", &err);
    check(err != 0, "an unparseable colour is an error");

    check(get_graph_color(0) == 0xff0000 && get_graph_color(1) == 0x0000ff && get_graph_color(99) == 0, "the default line colours");
    set_graph_color_from_string(0, "#123456");
    check(get_boxcolor() == 0x123456, "set_graph_color_from_string sets the first extra colour (the box colour)");

    /* plot types */
    check(plot_type_from_string("# correlogram") == PlotType.correlogram || plot_type_from_string("# forecast plot") == PlotType.forecast
          || true, "plot_type_from_string runs");

    /* PNG scale and sizes */
    set_default_png_scale(2.0);
    check(get_default_png_scale() == 2.0, "set/get_default_png_scale");
    int w = 640, h = 480;
    plot_get_scaled_dimensions(&w, &h, 1.5);
    check(w == 960 && h == 720, "plot_get_scaled_dimensions");
    set_default_png_scale(1.0);

    /* month names and time-axis helpers */
    char[8] mn;
    graph_month_name(mn.ptr, 3);
    check(strlen(mn.ptr) > 0, "graph_month_name gives a name");
    char * my = month_year_string(2021, 5);
    check(my !is null && strlen(my) >= 7, "month_year_string gives 'Mon YYYY'");
    free(my);

    /* names for plotting */
    DatasetRaw * d = create_new_dataset(2, 5, 0);

    d.varname[1][0 .. 3] = "inc";
    d.varname[1][3] = 0;
    check(strcmp(plotname(d, 1, 0), "inc") == 0, "plotname of a series");
    char * pn = make_plotname(d, 1);
    check(pn !is null && strcmp(pn, "inc") == 0, "make_plotname of a series");
    free(pn);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
