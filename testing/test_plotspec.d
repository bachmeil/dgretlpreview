/* Smoke tests for plotspec.d: creating and destroying a plot specification, the line-style name tables, key
   position strings, adding/deleting/cloning lines, labels and arrows, quote escaping, writing labels, arrows and
   user-defined lines to a stream, and the "recession bars" container.

   Not exercised: plotspec_print (needs the gnuplot support functions of graphing.c), fitted-line and loess
   handling, and the parsing of a bars file.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : FILE, _IOLBF, fclose, fgets, printf, rewind, setvbuf, stdout, tmpfile;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strlen, strstr;

import gretl.libgretlh;
import gretl.plotspec;
import gretl.plotspech : GpKeySpec, GpLineStyle, GptArrow, GptLabel, GptLine, GptSpec, PlotBars;

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

bool contains(string hay, string needle)
{
    import std.algorithm.searching : canFind;
    return hay.canFind(needle);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    GptSpec * spec = plotspec_new();

    check(spec !is null && spec.n_lines == 0 && spec.n_labels == 0 && spec.n_arrows == 0, "plotspec_new starts empty");

    /* style tables */
    check(strcmp(gp_line_style_display_name(GpLineStyle.linespoints), "lines/points") == 0, "display name of linespoints");
    check(gp_style_index_from_name("impulses") == GpLineStyle.impulses && gp_style_index_from_name("lp") == GpLineStyle.linespoints,
          "style index from names and abbreviations");
    check(gp_style_index_from_name("bogus") == GpLineStyle.auto_, "unknown style falls back to auto");
    check(gp_style_index_from_display_name("error bars") == GpLineStyle.errorbars, "style index from display name");

    /* key position strings */
    char * ks = gp_keyspec_string(cast(GpKeySpec) (GpKeySpec.left | GpKeySpec.top));
    check(ks !is null && strcmp(ks, " left top") == 0, "key spec string for left top");
    free(ks);
    check(gp_keyspec_from_string("left top") == cast(GpKeySpec) (GpKeySpec.left | GpKeySpec.top), "key spec parsed from a string");

    /* lines */
    check(plotspec_add_line(spec) == 0 && plotspec_add_line(spec) == 0 && spec.n_lines == 2, "two lines added");
    spec.lines[0].style = GpLineStyle.lines;
    spec.lines[1].style = GpLineStyle.points;
    GptLine * cl = plotspec_clone_lines(spec, null);
    int err = 0;
    cl = plotspec_clone_lines(spec, &err);
    check(err == 0 && cl !is null && cl[1].style == GpLineStyle.points, "lines are cloned");
    free(cl);
    check(plotspec_delete_line(spec, 0) == 0 && spec.n_lines == 1 && spec.lines[0].style == GpLineStyle.points, "deleting a line shifts the rest down");

    /* labels and arrows */
    check(plotspec_add_label(spec) == 0 && spec.n_labels == 1, "label added");
    spec.labels[0].text[0 .. 5] = "a\"b\0\0"[0 .. 5];
    spec.labels[0].pos[0] = 1.5;
    spec.labels[0].pos[1] = 2.5;
    check(plotspec_add_arrow(spec) == 0 && spec.n_arrows == 1, "arrow added");
    spec.arrows[0].x0 = 0;
    spec.arrows[0].y0 = 0;
    spec.arrows[0].x1 = 1;
    spec.arrows[0].y1 = 1;

    char * esc = escape_quotes("say \"hi\"");
    check(esc !is null && strcmp(esc, "say \\\"hi\\\"") == 0, "quotes are escaped with backslashes");
    free(esc);

    FILE * fp = tmpfile();
    print_plot_labelspec(&spec.labels[0], fp);
    string lab = slurp(fp);
    check(contains(lab, "set label") && contains(lab, "1.5") && contains(lab, "2.5"), "label written as a gnuplot command");
    fclose(fp);

    fp = tmpfile();
    print_plot_arrow(&spec.arrows[0], fp);
    check(contains(slurp(fp), "set arrow"), "arrow written as a gnuplot command");
    fclose(fp);

    check(plotspec_delete_label(spec, 0) == 0 && spec.n_labels == 0, "label deleted");
    check(plotspec_delete_arrow(spec, 0) == 0 && spec.n_arrows == 0, "arrow deleted");

    /* bars */
    spec.nbars = 3;
    check(plotspec_allocate_bars(spec) == 0 && spec.bars !is null && spec.bars.n == 3, "bars allocated for three start-stop pairs");
    check(plotspec_set_bar_info(spec, 0, 10.0, 20.0) == 0 && spec.bars.dx[0][0] == 10.0 && spec.bars.dx[0][1] == 20.0, "bar info set");
    plotspec_remove_bars(spec);
    check(spec.bars is null, "bars removed");

    plotspec_destroy(spec);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
