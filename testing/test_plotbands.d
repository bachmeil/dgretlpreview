/* Smoke tests for plotbands.d: the band-style names, band_info_new/free, the matrix comparison used to
   order bands, series creation for constant band widths, and the gnuplot text written for band
   lines, bars and titles.

   Not exercised: plot_with_band itself (needs the gnuplot support functions of graphing.c), the bundle-based
   band specifications and the legacy option parsing.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : FILE, _IOLBF, fclose, fgets, printf, rewind, setvbuf, stdout, tmpfile;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strlen;

import gretl.dataset : create_new_dataset;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.plotbands;
import gretl.plot_privh : GnuplotInfo;

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

    int err = 0;

    check(style_from_string("line", &err) == BandStyle.line && err == 0, "style 'line'");
    check(style_from_string("fill", &err) == BandStyle.fill && style_from_string("dash", &err) == BandStyle.dash
          && style_from_string("bars", &err) == BandStyle.bars && style_from_string("step", &err) == BandStyle.step && err == 0,
          "styles fill, dash, bars and step");
    style_from_string("bogus", &err);
    check(err != 0, "an unknown style is an error");

    BandInfo * bi = band_info_new();

    check(bi !is null && bi.factor == 1.0 && bi.style == BandStyle.line && bi.center == -1 && bi.width == -1, "band_info_new defaults: factor 1, line style, no series");
    free(bi);

    GretlMatrixRaw * a = gretl_matrix_alloc(2, 2);
    GretlMatrixRaw * b = gretl_matrix_alloc(2, 2);

    foreach (i; 0 .. 4) {
        a.val[i] = b.val[i] = i + 1.0;
    }
    check(band_matrices_differ(a, b) == 0, "equal matrices do not differ");
    b.val[3] = 9.0;
    check(band_matrices_differ(a, b) == 1, "a changed element is detected");
    gretl_matrix_free(a);
    gretl_matrix_free(b);

    /* constant series for a scalar band width: returns the series number, or 0 (the constant) for a width of 1 */
    DatasetRaw * d = create_new_dataset(3, 5, 0);
    int v_next = 1;

    d.t1 = 0;
    d.t2 = 4;
    int v = make_const_series(0.25, d, &v_next);
    check(v == 1 && v_next == 2 && d.Z[1][0] == 0.25 && d.Z[1][4] == 0.25, "make_const_series fills a new series with the constant");
    check(make_const_series(1.0, d, &v_next) == 0 && v_next == 2, "a width of 1 uses the constant series and adds nothing");

    /* gnuplot text */
    FILE * fp = tmpfile();
    gp_newline(0, fp);
    gp_newline(1, fp);
    string s = slurp(fp);
    check(s == "\n, \\\n", "gp_newline ends a plot element with a newline or a continuation");
    fclose(fp);

    BandInfo bb;
    bb.center = 2;
    bb.width = 3;
    bb.factor = 1.0;
    bb.style = BandStyle.line;
    bb.rgb[0] = 0;
    bb.title = "my band";

    fp = tmpfile();
    band_title(&bb, fp);
    check(slurp(fp) == "title 'my band' ", "band title");
    fclose(fp);

    bb.title = null;
    fp = tmpfile();
    band_title(&bb, fp);
    check(slurp(fp) == "notitle ", "no title");
    fclose(fp);

    fp = tmpfile();
    print_pm_lines(&bb, 1, "", 0, fp);
    s = slurp(fp);
    check(contains(s, "$data using 1:($2-1*$3)") && contains(s, "w lines") && contains(s, "($2+1*$3)"), "plus/minus lines use center -/+ factor * width");
    fclose(fp);

    bb.factor = 2.0;
    fp = tmpfile();
    print_pm_lines(&bb, 1, "", 0, fp);
    check(contains(slurp(fp), "2*$3"), "the multiplicative factor appears in the lines");
    fclose(fp);

    bb.style = BandStyle.step;
    fp = tmpfile();
    print_pm_lines(&bb, 1, "", 0, fp);
    check(contains(slurp(fp), "steps"), "step style gives steps");
    fclose(fp);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
