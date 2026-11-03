/* Smoke tests for modelprint.d: an OLS model is printed through printmodel() in plain text, RTF and CSV
   formats, and the buffers are checked for the coefficient names, the statistics headings and the format markers.
   Also tests the small helpers char_len and glyph_count, and the model predicates.

   Not exercised: TeX output (needs the unported tex_* helpers of printout.c), the printing of NLS/ARMA/GARCH/panel/discrete-choice/VAR-specific footers, model tables and
   the "alternate" (--alt) layouts.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strstr;

import gretl.dataset : create_new_dataset;
import gretl.estimate : lsq;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_model : clear_model;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new, gretl_print_set_format;
import gretl.gretl_prnh : PrnFormat, PrnType;
import gretl.libgretlh;
import gretl.modelprint;

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

bool has(const(char) * buf, string needle)
{
    import std.string : toStringz;
    return buf !is null && strstr(buf, toStringz(needle)) !is null;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum n = 14;
    immutable double[n] noise = [0.3, -0.2, 0.5, -0.6, 0.1, 0.4, -0.3, -0.2, -0.5, 0.2, 0.1, -0.1, 0.0, 0.2];
    DatasetRaw * d = create_new_dataset(3, n, 0);

    foreach (i, nm; ["const", "y", "x"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. n) {
        d.Z[0][t] = 1;
        d.Z[2][t] = t + 1.0;
        d.Z[1][t] = 1.0 + 2.0 * (t + 1.0) + noise[t];
    }
    int[4] list = [3, 1, 0, 2];
    Model m = lsq(list.ptr, d, GretlCmdIndex.ols, Gretlopt.OPT_NONE);

    check(m.errcode == 0 && m.ncoeff == 2, "OLS model estimated");

    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    check(err == 0 && prn !is null, "buffer Prn created");
    printmodel(&m, d, Gretlopt.OPT_NONE, prn);

    const(char) * buf = gretl_print_get_buffer(prn);

    check(has(buf, "Model 1") || has(buf, "OLS"), "plain text: model heading");
    check(has(buf, "coefficient") && has(buf, "std. error"), "plain text: coefficient table headings");
    check(has(buf, "const") && has(buf, " x"), "plain text: coefficient names");
    check(has(buf, "R-squared") && has(buf, "Sum squared resid"), "plain text: summary statistics");
    check(has(buf, "Akaike criterion"), "plain text: information criteria");
    gretl_print_destroy(prn);

    /* RTF */
    prn = gretl_print_new(PrnType.buffer, &err);
    gretl_print_set_format(prn, PrnFormat.rtf);
    printmodel(&m, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(has(buf, "\\par") || has(buf, "\\trowd") || has(buf, "\\intbl"), "RTF: control words present");
    gretl_print_destroy(prn);

    /* CSV */
    prn = gretl_print_new(PrnType.buffer, &err);
    gretl_print_set_format(prn, PrnFormat.csv);
    printmodel(&m, d, Gretlopt.OPT_NONE, prn);
    buf = gretl_print_get_buffer(prn);
    check(has(buf, "const") && has(buf, ","), "CSV: coefficient name and separators present");
    gretl_print_destroy(prn);

    /* quiet option prints nothing */
    prn = gretl_print_new(PrnType.buffer, &err);
    printmodel(&m, d, Gretlopt.OPT_Q, prn);
    buf = gretl_print_get_buffer(prn);
    check(buf is null || buf[0] == 0, "OPT_Q suppresses output");
    gretl_print_destroy(prn);

    /* helpers */
    check(char_len("abc") == 3, "char_len of ASCII");
    check(glyph_count("abc") == 3, "glyph_count of ASCII");
    check(glyph_count("\u4e2d\u6587") == 4, "glyph_count counts wide characters twice");
    check(binary_model(&m) == false && ordered_model(&m) == false, "OLS is not a binary or ordered model");
    check(XDIGITS(&m) > 0 && FDIGITS(&m) == 5, "digit settings for a non-MP model");

    clear_model(&m);
    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
