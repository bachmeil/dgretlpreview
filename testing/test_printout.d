/* Smoke tests for printout.d: the bufgets family (line-wise reading from a string buffer), digits settings,
   the small layout helpers (dashline, bufspace, print_centered), exponent fixing, number printing, print_smpl,
   list_series and printdata on a small dataset.

   Not exercised: frequency distributions and cross-tabulations, covariance-matrix and forecast printing, the
   RTF/TeX variants, and the series-by-column option of printdata with obs markers.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strcmp, strstr;

import gretl.dataset : create_new_dataset;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.printout;

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

bool same(const(char) * a, string b)
{
    import std.string : toStringz;
    return a !is null && strcmp(a, toStringz(b)) == 0;
}

bool has(const(char) * buf, string needle)
{
    import std.string : toStringz;
    return buf !is null && strstr(buf, toStringz(needle)) !is null;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    int err;

    /* bufgets family */
    const(char) * text = "first line\nsecond\nthird line here\n";
    char[64] line;

    check(bufgets_init(text) == 0, "bufgets_init");
    check(bufgets(line.ptr, line.sizeof, text) !is null && same(line.ptr, "first line\n"), "bufgets: first line");
    check(buftell(text) == 11, "buftell after one line");
    check(bufgets(line.ptr, line.sizeof, text) !is null && same(line.ptr, "second\n"), "bufgets: second line");
    check(bufgets_peek_line_length(text) == 16, "peek at the length of the next line");
    check(buf_back_lines(text, 1) == 0 && bufgets(line.ptr, line.sizeof, text) !is null && same(line.ptr, "second\n"),
          "buf_back_lines steps back one line");
    buf_rewind(text);
    check(bufgets(line.ptr, line.sizeof, text) !is null && same(line.ptr, "first line\n"), "buf_rewind");
    check(bufseek(text, 11) == 0 && bufgets(line.ptr, line.sizeof, text) !is null && same(line.ptr, "second\n"), "bufseek");
    bufgets_finalize(text);
    bufgets_cleanup();

    /* digits */
    int d0 = get_gretl_digits();
    set_gretl_digits(4);
    check(get_gretl_digits() == 4, "set/get_gretl_digits");
    set_gretl_digits(d0);

    /* layout helpers */
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    dashline(5, prn);
    check(has(gretl_print_get_buffer(prn), "-----"), "dashline");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    bufspace(3, prn);
    check(same(gretl_print_get_buffer(prn), "   "), "bufspace");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    print_centered("ab", 6, prn);
    check(has(gretl_print_get_buffer(prn), "ab") && strstr(gretl_print_get_buffer(prn), "ab") - gretl_print_get_buffer(prn) == 2,
          "print_centered pads on the left");
    gretl_print_destroy(prn);

    char[16] ex = 0;
    ex[0 .. 8] = "1.5e+007";
    gretl_fix_exponent(ex.ptr);
    check(same(ex.ptr, "1.5e+07"), "gretl_fix_exponent drops a leading zero");

    prn = gretl_print_new(PrnType.buffer, &err);
    gretl_print_value(2.5, prn);
    check(has(gretl_print_get_buffer(prn), "2.5"), "gretl_print_value");
    gretl_print_destroy(prn);

    /* dataset printing */
    enum n = 6;
    DatasetRaw * dset = create_new_dataset(3, n, 0);

    foreach (i, nm; ["const", "x", "y"]) {
        dset.varname[i][0 .. nm.length] = nm[];
        dset.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. n) {
        dset.Z[0][t] = 1;
        dset.Z[1][t] = t + 1.5;
        dset.Z[2][t] = 10.0 * (t + 1);
    }
    dset.t1 = 0;
    dset.t2 = n - 1;

    prn = gretl_print_new(PrnType.buffer, &err);
    list_series(dset, Gretlopt.OPT_NONE, prn);
    check(has(gretl_print_get_buffer(prn), "x") && has(gretl_print_get_buffer(prn), "y"), "list_series names the variables");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    print_smpl(dset, 0, Gretlopt.OPT_NONE, prn);
    check(gretl_print_get_buffer(prn) !is null && strlen2(gretl_print_get_buffer(prn)) > 0, "print_smpl writes something");
    gretl_print_destroy(prn);

    int[3] list = [2, 1, 2];
    prn = gretl_print_new(PrnType.buffer, &err);
    err = printdata(list.ptr, null, dset, Gretlopt.OPT_NONE, prn);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(err == 0 && has(buf, "x") && has(buf, "y") && has(buf, "60") && has(buf, "6.5"), "printdata of two series");
    gretl_print_destroy(prn);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

size_t strlen2(const(char) * s)
{
    import core.stdc.string : strlen;
    return strlen(s);
}
