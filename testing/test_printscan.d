/* Smoke tests for printscan.d: sprintf/printf with literal numeric arguments (do_sprintf_function, do_printf,
   do_printf_command), the format-chunk scanner and the error returns for malformed formats.

   Not exercised: arguments that need the genr parser (string literals, expressions, series, matrices, named scalars), sscanf
   (which stores into user variables) and obs-marker generation.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp;

import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.printscan;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    int err = 0;
    char * s = do_sprintf_function("%d items", "3", null, &err);

    check(err == 0 && same(s, "3 items"), "sprintf with an integer literal");
    free(s);

    s = do_sprintf_function("%.2f and %5.1f|", "2.5, 3.25", null, &err);
    check(err == 0 && same(s, "2.50 and   3.2|") || same(s, "2.50 and   3.3|"), "sprintf with two doubles");
    free(s);

    s = do_sprintf_function("100%%", "", null, &err);
    check(err == 0 && same(s, "100%"), "sprintf of an escaped percent");
    free(s);

    s = do_sprintf_function("%d %d", "1", null, &err);
    check(err != 0, "too few arguments is an error");
    err = 0;

    s = do_sprintf_function("%y", "1", null, &err);
    check(err != 0, "bad conversion character is an error");
    err = 0;

    int nchars = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);
    err = do_printf("x=%g\n", "4.5", null, prn, &nchars);
    const(char) * buf = gretl_print_get_buffer(prn);
    check(err == 0 && same(buf, "x=4.5\n") && nchars == 6, "printf to a buffer");
    gretl_print_destroy(prn);

    prn = gretl_print_new(PrnType.buffer, &err);
    err = do_printf_command("%d-%d", "7, 8", null, prn, 1);
    buf = gretl_print_get_buffer(prn);
    check(err == 0 && same(buf, "7-8"), "printf command with a quoted format");
    gretl_print_destroy(prn);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
