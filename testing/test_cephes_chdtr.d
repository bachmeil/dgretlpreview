/* Tests for cephes/chdtr.d: the chi-square distribution: chdtr (cdf), chdtrc (upper tail) and chdtri (inverse of the
   upper tail), against closed forms (1 and 2 degrees of freedom), tabulated critical values, the complement
   relation, the inversion round trip (chdtri uses the still-C igami), and the domain errors. A bit-for-bit
   comparison with the original C code over 330 combinations was run once when the module was ported (see
   log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : erf, exp, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.chdtr;

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

bool close(double a, double b, double tol = 1e-12)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* closed forms: df = 2 is exponential, df = 1 is erf(sqrt(x/2)) */
    bool ok = true;
    foreach (x; [0.05, 0.7, 2.0, 6.0, 15.0]) {
        ok = ok && close(chdtr(2, x), 1 - exp(-x / 2)) && close(chdtrc(2, x), exp(-x / 2));
        ok = ok && close(chdtr(1, x), erf(sqrt(x / 2)));
    }
    check(ok, "df = 1 and df = 2 match their closed forms");

    /* tabulated 5% critical values */
    check(fabs(chdtri(1, 0.05) - 3.8415) < 1e-4 && fabs(chdtri(5, 0.05) - 11.0705) < 1e-4 && fabs(chdtri(10, 0.05) - 18.3070) < 1e-4 && fabs(chdtri(30, 0.05) - 43.7730) < 1e-4, "5% critical values for 1, 5, 10 and 30 degrees of freedom");
    check(fabs(chdtrc(10, 18.307) - 0.05) < 1e-5, "chdtrc at a critical value");

    /* complement */
    check(close(chdtr(7, 5.5) + chdtrc(7, 5.5), 1.0, 1e-14) && chdtr(4, 0.0) == 0.0 && chdtrc(4, 0.0) == 1.0, "cdf + upper tail = 1; limits at 0");

    /* round trip */
    ok = true;
    foreach (df; [1.0, 2.0, 3.5, 10.0, 50.0, 400.0]) {
        foreach (y; [1e-6, 0.001, 0.05, 0.5, 0.95, 0.999]) {
            double x = chdtri(df, y);
            ok = ok && fabs(chdtrc(df, x) - y) <= 1e-9 * y + 1e-14;
        }
    }
    check(ok, "chdtrc(df, chdtri(df, y)) = y");

    /* domain errors */
    check(chdtr(0.5, 1.0) == 0.0 && chdtr(3, -1.0) == 0.0 && chdtrc(0.5, 1.0) == 0.0 && chdtrc(3, -1.0) == 0.0, "df below 1 and negative x give 0");
    check(chdtri(3, -0.1) == 0.0 && chdtri(3, 1.1) == 0.0 && chdtri(0.5, 0.5) == 0.0, "chdtri with y outside [0, 1] or df below 1 gives 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
