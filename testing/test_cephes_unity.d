/* Tests for cephes/unity.d: log(1 + x), exp(x) - 1 and cos(x) - 1 computed accurately near zero. cephes_log and
   cephes_exp are compared with libm's log1p and expm1, and cosm1 with -2 sin^2(x/2), in relative error over
   arguments from 1e-300 to a few units, including the ranges where the naive formulas lose all digits. A
   bit-for-bit comparison with the original C code over 205,218 arguments was run once when the module was ported
   (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : expm1, fabs, isnan, log1p, sin;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.unity;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tiny arguments: the first-order terms, where log(1 + x) and exp(x) - 1 in plain arithmetic give 0 or noise */
    check(cephes_log(1e-20) == 1e-20 && cephes_exp(1e-20) == 1e-20 && cephes_log(-1e-20) == -1e-20, "tiny x: log(1+x) = exp(x) - 1 = x");
    check(cosm1(1e-9) == -0.5e-18 && cosm1(0.0) == 0.0, "cosm1 for tiny x is -x^2/2 (plain cos(x) - 1 would give 0)");

    /* against libm over many magnitudes, both signs */
    double wl = 0, we = 0, wc = 0;
    for (double x = 1e-300; x < 4.0; x *= 1.31) {
        foreach (sgn; [1.0, -1.0]) {
            double y = sgn * x;
            if (y > -1.0) {
                double e = fabs(cephes_log(y) - log1p(y)) / fabs(log1p(y));
                if (e > wl) wl = e;
            }
            double e = fabs(cephes_exp(y) - expm1(y)) / fabs(expm1(y));
            if (e > we) we = e;
            if (x > 1e-100) {   // below this x^2 underflows to the denormals, where relative error is meaningless
                double r = -2.0 * sin(y / 2) * sin(y / 2);
                e = fabs(cosm1(y) - r) / fabs(r);
                if (e > wc) wc = e;
            }
        }
    }
    printf("  worst relative errors: log %g, exp %g, cosm1 %g\n", wl, we, wc);
    check(wl < 1e-15 && we < 1e-15 && wc < 1e-14, "accurate to a few units in the last place against libm over 1e-300 .. 4 (cosm1 from 1e-100)");

    /* branches: just inside and outside the rational ranges */
    check(fabs(cephes_log(0.4142) - log1p(0.4142)) < 1e-15 && fabs(cephes_log(0.42) - log1p(0.42)) < 1e-15 && fabs(cephes_exp(0.49) - expm1(0.49)) < 1e-15 && fabs(cephes_exp(0.51) - expm1(0.51)) < 1e-15, "continuity across the switch points");

    /* special values */
    check(cephes_exp(double.infinity) == double.infinity && cephes_exp(-double.infinity) == -1.0 && isnan(cephes_exp(double.nan)) != 0, "exp(x) - 1 at infinities and NaN");
    check(cephes_log(-1.0) == -double.infinity, "log(1 + x) at x = -1 is minus infinity");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
