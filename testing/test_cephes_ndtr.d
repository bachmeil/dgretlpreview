/* Tests for cephes/ndtr.d: the normal distribution function at tabulated points, symmetry, the tails (against
   libm's erfc), underflow to 0 and rounding to 1, and special values. A bit-for-bit comparison with the original
   C code over 42,658 arguments was run once when the module was ported (see log.md). The inverse (ndtri) is used to
   round-trip.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : erfc, fabs, isnan, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.ndtr;
import cephes.ndtri : ndtri;

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

bool close(double a, double b, double tol = 1e-13)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tabulated values */
    check(ndtr(0.0) == 0.5, "Phi(0) = 1/2");
    check(close(ndtr(1.0), 0.8413447460685429) && close(ndtr(-1.0), 0.15865525393145707), "Phi(+-1)");
    check(close(ndtr(1.96), 0.9750021048517795) && close(ndtr(-1.96), 0.024997895148220435), "Phi(+-1.96)");
    check(close(ndtr(3.0), 0.9986501019683699) && close(ndtr(-3.0), 0.0013498980316300946), "Phi(+-3)");

    /* against libm: Phi(x) = erfc(-x/sqrt 2)/2, over the whole range, in relative error (including the far tails) */
    double worst = 0;
    for (double x = -37.0; x < 8.0; x += 0.0731) {
        double r = 0.5 * erfc(-x / sqrt(2.0));
        double e = fabs(ndtr(x) - r) / r;
        if (e > worst) worst = e;
    }
    printf("  worst relative error against libm over -37..8: %g\n", worst);
    check(worst < 1e-11, "agrees with libm erfc over the lower tail and the body");

    /* symmetry, away from the upper tail where 1 - small loses digits */
    bool ok = true;
    for (double x = -6.0; x < 6.0; x += 0.113) {
        ok = ok && fabs(ndtr(x) + ndtr(-x) - 1.0) < 1e-15;
    }
    check(ok, "Phi(x) + Phi(-x) = 1");

    /* extreme arguments */
    check(ndtr(-40.0) == 0.0 && ndtr(40.0) == 1.0 && ndtr(-1e5) == 0.0 && ndtr(1e5) == 1.0, "underflow to 0 and rounding to 1");
    check(ndtr(double.infinity) == 1.0 && ndtr(-double.infinity) == 0.0 && isnan(ndtr(double.nan)), "infinities and NaN");

    /* round trip with the (C) inverse */
    ok = true;
    foreach (p; [1e-12, 1e-6, 0.01, 0.2, 0.5, 0.9, 0.999]) {
        ok = ok && fabs(ndtr(ndtri(p)) - p) <= 1e-12 * p;
    }
    check(ok, "ndtr(ndtri(p)) = p");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
