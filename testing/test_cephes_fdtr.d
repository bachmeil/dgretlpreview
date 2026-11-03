/* Tests for cephes/fdtr.d: the F distribution: cdf, upper tail and inverse, against closed forms (F(1,1) and
   F(2,2)), tabulated critical values, complement, inversion round trips (the inverse uses the still-C incbi),
   and the domain errors. A bit-for-bit comparison with the original C code over 3,630 combinations was run once
   when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : atan, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.fdtr;

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

    /* closed forms: F(2,2) cdf is x/(1+x), F(1,1) cdf is (2/pi) atan(sqrt x) */
    bool ok = true;
    foreach (x; [0.01, 0.3, 1.0, 4.0, 50.0]) {
        ok = ok && close(fdtr(2, 2, x), x / (1 + x)) && close(fdtr(1, 1, x), 2 / PI * atan(sqrt(x)), 1e-11);
        ok = ok && close(fdtrc(2, 2, x), 1 / (1 + x)) && close(fdtrc(1, 1, x), 1 - 2 / PI * atan(sqrt(x)), 1e-11);
    }
    check(ok, "F(2,2) and F(1,1) cdf and upper tail match their closed forms");

    check(close(fdtr(10, 10, 1.0), 0.5) && fdtr(3, 7, 0.0) == 0.0, "F(10,10) has median 1; cdf at 0 is 0");
    check(close(fdtr(4, 9, 2.7) + fdtrc(4, 9, 2.7), 1.0, 1e-14), "cdf + upper tail = 1");

    /* tabulated 5% critical values */
    check(fabs(fdtri(1, 1, 0.05) - 161.4476) < 1e-3 && fabs(fdtri(5, 5, 0.05) - 5.0503) < 1e-4 && fabs(fdtri(4, 10, 0.05) - 3.4780) < 1e-4, "critical values F(0.95; 1,1), (5,5) and (4,10)");

    /* round trips: fdtrc(a, b, fdtri(a, b, y)) = y, across the code's two solution branches */
    ok = true;
    foreach (a; [1.0, 3.0, 12.0, 100.0]) {
        foreach (b; [2.0, 7.0, 30.0, 500.0]) {
            foreach (y; [1e-5, 0.0005, 0.01, 0.2, 0.5, 0.9, 0.999]) {
                ok = ok && fabs(fdtrc(a, b, fdtri(a, b, y)) - y) <= 1e-9 * y + 1e-14;
            }
        }
    }
    check(ok, "fdtrc(fdtri(y)) = y");

    /* domain errors */
    check(fdtr(0.5, 3, 1.0) == 0.0 && fdtr(2, 0.9, 1.0) == 0.0 && fdtr(2, 3, -1.0) == 0.0 && fdtrc(0.5, 3, 1.0) == 0.0 && fdtrc(2, 3, -1.0) == 0.0, "degrees of freedom below 1 and negative x give 0");
    check(fdtri(2, 3, 0.0) == 0.0 && fdtri(2, 3, 1.5) == 0.0 && fdtri(0.5, 3, 0.5) == 0.0, "fdtri with y outside (0, 1] or degrees of freedom below 1 gives 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
