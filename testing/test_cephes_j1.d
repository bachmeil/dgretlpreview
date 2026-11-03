/* Tests for cephes/j1.d: the Bessel functions of order one, J1 and Y1, against libm's j1/y1 over the regimes
   (x <= 5 and the asymptotic range), tabulated values, oddness of J1 for |x| <= 5, the first zero, and the domain
   error of Y1. The C returns NaN for x < -5 (bugs.md #105); that behaviour is recorded here. A bit-for-bit
   comparison with the original C code over 531,802 arguments was run once when the module was ported (see
   log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, isnan;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.j1;

/* libm's j1 and y1, under other names so that they do not clash with the module names */
extern (C) {
    pragma(mangle, "j1") double libm_j1(double x);
    pragma(mangle, "y1") double libm_y1(double x);
}

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

    /* tabulated */
    check(cephes_j1(0.0) == 0.0 && close(cephes_j1(1.0), 0.4400505857449335) && close(cephes_j1(10.0), 0.04347274616886144), "J1 at 0, 1 and 10");
    check(close(cephes_y1(1.0), -0.7812128213002887) && close(cephes_y1(10.0), 0.24901542420695388), "Y1 at 1 and 10");

    /* odd in the range where the C reflects implicitly (|x| <= 5) */
    check(cephes_j1(-3.7) == -cephes_j1(3.7), "J1 is odd for |x| <= 5");

    /* against libm for positive x, and for negative x down to -5 */
    double wj = 0, wy = 0;
    for (double x = 0.001; x < 300.0; x += 0.0731) {
        double e = fabs(cephes_j1(x) - libm_j1(x));
        if (e > wj) wj = e;
        e = fabs(cephes_y1(x) - libm_y1(x));
        if (e > wy) wy = e;
    }
    for (double x = -5.0; x < 0.0; x += 0.0731) {
        double e = fabs(cephes_j1(x) - libm_j1(x));
        if (e > wj) wj = e;
    }
    printf("  worst absolute error against libm: J1 %g, Y1 %g\n", wj, wy);
    check(wj < 1e-14 && wy < 1e-14, "J1 and Y1 agree with libm");

    /* the first nonzero zero of J1, 3.8317059702075125 */
    check(fabs(cephes_j1(3.8317059702075125)) < 1e-15, "J1 vanishes at its first zero");

    /* the upstream behaviour for x < -5: NaN, where libm gives the (odd) value (bugs.md #105) */
    check(isnan(cephes_j1(-10.0)) != 0 && !isnan(libm_j1(-10.0)), "J1(-10) is NaN in the cephes code (known, bugs.md #105)");

    /* domain error of Y1 */
    check(cephes_y1(0.0) == -MAXNUM && cephes_y1(-2.0) == -MAXNUM, "Y1 at x <= 0 is -MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
