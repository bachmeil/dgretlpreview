/* Tests for cephes/j0.d: the Bessel functions of order zero, J0 and Y0, against libm's j0/y0 over the three
   regimes (small argument, x <= 5, asymptotic), tabulated values, the zeros, symmetry of J0, and the domain
   error of Y0. A bit-for-bit comparison with the original C code over 265,903 arguments was run once when the
   module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.j0;

/* libm's j0 and y0, under other names so that they do not clash with the module names */
extern (C) {
    pragma(mangle, "j0") double libm_j0(double x);
    pragma(mangle, "y0") double libm_y0(double x);
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
    check(cephes_j0(0.0) == 1.0 && close(cephes_j0(1.0), 0.7651976865579666) && close(cephes_j0(10.0), -0.2459357644513483), "J0 at 0, 1 and 10");
    check(close(cephes_y0(1.0), 0.08825696421567697) && close(cephes_y0(10.0), 0.05567116728359939), "Y0 at 1 and 10");

    /* symmetry */
    check(cephes_j0(-3.7) == cephes_j0(3.7), "J0 is even");

    /* against libm over the regimes; absolute error (the functions cross zero) */
    double wj = 0, wy = 0;
    for (double x = 0.001; x < 300.0; x += 0.0731) {
        double e = fabs(cephes_j0(x) - libm_j0(x));
        if (e > wj) wj = e;
        e = fabs(cephes_y0(x) - libm_y0(x));
        if (e > wy) wy = e;
    }
    printf("  worst absolute error against libm over 0.001..300: J0 %g, Y0 %g\n", wj, wy);
    check(wj < 1e-14 && wy < 1e-14, "J0 and Y0 agree with libm");

    /* the first zero of J0, 2.404825557695773 */
    check(fabs(cephes_j0(2.404825557695773)) < 1e-15, "J0 vanishes at its first zero");

    /* tiny argument branch: 1 - x^2/4 */
    check(cephes_j0(1.0e-6) == 1.0 - 1.0e-12 / 4.0, "J0 for x < 1e-5 is 1 - x^2/4");

    /* domain error of Y0 */
    check(cephes_y0(0.0) == -MAXNUM && cephes_y0(-1.0) == -MAXNUM, "Y0 at x <= 0 is -MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
