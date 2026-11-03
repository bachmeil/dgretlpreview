/* Tests for cephes/chbevl.d: the Chebyshev series evaluation against a direct summation with the Chebyshev
   polynomials T_k(y) = cos(k acos y) of the half argument y = x/2 (the convention of cephes: the series is
   evaluated at x/2, for x in [-2, 2]), with the coefficients stored highest order first and the order-zero term
   halved: result = sum_k c_k T_k(x/2) - c_0 / 2, c_0 being the last array element, small cases worked by hand, and a
   constant series. The users (i0, i1, k0, k1) are tested through the Bessel functions. A bit-for-bit comparison
   with the original C code over 60,000 random series was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : acos, cos, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.chbevl;

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

/* sum_{k=0}^{n-1} coef[n-1-k] T_k(x/2) - coef[n-1]/2 */
double direct(double x, const(double)[] coef)
{
    int n = cast(int) coef.length;
    double s = 0;
    foreach (k; 0 .. n) s += coef[n - 1 - k] * cos(k * acos(x / 2));
    return s - 0.5 * coef[n - 1];
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* hand cases: [0, 1] is the constant series c_0 = 1, which evaluates to 1/2 */
    double[2] c01 = [0.0, 1.0];
    check(chbevl(0.3, c01.ptr, 2) == 0.5 && chbevl(-1.6, c01.ptr, 2) == 0.5, "[0, 1] is the constant 1/2");

    /* [1, 0, 0] is T_2(x/2) = x^2/2 - 1 */
    double[3] c2 = [1.0, 0.0, 0.0];
    bool ok = true;
    foreach (x; [-2.0, -0.8, 0.0, 1.2, 2.0]) ok = ok && fabs(chbevl(x, c2.ptr, 3) - (x * x / 2 - 1)) < 1e-15;
    check(ok, "[1, 0, 0] is T_2(x/2)");

    /* [1, 0] is T_1(x/2) = x/2 (the constant term c_0 = 0) */
    double[2] c1 = [1.0, 0.0];
    check(chbevl(0.74, c1.ptr, 2) == 0.37, "[1, 0] is T_1(x/2)");

    /* random series against the direct sum */
    double[12] c = [0.3, -1.2, 0.7, 2.5, -0.1, 0.9, -3.3, 1.1, 0.4, -0.6, 1.7, 0.25];
    double worst = 0;
    foreach (n; [2, 3, 5, 8, 12]) {
        for (double x = -1.98; x < 2.0; x += 0.1463) {
            double e = fabs(chbevl(x, c.ptr, n) - direct(x, c[0 .. n]));
            if (e > worst) worst = e;
        }
    }
    printf("  worst absolute difference from the direct sum: %g\n", worst);
    check(worst < 1e-13, "agrees with direct summation of Chebyshev polynomials");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
