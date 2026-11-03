/* Tests for cephes/jn.d: the Bessel function of the first kind of integer order, J_n, against libm's jn over the
   documented range x <= 30 (the code degrades silently beyond it; that is recorded), the parity rules, the
   recurrence, and the argument below machine epsilon. A bit-for-bit comparison with the original C code over
   426,442 (n, x) pairs was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.jn;

extern (C) double jn(int n, double x);   // libm

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

    /* against libm: scaled absolute error (the functions cross zero) */
    double worst = 0;
    foreach (n; [0, 1, 2, 3, 4, 7, 15, 30, 60, 120, 250]) {
        for (double x = 0.05; x < 30.0; x *= 1.19) {
            double e = fabs(cephes_bessel_Jn(n, x) - jn(n, x)) / (1e-2 + fabs(jn(n, x)));
            if (e > worst) worst = e;
        }
    }
    printf("  worst scaled error against libm jn: %g\n", worst);
    check(worst < 1e-12, "J_n agrees with libm over n up to 250, x from 0.05 to 30 (the range of the code's accuracy table)");

    /* beyond the accuracy table (x > 30) the continued fraction with its fixed 53 terms no longer converges, and the
       result silently degrades: J_3(100) is off by about 1e-5 of the scale, J_250(300) by 10% */
    check(fabs(cephes_bessel_Jn(250, 300.0) - jn(250, 300.0)) / (1e-2 + fabs(jn(250, 300.0))) > 1e-3, "known limitation: inaccurate for large x (see log.md)");

    /* tabulated */
    check(fabs(cephes_bessel_Jn(3, 2.0) - 0.12894324947440206) < 1e-14 && fabs(cephes_bessel_Jn(10, 5.0) - 0.0014678026473) < 1e-12, "J_3(2) and J_10(5)");

    /* parity: J_{-n}(x) = (-1)^n J_n(x), J_n(-x) = (-1)^n J_n(x) */
    check(cephes_bessel_Jn(-3, 2.0) == -cephes_bessel_Jn(3, 2.0) && cephes_bessel_Jn(-4, 2.0) == cephes_bessel_Jn(4, 2.0) && cephes_bessel_Jn(3, -2.0) == -cephes_bessel_Jn(3, 2.0) && cephes_bessel_Jn(4, -2.0) == cephes_bessel_Jn(4, 2.0) && cephes_bessel_Jn(-3, -2.0) == cephes_bessel_Jn(3, 2.0), "parity in the order and in the argument");

    /* recurrence J_{n-1} + J_{n+1} = (2n/x) J_n */
    bool ok = true;
    foreach (n; [3, 8, 20, 50]) {
        foreach (x; [1.5, 6.0, 14.0, 25.0]) {
            double lhs = cephes_bessel_Jn(n - 1, x) + cephes_bessel_Jn(n + 1, x);
            double rhs = 2.0 * n / x * cephes_bessel_Jn(n, x);
            ok = ok && fabs(lhs - rhs) <= 1e-11 * (1e-3 + fabs(lhs) + fabs(rhs));
        }
    }
    check(ok, "J_{n-1} + J_{n+1} = (2n/x) J_n");

    /* tiny argument with n > 2 */
    check(cephes_bessel_Jn(5, 1e-20) == 0.0 && cephes_bessel_Jn(3, 0.0) == 0.0, "x below machine epsilon gives 0 for n > 2");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
