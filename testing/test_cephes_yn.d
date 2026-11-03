/* Tests for cephes/yn.d: the Bessel function of the second kind of integer order, Y_n, against libm's yn,
   the parity rule for negative orders, the Wronskian J_{n+1} Y_n - J_n Y_{n+1} = 2/(pi x) (with libm's jn), the
   recurrence, and the domain error. A bit-for-bit comparison with the original C code over 369,126 (n, x) pairs
   was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM, PI;
import cephes.yn;

extern (C) {
    double jn(int n, double x);   // libm
    double yn(int n, double x);   // libm
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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* against libm */
    double worst = 0;
    foreach (n; [0, 1, 2, 3, 5, 9, 16, 30, 50]) {
        for (double x = 0.1; x < 200.0; x *= 1.17) {
            double e = fabs(cephes_bessel_Yn(n, x) - yn(n, x)) / (1e-2 + fabs(yn(n, x)));
            if (e > worst) worst = e;
        }
    }
    printf("  worst scaled error against libm yn: %g\n", worst);
    check(worst < 1e-11, "Y_n agrees with libm over n up to 50, x from 0.1 to 200");

    /* tabulated */
    check(fabs(cephes_bessel_Yn(2, 1.0) - (-1.650682606816254)) < 1e-12 && fabs(cephes_bessel_Yn(3, 4.0) - (-0.182022115953)) < 1e-11, "Y_2(1) and Y_3(4)");

    /* parity: Y_{-n} = (-1)^n Y_n */
    check(cephes_bessel_Yn(-3, 2.0) == -cephes_bessel_Yn(3, 2.0) && cephes_bessel_Yn(-4, 2.0) == cephes_bessel_Yn(4, 2.0) && cephes_bessel_Yn(-1, 2.0) == -cephes_bessel_Yn(1, 2.0), "parity in the order");

    /* Wronskian */
    bool ok = true;
    foreach (n; [0, 1, 4, 10]) {
        foreach (x; [0.7, 3.0, 12.0, 40.0]) {
            double w = jn(n + 1, x) * cephes_bessel_Yn(n, x) - jn(n, x) * cephes_bessel_Yn(n + 1, x);
            ok = ok && fabs(w * PI * x / 2.0 - 1.0) < 1e-10;
        }
    }
    check(ok, "J_{n+1} Y_n - J_n Y_{n+1} = 2/(pi x)");

    /* recurrence */
    ok = true;
    foreach (n; [2, 6, 15]) {
        foreach (x; [1.0, 5.0, 30.0]) {
            double lhs = cephes_bessel_Yn(n - 1, x) + cephes_bessel_Yn(n + 1, x);
            double rhs = 2.0 * n / x * cephes_bessel_Yn(n, x);
            ok = ok && fabs(lhs - rhs) <= 1e-11 * (1e-3 + fabs(lhs) + fabs(rhs));
        }
    }
    check(ok, "Y_{n-1} + Y_{n+1} = (2n/x) Y_n");

    /* domain: n >= 2 with x <= 0 is a singularity; n = 0, 1 are handled by y0 and y1 */
    check(cephes_bessel_Yn(3, 0.0) == -MAXNUM && cephes_bessel_Yn(2, -1.0) == -MAXNUM && cephes_bessel_Yn(0, 0.0) == -MAXNUM, "x <= 0 gives -MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
