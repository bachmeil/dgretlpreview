/* Tests for cephes/kn.d: the modified Bessel function of the third kind of integer order, K_n(x), against
   tabulated values, the independent netlib routine (cephes.kv), the recurrence K_{n+1} = K_{n-1} + (2n/x) K_n,
   symmetry in the order, and the error cases (overflow for orders above 31, domain and singularity for x <= 0,
   underflow for large x). A bit-for-bit comparison with the original C code over 55,166 (n, x) pairs was run once
   when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.kn;
import cephes.kv : netlib_bessel_K;

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

    /* tabulated: K_0(1) = 0.42102443824070834, K_1(1) = 0.60190723019723457, K_2(1) = 1.6248388986351774 */
    check(close(cephes_bessel_Kn(0, 1.0), 0.42102443824070834) && close(cephes_bessel_Kn(1, 1.0), 0.60190723019723457) && close(cephes_bessel_Kn(2, 1.0), 1.6248388986351774), "K_0, K_1 and K_2 at 1");
    check(close(cephes_bessel_Kn(1, 5.0), 0.004044613445452164, 1e-11) && close(cephes_bessel_Kn(0, 0.1), 2.4270690247020164), "K_1(5) and K_0(0.1)");

    /* against the independent netlib routine: series range, asymptotic range, and orders up to 31 */
    double worst = 0;
    foreach (n; [0, 1, 2, 3, 5, 8, 13, 20, 31]) {
        for (double x = 0.05; x < 120.0; x *= 1.27) {
            double r = netlib_bessel_K(n, x, 1);
            if (r == 0.0 || r > 1e300) continue;
            double e = fabs(cephes_bessel_Kn(n, x) - r) / r;
            if (e > worst) worst = e;
        }
    }
    printf("  worst relative difference from netlib_bessel_K: %g\n", worst);
    /* the largest differences (4e-10) are just above the switch to the asymptotic series at x = 9.55, where the
       truncated series is least accurate (the header says it reaches 1.4e-17 only for x > 18.4) */
    check(worst < 1e-9, "agrees with the netlib routine over n <= 31 and x from 0.05 to 120");

    /* recurrence and symmetry */
    bool ok = true;
    foreach (n; [1, 2, 5, 10]) {
        foreach (x; [0.4, 2.0, 8.0, 15.0, 40.0]) {
            double lhs = cephes_bessel_Kn(n + 1, x);
            double rhs = cephes_bessel_Kn(n - 1, x) + 2.0 * n / x * cephes_bessel_Kn(n, x);
            ok = ok && close(lhs, rhs, 1e-9);
        }
    }
    check(ok, "K_{n+1} = K_{n-1} + (2n/x) K_n");
    check(cephes_bessel_Kn(-3, 2.2) == cephes_bessel_Kn(3, 2.2), "K_{-n} = K_n");

    /* errors */
    check(cephes_bessel_Kn(40, 2.0) == MAXNUM, "order above 31 overflows to MAXNUM");
    check(cephes_bessel_Kn(1, 0.0) == MAXNUM && cephes_bessel_Kn(1, -1.0) == MAXNUM, "x <= 0 gives MAXNUM");
    check(cephes_bessel_Kn(1, 800.0) == 0.0, "x above MAXLOG underflows to 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
