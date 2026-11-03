/* Tests for cephes/k1.d: the modified Bessel function of the third kind of order one, K1, and its exponentially
   scaled form k1e, against tabulated values, the independent netlib routine (cephes.kv), the scaling relation,
   the Wronskian I0 K1 + I1 K0 = 1/x (with the D I0, I1 and K0), and the domain errors. A bit-for-bit comparison
   with the original C code over 484,604 arguments was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.i0 : cephes_bessel_I0;
import cephes.i1 : cephes_bessel_I1;
import cephes.k0 : cephes_bessel_K0;
import cephes.k1;
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

bool close(double a, double b, double tol = 1e-13)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tabulated */
    check(close(cephes_bessel_K1(1.0), 0.6019072301972346) && close(cephes_bessel_K1(0.1), 9.853844780870606) && close(cephes_bessel_K1(5.0), 0.004044613445452164, 1e-12), "K1 at 0.1, 1 and 5");

    /* against the netlib routine, both Chebyshev ranges */
    double worst = 0;
    for (double x = 0.01; x < 600.0; x *= 1.19) {
        double r = netlib_bessel_K(1.0, x, 1);
        if (r == 0.0) continue;
        double e = fabs(cephes_bessel_K1(x) - r) / r;
        if (e > worst) worst = e;
    }
    printf("  worst relative difference from netlib_bessel_K(1, x): %g\n", worst);
    check(worst < 1e-9, "K1 agrees with the netlib routine over 0.01..600");

    /* scaling */
    check(close(k1e(0.7), exp(0.7) * cephes_bessel_K1(0.7), 1e-12) && close(k1e(6.0), exp(6.0) * cephes_bessel_K1(6.0), 1e-12), "k1e = exp(x) K1");
    check(k1e(700.0) > 0.0 && k1e(700.0) < 0.1, "k1e stays finite where K1 underflows");

    /* Wronskian: I0 K1 + I1 K0 = 1/x */
    bool ok = true;
    foreach (x; [0.3, 1.0, 1.9, 2.1, 6.0, 15.0]) {
        double w = cephes_bessel_I0(x) * cephes_bessel_K1(x) + cephes_bessel_I1(x) * cephes_bessel_K0(x);
        ok = ok && fabs(w * x - 1.0) < 1e-12;
    }
    check(ok, "I0 K1 + I1 K0 = 1/x with all four functions from the D port");

    /* domain */
    check(cephes_bessel_K1(0.0) == MAXNUM && cephes_bessel_K1(-1.0) == MAXNUM && k1e(0.0) == MAXNUM && k1e(-2.0) == MAXNUM, "x <= 0 gives MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
