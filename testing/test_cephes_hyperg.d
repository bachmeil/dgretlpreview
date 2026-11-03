/* Tests for cephes/hyperg.d: the confluent hypergeometric function 1F1(a; b; x) against closed forms
   (exponential, error function, terminating polynomials), Kummer's transformation, the power-series and
   asymptotic-series regimes (large |x|), the singular case, and the asymptotic 2F0 helper. A bit-for-bit
   comparison with the original C code over 8,008 parameter combinations was run once when the module was ported
   (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : erf, exp, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM, PI;
import cephes.hyperg;

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

bool close(double a, double b, double tol = 1e-11)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* 1F1(a; a; x) = exp(x) (the Kummer branch takes this one) */
    check(close(hyperg(1.7, 1.7, 2.0), exp(2.0)) && close(hyperg(2.5, 2.5, -3.0), exp(-3.0)), "1F1(a;a;x) = e^x");

    /* 1F1(1; 2; x) = (e^x - 1)/x, small and large |x| (power series and asymptotic series) */
    bool ok = true;
    foreach (x; [-60.0, -20.0, -0.5, 0.01, 0.5, 3.0, 40.0]) {
        ok = ok && close(hyperg(1, 2, x), (exp(x) - 1) / x, 1e-8);
    }
    check(ok, "1F1(1;2;x) = (e^x - 1)/x at small and at large |x|");
    /* the mid-range failures of the C code (bugs.md #103): 1F1(1;2;-3) comes out as 1/3 instead of 0.3167 and
       1F1(1;2;12) is off by 6e-6, because the asymptotic series, which terminates when 1+a-b = 0, is chosen on a
       misleading error estimate */
    check(fabs(hyperg(1, 2, -3.0) - (exp(-3.0) - 1) / -3.0) > 1e-2 && fabs(hyperg(1, 2, -3.0) - 1.0 / 3.0) < 1e-12, "known upstream inaccuracy at 1F1(1;2;-3), see bugs.md #103");

    /* 1F1(1; 3; x) = 2 (e^x - 1 - x)/x^2 */
    check(close(hyperg(1, 3, -20.0), 2 * (exp(-20.0) - 1 + 20.0) / 400.0, 1e-8) && close(hyperg(1, 3, 40.0), 2 * (exp(40.0) - 1 - 40.0) / 1600.0, 1e-11), "1F1(1;3;x) = 2(e^x - 1 - x)/x^2");

    /* 1F1(1/2; 3/2; -x^2) = sqrt(pi) erf(x) / (2x) */
    ok = true;
    foreach (x; [0.3, 1.0, 2.0]) {
        ok = ok && close(hyperg(0.5, 1.5, -x * x), sqrt(PI) * erf(x) / (2 * x), 1e-11);
    }
    check(ok && close(hyperg(0.5, 1.5, -16.0), sqrt(PI) * erf(4.0) / 8.0, 1e-7), "1F1(1/2;3/2;-x^2) = sqrt(pi) erf(x)/(2x)");

    /* terminating series: 1F1(-2; 3; x) = 1 - 2x/3 + x^2/12 */
    check(close(hyperg(-2, 3, 1.5), 1 - 2 * 1.5 / 3 + 1.5 * 1.5 / 12) && close(hyperg(-2, 3, -4.0), 1 + 8.0 / 3 + 16.0 / 12), "terminating polynomial");

    /* Kummer's transformation: 1F1(a;b;x) = e^x 1F1(b-a;b;-x) */
    check(close(hyperg(0.7, 2.9, 1.8), exp(1.8) * hyperg(2.2, 2.9, -1.8), 1e-12) && close(hyperg(1.3, 4.1, -6.0), exp(-6.0) * hyperg(2.8, 4.1, 6.0), 1e-11), "Kummer's transformation");

    /* singularity: b = 0 with a != 0 */
    check(hyperg(1.0, 0.0, 1.0) == MAXNUM, "b = 0 is a singularity: MAXNUM");

    /* the asymptotic 2F0 helper: 2F0(1, 1; -1/x) ~ x e^x E1(x), compare with the exact value 0.5963 at x = 10 */
    double err = 0;
    double v = hyp2f0(1.0, 1.0, -0.1, 1, &err);
    check(fabs(v - 0.91563333) < 1e-3 && err < 0.1, "hyp2f0(1,1;-1/10) is about x e^x E1(x) = 0.9156");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
