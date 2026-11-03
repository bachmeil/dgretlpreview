/* Tests for cephes/i0.d: the modified Bessel function of the first kind of order zero, I0, and its exponentially
   scaled form i0e, against the power series (an independent computation), tabulated values, evenness, the
   scaling relation, the relation I0'(x) = I1(x) (with the D I1) by a central difference, and large arguments. A
   bit-for-bit comparison with the original C code over 359,730 arguments was run once when the module was ported
   (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.i0;
import cephes.i1 : cephes_bessel_I1;

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

/* I0(x) = sum_k (x^2/4)^k / (k!)^2 */
double i0_series(double x)
{
    double q = x * x / 4, term = 1, sum = 1;
    for (int k = 1; k < 400; k++) {
        term *= q / (cast(double) k * k);
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    return sum;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tabulated */
    check(cephes_bessel_I0(0.0) == 1.0 && close(cephes_bessel_I0(1.0), 1.2660658777520084) && close(cephes_bessel_I0(10.0), 2815.716628466254, 1e-13), "I0 at 0, 1 and 10");

    /* against the power series over both Chebyshev ranges */
    double worst = 0;
    for (double x = 0.01; x < 60.0; x *= 1.23) {
        double e = fabs(cephes_bessel_I0(x) - i0_series(x)) / i0_series(x);
        if (e > worst) worst = e;
    }
    printf("  worst relative error against the power series on 0.01..60: %g\n", worst);
    check(worst < 1e-13, "I0 agrees with the power series");

    /* even */
    check(cephes_bessel_I0(-3.2) == cephes_bessel_I0(3.2) && i0e(-12.5) == i0e(12.5), "I0 and i0e are even");

    /* scaling */
    check(close(i0e(2.5), exp(-2.5) * cephes_bessel_I0(2.5), 1e-13) && close(i0e(20.0), exp(-20.0) * cephes_bessel_I0(20.0), 1e-13), "i0e = exp(-|x|) I0");
    check(i0e(500.0) > 0.0 && i0e(500.0) < 0.05, "i0e stays finite where I0 overflows");

    /* I0' = I1, by a central difference */
    bool ok = true;
    foreach (x; [0.5, 2.0, 7.0, 12.0, 25.0]) {
        double h = 1e-5 * x;
        double d = (cephes_bessel_I0(x + h) - cephes_bessel_I0(x - h)) / (2 * h);
        ok = ok && fabs(d / cephes_bessel_I1(x) - 1.0) < 1e-7;
    }
    check(ok, "the derivative of I0 is I1");

    /* large argument: I0(x) ~ exp(x)/sqrt(2 pi x) (1 + 1/(8x) + 9/(128 x^2)) */
    enum double PI = 3.14159265358979323846;
    double x = 200.0;
    check(fabs(i0e(x) / (1.0 / sqrt(2 * PI * x) * (1 + 1.0 / (8 * x) + 9.0 / (128 * x * x))) - 1.0) < 1e-6, "asymptotic form at x = 200");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
