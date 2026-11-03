/* Tests for cephes/i1.d: the modified Bessel function of the first kind of order one, I1, and its exponentially
   scaled form i1e, against the power series (an independent computation), tabulated values, oddness, the scaling
   relation, and large arguments. A bit-for-bit comparison with the original C
   code over 359,730 arguments was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.i1;

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

/* I1(x) = sum_k (x/2)^(2k+1) / (k! (k+1)!) */
double i1_series(double x)
{
    double h = x / 2, term = h, sum = h;
    for (int k = 1; k < 300; k++) {
        term *= h * h / (k * (k + 1.0));
        sum += term;
        if (fabs(term) < 1e-17 * fabs(sum)) break;
    }
    return sum;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tabulated */
    check(cephes_bessel_I1(0.0) == 0.0 && close(cephes_bessel_I1(1.0), 0.5651591039924851) && close(cephes_bessel_I1(10.0), 2670.988303701255, 1e-13), "I1 at 0, 1 and 10");

    /* against the power series over both Chebyshev ranges */
    double worst = 0;
    for (double x = 0.01; x < 60.0; x *= 1.23) {
        double e = fabs(cephes_bessel_I1(x) - i1_series(x)) / i1_series(x);
        if (e > worst) worst = e;
    }
    printf("  worst relative error against the power series on 0.01..60: %g\n", worst);
    check(worst < 1e-13, "I1 agrees with the power series");

    /* odd */
    check(cephes_bessel_I1(-3.2) == -cephes_bessel_I1(3.2) && i1e(-12.5) == -i1e(12.5), "I1 and i1e are odd");

    /* scaling */
    check(close(i1e(2.5), exp(-2.5) * cephes_bessel_I1(2.5), 1e-13) && close(i1e(20.0), exp(-20.0) * cephes_bessel_I1(20.0), 1e-13), "i1e = exp(-|x|) I1");
    check(i1e(500.0) > 0.0 && i1e(500.0) < 0.05, "i1e stays finite where I1 overflows");

    /* large argument: I1(x) ~ exp(x)/sqrt(2 pi x) (1 - 3/(8x)) */
    import core.stdc.math : sqrt;
    enum double PI = 3.14159265358979323846;
    double x = 200.0;
    check(fabs(i1e(x) / (1.0 / sqrt(2 * PI * x) * (1 - 3.0 / (8 * x) - 15.0 / (128 * x * x))) - 1.0) < 1e-6, "asymptotic form at x = 200");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
