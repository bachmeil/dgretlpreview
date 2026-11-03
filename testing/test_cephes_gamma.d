/* Tests for cephes/gamma.d: gamma and log-gamma at integers, half-integers, negative arguments, large and
   small arguments, the sign variable, the special values, and a comparison with libm's tgamma/lgamma over a
   grid of arguments. The incomplete gamma (igam, which calls lgam) is used as a consumer.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, isinf, isnan, lgamma, tgamma, sqrt, exp, log;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.gamma;
import cephes.igam : igam;


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
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

double relerr(double a, double b)
{
    return fabs(a - b) / fabs(b);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* integers: gamma(n) = (n - 1)! */
    check(cephes_gamma(1.0) == 1.0 && cephes_gamma(2.0) == 1.0 && close(cephes_gamma(5.0), 24.0, 1e-14) && close(cephes_gamma(10.0), 362880.0, 1e-14), "gamma at integers");
    check(close(cephes_gamma(0.5), sqrt(PI), 1e-14) && close(cephes_gamma(1.5), 0.5 * sqrt(PI), 1e-14), "gamma at half-integers");
    check(close(cephes_gamma(-0.5), -2.0 * sqrt(PI), 1e-14), "gamma(-0.5) is -2 sqrt(pi)");
    /* the sign variable is only updated in the large-negative branch (bugs.md #102): +1 for |x| <= 33 */
    check(cephes_sgngam == 1, "sgngam is +1 after gamma(-0.5), although gamma is negative there");
    cephes_gamma(-34.5);
    check(cephes_sgngam == -1 && tgamma(-34.5) < 0, "sgngam is -1 after gamma(-34.5)");
    check(close(cephes_gamma(-1.5), 4.0 / 3.0 * sqrt(PI), 1e-14) && get_cephes_sgngam() == 1, "gamma(-1.5) is positive, sgngam +1");

    /* poles and special values */
    check(isnan(cephes_gamma(0.0)) && isnan(cephes_gamma(-3.0)), "gamma at 0 and at negative integers is NaN");
    check(isnan(cephes_gamma(double.nan)) && cephes_gamma(double.infinity) == double.infinity && isnan(cephes_gamma(-double.infinity)), "gamma of NaN and of infinities");
    check(isinf(cephes_gamma(180.0)) || cephes_gamma(180.0) > 1e300, "gamma overflows past 171.6");

    /* small arguments use the 1/x expansion */
    check(close(cephes_gamma(1.0e-10), 1.0e10 - 0.5772156649015329, 1e-9), "gamma of a tiny argument");

    /* large negative arguments use the reflection formula; compare with libm */
    check(relerr(cephes_gamma(-33.5), tgamma(-33.5)) < 1e-13 && relerr(cephes_gamma(-40.25), tgamma(-40.25)) < 1e-13, "gamma at large negative arguments");

    /* grid comparison with libm */
    double worst = 0;
    for (double x = -30.9; x < 170.0; x += 0.37) {
        double r = tgamma(x);
        if (r == 0 || isnan(r) || isinf(r)) continue;
        double e = relerr(cephes_gamma(x), r);
        if (e > worst) worst = e;
    }
    printf("  worst relative error of gamma against libm: %g\n", worst);
    check(worst < 1e-12, "gamma agrees with libm over a grid");

    /* log gamma */
    check(lgam(1.0) == 0.0 && lgam(2.0) == 0.0 && close(lgam(10.0), log(362880.0), 1e-14), "lgam at 1, 2 and 10");
    check(close(lgam(0.5), 0.5 * log(PI), 1e-14) && cephes_sgngam == 1, "lgam(0.5)");
    check(close(lgam(-0.5), log(2.0 * sqrt(PI)), 1e-14) && cephes_sgngam == -1, "lgam(-0.5) and the sign of gamma(-0.5)");
    check(close(lgam(1000.0), lgamma(1000.0), 1e-13) && close(lgam(1.0e5), lgamma(1.0e5), 1e-13) && close(lgam(1.0e9), lgamma(1.0e9), 1e-13), "lgam for large arguments (Stirling branches)");
    check(close(lgam(-100.5), lgamma(-100.5), 1e-12), "lgam below -34 (reflection)");
    check(isinf(lgam(0.0)) && isinf(lgam(-5.0)) && isinf(lgam(double.infinity)) && isnan(lgam(double.nan)), "lgam at poles and special values");
    check(cephes_lgamma(7.3) == lgam(7.3), "cephes_lgamma is lgam");

    worst = 0;
    for (double x = -150.3; x < 500.0; x += 0.731) {
        double r = lgamma(x);
        if (isinf(r) || isnan(r) || fabs(r) < 1e-3) continue;
        double e = fabs(lgam(x) - r) / (1.0 + fabs(r));
        if (e > worst) worst = e;
    }
    printf("  worst error of lgam against libm: %g\n", worst);
    check(worst < 1e-12, "lgam agrees with libm over a grid");

    /* the incomplete gamma calls lgam: P(1, x) = 1 - exp(-x) */
    check(close(igam(1.0, 2.0), 1.0 - exp(-2.0), 1e-12) && close(igam(3.0, 0.5), 0.0143876739, 1e-8), "igam works through lgam");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
