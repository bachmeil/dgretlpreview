/* Tests for cephes/incbet.d: the incomplete beta integral against closed forms (a or b equal to 1, a = b,
   (2,2) and (1/2,1/2)), the reflection I_x(a,b) = 1 - I_{1-x}(b,a), symmetry at x = 1/2 for a = b, the large
   argument (logarithm) branch, monotonicity, and the domain cases. Student's t distribution (stdtr, which calls
   incbet) is used as a consumer. A bit-for-bit comparison with the original C code
   over 10,496 (a, b, x) combinations was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : asin, fabs, pow, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.incbet;
import cephes.stdtr : stdtr;

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

    /* boundaries and domain */
    check(incbet(2, 3, 0.0) == 0.0 && incbet(2, 3, 1.0) == 1.0, "x = 0 and x = 1");
    check(incbet(-1, 3, 0.5) == 0.0 && incbet(2, 0, 0.5) == 0.0 && incbet(2, 3, -0.5) == 0.0 && incbet(2, 3, 1.5) == 0.0, "domain errors give 0");

    /* closed forms */
    bool ok = true;
    foreach (x; [0.01, 0.2, 0.5, 0.8, 0.99]) {
        ok = ok && close(incbet(1, 1, x), x) && close(incbet(3.5, 1, x), pow(x, 3.5)) && close(incbet(1, 2.5, x), 1 - pow(1 - x, 2.5));
        ok = ok && close(incbet(2, 2, x), 3 * x * x - 2 * x * x * x);
        ok = ok && close(incbet(0.5, 0.5, x), 2 / PI * asin(sqrt(x)), 1e-11);
    }
    check(ok, "I_x(1,1), I_x(a,1), I_x(1,b), I_x(2,2) and I_x(1/2,1/2) match their closed forms");

    /* reflection, over parameters that use each branch (power series, both continued fractions) */
    ok = true;
    foreach (a; [0.3, 1.5, 4.0, 12.0, 60.0]) {
        foreach (b; [0.4, 2.0, 7.0, 25.0, 90.0]) {
            foreach (x; [0.05, 0.3, 0.5, 0.7, 0.95]) {
                ok = ok && fabs(incbet(a, b, x) + incbet(b, a, 1 - x) - 1.0) < 1e-12;
            }
        }
    }
    check(ok, "I_x(a,b) + I_{1-x}(b,a) = 1 on a grid");

    /* symmetry for a = b at x = 1/2, including the logarithm branch for large arguments */
    check(close(incbet(5, 5, 0.5), 0.5) && close(incbet(50, 50, 0.5), 0.5) && fabs(incbet(1000, 1000, 0.5) - 0.5) < 1e-9 && fabs(incbet(1e5, 1e5, 0.5) - 0.5) < 1e-8, "I_{1/2}(a,a) = 1/2, also for large a");

    /* monotone in x */
    ok = true;
    double prev = -1;
    for (double x = 0.001; x < 1.0; x += 0.0137) {
        double v = incbet(7.5, 3.2, x);
        ok = ok && v >= prev;
        prev = v;
    }
    check(ok, "increasing in x");

    /* a tabulated value: I_{0.4}(30, 40) is a binomial-type tail; check against the relation with the binomial
       distribution P(X >= 30) for X ~ Bin(69, 0.4) computed by direct summation */
    double binom_tail = 0;
    {
        import core.stdc.math : lgamma, exp, log;
        for (int k = 30; k <= 69; k++) {
            double lg = lgamma(70.0) - lgamma(k + 1.0) - lgamma(70.0 - k) + k * log(0.4) + (69 - k) * log(0.6);
            binom_tail += exp(lg);
        }
    }
    check(fabs(incbet(30, 40, 0.4) - binom_tail) < 1e-10, "I_x(k, n-k+1) is the binomial upper tail");

    /* the Student t distribution calls incbet: stdtr(1, t) is the Cauchy cdf 1/2 + atan(t)/pi */
    import core.stdc.math : atan;
    check(close(stdtr(1, 1.0), 0.75, 1e-12) && close(stdtr(1, -2.0), 0.5 + atan(-2.0) / PI, 1e-12), "stdtr uses incbet (Cauchy case)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
