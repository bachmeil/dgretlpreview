/* Tests for cephes/incbi.d: the inverse of the incomplete beta integral against closed forms, round trips
   incbet(a, b, incbi(a, b, y)) = y across the code's paths (approximation plus Newton, interval halving when a or b
   is at most 1, reversal for large y), symmetry, and the boundary values. The D bdtri (which calls incbi for
   the binomial) is used as an independent consumer. A bit-for-bit comparison with the original C code over
   3,136 combinations was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.bdtr : bdtri;
import cephes.incbet : incbet;
import cephes.incbi;


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

    /* boundaries */
    check(incbi(2, 3, 0.0) == 0.0 && incbi(2, 3, -0.1) == 0.0 && incbi(2, 3, 1.0) == 1.0 && incbi(2, 3, 1.5) == 1.0, "y <= 0 gives 0 and y >= 1 gives 1");

    /* closed forms: I_x(1,1) = x, I_x(2,1) = x^2, I_x(1,2) = 1 - (1-x)^2 */
    bool ok = true;
    foreach (y; [0.001, 0.1, 0.5, 0.9, 0.999]) {
        ok = ok && close(incbi(1, 1, y), y) && close(incbi(2, 1, y), sqrt(y)) && close(incbi(1, 2, y), 1 - sqrt(1 - y));
    }
    check(ok, "inverse closed forms for (1,1), (2,1) and (1,2)");

    /* the median of a symmetric beta is 1/2 */
    check(close(incbi(3, 3, 0.5), 0.5) && close(incbi(40, 40, 0.5), 0.5) && close(incbi(0.7, 0.7, 0.5), 0.5), "incbi(a, a, 1/2) = 1/2");

    /* round trips, for y from 1e-4 to 0.99. At the extremes (y = 1e-8 with a >= 150 and b <= 2.2, or y = 0.9999 with
       b = 0.3, where x is within a few ulp of 1) the original C code, too, only reaches 1e-8 to 1e-4 relative accuracy in
       y; those cases are left out. */
    double worst = 0;
    foreach (a; [0.2, 0.9, 1.0, 1.7, 5.0, 20.0, 150.0, 2000.0]) {
        foreach (b; [0.3, 1.0, 2.2, 9.0, 60.0, 700.0]) {
            foreach (y; [1e-4, 0.01, 0.2, 0.5, 0.8, 0.99]) {
                double x = incbi(a, b, y);
                double e = fabs(incbet(a, b, x) - y) / y;
                if (x > 0 && x < 1 && e > worst) worst = e;
            }
        }
    }
    printf("  worst relative round-trip error: %g\n", worst);
    check(worst < 1e-7, "incbet(a, b, incbi(a, b, y)) = y over a grid");

    /* reflection: incbi(a, b, y) = 1 - incbi(b, a, 1 - y) */
    ok = true;
    foreach (y; [0.05, 0.3, 0.6, 0.95]) {
        ok = ok && fabs(incbi(4.5, 11.0, y) - (1.0 - incbi(11.0, 4.5, 1.0 - y))) < 1e-9;
    }
    check(ok, "reflection incbi(a,b,y) = 1 - incbi(b,a,1-y)");

    /* consumer: the binomial success probability with P(X <= k) = y; for n = 1, k = 0, y = 0.25 this is 0.75 */
    check(close(bdtri(0, 1, 0.25), 0.75, 1e-10), "bdtri uses incbi (n = 1)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
