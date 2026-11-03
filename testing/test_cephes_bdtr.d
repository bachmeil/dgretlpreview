/* Tests for cephes/bdtr.d: the binomial distribution: bdtr (cdf), bdtrc (upper tail) and bdtri (inverse in p),
   against direct summation of the probabilities, closed forms for k = 0 and k = n, the complement relation, the
   round trip bdtr(k, n, bdtri(k, n, y)) = y, and the domain errors. A bit-for-bit comparison with the original C
   code over 3,366 combinations was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, lgamma, log, pow;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.bdtr;

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

/* P(X = j) for X ~ Bin(n, p) */
double pmf(int j, int n, double p)
{
    return exp(lgamma(n + 1.0) - lgamma(j + 1.0) - lgamma(n - j + 1.0) + j * log(p) + (n - j) * log(1 - p));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* direct summation */
    bool ok = true;
    foreach (n; [1, 5, 20, 60, 300]) {
        foreach (p; [0.05, 0.3, 0.5, 0.9]) {
            for (int k = 0; k < n; k += (n > 20 ? n / 6 : 1)) {
                double s = 0;
                foreach (j; 0 .. k + 1) s += pmf(j, n, p);
                ok = ok && fabs(bdtr(k, n, p) - s) <= 1e-10 * s + 1e-14;
                ok = ok && fabs(bdtrc(k, n, p) - (1 - s)) <= 1e-10 + 1e-10 * (1 - s);
            }
        }
    }
    check(ok, "bdtr and bdtrc agree with direct summation");

    /* closed forms and special cases */
    check(close(bdtr(0, 10, 0.3), pow(0.7, 10)) && bdtr(10, 10, 0.3) == 1.0 && bdtrc(10, 10, 0.3) == 0.0 && bdtrc(-1, 10, 0.3) == 1.0, "k = 0, k = n and k < 0");
    check(close(bdtrc(0, 10, 0.3), 1 - pow(0.7, 10)) && close(bdtrc(0, 10, 0.001), 1 - pow(0.999, 10), 1e-9), "upper tail at k = 0, both branches");
    check(close(bdtr(3, 8, 0.5) + bdtrc(3, 8, 0.5), 1.0, 1e-14), "bdtr + bdtrc = 1");

    /* inverse */
    ok = true;
    foreach (n; [1, 4, 25, 200]) {
        for (int k = 0; k < n; k += (n > 25 ? n / 5 : 1)) {
            foreach (y; [0.01, 0.2, 0.5, 0.8, 0.99]) {
                double p = bdtri(k, n, y);
                ok = ok && fabs(bdtr(k, n, p) - y) <= 1e-8;
            }
        }
    }
    check(ok, "bdtr(k, n, bdtri(k, n, y)) = y");
    check(close(bdtri(0, 1, 0.25), 0.75, 1e-10), "bdtri(0, 1, y) = 1 - y");

    /* domain errors */
    check(bdtr(3, 10, -0.1) == 0.0 && bdtr(3, 10, 1.2) == 0.0 && bdtr(-1, 10, 0.5) == 0.0 && bdtr(11, 10, 0.5) == 0.0, "bdtr domain errors give 0");
    check(bdtrc(11, 10, 0.5) == 0.0 && bdtrc(3, 10, 2.0) == 0.0, "bdtrc domain errors give 0");
    check(bdtri(3, 10, -0.1) == 0.0 && bdtri(3, 10, 1.1) == 0.0 && bdtri(10, 10, 0.5) == 0.0, "bdtri domain errors give 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
