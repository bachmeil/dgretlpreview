/* Tests for cephes/nbdtr.d: the negative binomial distribution: nbdtr (cdf), nbdtrc (upper tail) and nbdtri
   (inverse in p), against direct summation of the probabilities, the complement relation, the inversion round
   trip, the geometric case n = 1, and the domain errors. nbdtri is not called with n <= 0, which makes the C
   incbi loop forever (bugs.md #107). A bit-for-bit comparison with the original C code over 1,680 combinations
   was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, lgamma, log, pow;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.nbdtr;

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

/* P(X = j) for X the number of failures before the n-th success */
double pmf(int j, int n, double p)
{
    return exp(lgamma(n + j + 0.0) - lgamma(j + 1.0) - lgamma(n + 0.0) + n * log(p) + j * log(1 - p));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* direct summation */
    bool ok = true;
    foreach (n; [1, 2, 5, 12, 40]) {
        foreach (p; [0.1, 0.35, 0.6, 0.9]) {
            foreach (k; [0, 1, 4, 10, 30]) {
                double s = 0;
                foreach (j; 0 .. k + 1) s += pmf(j, n, p);
                ok = ok && fabs(nbdtr(k, n, p) - s) <= 1e-10 * s + 1e-14;
                ok = ok && fabs(nbdtrc(k, n, p) - (1 - s)) <= 1e-9;
            }
        }
    }
    check(ok, "nbdtr and nbdtrc agree with direct summation");

    /* geometric distribution, n = 1: P(X <= k) = 1 - (1-p)^(k+1) */
    check(close(nbdtr(3, 1, 0.3), 1 - pow(0.7, 4)) && close(nbdtrc(3, 1, 0.3), pow(0.7, 4)), "n = 1 is the geometric distribution");
    check(close(nbdtr(5, 4, 0.5) + nbdtrc(5, 4, 0.5), 1.0, 1e-13), "nbdtr + nbdtrc = 1");

    /* inverse */
    ok = true;
    foreach (n; [1, 3, 10, 50]) {
        foreach (k; [0, 2, 8, 40]) {
            foreach (y; [0.01, 0.2, 0.5, 0.8, 0.99]) {
                double p = nbdtri(k, n, y);
                ok = ok && fabs(nbdtr(k, n, p) - y) <= 1e-8;
            }
        }
    }
    check(ok, "nbdtr(k, n, nbdtri(k, n, y)) = y");

    /* domain errors */
    check(nbdtr(-1, 3, 0.5) == 0.0 && nbdtr(2, 3, -0.1) == 0.0 && nbdtr(2, 3, 1.1) == 0.0, "nbdtr domain errors give 0");
    check(nbdtrc(-1, 3, 0.5) == 0.0 && nbdtrc(2, 3, 1.5) == 0.0, "nbdtrc domain errors give 0");
    check(nbdtri(-1, 3, 0.5) == 0.0 && nbdtri(2, 3, -0.2) == 0.0 && nbdtri(2, 3, 1.2) == 0.0, "nbdtri domain errors give 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
