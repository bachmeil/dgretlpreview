/* Tests for cephes/pdtr.d: the Poisson distribution: pdtr (cdf), pdtrc (upper tail) and pdtri (inverse in the
   mean), against direct summation of the probabilities, the complement relation, the k = 0 closed form, the
   inversion round trip (pdtri uses the still-C igami), and the domain errors. A bit-for-bit comparison with the
   original C code over 306 combinations was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.pdtr;

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

    /* direct summation: P(X <= k) = exp(-m) sum_{j<=k} m^j / j! */
    bool ok = true;
    foreach (m; [0.1, 0.8, 3.0, 12.0, 40.0]) {
        foreach (k; [0, 1, 2, 5, 10, 30, 60]) {
            double term = exp(-m), s = term;
            foreach (j; 1 .. k + 1) { term *= m / j; s += term; }
            ok = ok && fabs(pdtr(k, m) - s) <= 1e-11 * s + 1e-15;
            ok = ok && fabs(pdtrc(k, m) - (1 - s)) <= 1e-10;
        }
    }
    check(ok, "pdtr and pdtrc agree with direct summation");

    /* k = 0 */
    check(close(pdtr(0, 2.0), exp(-2.0)) && close(pdtrc(0, 2.0), 1 - exp(-2.0)), "k = 0 closed forms");
    check(close(pdtr(4, 6.5) + pdtrc(4, 6.5), 1.0, 1e-14), "pdtr + pdtrc = 1");

    /* inverse: pdtr(k, pdtri(k, y)) = y */
    ok = true;
    foreach (k; [0, 1, 3, 9, 40]) {
        foreach (y; [0.01, 0.2, 0.5, 0.8, 0.99]) {
            double m = pdtri(k, y);
            ok = ok && fabs(pdtr(k, m) - y) <= 1e-9;
        }
    }
    check(ok, "pdtr(k, pdtri(k, y)) = y");

    /* domain errors */
    check(pdtr(-1, 2.0) == 0.0 && pdtr(3, 0.0) == 0.0 && pdtr(3, -1.0) == 0.0, "pdtr domain errors give 0");
    check(pdtrc(-1, 2.0) == 0.0 && pdtrc(3, 0.0) == 0.0, "pdtrc domain errors give 0");
    check(pdtri(-1, 0.5) == 0.0 && pdtri(2, -0.1) == 0.0 && pdtri(2, 1.0) == 0.0, "pdtri domain errors give 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
