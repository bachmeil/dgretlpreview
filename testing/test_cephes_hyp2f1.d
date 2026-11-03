/* Tests for cephes/hyp2f1.d: the Gauss hypergeometric function 2F1(a, b; c; x) against closed forms (logarithm,
   arcsine, binomial series, polynomial cases including the negative-integer cases that gretl's variant of the
   cephes code added), Gauss's theorem at x = 1, the transformations for x < -1 and x near 1, and the divergent
   case. A bit-for-bit comparison with the original C code over a grid of 42,240 parameter combinations was run
   once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : asin, fabs, log, pow, sqrt;
import core.stdc.float_ : DBL_MAX;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.hyp2f1;

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
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* trivial cases */
    check(hyp2f1(1.5, 2.5, 3.5, 0.0) == 1.0 && hyp2f1(0.0, 2.5, 3.5, 0.7) == 1.0 && hyp2f1(1.5, 0.0, 3.5, -0.7) == 1.0, "x = 0 and a or b zero give 1");

    /* 2F1(1,1;2;x) = -log(1-x)/x */
    bool ok = true;
    foreach (x; [-8.0, -2.5, -0.9, -0.3, 0.1, 0.5, 0.85, 0.95, 0.99]) {
        ok = ok && close(hyp2f1(1, 1, 2, x), -log(1 - x) / x, 1e-10);
    }
    check(ok, "2F1(1,1;2;x) = -log(1-x)/x, across the transformation ranges");

    /* 2F1(1/2,1/2;3/2;x^2) = asin(x)/x */
    ok = true;
    foreach (x; [0.1, 0.4, 0.7, 0.9]) {
        ok = ok && close(hyp2f1(0.5, 0.5, 1.5, x * x), asin(x) / x, 1e-11);
    }
    check(ok, "2F1(1/2,1/2;3/2;x^2) = asin(x)/x");

    /* 2F1(a,b;b;x) = (1-x)^(-a) */
    check(close(hyp2f1(1.7, 2.3, 2.3, 0.6), pow(0.4, -1.7)) && close(hyp2f1(2.3, 1.7, 2.3, -0.6), pow(1.6, -1.7)), "2F1(a,b;b;x) = (1-x)^-a");

    /* terminating cases: a negative integer makes the series a polynomial */
    check(close(hyp2f1(-2, 1.5, 2.5, 0.3), 1 - 2 * 1.5 * 0.3 / 2.5 + 1.5 * 2.5 * 0.09 / (2.5 * 3.5)), "negative-integer a: a polynomial");
    /* the awkward case in the file's own comment: hyp2f1(3, -2, 4, 0.99) = 1 - 1.5 x + 0.6 x^2 */
    check(close(hyp2f1(3, -2, 4, 0.99), 1 - 1.5 * 0.99 + 0.6 * 0.99 * 0.99), "hyp2f1(3, -2, 4, 0.99), the case the variant fixes");

    /* Gauss's theorem: 2F1(a,b;c;1) = Gamma(c)Gamma(c-a-b)/(Gamma(c-a)Gamma(c-b)); (1,2;5;1) = 24*1/(6*2) = 2 */
    check(close(hyp2f1(1, 2, 5, 1.0), 2.0, 1e-12), "Gauss's theorem at x = 1");

    /* divergence at x = 1 returns MAXNUM (with an error message) */
    check(hyp2f1(1, 1, 1.5, 1.0) == MAXNUM, "divergent case at x = 1 returns MAXNUM");

    /* c - a - b an integer, x near 1: the psi expansion; (1,1;2;0.95) */
    check(close(hyp2f1(1, 1, 2, 0.95), -log(0.05) / 0.95, 1e-10) && close(hyp2f1(1, 2, 4, 0.97), hyp2f1(2, 1, 4, 0.97), 1e-12), "x near 1 with integer c - a - b");

    /* Euler transformation and symmetry in a and b */
    check(close(hyp2f1(0.3, 1.2, 2.7, -0.4), hyp2f1(1.2, 0.3, 2.7, -0.4), 1e-13), "symmetric in a and b");
    check(close(hyp2f1(0.3, 1.2, 2.7, 0.4), pow(0.6, 2.7 - 0.3 - 1.2) * hyp2f1(2.7 - 0.3, 2.7 - 1.2, 2.7, 0.4), 1e-12), "Euler's transformation");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
