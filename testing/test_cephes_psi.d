/* Tests for cephes/psi.d: the digamma function at tabulated values (integers, halves, negative arguments), the
   recurrence psi(x + 1) = psi(x) + 1/x, the reflection formula, large arguments, and the poles. A bit-for-bit
   comparison with the original C code over 41,252 arguments was run once when the module was ported (see
   log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, log, tan;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM, PI;
import cephes.psi;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    enum double EULER = 0.5772156649015329;

    /* integers: psi(n) = -gamma + sum_{k<n} 1/k */
    check(close(psi(1.0), -EULER) && close(psi(2.0), 1 - EULER) && close(psi(5.0), 1 + 0.5 + 1.0 / 3 + 0.25 - EULER) && close(psi(10.0), 2.2517525890667211, 1e-13), "psi at positive integers up to 10");

    /* halves: psi(1/2) = -gamma - 2 ln 2 */
    check(close(psi(0.5), -EULER - 2 * log(2.0)) && close(psi(1.5), 2 - EULER - 2 * log(2.0)), "psi(1/2) and psi(3/2)");

    /* a larger argument */
    check(close(psi(100.0), 4.600161852738087, 1e-13), "psi(100)");

    /* recurrence */
    bool ok = true;
    foreach (x; [0.1, 0.9, 2.3, 7.7, 11.5, 40.2, 1000.5]) {
        ok = ok && fabs(psi(x + 1) - psi(x) - 1 / x) < 1e-12 * (1 + 1 / x);
    }
    check(ok, "psi(x + 1) = psi(x) + 1/x");

    /* negative arguments: psi(-1/2) = psi(1/2) + 2, and the reflection psi(1-x) - psi(x) = pi cot(pi x) */
    check(close(psi(-0.5), psi(0.5) + 2.0, 1e-13), "psi(-1/2)");
    ok = true;
    foreach (x; [0.2, 0.37, 0.8, 2.6, 5.9]) {
        ok = ok && fabs((psi(1 - x) - psi(x)) - PI / tan(PI * x)) < 1e-12 * (1 + fabs(PI / tan(PI * x)));
    }
    check(ok, "psi(1 - x) - psi(x) = pi cot(pi x)");
    ok = true;
    foreach (x; [-0.3, -1.7, -3.2, -8.9, -15.5]) {
        ok = ok && fabs(psi(x + 1) - psi(x) - 1 / x) < 1e-11 * (1 + 1 / fabs(x));
    }
    check(ok, "the recurrence also holds for negative arguments");

    /* large x: psi(x) ~ log x - 1/(2x) - 1/(12 x^2) */
    check(close(psi(1.0e6), log(1.0e6) - 0.5e-6 - 1.0 / (12 * 1.0e12), 1e-14), "large x");

    /* poles */
    check(psi(0.0) == MAXNUM && psi(-1.0) == MAXNUM && psi(-7.0) == MAXNUM, "poles at 0 and the negative integers give MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
