/* Tests for cephes/igam.d: the incomplete gamma integral: left tail igam and right tail igamc, against closed
   forms (a = 1, 1/2, 2 and integer a through the Poisson sum), the complement relation across the branches (power
   series and continued fraction), large arguments, underflow, and the domain cases. The inverse igami is used to
   round-trip. A bit-for-bit comparison with the original C code over 18,180 points was run once when the
   module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : erf, exp, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.igam;
import cephes.igami : igami;

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

    /* closed forms */
    bool ok = true;
    foreach (x; [0.01, 0.4, 1.0, 2.5, 7.0, 20.0]) {
        ok = ok && close(igam(1, x), 1 - exp(-x), 1e-12) && close(igamc(1, x), exp(-x), 1e-12);
        ok = ok && close(igam(0.5, x), erf(sqrt(x)), 1e-12);
        ok = ok && close(igamc(2, x), exp(-x) * (1 + x), 1e-12);
    }
    check(ok, "a = 1, 1/2 and 2 match their closed forms");

    /* the right tail at integer a is the Poisson sum: Q(n, x) = exp(-x) sum_{k<n} x^k/k! */
    ok = true;
    foreach (n; [3, 5, 10, 25]) {
        foreach (x; [0.5, 3.0, 12.0, 30.0]) {
            double term = 1.0, sum = 1.0;
            foreach (k; 1 .. n) { term *= x / k; sum += term; }
            ok = ok && close(igamc(n, x), exp(-x) * sum, 1e-11);
        }
    }
    check(ok, "igamc(n, x) is the Poisson sum");

    /* complement across the branches */
    ok = true;
    foreach (a; [0.3, 1.0, 4.5, 20.0, 150.0]) {
        foreach (x; [0.05, 0.9, 3.0, 10.0, 60.0, 200.0]) {
            ok = ok && fabs(igam(a, x) + igamc(a, x) - 1.0) < 1e-12;
        }
    }
    check(ok, "igam + igamc = 1");

    /* large arguments: P(a, a) is a little above 1/2 */
    check(fabs(igam(1e4, 1e4) - 0.5) < 0.01 && igam(1e4, 1e4) > 0.5 && fabs(igam(1e6, 1e6) - 0.5) < 0.001, "igam(a, a) is near 1/2 for large a");

    /* domain and underflow */
    check(igam(2, 0.0) == 0.0 && igam(2, -1.0) == 0.0 && igam(0.0, 1.0) == 0.0 && igam(-1.0, 1.0) == 0.0, "igam is 0 for x <= 0 or a <= 0");
    check(igamc(2, 0.0) == 1.0 && igamc(2, -1.0) == 1.0 && igamc(0.0, 1.0) == 1.0, "igamc is 1 for x <= 0 or a <= 0");
    check(igam(1000, 1e-3) == 0.0 && igamc(1e4, 1e5) == 0.0, "underflow to 0 in the left and right tails");

    /* round trip with the inverse of igamc */
    ok = true;
    foreach (a; [0.5, 2.0, 7.5, 30.0]) {
        foreach (p; [0.01, 0.3, 0.5, 0.9]) {
            double x = igami(a, p);
            ok = ok && fabs(igamc(a, x) - p) < 1e-9;
        }
    }
    check(ok, "igamc(a, igami(a, p)) = p");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
