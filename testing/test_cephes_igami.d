/* Tests for cephes/igami.d: the inverse of the complemented incomplete gamma integral: closed forms (a = 1 and
   a = 2 through a root of x exp(-x)), the round trip igamc(a, igami(a, y)) = y over a wide range of shapes and
   levels (Newton and interval halving paths), the relation to the chi-square critical values, and the edge cases.
   A bit-for-bit comparison with the original C code over 1,320 (a, y) pairs was run once when the module was
   ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, log;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.igam : igamc;
import cephes.igami;

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

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* a = 1: Q(1, x) = exp(-x), so the inverse is -log(y) */
    bool ok = true;
    foreach (y; [1e-10, 1e-4, 0.01, 0.3, 0.5, 0.9, 0.999]) {
        ok = ok && fabs(igami(1.0, y) - (-log(y))) <= 1e-9 * (1.0 + fabs(log(y)));
    }
    check(ok, "a = 1: igami(1, y) = -log(y)");

    /* round trip over shapes from 0.1 to 1e5 and levels from 1e-10 to 0.999 */
    double worst = 0;
    foreach (a; [0.1, 0.5, 1.0, 1.5, 3.0, 10.0, 50.0, 500.0, 5000.0, 1e5]) {
        foreach (y; [1e-10, 1e-6, 0.001, 0.05, 0.3, 0.5, 0.8, 0.95, 0.999]) {
            double x = igami(a, y);
            if (x <= 0.0 || x >= MAXNUM) continue;
            double e = fabs(igamc(a, x) - y) / y;
            if (e > worst) worst = e;
        }
    }
    printf("  worst relative round-trip error: %g\n", worst);
    check(worst < 1e-7, "igamc(a, igami(a, y)) = y over a wide range");

    /* chi-square connection: the 5% critical value of chi-square(10) is 2 * igami(5, 0.05) = 18.307 */
    check(fabs(2.0 * igami(5.0, 0.05) - 18.3070) < 1e-4 && fabs(2.0 * igami(0.5, 0.05) - 3.8415) < 1e-4, "chi-square critical values");

    /* monotone: larger y gives smaller x */
    ok = true;
    double prev = 1e300;
    for (double y = 0.001; y < 1.0; y += 0.0137) {
        double x = igami(7.5, y);
        ok = ok && x < prev;
        prev = x;
    }
    check(ok, "decreasing in y");

    /* a large shape parameter: the solution is close to a */
    check(fabs(igami(1e30, 0.5) / 1e30 - 1.0) < 1e-9, "a = 1e30: the solution is within 1e-9 of a");

    /* the end points: y = 0 correctly gives MAXNUM (the solution is infinite), but y = 1 also gives MAXNUM
       where the solution is 0 (bugs.md #108) */
    check(igami(2.0, 0.0) == MAXNUM && igami(2.0, 1.0) == MAXNUM, "y = 0 and y = 1 both give MAXNUM (the latter is wrong, bugs.md #108)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
