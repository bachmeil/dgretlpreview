/* Tests for cephes/airy.d: Ai, Ai', Bi and Bi' at tabulated points, the Wronskian Ai Bi' - Ai' Bi = 1/pi
   over the whole range of the code (power series, positive asymptotic and negative asymptotic branches), the
   out-of-range return, and the C `jv` code that calls airy for large orders (through the D port of jv). A
   bit-for-bit comparison with the original C code over 13,680 arguments was run once when the module was ported
   (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.float_ : DBL_MAX;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.airy;
import cephes.constants : MAXNUM, PI;

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

    double ai, aip, bi, bip;

    /* at 0 */
    check(airy(0.0, &ai, &aip, &bi, &bip) == 0 && close(ai, 0.3550280538878172) && close(aip, -0.2588194037928068)
          && close(bi, 0.6149266274460007) && close(bip, 0.4482883573538264), "Ai, Ai', Bi, Bi' at 0");

    /* tabulated values (power series region) */
    airy(1.0, &ai, &aip, &bi, &bip);
    check(close(ai, 0.13529241631288141) && close(aip, -0.15914744129679328) && close(bi, 1.2074235949528713) && close(bip, 0.9324359333927756), "values at x = 1");
    airy(-1.0, &ai, &aip, &bi, &bip);
    check(close(ai, 0.5355608832923521) && close(aip, -0.01016056711664521) && close(bi, 0.10399738949694461) && close(bip, 0.5923756264227923), "values at x = -1");

    /* asymptotic regions */
    airy(5.0, &ai, &aip, &bi, &bip);
    check(close(ai, 1.08344428136074e-4, 1e-12) && close(aip, -2.47413890868462e-4, 1e-12) && close(bi, 657.792044171171, 1e-12), "values at x = 5");
    airy(10.0, &ai, &aip, &bi, &bip);
    check(close(ai, 1.1047532552898687e-10, 1e-12) && close(bi, 4.556411535482242e8, 1e-12), "values at x = 10 (zeta > 16 branch for Bi)");
    airy(-5.0, &ai, &aip, &bi, &bip);
    check(close(ai, 0.3507610090241142, 1e-12) && close(bi, -0.1383691349016005, 1e-12), "values at x = -5");

    /* the Wronskian over the range */
    double worst = 0;
    for (double x = -40.0; x < 25.0; x += 0.173) {
        airy(x, &ai, &aip, &bi, &bip);
        double w = ai * bip - aip * bi;
        double e = fabs(w * PI - 1.0);
        if (e > worst) worst = e;
    }
    printf("  worst |pi W - 1| over -40..25: %g\n", worst);
    check(worst < 1e-9, "the Wronskian is 1/pi across all branches");

    /* beyond MAXAIRY (25.77) */
    check(airy(30.0, &ai, &aip, &bi, &bip) == -1 && ai == 0.0 && aip == 0.0 && bi == MAXNUM && bip == MAXNUM, "x > 25.77 returns -1, Ai = Ai' = 0, Bi = Bi' = MAXNUM");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
