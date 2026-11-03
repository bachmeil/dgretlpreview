/* Tests for cephes/ndtri.d: the inverse of the normal distribution function at tabulated quantiles, symmetry,
   the three approximation ranges (central, and the two tail ranges split at z = 8), round trips with ndtr over
   many orders of magnitude, and the domain errors. A bit-for-bit comparison with the original C code over 16,970
   arguments was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, isnan;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM;
import cephes.ndtr : ndtr;
import cephes.ndtri;

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

    /* tabulated quantiles */
    check(ndtri(0.5) == 0.0, "ndtri(1/2) = 0");
    check(close(ndtri(0.975), 1.959963984540054) && close(ndtri(0.95), 1.6448536269514722) && close(ndtri(0.995), 2.5758293035489004), "quantiles 0.975, 0.95 and 0.995");
    check(close(ndtri(0.001), -3.090232306167813) && close(ndtri(1e-10), -6.361340902404056, 1e-11) && close(ndtri(1e-20), -9.262340089798408, 1e-11), "lower-tail quantiles (z < 8 and z >= 8 branches)");

    /* symmetry (p not so small that 1 - p loses it) */
    bool ok = true;
    foreach (p; [1e-4, 0.01, 0.05, 0.2, 0.4]) {
        ok = ok && fabs(ndtri(p) + ndtri(1.0 - p)) <= 1e-11 * (1.0 + fabs(ndtri(p)));
    }
    check(ok, "ndtri(p) = -ndtri(1-p)");

    /* round trip with ndtr over a wide range, in relative error of p */
    double worst = 0;
    for (double p = 1e-300; p < 0.5; p *= 3.1) {
        double e = fabs(ndtr(ndtri(p)) - p) / p;
        if (e > worst) worst = e;
    }
    printf("  worst relative round-trip error over 1e-300..0.5: %g\n", worst);
    check(worst < 1e-9, "ndtr(ndtri(p)) = p from 1e-300 to 1/2");

    /* domain errors */
    check(ndtri(0.0) == -MAXNUM && ndtri(-0.5) == -MAXNUM && ndtri(1.0) == MAXNUM && ndtri(2.0) == MAXNUM, "p <= 0 gives -MAXNUM and p >= 1 gives MAXNUM (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
