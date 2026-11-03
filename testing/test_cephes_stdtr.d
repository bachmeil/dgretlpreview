/* Tests for cephes/stdtr.d: Student's t distribution function and its inverse, against closed forms (1 and 2
   degrees of freedom), tabulated critical values, symmetry, large degrees of freedom (approaching the normal),
   non-integer degrees of freedom, the inversion round trip, and the domain errors. A bit-for-bit comparison with
   the original C code over 39,150 combinations was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : atan, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM, PI;
import cephes.ndtr : ndtr;
import cephes.stdtr;

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

    /* t = 0 */
    check(stdtr(7, 0.0) == 0.5 && stdtr(2.5, 0.0) == 0.5, "median at 0");

    /* closed forms: Cauchy (k = 1) and k = 2 */
    bool ok = true;
    foreach (t; [-30.0, -3.0, -0.7, 0.4, 2.0, 11.0]) {
        ok = ok && close(stdtr(1, t), 0.5 + atan(t) / PI, 1e-12);
        ok = ok && close(stdtr(2, t), 0.5 + t / (2.0 * sqrt(2.0 + t * t)), 1e-12);
    }
    check(ok, "k = 1 (Cauchy) and k = 2 match their closed forms");

    /* tabulated two-sided 5% critical values: t_{0.975}(k) */
    check(fabs(stdtri(1, 0.975) - 12.7062) < 1e-4 && fabs(stdtri(5, 0.975) - 2.5706) < 1e-4 && fabs(stdtri(10, 0.975) - 2.2281) < 1e-4 && fabs(stdtri(30, 0.975) - 2.0423) < 1e-4, "t critical values for 1, 5, 10 and 30 degrees of freedom");
    check(fabs(stdtr(10, 2.2281) - 0.975) < 1e-5, "stdtr at a critical value");

    /* symmetry across odd, even, non-integer degrees of freedom and both code paths */
    ok = true;
    foreach (k; [1.0, 2.0, 3.0, 8.0, 15.0, 4.5, 21.3]) {
        foreach (t; [0.3, 1.5, 3.0, 8.0]) {
            ok = ok && fabs(stdtr(k, t) + stdtr(k, -t) - 1.0) < 1e-13;
        }
    }
    check(ok, "stdtr(t) + stdtr(-t) = 1");

    /* many degrees of freedom: close to the normal */
    check(fabs(stdtr(1e7, 1.96) - ndtr(1.96)) < 1e-6, "large k approaches the normal distribution");

    /* non-integer degrees of freedom through the incomplete beta: between neighbouring integer cases */
    double lo = stdtr(4, 1.2), mid = stdtr(4.5, 1.2), hi = stdtr(5, 1.2);
    check(mid > lo && mid < hi, "k = 4.5 lies between k = 4 and k = 5");

    /* inversion round trip, across the three branches of stdtri */
    ok = true;
    foreach (k; [1.0, 2.0, 5.0, 12.0, 4.5, 80.0]) {
        foreach (p; [1e-8, 0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.999, 1.0 - 1e-8]) {
            double t = stdtri(k, p);
            ok = ok && fabs(stdtr(k, t) - p) <= 1e-9 * (p < 0.5 ? p : 1.0 - p) + 1e-13;
        }
    }
    check(ok, "stdtr(k, stdtri(k, p)) = p");

    /* domain errors */
    check(stdtr(0.0, 1.0) == 0.0 && stdtr(-1.0, 1.0) == 0.0, "stdtr with k <= 0 gives 0");
    check(stdtri(5, 0.0) == 0.0 && stdtri(5, 1.0) == 0.0 && stdtri(0.0, 0.5) == 0.0, "stdtri with p outside (0, 1) or k <= 0 gives 0");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
