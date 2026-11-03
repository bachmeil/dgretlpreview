/* Tests for cephes/expx2.d: exp(x^2) and exp(-x^2) without the error amplification of squaring first, against
   the extended-precision expl(x*x), the improvement over the naive exp(x*x) for large x, symmetry in the sign of
   x, the overflow to infinity and the underflow to 0. A bit-for-bit comparison with the original C code over
   120,874 arguments was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, expl, fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.expx2;

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

    /* small and moderate arguments */
    check(expx2(0.0, 1) == 1.0 && expx2(0.0, -1) == 1.0, "x = 0");
    bool ok = true;
    foreach (x; [0.01, 0.3, 1.0, 2.5, 5.0]) {
        ok = ok && fabs(expx2(x, 1) / cast(double) expl(cast(real) x * x) - 1.0) < 3e-16;
        ok = ok && fabs(expx2(x, -1) / cast(double) expl(-cast(real) x * x) - 1.0) < 3e-16;
    }
    check(ok, "exp(x^2) and exp(-x^2) to within 1.5 units in the last place of the extended-precision value");

    /* large arguments, where the naive exp(x*x) has a relative error of x^2 * 2^-53 or more */
    double worst = 0, worstNaive = 0;
    for (double x = 15.0; x < 26.5; x += 0.37) {
        double ref_ = cast(double) expl(cast(real) x * x);
        double e = fabs(expx2(x, 1) - ref_) / ref_;
        double en = fabs(exp(x * x) - ref_) / ref_;
        if (e > worst) worst = e;
        if (en > worstNaive) worstNaive = en;
    }
    printf("  worst relative error for 15 <= x < 26.5: expx2 %g, naive exp(x*x) %g\n", worst, worstNaive);
    static if (real.sizeof > double.sizeof) {
        check(worst < 4e-16 && worstNaive > 5 * worst, "expx2 is accurate where exp(x*x) is not");
    } else {
        /* real is a double here (MSVC), so ref_ is the naive value and cannot judge anything */
        printf("  (no extended-precision type on this platform: the check against exp(x*x) is skipped)\n");
    }

    /* sign of the argument does not matter */
    check(expx2(-3.3, 1) == expx2(3.3, 1) && expx2(-3.3, -1) == expx2(3.3, -1), "even in x");

    /* limits: exp(x^2) overflows above x = sqrt(709.78) = 26.64; exp(-x^2) underflows to 0 */
    check(expx2(26.0, 1) < double.infinity && expx2(27.0, 1) == double.infinity && expx2(1e5, 1) == double.infinity, "overflow to infinity");
    check(expx2(30.0, -1) == 0.0 && expx2(1e5, -1) == 0.0, "underflow to 0 for exp(-x^2)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
