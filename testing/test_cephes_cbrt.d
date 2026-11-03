/* Tests for cephes/cbrt.d: the cephes cube root: perfect cubes, signs, zero, infinities and NaN, relative accuracy
   against libm's cbrtl (rounded to double) over many orders of magnitude (the cephes result differs from the
   correctly rounded one by a unit in the last place for about half of all arguments, which is why jv.d must use
   this one; the C `cbrt` symbol cannot be used as the reference here because the test is linked with libceph.a,
   whose cbrt.o would satisfy it), and cube(cbrt(x)) = x. A bit-for-bit comparison with the original C code over
   282,225 arguments was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : cbrtl, fabs, isnan;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.cbrt : cbrt;

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

    /* perfect cubes */
    bool ok = true;
    foreach (k; [1.0, 2.0, 3.0, 5.0, 10.0, 12.0, 100.0, 1000.0]) {
        ok = ok && fabs(cbrt(k * k * k) - k) <= 1e-15 * k && fabs(cbrt(-k * k * k) + k) <= 1e-15 * k;
    }
    check(ok, "cube roots of perfect cubes, both signs");
    check(cbrt(8.0) == 2.0 && cbrt(-27.0) == -3.0 && cbrt(1.0) == 1.0, "exact small cases");

    /* special values */
    check(cbrt(0.0) == 0.0 && cbrt(double.infinity) == double.infinity && cbrt(-double.infinity) == -double.infinity && isnan(cbrt(double.nan)) != 0, "zero, infinities and NaN");

    /* accuracy against libm over many orders of magnitude */
    double worst = 0;
    long ndiff = 0, n = 0;
    for (double x = 1e-300; x < 1e300; x *= 1.073) {
        double a = cbrt(x), b = cast(double) cbrtl(x);
        double e = fabs(a - b) / b;
        if (e > worst) worst = e;
        if (a != b) ndiff++;
        n++;
    }
    printf("  worst relative difference from libm: %g; %ld of %ld values differ from the correctly rounded cube root\n", worst, ndiff, n);
    check(worst < 5e-16 && ndiff > 0, "agrees with the correctly rounded cube root to within 2 units in the last place");

    /* cubing returns the argument */
    ok = true;
    foreach (x; [1e-10, 0.003, 0.5, 2.0, 7.3, 1000.0, 1e9, 3.2e20]) {
        double c = cbrt(x);
        ok = ok && fabs(c * c * c - x) <= 1e-15 * x;
    }
    check(ok, "cbrt(x)^3 = x");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
