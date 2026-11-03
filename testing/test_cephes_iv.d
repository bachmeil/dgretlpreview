/* Tests for cephes/iv.d: the modified Bessel function of the first kind of real order, I_v(x), against the closed
   forms for orders 3/2 and 5/2, the D I0 and I1 for orders 0 and 1, the recurrence I_{v-1} - I_{v+1} = (2v/x)
   I_v, symmetry and sign rules for negative integer orders and negative arguments, and the special values. The
   half-integer order 1/2 is where hyperg's asymptotic series is chosen wrongly (bugs.md #103): that failure is
   recorded. A bit-for-bit comparison with the original C code over 88,018 (v, x) pairs was run once when the module
   was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : cosh, fabs, sinh, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : MAXNUM, PI;
import cephes.i0 : cephes_bessel_I0;
import cephes.i1 : cephes_bessel_I1;
import cephes.iv;

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

    /* orders 0 and 1 against the D I0 and I1 */
    bool ok = true;
    foreach (x; [0.3, 1.0, 2.5, 6.0, 15.0, 30.0]) {
        ok = ok && close(cephes_bessel_Iv(0.0, x), cephes_bessel_I0(x), 1e-9) && close(cephes_bessel_Iv(1.0, x), cephes_bessel_I1(x), 1e-9);
    }
    check(ok, "I_0 and I_1 agree with the Chebyshev versions");

    /* closed forms: I_{3/2}(x) = sqrt(2/(pi x)) (cosh x - sinh x / x) */
    ok = true;
    foreach (x; [0.5, 1.0, 2.0, 5.0, 12.0, 30.0]) {
        double c = sqrt(2 / (PI * x));
        ok = ok && close(cephes_bessel_Iv(1.5, x), c * (cosh(x) - sinh(x) / x), 1e-10);
        ok = ok && close(cephes_bessel_Iv(2.5, x), c * ((1 + 3 / (x * x)) * sinh(x) - 3 / x * cosh(x)), 1e-9);
    }
    check(ok, "I_{3/2} and I_{5/2} match their closed forms");

    /* order 1/2: sqrt(2/(pi x)) sinh x is only reached for larger x; at x = 2 the C is off by 1.9% (bugs.md #103) */
    double ex = sqrt(2 / (PI * 2.0)) * sinh(2.0);
    check(fabs(cephes_bessel_Iv(0.5, 2.0) / ex - 1.0) > 1e-2 && close(cephes_bessel_Iv(0.5, 12.0), sqrt(2 / (PI * 12.0)) * sinh(12.0), 1e-9), "known upstream inaccuracy of I_{1/2}(2), see bugs.md #103");

    /* recurrence at non-integer orders */
    ok = true;
    foreach (v; [0.3, 1.7, 3.4]) {
        foreach (x; [1.0, 4.0, 9.0, 20.0]) {
            double lhs = cephes_bessel_Iv(v - 1, x) - cephes_bessel_Iv(v + 1, x);
            double rhs = 2 * v / x * cephes_bessel_Iv(v, x);
            ok = ok && fabs(lhs - rhs) <= 1e-8 * (fabs(lhs) + fabs(rhs));
        }
    }
    check(ok, "I_{v-1} - I_{v+1} = (2v/x) I_v");

    /* symmetry for negative integer orders, and the sign for negative x */
    check(close(cephes_bessel_Iv(-2.0, 3.0), cephes_bessel_Iv(2.0, 3.0), 1e-12) && close(cephes_bessel_Iv(3.0, -2.0), -cephes_bessel_Iv(3.0, 2.0), 1e-12) && close(cephes_bessel_Iv(2.0, -2.0), cephes_bessel_Iv(2.0, 2.0), 1e-12), "I_{-n} = I_n and I_n(-x) = (-1)^n I_n(x)");

    /* special values */
    check(cephes_bessel_Iv(0.0, 0.0) == 1.0 && cephes_bessel_Iv(2.5, 0.0) == 0.0 && cephes_bessel_Iv(-2.5, 0.0) == MAXNUM, "x = 0");
    check(cephes_bessel_Iv(0.5, -1.0) == 0.0, "negative x with a non-integer order gives 0 (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
