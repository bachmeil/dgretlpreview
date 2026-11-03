/* Tests for cephes/polevl.d: polynomial evaluation by Horner's rule, polevl (N + 1 coefficients, highest order
   first) and p1evl (leading coefficient 1 implied), against hand-worked polynomials and direct power sums. A
   bit-for-bit comparison with the original C code over 200,000 random evaluations was run once when the module was
   ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, pow;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.polevl;

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

    /* 2x^2 + 3x + 5 */
    double[3] c = [2.0, 3.0, 5.0];
    check(polevl(0.0, c.ptr, 2) == 5.0 && polevl(1.0, c.ptr, 2) == 10.0 && polevl(-2.0, c.ptr, 2) == 7.0 && polevl(0.5, c.ptr, 2) == 7.0, "polevl on 2x^2 + 3x + 5");

    /* x^2 + 3x + 5 with the leading 1 implied */
    double[2] d = [3.0, 5.0];
    check(p1evl(0.0, d.ptr, 2) == 5.0 && p1evl(1.0, d.ptr, 2) == 9.0 && p1evl(-2.0, d.ptr, 2) == 3.0 && p1evl(10.0, d.ptr, 2) == 135.0, "p1evl on x^2 + 3x + 5");

    /* a degree-one polynomial through polevl, and a higher degree against direct power sums */
    double[2] e = [4.0, -1.0];
    check(polevl(2.5, e.ptr, 1) == 9.0, "polevl on 4x - 1");

    double[8] f = [0.5, -1.5, 2.0, 0.25, -3.0, 1.0, 7.0, -0.125];
    bool ok = true;
    foreach (x; [-2.0, -0.7, 0.0, 0.3, 1.1, 1.9]) {
        double s = 0;
        foreach (k; 0 .. 8) s += f[k] * pow(x, 7 - k);
        ok = ok && fabs(polevl(x, f.ptr, 7) - s) <= 1e-13 * (1 + fabs(s));
        double s1 = pow(x, 8);
        foreach (k; 0 .. 8) s1 += f[k] * pow(x, 7 - k);
        ok = ok && fabs(p1evl(x, f.ptr, 8) - s1) <= 1e-13 * (1 + fabs(s1));
    }
    check(ok, "degree 7 polevl and degree 8 p1evl against direct power sums");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
