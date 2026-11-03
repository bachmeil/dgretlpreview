/* Tests for cephes/kv.d: the modified Bessel function of the third kind of real order, K_v(x), scaled and
   unscaled, against closed forms (half-integer orders), tabulated values, the recurrence
   K_{v+1} = K_{v-1} + (2v/x) K_v, symmetry in v, and the error cases. A comparison with the original C code
   (bit for bit over a grid of orders and arguments) was run once when the module was ported.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, isnan, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.kv;

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
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* tabulated: K_0(1) = 0.42102443824070834, K_1(1) = 0.60190723019723457, K_2(1) = 1.6248388986351774 */
    check(close(netlib_bessel_K(0.0, 1.0, 1), 0.42102443824070834) && close(netlib_bessel_K(1.0, 1.0, 1), 0.60190723019723457)
          && close(netlib_bessel_K(2.0, 1.0, 1), 1.6248388986351774), "K_0, K_1 and K_2 at 1");
    check(close(netlib_bessel_K(0.0, 0.1, 1), 2.4270690247020164) && close(netlib_bessel_K(1.0, 5.0, 1), 0.004044613445452164), "K_0(0.1) and K_1(5)");

    /* half-integer orders: K_{1/2}(x) = sqrt(pi/(2x)) e^-x and K_{3/2}(x) = K_{1/2}(x) (1 + 1/x) */
    foreach (x; [0.05, 0.7, 1.0, 3.0, 9.0, 40.0]) {
        double k12 = sqrt(PI / (2 * x)) * exp(-x);
        if (!(close(netlib_bessel_K(0.5, x, 1), k12, 1e-13) && close(netlib_bessel_K(1.5, x, 1), k12 * (1 + 1 / x), 1e-13))) {
            check(false, "half-integer orders at some x");
        }
    }
    check(true, "K_{1/2} and K_{3/2} match the closed forms on a grid");

    /* exponentially scaled values: K e^x */
    check(close(netlib_bessel_K(0.3, 2.5, 2), netlib_bessel_K(0.3, 2.5, 1) * exp(2.5), 1e-13) && close(netlib_bessel_K(1.7, 0.2, 2), netlib_bessel_K(1.7, 0.2, 1) * exp(0.2), 1e-13), "scaled values are K e^x");
    check(netlib_bessel_K(0.0, 800.0, 1) == 0.0 && netlib_bessel_K(0.0, 800.0, 2) > 0.0, "unscaled K underflows to 0 at x = 800, scaled does not");

    /* symmetry and recurrence */
    check(netlib_bessel_K(-2.3, 1.4, 1) == netlib_bessel_K(2.3, 1.4, 1), "K_{-v} = K_v");
    double v = 0.8, x = 2.2;
    check(close(netlib_bessel_K(v + 1, x, 1), netlib_bessel_K(0.2, x, 1) + netlib_bessel_K(v, x, 1) * 2 * v / x, 1e-12), "recurrence K_{v+1} = K_{v-1} + (2v/x) K_v (K_{-0.2} = K_{0.2})");

    /* errors and special values */
    check(netlib_bessel_K(1.0, 0.0, 1) == double.infinity, "K at x = 0 is infinite");
    check(isnan(netlib_bessel_K(1.0, -1.0, 1)) != 0, "negative x is a domain error (NaN)");
    check(isnan(netlib_bessel_K(double.nan, 1.0, 1)) && isnan(netlib_bessel_K(1.0, double.nan, 1)), "NaN in, NaN out");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
