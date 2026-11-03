/* Tests for cephes/jv.d: the Bessel function of the first kind of real order J_v(x), and Hankel's expansion.
   Integer orders are compared with libm's jn (including orders above 500, which use the large-order
   expansions); half-integer orders with their closed forms; other orders with the recurrence
   J_{v-1} + J_{v+1} = (2v/x) J_v; plus negative orders and arguments and the error case. A bit-for-bit
   comparison with the original C code over a large grid was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : cos, fabs, isnan, sin, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.jv;

extern (C) double jn(int n, double x);   // libm

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
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* integer orders against libm */
    double worst = 0;
    foreach (n; [0, 1, 2, 3, 7, 15, 30, 80, 200, 450]) {
        for (double x = 0.3; x < 400.0; x *= 1.31) {
            double r = jn(n, x);
            double e = fabs(cephes_bessel_Jv(n, x) - r) / (1e-3 + fabs(r));
            if (e > worst) worst = e;
        }
    }
    printf("  worst scaled error against libm jn (orders up to 450): %g\n", worst);
    version (Windows) {
        /* the jn of the Windows C library is not accurate for large orders (the reference is the less exact side) */
        check(worst < 1e-5, "integer orders agree with libm");
    } else {
        check(worst < 1e-9, "integer orders agree with libm");
    }

    /* orders above 500 use the uniform and transitional expansions */
    worst = 0;
    foreach (n; [510, 600, 800, 1500]) {
        foreach (f; [0.6, 0.9, 0.99, 1.0, 1.01, 1.1, 1.4]) {
            double x = n * f;
            double r = jn(n, x);
            double e = fabs(cephes_bessel_Jv(n, x) - r) / (1e-3 + fabs(r));
            if (e > worst) worst = e;
        }
    }
    printf("  worst scaled error for orders above 500: %g\n", worst);
    version (Windows) {
        printf("  (the jn of the Windows C library is too inaccurate for orders above 500: no check)\n");
    } else {
        check(worst < 1e-6, "orders above 500 agree with libm");
    }

    /* half-integer orders */
    bool ok = true;
    foreach (x; [0.4, 1.0, 3.0, 8.0, 25.0, 60.0]) {
        double c = sqrt(2 / (PI * x));
        ok = ok && close(cephes_bessel_Jv(0.5, x), c * sin(x), 1e-12);
        ok = ok && close(cephes_bessel_Jv(1.5, x), c * (sin(x) / x - cos(x)), 1e-11);
        ok = ok && close(cephes_bessel_Jv(-0.5, x), c * cos(x), 1e-12);
    }
    check(ok, "J_{1/2}, J_{3/2} and J_{-1/2} match their closed forms");

    /* recurrence at non-integer orders */
    ok = true;
    foreach (v; [0.3, 1.7, 5.2, 12.9, 33.4]) {
        foreach (x; [0.8, 4.0, 9.5, 31.0]) {
            double lhs = cephes_bessel_Jv(v - 1, x) + cephes_bessel_Jv(v + 1, x);
            double rhs = 2 * v / x * cephes_bessel_Jv(v, x);
            ok = ok && fabs(lhs - rhs) <= 1e-9 * (1e-3 + fabs(lhs) + fabs(rhs));
        }
    }
    check(ok, "the three-term recurrence holds at non-integer orders");

    /* integer order, negative order and argument: J_{-n} = (-1)^n J_n, J_n(-x) = (-1)^n J_n(x) */
    check(close(cephes_bessel_Jv(-3, 5.0), -cephes_bessel_Jv(3, 5.0)) && close(cephes_bessel_Jv(3, -5.0), -cephes_bessel_Jv(3, 5.0)) && close(cephes_bessel_Jv(-4, -5.0), cephes_bessel_Jv(4, 5.0)), "reflection for integer orders");

    /* small and zero arguments */
    check(cephes_bessel_Jv(2.5, 0.0) == 0.0 && cephes_bessel_Jv(0.0, 0.0) == 1.0, "J_v(0) is 0 for v > 0 and 1 for v = 0");
    check(close(cephes_bessel_Jv(7.2, 0.01), 0.0, 1e-9) && cephes_bessel_Jv(7.2, 0.01) > 0, "tiny positive value for small x and large order");

    /* Hankel's expansion at large x */
    check(close(cephes_hankel(0.0, 100.0), jn(0, 100.0), 1e-10) && close(cephes_hankel(2.0, 80.0), jn(2, 80.0), 1e-10), "Hankel's expansion at large x");

    /* domain error: negative argument with non-integer order */
    check(cephes_bessel_Jv(0.5, -1.0) == 0.0, "negative argument with a non-integer order gives 0 (with an error message)");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
