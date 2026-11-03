/* Tests for cephes/struve.d: the Struve function H_v(x) against numerical integration of its integral
   representation, Y_v(x) for half-integer orders (closed forms) and integer orders (libm's yn, through
   cephes_bessel_Yn), the negative half-integer case of struve (equal to a Bessel J), and the series
   helpers onef2 and threef0 on a case with a closed form. A bit-for-bit comparison with the original C code over
   75,180 points was run once when the module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : cos, exp, fabs, pow, sin, sqrt, tgamma;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.constants : PI;
import cephes.jv : cephes_bessel_Jv;
import cephes.struve;

extern (C) double yn(int n, double x);   // libm

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

bool close(double a, double b, double tol = 1e-10)
{
    return fabs(a - b) <= tol * (fabs(b) + 1e-300);
}

/* H_v(x) = 2 (x/2)^v / (sqrt(pi) Gamma(v + 1/2)) * integral_0^{pi/2} sin(x cos t) sin(t)^{2v} dt (v > -1/2),
   by Simpson's rule */
double struve_integral(double v, double x)
{
    enum N = 20000;
    double h = (PI / 2) / N;
    double s = 0;
    foreach (i; 0 .. N + 1) {
        double t = i * h;
        double f = sin(x * cos(t)) * pow(sin(t), 2 * v);
        double w = (i == 0 || i == N) ? 1 : (i % 2 ? 4 : 2);
        s += w * f;
    }
    s *= h / 3;
    return 2 * pow(x / 2, v) / (sqrt(PI) * tgamma(v + 0.5)) * s;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* H_v against the integral, series region (x < 18) and asymptotic region (x >= 18) */
    bool ok = true;
    double worst = 0;
    foreach (v; [0.0, 1.0, 2.0, 0.5, 1.5, 3.3]) {
        foreach (x; [0.3, 1.0, 4.0, 9.0, 17.0, 25.0, 40.0]) {
            if (v == 0.5 && x > 20.0) continue;   // see the next check (bugs.md #106)
            double e = fabs(struve(v, x) - struve_integral(v, x)) / (1e-3 + fabs(struve_integral(v, x)));
            if (e > worst) worst = e;
        }
    }
    printf("  worst scaled error of struve against the integral: %g\n", worst);   // ~2e-8, at x = 25 (series cancellation)
    check(worst < 1e-7, "H_v(x) agrees with its integral representation (v from 0 to 3.3, x from 0.3 to 40)");

    /* the failure of the C code at order 1/2 (bugs.md #106): for x > 30 it returns 0 instead of
       sqrt(2/(pi x)) (1 - cos x) = 0.2103 at x = 40, and at x = 25 is off by 1.4e-8 */
    check(struve(0.5, 40.0) == 0.0 && struve_integral(0.5, 40.0) > 0.2 && fabs(struve(0.5, 25.0) - struve_integral(0.5, 25.0)) < 1e-6, "known upstream failure: H_{1/2}(40) = 0, see bugs.md #106");

    /* tabulated: H_0(1) = 0.5686566270482879, H_1(1) = 0.19845699 */
    check(close(struve(0.0, 1.0), 0.5686566270482879, 1e-9), "H_0(1)");

    /* negative half-integer order: H_{-n-1/2}(x) = (-1)^n J_{n+1/2}(x) */
    check(close(struve(-0.5, 3.0), cephes_bessel_Jv(0.5, 3.0), 1e-12) && close(struve(-1.5, 3.0), -cephes_bessel_Jv(1.5, 3.0), 1e-12), "negative half-integer orders reduce to Bessel J");

    /* Y_v: half-integer closed forms, integer orders through cephes_bessel_Yn */
    ok = true;
    foreach (x; [0.5, 2.0, 7.0, 20.0]) {
        double c = sqrt(2 / (PI * x));
        ok = ok && close(cephes_bessel_Yv(0.5, x), -c * cos(x), 1e-11);
        ok = ok && close(cephes_bessel_Yv(1.5, x), -c * (cos(x) / x + sin(x)), 1e-11);
        ok = ok && close(cephes_bessel_Yv(2.0, x), yn(2, x), 1e-11);
        ok = ok && close(cephes_bessel_Yv(0.0, x), yn(0, x), 1e-11);
    }
    check(ok, "Y_v: half-integer closed forms and integer orders");

    /* onef2(a, b, c; x) with a = c: 1F2(a; b, a; x) = 0F1(; b; x), so 1F2(1; 1.5, 1; -t) = cos(2 sqrt(t)) *
       Gamma-free closed form for b = 1/2: 0F1(; 1/2; -z^2/4) = cos z; use b = 0.5, c = 1, a = 1 */
    double err = 0;
    double z = 2.0;
    double v = onef2(1.0, 0.5, 1.0, -z * z / 4, &err);
    check(close(v, cos(z), 1e-12) && err < 1e-10, "1F2(1; 1/2, 1; -z^2/4) = cos z, with a small error estimate");

    /* threef0 is an asymptotic series: with a small argument the first terms are accurate */
    double e3 = 0;
    double r3 = threef0(1.0, 0.5, 0.5, -0.001, &e3);
    check(fabs(r3 - (1 - 0.25 * 0.001 + 0)) < 1e-6, "3F0 with a small argument starts as 1 + (abc) x");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
