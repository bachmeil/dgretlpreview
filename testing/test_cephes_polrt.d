/* Tests for cephes/polrt.d: the complex roots of a real polynomial: real roots, a complex pair, a repeated
   root, zero roots, a higher degree polynomial built from known roots (checked by evaluating the polynomial at
   the returned roots), the return codes, and the use by gretl_matrix_polroots. A bit-for-bit comparison with the
   original C code over 4,000 random polynomials (roots and the modified coefficient array) was run once when the
   module was ported (see log.md).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : cos, fabs, sin, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import cephes.mconfh : Cmplx;
import cephes.polrt;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free, gretl_matrix_polroots;
import gretl.gretl_matrixh : GretlMatrixRaw;

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

/* |p(z)| for the polynomial with ascending coefficients c[0 .. m], by Horner's rule in complex arithmetic */
double absval(const(double)[] c, int m, Cmplx z)
{
    double re = c[m], im = 0;
    for (int k = m - 1; k >= 0; k--) {
        double nre = re * z.r - im * z.i + c[k];
        double nim = re * z.i + im * z.r;
        re = nre;
        im = nim;
    }
    return sqrt(re * re + im * im);
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    Cmplx[40] root;
    double[40] work;

    /* (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6 */
    double[4] c3 = [-6, 11, -6, 1];
    double[4] c3save = c3;
    int rc = polrt(c3.ptr, work.ptr, 3, root.ptr);
    bool ok = rc == 0;
    int found = 0;
    foreach (r; [1.0, 2.0, 3.0]) {
        foreach (i; 0 .. 3) if (fabs(root[i].r - r) < 1e-6 && root[i].i == 0) found++;
    }
    check(ok && found == 3, "real roots 1, 2, 3");

    /* x^2 + 1 */
    double[3] c2 = [1, 0, 1];
    rc = polrt(c2.ptr, work.ptr, 2, root.ptr);
    check(rc == 0 && fabs(root[0].r) < 1e-6 && fabs(fabs(root[0].i) - 1) < 1e-6 && root[0].i == -root[1].i && root[0].r == root[1].r, "complex pair +-i");

    /* (x - 2)^2 (x + 3) = x^3 - x^2 - 8x + 12: a repeated root, found to Newton's lower accuracy */
    double[4] cr = [12, -8, -1, 1];
    double[4] crsave = cr;
    rc = polrt(cr.ptr, work.ptr, 3, root.ptr);
    ok = rc == 0;
    foreach (i; 0 .. 3) ok = ok && absval(crsave[], 3, root[i]) < 1e-3;
    check(ok, "a repeated root: the roots nearly satisfy the polynomial");

    /* zero roots: x^2 (x - 1) = x^3 - x^2 */
    double[4] cz = [0, 0, -1, 1];
    rc = polrt(cz.ptr, work.ptr, 3, root.ptr);
    int zeros = 0, ones = 0;
    foreach (i; 0 .. 3) {
        if (root[i].r == 0 && root[i].i == 0) zeros++;
        if (fabs(root[i].r - 1) < 1e-6 && root[i].i == 0) ones++;
    }
    check(rc == 0 && zeros == 2 && ones == 1, "zero roots");

    /* degree 8, built from known roots: 1, -2, 0.5, 3 and the pairs -1 +- i, 2 +- 0.5i */
    double[9] p = 0;
    p[0] = 1;
    int deg = 0;
    void mulLinear(double a0, double a1)   // multiply by (a1 x + a0)
    {
        foreach_reverse (k; 0 .. deg + 2) {
            p[k] = a0 * (k <= deg ? p[k] : 0) + (k > 0 ? a1 * p[k - 1] : 0);
        }
        deg++;
    }
    foreach (r; [1.0, -2.0, 0.5, 3.0]) mulLinear(-r, 1);
    // (x - a)^2 + b^2 = x^2 - 2a x + (a^2 + b^2): multiply by it as two real factors
    foreach (ab; [[-1.0, 1.0], [2.0, 0.5]]) {
        double a = ab[0], b = ab[1];
        double[9] q = 0;
        foreach (k; 0 .. deg + 1) {
            q[k] += (a * a + b * b) * p[k];
            q[k + 1] += -2 * a * p[k];
            q[k + 2] += p[k];
        }
        deg += 2;
        p = q;
    }
    double[9] psave = p;
    rc = polrt(p.ptr, work.ptr, 8, root.ptr);
    ok = rc == 0 && deg == 8;
    double worst = 0;
    foreach (i; 0 .. 8) {
        double v = absval(psave[], 8, root[i]);
        if (v > worst) worst = v;
    }
    printf("  worst |p(root)| for the degree 8 polynomial: %g\n", worst);
    check(ok && worst < 1e-4, "degree 8 polynomial: the roots satisfy it");

    /* return codes */
    double[2] bad = [1, 1];
    check(polrt(bad.ptr, work.ptr, 0, root.ptr) == 1, "degree 0 returns 1");
    double[40] big = 1;
    check(polrt(big.ptr, work.ptr, 37, root.ptr) == 2, "degree above 36 returns 2");
    double[3] lead0 = [1, 2, 0];
    check(polrt(lead0.ptr, work.ptr, 2, root.ptr) == 4, "zero leading coefficient returns 4");

    /* through gretl_matrix_polroots: roots of 6 - 5x + x^2 are 2 and 3 */
    int err = 0;
    GretlMatrixRaw * a = gretl_matrix_alloc(3, 1);
    a.val[0] = 6; a.val[1] = -5; a.val[2] = 1;
    GretlMatrixRaw * r = gretl_matrix_polroots(a, 0, 0, &err);
    ok = err == 0 && r !is null && r.rows == 2;
    if (ok) {
        double r0 = r.val[0], r1 = r.val[1];
        ok = (fabs(r0 - 2) < 1e-6 && fabs(r1 - 3) < 1e-6) || (fabs(r0 - 3) < 1e-6 && fabs(r1 - 2) < 1e-6);
    }
    check(ok, "gretl_matrix_polroots finds 2 and 3");
    gretl_matrix_free(r);
    gretl_matrix_free(a);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
