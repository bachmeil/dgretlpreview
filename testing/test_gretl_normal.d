/* Smoke tests for gretl_normal.d: the inverse Mills ratio, the bivariate normal CDF (against closed forms), and
   the GHK simulator for multivariate normal probabilities (exact in one dimension, close to the product of the
   marginals for an identity covariance with a stratified set of uniform draws, error handling).

   Not exercised: the GHK score matrix (gretl_GHK2).

   Build and run with testing/run_tests.sh. */

import cephes.ndtr : ndtr;
import core.stdc.math : asin, fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_normal;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

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

bool close(double a, double b, double tol = 1e-8)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}


int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    check(close(invmills(0), 0.7978845608028654), "invmills(0) = sqrt(2/pi)");
    check(close(invmills(1), 1.5251352761609812, 1e-7), "invmills(1)");
    check(close(invmills(-1), 0.28759997093917836, 1e-7), "invmills(-1)");
    check(close(invmills(30), 30.0332, 1e-3) && invmills(-30) == 0, "invmills in the far tails");

    check(close(bvnorm_cdf(0.5, 0, 0), 1.0 / 3.0), "bvnorm_cdf(0.5, 0, 0) = 1/3");
    check(close(bvnorm_cdf(0, 0.3, -0.2), ndtr(0.3) * ndtr(-0.2)), "independent case is the product of the marginals");
    check(close(bvnorm_cdf(1, 0.3, -0.2), ndtr(-0.2)), "rho = 1: the smaller abscissa");
    check(close(bvnorm_cdf(-1, 1, -0.5), ndtr(1) - ndtr(-0.5)), "rho = -1");
    check(na(bvnorm_cdf(1.5, 0, 0)), "|rho| > 1 is NA");
    check(close(bvnorm_cdf(-0.7, 0.4, 1.1), bvnorm_cdf(-0.7, 1.1, 0.4)), "symmetry in the arguments");
    check(close(bvnorm_cdf(0.8, 40, 0.3), ndtr(0.3), 1e-6), "a very large abscissa leaves the other marginal");

    /* one dimension: P(-1 < 2z < 1) */
    int err = 0;
    GretlMatrixRaw * C = gretl_matrix_alloc(1, 1);
    GretlMatrixRaw * A = gretl_matrix_alloc(1, 1);
    GretlMatrixRaw * B = gretl_matrix_alloc(1, 1);
    GretlMatrixRaw * U = gretl_matrix_alloc(1, 10);

    C.val[0] = 2;
    A.val[0] = -1;
    B.val[0] = 1;
    foreach (i; 0 .. 10) U.val[i] = (i + 0.5) / 10;
    GretlMatrixRaw * P = gretl_GHK(C, A, B, U, &err);

    check(err == 0 && P !is null && P.rows == 1 && close(P.val[0], ndtr(0.5) - ndtr(-0.5), 1e-6), "GHK in one dimension is exact");
    gretl_matrix_free(P);

    /* bad bounds */
    B.val[0] = -2;
    P = gretl_GHK(C, A, B, U, &err);
    check(P is null && err == GretlError.data, "GHK with B < A is an error");
    err = 0;

    /* NaN bounds give NaN */
    A.val[0] = double.nan;
    B.val[0] = 1;
    P = gretl_GHK(C, A, B, U, &err);
    check(err == 0 && P !is null && P.val[0] != P.val[0], "GHK with a missing bound gives NaN");
    gretl_matrix_free(P);
    gretl_matrix_free(C);
    gretl_matrix_free(A);
    gretl_matrix_free(B);
    gretl_matrix_free(U);

    /* two independent dimensions (C = I), draws on a Halton-like grid */
    C = gretl_matrix_alloc(2, 2);
    A = gretl_matrix_alloc(1, 2);
    B = gretl_matrix_alloc(1, 2);
    U = gretl_matrix_alloc(2, 400);
    C.val[0] = 1; C.val[1] = 0; C.val[2] = 0; C.val[3] = 1;
    A.val[0] = -1e10; A.val[1] = -1e10;
    B.val[0] = 0.3; B.val[1] = -0.4;
    foreach (j; 0 .. 400) {
        U.val[2 * j] = (j + 0.5) / 400;
        double v = 0, f = 1.0 / 3;
        int k = j + 1;
        while (k > 0) { v += f * (k % 3); k /= 3; f /= 3; }
        U.val[2 * j + 1] = v;
    }
    P = gretl_GHK(C, A, B, U, &err);
    check(err == 0 && fabs(P.val[0] - ndtr(0.3) * ndtr(-0.4)) < 0.01, "GHK in two independent dimensions is near the product");
    gretl_matrix_free(P);
    gretl_matrix_free(C);
    gretl_matrix_free(A);
    gretl_matrix_free(B);
    gretl_matrix_free(U);

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
