/* Smoke tests for gretl_cmatrix.d: FFTs (real and complex), complex multiplication, inverse, rank, determinant,
   trace, transposes, add/sub, build/extract, elementwise operators, Kronecker/horizontal products, Cholesky, QR,
   eigenvalues, Schur, matrix exp/log, printing and the small helpers.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.string : strstr;

import gretl.gretl_cmatrix;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix : gretl_cmatrix_new, gretl_cmatrix_new0, gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;

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

bool close(double a, double b, double tol = 1e-9)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

bool zclose(Cmplx a, Cmplx b, double tol = 1e-9)
{
    return close(a.r, b.r, tol) && close(a.i, b.i, tol);
}

GretlMatrixRaw * cm(int r, int c, Cmplx[] v...)
{
    GretlMatrixRaw * m = gretl_cmatrix_new(r, c);

    m.z[0 .. v.length] = v[];
    return m;
}

GretlMatrixRaw * rm(int r, int c, double[] v...)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(r, c);

    m.val[0 .. v.length] = v[];
    return m;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);
    int err = 0;

    /* FFT of a real vector, and its inverse */
    GretlMatrixRaw * y = rm(4, 1, 1, 2, 3, 4);
    GretlMatrixRaw * f = gretl_matrix_fft(y, &err);

    check(err == 0 && f.is_complex && zclose(f.z[0], Cmplx(10, 0)) && zclose(f.z[1], Cmplx(-2, 2)) && zclose(f.z[2], Cmplx(-2, 0))
          && zclose(f.z[3], Cmplx(-2, -2)), "fft of 1..4");
    GretlMatrixRaw * b = gretl_matrix_ffti(f, &err);

    check(err == 0 && !b.is_complex && close(b.val[0], 1) && close(b.val[3], 4), "ffti of a Hermitian spectrum is the real input");
    gretl_matrix_free(b);
    GretlMatrixRaw * c2 = cm(2, 1, Cmplx(1, 1), Cmplx(0, 2));
    GretlMatrixRaw * c3 = gretl_cmatrix_fft(c2, 0, &err);

    check(err == 0 && zclose(c3.z[0], Cmplx(1, 3)) && zclose(c3.z[1], Cmplx(1, -1)), "complex fft");
    gretl_matrix_free(c3);
    gretl_matrix_free(c2);
    gretl_matrix_free(f);
    gretl_matrix_free(y);

    /* products, transposes, sums */
    GretlMatrixRaw * A = cm(2, 2, Cmplx(1, 1), Cmplx(0, 1), Cmplx(2, 0), Cmplx(1, -1));
    GretlMatrixRaw * B = cm(2, 2, Cmplx(1, 0), Cmplx(0, 0), Cmplx(0, 0), Cmplx(1, 0));
    GretlMatrixRaw * C = gretl_cmatrix_multiply(A, B, &err);

    check(err == 0 && zclose(C.z[0], A.z[0]) && zclose(C.z[3], A.z[3]), "multiply by the identity");
    gretl_matrix_free(C);
    GretlMatrixRaw * R = rm(2, 2, 1, 0, 0, 1);

    C = gretl_cmatrix_multiply(A, R, &err);
    check(err == 0 && zclose(C.z[1], A.z[1]), "multiply complex by real");
    gretl_matrix_free(C);
    C = gretl_cmatrix_AHB(A, B, &err);
    check(zclose(C.z[0], Cmplx(1, -1)) && zclose(C.z[2], Cmplx(0, -1)), "A^H B");
    gretl_matrix_free(C);
    C = gretl_ctrans(A, 1, &err);
    check(zclose(C.z[1], Cmplx(2, 0)) && zclose(C.z[2], Cmplx(0, -1)) && zclose(C.z[0], Cmplx(1, -1)), "conjugate transpose");
    gretl_matrix_free(C);
    C = gretl_cmatrix_add_sub(A, B, -1, &err);
    check(err == 0 && zclose(C.z[0], Cmplx(0, 1)) && zclose(C.z[1], Cmplx(0, 1)), "subtraction");
    gretl_matrix_free(C);
    GretlMatrixRaw * T = gretl_cmatrix_trace(A, &err);

    check(zclose(T.z[0], Cmplx(2, 0)), "trace");
    gretl_matrix_free(T);
    GretlMatrixRaw * Dt = gretl_cmatrix_determinant(A, 0, &err);

    /* det = (1+i)(1-i) - 2i = 2 - 2i */
    check(err == 0 && zclose(Dt.z[0], Cmplx(2, -2), 1e-8), "determinant via eigenvalues");
    gretl_matrix_free(Dt);

    /* inverse, rank */
    GretlMatrixRaw * Ai = gretl_cmatrix_inverse(A, &err);
    GretlMatrixRaw * P = gretl_cmatrix_multiply(A, Ai, &err);

    check(err == 0 && zclose(P.z[0], Cmplx(1, 0), 1e-8) && zclose(P.z[1], Cmplx(0, 0), 1e-8) && zclose(P.z[3], Cmplx(1, 0), 1e-8), "inverse");
    gretl_matrix_free(P);
    gretl_matrix_free(Ai);
    check(gretl_cmatrix_rank(A, &err) == 2, "rank of a nonsingular matrix");
    GretlMatrixRaw * S1 = cm(2, 2, Cmplx(1, 0), Cmplx(2, 0), Cmplx(1, 0), Cmplx(2, 0));

    check(gretl_cmatrix_rank(S1, &err) == 1, "rank of a rank-1 matrix");
    GretlMatrixRaw * Si = gretl_cmatrix_inverse(S1, &err);

    check(err == GretlError.singular && Si is null, "inverse of a singular matrix is an error");
    err = 0;
    gretl_matrix_free(S1);

    /* build, extract, parts, fill, diagonal */
    GretlMatrixRaw * re = rm(1, 2, 1, 2);
    GretlMatrixRaw * im = rm(1, 2, 3, 4);
    GretlMatrixRaw * Z = gretl_cmatrix_build(re, im, 0, 0, &err);

    check(err == 0 && zclose(Z.z[1], Cmplx(2, 4)), "build from real and imaginary parts");
    GretlMatrixRaw * X = gretl_cmatrix_extract(Z, 1, &err);

    check(X.val[0] == 3 && X.val[1] == 4, "extract the imaginary part");
    gretl_matrix_free(X);
    check(gretl_cmatrix_set_part(Z, null, 9, 0) == 0 && zclose(Z.z[0], Cmplx(9, 3)), "set the real part");
    check(gretl_cmatrix_fill(Z, Cmplx(1, 2)) == 0 && zclose(Z.z[1], Cmplx(1, 2)), "fill");
    gretl_matrix_free(Z);
    gretl_matrix_free(re);
    gretl_matrix_free(im);
    GretlMatrixRaw * E = gretl_cmatrix_new0(3, 3);

    check(gretl_cmatrix_set_diagonal(E, null, 5) == 0 && zclose(E.z[4], Cmplx(5, 0)) && zclose(E.z[1], Cmplx(0, 0)), "set_diagonal");
    gretl_matrix_free(E);
    check(matrix_is_complex(A) == 2 && matrix_is_complex(R) == 0, "matrix_is_complex");

    /* elementwise operators and statistics */
    GretlMatrixRaw * D = gretl_cmatrix_dot_op(A, A, '*', &err);

    check(err == 0 && zclose(D.z[0], Cmplx(0, 2)), "dot_op *");
    gretl_matrix_free(D);
    D = gretl_cmatrix_dot_op(A, A, '=', &err);
    check(zclose(D.z[1], Cmplx(1, 0)), "dot_op == gives 1");
    gretl_matrix_free(D);
    GretlMatrixRaw * V = gretl_cmatrix_vector_stat(A, GretlVecStat.sum, 0, 0, &err);

    check(err == 0 && zclose(V.z[0], Cmplx(1, 2)) && zclose(V.z[1], Cmplx(3, -1)), "column sums");
    gretl_matrix_free(V);

    /* Kronecker and horizontal product */
    GretlMatrixRaw * K = gretl_cmatrix_kronecker(A, B, &err);

    check(err == 0 && K.rows == 4 && K.cols == 4 && zclose(K.z[0], A.z[0]) && zclose(K.z[5], A.z[0]), "kronecker");
    gretl_matrix_free(K);
    K = gretl_cmatrix_hdprod(A, null, &err);
    check(err == 0 && K.rows == 2 && K.cols == 3, "hdprod of A with itself is r x p(p+1)/2");
    gretl_matrix_free(K);

    /* Cholesky, QR, eigenvalues */
    GretlMatrixRaw * H = cm(2, 2, Cmplx(2, 0), Cmplx(0, 1), Cmplx(0, -1), Cmplx(2, 0));
    GretlMatrixRaw * L = gretl_cmatrix_cholesky(H, &err);

    check(err == 0 && zclose(L.z[0], Cmplx(sqrt(2.0), 0)) && zclose(L.z[2], Cmplx(0, 0)), "cholesky of a Hermitian matrix");
    gretl_matrix_free(L);
    GretlMatrixRaw * Hc = gretl_matrix_alloc(2, 1);

    gretl_matrix_free(Hc);
    GretlMatrixRaw * Q = gretl_cmatrix_QR_decomp(A, null, &err);
    GretlMatrixRaw * QQ = gretl_cmatrix_AHB(Q, Q, &err);

    check(err == 0 && zclose(QQ.z[0], Cmplx(1, 0), 1e-8) && zclose(QQ.z[1], Cmplx(0, 0), 1e-8), "QR: Q has orthonormal columns");
    gretl_matrix_free(QQ);
    gretl_matrix_free(Q);
    GretlMatrixRaw * w = gretl_zgeev(H, null, null, &err);

    check(err == 0 && w.rows == 2 && close(w.z[0].r + w.z[1].r, 4), "zgeev eigenvalues of a Hermitian matrix sum to the trace");
    gretl_matrix_free(w);
    GretlMatrixRaw * ev = gretl_zheev(H, 0, &err);

    check(err == 0 && close(ev.val[0], 1) && close(ev.val[1], 3), "zheev eigenvalues");
    gretl_matrix_free(ev);
    gretl_matrix_free(H);

    /* matrix exponential and log of a diagonal matrix */
    GretlMatrixRaw * Dg = cm(2, 2, Cmplx(1, 0), Cmplx(0, 0), Cmplx(0, 0), Cmplx(2, 0));
    GretlMatrixRaw * Ex = gretl_cmatrix_exp(Dg, &err);

    check(err == 0 && close(Ex.z[0].r, 2.718281828459045, 1e-8) && close(Ex.z[3].r, 7.38905609893065, 1e-8), "matrix exp of a diagonal matrix");
    GretlMatrixRaw * Lg = gretl_matrix_log(Ex, &err);

    check(err == 0 && close(Lg.val[0], 1, 1e-8) && close(Lg.val[3], 2, 1e-8), "matrix log undoes exp (real result)");
    gretl_matrix_free(Lg);
    gretl_matrix_free(Ex);
    gretl_matrix_free(Dg);

    /* printing and errors */
    Prn * prn = gretl_print_new(PrnType.buffer, null);

    gretl_cmatrix_print(A, "A", prn);
    check(strstr(gretl_print_get_buffer(prn), "A (2 x 2)") !is null && strstr(gretl_print_get_buffer(prn), "i") !is null, "cmatrix_print");
    gretl_print_destroy(prn);
    prn = gretl_print_new(PrnType.buffer, null);
    gretl_cmatrix_printf(A, "%6.2f", prn);
    check(strstr(gretl_print_get_buffer(prn), "  1.00+ 1.00i") !is null || strstr(gretl_print_get_buffer(prn), "1.00") !is null, "cmatrix_printf");
    check(gretl_cmatrix_printf(A, "%d", prn) == GretlError.invarg, "cmatrix_printf accepts only floating-point formats");
    gretl_print_destroy(prn);
    check(gretl_cmatrix_rank(R, &err) == 0 && err == GretlError.invarg, "a real matrix is rejected by cmatrix_validate");
    err = 0;

    gretl_matrix_free(R);
    gretl_matrix_free(B);
    gretl_matrix_free(A);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
