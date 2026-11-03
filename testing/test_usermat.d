/* Smoke tests for usermat.d: matrix lookup by name (plain, stolen, copied), submatrix specifications (range
   selection, column selection, replacement and scalar assignment), vec/vech/unvech, in-place inversion and
   transposition, matrix OLS, the symmetric eigenproblem, and row/column names.

   Not exercised: selection by matrix or "diag"/"upper"/"lower" dummies, complex matrices, QR/SVD, RLS, the GHK
   simulator and the bundle/array name setters.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp;

import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.usermat;
import gretl.usermath : MSEL_MAX, MatrixSubspec, SelType;
import gretl.uservar : user_var_add;

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

GretlMatrixRaw * make(int r, int c, const(double)[] v)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(r, c);
    foreach (i, e; v) m.val[i] = e;
    return m;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* by name */
    GretlMatrixRaw * A = make(4, 3, [1, 2, 3, 4,  5, 6, 7, 8,  9, 10, 11, 12]);  /* column-major */
    check(user_var_add("A", GretlType.matrix, A) == 0 && get_matrix_by_name("A") is A, "get_matrix_by_name finds a user matrix");
    check(get_matrix_by_name("nonesuch") is null && get_matrix_by_name(null) is null && get_matrix_by_name("") is null, "unknown, null and empty names");
    int err = 0;
    GretlMatrixRaw * cp = get_matrix_copy_by_name("A", &err);
    check(err == 0 && cp !is null && cp !is A && cp.val[5] == 6, "get_matrix_copy_by_name returns a copy");
    gretl_matrix_free(cp);

    /* submatrix by range: rows 2:3, all columns */
    MatrixSubspec * spec = matrix_subspec_new();
    spec.ltype = SelType.range;
    spec.lsel.range[0] = 2;
    spec.lsel.range[1] = 3;
    spec.rtype = SelType.all;
    GretlMatrixRaw * S = matrix_get_submatrix(A, spec, 0, &err);
    check(err == 0 && S !is null && S.rows == 2 && S.cols == 3 && S.val[0] == 2 && S.val[1] == 3 && S.val[5] == 11, "rows 2:3 of a 4x3 matrix");
    gretl_matrix_free(S);
    free(spec);

    /* a single element: A[3,2] = 7 */
    spec = matrix_subspec_new();
    spec.ltype = SelType.range;
    spec.lsel.range[0] = 3;
    spec.lsel.range[1] = 3;
    spec.rtype = SelType.range;
    spec.rsel.range[0] = 2;
    spec.rsel.range[1] = 2;
    S = matrix_get_submatrix(A, spec, 0, &err);
    check(err == 0 && S !is null && S.rows == 1 && S.cols == 1 && S.val[0] == 7, "element selection");
    gretl_matrix_free(S);
    free(spec);

    /* replace columns 2:3 of a copy by a 4x2 block of -1 */
    GretlMatrixRaw * B = make(4, 3, [1, 2, 3, 4,  5, 6, 7, 8,  9, 10, 11, 12]);
    GretlMatrixRaw * R = make(4, 2, [-1, -1, -1, -1, -1, -1, -1, -1]);
    spec = matrix_subspec_new();
    spec.ltype = SelType.all;
    spec.rtype = SelType.range;
    spec.rsel.range[0] = 2;
    spec.rsel.range[1] = 3;
    err = matrix_replace_submatrix(B, R, spec);
    check(err == 0 && B.val[3] == 4 && B.val[4] == -1 && B.val[11] == -1, "replace columns 2:3");
    free(spec);

    /* scalar assignment to all of rows 1:2 */
    spec = matrix_subspec_new();
    spec.ltype = SelType.range;
    spec.lsel.range[0] = 1;
    spec.lsel.range[1] = 2;
    spec.rtype = SelType.all;
    err = assign_scalar_to_submatrix(B, null, 0.5, spec);
    check(err == 0 && B.val[0] == 0.5 && B.val[1] == 0.5 && B.val[2] == 3 && B.val[8] == 0.5, "assign_scalar_to_submatrix");
    free(spec);
    gretl_matrix_free(R);
    gretl_matrix_free(B);

    /* vec, vech, unvech */
    GretlMatrixRaw * Q = make(2, 2, [1, 2, 2, 3]);
    GretlMatrixRaw * v = user_matrix_vec(Q, &err);
    check(err == 0 && v !is null && v.rows == 4 && v.val[1] == 2 && v.val[3] == 3, "vec");
    GretlMatrixRaw * vh = user_matrix_vech(Q, 0, &err);
    check(err == 0 && vh !is null && vh.rows == 3 && vh.val[0] == 1 && vh.val[1] == 2 && vh.val[2] == 3, "vech");
    GretlMatrixRaw * uv = user_matrix_unvech(vh, NADBL, &err);
    check(err == 0 && uv !is null && uv.rows == 2 && uv.val[1] == 2 && uv.val[2] == 2 && uv.val[3] == 3, "unvech (NA diagonal) restores the symmetric matrix");
    gretl_matrix_free(v);
    gretl_matrix_free(vh);
    gretl_matrix_free(uv);

    /* in place */
    err = matrix_transpose_in_place(A);
    check(err == 0 && A.rows == 3 && A.cols == 4 && A.val[1] == 5 && A.val[3] == 2, "transpose in place");
    GretlMatrixRaw * I2 = make(2, 2, [2, 0, 0, 4]);
    err = matrix_invert_in_place(I2);
    check(err == 0 && close(I2.val[0], 0.5) && close(I2.val[3], 0.25), "invert in place");
    gretl_matrix_free(I2);

    /* OLS: y = 1 + 2x exactly */
    GretlMatrixRaw * Y = make(4, 1, [3, 5, 7, 9]);
    GretlMatrixRaw * X = make(4, 2, [1, 1, 1, 1, 1, 2, 3, 4]);
    GretlMatrixRaw * b = user_matrix_ols(Y, X, null, null, null, null, Gretlopt.OPT_NONE, &err);
    check(err == 0 && b !is null && close(b.val[0], 1.0, 1e-8) && close(b.val[1], 2.0, 1e-8), "user_matrix_ols recovers the coefficients");
    gretl_matrix_free(b);
    gretl_matrix_free(Y);
    gretl_matrix_free(X);

    /* symmetric eigenvalues, ascending */
    GretlMatrixRaw * E = user_matrix_eigensym(Q, null, &err);
    if (E !is null) {
        double tr = E.val[0] + E.val[1];
        check(close(tr, 4.0) && close(E.val[0] * E.val[1], 3.0 - 4.0), "eigenvalues of [[1,2],[2,3]]: trace 4, determinant -1");
    }
    gretl_matrix_free(E);
    gretl_matrix_free(Q);

    /* names */
    GretlMatrixRaw * N = make(2, 2, [1, 2, 3, 4]);
    err = umatrix_set_names_from_string(N, "r1 r2", 1, null);
    char * rn = user_matrix_get_row_name(N, 2, &err);
    check(err == 0 && rn !is null && strcmp(rn, "r2") == 0, "row names from a string");
    free(rn);
    err = 0;
    err = umatrix_set_names_from_string(N, "c1 c2", 0, null);
    char * cn = user_matrix_get_column_name(N, 1, &err);
    check(err == 0 && cn !is null && strcmp(cn, "c1") == 0, "column names from a string");
    free(cn);
    gretl_matrix_free(N);

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}

