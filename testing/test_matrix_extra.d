/* Smoke tests for matrix_extra.d: conversions between matrices and arrays/series/datasets, the data-subset
   functions, file reading and writing (text, csv, binary), printing, constructors, row/column cutting and masks,
   quadrature rules, 2-d convolution, which() and vector_from_strings().

   Not exercised: the model-based extractors (need an estimated model), matrix_mp_ols (plugin) and the complex
   print functions (gretl_cmatrix.c is not ported).

   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : _IOLBF, printf, remove, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strstr;

import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.graphingh : FitType;
import gretl.gretl_prnh : PrnType;
import gretl.libgretlh;
import gretl.matrix_extra;
import gretl.matrix_extrah : MMissingCode, QuadMethod;
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

bool close(double a, double b, double tol = 1e-9)
{
    return fabs(a - b) <= tol * (1.0 + fabs(b));
}

GretlMatrixRaw * mat(int r, int c, double[] v...)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(r, c);

    m.val[0 .. v.length] = v[];
    return m;
}

bool same(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b)
{
    if (a.rows != b.rows || a.cols != b.cols) return false;
    foreach (i; 0 .. a.rows * a.cols) {
        if (a.val[i] != b.val[i] && !(a.val[i] != a.val[i] && b.val[i] != b.val[i])) return false;
    }
    return true;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* conversions */
    double[5] arr = [1, double.nan, 3, 4, 5];
    GretlVector * v = gretl_vector_from_array(arr.ptr, 3, GretlMatrixMod.none);

    check(v.rows == 3 && v.val[0] == 1 && v.val[1] == 3 && v.val[2] == 4, "vector_from_array skips missing values");
    gretl_matrix_free(v);
    v = gretl_vector_from_array(arr.ptr, 2, GretlMatrixMod.square);
    check(v.val[0] == 1 && v.val[1] == 9, "vector_from_array with squares");
    gretl_matrix_free(v);
    v = gretl_vector_from_series(arr.ptr, 2, 4);
    check(v.rows == 3 && v.val[0] == 3 && gretl_matrix_get_t1(v) == 2 && gretl_matrix_get_t2(v) == 4, "vector_from_series records t1 and t2");
    gretl_matrix_free(v);
    check(gretl_vector_from_series(arr.ptr, 3, 1) is null, "vector_from_series with an empty range is NULL");

    double[2] c0 = [1, 2];
    double[2] c1 = [3, 4];
    const(double) *[2] cols = [c0.ptr, c1.ptr];
    GretlMatrixRaw * m = gretl_matrix_from_2d_array(cols.ptr, 2, 2);

    check(m.val[0] == 1 && m.val[1] == 2 && m.val[2] == 3 && m.val[3] == 4, "matrix_from_2d_array");
    int err = 0;
    double ** cc = gretl_matrix_get_columns(m, &err);

    check(err == 0 && cc[1] == m.val + 2 && cc[0][1] == 2, "matrix_get_columns points into the matrix");
    free(cc);
    double[2] row;

    check(gretl_matrix_row_to_array(m, 1, row.ptr) == 0 && row[0] == 2 && row[1] == 4, "matrix_row_to_array");
    check(gretl_matrix_row_to_array(m, 2, row.ptr) == GretlError.data, "matrix_row_to_array out of bounds");
    gretl_matrix_free(m);
    m = gretl_matrix_from_scalar(7.5);
    check(m.rows == 1 && m.cols == 1 && m.val[0] == 7.5, "matrix_from_scalar");
    gretl_matrix_free(m);

    /* data subsets */
    DatasetRaw * d = create_new_dataset(3, 5, 0);

    d.pd = 1;
    foreach (i, nm; ["const", "x", "y"]) {
        d.varname[i][0 .. nm.length] = nm[];
        d.varname[i][nm.length] = 0;
    }
    foreach (t; 0 .. 5) {
        d.Z[0][t] = 1;
        d.Z[1][t] = t + 1;
        d.Z[2][t] = 10 * (t + 1);
    }
    d.Z[2][2] = NADBL;
    int[3] list = [2, 1, 2];

    err = 0;
    m = gretl_matrix_data_subset(list.ptr, d, 0, 4, MMissingCode.ok, &err);
    check(err == 0 && m.rows == 5 && m.cols == 2 && na(m.val[7]) && m.val[4] == 5, "data_subset, missing OK");
    const(char) ** cn = gretl_matrix_get_colnames(m);
    check(cn !is null && strcmp(cn[0], "x") == 0 && strcmp(cn[1], "y") == 0, "data_subset attaches column names");
    gretl_matrix_free(m);
    m = gretl_matrix_data_subset(list.ptr, d, 0, 4, MMissingCode.skip, &err);
    check(err == 0 && m.rows == 4 && m.val[2] == 4, "data_subset, skip missing rows");
    gretl_matrix_free(m);
    m = gretl_matrix_data_subset(list.ptr, d, 0, 4, MMissingCode.error, &err);
    check(err == GretlError.missdata && m is null, "data_subset, missing is an error");
    err = 0;
    m = gretl_matrix_data_subset(list.ptr, d, 3, 4, MMissingCode.error, &err);
    check(err == 0 && m.rows == 2 && gretl_matrix_get_t1(m) == 3 && gretl_matrix_get_t2(m) == 4, "data_subset records the sample range");
    gretl_matrix_free(m);
    m = gretl_matrix_data_subset(null, d, -1, -1, MMissingCode.skip, &err);
    check(err == 0 && m.rows == 4 && m.cols == 2, "data_subset of the whole dataset, no list");
    gretl_matrix_free(m);
    char[5] mask = [0, 1, 0, 0, 1];

    m = gretl_matrix_data_subset_masked(list.ptr, d, 0, 4, mask.ptr, &err);
    check(err == GretlError.missdata, "data_subset_masked still fails on a missing value that is not masked");
    err = 0;
    mask[2] = 1;
    m = gretl_matrix_data_subset_masked(list.ptr, d, 0, 4, mask.ptr, &err);
    check(err == 0 && m.rows == 2 && m.val[0] == 1 && m.val[1] == 4, "data_subset_masked");
    gretl_matrix_free(m);
    GretlMatrixRaw * mm = mat(5, 1, 0, 1, 1, 0, 1);

    m = gretl_matrix_data_subset_special(list.ptr, d, mm, &err);
    check(err == 0 && m.rows == 3 && m.cols == 2 && m.val[0] == 2 && m.val[2] == 5, "data_subset_special");
    gretl_matrix_free(m);
    gretl_matrix_free(mm);

    /* datasets from matrices */
    {
        GretlMatrixRaw * X = mat(3, 2, 1, 2, 3, 4, 5, 6);
        DatasetRaw * ds = gretl_dataset_from_matrix(X, null, Gretlopt.OPT_NONE, &err);

        check(err == 0 && ds.v == 3 && ds.n == 3 && ds.Z[2][0] == 4 && ds.Z[2] != X.val + 3 && strcmp(ds.varname[1], "col1") == 0,
              "dataset_from_matrix copies and names the columns");
        destroy_dataset(ds);
        ds = gretl_dataset_from_matrix(X, null, Gretlopt.OPT_B, &err);
        check(ds.Z[2] == X.val + 3, "dataset_from_matrix with OPT_B borrows the data");
        free(ds.Z);
        int[2] cl = [1, 2];

        gretl_dataset_from_matrix(X, cl.ptr, Gretlopt.OPT_NONE, &err);
        check(err == 0, "dataset_from_matrix with a column list");
        cl[1] = 3;
        err = 0;
        gretl_dataset_from_matrix(X, cl.ptr, Gretlopt.OPT_NONE, &err);
        check(err == GretlError.data, "dataset_from_matrix with a bad column number");
        err = 0;
        GretlMatrixRaw * X2 = mat(2, 1, 1, 2);

        check(gretl_dataset_from_matrices(X, X2, &err) is null && err == GretlError.nonconf, "dataset_from_matrices needs equal row counts");
        gretl_matrix_free(X);
        gretl_matrix_free(X2);
        err = 0;
    }

    /* plot-fit matrices */
    {
        double[4] yy = [1, 2, 3, 4];
        double[4] xx = [1, 2, double.nan, 4];
        GretlMatrixRaw * py, pX;

        check(gretl_plotfit_matrices(yy.ptr, xx.ptr, FitType.quadratic, 0, 3, &py, &pX) == 0, "plotfit_matrices");
        check(py.rows == 3 && pX.rows == 3 && pX.cols == 3 && gretl_matrix_get(pX, 2, 2) == 16, "plotfit quadratic design matrix skips missing");
        gretl_matrix_free(py);
        gretl_matrix_free(pX);
        double[2] neg = [-1, 2];

        check(gretl_plotfit_matrices(neg.ptr, null, FitType.loglin, 0, 1, &py, &pX) == GretlError.data, "plotfit loglin rejects non-positive y");
    }

    /* files */
    {
        GretlMatrixRaw * A = mat(2, 3, 1.5, -2, 3, double.nan, 5, 6.25);
        const(char) * fn = "build/test_matrix_extra.mat";
        const(char) * fb = "build/test_matrix_extra.bin";
        const(char) * fc = "build/test_matrix_extra.csv";

        check(gretl_matrix_write_to_file(A, fn, 0) == 0, "write text");
        GretlMatrixRaw * B = gretl_matrix_read_from_file(fn, 0, &err);

        check(err == 0 && B !is null && same(A, B), "text file round trip (incl. NA)");
        gretl_matrix_free(B);
        check(gretl_matrix_write_to_file(A, fb, 0) == 0, "write binary");
        B = gretl_matrix_read_from_file(fb, 0, &err);
        check(err == 0 && B !is null && same(A, B), "binary file round trip");
        gretl_matrix_free(B);
        check(gretl_matrix_write_to_file(A, fc, 0) == 0, "write csv");
        remove(fn);
        remove(fb);
        remove(fc);
        B = gretl_matrix_read_from_file("build/no_such_dir/none.mat", 0, &err);
        check(B is null && err == GretlError.fopen, "reading a missing file is E_FOPEN");
        err = 0;
        gretl_matrix_free(A);
    }

    /* printing and constructors */
    {
        Prn * prn = gretl_print_new(PrnType.buffer, null);
        GretlMatrixRaw * A = mat(2, 2, 1, 2, 3, 4.5);

        gretl_matrix_print_to_prn(A, "A", prn);
        const(char) * buf = gretl_print_get_buffer(prn);

        check(strstr(buf, "A (2 x 2)") !is null && strstr(buf, "4.5") !is null, "print_to_prn: header and values");
        gretl_print_destroy(prn);

        prn = gretl_print_new(PrnType.buffer, null);
        gretl_matrix_print_with_format(A, "%8.3f", 0, 0, prn);
        buf = gretl_print_get_buffer(prn);
        check(strstr(buf, "   4.500") !is null && strstr(buf, "   1.000") !is null, "print_with_format");
        gretl_print_destroy(prn);

        prn = gretl_print_new(PrnType.buffer, null);
        gretl_matrix_print_to_prn(null, "N", prn);
        check(strstr(gretl_print_get_buffer(prn), "matrix is NULL") !is null, "print_to_prn of a NULL matrix");
        gretl_print_destroy(prn);

        char * s = gretl_matrix_write_constructor(A);

        check(strcmp(s, "{1,3;2,4.5}") == 0, "write_constructor");
        free(s);
        GretlMatrixRaw * big = gretl_matrix_alloc(5, 5);

        s = gretl_matrix_write_constructor(big);
        check(strcmp(s, "{}") == 0, "write_constructor of an oversize matrix is {}");
        free(s);
        gretl_matrix_free(big);
        gretl_matrix_free(A);
    }

    /* cutting and masks */
    {
        GretlMatrixRaw * A = mat(3, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9);
        char[3] mk = [0, 1, 0];

        check(gretl_matrix_cut_rows_cols(A, mk.ptr) == 0 && A.rows == 2 && A.cols == 2
              && A.val[0] == 1 && A.val[1] == 3 && A.val[2] == 7 && A.val[3] == 9, "cut_rows_cols");
        gretl_matrix_free(A);
        A = mat(3, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9);
        gretl_matrix_cut_rows(A, mk.ptr);
        check(A.rows == 2 && A.val[0] == 1 && A.val[1] == 3 && A.val[2] == 4 && A.val[3] == 6, "cut_rows");
        gretl_matrix_free(A);
        A = mat(3, 3, 1, 2, 3, 4, 5, 6, 7, 8, 9);
        gretl_matrix_cut_cols(A, mk.ptr);
        check(A.cols == 2 && A.val[2] == 3 && A.val[3] == 7 && A.val[5] == 9, "cut_cols");
        gretl_matrix_free(A);

        A = mat(3, 3, 1, 0, 3, 0, 0, 0, 7, 0, 9);
        char * zr = gretl_matrix_zero_row_mask(A, &err);

        check(zr !is null && zr[0] == 0 && zr[1] == 1 && zr[2] == 0, "zero_row_mask");
        free(zr);
        char * zc = gretl_matrix_zero_col_mask(A, &err);

        check(zc !is null && zc[0] == 0 && zc[1] == 1 && zc[2] == 0, "zero_col_mask");
        free(zc);
        char * zd = gretl_matrix_zero_diag_mask(A, &err);

        check(zd !is null && zd[1] == 1 && zd[0] == 0, "zero_diag_mask");
        free(zd);
        gretl_matrix_free(A);
        A = mat(3, 2, 1, 2, 3, 2, 4, 6);   /* second column = 2 * first */
        char * rk = gretl_matrix_rank_mask(A, &err);

        check(rk !is null && rk[0] == 0 && rk[1] == 1, "rank_mask flags a dependent column");
        free(rk);
        gretl_matrix_free(A);
    }

    /* quadrature */
    {
        GretlMatrixRaw * q = gretl_quadrule_matrix_new(2, QuadMethod.legendre, NADBL, NADBL, &err);

        check(err == 0 && close(q.val[0], -1 / sqrt(3.0)) && close(q.val[1], 1 / sqrt(3.0)) && close(q.val[2], 1) && close(q.val[3], 1),
              "Gauss-Legendre, n = 2");
        gretl_matrix_free(q);
        q = gretl_quadrule_matrix_new(2, QuadMethod.legendre, 0, 2, &err);
        check(close(q.val[0], 1 - 1 / sqrt(3.0)) && close(q.val[2], 1), "Gauss-Legendre rescaled to (0, 2)");
        gretl_matrix_free(q);
        q = gretl_quadrule_matrix_new(2, QuadMethod.ghermite, NADBL, NADBL, &err);
        check(close(q.val[1], 1 / sqrt(2.0)) && close(q.val[2], sqrt(3.141592653589793) / 2), "Gauss-Hermite, n = 2");
        gretl_matrix_free(q);
        q = gretl_quadrule_matrix_new(2, QuadMethod.laguerre, NADBL, NADBL, &err);
        check(close(q.val[0], 2 - sqrt(2.0)) && close(q.val[2] + q.val[3], 1), "Gauss-Laguerre, n = 2");
        gretl_matrix_free(q);
        q = gretl_quadrule_matrix_new(5, QuadMethod.ghermite, 0, 1, &err);
        double sw = 0;
        foreach (i; 0 .. 5) sw += q.val[5 + i];
        check(close(sw, 1), "normal-scaled Gauss-Hermite weights sum to 1");
        gretl_matrix_free(q);
        check(gretl_quadrule_matrix_new(3, 99, 0, 1, &err) is null && err == GretlError.data, "unknown method is E_DATA");
        err = 0;
        q = gretl_gauss_hermite_matrix_new(3, &err);
        check(err == 0 && q.rows == 3, "gauss_hermite_matrix_new");
        gretl_matrix_free(q);
    }

    /* convolution, which, strings */
    {
        GretlMatrixRaw * A = mat(1, 3, 1, 2, 3);
        GretlMatrixRaw * B = mat(1, 2, 1, 1);
        GretlMatrixRaw * C = gretl_matrix_2d_convolution(A, B, &err);

        check(C.rows == 1 && C.cols == 4 && C.val[0] == 1 && C.val[1] == 3 && C.val[2] == 5 && C.val[3] == 3, "2d_convolution");
        gretl_matrix_free(A);
        gretl_matrix_free(B);
        gretl_matrix_free(C);
        A = mat(1, 4, 0, 5, 0, 7);
        C = gretl_matrix_which(A, &err);
        check(C.rows == 2 && C.val[0] == 2 && C.val[1] == 4, "which() of a vector gives 1-based positions");
        gretl_matrix_free(C);
        gretl_matrix_free(A);
        A = mat(2, 2, 0, 1, 1, 0);
        C = gretl_matrix_which(A, &err);
        check(C.rows == 2 && C.cols == 2 && C.val[0] == 2 && C.val[2] == 1 && C.val[3] == 2, "which() of a matrix gives row, column pairs");
        gretl_matrix_free(C);
        gretl_matrix_free(A);

        char[8] s1 = "1.5", s2 = "abc", s3 = "4";
        char *[3] S = [s1.ptr, s2.ptr, s3.ptr];
        int nv;

        C = vector_from_strings(S.ptr, 3, "%m", &nv, &err);
        check(err == 0 && nv == 2 && C.rows == 3 && C.val[0] == 1.5 && na(C.val[1]) && C.val[2] == 4, "vector_from_strings");
        gretl_matrix_free(C);
        C = vector_from_strings(S.ptr, 3, "%2m", &nv, &err);
        check(C.rows == 2, "vector_from_strings with a row limit");
        gretl_matrix_free(C);
        check(vector_from_strings(S.ptr, 3, "%x", &nv, &err) is null && err == GretlError.invarg, "vector_from_strings rejects a bad format");
        err = 0;
    }

    /* covariance from a varlist */
    {
        int[3] cl = [2, 1, 2];
        d.Z[2][2] = 30;
        GretlMatrixRaw * means;
        GretlMatrixRaw * V = gretl_covariance_matrix_from_varlist(cl.ptr, d, &means, &err);

        check(err == 0 && close(V.val[0], 2.5) && close(V.val[1], 25) && close(V.val[3], 250) && close(means.val[1], 30),
              "covariance_matrix_from_varlist");
        gretl_matrix_free(V);
        gretl_matrix_free(means);
    }

    destroy_dataset(d);
    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
