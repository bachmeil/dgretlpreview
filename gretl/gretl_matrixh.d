/* Port of gretl_matrix.h */

module gretl.gretl_matrixh;


import gretl.gretl_matrix : MatrixInfo, gretl_matrix_alloc, gretl_matrix_free;
import gretl.libgretlh;

/* minimum value of diagonal element of R (as in X = QR) that counts
   as non-zero for the purpose of determining the rank of X */

enum R_DIAG_MIN = 1.0e-8;

enum GretlMatrixMod {
    none = 0,
    transpose,
    square,
    cumulate,
    decrement,
    ctransp
}

enum GretlMatrixStructure {
    square = 1,
    lower_triangular,
    upper_triangular,
    symmetric,
    hermitian,
    diagonal,
    identity,
    scalar,
}

enum GretlVecStat {
    sum,
    prod,
    mean
}

enum ConfType {
    none = 0,
    elements,
    a_colvec,
    b_colvec,
    a_rowvec,
    b_rowvec,
    a_scalar,
    b_scalar,
    ac_br,
    ar_bc
}

/* The basic libgretl matrix type; GretlVector is an alias
   that can be used for matrices with rows or cols = 1. */

struct GretlMatrixRaw {
    int rows;
    int cols;
    double * val;
    Cmplx * z;          /* was "complex" */
    int is_complex;
    /* private */
    MatrixInfo * info;
}

alias GretlVector = GretlMatrixRaw;

double gretl_matrix_get(const(GretlMatrixRaw) * m, int i, int j)
{
    return m.val[j * m.rows + i];
}

Cmplx gretl_cmatrix_get(const(GretlMatrixRaw) * m, int i, int j)
{
    return m.z[j * m.rows + i];
}

double gretl_vector_get(const(GretlMatrixRaw) * v, int i)
{
    return v.val[i];
}

void gretl_matrix_set(GretlMatrixRaw * m, int i, int j, double x)
{
    m.val[j * m.rows + i] = x;
}

void gretl_cmatrix_set(GretlMatrixRaw * m, int i, int j, Cmplx x)
{
    m.z[j * m.rows + i] = x;
}

void gretl_vector_set(GretlMatrixRaw * v, int i, double x)
{
    v.val[i] = x;
}

int gretl_matrix_cols(const(GretlMatrixRaw) * m)
{
    return (m == null) ? 0 : m.cols;
}

int gretl_matrix_rows(const(GretlMatrixRaw) * m)
{
    return (m == null) ? 0 : m.rows;
}

int gretl_vector_get_length(const(GretlMatrixRaw) * v)
{
    return (v == null) ? 0 :
        (v.cols == 1) ? v.rows :
        (v.rows == 1) ? v.cols : 0;
}

GretlMatrixRaw * gretl_vector_alloc(int i)
{
    return gretl_matrix_alloc(1, i);
}

GretlMatrixRaw * gretl_column_vector_alloc(int i)
{
    return gretl_matrix_alloc(i, 1);
}

void gretl_vector_free(GretlMatrixRaw * v)
{
    gretl_matrix_free(v);
}

bool gretl_matrix_is_scalar(const(GretlMatrixRaw) * m)
{
    return (m != null &&
            m.is_complex == 0 &&
            m.rows == 1 &&
            m.cols == 1);
}

bool gretl_matrix_is_cscalar(const(GretlMatrixRaw) * m)
{
    return (m != null &&
            m.is_complex &&
            m.rows == 1 &&
            m.cols == 1);
}

bool gretl_is_null_matrix(const(GretlMatrixRaw) * m)
{
    return (m == null || m.rows == 0 || m.cols == 0);
}

bool gretl_is_complex(const(GretlMatrixRaw) * m)
{
    return (m != null && m.is_complex == 1);
}
