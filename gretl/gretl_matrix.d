/* Port of gretl_matrix.c */

module gretl.gretl_matrix;


/*
 * Libgretl implements most of the matrix functionality that is
 * likely to be required in econometric calculation.  For basics
 * such as decomposition and inversion we use LAPACK as the
 * underlying engine.
 *
 * To get yourself a gretl matrix, use gretl_matrix_alloc() or
 * one of the more specialized constructors; to free such a
 * matrix use gretl_matrix_free().
 *
 * The SIMD (USE_SIMD) and OpenMP code paths of the C source are not
 * part of this port: only the plain serial code is.
 */

import cephes.mconfh : CephesCmplx = Cmplx;
import cephes.polrt : polrt;
import core.stdc.errno : errno;
import core.stdc.fenv : FE_ALL_EXCEPT, FE_DIVBYZERO, FE_INVALID, FE_OVERFLOW, feclearexcept, fetestexcept;
import core.stdc.stdarg : va_arg, va_end, va_list, va_start;
import core.stdc.stdio : FILE, fprintf, fputc, fputs, stderr, stdout;
import core.stdc.stdlib : atoi, calloc, free, getenv, malloc, qsort, realloc;
import std.math : isFinite, isNaN;
import core.stdc.string : memcpy, memset;
import std.algorithm.comparison : max, min;
import core.stdc.math : fabs, floor, log, log2, pow, sqrt;
import gretl.pvaluesh;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new_with_stream;
import gretl.clapack_complexh;
import gretl.clapack_doubleh;
import gretl.gretl_f2ch;
import gretl.glib;
import gretl.gretl_errorsh;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.libset : libset_get_bool;
import gretl.pvalues : x_factorial;
import gretl.libseth : SetKey;
import gretl.missingh;
import gretl.princomp_extra : matrix_is_incomplete, princomp_with_NAs;
import gretl.gretl_list : gretl_list_new, in_gretl_list;
import gretl.strutils : gretl_strdup, strings_array_dup, strings_array_free, strings_array_new, strings_array_reverse;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_set_from_errno, gretl_errmsg_sprintf;
import gretl.gretl_utils : blas_get_num_threads, blas_is_openblas, blas_is_threaded, blas_set_num_threads, count_distinct_values, doubles_array_free, doubles_array_new, gretl_compare_doubles, gretl_int_from_double, gretl_inverse_compare_doubles;
import gretl.describe : compare_xtab_rows, gretl_array_quantiles;
import gretl.gretl_cmatrix : gretl_cmatrix_print, real_to_complex_fill;
import gretl.random : gretl_rand_int_minmax, gretl_rand_normal, gretl_rand_uniform;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    /* from gretl_cmatrix.h; signature inferred from the call sites (that header was not provided) */
}

/* ---- end of declarations ---- */

/*
 * Libgretl implements most of the matrix functionality that is
 * likely to be required in econometric calculation.  For basics
 * such as decomposition and inversion we use LAPACK as the
 * underlying engine.
 *
 * To get yourself a gretl matrix, use gretl_matrix_alloc() or
 * one of the more specialized constructors; to free such a
 * matrix use gretl_matrix_free().
 */

struct GretlMatrixBlock {
    int n;
    double * val;
    GretlMatrixRaw ** matrix;
}

private bool gretl_is_vector(const(GretlMatrixRaw) * v)
{
    return (v.rows == 1 || v.cols == 1);
}

private bool matrix_is_scalar(const(GretlMatrixRaw) * m)
{
    return (m.rows == 1 && m.cols == 1);
}

private int mdx(const(GretlMatrixRaw) * a, int i, int j)
{
    return (j * a.rows + i);
}

private double matrix_transp_get(const(GretlMatrixRaw) * m, int i, int j)
{
    return m.val[i * m.rows + j];
}

private void matrix_transp_set(GretlMatrixRaw * m, int i, int j, double x)
{
    m.val[i * m.rows + j] = x;
}

private Cmplx cmatrix_transp_get(const(GretlMatrixRaw) * m, int i, int j)
{
    return m.z[i * m.rows + j];
}

private void cmatrix_transp_set(GretlMatrixRaw * m, int i, int j, Cmplx x)
{
    m.z[i * m.rows + j] = x;
}

private enum INFO_INVALID = cast(MatrixInfo *) 0xdeadbeef;

private bool is_block_matrix(const(GretlMatrixRaw) * m)
{
    return (m.info == INFO_INVALID);
}

private bool is_one_by_one(const(GretlMatrixRaw) * m)
{
    return (m.rows == 1 && m.cols == 1);
}

private bool no_metadata(const(GretlMatrixRaw) * m)
{
    return (m.info == null || is_block_matrix(m));
}

private void * mval_malloc(size_t sz)
{
    /* forestall "invalid reads" by OpenBLAS */
    return malloc(sz % 16 ? sz + 8 : sz);
}

private void * mval_realloc(void * ptr, size_t sz)
{
    /* comment as for mval_malloc() */
    return realloc(ptr, sz % 16 ? sz + 8 : sz);
}

private void mval_free(void * m)
{
    free(m);
}

/* Below: setting of the maximal value of K = the shared inner
   dimension in matrix multiplication for use of SIMD. Also
   setting of the minimum value of M x N for doing matrix
   addition and subtraction via SIMD. If these variables are
   set to -1 that disables SIMD by default (though the user
   can change that via the "set" command).

   NB: the SIMD code paths themselves are not part of this port, so
   these settings have no effect on the computations.
*/

private __gshared int simd_k_max = 8;   /* 2014-03-07: was -1 */
private __gshared int simd_mn_min = 16; /* 2014-03-07: was -1 */

void set_simd_k_max(int k)
{
    simd_k_max = k;
}

int get_simd_k_max()
{
    return simd_k_max;
}

void set_simd_mn_min(int mn)
{
    simd_mn_min = mn;
}

int get_simd_mn_min()
{
    return simd_mn_min;
}

private enum SVD_SMIN = 1.0e-9;

/* maybe experiment with these? */
private enum QR_RCOND_MIN  = 1.0e-14;
private enum QR_RCOND_WARN = 1.0e-07;

/* matrix metadata struct, not allocated by default */

struct MatrixInfo {
    int t1;
    int t2;
    char ** colnames;
    char ** rownames;
}

private enum NameFlags {
    colnames = 1 << 0,
    rownames = 1 << 1,
    reversed = 1 << 2
}

/* Central accounting for error in matrix allocation */

private __gshared int gretl_matrix_err;

/* get, and clear, the matrix error code */

int get_gretl_matrix_err()
{
    int ret = gretl_matrix_err;
    gretl_matrix_err = 0;
    return ret;
}

void clear_gretl_matrix_err()
{
    gretl_matrix_err = 0;
}

private void set_gretl_matrix_err(int err)
{
    if (gretl_matrix_err == 0) {
        gretl_matrix_err = err;
    }
}

private int wspace_fail(integer info, double w0)
{
    int iinfo = cast(int) info;

    fprintf(stderr, "gretl_matrix: workspace query failed: info = %d, "
            ~ "work[0] = %g\n", iinfo, w0);
    return GretlError.data;
}

/* An efficient means of allocating temporary storage for lapack
   operations: this should be used _only_ for temporary allocations
   that would ordinarily be freed before returning from the function
   in question.  In this mode we keep the chunk around for future use,
   expanding it as needed.
*/

/* NB: these are thread-local (D's default for module variables), which is
   what the C source's "omp threadprivate" gives in an OpenMP build. */
private void * lapack_mem_chunk;
private size_t lapack_mem_sz;

private void * lapack_malloc(size_t sz)
{
    void * mem = null;

    if (sz > lapack_mem_sz) {
        void * chunk = realloc(lapack_mem_chunk, sz);

        if (chunk != null) {
            lapack_mem_chunk = mem = chunk;
            lapack_mem_sz = sz;
        }
    } else {
        mem = lapack_mem_chunk;
    }

    return mem;
}

private void * lapack_realloc(void * p, size_t sz)
{
    return lapack_malloc(sz);
}

/*
 * lapack_mem_free:
 *
 * Cleanup function, called by libgretl_cleanup(). Frees
 * any memory that has been allocated internally as
 * temporary workspace for LAPACK functions.
 */

void lapack_mem_free()
{
    free(lapack_mem_chunk);
    lapack_mem_chunk = null;
    lapack_mem_sz = 0;
}

private void lapack_free(void * p)
{
    return;
}

private void math_err_init()
{
    errno = 0;
    feclearexcept(FE_ALL_EXCEPT);
}

/* the following is called after math operations only
   if @errno is found to be non-zero */

private int math_err_check(const(char) * msg, int errnum)
{
    int err = GretlError.data;

    if (!fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)) {
        /* we'll let "pure" underflow pass */
        fprintf(stderr, "warning: calculation underflow\n");
        err = 0;
    } else {
        gretl_errmsg_set_from_errno(msg, errnum);
    }
    feclearexcept(FE_ALL_EXCEPT);
    errno = 0;
    return err;
}

private int matrix_block_error(const(char) * f)
{
    fprintf(stderr, "CODING ERROR: illegal call to %s on "
            ~ "member of matrix block\n", f);
    return GretlError.data;
}

/* LIMIT_THREADS (dates from October 2024): this symbol has the effect of
   preventing multi-threading in calls to dpotr* and dgetr* on small- to
   moderate-sized matrices via openblas.  Experimentation has shown that
   using a single thread in this case is a good deal faster, particularly
   on Windows.
*/

/* Minimum matrix sizes for which we'll allow multi-threading, when
   not on MS Windows; as things stand we'll disallow multi-threading
   on Windows unconditionally.
*/
private enum POTRF_MT_MIN = 125;
private enum GETRF_MT_MIN = 80;

private void maybe_force_single(int n, int thresh, int * save_nt)
{
    int nt = blas_get_num_threads();
    int limit = (nt > 1 && n < thresh);

    if (limit) {
        blas_set_num_threads(1);
        *save_nt = nt;
    } else {
        *save_nt = 0;
    }
}

/*
 * gretl_matrix_alloc:
 * @rows: desired number of rows in matrix.
 * @cols: desired number of columns.
 *
 * Returns: pointer to a newly allocated gretl_matrix, or null
 * on failure.  Note that the actual data storage is not
 * initialized.
 */

GretlMatrixRaw * gretl_matrix_alloc(int rows, int cols)
{
    GretlMatrixRaw * m;
    size_t vsize = 0;
    guint64 chk;

    if (rows < 0 || cols < 0) {
        set_gretl_matrix_err(GretlError.invarg);
        fprintf(stderr, "gretl error: gretl_matrix_alloc: rows=%d, cols=%d\n",
                rows, cols);
        return null;
    }

    chk = rows * cast(guint64) cols;
    if (chk > cast(guint64) int.max) {
        set_gretl_matrix_err(GretlError.invarg);
        gretl_errmsg_set(_("Too many elements").ptr);
        return null;
    }

    m = cast(GretlMatrixRaw *) malloc(GretlMatrixRaw.sizeof);
    if (m == null) {
        set_gretl_matrix_err(GretlError.alloc);
        return null;
    }

    vsize = cast(size_t) chk * double.sizeof;
    if (vsize == 0) {
        m.val = null;
    } else {
        m.val = cast(double *) mval_malloc(vsize);
        if (m.val == null) {
            set_gretl_matrix_err(GretlError.alloc);
            free(m);
            return null;
        }
    }

    m.rows = rows;
    m.cols = cols;
    m.is_complex = 0;
    m.z = null;
    m.info = null;

    return m;
}

GretlMatrixRaw * gretl_cmatrix_new(int r, int c)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(2*r, c);

    if (m != null) {
        m.is_complex = 1;
        m.z = cast(Cmplx *) m.val;
        m.rows = r;
    }
    return m;
}

GretlMatrixRaw * gretl_cmatrix_new0(int r, int c)
{
    GretlMatrixRaw * m = gretl_zero_matrix_new(2*r, c);

    if (m != null) {
        m.is_complex = 1;
        m.z = cast(Cmplx *) m.val;
        m.rows = r;
    }
    return m;
}

GretlMatrixRaw * gretl_cmatrix_new1(int r, int c)
{
    GretlMatrixRaw * m = gretl_matrix_alloc(2*r, c);
    int i;

    if (m != null) {
        m.is_complex = 1;
        m.z = cast(Cmplx *) m.val;
        m.rows = r;
        for (i=0; i<r*c; i++) {
            m.z[i] = Cmplx(1.0, 0.0);
        }
    }
    return m;
}

GretlMatrixRaw * gretl_matching_matrix_new(int r, int c,
                                        const(GretlMatrixRaw) * m)
{
    if (m.is_complex) {
        return gretl_cmatrix_new(r, c);
    } else {
        return gretl_matrix_alloc(r, c);
    }
}

private int matrix_set_complex(GretlMatrixRaw * m, int c, int full)
{
    if (m == null) {
        return GretlError.invarg;
    } else if (c && !m.is_complex) {
        if (full && m.rows % 2 != 0) {
            return GretlError.invarg;
        } else {
            m.is_complex = 1;
            m.z = cast(Cmplx *) m.val;
            if (full) {
                m.rows /= 2;
            }
        }
    } else if (!c && m.is_complex) {
        m.is_complex = 0;
        m.z = null;
        if (full) {
            m.rows *= 2;
        }
    } else if (c && m.is_complex) {
        m.z = cast(Cmplx *) m.val;
    }

    return 0;
}

int gretl_matrix_set_complex(GretlMatrixRaw * m, int c)
{
    return matrix_set_complex(m, c, 0);
}

int gretl_matrix_set_complex_full(GretlMatrixRaw * m, int c)
{
    return matrix_set_complex(m, c, 1);
}

void gretl_matrix_block_destroy(GretlMatrixBlock * B)
{
    int i;

    if (B == null) {
        return;
    }

    if (B.matrix != null) {
        for (i=0; i<B.n; i++) {
            free(B.matrix[i]);
        }
        free(B.matrix);
    }

    free(B.val);
    free(B);
}

void gretl_matrix_block_zero(GretlMatrixBlock * B)
{
    if (B != null && B.matrix != null) {
        int i;

        for (i=0; i<B.n; i++) {
            gretl_matrix_zero(B.matrix[i]);
        }
    }
}

/* Create an array of n matrices.  The "..." should be filled with (at
   minimum) the number of rows and columns for the first matrix to
   create, which will be written to the location given by @pm.
   Following this there can be any number of triples of type
   (gretl_matrix **, int, int), representing the location of a matrix
   to be filled out and the desired number of rows and columns,
   respectively.  The argument list must be terminated by null.

   This is supposed to economize on calls to alloc() and free(), since
   we allocate a combined data block for all the matrices in the
   array.

   Matrices in this array should be destroyed by calling
   gretl_matrix_block_destroy() on the block -- do NOT call
   gretl_matrix_free on individual member-matrices.

   NB: this is a C-style variadic (extern (C)); callers must pass
   exactly the argument types described above (int, int, GretlMatrix **).
*/

extern (C) GretlMatrixBlock * gretl_matrix_block_new(GretlMatrixRaw ** pm, ...)
{
    va_list ap;
    GretlMatrixBlock * B;
    GretlMatrixRaw ** targ;
    GretlMatrixRaw * m;
    size_t vsize = 0;
    int i, err = 0;

    B = cast(GretlMatrixBlock *) malloc(GretlMatrixBlock.sizeof);
    if (B == null) {
        return null;
    }

    /* first pass: determine the number of
       (pointer, int, int) triples */
    va_start(ap, pm);
    for (i=1; ; i++) {
        va_arg!int(ap);
        va_arg!int(ap);
        targ = va_arg!(GretlMatrixRaw **)(ap);
        if (targ == null) {
            break;
        }
    }
    va_end(ap);

    /* initialize B */
    B.n = i;
    B.matrix = cast(GretlMatrixRaw **) malloc(B.n * (*B.matrix).sizeof);
    if (B.matrix == null) {
        free(B);
        return null;
    }

    /* null everything in case we fail */
    B.val = null;
    for (i=0; i<B.n; i++) {
        B.matrix[i] = null;
    }

    /* now allocate and initialize the matrices */
    for (i=0; i<B.n; i++) {
        B.matrix[i] = cast(GretlMatrixRaw *) malloc(GretlMatrixRaw.sizeof);
        if (B.matrix[i] == null) {
            gretl_matrix_block_destroy(B);
            return null;
        }
        B.matrix[i].info = INFO_INVALID;
        B.matrix[i].val = null;
        B.matrix[i].z = null;
        B.matrix[i].is_complex = 0;
    }

    /* second pass through arg list */

    va_start(ap, pm);

    for (i=0; i<B.n; i++) {
        m = B.matrix[i];
        if (i == 0) {
            *pm = m;
        } else {
            targ = va_arg!(GretlMatrixRaw **)(ap);
            *targ = m;
        }
        m.rows = va_arg!int(ap);
        m.cols = va_arg!int(ap);
        if (m.rows < 0 || m.cols < 0) {
            err = 1;
            break;
        }
        vsize += m.rows * m.cols;
    }

    va_end(ap);

    if (!err && vsize > 0) {
        /* allocate combined data block */
        B.val = cast(double *) malloc(vsize * (*B.val).sizeof);
        if (B.val == null) {
            err = 1;
        }
    }

    if (err) {
        gretl_matrix_block_destroy(B);
        B = null;
    } else {
        /* set the val pointers */
        double * val = B.val;
        int n;

        for (i=0; i<B.n; i++) {
            m = B.matrix[i];
            n = m.rows * m.cols;
            if (n > 0) {
                m.val = val;
                val += n;
            }
        }
    }

    return B;
}

int gretl_matrix_block_n_matrices(GretlMatrixBlock * B)
{
    return B == null ? 0 : B.n;
}

GretlMatrixRaw * gretl_matrix_block_get_matrix(GretlMatrixBlock * B,
                                            int i)
{
    if (B == null || i < 0 || i >= B.n) {
        return null;
    } else {
        return B.matrix[i];
    }
}

int gretl_matrix_is_nonfinite(const(GretlMatrixRaw) * m)
{
    if (m != null) {
        int i, n = m.rows * m.cols;

        for (i=0; i<n; i++) {
            if (!isFinite(m.val[i])) {
                return GretlError.nonfin;
            }
        }
    }

    return 0;
}

int gretl_matrix_get_structure(const(GretlMatrixRaw) * m)
{
    int ret = 0;

    if (gretl_is_null_matrix(m)) {
        return 0;
    }

    if (m.rows == m.cols) {
        ret = GretlMatrixStructure.square;
        if (m.rows == 1) {
            ret = GretlMatrixStructure.scalar;
        }
    }

    /* FIXME the tests in the next block of code are not applicable
       for complex matrices.
    */

    if (ret == GretlMatrixStructure.square) {
        guint8 uzero = 1;
        guint8 lzero = 1;
        guint8 symm = 1;
        guint8 udiag = 1;
        guint8 hermt = 0;
        int i, j;
        int k = 0;

        if (m.is_complex) {
            Cmplx z_0 = Cmplx(0.0, 0.0);
            Cmplx z_1 = Cmplx(1.0, 0.0);
            Cmplx zij;

            hermt = 1;
            for (j=0; j<m.cols; j++) {
                for (i=0; i<m.rows; i++) {
                    zij = m.z[k++];
                    if (j > i) {
                        if (uzero && zij != z_0) {
                            uzero = 0;
                        }
                    } else if (i > j) {
                        if (lzero && zij != z_0) {
                            lzero = 0;
                        }
                    } else if (i == j) {
                        if (udiag && zij != z_1) {
                            udiag = 0;
                        }
                    }
                    if (j != i && symm) {
                        if (zij != gretl_cmatrix_get(m,j,i)) {
                            symm = 0;
                        }
                    }
                    if (j != i && hermt) {
                        if (zij != conj(gretl_cmatrix_get(m,j,i))) {
                            hermt = 0;
                        }
                    }
                }
                if (!uzero && !lzero && !symm && !hermt) {
                    break;
                }
            }
        } else {
            double x;

            for (j=0; j<m.cols; j++) {
                for (i=0; i<m.rows; i++) {
                    x = m.val[k++];
                    if (j > i) {
                        if (uzero && x != 0.0) {
                            uzero = 0;
                        }
                    } else if (i > j) {
                        if (lzero && x != 0.0) {
                            lzero = 0;
                        }
                    } else if (i == j) {
                        if (udiag && x != 1.0) {
                            udiag = 0;
                        }
                    }
                    if (j != i && symm) {
                        if (x != gretl_matrix_get(m,j,i)) {
                            symm = 0;
                        }
                    }
                    if (!uzero && !lzero && !symm) {
                        break;
                    }
                }
                if (!uzero && !lzero && !symm) {
                    break;
                }
            }
        }

        if (udiag && uzero && lzero) {
            ret = GretlMatrixStructure.identity;
        } else if (uzero && lzero) {
            ret = GretlMatrixStructure.diagonal;
        } else if (uzero) {
            ret = GretlMatrixStructure.lower_triangular;
        } else if (lzero) {
            ret = GretlMatrixStructure.upper_triangular;
        } else if (symm) {
            ret = GretlMatrixStructure.symmetric;
        } else if (hermt) {
            ret = GretlMatrixStructure.hermitian;
        }
    }

    return ret;
}

private int matrix_is_triangular(const(GretlMatrixRaw) * m)
{
    double x;
    guint8 uzero = 1;
    guint8 lzero = 1;
    int i, j;
    int k = 0;

    if (m.cols != m.rows) {
        /* m must be square */
        return 0;
    }

    for (j=0; j<m.cols; j++) {
        for (i=0; i<m.rows; i++) {
            x = m.val[k++];
            if (j > i) {
                if (uzero && x != 0.0) {
                    uzero = 0;
                }
            } else if (i > j) {
                if (lzero && x != 0.0) {
                    lzero = 0;
                }
            }
            if (!uzero && !lzero) {
                break;
            }
        }
        if (!uzero && !lzero) {
            break;
        }
    }

    if (uzero) {
        return GretlMatrixStructure.lower_triangular;
    } else if (lzero) {
        return GretlMatrixStructure.upper_triangular;
    }

    return 0;
}

/*
 * gretl_matrix_reuse:
 * @m: matrix to reuse.
 * @rows: desired number of rows in "new" matrix, or -1
 * to leave the current value unchanged.
 * @cols: desired number of columns in "new" matrix, or -1
 * to leave the current value unchanged.
 *
 * An "experts only" memory-conservation trick. If @m is an
 * already-allocated gretl matrix, you can "resize" it by
 * specifying a new number of rows and columns.  This works
 * only if the product of @rows and @cols is less than or equal
 * to the product of the number of rows and columns in the
 * matrix as originally allocated; no actual reallocation of memory
 * is performed.  If you "reuse" with an excessive number of rows
 * or columns you will surely crash your program or smash the
 * stack. Note also that the matrix-pointer returned is not really
 * new, and when the matrix is to be freed, gretl_matrix_free()
 * should be applied only once.
 *
 * Returns: pointer to the "resized" gretl_matrix.
 */

GretlMatrixRaw * gretl_matrix_reuse(GretlMatrixRaw * m, int rows, int cols)
{
    int r = rows > 0 ? rows : m.rows;
    int c = cols > 0 ? cols : m.cols;

    m.rows = r;
    m.cols = c;

    return m;
}

/*
 * gretl_matrix_realloc:
 * @m: matrix to reallocate.
 * @rows: desired number of rows in "new" matrix.
 * @cols: desired number of columns in "new" matrix.
 *
 * Reallocates the storage in @m to the specified dimensions.
 *
 * Returns: 0 on success, GretlError.alloc on failure.
 */

int gretl_matrix_realloc(GretlMatrixRaw * m, int rows, int cols)
{
    int oldrows, oldcols;
    size_t n;
    guint64 chk;
    double * x = null;

    if (m == null) {
        return GretlError.data;
    }

    if (rows < 0 || cols < 0) {
        fprintf(stderr, "gretl error: gretl_matrix_realloc: rows=%d, cols=%d\n",
                rows, cols);
        return GretlError.invarg;
    }

    chk = rows * cast(guint64) cols;
    if (chk > cast(guint64) int.max) {
        gretl_errmsg_set(_("Too many elements").ptr);
        return GretlError.invarg;
    }

    if (rows == m.rows && cols == m.cols) {
        /* no-op */
        return 0;
    }

    n = rows * cols;

    if (n == m.rows * m.cols) {
        /* no need to reallocate storage */
        m.rows = rows;
        m.cols = cols;
        gretl_matrix_destroy_info(m);
        return 0;
    }

    if (is_block_matrix(m)) {
        matrix_block_error("gretl_matrix_realloc");
        return GretlError.data;
    }

    if (n == 0) {
        mval_free(m.val);
        m.val = null;
    } else {
        size_t sz = cast(size_t) chk * double.sizeof;

        if (m.is_complex) {
            x = cast(double *) mval_realloc(m.val, 2 * sz);
        } else {
            x = cast(double *) mval_realloc(m.val, sz);
        }
        if (x == null) {
            return GretlError.alloc;
        }
    }

    oldrows = m.rows;
    oldcols = m.cols;

    m.val = x;
    m.rows = rows;
    m.cols = cols;
    if (m.is_complex) {
        m.z = cast(Cmplx *) m.val;
    }

    if (m.info != null) {
        if (m.rows != oldrows && m.cols != oldcols) {
            gretl_matrix_destroy_info(m);
        } else if (m.rows != oldrows && m.info.rownames != null) {
            strings_array_free(m.info.rownames, oldrows);
            m.info.rownames = null;
        } else if (m.cols != oldcols && m.info.colnames != null) {
            strings_array_free(m.info.colnames, oldcols);
            m.info.colnames = null;
        }
    }

    return 0;
}

/*
 * gretl_matrix_init_full:
 * @m: matrix to be initialized.
 * @rows: number of rows.
 * @cols: number of columns.
 * @val: data array.
 *
 * Initializes @m to be @rows by @cols and have data @val. This
 * intended for use with automatic matrices declared "on
 * the stack". It is up to the user to ensure that the size of
 * @val is compatible with the @rows and @cols specification.
 */

GretlMatrixRaw * gretl_matrix_init_full(GretlMatrixRaw * m,
                                     int rows, int cols,
                                     double * val)
{
    m.rows = rows;
    m.cols = cols;
    m.val = val;
    m.info = null;
    m.is_complex = 0;
    m.z = null;
    return m;
}

/*
 * gretl_matrix_init:
 * @m: matrix to be initialized.
 *
 * Initializes @m to be zero by zero with null data.
 */

GretlMatrixRaw * gretl_matrix_init(GretlMatrixRaw * m)
{
    m.rows = m.cols = 0;
    m.val = null;
    m.info = null;
    m.is_complex = 0;
    m.z = null;
    return m;
}

/*
 * gretl_matrix_replace:
 * @pa: location of matrix to be replaced.
 * @b: replacement matrix.
 *
 * Frees the matrix at location @pa and substitutes @b.
 *
 * Returns: the replacement matrix.
 */

GretlMatrixRaw * gretl_matrix_replace(GretlMatrixRaw ** pa,
                                   GretlMatrixRaw * b)
{
    gretl_matrix_free(*pa);
    *pa = b;
    return b;
}

/*
 * gretl_matrix_replace_content:
 * @targ: matrix to receive new content.
 * @donor: matrix to donate content.
 *
 * Moves the content of @donor into @targ; @donor becomes
 * a null matrix in consequence.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_replace_content(GretlMatrixRaw * targ,
                                 GretlMatrixRaw * donor)
{
    if (is_block_matrix(targ) || is_block_matrix(donor)) {
        matrix_block_error("gretl_matrix_replace_content");
        return GretlError.data;
    } else {
        gretl_matrix_destroy_info(targ);
        free(targ.val);
        targ.rows = donor.rows;
        targ.cols = donor.cols;
        targ.val = donor.val;
        donor.val = null;
        gretl_matrix_set_complex(targ, donor.is_complex);
        return 0;
    }
}

/*
 * gretl_identity_matrix_new:
 * @n: desired number of rows and columns in the matrix.
 *
 * Returns: pointer to a newly allocated identity matrix, or null
 * on failure.
 */

GretlMatrixRaw * gretl_identity_matrix_new(int n)
{
    GretlMatrixRaw * m;
    int i, k;

    if (n < 0) {
        return null;
    } else if (n == 0) {
        return gretl_null_matrix_new();
    }

    m = gretl_matrix_alloc(n, n);

    if (m != null) {
        k = n * n;
        n++;
        for (i=0; i<k; i++) {
            m.val[i] = (i % n) ? 0.0 : 1.0;
        }
    }

    return m;
}

/*
 * gretl_DW_matrix_new:
 * @n: desired number of rows and columns in the matrix.
 *
 * Returns: pointer to a newly allocated Durbin-Watson matrix, or null
 * on failure.  This is a tridiagonal matrix with 2 on the leading
 * diagonal (apart from 1 at the ends) and -1 on the supra- and
 * infra-diagonals.
 */

GretlMatrixRaw * gretl_DW_matrix_new(int n)
{
    GretlMatrixRaw * m = gretl_zero_matrix_new(n, n);
    int i, j;

    if (m == null) {
        return null;
    }

    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if (j == i) {
                if (i == 0 || i == n-1) {
                    gretl_matrix_set(m, i, j, 1.0);
                } else {
                    gretl_matrix_set(m, i, j, 2.0);
                }
            } else if (j == i + 1 || i == j + 1) {
                gretl_matrix_set(m, i, j, -1.0);
            }
        }
    }

    return m;
}

private GretlMatrixRaw * gretl_filled_matrix_new(int r, int c,
                                              double val)
{
    GretlMatrixRaw * m = null;

    if (r < 0 || c < 0) {
        return null;
    } else if (r == 0 || c == 0) {
        m = gretl_null_matrix_new();
        if (m != null) {
            m.rows = r;
            m.cols = c;
        }
    } else {
        int i, n = r * c;

        m = gretl_matrix_alloc(r, c);
        if (m != null) {
            if (val == 0.0) {
                memset(m.val, 0, n * (*m.val).sizeof);
            } else {
                for (i=0; i<n; i++) {
                    m.val[i] = val;
                }
            }
        }
    }

    return m;
}

/*
 * gretl_zero_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated zero matrix, or null
 * on failure.
 */

GretlMatrixRaw * gretl_zero_matrix_new(int r, int c)
{
    return gretl_filled_matrix_new(r, c, 0.0);
}

/*
 * gretl_unit_matrix_new:
 * @r: desired number of rows in the matrix.
 * @c: desired number of columns in the matrix.
 *
 * Returns: pointer to a newly allocated matrix, all
 * of whose elements equal 1, or null on failure.
 */

GretlMatrixRaw * gretl_unit_matrix_new(int r, int c)
{
    return gretl_filled_matrix_new(r, c, 1.0);
}

/*
 * gretl_null_matrix_new:
 *
 * Returns: pointer to a newly allocated null matrix, or
 * null on failure.
 */

GretlMatrixRaw * gretl_null_matrix_new()
{
    GretlMatrixRaw * m = cast(GretlMatrixRaw *) malloc(GretlMatrixRaw.sizeof);

    if (m == null) {
        set_gretl_matrix_err(GretlError.alloc);
        return null;
    }

    gretl_matrix_init(m);

    return m;
}

/*
 * gretl_matrix_seq:
 * @start: first element.
 * @end: last element.
 * @step: positive step.
 * @err: location to recieve error code.
 *
 * Returns: pointer to a row vector, containing values from
 * @start to @end, in decreasing order if @start > @end --
 * or null on failure.
 */

GretlMatrixRaw * gretl_matrix_seq(double start, double end,
                               double step, int * err)
{
    GretlMatrixRaw * v;
    int reverse = (start > end);
    double k = start;
    int i, n = 0;

    if (step <= 0) {
        *err = GretlError.data;
        return null;
    }

    if (step == 1.0) {
        if (reverse) {
            n = cast(int) (start - end + 1);
            step = -step;
        } else {
            n = cast(int) (end - start + 1);
        }
    } else if (reverse) {
        step = -step;
        while (k >= end) {
            n++;
            k += step;
        }
    } else {
        while (k <= end) {
            n++;
            k += step;
        }
    }

    if (n == 0) {
        *err = GretlError.data;
        return null;
    }
    v = gretl_vector_alloc(n);

    if (v == null) {
        *err = GretlError.alloc;
    } else {
        k = start;
        if (step == 1.0) {
            for (i=0; i<n; i++) {
                v.val[i] = k++;
            }
        } else if (step == -1.0) {
            for (i=0; i<n; i++) {
                v.val[i] = k--;
            }
        } else {
            for (i=0; i<n; i++) {
                v.val[i] = k;
                k += step;
            }
        }
    }

    return v;
}

/*
 * gretl_matrix_fill:
 * @m: matrix to fill.
 * @x: value with which to fill.
 *
 * Sets all entries in @m to the value @x.
 */

void gretl_matrix_fill(GretlMatrixRaw * m, double x)
{
    if (m != null) {
        int i, n = m.rows * m.cols;

        if (m.is_complex) {
            for (i=0; i<n; i++) {
                m.z[i] = Cmplx(x, 0.0);
            }
        } else {
            for (i=0; i<n; i++) {
                m.val[i] = x;
            }
        }
    }
}

private GretlMatrixRaw *
gretl_matrix_copy_mod(const(GretlMatrixRaw) * m, int mod)
{
    GretlMatrixRaw * c;
    int rows, cols;
    int i, j;

    if (m == null) {
        return null;
    }

    if (mod == GretlMatrixMod.transpose) {
        rows = m.cols;
        cols = m.rows;
    } else {
        rows = m.rows;
        cols = m.cols;
    }

    c = gretl_matching_matrix_new(rows, cols, m);
    if (c == null) {
        return null;
    }

    if (mod == GretlMatrixMod.transpose) {
        int k = 0;

        if (m.is_complex) {
            /* we'll do the conjugate transpose */
            Cmplx mij;

            for (j=0; j<m.cols; j++) {
                for (i=0; i<m.rows; i++) {
                    mij = m.z[k++];
                    gretl_cmatrix_set(c, j, i, conj(mij));
                }
            }
        } else {
            /* a real matrix */
            if (min(m.cols, m.rows) > 1) {
                for (j=0; j<m.cols; j++) {
                    for (i=0; i<m.rows; i++) {
                        gretl_matrix_set(c, j, i, m.val[k++]);
                    }
                }
            } else {
                memcpy(c.val, m.val, rows * cols * (*m.val).sizeof);
            }
        }
    } else {
        /* not transposing */
        int n = rows * cols;

        if (m.is_complex) {
            memcpy(c.z, m.z, n * (*m.z).sizeof);
        } else {
            memcpy(c.val, m.val, n * (*m.val).sizeof);
        }
    }

    gretl_matrix_copy_info(c, m, mod);

    return c;
}

private GretlMatrixRaw * matrix_copy_plain(const(GretlMatrixRaw) * m)
{
    GretlMatrixRaw * c;

    if (m == null) {
        return null;
    }

    c = gretl_matching_matrix_new(m.rows, m.cols, m);

    if (c != null) {
        int n = c.rows * c.cols;

        if (m.is_complex) n *= 2;
        memcpy(c.val, m.val, n * (*m.val).sizeof);
    }

    return c;
}

/*
 * gretl_matrix_copy:
 * @m: source matrix to be copied.
 *
 * Returns: an allocated copy of matrix @m, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_copy(const(GretlMatrixRaw) * m)
{
    return gretl_matrix_copy_mod(m, GretlMatrixMod.none);
}

/*
 * gretl_matrix_copy_transpose:
 * @m: source matrix to be copied.
 *
 * Returns: an allocated copy of the tranpose of @m, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_copy_transpose(const(GretlMatrixRaw) * m)
{
    return gretl_matrix_copy_mod(m, GretlMatrixMod.transpose);
}

/* Relatively lightweight version of gretl_matrix_copy, for internal
   use when we just want a temporary copy of an original matrix as
   workspace, and we know that the original is not a null matrix.
*/

private GretlMatrixRaw * gretl_matrix_copy_tmp(const(GretlMatrixRaw) * a)
{
    size_t sz = a.rows * a.cols * double.sizeof;
    GretlMatrixRaw * b = cast(GretlMatrixRaw *) calloc(1, GretlMatrixRaw.sizeof);

    if (a.is_complex) sz *= 2;

    if (b != null && (b.val = cast(double *) mval_malloc(sz)) != null) {
        b.rows = a.rows;
        b.cols = a.cols;
        b.info = null;
        memcpy(b.val, a.val, sz);
        gretl_matrix_set_complex(b, a.is_complex);
    }

    return b;
}

/*
 * gretl_matrix_copy_row:
 * @dest: destination matrix.
 * @di: row to copy into.
 * @src: source matrix.
 * @si: row to copy from.
 *
 * Copies the values from row @si of @src into row @di
 * of @dest, provided @src and @dest have the same number
 * of columns.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_matrix_copy_row(GretlMatrixRaw * dest, int di,
                          const(GretlMatrixRaw) * src, int si)
{
    int err = 0;

    if (dest == null || src == null ||
        gretl_is_null_matrix(dest) ||
        gretl_is_null_matrix(src)) {
        err = GretlError.data;
    } else if (dest.cols != src.cols) {
        err = GretlError.nonconf;
    } else {
        double x;
        int j;

        for (j=0; j<src.cols; j++) {
            x = gretl_matrix_get(src, si, j);
            gretl_matrix_set(dest, di, j, x);
        }
    }

    return err;
}

/* Mechanism for copying row or column names from matrix
   @src to matrix @targ. If @sel is non-null it must be a
   boolean vector indicating that @targ holds a subset of
   the columns or rows of @src.
*/

private void maybe_preserve_names(GretlMatrixRaw * targ,
                                  const(GretlMatrixRaw) * src,
                                  NameFlags flags,
                                  const(GretlMatrixRaw) * sel)
{
    int cols = (flags & NameFlags.colnames);
    int reverse = (flags & NameFlags.reversed);
    char ** srcnames, S;
    int ns, nt, err = 0;

    if (no_metadata(src)) {
        return;
    } else if (cols) {
        srcnames = cast(char **) src.info.colnames;
        ns = src.cols;
        nt = targ.cols;
    } else {
        srcnames = cast(char **) src.info.rownames;
        ns = src.rows;
        nt = targ.rows;
    }

    if (srcnames == null || nt > ns) {
        return;
    } else if (nt < ns && sel == null) {
        return;
    }

    if (sel != null) {
        int i, n = gretl_vector_get_length(sel);
        int k = 0;

        for (i=0; i<n; i++) {
            k += (sel.val[i] != 0);
        }
        S = strings_array_new(k);
        if (S != null) {
            k = 0;
            for (i=0; i<n; i++) {
                if (sel.val[i] != 0) {
                    S[k++] = gretl_strdup(srcnames[i]);
                }
            }
        }
    } else if (reverse) {
        S = strings_array_reverse(srcnames, ns);
    } else {
        S = strings_array_dup(srcnames, ns);
    }

    if (S != null) {
        if (cols) {
            err = gretl_matrix_set_colnames(targ, S);
        } else {
            err = gretl_matrix_set_rownames(targ, S);
        }
        if (err) {
            strings_array_free(S, nt);
        }
    }
}

private void maybe_concat_names(GretlMatrixRaw * targ,
                                const(GretlMatrixRaw) * src1,
                                const(GretlMatrixRaw) * src2,
                                NameFlags flags)
{
    int cols = (flags & NameFlags.colnames);
    char ** srcnames1 = null;
    char ** srcnames2 = null;
    char ** S = null;
    int n1, ns, err = 0;

    if (no_metadata(src1) || no_metadata(src2)) {
        return;
    } else if (cols) {
        if (targ.cols == src1.cols + src2.cols) {
            srcnames1 = cast(char **) src1.info.colnames;
            srcnames2 = cast(char **) src2.info.colnames;
        }
    } else {
        if (targ.rows == src1.rows + src2.rows) {
            srcnames1 = cast(char **) src1.info.rownames;
            srcnames2 = cast(char **) src2.info.rownames;
        }
    }

    if (srcnames1 == null || srcnames2 == null) {
        return;
    }

    n1 = cols ? src1.cols : src1.rows;
    ns = cols ? targ.cols : targ.rows;
    S = strings_array_new(ns);

    if (S != null) {
        int i, j = 0, k = 0;

        for (i=0; i<ns; i++) {
            if (i < n1) {
                S[i] = gretl_strdup(srcnames1[j++]);
            } else {
                S[i] = gretl_strdup(srcnames2[k++]);
            }
        }
        if (cols) {
            err = gretl_matrix_set_colnames(targ, S);
        } else {
            err = gretl_matrix_set_rownames(targ, S);
        }
        if (err) {
            strings_array_free(S, ns);
        }
    }
}

/*
 * gretl_matrix_reverse_rows:
 * @m: source matrix whose rows are to be reversed.
 * @err: location to receive error code.
 *
 * Returns: a matrix with the same rows as @m, last to first.
 */

GretlMatrixRaw * gretl_matrix_reverse_rows(const(GretlMatrixRaw) * m,
                                        int * err)
{
    GretlMatrixRaw * ret;
    int i, j, r, c;

    if (m == null) {
        *err = GretlError.invarg;
        return null;
    } else if (gretl_is_null_matrix(m)) {
        return gretl_null_matrix_new();
    }

    r = m.rows;
    c = m.cols;

    ret = gretl_matching_matrix_new(r, c, m);

    if (ret == null) {
        *err = GretlError.alloc;
    } else {
        Cmplx z;
        double x;

        for (i=0; i<r; i++) {
            for (j=0; j<m.cols; j++) {
                if (m.is_complex) {
                    z = gretl_cmatrix_get(m, r-i-1, j);
                    gretl_cmatrix_set(ret, i, j, z);
                } else {
                    x = gretl_matrix_get(m, r-i-1, j);
                    gretl_matrix_set(ret, i, j, x);
                }
            }
        }
        maybe_preserve_names(ret, m, NameFlags.rownames | NameFlags.reversed, null);
        maybe_preserve_names(ret, m, NameFlags.colnames, null);
    }

    return ret;
}

/*
 * gretl_matrix_reverse_cols:
 * @m: source matrix whose columns are to be reversed.
 * @err: location to receive error code.
 *
 * Returns: a matrix with the same columns as @m, last to first.
 */

GretlMatrixRaw * gretl_matrix_reverse_cols(const(GretlMatrixRaw) * m,
                                        int * err)
{
    GretlMatrixRaw * ret;
    const(double) * x;
    double * y;
    size_t csize;
    int i, r, c;

    if (m == null) {
        *err = GretlError.invarg;
        return null;
    } else if (gretl_is_null_matrix(m)) {
        return gretl_null_matrix_new();
    }

    r = m.rows;
    c = m.cols;
    ret = gretl_matching_matrix_new(r, c, m);

    if (ret == null) {
        *err = GretlError.alloc;
    } else {
        if (m.is_complex) {
            r *= 2;
        }
        x = m.val;
        y = ret.val + r * (c-1);
        csize = r * (*x).sizeof;

        for (i=0; i<c; i++) {
            memcpy(y, x, csize);
            x += r;
            y -= r;
        }

        maybe_preserve_names(ret, m, NameFlags.colnames | NameFlags.reversed, null);
        maybe_preserve_names(ret, m, NameFlags.rownames, null);
    }

    return ret;
}

void gretl_matrix_destroy_info(GretlMatrixRaw * m)
{
    if (m != null && m.info != null && !is_block_matrix(m)) {
        strings_array_free(m.info.colnames, m.cols);
        strings_array_free(m.info.rownames, m.rows);
        free(m.info);
        m.info = null;
    }
}

/*
 * gretl_matrix_free:
 * @m: matrix to be freed.
 *
 * Frees the allocated storage in @m, then @m itself.
 */

void gretl_matrix_free(GretlMatrixRaw * m)
{
    if (m == null) return;

    if (is_block_matrix(m)) {
        matrix_block_error("gretl_matrix_free");
        return;
    }

    if (m.val != null) {
        mval_free(m.val);
    }

    if (m.info != null) {
        gretl_matrix_destroy_info(m);
    }

    free(m);
}

/*
 * gretl_matrix_zero:
 * @m: matrix to be set to zero.
 *
 * Sets all elements of @m to zero.
 */

void gretl_matrix_zero(GretlMatrixRaw * m)
{
    int i, n = m.rows * m.cols;

    for (i=0; i<n; i++) {
        m.val[i] = 0.0;
    }
}

/*
 * gretl_matrix_get_diagonal:
 * @m: input matrix.
 * @err: location to receive error code.
 *
 * Returns: a column vector containing the diagonal elements of
 * @m, otherwise null.  A non-zero value is assigned via @err
 * on failure.
 */

GretlMatrixRaw * gretl_matrix_get_diagonal(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * d = null;
    int i, n = 0;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
        d = gretl_null_matrix_new();
    } else {
        n = min(m.rows, m.cols);
        d = gretl_matching_matrix_new(n, 1, m);
    }

    if (d == null) {
        *err = GretlError.alloc;
    } else {
        if (m.is_complex) {
            for (i=0; i<n; i++) {
                d.z[i] = gretl_cmatrix_get(m, i, i);
            }
        } else {
            for (i=0; i<n; i++) {
                d.val[i] = gretl_matrix_get(m, i, i);
            }
        }
    }

    return d;
}

/*
 * gretl_matrix_set_diagonal:
 * @targ: target matrix.
 * @src: source vector (or null).
 * @x: (alternative) source scalar.
 *
 * Sets the diagonal elements of @targ using the elements of
 * @src, if non-null, or otherwise the constant value @x.
 * If @src is given it must be a vector of length equal to
 * that of the diagonal of @targ (that is, the minimum of
 * its rows and columns).
 *
 * Returns: 0 on success, error code on non-conformability.
 */

int gretl_matrix_set_diagonal(GretlMatrixRaw * targ,
                              const(GretlMatrixRaw) * src,
                              double x)
{
    int i, n, match = 0;
    int err = 0;

    if (gretl_is_null_matrix(targ) || targ.is_complex) {
        return GretlError.invarg;
    } else if (src != null && src.is_complex) {
        return GretlError.invarg;
    }

    n = min(targ.rows, targ.cols);

    if (src != null) {
        if (gretl_vector_get_length(src) == n) {
            match = 1;
        } else if (gretl_matrix_is_scalar(src)) {
            x = src.val[0];
            match = 2;
        }
    } else {
        match = 2;
    }

    if (match == 0) {
        err = GretlError.nonconf;
    } else {
        for (i=0; i<n; i++) {
            if (match == 1) {
                gretl_matrix_set(targ, i, i, src.val[i]);
            } else {
                gretl_matrix_set(targ, i, i, x);
            }
        }
    }

    return err;
}

/*
 * gretl_matrix_set_triangle:
 * @targ: target matrix.
 * @src: source vector (or null).
 * @x: (alternative) source scalar.
 * @upper: flag to set the upper part, the default
 * being to set the lower.
 *
 * Sets the lower or upper elements of the matrix
 * @targ using the elements of @src, if non-null, or
 * otherwise the constant value @x.
 *
 * If @src is given it must be a vector of length equal to
 * that of the number of infra- or supra-diagonal elements
 * of @targ.
 *
 * Returns: 0 on success, error code on non-conformability.
 */

int gretl_matrix_set_triangle(GretlMatrixRaw * targ,
                              const(GretlMatrixRaw) * src,
                              double x, int upper)
{
    int r, c, p, i, j, n;
    int lower = !upper;
    int match = 0;
    int err = 0;

    if (gretl_is_null_matrix(targ) || targ.is_complex) {
        return GretlError.invarg;
    } else if (src != null && src.is_complex) {
        return GretlError.invarg;
    }

    r = targ.rows;
    c = targ.cols;

    if ((c == 1 && upper) || (r == 1 && !upper)) {
        /* no such part */
        return GretlError.invarg;
    }

    p = min(r, c);
    n = (p * (p-1)) / 2;

    if (r > c && lower) {
        n += (r - c) * c;
    } else if (c > r && upper) {
        n += (c - r) * r;
    }

    if (src != null) {
        if (gretl_vector_get_length(src) == n) {
            match = 1;
        } else if (gretl_matrix_is_scalar(src)) {
            x = src.val[0];
            match = 2;
        }
    } else {
        match = 2;
    }

    if (match == 0) {
        err = GretlError.nonconf;
    } else {
        int jmin = upper ? 1 : 0;
        int jmax = upper ? c : r;
        int imin = upper ? 0 : 1;
        int imax = upper ? 1 : r;
        int k = 0;

        for (j=jmin; j<jmax; j++) {
            for (i=imin; i<imax; i++) {
                if (src != null) {
                    x = src.val[k++];
                }
                gretl_matrix_set(targ, i, j, x);
            }
            if (lower) {
                imin++;
            } else if (imax < r) {
                imax++;
            }
        }
    }

    return err;
}

/*
 * gretl_matrix_get_triangle:
 * @m: source matrix (real or complex).
 * @upper: flag to get the upper part, the default
 * being to get the lower part.
 * @err: location to receive error code.
 *
 * Returns: A column vector holding the vec of either the
 * infra- or supra-diagonal elements of @m, or null on failure.
 * Note that the "part" returned may not be an actual
 * triangle if @m is not square.
 */

GretlMatrixRaw * gretl_matrix_get_triangle(const(GretlMatrixRaw) * m,
                                        int upper, int * err)
{
    GretlMatrixRaw * ret = null;
    int r, c, p, n, i, j;
    int lower = !upper;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.invarg;
        return null;
    }

    r = m.rows;
    c = m.cols;

    if ((c == 1 && upper) || (r == 1 && !upper)) {
        /* no such part is available */
        *err = GretlError.invarg;
        return null;
    }

    p = min(r, c);
    n = (p * (p-1)) / 2;

    if (r > c && lower) {
        n += (r - c) * c;
    } else if (c > r && upper) {
        n += (c - r) * r;
    }

    ret = gretl_matching_matrix_new(n, 1, m);

    if (ret == null) {
        *err = GretlError.alloc;
    } else {
        int jmin = upper ? 1 : 0;
        int jmax = upper ? c : r;
        int imin = upper ? 0 : 1;
        int imax = upper ? 1 : r;
        int k = 0;

        for (j=jmin; j<jmax; j++) {
            for (i=imin; i<imax; i++) {
                if (m.is_complex) {
                    ret.z[k++] = gretl_cmatrix_get(m, i, j);
                } else {
                    ret.val[k++] = gretl_matrix_get(m, i, j);
                }
            }
            if (lower) {
                imin++;
            } else if (imax < r) {
                imax++;
            }
        }
    }

    return ret;
}

/*
 * gretl_matrix_get_row:
 * @m: input matrix.
 * @i: index of row to access.
 * @v: location to receive row values.
 *
 * Copies row @i of matrix @m into vector @v.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_matrix_get_row(const(GretlMatrixRaw) * m, int i, GretlVector * v)
{
    int j, nc = gretl_matrix_cols(m);

    if (gretl_vector_get_length(v) != nc) {
        return GretlError.nonconf;
    }

    for (j=0; j<nc; j++) {
        gretl_vector_set(v, j, gretl_matrix_get(m, i, j));
    }

    return 0;
}

/*
 * gretl_matrix_trace:
 * @m: square input matrix.
 *
 * Returns: the trace (sum of diagonal elements) of @m, if
 * @m is square, otherwise NADBL.
 */

double gretl_matrix_trace(const(GretlMatrixRaw) * m)
{
    double tr = 0.0;
    int i;

    if (gretl_is_null_matrix(m) || m.rows != m.cols) {
        return NADBL;
    }

    for (i=0; i<m.rows; i++) {
        tr += gretl_matrix_get(m, i, i);
    }

    return tr;
}

/*
 * gretl_matrix_random_fill:
 * @m: input matrix.
 * @dist: either DistCode.uniform or DistCode.normal.
 *
 * Fills @m with pseudo-random values from either the uniform
 * or the standard normal distribution.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_matrix_random_fill(GretlMatrixRaw * m, DistCode dist)
{
    int n;

    if (m == null || (dist != DistCode.uniform && dist != DistCode.normal)) {
        return 1;
    }

    n = m.rows * m.cols;

    if (n > 0) {
        if (dist == DistCode.normal) {
            gretl_rand_normal(m.val, 0, n - 1);
        } else if (dist == DistCode.uniform) {
            gretl_rand_uniform(m.val, 0, n - 1);
        }
    }

    return 0;
}

/* Solves a*x = b for triangular @a: on exit @b is overwritten
   by @x
*/

private int gretl_triangular_solve(const(GretlMatrixRaw) * a,
                                   GretlMatrixRaw * b,
                                   GretlMatrixMod mod,
                                   GretlMatrixStructure t)
{
    char uplo, transa;
    char side = 'L';
    char diag = 'N';
    double alpha = 1.0;
    integer m = b.rows;
    integer n = b.cols;

    uplo = (t == GretlMatrixStructure.lower_triangular) ? 'L' : 'U';
    transa = (mod == GretlMatrixMod.transpose) ? 'T' : 'N';

    dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha,
           cast(double *) a.val, &m, b.val, &m);

    return 0;
}

/* Fills the (column) p-vector @X with a single realization of p normal
   variates with a covariance structure specified by the lower-triangular
   p x p matrix L. If prec = 0, L is assumed to pertain to a covariance
   matrix, otherwise it's assumed to pertain to a precision matrix.
*/

int correlated_normal_vec(GretlVector * X, const(GretlMatrixRaw) * L,
                          int prec)
{
    gretl_matrix_random_fill(X, DistCode.normal);

    if (prec == 0) {
        /* covariance matrix: pre-multiply X by L */
        gretl_blas_dtrmm(L, X, "LLN");
    } else {
        /* precision matrix: do left division, L' \ X */
        gretl_triangular_solve(L, X, GretlMatrixMod.transpose,
                               GretlMatrixStructure.lower_triangular);
    }

    return 0;
}

/*
 * gretl_random_matrix_new:
 * @r: number of rows.
 * @c: number of columns.
 * @dist: either DistCode.uniform or DistCode.normal.
 *
 * Creates a new @r x @c matrix and fills it with pseudo-random
 * values from either the uniform or the standard normal
 * distribution.
 *
 * Returns: allocated matrix or null on failure.
 */

GretlMatrixRaw * gretl_random_matrix_new(int r, int c, DistCode dist)
{
    GretlMatrixRaw * m = null;

    if (dist != DistCode.uniform && dist != DistCode.normal) {
        return null;
    } else if (r < 0 || c < 0) {
        return null;
    } else if (r == 0 || c == 0) {
        m = gretl_null_matrix_new();
        if (m != null) {
            m.rows = r;
            m.cols = c;
        }
    } else {
        m = gretl_matrix_alloc(r, c);
        if (m != null) {
            if (dist == DistCode.normal) {
                gretl_rand_normal(m.val, 0, r * c - 1);
            } else if (dist == DistCode.uniform) {
                gretl_rand_uniform(m.val, 0, r * c - 1);
            }
        }
    }

    return m;
}

/*
 * gretl_vector_mean:
 * @v: input vector.
 *
 * Returns: the arithmetic mean of the elements of @v, or
 * NADBL on failure.
 */

double gretl_vector_mean(const(GretlVector) * v)
{
    double num = 0.0;
    int i, n, den = 0;

    if (gretl_is_null_matrix(v)) {
        return NADBL;
    }

    n = gretl_vector_get_length(v);
    if (n == 0) {
        return NADBL;
    }

    for (i=0; i<n; i++) {
        if (!na(v.val[i])) {
            num += v.val[i];
            den++;
        }
    }

    return (den > 0) ? (num / den) : NADBL;
}

/*
 * gretl_vector_variance:
 * @v: input vector.
 *
 * Returns: the variance of the elements of @v, or
 * NADBL on failure.
 */

double gretl_vector_variance(const(GretlVector) * v)
{
    double s2 = 0.0;
    double x, xbar = 0.0;
    int i, n, den = 0;

    if (gretl_is_null_matrix(v)) {
        return NADBL;
    }

    n = gretl_vector_get_length(v);
    if (n == 0) {
        return NADBL;
    }

    for (i=0; i<n; i++) {
        if (!na(v.val[i])) {
            xbar += v.val[i];
            den++;
        }
    }

    if (den == 0) {
        return NADBL;
    }

    xbar /= den;

    for (i=0; i<n; i++) {
        x = v.val[i];
        if (!na(x)) {
            x -= xbar;
            s2 += x * x;
        }
    }

    return s2 / den;
}

private int real_matrix_resample(GretlMatrixRaw * R, const(GretlMatrixRaw) * m)
{
    int i, j, k, t1, r = R.rows;
    int * z = cast(int *) malloc(r * int.sizeof);
    double x;

    if (z == null) {
        return GretlError.alloc;
    }

    /* generate r drawings from [0 .. r-1] */
    gretl_rand_int_minmax(z, r, 0, r - 1);

    /* sample from source matrix @m based on row indices */
    for (i=0; i<r; i++) {
        k = z[i] % m.rows;
        for (j=0; j<m.cols; j++) {
            x = gretl_matrix_get(m, k, j);
            gretl_matrix_set(R, i, j, x);
        }
    }

    t1 = gretl_matrix_get_t1(m);
    if (t1 > 0 && r <= m.rows) {
        gretl_matrix_set_t1(R, t1);
        gretl_matrix_set_t2(R, t1 + r - 1);
    }

    free(z);

    return 0;
}

/*
 * gretl_matrix_resample:
 * @m: input matrix.
 * @draws: number of draws (or 0 to use the number or rows
 * in @m).
 * @err: location to receive error code.
 *
 * Returns: a new matrix consisting of a random re-sampling
 * (with replacement) of the rows of @m, or null on
 * failure.
 */

GretlMatrixRaw * gretl_matrix_resample(const(GretlMatrixRaw) * m,
                                    int draws, int * err)
{
    GretlMatrixRaw * R = null;
    int r;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
        return null;
    } else if (m.is_complex) {
        *err = GretlError.cmplx;
        return null;
    }

    if (draws < 0) {
        *err = GretlError.invarg;
        return null;
    } else if (draws > 0) {
        r = draws;
    } else {
        r = m.rows;
    }

    R = gretl_matrix_alloc(r, m.cols);

    if (R == null) {
        *err = GretlError.alloc;
    } else {
        *err = real_matrix_resample(R, m);
    }

    return R;
}

int gretl_matrix_resample2(GretlMatrixRaw * targ,
                           const(GretlMatrixRaw) * src)
{
    if (gretl_is_null_matrix(targ) || gretl_is_null_matrix(src)) {
        return GretlError.data;
    } else if (targ.is_complex || src.is_complex) {
        return GretlError.cmplx;
    } else {
        return real_matrix_resample(targ, src);
    }
}

/*
 * gretl_matrix_block_resample:
 * @m: input matrix.
 * @blocklen: length of moving blocks.
 * @draws: number of draws (or 0 to use the rows of @m).
 * @err: location to receive error code.
 *
 * Returns: a new matrix consisting of a random re-sampling
 * (with replacement) of the rows of @m, using blocks of
 * contiguous rows of length @blocklen, or null on
 * failure.
 */

GretlMatrixRaw * gretl_matrix_block_resample(const(GretlMatrixRaw) * m,
                                          int blocklen, int draws,
                                          int * err)
{
    GretlMatrixRaw * R = null;
    int * z = null;
    double x;
    int b, n, s, r, rmax;
    int t1;
    int i, j, k;

    if (gretl_is_null_matrix(m) || blocklen <= 0 || draws < 0) {
        *err = GretlError.data;
        return null;
    } else if (m.is_complex) {
        *err = GretlError.cmplx;
        return null;
    }

    if (blocklen == 1) {
        return gretl_matrix_resample(m, draws, err);
    }

    r = draws > 0 ? draws : m.rows;

    /* Let n represent the number of blocks of @blocklen
       contiguous rows which we need to select; the
       last of these may not be fully used.
    */
    n = r / blocklen + (r % blocklen > 0);

    rmax = m.rows - blocklen;
    if (rmax < 0) {
        *err = GretlError.data;
        return null;
    }

    R = gretl_matrix_alloc(r, m.cols);
    z = cast(int *) malloc(n * (*z).sizeof);

    if (R == null || z == null) {
        gretl_matrix_free(R);
        free(z);
        *err = GretlError.alloc;
        return null;
    }

    /* generate n drawings from [0 .. rmax] */
    gretl_rand_int_minmax(z, n, 0, rmax);

    /* sample from source matrix based on block indices */
    i = 0;
    for (b=0; b<n; b++) {
        for (s=0; s<blocklen; s++) {
            if (i < r) {
                /* don't spill over the end */
                k = z[b] + s;
                for (j=0; j<m.cols; j++) {
                    x = gretl_matrix_get(m, k, j);
                    gretl_matrix_set(R, i, j, x);
                }
                i++;
            } else {
                break;
            }
        }
    }

    t1 = gretl_matrix_get_t1(m);
    if (t1 > 0 && r <= m.rows) {
        gretl_matrix_set_t1(R, t1);
        gretl_matrix_set_t2(R, t1 + r - 1);
    }

    free(z);

    return R;
}

/*
 * gretl_matrix_block_resample2:
 * @targ: target matrix.
 * @src: source matrix.
 * @blocklen: length of moving blocks.
 * @z: array of length XXX.
 *
 * An "in-place" version of gretl_matrix_block_resample().
 * It is assumed that @targ is a matrix of the same dimensions
 * as @src, that @blocklen is greater than 1, and that @z
 * is long enough to hold n integers, where n is the number
 * of rows in @src divided by @blocklen, rounded up to the
 * nearest integer.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_matrix_block_resample2(GretlMatrixRaw * targ,
                                 const(GretlMatrixRaw) * src,
                                 int blocklen,
                                 int * z)
{
    double x;
    int r = src.rows;
    int b, n, s, rmax;
    int i, j, k;

    n = r / blocklen + (r % blocklen > 0);

    rmax = r - blocklen;
    if (rmax < 0) {
        return GretlError.data;
    }

    /* generate n drawings from [0 .. rmax] */
    gretl_rand_int_minmax(z, n, 0, rmax);

    /* sample from source matrix based on block indices */
    i = 0;
    for (b=0; b<n; b++) {
        for (s=0; s<blocklen; s++) {
            if (i < r) {
                k = z[b] + s;
                for (j=0; j<src.cols; j++) {
                    x = gretl_matrix_get(src, k, j);
                    gretl_matrix_set(targ, i, j, x);
                }
                i++;
            } else {
                break;
            }
        }
    }

    return 0;
}

private int gretl_matrix_zero_triangle(GretlMatrixRaw * m, char t)
{
    int i, j;

    if (gretl_is_null_matrix(m)) {
        return GretlError.data;
    } else if (m.rows != m.cols) {
        return GretlError.nonconf;
    }

    if (t == 'U') {
        for (i=0; i<m.rows-1; i++) {
            for (j=i+1; j<m.cols; j++) {
                if (m.is_complex) {
                    gretl_cmatrix_set(m, i, j, Cmplx(0.0, 0.0));
                } else {
                    gretl_matrix_set(m, i, j, 0.0);
                }
            }
        }
    } else {
        for (i=1; i<m.rows; i++) {
            for (j=0; j<i; j++) {
                if (m.is_complex) {
                    gretl_cmatrix_set(m, i, j, Cmplx(0.0, 0.0));
                } else {
                    gretl_matrix_set(m, i, j, 0.0);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_zero_upper:
 * @m: square matrix to operate on.
 *
 * Sets the elements of @m outside of the lower triangle to zero.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_matrix_zero_upper(GretlMatrixRaw * m)
{
    return gretl_matrix_zero_triangle(m, 'U');
}

/*
 * gretl_matrix_zero_lower:
 * @m: square matrix to operate on.
 *
 * Sets the elements of @m outside of the upper triangle to zero.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_matrix_zero_lower(GretlMatrixRaw * m)
{
    return gretl_matrix_zero_triangle(m, 'L');
}

/*
 * gretl_matrix_multiply_by_scalar:
 * @m: matrix to operate on.
 * @x: scalar by which to multiply.
 *
 * Multiplies all elements of @m by @x.
 */

void gretl_matrix_multiply_by_scalar(GretlMatrixRaw * m, double x)
{
    int i, n = m.rows * m.cols;

    for (i=0; i<n; i++) {
        m.val[i] *= x;
    }
}

/*
 * gretl_matrix_divide_by_scalar:
 * @m: matrix to operate on.
 * @x: scalar by which to divide.
 *
 * Divides all elements of @m by @x.
 *
 * Returns: 0 on success, 1 if x = 0.
 */

int gretl_matrix_divide_by_scalar(GretlMatrixRaw * m, double x)
{
    if (x == 0.0) {
        return 1;
    } else {
        gretl_matrix_multiply_by_scalar(m, 1.0 / x);
        return 0;
    }
}

/*
 * gretl_matrix_switch_sign:
 * @m: matrix to operate on.
 *
 * Changes the sign of each element of @m.
 */

void gretl_matrix_switch_sign(GretlMatrixRaw * m)
{
    if (!gretl_is_null_matrix(m)) {
        int i, n = m.rows * m.cols;

        for (i=0; i<n; i++) {
            m.val[i] = -m.val[i];
        }
    }
}

/*
 * gretl_matrix_raise:
 * @m: matrix to operate on.
 * @x: exponent.
 *
 * Raises each element of @m to the power @x.
 */

void gretl_matrix_raise(GretlMatrixRaw * m, double x)
{
    if (!gretl_is_null_matrix(m)) {
        int i, n = m.rows * m.cols;

        for (i=0; i<n; i++) {
            m.val[i] = pow(m.val[i], x);
        }
    }
}

/* if x can be represented as an integer and is an exact power of 2,
   return the exact log_2 of x, otherwise use log() to compute an
   approximation.
*/

private double log_2(double x)
{
    int i, s;

    if (x <= 0) {
        return log(x);
    }

    if (floor(x) != x || x < 2 || x > cast(double) int.max) {
        return log2(x);
    }

    s = cast(int) floor(x);

    for (i=1; ; i++) {
        if (s % 2) {
            break;
        }
        s /= 2;
        if (s == 1) {
            return cast(double) i;
        }
    }

    return log2(x);
}

private double mexp_error_eps(int q)
{
    double x1, x2, x3, x4;
    double qf = x_factorial(q);

    x1 = pow(2.0, 3.0 - (q + q));
    x2 = qf * qf;
    x3 = x_factorial(2 * q);
    x4 = (2 * q + 1) * x3;
    x3 *= x4;

    return x1 * (x2 / x3);
}

/*
 * gretl_matrix_exp:
 * @m: square matrix to operate on.
 * @err: location to receive error code.
 *
 * Calculates the matrix exponential of @m, using algorithm
 * 11.3.1 from Golub and Van Loan, "Matrix Computations", 3e.
 *
 * Returns: the exponential, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_exp(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * A = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * N = null;
    GretlMatrixRaw * D = null;
    GretlMatrixRaw * W = null;
    double xa, c, j, delta = 1.0e-13;
    int q, k, n;

    if (gretl_is_null_matrix(m) || m.rows != m.cols) {
        *err = GretlError.data;
        return null;
    }

    n = m.rows;

    A = gretl_matrix_copy_tmp(m);
    X = gretl_identity_matrix_new(n);
    N = gretl_identity_matrix_new(n);
    D = gretl_identity_matrix_new(n);
    W = gretl_matrix_alloc(n, n);

    if (A == null || X == null || N == null ||
        D == null || W == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    xa = gretl_matrix_infinity_norm(A);

    j = floor(log_2(xa));
    if (j < 0) {
        j = 0;
    }

    gretl_matrix_divide_by_scalar(A, pow(2.0, j));

    for (q=1; q<16; q++) {
        c = mexp_error_eps(q);
        if (c * xa <= delta) {
            break;
        }
    }

    c = 1.0;

    for (k=1; k<=q; k++) {
        c *= (q - k + 1.0) / ((2.0*q - k + 1) * k);
        /* X = AX */
        gretl_matrix_multiply(A, X, W);
        gretl_matrix_copy_values(X, W);
        /* N = N + cX */
        gretl_matrix_multiply_by_scalar(W, c);
        gretl_matrix_add_to(N, W);
        /* D = D + (-1)^k cX */
        if (k % 2) {
            gretl_matrix_subtract_from(D, W);
        } else {
            gretl_matrix_add_to(D, W);
        }
    }

    /* solve DF = N for F */
    *err = gretl_LU_solve(D, N);

    if (!*err) {
        for (k=0; k<j; k++) {
            gretl_matrix_multiply(N, N, W);
            gretl_matrix_copy_values(N, W);
        }
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(X);
    gretl_matrix_free(D);
    gretl_matrix_free(W);

    if (*err) {
        gretl_matrix_free(N);
        N = null;
    }

    return N;
}

/*
 * gretl_matrix_polroots:
 * @a: vector of coefficients.
 * @force_complex: see below.
 * @legacy: see below.
 * @err: location to receive error code.
 *
 * Calculates the roots of the polynomial with coefficients
 * given by @a.  If the degree of the polynomial is p, then
 * @a should contain p + 1 coefficients in ascending order,
 * i.e. starting with the constant and ending with the
 * coefficient on x^p.
 *
 * As a transitional measure, if @legacy is non-zero, in the
 * case of complex roots the return vector is in old-style
 * format: real values in column 0 and imaginary parts in
 * column 1.
 *
 * Returns: by default, a regular p-vector if all the roots are real,
 * otherwise a p x 1 complex matrix. The @force_complex flag can
 * be used to produce a complex return value even if the imaginary
 * parts are all zero.
 */

GretlMatrixRaw * gretl_matrix_polroots(const(GretlMatrixRaw) * a,
                                    int force_complex,
                                    int legacy,
                                    int * err)
{
    GretlMatrixRaw * r = null;
    double * work = null;
    Cmplx * roots = null;
    double * xcof, cof;
    int i, m, order, polerr;

    m = gretl_vector_get_length(a);
    if (m < 2) {
        *err = GretlError.data;
        return null;
    }

    order = m - 1;
    work = cast(double *) malloc(2 * m * (*work).sizeof);
    roots = cast(Cmplx *) malloc(order * (*roots).sizeof);

    if (work == null || roots == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    xcof = work;
    cof = xcof + m;

    for (i=0; i<m; i++) {
        xcof[i] = a.val[i];
    }

    polerr = polrt(xcof, cof, order, cast(CephesCmplx *) roots);

    if (polerr) {
        *err = GretlError.data;
    } else {
        int allreal = !force_complex;
        Cmplx z;

        for (i=0; i<order && allreal; i++) {
            if (roots[i].i != 0) {
                allreal = 0;
            }
        }
        if (allreal) {
            r = gretl_matrix_alloc(order, 1);
        } else if (legacy) {
            r = gretl_matrix_alloc(order, 2);
        } else {
            r = gretl_cmatrix_new(order, 1);
        }
        if (r == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0; i<order; i++) {
                if (allreal) {
                    gretl_matrix_set(r, i, 0, roots[i].r);
                } else if (legacy) {
                    gretl_matrix_set(r, i, 0, roots[i].r);
                    gretl_matrix_set(r, i, 1, roots[i].i);
                } else {
                    z = Cmplx(roots[i].r, roots[i].i);
                    gretl_cmatrix_set(r, i, 0, z);
                }
            }
        }
    }

 bailout:

    free(work);
    free(roots);

    return r;
}

/*
 * gretl_vector_copy_values:
 * @targ: target vector.
 * @src: source vector.
 *
 * Copies the elements of @src into the corresponding elements
 * of @targ.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * two vectors are not of the same length.
 */

int gretl_vector_copy_values(GretlVector * targ,
                             const(GretlVector) * src)
{
    int n;

    if (src == null) {
        fprintf(stderr, "gretl_vector_copy_values: src is NULL\n");
        return GretlError.data;
    }

    if (targ == src) {
        /* no-op */
        return 0;
    }

    n = gretl_vector_get_length(src);

    if (gretl_vector_get_length(targ) != n) {
        return GretlError.nonconf;
    }

    if (n > 0) {
        memcpy(targ.val, src.val, n * (*targ.val).sizeof);
    }

    return 0;
}

/*
 * gretl_matrix_copy_values:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies the elements of @src into the corresponding elements
 * of @targ.
 *
 * Returns: 0 on successful completion, or
 * GretlError.nonconf if the two matrices are not
 * conformable for the operation.
 */

int gretl_matrix_copy_values(GretlMatrixRaw * targ,
                             const(GretlMatrixRaw) * src)
{
    int n;

    if (src == null) {
        fprintf(stderr, "gretl_matrix_copy_values: src is NULL\n");
        return GretlError.data;
    } else if (targ == src) {
        /* no-op */
        return 0;
    } else if (targ.is_complex + src.is_complex == 1) {
        return GretlError.mixed;
    }

    if (targ.rows != src.rows || targ.cols != src.cols) {
        fprintf(stderr, "gretl_matrix_copy_values: targ is %d x %d but src is %d x %d\n",
                targ.rows, targ.cols, src.rows, src.cols);
        return GretlError.nonconf;
    }

    n = src.rows * src.cols;
    if (n > 0) {
        if (src.is_complex) {
            n *= 2;
        }
        memcpy(targ.val, src.val, n * (*targ.val).sizeof);
    }

    return 0;
}

/*
 * gretl_matrix_copy_data:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies all data from @src to @targ: this does the same
 * as gretl_matrix_copy_values() but also transcribes
 * t1, t2, colnames and rownames information, if present.
 *
 * Returns: 0 on successful completion, non-zero code on
 * failure.
 */

int gretl_matrix_copy_data(GretlMatrixRaw * targ,
                           const(GretlMatrixRaw) * src)
{
    int err;

    err = gretl_matrix_copy_values(targ, src);

    if (!err) {
        err = gretl_matrix_copy_info(targ, src, GretlMatrixMod.none);
    }

    return err;
}

/*
 * gretl_matrix_copy_values_shaped:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Copies the elements of @src into @targ, column
 * by column.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * the two matrices do not contain the same number of
 * elements.
 */

int gretl_matrix_copy_values_shaped(GretlMatrixRaw * targ,
                                    const(GretlMatrixRaw) * src)
{
    int n = targ.rows * targ.cols;

    if (src.rows * src.cols != n) {
        fprintf(stderr, "gretl_matrix_copy_values_shaped: "
                ~ "targ is %d x %d but src is %d x %d\n",
                targ.rows, targ.cols, src.rows, src.cols);
        return GretlError.nonconf;
    }

    if (n > 0) {
        memcpy(targ.val, src.val, n * (*targ.val).sizeof);
    }

    return 0;
}

private int add_scalar_to_matrix(GretlMatrixRaw * targ, double x)
{
    int i, n = targ.rows * targ.cols;

    for (i=0; i<n; i++) {
        targ.val[i] += x;
    }

    return 0;
}

private int subtract_scalar_from_matrix(GretlMatrixRaw * targ, double x)
{
    int i, n = targ.rows * targ.cols;

    for (i=0; i<n; i++) {
        targ.val[i] -= x;
    }

    return 0;
}

/*
 * gretl_matrix_add_to:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Adds the elements of @src to the corresponding elements
 * of @targ.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * two matrices are not conformable for the operation.
 * In the special case where @src is in fact a scalar, the
 * operation always goes through OK, with the scalar being
 * added to each element of @targ.
 */

int
gretl_matrix_add_to(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int i, n;

    if (targ.rows != src.rows || targ.cols != src.cols) {
        if (matrix_is_scalar(src)) {
            return add_scalar_to_matrix(targ, src.val[0]);
        } else {
            fprintf(stderr, "gretl_matrix_add_to: adding %d x %d to %d x %d\n",
                    src.rows, src.cols, targ.rows, targ.cols);
            return GretlError.nonconf;
        }
    }

    n = src.rows * src.cols;

    for (i=0; i<n; i++) {
        targ.val[i] += src.val[i];
    }

    return 0;
}

/*
 * gretl_matrix_add:
 * @a: source matrix.
 * @b: source matrix.
 * @c: target matrix.
 *
 * Adds the elements of @a and @b, the result going to @c.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * matrices are not conformable for the operation.
 */

int
gretl_matrix_add(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b,
                 GretlMatrixRaw * c)
{
    int rows = a.rows, cols = a.cols;
    int i, n;

    if (a.is_complex || b.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_add\n");
        return GretlError.cmplx;
    } else if (b.rows != rows || c.rows != rows ||
        b.cols != cols || c.cols != cols) {
        fprintf(stderr, "gretl_matrix_add: non-conformable\n");
        return GretlError.nonconf;
    }

    n = rows * cols;

    for (i=0; i<n; i++) {
        c.val[i] = a.val[i] + b.val[i];
    }

    return 0;
}

/*
 * gretl_matrix_add_transpose_to:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Adds the elements of @src, transposed, to the corresponding
 * elements of @targ.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * two matrices are not conformable for the operation.
 */

int gretl_matrix_add_transpose_to(GretlMatrixRaw * targ,
                                  const(GretlMatrixRaw) * src)
{
    int i, j, k = 0;

    if (targ.is_complex || src.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_add_transpose_to\n");
        return GretlError.cmplx;
    } else if (targ.rows != src.cols || targ.cols != src.rows) {
        fprintf(stderr, "gretl_matrix_add_transpose_to: "
                ~ "adding %d x %d to %d x %d\n",
                src.cols, src.rows, targ.rows, targ.cols);
        return GretlError.nonconf;
    }

    /* note: the k index follows column-major order */
    for (i=0; i<src.rows; i++) {
        for (j=0; j<src.cols; j++) {
            targ.val[k++] += gretl_matrix_get(src, i, j);
        }
    }

    return 0;
}

/*
 * gretl_matrix_subtract_from:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Subtracts the elements of @src from the corresponding elements
 * of @targ.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * two matrices are not conformable for the operation.
 * In the special case where @src is in fact a scalar, the
 * operation always goes through OK, with the scalar being
 * subtracted from each element of @targ.
 */

int
gretl_matrix_subtract_from(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int i, n;

    if (targ.is_complex || src.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_subtract_from\n");
        return GretlError.cmplx;
    } else if (targ.rows != src.rows || targ.cols != src.cols) {
        if (matrix_is_scalar(src)) {
            return subtract_scalar_from_matrix(targ, src.val[0]);
        } else {
            return GretlError.nonconf;
        }
    }

    n = src.rows * src.cols;

    for (i=0; i<n; i++) {
        targ.val[i] -= src.val[i];
    }

    return 0;
}

/*
 * gretl_matrix_subtract:
 * @a: source_matrix.
 * @b: source matrix.
 * @c: target matrix.
 *
 * Subtracts the elements of @b from the corresponding elements
 * of @a, the result going to @c.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * matrices are not conformable for the operation.
 */

int
gretl_matrix_subtract(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b,
                      GretlMatrixRaw * c)
{
    int rows = a.rows, cols = a.cols;
    int i, n;

    if (a.is_complex || b.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_subtract\n");
        return GretlError.cmplx;
    } else if (b.rows != rows || c.rows != rows ||
        b.cols != cols || c.cols != cols) {
        fprintf(stderr, "gretl_matrix_subtract: non-conformable\n");
        return GretlError.nonconf;
    }

    n = rows * cols;

    for (i=0; i<n; i++) {
        c.val[i] = a.val[i] - b.val[i];
    }

    return 0;
}

/*
 * gretl_matrix_subtract_reversed:
 * @a: m x n matrix.
 * @b: m x n matrix.
 *
 * Operates on @b such that b_{ij} = a_{ij} - b_{ij}.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if the
 * two matrices are not conformable for the operation.
 */

int
gretl_matrix_subtract_reversed(const(GretlMatrixRaw) * a, GretlMatrixRaw * b)
{
    int i, n;

    if (a.rows != b.rows || a.cols != b.cols) {
        return GretlError.nonconf;
    }

    n = a.rows * a.cols;

    for (i=0; i<n; i++) {
        b.val[i] = a.val[i] - b.val[i];
    }

    return 0;
}

/*
 * gretl_matrix_I_minus:
 * @m: original square matrix, n x n.
 *
 * Rewrites @m as (I - m), where I is the n x n identity
 * matrix.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if @m is
 * not square.
 */

int gretl_matrix_I_minus(GretlMatrixRaw * m)
{
    double x;
    int n = m.rows;
    int i, j, k = 0;

    if (m.cols != n) {
        return GretlError.nonconf;
    }

    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            x = m.val[k];
            if (i == j) {
                m.val[k] = 1.0 - x;
            } else if (x != 0.0) {
                m.val[k] = -x;
            }
            k++;
        }
    }

    return 0;
}

/*
 * gretl_matrix_inscribe_I:
 * @m: original matrix.
 * @row: top row for insertion.
 * @col: leftmost row for insertion.
 * @n: dimension (rows and columns) of identity matrix.
 *
 * Writes an n x n identity matrix into matrix @m, the top left-hand
 * corner of the insertion being given by @row and @col (which are
 * 0-based).
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if an identity
 * matrix of the specified size cannot be fitted into @m at the
 * specified location.
 */

int gretl_matrix_inscribe_I(GretlMatrixRaw * m, int row, int col, int n)
{
    int i, j, cj, ri;

    if (n < 0) {
        return GretlError.invarg;
    } else if (n == 0) {
        /* no-op */
        return 0;
    }

    if (row < 0 || row + n > m.rows ||
        col < 0 || col + n > m.cols) {
        return GretlError.nonconf;
    }

    for (j=0; j<n; j++) {
        cj = col + j;
        for (i=0; i<n; i++) {
            ri = row + i;
            gretl_matrix_set(m, ri, cj, (i == j) ? 1.0 : 0.0);
        }
    }

    return 0;
}

/*
 * gretl_matrix_transpose_in_place:
 * @m: matrix to transpose.
 *
 * Tranposes @m in place.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_matrix_transpose_in_place(GretlMatrixRaw * m)
{
    int r = m.rows;
    int c = m.cols;
    int i, j;

    gretl_matrix_destroy_info(m);

    if (r == 1 || c == 1) {
        m.cols = r;
        m.rows = c;
        return 0;
    }

    if (r == c) {
        double mij, mji;
        int n = r - 1;

        for (i=0; i<n; i++) {
            for (j=i+1; j<c; j++) {
                mij = gretl_matrix_get(m, i, j);
                mji = gretl_matrix_get(m, j, i);
                gretl_matrix_set(m, i, j, mji);
                gretl_matrix_set(m, j, i, mij);
            }
        }
    } else {
        size_t sz = r * c * double.sizeof;
        double * val = cast(double *) mval_malloc(sz);
        int k = 0;

        if (val == null) {
            return GretlError.alloc;
        }

        memcpy(val, m.val, sz);
        m.rows = c;
        m.cols = r;
        for (j=0; j<c; j++) {
            for (i=0; i<r; i++) {
                gretl_matrix_set(m, j, i, val[k++]);
            }
        }

        mval_free(val);
    }

    return 0;
}

/*
 * gretl_matrix_transpose:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * Fills out @targ (which must be pre-allocated and of the right
 * dimensions) with the transpose of @src.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_matrix_transpose(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int r = src.rows;
    int c = src.cols;
    int i, j, k = 0;

    if (targ.rows != c || targ.cols != r) {
        return GretlError.nonconf;
    }

    for (j=0; j<c; j++) {
        for (i=0; i<r; i++) {
            gretl_matrix_set(targ, j, i, src.val[k++]);
        }
    }

    return 0;
}

/*
 * gretl_square_matrix_transpose:
 * @m: square matrix to operate on.
 *
 * Transposes the matrix @m.
 *
 * Returns: 0 on success, non-zero error code otherwise.
 */

int gretl_square_matrix_transpose(GretlMatrixRaw * m)
{
    double x, y;
    int mij, mji;
    int i, j;

    if (m.rows != m.cols) {
        fputs("gretl_square_matrix_transpose: matrix must be square\n",
              stderr);
        return 1;
    }

    for (i=0; i<m.rows-1; i++) {
        for (j=i+1; j<m.rows; j++) {
            mij = mdx(m,i,j);
            mji = mdx(m,j,i);
            x = m.val[mij];
            y = m.val[mji];
            m.val[mij] = y;
            m.val[mji] = x;
        }
    }

    return 0;
}

/*
 * gretl_matrix_xtr_symmetric:
 * @m: gretl_matrix.
 *
 * Computes the symmetric part of @m by averaging its off-diagonal
 * elements.
 */

void gretl_matrix_xtr_symmetric(GretlMatrixRaw * m)
{
    double x;
    int mij, mji;
    int i, j;

    for (i=0; i<m.rows; i++) {
        for (j=0; j<i; j++) {
            mij = mdx(m,i,j);
            mji = mdx(m,j,i);
            x = m.val[mij];
            x += m.val[mji];
            m.val[mij] = m.val[mji] = 0.5 * x;
        }
    }
}

/*
 * gretl_matrix_add_self_transpose:
 * @m: (square) matrix to operate on.
 *
 * Adds the transpose of @m to @m itself, yielding a symmetric
 * matrix.
 *
 * Returns: 0 on successful completion, or
 * 1 if the source matrix is not square.
 */

int gretl_matrix_add_self_transpose(GretlMatrixRaw * m)
{
    double x;
    int mij, mji;
    int i, j;

    if (m.rows != m.cols) {
        fputs("gretl_matrix_add_self_transpose: matrix must be square\n",
              stderr);
        return GretlError.nonconf;
    }

    for (i=0; i<m.rows; i++) {
        for (j=i; j<m.rows; j++) {
            mij = mdx(m,i,j);
            mji = mdx(m,j,i);
            x = m.val[mij];
            x += m.val[mji];
            m.val[mij] = m.val[mji] = x;
        }
    }

    return 0;
}

/*
 * gretl_matrix_vectorize:
 * @targ: target vector, (m * n) x 1.
 * @src: source matrix, m x n.
 *
 * Writes into @targ vec(@src), that is, a column vector
 * formed by stacking the columns of @src.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int
gretl_matrix_vectorize(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int n;

    if (gretl_is_null_matrix(src) || gretl_is_null_matrix(targ)) {
        return GretlError.data;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    n = src.rows * src.cols;

    if (targ.cols != 1 || targ.rows != n) {
        return GretlError.nonconf;
    }

    if (src.is_complex) {
        n *= 2;
    }
    memcpy(targ.val, src.val, n * (*src.val).sizeof);

    return 0;
}

/*
 * gretl_matrix_vectorize_new:
 * @m: matrix to be vectorized.
 *
 * Returns: a gretl column vector, vec(@m), or null on failure.
 */

GretlMatrixRaw * gretl_matrix_vectorize_new(const(GretlMatrixRaw) * m)
{
    GretlMatrixRaw * v;
    int n;

    if (gretl_is_null_matrix(m)) {
        return null;
    }

    n = m.rows * m.cols;

    v = gretl_matching_matrix_new(n, 1, m);

    if (v != null) {
        if (m.is_complex) {
            n *= 2;
        }
        memcpy(v.val, m.val, n * (*m.val).sizeof);
    }

    return v;
}

/*
 * gretl_matrix_unvectorize:
 * @targ: target matrix, m x n.
 * @src: source vector, (m * n) x 1.
 *
 * Writes successive blocks of length m from @src into
 * the successive columns of @targ (that is, performs the
 * inverse of the vec() operation).
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int
gretl_matrix_unvectorize(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int n;

    if (gretl_is_null_matrix(src) || gretl_is_null_matrix(targ)) {
        return GretlError.data;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    n = targ.rows * targ.cols;

    if (src.cols != 1 || src.rows != n) {
        return GretlError.nonconf;
    }

    if (src.is_complex) {
        n *= 2;
    }
    memcpy(targ.val, src.val, n * (*src.val).sizeof);

    return 0;
}

/*
 * gretl_matrix_vectorize_h:
 * @targ: target vector, length n * (n+1)/2.
 * @src: source square matrix, n x n.
 *
 * Writes into @targ vech(@src), that is, a vector
 * containing the lower-triangular elements of @src.
 * This is only useful for symmetric matrices, but for the
 * sake of performance we don't check for that.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int
gretl_matrix_vectorize_h(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int n = src.rows;
    int m = n * (n+1) / 2;
    int i, j, k;

    if (gretl_vector_get_length(targ) != m) {
        return GretlError.nonconf;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    k = 0;
    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            if (src.is_complex) {
                targ.z[k++] = gretl_cmatrix_get(src, i, j);
            } else {
                targ.val[k++] = gretl_matrix_get(src, i, j);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_vectorize_h_skip:
 * @targ: target vector, length n * (n-1)/2.
 * @src: source square matrix, n x n.
 *
 * Like gretl_matrix_vectorize_h() except that the diagonal
 * of @src is skipped in the vectorization.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int
gretl_matrix_vectorize_h_skip(GretlMatrixRaw * targ,
                              const(GretlMatrixRaw) * src)
{
    int n = src.rows;
    int m = n * (n-1) / 2;
    int i, j, k;

    if (gretl_vector_get_length(targ) != m) {
        return GretlError.nonconf;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    k = 0;
    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            if (i == j) {
                continue;
            }
            if (src.is_complex) {
                targ.z[k++] = gretl_cmatrix_get(src, i, j);
            } else {
                targ.val[k++] = gretl_matrix_get(src, i, j);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_unvectorize_h:
 * @targ: target matrix, n x n.
 * @src: source vector, m x 1.
 *
 * Rearranges successive blocks of decreasing length from @src into
 * the successive columns of @targ (that is, performs the
 * inverse of the vech() operation): @targ comes out symmetric.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int
gretl_matrix_unvectorize_h(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int m = gretl_vector_get_length(src);
    int n = targ.rows;
    Cmplx z;
    double x;
    int i, j, k;

    if (m == 0 || n * (n + 1) != 2 * m) {
        return GretlError.nonconf;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    k = 0;
    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            if (src.is_complex) {
                z = src.z[k++];
                gretl_cmatrix_set(targ, i, j, conj(z));
                gretl_cmatrix_set(targ, j, i, z);
            } else {
                x = src.val[k++];
                gretl_matrix_set(targ, i, j, x);
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_unvectorize_h_diag:
 * @targ: target matrix, n x n.
 * @src: source vector.
 * @diag: value for filling the diagonal.
 *
 * Like gretl_matrix_unvectorize_h(), but expects input @src in which
 * the diagonal elements are omitted, while @diag gives the value
 * to put in place of the omitted elements.
 *
 * Returns: 0 on successful completion, or GretlError.nonconf if
 * @targ is not correctly dimensioned.
 */

int gretl_matrix_unvectorize_h_diag(GretlMatrixRaw * targ,
                                    const(GretlMatrixRaw) * src,
                                    double diag)
{
    int m = gretl_vector_get_length(src);
    int n = targ.rows;
    Cmplx z;
    double x;
    int i, j, k;

    if (m == 0 || n * (n - 1) != 2 * m) {
        return GretlError.nonconf;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    }

    k = 0;
    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            if (i == j) {
                if (src.is_complex) {
                    z = Cmplx(diag, 0);
                    gretl_cmatrix_set(targ, i, j, z);
                } else {
                    gretl_matrix_set(targ, i, j, diag);
                }
            } else if (src.is_complex) {
                z = src.z[k++];
                gretl_cmatrix_set(targ, i, j, conj(z));
                gretl_cmatrix_set(targ, j, i, z);
            } else {
                x = src.val[k++];
                gretl_matrix_set(targ, i, j, x);
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_inscribe_matrix:
 * @targ: target matrix.
 * @src: source matrix.
 * @row: row offset for insertion (0-based).
 * @col: column offset for insertion.
 * @mod: GretlMatrixMod.none, GretlMatrixMod.transpose or
 * GretlMatrixMod.cumulate.
 *
 * Writes @src into @targ, starting at offset @row, @col.
 * The @targ matrix must be large enough to sustain the
 * inscription of @src at the specified point.  If @mod
 * is GretlMatrixMod.transpose it is in fact the transpose of
 * @src that is written into @targ. Or if @mod is
 * GretlMatrixMod.cumulate the incoming values from @src are
 * added to those already present in @targ, rather than
 * replacing them.
 *
 * Returns: 0 on success, GretlError.nonconf if the matrices are
 * not conformable for the operation.
 */

int gretl_matrix_inscribe_matrix(GretlMatrixRaw * targ,
                                 const(GretlMatrixRaw) * src,
                                 int row, int col,
                                 GretlMatrixMod mod)
{
    int m = (mod == GretlMatrixMod.transpose) ? src.cols : src.rows;
    int n = (mod == GretlMatrixMod.transpose) ? src.rows : src.cols;
    double * tval = null;
    const(double) * sval = null;
    size_t sz = 0;
    Cmplx z;
    double x;
    int i, j, ri, cj;

    if (row < 0 || col < 0) {
        return GretlError.nonconf;
    } else if (targ.is_complex + src.is_complex == 1) {
        return GretlError.mixed;
    } else if (row + m > targ.rows ||
        col + n > targ.cols) {
        fprintf(stderr, "gretl_matrix_inscribe_matrix: out of bounds\n");
        return GretlError.nonconf;
    }

    if (!targ.is_complex && mod == GretlMatrixMod.none) {
        /* we'll use memcpy below */
        tval = targ.val + col * targ.rows + row;
        sval = src.val;
        sz = src.rows * (*sval).sizeof;
    }

    for (j=0; j<n; j++) {
        cj = col + j;
        if (tval != null) {
            /* do the column in one go */
            memcpy(tval, sval, sz);
            tval += targ.rows;
            sval += src.rows;
            continue;
        }
        for (i=0; i<m; i++) {
            /* proceed row by row */
            ri = row + i;
            if (src.is_complex) {
                if (mod == GretlMatrixMod.transpose) {
                    z = cmatrix_transp_get(src, i, j);
                } else {
                    z = gretl_cmatrix_get(src, i, j);
                    if (mod == GretlMatrixMod.cumulate) {
                        z += gretl_cmatrix_get(targ, ri, cj);
                    }
                }
                gretl_cmatrix_set(targ, ri, cj, z);
            } else {
                if (mod == GretlMatrixMod.transpose) {
                    x = matrix_transp_get(src, i, j);
                } else {
                    /* we must be cumulating */
                    x = gretl_matrix_get(targ, ri, cj) +
                        gretl_matrix_get(src, i, j);
                }
                gretl_matrix_set(targ, ri, cj, x);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_extract_matrix:
 * @targ: target matrix.
 * @src: source matrix.
 * @row: row offset for extraction (0-based).
 * @col: column offset for extraction.
 * @mod: either GretlMatrixMod.transpose or GretlMatrixMod.none.
 *
 * Writes into @targ a sub-matrix of @src, taken from the
 * offset @row, @col.  The @src matrix must be large enough
 * to provide a sub-matrix of the dimensions of @targ.
 * If @mod is GretlMatrixMod.transpose it is in fact the transpose
 * of the sub-matrix that that is written into @targ.
 *
 * Returns: 0 on success, GretlError.nonconf if the matrices are
 * not conformable for the operation.
 */

int gretl_matrix_extract_matrix(GretlMatrixRaw * targ,
                                const(GretlMatrixRaw) * src,
                                int row, int col,
                                GretlMatrixMod mod)
{
    int m = (mod == GretlMatrixMod.transpose) ? targ.cols : targ.rows;
    int n = (mod == GretlMatrixMod.transpose) ? targ.rows : targ.cols;
    Cmplx z;
    double x;
    int i, j, si, sj;

    if (row < 0 || col < 0) {
        return GretlError.nonconf;
    } else if (src.is_complex + targ.is_complex == 1) {
        return GretlError.mixed;
    } else if (row >= src.rows) {
        fprintf(stderr, "extract_matrix: requested starting row=%d, but "
                ~ "src has %d rows\n", row, src.rows);
        return GretlError.nonconf;
    } else if (col >= src.cols) {
        fprintf(stderr, "extract_matrix: requested starting col=%d, but "
                ~ "src has %d cols\n", col, src.cols);
        return GretlError.nonconf;
    } else if (row + m > src.rows || col + n > src.cols) {
        fprintf(stderr, "gretl_matrix_extract_matrix: out of bounds\n");
        return GretlError.nonconf;
    }

    si = row;
    for (i=0; i<m; i++) {
        sj = col;
        for (j=0; j<n; j++) {
            if (src.is_complex) {
                z = gretl_cmatrix_get(src, si, sj++);
                if (mod == GretlMatrixMod.transpose) {
                    cmatrix_transp_set(targ, i, j, z);
                } else {
                    gretl_cmatrix_set(targ, i, j, z);
                }
            } else {
                x = gretl_matrix_get(src, si, sj++);
                if (mod == GretlMatrixMod.transpose) {
                    matrix_transp_set(targ, i, j, x);
                } else {
                    gretl_matrix_set(targ, i, j, x);
                }
            }
        }
        si++;
    }

    return 0;
}

/*
 * gretl_matrix_steal_data:
 * @m: matrix to operate on.
 *
 * "Steals" the allocated data from @m, which is left with a
 * null data pointer.
 *
 * Returns: a pointer to the "stolen" data.
 */

double * gretl_matrix_steal_data(GretlMatrixRaw * m)
{
    double * vals = null;

    if (m != null) {
        if (is_block_matrix(m)) {
            matrix_block_error("gretl_matrix_steal_data");
            return null;
        }
        vals = m.val;
        m.val = null;
        m.z = null;
    }

    return vals;
}

private void real_gretl_matrix_print(const(GretlMatrixRaw) * m,
                                     const(char) * msg,
                                     FILE * fp)
{
    const(char) * fmt = "%#12.5g ";
    char * envstr;
    int i, j;

    if (m == null) {
        if (msg != null && *msg != '\0') {
            fprintf(fp, "%s: matrix is NULL\n", msg);
        } else {
            fputs("matrix is NULL\n", fp);
        }
        return;
    } else if (m.val == null) {
        if (msg != null && *msg != '\0') {
            fprintf(fp, "%s: matrix is empty, %d x %d\n", msg,
                    m.rows, m.cols);
        } else {
            fprintf(fp, "matrix is empty: %d x %d\n",
                    m.rows, m.cols);
        }
        return;
    }

    if (m.is_complex) {
        Prn * prn = gretl_print_new_with_stream(fp);

        if (prn != null) {
            gretl_cmatrix_print(m, msg, prn);
            gretl_print_destroy(prn);
        }
        return;
    }

    envstr = getenv("GRETL_MATRIX_DEBUG");
    if (envstr != null && atoi(envstr) > 0) {
        fmt = "%#22.15g ";
    } else {
        envstr = getenv("GRETL_MATRIX_PRINT6");
        if (envstr != null && atoi(envstr) > 0) {
            fmt = "%#12.6g ";
        }
    }

    if (msg != null && *msg != '\0') {
        fprintf(fp, "%s (%d x %d)", msg, m.rows, m.cols);
        if (is_block_matrix(m)) {
            fprintf(fp, " (part of matrix block)\n\n");
        } else if (gretl_matrix_is_dated(m)) {
            int mt1 = gretl_matrix_get_t1(m);
            int mt2 = gretl_matrix_get_t2(m);

            fprintf(fp, " [t1 = %d, t2 = %d]\n\n", mt1 + 1, mt2 + 1);
        } else {
            fputs("\n\n", fp);
        }
    }

    for (i=0; i<m.rows; i++) {
        for (j=0; j<m.cols; j++) {
            fprintf(fp, fmt, gretl_matrix_get(m, i, j));
        }
        fputc('\n', fp);
    }

    fputc('\n', fp);
}

/*
 * gretl_matrix_print:
 * @m: matrix.
 * @msg: message to print with matrix, or null.
 *
 * Prints the given matrix to stderr.
 */

void gretl_matrix_print(const(GretlMatrixRaw) * m, const(char) * msg)
{
    real_gretl_matrix_print(m, msg, stderr);
}

/*
 * gretl_matrix_print2:
 * @m: matrix.
 * @msg: message to print with matrix, or null.
 *
 * Prints the given matrix to stdout.
 */

void gretl_matrix_print2(const(GretlMatrixRaw) * m, const(char) * msg)
{
    real_gretl_matrix_print(m, msg, stdout);
}

private enum DEFAULT_EQTOL = 1.0e-9; /* 2014-08-05: was 1.5e-12 */

private __gshared double eq_tol = DEFAULT_EQTOL;

/*
 * gretl_matrix_set_equals_tolerance:
 * @tol: tolerance value.
 *
 * Sets the tolerance for judging whether or not a matrix is symmetric
 * (see gretl_matrix_is_symmetric() and also
 * gretl_invert_symmetric_matrix()).  The tolerance is the maximum
 * relative difference between corresponding off-diagonal elements that
 * is acceptable in a supposedly "symmetric" matrix.  The default
 * value is 1.0e-9.
 */

void gretl_matrix_set_equals_tolerance(double tol)
{
    eq_tol = tol;
}

/*
 * gretl_matrix_unset_equals_tolerance:
 *
 * Sets the tolerance for judging whether or not a matrix is symmetric
 * to its default value.  See also gretl_matrix_set_equals_tolerance().
 */

void gretl_matrix_unset_equals_tolerance()
{
    eq_tol = DEFAULT_EQTOL;
}

private double sneq_reldiff(double x, double y, double ad)
{
    double rd;

    if (x == 0.0) {
        rd = y;
    } else if (y == 0.0) {
        rd = x;
    } else if (x > y) {
        rd = ad / y;
    } else {
        rd = ad / x;
    }

    return fabs(rd);
}

private int real_gretl_matrix_is_symmetric(const(GretlMatrixRaw) * m,
                                           int verbose)
{
    double x, y, ad, rd;
    int i, j;

    if (gretl_is_null_matrix(m)) {
        return 0;
    }

    for (i=1; i<m.rows; i++) {
        for (j=0; j<i; j++) {
            x = gretl_matrix_get(m, i, j);
            y = gretl_matrix_get(m, j, i);
            ad = fabs(y - x);
            if (ad < 1.0e-12) {
                continue;
            }
            rd = sneq_reldiff(x, y, ad);
            if (rd > eq_tol) {
                if (verbose) {
                    fprintf(stderr, "M(%d,%d) = %.16g but M(%d,%d) = %.16g\n"
                            ~ " reldiff = %g\n", i, j, x, j, i, y, rd);
                    if (m.rows < 100) {
                        gretl_matrix_print(m, "gretl_matrix_is_symmetric()");
                    }
                }
                return 0;
            }
        }
    }

    return 1;
}

/*
 * gretl_matrix_is_symmetric:
 * @m: gretl_matrix.
 *
 * Returns: 1 if @m is symmetric (with a small relative tolerance
 * for asymmetry), otherwise 0.
 */

int gretl_matrix_is_symmetric(const(GretlMatrixRaw) * m)
{
    return real_gretl_matrix_is_symmetric(m, 0);
}

/*
 * gretl_matrix_is_idempotent:
 * @m: gretl_matrix.
 * @tol: numerical tolerance
 *
 * Returns: 1 if @m is idempotent, otherwise 0.
 */

int gretl_matrix_is_idempotent(const(GretlMatrixRaw) * m, double tol)
{
    GretlMatrixRaw * b;
    int k, ret, err;

    if (gretl_is_null_matrix(m)) {
        return 0;
    }

    k = m.rows;

    if (m.cols != k) {
        return 0;
    }

    b = gretl_matrix_alloc(k, k);
    if (b == null) {
        return 0;
    }

    gretl_matrix_multiply(m, m, b);
    ret = gretl_matrices_are_equal(m, b, tol, &err);
    gretl_matrix_free(b);

    return ret;
}

/*
 * gretl_matrix_infinity_norm:
 * @m: gretl_matrix.
 *
 * Returns: the infinity-norm of @m (the maximum value across
 * the rows of @m of the sum of the absolute values of
 * the elements in the given row).
 */

double gretl_matrix_infinity_norm(const(GretlMatrixRaw) * m)
{
    double rsum, rmax = 0.0;
    int i, j;

    if (gretl_is_null_matrix(m)) {
        return NADBL;
    }

    for (i=0; i<m.rows; i++) {
        rsum = 0.0;
        for (j=0; j<m.cols; j++) {
            if (m.is_complex) {
                rsum += cabs(gretl_cmatrix_get(m, i, j));
            } else {
                rsum += fabs(gretl_matrix_get(m, i, j));
            }
        }
        if (rsum > rmax) {
            rmax = rsum;
        }
    }

    return rmax;
}

/*
 * gretl_matrix_one_norm:
 * @m: gretl_matrix.
 *
 * Returns: the 1-norm of @m (the maximum value across
 * the columns of @m of the sum of the absolute values of
 * the elements in the given column).
 */

double gretl_matrix_one_norm(const(GretlMatrixRaw) * m)
{
    double csum, cmax = 0.0;
    int i, j;

    if (gretl_is_null_matrix(m)) {
        return NADBL;
    }

    for (j=0; j<m.cols; j++) {
        csum = 0.0;
        for (i=0; i<m.rows; i++) {
            if (m.is_complex) {
                csum += cabs(gretl_cmatrix_get(m, i, j));
            } else {
                csum += fabs(gretl_matrix_get(m, i, j));
            }
        }
        if (csum > cmax) {
            cmax = csum;
        }
    }

    return cmax;
}

/*
 * gretl_vcv_log_determinant:
 * @m: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the log determinant of the symmetric positive-definite
 * matrix @m using Cholesky decomposition.
 *
 * Returns: the log determinant, or NADBL on failure.
 */

double gretl_vcv_log_determinant(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * a = null;
    char uplo = 'L';
    integer n, info;
    double det = NADBL;
    int i;

    if (gretl_is_null_matrix(m)) {
        return NADBL;
    }

    n = m.rows;

    if (m.rows != m.cols) {
        fputs("gretl_vcv_log_determinant: matrix must be square\n", stderr);
        *err = GretlError.invarg;
        return det;
    }

    if (!real_gretl_matrix_is_symmetric(m, 1)) {
        fputs("gretl_vcv_log_determinant: matrix is not symmetric\n", stderr);
        *err = GretlError.invarg;
        return det;
    }

    a = gretl_matrix_copy_tmp(m);
    if (a == null) {
        fputs("gretl_vcv_log_determinant: out of memory\n", stderr);
        *err = GretlError.alloc;
        return det;
    }

    dpotrf_(&uplo, &n, a.val, &n, &info);

    if (info != 0) {
        if (info > 0) {
            *err = GretlError.notpd;
        } else {
            fputs("gretl_vcv_log_determinant: illegal argument to dpotrf\n",
                  stderr);
            *err = GretlError.invarg;
        }
    } else {
        double x;

        det = 1.0;
        for (i=0; i<n; i++) {
            x = gretl_matrix_get(a, i, i);
            det *= x * x;
        }
        det = log(det);
    }

    gretl_matrix_free(a);

    return det;
}

/* This is really only necessary when using OpenBLAS: when the matrix
   under analysis contains NaN values the pivot values calculated by
   dgetrf() can go out of bounds. We could flag an error here if we
   find an out-of-bounds value, but the downside of that is that we'd
   get different results when using OpenBLAS versus netlib lapack/blas
   (since netlib returns a NaN matrix rather than erroring out).

   This check added 2015-12-24, required for OpenBLAS 0.2.16.dev and
   earlier.
*/

private void pivot_check(integer * ipiv, int n)
{
    int i;

    for (i=0; i<n; i++) {
        if (ipiv[i] > n) {
            /* clamp the bad value to avoid a crash */
            fprintf(stderr, "pivot_check: clamped bad ipiv[%d] = %d\n",
                    i, ipiv[i]);
            ipiv[i] = n;
        }
    }
}

/* Calculate the determinant of @a using LU factorization.  If logdet
   != 0 and absval == 0, return the log of the determinant, or NA if
   the determinant is non-positive.  if logdet != 0 and absval != 0,
   return the log of the absolute value of the determinant. Otherwise
   return the determinant itself.
*/

private double gretl_LU_determinant(GretlMatrixRaw * a, int logdet,
                                    int absval, int * err)
{
    integer n, info;
    integer * ipiv;
    double det;
    int i;

    if (gretl_is_null_matrix(a)) {
        *err = GretlError.data;
        return NADBL;
    }

    *err = 0;
    n = a.rows;
    if (a.cols != n) {
        fputs("gretl_LU_determinant: matrix must be square\n", stderr);
        *err = GretlError.nonconf;
        return NADBL;
    }

    if (n == 1) {
        /* simple 1 x 1 case */
        det = a.val[0];
        if (logdet) {
            if (det > 0) {
                return log(det);
            } else if (det < 0) {
                return absval ? log(-det) : NADBL;
            } else {
                return NADBL;
            }
        } else {
            return det;
        }
    }

    ipiv = cast(integer *) malloc(n * (*ipiv).sizeof);
    if (ipiv == null) {
        *err = GretlError.alloc;
        return NADBL;
    }

    dgetrf_(&n, &n, a.val, &n, ipiv, &info);

    if (info > 0) {
        if (logdet) {
            return NADBL;
        } else {
            return 0;
        }
    } else if (info < 0) {
        fprintf(stderr, "gretl_LU_determinant: dgetrf gave info = %d\n",
                cast(int) info);
        free(ipiv);
        *err = GretlError.data;
        return NADBL;
    } else {
        pivot_check(ipiv, n);
    }

    if (logdet) {
        int negcount = 0;
        double di;

        /* Note: we're better off here taking logs and adding, rather
           than multiplying terms then taking the log of the product.
           In this way we can get a finite result for the log-determinant
           of a matrix whose determinant is numerically "infinite" --
           up to a point.
        */
        det = 0.0;
        for (i=0; i<n; i++) {
            di = gretl_matrix_get(a, i, i);
            if (di == 0.0) {
                fputs("gretl_matrix_log_determinant: determinant = 0\n", stderr);
                det = NADBL;
                break;
            }
            if (ipiv[i] != i + 1) {
                di = -di;
            }
            if (di < 0) {
                di = -di;
                negcount++;
            }
            det += log(di);
        }
        if (!absval && negcount % 2) {
            /* got a negative value: try calculating the determinant
               itself, for reporting
            */
            double d = 1.0;

            for (i=0; i<n; i++) {
                di = gretl_matrix_get(a, i, i);
                if (ipiv[i] != i + 1) {
                    di = -di;
                }
                d *= di;
            }
            fprintf(stderr, "gretl_matrix_log_determinant: determinant is < 0 (%g)\n", d);
            det = NADBL;
        }
    } else {
        /* plain determinant */
        det = 1.0;
        for (i=0; i<n; i++) {
            if (ipiv[i] != i + 1) {
                det = -det;
            }
            det *= gretl_matrix_get(a, i, i);
        }
    }

    free(ipiv);

    return det;
}

private double det_22(const(double) * a, int * err)
{
    return a[0]*a[3] - a[1]*a[2];
}

private double det_33(const(double) * a, int * err)
{
    double d = a[0]*a[4]*a[8] - a[0]*a[7]*a[5]
        + a[3]*a[7]*a[2] - a[3]*a[1]*a[8]
        + a[6]*a[1]*a[5] - a[6]*a[4]*a[2];

    return d;
}

/*
 * gretl_matrix_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the determinant of the square matrix @a using the LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 *
 * Returns: the determinant, or NADBL on failure.
 */

double gretl_matrix_determinant(GretlMatrixRaw * a, int * err)
{
    if (a != null) {
        if (a.rows == 2 && a.cols == 2) {
            return det_22(a.val, err);
        } else if (a.rows == 3 && a.cols == 3) {
            return det_33(a.val, err);
        }
    }

    return gretl_LU_determinant(a, 0, 0, err);
}

/*
 * gretl_matrix_log_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the log of the determinant of the square matrix @a using LU
 * factorization.  Matrix @a is not preserved: it is overwritten
 * by the factorization.
 *
 * Returns: the determinant, or NADBL on failure.
 */

double gretl_matrix_log_determinant(GretlMatrixRaw * a, int * err)
{
    return gretl_LU_determinant(a, 1, 0, err);
}

/*
 * gretl_matrix_log_abs_determinant:
 * @a: gretl_matrix.
 * @err: location to receive error code.
 *
 * Compute the log of the absolute value of the determinant of the
 * square matrix @a using LU factorization.  Matrix @a is not
 * preserved: it is overwritten by the factorization.
 *
 * Returns: the determinant, or NADBL on failure.
 */

double gretl_matrix_log_abs_determinant(GretlMatrixRaw * a, int * err)
{
    return gretl_LU_determinant(a, 1, 1, err);
}

private void matrix_grab_content(GretlMatrixRaw * targ, GretlMatrixRaw * src)
{
    targ.rows = src.rows;
    targ.cols = src.cols;
    targ.is_complex = src.is_complex;

    mval_free(targ.val);
    targ.val = src.val;
    targ.z = src.z;
    src.val = null;
    src.z = null;

    gretl_matrix_destroy_info(targ);
    targ.info = src.info;
    src.info = null;
}

/* Least squares solution using QR, with column pivoting and detection
   of rank deficiency, using lapack dgelsy, or zgelsy if the arguments
   are complex.
*/

private int QR_solve(GretlMatrixRaw * A, GretlMatrixRaw * B)
{
    integer m, n, lda, nrhs;
    integer rank, info = 0;
    integer lwork = -1;
    integer * jpvt = null;
    double * work = null;
    double * rwork = null;
    double rcond = QR_RCOND_MIN;
    int zfunc = 0;
    int wsz = 1;
    int i, err = 0;

    lda = m = A.rows;
    n = A.cols;
    nrhs = B.cols;

    if (m > n && is_block_matrix(B)) {
        matrix_block_error("QR solve");
        return GretlError.data;
    }

    if (n > m || B.rows != m) {
        return GretlError.nonconf;
    }

    if (A.is_complex) {
        if (!B.is_complex) {
            return GretlError.invarg;
        }
        zfunc = 1;
        wsz = 2;
    }

    jpvt = cast(integer *) malloc(n * (*jpvt).sizeof);
    work = cast(double *) lapack_malloc(wsz * (*work).sizeof);
    if (jpvt == null || work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (zfunc) {
        rwork = cast(double *) malloc(2 * n * (*rwork).sizeof);
        if (rwork == null) {
            err = GretlError.alloc;
            goto bailout;
        }
    }

    for (i=0; i<n; i++) {
        jpvt[i] = 0;
    }

    /* workspace query */
    if (zfunc) {
        zgelsy_(&m, &n, &nrhs, cast(Cmplx *) A.z, &lda, cast(Cmplx *) B.z, &lda,
                jpvt, &rcond, &rank, cast(Cmplx *) work, &lwork, rwork, &info);
    } else {
        dgelsy_(&m, &n, &nrhs, A.val, &lda, B.val, &lda,
                jpvt, &rcond, &rank, work, &lwork, &info);
    }
    if (info != 0) {
        fprintf(stderr, "gelsy: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

    /* optimally sized work array */
    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, cast(size_t) lwork * wsz * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* run actual computation */
    if (zfunc) {
        zgelsy_(&m, &n, &nrhs, cast(Cmplx *) A.z, &lda, cast(Cmplx *) B.z, &lda,
                jpvt, &rcond, &rank, cast(Cmplx *) work, &lwork, rwork, &info);
    } else {
        dgelsy_(&m, &n, &nrhs, A.val, &lda, B.val, &lda,
                jpvt, &rcond, &rank, work, &lwork, &info);
    }
    if (info != 0) {
        fprintf(stderr, "gelsy: info = %d\n", cast(int) info);
        err = 1;
    } else if (rank < n) {
        fprintf(stderr, "gelsy: cols(A) = %d, rank(A) = %d\n",
                A.cols, rank);
    }

    if (!err && m > n) {
        GretlMatrixRaw * C;

        C = gretl_matrix_trim_rows(B, 0, m - n, &err);
        if (!err) {
            matrix_grab_content(B, C);
            gretl_matrix_free(C);
        }
    }

 bailout:

    free(jpvt);
    lapack_free(work);
    if (zfunc) {
        free(rwork);
    }

    return err;
}

/*
 * gretl_LU_solve_invert:
 * @a: square matrix.
 * @b: matrix.
 *
 * Solves ax = b for the unknown x, using LU decomposition,
 * then proceeds to use the decomposition to invert @a. Calls
 * the LAPACK functions dgetrf(), dgetrs() and dgetri(); the
 * decomposition proceeds via partial pivoting with row
 * interchanges.
 *
 * On exit, @b is replaced by the solution and @a is replaced
 * by its inverse.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_LU_solve_invert(GretlMatrixRaw * a, GretlMatrixRaw * b)
{
    char trans = 'N';
    integer info;
    integer n, ldb, nrhs = 1;
    integer lwork = -1;
    double * work = null;
    integer * ipiv;
    int err = 0;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        a.rows != a.cols) {
        return GretlError.data;
    }

    n = a.rows;

    if (b.cols == 1) {
        ldb = b.rows;
    } else if (b.rows == 1) {
        ldb = b.cols;
    } else {
        nrhs = b.cols;
        ldb = b.rows;
    }

    ipiv = cast(integer *) malloc(n * (*ipiv).sizeof);
    if (ipiv == null) {
        return GretlError.alloc;
    }

    dgetrf_(&n, &n, a.val, &n, ipiv, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_LU_solve_invert: dgetrf gave info = %d\n",
                cast(int) info);
        err = (info < 0) ? GretlError.data : GretlError.singular;
    } else {
        pivot_check(ipiv, n);
    }

    if (!err) {
        dgetrs_(&trans, &n, &nrhs, a.val, &n, ipiv, b.val, &ldb, &info);
        if (info != 0) {
            fprintf(stderr, "gretl_LU_solve_invert: dgetrs gave info = %d\n",
                    cast(int) info);
            err = GretlError.data;
        }
    }

    if (!err) {
        work = cast(double *) lapack_malloc((*work).sizeof);
        if (work == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        dgetri_(&n, a.val, &n, ipiv, work, &lwork, &info);
        if (info != 0) {
            err = wspace_fail(info, work[0]);
        } else {
            lwork = cast(integer) work[0];
            work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
            if (work == null) {
                err = GretlError.alloc;
            }
        }
    }

    if (!err) {
        dgetri_(&n, a.val, &n, ipiv, work, &lwork, &info);
        if (info != 0) {
            fprintf(stderr, "gretl_LU_solve_invert: dgetri gave info = %d\n",
                    cast(int) info);
            err = GretlError.data;
        }
    }

    free(ipiv);
    lapack_free(work);

    return err;
}

/*
 * gretl_LU_solve:
 * @a: square matrix.
 * @b: matrix.
 *
 * Solves ax = b for the unknown x, via LU decomposition
 * using partial pivoting with row interchanges.
 * On exit, @b is replaced by the solution and @a is replaced
 * by its decomposition. Calls the LAPACK functions dgetrf()
 * and dgetrs().
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_LU_solve(GretlMatrixRaw * a, GretlMatrixRaw * b)
{
    char trans = 'N';
    integer info;
    integer n, ldb, nrhs = 1;
    integer * ipiv;
    int zfunc = 0;
    int debug_ = 0;
    int err = 0;

    /* (BLASDEBUG is 1 in the C source, so this block is live) */
    const(char) * envstr = getenv("GRETL_MATRIX_DEBUG");

    debug_ = (envstr != null && atoi(envstr) > 0);

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        a.rows != a.cols) {
        return GretlError.data;
    }

    zfunc = a.is_complex;
    if (zfunc && !b.is_complex) {
        return GretlError.invarg;
    }

    if (debug_) {
        fputs("gretl_LU_solve\n", stderr);
        gretl_matrix_print(a, "a, on input");
        gretl_matrix_print(b, "b, on input");
    }

    n = a.cols;

    if (b.cols == 1) {
        ldb = b.rows;
    } else if (b.rows == 1) {
        ldb = b.cols;
    } else {
        nrhs = b.cols;
        ldb = b.rows;
    }

    ipiv = cast(integer *) malloc(n * (*ipiv).sizeof);
    if (ipiv == null) {
        return GretlError.alloc;
    }

    if (zfunc) {
        zgetrf_(&n, &n, cast(Cmplx *) a.val, &n, ipiv, &info);
    } else {
        dgetrf_(&n, &n, a.val, &n, ipiv, &info);
    }

    if (info != 0) {
        fprintf(stderr, "gretl_LU_solve: getrf gave info = %d\n",
                cast(int) info);
        err = (info < 0) ? GretlError.data : GretlError.singular;
    } else {
        pivot_check(ipiv, n);
    }

    if (!err) {
        if (zfunc) {
            zgetrs_(&trans, &n, &nrhs, cast(Cmplx *) a.val, &n, ipiv,
                    cast(Cmplx *) b.val, &ldb, &info);
        } else {
            dgetrs_(&trans, &n, &nrhs, a.val, &n, ipiv, b.val, &ldb, &info);
        }
        if (info != 0) {
            fprintf(stderr, "gretl_LU_solve: dgetrs gave info = %d\n",
                    cast(int) info);
            err = GretlError.data;
        }
    }

    if (debug_) {
        gretl_matrix_print(a, "a, on return");
        gretl_matrix_print(b, "b, on return");
        fprintf(stderr, "err, on return = %d\n", err);
    }

    free(ipiv);

    return err;
}

/*
 * gretl_matrix_solve:
 * @a: m x n matrix, with m >= n.
 * @b: gretl_matrix.
 *
 * Solves ax = b for the unknown x. If @a is square the method
 * is LU decomposition, on which see also gretl_LU_solve().
 * If m > n the QR decomposition is used to find the least
 * squares solution.
 *
 * On exit, @b is replaced by the solution and @a is replaced
 * by its decomposition.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

private int gretl_matrix_solve(GretlMatrixRaw * a, GretlMatrixRaw * b)
{
    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        return GretlError.data;
    }

    if (a.rows == a.cols) {
        return gretl_LU_solve(a, b);
    } else if (a.rows > a.cols) {
        return QR_solve(a, b);
    } else {
        return GretlError.data;
    }
}

private enum CHOL_TINY  = 8.0e-09;
private enum CHOL_SMALL = 1.0e-08;

/* Native Cholesky decomposition and back-solution, with a check for
   excessively small diagonal elements.  The matrix @a is the vech
   of X'X, and b is the right-hand side on entry, the solution on
   successful exit.  This may not be as efficient as lapack's
   dpotrf/dpotrs, but the advantage is that this function flags
   near-singularity whereas the LAPACK functions generate an error
   condition only on outright singularity.

   Note that the matrix @a is overwritten.
*/

private int native_cholesky_decomp_solve(GretlMatrixRaw * a, GretlMatrixRaw * b)
{
    double * xtx = a.val;
    double * xty = b.val;
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;
    int nc = b.rows;

    if (xtx[0] <= 0.0) {
        fprintf(stderr, "%s %d: xtx <= 0.0\n", __FILE__.ptr, __LINE__);
        return GretlError.nan;
    }

    e = 1.0 / sqrt(xtx[0]);
    xtx[0] = e;
    xty[0] *= e;
    for (i=1; i<nc; i++) {
        xtx[i] *= e;
    }

    kk = nc;

    for (j=1; j<nc; j++) {
        /* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j;

        for (l=1; l<=jm1; l++) {
            xx = xtx[k];
            d1 += xx * xty[l-1];
            d += xx * xx;
            k += nc-l;
        }

        d2 = xtx[kk] - d;
        test = d2 / xtx[kk];

        /* check for effective singularity */
        if (test < CHOL_TINY) {
            fprintf(stderr, "cholesky: test[%d] = %g\n", j, test);
            return GretlError.singular;
        } else if (test < CHOL_SMALL) {
            fprintf(stderr, "cholesky: test[%d] = %g\n", j, test);
        }

        e = 1 / sqrt(d2);
        xtx[kk] = e;
        xty[j] = (xty[j] - d1) * e;

        /* off-diagonal elements */
        for (i=j+1; i<nc; i++) {
            kk++;
            d = 0.0;
            k = j;
            for (l=1; l<=jm1; l++) {
                d += xtx[k] * xtx[k-j+i];
                k += nc - l;
            }
            xtx[kk] = (xtx[kk] - d) * e;
        }
        kk++;
    }

    kk--;

    /* back-solve for the coefficients, into b */

    xty[nc-1] *= xtx[kk];

    for (j=nc-2; j>=0; j--) {
        d = xty[j];
        for (i=nc-1; i>j; i--) {
            d -= xty[i] * xtx[--kk];
        }
        xty[j] = d * xtx[--kk];
    }

    for (j=0; j<nc; j++) {
        if (isNaN(xty[j])) {
            fprintf(stderr, "%s %d: coeff %d is NaN\n", __FILE__.ptr, __LINE__, j);
            return GretlError.nan;
        }
    }

    return 0;
}

private enum CHOL_RCOND_MIN = 1.0e-6;

/*
 * gretl_cholesky_decomp_solve:
 * @a: symmetric positive-definite matrix.
 * @b: matrix 'x' on input, solution 'b' on output.
 *
 * Solves ax = b for the unknown matrix x, using Cholesky decomposition
 * via the LAPACK functions dpotrf() and dpotrs().
 *
 * On exit, @b is replaced by the solution and @a is replaced by its
 * Cholesky decomposition.
 *
 * Returns: 0 on successful completion, or non-zero code on error.
 */

int gretl_cholesky_decomp_solve(GretlMatrixRaw * a, GretlMatrixRaw * b)
{
    integer n, m, info = 0;
    double rcond;
    double * work = null;
    integer * iwork = null;
    char diag = 'N';
    char norm = '1';
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        return GretlError.data;
    }

    n = a.cols;
    m = b.cols;

    dpotrf_(&uplo, &n, a.val, &n, &info);
    if (info != 0) {
        fprintf(stderr, "gretl_cholesky_decomp_solve: "
                ~ "dpotrf failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
        err = (info > 0) ? GretlError.notpd : GretlError.data;
    }

    if (!err) {
        work = cast(double *) lapack_malloc(3 * n * (*work).sizeof);
        iwork = cast(integer *) malloc(n * (*iwork).sizeof);
        if (work == null || iwork == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        dtrcon_(&norm, &uplo, &diag, &n, a.val, &n, &rcond, work, iwork, &info);
        if (rcond < CHOL_RCOND_MIN) {
            err = GretlError.singular;
        }
    }

    if (!err) {
        dpotrs_(&uplo, &n, &m, a.val, &n, b.val, &n, &info);
        if (info != 0) {
            fprintf(stderr, "gretl_cholesky_decomp_solve:\n"
                    ~ " dpotrs failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
            err = GretlError.singular;
        }
    }

    lapack_free(work);
    free(iwork);

    return err;
}

/*
 * gretl_cholesky_solve:
 * @a: Lower triangular Cholesky factor of symmetric p.d. matrix
 * @b: right-hand side vector.
 *
 * Solves ax = b for the unknown x, using the precomputed
 * Cholesky factor in @a. On exit, @b is replaced by the solution.
 *
 * Returns: 0 on successful completion, or non-zero code on error.
 */

int gretl_cholesky_solve(const(GretlMatrixRaw) * a, GretlVector * b)
{
    integer n, info, one = 1;
    char uplo = 'L';

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        return GretlError.data;
    }

    n = a.cols;

    dpotrs_(&uplo, &n, &one, cast(double *) a.val, &n, b.val, &n, &info);
    if (info != 0) {
        fprintf(stderr, "gretl_cholesky_solve:\n"
                ~ " dpotrs failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
        return GretlError.singular;
    }

    return 0;
}

/*
 * gretl_cholesky_invert:
 * @a: Cholesky-decomposed symmetric positive-definite matrix.
 *
 * Inverts @a, which must contain the pre-computed Cholesky
 * decomposition of a p.d. matrix, as may be obtained using
 * gretl_cholesky_decomp_solve(). For speed there's no error
 * checking of the input -- the caller should make sure it's OK.
 *
 * Returns: 0 on successful completion, or non-zero code on error.
 */

int gretl_cholesky_invert(GretlMatrixRaw * a)
{
    integer info, n = a.cols;
    char uplo = 'L';
    int err = 0;

    dpotri_(&uplo, &n, a.val, &n, &info);

    if (info != 0) {
        err = GretlError.singular;
        fprintf(stderr, "gretl_cholesky_invert:\n"
                ~ " dpotri failed with info = %d\n", cast(int) info);
    } else {
        gretl_matrix_mirror(a, uplo);
    }

    return err;
}

/* translation to C of tsld1.f in the netlib toeplitz package,
   code as of 07/23/82; see http://www.netlib.org/toeplitz/

   tsld1 solves the double precision linear system
   A * x = b for Toeplitz matrix A

   on entry:
     a1     double precision(m), the first row of A
     a2     double precision(m - 1), the first column of A
            beginning with the second element
      b     double precision(m), the right hand side vector
     c1     double precision(m - 1), workspace
     c2     double precision(m - 1), workspace
      m     integer, order of the matrix A
     pd     integer, flag if positive definiteness is assumed

     (c1 and c2 are internalized below)

   on exit:
      x     double precision(m), the solution vector
      dt    pointer to double, determinant if pd == 0, log determinant otherwise
*/

private enum TOEPLITZ_SMALL = 1.0e-30;

private int tsld1(const(double) * a1, const(double) * a2,
                  const(double) * b, double * x, double * dt,
                  int m, int pd)
{
    double r1, r2, r3, r5, r6;
    int n, i, n1;
    double * c1 = null;
    double * c2 = null;
    int err = 0;

    if (fabs(a1[0]) < TOEPLITZ_SMALL) {
        return GretlError.singular;
    }

    /* solve the system with principal minor of order 1 */

    r1 = a1[0];
    if (dt != null) {
        *dt = pd ? log(r1) : r1;
    }

    x[0] = b[0] / r1;
    if (m == 1) {
        return 0;
    }

    c1 = cast(double *) malloc((m-1) * (*c1).sizeof);
    c2 = cast(double *) malloc((m-1) * (*c2).sizeof);

    if (c1 == null || c2 == null) {
        free(c1);
        free(c2);
        return GretlError.alloc;
    }

    r2 = 0.0;

    /* recurrent process for solving the system for
       order = 2 to m */

    for (n=1; n<m; n++) {
        /* compute multiples of the first and last columns of
           the inverse of the principal minor of order n + 1
        */
        n1 = n - 1;
        r5 = a2[n1];
        r6 = a1[n];
        if (n > 1) {
            c1[n1] = r2;
            for (i=0; i<n1; i++) {
                r5 += a2[i] * c1[n1-i];
                r6 += a1[i+1] * c2[i];
            }
        }

        r2 = -r5 / r1;
        r3 = -r6 / r1;
        r1 += r5 * r3;

        if (fabs(r1) < TOEPLITZ_SMALL) {
            err = GretlError.singular;
            goto bailout;
        } else if (dt != null) {
            if (pd) {
                *dt += log(r1);
            } else {
                *dt *= r1;
            }
        }

        if (n > 1) {
            r6 = c2[0];
            c2[n1] = 0.0;
            for (i=1; i<n; i++) {
                r5 = c2[i];
                c2[i] = c1[i] * r3 + r6;
                c1[i] += r6 * r2;
                r6 = r5;
            }
        }
        c2[0] = r3;

        /* compute the solution of the system with
           principal minor of order n + 1
        */
        r5 = 0.0;
        for (i=0; i<n; i++) {
            r5 += a2[i] * x[n1-i];
        }

        r6 = (b[n] - r5) / r1;
        for (i=0; i<n; i++) {
            x[i] += c2[i] * r6;
        }
        x[n] = r6;
    }

 bailout:

    if (err && dt != null) {
        *dt = NADBL;
    }

    free(c1);
    free(c2);

    return err;
}

/*
 * gretl_toeplitz_solve:
 * @c: Toeplitz column.
 * @r: Toeplitz row (null implies that r = c).
 * @b: right-hand side vector.
 * @pd: flag for positive definitness assumed
 * @det: on exit, determinant if pd==0, log determinant otherwise.
 * @err: location to receive error code.
 *
 * Solves Tx = b for the unknown vector x, where T is a Toeplitz
 * matrix, that is (zero-based)
 *
 * T_{ij} = c_{i-j} for i <= j
 * T_{ij} = r_{i-j} for i > j
 *
 * Note that c[0] should equal r[0].
 *
 * Returns: a newly allocated vector, containing the solution, x.
 */

GretlVector * gretl_toeplitz_solve(const(GretlVector) * c,
                                   const(GretlVector) * r,
                                   const(GretlVector) * b,
                                   int pd,
                                   double * det,
                                   int * err)
{
    int m = gretl_vector_get_length(c);
    GretlMatrixRaw * y = null;

    if (det != null) {
        *det = NADBL;
    }

    if (gretl_is_complex(c) || gretl_is_complex(b) ||
        (r != null && gretl_is_complex(r))) {
        fprintf(stderr, "E_CMPLX in gretl_toeplitz_solve\n");
        *err = GretlError.cmplx;
        return null;
    }

    /* a few sanity checks */

    if (m == 0 || m != gretl_vector_get_length(b) ||
        (r != null && m != gretl_vector_get_length(r))) {
        *err = GretlError.nonconf;
        return null;
    }

    if (r != null && r.val[0] != c.val[0]) {
        *err = GretlError.data;
        return null;
    }

    y = gretl_column_vector_alloc(m);

    if (y == null) {
        *err = GretlError.alloc;
    } else {
        /* invoke gretlized netlib routine */
        const(double) * a1 = r != null ? r.val : c.val;
        const(double) * a2 = c.val + 1;

        *err = tsld1(a1, a2, b.val, y.val, det, m, pd);
        if (*err) {
            gretl_matrix_free(y);
            y = null;
        }
    }

    return y;
}

private void gretl_matrix_cum(GretlMatrixRaw * m, int i, int j, double x)
{
    m.val[j * m.rows + i] += x;
}

/* FIXME set this to a positive value under macOS, to take advantage
   of VecLib?
*/

private __gshared int blas_mnk_min = -1;

/*
 * set_blas_mnk_min:
 * @mnk: value to set.
 *
 * By default all matrix multiplication within libgretl is
 * done using native code, but there is the possibility of
 * farming out multiplication to the BLAS.
 *
 * When multiplying an m x n matrix into an n x k matrix
 * libgretl finds the product of the dimensions, m*n*k,
 * and compares this with an internal threshhold variable,
 * blas_mnk_min. If and only if blas_mnk_min >= 0 and
 * n*m*k >= blas_mnk_min, then we use the BLAS.
 */

void set_blas_mnk_min(int mnk)
{
    blas_mnk_min = mnk;
}

/*
 * get_blas_mnk_min:
 *
 * Returns: the value of the internal variable blas_mnk_min.
 * See set_blas_mnk_min().
 */

int get_blas_mnk_min()
{
    return blas_mnk_min;
}

private int use_blas(int m, int n, int k)
{
    if (blas_mnk_min >= 0) {
        guint64 mnk = cast(guint64) m * n * k;

        return mnk >= cast(guint64) blas_mnk_min;
    } else {
        return 0;
    }
}

private void gretl_blas_dsyrk(const(GretlMatrixRaw) * a, int atr,
                              GretlMatrixRaw * c, GretlMatrixMod cmod)
{
    char uplo = 'U';
    char tr = (atr) ? 'T' : 'N';
    integer n = c.rows;
    integer k = (atr) ? a.rows : a.cols;
    integer lda = a.rows;
    double x, alpha = 1.0, beta = 0.0;
    int i, j;

    if (cmod == GretlMatrixMod.cumulate) {
        beta = 1.0;
    } else if (cmod == GretlMatrixMod.decrement) {
        alpha = -1.0;
        beta = 1.0;
    }

    dsyrk_(&uplo, &tr, &n, &k, &alpha, a.val, &lda,
           &beta, c.val, &n);

    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            x = gretl_matrix_get(c, i, j);
            gretl_matrix_set(c, j, i, x);
        }
    }
}

private void gretl_st_result(GretlMatrixRaw * c, int i, int j, double x,
                             GretlMatrixMod m)
{
    if (m == GretlMatrixMod.cumulate) {
        c.val[j * c.rows + i] += x;
        if (i != j) c.val[i * c.rows + j] += x;
    } else if (m == GretlMatrixMod.decrement) {
        c.val[j * c.rows + i] -= x;
        if (i != j) c.val[i * c.rows + j] -= x;
    } else {
        gretl_matrix_set(c, i, j, x);
        gretl_matrix_set(c, j, i, x);
    }
}

private int
matrix_multiply_self_transpose(const(GretlMatrixRaw) * a, int atr,
                               GretlMatrixRaw * c, GretlMatrixMod cmod)
{
    int i, j, k;
    int nc = (atr) ? a.cols : a.rows;
    int nr = (atr) ? a.rows : a.cols;
    int idx1, idx2;
    double x;

    if (c.rows != nc) {
        return GretlError.nonconf;
    }

    if (use_blas(nc, nc, nr)) {
        gretl_blas_dsyrk(a, atr, c, cmod);
        return 0;
    }

    if (c.rows == 1) {
        k = a.cols * a.rows;
        if (cmod != GretlMatrixMod.cumulate) {
            c.val[0] = 0.0;
        }
        for (i=0; i<k; i++) {
            c.val[0] += a.val[i] * a.val[i];
        }
        return 0;
    }

    if (atr) {
        for (i=0; i<nc; i++) {
            for (j=i; j<nc; j++) {
                idx1 = i * a.rows;
                idx2 = j * a.rows;
                x = 0.0;
                for (k=0; k<nr; k++) {
                    x += a.val[idx1++] * a.val[idx2++];
                }
                gretl_st_result(c,i,j,x,cmod);
            }
        }
    } else {
        for (i=0; i<nc; i++) {
            for (j=i; j<nc; j++) {
                idx1 = i;
                idx2 = j;
                x = 0.0;
                for (k=0; k<nr; k++) {
                    x += a.val[idx1] * a.val[idx2];
                    idx1 += a.rows;
                    idx2 += a.rows;
                }
                gretl_st_result(c,i,j,x,cmod);
            }
        }
    }

    return 0;
}

private int
matrix_multiply_self_transpose_single(const(GretlMatrixRaw) * a,
                                      int atr,
                                      GretlMatrixRaw * c,
                                      GretlMatrixMod cmod)
{
    int i, j, k;
    int nc = (atr) ? a.cols : a.rows;
    int nr = (atr) ? a.rows : a.cols;
    int idx1, idx2;
    double x;

    if (c.rows != nc) {
        return GretlError.nonconf;
    }

    if (c.rows == 1) {
        k = a.cols * a.rows;
        if (cmod != GretlMatrixMod.cumulate) {
            c.val[0] = 0.0;
        }
        for (i=0; i<k; i++) {
            c.val[0] += a.val[i] * a.val[i];
        }
        return 0;
    }

    if (atr) {
        for (i=0; i<nc; i++) {
            for (j=i; j<nc; j++) {
                idx1 = i * a.rows;
                idx2 = j * a.rows;
                x = 0.0;
                for (k=0; k<nr; k++) {
                    x += a.val[idx1++] * a.val[idx2++];
                }
                gretl_st_result(c,i,j,x,cmod);
            }
        }
    } else {
        for (i=0; i<nc; i++) {
            for (j=i; j<nc; j++) {
                idx1 = i;
                idx2 = j;
                x = 0.0;
                for (k=0; k<nr; k++) {
                    x += a.val[idx1] * a.val[idx2];
                    idx1 += a.rows;
                    idx2 += a.rows;
                }
                gretl_st_result(c,i,j,x,cmod);
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_XTX_new:
 * @X: matrix to process.
 *
 * Returns: a newly allocated matrix containing X'X, or
 * null on error.
*/

GretlMatrixRaw * gretl_matrix_XTX_new(const(GretlMatrixRaw) * X)
{
    GretlMatrixRaw * XTX = null;

    if (!gretl_is_null_matrix(X)) {
        XTX = gretl_matrix_alloc(X.cols, X.cols);
    }

    if (XTX == null) {
        fprintf(stderr, "gretl_matrix_XTX_new: %d x %d is too big\n",
                X.cols, X.cols);
    } else {
        matrix_multiply_self_transpose(X, 1, XTX, GretlMatrixMod.none);
        maybe_preserve_names(XTX, X, NameFlags.colnames, null);
    }

    return XTX;
}

/*
 * gretl_matrix_packed_XTX_new:
 * @X: matrix to process.
 * @nasty: location to receive warning if any diagonal element
 * is less than double.epsilon.
 *
 * Performs the multiplication X'X producing a packed result,
 * that is, the vech() of the full solution.
 *
 * Returns: a newly allocated column vector containing the
 * unique elements of X'X, or null on error.
*/

private GretlMatrixRaw * gretl_matrix_packed_XTX_new(const(GretlMatrixRaw) * X,
                                                  int * nasty)
{
    GretlMatrixRaw * XTX = null;
    double x;
    int i, j, k, nc, nr, n;

    if (gretl_is_null_matrix(X)) {
        return null;
    }

    nc = X.cols;
    nr = X.rows;
    n = (nc * nc + nc) / 2;
    XTX = gretl_matrix_alloc(n, 1);

    if (XTX == null) {
        return null;
    }

    n = 0;
    for (i=0; i<nc; i++) {
        for (j=i; j<nc; j++) {
            x = 0.0;
            for (k=0; k<nr; k++) {
                x += X.val[i*nr+k] * X.val[j*nr+k];
            }
            if (i == j && x < double.epsilon) {
                *nasty = 1;
            }
            XTX.val[n++] = x;
        }
    }

    return XTX;
}

void gretl_blas_dgemm(const(GretlMatrixRaw) * a, int atr,
                      const(GretlMatrixRaw) * b, int btr,
                      GretlMatrixRaw * c, GretlMatrixMod cmod,
                      int m, int n, int k)
{
    char TransA = atr ? 'T' : 'N';
    char TransB = btr ? 'T' : 'N';
    double alpha = 1.0, beta = 0.0;

    if (cmod == GretlMatrixMod.cumulate) {
        beta = 1.0;
    } else if (cmod == GretlMatrixMod.decrement) {
        alpha = -1.0;
        beta = 1.0;
    }

    dgemm_(&TransA, &TransB, &m, &n, &k,
           &alpha, a.val, &a.rows, b.val, &b.rows, &beta,
           c.val, &c.rows);
}

void gretl_blas_dsymm(const(GretlMatrixRaw) * a, int asecond,
                      const(GretlMatrixRaw) * b, int upper,
                      GretlMatrixRaw * c, GretlMatrixMod cmod,
                      int m, int n)
{
    double alpha = 1.0, beta = 0.0;
    char side = asecond ? 'R' : 'L';
    char uplo = upper ? 'U' : 'L';

    if (cmod == GretlMatrixMod.cumulate) {
        beta = 1.0;
    } else if (cmod == GretlMatrixMod.decrement) {
        alpha = -1.0;
        beta = 1.0;
    }

    dsymm_(&side, &uplo, &m, &n, &alpha, a.val, &a.rows,
           b.val, &b.rows, &beta, c.val, &c.rows);
}

/*
  BLAS dtrmm performs one of the matrix-matrix operations

    B <- alpha * op(A) * B (if side = 'L') or
    B <- alpha * B * op(A) (if side = 'R')

  where alpha is a scalar, B is an m x n matrix, A is an upper
  (if uplo = 'U') or lower (if uplo = 'L') triangular matrix,
  and op(A) is either A (if transa = 'N') or A' (transa = 'T').

  The @flags argument here must contain side, uplo and transa.
*/

void gretl_blas_dtrmm(const(GretlMatrixRaw) * a,
                      GretlMatrixRaw * b,
                      const(char) * flags)
{
    char side = flags[0];
    char uplo = flags[1];
    char transa = flags[2];
    char diag = 'N';
    double alpha = 1.0;
    integer m = b.rows;
    integer n = b.cols;
    integer lda = a.rows;

    dtrmm_(&side, &uplo, &transa, &diag, &m, &n,
           &alpha, cast(double *) a.val, &lda, b.val, &m);
}

/* below: a native C re-write of netlib BLAS dgemm.f: note that
   for gretl's purposes we do not support values of 'beta'
   other than 0 or 1
*/

private void gretl_dgemm(const(GretlMatrixRaw) * a, int atr,
                         const(GretlMatrixRaw) * b, int btr,
                         GretlMatrixRaw * c, GretlMatrixMod cmod,
                         int m, int n, int k)
{
    const(double) * A = a.val;
    const(double) * B = b.val;
    double * C = c.val;
    double x, alpha = 1.0;
    int beta = 0;
    int ar = a.rows;
    int br = b.rows;
    int cr = c.rows;
    int i, j, l;

    if (cmod == GretlMatrixMod.cumulate) {
        beta = 1;
    } else if (cmod == GretlMatrixMod.decrement) {
        alpha = -1.0;
        beta = 1;
    }

    if (!btr) {
        if (!atr) {
            /* C := alpha*A*B + beta*C */
            for (j=0; j<n; j++) {
                if (beta == 0) {
                    for (i=0; i<m; i++) {
                        C[j*cr+i] = 0.0;
                    }
                }
                for (l=0; l<k; l++) {
                    if (B[j*br+l] != 0.0) {
                        x = alpha * B[j*br+l];
                        for (i=0; i<m; i++) {
                            C[j*cr+i] += x * A[l*ar+i];
                        }
                    }
                }
            }
        } else {
            /* C := alpha*A'*B + beta*C */
            for (j=0; j<n; j++) {
                for (i=0; i<m; i++) {
                    x = 0.0;
                    for (l=0; l<k; l++) {
                        x += A[i*ar+l] * B[j*br+l];
                    }
                    if (beta == 0) {
                        C[j*cr+i] = alpha * x;
                    } else {
                        C[j*cr+i] += alpha * x;
                    }
                }
            }
        }
    } else {
        if (!atr) {
            /* C := alpha*A*B' + beta*C */
            for (j=0; j<n; j++) {
                if (beta == 0) {
                    for (i=0; i<m; i++) {
                        C[j*cr+i] = 0.0;
                    }
                }
                for (l=0; l<k; l++) {
                    if (B[l*br+j] != 0.0) {
                        x = alpha * B[l*br+j];
                        for (i=0; i<m; i++) {
                            C[j*cr+i] += x * A[l*ar+i];
                        }
                    }
                }
            }
        } else {
            /* C := alpha*A'*B' + beta*C */
            for (j=0; j<n; j++) {
                for (i=0; i<m; i++) {
                    x = 0.0;
                    for (l=0; l<k; l++) {
                        x += A[i*ar+l] * B[l*br+j];
                    }
                    if (beta == 0) {
                        C[j*cr+i] = alpha * x;
                    } else {
                        C[j*cr+i] += alpha * x;
                    }
                }
            }
        }
    }
}

private int
matmul_mod_w_scalar(double x, const(GretlMatrixRaw) * m, int mtr,
                    GretlMatrixRaw * c, GretlMatrixMod cmod)
{
    int cr = mtr ? m.cols : m.rows;
    int cc = mtr ? m.rows : m.cols;

    if (c.rows != cr || c.cols != cc) {
        return GretlError.nonconf;
    }

    if (mtr) {
        double xm, cij;
        int i, j, k = 0;

        for (i=0; i<cr; i++) {
            for (j=0; j<cc; j++) {
                xm = x * m.val[k++];
                if (cmod == GretlMatrixMod.cumulate) {
                    cij = gretl_matrix_get(c, i, j) + xm;
                } else if (cmod == GretlMatrixMod.decrement) {
                    cij = gretl_matrix_get(c, i, j) - xm;
                } else {
                    cij = xm;
                }
                gretl_matrix_set(c, i, j, cij);
            }
        }
    } else {
        double xm;
        int i, n = cr * cc;

        for (i=0; i<n; i++) {
            xm = x * m.val[i];
            if (cmod == GretlMatrixMod.cumulate) {
                c.val[i] += xm;
            } else if (cmod == GretlMatrixMod.decrement) {
                c.val[i] -= xm;
            } else {
                c.val[i] = xm;
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_multiply_mod:
 * @a: left-hand matrix.
 * @amod: modifier: GretlMatrixMod.none or GretlMatrixMod.transpose.
 * @b: right-hand matrix.
 * @bmod: modifier: GretlMatrixMod.none or GretlMatrixMod.transpose.
 * @c: matrix to hold the product.
 * @cmod: modifier: GretlMatrixMod.none or GretlMatrixMod.cumulate to
 * add the result to the existing value of @c, or
 * GretlMatrixMod.decrement to subtract from the existing value
 * of @c.
 *
 * Multiplies @a (or a-transpose) into @b (or b transpose),
 * with the result written into @c.
 *
 * Returns: 0 on success; non-zero error code on
 * failure.
*/

int gretl_matrix_multiply_mod(const(GretlMatrixRaw) * a, GretlMatrixMod amod,
                              const(GretlMatrixRaw) * b, GretlMatrixMod bmod,
                              GretlMatrixRaw * c, GretlMatrixMod cmod)
{
    const int atr = (amod == GretlMatrixMod.transpose);
    const int btr = (bmod == GretlMatrixMod.transpose);
    int lrows, lcols;
    int rrows, rcols;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        gretl_is_null_matrix(c)) {
        return GretlError.data;
    }

    if (a == c || b == c) {
        fputs("gretl_matrix_multiply:\n product matrix must be "
              ~ "distinct from both input matrices\n", stderr);
        fprintf(stderr, "a = %p, b = %p, c = %p\n",
                cast(void *) a, cast(void *) b, cast(void *) c);
        return 1;
    }

    if (a == b && atr != btr && c.rows == c.cols) {
        return matrix_multiply_self_transpose(a, atr, c, cmod);
    }

    if (a.rows == 1 && a.cols == 1) {
        return matmul_mod_w_scalar(a.val[0], b, btr, c, cmod);
    } else if (b.rows == 1 && b.cols == 1) {
        return matmul_mod_w_scalar(b.val[0], a, atr, c, cmod);
    }

    lrows = (atr) ? a.cols : a.rows;
    lcols = (atr) ? a.rows : a.cols;
    rrows = (btr) ? b.cols : b.rows;
    rcols = (btr) ? b.rows : b.cols;

    if (lcols != rrows) {
        fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
        fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
                lrows, lcols, rrows, rcols, c.rows, c.cols);
        return GretlError.nonconf;
    }

    if (c.rows != lrows || c.cols != rcols) {
        fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
        fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
                lrows, lcols, rrows, rcols, c.rows, c.cols);
        return GretlError.nonconf;
    }

    if (use_blas(lrows, rcols, lcols)) {
        gretl_blas_dgemm(a, atr, b, btr, c, cmod, lrows, rcols, lcols);
    } else {
        gretl_dgemm(a, atr, b, btr, c, cmod, lrows, rcols, lcols);
    }

    return 0;
}

/* single-threaded version of gretl_matrix_multiply_mod()
   for use when we're performing some larger task under
   OMP.
*/

int gretl_matrix_multiply_mod_single(const(GretlMatrixRaw) * a,
                                     GretlMatrixMod amod,
                                     const(GretlMatrixRaw) * b,
                                     GretlMatrixMod bmod,
                                     GretlMatrixRaw * c,
                                     GretlMatrixMod cmod)
{
    const int atr = (amod == GretlMatrixMod.transpose);
    const int btr = (bmod == GretlMatrixMod.transpose);
    int lrows, lcols;
    int rrows, rcols;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        gretl_is_null_matrix(c)) {
        return GretlError.data;
    }

    if (a == c || b == c) {
        fputs("gretl_matrix_multiply:\n product matrix must be "
              ~ "distinct from both input matrices\n", stderr);
        fprintf(stderr, "a = %p, b = %p, c = %p\n",
                cast(void *) a, cast(void *) b, cast(void *) c);
        return 1;
    }

    if (a == b && atr != btr && c.rows == c.cols) {
        return matrix_multiply_self_transpose_single(a, atr, c, cmod);
    }

    if (a.rows == 1 && a.cols == 1) {
        return matmul_mod_w_scalar(a.val[0], b, btr, c, cmod);
    } else if (b.rows == 1 && b.cols == 1) {
        return matmul_mod_w_scalar(b.val[0], a, atr, c, cmod);
    }

    lrows = (atr) ? a.cols : a.rows;
    lcols = (atr) ? a.rows : a.cols;
    rrows = (btr) ? b.cols : b.rows;
    rcols = (btr) ? b.rows : b.cols;

    if (lcols != rrows) {
        fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
        fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
                lrows, lcols, rrows, rcols, c.rows, c.cols);
        return GretlError.nonconf;
    }

    if (c.rows != lrows || c.cols != rcols) {
        fputs("gretl_matrix_multiply_mod: matrices not conformable\n", stderr);
        fprintf(stderr, " Requested (%d x %d) * (%d x %d) = (%d x %d)\n",
                lrows, lcols, rrows, rcols, c.rows, c.cols);
        return GretlError.nonconf;
    }

    gretl_dgemm(a, atr, b, btr, c, cmod, lrows, rcols, lcols);

    return 0;
}

/*
 * gretl_matrix_I_kronecker:
 * @p: dimension of left-hand identity matrix.
 * @B: right-hand matrix, r x s.
 * @K: target matrix, (p * r) x (p * s).
 *
 * Writes the Kronecker product of the identity matrix
 * of order @p and @B into @K.
 *
 * Returns: 0 on success, GretlError.nonconf if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_I_kronecker(int p, const(GretlMatrixRaw) * B,
                         GretlMatrixRaw * K)
{
    double bkl;
    int r, s, d;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(B)) {
        return GretlError.data;
    }

    r = B.rows;
    s = B.cols;

    if (K.rows != p * r || K.cols != p * s) {
        return GretlError.nonconf;
    }

    for (i=0; i<p; i++) {
        ioff = i * r;
        for (j=0; j<p; j++) {
            /* block ij is an r * s matrix, I_{ij} * B */
            d = (i == j);
            joff = j * s;
            for (k=0; k<r; k++) {
                Ki = ioff + k;
                for (l=0; l<s; l++) {
                    Kj = joff + l;
                    if (d) {
                        bkl = gretl_matrix_get(B, k, l);
                        gretl_matrix_set(K, Ki, Kj, bkl);
                    } else {
                        gretl_matrix_set(K, Ki, Kj, 0.0);
                    }
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_I_kronecker_new:
 * @p: dimension of left-hand identity matrix.
 * @B: right-hand matrix, r x s.
 * @err: location to receive error code.
 *
 * Writes the Kronecker product of the identity matrix
 * of order @r and @B into a newly allocated matrix.
 *
 * Returns: the new matrix, or null on failure.
 */

GretlMatrixRaw *
gretl_matrix_I_kronecker_new(int p, const(GretlMatrixRaw) * B, int * err)
{
    GretlMatrixRaw * K;

    if (gretl_is_null_matrix(B)) {
        *err = GretlError.data;
        return null;
    }

    K = gretl_matrix_alloc(p * B.rows, p * B.cols);

    if (K == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_I_kronecker(p, B, K);
    }

    return K;
}

/*
 * gretl_matrix_kronecker_I:
 * @A: left-hand matrix, p x q.
 * @r: dimension of right-hand identity matrix.
 * @K: target matrix, (p * r) x (q * r).
 *
 * Writes the Kronecker product of @A and the identity
 * matrix of order @r into @K.
 *
 * Returns: 0 on success, GretlError.nonconf if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_kronecker_I(const(GretlMatrixRaw) * A, int r,
                         GretlMatrixRaw * K)
{
    double aij;
    int p, q;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(A)) {
        return GretlError.data;
    }

    p = A.rows;
    q = A.cols;

    if (K.rows != p * r || K.cols != q * r) {
        return GretlError.nonconf;
    }

    for (i=0; i<p; i++) {
        ioff = i * r;
        for (j=0; j<q; j++) {
            /* block ij is an r * r matrix, a_{ij} * I_r */
            aij = gretl_matrix_get(A, i, j);
            joff = j * r;
            for (k=0; k<r; k++) {
                Ki = ioff + k;
                for (l=0; l<r; l++) {
                    Kj = joff + l;
                    if (k == l) {
                        gretl_matrix_set(K, Ki, Kj, aij);
                    } else {
                        gretl_matrix_set(K, Ki, Kj, 0.0);
                    }
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_kronecker_I_new:
 * @A: left-hand matrix, p x q.
 * @r: dimension of right-hand identity matrix.
 * @err: location to receive error code.
 *
 * Writes into a newl allocated matrix the Kronecker
 * product of @A and the identity matrix of order @r.
 *
 * Returns: the new matrix, or null on failure.
 */

GretlMatrixRaw *
gretl_matrix_kronecker_I_new(const(GretlMatrixRaw) * A, int r, int * err)
{
    GretlMatrixRaw * K;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return null;
    }

    K = gretl_matrix_alloc(A.rows * r, A.cols * r);

    if (K == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_kronecker_I(A, r, K);
    }

    return K;
}

/*
 * gretl_matrix_kronecker_product:
 * @A: left-hand matrix, p x q.
 * @B: right-hand matrix, r x s.
 * @K: target matrix, (p * r) x (q * s).
 *
 * Writes the Kronecker product of @A and @B into @K.
 *
 * Returns: 0 on success, GretlError.nonconf if matrix @K is
 * not correctly dimensioned for the operation.
 */

int
gretl_matrix_kronecker_product(const(GretlMatrixRaw) * A,
                               const(GretlMatrixRaw) * B,
                               GretlMatrixRaw * K)
{
    double aij, bkl;
    int p, q, r, s;
    int i, j, k, l;
    int ioff, joff;
    int Ki, Kj;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(B) ||
        gretl_is_null_matrix(K)) {
        return GretlError.data;
    }

    p = A.rows;
    q = A.cols;
    r = B.rows;
    s = B.cols;

    if (K.rows != p * r || K.cols != q * s) {
        return GretlError.nonconf;
    }

    for (i=0; i<p; i++) {
        ioff = i * r;
        for (j=0; j<q; j++) {
            /* block ij is an r * s matrix, a_{ij} * B */
            aij = gretl_matrix_get(A, i, j);
            joff = j * s;
            for (k=0; k<r; k++) {
                Ki = ioff + k;
                for (l=0; l<s; l++) {
                    bkl = gretl_matrix_get(B, k, l);
                    Kj = joff + l;
                    gretl_matrix_set(K, Ki, Kj, aij * bkl);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_kronecker_product_new:
 * @A: left-hand matrix, p x q.
 * @B: right-hand matrix, r x s.
 * @err: location to receive error code.
 *
 * Returns: A newly allocated (p * r) x (q * s) matrix which
 * is the Kronecker product of matrices @A and @B, or null
 * on failure.
 */

GretlMatrixRaw *
gretl_matrix_kronecker_product_new(const(GretlMatrixRaw) * A,
                                   const(GretlMatrixRaw) * B,
                                   int * err)
{
    GretlMatrixRaw * K;
    guint64 chk;
    int p, q, r, s;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
        *err = GretlError.data;
        return null;
    }

    p = A.rows;
    q = A.cols;
    r = B.rows;
    s = B.cols;

    chk = cast(guint64) p * r * q * s;
    if (chk > int.max) {
        *err = GretlError.invarg;
        gretl_errmsg_set(_("Too many elements").ptr);
        return null;
    }

    K = gretl_matrix_alloc(p * r, q * s);

    if (K == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_kronecker_product(A, B, K);
    }

    return K;
}

/*
 * gretl_matrix_hdproduct:
 * @A: left-hand matrix, r x p.
 * @B: right-hand matrix, r x q or null.
 * @C: target matrix, r x (p * q).
 *
 * Writes into @C the horizontal direct product of @A and @B.
 * That is, $C_i' = A_i' \otimes B_i'$ (in TeX notation). If @B
 * is null, then it's understood to be equal to @A.
 *
 * Returns: 0 on success, GretlError.nonconf if @A and @B have different
 * numbers of rows or matrix @C is not correctly dimensioned for the
 * operation.
 */

int gretl_matrix_hdproduct(const(GretlMatrixRaw) * A,
                           const(GretlMatrixRaw) * B,
                           GretlMatrixRaw * C)
{
    double aij, bik;
    int r, p, q;
    int i, j, k;
    int ndx, retcols;
    int do_symmetric;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(C)) {
        return GretlError.data;
    }

    r = A.rows;
    p = A.cols;
    do_symmetric = gretl_is_null_matrix(B);

    if (do_symmetric) {
        q = p;
        retcols = p * (p+1) / 2;
        if (C.rows != r || C.cols != retcols) {
            return GretlError.nonconf;
        }
    } else {
        q = B.cols;
        retcols = p * q;
        if (B.rows != r || C.rows != r || C.cols != retcols) {
            return GretlError.nonconf;
        }
    }

    for (i=0; i<r; i++) {
        ndx = 0;
        for (j=0; j<p; j++) {
            aij = gretl_matrix_get(A, i, j);
            if (do_symmetric) {
                for (k=j; k<q; k++) {
                    bik = gretl_matrix_get(A, i, k);
                    gretl_matrix_set(C, i, ndx++, aij*bik);
                }
            } else if (aij != 0.0) {
                ndx = j * q;
                for (k=0; k<q; k++) {
                    bik = gretl_matrix_get(B, i, k);
                    gretl_matrix_set(C, i, ndx + k, aij*bik);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_hdproduct_new:
 * @A: left-hand matrix, r x p.
 * @B: right-hand matrix, r x q or null.
 * @err: location to receive error code.
 *
 * If @B is null, then it is implicitly taken as equal to @A; in this case,
 * the returned matrix only contains the non-redundant elements; therefore,
 * it has ncols = p*(p+1)/2 elements. Otherwise, all the products are computed
 * and ncols = p*q.
 *
 * Returns: newly allocated r x ncols matrix which is the horizontal
 * direct product of matrices @A and @B, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_hdproduct_new(const(GretlMatrixRaw) * A,
                                         const(GretlMatrixRaw) * B,
                                         int * err)
{
    GretlMatrixRaw * K = null;
    int r, p, q, ncols;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
    } else if (gretl_is_complex(A) && gretl_is_null_matrix(B)) {
        *err = GretlError.data;
    } else if (gretl_is_complex(A) || gretl_is_complex(B)) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_hdproduct_new\n");
        *err = GretlError.cmplx;
    }

    if (*err) {
        return null;
    }

    r = A.rows;
    p = A.cols;

    if (gretl_is_null_matrix(B)) {
        q = A.cols;
        ncols = p * (p+1) / 2;
    } else {
        if (B.rows != r) {
            *err = GretlError.nonconf;
        } else {
            q = B.cols;
            ncols = p * q;
        }
    }

    if (!*err) {
        K = gretl_zero_matrix_new(r, ncols);
        if (K == null) {
            *err = GretlError.alloc;
        } else {
            gretl_matrix_hdproduct(A, B, K);
        }
    }

    return K;
}

/* returns the sequence of bits in the binary expansion of @s */

private char * binary_expansion(int s, int * t, int * pow2)
{
    char * bits = null;
    double l2 = log_2(s);
    int k = cast(int) floor(l2);

    if (l2 == k) {
        *pow2 = 1;
    }

    *t = k;
    bits = cast(char *) calloc(k + 1, 1);

    if (bits != null) {
        while (1) {
            bits[k] = 1;
            s = cast(int) (s - pow(2.0, k));
            if (s == 0) {
                break;
            }
            l2 = log_2(s);
            k = cast(int) floor(l2);
        }
    }

    return bits;
}

private GretlMatrixRaw * matrix_frac_pow(const(GretlMatrixRaw) * m,
                                      double a, int * err)
{
    GretlMatrixRaw * ret;
    GretlMatrixRaw * tmp;
    GretlMatrixRaw * lam;
    double eps = 1.0e-12;
    int n = m.rows;

    tmp = gretl_matrix_copy(m);
    ret = gretl_matrix_alloc(n, n);

    if (tmp == null || ret == null) {
        gretl_matrix_free(tmp);
        gretl_matrix_free(ret);
        *err = GretlError.alloc;
        return null;
    }

    lam = gretl_symmetric_matrix_eigenvals(tmp, 1, err);

    if (!*err) {
        if (lam.val[0] < -eps) {
            /* be a little lenient with positive
               semidefinite matrices */
            *err = GretlError.notpd;
        } else if (lam.val[0] < eps && a < 0) {
            /* but don't allow negative exponents if @m
               is singular */
            *err = GretlError.invarg;
        } else {
            double x, y, a2 = a/2;
            int i, j;

            for (j=0; j<n; j++) {
                y = pow(fabs(lam.val[j]), a2);
                for (i=0; i<n; i++) {
                    x = gretl_matrix_get(tmp, i, j);
                    gretl_matrix_set(tmp, i, j, x * y);
                }
            }
            matrix_multiply_self_transpose(tmp, 0, ret, GretlMatrixMod.none);
        }
    }

    gretl_matrix_free(lam);
    gretl_matrix_free(tmp);

    if (*err) {
        gretl_matrix_free(ret);
        ret = null;
    }

    return ret;
}

private GretlMatrixRaw * matrix_int_pow(const(GretlMatrixRaw) * A,
                                     int s, int * err)
{
    GretlMatrixRaw * B = null;
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * W = null;
    char * bits = null;
    int n, t, pow2 = 0;

    n = A.rows;

    if (s < 2) {
        B = (s == 0) ? gretl_identity_matrix_new(n) :
            gretl_matrix_copy(A);
        if (B == null) {
            *err = GretlError.alloc;
        }
        return B;
    }

    bits = binary_expansion(s, &t, &pow2);
    if (bits == null) {
        *err = GretlError.alloc;
        return null;
    }

    B = gretl_matrix_copy_tmp(A);
    C = gretl_matrix_alloc(n, n);

    if (!pow2) {
        W = gretl_matrix_alloc(n, n);
    }

    if (B == null || C == null || (W == null && !pow2)) {
        gretl_matrix_free(C);
        C = null;
        *err = GretlError.alloc;
    }

    if (!*err) {
        int q = 0;

        while (bits[q] == 0) {
            /* B = B^2 */
            gretl_matrix_multiply(B, B, C);
            gretl_matrix_copy_values(B, C);
            q++;
        }

        if (!pow2) {
            /* more work needed */
            int k;

            gretl_matrix_copy_values(C, B);

            for (k=q+1; k<=t; k++) {
                /* B = B^2 */
                gretl_matrix_multiply(B, B, W);
                gretl_matrix_copy_values(B, W);
                if (bits[k]) {
                    /* C = CB */
                    gretl_matrix_multiply(C, B, W);
                    gretl_matrix_copy_values(C, W);
                }
            }
        }
    }

    gretl_matrix_free(B);
    gretl_matrix_free(W);
    free(bits);

    return C;
}

/*
 * gretl_matrix_pow:
 * @A: square source matrix.
 * @s: exponent.
 * @err: location to receive error code.
 *
 * Calculates the matrix A^s. If @s is a non-negative integer
 * Golub and Van Loan's Algorithm 11.2.2 ("Binary Powering")
 * is used. Otherwise @A must be positive definite, and the
 * power is computed via the eigen-decomposition of @A.
 *
 * Returns: allocated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_pow(const(GretlMatrixRaw) * A,
                               double s, int * err)
{
    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
    } else if (A.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_pow\n");
        *err = GretlError.cmplx;
    } else if (A.rows != A.cols) {
        *err = GretlError.nonconf;
    }

    if (*err) {
        return null;
    } else if (s != floor(s) || s < 0) {
        return matrix_frac_pow(A, s, err);
    } else {
        int k = gretl_int_from_double(s, err);

        if (*err) {
            return null;
        } else {
            return matrix_int_pow(A, k, err);
        }
    }
}

/*
 * gretl_vector_dot_product:
 * @a: first vector.
 * @b: second vector.
 * @err: pointer to receive error code (zero on success,
 * non-zero on failure), or null.
 *
 * Returns: The dot (scalar) product of @a and @b, or NADBL on
 * failure.
 */

double gretl_vector_dot_product(const(GretlVector) * a,
                                const(GretlVector) * b,
                                int * err)
{
    int i, dima, dimb;
    double dp = 0.0;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        return NADBL;
    } else if (a.is_complex || b.is_complex) {
        *err = GretlError.cmplx;
        return NADBL;
    }

    dima = (a.rows > 1) ? a.rows : a.cols;
    dimb = (b.rows > 1) ? b.rows : b.cols;

    if (!gretl_is_vector(a) || !gretl_is_vector(b) || dima != dimb) {
        if (err != null) {
            *err = GretlError.nonconf;
        }
        dp = NADBL;
    } else {
        for (i=0; i<dima; i++) {
            dp += a.val[i] * b.val[i];
        }
    }

    return dp;
}

/*
 * gretl_matrix_dot_product:
 * @a: left-hand matrix.
 * @amod: modifier: GretlMatrixMod.none or GretlMatrixMod.transpose.
 * @b: right-hand matrix.
 * @bmod: modifier: GretlMatrixMod.none or GretlMatrixMod.transpose.
 * @err: pointer to receive error code (zero on success,
 * non-zero on failure), or null.
 *
 * Returns: The dot (scalar) product of @a (or @a-transpose) and
 * @b (or @b-transpose), or NADBL on failure.
 */

double gretl_matrix_dot_product(const(GretlMatrixRaw) * a,
                                GretlMatrixMod amod,
                                const(GretlMatrixRaw) * b,
                                GretlMatrixMod bmod,
                                int * err)
{
    GretlMatrixRaw * c = null;
    double ret = NADBL;
    int myerr = 0;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        return NADBL;
    }

    if (gretl_is_vector(a) && gretl_is_vector(b)) {
        return gretl_vector_dot_product(a, b, err);
    }

    c = gretl_matrix_alloc(1, 1);
    if (c == null) {
        myerr = GretlError.alloc;
    }

    if (!myerr) {
        myerr = gretl_matrix_multiply_mod(a, amod, b, bmod,
                                          c, GretlMatrixMod.none);
    }

    if (!myerr) {
        ret = c.val[0];
    }

    gretl_matrix_free(c);

    if (err != null) {
        *err = myerr;
    }

    return ret;
}

/*
 * dot_operator_conf:
 * @A: first matrix.
 * @B: second matrix.
 * @r: pointer to rows of result.
 * @c: pointer to columns of result.
 *
 * Used to establish the validity of a "dot operation" such as A .* B
 * or A .+ B and, if the operation is valid, the dimensions of the
 * result.
 *
 * Returns: a numeric code identifying the convention to be used,
 * where ConfType.none indicates non-conformability.
 */

ConfType dot_operator_conf(const(GretlMatrixRaw) * A,
                           const(GretlMatrixRaw) * B,
                           int * r, int * c)
{
    int ra = A.rows;
    int rb = B.rows;
    int ca = A.cols;
    int cb = B.cols;
    int confr = (ra == rb);
    int confc = (ca == cb);
    int colva = (ca == 1);
    int colvb = (cb == 1);
    int rowva = (ra == 1);
    int rowvb = (rb == 1);
    ConfType ret = ConfType.none;

    if (confr && confc) {
        /* element-by-element operation */
        ret = ConfType.elements;
        *r = ra;
        *c = ca;
    } else if (confr && colva) {
        /* rows match; A is a column vector */
        ret = ConfType.a_colvec;
        *r = ra;
        *c = (colva) ? cb : ca;
    } else if (confr && colvb) {
        /* rows match; B is a column vector */
        ret = ConfType.b_colvec;
        *r = ra;
        *c = (colva) ? cb : ca;
    } else if (confc && rowva) {
        /* columns match; A is a row vector */
        ret = ConfType.a_rowvec;
        *r = (rowva) ? rb : ra;
        *c = ca;
    } else if (confc && rowvb) {
        /* columns match; B is a row vector */
        ret = ConfType.b_rowvec;
        *r = (rowva) ? rb : ra;
        *c = ca;
    } else if (ra == 1 && ca == 1) {
        /* A is a scalar in disguise */
        ret = ConfType.a_scalar;
        *r = rb;
        *c = cb;
    } else if (rb == 1 && cb == 1) {
        /* B is a scalar in disguise */
        ret = ConfType.b_scalar;
        *r = ra;
        *c = ca;
    } else if (colva && rowvb) {
        /* A is a column and B is a row */
        ret = ConfType.ac_br;
        *r = ra;
        *c = cb;
    } else if (rowva && colvb) {
        /* A is a row and B is a column */
        ret = ConfType.ar_bc;
        *r = rb;
        *c = ca;
    }

    return ret;
}

/* give an OPENMP parallelization? */

private void vec_x_op_vec_y(double * z, const(double) * x,
                            const(double) * y, int n,
                            int op)
{
    int i;

    switch (op) {
    case '*':
        for (i=0; i<n; i++) {
            z[i] = x[i] * y[i];
        }
        break;
    case '/':
        for (i=0; i<n; i++) {
            z[i] = x[i] / y[i];
        }
        break;
    case '+':
        for (i=0; i<n; i++) {
            z[i] = x[i] + y[i];
        }
        break;
    case '-':
        for (i=0; i<n; i++) {
            z[i] = x[i] - y[i];
        }
        break;
    case '^':
        for (i=0; i<n; i++) {
            z[i] = pow(x[i], y[i]);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = x[i] == y[i];
        }
        break;
    case '>':
        for (i=0; i<n; i++) {
            z[i] = x[i] > y[i];
        }
        break;
    case '<':
        for (i=0; i<n; i++) {
            z[i] = x[i] < y[i];
        }
        break;
    case ']':
        for (i=0; i<n; i++) {
            z[i] = x[i] >= y[i];
        }
        break;
    case '[':
        for (i=0; i<n; i++) {
            z[i] = x[i] <= y[i];
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = x[i] != y[i];
        }
        break;
    default:
        break;
    }
}

private void vec_x_op_y(double * z, const(double) * x,
                        double y, int n, int op)
{
    int i;

    switch (op) {
    case '*':
        for (i=0; i<n; i++) {
            z[i] = x[i] * y;
        }
        break;
    case '/':
        for (i=0; i<n; i++) {
            z[i] = x[i] / y;
        }
        break;
    case '+':
        for (i=0; i<n; i++) {
            z[i] = x[i] + y;
        }
        break;
    case '-':
        for (i=0; i<n; i++) {
            z[i] = x[i] - y;
        }
        break;
    case '^':
        for (i=0; i<n; i++) {
            z[i] = pow(x[i], y);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = x[i] == y;
        }
        break;
    case '>':
        for (i=0; i<n; i++) {
            z[i] = x[i] > y;
        }
        break;
    case '<':
        for (i=0; i<n; i++) {
            z[i] = x[i] < y;
        }
        break;
    case ']':
        for (i=0; i<n; i++) {
            z[i] = x[i] >= y;
        }
        break;
    case '[':
        for (i=0; i<n; i++) {
            z[i] = x[i] <= y;
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = x[i] != y;
        }
        break;
    default:
        break;
    }
}

private void x_op_vec_y(double * z, double x,
                        const(double) * y, int n,
                        int op)
{
    int i;

    switch (op) {
    case '*':
        for (i=0; i<n; i++) {
            z[i] = x * y[i];
        }
        break;
    case '/':
        for (i=0; i<n; i++) {
            z[i] = x / y[i];
        }
        break;
    case '+':
        for (i=0; i<n; i++) {
            z[i] = x + y[i];
        }
        break;
    case '-':
        for (i=0; i<n; i++) {
            z[i] = x - y[i];
        }
        break;
    case '^':
        for (i=0; i<n; i++) {
            z[i] = pow(x, y[i]);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = x == y[i];
        }
        break;
    case '>':
        for (i=0; i<n; i++) {
            z[i] = x > y[i];
        }
        break;
    case '<':
        for (i=0; i<n; i++) {
            z[i] = x < y[i];
        }
        break;
    case ']':
        for (i=0; i<n; i++) {
            z[i] = x >= y[i];
        }
        break;
    case '[':
        for (i=0; i<n; i++) {
            z[i] = x <= y[i];
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = x != y[i];
        }
        break;
    default:
        break;
    }
}

private double x_op_y(double x, double y, int op)
{
    switch (op) {
    case '*':
        return x * y;
    case '/':
        return x / y;
    case '+':
        return x + y;
    case '-':
        return x - y;
    case '^':
        return pow(x, y);
    case '=':
        return x == y;
    case '>':
        return x > y;
    case '<':
        return x < y;
    case ']':
        return x >= y;
    case '[':
        return x <= y;
    case '!':
        return x != y;
    default:
        return 0;
    }
}

/*
 * gretl_matrix_dot_op:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @op: operator.
 * @err: location to receive error code.
 *
 * @op should be one of the standard ASCII representations of
 * arithmetical operators:
 *
 * '*' : multiply, '/' divide
 * '+' : add, '-' subtract
 * '^': exponentiate
 *
 * or one of these test symbols:
 *
 * '=' : is equal to
 * '>' : is greater than, '<' is less than
 * ']' : is greater than or equal to
 * '[' : is less than or equal to
 * '!' : is not equal to
 *
 * Returns: a new matrix, each of whose elements is the result
 * of (x @op y), where x and y are the  corresponding elements of
 * the matrices @a and @b (or null on failure).
 */

GretlMatrixRaw * gretl_matrix_dot_op(const(GretlMatrixRaw) * a,
                                  const(GretlMatrixRaw) * b,
                                  int op, int * err)
{
    ConfType conftype;
    GretlMatrixRaw * c = null;
    double x, y;
    int nr, nc;
    int i, j, off;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        *err = GretlError.data;
        return null;
    }

    conftype = dot_operator_conf(a, b, &nr, &nc);

    if (conftype == ConfType.none) {
        fputs("gretl_matrix_dot_op: matrices not conformable\n", stderr);
        fprintf(stderr, " op = '%c', A is %d x %d, B is %d x %d\n",
                cast(char) op, a.rows, a.cols, b.rows, b.cols);
        *err = GretlError.nonconf;
        return null;
    }

    c = gretl_matrix_alloc(nr, nc);
    if (c == null) {
        *err = GretlError.alloc;
        return null;
    }

    math_err_init();

    switch (conftype) {
    case ConfType.elements:
        vec_x_op_vec_y(c.val, a.val, b.val, nr*nc, op);
        break;
    case ConfType.a_colvec:
        for (i=0; i<nr; i++) {
            x = a.val[i];
            for (j=0; j<nc; j++) {
                y = gretl_matrix_get(b, i, j);
                y = x_op_y(x, y, op);
                gretl_matrix_set(c, i, j, y);
            }
        }
        break;
    case ConfType.b_colvec:
        for (i=0; i<nr; i++) {
            y = b.val[i];
            for (j=0; j<nc; j++) {
                x = gretl_matrix_get(a, i, j);
                x = x_op_y(x, y, op);
                gretl_matrix_set(c, i, j, x);
            }
        }
        break;
    case ConfType.a_rowvec:
        off = 0;
        for (j=0; j<nc; j++) {
            x = a.val[j];
            x_op_vec_y(c.val + off, x, b.val + off, nr, op);
            off += nr;
        }
        break;
    case ConfType.b_rowvec:
        off = 0;
        for (j=0; j<nc; j++) {
            y = b.val[j];
            vec_x_op_y(c.val + off, a.val + off, y, nr, op);
            off += nr;
        }
        break;
    case ConfType.a_scalar:
        x_op_vec_y(c.val, a.val[0], b.val, nr*nc, op);
        break;
    case ConfType.b_scalar:
        vec_x_op_y(c.val, a.val, b.val[0], nr*nc, op);
        break;
    case ConfType.ac_br:
        for (i=0; i<nr; i++) {
            x = a.val[i];
            for (j=0; j<nc; j++) {
                y = b.val[j];
                y = x_op_y(x, y, op);
                gretl_matrix_set(c, i, j, y);
            }
        }
        break;
    case ConfType.ar_bc:
        for (j=0; j<nc; j++) {
            x = a.val[j];
            for (i=0; i<nr; i++) {
                y = b.val[i];
                y = x_op_y(x, y, op);
                gretl_matrix_set(c, i, j, y);
            }
        }
        break;
    default: /* hush a warning */
        break;
    }

    if (errno) {
        *err = math_err_check("gretl_matrix_dot_op", errno);
        if (*err) {
            gretl_matrix_free(c);
            c = null;
        }
    }

    return c;
}

/* Multiplication or division for complex matrices in the old
   gretl representation, with real parts in the first column
   and imaginary parts (if present) in the second.
*/

private GretlMatrixRaw *
gretl_matrix_complex_muldiv(const(GretlMatrixRaw) * a,
                            const(GretlMatrixRaw) * b,
                            int multiply,
                            int force_complex,
                            int * err)
{
    GretlMatrixRaw * c = null;
    const(double) * ar, ai;
    const(double) * br, bi;
    double * cr, ci;
    Cmplx az, bz, cz;
    int m, n, p, q;
    int i, izero = 1;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        *err = GretlError.data;
        return null;
    }

    m = a.rows;
    n = a.cols;
    p = b.rows;
    q = b.cols;

    if (m != p) {
        *err = GretlError.nonconf;
        return null;
    }

    if ((n != 1 && n != 2) || (q != 1 && q != 2)) {
        *err = GretlError.nonconf;
        return null;
    }

    if (force_complex) {
        p = 2;
    } else {
        p = (n == 1 && q == 1) ? 1 : 2;
    }

    c = gretl_matrix_alloc(m, p);
    if (c == null) {
        *err = GretlError.alloc;
        return null;
    }

    math_err_init();

    ar = a.val;
    ai = (a.cols == 2) ? ar + m : null;

    br = b.val;
    bi = (b.cols == 2) ? br + m : null;

    cr = c.val;
    ci = (c.cols == 2) ? cr + m : null;

    for (i=0; i<m; i++) {
        az = ai == null ? Cmplx(ar[i], 0.0) : Cmplx(ar[i], ai[i]);
        bz = bi == null ? Cmplx(br[i], 0.0) : Cmplx(br[i], bi[i]);
        if (multiply) {
            cz = az * bz;
        } else {
            cz = az / bz;
        }
        cr[i] = cz.r;
        if (ci != null) {
            ci[i] = cz.i;
            if (ci[i] != 0.0) {
                izero = 0;
            }
        }
    }

    if (errno) {
        *err = math_err_check("gretl_matrix_complex_muldiv", errno);
        if (*err) {
            gretl_matrix_free(c);
            c = null;
        }
    }

    if (!*err && !force_complex && c.cols == 2 && izero) {
        /* drop the all-zero imaginary part */
        *err = gretl_matrix_realloc(c, c.rows, 1);
        if (*err) {
            gretl_matrix_free(c);
            c = null;
        }
    }

    return c;
}

/*
 * gretl_matrix_complex_multiply:
 * @a: m x (1 or 2) matrix.
 * @b: m x (1 or 2) matrix.
 * @force_complex: see below.
 * @err: location to receive error code.
 *
 * Computes the complex product of @a and @b.  The first
 * column in these matrices is assumed to contain real
 * values, and the second column (if present) imaginary
 * coefficients.
 *
 * Returns: a matrix with the result of the multiplication
 * of the two vectors of complex numbers. If both @a and @b have no
 * imaginary part and the @force_complex flag is zero, the return
 * value will be m x 1, otherwise it will be m x 2.
 */

GretlMatrixRaw * gretl_matrix_complex_multiply(const(GretlMatrixRaw) * a,
                                            const(GretlMatrixRaw) * b,
                                            int force_complex,
                                            int * err)
{
    return gretl_matrix_complex_muldiv(a, b, 1, force_complex, err);
}

/*
 * gretl_matrix_complex_divide:
 * @a: m x (1 or 2) matrix.
 * @b: m x (1 or 2) matrix.
 * @force_complex: see below.
 * @err: location to receive error code.
 *
 * Computes the complex division of @a over @b.  The first
 * column in these matrices is assumed to contain real
 * values, and the second column (if present) imaginary
 * coefficients.
 *
 * Returns: a matrix with the result of the division of the
 * two vectors of complex numbers. If both @a and @b have no
 * imaginary part and the @force_complex flag is zero, the return
 * value will be m x 1, otherwise it will be m x 2.
 */

GretlMatrixRaw * gretl_matrix_complex_divide(const(GretlMatrixRaw) * a,
                                          const(GretlMatrixRaw) * b,
                                          int force_complex,
                                          int * err)
{
    return gretl_matrix_complex_muldiv(a, b, 0, force_complex, err);
}

/*
 * gretl_rmatrix_vector_stat:
 * @m: source matrix.
 * @vs: the required statistic or quantity: sum, product or mean.
 * @rowwise: if non-zero go by rows, otherwise go by columns.
 * @skip_na: ignore missing values.
 * @err: location to receive error code.
 *
 * Returns: a row vector or column vector containing the sums,
 * products or means of the columns or rows of @m. See also
 * gretl_cmatrix_vector_stat() for the complex variant.
 */

GretlMatrixRaw * gretl_rmatrix_vector_stat(const(GretlMatrixRaw) * m,
                                        GretlVecStat vs,
                                        int rowwise,
                                        int skip_na,
                                        int * err)
{
    GretlMatrixRaw * ret;
    double x, mij;
    int r, c, i, j;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
        return null;
    }

    r = rowwise ? m.rows : 1;
    c = rowwise ? 1 : m.cols;

    ret = gretl_matrix_alloc(r, c);
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (rowwise) {
        /* by rows */
        double x0 = vs == GretlVecStat.prod ? 1 : 0;
        int ncols;

        for (i=0; i<m.rows; i++) {
            ncols = m.cols;
            x = x0;
            for (j=0; j<m.cols; j++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    ncols--;
                    continue;
                } else if (vs == GretlVecStat.prod) {
                    x *= mij;
                } else {
                    x += mij;
                }
            }
            if (ncols == 0) {
                x = NADBL;
            } else if (vs == GretlVecStat.mean) {
                x /= ncols;
            }
            gretl_matrix_set(ret, i, 0, x);
        }
    } else {
        /* by columns */
        double x0 = vs == GretlVecStat.prod ? 1 : 0;
        int nrows;

        for (j=0; j<m.cols; j++) {
            nrows = m.rows;
            x = x0;
            for (i=0; i<m.rows; i++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    nrows--;
                    continue;
                } else if (vs == GretlVecStat.prod) {
                    x *= mij;
                } else {
                    x += mij;
                }
            }
            if (nrows == 0) {
                x = NADBL;
            } else if (vs == GretlVecStat.mean) {
                x /= nrows;
            }
            gretl_matrix_set(ret, 0, j, x);
        }
    }

    if (rowwise) {
        maybe_preserve_names(ret, m, NameFlags.rownames, null);
    } else {
        maybe_preserve_names(ret, m, NameFlags.colnames, null);
    }

    return ret;
}

/*
 * gretl_matrix_column_sd:
 * @m: source matrix.
 * @df: degrees of freedom for standard deviations.
 * @skip_na: ignore missing values.
 * @err: location to receive error code.
 *
 * Returns: a row vector containing the standard deviations of
 * the columns of @m, or null on failure. If @df is positive
 * it is used as the divisor when calculating the column
 * variance, otherwise the divisor is the number of rows in
 * @m.
 */

GretlMatrixRaw * gretl_matrix_column_sd(const(GretlMatrixRaw) * m,
                                     int df, int skip_na,
                                     int * err)
{
    GretlMatrixRaw * s;
    double mij, xbar, dev, v;
    int i, j, nrows;
    int den, dfc = 0;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
        return null;
    } else if (m.is_complex) {
        *err = GretlError.cmplx;
        return null;
    }

    s = gretl_matrix_alloc(1, m.cols);

    if (s == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (df > 0) {
        /* respect the df supplied by the caller, taking it
           as implicitly defining the df correction, dfc
        */
        dfc = m.rows - df;
    }

    for (j=0; j<m.cols; j++) {
        nrows = m.rows;
        xbar = 0.0;
        for (i=0; i<m.rows; i++) {
            mij = gretl_matrix_get(m, i, j);
            if (na(mij) && skip_na) {
                nrows--;
            } else {
                xbar += mij;
            }
        }
        den = nrows - dfc;
        if (nrows < 2 || den <= 0 || na(xbar)) {
            s.val[j] = NADBL;
        } else {
            xbar /= nrows;
            v = 0.0;
            for (i=0; i<m.rows; i++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    ; /* skip */
                } else {
                    dev = mij - xbar;
                    v += dev * dev;
                }
            }
            s.val[j] = sqrt(v / den);
        }
    }

    return s;
}

/*
 * gretl_matrix_demean_by_row:
 * @m: matrix on which to operate.
 *
 * For each row of @m, subtracts the row mean from each
 * element on the row.
 */

void gretl_matrix_demean_by_row(GretlMatrixRaw * m)
{
    double x, rmean;
    int i, j;

    for (i=0; i<m.rows; i++) {
        rmean = 0;
        for (j=0; j<m.cols; j++) {
            rmean += gretl_matrix_get(m, i, j);
        }
        rmean /= m.cols;
        for (j=0; j<m.cols; j++) {
            x = gretl_matrix_get(m, i, j);
            gretl_matrix_set(m, i, j, x - rmean);
        }
    }
}

/*
 * gretl_matrix_center:
 * @m: matrix on which to operate.
 *
 * Subtracts the column mean from each column of @m.
 */

int gretl_matrix_center(GretlMatrixRaw * m, int skip_na)
{
    double mij, xbar;
    int i, j, nrows;

    for (j=0; j<m.cols; j++) {
        nrows = m.rows;
        xbar = 0;
        for (i=0; i<m.rows; i++) {
            mij = gretl_matrix_get(m, i, j);
            if (na(mij) && skip_na) {
                nrows--;
            } else {
                xbar += mij;
            }
        }
        if (nrows > 0) {
            xbar /= nrows;
            for (i=0; i<m.rows; i++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    ; /* skip */
                } else {
                    gretl_matrix_set(m, i, j, mij - xbar);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_standardize:
 * @m: matrix on which to operate.
 * @dfcorr: degrees of freedom correction.
 * @skip_na: ignore missing values.
 *
 * Subtracts the column mean from each column of @m and
 * divides by the column standard deviation, using @dfcorr
 * as degrees of freedom correction (0 for MLE).
 */

int gretl_matrix_standardize(GretlMatrixRaw * m, int dfcorr, int skip_na)
{
    double mij, x, xbar, sdc, den;
    int i, j, nrows;

    if (m.rows < 2) {
        return GretlError.toofew;
    }

    for (j=0; j<m.cols; j++) {
        nrows = m.rows;
        xbar = sdc = 0;
        for (i=0; i<m.rows; i++) {
            mij = gretl_matrix_get(m, i, j);
            if (na(mij) && skip_na) {
                nrows--;
            } else {
                xbar += mij;
            }
        }
        if (nrows == 0) {
            xbar = NADBL;
        } else {
            xbar /= nrows;
            for (i=0; i<m.rows; i++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    ; /* skip */
                } else {
                    x = mij - xbar;
                    gretl_matrix_set(m, i, j, x);
                    sdc += x * x;
                }
            }
        }
        den = nrows - dfcorr;
        if (nrows < 2 || den <= 0 || na(xbar)) {
            for (i=0; i<m.rows; i++) {
                gretl_matrix_set(m, i, j, NADBL);
            }
        } else {
            sdc = sqrt(sdc / den);
            for (i=0; i<m.rows; i++) {
                mij = gretl_matrix_get(m, i, j);
                if (na(mij) && skip_na) {
                    ; /* skip */
                } else {
                    gretl_matrix_set(m, i, j, mij / sdc);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_quantiles:
 * @m: matrix on which to operate.
 * @p: vector of desired quantiles.
 * @err: location to receive error code.
 *
 * Returns: a matrix containing the @p quantiles
 * of the columns of @m, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_quantiles(const(GretlMatrixRaw) * m,
                                     const(GretlMatrixRaw) * p,
                                     int * err)
{
    GretlMatrixRaw * qvals;
    const(double) * mval;
    double * a, q;
    int i, j, k;
    int n, plen;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.invarg;
        return null;
    }

    plen = gretl_vector_get_length(p);

    if (plen == 0) {
        *err = GretlError.invarg;
        return null;
    }

    for (i=0; i<plen; i++) {
        if (p.val[i] <= 0 || p.val[i] >= 1 || na(p.val[i])) {
            *err = GretlError.invarg;
            return null;
        }
    }

    qvals = gretl_matrix_alloc(plen, m.cols);
    if (qvals == null) {
        *err = GretlError.alloc;
        return null;
    }

    n = m.rows;
    a = cast(double *) malloc(n * (*a).sizeof);
    q = cast(double *) malloc(plen * (*q).sizeof);

    if (a == null || q == null) {
        *err = GretlError.alloc;
        gretl_matrix_free(qvals);
        free(a);
        free(q);
        return null;
    }

    mval = m.val;

    for (j=0; j<m.cols && !*err; j++) {
        k = 0;
        for (i=0; i<n; i++) {
            if (!na(mval[i])) {
                a[k++] = mval[i];
            }
        }
        memcpy(q, p.val, plen * (*q).sizeof);
        if (k == 0) {
            for (i=0; i<plen; i++) {
                gretl_matrix_set(qvals, i, j, NADBL);
            }
        } else {
            *err = gretl_array_quantiles(a, k, q, plen);
            if (!*err) {
                for (i=0; i<plen; i++) {
                    gretl_matrix_set(qvals, i, j, q[i]);
                }
            }
        }
        mval += n;
    }

    if (*err) {
        gretl_matrix_free(qvals);
        qvals = null;
    }

    free(a);
    free(q);

    return qvals;
}

/*
 * gretl_matrix_multiply:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @c: matrix to hold the product.
 *
 * Multiplies @a into @b, with the result written into @c.
 *
 * Returns: 0 on success; non-zero error code on
 * failure.
 */

int gretl_matrix_multiply(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b,
                          GretlMatrixRaw * c)
{
    int err = 0;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        gretl_is_null_matrix(c)) {
        return GretlError.data;
    }

    if (matrix_is_scalar(a)) {
        err = gretl_matrix_copy_values(c, b);
        if (!err) {
            gretl_matrix_multiply_by_scalar(c, a.val[0]);
        }
    } else if (matrix_is_scalar(b)) {
        err = gretl_matrix_copy_values(c, a);
        if (!err) {
            gretl_matrix_multiply_by_scalar(c, b.val[0]);
        }
    } else {
        err = gretl_matrix_multiply_mod(a, GretlMatrixMod.none,
                                        b, GretlMatrixMod.none,
                                        c, GretlMatrixMod.none);
    }

    return err;
}

int gretl_matrix_multiply_single(const(GretlMatrixRaw) * a,
                                 const(GretlMatrixRaw) * b,
                                 GretlMatrixRaw * c)
{
    int err = 0;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b) ||
        gretl_is_null_matrix(c)) {
        return GretlError.data;
    }

    if (matrix_is_scalar(a)) {
        err = gretl_matrix_copy_values(c, b);
        if (!err) {
            gretl_matrix_multiply_by_scalar(c, a.val[0]);
        }
    } else if (matrix_is_scalar(b)) {
        err = gretl_matrix_copy_values(c, a);
        if (!err) {
            gretl_matrix_multiply_by_scalar(c, b.val[0]);
        }
    } else {
        err = gretl_matrix_multiply_mod_single(a, GretlMatrixMod.none,
                                               b, GretlMatrixMod.none,
                                               c, GretlMatrixMod.none);
    }

    return err;
}

/*
 * gretl_matrix_multiply_new:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @err: location for error code.
 *
 * Multiplies @a into @b, with the result written into a newly
 * allocated matrix.
 *
 * Returns: matrix product on success, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_multiply_new(const(GretlMatrixRaw) * a,
                                        const(GretlMatrixRaw) * b,
                                        int * err)
{
    GretlMatrixRaw * c;

    if (gretl_is_null_matrix(a) || gretl_is_null_matrix(b)) {
        *err = GretlError.data;
        return null;
    }

    if (a.cols != b.rows) {
        fprintf(stderr, "gretl_matrix_multiply_new: requested (%d x %d) * (%d x %d)\n",
                a.rows, a.cols, b.rows, b.cols);
        *err = GretlError.nonconf;
        return null;
    }

    c = gretl_matrix_alloc(a.rows, b.cols);
    if (c == null) {
        *err = GretlError.alloc;
        return null;
    }

    *err = gretl_matrix_multiply_mod(a, GretlMatrixMod.none,
                                     b, GretlMatrixMod.none,
                                     c, GretlMatrixMod.none);

    if (*err) {
        gretl_matrix_free(c);
        c = null;
    }

    return c;
}

private int matrix_divide_by_scalmat(GretlMatrixRaw * num,
                                     const(GretlMatrixRaw) * den)
{
    int i, n = num.rows * num.cols;

    if (num.is_complex) {
        Cmplx zden;

        zden = den.is_complex ? den.z[0] : Cmplx(den.val[0], 0.0);
        for (i=0; i<n; i++) {
            num.z[i] /= zden;
        }
    } else {
        if (den.is_complex) {
            return GretlError.types;
        } else {
            for (i=0; i<n; i++) {
                num.val[i] /= den.val[0];
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_divide:
 * @a: left-hand matrix.
 * @b: right-hand matrix.
 * @mod: GretlMatrixMod.none for left division, or
 * GretlMatrixMod.transpose for right division.
 * @err: location to receive error code.
 *
 * Follows the semantics of Matlab/Octave for left and right
 * matrix "division". In left division, A \ B is in principle
 * equivalent to A^{-1} * B, and in right division A / B is
 * in principle equivalent to A * B^{-1}, but the result is
 * obtained without explicit computation of the inverse.
 *
 * In left division @a and @b must have the same number of
 * rows; in right division they must have the same number
 * of columns.
 *
 * Returns: the "quotient" matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_divide(const(GretlMatrixRaw) * a,
                                  const(GretlMatrixRaw) * b,
                                  GretlMatrixMod mod,
                                  int * err)
{
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * AT = null;
    GretlMatrixRaw * BT = null;
    GretlMatrixRaw * Tmp;

    if (gretl_is_null_matrix(a) ||
        gretl_is_null_matrix(b)) {
        *err = GretlError.data;
        return null;
    }

    /* detect and handle scalar cases */
    if (mod == GretlMatrixMod.none && is_one_by_one(a)) {
        Q = gretl_matrix_copy(b);
        if (Q == null) {
            *err = GretlError.alloc;
        } else {
            *err = matrix_divide_by_scalmat(Q, a);
        }
        return Q;
    } else if (mod == GretlMatrixMod.transpose && is_one_by_one(b)) {
        Q = gretl_matrix_copy(a);
        if (Q == null) {
            *err = GretlError.alloc;
        } else {
            *err = matrix_divide_by_scalmat(Q, b);
        }
        return Q;
    }

    if (mod == GretlMatrixMod.none && a.rows != b.rows) {
        *err = GretlError.nonconf;
    } else if (mod == GretlMatrixMod.transpose && a.cols != b.cols) {
        *err = GretlError.nonconf;
    }

    if (*err) {
        return Q;
    }

    /* USE_TRISOLVE is 1 in the C source */
    /* 2024-09-06 */
    if (mod == GretlMatrixMod.none) {
        int tri = matrix_is_triangular(a);

        if (tri) {
            Q = gretl_matrix_copy(b);
            if (Q == null) {
                *err = GretlError.alloc;
            } else {
                gretl_triangular_solve(a, Q, mod, cast(GretlMatrixStructure) tri);
            }
            return Q;
        }
    } else {
        /* the following may not be worth doing? */
        int tri = matrix_is_triangular(b);

        if (tri) {
            Q = gretl_matrix_copy_transpose(a);
            if (Q == null) {
                *err = GretlError.alloc;
            } else {
                gretl_triangular_solve(b, Q, mod, cast(GretlMatrixStructure) tri);
                gretl_matrix_transpose_in_place(Q);
            }
            return Q;
        }
    }

    if (mod == GretlMatrixMod.transpose) {
        AT = gretl_matrix_copy_transpose(b);
        BT = gretl_matrix_copy_transpose(a);
        if (AT == null || BT == null) {
            *err = GretlError.alloc;
            goto bailout;
        } else {
            a = AT;
            b = BT;
        }
    }

    Q = gretl_matrix_copy(b);
    if (Q == null) {
        *err = GretlError.alloc;
    } else {
        Tmp = gretl_matrix_copy(a);
        if (Tmp == null) {
            *err = GretlError.alloc;
        } else {
            *err = gretl_matrix_solve(Tmp, Q);
            gretl_matrix_free(Tmp);
        }
    }

    if (mod == GretlMatrixMod.transpose && *err == 0) {
        Tmp = Q;
        Q = gretl_matrix_copy_transpose(Tmp);
        if (Q == null) {
            *err = GretlError.alloc;
        }
        gretl_matrix_free(Tmp);
    }

 bailout:

    if (mod == GretlMatrixMod.transpose) {
        gretl_matrix_free(AT);
        gretl_matrix_free(BT);
    }

    if (*err && Q != null) {
        gretl_matrix_free(Q);
        Q = null;
    }

    return Q;
}

/*
 * gretl_general_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 *
 * Estimates the reciprocal condition number of the general
 * real matrix @m (in the 1-norm), using the LAPACK
 * functions dgetrf() and dgecon().
 *
 * Returns: the estimate, or NADBL on failure to allocate memory.
 */

private double gretl_general_matrix_rcond(const(GretlMatrixRaw) * A,
                                          int * err)
{
    GretlMatrixRaw * a = null;
    char norm = '1';
    integer m, n, lda, info;
    integer * iwork = null;
    integer * ipiv = null;
    double * work = null;
    double rcond = NADBL;

    *err = 0;

    if (gretl_is_null_matrix(A)) {
        return NADBL;
    }

    m = A.rows;
    n = A.cols;
    lda = A.rows;

    a = gretl_matrix_copy_tmp(A);
    work = cast(double *) malloc((4 * n) * (*work).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);
    ipiv = cast(integer *) malloc(min(m, n) * (*ipiv).sizeof);

    if (a == null || work == null || iwork == null || ipiv == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    dgetrf_(&m, &n, a.val, &lda, ipiv, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_general_matrix_rcond:\n"
                ~ " dgetrf failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
        gretl_matrix_print(A, "A in rcond");
        *err = GretlError.data;
        rcond = NADBL;
    } else {
        pivot_check(ipiv, min(m, n));
    }

    if (!*err) {
        double anorm = gretl_matrix_one_norm(A);

        dgecon_(&norm, &n, a.val, &lda, &anorm, &rcond, work, iwork, &info);
        if (info != 0) {
            *err = GretlError.data;
            rcond = NADBL;
        }
    }

 bailout:

    free(work);
    free(iwork);
    free(ipiv);
    gretl_matrix_free(a);

    return rcond;
}

/*
 * gretl_symmetric_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 *
 * Estimates the reciprocal condition number of the real symmetric
 * positive definite matrix @m (in the 1-norm), using the LAPACK
 * functions dpotrf() and dpocon().
 *
 * Returns: the estimate, or NADBL on failure to allocate memory.
 */

double gretl_symmetric_matrix_rcond(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * a = null;
    char uplo = 'L';
    integer n, lda;
    integer info;
    integer * iwork = null;
    double * work = null;
    double rcond = NADBL;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
        return NADBL;
    }

    n = m.rows;
    lda = m.rows;

    a = gretl_matrix_copy_tmp(m);
    work = cast(double *) malloc((3 * n) * (*work).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);

    if (a == null || work == null || iwork == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    dpotrf_(&uplo, &n, a.val, &n, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_symmetric_matrix_rcond: "
                ~ "dpotrf failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
        rcond = 0.0;
    } else {
        double anorm = gretl_matrix_one_norm(m);

        dpocon_(&uplo, &n, a.val, &lda, &anorm, &rcond, work, iwork, &info);
        if (info != 0) {
            *err = 1;
            rcond = NADBL;
        }
    }

 bailout:

    free(work);
    free(iwork);
    gretl_matrix_free(a);

    return rcond;
}

/*
 * gretl_matrix_rcond:
 * @m: matrix to examine.
 * @err: location to receive error code.
 *
 * Estimates the reciprocal condition number of the real
 * matrix @m (in the 1-norm).
 *
 * Returns: the estimate, or NADBL on failure to allocate memory.
 */

double gretl_matrix_rcond(const(GretlMatrixRaw) * m, int * err)
{
    return gretl_general_matrix_rcond(m, err);
}

/*
 * gretl_matrix_cond_index:
 * @m: matrix to examine.
 * @err: location to receive error code.
 *
 * Estimates the condition number (a la Belsley) of the real
 * matrix @m.
 *
 * Returns: the estimate, or NADBL on failure.
 */

double gretl_matrix_cond_index(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * X, XX, v;
    double xij, den, cidx = NADBL;
    int i, j, r, c;

    if (gretl_is_null_matrix(m)) {
        return NADBL;
    }

    r = m.rows;
    c = m.cols;

    X = gretl_matrix_alloc(r, c);
    XX = gretl_matrix_alloc(c, c);

    if (X == null || XX == null) {
        gretl_matrix_free(X);
        gretl_matrix_free(XX);
        *err = GretlError.alloc;
        return NADBL;
    }

    /* normalize columns of @m into X */
    for (j=0; j<c; j++) {
        den = 0.0;
        for (i=0; i<r; i++) {
            xij = gretl_matrix_get(m, i, j);
            den += xij * xij;
        }
        den = sqrt(den);
        for (i=0; i<r; i++) {
            xij = gretl_matrix_get(m, i, j);
            gretl_matrix_set(X, i, j, xij / den);
        }
    }

    /* form X'X */
    gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                              X, GretlMatrixMod.none,
                              XX, GretlMatrixMod.none);

    v = gretl_symmetric_matrix_eigenvals(XX, 0, err);

    if (!*err) {
        cidx = sqrt(v.val[c-1] / v.val[0]);
    }

    gretl_matrix_free(X);
    gretl_matrix_free(XX);
    gretl_matrix_free(v);

    return cidx;
}

/*
 * gretl_matrix_cholesky_decomp:
 * @a: matrix to operate on.
 *
 * Computes the Cholesky factorization of the symmetric,
 * positive definite matrix @a.  On exit the lower triangle of
 * @a is replaced by the factor L, as in a = LL', and the
 * upper triangle is set to zero.  Uses the lapack function
 * dpotrf.
 *
 * Returns: 0 on success; 1 on failure.
 */

int gretl_matrix_cholesky_decomp(GretlMatrixRaw * a)
{
    char uplo = 'L';
    integer n;
    integer info;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    n = a.rows;
    if (a.cols != n) {
        return GretlError.nonconf;
    }

    dpotrf_(&uplo, &n, a.val, &n, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_matrix_cholesky_decomp: info = %d\n",
                cast(int) info);
        err = (info > 0) ? GretlError.notpd : GretlError.data;
    } else {
        gretl_matrix_zero_upper(a);
    }

    return err;
}

/*
 * cholesky_factor_of_inverse:
 * @a: matrix to operate on.
 *
 * Computes the lower-triangular Cholesky factorization of
 * the inverse of a symmetric, positive definite matrix @A.
 * On exit the lower triangle of @A is replaced by the factor
 * L, as in LL' = A^{-1}, and the upper triangle is set to zero.
 * Uses the lapack functions dpotrf and dpotri.
 *
 * Returns: 0 on success; non-zero on failure.
 */

int cholesky_factor_of_inverse(GretlMatrixRaw * a)
{
    int save_nt = 0;
    integer n, info;
    char uplo = 'L';
    int err = 0;

    n = a.cols;
    if (n == 1) {
        a.val[0] = 1.0 / a.val[0];
        return 0;
    }

    if (blas_is_openblas()) {
        maybe_force_single(n, POTRF_MT_MIN, &save_nt);
    }

    /* obtain Cholesky factor of @a */
    dpotrf_(&uplo, &n, a.val, &n, &info);

    if (info == 0) {
        /* obtain inverse of @a */
        dpotri_(&uplo, &n, a.val, &n, &info);
    }

    if (info == 0) {
        /* obtain Cholesky factor of inverse */
        dpotrf_(&uplo, &n, a.val, &n, &info);
    }

    if (info == 0) {
        /* zero the upper triangle of @a */
        int i, j;

        for (j=1; j<n; j++) {
            for (i=0; i<j; i++) {
                gretl_matrix_set(a, i, j, 0.0);
            }
        }
    } else {
        err = (info > 0) ? GretlError.notpd : GretlError.data;
    }

    if (save_nt > 1) {
        blas_set_num_threads(save_nt);
    }

    return err;
}

private int process_psd_root(GretlMatrixRaw * L,
                             const(GretlMatrixRaw) * A,
                             integer rank,
                             integer * piv)
{
    GretlMatrixRaw * LL = null;
    double toler = 1.0e-8;
    int i, j, n = L.rows;
    int err = 0;

    LL = gretl_matrix_alloc(n, n);

    if (LL == null) {
        err = GretlError.alloc;
    } else {
        /* form LL' and compare with A to see if @L is
           really a viable factor
        */
        double dj, dmax = 0;

        gretl_matrix_multiply_mod(L, GretlMatrixMod.none,
                                  L, GretlMatrixMod.transpose,
                                  LL, GretlMatrixMod.none);
        for (j=0; j<n; j++) {
            dj = 0.0;
            for (i=0; i<n; i++) {
                dj += fabs(gretl_matrix_get(LL, i, j) - gretl_matrix_get(A, i, j));
            }
            if (dj > dmax) {
                dmax = dj;
            }
        }

        if (dmax > toler) {
            gretl_errmsg_sprintf("psdroot: norm-test of %g exceeds tolerance (%g)",
                                 dmax, toler);
            err = GretlError.data;
        }

        gretl_matrix_free(LL);
    }

    return err;
}

/* PSD cholesky-type factor via the simple algorithm from
   Golub and Van Loan.
*/

private int real_psd_root(GretlMatrixRaw * a, const(GretlMatrixRaw) * a0)
{
    double d, x1, x2, x3;
    int i, j, k, n = a.rows;
    int err = 0;

    /* Golub and Van Loan, algorithm 4.2.11 */

    for (k=0; k<n && !err; k++) {
        d = gretl_matrix_get(a, k, k);
        if (d > 0) {
            d = sqrt(d);
            gretl_matrix_set(a, k, k, d);
            for (i=k+1; i<n; i++) {
                x1 = gretl_matrix_get(a, i, k);
                gretl_matrix_set(a, i, k, x1 / d);
            }
            for (j=k+1; j<n; j++) {
                x1 = gretl_matrix_get(a, j, k);
                for (i=j; i<n; i++) {
                    x2 = gretl_matrix_get(a, i, j);
                    x3 = gretl_matrix_get(a, i, k);
                    gretl_matrix_set(a, i, j, x2 - x3 * x1);
                }
            }
        } else {
            if (a0 == null && d < -1.0e-8) {
                /* Since we can't perform the check against a0, we'll
                   reject a matrix that has a "significantly" negative
                   diagonal element.
                */
                err = GretlError.data;
            }
            for (i=k; i<n; i++) {
                gretl_matrix_set(a, i, k, 0.0);
            }
        }
    }

    gretl_matrix_zero_triangle(a, 'U');

    if (!err && a0 != null) {
        err = process_psd_root(a, a0, 0, null);
    }

    return err;
}

/*
 * gretl_matrix_psd_root:
 * @a: matrix to operate on.
 * @check: if non-zero, perform a test for psd status.
 *
 * Computes the LL' factorization of the symmetric,
 * positive semidefinite matrix @a.  On successful exit
 * the lower triangle of @a is replaced by the factor L
 * and the upper triangle is set to zero.
 *
 * Returns: 0 on success; non-zero on failure.
 */

int gretl_matrix_psd_root(GretlMatrixRaw * a, int check)
{
    GretlMatrixRaw * a0 = null;
    int err = 0;

    if (gretl_is_null_matrix(a) || a.rows != a.cols) {
        return GretlError.invarg;
    }

    if (check) {
        /* make a copy of @a so we can test for its
           supposed psd attribute
        */
        a0 = gretl_matrix_copy(a);
        if (a0 == null) {
            return GretlError.alloc;
        }
    }

    err = real_psd_root(a, a0);
    gretl_matrix_free(a0);

    return err;
}

int gretl_matrix_QR_pivot_decomp(GretlMatrixRaw * M,
                                 GretlMatrixRaw * R,
                                 GretlMatrixRaw * P)
{
    integer m = M.rows;
    integer n = M.cols;
    integer info = 0;
    integer lwork = -1;
    integer lda = m;
    integer * iwork = null;
    double * tau = null;
    double * work = null;
    integer * jpvt = null;
    int i, j;
    int err = 0;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns").ptr);
        return GretlError.nonconf;
    }

    if (R != null && (R.rows != n || R.cols != n)) {
        return GretlError.nonconf;
    }
    if (P != null && (P.rows != 1 || P.cols != n)) {
        return GretlError.nonconf;
    }

    /* dim of tau is min (m, n) */
    tau = cast(double *) malloc(n * (*tau).sizeof);
    work = cast(double *) lapack_malloc((*work).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);

    if (tau == null || work == null || iwork == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* pivot array */
    jpvt = cast(integer *) malloc(n * (*jpvt).sizeof);
    if (jpvt == null) {
        err = GretlError.alloc;
        goto bailout;
    } else {
        /* pin the first column in place */
        jpvt[0] = 1;
        /* but allow permutation of the others */
        for (i=1; i<n; i++) {
            jpvt[i] = 0;
        }
    }

    /* workspace size query */
    dgeqp3_(&m, &n, M.val, &lda, jpvt, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dgeqp3: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

    /* optimally sized work array */
    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, cast(size_t) lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* run actual QR factorization */
    dgeqp3_(&m, &n, M.val, &lda, jpvt, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dgeqp3: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

    if (R != null) {
        /* copy the upper triangular R out of M */
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    gretl_matrix_set(R, i, j,
                                     gretl_matrix_get(M, i, j));
                } else {
                    gretl_matrix_set(R, i, j, 0.0);
                }
            }
        }
    }

    /* obtain the real "Q" matrix (in M) */
    dorgqr_(&m, &n, &n, M.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dorgqr: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

 bailout:

    free(tau);
    lapack_free(work);
    free(iwork);

    if (P != null) {
        for (i=0; i<n; i++) {
            P.val[i] = jpvt[i];
        }
    }

    free(jpvt);

    return err;
}

/*
 * gretl_matrix_QR_decomp:
 * @M: m x n matrix to be decomposed.
 * @R: n x n matrix into which to write R, as in M = Q * R,
 * or null if this is not wanted.
 *
 * Computes the QR factorization of @M.  On successful exit
 * the matrix @M holds Q, and, if @R is not null, the upper
 * triangle of @R holds R.  Uses the LAPACK functions
 * dgeqrf() and dorgqr().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_matrix_QR_decomp(GretlMatrixRaw * M, GretlMatrixRaw * R)
{
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    double * tau = null;
    double * work = null;
    int i, j;
    int err = 0;

    if (gretl_is_null_matrix(M)) {
        return GretlError.data;
    }

    lda = m = M.rows;
    n = M.cols;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns").ptr);
        return GretlError.nonconf;
    }

    if (R != null && (R.rows != n || R.cols != n)) {
        return GretlError.nonconf;
    }

    /* dim of tau is min (m, n) */
    tau = cast(double *) malloc(n * (*tau).sizeof);
    work = cast(double *) lapack_malloc((*work).sizeof);

    if (tau == null || work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* workspace size query */
    dgeqrf_(&m, &n, M.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dgeqrf: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

    /* optimally sized work array */
    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, cast(size_t) lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* run actual QR factorization */
    dgeqrf_(&m, &n, M.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dgeqrf: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

    if (R != null) {
        /* copy the upper triangular R out of M */
        double x;

        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    x = gretl_matrix_get(M, i, j);
                    gretl_matrix_set(R, i, j, x);
                } else {
                    gretl_matrix_set(R, i, j, 0.0);
                }
            }
        }
    }

    /* obtain the real "Q" matrix (in M) */
    dorgqr_(&m, &n, &n, M.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dorgqr: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    }

 bailout:

    free(tau);
    lapack_free(work);

    return err;
}

private int get_R_rank(const(GretlMatrixRaw) * R)
{
    double d;
    int i, rank = R.rows;

    for (i=0; i<R.rows; i++) {
        d = gretl_matrix_get(R, i, i);
        if (fabs(d) < R_DIAG_MIN || na(d)) {
            rank--;
        }
    }

    return rank;
}

/*
 * gretl_triangular_matrix_rcond:
 * @A: triangular matrix.
 * @uplo: 'U' for upper triangular input, 'L' for lower.
 * @diag: 'N' for non-unit triangular, 'U' for unit triangular.
 *
 * Returns: the reciprocal condition number of @R in the 1-norm,
 * or NADBL on failure.
 */

double gretl_triangular_matrix_rcond(const(GretlMatrixRaw) * A,
                                     char uplo, char diag)
{
    integer * iwork = null;
    double * work = null;
    integer n, info = 0;
    double rcond = NADBL;
    char norm = '1';

    n = A.rows;
    work = cast(double *) lapack_malloc(3 * n * (*work).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);

    if (work != null && iwork != null) {
        dtrcon_(&norm, &uplo, &diag, &n, cast(double *) A.val, &n, &rcond, work,
                iwork, &info);
        if (info != 0) {
            fprintf(stderr, "dtrcon: info = %d\n", cast(int) info);
            rcond = NADBL;
        }
    }

    lapack_free(work);
    free(iwork);

    return rcond;
}

/*
 * gretl_check_QR_rank:
 * @R: matrix R from QR decomposition.
 * @err: location to receive error code.
 * @rcnd: location to receive reciprocal condition number.
 *
 * Checks the reciprocal condition number of R and calculates
 * the rank of the matrix QR.  If @rcnd is not null it receives
 * the reciprocal condition number.
 *
 * Returns: on success, the rank of QR.
 */

int gretl_check_QR_rank(const(GretlMatrixRaw) * R, int * err, double * rcnd)
{
    integer * iwork = null;
    double * work = null;
    integer n, info = 0;
    double rcond;
    char uplo = 'U';
    char diag = 'N';
    char norm = '1';
    int rank;

    if (gretl_is_null_matrix(R)) {
        *err = GretlError.data;
        return 0;
    }

    *err = 0;

    rank = n = R.rows;
    work = cast(double *) lapack_malloc(3 * n * (*work).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);

    if (work == null || iwork == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    dtrcon_(&norm, &uplo, &diag, &n, cast(double *) R.val, &n, &rcond, work,
            iwork, &info);

    if (info != 0) {
        fprintf(stderr, "dtrcon: info = %d\n", cast(int) info);
        *err = 1;
        goto bailout;
    }

    if (rcond < QR_RCOND_MIN) {
        rank = get_R_rank(R);
    }

    if (rcnd != null) {
        *rcnd = rcond;
    }

 bailout:

    lapack_free(work);
    free(iwork);

    return rank;
}

private double svd_smin(const(GretlMatrixRaw) * a, double smax, double eps)
{

    int dmax = (a.rows > a.cols) ? a.rows : a.cols;
    double actual_eps = (na(eps) || eps <= 0) ? 2.20e-16 : eps;

    /* numpy and Matlab use the "Numerical Recipes" convention
       by which eps should be machine epsilon (for 8-byte reals,
       2.20e-16).
    */

    return dmax * actual_eps * smax;
}

/*
 * gretl_matrix_rank:
 * @a: matrix to examine.
 * @eps: value below which a singular value is taken to
 * be zero. Give 0 or NADBL for automatic use of machine
 * epsilon.
 * @err: location to receive error code on failure.
 *
 * Computes the rank of @a via its SV decomposition.
 *
 * Returns: the rank of @a, or 0 on failure.
 */

int gretl_matrix_rank(const(GretlMatrixRaw) * a, double eps,
                      int * err)
{
    GretlMatrixRaw * s = null;
    int i, k, rank = 0;

    if (gretl_is_null_matrix(a)) {
        return 0;
    }

    k = (a.rows < a.cols) ? a.rows : a.cols;

    if (a.rows > 4 * k || a.cols > 4 * k) {
        GretlMatrixRaw * b = gretl_matrix_alloc(k, k);
        GretlMatrixMod mod1, mod2;

        mod1 = a.rows > k ? GretlMatrixMod.transpose : GretlMatrixMod.none;
        mod2 = a.cols > k ? GretlMatrixMod.transpose : GretlMatrixMod.none;
        gretl_matrix_multiply_mod(a, mod1, a, mod2, b, GretlMatrixMod.none);
        *err = real_gretl_matrix_SVD(b, null, &s, null, 0);
        gretl_matrix_free(b);
    } else {
        *err = real_gretl_matrix_SVD(a, null, &s, null, 0);
    }

    if (!*err) {
        double smin = svd_smin(a, s.val[0], eps);

        for (i=0; i<k; i++) {
            if (s.val[i] > smin) {
                rank++;
            }
        }
    }

    gretl_matrix_free(s);

    return rank;
}

/*
 * gretl_invert_triangular_matrix:
 * @a: triangular matrix to invert.
 * @uplo: 'L' for lower triangular @a, 'U' for upper.
 *
 * Computes the inverse of a triangular matrix.  On exit
 * @a is overwritten with the inverse.  Uses the lapack
 * function dtrtri.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_triangular_matrix(GretlMatrixRaw * a, char uplo)
{
    char diag = 'N';
    integer n, info = 0;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    n = a.rows;

    if (a.rows != a.cols) {
        return GretlError.nonconf;
    }

    dtrtri_(&uplo, &diag, &n, a.val, &n, &info);

    if (info < 0) {
        fprintf(stderr, "dtrtri: info = %d\n", cast(int) info);
        err = GretlError.data;
    } else if (info > 0) {
        err = GretlError.singular;
    }

    return err;
}

/*
 * gretl_invert_general_matrix:
 * @a: matrix to invert.
 *
 * Computes the inverse of a general matrix using LU
 * factorization.  On exit @a is overwritten with the inverse.
 * Uses the LAPACK functions dgetrf() and dgetri().
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_general_matrix(GretlMatrixRaw * a)
{
    int save_nt = 0;
    integer n;
    integer info;
    integer lwork;
    integer * ipiv;
    double * work;
    int err = 0;

    if (gretl_is_null_matrix(a) || (a.rows != a.cols)) {
        return GretlError.data;
    }

    n = a.rows;

    ipiv = cast(integer *) malloc(n * (*ipiv).sizeof);
    if (ipiv == null) {
        return GretlError.alloc;
    }

    work = cast(double *) lapack_malloc((*work).sizeof);
    if (work == null) {
        free(ipiv);
        return GretlError.alloc;
    }

    if (blas_is_openblas()) {
        maybe_force_single(n, GETRF_MT_MIN, &save_nt);
    }

    dgetrf_(&n, &n, a.val, &n, ipiv, &info);

    if (info != 0) {
        fprintf(stderr, "dgetrf: matrix is singular (info=%d)\n", info);
        err = GretlError.singular;
        goto bailout;
    } else {
        pivot_check(ipiv, n);
    }

    lwork = -1;
    dgetri_(&n, a.val, &n, ipiv, work, &lwork, &info);

    if (info != 0 || work[0] <= 0.0) {
        err = wspace_fail(info, work[0]);
        goto bailout;
    }

    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    dgetri_(&n, a.val, &n, ipiv, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dgetri: matrix is singular (info = %d)\n", info);
        err = GretlError.singular;
    }

 bailout:

    if (save_nt > 1) {
        blas_set_num_threads(save_nt);
    }

    lapack_free(work);
    free(ipiv);

    return err;
}

/*
 * gretl_matrix_mirror:
 * @m: matrix to expand.
 * @uplo: 'L' or 'U'.
 *
 * If @uplo = 'L', copy the lower triangle of @m into
 * the upper triangle; or if @uplo = 'U' copy the upper
 * triangle into the lower, in either case producing a
 * symmetric result.
 *
 * Returns: 0 on success; non-zero error code if @m is
 * not square.
 */

int gretl_matrix_mirror(GretlMatrixRaw * m, char uplo)
{
    int i, j, n;
    double x;

    if (m.cols != m.rows) {
        fputs("gretl_matrix_mirror: input is not square\n",
              stderr);
        return 1;
    }

    n = m.rows;

    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            if (uplo == 'U') {
                x = gretl_matrix_get(m, i, j);
                gretl_matrix_set(m, j, i, x);
            } else {
                x = gretl_matrix_get(m, j, i);
                gretl_matrix_set(m, i, j, x);
            }
        }
    }

    return 0;
}

/*
 * gretl_invert_diagonal_matrix:
 * @a: matrix to invert.
 *
 * Computes the inverse of a diagonal matrix.
 * On exit @a is overwritten with the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_diagonal_matrix(GretlMatrixRaw * a)
{
    double x;
    int i;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    if (a.cols != a.rows) {
        fputs("gretl_invert_diagonal_matrix: input is not square\n",
              stderr);
        return GretlError.nonconf;
    }

    for (i=0; i<a.rows; i++) {
        if (gretl_matrix_get(a, i, i) == 0.0) {
            return GretlError.singular;
        }
    }

    for (i=0; i<a.rows; i++) {
        x = gretl_matrix_get(a, i, i);
        gretl_matrix_set(a, i, i, 1.0 / x);
    }

    return 0;
}

/*
 * gretl_invert_matrix:
 * @a: matrix to invert.
 *
 * Computes the inverse of matrix @a: on exit @a is
 * overwritten with the inverse.  If @a is diagonal, triangular
 * or symmetric, appropriate simple inversion routines
 * are called.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_matrix(GretlMatrixRaw * a)
{
    int s, err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    s = gretl_matrix_get_structure(a);

    if (s == GretlMatrixStructure.identity) {
        return 0;
    } else if (s == GretlMatrixStructure.diagonal) {
        err = gretl_invert_diagonal_matrix(a);
    } else if (s == GretlMatrixStructure.symmetric) {
        err = real_invert_symmetric_matrix(a, 0, 1, null);
        if (err) {
            err = gretl_invert_symmetric_indef_matrix(a);
        }
    } else if (s == GretlMatrixStructure.lower_triangular) {
        err = gretl_invert_triangular_matrix(a, 'L');
    } else if (s == GretlMatrixStructure.upper_triangular) {
        err = gretl_invert_triangular_matrix(a, 'U');
    } else if (s >= GretlMatrixStructure.square) {
        err = gretl_invert_general_matrix(a);
    } else {
        err = GretlError.nonconf;
    }

    return err;
}

private enum RS_RCOND_MIN = 1.0e-15;

/*
 * gretl_invert_symmetric_indef_matrix:
 * @a: matrix to invert.
 *
 * Computes the inverse of a real symmetric matrix via the
 * Bunch-Kaufman diagonal pivoting method.  Uses the LAPACK
 * functions dsytrf() and dsytri().  On exit @a is overwritten
 * with the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_symmetric_indef_matrix(GretlMatrixRaw * a)
{
    char uplo = 'U';
    integer n, info;
    integer * ipiv;
    integer * iwork;
    integer lwork = -1;
    double anorm, rcond;
    double * work;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    if (a.cols != a.rows) {
        fputs("gretl_invert_symmetric_indef_matrix: input is not square\n",
              stderr);
        return GretlError.nonconf;
    }

    n = a.rows;
    ipiv = cast(integer *) malloc(n * (*ipiv).sizeof);
    iwork = cast(integer *) malloc(n * (*iwork).sizeof);
    work = cast(double *) lapack_malloc((*work).sizeof);

    if (ipiv == null || iwork == null || work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    anorm = gretl_matrix_one_norm(a);

    /* workspace query */
    dsytrf_(&uplo, &n, a.val, &n, ipiv, work, &lwork, &info);
    if (info != 0 || work[0] <= 0.0) {
        err = wspace_fail(info, work[0]);
        goto bailout;
    }

    lwork = cast(integer) work[0];

    if (lwork < 2 * n) {
        lwork = 2 * n;
    }

    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* decompose */
    dsytrf_(&uplo, &n, a.val, &n, ipiv, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dsytrf: matrix is singular\n");
        err = GretlError.singular;
        goto bailout;
    }

    /* check condition number */
    dsycon_(&uplo, &n, a.val, &n, ipiv, &anorm, &rcond,
            work, iwork, &info);
    if (info != 0) {
        fprintf(stderr, "dsycon: info = %d\n", cast(int) info);
        err = 1;
        goto bailout;
    } else if (rcond < RS_RCOND_MIN) {
        fprintf(stderr, "dsycon: rcond = %g\n", rcond);
        err = GretlError.singular;
        goto bailout;
    }

    /* invert */
    dsytri_(&uplo, &n, a.val, &n, ipiv, work, &info);

 bailout:

    lapack_free(work);
    free(ipiv);
    free(iwork);

    if (!err) {
        if (info != 0) {
            fputs("dsytri: matrix is singular\n", stderr);
            err = GretlError.singular;
        } else {
            gretl_matrix_mirror(a, uplo);
        }
    }

    return err;
}

private int real_invert_symmetric_matrix(GretlMatrixRaw * a,
                                         int symmcheck,
                                         int preserve,
                                         double * ldet)
{
    int save_nt = 0;
    integer n, info;
    double * aval = null;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    if (a.cols != a.rows) {
        fputs("invert_symmetric_matrix: input is not square\n",
              stderr);
        return GretlError.nonconf;
    }

    n = a.cols;

    if (n == 1) {
        if (ldet != null) {
            *ldet = log(a.val[0]);
        }
        a.val[0] = 1.0 / a.val[0];
        return 0;
    }

    if (symmcheck && !real_gretl_matrix_is_symmetric(a, 1)) {
        fputs("invert_symmetric_matrix: matrix is not symmetric\n", stderr);
        return GretlError.notpd;
    }

    if (preserve) {
        /* back-up, just in case */
        size_t bytes = n * n * (*aval).sizeof;

        aval = cast(double *) lapack_malloc(bytes);
        if (aval == null) {
            return GretlError.alloc;
        } else {
            memcpy(aval, a.val, bytes);
        }
    }

    if (blas_is_openblas()) {
        maybe_force_single(n, POTRF_MT_MIN, &save_nt);
    }

    dpotrf_(&uplo, &n, a.val, &n, &info);

    if (info != 0) {
        err = (info > 0) ? GretlError.notpd : GretlError.data;
        if (err == GretlError.data) {
            fprintf(stderr, "invert_symmetric_matrix: "
                    ~ "dpotrf failed with info = %d (n = %d)\n",
                    cast(int) info, cast(int) n);
        }
    }

    if (!err && ldet != null) {
        double x = 0.0;
        int i;

        for (i=0; i<n; i++) {
            x += log(gretl_matrix_get(a,i,i));
        }
        *ldet = 2.0 * x;
    }

    if (!err) {
        dpotri_(&uplo, &n, a.val, &n, &info);
        if (info != 0) {
            err = GretlError.notpd;
        } else {
            gretl_matrix_mirror(a, uplo);
        }
    }

    if (save_nt > 1) {
        blas_set_num_threads(save_nt);
    }

    if (err && preserve) {
        memcpy(a.val, aval, n * n * (*aval).sizeof);
        if (getenv("GRETL_MATRIX_DEBUG")) {
            gretl_matrix_print(a, "input matrix");
        }
    }

    if (aval != null) {
        lapack_free(aval);
    }

    return err;
}

/*
 * gretl_invert_symmetric_matrix:
 * @a: matrix to invert.
 *
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization.  On exit @a is overwritten with
 * the inverse. Uses the LAPACK functions dpotrf() and dpotri().
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_symmetric_matrix(GretlMatrixRaw * a)
{
    /* note: the second '1' for @preserve is debatable */
    return real_invert_symmetric_matrix(a, 1, 1, null);
}

/*
 * gretl_invpd:
 * @a: matrix to invert.
 *
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization.  On exit @a is overwritten with
 * the inverse. Uses the LAPACK functions dpotrf() and dpotri().
 * Little checking is done, for speed: we assume the caller
 * knows what he's doing.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invpd(GretlMatrixRaw * a)
{
    return real_invert_symmetric_matrix(a, 0, 0, null);
}

/*
 * gretl_inverse_from_cholesky_decomp:
 * @targ: matrix to hold inverse.
 * @src: Cholesky-decomposed matrix.
 *
 * Computes in @targ the inverse of a symmetric positive definite
 * matrix, on the assumption that the original matrix this has already
 * been Cholesky-decomposed in @src.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_inverse_from_cholesky_decomp(GretlMatrixRaw * targ,
                                       const(GretlMatrixRaw) * src)
{
    integer info, n;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(targ) || gretl_is_null_matrix(src)) {
        return GretlError.data;
    }

    n = src.cols;

    if (n != src.rows || n != targ.rows || n != targ.cols) {
        return GretlError.nonconf;
    }

    memcpy(targ.val, src.val, n * n * (*src.val).sizeof);

    dpotri_(&uplo, &n, targ.val, &n, &info);

    if (info != 0) {
        err = GretlError.singular;
        fprintf(stderr, "invert_symmetric_matrix:\n"
                ~ " dpotri failed with info = %d\n", cast(int) info);
    } else {
        gretl_matrix_mirror(targ, uplo);
    }

    return err;
}

/*
 * gretl_invert_symmetric_matrix2:
 * @a: matrix to invert.
 * @ldet: location to receive log determinant, or null.
 *
 * Computes the inverse of a symmetric positive definite matrix
 * using Cholesky factorization, computing the log-determinant
 * in the process.  On exit @a is overwritten with the inverse
 * and if @ldet is not null the log-determinant is written to
 * that location.  Uses the LAPACK functions dpotrf() and dpotri().
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_symmetric_matrix2(GretlMatrixRaw * a, double * ldet)
{
    return real_invert_symmetric_matrix(a, 0, 0, ldet);
}

/*
 * gretl_invert_packed_symmetric_matrix:
 * @v: symmetric matrix in vech form (lower triangle packed
 * as a column vector).
 *
 * Computes the inverse of a symmetric positive definite matrix,
 * stored in vech form, using Cholesky factorization.  On exit
 * @v is overwritten with the lower triangle of the inverse.
 * Uses the LAPACK functions dpptrf() and dpptri().
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_invert_packed_symmetric_matrix(GretlMatrixRaw * v)
{
    GretlMatrixRaw * vcpy = null;
    integer info, n;
    char uplo = 'L';
    int err = 0;

    if (gretl_is_null_matrix(v)) {
        return GretlError.data;
    }

    if (v.cols != 1) {
        fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
                ~ " matrix is not in vech form\n");
        return GretlError.data;
    }

    if (v.rows == 1) {
        v.val[0] = 1.0 / v.val[0];
        return 0;
    }

    if (v.rows < 100) {
        vcpy = gretl_matrix_copy_tmp(v);
    }

    n = cast(integer) ((sqrt(1.0 + 8.0 * v.rows) - 1.0) / 2.0);

    dpptrf_(&uplo, &n, v.val, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
                ~ " dpptrf failed with info = %d (n = %d)\n", cast(int) info, cast(int) n);
        if (info > 0) {
            fputs(" matrix is not positive definite\n", stderr);
            err = GretlError.notpd;
        } else {
            err = GretlError.data;
        }
        if (vcpy != null) {
            gretl_matrix_print(vcpy, "input matrix");
        }
        return err;
    }

    dpptri_(&uplo, &n, v.val, &info);

    if (info != 0) {
        err = GretlError.singular;
        fprintf(stderr, "gretl_invert_packed_symmetric_matrix:\n"
                ~ " dpptri failed with info = %d\n", cast(int) info);

    }

    gretl_matrix_free(vcpy);

    return err;
}

/* Below: we need an n x n complex matrix to return eigenvectors to the
   caller. The input matrix @m may or may not be correctly sized to do
   this job -- if not, we have to allocate a new matrix, which we pass
   pack via @pev. In addition we need a temporary n x n real matrix to
   pass to dgeev_ to hold a compressed representation of the
   eigenvectors; this we allocate unconditionally and pass back via
   @pec.
*/

private int dgeev_eigvecs_alloc(GretlMatrixRaw * m,
                                GretlMatrixRaw ** pev,
                                GretlMatrixRaw ** pec,
                                int n)
{
    GretlMatrixRaw * ev = null;
    GretlMatrixRaw * ec = null;
    int mrc = m.rows * m.cols;
    int dim = n * n;

    if (m.is_complex && mrc == dim) {
        /* OK */
        m.rows = m.cols = n;
    } else if (!m.is_complex && mrc == 2*dim) {
        /* @m can be converted */
        m.rows = 2*n;
        m.cols = n;
        matrix_set_complex(m, 1, 1);
    } else {
        /* we have to allocate */
        ev = gretl_cmatrix_new0(n, n);
        if (ev == null) {
            return GretlError.alloc;
        }
    }

    ec = gretl_matrix_alloc(n, n);
    if (ec == null) {
        gretl_matrix_free(ev);
        return GretlError.alloc;
    }

    *pev = ev;
    *pec = ec;

    return 0;
}

/* Transcribe from compact representation of eigenvectors in the n x n
   real matrix @src to the "new-style" n x n complex matrix @targ. What
   happens for each column depends on whether the associated eigenvalue
   is real or a member of a conjugate pair. The arrays @wr and @wi hold
   the real and imaginary parts of the eigenvalues, respectively.
*/

private void dgeev_eigvecs_transcribe(GretlMatrixRaw * targ,
                                      GretlMatrixRaw * src,
                                      double * wr, double * wi)
{
    double re, im;
    int i, j, isreal;
    int n = src.rows;

    for (j=0; j<n; j++) {
        isreal = (wi[j] == 0);
        for (i=0; i<n; i++) {
            re = gretl_matrix_get(src, i, j);
            if (isreal) {
                /* lambda(j) is real */
                gretl_cmatrix_set(targ, i, j, Cmplx(re, 0.0));
            } else {
                /* lambda(j) and lambda(j+1) are a conjugate pair */
                im = gretl_matrix_get(src, i, j+1);
                gretl_cmatrix_set(targ, i, j, Cmplx(re, im));
                gretl_cmatrix_set(targ, i, j+1, Cmplx(re, -im));
            }
        }
        if (!isreal) {
            j++;
        }
    }
}

private GretlMatrixRaw * eigen_trivial(const(GretlMatrixRaw) * A,
                                    GretlMatrixRaw * VR,
                                    GretlMatrixRaw * VL)
{
    GretlMatrixRaw * ret = gretl_matrix_copy(A);

    if (VR != null || VL != null) {
        GretlMatrixRaw *[2] targ = [VR, VL];
        GretlMatrixRaw * one;
        int i;

        for (i=0; i<2; i++) {
            if (targ[i] != null) {
                one = gretl_matrix_alloc(1, 1);
                one.val[0] = 1.0;
                gretl_matrix_replace_content(targ[i], one);
                gretl_matrix_free(one);
            }
        }
    }

    return ret;
}

/* convert dgeev eigenvalues to cmatrix format */

private void eigenvals_to_cmatrix(GretlMatrixRaw * lam,
                                  double * a, int n)
{
    int i, k = 0;

    for (i=0; i<2*n; i++) {
        a[i] = lam.val[i];
    }
    for (i=0; i<n; i++) {
        lam.val[k++] = a[i];
        lam.val[k++] = a[i+n];
    }
    lam.cols = 1;
    matrix_set_complex(lam, 1, 0);
}

private void maybe_eigen_trim(GretlMatrixRaw * lam)
{
    double * lv = lam.val + lam.rows;
    int i;

    for (i=0; i<lam.rows; i++) {
        if (lv[i] != 0.0) {
            return;
        }
    }

    /* drop the second column */
    gretl_matrix_reuse(lam, -1, 1);
}

/* The @legacy case supports some non-updated usage in var.c, kalman.c
   and genfuncs.c. But we know that @VR and @VL will be null in such
   cases: only eigenvalues are wanted.
*/

private GretlMatrixRaw * real_gretl_dgeev(const(GretlMatrixRaw) * A,
                                       GretlMatrixRaw * VR,
                                       GretlMatrixRaw * VL,
                                       int legacy,
                                       int * err)
{
    GretlMatrixRaw * ret = null;
    GretlMatrixRaw * Acpy = null;
    GretlMatrixRaw * Ltmp = null;
    GretlMatrixRaw * Rtmp = null;
    GretlMatrixRaw * VLz = null;
    GretlMatrixRaw * VRz = null;
    integer n, info, lwork;
    integer ldvl, ldvr;
    double * wr, wi;
    double * a = null;
    double * work = null;
    double * vl = null, vr = null;
    char jobvl = VL != null ? 'V' : 'N';
    char jobvr = VR != null ? 'V' : 'N';

    if (gretl_is_null_matrix(A) || A.rows != A.cols) {
        *err = GretlError.invarg;
        return null;
    }

    n = A.rows;
    if (n == 1) {
        /* Dispatch the scalar case, hence ensuring that
           A has at least two columns, which is useful
           to know below.
        */
        return eigen_trivial(A, VR, VL);
    }

    ldvl = VL != null ? n : 1;
    ldvr = VR != null ? n : 1;

    /* we need a copy of @A, which gets overwritten */
    Acpy = gretl_matrix_copy(A);
    if (Acpy == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    a = Acpy.val;

    if (VL != null) {
        /* left eigenvectors wanted */
        *err = dgeev_eigvecs_alloc(VL, &VLz, &Ltmp, n);
        if (*err) {
            goto bailout;
        }
        vl = Ltmp.val;
    }

    if (VR != null) {
        /* right eigenvectors wanted */
        *err = dgeev_eigvecs_alloc(VR, &VRz, &Rtmp, n);
        if (*err) {
            goto bailout;
        }
        vr = Rtmp.val;
    }

    work = cast(double *) lapack_malloc((*work).sizeof);
    ret = gretl_zero_matrix_new(n, 2);
    if (work == null || ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    wr = ret.val; /* eigenvalues, real part */
    wi = wr + n;  /* eigenvalues, imag. part */

    /* get optimal workspace size */
    lwork = -1;
    dgeev_(&jobvl, &jobvr, &n, a, &n, wr, wi, vl, &ldvl,
           vr, &ldvr, work, &lwork, &info);
    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* do the actual decomposition */
    dgeev_(&jobvl, &jobvr, &n, a, &n, wr, wi, vl, &ldvl,
           vr, &ldvr, work, &lwork, &info);

    if (info != 0) {
        fprintf(stderr, "dgeev: info = %d\n", info);
        *err = GretlError.invarg;
    } else {
        if (VL != null) {
            if (VLz != null) {
                dgeev_eigvecs_transcribe(VLz, Ltmp, wr, wi);
                gretl_matrix_replace_content(VL, VLz);
            } else {
                dgeev_eigvecs_transcribe(VL, Ltmp, wr, wi);
            }
        }
        if (VR != null) {
            if (VRz != null) {
                dgeev_eigvecs_transcribe(VRz, Rtmp, wr, wi);
                gretl_matrix_replace_content(VR, VRz);
            } else {
                dgeev_eigvecs_transcribe(VR, Rtmp, wr, wi);
            }
        }
    }

 bailout:

    if (*err) {
        gretl_matrix_free(ret);
        ret = null;
    } else if (legacy) {
        maybe_eigen_trim(ret);
    } else {
        eigenvals_to_cmatrix(ret, a, n);
    }

    lapack_free(work);
    gretl_matrix_free(Acpy);
    gretl_matrix_free(Ltmp);
    gretl_matrix_free(Rtmp);
    gretl_matrix_free(VLz);
    gretl_matrix_free(VRz);

    return ret;
}

GretlMatrixRaw * gretl_dgeev(const(GretlMatrixRaw) * A,
                          GretlMatrixRaw * VR,
                          GretlMatrixRaw * VL,
                          int * err)
{
    return real_gretl_dgeev(A, VR, VL, 0, err);
}

/*
 * legacy_general_matrix_eigenvals:
 * @m: square matrix on which to operate.
 * @err: location to receive error code.
 *
 * Computes the eigenvalues of the general matrix @m.
 *
 * Returns: allocated matrix containing the eigenvalues, or null
 * on failure.  The returned matrix, on successful completion,
 * is n x 2 (where n = the number of rows and columns in the
 * matrix @m); the first column holds the real parts of
 * the eigenvalues of @m and the second the imaginary parts.
 */

GretlMatrixRaw *
legacy_general_matrix_eigenvals(const(GretlMatrixRaw) * m, int * err)
{
    return real_gretl_dgeev(m, null, null, 1, err);
}

/*
 * gretl_symmetric_eigen_sort:
 * @evals: array of real eigenvalues from symmetric matrix.
 * @evecs: matrix of eigenvectors.
 * @rank: desired number of columns in output.
 *
 * Sorts the eigenvalues in @evals from largest to smallest, and
 * rearranges the columns in @evecs correspondingly.  If @rank is
 * greater than zero and less than the number of columns in @evecs,
 * then on output @evecs is shrunk so that it contains only the
 * columns associated with the largest @rank eigenvalues.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_symmetric_eigen_sort(GretlMatrixRaw * evals,
                               GretlMatrixRaw * evecs,
                               int rank)
{
    double * tmp = null;
    int n, m, err = 0;

    n = gretl_vector_get_length(evals);
    if (n == 0) {
        return GretlError.data;
    }

    if (evecs != null && (evecs.rows != n || evecs.cols != n)) {
        return GretlError.data;
    }

    if (rank <= 0) {
        rank = n;
    }
    m = n / 2;

    if (evecs != null && rank >= m) {
        /* we'll need some temporary storage for
           swapping eigenvectors
        */
        tmp = cast(double *) malloc(n * (*tmp).sizeof);
        if (tmp == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        int i, j, k;
        double x;

        /* reverse the eigenvalues in @evals */
        k = n - 1;
        for (i=0; i<m; i++) {
            x = evals.val[i];
            evals.val[i] = evals.val[k];
            evals.val[k] = x;
            k--;
        }

        if (evecs != null) {
            size_t colsize = n * (*tmp).sizeof;
            double * colj = evecs.val;
            double * colk = evecs.val + (n-1)*n;

            if (rank < m) {
                /* we just have to copy the last @rank cols
                   to the front in reverse order
                */
                m = rank;
            }

            for (j=0; j<m; j++) {
                if (tmp == null) {
                    /* col k -> col j */
                    memcpy(colj, colk, colsize);
                } else {
                    /* col j -> tmp */
                    memcpy(tmp, colj, colsize);
                    /* col k -> col j */
                    memcpy(colj, colk, colsize);
                    /* tmp -> col k */
                    memcpy(colk, tmp, colsize);
                }
                colj += n;
                colk -= n;
            }
            /* and "shrink" @evecs, if wanted */
            if (rank < n) {
                evecs.cols = rank;
            }
        }
    }

    free(tmp);

    return err;
}

private GretlMatrixRaw * eigensym_rrr(GretlMatrixRaw * m,
                                   int eigenvecs,
                                   int * err)
{
    integer n, info, lwork, liwork;
    integer nv, ldz = 1;
    double vl = 0, vu = 0;
    int il = 1, iu = m.rows;
    GretlMatrixRaw * evals = null;
    double * z = null;
    double * work = null;
    double * w = null;
    integer * iwork = null;
    integer * isuppz = null;
    char jobz = eigenvecs ? 'V' : 'N';
    double abstol = 0;
    char range = 'A';
    char uplo = 'U';

    /* Note: @vl and @vu, @il and @iu, are required to work around
       buggy implementations of dsyevr, which reference these terms
       even when they're not supposed to. @vl and @vu were added in
       October 2020 as a fix for Apple's libLAPACK.dylib; @il and @iu
       were added to prevent a crash in libopenblas from
       MSYS2/clangarm64 in May 2025.
    */

    n = m.rows;

    work = cast(double *) lapack_malloc((*work).sizeof);
    iwork = cast(integer *) malloc((*iwork).sizeof);
    if (work == null || iwork == null) {
        *err = GretlError.alloc;
        return null;
    }

    evals = gretl_column_vector_alloc(n);
    if (evals == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (eigenvecs) {
        z = cast(double *) malloc(n * n * (*z).sizeof);
        isuppz = cast(integer *) malloc(2 * n * (*isuppz).sizeof);
        if (z == null || isuppz == null) {
            *err = GretlError.alloc;
            goto bailout;
        }
        ldz = n;
    }

    w = evals.val;

    lwork = liwork = -1; /* find optimal workspace size */
    dsyevr_(&jobz, &range, &uplo, &n, m.val, &n,
            &vl, &vu, &il, &iu, &abstol, &nv, w,
            z, &ldz, isuppz, work, &lwork, iwork,
            &liwork, &info);

    if (info != 0 || work[0] <= 0.0) {
        *err = wspace_fail(info, work[0]);
        goto bailout;
    }

    lwork = cast(integer) work[0];
    liwork = iwork[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    iwork = cast(integer *) realloc(iwork, liwork * (*iwork).sizeof);
    if (work == null || iwork == null) {
        *err = GretlError.alloc;
    }

    if (!*err) {
        dsyevr_(&jobz, &range, &uplo, &n, m.val, &n,
                &vl, &vu, &il, &iu, &abstol, &nv, w,
                z, &ldz, isuppz, work, &lwork, iwork,
                &liwork, &info);
        if (info != 0) {
            fprintf(stderr, "dsyevr: info = %d\n", info);
            *err = GretlError.data;
        }
    }

    if (!*err && eigenvecs) {
        memcpy(m.val, z, n*n * (*z).sizeof);
    }

 bailout:

    lapack_free(work);
    free(iwork);
    free(isuppz);
    free(z);

    if (*err && evals != null) {
        gretl_matrix_free(evals);
        evals = null;
    }

    return evals;
}

private int basic_eigensym_work(double * mval, double * w,
                                int dim, int eigvecs)
{
    double * work = null;
    integer lwork = -1;
    integer info = 0;
    integer n = dim;
    char jobz = eigvecs ? 'V' : 'N';
    char uplo = 'U';
    int err = 0;

    work = cast(double *) lapack_malloc((*work).sizeof);
    if (work == null) {
        return GretlError.alloc;
    }

    dsyev_(&jobz, &uplo, &n, mval, &n, w, work, &lwork, &info);
    if (info != 0 || work[0] <= 0.0) {
        return wspace_fail(info, work[0]);
    }

    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        return GretlError.alloc;
    }

    dsyev_(&jobz, &uplo, &n, mval, &n, w, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "dsyev: info = %d\n", info);
        err = GretlError.invarg;
    }

    lapack_free(work);

    return err;
}

private GretlMatrixRaw * eigensym_standard(GretlMatrixRaw * m,
                                        int eigenvecs,
                                        int * err)
{
    GretlMatrixRaw * evals = null;

    evals = gretl_column_vector_alloc(m.rows);
    if (evals == null) {
        *err = GretlError.alloc;
    } else {
        *err = basic_eigensym_work(m.val, evals.val, m.rows,
                                   eigenvecs);
    }

    if (*err && evals != null) {
        gretl_matrix_free(evals);
        evals = null;
    }

    return evals;
}

/*
 * gretl_symmetric_matrix_eigenvals:
 * @m: n x n matrix to operate on.
 * @eigenvecs: non-zero to calculate eigenvectors, 0 to omit.
 * @err: location to receive error code.
 *
 * Computes the eigenvalues of the real symmetric matrix @m.
 * If @eigenvecs is non-zero, also compute the orthonormal
 * eigenvectors of @m, which are stored in @m. Uses the lapack
 * function dsyevr(), or dsyev() for small matrices.
 *
 * Returns: n x 1 matrix containing the eigenvalues in ascending
 * order, or null on failure.
 */

GretlMatrixRaw *
gretl_symmetric_matrix_eigenvals(GretlMatrixRaw * m, int eigenvecs, int * err)
{
    GretlMatrixRaw * ret = null;
    static __gshared int ev_ver = 0; /* questionable? */
    int save_nt = 0;

    *err = 0;

    if (gretl_is_null_matrix(m) || m.rows != m.cols) {
        /* If we're not actually testing for symmetry, we must
           at least test for squareness, on pain of crashing.
        */
        *err = GretlError.invarg;
        return null;
    }

    if (blas_is_threaded()) {
        save_nt = blas_get_num_threads();
        if (save_nt > 1) {
            blas_set_num_threads(1);
        }
    }

    if (ev_ver == 0) {
        char * s = getenv("GRETL_OLD_EV");

        ev_ver = s != null ? 1 : 2;
    }

    if (m.rows < 10 || ev_ver == 1) {
        ret = eigensym_standard(m, eigenvecs, err);
    } else {
        ret = eigensym_rrr(m, eigenvecs, err);
    }

    if (blas_is_threaded() && save_nt > 1) {
        blas_set_num_threads(save_nt);
    }

    return ret;
}

private GretlMatrixRaw *
real_symm_eigenvals_descending(GretlMatrixRaw * m,
                               int eigenvecs,
                               int rank,
                               int * err)
{
    GretlMatrixRaw * v =
        gretl_symmetric_matrix_eigenvals(m, eigenvecs, err);

    if (!*err) {
        m = eigenvecs ? m : null;
        *err = gretl_symmetric_eigen_sort(v, m, rank);
    }

    if (*err && v != null) {
        gretl_matrix_free(v);
        v = null;
    }

    return v;
}

/*
 * gretl_symm_matrix_eigenvals_descending:
 * @m: n x n matrix to operate on.
 * @eigenvecs: non-zero to calculate eigenvectors, 0 to omit.
 * @err: location to receive error code.
 *
 * Computes the eigenvalues of the real symmetric matrix @m.
 * If @eigenvecs is non-zero, also compute the orthonormal
 * eigenvectors of @m, which are stored in @m. Uses the lapack
 * function dsyev().
 *
 * Returns: n x 1 matrix containing the eigenvalues in descending
 * order, or null on failure.
 */

GretlMatrixRaw *
gretl_symm_matrix_eigenvals_descending(GretlMatrixRaw * m,
                                       int eigenvecs,
                                       int * err)
{
    return real_symm_eigenvals_descending(m, eigenvecs,
                                          0, err);
}

private double get_extreme_eigenvalue(GretlMatrixRaw * m, int getmax,
                                      int * err)
{
    double ev = double.nan;
    GretlMatrixRaw * v;

    v = gretl_symmetric_matrix_eigenvals(m, 0, err);

    if (!*err) {
        int n = gretl_vector_get_length(v);

        /* the eigenvalues, from lapack's dsyev(),
           are in ascending order */

        if (getmax) {
            ev = v.val[n-1];
        } else {
            ev = v.val[0];
        }

        gretl_matrix_free(v);
    }

    if (*err == 0 || *err == 1) {
        /* reconstitute full matrix */
        gretl_matrix_mirror(m, 'L');
    }

    return ev;
}

/*
 * gretl_symm_matrix_lambda_min:
 * @m: n x n matrix to operate on.
 * @err: location to receive error code.
 *
 * Returns: the minimum eigenvalue of the real symmetric matrix @m,
 * or NaN on error.
 */

double gretl_symm_matrix_lambda_min(const(GretlMatrixRaw) * m, int * err)
{
    return get_extreme_eigenvalue(cast(GretlMatrixRaw *) m, 0, err);
}

/*
 * gretl_symm_matrix_lambda_max:
 * @m: n x n matrix to operate on.
 * @err: location to receive error code.
 *
 * Returns: the maximum eigenvalue of the real symmetric matrix @m,
 * or NaN on error.
 */

double gretl_symm_matrix_lambda_max(const(GretlMatrixRaw) * m, int * err)
{
    return get_extreme_eigenvalue(cast(GretlMatrixRaw *) m, 1, err);
}

private int gensymm_conformable(const(GretlMatrixRaw) * A,
                                const(GretlMatrixRaw) * B)
{
    if (!real_gretl_matrix_is_symmetric(A, 1)) {
        fputs("gretl_gensymm_eigenvals: matrix A is not symmetric\n",
              stderr);
        return 0;
    }

    if (!real_gretl_matrix_is_symmetric(B, 1)) {
        fputs("gretl_gensymm_eigenvals: matrix B is not symmetric\n",
              stderr);
        return 0;
    }

    if (B.rows != A.rows) {
        fputs("gretl_gensymm_eigenvals: matrices A and B have different size\n",
              stderr);
        return 0;
    }

    return 1;
}

/*
 * gretl_gensymm_eigenvals:
 * @A: symmetric matrix.
 * @B: symmetric positive definite matrix.
 * @V: matrix to hold the generalized eigenvectors, or null if
 * these are not required.
 * @err: location to receive error code.
 *
 * Solves the generalized eigenvalue problem
 * | A - \lambda B | = 0 , where both A and B are symmetric
 * and B is positive definite.
 *
 * Returns: allocated storage containing the eigenvalues, in
 * ascending order, or null on failure.
 */

GretlMatrixRaw * gretl_gensymm_eigenvals(const(GretlMatrixRaw) * A,
                                      const(GretlMatrixRaw) * B,
                                      GretlMatrixRaw * V,
                                      int * err)
{
    GretlMatrixRaw * K = null;
    GretlMatrixRaw * tmp = null;
    GretlMatrixRaw * evals = null;
    int n;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
        *err = GretlError.data;
        return null;
    }

    if (!gensymm_conformable(A, B)) {
        *err = GretlError.nonconf;
        return null;
    }

    n = A.rows;
    K = gretl_matrix_copy_tmp(B);
    tmp = gretl_matrix_alloc(n, n);

    if (K == null || tmp == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    *err = gretl_matrix_cholesky_decomp(K);
    if (*err) {
        fputs("gretl_gensymm_eigenvals: matrix B not p.d.\n",
              stderr);
        *err = GretlError.nonconf;
        goto bailout;
    }

    *err = gretl_invert_triangular_matrix(K, 'L');
    if (*err) {
        fputs("gretl_gensymm_eigenvals: matrix B only p.s.d.\n",
              stderr);
        *err = GretlError.nonconf;
        goto bailout;
    }

    gretl_matrix_qform(K, GretlMatrixMod.none, A, tmp, GretlMatrixMod.none);

    evals = gretl_symmetric_matrix_eigenvals(tmp, 1, err);
    if (*err) {
        goto bailout;
    }

    if (V != null) {
        *err = gretl_matrix_multiply_mod(K, GretlMatrixMod.transpose,
                                         tmp, GretlMatrixMod.none,
                                         V, GretlMatrixMod.none);
    }

 bailout:

    gretl_matrix_free(K);
    gretl_matrix_free(tmp);

    if (*err && evals != null) {
        gretl_matrix_free(evals);
        evals = null;
    }

    return evals;
}

/* Compute SVD via eigen-decomposition for the case where
   @X is "tall": more rows than columns.
*/

private int tall_SVD(const(GretlMatrixRaw) * X,
                     GretlMatrixRaw ** pU,
                     GretlMatrixRaw ** psv,
                     GretlMatrixRaw ** pVt)
{
    GretlMatrixRaw * XTX;
    GretlMatrixRaw * sv;
    GretlMatrixRaw * lam = null;
    GretlMatrixRaw * U = null;
    GretlMatrixRaw * Vt = null;
    GretlMatrixRaw * Vl = null;
    double lj, vij;
    int vecs, c = X.cols;
    int i, j, jj;
    int err = 0;

    XTX = gretl_matrix_alloc(c, c);
    sv = gretl_matrix_alloc(1, c);
    if (XTX == null || sv == null) {
        return GretlError.alloc;
    }

    vecs = pU != null || pVt != null;

    gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                              X, GretlMatrixMod.none,
                              XTX, GretlMatrixMod.none);

    lam = gretl_symmetric_matrix_eigenvals(XTX, vecs, &err);
    if (!err) {
        for (i=0; i<c; i++) {
            lj = lam.val[c-i-1];
            if (lj < 0) {
                err = GretlError.singular;
                break;
            }
            sv.val[i] = sqrt(lj);
        }
    }

    if (err) {
        gretl_matrix_free(XTX);
        gretl_matrix_free(sv);
        return err;
    }

    if (pVt != null) {
        Vt = gretl_matrix_alloc(c, c);
        if (Vt == null) {
            err = GretlError.alloc;
        } else {
            for (j=0; j<c; j++) {
                jj = c - j - 1;
                for (i=0; i<c; i++) {
                    vij = gretl_matrix_get(XTX, i, j);
                    gretl_matrix_set(Vt, jj, i, vij);
                }
            }
        }
    }

    if (!err && pU != null) {
        U = gretl_matrix_alloc(X.rows, c);
        Vl = gretl_matrix_alloc(c, c);
        if (U == null || Vl == null) {
            err = GretlError.alloc;
        } else {
            for (j=0; j<c; j++) {
                jj = c - j - 1;
                for (i=0; i<c; i++) {
                    vij = gretl_matrix_get(XTX, i, jj);
                    gretl_matrix_set(Vl, i, j, vij / sv.val[j]);
                }
            }
            gretl_matrix_multiply(X, Vl, U);
        }
    }

    if (psv != null) {
        *psv = sv;
        sv = null;
    }
    if (pU != null) {
        *pU = U;
        U = null;
    }
    if (pVt != null) {
        *pVt = Vt;
        Vt = null;
    }

    gretl_matrix_free(XTX);
    gretl_matrix_free(sv);
    gretl_matrix_free(lam);
    gretl_matrix_free(U);
    gretl_matrix_free(Vt);
    gretl_matrix_free(Vl);

    return err;
}

private int real_gretl_matrix_SVD(const(GretlMatrixRaw) * x,
                                  GretlMatrixRaw ** pu,
                                  GretlVector ** ps,
                                  GretlMatrixRaw ** pvt,
                                  int full)
{
    integer m, n, lda;
    integer ldu = 1, ldvt = 1;
    integer lwork = -1;
    integer * iwork = null;
    integer info;
    GretlMatrixRaw * a = null;
    GretlMatrixRaw * s = null;
    GretlMatrixRaw * u = null;
    GretlMatrixRaw * vt = null;
    char jobu = 'N', jobvt = 'N';
    char jobz = 'N';
    double xu, xvt;
    double * uval = &xu, vtval = &xvt;
    double * work = null;
    int k, dnc;
    int err = 0;

    a = gretl_matrix_copy_tmp(x);
    if (a == null) {
        return GretlError.alloc;
    }

    lda = m = x.rows;
    n = x.cols;
    k = (m < n) ? m : n;
    dnc = k > 20;

    s = gretl_vector_alloc(k);
    if (s == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (dnc) {
        /* divide and conquer */
        if (pu != null || pvt != null) {
            int ucols = full ? m : k;

            ldu = m;
            ldvt = full ? n : k;
            u = gretl_matrix_alloc(ldu, ucols);
            vt = gretl_matrix_alloc(ldvt, n);
            if (u == null || vt == null) {
                err = GretlError.alloc;
                goto bailout;
            } else {
                uval = u.val;
                vtval = vt.val;
                jobz = full ? 'A' : 'S';
            }
        }

        work = cast(double *) lapack_malloc((*work).sizeof);
        iwork = cast(integer *) malloc(8 * k * (*iwork).sizeof);
        if (work == null || iwork == null) {
            err = GretlError.alloc;
            goto bailout;
        }

        /* workspace query */
        dgesdd_(&jobz, &m, &n, a.val, &lda, s.val, uval, &ldu,
                vtval, &ldvt, work, &lwork, iwork, &info);
    } else {
        /* vanilla SVD computation */
        if (pu != null) {
            ldu = m;
            if (full) {
                u = gretl_matrix_alloc(ldu, m);
            } else {
                u = gretl_matrix_alloc(ldu, k);
            }
            if (u == null) {
                err = GretlError.alloc;
                goto bailout;
            } else {
                uval = u.val;
                jobu = full ? 'A' : 'S';
            }
        }
        if (pvt != null) {
            ldvt = full ? n : k;
            vt = gretl_matrix_alloc(ldvt, n);
            if (vt == null) {
                err = GretlError.alloc;
                goto bailout;
            } else {
                vtval = vt.val;
                jobvt = full ? 'A' : 'S';
            }
        }

        work = cast(double *) lapack_malloc((*work).sizeof);
        if (work == null) {
            err = GretlError.alloc;
            goto bailout;
        }

        /* workspace query */
        dgesvd_(&jobu, &jobvt, &m, &n, a.val, &lda, s.val, uval, &ldu,
                vtval, &ldvt, work, &lwork, &info);
    }

    if (info != 0 || work[0] <= 0.0) {
        err = wspace_fail(info, work[0]);
        goto bailout;
    }

    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* actual computation */
    if (dnc) {
        dgesdd_(&jobz, &m, &n, a.val, &lda, s.val, uval, &ldu,
                vtval, &ldvt, work, &lwork, iwork, &info);
    } else {
        dgesvd_(&jobu, &jobvt, &m, &n, a.val, &lda, s.val, uval, &ldu,
                vtval, &ldvt, work, &lwork, &info);
    }

    if (info != 0) {
        fprintf(stderr, "gretl_matrix_SVD: info = %d\n", cast(int) info);
        err = GretlError.data;
        goto bailout;
    }

    if (ps != null) {
        *ps = s;
        s = null;
    }
    if (pu != null) {
        *pu = u;
        u = null;
    }
    if (pvt != null) {
        *pvt = vt;
        vt = null;
    }

 bailout:

    lapack_free(work);
    free(iwork);
    gretl_matrix_free(a);
    gretl_matrix_free(s);
    gretl_matrix_free(u);
    gretl_matrix_free(vt);

    return err;
}

/*
 * gretl_matrix_SVD:
 * @x: m x n matrix to decompose.
 * @pu: location for matrix U, or null if not wanted.
 * @ps: location for vector of singular values, or null if not wanted.
 * @pvt: location for matrix V (transposed), or null if not wanted.
 * @full: if U and/or V are to be computed, a non-zero value flags
 * production of "full-size" U (m x m) and/or V (n x n). Otherwise U
 * will be m x min(m,n) and and V' will be min(m,n) x n. Note that
 * this flag matters only if @x is not square.
 *
 * Computes SVD factorization of a general matrix using one of the
 * the lapack functions dgesvd() or dgesdd(). A = U * diag(s) * Vt.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_SVD(const(GretlMatrixRaw) * x, GretlMatrixRaw ** pu,
                     GretlVector ** ps, GretlMatrixRaw ** pvt,
                     int full)
{
    int err = 0;

    if (pu == null && ps == null && pvt == null) {
        /* no-op */
        return 0;
    } else if (gretl_is_null_matrix(x)) {
        return GretlError.data;
    }

    if (!full && x.rows > x.cols && getenv("GRETL_REAL_SVD") == null) {
        /* The "tall" variant is very fast, but not at all
           accurate for near-singular matrices. If @x is
           too close to singular this will be flagged by an
           error code of GretlError.singular from tall_SVD(), in which
           case we'll proceed to try "real" SVD; any other
           error will be treated as fatal.
        */
        err = tall_SVD(x, pu, ps, pvt);
        if (err != GretlError.singular) {
            /* either OK or fatal error */
            return err;
        }
    }

    return real_gretl_matrix_SVD(x, pu, ps, pvt, full);
}

/*
 * gretl_matrix_SVD_johansen_solve:
 * @R0: T x p matrix of residuals.
 * @R1: T x p1 matrix of residuals.
 * @evals: vector to receive eigenvals, or null if not wanted.
 * @B: matrix to hold \beta, or null if not wanted.
 * @A: matrix to hold \alpha, or null if not wanted.
 * @jrank: cointegration rank, <= p.
 *
 * Solves the Johansen generalized eigenvalue problem via
 * SVD decomposition.  See J. A. Doornik and R. J. O'Brien,
 * "Numerically stable cointegration analysis", Computational
 * Statistics and Data Analysis, 41 (2002), pp. 185-193,
 * Algorithm 4.
 *
 * If @B is non-null it should be p1 x p on input; it will
 * be trimmed to p1 x @jrank on output if @jrank < p.
 * If @A is non-null it should be p x p on input; it will
 * be trimmed to p x @jrank on output if @jrank < p.
 * @evals should be a vector of length @jrank.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_SVD_johansen_solve(const(GretlMatrixRaw) * R0,
                                    const(GretlMatrixRaw) * R1,
                                    GretlMatrixRaw * evals,
                                    GretlMatrixRaw * B,
                                    GretlMatrixRaw * A,
                                    int jrank)
{
    GretlMatrixRaw * U0 = null;
    GretlMatrixRaw * U1 = null;
    GretlMatrixRaw * Uz = null;
    GretlMatrixRaw * S1 = null;
    GretlMatrixRaw * Sz = null;
    GretlMatrixRaw * V1 = null;
    GretlMatrixRaw * Z = null;
    int T = R0.rows;
    int p = R0.cols;
    int p1 = R1.cols;
    int r, err;

    if (evals == null && B == null && A == null) {
        /* no-op */
        return 0;
    }

    r = (jrank == 0) ? p : jrank;

    if (r < 1 || r > p) {
        fprintf(stderr, "Johansen SVD: r is wrong (%d)\n", r);
        return GretlError.nonconf;
    }

    if (evals != null && gretl_vector_get_length(evals) < r) {
        fprintf(stderr, "Johansen SVD: evals is too short\n");
        return GretlError.nonconf;
    }

    if (B != null && (B.rows != p1 || B.cols != p)) {
        fprintf(stderr, "Johansen SVD: B is wrong size\n");
        return GretlError.nonconf;
    }

    if (A != null && (A.rows != p || A.cols != p)) {
        fprintf(stderr, "Johansen SVD: A is wrong size\n");
        return GretlError.nonconf;
    }

    err = real_gretl_matrix_SVD(R0, &U0, null, null, 0);

    if (!err) {
        err = real_gretl_matrix_SVD(R1, &U1, &S1, &V1, 0);
    }

    if (!err) {
        Z = gretl_matrix_alloc(p1, p);
        if (Z == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_matrix_multiply_mod(U1, GretlMatrixMod.transpose,
                                            U0, GretlMatrixMod.none,
                                            Z, GretlMatrixMod.none);
        }
    }

    if (!err) {
        err = real_gretl_matrix_SVD(Z, &Uz, &Sz, null, 0);
    }

    if (!err) {
        double x, si;
        int i, j;

        if (evals != null) {
            for (i=0; i<r; i++) {
                evals.val[i] = Sz.val[i] * Sz.val[i];
            }
        }

        if (B != null) {
            /* \hat{\beta} = T^{1/2} V_1 {\Sigma_1}^{-1} U_z */
            for (i=0; i<p1; i++) {
                si = S1.val[i];
                for (j=0; j<p1; j++) {
                    if (si > SVD_SMIN) {
                        x = gretl_matrix_get(V1, i, j);
                        gretl_matrix_set(V1, i, j, x / si);
                    } else {
                        gretl_matrix_set(V1, i, j, 0);
                    }
                }
            }

            gretl_matrix_multiply_mod(V1, GretlMatrixMod.transpose,
                                      Uz, GretlMatrixMod.none,
                                      B, GretlMatrixMod.none);
            gretl_matrix_multiply_by_scalar(B, sqrt(cast(double) T));
            if (r < p) {
                gretl_matrix_reuse(B, -1, r);
            }
        }

        if (A != null) {
            /* \hat{\alpha} = T^{-1/2} R_0' U_1 U_z */
            gretl_matrix_reuse(Z, p, p1);
            gretl_matrix_multiply_mod(R0, GretlMatrixMod.transpose,
                                      U1, GretlMatrixMod.none,
                                      Z, GretlMatrixMod.none);
            gretl_matrix_multiply(Z, Uz, A);
            gretl_matrix_divide_by_scalar(A, sqrt(cast(double) T));
            if (r < p) {
                gretl_matrix_reuse(A, -1, r);
            }
        }
    }

    gretl_matrix_free(U0);
    gretl_matrix_free(U1);
    gretl_matrix_free(Uz);
    gretl_matrix_free(S1);
    gretl_matrix_free(Sz);
    gretl_matrix_free(V1);
    gretl_matrix_free(Z);

    return err;
}

/* just remove ugliness for printing (in the C source this is the
   non-OLD_NULLSPACE variant; OLD_NULLSPACE is 0) */

private void normalize_nullspace(GretlMatrixRaw * M)
{
    int i, k = M.rows * M.cols;

    for (i=0; i<k; i++) {
        if (M.val[i] == -0) {
            M.val[i] = 0;
        }
    }
}

/*
 * gretl_matrix_right_nullspace:
 * @M: matrix to operate on.
 * @err: location to receive error code.
 *
 * Given an m x n matrix @M, construct a conformable matrix
 * R such that MR = 0 (that is, all the columns of R are
 * orthogonal to the space spanned by the rows of @M).
 *
 * Returns: the allocated matrix R, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_right_nullspace(const(GretlMatrixRaw) * M, int * err)
{
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * S = null;
    int i, j, k;

    if (gretl_is_null_matrix(M)) {
        *err = GretlError.data;
        return null;
    }

    /* we'll need the full SVD here */
    *err = real_gretl_matrix_SVD(M, null, &S, &V, 1);

    if (!*err) {
        char E = 'E';
        int m = M.rows;
        int n = M.cols;
        int r = min(m, n);
        int sz = max(m, n);
        double x, eps = dlamch_(&E);
        double smin = sz * S.val[0] * eps;

        /* rank plus nullity = n */
        k = n;
        for (i=0; i<r; i++) {
            if (S.val[i] > smin) {
                k--;
            }
        }

        if (k == 0) {
            R = gretl_null_matrix_new();
        } else {
            R = gretl_matrix_alloc(n, k);
        }

        if (R == null) {
            *err = GretlError.alloc;
        } else if (k > 0) {
            for (i=0; i<n; i++) {
                for (j=0; j<k; j++) {
                    x = gretl_matrix_get(V, j + n - k, i);
                    gretl_matrix_set(R, i, j, x);
                }
            }
            normalize_nullspace(R);
        }
    }

    gretl_matrix_free(S);
    gretl_matrix_free(V);

    return R;
}

/*
 * gretl_matrix_left_nullspace:
 * @M: matrix to operate on.
 * @mod: GretlMatrixMod.none or GretlMatrixMod.transpose
 * @err: location to receive error code.
 *
 * Given an m x n matrix @M, construct a conformable matrix
 * L such that LM = 0 (that is, all the columns of @M are
 * orthogonal to the space spanned by the rows of L).
 *
 * Returns: the allocated matrix L, or if @mod is
 * GretlMatrixMod.transpose, L', or null on failure.
 */

GretlMatrixRaw * gretl_matrix_left_nullspace(const(GretlMatrixRaw) * M,
                                          GretlMatrixMod mod,
                                          int * err)
{
    GretlMatrixRaw * Tmp = null;
    GretlMatrixRaw * L = null;

    if (gretl_is_null_matrix(M)) {
        *err = GretlError.data;
        return null;
    }

    Tmp = gretl_matrix_copy_transpose(M);
    if (Tmp == null) {
        *err = GretlError.alloc;
        return null;
    }

    L = gretl_matrix_right_nullspace(Tmp, err);
    gretl_matrix_free(Tmp);

    if (!*err && mod == GretlMatrixMod.transpose) {
        Tmp = gretl_matrix_copy_transpose(L);
        if (Tmp == null) {
            *err = GretlError.alloc;
        } else {
            gretl_matrix_free(L);
            L = Tmp;
        }
    }

    return L;
}

private bool true_null_matrix(const(GretlMatrixRaw) * a)
{
    return (a.rows == 0 && a.cols == 0);
}

/*
 * gretl_matrix_row_concat:
 * @a: upper source matrix (m x n).
 * @b: lower source matrix (p x n).
 * @err: location to receive error code.
 *
 * Returns: newly allocated matrix ((m+p) x n) that results from
 * the row-wise concatenation of @a and @b, or null on failure.
 */

GretlMatrixRaw *
gretl_matrix_row_concat(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b,
                        int * err)
{
    GretlMatrixRaw * c = null;

    if (a == null || b == null) {
        *err = GretlError.data;
    } else if (true_null_matrix(a)) {
        c = gretl_matrix_copy(b);
        goto finish;
    } else if (true_null_matrix(b)) {
        c = gretl_matrix_copy(a);
        goto finish;
    }

    if (!*err) {
        int cmplx_a = a.is_complex;
        int cmplx_b = b.is_complex;
        int cmplx_c = cmplx_a || cmplx_b;
        int scalar_a = 0;
        int scalar_b = 0;
        Cmplx z;
        double x;
        int cr, cc;
        int i, j, k;

        if (matrix_is_scalar(a) && b.cols > 1) {
            /* 2024-11-29: was b->cols != 1 above */
            scalar_a = 1;
            cr = b.rows + 1;
            cc = b.cols;
        } else if (matrix_is_scalar(b) && a.cols > 1) {
            /* 2024-11-29: was a->cols != 1 above */
            scalar_b = 1;
            cr = a.rows + 1;
            cc = a.cols;
        } else if (a.cols != b.cols) {
            *err = GretlError.nonconf;
            return null;
        } else {
            cr = a.rows + b.rows;
            cc = a.cols;
        }

        if (cr == 0 && cc == 0) {
            c = gretl_null_matrix_new();
        } else if (cmplx_c) {
            c = gretl_cmatrix_new(cr, cc);
        } else {
            c = gretl_matrix_alloc(cr, cc);
        }
        if (c == null) {
            *err = GretlError.alloc;
            return null;
        } else if (cr == 0) {
            return c;
        }

        if (scalar_a) {
            for (j=0; j<b.cols; j++) {
                if (cmplx_c) {
                    z = cmplx_a ? a.z[0] : Cmplx(a.val[0], 0.0);
                    gretl_cmatrix_set(c, 0, j, z);
                } else {
                    gretl_matrix_set(c, 0, j, a.val[0]);
                }
            }
        } else {
            for (i=0; i<a.rows; i++) {
                for (j=0; j<a.cols; j++) {
                    if (cmplx_c) {
                        z = cmplx_a ? gretl_cmatrix_get(a, i, j) :
                            Cmplx(gretl_matrix_get(a, i, j), 0.0);
                        gretl_cmatrix_set(c, i, j, z);
                    } else {
                        x = gretl_matrix_get(a, i, j);
                        gretl_matrix_set(c, i, j, x);
                    }
                }
            }
        }

        k = a.rows;
        if (scalar_b) {
            for (j=0; j<a.cols; j++) {
                if (cmplx_c) {
                    z = cmplx_b ? b.z[0] : Cmplx(b.val[0], 0.0);
                    gretl_cmatrix_set(c, k, j, z);
                } else {
                    gretl_matrix_set(c, k, j, b.val[0]);
                }
            }
        } else {
            for (i=0; i<b.rows; i++) {
                for (j=0; j<b.cols; j++) {
                    if (cmplx_c) {
                        z = cmplx_b ? gretl_cmatrix_get(b, i, j) :
                            Cmplx(gretl_matrix_get(b, i, j), 0.0);
                        gretl_cmatrix_set(c, k, j, z);
                    } else {
                        x = gretl_matrix_get(b, i, j);
                        gretl_matrix_set(c, k, j, x);
                    }
                }
                k++;
            }
        }
    }

 finish:

    if (!*err) {
        if (c == null) {
            *err = GretlError.alloc;
        } else {
            maybe_preserve_names(c, a, NameFlags.colnames, null);
            maybe_concat_names(c, a, b, NameFlags.rownames);
        }
    }

    return c;
}

/*
 * gretl_matrix_col_concat:
 * @a: left-hand source matrix (m x n).
 * @b: right-hand source matrix (m x p).
 * @err: location to receive error code.
 *
 * Returns: newly allocated matrix (m x (n+p)) that results from
 * the column-wise concatenation of @a and @b, or null on failure.
 */

GretlMatrixRaw *
gretl_matrix_col_concat(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b,
                        int * err)
{
    GretlMatrixRaw * c = null;

    if (a == null || b == null) {
        *err = GretlError.data;
    } else if (true_null_matrix(a)) {
        c = gretl_matrix_copy(b);
        goto finish;
    } else if (true_null_matrix(b)) {
        c = gretl_matrix_copy(a);
        goto finish;
    }

    if (!*err) {
        int cmplx_a = a.is_complex;
        int cmplx_b = b.is_complex;
        int cmplx_c = cmplx_a || cmplx_b;
        int scalar_a = 0;
        int scalar_b = 0;
        int n_a = a.rows * a.cols;
        int n_b = b.rows * b.cols;
        size_t xsize = double.sizeof;
        size_t zsize = Cmplx.sizeof;
        Cmplx z;
        int i, cr, cc;

        if (matrix_is_scalar(a) && b.rows > 1) {
            scalar_a = 1;
            cr = b.rows;
            cc = b.cols + 1;
        } else if (matrix_is_scalar(b) && a.rows > 1) {
            scalar_b = 1;
            cr = a.rows;
            cc = a.cols + 1;
        } else if (a.rows != b.rows) {
            *err = GretlError.nonconf;
            return null;
        } else {
            cr = a.rows;
            cc = a.cols + b.cols;
        }

        if (cr == 0 && cc == 0) {
            c = gretl_null_matrix_new();
        } else if (cmplx_c) {
            c = gretl_cmatrix_new(cr, cc);
        } else {
            c = gretl_matrix_alloc(cr, cc);
        }
        if (c == null) {
            *err = GretlError.alloc;
            return null;
        } else if (cr == 0) {
            return c;
        }

        if (scalar_a) {
            if (!cmplx_c) {
                memcpy(c.val + b.rows, b.val, n_b * xsize);
            } else if (cmplx_b) {
                memcpy(c.z + b.rows, b.z, n_b * zsize);
            } else {
                real_to_complex_fill(c, b, 0, 1);
            }
            for (i=0; i<b.rows; i++) {
                if (cmplx_c) {
                    z = cmplx_a ? a.z[0] : Cmplx(a.val[0], 0.0);
                    gretl_cmatrix_set(c, i, 0, z);
                } else {
                    gretl_matrix_set(c, i, 0, a.val[0]);
                }
            }
        } else if (scalar_b) {
            if (!cmplx_c) {
                memcpy(c.val, a.val, n_a * xsize);
            } else if (cmplx_a) {
                memcpy(c.z, a.z, n_a * zsize);
            } else {
                real_to_complex_fill(c, a, 0, 0);
            }
            for (i=0; i<a.rows; i++) {
                if (cmplx_c) {
                    z = cmplx_b ? b.z[0] : Cmplx(b.val[0], 0.0);
                    gretl_cmatrix_set(c, i, a.cols, z);
                } else {
                    gretl_matrix_set(c, i, a.cols, b.val[0]);
                }
            }
        } else {
            /* neither @a nor @b is scalar */
            if (!cmplx_c) {
                memcpy(c.val, a.val, n_a * xsize);
                memcpy(c.val + n_a, b.val, n_b * xsize);
            } else {
                if (cmplx_a) {
                    memcpy(c.z, a.z, n_a * zsize);
                } else {
                    real_to_complex_fill(c, a, 0, 0);
                }
                if (cmplx_b) {
                    memcpy(c.z + n_a, b.z, n_b * zsize);
                } else {
                    real_to_complex_fill(c, b, 0, a.cols);
                }
            }
        }
    }

 finish:

    if (!*err) {
        if (c == null) {
            *err = GretlError.alloc;
        } else {
            maybe_preserve_names(c, a, NameFlags.rownames, null);
            maybe_concat_names(c, a, b, NameFlags.colnames);
        }
    }

    return c;
}

/*
 * gretl_matrix_direct_sum:
 * @a: top left matrix.
 * @b: bottom right matrix.
 * @err: location to receive error code.
 *
 * Returns: a new matrix containing the direct sum of @a and
 * @b, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_direct_sum(const(GretlMatrixRaw) * a,
                                      const(GretlMatrixRaw) * b,
                                      int * err)
{
    GretlMatrixRaw * c = null;

    if (gretl_is_null_matrix(a) && gretl_is_null_matrix(b)) {
        c = gretl_null_matrix_new();
    } else if (a.is_complex + b.is_complex == 1) {
        *err = GretlError.mixed;
    } else if (gretl_is_null_matrix(a)) {
        c = gretl_matrix_copy(b);
    } else if (gretl_is_null_matrix(b)) {
        c = gretl_matrix_copy(a);
    } else {
        int m = a.rows + b.rows;
        int n = a.cols + b.cols;
        int i, j, k;
        Cmplx z;
        double x;

        if (a.is_complex) {
            c = gretl_cmatrix_new0(m, n);
        } else {
            c = gretl_zero_matrix_new(m, n);
        }

        if (c != null) {
            for (i=0; i<a.rows; i++) {
                for (j=0; j<a.cols; j++) {
                    if (a.is_complex) {
                        z = gretl_cmatrix_get(a, i, j);
                        gretl_cmatrix_set(c, i, j, z);
                    } else {
                        x = gretl_matrix_get(a, i, j);
                        gretl_matrix_set(c, i, j, x);
                    }
                }
            }
            for (i=0; i<b.rows; i++) {
                k = i + a.rows;
                for (j=0; j<b.cols; j++) {
                    if (a.is_complex) {
                        z = gretl_cmatrix_get(b, i, j);
                        gretl_cmatrix_set(c, k, j + a.cols, z);
                    } else {
                        x = gretl_matrix_get(b, i, j);
                        gretl_matrix_set(c, k, j + a.cols, x);
                    }
                }
            }
        }
    }

    if (!*err && c == null) {
        *err = GretlError.alloc;
    }

    return c;
}

/*
 * gretl_matrix_inplace_colcat:
 * @a: matrix to be enlarged (m x n).
 * @b: matrix from which columns should be added (m x p).
 * @mask: char array, of length p, with 1s in positions
 * corresponding to columns of @b that are to be added
 * to @a, 0s elsewhere; or null to add all columns
 * of @b.
 *
 * Concatenates onto @a the selected columns of @b, if the
 * two matrices are conformable.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_matrix_inplace_colcat(GretlMatrixRaw * a,
                                const(GretlMatrixRaw) * b,
                                const(char) * mask)
{
    double x;
    int addc;
    int i, j, k;

    if (a == null || b == null) {
        return GretlError.data;
    } else if (a.is_complex || b.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_inplace_colcat\n");
        return GretlError.cmplx;
    } else if (a.rows != b.rows) {
        return GretlError.nonconf;
    }

    if (mask == null) {
        addc = b.cols;
    } else {
        addc = 0;
        for (j=0; j<b.cols; j++) {
            if (mask[j]) addc++;
        }
        if (addc == 0) {
            return 0;
        }
    }

    k = a.cols;

    if (gretl_matrix_realloc(a, a.rows, k + addc)) {
        return GretlError.alloc;
    }

    if (mask == null) {
        size_t bsize = b.rows * b.cols * (*b.val).sizeof;

        memcpy(a.val + a.rows * k, b.val, bsize);
    } else {
        for (j=0; j<b.cols; j++) {
            if (mask[j]) {
                for (i=0; i<b.rows; i++) {
                    x = gretl_matrix_get(b, i, j);
                    gretl_matrix_set(a, i, k, x);
                }
                k++;
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_cumcol:
 * @m: source matrix.
 * @err: error code.
 *
 * Returns: a matrix of the same dimensions as @m, containing
 * the cumulated columns of @m.
 */

GretlMatrixRaw * gretl_matrix_cumcol(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * a;
    int t, i;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
        return null;
    }

    a = gretl_matching_matrix_new(m.rows, m.cols, m);

    if (a == null) {
        *err = GretlError.alloc;
    } else if (a.is_complex) {
        Cmplx z;

        for (i=0; i<m.cols; i++) {
            z = Cmplx(0, 0);
            for (t=0; t<m.rows; t++) {
                z += gretl_cmatrix_get(m, t, i);
                gretl_cmatrix_set(a, t, i, z);
            }
        }
    } else {
        double x;

        for (i=0; i<m.cols; i++) {
            x = 0;
            for (t=0; t<m.rows; t++) {
                x += gretl_matrix_get(m, t, i);
                gretl_matrix_set(a, t, i, x);
            }
        }
    }

    return a;
}

/*
 * gretl_matrix_diffcol:
 * @m: source matrix.
 * @missval: value to represent missing observations.
 * @err: error code.
 *
 * Returns: a matrix of the same dimensions as @m, containing
 * @missval in the first row and the difference between consecutive
 * rows of @m afterwards.
 */

GretlMatrixRaw * gretl_matrix_diffcol(const(GretlMatrixRaw) * m,
                                   double missval, int * err)
{
    GretlMatrixRaw * a;
    int t, i;

    *err = 0;

    if (gretl_is_null_matrix(m)) {
        return null;
    }

    a = gretl_matching_matrix_new(m.rows, m.cols, m);

    if (a == null) {
        *err = GretlError.alloc;
    } else if (a.is_complex) {
        Cmplx z, zlag;

        for (i=0; i<m.cols; i++) {
            gretl_cmatrix_set(a, 0, i, Cmplx(missval, 0.0));
        }
        for (i=0; i<m.cols; i++) {
            zlag = gretl_cmatrix_get(m, 0, i);
            for (t=1; t<m.rows; t++) {
                z = gretl_cmatrix_get(m, t, i);
                gretl_cmatrix_set(a, t, i, z - zlag);
                zlag = z;
            }
        }
    } else {
        double x, xlag;

        for (i=0; i<m.cols; i++) {
            gretl_matrix_set(a, 0, i, missval);
        }
        for (i=0; i<m.cols; i++) {
            xlag = gretl_matrix_get(m, 0, i);
            for (t=1; t<m.rows; t++) {
                x = gretl_matrix_get(m, t, i);
                gretl_matrix_set(a, t, i, x - xlag);
                xlag = x;
            }
        }
    }

    return a;
}

/*
 * gretl_matrix_lag:
 * @m: source matrix.
 * @k: vector of lag orders (> 0 for lags, < 0 for leads).
 * @opt: use OPT_L to arrange multiple lags by lag rather than by variable.
 * @missval: value to represent missing observations.
 *
 * Returns: A matrix of the same dimensions as @m, containing lags
 * of the variables in the columns of @m, with missing values set
 * to @missval.
 */

GretlMatrixRaw * gretl_matrix_lag(const(GretlMatrixRaw) * m,
                               const(GretlVector) * k,
                               Gretlopt opt,
                               double missval)
{
    GretlMatrixRaw * a;
    double x;
    int nl = gretl_vector_get_length(k);
    int s, t, i, j, n, kj;

    if (gretl_is_null_matrix(m) || nl == 0 || m.is_complex) {
        return null;
    }

    a = gretl_matrix_alloc(m.rows, m.cols * nl);
    if (a == null) {
        return null;
    }

    if (opt & Gretlopt.OPT_L) {
        /* by lag */
        n = 0;
        for (j=0; j<nl; j++) {
            kj = cast(int) gretl_vector_get(k, j);
            for (t=0; t<m.rows; t++) {
                s = t - kj;
                if (s < 0 || s >= m.rows) {
                    for (i=0; i<m.cols; i++) {
                        gretl_matrix_set(a, t, n+i, missval);
                    }
                } else {
                    for (i=0; i<m.cols; i++) {
                        x = gretl_matrix_get(m, s, i);
                        gretl_matrix_set(a, t, n+i, x);
                    }
                }
            }
            n += m.cols;
        }
    } else {
        /* by variable */
        n = 0;
        for (i=0; i<m.cols; i++) {
            for (j=0; j<nl; j++) {
                kj = cast(int) gretl_vector_get(k, j);
                for (t=0; t<m.rows; t++) {
                    s = t - kj;
                    if (s < 0 || s >= m.rows) {
                        gretl_matrix_set(a, t, n+j, missval);
                    } else {
                        x = gretl_matrix_get(m, s, i);
                        gretl_matrix_set(a, t, n+j, x);
                    }
                }
            }
            n += nl;
        }
    }

    return a;
}

/*
 * gretl_matrix_inplace_lag:
 * @targ: target matrix.
 * @src: source matrix.
 * @k: lag order (> 0 for lags, < 0 for leads).
 *
 * Fills out @targ (if it is of the correct dimensions),
 * with (columnwise) lags of @src, using 0 for missing
 * values.
 *
 * Returns: 0 on success, non-zero code otherwise.
 */

int gretl_matrix_inplace_lag(GretlMatrixRaw * targ,
                             const(GretlMatrixRaw) * src,
                             int k)
{
    int m, n;
    double x;
    int s, t, i;

    if (gretl_is_null_matrix(targ) || gretl_is_null_matrix(src)) {
        return GretlError.data;
    }

    m = src.rows;
    n = src.cols;

    if (targ.rows != m || targ.cols != n) {
        return GretlError.nonconf;
    }

    for (t=0; t<m; t++) {
        s = t - k;
        if (s < 0 || s >= m) {
            for (i=0; i<n; i++) {
                gretl_matrix_set(targ, t, i, 0.0);
            }
        } else {
            for (i=0; i<n; i++) {
                x = gretl_matrix_get(src, s, i);
                gretl_matrix_set(targ, t, i, x);
            }
        }
    }

    return 0;
}

/* The most common use-case here will be updating the t1 and t2 members
   of targ's info based on src's info: this naturally arises when a new
   m x n matrix is generated and its content is assigned to an existing
   m x n matrix. We'll also handle the passing on of rownames and/or
   colnames to the copy, taking transposition into account.
*/

private int gretl_matrix_copy_info(GretlMatrixRaw * targ,
                                   const(GretlMatrixRaw) * src,
                                   int mod)
{
    int err = 0;

    if (is_block_matrix(targ) || is_block_matrix(src)) {
        return GretlError.data;
    }

    if (src.info == null || src.is_complex) {
        if (targ.info != null) {
            gretl_matrix_destroy_info(targ);
        }
        return 0;
    }

    if (targ.info == null) {
        targ.info = cast(MatrixInfo *) malloc((*targ.info).sizeof);
        if (targ.info == null) {
            return GretlError.alloc;
        }
    } else {
        strings_array_free(targ.info.colnames, targ.cols);
        strings_array_free(targ.info.rownames, targ.rows);
    }

    targ.info.colnames = null;
    targ.info.rownames = null;

    if (mod == GretlMatrixMod.transpose) {
        if (src.info.colnames != null) {
            targ.info.rownames = strings_array_dup(cast(char **) src.info.colnames,
                                                   src.cols);
            if (targ.info.rownames == null) {
                err = GretlError.alloc;
            }
        }
        if (!err && src.info.rownames != null) {
            targ.info.colnames = strings_array_dup(cast(char **) src.info.rownames,
                                                   src.rows);
            if (targ.info.colnames == null) {
                err = GretlError.alloc;
            }
        }
    } else {
        targ.info.t1 = src.info.t1;
        targ.info.t2 = src.info.t2;
        if (src.info.colnames != null) {
            targ.info.colnames = strings_array_dup(cast(char **) src.info.colnames,
                                                   src.cols);
            if (targ.info.colnames == null) {
                err = GretlError.alloc;
            }
        }
        if (!err && src.info.rownames != null) {
            targ.info.rownames = strings_array_dup(cast(char **) src.info.rownames,
                                                   src.rows);
            if (targ.info.rownames == null) {
                err = GretlError.alloc;
            }
        }
    }

    return err;
}

private int gretl_matrix_add_info(GretlMatrixRaw * m)
{
    m.info = cast(MatrixInfo *) malloc((*m.info).sizeof);

    if (m.info == null) {
        return GretlError.alloc;
    } else {
        m.info.t1 = 0;
        m.info.t2 = 0;
        m.info.colnames = null;
        m.info.rownames = null;
        return 0;
    }
}

/*
 * gretl_matrix_set_t1:
 * @m: matrix to operate on.
 * @t: integer value to set.
 *
 * Sets an integer value on @m, which can be retrieved using
 * gretl_matrix_get_t1().
 *
 * Returns: 0 on success, non-ero on error.
 */

int gretl_matrix_set_t1(GretlMatrixRaw * m, int t)
{
    if (m == null) {
        return GretlError.data;
    } else if (is_block_matrix(m)) {
        return matrix_block_error("gretl_matrix_set_t1");
    } else if (m.info == null && gretl_matrix_add_info(m)) {
        return GretlError.alloc;
    }

    m.info.t1 = t;

    return 0;
}

/*
 * gretl_matrix_set_t2:
 * @m: matrix to operate on.
 * @t: integer value to set.
 *
 * Sets an integer value on @m, which can be retrieved using
 * gretl_matrix_get_t2().
 *
 * Returns: 0 on success, non-ero on error.
 */

int gretl_matrix_set_t2(GretlMatrixRaw * m, int t)
{
    if (m == null) {
        return GretlError.data;
    } else if (is_block_matrix(m)) {
        return matrix_block_error("gretl_matrix_set_t2");
    } else if (m.info == null && gretl_matrix_add_info(m)) {
        return GretlError.alloc;
    }

    m.info.t2 = t;

    return 0;
}

/*
 * gretl_matrix_get_t1:
 * @m: matrix to read from.
 *
 * Returns: the integer that has been set on @m using
 * gretl_matrix_set_t1(), or zero if no such value has
 * been set.
 */

int gretl_matrix_get_t1(const(GretlMatrixRaw) * m)
{
    if (m != null && !is_block_matrix(m) && m.info != null) {
        return m.info.t1;
    } else {
        return 0;
    }
}

/*
 * gretl_matrix_get_t2:
 * @m: matrix to read from.
 *
 * Returns: the integer that has been set on @m using
 * gretl_matrix_set_t2(), or zero if no such value has
 * been set.
 */

int gretl_matrix_get_t2(const(GretlMatrixRaw) * m)
{
    if (m != null && !is_block_matrix(m) && m.info != null) {
        return m.info.t2;
    } else {
        return 0;
    }
}

/*
 * gretl_matrix_is_dated:
 * @m: matrix to examine.
 *
 * Returns: 1 if matrix @m has integer indices recorded
 * via gretl_matrix_set_t1() and gretl_matrix_set_t2(),
 * such that t1 >= 0 and t2 > t1, otherwise zero.
 */

int gretl_matrix_is_dated(const(GretlMatrixRaw) * m)
{
    if (m != null && !is_block_matrix(m) && m.info != null) {
        return (m.info.t1 >= 0 && (m.info.t2 > m.info.t1));
    } else {
        return 0;
    }
}

/*
 * gretl_matrix_copy_dates:
 * @target: target matrix.
 * @src: source_matrix.
 *
 * Copies the "date" elements t1 and t2 from @src to @targ,
 * provided that @src has dates set and the two matrices have
 * the same number of rows.

 * Returns: 0 on successful completion or an error code
 * on failure.
 */

int gretl_matrix_copy_dates(GretlMatrixRaw * targ,
                            const(GretlMatrixRaw) * src)
{
    int err;

    if (targ.rows != src.rows) {
        err = GretlError.nonconf;
    } else if (!gretl_matrix_is_dated(src)) {
        err = GretlError.data;
    } else {
        err = gretl_matrix_set_t1(targ, src.info.t1);
    }
    if (!err) {
        gretl_matrix_set_t2(targ, src.info.t2);
    }

    return err;
}

private int
get_SVD_ols_vcv(const(GretlMatrixRaw) * A, const(GretlMatrixRaw) * B,
                const(double) * s, GretlMatrixRaw * V, double * s2)
{
    double aik, ajk, vij;
    int m = A.cols;
    int i, j, k;

    /* Get X'X{-1}, based on the work done by the SV decomp:
       reciprocals of the squares of the (positive) singular values,
       premultiplied by V and postmultiplied by V-transpose
    */
    for (i=0; i<m; i++) {
        for (j=i; j<m; j++) {
            vij = 0.0;
            for (k=0; k<m; k++) {
                if (s[k] > 0.0) {
                    aik = gretl_matrix_get(A, k, i);
                    ajk = gretl_matrix_get(A, k, j);
                    vij += aik * ajk / (s[k] * s[k]);
                }
            }
            gretl_matrix_set(V, i, j, vij);
            if (j != i) {
                gretl_matrix_set(V, j, i, vij);
            }
        }
    }

    if (s2 != null) {
        double sigma2 = 0.0;
        int T = A.rows;

        for (i=m; i<T; i++) {
            sigma2 += B.val[i] * B.val[i];
        }
        sigma2 /= T - m;
        gretl_matrix_multiply_by_scalar(V, sigma2);
        *s2 = sigma2;
    }

    return 0;
}

private double
get_ols_error_variance(const(GretlVector) * y, const(GretlMatrixRaw) * X,
                       const(GretlVector) * b, int nr)
{
    double u, s2 = 0.0;
    int k = X.cols;  /* number of regressors */
    int n = X.rows;  /* number of observations */
    int i, j;

    for (i=0; i<n; i++) {
        u = y.val[i];
        for (j=0; j<k; j++) {
            u -= gretl_matrix_get(X, i, j) * b.val[j];
        }
        s2 += u * u;
    }

    s2 /= (n - k + nr); /* nr = number of restrictions */

    return s2;
}

private int get_ols_vcv(GretlMatrixRaw * V, double * s2)
{
    if (gretl_invert_general_matrix(V)) {
        gretl_matrix_print(V, "get_ols_vcv: inversion failed");
        return 1;
    }

    if (s2 != null) {
        gretl_matrix_multiply_by_scalar(V, *s2);
    }

    return 0;
}

private void
get_ols_uhat(const(GretlVector) * y, const(GretlMatrixRaw) * X,
             const(GretlVector) * b, GretlVector * uhat)
{
    int ncoeff = gretl_vector_get_length(b);
    int n = gretl_vector_get_length(uhat);
    int i, j;
    double uh;

    for (i=0; i<n; i++) {
        uh = y.val[i];
        for (j=0; j<ncoeff; j++) {
            uh -= b.val[j] * gretl_matrix_get(X, i, j);
        }
        uhat.val[i] = uh;
    }
}

private int svd_ols_work(GretlMatrixRaw * A,
                         GretlMatrixRaw * B,
                         double * s,
                         int use_dc)
{
    double * work = null;
    double rcond = 0.0;
    integer m, n, nrhs;
    integer lda, ldb;
    integer lwork = -1;
    integer liwork = 0;
    integer rank;
    integer info;
    integer * iwork = null;
    int err = 0;

    work = cast(double *) lapack_malloc((*work).sizeof);
    if (work == null) {
        return GretlError.alloc;
    }

    lda = ldb = m = A.rows;
    n = A.cols;
    nrhs = B.cols;

    /* workspace query */
    if (use_dc) {
        dgelsd_(&m, &n, &nrhs, A.val, &lda, B.val, &ldb, s, &rcond,
                &rank, work, &lwork, &liwork, &info);
    } else {
        dgelss_(&m, &n, &nrhs, A.val, &lda, B.val, &ldb, s, &rcond,
                &rank, work, &lwork, &info);
    }

    if (info != 0 || work[0] <= 0.0) {
        return wspace_fail(info, work[0]);
    }

    lwork = cast(integer) work[0];
    work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        return GretlError.alloc;
    }

    if (use_dc) {
        iwork = cast(integer *) malloc(liwork * (*iwork).sizeof);
        if (iwork == null) {
            return GretlError.alloc;
        }
    }

    /* get actual solution */
    if (use_dc) {
        dgelsd_(&m, &n, &nrhs, A.val, &lda, B.val, &ldb, s, &rcond,
                &rank, work, &lwork, iwork, &info);
    } else {
        dgelss_(&m, &n, &nrhs, A.val, &lda, B.val, &ldb, s, &rcond,
                &rank, work, &lwork, &info);
    }

    if (info != 0) {
        fprintf(stderr, "svd_ols_work: got info = %d (with use_dc = %d)\n",
                info, use_dc);
        err = GretlError.noconv;
    } else if (rank < n) {
        fprintf(stderr, "svd_ols_work:\n"
                ~ " data matrix X (%d x %d) has column rank %d\n",
                m, n, cast(int) rank);
    }

    lapack_free(work);
    free(iwork);

    return err;
}

/*
 * gretl_matrix_SVD_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or null if this is not needed.
 * @uhat: vector to hold the regression residuals, or null if
 * these are not needed.
 * @s2: pointer to receive residual variance, or null.  Note:
 * if @s2 is null, the vcv estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using SVD decomposition, and puts the
 * coefficient estimates in @b.  Optionally, calculates the
 * covariance matrix in @vcv and the residuals in @uhat.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_SVD_ols(const(GretlVector) * y, const(GretlMatrixRaw) * X,
                         GretlVector * b, GretlMatrixRaw * vcv,
                         GretlVector * uhat, double * s2)
{
    GretlVector * A = null;
    GretlMatrixRaw * B = null;
    double * s = null;
    int k, use_dc = 0;
    int err = 0;

    if (gretl_is_null_matrix(y) ||
        gretl_is_null_matrix(X) ||
        gretl_is_null_matrix(b)) {
        return GretlError.data;
    }

    /* (the PREFER_DGELSD variant, 0 in the C source, is not ported) */

    k = X.cols;

    if (gretl_vector_get_length(b) != k) {
        return GretlError.nonconf;
    }

    A = gretl_matrix_copy_tmp(X);
    B = gretl_matrix_copy_tmp(y);

    if (A == null || B == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* for singular values of A */
    s = cast(double *) malloc(k * (*s).sizeof);
    if (s == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = svd_ols_work(A, B, s, use_dc);

    if (!err) {
        int i;

        for (i=0; i<k; i++) {
            b.val[i] = B.val[i];
        }
        if (vcv != null) {
            err = get_SVD_ols_vcv(A, B, s, vcv, s2);
        }
        if (uhat != null) {
            get_ols_uhat(y, X, b, uhat);
        }
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(B);
    free(s);

    return err;
}

/*
 * gretl_matrix_multi_SVD_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates, or null.
 * @E: T x g matrix to hold the regression residuals, or null if these are
 * not needed.
 * @XTXi: location to receive (X'X)^{-1}, or null if this is not needed.
 *
 * Computes OLS estimates using SVD decomposition, and puts the
 * coefficient estimates in @B.  Optionally, calculates the
 * residuals in @E, (X'X)^{-1} in @XTXi.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_multi_SVD_ols(const(GretlMatrixRaw) * Y,
                               const(GretlMatrixRaw) * X,
                               GretlMatrixRaw * B,
                               GretlMatrixRaw * E,
                               GretlMatrixRaw ** XTXi)
{
    int g, k, T;
    GretlMatrixRaw * A = null;
    GretlMatrixRaw * C = null;
    double * s = null;
    int free_B = 0;
    int use_dc = 0;
    int err = 0;

    if (gretl_is_null_matrix(Y) ||
        gretl_is_null_matrix(X)) {
        return GretlError.data;
    }

    /* (the PREFER_DGELSD variant, 0 in the C source, is not ported) */

    g = Y.cols;
    k = X.cols;
    T = X.rows;

    if (B == null) {
        B = gretl_matrix_alloc(k, g);
        if (B == null) {
            return GretlError.alloc;
        }
        free_B = 1;
    }

    if (B.rows != k || B.cols != g) {
        err = GretlError.nonconf;
    } else if (Y.rows != T) {
        err = GretlError.nonconf;
    } else if (E != null && (E.cols != g || E.rows != T)) {
        err = GretlError.nonconf;
    } else if (k > T) {
        err = GretlError.df;
    }

    A = gretl_matrix_copy_tmp(X);
    C = gretl_matrix_copy_tmp(Y);

    if (A == null || C == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* for singular values of A */
    s = cast(double *) malloc(k * (*s).sizeof);
    if (s == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = svd_ols_work(A, C, s, use_dc);

    if (!err) {
        /* coeffs: extract the first k rows from @C */
        double bij;
        int i, j;

        for (i=0; i<k; i++) {
            for (j=0; j<g; j++) {
                bij = gretl_matrix_get(C, i, j);
                gretl_matrix_set(B, i, j, bij);
            }
        }
    }

    if (!err && E != null) {
        /* compute residuals, if wanted */
        int i, imax = E.rows * E.cols;

        gretl_matrix_multiply(X, B, E);
        for (i=0; i<imax; i++) {
            E.val[i] = Y.val[i] - E.val[i];
        }
    }

    if (!err && XTXi != null) {
        /* build (X'X)^{-1}, if wanted */
        *XTXi = gretl_matrix_alloc(k, k);
        if (*XTXi == null) {
            err = GretlError.alloc;
        } else {
            err = get_SVD_ols_vcv(A, C, s, *XTXi, null);
        }
    }

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(C);
    free(s);

    if (free_B) {
        gretl_matrix_free(B);
    }

    return err;
}

/*
 * gretl_matrix_moore_penrose:
 * @a: m x n matrix.
 *
 * Computes the generalized inverse of matrix @a via its SVD
 * factorization, with the help of the lapack function
 * dgesvd.  On exit the original matrix is overwritten by
 * the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_moore_penrose(GretlMatrixRaw * A, double tol)
{
    GretlMatrixRaw * U = null;
    GretlMatrixRaw * S = null;
    GretlMatrixRaw * VT = null;
    int err = 0;

    if (gretl_is_null_matrix(A)) {
        return GretlError.data;
    }

    err = real_gretl_matrix_SVD(A, &U, &S, &VT, 0);
    tol = svd_smin(A, S.val[0], tol);

    if (!err) {
        GretlMatrixRaw * Vsel = null;
        int nsv = min(A.rows, A.cols);
        int i, j, k = 0;
        double x;

        for (i=0; i<nsv; i++) {
            if (S.val[i] > tol) {
                k++;
            }
        }

        if (k < VT.rows) {
            Vsel = gretl_matrix_alloc(k, VT.cols);
            if (Vsel == null) {
                err = GretlError.alloc;
                goto bailout;
            }
            for (j=0; j<VT.cols; j++) {
                for (i=0; i<k; i++) {
                    x = gretl_matrix_get(VT, i, j);
                    gretl_matrix_set(Vsel, i, j, x);
                }
            }
        }

        /* U <- U .* S^{-1}, for S[j] > min */
        for (i=0; i<U.rows; i++) {
            for (j=0; j<k; j++) {
                x = gretl_matrix_get(U, i, j);
                gretl_matrix_set(U, i, j, x / S.val[j]);
            }
        }
        if (k < U.cols) {
            gretl_matrix_reuse(U, -1, k);
        }

        err = gretl_matrix_multiply_mod(U, GretlMatrixMod.none,
                                        Vsel != null ? Vsel : VT,
                                        GretlMatrixMod.none,
                                        A, GretlMatrixMod.none);
        if (!err) {
            gretl_matrix_transpose_in_place(A);
        }
        gretl_matrix_free(Vsel);
    }

 bailout:

    gretl_matrix_free(U);
    gretl_matrix_free(S);
    gretl_matrix_free(VT);

    return err;
}

/*
 * gretl_SVD_invert_matrix:
 * @a: n x n matrix to invert.
 *
 * Computes the inverse (or generalized inverse) of a general square
 * matrix using SVD factorization, with the help of the lapack function
 * dgesvd.  If any of the singular values of @a are less than 1.0e-9
 * the Moore-Penrose generalized inverse is computed instead of the
 * standard inverse.  On exit the original matrix is overwritten by
 * the inverse.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_SVD_invert_matrix(GretlMatrixRaw * a)
{
    GretlMatrixRaw * u = null;
    GretlMatrixRaw * s = null;
    GretlMatrixRaw * vt = null;
    double x;
    int i, j, n;
    int rank = 0;
    int err = 0;

    if (gretl_is_null_matrix(a)) {
        return GretlError.data;
    }

    if (a.rows != a.cols) {
        err = GretlError.nonconf;
        goto bailout;
    }

    n = a.rows;

    /* a = USV' ; a^{-1} = VWU' where W holds inverse of diag elements of S */

    err = real_gretl_matrix_SVD(a, &u, &s, &vt, 0);

    if (!err) {
        double smin = svd_smin(a, s.val[0], NADBL);

        for (i=0; i<n; i++) {
            if (s.val[i] > smin) {
                rank++;
            } else {
                break;
            }
        }

        if (rank < n) {
            GretlMatrixRaw * vt2;

            fprintf(stderr, "gretl_SVD_invert_matrix: rank = %d (dim = %d)\n",
                    rank, n);
            fputs("Warning: computing Moore-Penrose generalized inverse\n", stderr);

            vt2 = gretl_matrix_alloc(rank, n);
            if (vt2 == null) {
                err = GretlError.alloc;
                goto bailout;
            }
            for (i=0; i<rank; i++) {
                for (j=0; j<n; j++) {
                    x = gretl_matrix_get(vt, i, j);
                    gretl_matrix_set(vt2, i, j, x);
                }
            }
            gretl_matrix_free(vt);
            vt = vt2;
            gretl_matrix_reuse(u, n, rank);
        }
    }

    if (!err) {
        /* invert singular values */
        for (j=0; j<rank; j++) {
            for (i=0; i<n; i++) {
                x = gretl_matrix_get(u, i, j);
                gretl_matrix_set(u, i, j, x / s.val[j]);
            }
        }
        err = gretl_matrix_multiply_mod(vt, GretlMatrixMod.transpose,
                                        u, GretlMatrixMod.transpose,
                                        a, GretlMatrixMod.none);
    }

 bailout:

    gretl_matrix_free(u);
    gretl_matrix_free(s);
    gretl_matrix_free(vt);

    return err;
}

/*
 * gretl_matrix_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @b: vector to hold coefficient estimates, or null if this
 * is not needed.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or null if this is not needed.
 * @uhat: vector to hold the regression residuals, or null if
 * these are not needed.
 * @s2: pointer to receive residual variance, or null.  Note:
 * if @s2 is null, the "vcv" estimate will be plain (X'X)^{-1}.
 *
 * Computes OLS estimates using Cholesky factorization by default,
 * but with a fallback to QR decomposition if the data are highly
 * ill-conditioned, and puts the coefficient estimates in @b.
 * Optionally, calculates the covariance matrix in @vcv and the
 * residuals in @uhat.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_ols(const(GretlVector) * y,
                     const(GretlMatrixRaw) * X,
                     GretlVector * b,
                     GretlMatrixRaw * vcv,
                     GretlVector * uhat,
                     double * s2)
{
    GretlMatrixRaw * XTX = null;
    GretlMatrixRaw * c = null;
    int use_lapack = 0;
    int try_QR = 0;
    int nasty = 0;
    int k, T, err = 0;

    if (gretl_is_null_matrix(y) ||
        gretl_is_null_matrix(X)) {
        return GretlError.data;
    }

    k = X.cols;
    T = X.rows;

    /* check dimensions */
    if ((b != null && gretl_vector_get_length(b) != k) ||
        gretl_vector_get_length(y) != T) {
        return GretlError.nonconf;
    } else if (T < k) {
        return GretlError.df;
    } else if (vcv != null && (vcv.rows != k || vcv.cols != k)) {
        return GretlError.nonconf;
    }

    if (b == null) {
        /* we'll handle the case where @b is not wanted */
        c = gretl_matrix_alloc(X.cols, 1);
        b = c;
    }

    if (libset_get_bool(SetKey.use_svd)) {
        /* this call is deferred in case @b is null on input */
        err = gretl_matrix_SVD_ols(y, X, b, vcv, uhat, s2);
        goto finish;
    }

    /* It should be worth using lapack's Cholesky routine if the input
       is big enough, but what's the best definition of "big enough"?
    */
    if (k >= 50 || (T >= 250 && k >= 30)) {
        use_lapack = 1;
        XTX = gretl_matrix_XTX_new(X);
    } else {
        /* using gretl's native Cholesky */
        XTX = gretl_matrix_packed_XTX_new(X, &nasty);
    }

    if (XTX == null) {
        err = GretlError.alloc;
        goto finish;
    }

    if (use_lapack || !nasty) {
        /* preliminary Cholesky step: shouldn't fail */
        gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                                  y, GretlMatrixMod.none,
                                  b, GretlMatrixMod.none);
    }

    if (use_lapack) {
        if (vcv != null) {
            /* we'll want this even if we switch to QR */
            gretl_matrix_copy_values(vcv, XTX);
        }
        err = gretl_cholesky_decomp_solve(XTX, b);
        if (err) {
            try_QR = 1;
        }
    } else {
        if (vcv != null) {
            /* we'll want this even if the next step fails */
            gretl_matrix_unvectorize_h(vcv, XTX);
        }
        if (!nasty) {
            err = native_cholesky_decomp_solve(XTX, b);
        }
        if (nasty || err == GretlError.singular) {
            try_QR = 1;
        }
    }

    gretl_matrix_free(XTX);

    if (try_QR) {
        fprintf(stderr, "gretl_matrix_ols: switching to QR decomp\n");
        err = gretl_matrix_QR_ols(y, X, b, null, null, null);
    }

    if (!err) {
        if (s2 != null) {
            *s2 = get_ols_error_variance(y, X, b, 0);
        }
        if (vcv != null) {
            err = get_ols_vcv(vcv, s2);
        }
        if (uhat != null) {
            get_ols_uhat(y, X, b, uhat);
        }
    }

 finish:

    gretl_matrix_free(c);

    return err;
}

/*
 * gretl_matrix_multi_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates, or null.
 * @E: T x g matrix to hold the regression residuals, or null if these are
 * not needed.
 * @XTXi: location to receive (X'X)^{-1} on output, or null if this is not
 * needed.
 *
 * Computes OLS estimates using Cholesky decomposition by default, but
 * with a fallback to QR decomposition if the data are highly
 * ill-conditioned, and puts the coefficient estimates in @B.
 * Optionally, calculates the residuals in @E.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_multi_ols(const(GretlMatrixRaw) * Y,
                           const(GretlMatrixRaw) * X,
                           GretlMatrixRaw * B,
                           GretlMatrixRaw * E,
                           GretlMatrixRaw ** XTXi)
{
    GretlMatrixRaw * XTX = null;
    int g, T, k;
    int free_B = 0;
    int nasty = 0;
    int err = 0;

    if (libset_get_bool(SetKey.use_svd)) {
        return gretl_matrix_multi_SVD_ols(Y, X, B, E, XTXi);
    }

    if (gretl_is_null_matrix(Y) ||
        gretl_is_null_matrix(X)) {
        return GretlError.invarg;
    }

    g = Y.cols;
    T = X.rows;
    k = X.cols;

    if (B == null) {
        /* create a throw-away B */
        B = gretl_matrix_alloc(k, g);
        if (B == null) {
            return GretlError.alloc;
        }
        free_B = 1;
    }

    if (B.rows != k || B.cols != g) {
        fprintf(stderr, "gretl_matrix_multi_ols: B is %d x %d, should be %d x %d\n",
                B.rows, B.cols, k, g);
        err = GretlError.nonconf;
    } else if (Y.rows != T) {
        fprintf(stderr, "gretl_matrix_multi_ols: Y has %d rows, should have %d\n",
                Y.rows, T);
        err = GretlError.nonconf;
    } else if (E != null && (E.rows != T || E.cols != g)) {
        fprintf(stderr, "gretl_matrix_multi_ols: E is %d x %d, should be %d x %d\n",
                E.rows, E.cols, T, g);
        err = GretlError.nonconf;
    } else if (k > T) {
        err = GretlError.df;
    }

    if (!err) {
        XTX = gretl_matrix_XTX_new(X);
        if (XTX == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                                        Y, GretlMatrixMod.none,
                                        B, GretlMatrixMod.none);
    }

    if (!err) {
        err = nasty = gretl_cholesky_decomp_solve(XTX, B);
        if (err == GretlError.singular) {
            fprintf(stderr, "gretl_matrix_multi_ols: switching to QR decomp\n");
            err = gretl_matrix_QR_ols(Y, X, B, E, XTXi, null);
        }
    }

    if (!err && !nasty && E != null) {
        gretl_matrix_copy_values(E, Y);
        gretl_matrix_multiply_mod(X, GretlMatrixMod.none,
                                  B, GretlMatrixMod.none,
                                  E, GretlMatrixMod.decrement);
    }

    if (!err && !nasty && XTXi != null) {
        integer info = 0, ik = k;
        char uplo = 'L';

        dpotri_(&uplo, &ik, XTX.val, &ik, &info);
        gretl_matrix_mirror(XTX, uplo);
        *XTXi = XTX;
    } else {
        gretl_matrix_free(XTX);
    }

    if (free_B) {
        gretl_matrix_free(B);
    }

    return err;
}

/* helper function for gretl_matrix_factorized_ols() */

private int * get_factor_sorted_data(const(GretlMatrixRaw) * Y,
                                     const(GretlMatrixRaw) * X,
                                     const(GretlVector) * f,
                                     GretlMatrixRaw * fY,
                                     GretlMatrixRaw * fX,
                                     GretlMatrixRaw * sf,
                                     int * err)
{
    int * ret = null;
    struct RSort {
        double x;
        int row;
    }
    RSort * rs;
    double d;
    int T = Y.rows;
    int g = Y.cols;
    int k = X.cols;
    int i, t;

    rs = cast(RSort *) malloc(T * (*rs).sizeof);
    if (rs == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (t=0; t<T; t++) {
        rs[t].x = f.val[t];
        rs[t].row = t;
    }

    qsort(rs, T, (*rs).sizeof, &compare_values);

    ret = cast(int *) malloc(T * (*ret).sizeof);

    /* transcribe Y and X in sort order: put t in the inner loop to
       mitigate cache misses
    */
    for (i=0; i<g; i++) {
        for (t=0; t<T; t++) {
            d = gretl_matrix_get(Y, rs[t].row, i);
            gretl_matrix_set(fY, t, i, d);
        }
    }
    for (i=0; i<k; i++) {
        for (t=0; t<T; t++) {
            d = gretl_matrix_get(X, rs[t].row, i);
            gretl_matrix_set(fX, t, i, d);
        }
    }
    for (t=0; t<T; t++) {
        sf.val[t] = rs[t].x;
        ret[t] = rs[t].row;
    }

    free(rs);

    return ret;
}

private int vector_is_sorted(const(GretlVector) * f, int T)
{
    int t;

    for (t=1; t<T; t++) {
        if (f.val[t] < f.val[t-1]) {
            return 0;
        }
    }

    return 1;
}

private int reorder_rows(GretlMatrixRaw * m, int * r)
{
    GretlMatrixRaw * tmp = gretl_matrix_copy(m);
    double mti;
    int T = m.rows;
    int k = m.cols;
    int i, t;

    for (t=0; t<T; t++) {
        for (i=0; i<k; i++) {
            mti = gretl_matrix_get(tmp, t, i);
            gretl_matrix_set(m, r[t], i, mti);
        }
    }

    gretl_matrix_free(tmp);

    return 0;
}

/*
 * gretl_matrix_factorized_ols:
 * @Y: T x g matrix of dependent variable vector.
 * @X: T x k matrix of independent variables.
 * @fac: discrete factor vector.
 * @B: matrix to hold coefficient estimates, or null if this
 * is not needed.
 * @V: matrix to hold the covariance matrix of the coefficients,
 * or null if this is not needed.
 * @U: matrix to hold the regression residuals, or null if
 * these are not needed.
 * @A: matrix to hold the "fixed effects" (per factor-value constants),
 * or null if these are not needed.
 *
 * Works like gretl_matrix_ols() except that prior to the OLS
 * computation both @Y and @X are demeaned, using group means
 * based on the values of @fac. This is equivalent to regressing
 * @Y on @X along with a set of dummy variables coding for the
 * value of @fac, but potentially much more efficient if @fac
 * has many distinct values.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_factorized_ols(const(GretlMatrixRaw) * Y,
                                const(GretlMatrixRaw) * X,
                                const(GretlVector) * fac,
                                GretlMatrixRaw * B,
                                GretlMatrixRaw * V,
                                GretlVector * U,
                                GretlVector * A)
{
    GretlVector * fY = null;
    GretlMatrixRaw * fX = null;
    GretlMatrixRaw * sf = null;
    int * order = null;
    double * mean = null;
    double * ymean;
    double * xmean;
    double f, s2;
    double yti, xti;
    int T = X.rows;
    int k = X.cols;
    int g = Y.cols;
    int nm = g + k;
    int nfvals = 0;
    int do_sort = 0;
    int i, n, t, t0, s;
    int err = 0;

    if (Y.rows != T || fac.rows != T) {
        return GretlError.invarg;
    }

    for (t=0; t<T; t++) {
        f = fac.val[t];
        if (na(f) || f != floor(f)) {
            gretl_errmsg_set(_("The factor variable must be integer-valued").ptr);
            return GretlError.invarg;
        }
    }

    if (vector_is_sorted(fac, T)) {
        /* storage for de-meaned Y and X */
        fY = gretl_matrix_copy(Y);
        fX = gretl_matrix_copy(X);
        sf = cast(GretlMatrixRaw *) fac;
    } else {
        /* storage for fac-sorted and de-meaned Y and X, plus sorted fac */
        fY = gretl_matrix_alloc(T, g);
        fX = gretl_matrix_alloc(T, k);
        sf = gretl_matrix_alloc(T, 1);
        do_sort = 1;
    }

    if (fY == null || fX == null || sf == null) {
        err = GretlError.alloc;
    } else if (do_sort) {
        order = get_factor_sorted_data(Y, X, fac, fY, fX, sf, &err);
    }

    if (!err) {
        mean = cast(double *) malloc(nm * (*mean).sizeof);
        if (mean == null) {
            err = GretlError.alloc;
        } else {
            ymean = mean;
            xmean = ymean + g;
            for (i=0; i<nm; i++) {
                mean[i] = 0.0;
            }
        }
    }

    if (err) {
        goto bailout;
    }

    nfvals = t0 = 0;
    for (t=0; t<T; t++) {
        for (i=0; i<g; i++) {
            yti = gretl_matrix_get(fY, t, i);
            if (na(yti)) {
                err = GretlError.missdata;
                break;
            }
            ymean[i] += yti;
        }
        for (i=0; i<k; i++) {
            xti = gretl_matrix_get(fX, t, i);
            if (na(xti)) {
                err = GretlError.missdata;
                break;
            }
            xmean[i] += xti;
        }
        if (err) {
            break;
        }
        if ((t < T-1 && sf.val[t+1] != sf.val[t]) || t == T-1) {
            /* finish the current factor value */
            nfvals++;
            n = t - t0 + 1;
            for (i=0; i<g; i++) {
                ymean[i] /= n;
                for (s=t0; s<=t; s++) {
                    yti = gretl_matrix_get(fY, s, i);
                    gretl_matrix_set(fY, s, i, yti - ymean[i]);
                }
                ymean[i] = 0.0;
            } for (i=0; i<k; i++) {
                xmean[i] /= n;
                for (s=t0; s<=t; s++) {
                    xti = gretl_matrix_get(fX, s, i);
                    gretl_matrix_set(fX, s, i, xti - xmean[i]);
                }
                xmean[i] = 0.0;
            }
            t0 = t + 1;
        }
    }

    if (err) {
        goto bailout;
    }

    if (g == 1) {
        double * ps2 = V != null ? &s2 : null;

        err = gretl_matrix_ols(fY, fX, B, V, U, ps2);
    } else if (V == null) {
        err = gretl_matrix_multi_ols(fY, fX, B, U, null);
    } else {
        GretlMatrixRaw * Vtmp = null;

        err = gretl_matrix_multi_ols(fY, fX, B, U, &Vtmp);
        if (!err && Vtmp != null) {
            gretl_matrix_replace_content(V, Vtmp);
            gretl_matrix_free(Vtmp);
        }
    }

    if (U != null && order != null) {
        reorder_rows(U, order);
    }

    if (!err && A != null) {
        /* compute "fixed-effects" */
        double aij, aijt, xth;
        int j, h;

        err = gretl_matrix_realloc(A, nfvals, g);
        if (err) {
            goto bailout;
        }
        gretl_matrix_zero(A);

        i = t0 = 0;
        for (t=0; t<T; t++) {
            s = do_sort ? order[t] : t;
            for (j=0; j<g; j++) {
                aijt = gretl_matrix_get(Y, s, j);
                for (h=0; h<k; h++) {
                    xth = gretl_matrix_get(X, s, h);
                    aijt -= xth * gretl_matrix_get(B, h, j);
                }
                /* cumulate aij */
                aij = gretl_matrix_get(A, i, j);
                gretl_matrix_set(A, i, j, aij + aijt);
            }
            if ((t < T-1 && sf.val[t+1] != sf.val[t]) || t == T-1) {
                /* finish the current factor value */
                n = t - t0 + 1;
                for (j=0; j<g; j++) {
                    aij = gretl_matrix_get(A, i, j);
                    gretl_matrix_set(A, i, j, aij / n);
                }
                t0 = t + 1;
                i++;
            }
        }
    }

 bailout:

    gretl_matrix_free(fY);
    gretl_matrix_free(fX);
    free(mean);
    free(order);
    if (sf != fac) {
        gretl_matrix_free(sf);
    }

    return err;
}

/* Construct W = X'X augmented by R and R': we need this if we're
   calculating the covariance matrix for restricted least squares.
*/

private GretlMatrixRaw * build_augmented_XTX(const(GretlMatrixRaw) * X,
                                          const(GretlMatrixRaw) * R,
                                          int * err)
{
    GretlMatrixRaw * XTX, W;
    int k = X.cols;
    int kW = k + R.rows;

    XTX = gretl_matrix_XTX_new(X);
    W = gretl_zero_matrix_new(kW, kW);

    if (XTX == null || W == null) {
        gretl_matrix_free(XTX);
        gretl_matrix_free(W);
        *err = GretlError.alloc;
        return null;
    }

    if (!*err) {
        double x;
        int i, j;

        for (i=0; i<k; i++) {
            for (j=0; j<k; j++) {
                x = gretl_matrix_get(XTX, i, j);
                gretl_matrix_set(W, i, j, x);
            }
        }
        for (i=0; i<R.rows; i++) {
            for (j=0; j<R.cols; j++) {
                x = gretl_matrix_get(R, i, j);
                gretl_matrix_set(W, i+k, j, x);
                gretl_matrix_set(W, j, i+k, x);
            }
        }
    }

    gretl_matrix_free(XTX);

    return W;
}

/* constrained least squares using lapack:

    minimize || y - X*b ||_2  subject to R*b = q
*/

private int gretl_matrix_gglse(const(GretlVector) * y,
                               const(GretlMatrixRaw) * X,
                               const(GretlMatrixRaw) * R,
                               const(GretlVector) * q,
                               GretlVector * b)
{
    GretlMatrixRaw * A, B, c, d;
    integer info;
    integer m = X.rows;
    integer n = X.cols;
    integer p = R.rows;
    integer lwork = -1;
    double * work;
    int err = 0;

    /* note: all the input matrices get overwritten */
    A = gretl_matrix_copy(X);
    B = gretl_matrix_copy(R);
    c = gretl_matrix_copy(y);

    if (q != null) {
        d = gretl_matrix_copy(q);
    } else {
        d = gretl_zero_matrix_new(p, 1);
    }

    work = cast(double *) lapack_malloc((*work).sizeof);

    if (A == null || B == null || c == null ||
        d == null || work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* determine optimal workspace */
    dgglse_(&m, &n, &p, A.val, &m, B.val, &p, c.val,
            d.val, b.val, work, &lwork, &info);

    if (info != 0) {
        err = wspace_fail(info, work[0]);
    } else {
        lwork = cast(integer) work[0];
        work = cast(double *) lapack_realloc(work, lwork * (*work).sizeof);
        if (work == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        /* do constrained calculation */
        dgglse_(&m, &n, &p, A.val, &m, B.val, &p, c.val,
                d.val, b.val, work, &lwork, &info);
        if (info != 0) {
            fprintf(stderr, "dgglse gave info = %d\n", cast(int) info);
            err = (info < 0) ? GretlError.data : GretlError.singular;
        }
    }

    lapack_free(work);

 bailout:

    gretl_matrix_free(A);
    gretl_matrix_free(B);
    gretl_matrix_free(c);
    gretl_matrix_free(d);

    return err;
}

/* get_exact_list: constructs a list of the parameter positions
   (columns of the @R matrix) corresponding to unitary restrictions --
   that is, restrictions that stipulate a definite numerical value for
   a given parameter. We want this so that we can set the variance of
   such parameters (and also their covariances with other parameter
   estimates) to exactly zero.
*/

private int * get_exact_list(const(GretlMatrixRaw) * R)
{
    int * list = null;
    int n, n_exact = 0;
    int i, j;

    for (i=0; i<R.rows; i++) {
        n = 0;
        for (j=0; j<R.cols && n<2; j++) {
            if (gretl_matrix_get(R, i, j) != 0.0) {
                n++;
            }
        }
        if (n == 1) {
            n_exact++;
        }
    }

    if (n_exact > 0) {
        list = gretl_list_new(n_exact);
    }

    if (list != null) {
        int col = 0, k = 1;

        for (i=0; i<R.rows && k<=n_exact; i++) {
            n = 0;
            for (j=0; j<R.cols && n<2; j++) {
                if (gretl_matrix_get(R, i, j) != 0.0) {
                    col = j;
                    n++;
                }
            }
            if (n == 1) {
                list[k++] = col;
            }
        }
    }

    return list;
}

/*
 * gretl_matrix_restricted_ols:
 * @y: dependent variable vector.
 * @X: matrix of independent variables.
 * @R: left-hand restriction matrix, as in Rb = q.
 * @q: right-hand restriction vector or null.
 * @b: vector to hold coefficient estimates.
 * @vcv: matrix to hold the covariance matrix of the coefficients,
 * or null if this is not needed.
 * @uhat: vector to hold residuals, if wanted.
 * @s2: pointer to receive residual variance, or null.  If vcv is non-null
 * and s2 is null, the vcv estimate is just W^{-1}.
 *
 * Computes OLS estimates restricted by R and q, using the lapack
 * function dgglse(), and puts the coefficient estimates in @b.
 * Optionally, calculates the covariance matrix in @vcv. If @q is
 * null this is taken as equivalent to a zero vector.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int
gretl_matrix_restricted_ols(const(GretlVector) * y, const(GretlMatrixRaw) * X,
                            const(GretlMatrixRaw) * R, const(GretlVector) * q,
                            GretlVector * b, GretlMatrixRaw * vcv,
                            GretlVector * uhat, double * s2)
{
    GretlMatrixRaw * W = null;
    double x;
    int k = X.cols;
    int nr = R.rows;
    int i, j, err = 0;

    if (gretl_vector_get_length(b) != k) {
        fprintf(stderr, "gretl_matrix_restricted_ols: "
                ~ "b should be a %d-vector\n", k);
        err = GretlError.nonconf;
    }

    if (!err && vcv != null) {
        W = build_augmented_XTX(X, R, &err);
    }

    if (!err) {
        err = gretl_matrix_gglse(y, X, R, q, b);
    }

    if (!err) {
        if (s2 != null) {
            *s2 = get_ols_error_variance(y, X, b, nr);
        }

        if (W != null) {
            int * exlist = null;

            err = get_ols_vcv(W, s2);

            if (!err) {
                for (i=0; i<k; i++) {
                    for (j=0; j<k; j++) {
                        x = gretl_matrix_get(W, i, j);
                        gretl_matrix_set(vcv, i, j, x);
                    }
                }
                exlist = get_exact_list(R);
            }

            if (exlist != null) {
                int p;

                for (p=0; p<k; p++) {
                    if (in_gretl_list(exlist, p)) {
                        for (i=0; i<k; i++) {
                            gretl_matrix_set(vcv, i, p, 0.0);
                        }
                        for (j=0; j<k; j++) {
                            gretl_matrix_set(vcv, p, j, 0.0);
                        }
                    }
                }
                free(exlist);
            }
        }

        if (uhat != null) {
            get_ols_uhat(y, X, b, uhat);
        }
    }

    if (W != null) {
        gretl_matrix_free(W);
    }

    return err;
}

/*
 * gretl_matrix_restricted_multi_ols:
 * @Y: dependent variable matrix, T x g.
 * @X: matrix of independent variables, T x k.
 * @R: left-hand restriction matrix, as in RB = q.
 * @q: right-hand restriction matrix.
 * @B: matrix to hold coefficient estimates, k x g.
 * @U: matrix to hold residuals (T x g), if wanted.
 * @pW: pointer to matrix to hold the RLS counterpart to (X'X)^{-1},
 * if wanted.
 *
 * Computes LS estimates restricted by @R and @q, putting the
 * coefficient estimates into @B. The @R matrix must have
 * g*k columns; each row represents a linear restriction;
 * @q must be a column vector of length equal to the
 * number of rows of @R.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int
gretl_matrix_restricted_multi_ols(const(GretlMatrixRaw) * Y,
                                  const(GretlMatrixRaw) * X,
                                  const(GretlMatrixRaw) * R,
                                  const(GretlMatrixRaw) * q,
                                  GretlMatrixRaw * B,
                                  GretlMatrixRaw * U,
                                  GretlMatrixRaw ** pW)
{
    GretlMatrixBlock * M;
    GretlMatrixRaw * XTX, RXR, XYq;
    GretlMatrixRaw * Yi, XYi;
    GretlMatrixRaw * V = null;
    int T = Y.rows;  /* sample length */
    int k = X.cols;  /* number of regressors */
    int g = Y.cols;  /* number of dependent vars */
    int nc = g * k;  /* total coefficients */
    int nr = R.rows; /* number of restrictions */
    int p = nc + nr; /* coeffs plus restrictions */
    int dsize, offset;
    int i, r, err = 0;

    if (X.rows != T) {
        return GretlError.nonconf;
    } else if (B.rows != k || B.cols != g) {
        return GretlError.nonconf;
    } else if (R.cols != nc || q.rows != nr || q.cols != 1) {
        return GretlError.nonconf;
    } else if (U != null && (U.rows != T || U.cols != g)) {
        return GretlError.nonconf;
    }

    M = gretl_matrix_block_new(&XTX, k, k,
                               &RXR, p, p,
                               &XYq, p, 1,
                               &Yi,  T, 1,
                               &XYi, k, 1,
                               cast(GretlMatrixRaw **) null);
    if (M == null) {
        return GretlError.alloc;
    }

    gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                              X, GretlMatrixMod.none,
                              XTX, GretlMatrixMod.none);

    gretl_matrix_zero(RXR);

    dsize = cast(int) (T * double.sizeof);
    offset = r = 0;

    /* Form the "big" X'X and X'y matrices, bordered by the
       restriction:

       "RXR" = [(I_g ** X'X) ~ R'] | (R ~ 0)
       "XYq" = vec(X'Y) | q
    */

    for (i=0; i<g; i++) {
        gretl_matrix_inscribe_matrix(RXR, XTX, r, r,
                                     GretlMatrixMod.none);
        memcpy(Yi.val, Y.val + offset, dsize);
        gretl_matrix_multiply_mod(X, GretlMatrixMod.transpose,
                                  Yi, GretlMatrixMod.none,
                                  XYi, GretlMatrixMod.none);
        gretl_matrix_inscribe_matrix(XYq, XYi, r, 0,
                                     GretlMatrixMod.none);
        r += k;
        offset += T;
    }

    gretl_matrix_inscribe_matrix(RXR, R, r, 0, GretlMatrixMod.none);
    gretl_matrix_inscribe_matrix(RXR, R, 0, nc, GretlMatrixMod.transpose);
    gretl_matrix_inscribe_matrix(XYq, q, nc, 0, GretlMatrixMod.none);

    if (pW != null) {
        /* keep a copy of @V for inversion */
        V = gretl_matrix_copy(RXR);
        if (V == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        /* solve for stacked coeff vector in XYq */
        err = gretl_LU_solve(RXR, XYq);
        if (!err) {
            /* transcribe to B */
            dsize = cast(int) (nc * double.sizeof);
            memcpy(B.val, XYq.val, dsize);
        }
    }

    if (!err && U != null) {
        /* compute residuals */
        gretl_matrix_copy_values(U, Y);
        gretl_matrix_multiply_mod(X, GretlMatrixMod.none,
                                  B, GretlMatrixMod.none,
                                  U, GretlMatrixMod.decrement);
    }

    if (!err && pW != null) {
        /* compute variance-related matrix */
        err = gretl_invert_general_matrix(V);
        if (!err) {
            *pW = gretl_matrix_alloc(nc, nc);
            if (*pW == null) {
                err = GretlError.alloc;
            } else {
                double wij;
                int j;

                for (j=0; j<nc; j++) {
                    for (i=0; i<nc; i++) {
                        wij = gretl_matrix_get(V, i, j);
                        gretl_matrix_set(*pW, i, j, wij);
                    }
                }
            }
        }
    }

    gretl_matrix_block_destroy(M);
    gretl_matrix_free(V);

    return err;
}

private int QR_OLS_work(GretlMatrixRaw * Q, GretlMatrixRaw * R)
{
    integer k = gretl_matrix_rows(R);
    int r, err;

    /* basic decomposition */
    err = gretl_matrix_QR_decomp(Q, R);
    if (err) {
        return err;
    }

    /* check rank of QR */
    r = gretl_check_QR_rank(R, &err, null);
    if (err) {
        return err;
    }

    if (r < k) {
        err = GretlError.singular;
    } else {
        /* invert R */
        char uplo = 'U';
        char diag = 'N';
        integer info = 0;

        dtrtri_(&uplo, &diag, &k, R.val, &k, &info);
        if (info != 0) {
            fprintf(stderr, "dtrtri: info = %d\n", cast(int) info);
            err = 1;
        }
    }

    return err;
}

/*
 * gretl_matrix_QR_ols:
 * @y: T x g matrix of dependent variables.
 * @X: T x k matrix of independent variables.
 * @B: k x g matrix to hold coefficient estimates.
 * @E: T x g matrix to hold the regression residuals, or null if these are
 * not needed.
 * @XTXi: location to receive (X'X)^{-1}, or null if this is not needed.
 * @Qout: location to receive Q on output, or null.
 *
 * Computes OLS estimates using QR decomposition, and puts the
 * coefficient estimates in @B.  Optionally, calculates the
 * residuals in @E, (X'X)^{-1} in @XTXi, and/or the matrix Q
 * in @Qout.
 *
 * Returns: 0 on success, non-zero error code on failure.
 */

int gretl_matrix_QR_ols(const(GretlMatrixRaw) * Y,
                        const(GretlMatrixRaw) * X,
                        GretlMatrixRaw * B,
                        GretlMatrixRaw * E,
                        GretlMatrixRaw ** XTXi,
                        GretlMatrixRaw ** Qout)
{
    int g = Y.cols;
    int k = X.cols;
    int T = X.rows;
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * G = null;
    int err = 0;

    if (B.rows != k || B.cols != g) {
        err = GretlError.nonconf;
    } else if (Y.rows != T) {
        err = GretlError.nonconf;
    } else if (E != null && (E.cols != g || E.rows != T)) {
        err = GretlError.nonconf;
    } else if (k > T) {
        err = GretlError.df;
    }

    if (!err) {
        Q = gretl_matrix_copy(X);
        R = gretl_matrix_alloc(k, k);
        G = gretl_matrix_alloc(k, g);
        if (Q == null || R == null || G == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        err = QR_OLS_work(Q, R);
    }

    if (!err) {
        /* make "G" into gamma-hat */
        gretl_matrix_multiply_mod(Q, GretlMatrixMod.transpose,
                                  Y, GretlMatrixMod.none,
                                  G, GretlMatrixMod.none);
    }

    if (!err) {
        /* OLS coefficients */
        gretl_matrix_multiply(R, G, B);
    }

    if (!err && E != null) {
        /* compute residuals */
        int i, imax = E.rows * E.cols;

        gretl_matrix_multiply(X, B, E);
        for (i=0; i<imax; i++) {
            E.val[i] = Y.val[i] - E.val[i];
        }
    }

    /* create (X'X)^{-1} = RR' */
    if (!err && XTXi != null) {
        *XTXi = gretl_matrix_alloc(k, k);
        if (*XTXi == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_multiply_mod(R, GretlMatrixMod.none,
                                      R, GretlMatrixMod.transpose,
                                      *XTXi, GretlMatrixMod.none);
        }
    }

    if (!err && Qout != null) {
        *Qout = Q;
    } else {
        gretl_matrix_free(Q);
    }

    gretl_matrix_free(R);
    gretl_matrix_free(G);

    return err;
}

/*
 * gretl_matrix_r_squared:
 * @y: dependent variable, T-vector.
 * @X: independent variables matrix, T x k.
 * @b: coefficients, k-vector.
 * @err: location to receive error code.
 *
 * Returns: the unadjusted R-squared, based on the regression
 * represented by @y, @X and @b, or NADBL on failure.
 */

double gretl_matrix_r_squared(const(GretlMatrixRaw) * y,
                              const(GretlMatrixRaw) * X,
                              const(GretlMatrixRaw) * b,
                              int * err)
{
    double ess = 0.0, tss = 0.0;
    double xx, ybar;
    int i, j;

    if (gretl_vector_get_length(y) != X.rows ||
        gretl_vector_get_length(b) != X.cols) {
        *err = GretlError.nonconf;
        return NADBL;
    }

    ybar = gretl_vector_mean(y);

    for (i=0; i<X.rows; i++) {
        xx = y.val[i];
        for (j=0; j<X.cols; j++) {
            xx -= b.val[j] * gretl_matrix_get(X, i, j);
        }
        ess += xx * xx;
        xx = y.val[i] - ybar;
        tss += xx * xx;
    }

    return 1.0 - ess / tss;
}

/*
 * gretl_matrix_columwise_product:
 * @A: T x k matrix.
 * @B: T x n matrix.
 * @S: k x n selection matrix, or null.
 * @C: T x (k*n) matrix to hold the product (but see below).
 *
 * If @S is null, computes a columnwise product in k blocks, each
 * of n columns. The first block consists of the Hadamard product
 * of the first column of @A and the matrix @B, the second block
 * holds the Hadamard product of the second column of @A and
 * matrix @B, and so on.
 *
 * A non-null @S matrix may be used to filter the column-pairs for
 * multiplication: @C will include the product of column i of @A
 * and column j of @B if and only if the i, j element of @S is
 * non-zero. In this case @C should have a number of columns
 * equal to the number of non-zero elements of @S.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_columnwise_product(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * B,
                                    const(GretlMatrixRaw) * S,
                                    GretlMatrixRaw * C)
{
    int k, n, T;
    double x, y;
    int i, j, t, p;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(B) ||
        gretl_is_null_matrix(C)) {
        return GretlError.data;
    }

    k = A.cols;
    n = B.cols;
    T = A.rows;

    if (B.rows != T || C.rows != T) {
        return GretlError.nonconf;
    }

    if (S != null) {
        if (S.rows != k || S.cols != n) {
            return GretlError.nonconf;
        } else {
            int c = 0;

            for (i=0; i<k*n; i++) {
                if (S.val[i] != 0) {
                    c++;
                }
            }
            if (C.cols != c) {
                return GretlError.nonconf;
            }
        }
    } else if (C.cols != k * n) {
        return GretlError.nonconf;
    }

    p = 0;
    for (i=0; i<k; i++) {
        for (j=0; j<n; j++) {
            if (S == null || gretl_matrix_get(S, i, j) != 0) {
                for (t=0; t<T; t++) {
                    x = gretl_matrix_get(A, t, i);
                    y = gretl_matrix_get(B, t, j);
                    gretl_matrix_set(C, t, p, x * y);
                }
                p++;
            }
        }
    }

    return 0;
}

private int alt_qform(const(GretlMatrixRaw) * A, GretlMatrixMod amod,
                      const(GretlMatrixRaw) * X, GretlMatrixRaw * C,
                      GretlMatrixMod cmod)
{
    GretlMatrixRaw * Tmp;
    int r = (amod) ? A.cols : A.rows;

    Tmp = gretl_matrix_alloc(r, X.cols);
    if (Tmp == null) {
        return GretlError.alloc;
    }

    if (amod == GretlMatrixMod.transpose) {
        /* A' * X * A */
        gretl_matrix_multiply_mod(A, GretlMatrixMod.transpose,
                                  X, GretlMatrixMod.none,
                                  Tmp, GretlMatrixMod.none);
        gretl_matrix_multiply_mod(Tmp, GretlMatrixMod.none,
                                  A, GretlMatrixMod.none,
                                  C, cmod);
    } else {
        /* A * X * A' */
        gretl_matrix_multiply(A, X, Tmp);
        gretl_matrix_multiply_mod(Tmp, GretlMatrixMod.none,
                                  A, GretlMatrixMod.transpose,
                                  C, cmod);
    }

    gretl_matrix_xtr_symmetric(C);
    gretl_matrix_free(Tmp);

    return 0;
}

/*
 * gretl_matrix_qform:
 * @A: m * k matrix or k * m matrix, depending on @amod.
 * @amod: GretlMatrixMod.none or GretlMatrixMod.transpose: in the first
 * case @A should be m * k; in the second, k * m;
 * @X: symmetric k * k matrix.
 * @C: matrix to hold the product.
 * @cmod: modifier: GretlMatrixMod.none, or GretlMatrixMod.cumulate to
 * add the result to the existing value of @C, or
 * GretlMatrixMod.decrement to subtract from the existing value
 * of @C.
 *
 * Computes either A * X * A' (if amod = GretlMatrixMod.none) or
 * A' * X * A (if amod = GretlMatrixMod.transpose), with the result
 * written into @C.  The matrix @X must be symmetric, but this
 * is not checked, to save time.  If you are in doubt on this
 * point you can call gretl_matrix_is_symmetric() first.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

private enum QFORM_SMALL = 1.0e-21;

int gretl_matrix_qform(const(GretlMatrixRaw) * A, GretlMatrixMod amod,
                       const(GretlMatrixRaw) * X, GretlMatrixRaw * C,
                       GretlMatrixMod cmod)
{
    int i, j, ii, jj;
    double xi, xj, xij, xx, cij;
    int m, k;
    guint64 N;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(X) ||
        gretl_is_null_matrix(C)) {
        return GretlError.data;
    } else if (A.is_complex || X.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_qform\n");
        if (A.is_complex) fprintf(stderr, "\touter is complex\n");
        if (X.is_complex) fprintf(stderr, "\tinner is complex\n");
        return GretlError.cmplx;
    }

    m = (amod) ? A.cols : A.rows;
    k = (amod) ? A.rows : A.cols;

    if (k != X.rows) {
        fprintf(stderr, "gretl_matrix_qform: %s is (%d x %d) but X is (%d x %d)\n",
                ((amod) ? "A'" : "A").ptr, m, k, X.rows, X.cols);
        return GretlError.nonconf;
    }

    if (C.rows != m || C.cols != m) {
        fputs("gretl_matrix_qform: destination matrix not conformable\n", stderr);
        return GretlError.nonconf;
    }

    N = cast(guint64) m * m * k * k;

    if (N > 80000) {
        /* take advantage of optimized matrix multiplication */
        return alt_qform(A, amod, X, C, cmod);
    }

    if (amod) {
        for (i=0; i<m; i++) {
            for (j=i; j<m; j++) {
                xx = 0.0;
                for (ii=0; ii<k; ii++) {
                    xi = gretl_matrix_get(A,ii,i);
                    if (fabs(xi) > QFORM_SMALL) {
                        for (jj=0; jj<k; jj++) {
                            xj = gretl_matrix_get(A,jj,j);
                            xij = gretl_matrix_get(X,ii,jj);
                            xx += xij * xi * xj;
                        }
                    }
                }
                if (cmod == GretlMatrixMod.cumulate) {
                    cij = gretl_matrix_get(C, i, j) + xx;
                } else if (cmod == GretlMatrixMod.decrement) {
                    cij = gretl_matrix_get(C, i, j) - xx;
                } else {
                    cij = xx;
                }
                gretl_matrix_set(C, i, j, cij);
                if (j != i) {
                    gretl_matrix_set(C, j, i, cij);
                }
            }
        }
    } else {
        for (i=0; i<m; i++) {
            for (j=i; j<m; j++) {
                xx = 0.0;
                for (ii=0; ii<k; ii++) {
                    xi = gretl_matrix_get(A,i,ii);
                    if (fabs(xi) > QFORM_SMALL) {
                        for (jj=0; jj<k; jj++) {
                            xj = gretl_matrix_get(A,j,jj);
                            xij = gretl_matrix_get(X,ii,jj);
                            xx += xij * xi * xj;
                        }
                    }
                }
                if (cmod == GretlMatrixMod.cumulate) {
                    cij = gretl_matrix_get(C, i, j) + xx;
                } else if (cmod == GretlMatrixMod.decrement) {
                    cij = gretl_matrix_get(C, i, j) - xx;
                } else {
                    cij = xx;
                }
                gretl_matrix_set(C, i, j, cij);
                if (j != i) {
                    gretl_matrix_set(C, j, i, cij);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_matrix_diag_qform:
 * @A: m * k matrix or k * m matrix, depending on @amod.
 * @amod: GretlMatrixMod.none or GretlMatrixMod.transpose: in the first
 * case @A should be m * k; in the second, k * m;
 * @d: k-element vector.
 * @C: matrix to hold the product.
 * @cmod: modifier: GretlMatrixMod.none, or GretlMatrixMod.cumulate to
 * add the result to the existing value of @C, or
 * GretlMatrixMod.decrement to subtract from the existing value of
 * @C.
 *
 * Computes either A * md * A' (if amod = GretlMatrixMod.none) or A' *
 * md * A (if amod = GretlMatrixMod.transpose), where md is a diagonal
 * matrix holding the vector @d. The result is written into @C.
 *
 * Returns: 0 on success; non-zero error code on failure.
 */

int gretl_matrix_diag_qform(const(GretlMatrixRaw) * A, GretlMatrixMod amod,
                            const(GretlVector) * d, GretlMatrixRaw * C,
                            GretlMatrixMod cmod)
{
    int i, j, k;
    double x, xi, xj, cij;
    int r, c, ld;
    int dim, hybrid_min = 9500;
    int use_hybrid = 1;

    if (gretl_is_null_matrix(A) ||
        gretl_is_null_matrix(d) ||
        gretl_is_null_matrix(C)) {
        return GretlError.data;
    } else if (A.is_complex || d.is_complex) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_diag_qform\n");
        if (A.is_complex) fprintf(stderr, "\touter is complex\n");
        if (d.is_complex) fprintf(stderr, "\tinner is complex\n");
        return GretlError.cmplx;
    }

    r = (amod) ? A.cols : A.rows;
    c = (amod) ? A.rows : A.cols;
    ld = gretl_vector_get_length(d);

    if (c != ld) {
        fprintf(stderr, "gretl_matrix_diag_qform: %s is (%d x %d) but d has %d elements\n",
                ((amod) ? "A'" : "A").ptr, r, c, ld);
        return GretlError.nonconf;
    }

    if (C.rows != r || C.cols != r) {
        fputs("gretl_matrix_diag_qform: destination matrix not conformable\n", stderr);
        return GretlError.nonconf;
    }

    if (hybrid_min > 0) {
        char * s = getenv("HYBRID_MIN");

        if (s != null) {
            hybrid_min = atoi(s);
        }
        /* condition on the number of flops required */
        dim = cast(int) (3 * c * 0.5 * (r+1) * r);
        use_hybrid = (dim > hybrid_min);
    }

    if (use_hybrid) {
        /* hybrid of special code and optimized matrix multiplication */
        GretlMatrixRaw * AD = gretl_matrix_alloc(A.rows, A.cols);

        if (AD != null) {
            k = 0;
            if (amod) {
                /* here AD is actually <d> A */
                for (i=0; i<A.rows; i++) {
                    x = d.val[i];
                    k = i;
                    for (j=0; j<A.cols; j++) {
                        AD.val[k] = A.val[k] * x;
                        k += A.rows;
                    }
                }
                gretl_matrix_multiply_mod(AD, GretlMatrixMod.transpose,
                                          A, GretlMatrixMod.none,
                                          C, cmod);
            } else {
                /* here AD is A <d> */
                for (i=0; i<A.cols; i++) {
                    x = d.val[i];
                    for (j=0; j<A.rows; j++) {
                        AD.val[k] = A.val[k] * x;
                        k++;
                    }
                }
                gretl_matrix_multiply_mod(AD, GretlMatrixMod.none,
                                          A, GretlMatrixMod.transpose,
                                          C, cmod);
            }
            gretl_matrix_free(AD);
        }
    } else {
        /* fully specialized code which works well for small input */
        for (i=0; i<r; i++) {
            for (j=0; j<=i; j++) {
                x = 0.0;
                for (k=0; k<c; k++) {
                    if (amod) {
                        xi = gretl_matrix_get(A,k,i);
                        xj = gretl_matrix_get(A,k,j);
                    } else {
                        xi = gretl_matrix_get(A,i,k);
                        xj = gretl_matrix_get(A,j,k);
                    }
                    x += d.val[k] * xi * xj;
                }
                if (cmod == GretlMatrixMod.cumulate) {
                    cij = gretl_matrix_get(C, i, j) + x;
                } else if (cmod == GretlMatrixMod.decrement) {
                    cij = gretl_matrix_get(C, i, j) - x;
                } else {
                    cij = x;
                }
                gretl_matrix_set(C, i, j, cij);
                if (j != i) {
                    gretl_matrix_set(C, j, i, cij);
                }
            }
        }
    }

    return 0;
}

/*
 * gretl_scalar_qform:
 * @b: k-vector.
 * @X: symmetric k x k matrix.
 * @err: pointer to receive error code.
 *
 * Computes the scalar product bXb', or b'Xb if @b is a column
 * vector. The content of @err is set to a non-zero code on
 * failure.
 *
 * Returns: the scalar product, or NADBL on failure.
 */

double gretl_scalar_qform(const(GretlVector) * b,
                          const(GretlMatrixRaw) * X,
                          int * err)
{
    double tmp, ret = 0.0;
    int i, j, k, p;

    if (gretl_is_null_matrix(b) || gretl_is_null_matrix(X)) {
        *err = GretlError.data;
        return NADBL;
    }

    k = gretl_vector_get_length(b);

    if (k == 0 || X.rows != k || X.cols != k) {
        *err = GretlError.nonconf;
        return NADBL;
    }

    p = 0;
    for (j=0; j<k; j++) {
        tmp = 0.0;
        for (i=0; i<k; i++) {
            tmp += b.val[i] * X.val[p++];
        }
        ret += tmp * b.val[j];
    }

    return ret;
}

/*
 * gretl_matrix_diagonal_sandwich:
 * @d: k-vector.
 * @X: k * k matrix.
 * @DXD: target k * k matrix.
 *
 * Computes in @DXD (which must be pre-allocated), the matrix
 * product D * X * D, where D is a diagonal matrix with elements
 * given by the vector @d.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_matrix_diagonal_sandwich(const(GretlVector) * d,
                                   const(GretlMatrixRaw) * X,
                                   GretlMatrixRaw * DXD)
{
    int dim = (d.rows == 1) ? d.cols : d.rows;
    double x, xij;
    int i, j, err = 0;

    if (dim != X.rows || dim != X.cols ||
        dim != DXD.rows || dim != DXD.cols) {
        err = GretlError.nonconf;
    } else {
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                xij = gretl_matrix_get(X, i, j);
                x = xij * d.val[i] * d.val[j];
                gretl_matrix_set(DXD, i, j, x);
            }
        }
    }

    return err;
}

/*
 * gretl_is_identity_matrix:
 * @m: matrix to examine.
 *
 * Returns: 1 if @m is an identity matrix, 0 otherwise.
 */

int gretl_is_identity_matrix(const(GretlMatrixRaw) * m)
{
    double x;
    int i, j;

    if (gretl_is_null_matrix(m)) {
        return 0;
    } else if (m.rows != m.cols) {
        return 0;
    }

    for (j=0; j<m.cols; j++) {
        for (i=0; i<m.rows; i++) {
            x = gretl_matrix_get(m, i, j);
            if (i == j && x != 1.0) return 0;
            if (i != j && x != 0.0) return 0;
        }
    }

    return 1;
}

/*
 * gretl_is_zero_matrix:
 * @m: matrix to examine.
 *
 * Returns: 1 if @m is a zero matrix, 0 otherwise.
 */

int gretl_is_zero_matrix(const(GretlMatrixRaw) * m)
{
    int i, n;

    if (gretl_is_null_matrix(m)) {
        return 0;
    }

    n = m.rows * m.cols;

    for (i=0; i<n; i++) {
        if (m.val[i] != 0.0) {
            return 0;
        }
    }

    return 1;
}

/*
 * gretl_matrix_isfinite:
 * @m: matrix to examine.
 * @err: location to receive error code.
 *
 * Returns: a matrix with 1s in positions corresponding to
 * finite elements of @m, zeros otherwise.
 */

GretlMatrixRaw * gretl_matrix_isfinite(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * f;

    if (m == null) {
        *err = GretlError.data;
        return null;
    }

    f = gretl_matrix_alloc(m.rows, m.cols);

    if (f == null) {
        *err = GretlError.alloc;
    } else {
        int i, n = m.rows * m.cols;

        for (i=0; i<n; i++) {
            f.val[i] = (na(m.val[i])) ? 0 : 1;
        }
    }

    return f;
}

/*
 * gretl_matrices_are_equal:
 * @a: first matrix in comparison.
 * @b: second matrix in comparison.
 * @tol: numerical tolerance.
 * @err: location to receive error code.
 *
 * Returns: 1 if the matrices @a and @b compare equal, 0 if they
 * differ, and -1 if the comparison is invalid, in which case
 * GretlError.nonconf is written to @err.
 */

int gretl_matrices_are_equal(const(GretlMatrixRaw) * a,
                             const(GretlMatrixRaw) * b,
                             double tol, int * err)
{
    double ax, bx;
    int nas;
    int i, n;

    if (a == null || b == null) {
        *err = GretlError.data;
        return -1;
    }

    if (a == b) {
        return 1;
    } else if (a.rows != b.rows || a.cols != b.cols) {
        *err = GretlError.nonconf;
        return -1;
    } else if (a.is_complex + b.is_complex > 0) {
        *err = GretlError.cmplx;
        return -1;
    }

    n = a.rows * a.cols;

    for (i=0; i<n; i++) {
        ax = a.val[i];
        bx = b.val[i];
        nas = na(ax) + na(bx);
        if (nas == 1) {
            return 0;
        } else if (nas == 0 && fabs(ax - bx) > tol) {
            int j = i / a.rows;
            int k = i % a.rows;

            fprintf(stderr, "gretl_matrices_are_equal:\n "
                    ~ "at row %d, col %d (1-based)\n aij = %.15g but bij = %.15g\n",
                    k + 1, j + 1, ax, bx);
            return 0;
        }
    }

    return 1;
}

/*
 * gretl_covariance_matrix:
 * @m: (n x k) matrix containing n observations on each of k
 * variables.
 * @corr: flag for computing correlations.
 * @dfc: degrees of freedom correction: use 1 for sample
 * variance, 0 for MLE.
 * @err: pointer to receive non-zero error code in case of
 * failure, or null.
 *
 * Returns: the covariance matrix of variables in the columns of
 * @m, or the correlation matrix if @corr is non-zero.
 */

GretlMatrixRaw * gretl_covariance_matrix(const(GretlMatrixRaw) * m,
                                      int corr, int dfc,
                                      int * err)
{
    GretlMatrixRaw * D = null;
    GretlMatrixRaw * V = null;

    if (gretl_is_null_matrix(m) || dfc < 0) {
        *err = GretlError.invarg;
    } else if (m.rows < 2 || dfc >= m.rows) {
        *err = GretlError.toofew;
    } else {
        D = gretl_matrix_copy(m);
        if (D == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        return null;
    }

    if (corr) {
        gretl_matrix_standardize(D, dfc, 0);
    } else {
        gretl_matrix_center(D, 0);
    }

    V = gretl_matrix_XTX_new(D);
    if (V == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_divide_by_scalar(V, m.rows - dfc);
    }

    gretl_matrix_free(D);

    return V;
}

private GretlMatrixRaw *
real_gretl_matrix_values(const(double) * x, int n,
                         Gretlopt opt, int * n_vals,
                         int * missvals, int * err)
{
    GretlMatrixRaw * v = null;
    double * sorted = null;
    double last;
    int m = 0;
    int i, k;

    sorted = cast(double *) malloc(n * (*sorted).sizeof);
    if (sorted == null) {
        *err = GretlError.alloc;
        return null;
    }

    k = 0;
    for (i=0; i<n; i++) {
        if (!na(x[i])) {
            sorted[k++] = x[i];
        }
    }

    if (k == 0) {
        if (n_vals == null) {
            v = gretl_null_matrix_new();
        }
        goto bailout;
    } else if (missvals != null) {
        *missvals = n - k;
    }

    qsort(sorted, k, (*sorted).sizeof, &gretl_compare_doubles);
    m = count_distinct_values(sorted, k);

    if (n_vals != null) {
        /* the caller just wants the count */
        *n_vals = m;
        goto bailout;
    }

    v = gretl_column_vector_alloc(m);
    if (v == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (opt & Gretlopt.OPT_S) {
        /* sorted */
        v.val[0] = last = sorted[0];
        for (i=1, m=1; i<k; i++) {
            if (sorted[i] != last) {
                last = sorted[i];
                v.val[m++] = sorted[i];
            }
        }
    } else {
        /* unsorted */
        int j, add;

        for (i=0, m=0; i<n; i++) {
            if (!na(x[i])) {
                add = 1;
                for (j=0; j<m; j++) {
                    if (v.val[j] == x[i]) {
                        add = 0;
                        break;
                    }
                }
                if (add) {
                    v.val[m++] = x[i];
                }
            }
        }
    }

 bailout:

    free(sorted);

    return v;
}

/*
 * gretl_matrix_values:
 * @x: array to process.
 * @n: length of array.
 * @opt: if OPT_S the array of values will be sorted, otherwise
 * given in order of occurrence.
 * @err: location to receive error code.
 *
 * Returns: an allocated matrix containing the distinct
 * values in array @x, skipping any missing values, or
 * null on failure.
 */

GretlMatrixRaw * gretl_matrix_values(const(double) * x, int n,
                                  Gretlopt opt, int * err)
{
    return real_gretl_matrix_values(x, n, opt, null, null, err);
}

/*
 * gretl_matrix_values_full:
 * @x: array to process.
 * @n: length of array.
 * @opt: if OPT_S the array of values will be sorted, otherwise
 * given in order of occurrence.
 * @missvals: location to receive count of missing values.
 * @err: location to receive error code.
 *
 * Returns: an allocated matrix containing the distinct
 * values in array @x, skipping any missing values, or
 * null on failure.
 */

GretlMatrixRaw * gretl_matrix_values_full(const(double) * x, int n,
                                       Gretlopt opt, int * missvals,
                                       int * err)
{
    return real_gretl_matrix_values(x, n, opt, null, missvals, err);
}

/*
 * gretl_matrix_n_values:
 * @x: array to process.
 * @n: length of array.
 * @err: location to receive error code.
 *
 * Returns: the number of distinct non-missing values in array @x,
 * or 0 on failure.
 */

int gretl_matrix_n_values(const(double) * x, int n, int * err)
{
    int n_vals = 0;

    real_gretl_matrix_values(x, n, Gretlopt.OPT_NONE, &n_vals, null, err);

    return n_vals;
}

/*
 * gretl_matrix_shape:
 * @A: array to process.
 * @r: rows of target matrix.
 * @c: columns of target matrix.
 *
 * Creates an (r x c) matrix containing the re-arranged
 * values of A.  Elements are read from A by column and
 * written into the target, also by column.  If A contains
 * less elements than n = r*c, they are repeated cyclically;
 * if A has more elements, only the first n are used.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_shape(const(GretlMatrixRaw) * A,
                                 int r, int c, int * err)
{
    GretlMatrixRaw * B = null;

    if (gretl_is_null_matrix(A) || r < 0 || c < 0) {
        *err = GretlError.invarg;
        return null;
    }

    if (r == 0 && c == 0) {
        return gretl_null_matrix_new();
    }

    if (A.is_complex) {
        B = gretl_cmatrix_new(r, c);
    } else {
        B = gretl_matrix_alloc(r, c);
    }

    if (B == null) {
        *err = GretlError.alloc;
    } else {
        int nA = A.rows * A.cols;
        int nB = r * c;
        int i, k = 0;

        if (A.is_complex) {
            nA *= 2;
            nB *= 2;
        }

        k = 0;
        for (i=0; i<nB; i++) {
            B.val[i] = A.val[k++];
            if (k == nA) {
                k = 0;
            }
        }
    }

    return B;
}

/*
 * gretl_matrix_trim_rows:
 * @A: array to process.
 * @ttop: rows to trim at top.
 * @tbot: rows to trim at bottom.
 * @err: location to receive error code.
 *
 * Creates a new matrix which is a copy of @A with @ttop rows
 * trimmed from the top and @tbot rows trimmed from the
 * bottom.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_trim_rows(const(GretlMatrixRaw) * A,
                                     int ttop, int tbot,
                                     int * err)
{
    GretlMatrixRaw * B;
    Cmplx z;
    double x;
    int i, j, m;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return null;
    }

    m = A.rows - (ttop + tbot);

    if (ttop < 0 || tbot < 0 || m <= 0) {
        *err = GretlError.data;
        return null;
    }

    B = gretl_matching_matrix_new(m, A.cols, A);

    if (B == null) {
        *err = GretlError.alloc;
    } else {
        for (j=0; j<A.cols; j++) {
            for (i=0; i<m; i++) {
                if (A.is_complex) {
                    z = gretl_cmatrix_get(A, i + ttop, j);
                    gretl_cmatrix_set(B, i, j, z);
                } else {
                    x = gretl_matrix_get(A, i + ttop, j);
                    gretl_matrix_set(B, i, j, x);
                }
            }
        }
    }

    return B;
}

/*
 * gretl_matrix_minmax:
 * @A: m x n matrix to examine.
 * @mm: 0 for minima, 1 for maxima.
 * @rc: 0 for row-wise, 1 for column-wise.
 * @idx: 0 for values, 1 for indices.
 * @skip_na: ignore missing values.
 * @err: location to receive error code.
 *
 * Creates a matrix holding the row or column mimima or
 * maxima from @A, either as values or as location indices.
 * For example, if @mm = 0, @rc = 0, and @idx = 0, the
 * created matrix is m x 1 and holds the values of the row
 * minima.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_minmax(const(GretlMatrixRaw) * A,
                                  int mm, int rc, int idx,
                                  int skip_na, int * err)
{
    GretlMatrixRaw * B;
    double d, x;
    int i, j, k;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return null;
    }

    if (rc == 0) {
        B = gretl_zero_matrix_new(A.rows, 1);
    } else {
        B = gretl_zero_matrix_new(1, A.cols);
    }

    if (B == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (rc == 0) {
        /* going by rows */
        for (i=0; i<A.rows; i++) {
            int valid_cols = A.cols;

            d = mm > 0 ? -double.max : double.max;
            k = 0;
            for (j=0; j<A.cols; j++) {
                x = gretl_matrix_get(A, i, j);
                if (na(x)) {
                    if (skip_na) {
                        valid_cols--;
                        continue;
                    } else {
                        B.val[i] = NADBL;
                        break;
                    }
                } else if (mm > 0) {
                    /* looking for max */
                    if (x > d) {
                        d = x;
                        k = j;
                    }
                } else {
                    /* looking for min */
                    if (x < d) {
                        d = x;
                        k = j;
                    }
                }
            }
            if (valid_cols == 0) {
                B.val[i] = NADBL;
            }
            if (!na(B.val[i])) {
                B.val[i] = idx ? cast(double) k + 1 : d;
            }
        }
    } else {
        /* going by columns */
        for (j=0; j<A.cols; j++) {
            int valid_rows = A.rows;

            d = mm > 0 ? -double.max : double.max;
            k = 0;
            for (i=0; i<A.rows; i++) {
                x = gretl_matrix_get(A, i, j);
                if (na(x)) {
                    if (skip_na) {
                        valid_rows--;
                        continue;
                    } else {
                        B.val[j] = NADBL;
                        break;
                    }
                } else if (mm > 0) {
                    /* looking for max */
                    if (x > d) {
                        d = x;
                        k = i;
                    }
                } else {
                    /* looking for min */
                    if (x < d) {
                        d = x;
                        k = i;
                    }
                }
            }
            if (valid_rows == 0) {
                B.val[j] = NADBL;
            }
            if (!na(B.val[j])) {
                B.val[j] = idx ? cast(double) k + 1 : d;
            }
        }
    }

    return B;
}

/*
 * gretl_matrix_global_minmax:
 * @A: matrix to examine.
 * @mm: 0 for minimum, 1 for maximum.
 * @err: location to receive error code.
 *
 * Returns: the smallest or greatest element of @A,
 * ignoring NaNs but not infinities (?), or NADBL on
 * failure.
 */

double gretl_matrix_global_minmax(const(GretlMatrixRaw) * A,
                                  int mm, int * err)
{
    double x, ret = NADBL;
    int i, n, started = 0;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return NADBL;
    }

    n = A.rows * A.cols;

    for (i=0; i<n; i++) {
        x = A.val[i];
        if (isNaN(x)) {
            /* (MINMAX_IGNORE_NAN is 1 in the C source) */
            continue;
        } else if (!started) {
            ret = x;
            started = 1;
        } else if ((mm == 0 && x < ret) ||
                   (mm == 1 && x > ret)) {
            ret = x;
        }
    }

    return ret;
}

/*
 * gretl_matrix_global_sum:
 * @A: matrix to examine.
 * @err: location to receive error code.
 *
 * Returns: the sum of the elements of @A,
 * or NADBL on failure.
 */

double gretl_matrix_global_sum(const(GretlMatrixRaw) * A,
                               int * err)
{
    double ret = 0.0;
    int i, n;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return NADBL;
    } else if (A.is_complex) {
        *err = GretlError.cmplx;
        return NADBL;
    }

    n = A.rows * A.cols;

    for (i=0; i<n; i++) {
        ret += A.val[i];
        if (isNaN(ret)) {
            break;
        }
    }

    return ret;
}

GretlMatrixRaw * real_gretl_matrix_pca(const(GretlMatrixRaw) * X,
                                            int p, Gretlopt opt,
                                            int * err)
{
    GretlMatrixRaw * D = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * P = null;
    GretlMatrixRaw * e;

    if (gretl_is_null_matrix(X)) {
        *err = GretlError.data;
        return null;
    }

    if (p <= 0 || p > X.cols) {
        *err = GretlError.invarg;
        return null;
    } else if (X.rows < 2) {
        *err = GretlError.toofew;
        return null;
    } else if (X.is_complex) {
        *err = GretlError.cmplx;
        return null;
    }

    D = gretl_matrix_copy(X);
    if (D == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (opt & Gretlopt.OPT_V) {
        /* use covariance matrix */
        gretl_matrix_center(D, 0);
    } else {
        /* use correlation matrix */
        gretl_matrix_standardize(D, 1, 0);
    }

    V = gretl_matrix_XTX_new(D);
    if (V == null) {
        *err = GretlError.alloc;
    } else {
        /* note: we don't need the eigenvalues of V, but if we
           don't grab and then free the return value below,
           we'll leak a gretl_matrix
        */
        e = real_symm_eigenvals_descending(V, 1, p, err);
        gretl_matrix_free(e);
    }

    if (!*err) {
        P = gretl_matrix_multiply_new(D, V, err);
    }

    gretl_matrix_free(D);
    gretl_matrix_free(V);

    return P;
}

/* NOTE: the C source has `#include "princomp_extra.c"` at this point; in
   the port that file is the separate module princomp_extra.d (imported at
   the top). real_gretl_matrix_pca is public so that module can call it. */

/*
 * gretl_matrix_pca:
 * @X: T x k data matrix.
 * @p: number of principal components to return: 0 < p <= k.
 * @opt: if OPT_V, use the covariance matrix rather than the
 * correlation matrix as basis.
 * @err: location to receive error code.
 *
 * Carries out a Principal Components analysis of @X and
 * returns the first @p components: the component corresponding
 * to the largest eigenvalue of the correlation matrix of @X
 * is placed in column 1, and so on.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_pca(const(GretlMatrixRaw) * X, int p,
                               Gretlopt opt, int * err)
{
    if (matrix_is_incomplete(X)) {
        return princomp_with_NAs(X, p, opt, err);
    } else {
        return real_gretl_matrix_pca(X, p, opt, err);
    }
}

private bool complete_obs(const(double) * x, const(double) * y, int i)
{
    return (!na(x[i]) && !na(y[i]));
}

private int ok_xy_count(const(double) * x,
                        const(double) * y,
                        int n, int * err)
{
    int i, n_ok = 0;

    for (i=0; i<n; i++) {
        if (complete_obs(x, y, i)) {
            n_ok++;
        }
    }

    if (!*err && n_ok < 2) {
        *err = GretlError.toofew;
    }

    return n_ok;
}

private void make_matrix_xtab(double ** X, int n,
                              const(GretlMatrixRaw) * vx,
                              const(GretlMatrixRaw) * vy,
                              GretlMatrixRaw * tab)
{
    double xr, xc;
    int rndx, cndx;
    int counter, i;

    qsort(X, n, (*X).sizeof, &compare_xtab_rows);

    /* compute frequencies by going through sorted X */

    counter = rndx = cndx = 0;
    xr = vx.val[0];
    xc = vy.val[0];

    for (i=0; i<n; i++) {
        while (X[i][0] > xr) {
            /* skip row */
            gretl_matrix_set(tab, rndx, cndx, cast(double) counter);
            counter = 0;
            xr = vx.val[++rndx];
            cndx = 0;
            xc = vy.val[0];
        }
        while (X[i][1] > xc) {
            /* skip column */
            gretl_matrix_set(tab, rndx, cndx, cast(double) counter);
            counter = 0;
            xc = vy.val[++cndx];
        }
        counter++;
    }

    gretl_matrix_set(tab, rndx, cndx, cast(double) counter);
}

/*
 * gretl_matrix_xtab:
 * @x: data array.
 * @y: data array.
 * @n: length of the two arrays.
 * @err: location to receive error code.
 *
 * Computes the cross tabulation of the values contained in the
 * arrays @x (by row) and @y (by column). These should generally
 * be discrete values otherwise the cross-tabulation may be
 * very large and uninformative.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_xtab(const(double) * x,
                                const(double) * y,
                                int n, int * err)
{
    GretlMatrixRaw * tab = null;
    GretlMatrixRaw * vx = null;
    GretlMatrixRaw * vy = null;
    int complete;
    int i, j, n_ok;

    *err = 0;

    n_ok = ok_xy_count(x, y, n, err);
    if (*err) {
        return null;
    }

    complete = (n_ok == n);

    if (complete) {
        /* no need to dodge missing values */
        vx = gretl_matrix_values(x, n, Gretlopt.OPT_S, err);
        if (!*err) {
            vy = gretl_matrix_values(y, n, Gretlopt.OPT_S, err);
        }
    } else {
        double * tmp = cast(double *) malloc(2 * n_ok * double.sizeof);
        double * okx, oky;

        if (tmp == null) {
            *err = GretlError.alloc;
        } else {
            okx = tmp;
            oky = tmp + n_ok;
            for (i=0, j=0; i<n; i++) {
                if (complete_obs(x, y, i)) {
                    okx[j] = x[i];
                    oky[j] = y[i];
                    j++;
                }
            }
            vx = gretl_matrix_values(okx, n_ok, Gretlopt.OPT_S, err);
            if (!*err) {
                vy = gretl_matrix_values(oky, n_ok, Gretlopt.OPT_S, err);
            }
            free(tmp);
        }
    }

    if (!*err) {
        tab = gretl_zero_matrix_new(vx.rows, vy.rows);
        if (tab == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        double ** X = doubles_array_new(n_ok, 2);

        if (X == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0, j=0; i<n; i++) {
                if (complete || complete_obs(x, y, i)) {
                    X[j][0] = x[i];
                    X[j][1] = y[i];
                    j++;
                }
            }
        }
        make_matrix_xtab(X, n_ok, vx, vy, tab);
        doubles_array_free(X, n_ok);
    }

    gretl_matrix_free(vx);
    gretl_matrix_free(vy);

    return tab;
}

/*
 * gretl_matrix_bool_sel:
 * @A: matrix.
 * @sel: selection vector.
 * @rowsel: row/column mode selector.
 * @err: location to receive error code.
 *
 * If @rowsel = 1, constructs a matrix which contains the rows
 * of A corresponding to non-zero values in the vector @sel;
 * if @rowsel = 0, does the same thing but column-wise.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_bool_sel(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * sel,
                                    int rowsel, int * err)
{
    GretlMatrixRaw * ret = null;
    int nonzero = 0;
    int ra, ca, rs, cs;
    int rret, cret;
    int i, j, k, n;
    double x;

    *err = 0;

    if (gretl_is_null_matrix(A)) {
        return gretl_null_matrix_new();
    } else if (sel.is_complex) {
        *err = GretlError.invarg;
        return null;
    }

    ra = A.rows;
    ca = A.cols;
    rs = sel.rows;
    cs = sel.cols;

    /* check dimensions */
    if (rowsel) {
        if ((ra != rs) || (cs > 1)) {
            *err = GretlError.nonconf;
            return null;
        }
    } else {
        if ((ca != cs) || (rs > 1)) {
            *err = GretlError.nonconf;
            return null;
        }
    }

    /* count nonzeros */
    n = (rowsel) ? rs : cs ;
    for (i=0; i<n; i++) {
        x = gretl_vector_get(sel, i);
        if (na(x)) {
            *err = GretlError.missdata;
            return null;
        } else if (x != 0) {
            nonzero++;
        }
    }

    /* check for extreme cases */
    if (nonzero == n) {
        ret = gretl_matrix_copy(A);
        goto bailout;
    } else if (nonzero == 0) {
        ret = gretl_null_matrix_new();
        goto bailout;
    }

    rret = rowsel ? nonzero : ra;
    cret = rowsel ? ca : nonzero;

    if (A.is_complex) {
        ret = gretl_cmatrix_new(rret, cret);
    } else {
        ret = gretl_matrix_alloc(rret, cret);
    }
    if (ret == null) {
        goto bailout;
    }

    /* copy selected row/columns */
    if (rowsel) {
        /* selection of rows */
        Cmplx z;

        k = 0;
        for (i=0; i<ra; i++) {
            if (gretl_vector_get(sel, i) != 0) {
                for (j=0; j<ca; j++) {
                    if (A.is_complex) {
                        z = gretl_cmatrix_get(A, i, j);
                        gretl_cmatrix_set(ret, k, j, z);
                    } else {
                        x = gretl_matrix_get(A, i, j);
                        gretl_matrix_set(ret, k, j, x);
                    }
                }
                k++;
            }
        }
    } else {
        /* selection of columns */
        double * targ = ret.val;
        const(double) * src = A.val;
        int rdim = A.is_complex ? 2 : 1;
        size_t colsize = ra * rdim * (*src).sizeof;

        for (j=0; j<ca; j++) {
            if (gretl_vector_get(sel, j) != 0) {
                memcpy(targ, src, colsize);
                targ += ra * rdim;
            }
            src += ra * rdim;
        }
    }

    if (rowsel) {
        maybe_preserve_names(ret, A, NameFlags.rownames, sel);
        maybe_preserve_names(ret, A, NameFlags.colnames, null);
    } else {
        maybe_preserve_names(ret, A, NameFlags.colnames, sel);
        maybe_preserve_names(ret, A, NameFlags.rownames, null);
    }

 bailout:

    if (ret == null) {
        *err = GretlError.alloc;
    }

    return ret;
}

private int unstable_comp(double a, double b)
{
    int ret = 0;

    if (isNaN(a) || isNaN(b)) {
        if (!isNaN(a)) {
            ret = -1;
        } else if (!isNaN(b)) {
            ret = 1;
        }
    } else {
        ret = cast(int) (a > b) - cast(int) (a < b);
    }

    return ret;
}

/* (extern (C) so these can be handed to qsort) */

private extern (C) int compare_values(const(void) * a, const(void) * b)
{
    const(double) * da = cast(const(double) *) a;
    const(double) * db = cast(const(double) *) b;
    int ret = unstable_comp(*da, *db);

    if (ret == 0) {
        /* ensure stable sort */
        ret = (cast(const(ubyte) *) a) - (cast(const(ubyte) *) b) > 0 ? 1 : -1;
    }

    return ret;
}

private extern (C) int inverse_compare_values(const(void) * a, const(void) * b)
{
    const(double) * da = cast(const(double) *) a;
    const(double) * db = cast(const(double) *) b;
    int ret = unstable_comp(*db, *da);

    if (ret == 0) {
        /* ensure stable sort */
        ret = (cast(const(ubyte) *) a) - (cast(const(ubyte) *) b) > 0 ? 1 : -1;
    }

    return ret;
}

/*
 * gretl_matrix_sort_by_column:
 * @m: matrix.
 * @k: column by which to sort.
 * @err: location to receive error code.
 *
 * Produces a matrix which contains the rows of @m, re-
 * ordered by increasing value of the elements in column
 * @k.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_sort_by_column(const(GretlMatrixRaw) * m,
                                          int k, int * err)
{
    struct RSort {
        double x;
        int row;
    }
    RSort * rs;
    GretlMatrixRaw * a;
    double x;
    int i, j;

    if (gretl_is_null_matrix(m) || k < 0 || k >= m.cols) {
        *err = GretlError.data;
        return null;
    }

    rs = cast(RSort *) malloc(m.rows * (*rs).sizeof);
    if (rs == null) {
        *err = GretlError.alloc;
        return null;
    }

    a = gretl_matrix_copy(m);
    if (a == null) {
        free(rs);
        *err = GretlError.alloc;
        return null;
    }

    for (i=0; i<m.rows; i++) {
        rs[i].x = gretl_matrix_get(m, i, k);
        rs[i].row = i;
    }

    qsort(rs, m.rows, (*rs).sizeof, &compare_values);

    for (j=0; j<m.cols; j++) {
        for (i=0; i<m.rows; i++) {
            x = gretl_matrix_get(m, rs[i].row, j);
            gretl_matrix_set(a, i, j, x);
        }
    }

    if (a.info != null && a.info.rownames != null) {
        char ** S = cast(char **) malloc(a.rows * (char *).sizeof);

        if (S != null) {
            for (i=0; i<a.rows; i++) {
                S[i] = a.info.rownames[i];
            }
            for (i=0; i<a.rows; i++) {
                a.info.rownames[i] = S[rs[i].row];
            }
            free(S);
        }
    }

    free(rs);

    return a;
}

private struct RSortX {
    double * x;
    int n;
    int row;
}

private extern (C) int compare_rsortx(const(void) * a, const(void) * b)
{
    const(RSortX) * rsa = cast(const(RSortX) *) a;
    const(RSortX) * rsb = cast(const(RSortX) *) b;
    int i, ret = 0;

    for (i=0; i<rsa.n; i++) {
        ret = unstable_comp(rsa.x[i], rsb.x[i]);
        if (ret != 0) {
            break;
        } else {
            /* ensure stable sort and continue */
            ret = (cast(const(ubyte) *) a) - (cast(const(ubyte) *) b) > 0 ? 1 : -1;
        }
    }

    return ret;
}

GretlMatrixRaw * gretl_matrix_sort_by_columns(const(GretlMatrixRaw) * m,
                                           int * cols, int * err)
{
    RSortX * rsx = null;
    GretlMatrixRaw * a = null;
    double * rvals = null;
    double * rptr = null;
    double x;
    int ns, i, j;

    if (gretl_is_null_matrix(m)) {
        *err = GretlError.data;
        return null;
    }

    if (cols == null || cols[0] == 0) {
        return gretl_matrix_copy(m);
    }

    /* check out the @cols argument */
    ns = cols[0];
    if (ns > m.cols) {
        *err = GretlError.invarg;
    } else {
        for (i=1; i<=ns; i++) {
            if (cols[i] < 0 || cols[i] >= m.cols) {
                *err = GretlError.invarg;
                break;
            }
        }
    }
    if (*err) {
        return null;
    }

    rsx = cast(RSortX *) malloc(m.rows * (*rsx).sizeof);
    rvals = cast(double *) malloc(m.rows * ns * (*rvals).sizeof);
    if (rsx == null || rvals == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    a = gretl_matrix_copy(m);
    if (a == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* fill the rsortx structs */
    rptr = rvals;
    for (i=0; i<m.rows; i++) {
        rsx[i].x = rptr;
        for (j=0; j<ns; j++) {
            rsx[i].x[j] = gretl_matrix_get(m, i, cols[j+1]);
        }
        rsx[i].n = ns;
        rsx[i].row = i;
        rptr += ns;
    }

    qsort(rsx, m.rows, (*rsx).sizeof, &compare_rsortx);

    for (j=0; j<m.cols; j++) {
        for (i=0; i<m.rows; i++) {
            x = gretl_matrix_get(m, rsx[i].row, j);
            gretl_matrix_set(a, i, j, x);
        }
    }

    if (a.info != null && a.info.rownames != null) {
        char ** S = cast(char **) malloc(a.rows * (char *).sizeof);

        if (S != null) {
            for (i=0; i<a.rows; i++) {
                S[i] = a.info.rownames[i];
            }
            for (i=0; i<a.rows; i++) {
                a.info.rownames[i] = S[rsx[i].row];
            }
            free(S);
        }
    }

 bailout:

    free(rsx);
    free(rvals);

    return a;
}

private bool has_colnames(const(GretlMatrixRaw) * m)
{
    return (m != null && !is_block_matrix(m) &&
            m.info != null && m.info.colnames != null);
}

private bool has_rownames(const(GretlMatrixRaw) * m)
{
    return (m != null && !is_block_matrix(m) &&
            m.info != null && m.info.rownames != null);
}

private struct NamedVal {
    double x;
    const(char) * s;
}

private NamedVal *
make_named_vals(const(GretlMatrixRaw) * m, char ** S, int n)
{
    NamedVal * nv = cast(NamedVal *) malloc(n * NamedVal.sizeof);

    if (nv != null) {
        int i;

        for (i=0; i<n; i++) {
            nv[i].x = m.val[i];
            nv[i].s = S[i];
        }
    }

    return nv;
}

private int vector_copy_marginal_names(GretlVector * v,
                                       NamedVal * nv,
                                       int n)
{
    int err = gretl_matrix_add_info(v);

    /* note: we assume v->info is null on entry */

    if (!err) {
        char *** pS;
        int i;

        pS = v.cols > 1 ? &v.info.colnames : &v.info.rownames;
        *pS = strings_array_new(n);
        if (*pS != null) {
            for (i=0; i<n; i++) {
                (*pS)[i] = gretl_strdup(nv[i].s);
            }
        } else {
            err = GretlError.alloc;
        }
    }

    return err;
}

GretlMatrixRaw * gretl_vector_sort(const(GretlMatrixRaw) * v,
                                int descending,
                                int * err)
{
    int n = gretl_vector_get_length(v);
    GretlMatrixRaw * vs = null;

    if (n == 0) {
        *err = GretlError.types;
        return null;
    }

    vs = matrix_copy_plain(v);

    if (vs == null) {
        *err = GretlError.alloc;
    } else {
        NamedVal * nvals = null;
        char ** S = null;

        if (v.cols > 1 && has_colnames(v)) {
            S = cast(char **) v.info.colnames;
        } else if (v.rows > 1 && has_rownames(v)) {
            S = cast(char **) v.info.rownames;
        }

        if (S != null) {
            nvals = make_named_vals(v, S, n);
            if (nvals == null) {
                *err = GretlError.alloc;
            }
        }

        if (nvals != null) {
            int i;

            qsort(nvals, n, (*nvals).sizeof, descending ?
                  &inverse_compare_values : &compare_values);
            for (i=0; i<n; i++) {
                vs.val[i] = nvals[i].x;
            }
            vector_copy_marginal_names(vs, nvals, n);
            free(nvals);
        } else if (!*err) {
            double * x = vs.val;

            qsort(x, n, (*x).sizeof, descending ? &gretl_inverse_compare_doubles :
                  &gretl_compare_doubles);
        }
    }

    return vs;
}

/* Calculate X(t)-transpose * X(t-lag) */

private void xtxlag(GretlMatrixRaw * wt, const(GretlMatrixRaw) * X,
                    int n, int t, int lag)
{
    double xi, xj;
    int i, j;

    for (i=0; i<n; i++) {
        xi = gretl_matrix_get(X, t, i);
        for (j=0; j<n; j++) {
            xj = gretl_matrix_get(X, t - lag, j);
            gretl_matrix_set(wt, i, j, xi * xj);
        }
    }
}

/*
 * gretl_matrix_covariogram:
 * @X: T x k matrix (typically containing regressors).
 * @u: T-vector (typically containing residuals), or null.
 * @w: (p+1)-vector of weights, or null.
 * @p: lag order >= 0.
 * @err: location to receive error code.
 *
 * Produces the matrix covariogram,
 *
 * \sum_{j=-p}^{p} \sum_j w_{|j|} (X_t' u_t u_{t-j} X_{t-j})
 *
 * If @u is not given the u terms are omitted, and if @w
 * is not given, all the weights are 1.0.
 *
 * Returns: the generated matrix, or null on failure.
 */

GretlMatrixRaw * gretl_matrix_covariogram(const(GretlMatrixRaw) * X,
                                       const(GretlMatrixRaw) * u,
                                       const(GretlMatrixRaw) * w,
                                       int p, int * err)
{
    GretlMatrixRaw * V;
    GretlMatrixRaw * G;
    GretlMatrixRaw * xtj;
    double uu;
    int j, k, t, T;

    if (gretl_is_null_matrix(X)) {
        return null;
    }

    if (gretl_is_complex(X) ||
        gretl_is_complex(u) ||
        gretl_is_complex(w)) {
        fprintf(stderr, "E_CMPLX in gretl_matrix_covariogram\n");
        *err = GretlError.cmplx;
        return null;
    }

    k = X.cols;
    T = X.rows;

    if (u != null && gretl_vector_get_length(u) != T) {
        *err = GretlError.nonconf;
        return null;
    }

    if (p < 0 || p > T) {
        *err = GretlError.nonconf;
        return null;
    }

    if (w != null && gretl_vector_get_length(w) != p + 1) {
        *err = GretlError.nonconf;
        return null;
    }

    V = gretl_zero_matrix_new(k, k);
    xtj = gretl_matrix_alloc(k, k);
    G = gretl_matrix_alloc(k, k);

    if (V == null || G == null || xtj == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    for (j=0; j<=p; j++) {
        gretl_matrix_zero(G);
        for (t=j; t<T; t++) {
            xtxlag(xtj, X, k, t, j);
            if (u != null) {
                uu = u.val[t] * u.val[t-j];
                gretl_matrix_multiply_by_scalar(xtj, uu);
            }
            gretl_matrix_add_to(G, xtj);
        }
        if (j > 0) {
            gretl_matrix_add_self_transpose(G);
        }
        if (w != null) {
            gretl_matrix_multiply_by_scalar(G, w.val[j]);
        }
        gretl_matrix_add_to(V, G);
    }

 bailout:

    gretl_matrix_free(G);
    gretl_matrix_free(xtj);

    if (*err) {
        gretl_matrix_free(V);
        V = null;
    }

    return V;
}

/*
 * gretl_matrix_GG_inverse:
 * @G: T x k source matrix.
 * @err: location to receive error code.
 *
 * Multiples G' into G and inverts the result. A shortcut
 * function intended for producing an approximation to
 * the Hessian given a gradient matrix.
 *
 * Returns: the newly allocated k x k inverse on success,
 * or null on error.
 */

GretlMatrixRaw * gretl_matrix_GG_inverse(const(GretlMatrixRaw) * G, int * err)
{
    GretlMatrixRaw * H = null;
    int k = G.cols;

    H = gretl_matrix_alloc(k, k);
    if (H == null) {
        *err = GretlError.alloc;
        return null;
    }

    gretl_matrix_multiply_mod(G, GretlMatrixMod.transpose,
                              G, GretlMatrixMod.none,
                              H, GretlMatrixMod.none);

    *err = gretl_invert_symmetric_matrix(H);

    if (*err) {
        fprintf(stderr, "gretl_matrix_GG_inverse: H not pd\n");
        gretl_matrix_free(H);
        H = null;
    }

    return H;
}

/*
 * gretl_matrix_commute:
 * @A: source matrix.
 * @r: row dimension.
 * @c: column dimension.
 * @pre: premultiply (Boolean flag).
 * @add_id: add identity matrix (Boolean flag).
 * @err: location to receive error code.
 *
 * It is assumed that @A is a matrix with (@r*@c) rows, so each of its
 * columns can be seen as the vectorization of an (@r x @c)
 * matrix. Each column of the output matrix contains the vectorization
 * of the transpose of the corresponding column of @A. This is
 * equivalent to premultiplying @A by the so-called "commutation
 * matrix" $K_{r,c}$. If the @add_id flag is non-zero, then @A is
 * added to the output matrix, so that @A is premultiplied by (I +
 * K_{r,c}) if @pre is nonzero, postmultiplied if @pre is 0.
 *
 * See eg Magnus and Neudecker (1988), "Matrix Differential Calculus
 * with Applications in Statistics and Econometrics".
 */

GretlMatrixRaw * gretl_matrix_commute(GretlMatrixRaw * A, int r, int c,
                                   int pre, int add_id, int * err)
{
    /* dim0 is the dimension on which the swapping has to happen; dim1
       is the other one */
    int dim0 = r * c;
    int dim1 = pre ? A.cols : A.rows;
    int * indices;
    GretlMatrixRaw * ret;

    /* dimension check */
    int dim_ok = pre ? (dim0 == A.rows) : (dim0 == A.cols);
    if (!dim_ok) {
        *err = GretlError.nonconf;
        return null;
    }

    indices = cast(int *) malloc(dim0 * (*indices).sizeof);
    if (indices == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (add_id) {
        ret = gretl_matrix_copy(A);
    } else {
        ret = gretl_zero_matrix_new(A.rows, A.cols);
    }

    if (ret == null) {
        *err = GretlError.alloc;
    } else {
        int i, j, h, k = 0;
        double x;

        for (i=0; i<r; i++) {
            for (j=0; j<c; j++) {
                indices[k++] = j*r + i;
            }
        }

        k = 0;
        if (pre) {
            for (j=0; j<dim1; j++) {
                for (i=0; i<dim0; i++) {
                    h = indices[i];
                    x = gretl_matrix_get(A, h, j);
                    ret.val[k++] += x;
                }
            }
        } else {
            for (j=0; j<dim0; j++) {
                h = indices[j];
                for (i=0; i<dim1; i++) {
                    x = gretl_matrix_get(A, i, h);
                    ret.val[k++] += x;
                }
            }
        }
    }

    free(indices);

    return ret;
}

/*
 * gretl_matrix_transcribe_obs_info:
 * @targ: target matrix.
 * @src: source matrix.
 *
 * If @targ and @src have the same number of rows, and if
 * the rows of @src are identified by observation stamps
 * while those of @targ are not so identified, copy the
 * stamp information across to @targ.  (Or if the given
 * conditions are not satified, do nothing.)
 */

void gretl_matrix_transcribe_obs_info(GretlMatrixRaw * targ,
                                      const(GretlMatrixRaw) * src)
{
    if (targ.rows == src.rows &&
        src.info != null && targ.info == null) {
        gretl_matrix_set_t1(targ, src.info.t1);
        gretl_matrix_set_t2(targ, src.info.t2);
    }
}

private GretlMatrixRaw * reorder_A(const(GretlMatrixRaw) * A,
                                int n, int np, int * err)
{
    GretlMatrixRaw * B;
    int p = np / n;

    B = gretl_matrix_alloc(np, n);

    if (B == null) {
        *err = GretlError.alloc;
        return null;
    } else {
        int i, j, k;
        int from, to;
        double x, y;

        for (j=0; j<n; j++) {
            for (k=0; k<=p/2; k++) {
                from = k*n;
                to = n*(p-k-1);
                for (i=0; i<n; i++) {
                    x = gretl_matrix_get(A, j, from + i);
                    y = gretl_matrix_get(A, j, to + i);
                    gretl_matrix_set(B, to + i, j, x);
                    gretl_matrix_set(B, from + i, j, y);
                }
            }
        }
    }

    return B;
}

/*
 * gretl_matrix_varsimul:
 * @A: n x np coefficient matrix.
 * @U: T x n data matrix.
 * @x0: p x n matrix for initialization.
 * @err: location to receive error code.
 *
 * Simulates a p-order n-variable VAR:
 * x_t = \sum A_i x_{t-i} + u_t
 *
 * The A_i matrices must be stacked horizontally into the @A
 * argument, that is: A = A_1 ~ A_2 ~ A_p. The u_t vectors are
 * contained (as rows) in @U. Initial values are in @x0.
 *
 * Note the that the arrangement of the @A matrix is somewhat
 * sub-optimal computationally, since its elements have to be
 * reordered by the function reorder_A (see above). However, the
 * present form is more intuitive for a human being, and that's
 * what counts.
 *
 * Returns: a newly allocated T+p x n matrix on success, whose t-th
 * row is (x_t)', or null on error.
 */

GretlMatrixRaw * gretl_matrix_varsimul(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * U,
                                    const(GretlMatrixRaw) * x0,
                                    int * err)
{
    GretlMatrixRaw * A2, X, UT;
    GretlVector xt, xtlag, ut;
    double x;
    int p = x0.rows;
    int n = x0.cols;
    int np = n * p;
    int T = p + U.rows;
    int t, i;

    if (A.rows != n || A.cols != np || U.cols != n) {
        *err = GretlError.nonconf;
        return null;
    }

    A2 = reorder_A(A, n, np, err);
    X = gretl_matrix_alloc(n, T);
    UT = gretl_matrix_copy_transpose(U);

    if (X == null || A2 == null || UT == null) {
        *err = GretlError.alloc;
        gretl_matrix_free(A2);
        gretl_matrix_free(X);
        gretl_matrix_free(UT);
        return null;
    }

    for (t=0; t<p; t++) {
        for (i=0; i<n; i++) {
            x = gretl_matrix_get(x0, t, i);
            gretl_matrix_set(X, i, t, x);
        }
    }

    gretl_matrix_init_full(&xt, 1, n, X.val + np);
    gretl_matrix_init_full(&ut, 1, n, UT.val);
    gretl_matrix_init_full(&xtlag, 1, np, X.val);

    for (t=p; t<T; t++) {
        gretl_matrix_multiply(&xtlag, A2, &xt);
        gretl_matrix_add_to(&xt, &ut);
        xt.val += n;
        xtlag.val += n;
        ut.val += n;
    }

    *err = gretl_matrix_transpose_in_place(X);

    if (!*err) {
        /* set dates on output matrix if possible */
        int t1 = gretl_matrix_get_t1(U) - p;

        if (t1 > 0) {
            gretl_matrix_set_t1(X, t1);
            gretl_matrix_set_t2(X, t1 + T - 1);
        }
    }

    gretl_matrix_free(A2);
    gretl_matrix_free(UT);

    return X;
}

/*
 * gretl_matrix_set_colnames:
 * @m: target matrix.
 * @S: array of strings.
 *
 * Sets an array of strings on @m which can be retrieved
 * using gretl_matrix_get_colnames(). Note that @S must
 * contain as many strings as @m has columns. The matrix
 * takes ownership of @S, which should be allocated and
 * not subsequently touched by the caller.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_matrix_set_colnames(GretlMatrixRaw * m, char ** S)
{
    if (m == null) {
        return GretlError.data;
    } else if (is_block_matrix(m)) {
        return matrix_block_error("gretl_matrix_set_colnames");
    } else if (S != null && m.info == null &&
               gretl_matrix_add_info(m)) {
        return GretlError.alloc;
    }

    if (m.info != null) {
        if (m.info.colnames != null) {
            strings_array_free(m.info.colnames, m.cols);
        }
        m.info.colnames = S;
    }

    return 0;
}

/*
 * gretl_matrix_set_rownames:
 * @m: target matrix.
 * @S: array of strings.
 *
 * Sets an array of strings on @m which can be retrieved
 * using gretl_matrix_get_rownames(). Note that @S must
 * contain as many strings as @m has rows. The matrix
 * takes ownership of @S, which should be allocated and
 * not subsequently touched by the caller.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_matrix_set_rownames(GretlMatrixRaw * m, char ** S)
{
    if (m == null) {
        return GretlError.data;
    } else if (is_block_matrix(m)) {
        return matrix_block_error("gretl_matrix_set_rownames");
    } else if (S != null && m.info == null &&
               gretl_matrix_add_info(m)) {
        return GretlError.alloc;
    }

    if (m.info != null) {
        if (m.info.rownames != null) {
            strings_array_free(m.info.rownames, m.rows);
        }
        m.info.rownames = S;
    }

    return 0;
}

/*
 * gretl_matrix_get_colnames:
 * @m: matrix
 *
 * Returns: The array of strings set on @m using
 * gretl_matrix_set_colnames(), or null if no such
 * strings have been set. The returned array will
 * contain as many strings as @m has columns.
 */

const(char) ** gretl_matrix_get_colnames(const(GretlMatrixRaw) * m)
{
    if (has_colnames(m)) {
        return cast(const(char) **) m.info.colnames;
    } else {
        return null;
    }
}

/*
 * gretl_matrix_get_rownames:
 * @m: matrix
 *
 * Returns:The array of strings set on @m using
 * gretl_matrix_set_rownames(), or null if no such
 * strings have been set. The returned array will
 * contain as many strings as @m has rows.
 */

const(char) ** gretl_matrix_get_rownames(const(GretlMatrixRaw) * m)
{
    if (has_rownames(m)) {
        return cast(const(char) **) m.info.rownames;
    } else {
        return null;
    }
}
