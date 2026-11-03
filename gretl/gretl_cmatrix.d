/* Port of lib/src/gretl_cmatrix.c: operations on complex matrices (the header gretl_cmatrix.h is
   gretl_cmatrixh.d), via LAPACK's complex routines and FFTW.

   Differences from the C:
   - C99 `double complex` is the Cmplx struct (libgretlh.d), so arithmetic is done with its operators; the
     elementwise complex functions (cpow, clog, cexp) are the C99 ones, re-declared with Cmplx and pragma(mangle);
   - the comparison operators of cmatrix_dot_op store 0/1 as complex values, through the helper cbool();
   - the function arguments of apply_cmatrix_dfunc/cfunc are the extern (C) function pointer types CDFunc and
     CCFunc;
   - FFTW is declared here with extern (C) (fftw_complex is Cmplx, an fftw_plan is an opaque pointer);
   - the `#if 0` code (gretl_zgemm_full) and the unused CDEC5 alternative are not ported.
   Dependencies not yet ported are declared extern (C) below, in the style of the other ports. */

module gretl.gretl_cmatrix;


import core.stdc.errno : errno;
import core.stdc.math : floor, fabs;
import core.stdc.stdio : fprintf, fputs, stderr;
import core.stdc.stdlib : free, malloc, realloc;
import core.stdc.string : memcpy, strlen;

import std.algorithm : max, min;

import gretl.clapack_complexh;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_errmsg_set_from_errno;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_f2ch : integer;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.gmem : g_free;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;

alias FftwPlan = void *;

extern (C) {
    void * fftw_malloc(size_t n);
    void fftw_free(void * p);
    FftwPlan fftw_plan_dft_r2c_1d(int n, double * in_, Cmplx * out_, uint flags);
    FftwPlan fftw_plan_dft_c2r_1d(int n, Cmplx * in_, double * out_, uint flags);
    FftwPlan fftw_plan_dft_1d(int n, Cmplx * in_, Cmplx * out_, int sign, uint flags);
    void fftw_execute(FftwPlan p);
    void fftw_destroy_plan(FftwPlan p);

    pragma(mangle, "cpow") Cmplx c_cpow(Cmplx x, Cmplx y);
    pragma(mangle, "clog") Cmplx c_clog(Cmplx z);
    pragma(mangle, "cexp") Cmplx c_cexp(Cmplx z);
}

enum FFTW_FORWARD = -1;
enum FFTW_BACKWARD = 1;
enum FFTW_ESTIMATE = 1U << 6;

bool cscalar(const(GretlMatrixRaw) * m)
{
    return m.rows == 1 && m.cols == 1;
}

bool cna(Cmplx z)
{
    return na(z.r) || na(z.i);
}

/* helper function for fftw-based real FFT functions */

int fft_allocate(double ** ffx, Cmplx ** ffz,
                 GretlMatrixRaw ** ret, int r, int c,
                 int inverse)
{
    /* real workspace, per series */
    *ffx = cast(double *) fftw_malloc(r * (**ffx).sizeof);
    if (*ffx == null) {
        return GretlError.alloc;
    }

    /* complex workspace, per series */
    *ffz = cast(Cmplx *) fftw_malloc((r/2 + 1 + r % 2) * (**ffz).sizeof);
    if (*ffz == null) {
        free(*ffx);
        return GretlError.alloc;
    }

    /* matrix to hold output */
    if (!inverse) {
        *ret = gretl_cmatrix_new(r, c);
    } else {
        *ret = gretl_matrix_alloc(r, c);
    }
    if (*ret == null) {
        free(*ffx);
        free(*ffz);
        return GretlError.alloc;
    }

    return 0;
}

/* FFT for real input -> complex output and FFTI for Hermetian input ->
   real output.  Both old and new-style complex formats are supported.
*/

GretlMatrixRaw * real_matrix_fft(const(GretlMatrixRaw) * y,
                              int inverse, int * err)
{
    GretlMatrixRaw * ret = null;
    FftwPlan p = null;
    double * ffx = null;
    Cmplx * ffz = null;
    int r, c, m, odd;
    int i, j;

    if (y.rows < 2) {
        *err = GretlError.data;
        return null;
    }

    r = y.rows;
    c = y.cols;
    m = r / 2;
    odd = r % 2;

    *err = fft_allocate(&ffx, &ffz, &ret, r, c, inverse);
    if (*err) {
        return null;
    }

    for (j=0; j<c; j++) {
        /* load the data */
        if (inverse) {
            for (i=0; i<=m+odd; i++) {
                ffz[i] = gretl_cmatrix_get(y, i, j);
            }
        } else {
            /* going in the real -> complex direction */
            for (i=0; i<r; i++) {
                ffx[i] = gretl_matrix_get(y, i, j);
            }
        }

        if (j == 0) {
            /* make the plan just once */
            if (inverse) {
                p = fftw_plan_dft_c2r_1d(r, ffz, ffx, FFTW_ESTIMATE);
            } else {
                p = fftw_plan_dft_r2c_1d(r, ffx, ffz, FFTW_ESTIMATE);
            }
        }

        /* run the transform */
        fftw_execute(p);

        /* transcribe the result */
        if (inverse) {
            for (i=0; i<r; i++) {
                gretl_matrix_set(ret, i, j, ffx[i] / r);
            }
        } else {
            for (i=0; i<=m+odd; i++) {
                gretl_cmatrix_set(ret, i, j, ffz[i]);
            }
            for (i=m; i>0; i--) {
                gretl_cmatrix_set(ret, r-i, j, conj(ffz[i]));
            }
        }
    }

    fftw_destroy_plan(p);
    fftw_free(ffz);
    fftw_free(ffx);

    return ret;
}

int row_is_real(const(GretlMatrixRaw) * y, int i)
{
    Cmplx z;
    int j;

    for (j=0; j<y.cols; j++) {
        z = gretl_cmatrix_get(y, i, j);
        if (z.i != 0) {
            return 0;
        }
    }

    return 1;
}

int rows_are_conjugate(const(GretlMatrixRaw) * y, int i, int k)
{
    Cmplx z1, z2;
    int j;

    for (j=0; j<y.cols; j++) {
        z1 = gretl_cmatrix_get(y, i, j);
        z2 = gretl_cmatrix_get(y, k, j);
        if (z1 != conj(z2)) {
            return 0;
        }
    }

    return 1;
}

int fft_is_hermitian(const(GretlMatrixRaw) * y)
{
    if (!row_is_real(y, 0)) {
        return 0;
    } else {
        int i, k, n2 = y.rows / 2;

        if (y.rows % 2 == 0 && !row_is_real(y, n2)) {
            /* In the Hermitian case the middle row of the
               odd-numbered remainder, on excluding the first
               row, must be real.
            */
            return 0;
        }
        for (i=1; i<n2; i++) {
            k = y.rows - i;
            if (!rows_are_conjugate(y, i, k)) {
                return 0;
            }
        }
    }

    return 1;
}

GretlMatrixRaw * gretl_matrix_fft(const(GretlMatrixRaw) * y, int * err)
{
    return real_matrix_fft(y, 0, err);
}

GretlMatrixRaw * gretl_matrix_ffti(const(GretlMatrixRaw) * y, int * err)
{
    if (!y.is_complex) {
        *err = GretlError.types;
        return null;
    } else if (fft_is_hermitian(y)) {
        /* its inverse should be real */
        return real_matrix_fft(y, 1, err);
    } else {
        return gretl_cmatrix_fft(y, 1, err);
    }
}

int cmatrix_validate(const(GretlMatrixRaw) * m, int square)
{
    int ret = 1;

    if (gretl_is_null_matrix(m)) {
        ret = 0;
    } else if (!m.is_complex || m.z == null) {
        ret = 0;
    } else if (square && m.rows != m.cols) {
        ret = 0;
    }

    if (!ret) {
        fprintf(stderr, "cmatrix_validate: failed\n".ptr);
    }

    return ret;
}

/* Construct a complex matrix, with all-zero imaginary part,
   from real matrix @A.
*/

GretlMatrixRaw * complex_from_real(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * C = null;

    if (gretl_is_null_matrix(A)) {
        *err = GretlError.data;
        return null;
    }

    C = gretl_cmatrix_new0(A.rows, A.cols);

    if (C == null) {
        *err = GretlError.alloc;
    } else {
        int i, n = A.rows * A.cols;

        for (i=0; i<n; i++) {
            C.z[i] = Cmplx(A.val[i], 0);
        }
    }

    return C;
}

void * oversize_error(int * err)
{
    *err = GretlError.invarg;
    gretl_errmsg_set("Too many elements");
    return null;
}

/* Variant of zgemm: simplified version of gretl_zgemm_full which
   allocates the product matrix, C.
*/

GretlMatrixRaw * gretl_zgemm(const(GretlMatrixRaw) * A,
                          char transa,
                          const(GretlMatrixRaw) * B,
                          char transb,
                          int * err)
{
    GretlMatrixRaw * C;
    Cmplx alpha = Cmplx(1, 0);
    Cmplx beta = Cmplx(0, 0);
    integer lda, ldb, ldc;
    integer m, n, k, kb;
    ulong chk;

    if (!cmatrix_validate(A, 0) || !cmatrix_validate(B, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    if (transa == 'C') {
        m = A.cols; /* rows of op(A) */
        k = A.rows; /* cols of op(A) */
    } else {
        m = A.rows; /* complex rows of A */
        k = A.cols; /* cols of A */
    }

    if (transb == 'C') {
        n = B.rows;  /* columns of op(B) */
        kb = B.cols; /* rows of op(B) */
    } else {
        n = B.cols;
        kb = B.rows;
    }

    if (k != kb) {
        *err = GretlError.nonconf;
        return null;
    }

    chk = m * cast(ulong) n;
    if (chk > cast(ulong) int.max) {
        return cast(GretlMatrixRaw *) oversize_error(err);
    }

    C = gretl_cmatrix_new(m, n);
    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    lda = A.rows;
    ldb = B.rows;
    ldc = C.rows;

    zgemm_(&transa, &transb, &m, &n, &k, &alpha,
           cast(Cmplx *) A.val, &lda, cast(Cmplx *) B.val, &ldb,
           &beta, cast(Cmplx *) C.val, &ldc);

    return C;
}

GretlMatrixRaw * real_cmatrix_multiply(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * B,
                                    char transa, int * err)
{
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * T = null;

    if (A.is_complex && B.is_complex) {
        if (transa == 'N' && (cscalar(A) || cscalar(B))) {
            return gretl_cmatrix_dot_op(A, B, '*', err);
        } else {
            C = gretl_zgemm(A, transa, B, 'N', err);
        }
    } else if (A.is_complex) {
        /* case of real B */
        T = complex_from_real(B, err);
        if (T != null) {
            C = gretl_zgemm(A, transa, T, 'N', err);
        }
    } else if (B.is_complex) {
        /* case of real A */
        T = complex_from_real(A, err);
        if (T != null) {
            C = gretl_zgemm(T, transa, B, 'N', err);
        }
    } else {
        *err = GretlError.types;
    }

    gretl_matrix_free(T);

    return C;
}

/* Multiplication of @A times @B, where we know that at least one of
   them is complex; the other, if it's not complex, must be converted to
   a complex matrix with a zero imaginary part first.
*/

GretlMatrixRaw * gretl_cmatrix_multiply(const(GretlMatrixRaw) * A,
                                     const(GretlMatrixRaw) * B,
                                     int * err)
{
    return real_cmatrix_multiply(A, B, 'N', err);
}

/* Returns (conjugate transpose of A, or A^H) times B, allowing for the
   possibility that either A or B (but not both!) may be a real matrix
   on input.
*/

GretlMatrixRaw * gretl_cmatrix_AHB(const(GretlMatrixRaw) * A,
                                const(GretlMatrixRaw) * B,
                                int * err)
{
    return real_cmatrix_multiply(A, B, 'C', err);
}

/* Eigen decomposition of complex (Hermitian) matrix using LAPACK's
   zheev()
*/

GretlMatrixRaw * gretl_zheev(GretlMatrixRaw * A, int eigenvecs, int * err)
{
    GretlMatrixRaw * evals = null;
    integer n, info, lwork;
    double * w = null;
    double * rwork = null;
    Cmplx * a = null;
    Cmplx * work = null;
    Cmplx wsz;
    char jobz = eigenvecs ? 'V' : 'N';
    char uplo = 'U';

    if (!cmatrix_validate(A, 1)) {
        *err = GretlError.invarg;
        return null;
    }

    n = A.rows;
    evals = gretl_matrix_alloc(n, 1);
    if (evals == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    w = evals.val;
    a = cast(Cmplx *) A.val;

    /* get optimal workspace size */
    lwork = -1;
    zheev_(&jobz, &uplo, &n, a, &n, w, &wsz, &lwork, rwork, &info);

    lwork = cast(integer) wsz.r;
    work = cast(Cmplx *) malloc(lwork * (*work).sizeof);
    rwork = cast(double *) malloc((3 * n - 2) * (*rwork).sizeof);
    if (work == null || rwork == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* do the actual decomposition */
    zheev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zheev: info = %d\n".ptr, info);
        *err = GretlError.data;
    }

 bailout:

    free(rwork);
    free(work);

    if (*err) {
        gretl_matrix_free(evals);
        evals = null;
    }

    return evals;
}

int ensure_aux_cmatrix(GretlMatrixRaw * m,
                       GretlMatrixRaw ** paux,
                       Cmplx ** ppz,
                       int r, int c)
{
    GretlMatrixRaw * aux = null;
    int mrc, dim = r * c;

    mrc = (m == null)? 0 : m.rows * m.cols;

    /* We need an r x c complex matrix for output:
       can we just use @m?
    */
    if (mrc == dim && m.is_complex) {
        /* OK, reusable complex matrix */
        m.rows = r;
        m.cols = c;
        if (ppz != null) {
            *ppz = cast(Cmplx *) m.val;
        }
    } else if (mrc == 2*dim && !m.is_complex) {
        /* OK, reusable real matrix */
        m.rows = 2*r;
        m.cols = c;
        gretl_matrix_set_complex_full(m, 1);
        if (ppz != null) {
            *ppz = cast(Cmplx *) m.val;
        }
    } else {
        /* nope, have to allocate anew */
        aux = gretl_cmatrix_new0(r, c);
        if (aux == null) {
            return GretlError.alloc;
        } else {
            *paux = aux;
            if (ppz != null) {
                *ppz = cast(Cmplx *) aux.val;
            }
        }
    }

    return 0;
}

int ensure_aux_matrix(GretlMatrixRaw * m,
                      GretlMatrixRaw ** paux,
                      int r, int c)
{
    GretlMatrixRaw * aux = null;
    int mrc, dim = r * c;

    mrc = (m == null)? 0 : m.rows * m.cols;

    if (mrc == dim && !m.is_complex) {
        /* OK, reusable real matrix */
        m.rows = r;
        m.cols = c;
    } else {
        /* have to allocate anew */
        aux = gretl_matrix_alloc(r, c);
        if (aux == null) {
            return GretlError.alloc;
        } else {
            *paux = aux;
        }
    }

    return 0;
}

/* Eigen decomposition of complex (non-Hermitian) matrix using
   LAPACK's zgeev() */

GretlMatrixRaw * gretl_zgeev(const(GretlMatrixRaw) * A,
                          GretlMatrixRaw * VR,
                          GretlMatrixRaw * VL,
                          int * err)
{
    GretlMatrixRaw * ret = null;
    GretlMatrixRaw * Acpy = null;
    GretlMatrixRaw * Ltmp = null;
    GretlMatrixRaw * Rtmp = null;
    integer n, info, lwork;
    integer ldvl, ldvr;
    double * w = null;
    double * rwork = null;
    Cmplx * a = null;
    Cmplx * work = null;
    Cmplx * vl = null;
    Cmplx * vr = null;
    Cmplx wsz;
    char jobvl = VL != null ? 'V' : 'N';
    char jobvr = VR != null ? 'V' : 'N';

    if (!cmatrix_validate(A, 1)) {
        *err = GretlError.invarg;
        return null;
    }

    n = A.rows;
    ldvl = VL != null ? n : 1;
    ldvr = VR != null ? n : 1;

    /* we need a copy of @A, which gets overwritten */
    Acpy = gretl_matrix_copy(A);
    if (Acpy == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    a = cast(Cmplx *) Acpy.val;

    if (VL != null) {
        /* left eigenvectors wanted */
        *err = ensure_aux_cmatrix(VL, &Ltmp, &vl, n, n);
        if (*err) {
            goto bailout;
        }
    }

    if (VR != null) {
        /* right eigenvectors wanted */
        *err = ensure_aux_cmatrix(VR, &Rtmp, &vr, n, n);
        if (*err) {
            goto bailout;
        }
    }

    rwork = cast(double *) malloc(2 * n * (*rwork).sizeof);
    ret = gretl_cmatrix_new(n, 1);
    if (rwork == null || ret == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    w = ret.val;

    /* get optimal workspace size */
    lwork = -1;
    zgeev_(&jobvl, &jobvr, &n, a, &n, w, vl, &ldvl, vr, &ldvr,
           &wsz, &lwork, rwork, &info);
    lwork = cast(integer) wsz.r;
    work = cast(Cmplx *) malloc(lwork * (*work).sizeof);
    if (work == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* do the actual decomposition */
    zgeev_(&jobvl, &jobvr, &n, a, &n, w, vl, &ldvl, vr, &ldvr,
           work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeev: info = %d\n".ptr, info);
        *err = GretlError.data;
    } else {
        if (Ltmp != null) {
            gretl_matrix_replace_content(VL, Ltmp);
        }
        if (Rtmp != null) {
            gretl_matrix_replace_content(VR, Rtmp);
        }
    }

 bailout:

    free(rwork);
    free(work);
    gretl_matrix_free(Acpy);
    gretl_matrix_free(Ltmp);
    gretl_matrix_free(Rtmp);

    if (*err) {
        gretl_matrix_free(ret);
        ret = null;
    }

    return ret;
}

/* Schur factorization of @A, with optional assignment of the
   matrix of Schur vectors to @Z and/or the eigenvalues of @A
   to @W, via LAPACK zgees().
*/

GretlMatrixRaw * gretl_zgees(const(GretlMatrixRaw) * A,
                          GretlMatrixRaw * Z,
                          GretlMatrixRaw * W,
                          int * err)
{
    GretlMatrixRaw * ret = null;
    GretlMatrixRaw * Ztmp = null;
    GretlMatrixRaw * Wtmp = null;
    integer n, info, lwork;
    integer ldvs;
    double * rwork = null;
    Cmplx * a = null;
    Cmplx * work = null;
    Cmplx * vs = null;
    Cmplx * w = null;
    char jobvs = Z != null ? 'V' : 'N';
    char srt = 'N';
    integer sdim = 0;

    if (!cmatrix_validate(A, 1)) {
        *err = GretlError.invarg;
        return null;
    }

    n = A.rows;
    ldvs = Z != null ? n : 1;

    /* we need a copy of @A, which gets overwritten */
    ret = gretl_matrix_copy(A);
    if (ret == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    a = ret.z;

    if (Z != null) {
        /* Schur vectors wanted */
        *err = ensure_aux_cmatrix(Z, &Ztmp, &vs, n, n);
        if (*err) {
            goto bailout;
        }
    }

    /* Eigenvalues vector: seems like we need this,
       regardless of whether @W is passed?
    */
    *err = ensure_aux_cmatrix(W, &Wtmp, &w, n, 1);
    if (*err) {
        goto bailout;
    }

    work = cast(Cmplx *) malloc((*work).sizeof);
    rwork = cast(double *) malloc(n * (*rwork).sizeof);
    if (work == null || rwork == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* get optimal workspace size */
    lwork = -1;
    zgees_(&jobvs, &srt, null, &n, a, &n, &sdim, w, vs, &ldvs,
           work, &lwork, rwork, null, &info);
    lwork = cast(integer) work[0].r;
    work = cast(Cmplx *) realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* do the actual factorization */
    zgees_(&jobvs, &srt, null, &n, a, &n, &sdim, w, vs, &ldvs,
           work, &lwork, rwork, null, &info);
    if (info != 0) {
        *err = GretlError.data;
    } else {
        if (Ztmp != null) {
            gretl_matrix_replace_content(Z, Ztmp);
        }
        if (W != null && Wtmp != null) {
            gretl_matrix_replace_content(W, Wtmp);
        }
    }

 bailout:

    free(rwork);
    free(work);
    gretl_matrix_free(Ztmp);
    gretl_matrix_free(Wtmp);

    if (*err) {
        gretl_matrix_free(ret);
        ret = null;
    }

    return ret;
}

/* SVD of a complex matrix via the LAPACK function zgesvd() */

int gretl_cmatrix_SVD(const(GretlMatrixRaw) * x, GretlMatrixRaw ** pu,
                      GretlVector ** ps, GretlMatrixRaw ** pvt,
                      int full)
{
    integer m, n, lda;
    integer ldu = 1, ldvt = 1;
    integer lwork = -1L;
    double * rwork = null;
    integer info;
    GretlMatrixRaw * s = null;
    GretlMatrixRaw * u = null;
    GretlMatrixRaw * vt = null;
    Cmplx * az;
    char jobu = 'N', jobvt = 'N';
    Cmplx zu, zvt;
    Cmplx * uval = &zu;
    Cmplx * vtval = &zvt;
    Cmplx * work = null;
    int xsize, k, err = 0;

    if (pu == null && ps == null && pvt == null) {
        /* no-op */
        return 0;
    }

    if (!cmatrix_validate(x, 0)) {
        return GretlError.invarg;
    }

    lda = m = x.rows;
    n = x.cols;
    xsize = lda * n;

    az = cast(Cmplx *) malloc(xsize * (*az).sizeof);
    if (az == null) {
        return GretlError.alloc;
    }
    memcpy(az, x.val, xsize * (*az).sizeof);

    k = min(m, n);

    s = gretl_vector_alloc(k);
    if (s == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (pu != null) {
        ldu = m;
        if (full) {
            u = gretl_cmatrix_new(ldu, m);
        } else {
            u = gretl_cmatrix_new(ldu, k);
        }
        if (u == null) {
            err = GretlError.alloc;
            goto bailout;
        } else {
            uval = cast(Cmplx *) u.val;
            jobu = full ? 'A' : 'S';
        }
    }

    if (pvt != null) {
        ldvt = full ? n : k;
        vt = gretl_cmatrix_new(ldvt, n);
        if (vt == null) {
            err = GretlError.alloc;
            goto bailout;
        } else {
            vtval = cast(Cmplx *) vt.val;
            jobvt = full ? 'A' : 'S';
        }
    }

    work = cast(Cmplx *) malloc((*work).sizeof);
    rwork = cast(double *) malloc(5 * k * (*rwork).sizeof);
    if (work == null || rwork == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* workspace query */
    lwork = -1;
    zgesvd_(&jobu, &jobvt, &m, &n, az, &lda, s.val, uval, &ldu,
            vtval, &ldvt, work, &lwork, rwork, &info);

    if (info != 0 || work[0].r <= 0.0) {
        fprintf(stderr, "zgesvd: workspace query failed\n".ptr);
        err = GretlError.data;
        goto bailout;
    }

    lwork = cast(integer) work[0].r;
    work = cast(Cmplx *) realloc(work, lwork * (*work).sizeof);
    if (work == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* actual computation */
    zgesvd_(&jobu, &jobvt, &m, &n, az, &lda, s.val, uval, &ldu,
            vtval, &ldvt, work, &lwork, rwork, &info);

    if (info != 0) {
        fprintf(stderr, "gretl_cmatrix_SVD: info = %d\n".ptr, cast(int) info);
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

    free(az);
    free(work);
    free(rwork);
    gretl_matrix_free(s);
    gretl_matrix_free(u);
    gretl_matrix_free(vt);

    return err;
}

double csvd_smin(const(GretlMatrixRaw) * a, double smax)
{
    const double macheps = 2.20e-16;
    int dmax = (a.rows > a.cols)? a.rows : a.cols;

    /* as per numpy, Matlab (2015-09-28) */
    return dmax * macheps * smax;
}

int gretl_cmatrix_rank(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * S = null;
    int i, k, rank = 0;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return 0;
    }

    k = min(A.rows, A.cols);

    if (A.rows > 4 * k || A.cols > 4 * k) {
        char trans1 = A.rows > k ? 'C' : 'N';
        char trans2 = A.cols > k ? 'C' : 'N';
        GretlMatrixRaw * B;

        B = gretl_zgemm(A, trans1, A, trans2, err);
        if (!*err) {
            *err = gretl_cmatrix_SVD(B, null, &S, null, 1);
        }
        gretl_matrix_free(B);
    } else {
        *err = gretl_cmatrix_SVD(A, null, &S, null, 1);
    }

    if (!*err) {
        double smin = csvd_smin(A, S.val[0]);

        for (i=0; i<k; i++) {
            if (S.val[i] > smin) {
                rank++;
            }
        }
    }

    gretl_matrix_free(S);

    return rank;
}

/* Inverse of a complex matrix via its SVD. If the  @generalized
   flag is non-zero the generalized inverse is produced in case
   of rank deficiency, othewise an error is flagged if @A is not
   of full rank.
*/

GretlMatrixRaw * cmatrix_SVD_inverse(const(GretlMatrixRaw) * A,
                                  int generalized,
                                  int * err)
{
    GretlMatrixRaw * U = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * s = null;
    GretlMatrixRaw * Vt = null;
    GretlMatrixRaw * ret = null;

    *err = gretl_cmatrix_SVD(A, &U, &s, &V, 1);

    if (!*err) {
        double smin = csvd_smin(A, s.val[0]);
        Cmplx vij;
        int i, j, h = 0;

        for (i=0; i<s.cols; i++) {
            h += s.val[i] > smin;
        }
        if (!generalized && h < s.cols) {
            gretl_errmsg_set(_("Matrix is singular").ptr);
            *err = GretlError.singular;
            goto bailout;
        }
        Vt = gretl_ctrans(V, 1, err);
        if (!*err) {
            for (j=0; j<h; j++) {
                for (i=0; i<Vt.rows; i++) {
                    vij = gretl_cmatrix_get(Vt, i, j);
                    gretl_cmatrix_set(Vt, i, j, vij / Cmplx(s.val[j], 0));
                }
            }
        }
        if (!*err) {
            Vt.cols = U.cols = h;
            ret = gretl_zgemm(Vt, 'N', U, 'C', err);
        }
    }

 bailout:

    gretl_matrix_free(U);
    gretl_matrix_free(V);
    gretl_matrix_free(s);
    gretl_matrix_free(Vt);

    return ret;
}

GretlMatrixRaw * gretl_cmatrix_ginv(const(GretlMatrixRaw) * A, int * err)
{
    return cmatrix_SVD_inverse(A, 1, err);
}

GretlMatrixRaw * gretl_cmatrix_inverse(const(GretlMatrixRaw) * A, int * err)
{
    return cmatrix_SVD_inverse(A, 0, err);
}

/* Horizontal direct product of complex @A and @B, or of
   @A with itself if B is NULL.
*/

GretlMatrixRaw * real_cmatrix_hdp(const(GretlMatrixRaw) * A,
                               const(GretlMatrixRaw) * B,
                               int * err)
{
    GretlMatrixRaw * C = null;
    Cmplx aij, bik;
    int do_symmetric = 0;
    ulong chk;
    int i, j, k, ndx;
    int r, p, q, ccols;

    if (!cmatrix_validate(A,0)) {
        *err = GretlError.invarg;
    } else if (B != null) {
        if (!cmatrix_validate(B,0)) {
            *err = GretlError.invarg;
        } else if (B.rows != A.rows) {
            *err = GretlError.nonconf;
        }
    } else {
        do_symmetric = 1;
    }

    if (*err) {
        return null;
    }

    r = A.rows;
    p = A.cols;

    if (do_symmetric) {
        q = p;
        chk = cast(ulong) p * (p+1) / 2;
        if (chk > cast(ulong) int.max) {
            return cast(GretlMatrixRaw *) oversize_error(err);
        }
        ccols = p * (p+1) / 2;
    } else {
        q = B.cols;
        chk = cast(ulong) p * q;
        if (chk > cast(ulong) int.max) {
            return cast(GretlMatrixRaw *) oversize_error(err);
        }
        ccols = p * q;
    }

    chk = r * cast(ulong) ccols;
    if (chk > cast(ulong) int.max) {
        return cast(GretlMatrixRaw *) oversize_error(err);
    }

    C = gretl_cmatrix_new0(r, ccols);

    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=0; i<r; i++) {
        ndx = 0;
        for (j=0; j<p; j++) {
            aij = gretl_cmatrix_get(A, i, j);
            if (do_symmetric) {
                for (k=j; k<q; k++) {
                    bik = gretl_cmatrix_get(A, i, k);
                    gretl_cmatrix_set(C, i, ndx++, aij*conj(bik));
                }
            } else if (aij != Cmplx(0, 0)) {
                ndx = j * q;
                for (k=0; k<q; k++) {
                    bik = gretl_cmatrix_get(B, i, k);
                    gretl_cmatrix_set(C, i, ndx + k, aij*bik);
                }
            }
        }
    }

    return C;
}

/* Kronecker product of complex @A and @B */

GretlMatrixRaw * real_cmatrix_kron(const(GretlMatrixRaw) * A,
                                const(GretlMatrixRaw) * B,
                                int * err)
{
    GretlMatrixRaw * K = null;
    ulong chk;
    int p, q, r, s;

    if (!cmatrix_validate(A,0) || !cmatrix_validate(B,0)) {
        *err = GretlError.invarg;
        return null;
    }

    p = A.rows;
    q = A.cols;
    r = B.rows;
    s = B.cols;

    chk = cast(ulong) p * r * q * s;
    if (chk > cast(ulong) int.max) {
        return cast(GretlMatrixRaw *) oversize_error(err);
    }

    K = gretl_cmatrix_new0(p * r, q * s);

    if (K == null) {
        *err = GretlError.alloc;
    } else {
        Cmplx x, aij, bkl;
        int i, j, k, l;
        int ioff, joff;
        int Ki, Kj;

        for (i=0; i<p; i++) {
            ioff = i * r;
            for (j=0; j<q; j++) {
                /* block ij is an r * s matrix, a_{ij} * B */
                aij = gretl_cmatrix_get(A, i, j);
                joff = j * s;
                for (k=0; k<r; k++) {
                    Ki = ioff + k;
                    for (l=0; l<s; l++) {
                        bkl = gretl_cmatrix_get(B, k, l);
                        Kj = joff + l;
                        x = aij * bkl;
                        gretl_cmatrix_set(K, Ki, Kj, x);
                    }
                }
            }
        }
    }

    return K;
}

GretlMatrixRaw * gretl_cmatrix_kronlike(const(GretlMatrixRaw) * A,
                                     const(GretlMatrixRaw) * B,
                                     int hdp, int * err)
{
    GretlMatrixRaw * L = cast(GretlMatrixRaw *) A;
    GretlMatrixRaw * R = cast(GretlMatrixRaw *) B;
    GretlMatrixRaw * C = null;

    if (A.is_complex && hdp && B == null) {
        ; /* OK */
    } else if (A.is_complex && B.is_complex) {
        ; /* OK */
    } else if (A.is_complex) {
        R = complex_from_real(B, err);
    } else if (B.is_complex) {
        L = complex_from_real(A, err);
    } else {
        *err = GretlError.types;
    }

    if (!*err) {
        if (hdp) {
            C = real_cmatrix_hdp(L, R, err);
        } else {
            C = real_cmatrix_kron(L, R, err);
        }
    }

    if (L != A) gretl_matrix_free(L);
    if (R != B) gretl_matrix_free(R);

    return C;
}

GretlMatrixRaw * gretl_cmatrix_hdprod(const(GretlMatrixRaw) * A,
                                   const(GretlMatrixRaw) * B,
                                   int * err)
{
    return gretl_cmatrix_kronlike(A, B, 1, err);
}

GretlMatrixRaw * gretl_cmatrix_kronecker(const(GretlMatrixRaw) * A,
                                      const(GretlMatrixRaw) * B,
                                      int * err)
{
    return gretl_cmatrix_kronlike(A, B, 0, err);
}

/* Complex FFT (or inverse) via fftw */

GretlMatrixRaw * gretl_cmatrix_fft(const(GretlMatrixRaw) * A,
                                int inverse, int * err)
{
    GretlMatrixRaw * B = null;
    Cmplx * tmp;
    Cmplx * ptr;
    FftwPlan p;
    int sign;
    int r, c, j;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    B = gretl_matrix_copy(A);
    if (B == null) {
        *err = GretlError.alloc;
        return null;
    }

    r = A.rows;
    c = A.cols;

    tmp = cast(Cmplx *) B.val;
    sign = inverse ? FFTW_BACKWARD : FFTW_FORWARD;

    ptr = tmp;
    for (j=0; j<c; j++) {
        p = fftw_plan_dft_1d(r, ptr, ptr, sign, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);
        /* advance pointer to next column */
        ptr += r;
    }

    if (inverse) {
        /* "FFTW computes an unnormalized transform: computing a
            forward followed by a backward transform (or vice versa)
            will result in the original data multiplied by the size of
            the transform (the product of the dimensions)."
            So should we do the following?
        */
        for (j=0; j<r*c; j++) {
            tmp[j] /= Cmplx(r, 0);
        }
    }

    return B;
}

enum CDEC5 = true; /* five decimal places in default format */

int gretl_cmatrix_print_range(const(GretlMatrixRaw) * A,
                              const(char) * name,
                              int rmin, int rmax,
                              Prn * prn)
{
    Cmplx aij;
    double re, im, ai, xmax;
    int alt_default = 0;
    int all_ints = 1;
    int intmax = 1000;
    int altmin = 1000;
    int zwidth = 0;
    int rdecs, idecs;
    char pm;
    int r, c, i, j;

    if (!cmatrix_validate(A, 0)) {
        return GretlError.invarg;
    }

    r = A.rows;
    c = A.cols;

    if (rmin < 0) rmin = 0;
    if (rmax < 0) rmax = r;

    xmax = 0;

    for (j=0; j<c; j++) {
        for (i=rmin; i<rmax; i++) {
            aij = gretl_cmatrix_get(A, i, j);
            re = aij.r;
            im = aij.i;
            if (all_ints && (floor(re) != re || floor(im) != im)) {
                all_ints = 0;
            }
            re = max(fabs(re), fabs(im));
            if (re > xmax) {
                xmax = re;
            }
            if (all_ints && xmax >= intmax) {
                all_ints = 0;
            }
            if (!all_ints && xmax >= altmin) {
                break;
            }
        }
        if (!all_ints && xmax >= altmin) {
            break;
        }
    }

    if (all_ints) {
        /* apply a more compact format */
        zwidth = 2 + (xmax >= 10) + (xmax >= 100);
    } else if (xmax >= altmin) {
        /* we'll want a different default format */
        alt_default = 1;
    }

    if (name != null && *name != '\0') {
        if (r > 1 || c > 1) {
            pprintf(prn, "%s (%d x %d)".ptr, name, r, c);
            pputs(prn, "\n\n");
        } else {
            pprintf(prn, "%s: ".ptr, name);
        }
    }

    for (i=rmin; i<rmax; i++) {
        for (j=0; j<c; j++) {
            aij = gretl_cmatrix_get(A, i, j);
            re = aij.r;
            im = aij.i;
            pm = (im >= -0) ? '+' : '-';
            ai = fabs(im);
            if (zwidth > 0) {
                pprintf(prn, "%*g %c %*gi".ptr, zwidth, re, pm, zwidth-1, ai);
            } else if (alt_default) {
                pprintf(prn, "%# 9.4e %c %#8.4ei".ptr, re, pm, ai);
            } else {
                rdecs = re <= -10 ? 4 : 5;
                idecs = ai >= 10 ? 4 : 5;
                pprintf(prn, "%8.*f %c %7.*fi".ptr, rdecs, re, pm, idecs, ai);
            }
            if (j < c - 1) {
                pputs(prn, "  ");
            }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    return 0;
}

int gretl_cmatrix_print(const(GretlMatrixRaw) * A,
                        const(char) * name,
                        Prn * prn)
{
    return gretl_cmatrix_print_range(A, name, -1, -1, prn);
}

int gretl_cmatrix_printf(const(GretlMatrixRaw) * A,
                         const(char) * fmt,
                         Prn * prn)
{
    Cmplx aij;
    double re, im;
    char * fmtstr = null;
    char[3] s = [' ', ' ', 0];
    int r, c, i, j;

    if (!cmatrix_validate(A, 0)) {
        return GretlError.invarg;
    }

    if (fmt == null) {
        fmt = "%7.4f";
    } else {
        /* we only accept floating-point formats */
        char lc = fmt[strlen(fmt) - 1];

        if (lc != 'f' && lc != 'g') {
            return GretlError.invarg;
        }
    }

    fmtstr = g_strdup_printf("%s%%s%si".ptr, fmt, fmt);

    r = A.rows;
    c = A.cols;

    for (i=0; i<r; i++) {
        for (j=0; j<c; j++) {
            aij = gretl_cmatrix_get(A, i, j);
            re = aij.r;
            im = aij.i;
            s[1] = (im >= 0) ? '+' : '-';
            pprintf(prn, fmtstr, re, s.ptr, fabs(im));
            if (j < c - 1) {
                pputs(prn, "  ");
            }
        }
        pputc(prn, '\n');
    }
    pputc(prn, '\n');

    g_free(fmtstr);

    return 0;
}

/* Compose a complex matrix from its real and imaginary
   components. If @Re is NULL the result will have a
   constant real part given by @x; and if @Im is NULL it
   will have a constant imaginary part given by @y. If
   both matrix arguments are non-NULL they must be of
   the same dimensions.
*/

GretlMatrixRaw * gretl_cmatrix_build(const(GretlMatrixRaw) * Re,
                                  const(GretlMatrixRaw) * Im,
                                  double x, double y,
                                  int * err)
{
    GretlMatrixRaw * C = null;
    int r = 1, c = 1;
    int i, n;

    if (Re != null) {
        if (Re.is_complex) {
            *err = GretlError.invarg;
        } else {
            r = Re.rows;
            c = Re.cols;
        }
    }
    if (!*err && Im != null) {
        if (Im.is_complex) {
            *err = GretlError.invarg;
        } else if (Re != null) {
            if (Im.rows != r || Im.cols != c) {
                *err = GretlError.nonconf;
            }
        } else {
            r = Im.rows;
            c = Im.cols;
        }
    }

    if (!*err) {
        n = r * c;
        C = gretl_cmatrix_new(r, c);
        if (C == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        for (i=0; i<n; i++) {
            if (Re != null) {
                x = Re.val[i];
            }
            if (Im != null) {
                y = Im.val[i];
            }
            C.z[i] = Cmplx(x, y);
        }
    }

    return C;
}

/* Extract the real part of @A if @im = 0 or the imaginary part
   if @im is non-zero.
*/

GretlMatrixRaw * gretl_cmatrix_extract(const(GretlMatrixRaw) * A,
                                    int im, int * err)
{
    GretlMatrixRaw * B = null;
    int i, n;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    B = gretl_matrix_alloc(A.rows, A.cols);

    if (B == null) {
        *err = GretlError.alloc;
    } else {
        n = A.rows * A.cols;
        for (i=0; i<n; i++) {
            B.val[i] = im ? A.z[i].i : A.z[i].r;
        }
    }

    return B;
}

/* Return [conjugate] transpose of complex matrix @A */

GretlMatrixRaw * gretl_ctrans(const(GretlMatrixRaw) * A,
                           int conjugate, int * err)
{
    GretlMatrixRaw * C = null;
    Cmplx aij;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    C = gretl_cmatrix_new(A.cols, A.rows);
    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (j=0; j<A.cols; j++) {
        for (i=0; i<A.rows; i++) {
            aij = gretl_cmatrix_get(A, i, j);
            if (conjugate) {
                aij = conj(aij);
            }
            gretl_cmatrix_set(C, j, i, aij);
        }
    }

    return C;
}

/* Convert complex matrix @A to its conjugate transpose */

int gretl_ctrans_in_place(GretlMatrixRaw * A)
{
    GretlMatrixRaw * C = null;
    Cmplx aij;
    int i, j, n;
    int err = 0;

    if (!cmatrix_validate(A, 0)) {
        return GretlError.invarg;
    }

    /* temporary matrix */
    C = gretl_cmatrix_new(A.cols, A.rows);
    if (C == null) {
        return GretlError.alloc;
    }

    for (j=0; j<A.cols; j++) {
        for (i=0; i<A.rows; i++) {
            aij = gretl_cmatrix_get(A, i, j);
            gretl_cmatrix_set(C, j, i, conj(aij));
        }
    }

    /* now rejig @A */
    A.rows = C.rows;
    A.cols = C.cols;
    n = C.rows * C.cols;
    memcpy(A.z, C.z, n * (*A.z).sizeof);
    gretl_matrix_destroy_info(A);

    gretl_matrix_free(C);

    return err;
}

/* Addition or subtraction of matrices: handle the case
   where one operand is complex and the other is real.
   In addition handle the case where one of them is a
   scalar matrix, real or complex.
*/

GretlMatrixRaw * gretl_cmatrix_add_sub(const(GretlMatrixRaw) * A,
                                    const(GretlMatrixRaw) * B,
                                    int sgn, int * err)
{
    GretlMatrixRaw * C = null;
    Cmplx aval = Cmplx(0, 0);
    Cmplx bval = Cmplx(0, 0);
    int cr = A.rows;
    int cc = A.cols;
    int a_scalar = 0;
    int b_scalar = 0;

    if (A.is_complex && B.rows == 1 && B.cols == 1) {
        b_scalar = 1;
        bval = B.is_complex ? B.z[0] : Cmplx(B.val[0], 0);
        goto allocate;
    } else if (B.is_complex && A.rows == 1 && A.cols == 1) {
        cr = B.rows;
        cc = B.cols;
        a_scalar = 1;
        aval = A.is_complex ? A.z[0] : Cmplx(A.val[0], 0);
        goto allocate;
    }

    if (B.cols != A.cols) {
        *err = GretlError.nonconf;
        return null;
    }

    if (A.is_complex && B.is_complex) {
        /* both complex */
        if (B.rows != cr) {
            *err = GretlError.nonconf;
        }
    } else if (A.is_complex) {
        /* A complex, B real */
        if (B.rows != cr) {
            *err = GretlError.nonconf;
        }
    } else {
        /* A real, B complex */
        cr = B.rows;
        if (A.rows != cr) {
            *err = GretlError.nonconf;
        }
    }

 allocate:

    if (!*err) {
        C = gretl_cmatrix_new(cr, cc);
        if (C == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        int i, n = cc * cr;

        if (b_scalar) {
            for (i=0; i<n; i++) {
                C.z[i] = sgn < 0 ? A.z[i] - bval : A.z[i] + bval;
            }
        } else if (a_scalar) {
            for (i=0; i<n; i++) {
                C.z[i] = sgn < 0 ? aval - B.z[i] : aval + B.z[i];
            }
        } else if (A.is_complex && B.is_complex) {
            for (i=0; i<n; i++) {
                C.z[i] = sgn < 0 ? A.z[i] - B.z[i] : A.z[i] + B.z[i];
            }
        } else if (A.is_complex) {
            for (i=0; i<n; i++) {
                Cmplx bi = Cmplx(B.val[i], 0);

                C.z[i] = sgn < 0 ? A.z[i] - bi : A.z[i] + bi;
            }
        } else {
            for (i=0; i<n; i++) {
                Cmplx ai = Cmplx(A.val[i], 0);

                C.z[i] = sgn < 0 ? ai - B.z[i] : ai + B.z[i];
            }
        }
    }

    return C;
}

double gretl_cquad(Cmplx z)
{
    double zr = z.r;
    double zi = z.i;

    return zr*zr + zi*zi;
}

/* function pointer types for the functions that are applied
   elementwise to complex matrices */

alias CDFunc = extern (C) double function(Cmplx);
alias CCFunc = extern (C) Cmplx function(Cmplx);

/* Apply a function which maps from complex to real:
   creal, cimag, carg, cmod.
*/

int apply_cmatrix_dfunc(GretlMatrixRaw * targ,
                        const(GretlMatrixRaw) * src,
                        CDFunc dfunc)
{
    int n = src.cols * src.rows;
    int i, err = 0;

    if (!cmatrix_validate(src, 0) || targ.is_complex) {
        return GretlError.invarg;
    }

    errno = 0;

    for (i=0; i<n; i++) {
        targ.val[i] = dfunc(src.z[i]);
    }

    if (errno) {
        gretl_errmsg_set_from_errno(null, errno);
        err = GretlError.data;
    }

    return err;
}

/* Apply a function which maps from complex to complex;
   includes trigonometric functions, log, exp.
*/

int apply_cmatrix_cfunc(GretlMatrixRaw * targ,
                        const(GretlMatrixRaw) * src,
                        CCFunc cfunc)
{
    int n = src.cols * src.rows;
    int i, err = 0;

    if (!cmatrix_validate(src, 0) || !cmatrix_validate(targ, 0)) {
        return GretlError.invarg;
    }

    errno = 0;

    for (i=0; i<n; i++) {
        targ.z[i] = cfunc(src.z[i]);
    }

    if (errno) {
        gretl_errmsg_set_from_errno(null, errno);
        err = GretlError.data;
    }

    return err;
}

int apply_cmatrix_unary_op(GretlMatrixRaw * targ,
                           const(GretlMatrixRaw) * src,
                           int op)
{
    int i, n = src.cols * src.rows;
    int err = 0;

    if (!cmatrix_validate(src, 0) || !cmatrix_validate(targ, 0)) {
        return GretlError.invarg;
    }

    for (i=0; i<n && !err; i++) {
        if (op == 1) {
            /* U_NEG */
            targ.z[i] = Cmplx(-src.z[i].r, -src.z[i].i);
        } else if (op == 2) {
            /* U_POS */
            targ.z[i] = src.z[i];
        } else if (op == 3) {
            /* U_NOT */
            targ.z[i] = Cmplx(src.z[i] == Cmplx(0, 0), 0);
        } else {
            err = GretlError.invarg;
        }
    }

    return err;
}

GretlMatrixRaw * gretl_cmatrix_from_scalar(Cmplx z, int * err)
{
    GretlMatrixRaw * ret = gretl_cmatrix_new(1, 1);

    if (ret == null) {
        *err = GretlError.alloc;
    } else {
        ret.z[0] = z;
    }

    return ret;
}

/* Determinant via eigenvalues */

GretlMatrixRaw * gretl_cmatrix_determinant(const(GretlMatrixRaw) * X,
                                        int ldet, int * err)
{
    GretlMatrixRaw * E = null;
    GretlMatrixRaw * ret = null;

    if (!cmatrix_validate(X, 1)) {
        *err = GretlError.invarg;
        return ret;
    }

    E = gretl_zgeev(X, null, null, err);

    if (E != null) {
        Cmplx cret = Cmplx(1, 0);
        int i;

        for (i=0; i<E.rows; i++) {
            cret *= E.z[i];
        }
        gretl_matrix_free(E);
        if (ldet) {
            cret = c_clog(cret);
        }
        ret = gretl_cmatrix_from_scalar(cret, err);
    }

    return ret;
}

GretlMatrixRaw * gretl_cmatrix_trace(const(GretlMatrixRaw) * X, int * err)
{
    GretlMatrixRaw * ret = null;

    if (!cmatrix_validate(X, 1)) {
        *err = GretlError.invarg;
    } else {
        Cmplx tr = Cmplx(0, 0);
        int i;

        for (i=0; i<X.rows; i++) {
            tr += gretl_cmatrix_get(X, i, i);
        }
        ret = gretl_cmatrix_from_scalar(tr, err);
    }

    return ret;
}

/* Set the diagonal of complex matrix @targ using either
   @src (if not NULL) or @x. In the first case @src can
   be either a complex vector of the right length, or a
   real vector, or a complex scalar.
*/

int gretl_cmatrix_set_diagonal(GretlMatrixRaw * targ,
                               const(GretlMatrixRaw) * src,
                               double x)
{
    Cmplx zi = Cmplx(0, 0);
    int d, i;
    int match = 0;
    int err = 0;

    if (!cmatrix_validate(targ, 0)) {
        return GretlError.invarg;
    }

    d = min(targ.rows, targ.cols);

    if (src != null) {
        if (src.is_complex) {
            if (gretl_vector_get_length(src) == d) {
                /* conformable complex vector */
                match = 1;
            } else if (cscalar(src)) {
                /* complex scalar */
                zi = src.z[0];
                match = 2;
            }
        } else if (gretl_vector_get_length(src) == d) {
            /* conformable real vector */
            match = 3;
        }
    } else {
        /* use real scalar, @x */
        zi = Cmplx(x, 0);
        match = 4;
    }

    if (match == 0) {
        return GretlError.nonconf;
    }

    for (i=0; i<d; i++) {
        if (match == 1) {
            gretl_cmatrix_set(targ, i, i, src.z[i]);
        } else if (match == 3) {
            gretl_cmatrix_set(targ, i, i, Cmplx(src.val[i], 0));
        } else {
            gretl_cmatrix_set(targ, i, i, zi);
        }
    }

    return err;
}

/* Set the lower or upper part of square complex matrix
   @targ using either @src (if not NULL) or @x. In the first
   case @src can be either a complex vector of the right length,
   or a real vector, or a complex scalar.
*/

int gretl_cmatrix_set_triangle(GretlMatrixRaw * targ,
                               const(GretlMatrixRaw) * src,
                               double x, int upper)
{
    Cmplx zi = Cmplx(0, 0);
    int r, c, p, i, j, n;
    int lower = !upper;
    int match = 0;
    int err = 0;

    if (!cmatrix_validate(targ, 0)) {
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
        if (src.is_complex) {
            if (gretl_vector_get_length(src) == n) {
                /* conformable complex vector */
                match = 1;
            } else if (cscalar(src)) {
                /* complex scalar */
                zi = src.z[0];
                match = 2;
            }
        } else if (gretl_vector_get_length(src) == n) {
            /* conformable real vector */
            match = 3;
        }
    } else {
        /* use real scalar, @x */
        zi = Cmplx(x, 0);
        match = 2;
    }

    if (match == 0) {
        return GretlError.nonconf;
    } else {
        int jmin = upper ? 1 : 0;
        int jmax = upper ? c : r;
        int imin = upper ? 0 : 1;
        int imax = upper ? 1 : r;
        int k = 0;

        for (j=jmin; j<jmax; j++) {
            for (i=imin; i<imax; i++) {
                if (match == 1) {
                    gretl_cmatrix_set(targ, i, j, src.z[k++]);
                } else if (match == 3) {
                    gretl_cmatrix_set(targ, i, j, Cmplx(src.val[k++], 0));
                } else {
                    gretl_cmatrix_set(targ, i, j, zi);
                }
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

/* switch between "legacy" and new representations of a
   complex matrix */

GretlMatrixRaw * gretl_cmatrix_switch(const(GretlMatrixRaw) * m,
                                   int to_new, int * err)
{
    GretlMatrixRaw * ret = null;
    int r, c, i, j, k;

    if (gretl_is_null_matrix(m)) {
        return gretl_null_matrix_new();
    }

    if ((to_new && m.is_complex) ||
        (!to_new && !cmatrix_validate(m, 0))) {
        *err = GretlError.invarg;
        return null;
    }

    r = m.rows;
    c = to_new ? m.cols / 2 : m.cols * 2;

    if (to_new) {
        ret = gretl_cmatrix_new(r, c);
    } else {
        ret = gretl_matrix_alloc(r, c);
    }
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    k = 0;

    if (to_new) {
        /* make re and im components contiguous */
        double xr, xi;

        for (j=0; j<ret.cols; j++) {
            for (i=0; i<m.rows; i++) {
                xr = gretl_matrix_get(m, i, k);
                xi = gretl_matrix_get(m, i, k+1);
                gretl_cmatrix_set(ret, i, j, Cmplx(xr, xi));
            }
            k += 2;
        }
    } else {
        /* put re and im components in adjacent columns */
        Cmplx mij;

        for (j=0; j<m.cols; j++) {
            for (i=0; i<ret.rows; i++) {
                mij = gretl_cmatrix_get(m, i, j);
                gretl_matrix_set(ret, i, k, mij.r);
                gretl_matrix_set(ret, i, k+1, mij.i);
            }
            k += 2;
        }
    }

    return ret;
}

/* Compute column or row sums, means or products */

GretlMatrixRaw * gretl_cmatrix_vector_stat(const(GretlMatrixRaw) * m,
                                        GretlVecStat vs, int rowwise,
                                        int skip_na, int * err)
{
    Cmplx z;
    Cmplx mij;
    GretlMatrixRaw * ret;
    int product;
    int r, c, i, j;

    if (!cmatrix_validate(m, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    r = rowwise ? m.rows : 1;
    c = rowwise ? 1 : m.cols;

    ret = gretl_cmatrix_new(r, c);
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    product = (vs == GretlVecStat.prod);

    if (rowwise) {
        /* by rows */
        Cmplx z0 = product ? Cmplx(1, 0) : Cmplx(0, 0);
        int ncols;

        for (i=0; i<m.rows; i++) {
            ncols = m.cols;
            z = z0;
            for (j=0; j<m.cols; j++) {
                mij = gretl_cmatrix_get(m, i, j);
                if (cna(mij) && skip_na) {
                    ncols--;
                    continue;
                } else if (product) {
                    z *= mij;
                } else {
                    z += mij;
                }
            }
            if (ncols == 0) {
                z = Cmplx(NADBL, 0);
            } else if (vs == GretlVecStat.mean) {
                z /= Cmplx(ncols, 0);
            }
            gretl_cmatrix_set(ret, i, 0, z);
        }
    } else {
        /* by columns */
        Cmplx z0 = product ? Cmplx(1, 0) : Cmplx(0, 0);
        int nrows;

        for (j=0; j<m.cols; j++) {
            nrows = m.rows;
            z = z0;
            for (i=0; i<m.rows; i++) {
                mij = gretl_cmatrix_get(m, i, j);
                if (cna(mij) && skip_na) {
                    nrows--;
                    continue;
                } else if (product) {
                    z *= mij;
                } else {
                    z += mij;
                }
            }
            if (nrows == 0) {
                z = Cmplx(NADBL, 0);
            } else if (vs == GretlVecStat.mean) {
                z /= Cmplx(nrows, 0);
            }
            gretl_cmatrix_set(ret, 0, j, z);
        }
    }

    return ret;
}

int gretl_cmatrix_fill(GretlMatrixRaw * m, Cmplx z)
{
    if (!cmatrix_validate(m, 0)) {
        return GretlError.types;
    } else {
        int i, n = m.cols * m.rows;

        for (i=0; i<n; i++) {
            m.z[i] = z;
        }
        return 0;
    }
}

/* the boolean results of the comparison operators are
   stored as complex 0/1 values */

Cmplx cbool(bool b)
{
    return Cmplx(b ? 1 : 0, 0);
}

void vec_x_op_vec_y(Cmplx * z,
                    const(Cmplx) * x,
                    const(Cmplx) * y,
                    int n, int op)
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
            z[i] = c_cpow(x[i], y[i]);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = cbool(x[i] == y[i]);
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = cbool(x[i] != y[i]);
        }
        break;
    default:
        break;
    }
}

void vec_x_op_y(Cmplx * z,
                const(Cmplx) * x,
                Cmplx y,
                int n, int op)
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
            z[i] = c_cpow(x[i], y);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = cbool(x[i] == y);
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = cbool(x[i] != y);
        }
        break;
    default:
        break;
    }
}

void x_op_vec_y(Cmplx * z,
                Cmplx x,
                const(Cmplx) * y,
                int n, int op)
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
            z[i] = c_cpow(x, y[i]);
        }
        break;
    case '=':
        for (i=0; i<n; i++) {
            z[i] = cbool(x == y[i]);
        }
        break;
    case '!':
        for (i=0; i<n; i++) {
            z[i] = cbool(x != y[i]);
        }
        break;
    default:
        break;
    }
}

Cmplx x_op_y(Cmplx x, Cmplx y, int op)
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
        return c_cpow(x, y);
    case '=':
        return cbool(x == y);
    case '!':
        return cbool(x != y);
    default:
        return Cmplx(0, 0);
    }
}

GretlMatrixRaw * cmatrix_dot_op(const(GretlMatrixRaw) * A,
                             const(GretlMatrixRaw) * B,
                             int op, int * err)
{
    GretlMatrixRaw * C = null;
    Cmplx x, y;
    int nr, nc;
    int conftype;
    int i, j, off;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
        *err = GretlError.data;
        return null;
    }

    conftype = dot_operator_conf(A, B, &nr, &nc);

    if (conftype == ConfType.none) {
        fputs("gretl_cmatrix_dot_op: matrices not conformable\n", stderr);
        fprintf(stderr, " op = '%c', A is %d x %d, B is %d x %d\n".ptr,
                cast(char) op, A.rows/2, A.cols, B.rows/2, B.cols);
        *err = GretlError.nonconf;
        return null;
    }

    C = gretl_cmatrix_new(nr, nc);
    if (C == null) {
        *err = GretlError.alloc;
        return null;
    }

    switch (conftype) {
    case ConfType.elements:
        vec_x_op_vec_y(C.z, A.z, B.z, nr*nc, op);
        break;
    case ConfType.a_colvec:
        for (i=0; i<nr; i++) {
            x = A.z[i];
            for (j=0; j<nc; j++) {
                y = gretl_cmatrix_get(B, i, j);
                y = x_op_y(x, y, op);
                gretl_cmatrix_set(C, i, j, y);
            }
        }
        break;
    case ConfType.b_colvec:
        for (i=0; i<nr; i++) {
            y = B.z[i];
            for (j=0; j<nc; j++) {
                x = gretl_cmatrix_get(A, i, j);
                x = x_op_y(x, y, op);
                gretl_cmatrix_set(C, i, j, x);
            }
        }
        break;
    case ConfType.a_rowvec:
        off = 0;
        for (j=0; j<nc; j++) {
            x = A.z[j];
            x_op_vec_y(C.z + off, x, B.z + off, nr, op);
            off += nr;
        }
        break;
    case ConfType.b_rowvec:
        off = 0;
        for (j=0; j<nc; j++) {
            y = B.z[j];
            vec_x_op_y(C.z + off, A.z + off, y, nr, op);
            off += nr;
        }
        break;
    case ConfType.a_scalar:
        x_op_vec_y(C.z, A.z[0], B.z, nr*nc, op);
        break;
    case ConfType.b_scalar:
        vec_x_op_y(C.z, A.z, B.z[0], nr*nc, op);
        break;
    case ConfType.ac_br:
        for (i=0; i<nr; i++) {
            x = A.z[i];
            for (j=0; j<nc; j++) {
                y = B.z[j];
                y = x_op_y(x, y, op);
                gretl_cmatrix_set(C, i, j, y);
            }
        }
        break;
    case ConfType.ar_bc:
        for (j=0; j<nc; j++) {
            x = A.z[j];
            for (i=0; i<nr; i++) {
                y = B.z[i];
                y = x_op_y(x, y, op);
                gretl_cmatrix_set(C, i, j, y);
            }
        }
        break;
    default: /* hush a warning */
        break;
    }

    if (errno) {
        if (*err) {
            gretl_matrix_free(C);
            C = null;
        }
    }

    return C;
}

GretlMatrixRaw * gretl_cmatrix_dot_op(const(GretlMatrixRaw) * A,
                                   const(GretlMatrixRaw) * B,
                                   int op, int * err)
{
    GretlMatrixRaw * T = null;
    GretlMatrixRaw * C = null;

    if (A.is_complex && B.is_complex) {
        C = cmatrix_dot_op(A, B, op, err);
    } else if (A.is_complex) {
        /* case of real B */
        T = complex_from_real(B, err);
        if (T != null) {
            C = cmatrix_dot_op(A, T, op, err);
        }
    } else if (B.is_complex) {
        /* case of real A */
        T = complex_from_real(A, err);
        if (T != null) {
            C = cmatrix_dot_op(T, B, op, err);
        }
    } else {
        *err = GretlError.types;
    }

    gretl_matrix_free(T);

    return C;
}

GretlMatrixRaw * gretl_cmatrix_divide(const(GretlMatrixRaw) * A,
                                   const(GretlMatrixRaw) * B,
                                   GretlMatrixMod mod,
                                   int * err)
{
    GretlMatrixRaw * T = null;
    GretlMatrixRaw * C = null;

    if (A.is_complex && B.is_complex) {
        C = gretl_matrix_divide(A, B, mod, err);
    } else if (A.is_complex) {
        /* case of real B */
        T = complex_from_real(B, err);
        if (T != null) {
            C = gretl_matrix_divide(A, T, mod, err);
        }
    } else if (B.is_complex) {
        /* case of real A */
        T = complex_from_real(A, err);
        if (T != null) {
            C = gretl_matrix_divide(T, B, mod, err);
        }
    } else {
        *err = GretlError.types;
    }

    gretl_matrix_free(T);

    return C;
}

GretlMatrixRaw * cmatrix_get_element(const(GretlMatrixRaw) * M, int i, int * err)
{
    GretlMatrixRaw * ret = null;

    if (M == null) {
        *err = GretlError.data;
    } else if (i < 0 || i >= M.rows * M.cols) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
        *err = GretlError.invarg;
    } else {
        ret = gretl_cmatrix_from_scalar(M.z[i], err);
    }

    return ret;
}

/* Set either the real or the imaginary part of @targ,
   using @src if non-NULL or @x otherwise.
*/

int gretl_cmatrix_set_part(GretlMatrixRaw * targ,
                           const(GretlMatrixRaw) * src,
                           double x, int im)
{
    Cmplx z;
    int i, j;

    if (!cmatrix_validate(targ, 0)) {
        return GretlError.types;
    } else if (src != null) {
        if (gretl_is_null_matrix(src) ||
            src.rows != targ.rows ||
            src.cols != targ.cols ||
            src.is_complex) {
            return GretlError.invarg;
        }
    }

    for (j=0; j<targ.cols; j++) {
        for (i=0; i<targ.rows; i++) {
            if (src != null) {
                x = gretl_matrix_get(src, i, j);
            }
            z = gretl_cmatrix_get(targ, i, j);
            if (im) {
                z = Cmplx(z.r, x);
            } else {
                z = Cmplx(x, z.i);
            }
            gretl_cmatrix_set(targ, i, j, z);
        }
    }

    return 0;
}

int imag_part_zero(const(GretlMatrixRaw) * A)
{
    int i, n = A.rows * A.cols;
    double tol = 1.0e-15; /* ?? */

    for (i=0; i<n; i++) {
        if (fabs(A.z[i].i) > tol) {
            fprintf(stderr, "imag_part_zero? no, got %g\n".ptr,
                    A.z[i].i);
            return 0;
        }
    }
    return 1;
}

/* Matrix logarithm, for a diagonalizable matrix */

GretlMatrixRaw * gretl_matrix_log(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * Vi = null;
    GretlMatrixRaw * w = null;
    GretlMatrixRaw * R = null;
    int i;

    if (gretl_is_null_matrix(A) || A.rows != A.cols) {
        /* we require a square matrix, real or complex */
        *err = GretlError.invarg;
        return null;
    }

    if (A.is_complex) {
        C = gretl_matrix_copy(A);
        if (C == null) {
            *err = GretlError.alloc;
        }
    } else {
        C = gretl_cmatrix_build(A, null, 0, 0, err);
    }

    if (!*err) {
        V = gretl_cmatrix_new(A.rows, A.rows);
        if (V == null) {
            *err = GretlError.alloc;
        } else {
            w = gretl_zgeev(C, V, null, err);
        }
    }

    if (!*err) {
        /* The following step will fail if
           @A is not diagonalizable
        */
        Vi = gretl_cmatrix_inverse(V, err);
    }

    if (!*err) {
        for (i=0; i<A.rows; i++) {
            w.z[i] = c_clog(w.z[i]);
        }
        R = cmatrix_dot_op(w, Vi, '*', err);
        if (!*err) {
            gretl_matrix_free(Vi);
            Vi = gretl_cmatrix_multiply(V, R, err);
        }
    }

    gretl_matrix_free(C);

    if (!*err) {
        if (imag_part_zero(Vi)) {
            /* should we do this? */
            C = gretl_cmatrix_extract(Vi, 0, err);
        } else {
            C = Vi;
            Vi = null;
        }
    } else {
        C = null;
    }

    gretl_matrix_free(V);
    gretl_matrix_free(Vi);
    gretl_matrix_free(w);
    gretl_matrix_free(R);

    return C;
}

/* Matrix exponential for a diagonalizable complex matrix */

GretlMatrixRaw * gretl_cmatrix_exp(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * Vi = null;
    GretlMatrixRaw * w = null;
    GretlMatrixRaw * R = null;
    int i;

    if (!cmatrix_validate(A, 1)) {
        *err = GretlError.invarg;
        return null;
    }

    C = gretl_matrix_copy(A);
    if (C == null) {
        *err = GretlError.alloc;
    }

    if (!*err) {
        V = gretl_cmatrix_new(A.rows, A.rows);
        if (V == null) {
            *err = GretlError.alloc;
        } else {
            w = gretl_zgeev(C, V, null, err);
        }
    }

    if (!*err) {
        /* The following step will fail if
           @A is not diagonalizable
        */
        Vi = gretl_cmatrix_inverse(V, err);
    }

    gretl_matrix_free(C);
    C = null;

    if (!*err) {
        for (i=0; i<A.rows; i++) {
            w.z[i] = c_cexp(w.z[i]);
        }
        R = cmatrix_dot_op(w, Vi, '*', err);
        if (!*err) {
            C = gretl_cmatrix_multiply(V, R, err);
        }
    }

    gretl_matrix_free(V);
    gretl_matrix_free(Vi);
    gretl_matrix_free(w);
    gretl_matrix_free(R);

    return C;
}

int matrix_is_hermitian(const(GretlMatrixRaw) * z)
{
    Cmplx zij, zji;
    int i, j, n = z.rows;

    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            zij = gretl_cmatrix_get(z, i, j);
            zji = gretl_cmatrix_get(z, j, i);
            if (zji != conj(zij)) {
                return 0;
            }
        }
    }
    return 1;
}

/* Cholesky decomposition of Hermitian matrix using LAPACK
   function zpotrf()
*/

GretlMatrixRaw * gretl_cmatrix_cholesky(const(GretlMatrixRaw) * A, int * err)
{
    GretlMatrixRaw * C = null;
    char uplo = 'L';
    integer n, info;

    if (!cmatrix_validate(A, 1)) {
        *err = GretlError.invarg;
    } else if (!matrix_is_hermitian(A)) {
        gretl_errmsg_set(_("Matrix is not Hermitian").ptr);
        *err = GretlError.invarg;
    } else {
        C = gretl_matrix_copy(A);
        if (C == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        return null;
    }

    n = A.rows;
    zpotrf_(&uplo, &n, cast(Cmplx *) C.val, &n, &info);
    if (info != 0) {
        fprintf(stderr, "gretl_cmatrix_cholesky: zpotrf failed with info = %d (n = %d)\n".ptr,
                cast(int) info, cast(int) n);
        *err = (info > 0)? GretlError.notpd : GretlError.data;
    }

    if (*err) {
        gretl_matrix_free(C);
        C = null;
    } else {
        gretl_matrix_zero_upper(C);
    }

    return C;
}

GretlMatrixRaw * gretl_cmatrix_QR_pivot_decomp(const(GretlMatrixRaw) * A,
                                            GretlMatrixRaw * R,
                                            GretlMatrixRaw * P,
                                            int * err)
{
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * Rtmp = null;
    GretlMatrixRaw * Ptmp = null;
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    integer * jpvt = null;
    Cmplx * tau = null;
    Cmplx * work = null;
    double * rwork = null;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    lda = m = A.rows;
    n = A.cols;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns").ptr);
        *err = GretlError.nonconf;
        return null;
    }

    Q = gretl_matrix_copy(A);
    if (Q == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (R != null) {
        *err = ensure_aux_cmatrix(R, &Rtmp, null, n, n);
        if (*err) {
            goto bailout;
        } else if (Rtmp != null) {
            gretl_matrix_replace_content(R, Rtmp);
        }
    }
    if (P != null) {
        *err = ensure_aux_matrix(P, &Ptmp, 1, n);
        if (*err) {
            goto bailout;
        } else if (Ptmp != null) {
            gretl_matrix_replace_content(P, Ptmp);
        }
    }

    /* dim of tau is min (m, n) */
    tau = cast(Cmplx *) malloc(n * (*tau).sizeof);
    work = cast(Cmplx *) malloc((*work).sizeof);
    rwork = cast(double *) malloc(2 * n * (*rwork).sizeof);
    if (tau == null || work == null || rwork == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* pivot array */
    jpvt = cast(integer *) malloc(n * (*jpvt).sizeof);
    if (jpvt == null) {
        *err = GretlError.alloc;
        goto bailout;
    } else {
        for (i=0; i<n; i++) {
            jpvt[i] = 0;
        }
    }

    /* workspace size query */
    zgeqp3_(&m, &n, cast(Cmplx *) Q.val, &lda, jpvt, tau, work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqp3: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
    } else {
        /* optimally sized work array */
        lwork = cast(integer) work[0].r;
        work = cast(Cmplx *) realloc(work, cast(size_t) lwork * (*work).sizeof);
        if (work == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        goto bailout;
    }

    /* run actual QR factorization */
    zgeqp3_(&m, &n, cast(Cmplx *) Q.val, &lda, jpvt, tau, work, &lwork, rwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqp3: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
        goto bailout;
    }

    if (R != null) {
        /* copy the upper triangular R out of Q */
        Cmplx z;

        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    z = gretl_cmatrix_get(Q, i, j);
                    gretl_cmatrix_set(R, i, j, z);
                } else {
                    gretl_cmatrix_set(R, i, j, Cmplx(0, 0));
                }
            }
        }
    }

    /* turn Q into "the real" Q, with the help of tau */
    zungqr_(&m, &n, &n, cast(Cmplx *) Q.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zungqr: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
    }

    if (P != null) {
        for (i=0; i<n; i++) {
            P.val[i] = jpvt[i];
        }
    }

 bailout:

    free(tau);
    free(work);
    free(rwork);
    free(jpvt);
    gretl_matrix_free(Rtmp);
    gretl_matrix_free(Ptmp);

    if (*err) {
        gretl_matrix_free(Q);
        Q = null;
    }

    return Q;
}

/* QR decomposition of complex matrix using LAPACK functions
   zgeqrf() and zungqr().
*/

GretlMatrixRaw * gretl_cmatrix_QR_decomp(const(GretlMatrixRaw) * A,
                                      GretlMatrixRaw * R,
                                      int * err)
{
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * Rtmp = null;
    integer m, n, lda;
    integer info = 0;
    integer lwork = -1;
    Cmplx * tau = null;
    Cmplx * work = null;
    int i, j;

    if (!cmatrix_validate(A, 0)) {
        *err = GretlError.invarg;
        return null;
    }

    lda = m = A.rows;
    n = A.cols;

    if (n > m) {
        gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns").ptr);
        *err = GretlError.nonconf;
        return null;
    }

    Q = gretl_matrix_copy(A);
    if (Q == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (R != null) {
        *err = ensure_aux_cmatrix(R, &Rtmp, null, n, n);
        if (*err) {
            goto bailout;
        } else if (Rtmp != null) {
            gretl_matrix_replace_content(R, Rtmp);
        }
    }

    /* dim of tau is min (m, n) */
    tau = cast(Cmplx *) malloc(n * (*tau).sizeof);
    work = cast(Cmplx *) malloc((*work).sizeof);
    if (tau == null || work == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    /* workspace size query */
    zgeqrf_(&m, &n, cast(Cmplx *) Q.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqrf: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
    } else {
        /* optimally sized work array */
        lwork = cast(integer) work[0].r;
        work = cast(Cmplx *) realloc(work, cast(size_t) lwork * (*work).sizeof);
        if (work == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        goto bailout;
    }

    /* run actual QR factorization */
    zgeqrf_(&m, &n, cast(Cmplx *) Q.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zgeqrf: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
        goto bailout;
    }

    if (R != null) {
        /* copy the upper triangular R out of Q */
        Cmplx z;

        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (i <= j) {
                    z = gretl_cmatrix_get(Q, i, j);
                    gretl_cmatrix_set(R, i, j, z);
                } else {
                    gretl_cmatrix_set(R, i, j, Cmplx(0, 0));
                }
            }
        }
    }

    /* turn Q into "the real" Q, with the help of tau */
    zungqr_(&m, &n, &n, cast(Cmplx *) Q.val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        fprintf(stderr, "zungqr: info = %d\n".ptr, cast(int) info);
        *err = GretlError.data;
    }

 bailout:

    free(tau);
    free(work);
    gretl_matrix_free(Rtmp);

    if (*err) {
        gretl_matrix_free(Q);
        Q = null;
    }

    return Q;
}

/* Copy data from real matrix @src to complex matrix *targ,
   starting at offset (@r0, @c0) into @targ. @src is treated
   as a complex matrix with zero imaginary part.
*/

void real_to_complex_fill(GretlMatrixRaw * targ,
                          const(GretlMatrixRaw) * src,
                          int r0, int c0)
{
    Cmplx z;
    int i, j, r, c;

    for (j=0, c=c0; j<src.cols && c<targ.cols; j++, c++) {
        for (i=0, r=r0; i<src.rows && r<targ.rows; i++, r++) {
            z = Cmplx(gretl_matrix_get(src, i, j), 0);
            gretl_cmatrix_set(targ, r, c, z);
        }
    }
}

/* Tests whether a matrix has the is_complex flag set, and also
   whether it's "really complex" (has a non-zero imaginary part).
*/

int matrix_is_complex(const(GretlMatrixRaw) * M)
{
    if (M == null || !M.is_complex) {
        return 0;
    } else {
        int i, n = M.rows * M.cols;
        int ret = 1;

        for (i=0; i<n; i++) {
            if (M.z[i].i != 0) {
                ret = 2;
                break;
            }
        }
        return ret;
    }
}
