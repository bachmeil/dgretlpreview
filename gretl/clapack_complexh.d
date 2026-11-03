/* Port of clapack_complex.h */

module gretl.clapack_complexh;


/* LAPACK subroutines: double-precision complex versions only.

   As in clapack_doubleh.d these are declarations of external Fortran
   routines, kept as extern (C) declarations. Complex arrays are passed as
   Cmplx * (the C "cmplx" type: two doubles, same layout as C complex double). */

import gretl.gretl_f2ch;
import gretl.libgretlh;

extern (C) {

void zheev_(const(char) * jobz, const(char) * uplo, integer * n, Cmplx * a,
            integer * lda, double * w, Cmplx * work, integer * lwork,
            double * rwork, integer * info);

void zgeev_(const(char) * jobvl, const(char) * jobvr, integer * n, Cmplx * a,
            integer * lda, double * w, Cmplx * vl, integer * ldvl, Cmplx * vr,
            integer * ldvr, Cmplx * work, integer * lwork, double * rwork,
            integer * info);

void zgelsy_(integer * m, integer * n, integer * nrhs, Cmplx * a,
             integer * lda, Cmplx * b, integer * ldb, integer * jpvt,
             double * rcond, integer * rank, Cmplx * work, integer * lwork,
             double * rwork, integer * info);

void zgetrs_(char * trans, integer * n, integer * nrhs, Cmplx * a,
             integer * lda, integer * ipiv, Cmplx * b, integer * ldb,
             integer * info);

void zgetrf_(integer * m, integer * n, Cmplx * a, integer * lda,
             integer * ipiv, integer * info);

void zgetri_(integer * n, Cmplx * a, integer * lda, integer * ipiv,
             Cmplx * work, integer * lwork, integer * info);

void zgemm_(const(char) * transa, const(char) * transb, integer * m,
            integer * n, integer * k, Cmplx * alpha, Cmplx * a, integer * lda,
            Cmplx * b, integer * ldb, Cmplx * beta, Cmplx * c, integer * ldc);

void zsyrk_(const(char) * uplo, const(char) * trans, integer * n, integer * k,
            Cmplx * alpha, Cmplx * a, integer * lda, Cmplx * beta, Cmplx * c,
            integer * ldc);

void zgesvd_(const(char) * jobu, const(char) * jobvt, integer * m,
             integer * n, Cmplx * a, integer * lda, double * s, Cmplx * u,
             integer * ldu, Cmplx * vt, integer * ldvt, Cmplx * work,
             integer * lwork, double * rwork, integer * info);

void zgees_(const(char) * jobvs, const(char) * srt, void * select,
            integer * n, Cmplx * a, integer * lda, integer * sdim, Cmplx * w,
            Cmplx * vs, integer * ldvs, Cmplx * work, integer * lwork,
            double * rwork, integer * bwork, integer * info);

void zpotrf_(const(char) * uplo, integer * n, Cmplx * a, integer * lda,
             integer * info);

void zgeqrf_(integer * m, integer * n, Cmplx * a, integer * lda, Cmplx * tau,
             Cmplx * work, integer * lwork, integer * info);

void zungqr_(integer * m, integer * n, integer * k, Cmplx * a, integer * lda,
             Cmplx * tau, Cmplx * work, integer * lwork, integer * info);

void zgeqp3_(integer * m, integer * n, Cmplx * a, integer * lda,
             integer * jpvt, Cmplx * tau, Cmplx * work, integer * lwork,
             double * rwork, integer * info);

}
