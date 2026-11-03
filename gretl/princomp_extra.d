/* Port of princomp_extra.c */

module gretl.princomp_extra;


/*
  Note 2024-11-18: This file is (hopefully) a temporary hack, an
  intermediate step on the way to a proper renovation of gretl's
  princomp() function.

  NB: in the C source this file is #included into gretl_matrix.c; here it
  is its own module, importing from gretl_matrix (which imports back the
  two public functions below).
*/

import core.stdc.math : fabs, sqrt;
import core.stdc.stdio : fprintf, stderr;
import core.stdc.stdlib : free, malloc;

import gretl.glib;
import gretl.gretl_errorsh;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_copy,
    gretl_matrix_copy_values, gretl_matrix_free, gretl_matrix_multi_ols,
    gretl_matrix_multiply, gretl_zero_matrix_new, real_gretl_matrix_pca;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.libset : libset_get_user_tolerance;
import gretl.libseth : SetKey;
import gretl.missingh;
import gretl.gretl_errors : gretl_errmsg_sprintf;

/* ---- declarations of not-yet-ported items used below ---- */

extern (C) {
}

/* ---- end of declarations ---- */

int matrix_is_incomplete(const(GretlMatrixRaw) * m)
{
    int i, n = m.rows * m.cols;

    for (i=0; i<n; i++) {
        if (na(m.val[i])) {
            return 1;
        }
    }

    return 0;
}

private gint8 * makemask(const(GretlMatrixRaw) * X, int * nfull, int * ndrop)
{
    gint8 * ret;
    double xij;
    int i, j;
    int nna;

    ret = cast(gint8 *) malloc(X.rows * (*ret).sizeof);

    for (i=0; i<X.rows; i++) {
        nna = 0;
        for (j=0; j<X.cols; j++) {
            xij = gretl_matrix_get(X, i, j);
            if (na(xij)) {
                nna++;
            }
        }
        if (nna == 0) {
            ret[i] = 1;
            *nfull += 1;
        } else if (nna == X.cols) {
            ret[i] = -1;
            *ndrop += 1;
        } else {
            ret[i] = 0;
        }
    }

    return ret;
}

private GretlMatrixRaw * select_rows_from(const(GretlMatrixRaw) * M,
                                       int nr,
                                       const(gint8) * mask)
{
    GretlMatrixRaw * X = gretl_matrix_alloc(nr, M.cols);
    int i, j, k = 0;
    double mij;

    for (i=0; i<M.rows; i++) {
        if (mask[i] > 0) {
            for (j=0; j<M.cols; j++) {
                mij = gretl_matrix_get(M, i, j);
                gretl_matrix_set(X, k, j, mij);
            }
            k++;
        }
    }

    return X;
}

private void F_insert_P(GretlMatrixRaw * F,
                        const(GretlMatrixRaw) * P,
                        const(gint8) * mask)
{
    double pkj;
    int i, j, k = 0;

    for (i=0; i<F.rows; i++) {
        if (mask[i] > 0) {
            for (j=0; j<F.cols; j++) {
                pkj = gretl_matrix_get(P, k, j);
                gretl_matrix_set(F, i, j, pkj);
            }
            k++;
        }
    }
}

/* Copy from @a to @targ unless an @a element is missing,
   in which case copy that element from @b instead.
*/

private void alt_transcribe_values(GretlMatrixRaw * targ,
                                   const(GretlMatrixRaw) * a,
                                   const(GretlMatrixRaw) * b)
{
    int i, n = targ.rows * targ.cols;

    for (i=0; i<n; i++) {
        targ.val[i] = na(a.val[i]) ? b.val[i] : a.val[i];
    }
}

private double max_abs_value(const(GretlMatrixRaw) * E)
{
    int i, n = E.rows * E.cols;
    double ai, ret = 0;

    for (i=0; i<n; i++) {
        ai = fabs(E.val[i]);
        if (ai > ret) {
            ret = ai;
        }
    }

    return ret;
}

private int demean_with_NAs(GretlMatrixRaw * X, int stdize)
{
    double * Xj = X.val;
    double cmean, csd;
    int i, j, Tj;
    int err = 0;

    for (j=0; j<X.cols; j++) {
        cmean = csd = 0.0;
        Tj = 0;
        for (i=0; i<X.rows; i++) {
            if (!na(Xj[i])) {
                cmean += Xj[i];
                Tj++;
            }
        }
        if (Tj == 0) {
            gretl_errmsg_sprintf("Column %d is entirely missing", j+1);
            err = GretlError.missdata;
            break;
        }
        cmean /= Tj;
        for (i=0; i<X.rows; i++) {
            if (!na(Xj[i])) {
                Xj[i] -= cmean;
                if (stdize) {
                    csd += Xj[i] * Xj[i];
                }
            }
        }
        if (stdize) {
            csd = sqrt(csd / Tj);
            for (i=0; i<X.rows; i++) {
                if (!na(Xj[i])) {
                    Xj[i] /= csd;
                }
            }
        }
        Xj += X.rows;
    }

    return err;
}

GretlMatrixRaw * princomp_with_NAs(const(GretlMatrixRaw) * M,
                                int p, Gretlopt opt,
                                int * err)
{
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * Tmp = null;
    GretlMatrixRaw * F = null;
    GretlMatrixRaw * L = null;
    GretlMatrixRaw * FL = null;
    GretlMatrixRaw * F0 = null;
    GretlMatrixRaw * E = null;
    Gretlopt pca_opt = Gretlopt.OPT_NONE;
    double tol = libset_get_user_tolerance(SetKey.bhhh_toler);
    double crit = NADBL;
    gint8 * mask = null;
    int maxiter = 2048;
    int stdize = 1;
    int conv = 0;
    int nfull = 0;
    int ndrop = 0;
    int iter = 0;

    if (opt & Gretlopt.OPT_C) {
        stdize = 0;
        pca_opt |= Gretlopt.OPT_V; /* FIXME? */
    }

    mask = makemask(M, &nfull, &ndrop);

    if (p <= 0 || p > M.cols) {
        gretl_errmsg_sprintf("%d: invalid number of components", p);
        *err = GretlError.invarg;
        return null;
    }

    if (ndrop > 0) {
        *err = GretlError.missdata;
        gretl_errmsg_sprintf("%d rows were entirely missing", ndrop);
    } else if (nfull < p) {
        *err = GretlError.toofew;
        gretl_errmsg_sprintf("%d complete rows are needed but only %d were found",
                             p, nfull);
    } else {
        X = gretl_matrix_copy(M);
        if (X == null) {
            *err = GretlError.alloc;
        }
    }

    if (*err) {
        free(mask);
        return null;
    }

    /* demean, and possibly standardize */
    *err = demean_with_NAs(X, stdize);

    if (!*err) {
        /* initialize loadings via PCA on complete observations */
        GretlMatrixRaw * P;

        Tmp = select_rows_from(X, nfull, mask);
        F = gretl_zero_matrix_new(X.rows, p);
        P = real_gretl_matrix_pca(Tmp, p, pca_opt, err);
        F_insert_P(F, P, mask);
        L = gretl_zero_matrix_new(P.cols, Tmp.cols);
        *err = gretl_matrix_multi_ols(Tmp, P, L, null, null);
        gretl_matrix_free(P);
        gretl_matrix_free(Tmp);
        Tmp = null;
    }

    if (!*err) {
        /* additional matrix allocation */
        Tmp = gretl_matrix_alloc(X.rows, X.cols);
        F0 = gretl_matrix_alloc(F.rows, F.cols);
        FL = gretl_matrix_alloc(X.rows, X.cols);
        E = gretl_matrix_alloc(X.rows, F.cols);
    }

    /* EM iteration */
    while (!conv && !*err && iter < maxiter) {
        gretl_matrix_copy_values(F0, F);
        gretl_matrix_multiply(F, L, FL);
        alt_transcribe_values(Tmp, X, FL);
        gretl_matrix_free(F);
        F = real_gretl_matrix_pca(Tmp, p, pca_opt, err);
        if (!*err) {
            *err = gretl_matrix_multi_ols(F, F0, null, E, null);
        }
        if (*err) {
            break;
        }
        crit = max_abs_value(E);
        if (crit < tol) {
            conv = 1;
        }
        *err = gretl_matrix_multi_ols(Tmp, F, L, null, null);
        iter++;
    }

    if (!*err && !conv) {
        fprintf(stderr, "PC EM algo didn't converge: crit = %g\n", crit);
        *err = GretlError.noconv;
    }
    if (*err && F != null) {
        gretl_matrix_free(F);
        F = null;
    }

    gretl_matrix_free(X);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(FL);
    gretl_matrix_free(F0);
    gretl_matrix_free(E);
    free(mask);

    return F;
}
