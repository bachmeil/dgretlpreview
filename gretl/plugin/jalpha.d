/* Port of plugin/jalpha.c: computations concerned with testing restrictions on the alpha (adjustments) matrix in a VECM.
   Handles only homogeneous, common restrictions. Notation based on S. Johansen, "Likelihood-Based Inference in Cointegrated
   Vector Autoregressive Models" (Oxford, 1995). See section 8.2.1, "The same restriction on all \alpha", p. 124ff.

   In the C this file is part of the same plugin as johansen.c and jrestrict.c; here the functions are ordinary D functions of
   this module (vecm_alpha_test is called from gretl.plugin.johansen, and the other two modules import each other: the C
   prototypes of jprivate.h become imports). The debugging output under ADEBUG is left out. */

module gretl.plugin.jalpha;

import core.stdc.math : fabs;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn;
import gretl.gretl_restrict : rset_get_Ra_matrix;
import gretl.gretl_restricth : GretlRestriction;
import gretl.johansenh : JohansenInfo;
import gretl.libgretlh;
import gretl.libset : libset_get_int;
import gretl.libseth : SetKey, VECMnorm;
import gretl.plugin.johansen : VBeta, johansen_LR_calc, johansen_ll_calc, print_beta_alpha_Pi;
import gretl.varh : GretlVar;

/* \bar{A} = A(A'A)^{-1}, where A is the orthogonal complement of the
   restriction R imposed on \alpha in the form R\alpha = 0
*/

private GretlMatrixRaw * make_A_bar(const(GretlMatrixRaw) * R,
                                    GretlMatrixRaw ** pA,
                                    int * err)
{
    GretlMatrixRaw * A = null;
    GretlMatrixRaw * Abar = null;
    GretlMatrixRaw * Tmp = null;

    A = gretl_matrix_right_nullspace(R, err);
    if (*err) {
        return null;
    }

    Tmp = gretl_matrix_alloc(A.cols, A.cols);
    if (Tmp == null) {
        *err = GretlError.alloc;
        gretl_matrix_free(A);
        return null;
    }

    gretl_matrix_multiply_mod(A, GretlMatrixMod.transpose,
                              A, GretlMatrixMod.none,
                              Tmp, GretlMatrixMod.none);

    *err = gretl_invert_symmetric_matrix(Tmp);

    if (!*err) {
        Abar = gretl_matrix_multiply_new(A, Tmp, err);
    }

    if (!*err && pA != null) {
        *pA = A;
    } else {
        gretl_matrix_free(A);
    }

    gretl_matrix_free(Tmp);

    return Abar;
}

/* "AS00" is the counterpart to S_{00}^{-1} in the unconstrained
   case.  It is given by

   \bar{A} (\bar{A}'S_{00.A_{\perp}}\bar{A})^{-1} \bar{A}'

*/

private GretlMatrixRaw * make_AS00(const(GretlMatrixRaw) * R,
                                   const(GretlMatrixRaw) * S00a,
                                   GretlMatrixRaw ** pA,
                                   int * err)
{
    GretlMatrixRaw * AS00 = null;
    GretlMatrixRaw * Abar = null;
    GretlMatrixRaw * Tmp = null;

    Abar = make_A_bar(R, pA, err);
    if (*err) {
        return null;
    }

    Tmp = gretl_matrix_alloc(Abar.cols, Abar.cols);
    if (Tmp == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    *err = gretl_matrix_qform(Abar, GretlMatrixMod.transpose,
                              S00a, Tmp, GretlMatrixMod.none);

    if (!*err) {
        *err = gretl_invert_symmetric_matrix(Tmp);
    }

    if (!*err) {
        AS00 = gretl_matrix_alloc(Abar.rows, Abar.rows);
        if (AS00 == null) {
            *err = GretlError.alloc;
            goto bailout;
        }
    }

    if (!*err) {
        *err = gretl_matrix_qform(Abar, GretlMatrixMod.none, Tmp,
                                  AS00, GretlMatrixMod.none);
    }

 bailout:

    gretl_matrix_free(Tmp);
    gretl_matrix_free(Abar);

    return AS00;
}

/* Solve the eigenvalue problem

   |\tilde{\lambda}S_{11.A_{\perp}} - S_{10.A_{\perp}}
      \bar{A}(\bar{A}' S_{00.A_{\perp}} \bar{A})^{-1}
        \bar{A}'S_{01.A_{\perp}}| = 0

   (Johansen, p. 126)
*/

int alt_get_eigenvalues(GretlMatrixRaw * AS00,
                        const(GretlMatrixRaw) * S01a,
                        const(GretlMatrixRaw) * S11a,
                        GretlMatrixRaw * M,
                        GretlMatrixRaw ** evals,
                        GretlMatrixRaw * Tmp,
                        int rank)
{
    int err = 0;

    gretl_matrix_qform(S01a, GretlMatrixMod.transpose,
                       AS00, Tmp, GretlMatrixMod.none);

    *evals = gretl_gensymm_eigenvals(Tmp, S11a, M, &err);

    if (!err) {
        err = gretl_symmetric_eigen_sort(*evals, M, rank);
    }

    return err;
}

private int alpha_compute_alpha(JohansenInfo * jv,
                                const(GretlMatrixRaw) * S11a,
                                const(GretlMatrixRaw) * S01a)
{
    const(GretlMatrixRaw) * B = jv.Beta;
    GretlMatrixRaw * BSB = null;
    GretlMatrixRaw * Tmp = null;
    int err = 0;

    BSB = gretl_matrix_alloc(B.cols, B.cols);
    Tmp = gretl_matrix_alloc(B.rows, B.cols);

    if (BSB == null || Tmp == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = gretl_matrix_qform(B, GretlMatrixMod.transpose, S11a,
                                 BSB, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(BSB);
    }

    if (!err) {
        gretl_matrix_multiply(B, BSB, Tmp);
        gretl_matrix_multiply(S01a, Tmp, jv.Alpha);
    }

    gretl_matrix_free(BSB);
    gretl_matrix_free(Tmp);

    return err;
}

private int alpha_calc_full(GretlVar * jvar,
                            const(GretlMatrixRaw) * M,
                            const(GretlMatrixRaw) * S11a,
                            const(GretlMatrixRaw) * S01a)
{
    JohansenInfo * jv = jvar.jinfo;
    GretlMatrixRaw * B = null;
    int vnorm = libset_get_int(SetKey.vecm_norm);

    if (jv.Beta == null) {
        jv.Beta = gretl_matrix_copy(M);
    } else {
        gretl_matrix_copy_values(jv.Beta, M);
    }

    B = jv.Beta;
    if (B == null) {
        return GretlError.alloc;
    }

    /* if we're going to normalize beta (?) then do it
       before computing alpha */

    if (vnorm == VECMnorm.diag || vnorm == VECMnorm.first) {
        double x, den;
        int i, j, row;

        for (j=0; j<B.cols; j++) {
            row = (vnorm == VECMnorm.diag)? j : 0;
            den = gretl_matrix_get(B, row, j);
            if (den != 0.0) {
                for (i=0; i<B.rows; i++) {
                    x = gretl_matrix_get(B, i, j);
                    gretl_matrix_set(B, i, j, x / den);
                }
            }
        }
    }

    return alpha_compute_alpha(jv, S11a, S01a);
}

private void set_true_zeros(GretlMatrixRaw * m)
{
    int i, n = m.rows * m.cols;

    for (i=0; i<n; i++) {
        if (fabs(m.val[i]) < 3.0e-19) {
            m.val[i] = 0;
        }
    }
}

/* See Johansen's 1995 book, section 8.2.1, "The same restriction on
   all \alpha"
*/

int vecm_alpha_test(GretlVar * jvar,
                    GretlRestriction * rset,
                    const(DatasetRaw) * dset,
                    Gretlopt opt,
                    Prn * prn)
{
    const(GretlMatrixRaw) * Ap = rset_get_Ra_matrix(rset);
    const(GretlMatrixRaw) * S00 = jvar.jinfo.S00;
    GretlMatrixRaw * S01 = jvar.jinfo.S01;
    GretlMatrixRaw * S11 = jvar.jinfo.S11;
    GretlMatrixRaw * ASA = null;
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * Tmp = null;
    GretlMatrixRaw * S00a = null;
    GretlMatrixRaw * S11a = null;
    GretlMatrixRaw * S01a = null;
    int r = jvar.jinfo.rank;
    int p = jvar.neqns;
    int m = S11.rows;
    int err = 0;

    clear_gretl_matrix_err();

    ASA = gretl_matrix_alloc(Ap.rows, Ap.rows);
    C = gretl_matrix_alloc(p, p);
    Tmp = gretl_matrix_alloc(m, m);
    S00a = gretl_zero_matrix_new(p, p);
    S11a = gretl_zero_matrix_new(m, m);
    S01a = gretl_zero_matrix_new(p, m);

    err = get_gretl_matrix_err();
    if (err) {
        goto bailout;
    }

    /* form A\perp' S_{00} A_\perp in "ASA" */
    gretl_matrix_qform(Ap, GretlMatrixMod.none, S00,
                       ASA, GretlMatrixMod.none);

    /* "ASA" <- (A\perp' S_{00} A_\perp)^{-1} */
    gretl_invert_symmetric_matrix(ASA);

    /* C = A\perp' (A\perp' S_{00} A_\perp)^{-1} A_\perp */
    gretl_matrix_qform(Ap, GretlMatrixMod.transpose, ASA,
                       C, GretlMatrixMod.none);

    /* Johansen, p. 124: compute the A_{\perp}-transformed
       moment matrices */

    /* S_{00.A_{\perp}} */
    gretl_matrix_qform(S00, GretlMatrixMod.transpose, C,
                       S00a, GretlMatrixMod.none);
    gretl_matrix_subtract_reversed(S00, S00a);

    /* S_{11.A_{\perp}} */
    gretl_matrix_qform(S01, GretlMatrixMod.transpose, C,
                       S11a, GretlMatrixMod.none);
    gretl_matrix_subtract_reversed(S11, S11a);

    /* S_{01.A_{\perp}} */
    gretl_matrix_reuse(Tmp, p, p);
    gretl_matrix_multiply(S00, C, Tmp);
    gretl_matrix_multiply(Tmp, S01, S01a);
    gretl_matrix_subtract_reversed(S01, S01a);

    set_true_zeros(S00a);
    set_true_zeros(S01a);

    if (!err) {
        /* do the eigenvalues thing */
        GretlMatrixRaw * A = null;
        GretlMatrixRaw * AS00 = null;
        GretlMatrixRaw * M = null;
        GretlMatrixRaw * evals = null;

        AS00 = make_AS00(Ap, S00a, &A, &err);

        if (!err) {
            M = gretl_matrix_alloc(m, m);
            if (M == null) {
                err = GretlError.alloc;
            }
        }

        if (!err) {
            gretl_matrix_reuse(Tmp, m, m);
            err = alt_get_eigenvalues(AS00, S01a, S11a,
                                      M, &evals, Tmp, r);
        }

        if (!err) {
            if (opt & Gretlopt.OPT_F) {
                johansen_ll_calc(jvar, evals);
                jvar.jinfo.lrdf = r * (p - A.cols);
            } else if (!(opt & Gretlopt.OPT_B)) {
                err = johansen_LR_calc(jvar, evals, A, rset, VBeta.alpha, prn);
            }
        }

        if (!err && (opt & (Gretlopt.OPT_V | Gretlopt.OPT_F))) {
            err = alpha_calc_full(jvar, M, S11a, S01a);
        }

        if (!err && (opt & Gretlopt.OPT_V)) {
            print_beta_alpha_Pi(jvar, dset, prn);
        }

        gretl_matrix_free(evals);
        gretl_matrix_free(M);
        gretl_matrix_free(AS00);
        gretl_matrix_free(A);
    }

 bailout:

    gretl_matrix_free(ASA);
    gretl_matrix_free(C);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(S00a);
    gretl_matrix_free(S11a);
    gretl_matrix_free(S01a);

    return err;
}
