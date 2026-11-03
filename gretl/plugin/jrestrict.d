/* Port of plugin/jrestrict.c: general restrictions on beta (and alpha) in the context of a VECM: identification checks,
   the switching algorithm and (L)BFGS maximization of the restricted likelihood, variances from the information matrix,
   output. Based on an ox program by Jack, June 2007; rendered into C by Allin.

   In the C this file is part of the same plugin as johansen.c and jalpha.c; here general_vecm_analysis, the entry point used by
   gretl.plugin.johansen, is an ordinary D function of this module (the C prototype of jprivate.h becomes an import). The
   typedefs gradhelper_, switcher_ and Jwrap_ are the structs GradHelper, Switcher and Jwrap, the flag macros are the enum JFlags
   with small functions, the anonymous enum OMEGA_ONLY/OMEGA_PLUS is OmegaCode, the BFGS callbacks Jloglik and Jgradient are
   extern (C) as BFGS_CRIT_FUNC and BFGS_GRAD_FUNC require, a C char used as a signed flag (in check_for_scaling) is a byte, and
   the debugging output under JDEBUG is left out. */

module gretl.plugin.jrestrict;

import core.stdc.math : fabs, floor, log, sqrt;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import std.algorithm : max, min;

import gretl.gretl_bfgs;
import gretl.gretl_bfgsh : BFGS_CRIT_FUNC;
import gretl.nlsh : CriterionType;
import gretl.gretl_list : gretl_list_new, in_gretl_list;
import gretl.pvaluesh : DistCode;
import gretl.random;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_restrict : gretl_restriction_get_options, rset_VECM_acols, rset_add_results, rset_VECM_bcols, rset_get_R_matrix, rset_get_Ra_matrix, rset_get_q_matrix, rset_get_qa_matrix;
import gretl.gretl_restricth : GretlRestriction;
import gretl.johansenh : JohansenCode, JohansenInfo, auto_restr, jcode, jrank;
import gretl.libgretlh;
import gretl.libset : get_initvals, libset_get_int, n_initvals;
import gretl.libseth : SetKey;
import gretl.missingh : NADBL, na;
import gretl.printout : bufspace, maybe_trim_varname;
import gretl.pvalues : chisq_cdf_comp;
import gretl.printouth : SHORTLEN;
import gretl.varh;
import gretl.varprint : vecm_beta_varname;

enum JFlags {
    use_lbfgs   = 1 << 0,
    cross_alpha = 1 << 1
}

struct Jwrap {
    ubyte flags;    /* status bitflags */
    int T;          /* length of time series used */
    int p;          /* number of equations */
    int p1;         /* number of rows in beta (>= p) */
    int r;          /* cointegrating rank of VECM */
    int blen;       /* number of unrestricted coefficients in beta */
    int alen;       /* number of unrestricted coefficients in alpha */
    int df;         /* degrees of freedom for LR test */
    int jr;         /* rank of Jacobian */
    double llk;     /* constant in log-likelihood */
    double ll;      /* log-likelihood */

    /* moment matrices and copies */
    const(GretlMatrixRaw) * R0;
    const(GretlMatrixRaw) * R1;
    const(GretlMatrixRaw) * S00;
    const(GretlMatrixRaw) * S01;
    const(GretlMatrixRaw) * S11;

    GretlMatrixRaw * S00i;

    /* restrictions on beta */
    GretlMatrixRaw * H;
    GretlMatrixRaw * h;

    /* restrictions on beta, incl. normalizations */
    GretlMatrixRaw * H0;

    /* augmented beta restriction matrices, if needed */
    GretlMatrixRaw * bigR;
    GretlMatrixRaw * bigq;

    /* homogeneous restrictions on alpha */
    GretlMatrixRaw * G;

    /* coefficients and variances */
    GretlMatrixRaw * beta;
    GretlMatrixRaw * alpha;
    GretlMatrixRaw * Pi;
    GretlMatrixRaw * lsPi;
    GretlMatrixRaw * Omega;
    GretlMatrixRaw * iOmega;
    GretlMatrixRaw * Vb;
    GretlMatrixRaw * Va;
    GretlMatrixRaw * bse;
    GretlMatrixRaw * ase;

    /* information matrix diagonal blocks */
    GretlMatrixRaw * I00;
    GretlMatrixRaw * I11;

    /* free parameter vector (BFGS or jittering only) */
    GretlMatrixRaw * theta;

    /* temp storage for beta calculation */
    GretlMatrixRaw * phi;

    /* alpha calculation */
    GretlMatrixRaw * psi;

    /* temporary storage */
    GretlMatrixRaw * qf1;
    GretlMatrixRaw * qf2;
    GretlMatrixRaw * Tmppp;
    GretlMatrixRaw * Tmprp;

    /* for analytical derivatives in LBFGS */
    GradHelper * ghelper;

    /* for keeping track of column scaling */
    int * normrow;
    int * normcol;
    double * normval;
}

private bool do_scaling(const(Jwrap) * J) { return J.normrow != null; }
private bool using_bfgs(const(Jwrap) * J) { return (J.flags & JFlags.use_lbfgs) != 0; }
private bool cross_alpha(const(Jwrap) * J) { return (J.flags & JFlags.cross_alpha) != 0; }

private int make_moment_matrices(Jwrap * J, const(GretlVar) * jvar)
{
    GretlMatrixRaw * Tmp = null;
    int err = 0;

    J.R0 = jvar.jinfo.R0;
    J.R1 = jvar.jinfo.R1;

    J.S00 = jvar.jinfo.S00;
    J.S01 = jvar.jinfo.S01;
    J.S11 = jvar.jinfo.S11;

    J.S00i = gretl_matrix_copy(J.S00);

    if (J.S00i == null) {
        return GretlError.alloc;
    }

    J.p1 = J.S01.cols;

    Tmp = gretl_matrix_alloc(J.p1, J.p1);

    if (Tmp == null) {
        return GretlError.alloc;
    }

    if (!err) {
        gretl_matrix_copy_values(J.S00i, J.S00);
        err = gretl_invert_symmetric_matrix(J.S00i);
    }

    if (!err) {
        err = gretl_matrix_qform(J.S01, GretlMatrixMod.transpose,
                                 J.S00i, Tmp, GretlMatrixMod.none);
    }

    if (!err) {
        /* allocate beta, alpha, etc. while we're at it:
           but are these sizes always right? */
        J.beta = gretl_zero_matrix_new(J.p1, J.r);
        J.alpha = gretl_zero_matrix_new(J.p, J.r);
        J.Pi = gretl_matrix_alloc(J.p, J.p1);
        J.Omega = gretl_matrix_alloc(J.p, J.p);
        if (J.beta == null || J.alpha == null ||
            J.Pi == null || J.Omega == null) {
            err = GretlError.alloc;
        } else {
            J.blen = J.p1 * J.r;
            J.alen = J.p * J.r;
        }
    }

    gretl_matrix_free(Tmp);

    return err;
}

private int J_alloc_aux(Jwrap * J)
{
    clear_gretl_matrix_err();

    J.qf1 = gretl_matrix_alloc(J.r, J.r);
    J.qf2 = gretl_matrix_alloc(J.r, J.r);
    J.Tmppp = gretl_matrix_alloc(J.p, J.p);
    J.Tmprp = gretl_matrix_alloc(J.r, J.p);

    return get_gretl_matrix_err();
}

private int allocate_info_blocks_etc(Jwrap * J)
{
    int err = 0;

    if (J.G != null) {
        J.I00 = gretl_matrix_alloc(J.alen, J.alen);
        if (J.I00 == null) {
            err = GretlError.alloc;
        }
    }

    if (J.blen > 0 && !err) {
        J.I11 = gretl_matrix_alloc(J.blen, J.blen);
        if (J.I11 == null) {
            err = GretlError.alloc;
        }
    }

    if (J.iOmega == null && !err) {
        J.iOmega = gretl_matrix_alloc(J.p, J.p);
        if (J.iOmega == null) {
            err = GretlError.alloc;
        }
    }

    return err;
}

private void jwrap_destroy(Jwrap * J)
{
    gretl_matrix_free(J.S00i);

    gretl_matrix_free(J.G);

    if (J.H0 == J.H) {
        J.H0 = null;
    }

    gretl_matrix_free(J.H);
    gretl_matrix_free(J.h);
    gretl_matrix_free(J.H0);

    gretl_matrix_free(J.bigR);
    gretl_matrix_free(J.bigq);

    gretl_matrix_free(J.beta);
    gretl_matrix_free(J.alpha);
    gretl_matrix_free(J.Pi);
    gretl_matrix_free(J.lsPi);
    gretl_matrix_free(J.Omega);
    gretl_matrix_free(J.iOmega);
    gretl_matrix_free(J.Vb);
    gretl_matrix_free(J.Va);
    gretl_matrix_free(J.bse);
    gretl_matrix_free(J.ase);

    gretl_matrix_free(J.I00);
    gretl_matrix_free(J.I11);

    gretl_matrix_free(J.phi);
    gretl_matrix_free(J.psi);
    gretl_matrix_free(J.theta);

    gretl_matrix_free(J.qf1);
    gretl_matrix_free(J.qf2);
    gretl_matrix_free(J.Tmppp);
    gretl_matrix_free(J.Tmprp);

    gradhelper_free(J.ghelper);

    free(J.normrow);
    free(J.normcol);
    free(J.normval);

    free(J);
}

private Jwrap * jwrap_new(const(GretlVar) * jvar, Gretlopt opt, int * err)
{
    Jwrap * J = cast(Jwrap *) malloc(Jwrap.sizeof);

    if (J == null) {
        return null;
    }

    J.flags = 0;
    if (opt & Gretlopt.OPT_L) {
        J.flags |= JFlags.use_lbfgs;
    }

    J.T = jvar.T;
    J.p = jvar.neqns;
    J.r = jrank(jvar);
    J.blen = 0;
    J.alen = 0;
    J.df = 0;
    J.jr = 0;

    J.llk = J.T * 0.5 * (J.p * (1.0 + LN_2_PI));
    J.ll = NADBL;

    J.R0 = J.R1 = null;

    J.S00 = null;
    J.S01 = null;
    J.S11 = null;

    J.S00i = null;

    J.H = null;
    J.h = null;
    J.H0 = null;

    J.bigR = null;
    J.bigq = null;

    J.G = null;

    J.beta = null;
    J.alpha = null;
    J.Pi = null;
    J.lsPi = null;
    J.Omega = null;
    J.iOmega = null;
    J.Vb = null;
    J.Va = null;
    J.bse = null;
    J.ase = null;

    J.I00 = null;
    J.I11 = null;

    J.phi = null;
    J.psi = null;
    J.theta = null;

    J.qf1 = null;
    J.qf2 = null;
    J.Tmppp = null;
    J.Tmprp = null;

    J.ghelper = null;

    J.normrow = null;
    J.normcol = null;
    J.normval = null;

    *err = make_moment_matrices(J, jvar);

    if (!*err) {
        *err = J_alloc_aux(J);
    }

    if (*err) {
        jwrap_destroy(J);
        J = null;
    }

    return J;
}

/* vectorize src into targ */

private void vec_simple(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int i, n = src.rows * src.cols;

    for (i=0; i<n; i++) {
        targ.val[i] = src.val[i];
    }
}

/* vectorize the transpose of src into targ */

private void vec_transpose(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int i, j, k = 0;

    for (i=0; i<src.rows; i++) {
        for (j=0; j<src.cols; j++) {
            targ.val[k++] = gretl_matrix_get(src, i, j);
        }
    }
}

/* unvectorize src into targ, transposing */

private void unvec_transpose(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src)
{
    int i, j, k = 0;

    for (i=0; i<targ.rows; i++) {
        for (j=0; j<targ.cols; j++) {
            gretl_matrix_set(targ, i, j, src.val[k++]);
        }
    }
}

/* produce the \beta matrix in the case where it is fully
   constrained (and therefore no estimation is needed)
*/

private int solve_for_beta(Jwrap * J,
                           const(GretlMatrixRaw) * R,
                           const(GretlMatrixRaw) * q)
{
    GretlMatrixRaw * b = gretl_matrix_copy(q);
    int err = 0;

    if (b == null) {
        return GretlError.alloc;
    }

    if (!gretl_is_identity_matrix(R)) {
        GretlMatrixRaw * Rcpy = gretl_matrix_copy(R);

        if (Rcpy == null) {
            err = GretlError.alloc;
        } else {
            err = gretl_LU_solve(Rcpy, b);
            gretl_matrix_free(Rcpy);
        }
    }

    if (!err) {
        err = gretl_matrix_copy_values_shaped(J.beta, b);
    }

    if (!err) {
        J.blen = 0;
    }

    gretl_matrix_free(b);

    return err;
}

private int psi_from_alpha(Jwrap * J)
{
    GretlMatrixRaw * GG = null;
    GretlMatrixRaw * GGG = null;
    int err = 0;

    if (J.G == null) {
        /* just vectorize alpha' into psi */
        vec_transpose(J.psi, J.alpha);
        return 0;
    }

    GG = gretl_matrix_alloc(J.G.cols, J.G.cols);
    GGG = gretl_matrix_alloc(J.G.cols, J.G.rows);

    if (GG == null || GGG == null) {
        gretl_matrix_free(GG);
        gretl_matrix_free(GGG);
        return GretlError.alloc;
    }

    gretl_matrix_multiply_mod(J.G, GretlMatrixMod.transpose,
                              J.G, GretlMatrixMod.none,
                              GG, GretlMatrixMod.none);

    err = gretl_invert_symmetric_matrix(GG);
    if (err) {
        goto bailout;
    }

    gretl_matrix_multiply_mod(GG, GretlMatrixMod.none,
                              J.G, GretlMatrixMod.transpose,
                              GGG, GretlMatrixMod.none);

    gretl_matrix_reuse(J.Tmprp, J.r * J.p, 1);

    /* vectorize alpha' into Tmprp */
    vec_transpose(J.Tmprp, J.alpha);

    gretl_matrix_multiply(GGG, J.Tmprp, J.psi);
    gretl_matrix_reuse(J.Tmprp, J.r, J.p);

 bailout:

    gretl_matrix_free(GG);
    gretl_matrix_free(GGG);

    return err;
}

struct Switcher {
    GretlMatrixBlock * B; /* holder for the following */
    GretlMatrixRaw * K1;      /* holds kronecker product */
    GretlMatrixRaw * K2;      /* holds kronecker product */
    GretlMatrixRaw * TmpL;    /* used only in \phi calculation */
    GretlMatrixRaw * TmpR;    /* used only in \phi calculation */
    GretlMatrixRaw * TmpR1;   /* used when \alpha is restricted */
    GretlMatrixRaw * Tmprp1;  /* r x p1 temp */
}

private void switcher_free(Switcher * s)
{
    gretl_matrix_block_destroy(s.B);

    gretl_matrix_free(s.TmpL);
    gretl_matrix_free(s.TmpR1);
}

private int make_lsPi(Jwrap * J)
{
    GretlMatrixRaw * S11i;
    int err = 0;

    if (J.lsPi != null) {
        /* already done */
        return 0;
    }

    S11i = gretl_matrix_copy(J.S11);
    if (S11i == null) {
        return GretlError.alloc;
    }

    J.lsPi = gretl_matrix_alloc(J.p1, J.p);
    if (J.lsPi == null) {
        gretl_matrix_free(S11i);
        return GretlError.alloc;
    }

    err = gretl_invert_symmetric_matrix(S11i);

    if (!err) {
        gretl_matrix_multiply_mod(S11i, GretlMatrixMod.none,
                                  J.S01, GretlMatrixMod.transpose,
                                  J.lsPi, GretlMatrixMod.none);
        /* make into vec(\Pi'_{LS}) */
        gretl_matrix_reuse(J.lsPi, J.p1 * J.p, 1);
    }

    gretl_matrix_free(S11i);

    return err;
}

private int switcher_init(Switcher * s, Jwrap * J)
{
    int err = 0;

    s.TmpL = null;
    s.TmpR1 = null;

    s.B = gretl_matrix_block_new(&s.K1,   J.p1 * J.r, J.p1 * J.r,
                                 &s.K2,   J.p1 * J.r, J.p1 * J.p1,
                                 &s.TmpR, J.p * J.p1, 1,
                                 &s.Tmprp1, J.r, J.p1,
                                 null);
    if (s.B == null) {
        return GretlError.alloc;
    }

    if (J.blen > 0) {
        s.TmpL = gretl_matrix_alloc(J.blen, J.p * J.p1);
        if (s.TmpL == null) {
            return GretlError.alloc;
        }
    }

    if (J.G != null) {
        s.TmpR1 = gretl_matrix_alloc(J.alen, J.p * J.p1);
        if (s.TmpR1 == null) {
            return GretlError.alloc;
        }
    }

    err = make_lsPi(J);

    if (!err) {
        err = allocate_info_blocks_etc(J);
    }

    return err;
}

/* The following functions implement the switching algorithm
   as set out in Boswijk and Doornik, 2004, p. 455.
*/

private void alpha_from_psi(Jwrap * J)
{
    int i, j, k = 0;

    if (J.G != null) {
        gretl_matrix_reuse(J.Tmprp, J.p * J.r, 1);
        gretl_matrix_multiply(J.G, J.psi, J.Tmprp);
        for (i=0; i<J.p; i++) {
            for (j=0; j<J.r; j++) {
                gretl_matrix_set(J.alpha, i, j, J.Tmprp.val[k++]);
            }
        }
        gretl_matrix_reuse(J.Tmprp, J.r, J.p);
    } else {
        for (i=0; i<J.p; i++) {
            for (j=0; j<J.r; j++) {
                gretl_matrix_set(J.alpha, i, j, J.psi.val[k++]);
            }
        }
    }
}

/*
   Update \psi using:

   [ G'(\Omega^{-1} \otimes \beta'S_{11}\beta)G ]^{-1}
   \times G'(\Omega^{-1} \otimes \beta'S_{11}) vec(Pi'_{LS})

   or, in case \alpha is not restricted, just compute the
   ML alpha conditional on \beta in the usual Johansen
   manner.
*/

private int update_psi(Jwrap * J, Switcher * s)
{
    int err = 0;

    if (J.G == null) {
        J_compute_alpha(J);
        return 0;
    }

    gretl_matrix_reuse(s.K1, J.p * J.r, J.p * J.r);
    gretl_matrix_reuse(s.K2, J.p * J.r, J.p * J.p1);

    /* left-hand chunk */
    gretl_matrix_qform(J.beta, GretlMatrixMod.transpose, J.S11,
                       J.qf1, GretlMatrixMod.none);
    gretl_matrix_kronecker_product(J.iOmega, J.qf1, s.K1);
    gretl_matrix_qform(J.G, GretlMatrixMod.transpose, s.K1,
                       J.I00, GretlMatrixMod.none);

    /* right-hand chunk */
    gretl_matrix_multiply_mod(J.beta, GretlMatrixMod.transpose,
                              J.S11, GretlMatrixMod.none,
                              s.Tmprp1, GretlMatrixMod.none);
    gretl_matrix_kronecker_product(J.iOmega, s.Tmprp1, s.K2);
    gretl_matrix_multiply_mod(J.G, GretlMatrixMod.transpose,
                              s.K2, GretlMatrixMod.none,
                              s.TmpR1, GretlMatrixMod.none);
    gretl_matrix_multiply(s.TmpR1, J.lsPi, J.psi);

    /* combine */
    err = gretl_cholesky_decomp_solve(J.I00, J.psi);

    /* update alpha */
    if (!err) {
        alpha_from_psi(J);
    }

    return err;
}

/* update \beta based on \phi vector */

private void beta_from_phi(Jwrap * J)
{
    if (J.blen == 0) {
        return;
    }

    if (J.H != null) {
        gretl_matrix_reuse(J.beta, J.p1 * J.r, 1);
        gretl_matrix_multiply(J.H, J.phi, J.beta);
        if (!gretl_is_zero_matrix(J.h)) {
            gretl_matrix_add_to(J.beta, J.h);
        }
        gretl_matrix_reuse(J.beta, J.p1, J.r);
    } else {
        gretl_matrix_copy_values_shaped(J.beta, J.phi);
    }
}

/*
    Update \phi using:

    [H'(\alpha'\Omega^{-1}\alpha \otimes S_{11})H]^{-1} \times
       H'(\alpha'\Omega^{-1} \otimes S_{11}) \times
         [vec(Pi'_{LS}) - (\alpha \otimes I_{p1})h_0]
 */

private int update_phi(Jwrap * J, Switcher * s)
{
    int err = 0;

    if (J.blen == 0) {
        return 0;
    }

    gretl_matrix_reuse(s.K1, J.r * J.p1, J.r * J.p1);
    gretl_matrix_reuse(s.K2, J.r * J.p1, J.p * J.p1);

    /* first big inverse */
    gretl_matrix_qform(J.alpha, GretlMatrixMod.transpose, J.iOmega,
                       J.qf1, GretlMatrixMod.none);
    gretl_matrix_kronecker_product(J.qf1, J.S11, s.K1);
    if (J.H != null) {
        gretl_matrix_qform(J.H, GretlMatrixMod.transpose, s.K1,
                           J.I11, GretlMatrixMod.none);
    } else {
        gretl_matrix_copy_values(J.I11, s.K1);
    }

    /* second chunk */
    gretl_matrix_multiply_mod(J.alpha, GretlMatrixMod.transpose,
                              J.iOmega, GretlMatrixMod.none,
                              J.Tmprp, GretlMatrixMod.none);
    gretl_matrix_kronecker_product(J.Tmprp, J.S11, s.K2);
    if (J.H != null) {
        gretl_matrix_multiply_mod(J.H, GretlMatrixMod.transpose,
                                  s.K2, GretlMatrixMod.none,
                                  s.TmpL, GretlMatrixMod.none);
    } else {
        gretl_matrix_copy_values(s.TmpL, s.K2);
    }

    /* combine first and second chunks */
    err = gretl_cholesky_decomp_solve(J.I11, s.TmpL);
    if (err) {
        fprintf(stderr, "cholesky decomp failed in update_phi\n");
    }

    if (!err) {
        /* right-hand chunk */
        gretl_matrix_copy_values(s.TmpR, J.lsPi);
        if (J.h != null && !gretl_is_zero_matrix(J.h)) {
            gretl_matrix_reuse(s.K2, J.p * J.p1, J.r * J.p1);
            gretl_matrix_kronecker_I(J.alpha, J.p1, s.K2);
            gretl_matrix_multiply_mod(s.K2, GretlMatrixMod.none,
                                      J.h, GretlMatrixMod.none,
                                      s.TmpR, GretlMatrixMod.decrement);
        }

        /* combine */
        gretl_matrix_multiply(s.TmpL, s.TmpR, J.phi);

        /* update beta */
        beta_from_phi(J);
    }

    return err;
}

enum OmegaCode {
    only,
    plus
}

/*
    Update Omega using:

    S_{00} - S_{01} \beta\alpha' - \alpha\beta' S_{10} +
       \alpha\beta' S_{11} \beta\alpha'

    then invert into iOmega if wanted.
 */

private int make_Omega(Jwrap * J, int code)
{
    int err = 0;

    gretl_matrix_copy_values(J.Omega, J.S00);

    gretl_matrix_multiply_mod(J.alpha, GretlMatrixMod.none,
                              J.beta, GretlMatrixMod.transpose,
                              J.Pi, GretlMatrixMod.none);

    gretl_matrix_multiply_mod(J.S01, GretlMatrixMod.none,
                              J.Pi, GretlMatrixMod.transpose,
                              J.Tmppp, GretlMatrixMod.none);

    gretl_matrix_add_self_transpose(J.Tmppp);
    gretl_matrix_subtract_from(J.Omega, J.Tmppp);

    gretl_matrix_qform(J.Pi, GretlMatrixMod.none, J.S11,
                       J.Omega, GretlMatrixMod.cumulate);

    if (code == OmegaCode.plus) {
        gretl_matrix_copy_values(J.iOmega, J.Omega);
        err = gretl_invert_symmetric_matrix(J.iOmega);
    }

    return err;
}

/* reduced version of log-likelihood calculation for
   switching algorithm loop */

private int switcher_ll(Jwrap * J)
{
    int err = 0;

    gretl_matrix_copy_values(J.Tmppp, J.Omega);
    J.ll = gretl_matrix_log_determinant(J.Tmppp, &err);
    if (!err) {
        J.ll *= -J.T / 2.0;
    }

    return err;
}

private int form_I00_inverse(Jwrap * J)
{
    GretlMatrixRaw * K = null;
    int err = 0;

    if (J.I00 == null) {
        J.I00 = gretl_matrix_alloc(J.alen, J.alen);
        if (J.I00 == null) {
            return GretlError.alloc;
        }
    }

    gretl_matrix_qform(J.beta, GretlMatrixMod.transpose, J.S11,
                       J.qf1, GretlMatrixMod.none);

    K = gretl_matrix_kronecker_product_new(J.iOmega, J.qf1, &err);

    if (!err) {
        if (J.G != null) {
            gretl_matrix_qform(J.G, GretlMatrixMod.transpose, K,
                               J.I00, GretlMatrixMod.none);
        } else {
            gretl_matrix_copy_values(J.I00, K);
        }
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(J.I00);
    }

    gretl_matrix_free(K);

    return err;
}

private int form_I11_inverse(Jwrap * J)
{
    GretlMatrixRaw * K = null;
    int err = 0;

    if (J.blen == 0) {
        return 0;
    }

    if (J.I11 == null) {
        J.I11 = gretl_matrix_alloc(J.blen, J.blen);
        if (J.I11 == null) {
            return GretlError.alloc;
        }
    }

    gretl_matrix_qform(J.alpha, GretlMatrixMod.transpose, J.iOmega,
                       J.qf1, GretlMatrixMod.none);

    K = gretl_matrix_kronecker_product_new(J.qf1, J.S11, &err);

    if (!err) {
        const(GretlMatrixRaw) * H = (do_scaling(J))? J.H0 : J.H;

        if (H != null) {
            gretl_matrix_qform(H, GretlMatrixMod.transpose, K,
                               J.I11, GretlMatrixMod.none);
        } else {
            gretl_matrix_copy_values(J.I11, K);
        }
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(J.I11);
    }

    gretl_matrix_free(K);

    return err;
}

private int add_Omega_inverse(Jwrap * J)
{
    int err = 0;

    J.iOmega = gretl_matrix_copy(J.Omega);

    if (J.iOmega == null) {
        err = GretlError.alloc;
    } else {
        err = gretl_invert_symmetric_matrix(J.iOmega);
    }

    return err;
}

private void restricted_vecm_set_df(Jwrap * J, GretlVar * v)
{
    int K = 0;
    double c;

    /* lagged differences */
    K += v.order * J.p;

    /* deterministic stuff */
    K += v.jinfo.seasonals + v.ifc;
    if (jcode(v) == JohansenCode.unrest_trend) {
        K++;
    }

    /* exogenous vars? */
    if (v.xlist != null) {
        K += v.xlist[0];
    }

    K *= J.p;

    /* add number of free beta and alpha terms */
    K += J.blen + J.alen;

    c = K / cast(double) J.p;
    v.df = v.T - cast(int) floor(c);
}

/* Use the information matrix to compute standard errors for
   beta and alpha.
*/

private int variance_from_info_matrix(Jwrap * J, GretlVar * v)
{
    int npar = J.alen + J.blen;
    int err = 0;

    if (J.jr < npar) {
        /* model is not fully identified */
        return 0;
    }

    /* preliminary */

    restricted_vecm_set_df(J, v);

    if (J.iOmega == null) {
        err = add_Omega_inverse(J);
        if (err) return err;
    }

    /* variance of beta */

    if (J.blen == 0) {
        int nb = J.p1 * J.r;

        J.Vb = gretl_zero_matrix_new(nb, nb);
        J.bse = gretl_zero_matrix_new(J.p1, J.r);

        if (J.Vb == null || J.bse == null) {
            err = GretlError.alloc;
        }

        goto beta_done;
    }

    if (do_scaling(J)) {
        gretl_matrix_free(J.I11);
        J.I11 = null;
    }

    err = form_I11_inverse(J);
    if (err) return err;

    gretl_matrix_divide_by_scalar(J.I11, v.df);

    if (J.H0 != null) {
        int nb = J.r * J.p1;

        J.Vb = gretl_matrix_alloc(nb, nb);
        if (J.Vb == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_qform(J.H0, GretlMatrixMod.none, J.I11,
                               J.Vb, GretlMatrixMod.none);
        }
    } else {
        J.Vb = J.I11;
        J.I11 = null;
    }

    if (!err) {
        err = make_beta_se(J);
    }

 beta_done:

    if (err) {
        return err;
    }

    /* variance of alpha */

    err = form_I00_inverse(J);
    if (err) return err;

    gretl_matrix_divide_by_scalar(J.I00, v.df);

    if (J.G != null) {
        int na = J.r * J.p;

        J.Va = gretl_matrix_alloc(na, na);
        if (J.Va == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_qform(J.G, GretlMatrixMod.none, J.I00,
                               J.Va, GretlMatrixMod.none);
        }
    } else {
        J.Va = J.I00;
        J.I00 = null;
    }

    if (!err) {
        err = make_alpha_se(J);
    }

    return err;
}

/* preliminary test for whether we're using manual
   initialization with the switching algorithm: if so,
   we won't try to extract column scale factors
*/

private int doing_manual_init(Jwrap * J)
{
    int ulen = n_initvals();
    int ret = 0;

    if (ulen > 0) {
        int psilen = (J.G == null)? 0 : J.alen;

        if (ulen == J.blen + psilen) {
            ret = 1;
        } else if (ulen == J.blen + J.alen) {
            ret = 1;
        } else if (ulen == (J.p1 + J.p) * J.r) {
            ret = 1;
        }
    }

    return ret;
}

private int user_switch_init(Jwrap * J, int * uinit)
{
    GretlMatrixRaw * U;
    int psilen;
    int ulen, i, k;
    int err = 0;

    U = get_initvals();
    if (U == null) {
        return 0;
    }

    psilen = (J.G == null)? 0 : J.alen;
    ulen = gretl_vector_get_length(U);

    if (ulen == J.blen + psilen) {
        /* we're given phi and psi? */
        k = 0;
        for (i=0; i<J.blen; i++) {
            J.phi.val[i] = U.val[k++];
        }
        if (psilen == 0) {
            beta_from_phi(J);
            J_compute_alpha(J);
        } else {
            for (i=0; i<psilen; i++) {
                J.psi.val[i] = U.val[k++];
            }
        }
    } else if (ulen == J.blen + J.alen) {
        /* we're given phi and (redundant) psi? */
        k = 0;
        for (i=0; i<J.blen; i++) {
            J.phi.val[i] = U.val[k++];
        }
        for (i=0; i<J.alen; i++) {
            J.psi.val[i] = U.val[k++];
        }
        alpha_from_psi(J);
    } else if (ulen == (J.p1 + J.p) * J.r) {
        /* we're given vec(beta), vec(alpha)? */
        int n = J.p1 * J.r;

        k = 0;
        for (i=0; i<n; i++) {
            J.beta.val[i] = U.val[k++];
        }

        n = J.p * J.r;
        for (i=0; i<n; i++) {
            J.alpha.val[i] = U.val[k++];
        }

        err = phi_from_beta(J);
        if (!err) {
            err = psi_from_alpha(J);
        }
    } else {
        gretl_errmsg_set("Invalid initialization given");
        err = GretlError.data;
    }

    gretl_matrix_free(U);

    *uinit = 1;

    return err;
}

/* driver function for switching algorithm */

private int switchit(Jwrap * J, Gretlopt opt, Prn * prn)
{
    Switcher s;
    double lldiff = NADBL;
    double llbak = -1.0e+200;
    double eps1 = 0.0001;
    double eps2 = 0.005;
    double stol = 0.01 * eps1;
    double wtol = 0.01 * eps2;
    int j, jmax = 50000;
    int wcount = 0;
    int uinit = 0;
    int conv = 0;
    int err;

    /* old: tol = 2.0e-11 */

    err = user_switch_init(J, &uinit);
    if (err) {
        return err;
    }

    err = switcher_init(&s, J);

    if (!err && J.H != null) {
        beta_from_phi(J);
    }

    if (!err && J.G != null) {
        alpha_from_psi(J);
    }

    if (!err) {
        /* initialize Omega */
        err = make_Omega(J, OmegaCode.plus);
    }

    for (j=0; j<jmax && !err; j++) {
        err = update_phi(J, &s);
        if (!err) {
            err = update_psi(J, &s);
        }
        if (!err) {
            err = make_Omega(J, OmegaCode.plus);
        }
        if (!err) {
            err = switcher_ll(J);
        }

        if (!err && j > 1) {
            lldiff = J.ll - llbak;
            if (lldiff < stol) {
                conv = 1;
                break;
            } else if (lldiff < wtol) {
                if (++wcount == 4) {
                    conv = 2;
                    break;
                }
            } else {
                wcount = 0;
            }
        }

        llbak = J.ll;
    }

    if (err) {
        pputs(prn, _("Switching algorithm: failed").ptr);
        pputc(prn, '\n');
    } else if (!(opt & Gretlopt.OPT_S)) {
        /* not silent */
        pprintf(prn, _("Switching algorithm: %d iterations").ptr, j);
        if (uinit) {
            pprintf(prn, " (%s)", _("user-supplied initial values").ptr);
        }
        pputc(prn, '\n');

        if (!na(J.ll) && !na(lldiff)) {
            pprintf(prn, " -(T/2)log|Omega| = %.8g, lldiff = %g\n", J.ll, lldiff);
        }

        if (conv == 1) {
            pprintf(prn, "%s\n", _("Strong convergence").ptr);
        } else if (conv == 2) {
            pprintf(prn, "%s\n", _("Weak convergence").ptr);
        }
    }

    if (!err && !conv) {
        err = GretlError.noconv;
    } else if (j > 1500) {
        pputs(prn, "*** warning: large number of iterations may "
              ~ "indicate a problem");
        pputc(prn, '\n');
    }

    if (!(opt & Gretlopt.OPT_S)) {
        pputc(prn, '\n');
    }

    if (!err) {
        J.ll -= J.llk;
    }

    switcher_free(&s);

    return err;
}

/* re-scale beta and alpha based on the saved scaling
   information
*/

private int replace_col_scaling(Jwrap * J)
{
    double x, s;
    int i, j, k, ii;
    int err = 0;

    for (i=0; i<J.normcol[0]; i++) {
        j = J.normcol[i+1];
        if (j < 0) {
            continue;
        }
        k = j / J.p1;
        s = gretl_matrix_get(J.beta, j % J.p1, k) / J.normval[i];
        for (ii=0; ii<J.p1; ii++) {
            x = gretl_matrix_get(J.beta, ii, k);
            if (x != 0) {
                gretl_matrix_set(J.beta, ii, k, x / s);
            }
        }
        for (ii=0; ii<J.p; ii++) {
            x = gretl_matrix_get(J.alpha, ii, k);
            if (x != 0) {
                gretl_matrix_set(J.alpha, ii, k, x * s);
            }
        }
        J.blen -= 1;
    }

    return err;
}

private int reset_null_H(Jwrap * J)
{
    J.H0 = J.H;
    gretl_matrix_free(J.h);
    J.H = J.h = null;

    J.blen = J.r * J.p1;

    gretl_matrix_free(J.phi);
    J.phi = gretl_zero_matrix_new(J.blen, 1);
    if (J.phi == null) {
        return GretlError.alloc;
    }

    return 0;
}

/* recalculate the H matrix after removing column
   scalings */

private int reset_H(Jwrap * J, const(GretlMatrixRaw) * R,
                    const(GretlMatrixRaw) * q)
{
    int err;

    J.H0 = J.H;
    gretl_matrix_free(J.h);
    J.h = null;

    err = real_set_up_H(J, R, q);

    if (!err) {
        gretl_matrix_free(J.phi);
        J.phi = gretl_zero_matrix_new(J.blen, 1);
        if (J.phi == null) {
            err = GretlError.alloc;
        }
    }

    return err;
}

private const(GretlMatrixRaw) * J_get_R_matrix(Jwrap * J, const(GretlRestriction) * rset)
{
    if (J.bigR != null) {
        return J.bigR;
    } else {
        return rset_get_R_matrix(rset);
    }
}

private const(GretlMatrixRaw) * J_get_q_matrix(Jwrap * J, const(GretlRestriction) * rset)
{
    if (J.bigq != null) {
        return J.bigq;
    } else {
        return rset_get_q_matrix(rset);
    }
}

private int get_master_col(Jwrap * J, int pos, int * mpos)
{
    int i, k = -J.normcol[pos];
    int ci, bcol = k / J.p1;

    for (i=0; i<J.normcol[0]; i++) {
        ci = J.normcol[i+1];
        if (ci >= 0 && ci / J.p1 == bcol) {
            *mpos = i + 1;
            return ci;
        }
    }

    return -1;
}

/* Take a row from the original R matrix and re-write it as a
   homogenous restriction referencing the "master" value which
   was used as a basis for scaling the relevant column of beta.
   For example, if we found the scaling b[2,1] = 1, and then the
   following restriction b[2,2] = -1, we remove the first
   restriction (temporarily), and we replace the second with
   b[1,1] + b[2,1] = 0.
*/

private int homogenize_R_line(GretlMatrixRaw * R,
                              GretlMatrixRaw * q,
                              Jwrap * J, int pos,
                              int k)
{
    int j, j0, j1, mpos;
    double x;

    j0 = get_master_col(J, pos, &mpos);
    if (j0 < 0) {
        return GretlError.data;
    }

    j1 = -J.normcol[pos];
    x = -J.normval[mpos-1] / J.normval[pos-1];

    for (j=0; j<R.cols; j++) {
        if (j == j0) {
            gretl_matrix_set(R, k, j, 1);
        } else if (j == j1) {
            gretl_matrix_set(R, k, j, x);
        } else {
            gretl_matrix_set(R, k, j, 0);
        }
    }

    gretl_vector_set(q, k, 0);

    return 0;
}

private void norm_destroy(Jwrap * J)
{
    free(J.normrow);
    free(J.normcol);
    free(J.normval);

    J.normrow = null;
    J.normcol = null;
    J.normval = null;
}

private int norm_allocate(Jwrap * J, const(GretlMatrixRaw) * R)
{
    int n = R.rows;
    int err = 0;

    J.normrow = gretl_list_new(n);
    if (J.normrow == null) err = GretlError.alloc;

    if (!err) {
        J.normcol = gretl_list_new(n);
        if (J.normcol == null) err = GretlError.alloc;
    }

    if (!err) {
        J.normval = cast(double *) malloc(n * (*J.normval).sizeof);
        if (J.normval == null) err = GretlError.alloc;
    }

    if (err) {
        norm_destroy(J);
    } else {
        J.normrow[0] = 0;
        J.normcol[0] = 0;
    }

    return err;
}

private void list_push(int * list, int k)
{
    list[0] += 1;
    list[list[0]] = k;
}

/* test whether or not a row of the R matrix (taking
   into account the corresponding q value) represents
   a scaling of a column of beta
*/

private int is_scaling_row(const(GretlMatrixRaw) * R,
                           const(GretlMatrixRaw) * q,
                           int i, int * col,
                           double * sval)
{
    int j, c = 0;
    double rij;

    for (j=0; j<R.cols; j++) {
        rij = gretl_matrix_get(R, i, j);
        if (rij != 0 && q.val[i] != 0) {
            c++;
            if (c > 1) {
                break;
            } else {
                *sval = q.val[i] / rij;
                *col = j;
            }
        }
    }

    return c == 1;
}

/* Check for restrictions which make scale removal impossible,
   for each column of beta.  Record a -1 in the appropriate
   place in the @sc array if scaling won't work for that column.

   We consider scale removal infeasible for a given beta column
   if (a) we find a non-homogeneous restriction involving two
   or more beta elements (e.g. b[1,1] - b[1,3] = 2), or (b) we
   find any cross-column restrictions.

   Return 1 if there are no scaleable columns on this test, else 0.
*/

private int screen_beta_cols(Jwrap * J,
                             const(GretlMatrixRaw) * R,
                             const(GretlMatrixRaw) * q,
                             byte * sc)
{
    double x;
    int ns = 0;
    int i, j, k, m, c;

    for (i=0; i<R.rows; i++) {
        k = -1;
        c = 0;
        for (j=0; j<R.cols; j++) {
            x = gretl_matrix_get(R, i, j);
            if (x != 0) {
                c++;
                if (k < 0) {
                    k = j / J.p1;
                } else {
                    m = j / J.p1;
                    if (m != k) {
                        /* cross-column restriction */
                        sc[k] = sc[m] = -1;
                    }
                }
            }
        }
        if (c > 1 && q.val[i] != 0) {
            sc[k] = -1;
        }
    }

    for (i=0; i<J.r; i++) {
        if (sc[i] == 0) ns++;
    }

    return (ns == 0);
}

/* Check for restrictions (rows of R) which are just scalings of a
   given column of beta.  There can only be one such scaling per
   beta column.
*/

private int check_for_scaling(Jwrap * J,
                              const(GretlMatrixRaw) * R,
                              const(GretlMatrixRaw) * q)
{
    byte * sc = cast(byte *) calloc(J.r, 1);
    double sval = 0;
    int i, k, rcol = 0;

    if (sc == null) {
        return GretlError.alloc;
    }

    /* first pass: screen for conditions which make
       scale-removal impossible */
    if (screen_beta_cols(J, R, q, sc)) {
        free(sc);
        return 0;
    }

    if (norm_allocate(J, R)) {
        return GretlError.alloc;
    }

    for (i=0; i<R.rows; i++) {
        if (is_scaling_row(R, q, i, &rcol, &sval)) {
            k = rcol / J.p1;
            if (sc[k] >= 0) {
                /* column is scaleable */
                J.normval[J.normrow[0]] = sval;
                list_push(J.normrow, i);
                if (sc[k] == 0) {
                    /* beta column k is not yet scaled */
                    list_push(J.normcol, rcol);
                } else {
                    /* negative rcol entry indicates a "follower" */
                    list_push(J.normcol, -rcol);
                }
                sc[k] += 1;
            }
        }
    }

    if (J.normrow[0] == 0) {
        /* no (practicable) scalings found */
        norm_destroy(J);
    }

    free(sc);

    return 0;
}

private int just_identified(Jwrap * J)
{
    int npar = J.alen + J.blen;

    return (J.jr >= npar && J.df == 0);
}

private int maybe_remove_col_scaling(Jwrap * J,
                                     const(GretlRestriction) * rset)
{
    const(GretlMatrixRaw) * R = J_get_R_matrix(J, rset);
    const(GretlMatrixRaw) * q = J_get_q_matrix(J, rset);
    GretlMatrixRaw * Rr = null;
    GretlMatrixRaw * qr = null;
    double x;
    int i, j, k, rr, pos;
    int err = 0;

    if (R == null || q == null || gretl_is_zero_matrix(q)) {
        return 0;
    }

    if (using_bfgs(J) || cross_alpha(J) || doing_manual_init(J)) {
        return 0;
    }

    /* ?? */
    if (just_identified(J)) {
        return 0;
    }

    err = check_for_scaling(J, R, q);
    if (err || !do_scaling(J)) {
        return err;
    }

    rr = R.rows;
    for (i=0; i<J.normcol[0]; i++) {
        if (J.normcol[i+1] >= 0) {
            rr--;
        }
    }

    if (rr == 0) {
        /* all the restrictions were column scalings */
        return reset_null_H(J);
    }

    Rr = gretl_matrix_alloc(rr, R.cols);
    qr = gretl_column_vector_alloc(rr);

    if (Rr == null || qr == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    k = 0;
    for (i=0; i<R.rows; i++) {
        pos = in_gretl_list(J.normrow, i);
        if (pos > 0 && J.normcol[pos] < 0) {
            /* convert to homogeneous restriction */
            homogenize_R_line(Rr, qr, J, pos, k);
            k++;
        } else if (pos == 0) {
            /* pass the row through unmodified */
            for (j=0; j<R.cols; j++) {
                x = gretl_matrix_get(R, i, j);
                gretl_matrix_set(Rr, k, j, x);
                x = gretl_vector_get(q, i);
                gretl_vector_set(qr, k, x);
            }
            k++;
        }
    }

    err = reset_H(J, Rr, qr);

 bailout:

    gretl_matrix_free(Rr);
    gretl_matrix_free(qr);

    return err;
}

/*
   J = [(I_p \otimes \beta)G : (\alpha \otimes I_{p1})H]

   Boswijk and Doornik (2004), equation (40), page 455.
*/

private int check_jacobian(Jwrap * J)
{
    GretlMatrixRaw * L = null;
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * Jac = null;
    int err = 0;

    /* form beta and alpha, based on randomized theta */

    if (J.H != null) {
        GretlMatrixRaw * phi;

        if (J.H.rows != J.p1 * J.r) {
            fprintf(stderr, "?? matrices wrongly sized in check_jacobian\n");
            return GretlError.nonconf;
        }

        phi = gretl_column_vector_alloc(J.H.cols);
        if (phi == null) {
            return GretlError.alloc;
        }

        gretl_matrix_random_fill(phi, DistCode.normal);
        gretl_matrix_reuse(J.beta, J.p1 * J.r, 1);
        err = gretl_matrix_multiply(J.H, phi, J.beta);
        if (!err && J.h != null) {
            err = gretl_matrix_add_to(J.beta, J.h);
        }
        gretl_matrix_reuse(J.beta, J.p1, J.r);
        gretl_matrix_free(phi);
    } else if (J.blen > 0) {
        gretl_matrix_random_fill(J.beta, DistCode.normal);
    }

    if (!err && J.G != null) {
        GretlMatrixRaw * psi = gretl_column_vector_alloc(J.G.cols);
        GretlMatrixRaw * avec = gretl_column_vector_alloc(J.p * J.r);

        if (psi == null || avec == null) {
            err = GretlError.alloc;
            goto bailout;
        }

        gretl_matrix_random_fill(psi, DistCode.normal);
        gretl_matrix_multiply(J.G, psi, avec);
        unvec_transpose(J.alpha, avec);
        gretl_matrix_free(psi);
        gretl_matrix_free(avec);
    } else if (!err) {
        J_compute_alpha(J);
    }

    if (!err) {
        L = gretl_matrix_I_kronecker_new(J.p, J.beta, &err);
    }

    if (!err) {
        R = gretl_matrix_kronecker_I_new(J.alpha, J.p1, &err);
    }

    if (!err && J.G != null) {
        /* alpha is restricted */
        GretlMatrixRaw * LG = gretl_matrix_multiply_new(L, J.G, &err);

        if (!err) {
            gretl_matrix_replace(&L, LG);
        }
    }

    if (!err && J.H != null) {
        /* beta is restricted */
        GretlMatrixRaw * RH = gretl_matrix_multiply_new(R, J.H, &err);

        if (!err) {
            gretl_matrix_replace(&R, RH);
        }
    }

    if (!err) {
        Jac = gretl_matrix_col_concat(L, R, &err);
    }

    if (!err) {
        J.jr = gretl_matrix_rank(Jac, NADBL, &err);
    }

 bailout:

    gretl_matrix_free(L);
    gretl_matrix_free(R);
    gretl_matrix_free(Jac);

    return err;
}

/* Doornik's approach to checking identification */

private int vecm_id_check(Jwrap * J, GretlVar * jvar,
                          Gretlopt opt, Prn * prn)
{
    int npar = J.alen + J.blen;
    int silent = (opt & Gretlopt.OPT_S)? 1 : 0;
    int err = 0;

    err = check_jacobian(J);

    if (!err) {
        if (!silent) {
            pputc(prn, '\n');
            pprintf(prn, _("Rank of Jacobian = %d, number of free "
                           ~ "parameters = %d\n").ptr, J.jr, npar);
        }
        if (J.jr < npar) {
            if (!silent) {
                pputs(prn, _("Model is not fully identified\n").ptr);
            }
            if (using_bfgs(J)) {
                err = GretlError.noident;
            }
        } else if (!silent) {
            pputs(prn, _("Model is fully identified\n").ptr);
        }

        J.df = (J.p + J.p1 - J.r) * J.r - J.jr;
        if (!silent) {
            pprintf(prn, _("Based on Jacobian, df = %d\n").ptr, J.df);
        }

        if (!(opt & Gretlopt.OPT_B)) {
            /* system was subject to a prior restriction? */
            if (jvar.jinfo.lrdf > 0) {
                J.df -= jvar.jinfo.lrdf;
                if (!silent) {
                    pprintf(prn, _("Allowing for prior restriction, df = %d\n").ptr,
                            J.df);
                }
            }
        }
    }

    return err;
}

/* common alpha restriction needs to be replicated */

private int G_from_expanded_R(Jwrap * J, const(GretlMatrixRaw) * R)
{
    GretlMatrixRaw * Rtmp = null;
    GretlMatrixRaw * Rcpy = null;
    double x;
    int c0, c1 = 0;
    int i, j, k;
    int err = 0;

    Rtmp = gretl_matrix_I_kronecker_new(J.r, R, &err);

    if (!err) {
        Rcpy = gretl_matrix_copy(Rtmp);
        if (Rcpy == null) {
            err = GretlError.alloc;
        }
    }

    for (i=0; i<J.p; i++) {
        for (j=0; j<J.r; j++) {
            /* write col of orig R into col of remapped R */
            c0 = i + j * J.p;
            for (k=0; k<Rcpy.rows; k++) {
                x = gretl_matrix_get(Rcpy, k, c0);
                gretl_matrix_set(Rtmp, k, c1, x);
            }
            c1++;
        }
    }

    if (!err) {
        J.G = gretl_matrix_right_nullspace(Rtmp, &err);
    }

    gretl_matrix_free(Rtmp);
    gretl_matrix_free(Rcpy);

    return err;
}

enum NSMIN = 1.0e-16;

private void jrestrict_normalize_nullspace(GretlMatrixRaw * A)
{
    int i, j, nz, nm1, pos;
    double x, aij;

    if (A.cols == 1) {
        x = 0;
        for (i=0; i<A.rows; i++) {
            if (fabs(A.val[i]) > x) {
                x = A.val[i];
            }
        }
        for (i=0; i<A.rows; i++) {
            aij = A.val[i] / x;
            A.val[i] = fabs(aij) < NSMIN ? 0 : aij;
        }
    }

    for (i=0; i<A.rows; i++) {
        nz = nm1 = pos = 0;
        for (j=0; j<A.cols; j++) {
            aij = gretl_matrix_get(A, i, j);
            if (aij == 0.0) {
                nz++;
            } else if (aij == -1.0) {
                nm1++;
                pos = j;
            } else {
                break;
            }
        }
        if (nz == A.cols - 1 && nm1 == 1) {
            /* got a solitary -1 on row i */
            gretl_matrix_set(A, i, pos, 1.0);
        }
    }
}

/* user-supplied R needs to be remapped for conformity with
   vec(\alpha') = G\psi */

private int G_from_remapped_R(Jwrap * J, const(GretlMatrixRaw) * R)
{
    GretlMatrixRaw * Rtmp = null;
    double x;
    int c0, c1 = 0;
    int i, j, k;
    int err = 0;

    Rtmp = gretl_matrix_copy(R);
    if (Rtmp == null) {
        return GretlError.alloc;
    }

    for (i=0; i<J.p; i++) {
        for (j=0; j<J.r; j++) {
            /* write col of orig R into col of remapped R */
            c0 = i + j * J.p;
            for (k=0; k<R.rows; k++) {
                x = gretl_matrix_get(R, k, c0);
                gretl_matrix_set(Rtmp, k, c1, x);
            }
            c1++;
        }
    }

    J.G = gretl_matrix_right_nullspace(Rtmp, &err);
    gretl_matrix_free(Rtmp);

    return err;
}

/* Old comment: "detect the presence of a cross-equation alpha
   restriction". As of 2018-05-10 I'm unsure about this function.
   What it actually does is detect the case where @Ra, as a whole,
   places restrictions on more than one column of alpha, that is,
   on more than one EC term.

   @Ra will have as many rows as there are restrictions on
   alpha and as many columns as elements in vec(alpha). When
   we integer-divide the column index j by J->p (the number
   of equations in the system) we get the 0-based index of the
   associated alpha column. The code looks for non-zero
   elements in @Ra: on finding the first one it records the
   alpha column (@acol); if it later finds a non-zero entry
   associated with a different alpha column it records the
   "cross_alpha" flag and returns 1.

   My puzzlement now is whether we're really seeking "cross-
   equation alpha restrictions" -- if so, the function seems
   wrong... Note that the only actual use of the "cross_alpha"
   flag is to rule out scale removal in restricted estimation,
   so the real question is, what feature of @Ra ought to trigger
   that limitation?
*/

private int cross_alpha_check(Jwrap * J, const(GretlMatrixRaw) * Ra)
{
    double x;
    int i, j, acol = -1;

    for (i=0; i<Ra.rows; i++) {
        for (j=0; j<Ra.cols; j++) {
            x = gretl_matrix_get(Ra, i, j);
            if (x != 0) {
                if (acol < 0) {
                    /* first non-zero element found */
                    acol = j / J.p;
                } else if (j / J.p != acol) {
                    /* this non-zero element pertains
                       to a different alpha column
                    */
                    J.flags |= JFlags.cross_alpha;
                    return 1;
                }
            }
        }
    }

    return 0;
}

/* set up restriction matrix, G, for alpha */

private int set_up_G(Jwrap * J, const(GretlRestriction) * rset)
{
    const(GretlMatrixRaw) * Ra = rset_get_Ra_matrix(rset);
    const(GretlMatrixRaw) * qa = rset_get_qa_matrix(rset);
    int err = 0;

    if (qa != null && !gretl_is_zero_matrix(qa)) {
        gretl_errmsg_set("alpha restriction is not homogeneous: "
                         ~ "not supported");
        return GretlError.notimp;
    }

    if (J.r > 1) {
        if (Ra.cols == J.p) {
            /* got a common alpha restriction */
            err = G_from_expanded_R(J, Ra);
        } else if (Ra.cols > J.p) {
            err = G_from_remapped_R(J, Ra);
        }
    } else {
        J.G = gretl_matrix_right_nullspace(Ra, &err);
    }

    if (!err) {
        J.alen = J.G.cols;
    }

    if (!err) {
        jrestrict_normalize_nullspace(J.G);
    }

    if (!err) {
        cross_alpha_check(J, Ra);
    }

    return err;
}

private int real_set_up_H(Jwrap * J, const(GretlMatrixRaw) * R,
                          const(GretlMatrixRaw) * q)
{
    GretlMatrixRaw * RRT = null;
    GretlMatrixRaw * Tmp = null;
    int err = 0;

    J.H = gretl_matrix_right_nullspace(R, &err);

    if (err) {
        return err;
    }

    J.blen = J.H.cols;

    if (q == null || gretl_is_zero_matrix(q)) {
        /* implicitly or explicitly zero */
        J.h = gretl_zero_matrix_new(R.cols, 1);
        return (J.h == null)? GretlError.alloc : 0;
    }

    /* now for h, for non-zero q */

    RRT = gretl_matrix_alloc(R.rows, R.rows);
    Tmp = gretl_matrix_alloc(R.cols, R.rows);
    if (RRT == null || Tmp == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(R, GretlMatrixMod.none,
                                        R, GretlMatrixMod.transpose,
                                        RRT, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(RRT);
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(R, GretlMatrixMod.transpose,
                                        RRT, GretlMatrixMod.none,
                                        Tmp, GretlMatrixMod.none);
    }

    if (!err) {
        J.h = gretl_matrix_multiply_new(Tmp, q, &err);
    }

    gretl_matrix_free(RRT);
    gretl_matrix_free(Tmp);

    return err;
}

private GretlMatrixRaw * replicate_q(int r, const(GretlMatrixRaw) * q0,
                                     int * err)
{
    GretlMatrixRaw * q;
    int k = q0.rows * r;
    int i, j;

    q = gretl_column_vector_alloc(k);

    if (q == null) {
        *err = GretlError.alloc;
    } else {
        k = 0;
        for (j=0; j<r; j++) {
            for (i=0; i<q0.rows; i++) {
                q.val[k++] = q0.val[i];
            }
        }
    }

    return q;
}

/* set up restriction matrices, H and h, for beta */

private int set_up_H(Jwrap * J, const(GretlRestriction) * rset)
{
    const(GretlMatrixRaw) * R = rset_get_R_matrix(rset);
    const(GretlMatrixRaw) * q = rset_get_q_matrix(rset);
    int qzero, err = 0;

    if (R.rows == J.p1 * J.r) {
        /* number of restrictions = total betas */
        return solve_for_beta(J, R, q);
    }

    qzero = (q == null || gretl_is_zero_matrix(q));

    if (J.r > 1 && R.cols == J.p1) {
        /* common beta restriction */
        J.bigR = gretl_matrix_I_kronecker_new(J.r, R, &err);
        if (!err) {
            R = J.bigR;
            if (!qzero) {
                J.bigq = replicate_q(J.r, q, &err);
                if (!err) {
                    q = J.bigq;
                }
            }
        }
    }

    if (!err) {
        err = real_set_up_H(J, R, q);
    }

    if (!err) {
        J.H0 = J.H;
    }

    return err;
}

/* See H. Peter Boswijk, "Identifiability of Cointegrated Systems",
   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.500.8189
   We're assuming here that vec(\beta) = H\psi + h_0, for
   non-zero h_0.
*/

private int phi_init_nonhomog(Jwrap * J)
{
    GretlMatrixRaw * BP = null;
    GretlMatrixRaw * IBP = null;
    GretlMatrixRaw * IBPH = null;
    GretlMatrixRaw * IBPh = null;
    GretlMatrixRaw * e = null;
    int nzero = J.p1 - J.r;
    int err = 0;

    if (J.h == null || gretl_is_zero_matrix(J.h) ||
        nzero == 0 || J.blen == 0) {
        /* shouldn't be here! */
        return 0;
    }

    BP = gretl_matrix_alloc(J.p1, J.p1);
    IBPH = gretl_matrix_alloc(J.r * nzero, J.blen);
    IBPh = gretl_column_vector_alloc(J.r * nzero);

    if (BP == null || IBPH == null || IBPh == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    gretl_matrix_multiply_mod(J.beta, GretlMatrixMod.none,
                              J.beta, GretlMatrixMod.transpose,
                              BP, GretlMatrixMod.none);

    e = gretl_symmetric_matrix_eigenvals(BP, 1, &err);
    if (!err) {
        /* truncate @BP to just the eigenvectors corresponding
           to the @nzero (in principle) zero eigenvalues
        */
        BP.cols = nzero;
    }

    if (!err) {
        IBP = gretl_matrix_I_kronecker_new(J.r, BP, &err);
    }

    if (!err) {
        gretl_matrix_multiply_mod(IBP, GretlMatrixMod.transpose,
                                  J.H, GretlMatrixMod.none,
                                  IBPH, GretlMatrixMod.none);
        gretl_matrix_multiply_mod(IBP, GretlMatrixMod.transpose,
                                  J.h, GretlMatrixMod.none,
                                  IBPh, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_matrix_moore_penrose(IBPH, NADBL);
    }

    if (!err) {
        gretl_matrix_multiply(IBPH, IBPh, J.phi);
        if (gretl_is_zero_matrix(J.phi)) {
            fprintf(stderr, "Got a zero vector for phi\n");
            gretl_matrix_print(IBPh, "IBPh");
        } else {
            gretl_matrix_switch_sign(J.phi);
        }
    }

 bailout:

    gretl_matrix_free(BP);
    gretl_matrix_free(IBP);
    gretl_matrix_free(IBPH);
    gretl_matrix_free(IBPh);
    gretl_matrix_free(e);

    return err;
}

/* initialization of \phi in case the restriction on beta
   is homogeneous */

private int phi_init_homog(Jwrap * J)
{
    GretlMatrixRaw * HH = null;
    GretlMatrixRaw * Hb = null;
    GretlMatrixRaw * b = null;
    int nbeta = J.p1 * J.r;
    int i, err = 0;

    b = gretl_matrix_copy(J.beta);
    HH = gretl_matrix_alloc(J.blen, J.blen);
    Hb = gretl_column_vector_alloc(J.blen);

    if (b == null || HH == null || Hb == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    gretl_matrix_reuse(b, nbeta, 1);

    err = gretl_matrix_multiply_mod(J.H, GretlMatrixMod.transpose,
                                    J.H, GretlMatrixMod.none,
                                    HH, GretlMatrixMod.none);

    if (!err) {
        err = gretl_invert_symmetric_matrix(HH);
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(J.H, GretlMatrixMod.transpose,
                                        b, GretlMatrixMod.none,
                                        Hb, GretlMatrixMod.none);
    }

    if (!err) {
        gretl_matrix_reuse(b, Hb.rows, 1);
        err = gretl_matrix_multiply(HH, Hb, b);
    }

    if (!err) {
        for (i=0; i<b.rows; i++) {
            J.phi.val[i] = b.val[i];
        }
    }

 bailout:

    gretl_matrix_free(HH);
    gretl_matrix_free(Hb);
    gretl_matrix_free(b);

    return err;
}

private int phi_from_beta(Jwrap * J)
{
    int err = 0;

    if (J.H == null) {
        /* just vectorize beta into phi */
        vec_simple(J.phi, J.beta);
    } else if (gretl_is_zero_matrix(J.h)) {
        err = phi_init_homog(J);
    } else {
        /* solve for \phi a la Boswijk */
        err = phi_init_nonhomog(J);
    }

    return err;
}

/* solve the unrestricted eigenvalue problem for \beta */

private int case0(Jwrap * J)
{
    GretlMatrixRaw * B = null;
    int err = 0;

    B = gretl_matrix_alloc(J.p1, J.p);
    if (B == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = gretl_matrix_SVD_johansen_solve(J.R0, J.R1,
                                              null, B, null,
                                              J.r);
    }

    if (!err) {
        gretl_matrix_copy_values(J.beta, B);
    }

    gretl_matrix_free(B);

    return err;
}

/*
   Initialize beta, or the free parameters associated with
   beta in case beta is restricted.
*/

private int beta_init(Jwrap * J)
{
    int err;

    err = case0(J);

    if (!err) {
        err = phi_from_beta(J);
    }

    if (!err && J.H != null) {
        beta_from_phi(J);
    }

    return err;
}

/* do initial computation of alpha based on beta (unless this is
   already done); then, if relevant, impose alpha restriction and make
   psi.
*/

private int alpha_init(Jwrap * J)
{
    int err = 0;

    err = J_compute_alpha(J);
    if (err) {
        return err;
    }

    if (!err) {
        err = psi_from_alpha(J);
    }

    return err;
}

private int make_beta_se(Jwrap * J)
{
    double x;
    int i, j;

    J.bse = gretl_matrix_alloc(J.p1, J.r);
    if (J.bse == null) {
        return GretlError.alloc;
    }

    for (i=0; i<J.Vb.rows; i++) {
        for (j=0; j<J.Vb.cols; j++) {
            x = gretl_matrix_get(J.Vb, i, j);
            if (fabs(x) < 1.0e-20) {
                gretl_matrix_set(J.Vb, i, j, 0);
                x = 0;
            }
            if (j == i) {
                J.bse.val[i] = sqrt(x);
            }
        }
    }

    return 0;
}

private int make_alpha_se(Jwrap * J)
{
    double x;
    int i, j, k;

    J.ase = gretl_matrix_alloc(J.p, J.r);
    if (J.ase == null) {
        return GretlError.alloc;
    }

    k = 0;
    for (i=0; i<J.ase.rows; i++) {
        for (j=0; j<J.ase.cols; j++) {
            x = gretl_matrix_get(J.Va, k, k);
            gretl_matrix_set(J.ase, i, j, sqrt(x));
            k++;
        }
    }

    return 0;
}

/* called by BFGS/simann callback (only) */

private int real_compute_ll(Jwrap * J)
{
    int err = 0;

    if (J.G == null) {
        err = J_compute_alpha(J);
    } else {
        alpha_from_psi(J);
    }

    if (!err) {
        err = make_Omega(J, OmegaCode.only);
    }

    if (!err) {
        gretl_matrix_copy_values(J.Tmppp, J.Omega);
        J.ll = gretl_matrix_log_determinant(J.Tmppp, &err);
    }

    if (!err) {
        J.ll *= -J.T * 0.5;
        J.ll -= J.llk;
    }

    return err;
}

private void sync_with_theta(Jwrap * J, const(double) * theta)
{
    int i, k = 0;

    if (J.H != null) {
        for (i=0; i<J.blen; i++) {
            J.phi.val[i] = theta[k++];
        }
    }

    if (J.G != null) {
        for (i=0; i<J.alen; i++) {
            J.psi.val[i] = theta[k++];
        }
    }
}

/* BFGS (and simann) callback function */

private extern (C) double Jloglik(const(double) * theta, void * data)
{
    Jwrap * J = cast(Jwrap *) data;

    sync_with_theta(J, theta);

    if (J.blen > 0) {
        beta_from_phi(J);
    }

    real_compute_ll(J);

    return J.ll;
}

struct GradHelper {
    GretlMatrixRaw * vPiT;
    GretlMatrixRaw * dPi;
    GretlMatrixRaw * BS;
    GretlMatrixRaw * K;
    GretlMatrixRaw * LK;
    GretlMatrixRaw * psigrad;
    GretlMatrixRaw * phigrad;
}

private void gradhelper_free(GradHelper * g)
{
    if (g == null) {
        return;
    }

    gretl_matrix_free(g.vPiT);
    gretl_matrix_free(g.dPi);
    gretl_matrix_free(g.BS);
    gretl_matrix_free(g.K);
    gretl_matrix_free(g.LK);
    gretl_matrix_free(g.psigrad);
    gretl_matrix_free(g.phigrad);

    free(g);
}

private int add_gradhelper(Jwrap * J)
{
    GradHelper * g = cast(GradHelper *) malloc(GradHelper.sizeof);
    int LKr, err;

    if (g == null) {
        return GretlError.alloc;
    }

    err = make_lsPi(J);
    if (err) {
        free(g);
        return err;
    }

    J.iOmega = gretl_matrix_alloc(J.p, J.p);
    if (J.iOmega == null) {
        free(g);
        return GretlError.alloc;
    }

    g.vPiT = null;
    g.dPi = null;
    g.K = null;
    g.LK = null;
    g.BS = null;
    g.psigrad = null;
    g.phigrad = null;

    if (J.blen == 0) {
        LKr = J.alen;
    } else if (J.G == null) {
        LKr = J.blen;
    } else if (J.H == null) {
        LKr = J.alen;
    } else {
        LKr = (J.blen > J.alen)? J.blen : J.alen;
    }

    clear_gretl_matrix_err();

    g.vPiT = gretl_matrix_alloc(J.p * J.p1, 1);
    g.dPi = gretl_matrix_alloc(J.p * J.p1, 1);
    g.K = gretl_matrix_alloc(J.p1 * J.r, J.p * J.p1);
    g.LK = gretl_matrix_alloc(LKr, J.p * J.p1);

    if (J.G != null) {
        g.BS = gretl_matrix_alloc(J.r, J.p1);
        g.psigrad = gretl_matrix_alloc(J.alen, 1);
    }

    if (J.blen > 0) {
        g.phigrad = gretl_matrix_alloc(J.blen, 1);
    }

    err = get_gretl_matrix_err();
    if (err) {
        gradhelper_free(g);
    } else {
        J.ghelper = g;
    }

    return err;
}

/* optional BFGS callback */

private extern (C) int Jgradient(double * theta, double * gr, int n,
                                 BFGS_CRIT_FUNC f, void * data)
{
    Jwrap * J = cast(Jwrap *) data;
    GradHelper * g = J.ghelper;
    int i, k, err = 0;

    sync_with_theta(J, theta);

    if (J.blen > 0) {
        beta_from_phi(J);
    }

    if (J.G == null) {
        err = J_compute_alpha(J);
    } else {
        alpha_from_psi(J);
    }

    if (!err) {
        err = make_Omega(J, OmegaCode.plus);
    }

    /* form vec(\Pi_{LS}' - \Pi') */

    vec_transpose(g.vPiT, J.Pi);
    gretl_matrix_copy_values(g.dPi, J.lsPi);
    gretl_matrix_subtract_from(g.dPi, g.vPiT);

    /* form psi component of gradient, unless alpha is unrestricted */

    if (J.G != null) {
        gretl_matrix_multiply_mod(J.beta, GretlMatrixMod.transpose,
                                  J.S11, GretlMatrixMod.none,
                                  g.BS, GretlMatrixMod.none);
        gretl_matrix_reuse(g.K, J.p * J.r, -1);
        gretl_matrix_kronecker_product(J.iOmega, g.BS, g.K);

        gretl_matrix_reuse(g.LK, J.alen, -1);
        gretl_matrix_multiply_mod(J.G, GretlMatrixMod.transpose,
                                  g.K, GretlMatrixMod.none,
                                  g.LK, GretlMatrixMod.none);
        gretl_matrix_multiply(g.LK, g.dPi, g.psigrad);
        gretl_matrix_multiply_by_scalar(g.psigrad, J.T);
    }

    /* form phi component, unless beta is fixed */

    if (J.blen > 0) {
        gretl_matrix_multiply_mod(J.alpha, GretlMatrixMod.transpose,
                                  J.iOmega, GretlMatrixMod.none,
                                  J.Tmprp, GretlMatrixMod.none);
        gretl_matrix_reuse(g.K, J.p1 * J.r, -1);
        gretl_matrix_kronecker_product(J.Tmprp, J.S11, g.K);

        if (J.H != null) {
            gretl_matrix_reuse(g.LK, J.blen, -1);
            gretl_matrix_multiply_mod(J.H, GretlMatrixMod.transpose,
                                      g.K, GretlMatrixMod.none,
                                      g.LK, GretlMatrixMod.none);
            gretl_matrix_multiply(g.LK, g.dPi, g.phigrad);
        } else {
            gretl_matrix_multiply(g.K, g.dPi, g.phigrad);
        }

        gretl_matrix_multiply_by_scalar(g.phigrad, J.T);
    }

    /* transcribe results to the gr array */
    k = 0;
    for (i=0; i<J.blen; i++) {
        gr[k++] = g.phigrad.val[i];
    }
    if (J.G != null) {
        for (i=0; i<J.alen; i++) {
            gr[k++] = g.psigrad.val[i];
        }
    }

    return err;
}

private int max_beta_vname(char * tmp,
                           GretlVar * v,
                           const(DatasetRaw) * dset)
{
    int r = gretl_matrix_rows(v.jinfo.Beta);
    const(char) * s;
    int i, ni, n = 0;

    for (i=0; i<r; i++) {
        s = vecm_beta_varname(tmp, v, dset, i);
        ni = cast(int) strlen(s);
        if (ni > n) {
            n = ni;
        }
    }

    return n;
}

enum VECM_WIDTH = 13;

private int printres(Jwrap * J, GretlVar * jvar,
                     GretlRestriction * rset,
                     const(DatasetRaw) * dset,
                     Prn * prn)
{
    const(GretlMatrixRaw) * c = J.beta;
    const(GretlMatrixRaw) * sd = J.bse;
    char[16] s;
    char[SHORTLEN] vname;
    char[16] namefmt;
    int nwid, sdshow;
    int i, j;

    if (J.df > 0) {
        pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g\n").ptr, jvar.ll);
        pprintf(prn, _("Restricted loglikelihood (lr) = %.8g\n").ptr, J.ll);
    } else {
        char[32] ll0, ll1;

        pprintf(prn, "%s = %.8g\n", _("loglikelihood").ptr, J.ll);
        sprintf(ll0.ptr, "%.8g", jvar.ll);
        sprintf(ll1.ptr, "%.8g", J.ll);
        if (strcmp(ll0.ptr, ll1.ptr)) {
            pprintf(prn, "*** warning: should equal %.8g\n", jvar.ll);
        }
    }

    if (J.df > 0) {
        double x = 2.0 * (jvar.ll - J.ll);
        double pv = chisq_cdf_comp(J.df, x);

        pprintf(prn, "2 * (lu - lr) = %g\n", x);
        pprintf(prn, "P(%s(%d) > %g) = %g\n", _("Chi-square").ptr, J.df, x, pv);
        rset_add_results(rset, x, pv, J.ll);
    }

    sdshow = (sd != null && !gretl_is_zero_matrix(sd));

    pputc(prn, '\n');
    pputs(prn, _("Cointegrating vectors").ptr);
    if (sdshow) {
        pprintf(prn, " (%s)", _("standard errors in parentheses").ptr);
    }
    pputs(prn, "\n\n");

    nwid = max_beta_vname(vname.ptr, jvar, dset) + 1;
    sprintf(namefmt.ptr, "%%-%ds", nwid);

    for (i=0; i<J.p1; i++) {
        vecm_beta_varname(vname.ptr, jvar, dset, i);
        pprintf(prn, namefmt.ptr, vname.ptr);

        for (j=0; j<J.r; j++) {
            pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
        }
        pputc(prn, '\n');

        if (sdshow) {
            bufspace(nwid + 1, prn);
            for (j=0; j<J.r; j++) {
                sprintf(s.ptr, "(%#.5g)", gretl_matrix_get(sd, i, j));
                pprintf(prn, "%12s ", s.ptr);
            }
            pputc(prn, '\n');
        }
    }

    c = J.alpha;
    sd = J.ase;

    sdshow = (sd != null && !gretl_is_zero_matrix(sd));

    pputc(prn, '\n');
    pputs(prn, _("alpha (adjustment vectors)").ptr);
    if (sdshow) {
        pprintf(prn, " (%s)", _("standard errors in parentheses").ptr);
    }
    pputs(prn, "\n\n");

    for (i=0; i<J.p; i++) {
        maybe_trim_varname(vname.ptr, dset.varname[jvar.ylist[i+1]]);
        pprintf(prn, namefmt.ptr, vname.ptr);
        for (j=0; j<J.r; j++) {
            pprintf(prn, "%#12.5g ", gretl_matrix_get(c, i, j));
        }
        pputc(prn, '\n');

        if (sdshow) {
            bufspace(nwid + 1, prn);
            for (j=0; j<J.r; j++) {
                sprintf(s.ptr, "(%#.5g)", gretl_matrix_get(sd, i, j));
                pprintf(prn, "%12s ", s.ptr);
            }
            pputc(prn, '\n');
        }
    }

    pputc(prn, '\n');

    return 0;
}

/* solve for unrestricted alpha conditional on beta */

private int J_compute_alpha(Jwrap * J)
{
    GretlMatrixRaw * S01b;
    int err = 0;

    S01b = gretl_matrix_reuse(J.Tmprp, J.p, J.r);
    gretl_matrix_multiply(J.S01, J.beta, S01b);

    if (!err) {
        err = gretl_matrix_qform(J.beta, GretlMatrixMod.transpose,
                                 J.S11, J.qf1, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(J.qf1);
        if (err) {
            gretl_matrix_print(J.qf1, "J->qf1: couldn't invert");
        }
    }

    if (!err) {
        gretl_matrix_multiply(S01b, J.qf1, J.alpha);
    }

    gretl_matrix_reuse(J.Tmprp, J.r, J.p);

    return err;
}

private int allocate_phi(Jwrap * J)
{
    int n = (J.H != null)? J.H.cols : J.p1 * J.r;

    J.phi = gretl_zero_matrix_new(n, 1);

    return (J.phi == null)? GretlError.alloc : 0;
}

private int allocate_psi(Jwrap * J)
{
    int n = (J.G != null)? J.G.cols : J.p * J.r;

    J.psi = gretl_zero_matrix_new(n, 1);

    return (J.psi == null)? GretlError.alloc : 0;
}

/* Needed only when using BFGS (or when initializing using simulated
   annealing): concatenate phi and psi to form the full vector of free
   parameters, theta.  But note: if alpha is unrestricted, we'll just
   compute alpha conditional on beta, and so we don't append phi to
   the free params in theta.
*/

private int make_theta(Jwrap * J)
{
    int i, k = 0;

    if (J.H != null) {
        k += J.blen;
    }

    if (J.G != null) {
        k += J.alen;
    }

    if (k == 0) {
        fprintf(stderr, "jrestrict: make_theta: k = 0\n");
        return GretlError.data;
    }

    J.theta = gretl_column_vector_alloc(k);
    if (J.theta == null) {
        return GretlError.alloc;
    }

    k = 0;

    if (J.H != null) {
        for (i=0; i<J.blen; i++) {
            J.theta.val[k++] = J.phi.val[i];
        }
    }

    if (J.G != null) {
        for (i=0; i<J.alen; i++) {
            J.theta.val[k++] = J.psi.val[i];
        }
    }

    return 0;
}

/* write info from the temporary Jwrap structure into
   the "permanent" GRETL_VAR structure */

private void transcribe_to_jvar(Jwrap * J, GretlVar * jvar,
                                Gretlopt opt)
{
    if (!(opt & Gretlopt.OPT_B)) {
        jvar.jinfo.ll0 = jvar.ll;
        jvar.ll = J.ll;
        jvar.jinfo.lrdf += J.df; /* ? */
    }

    gretl_matrix_replace(&jvar.S, J.Omega);
    J.Omega = null;

    gretl_matrix_replace(&jvar.jinfo.Beta, J.beta);
    J.beta = null;

    gretl_matrix_replace(&jvar.jinfo.Alpha, J.alpha);
    J.alpha = null;

    if (opt & Gretlopt.OPT_B) {
        /* bootstrapping: we don't need the variances */
        return;
    }

    gretl_matrix_replace(&jvar.jinfo.Bvar, J.Vb);
    J.Vb = null;

    gretl_matrix_replace(&jvar.jinfo.Bse, J.bse);
    J.bse = null;

    gretl_matrix_replace(&jvar.jinfo.Ase, J.ase);
    J.ase = null;
}

private int do_bfgs(Jwrap * J, Gretlopt opt, Prn * prn)
{
    int err = add_gradhelper(J);

    if (!err) {
        int maxit = 4000, fncount = 0, grcount = 0;
        int nn = J.theta.rows;
        double toler = 1.0e-11;

        pputs(prn, "LBFGS: using analytical derivatives\n\n");

        err = BFGS_max(J.theta.val, nn, maxit, toler,
                       &fncount, &grcount, &Jloglik, CriterionType.loglik,
                       &Jgradient, J, null, opt | Gretlopt.OPT_L, prn);
    }

    return err;
}

/*
   OPT_L: use LBFGS approach instead of switching algorithm.
   OPT_F: doing full estimation of restricted system, not just
          testing the restriction.
   OPT_J: ("jitter") use simulated annealing in initialization.
   OPT_B: bootstrapping impulse response
*/

int general_vecm_analysis(GretlVar * jvar,
                          GretlRestriction * rset,
                          const(DatasetRaw) * dset,
                          Prn * prn)
{
    Jwrap * J = null;
    Gretlopt opt = gretl_restriction_get_options(rset);
    int full = (opt & Gretlopt.OPT_F) != 0;
    int do_simann = (opt & Gretlopt.OPT_J) != 0;
    int err = 0;

    J = jwrap_new(jvar, opt, &err);
    if (err) {
        return err;
    }

    /* Identification analysis; notation follows quite closely
       Boswijk and Doornik (2004) */

    /* do we have constraints on beta? If so, set up the H matrix */

    if (rset_VECM_bcols(rset) > 0) {
        err = set_up_H(J, rset);
    }

    if (!err) {
        err = allocate_phi(J);
    }

    /* do we have constraints on alpha? if so, set up the G matrix */

    if (!err && rset_VECM_acols(rset) > 0) {
        err = set_up_G(J, rset);
    }

    if (!err) {
        err = allocate_psi(J);
    }

    if (!err) {
        /* check identification via the rank of the Jacobian and
           assign to J->df the value

           (rows_alpha + rows_beta - rank)*rank - Jacobian_rank

           if this is zero, there's no need for numerical optimization.
        */
        err = vecm_id_check(J, jvar, opt, prn);
        if (J.df == 0) {
            J.flags &= cast(ubyte) ~JFlags.use_lbfgs;
            do_simann = 0;
        }
    }

    if (!err && J.blen == 0 && J.G == null) {
        /* beta doesn't need to be estimated */
        J.df = (J.p1 - J.r) * J.r;
        err = real_compute_ll(J);
        goto skipest;
    }

    if (!err && !(opt & Gretlopt.OPT_N)) {
        err = maybe_remove_col_scaling(J, rset);
    }

    if (!err) {
        err = beta_init(J);
    }

    if (!err && (!using_bfgs(J) || J.G != null)) {
        err = alpha_init(J);
    }

    if (!err && (using_bfgs(J) || do_simann)) {
        err = make_theta(J);
    }

    if (!err && do_simann) {
        int n = gretl_vector_get_length(J.theta);

        err = gretl_simann(J.theta.val, n, 4096,
                           &Jloglik, J, opt, prn);
        if (!err) {
            sync_with_theta(J, J.theta.val);
        }
    }

    if (!err) {
        if (using_bfgs(J)) {
            err = do_bfgs(J, opt, prn);
        } else {
            err = switchit(J, opt, prn);
        }
    }

    if (!err && do_scaling(J)) {
        err = replace_col_scaling(J);
    }

 skipest:

    if (!err && !(opt & Gretlopt.OPT_B)) {
        err = variance_from_info_matrix(J, jvar);
    }

    if (!err) {
        if (full) {
            transcribe_to_jvar(J, jvar, Gretlopt.OPT_F);
        } else if (opt & Gretlopt.OPT_B) {
            transcribe_to_jvar(J, jvar, Gretlopt.OPT_B);
        } else {
            printres(J, jvar, rset, dset, prn);
        }
    }

    jwrap_destroy(J);

    return err;
}
