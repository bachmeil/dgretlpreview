/* Port of lib/src/irfboot.c: bootstrap confidence intervals for VAR and VECM impulse responses (residual
   resampling, re-estimation, percentile and studentized intervals). */

module gretl.irfboot;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.gretl_restricth : GretlRestriction;
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;

import gretl.dataset : copy_dataset_obs_info, create_auxiliary_dataset, destroy_dataset;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_error_clear;
import gretl.gretl_list : gretl_list_copy;
import gretl.gretl_utils : gretl_compare_doubles;
import gretl.libset : libset_get_int;
import gretl.libseth : SetKey;
import gretl.matrix_extra : gretl_matrix_row_to_array;
import gretl.plugins : get_plugin_function;
import gretl.var : VAR_write_A_matrix, decompanionize, gretl_VAR_clear, gretl_VAR_do_error_decomp, gretl_VAR_free, johansen_stage_1, n_restricted_terms, reorder_responses;
import gretl.gretl_array : gretl_array_destroy, gretl_array_get_data, gretl_matrix_array_sized;
import gretl.random : gretl_rand_int_max;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


struct IrfBoot {
    int ncoeff; /* number of coefficients per equation */
    int horizon; /* horizon for impulse responses */
    int iters; /* number of iterations */
    int nresp; /* number of responses to calculate */
    GretlMatrixBlock * MB; /* wrapper for some of the following */
    GretlMatrixRaw * rE; /* matrix of resampled original residuals */
    GretlMatrixRaw * Xt; /* row t of X matrix */
    GretlMatrixRaw * Yt; /* yhat at t */
    GretlMatrixRaw * Et; /* residuals at t */
    GretlMatrixRaw * rtmp; /* temporary storage for responses */
    GretlMatrixRaw * wk1; /* workspace */
    GretlMatrixRaw * wk2; /* workspace */
    GretlMatrixRaw * resp; /* impulse response matrix */
    GretlArray * aresp; /* array variant of @resp */
    GretlMatrixRaw * C0; /* initial coefficient estimates (VECM only) */
    int * sample; /* resampling array */
    DatasetRaw * dset; /* dummy dataset for levels (VECM only) */
};

void irf_boot_free(IrfBoot * b)
{
    if (b == null) {
    return;
    }

    gretl_matrix_block_destroy(b.MB);

    gretl_matrix_free(b.Xt);
    gretl_matrix_free(b.Yt);
    gretl_matrix_free(b.Et);
    gretl_matrix_free(b.C0);

    if (b.dset != null) {
    destroy_dataset(b.dset);
    }

    if (b.aresp != null) {
    gretl_array_destroy(b.aresp);
    }

    free(b.sample);
    free(b);
}

int boot_allocate(IrfBoot * b, GretlVar * v)
{
    int n = v.neqns;
    int np = n * levels_order(v);
    int err = 0;

    if (b.nresp > 1) {
    b.aresp = gretl_matrix_array_sized(b.nresp,
                        b.horizon, b.iters,
                        &err);
    if (!err) {
        b.MB = gretl_matrix_block_new(&b.rtmp, np, n,
                       &b.wk1,   n, n,
                       &b.wk2,  np, n,
                       &b.rE, v.T, n,
                       null);
    }
    } else {
    b.MB = gretl_matrix_block_new(&b.rtmp, np, n,
                       &b.wk1,   n, n,
                       &b.wk2,  np, n,
                       &b.rE, v.T, n,
                       &b.resp, b.horizon, b.iters,
                       null);
    }

    if (b.MB == null) {
    return GretlError.alloc;
    }

    b.sample = cast(typeof(b.sample)) malloc(v.T * (*b.sample).sizeof);
    if (b.sample == null) {
    return GretlError.alloc;
    }

    if (v.ci == GretlCmdIndex.var) {
    /* not needed for VECM */
    b.Xt = gretl_matrix_alloc(1, v.X.cols);
    b.Yt = gretl_matrix_alloc(1, v.neqns);
    b.Et = gretl_matrix_alloc(1, v.neqns);

    if (b.Xt == null || b.Yt == null || b.Et == null) {
        return GretlError.alloc;
    }
    }

    return 0;
}

IrfBoot * irf_boot_new(GretlVar * var, int periods, int nresp)
{
    IrfBoot *b;
    int err = 0;

    b = cast(typeof(b)) malloc((*b).sizeof);
    if (b == null) {
    return null;
    }

    b.MB = null;
    b.Xt = null;
    b.Yt = null;
    b.Et = null;
    b.C0 = null;
    b.sample = null;
    b.dset = null;

    b.nresp = nresp;
    b.resp = null;
    b.aresp = null;

    b.horizon = periods;
    b.iters = libset_get_int(SetKey.boot_iters);

    if (var.jinfo != null) {
    b.ncoeff = var.ncoeff + (var.neqns - jrank(var)) +
        n_restricted_terms(var); /* FIXME? */
    } else {
    b.ncoeff = var.ncoeff;
    }

    err = boot_allocate(b, var);

    if (err) {
    irf_boot_free(b);
    b = null;
    }

    return b;
}

/* Copy from @src into @rmax rows in @targ, starting
   at the row given by @offset; it is assumed that
   the two matrices have the same number of columns
   although @targ has more rows.
*/

void copy_part(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int offset, int rmax)
{
    double * tval = targ.val + offset;
    double * sval = cast(double *) src.val;
    size_t sz = rmax * (*tval).sizeof;
    int j;

    for (j=0; j<targ.cols; j++) {
    memcpy(tval, sval, sz);
    tval += targ.rows;
    sval += src.rows;
    }
}

int recalculate_impulse_responses(IrfBoot * b, GretlVar * var, int targ, int shock, int iter)
{
    GretlMatrixRaw * C = var.C;
    GretlMatrixRaw * resp = b.resp;
    double rij;
    int t, i, j, k;
    int n, nxr;
    int err = 0;

    if (var.ord != null) {
    C = reorder_responses(var, &err);
    if (err) {
        return err;
    }
    }

    n = C.cols;
    /* note A is transposed here */
    nxr = var.A.rows - n;
    gretl_matrix_zero(b.rtmp);

    for (t=0; t<b.horizon; t++) {
    if (t == 0) {
        /* initial estimated responses */
        copy_part(b.rtmp, C, 0, n);
    } else {
        /* calculate further estimated responses */
        gretl_matrix_multiply_mod(var.A, GretlMatrixMod.transpose,
                      b.rtmp, GretlMatrixMod.none,
                      b.wk1, GretlMatrixMod.none);
        copy_part(b.wk2, b.wk1, 0, n);
        copy_part(b.wk2, b.rtmp, n, nxr);
        gretl_matrix_copy_values(b.rtmp, b.wk2);
    }
    if (targ >= 0 && shock >= 0) {
        rij = gretl_matrix_get(b.rtmp, targ, shock);
        gretl_matrix_set(resp, t, iter, rij);
    } else if (shock >= 0) {
        /* all targets for one shock */
        for (i=0; i<n; i++) {
        resp = cast(GretlMatrixRaw *) gretl_array_get_data(b.aresp, i);
        rij = gretl_matrix_get(b.rtmp, i, shock);
        gretl_matrix_set(resp, t, iter, rij);
        }
    } else if (targ >= 0) {
        /* all shocks for one target */
        for (j=0; j<n; j++) {
        resp = cast(GretlMatrixRaw *) gretl_array_get_data(b.aresp, j);
        rij = gretl_matrix_get(b.rtmp, targ, j);
        gretl_matrix_set(resp, t, iter, rij);
        }
    } else {
        /* all shocks, all targets */
        k = 0;
        for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            resp = cast(GretlMatrixRaw *) gretl_array_get_data(b.aresp, k++);
            rij = gretl_matrix_get(b.rtmp, i, j);
            gretl_matrix_set(resp, t, iter, rij);
        }
        }
    }
    }

    if (C != var.C) {
    gretl_matrix_free(C);
    }

    return err;
}

void maybe_resize_vecm_matrices(GretlVar * v)
{
    int nc0 = (v.order * v.neqns) + v.ifc + v.jinfo.seasonals;

    if (v.xlist != null) {
    nc0 += v.xlist[0];
    }

    if (v.detflags & DetFlags.trend) {
    nc0++;
    }

    if (v.X != null && v.X.cols > nc0) {
    /* for stage 1: skip the extra EC-term columns in X */
    gretl_matrix_reuse(v.X, -1, nc0);
    gretl_matrix_reuse(v.B, nc0, -1);
    }
}

/* In re-estimation of VAR or VECM we'll tolerate only so many cases
   of near-perfect collinearity. These can surely arise by chance but
   if their frequency is too high this would seem to indicate an
   inherent problem.
*/

enum MAXSING = 0.10;

int irf_fatal(int err, IrfBoot * b, int iter, int scount)
{
    if (err != GretlError.singular || iter == 0) {
    return 1;
    } else if ((scount / cast(double) b.iters) > MAXSING) {
    return 1;
    } else {
    return 0;
    }
}

int re_estimate_VECM(IrfBoot * b, GretlVar * v, int targ, int shock, int iter, int scount)
{
    alias JbrFn = extern (C) int function(GretlVar *, const(DatasetRaw) *);
    __gshared JbrFn jbr = null;
    int err = 0;

    gretl_error_clear();

    if (jbr == null) {
    /* open the Johansen plugin */
    jbr = cast(JbrFn) get_plugin_function("johansen_boot_round");
    if (jbr == null) {
        return GretlError.fopen;
    }
    }

    /* The various VECM matrices may need to be re-set to the sizes
       expected by johansen_stage_1() */
    maybe_resize_vecm_matrices(v);

    err = johansen_stage_1(v, b.dset, Gretlopt.OPT_NONE, null);

    if (!err) {
    /* call the plugin function */
    err = jbr(v, b.dset);
    }

    if (!err) {
    err = gretl_VAR_do_error_decomp(v.S, v.C, v.ord);
    }

    if (!err) {
    err = recalculate_impulse_responses(b, v, targ, shock, iter);
    }

    return err;
}

int re_estimate_VAR(IrfBoot * b, GretlVar * v, int targ, int shock, int iter)
{
    int err;

    err = gretl_matrix_multi_ols(v.Y, v.X, v.B, v.E, null);

    if (!err) {
    VAR_write_A_matrix(v, GretlMatrixMod.transpose);
    }

    if (!err) {
    gretl_matrix_multiply_mod(v.E, GretlMatrixMod.transpose,
                  v.E, GretlMatrixMod.none,
                  v.S, GretlMatrixMod.none);
    gretl_matrix_divide_by_scalar(v.S, v.T);
    err = gretl_VAR_do_error_decomp(v.S, v.C, v.ord);
    }

    if (!err) {
    err = recalculate_impulse_responses(b, v, targ, shock, iter);
    }

    return err;
}

/* VECM: make a matrix, rbeta, containing the coefficients on the
   "restricted" constant or trend and/or restricted exogenous
   variables, if this is required.
*/

GretlMatrixRaw * make_restricted_coeff_matrix(GretlVar * v)
{
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * rbeta = null;
    int rank = v.jinfo.rank;
    int nr = n_restricted_terms(v);
    double x;
    int j, k, err = 0;

    b = gretl_matrix_alloc(rank, nr);
    if (b == null) {
    return null;
    }

    for (j=0; j<rank; j++) {
    /* extract the last row(s) of \beta, transposed */
    for (k=0; k<nr; k++) {
        x = gretl_matrix_get(v.jinfo.Beta, v.neqns + k, j);
        gretl_matrix_set(b, j, k, x);
    }
    }

    rbeta = gretl_matrix_multiply_new(v.jinfo.Alpha, b, &err);

    gretl_matrix_free(b);

    return rbeta;
}

/* VECM: make a record of all the original coefficients in the VAR
   representation of the VECM, so we have these on hand in convenient
   form for recomputing the dataset after resampling the residuals,
   or for forecasting.

   We get these coefficients from three sources: (1) the
   companion matrix, A; (2) the VECM models (unrestricted constant
   and/or trend, seasonals if applicable); and (3) the implied
   coefficient on a restricted constant or trend.

   Note that the construction here depends on the order in which
   variables are added to the dataset; that order can't be changed
   without breaking stuff here.
*/

GretlMatrixRaw * VAR_coeff_matrix_from_VECM(GretlVar * var, int Atrans)
{
    GretlMatrixRaw * C0 = null;
    GretlMatrixRaw * rbeta = null;
    GretlMatrixRaw * B = var.B;
    int order = var.order + 1;
    int nexo = (var.xlist != null)? var.xlist[0] : 0;
    int utrend = jcode(var) == JohansenCode.unrest_trend;
    int ndelta = var.order * var.neqns;
    int nseas = var.jinfo.seasonals;
    int nr = n_restricted_terms(var);
    int ncoeff, nb;
    int X0, S0, T0;
    double aij;
    int i, j, k;

    /* total coeffs in VAR representation */
    ncoeff = var.ncoeff + (var.neqns - var.jinfo.rank) + nr;

    /* Is var->B required? And if so how many B rows do we need? */
    nb = var.ifc + nexo + nseas + utrend;

    if (nb > 0) {
    if (B == null) {
        gretl_errmsg_set("VAR coefficient matrix B is missing!");
        return null;
    } else if (B.rows != nb && B.rows != var.ncoeff) {
        gretl_errmsg_set("VAR coefficient matrix B is of wrong size!");
        fprintf(stderr, "B should have %d or %d rows, but has %d\n".ptr,
            nb, var.ncoeff, B.rows);
        gretl_matrix_print(B, "var->B");
        return null;
    } else {
        nb = B.rows;
    }
    }

    if (nr > 0) {
    rbeta = make_restricted_coeff_matrix(var);
    if (rbeta == null) {
        return null;
    }
    }

    C0 = gretl_matrix_alloc(var.neqns, ncoeff);
    if (C0 == null) {
    gretl_matrix_free(rbeta);
    return null;
    }

    /* position of first exog var coeff in VECM model */
    X0 = var.ifc + ndelta;
    /* position of first seasonal coeff in VECM model */
    S0 = X0 + nexo;
    /* position of trend coeff in VECM model */
    T0 = S0 + nseas;

    for (i=0; i<var.neqns; i++) {
    const(double) * coeff = null;
    int col = 0;

    if (nb > 0) {
        coeff = B.val + i * nb;
    }

    /* constant, if present */
    if (var.ifc) {
        gretl_matrix_set(C0, i, col++, coeff[0]);
    }

    /* endogenous vars: use companion matrix */
    for (j=0; j<var.neqns; j++) {
        for (k=0; k<order; k++) {
        if (Atrans) {
            aij = gretl_matrix_get(var.A, k * var.neqns + j, i);
        } else {
            aij = gretl_matrix_get(var.A, i, k * var.neqns + j);
        }
        gretl_matrix_set(C0, i, col++, aij);
        }
    }

    /* exogenous vars */
    for (j=0; j<nexo; j++) {
        gretl_matrix_set(C0, i, col++, coeff[X0+j]);
    }

    /* seasonals, if present */
    for (j=0; j<nseas; j++) {
        gretl_matrix_set(C0, i, col++, coeff[S0+j]);
    }

    /* unrestricted trend, if present */
    if (utrend) {
        gretl_matrix_set(C0, i, col++, coeff[T0]);
    }

    /* additional restricted term(s), if present */
    if (rbeta != null) {
        for (j=0; j<nr; j++) {
        aij = gretl_matrix_get(rbeta, i, j);
        gretl_matrix_set(C0, i, col++, aij);
        }
    }
    }

    if (rbeta != null) {
    gretl_matrix_free(rbeta);
    }

    return C0;
}

/* VECM: copy levels of Y vars from main Z into temporary dataset (for
   passing to Johansen stage 1).  If there are restricted exogenous
   vars they must be included in the temp dataset too.
*/

int init_VECM_dataset(IrfBoot * b, GretlVar * var, const(DatasetRaw) * dset)
{
    int nv = var.neqns + 1;
    int i, j, vi, t;

    if (var.rlist != null) {
    nv += var.rlist[0];
    }

    b.dset = create_auxiliary_dataset(nv, dset.n, Gretlopt.OPT_NONE);
    if (b.dset == null) {
    return GretlError.alloc;
    }

    copy_dataset_obs_info(b.dset, dset);
    b.dset.t1 = dset.t1;
    b.dset.t2 = dset.t2;

    /* copy levels of Y into boot->Z and adjust ylist */
    for (i=0, j=1; i<var.neqns; i++, j++) {
    vi = var.ylist[i+1];
    for (t=0; t<dset.n; t++) {
        b.dset.Z[j][t] = dset.Z[vi][t];
    }
    var.ylist[j] = j;
    }

    if (var.rlist != null) {
    /* copy restricted X into boot->Z and adjust rlist */
    for (i=1; i<=var.rlist[0]; i++, j++) {
        vi = var.rlist[i];
        for (t=0; t<dset.n; t++) {
        b.dset.Z[j][t] = dset.Z[vi][t];
        }
        var.rlist[i] = j;
    }
    }

    return 0;
}

/* Compute VECM Y (in levels) using the VAR representation. */

void compute_VECM_dataset(IrfBoot * b, GretlVar * var, int iter)
{
    int order = var.order + 1;
    int nexo = (var.xlist != null)? var.xlist[0] : 0;
    int nseas = var.jinfo.seasonals;
    double cij, eti, xti;
    int i, j, k, s, t;

    for (t=var.t1, s=0; t<=var.t2; t++, s++) {
    for (i=0; i<var.neqns; i++) {
        double bti = 0.0;
        int xcol = var.ifc + var.neqns * var.order;
        int col = 0;

        /* unrestricted constant, if present */
        if (var.ifc) {
        cij = gretl_matrix_get(b.C0, i, col++);
        bti += cij;
        }

        /* lags of endogenous vars */
        for (j=1; j<=var.neqns; j++) {
        for (k=1; k<=order; k++) {
            cij = gretl_matrix_get(b.C0, i, col++);
            bti += cij * b.dset.Z[j][t-k];
        }
        }

        /* exogenous vars, if present */
        for (j=0; j<nexo; j++) {
        cij = gretl_matrix_get(b.C0, i, col++);
        xti = gretl_matrix_get(var.X, s, xcol++);
        bti += cij * xti;
        }

        /* seasonals, if present */
        for (j=0; j<nseas; j++) {
        cij = gretl_matrix_get(b.C0, i, col++);
        xti = gretl_matrix_get(var.X, s, xcol++);
        bti += cij * xti;
        }

        if (jcode(var) == JohansenCode.unrest_trend) {
        /* unrestricted trend */
        cij = gretl_matrix_get(b.C0, i, col++);
        bti += cij * (t + 1);
        } else if (jcode(var) == JohansenCode.rest_const) {
        /* restricted constant */
        bti += gretl_matrix_get(b.C0, i, col++);
        } else if (jcode(var) == JohansenCode.rest_trend) {
        /* restricted trend */
        cij = gretl_matrix_get(b.C0, i, col++);
        bti += cij * t;
        }

        /* restricted exogenous vars, if present */
        if (var.rlist != null) {
        for (j=1; j<=var.rlist[0]; j++) {
            cij = gretl_matrix_get(b.C0, i, col++);
            bti += cij * b.dset.Z[j+var.neqns][t-1];
        }
        }

        /* set level of Y(t, i) to fitted value + re-sampled error */
        eti = gretl_matrix_get(b.rE, s, i);
        b.dset.Z[i+1][t] = bti + eti;
    }
    }

    if (var.X != null) {
    /* now rewrite lagged differences into X matrix */
    k = var.ifc;
    for (i=1; i<=var.neqns; i++) {
        for (j=1; j<=var.order; j++) {
        s = 0;
        for (t=var.t1; t<=var.t2; t++) {
            xti = b.dset.Z[i][t-j] - b.dset.Z[i][t-j-1];
            gretl_matrix_set(var.X, s++, k, xti);
        }
        k++;
        }
    }
    }
}

void grab_row_from(GretlMatrixRaw * targ, GretlMatrixRaw * src, int t)
{
    int j;

    for (j=0; j<src.cols; j++) {
    targ.val[j] = gretl_matrix_get(src, t, j);
    }
}

void stick_row_into(GretlMatrixRaw * targ, GretlMatrixRaw * src, int t)
{
    int j;

    for (j=0; j<src.cols; j++) {
    gretl_matrix_set(targ, t, j, src.val[j]);
    }
}

/* (Re-)fill the bootstrap dataset with artificial data, based on the
   re-sampled residuals from the original VAR (case of simple VAR, not
   VECM).
*/

void compute_VAR_dataset(IrfBoot * b, GretlVar * var, GretlVar * vbak, int iter)
{
    double x;
    int i, j, k, t;
    int nl = var_n_lags(var);

    for (t=0; t<var.T; t++) {
    /* extract row of var->X at t */
    grab_row_from(b.Xt, var.X, t);

    /* multiply Xt into original coeff matrix, forming Yt */
    gretl_matrix_multiply(b.Xt, vbak.B, b.Yt);

    /* extract resampled residuals at t */
    grab_row_from(b.Et, b.rE, t);

    /* add resampled residual to Yt */
    gretl_matrix_add_to(b.Yt, b.Et);

    /* write into big Y matrix */
    stick_row_into(var.Y, b.Yt, t);

    /* revise lagged Y columns in X */
    for (i=0; i<var.neqns; i++) {
        k = var.ifc + i*nl;
        x = b.Yt.val[i];
        for (j=1; j<=nl && t+j < var.T; j++) {
        gretl_matrix_set(var.X, t+j, k++, x);
        }
    }
    }

}

/* Resample the original VAR or VECM residuals, stored in
   vbak->E, writing the new_ sample into b->rE.

   Note the facility to "resample" _without_ actually changing the
   order, if BDEBUG > 1.  This is useful for checking that the VarPrint.irf
   bootstrap rounds are idempotent: we should then get exactly the
   same set of responses as in the original estimation of the
   VAR/VECM.
*/

void irf_resample_resids(IrfBoot * b, GretlVar * vbak)
{
    double eti;
    int i, t;

    /* construct sampling array */
    for (t=0; t<vbak.T; t++) {
    b.sample[t] = gretl_rand_int_max(vbak.T);
    }

    /* draw from the original residuals */
    for (i=0; i<vbak.neqns; i++) {
    for (t=0; t<vbak.T; t++) {
        eti = gretl_matrix_get(vbak.E, b.sample[t], i);
        gretl_matrix_set(b.rE, t, i, eti);
    }
    }
}

int irf_boot_quantiles(IrfBoot * b, GretlMatrixRaw * R, double alpha, GretlMatrixRaw * point)
{
    double r0;
    double * ri;
    int nr = b.iters + 1;
    int i, ilo, ihi;

    ri = cast(typeof(ri)) malloc(nr * (*ri).sizeof);
    if (ri == null) {
    return GretlError.alloc;
    }

    ilo = cast(int) (nr * alpha / 2.0);
    ihi = cast(int) (nr * (1.0 - alpha / 2.0));

    if (b.nresp > 1) {
    GretlMatrixRaw * resp;
    int j, Rcol = 0;

    for (j=0; j<b.nresp; j++) {
        resp = cast(GretlMatrixRaw *) gretl_array_get_data(b.aresp, j);
        for (i=0; i<b.horizon; i++) {
        r0 = gretl_matrix_get(point, i, j);
        gretl_matrix_row_to_array(resp, i, ri);
        ri[b.iters] = r0;
        qsort(ri, nr, (*ri).sizeof, &gretl_compare_doubles);
        gretl_matrix_set(R, i, Rcol, r0);
        gretl_matrix_set(R, i, Rcol+1, ri[ilo-1]);
        gretl_matrix_set(R, i, Rcol+2, ri[ihi-1]);
        }
        Rcol += 3;
    }
    } else {
    for (i=0; i<b.horizon; i++) {
        r0 = point.val[i];
        gretl_matrix_row_to_array(b.resp, i, ri);
        ri[b.iters] = r0;
        qsort(ri, nr, (*ri).sizeof, &gretl_compare_doubles);
        gretl_matrix_set(R, i, 0, r0);
        gretl_matrix_set(R, i, 1, ri[ilo-1]);
        gretl_matrix_set(R, i, 2, ri[ihi-1]);
    }
    }

    free(ri);

    return 0;
}

void restore_VAR_data(GretlVar * v, GretlVar * vbak)
{
    gretl_matrix_replace(&v.Y, vbak.Y);
    gretl_matrix_replace(&v.X, vbak.X);
    gretl_matrix_replace(&v.B, vbak.B);
    gretl_matrix_replace(&v.XTX, vbak.XTX);
    gretl_matrix_replace(&v.A, vbak.A);
    gretl_matrix_replace(&v.E, vbak.E);
    gretl_matrix_replace(&v.C, vbak.C);
    gretl_matrix_replace(&v.S, vbak.S);

    if (vbak.jinfo != null) {
    gretl_matrix_replace(&v.jinfo.R0, vbak.jinfo.R0);
    gretl_matrix_replace(&v.jinfo.R1, vbak.jinfo.R1);
    gretl_matrix_replace(&v.jinfo.S00, vbak.jinfo.S00);
    gretl_matrix_replace(&v.jinfo.S11, vbak.jinfo.S11);
    gretl_matrix_replace(&v.jinfo.S01, vbak.jinfo.S01);
    gretl_matrix_replace(&v.jinfo.Beta, vbak.jinfo.Beta);
    gretl_matrix_replace(&v.jinfo.Alpha, vbak.jinfo.Alpha);

    if (vbak.ylist != null) {
        free(v.ylist);
        v.ylist = vbak.ylist;
    }

    if (vbak.rlist != null) {
        free(v.rlist);
        v.rlist = vbak.rlist;
    }

    free(vbak.jinfo);
    }

    free(vbak);
}

/* make a back-up copy of the VAR data that we'll want to restore on
   exit from the VarPrint.irf bootstrap operation */

void matrix_swap_copy(GretlMatrixRaw ** targ, GretlMatrixRaw ** src)
{
    *targ = *src;
    *src = gretl_matrix_copy(*targ);
}

int list_swap_copy(int ** targ, int ** src)
{
    *targ = *src;
    *src = gretl_list_copy(*targ);
    return *src == null ? GretlError.alloc : 0;
}

GretlVar * back_up_VAR(GretlVar * v)
{
    GretlVar * vbak;
    int err = 0;

    vbak = cast(typeof(vbak)) malloc((*vbak).sizeof);
    if (vbak == null) {
    return null;
    }

    gretl_VAR_clear(vbak);
    clear_gretl_matrix_err();

    if (v.Y != null) {
    matrix_swap_copy(&vbak.Y, &v.Y);
    }
    if (v.B != null) {
    matrix_swap_copy(&vbak.B, &v.B);
    }

    if (v.X != null) {
    if (v.xcols > v.X.cols) {
        int save_k = v.X.cols;

        vbak.X = v.X;
        gretl_matrix_reuse(vbak.X, -1, v.xcols);
        v.X = gretl_matrix_copy(vbak.X);
        gretl_matrix_reuse(vbak.X, -1, save_k);
        gretl_matrix_reuse(v.X, -1, save_k);
    } else {
        matrix_swap_copy(&vbak.X, &v.X);
    }
    }

    if (v.XTX != null) {
    matrix_swap_copy(&vbak.XTX, &v.XTX);
    }
    matrix_swap_copy(&vbak.E, &v.E);
    matrix_swap_copy(&vbak.C, &v.C);
    matrix_swap_copy(&vbak.S, &v.S);

    /* @A is a special case */
    vbak.A = v.A;
    v.A = decompanionize(vbak.A, v.neqns, GretlMatrixMod.transpose);

    err = get_gretl_matrix_err() || v.A == null;

    if (!err && v.jinfo != null) {
    vbak.jinfo = cast(typeof(vbak.jinfo)) calloc(1, (*vbak.jinfo).sizeof);
    if (vbak.jinfo == null) {
        err = GretlError.alloc;
    } else {
        err = list_swap_copy(&vbak.ylist, &v.ylist);
    }
    }

    if (!err && v.rlist != null) {
    err = list_swap_copy(&vbak.rlist, &v.rlist);
    }

    if (!err && vbak.jinfo != null) {
    matrix_swap_copy(&vbak.jinfo.R0, &v.jinfo.R0);
    matrix_swap_copy(&vbak.jinfo.R1, &v.jinfo.R1);
    matrix_swap_copy(&vbak.jinfo.S00, &v.jinfo.S00);
    matrix_swap_copy(&vbak.jinfo.S11, &v.jinfo.S11);
    matrix_swap_copy(&vbak.jinfo.S01, &v.jinfo.S01);
    matrix_swap_copy(&vbak.jinfo.Beta, &v.jinfo.Beta);
    matrix_swap_copy(&vbak.jinfo.Alpha, &v.jinfo.Alpha);
    err = get_gretl_matrix_err();
    }

    if (err) {
    restore_VAR_data(v, vbak);
    gretl_VAR_free(vbak);
    vbak = null;
    } else {
    vbak.T = v.T;
    vbak.neqns = v.neqns;
    vbak.order = v.order;
    }

    return vbak;
}

/* public bootstrapping function, called from var.c */

GretlMatrixRaw * irf_bootstrap(GretlVar * var, int targ, int shock, int periods, double alpha, GretlMatrixRaw * point, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * R = null; /* the return value */
    GretlVar * vbak = null;
    IrfBoot *boot = null;
    int scount = 0;
    int nresp = 1;
    int iter;

    if (targ < 0 && shock < 0) {
    /* compute all */
    nresp = var.neqns * var.neqns;
    } else if (targ < 0 || shock < 0) {
    nresp = var.neqns;
    }

    R = gretl_zero_matrix_new(periods, 3*nresp);
    if (R == null) {
    *err = GretlError.alloc;
    return null;
    }

    vbak = back_up_VAR(var);
    if (vbak == null) {
    *err = GretlError.alloc;
    gretl_matrix_free(R);
    return null;
    }

    boot = irf_boot_new(var, periods, nresp);
    if (boot == null) {
    *err = GretlError.alloc;
    goto bailout;
    }

    if (var.ci == GretlCmdIndex.vecm) {
    boot.C0 = VAR_coeff_matrix_from_VECM(var, 1);
    if (boot.C0 == null) {
        *err = GretlError.alloc;
    } else {
        *err = init_VECM_dataset(boot, var, dset);
    }
    }

    for (iter=0; iter<boot.iters && !*err; iter++) {
    irf_resample_resids(boot, vbak);
    if (var.ci == GretlCmdIndex.vecm) {
        compute_VECM_dataset(boot, var, iter);
        *err = re_estimate_VECM(boot, var, targ, shock, iter, scount);
    } else {
        compute_VAR_dataset(boot, var, vbak, iter);
        *err = re_estimate_VAR(boot, var, targ, shock, iter);
    }
    if (*err && !irf_fatal(*err, boot, iter, scount)) {
        /* excessive collinearity: try again, unless this is
           becoming a serious habit
        */
        scount++;
        iter--;
        *err = 0;
    }
    }

    if (*err && scount / cast(double) boot.iters >= MAXSING) {
    gretl_errmsg_set("Excessive collinearity in resampled datasets");
    }

    if (!*err) {
    *err = irf_boot_quantiles(boot, R, alpha, point);
    }

    irf_boot_free(boot);

 bailout:

    restore_VAR_data(var, vbak);

    if (*err) {
    gretl_matrix_free(R);
    R = null;
    }

    return R;
}

