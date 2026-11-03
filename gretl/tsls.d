/* Port of lib/src/tsls.c: two-stage least squares (and LIML helpers): instrument and regressor list handling,
   first-stage regressions, second-stage estimation, Sargan/Basmann over-identification and Hausman tests. The
   header tsls.h has only prototypes. */

module gretl.tsls;


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

import gretl.dataset : dataset_add_series, dataset_drop_last_variables;
import gretl.dataseth : dataset_is_time_series;
import gretl.describe : gretl_corr_rsq;
import gretl.estimate : dwstat, lsq, rhohat;
import gretl.genparseh : ConstId;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_sprintf;
import gretl.gretl_list : gretl_list_add, gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_insert_list, gretl_list_new, gretl_list_separator_position, gretl_list_split_on_separator, gretl_list_sublist, gretl_lists_join_with_separator, in_gretl_list, reglist_check_for_const;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : ModelTest, clear_model, gretl_model_detach_data_item, gretl_model_get_data, gretl_model_get_data_full, gretl_model_get_list, gretl_model_init, gretl_model_set_data, gretl_model_set_double, gretl_model_set_int, gretl_model_set_list_as_data, gretl_model_set_matrix_as_data, maybe_add_test_to_model, mle_criteria, model_count_minus, model_test_new, model_test_set_dfn, model_test_set_pvalue, model_test_set_teststat, model_test_set_value;
import gretl.gretl_modelh : ModelTestType;
import gretl.matrix_extra : gretl_matrix_cut_rows_cols, gretl_matrix_data_subset_masked, gretl_matrix_rank_mask;
import gretl.missing : model_has_missing_obs, model_missing;
import gretl.options : get_optval_string, incompatible_options;
import gretl.pvalues : chisq_cdf_comp;
import gretl.qr_estimate : qr_tsls_vcv;
import gretl.compare : wald_omit_F, wald_omit_chisq;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


struct IvInfo {
    Model * pmod; /* pointer to the model to be returned */
    GretlMatrixRaw * Q; /* matrix for creating "hat" variables */
    GretlMatrixRaw * g; /* auxiliary matrix */
    char * missmask; /* missing observations mask */
    int * reglist; /* full list of regressors */
    int * instlist; /* list of instruments */
    int * hatlist; /* list of fitted values of endogenous regressors */
    int * s2list; /* regression list for second stage */
    int * exolist; /* list of exogenous regressors */
    int * endolist; /* list of endogenous regressors */
    int * idroplist; /* list of dropped instruments */
    int overid; /* over-identification rank */
};

IvInfo * iv_info_new(Model * pmod)
{
    IvInfo * ivi = cast(IvInfo *) malloc(IvInfo.sizeof);

    if (ivi != null) {
    ivi.pmod = pmod;
    ivi.Q = ivi.g = null;
    ivi.missmask = null;
    ivi.reglist = null;
    ivi.instlist = null;
    ivi.hatlist = null;
    ivi.s2list = null;
    ivi.exolist = null;
    ivi.endolist = null;
    ivi.idroplist = null;
    ivi.overid = 0;
    }

    return ivi;
}

void iv_info_destroy(IvInfo * ivi)
{
    if (ivi == null) {
    return;
    }

    gretl_matrix_free(ivi.Q);
    gretl_matrix_free(ivi.g);

    free(ivi.missmask);
    free(ivi.reglist);
    free(ivi.instlist);
    free(ivi.hatlist);
    free(ivi.s2list);
    free(ivi.exolist);
    free(ivi.endolist);
    free(ivi.idroplist);

    free(ivi);
}

void tsls_omitzero(int * list, const(char) * mask, const(DatasetRaw) * dset)
{
    int i, v, t, allzero;

    for (i=list[0]; i>1; i--) {
        v = list[i];
        allzero = 1;
        for (t=dset.t1; t<=dset.t2; t++) {
            if (mask != null && mask[t - dset.t1]) {
                continue;
            }
            if (dset.Z[v][t] != 0.0) {
                allzero = 0;
                break;
            }
        }
        if (allzero) {
            gretl_list_delete_at_pos(list, i);
        }
    }
}

/**
 * tsls_free_data:
 * @pmod: model to operate on.
 *
 * Frees the first-stage data attached to a two-stage least squares
 * model that was estimated using %OPT_E (part of a system).
 */

void tsls_free_data(const(Model) * pmod)
{
    GretlMatrixRaw * endog = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "endog");
    double ** X = cast(double **) gretl_model_get_data(pmod, "tslsX");
    int i, m = 0;

    if (endog != null && X != null) {
        for (i=0; i<pmod.ncoeff; i++) {
            if (endog.val[i] != 0) m++;
        }
        for (i=0; i<m; i++) {
            free(X[i]);
        }
    }
}

/* when tsls is called for initial estimation of an equation
   within a system of equations, save the first-stage fitted
   values on the model under the key "tslsX".
*/

int tsls_save_data(IvInfo * ivi, DatasetRaw * dset)
{
    Model * pmod = ivi.pmod;
    GretlMatrixRaw * endog = null;
    double ** X = null;
    size_t Xsize, xs_old;
    int recycle_X = 0;
    int recycle_e = 0;
    int i, v, err = 0;

    Xsize = ivi.hatlist[0] * (*X).sizeof;

    /* re-use old pointers if applicable */

    X = cast(double **) gretl_model_get_data_full(pmod, "tslsX", null, null, &xs_old);
    if (X != null) {
        if (Xsize == xs_old) {
            recycle_X = 1;
        } else {
            tsls_free_data(pmod);
            gretl_model_detach_data_item(pmod, "tslsX");
            free(X);
            X = null;
        }
    }

    if (!recycle_X && Xsize > 0) {
        X = cast(typeof(X)) malloc(Xsize);
        if (X == null) {
            return GretlError.alloc;
        }
    }

    endog = cast(typeof(endog)) gretl_model_get_data(pmod, "endog");
    if (endog != null) {
        if (gretl_vector_get_length(endog) == pmod.ncoeff) {
            recycle_e = 1;
        } else {
            gretl_model_detach_data_item(pmod, "endog");
            gretl_matrix_free(endog);
            endog = null;
        }
    }

    if (!recycle_e) {
        endog = gretl_matrix_alloc(pmod.ncoeff, 1);
        if (endog == null) {
            free(X);
            return GretlError.alloc;
        }
    }

    /* steal the appropriate columns from Z */
    for (i=1; i<=ivi.hatlist[0]; i++) {
        v = ivi.hatlist[i];
        X[i-1] = dset.Z[v];
        dset.Z[v] = null;
    }

    for (i=0; i<pmod.ncoeff; i++) {
        v = pmod.list[i+2];
        endog.val[i] = !in_gretl_list(ivi.exolist, v);
    }

     /* now attach X and endog to the model */
    if (!recycle_X && X != null) {
        gretl_model_set_data(pmod, "tslsX", X, GretlType.double_array,
                             Xsize);
    }
    if (!recycle_e && endog != null) {
        gretl_model_set_matrix_as_data(pmod, "endog", endog);
    }

    return err;
}

void add_tsls_var(int * list, int v, Gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i;

    if (opt & Gretlopt.OPT_T) {
        /* add as instrument */
        list[0] += 1;
        list[list[0]] = v;
    } else if (opt & Gretlopt.OPT_B) {
        /* add as exogenous regressor */
        list[0] += 2;
        for (i=list[0]-1; i>pos; i--) {
            list[i] = list[i-1];
        }
        list[pos] = v;
        list[list[0]] = v;
    } else {
        /* add as endogenous regressor */
        list[0] += 1;
        for (i=list[0]; i>pos; i--) {
            list[i] = list[i-1];
        }
        list[pos] = v;
    }
}

void delete_tsls_var(int * list, int v, Gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, j;

    if (opt & Gretlopt.OPT_T) {
        /* delete as instrument */
        for (i=pos+1; i<=list[0]; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
                break;
            }
        }
    } else if (opt & Gretlopt.OPT_B) {
        /* delete from both sub-lists */
        for (i=2; i<=list[0]; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
            }
        }
    } else {
        /* delete from regressor list */
        for (i=2; i<pos; i++) {
            if (list[i] == v) {
                for (j=i; j<list[0]; j++) {
                    list[j] = list[j + 1];
                }
                list[0] -= 1;
                break;
            }
        }
    }
}

int in_ivreg_list(const(int) * list, int v, Gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, imin, imax;

    if (pos == 0) {
        return 0;
    }

    if (opt & Gretlopt.OPT_T) {
        imin = pos + 1;
        imax = list[0];
    } else if (opt & Gretlopt.OPT_B) {
        imin = 2;
        imax = list[0];
    } else {
        imin = 2;
        imax = pos - 1;
    }

    for (i=imin; i<=imax; i++) {
        if (list[i] == v) {
            return i;
        }
    }

    return 0;
}

/**
 * ivreg_list_omit:
 * @orig: list specifying original IV model.
 * @drop: list of variables to be omitted.
 * @opt: may include %OPT_T (omit variable only as instrument)
 * or %OPT_B (omit both as regressor and instrument).  The default
 * is to omit (only) as regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new IVREG specification list, by first copying @orig
 * then deleting from the copy the variables found in @drop.
 *
 * Returns: the new, reduced list or %NULL on error.
 */

int * ivreg_list_omit(const(int) * orig, const(int) * drop, Gretlopt opt, int * err)
{
    int * newlist;
    int i;

    *err = incompatible_options(opt, Gretlopt.OPT_T | Gretlopt.OPT_B);
    if (*err) {
        return null;
    }

    newlist = gretl_list_copy(orig);

    for (i=1; i<=drop[0]; i++) {
        if (in_ivreg_list(orig, drop[i], opt)) {
            delete_tsls_var(newlist, drop[i], opt);
        } else {
            *err = GretlError.unspec;
        }
    }

    if (*err) {
        free(newlist);
        newlist = null;
    }

    return newlist;
}

/**
 * ivreg_list_add:
 * @orig: list specifying original IV model.
 * @add: list of variables to be added.
 * @opt: may include %OPT_T (add variable only as instrument)
 * or %OPT_B (add both as regressor and instrument).  The default
 * is to add as an endogenous regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new IVREG specification list, by copying @orig
 * and adding the variables found in @add.
 *
 * Returns: the new, augmented list or %NULL on error.
 */

int * ivreg_list_add(const(int) * orig, const(int) * add, Gretlopt opt, int * err)
{
    int norig = orig[0];
    int nadd = add[0];
    int * newlist;
    int i;

    *err = incompatible_options(opt, Gretlopt.OPT_T | Gretlopt.OPT_B);
    if (*err) {
        return null;
    }

    if (opt & Gretlopt.OPT_B) {
        nadd *= 2;
    }

    newlist = gretl_list_new(norig + nadd);

    for (i=0; i<=norig; i++) {
        newlist[i] = orig[i];
    }

    for (i=1; i<=add[0]; i++) {
        if (in_ivreg_list(orig, add[i], opt)) {
            *err = GretlError.adddup;
        } else {
            add_tsls_var(newlist, add[i], opt);
        }
    }

    if (*err) {
        free(newlist);
        newlist = null;
    }

    return newlist;
}

/* prepend constant to list of instruments */

int zlist_prepend_const(int ** pzlist)
{
    int[2] clist = [1, 0];

    return gretl_list_insert_list(pzlist, clist.ptr, 1);
}

/**
 * tsls_make_endolist:
 * @reglist: regression specification.
 * @instlist: predetermined and exogenous variables.
 * @err: location to receive error code.
 *
 * Determines which variables in @reglist, when checked against the
 * predetermined and exogenous vars in @instlist, need to be
 * instrumented, and populates the returned list accordingly.
 *
 * Returns: allocated list of variables to be instrumented or
 * NULL if there are no such variables.
 */

int * tsls_make_endolist(const(int) * reglist, const(int) * instlist, int * err)
{
    int * ret = null;
    int i, nendo = 0;

    for (i=2; i<=reglist[0]; i++) {
        if (!in_gretl_list(instlist, reglist[i])) {
        nendo++;
    }
    }

    if (nendo == 0) {
    return ret;
    }

    ret = gretl_list_new(nendo);
    if (ret == null) {
    *err = GretlError.alloc;
    } else {
    int vi, j = 1;

    for (i=2; i<=reglist[0]; i++) {
        vi = reglist[i];
        if (!in_gretl_list(instlist, vi)) {
        ret[j++] = vi;
        }
    }
    }

    return ret;
}

/* Fill the residuals matrix for tsls likelihood calculation.
   We run a multivariate regression of the dependent variable
   and the endogenous regressors on all the instruments.
*/

int fill_E_matrix(GretlMatrixRaw * E, IvInfo * ivi, const(DatasetRaw) * dset)
{
    Model * pmod = ivi.pmod;
    GretlMatrixRaw * Y = null;
    GretlMatrixRaw * X = null;
    int ny;
    int * ylist;
    int i, err = 0;

    /* form list of LHS series */
    ny = ivi.endolist[0] + 1;
    ylist = gretl_list_new(ny);
    ylist[1] = ivi.reglist[1];
    for (i=2; i<=ny; i++) {
        ylist[i] = ivi.endolist[i-1];
    }

    /* LHS: depvar and endogenous regressors */
    Y = gretl_matrix_data_subset_masked(ylist, dset, pmod.t1, pmod.t2,
                                        ivi.missmask, &err);
    if (!err) {
        /* RHS: all instruments */
        X = gretl_matrix_data_subset_masked(ivi.instlist, dset,
                        pmod.t1, pmod.t2,
                                            ivi.missmask, &err);
    }

    if (!err) {
        /* populate the residuals matrix, E */
        err = gretl_matrix_multi_ols(Y, X, null, E, null);
    }

    gretl_matrix_free(Y);
    gretl_matrix_free(X);
    free(ylist);

    return err;
}

/* calculate the log determinant */

double tsls_get_ldet(GretlMatrixRaw * S, int * err)
{
    double ret = NADBL;
    char * mask;

    /* allow for the possibility that S is rank-deficient? */
    mask = gretl_matrix_rank_mask(S, err);
    if (mask != null) {
        *err = gretl_matrix_cut_rows_cols(S, mask);
        free(mask);
    }
    if (!*err) {
        ret = gretl_matrix_log_determinant(S, err);
    }

    return ret;
}

/* calculate the maximized log-likelihood for a just-identified
   model */

int tsls_loglik(IvInfo * ivi, DatasetRaw * dset)
{
    GretlMatrixRaw * E, S;
    int T = ivi.pmod.nobs;
    int k = 1 + ivi.endolist[0];
    int err = 0;

    ivi.pmod.lnL = NADBL;

    E = gretl_matrix_alloc(T, k);
    S = gretl_matrix_alloc(k, k);

    if (E == null || S == null) {
        err = GretlError.alloc;
    } else {
        err = fill_E_matrix(E, ivi, dset);
    }

    if (!err) {
        err = gretl_matrix_multiply_mod(E, GretlMatrixMod.transpose,
                                        E, GretlMatrixMod.none,
                                        S, GretlMatrixMod.none);
    }

    if (!err) {
        double ldet = tsls_get_ldet(S, &err);

        if (na(ldet)) {
            ivi.pmod.lnL = NADBL;
        } else {
            /* Davidson and MacKinnon, ETM, p. 538, taking
               kappa = 1 for the just-identified case
            */
            ivi.pmod.lnL = -(T / 2.0) * (ConstId.pi + ldet);
        }
    }

    mle_criteria(ivi.pmod, 0);

    gretl_matrix_free(E);
    gretl_matrix_free(S);

    return err;
}

/* perform the Sargan overidentification test for an IV model */

int ivreg_sargan_test(IvInfo * ivi, DatasetRaw * dset)
{
    Model * pmod = ivi.pmod;
    Model smod;
    int i, t, err = 0;
    int * OT_list = null;
    int nv = dset.v;

    if (pmod.list[0] == 2 && pmod.list[2] == 0) {
        /* degenerate model with const only */
        return 0;
    }

    err = dataset_add_series(dset, 1);
    if (err) {
        return err;
    }

    /* add the @pmod residual to the dataset */
    for (t=pmod.t1; t<=pmod.t2; t++) {
        dset.Z[nv][t] = pmod.uhat[t];
    }

    OT_list = gretl_list_new(ivi.instlist[0] + 1);
    if (OT_list == null) {
        dataset_drop_last_variables(dset, 1);
        return GretlError.alloc;
    }

    OT_list[1] = nv;
    for (i=2; i<=OT_list[0]; i++) {
        OT_list[i] = ivi.instlist[i-1];
    }

    /* OLS: the dependent variable is the residual from @pmod;
       the regressors are the (full set of) instruments
    */
    smod = lsq(OT_list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);

    if (smod.errcode) {
        fprintf(stderr, "ivreg_sargan_test: smod.errcode = %d\n".ptr, smod.errcode);
        err = smod.errcode;
    } else {
        ModelTest *test = model_test_new(ModelTestType.sargan);
        double OTest = smod.rsq * smod.nobs;

        if (test != null) {
            model_test_set_teststat(test, TestStats.lm);
            model_test_set_dfn(test, ivi.overid);
            model_test_set_value(test, OTest);
            model_test_set_pvalue(test, chisq_cdf_comp(ivi.overid, OTest));
            maybe_add_test_to_model(pmod, test);
        }
    }

    clear_model(&smod);
    dataset_drop_last_variables(dset, 1);
    free(OT_list);

    return err;
}

int hausman_drop_check(Model * hmod, int * hatlist, int * endolist, DatasetRaw * dset)
{
    int * dlist = gretl_model_get_list(hmod, "droplist");
    int n_drop = 0;

    if (dlist != null) {
        int i, pos;

        for (i=1; i<=dlist[0]; i++) {
            pos = in_gretl_list(hatlist, dlist[i]);
            if (pos > 0 && pos <= endolist[0]) {
                n_drop++;
                gretl_warnmsg_sprintf("%s does not seem to be endogenous".ptr,
                                      dset.varname[endolist[pos]]);
            }
        }
    }

    return n_drop;
}

double tsls_robust_hausman_fixup(const(int) * hatlist, Model * pmod, int * df)
{
    int * wlist = gretl_list_copy(hatlist);
    double H = NADBL;

    if (wlist != null) {
        int i;

        for (i=wlist[0]; i>0; i--) {
            if (!in_gretl_list(pmod.list, wlist[i])) {
                gretl_list_delete_at_pos(wlist, i);
            }
        }
        H = wald_omit_chisq(wlist, pmod);
        *df = wlist[0];
        free(wlist);
    }

    return H;
}

void add_hausman_to_model(Model * pmod, double HTest, int df)
{
    ModelTest *test = model_test_new(ModelTestType.iv_hausman);

    if (test != null) {
        model_test_set_teststat(test, TestStats.wald_chisq);
        model_test_set_dfn(test, df);
        model_test_set_value(test, HTest);
        model_test_set_pvalue(test, chisq_cdf_comp(df, HTest));
        maybe_add_test_to_model(pmod, test);
    }
}

/* Hausman test via the matrix method (both plain and robust).  See
   Alecos Papadopoulos, 'A New Matrix Statistic for the Hausman
   Endogeneity Test under Heteroskedasticity', Econometrics 11:23,
   2023. https://doi.org/10.3390/econometrics11040023
*/

int tsls_matrix_hausman_test(IvInfo * ivi, Gretlopt opt, DatasetRaw * dset)
{
    const(char) * mask = ivi.missmask;
    Model * pmod = ivi.pmod;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X1 = null;
    GretlMatrixRaw * Z = null;
    GretlMatrixRaw * u = null;
    GretlMatrixRaw * E = null;
    GretlMatrixRaw * X1u = null;
    GretlMatrixRaw * A = null;
    int * xlist = null;
    int i;
    int[2] ylist;
    double q;
    int err = 0;

    /* Step 1: obtain residuals from estimation of the specified
       model via plain OLS.
    */

    xlist = gretl_list_new(ivi.reglist[0] - 1);
    if (xlist == null) {
    err = GretlError.alloc;
    goto bailout;
    }
    for (i=1; i<=xlist[0]; i++) {
    xlist[i] = ivi.reglist[i+1];
    }
    ylist[0] = 1;
    ylist[1] = ivi.reglist[1];

    /* all regressors */
    X = gretl_matrix_data_subset_masked(xlist, dset, pmod.t1, pmod.t2,
                    mask, &err);
    if (!err) {
    y = gretl_matrix_data_subset_masked(ylist.ptr, dset, pmod.t1, pmod.t2,
                        mask, &err);
    }
    if (!err) {
    /* put OLS residuals into @u */
    u = gretl_matrix_alloc(X.rows, 1);
    gretl_matrix_ols(y, X, null, null, u, null);
    }

    /* Step 2: fill @X1 with \hat{X1} */
    X1 = gretl_matrix_data_subset_masked(ivi.hatlist, dset, pmod.t1, pmod.t2,
                     mask, &err);

    if (!err) {
    /* some more required allocation */
    X1u = gretl_matrix_alloc(X1.cols, 1);
    E = gretl_matrix_alloc(X1.rows, X1.cols);
    A = gretl_matrix_alloc(X1.cols, X1.cols);
    if (X1u == null || E == null || A == null) {
        err = GretlError.alloc;
    }
    }

    if (err) {
    goto bailout;
    }

    /* Step 3: form X1u = \hat{X}_1'\hat{u} */
    gretl_matrix_multiply_mod(X1, GretlMatrixMod.transpose,
                  u, GretlMatrixMod.none,
                  X1u, GretlMatrixMod.none);

    /* Step 4: form E = M_x \hat{X}_1 */
    err = gretl_matrix_multi_ols(X1, X, null, E, null);

    if (!err) {
    if (opt & Gretlopt.OPT_R) {
        /* Step 5 (robust): A = qform(M_x \hat{X_1}', diag(u.^2)) */
        double ui;

        for (i=0; i<u.rows; i++) {
        ui = u.val[i];
        u.val[i] = ui * ui;
        }
        gretl_matrix_diag_qform(E, GretlMatrixMod.transpose, u,
                    A, GretlMatrixMod.none);
    } else {
        /* Step 5 (plain): A = \hat{X_1}' M_x \hat{X_1} */
        double s2 = 0;

        for (i=0; i<u.rows; i++) {
        s2 += u.val[i] * u.val[i];
        }
        s2 /= u.rows;
        gretl_matrix_multiply_mod(X1, GretlMatrixMod.transpose,
                      E, GretlMatrixMod.none,
                      A, GretlMatrixMod.none);
        gretl_matrix_multiply_by_scalar(A, s2);
    }
    err = gretl_invert_symmetric_matrix(A);
    if (err) {
        gretl_errmsg_set("Hausman test matrix could not be inverted: this indicates invalid instruments");
    }
    }

    if (err) {
    q = NADBL;
    } else {
    q = gretl_scalar_qform(X1u, A, &err);
    if (!err) {
        add_hausman_to_model(pmod, q, X1.cols);
    }
    }

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(X1);
    gretl_matrix_free(Z);
    gretl_matrix_free(u);
    gretl_matrix_free(E);
    gretl_matrix_free(X1u);
    gretl_matrix_free(A);
    free(xlist);

    return err;
}

/* Hausman test via the regression method */

int tsls_hausman_test(IvInfo * ivi, Gretlopt opt, DatasetRaw * dset)
{
    Model hmod;
    double URSS;
    Gretlopt uopt;
    int * HT_list = null;
    int t, df;
    int nv = dset.v;
    int n_drop = 0;
    int ku = 0;
    int err = 0;

    /* create regression list for unrestricted model by appending
       @hatlist to @reglist */
    HT_list = gretl_list_add(ivi.reglist, ivi.hatlist, &err);
    if (err == GretlError.noadd) {
        err = 0; /* more of a no-op than an error */
    } else if (err) {
        return err;
    }

    /* option for unrestricted regression */
    uopt = (opt & Gretlopt.OPT_R)? (Gretlopt.OPT_A | Gretlopt.OPT_R) : Gretlopt.OPT_A;

    /* estimate the unrestricted model */
    hmod = lsq(HT_list, dset, GretlCmdIndex.ols, uopt);
    if (hmod.errcode) {
        fprintf(stderr, "tsls_hausman_test: hmod.errcode (U) = %d\n".ptr, hmod.errcode);
        err = hmod.errcode;
        goto bailout;
    }

    n_drop = hausman_drop_check(&hmod, ivi.hatlist, ivi.endolist, dset);

    if (hmod.dfd == 0) {
        /* perfect fit, can't do test */
        goto bailout;
    }

    if (opt & Gretlopt.OPT_R) {
        /* 2023-08-31: if the --robust option was in force for
           tsls estimation, compute the Hausman test as a Wald
           test on the robust variance matrix for the U model.
        */
        double HTest = NADBL;

        if (n_drop > 0) {
            HTest = tsls_robust_hausman_fixup(ivi.hatlist, &hmod, &df);
        } else {
            HTest = wald_omit_chisq(ivi.hatlist, &hmod);
            df = ivi.hatlist[0];
        }
        if (!na(HTest)) {
            add_hausman_to_model(ivi.pmod, HTest, df);
        }
        goto bailout;
    }

    /* record U-model info */
    URSS = hmod.ess;
    ku = hmod.ncoeff;

    /* add fitted values from unrestricted model to dataset */
    err = dataset_add_series(dset, 1);
    if (err) {
        goto bailout;
    } else {
        for (t=hmod.t1; t<=hmod.t2; t++) {
            dset.Z[nv][t] = hmod.yhat[t];
        }
    }

    clear_model(&hmod);

    /* redefine HT_list: dep var is the fitted value from the
       U-model, regressors are from @reglist */
    free(HT_list);
    HT_list = gretl_list_copy(ivi.reglist);
    if (HT_list == null) {
        err = GretlError.alloc;
        goto bailout;
    }
    HT_list[1] = nv;

    /* regress the U fitted values on the regressors of the restricted
       model (so the sum of squares equals RRSS-URSS)
    */
    hmod = lsq(HT_list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);

    if (hmod.errcode) {
        fprintf(stderr, "tsls_hausman_test: hmod.errcode (D) = %d\n".ptr, hmod.errcode);
        err = hmod.errcode;
    } else if ((df = ku - hmod.ncoeff) > 0) {
        /* Degrees of freedom for the Hausman test need not be equal to the
           number of added regressors, since some of them may not really be
           endogenous.
        */
        double DRSS = hmod.ess;
        double HTest = (DRSS/URSS) * hmod.nobs;

        add_hausman_to_model(ivi.pmod, HTest, df);
    }

 bailout:

    dataset_drop_last_variables(dset, dset.v - nv);
    clear_model(&hmod);
    free(HT_list);

    return err;
}

/* form matrix of instruments and perform QR decomposition
   of this matrix; attach return Q to @ivi; return 0 on
   success, error code on error.
*/

int tsls_add_Q(IvInfo * ivi, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * R = null;
    int rank, ndrop = 0;
    int * dlist = null;
    double test;
    int i, j, k;
    int err = 0;

    ivi.Q = gretl_matrix_data_subset_masked(ivi.instlist, dset,
                         dset.t1, dset.t2,
                         ivi.missmask, &err);
    if (err) {
        return err;
    }

    k = gretl_matrix_cols(ivi.Q);

    if (k > ivi.Q.rows) {
        /* can't do QR decomp! */
        gretl_errmsg_sprintf(_("Number of instruments (%d) exceeds the number of observations (%d)").ptr,
                             k, ivi.Q.rows);
        err = GretlError.df;
        goto bailout;
    }

    R = gretl_matrix_alloc(k, k);
    if (R == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    err = gretl_matrix_QR_decomp(ivi.Q, R);
    if (err) {
        goto bailout;
    }

    rank = gretl_check_QR_rank(R, &err, null);

    if (err) {
        goto bailout;
    } else if (rank < k) {
        fprintf(stderr, "tsls_Q: k = %d, rank = %d\n".ptr, k, rank);
        ndrop = k - rank;
    }

    if (ndrop > 0) {
        dlist = gretl_list_new(ndrop);
        if (dlist != null) {
            dlist[0] = 0;
        }
        j = 1;
        for (i=0; i<k; i++) {
            test = gretl_matrix_get(R, i, i);
            if (fabs(test) < R_DIAG_MIN) {
                if (dlist != null) {
                    dlist[0] += 1;
                    dlist[dlist[0]] = ivi.instlist[j];
                }
                fprintf(stderr, "tsls_Q: dropping redundant instrument %d (%s)\n".ptr,
                        ivi.instlist[j], dset.varname[ivi.instlist[j]]);
                gretl_list_delete_at_pos(ivi.instlist, j--);
            }
            j++;
        }

        k = ivi.instlist[0];
        gretl_matrix_free(ivi.Q);
        ivi.Q = gretl_matrix_data_subset_masked(ivi.instlist, dset,
                         dset.t1, dset.t2,
                         ivi.missmask, &err);
        if (!err) {
            R = gretl_matrix_reuse(R, k, k);
            err = gretl_matrix_QR_decomp(ivi.Q, R);
        }
    }

 bailout:

    gretl_matrix_free(R);

    if (err) {
        free(dlist);
        gretl_matrix_free(ivi.Q);
        ivi.Q = null;
    } else {
        ivi.idroplist = dlist;
    }

    return err;
}

int tsls_form_xhat(IvInfo * ivi, DatasetRaw * dset, int v0, int v1)
{
    const(char) * mask = ivi.missmask;
    const(double) * x = dset.Z[v0];
    double * xhat = dset.Z[v1];
    int k = gretl_matrix_cols(ivi.Q);
    int allzero = 1;
    int i, t, s = 0;
    int err = 0;

    /* form g = Q'x */
    for (i=0; i<k; i++) {
        ivi.g.val[i] = 0.0;
        s = 0;
        for (t=dset.t1; t<=dset.t2; t++) {
            if (mask != null && mask[t - dset.t1]) {
                continue;
            }
            ivi.g.val[i] += gretl_matrix_get(ivi.Q, s++, i) * x[t];
        }
    }

    /* form \hat{x} = Qg = QQ'x */
    s = 0;
    for (t=dset.t1; t<=dset.t2; t++) {
        if (mask != null && mask[t - dset.t1]) {
            xhat[t] = NADBL;
            continue;
        }
        xhat[t] = 0.0;
        for (i=0; i<k; i++) {
            xhat[t] += gretl_matrix_get(ivi.Q, s, i) * ivi.g.val[i];
        }
        if (xhat[t] != 0) {
            allzero = 0;
        }
        s++;
    }

    if (allzero) {
        gretl_errmsg_sprintf(_("The first-stage fitted values for %s are all zero").ptr,
                             dset.varname[v0]);
        err = GretlError.data;
    } else {
        /* name the fitted series according to the original */
        strcpy(dset.varname[v1], "h_");
        strncat(dset.varname[v1], dset.varname[v0], VNAMELEN - 3);
    }

    return err;
}

void tsls_residuals(IvInfo * ivi, Gretlopt opt, const(DatasetRaw) * dset)
{
    Model * pmod = ivi.pmod;
    double yh, sigma0 = pmod.sigma;
    int yno = ivi.reglist[1];
    int i, t;

    pmod.ess = 0.0;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        yh = 0.0;
        for (i=0; i<pmod.ncoeff; i++) {
            yh += pmod.coeff[i] * dset.Z[ivi.reglist[i+2]][t];
        }
        pmod.yhat[t] = yh;
        pmod.uhat[t] = dset.Z[yno][t] - yh;
        pmod.ess += pmod.uhat[t] * pmod.uhat[t];
    }

    if (pmod.ess <= 0.0) {
        pmod.sigma = 0.0;
    } else {
        int den = (opt & Gretlopt.OPT_N)? pmod.nobs : pmod.dfd;

        pmod.sigma = sqrt(pmod.ess / den);
    }

    if (sigma0 > 0.0) {
        double corr = pmod.sigma / sigma0;

        for (i=0; i<pmod.ncoeff; i++) {
            pmod.sderr[i] *= corr;
        }
    }
}

void tsls_extra_stats(IvInfo * ivi, Gretlopt opt, const(DatasetRaw) * dset)
{
    Model * pmod = ivi.pmod;
    int yno = ivi.reglist[1];
    double r;

    /* squared correlation between y and yhat */
    pmod.rsq = gretl_corr_rsq(pmod.t1, pmod.t2, dset.Z[yno], pmod.yhat);

    r = 1 - pmod.rsq;
    pmod.adjrsq = 1.0 - (r * (pmod.nobs - 1) / pmod.dfd);

    pmod.fstt = pmod.chisq = NADBL;

    if (pmod.ncoeff > 1) {
        if ((opt & Gretlopt.OPT_N) || (ivi.overid && pmod.ncoeff == 2)) {
            pmod.chisq = wald_omit_chisq(null, pmod);
        } else {
            pmod.fstt = wald_omit_F(null, pmod);
        }
    }

    /* tsls is not an ML estimator unless it's exactly identified,
       and in that case we require a special calculation */
    pmod.lnL = NADBL;
    pmod.criterion[InfoCriterion.aic] = NADBL;
    pmod.criterion[InfoCriterion.bic] = NADBL;
    pmod.criterion[InfoCriterion.hqc] = NADBL;

    if (dataset_is_time_series(dset) && !model_has_missing_obs(pmod)) {
        /* time series, no missing obs within sample range */
        pmod.rho = rhohat(1, pmod.t1, pmod.t2, pmod.uhat);
        pmod.dw = dwstat(1, pmod, dset);
    } else {
        pmod.rho = pmod.dw = NADBL;
    }
}

void tsls_recreate_full_list(IvInfo * ivi)
{
    int * full_list;

    if (ivi.instlist != null && ivi.instlist[0] > 0) {
        full_list = gretl_lists_join_with_separator(ivi.reglist, ivi.instlist);
    } else {
        full_list = gretl_list_copy(ivi.reglist);
    }

    if (full_list == null) {
        ivi.pmod.errcode = GretlError.alloc;
    } else {
        free(ivi.pmod.list);
        ivi.pmod.list = full_list;
    free(ivi.pmod.xlist); /* don't leak memory */
    ivi.pmod.xlist = gretl_list_sublist(ivi.reglist, 2, -1);
    }
}

void replace_list_element(int * list, int targ, int repl)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == targ) {
            list[i] = repl;
            break;
        }
    }
}

int reglist_remove_redundant_vars(IvInfo * ivi)
{
    int * dlist = gretl_model_get_list(ivi.pmod, "droplist");
    int i, pos;

    if (dlist == null) {
        return 1;
    }

    for (i=1; i<=dlist[0]; i++) {
        pos = in_gretl_list(ivi.s2list, dlist[i]);
        if (pos > 1) {
            gretl_list_delete_at_pos(ivi.reglist, pos);
        }
        if (ivi.endolist != null) {
            /* note: if endolist is NULL, so is hatlist */
            pos = in_gretl_list(ivi.hatlist, dlist[i]);
            if (pos > 1) {
                /* First replace dlist[i] with the ID of the original
                   regressor from @endolist, in place of the "hatlist"
                   variable (which will be deleted when the model is
                   returned), so that the printout of regressors that are
                   dropped due to exact collinearity will show the
                   appropriate names.

                   Second, delete the redundant term from both endolist
                   and hatlist, so that subsequent calculations will not
                   get messed up.
                */
                dlist[i] = ivi.endolist[pos];
                gretl_list_delete_at_pos(ivi.endolist, pos);
                gretl_list_delete_at_pos(ivi.hatlist, pos);
            }
        }
    }

    return 0;
}

/* Stock-Yogo: compute test statistic for weakness of instruments.
   This is the Cragg-Donald statistic: the minimum eigenvalue of the
   matrix defined below (which reduces to the first-stage F-statistic
   when there is just one endogenous regressor).  See Stock and Yogo,
   "Testing for Weak Instruments in Linear IV Regressions" (originally
   NBER Technical Working Paper 284; revised February 2003), at
   http://ksghome.harvard.edu/~JStock/pdf/rfa_6.pdf
*/

int compute_stock_yogo(IvInfo * ivi, const(DatasetRaw) * dset)
{
    const(double) ** dZ = cast(const(double) **) dset.Z;
    Model * pmod = ivi.pmod;
    GretlMatrixBlock * B = null;
    GretlMatrixBlock * B2 = null;
    GretlMatrixRaw * G, S, Y;
    GretlMatrixRaw * Ya, YPY;
    GretlMatrixRaw * X, Z, E = null;
    int T = pmod.nobs;
    int n = ivi.endolist[0];
    int K1 = pmod.ncoeff - n;
    int K2 = ivi.instlist[0] - K1;
    int i, s, vi, vj, t;
    int err = 0;

    B = gretl_matrix_block_new(&G, n, n,
                               &S, n, n,
                               &Y, T, n,
                               &Ya, T, n,
                               &YPY, n, n,
                               null);
    if (B == null) {
        return GretlError.alloc;
    }

    if (K1 > 0) {
        /* we have some exogenous regressors */
        GretlMatrixRaw * M;
        int j, ix = 0, iz = 0;

        B2 = gretl_matrix_block_new(&X, T, K1,
                                    &Z, T, K2,
                                    &E, T, K2,
                                    null);
        if (B2 == null) {
            gretl_matrix_block_destroy(B);
            return GretlError.alloc;
        }

        /* form the X and Z matrices */
        for (i=1; i<=ivi.instlist[0]; i++) {
            vi = ivi.instlist[i];
            if (in_gretl_list(ivi.pmod.list, vi)) {
                /* included exogenous var -> X */
                M = X;
                j = ix++;
            } else {
                /* excluded exogenous var -> Z */
                M = Z;
                j = iz++;
            }
            s = 0;
            for (t=pmod.t1; t<=pmod.t2; t++) {
                if (!na(pmod.uhat[t])) {
                    gretl_matrix_set(M, s++, j, dZ[vi][t]);
                }
            }
        }
    }

    /* form Y matrix */
    for (i=0; i<n; i++) {
        vi = ivi.endolist[i+1];
        s = 0;
        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (!na(pmod.uhat[t])) {
                gretl_matrix_set(Y, s++, i, dZ[vi][t]);
            }
        }
    }

    if (K1 == 0) {
        /* form Y-hat matrix in Ya */
        for (i=0; i<n; i++) {
            vi = ivi.hatlist[i+1];
            s = 0;
            for (t=pmod.t1; t<=pmod.t2; t++) {
                if (!na(pmod.uhat[t])) {
                    gretl_matrix_set(Ya, s++, i, dZ[vi][t]);
                }
            }
        }
    } else {
        int bdim = (K1 > K2)? K1 : K2;
        GretlMatrixRaw * b;

        /* the leading dimension of b is the number of regressors */
        b = gretl_matrix_alloc(bdim, K2);
        if (b == null) {
            err = GretlError.alloc;
        }

        if (!err) {
            /* partial X out of Y: Y <- M_x Y */
            gretl_matrix_reuse(b, K1, n);
            gretl_matrix_reuse(E, -1, n);
            err = gretl_matrix_multi_ols(Y, X, b, E, null);
            gretl_matrix_copy_values(Y, E);
        }

        if (!err) {
            /* partial X out of Z: Z <- M_x Z */
            gretl_matrix_reuse(b, K1, K2);
            gretl_matrix_reuse(E, -1, K2);
            err = gretl_matrix_multi_ols(Z, X, b, E, null);
            gretl_matrix_copy_values(Z, E);
        }

        if (!err) {
            /* form projection P_z Y, in Ya */
            gretl_matrix_reuse(b, K2, n);
            gretl_matrix_reuse(E, -1, n);
            err = gretl_matrix_multi_ols(Y, Z, b, E, null);
            gretl_matrix_copy_values(Ya, Y);
            gretl_matrix_subtract_from(Ya, E);
        }

        gretl_matrix_free(b);
    }

    if (!err) {
        /* form Y' P_z Y */
        err = gretl_matrix_multiply_mod(Y, GretlMatrixMod.transpose,
                                        Ya, GretlMatrixMod.none,
                                        YPY, GretlMatrixMod.none);
    }

    /* now write first-stage residuals into Ya */
    for (i=0; i<n && !err; i++) {
        vi = ivi.endolist[i+1];
        vj = ivi.hatlist[i+1];
        s = 0;
        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (!na(pmod.uhat[t])) {
                gretl_matrix_set(Ya, s++, i, dZ[vi][t] - dZ[vj][t]);
            }
        }
    }

    if (!err) {
        /* form S = \hat{\Sigma}_{vv} */
        err = gretl_matrix_multiply_mod(Y, GretlMatrixMod.transpose,
                                        Ya, GretlMatrixMod.none,
                                        S, GretlMatrixMod.none);
        if (!err) {
            gretl_matrix_xtr_symmetric(S);
            gretl_matrix_divide_by_scalar(S, T - K1 - K2);
            for (i=0; i<S.rows; i++) {
                if (gretl_matrix_get(S, i, i) <= 0) {
                    err = GretlError.singular;
                    break;
                }
            }
        }
    }

    if (!err) {
        double rc = gretl_symmetric_matrix_rcond(S, &err);

        if (!err && (na(rc) || rc < 1.0e-7)) {
            err = GretlError.singular;
        }
    }

    if (!err) {
        /* S^{-1/2}: invert S and Cholesky-decompose */
        err = gretl_invert_symmetric_matrix(S);
        if (!err) {
            err = gretl_matrix_cholesky_decomp(S);
        }
    }

    if (!err) {
        /* finally, form big sandwich */
        err = gretl_matrix_qform(S, GretlMatrixMod.transpose,
                                 YPY, G, GretlMatrixMod.none);
        if (!err) {
            gretl_matrix_divide_by_scalar(G, K2);
        }
    }

    if (!err) {
        double gmin = gretl_symm_matrix_lambda_min(G, &err);

        if (!err) {
            gretl_model_set_double(pmod, "gmin", gmin);
            gretl_model_set_int(pmod, "n", n);
            gretl_model_set_int(pmod, "K2", K2);
        }
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(B2);

    return err;
}

/* Calculate first-stage F-test as diagnostic for weak instruments.
   This variant can only handle one endogenous regressor, but it
   supports robust estimation if the model in question uses a robust
   covariance estimator.
*/

int compute_first_stage_F(IvInfo * ivi, Gretlopt opt, DatasetRaw * dset)
{
    Model mod1;
    int * list1 = null;
    int * flist = null;
    double F;
    int i, ev, vi;
    int err = 0;

    gretl_model_init(&mod1, dset);
    ev = ivi.endolist[1];

    /* The regression list, list1, has endogenous regressor ev as the
       dependent variable and includes all instruments; the list for
       omission, flist, includes all the excluded instruments.
    */

    list1 = gretl_list_append_term(&list1, ev);
    if (list1 == null) {
        return GretlError.alloc;
    }

    for (i=1; i<=ivi.instlist[0]; i++) {
        vi = ivi.instlist[i];
        list1 = gretl_list_append_term(&list1, vi);
        if (list1 == null) {
            err = GretlError.alloc;
            break;
        }
        if (!in_gretl_list(ivi.reglist, vi)) {
            flist = gretl_list_append_term(&flist, vi);
            if (flist == null) {
                err = GretlError.alloc;
                break;
            }
        }
    }

    if (!err) {
        Gretlopt myopt = (opt & Gretlopt.OPT_R)? Gretlopt.OPT_R : Gretlopt.OPT_NONE;

        mod1 = lsq(list1, dset, GretlCmdIndex.ols, myopt | Gretlopt.OPT_A);
        err = mod1.errcode;
        if (err) {
            fprintf(stderr, "compute_first_stage F: lsq failed\n".ptr);
        }
    }

    if (!err) {
        F = wald_omit_F(flist, &mod1);
        if (na(F)) {
            err = 1;
        }
    }

    if (!err) {
        gretl_model_set_double(ivi.pmod, "stage1-F", F);
        gretl_model_set_int(ivi.pmod, "stage1-dfn", flist[0]);
        gretl_model_set_int(ivi.pmod, "stage1-dfd", mod1.dfd);
        if (!(opt & Gretlopt.OPT_R)) {
            /* flag lookup of Stock-Yogo critical values */
            gretl_model_set_double(ivi.pmod, "gmin", F);
        }
    }

    clear_model(&mod1);
    free(list1);
    free(flist);

    return err;
}

int tsls_adjust_sample(const(int) * list, DatasetRaw * dset, char ** pmask)
{
    int i, t, t1min = dset.t1, t2max = dset.t2;
    char * mask = null;
    int T, vi, missobs;
    int err = 0;

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
        missobs = 0;
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset.Z[vi][t])) {
                missobs = 1;
                break;
            }
        }
        if (missobs) {
            t1min++;
        } else {
            break;
        }
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
        missobs = 0;
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset.Z[vi][t])) {
                missobs = 1;
                break;
            }
        }
        if (missobs) {
            t2max--;
        } else {
            break;
        }
    }

    /* count missing values in mid-range of data */
    missobs = 0;
    for (t=t1min; t<=t2max; t++) {
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == 0 || vi == LISTSEP) {
                continue;
            }
            if (na(dset.Z[vi][t])) {
                missobs++;
                break;
            }
        }
    }

    T = t2max - t1min + 1;

    if (missobs == T) {
        err = GretlError.missdata;
    } else if (missobs > 0) {
        mask = cast(typeof(mask)) calloc(T, 1); /* all NUL bytes */
        if (mask == null) {
            err = GretlError.alloc;
        } else {
            for (t=t1min; t<=t2max; t++) {
                for (i=1; i<=list[0]; i++) {
                    vi = list[i];
                    if (vi == 0 || vi == LISTSEP) {
                        continue;
                    }
                    if (na(dset.Z[vi][t])) {
                        mask[t - t1min] = 1;
                        break;
                    }
                }
            }
        }
    }

    if (!err) {
        dset.t1 = t1min;
        dset.t2 = t2max;
    }

    *pmask = mask;

    return err;
}

/**
 * ivreg_process_lists:
 * @list: original specification list.
 * @reglist: location to receive regression list.
 * @instlist: location to receive list of instruments.
 *
 * Split the incoming list into its two components and perform
 * some basic checks; if the checks fail the two created
 * lists are destroyed.
 *
 * Returns: 0, on success, non-zero error code on failure.
 */

int ivreg_process_lists(const(int) * list, int ** reglist, int ** instlist)
{
    int * rlist = null;
    int * zlist = null;
    int i, oc, err;

    err = gretl_list_split_on_separator(list, &rlist, &zlist);
    if (err) {
        return err;
    }

    if (rlist[0] < 2 || zlist == null || zlist[0] < 1) {
    /* The regression list must have at least two members, and the
       instrument list must not be empty
    */
        err = GretlError.args;
    } else {
        for (i=1; i<=zlist[0] && !err; i++) {
            if (zlist[i] == rlist[1]) {
                gretl_errmsg_set(_("You can't use the dependent variable as an instrument").ptr);
                err = GretlError.data;
            }
        }
    }

    if (!err && in_gretl_list(rlist, 0) > 1 && !in_gretl_list(zlist, 0)) {
    /* Don't treat the constant as if it were endogenous: if
       it's in the regressor list ensure it's also in the
       instrument list.
    */
    err = zlist_prepend_const(&zlist);
    }

    if (!err && (oc = zlist[0] - rlist[0] + 1) < 0) {
    gretl_errmsg_sprintf(_("The order condition for identification is not satisfied.\nAt least %d more instruments are needed.").ptr, -oc);
    err = GretlError.data;
    }

    if (err) {
    free(rlist);
    free(zlist);
    } else {
    *reglist = rlist;
    *instlist = zlist;
    }

    return err;
}

/**
 * tsls:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: may contain %OPT_R for robust VCV, %OPT_N for no df correction,
 * %OPT_A if this is an auxiliary reression, %OPT_E if we're
 * estimating one equation within a system, %OPT_H to
 * add "hatlist" to model even under %OPT_E, %OPT_X to skip
 * the set of tests that are usually calculated. Also %OPT_C
 * for clustered standard errors.
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.  If %OPT_E is given, fitted values from the first-stage
 * regressions are saved with the model.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model tsls(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    Model tsls;
    IvInfo *ivi = null;
    int orig_t1 = dset.t1, orig_t2 = dset.t2;
    int orig_nvar = dset.v;
    int sysest = (opt & Gretlopt.OPT_E);
    int no_tests = 0;
    int no_hausman = 0;
    int nendo = 0;
    int nreg = 0;
    int i, err = 0;

    /* initialize model in case we bail out early on error */
    gretl_model_init(&tsls, dset);

    /* allocate container struct */
    ivi = iv_info_new(&tsls);
    if (ivi == null) {
    tsls.errcode = GretlError.alloc;
    return tsls;
    }

    if (opt & Gretlopt.OPT_X) {
        /* supplementary tests wanted? */
        const(char) * s = get_optval_string(GretlCmdIndex.ivreg, Gretlopt.OPT_X);

        if (s != null && !strcmp(s, "hausman")) {
            /* skip the hausman test only */
            no_hausman = 1;
        } else {
            /* skip all extra tests */
            no_tests = 1;
        }
    }

    /* reglist: dependent var plus list of regressors
       instlist: list of instruments
    */
    tsls.errcode = ivreg_process_lists(list, &ivi.reglist, &ivi.instlist);
    if (tsls.errcode) {
    iv_info_destroy(ivi);
        return tsls;
    }

    /* adjust sample range for missing observations */
    err = tsls_adjust_sample(list, dset, &ivi.missmask);

    if (!err) {
        /* allocate second stage regression list */
        ivi.s2list = gretl_list_copy(ivi.reglist);
        if (ivi.s2list == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        goto bailout;
    }

    /* in case we drop any instruments as redundant */
    if (sysest) {
        ivi.exolist = gretl_list_copy(ivi.instlist);
    }

    /* drop any vars that are all zero, and reshuffle the constant
       into first position among the independent vars
    */
    tsls_omitzero(ivi.reglist, ivi.missmask, dset);
    reglist_check_for_const(ivi.reglist, dset);
    tsls_omitzero(ivi.instlist, ivi.missmask, dset);

    /* Determine the list of variables (endolist) for which we need to
       obtain fitted values in the first stage.  Note that we accept
       the condition where there are no such variables (and hence
       TSLS is in fact redundant); this may sometimes be required
       when tsls is called in the context of system estimation.
    */
    ivi.endolist = tsls_make_endolist(ivi.reglist,
                       ivi.instlist,
                       &err);
    if (err) {
        goto bailout;
    }

    if (ivi.endolist != null) {
        nendo = ivi.endolist[0];
        ivi.hatlist = gretl_list_copy(ivi.endolist);
        if (ivi.hatlist == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        err = tsls_add_Q(ivi, dset);
    }

    if (err) {
        goto bailout;
    }

    /* check (again) for order condition: we do this after tsls_Q in case
       any instruments get dropped at that stage
    */
    nreg = ivi.reglist[0] - 1;
    ivi.overid = ivi.instlist[0] - nreg;
    if (ivi.overid < 0) {
        gretl_errmsg_sprintf(_("The order condition for identification is not satisfied.\nAt least %d more instruments are needed.").ptr, -ivi.overid);
        err = GretlError.data;
        goto bailout;
    }

    /* allocate storage for fitted vars (etc.), if needed */
    if (nendo > 0) {
        err = dataset_add_series(dset, nendo);
        if (!err) {
            ivi.g = gretl_vector_alloc(ivi.Q.cols);
            if (ivi.g == null) {
                err = GretlError.alloc;
            }
        }
        if (err) {
            goto bailout;
        }
    }

    /* Deal with the variables for which instruments are needed: loop
       across the list of variables to be instrumented (endolist),
       form the fitted values as QQ'x_i, and add these fitted values
       into the data array, dset->Z.
    */

    for (i=0; i<nendo; i++) {
        int v0 = ivi.endolist[i+1];
        int v1 = orig_nvar + i;

        /* form xhat_i = QQ'x_i */
        err = tsls_form_xhat(ivi, dset, v0, v1);
        if (err) {
            goto bailout;
        }

        /* substitute v1 into the right place in the second-stage
           regression list */
        replace_list_element(ivi.s2list, v0, v1);

        /* and update hatlist */
        ivi.hatlist[i+1] = v1;
    }

    /* second-stage regression */
    tsls = lsq(ivi.s2list, dset, GretlCmdIndex.ols, (sysest)? Gretlopt.OPT_Z : Gretlopt.OPT_NONE);
    if (tsls.errcode) {
        fprintf(stderr, "tsls, stage 2: tsls.errcode = %d\n".ptr, tsls.errcode);
        goto bailout;
    }

    if (tsls.list[0] < ivi.s2list[0]) {
        /* Were collinear regressors dropped? If so, adjustments are needed */
        ivi.overid += ivi.s2list[0] - tsls.list[0];
        reglist_remove_redundant_vars(ivi);
        if (ivi.endolist != null) {
            nendo = ivi.endolist[0];
        }
    }

    if (ivi.instlist != null && ivi.instlist[0] > 0) {
        /* record the number of instruments used */
        gretl_model_set_int(&tsls, "ninst", ivi.instlist[0]);
    }

    if (!no_tests) {
        if (!sysest || (opt & Gretlopt.OPT_H)) {
            if (nendo == 1) {
                /* handles robust estimation, for single endogenous regressor */
                compute_first_stage_F(ivi, opt, dset);
            } else if (nendo > 0 && !(opt & Gretlopt.OPT_R)) {
                /* at present, only handles case of i.i.d. errors */
                compute_stock_yogo(ivi, dset);
            }
        }
    }

    if (nendo > 0) {
        /* special: we need to use the original RHS vars to compute
           residuals and associated statistics */
        tsls_residuals(ivi, opt, dset);
    }

    if (opt & Gretlopt.OPT_C) {
        /* clustered implies robust */
        opt |= Gretlopt.OPT_R;
    }

    if (opt & Gretlopt.OPT_R) {
        /* robust standard errors called for */
        qr_tsls_vcv(&tsls, dset, opt);
        if (tsls.errcode) {
            fprintf(stderr, "qr_tsls_vcv: err = %d\n".ptr, tsls.errcode);
            goto bailout;
        }
    }

    if (nendo > 0) {
        /* compute additional statistics (R^2, F, etc.) */
        tsls_extra_stats(ivi, opt, dset);
    }

    if (!sysest && !no_tests) {
        if (!no_hausman && nendo > 0 && ivi.hatlist != null) {
        if (opt & Gretlopt.OPT_M) {
        tsls_matrix_hausman_test(ivi, opt, dset);
        } else {
        tsls_hausman_test(ivi, opt, dset);
        }
        }
        if (ivi.overid > 0) {
            ivreg_sargan_test(ivi, dset);
        } else if (nendo > 0) {
            tsls_loglik(ivi, dset);
        }
    }

    /* set command code on the model */
    tsls.ci = GretlCmdIndex.ivreg;

    /* write the full 2SLS list (dep. var. and regressors, followed by
       instruments, possibly purged of redundant elements) into the model
       for future reference
    */
    tsls_recreate_full_list(ivi);

    if (ivi.idroplist != null) {
        gretl_model_set_list_as_data(&tsls, "inst_droplist", ivi.idroplist);
    ivi.idroplist = null; /* don't free */
    }

 bailout:

    if (err && !tsls.errcode) {
        tsls.errcode = err;
    }

    if (!tsls.errcode && nendo > 0) {
        if (sysest) {
            /* save first-stage fitted values */
            tsls_save_data(ivi, dset);
        }
        if (!sysest || (opt & Gretlopt.OPT_H)) {
            /* save list of endogenous regressors on model */
            gretl_model_set_list_as_data(&tsls, "endolist", ivi.endolist);
            ivi.endolist = null; /* model takes ownership */
            gretl_model_set_list_as_data(&tsls, "instlist", ivi.instlist);
            ivi.instlist = null; /* ditto */
        }
    }

    /* delete any first-stage fitted values from dataset */
    dataset_drop_last_variables(dset, dset.v - orig_nvar);

    iv_info_destroy(ivi);

    if ((opt & Gretlopt.OPT_A) || tsls.errcode) {
        model_count_minus(&tsls); /* OK? */
    }

    if (!tsls.errcode && !(opt & Gretlopt.OPT_N)) {
        gretl_model_set_int(&tsls, "dfcorr", 1);
    }

    /* restore original sample range */
    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    return tsls;
}

