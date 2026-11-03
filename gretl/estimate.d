/* Port of estimate.c */

module gretl.estimate;


import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : exp, fabs, isnan, log, sqrt;
import core.stdc.stdio : FILE, fprintf, fputs, sprintf, stderr;
import core.stdc.stdlib : atoi, free, malloc;
import core.stdc.string : memcpy, strcpy, strstr;

import gretl.gretl_commandsh;
import gretl.dataset : dataset_add_series, dataset_drop_last_variables, series_get_lag, series_get_parent_id, series_is_log;
import gretl.gretl_errorsh;
import gretl.glib;
import gretl.gstrfuncs;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_free, gretl_matrix_ols,
    gretl_matrix_r_squared;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrixh;
import gretl.gretl_model : ModelTest, clear_model, gretl_model_add_arinfo, gretl_model_destroy_data_item, gretl_model_get_double, gretl_model_get_double_default, gretl_model_get_int, gretl_model_init, gretl_model_set_data, gretl_model_set_double, gretl_model_set_int, gretl_model_set_list_as_data, gretl_model_set_string_as_data, gretl_model_test_print_direct, impose_model_smpl, maybe_add_test_to_model, mle_criteria, model_test_free, model_test_new, model_test_set_dfn, model_test_set_opt, model_test_set_order, model_test_set_pvalue, model_test_set_teststat, model_test_set_value, set_model_id;
import gretl.gretl_utilsh : floatneq;
import gretl.dataseth : dataset_is_cross_section, dataset_is_time_series, sample_size;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.libgretlh;
import gretl.libset : libset_get_bool, libset_get_int, libset_set_int;
import gretl.libseth : HACopts, SetKey;
import gretl.missingh;
import gretl.missing_privateh;
import gretl.pvalues : chisq_cdf_comp, normal_cdf;
import gretl.graphingh : GptFlags, PlotType;
import gretl.compareh : ModelAuxCode;
import gretl.gretl_modelh : ModelTestType;
import gretl.gretl_list : gretl_list_append_term, gretl_list_const_pos, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_has_separator, gretl_list_new, gretl_list_separator_position, gretl_list_split_on_separator, gretl_list_sublist, gretl_null_list, list_members_replaced, reglist_check_for_const;
import gretl.strutils : gretl_strdup, strings_array_free, strings_array_new_with_length;
import gretl.options : incompatible_options, option_prereq_missing, options_incompatible_with;
import gretl.gretl_errors : gretl_errmsg_append, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear;
import gretl.missing : apply_reference_missmask, list_adjust_sample, model_adjust_sample, model_missing, reference_missmask_present;
import gretl.transforms : loggenr, xpxgenr;
import gretl.gretl_utils : gretl_isdummy, gretl_ispositive, ijton, ls_criteria, record_test_result;
import gretl.qr_estimate : gretl_qr_regress, lapack_cholesky_regress;
import gretl.gretl_paths : gretl_x12_arima;
import gretl.plugins : get_plugin_function;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.describe : gretl_corr, gretl_corr_rsq, gretl_mean, gretl_stddev;
import gretl.gretl_panel : panel_wls_by_unit, real_panel_model;
import gretl.modelprint : print_coeffs, printmodel;
import gretl.printout : bufspace;
import gretl.tsls : tsls;
import gretl.nls : ivreg_via_gmm;
import gretl.system : single_equation_liml;
import gretl.graphing : finalize_plot_input_file, open_plot_input_file;

/* plugin entry points, fetched at run time with get_plugin_function()
   (the plugins are still C, hence extern (C)) */

private alias LadDriver = extern (C) int function(Model *, DatasetRaw *, Gretlopt);
private alias RqDriver = extern (C) int function(const(GretlMatrixRaw) *, Model *,
                                                 DatasetRaw *, Gretlopt, Prn *);
private alias ArmaModelFn = extern (C) Model function(const(int) *, const(int) *,
                                                      DatasetRaw *, Gretlopt, Prn *);
private alias ArmaX12ModelFn = extern (C) Model function(const(int) *, const(int) *,
                                                         const(DatasetRaw) *, int,
                                                         Gretlopt, Prn *);
private alias ArmaSelectFn = extern (C) int function(const(int) *, const(int) *,
                                                     DatasetRaw *, Gretlopt, Prn *);
private alias GarchModelFn = extern (C) Model function(const(int) *, DatasetRaw *,
                                                       Prn *, Gretlopt);
private alias MplsqFn = extern (C) int function(const(int) *, const(int) *,
                                                const(int) *, DatasetRaw *, Model *,
                                                Gretlopt);
private alias DpdEstimateFn = extern (C) Model function(const(int) *, const(int) *,
                                                        const(char) *, const(DatasetRaw) *,
                                                        Gretlopt, Prn *);
private alias GretlAnovaFn = extern (C) int function(const(int) *, const(DatasetRaw) *,
                                                     Gretlopt, Prn *);

/* ---- declarations of not-yet-ported libgretl items used below ----
   Replace each with an import when its defining file is ported.
   The prototypes are taken from the libgretl headers (missing_private.h,
   gretl_model.h, gretl_list.h, describe.h, dataset.h, options.h,
   qr_estimate.h, tsls.h, system.h, gretl_panel.h, estim_private.h, ...). */

extern (C) {
}

/* ---- end of declarations ---- */

private enum TINY = 8.0e-09; /* was 2.1e-09 (last changed 2007/01/20) */
private enum SMALL = 2.0e-08; /* threshold for printing a warning for collinearity */
private enum YBARZERO = 0.5e-14; /* threshold for treating mean of dependent
                                    variable as effectively zero */
private enum ESSZERO = 1.0e-22; /* threshold for considering a tiny error-sum-of-
                                   squares value to be effectively zero */

/* compute statistics for the dependent variable in a model */

private void model_depvar_stats(Model * pmod, DatasetRaw * dset)
{
    double xx, sum = 0.0;
    int yno = pmod.list[1];
    int t, zw = 0;

    if (pmod.ci == GretlCmdIndex.wls && gretl_model_get_int(pmod, "wt_zeros")) {
        /* WLS with some zero weights */
        zw = 1;
    }

    pmod.ybar = pmod.sdy = NADBL;

    if (pmod.nobs <= 0) {
        return;
    }

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (zw && dset.Z[pmod.nwt][t] == 0.0) {
            continue;
        }
        if (!model_missing(pmod, t)) {
            sum += dset.Z[yno][t];
        }
    }

    pmod.ybar = sum / pmod.nobs;

    sum = 0.0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (zw && dset.Z[pmod.nwt][t] == 0.0) {
            continue;
        }
        if (!model_missing(pmod, t)) {
            sum += (dset.Z[yno][t] - pmod.ybar);
        }
    }

    pmod.ybar = pmod.ybar + sum / pmod.nobs;

    if (fabs(pmod.ybar) < YBARZERO) {
        pmod.ybar = 0.0;
    }

    sum = 0.0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (zw && dset.Z[pmod.nwt][t] == 0.0) {
            continue;
        }
        if (!model_missing(pmod, t)) {
            xx = dset.Z[yno][t] - pmod.ybar;
            sum += xx * xx;
        }
    }

    sum = (pmod.nobs > 1)? sum / (pmod.nobs - 1) : 0.0;

    pmod.sdy = (sum >= 0)? sqrt(sum) : NADBL;
}

/* determine the degrees of freedom for a model */

private int get_model_df(Model * pmod)
{
    int err = 0;

    pmod.ncoeff = pmod.list[0] - 1;
    pmod.dfd = pmod.nobs - pmod.ncoeff;

    if (pmod.dfd < 0) {
        pmod.errcode = GretlError.df;
        gretl_errmsg_sprintf(_("No. of obs (%d) is less than no. "
                               ~ "of parameters (%d)").ptr, pmod.nobs,
                             pmod.ncoeff);
        err = 1;
    } else {
        pmod.dfn = pmod.ncoeff - pmod.ifc;
    }

    return err;
}

private int transcribe_ld_vcv(Model * targ, Model * src)
{
    int nv = targ.ncoeff;
    int nxpx = (nv * nv + nv) / 2;
    int i, j, err = 0;

    err = makevcv(src, src.sigma);

    if (!err && targ.vcv == null) {
        targ.vcv = cast(double *) malloc(nxpx * double.sizeof);
        if (targ.vcv == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        for (i=0; i<nv; i++) {
            for (j=i; j<nv; j++) {
                targ.vcv[ijton(i, j, nv)] =
                    src.vcv[ijton(i, j, src.ncoeff)];
            }
        }
    }

    return err;
}

/* Calculate consistent standard errors (and VCV matrix) when doing
   AR(1) estimation of a model with lagged dependent variable.  See
   Ramanathan, Introductory Econometrics, 5e, p. 450.
*/

private int ldepvar_std_errors(Model * pmod, DatasetRaw * dset)
{
    Model emod;
    const(double) * x;
    int orig_t1 = dset.t1;
    int orig_t2 = dset.t2;
    double rho = gretl_model_get_double(pmod, "rho_gls");
    int origv = dset.v;
    int vnew = pmod.list[0] + 1 - pmod.ifc;
    int * list;
    int vi, vm;
    int i, t;
    int err = 0;

    list = gretl_list_new(vnew + pmod.ifc);
    if (list == null) {
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    err = dataset_add_series(dset, vnew);
    if (err) {
        free(list);
        pmod.errcode = GretlError.alloc;
        return 1;
    }

    vi = origv;

    /* dependent var is residual from original model */
    for (t=0; t<dset.n; t++) {
        dset.Z[vi][t] = pmod.uhat[t];
    }
    strcpy(dset.varname[vi], "eps");
    list[1] = vi++;

    /* indep vars are rho-differenced vars from original model */
    for (i=2; i<=pmod.list[0]; i++) {
        vm = pmod.list[i];
        if (vm == 0) {
            list[i] = 0;
            continue;
        }
        sprintf(dset.varname[vi], "%.6s_r", dset.varname[vm]);
        x = dset.Z[vm];
        for (t=0; t<dset.n; t++) {
            if (t == 0 || na(x[t]) || na(x[t-1])) {
                dset.Z[vi][t] = NADBL;
            } else {
                dset.Z[vi][t] = x[t] - rho * x[t-1];
            }
        }
        list[i] = vi++;
    }

    /* last indep var is lagged u-hat */
    for (t=0; t<dset.n; t++) {
        if (t == 0) {
            dset.Z[vi][t] = NADBL;
        } else {
            dset.Z[vi][t] = dset.Z[pmod.list[1]][t-1];
        }
        if (na(dset.Z[vi][t])) {
            continue;
        }
        for (i=0; i<pmod.ncoeff; i++) {
            x = dset.Z[pmod.list[i+2]];
            if (na(x[t-1])) {
                dset.Z[vi][t] = NADBL;
                break;
            } else {
                dset.Z[vi][t] -= pmod.coeff[i] * x[t-1];
            }
        }
    }

    list[list[0]] = vi;
    strcpy(dset.varname[vi], "uhat_1");

    dset.t1 = pmod.t1;
    dset.t2 = pmod.t2;

    emod = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (emod.errcode) {
        err = emod.errcode;
    } else {
        for (i=0; i<pmod.ncoeff; i++) {
            pmod.sderr[i] = emod.sderr[i];
        }

        err = transcribe_ld_vcv(pmod, &emod);
    }

    clear_model(&emod);

    free(list);
    dataset_drop_last_variables(dset, vnew);

    dset.t1 = orig_t1;
    dset.t2 = orig_t2;

    if (err && !pmod.errcode) {
        pmod.errcode = err;
    }

    return err;
}

/* special computation of statistics for autoregressive models */

private int compute_ar_stats(Model * pmod, const(DatasetRaw) * dset,
                             double rho)
{
    int i, v, t, yno = pmod.list[1];
    int pwe = (pmod.opt & Gretlopt.OPT_P) != 0;
    double x, pw1 = 0.0;

    if (pwe) {
        pw1 = sqrt(1.0 - rho * rho);
    }

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (t == pmod.t1 && pwe) {
            x = pw1 * dset.Z[yno][t];
            for (i=pmod.ifc; i<pmod.ncoeff; i++) {
                v = pmod.list[i+2];
                x -= pmod.coeff[i] * pw1 * dset.Z[v][t];
            }
            if (pmod.ifc) {
                x -= pw1 * pmod.coeff[0];
            }
        } else {
            x = dset.Z[yno][t] - rho*dset.Z[yno][t-1];
            for (i=0; i<pmod.ncoeff; i++) {
                v = pmod.list[i+2];
                x -= pmod.coeff[i] * (dset.Z[v][t] - rho*dset.Z[v][t-1]);
            }
        }
        pmod.uhat[t] = x;
        pmod.yhat[t] = dset.Z[yno][t] - x;
    }

    pmod.rsq = gretl_corr_rsq(pmod.t1, pmod.t2, dset.Z[yno], pmod.yhat);
    pmod.adjrsq =
        1.0 - ((1.0 - pmod.rsq) * (pmod.t2 - pmod.t1) /
               cast(double) pmod.dfd);

    return 0;
}

/* calculation of WLS stats (in agreement with GNU R) */

private void get_wls_stats(Model * pmod, const(DatasetRaw) * dset,
                           Gretlopt opt)
{
    const(double) * y = dset.Z[pmod.list[1]];
    const(double) * w = dset.Z[pmod.nwt];
    double wmean = 0.0, wsum = 0.0;
    double dy, tss = 0;
    int t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!model_missing(pmod, t) && w[t] > 0) {
            wmean += w[t] * y[t];
            wsum += w[t];
        }
    }

    wmean /= wsum;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!model_missing(pmod, t) && w[t] > 0) {
            dy = y[t] - wmean;
            tss += w[t] * dy * dy;
        }
    }

    if (!(opt & Gretlopt.OPT_R)) {
        pmod.fstt = ((tss - pmod.ess) * pmod.dfd) / (pmod.dfn * pmod.ess);
    }
    pmod.rsq = (1 - (pmod.ess / tss));
    pmod.adjrsq = 1 - ((1 - pmod.rsq) * (pmod.nobs - 1)/pmod.dfd);
}

private void fix_wls_values(Model * pmod, const(DatasetRaw) * dset)
{
    int dwt = gretl_model_get_int(pmod, "wt_dummy");
    double yh, ess_orig = 0.0;
    double sw, sigma_orig;
    int t, i, v;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        if (dset.Z[pmod.nwt][t] == 0.0) {
            yh = 0.0;
            for (i=0; i<pmod.ncoeff; i++) {
                v = pmod.list[i+2];
                yh += pmod.coeff[i] * dset.Z[v][t];
            }
            pmod.yhat[t] = yh;
            v = pmod.list[1];
            pmod.uhat[t] = dset.Z[v][t] - yh;
        } else if (!dwt) {
            sw = sqrt(dset.Z[pmod.nwt][t]);
            pmod.yhat[t] /= sw;
            pmod.uhat[t] /= sw;
            ess_orig += pmod.uhat[t] * pmod.uhat[t];
        }
    }

    if (!dwt) {
        sigma_orig = sqrt(ess_orig / pmod.dfd);
        gretl_model_set_double(pmod, "ess_orig", ess_orig);
        gretl_model_set_double(pmod, "sigma_orig", sigma_orig);
    }
}

/* drop the weight var from the list of regressors (WLS) */

private void dropwt(Model * pmod)
{
    int i;

    pmod.list[0] -= 1;
    for (i=1; i<=pmod.list[0]; i++) {
        pmod.list[i] = pmod.list[i+1];
    }
}

private int model_missval_count(const(Model) * pmod)
{
    int mc = 0;

    if (pmod.missmask != null) {
        int t;

        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (pmod.missmask[t] == '1') {
                mc++;
            }
        }
    }

    return mc;
}

private int wls_usable_obs(Model * pmod, const(DatasetRaw) * dset)
{
    int t, n = pmod.nobs;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (pmod.missmask[t] == '1') {
            n--;
        } else if (dset.Z[pmod.nwt][t] == 0) {
            n--;
        }
    }

    return n;
}

private int lsq_check_for_missing_obs(Model * pmod, Gretlopt opt, DatasetRaw * dset,
                                      int * misst)
{
    int reject_missing = 0;
    int missv = 0;

    if (reference_missmask_present()) {
        int err = apply_reference_missmask(pmod);

        /* If there was a reference mask present, it was put there
           as part of a hypothesis test on some original model, and
           it should be respected in estimation of this model
        */
        if (err) {
            pmod.errcode = GretlError.alloc;
            return 1;
        } else {
            return 0;
        }
    }

    if (opt & Gretlopt.OPT_M) {
        reject_missing = 1;
    } else if (libset_get_int(SetKey.hac_missvals) == HACopts.refuse) {
        /* we won't do HAC VCV with embedded missing obs */
        if ((opt & Gretlopt.OPT_R) && dataset_is_time_series(dset) &&
            !libset_get_bool(SetKey.force_hc)) {
            reject_missing = 2;
        }
    }

    if (reject_missing) {
        /* reject missing obs within adjusted sample */
        missv = model_adjust_sample(pmod, dset, misst);
        if (reject_missing == 2) {
            gretl_errmsg_append(_("HAC standard errors not available").ptr, 0);
        }
    } else {
        /* we'll try to work around missing obs */
        missv = model_adjust_sample(pmod, dset, null);
    }

    return missv;
}

private int check_for_lags(Model * pmod, const(DatasetRaw) * dset)
{
    int pmax = 0;
    int ldv = model_lags(pmod.list, dset, &pmax);

    if (pmax > 0) {
        gretl_model_set_int(pmod, "maxlag", pmax);
    } else if (gretl_model_get_int(pmod, "maxlag")) {
        gretl_model_destroy_data_item(pmod, "maxlag");
    }

    if (ldv) {
        gretl_model_set_int(pmod, "ldepvar", ldv);
    } else if (gretl_model_get_int(pmod, "ldepvar")) {
        gretl_model_destroy_data_item(pmod, "ldepvar");
    }

    return ldv;
}

private void log_depvar_ll(Model * pmod, const(DatasetRaw) * dset)
{
    char[VNAMELEN] parent;

    if (series_is_log(dset, pmod.list[1], parent.ptr)) {
        double jll = pmod.lnL;
        int t;

        for (t=0; t<dset.n; t++) {
            if (!na(pmod.uhat[t])) {
                jll -= dset.Z[pmod.list[1]][t];
            }
        }
        gretl_model_set_double(pmod, "jll", jll);
        gretl_model_set_string_as_data(pmod,
                                       "log-parent",
                                       gretl_strdup(parent.ptr));
    }
}

/* check_weight_var(): Here we're screening out certain errors to do
   with the weight series for WLS, and also checking for 0/1 dummy
   weights.

   * Negative weights are not permitted: we flag an error if any are
     found.
   * 0/1 weights are OK, but otherwise zero weights are not accepted
     unless the --allow-zeros option (OPT_Z) is passed to the "wls"
     command.
*/

private int check_weight_var(Model * pmod, const(double) * w, int * effobs,
                             Gretlopt opt)
{
    string all_zeros =
        N_("Weight variable is all zeros, aborting regression");
    string some_neg =
        N_("Weight variable contains negative values");
    string bad_zeros =
        N_("Weight variable is not a dummy but contains zeros");
    int zeros = 0, nobs = 0;
    int t, is_dummy = 1;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (w[t] < 0.0) {
            gretl_errmsg_set(_(some_neg).ptr);
            pmod.errcode = GretlError.data;
            return 1;
        } else if (w[t] > 0.0) {
            nobs++;
            if (w[t] != 1.0) {
                is_dummy = 0;
            }
        } else if (w[t] == 0.0) {
            zeros++;
        }
    }

    if (nobs == 0) {
        gretl_errmsg_set(_(all_zeros).ptr);
        pmod.errcode = GretlError.data;
        return 1;
    }

    /* note: this @effobs value does not take into account missing
       values of dep var or regressors, and so is potentially
       misleading! (FIXME)
    */
    *effobs = nobs;

    if (is_dummy) {
        /* the weight var is a 0/1 dummy */
        gretl_model_set_int(pmod, "wt_dummy", 1);
        gretl_model_set_int(pmod, "wt_zeros", zeros);
    } else if (zeros > 0) {
        if (opt & Gretlopt.OPT_Z) {
            gretl_model_set_int(pmod, "wt_zeros", zeros);
        } else {
            gretl_errmsg_set(_(bad_zeros).ptr);
            pmod.errcode = GretlError.data;
            return 1;
        }
    }

    return 0;
}

void maybe_shift_ldepvar(Model * pmod, DatasetRaw * dset)
{
    if (gretl_model_get_int(pmod, "ldepvar")) {
        check_for_lags(pmod, dset);
    }
}

/*
 * XTX_XTy:
 * @list: list of variables in model.
 * @t1: starting observation.
 * @t2: ending observation.
 * @dset: pointer to dataset.
 * @nwt: ID number of variable used as weight, or 0.
 * @rho: quasi-differencing coefficent, or 0.0;
 * @pwe: if non-zero, use Prais-Winsten for first observation.
 * @xpx: on output, X'X matrix as lower triangle, stacked by columns.
 * @xpy: on output, X'y vector (but see below).
 * @ysum: location to recieve \sum y_i, or NULL.
 * @ypy: location to recieve (scalar) y'y, or NULL.
 * @mask: missing observations mask, or NULL.
 *
 * Calculates X'X and X'y, with various possible transformations
 * of the original data, depending on @nwt, @rho and @pwe.
 * If X'y is not required, @xpy can be given as NULL.
 *
 * Note: the y-related pointer arguments @xpy, @ysum, and @ypy
 * form a "package": either all should be given, or all should
 * be NULL.
 *
 * Returns: 0 on success, non-zero on error.
 */

private int XTX_XTy(const(int) * list, int t1, int t2,
                    const(DatasetRaw) * dset, int nwt,
                    double rho, int pwe,
                    double * xpx, double * xpy,
                    double * ysum, double * ypy,
                    const(char) * mask)
{
    int yno = list[1];
    int lmin = (xpy != null)? 2 : 1;
    int lmax = list[0];
    int qdiff = (rho != 0.0);
    const(double) * y = null;
    const(double) * w = null;
    const(double) * xi = null;
    const(double) * xj = null;
    double x, pw1;
    int i, j, t, m;
    int err = 0;

    /* Prais-Winsten term */
    if (qdiff && pwe) {
        pw1 = sqrt(1.0 - rho * rho);
    } else {
        pwe = 0;
        pw1 = 0.0;
    }

    /* dependent variable */
    y = dset.Z[yno];

    if (nwt) {
        /* weight variable */
        w = dset.Z[nwt];
    }

    if (xpy != null) {
        *ysum = *ypy = 0.0;

        for (t=t1; t<=t2; t++) {
            if (masked(mask, t)) {
                continue;
            }
            x = y[t];
            if (qdiff) {
                if (pwe && t == t1) {
                    x = pw1 * y[t];
                } else {
                    x -= rho * y[t-1];
                }
            } else if (nwt) {
                x *= sqrt(w[t]);
            }
            *ysum += x;
            *ypy += x * x;
        }

        if (*ypy <= 0.0) {
            /* error condition */
            return yno;
        }
    }

    m = 0;

    if (qdiff) {
        /* quasi-difference the data */
        for (i=lmin; i<=lmax; i++) {
            xi = dset.Z[list[i]];
            for (j=i; j<=lmax; j++) {
                xj = dset.Z[list[j]];
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (pwe && t == t1) {
                        x += pw1 * xi[t1] * pw1 * xj[t];
                    } else {
                        x += (xi[t] - rho * xi[t-1]) *
                            (xj[t] - rho * xj[t-1]);
                    }
                }
                if (i == j && x < DBL_EPSILON)  {
                    return GretlError.singular;
                }
                xpx[m++] = x;
            }
            if (xpy != null) {
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (pwe && t == t1) {
                        x += pw1 * y[t] * pw1 * xi[t];
                    } else {
                        x += (y[t] - rho * y[t-1]) *
                            (xi[t] - rho * xi[t-1]);
                    }
                }
                xpy[i-2] = x;
            }
        }
    } else if (nwt) {
        /* weight the data */
        for (i=lmin; i<=lmax; i++) {
            xi = dset.Z[list[i]];
            for (j=i; j<=lmax; j++) {
                xj = dset.Z[list[j]];
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (!masked(mask, t)) {
                        x += w[t] * xi[t] * xj[t];
                    }
                }
                if (i == j && x < DBL_EPSILON) {
                    return GretlError.singular;
                }
                xpx[m++] = x;
            }
            if (xpy != null) {
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (!masked(mask, t)) {
                        x += w[t] * y[t] * xi[t];
                    }
                }
                xpy[i-2] = x;
            }
        }
    } else if (mask != null) {
        /* plain data, but with missing obs mask */
        for (i=lmin; i<=lmax; i++) {
            xi = dset.Z[list[i]];
            for (j=i; j<=lmax; j++) {
                xj = dset.Z[list[j]];
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (!masked(mask, t)) {
                        x += xi[t] * xj[t];
                    }
                }
                if (i == j && x < DBL_EPSILON) {
                    return GretlError.singular;
                }
                xpx[m++] = x;
            }
            if (xpy != null) {
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    if (!masked(mask, t)) {
                        x += y[t] * xi[t];
                    }
                }
                xpy[i-2] = x;
            }
        }
    } else {
        /* plain data, no missing obs mask */
        for (i=lmin; i<=lmax; i++) {
            xi = dset.Z[list[i]];
            for (j=i; j<=lmax; j++) {
                xj = dset.Z[list[j]];
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    x += xi[t] * xj[t];
                }
                if (i == j && x < DBL_EPSILON) {
                    return GretlError.singular;
                }
                xpx[m++] = x;
            }
            if (xpy != null) {
                x = 0.0;
                for (t=t1; t<=t2; t++) {
                    x += y[t] * xi[t];
                }
                xpy[i-2] = x;
            }
        }
    }

    return err;
}

private int native_cholesky_regress(Model * pmod,
                                    const(DatasetRaw) * dset,
                                    Gretlopt opt)
{
    double rho, ysum = 0.0, ypy = 0.0;
    double * xpy;
    int k = pmod.ncoeff;
    int nxpx = k * (k + 1) / 2;
    int pwe = (pmod.opt & Gretlopt.OPT_P) != 0;
    int i;

    if (nxpx == 0) {
        fprintf(stderr, "problem: nxpx = 0\n");
        pmod.errcode = GretlError.data;
        return pmod.errcode;
    }

    xpy = cast(double *) malloc(k * double.sizeof);
    if (xpy == null) {
        pmod.errcode = GretlError.alloc;
        return pmod.errcode;
    }

    pmod.xpx = cast(double *) malloc(nxpx * double.sizeof);
    pmod.coeff = cast(double *) malloc(k * double.sizeof);
    pmod.sderr = cast(double *) malloc(k * double.sizeof);

    if (pmod.yhat == null) {
        pmod.yhat = cast(double *) malloc(pmod.full_n * double.sizeof);
    }

    if (pmod.uhat == null) {
        pmod.uhat = cast(double *) malloc(pmod.full_n * double.sizeof);
    }

    if (pmod.xpx == null || pmod.coeff == null ||
        pmod.sderr == null || pmod.yhat == null ||
        pmod.uhat == null) {
        free(xpy);
        pmod.errcode = GretlError.alloc;
        return pmod.errcode;
    }

    for (i=0; i<k; i++) {
        xpy[i] = 0.0;
    }
    for (i=0; i<nxpx; i++) {
        pmod.xpx[i] = 0.0;
    }

    rho = gretl_model_get_double_default(pmod, "rho_gls", 0.0);

    /* calculate regression results, Cholesky style */
    pmod.errcode = XTX_XTy(pmod.list, pmod.t1, pmod.t2, dset,
                           pmod.nwt, rho, pwe, pmod.xpx, xpy,
                           &ysum, &ypy, pmod.missmask);

    if (!pmod.errcode) {
        regress(pmod, xpy, ysum, ypy, dset, rho, opt);
    }

    free(xpy);

    return pmod.errcode;
}

private int gretl_null_regress(Model * pmod, const(DatasetRaw) * dset)
{
    double yt;
    int t;

    if (pmod.yhat == null) {
        pmod.yhat = cast(double *) malloc(pmod.full_n * double.sizeof);
    }

    if (pmod.uhat == null) {
        pmod.uhat = cast(double *) malloc(pmod.full_n * double.sizeof);
    }

    if (pmod.yhat == null || pmod.uhat == null) {
        pmod.errcode = GretlError.alloc;
        return pmod.errcode;
    }

    pmod.nobs = 0;
    pmod.ifc = 0;
    pmod.ess = 0.0;
    pmod.rsq = pmod.adjrsq = 0.0;

    for (t=0; t<pmod.full_n; t++) {
        yt = dset.Z[pmod.list[1]][t];
        if (t < pmod.t1 || t > pmod.t2 || na(yt)) {
            pmod.uhat[t] = pmod.yhat[t] = NADBL;
        } else {
            pmod.uhat[t] = yt;
            pmod.yhat[t] = 0.0;
            pmod.ess += yt * yt;
            pmod.nobs += 1;
        }
    }

    if (pmod.ess == 0) {
        pmod.sigma = 0.0;
    } else if (pmod.nobs > 1) {
        pmod.sigma = sqrt(pmod.ess / (pmod.nobs - 1));
    } else {
        pmod.errcode = GretlError.data;
    }

    if (pmod.errcode == 0) {
        ls_criteria(pmod);
    }

    return pmod.errcode;
}

/* check whether the variable with ID @yno is all zeros;
   return 1 if so, otherwise return 0 */

private int depvar_zero(const(Model) * pmod, int yno, const(double *) * Z)
{
    double yt;
    int t, ret = 1;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        yt = Z[yno][t];
        if (pmod.nwt) {
            yt *= Z[pmod.nwt][t];
        }
        if (yt != 0.0) {
            ret = 0;
            break;
        }
    }

    return ret;
}

private int cholesky_regress(Model * pmod, const(DatasetRaw) * dset,
                             Gretlopt opt)
{
    int T = pmod.t2 - pmod.t1 + 1;
    int k = pmod.list[0] - 1;

    if (k >= 50 || (T >= 250 && k >= 30)) {
        return lapack_cholesky_regress(pmod, dset, opt);
    } else {
        return native_cholesky_regress(pmod, dset, opt);
    }
}

/* limited freeing of elements before passing a model
   on for QR estimation in the case of (near-)singularity
*/

private void model_free_storage(Model * pmod)
{
    free(pmod.xpx);
    free(pmod.coeff);
    free(pmod.sderr);

    pmod.xpx = null;
    pmod.coeff = null;
    pmod.sderr = null;
}

/*
 * ar1_lsq:
 * @list: model specification.
 * @dset: dataset struct.
 * @ci: command index, e.g. OLS, AR1.
 * @opt: option flags (see lsq and ar1_model).
 * @rho: coefficient for quasi-differencing the data, or 0.
 *
 * Nests lsq() and ar1_model(); is also used internally with ci == OLS
 * but rho != 0.0 when estimating rho via, e.g., Cochrane-Orcutt
 * iteration.
 *
 * Returns: model struct containing the estimates.
 */

private Model ar1_lsq(const(int) * list, DatasetRaw * dset,
                      GretlCmdIndex ci, Gretlopt opt,
                      double rho)
{
    Model mdl;
    int wtdobs = 0;
    int missv = 0, misst = 0;
    int pwe = (opt & Gretlopt.OPT_P) != 0;
    int nullmod = 0, ldv = 0;
    int save_hc = -1;
    int yno, i;

    gretl_model_init(&mdl, dset);

    if (list == null || dset == null || dset.Z == null) {
        fprintf(stderr, "E_DATA: lsq: list = %p, dset = %p\n",
                cast(void *) list, cast(void *) dset);
        mdl.errcode = GretlError.data;
        return mdl;
    }

    if (ci == GretlCmdIndex.ols &&
        incompatible_options(opt, Gretlopt.OPT_C | Gretlopt.OPT_J)) {
        /* at present --cluster and --jackknife can't be combined */
        mdl.errcode = GretlError.badopt;
        return mdl;
    }

    if (opt & Gretlopt.OPT_C) {
        /* cluster option implies robust */
        opt |= Gretlopt.OPT_R;
    }

    if (ci == GretlCmdIndex.ar1) {
        if (opt & Gretlopt.OPT_P) {
            /* OPT_P: Prais-Winsten */
            mdl.opt |= Gretlopt.OPT_P;
        } else if (opt & Gretlopt.OPT_H) {
            /* OPT_H: Hildreth-Lu */
            mdl.opt |= Gretlopt.OPT_H;
        }
    } else if (rho != 0.0 && (opt & Gretlopt.OPT_P)) {
        /* estimation of rho in progress */
        mdl.opt |= Gretlopt.OPT_P;
    }

    if (list[0] == 1 && ci == GretlCmdIndex.ols && (opt & Gretlopt.OPT_U)) {
        /* null model OK */
        nullmod = 1;
    } else if (list[0] == 1 || dset.v == 1) {
        fprintf(stderr, "E_DATA: lsq: list[0] = %d, dset->v = %d\n",
                list[0], dset.v);
        mdl.errcode = GretlError.data;
        return mdl;
    }

    /* preserve a copy of the list supplied, for future reference */
    mdl.list = gretl_list_copy(list);
    if (mdl.list == null) {
        mdl.errcode = GretlError.alloc;
        return mdl;
    }

    mdl.t1 = dset.t1;
    mdl.t2 = dset.t2;
    mdl.full_n = dset.n;
    mdl.ci = ci;

    /* Doing weighted least squares? */
    if (ci == GretlCmdIndex.wls) {
        check_weight_var(&mdl, dset.Z[mdl.list[1]], &wtdobs, opt);
        if (mdl.errcode) {
            return mdl;
        }
        mdl.nwt = mdl.list[1];
    } else {
        mdl.nwt = 0;
    }

    /* sanity check */
    if (mdl.t1 < 0 || mdl.t2 > dset.n - 1) {
        mdl.errcode = GretlError.data;
        goto lsq_abort;
    }

    /* adjust sample range and check for missing obs: this
       may set the model errcode */
    missv = lsq_check_for_missing_obs(&mdl, opt, dset, &misst);
    if (mdl.errcode) {
        goto lsq_abort;
    }

    /* react to presence of unhandled missing obs */
    if (missv) {
        gretl_errmsg_sprintf(_("Missing value encountered for "
                               ~ "variable %d, obs %d").ptr,
                             missv, misst);
        mdl.errcode = GretlError.data;
        return mdl;
    }

    if (ci == GretlCmdIndex.wls) {
        /* note: remove the WLS weight series from regression list */
        dropwt(&mdl);
    }

    yno = mdl.list[1];

    /* check for unknown vars in list */
    for (i=1; i<=mdl.list[0]; i++) {
        if (mdl.list[i] > dset.v - 1) {
            mdl.errcode = GretlError.unkvar;
            goto lsq_abort;
        }
    }

    /* check for zero dependent var */
    if (depvar_zero(&mdl, yno, cast(const(double *) *) dset.Z)) {
        mdl.errcode = GretlError.zero;
        goto lsq_abort;
    }

    /* drop any vars that are all zero and repack the list */
    omitzero(&mdl, dset, opt);

    /* if regressor list contains a constant, record this fact and
       place it first among the regressors */
    mdl.ifc = reglist_check_for_const(mdl.list, dset);

    /* Check for presence of lagged dependent variable?
       (Don't bother if this is an auxiliary regression.)
    */
    if (!(opt & Gretlopt.OPT_A) && !dataset_is_cross_section(dset)) {
        ldv = check_for_lags(&mdl, dset);
    }

    /* AR1: advance the starting observation by one? */
    if (rho != 0.0 && !pwe) {
        mdl.t1 += 1;
    }

    mdl.ncoeff = mdl.list[0] - 1;
    if (wtdobs > 0 && wtdobs <= mdl.t2 - mdl.t1 + 1 && mdl.missmask == null) {
        /* FIXME: this is WLS-specific, but is it right? */
        mdl.nobs = wtdobs;
    } else {
        mdl.nobs = mdl.t2 - mdl.t1 + 1;
        if (mdl.missmask != null) {
            if (mdl.nwt) {
                mdl.nobs = wls_usable_obs(&mdl, dset);
            } else {
                mdl.nobs -= model_missval_count(&mdl);
            }
        }
    }

    /* check degrees of freedom */
    if (get_model_df(&mdl)) {
        goto lsq_abort;
    }

    /* if df correction is not wanted, record this fact */
    if (opt & Gretlopt.OPT_N) {
        mdl.opt |= Gretlopt.OPT_N;
    }

    if (dataset_is_time_series(dset)) {
        opt |= Gretlopt.OPT_T;
    }

    /* remove Durbin-Watson p-value flag if not appropriate */
    if (ldv || !(opt & Gretlopt.OPT_T)) {
        opt &= ~Gretlopt.OPT_I;
    }

    if (opt & Gretlopt.OPT_J) {
        /* --jackknife */
        save_hc = libset_get_int(SetKey.hc_version);
        libset_set_int(SetKey.hc_version, 4);
        opt |= Gretlopt.OPT_R;
        mdl.opt |= (Gretlopt.OPT_J | Gretlopt.OPT_R);
    }

    if (rho != 0.0) {
        gretl_model_set_double(&mdl, "rho_gls", rho);
    }

    if (nullmod) {
        gretl_null_regress(&mdl, dset);
    } else if (libset_get_int(SetKey.use_qr) > 0) {
        gretl_qr_regress(&mdl, dset, opt);
    } else if (opt & (Gretlopt.OPT_R | Gretlopt.OPT_I)) {
        gretl_qr_regress(&mdl, dset, opt);
    } else {
        cholesky_regress(&mdl, dset, opt);
        if (mdl.errcode == GretlError.singular) {
            /* near-perfect collinearity is better handled by QR */
            model_free_storage(&mdl);
            gretl_qr_regress(&mdl, dset, opt);
        }
    }

    if (save_hc >= 0) {
        libset_set_int(SetKey.hc_version, save_hc);
    }

    if (mdl.errcode) {
        goto lsq_abort;
    }

    /* get the mean and sd of dep. var. and make available */
    model_depvar_stats(&mdl, dset);

    /* Doing an autoregressive procedure? */
    if (ci == GretlCmdIndex.ar1) {
        if (compute_ar_stats(&mdl, dset, rho)) {
            goto lsq_abort;
        }
        if (ldv) {
            ldepvar_std_errors(&mdl, dset);
            if (mdl.errcode) {
                goto lsq_abort;
            }
        }
        if ((opt & Gretlopt.OPT_H) && (opt & Gretlopt.OPT_B)) {
            gretl_model_set_int(&mdl, "no-corc", 1);
        }
        if (gretl_model_get_int(&mdl, "maxlag") < 1) {
            gretl_model_set_int(&mdl, "maxlag", 1);
        }
    }

    /* weighted least squares: fix yhat and uhat; add calculation of
       ESS and sigma based on unweighted data
    */
    if (ci == GretlCmdIndex.wls) {
        get_wls_stats(&mdl, dset, opt);
        fix_wls_values(&mdl, dset);
    }

    if (mdl.missmask == null && (opt & Gretlopt.OPT_T) && !(opt & Gretlopt.OPT_I)) {
        mdl.rho = rhohat(1, mdl.t1, mdl.t2, mdl.uhat);
        mdl.dw = dwstat(1, &mdl, dset);
    } else if (!(opt & Gretlopt.OPT_I)) {
        mdl.rho = mdl.dw = NADBL;
    }

    /* special case: degenerate model */
    if (mdl.ncoeff == 1 && mdl.ifc) {
        mdl.rsq = mdl.adjrsq = 0.0;
        mdl.fstt = NADBL;
    }

    /* Generate model selection statistics */
    if (ci == GretlCmdIndex.ar1) {
        mdl.lnL = NADBL;
        mle_criteria(&mdl, 0);
    } else {
        ls_criteria(&mdl);
    }
    if (!(opt & Gretlopt.OPT_A) && !na(mdl.lnL)) {
        log_depvar_ll(&mdl, dset);
    }

    free(mdl.xlist);
    mdl.xlist = gretl_list_sublist(mdl.list, 2, -1);

 lsq_abort:

    if (!(opt & (Gretlopt.OPT_A | Gretlopt.OPT_S))) {
        /* if it's not an auxiliary regression, or part of
           a system, set an ID number on the model
        */
        set_model_id(&mdl, opt);
    }

    return mdl;
}

/**
 * lsq:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @ci: one of the command indices in #LSQ_MODEL.
 * @opt: option flags: zero or more of the following --
 *   OPT_R: compute robust standard errors;
 *   OPT_A: treat as auxiliary regression (don't bother checking
 *     for presence of lagged dependent var, don't augment model count);
 *   OPT_N: don't use degrees of freedom correction for standard
 *      error of regression;
 *   OPT_M: reject missing observations within sample range;
 *   OPT_Z: (internal use) suppress the automatic elimination of
 *      perfectly collinear variables.
 *   OPT_X: compute "variance matrix" as just (X'X)^{-1}
 *   OPT_I: compute Durbin-Watson p-value.
 *   OPT_U: treat null model as OK.
 *   OPT_P: (ar1 only): use Prais-Winsten.
 *   OPT_S: estimation as part of system.
 *
 * Computes least squares estimates of the model specified by @list,
 * using an estimator determined by the value of @ci.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model lsq(const(int) * list, DatasetRaw * dset, GretlCmdIndex ci,
          Gretlopt opt)
{
    return ar1_lsq(list, dset, ci, opt, 0.0);
}

/**
 * ar1_model:
 * @list: model specification.
 * @dset: dataset struct.
 * @opt: option flags: may include OPT_H to use Hildreth-Lu,
 *       OPT_P to use Prais-Winsten, OPT_B to suppress Cochrane-Orcutt
 *       fine-tuning of Hildreth-Lu results, OPT_G to generate
 *       a gnuplot graph of the search in Hildreth-Lu case.
 * @prn: gretl printer.
 *
 * Returns: model struct containing the estimates.
 */

Model ar1_model(const(int) * list, DatasetRaw * dset,
                Gretlopt opt, Prn * prn)
{
    Model mdl;
    double rho;
    int err = 0;

    rho = estimate_rho(list, dset, opt, prn, &err);

    if (err) {
        gretl_model_init(&mdl, null);
        mdl.errcode = err;
        return mdl;
    } else {
        return ar1_lsq(list, dset, GretlCmdIndex.ar1, opt, rho);
    }
}

private int make_ess(Model * pmod, const(DatasetRaw) * dset)
{
    int i, t, yno = pmod.list[1], l0 = pmod.list[0];
    int nwt = pmod.nwt;
    double yhat, ut;

    pmod.ess = 0.0;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (nwt && dset.Z[nwt][t] == 0.0) {
            continue;
        }
        if (model_missing(pmod, t)) {
            continue;
        }
        yhat = 0.0;
        for (i=2; i<=l0; i++) {
            yhat += pmod.coeff[i-2] * dset.Z[pmod.list[i]][t];
        }
        ut = dset.Z[yno][t] - yhat;
        if (nwt) {
            ut *= sqrt(dset.Z[nwt][t]);
        }
        pmod.ess += ut * ut;
    }

    return 0;
}

/* The heuristic used here is that the model effectively
   contains a constant or intercept if the means of y and
   yhat are the same, or in other words the mean residual
   is zero -- but allowing for some numerical slop.
*/

int check_for_effective_const(Model * pmod, const(DatasetRaw) * dset)
{
    double ubar = 0.0, ybar = 0.0;
    const(double) * y = dset.Z[pmod.list[1]];
    int t, ret = 0;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!na(pmod.uhat[t])) {
            ubar += pmod.uhat[t];
            ybar += y[t];
        }
    }

    ubar = fabs(ubar / pmod.nobs);
    ybar = fabs(ybar / pmod.nobs);

    /* the calibration below is debatable: watch out for
       "wrong" results */

    if (ubar < 5.0e-7) {
        /* absolute value of mean residual small enough? */
        ret = 1;
    } else if (ubar < 0.01 && ybar > 1.0e-20) {
        /* try scaled variant? */
        ret = ubar / ybar < 2.0e-15;
    }

    if (ret) {
        gretl_model_set_int(pmod, "effconst", 1);
        pmod.dfn -= 1;
    } else if (gretl_model_get_int(pmod, "effconst")) {
        gretl_model_set_int(pmod, "effconst", 0);
        pmod.dfn += 1;
    }

    return ret;
}

private void uncentered_r_squared(Model * pmod, const(DatasetRaw) * dset)
{
    double y0 = dset.Z[pmod.list[1]][pmod.t1];

    /* special computation for the case of TSS = 0, i.e.
       the dependent variable is a constant */

    if (y0 != 0) {
        double den = pmod.nobs * y0 * y0;

        pmod.rsq = 1 - (pmod.ess / den);
        gretl_model_set_int(pmod, "uncentered", 1);
    }
}

private void compute_r_squared(Model * pmod, const(DatasetRaw) * dset,
                               int * ifc)
{
    int yno = pmod.list[1];

    pmod.rsq = 1.0 - (pmod.ess / pmod.tss);

    if (*ifc == 0) {
        *ifc = check_for_effective_const(pmod, dset);
    }

    if (pmod.dfd > 0) {
        double yt, den = 0.0;

        if (*ifc) {
            den = pmod.tss * pmod.dfd;
            pmod.adjrsq = 1 - (pmod.ess * (pmod.nobs - 1) / den);
        } else {
            /* model does not contain a constant */
            int t;

            for (t=pmod.t1; t<=pmod.t2; t++) {
                if (!na(pmod.yhat[t])) {
                    yt = dset.Z[yno][t];
                    den += yt * yt;
                }
            }

            /* make the centered R^2 available for output */
            gretl_model_set_double(pmod, "centered-R2", pmod.rsq);

            /* but flag the fact that the reported R-squared is
               in fact uncentered */
            gretl_model_set_int(pmod, "uncentered", 1);

            /* and make the "official" figure the uncentered R^2,
               as per NIST, R, Stata, SPSS, ... */
            pmod.rsq = 1 - pmod.ess / den;
            pmod.adjrsq =
                1.0 - ((1.0 - pmod.rsq) * (pmod.nobs - 1.0) / pmod.dfd);
        }
    }

    if (pmod.rsq < 0.0) {
        pmod.rsq = 0.0;
    }
}

/*
  diaginv: solves for the diagonal elements of X'X inverse.

  @xpx = Cholesky-decomposed X'X
  @tmp = work array, length >= nv
  @diag = diagonal elements of {X'X}^{-1} (output)
  @nv = number of regression coefficients = length of diag
*/

private void diaginv(const(double) * xpx, double * tmp, double * diag,
                     int nv)
{
    int kk, l, m, k, i, j;
    double d, e;

    kk = 0;

    for (l=0; l<nv-1; l++) {
        d = xpx[kk];
        tmp[l] = d;
        e = d * d;
        m = 0;
        if (l > 0) {
            for (j=1; j<=l; j++) {
                m += nv - j;
            }
        }
        for (i=l+1; i<nv; i++) {
            d = 0.0;
            k = i + m;
            for (j=l; j<i; j++) {
                d += tmp[j] * xpx[k];
                k += nv - (j+1);
            }
            d = -d * xpx[k];
            tmp[i] = d;
            e += d * d;
        }
        kk += nv - l;
        diag[l] = e;
    }

    diag[nv-1] = xpx[kk] * xpx[kk];
}

/*
  cholbeta: in-place Cholesky decomposition of X'X (pmod->xpx) (lower
  triangular matrix stacked in columns) plus solution for the
  least-squares coefficient estimates.

  pmod->xpx = X'X on input and Cholesky decomposition on output
  xpy = the X'y vector on input and Cholesky-transformed t
        vector on output
  rss = location to receive regression sum of squares

  The number of floating-point operations is basically 3.5 * nv^2
  plus (nv^3) / 3.
*/

private int cholbeta(Model * pmod, double * xpy, double * rss)
{
    int i, j, k, kk, l, jm1;
    double e, d, d1, d2, test, xx;
    double * xpx = pmod.xpx;
    double * b = pmod.coeff;
    int nc = pmod.ncoeff;

    if (xpx[0] <= 0.0) {
        fprintf(stderr, "%s %d: xpx <= 0.0\n", __FILE__.ptr, __LINE__);
        return GretlError.nan;
    }

    e = 1.0 / sqrt(xpx[0]);
    xpx[0] = e;
    xpy[0] *= e;
    for (i=1; i<nc; i++) {
        xpx[i] *= e;
    }

    kk = nc;

    for (j=1; j<nc; j++) {
        /* diagonal elements */
        d = d1 = 0.0;
        k = jm1 = j;

        for (l=1; l<=jm1; l++) {
            xx = xpx[k];
            d1 += xx * xpy[l-1];
            d += xx * xx;
            k += nc - l;
        }

        d2 = xpx[kk] - d;
        test = d2 / xpx[kk];

        /* check for singularity */
        if (test < TINY) {
            *rss = -1.0;
            return GretlError.singular;
        } else if (test < SMALL) {
            gretl_model_set_int(pmod, "near-singular", 1);
        }

        e = 1 / sqrt(d2);
        xpx[kk] = e;
        xpy[j] = (xpy[j] - d1) * e;

        /* off-diagonal elements */
        for (i=j+1; i<nc; i++) {
            kk++;
            d = 0.0;
            k = j;
            for (l=1; l<=jm1; l++) {
                d += xpx[k] * xpx[k-j+i];
                k += nc - l;
            }
            xpx[kk] = (xpx[kk] - d) * e;
        }
        kk++;
    }

    kk--;

    /* calculate regression sum of squares */
    d = 0.0;
    for (j=0; j<nc; j++) {
        d += xpy[j] * xpy[j];
    }
    *rss = d;

    /* back-solve for the coefficients */

    b[nc-1] = xpy[nc-1] * xpx[kk];

    for (j=nc-2; j>=0; j--) {
        d = xpy[j];
        for (i=nc-1; i>j; i--) {
            d -= b[i] * xpx[--kk];
        }
        b[j] = d * xpx[--kk];
    }

    for (j=0; j<nc; j++) {
        if (isnan(b[j])) {
            fprintf(stderr, "%s %d: coeff %d is NaN\n", __FILE__.ptr, __LINE__, j);
            return GretlError.nan;
        }
    }

    return 0;
}

/* compute fitted values and residuals */

private int hatvars(Model * pmod, const(DatasetRaw) * dset)
{
    int yno = pmod.list[1];
    int xno, i, t;
    double x, y;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (model_missing(pmod, t)) {
            continue;
        }
        pmod.yhat[t] = 0.0;
        for (i=0; i<pmod.ncoeff; i++) {
            xno = pmod.list[i+2];
            x = dset.Z[xno][t];
            if (pmod.nwt) {
                x *= sqrt(dset.Z[pmod.nwt][t]);
            }
            pmod.yhat[t] += pmod.coeff[i] * x;
        }
        y = dset.Z[yno][t];
        if (pmod.nwt) {
            y *= sqrt(dset.Z[pmod.nwt][t]);
        }
        pmod.uhat[t] = y - pmod.yhat[t];
    }

    return 0;
}

/*
  regress: takes X'X (pmod->xpx) and X'y (@xpy) and
  computes OLS estimates and associated statistics.

  n = total number of observations per series in data set
  pmod->ifc = 1 if constant is present in model, else = 0

  ess = error sum of squares
  sigma = standard error of regression
  fstt = F-statistic
  pmod->coeff = array of regression coefficients
  pmod->sderr = corresponding array of standard errors
*/

private void regress(Model * pmod, double * xpy,
                     double ysum, double ypy,
                     const(DatasetRaw) * dset,
                     double rho, Gretlopt opt)
{
    int ifc = pmod.ifc;
    double zz, rss = 0.0;
    double s2 = 0.0;
    double * diag = null;
    int i, err = 0;

    for (i=0; i<pmod.full_n; i++) {
        pmod.yhat[i] = pmod.uhat[i] = NADBL;
    }

    zz = ysum * ysum / pmod.nobs;
    pmod.tss = ypy - zz;

    /* Cholesky-decompose X'X and find the coefficients */
    err = cholbeta(pmod, xpy, &rss);
    if (err) {
        pmod.errcode = err;
        return;
    }

    if (rho != 0.0) {
        pmod.ess = ypy - rss;
    } else {
        make_ess(pmod, dset);
        rss = ypy - pmod.ess;
    }

    if (fabs(pmod.ess) < ESSZERO) {
        pmod.ess = 0.0;
    } else if (pmod.ess < 0.0) {
        gretl_errmsg_sprintf(_("Error sum of squares (%g) is not > 0").ptr,
                             pmod.ess);
        return;
    }

    if (pmod.dfd == 0) {
        pmod.sigma = 0.0;
        pmod.adjrsq = NADBL;
    } else {
        if (pmod.opt & Gretlopt.OPT_N) {
            /* no-df-corr */
            s2 = pmod.ess / pmod.nobs;
        } else {
            s2 = pmod.ess / pmod.dfd;
        }
        pmod.sigma = sqrt(s2);
    }

    if (pmod.tss < DBL_EPSILON) {
        pmod.rsq = pmod.adjrsq = NADBL;
    }

    hatvars(pmod, dset);

    if (pmod.errcode) return;

    if (pmod.tss > 0.0) {
        compute_r_squared(pmod, dset, &ifc);
    } else if (pmod.tss == 0.0 && !(opt & Gretlopt.OPT_A)) {
        uncentered_r_squared(pmod, dset);
    }

    if (s2 <= 0.0 || pmod.dfd == 0 || pmod.dfn == 0) {
        pmod.fstt = NADBL;
    } else if (pmod.rsq == 1.0) {
        pmod.fstt = NADBL;
    } else if (pmod.opt & Gretlopt.OPT_N) {
        /* for consistency with no-df-corr */
        pmod.fstt = NADBL;
        pmod.chisq = (rss - zz * ifc) / s2;
    } else {
        pmod.fstt = (rss - zz * ifc) / (s2 * pmod.dfn);
        if (pmod.fstt < 0.0) {
            pmod.fstt = 0.0;
        }
    }

    diag = cast(double *) malloc(pmod.ncoeff * double.sizeof);
    if (diag == null) {
        pmod.errcode = GretlError.alloc;
        return;
    }

    /* Note: 'xpy' is used purely as workspace in diaginv() below, as
       a matter of economy/convenience; no values in this array are
       used before they are written.
    */
    diaginv(pmod.xpx, xpy, diag, pmod.ncoeff);

    for (i=0; i<pmod.ncoeff; i++) {
        if (diag[i] >= 0.0) {
            pmod.sderr[i] = pmod.sigma * sqrt(diag[i]);
        } else {
            pmod.sderr[i] = 0.0;
        }
    }

    free(diag);
}

/**
 * makevcv:
 * @pmod: pointer to model.
 * @sigma: square root of error variance, or 1.0 to
 * produce just X'X^{-1}.
 *
 * Inverts the Cholesky-decomposed X'X and computes the
 * coefficient covariance matrix.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int makevcv(Model * pmod, double sigma)
{
    int dec, mst, kk, i, j, kj, icnt, m, k, l = 0;
    int nv, nxpx;
    double d;

    if (pmod.vcv != null) {
        /* already done */
        return 0;
    }

    if (pmod.xpx == null) {
        /* raw material not available */
        return GretlError.badstat;
    }

    nv = pmod.ncoeff;
    nxpx = (nv * nv + nv) / 2;
    mst = nxpx;
    kk = nxpx - 1;

    pmod.vcv = cast(double *) malloc(nxpx * double.sizeof);
    if (pmod.vcv == null) {
        return GretlError.alloc;
    }

    for (i=0; i<nv; i++) {
        mst -= i;
        /* find diagonal element */
        d = pmod.xpx[kk];
        if (i > 0) {
            for (j=kk+1; j<=kk+i; j++) {
                d -= pmod.xpx[j] * pmod.vcv[j];
            }
        }
        pmod.vcv[kk] = d * pmod.xpx[kk];
        /* find off-diagonal elements indexed by kj */
        kj = kk;
        kk = kk - i - 2;
        if (i > nv - 2) {
            continue;
        }
        for (j=i+1; j<nv; j++) {
            icnt = i+1;
            kj -= j;
            d = 0.0;
            m = mst + 1;
            for (k=0; k<=j-1; k++) {
                if (icnt > 0) {
                    dec = 1;
                    icnt--;
                } else {
                    dec = k;
                }
                m -= dec;
                l = kj + i - k;
                d += pmod.vcv[m-1] * pmod.xpx[l];
            }
            pmod.vcv[kj] = -d * pmod.xpx[l-1];
        }
    }

    if (pmod.ci == GretlCmdIndex.logit || pmod.ci == GretlCmdIndex.probit) {
        sigma = 1.0;
    }

    if (sigma != 1.0) {
        double s2 = sigma * sigma;

        for (i=0; i<nxpx; i++) {
            pmod.vcv[i] *= s2;
        }
    }

    return 0;
}

/**
 * dwstat:
 * @order: order of autoregression (usually 1).
 * @pmod: pointer to model.
 * @dset: pointer to dataset.
 *
 * Computes the Durbin-Watson statistic for @pmod.
 *
 * Returns: the D-W value, or #NADBL on error.
 */

double dwstat(int order, Model * pmod, const(DatasetRaw) * dset)
{
    const(double) * w = null;
    double ut, u1;
    double num = 0.0;
    double den = 0.0;
    int t, t1;

    if (pmod.ess <= 0.0) {
        return NADBL;
    }

    t1 = pmod.t1 + order;

    if (pmod.nwt) {
        w = dset.Z[pmod.nwt];
        ut = pmod.uhat[t1 - 1];
        if (!na(ut)) {
            den += ut * ut;
        }
    } else {
        den = pmod.ess;
    }

    for (t=t1; t<=pmod.t2; t++)  {
        ut = pmod.uhat[t];
        u1 = pmod.uhat[t-1];
        if (na(ut) || na(u1) ||
            (w != null && (w[t] == 0.0 || w[t-1] == 0.0))) {
            continue;
        }
        num += (ut - u1) * (ut - u1);
        if (pmod.nwt) {
            den += ut * ut;
        }
    }

    return num / den;
}

/* altrho: alternative calculation of rho */

private double altrho(int order, int t1, int t2, const(double) * uhat)
{
    int t, n, len = t2 - (t1 + order) + 1;
    double * ut, u1;
    double uht, uh1, rho;

    ut = cast(double *) malloc(2 * len * double.sizeof);
    if (ut == null) {
        return NADBL;
    }

    u1 = ut + len;
    n = 0;

    for (t=t1+order; t<=t2; t++) {
        uht = uhat[t];
        uh1 = (t > 0)? uhat[t-1] : NADBL;
        if (!na(uht) && !na(uh1)) {
            ut[n] = uht;
            u1[n] = uh1;
            n++;
        }
    }

    rho = gretl_corr(0, n - 1, ut, u1, null);

    free(ut);

    return rho;
}

/**
 * rhohat:
 * @order: order of autoregression, usually 1.
 * @t1: start of sample range.
 * @t2: end of sample range.
 * @uhat: array of regression residuals.
 *
 * Computes the first order serial correlation coefficient
 * for @uhat, over the range @t1 to @t2.
 *
 * Returns: the \hat{rho} value, or #NADBL on error.
 */

double rhohat(int order, int t1, int t2, const(double) * uhat)
{
    double ut, u1, num = 0.0, den = 0.0;
    double rho;
    int t;

    for (t=t1+order; t<=t2; t++) {
        ut = uhat[t];
        u1 = uhat[t-1];
        if (na(ut) || na(u1)) {
            continue;
        }
        num += ut * u1;
        den += u1 * u1;
    }

    if (den < DBL_EPSILON) {
        return NADBL;
    }

    rho = num / den;

    if (rho > 1.0 || rho < -1.0) {
        /* out of bounds, try again */
        rho = altrho(order, t1, t2, uhat);
    }

    return rho;
}

private int hilu_plot(double * ssr, double * rho, int n)
{
    FILE * fp;
    int i, err = 0;

    fp = open_plot_input_file(PlotType.regular, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    fputs("# hildreth-lu\n", fp);
    fputs("set xlabel 'rho'\n", fp);

    fprintf(fp, "set ylabel '%s'\n", _("ESS").ptr);

    fputs("set nokey\n", fp);
    fputs("set xrange [-1.0:1.0]\n", fp);
    fputs("plot '-' using 1:2 w impulses\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        fprintf(fp, "%g %g\n", rho[i], ssr[i]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* calculation of rhohat for use with CORC/PWE */

private double autores(Model * pmod, DatasetRaw * dset, Gretlopt opt)
{
    double * uhat = pmod.uhat;
    double ut, num = 0.0, den = 0.0;
    int yno = pmod.list[1];
    int t1 = pmod.t1;
    int t, i, v;

    if (!(opt & Gretlopt.OPT_P)) {
        /* not using Prais-Winsten */
        t1--;
    }

    for (t=t1; t<=pmod.t2; t++) {
        ut = dset.Z[yno][t];
        for (i=0; i<pmod.ncoeff; i++) {
            v = pmod.list[i+2];
            ut -= pmod.coeff[i] * dset.Z[v][t];
        }
        uhat[t] = ut;
        if (t > t1) {
            num += uhat[t] * uhat[t-1];
            den += uhat[t-1] * uhat[t-1];
        }
    }

    return num / den;
}

private int rho_val_ends_row(int k)
{
    int i, endval = -80; /* rho = -0.80 */

    for (i=0; i<7; i++) {
        if (k == endval) {
            return 1;
        }
        endval += 30;
    }

    return 0;
}

private void hilu_print_iter(double r, double ess, int iter,
                             double * ssr, double * rh,
                             int * ip, Prn * prn)
{
    char[8] num;
    int k;

    sprintf(num.ptr, "%.1f", 100 * r);
    k = atoi(num.ptr);

    /* we'll not print all 199 rho values */

    if (k == -99 || k == 99 || k % 10 == 0) {
        pprintf(prn, " %5.2f %#12.6g", r, ess);
        if (rho_val_ends_row(k)) {
            pputc(prn, '\n');
        } else {
            bufspace(3, prn);
        }
    }
}

private void hilu_header(Prn * prn)
{
    int i;

    pputs(prn, "\n   ");

    for (i=0; i<3; i++) {
        pputs(prn, "rho");
        bufspace(10, prn);
        pputs(prn, _("ESS").ptr);
        if (i < 2) {
            pputs(prn, "      ");
        }
    }

    pputc(prn, '\n');
}

private double hilu_search(const(int) * list, DatasetRaw * dset,
                           Model * pmod, Gretlopt opt,
                           Prn * prn)
{
    double * rh = null, ssr = null;
    double essmin = 1.0e200;
    double ess, r, hl_rho = 0;
    int quiet = (opt & Gretlopt.OPT_Q) != 0;
    int do_plot = (opt & Gretlopt.OPT_G) != 0;
    int iter = 0, i = 0;
    int nplot = 0;

    if (!quiet) {
        /* allocate arrays for recording plot info */
        nplot = do_plot ? 199 : 21;
        ssr = cast(double *) malloc(2 * nplot * double.sizeof);
        if (ssr == null) {
            pmod.errcode = GretlError.alloc;
            return NADBL;
        }
        rh = ssr + nplot;
        /* and print header */
        hilu_header(prn);
    }

    for (r = -0.99; r < 1.0; r += .01, iter++) {
        clear_model(pmod);

        *pmod = ar1_lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A, r);
        if (pmod.errcode) {
            free(ssr);
            return NADBL;
        }

        ess = pmod.ess;

        if (!quiet) {
            hilu_print_iter(r, ess, iter, ssr, rh, &i, prn);
            if (do_plot) {
                /* record full data for plotting */
                ssr[i] = ess;
                rh[i++] = r;
            }
        }

        if (iter == 0 || ess < essmin) {
            essmin = ess;
            hl_rho = r;
        }
    } /* end of basic iteration */

    if (hl_rho > 0.989) {
        /* try exploring this funny region? */
        for (r = 0.99; r <= 0.999; r += .001) {
            clear_model(pmod);
            *pmod = ar1_lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A, r);
            if (pmod.errcode) {
                free(ssr);
                return NADBL;
            }
            ess = pmod.ess;
            if (ess < essmin) {
                essmin = ess;
                hl_rho = r;
            }
        }
    }

    if (hl_rho > 0.9989) {
        /* this even funnier one? */
        for (r = 0.9991; r <= 0.9999; r += .0001) {
            clear_model(pmod);
            *pmod = ar1_lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A, r);
            if (pmod.errcode) {
                free(ssr);
                return NADBL;
            }
            ess = pmod.ess;
            if (ess < essmin) {
                essmin = ess;
                hl_rho = r;
            }
        }
    }

    if (!quiet) {
        pprintf(prn, _("\n\nESS is minimized for rho = %g\n\n").ptr, hl_rho);
        if (do_plot) {
            hilu_plot(ssr, rh, nplot);
        }
    }

    free(ssr);

    return hl_rho;
}

private void corc_header(Gretlopt opt, Prn * prn)
{
    if (opt & Gretlopt.OPT_H) {
        pputs(prn, _("\nFine-tune rho using the CORC procedure...\n\n").ptr);
    } else {
        pputs(prn, _("\nPerforming iterative calculation of rho...\n\n").ptr);
    }

    pputs(prn, _("                 ITER       RHO        ESS").ptr);
    pputc(prn, '\n');
}

/* Unlike plain CORC, PWE needs the term sqrt(1 - rho^2), so
   we really cannot proceed if |rho| >= 1.0, even just to see
   where rho ends up, as we now do with CORC.
*/

private int pwe_error(double rho, int quiet, int iter, Prn * prn)
{
    gretl_errmsg_sprintf(_("Prais-Winsten: invalid rho value %g "
                           ~ "(explosive error)").ptr, rho);
    if (!quiet) {
        pprintf(prn, "          %10d %12.5f", iter, rho);
        pprintf(prn, "        NA\n");
    }

    return GretlError.noconv;
}

/*
 * estimate_rho:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: option flags: may include OPT_H to use Hildreth-Lu,
 *       OPT_P to use Prais-Winsten, OPT_B to suppress Cochrane-Orcutt
 *       fine-tuning of Hildreth-Lu results, OPT_G to generate
 *       a gnuplot graph of the search in Hildreth-Lu case.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Estimate the quasi-differencing coefficient for use with the
 * Cochrane-Orcutt, Hildreth-Lu or Prais-Winsten procedures for
 * handling first-order serial correlation. Print a trace of the
 * search for rho.
 *
 * Returns: rho estimate on successful completion, #NADBL on error.
 */

private double estimate_rho(const(int) * list, DatasetRaw * dset,
                            Gretlopt opt, Prn * prn, int * err)
{
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int quiet = (opt & Gretlopt.OPT_Q) != 0;
    int pwe = (opt & Gretlopt.OPT_P) != 0;
    double rho = 0.0;
    Model armod;

    *err = list_adjust_sample(list, &dset.t1, &dset.t2,
                              dset, null);
    if (*err) {
        goto bailout;
    }

    gretl_model_init(&armod, dset);

    if (opt & Gretlopt.OPT_H) {
        /* Do Hildreth-Lu first */
        rho = hilu_search(list, dset, &armod, opt, prn);
    } else {
        /* Initialize Cochrane-Orcutt (or Prais-Winsten) */
        armod = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
        if (!armod.errcode && armod.dfd == 0) {
            armod.errcode = GretlError.df;
        }
        if (!armod.errcode) {
            rho = armod.rho;
        }
    }

    if (armod.errcode) {
        *err = armod.errcode;
    } else if (!(opt & Gretlopt.OPT_H) || !(opt & Gretlopt.OPT_B)) {
        /* either not Hildreth-Lu, or Hildreth-Lu but CORC
           fine-tuning is wanted (has not been suppressed
           via OPT_B)
        */
        Gretlopt lsqopt = pwe ? (Gretlopt.OPT_A | Gretlopt.OPT_P) : Gretlopt.OPT_A;
        double edsmall = 5.0e-8;
        double ess0 = armod.ess;
        double rho0 = rho;
        double rdsmall = (opt & Gretlopt.OPT_L)? 0.001 : 1.0e-6;
        int iter = 0, itermax = 100;
        int converged = 0;

        if (!quiet) {
            corc_header(opt, prn);
        }

        while (++iter <= itermax && !converged) {
            clear_model(&armod);
            armod = ar1_lsq(list, dset, GretlCmdIndex.ols, lsqopt, rho);
            if ((*err = armod.errcode) != 0) {
                break;
            }

            if (!quiet) {
                pprintf(prn, "          %10d %12.5f", iter, rho);
                pprintf(prn, "   %#.6g\n", armod.ess);
            }

            rho = autores(&armod, dset, opt);

            if (pwe && fabs(rho) >= 1.0) {
                *err = pwe_error(rho, quiet, iter + 1, prn);
                break;
            }

            if (armod.ess == 0.0) {
                converged = 1;
            } else if ((ess0 - armod.ess) / ess0 < edsmall) {
                converged = 1;
            } else if (fabs(rho - rho0) < rdsmall) {
                converged = 1;
            }

            ess0 = armod.ess;
            rho0 = rho;
        } /* end Cochrane-Orcutt iteration */

        if (converged && (rho >= 1.0 || rho <= -1.0)) {
            gretl_errmsg_sprintf(_("%s: rho = %g, error is unbounded").ptr,
                                 "ar1".ptr, rho);
            converged = 0;
        }

        if (!converged && !*err) {
            *err = GretlError.noconv;
        }

        if (!quiet) {
            if (*err) {
                pputc(prn, '\n');
            } else {
                pprintf(prn, "          %10d %12.5f", iter, rho);
                pprintf(prn, "   %#.6g\n", armod.ess);
            }
        }
    }

    clear_model(&armod);

 bailout:

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    if (*err) {
        rho = NADBL;
    }

    return rho;
}

/**
 * augment_regression_list:
 * @orig: list giving original regression specification.
 * @aux: either %AUX_SQ, %AUX_LOG or %AUX_WHITE.
 * @dset: dataset struct.
 * @err: location to receive error code.
 *
 * Augment the regression list @orig with auxiliary terms.  If @aux
 * is %AUX_SQ add the squares of the original regressors; if @aux
 * is %AUX_WHITE add squares and cross-products, or if @aux is
 * %AUX_LOG add the natural logs of the original regressors.
 * If they are not already present, these variables are added
 * to the data array.
 *
 * Returns: the augmented list, or NULL on failure.
 */

int * augment_regression_list(const(int) * orig, int aux,
                              DatasetRaw * dset, int * err)
{
    int * list;
    int listlen;
    int v_prev = dset.v;
    int k = orig[0];
    int cnum = 0;
    int i;

    if (aux == ModelAuxCode.white) {
        int cpos = gretl_list_const_pos(orig, 2, dset);
        int nt, trv = orig[0] - 1;

        if (cpos > 0) {
            trv--;
            cnum = orig[cpos];
        }
        nt = (trv * trv + trv) / 2;
        listlen = orig[0] + nt + 1;
    } else {
        listlen = 2 * orig[0];
    }

    list = cast(int *) malloc(listlen * int.sizeof);
    if (list == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* transcribe original list */
    for (i=0; i<=orig[0]; i++) {
        list[i] = orig[i];
    }

    /* add squares, cross-products or logs of independent vars */

    for (i=2; i<=orig[0]; i++) {
        int vnew, vi = orig[i];

        if (vi == 0) {
            continue;
        }

        if (aux == ModelAuxCode.sq || aux == ModelAuxCode.white) {
            vnew = xpxgenr(vi, vi, dset);
            if (vnew > 0) {
                list[++k] = vnew;
            }
            if (aux == ModelAuxCode.white) {
                int j, vj;

                for (j=i+1; j<=orig[0]; j++) {
                    vj = orig[j];
                    if (vj == cnum) {
                        continue;
                    }
                    vnew = xpxgenr(vi, vj, dset);
                    if (vnew > 0) {
                        if (vnew >= v_prev) {
                            /* ensure uniqueness of generated varnames */
                            sprintf(dset.varname[vnew], "X%d_X%d", i-1, j-1);
                        }
                        list[++k] = vnew;
                    }
                }
            }
        } else if (aux == ModelAuxCode.log) {
            /* don't try to log series with non-positive values */
            if (gretl_ispositive(dset.t1, dset.t2, dset.Z[vi], 1)) {
                vnew = loggenr(vi, dset);
                if (vnew > 0) {
                    list[++k] = vnew;
                }
            }
        }
    }

    list[0] = k;

    return list;
}

/* For observation s, see if the regression list contains a variable
   that has a single non-zero value at that particular observation.
   We run this check on observations showing an OLS residual of zero.
*/

private int observation_is_dummied(const(Model) * pmod,
                                   int * list, const(double *) * Z,
                                   int s)
{
    int i, t, v;
    int ret = 0;

    for (i=list[0]; i>=2; i--) {
        v = list[i];
        if (v == 0) {
            continue;
        }
        ret = 1;
        for (t=pmod.t1; t<=pmod.t2; t++) {
            if ((t == s && Z[v][t] == 0.0) || (t != s && Z[v][t] != 0.0)) {
                ret = 0;
                break;
            }
        }
        if (ret) {
            gretl_list_delete_at_pos(list, i);
            break;
        }
    }

    return ret;
}

private int * copy_ensure_const(const(int) * orig)
{
    int * list = gretl_list_new(orig[0] + 1);
    int i;

    if (list != null) {
        list[1] = orig[1];
        list[2] = 0;
        for (i=3; i<=list[0]; i++) {
            list[i] = orig[i-1];
        }
    }

    return list;
}

/* get_hsk_weights: take the residuals from the model @pmod, square
   them and take logs; find the fitted values for this series using an
   auxiliary regression including the original independent variables
   (and their squares, if not given OPT_N); exponentiate the fitted
   values; and add the resulting series to the data set.
*/

private int get_hsk_weights(Model * pmod, DatasetRaw * dset, Gretlopt opt)
{
    int oldv = dset.v;
    int t, t1 = dset.t1, t2 = dset.t2;
    int * lcpy = null;
    int * auxlist = null;
    int err = 0, shrink = 0;
    double xx;
    Model aux;

    if (pmod.ifc) {
        lcpy = gretl_list_copy(pmod.list);
    } else {
        lcpy = copy_ensure_const(pmod.list);
    }

    if (lcpy == null) {
        return GretlError.alloc;
    }

    /* allocate space for an additional variable */
    if (dataset_add_series(dset, 1)) {
        free(lcpy);
        return GretlError.alloc;
    }

    /* add transformed residuals to data set */
    for (t=0; t<dset.n; t++) {
        xx = pmod.uhat[t];
        if (na(xx)) {
            dset.Z[oldv][t] = NADBL;
        } else if (xx == 0.0) {
            if (observation_is_dummied(pmod, lcpy, cast(const(double *) *) dset.Z, t)) {
                dset.Z[oldv][t] = NADBL;
            } else {
                fprintf(stderr, "hsk: got a zero residual, could be a problem!\n");
                dset.Z[oldv][t] = -1.0e16; /* ?? */
            }
        } else {
            dset.Z[oldv][t] = log(xx * xx);
        }
    }

    if (opt & Gretlopt.OPT_N) {
        /* --no-squares option */
        auxlist = lcpy;
    } else {
        /* build regression list, adding the squares of the original
           independent vars */
        auxlist = augment_regression_list(lcpy, ModelAuxCode.sq, dset, &err);
        if (err) {
            return err;
        }
    }

    auxlist[1] = oldv; /* the newly added log(uhat-squared) */

    dset.t1 = pmod.t1;
    dset.t2 = pmod.t2;

    aux = lsq(auxlist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = aux.errcode;
    if (err) {
        shrink = dset.v - oldv;
    } else {
        /* write into the data set the required weights */
        for (t=aux.t1; t<=aux.t2; t++) {
            xx = aux.yhat[t];
            if (na(xx)) {
                dset.Z[oldv][t] = NADBL;
            } else {
                dset.Z[oldv][t] = 1.0 / exp(xx);
            }
        }
        shrink = dset.v - oldv - 1;
    }

    dset.t1 = t1;
    dset.t2 = t2;

    clear_model(&aux);

    if (shrink > 0) {
        dataset_drop_last_variables(dset, shrink);
    }

    if (auxlist != lcpy) {
        free(auxlist);
    }
    free(lcpy);

    return err;
}

/**
 * hsk_model:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: may include OPT_N to suppress use of squares
 * of regressors in the auxiliary regression.
 *
 * Estimate the model given in @list using a correction for
 * heteroskedasticity.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model hsk_model(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    int i, err;
    int orig_nvar = dset.v;
    int * hsklist;
    Model hsk;

    if (dset.Z == null) {
        hsk.errcode = GretlError.data;
        return hsk;
    }

    gretl_error_clear();

    /* run initial OLS */
    hsk = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (hsk.errcode) {
        return hsk;
    }

    /* form weights based on the OLS residuals */
    err = get_hsk_weights(&hsk, dset, opt);
    if (err) {
        hsk.errcode = err;
        return hsk;
    }

    /* allocate regression list for weighted least squares */
    hsklist = gretl_list_new(list[0] + 1);
    if (hsklist == null) {
        hsk.errcode = GretlError.alloc;
        return hsk;
    }

    /* the last variable in the dataset will be the weight var */
    hsklist[1] = dset.v - 1;

    /* put the original dependent variable in at position 2 */
    hsklist[2] = list[1];

    /* add the original independent vars into the WLS list */
    for (i=3; i<=hsklist[0]; i++) {
        hsklist[i] = list[i-1];
    }

    clear_model(&hsk);
    hsk = lsq(hsklist, dset, GretlCmdIndex.wls, Gretlopt.OPT_NONE);
    hsk.ci = GretlCmdIndex.hsk;

    if (opt & Gretlopt.OPT_N) {
        gretl_model_set_int(&hsk, "no-squares", 1);
        hsk.opt |= Gretlopt.OPT_N;
    }

    if (!hsk.errcode) {
        hsk.xlist = gretl_list_sublist(hsk.list, 2, -1);
    }

    dataset_drop_last_variables(dset, dset.v - orig_nvar);

    free(hsklist);

    return hsk;
}

private void print_HET_1(double z, double pval, Prn * prn)
{
    pprintf(prn, "\n%s\n", _("Pesaran-Taylor test for heteroskedasticity").ptr);
    pprintf(prn, "\n%s: HET_1 = %f,\n", _("Test statistic").ptr, z);
    pprintf(prn, "%s = 2 * P(z > %f) = %.3g\n\n",
            _("with p-value").ptr, z, pval);
}

private void print_whites_test(double LM, int df, double pval,
                               Gretlopt opt, Prn * prn)
{
    if (opt & Gretlopt.OPT_B) {
        pprintf(prn, "\n%s\n", _("Breusch-Pagan test for heteroskedasticity").ptr);
        if (opt & Gretlopt.OPT_R) {
            pprintf(prn, "%s\n", _("with Koenker robust variance estimator").ptr);
        }
        pprintf(prn, "\n%s: LM = %f,\n", _("Test statistic").ptr, LM);
    } else {
        pprintf(prn, "\n%s", _("White's test for heteroskedasticity").ptr);
        if (opt & Gretlopt.OPT_X) {
            pprintf(prn, " (%s)\n", _("squares only").ptr);
        } else {
            pputc(prn, '\n');
        }
        pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic").ptr, LM);
    }

    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n",
            _("with p-value").ptr, _("Chi-square").ptr,
            df, LM, pval);
}

/* For White's test, see if we have enough degrees of freedom
   to add squares -- and if so, if we have enough df to add
   cross-products also.
*/

private int get_whites_aux(const(Model) * pmod, const(double *) * Z)
{
    int aux = ModelAuxCode.white;
    int rem = pmod.ncoeff - pmod.ifc - 1;
    int nsq = 0, nxpx = 0;
    int i, v;

    for (i=2; i<=pmod.list[0]; i++) {
        v = pmod.list[i];
        if (v > 0) {
            if (!gretl_isdummy(pmod.t1, pmod.t2, Z[v])) {
                nsq++;
            }
            nxpx += rem--;
        }
    }

    if (pmod.dfd - nsq < 1) {
        aux = ModelAuxCode.none;
    } else if (pmod.dfd - nsq - nxpx < 1) {
        aux = ModelAuxCode.sq;
    }

    return aux;
}

/**
 * tsls_hetero_test:
 * @pmod: pointer to model.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save results to model; OPT_Q
 * means don't print the auxiliary regression.
 * @prn: gretl printing struct.
 *
 * Runs Pesaran and Taylor's (1999) HET_1 test for heteroskedasticity
 * on the given tsls model. The statistic is just a t-statistic, so
 * under the null it is distributed as a standard normal.
 *
 * Returns: 0 on successful completion, error code on error.
 */

private int tsls_hetero_test(Model * pmod, DatasetRaw * dset,
                             Gretlopt opt, Prn * prn)
{
    int pos, newv = dset.v;
    int * ptlist = null, testlist = null;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    Model ptmod;
    double x;
    int i, h, t;
    int err = 0;

    pos = gretl_list_separator_position(pmod.list);
    h = pmod.list[0] - pos;

    ptlist = gretl_list_new(h + 1);
    testlist = gretl_list_new(3);

    if (ptlist == null || testlist == null) {
        free(ptlist);
        free(testlist);
        return GretlError.alloc;
    }

    ptlist[1] = pmod.list[1];
    for (i=2; i<=ptlist[0]; i++) {
        ptlist[i] = pmod.list[i + pos - 1];
    }

    testlist[1] = newv;
    testlist[2] = 0;
    testlist[3] = newv + 1;

    /* reduced form: regress the original dependent variable on all of
       the instruments from the original model
    */
    ptmod = lsq(ptlist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = ptmod.errcode;
    if (err) {
        goto bailout;
    }

    /* add two series: (i) the squares of the residuals from the
       original model and (ii) the squares of the fitted values from
       the auxiliary regression above
     */
    err = dataset_add_series(dset, 2);
    if (err) {
        clear_model(&ptmod);
        goto bailout;
    }

    strcpy(dset.varname[newv+1], "yhat^2");

    for (t=pmod.t1; t<=pmod.t2; t++) {
        x = pmod.uhat[t];
        dset.Z[newv][t] = x*x;   /* squared residual */
        x = ptmod.yhat[t];
        dset.Z[newv+1][t] = x*x; /* squared fitted value */
    }

    clear_model(&ptmod);

    dset.t1 = pmod.t1;
    dset.t2 = pmod.t2;

    /* regress the squared residuals on the squared fitted
       values from the reduced-form auxiliary regression
    */
    ptmod = lsq(testlist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = ptmod.errcode;

    if (!err) {
        double z = fabs(ptmod.coeff[1]) / ptmod.sderr[1];
        double pval = 2.0 * (1 - normal_cdf(z));

        if (opt & Gretlopt.OPT_Q) {
            print_HET_1(z, pval, prn);
        } else {
            ptmod.aux = ModelAuxCode.het_1;
            printmodel(&ptmod, dset, Gretlopt.OPT_NONE, prn);
        }

        if (opt & Gretlopt.OPT_S) {
            ModelTest * test = model_test_new(ModelTestType.het_1);

            if (test != null) {
                model_test_set_teststat(test, cast(ubyte) TestStats.z);
                model_test_set_value(test, z);
                model_test_set_pvalue(test, pval);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(z, pval);
    }

    clear_model(&ptmod);

    dataset_drop_last_variables(dset, 2);

 bailout:

    free(ptlist);
    free(testlist);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

/* Compute LM statistic as per Breusch and Pagan (Econometrica, 1979),
   with the option to use the robust variance estimator proposed by
   Koenker (Journal of Econometrics, 1981).
*/

private double get_BP_LM(Model * pmod, int * list, Model * aux,
                         DatasetRaw * dset, Gretlopt opt, int * err)
{
    double s2, u2t, gt;
    double V = 0.0, LM = NADBL;
    int t, v = list[1];

    /* note: no df adjustment */
    s2 = pmod.ess / pmod.nobs;

    if (opt & Gretlopt.OPT_R) {
        /* calculate robust variance estimate a la Koenker */
        for (t=pmod.t1; t<=pmod.t2; t++) {
            if (!na(pmod.uhat[t])) {
                u2t = pmod.uhat[t] * pmod.uhat[t];
                V += (u2t - s2) * (u2t - s2);
            }
        }
        V /= pmod.nobs;
    }

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (na(pmod.uhat[t])) {
            dset.Z[v][t] = NADBL;
        } else {
            u2t = pmod.uhat[t] * pmod.uhat[t];
            gt = (opt & Gretlopt.OPT_R)? (u2t - s2) : (u2t / s2);
            dset.Z[v][t] = gt;
        }
    }

    *aux = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    *err = aux.errcode;

    if (!*err) {
        double RSS = aux.tss - aux.ess;

        if (RSS < 0) {
            *err = GretlError.data;
        } else {
            if (opt & Gretlopt.OPT_R) {
                aux.opt |= Gretlopt.OPT_R;
                LM = RSS / V;
            } else {
                LM = .5 * RSS;
            }
            gretl_model_set_double(aux, "BPLM", LM);
            aux.aux = ModelAuxCode.bp;
        }
    }

    return LM;
}

private int * ensure_const_copy_list(const(Model) * pmod, int * err)
{
    int * list = null;

    if (pmod.ifc) {
        list = gretl_list_copy(pmod.list);
    } else {
        int i, nl = pmod.list[0] + 1;

        list = gretl_list_new(nl);
        if (list != null) {
            list[1] = pmod.list[1];
            list[2] = 0; /* insert const */
            for (i=3; i<=nl; i++) {
                list[i] = pmod.list[i-1];
            }
        }
    }

    if (list == null) {
        *err = GretlError.alloc;
    }

    return list;
}

private int * ensure_const_augment_list(const(Model) * pmod, int aux,
                                        DatasetRaw * dset, int * err)
{
    int * list = null;

    if (pmod.ifc) {
        list = augment_regression_list(pmod.list, aux, dset, err);
    } else {
        int * tmp = ensure_const_copy_list(pmod, err);

        if (!*err) {
            list = augment_regression_list(tmp, aux, dset, err);
            free(tmp);
        }
    }

    return list;
}

/**
 * whites_test:
 * @pmod: pointer to model.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save results to model; OPT_Q
 * means don't print the auxiliary regression; OPT_B means
 * do the simpler Breusch-Pagan variant; OPT_I means run
 * silently.
 * @prn: gretl printing struct.
 *
 * Runs White's test for heteroskedasticity on the given model.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int whites_test(Model * pmod, DatasetRaw * dset,
                Gretlopt opt, Prn * prn)
{
    int BP = (opt & Gretlopt.OPT_B) != 0;
    int aux = ModelAuxCode.none;
    int v = dset.v;
    int * list = null;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    double zz, LM = 0;
    Model white;
    int t, err = 0;

    if (pmod.ci == GretlCmdIndex.ivreg) {
        return tsls_hetero_test(pmod, dset, opt, prn);
    }

    if (pmod.list == null || gretl_list_has_separator(pmod.list)) {
        return GretlError.notimp;
    }

    if (pmod.ci == GretlCmdIndex.nls || pmod.ci == GretlCmdIndex.mle ||
        pmod.ci == GretlCmdIndex.gmm || pmod.ci == GretlCmdIndex.arma ||
        pmod.ci == GretlCmdIndex.logistic) {
        return GretlError.notimp;
    }

    if ((err = list_members_replaced(pmod, dset)) != 0) {
        return err;
    }

    impose_model_smpl(pmod, dset);

    /* what can we do, with the degrees of freedom available? */
    if (!BP) {
        if (opt & Gretlopt.OPT_X) {
            aux = ModelAuxCode.sq;
        } else {
            aux = get_whites_aux(pmod, cast(const(double *) *) dset.Z);
            if (aux == ModelAuxCode.none) {
                dset.t1 = save_t1;
                dset.t2 = save_t2;
                return GretlError.df;
            }
        }
    }

    gretl_model_init(&white, dset);

    /* make space in data set */
    if (dataset_add_series(dset, 1)) {
        err = GretlError.alloc;
    }

    if (!err) {
        /* get residuals, square and add to data set */
        for (t=0; t<dset.n; t++) {
            zz = pmod.uhat[t];
            if (na(zz)) {
                dset.Z[v][t] = NADBL;
            } else {
                dset.Z[v][t] = zz * zz;
            }
        }
        strcpy(dset.varname[v], "uhatsq");
    }

    if (!err) {
        if (BP) {
            list = ensure_const_copy_list(pmod, &err);
        } else {
            /* build aux regression list, adding squares and
               (possibly) cross-products of the original
               independent vars */
            list = ensure_const_augment_list(pmod, aux, dset, &err);
        }
        if (!err) {
            list[1] = v; /* the newly added uhat-squared */
        }
    }

    if (!err) {
        /* run auxiliary regression */
        if (BP) {
            LM = get_BP_LM(pmod, list, &white, dset, opt, &err);
        } else {
            white = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
            err = white.errcode;
            if (!err) {
                LM = white.rsq * white.nobs;
                white.aux = ModelAuxCode.white;
            }
        }
    }

    if (!err) {
        int df = white.ncoeff - 1;
        double pval = chisq_cdf_comp(df, LM);
        Gretlopt testopt = Gretlopt.OPT_NONE;

        if (BP) {
            if (opt & Gretlopt.OPT_R) {
                /* Koenker robust variant */
                testopt = Gretlopt.OPT_R;
            }
        } else if (opt & Gretlopt.OPT_X) {
            /* squares only */
            testopt = Gretlopt.OPT_X;
        }

        if (!(opt & Gretlopt.OPT_I)) {
            if (opt & Gretlopt.OPT_Q) {
                print_whites_test(LM, df, pval, opt, prn);
            } else {
                white.opt |= testopt;
                printmodel(&white, dset, Gretlopt.OPT_NONE, prn);
            }
        }

        if (opt & Gretlopt.OPT_S) {
            ModelTestType mt = BP? ModelTestType.bp : ModelTestType.whites;
            ModelTest * test = model_test_new(mt);

            if (test != null) {
                model_test_set_teststat(test, cast(ubyte) TestStats.lm);
                model_test_set_dfn(test, df);
                model_test_set_value(test, LM);
                model_test_set_pvalue(test, pval);
                model_test_set_opt(test, testopt);
                maybe_add_test_to_model(pmod, test);
            }
        }

        record_test_result(LM, pval);
    }

    clear_model(&white);
    dataset_drop_last_variables(dset, dset.v - v);
    free(list);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

private int ar_list_max(const(int) * list)
{
    int i, lmax = 0;

    for (i=1; i<=list[0]; i++) {
        if (list[i] > lmax) {
            lmax = list[i];
        }
    }

    return lmax;
}

/**
 * ar_model:
 * @list: list of lags plus dependent variable and list of regressors.
 * @dset: dataset struct.
 * @opt: may contain OPT_O to print covariance matrix.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the generalized
 * Cochrane-Orcutt procedure for autoregressive errors.
 *
 * Returns: #MODEL struct containing the results.
 */

Model ar_model(const(int) * list, DatasetRaw * dset,
               Gretlopt opt, Prn * prn)
{
    double diff, ess, tss, xx;
    int i, j, vi, t, t1, t2, vc, yno, ryno, iter;
    int k, pmax, v = dset.v;
    int * arlist = null, rholist = null;
    int * reglist = null, reglist2 = null;
    int cpos, nvadd;
    Model ar, rhomod;

    gretl_error_clear();
    gretl_model_init(&ar, dset);
    gretl_model_init(&rhomod, dset);

    ar.errcode = gretl_list_split_on_separator(list, &arlist, &reglist);

    if (!ar.errcode && reglist != null && reglist[0] == 1) {
        gretl_errmsg_set(_("At least one regressor is required").ptr);
        ar.errcode = GretlError.args;
        goto bailout;
    }

    if (!ar.errcode && arlist[0] == 1 && arlist[1] == 1) {
        /* special case: ar 1 ; ... => use AR1 apparatus */
        ar = ar1_model(reglist, dset, opt, prn);
        goto bailout;
    }

    if (!ar.errcode) {
        reglist2 = gretl_list_new(reglist[0]);
        rholist = gretl_list_new(arlist[0] + 1);
        if (reglist2 == null || rholist == null) {
            ar.errcode = GretlError.alloc;
        }
    }

    if (ar.errcode) {
        goto bailout;
    }

    pmax = ar_list_max(arlist);
    cpos = reglist_check_for_const(reglist, dset);

    /* first pass: estimate model via OLS: use OPT_M to generate an
       error in case of missing values within sample range
    */
    ar = lsq(reglist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_M);
    if (ar.errcode) {
        goto bailout;
    }

    nvadd = arlist[0] + 1 + reglist[0];

    /* allocate space for the uhat terms and transformed data */
    if (dataset_add_series(dset, nvadd)) {
        ar.errcode = GretlError.alloc;
        goto bailout;
    }

    yno = reglist[1];
    t1 = ar.t1;
    t2 = ar.t2;
    rholist[1] = v;

    pprintf(prn, "%s\n\n", _("Generalized Cochrane-Orcutt estimation").ptr);
    bufspace(17, prn);
    /* xgettext:no-c-format */
    pputs(prn, _("ITER             ESS           % CHANGE").ptr);
    pputs(prn, "\n\n");

    /* now loop while ess is changing */
    diff = 1.0e6;
    ess = 0.0;

    for (iter = 1; iter <= 20 && diff > 0.005; iter++) {
        for (t=0; t<dset.n; t++) {
            if (t < t1 || t > t2) {
                dset.Z[v][t] = NADBL;
            } else {
                /* special computation of uhat */
                xx = dset.Z[yno][t];
                for (j=0; j<reglist[0]-1; j++) {
                    xx -= ar.coeff[j] * dset.Z[reglist[j+2]][t];
                }
                dset.Z[v][t] = xx;
            }
        }
        for (i=1; i<=arlist[0]; i++) {
            k = arlist[i];
            rholist[1+i] = v + i;
            for (t=0; t<dset.n; t++) {
                if (t < t1 + k || t > t2) {
                    dset.Z[v+i][t] = NADBL;
                } else {
                    dset.Z[v+i][t] = dset.Z[v][t-k];
                }
            }
        }

        /* estimate the rho terms */
        if (iter > 1) {
            clear_model(&rhomod);
        }
        rhomod = lsq(rholist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);

        /* and rho-transform the data */
        ryno = vc = v + i;
        for (i=1; i<=reglist[0]; i++) {
            vi = reglist[i];
            for (t=0; t<dset.n; t++) {
                if (t < t1 + pmax || t > t2) {
                    dset.Z[vc][t] = NADBL;
                } else {
                    xx = dset.Z[vi][t];
                    for (j=1; j<=arlist[0]; j++) {
                        k = arlist[j];
                        xx -= rhomod.coeff[j-1] * dset.Z[vi][t-k];
                    }
                    dset.Z[vc][t] = xx;
                }
            }
            reglist2[i] = vc++;
        }

        /* estimate the transformed model */
        clear_model(&ar);
        ar = lsq(reglist2, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);

        if (iter > 1) {
            diff = fabs(100 * (ar.ess - ess) / ess);
        }

        ess = ar.ess;
        pprintf(prn, "%16c%3d %20f ", ' ', iter, ess);

        if (iter > 1) {
            pprintf(prn, "%13.3f\n", diff);
        } else {
            pputc(prn, '\n');
        }
    } /* end "ess changing" loop */

    for (i=0; i<=reglist[0]; i++) {
        ar.list[i] = reglist[i];
    }
    if (cpos > 0) {
        ar.ifc = 1;
    }
    if (ar.ifc) {
        if (!gretl_model_get_int(&ar, "effconst")) {
            ar.dfn -= 1;
        }
    }
    ar.ci = GretlCmdIndex.ar;

    /* special computation of fitted values */
    for (t=t1; t<=t2; t++) {
        xx = 0.0;
        for (j=2; j<=reglist[0]; j++) {
            xx += ar.coeff[j-2] * dset.Z[reglist[j]][t];
        }
        ar.uhat[t] = dset.Z[yno][t] - xx;
        for (j=1; j<=arlist[0]; j++) {
            if (t - t1 >= arlist[j]) {
                k = arlist[j];
                xx += rhomod.coeff[j-1] * ar.uhat[t-k];
            }
        }
        ar.yhat[t] = xx;
    }

    for (t=t1; t<=t2; t++) {
        ar.uhat[t] = dset.Z[yno][t] - ar.yhat[t];
    }

    ar.rsq = gretl_corr_rsq(ar.t1, ar.t2, dset.Z[reglist[1]], ar.yhat);
    ar.adjrsq = 1.0 - ((1.0 - ar.rsq) * (ar.nobs - 1.0) / ar.dfd);

    /* special computation of TSS */
    xx = gretl_mean(ar.t1, ar.t2, dset.Z[ryno]);
    tss = 0.0;
    for (t=ar.t1; t<=ar.t2; t++) {
        tss += (dset.Z[ryno][t] - xx) * (dset.Z[ryno][t] - xx);
    }
    if (ar.dfn == 0) {
        ar.fstt = NADBL;
    } else {
        ar.fstt = ar.dfd * (tss - ar.ess) / (ar.dfn * ar.ess);
    }
    ar.lnL = NADBL;
    mle_criteria(&ar, 0);
    ar.dw = dwstat(pmax, &ar, dset);
    ar.rho = rhohat(pmax, ar.t1, ar.t2, ar.uhat);

    dataset_drop_last_variables(dset, nvadd);

    if (gretl_model_add_arinfo(&ar, pmax)) {
        ar.errcode = GretlError.alloc;
    } else {
        double * y = dset.Z[reglist[1]];

        for (i=0; i<=arlist[0]; i++) {
            ar.arinfo.arlist[i] = arlist[i];
            if (i >= 1) {
                ar.arinfo.rho[i-1] = rhomod.coeff[i-1];
                ar.arinfo.sderr[i-1] = rhomod.sderr[i-1];
            }
        }
        ar.ybar = gretl_mean(ar.t1, ar.t2, y);
        ar.sdy = gretl_stddev(ar.t1, ar.t2, y);
    }
    clear_model(&rhomod);

    if (gretl_model_get_int(&ar, "maxlag") < pmax) {
        gretl_model_set_int(&ar, "maxlag", pmax);
    }

    set_model_id(&ar, opt);
    if (!ar.errcode) {
        ar.xlist = gretl_list_sublist(ar.list, 2, -1);
    }

 bailout:

    free(reglist);
    free(reglist2);
    free(rholist);
    free(arlist);

    return ar;
}

private int modelvar_iszero(const(Model) * pmod, const(double) * x)
{
    int t;

    for (t=pmod.t1; t<=pmod.t2; t++) {
        if (!model_missing(pmod, t) && floatneq(x[t], 0.0)) {
            return 0;
        }
    }

    return 1;
}

/* From the position of the first regressor to end of list, omits
   variables with all-zero observations and repacks the rest of
   them. If the regression in question is not just an auxiliary
   one (OPT_A), records the list of dropped variables, if any,
   on @pmod.
*/

private void omitzero(Model * pmod, const(DatasetRaw) * dset,
                      Gretlopt opt)
{
    int * zlist = null;
    int i, v, offset;
    double x = 0.0;

    offset = (pmod.ci == GretlCmdIndex.wls)? 3 : 2;

    if (!(opt & Gretlopt.OPT_A)) {
        zlist = gretl_null_list();
    }

    for (i=pmod.list[0]; i>=offset; i--) {
        v = pmod.list[i];
        if (modelvar_iszero(pmod, dset.Z[v])) {
            if (zlist != null) {
                gretl_list_append_term(&zlist, v);
            }
            fprintf(stderr, "Deleting var %d (%s) at list pos %d: all zero\n",
                    v, dset.varname[v], i);
            gretl_list_delete_at_pos(pmod.list, i);
        }
    }

    if (pmod.nwt) {
        int t, wtzero;

        for (i=pmod.list[0]; i>=offset; i--) {
            v = pmod.list[i];
            wtzero = 1;
            for (t=pmod.t1; t<=pmod.t2; t++) {
                x = dset.Z[v][t] * dset.Z[pmod.nwt][t];
                if (!model_missing(pmod, t) && floatneq(x, 0.0)) {
                    wtzero = 0;
                    break;
                }
            }
            if (wtzero) {
                if (zlist != null) {
                    gretl_list_append_term(&zlist, v);
                }
                gretl_list_delete_at_pos(pmod.list, i);
            }
        }
    }

    if (zlist != null) {
        if (zlist[0] > 0) {
            gretl_model_set_list_as_data(pmod, "zerolist", zlist);
        } else {
            free(zlist);
        }
    }
}

/* model_lags: attempt to detect presence of a lagged dependent
   variable among the regressors -- if found, return its position
   in @list, otherwise return 0. Also write into *pmax the maximum
   lag among the regressors.

   FIXME: use inside a function, when the name of the dependent
   variable may be changed if it was a function argument?
*/

private int model_lags(const(int) * list, const(DatasetRaw) * dset, int * pmax)
{
    int xno, yno = list[1];
    int ilag, maxlag = 0;
    int i, ret = 0;

    for (i=2; i<=list[0]; i++) {
        xno = list[i];
        if (xno == LISTSEP) {
            break;
        }
        ilag = series_get_lag(dset, xno);
        if (ilag > 0) {
            if (ret == 0 && series_get_parent_id(dset, xno) == yno) {
                ret = i;
            }
            if (ilag > maxlag) {
                maxlag = ilag;
            }
        }
    }

    *pmax = maxlag;

    return ret;
}

private void arch_test_print_simple(int order, double LM,
                                    double pval, Gretlopt opt,
                                    Prn * prn)
{
    if (opt & Gretlopt.OPT_M) {
        /* multi-equation system */
        pprintf(prn, "%s: TR^2 = %f,\n", _("Test statistic").ptr, LM);
    } else {
        pprintf(prn, "\n%s %d\n", _("Test for ARCH of order").ptr, order);
        pprintf(prn, "\n%s: TR^2 = %f,\n", _("Test statistic").ptr, LM);
    }
    pprintf(prn, "%s = P(%s(%d) > %f) = %f\n\n",
            _("with p-value").ptr, _("Chi-square").ptr,
            order, LM, pval);
}

private void print_arch_regression(const(GretlMatrixRaw) * b,
                                   const(GretlMatrixRaw) * V,
                                   int T, int order,
                                   Gretlopt opt, Prn * prn)
{
    int i, k = order + 1;
    double * se = cast(double *) malloc(k * double.sizeof);
    char ** names;

    names = strings_array_new_with_length(k, 16);

    if (se != null && names != null) {
        if (!(opt & Gretlopt.OPT_M)) {
            /* not multi-equation */
            pputc(prn, '\n');
            pprintf(prn, _("Test for ARCH of order %d").ptr, order);
            pputs(prn, "\n\n");
        }

        for (i=0; i<k; i++) {
            se[i] = sqrt(gretl_matrix_get(V, i, i));
            sprintf(names[i], "alpha(%d)", i);
        }

        /* is a df correction wanted here? */
        print_coeffs(b.val, se, cast(const(char) **) names,
                     k, T - k, GretlCmdIndex.arch, prn);
    }

    free(se);
    strings_array_free(names, k);
}

private void arch_test_save_or_print(const(GretlMatrixRaw) * b, const(GretlMatrixRaw) * V,
                                     int T, int order, double rsq, Model * pmod,
                                     Gretlopt opt, Prn * prn)
{
    double LM = T * rsq;
    double pv = chisq_cdf_comp(order, LM);

    record_test_result(LM, pv);

    if (!(opt & Gretlopt.OPT_I)) {
        if (V != null) {
            /* V will be NULL if --quiet is in force */
            print_arch_regression(b, V, T, order, opt, prn);
        }
        if (opt & Gretlopt.OPT_Q) {
            arch_test_print_simple(order, LM, pv, opt, prn);
        }
    }

    if ((opt & Gretlopt.OPT_S) || !(opt & Gretlopt.OPT_Q)) {
        ModelTest * test = model_test_new(ModelTestType.arch);

        if (test != null) {
            int heading = (V == null);

            model_test_set_teststat(test, cast(ubyte) TestStats.lm);
            model_test_set_order(test, order);
            model_test_set_dfn(test, order);
            model_test_set_value(test, LM);
            model_test_set_pvalue(test, pv);

            if (!(opt & Gretlopt.OPT_Q)) {
                gretl_model_test_print_direct(test, heading, prn);
            }

            if (pmod != null && (opt & Gretlopt.OPT_S)) {
                maybe_add_test_to_model(pmod, test);
            } else {
                model_test_free(test);
            }
        }
    }
}

private int real_arch_test(const(double) * u, int T, int order,
                           Model * pmod, Gretlopt opt,
                           Prn * prn)
{
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * V = null;
    int i, k, s, t;
    double x, s2, rsq;
    double * ps2 = null;
    int err = 0;

    gretl_error_clear();

    if (order < 1 || order > T - 1) {
        gretl_errmsg_sprintf(_("Invalid lag order for arch (%d)").ptr, order);
        return GretlError.data;
    }

    T -= order;
    k = order + 1;

    X = gretl_matrix_alloc(T, k);
    y = gretl_column_vector_alloc(T);
    b = gretl_column_vector_alloc(k);

    if (X == null || y == null || b == null) {
        gretl_matrix_free(X);
        gretl_matrix_free(y);
        gretl_matrix_free(b);
        return GretlError.alloc;
    }

    if (!(opt & Gretlopt.OPT_Q)) {
        V = gretl_matrix_alloc(k, k);
        if (V == null) {
            err = GretlError.alloc;
            goto bailout;
        }
        ps2 = &s2;
    }

    /* fill out the matrices with squared residuals
       and lags of same */

    for (i=0; i<k; i++) {
        for (t=0; t<T; t++) {
            s = t + order;
            if (i == 0) {
                x = u[s];
                gretl_vector_set(y, t, x * x);
                gretl_matrix_set(X, t, i, 1.0);
            } else {
                x = u[s - i];
                gretl_matrix_set(X, t, i, x * x);
            }
        }
    }

    err = gretl_matrix_ols(y, X, b, V, null, ps2);

    if (!err) {
        rsq = gretl_matrix_r_squared(y, X, b, &err);
    }

    if (!err) {
        arch_test_save_or_print(b, V, T, order, rsq, pmod,
                                opt, prn);
    }

 bailout:

    gretl_matrix_free(X);
    gretl_matrix_free(y);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    return err;
}

/**
 * arch_test:
 * @pmod: model to be tested.
 * @order: lag order for ARCH process.
 * @dset: dataset struct.
 * @opt: if flags include OPT_S, save test results to model;
 * if OPT_Q, be less verbose; if OPT_I, be silent.
 * @prn: gretl printing struct.
 *
 * Tests @pmod for AutoRegressive Conditional Heteroskedasticity.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int arch_test(Model * pmod, int order, const(DatasetRaw) * dset,
              Gretlopt opt, Prn * prn)
{
    int err;

    if (pmod.missmask != null) {
        err = GretlError.missdata;
    } else {
        const(double) * u = pmod.uhat + pmod.t1;

        if (order == 0) {
            /* use data frequency as default lag order */
            order = dset.pd;
        }

        err = real_arch_test(u, pmod.nobs, order, pmod,
                             opt, prn);
    }

    return err;
}

int array_arch_test(const(double) * u, int n, int order,
                    Gretlopt opt, Prn * prn)
{
    return real_arch_test(u, n, order, null, opt, prn);
}

/**
 * arch_model:
 * @list: dependent variable plus list of regressors.
 * @order: lag order for ARCH process.
 * @dset: dataset struct.
 * @opt: may contain OPT_O to print covariance matrix.
 *
 * Estimate the model given in @list via weighted least squares,
 * with the weights based on the predicted error variances from
 * an auxiliary regression of the squared residuals on their lagged
 * values.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model arch_model(const(int) * list, int order, DatasetRaw * dset,
                 Gretlopt opt)
{
    Model amod;
    int * wlist = null, alist = null;
    int T = sample_size(dset);
    int oldv = dset.v;
    int i, t, nwt, k, n = dset.n;
    double * a = null;
    double * se = null;
    double xx;

    gretl_error_clear();
    gretl_model_init(&amod, dset);

    if (order == 0) {
        /* use data frequency as default lag order */
        order = dset.pd;
    }

    if (order < 1 || order > T - list[0]) {
        amod.errcode = GretlError.unspec;
        gretl_errmsg_sprintf(_("Invalid lag order for arch (%d)").ptr, order);
        return amod;
    }

    if (dataset_add_series(dset, order + 1)) {
        amod.errcode = GretlError.alloc;
    } else {
        alist = gretl_list_new(order + 2);
        if (alist == null) {
            amod.errcode = GretlError.alloc;
        }
    }

    if (amod.errcode) {
        goto bailout;
    }

    /* start list for aux regression */
    alist[1] = dset.v - order - 1;
    alist[2] = 0;

    /* run initial OLS and get squared residuals */
    amod = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_M);
    if (amod.errcode) {
        goto bailout;
    }

    k = dset.v - order - 1;
    strcpy(dset.varname[k], "utsq");
    for (t=0; t<n; t++) {
        dset.Z[k][t] = NADBL;
    }
    for (t=amod.t1; t<=amod.t2; t++) {
        xx = amod.uhat[t];
        dset.Z[k][t] = xx * xx;
    }
    /* also lags of squared resids */
    for (i=1; i<=order; i++) {
        k =  dset.v - order + i - 1;
        alist[i+2] = k;
        sprintf(dset.varname[k], "utsq_%d", i);
        for (t=0; t<n; t++) {
            dset.Z[k][t] = NADBL;
        }
        for (t=amod.t1+i; t<=amod.t2; t++) {
            dset.Z[k][t] = dset.Z[alist[1]][t-i];
        }
    }

    /* run auxiliary regression */
    clear_model(&amod);
    amod = lsq(alist, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (amod.errcode) {
        goto bailout;
    }

    /* steal the coefficients for reference */
    a = amod.coeff;
    amod.coeff = null;
    se = amod.sderr;
    amod.sderr = null;

    /* prepare for weighted estimation */
    wlist = gretl_list_new(list[0] + 1);

    if (wlist == null) {
        amod.errcode = GretlError.alloc;
    } else {
        /* construct the weight variable */
        nwt = wlist[1] = dset.v - 1;
        strcpy(dset.varname[nwt], "1/sigma");

        for (i=2; i<=wlist[0]; i++) {
            wlist[i] = list[i-1];
        }

        k = dset.v - order - 1;

        for (t=amod.t1; t<=amod.t2; t++) {
            xx = amod.yhat[t];
            if (xx <= 0.0) {
                xx = dset.Z[k][t];
            }
            dset.Z[nwt][t] = 1.0 / xx; /* FIXME is this right? */
        }

        clear_model(&amod);
        amod = lsq(wlist, dset, GretlCmdIndex.wls, Gretlopt.OPT_NONE);
        amod.ci = GretlCmdIndex.arch;

        if (!amod.errcode) {
            gretl_model_set_int(&amod, "arch_order", order);
            gretl_model_set_data(&amod, "arch_coeff", a,
                                 GretlType.double_array,
                                 (order + 1) * double.sizeof);
            gretl_model_set_data(&amod, "arch_sderr", se,
                                 GretlType.double_array,
                                 (order + 1) * double.sizeof);
            amod.xlist = gretl_list_sublist(amod.list, 2, -1);
        }
    }

 bailout:

    if (alist != null) free(alist);
    if (wlist != null) free(wlist);

    dataset_drop_last_variables(dset, dset.v - oldv);

    return amod;
}

/**
 * lad_model:
 * @list: dependent variable plus list of regressors.
 * @dset: dataset struct.
 * @opt: may include OPT_Q for quiet operation.
 *
 * Estimate the model given in @list using the method of Least
 * Absolute Deviation (LAD).
 *
 * Returns: a #MODEL struct containing the estimates.
 */

Model lad_model(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    Model mod;
    LadDriver lad_driver;

    /* run an initial OLS to "set the model up" and check for errors.
       the lad_driver function will overwrite the coefficients etc.
    */
    mod = lsq(list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (mod.errcode) {
        return mod;
    }

    lad_driver = cast(typeof(lad_driver)) get_plugin_function("lad_driver");
    if (lad_driver == null) {
        mod.errcode = GretlError.fopen;
        return mod;
    }

    (*lad_driver)(&mod, dset, opt);
    set_model_id(&mod, opt);

    return mod;
}

/**
 * quantreg:
 * @tau: vector containing one or more quantile values, in the range
 * 0.01 to 0.99.
 * @list: model specification: dependent var and regressors.
 * @dset: dataset struct.
 * @opt: may contain OPT_R for robust standard errors,
 * OPT_I to produce confidence intervals.
 * @prn: gretl printing struct.
 *
 * Estimate the model given in @list using the method of
 * quantile regression.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model quantreg(const(GretlMatrixRaw) * tau, const(int) * list,
               DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model mod;
    RqDriver rq_driver;
    Gretlopt olsopt = (Gretlopt.OPT_A | Gretlopt.OPT_M);

    /* Run an initial OLS to "set the model up" and check for errors.
       the driver function will overwrite the coefficients, etc.
       For now we make life easier by rejecting within-sample missing
       values (OPT_M).
    */

    if (opt & Gretlopt.OPT_R) {
        olsopt |= Gretlopt.OPT_R;
    }

    mod = lsq(list, dset, GretlCmdIndex.ols, olsopt);

    if (mod.errcode) {
        return mod;
    }

    rq_driver = cast(typeof(rq_driver)) get_plugin_function("rq_driver");

    if (rq_driver == null) {
        mod.errcode = GretlError.fopen;
        return mod;
    }

    (*rq_driver)(tau, &mod, dset, opt, prn);
    set_model_id(&mod, opt);

    return mod;
}

/* Replaces glib_grab_output(), which used g_spawn_sync() to run
   @prog (searched for in PATH, no arguments) and capture its stdout.
   Here that is done with std.process.execute. As in the C code, *sout
   is left null if the program could not be run or exited with a
   non-zero status; otherwise it is a malloc'd, NUL-terminated copy
   of the program's output.
*/

private void glib_grab_output(const(char) * prog,
                              char ** sout)
{
    import std.process : execute, ProcessException;
    import std.string : fromStringz;

    *sout = null;

    try {
        auto res = execute([cast(string) fromStringz(prog)]);

        if (res.status == 0) {
            char * s = cast(char *) malloc(res.output.length + 1);

            if (s != null) {
                memcpy(s, res.output.ptr, res.output.length);
                s[res.output.length] = '\0';
                *sout = s;
            }
        }
    } catch (ProcessException e) {
        /* could not spawn the program */
    }
}

/*
 * get_x13as_maxpd:
 *
 * Retrieve the highest data frequency handled by X-13ARIMA,
 * which may vary depending on how the program was built.
 * This may be relevant when executing the gretl arima
 * command with the option to use X-13ARIMA.
 */

int get_x13as_maxpd()
{
    __gshared static int n;

    if (n == 0) {
        const(char) * x12a = gretl_x12_arima();
        char * sout = null;

        glib_grab_output(x12a, &sout);
        if (sout != null) {
            char * p = strstr(sout, "PSP = ");

            if (p != null) {
                n = atoi(p + 6);
            }
            free(sout);
        }
        if (n <= 0) {
            n = 12;
        }
    }

    return n;
}

/**
 * arma:
 * @list: AR and MA orders, dependent variable, and any exogenous
 * regressors.
 * @pqlags: list giving specific non-seasonal AR/MA lags (or NULL).
 * @dset: dataset struct.
 * @opt: option flags.
 * @PRN: for printing details of iterations (or NULL).
 *
 * Calculates ARMA estimates. By default the estimator is
 * exact ML, implemented via gretl's Kalman filter in
 * conjunction with the BFGS maximizer. Other options:
 * if @opt includes OPT_L we use L-BFGS-B rather than standard
 * BFGS; OPT_C (incompatible with OPT_L) calls for use of
 * conditional ML (BHHH algorithm); and OPT_X requests estimation via
 * X-12-ARIMA rather than native code.
 *
 * If @pqlags is non-NULL it should take the form of two
 * gretl sub-lists joined by #LISTSEP: the first sub-list holds
 * a set of AR lags and the second a list of MA lags. If only
 * AR lags are specified, a single list may be given; if only
 * MA lags are specified, use #LISTSEP with a null first sub-list.
 *
 * If the model specification does not include a constant
 * this is added automatically, unless @opt includes OPT_N
 * ("no constant").
 *
 * When the estimation method is native exact ML, two (mutually
 * exclusive) flags in @opt may be used to control the estimator
 * of the covariance matrix: OPT_G specifies use of the outer
 * product of the gradient (OPG), while OPT_H specifies use of
 * the (numerical) Hessian. If neither of these flags is given, the
 * default is to use the Hessian by preference, but to fall back
 * to OPG if computation of the numerical Hessian fails. These
 * flags are ignored if estimation is not via native exact ML.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model arma(const(int) * list, const(int) * pqlags,
           DatasetRaw * dset, Gretlopt opt,
           Prn * prn)
{
    ArmaModelFn arma_model;
    ArmaX12ModelFn arma_x12_model;
    Model armod;
    int err = 0;

    gretl_model_init(&armod, dset);

    /* --opg and --hessian cannot be combined */
    err = incompatible_options(opt, Gretlopt.OPT_G | Gretlopt.OPT_H);
    if (!err) {
        /* --x-13arima cannot be combined with --robust, --stdx or
           --lagselect. (Note that --stdx is not documented.)
        */
        err = options_incompatible_with(opt, Gretlopt.OPT_X,
                                        Gretlopt.OPT_R | Gretlopt.OPT_S | Gretlopt.OPT_Z);
    }
    if (err) {
        armod.errcode = err;
        return armod;
    }

    if (opt & Gretlopt.OPT_X) {
        int pdmax = get_x13as_maxpd();

        if ((dset.t2 - dset.t1 + 1) > pdmax * 50) {
            gretl_errmsg_sprintf(_("X-13ARIMA can't handle more than %d observations.\n"
                                   ~ "Please select a smaller sample.").ptr, pdmax * 50);
            armod.errcode = GretlError.data;
            return armod;
        }
        arma_x12_model = cast(typeof(arma_x12_model)) get_plugin_function("arma_x12_model");
        if (arma_x12_model == null) {
            err = GretlError.fopen;
        } else {
            armod = (*arma_x12_model)(list, pqlags, dset, pdmax, opt, prn);
        }
    } else if (opt & Gretlopt.OPT_Z) {
        ArmaSelectFn arma_select;

        arma_select = cast(typeof(arma_select)) get_plugin_function("arma_select");
        if (arma_select == null) {
            err = GretlError.fopen;
        } else {
            err = (*arma_select)(list, pqlags, dset, opt, prn);
        }
    } else {
        arma_model = cast(typeof(arma_model)) get_plugin_function("arma_model");
        if (arma_model == null) {
            err = GretlError.fopen;
        } else {
            armod = (*arma_model)(list, pqlags, dset, opt, prn);
        }
    }

    if (err) {
        armod.errcode = err;
        return armod;
    }

    if (!(opt & Gretlopt.OPT_Z)) {
        set_model_id(&armod, opt);
    }

    return armod;
}

/**
 * garch:
 * @list: dependent variable plus arch and garch orders.
 * @dset: dataset struct.
 * @opt: can specify robust standard errors and VCV.
 * @prn: for printing details of iterations (or NULL).
 *
 * Calculate GARCH estimates.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model garch(const(int) * list, DatasetRaw * dset, Gretlopt opt,
            Prn * prn)
{
    Model mod;
    Prn * myprn;
    GarchModelFn garch_model;

    gretl_error_clear();

    garch_model = cast(typeof(garch_model)) get_plugin_function("garch_model");

    if (garch_model == null) {
        gretl_model_init(&mod, dset);
        mod.errcode = GretlError.fopen;
        return mod;
    }

    if (opt & Gretlopt.OPT_V) {
        myprn = prn;
    } else {
        myprn = null;
    }

    mod = (*garch_model)(list, dset, myprn, opt);
    set_model_id(&mod, opt);

    return mod;
}

/**
 * mp_ols:
 * @list: specification of variables to use.
 * @dset: pointer to dataset struct.
 * @opt: may include OPT_Q for quiet operation, OPT_S for
 * simple output.
 *
 * Estimate an OLS model using multiple-precision arithmetic
 * via the GMP library.
 *
 * Returns: a #MODEL struct containing the estimates.
 */

Model mp_ols(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    MplsqFn mplsq;
    Model mpmod;

    gretl_model_init(&mpmod, dset);

    mplsq = cast(typeof(mplsq)) get_plugin_function("mplsq");
    if (mplsq == null) {
        mpmod.errcode = 1;
        return mpmod;
    }

    if (gretl_list_has_separator(list)) {
        int * base = null;
        int * poly = null;

        mpmod.errcode = gretl_list_split_on_separator(list, &base, &poly);
        if (mpmod.errcode == 0 && (base == null || poly == null)) {
            mpmod.errcode = GretlError.args;
        } else {
            mpmod.errcode = (*mplsq)(base, poly, null, dset,
                                     &mpmod, Gretlopt.OPT_X);
        }
        free(base);
        free(poly);
    } else {
        mpmod.errcode = (*mplsq)(list, null, null, dset,
                                 &mpmod, Gretlopt.OPT_X);
    }

    set_model_id(&mpmod, opt);

    return mpmod;
}

private int check_panel_options(Gretlopt opt)
{
    if ((opt & Gretlopt.OPT_U) && (opt & Gretlopt.OPT_H)) {
        /* can't specify random effects + weighted least squares */
        return GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_T) && !(opt & Gretlopt.OPT_H)) {
        /* iterate option requires weighted least squares option */
        return GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_E) && !(opt & Gretlopt.OPT_U)) {
        /* the Nerlove option requires random effects */
        return GretlError.badopt;
    } else if (incompatible_options(opt, Gretlopt.OPT_B | Gretlopt.OPT_U | Gretlopt.OPT_P)) {
        /* mutually exclusive estimator requests */
        return GretlError.badopt;
    } else if (opt & Gretlopt.OPT_J) {
        /* jackknife option not OK for panel data, at present */
        gretl_errmsg_set("The --jackknife option is not supported for panel data".ptr);
        return GretlError.badopt;
    }

    return 0;
}

/**
 * panel_model:
 * @list: regression list (dependent variable plus independent
 * variables).
 * @dset: dataset struct.
 * @opt: can include OPT_Q (quiet estimation), OPT_U (random
 * effects model), OPT_H (weights based on the error variance
 * for the respective cross-sectional units), OPT_I (iterate,
 * only available in conjunction with OPT_H).
 * @prn: printing struct (or NULL).
 *
 * Calculate estimates for a panel dataset, using fixed
 * effects (the default), random effects, or weighted
 * least squares based on the respective variances for the
 * cross-sectional units.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model panel_model(const(int) * list, DatasetRaw * dset,
                  Gretlopt opt, Prn * prn)
{
    Model mod;

    if (check_panel_options(opt)) {
        gretl_model_init(&mod, dset);
        mod.errcode = GretlError.badopt;
    } else if (opt & Gretlopt.OPT_H) {
        mod = panel_wls_by_unit(list, dset, opt, prn);
    } else {
        mod = real_panel_model(list, dset, opt, prn);
    }

    return mod;
}

/**
 * ivreg:
 * @list: regression list; see the documentation for the "tsls"
 * command.
 * @dset: dataset struct.
 * @opt: option flags.
 *
 * Calculate IV estimates for a linear model, by default via two-stage
 * least squares. The option flags can include OPT_Q for quiet
 * operation, and either OPT_L to specify estimation via Limited
 * Information Maximum Likelihood or OPT_G for estimation via
 * Generalized method of Moments.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model ivreg(const(int) * list, DatasetRaw * dset, Gretlopt opt)
{
    Model mod;
    int err;

    gretl_error_clear();

    /* can't have both LIML and GMM options */
    err = incompatible_options(opt, Gretlopt.OPT_G | Gretlopt.OPT_L);

    if (!err) {
        /* two-step, iterate and weights options are GMM-only */
        err = option_prereq_missing(opt, Gretlopt.OPT_T | Gretlopt.OPT_I | Gretlopt.OPT_H,
                                    Gretlopt.OPT_G);
    }

    if (err) {
        gretl_model_init(&mod, dset);
        mod.errcode = err;
        return mod;
    }

    if (opt & Gretlopt.OPT_L) {
        mod = single_equation_liml(list, dset, opt);
    } else if (opt & Gretlopt.OPT_G) {
        mod = ivreg_via_gmm(list, dset, opt);
    } else {
        mod = tsls(list, dset, opt);
    }

    return mod;
}

/**
 * dpd_model:
 * @list: regression list.
 * @laglist: list of specific lags of the dependent variable, or
 * NULL.
 * @ispec: may contain additional instrument specification.
 * @dset: dataset struct.
 * @opt: may include... see options.c.
 * @prn: printing struct.
 *
 * Implements the "dpanel" command.
 *
 * Returns: a #MODEL struct, containing the estimates.
 */

Model dpd_model(const(int) * list, const(int) * laglist,
                const(char) * ispec, const(DatasetRaw) * dset,
                Gretlopt opt, Prn * prn)
{
    DpdEstimateFn dpd_estimate;
    Model mod;

    gretl_model_init(&mod, dset);

    dpd_estimate = cast(typeof(dpd_estimate)) get_plugin_function("dpd_estimate");
    if (dpd_estimate == null) {
        mod.errcode = 1;
        return mod;
    }

    mod = (*dpd_estimate)(list, laglist, ispec, dset, opt, prn);
    set_model_id(&mod, opt);

    return mod;
}

/**
 * anova:
 * @list: must contain the response and treatment variables.
 * @dset: dataset struct.
 * @opt: unused at present.
 * @prn: printing struct.
 *
 * Does one-way or two-way Analysis of Variance (prints table and F-test).
 *
 * Returns: 0 on success, non-zero on failure.
*/

int anova(const(int) * list, const(DatasetRaw) * dset,
          Gretlopt opt, Prn * prn)
{
    GretlAnovaFn gretl_anova;

    gretl_anova = cast(typeof(gretl_anova)) get_plugin_function("gretl_anova");
    if (gretl_anova == null) {
        return 1;
    }

    return (*gretl_anova)(list, dset, opt, prn);
}
