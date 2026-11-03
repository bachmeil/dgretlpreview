/* Port of plugin/arma.c, with plugin/arma_common.c (which arma.c includes): estimation of ARMA and ARIMA models (arma_model), and
   selection of the orders by information criteria (arma_select, gui_arma_select).

   In the C the files kalman_arma.c, as197.c, as154.c and as_driver.c are also included by arma.c; kalman_arma.c is the module
   gretl.plugin.kalman_arma (ported), the AS 154 / AS 197 driver is gretl.plugin.as_driver (ported; AS 154 itself, as154.c, is not in the tree).
   The initialization functions of arma_init.c are in gretl.plugin.arma_init; bhhh_arma of arma_cml.c is in gretl.plugin.arma_cml.

   The functions that kalman_arma.c uses from arma.c (form_arma_y_vector, form_arma_X_matrix, arma_QML_vcv, arma_OPG_vcv,
   arma_use_opg, arima_ydiff_only) are extern (C) and public here, as are maybe_correct_MA and arma_model_add_roots, flip_poly and
   so on (the prototypes of arma_priv.h). Debug output under ARMA_DEBUG, SAMPLE_DEBUG and the like (all off in the C) is left out.
   The #ifdef X12A_CODE branch (not defined in the plugin) is dropped. */

module gretl.plugin.arma;

import cephes.polrt : polrt;
import gretl.dataset : dataset_add_series, dataset_drop_last_variables;
import gretl.describe : gretl_corr_rsq, gretl_mean, gretl_moments, gretl_stddev;
import gretl.gretl_bfgs : numerical_score_matrix;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_list_copy, gretl_list_delete_at_pos, gretl_list_has_separator, gretl_list_new, gretl_list_separator_position, gretl_list_split_on_separator;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_invert_symmetric_matrix, gretl_matrix_XTX_new, gretl_matrix_alloc, gretl_matrix_divide_by_scalar, gretl_matrix_free, gretl_matrix_init, gretl_matrix_multiply_by_scalar, gretl_matrix_ols, gretl_matrix_set_colnames, gretl_matrix_set_t1, gretl_matrix_set_t2, gretl_symmetric_matrix_rcond;
import gretl.gretl_model : arima_delta_coeffs, clear_model, clear_model_xpx, gretl_model_add_QML_vcv, gretl_model_add_arma_varnames, gretl_model_allocate_storage, gretl_model_init, gretl_model_set_data, gretl_model_set_double, gretl_model_set_int, gretl_model_set_list_as_data, gretl_model_set_string_as_data, gretl_model_smpl_init, gretl_model_write_vcv, mle_criteria;
import gretl.gretl_modelh : ArmaFlags;
import gretl.gretl_prn : pprintf, pputc, pputs;
import gretl.gretl_utils : doubles_array_free, record_matrix_test_result, true_const;
import gretl.libset : get_initvals, libset_get_int, n_initvals;
import gretl.matrix_extra : gretl_matrix_data_subset, gretl_vector_from_series;
import gretl.missingh : NADBL, na;
import gretl.modelprint : arma_extra_info;
import gretl.options : incompatible_options, options_incompatible_with;
import gretl.printout : dashline;
import gretl.strutils : gretl_strdup, strings_array_new;

import gretl.gretl_prn : Prn;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc, qsort;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strchr;

import gretl.plugin.arma_priv;
import gretl.plugin.as_driver : as154_llt_callback, as197_llt_callback, as_arma;
import gretl.plugin.arma_cml : bhhh_arma;
import gretl.plugin.arma_init : ar_arma_init, arma_by_ls, arma_by_simple_ols, flip_poly, hr_arma_init, maybe_set_yscale;
import gretl.matrix_extrah : MMissingCode;
import cephes.mconfh : CephesCmplx = Cmplx;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.libseth : SetKey;
import gretl.gretl_modelh : MLVCVType;
import gretl.plugin.kalman_arma : kalman_arma, kalman_arma_ll, kalman_arma_llt_callback;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;










enum MAX_ARMA_ORDER = 128;
enum MAX_ARIMA_DIFF = 2;


private void arma_info_init(ArmaInfo * ainfo, Gretlopt opt, const(int) * pqspec, const(DatasetRaw) * dset)
{
    ainfo.yno = 0;
    ainfo.flags = 0;
    ainfo.pflags = 0;
    ainfo.init_ = cast(IniMethod) 0;
    ainfo.alist = null;

    if (opt & Gretlopt.OPT_X) {
	/* we got --x-12-arima */
	ainfo.flags |= ArmaFlags.x12a;
    }

    if (!(opt & Gretlopt.OPT_C)) {
	/* we didn't get --conditional */
	ainfo.flags |= ArmaFlags.exact;
    }

    ainfo.ll = NADBL;

    ainfo.pqspec = pqspec;
    ainfo.pmask = null;
    ainfo.qmask = null;

    ainfo.p = 0;
    ainfo.d = 0;
    ainfo.q = 0;
    ainfo.P = 0;
    ainfo.D = 0;
    ainfo.Q = 0;

    ainfo.np = 0;
    ainfo.nq = 0;

    ainfo.maxlag = 0;
    ainfo.ifc = 0;
    ainfo.nexo = 0;
    ainfo.nc = 0;

    ainfo.t1 = dset.t1;
    ainfo.t2 = dset.t2;
    ainfo.pd = dset.pd;
    ainfo.T = 0;
    ainfo.r0 = 0;

    ainfo.fncount = 0;
    ainfo.grcount = 0;

    ainfo.y = null;
    ainfo.e = null;
    ainfo.Z = null;
    ainfo.yscale = 1.0;
    ainfo.yshift = 0.0;

    ainfo.xlist = null;
    ainfo.misslist = null;
    ainfo.xstats = null;
    ainfo.dX = null;
    ainfo.G = null;
    ainfo.V = null;

    ainfo.n_aux = 0;
    ainfo.aux = null;

    ainfo.prn = null;
}

private void arma_info_cleanup(ArmaInfo * ainfo)
{
    free(ainfo.alist);
    free(ainfo.pmask);
    free(ainfo.qmask);
    free(ainfo.e);
    free(ainfo.Z);
    free(ainfo.xlist);
    free(ainfo.misslist);

    gretl_matrix_free(ainfo.xstats);
    gretl_matrix_free(ainfo.dX);
    gretl_matrix_free(ainfo.G);
    gretl_matrix_free(ainfo.V);

    if (arima_ydiff(ainfo)) {
	free(ainfo.y);
    }

    doubles_array_free(ainfo.aux, ainfo.n_aux);
}

enum ArmaMask {
    AR_MASK,
    MA_MASK
}

/* Create a mask for skipping certain intermediate lags,
   AR or MA.  This function also sets ainfo.np and ainfo.nq,
   which record the actual number of non-seasonal AR and MA
   lags used.
*/

private char * mask_from_list(const(int) * list, ArmaInfo * ainfo, int m, int * err)
{
    int mlen = (m == ArmaMask.AR_MASK)? ainfo.p : ainfo.q;
    int nv = 0;
    int nmax = 0;
    char * mask;
    int i;
    int k;

    mask = cast(char *) malloc(mlen + 1);
    if (mask == null) {
	*err = GretlError.alloc;
	return null;
    }

    for (i=0; i<mlen; i++) {
	mask[i] = '0';
    }
    mask[mlen] = '\0';

    for (i=1; i<=list[0]; i++) {
	k = list[i];
	if (k > 0) {
	    mask[k-1] = '1';
	    nv++;
	    if (k > nmax) {
		nmax = k;
	    }
	}
    }

    if (m == ArmaMask.AR_MASK) {
	ainfo.p = nmax;
	ainfo.np = nv;
    } else {
	ainfo.q = nmax;
	ainfo.nq = nv;
    }

    if (nv == 0) {
	free(mask);
	mask = null;
    }

    return mask;
}

private int arma_make_masks(ArmaInfo * ainfo)
{
    int * plist = null;
    int * qlist = null;
    int err = 0;

    if (ainfo.pqspec != null) {
	if (gretl_list_has_separator(ainfo.pqspec)) {
	    gretl_list_split_on_separator(ainfo.pqspec, &plist, &qlist);
	} else {
	    plist = gretl_list_copy(ainfo.pqspec);
	}
    }

    if (ainfo.p > 0) {
	ainfo.np = ainfo.p;
	if (plist != null && plist[0] > 0) {
	    ainfo.pmask = mask_from_list(plist, ainfo, ArmaMask.AR_MASK, &err);
	}
    }

    if (ainfo.q > 0 && !err) {
	ainfo.nq = ainfo.q;
	if (qlist != null && qlist[0] > 0) {
	    ainfo.qmask = mask_from_list(qlist, ainfo, ArmaMask.MA_MASK, &err);
	}
    }

    free(plist);
    free(qlist);

    return err;
}

extern (C) int arma_list_y_position(ArmaInfo * ainfo)
{
    int ypos;

    if (arma_is_arima(ainfo)) {
	ypos = arma_has_seasonal(ainfo) ? 9 : 5;
    } else {
	ypos = arma_has_seasonal(ainfo) ? 7 : 4;
    }

    return ypos;
}

private int arima_integrate(double * dx, const(double) * x, int t1, int t2, int d, int D, int s)
{
    double * ix;
    int * c;
    int k = d + s * D;
    int i;
    int t;

    ix = cast(double *) malloc((t2 + 1) * (*ix).sizeof);
    if (ix == null) {
	return GretlError.alloc;
    }

    c = arima_delta_coeffs(d, D, s);
    if (c == null) {
	free(ix);
	return GretlError.alloc;
    }

    for (t=0; t<t1; t++) {
	ix[t] = 0.0;
    }

    for (t=t1; t<=t2; t++) {
	ix[t] = dx[t];
	for (i=0; i<k; i++) {
	    if (c[i] != 0) {
		ix[t] += c[i] * x[t-i-1];
	    }
	}
    }

    /* transcribe integrated result back into "dx" */
    for (t=0; t<=t2; t++) {
	if (t < t1) {
	    dx[t] = NADBL;
	} else {
	    dx[t] = ix[t];
	}
    }

    free(ix);
    free(c);

    return 0;
}

private void ainfo_data_to_model(ArmaInfo * ainfo, Model * pmod)
{
    pmod.ifc = ainfo.ifc;
    pmod.dfn = ainfo.nc - pmod.ifc;
    pmod.dfd = pmod.nobs - pmod.dfn;
    pmod.ncoeff = ainfo.nc;

    if (arma_has_seasonal(ainfo)) {
	gretl_model_set_int(pmod, "arma_P", ainfo.P);
	gretl_model_set_int(pmod, "arma_Q", ainfo.Q);
	gretl_model_set_int(pmod, "arma_pd", ainfo.pd);
    }

    if (ainfo.d > 0 || ainfo.D > 0) {
	gretl_model_set_int(pmod, "arima_d", ainfo.d);
	gretl_model_set_int(pmod, "arima_D", ainfo.D);
    }

    if (ainfo.nexo > 0) {
	gretl_model_set_int(pmod, "armax", 1);
    }

    if (ainfo.pmask != null) {
	gretl_model_set_string_as_data(pmod, "pmask",
				       gretl_strdup(ainfo.pmask));
    }

    if (ainfo.qmask != null) {
	gretl_model_set_string_as_data(pmod, "qmask",
				       gretl_strdup(ainfo.qmask));
    }
}

private void arma_depvar_stats(Model * pmod, ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    if (arma_is_arima(ainfo) && !arima_ydiff(ainfo)) {
	/* calculate differenced y for stats */
	int d = ainfo.d;
	int D = ainfo.D;
	int T = pmod.t2 - pmod.t1 + 1;
	double * dy = cast(double *) malloc(T * double.sizeof);
	int *delta = arima_delta_coeffs(d, D, ainfo.pd);

	if (dy != null && delta != null) {
	    int k = d + ainfo.pd * D;

	    real_arima_difference_series(dy, dset.Z[ainfo.yno],
					 pmod.t1, pmod.t2, delta, k);
	    pmod.ybar = gretl_mean(0, T - 1, dy);
	    pmod.sdy = gretl_stddev(0, T - 1, dy);
	}
	free(dy);
	free(delta);
    } else {
	pmod.ybar = gretl_mean(pmod.t1, pmod.t2, ainfo.y);
	pmod.sdy = gretl_stddev(pmod.t1, pmod.t2, ainfo.y);
    }
}

private void handle_null_model(Model * pmod, ArmaInfo * ainfo)
{
    int full_n = pmod.full_n;

    pmod.ncoeff = 1;
    pmod.full_n = 0;
    pmod.errcode = gretl_model_allocate_storage(pmod);
    pmod.full_n = full_n;

    if (!pmod.errcode) {
	gretl_model_set_int(pmod, "null-model", 1);
	pmod.coeff[0] = 0.0;
	pmod.sigma = pmod.sdy;
    }
}


/* write the various statistics from ARMA estimation into
   a gretl MODEL struct */

extern (C) void write_arma_model_stats(Model * pmod, ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    double mean_error;
    int do_criteria = 1;
    int t;

    pmod.ci = GretlCmdIndex.arma;
    ainfo_data_to_model(ainfo, pmod);

    free(pmod.list);
    pmod.list = gretl_list_copy(ainfo.alist);

    if (!arma_least_squares(ainfo)) {
	arma_depvar_stats(pmod, ainfo, dset);
    }

    mean_error = pmod.ess = 0.0;

    for (t=pmod.t1; t<=pmod.t2; t++) {
	if (!na(ainfo.y[t]) && !na(pmod.uhat[t])) {
	    pmod.yhat[t] = ainfo.y[t] - pmod.uhat[t];
	    pmod.ess += pmod.uhat[t] * pmod.uhat[t];
	    mean_error += pmod.uhat[t];
	}
    }

    if (arma_is_arima(ainfo) && arima_ydiff(ainfo)) {
	arima_integrate(pmod.yhat, dset.Z[ainfo.yno],
			pmod.t1, pmod.t2,
			ainfo.d, ainfo.D, ainfo.pd);
    }

    mean_error /= pmod.nobs;
    if (arma_least_squares(ainfo) && pmod.ifc &&
	mean_error < 1.0e-15) {
	mean_error = 0.0;
    }
    gretl_model_set_double(pmod, "mean_error", mean_error);

    if (na(pmod.sigma)) {
	/* in x12a or native exact cases this is already done */
	pmod.sigma = sqrt(pmod.ess / pmod.nobs);
    }

    /* R-squared as suggested by J\"org Breitung */
    pmod.rsq = gretl_corr_rsq(pmod.t1, pmod.t2, dset.Z[ainfo.yno], pmod.yhat);
    pmod.adjrsq = 1.0 - ((1.0 - pmod.rsq) * (pmod.t2 - pmod.t1) /
			  cast(double) pmod.dfd);
    pmod.fstt = pmod.chisq = NADBL;
    pmod.tss = NADBL;

    if (arma_least_squares(ainfo)) {
	/* not applicable */
	do_criteria = 0;
    } else if (arma_by_x12a(ainfo) && !na(pmod.criterion[InfoCriterion.aic])) {
	/* already given by x12a */
	do_criteria = 0;
    }

    if (do_criteria) {
	mle_criteria(pmod, 1);
    }

    if (!pmod.errcode && pmod.ncoeff == 0) {
	handle_null_model(pmod, ainfo);
    }

    if (!pmod.errcode) {
	gretl_model_add_arma_varnames(pmod, dset, ainfo.yno,
				      ainfo.p, ainfo.q,
				      ainfo.pmask, ainfo.qmask,
				      ainfo.P, ainfo.Q,
				      ainfo.nexo);
    }
}

private void calc_max_lag(ArmaInfo * ainfo)
{
    if (arma_exact_ml(ainfo)) {
	ainfo.maxlag = ainfo.d + ainfo.D * ainfo.pd;
    } else {
	/* conditional ML */
	int pmax = ainfo.p + ainfo.P * ainfo.pd;
	int dmax = ainfo.d + ainfo.D * ainfo.pd;

	ainfo.maxlag = pmax + dmax;
    }

}

private int arma_adjust_sample(ArmaInfo * ainfo, const(DatasetRaw) * dset, int * missv, int * misst)
{
    int * list = ainfo.alist;
    int ypos = arma_list_y_position(ainfo);
    int t0;
    int t1 = dset.t1;
    int t2 = dset.t2;
    int i;
    int vi;
    int vlmax;
    int k;
    int t;
    int missing;
    int err = 0;


    t0 = t1 - ainfo.maxlag;
    if (t0 < 0) {
	t1 -= t0;
    }

    /* list position of last var to check for lags */
    if (arma_xdiff(ainfo)) {
	vlmax = list[0];
    } else {
	vlmax = ypos;
    }

    /* advance the starting point if need be */
    for (t=t1; t<=t2; t++) {
	missing = 0;
	for (i=ypos; i<=list[0] && !missing; i++) {
	    vi = list[i];
	    if (na(dset.Z[vi][t])) {
		/* current value missing */
		missing = 1;
	    }
	    if (i <= vlmax) {
		for (k=1; k<=ainfo.maxlag && !missing; k++) {
		    if (na(dset.Z[vi][t-k])) {
			/* lagged value missing */
			missing = 1;
		    }
		}
	    }
	}
	if (missing) {
	    t1++;
	} else {
	    break;
	}
    }

    /* retard the ending point if need be */
    for (t=t2; t>=t1; t--) {
	missing = 0;
	for (i=ypos; i<=list[0] && !missing; i++) {
	    vi = list[i];
	    if (na(dset.Z[vi][t])) {
		missing = 1;
	    }
	}
	if (missing) {
	    t2--;
	} else {
	    break;
	}
    }

    if (t2 < t1) {
	gretl_errmsg_set(_("No usable data were found").ptr);
	return GretlError.missdata;
    }

    missing = 0;

    /* check for missing obs within the adjusted sample range */
    for (t=t1; t<t2; t++) {
	int tmiss = 0;

	for (i=ypos; i<=list[0]; i++) {
	    vi = list[i];
	    if (na(dset.Z[vi][t])) {
		if (missv != null && misst != null && *missv == 0) {
		    /* record info on first missing obs */
		    *missv = vi;
		    *misst = t + 1;
		}
		tmiss = 1;
	    }
	}
	if (tmiss) {
	    missing++;
	}
    }

    if (missing > 0 && !arma_na_ok(ainfo)) {
	err = GretlError.missdata;
    }

    if (!err) {
	ainfo.fullT = t2 - t1 + 1;
	ainfo.T = ainfo.fullT - missing;
	if (ainfo.T <= ainfo.nc) {
	    /* insufficient observations */
	    err = GretlError.df;
	}
    }

    if (!err) {
	ainfo.t1 = t1;
	ainfo.t2 = t2;
    }

    return err;
}

/* remove the intercept from list of regressors */

private int arma_remove_const(ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    int * list = ainfo.alist;
    int seasonal = arma_has_seasonal(ainfo);
    int diffs = arma_is_arima(ainfo);
    int xstart;
    int ret = 0;
    int i;
    int j;

    if (diffs) {
	xstart = (seasonal)? 10 : 6;
    } else {
	xstart = (seasonal)? 8 : 5;
    }

    for (i=xstart; i<=list[0]; i++) {
	if (list[i] == 0 || true_const(list[i], dset)) {
	    for (j=i; j<list[0]; j++) {
		list[j] = list[j+1];
	    }
	    list[0] -= 1;
	    ret = 1;
	    break;
	}
    }

    return ret;
}

private int check_arma_sep(ArmaInfo * ainfo, int sep1)
{
    int * list = ainfo.alist;
    int sep2 = (sep1 == 3)? 6 : 8;
    int i;
    int err = 0;

    for (i=sep1+1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    if (i == sep2) {
		/* there's a second list separator in the right place:
		   we've got a seasonal specification */
		set_arma_has_seasonal(ainfo);
	    } else {
		err = 1;
	    }
	}
    }

    if (!err && sep1 == 4) {
	/* check for apparent but not "real" arima spec */
	if (arma_has_seasonal(ainfo)) {
	    if (list[2] == 0 && list[6] == 0) {
		gretl_list_delete_at_pos(list, 2);
		gretl_list_delete_at_pos(list, 5);
		unset_arma_is_arima(ainfo);
	    }
	} else {
	    if (list[2] == 0) {
		gretl_list_delete_at_pos(list, 2);
		unset_arma_is_arima(ainfo);
	    }
	}
    }

    if (!err && arma_has_seasonal(ainfo) && ainfo.pd == 1) {
        gretl_errmsg_set("Invalid seasonal specification: data frequency = 1");
        err = GretlError.pdwrong;
    }

    return err;
}

private int arma_add_xlist(ArmaInfo * ainfo, int ypos)
{
    int i;
    int err = 0;

    ainfo.xlist = gretl_list_new(ainfo.nexo);

    if (ainfo.xlist == null) {
	err = GretlError.alloc;
    } else {
	for (i=1; i<=ainfo.nexo; i++) {
	    ainfo.xlist[i] = ainfo.alist[ypos + i];
	}
    }

    return err;
}

int count_arma_coeffs(const(ArmaInfo) * a)
{
    return a.ifc + a.np + a.nq + a.P + a.Q + a.nexo;
}

private int check_arma_list(ArmaInfo * ainfo, const(DatasetRaw) * dset, Gretlopt opt)
{
    int * list = ainfo.alist;
    int ypos = arma_has_seasonal(ainfo) ? 7 : 4;
    int armax = (list[0] > ypos);
    int hadconst = 0;
    int err = 0;

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	/* non-seasonal AR order out of bounds */
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARMA_ORDER) {
	/* non-seasonal MA order out of bounds */
	err = 1;
    }

    if (!err) {
	ainfo.p = list[1];
	ainfo.q = list[2];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 7) {
	    err = 1;
	} else if (list[4] < 0 || list[4] > MAX_ARMA_ORDER) {
	    /* seasonal AR order out of bounds */
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    /* seasonal MA order out of bounds */
	    err = 1;
	}
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo.P = list[4];
	ainfo.Q = list[5];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.0  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.0
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(ainfo, dset);
	}
	if ((opt & Gretlopt.OPT_N) || (armax && !hadconst)) {
	    ; /* no constant present */
	} else {
	    ainfo.ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command").ptr);
    } else {
	ainfo.nexo = list[0] - ypos;
	ainfo.nc = count_arma_coeffs(ainfo);
	ainfo.yno = list[ypos];
	if (ainfo.nexo > 0) {
	    err = arma_add_xlist(ainfo, ypos);
	}
    }

    return err;
}

private int check_arima_list(ArmaInfo * ainfo, const(DatasetRaw) * dset, Gretlopt opt)
{
    int * list = ainfo.alist;
    int ypos = arma_has_seasonal(ainfo) ? 9 : 5;
    int armax = (list[0] > ypos);
    int hadconst = 0;
    int err = 0;

    if (list[1] < 0 || list[1] > MAX_ARMA_ORDER) {
	err = 1;
    } else if (list[2] < 0 || list[2] > MAX_ARIMA_DIFF) {
	err = 1;
    } else if (list[3] < 0 || list[3] > MAX_ARMA_ORDER) {
	err = 1;
    }

    if (!err) {
	ainfo.p = list[1];
	ainfo.d = list[2];
	ainfo.q = list[3];
    }

    if (!err && arma_has_seasonal(ainfo)) {
	if (list[0] < 9) {
	    err = 1;
	} else if (list[5] < 0 || list[5] > MAX_ARMA_ORDER) {
	    err = 1;
	} else if (list[6] < 0 || list[6] > MAX_ARIMA_DIFF) {
	    err = 1;
	} else if (list[7] < 0 || list[7] > MAX_ARMA_ORDER) {
	    err = 1;
	}
    }

    if (!err && arma_has_seasonal(ainfo)) {
	ainfo.P = list[5];
	ainfo.D = list[6];
	ainfo.Q = list[7];
    }

    /* now that we have p and q we can check for masked lags */

    if (!err) {
	err = arma_make_masks(ainfo);
    }

    /* If there's an explicit constant in the list here, we'll remove
       it, since it is added implicitly later.  But if we're supplied
       with OPT_N (meaning: no intercept) we'll flag this by
       setting ifc = 0.0  Also, if the user gave an armax list
       (specifying regressors) we'll respect the absence of a constant
       from that list by setting ifc = 0.0
    */

    if (!err) {
	if (armax) {
	    hadconst = arma_remove_const(ainfo, dset);
	}
	if ((opt & Gretlopt.OPT_N) || (armax && !hadconst)) {
	    ;
	} else {
	    ainfo.ifc = 1;
	}
    }

    if (err) {
	gretl_errmsg_set(_("Error in arma command").ptr);
    } else {
	ainfo.nexo = list[0] - ypos;
	ainfo.nc = count_arma_coeffs(ainfo);
	ainfo.yno = list[ypos];
	if (ainfo.nexo > 0) {
	    err = arma_add_xlist(ainfo, ypos);
	}
    }

    return err;
}

private int arma_check_list(ArmaInfo * ainfo, const(DatasetRaw) * dset, Gretlopt opt)
{
    int * list = ainfo.alist;
    int sep1 = gretl_list_separator_position(list);
    int err = 0;

    if (sep1 == 3) {
	if (list[0] < 4) {
	    err = GretlError.parse;
	}
    } else if (sep1 == 4) {
	if (list[0] < 5) {
	    err = GretlError.parse;
	} else {
	    set_arma_is_arima(ainfo);
	}
    } else {
	err = GretlError.parse;
    }

    if (!err) {
	err = check_arma_sep(ainfo, sep1);
    }

    if (!err) {
	if (arma_is_arima(ainfo)) {
	    /* check for arima spec */
	    err = check_arima_list(ainfo, dset, opt);
	} else {
	    /* check for simple arma spec */
	    err = check_arma_list(ainfo, dset, opt);
	    /* catch null model */
	    if (ainfo.nc == 0) {
		err = GretlError.args;
	    }
	}
    }

    return err;
}

private void real_arima_difference_series(double * dx, const(double) * x, int t1, int t2, int * delta, int k)
{
    int i;
    int p;
    int t;
    int s = 0;

    for (t=t1; t<=t2; t++) {
	dx[s] = x[t];
	for (i=0; i<k && !na(dx[s]); i++) {
	    if (delta[i] != 0) {
		p = t - i - 1;
		if (p < 0 || na(x[p])) {
		    dx[s]  = NADBL;
		} else {
		    dx[s] -= delta[i] * x[p];
		}
	    }
	}
	s++;
    }
}

private int transcribe_extra_info(ArmaInfo * ainfo, Model * armod)
{
    int *ainfo_list = gretl_list_new(9);
    int err = 0;

    if (ainfo_list == null) {
	armod.errcode = err = GretlError.alloc;
    } else {
	/* wrap up a bunch of relevant integers */
	ainfo_list[1] = ainfo.p;
	ainfo_list[2] = ainfo.q;
	ainfo_list[3] = ainfo.P;
	ainfo_list[4] = ainfo.Q;
	ainfo_list[5] = ainfo.np;
	ainfo_list[6] = ainfo.nq;
	ainfo_list[7] = ainfo.d;
	ainfo_list[8] = ainfo.D;
	ainfo_list[9] = ainfo.pd;
	err = gretl_model_set_list_as_data(armod, "ainfo", ainfo_list);
    }

    if (!err && arma_xdiff(ainfo)) {
	gretl_model_set_int(armod, "xdiff", 1);
    }

    if (ainfo.ifc || ainfo.nexo > 0) {
	int i;

	armod.xlist = gretl_list_new(ainfo.ifc + ainfo.nexo);
	if (ainfo.ifc) {
	    armod.xlist[1] = 0;
	}
	for (i=0; i<ainfo.nexo; i++) {
	    armod.xlist[i+ainfo.ifc+1] = ainfo.xlist[i+1];
	}
    }

    return err;
}


/* Add to the ainfo struct a full-length series y holding
   the differenced version of the dependent variable.
   If the "xdiff" flag is set on ainfo, in addition
   create a matrix dX holding the differenced regressors;
   in that case the time-series length of dX depends on
   the @fullX flag -- if fullX = 0, this equals
   ainfo.T but if fullX = 0 it equals ainfo.t2 + 1.0
*/

extern (C) int arima_difference(ArmaInfo * ainfo, const(DatasetRaw) * dset, int fullX)
{
    const(double) * y = dset.Z[ainfo.yno];
    double * dy = null;
    int * delta = null;
    int s = ainfo.pd;
    int k;
    int t;
    int t1 = 0;
    int err = 0;


    /* note: dy is a full length series (dset.n) */

    dy = cast(double *) malloc(dset.n * (*dy).sizeof);
    if (dy == null) {
	return GretlError.alloc;
    }

    delta = arima_delta_coeffs(ainfo.d, ainfo.D, s);
    if (delta == null) {
	free(dy);
	return GretlError.alloc;
    }

    for (t=0; t<dset.n; t++) {
	dy[t] = NADBL;
    }

    for (t=0; t<dset.n; t++) {
	if (na(y[t])) {
	    t1++;
	} else {
	    break;
	}
    }

    t1 += ainfo.d + ainfo.D * s;
    k = ainfo.d + s * ainfo.D;

    real_arima_difference_series(dy + t1, y, t1, ainfo.t2, delta, k);


    ainfo.y = dy;
    set_arima_ydiff(ainfo);

    if (arma_xdiff(ainfo)) {
	/* also difference the ARIMAX regressors */
	int xt1 = ainfo.t1;
	int xT = ainfo.T;

	if (fullX) {
	    xt1 = 0;
	    xT = ainfo.t2 + 1;
	}

	ainfo.dX = gretl_matrix_alloc(xT, ainfo.nexo);

	if (ainfo.dX == null) {
	    err = GretlError.alloc;
	} else {
	    double * val = ainfo.dX.val;
	    int i;
	    int vi;

	    for (i=0; i<ainfo.nexo; i++) {
		vi = ainfo.xlist[i+1];
		real_arima_difference_series(val, dset.Z[vi], xt1,
					     ainfo.t2, delta, k);
		val += xT;
	    }
	}
    }

    free(delta);

    return err;
}

extern (C) void arima_difference_undo(ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    free(ainfo.y);
    ainfo.y = cast(double *) dset.Z[ainfo.yno];

    if (ainfo.dX != null) {
	gretl_matrix_free(ainfo.dX);
	ainfo.dX = null;
    }

    unset_arima_ydiff(ainfo);
}


/* ---- arma.c proper ---- */

enum ArmaListPos {
    p_, q_, P_, Q_, d_, D_, y_
}


extern (C) int maybe_correct_MA(ArmaInfo * ainfo, double * theta, double * Theta)
{
    int err = 0;

    if (ainfo.q > 0) {
	err = flip_poly(theta, ainfo, 0, 0);
    }
    if (!err && ainfo.Q > 0) {
	err = flip_poly(Theta, ainfo, 0, 1);
    }
    return err;
}

/*
  Given an ARMA process $A(L)B(L) y_t = C(L)D(L) \epsilon_t$, finds the
  roots of the four polynomials -- or just two polynomials if seasonal
  AR and MA effects, B(L) and D(L) are not present -- and attaches
  this information to the ARMA model.

  pmod: MODEL pointer to which the roots info should be attached.

  ainfo: gives various pieces of information on the ARMA model,
  including seasonal and non-seasonal AR and MA orders.

  coeff: ifc + p + q + P + Q vector of coefficients (if an intercept
  is present it is element 0 and is ignored)

  returns: zero on success, non-zero on failure
*/

extern (C) int arma_model_add_roots(Model * pmod, ArmaInfo * ainfo, const(double) * coeff)
{
    const(double) * phi = coeff + ainfo.ifc;
    const(double) * Phi = phi + ainfo.np;
    const(double) * theta = Phi + ainfo.P;
    const(double) * Theta = theta + ainfo.nq;
    int nr = ainfo.p + ainfo.P + ainfo.q + ainfo.Q;
    int pmax;
    int qmax;
    int lmax;
    double * temp = null;
    double * tmp2 = null;
    CephesCmplx * rptr;
    CephesCmplx * roots = null;
    int i;
    int k;
    int cerr = 0;

    pmax = (ainfo.p > ainfo.P)? ainfo.p : ainfo.P;
    qmax = (ainfo.q > ainfo.Q)? ainfo.q : ainfo.Q;
    lmax = (pmax > qmax)? pmax : qmax;

    if (pmax == 0 && qmax == 0) {
	return 0;
    }

    temp = cast(double *) malloc((lmax + 1) * (*temp).sizeof);
    tmp2 = cast(double *) malloc((lmax + 1) * (*tmp2).sizeof);
    roots = cast(CephesCmplx *) malloc(nr * (*roots).sizeof);

    if (temp == null || tmp2 == null || roots == null) {
	free(temp);
	free(tmp2);
	free(roots);
	return GretlError.alloc;
    }

    temp[0] = 1.0;
    rptr = roots;

    if (ainfo.p > 0) {
	/* A(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo.p; i++) {
	    if (AR_included(ainfo, i)) {
		temp[i+1] = -phi[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}
	cerr = polrt(temp, tmp2, ainfo.p, rptr);
	rptr += ainfo.p;
    }

    if (!cerr && ainfo.P > 0) {
	/* B(L), seasonal */
	for (i=0; i<ainfo.P; i++) {
	    temp[i+1] = -Phi[i];
	}
	cerr = polrt(temp, tmp2, ainfo.P, rptr);
	rptr += ainfo.P;
    }

    if (!cerr && ainfo.q > 0) {
	/* C(L), non-seasonal */
	k = 0;
	for (i=0; i<ainfo.q; i++) {
	    if (MA_included(ainfo, i)) {
		temp[i+1] = theta[k++];
	    } else {
		temp[i+1] = 0;
	    }
	}
	cerr = polrt(temp, tmp2, ainfo.q, rptr);
	rptr += ainfo.q;
    }

    if (!cerr && ainfo.Q > 0) {
	/* D(L), seasonal */
	for (i=0; i<ainfo.Q; i++) {
	    temp[i+1] = Theta[i];
	}
	cerr = polrt(temp, tmp2, ainfo.Q, rptr);
    }

    free(temp);
    free(tmp2);

    if (cerr) {
	free(roots);
    } else {
	gretl_model_set_data(pmod, "roots", roots, GretlType.cmplx_array,
			     nr * (*roots).sizeof);
    }

    return 0;
}

/* add covariance matrix and standard errors based on Outer Product of
   Gradient
*/

extern (C) int arma_OPG_vcv(Model * pmod, void * data, int algo, double * b, double s2, int k, int T, Prn * prn)
{
    GretlMatrixRaw * G = null;
    GretlMatrixRaw * V = null;
    int err = 0;

    if (algo == 154) {
	G = numerical_score_matrix(b, T, k, &as154_llt_callback,
				   data, &err);
    } else if (algo == 197) {
	G = numerical_score_matrix(b, T, k, &as197_llt_callback,
				   data, &err);
    } else {
	G = numerical_score_matrix(b, T, k, &kalman_arma_llt_callback,
				   data, &err);
    }

    if (!err) {
	V = gretl_matrix_XTX_new(G);
	if (V == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err) {
	double rcond = gretl_symmetric_matrix_rcond(V, &err);

	if (!err && rcond < 1.0E-10) {
	    pprintf(prn, _("OPG: rcond = %g; will try Hessian\n").ptr, rcond);
	    err = 1;
	}
    }

    if (!err) {
	err = gretl_invert_symmetric_matrix(V);
    }

    if (!err) {
	gretl_matrix_multiply_by_scalar(V, s2);
	err = gretl_model_write_vcv(pmod, V);
    }

    gretl_matrix_free(G);
    gretl_matrix_free(V);

    return err;
}

extern (C) int arma_QML_vcv(Model * pmod, GretlMatrixRaw * H, void * data, int algo, double * b, double s2, int k, int T, Prn * prn)
{
    GretlMatrixRaw * G;
    int err = 0;

    if (algo == 154) {
	G = numerical_score_matrix(b, T, k, &as154_llt_callback,
				   data, &err);
    } else if (algo == 197) {
	G = numerical_score_matrix(b, T, k, &as197_llt_callback,
				   data, &err);
    } else {
	G = numerical_score_matrix(b, T, k, &kalman_arma_llt_callback,
				   data, &err);
    }

    if (!err) {
	gretl_matrix_divide_by_scalar(G, sqrt(s2));
	err = gretl_model_add_QML_vcv(pmod, GretlCmdIndex.arma, H, G,
				      null, Gretlopt.OPT_NONE, null);
    }

    gretl_matrix_free(G);

    return err;
}

extern (C) int arima_ydiff_only(ArmaInfo * ainfo)
{
    if ((ainfo.d > 0 || ainfo.D > 0) &&
	ainfo.nexo > 0 && !arma_xdiff(ainfo)) {
	return 1;
    } else {
	return 0;
    }
}

extern (C) int arma_use_opg(Gretlopt opt)
{
    int ret = 0; /* use of the Hessian is the default */

    if (opt & Gretlopt.OPT_G) {
	ret = 1;
    } else if (libset_get_int(SetKey.arma_vcv) == MLVCVType.op) {
	ret = 1;
    }

    return ret;
}

private void kalman_rescale_y(GretlMatrixRaw * y, ArmaInfo * ainfo)
{
    int i;


    for (i=0; i<y.rows; i++) {
	if (!isnan(y.val[i])) {
	    y.val[i] -= ainfo.yshift;
	    y.val[i] *= ainfo.yscale;
	}
    }
}

/* for Kalman: convert from full-length y series to
   y vector of length ainfo.T */

extern (C) GretlMatrixRaw * form_arma_y_vector(ArmaInfo * ainfo, int * err)
{
    GretlMatrixRaw * yvec;

    yvec = gretl_vector_from_series(ainfo.y, ainfo.t1, ainfo.t2);

    if (yvec == null) {
	*err = GretlError.alloc;
    } else {
	if (ainfo.yscale != 1.0) {
	    kalman_rescale_y(yvec, ainfo);
	}
    }

    return yvec;
}

extern (C) GretlMatrixRaw * form_arma_X_matrix(ArmaInfo * ainfo, const(DatasetRaw) * dset, int * err)
{
    GretlMatrixRaw * X;
    int missop;


    if (arma_na_ok(ainfo)) {
	missop = MMissingCode.ok;
    } else {
	missop = MMissingCode.error;
    }

    X = gretl_matrix_data_subset(ainfo.xlist, dset,
				 ainfo.t1, ainfo.t2,
				 missop, err);


    return X;
}

bool y_missing(double y)
{
    return na(y) || isnan(y);
}

/* support for AS 154 and AS 197 */



private void arma_init_message(ArmaInfo * ainfo)
{
    pprintf(ainfo.prn, "\n%s: ", _("ARMA initialization").ptr);

    if (ainfo.init_ == IniMethod.user) {
	pprintf(ainfo.prn, "%s\n\n", _("user-specified values").ptr);
    } else if (ainfo.init_ == IniMethod.hr) {
	pprintf(ainfo.prn, "%s\n\n", _("Hannan-Rissanen method").ptr);
    } else if (ainfo.init_ == IniMethod.small) {
	pprintf(ainfo.prn, "%s\n\n", _("small MA values").ptr);
    } else if (ainfo.init_ == IniMethod.nls) {
	pprintf(ainfo.prn, "%s\n\n", _("using nonlinear AR model").ptr);
    } else if (ainfo.init_ == IniMethod.ols) {
	pprintf(ainfo.prn, "%s\n\n", _("using linear AR model").ptr);
    }
}

private int user_arma_init(double * coeff, ArmaInfo * ainfo)
{
    int i, nc = n_initvals();

    if (nc == 0) {
	return 0;
    } else if (nc < ainfo.nc) {
	pprintf(ainfo.prn, _("ARMA initialization: need %d coeffs but got %d\n").ptr,
		ainfo.nc, nc);
	return GretlError.data;
    }

    if (arma_exact_ml(ainfo)) {
	/* user-specified initializer is handled within BFGSmax */
	for (i=0; i<ainfo.nc; i++) {
	    coeff[i] = 0.0;
	}
    } else {
	GretlMatrixRaw * m = get_initvals();

	for (i=0; i<ainfo.nc; i++) {
	    coeff[i] = gretl_vector_get(m, i);
	}
	gretl_matrix_free(m);
    }

    ainfo.init_ = IniMethod.user;

    return 0;
}

/* Should we try Hannan-Rissanen initialization of ARMA
   coefficients? */

private int prefer_hr_init(ArmaInfo * ainfo)
{
    int ret = 0;

    if (ainfo.q > 1 || ainfo.Q > 0) {
	ret = 1;
	if (arma_xdiff(ainfo)) {
	    /* don't use for ARIMAX (yet?) */
	    ret = 0;
	} else if (ainfo.T < 100) {
	    /* unlikely to work well with small sample */
	    ret = 0;
	} else if (ainfo.p > 0 && ainfo.P > 0) {
	    /* not sure about this: HR catches the MA terms, but NLS
	       handles the seasonal/non-seasonal AR interactions
	       better?
	    */
	    ret = 0;
	} else if ((ainfo.P > 0 && ainfo.p >= ainfo.pd) ||
		   (ainfo.Q > 0 && ainfo.q >= ainfo.pd)) {
	    /* overlapping seasonal/non-seasonal orders screw things up */
	    ret = 0;
	} else if (ret && arma_exact_ml(ainfo)) {
	    /* screen for cases where we'll use NLS */
	    if (ainfo.P > 0) {
		ret = 0;
	    } else if (ainfo.p + ainfo.P > 0 && ainfo.nexo > 0) {
		ret = 0;
	    } else if (ainfo.Q > 0 && arma_missvals(ainfo)) {
		ret = 0;
	    }
	}
    }


    return ret;
}

/* estimate an ARIMA (0,d,0) x (0,D,0) model via OLS */

private int arima_by_ls(const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod)
{
    GretlMatrixRaw * X;
    GretlMatrixRaw * b;
    GretlMatrixRaw * u;
    GretlMatrixRaw * V;
    double x;
    double s2;
    int i;
    int t;
    int k = ainfo.dX.cols;
    int err = 0;

    if (ainfo.ifc) {
	/* the constant will not have been included in ainfo.dX */
	X = gretl_matrix_alloc(ainfo.T, k + 1);
	if (X == null) {
	    return GretlError.alloc;
	}
	for (i=0; i<=k; i++) {
	    for (t=0; t<ainfo.T; t++) {
		if (i == 0) {
		    gretl_matrix_set(X, t, i, 1.0);
		} else {
		    x = gretl_matrix_get(ainfo.dX, t, i-1);
		    gretl_matrix_set(X, t, i, x);
		}
	    }
	}
	k++;
    } else {
	X = ainfo.dX;
    }

    b = gretl_column_vector_alloc(k);
    u = gretl_column_vector_alloc(ainfo.T);
    V = gretl_matrix_alloc(k, k);

    if (b == null || u == null || V == null) {
	err = GretlError.alloc;
    } else {
	GretlMatrixRaw y;

	gretl_matrix_init(&y);
	y.rows = ainfo.T;
	y.cols = 1;
	y.val = ainfo.y + ainfo.t1;
	gretl_matrix_set_t1(&y, ainfo.t1);
	gretl_matrix_set_t2(&y, ainfo.t2);

	err = gretl_matrix_ols(&y, X, b, V, u, &s2);
    }

    if (!err) {
	pmod.ncoeff = k;
	pmod.full_n = dset.n;
	err = gretl_model_allocate_storage(pmod);
    }

    if (!err) {
	for (i=0; i<k; i++) {
	    pmod.coeff[i] = b.val[i];
	}
	for (t=0; t<ainfo.T; t++) {
	    pmod.uhat[t + ainfo.t1] = u.val[t];
	}
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
    	pmod.ybar = gretl_mean(ainfo.t1, ainfo.t2, ainfo.y);
	pmod.sdy = gretl_stddev(ainfo.t1, ainfo.t2, ainfo.y);
	pmod.nobs = ainfo.T;
    }

    gretl_matrix_free(b);
    gretl_matrix_free(u);
    gretl_matrix_free(V);

    if (X != ainfo.dX) {
	gretl_matrix_free(X);
    }

    return err;
}

bool simple_OLS_case(const(ArmaInfo) * a)
{
    return a.p == 0 && a.P == 0 && a.q == 0 && a.Q == 0 && a.d == 0 && a.D == 0;
}

/* calculate info criteria for compatibility with ML? */

private int arma_via_OLS(ArmaInfo * ainfo, const(double) * coeff, const(DatasetRaw) * dset, Model * pmod)
{
    int err = 0;

    ainfo.flags |= ArmaFlags.ls;

    if (arma_xdiff(ainfo)) {
	err = arima_by_ls(dset, ainfo, pmod);
    } else if (simple_OLS_case(ainfo)) {
        err = arma_by_simple_ols(coeff, dset, ainfo, pmod);
    } else {
	err = arma_by_ls(coeff, dset, ainfo, pmod);
    }

    if (!err) {
	ArmaFlags f = arma_exact_ml(ainfo) ? ArmaFlags.ols : ArmaFlags.ls;

	pmod.t1 = ainfo.t1;
	pmod.t2 = ainfo.t2;
	pmod.full_n = dset.n;
	write_arma_model_stats(pmod, ainfo, dset);
	if (arma_exact_ml(ainfo)) {
	    /* In the case of ainfo.nc == 0 (no coefficients
	       actually estimated), pmod.ncoeff will be 1, since
	       we add a dummy constant with value 0.0 That "1"
	       will account for the variance estimate so that
	       addk should be zero. Otherwise we add 1 for the
	       variance estimate.
	    */
	    int addk = ainfo.nc == 0 ? 0 : 1;

	    mle_criteria(pmod, addk);
	} else {
	    arma_model_add_roots(pmod, ainfo, pmod.coeff);
	}
	gretl_model_set_int(pmod, "arma_flags", f);
    }

    if (!err && pmod.errcode) {
	err = pmod.errcode;
    }

    return err;
}

/* Set flag to indicate differencing of exogenous regressors, in the
   case of an ARIMAX model using native exact ML -- unless this is
   forbidden by OPT_Y (--y-diff-only).  Note that we don't do this
   when we're using conditional ML (BHHH).
*/

private void maybe_set_xdiff_flag(ArmaInfo * ainfo, Gretlopt opt)
{
    if (arma_exact_ml(ainfo) &&
	(ainfo.d > 0 || ainfo.D > 0) &&
	ainfo.nexo > 0 && !(opt & Gretlopt.OPT_Y)) {
	ainfo.pflags |= PrivFlags.xdiff;
    }
}

/* Respond to OPT_S (--stdx): standardize exogenous regressors */

private int arma_standardize_x(ArmaInfo * ainfo, DatasetRaw * dset)
{
    int orig_v = dset.v;
    int err = 0;

    ainfo.xstats = gretl_matrix_alloc(ainfo.nexo, 2);
    if (ainfo.xstats == null) {
	return GretlError.alloc;
    }

    err = dataset_add_series(dset, ainfo.nexo);

    if (!err) {
	double xbar;
	double sdx;
	int i;
	int vi;
	int vj;
	int t;

	for (i=0; i<ainfo.nexo && !err; i++) {
	    vi = ainfo.xlist[i+1];
	    err = gretl_moments(ainfo.t1, ainfo.t2, dset.Z[vi],
				null, &xbar, &sdx, null, null, 1);
	    if (!err) {
		vj = orig_v + i;
		for (t=0; t<dset.n; t++) {
		    dset.Z[vj][t] = (dset.Z[vi][t] - xbar) / sdx;
		}
		/* replace x-ref with standardized version */
		ainfo.xlist[i+1] = vj;
		/* and record the stats used */
		gretl_matrix_set(ainfo.xstats, i, 0, xbar);
		gretl_matrix_set(ainfo.xstats, i, 1, sdx);
	    }
	}
    }

    if (!err) {
	set_arma_stdx(ainfo);
    }

    return err;
}

/* Set flag to allow NAs within the sample range for an
   ARMA model using native exact ML, or for one that can
   be estimated via OLS (with no MA term).
*/

private void maybe_allow_missvals(ArmaInfo * ainfo)
{
    if (arma_exact_ml(ainfo)) {
	ainfo.pflags |= PrivFlags.naok;
    } else if (ainfo.q == 0 && ainfo.Q == 0) {
	ainfo.pflags |= PrivFlags.naok;
    }
}

private int check_arma_options(Gretlopt opt)
{
    int err;

    /* can't specify LBFGS or --robust with conditional ML */
    err = options_incompatible_with(opt, Gretlopt.OPT_C, Gretlopt.OPT_L | Gretlopt.OPT_R);

    if (!err) {
	/* nor more than one of AS 154, CML, Kalman */
	err = incompatible_options(opt, Gretlopt.OPT_A | Gretlopt.OPT_C | Gretlopt.OPT_K);
    }

    return err;
}

private int arma_precheck(const(int) * list, const(int) * pqspec, Gretlopt opt, DatasetRaw * dset, ArmaInfo * painfo, int * dim)
{
    ArmaInfo ainfo;
    int err;

    memset(&ainfo, 0, ArmaInfo.sizeof);

    arma_info_init(&ainfo, opt, pqspec, dset);
    err = check_arma_options(opt);
    if (!err) {
	ainfo.alist = gretl_list_copy(list);
	if (ainfo.alist == null) {
	    err = GretlError.alloc;
	}
    }
    if (!err) {
	err = arma_check_list(&ainfo, dset, opt);
    }
    if (!err) {
	if (painfo != null) {
	    *painfo = ainfo;
	} else if (dim != null) {
	    /* just transcribe the basics */
	    dim[ArmaListPos.p_] = ainfo.p;
	    dim[ArmaListPos.q_] = ainfo.q;
	    dim[ArmaListPos.P_] = ainfo.P;
	    dim[ArmaListPos.Q_] = ainfo.Q;
	    dim[ArmaListPos.d_] = ainfo.d;
	    dim[ArmaListPos.D_] = ainfo.D;
	    dim[ArmaListPos.y_] = ainfo.yno;
	}
	if (painfo == null) {
	    free(ainfo.alist);
	}
    }

    return err;
}

extern (C) Model arma_model(const(int) * list, const(int) * pqspec, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    double * coeff = null;
    Model armod;
    ArmaInfo ainfo_s;
    ArmaInfo * ainfo;
    int missv = 0;
    int misst = 0;
    int orig_v = dset.v;
    int err;

    gretl_model_init(&armod, dset);

    err = arma_precheck(list, pqspec, opt, dset, &ainfo_s, null);
    if (err) {
	armod.errcode = err;
	return armod;
    }

    ainfo = &ainfo_s;
    if (opt & Gretlopt.OPT_V) {
	ainfo.prn = prn;
    }

    if (!err) {
	/* calculate maximum lag */
	maybe_set_xdiff_flag(ainfo, opt);
	calc_max_lag(ainfo);
    }

    if (!err) {
	/* adjust sample range if need be */
	maybe_allow_missvals(ainfo);
	err = arma_adjust_sample(ainfo, dset, &missv, &misst);
	if (err) {
	    if (missv > 0 && misst > 0) {
		gretl_errmsg_sprintf(_("Missing value encountered for " ~ "variable %d, obs %d").ptr, missv, misst);
	    }
	} else if (missv > 0) {
	    set_arma_missvals(ainfo);
	}
    }

    if (!err && ainfo.nexo > 0 && arma_exact_ml(ainfo)) {
	/* FIXME check conditionality more rigorously */
	if ((opt & Gretlopt.OPT_S) && !arma_xdiff(ainfo)) {
	    /* --stdx */
	    err = arma_standardize_x(ainfo, dset);
	}
    }

    if (!err) {
	/* allocate initial coefficient vector */
	coeff = cast(double *) malloc(ainfo.nc * (*coeff).sizeof);
	if (coeff == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err) {
	/* organize the dependent variable */
	ainfo.y = cast(double *) dset.Z[ainfo.yno];
	if (ainfo.d > 0 || ainfo.D > 0) {
	    if (arma_missvals(ainfo)) {
		/* for now: insist on native Kalman, since only it
		   handles the levels formulation of ARIMA
		*/
		opt &= ~Gretlopt.OPT_A;
		opt |= Gretlopt.OPT_K;
		set_arima_levels(ainfo);
	    } else {
		/* note: this replaces ainfo.y */
		err = arima_difference(ainfo, dset, 0);
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    if (ainfo.p == 0 && ainfo.P == 0 &&
	ainfo.q == 0 && ainfo.Q == 0 &&
	arma_exact_ml(ainfo) && !arma_missvals(ainfo)) {
	/* pure "I" model, no NAs: OLS provides the MLE */
	err = arma_via_OLS(ainfo, null, dset, &armod);
	goto bailout; /* estimation handled */
    }

    /* start initialization of the coefficients */

    /* first see if the user specified some values */
    err = user_arma_init(coeff, ainfo);
    if (err) {
	goto bailout;
    }

    if (!arma_exact_ml(ainfo) && ainfo.q == 0 && ainfo.Q == 0) {
	/* for a pure AR model, the conditional MLE is least
	   squares (OLS or NLS); in the NLS case a user-specified
	   initializer may be useful, if present
	*/
	const(double) * b = ainfo.init_ ? coeff : null;

	err = arma_via_OLS(ainfo, b, dset, &armod);
	goto bailout; /* estimation handled */
    }

    /* see if it may be helpful to scale the dependent
       variable, if we're doing exact ML */
    if (arma_exact_ml(ainfo) && ainfo.ifc && ainfo.init_ != IniMethod.user) {
	maybe_set_yscale(ainfo);
    }

    /* try Hannan-Rissanen init, if suitable */
    if (!ainfo.init_ && prefer_hr_init(ainfo)) {
	hr_arma_init(coeff, dset, ainfo);
    }

    /* initialize via AR model by OLS or NLS, adding minimal
       MA coefficients if needed: this is the fallback if
       Hannan-Rissanen fails, but also the default if the
       conditions of applicability of H-R are not met
    */
    if (!err && !ainfo.init_) {
	err = ar_arma_init(coeff, dset, ainfo, &armod, opt);
    }

    if (ainfo.prn != null && ainfo.init_) {
	arma_init_message(ainfo);
    }

    if (!err) {
	clear_model_xpx(&armod);
	if (arma_exact_ml(ainfo)) {
	    if (opt & Gretlopt.OPT_K) {
		err = kalman_arma(coeff, dset, ainfo, &armod, opt);
	    } else {
		err = as_arma(coeff, dset, ainfo, &armod, opt);
	    }
	} else {
	    err = bhhh_arma(coeff, dset, ainfo, &armod, opt);
	}
    }

 bailout:

    if (err && !armod.errcode) {
	armod.errcode = err;
    }

    if (!err) {
	transcribe_extra_info(ainfo, &armod);
    }

    if (!armod.errcode) {
	gretl_model_smpl_init(&armod, dset);
    }

    free(coeff);
    arma_info_cleanup(ainfo);

    if (dset.v > orig_v) {
	dataset_drop_last_variables(dset, dset.v - orig_v);
    }

    return armod;
}

struct ArmaSel {
    GretlMatrixRaw * m;
    int * list;
    int[7] dim;
    double[3] mincrit;
    int[3] minidx;
    int seasonal;
    int qpos; /* position of 'q' in command list */
    int Ppos; /* position of 'P' in command list */
    int Qpos; /* position of 'Q' in command list */
    int width; /* width to print criteria */
    int ndec; /* number of decimal places */
    int nbad; /* count of non-estimated specifications */
};



/* Allocate the matrix required to hold arma lag selection
   information, and add column names. The dimensions of this
   matrix depend on whether seasonal terms should be
   included.
*/

private int add_select_matrix(ArmaSel * asel)
{
    int i;
    int cols;
    int rows = 0;
    int err = 0;

    for (i=0; i<4; i++) {
	if (asel.dim[i] > 0) {
	    if (rows > 0) {
		rows *= asel.dim[i] + 1;
	    } else {
		rows = asel.dim[i] + 1;
	    }
	}
    }

    if (rows == 0) {
	gretl_errmsg_set("arima --select: search space is null");
	return GretlError.data;
    }

    cols = asel.seasonal ? 8 : 6;
    asel.m = gretl_matrix_alloc(rows, cols);

    if (asel.m == null) {
	err = GretlError.alloc;
    } else {
	const(char) *[] S1 = ["p", "q", "AIC", "BIC", "HQC", "loglik"];
	const(char) *[] S2 = ["p", "q", "P", "Q", "AIC", "BIC", "HQC", "loglik"];
        char **S = strings_array_new(cols);

        if (S != null) {
	    for (i=0; i<cols; i++) {
		S[i] = gretl_strdup(cols == 6 ? S1[i] : S2[i]);
	    }
            gretl_matrix_set_colnames(asel.m, S);
        }
    }

    return err;
}

/* Output a header for the arma lag-selection matrix, including
   information gleaned from the last model estimated, @amod.
*/

private void arma_select_header(ArmaSel * asel, Model * amod, DatasetRaw * dset, Prn * prn)
{
    char[8] word;
    char[32] s1;
    char[32] s2;
    int w = asel.width;
    int * num = asel.dim.ptr;

    if (num[ArmaListPos.d_] || num[ArmaListPos.D_]) {
	strcpy(word.ptr, "ARIMA".ptr);
    } else {
	strcpy(word.ptr, "ARMA".ptr);
    }
    s2[0] = '\0';

    if (num[ArmaListPos.d_]) {
	sprintf(s1.ptr, "(p, %d, q)", num[ArmaListPos.d_]);
    } else {
	strcpy(s1.ptr, "(p, q)".ptr);
    }
    if (num[ArmaListPos.D_]) {
	if (asel.seasonal) {
	    sprintf(s2.ptr, "(P, %d, Q)", num[ArmaListPos.D_]);
	} else {
	    sprintf(s2.ptr, "(%d, %d, %d)", num[ArmaListPos.P_], num[ArmaListPos.D_], num[ArmaListPos.Q_]);
	}
    } else if (num[ArmaListPos.P_] || num[ArmaListPos.Q_]) {
	strcpy(s2.ptr, "(P, Q)".ptr);
    }

    pputc(prn, '\n');
    arma_extra_info(amod, prn);
    pprintf(prn, "Dependent variable %s, T = %d\n", dset.varname[num[ArmaListPos.y_]],
	    amod.nobs);
    pprintf(prn, "Criteria for %s%s%s specifications\n", word.ptr, s1.ptr, s2.ptr);
    if (asel.seasonal) {
	dashline(66, prn);
	pprintf(prn, " p, q, P, Q %*s %*s %*s %*s\n", w, "AIC".ptr, w+1,
		"BIC".ptr, w+1, "HQC".ptr, w+1, "loglik".ptr);
	dashline(66, prn);
    } else {
	dashline(60, prn);
	pprintf(prn, " p, q %*s %*s %*s %*s\n", w, "AIC".ptr, w+1, "BIC".ptr,
		w+1, "HQC".ptr, w+1, "loglik".ptr);
	dashline(60, prn);
    }
}

/* Notes to be printed under the lag-selection matrix */

private void arma_select_footer(ArmaSel * asel, Prn * prn)
{
    pputc(prn, '\n');
    pputs(prn, "'*' indicates best, per criterion\n");
    if (asel.nbad > 0) {
	pputs(prn, "'NA' indicates that a specification could not be estimated\n");
    }
    pputs(prn, "Log-likelihood ('loglik') is provided for reference\n");
    pputc(prn, '\n');
}

/* Note: asel.dim should already be filled out when arma_sel_init()
   is called
*/

private void arma_sel_init(ArmaSel * asel, const(int) * list)
{
    int i;
    int s1;

    for (i=0; i<3; i++) {
	asel.mincrit[i] = NADBL;
	asel.minidx[i] = 0;
    }

    asel.Ppos = asel.Qpos = 0;
    asel.seasonal = 0;
    asel.m = null;
    asel.list = null;
    asel.width = 12;
    asel.ndec = 4;
    asel.nbad = 0;

    s1 = gretl_list_separator_position(list);
    asel.qpos = (s1 == 4)? 3 : 2;

    if (asel.dim[ArmaListPos.P_] || asel.dim[ArmaListPos.Q_]) {
	asel.seasonal = 1;
	asel.Ppos = s1 + 1;
	for (i=s1+1; i<=list[0]; i++) {
	    if (list[i] == LISTSEP) {
		asel.Qpos = i-1;
		break;
	    }
	}
    }
}

private void record_criterion(ArmaSel * asel, Model * amod, int i, int j, int k)
{
    if (j < 3) {
        /* one of the three info criteria */
	double cij = amod.criterion[j];

	gretl_matrix_set(asel.m, i, k, cij);
	if (na(asel.mincrit[j])) {
	    asel.mincrit[j] = cij;
	    asel.minidx[j] = i;
	} else if (cij < asel.mincrit[j]) {
	    asel.mincrit[j] = cij;
	    asel.minidx[j] = i;
	}
    } else {
        /* last column: log-likelihood */
	gretl_matrix_set(asel.m, i, k, amod.lnL);
    }
}

private void fill_arma_sel_matrix(ArmaSel * asel, Model * amod, Gretlopt opt, DatasetRaw * dset)
{
    int p;
    int q;
    int P;
    int Q;
    int j;
    int i = 0;

    for (p=0; p<=asel.dim[ArmaListPos.p_]; p++) {
	asel.list[1] = p;
        for (q=0; q<=asel.dim[ArmaListPos.q_]; q++) {
	    asel.list[asel.qpos] = q;
	    if (asel.seasonal) {
		/* include seasonal terms */
		for (P=0; P<=asel.dim[ArmaListPos.P_]; P++) {
		    asel.list[asel.Ppos] = P;
		    for (Q=0; Q<=asel.dim[ArmaListPos.Q_]; Q++) {
			asel.list[asel.Qpos] = Q;
			gretl_matrix_set(asel.m, i, 0, p);
			gretl_matrix_set(asel.m, i, 1, q);
			gretl_matrix_set(asel.m, i, 2, P);
			gretl_matrix_set(asel.m, i, 3, Q);
			*amod = arma_model(asel.list, null, dset, opt, null);
			for (j=0; j<4; j++) {
			    if (amod.errcode) {
				gretl_matrix_set(asel.m, i, j+4, NADBL);
				asel.nbad += 1;
			    } else {
				record_criterion(asel, amod, i, j, j+4);
			    }
			}
			if (i < asel.m.rows - 1) {
			    clear_model(amod);
			}
			i++;
		    }
		}
	    } else {
		/* exclude seasonal terms */
		gretl_matrix_set(asel.m, i, 0, p);
		gretl_matrix_set(asel.m, i, 1, q);
		*amod = arma_model(asel.list, null, dset, opt, null);
		for (j=0; j<4; j++) {
		    if (amod.errcode) {
			gretl_matrix_set(asel.m, i, j+2, NADBL);
			asel.nbad += 1;
		    } else {
			record_criterion(asel, amod, i, j, j+2);
		    }
		}
		if (i < asel.m.rows - 1) {
		    clear_model(amod);
		}
		i++;
	    }
        }
    }
}

/* pretty-printer for arma lag selection matrix */

private void print_arma_sel_matrix(ArmaSel * asel, Prn * prn)
{
    int w;
    int i;
    int j;
    int p;
    int q;
    int P;
    int Q;
    double cij;

    i = 0;
    for (p=0; p<=asel.dim[ArmaListPos.p_]; p++) {
	for (q=0; q<=asel.dim[ArmaListPos.q_]; q++) {
	    if (asel.seasonal) {
		for (P=0; P<=asel.dim[ArmaListPos.P_]; P++) {
		    for (Q=0; Q<=asel.dim[ArmaListPos.Q_]; Q++) {
			pprintf(prn, "%2d,%2d,%2d,%2d", p, q, P, Q);
			for (j=0; j<4; j++) {
			    cij = gretl_matrix_get(asel.m, i, j+4);
			    if (na(cij)) {
				w = j==0 ? asel.width : asel.width + 1;
				pprintf(prn, " %*s", w, "NA".ptr);
			    } else {
				pprintf(prn, " %*.*f", asel.width, asel.ndec, cij);
				if (j < 3) {
				    pputc(prn, i == asel.minidx[j] ? '*' : ' ');
				}
			    }
			}
			pputc(prn, '\n');
			i++;
		    }
		}
	    } else {
		pprintf(prn, "%2d,%2d", p, q);
		for (j=0; j<4; j++) {
		    cij = gretl_matrix_get(asel.m, i, j+2);
		    if (na(cij)) {
			w = j==0 ? asel.width : asel.width + 1;
			pprintf(prn, " %*s", w, "NA".ptr);
		    } else {
			pprintf(prn, " %*.*f", asel.width, asel.ndec, cij);
			if (j < 3) {
			    pputc(prn, i == asel.minidx[j] ? '*' : ' ');
			}
		    }
		}
		pputc(prn, '\n');
		i++;
	    }
	}
    }
}

/* Specific to conditional ML estimation: ensure a common sample by
   first estimating the biggest model and restricting the sample range
   to the observations usable in that case.
*/

private int set_common_sample(ArmaSel * asel, Model * amod, Gretlopt opt, DatasetRaw * dset)
{
    int err = 0;

    /* set all orders to max */
    asel.list[1] = asel.dim[ArmaListPos.p_];
    asel.list[asel.qpos] = asel.dim[ArmaListPos.q_];
    if (asel.seasonal) {
	asel.list[asel.Ppos] = asel.dim[ArmaListPos.P_];
	asel.list[asel.Qpos] = asel.dim[ArmaListPos.Q_];
    }

    *amod = arma_model(asel.list, null, dset, opt, null);
    if (amod.errcode) {
	gretl_errmsg_set("Error estimating biggest spec via CML");
	err = amod.errcode;
    } else {
	dset.t1 = amod.t1;
	dset.t2 = amod.t2;
    }

    clear_model(amod);

    return err;
}

/* Provides built-in means of selecting the AR and MA orders of an
   AR(I)MA model via Information Criteria. We accept a general
   specification of the form

   (p d q)(P D Q)

   and search over p, q, P and Q, taking their values on input as
   maxima; d and D remain clamped at their input values.
*/

private int real_arma_select(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn, GretlMatrixRaw ** pm)
{
    Model amod;
    ArmaSel asel;

    memset(&asel, 0, ArmaSel.sizeof);
    Gretlopt aopt = opt | Gretlopt.OPT_Q;
    int save_t1 = dset.t1;
    int save_t2 = dset.t2;
    int err = 0;

    err = arma_precheck(list, null, opt, dset, null, asel.dim.ptr);
    if (err) {
	return err;
    }

    arma_sel_init(&asel, list);
    err = add_select_matrix(&asel);
    if (err) {
        return err;
    }

    /* to preserve const property of @list */
    asel.list = gretl_list_copy(list);


    if (opt & Gretlopt.OPT_C) {
	/* conditional ML: ensure a common sample */
	err = set_common_sample(&asel, &amod, aopt, dset);
	if (err) {
	    gretl_matrix_free(asel.m);
	    free(asel.list);
	    return err;
	}
    }

    /* fill the results matrix */
    fill_arma_sel_matrix(&asel, &amod, aopt, dset);

    if (!(opt & Gretlopt.OPT_Q)) {
        /* print as table */
	arma_select_header(&asel, &amod, dset, prn);
	print_arma_sel_matrix(&asel, prn);
	arma_select_footer(&asel, prn);
    }

    free(asel.list);

    if (asel.m != null) {
	if (pm != null) {
	    *pm = asel.m;
	} else {
	    record_matrix_test_result(asel.m, null);
	}
    }

    clear_model(&amod);

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    return err;
}

extern (C) int arma_select(const(int) * list, const(int) * pqspec, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    if (pqspec != null) {
        return GretlError.data;
    } else {
	return real_arma_select(list, dset, opt, prn, null);
    }
}

extern (C) int gui_arma_select(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn, GretlMatrixRaw ** pm)
{
    return real_arma_select(list, dset, opt, prn, pm);
}
