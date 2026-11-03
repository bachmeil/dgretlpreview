/* Port of plugin/garch.c: the GARCH plugin driver (Fiorentini, Calzolari and Panattoni mixed-gradient algorithm): list handling, the
   OLS pre-estimation and initialization (also from an ARMA fit of the squared residuals), scaling of the dependent variable, the call of the
   estimator (garch_estimate of fcp.c, or garch_estimate_mod of mod_garch.c with the BFGS-based variants), and the model statistics and
   parameter names. garch_model is the entry point of the plugin.

   garch.h only declares garch_estimate, garch_analytical_hessian (both in gretl.plugin.fcp) and garch_estimate_mod (gretl.plugin.mod_garch). The debug output
   (VPARM_DEBUG, off) and the code under #if 0 are not in the port; GARCH_AUTOCORR_TEST and GARCH_SCALE_SIGMA are on, as in the C. */

module gretl.plugin.garch;

import gretl.compare : autocorr_test;
import gretl.dataset : dataset_add_series, dataset_drop_last_variables;
import gretl.estimate : arma, lsq;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_list_copy, gretl_list_new, gretl_list_sublist;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_matrix_free, gretl_zero_matrix_new;
import gretl.gretl_model : clear_model, clear_model_xpx, gretl_model_allocate_param_names, gretl_model_destroy_tests, gretl_model_get_data, gretl_model_init, gretl_model_set_data, gretl_model_set_double, gretl_model_set_int, gretl_model_set_param_name, gretl_model_set_vcv_info, gretl_model_write_coeffs, gretl_model_write_vcv, mle_criteria, model_count_minus;
import gretl.gretl_prn : pprintf, pputc, pputs;
import gretl.gretl_utils : doubles_array_free, doubles_array_new, get_last_pvalue, get_last_test_statistic;
import gretl.libset : get_initvals, libset_get_int, n_initvals;
import gretl.missingh : NADBL, na;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import gretl.gretl_prn : Prn;
import gretl.plugin.fcp : garch_estimate;
import gretl.plugin.mod_garch : garch_estimate_mod;
import core.stdc.string : strcpy;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.libseth : SetKey;
import gretl.gretl_modelh : MLVCVType, VCVMajorType;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;




/* GARCH plugin for gretl using the Fiorentini, Calzolari and
   Panattoni mixed-gradient algorithm.
*/




enum PQ_MAX = 7;               /* max sum of GARCH p and q */
enum GARCH_PARAM_MAX = 0.999;

private void add_garch_varnames(Model * pmod, const(DatasetRaw) * dset, const(int) * list)
{
    char[24] tmp;
    int p = list[1]; /* GARCH beta terms */
    int q = list[2]; /* ARCH alpha terms > 0 */
    int r = list[0] - 4; /* regressors */
    int np = 1 + p + q + r; /* the "1" is for alpha(0) */
    int i;
    int j;

    free(pmod.list);
    pmod.list = gretl_list_copy(list);

    gretl_model_allocate_param_names(pmod, np);
    if (pmod.errcode) {
	return;
    }

    j = 0;

    for (i=0; i<r; i++) {
	gretl_model_set_param_name(pmod, j++, dset.varname[pmod.list[5+i]]);
    }

    gretl_model_set_param_name(pmod, j++, "alpha(0)");

    for (i=0; i<q; i++) {
	sprintf(tmp.ptr, "alpha(%d)".ptr, i + 1);
	gretl_model_set_param_name(pmod, j++, tmp.ptr);
    }

    for (i=0; i<p; i++) {
	sprintf(tmp.ptr, "beta(%d)".ptr, i + 1);
	gretl_model_set_param_name(pmod, j++, tmp.ptr);
    }
}

private void rescale_results(double * theta, GretlMatrixRaw * V, double scale, int npar, int nc)
{
    double vij;
    double sfi;
    double sf;
    double sc2 = scale * scale;
    int i;
    int j;

    for (i=0; i<nc; i++) {
	theta[i] *= scale;
    }

    theta[nc] *= sc2;

    for (i=0; i<npar; i++) {
	sfi = (i < nc)? scale : (i == nc)? sc2 : 1.0;
	for (j=0; j<=i; j++) {
	    sf = (j < nc)? scale * sfi : (j == nc)? sc2 * sfi : sfi;
	    vij = gretl_matrix_get(V, i, j) * sf;
	    gretl_matrix_set(V, i, j, vij);
	    gretl_matrix_set(V, j, i, vij);
	}
    }
}

private int write_garch_stats(Model * pmod, const(int) * list, const(DatasetRaw) * dset, double * theta, GretlMatrixRaw * V, int p, int q, double scale, const(double) * e, const(double) * h, int npar, int nc, int pad, int ifc, Prn * prn)
{
    double * garch_h;
    double den;
    const(double) * vcoef;
    int ynum = list[4];
    int nvp = list[1] + list[2];
    int xvars = list[0] - 4;
    int i;
    int err;

    err = gretl_model_set_int(pmod, "garch_p", p);
    if (!err) {
	err = gretl_model_set_int(pmod, "garch_q", q);
    }

    if (!err) {
	if (scale != 1.0) {
	    rescale_results(theta, V, scale, npar, nc);
	}
	err = gretl_model_write_coeffs(pmod, theta, npar);
    }

    if (err) {
	return err;
    }

    gretl_model_write_vcv(pmod, V);

    /* verbose? */
    if (prn != null) {
	for (i=0; i<npar; i++) {
	    pprintf(prn, "theta[%d]: %#14.6g (%#0.6g)\n", i, theta[i],
		    pmod.sderr[i]);
	}
	pputc(prn, '\n');
    }

    pmod.ess = 0.0;
    for (i=pmod.t1; i<=pmod.t2; i++) {
	pmod.uhat[i] = e[i + pad] * scale;
	pmod.ess += pmod.uhat[i] * pmod.uhat[i];
	pmod.yhat[i] = dset.Z[ynum][i] * scale - pmod.uhat[i];
    }

    vcoef = pmod.coeff + xvars;

    /* set sigma to its unconditional or steady-state value */
    den = 1.0;
    for (i=1; i<=nvp; i++) {
	den -= vcoef[i];
    }
    pmod.sigma = sqrt(vcoef[0] / den);

    pmod.adjrsq = NADBL;
    pmod.fstt = NADBL;

    mle_criteria(pmod, 1);

    pmod.ci = GretlCmdIndex.garch;
    pmod.ifc = ifc;

    add_garch_varnames(pmod, dset, list);

    /* add predicted error variance to model */
    garch_h = cast(double *) malloc(dset.n * (*garch_h).sizeof);
    if (garch_h != null) {
	for (i=0; i<dset.n; i++) {
	    if (i < pmod.t1 || i > pmod.t2) {
		garch_h[i] = NADBL;
	    } else {
		garch_h[i] = h[i + pad] * scale * scale;
	    }
	}
	gretl_model_set_data(pmod, "garch_h", garch_h,
			     GretlType.double_array,
			     dset.n * (*garch_h).sizeof);
    }

    return err;
}

private int make_garch_dataset(const(int) * list, const(DatasetRaw) * dset, int bign, int pad, int nx, double ** py, double *** pX)
{
    double * y = null;
    double ** X = null;
    int vx;
    int vy = list[4];
    int i;
    int k;
    int s;
    int t;

    /* If pad > 0 we have to create a newly allocated, padded
       dataset.  Otherwise we can use a virtual dataset, made
       up of pointers into the original dataset, Z.
    */

    if (pad > 0) {
	y = cast(double *) malloc(bign * (*y).sizeof);
	if (y == null) {
	    return GretlError.alloc;
	}
	*py = y;
    }

    if (nx > 0) {
	if (pad) {
	    X = doubles_array_new(nx, bign);
	} else {
	    X = cast(double* *) malloc(nx * (*X).sizeof);
	}

	if (X == null) {
	    free(y);
	    *py = null;
	    return GretlError.alloc;
	}
    }

    if (pad > 0) {
	/* build padded dataset */
	for (t=0; t<bign; t++) {
	    if (t < pad) {
		y[t] = 0.0;
		for (i=0; i<nx; i++) {
		    X[i][t] = 0.0;
		}
	    } else {
		s = t - pad;
		y[t] = dset.Z[vy][s];
		k = 5;
		for (i=0; i<nx; i++) {
		    vx = list[k++];
		    X[i][t] = dset.Z[vx][s];
		}
	    }
	}
    } else {
	/* build virtual dataset */
	*py = cast(double *) dset.Z[vy];
	k = 5;
	for (i=0; i<nx; i++) {
	    vx = list[k++];
	    X[i] = cast(double *) dset.Z[vx];
	}
    }

    *pX = X;

    return 0;
}

private int get_vopt(int robust)
{
    int vopt = libset_get_int(SetKey.garch_vcv);
    int ropt = libset_get_int(SetKey.garch_alt_vcv);

    /* The defaults: QML if "robust" option is in force,
       otherwise negative Hessian */
    if (vopt == MLVCVType.unset) {
	if (robust) {
	    if (ropt == MLVCVType.unset) {
		vopt = MLVCVType.qml;
	    } else {
		vopt = ropt;
	    }
	} else {
	    vopt = MLVCVType.hessian;
	}
    }

    return vopt;
}

private void garch_print_init(const(double) * theta, int k, int p, int q, int manual, Prn * prn)
{
    int i;
    int j = 0;

    pputc(prn, '\n');

    if (manual) {
	pputs(prn, _("Manual initialization of parameters").ptr);
    } else {
	pputs(prn, _("Automatic initialization of parameters").ptr);
    }

    pprintf(prn, "\n\n %s:\n", _("Regression coefficients").ptr);

    for (i=0; i<k; i++) {
	pprintf(prn, "  theta[%d] = %g\n", i, theta[j++]);
    }

    pprintf(prn, "\n %s:\n", _("Variance parameters").ptr);

    pprintf(prn, "  alpha[0] = %g\n", theta[j++]);
    for (i=0; i<q; i++) {
	pprintf(prn, "  alpha[%d] = %g\n", i+1, theta[j++]);
    }
    for (i=0; i<p; i++) {
	pprintf(prn, "   beta[%d] = %g\n", i, theta[j++]);
    }

    pputc(prn, '\n');
}

/* pick up any manually set initial values (if these
   have been set via "set initvals")
*/

private int garch_manual_init(double * theta, int k, int p, int q, Gretlopt opt, Prn * prn)
{
    int mlen = n_initvals();
    int n = k + p + q + 1;

    if (mlen != n) {
	if (mlen > 0) {
	    fprintf(stderr, "Number of initvals = %d, but we want %d " ~
			"values for GARCH\n", mlen, n);
	}
	/* initialization not done */
	return 0;
    }

    /* if we're _not_ using FCP, the following is handled
       within the BFGS routine */

    if (opt & Gretlopt.OPT_F) {
	GretlMatrixRaw * m = get_initvals();
	int i;

	/* order: coeffs on regressors; variance params */
	for (i=0; i<n; i++) {
	    theta[i] = m.val[i];
	}

	garch_print_init(theta, k, p, q, 1, prn);
	gretl_matrix_free(m);
    }

    return 1;
}

private int garch_peek_manual_init(int k, int p, int q)
{
    return k + p + q + 1 == n_initvals();
}

private int garch_driver(const(int) * list, double scale, const(DatasetRaw) * dset, Model * pmod, double * vparm, int ifc, Gretlopt opt, Prn * prn)
{
    int t1 = pmod.t1;
    int t2 = pmod.t2;
    int nc = pmod.ncoeff;
    int p = list[1];
    int q = list[2];
    double * y = null;
    double **X = null;
    double * h = null;
    double * e = null;
    double * e2 = null;
    double * theta = null;
    double ll = NADBL;
    GretlMatrixRaw * V = null;
    int fnc = 0;
    int grc = 0;
    int iters = 0;
    int nobs;
    int maxlag;
    int bign;
    int pad = 0;
    int i;
    int npar;
    int vopt;
    int err = 0;

    vopt = get_vopt(opt & Gretlopt.OPT_R);

    maxlag = (p > q)? p : q;
    npar = nc + p + q + 1;

    nobs = t2 + 1; /* number of obs in full dataset */

    if (maxlag > t1) {
	/* need to pad data series at start */
	pad = maxlag - t1;
    }

    /* length of series to pass to garch_estimate */
    bign = nobs + pad;

    e = cast(double *) malloc(bign * (*e).sizeof);
    e2 = cast(double *) malloc(bign * (*e2).sizeof);
    h = cast(double *) malloc(bign * (*h).sizeof);
    if (e == null || e2 == null || h == null) {
	err = GretlError.alloc;
	goto bailout;
    }

    for (i=0; i<bign; i++) {
	e[i] = e2[i] = h[i] = 0.0;
    }

    theta = cast(double *) malloc(npar * (*theta).sizeof);
    if (theta == null) {
	err = GretlError.alloc;
	goto bailout;
    }

    V = gretl_zero_matrix_new(npar, npar);
    if (V == null) {
	err = GretlError.alloc;
	goto bailout;
    }

    /* create dataset for garch estimation */
    err = make_garch_dataset(list, dset, bign, pad, nc, &y, &X);
    if (err) {
	goto bailout;
    }

    if (!garch_manual_init(theta, nc, p, q, opt, prn)) {
	/* initial coefficients from OLS */
	for (i=0; i<nc; i++) {
	    theta[i] = pmod.coeff[i];
	}
	/* initialize variance parameters */
	for (i=0; i<p+q+1; i++) {
	    theta[i+nc] = vparm[i];
	}
	if (opt & Gretlopt.OPT_V) {
	    garch_print_init(theta, nc, p, q, 0, prn);
	}
    }

    if (opt & Gretlopt.OPT_F) {
	/* --fcp */
	err = garch_estimate(y, cast(const(double) **) X,
			     t1 + pad, t2 + pad, bign, nc,
			     p, q, theta, V, e, e2, h,
			     scale, &ll, &iters, vopt, prn);
    } else {
	err = garch_estimate_mod(y, cast(const(double) **) X,
				 t1 + pad, t2 + pad, bign, nc,
				 p, q, theta, V, e, e2, h,
				 scale, &ll, &fnc, &grc, vopt, prn);
    }

    if (!err) {
	pmod.lnL = ll;
	write_garch_stats(pmod, list, dset, theta, V, p, q,
			  scale, e, h, npar, nc, pad, ifc, prn);
	if (iters > 0) {
	    gretl_model_set_int(pmod, "iters", iters);
	} else if (grc > 0) {
	    gretl_model_set_int(pmod, "fncount", fnc);
	    gretl_model_set_int(pmod, "grcount", grc);
	} else {
	    gretl_model_set_int(pmod, "iters", fnc);
	}
	gretl_model_set_vcv_info(pmod, VCVMajorType.ml, vopt);
	if (opt & Gretlopt.OPT_F) {
	    pmod.opt |= Gretlopt.OPT_F;
	}
    }

 bailout:

    free(e);
    free(e2);
    free(h);
    free(theta);
    gretl_matrix_free(V);

    if (pad > 0) {
	free(y);
	doubles_array_free(X, nc);
    } else {
	free(X);
    }

    if (err && !pmod.errcode) {
	pmod.errcode = err;
    }

    return err;
}

private int add_uhat_squared(const(Model) * pmod, double scale, DatasetRaw * dset)
{
    int t;
    int v = dset.v;

    if (dataset_add_series(dset, 1)) {
	return GretlError.alloc;
    }

    for (t=0; t<dset.n; t++) {
	double u = pmod.uhat[t];

	if (na(u)) {
	    dset.Z[v][t] = NADBL;
	} else {
	    u /= scale;
	    dset.Z[v][t] = u * u;
	}
    }

    strcpy(dset.varname[v], "uhat2");

    return 0;
}

/*
  p and q are the GARCH orders
  ao = max(q,p) is the ar order
  mo = q is the ma order

  it is assumed that armapar contains the arma parameters
  in the following order:
  armapar[0] : intercept
  armapar[1..ao] : ar terms
  armapar[ao+1..ao+mo] : ma terms
*/

private void garchpar_from_armapar(const(double) * armapar, int q, int p, double * vparm)
{
    double x;
    double sum_ab = 0.0;
    int ao = (p > q)? p : q;
    int mo = q;
    int i;


    for (i=1; i<=p; i++) {
	x = 0.0;
	if (i <= ao) {
	    x += armapar[i];
	}
	if (i<=mo) {
	    x += armapar[p+i];
	}
	vparm[i] = (x < 0.0)? 0.01 : x;
	sum_ab += vparm[i];
    }

    for (i=1; i<=q; i++) {
	x = armapar[p+i];
	vparm[p+i] = (x > 0.0)? 0.0001 : -x;
	sum_ab += vparm[p+i];
    }


    if (sum_ab > GARCH_PARAM_MAX) {
	for (i=1; i<=p+q; i++) {
	    vparm[i] *= GARCH_PARAM_MAX / sum_ab;
	}
	sum_ab = GARCH_PARAM_MAX;
    }

    vparm[0] = armapar[0];
}

private int garch_init_by_arma(const(Model) * pmod, const(int) * glist, DatasetRaw * dset, double scale, double * vparm)
{
    int p = glist[1];
    int q = glist[2];
    int v = dset.v;
    int * list = null;
    int err = 0;

    /* for now we'll try this only for GARCH up to (2,2) */
    if (q > 2 || p > 2) {
 	return 0;
    }

    /* add OLS uhat squared to dataset */
    if (add_uhat_squared(pmod, scale, dset)) {
	return GretlError.alloc;
    }

    list = gretl_list_copy(glist);

    if (list == null) {
	err = GretlError.alloc;
    } else {
	Model amod;
	int i;

	list[1] = (q > p)? q : p;
	list[2] = q;
	/* dep var is squared OLS residual: last var added */
	list[4] = v;

	amod = arma(list, null, dset, Gretlopt.OPT_C, null);
	err = amod.errcode;
	if (!err) {
	    model_count_minus(&amod);
	    garchpar_from_armapar(amod.coeff, p, q, vparm);
	    for (i=0; i<q+p+1; i++) {
		fprintf(stderr, "from ARMA: vparm_init[%d] = %#12.6g\n", i,
			vparm[i]);
	    }
	}
	clear_model(&amod);
    }

    dataset_drop_last_variables(dset, dset.v - v);
    free(list);

    return err;
}

/* XPOS is the list position of the first regressor (if any).
   Note that the garch list structure is:

    0 1 2   3   4  5  6 ...
    # p q <sep> y x0 x1 ...

*/
enum XPOS = 5;

private int * get_garch_list(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt, int * ifc, int * err)
{
    int * glist = null;
    int i;
    int p;
    int q;
    int cpos = 0;
    int add0 = 0;

    /* is the list well-formed? */
    if (list[0] < 4 || list[1] == LISTSEP ||
	list[2] == LISTSEP || list[3] != LISTSEP) {
	*err = GretlError.parse;
	return null;
    }

    p = list[1];
    q = list[2];

    *err = 0;

    /* negative orders don't make sense */
    if (p < 0 || q < 0) {
	gretl_errmsg_set(_("GARCH: neither p nor q can be negative").ptr);
	*err = GretlError.data;
	return null;
    }

    /* rule out pure AR in variance: the model is unidentified */
    if (p > 0 && q == 0) {
	gretl_errmsg_set(_("GARCH: p > 0 and q = 0: the model is unidentified").ptr);
	*err = GretlError.data;
	return null;
    }

    /* rule out excessive total GARCH-iness */
    if (p + q > PQ_MAX) {
	gretl_errmsg_sprintf(_("GARCH: p + q must not exceed %d").ptr, PQ_MAX);
	*err = GretlError.data;
	return null;
    }

    /* check for presence of constant among regressors */
    for (i=XPOS; i<=list[0]; i++) {
	if (list[i] == 0) {
	    /* got the constant: OK */
	    cpos = i;
	    break;
	}
    }

    /* OPT_N means don't auto-add a constant */
    if (cpos == 0 && !(opt & Gretlopt.OPT_N)) {
	add0 = 1;
    }

    *ifc = (cpos > 0 || add0);

    glist = gretl_list_new(list[0] + add0);

    if (glist == null) {
	*err = GretlError.alloc;
    } else {
	int j = 1;

	/* transcribe first portion of original list */
	for (i=1; i<XPOS; i++) {
	    glist[j++] = list[i];
	}

	if (add0 || (cpos > 0 && cpos != XPOS)) {
	    /* insert constant here if not already present,
	       or if originally placed later */
	    glist[j++] = 0;
	}

	/* transcribe the original regressors, if any */
	for (i=XPOS; i<=list[0]; i++) {
	    if (i == XPOS || list[i] != 0) {
		glist[j++] = list[i];
	    }
	}
    }

    return glist;
}



int garch_pretest(Model * pmod, DatasetRaw * dset, double * LMF, double * pvF)
{
    int err;

    err = autocorr_test(pmod, dset.pd, dset,
			Gretlopt.OPT_S | Gretlopt.OPT_Q, null);

    if (!err) {
	*LMF = get_last_test_statistic();
	*pvF = get_last_pvalue();
    }

    return err;
}

private void autocorr_message(double LMF, double pvF, int order, Prn * prn)
{
    if (!na(LMF) && pvF < 0.05) {
	pputs(prn, _("\nConvergence was not reached.  One possible reason " ~
			"for this is\nautocorrelation in the error term.\n").ptr);
	pprintf(prn, _("After estimating the model by OLS, the following result " ~
			"was\nobtained for a test of autocorrelation of order %d:\n").ptr,
		order);
	pprintf(prn, "LMF = %g, with p-value %g\n", LMF, pvF);
    }
}




private double garch_scale_sigma(int yno, Model * pmod, DatasetRaw * dset)
{
    double scale = pmod.sigma;
    int i;

    for (i=0; i<dset.n; i++) {
	if (!na(dset.Z[yno][i])) {
	    dset.Z[yno][i] /= scale;
	}
    }

    for (i=0; i<pmod.ncoeff; i++) {
	pmod.coeff[i] /= scale;
    }

    pmod.ess /= scale * scale;
    pmod.sigma = 1.0;

    return scale;
}

private void garch_undo_scaling(int yno, double scale, DatasetRaw * dset)
{
    int t;

    if (scale != 1.0) {
	for (t=0; t<dset.n; t++) {
	    if (!na(dset.Z[yno][t])) {
		dset.Z[yno][t] *= scale;
	    }
	}
    }
}


/* default variance parameter initialization */

private void garch_vparm_init(const(int) * list, double sigma, double * vparm)
{
    int i;
    int q = list[1];
    int p = list[2];
    double den = 1.0;
    double tmp = (q>0) ? 0.2 : 0.8;

    if (p > 0) {
	for (i=1; i<=p; i++) {
	    vparm[i] = tmp / p;
	    den -= vparm[i];
	}
    }

    if (q > 0) {
	for (i=p+1; i<=p+q; i++) {
	    vparm[i] = 0.7 / q;
	    den -= vparm[i];
	}
    }

    vparm[0] = sigma * sigma * den;
}

/* make regression list for initial OLS: we skip three terms
   from the GARCH list, namely p, q and the separator
   that follows.
*/

private int * make_ols_list(const(int) * list, int * err)
{
    int * olist;
    int i;

    olist = gretl_list_new(list[0] - 3);

    if (olist == null) {
	*err = GretlError.alloc;
    } else {
	for (i=4; i<=list[0]; i++) {
	    olist[i-3] = list[i];
	}
    }

    return olist;
}

private Model garch_run_ols(const(int) * list, DatasetRaw * dset, Prn * prn)
{
    int * ols_list;
    Model model;
    int err = 0;

    ols_list = make_ols_list(list, &err);
    if (err) {
	gretl_model_init(&model, null);
	model.errcode = err;
	return model;
    }

    model = lsq(ols_list, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_M | Gretlopt.OPT_U);


    free(ols_list);

    if (!model.errcode) {
	free(model.xlist);
	model.xlist = gretl_list_sublist(model.list, 2, -1);
	clear_model_xpx(&model);
    }

    return model;
}

private void clean_dropped_vars(Model mod, int * list)
{
    if (list[0] - mod.list[0] > 3) {
	int i;

        list[0] = mod.list[0] + 3;
	for (i=4; i<=list[0]; i++) {
	    list[i] = mod.list[i-3];
	}
    }
}

private void garch_add_lr_test(Model * pmod, double llr, const(int) * list)
{
    if (!na(pmod.lnL) && llr <= pmod.lnL) {
	double LR = 2.0 * (pmod.lnL - llr);
	int LRdf = list[1] + list[2];

	gretl_model_set_double(pmod, "garch_LR", LR);
	gretl_model_set_int(pmod, "garch_LR_df", LRdf);
    }
}

private void garch_standardize_residuals(Model * pmod)
{
    double * h = cast(double *) gretl_model_get_data(pmod, "garch_h");

    if (h != null) {
	int t;

	for (t=pmod.t1; t<=pmod.t2; t++) {
	    pmod.uhat[t] /= sqrt(h[t]);
	}
	pmod.opt |= Gretlopt.OPT_Z;
    }
}

/* the driver function for the plugin */

extern (C) Model garch_model(const(int) * cmdlist, DatasetRaw * dset, Prn * prn, Gretlopt opt)
{
    Model model;
    int * list = null;
    double[PQ_MAX+1] vparm = 0;
    double LMF = NADBL;
    double pvF = NADBL;
    double llr = NADBL;
    double scale = 1.0;
    int ols_T;
    int ifc;
    int yno = 0;
    int have_init = 0;
    int err = 0;

    list = get_garch_list(cmdlist, dset, opt, &ifc, &err);
    if (err) {
	gretl_model_init(&model, null);
	model.errcode = err;
	return model;
    }

    /* run initial OLS */
    model = garch_run_ols(list, dset, prn);
    if (model.errcode) {
	free(list);
	return model;
    }

    clean_dropped_vars(model, list);
    llr = model.lnL;
    ols_T = model.nobs;

    have_init = garch_peek_manual_init(model.ncoeff, list[1], list[2]);

    if (prn != null) {
	garch_pretest(&model, dset, &LMF, &pvF);
    }

    if (!have_init) {
	yno = list[4];
	scale = garch_scale_sigma(yno, &model, dset);
    }

    if (!have_init) {
	/* variance parameter initialization */
	garch_vparm_init(list, model.sigma, vparm.ptr);
	if (opt & Gretlopt.OPT_A) {
	    /* "--arma-init": try initializing params via ARMA */
	    garch_init_by_arma(&model, list, dset, scale, vparm.ptr);
	}
    }

    garch_driver(list, scale, dset, &model, vparm.ptr,
		 ifc, opt, prn);

    garch_undo_scaling(yno, scale, dset);

    if (!na(LMF)) {
	if (model.errcode == GretlError.noconv) {
	    autocorr_message(LMF, pvF, dset.pd, prn);
	} else {
	    gretl_model_destroy_tests(&model);
	}
    }

    if (!model.errcode) {
	if (opt & Gretlopt.OPT_Z) {
	    garch_standardize_residuals(&model);
	}
	if (!na(llr) && ols_T == model.nobs) {
	    garch_add_lr_test(&model, llr, cmdlist);
	}
    }

    free(list);

    return model;
}
