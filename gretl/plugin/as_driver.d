/* Port of plugin/as_driver.c: exact ML for ARMA models by the AS 154 (Gardner et al.) or AS 197 (Melard, gretl.plugin.as197) algorithms:
   the shared workspace (AsInfo, C as_info), the likelihood iterations and the BFGS driver as_arma, which arma.c calls, the per-observation
   likelihood callbacks used for the OPG matrix, and the finishing (undoing of y scaling, coefficient and model statistics).

   In the C this file and as197.c / as154.c are #included by arma.c; here AS 197 is gretl.plugin.as197 and AS 154 itself (starma, karma)
   is gretl.plugin.as154. as_arma, as154_llt_callback and as197_llt_callback are
   extern (C) and public. */

module gretl.plugin.as_driver;

import gretl.gretl_bfgs : BFGS_defaults, BFGS_max, numerical_hessian_inverse;
import gretl.gretl_matrix : gretl_identity_matrix_new, gretl_matrix_alloc, gretl_matrix_free, gretl_matrix_qform;
import gretl.gretl_model : gretl_model_allocate_storage, gretl_model_set_int, gretl_model_set_vcv_info, gretl_model_write_vcv;
import gretl.gretl_utils : copyvec;
import gretl.matrix_extra : gretl_vcv_matrix_from_model;
import gretl.missingh : NADBL, na;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import gretl.plugin.arma_priv;
import gretl.plugin.arma : arima_ydiff_only, arma_OPG_vcv, arma_QML_vcv, arma_model_add_roots, arma_use_opg, form_arma_X_matrix, form_arma_y_vector, maybe_correct_MA, write_arma_model_stats;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_modelh : MLVCVType, VCVMajorType;
import gretl.nlsh : CriterionType;
import gretl.gretl_bfgsh : BFGS_CRIT_FUNC;
import gretl.plugin.as197 : flikam;
import gretl.plugin.as154 : karma, starma;
import gretl.gretl_prn : Prn;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_modelh : ArmaFlags;


/* shared between AS 154 and AS 197 */

struct AsInfo {
    int algo;
    int p;
    int P;
    int q;
    int Q;
    int pd;
    int plen;
    int qlen;
    int r;
    int rp1; /* AS 197 */
    int np;
    int nrbar; /* AS 154 */
    int ifault;
    int n;
    int ok_n;
    int ifc;
    /* AR and MA coeffs */
    double * phi;
    double * theta;
     /* dependent var and forecast errors */
    double * y;
    double * y0;
    double * e;
    /* AS 197 workspace */
    double * vw;
    double * vl;
    double * vk;
    /* AS 154 workspace */
    double * A;
    double * P0;
    double * V;
    double * thetab;
    double * xnext;
    double * xrow;
    double * rbar;
    double * evec;
    /* components of likelihood */
    double sumsq;
    double fact;
    double sumlog;
    double toler; /* tolerance for switching to fast iterations */
    double loglik;
    BFGS_CRIT_FUNC cfunc;
    int ma_check;
    int iupd; /* specific to AS 154 */
    int ncalls;
    ArmaInfo *ai;
    GretlMatrixRaw * X;
    int free_X;
};

private int as_154_alloc(AsInfo * as)
{
    int err = 0;

    /* unused pointers specific to AS 197 */
    as.vw = as.vl = as.vk = null;

    as.phi =   cast(double *) malloc(as.r * (*as.phi).sizeof);
    as.theta = cast(double *) malloc(as.r * (*as.theta).sizeof);
    as.A =     cast(double *) malloc(as.r * (*as.A).sizeof);
    as.P0 =    cast(double *) malloc(as.np * (*as.P0).sizeof);
    as.V =     cast(double *) malloc(as.np * (*as.V).sizeof);
    as.e =     cast(double *) malloc(as.n * (*as.e).sizeof);
    as.evec =  cast(double *) malloc(as.r * (*as.evec).sizeof);

    if (as.phi == null || as.theta == null || as.A == null ||
	as.P0 == null || as.V == null || as.e == null ||
	as.evec == null) {
	err = GretlError.alloc;
    }

    if (!err) {
	int worklen = 3 * as.np + as.nrbar;

	as.thetab = cast(double *) malloc(worklen * (*as.thetab).sizeof);
	if (as.thetab == null) {
	    err = GretlError.alloc;
	} else {
	    as.xnext = as.thetab + as.np;
	    as.xrow = as.xnext + as.np;
	    as.rbar = as.xrow + as.np;
	}
    }

    return err;
}

private int as_197_alloc(AsInfo * as)
{
    int err = 0;

    /* unused pointers specific to AS 154 */
    as.A = as.P0 = as.V = as.evec = null;
    as.thetab = as.xnext = as.xrow = as.rbar = null;

    as.phi = as.theta = null;

    if (as.plen > 0) {
	as.phi = cast(double *) malloc(as.plen * (*as.phi).sizeof);
	if (as.phi == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err && as.qlen > 0) {
	as.theta = cast(double *) malloc(as.qlen * (*as.theta).sizeof);
	if (as.theta == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err) {
	int worklen = as.n + 2*as.rp1 + as.r;

	as.e =  cast(double *) malloc(worklen * (*as.e).sizeof);
	if (as.e == null) {
	    err = GretlError.alloc;
	} else {
	    as.vw = as.e + as.n;
	    as.vl = as.vw + as.rp1;
	    as.vk = as.vl + as.rp1;
	}
    }

    return err;
}

private int as_info_init(AsInfo * as, int algo, ArmaInfo * ai, double toler)
{
    int err = 0;

    as.algo = algo;
    as.ai = ai; /* create accessor */

    /* convenience copies of @ai integer values */
    as.p = ai.p;
    as.P = ai.P;
    as.q = ai.q;
    as.Q = ai.Q;
    as.pd = ai.pd;
    as.n = ai.fullT;
    as.ok_n = ai.T;
    as.ifc = ai.ifc;

    as.plen = as.p + as.pd * as.P;
    as.qlen = as.q + as.pd * as.Q;
    as.r = (as.plen > as.qlen + 1)? as.plen : as.qlen + 1;

    if (algo == 154) {
	as.np = as.r * (as.r + 1)/2;
	as.nrbar = as.np * (as.np - 1)/2;
	as.rp1 = 0;
    } else {
	as.rp1 = as.r + 1;
	as.np = as.nrbar = 0;
    }

    as.y = as.y0 = null; /* later! */
    as.X = null; /* later too */
    as.free_X = 0;

    if (algo == 154) {
	err = as_154_alloc(as);
    } else {
	err = as_197_alloc(as);
    }

    if (!err) {
	as.toler = toler;
	as.loglik = NADBL;
	as.ifault = 0;
	as.ma_check = 0;
	as.iupd = 1; /* AS 154: FIXME AR(1) */
	as.ncalls = 0;
    }

    return err;
}

private void as_info_free(AsInfo * as)
{
    free(as.phi);
    free(as.theta);
    free(as.e);
    free(as.y0);

    if (as.algo == 154) {
	free(as.A);
	free(as.P0);
	free(as.V);
	free(as.evec);
	free(as.thetab);
    }

    if (as.free_X) {
	gretl_matrix_free(as.X);
    }
}

private void as_write_big_phi(const(double) * b, AsInfo * as)
{
    const(double) * bs = b + as.ai.np;
    double x;
    double y;
    int i;
    int j;
    int k;
    int ii;

    for (i=0; i<as.plen; i++) {
	as.phi[i] = 0.0;
    }

    for (j=-1; j<as.P; j++) {
	x = (j < 0)? -1 : bs[j];
	k = 0;
        for (i=-1; i<as.p; i++) {
	    if (i < 0) {
		y = -1;
	    } else if (AR_included(as.ai, i)) {
		y = b[k++];
	    } else {
		y = 0.0;
	    }
            ii = (j+1) * as.pd + (i+1);
	    if (ii > 0) {
		as.phi[ii-1] -= x * y;
	    }
        }
    }
}

private void as_write_big_theta(const(double) * b, AsInfo * as)
{
    const(double) * bs = b + as.ai.nq;
    double x;
    double y;
    int i;
    int j;
    int k;
    int ii;

    for (i=0; i<as.qlen; i++) {
	as.theta[i] = 0.0;
    }

    for (j=-1; j<as.Q; j++) {
	x = (j < 0)? 1 : bs[j];
	k = 0;
        for (i=-1; i<as.q; i++) {
	    if (i < 0) {
		y = 1;
	    } else if (MA_included(as.ai, i)) {
		y = b[k++];
	    } else {
		y = 0.0;
	    }
            ii = (j+1) * as.pd + (i+1);
	    if (ii > 0) {
		as.theta[ii-1] += x * y;
	    }
        }
    }
}

private void as_fill_arrays(AsInfo * as, const(double) * b)
{
    int np = as.ai.np + as.P;
    int nq = as.ai.nq + as.Q;
    double mu = 0.0;
    int i;
    int j;

    if (as.ifc) {
	mu = b[0];
	if (as.ai.nexo == 0) {
	    /* just subtract the constant */
	    for (i=0; i<as.n; i++) {
		as.y[i] = as.y0[i];
		if (!isnan(as.y0[i])) {
		    as.y[i] -= mu;
		}
	    }
	}
	b++;
    }

    if (as.P > 0) {
	as_write_big_phi(b, as);
    } else if (as.p > 0) {
	j = 0;
	for (i=0; i<as.p; i++) {
	    if (AR_included(as.ai, i)) {
		as.phi[i] = b[j++];
	    } else {
		as.phi[i] = 0.0;
	    }
	}
    }
    b += np;

    if (as.Q > 0) {
	as_write_big_theta(b, as);
    } else if (as.q > 0) {
	j = 0;
	for (i=0; i<as.q; i++) {
	    if (MA_included(as.ai, i)) {
		as.theta[i] = b[j++];
	    } else {
		as.theta[i] = 0.0;
	    }
	}
    }
    b += nq;

    if (as.ai.nexo > 0) {
	/* subtract the regression effect */
	double xij;

	for (i=0; i<as.n; i++) {
	    as.y[i] = as.y0[i];
	    if (!isnan(as.y[i])) {
		if (as.ifc) {
		    as.y[i] -= mu;
		}
		for (j=0; j<as.ai.nexo; j++) {
		    xij = gretl_matrix_get(as.X, i, j);
		    as.y[i] -= xij * b[j];
		}
	    }
	}
    }
}

/* full ARMA loglikelihood */

private double as_loglikelihood(const(AsInfo) * as)
{
    double ll1 = 1.0 + LN_2_PI + log(as.sumsq / as.ok_n);

    if (as.algo == 154) {
	return -0.5 * (as.ok_n * ll1 + as.sumlog);
    } else {
	return -0.5 * as.ok_n * (ll1 + log(as.fact));
    }
}

extern (C) double as197_iteration(const(double) * b, void * data)
{
    AsInfo * as = cast(AsInfo *) data;
    double crit = NADBL;
    /* number of actually included AR terms */
    int np = as.ai.np + as.P;

    as.ncalls += 1;

    if (as.ma_check) {
	/* check that MA term(s) are within bounds */
	double *theta = cast(double *) b + as.ifc + np;
	double * Theta = theta + as.ai.nq;

	if (maybe_correct_MA(as.ai, theta, Theta)) {
	    return NADBL;
	}
    }

    as_fill_arrays(as, b);

    as.ifault = flikam(as.phi, as.plen, as.theta, as.qlen,
			as.y, as.e, as.n, &as.sumsq, &as.fact,
			as.vw, as.vl, as.rp1, as.vk, as.r,
			as.toler);

    if (as.ifault > 0) {
	if (as.ifault == 5) {
	    ; // fputs("flikam: (near) non-stationarity\n", stderr);
	} else {
	    fprintf(stderr, "flikam: ifault = %d\n", as.ifault);
	}
	return NADBL;
    }

    if (isnan(as.sumsq) || isnan(as.fact)) {
	; /* leave crit as NA */
    } else {
	as.loglik = crit = as_loglikelihood(as);
    }

    return crit;
}

extern (C) double as154_iteration(const(double) * b, void * data)
{
    AsInfo * as = cast(AsInfo *) data;
    double crit = NADBL;
    /* number of actually included AR terms */
    int np = as.ai.np + as.P;
    int nit = 0;

    if (as.ma_check) {
	/* check that MA term(s) are within bounds */
	double *theta = cast(double *) b + as.ifc + np;
	double * Theta = theta + as.ai.nq;

	if (maybe_correct_MA(as.ai, theta, Theta)) {
	    return NADBL;
	}
    }

    as_fill_arrays(as, b);

    as.ifault = starma(as.plen, as.qlen, as.r, as.np,
			as.phi, as.theta, as.A, as.P0, as.V,
			as.thetab, as.xnext, as.xrow,
			as.rbar, as.nrbar);

    if (as.ifault) {
	fprintf(stderr, "starma: ifault = %d\n", as.ifault);
	return NADBL;
    }

    /* initialization required */
    as.sumlog = as.sumsq = 0;

    karma(as.plen, as.qlen, as.r, as.np,
	  as.phi, as.theta, as.A, as.P0, as.V,
	  as.n, as.y, as.e,
	  &as.sumlog, &as.sumsq, as.iupd,
	  as.toler, as.evec, &nit);

    if (isnan(as.sumlog) || isnan(as.sumsq) || as.sumsq <= 0) {
	; // fprintf(stderr, "karma: got NaNs, nit = %d\n", nit);
    } else {
	as.loglik = crit = as_loglikelihood(as);
    }

    return crit;
}

extern (C) const(double) * as197_llt_callback(const(double) * b, int i, void * data)
{
    AsInfo * as = cast(AsInfo *) data;
    int err;

    as_fill_arrays(as, b);
    err = flikam(as.phi, as.plen, as.theta, as.qlen,
		 as.y, as.e, as.n, &as.sumsq, &as.fact,
		 as.vw, as.vl, as.rp1, as.vk, as.r,
		 as.toler);

    return (err)? null : as.e;
}

extern (C) const(double) * as154_llt_callback(const(double) * b, int i, void * data)
{
    AsInfo * as = cast(AsInfo *) data;
    int err = 0;
    int nit = 0;

    as_fill_arrays(as, b);
    as.ifault = starma(as.plen, as.qlen, as.r, as.np,
			as.phi, as.theta, as.A, as.P0, as.V,
			as.thetab, as.xnext, as.xrow,
			as.rbar, as.nrbar);

    as.sumlog = as.sumsq = 0;
    karma(as.plen, as.qlen, as.r, as.np,
	  as.phi, as.theta, as.A, as.P0, as.V,
	  as.n, as.y, as.e,
	  &as.sumlog, &as.sumsq, as.iupd,
	  as.toler, as.evec, &nit);

    if (isnan(as.sumlog) || isnan(as.sumsq) || as.sumsq <= 0) {
	fprintf(stderr, "as154_llt_callback: failed\n");
	err = GretlError.nan;
    }

    return (err)? null : as.e;
}

/* Undo scalings, allowing for standardization of exogenous
   regressors. At present we come here only if such
   standardization has been done.
*/

private int unscramble_scalings(ArmaInfo * ainfo, Model * pmod)
{
    GretlMatrixRaw * difmat;
    GretlMatrixRaw * V;
    GretlMatrixRaw * V0;
    double * xbar;
    double * sdx;
    double * b;
    double sdy;
    int nc = ainfo.nc;
    int nx = ainfo.nexo;
    int i;
    int j;
    int xpos;
    int err = 0;

    difmat = gretl_identity_matrix_new(nc);
    V = gretl_matrix_alloc(nc, nc);

    if (difmat == null || V == null) {
	gretl_matrix_free(difmat);
	gretl_matrix_free(V);
	return GretlError.alloc;
    }

    sdy = 1 / ainfo.yscale;
    b = pmod.coeff;
    V0 = gretl_vcv_matrix_from_model(pmod, null, &err);
    if (err) {
	goto bailout;
    }

    xbar = ainfo.xstats.val;
    sdx = xbar + nx;
    xpos = nc - nx;
    for (i=xpos, j=0; i<nc; i++, j++) {
	b[i] = sdy * b[i] / sdx[j];
    }

    if (ainfo.ifc && ainfo.yscale != 1.0) {
	b[0] /= ainfo.yscale;
	b[0] += ainfo.yshift;
    }

    difmat.val[0] = sdy;

    for (i=0; i<nx; i++) {
	j = xpos + i;
	gretl_matrix_set(difmat, j, j, sdy / sdx[i]);
    }
    for (i=0; i<nx; i++) {
	if (ainfo.ifc) {
	    b[0] -= xbar[i] * b[xpos+i];
	}
	/* the following line also conditional on ifc? */
	gretl_matrix_set(difmat, 0, xpos+i, -sdy * xbar[i] / sdx[i]);
    }

    err = gretl_matrix_qform(difmat, GretlMatrixMod.none,
			     V0, V, GretlMatrixMod.none);
    if (!err) {
	err = gretl_model_write_vcv(pmod, V);
    }

 bailout:

    gretl_matrix_free(difmat);
    gretl_matrix_free(V);
    gretl_matrix_free(V0);

    return err;
}

/* Undo y scaling, on the assumption that standardization of
   exogenous regressors has NOT been done.
*/

private int as_undo_y_scaling(ArmaInfo * ainfo, GretlMatrixRaw * y, double * b, AsInfo * as)
{
    double *beta = b + ainfo.ifc + ainfo.np + ainfo.P +
	ainfo.nq + ainfo.Q;
    double lnl;
    int i;
    int t;
    int err = 0;

    if (ainfo.ifc) {
	b[0] /= ainfo.yscale;
	b[0] += ainfo.yshift;
    }

    for (i=0; i<ainfo.nexo; i++) {
	beta[i] /= ainfo.yscale;
    }

    for (t=0; t<ainfo.fullT; t++) {
	/* restore original y for loglikelihood calculation */
	if (!isnan(as.y[t])) {
	    as.y[t] /= ainfo.yscale;
	    as.y[t] += ainfo.yshift;
	    if (as.y0 != null) {
		as.y0[t] /= ainfo.yscale;
		as.y0[t] += ainfo.yshift;
	    }
	}
    }

    lnl = as.cfunc(b, as);

    if (na(lnl)) {
	err = 1;
    }

    return err;
}

private int as_arma_finish(Model * pmod, ArmaInfo * ainfo, const(DatasetRaw) * dset, AsInfo * as, double * b, Gretlopt opt, Prn * prn)
{
    int i;
    int t;
    int k = ainfo.nc;
    int do_opg = arma_use_opg(opt);
    int QML = (opt & Gretlopt.OPT_R);
    int vcv_err = 0;
    double s2;
    int err;

    pmod.t1 = ainfo.t1;
    pmod.t2 = ainfo.t2;
    pmod.nobs = ainfo.T;
    pmod.ncoeff = ainfo.nc;
    pmod.full_n = dset.n;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    for (i=0; i<k; i++) {
	pmod.coeff[i] = b[i];
    }

    s2 = 0.0;
    i = 0;
    for (t=pmod.t1; t<=pmod.t2; t++) {
	if (isnan(as.e[i])) {
	    pmod.uhat[t] = NADBL;
	} else {
	    s2 += as.e[i] * as.e[i];
	    pmod.uhat[t] = as.e[i];
	}
	i++;
    }

    s2 /= ainfo.T;
    pmod.sigma = sqrt(s2);
    pmod.lnL = as.loglik;

    /* configure for computing variance matrix */
    as.ma_check = 0;

    if (!do_opg) {
	/* base covariance matrix on Hessian (perhaps QML) */
	GretlMatrixRaw * Hinv;
	double d = 0.0; /* adjust? */

	Hinv = numerical_hessian_inverse(b, ainfo.nc, as.cfunc,
					 as, d, &vcv_err);
	if (!vcv_err) {
	    if (QML) {
		vcv_err = arma_QML_vcv(pmod, Hinv, as, as.algo, b, s2,
				       k, ainfo.T, prn);
	    } else {
		err = gretl_model_write_vcv(pmod, Hinv);
		if (!err) {
		    gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.hessian);
		}
	    }
	} else if (!(opt & Gretlopt.OPT_H)) {
	    /* fallback when Hessian not explicitly requested */
	    do_opg = 1;
	    gretl_model_set_int(pmod, "hess-error", 1);
	}
	gretl_matrix_free(Hinv);
    }

    if (do_opg) {
	vcv_err = arma_OPG_vcv(pmod, as, as.algo, b, s2, k, ainfo.T, prn);
	if (!vcv_err) {
	    gretl_model_set_vcv_info(pmod, VCVMajorType.ml, MLVCVType.op);
	    pmod.opt |= Gretlopt.OPT_G;
	}
    }

    if (!err && arma_stdx(ainfo) && ainfo.yscale != 1.0) {
	pmod.lnL -= ainfo.T * log(1/ainfo.yscale);
    }

    if (!err) {
	write_arma_model_stats(pmod, ainfo, dset);
	arma_model_add_roots(pmod, ainfo, b);
	gretl_model_set_int(pmod, "arma_flags", ArmaFlags.exact);
	gretl_model_set_int(pmod, "as_algo", as.algo);
	if (arima_ydiff_only(ainfo)) {
	    pmod.opt |= Gretlopt.OPT_Y;
	}
    }

    if (arma_stdx(ainfo)) {
	/* maybe: ainfo.yscale != 1 || arma_stdx(ainfo) */
	unscramble_scalings(ainfo, pmod);
    }

    return err;
}

/* As of 2018-04, the AS 197 implementation for gretl
   can't properly handle missing values within the sample
   range; all other "special cases" should be OK.
*/

private int as197_ok(ArmaInfo * ainfo)
{
    return arma_missvals(ainfo) ? 0 : 1;
}

extern (C) int as_arma(const(double) * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod, Gretlopt opt)
{
    AsInfo as;
    GretlMatrixRaw * y = null;
    double * b = null;
    double toler = -1.0;
    int algo;
    int err = 0;

    memset(&as, 0, AsInfo.sizeof);

    if (opt & Gretlopt.OPT_A) {
	/* --as154 */
	algo = 154;
    } else {
	/* prefer AS 197 if it's usable */
	algo = as197_ok(ainfo) ? 197 : 154;
    }

    err = as_info_init(&as, algo, ainfo, toler);
    if (err) {
	return err;
    }

    b = copyvec(coeff, ainfo.nc);
    if (b == null) {
	return GretlError.alloc;
    }

    y = form_arma_y_vector(ainfo, &err);
    if (!err) {
	as.y = y.val;
	if (as.ifc || ainfo.nexo > 0) {
	    as.y0 = copyvec(as.y, as.n);
	    if (as.y0 == null) {
		err = GretlError.alloc;
	    }
	}
    }

    if (!err && ainfo.nexo > 0) {
	if (ainfo.dX != null) {
	    as.X = ainfo.dX;
	} else {
	    as.X = form_arma_X_matrix(ainfo, dset, &err);
	    as.free_X = 1;
	}
    }

    if (!err) {
	/* maximize loglikelihood via BFGS */
	Gretlopt maxopt = opt | (Gretlopt.OPT_A | Gretlopt.OPT_U);
	int maxit;
	double tol;                  /* the C redeclares toler here (shadowing) */

	if (as.algo == 197) {
	    as.cfunc = &as197_iteration;
	    if (as.n > 2000) {
		/* try to avoid slowdown on big samples? */
		as.toler = 0.0001;
	    }
	} else {
	    /* AS 154 */
	    as.cfunc = &as154_iteration;
	}

	if (as.q > 0 || as.Q > 0) {
	    as.ma_check = 1;
	}

	BFGS_defaults(&maxit, &tol, GretlCmdIndex.arma);

	err = BFGS_max(b, ainfo.nc, maxit, tol,
		       &ainfo.fncount, &ainfo.grcount,
		       as.cfunc, CriterionType.loglik, null, &as, null,
		       maxopt, ainfo.prn);

	if (!err) {
	    if (ainfo.yscale != 1.0 && !arma_stdx(ainfo)) {
		/* note: this implies recalculation of loglik */
		as_undo_y_scaling(ainfo, y, b, &as);
	    }
	    gretl_model_set_int(pmod, "fncount", ainfo.fncount);
	    gretl_model_set_int(pmod, "grcount", ainfo.grcount);
	    err = as_arma_finish(pmod, ainfo, dset, &as, b,
				 opt, ainfo.prn);
	}
    }

    if (err && !pmod.errcode) {
	pmod.errcode = err;
    }

    as_info_free(&as);
    gretl_matrix_free(y);
    free(b);

    return err;
}
