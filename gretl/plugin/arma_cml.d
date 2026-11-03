/* Port of plugin/arma_cml.c: ARMA estimation via conditional ML (BHHH method, gretl.bhhh_max): the analytical score, the conditional
   forecast errors and the BHHH callback, and bhhh_arma, the entry point arma.c calls.

   The function adjust_score_t1 (inside #if 0 in the C, never called) and the debug output (CML_DEBUG, ARMA_DEBUG: off) are not in the port. */

module gretl.plugin.arma_cml;

import gretl.bhhh_max : bhhh_max;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_zero_matrix_new;
import gretl.gretl_model : gretl_model_allocate_storage, gretl_model_set_int, gretl_model_write_vcv;
import gretl.gretl_prn : pputs;
import gretl.gretl_utils : doubles_array_new0;
import gretl.libset : libset_get_double;
import gretl.missingh : NADBL, na;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import gretl.plugin.arma_priv;
import gretl.plugin.arma : arma_list_y_position, arma_model_add_roots, maybe_correct_MA, write_arma_model_stats;
import gretl.libseth : SetKey;
import core.stdc.stdio : fputs;
import gretl.gretl_prn : Prn;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_modelh : ArmaFlags;



/* ARMA estimation via conditional ML (BHHH method) */



/* ln(sqrt(2*pi)) + 0.5 */
enum LN_SQRT_2_PI_P5 = 1.41893853320467274178;

private void do_MA_partials(double * drv, ArmaInfo * ainfo, const(double) * theta, const(double) * Theta, int t)
{
    int i;
    int j;
    int k;
    int p;
    int s;

    k = 0;
    for (i=0; i<ainfo.q; i++) {
	if (MA_included(ainfo, i)) {
	    p = i + 1;
	    if (t - p >= 0) {
		drv[0] -= theta[k] * drv[p];
	    }
	    k++;
	}
    }

    for (j=0; j<ainfo.Q; j++) {
	p = (j + 1) * ainfo.pd;
	if (t - p >= 0) {
	    drv[0] -= Theta[j] * drv[p];
	    k = 0;
	    for (i=0; i<ainfo.q; i++) {
		if (MA_included(ainfo, i)) {
		    s = p + i + 1;
		    if (t - s >= 0) {
			drv[0] -= Theta[j] * theta[k] * drv[s];
		    }
		    k++;
		}
	    }
	}
    }
}

/* for each of the arrays of derivatives, shuffle each
   value one place up */

private void push_derivs(ArmaInfo * ainfo, double ** de, int dlen)
{
    int i;
    int j;

    for (i=0; i<ainfo.n_aux; i++) {
	for (j=dlen-1; j>0; j--) {
	    de[i][j] = de[i][j-1];
	}
	de[i][0] = 0.0;
    }
}

private void zero_derivs(ArmaInfo * ainfo, double ** de, int dlen)
{
    int i;
    int j;

    for (i=0; i<ainfo.n_aux; i++) {
	for (j=0; j<dlen; j++) {
	    de[i][j] = 0.0;
	}
    }
}


private int arma_analytical_score(ArmaInfo * ainfo, const(double) * y, const(double) ** X, const(double) * phi, const(double) * Phi, const(double) * theta, const(double) * Theta, double s2, GretlMatrixRaw * G)
{
    /* forecast errors */
    const(double) * e = ainfo.e;
    /* pointers to blocks of derivatives (workspace) */
    double **de = ainfo.aux;
    double **de_a =    de + ainfo.ifc;
    double **de_sa = de_a + ainfo.np;
    double **de_m = de_sa + ainfo.P;
    double **de_sm = de_m + ainfo.nq;
    double **de_r = de_sm + ainfo.Q;
    int dlen = 1 + ainfo.q + ainfo.pd * ainfo.Q;
    double x;
    double Gsi;
    int t;
    int gt;
    int t1 = ainfo.t1;
    int i;
    int j;
    int k;
    int p;
    int s;
    int err = 0;

    zero_derivs(ainfo, de, dlen);


    for (t=t1, gt=0; t<=ainfo.t2 && !err; t++, gt++) {

	/* the constant term (de_0) */
	if (ainfo.ifc) {
	    de[0][0] = -1.0;
	    do_MA_partials(de[0], ainfo, theta, Theta, t);
	}

	/* non-seasonal AR terms (de_a) */
	k = 0;
	for (i=0; i<ainfo.p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    p = i + 1;
	    if (t - p >= 0) {
		de_a[k][0] = -y[t-p];
		/* cross-partial with seasonal AR */
		for (j=0; j<ainfo.P; j++) {
		    s = p + (j + 1) * ainfo.pd;
		    if (t - s >= 0) {
			de_a[k][0] += Phi[j] * y[t-s];
		    }
		}
		do_MA_partials(de_a[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal AR terms (de_sa) */
	for (j=0; j<ainfo.P; j++) {
	    p = (j + 1) * ainfo.pd;
	    if (t - p >= 0) {
		de_sa[j][0] = -y[t-p];
		/* cross-partial with non-seasonal AR */
		k = 0;
		for (i=0; i<ainfo.p; i++) {
		    if (AR_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sa[j][0] += phi[k] * y[t-s];
			}
			k++;
		    }
		}
		do_MA_partials(de_sa[j], ainfo, theta, Theta, t);
	    }
	}

	/* non-seasonal MA terms (de_m) */
	k = 0;
	for (i=0; i<ainfo.q; i++) {
	    if (!MA_included(ainfo, i)) {
		continue;
	    }
	    p = i + 1;
	    if (t - p >= 0) {
		de_m[k][0] = -e[t-p];
		/* cross-partial with seasonal MA */
		for (j=0; j<ainfo.Q; j++) {
		    s = p + (j + 1) * ainfo.pd;
		    if (t - s >= 0) {
			de_m[k][0] -= Theta[j] * e[t-s];
		    }
		}
		do_MA_partials(de_m[k], ainfo, theta, Theta, t);
	    }
	    k++;
	}

	/* seasonal MA terms (de_sm) */
	for (j=0; j<ainfo.Q; j++) {
	    p = (j + 1) * ainfo.pd;
	    if (t - p >= 0) {
		de_sm[j][0] = -e[t-p];
		/* cross-partial with non-seasonal MA */
		k = 0;
		for (i=0; i<ainfo.q; i++) {
		    if (MA_included(ainfo, i)) {
			s = p + i + 1;
			if (t - s >= 0) {
			    de_sm[j][0] -= theta[k] * e[t-s];
			}
			k++;
		    }
		}
		do_MA_partials(de_sm[j], ainfo, theta, Theta, t);
	    }
	}

	/* exogenous regressors (de_r) */
	for (j=0; j<ainfo.nexo; j++) {
	    de_r[j][0] = -X[j][t];
	    do_MA_partials(de_r[j], ainfo, theta, Theta, t);
	}

	/* update gradient matrix */
	x = e[t] / s2; /* sqrt(s2)? does it matter? */
	for (i=0; i<ainfo.nc; i++) {
	    Gsi = -de[i][0] * x;
	    if (na(Gsi)) {
		fprintf(stderr, "arma score, bad value at t=%d, i=%d\n", t, i);
		err = GretlError.nan;
		break;
	    }
	    gretl_matrix_set(G, gt, i, Gsi);
	}

	push_derivs(ainfo, de, dlen);
    }

    return err;
}

private int conditional_arma_forecast_errors(ArmaInfo * ainfo, const(double) * y, const(double) ** X, double b0, const(double) * phi, const(double) * Phi, const(double) * theta, const(double) * Theta, const(double) * beta, double * s2)
{
    double * e = ainfo.e;
    int t1 = ainfo.t1;
    int i;
    int j;
    int k;
    int s;
    int t;
    int p;

    *s2 = 0.0;

    for (t=t1; t<=ainfo.t2; t++) {
	e[t] = y[t];

	/* intercept */
	if (ainfo.ifc) {
	    e[t] -= b0;
	}

	/* non-seasonal AR component */
	k = 0;
	for (i=0; i<ainfo.p; i++) {
	    if (AR_included(ainfo, i)) {
		s = t - (i + 1);
		e[t] -= phi[k++] * y[s];
	    }
	}

	/* seasonal AR component plus interactions */
	for (j=0; j<ainfo.P; j++) {
	    s = t - (j + 1) * ainfo.pd;
	    e[t] -= Phi[j] * y[s];
	    k = 0;
	    for (i=0; i<ainfo.p; i++) {
		if (AR_included(ainfo, i)) {
		    p = s - (i + 1);
		    e[t] += Phi[j] * phi[k++] * y[p];
		}
	    }
	}

	/* non-seasonal MA component */
	k = 0;
	for (i=0; i<ainfo.q; i++) {
	    if (MA_included(ainfo, i)) {
		s = t - (i + 1);
		if (s >= ainfo.t1) {
		    e[t] -= theta[k] * e[s];
		}
		k++;
	    }
	}

	/* seasonal MA component plus interactions */
	for (j=0; j<ainfo.Q; j++) {
	    s = t - (j + 1) * ainfo.pd;
	    if (s >= ainfo.t1) {
		e[t] -= Theta[j] * e[s];
		k = 0;
		for (i=0; i<ainfo.q; i++) {
		    if (MA_included(ainfo, i)) {
			p = s - (i + 1);
			if (p >= ainfo.t1) {
			    e[t] -= Theta[j] * theta[k] * e[p];
			}
			k++;
		    }
		}
	    }
	}

	/* exogenous regressors */
	for (i=0; i<ainfo.nexo; i++) {
	    e[t] -= beta[i] * X[i][t];
	}

	*s2 += e[t] * e[t];
    }

    return 0;
}

/* Calculate ARMA log-likelihood.  This function is passed to the
   bhhh_max() routine as a callback. */

extern (C) double bhhh_arma_callback(double * coeff, GretlMatrixRaw * G, void * data, int do_score, int * err)
{
    ArmaInfo *ainfo = cast(ArmaInfo *) data;
    /* pointers to blocks of data */
    const(double) * y = ainfo.Z[0];
    const(double) ** X = ainfo.Z + 1;
    /* pointers to blocks of coefficients */
    double * phi = coeff + ainfo.ifc;
    double * Phi = phi + ainfo.np;
    double * theta = Phi + ainfo.P;
    double * Theta = theta + ainfo.nq;
    double * beta = Theta + ainfo.Q;
    double ll;
    double s2 = 0.0;
    int any_ma = ainfo.q > 0 || ainfo.Q > 0;

    *err = 0;


    if (any_ma && maybe_correct_MA(ainfo, theta, Theta)) {
	pputs(ainfo.prn, "arma: MA estimate(s) out of bounds\n");
	fputs("bhhh_arma_callback: MA estimate(s) out of bounds\n".ptr, stderr);
	*err = GretlError.noconv;
	return NADBL;
    }

    conditional_arma_forecast_errors(ainfo, y, X, coeff[0],
				     phi, Phi, theta, Theta,
				     beta, &s2);

    /* error variance and log-likelihood */
    s2 /= cast(double) ainfo.T;
    ll = -ainfo.T * (0.5 * log(s2) + LN_SQRT_2_PI_P5);

    if (isnan(ll)) {
	*err = GretlError.nan;
    }

    if (do_score) {
	ainfo.ll = ll;
	arma_analytical_score(ainfo, y, X,
			      phi, Phi, theta, Theta,
			      s2, G);
    }

    return ll;
}

/* construct a "virtual dataset" in the form of a set of pointers into
   the main dataset: this will be passed to the bhhh_max function.
   The dependent variable is put in position 0; following this are the
   independent variables.
*/

private const(double) ** make_arma_Z(ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    const(double) ** aZ;
    int * list = ainfo.alist;
    int ypos = arma_list_y_position(ainfo);
    int nx = list[0] - ypos;
    int v;
    int i;


    aZ = cast(const(double) **) malloc((nx + 1) * (*aZ).sizeof);
    if (aZ == null) {
	return null;
    }

    /* the dependent variable */
    if (ainfo.y != null) {
	aZ[0] = ainfo.y;
    } else {
	aZ[0] = dset.Z[list[ypos]];
    }

    /* the independent variables */
    for (i=1; i<=nx; i++) {
	v = list[i + ypos];
	aZ[i] = dset.Z[v];
    }

    return aZ;
}

/* add extra OPG-related stuff to the arma info struct */

private int set_up_arma_OPG_info(ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    /* array length needed for derivatives */
    int nd = 1 + ainfo.q + ainfo.pd * ainfo.Q;
    /* number of derivatives */
    int k = ainfo.nc;
    int err = 0;

    /* construct virtual dataset for dep var, real regressors */
    ainfo.Z = make_arma_Z(ainfo, dset);
    if (ainfo.Z == null) {
	err = GretlError.alloc;
    }

    if (!err) {
	/* allocate gradient matrix */
	ainfo.G = gretl_zero_matrix_new(ainfo.T, k);
	if (ainfo.G == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err && !arma_exact_ml(ainfo)) {
	/* allocate covariance matrix */
	ainfo.V = gretl_matrix_alloc(k, k);
	if (ainfo.V == null) {
	    err = GretlError.alloc;
	}
    }

    if (!err) {
	/* forecast errors array */
	ainfo.e = cast(double *) malloc((ainfo.t2 + 1) * (*ainfo.e).sizeof);
	if (ainfo.e == null) {
	    err = GretlError.alloc;
	} else {
	    int t;

	    for (t=0; t<=ainfo.t2; t++) {
		ainfo.e[t] = 0.0;
	    }
	}
    }

    if (!err) {
	/* derivatives arrays */
	ainfo.aux = doubles_array_new0(k, nd);
	if (ainfo.aux == null) {
	    err = GretlError.alloc;
	} else {
	    ainfo.n_aux = k;
	}
    }

    return err;
}

/* retrieve results specific to bhhh procedure */

private int conditional_arma_model_prep(Model * pmod, ArmaInfo * ainfo, double * theta)
{
    int i;
    int t;
    int err;

    pmod.t1 = ainfo.t1;
    pmod.t2 = ainfo.t2;
    pmod.nobs = pmod.t2 - pmod.t1 + 1;
    pmod.ncoeff = ainfo.nc;

    err = gretl_model_allocate_storage(pmod);
    if (err) {
	return err;
    }

    pmod.lnL = ainfo.ll;
    pmod.sigma = NADBL; /* will be replaced */

    for (i=0; i<pmod.ncoeff; i++) {
	pmod.coeff[i] = theta[i];
    }

    for (t=pmod.t1; t<=pmod.t2; t++) {
	pmod.uhat[t] = ainfo.e[t];
    }

    err = gretl_model_write_vcv(pmod, ainfo.V);

    return err;
}

extern (C) int bhhh_arma(double * theta, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod, Gretlopt opt)
{
    Gretlopt bhhh_opt = Gretlopt.OPT_NONE;
    double tol = libset_get_double(SetKey.bhhh_toler);
    int err = 0;

    err = set_up_arma_OPG_info(ainfo, dset);
    if (err) {
	pmod.errcode = err;
	return err;
    }

    if (opt & Gretlopt.OPT_V) {
	bhhh_opt |= Gretlopt.OPT_V;
    }

    err = bhhh_max(theta, ainfo.nc, ainfo.G,
		   &bhhh_arma_callback, tol,
		   &ainfo.fncount, &ainfo.grcount,
		   ainfo, ainfo.V, bhhh_opt, ainfo.prn);

    if (err) {
	fprintf(stderr, "arma: bhhh_max returned %d\n", err);
    } else {
	pmod.full_n = dset.n;
	err = conditional_arma_model_prep(pmod, ainfo, theta);
    }

    if (!err) {
	gretl_model_set_int(pmod, "fncount", ainfo.fncount);
	gretl_model_set_int(pmod, "grcount", ainfo.grcount);
	write_arma_model_stats(pmod, ainfo, dset);
	arma_model_add_roots(pmod, ainfo, theta);
    }

    if (err && !pmod.errcode) {
	pmod.errcode = err;
    }

    return pmod.errcode;
}
