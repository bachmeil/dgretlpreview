/* Port of plugin/arma_init.c: initialization of ARMA estimation (Hannan-Rissanen, and an AR model by OLS or NLS with minimal MA
   coefficients added), the OLS and least-squares estimators for pure AR models (arma_by_simple_ols, arma_by_ls), the scaling of y for
   exact ML (maybe_set_yscale), the transformation of the constant (transform_arma_const) and the reflection of roots inside the unit
   circle (flip_poly, used by maybe_correct_MA in gretl.plugin.arma).

   The entry points that arma_priv.h declares are extern (C) and public here; the debug output under AINIT_DEBUG (off in the C) is left out. */

module gretl.plugin.arma_init;

import gretl.dataset : create_auxiliary_dataset, dataset_add_series, destroy_dataset;
import gretl.describe : gretl_mean, gretl_moments;
import gretl.estimate : lsq;
import gretl.gretl_list : gretl_list_append_term, gretl_list_new;
import gretl.gretl_matrix : gretl_matrix_alloc, gretl_matrix_complex_divide, gretl_matrix_complex_multiply, gretl_matrix_free, gretl_matrix_multiply_by_scalar, gretl_matrix_polroots, gretl_matrix_shape, gretl_matrix_transpose_in_place, gretl_matrix_zero, gretl_zero_matrix_new;
import gretl.gretl_model : clear_model, gretl_model_init;
import gretl.gretl_prn : pprintf, pputc;
import gretl.gretl_utils : copyvec;
import gretl.libset : libset_get_double, libset_set_double;
import gretl.missingh : NADBL, na;
import gretl.nls : nlspec_new, aux_nlspec_add_param_list, model_from_nlspec, nlspec_destroy, nlspec_set_regression_function, nlspec_set_t1_t2;
import gretl.strutils : strings_array_free, strings_array_new_with_length;
import gretl.uservar : set_auxiliary_scalars, unset_auxiliary_scalars;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen, strcpy, strncpy, strcat, strchr;

import gretl.plugin.arma_priv;
import gretl.plugin.arma : arima_difference, arima_difference_undo;
import gretl.nlspech : Nlspec;
import gretl.libseth : SetKey;
import gretl.gretl_prn : Prn;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_commandsh : GretlCmdIndex;
import cephes.mconfh : CephesCmplx = Cmplx;
import gretl.matrix_extrah : MMissingCode;
import gretl.gretl_modelh : ArmaFlags;





/* Given an estimate of the ARMA constant via OLS, convert to the form
   wanted for initializing the Kalman filter.  Note: the @b array
   goes: const, phi, Phi, theta, Theta, beta.
*/

extern (C) void transform_arma_const(double * b, ArmaInfo * ainfo)
{
    const(double) * phi = b + 1;
    const(double) * Phi = phi + ainfo.np;
    double narfac = 1.0; /* nonseasonal AR factor */
    double sarfac = 1.0; /* seasonal AR factor */
    int i;
    int k = 0;

    if (ainfo.np == 0 && ainfo.P == 0) {
	/* nothing to be done */
	return;
    }


    for (i=0; i<ainfo.p; i++) {
	if (AR_included(ainfo, i)) {
	    narfac -= phi[k++];
	}
    }

    for (i=0; i<ainfo.P; i++) {
	sarfac -= Phi[i];
    }

    b[0] /= (narfac * sarfac);

}

private int init_transform_const(ArmaInfo * ainfo)
{
    return ainfo.ifc && arma_exact_ml(ainfo);
}

extern (C) void maybe_set_yscale(ArmaInfo * ainfo)
{
    double ybar;
    double sdy;
    int err;

    if (arima_levels(ainfo)) {
	/* not sure about this clause! */
	ybar = gretl_mean(ainfo.t1, ainfo.t2, ainfo.y);
	if (fabs(ybar) > 250) {
	    set_arma_avg_ll(ainfo);
	}
	return;
    }

    err = gretl_moments(ainfo.t1, ainfo.t2, ainfo.y,
			null, &ybar, &sdy, null, null, 1);

    if (!err && sdy > 0) {
	/* try to catch cases where (a) the overall scale is
	   too big or (b) the coefficient of variation is
	   too small: in such cases set up conversion to
	   (y - ybar) / sdy + 1 = (y - (ybar - sdy)) * 1/sdy.
	*/
	double abs_ybar = fabs(ybar);
	double hi = 200;
	double lo = 0.01;
	// double hi = 1.25, lo = 0.75;

	if (abs_ybar > hi || abs_ybar < lo || sdy/abs_ybar < lo) {
	    ainfo.yshift = ybar - sdy; /* subtract */
	    ainfo.yscale = 1 / sdy;    /* multiply */
	}
    }

    if (!err && ainfo.prn != null && ainfo.yscale != 1.0) {
	pputc(ainfo.prn, '\n');
	pprintf(ainfo.prn, _("Shifting y by %g, scaling by %g\n").ptr,
		ainfo.yshift, ainfo.yscale);
    }
}

bool apply_yscaling(const(ArmaInfo) * a, double x)
{
    return arma_exact_ml(a) && !na(x);
}

enum HR_MINLAGS = 16;

/* complex inversion of @z */

private GretlMatrixRaw * cinv(GretlMatrixRaw * z)
{
    GretlMatrixRaw * tmp;
    GretlMatrixRaw * ret = null;
    int n = z.rows;
    int i;
    int err = 0;

    tmp = gretl_zero_matrix_new(n, 2);
    for (i=0; i<n; i++) {
	tmp.val[i] = 1.0;
    }
    ret = gretl_matrix_complex_divide(tmp, z, 1, &err);
    gretl_matrix_free(tmp);

    return ret;
}

private void copy_row(GretlMatrixRaw * targ, int it, const(GretlMatrixRaw) * src, int is_, int neg)
{
    double d;
    int j;

    for (j=0; j<src.cols; j++) {
	d = gretl_matrix_get(src, is_, j);
	gretl_matrix_set(targ, it, j, neg ? -d : d);
    }
}

/* computes a polynomial from its roots: @r is assumed
   to be n x 2 (complex)
*/

private GretlMatrixRaw * pol_from_roots(const(GretlMatrixRaw) * r)
{
    GretlMatrixRaw * tmp;
    GretlMatrixRaw * ret = null;
    int n = r.rows;
    int err = 0;

    tmp = gretl_matrix_alloc(1, 2);

    if (n == 0) {
	tmp.val[0] = 1;
	tmp.val[1] = 0;
	ret = tmp;
    } else {
	copy_row(tmp, 0, r, n-1, 0);
	if (tmp.val[0] == 0 && tmp.val[1] == 0) {
	    tmp.val[0] = tmp.val[1] = NADBL;
	    ret = tmp;
        } else {
	    GretlMatrixRaw * ix = cinv(tmp);
	    int i;

            if (n == 1) {
		/* hansl: ret = {1,0} | -ix */
		ret = gretl_zero_matrix_new(ix.rows + 1, 2);
		ret.val[0] = 1;
		for (i=0; i<ix.rows; i++) {
		    copy_row(ret, i+1, ix, i, 1);
		}
            } else {
		GretlMatrixRaw * rslice; /* hansl: = r[1:n-1,] */
		GretlMatrixRaw * tmp1;
		GretlMatrixRaw * tmp2;
		double d0;
		double d1;

		rslice = gretl_matrix_alloc(n-1, 2);
		for (i=0; i<rslice.rows; i++) {
		    copy_row(rslice, i, r, i, 0);
		}
		gretl_matrix_free(tmp);
                tmp = pol_from_roots(rslice);
		/* hansl: ret = tmp | {0,0} */
		ret = gretl_zero_matrix_new(tmp.rows + 1, 2);
		for (i=0; i<tmp.rows; i++) {
		    copy_row(ret, i, tmp, i, 0);
		}
		/* hansl: ix = transp(mshape(ix, 2, n)) */
		tmp1 = gretl_matrix_shape(ix, 2, n, &err);
		gretl_matrix_transpose_in_place(tmp1);
		/* hansl: tmp = cmult(tmp1, -ix) */
		gretl_matrix_multiply_by_scalar(tmp1, -1.0);
		tmp2 = gretl_matrix_complex_multiply(tmp, tmp1, 1, &err);
		/* hansl: ret[2:,] += tmp */
		for (i=1; i<ret.rows; i++) {
		    d0 = gretl_matrix_get(ret, i, 0);
		    d0 += gretl_matrix_get(tmp2, i-1, 0);
		    d1 = gretl_matrix_get(ret, i, 1);
		    d1 += gretl_matrix_get(tmp2, i-1, 1);
		    gretl_matrix_set(ret, i, 0, d0);
		    gretl_matrix_set(ret, i, 1, d1);
		}
		gretl_matrix_free(tmp1);
		gretl_matrix_free(tmp2);
		gretl_matrix_free(rslice);
	    }
	    gretl_matrix_free(ix);
	}
    }

    if (ret == tmp) {
	tmp = null;
    }
    gretl_matrix_free(tmp);

    return ret;
}

/* for non-seasonal "gappy" coeff vector: expand
   by inserting zeros as needed */

private GretlMatrixRaw * poly_from_coeff(const(double) * coeff, const(char) * mask, int n, int ar)
{
    GretlMatrixRaw * ret;
    int i;
    int k = 0;

    ret = gretl_zero_matrix_new(n + 1, 1);
    ret.val[0] = 1.0;

    for (i=0; i<n; i++) {
        if (mask[i] == '1') {
	    ret.val[i+1] = ar ? -coeff[k] : coeff[k];
	    k++;
        }
    }

    return ret;
}

/* checks if the polynomial given by @coeff is fundamental,
   and modifies it if that is not the case.
*/

extern (C) int flip_poly(double * coeff, ArmaInfo * ainfo, int ar, int seasonal)
{
    GretlMatrixRaw * tmp = null;
    GretlMatrixRaw * r = null;
    const(char) * mask;
    double re;
    double im;
    int n;
    int n_inside = 0;
    int i;
    int err = 0;

    if (ar) {
	n = seasonal ? ainfo.P : ainfo.p;
	mask = seasonal ? null : ainfo.pmask;
    } else {
	n = seasonal ? ainfo.Q : ainfo.q;
	mask = seasonal ? null : ainfo.qmask;
    }

    if (mask == null) {
	/* no expansion needed */
	tmp = gretl_matrix_alloc(n + 1, 1);
	tmp.val[0] = 1.0;
	for (i=1; i<=n; i++) {
	    tmp.val[i] = ar ? -coeff[i-1] : coeff[i-1];
	}
    } else {
	/* expand to handle gappiness */
	tmp = poly_from_coeff(coeff, mask, n, ar);
    }

    /* force legacy form of complex output, for now */
    r = gretl_matrix_polroots(tmp, 1, 1, &err);

    if (err) {
	goto bailout;
    }

    gretl_matrix_zero(tmp);
    for (i=0; i<r.rows; i++) {
	re = gretl_matrix_get(r, i, 0);
	im = gretl_matrix_get(r, i, 1);
	if (re*re + im*im < 1.0) {
	    /* record row reference */
	    tmp.val[i] = 1;
	    n_inside++;
	}
    }

    if (n_inside > 0) {
	GretlMatrixRaw * rfix;
	GretlMatrixRaw * ifix;
	double ci1;
	int k = 0;

	/* compose sub-matrix */
	rfix = gretl_matrix_alloc(n_inside, 2);
	for (i=0; i<r.rows; i++) {
	    if (tmp.val[i] == 1) {
		copy_row(rfix, k++, r, i, 0);
	    }
	}
	/* complex inversion */
	ifix = cinv(rfix);
	if (ifix == null) {
	    err = GretlError.data;
	    goto bailout;
	}
	/* replace the inverted portion of r */
	k = 0;
	for (i=0; i<r.rows; i++) {
	    if (tmp.val[i] == 1) {
		copy_row(r, i, ifix, k++, 0);
	    }
	}
	gretl_matrix_free(tmp);
        tmp = pol_from_roots(r);
	if (tmp == null) {
	    err = GretlError.data;
	    goto bailout;
	}
	if (mask != null) {
	    /* shrink to coeff */
	    k = 0;
	    for (i=0; i<n; i++) {
		if (mask[i] == '1') {
		    ci1 = tmp.val[i+1];
		    coeff[k++] = ar ? -ci1 : ci1;
		}
	    }
	} else {
	    /* just copy to coeff */
	    for (i=0; i<n; i++) {
		ci1 = tmp.val[i+1];
		coeff[i] = ar ? -ci1 : ci1;
	    }
	}
	gretl_matrix_free(rfix);
	gretl_matrix_free(ifix);
    }

 bailout:

    gretl_matrix_free(r);
    gretl_matrix_free(tmp);

    return err;
}

/* @pmod.coeff contains coefficients from step 2 of
   the H-R procedure, in the order: intercept, exogenous
   vars, phi, Phi, theta, Theta (in each case, if present).
   The array @b has to be filled in the order: intercept,
   phi, Phi, theta, Theta, exogenous vars.
*/

private int hr_transcribe_coeffs(ArmaInfo * ainfo, Model * pmod, double * b)
{
    double * phi = null;
    double * Phi = null;
    double * theta = null;
    double * Theta = null;
    int j = ainfo.nexo + ainfo.ifc;
    int i;
    int k = 0;
    int err = 0;

    if (ainfo.ifc) {
	b[0] = pmod.coeff[0];
	if (arma_xdiff(ainfo)) {
	    b[0] /= ainfo.T;
	}
	k = 1;
    }

    phi = b + k;
    for (i=0; i<ainfo.p; i++) {
	if (AR_included(ainfo, i)) {
	    b[k++] = pmod.coeff[j++];
	}
    }

    Phi = b + k;
    for (i=0; i<ainfo.P; i++) {
	b[k++] = pmod.coeff[j];
	j += ainfo.np + 1; /* assumes ainfo.p < pd */
    }

    theta = b + k;
    for (i=0; i<ainfo.q; i++) {
	if (MA_included(ainfo, i)) {
	    b[k++] = pmod.coeff[j++];
	}
    }

    Theta = b + k;
    for (i=0; i<ainfo.Q; i++) {
	b[k++] = pmod.coeff[j];
	j += ainfo.nq + 1; /* assumes ainfo.q < pd */
    }

    j = ainfo.ifc;
    for (i=0; i<ainfo.nexo; i++) {
	b[k++] = pmod.coeff[j++];
    }

    if (ainfo.p > 0) {
	flip_poly(phi, ainfo, 1, 0);
    }
    if (ainfo.P > 0) {
	flip_poly(Phi, ainfo, 1, 1);
    }
    if (ainfo.q > 0) {
	flip_poly(theta, ainfo, 0, 0);
    }
    if (ainfo.Q > 0) {
	flip_poly(Theta, ainfo, 0, 1);
    }

    return err;
}

private int pre_sample_count(ArmaInfo * ainfo, const(double) * y, int maxlag)
{
    int t;
    int n = 0;

    for (t=ainfo.t1-1; t>=0; t--) {
	if (!na(y[t])) {
	    n++;
	    if (n == maxlag) {
		break;
	    }
	} else {
	    break;
	}
    }

    return n;
}

private double * prescale_y(double * y, ArmaInfo * ainfo, int n)
{
    double *ys = copyvec(y, n);
    int t;

    for (t=0; t<n; t++) {
	if (!na(ys[t])) {
	    ys[t] -= ainfo.yshift;
	    ys[t] *= ainfo.yscale;
	}
    }

    return ys;
}

/* Hannan-Rissanen ARMA initialization via two OLS passes. In the
   first pass we run an OLS regression of y on the exogenous vars plus
   a certain (biggish) number of lags. In the second we estimate the
   ARMA model by OLS, substituting innovations and corresponding lags
   with the first-pass residuals.
*/

private int real_hr_arma_init(double * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Prn * prn)
{
    double * y;
    double * dest;
    double * src;
    DatasetRaw * aset = null;
    Model armod;
    int * hrlist = null;
    char * done = null;
    int maxp2p;
    int maxp2q;
    int p1lags;
    int nv;
    int nv1;
    int nv2;
    int np2p;
    int np2q;
    int m;
    int xpos;
    int pos;
    int mapos;
    size_t datalen;
    int free_y = 0;
    int i;
    int j;
    int T;
    int t1;
    int err = 0;

    /* the dependent variable (of length dset.n) */
    if (arma_xdiff(ainfo)) {
	/* for initialization, use the level of y */
	y = cast(double *) dset.Z[ainfo.yno];
    } else {
	y = ainfo.y;
    }

    /* the greatest AR lag-length in pass 2 */
    maxp2p = ainfo.p + ainfo.pd * ainfo.P;
    /* the greatest MA lag-length in pass 2 */
    maxp2q = ainfo.q + ainfo.pd * ainfo.Q;

    /* the actual number of lags in pass 2 */
    np2p = ainfo.np + ainfo.P + ainfo.p * ainfo.P;
    np2q = ainfo.nq + ainfo.Q + ainfo.q * ainfo.Q;
    nv2 = 2 + ainfo.nexo + np2p + np2q;

    /* do we need more AR lags for H-R? */
    p1lags = maxp2p < HR_MINLAGS ? HR_MINLAGS : maxp2p;
    nv1 = 2 + ainfo.nexo + p1lags;

    /* how many variables do we need to allocate for? */
    nv = nv1 > nv2 ? nv1 : nv2;

    /* and how many observations? */
    T = ainfo.T - p1lags;
    if (ainfo.t1 > 0) {
	/* use non-missing pre-sample obs? */
	T += pre_sample_count(ainfo, y, p1lags);
    }
    /* benchmark position for reading from dset */
    t1 = ainfo.t1 + ainfo.T - T;

    /* allocate sufficient storage for both passes */
    aset = create_auxiliary_dataset(nv, T, Gretlopt.OPT_NONE);
    if (aset == null) {
	return GretlError.alloc;
    }


    if (ainfo.yscale != 1.0) {
	y = prescale_y(y, ainfo, dset.n);
	if (y == null) {
	    err = GretlError.alloc;
	    goto bailout;
	} else {
	    free_y = 1;
	}
    }

    /* in case we fail before estimating a model */
    gretl_model_init(&armod, dset);

    /* regression list */
    hrlist = gretl_list_new(nv);
    if (hrlist == null) {
	err = GretlError.alloc;
	goto bailout;
    } else {
	hrlist[1] = 1;
	hrlist[2] = 0;
	for (i=2; i<nv; i++) {
	    hrlist[i+1] = i;
	}
    }

    /* adjust the list for pass 1 */
    hrlist[0] = nv1;

    /* recorder array, etc. */
    done = cast(char *) calloc(p1lags, 1);
    datalen = T * double.sizeof;

    /* dependent var */
    strcpy(aset.varname[1], "y");
    memcpy(aset.Z[1], y + t1, datalen);
    /* exogenous vars */
    pos = 2;
    for (i=0; i<ainfo.nexo; i++) {
	sprintf(aset.varname[pos], "x%d", i);
	xpos = ainfo.xlist[i+1];
	memcpy(aset.Z[pos], dset.Z[xpos] + t1, datalen);
	pos++;
    }
    /* pass2 non-seasonal AR lags first */
    for (i=1; i<=ainfo.p; i++) {
	if (AR_included(ainfo, i-1)) {
	    sprintf(aset.varname[pos], "y_%d", i);
	    memcpy(aset.Z[pos], y + t1 - i, datalen);
	    done[i-1] = 1;
	    pos++;
	}
    }
    /* then pass2 seasonal AR lags */
    for (j=1; j<=ainfo.P; j++) {
	m = j * ainfo.pd;
	sprintf(aset.varname[pos], "y_%d", m);
	memcpy(aset.Z[pos], y + t1 - m, datalen);
	done[m-1] = 1;
	pos++;
    }
    /* then pass2 AR interactions */
    for (j=1; j<=ainfo.P; j++) {
	for (i=1; i<=ainfo.p; i++) {
	    if (AR_included(ainfo, i-1)) {
		m = j * ainfo.pd + i;
		sprintf(aset.varname[pos], "y_%d", m);
		memcpy(aset.Z[pos], y + t1 - m, datalen);
		done[m-1] = 1;
		pos++;
	    }
	}
    }
    mapos = pos; /* insertion point for pass2 MA lags */
    /* then any "extra" AR lags for pass 1 only */
    for (i=1; i<=p1lags; i++) {
	if (!done[i-1]) {
	    sprintf(aset.varname[pos], "y_%d", i);
	    memcpy(aset.Z[pos], y + t1 - i, datalen);
	    pos++;
	}
    }

    /* pass 1 estimation (FIXME constant?) */
    armod = lsq(hrlist, aset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    if (armod.errcode) {
	err = armod.errcode;
	goto bailout;
    }


    /* revise hrlist if needed (FIXME constant?) */
    hrlist[0] = nv2;
    /* position for insertion of MA terms: at least some of
       these will likely overwrite AR terms that were wanted
       only for pass 1
    */
    pos = mapos;

    /* pass 2 sample: skip the leading observations for
       which we don't have lagged residuals from pass 1
    */
    aset.t1 = maxp2q;
    datalen = (T - maxp2q) * double.sizeof;

    for (i=1; i<=ainfo.q; i++) {
	if (MA_included(ainfo, i-1)) {
	    sprintf(aset.varname[pos], "e_%d", i);
	    dest = aset.Z[pos] + aset.t1;
	    src = armod.uhat + aset.t1 - i;
	    memcpy(dest, src, datalen);
	    pos++;
	}
    }
    for (j=1; j<=ainfo.Q; j++) {
	m = j * ainfo.pd;
	sprintf(aset.varname[pos], "e_%d", m);
	dest = aset.Z[pos] + aset.t1;
	src = armod.uhat + aset.t1 - m;
	memcpy(dest, src, datalen);
	pos++;
    }
    for (j=1; j<=ainfo.Q; j++) {
	for (i=1; i<=ainfo.q; i++) {
	    if (MA_included(ainfo, i-1)) {
		m = j * ainfo.pd + i;
		sprintf(aset.varname[pos], "e_%d", m);
		dest = aset.Z[pos] + aset.t1;
		src = armod.uhat + aset.t1 - m;
		memcpy(dest, src, datalen);
		pos++;
	    }
	}
    }

    /* pass 2 estimation */
    clear_model(&armod);
    armod = lsq(hrlist, aset, GretlCmdIndex.ols, Gretlopt.OPT_A);

    if (armod.errcode) {
	err = armod.errcode;
    } else {
	err = hr_transcribe_coeffs(ainfo, &armod, coeff);
	if (!err && ainfo.nexo == 0 && init_transform_const(ainfo)) {
	    transform_arma_const(coeff, ainfo);
	}
    }


 bailout:

    free(hrlist);
    free(done);
    destroy_dataset(aset);
    clear_model(&armod);
    if (free_y) {
	free(y);
    }

    if (!err) {
	ainfo.init_ = IniMethod.hr;
    }

    return err;
}

/* Do we have enough observations to do Hannan-Rissanen? */

private int hr_df_check(ArmaInfo * ainfo, const(DatasetRaw) * dset)
{
    int nobs = ainfo.T;
    int nlags = (ainfo.P + ainfo.Q) * dset.pd;
    int ncoeff;
    int df;
    int ok = 1;

    if (nlags < HR_MINLAGS) {
	nlags = HR_MINLAGS;
    }

    ncoeff = nlags + ainfo.nexo + ainfo.ifc;
    nobs -= nlags;
    df = nobs - ncoeff;

    if (df < 1) {
	ok = 0;
    }


    return ok;
}

extern (C) int hr_arma_init(double * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo)
{
    int ok = hr_df_check(ainfo, dset);
    int err = 0;

    if (ok) {
	err = real_hr_arma_init(coeff, dset, ainfo, ainfo.prn);
    }


    return err;
}

private double get_y_mean(ArmaInfo * ainfo)
{
    double ysum = 0.0;
    int t;
    int T = 0;

    for (t=ainfo.t1; t<=ainfo.t2; t++) {
	if (!na(ainfo.y[t])) {
	    if (ainfo.yscale != 1.0) {
		ysum += (ainfo.y[t] - ainfo.yshift) * ainfo.yscale;
	    } else {
		ysum += ainfo.y[t];
	    }
	    T++;
	}
    }

    return ysum / T;
}

enum MA_SMALL = 0.0001;

/* Transcribe coeffs from the OLS or NLS model used for initializing,
   into the array @b that will be passed to the maximizer. While
   we're at it, check that the AR polynomial(s) are in the
   stationary zone.
*/

private void ar_init_transcribe_coeffs(ArmaInfo * ainfo, Model * pmod, double * b)
{
    int q0 = ainfo.ifc + ainfo.np + ainfo.P;
    int totq = ainfo.nq + ainfo.Q;
    int i;
    int j = 0;

    for (i=0; i<pmod.ncoeff; i++) {
	if (i == q0 && totq > 0) {
	    /* reserve space for MA terms */
	    j += totq;
	}
	if (j < ainfo.nc) {
	    b[j++] = pmod.coeff[i];
	}
    }

    if (arma_xdiff(ainfo) && ainfo.ifc) {
	/* is this a good idea? */
	b[0] /= ainfo.T;
    }

    /* insert near-zeros for MA terms */
    for (i=0; i<totq; i++) {
	b[q0 + i] = MA_SMALL;
    }

    /* stationarity check/fix */
    if (ainfo.p > 0) {
	flip_poly(b + ainfo.ifc, ainfo, 1, 0);
    }
    if (ainfo.P > 0) {
	flip_poly(b + ainfo.ifc + ainfo.np, ainfo, 1, 1);
    }
}

/* compose variable names for temporary dataset */

private void arma_init_add_varnames(ArmaInfo * ainfo, int ptotal, int narmax, DatasetRaw * aset)
{
    int i;
    int j;
    int k;
    int kx;
    int ky;
    int lag;
    int k0 = 2;

    strcpy(aset.varname[1], "y");

    k = k0;
    kx = ptotal + ainfo.nexo + k0;

    for (i=0; i<ainfo.p; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(aset.varname[k++], "y_%d", lag);
	    for (j=0; j<narmax; j++) {
		sprintf(aset.varname[kx++], "x%d_%d", j+1, lag);
	    }
	}
    }

    ky = ainfo.np + ainfo.P + k0;

    for (j=0; j<ainfo.P; j++) {
	lag = (j + 1) * ainfo.pd;
	k = k0 + ainfo.np + j;
	sprintf(aset.varname[k], "y_%d", lag);
	for (i=0; i<narmax; i++) {
	    sprintf(aset.varname[kx++], "x%d_%d", i+1, lag);
	}
	for (i=0; i<ainfo.p; i++) {
	    if (AR_included(ainfo, i)) {
		lag = (j + 1) * ainfo.pd + (i + 1);
		sprintf(aset.varname[ky++], "y_%d", lag);
		for (k=0; k<narmax; k++) {
		    sprintf(aset.varname[kx++], "x%d_%d", k+1, lag);
		}
	    }
	}
    }

    kx = ptotal + k0;

    for (i=0; i<ainfo.nexo; i++) {
	sprintf(aset.varname[kx++], "x%d", i+1);
    }
}

/* experimental: when initializing an AR(I)MA model via
   NLS, work around interior NAs by adding observation-
   specific dummies to the dataset
*/

private int arma_init_add_dummies(ArmaInfo * ainfo, DatasetRaw * dset)
{
    int * misslist = null;
    int t1 = dset.t1;
    int i;
    int t;
    int err = 0;

    /* if we have a block of leading NAs, skip it */

    for (t=t1; t<=dset.t2 && !err; t++) {
	int miss = 0;

	for (i=1; i<dset.v; i++) {
	    if (na(dset.Z[i][t])) {
		miss = 1;
		break;
	    }
	}
	if (miss) {
	    t1++;
	} else {
	    break;
	}
    }

    /* form list of observation indices of interior NAs */

    for (t=t1; t<=dset.t2 && !err; t++) {
	for (i=1; i<dset.v; i++) {
	    if (na(dset.Z[i][t])) {
		misslist = gretl_list_append_term(&misslist, t);
		if (misslist == null) {
		    err = GretlError.alloc;
		}
		break;
	    }
	}
    }


    if (misslist != null) {
	/* For each observation with any missing values, add
	   a specific dummy and zero out the missing data.
	*/
	int origv = dset.v;
	int j;
	int v;
	int nd = misslist[0];

	err = dataset_add_series(dset, nd);
	if (!err) {
	    for (i=1; i<=misslist[0]; i++) {
		v = origv + i - 1;
		t = misslist[i];
		sprintf(dset.varname[v], "d%d", i);
		dset.Z[v][t] = 1.0;
		for (j=1; j<origv; j++) {
		    if (na(dset.Z[j][t])) {
			dset.Z[j][t] = 0.0;
		    }
		}
	    }
	}
    }

    ainfo.misslist = misslist;

    return err;
}

/* X, if non-null, holds the differenced regressors */

private double get_xti(const(DatasetRaw) * dset, int i, int t, const(int) * xlist, const(GretlMatrixRaw) * X)
{
    if (X != null) {
	return gretl_matrix_get(X, t, i);
    } else {
        return dset.Z[xlist[i+1]][t];
    }
}

/* Build temporary dataset including lagged vars: if we're doing exact
   ML on an ARMAX model we need lags of the exogenous variables as
   well as lags of y_t.  Note that the auxiliary dataset has "t = 0"
   at an offset of ainfo.t1 into the "real", external dataset.
*/

private int arma_init_build_dataset(ArmaInfo * ainfo, int ptotal, int narmax, const(int) * list, const(DatasetRaw) * dset, DatasetRaw * aset, int nonlin)
{
    double **aZ = aset.Z;
    const(double) * y;
    const(GretlMatrixRaw) * X = null;
    const(int) * xlist = ainfo.xlist;
    int i;
    int j;
    int k;
    int kx;
    int ky;
    int t;
    int s;
    int k0 = 2;
    int undo_diff = 0;
    int err = 0;

    if (arima_levels(ainfo)) {
	/* we'll need differences for initialization */
	err = arima_difference(ainfo, dset, 1);
	if (err) {
	    return err;
	}
	undo_diff = 1;
	y = ainfo.y;
	X = ainfo.dX;
    } else if (arma_xdiff(ainfo)) {
	/* run init in levels (FIXME?) */
	y = dset.Z[ainfo.yno];
    } else {
	y = ainfo.y;
    }

    /* add variable names to auxiliary dataset */
    arma_init_add_varnames(ainfo, ptotal, narmax, aset);

    for (t=0; t<aset.n; t++) {
	int realt = t + ainfo.t1;
	int miss = 0;

	if (apply_yscaling(ainfo, y[realt])) {
	    aZ[1][t] = (y[realt] - ainfo.yshift) * ainfo.yscale;
	} else {
	    aZ[1][t] = y[realt];
	}

	k = k0;
	kx = ptotal + ainfo.nexo + k0;

	for (i=0; i<ainfo.p; i++) {
	    if (!AR_included(ainfo, i)) {
		continue;
	    }
	    s = realt - (i + 1);
	    if (s < 0) {
		miss = 1;
		aZ[k++][t] = NADBL;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		if (apply_yscaling(ainfo, y[s])) {
		    aZ[k][t] -= ainfo.yshift;
		    aZ[k][t] *= ainfo.yscale;
		}
		k++;
		for (j=0; j<narmax; j++) {
		    aZ[kx++][t] = get_xti(dset, j, s, xlist, X);
		}
	    }
	}

	ky = ainfo.np + ainfo.P + k0;

	for (j=0; j<ainfo.P; j++) {
	    s = realt - (j + 1) * ainfo.pd;
	    k = ainfo.np + k0 + j;
	    if (s < 0) {
		miss = 1;
		aZ[k][t] = NADBL;
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = NADBL;
		}
	    } else {
		aZ[k][t] = y[s];
		if (apply_yscaling(ainfo, y[s])) {
		    aZ[k][t] -= ainfo.yshift;
		    aZ[k][t] *= ainfo.yscale;
		}
		for (k=0; k<narmax; k++) {
		    aZ[kx++][t] = get_xti(dset, k, s, xlist, X);
		}
	    }
	    for (i=0; i<ainfo.p; i++) {
		if (!AR_included(ainfo, i)) {
		    continue;
		}
		s = realt - ((j + 1) * ainfo.pd + (i + 1));
		if (s < 0) {
		    miss = 1;
		    aZ[ky++][t] = NADBL;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = NADBL;
		    }
		} else {
		    aZ[ky][t] = y[s];
		    if (apply_yscaling(ainfo, y[s])) {
			aZ[ky][t] -= ainfo.yshift;
			aZ[ky][t] *= ainfo.yscale;
		    }
		    ky++;
		    for (k=0; k<narmax; k++) {
			aZ[kx++][t] = get_xti(dset, k, s, xlist, X);
		    }
		}
	    }
	}

	kx = ptotal + k0;

	for (i=0; i<ainfo.nexo; i++) {
	    aZ[kx++][t] = get_xti(dset, i, realt, xlist, X);
	}

	if (miss) {
	    aset.t1 = t + 1;
	}
    }

    if (nonlin && arma_missvals(ainfo)) {
	err = arma_init_add_dummies(ainfo, aset);
    }

    if (undo_diff) {
	arima_difference_undo(ainfo, dset);
    }


    return err;
}

private void nls_kickstart(Model * pmod, DatasetRaw * dset, double * b0, double * by1)
{
    int[4] list;

    if (b0 != null) {
	list[0] = 3;
	list[1] = 1;
	list[2] = 0;
	list[3] = 2;
    } else {
	list[0] = 2;
	list[1] = 1;
	list[2] = 2;
    }

    *pmod = lsq(list.ptr, dset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_Z);

    if (!pmod.errcode) {
	if (b0 != null) {
	    *b0 = pmod.coeff[0];
	    *by1 = pmod.coeff[1];
	} else {
	    *by1 = pmod.coeff[0];
	}
	if (*by1 >= 1.0) {
	   *by1 = 0.95;
	}
    }

    clear_model(pmod);
}

private int add_to_spec(char * targ, const(char) * src)
{
    if (strlen(src) + strlen(targ) > MAXLINE - 1) {
	return 1;
    } else {
	strcat(targ, src);
	return 0;
    }
}

/* for ARMAX: write the component of the NLS specification
   that takes the form (y_{t-i} - X_{t-i} \beta)
*/

private int y_Xb_at_lag(char * spec, ArmaInfo * ainfo, int narmax, int lag)
{
    char[32] chunk;
    int i;
    int nt;
    int err = 0;

    if (narmax == 0) {
	sprintf(chunk.ptr, "y_%d".ptr, lag);
	return add_to_spec(spec, chunk.ptr);
    }

    nt = ainfo.ifc + narmax;

    sprintf(chunk.ptr, "(y_%d-".ptr, lag);

    if (nt > 1) {
	strcat(chunk.ptr, "(".ptr);
    }

    if (ainfo.ifc) {
	strcat(chunk.ptr, "_b0".ptr);
    }

    err = add_to_spec(spec, chunk.ptr);

    for (i=0; i<narmax && !err; i++) {
	if (ainfo.ifc || i > 0) {
	    err += add_to_spec(spec, "+");
	}
	sprintf(chunk.ptr, "_b%d*x%d_%d".ptr, i+1, i+1, lag);
	err += add_to_spec(spec, chunk.ptr);
    }

    if (nt > 1) {
	err += add_to_spec(spec, "))");
    } else {
	err += add_to_spec(spec, ")");
    }

    return err;
}

private int arma_get_nls_model(Model * amod, ArmaInfo * ainfo, int narmax, const(double) * coeff, DatasetRaw * dset, Prn * prn)
{
    Gretlopt nlsopt = Gretlopt.OPT_A;
    char[MAXLINE] fnstr;
    char[32] term;
    Nlspec * spec;
    double * parms = null;
    char **pnames = null;
    double * b0 = null;
    double * by1 = null;
    int nparam;
    int lag;
    int i;
    int j;
    int k;
    int err = 0;

    spec = nlspec_new(GretlCmdIndex.nls, dset);
    if (spec == null) {
	return GretlError.alloc;
    }

    if (arma_least_squares(ainfo)) {
	/* respect verbose option */
	if (prn != null) {
	    nlsopt |= Gretlopt.OPT_V;
	}
    } else {
	/* don't bother with standard errors */
	nlsopt |= Gretlopt.OPT_C;
    }

    nlspec_set_t1_t2(spec, 0, ainfo.T - 1);

    nparam = ainfo.ifc + ainfo.np + ainfo.P + ainfo.nexo;

    if (ainfo.misslist != null) {
	nparam += ainfo.misslist[0];
    }

    parms = cast(double *) malloc(nparam * (*parms).sizeof);
    if (parms == null) {
	err = GretlError.alloc;
	goto bailout;
    }

    pnames = strings_array_new_with_length(nparam, VNAMELEN);
    if (pnames == null) {
	err = GretlError.alloc;
	goto bailout;
    }

    /* make names for the parameters; construct the param list;
       and do some rudimentary fall-back initialization */

    for (i=0; i<nparam; i++) {
	parms[i] = 0.0;
    }

    k = 0;

    if (ainfo.ifc) {
	if (coeff != null) {
	    parms[k] = coeff[k];
	} else {
	    parms[k] = gretl_mean(0, dset.n - 1, dset.Z[1]);
	}
	b0 = &parms[k];
	strcpy(pnames[k++], "_b0");
    }

    for (i=0; i<ainfo.p; i++) {
	if (AR_included(ainfo, i)) {
	    if (by1 == null) {
		by1 = &parms[k];
		if (coeff == null) {
		    parms[k] = 0.1;
		}
	    }
	    if (coeff != null) {
		parms[k] = coeff[k];
	    }
	    sprintf(pnames[k++], "_phi%d", i+1);
	}
    }

    for (i=0; i<ainfo.P; i++) {
	if (by1 == null) {
	    by1 = &parms[k];
	    if (coeff == null) {
		parms[k] = 0.1;
	    }
	}
	if (coeff != null) {
	    parms[k] = coeff[k];
	}
	sprintf(pnames[k++], "_Phi%d", i+1);
    }

    for (i=0; i<ainfo.nexo; i++) {
	if (coeff != null) {
	    parms[k] = coeff[k];
	}
	sprintf(pnames[k++], "_b%d", i+1);
    }

    if (ainfo.misslist != null) {
	for (i=1; i<=ainfo.misslist[0]; i++) {
	    j = ainfo.misslist[i];
	    parms[k] = dset.Z[1][j];
	    sprintf(pnames[k++], "_c%d", i);
	}
    }

    /* construct NLS specification string */

    strcpy(fnstr.ptr, "y=".ptr);

    if (ainfo.ifc) {
	strcat(fnstr.ptr, "_b0".ptr);
    } else {
	strcat(fnstr.ptr, "0".ptr);
    }

    for (i=0; i<ainfo.p && !err; i++) {
	if (AR_included(ainfo, i)) {
	    lag = i + 1;
	    sprintf(term.ptr, "+_phi%d*".ptr, lag);
	    err = add_to_spec(fnstr.ptr, term.ptr);
	    if (!err) {
		err = y_Xb_at_lag(fnstr.ptr, ainfo, narmax, lag);
	    }
	}
    }

    for (j=0; j<ainfo.P && !err; j++) {
	sprintf(term.ptr, "+_Phi%d*".ptr, j+1);
	strcat(fnstr.ptr, term.ptr);
	lag = (j + 1) * ainfo.pd;
	y_Xb_at_lag(fnstr.ptr, ainfo, narmax, lag);
	for (i=0; i<ainfo.p; i++) {
	    if (AR_included(ainfo, i)) {
		sprintf(term.ptr, "-_phi%d*_Phi%d*".ptr, i+1, j+1);
		err = add_to_spec(fnstr.ptr, term.ptr);
		if (!err) {
		    lag = (j+1) * ainfo.pd + (i+1);
		    y_Xb_at_lag(fnstr.ptr, ainfo, narmax, lag);
		}
	    }
	}
    }

    for (i=0; i<ainfo.nexo && !err; i++) {
	sprintf(term.ptr, "+_b%d*x%d".ptr, i+1, i+1);
	err = add_to_spec(fnstr.ptr, term.ptr);
    }

    if (!err && ainfo.misslist != null) {
	for (i=1; i<=ainfo.misslist[0]; i++) {
	    sprintf(term.ptr, "+_c%d*d%d".ptr, i, i);
	    err = add_to_spec(fnstr.ptr, term.ptr);
	}
    }

    if (!err) {
	if (coeff == null) {
	    nls_kickstart(amod, dset, b0, by1);
	}

	err = nlspec_set_regression_function(spec, fnstr.ptr, dset);
    }

    if (!err) {
	double save_tol = libset_get_double(SetKey.nls_toler);

	libset_set_double(SetKey.nls_toler, 1.0e-5);
	set_auxiliary_scalars();
	err = aux_nlspec_add_param_list(spec, nparam, parms, pnames);
	if (!err) {
	    *amod = model_from_nlspec(spec, dset, nlsopt, prn);
	    err = amod.errcode;
	}
	unset_auxiliary_scalars();
	libset_set_double(SetKey.nls_toler, save_tol);
    }

 bailout:

    nlspec_destroy(spec);
    free(parms);
    strings_array_free(pnames, nparam);

    return err;
}

/* compose the regression list for the case where we're initializing
   ARMA via plain OLS (not NLS)
*/

private int * make_ar_ols_list(ArmaInfo * ainfo, int av)
{
    int *list = gretl_list_new(av);
    int i;
    int k;
    int vi;

    if (list == null) {
	return null;
    }

    list[1] = 1;

    if (ainfo.ifc) {
	list[2] = 0;
	k = 3;
    } else {
	list[0] -= 1;
	k = 2;
    }

    /* allow for const and y */
    vi = 2;

    for (i=0; i<ainfo.p; i++) {
	if (AR_included(ainfo, i)) {
	    list[k++] = vi++;
	}
    }

    for (i=0; i<ainfo.P; i++) {
	list[k++] = vi++;
    }

    for (i=0; i<ainfo.nexo; i++) {
	list[k++] = vi++;
    }

    return list;
}

/* Apply least squares to get initial values for the AR coefficients,
   either OLS or NLS.  We use NLS if there is nonlinearity due to
   either (a) the presence of both a seasonal and a non-seasonal AR
   component or (b) the presence of exogenous variables in the context
   of a non-zero AR order, where estimation will be via exact ML.

   In this initialization any MA coefficients are simply set to
   "near-zero" (MA_SMALL).
*/

extern (C) int ar_arma_init(double * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod, Gretlopt opt)
{
    int * list = ainfo.alist;
    int nmixed = ainfo.np * ainfo.P;
    int ptotal = ainfo.np + ainfo.P + nmixed;
    int av = ptotal + ainfo.nexo + 2;
    DatasetRaw * aset = null;
    int * arlist = null;
    Model armod;
    int narmax;
    int nonlin = 0;
    int i;
    int err = 0;


    if (ptotal == 0 && ainfo.nexo == 0 && !ainfo.ifc) {
	/* special case of pure MA model */
	for (i=0; i<ainfo.nq + ainfo.Q; i++) {
	    coeff[i] = MA_SMALL;
	}
	ainfo.init_ = IniMethod.small;
	return 0;
    }

    gretl_model_init(&armod, dset);

    narmax = arma_exact_ml(ainfo) ? ainfo.nexo : 0;
    if (narmax > 0 && ptotal > 0) {
	/* ARMAX-induced lags of exog vars */
	av += ainfo.nexo * ptotal;
    }

    if (ptotal == 0 && ainfo.nexo == 0 && ainfo.ifc) {
	/* straight MA model with constant */
	coeff[0] = get_y_mean(ainfo);
	for (i=1; i<=ainfo.nq + ainfo.Q; i++) {
	    coeff[i] = MA_SMALL;
	}
	ainfo.init_ = IniMethod.small;
	return 0;
    }

    aset = create_auxiliary_dataset(av, ainfo.fullT, Gretlopt.OPT_NONE);
    if (aset == null) {
	return GretlError.alloc;
    }

    if (ptotal > 0 && (narmax > 0 || nmixed > 0)) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	arlist = make_ar_ols_list(ainfo, av);
    }

    /* build temporary dataset, dset . aset */
    arma_init_build_dataset(ainfo, ptotal, narmax, list,
			    dset, aset, nonlin);

    if (nonlin) {
	Prn * dprn = null;

	err = arma_get_nls_model(&armod, ainfo, narmax, null, aset,
				 dprn);
    } else {
	armod = lsq(arlist, aset, GretlCmdIndex.ols, Gretlopt.OPT_A | Gretlopt.OPT_Z);
	err = armod.errcode;
    }


    if (!err) {
	ar_init_transcribe_coeffs(ainfo, &armod, coeff);
    }

    /* handle the case where we need to translate from an
       estimate of the regression constant to the
       unconditional mean of y_t
    */
    if (!err && (!nonlin || ainfo.nexo == 0) &&
	init_transform_const(ainfo)) {
	transform_arma_const(coeff, ainfo);
    }

    if (!err) {
	ainfo.init_ = nonlin ? IniMethod.nls : IniMethod.ols;
    }

    /* clean up */
    clear_model(&armod);
    destroy_dataset(aset);
    free(arlist);

    return err;
}

private int * plain_ols_list(ArmaInfo * ainfo)
{
    int * list;
    int nl = 1; /* y */
    int i;
    int j;

    if (ainfo.ifc) {
        nl++;
    }
    if (ainfo.xlist != null) {
        nl += ainfo.xlist[0];
    }
    list = gretl_list_new(nl);
    i = 1;
    list[i++] = ainfo.yno;
    if (ainfo.ifc) {
        list[i++] = 0;
    }
    if (ainfo.xlist != null) {
        for (j=1; j<=ainfo.xlist[0]; j++) {
            list[i++] = ainfo.xlist[j];
        }
    }

    return list;
}

private Gretlopt arma_ols_opt(ArmaInfo * ainfo)
{
    Gretlopt opt = Gretlopt.OPT_A | Gretlopt.OPT_Z;

    if (ainfo.nc == 0) {
        opt |= Gretlopt.OPT_U;
    }

    return opt;
}

extern (C) int arma_by_simple_ols(const(double) * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod)
{
    int * arlist;
    Gretlopt opt;

    arlist = plain_ols_list(ainfo);
    opt = arma_ols_opt(ainfo);
    *pmod = lsq(arlist, cast(DatasetRaw *) dset, GretlCmdIndex.ols, opt);

    return pmod.errcode;
}

extern (C) int arma_by_ls(const(double) * coeff, const(DatasetRaw) * dset, ArmaInfo * ainfo, Model * pmod)
{
    Prn * prn = ainfo.prn;
    int * list = ainfo.alist;
    int nmixed = ainfo.np * ainfo.P;
    int ptotal = ainfo.np + ainfo.P + nmixed;
    int av = ptotal + ainfo.nexo + 2;
    DatasetRaw * aset = null;
    int * arlist = null;
    int nonlin = 0;

    aset = create_auxiliary_dataset(av, ainfo.fullT, Gretlopt.OPT_NONE);
    if (aset == null) {
	return GretlError.alloc;
    }

    if (ptotal > 0 && nmixed > 0) {
	/* we'll have to use NLS */
	nonlin = 1;
    } else {
	/* OLS: need regression list */
	arlist = make_ar_ols_list(ainfo, av);
    }

    /* build temporary dataset */
    arma_init_build_dataset(ainfo, ptotal, 0, list,
			    dset, aset, nonlin);

    if (nonlin) {
	pmod.errcode = arma_get_nls_model(pmod, ainfo, 0, coeff, aset,
					   prn);
    } else {
	Gretlopt opt = arma_ols_opt(ainfo);

	*pmod = lsq(arlist, aset, GretlCmdIndex.ols, opt);
    }

    /* clean up */
    free(arlist);
    destroy_dataset(aset);

    if (!pmod.errcode && pmod.full_n < dset.n) {
	/* the model series are short */
	double * uhat = cast(double *) malloc(dset.n * double.sizeof);
	double * yhat = cast(double *) malloc(dset.n * double.sizeof);
	int s;
	int t;

	if (uhat == null || yhat == null) {
	    free(uhat);
	    free(yhat);
	    pmod.errcode = GretlError.alloc;
	} else {
	    for (t=0; t<dset.n; t++) {
		uhat[t] = yhat[t] = NADBL;
	    }
	    t = ainfo.t1;
	    for (s=0; s<pmod.full_n; s++, t++) {
		uhat[t] = pmod.uhat[s];
		yhat[t] = pmod.yhat[s];
	    }
	    free(pmod.uhat);
	    pmod.uhat = uhat;
	    free(pmod.yhat);
	    pmod.yhat = yhat;
	}
    }

    return pmod.errcode;
}
