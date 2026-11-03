/* Port of plugin/fcp.c: GARCH(p, q) estimation of a linear equation by the mixed-gradient algorithm of Fiorentini, Calzolari and
   Panattoni (Journal of Applied Econometrics, 1996); see Bollerslev, JoE 31 (1986), 307-327. garch_estimate is the estimator that
   gretl.plugin.garch calls with the option --fcp; garch_analytical_hessian is the analytical Hessian used by the other estimator
   (mod_garch.c). The debug output (FDEBUG, off in the C) is not in the port.

   The static variables of garch_info_matrix and garch_hessian (ll1, fs) are function-local __gshared variables, as in the C. */

module gretl.plugin.fcp;

import gretl.gretl_matrix : gretl_invert_symmetric_indef_matrix, gretl_invert_symmetric_matrix, gretl_matrix_alloc, gretl_matrix_copy_values, gretl_matrix_free, gretl_matrix_qform, gretl_matrix_switch_sign, gretl_matrix_zero, gretl_zero_matrix_new;
import gretl.gretl_prn : pprintf, pputc, pputs;
import gretl.gretl_utils : doubles_array_free, doubles_array_new;
import gretl.missingh : na;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import gretl.gretl_prn : Prn;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_modelh : MLVCVType;


/*
    The functions here are derived from the Fortran in VSGARCMX.FOR,
    from the Fiorentini, Calzolari and Panattoni GARCH implementation
    (Journal of Applied Econometrics, 1996), the original notice for
    which is reproduced infra.  The relevant code was extracted from
    its original Monte Carlo context, translated to C using f2c,
    then edited extensively to turn it into more idiomatic C.

    I have also modified the comments at certain places in the code,
    where reference is made to the equations in the FCP paper.  The
    comments in the Fortran apparently pertained to a draft of the
    paper; I have updated them relative to the paper as published in
    JAE, 1996, pp. 399-417.0

    Allin Cottrell, Wake Forest University, March 2004.0
    Updated and re-worked, June 2007.0
    Further minor modifications, September 2008.0
*/

/* Gabriele FIORENTINI, Giorgio CALZOLARI, Lorenzo PANATTONI
   Journal of Applied Econometrics, 1996

   mixed gradient algorithm

   garch(p,q) estimates of a linear equation

   See BOLLERSLEV, JOE 31(1986), 307-327.0
*/



enum SMALL_HT = 1.0e-7;
enum S1MIN = 1.0e-10;
enum SUMGRMAX = 1.0e-4;

enum FcpCode {
    FCP_FULL, /* doing the full job of GARCH estimation */
    FCP_HESS  /* just calculating the Hessian */
}


struct FcpInfo {
    int nc;
    int t1;
    int t2;
    int T;
    int p;
    int q;
    int npar;
    double scale;

    /* supplied arrays */
    const(double) * y;
    const(double) ** X;
    double * theta;
    double * e;
    double * e2;
    double * h;

    /* temporary arrays */
    double * grad;
    double * parpre;
    double * gg;
    double * step;
    double * zt;
    double * asum2;
    double **dhdp;
    double ***H;

    GretlMatrixRaw * V;
};

private void free_H(double *** H, int np)
{
    int i;
    int j;

    for (i=0; i<np; i++) {
	if (H[i] != null) {
	    for (j=0; j<np; j++) {
		free(H[i][j]);
	    }
	    free(H[i]);
	}
    }

    free(H);
}

private double *** allocate_H(int np, int p, int q)
{
    int i, j, lag = (p > q)? p : q;
    double ***H;

    H = cast(double** *) malloc(np * (*H).sizeof);
    if (H == null) {
	return null;
    }

    for (i=0; i<np; i++) {
	H[i] = null;
    }

    for (i=0; i<np; i++) {
	H[i] = cast(double **) malloc(np * (double *).sizeof);
	if (H[i] == null) {
	    goto bailout;
	}
	for (j=0; j<np; j++) {
	    H[i][j] = null;
	}
	for (j=0; j<np; j++) {
	    H[i][j] = cast(double *) malloc((lag + 1) * double.sizeof);
	    if (H[i][j] == null) {
		goto bailout;
	    }
	}
    }

    return H;

 bailout:

    free_H(H, np);
    return null;
}

private int fcp_allocate(FcpInfo * f, int code)
{
    f.zt = cast(double *) malloc((f.p + f.q + 1) * (*f.zt).sizeof);
    f.asum2 = cast(double *) malloc(f.nc * (*f.asum2).sizeof);
    f.grad = cast(double *) malloc(f.npar * (*f.grad).sizeof);
    if (f.zt == null || f.asum2 == null || f.grad == null) {
	return GretlError.alloc;
    }

    if (code == FcpCode.FCP_FULL) {
	f.parpre = cast(double *) malloc(f.npar * (*f.parpre).sizeof);
	f.gg = cast(double *) malloc(f.npar * (*f.gg).sizeof);
	f.step = cast(double *) malloc(f.npar * (*f.step).sizeof);
	if (f.parpre == null || f.gg == null || f.step == null) {
	    return GretlError.alloc;
	}
    }

    f.dhdp = doubles_array_new(f.npar, f.T);
    if (f.dhdp == null) {
	return GretlError.alloc;
    }

    f.V = gretl_zero_matrix_new(f.npar, f.npar);
    if (f.V == null) {
	return GretlError.alloc;
    }

    f.H = allocate_H(f.npar, f.p, f.q);
    if (f.H == null) {
	return GretlError.alloc;
    }

    return 0;
}

private void fcpinfo_destroy(FcpInfo * f)
{
    free(f.grad);
    free(f.parpre);
    free(f.gg);
    free(f.step);
    free(f.zt);
    free(f.asum2);

    doubles_array_free(f.dhdp, f.npar);
    gretl_matrix_free(f.V);
    free_H(f.H, f.npar);

    free(f);
}

private FcpInfo * fcpinfo_new(int q, int p, int t1, int t2, int T, const(double) * y, const(double) ** X, int nc, double * theta, double * e, double * e2, double * h, double scale, int code)
{
    FcpInfo * f = cast(FcpInfo *) malloc(FcpInfo.sizeof);

    if (f == null) {
	return null;
    }

    f.theta = theta;

    f.grad = null;
    f.parpre = null;
    f.gg = null;
    f.step = null;
    f.zt = null;
    f.asum2 = null;
    f.dhdp = null;
    f.V = null;

    f.nc = nc;
    f.t1 = t1;
    f.t2 = t2;
    f.T = T,
    f.p = p;
    f.q = q;
    f.y = y;
    f.X = X;
    f.e = e;
    f.e2 = e2;
    f.h = h;
    f.scale = scale;

    f.npar = f.nc + 1 + q + p;

    if (fcp_allocate(f, code)) {
	fcpinfo_destroy(f);
	f = null;
    }

    return f;
}

private void print_iter_val(double x, int i, int k, Prn * prn)
{
    if (na(x)) {
	pprintf(prn, "%-12s", "NA".ptr);
    } else {
	pprintf(prn, "%#12.5g", x);
    }
    if (i && i % 6 == 5 && i < k-1) {
	pprintf(prn, "\n%12s", " ".ptr);
    }
}

private void garch_iter_info(FcpInfo * f, int iter, double ll, int hess, Prn * prn)
{
    double x;
    int i;

    pprintf(prn, "\n*** %s %d%s\n", "iteration".ptr, iter + 1,
	    (hess)? _(" (using Hessian)").ptr : " (using Information Matrix)".ptr);


    pputs(prn, _("Parameters: ").ptr);

    for (i=0; i<f.npar; i++) {
	x = f.parpre[i];
	if (i < f.nc) {
	    x *= f.scale;
	} else if (i == f.nc) {
	    x *= f.scale * f.scale;
	}
	print_iter_val(x, i, f.npar, prn);
    }

    pputc(prn, '\n');
    pputs(prn, _("Gradients:  ").ptr);

    for (i=0; i<f.npar; i++) {
	print_iter_val(f.grad[i], i, f.npar, prn);
    }

    pprintf(prn, "\nll = %f\n", ll);
}

private double get_yhat(const(double) ** X, int n, int t, const(double) * a)
{
    double yhat = 0.0;
    int i;

    for (i=0; i<n; i++) {
	yhat += a[i] * X[i][t];
    }

    return yhat;
}

/* Compute the GARCH log-likelihood.  Params are passed in f.theta;
   e, e2 and ht are computed here (e2 holds squared residuals).
*/

private double garch_ll(FcpInfo * f)
{
    int t1 = f.t1;
    int t2 = f.t2;
    int p = f.p;
    int q = f.q;
    int nc = f.nc;
    int i;
    int t;
    int lag;
    int n = t2 - t1 + 1;
    double uncvar;
    double ll;
    double scale2;
    double hts;

    const(double) * alpha = f.theta + nc + 1;
    const(double) * beta = alpha + q;


    /* Compute residuals, squared residuals, and unconditional
       variance over the real estimation period
    */
    uncvar = 0.0;
    for (t = t1; t <= t2; t++) {
	f.e[t] = f.y[t] - get_yhat(f.X, f.nc, t, f.theta);
	f.e2[t] = f.e[t] * f.e[t];
	uncvar += f.e2[t];
    }
    uncvar /= n;


    /*
       We use sample unconditional variance as the starting value
       (at time 0, -1, -2, etc.) for the squared residuals and ht;
       we use 0 as the starting value for residuals.
    */

    lag = (p > q)? p : q;

    for (t = t1-lag; t < t1; ++t) {
	f.e[t] = 0.0;
	f.e2[t] = f.h[t] = uncvar;
    }

    for (t=t1; t<=t2; t++) {
	f.h[t] = f.theta[nc];
	for (i=1; i<=q; i++) {
	    f.h[t] += f.e2[t-i] * alpha[i-1];
	}
	for (i=1; i<=p; i++) {
	    f.h[t] += f.h[t-i] * beta[i-1];
	}
	/* arbitrary */
	if (f.h[t] <= 0.0) {
	    f.h[t] = SMALL_HT;
	}
    }

    ll = 0.0;
    scale2 = f.scale * f.scale;


    for (t=t1; t<=t2; t++) {
	hts = f.h[t] * scale2;
	ll -= 0.5 * log(hts) + 0.5 * f.e2[t] / f.h[t] + LN_SQRT_2_PI;
    }

    return ll;
}

/* combined setup for OP matrix, information matrix and Hessian */

private int vcv_setup(FcpInfo * f, GretlMatrixRaw * V, int code)
{
    int i;
    int j;
    int k;
    int t;
    int n;
    int lag;
    int nvpar;
    double x;

    int t1 = f.t1;
    int t2 = f.t2;
    int p = f.p;
    int q = f.q;
    int nc = f.nc;
    int npar = f.npar;

    double **dhdp = f.dhdp;
    const double **g = f.X;
    double ***H = f.H;
    double * e = f.e;
    double * e2 = f.e2;
    double * h = f.h;
    double * zt = f.zt;
    double * asum2 = f.asum2;

    const(double) * alpha = f.theta + nc + 1;
    const(double) * beta = alpha + q;


    /* some useful abbreviations */
    nvpar = 1 + q + p;
    lag = (p > q)? p : q;
    n = t2 - t1 + 1;


    /* Begin computation of dhtdp wrt the variance parameters; we
       start computing derivatives of starting values; for ht starting
       values are obtained from the unconditional variance of the
       residuals.
     */

    for (k=1; k<=p; k++) {
	for (i=0; i<nvpar; i++) {
	    dhdp[nc+i][t1-k] = 0.0;
	    if (H != null) {
		/* hessian only */
		for (j=0; j<nvpar; j++) {
		    H[nc+i][nc+j][k] = 0.0;
		}
	    }
	}
    }

    for (t=t1; t<=t2; t++) {
	/* fill in zt at time t (see p. 401) */
	zt[0] = 1.0;
	for (i=1; i<=q; i++) {
	    zt[i] = e2[t-i];
	}
	for (i=1; i<=p; i++) {
	    zt[q+i] = h[t-i];
	}

	/* Fill in dhtdp at time t, part relative to variance parameters
	   (eq. 7, p. 402) */
	for (i=0; i<nvpar; i++) {
	    dhdp[nc+i][t] = zt[i];
	    for (j=1; j<=p; j++) {
		dhdp[nc+i][t] += dhdp[nc+i][t-j] * beta[j-1];
	    }
	}
    }

    /* Build matrix dhtdp, block for regression coefficients
       (eq. 13). We use 0 as starting value (time 0, -1, etc.) for the
       derivatives of ht wrt the coefficients; we also use 0 as
       starting values for the residuals.
    */

    /* building blocks for pre-sample terms */
    for (i=0; i<nc; i++) {
	asum2[i] = 0.0;
	for (t=t1; t<=t2; t++) {
	    asum2[i] -= e[t] * 2.0 * g[i][t];
	}
	asum2[i] /= n;
    }

    /* pre-sample range */
    for (t=t1-lag; t<t1; t++) {
	for (i=0; i<nc; i++) {
	    dhdp[i][t] = asum2[i];
	}
    }

    /* actual sample range */
    for (t=t1; t<=t2; t++) {
	for (i=0; i<nc; i++) {
	    dhdp[i][t] = 0.0;
	    for (j=1; j<=q; j++) {
		if (t - q < t1) {
		    dhdp[i][t] += alpha[j-1] * asum2[i];
		} else {
		    dhdp[i][t] -= alpha[j-1] * 2.0 * g[i][t-j] * e[t-j];
		}
	    }
	    for (j=1; j<=p; j++) {
		dhdp[i][t] += dhdp[i][t-j] * beta[j-1];
	    }
	}
    }

    /* Initialize gradient and vcv */
    for (i=0; i<npar; i++) {
	f.grad[i] = 0.0;
    }
    gretl_matrix_zero(V);

    for (t=t1; t<=t2; t++) {
	double r_h = e[t] / h[t];
	double r2_h = e[t] * r_h;
	double aa;
	double bb;
	int ncj;
	int nci;

	/*
	   First part, relative to regression coefficients (eq. 10, p. 402)
 	*/
	for (i=0; i<nc; i++) {
	    aa = r_h * g[i][t] + 0.5 / h[t] * dhdp[i][t] * (r2_h - 1.0);
	    f.grad[i] += aa;
	    if (code == MLVCVType.op) {
		for (j=0; j<=i; j++) {
		    bb = r_h * g[j][t] + 0.5 / h[t] * dhdp[j][t] * (r2_h - 1.0);
		    x = gretl_matrix_get(V, i, j);
		    x += aa * bb;
		    gretl_matrix_set(V, i, j, x);
		    gretl_matrix_set(V, j, i, x);
		}
		for (j=0; j<nvpar; j++) {
		    ncj = nc + j;
		    x = gretl_matrix_get(V, i, ncj);
		    x += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    gretl_matrix_set(V, i, ncj, x);
		    gretl_matrix_set(V, ncj, i, x);
		}
	    }
	}

	/*
	   Second part, relative to variance parameters (eq. 6, p. 401)
	*/
	for (i=0; i<nvpar; i++) {
	    nci = nc + i;
	    aa = 0.5 / h[t] * dhdp[nci][t] * (r2_h - 1.0);
	    f.grad[nci] += aa;
	    if (code == MLVCVType.op) {
		for (j=0; j<=i; j++) {
		    ncj = nc + j;
		    x = gretl_matrix_get(V, nci, ncj);
		    x += aa * 0.5 / h[t] * dhdp[ncj][t] * (r2_h - 1.0);
		    gretl_matrix_set(V, nci, ncj, x);
		    gretl_matrix_set(V, ncj, nci, x);
		}
	    }
	}
    }

    if (code == MLVCVType.im) {
	for (t=t1; t<=t2; t++) {
	    double ht2 = h[t] * h[t];

	    /* Part relative to the coefficients (eq. 30, p. 406).
	       Since we take the expected value, only the first two terms
	       remain.
	    */
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    x = gretl_matrix_get(V, i, j);
		    x -= g[i][t] * g[j][t] / h[t] + 0.5 * dhdp[i][t] * dhdp[j][t] / ht2;
		    gretl_matrix_set(V, i, j, x);
		}
	    }

	    /* Part relative to the variance parameters (eq. 29, p. 406).
	       Since we take the expected value, only the second term
	       remains.
	    */
	    for (i=nc; i<npar; i++) {
		for (j=nc; j<npar; j++) {
		    x = gretl_matrix_get(V, i, j);
		    x -= 0.5 * dhdp[i][t] * dhdp[j][t] / ht2;
		    gretl_matrix_set(V, i, j, x);
		}
	    }
	}
    }

    if (code != MLVCVType.hessian) {
	/* we're done */
	return 0;
    }

    /* Initial values in dhdpdp (here, "H") are 2/t x'x, and zero for
       off-diagonal (mixed) blocks. */

    for (k=0; k<lag; k++) {
	for (i=0; i<nc; i++) {
	    for (j=0; j<nc; j++) {
		H[i][j][k+1] = 0.0;
	    }
	}
	for (t=t1; t<=t2; t++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    H[i][j][k+1] += 2.0 * g[i][t] * g[j][t] / n;
		}
	    }
	}
	for (i=0; i<nc; i++) {
	    for (j=0; j<nvpar; j++) {
		/* mod. by AC: zero _all_ mixed entries */
		H[i][nc+j][k+1] = H[nc+j][i][k+1] = 0.0;
	    }
	}
    }

    /* Now we fill out the full Hessian */

    for (t=t1; t<=t2; ++t) {
	double r_h = e[t] / h[t];
	double r2_h = r_h * e[t];
	double r2_h3 = r2_h / (h[t] * h[t]);
	double u_h2 = 1.0 / (h[t] * h[t]);

	for (i=0; i<npar; i++) {
	    for (j=0; j<npar; j++) {
		H[i][j][0] = 0.0;
	    }
	}

	if (lag <= 0) {
	    goto lag0;
	}

	for (k=1; k<=q; k++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    if (t - q < t1) {
			H[i][j][0] += H[i][j][q] * alpha[k-1];
		    } else {
			H[i][j][0] += 2.0 *
			    g[i][t-k] * g[j][t-k] * alpha[k-1];
		    }
		}
	    }
	}

	for (k=1; k<=p; k++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nc; j++) {
		    H[i][j][0] += H[i][j][k] * beta[k-1];
		}
	    }
	}

	for (i=0; i<nc; i++) {
	    for (k=1; k<=q; k++) {
		if (t - q < t1) {
		    H[i][nc+k][0] += asum2[i];
		} else {
		    H[i][nc+k][0] -= 2.0 * g[i][t-k] * e[t-k];
		}
	    }
	    for (k=1; k<=p; k++) {
		H[i][nc+q+k][0] += dhdp[i][t-k];
	    }
	}

	for (k=1; k<=p; k++) {
	    for (i=0; i<nc; i++) {
		for (j=0; j<nvpar; j++) {
		    H[i][nc+j][0] += H[i][nc+j][k] * beta[k-1];
		}
	    }
	}

    lag0:

	/* Part relative to the coefficients (eq. 15, p. 403).
	   Since we take the expected value, only the first two terms
	   remain.
 	*/

	for (i=0; i<nc; i++) {
	    for (j=0; j<nc; j++) {
		x = gretl_matrix_get(V, i, j);
		x = x - g[i][t] * g[j][t] / h[t]
		    - 0.5 * r2_h3 * dhdp[i][t] * dhdp[j][t]
		    - (r_h * g[j][t] * dhdp[i][t]) / h[t]
		    - (r_h * g[i][t] * dhdp[j][t]) / h[t]
		    + 0.5 * (r2_h - 1.0) *
		    (H[i][j][0] / h[t] - dhdp[i][t]
		     * dhdp[j][t] / (h[t] * h[t]));
		gretl_matrix_set(V, i, j, x);
	    }
	}

	/* Part relative to the variance parameters (eq. 14, p. 403).
	   Since we take the expected value, only the second term
	   remains.
 	*/

	if (p > 0) {
	    for (i=0; i<nvpar; i++) {
		for (j=1; j<=p; j++) {
		    H[nc+i][nc+q+j][0] += dhdp[nc+i][t-j];
		}
	    }
	    for (i=1; i<=p; i++) {
		for (j=0; j<nvpar; j++) {
		    H[nc+q+i][nc+j][0] += dhdp[nc+j][t-i];
		}
	    }
	    for (k=1; k<=p; k++) {
		for (i=0; i<nvpar; i++) {
		    for (j=0; j<nvpar; j++) {
			H[nc+i][nc+j][0] += H[nc+i][nc+j][k] * beta[k-1];
		    }
		}
	    }
	}

	for (i=nc; i<npar; i++) {
	    for (j=nc; j<npar; j++) {
		x = gretl_matrix_get(V, i, j);
		x = x + 0.5 * u_h2 * dhdp[i][t] * dhdp[j][t]
		    - r2_h3 * dhdp[i][t] * dhdp[j][t]
		    + 0.5 * (r2_h - 1.0) / h[t] * H[i][j][0];
		gretl_matrix_set(V, i, j, x);
	    }
	}

	/* top-right mixed part (eq. 17, p. 403) */
	for (i=0; i<nc; i++) {
	    for (j=0; j<nvpar; j++) {
		x = gretl_matrix_get(V, i, nc+j);
		x = x - g[i][t] * r_h * dhdp[nc+j][t] / h[t]
		    - 0.5 * (r2_h - 1.0) * dhdp[nc+j][t] * dhdp[i][t] / (h[t] * h[t])
		    + 0.5 * (r2_h - 1.0) * H[i][nc+j][0] / h[t]
		    - 0.5 * r2_h * u_h2 * dhdp[i][t] * dhdp[nc+j][t];
		gretl_matrix_set(V, i, nc+j, x);
		/* and bottom left too */
		gretl_matrix_set(V, nc+j, i, x);
	    }
	}

	/* before quitting time t, tidy up dhdpdp */
	for (k=0; k<lag; k++) {
	    for (i=0; i<npar; i++) {
		for (j=0; j<npar; j++) {
		    H[i][j][lag-k] = H[i][j][lag-k-1];
		}
	    }
	}
    }

    return 0;
} /* vcv_setup */

/* Update theta using the current step vector.  Check that the
   conditional variance parameters are admissible.
 */

private void update_theta(FcpInfo * f, double d)
{
    int i;
    int nv = f.p + f.q + 1;
    double * a = f.theta + f.nc;
    double sum = 0.0;

    for (i=0; i<f.npar; i++) {
	f.theta[i] = f.gg[i] + f.step[i] * d;
    }

    if (a[0] <= 0.0) {
	a[0] = SMALL_HT;
    }

    for (i=1; i<nv; i++) {
	if (a[i] < 0.0) {
	    a[i] = 0.0;
	}
	sum += a[i];
    }

    if (sum > 1.0) {
	for (i=1; i<nv; i++) {
	    a[i] /= sum;
	}
    }
}

/* calculate the step for the new coefficients */

private double step_calc(FcpInfo * f, const(GretlMatrixRaw) * V, double * ds, double * ps2)
{
    double s1 = 0.0;
    double s2 = 0.0;
    double relstep;
    int i;
    int j;

    for (i=0; i<f.npar; i++) {
	s1 += f.theta[i] * f.theta[i];
	f.gg[i] = f.theta[i];
	f.step[i] = 0.0;
	for (j=0; j<f.npar; j++) {
	    f.step[i] -= f.grad[j] * gretl_matrix_get(V, i, j);
	}
	s2 += f.step[i] * f.step[i];
    }

    if (s1 == 0.0) {
	s1 = S1MIN;
    }

    relstep = sqrt(s2 / s1);
    s2 = sqrt(s2);

    for (i=0; i<f.npar; i++) {
	f.step[i] /= s2;
    }

    *ds = *ps2 = s2;

    return relstep;
}

/* Common iteration routine called when adjusting the parameters,
   either using the information matrix or the Hessian.  The flow
   control here is fairly inscrutable for ordinary mortals, but
   it is very effective!
*/

private void fcp_iterate(FcpInfo * f, GretlMatrixRaw * V, double * pll1, double * pfs, double toler, int count)
{
    double ll2;
    double ll3;
    double ll1 = *pll1;
    double fs = *pfs;
    double d0;
    double d1;
    double d2;
    double d3;
    double d12;
    double d31;
    double d23;
    double d12s;
    double d31s;
    double d23s;
    double di;
    double dm;
    double ds;
    double dmin;
    double dmax;
    double a1s;
    double a2s;
    double a3s;
    double s2;
    double bigd;
    double step;
    int nexp;
    int ncall = 0;

    step = step_calc(f, V, &ds, &s2);

    if (step <= toler) {
	update_theta(f, ds);
	*pll1 = fs;
	*pfs = -fs;
	return;
    }

    nexp = count / 5;
    if (nexp > 5) {
	nexp = 5;
    }

    d0 = s2 / pow(2.0, nexp);
    dmin = d0 * 0.001;
    dmax = d0 * 4.0;

    if (count == 1) {
	ll1 = -garch_ll(f);
    }

    update_theta(f, d0);
    ll2 = -garch_ll(f);

    if (ll2 > ll1) {
	d1 = -d0;
	d2 = 0.0;
	d3 = d0;
	ll3 = ll2;
	ll2 = ll1;
	update_theta(f, d1);
	ll1 = -garch_ll(f);
    } else {
	d1 = 0.0;
	d2 = d0;
	d3 = d0 + d0;
	update_theta(f, d3);
	ll3 = -garch_ll(f);
    }

    while (1) {
	d23 = d2 - d3;
	d31 = d3 - d1;
	d12 = d1 - d2;
	di = d23 * ll1 + d31 * ll2 + d12 * ll3;
	bigd = di * -2.0 / (d23 * d31 * d12);
	if (bigd > 0.0) {
	    goto LLS;
	}
	if (ll3 <= ll1) {
	    goto LL3;
	}

    LL1:
	d3 = d2;
	d2 = d1;
	ll3 = ll2;
	ll2 = ll1;
	d1 -= dmax;
	update_theta(f, d1);
	ll1 = -garch_ll(f);
	if (++ncall > 100) {
	    break;
	}
	continue;

    LL3:
	d1 = d2;
	d2 = d3;
	ll1 = ll2;
	ll2 = ll3;
	d3 += dmax;
	update_theta(f, d3);
	ll3 = -garch_ll(f);
	if (++ncall > 100) {
	    break;
	}
	continue;

    LLS:
	d23s = d23 * (d2 + d3);
	d31s = d31 * (d3 + d1);
	d12s = d12 * (d1 + d2);
	ds = (d23s * ll1 + d31s * ll2 + d12s * ll3) * 0.5 / di;
	update_theta(f, ds);
	fs = -garch_ll(f);
	if (++ncall > 100) {
	    break;
	}

	a1s = fabs(d1 - ds);
	a2s = fabs(d2 - ds);
	a3s = fabs(d3 - ds);

	dm = (a1s > a3s)? a3s : a1s;
	if (dm > dmax) {
	    if (ds < d1 - dmax) {
		goto LL1;
	    }
	    if (ds > d3 + dmax) {
		goto LL3;
	    }
	}

	if (a1s < dmin || a2s < dmin || a3s < dmin) {
	    break;
	}

	if (ll1 < ll2 || ll1 < ll3) {
	    if (ll2 < ll3 || ll2 < ll1) {
		d3 = ds;
		ll3 = fs;
	    } else {
		d2 = ds;
		ll2 = fs;
	    }
	} else {
	    d1 = ds;
	    ll1 = fs;
	}

	while (d1 > d2 || d2 > d3) {
	    double tmp;

	    if (d2 > d3) {
		tmp = d2;
		d2 = d3;
		d3 = tmp;
		tmp = ll2;
		ll2 = ll3;
		ll3 = tmp;
	    }
	    if (d1 > d2) {
		tmp = d1;
		d1 = d2;
		d2 = tmp;
		tmp = ll1;
		ll1 = ll2;
		ll2 = tmp;
	    }
	}
    }

    if (fs > ll1) {
	fs = ll1;
	ds = d1;
    }
    if (fs > ll2) {
	fs = ll2;
	ds = d2;
    }
    if (fs > ll3) {
	fs = ll3;
	ds = d3;
    }

    update_theta(f, ds);

    *pll1 = fs;
    *pfs = -fs;
}

private int garch_info_matrix(FcpInfo * f, GretlMatrixRaw * V, double toler, int * count)
{
    static __gshared double ll1 = 0.0;
    static __gshared double fs = 0.0;
    int err;

    vcv_setup(f, V, MLVCVType.im);

    if (count != null) {
	*count += 1;
    }

    err = gretl_invert_symmetric_indef_matrix(V);
    if (err) {
	fprintf(stderr, "garch_info_matrix: matrix inversion failed\n");
	return err;
    }

    if (count != null) {
	/* not just calculating vcv at convergence */
	fcp_iterate(f, V, &ll1, &fs, toler, *count);
    }

    gretl_matrix_switch_sign(V);

    return err;
}

private int garch_hessian(FcpInfo * f, GretlMatrixRaw * V, double toler, int * count)
{
    static __gshared double ll1 = 0.0;
    static __gshared double fs = 0.0;
    int i;
    int sign_done = 0;
    int err;

    vcv_setup(f, V, MLVCVType.hessian);

    if (count != null) {
	*count += 1;
    }

    if (toler == 0.0) {
	for (i=0; i<V.rows; i++) {
	    if (gretl_matrix_get(V, i, i) < 0.0) {
		gretl_matrix_switch_sign(V);
		sign_done = 1;
		break;
	    }
	}
    }

    if (toler == 0.0) {
	err = gretl_invert_symmetric_matrix(V);
    } else {
	err = gretl_invert_symmetric_indef_matrix(V);
    }

    if (err) {
	fprintf(stderr, "garch_hessian: matrix inversion failed\n");
	return err;
    }

    if (count != null) {
	fcp_iterate(f, V, &ll1, &fs, toler, *count);
    }

    if (!sign_done) {
	gretl_matrix_switch_sign(V);
    }

    return err;
}

private int make_garch_vcv(FcpInfo * f, const(GretlMatrixRaw) * ihess, GretlMatrixRaw * V, int vopt)
{
    GretlMatrixRaw * OP = null;
    GretlMatrixRaw * iinfo = null;
    int k = f.npar;
    int err = 0;

    /* OP and robust variants need OP matrix */
    if (vopt == MLVCVType.op || vopt == MLVCVType.qml || vopt == MLVCVType.bw) {
	OP = gretl_matrix_alloc(k, k);
	if (OP == null) {
	    err = GretlError.alloc;
	    goto bailout;
	}

	vcv_setup(f, OP, MLVCVType.op);

	if (vopt == MLVCVType.op) {
	    gretl_matrix_copy_values(V, OP);
	    err = gretl_invert_symmetric_matrix(V);
	}
    }

    /* IM and BW variants need the info matrix */
    if (vopt == MLVCVType.im || vopt == MLVCVType.bw) {
	iinfo = gretl_matrix_alloc(k, k);
	if (iinfo == null) {
	    err = GretlError.alloc;
	    goto bailout;
	}

	garch_info_matrix(f, iinfo, 0.0, null);

	if (vopt == MLVCVType.im) {
	    gretl_matrix_copy_values(V, iinfo);
	} else {
	    /* Bollerslev-Wooldridge */
	    gretl_matrix_qform(iinfo, GretlMatrixMod.none,
			       OP, V, GretlMatrixMod.none);
	}
    } else if (vopt == MLVCVType.qml) {
	gretl_matrix_qform(ihess, GretlMatrixMod.none,
			   OP, V, GretlMatrixMod.none);
    } else if (vopt == MLVCVType.hessian) {
	gretl_matrix_copy_values(V, ihess);
    }

 bailout:

    gretl_matrix_free(OP);
    gretl_matrix_free(iinfo);

    return err;
}

private int converged(FcpInfo * f, double tol)
{
    double s1 = 0.0;
    double s2 = 0.0;
    double pdiff;
    int i;

    for (i=0; i<f.npar; i++) {
	s1 += f.parpre[i] * f.parpre[i];
	pdiff = f.theta[i] - f.parpre[i];
	s2 += pdiff * pdiff;
    }

    if (s1 == 0.0) {
	s1 = S1MIN;
    }

    return s2 / s1 <= tol * tol;
}

/*
   Parameters to garch_estimate()

   y:     the dependent variable
   X:     the independent regressors (including the constant)
   t1:    beginning of sample relative to the arrays y, X
   t2:    end of sample
   nobs:  total number of observations in y, X
   nc:    number of columns in X
   p:     number of 'beta' variance params
   q:     number of 'alpha' variance params (excluding the constant)
   theta: full parameter vector, pre-initialized; on output, the
          estimates at convergence.
   V:     covariance matrix of parameters (all 0 on input)
   e:     vector of 0's on input, residuals on output
   e2:    storage for squared residuals on output
   h:     storage for conditional variances on output
   scale: factor used to scale the dependent variable
   pll:   location to receive log-likelihood on output
   iters: 0 on input, holds number of iterations on output
   vopt:  code indicating which version of the covariance
          matrix to compute in V
   prn:   print handle for info on iterations etc.

*/

extern (C) int garch_estimate(const(double) * y, const(double) ** X, int t1, int t2, int nobs, int nc, int p, int q, double * theta, GretlMatrixRaw * V, double * e, double * e2, double * h, double scale, double * pll, int * iters, int vopt, Prn * prn)
{
    FcpInfo *f;
    int it1;
    int it2;
    int ittot;
    int count = 0;
    int npar = nc + 1 + p + q;
    double tol1 = 0.05; /* tolerance when using info matrix */
    double tol2 = 1e-8; /* tolerance when using Hessian */
    double ll;
    double sumgra;
    int i;
    int err = 0;

    f = fcpinfo_new(q, p, t1, t2, nobs, y, X, nc,
		    theta, e, e2, h, scale, FcpCode.FCP_FULL);
    if (f == null) {
	return GretlError.alloc;
    }

    /* Step 1: iterate to a first approximation using the info matrix */

    for (it1=0; it1<100; it1++) {
	ll = garch_ll(f);
	for (i=0; i<npar; i++) {
	    f.parpre[i] = f.theta[i];
	}

	err = garch_info_matrix(f, f.V, tol1, &count);
	if (err) {
	    goto garch_exit;
	}

	garch_iter_info(f, it1, ll, 0, prn);

	if (converged(f, tol1)) {
	    break;
	}
    }

    ittot = it1 + 1;

    /* Step 2: fine-tune using the Hessian */

    for (it2=0; it2<100; it2++) {
	ll = garch_ll(f);
	for (i=0; i<npar; i++) {
	    f.parpre[i] = f.theta[i];
	}

	err = garch_hessian(f, f.V, tol2, &it2);
	if (err) {
	    goto garch_exit;
	}

	garch_iter_info(f, ittot++, ll, 1, prn);

	if (converged(f, tol2)) {
	    break;
	}
    }

    *iters = ittot;

    sumgra = 0.0;
    for (i=0; i<npar; i++) {
	sumgra += f.grad[i] * f.grad[i];
    }

    if (sumgra >= SUMGRMAX) {
	pprintf(prn, _("\nParameters and gradients at iteration %d:\n\n").ptr,
		ittot);
	for (i=0; i<f.npar; i++) {
	    pprintf(prn, "%12.6f (%9.6f)\n", f.theta[i], f.grad[i]);
	}
	pprintf(prn, _("\nSum of squared gradients = %0.9g (should be less " ~
			"than %g)\n").ptr, sumgra, SUMGRMAX);
	if (!err) {
	    err = GretlError.noconv;
	}
    } else {
	pprintf(prn, _("\nFull Hessian convergence at iteration %d, " ~
			"tol = %0.9g\n\n").ptr, ittot, tol2);
	*pll = ll;
    }

    if (!err) {
	/* build the desired VCV variant */
	err = make_garch_vcv(f, f.V, V, vopt);
    }

 garch_exit:

    fcpinfo_destroy(f);

    return err;
}

extern (C) GretlMatrixRaw * garch_analytical_hessian(const(double) * y, const(double) ** X, int t1, int t2, int nobs, int nc, int p, int q, double * theta, double * e, double * e2, double * h, double scale, int * err)
{
    GretlMatrixRaw * H = null;
    FcpInfo *f;

    f = fcpinfo_new(q, p, t1, t2, nobs, y, X, nc,
		    theta, e, e2, h, scale, FcpCode.FCP_HESS);
    if (f == null) {
	*err = GretlError.alloc;
	return null;
    }

    *err = garch_hessian(f, f.V, 0.0, null);

    if (!*err) {
	H = f.V;
	f.V = null;
    }

    fcpinfo_destroy(f);

    return H;
}
