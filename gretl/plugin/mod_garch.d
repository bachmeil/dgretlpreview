/* Port of plugin/mod_garch.c: modular GARCH routines by Jack Lucchetti (October 2006), meant to replace fcp.c seamlessly: the log-likelihood
   of the Gaussian GARCH model through its residual and variance recursions with their derivatives (garch_etht), the analytical score,
   the covariance matrix variants (inverse of the information matrix, OPG, QML, Bollerslev-Wooldridge; the analytical Hessian comes
   from gretl.plugin.fcp) and the maximization by BFGS (garch_estimate_mod, the default estimator of gretl.plugin.garch).

   The functions passed to BFGS_max (garch_loglik, garch_score) are extern (C). The unused function mark (inside #if 0) and the debug output
   (GDEBUG, MOD_DEBUG, off in the C) are not in the port. The anonymous enums are InitVar and GarchDist; the struct garch_container is
   GarchContainer. */

module gretl.plugin.mod_garch;

import gretl.gretl_bfgs : BFGS_defaults, BFGS_max, newton_raphson_max;
import gretl.gretl_matrix : gretl_invert_symmetric_matrix, gretl_matrix_alloc, gretl_matrix_copy_values, gretl_matrix_free, gretl_matrix_qform;
import gretl.gretl_utils : doubles_array_free, doubles_array_new;
import gretl.libset : libset_get_int;
import gretl.missingh : NADBL, na;

import gretl.gretl_bfgsh : BFGS_CRIT_FUNC;

import core.stdc.math : fabs, floor, isnan, log, sqrt, pow, exp;
import core.stdc.stdlib : free, malloc, calloc, realloc;
import core.stdc.stdio : fprintf, printf, sprintf, stderr;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import gretl.gretl_prn : Prn;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.libseth : OptimCode, SetKey;
import gretl.nlsh : CriterionType;
import gretl.libgretlh;
import gretl.gretl_matrixh;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_modelh : MLVCVType;
import gretl.plugin.fcp : garch_analytical_hessian;


/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/*
  Modular GARCH routines by Jack Lucchetti, October 2006. For the 
  moment, meant to replace seamlessly fcp.c, but syntax should evolve in 
  the future.
*/



enum InitVar {
    INIT_VAR_THEO,
    INIT_VAR_OLS,
    INIT_VAR_RESID
}

enum GarchDist {
    DIST_NORM,
    DIST_T
}


struct GarchContainer {
    const(double) * y; /* dependent variable */
    const(double) ** X;    /* regressors (constant excluded) */
    int t1; /* beginning of sample */
    int t2; /* end of sample */
    int nobs; /* number of observations */
    int ncm; /* number of regressors, including the constant */
    int p; /* GARCH p */
    int q; /* GARCH q */
    int k; /* total number of parameters */
    int init_; /* h0 initialisation method (C name: init) */
    int distrib; /* innovations distribution (only Gaussian for now) */
    double * e; /* residuals */
    double * e2; /* squared residuals */
    double * h; /* conditional variance */
    double **score_e;    /* derivatives of the residuals wrt the parameters */
    double **score_h;    /* derivatives of the variances wrt the parameters */
    double **blockglue;  /* derivatives of the loglik wrt residuals and variances */
    double **G;          /* score matrix */
    double * tot_score; /* score vector (sum of G) */
    double scale; /* scale factor for dependent var */
    int boundcheck; /* enable bounds check */
};


private void free_eh_derivs(GarchContainer * DH)
{
    doubles_array_free(DH.score_e, DH.k);
    doubles_array_free(DH.score_h, DH.k);
    doubles_array_free(DH.G, DH.k);
    doubles_array_free(DH.blockglue, 2);
}

private int allocate_eh_derivs(GarchContainer * DH)
{
    int k = DH.k;
    int n = DH.nobs;
    int err = 0;

    DH.score_e = doubles_array_new(k, n);
    DH.score_h = doubles_array_new(k, n);
    DH.G = doubles_array_new(k, n);
    DH.blockglue = doubles_array_new(2, n);

    if (DH.score_e == null ||
	DH.score_h == null ||
	DH.G == null ||
	DH.blockglue == null) {
	free_eh_derivs(DH);
	err = GretlError.alloc;
    } 

    return err;
}

private GarchContainer * garch_container_new(const(double) * y, const(double) ** X, int t1, int t2, int nobs, int nc, int p, int q, int init_method, double * e, double * e2, double * h, double scale)
{
    GarchContainer * DH = cast(GarchContainer *) malloc(GarchContainer.sizeof);

    if (DH == null) {
	return null;
    }

    DH.y = y;
    DH.X = X;
    DH.t1 = t1;
    DH.t2 = t2;
    DH.nobs = nobs;
    DH.ncm = nc;
    DH.p = p;
    DH.q = q;
    DH.init_ = init_method;
    DH.e = e;
    DH.e2 = e2;
    DH.h = h;
    DH.k = nc + 1 + p + q;
    DH.boundcheck = 1;
    DH.scale = scale;

    DH.score_e = null;
    DH.score_h = null;
    DH.G = null;
    DH.blockglue = null;

    if (allocate_eh_derivs(DH)) {
	free(DH);
	DH = null;
    }

    return DH;
}

private void garch_container_destroy(GarchContainer * DH)
{
    free_eh_derivs(DH);
    free(DH);
}

/* *ARCH log-likelihood for Gaussian innovations */

private double normal_ll(GarchContainer * DH)
{
    double e2t;
    double ht;
    double ll = 0.0;
    int t;

    for (t=DH.t1; t<=DH.t2; t++) {
	e2t = DH.e2[t];
	ht = DH.h[t];
	if (na(e2t) || na(ht)) {
	    return NADBL;
	}
	ll -= log(ht) + e2t / ht;
    }

    ll *= 0.5;
    ll -= (DH.t2 - DH.t1 + 1) * LN_SQRT_2_PI;

    return ll;
} 

private void normal_score(GarchContainer * DH)
{
    double ut;
    int t;

    for (t=DH.t1; t<=DH.t2; t++) {
	DH.blockglue[0][t] = ut = -DH.e[t] / DH.h[t];
	DH.blockglue[1][t] = 0.5 * (ut * ut - 1.0 / DH.h[t]);
    }
} 

private int params_in_bounds(const(double) * par, int ncm, int k)
{
    int ok = (par[ncm] >= 0.0); /* omega */
    double sum = 0.0;
    int i;

    for (i=ncm+1; i<k && ok; i++) {
	ok &= (par[i] >= 1.0e-12);
	sum += par[i];
    }

    ok &= (sum <= 1.0);

    return ok;
}

/* Compute the GARCH quantities */

private int garch_etht(const(double) * par, GarchContainer * DH)
{
    double **dedq = DH.score_e;
    double **dhdq = DH.score_h;
    int ncm = DH.ncm;
    int t1 = DH.t1;
    int t2 = DH.t2;
    int p = DH.p;
    int q = DH.q;
    int maxlag = (p > q)? p : q;
    int t0 = t1 - maxlag;
    int t;
    int T = t2 - t1 + 1;
    double et;
    double ht;
    double tmp;
    double h0 = 0.0;
    double dh0;
    double u_var = 0.0;
    int i;
    int j;
    int k;

    /* check for nonnegative params */

    if (DH.boundcheck && !params_in_bounds(par, ncm, DH.k)) {
	return GretlError.data;
    }

    /* compute residuals */

    tmp = 0.0;
    for (t=t0; t<=t2; t++) {
	if (t < t1) {
	    et = 0.0;
	} else {
	    et = DH.y[t];
	    for (i=0; i<ncm; i++) {
		et -= DH.X[i][t] * par[i];
	    }
	    DH.e[t] = et;
	    DH.e2[t] = et * et;
	    tmp += DH.e2[t]; 
	}
    }

    for (t=t0; t<t1; t++) {
	for (i=0; i<DH.k; i++) {
	    dedq[i][t] = 0.0;
	}
    }
	
    /* h0 and derivatives */

    if (DH.init_ == InitVar.INIT_VAR_OLS) {
	h0 = 1.0;
    } else if (DH.init_ == InitVar.INIT_VAR_RESID) {
	h0 = tmp / T;
    } else if (DH.init_ == InitVar.INIT_VAR_THEO) {
	tmp = 1.0;
	for (i=ncm+1; i<DH.k; i++) {
	    tmp -= par[i];
	}
	h0 = u_var = par[ncm] / tmp;
    }


    for (t=t0; t<t1; t++) {
	DH.h[t] = DH.e2[t] = h0;
    }

    if (DH.init_ == InitVar.INIT_VAR_OLS) {
	for (t=t0; t<t1; t++) {
	    for (i=0; i<DH.k; i++) {
		dhdq[i][t] = 0.0;
	    }
	}
    } else if (DH.init_ == InitVar.INIT_VAR_RESID) {
	for (i=0; i<ncm; i++) {
	    dh0 = 0.0;
	    for (t=t1; t<=t2; t++) {
		dh0 -= DH.e[t] * DH.X[i][t];
	    }
	    for (t=t0; t<t1; t++) {
		dhdq[i][t] = dh0 * 2.0 / T;
	    }
	}
	for (t=t0; t<t1; t++) {
	    for (i=ncm; i<DH.k; i++) {
		dhdq[i][t] = 0.0;
	    }
	}
    } else if (DH.init_ == InitVar.INIT_VAR_THEO) {
	for (t=t0; t<t1; t++) {
	    for (i=0; i<ncm; i++) {
		dhdq[i][t] = 0.0;
	    }
	}
	dh0 = u_var / par[ncm];
	for (t=t0; t<t1; t++) {
	    dhdq[ncm][t] = dh0;
	}
	dh0 *= u_var;
	for (t=t0; t<t1; t++) {
	    for (i=ncm+1; i<DH.k; i++) {
		dhdq[i][t] = dh0;
	    }
	}
    }

    /* in-sample loop */

    for (t=t1; t<=t2; t++) {
	ht = par[ncm]; /* omega */

	for (i=1; i<=q; i++) {
	    ht += DH.e2[t-i] * par[ncm+i];
	}

	for (i=1; i<=p; i++) {
	    ht += DH.h[t-i] * par[ncm+i+q];
	}
	
	DH.h[t] = ht;
	    
	/* regressors in mean equation */
	for (i=0; i<ncm; i++) {
	    dedq[i][t] = -(DH.X[i][t]);
	    k = ncm;
	    dhdq[i][t] = 0.0;
	    for (j=1; j<=q; j++) {
		if (t - q < t1 && DH.init_ == InitVar.INIT_VAR_RESID) { 
		    /* add INIT_THEO here */
		    dhdq[i][t] += par[k+j] * dhdq[i][t1-1];
		} else {	
		    dhdq[i][t] += 2.0 * par[k+j] * DH.e[t-j] * dedq[i][t-j];
		}
	    }
	}
	    
	/* garch params: omega */
	dedq[ncm][t] = 0.0;
	dhdq[ncm][t] = 1.0;
	if (t - p < t1 && DH.init_ == InitVar.INIT_VAR_THEO) {
	    for (i=1; i<=p; i++) {
		dhdq[ncm][t] += par[ncm+i] * dhdq[ncm][t1-1];
	    }
	}
	    
	/* garch params: alphas */
	k = ncm + 1;
	for (i=1; i<=q; i++) {
	    dedq[k][t] = 0.0;
	    dhdq[k][t] = DH.e2[t-i];
	    if (t - p < t1 && DH.init_ == InitVar.INIT_VAR_THEO) {
		for (j=0; j<p; j++) {
		    dhdq[k][t] += par[k+j] * dhdq[k][t1-1];
		}
	    }
	    k++;
	}
	    
	/* garch params: betas */
	k = ncm + q + 1;
	for (i=1; i<=p; i++) {
	    dedq[k][t] = 0.0;
	    dhdq[k][t] = DH.h[t-i];
	    if (t - p < t1 && DH.init_ == InitVar.INIT_VAR_THEO) {
		for (j=0; j<p; j++) {
		    dhdq[k][t] += par[k+j-p] * dhdq[k][t1-1];
		}
	    }
	    k++;
	}
	    
	/* "real" recursive part */
	for (i=0; i<DH.k; i++) {
	    k = ncm + q + 1;
	    for (j=1; j<=p; j++) {
		dhdq[i][t] += par[k++] * dhdq[i][t-j];
	    }
	}
    }	
    

    return 0;
} 

extern (C) double garch_loglik(const(double) * theta, void * ptr)
{
    GarchContainer *DH = cast(GarchContainer *) ptr;
    double ll = NADBL;
    int err;

    err = garch_etht(theta, DH);
    if (!err) {
	ll = normal_ll(DH);
    }

    return ll;
}

private int score_fill_matrices(const(double) * theta, GarchContainer * DH)
{
    int i;
    int t;
    int err;

    err = garch_etht(theta, DH);
    if (err) {
	return err;
    }

    normal_score(DH);

    for (t=DH.t1; t<=DH.t2; t++) {
	for (i=0; i<DH.k; i++) {
	    DH.G[i][t] = DH.score_e[i][t] * DH.blockglue[0][t] + 
		DH.score_h[i][t] * DH.blockglue[1][t];
	}
    }

    return err;
}

extern (C) int garch_score(double * theta, double * s, int npar, BFGS_CRIT_FUNC ll, void * ptr)
{
    GarchContainer *DH = cast(GarchContainer *) ptr;
    double tmp;
    int t;
    int i;
    int err;

    err = score_fill_matrices(theta, DH);
    if (err) {
	return err;
    }

    for (i=0; i<npar; i++) {
	tmp = 0.0;
	for (t=DH.t1;t<=DH.t2; t++) {
	    tmp += DH.G[i][t];
	}
	s[i] = tmp;
    }

    return err;
}

private GretlMatrixRaw * garch_iinfo(GarchContainer * DH, int * err)
{
    GretlMatrixRaw * info;
    double **tmp_info;
    double tmpi;
    double tmpj;
    double tmpx1;
    double tmpx2;
    double x;
    int ncm = DH.ncm;
    int i;
    int j;
    int t;

    info = gretl_matrix_alloc(DH.k, DH.k);
    if (info == null) {
	*err = GretlError.alloc;
	return null;
    }

    tmp_info = doubles_array_new(DH.k, DH.k);
    if (tmp_info == null) {
	*err = GretlError.alloc;
	gretl_matrix_free(info);
	return null;
    }

    for (i=0; i<DH.k; i++) {
	for (j=0; j<=i; j++) {
	    tmp_info[i][j] = 0.0;
	}
    }

    for (t=DH.t1; t<=DH.t2; t++) {
	for (i=0; i<ncm; i++) {
	    tmpi = DH.score_h[i][t] / DH.h[t];
	    tmpx1 = 2.0 * DH.X[i][t];
	    tmpx1 /= DH.h[t];
	    for (j=0; j<=i; j++) {
		tmpj = DH.score_h[j][t] / DH.h[t];
		tmpx2 = DH.X[j][t];
		x = tmpx1 * tmpx2 + tmpi * tmpj;
		tmp_info[i][j] += x;
	    }
	}

	for (i=ncm; i<DH.k; i++) {
	    tmpi = DH.score_h[i][t] / DH.h[t];
	    for (j=ncm; j<=i; j++) {
		tmpj = DH.score_h[j][t] / DH.h[t];
		x = tmpi * tmpj;
		tmp_info[i][j] += x;
	    }
	}
    }

    for (i=0; i<DH.k; i++) {
	for (j=0; j<=i; j++) {
	    gretl_matrix_set(info, i, j, 0.5 * tmp_info[i][j]);
	    if (j < i) {
		gretl_matrix_set(info, j, i, 0.5 * tmp_info[i][j]);
	    }
	}
    }

    doubles_array_free(tmp_info, DH.k);


    gretl_invert_symmetric_matrix(info);


    return info;
}

private GretlMatrixRaw * garch_opg(GarchContainer * DH, int * err)
{
    GretlMatrixRaw * GG;
    double **tmp_GG;
    double tmpi;
    double x;
    int t;
    int i;
    int j;

    GG = gretl_matrix_alloc(DH.k, DH.k);
    if (GG == null) {
	*err = GretlError.alloc;
	return null;
    }

    tmp_GG = doubles_array_new(DH.k, DH.k);
    if (tmp_GG == null) {
	*err = GretlError.alloc;
	gretl_matrix_free(GG);
	return null;
    }

    for (i=0; i<DH.k; i++) {
	for (j=0; j<DH.k; j++) {
	    tmp_GG[i][j] = 0.0;
	}
    }

    for (t=DH.t1; t<=DH.t2; t++) {
	for (i=0; i<DH.k; i++) {
	    tmpi = DH.G[i][t];
	    for (j=0; j<=i; j++) {
		x = tmpi * DH.G[j][t]; 
		tmp_GG[i][j] += x;
	    }
	}
    }

    for (i=0; i<DH.k; i++) {
	for (j=0; j<=i; j++) {
	    gretl_matrix_set(GG, i, j, tmp_GG[i][j]);
	    if (j < i) {
		gretl_matrix_set(GG, j, i, tmp_GG[i][j]);
	    }
	}
    }

    doubles_array_free(tmp_GG, DH.k);


    return GG;
}

private GretlMatrixRaw * garch_ihess(GarchContainer * DH, double * theta, int * err)
{
    return garch_analytical_hessian(DH.y, DH.X, DH.t1, DH.t2, DH.nobs, 
				    DH.ncm, DH.p, DH.q, theta, DH.e, 
				    DH.e2, DH.h, DH.scale, err);
}

private int garch_covariance_matrix(int vopt, double * theta, GarchContainer * DH, GretlMatrixRaw * V)
{
    GretlMatrixRaw * GG = null;
    GretlMatrixRaw * iinfo = null;
    GretlMatrixRaw * invhess = null;
    int err = 0;

    if (vopt == MLVCVType.op || vopt == MLVCVType.qml || vopt == MLVCVType.bw) { 
	/* GG' needed */
	GG = garch_opg(DH, &err);
    }

    if (vopt == MLVCVType.im || vopt == MLVCVType.bw) { 
	/* information matrix needed */
	iinfo = garch_iinfo(DH, &err);
    }

    if (vopt == MLVCVType.qml || vopt == MLVCVType.hessian) { 
	/* Hessian matrix needed */
	invhess = garch_ihess(DH, theta, &err);
    }

    if (err) {
	goto bailout;
    }

    switch (vopt) {
    case MLVCVType.hessian:
	gretl_matrix_copy_values(V, invhess);
	break;
    case MLVCVType.im:
	gretl_matrix_copy_values(V, iinfo);
	break;
    case MLVCVType.op:
	gretl_matrix_copy_values(V, GG);
	err = gretl_invert_symmetric_matrix(V);
	break;
    case MLVCVType.bw:
	gretl_matrix_qform(iinfo, GretlMatrixMod.none, GG,
			   V, GretlMatrixMod.none);
	break;
    case MLVCVType.qml:
	gretl_matrix_qform(invhess, GretlMatrixMod.none, GG,
			   V, GretlMatrixMod.none);
	break;
    default: 
	break;
    }


 bailout:

    gretl_matrix_free(GG);
    gretl_matrix_free(iinfo);
    gretl_matrix_free(invhess);

    return err;
}


/*
   Parameters to garch_estimate_mod()

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
   fncount: 0 on input, number of function evaluations on output
   fncount: 0 on input, number of gradient evaluations on output
   vopt:  code indicating which version of the covariance
          matrix to compute in V
   prn:   print handle for info on iterations etc.

*/

extern (C) int garch_estimate_mod(const(double) * y, const(double) ** X, int t1, int t2, int nobs, int nc, int p, int q, double * theta, GretlMatrixRaw * V, double * e, double * e2, double * h, double scale, double * pll, int * fncount, int * grcount, int vopt, Prn * prn)
{
    GarchContainer *DH;
    Gretlopt opt;
    double toler;
    int npar;
    int maxit;
    int use_newton = 0;
    int err = 0;

    DH = garch_container_new(y, X, t1, t2, nobs, nc, p, q, 
			     InitVar.INIT_VAR_RESID, e, e2, h, 
			     scale);

    if (DH == null) {
	return GretlError.alloc;
    }

    npar = nc + 1 + p + q; 

    if (libset_get_int(SetKey.gretl_optim) == OptimCode.newton) {
	use_newton = 1;
    }

    BFGS_defaults(&maxit, &toler, GretlCmdIndex.garch);
    opt = (prn != null)? (Gretlopt.OPT_U | Gretlopt.OPT_V) : Gretlopt.OPT_U;

    if (use_newton) {
	double crittol = 1.0e-7;
	double gradtol = 1.0e-7;

	maxit = 100;
	err = newton_raphson_max(theta, npar, maxit, 
				 crittol, gradtol, fncount, 
				 CriterionType.loglik, &garch_loglik, 
				 &garch_score, null, DH,
				 opt, prn);
    } else {
	err = BFGS_max(theta, npar, maxit, toler, 
		       fncount, grcount, &garch_loglik, CriterionType.loglik,
		       &garch_score, DH, null, opt, prn);
    }


    if (err) {
	*pll = NADBL;
    } else {
	*pll = garch_loglik(theta, DH) - (t2 - t1 + 1) * log(scale);
    }
    

    if (!err) {
	err = garch_covariance_matrix(vopt, theta, DH, V);
    }

    garch_container_destroy(DH);

    return err;
}
