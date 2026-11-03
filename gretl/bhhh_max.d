/* Port of lib/src/bhhh_max.c: maximum likelihood by the BHHH method (iteration of the outer-product-of-the-gradient
   regression with a line search; the header bhhh_max.h is bhhh_maxh.d).

   Differences from the C:
   - the BHHH_DEBUG code (compiled out in the C) is not ported;
   - the C99 `restrict` qualifiers are dropped. */

module gretl.bhhh_max;


import core.stdc.stdio : fprintf, stderr;
import core.stdc.stdlib : free, malloc;

import gretl.bhhh_maxh : BHHH_FUNC;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pputs;
import gretl.libgretlh;
import gretl.libseth : SetKey;
import gretl.libset : libset_get_int;
import gretl.missingh : NADBL, na;
import gretl.nlsh : CriterionType;
import gretl.printout : print_iter_info;

/* unported dependency (nls.c): per-file declaration from the 2026c headers */
extern (C) {
}

/* just experimental at this point: call @func to compute the
   per-observation contributions to the log-likelihood in @L,
   and use that info to compute @G numerically.
*/

double bhhh_numeric(double * b, int k,
                    GretlMatrixRaw * L,
                    GretlMatrixRaw * G,
                    BHHH_FUNC func, void * data,
                    int * err)
{
    double * Lval = L.val;
    const double h = 1.0e-8;
    double gti, ret;
    int i, t, T = G.rows;

    /* loglik at current evaluation point */
    ret = func(b, L, data, 0, err);

    for (i=0; i<k && !*err; i++) {
        double ll, bi0 = b[i];

        b[i] = bi0 - h;
        ll = func(b, L, data, 0, err);

        if (na(ll)) {
            ret = NADBL;
            break;
        }

        for (t=0; t<T; t++) {
            gretl_matrix_set(G, t, i, Lval[t]);
        }

        b[i] = bi0 + h;
        ll = func(b, L, data, 0, err);

        b[i] = bi0;

        if (na(ll)) {
            ret = NADBL;
            break;
        }

        for (t=0; t<T; t++) {
            gti = gretl_matrix_get(G, t, i);
            gti = (Lval[t] - gti) / (2.0 * h);
            gretl_matrix_set(G, t, i, gti);
        }
    }

    return ret;
}

/* in case we want to compute the score numerically */

GretlMatrixRaw * make_score_matrix(GretlMatrixRaw * L, int k, int * err)
{
    int T = gretl_matrix_rows(L);
    GretlMatrixRaw * G = null;

    if (T == 0) {
        *err = GretlError.data;
    } else if (T <= k) {
        *err = GretlError.df;
    } else {
        G = gretl_zero_matrix_new(T, k);
        if (G == null) {
            *err = GretlError.alloc;
        }
    }

    return G;
}

/**
 * bhhh_max:
 * @theta: array of adjustable coefficients.
 * @k: number of elements in @theta.
 * @M: T x @k matrix to hold the gradient.
 * @callback: pointer to function for calculating log-likelihood and
 * score.
 * @toler: tolerance for convergence.
 * @fncount: location to receive count of function evaluations.
 * @grcount: location to receive count of gradient evaluations.
 * @data: pointer to be passed to @loglik.
 * @V: matrix in which to store covariance, or %NULL.
 * @opt: can include %OPT_V for verbose output.
 * @prn: printing struct for iteration info (or %NULL).
 *
 * Maximize likelihood using the BHHH method, implemented via
 * iteration of the Outer Product of the Gradient (OPG) regression
 * with line search.
 *
 * @callback is called to calculate the loglikelihood for the model
 * in question.  The parameters passed to this function are:
 * (1) the current array of estimated coefficients; (2) @G;
 * (3) the @data pointer; (4) an integer indicator that is 1 if
 * the gradient should be calculated in @G, 0 if only the
 * loglikelihood is needed; and (5) an int pointer to receive
 * an error code. The return value from @callback should be the
 * loglikelihood (or #NADBL on error).
 *
 * For an example of the use of such a function, see arma.c in the
 * %plugin directory of the gretl source.
 *
 * Returns: 0 on successful completion, non-zero error code otherwise.
 */

int bhhh_max(double * theta, int k,
             GretlMatrixRaw * M,
             BHHH_FUNC callback,
             double toler,
             int * fncount, int * grcount,
             void * data,
             GretlMatrixRaw * V,
             Gretlopt opt,
             Prn * prn)
{
    GretlMatrixRaw * c = null;
    GretlMatrixRaw * g = null;
    GretlMatrixRaw * G = null;
    double * delta = null;
    double * ctemp = null;
    double * grad;
    int itermax, iter = 0;
    int fcount = 0;
    int gcount = 0;
    double minstep = 1.0e-06;
    double crit = 1.0;
    double stepsize = 0.25;
    double shrink = 0.5;
    double ll2, ll = 0.0;
    int numeric = 0;
    int i, T, err = 0;

    if (opt & Gretlopt.OPT_N) {
        /* numerical score (not ready) */
        G = make_score_matrix(M, k, &err);
        if (err) {
            return err;
        }
        numeric = 1;
    } else {
        /* analytical score */
        if (gretl_matrix_cols(M) != k) {
            return GretlError.nonconf;
        } else {
            G = M;
        }
    }

    T = gretl_matrix_rows(G);
    c = gretl_unit_matrix_new(T, 1);
    g = gretl_column_vector_alloc(k);
    delta = cast(double *) malloc(k * (*delta).sizeof);
    ctemp = cast(double *) malloc(k * (*ctemp).sizeof);

    if (c == null || g == null || delta == null || ctemp == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    if (opt & Gretlopt.OPT_I) {
        /* using BHHH for initialization of exact ML */
        minstep = 1.0e-03;
        shrink = 0.4;
    }

    itermax = libset_get_int(SetKey.bhhh_maxiter);
    grad = g.val;

    while (crit > toler && iter++ < itermax) {
        /* compute loglikelihood and score */
        if (numeric) {
            ll = bhhh_numeric(theta, k, M, G, callback, data, &err);
        } else {
            ll = callback(theta, G, data, 1, &err);
        }
        fcount++;
        gcount++;

        if (err) {
            pputs(prn, "Error calculating log-likelihood\n");
            break;
        }

        /* BHHH via OPG regression */
        err = gretl_matrix_ols(c, G, g, null, null, null);
        if (err) {
            fprintf(stderr, "BHHH OLS error code = %d\n".ptr, err);
            break;
        }
        for (i=0; i<k; i++) {
            delta[i] = grad[i] * stepsize;
            ctemp[i] = theta[i] + delta[i];
        }

        /* see if we've gone up...  (0 means don't compute score) */
        ll2 = callback(ctemp, G, data, 0, &err);
        fcount++;
        while (err || ll2 < ll) {
            /* ... or if not, halve the steplength */
            stepsize *= shrink;
            if (stepsize < minstep) {
                fprintf(stderr, "BHHH: hit minimum step size %g\n".ptr, minstep);
                err = GretlError.noconv;
                break;
            }
            for (i=0; i<k; i++) {
                delta[i] *= shrink;
                ctemp[i] = theta[i] + delta[i];
            }
            ll2 = callback(ctemp, G, data, 0, &err);
            fcount++;
        }

        if (err) break;

        /* print iteration info, if wanted */
        if (opt & Gretlopt.OPT_V) {
            print_iter_info(iter, ll, CriterionType.loglik, k, theta, grad,
                            stepsize, prn);
        }

        /* actually update parameter estimates */
        for (i=0; i<k; i++) {
            theta[i] = ctemp[i];
        }
        /* double the steplength? (was < 4.0 below) */
        if (stepsize < 1.0) {
            stepsize *= 2.0;
        }
        crit = ll2 - ll;
    }

    *fncount = fcount;
    *grcount = gcount;

    if ((opt & Gretlopt.OPT_V) && !(opt & Gretlopt.OPT_I)) {
        print_iter_info(-1, ll, CriterionType.loglik, k, theta, grad,
                        stepsize, prn);
    }

    if (crit > toler && !err) {
        err = GretlError.noconv;
    }

    if (err) {
        fprintf(stderr, "bhhh_max: iters = %d, crit = %g, tol = %g, err = %d\n".ptr,
                iter, crit, toler, err);
    } else if (V != null) {
        /* run OPG once more, to get VCV */
        double s2 = 0.0;

        err = gretl_matrix_ols(c, G, g, V, null, &s2);
        if (!err) {
            /* undo the default df correction */
            gretl_matrix_multiply_by_scalar(V, (T-k) / cast(double) T);
        }
    }

 bailout:

    gretl_matrix_free(c);
    gretl_matrix_free(g);
    if (G != M) {
        gretl_matrix_free(G);
    }
    free(delta);
    free(ctemp);

    return err;
}
