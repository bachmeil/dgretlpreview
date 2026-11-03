/* Port of lib/src/vartest.c: tests on estimated VARs and VECMs: multivariate normality, autocorrelation and ARCH
   tests, lag-length criteria, omission tests, Wald tests and the helpers shared with var.c. The header vartest.h is
   vartesth.d. */

module gretl.vartest;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno, ERANGE;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fputs, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, realloc;
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
import gretl.vartesth : N_IVALS, VAR_S_DFCORR;

private alias XUC = const(xmlChar) *;

import gretl.describe : ljung_box, multivariate_normality_test;
import gretl.estimate : arch_test;
import gretl.genparseh : ConstId;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_list : gretl_list_const_pos, gretl_list_copy, gretl_list_diff_new, gretl_list_omit, in_gretl_list;
import gretl.gretl_model : gretl_model_set_hac_vcv_info, gretl_model_set_vcv_info, gretl_model_write_vcv;
import gretl.gretl_modelh : VCVInfo, VCVMajorType;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : get_last_pvalue, get_last_test_statistic, record_matrix_test_result, record_test_result;
import gretl.libset : libset_get_int;
import gretl.libseth : SetKey;
import gretl.options : get_optval_int;
import gretl.printout : bufspace;
import gretl.pvalues : chisq_cdf_comp, snedecor_cdf_comp;
import gretl.qr_estimate : HAC_XOX;
import gretl.strutils : gretl_strdup, strings_array_new;
import gretl.var : VAR_fill_X, VAR_list_composite, gretl_VAR, gretl_VAR_get_residual_matrix, gretl_VAR_ldet;
import gretl.varprint : gretl_VAR_print;
import gretl.compare : print_add_omit_null;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}
int gretl_VAR_normality_test(GretlVar * var, Gretlopt opt, Prn * prn)
{
    int err = 0;

    if (var.E == null || var.S == null) {
    err = 1;
    } else {
    err = multivariate_normality_test(var.E, var.S,
                      opt, prn);
    }

    return err;
}

/* Multivariate version of Breusch-Godfrey test, as per
   H. Lutkepohl, New Introduction to Multiple Time Series
   Analysis (Springer, 2005) section 4.4, pp. 173-4.
*/

int multivariate_autocorr_test(GretlVar * var, int H, int autoH, Gretlopt opt, Prn * prn)
{
    const(GretlMatrixRaw) * U;
    char[32] Fspec;
    GretlMatrixRaw * tests, pvals;
    GretlMatrixRaw * B = null;
    GretlMatrixRaw * E = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * S = null;
    double * targ;
    double ulag, s, FRao;
    double detUU, detEE;
    double N, dfn, dfd;
    int quiet = (opt & Gretlopt.OPT_Q);
    int h, h2, K = var.neqns;
    int i, t, nx, K2 = K * K;
    int g, lagcol, hw = 0;
    int p = var.order;
    int T = var.T;
    int err = 0;

    /* we need var->X and var->S below */
    if (var == null || var.X == null || var.S == null) {
    return GretlError.data;
    }

    /* and also the matrix of residuals */
    U = gretl_VAR_get_residual_matrix(var);
    if (U == null) {
    return GretlError.data;
    }

    /* how many regressors are we going to have at max? */
    g = var.ncoeff;
    if (!var.ifc) {
    /* we'll have to add an intercept */
    g++;
    }
    nx = g + K * H;

    if (nx >= T) {
    /* not enough data to do this? */
    if (autoH) {
        H = cast(int) floor((T - 1 - g) / cast(double) K);
        if (H <= 0) {
        return GretlError.toofew;
        } else {
        nx = g + K * H;
        }
    } else {
        return GretlError.toofew;
    }
    }

    tests = gretl_column_vector_alloc(H);
    pvals = gretl_column_vector_alloc(H);

    if (tests == null || pvals == null) {
    return GretlError.alloc;
    }

    /* calculate determinant of cross-equation Sigma */
    S = gretl_matrix_copy(var.S);
    if (S == null) {
    err = GretlError.alloc;
    } else {
    detUU = gretl_matrix_determinant(S, &err);
    }

    if (err) {
    goto bailout;
    }

    /* allocate E, B, X to max size */
    E = gretl_matrix_alloc(T, K);
    B = gretl_matrix_alloc(nx, K);
    X = gretl_matrix_alloc(T, nx);

    if (E == null || B == null || X == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    if (var.ifc) {
    targ = X.val;
    } else {
    /* insert constant */
    for (t=0; t<T; t++) {
        X.val[t] = 1.0;
    }
    targ = X.val + T;
    }

    /* copy the original regressors into X */
    memcpy(targ, var.X.val, var.ncoeff * T * double.sizeof);

    if (!quiet) {
    pprintf(prn, "%s %d\n".ptr, _("Test for autocorrelation of order up to").ptr, H);
    hw = cast(int) ceil(log10(H));
    pputc(prn, '\n');
    bufspace(9 + hw, prn);
    pputs(prn, "Rao F   Approx dist.  p-value\n");
    }

    for (h=1; h<=H && !err; h++) {
    nx = g + K * h;
    gretl_matrix_reuse(B, nx, K);
    gretl_matrix_reuse(X, T, nx);
    h2 = h * h;
    /* add next lag of U to the X matrix */
    for (i=0; i<K; i++) {
        lagcol = g + K*(h-1) + i;
        for (t=0; t<T; t++) {
        ulag = (t - h < 0)? 0.0 : gretl_matrix_get(U, t-h, i);
        gretl_matrix_set(X, t, lagcol, ulag);
        }
    }
    err = gretl_matrix_multi_SVD_ols(U, X, B, E, null);
    if (!err) {
        err = gretl_matrix_multiply_mod(E, GretlMatrixMod.transpose,
                        E, GretlMatrixMod.none,
                        S, GretlMatrixMod.none);
    }
    if (err) {
        break;
    }
    gretl_matrix_divide_by_scalar(S, T);
    /* calculate the test statistic */
    s = sqrt((pow(K, 4) * h2 - 4.0) / (K2 + K2 * h2 - 5.0));
    N = T - K*p - 1 - K*h - (K - K*h + 1) / 2.0;
    dfn = K2 * h;
    dfd = floor(N*s - 0.5*(K2 * h) + 1);
    detEE = gretl_matrix_determinant(S, &err);
    if (!err) {
        FRao = pow(detUU / detEE, 1/s) - 1.0;
        FRao *= dfd / dfn;
        tests.val[h-1] = FRao;
        pvals.val[h-1] = snedecor_cdf_comp(dfn, dfd, FRao);
        if (!quiet) {
        sprintf(Fspec.ptr, "F(%g, %g)".ptr, dfn, dfd);
        pprintf(prn, "lag %*d %9.3f    %-13s %5.4f\n".ptr, hw, h,
            FRao, Fspec.ptr, pvals.val[h-1]);
        }
    }
    }

    if (!quiet) {
    pputc(prn, '\n');
    }

 bailout:

    gretl_matrix_free(E);
    gretl_matrix_free(B);
    gretl_matrix_free(X);
    gretl_matrix_free(S);

    if (!err) {
    record_matrix_test_result(tests, pvals);
    } else {
    gretl_matrix_free(tests);
    gretl_matrix_free(pvals);
    }

    return err;
}

int univariate_autocorr_test(GretlVar * var, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    Model * pmod;
    GretlMatrixRaw * tests, pvals;
    double lb;
    int i, quiet = (opt & Gretlopt.OPT_Q);
    int err = 0;

    tests = gretl_column_vector_alloc(var.neqns);
    pvals = gretl_column_vector_alloc(var.neqns);

    if (tests == null || pvals == null) {
    err = GretlError.alloc;
    }

    if (!quiet) {
    pprintf(prn, "%s %d\n\n".ptr, _("Test for autocorrelation of order").ptr,
        order);
    }

    for (i=0; i<var.neqns && !err; i++) {
    pmod = var.models[i];
    if (!quiet) {
        pprintf(prn, "%s %d:\n".ptr, _("Equation").ptr, i + 1);
    }
    lb = ljung_box(order, pmod.t1, pmod.t2, pmod.uhat, &err);
    if (!err) {
        tests.val[i] = lb;
        pvals.val[i] = chisq_cdf_comp(order, lb);
        if (!(opt & Gretlopt.OPT_Q)) {
        pprintf(prn, "Ljung-Box Q' = %g %s = P(%s(%d) > %g) = %.3g\n".ptr,
            lb, _("with p-value").ptr, _("Chi-square").ptr, order,
            lb, pvals.val[i]);
        pputc(prn, '\n');
        }
    }
    }

    if (!err) {
    record_matrix_test_result(tests, pvals);
    } else {
    gretl_matrix_free(tests);
    gretl_matrix_free(pvals);
    }

    return err;
}

int gretl_VAR_autocorrelation_test(GretlVar * var, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int h = order;
    int err;

    if (order == 0) {
    h = dset.pd;
    }

    if (opt & Gretlopt.OPT_U) {
    err = univariate_autocorr_test(var, h, dset, opt, prn);
    } else {
    int autoH = (order == 0);

    err = multivariate_autocorr_test(var, h, autoH, opt, prn);
    }

    return err;
}

/* write the vech of @src into row @t of @targ */

void vech_into_row(GretlMatrixRaw * targ, int t, const(GretlMatrixRaw) * src)
{
    double vij;
    int j, i, k = 0;

    for (j=0; j<src.cols; j++) {
    for (i=0; i<=j; i++) {
        vij = gretl_matrix_get(src, i, j);
        gretl_matrix_set(targ, t, k++, vij);
    }
    }
}

/* write from @src into @targ, with a lag of @h,
   starting at column @inicol and implicitly
   discarding rows of @src that would go off the
   end of @targ
*/

void inscribe_lag(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, int h, int inicol)
{
    double xij;
    int j, t, k = inicol;

    for (j=0; j<src.cols; j++) {
    for (t=h; t<targ.rows; t++) {
        xij = gretl_matrix_get(src, t-h, j);
        gretl_matrix_set(targ, t, k, xij);
    }
    k++;
    }
}

int multivariate_arch_test(GretlVar * var, int H, int autoH, Gretlopt opt, Prn * prn)
{
    const(GretlMatrixRaw) * U;
    GretlMatrixBlock * Bk;
    GretlMatrixRaw * tests, pvals;
    GretlMatrixRaw * vU, B, X;
    GretlMatrixRaw * E, S, SS;
    GretlMatrixRaw * iC0, ut, uu;
    double df, LM, K24;
    int quiet = (opt & Gretlopt.OPT_Q);
    int h, K = var.neqns;
    int i, t, nx, KK1;
    int vdim, inicol;
    int hw = 0;
    int T = var.T;
    int err = 0;

    if (var == null) {
    return GretlError.data;
    }

    U = gretl_VAR_get_residual_matrix(var);
    if (U == null) {
    return GretlError.data;
    }

    /* dimension of vech of variance */
    vdim = (K * (K + 1)) / 2;

    /* how many regressors are we going to have at max? */
    nx = 1 + vdim * H;

    if (nx >= T) {
    /* not enough data to do this? */
    if (autoH) {
        H = cast(int) floor((T - 1) / cast(double) vdim);
        if (H <= 0) {
        return GretlError.toofew;
        } else {
        nx = 1 + vdim * H;
        }
    } else {
        return GretlError.toofew;
    }
    }

    tests = gretl_column_vector_alloc(H);
    pvals = gretl_column_vector_alloc(H);

    if (tests == null || pvals == null) {
    return GretlError.alloc;
    }

    Bk = gretl_matrix_block_new(&vU, T, vdim,
                &B, nx, vdim,
                &E, T, vdim,
                &X, T, nx,
                &ut, 1, K,
                &uu, K, K,
                &S, vdim, vdim,
                &SS, vdim, vdim,
                &iC0, vdim, vdim,
                null);
    if (Bk == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    /* construct LHS (vU) */
    for (t=0; t<T; t++) {
    for (i=0; i<K; i++) {
        ut.val[i] = gretl_matrix_get(U, t, i);
    }
    gretl_matrix_multiply_mod(ut, GretlMatrixMod.transpose,
                  ut, GretlMatrixMod.none,
                  uu, GretlMatrixMod.none);
    vech_into_row(vU, t, uu);
    }

    /* initial setup of RHS (X) */
    gretl_matrix_zero(X);
    for (t=0; t<T; t++) {
    X.val[t] = 1.0;
    }

    /* useful constants */
    KK1 = K * (K+1);
    K24 = KK1 * KK1 / 4.0;

    if (!quiet) {
    pprintf(prn, "%s %d\n".ptr, _("Test for ARCH of order up to").ptr, H);
    hw = cast(int) ceil(log10(H));
    pputc(prn, '\n');
    bufspace(10 + hw, prn);
    pputs(prn, "LM       df     p-value\n");
    }

    inicol = 1;

    for (h=0; h<=H && !err; h++) {
    nx = 1 + vdim * h;
    gretl_matrix_reuse(B, nx, vdim);
    gretl_matrix_reuse(X, T, nx);
    if (h > 0) {
        /* add next lag of vU to the X matrix */
        inscribe_lag(X, vU, h, inicol);
        inicol += vdim;
    }
    /* run aux regression */
    err = gretl_matrix_multi_SVD_ols(vU, X, B, E, null);
    if (!err) {
        err = gretl_matrix_multiply_mod(E, GretlMatrixMod.transpose,
                        E, GretlMatrixMod.none,
                        S, GretlMatrixMod.none);
    }
    if (err) {
        break;
    }
    gretl_matrix_divide_by_scalar(S, T);
    if (h == 0) {
        /* record baseline: inv(Cov0) */
        gretl_matrix_copy_values(iC0, S);
        err = gretl_invert_symmetric_matrix(iC0);
    } else {
        /* calculate the test statistic */
        gretl_matrix_multiply(S, iC0, SS);
        LM = T * (KK1 / 2.0 - gretl_matrix_trace(SS));
        df = floor(h * K24);
        tests.val[h-1] = LM;
        pvals.val[h-1] = chisq_cdf_comp(df, LM);
        if (!quiet) {
        pprintf(prn, "lag %*d %9.3f %6d %11.4f\n".ptr, hw, h,
            LM, cast(int) df, pvals.val[h-1]);
        }
    }
    }

    if (!quiet) {
    pputc(prn, '\n');
    }

 bailout:

    gretl_matrix_block_destroy(Bk);

    if (!err) {
    record_matrix_test_result(tests, pvals);
    } else {
    gretl_matrix_free(tests);
    gretl_matrix_free(pvals);
    }

    return err;
}

int univariate_arch_test(GretlVar * var, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * tests, pvals;
    int i, err = 0;

    tests = gretl_column_vector_alloc(var.neqns);
    pvals = gretl_column_vector_alloc(var.neqns);

    if (tests == null || pvals == null) {
    err = GretlError.alloc;
    } else if (!(opt & Gretlopt.OPT_I)) {
    pprintf(prn, "%s %d\n\n".ptr, _("Test for ARCH of order").ptr, order);
    }

    for (i=0; i<var.neqns && !err; i++) {
    if (!(opt & Gretlopt.OPT_I)) {
        pprintf(prn, "%s %d:\n".ptr, _("Equation").ptr, i + 1);
    }
    /* add OPT_M for multi-equation output */
    err = arch_test(var.models[i], order, dset, opt | Gretlopt.OPT_M, prn);
    if (!err) {
        tests.val[i] = get_last_test_statistic();
        pvals.val[i] = get_last_pvalue();
    }
    }

    if (!err) {
    record_matrix_test_result(tests, pvals);
    } else {
    gretl_matrix_free(tests);
    gretl_matrix_free(pvals);
    }

    return err;
}

int gretl_VAR_arch_test(GretlVar * var, int order, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    int h = order;
    int err;

    if (order == 0) {
    h = dset.pd;
    }

    if (opt & Gretlopt.OPT_U) {
    err = univariate_arch_test(var, h, dset, opt, prn);
    } else {
    int autoH = (order == 0);

    err = multivariate_arch_test(var, h, autoH, opt, prn);
    }

    return err;
}

void form_C0j(GretlVar * var, GretlMatrixRaw * C0j, GretlMatrixRaw * et, GretlMatrixRaw * ej, int j)
{
    int i, t;

    gretl_matrix_zero(C0j);

    for (t=j; t<var.T; t++) {
    /* load e_t and e_{t-j} */
    for (i=0; i<var.neqns; i++) {
        et.val[i] = gretl_matrix_get(var.E, t, i);
        ej.val[i] = gretl_matrix_get(var.E, t-j, i);
    }
    /* add e_t' * e_{t-j} */
    gretl_matrix_multiply_mod(et, GretlMatrixMod.transpose,
                  ej, GretlMatrixMod.none,
                  C0j, GretlMatrixMod.cumulate);
    }

    gretl_matrix_divide_by_scalar(C0j, var.T);
}

/* See S. Johansen, Likelihood-Based Inference in Cointegrated
   Vector Autoregressive Models, 1995, pp. 21-22
*/

int VAR_portmanteau_test(GretlVar * var)
{
    GretlMatrixBlock * B;
    GretlMatrixRaw * C00, C0j;
    GretlMatrixRaw * et, ej;
    GretlMatrixRaw * L, R, Tmp;
    int k, n = var.neqns;
    double trj, LB = 0.0;
    int s, j, err = 0;

    var.LB = NADBL;
    var.LBs = 0;

    /* we'll do this only for unrestricted VARs */
    if (var.ci == GretlCmdIndex.vecm && jrank(var) < var.neqns) {
    return 0;
    }

    /* any guidance on the order for this test? */
    s = var.T / 4;
    if (s > 48) s = 48;

    k = levels_order(var);
    if (s - k <= 0) {
    /* no degrees of freedom */
    return 0;
    }

    B = gretl_matrix_block_new(&C00, n, n,
                   &C0j, n, n,
                   &et,  1, n,
                   &ej,  1, n,
                   &L,   n, n,
                   &R,   n, n,
                   &Tmp, n, n,
                   null);

    if (B == null) {
    return GretlError.alloc;
    }

    form_C0j(var, C00, et, ej, 0);
    err = gretl_invert_symmetric_matrix(C00);

    for (j=1; j<=s && !err; j++) {
    form_C0j(var, C0j, et, ej, j);
    gretl_matrix_multiply(C0j, C00, L);
    gretl_matrix_multiply_mod(C0j, GretlMatrixMod.transpose,
                  C00, GretlMatrixMod.none,
                  R, GretlMatrixMod.none);
    gretl_matrix_multiply(L, R, Tmp);
    trj = gretl_matrix_trace(Tmp);
    LB += (1.0 / (var.T - j)) * trj;
    }

    if (!err) {
    LB *= var.T * (var.T + 2);
    var.LB = LB;
    var.LBs = s;
    }

    gretl_matrix_block_destroy(B);

    return err;
}

int VAR_LR_lag_test(GretlVar * var, const(GretlMatrixRaw) * E)
{
    double test_ldet;
    int err = 0;

    test_ldet = gretl_VAR_ldet(var, E, &err);

    if (!err) {
    double ll, AIC, BIC, HQC;
    int T = var.T;
    int g = var.neqns;
    int m = var.ncoeff - g;
    int k = g * m;

    var.LR = T * (test_ldet - var.ldet);

    ll = -(g * T / 2.0) * (ConstId.pi + 1) - (T / 2.0) * test_ldet;
    AIC = (-2.0 * ll + 2.0 * k) / T;
    BIC = (-2.0 * ll + k * log(T)) / T;
    HQC = (-2.0 * ll + 2.0 * k * log(log(T))) / T;
    var.Ivals[0] = AIC;
    var.Ivals[1] = BIC;
    var.Ivals[2] = HQC;
    }

    return err;
}

void gretl_VAR_print_lagsel(GretlMatrixRaw * seltab, int * best_row, Prn * prn)
{
    int maxlag, minlag, nrows;
    double x;
    int i, j;

    minlag = cast(int) gretl_matrix_get(seltab, 0, 0);
    nrows = gretl_matrix_rows(seltab);
    maxlag = nrows + minlag - 1;

    pprintf(prn, _("VAR system, maximum lag order %d").ptr, maxlag);
    pputs(prn, "\n\n");

    pputs(prn, _("The asterisks below indicate the best (that is, minimized) values\nof the respective information criteria, AIC = Akaike criterion,\nBIC = Schwarz Bayesian criterion and HQC = Hannan-Quinn criterion.").ptr);
    pputs(prn, "\n\n");

    pputs(prn, _("lags        loglik    p(LR)       AIC          BIC          HQC").ptr);
    pputs(prn, "\n\n");

    for (i=0; i<nrows; i++) {
    pprintf(prn, "%4d".ptr, minlag + i);
    x = gretl_matrix_get(seltab, i, 1);
    pprintf(prn, "%14.5f".ptr, x);
    x = gretl_matrix_get(seltab, i, 2);
    pprintf(prn, "%9.5f".ptr, x);
    for (j=0; j<N_IVALS; j++) {
        x = gretl_matrix_get(seltab, i, 3+j);
        pprintf(prn, "%12.6f".ptr, x);
        if (i == best_row[j]) {
        pputc(prn, '*');
        } else {
        pputc(prn, ' ');
        }
    }
    pputc(prn, '\n');
    }
    pputc(prn, '\n');
}

int lagsel_get_min_lag(int p, int * err)
{
    int m = get_optval_int(GretlCmdIndex.var, Gretlopt.OPT_M, err);

    if (!*err && (m < 0 || m > p)) {
    *err = GretlError.invarg;
    }

    return m;
}

int lagsel_run_VAR(GretlVar * var, const(DatasetRaw) * dset, GretlMatrixRaw * E, int cols0, int j, int use_QR, double * ll)
{
    int jxcols = cols0 + j * var.neqns;
    int n = var.neqns;
    int err = 0;

    if (jxcols == 0) {
    gretl_matrix_copy_values(E, var.Y);
    } else {
    VAR_fill_X(var, j, dset);
    gretl_matrix_reuse(var.X, var.T, jxcols);
    gretl_matrix_reuse(var.B, jxcols, n);
    if (use_QR) {
        err = gretl_matrix_QR_ols(var.Y, var.X, var.B,
                      E, null, null);
    } else {
        err = gretl_matrix_multi_ols(var.Y, var.X, var.B,
                     E, null);
    }
    }

    if (!err) {
    double ldet = gretl_VAR_ldet(var, E, &err);
    double T2 = var.T / 2.0;

    if (!err) {
        *ll = -(n * T2) * (ConstId.pi + 1) - T2 * ldet;
    }
    }

    return err;
}

/* apparatus for selecting the optimal lag length for a VAR */

int VAR_do_lagsel(GretlVar * var, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn, GretlMatrixRaw ** pm)
{
    GretlMatrixRaw * selmat = null;
    GretlMatrixRaw * E = null;
    int p = var.order;
    int r = cast(int) (p - 1);
    int T = var.T;
    int n = var.neqns;
    /* initialize the "best" at the longest lag */
    double[N_IVALS] best = [var.AIC, var.BIC, var.HQC];
    int[N_IVALS] best_row = [r, r, r];
    double[N_IVALS] crit;
    double LRtest;
    double ll, ll0 = NADBL;
    int cols0, minlag = 1;
    int nrows, ncols;
    int use_QR = 0;
    int j, m = 0;
    int err = 0;

    /* number of cols in X that are not Y lags */
    cols0 = var.ncoeff - p * n;

    if (opt & Gretlopt.OPT_M) {
    minlag = lagsel_get_min_lag(p, &err);
    }
    if (p < minlag + 1) {
    return 0;
    }

    E = gretl_matrix_alloc(T, n);
    if (E == null) {
    err = GretlError.alloc;
    } else {
    nrows = p - minlag + 1;
    ncols = 3 + N_IVALS;
    selmat = gretl_matrix_alloc(nrows, ncols);
    if (selmat == null) {
        err = GretlError.alloc;
    }
    }

    if (err) {
    goto bailout;
    }

    if (getenv("VAR_USE_QR") != null) {
    use_QR = 1;
    }

    /* prep for the first likelihood-ratio test: 0-order VAR */
    err = lagsel_run_VAR(var, dset, E, cols0, 0, use_QR, &ll0);

    for (j=minlag; j<p && !err; j++) {
    err = lagsel_run_VAR(var, dset, E, cols0, j, use_QR, &ll);
    if (!err) {
        int q = var.ncoeff - (n * (p - j));
        int c, k = n * q;

        crit[0] = (-2.0 * ll + 2.0 * k) / T;               /* AIC */
        crit[1] = (-2.0 * ll + k * log(T)) / T;            /* BIC */
        crit[2] = (-2.0 * ll + 2.0 * k * log(log(T))) / T; /* HQC */

        gretl_matrix_set(selmat, m, 0, cast(double) j);
        gretl_matrix_set(selmat, m, 1, ll);
        if (j == minlag) {
        LRtest = 2.0 * (ll - ll0);
        } else {
        LRtest = 2.0 * (ll - gretl_matrix_get(selmat, m-1, 1));
        }
        gretl_matrix_set(selmat, m, 2, chisq_cdf_comp(n * n, LRtest));
        for (c=0; c<N_IVALS; c++) {
        gretl_matrix_set(selmat, m, 3+c, crit[c]);
        if (crit[c] < best[c]) {
            best[c] = crit[c];
            best_row[c] = m;
        }
        }
        m++; /* increment table row */
    }
    }

    if (!err) {
    /* append results from the incoming VAR */
    gretl_matrix_set(selmat, m, 0, cast(double) p);
    gretl_matrix_set(selmat, m, 1, var.ll);
    LRtest = 2.0 * (var.ll - gretl_matrix_get(selmat, m - 1, 1));
    gretl_matrix_set(selmat, m, 2, chisq_cdf_comp(n * n, LRtest));
    gretl_matrix_set(selmat, m, 3, var.AIC);
    gretl_matrix_set(selmat, m, 4, var.BIC);
    gretl_matrix_set(selmat, m, 5, var.HQC);
    }

    if (!err) {
    /* finalize and ship out the matrix */
    static const(char) *[6] hds = ["lags".ptr, "loglik".ptr, "p(LR)".ptr, "AIC".ptr, "BIC".ptr, "HQC".ptr];
    char ** S = strings_array_new(ncols);
    int i;

    if (S != null) {
        for (i=0; i<ncols; i++) {
        S[i] = gretl_strdup(hds[i]);
        }
        gretl_matrix_set_colnames(selmat, S);
    }
    if (!(opt & Gretlopt.OPT_S)) {
        gretl_VAR_print_lagsel(selmat, best_row.ptr, prn);
    }
    if (pm != null) {
        *pm = selmat;
    } else {
        record_matrix_test_result(selmat, null);
    }
    selmat = null;
    }

 bailout:

    gretl_matrix_free(selmat);
    gretl_matrix_free(E);

    return err;
}

GretlMatrixRaw * VAR_get_hvec(const(GretlMatrixRaw) * X, const(GretlMatrixRaw) * XTX, int * err)
{
    GretlMatrixRaw * hvec;
    GretlMatrixRaw * Xt;
    int t, T = X.rows;
    int k = X.cols;
    double ht;

    Xt = gretl_matrix_alloc(1, k);
    hvec = gretl_column_vector_alloc(T);

    if (Xt == null || hvec == null) {
    gretl_matrix_free(Xt);
    gretl_matrix_free(hvec);
    *err = GretlError.alloc;
    return null;
    }

    for (t=0; t<T && !*err; t++) {
    gretl_matrix_get_row(X, t, Xt);
    ht = gretl_scalar_qform(Xt, XTX, err);
    gretl_vector_set(hvec, t, ht);
    }

    gretl_matrix_free(Xt);

    return hvec;
}

GretlMatrixRaw * var_hac_xox(GretlVar * var, int k, VCVInfo * vi, int * err)
{
    GretlMatrixRaw * XOX = null;
    GretlMatrixRaw * uk;
    int t;

    uk = gretl_column_vector_alloc(var.T);

    if (uk == null) {
    *err = GretlError.alloc;
    } else {
    for (t=0; t<var.T; t++) {
        uk.val[t] = gretl_matrix_get(var.E, t, k);
    }
    XOX = HAC_XOX(var.X, uk, vi, 0, err);
    gretl_matrix_free(uk);
    }

    return XOX;
}

GretlMatrixRaw * var_hc_xox(GretlVar * var, int k, int hcv, int * err)
{
    GretlMatrixRaw * XOX;
    GretlMatrixRaw * hvec = null;
    int T = var.T;
    int g = var.ncoeff;

    XOX = gretl_matrix_alloc(g, g);

    if (XOX == null) {
    *err = GretlError.alloc;
    } else if (hcv > 1) {
    hvec = VAR_get_hvec(var.X, var.XTX, err);
    }

    if (!*err) {
    /* form X' \Omega X */
    double xti, xtj, xij;
    double utk, u2, ht;
    int i, j, t;

    for (i=0; i<g; i++) {
        for (j=i; j<g; j++) {
        xij = 0.0;
        for (t=0; t<T; t++) {
            utk = gretl_matrix_get(var.E, t, k);
            u2 = utk * utk;
            if (hcv > 1) {
            ht = gretl_vector_get(hvec, t);
            u2 /= 1.0 - ht;
            if (hcv > 2) {
                u2 /= 1.0 - ht;
            }
            }
            xti = gretl_matrix_get(var.X, t, i);
            xtj = gretl_matrix_get(var.X, t, j);
            xij += u2 * xti * xtj;
        }
        if (hcv == 1) {
            xij *= cast(double) T / (T - g);
        }
        gretl_matrix_set(XOX, i, j, xij);
        if (i != j) {
            gretl_matrix_set(XOX, j, i, xij);
        }
        }
    }
    }

    gretl_matrix_free(hvec);

    return XOX;
}

/* (X'X)^{-1} * X'\Omega X * (X'X)^{-1} */

int VAR_robust_vcv(GretlVar * var, GretlMatrixRaw * V, Model * pmod, int hcv, int k)
{
    GretlMatrixRaw * XOX = null;
    VCVInfo vi = {0};
    int err = 0;

    if (var.robust == VarRobust.hac) {
    XOX = var_hac_xox(var, k, &vi, &err);
    } else {
    XOX = var_hc_xox(var, k, hcv, &err);
    }

    if (!err) {
    gretl_matrix_qform(var.XTX, GretlMatrixMod.transpose, XOX,
               V, GretlMatrixMod.none);

    if (var.robust == VarRobust.hac) {
        gretl_model_set_hac_vcv_info(pmod, vi.vmin,
                     vi.order, vi.flags,
                     vi.bw);
    } else {
        gretl_model_set_vcv_info(pmod, VCVMajorType.hc, hcv);
    }
    }

    gretl_matrix_free(XOX);

    return err;
}

/* Run the various per-equation omit tests (all lags of each var in
   turn, last lag of all vars) using the Wald method.  We also
   add the standard errors to the models here, since we have the
   per-equation covariance matrices to hand.
*/

int VAR_wald_omit_tests(GretlVar * var)
{
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * C = null;
    GretlVector * b = null;
    int hcv = libset_get_int(SetKey.hc_version);
    int p = (var.lags != null)? var.lags[0] : var.order;
    int n = var.neqns;
    int g = var.ncoeff;
    int dim = (p > n)? p : n;
    int i, j, k, m = 0;
    int any_F_err = 0;
    int err = 0;

    if (var.ifc && var.robust && g - 1 > dim) {
    /* need bigger arrays for robust overall F-test */
    dim = g - 1;
    }

    V = gretl_matrix_alloc(g, g);
    C = gretl_matrix_alloc(dim, dim);
    b = gretl_column_vector_alloc(dim);

    if (V == null || C == null || b == null) {
    return GretlError.alloc;
    }

    for (i=0; i<var.neqns && !err; i++) {
    Model * pmod = var.models[i];
    int ii, jj, jpos, ipos = var.ifc;
    double w, vij;
    int F_err = 0;

    if (var.robust) {
        err = VAR_robust_vcv(var, V, pmod, hcv, i);
    } else {
        gretl_matrix_copy_values(V, var.XTX);
        gretl_matrix_multiply_by_scalar(V, pmod.sigma * pmod.sigma);
    }

    if (err) {
        break;
    }

    /* set the per-model VCV */
    gretl_model_write_vcv(pmod, V);

    gretl_matrix_reuse(C, p, p);
    gretl_matrix_reuse(b, p, 1);

    /* exclusion of each var, all lags */
    if (p > 0) {
        for (j=0; j<var.neqns; j++) {
        double wj = NADBL;

        gretl_matrix_extract_matrix(C, V, ipos, ipos, GretlMatrixMod.none);
        for (k=0; k<p; k++) {
            b.val[k] = pmod.coeff[k + ipos];
        }
        F_err = gretl_invert_symmetric_matrix(C);
        if (!F_err) {
            wj = gretl_scalar_qform(b, C, &F_err);
        }
        if (F_err) {
            any_F_err = 1;
            var.Fvals[m++] = NADBL;
        } else {
            var.Fvals[m++] = wj / p;
        }
        ipos += p;
        }
    }

    /* exclusion of last lag, all vars? */
    if (p > 0) {
        gretl_matrix_reuse(C, n, n);
        gretl_matrix_reuse(b, n, 1);

        ipos = var.ifc + p - 1;
        for (ii=0; ii<n; ii++) {
        jpos = var.ifc + p - 1;
        for (jj=0; jj<n; jj++) {
            vij = gretl_matrix_get(V, ipos, jpos);
            gretl_matrix_set(C, ii, jj, vij);
            jpos += p;
        }
        b.val[ii] = pmod.coeff[ipos];
        ipos += p;
        }

        F_err = gretl_invert_symmetric_matrix(C);
        if (!F_err) {
        w = gretl_scalar_qform(b, C, &F_err);
        }
        if (F_err) {
        any_F_err = 1;
        var.Fvals[m++] = NADBL;
        } else {
        var.Fvals[m++] = w / n;
        }
    }

    /* exclusion of all but const? */
    if (var.ifc && var.robust) {
        gretl_matrix_reuse(C, g-1, g-1);
        gretl_matrix_reuse(b, g-1, 1);

        gretl_matrix_extract_matrix(C, V, 1, 1, GretlMatrixMod.none);
        for (k=0; k<g-1; k++) {
        b.val[k] = pmod.coeff[k+1];
        }
        F_err = gretl_invert_symmetric_matrix(C);
        if (!F_err) {
        w = gretl_scalar_qform(b, C, &F_err);
        }
        if (F_err) {
        any_F_err = 1;
        pmod.fstt = NADBL;
        } else {
        pmod.fstt = w / (g-1);
        }
    }
    }

    gretl_matrix_free(V);
    gretl_matrix_free(C);
    gretl_matrix_free(b);

    if (!err && any_F_err) {
    fprintf(stderr, "*** Warning: some F-tests could not be computed\n".ptr);
    }

    return err;
}

const(int) * gretl_VAR_get_exo_list(GretlVar * var)
{
    return (var == null)? null : var.xlist;
}

const(int) * gretl_VAR_get_endo_list(GretlVar * var)
{
    return (var == null)? null : var.ylist;
}

/* Based on the specification stored in the original VAR struct,
   plus a new_ exogenous list, constitute a full VAR list.
*/

int * build_VAR_list(GretlVar * var, int * exolist, int * err)
{
    return VAR_list_composite(var.ylist, exolist, var.rlist);
}

int gretl_VAR_real_omit_test(GretlVar * orig, GretlVar * new_, const(DatasetRaw) * dset, Gretlopt opt, Prn * prn)
{
    int * omitlist = null;
    double LR, pval;
    int df, nr = 0;
    int err = 0;

    if (orig.xlist != null) {
    if (new_.xlist == null) {
        omitlist = gretl_list_copy(orig.xlist);
    } else {
        omitlist = gretl_list_diff_new(orig.xlist, new_.xlist, 1);
    }
    if (omitlist == null) {
        return GretlError.alloc;
    }
    nr = omitlist[0];
    }
    if (opt & Gretlopt.OPT_E) {
    /* omitting seasonals */
    nr += dset.pd + 1;
    }
    if (opt & Gretlopt.OPT_T) {
    /* omitting trend */
    nr++;
    }

    LR = orig.T * (new_.ldet - orig.ldet);
    df = orig.neqns * nr;
    pval = chisq_cdf_comp(df, LR);

    record_test_result(LR, pval);

    pprintf(prn, "%s:\n".ptr, _("Test on the original VAR").ptr);

    print_add_omit_null(omitlist, dset, opt | Gretlopt.OPT_S, prn);

    pprintf(prn, "  %s: %s(%d) = %g, ".ptr, _("LR test").ptr,
        _("Chi-square").ptr, df, LR);
    pprintf(prn, _("with p-value = %g\n").ptr, pval);

    free(omitlist);

    return err;
}

int VAR_omit_check(GretlVar * var, const(int) * omitlist, Gretlopt opt)
{
    int nx = omitlist == null ? 0 : omitlist[0];
    int err = 0;

    if (nx == 0 && !(opt & (Gretlopt.OPT_T | Gretlopt.OPT_E))) {
    /* nothing to be omitted */
    err = GretlError.noomit;
    } else if ((opt & Gretlopt.OPT_T) && !(var.detflags & DetFlags.trend)) {
    /* can't have the --trend option with no auto-trend */
    err = GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_E) && !(var.detflags & DetFlags.seas)) {
    /* can't have the --seasonals option with no auto-seasonals */
    err = GretlError.badopt;
    }

    return err;
}

/**
 * gretl_VAR_omit_test:
 * @var: pointer to original VAR.
 * @omitlist: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Re-estimates a given VAR after removing the variables
 * specified in @omitlist, and reports per-equation F-tests
 * and system-wide LR tests for the null hypothesis that
 * the omitted variables have zero parameters.
 *
 * Returns: restricted VAR on success, %NULL on error.
 */

GretlVar * gretl_VAR_omit_test(GretlVar * var, const(int) * omitlist, DatasetRaw * dset, Gretlopt opt, Prn * prn, int * err)
{
    GretlVar * vnew = null;
    Gretlopt varopt = Gretlopt.OPT_NONE;
    int smpl_t1 = dset.t1;
    int smpl_t2 = dset.t2;
    int * tmplist = null;
    int * varlist = null;
    int c1 = 0;

    if (var == null) {
    *err = GretlError.data;
    return null;
    }

    *err = VAR_omit_check(var, omitlist, opt);
    if (*err) {
    return null;
    }

    if (var.ifc) {
    c1 = !gretl_list_const_pos(omitlist, 1, dset);
    }

    if (omitlist != null && omitlist[0] > 0) {
    /* create reduced exogenous vars list for test VAR */
    tmplist = gretl_list_omit(var.xlist, omitlist, 1, err);
    if (tmplist == null) {
        goto bailout;
    }
    } else if (var.xlist != null) {
    tmplist = gretl_list_copy(var.xlist);
    if (tmplist == null) {
        goto bailout;
    }
    }

    /* create full input VAR list for test VAR */
    varlist = build_VAR_list(var, tmplist, err);
    if (varlist == null) {
    goto bailout;
    }

    if ((var.detflags & DetFlags.seas) && !(opt & Gretlopt.OPT_E)) {
    varopt |= Gretlopt.OPT_D;
    }

    if ((var.detflags & DetFlags.trend) && !(opt & Gretlopt.OPT_T)) {
    varopt |= Gretlopt.OPT_T;
    }

    /* If the original VAR did not include a constant, we need to
       pass OPT_N to the test VAR to suppress the constant.
       We also need to pass OPT_N in case the constant was
       present originally but is now to be omitted.
    */
    if (var.ifc == 0 || c1 == 0) {
    varopt |= Gretlopt.OPT_N;
    }

    /* impose as sample range the estimation range of the
       original VAR */
    dset.t1 = var.t1;
    dset.t2 = var.t2;

    vnew = gretl_VAR(var.order, var.lags, varlist, dset,
             varopt, null, err);

    if (vnew != null) {
    /* do the actual test(s) */
    *err = gretl_VAR_real_omit_test(var, vnew, dset, opt, prn);
    if (!*err && prn != null) {
        gretl_VAR_print(vnew, dset, Gretlopt.OPT_NONE, prn);
    }
    }

    /* put back into dset what was there on input */
    dset.t1 = smpl_t1;
    dset.t2 = smpl_t2;

 bailout:

    free(tmplist);
    free(varlist);

    return vnew;
}

int set_wald_R_ones(GretlMatrixRaw * R, int k, int r0, int neq, int stride)
{
    if (k < 0) {
    return GretlError.data;
    } else {
    int i, rmax = r0 + neq;

    for (i=r0; i<rmax; i++) {
        gretl_matrix_set(R, i, k, 1.0);
        k += stride;
    }
    return 0;
    }
}

/**
 * gretl_VAR_wald_omit_test:
 * @var: pointer to original.
 * @omitlist: list of variables to omit from original model.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: gretl printing struct.
 *
 * Runs a Wald test for the joint omission of the exogenous
 * variables specified in @omitlist.
 *
 * Returns: 0 on successful completion, non-zero error code
 * otherwise.
 */

int gretl_VAR_wald_omit_test(GretlVar * var, const(int) * omitlist, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * B = null;
    GretlMatrixRaw * S = null;
    GretlMatrixRaw * V = null;
    GretlMatrixRaw * R = null;
    GretlMatrixRaw * RB = null;
    GretlMatrixRaw * RVR = null;
    double x, test = NADBL;
    int i, j, neq, eqk, row0;
    int pos, nr = 0, nx_omit = 0;
    int err = 0;

    if (var == null || var.B == null ||
    var.S == null || var.XTX == null) {
    return GretlError.data;
    }

    err = VAR_omit_check(var, omitlist, opt);
    if (err) {
    return err;
    }

    if (omitlist != null) {
    nx_omit = omitlist[0];
    }

    B = gretl_matrix_vectorize_new(var.B);
    if (B == null) {
    return GretlError.alloc;
    }

    neq = var.neqns;

    if (nx_omit > 0) {
    /* omitting user-specified exogenous vars */
    nr += nx_omit * neq;
    }
    if (opt & Gretlopt.OPT_T) {
    /* omit trend */
    nr += neq;
    }
    if (opt & Gretlopt.OPT_E) {
    /* omit seasonals */
    nr += (dset.pd - 1) * neq;
    }

    S = gretl_zero_matrix_new(neq, neq);
    R = gretl_zero_matrix_new(nr, B.rows);
    RB = gretl_matrix_alloc(nr, 1);
    RVR = gretl_matrix_alloc(nr, nr);

    if (S == null || R == null || RB == null || RVR == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    for (i=0; i<neq; i++) {
    x = gretl_matrix_get(var.S, i, i);
    gretl_matrix_set(S, i, i, x);
    }

    V = gretl_matrix_kronecker_product_new(S, var.XTX, &err);
    if (err) {
    goto bailout;
    }

    eqk = var.B.rows;
    pos = var.ifc + neq * var.order;
    row0 = 0;

    if (nx_omit > 0) {
    /* @omitlist holds ID numbers of exog variables */
    int vi;

    for (i=1; i<=omitlist[0]; i++) {
        vi = omitlist[i];
        if (vi == 0) {
        j = 0;
        } else {
        j = in_gretl_list(var.xlist, vi);
        j += pos - 1;
        }
        set_wald_R_ones(R, j, row0, neq, eqk);
        row0 += neq;
    }
    }
    if (opt & Gretlopt.OPT_E) {
    /* omit auto-seasonals */
    for (i=1; i<dset.pd; i++) {
        j = pos + i - 1;
        set_wald_R_ones(R, j, row0, neq, eqk);
        row0 += neq;
    }
    }
    if (opt & Gretlopt.OPT_T) {
    /* omit auto-trend (comes as last param) */
    j = eqk - 1;
    set_wald_R_ones(R, j, row0, neq, eqk);
    row0 += neq;
    }

    gretl_matrix_multiply(R, B, RB);
    gretl_matrix_qform(R, GretlMatrixMod.none, V,
               RVR, GretlMatrixMod.none);
    err = gretl_invert_symmetric_matrix(RVR);

    if (!err) {
    test = gretl_scalar_qform(RB, RVR, &err);
    if (!err) {
        double pval = chisq_cdf_comp(nr, test);

        record_test_result(test, pval);

        if (!(opt & Gretlopt.OPT_I)) {
        /* not silent */
        pprintf(prn, "%s:\n".ptr, _("Test on VAR").ptr);
        print_add_omit_null(omitlist, dset, opt | Gretlopt.OPT_S, prn);
        pprintf(prn, "  %s: %s(%d) = %g, %s %g\n\n".ptr,  _("Wald test").ptr,
            _("Chi-square").ptr, nr, test,
            _("p-value").ptr, pval);
        }
    }
    }

 bailout:

    gretl_matrix_free(B);
    gretl_matrix_free(S);
    gretl_matrix_free(V);
    gretl_matrix_free(R);
    gretl_matrix_free(RB);
    gretl_matrix_free(RVR);

    return err;
}

