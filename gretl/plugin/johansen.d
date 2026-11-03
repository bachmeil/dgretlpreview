/* Port of plugin/johansen.c: the Johansen cointegration procedure and VECM estimation: the trace-test p-values by J. Doornik's
   gamma approximation (trace_pvalue and the helpers for the other tests), the eigen-analysis, the cointegration test with its
   output, estimation of a VECM unrestricted or subject to "simple" restrictions on beta or alpha (the general restrictions are
   handled by jrestrict.d and jalpha.d), beta and alpha normalization and
   variances, and the test of restrictions on beta.

   In the C this is part of a gretl plugin loaded with get_plugin_function("johansen_estimate") and so on; here the four entry
   points trace_pvalue, johansen_coint_test, johansen_estimate, johansen_boot_round and vecm_test_restriction are extern (C) and
   plugins.d returns their addresses for those names (the types that var.d, irfboot.d and pvalues.d cast them to). The
   line-by-line translation changes nothing but the notation of the 2-d arrays (double[6][5] is 5 rows of 6) and the macros
   that became functions or enum members (the flags of VECM_estimate_full are VecmFlags, V_ALPHA and V_BETA are VBeta); the
   debugging output under JDEBUG is left out. */

module gretl.plugin.johansen;

import core.stdc.math : exp, fabs, floor, log, sqrt;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : free, qsort;
import core.stdc.string : memmove, strlen, strspn;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.gretl_array : gretl_array_get_data, gretl_array_new, gretl_array_set_data;
import gretl.gretl_list : in_gretl_list;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : gretl_inverse_compare_doubles, record_matrix_test_result;
import gretl.gretl_utilsh : floateq;
import gretl.gretl_restrict : rset_from_VECM, rset_record_LR_result, rset_VECM_acols, rset_VECM_bcols, rset_add_results, rset_get_R_matrix, rset_get_Ra_matrix, rset_get_q_matrix, rset_get_qa_matrix;
import gretl.gretl_restricth : GretlRestriction;
import gretl.johansenh : JohansenCode, JohansenInfo, auto_restr, jcode, jrank;
import gretl.libgretlh;
import gretl.libset : libset_get_int;
import gretl.matrix_extra : gretl_matrix_print_to_prn;
import gretl.libseth : SetKey, VECMnorm;
import gretl.missingh : NADBL, na;
import gretl.plugin.jalpha : vecm_alpha_test;
import gretl.plugin.jrestrict : general_vecm_analysis;
import gretl.printout : max_namelen_in_list, maybe_trim_varname;
import gretl.printouth : SHORTLEN;
import gretl.pvalues : chisq_cdf_comp, gamma_cdf_comp;
import gretl.var : gretl_VAR_do_error_decomp, gretl_VECM_rank, n_restricted_terms, transcribe_VAR_models;
import gretl.vartest : VAR_portmanteau_test;
import gretl.varh;
import gretl.varprint : print_Johansen_test_case, vecm_beta_varname;

/* coefficient matrices for the trace test */

immutable double[6][5] trace_m_coef = [
    /* N^2  N       1   N==1   N==2  N^1/2 */
    [ 2, -1.00,  0.07,  0.07,     0,     0 ],
    [ 2,  2.01,  0.00,  0.06,  0.05,     0 ],
    [ 2,  1.05, -1.55, -0.50, -0.23,     0 ],
    [ 2,  4.05,  0.50, -0.23, -0.07,     0 ],
    [ 2,  2.85, -5.10, -0.10, -0.06,  1.35 ]
];


immutable double[6][5] trace_v_coef = [
    [ 3, -0.33, -0.55,  0.0,  0.00,  0 ],
    [ 3,  3.60,  0.75, -0.4, -0.30,  0 ],
    [ 3,  1.80,  0.00, -2.8, -1.10,  0 ],
    [ 3,  5.70,  3.20, -1.3, -0.50,  0 ],
    [ 3,  4.00,  0.80, -5.8, -2.66,  0 ]
];

/* trace-test sample size correction: log(m)-log(m-hat) */

immutable double[7][5] trace_m_corr = [
    /* sqrt(n)/T n/T  n^2/T^2  n==1/T     n==1     n==2      n==3 */
    [ -0.101,  0.499,  0.896, -0.562,  0.00229, 0.00662,        0 ],
    [      0,  0.465,  0.984, -0.273, -0.00244,       0,        0 ],
    [  0.134,  0.422,   1.02,   2.17, -0.00182,       0, -0.00321 ],
    [ 0.0252,  0.448,   1.09, -0.353,        0,       0,        0 ],
    [ -0.819,  0.615,  0.896,   2.43,  0.00149,       0,        0 ]
];

/* trace-test sample size correction: log(v)-log(v-hat) */

immutable double[7][5] trace_v_corr = [
    [ -0.204,  0.980,  3.11,  -2.14,   0.0499,  -0.0103, -0.00902 ],
    [  0.224,  0.863,  3.38, -0.807,        0,        0,  -0.0091 ],
    [  0.422,  0.734,  3.76,   4.32, -0.00606,        0, -0.00718 ],
    [      0,  0.836,  3.99,  -1.33, -0.00298, -0.00139, -0.00268 ],
    [  -1.29,   1.01,  3.92,   4.67,  0.00484, -0.00127,  -0.0199 ]
];

/* coefficient matrices for the lambdamax test */

immutable double[5][5] maxev_m_coef = [
    /*   N        1        N==1       N==2     N^1/2 */
    [ 6.0019, -2.75580,  0.67185,  0.114900, -2.77640 ],
    [ 5.9498,  0.43402,  0.04836,  0.018198, -2.36690 ],
    [ 5.8271, -1.64870, -1.61180, -0.259490, -1.56660 ],
    [ 5.8658,  2.55950, -0.34443, -0.077991, -1.75520 ],
    [ 5.6364, -0.90531, -3.51660, -0.479660, -0.21447 ]
];

immutable double[5][5] maxev_v_coef = [
    [ 1.8806, -15.499,  1.11360,  0.070508, 14.714 ],
    [ 2.2231, -7.9064,  0.58592, -0.034324, 12.058 ],
    [ 2.0785, -9.7846, -3.36800, -0.245280, 13.074 ],
    [ 1.9955, -5.5428,  1.24250,  0.419490, 12.841 ],
    [ 2.0899, -5.3303, -7.15230, -0.252600, 12.393 ]
];

private void fill_x_asy_array(double * x, int n)
{
    x[0] = n * n;
    x[1] = n;
    x[2] = 1.0;
    x[3] = (n == 1)? 1.0 : 0.0;
    x[4] = (n == 2)? 1.0 : 0.0;
    x[5] = sqrt(cast(double) n);
}

private void fill_x_corr_array(double * x, int n, int T)
{
    x[0] = sqrt(cast(double) n) / T;
    x[1] = n / cast(double) T;
    x[2] = (n * n) / cast(double) (T * T);
    x[3] = (n == 1)? 1.0 / T : 0.0;
    x[4] = (n == 1)? 1.0 : 0.0;
    x[5] = (n == 2)? 1.0 : 0.0;
    x[6] = (n == 3)? 1.0 : 0.0;;
}

/* public function accessible via gretl's plugin apparatus

   In a script, do: pv = pvalue(J, n, det, T, trace)

   Use T = 0 for the plain asymptotic result.
*/

extern (C) double trace_pvalue(double trace, int n, int det, int T)
{
    if (det < JohansenCode.no_const || det > JohansenCode.unrest_trend || n < 1) {
        return NADBL;
    } else {
        const(double) * tracem = trace_m_coef[det].ptr;
        const(double) * tracev = trace_v_coef[det].ptr;
        double mt = 0, vt = 0;
        double[7] x;
        int i;

        fill_x_asy_array(x.ptr, n);

        for (i=0; i<6; i++) {
            mt += x[i] * tracem[i];
            vt += x[i] * tracev[i];
        }

        if (T > 0 && T < 10000) {
            const(double) * mcorr = trace_m_corr[det].ptr;
            const(double) * vcorr = trace_v_corr[det].ptr;
            double mtcorr = 0, vtcorr = 0;

            fill_x_corr_array(x.ptr, n, T);

            for (i=0; i<7; i++) {
                mtcorr += x[i] * mcorr[i];
                vtcorr += x[i] * vcorr[i];
            }

            mt *= exp(mtcorr);
            vt *= exp(vtcorr);
        }

        return gamma_cdf_comp(mt, vt, trace, 2);
    }
}

/*
   Asymptotic p-values for Johansen's likelihood ratio tests
   computed using J. Doornik's gamma approximation --
   see "Approximations to the Asymptotic Distributions of
   Cointegration Tests", Journal of Economic Surveys,
   12/5, 1998, pp. 573-593.

   trace, lmax: trace and lambdamax statistics

   det: index of setup of deterministic regressors
        J_NO_CONST     = no constant
        J_REST_CONST   = restricted constant
        J_UNREST_CONST = unrestricted constant
        J_REST_TREND   = restricted trend
        J_UNREST_TREND = unrestricted trend

   n: (>= 1), the number of potentially cointegrated
   variables minus the cointegration rank under the null

   pval: on output, array of p-values for the two tests
*/

private int gamma_LR_asy_pvals(double trace, double lmax, JohansenCode det,
                               int n, double * pval)
{
    const(double) * tracem = trace_m_coef[det].ptr;
    const(double) * tracev = trace_v_coef[det].ptr;
    const(double) * maxevm = maxev_m_coef[det].ptr;
    const(double) * maxevv = maxev_v_coef[det].ptr;
    double mt = 0, vt = 0;
    double ml = 0, vl = 0;
    double[6] x;
    int i;

    fill_x_asy_array(x.ptr, n);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
        if (i > 0) {
            ml += x[i] * maxevm[i-1];
            vl += x[i] * maxevv[i-1];
        }
    }

    pval[0] = gamma_cdf_comp(mt, vt, trace, 2);
    pval[1] = gamma_cdf_comp(ml, vl, lmax, 2);

    return 0;
}

/* Variant of gamma_LR_asy_pvals() with sample-size correction.
   Only the trace test is handled.
*/

private double gamma_LR_T_pval(double trace, JohansenCode det, int n, int T)
{
    const(double) * tracem = trace_m_coef[det].ptr;
    const(double) * tracev = trace_v_coef[det].ptr;
    const(double) * mcorr = trace_m_corr[det].ptr;
    const(double) * vcorr = trace_v_corr[det].ptr;
    double mt = 0, vt = 0;
    double mtcorr = 0, vtcorr = 0;
    double[7] x;
    int i;

    fill_x_asy_array(x.ptr, n);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
    }

    fill_x_corr_array(x.ptr, n, T);

    for (i=0; i<7; i++) {
        mtcorr += x[i] * mcorr[i];
        vtcorr += x[i] * vcorr[i];
    }

    mt = exp(log(mt) + mtcorr);
    vt = exp(log(vt) + vtcorr);

    return gamma_cdf_comp(mt, vt, trace, 2);
}

/*
  @n1: number of potentially cointegrated variables
  @n2: number of variables that are conditioned upon
  @r:  cointegrating rank under H0.
  @T:  sample size minus number of parameters, or
       0 for the asymptotic result

  The Doornik gamma-approximation approach with correction for the
  case of a "partial system" as discussed in Harbo, Johansen, Nielsen
  and Rahbek, "Asymptotic Inference on Cointegrating Rank in Partial
  Systems" (Journal of Business and Economic Statistics 16/4, October
  1998). The critical values were re-simulated in gretl with 50,000
  replications.
*/

private double gamma_harbo_trace_pval(double trace, JohansenCode det,
                                      int n1, int n2, int r, int T)
{
    const(double) * tracem = trace_m_coef[det].ptr;
    const(double) * tracev = trace_v_coef[det].ptr;
    double mt = 0, vt = 0;
    double cov;
    double[7] x;
    int n = n1 + n2;
    int i, m = n - r;

    fill_x_asy_array(x.ptr, m);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
    }

    if (T > 0) {
        const(double) * mcorr = trace_m_corr[det].ptr;
        const(double) * vcorr = trace_v_corr[det].ptr;
        double mtcorr = 0, vtcorr = 0;

        fill_x_corr_array(x.ptr, m, T);

        for (i=0; i<7; i++) {
            mtcorr += x[i] * mcorr[i];
            vtcorr += x[i] * vcorr[i];
        }

        mt = exp(log(mt) + mtcorr);
        vt = exp(log(vt) + vtcorr);
    }

    if (det == JohansenCode.rest_trend) {
        cov = -1.35;
    } else if (det == JohansenCode.rest_const) {
        cov = -1.066;
    } else {
        cov = -1.270;
    }

    mt *= (n1 - r) / cast(double) (n - r);
    vt *= (n1 - r) / cast(double) (n - r);
    vt -= (n - n1) * (n1 - r) * cov;

    return gamma_cdf_comp(mt, vt, trace, 2);
}

/* Remove a possible excess zero from the end of a floating point
   number printed to the given precision p (working around a bug
   in the C library).
*/

private void fix_xstr(char * s, int p)
{
    int n = cast(int) strlen(s);

    if (n > p && strspn(s + n - p, "0") == p) {
        s[n-1] = 0;
    }
}

enum ABMIN = 1.0e-15;

enum VBeta {
    alpha,
    beta
}

/* for cointegration test: print cointegrating vectors or adjustments,
   either "raw" or re-scaled */

private void print_beta_or_alpha(const(GretlVar) * jvar, int k,
                                 const(DatasetRaw) * dset, Prn * prn,
                                 int job, int rescale)
{
    JohansenInfo * jv = cast(JohansenInfo *) jvar.jinfo;
    GretlMatrixRaw * c = (job == VBeta.beta)? jv.Beta : jv.Alpha;
    int rows = gretl_matrix_rows(c);
    int vnorm = libset_get_int(SetKey.vecm_norm);
    char[32] xstr;
    char[SHORTLEN] tmp;
    int n, namelen = 8;
    int i, j, row;
    double x, y;

    if (vnorm == VECMnorm.none && rescale) {
        return;
    }

    if (rescale) {
        pprintf(prn, "\n%s\n", (job == VBeta.beta)?
                _("renormalized beta").ptr :
                _("renormalized alpha").ptr);
    } else {
        pprintf(prn, "\n%s\n", (job == VBeta.beta)?
                _("beta (cointegrating vectors)").ptr :
                _("alpha (adjustment vectors)").ptr);
    }

    for (i=0; i<rows; i++) {
        vecm_beta_varname(tmp.ptr, cast(GretlVar *) jvar, dset, i);
        n = cast(int) strlen(tmp.ptr);
        if (n > namelen) {
            namelen = n;
        }
    }

    for (i=0; i<rows; i++) {
        vecm_beta_varname(tmp.ptr, cast(GretlVar *) jvar, dset, i);
        pprintf(prn, "%-*s", namelen + 2, tmp.ptr);
        for (j=0; j<k; j++) {
            x = gretl_matrix_get(c, i, j);
            if (rescale) {
                row = (vnorm == VECMnorm.first)? 0 : j;
                y = gretl_matrix_get(jv.Beta, row, j);
                if (job == VBeta.beta) {
                    x /= y;
                } else {
                    x *= y;
                }
            }
            if (x == -0.0 || fabs(x) < ABMIN) {
                x = 0.0;
            }
            sprintf(xstr.ptr, "%#.5g", x);
            fix_xstr(xstr.ptr, 5);
            pprintf(prn, "%12s ", xstr.ptr);
        }
        pputc(prn, '\n');
    }
}

/* Calculate \alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554.  Not needed when doing a VECM (without restrictions
   on alpha), in which case we get \alpha via VECM_estimate_full()
   below.  We use this to support verbose output when testing a
   restriction on \beta.
*/

private int compute_alpha(JohansenInfo * jv)
{
    const(GretlMatrixRaw) * B = jv.Beta;
    GretlMatrixRaw * alpha = null;
    GretlMatrixRaw * BSB = null;
    GretlMatrixRaw * Tmp = null;
    int p = jv.S01.rows;
    int r = B.cols;
    int err = 0;

    BSB = gretl_matrix_alloc(r, r);
    Tmp = gretl_matrix_alloc(B.rows, r);
    alpha = gretl_matrix_alloc(p, r);

    if (BSB == null || Tmp == null || alpha == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = gretl_matrix_qform(B, GretlMatrixMod.transpose, jv.S11,
                                 BSB, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(BSB);
    }

    if (!err) {
        gretl_matrix_multiply(B, BSB, Tmp);
        gretl_matrix_multiply(jv.S01, Tmp, alpha);
    }

    gretl_matrix_free(BSB);
    gretl_matrix_free(Tmp);

    if (!err) {
        gretl_matrix_replace(&jv.Alpha, alpha);
    } else {
        gretl_matrix_free(alpha);
    }

    return err;
}

/* print the long-run matrix, \alpha \beta' */

private int print_long_run_matrix(const(GretlVar) * jvar,
                                  const(DatasetRaw) * dset,
                                  Prn * prn)
{
    JohansenInfo * jv = cast(JohansenInfo *) jvar.jinfo;
    GretlMatrixRaw * Pi;
    char[SHORTLEN] tmp;
    const(char) * vname;
    int firstlen = 10;
    int namelen;
    double x;
    int i, j;

    Pi = gretl_matrix_alloc(jv.Alpha.rows, jv.Beta.rows);
    if (Pi == null) {
        return GretlError.alloc;
    }

    namelen = max_namelen_in_list(jvar.ylist, dset);
    if (firstlen <= namelen) {
        firstlen = namelen + 1;
    }
    if (namelen < 12) {
        namelen = 12;
    }

    gretl_matrix_multiply_mod(jv.Alpha, GretlMatrixMod.none,
                              jv.Beta, GretlMatrixMod.transpose,
                              Pi, GretlMatrixMod.none);

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')").ptr);

    vname = dset.varname[jvar.ylist[1]];
    maybe_trim_varname(tmp.ptr, vname);

    pprintf(prn, "%*s", firstlen + namelen, tmp.ptr);
    for (j=1; j<Pi.cols; j++) {
        vecm_beta_varname(tmp.ptr, cast(GretlVar *) jvar, dset, j);
        pprintf(prn, "%*s", namelen + 1, tmp.ptr);
    }

    pputc(prn, '\n');

    for (i=0; i<Pi.rows; i++) {
        vname = dset.varname[jvar.ylist[i+1]];
        maybe_trim_varname(tmp.ptr, vname);
        pprintf(prn, "%-*s", firstlen, tmp.ptr);
        for (j=0; j<Pi.cols; j++) {
            x = gretl_matrix_get(Pi, i, j);
            if (fabs(x) < 0.5e-14) {
                x = 0.0;
            }
            pprintf(prn, "%#*.5g ", namelen, x);
        }
        pputc(prn, '\n');
    }

    pputc(prn, '\n');

    gretl_matrix_free(Pi);

    return 0;
}

/* Compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix.
*/

private int compute_omega(GretlVar * vecm)
{
    if (vecm.S == null) {
        vecm.S = gretl_matrix_alloc(vecm.neqns, vecm.neqns);
        if (vecm.S == null) {
            return GretlError.alloc;
        }
    }

    gretl_matrix_multiply_mod(vecm.E, GretlMatrixMod.transpose,
                              vecm.E, GretlMatrixMod.none,
                              vecm.S, GretlMatrixMod.none);

    gretl_matrix_divide_by_scalar(vecm.S, vecm.T);

    return 0;
}

private void gretl_matrix_I(GretlMatrixRaw * A, int n)
{
    int i;

    gretl_matrix_zero(A);
    for (i=0; i<n; i++) {
        gretl_matrix_set(A, i, i, 1.0);
    }
}

private bool lag_wanted(const(GretlVar) * v, int i)
{
    return v.lags == null || in_gretl_list(v.lags, i);
}

/* After doing OLS estimation of the VECM conditional on \beta: copy
   the coefficients on the lagged differences (i.e. form the \Gamma
   matrices) so we can compute the VAR representation */

private void copy_coeffs_to_Gamma(GretlVar * vecm, GretlArray * G)
{
    GretlMatrixRaw * Gk;
    int nl = var_n_lags(vecm);
    int i, j, k, h;
    double x;

    for (i=0; i<vecm.neqns; i++) {
        for (k=0; k<vecm.order; k++) {
            Gk = cast(GretlMatrixRaw *) gretl_array_get_data(G, k);
            if (!lag_wanted(vecm, k+1)) {
                gretl_matrix_zero(Gk);
                continue;
            }
            h = k + vecm.ifc;
            /* successive lags (distinct \Gamma_i matrices) */
            for (j=0; j<vecm.neqns; j++) {
                /* successive \Delta x_j */
                x = gretl_matrix_get(vecm.B, h, i);
                gretl_matrix_set(Gk, i, j, x);
                h += nl;
            }
        }
    }
}

/* \Pi, as will be used in forming the VAR representation */

private void form_Pi(GretlVar * v, GretlMatrixRaw * Pi)
{
    gretl_matrix_multiply_mod(v.jinfo.Alpha, GretlMatrixMod.none,
                              v.jinfo.Beta, GretlMatrixMod.transpose,
                              Pi, GretlMatrixMod.none);
}

/* After doing OLS estimation of the VECM conditional on \beta:
   copy the coefficients on the EC terms (\beta' X) into the \alpha
   matrix.
*/

private int OLS_to_alpha(GretlVar * v)
{
    int rank = v.jinfo.rank;
    int pos = v.ncoeff - rank;
    double x;
    int i, j;

    for (i=0; i<v.neqns; i++) {
        for (j=0; j<rank; j++) {
            x = gretl_matrix_get(v.B, pos + j, i);
            gretl_matrix_set(v.jinfo.Alpha, i, j, x);
        }
    }

    return 0;
}

/* flags for controlling "full" estimation of VECM */

enum VecmFlags {
    alpha_restricted = 1 << 0,
    bootstrapping    = 1 << 1
}

private bool bootstrap(int f) { return (f & VecmFlags.bootstrapping) != 0; }
private bool alpha_restricted(int f) { return (f & VecmFlags.alpha_restricted) != 0; }
private bool estimate_alpha(int f) { return !(f & VecmFlags.alpha_restricted); }

/* VAR representation: transcribe the coefficient matrix A_i (for lag
   i) into its place in the full VAR coefficient matrix, A
*/

private void add_Ai_to_VAR_A(GretlMatrixRaw * Ai, GretlVar * vecm,
                             int k, int flags)
{
    int i, j, offset = k * vecm.neqns;
    int tr = bootstrap(flags);
    double x;

    for (i=0; i<vecm.neqns; i++) {
        for (j=0; j<vecm.neqns; j++) {
            x = gretl_matrix_get(Ai, i, j);
            if (tr) {
                gretl_matrix_set(vecm.A, j + offset, i, x);
            } else {
                gretl_matrix_set(vecm.A, i, j + offset, x);
            }
        }
    }
}

/* write pre-computed ML alpha into model structs */

private void transcribe_alpha(GretlVar * v)
{
    Model * pmod;
    double aij, sij = NADBL;
    int r = jrank(v);
    int k = v.ncoeff - r;
    int i, j;

    for (i=0; i<v.neqns; i++) {
        pmod = v.models[i];
        for (j=0; j<r; j++) {
            aij = gretl_matrix_get(v.jinfo.Alpha, i, j);
            if (v.jinfo.Ase != null) {
                sij = gretl_matrix_get(v.jinfo.Ase, i, j);
            }
            pmod.coeff[k+j] = aij;
            pmod.sderr[k+j] = sij;
        }
    }
}

/* The X (data) and B (coefficient) matrices may need expanding
   to take account of the EC terms */

private int vecm_check_size(GretlVar * v, int flags)
{
    int xc = (v.X != null)? v.X.cols : 0;
    int err = 0;

    if (bootstrap(flags)) {
        /* We're re-visiting this function: v->ncoeff will
           already be full-size, as will the v->X matrix.
        */
        int dim = v.ncoeff;

        if (alpha_restricted(flags)) {
            dim -= gretl_VECM_rank(v);
        }

        if (dim > 0 && (v.X == null || v.B == null)) {
            gretl_errmsg_set("vecm_check_size: X and/or B wrong");
            return GretlError.data;
        } else if (dim == 0) {
            /* FIXME! */
            if (v.X != null) v.X.cols = 0;
            if (v.B != null) v.B.rows = 0;
        } else {
            v.X.cols = dim;
            v.B.rows = dim;
        }
        return 0;
    }

    v.ncoeff += jrank(v);

    if (estimate_alpha(flags)) {
        xc += jrank(v);
    }

    if (xc == 0) {
        /* nothing to be done */
        return 0;
    }

    if (v.X == null) {
        v.X = gretl_matrix_alloc(v.T, xc);
        if (v.X == null) {
            err = GretlError.alloc;
        } else {
            /* record allocated size */
            v.xcols = xc;
        }
    } else if (v.X.cols < xc) {
        if (v.xcols == xc) {
            gretl_matrix_reuse(v.X, -1, xc);
        } else {
            err = gretl_matrix_realloc(v.X, v.T, xc);
            if (!err) {
                /* record revised full size */
                v.xcols = xc;
            }
        }
    }

    if (!err) {
        int rows = v.ncoeff; /* 2021-01-06: was xc */

        if (v.B == null) {
            v.B = gretl_matrix_alloc(rows, v.neqns);
            if (v.B == null) {
                err = GretlError.alloc;
            }
        } else if (v.B.rows < rows) {
            err = gretl_matrix_realloc(v.B, rows, v.neqns);
        }
    }

    return err;
}

/* For estimating both alpha and Gamma: add the EC terms into
   the X data matrix */

private int add_EC_terms_to_X(GretlVar * v, GretlMatrixRaw * X,
                              const(DatasetRaw) * dset)
{
    const(GretlMatrixRaw) * B = v.jinfo.Beta;
    int rank = jrank(v);
    int k, k0 = v.ncoeff - rank;
    double xt, bxt, bij;
    int i, ii, j, s, t;
    int err = 0;

    for (j=0, k=k0; j<rank; j++, k++) {
        for (t=v.t1, s=0; t<=v.t2; t++, s++) {
            bxt = 0.0;
            ii = 0;

            /* beta * X(t-1) */
            for (i=0; i<v.neqns; i++) {
                xt = dset.Z[v.ylist[i+1]][t-1];
                bij = gretl_matrix_get(B, ii++, j);
                bxt += bij * xt;
            }

            /* restricted const or trend */
            if (auto_restr(v)) {
                bij = gretl_matrix_get(B, ii++, j);
                if (jcode(v) == JohansenCode.rest_trend) {
                    bij *= t;
                }
                bxt += bij;
            }

            /* restricted exog vars */
            if (v.rlist != null) {
                for (i=0; i<v.rlist[0]; i++) {
                    xt = dset.Z[v.rlist[i+1]][t]; /* was t-1 */
                    bij = gretl_matrix_get(B, ii++, j);
                    bxt += bij * xt;
                }
            }

            gretl_matrix_set(X, s, k, bxt);
        }
    }

    return err;
}

/* preparing for OLS conditional on beta: construct the
   appropriate dependent variable matrix, Y
*/

private int make_vecm_Y(GretlVar * v, const(DatasetRaw) * dset,
                        const(GretlMatrixRaw) * Pi)
{
    int i, s, t, vi, vj;
    double pij, yti, xti;
    int err = 0;

    if (v.Y == null) {
        fprintf(stderr, "make_vecm_Y: v->Y is NULL\n");
        return GretlError.data;
    }

    if (Pi == null) {
        /* "Y" is composed of plain DYt */
        for (i=0; i<v.neqns; i++) {
            vi = v.ylist[i+1];
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                yti = dset.Z[vi][t] - dset.Z[vi][t-1];
                gretl_matrix_set(v.Y, s++, i, yti);
            }
        }
    } else {
        /* net out \alpha: "Y" = DY_t - \Pi Y*_t */
        int j, k, wexo, p1 = v.jinfo.Beta.rows;

        for (i=0; i<v.neqns; i++) {
            wexo = 1;
            vi = v.ylist[i+1];
            s = 0;
            for (t=v.t1; t<=v.t2; t++) {
                /* first difference */
                yti = dset.Z[vi][t] - dset.Z[vi][t-1];
                for (j=0; j<p1; j++) {
                    pij = gretl_matrix_get(Pi, i, j);
                    if (pij != 0.0) {
                        if (j < v.neqns) {
                            /* lagged Y level */
                            wexo = 0;
                            vj = v.ylist[j+1];
                            xti = dset.Z[vj][t-1];
                        } else if (j == v.neqns && auto_restr(v)) {
                            xti = (jcode(v) == JohansenCode.rest_trend)? t : 1;
                        } else {
                            k = j - v.ylist[0] - auto_restr(v) + 1;
                            vj = v.rlist[k];
                            xti = dset.Z[vj][t]; /* was t-1 */
                        }
                        yti -= pij * xti;
                    }
                }
                gretl_matrix_set(v.Y, s++, i, yti);
            }
            if (wexo) {
                fprintf(stderr, "make_vecm_Y: var %d is weakly exogenous\n", i);
            }
        }
    }

    return err;
}

/* As in PcGive: df = T - c, where c is the "average number of
   estimated parameters per equation, rounded towards zero".
   Note that this function is _not_ invoked in the case of
   "general" restriction on alpha/beta; in that case the df
   calculation is handled specially.
*/

private void vecm_set_df(GretlVar * v, const(GretlMatrixRaw) * H,
                         const(GretlMatrixRaw) * R)
{
    int p = v.neqns;
    int r = v.jinfo.rank;
    int p1 = v.jinfo.Beta.rows;
    int Kpi, K = 0;
    double c;

    if (r == 0) {
        /* full-rank VAR */
        r = p;
    }

    /* lagged differences */
    K += var_n_lags(v) * p;

    /* deterministic stuff */
    K += v.jinfo.seasonals;
    if (jcode(v) > JohansenCode.rest_const) {
        K++;
    }
    if (jcode(v) == JohansenCode.unrest_trend) {
        K++;
    }

    /* unrestricted exogenous vars? */
    if (v.xlist != null) {
        K += v.xlist[0];
    }

    K *= p;

    /* df loss in \Pi = \alpha \beta': see Johansen (1995) chapter 7 */
    Kpi = (p + p1 - r) * r;
    if (H != null || R != null) {
        Kpi -= v.jinfo.lrdf;
    }
    K += Kpi;

    c = K / cast(double) p;
    v.df = v.T - cast(int) floor(c);
}

/* In the case where the EC coefficients were not estimated via
   OLS but precomputed in v->jinfo->Alpha, augment the B matrix
   with the alpha values. We need to add a number of rows equal
   to v->ncoeff minus the current B->rows, then transcribe the
   transpose of v->jinfo->Alpha into the new rows.
*/

private int var_B_insert_alpha(GretlVar * v)
{
    GretlMatrixRaw * A = v.jinfo.Alpha;
    GretlMatrixRaw * B = v.B;
    double aji;
    double * dest, src;
    int r_orig = B.rows;
    int r_full = v.ncoeff;
    size_t sz;
    int i, j;

    /* sanity check */
    if (r_full - r_orig <= 0 || A.rows != B.cols) {
        fprintf(stderr, "var_B_insert_alpha: matrix sizes wrong!\n");
        return GretlError.data;
    }

    sz = r_orig * (*B.val).sizeof;

    for (j=B.cols-1; j>0; j--) {
        dest = B.val + r_full * j;
        src = B.val + r_orig * j;
        memmove(dest, src, sz);
    }

    gretl_matrix_reuse(B, r_full, B.cols);

    for (j=0; j<B.cols; j++) {
        for (i=0; i<A.cols; i++) {
            aji = gretl_matrix_get(A, j, i);
            gretl_matrix_set(B, r_orig + i, j, aji);
        }
    }

    return 0;
}

/* The following is designed to accommodate the case where alpha is
   restricted, in which case we can't just run OLS conditional on
   beta.

   DY_t = \Pi Y*_{t-1} + \sum_{i=1}^{k-1}\Gamma_i DY_{t-1} + ...

   Subtract \Pi Y*_t from both sides, call DY_t - \Pi Y*_t "Y",
   and call the lagged DYs "X".  Regress "Y" on "X" to find
   estimates of the \Gammas, etc.

   But the function also handles the case where \alpha will
   be estimated along with \Gamma.
*/

private int VECM_estimate_full(GretlVar * v, const(GretlRestriction) * rset,
                               const(DatasetRaw) * dset, int flags)
{
    GretlMatrixRaw * beta = v.jinfo.Beta;
    GretlMatrixRaw * Pi = null;
    GretlMatrixRaw * Ai = null;
    GretlArray * G = null;
    int order = v.order;
    int xc, n = v.neqns;
    int use_XTX = 0;
    int fix_Y = 0;
    int i, err;

    if (alpha_restricted(flags) && v.jinfo.Alpha == null) {
        /* error: alpha must be pre-computed */
        fprintf(stderr, " alpha must be pre-computed, but Alpha is NULL\n");
        return GretlError.data;
    }

    err = vecm_check_size(v, flags);
    if (err) {
        fprintf(stderr, " error %d from vecm_check_size\n", err);
        return err;
    }

    xc = (v.X != null)? v.X.cols : 0;

    Pi = gretl_matrix_alloc(n, beta.rows);
    Ai = gretl_matrix_alloc(n, n);
    if (Pi == null || Ai == null) {
        err = GretlError.alloc;
    }

    if (!err && order > 0) {
        G = gretl_array_new(GretlType.matrices, order, &err);
        if (G == null) {
            err = GretlError.alloc;
        } else {
            for (i=0; i<order; i++) {
                gretl_array_set_data(G, i, gretl_matrix_alloc(n, n));
            }
        }
    }

    if (!err && alpha_restricted(flags)) {
        form_Pi(v, Pi);
    }

    if (!err && v.Y != null) {
        if (alpha_restricted(flags)) {
            fix_Y = 1;
            err = make_vecm_Y(v, dset, Pi);
        } else {
            err = make_vecm_Y(v, dset, null);
        }
        if (err) {
            fprintf(stderr, " error %d from make_vecm_Y\n", err);
        }
    }

    if (!err && estimate_alpha(flags)) {
        err = add_EC_terms_to_X(v, v.X, dset);
        if (err) {
            fprintf(stderr, " error %d from add_EC_terms_to_X\n", err);
        }
    }

    if (!err) {
        if (xc > 0) {
            /* run the regressions */
            if (v.B.rows > xc) {
                gretl_matrix_reuse(v.B, xc, -1);
            }
            if (bootstrap(flags)) {
                err = gretl_matrix_multi_ols(v.Y, v.X, v.B, v.E, null);
            } else {
                if (v.XTX != null) {
                    gretl_matrix_replace(&v.XTX, null);
                }
                err = gretl_matrix_multi_SVD_ols(v.Y, v.X, v.B, v.E, &v.XTX);
                use_XTX = 1;
            }
            if (err) {
                fprintf(stderr, " error %d from matrix_ols (xc = %d)\n", err, xc);
            }
        } else if (v.Y != null) {
            /* nothing to estimate, with alpha already in hand */
            gretl_matrix_copy_values(v.E, v.Y);
        } else {
            gretl_matrix_zero(v.E);
        }
    }

    if (!err && order > 0) {
        copy_coeffs_to_Gamma(v, G);
    }

    if (!err) {
        if (estimate_alpha(flags)) {
            err = OLS_to_alpha(v);
            if (!err) {
                form_Pi(v, Pi);
            }
        } else if (xc < v.ncoeff && v.B != null) {
            /* transcribe EC terms to v->B */
            var_B_insert_alpha(v);
        }
    }

    if (err) {
        goto bailout;
    }

    if (Pi.cols > n) {
        gretl_matrix_reuse(Pi, -1, n);
    }

    if (order == 0) {
        gretl_matrix_I(Ai, n);
        gretl_matrix_add_to(Ai, Pi);
        add_Ai_to_VAR_A(Ai, v, 0, flags);
    } else {
        GretlMatrixRaw * Gm;

        for (i=0; i<=order; i++) {
            if (i == 0) {
                gretl_matrix_I(Ai, n);
                gretl_matrix_add_to(Ai, Pi);
                Gm = cast(GretlMatrixRaw *) gretl_array_get_data(G, i);
                gretl_matrix_add_to(Ai, Gm);
            } else if (i == order) {
                gretl_matrix_zero(Ai);
                Gm = cast(GretlMatrixRaw *) gretl_array_get_data(G, i-1);
                gretl_matrix_subtract_from(Ai, Gm);
            } else {
                Gm = cast(GretlMatrixRaw *) gretl_array_get_data(G, i);
                gretl_matrix_copy_values(Ai, Gm);
                Gm = cast(GretlMatrixRaw *) gretl_array_get_data(G, i-1);
                gretl_matrix_subtract_from(Ai, Gm);
            }
            add_Ai_to_VAR_A(Ai, v, i, flags);
        }
    }

    if (!err && !bootstrap(flags)) {
        const(GretlMatrixRaw) * XTX = null;

        if (use_XTX || estimate_alpha(flags)) {
            XTX = v.XTX;
        }
        if (fix_Y) {
            /* ensure "Y" holds plain differences */
            make_vecm_Y(v, dset, null);
        }
        transcribe_VAR_models(v, dset, XTX);
        if (alpha_restricted(flags)) {
            transcribe_alpha(v);
        }
    }

    if (!err) {
        /* see Johansen (1995, p. 45) */
        GretlMatrixRaw * Tmp = gretl_identity_matrix_new(n);
        GretlMatrixRaw * Gi;

        for (i=0; i<order; i++) {
            Gi = cast(GretlMatrixRaw *) gretl_array_get_data(G, i);
            gretl_matrix_subtract_from(Tmp, Gi);
        }
        v.jinfo.Gamma = Tmp;
    }

 bailout:

    gretl_matrix_free(Pi);
    gretl_matrix_free(Ai);
    if (G != null) {
        v.jinfo.G = G;
    }

    return err;
}

/* Print both "raw" and re-scaled versions of the beta and alpha
   matrices (cointegrating vectors and vectors of adjustments
   respectively).
*/

private int print_beta_and_alpha(const(GretlVar) * jvar, GretlMatrixRaw * evals, int h,
                                 const(DatasetRaw) * dset, Prn * prn)
{
    int i, err = 0;

    pputs(prn, _("eigenvalue").ptr);
    for (i=0; i<h; i++) {
        pprintf(prn, "%#12.5g ", evals.val[i]);
    }
    pputc(prn, '\n');

    /* "raw" vectors */
    print_beta_or_alpha(jvar, h, dset, prn, VBeta.beta, 0);
    print_beta_or_alpha(jvar, h, dset, prn, VBeta.alpha, 0);

    /* re-scaled versions */
    print_beta_or_alpha(jvar, h, dset, prn, VBeta.beta, 1);
    print_beta_or_alpha(jvar, h, dset, prn, VBeta.alpha, 1);

    pputc(prn, '\n');

    return err;
}

/*
   renormalize \beta such that its uppermost submatrix of
   size rank * rank is the identity matrix:

   \beta' = [ I | *free elements* ]
*/

private int phillips_normalize_beta(GretlVar * vecm)
{
    GretlMatrixRaw * c = null;
    GretlMatrixRaw * beta_c = null;

    int r = jrank(vecm);
    int n = gretl_matrix_rows(vecm.jinfo.Beta);
    int i, j, err = 0;

    double x;

    c = gretl_matrix_alloc(r, r);
    beta_c = gretl_matrix_alloc(n, r);
    if (c == null || beta_c == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    for (i=0; i<r; i++) {
        for (j=0; j<r; j++) {
            x = gretl_matrix_get(vecm.jinfo.Beta, i, j);
            gretl_matrix_set(c, i, j, x);
        }
    }

    /* form \beta_c = \beta c^{-1} */
    err = gretl_invert_general_matrix(c);
    if (err) {
        fprintf(stderr, "phillips_normalize_beta: c is singular\n");
        goto bailout;
    }

    gretl_matrix_multiply(vecm.jinfo.Beta, c, beta_c);

    /* correct rounding error: set true zeros in \beta_c */
    for (i=0; i<n; i++) {
        for (j=0; j<r; j++) {
            if (i >= r) {
                if (gretl_matrix_get(beta_c, i, j) == -0) {
                    gretl_matrix_set(beta_c, i, j, 0);
                }
            } else if (i == j) {
                gretl_matrix_set(beta_c, i, j, 1.0);
            } else {
                gretl_matrix_set(beta_c, i, j, 0.0);
            }
        }
    }

    gretl_matrix_copy_values(vecm.jinfo.Beta, beta_c);

 bailout:

    gretl_matrix_free(c);
    gretl_matrix_free(beta_c);

    return err;
}

private int col_normalize_beta(GretlVar * vecm, int vnorm)
{
    GretlMatrixRaw * B = vecm.jinfo.Beta;
    double x, den;
    int i, j, row;

    for (j=0; j<B.cols; j++) {
        row = (vnorm == VECMnorm.diag)? j : 0;
        den = gretl_matrix_get(B, row, j);
        if (den != 0.0) {
            for (i=0; i<B.rows; i++) {
                x = gretl_matrix_get(B, i, j);
                gretl_matrix_set(B, i, j, x / den);
            }
        }
    }

    return 0;
}

private int normalize_beta(GretlVar * vecm, const(GretlMatrixRaw) * H,
                           int * do_stderrs)
{
    int vnorm = libset_get_int(SetKey.vecm_norm);

    if (vnorm == VECMnorm.none) {
        if (do_stderrs != null) {
            *do_stderrs = 0;
        }
        return 0;
    }

    if (H == null) {
        if (vnorm == VECMnorm.phillips) {
            return phillips_normalize_beta(vecm);
        } else {
            if (do_stderrs != null) {
                *do_stderrs = 0;
            }
            return col_normalize_beta(vecm, vnorm);
        }
    } else {
        GretlMatrixRaw * B = vecm.jinfo.Beta;

        if (B.cols == 1) {
            double den = B.val[0];
            int i;

            if (den != 0.0) {
                for (i=0; i<B.rows; i++) {
                    if (B.val[i] != 0) {
                        B.val[i] /= den;
                    }
                }
            }
        }
    }

    return 0;
}

private int restricted_beta_se(GretlVar * v, int r, int p1)
{
    double x;
    int i;

    v.jinfo.Bse = gretl_matrix_alloc(p1, r);
    if (v.jinfo.Bse == null) {
        return GretlError.alloc;
    }

    for (i=0; i<v.jinfo.Bvar.rows; i++) {
        x = gretl_matrix_get(v.jinfo.Bvar, i, i);
        v.jinfo.Bse.val[i] = sqrt(x);
    }

    return 0;
}

private int restricted_beta_variance(GretlVar * vecm,
                                     GretlMatrixRaw * Hin)
{
    GretlMatrixRaw * H = null;
    GretlMatrixRaw * O = null;
    GretlMatrixRaw * aOa = null;
    GretlMatrixRaw * K = null;
    GretlMatrixRaw * Vphi = null;
    int r = jrank(vecm);
    int p = vecm.neqns;
    int p1 = p + n_restricted_terms(vecm);
    int nb = r * p1;
    int freeH = 0;
    int err = 0;

    if (r > 1) {
        H = gretl_matrix_I_kronecker_new(r, Hin, &err);
        if (err) {
            return err;
        }
        freeH = 1;
    } else {
        H = Hin;
    }

    clear_gretl_matrix_err();

    O = gretl_matrix_copy(vecm.S);
    aOa = gretl_matrix_alloc(r, r);
    K = gretl_matrix_alloc(nb, nb);
    Vphi = gretl_matrix_alloc(H.cols, H.cols);

    err = get_gretl_matrix_err();
    if (err) {
        goto bailout;
    }

    err = gretl_invert_symmetric_matrix(O);

    if (!err) {
        err = gretl_matrix_qform(vecm.jinfo.Alpha, GretlMatrixMod.transpose,
                                 O, aOa, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_matrix_kronecker_product(aOa, vecm.jinfo.S11, K);
    }

    if (!err) {
        gretl_matrix_qform(H, GretlMatrixMod.transpose, K,
                           Vphi, GretlMatrixMod.none);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(Vphi);
    }

    if (!err) {
        gretl_matrix_divide_by_scalar(Vphi, vecm.df);
    }

    if (!err) {
        vecm.jinfo.Bvar = gretl_matrix_alloc(nb, nb);
        if (vecm.jinfo.Bvar == null) {
            err = GretlError.alloc;
        } else {
            gretl_matrix_qform(H, GretlMatrixMod.none, Vphi,
                               vecm.jinfo.Bvar,
                               GretlMatrixMod.none);
        }
    }

    if (!err) {
        err = restricted_beta_se(vecm, r, p1);
    }

 bailout:

    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(K);
    gretl_matrix_free(Vphi);

    if (freeH) {
        gretl_matrix_free(H);
    }

    return err;
}

/* VECM: compute the variance of the estimator of \beta, after doing
   Phillips normalization */

private int beta_variance(GretlVar * vecm)
{
    GretlMatrixRaw * O = null;
    GretlMatrixRaw * aOa = null;
    GretlMatrixRaw * HSH = null;
    double x;
    int r = jrank(vecm);
    int p1 = gretl_matrix_rows(vecm.jinfo.Beta);
    int nh = p1 - r;
    int i, j, k, err = 0;

    O = gretl_matrix_copy(vecm.S);
    aOa = gretl_matrix_alloc(r, r);
    HSH = gretl_matrix_alloc(nh, nh);

    if (O == null || aOa == null || HSH == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* compute \alpha' \Omega^{-1} \alpha */

    err = gretl_invert_symmetric_matrix(O);
    if (err) {
        goto bailout;
    }

    gretl_matrix_qform(vecm.jinfo.Alpha, GretlMatrixMod.transpose, O,
                       aOa, GretlMatrixMod.none);

    /* form H'*S11*H: just keep the south-east corner */

    for (i=r; i<p1; i++) {
        for (j=r; j<p1; j++) {
            x = gretl_matrix_get(vecm.jinfo.S11, i, j);
            gretl_matrix_set(HSH, i - r, j - r, x);
        }
    }

    vecm.jinfo.Bvar = gretl_matrix_kronecker_product_new(aOa, HSH, &err);
    if (err) {
        goto bailout;
    }

    err = gretl_invert_symmetric_matrix(vecm.jinfo.Bvar);
    if (err) {
        goto bailout;
    }

    vecm.jinfo.Bse = gretl_zero_matrix_new(p1, r);
    if (vecm.jinfo.Bse == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* note use of df rather than T here */
    gretl_matrix_divide_by_scalar(vecm.jinfo.Bvar, vecm.df);

    k = 0;
    for (j=0; j<r; j++) {
        /* cointegrating vector j */
        for (i=r; i<p1; i++) {
            x = gretl_matrix_get(vecm.jinfo.Bvar, k, k);
            gretl_matrix_set(vecm.jinfo.Bse, i, j, sqrt(x));
            k++;
        }
    }

 bailout:

    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(HSH);

    return err;
}

int johansen_ll_calc(GretlVar * jvar, const(GretlMatrixRaw) * evals)
{
    GretlMatrixRaw * S00;
    int n = jvar.neqns;
    int r = jrank(jvar);
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar.jinfo.S00);

    if (S00 == null) {
        err = GretlError.alloc;
    } else {
        double ldet = gretl_matrix_log_determinant(S00, &err);

        jvar.ll = n * (1.0 + LN_2_PI) + ldet;
        for (i=0; i<h; i++) {
            jvar.ll += log(1.0 - evals.val[i]);
        }
        jvar.ll *= -(jvar.T / 2.0);
        gretl_matrix_free(S00);
    }

    return err;
}

private int vecm_ll_stats(GretlVar * vecm)
{
    int T = vecm.T;
    int g = vecm.neqns;
    int k = g * (vecm.order + 1); /* FIXME gappy */
    int err = 0;

    vecm.ldet = gretl_vcv_log_determinant(vecm.S, &err);
    if (err) {
        return err;
    }

    k += vecm.jinfo.seasonals;

    /* FIXME: is the following right for k? */

    if (jcode(vecm) >= JohansenCode.unrest_const) {
        k++;
    }
    if (jcode(vecm) == JohansenCode.unrest_trend) {
        k++;
    }
    if (vecm.xlist != null) {
        k += vecm.xlist[0];
    }

    k *= g;

    vecm.AIC = (-2.0 * vecm.ll + 2.0 * k) / T;
    vecm.BIC = (-2.0 * vecm.ll + log(cast(double) T) * k) / T;
    vecm.HQC = (-2.0 * vecm.ll + 2.0 * log(log(cast(double) T)) * k) / T;

    VAR_portmanteau_test(vecm);

    return 0;
}

private void coint_test_print_exog(GretlVar * jvar,
                                   const(DatasetRaw) * dset,
                                   Prn * prn)
{
    int i, vi;

    if (jvar.xlist != null && jvar.xlist[0] > 0) {
        pprintf(prn, "\n%s: ", _("Exogenous regressor(s)").ptr);
        for (i=1; i<=jvar.xlist[0]; i++) {
            vi = jvar.xlist[i];
            pprintf(prn, "%s ", dset.varname[vi]);
        }
    }
}

private int compute_coint_test(GretlVar * jvar, const(DatasetRaw) * dset,
                               Gretlopt opt, Prn * prn)
{
    GretlMatrixRaw * evals = jvar.jinfo.evals;
    GretlMatrixRaw * tests;
    GretlMatrixRaw * pvals;
    int T = jvar.T;
    int n = jvar.neqns;
    int nexo = 0, nrexo = 0;
    int partial = 0;
    double llc, trace, lmax;
    double cumeig = 0.0;
    JohansenCode jcase;
    int i;

    tests = gretl_matrix_alloc(n, 2);
    pvals = gretl_matrix_alloc(n, 2);

    if (tests == null || pvals == null) {
        gretl_matrix_free(tests);
        gretl_matrix_free(pvals);
        return GretlError.alloc;
    }

    for (i=n-1; i>=0; i--) {
        lmax = -T * log(1.0 - evals.val[i]);
        cumeig += lmax;
        trace = cumeig;
        gretl_matrix_set(tests, i, 0, trace);
        gretl_matrix_set(tests, i, 1, lmax);
    }

    if (jvar.xlist != null) {
        nexo = jvar.xlist[0];
    }

    if (jvar.rlist != null) {
        nrexo = jvar.rlist[0];
    }

    jcase = jcode(jvar);

    if (nrexo > 0 && jcase != JohansenCode.unrest_trend) {
        /* use gamma-approx for partial system */
        partial = 1;
    }

    print_Johansen_test_case(jcase, prn);
    if (nexo > 0) {
        coint_test_print_exog(jvar, dset, prn);
    }
    pputc(prn, '\n');

    llc = (1.0 + LN_2_PI) * jvar.T;
    pprintf(prn, "\n%s = %g (%s: %g)\n", _("Log-likelihood").ptr,
            jvar.ll + llc, _("including constant term").ptr, jvar.ll);

    /* make sure the df is calculated */
    vecm_set_df(jvar, null, null);

    if (partial || nexo > 0) {
        pputc(prn, '\n');
        pputs(prn, _("Cointegration tests, ignoring exogenous variables").ptr);
    }

    pprintf(prn, "\n%s %s %s  %s  %s  %s\n", _("Rank").ptr, _("Eigenvalue").ptr,
            _("Trace test").ptr, _("p-value").ptr,
            _("Lmax test").ptr, _("p-value").ptr);

    for (i=0; i<n; i++) {
        double[2] pv = 0;

        trace = gretl_matrix_get(tests, i, 0);
        lmax = gretl_matrix_get(tests, i, 1);
        gamma_LR_asy_pvals(trace, lmax, jcase, n - i, pv.ptr);
        pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]%#11.5g [%6.4f]\n",
                i, evals.val[i], trace, pv[0], lmax, pv[1]);
        gretl_matrix_set(pvals, i, 0, pv[0]);
        gretl_matrix_set(pvals, i, 1, pv[1]);
    }

    if (partial) {
        pputc(prn, '\n');
        pprintf(prn, _("Cointegration tests conditional on %d I(1) variable(s)").ptr,
                nrexo);
        pprintf(prn, "\n%s %s %s  %s  %s\n", _("Rank").ptr, _("Eigenvalue").ptr,
                _("Trace test").ptr, _("p-value").ptr, _("pval(T)").ptr);

        for (i=0; i<n; i++) {
            trace = gretl_matrix_get(tests, i, 0);
            double pv;
            JohansenCode det = jcase;

            /* bodge for comparability with CATS */
            if (jcase == JohansenCode.unrest_const) {
                det = JohansenCode.rest_trend;
            }
            pv = gamma_harbo_trace_pval(trace, det, n, nrexo, i, 0);
            pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]", i,
                    evals.val[i], trace, pv);
            gretl_matrix_set(pvals, i, 0, pv);
            gretl_matrix_set(pvals, i, 1, NADBL);
            pv = gamma_harbo_trace_pval(trace, det, n, nrexo, i, jvar.df);
            pprintf(prn, " [%6.4f]\n", pv);
        }
    } else {
        /* !partial */
        double pv;

        pputc(prn, '\n');
        pprintf(prn, _("Corrected for sample size (df = %d)").ptr, jvar.df);
        pprintf(prn, "\n%s %s %s\n", _("Rank").ptr, _("Trace test").ptr,
                _("p-value").ptr);
        for (i=0; i<n; i++) {
            trace = gretl_matrix_get(tests, i, 0);
            pv = gamma_LR_T_pval(trace, jcase, n - i, jvar.df);
            pprintf(prn, "%4d%#11.5g [%6.4f]\n", i, trace, pv);
            if (!(opt & Gretlopt.OPT_Y)) {
                /* not asymptotic: prefer sample-corrected p-value */
                gretl_matrix_set(pvals, i, 0, pv);
            }
        }
    }

    pputc(prn, '\n');

    if (nexo > 0 || nrexo > 0) {
        if (partial && jcase == JohansenCode.unrest_const) {
            pputs(prn, _("Warning: the p-values shown are for the "
                  ~ "case of\na restricted trend").ptr);
            pputs(prn, "\n\n");
        } else if (!partial) {
            pputs(prn, _("Note: in general, the test statistics above "
                         ~ "are valid only in the\nabsence of additional "
                         ~ "regressors.").ptr);
            pputs(prn, "\n\n");
        }
    }

    record_matrix_test_result(tests, pvals);

    return 0;
}

private int johansen_get_eigenvalues(GretlMatrixRaw * S00,
                                     const(GretlMatrixRaw) * S01,
                                     const(GretlMatrixRaw) * S11,
                                     GretlMatrixRaw ** M,
                                     GretlMatrixRaw ** evals,
                                     int rank)
{
    GretlMatrixRaw * Tmp = null;
    int n = S11.cols;
    int err;

    err = gretl_invert_symmetric_matrix(S00);
    if (err) {
        return err;
    }

    Tmp = gretl_matrix_alloc(n, n);
    if (Tmp == null) {
        return GretlError.alloc;
    }

    *M = gretl_matrix_alloc(n, n);
    if (*M == null) {
        gretl_matrix_free(Tmp);
        return GretlError.alloc;
    }

    gretl_matrix_qform(S01, GretlMatrixMod.transpose,
                       S00, Tmp, GretlMatrixMod.none);

    *evals = gretl_gensymm_eigenvals(Tmp, S11, *M, &err);

    if (!err) {
        err = gretl_symmetric_eigen_sort(*evals, *M, rank);
    }

    gretl_matrix_free(Tmp);

    return err;
}

/* Public entry point for cointegration test */

extern (C) int johansen_coint_test(GretlVar * jvar, const(DatasetRaw) * dset,
                        Gretlopt opt, Prn * prn)
{
    int p1 = jvar.jinfo.R1.cols;
    int p = jvar.neqns;
    int err = 0;

    jvar.jinfo.Beta = gretl_matrix_alloc(p1, p);
    jvar.jinfo.Alpha = gretl_matrix_alloc(p, p);
    jvar.jinfo.evals = gretl_column_vector_alloc(p);

    if (jvar.jinfo.Beta == null ||
        jvar.jinfo.Alpha == null ||
        jvar.jinfo.evals == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = gretl_matrix_SVD_johansen_solve(jvar.jinfo.R0,
                                              jvar.jinfo.R1,
                                              jvar.jinfo.evals,
                                              jvar.jinfo.Beta,
                                              jvar.jinfo.Alpha, 0);
    }

    if (err) {
        pputs(prn, _("Failed to find eigenvalues\n").ptr);
    } else {
        johansen_ll_calc(jvar, jvar.jinfo.evals);
        compute_coint_test(jvar, dset, opt, prn);

        if (!(opt & Gretlopt.OPT_Q)) {
            print_beta_and_alpha(jvar, jvar.jinfo.evals, p,
                                 dset, prn);
            print_long_run_matrix(jvar, dset, prn);
        }
    }

    return err;
}

private void set_beta_test_df(GretlVar * jvar, const(GretlMatrixRaw) * H)
{
    int r = jrank(jvar);
    int p1 = jvar.jinfo.Beta.rows;
    int s = H.cols;

    jvar.jinfo.lrdf = r * (p1 - s);
}

/* Likelihood ratio test calculation, for restriction on
   an existing VECM -- which is not modified */

int johansen_LR_calc(const(GretlVar) * jvar, const(GretlMatrixRaw) * evals,
                     const(GretlMatrixRaw) * H, GretlRestriction * rset,
                     int job, Prn * prn)
{
    GretlMatrixRaw * S00;
    double llr = 0.0;
    double ldet = 0.0;
    double T_2 = cast(double) jvar.T / 2.0;
    int n = jvar.neqns;
    int r = jrank(jvar);
    int s = H.cols;
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar.jinfo.S00);

    if (S00 == null) {
        err = GretlError.alloc;
    } else {
        ldet = gretl_matrix_log_determinant(S00, &err);
    }

    if (!err) {
        llr = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
        for (i=0; i<h; i++) {
            pprintf(prn, _("eigenvalue %d = %g\n").ptr, i+1, evals.val[i]);
            llr -= T_2 * log(1.0 - evals.val[i]);
        }
        pputc(prn, '\n');
    }

    if (S00 != null) {
        gretl_matrix_free(S00);
    }

    if (!err) {
        double x = 2.0 * (jvar.ll - llr);
        int nb = gretl_matrix_rows(jvar.jinfo.Beta);
        int df;

        if (job == VBeta.beta) {
            df = h * (nb - s);
        } else {
            df = h * (n - s);
        }

        /* allow for possible prior restriction */
        df -= jvar.jinfo.lrdf;

        pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g\n").ptr, jvar.ll);
        pprintf(prn, _("Restricted loglikelihood (lr) = %.8g\n").ptr, llr);
        pprintf(prn, "2 * (lu - lr) = %g\n", x);

        if (df > 0) {
            double pv = chisq_cdf_comp(df, x);

            if (jvar.jinfo.lrdf > 0) {
                pprintf(prn, _("Allowing for prior restriction, df = %d\n").ptr, df);
            }
            pprintf(prn, "P(%s(%d) > %g) = %g\n", _("Chi-square").ptr, df, x, pv);
            rset_add_results(rset, x, pv, llr);
        }
    }

    return err;
}

enum NSMIN = 1.0e-16;

private GretlMatrixRaw * johansen_nullspace(const(GretlMatrixRaw) * R,
                                            int * err)
{
    GretlMatrixRaw * H = gretl_matrix_right_nullspace(R, err);

    if (!*err && H.cols == 1) {
        /* this normalization used to be applied in
           gretl_matrix.c, prior to 2023-05-08
        */
        double hij, x = 0;
        int i;

        for (i=0; i<H.rows; i++) {
            if (fabs(H.val[i]) > x) {
                x = H.val[i];
            }
        }
        for (i=0; i<H.rows; i++) {
            hij = H.val[i] / x;
            H.val[i] = fabs(hij) < NSMIN ? 0 : hij;
        }
    }

    return H;
}

/* see Johansen (1995), section 7.2: compute the p x s
   direct restriction matrix H = R_\perp and form the
   following two intermediate matrices:

   H' * S_{11} * H  in "S11"
   S_{01} * H       in "S01"
*/

private int prep_beta_restriction(const(GretlVar) * jvar,
                                  const(GretlMatrixRaw) * R,
                                  GretlMatrixRaw ** S01,
                                  GretlMatrixRaw ** S11,
                                  GretlMatrixRaw ** pH)
{
    GretlMatrixRaw * H;
    int s, p = jvar.neqns;
    int err = 0;

    if (R == null) {
        return GretlError.data;
    }

    H = johansen_nullspace(R, &err);
    if (err) {
        return err;
    }

    *pH = H;
    s = gretl_matrix_cols(H);

    *S11 = gretl_matrix_alloc(s, s);
    *S01 = gretl_matrix_alloc(p, s);
    if (*S11 == null || *S01 == null) {
        return GretlError.alloc;
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GretlMatrixMod.transpose,
                             jvar.jinfo.S11,
                             *S11, GretlMatrixMod.none);

    if (!err) {
        /* S01 <- S01H */
        err = gretl_matrix_multiply(jvar.jinfo.S01, H, *S01);
    }

    return err;
}

/* check whether the restriction on @jvar in @rset is
   a simple restriction on \beta */

private int simple_beta_restriction(const(GretlVar) * jvar,
                                    const(GretlRestriction) * rset)
{
    int ret = 0;

    if (rset_VECM_acols(rset) == 0) {
        const(GretlMatrixRaw) * R = rset_get_R_matrix(rset);
        const(GretlMatrixRaw) * q = rset_get_q_matrix(rset);
        int rcols = jvar.neqns + n_restricted_terms(cast(GretlVar *) jvar);

        ret = 1;

        if (!gretl_is_zero_matrix(q)) {
            /* non-homogeneous */
            ret = 0;
        } else if (R.cols > rcols) {
            /* not common to all columns */
            ret = 0;
        }
    }

    return ret;
}

/* check whether the restriction on @jvar in @rset is
   a simple restriction on \alpha */

private int simple_alpha_restriction(const(GretlVar) * jvar,
                                     const(GretlRestriction) * rset)
{
    int ret = 0;

    if (rset_VECM_bcols(rset) == 0) {
        const(GretlMatrixRaw) * Ra = rset_get_Ra_matrix(rset);
        const(GretlMatrixRaw) * qa = rset_get_qa_matrix(rset);

        ret = 1;

        if (!gretl_is_zero_matrix(qa)) {
            /* non-homogeneous */
            ret = 0;
        } else if (Ra.cols > jvar.neqns) {
            /* not common to all columns */
            ret = 0;
        }
    }

    return ret;
}

private int transcribe_restriction_matrices(const(GretlVar) * jvar,
                                            const(GretlRestriction) * rset)
{
    JohansenInfo * jv = cast(JohansenInfo *) jvar.jinfo;
    int err = 0;

    if (rset_VECM_bcols(rset) > 0) {
        const(GretlMatrixRaw) * R = rset_get_R_matrix(rset);
        const(GretlMatrixRaw) * q = rset_get_q_matrix(rset);

        if (R != jv.R) {
            gretl_matrix_replace(&jv.R, gretl_matrix_copy(R));
        }

        if (q != jv.q) {
            gretl_matrix_replace(&jv.q, gretl_matrix_copy(q));
        }

        if (jv.R == null ||
            (q != null && jv.q == null)) {
            err = GretlError.alloc;
        }
    }

    if (!err && rset_VECM_acols(rset) > 0) {
        const(GretlMatrixRaw) * Ra = rset_get_Ra_matrix(rset);
        const(GretlMatrixRaw) * qa = rset_get_qa_matrix(rset);

        if (Ra != jv.Ra) {
            gretl_matrix_replace(&jv.Ra, gretl_matrix_copy(Ra));
        }

        if (qa != jv.qa) {
            gretl_matrix_replace(&jv.qa, gretl_matrix_copy(qa));
        }

        if (jv.Ra == null ||
            (qa != null && jv.qa == null)) {
            err = GretlError.alloc;
        }
    }

    return err;
}

/* driver for VECM estimation subject to "general" restrictions
   on beta and/or alpha */

private int j_general_restrict(GretlVar * jvar, GretlRestriction * rset,
                               const(DatasetRaw) * dset, int flags,
                               Prn * prn)
{
    int err = 0;

    err = general_vecm_analysis(jvar, rset, dset, prn);

    if (!err) {
        if (rset_VECM_acols(rset) > 0) {
            flags |= VecmFlags.alpha_restricted;
        }
        err = VECM_estimate_full(jvar, rset, dset, flags);
    }

    if (!err) {
        err = gretl_VAR_do_error_decomp(jvar.S, jvar.C, null);
    }

    if (!err && !bootstrap(flags)) {
        /* FIXME 'k' for AIC etc? */
        err = vecm_ll_stats(jvar);
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
    }

    return err;
}

/* Obtain the unrestricted log-likelihood for running the LR test.  We
   need do this (only) in the context where we're doing full
   estimation of a restricted system.  The prior system, relative to
   which the (new) restriction is defined, may have been restricted
   already, in which case the unrestricted ll is not available for
   comparison.

   This function is low-budget in that we don't bother with the
   eigenvectors, just the eigenvalues.
*/

private int get_unrestricted_ll(GretlVar * jvar)
{
    GretlMatrixRaw * S00 = null;
    GretlMatrixRaw * Tmp = null;
    GretlMatrixRaw * e = null;
    double ldet;
    int n1 = jvar.jinfo.S11.cols;
    int n = jvar.neqns;
    int r = jrank(jvar);
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar.jinfo.S00);
    if (S00 == null) {
        return GretlError.alloc;
    }

    Tmp = gretl_matrix_alloc(n1, n1);
    if (Tmp == null) {
        gretl_matrix_free(S00);
        return GretlError.alloc;
    }

    err = gretl_invert_symmetric_matrix(S00);

    if (!err) {
        gretl_matrix_qform(jvar.jinfo.S01, GretlMatrixMod.transpose,
                           S00, Tmp, GretlMatrixMod.none);
        e = gretl_gensymm_eigenvals(Tmp, jvar.jinfo.S11, null, &err);
    }

    if (!err) {
        gretl_matrix_copy_values(S00, jvar.jinfo.S00);
        ldet = gretl_matrix_log_determinant(S00, &err);
    }

    if (!err) {
        qsort(e.val, n1, (*e.val).sizeof, &gretl_inverse_compare_doubles);
        jvar.jinfo.ll0 = n * (1.0 + LN_2_PI) + ldet;
        for (i=0; i<r; i++) {
            jvar.jinfo.ll0 += log(1.0 - e.val[i]);
        }
        jvar.jinfo.ll0 *= -(jvar.T / 2.0);
    }

    gretl_matrix_free(S00);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(e);

    return err;
}

/* common finalization for estimation subject to simple beta
   restriction, simple alpha restriction, or no restriction.
*/

private int vecm_finalize(GretlVar * jvar, GretlMatrixRaw * H,
                          const(GretlMatrixRaw) * Ra,
                          const(DatasetRaw) * dset,
                          int flags)
{
    int do_stderrs = jrank(jvar) < jvar.neqns;
    int err = 0;

    if (Ra == null) {
        /* alpha unrestricted */
        err = normalize_beta(jvar, H, &do_stderrs);
    } else {
        do_stderrs = 0;
    }

    if (!err) {
        if (bootstrap(flags)) {
            do_stderrs = 0;
        } else {
            vecm_set_df(jvar, H, Ra);
        }
        err = VECM_estimate_full(jvar, null, dset, flags);
    }

    if (!err) {
        err = compute_omega(jvar);
    }

    if (!err && do_stderrs) {
        if (H != null) {
            err = restricted_beta_variance(jvar, H);
        } else {
            err = beta_variance(jvar);
        }
    }

    if (!err) {
        err = gretl_VAR_do_error_decomp(jvar.S, jvar.C, null);
    }

    if (!err && !bootstrap(flags)) {
        err = vecm_ll_stats(jvar);
    }

    return err;
}

/* estimation subject to "simple" restriction on alpha */

private int est_simple_alpha_restr(GretlVar * jvar, GretlRestriction * rset,
                                   const(DatasetRaw) * dset, int flags,
                                   Prn * prn)
{
    const(GretlMatrixRaw) * Ra = rset_get_Ra_matrix(rset);
    Gretlopt opt = Gretlopt.OPT_F;
    int err = 0;

    if (bootstrap(flags)) {
        opt |= Gretlopt.OPT_B;
    } else {
        err = get_unrestricted_ll(jvar);
    }

    if (!err) {
        err = vecm_alpha_test(jvar, rset, dset, opt, prn);
    }

    if (!err) {
        err = vecm_finalize(jvar, null, Ra, dset,
                            flags | VecmFlags.alpha_restricted);
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
    }

    return err;
}

/* full estimation subject to "simple" restriction on beta */

private int est_simple_beta_restr(GretlVar * jvar, const(GretlRestriction) * rset,
                                  const(DatasetRaw) * dset, int flags)
{
    const(GretlMatrixRaw) * R = null;
    GretlMatrixRaw * H = null;
    GretlMatrixRaw * M = null;
    GretlMatrixRaw * S00 = null;
    GretlMatrixRaw * S01 = null;
    GretlMatrixRaw * S11 = null;
    GretlMatrixRaw * evals = null;
    int r = jrank(jvar);
    int err = 0;

    if (!bootstrap(flags)) {
        err = get_unrestricted_ll(jvar);
    }

    if (!err) {
        R = rset_get_R_matrix(rset);
        err = prep_beta_restriction(jvar, R, &S01, &S11, &H);
    }

    if (!err) {
        S00 = gretl_matrix_copy(jvar.jinfo.S00);
        if (S00 == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00, S01, S11,
                                       &M, &evals, r);
    }

    if (!err) {
        if (jvar.jinfo.Beta == null) {
            jvar.jinfo.Beta = gretl_matrix_multiply_new(H, M, &err);
        } else {
            err = gretl_matrix_multiply(H, M, jvar.jinfo.Beta);
        }
    }

    if (!err && !bootstrap(flags)) {
        set_beta_test_df(jvar, H);
        err = johansen_ll_calc(jvar, evals);
    }

    if (!err) {
        err = vecm_finalize(jvar, H, null, dset, flags);
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
    }

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);
    gretl_matrix_free(S00);
    gretl_matrix_free(S01);
    gretl_matrix_free(S11);

    return err;
}

/* "unrestricted" VECM estimation */

private int j_estimate_unrestr(GretlVar * jvar,
                               const(DatasetRaw) * dset)
{
    GretlMatrixRaw * S00 = null;
    GretlMatrixRaw * evals = null;
    int r = jrank(jvar);
    int err = 0;

    S00 = gretl_matrix_copy(jvar.jinfo.S00);
    if (S00 == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00,
                                       jvar.jinfo.S01,
                                       jvar.jinfo.S11,
                                       &jvar.jinfo.Beta,
                                       &evals, r);
    }

    if (!err) {
        err = johansen_ll_calc(jvar, evals);
    }

    if (!err) {
        err = vecm_finalize(jvar, null, null, dset, 0);
    }

    if (!err) {
        if (evals.rows > r) {
            /* include only the non-zero values */
            evals.rows = r;
        }
        jvar.jinfo.evals = evals;
    } else {
        gretl_matrix_free(evals);
    }

    gretl_matrix_free(S00);

    return err;
}

/* Here we prep the system with the initial eigen-analysis, then
   basically hand over to jrestrict.c
*/

private int j_estimate_general(GretlVar * jvar, GretlRestriction * rset,
                               const(DatasetRaw) * dset, int flags,
                               Prn * prn)
{
    GretlMatrixRaw * S00 = null;
    GretlMatrixRaw * evals = null;
    GretlMatrixRaw * M = null;
    int r = jrank(jvar);
    int err = 0;

    S00 = gretl_matrix_copy(jvar.jinfo.S00);
    if (S00 == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00,
                                       jvar.jinfo.S01,
                                       jvar.jinfo.S11,
                                       &M, &evals, r);
    }

    if (bootstrap(flags)) {
        if (M != null) {
            gretl_matrix_copy_values(jvar.jinfo.Beta, M);
            gretl_matrix_free(M);
        }
    } else {
        jvar.jinfo.Beta = M;
    }

    if (!err && !bootstrap(flags)) {
        err = johansen_ll_calc(jvar, evals);
    }

    gretl_matrix_free(S00);
    gretl_matrix_free(evals);

    if (!err) {
        err = j_general_restrict(jvar, rset, dset, flags, prn);
    }

    return err;
}

private int real_johansen_estimate(GretlVar * jvar,
                                   GretlRestriction * rset,
                                   const(DatasetRaw) * dset,
                                   int flags,
                                   Prn * prn)
{
    int err = 0;

    if (rset == null) {
        err = j_estimate_unrestr(jvar, dset);
    } else if (simple_beta_restriction(jvar, rset)) {
        err = est_simple_beta_restr(jvar, rset, dset, flags);
    } else if (simple_alpha_restriction(jvar, rset)) {
        err = est_simple_alpha_restr(jvar, rset, dset, flags, prn);
    } else {
        err = j_estimate_general(jvar, rset, dset, flags, prn);
    }

    return err;
}

/* Public entry point for VECM estimation, restricted or not */

extern (C) int johansen_estimate(GretlVar * jvar, GretlRestriction * rset,
                      const(DatasetRaw) * dset, Prn * prn)
{
    return real_johansen_estimate(jvar, rset, dset, 0, prn);
}

/* bootstrap round for "unrestricted" VECM: take some
   shortcuts
*/

private int boot_round_simple(GretlVar * jvar, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * M = null;
    GretlMatrixRaw * evals = null;
    int err;

    err = johansen_get_eigenvalues(jvar.jinfo.S00, jvar.jinfo.S01,
                                   jvar.jinfo.S11, &M, &evals,
                                   jrank(jvar));

    if (!err) {
        gretl_matrix_copy_values(jvar.jinfo.Beta, M);
        err = normalize_beta(jvar, null, null);
        if (!err) {
            err = VECM_estimate_full(jvar, null, dset,
                                     VecmFlags.bootstrapping);
        }
        if (!err) {
            err = compute_omega(jvar);
        }
    }

    gretl_matrix_free(M);
    gretl_matrix_free(evals);

    return err;
}

/* Simplified version of the Johansen procedure, to be called in
   the process of computing bootstrap confidence intervals for
   impulse response functions.  We just have to do enough to
   generate the VAR representation.
*/

extern (C) int johansen_boot_round(GretlVar * jvar, const(DatasetRaw) * dset)
{
    GretlRestriction * rset;
    int err = 0;

    rset = rset_from_VECM(jvar, &err);
    if (err) {
        return err;
    }

    if (rset != null) {
        /* experimental */
        err = real_johansen_estimate(jvar, rset, dset,
                                     VecmFlags.bootstrapping, null);
        free(rset);
    } else {
        err = boot_round_simple(jvar, dset);
    }

    return err;
}

void print_beta_alpha_Pi(const(GretlVar) * jvar,
                         const(DatasetRaw) * dset,
                         Prn * prn)
{
    int r = jrank(jvar);

    print_beta_or_alpha(jvar, r, dset, prn, VBeta.beta, 0);
    print_beta_or_alpha(jvar, r, dset, prn, VBeta.alpha, 0);
    pputc(prn, '\n');
    print_long_run_matrix(jvar, dset, prn);
}

/* compute and print beta, alpha and alpha*beta', in the context where
   we've tested a (common, homogeneous) restriction on beta,
   represented by H, and verbose output has been requested.
*/

private int show_beta_alpha_etc(const(GretlVar) * jvar,
                                const(GretlMatrixRaw) * H,
                                const(GretlMatrixRaw) * M,
                                const(DatasetRaw) * dset,
                                Prn * prn)
{
    JohansenInfo * jv = cast(JohansenInfo *) jvar.jinfo;
    int err = 0;

    gretl_matrix_multiply_mod(H, GretlMatrixMod.none,
                              M, GretlMatrixMod.none,
                              jv.Beta, GretlMatrixMod.none);

    if (jv.rank == 1) {
        /* and if r > 1? */
        double den = jv.Beta.val[0];

        if (!floateq(den, 0.0)) {
            gretl_matrix_divide_by_scalar(jv.Beta, den);
        }
    }

    if (!err) {
        err = compute_alpha(jv);
    }

    if (!err) {
        print_beta_alpha_Pi(jvar, dset, prn);
    }

    return err;
}

/* test for a common, homogeneous restriction on beta (only) */

private int vecm_beta_test(GretlVar * jvar,
                           GretlRestriction * rset,
                           const(DatasetRaw) * dset,
                           Gretlopt opt,
                           Prn * prn)
{
    const(GretlMatrixRaw) * R;
    GretlMatrixRaw * H = null;
    GretlMatrixRaw * M = null;
    GretlMatrixRaw * S11 = null;
    GretlMatrixRaw * S01 = null;
    GretlMatrixRaw * S00 = null;
    GretlMatrixRaw * evals = null;
    int verbose = (opt & Gretlopt.OPT_V) != 0;
    int s, p, r;
    int err = 0;

    R = rset_get_R_matrix(rset);
    H = johansen_nullspace(R, &err);

    if (err) {
        return err;
    }

    p = jvar.neqns;
    r = jrank(jvar);
    s = gretl_matrix_cols(H);

    S11 = gretl_matrix_alloc(s, s);
    S01 = gretl_matrix_alloc(p, s);
    S00 = gretl_matrix_copy(jvar.jinfo.S00);

    if (S11 == null || S01 == null || S00 == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    pprintf(prn, "\n%s\n\n",
            _("Test of restrictions on cointegrating relations").ptr);

    if (verbose) {
        gretl_matrix_print_to_prn(H, "Restriction matrix, H", prn);
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GretlMatrixMod.transpose,
                             jvar.jinfo.S11, S11,
                             GretlMatrixMod.none);

    if (verbose) {
        gretl_matrix_print_to_prn(S11, "H'*S11*H", prn);
    }

    if (!err) {
        /* S01 <- S01*H */
        err = gretl_matrix_multiply(jvar.jinfo.S01, H, S01);
    }

    if (verbose) {
        gretl_matrix_print_to_prn(S01, "S01*H", prn);
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00, S01, S11,
                                       &M, &evals, r);
    }

    if (!err) {
        if (verbose) {
            gretl_matrix_print_to_prn(M, "M", prn);
        }
        johansen_LR_calc(jvar, evals, H, rset, VBeta.beta, prn);
    }

    if (!err && verbose) {
        show_beta_alpha_etc(jvar, H, M, dset, prn);
    }

 bailout:

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);
    gretl_matrix_free(S00);
    gretl_matrix_free(S11);
    gretl_matrix_free(S01);

    return err;
}

/* Test of linear restrictions on the cointegrating relations in a
   VECM.  If the restrictions are "simple" (homogeneous and in common)
   we do the test using the eigen-system approach.  If they are
   "general" restrictions we hand off to the specialized machinery in
   jrestrict.c.
*/

extern (C) int vecm_test_restriction(GretlVar * jvar,
                          GretlRestriction * rset,
                          const(DatasetRaw) * dset,
                          Gretlopt opt,
                          Prn * prn)
{
    GretlMatrixRaw * B0 = null;
    GretlMatrixRaw * A0 = null;
    Prn * vprn;
    int err = 0;

    B0 = gretl_matrix_copy(jvar.jinfo.Beta);
    A0 = gretl_matrix_copy(jvar.jinfo.Alpha);

    if (B0 == null || A0 == null) {
        return GretlError.alloc;
    }

    vprn = (opt & Gretlopt.OPT_S)? null : prn;

    if (simple_beta_restriction(jvar, rset)) {
        err = vecm_beta_test(jvar, rset, dset, opt, vprn);
    } else if (simple_alpha_restriction(jvar, rset)) {
        err = vecm_alpha_test(jvar, rset, dset, opt, vprn);
    } else {
        err = general_vecm_analysis(jvar, rset, dset, vprn);
    }

    if (!err) {
        rset_record_LR_result(rset);
    }

    /* restore orginal Beta, Alpha on exit */

    gretl_matrix_replace(&jvar.jinfo.Beta, B0);
    gretl_matrix_replace(&jvar.jinfo.Alpha, A0);

    return err;
}
