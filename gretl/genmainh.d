/* Port of genmain.h */

module gretl.genmainh;


enum RetrievalIndex {
    nobs = 1,  /* number of observations in current sample range */
    nvars,  /* number of variables in dataset (including the constant) */
    pd,  /* periodicity of dataset */
    panel_pd,  /* pd in time dimension of panel */
    t1,  /* start of current sample range */
    t2,  /* end of current sample range */
    tmax,  /* maximum value for t2 */
    datatype,  /* dataset structure (x-section, time-series, panel) */
    windows,  /* running on MS Windows (1) or not (0) */
    version_,  /* gretl version number */
    errno,  /* internal gretl error code */
    seed,  /* RNG seed */
    huge,  /* conventional "huge" number, like eg 1.0e100 */
    dset_max,  /* separator */
    test_lnl,  /* log-likelihood from last test (if applicable) */
    stopwatch,  /* stopwatch */
    test_brk,  /* obs at which break occurs (QLR test) */
    loglevel,  /* current loglevel for errors/warnings/etc */
    logstamp,  /* logger shows timestamp? (0/1) */
    scalar_max,  /* separator: scalars vs series */
    index,  /* 1-based observations index */
    time,  /* 1-based time index */
    punit,  /* 1-based panel unit index */
    obsmaj,  /* major component of observation (e.g. year) */
    obsmin,  /* minor component of observation (e.g. quarter, month) */
    obsmic,  /* micro component of observation (e.g. day) */
    dates,  /* ISO 8601 "basic" dates series */
    series_max,  /* separator: series vs variants */
    test_stat,  /* last test statistic(s) (scalar or matrix) */
    test_pval,  /* last test p-value(s) (scalar or matrix) */
    now,  /* current date/time (matrix) */
    result,  /* result of a "result-compatible" command */
    grpnames,  /* panel group names (array of strings) */
    pngfont,  /* name of font selected for plots */
    mapfile,  /* name of current map file, if any */
    map,  /* current map (if any) as bundle */
    max
}

/* Symbolic names for variables that can be retrieved
   from gretl models. */
enum ModelDataIndex {
    ess = RetrievalIndex.max + 1,  /* error sum of squares */
    t,  /* observations used */
    rsq,  /* R-squared */
    sigma,  /* standard error of residuals */
    df,  /* degrees of freedom */
    ncoeff,  /* total number of estimated coefficients */
    lnl,  /* log-likelihood */
    gmmcrit,  /* GMM criterion */
    aic,  /* Akaike info criterion */
    bic,  /* Bayesian info criterion */
    hqc,  /* Hannan-Quinn criterion */
    trsq,  /* T * R-squared, last model */
    dw,  /* Durbin-Watson statistic */
    dwpval,  /* Durbin-Watson p-value, last model */
    fstt,  /* overall F-statistic, last model */
    chisq,  /* overall chi-square stat, last model */
    diagtest,  /* system test for diagonal covariance matrix */
    diagpval,  /* p-value for the above */
    scalar_max,  /*** SEPARATOR, scalars/series ***/
    uhat,  /* residuals */
    yhat,  /* fitted values */
    llt,  /* per-observation loglikelihood */
    ahat,  /* individual effects in panel model */
    h,  /* GARCH predicted variances */
    sample,  /* observations used in estimation */
    uhat2,  /* squared residuals */
    series_max,  /*** SEPARATOR, series/matrices ***/
    coeff,  /* parameter estimates */
    se,  /* parameter standard errors */
    vcv,  /* parameter covariance matrix */
    rho,  /* autoregressive coefficients */
    compan,  /* VAR companion matrix */
    xtxinv,  /* VARs, VECMs: X'X^{-1} */
    vecg,  /* VECMs: the Gamma matrices */
    evals,  /* VECMs: eigenvalues */
    jalpha,  /* Johansen's alpha */
    jbeta,  /* Johansen's beta */
    jvbeta,  /* Covariance matrix for Johansen's normalized beta */
    js00,  /* VECM residual covariance matrix (1st differences) */
    js11,  /* VECM residual covariance matrix (levels) */
    js01,  /* VECM residual cross-product matrix */
    hausman,  /* Hausman test after tsls or fixed effects */
    sargan,  /* Sargan over-identification test after tsls */
    sysgam,  /* Parameter matrix Gamma (simultaneous systems) */
    sysa,  /* Parameter matrix A (simultaneous systems) */
    sysb,  /* Parameter matrix B (simultaneous systems) */
    fcast,  /* last forecast generated via fcast command */
    fcse,  /* standard errors associated with M_FCAST */
    coeff_ci,  /* (asymmetric) confidence intervals for coeffs */
    ehat,  /* ARMA: vector of estimated innovations */
    pmanteau,  /* VAR portmanteau test plus 's' value */
    oddsratios,  /* binary logit odds ratios */
    matrix_max,  /*** SEPARATOR, end of matrices ***/
    ec,  /* VECM error-correction terms */
    vma,  /* VARs, VECMs: vector moving average representation */
    fevd,  /* VAR variance decomposition */
    mnlprobs,  /* case probabilities for multinomial logit (legacy) */
    allprobs,  /* outcome probabilities for ordered and multinomial models */
    mbuild_max,  /*** SEPARATOR, end of matrix-builders ***/
    xlist,  /* list of regressors */
    ylist,  /* list of endogenous variables */
    list_max,  /*** SEPARATOR, end of lists ***/
    command,  /* model command word */
    depvar,  /* name of dependent variable */
    parnames,  /* array of parameter names */
    max  /* sentinel */
}

enum BundleDataIndex {
    model = ModelDataIndex.max + 1,  /* last model as bundle */
    system,  /* last VAR/VECM/system as bundle */
    sysinfo  /* system information */
}

bool model_data_scalar(int i) {
    return (i > RetrievalIndex.max && i < ModelDataIndex.scalar_max);
}

bool model_data_series(int i) {
    return (i > ModelDataIndex.scalar_max && i < ModelDataIndex.series_max);
}

bool model_data_matrix(int i) {
    return (i > ModelDataIndex.series_max && i < ModelDataIndex.matrix_max);
}

bool model_data_matrix_builder(int i) {
    return (i > ModelDataIndex.matrix_max && i < ModelDataIndex.mbuild_max);
}

bool model_data_list(int i) {
    return (i > ModelDataIndex.mbuild_max && i < ModelDataIndex.list_max);
}

bool model_data_string(int i) {
    return (i > ModelDataIndex.list_max && i < ModelDataIndex.parnames);
}

bool model_data_array(int i) {
    return (i == ModelDataIndex.parnames);
}

/* "typedef struct parser_ GENERATOR": the public name for the parser struct (defined in genparseh.d) */
import gretl.genparseh : Parser;

alias Generator = Parser;
