/* Port of gretl_model.h */

module gretl.gretl_modelh;


import gretl.gretl_commandsh;
import gretl.gretl_model : gretl_model_get_data, gretl_model_get_int;
import gretl.libgretlh;

enum ArmaFlags {
    x12a  = 1 << 0, /* using X-12-ARIMA (or X-13) to generate estimates */
    exact = 1 << 1, /* using exact ML */
    ls    = 1 << 2, /* using conditional ML, and O/NLS == CML */
    ols   = 1 << 3  /* OLS == MLE */
}

struct CoeffIntervals {
    int asy;
    int ncoeff;
    double alpha;
    double t;
    char ** names;
    double * coeff;
    double * maxerr;
    int df;
    Gretlopt opt;
}

enum VCVMajorType {
    classical,
    hc,
    hac,
    ml,
    panel,
    rq,
    cluster
}

enum MLVCVType {
    unset,
    hessian,
    im,
    op,
    qml,
    bw,
    hac,
    vcvmax
}

enum HACKernel {
    bartlett,
    parzen,
    qs,
    max
}

enum PanelVCVType {
    hac,  /* clustered by individual (Arellano) */
    bk,   /* Beck-Katz PCSE */
    time, /* clustered by period */
    dk,   /* Driscoll-Kraay SCC */
    both  /* clustered by both unit and period */
}

enum RQVCVType {
    asy,  /* asymptotic */
    nid   /* sandwich */
}

enum VCVFlags {
    prewhiten = 1
}

struct VCVInfo {
    int vmaj;        /* general type of VCV (see VCVMajorType) */
    int vmin;        /* variant of general type */
    int order;       /* for use with HAC */
    VCVFlags flags;  /* includes prewhitening */
    double bw;       /* for use with QS HAC kernel */
    char * cv1;      /* name of (first) cluster variable */
    char * cv2;      /* name of second cluster var */
}

/* single-equation model commands */

bool MODEL_COMMAND(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.ar ||
            c == GretlCmdIndex.ar1 ||
            c == GretlCmdIndex.arch ||
            c == GretlCmdIndex.arma ||
            c == GretlCmdIndex.dpanel ||
            c == GretlCmdIndex.duration ||
            c == GretlCmdIndex.garch ||
            c == GretlCmdIndex.gmm ||
            c == GretlCmdIndex.heckit ||
            c == GretlCmdIndex.hsk ||
            c == GretlCmdIndex.intreg ||
            c == GretlCmdIndex.ivreg ||
            c == GretlCmdIndex.lad ||
            c == GretlCmdIndex.logistic ||
            c == GretlCmdIndex.logit ||
            c == GretlCmdIndex.midasreg ||
            c == GretlCmdIndex.mle ||
            c == GretlCmdIndex.mpols ||
            c == GretlCmdIndex.negbin ||
            c == GretlCmdIndex.nls ||
            c == GretlCmdIndex.ols ||
            c == GretlCmdIndex.panel ||
            c == GretlCmdIndex.poisson ||
            c == GretlCmdIndex.probit ||
            c == GretlCmdIndex.biprobit ||
            c == GretlCmdIndex.quantreg ||
            c == GretlCmdIndex.tobit ||
            c == GretlCmdIndex.wls);
}

bool AR_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.ar ||
            c == GretlCmdIndex.ar1 ||
            c == GretlCmdIndex.arma ||
            c == GretlCmdIndex.garch);
}

bool SIMPLE_AR_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.ar ||
            c == GretlCmdIndex.ar1);
}

bool ML_ESTIMATOR(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.arma ||
            c == GretlCmdIndex.duration ||
            c == GretlCmdIndex.garch ||
            c == GretlCmdIndex.heckit ||
            c == GretlCmdIndex.logit ||
            c == GretlCmdIndex.mle ||
            c == GretlCmdIndex.negbin ||
            c == GretlCmdIndex.poisson ||
            c == GretlCmdIndex.probit ||
            c == GretlCmdIndex.biprobit ||
            c == GretlCmdIndex.tobit);
}

bool LIMDEP(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.logit ||
            c == GretlCmdIndex.probit ||
            c == GretlCmdIndex.tobit ||
            c == GretlCmdIndex.intreg);
}

bool COUNT_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.poisson ||
            c == GretlCmdIndex.negbin);
}

bool LSQ_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.ar1 ||
            c == GretlCmdIndex.hsk ||
            c == GretlCmdIndex.ols ||
            c == GretlCmdIndex.wls);
}

bool ASYMPTOTIC_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.arma ||
            c == GretlCmdIndex.dpanel ||
            c == GretlCmdIndex.duration ||
            c == GretlCmdIndex.garch ||
            c == GretlCmdIndex.gmm ||
            c == GretlCmdIndex.heckit ||
            c == GretlCmdIndex.intreg ||
            c == GretlCmdIndex.ivreg ||
            c == GretlCmdIndex.logit ||
            c == GretlCmdIndex.mle ||
            c == GretlCmdIndex.negbin ||
            c == GretlCmdIndex.poisson ||
            c == GretlCmdIndex.probit ||
            c == GretlCmdIndex.tobit ||
            c == GretlCmdIndex.biprobit);
}

bool EQN_SYSTEM_COMMAND(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.var ||
            c == GretlCmdIndex.vecm ||
            c == GretlCmdIndex.system);
}

/* model where the specification is not based on a list
   of variables */
bool NONLIST_MODEL(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.nls ||
            c == GretlCmdIndex.mle ||
            c == GretlCmdIndex.gmm ||
            c == GretlCmdIndex.midasreg);
}

bool is_model_ref_cmd(GretlCmdIndex c)
{
    return (c == GretlCmdIndex.add ||
            c == GretlCmdIndex.arch ||
            c == GretlCmdIndex.chow ||
            c == GretlCmdIndex.cusum ||
            c == GretlCmdIndex.fcast ||
            c == GretlCmdIndex.leverage ||
            c == GretlCmdIndex.modtest ||
            c == GretlCmdIndex.omit ||
            c == GretlCmdIndex.restrict ||
            c == GretlCmdIndex.vif);
}

bool RQ_SPECIAL_MODEL(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.lad || m.ci == GretlCmdIndex.quantreg) &&
            null != gretl_model_get_data(m, "rq_tauvec"));
}

bool POOLED_MODEL(const(Model) * m)
{
    return ((m.ci == GretlCmdIndex.ols || m.ci == GretlCmdIndex.panel) &&
            gretl_model_get_int(m, "pooled"));
}

enum ModelTestType {
    add,
    arch,
    autocorr,
    chow,
    cusum,
    qlr,
    groupwise,
    logs,
    normal,
    omit,
    reset,
    squares,
    whites,
    sargan,
    iv_hausman,
    panel_hausman,
    panel_f,
    panel_bp,
    panel_timedum,
    panel_ar,
    het_1,
    bp,
    chowdum,
    comfac,
    indep,
    re,
    within_f,
    panel_welch,
    re_wald,
    xdepend,
    bds,
    max
}
