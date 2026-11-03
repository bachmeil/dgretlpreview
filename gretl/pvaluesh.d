/* Port of pvalues.h */

module gretl.pvaluesh;


enum DistCode {
    none = 0,
    uniform,
    udiscrt,
    normal,
    student,
    chisq,
    snedecor,
    binomial,
    poisson,
    expon,
    weibull,
    gamma,
    ged,
    laplace,
    beta,
    dw,
    binorm,
    johansen,
    betabin,
    nc_chisq,
    nc_f,
    nc_t,
    logistic,
    dirichlet,
    discrete,
    tnormal
}
