/* Port of compare.h */

module gretl.compareh;


enum SeriesSaveCode {
    leverage  = 1 << 0,
    influence = 1 << 1,
    dffits    = 1 << 2,
    trend     = 1 << 3,
    cycle     = 1 << 4
}

/* Symbolic names to keep track of auxiliary regression models,
   which are estimated either for the purpose of carrying out
   some sort of diagnostic test or which form part of a
   multi-equation system. */
enum ModelAuxCode {
    none,       /* not an auxiliary regression */
    sq,         /* nonlinearity test (squared terms) */
    log,        /* nonlinearity test (log terms) */
    chow,       /* Chow test */
    add,        /* LM test regression for added variables */
    ar,         /* autocorrelation test */
    arch,       /* ARCH test */
    white,      /* heteroskedasticity (White's test) */
    coint,      /* cointegration test */
    df,         /* Dickey-Fuller test */
    adf,        /* augmented Dickey-Fuller test */
    kpss,       /* KPSS unit-root test */
    omit,       /* unused */
    reset,      /* Ramsey's RESET */
    sys,        /* single equation from multivariate system */
    var,        /* single equation from VAR system */
    vecm,       /* single equation from VECM system */
    johansen,   /* Johansen cointegration test */
    groupwise,  /* test for groupwise heteroskedasticity */
    het_1,      /* Pesaran-Taylor HET_1 test */
    bp,         /* Breusch-Pagan heteroskedastcity test */
    aux,        /* auxiliary regression not otherwise specified */
    comfac,     /* common factor test */
    biprob      /* biprobit initializer */
}
