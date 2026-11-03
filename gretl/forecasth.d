/* Port of the FITRESID struct in lib/src/forecast.h (the prototypes are for forecast.c, now forecast.d). */

module gretl.forecasth;


import gretl.libgretlh;
import gretl.printouth : SHORTLEN;

struct FitResid {
    int model_ID;   /* ID of model on which forecast is based */
    int asymp;      /* 0/1 flag for asymptotic estimator */
    int std;        /* 0/1 flag for standardized residuals */
    int model_t1;   /* start of model estimation range */
    FcastMethod method; /* one of the FcastMethod options */
    double * actual; /* array of values of dependent variable */
    double * fitted; /* array of fitted values */
    double * resid;  /* array of residuals */
    double * sderr;  /* array of forecast standard errors (or NULL) */
    double sigma;   /* standard error of regression */
    double alpha;   /* for confidence intervals */
    double a0;      /* correction factor for log-to-level forecasts */
    int a0meth;     /* method used for computing a0 (1 = normal, 2 = smear) */
    int pmax;       /* if positive, suggested number of decimal places
                       for use in printing */
    int df;         /* degrees of freedom for model */
    int t0;         /* start of pre-forecast data range */
    int t1;         /* start of forecast range */
    int t2;         /* end of forecast range */
    int k;          /* number of steps ahead (method = FC_KSTEP only) */
    int nobs;       /* length of the arrays actual, fitted, resid */
    char[VNAMELEN] depvar; /* full name of dependent variable */
    char[SHORTLEN] ylabel; /* name to display in printout or plot */
    Gretlopt opt;   /* option flags applied to forecast */
}
