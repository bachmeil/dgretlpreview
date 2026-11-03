/* Port of lib/src/johansen.h: the Johansen (VECM) information struct and the macros that inspect it. The
   prototypes in the C header belong to johansen.c/var.c. The enumeration prefix J_ is dropped, and the macros jcode,
   jrank and auto_restr are functions. */

module gretl.johansenh;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.varh : GretlVar;

enum JohansenCode {
    no_const = 0,
    rest_const,
    unrest_const,
    rest_trend,
    unrest_trend
}

struct JohansenInfo {
    int ID;               /* for identifying saved vars */
    JohansenCode code;    /* see enumeration above */
    int rank;             /* if specified, chosen cointegration rank, else 0 */
    int seasonals;        /* number of seasonal dummies included */
    GretlMatrixRaw * R0;     /* residuals, VAR in differences */
    GretlMatrixRaw * R1;     /* residuals, second regressions */
    GretlMatrixRaw * S00;    /* cross-products of residuals */
    GretlMatrixRaw * S11;    /* cross-products of residuals */
    GretlMatrixRaw * S01;    /* cross-products of residuals */
    GretlMatrixRaw * evals;  /* vector of eigenvalues */
    GretlMatrixRaw * Beta;   /* matrix of eigenvectors */
    GretlMatrixRaw * Alpha;  /* matrix of adjustments */
    GretlMatrixRaw * Gamma;  /* I - sum of \Gamma_i matrices */
    GretlMatrixRaw * JC;     /* Johansen's 'C' (Lütkepohl's $\Xi$) */
    GretlMatrixRaw * JVC;    /* Variance of vec(C) */
    GretlMatrixRaw * Cse;    /* Standard errors of C */
    GretlMatrixRaw * Bvar;   /* variance matrix of beta */
    GretlMatrixRaw * Bse;    /* standard errors of beta */
    GretlMatrixRaw * Ase;    /* standard errors of alpha */
    GretlMatrixRaw * R;      /* beta-restriction LHS matrix */
    GretlMatrixRaw * q;      /* beta-restrictions RHS matrix */
    GretlMatrixRaw * Ra;     /* alpha-restriction LHS matrix */
    GretlMatrixRaw * qa;     /* alpha-restrictions RHS matrix */
    GretlMatrixRaw * YY;     /* double-size Y matrix */
    GretlMatrixRaw * RR;     /* double-size residuals matrix */
    GretlMatrixRaw * BB;     /* double-size coefficient matrix */
    GretlArray * G;       /* array of \Gamma_i matrices */
    double ll0;           /* unrestricted log-likelihood */
    int lrdf;             /* df for likelihood ratio test */
    double prior_ll;      /* ll for prior model in restriction sequence */
    int prior_df;         /* df for prior model in restriction sequence */
}

JohansenCode jcode(const(GretlVar) * v)
{
    return (v.jinfo == null)? JohansenCode.no_const : v.jinfo.code;
}

int jrank(const(GretlVar) * v)
{
    return (v.jinfo == null)? 0 : v.jinfo.rank;
}

/* vecm contains an "automatic" restricted term */
bool auto_restr(const(GretlVar) * v)
{
    return v.jinfo != null && (v.jinfo.code == JohansenCode.rest_const || v.jinfo.code == JohansenCode.rest_trend);
}
