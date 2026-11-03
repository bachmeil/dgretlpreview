/* Port of the types in lib/src/var.h (the prototypes are for var.c, which is not ported yet; the C functions are
   declared extern (C) where used). The macros var_n_lags, var_max_lag and levels_order are functions. */

module gretl.varh;


import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.johansenh : JohansenInfo;
import gretl.libgretlh;

enum DetFlags {
    const_ = 1 << 0,
    trend  = 1 << 1,
    seas   = 1 << 2
}

enum VarRobust {
    hc = 1,
    hac
}

struct GretlVar {
    int ci;              /* command index (VAR or VECM) */
    int refcount;        /* for saving/deleting */
    int err;             /* error code */
    int neqns;           /* number of equations in system */
    int order;           /* maximum lag order */
    int t1;              /* starting observation */
    int t2;              /* ending observation */
    int T;               /* number of observations */
    int df;              /* T - average coeffs per equation */
    int ifc;             /* equations include a constant (1) or not (0) */
    int ncoeff;          /* total coefficients per equation */
    int * lags;          /* list of specific lags */
    int * ylist;         /* list of stochastic vars */
    int * xlist;         /* list of exogenous variables */
    int * rlist;         /* restricted exogenous variables (VECM only) */
    int detflags;        /* record of automatic deterministic vars added */
    int robust;          /* computing robust std errors? */
    int xcols;           /* full column size of X matrix (VECM special) */
    GretlMatrixRaw * Y;     /* matrix of dependent variables */
    GretlMatrixRaw * X;     /* matrix of independent variables */
    GretlMatrixRaw * B;     /* basic coefficient matrix */
    GretlMatrixRaw * XTX;   /* X'X inverse */
    GretlMatrixRaw * A;     /* augmented coefficient matrix (companion form) */
    GretlMatrixRaw * L;     /* lambda: inverse roots of A(L) polynomial */
    GretlMatrixRaw * E;     /* residuals matrix */
    GretlMatrixRaw * C;     /* Cholesky-decomposed covariance matrix */
    GretlMatrixRaw * S;     /* cross-equation variance matrix */
    GretlMatrixRaw * F;     /* optional forecast matrix */
    GretlMatrixRaw * V;     /* full parameter covariance matrix */
    GretlMatrixRaw * ord;   /* optional Cholesky-ordering vector */
    Model ** models;     /* pointers to individual equation estimates */
    double * Fvals;      /* hold results of F-tests */
    double * Ivals;      /* hold results of info criteria comparisons */
    double ldet;         /* log-determinant of S */
    double ll;           /* log-likelihood */
    double AIC;          /* Akaike criterion */
    double BIC;          /* Bayesian criterion */
    double HQC;          /* Hannan-Quinn criterion */
    double LR;           /* for likelihood-ratio testing */
    double LB;           /* Ljung-Box (Portmanteau) test statistic */
    int LBs;             /* order for Portmanteau test */
    JohansenInfo * jinfo; /* extra information for VECMs */
    char * name;         /* for use in session management */
}

int var_n_lags(const(GretlVar) * v)
{
    return (v.lags != null)? v.lags[0] : v.order;
}

int var_max_lag(const(GretlVar) * v)
{
    return (v.lags != null)? v.lags[v.lags[0]] : v.order;
}

int levels_order(const(GretlVar) * v)
{
    return v.order + (v.ci == GretlCmdIndex.vecm);
}
