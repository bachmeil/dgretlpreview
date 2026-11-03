/* Port of the types in lib/src/system.h (the prototypes belong to system.d). The members of the GretlSystemMethods enum have the
   SYS_METHOD_ prefix dropped; "3sls" is not an identifier, so it is `sls3`. The anonymous flags enum is
   SystemFlags. */

module gretl.systemh;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;

enum GretlSystemMethods {
    sur = 0,
    sls3,
    fiml,
    liml,
    ols,
    tsls,
    wls,
    max
}

enum SystemFlags {
    /* note: allow for obsoleted flags in saved sessions --
       treat this enumeration as append-only */
    dfcorr      = 1 << 2,
    vcv_geomean = 1 << 3,
    restrict    = 1 << 4,
    iterate     = 1 << 5,
    saveit      = 1 << 6,
    liml1       = 1 << 7,
    quiet       = 1 << 8,
    robust      = 1 << 9
}

/* the helper structs below are defined in system.c in the C (opaque in system.h); they are declared here because
   EquationSystem refers to them and system.d needs the layouts */
struct IdAtom {
    int op;         /* operator (plus or minus) */
    int varnum;     /* ID number of variable to right of operator */
}

struct Identity {
    int n_atoms;    /* number of "atomic" elements in identity */
    int depvar;     /* LHS variable in indentity */
    IdAtom * atoms; /* pointer to RHS "atoms" */
}

struct Predet {
    int id;         /* ID number of lag variable */
    int src;        /* ID number of "parent" */
    int lag;        /* lag order */
}

struct LimlData {
    double * lmin;  /* min. eigenvalues */
    double * ll;    /* per-equation log likelihood */
    int * idf;      /* overidentification df */
}

struct EquationSystem {
    char * name;                /* user-specified name for system, or NULL */
    int refcount;               /* for saving/deleting */
    int fd;                     /* function execution depth */
    int t1;                     /* starting observation number */
    int t2;                     /* ending observation number */
    int smpl_t1;                /* first obs in sample range */
    int smpl_t2;                /* last obs in sample range */
    int T;                      /* number of observations per equation */
    int df;                     /* T - average coeffs per equation */
    int method;                 /* estimation method */
    int neqns;                  /* number of stochastic equations */
    int nidents;                /* number of identities */
    int order;                  /* max lag of endogenous variable */
    int iters;                  /* number of iterations taken */
    int flags;                  /* to record options */
    double ll;                  /* log-likelihood (restricted) */
    double llu;                 /* unrestricted log-likelihood */
    double X2;                  /* chi-square test value */
    double ess;                 /* total error sum of squares */
    double diag_test;           /* test stat for diagonal covariance matrix */
    double bdiff;               /* summary stat for change in coefficients */
    double ldet;                /* log-determinant of covariance matrix */
    int ** lists;               /* regression lists for stochastic equations */
    int ** tslists;             /* back-up of TSLS-style lists */
    int * ylist;                /* list of endogenous variables */
    int * ilist;                /* list of instruments */
    int * xlist;                /* list of truly exogenous variables */
    int * plist;                /* list of predetermined variables */
    int * biglist;              /* list of all variables, for data checking */
    Predet * pre_vars;          /* array of info on predetermined regressors */
    Identity ** idents;         /* set of identities */
    GretlMatrixRaw * b;            /* coefficient estimates */
    GretlMatrixRaw * vcv;          /* covariance matrix of coefficients */
    GretlMatrixRaw * S;            /* cross-equation covariance matrix */
    GretlMatrixRaw * R;            /* LHS of any linear restrictions */
    GretlMatrixRaw * q;            /* RHS of any linear restrictions */
    GretlMatrixRaw * E;            /* residuals, all equations */
    GretlMatrixRaw * yhat;         /* fitted values, all equations */
    GretlMatrixRaw * Gamma;        /* structural form Gamma matrix (endog + identities)*/
    GretlMatrixRaw * B;            /* structural form B matrix (exogenous) */
    GretlMatrixRaw * A;            /* structural form A matrix (lagged endogenous) */
    GretlMatrixRaw * F;            /* forecast matrix */
    GretlMatrixRaw * Sr;           /* reduced-form error covariance matrix */
    Model ** models;            /* set of pointers to per-equation models */
    LimlData * ldata;           /* extra info from LIML estimation */
}
