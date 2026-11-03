/* Port of lib/src/nlspec.h: the private header shared by nls.c and gmm.c, with the nonlinear-model specification
   struct. The C type `nlspec` (opaque in nls.h) is `Nlspec`; `parm` (defined in nls.c) is declared here with
   the struct that holds an array of them, and `ocset` (gmm.c) is opaque. The NL_ prefix of the flags is
   dropped, and the macros of nls.c that inspect a spec (scalar_param, matrix_deriv, numeric_mode, ...) are
   functions here. The prototypes belong to nls.c and gmm.c. */

module gretl.nlspech;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.genparseh : Parser;
import gretl.gretl_prn : Prn;
import gretl.libgretlh;

enum NlFlags {
    analytical  = 1 << 0,
    autoreg     = 1 << 1,
    ahess       = 1 << 2,
    newton      = 1 << 3,
    smallstep   = 1 << 4,
    names_array = 1 << 5
}

struct Parm {
    char[VNAMELEN] name;   /* name of parameter */
    GretlBundle * bundle;  /* parent bundle, if applicable */
    GretlType type;        /* type of parameter (scalar or matrix) */
    int dvtype;            /* type of derivative for parameter */
    char * deriv;          /* string representation of derivative of regression
                              function with respect to param (or NULL) */
    int dnum;              /* ID number of series holding the derivative, or 0 */
    int nc;                /* number of individual coefficients associated
                              with the parameter */
    char[VNAMELEN] dname;  /* name of variable holding the derivative */
    Parser * dgenr;     /* generator for derivative */
    GretlMatrixRaw * mat;     /* pointer to matrix parameter */
}

struct Ocset;   /* opaque: orthogonality information, defined in gmm.c */

struct Nlspec {
    int ci;             /* NLS, MLE or GMM */
    int generr;         /* error from genr */
    NlFlags flags;      /* numeric or analytic derivatives, etc. */
    Gretlopt opt;       /* can include OPT_V for verbose output; if ci = MLE
                           can also include OPT_H (Hessian) or OPT_R (QML)
                           to control the estimator of the variance matrix;
                           can also include OPT_N to force use of
                           numerical derivatives. */
    int dv;             /* ID number of dependent variable (NLS) */
    int lhtype;         /* type of the LHS variable */
    char[VNAMELEN] lhname; /* name of LHS var in criterion function */
    char[VNAMELEN] hname;  /* name of Hessian matrix, if present */
    char * parnames;    /* user-set names for parameters */
    int lhv;            /* ID number of LHS series in function being
                           minimized or maximized... */
    GretlMatrixRaw * lvec; /* or LHS vector */
    char * nlfunc;      /* string representation of function,
                           expressed in terms of the residuals (NLS,
                           GMM) or the log-likelihood (MLE) */
    int nparam;         /* number of parameters */
    int ncoeff;         /* number of coefficients (allows for vector params) */
    int nmat;           /* number of matrix parameters */
    int naux;           /* number of auxiliary commands */
    int ngenrs;         /* number of variable-generating formulae */
    int iters;          /* number of iterations performed */
    int fncount;        /* number of function evaluations (BFGS) */
    int grcount;        /* number of gradient evaluations (BFGS) */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int real_t1;        /* real starting observation (if sub-sampled) */
    int real_t2;        /* real ending observation (if sub-sampled) */
    int nobs;           /* number of observations used */
    double crit;        /* criterion (minimand or maximand) */
    double tol;         /* tolerance for stopping iteration */
    Parm * params;      /* array of information on function parameters */
    double * fvec;      /* function vector */
    double * coeff;     /* coefficient estimates */
    GretlMatrixRaw * J;    /* Jacobian matrix */
    GretlMatrixRaw * Hinv; /* negative inverse of Hessian */
    GretlMatrixRaw * bounds; /* param bounds, specific to GMM with --lbfgs */
    char ** aux;        /* auxiliary commands */
    char * hesscall;    /* function call for Hessian */
    Parser ** genrs; /* variable-generation pointers */
    Parser * hgen;   /* generator for Hessian */
    DatasetRaw * dset;     /* pointer to dataset */
    Prn * prn;          /* printing aparatus */
    Ocset * oc;         /* orthogonality info (GMM) */
    char * missmask;    /* mask for missing observations */
}

bool scalar_param(const(Nlspec) * s, int i)
{
    return s.params[i].type == GretlType.double_;
}

bool matrix_deriv(const(Nlspec) * s, int i)
{
    return s.params[i].dvtype == GretlType.matrix;
}

bool scalar_deriv(const(Nlspec) * s, int i)
{
    return s.params[i].dvtype == GretlType.double_;
}

bool numeric_mode(const(Nlspec) * s)
{
    return !(s.flags & NlFlags.analytical);
}

bool analytic_mode(const(Nlspec) * s)
{
    return (s.flags & NlFlags.analytical) != 0;
}

bool scalar_loglik(const(Nlspec) * s)
{
    return s.lhtype == GretlType.double_;
}

bool suppress_grad_check(const(Nlspec) * s)
{
    return (s.opt & Gretlopt.OPT_S) != 0;
}
