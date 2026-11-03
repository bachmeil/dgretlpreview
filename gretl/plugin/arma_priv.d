/* Port of plugin/arma_priv.h: the private types, flags and helpers shared by the files of the ARMA plugin (arma.c, arma_common.c,
   arma_init.c, arma_cml.c, kalman_arma.c).

   In the C the macros arma_has_seasonal() and so on are macros on the struct; here they are small functions. The flag fields
   are plain ints holding PrivFlags / ArmaFlags values (as in the Kalman struct), so that the bit operations need no casts.

   The functions the header declares are defined, extern (C), in gretl.plugin.arma (maybe_correct_MA, arma_model_add_roots,
   write_arma_model_stats, arma_list_y_position, arima_difference, arima_difference_undo), in gretl.plugin.arma_init (flip_poly,
   maybe_set_yscale, hr_arma_init, ar_arma_init, arma_by_simple_ols, arma_by_ls, transform_arma_const) and in
   gretl.plugin.arma_cml (bhhh_arma). */

module gretl.plugin.arma_priv;

import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.gretl_modelh : ArmaFlags;
import gretl.gretl_prn : Prn;
import gretl.libgretlh : DatasetRaw, Gretlopt, Model;

enum PrivFlags {
    seas  = 1 << 0, /* includes seasonal component */
    dspec = 1 << 1, /* input list includes differences */
    xdiff = 1 << 2, /* ARIMA: exogenous regressors are differenced */
    lbfgs = 1 << 3, /* using L-BFGS-B with native exact ML */
    vech  = 1 << 4, /* using vech representation when computing variance matrix of state for Kalman filter */
    naok  = 1 << 5, /* allow missing observations */
    nas   = 1 << 6, /* sample contains NAs */
    lev   = 1 << 7, /* doing ARIMA via levels formulation */
    ydiff = 1 << 8, /* ainfo->y contains differenced y */
    avgll = 1 << 9, /* passing average likelihood option to Kalman */
    stdx  = 1 << 10 /* exogenous regressors standardized */
}

enum IniMethod {
    user = 1,
    hr,
    small,
    nls,
    ols
}

struct ArmaInfo {
    int yno;            /* ID of dependent variable */
    int flags;          /* specification flags (ArmaFlags) */
    int pflags;         /* "private" flags for estimation (PrivFlags) */
    int * alist;        /* copy of incoming list */
    const(int) * pqspec;/* auxiliary list with specific AR, MA lags */
    char * pmask;       /* specific AR lags included */
    char * qmask;       /* specific MA lags included */
    double ll;          /* log-likelihood */
    IniMethod init_;    /* initialization method (C name: init) */
    int ifc;            /* specification includes a constant? */
    int p;              /* max non-seasonal AR order */
    int d;              /* non-seasonal difference */
    int q;              /* max non-seasonal MA order */
    int P;              /* seasonal AR order */
    int D;              /* seasonal difference */
    int Q;              /* seasonal MA order */
    int np;             /* total non-seasonal AR lags */
    int nq;             /* total non-seasonal MA lags */
    int maxlag;         /* longest lag in model */
    int nexo;           /* number of other regressors (ARMAX) */
    int nc;             /* total number of coefficients */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int pd;             /* periodicity of data */
    int T;              /* number of valid observations in sample */
    int fullT;          /* total obs (possibly including interior NAs) */
    int r0;             /* size of ARMA state vector */
    int fncount;        /* count of function evaluations */
    int grcount;        /* count of gradient evaluations */
    double * y;         /* dependent variable (possibly differenced) */
    double * e;         /* forecast errors */
    const(double) ** Z; /* virtual dataset */
    double yscale;      /* scale factor for y */
    double yshift;      /* shift factor for y */
    int * xlist;        /* list of regressors (ARMAX) */
    int * misslist;     /* list of missing observations */
    GretlMatrixRaw * xstats; /* mean and std dev of regressors */
    GretlMatrixRaw * dX;     /* differenced regressors (ARIMAX) */
    GretlMatrixRaw * G;      /* score matrix */
    GretlMatrixRaw * V;      /* covariance matrix */
    int n_aux;          /* number of auxiliary arrays */
    double ** aux;      /* auxiliary arrays */
    Prn * prn;          /* verbose printer */
}

bool arma_by_x12a(const(ArmaInfo) * a)       { return (a.flags & ArmaFlags.x12a) != 0; }
bool arma_exact_ml(const(ArmaInfo) * a)      { return (a.flags & ArmaFlags.exact) != 0; }
bool arma_least_squares(const(ArmaInfo) * a) { return (a.flags & ArmaFlags.ls) != 0; }

void set_arma_least_squares(ArmaInfo * a)    { a.flags |= ArmaFlags.ls; }

bool arma_has_seasonal(const(ArmaInfo) * a)  { return (a.pflags & PrivFlags.seas) != 0; }
bool arma_is_arima(const(ArmaInfo) * a)      { return (a.pflags & PrivFlags.dspec) != 0; }
bool arma_xdiff(const(ArmaInfo) * a)         { return (a.pflags & PrivFlags.xdiff) != 0; }
bool arma_lbfgs(const(ArmaInfo) * a)         { return (a.pflags & PrivFlags.lbfgs) != 0; }
bool arma_using_vech(const(ArmaInfo) * a)    { return (a.pflags & PrivFlags.vech) != 0; }
bool arma_na_ok(const(ArmaInfo) * a)         { return (a.pflags & PrivFlags.naok) != 0; }
bool arma_missvals(const(ArmaInfo) * a)      { return (a.pflags & PrivFlags.nas) != 0; }
bool arima_levels(const(ArmaInfo) * a)       { return (a.pflags & PrivFlags.lev) != 0; }
bool arima_ydiff(const(ArmaInfo) * a)        { return (a.pflags & PrivFlags.ydiff) != 0; }
bool arma_avg_ll(const(ArmaInfo) * a)        { return (a.pflags & PrivFlags.avgll) != 0; }
bool arma_stdx(const(ArmaInfo) * a)          { return (a.pflags & PrivFlags.stdx) != 0; }

void set_arma_has_seasonal(ArmaInfo * a)  { a.pflags |= PrivFlags.seas; }
void set_arma_is_arima(ArmaInfo * a)      { a.pflags |= PrivFlags.dspec; }
void unset_arma_is_arima(ArmaInfo * a)    { a.pflags &= ~PrivFlags.dspec; }
void set_arma_use_vech(ArmaInfo * a)      { a.pflags |= PrivFlags.vech; }
void set_arma_na_ok(ArmaInfo * a)         { a.pflags |= PrivFlags.naok; }
void set_arma_missvals(ArmaInfo * a)      { a.pflags |= PrivFlags.nas; }
void set_arima_levels(ArmaInfo * a)       { a.pflags |= PrivFlags.lev; }
void set_arma_avg_ll(ArmaInfo * a)        { a.pflags |= PrivFlags.avgll; }
void set_arima_ydiff(ArmaInfo * a)        { a.pflags |= PrivFlags.ydiff; }
void unset_arima_ydiff(ArmaInfo * a)      { a.pflags &= ~PrivFlags.ydiff; }
void set_arma_stdx(ArmaInfo * a)          { a.pflags |= PrivFlags.stdx; }

bool AR_included(const(ArmaInfo) * a, int i) { return a.pmask == null || a.pmask[i] == '1'; }
bool MA_included(const(ArmaInfo) * a, int i) { return a.qmask == null || a.qmask[i] == '1'; }
