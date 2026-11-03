/* Port of libgretl.h */

module gretl.libgretlh;


import std.math;
import gretl.ghashtable : GHashTable;
import gretl.gretl_bundleh : BundleType;
import gretl.glib;
import gretl.gretl_model : ModelDataItem, ModelTest;
import gretl.varinfo_privh : VarInfo;

enum GLIB_VERSION_MIN_REQUIRED = GLIB_VERSION_2_28;
enum GLIB_VERSION_MAX_ALLOWED = GLIB_VERSION_CUR_STABLE;

/* On MS Windows FILENAME_MAX is 260; on Linux it's typically
   4096; and on macOS we read that it's typically 1024. Here
   we define it to a value that should be big enough for
   most purposes but is not too wasteful of stack space. */
enum FILENAME_MAX = 1024;

enum YMD_READ_FMT = "%d-%d-%d";
enum YMD_WRITE_FMT = "%d-%02d-%02d";
enum YMD_WRITE_Y2_FMT = "%02d-%02d-%02d";
enum YMD_WRITE_Y4_FMT = "%04d-%02d-%02d";

/* The non-NLS versions of the translation macros */
string _(string s) { return s; }
string N_(string s) { return s; }
string gettext_noop(string s) { return s; }

enum MAXLINE = 131072; /* maximum length of command line */
enum MAXLABEL = 128; /* maximum length of descriptive labels for variables */
enum MAXLEN = 512; /* max length of regular "long" strings */
enum MAXDISP = 32; /* max length of "display names" for variables */
enum VNAMELEN = 32; /* space allocated for var names (including termination) */
enum OBSLEN = 16; /* space allocated for obs strings (including termination) */

enum M_PI = 3.1415926535897932384626432;
enum M_2PI = 6.2831853071795864769252864;

enum SQRT_2_PI = 2.506628274631000502415765284811;
enum LN_PI = 1.144729885849400174143427351353;
enum LN_2_PI = 1.837877066409345483560659472811;
enum LN_SQRT_2_PI = 0.918938533204672741780329736406;

enum PMAX_NOT_AVAILABLE = 666;

/* numbers smaller than the given limit will print as zero */
double screen_zero(double x) {
    return (fabs(x) > 1.0e-13) ? x : 0.0;
}

enum GretlType {
    none,
    bool_,
    int_,
    uint32,
    uint64,
    obs,
    list,
    double_,
    int_array,
    double_array,
    string_,
    cmplx_array,
    series,
    matrix,
    struct_,
    scalar_ref,
    series_ref,
    matrix_ref,
    string_ref,
    list_ref,
    useries,
    date,
    bundle,
    bundle_ref,
    array,
    array_ref,
    strings,
    matrices,
    bundles,
    lists,
    arrays,
    strings_ref,
    matrices_ref,
    bundles_ref,
    lists_ref,
    arrays_ref,
    void_,
    numeric,
    any
}

bool gretl_scalar_type(GretlType t) {
    return (t == GretlType.bool_ ||
            t == GretlType.int_ ||
            t == GretlType.uint32 ||
            t == GretlType.uint64 ||
            t == GretlType.obs ||
            t == GretlType.double_);
}

bool gretl_ref_type(GretlType t) {
    return (t == GretlType.scalar_ref ||
            t == GretlType.series_ref ||
            t == GretlType.matrix_ref ||
            t == GretlType.bundle_ref ||
            t == GretlType.string_ref ||
            t == GretlType.strings_ref ||
            t == GretlType.matrices_ref ||
            t == GretlType.bundles_ref ||
            t == GretlType.lists_ref ||
            t == GretlType.arrays_ref);
}

bool gretl_array_type(GretlType t) {
    return (t == GretlType.strings ||
            t == GretlType.matrices ||
            t == GretlType.bundles ||
            t == GretlType.lists ||
            t == GretlType.arrays);
}

enum TsCodes {
    cross_section,
    time_series,
    stacked_time_series,
    stacked_cross_section,
    panel_unknown,
    panel_side_by_side,
    special_time_series,
    structure_unknown
}

enum ProgressFlags {
    none,
    total,
    load_init,
    save_init,
    font_init,
    updater_init,
    finish
}

enum ProgressReturnFlags {
    ok,
    done,
    canceled
}

enum TestStats {
    none,
    normal_chisq,
    lm,
    f,
    lmf,
    harvey_collier,
    reset,
    lr,
    wald_chisq,
    sup_wald,
    z,
    student,
    lb_chisq,
    wf
}

enum Gretlopt {
    OPT_NONE = 0,
    OPT_A = 1 <<  0,
    OPT_B = 1 <<  1,
    OPT_C = 1 <<  2,
    OPT_D = 1 <<  3,
    OPT_E = 1 <<  4,
    OPT_F = 1 <<  5,
    OPT_G = 1 <<  6,
    OPT_H = 1 <<  7,
    OPT_I = 1 <<  8,
    OPT_J = 1 <<  9,
    OPT_K = 1 << 10,
    OPT_L = 1 << 11,
    OPT_M = 1 << 12,
    OPT_N = 1 << 13,
    OPT_O = 1 << 14,
    OPT_P = 1 << 15,
    OPT_Q = 1 << 16,
    OPT_R = 1 << 17,
    OPT_S = 1 << 18,
    OPT_T = 1 << 19,
    OPT_U = 1 << 20,
    OPT_V = 1 << 21,
    OPT_W = 1 << 22,
    OPT_X = 1 << 23,
    OPT_Y = 1 << 24,
    OPT_Z = 1 << 25,
    OPT_a = 1 << 26,
    OPT_b = 1 << 27,
    OPT_i = 1 << 28,
    OPT_UNSET = 1 << 30
}

enum GretlOp {
    eq  = '=',
    gt  = '>',
    lt  = '<',
    neq = 21,
    gte = 22,
    lte = 23
}

enum InfoCriterion {
    aic,
    bic,
    hqc,
    max
}

enum FcastMethod {
    static_,
    dynamic,
    auto_,
    kstep
}

enum FcastFlags {
    auto_ok      = 1 << 0,
    dynamic_ok   = 1 << 1,
    addobs_ok    = 1 << 2,
    integrate_ok = 1 << 3,
    expon_ok     = 1 << 4
}

struct Cmplx {
    double r;
    double i;

    /* Arithmetic, standing in for C99 `double _Complex`. Division uses
       Smith's algorithm; unlike C's __divdc3 it does not try to recover
       infinities/NaNs from NaN results. */
    Cmplx opBinary(string op)(Cmplx rhs) const
    {
        static if (op == "+") {
            return Cmplx(r + rhs.r, i + rhs.i);
        } else static if (op == "-") {
            return Cmplx(r - rhs.r, i - rhs.i);
        } else static if (op == "*") {
            return Cmplx(r * rhs.r - i * rhs.i, r * rhs.i + i * rhs.r);
        } else static if (op == "/") {
            if (fabs(rhs.r) >= fabs(rhs.i)) {
                double q = rhs.i / rhs.r;
                double d = rhs.r + rhs.i * q;
                return Cmplx((r + i * q) / d, (i - r * q) / d);
            } else {
                double q = rhs.r / rhs.i;
                double d = rhs.r * q + rhs.i;
                return Cmplx((r * q + i) / d, (i * q - r) / d);
            }
        } else {
            static assert(0, "Cmplx: operator " ~ op ~ " not supported");
        }
    }

    ref Cmplx opOpAssign(string op)(Cmplx rhs)
    {
        this = this.opBinary!op(rhs);
        return this;
    }
}

Cmplx conj(Cmplx z)
{
    return Cmplx(z.r, -z.i);
}

double cabs(Cmplx z)
{
    return hypot(z.r, z.i);
}

/* opaque types declared in libgretl.h */
/* struct gretl_bundle_ is private to gretl_bundle.c in the C (opaque in gretl_bundle.h); the layout is declared here so
   that gretl_bundle.d can use it, while the other modules only hold pointers */
struct GretlBundle {
    BundleType type;    /* see enum in gretl_bundleh.d */
    GHashTable * ht;    /* holds key/value pairs */
    char * creator;     /* name of function that built the bundle */
    void * data;        /* holds pointer to struct for some uses */
}
/* struct gretl_array_ is private to gretl_array.c in the C (opaque in gretl_array.h); the layout is declared here so
   that gretl_array.d can use it, while the other modules only hold pointers */
struct GretlArray {
    GretlType type;  /* type of data */
    int n;           /* number of elements */
    void ** data;    /* actual data array */
    double * mdata;  /* for matrix block */
}
public import gretl.gretl_string_table : GretlStringTable;

/* gretl data set */
struct DatasetRaw {
    int v;              /* number of data series, including const */
    int n;              /* number of observations */
    int pd;             /* periodicity or frequency of data */
    int structure;      /* time series, cross section or whatever */
    double sd0;         /* floating-point representation of stobs */
    int t1, t2;         /* start and end of current sample */
    char[OBSLEN] stobs;  /* string representation of starting obs */
    char[OBSLEN] endobs; /* string representation of ending obs */
    double ** Z;        /* two-dimensional data array */
    char ** varname;    /* array of names of series */
    VarInfo ** varinfo; /* array of metadata per series */
    char markers;       /* NO_MARKERS (0), REGULAR MARKERS or DAILY_DATE_STRINGS */
    char modflag;       /* binary flag for dataset modified or not */
    char ** S;          /* to hold observation markers */
    char * descrip;     /* to hold info on data sources, etc. */
    char * submask;     /* subsampling mask */
    char * restriction; /* record of sub-sampling restriction */
    char * padmask;     /* record of padding to re-balance panel data */
    char * mapfile;     /* name of associated map (polygons) file, if any */
    ulong rseed;        /* resampling seed */
    int auxiliary;      /* 0 for regular dataset, 1 for auxiliary dataset */
    int n_varinfo;      /* # of named series with metadata, if < @v above */
    char * pangrps;     /* panel-only: name of series holding group names */
    int panel_pd;       /* panel-only: panel time-series frequency */
    double panel_sd0;   /* panel-only: time-series start */
}

struct VMatrix {
    int ci;
    int dim;
    int t1, t2;
    int nmin, nmax;
    int ncrit; /* for assessing critical values */
    char ** names;
    double * vec;
    double * xbar;
    double * ssx;
    int * list;
    int missing;
}

struct Sample {
    int t1;
    int t2;
    ulong rseed;
}

struct ArInfo {
    int * arlist;          /* list of autoregressive lags */
    double * rho;          /* array of autoreg. coeffs. */
    double * sderr;        /* and their standard errors */
}

/* struct to hold model results */
struct Model {
    int ID;                      /* ID number for model */
    int refcount;                /* for saving/deleting */
    int ci;                      /* "command index" -- estimation method */
    Gretlopt opt;                /* record of options */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    char * submask;              /* keep track of sub-sample in force
                                    when model was estimated */
    char * missmask;             /* missing observations mask */
    Sample smpl;                 /* numeric start and end of current sample
                                    when model was estimated */
    int full_n;                  /* full length of dataset on estimation */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int * list;                  /* list of variables by ID number */
    int * xlist;                 /* list of regressors */
    int ifc;                     /* = 1 if the equation includes a constant,
                                    else = 0 */
    int nwt;                     /* ID number of the weight variable (WLS) */
    int aux;                     /* code representing the sort of
                                    auxiliary regression this is (or not) */
    double * coeff;              /* array of coefficient estimates */
    double * sderr;              /* array of estimated std. errors */
    double * uhat;               /* regression residuals */
    double * yhat;               /* fitted values from regression */
    double * xpx;                /* X'X matrix, in packed form */
    double * vcv;                /* VCV matrix for coefficient estimates */
    double ess, tss;             /* Error and Total Sums of Squares */
    double sigma;                /* Standard error of regression */
    double rsq, adjrsq;          /* Unadjusted and adjusted R^2 */
    double fstt;                 /* overall F-statistic */
    double chisq;                /* overall chi-square statistic */
    double lnL;                  /* log-likelihood */
    double ybar, sdy;            /* mean and std. dev. of dependent var. */
    double[InfoCriterion.max] criterion; /* array of model selection statistics */
    double dw, rho;              /* Durbin-Watson stat. and estimated 1st
                                    order autocorrelation coefficient */
    ArInfo * arinfo;             /* pointer to struct to hold special info for
                                    autoregressive model */
    int errcode;                 /* Error code in case of failure */
    char * name;                 /* for use in GUI */
    char * depvar;               /* name of dependent var in special cases */
    int nparams;                 /* number of named model parameters */
    char ** params;              /* for named model parameters */
    gint64 esttime;              /* time of estimation */
    int ntests;                  /* number of attached test results */
    ModelTest * tests;           /* attached hypothesis test results */
    DatasetRaw * dataset;           /* for handling models estimated on a
                                    sub-sampled portion of the dataset */
    int n_data_items;            /* number of extra data items */
    ModelDataItem ** data_items; /* pointer to additional data */
}

/* pointer flavour of the (identity) translation function, for the arrays of N_() strings of the C */
const(char) * _(const(char) * s) { return s; }
