/* Port of lib/src/gretl_restrict.h: the restriction-set type. In C it is opaque outside gretl_restrict.c; the
   struct is defined here (together with its row type) so that gretl_restrict.d can use it, while the other
   modules only handle pointers to it. The prototypes of the C header belong to gretl_restrict.d. */

module gretl.gretl_restricth;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;
import gretl.objstackh : GretlObjType;

struct RRow {
    int nterms;      /* number of terms in restriction */
    double * mult;   /* array of numerical multipliers on coeffs */
    int * eq;        /* array of equation numbers (for multi-equation case) */
    int * bnum;      /* array of coeff numbers */
    char * letter;   /* array of coeff letters */
    double rhs;      /* numerical value on right-hand side */
}

struct GretlRestriction {
    int g;                   /* number of restrictions (rows of R) */
    int gmax;                /* max. possible restrictions */
    int bmulti;              /* flag for case where coeffs need two indices */
    int amulti;              /* VECM only */
    int gb, ga;              /* restrictions on beta, alpha (VECM only) */
    int bcols, acols;        /* VECM only */
    int vecm;                /* pertains to VECM beta or alpha? */
    GretlMatrixRaw * R;         /* LHS restriction matrix */
    GretlMatrixRaw * q;         /* RHS restriction matrix */
    GretlMatrixRaw * Ra;        /* second LHS restriction matrix */
    GretlMatrixRaw * qa;        /* second RHS restriction matrix */
    char * mask;             /* selection mask for coeffs */
    RRow ** rows;            /* "atomic" restrictions */
    void * obj;              /* pointer to model, system */
    GretlObjType otype;      /* type of model, system */
    Gretlopt opt;            /* OPT_C is used for "coeffsum" */
    char * rfunc;            /* name of nonlinear restriction function */
    double test;             /* test statistic */
    double pval;             /* p-value of test statistic */
    double lnl;              /* loglikelihood */
    double bsum;             /* sum of coefficients */
    double bse;              /* standard error of bsum */
    int code;                /* identifier for statistic type */
}
