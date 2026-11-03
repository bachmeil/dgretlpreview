/* Port of describe.h */

module gretl.describeh;


import gretl.libgretlh;
import gretl.pvaluesh;

struct Summary {
    Gretlopt opt;
    int n;
    int weight_var;
    int * misscount;
    int * list;
    int * ival;
    double * stats;
    double * mean;
    double * median;
    double * sd;
    double * skew;
    double * xkurt;
    double * low;
    double * high;
    double * cv;
    double * perc05;
    double * perc95;
    double * iqr;
    double sw;
    double sb;
}

struct FreqDist {
    char[VNAMELEN] varname;  /* for ID purposes */
    char[MAXDISP] gname;     /* name to use in plots */
    int discrete;            /* 1 if series contains integers */
    int strvals;             /* 1 if series is string-valued */
    DistCode dist;           /* code for theoretical distribution */
    int numbins;             /* number of bins or intervals */
    double xbar, sdx;        /* mean and std dev of variable */
    double * midpt;          /* array of midpoints of intervals */
    double * endpt;          /* array of endpoints of intervals */
    char ** S;               /* array of string values */
    int * f;                 /* frequencies */
    double test;             /* Chi-squared statistic for testing for a
                                Gaussian distribution, or z statistic
                                for testing for Gamma dist. */
    int n;                   /* included observation */
    int t1, t2;              /* sample limits */
}

struct Xtab {
    char * rname;          /* identifier of rows series */
    char * cname;          /* identifier of cols series */
    char ** Sr;            /* array of row string values */
    char ** Sc;            /* array of column string values */
    int rows, cols;        /* dimensions of table */
    double * rval, cval;   /* row and column numeric values */
    int * rtotal, ctotal;  /* marginal totals */
    int ** f;              /* array of frequencies */
    int n, missing;        /* observation counts */
    int t1, t2;            /* sample limits */
    int rstrs;             /* row string-valued flag */
    int cstrs;             /* col string-valued flag */
}
