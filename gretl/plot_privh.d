/* Port of lib/src/plot_priv.h: the private header shared by graphing.c and plotbands.c (the prototypes belong to
   those files). BPMode loses its BP_ prefix, `gnuplot_info` is `GnuplotInfo`, and GPNA is a string enum. */

module gretl.plot_privh;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.graphingh : FitType, GptFlags;

enum BPMode {
    regular,
    blockmat
}

enum GPNA = "NaN";

struct GnuplotInfo {
    GptFlags flags;
    FitType fit;
    int * list;
    int t1;
    int t2;
    double xrange;
    char[64] xtics;
    char[16] xfmt;
    char[16] yfmt;
    const(char) * yformula;
    const(char) * facname;
    const(double) * x;
    GretlMatrixRaw * fvals; /* factor values */
    int n_fvals;         /* number of factor values */
    int fid;             /* factor variable ID */
    int * flist;         /* lisy of y-vars to be factorized */
    int * withlist;      /* for recording 'with' gnuplot options */
    int * y2vars;        /* record of variables for second y axis */
    int band;
    double ybase;
}
