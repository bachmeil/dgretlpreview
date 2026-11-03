/* Port of the struct in lib/src/gretl_tdisagg.h (the function
   prototypes there get no D module) */

module gretl.gretl_tdisaggh;


import gretl.libgretlh : GretlType;
import gretl.gretl_matrixh : GretlMatrixRaw;

struct TdisaggInfo {
    int ynum, xnum;            /* series IDs, if applicable */
    const(int) * ylist;        /* series lists, if applicable */
    const(int) * xlist;
    const(double) * yval;      /* series values */
    const(double) * xval;
    GretlMatrixRaw * Y;           /* matrix values */
    GretlMatrixRaw * X;
    int xmidas;                /* 0/1: X given as a MIDAS list? */
    int efac;                  /* expansion factor */
    int extmax;                /* max # of extrapolated obs */
    GretlType targ;            /* target return type */
}
