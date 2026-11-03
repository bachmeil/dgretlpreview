/* Port of minpack/minpack.h: the callback typedefs and the min / max macros (as functions). The prototypes of the routines are not
   here: the routines are defined, extern (C), in gretl.minpack.enorm, fdjac2, lmder, lmdif, lmpar, qrfac and qrsolv (setulb_ in lbfgsb). */

module gretl.minpackh;


alias S_fp6 = extern (C) int function(int, int, double *, double *, int *, void *);
alias S_fp8 = extern (C) int function(int, int, double *, double *, double *, int, int *, void *);

/* the min and max macros of minpack.h */
T min(T)(T a, T b) { return (a <= b) ? a : b; }
T max(T)(T a, T b) { return (a >= b) ? a : b; }
