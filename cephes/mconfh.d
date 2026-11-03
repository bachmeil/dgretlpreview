/* Port of cephes/mconf.h: the common include file of the cephes routines. What it contributes to the D
   modules is the package name and version, the complex-number struct used by polrt, and the build
   configuration that the C selects with preprocessor symbols (see the notes below).

   Not ported:
   - the prototypes of mtherr() and mtherr_with_arg(): they belong to the port of mtherr.c;
   - the includes of config.h, math.h, libprob.h (cephes.libprobh) and cephes.h: D modules import what they use;
   - the fallback macros isnan() and isfinite(): D has them in core.stdc.math (the C only defines them when
     math.h lacks them);
   - the long double complex type cmplxl (HAVE_LONG_DOUBLE): nothing in the cephes sources uses it;
   - the CMPLX marker define (also unused).

   Configuration the C selects here, which the ports of the other files follow: MINUSZERO is defined (const.c
   then sets NEGZERO to -0.0), and DEC, DENORMAL, UNK, IBMPC and MIEEE are not tested by any cephes source file
   other than as the unused DEC branch of incbet.c (DEC is not defined, so that branch is dropped). */

module cephes.mconfh;

/* Name of package */
enum PROB_PACKAGE = "cephes";

/* Version number of package */
enum VERSION = "2.7";

/* complex number (the "cmplx" of the C; used by polrt) */
struct Cmplx {
    double r;
    double i;
}
