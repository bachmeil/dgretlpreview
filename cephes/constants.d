/* Port of cephes/const.c: the global mathematical constants and machine-arithmetic parameters of the cephes
   routines (MACHEP, MAXNUM, MAXLOG, MINLOG, PI, PIO2, PIO4, SQRT2, SQRTH, LOG2E, SQ2OPI, LOGE2, LOGSQ2, THPIO4,
   TWOOPI, NEGZERO, UFLOWTHRESH). The module is named constants (const is a D keyword, so const.c cannot keep its name).

   Differences from the C:
   - the values are the "Beebe" set (the C defines BEEBE: values found by running the Beebe floating-point test
     suite), and the unused alternative set is dropped;
   - INFINITY and NAN are not defined here: the C only defines them when math.h lacks them, and D has them in
     core.stdc.math;
   - NEGZERO is -0.0 (the C defines MINUSZERO, see mconfh.d);
   - the constants keep the C names. They stay ordinary (__gshared) variables, as in the C, and not immutable or
     enum constants: D may fold a constant with a literal initialiser at higher than double precision, so that,
     for instance, an immutable MACHEP does not compare equal to DBL_EPSILON, whereas the C (and a variable)
     holds the value rounded to double. Nothing in the cephes sources assigns to them. */

module cephes.constants;

/* values gathered by running the Beebe floating-point test suite */
__gshared double MACHEP      = 2.22044604925031308e-16;
__gshared double UFLOWTHRESH = 2.22507385850720188e-308;
__gshared double MAXNUM      = 1.79769313486231571e+308;

__gshared double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
__gshared double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */

__gshared double PI     =  3.14159265358979323846;       /* pi */
__gshared double PIO2   =  1.57079632679489661923;       /* pi/2 */
__gshared double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
__gshared double SQRT2  =  1.41421356237309504880;       /* sqrt(2) */
__gshared double SQRTH  =  7.07106781186547524401E-1;    /* sqrt(2)/2 */
__gshared double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
__gshared double SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
__gshared double LOGE2  =  6.93147180559945309417E-1;    /* log(2) */
__gshared double LOGSQ2 =  3.46573590279972654709E-1;    /* log(2)/2 */
__gshared double THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
__gshared double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */

__gshared double NEGZERO = -0.0;
