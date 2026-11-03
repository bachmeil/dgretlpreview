/* Port of cephes/iv.c: the modified Bessel function of the first kind of real order, I_v(x) (cephes_bessel_Iv),
   in terms of the confluent hypergeometric function: I_v(x) = (x/2)^v exp(-x) 1F1(v + 1/2; 2v + 1; 2x) /
   Gamma(v + 1). For negative x the order must be an integer.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier.

   Only genfuncs.d calls cephes_bessel_Iv, so it has D linkage. hyperg, cephes_gamma and mtherr come from the
   ported modules. */

module cephes.iv;

import core.stdc.math : exp, fabs, floor, log;

import cephes.constants : MAXNUM;
import cephes.gamma : cephes_gamma;
import cephes.hyperg : hyperg;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double cephes_bessel_Iv(double v, double x)
{
    double ax, t = floor(v);
    int sign;

    /* If v is a negative integer, invoke symmetry */
    if (v < 0.0 && t == v) {
	t = v = -v;
    }

    /* If x is negative, require v to be an integer */
    sign = 1;
    if (x < 0.0) {
	if (t != v) {
	    mtherr("iv", CephesError.domain);
	    return 0.0;
	}
	if (v != 2.0 * floor(v/2.0))
	    sign = -1;
    }

    /* Avoid logarithm singularity */
    if (x == 0.0) {
	if (v == 0.0)
	    return 1.0;
	if (v < 0.0) {
	    mtherr("iv", CephesError.overflow);
	    return MAXNUM;
	} else
	    return 0.0;
    }

    ax = fabs(x);
    t = v * log(0.5 * ax) - x;
    t = sign * exp(t) / cephes_gamma(v + 1.0);
    ax = v + 0.5;

    return t * hyperg(ax, 2*ax, 2*x);
}

