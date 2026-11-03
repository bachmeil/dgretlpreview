/* Port of cephes/igami.c: the inverse of the complemented incomplete gamma integral, igami(a, y0): the x for which
   igamc(a, x) = y0. It starts from an approximation based on the inverse normal distribution (Wilson-Hilferty),
   refines by Newton's method, and falls back to interval halving if Newton does not converge; the solution is
   returned as MAXNUM when the initial approximation is not finite.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   The gotos of the C (ihalve, done) are kept. No C file calls igami (this was the last C file of the cephes
   library), so it has D linkage. igamc, lgam, ndtri and mtherr come from the ported modules. */

module cephes.igami;

import core.stdc.math : exp, fabs, isfinite, log, sqrt;

import cephes.constants : MACHEP, MAXLOG, MAXNUM;
import cephes.gamma : lgam;
import cephes.igam : igamc;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;
import cephes.ndtri : ndtri;

double igami(double a, double y0)
{
    double x0, x1, x, yl, yh, y, d, lgm, dithresh;
    int i, dir;

    /* bound the solution */
    x0 = MAXNUM;
    yl = 0;
    x1 = 0;
    yh = 1.0;
    dithresh = 5.0 * MACHEP;

    /* approximation to inverse function */
    d = 1.0/(9.0*a);
    y = (1.0 - d - ndtri(y0) * sqrt(d));
    x = a * y * y * y;

    if (!isfinite(x)) {
	/* this is not going to work */
	return MAXNUM;
    }

    lgm = lgam(a);

    for (i=0; i<10; i++) {
	if (x > x0 || x < x1)
	    goto ihalve;
	y = igamc(a,x);
	if (y < yl || y > yh)
	    goto ihalve;
	if (y < y0) {
	    x0 = x;
	    yl = y;
	} else {
	    x1 = x;
	    yh = y;
	}
	/* compute the derivative of the function at this point */
	d = (a - 1.0) * log(x) - x - lgm;
	if (d < -MAXLOG)
	    goto ihalve;
	d = -exp(d);
	/* compute the step to the next approximation of x */
	d = (y - y0)/d;
	if (fabs(d/x) < MACHEP)
	    goto done;
	x = x - d;
    }

    /* Resort to interval halving if Newton iteration did not converge. */
 ihalve:

    d = 0.0625;
    if (x0 == MAXNUM) {
	if (x <= 0.0)
	    x = 1.0;
	while (x0 == MAXNUM && isfinite(d) && isfinite(x)) {
	    x = (1.0 + d) * x;
	    y = igamc(a, x);
	    if (y < y0) {
		x0 = x;
		yl = y;
		break;
	    }
	    d = d + d;
	}
    }
    d = 0.5;
    dir = 0;

    for (i=0; i<400; i++) {
	x = x1  +  d * (x0 - x1);
	y = igamc(a, x);
	lgm = (x0 - x1)/(x1 + x0);
	if (fabs(lgm) < dithresh)
	    break;
	lgm = (y - y0)/y0;
	if (fabs(lgm) < dithresh)
	    break;
	if (x <= 0.0)
	    break;
	if (y >= y0) {
	    x1 = x;
	    yh = y;
	    if (dir < 0) {
		dir = 0;
		d = 0.5;
	    }
	    else if (dir > 1)
		d = 0.5 * d + 0.5; 
	    else
		d = (y0 - yl)/(yh - yl);
	    dir += 1;
	} else {
	    x0 = x;
	    yl = y;
	    if (dir > 0) {
		dir = 0;
		d = 0.5;
	    } else if (dir < -1) {
		d = 0.5 * d;
	    } else {
		d = (y0 - yl)/(yh - yl);
	    }
	    dir -= 1;
	}
    }

    if (x == 0.0)
	mtherr("igami", CephesError.underflow);

 done:
    return x;
}
