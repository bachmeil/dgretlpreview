/* Port of cephes/stdtr.c: Student's t distribution. stdtr(k, t) is the integral of the t density with k > 0
   degrees of freedom from minus infinity to t (for integer k by the finite trigonometric series of Abramowitz and
   Stegun 26.7.3-4, otherwise through the incomplete beta integral), and stdtri(k, p) is its inverse (the t for
   which stdtr(k, t) = p).

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. incbet and incbi come from the ported modules; mtherr
   comes from cephes.mtherr. */

module cephes.stdtr;

import core.stdc.math : atan, fabs, sqrt;

import cephes.constants : MACHEP, MAXNUM, PI;
import cephes.incbet : incbet;
import cephes.incbi : incbi;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double stdtr(double rk, double t)
{
    double x, z, f, tz, p, xsqk;
    int k, j, is_int;

    if (rk <= 0) {
	mtherr("stdtr", CephesError.domain);
	return 0.0;
    } else {
	k = cast(int) rk;
	is_int = (rk - k) == 0.0;
    }

    if (t == 0) {
	return 0.5;
    }

    if (t < -2.0) {
	z = rk / (rk + t * t);
	p = 0.5 * incbet(0.5*rk, 0.5, z);
	return p;
    }

    if (is_int) {
	/* integer df */

	/*	compute integral from -t to +t */

	if (t < 0) {
	    x = -t;
	} else {
	    x = t;
	}

	z = 1.0 + (x * x)/rk;

	/* test if k is odd or even */
	if ((k & 1) != 0) {

	    /* computation for odd k */

	    xsqk = x/sqrt(rk);
	    p = atan(xsqk);
	    if (k > 1) {
		f = 1.0;
		tz = 1.0;
		j = 3;
		while ((j <= (k-2)) && ((tz/f) > MACHEP)) {
		    tz *= (j - 1)/(z * j);
		    f += tz;
		    j += 2;
		}
		p += f * xsqk/z;
	    }
	    p *= 2.0/PI;
	} else {

	    /* computation for even k */

	    f = 1.0;
	    tz = 1.0;
	    j = 2;

	    while ((j <= (k-2)) && ((tz/f) > MACHEP)) {
		tz *= (j - 1)/(z * j);
		f += tz;
		j += 2;
	    }
	    p = f * x/sqrt(z*rk);
	}

	/* common exit */
	
	if (t < 0) {
	    p = -p;	/* note destruction of relative accuracy */
	}

	p = 0.5 + 0.5 * p;

    } else {
	/* non-integer df */
	z = rk / (rk + t * t);
	p = 0.5 * incbet(0.5*rk, 0.5, z);
	if (t > 0) {
	    p = 1 - p;
	} 
    }


    return p;
}

double stdtri(double rk, double p)
{
    double t, z;
    int rflg;

    if (rk <= 0 || p <= 0.0 || p >= 1.0) {
	mtherr("stdtri", CephesError.domain);
	return 0.0;
    }

    if (p > 0.25 && p < 0.75) {
	if (p == 0.5) {
	    return 0.0;
	}
	z = 1.0 - 2.0 * p;
	z = incbi(0.5, 0.5*rk, fabs(z));
	t = sqrt(rk * z/(1.0-z));
	if (p < 0.5) {
	    t = -t;
	}
	return t;
    }

    rflg = -1;

    if (p >= 0.5) {
	p = 1.0 - p;
	rflg = 1;
    }

    z = incbi(0.5*rk, 0.5, 2.0*p);

    if (MAXNUM * z < rk) {
	return rflg * MAXNUM;
    }

    t = sqrt(rk/z - rk);

    return rflg * t;
}
