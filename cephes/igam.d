/* Port of cephes/igam.c: the incomplete gamma integral. igam(a, x) is the left tail, the integral of
   t^(a-1) exp(-t) from 0 to x divided by Gamma(a) (by a power series for x <= 1 or x <= a, otherwise as 1 -
   igamc), and igamc(a, x) the right tail, the complementary integral (by a continued fraction for x >= 1 and
   x >= a, otherwise as 1 - igam).

   Cephes Math Library Release 2.8, June 2000, Copyright 1985, 1987, 2000 by Stephen L. Moshier.

   igam and igamc have D linkage. The `#if 0` underflow message in igamc is dropped. big and biginv are __gshared
   variables, as in the C. lgam comes from cephes.gamma; mtherr comes from cephes.mtherr. */

module cephes.igam;

import core.stdc.math : exp, fabs, log;

import cephes.constants : MACHEP, MAXLOG;
import cephes.gamma : lgam;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

private __gshared double big = 4.503599627370496e15;
private __gshared double biginv =  2.22044604925031308085e-16;

double igamc(double a, double x)
{
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    if (x <= 0 || a <= 0) {
	return 1.0;
    }

    if (x < 1.0 || x < a) {
	return 1.0 - igam(a,x);
    }

    ax = a * log(x) - x - lgam(a);

    if (ax < -MAXLOG) {
	return 0.0;
    }

    ax = exp(ax);

    /* continued fraction */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1/qkm1;

    do {
	c += 1.0;
	y += 1.0;
	z += 2.0;
	yc = y * c;
	pk = pkm1 * z  -  pkm2 * yc;
	qk = qkm1 * z  -  qkm2 * yc;
	if (qk != 0) {
	    r = pk/qk;
	    t = fabs((ans - r)/r);
	    ans = r;
	} else {
	    t = 1.0;
	}
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if (fabs(pk) > big) {
	    pkm2 *= biginv;
	    pkm1 *= biginv;
	    qkm2 *= biginv;
	    qkm1 *= biginv;
	}
    } while (t > MACHEP);

    return ans * ax;
}

/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */

double igam(double a, double x)
{
    double ans, ax, c, r;

    if (x <= 0 || a <= 0) {
	return 0.0;
    }

    if (x > 1.0 && x > a) {
	return 1.0 - igamc(a,x);
    }

    /* Compute x**a * exp(-x) / gamma(a) */
    ax = a * log(x) - x - lgam(a);
    if (ax < -MAXLOG) {
	mtherr("igam", CephesError.underflow);
	return 0.0;
    }

    ax = exp(ax);

    /* power series */
    r = a;
    c = 1.0;
    ans = 1.0;

    do {
	r += 1.0;
	c *= x/r;
	ans += c;
    } while (c/ans > MACHEP);

    return ans * ax/a;
}
