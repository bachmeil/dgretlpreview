/* Port of cephes/bdtr.c: the binomial distribution. bdtr(k, n, p) is the sum of the first k + 1 terms of the
   binomial distribution with n trials and success probability p (the area from 0 to k), bdtrc(k, n, p) the sum of
   the terms k + 1 to n, and bdtri(k, n, y) the p for which bdtr(k, n, p) = y. Computed through the incomplete beta
   integral and its inverse.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. incbet and incbi come from the ported modules. mtherr
   comes from cephes.mtherr. cephes_exp and cephes_log come from cephes.unity. */

module cephes.bdtr;

import core.stdc.math : pow;

import cephes.incbet : incbet;
import cephes.incbi : incbi;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;
import cephes.unity : cephes_exp, cephes_log;

double bdtrc(int k, int n, double p)
{
    double dk, dn;

    if (p < 0.0 || p > 1.0)
	goto domerr;
    if (k < 0)
	return 1.0;

    if (n < k) {
    domerr:
	mtherr("bdtrc", CephesError.domain);
	return 0.0;
    }

    if (k == n)
	return 0.0;
    dn = n - k;
    if (k == 0) {
	if (p < .01)
	    dk = -cephes_exp(dn * cephes_log(-p));
	else
	    dk = 1.0 - pow(1.0 - p, dn);
    } else {
	dk = k + 1;
	dk = incbet(dk, dn, p);
    }

    return dk;
}

double bdtr(int k, int n, double p)
{
    double dk, dn;

    if (p < 0.0 || p > 1.0)
	goto domerr;
    if (k < 0 || n < k) {
    domerr:
	mtherr("bdtr", CephesError.domain);
	return 0.0;
    }

    if (k == n)
	return 1.0;

    dn = n - k;
    if (k == 0) {
	dk = pow(1.0 - p, dn);
    } else {
	dk = k + 1;
	dk = incbet(dn, dk, 1.0 - p);
    }

    return dk;
}

double bdtri(int k, int n, double y)
{
    double dk, dn, p;

    if (y < 0.0 || y > 1.0)
	goto domerr;
    if (k < 0 || n <= k) {
    domerr:
	mtherr("bdtri", CephesError.domain);
	return 0.0;
    }

    dn = n - k;
    if (k == 0) {
	if (y > 0.8)
	    p = -cephes_exp(cephes_log(y - 1.0) / dn);
	else
	    p = 1.0 - pow(y, 1.0 / dn);
    } else {
	dk = k + 1;
	p = incbet(dn, dk, 0.5);
	if (p > 0.5)
	    p = incbi(dk, dn, 1.0 - y);
	else
	    p = 1.0 - incbi(dn, dk, y);
    }

    return p;
}
