/* Port of cephes/jn.c: the Bessel function of the first kind of integer order, J_n(x) (cephes_bessel_Jn), by
   backward recurrence from a continued fraction for J_(n+1)/J_n, normalised with J0 or J1. Orders 0, 1 and 2 are
   computed directly from cephes_j0 and cephes_j1; negative orders and arguments use the parity.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 2000 by Stephen L. Moshier.

   Nothing calls cephes_bessel_Jn (the real-order cephes_bessel_Jv handles integer orders itself), so it has D
   linkage. cephes_j0 and cephes_j1 come from cephes.j0 and cephes.j1. */

module cephes.jn;

import core.stdc.math : fabs;

import cephes.constants : MACHEP;
import cephes.j0 : cephes_j0;
import cephes.j1 : cephes_j1;

double cephes_bessel_Jn(int n, double x)
{
    double pkm2, pkm1, pk, xk, r, ans;
    int k, sign;

    if (n < 0) {
	n = -n;
	if ((n & 1) == 0) /* -1**n */
	    sign = 1;
	else
	    sign = -1;
    } else {
	sign = 1;
    }

    if (x < 0.0) {
	if (n & 1)
	    sign = -sign;
	x = -x;
    }

    if (n == 0)
	return sign * cephes_j0(x);
    if (n == 1)
	return sign * cephes_j1(x);
    if (n == 2)
	return sign * (2.0 * cephes_j1(x) / x - cephes_j0(x));

    if (x < MACHEP)
	return 0.0;

    /* continued fraction */
    k = 53;

    pk = 2 * (n + k);
    ans = pk;
    xk = x * x;

    do {
	pk -= 2.0;
	ans = pk - (xk/ans);
    } while (--k > 0);

    ans = x/ans;

    /* backward recurrence */

    pk = 1.0;
    pkm1 = 1.0/ans;
    k = n-1;
    r = 2 * k;

    do {
	pkm2 = (pkm1 * r - pk * x) / x;
	pk = pkm1;
	pkm1 = pkm2;
	r -= 2.0;
    } while (--k > 0);

    if (fabs(pk) > fabs(pkm1))
	ans = cephes_j1(x)/pk;
    else
	ans = cephes_j0(x)/pkm1;

    return sign * ans;
}
