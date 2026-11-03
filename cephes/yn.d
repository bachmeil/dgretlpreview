/* Port of cephes/yn.c: the Bessel function of the second kind of integer order, Y_n(x) (cephes_bessel_Yn), by
   forward recurrence in n from Y0 and Y1 (cephes_y0, cephes_y1); negative orders use the parity. Returns -MAXNUM,
   with an error message, for x <= 0 and n >= 2.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 2000 by Stephen L. Moshier.

   Only struve.d calls cephes_bessel_Yn, so it has D linkage. cephes_y0 and cephes_y1 come from cephes.j0 and
   cephes.j1; mtherr from cephes.mtherr. */

module cephes.yn;

import cephes.constants : MAXNUM;
import cephes.j0 : cephes_y0;
import cephes.j1 : cephes_y1;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double cephes_bessel_Yn(int n, double x)
{
    double an, anm1, anm2, r;
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

    if (n == 0)
	return sign * cephes_y0(x);
    if (n == 1)
	return sign * cephes_y1(x);

    /* test for overflow */
    if (x <= 0.0) {
	mtherr("yn", CephesError.sing);
	return -MAXNUM;
    }

    /* forward recurrence on n */

    anm2 = cephes_y0(x);
    anm1 = cephes_y1(x);
    k = 1;
    r = 2 * k;
    do {
	an = r * anm1 / x  -  anm2;
	anm2 = anm1;
	anm1 = an;
	r += 2.0;
	++k;
    } while (k < n);

    return sign * an;
}
