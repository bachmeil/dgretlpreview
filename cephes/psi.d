/* Port of cephes/psi.c: the psi (digamma) function, psi(x) = d/dx log Gamma(x), by the asymptotic series after
   recurrence up to x = 10, with the exact finite sum for positive integers up to 10 and the reflection formula for
   x <= 0 (singular at the non-positive integers, where MAXNUM is returned with an error message).

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier.

   Only D code (hyp2f1.d and pvalues.d) calls psi and no C file does, so it has D linkage. The coefficient array
   is private __gshared. polevl comes from cephes.polevl. */

module cephes.psi;

import core.stdc.math : floor, log, tan;

import cephes.constants : MAXNUM, PI;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr_with_arg;
import cephes.polevl : polevl;

private __gshared double[7] A = [
 8.33333333333333333333E-2,
-2.10927960927960927961E-2,
 7.57575757575757575758E-3,
-4.16666666666666666667E-3,
 3.96825396825396825397E-3,
-8.33333333333333333333E-3,
 8.33333333333333333333E-2];

private enum EUL = 0.57721566490153286061;

double psi(double x)
{
    double p, q, nz, s, w, y, z;
    int i, n, negative = 0;

    nz = 0.0;

    if (x <= 0.0) {
	negative = 1;
	q = x;
	p = floor(q);
	if (p == q) {
	    mtherr_with_arg("psi", CephesError.sing, x);
	    return MAXNUM;
	}
	/* Remove the zeros of tan(PI x) by subtracting the nearest
	   integer from x
	 */
	nz = q - p;
	if (nz != 0.5) {
	    if (nz > 0.5) {
		p += 1.0;
		nz = q - p;
	    }
	    nz = PI/tan(PI*nz);
	} else {
	    nz = 0.0;
	}
	x = 1.0 - x;
    }

    /* check for positive integer up to 10 */
    if (x <= 10.0 && x == floor(x)) {
	y = 0.0;
	n = cast(int) x;
	for (i=1; i<n; i++) {
	    w = i;
	    y += 1.0/w;
	}
	y -= EUL;
	goto done;
    }

    s = x;
    w = 0.0;
    while (s < 10.0) {
	w += 1.0/s;
	s += 1.0;
    }

    if (s < 1.0e17) {
	z = 1.0/(s * s);
	y = z * polevl(z, A.ptr, 6);
    } else {
	y = 0.0;
    }

    y = log(s) - (0.5/s) - y - w;

 done:

    if (negative) {
	y -= nz;
    }

    return y;
}
