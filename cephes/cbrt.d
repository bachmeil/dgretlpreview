/* Port of cephes/cbrt.c: the cube root, cbrt(x), by a polynomial approximation of the cube root of the mantissa
   (peak relative error 9.2e-6), the exponent divided by three, and two Newton iterations. This is the cephes
   function, not libm's cbrt (the Bessel code in jv.d uses it, and its results depend on the difference in the
   last digit; see the jv entry in log.md).

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1991, 2000 by Stephen L. Moshier.

   Only jv.d calls it, so it has D linkage (a D function named cbrt does not clash with libm's). The four constants
   are __gshared variables, as in the C. */

module cephes.cbrt;

import core.stdc.math : frexp, isfinite, isnan, ldexp;

private __gshared double CBRT2 = 1.2599210498948731647672;
private __gshared double CBRT4 = 1.5874010519681994747517;
private __gshared double CBRT2I = 0.79370052598409973737585;
private __gshared double CBRT4I = 0.62996052494743658238361;

double cbrt(double x)
{
    int e, rem, sign;
    double z;

    if (isnan(x))
	return x;

    if (!isfinite(x))
	return x;

    if (x == 0)
	return x;

    if (x > 0) {
	sign = 1;
    } else {
	sign = -1;
	x = -x;
    }

    z = x;
    /* extract power of 2, leaving
     * mantissa between 0.5 and 1
     */
    x = frexp(x, &e);

    /* Approximate cube root of number between .5 and 1,
     * peak relative error = 9.2e-6
     */
    x = (((-1.3466110473359520655053e-1  * x
	   + 5.4664601366395524503440e-1) * x
	  - 9.5438224771509446525043e-1) * x
	 + 1.1399983354717293273738e0 ) * x
	+ 4.0238979564544752126924e-1;

    /* exponent divided by 3 */
    if (e >= 0) {
	rem = e;
	e /= 3;
	rem -= 3*e;
	if (rem == 1)
	    x *= CBRT2;
	else if (rem == 2)
	    x *= CBRT4;
    } else {
	/* argument less than 1 */
	e = -e;
	rem = e;
	e /= 3;
	rem -= 3*e;
	if (rem == 1)
	    x *= CBRT2I;
	else if (rem == 2)
	    x *= CBRT4I;
	e = -e;
    }

    /* multiply by power of 2 */
    x = ldexp(x, e);

    /* Newton iteration */
    x -= (x - (z/(x*x))) * 0.33333333333333333333;
    x -= (x - (z/(x*x))) * 0.33333333333333333333;

    return (sign < 0)? -x : x;
}
