/* Port of cephes/expx2.c: exp(x^2) (sign = 1) or exp(-x^2) (sign = -1) computed without the error amplification
   of squaring first: x is split into an exact multiple m of 1/128 and a small residual f, x^2 = m^2 + (2mf + f^2)
   with m^2 exact, and the result is exp(m^2) exp(2mf + f^2). Returns infinity when x^2 exceeds the largest
   argument of exp.

   Cephes Math Library Release 2.9, June 2000, Copyright 2000 by Stephen L. Moshier.

   Only ndtr.d calls expx2 (cephes_erfc does, see bugs.md #104), so it has D linkage. M and MINV are macros in the C
   (exact powers of two) and enums here. */

module cephes.expx2;

import core.stdc.math : exp, fabs, floor;

import cephes.constants : MAXLOG;

private enum M = 128.0;
private enum MINV = .0078125;

double expx2(double x, int sign)
{
    double u, u1, m, f;

    x = fabs(x);

    if (sign < 0) {
	x = -x;
    }

    /* Represent x as an exact multiple of M plus a residual.
       M is a power of 2 chosen so that exp(m * m) does not overflow
       or underflow and so that |x - m| is small.  */
    m = MINV * floor(M * x + 0.5);
    f = x - m;

    /* x^2 = m^2 + 2mf + f^2 */
    u = m * m;
    u1 = 2 * m * f + f * f;

    if (sign < 0) {
	u = -u;
	u1 = -u1;
    }

    if (u + u1 > MAXLOG) {
	return double.infinity;
    }

    /* u is exact, u1 is small.  */
    u = exp(u) * exp(u1);

    return u;
}
