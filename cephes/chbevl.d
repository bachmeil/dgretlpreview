/* Port of cephes/chbevl.c: the evaluation of a Chebyshev series, chbevl(x, array, n), the sum of n coefficients
   (highest order first) times the Chebyshev polynomials of the argument x, by Clenshaw's recurrence. The
   Bessel-function approximations (i0, i1, k0, k1) are built on it. The coefficient array must have n >= 2 entries:
   the do-while loop of the C decrements its counter before testing it, so for n = 1 or n = 0 it never stops (and
   reads past the array).

   Cephes Math Library Release 2.0, April 1987, Copyright 1985, 1987 by Stephen L. Moshier.

   The pointer walking of the C (*p++) is kept. Only D code (i0.d, i1.d, k0.d, k1.d) calls chbevl, so it has D
   linkage. */

module cephes.chbevl;

double chbevl(double x, const(double) * array, int n)
{
    double b0, b1, b2;
    const(double) * p = array;
    int i;

    b0 = *p++;
    b1 = 0.0;
    i = n - 1;

    do {
	b2 = b1;
	b1 = b0;
	b0 = x * b1 - b2 + *p++;
    } while (--i);

    return 0.5 * (b0-b2);
}
