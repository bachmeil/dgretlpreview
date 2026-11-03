/* Port of cephes/k1.c: the modified Bessel function of the third kind of order one, K1(x) (cephes_bessel_K1), and
   its exponentially scaled form exp(x) K1(x) (k1e), by Chebyshev expansions on 0 < x <= 2 (for x (K1(x) - log(x/2)
   I1(x))) and on 2 < x < infinity (for exp(x) sqrt(x) K1(x)). Both return MAXNUM, with an error message, for
   x <= 0.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 2000 by Stephen L. Moshier.

   Only genfuncs.d calls cephes_bessel_K1 and k1e has no caller, so both have D linkage. The coefficient arrays
   are private __gshared. cephes_bessel_I1 comes from cephes.i1; mtherr comes from cephes.mtherr; chbevl comes
   from cephes.chbevl. */

module cephes.k1;

import core.stdc.math : exp, log, sqrt;

import cephes.chbevl : chbevl;
import cephes.constants : MAXNUM;
import cephes.i1 : cephes_bessel_I1;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

private __gshared double[11] A = [
-7.02386347938628759343E-18,
-2.42744985051936593393E-15,
-6.66690169419932900609E-13,
-1.41148839263352776110E-10,
-2.21338763073472585583E-8,
-2.43340614156596823496E-6,
-1.73028895751305206302E-4,
-6.97572385963986435018E-3,
-1.22611180822657148235E-1,
-3.53155960776544875667E-1,
 1.52530022733894777053E0];

/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */

private __gshared double[25] B = [
-5.75674448366501715755E-18,
 1.79405087314755922667E-17,
-5.68946255844285935196E-17,
 1.83809354436663880070E-16,
-6.05704724837331885336E-16,
 2.03870316562433424052E-15,
-7.01983709041831346144E-15,
 2.47715442448130437068E-14,
-8.97670518232499435011E-14,
 3.34841966607842919884E-13,
-1.28917396095102890680E-12,
 5.13963967348173025100E-12,
-2.12996783842756842877E-11,
 9.21831518760500529508E-11,
-4.19035475934189648750E-10,
 2.01504975519703286596E-9,
-1.03457624656780970260E-8,
 5.74108412545004946722E-8,
-3.50196060308781257119E-7,
 2.40648494783721712015E-6,
-1.93619797416608296024E-5,
 1.95215518471351631108E-4,
-2.85781685962277938680E-3,
 1.03923736576817238437E-1,
 2.72062619048444266945E0];

double cephes_bessel_K1(double x)
{
    double y, z;

    z = 0.5 * x;
    if (z <= 0.0) {
	mtherr("k1", CephesError.domain);
	return MAXNUM;
    }

    if (x <= 2.0) {
	y = x * x - 2.0;
	y =  log(z) * cephes_bessel_I1(x) + chbevl(y, A.ptr, 11) / x;
	return y;
    }

    return exp(-x) * chbevl(8.0/x - 2.0, B.ptr, 25) / sqrt(x);
}

double k1e(double x)
{
    double y;

    if (x <= 0.0) {
	mtherr("k1e", CephesError.domain);
	return MAXNUM;
    }

    if (x <= 2.0) {
	y = x * x - 2.0;
	y =  log(0.5 * x) * cephes_bessel_I1(x) + chbevl(y, A.ptr, 11) / x;
	return y * exp(x);
    }

    return chbevl(8.0/x - 2.0, B.ptr, 25) / sqrt(x);
}
