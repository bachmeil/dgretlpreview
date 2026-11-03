/* Port of cephes/ndtri.c: the inverse of the normal distribution function, ndtri(y): the x for which the area
   under the Gaussian density from minus infinity to x is y. Rational approximations for y between exp(-2) and
   1 - exp(-2) (in y - 1/2), and for smaller tails in z = sqrt(-2 log y) split at z = 8; the upper tail uses the
   symmetry.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier.

   ndtri has D linkage. The coefficient arrays and s2pi are __gshared variables, as in the C. mtherr_with_arg
   comes from cephes.mtherr; polevl and p1evl come from cephes.polevl. */

module cephes.ndtri;

import core.stdc.math : log, sqrt;

import cephes.constants : MAXNUM;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr_with_arg;
import cephes.polevl : polevl, p1evl;

/* sqrt(2pi) */
private __gshared double s2pi = 2.50662827463100050242E0;

/* approximation for 0 <= |y - 0.5| <= 3/8 */
private __gshared double[5] P0 = [
    -5.99633501014107895267E1,
    9.80010754185999661536E1,
    -5.66762857469070293439E1,
    1.39312609387279679503E1,
    -1.23916583867381258016E0];

private __gshared double[8] Q0 = [
    /* 1.00000000000000000000E0,*/
    1.95448858338141759834E0,
    4.67627912898881538453E0,
    8.63602421390890590575E1,
    -2.25462687854119370527E2,
    2.00260212380060660359E2,
    -8.20372256168333339912E1,
    1.59056225126211695515E1,
    -1.18331621121330003142E0];

/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
 * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
 */

private __gshared double[9] P1 = [
    4.05544892305962419923E0,
    3.15251094599893866154E1,
    5.71628192246421288162E1,
    4.40805073893200834700E1,
    1.46849561928858024014E1,
    2.18663306850790267539E0,
    -1.40256079171354495875E-1,
    -3.50424626827848203418E-2,
    -8.57456785154685413611E-4];

private __gshared double[8] Q1 = [
    /*  1.00000000000000000000E0,*/
    1.57799883256466749731E1,
    4.53907635128879210584E1,
    4.13172038254672030440E1,
    1.50425385692907503408E1,
    2.50464946208309415979E0,
    -1.42182922854787788574E-1,
    -3.80806407691578277194E-2,
    -9.33259480895457427372E-4];

/* Approximation for interval z = sqrt(-2 log y) between 8 and 64
 * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
 */

private __gshared double[9] P2 = [
    3.23774891776946035970E0,
    6.91522889068984211695E0,
    3.93881025292474443415E0,
    1.33303460815807542389E0,
    2.01485389549179081538E-1,
    1.23716634817820021358E-2,
    3.01581553508235416007E-4,
    2.65806974686737550832E-6,
    6.23974539184983293730E-9];

private __gshared double[8] Q2 = [
    /*  1.00000000000000000000E0,*/
    6.02427039364742014255E0,
    3.67983563856160859403E0,
    1.37702099489081330271E0,
    2.16236993594496635890E-1,
    1.34204006088543189037E-2,
    3.28014464682127739104E-4,
    2.89247864745380683936E-6,
    6.79019408009981274425E-9];

double ndtri(double y0)
{
    double expm2 = 0.13533528323661269189; /* exp(-2) */
    double x, y, z, y2, x0, x1;
    int code;

    if (y0 <= 0.0) {
	mtherr_with_arg("ndtri", CephesError.domain, y0);
	return -MAXNUM;
    }

    if (y0 >= 1.0) {
	mtherr_with_arg("ndtri", CephesError.domain, y0);
	return MAXNUM;
    }

    code = 1;
    y = y0;

    if (y > (1.0 - expm2)) {
	y = 1.0 - y;
	code = 0;
    }

    if (y > expm2) {
	y = y - 0.5;
	y2 = y * y;
	x = y + y * (y2 * polevl(y2, P0.ptr, 4)/p1evl(y2, Q0.ptr, 8));
	x = x * s2pi; 
	return x;
    }

    x = sqrt( -2.0 * log(y) );
    x0 = x - log(x)/x;

    z = 1.0/x;

    if (x < 8.0) { 
	/* y > exp(-32) = 1.2664165549e-14 */
	x1 = z * polevl(z, P1.ptr, 8) / p1evl(z, Q1.ptr, 8);
    } else {
	x1 = z * polevl(z, P2.ptr, 8) / p1evl(z, Q2.ptr, 8);
    }

    x = x0 - x1;

    if (code != 0) {
	x = -x;
    }

    return x;
}
