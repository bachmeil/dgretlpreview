/* Port of cephes/unity.c: functions accurate near unity: cephes_log(x) = log(1 + x) (log1p), cephes_exp(x) =
   exp(x) - 1 (expm1) and cosm1(x) = cos(x) - 1, each by a rational approximation on the interval where the plain
   computation would lose digits and by the library function elsewhere. (cephes_log and cephes_exp are gretl's
   renamings of the cephes log1p and expm1.)

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier.

   Only bdtr.d calls cephes_exp and cephes_log, and nothing calls cosm1, so all three have D linkage. SQRTH and
   SQRT2 are macros in this file (literals, shadowing the constants of constants.d) and are enums here; PIO4 is
   the cephes constant. The coefficient arrays are private __gshared. polevl and p1evl come from cephes.polevl. */

module cephes.unity;

import core.stdc.math : cos, exp, isfinite, isnan, log;

import cephes.constants : PIO4;
import cephes.polevl : polevl, p1evl;

/* log1p(x) = log(1 + x)  */

/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 * Theoretical peak relative error = 2.32e-20
 */

private __gshared double[7] LP = [
 4.5270000862445199635215E-5,
 4.9854102823193375972212E-1,
 6.5787325942061044846969E0,
 2.9911919328553073277375E1,
 6.0949667980987787057556E1,
 5.7112963590585538103336E1,
 2.0039553499201281259648E1];

private __gshared double[6] LQ = [
/* 1.0000000000000000000000E0,*/
 1.5062909083469192043167E1,
 8.3047565967967209469434E1,
 2.2176239823732856465394E2,
 3.0909872225312059774938E2,
 2.1642788614495947685003E2,
 6.0118660497603843919306E1];

private enum SQRTH = 0.70710678118654752440;
private enum SQRT2 = 1.41421356237309504880;

double cephes_log(double x)
{
    double z = 1.0 + x;

    if (z < SQRTH || z > SQRT2) {
	return log(z);
    }

    z = x*x;
    z = -0.5 * z + x * (z * polevl(x, LP.ptr, 6) / p1evl(x, LQ.ptr, 6));

    return x + z;
}

/* expm1(x) = exp(x) - 1  */

/*  e^x =  1 + 2x P(x^2)/( Q(x^2) - P(x^2) )
 * -0.5 <= x <= 0.5
 */

private __gshared double[3] EP = [
 1.2617719307481059087798E-4,
 3.0299440770744196129956E-2,
 9.9999999999999999991025E-1];

private __gshared double[4] EQ = [
 3.0019850513866445504159E-6,
 2.5244834034968410419224E-3,
 2.2726554820815502876593E-1,
 2.0000000000000000000897E0];

double cephes_exp(double x)
{
    double r, xx;

    if (isnan(x)) {
	return x;
    }

    if (!isfinite(x)) {
	return (x < 0)? -1.0 : x;
    }

    if (x < -0.5 || x > 0.5) {
	return exp(x) - 1.0;
    }

    xx = x * x;
    r = x * polevl(xx, EP.ptr, 2);
    r = r/(polevl(xx, EQ.ptr, 3) - r);

    return r + r;
}

/* cosm1(x) = cos(x) - 1  */

private __gshared double[7] coscof = [
 4.7377507964246204691685E-14,
-1.1470284843425359765671E-11,
 2.0876754287081521758361E-9,
-2.7557319214999787979814E-7,
 2.4801587301570552304991E-5,
-1.3888888888888872993737E-3,
 4.1666666666666666609054E-2];

double cosm1(double x)
{
    double xx;

    if (x < -PIO4 || x > PIO4) {
	return cos(x) - 1.0;
    }

    xx = x * x;
    xx = -0.5*xx + xx * xx * polevl(xx, coscof.ptr, 6);

    return xx;
}
