/* Port of cephes/gdtr.c: the gamma distribution. gdtr(a, b, x) is the integral from 0 to x of a^b t^(b-1)
   exp(-a t) / Gamma(b) (the gamma density with scale 1/a and shape b), computed as the incomplete gamma integral
   igam(b, a x); gdtrc(a, b, x) is the integral from x to infinity, igamc(b, a x).

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. igam and igamc come from cephes.igam and mtherr from
   cephes.mtherr. */

module cephes.gdtr;

import cephes.igam : igam, igamc;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double gdtr(double a, double b, double x)
{
    if (x < 0.0) {
	mtherr("gdtr", CephesError.domain);
	return 0.0;
    }
    return igam(b, a * x);
}

double gdtrc(double a, double b, double x)
{
    if (x < 0.0) {
	mtherr("gdtrc", CephesError.domain);
	return 0.0;
    }
    return igamc(b, a * x);
}
