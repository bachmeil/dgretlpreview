/* Port of cephes/chdtr.c: the chi-square distribution. chdtr(df, x) is the area under the chi-square density with
   df degrees of freedom from 0 to x, chdtrc(df, x) the area from x to infinity, and chdtri(df, y) the x for which
   the upper tail area is y. Computed through the incomplete gamma integral and its inverse.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. igam and igamc come from cephes.igam. mtherr comes
   from cephes.mtherr. igami comes from cephes.igami. */

module cephes.chdtr;

import cephes.igam : igam, igamc;
import cephes.igami : igami;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double chdtrc(double df, double x)
{
    if (x < 0.0 || df < 1.0) {
	mtherr("chdtrc", CephesError.domain);
	return 0.0;
    }

    return igamc(df/2.0, x/2.0);
}

double chdtr(double df, double x)
{
    if (x < 0.0 || df < 1.0) {
	mtherr("chdtr", CephesError.domain);
	return 0.0;
    }

    return igam(df/2.0, x/2.0);
}

double chdtri(double df, double y)
{
    double x;

    if (y < 0.0 || y > 1.0 || df < 1.0) {
	mtherr("chdtri", CephesError.domain);
	return 0.0;
    }

    x = igami(0.5 * df, y);

    return 2.0 * x;
}
