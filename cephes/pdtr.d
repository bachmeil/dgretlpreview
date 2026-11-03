/* Port of cephes/pdtr.c: the Poisson distribution. pdtr(k, m) is the sum of the first k + 1 terms of the Poisson
   distribution with mean m (computed as the right tail of the incomplete gamma integral at k + 1), pdtrc(k, m) the
   sum of the terms k + 1 to infinity, and pdtri(k, y) the mean m for which pdtr(k, m) = y (by the inverse gamma
   integral).

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. igam and igamc come from cephes.igam and mtherr from
   cephes.mtherr. igami comes from cephes.igami. */

module cephes.pdtr;

import cephes.igam : igam, igamc;
import cephes.igami : igami;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double pdtrc(int k, double m)
{
    double v;

    if ((k < 0) || (m <= 0.0)) {
	mtherr("pdtrc", CephesError.domain);
	return 0.0;
    }

    v = k+1;
    return igam(v, m);
}

double pdtr(int k, double m)
{
    double v;

    if ((k < 0) || (m <= 0.0)) {
	mtherr("pdtr", CephesError.domain);
	return 0.0;
    }

    v = k+1;
    return igamc(v, m);
}

double pdtri(int k, double y)
{
    double v;

    if ((k < 0) || (y < 0.0) || (y >= 1.0)) {
	mtherr("pdtri", CephesError.domain);
	return 0.0;
    }

    v = k+1;
    v = igami(v, y);
    return v;
}
