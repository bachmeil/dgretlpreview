/* Port of cephes/nbdtr.c: the negative binomial distribution. nbdtr(k, n, p) is the sum of the terms 0 through k of
   the negative binomial distribution with n successes and success probability p, nbdtrc(k, n, p) the sum of the
   terms k + 1 to infinity, and nbdtri(k, n, p) the success probability for which nbdtr(k, n, .) equals p. Computed
   through the incomplete beta integral and its inverse.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No file calls these functions, so they have D linkage. incbet and incbi come from the ported modules; mtherr
   comes from cephes.mtherr. */

module cephes.nbdtr;

import cephes.incbet : incbet;
import cephes.incbi : incbi;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double nbdtrc(int k, int n, double p)
{
    double dk, dn;

    if ((p < 0.0) || (p > 1.0)) {
	goto domerr;
    }

    if (k < 0) {
    domerr:
	mtherr("nbdtr", CephesError.domain);
	return 0.0;
    }

    dk = k+1;
    dn = n;
    return incbet(dk, dn, 1.0 - p);
}

double nbdtr(int k, int n, double p)
{
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
	goto domerr;

    if (k < 0) {
    domerr:
	mtherr("nbdtr", CephesError.domain);
	return 0.0;
    }

    dk = k+1;
    dn = n;
    return incbet(dn, dk, p);
}

double nbdtri(int k, int n, double p)
{
    double dk, dn, w;

    if ((p < 0.0) || (p > 1.0))
	goto domerr;

    if (k < 0) {
    domerr:
	mtherr("nbdtri", CephesError.domain);
	return 0.0;
    }

    dk = k+1;
    dn = n;
    w = incbi(dn, dk, p);

    return w;
}
