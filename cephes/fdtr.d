/* Port of cephes/fdtr.c: the F distribution. fdtr(a, b, x) is the area under the F density (Snedecor's density,
   the variance-ratio density) from 0 to x with a and b degrees of freedom, fdtrc(a, b, x) the area from x to
   infinity, and fdtri(a, b, y) the x for which that upper tail area is y. All three are computed through the
   incomplete beta integral.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls these functions, so they have D linkage. incbet comes from cephes.incbet. incbi comes from
   cephes.incbi. mtherr comes from cephes.mtherr. */

module cephes.fdtr;

import cephes.incbet : incbet;
import cephes.incbi : incbi;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr;

double fdtrc(double a, double b, double x)
{
    double w;

    if (a < 1 || b < 1 || x < 0.0) {
        mtherr("fdtrc", CephesError.domain);
        return 0.0;
    }

    w = b / (b + a * x);

    return incbet(0.5*b, 0.5*a, w);
}

double fdtr(double a, double b, double x)
{
    double w;

    if (a < 1 || b < 1 || x < 0.0) {
        mtherr("fdtr", CephesError.domain);
        return 0.0;
    }

    w = a * x;
    w = w / (b + w);

    return incbet(0.5*a, 0.5*b, w);
}

double fdtri(double a, double b, double y)
{
    double w, x;

    if (a < 1 || b < 1 || y <= 0.0 || y > 1.0) {
        mtherr("fdtri", CephesError.domain);
        return 0.0;
    }

    /* Compute probability for x = 0.5 */
    w = incbet(0.5*b, 0.5*a, 0.5);

    /* If that is greater than y, then the solution w < .5.
       Otherwise, solve at 1-y to remove cancellation in (b - b*w).
    */
    if (w > y || y < 0.001) {
        w = incbi(0.5*b, 0.5*a, y);
        x = (b - b*w)/(a*w);
    } else {
        w = incbi(0.5*a, 0.5*b, 1.0-y);
        x = b*w/(a*(1.0-w));
    }

    return x;
}
