/* Port of cephes/btdtr.c: the beta distribution function btdtr(a, b, x), the area under the beta density with
   parameters a and b from 0 to x; it is the incomplete beta integral, incbet.

   Cephes Math Library Release 2.8, June 2000, Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier.

   No C file calls it, so it has D linkage. incbet comes from cephes.incbet. */

module cephes.btdtr;

import cephes.incbet : incbet;

double btdtr(double a, double b, double x)
{
    return incbet(a, b, x);
}
