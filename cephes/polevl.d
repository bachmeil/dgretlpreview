/* Port of cephes/polevl.c: polynomial evaluation by Horner's rule. polevl(x, coef, N) evaluates the polynomial
   of degree N whose N + 1 coefficients are stored highest order first; p1evl(x, coef, N) does the same for a
   polynomial whose leading coefficient is 1 and not stored (so it has N coefficients and degree N). Used by almost
   every cephes function with a rational or polynomial approximation.

   Cephes Math Library Release 2.1, December 1988, Copyright 1984, 1987, 1988 by Stephen L. Moshier.

   The pointer walking of the C (*p++) is kept, and with it the do-while loops that decrement their counter before
   testing it: polevl with N = 0 and p1evl with N = 1 (or less) never stop and read past the array. All callers pass
   N >= 1 for polevl and N >= 2 for p1evl. The coefficient parameters are const(double) *, where the C had double *.
   Only D code calls these functions, so they have D linkage. */

module cephes.polevl;

double polevl(double x, const(double) * coef, int N)
{
    double ans;
    int i;
    const(double) * p;

    p = coef;
    ans = *p++;
    i = N;

    do {
        ans = ans * x + *p++;
    } while (--i);

    return ans;
}

double p1evl(double x, const(double) * coef, int N)
{
    double ans;
    const(double) * p;
    int i;

    p = coef;
    ans = x + *p++;
    i = N-1;

    do {
        ans = ans * x + *p++;
    } while (--i);

    return ans;
}
