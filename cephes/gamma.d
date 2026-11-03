/* Port of cephes/gamma.c: the gamma function (cephes_gamma; the cephes gamma() renamed) and the natural
   logarithm of the absolute value of the gamma function (lgam, with the alias cephes_lgamma). The sign of the
   gamma function is left in the global cephes_sgngam, read through get_cephes_sgngam().

   Arguments |x| <= 34 are reduced by recurrence to a rational approximation on (2,3); large arguments use
   Stirling's formula and large negative ones a reflection formula. Cephes Math Library Release 2.8, June 2000,
   Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier.

   The functions and the sign variable have D linkage. mtherr_with_arg comes from cephes.mtherr. polevl and p1evl
   come from cephes.polevl. The scalar constants are __gshared variables, not enums, so that they hold doubles as
   in the C (see constants.d). */

module cephes.gamma;

import core.stdc.math : exp, fabs, floor, isfinite, isnan, log, pow, sin;

import cephes.constants : PI;
import cephes.libprobh : CephesError;
import cephes.mtherr : mtherr_with_arg;
import cephes.polevl : polevl, p1evl;

private immutable double[7] P = [
    1.60119522476751861407E-4,
    1.19135147006586384913E-3,
    1.04213797561761569935E-2,
    4.76367800457137231464E-2,
    2.07448227648435975150E-1,
    4.94214826801497100753E-1,
    9.99999999999999996796E-1
];

private immutable double[8] Q = [
    -2.31581873324120129819E-5,
    5.39605580493303397842E-4,
    -4.45641913851797240494E-3,
    1.18139785222060435552E-2,
    3.58236398605498653373E-2,
    -2.34591795718243348568E-1,
    7.14304917030273074085E-2,
    1.00000000000000000320E0
];

private enum MAXGAM = 171.624376956302725;

private __gshared double LOGPI = 1.14472988584940017414;

/* Stirling's formula for the gamma function */

private immutable double[5] STIR = [
    7.87311395793093628397E-4,
    -2.29549961613378126380E-4,
    -2.68132617805781232825E-3,
    3.47222221605458667310E-3,
    8.33333333333482257126E-2,
];

private enum MAXSTIR = 143.01608;

private __gshared double SQTPI = 2.50662827463100050242E0;

__gshared int cephes_sgngam = 0;

/* accessor for sgngam */

int get_cephes_sgngam()
{
    return cephes_sgngam;
}

/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */

private double stirf(double x)
{
    double y, w, v;

    w = 1.0/x;
    w = 1.0 + w * polevl(w, STIR.ptr, 4);
    y = exp(x);

    if (x > MAXSTIR) {
        /* Avoid overflow in pow() */
        v = pow(x, 0.5 * x - 0.25);
        y = v * (v / y);
    } else {
        y = pow(x, x - 0.5) / y;
    }

    y = SQTPI * y * w;

    return y;
}

double cephes_gamma(double x)
{
    double p, q, z;
    int i;

    cephes_sgngam = 1;

    if (isnan(x)) {
        return x;
    }

    if (!isfinite(x)) {
        return (x < 0)? double.nan : x;
    }

    q = fabs(x);

    if (q > 33.0) {
        if (x < 0.0) {
            p = floor(q);
            if (p == q) {
                goto gamnan;
            }
            i = cast(int) p;
            if ((i & 1) == 0) {
                cephes_sgngam = -1;
            }
            z = q - p;
            if (z > 0.5) {
                p += 1.0;
                z = q - p;
            }
            z = q * sin(PI * z);
            if (z == 0.0) {
                return cephes_sgngam * double.infinity;
            }
            z = fabs(z);
            z = PI/(z * stirf(q));
        } else {
            z = stirf(x);
        }
        return cephes_sgngam * z;
    }

    z = 1.0;
    while (x >= 3.0) {
        x -= 1.0;
        z *= x;
    }

    while (x < 0.0) {
        if (x > -1.0e-9) {
            goto small;
        }
        z /= x;
        x += 1.0;
    }

    while (x < 2.0) {
        if (x < 1.0e-9) {
            goto small;
        }
        z /= x;
        x += 1.0;
    }

    if (x == 2.0) {
        return z;
    }

    x -= 2.0;
    p = polevl(x, P.ptr, 6);
    q = polevl(x, Q.ptr, 7);

    return z * p / q;

 small:

    if (x == 0.0) {
        goto gamnan;
    } else {
        return z/((1.0 + 0.5772156649015329 * x) * x);
    }

 gamnan:
    mtherr_with_arg("gamma", CephesError.domain, x);
    return double.nan;
}

/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */

private immutable double[5] A = [
    8.11614167470508450300E-4,
    -5.95061904284301438324E-4,
    7.93650340457716943945E-4,
    -2.77777777730099687205E-3,
    8.33333333333331927722E-2
];

private immutable double[6] B = [
    -1.37825152569120859100E3,
    -3.88016315134637840924E4,
    -3.31612992738871184744E5,
    -1.16237097492762307383E6,
    -1.72173700820839662146E6,
    -8.53555664245765465627E5
];

private immutable double[6] C = [
    /* 1.00000000000000000000E0, */
    -3.51815701436523470549E2,
    -1.70642106651881159223E4,
    -2.20528590553854454839E5,
    -1.13933444367982507207E6,
    -2.53252307177582951285E6,
    -2.01889141433532773231E6
];

private __gshared double LS2PI = 0.91893853320467274178;

private enum MAXLGM = 2.556348e305;

/* Logarithm of gamma function */

double lgam(double x)
{
    double p, q, u, w, z;
    int i;

    cephes_sgngam = 1;

    if (isnan(x)) {
        return x;
    }

    if (!isfinite(x)) {
        return double.infinity;
    }

    if (x < -34.0) {
        q = -x;
        w = lgam(q); /* note this modifies cephes_sgngam! */
        p = floor(q);
        if (p == q) {
            goto lgsing;
        }
        i = cast(int) p;
        if ((i & 1) == 0) {
            cephes_sgngam = -1;
        } else {
            cephes_sgngam = 1;
        }
        z = q - p;
        if (z > 0.5) {
            p += 1.0;
            z = p - q;
        }
        z = q * sin(PI * z);
        if (z == 0.0) {
            goto lgsing;
        }
        z = LOGPI - log(z) - w;
        return z;
    }

    if (x < 13.0) {
        z = 1.0;
        p = 0.0;
        u = x;
        while (u >= 3.0) {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while (u < 2.0) {
            if (u == 0.0) {
                goto lgsing;
            }
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if (z < 0.0) {
            cephes_sgngam = -1;
            z = -z;
        } else {
            cephes_sgngam = 1;
        }
        if (u == 2.0) {
            return log(z);
        }
        p -= 2.0;
        x = x + p;
        p = x * polevl(x, B.ptr, 5) / p1evl(x, C.ptr, 6);
        return log(z) + p;
    }

    if (x > MAXLGM) {
        return cephes_sgngam * double.infinity;
    }

    q = (x - 0.5) * log(x) - x + LS2PI;
    if (x > 1.0e8) {
        return q;
    }

    p = 1.0/(x*x);
    if (x >= 1000.0) {
        q += ((7.9365079365079365079365e-4 * p
               - 2.7777777777777777777778e-3) * p
              + 0.0833333333333333333333) / x;
    } else {
        q += polevl(p, A.ptr, 4) / x;
    }

    return q;

 lgsing:
    mtherr_with_arg("lgam", CephesError.sing, x);
    return double.infinity;
}

double cephes_lgamma(double x)
{
    return lgam(x);
}
