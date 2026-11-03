/* Port of pvalues.c */

module gretl.pvalues;


import core.stdc.ctype : tolower;
import core.stdc.errno : errno;
import core.stdc.math : ceil, erfc, exp, fabs, floor, fmin, FP_INFINITE, FP_NAN,
    log, log2, pow, sqrt, trunc;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : free, malloc;
import core.stdc.string : strcat, strcmp, strcpy, strlen, strncmp;

import gretl.gretl_errorsh;
import gretl.gretl_matrix : gretl_matrix_alloc;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.libgretlh;
import cephes.mtherr : get_cephes_errno;
import cephes.ndtr : ndtr;
import cephes.stdtr : stdtr, stdtri;
import cephes.ndtri : ndtri;
import cephes.pdtr : pdtr, pdtrc;
import cephes.psi : psi;
import cephes.gamma : cephes_gamma, cephes_lgamma;
import cephes.bdtr : bdtr, bdtrc;
import cephes.btdtr : btdtr;
import cephes.chdtr : chdtr, chdtrc, chdtri;
import cephes.fdtr : fdtr, fdtrc, fdtri;
import cephes.gdtr : gdtr, gdtrc;
import cephes.igami : igami;
import cephes.hyp2f1 : hyp2f1;
import cephes.incbet : incbet;
import cephes.libprobh : CephesError;
import gretl.missingh;
import gretl.pvaluesh : DistCode;
import gretl.strutils : gretl_string_split, strings_array_free;
import gretl.plugins : get_plugin_function;
import gretl.genmain : generate_scalar;
import gretl.genfuncs : gretl_bessel, logistic_cdf;
import gretl.random : gretl_rand_F, gretl_rand_GED, gretl_rand_beta, gretl_rand_beta_binomial, gretl_rand_binomial, gretl_rand_chisq, gretl_rand_exponential, gretl_rand_gamma, gretl_rand_laplace, gretl_rand_logistic, gretl_rand_normal_full, gretl_rand_poisson, gretl_rand_student, gretl_rand_tnormal, gretl_rand_uniform_int_minmax, gretl_rand_uniform_minmax, gretl_rand_weibull;

/* plugin entry points, fetched at run time with get_plugin_function()
   (the plugins are still C, hence extern (C)) */

private alias TracePvalueFn = extern (C) double function(double, int, int, int);
private alias DwLookupFn = extern (C) int function(int, int, GretlMatrixRaw **);

private enum NORM_CDF_MAX = 0.9999999999999999;

/**
 * gammafun:
 * @x: argument.
 *
 * Returns: the gamma function of @x, or #NADBL on failure.
 */

double gammafun(double x)
{
    double ret = cephes_gamma(x);

    if (get_cephes_errno()) {
        ret = NADBL;
    }

    return ret;
}

/**
 * lngamma:
 * @x: argument.
 *
 * Returns: the log gamma function of @x, or #NADBL on failure.
 */

double lngamma(double x)
{
    double ret = cephes_lgamma(x);

    if (get_cephes_errno()) {
        ret = NADBL;
    }

    return ret;
}

/**
 * lnmgamma:
 * @x: argument
 * @p: shape parameter.
 *
 * Returns: the log multivariate gamma function of @x, or #NADBL
 * on failure. Note that the shape parameter @p must be at least
 * 1 and no greater than 2 times the argument @x, otherwise
 * NODBL is returned.
 */

double lnmgamma(double x, int p)
{
    double ret = NADBL;

    if (p >= 1 && p <= 2*x) {
        int i;

        ret = lngamma(x);
        for (i=1; i<p; i++) {
            ret += lngamma(x - 0.5*i);
        }
        ret += p*(p-1) * LN_PI / 4.0;
    }

    return ret;
}

/**
 * digamma:
 * @x: argument.
 *
 * Returns: the digamma (or Psi) function of @x, or #NADBL on failure.
 */

double digamma(double x)
{
    double ret = psi(x);

    if (get_cephes_errno()) {
        ret = NADBL;
    }

    return ret;
}

/**
 * trigamma:
 * @x: argument.
 *
 * Returns: the trigamma function of @x, or #NADBL on failure. The code
 * is adapted from
 * https://people.sc.fsu.edu/~jburkardt/f_src/asa121/asa121.html
 * See BE Schneider, Algorithm AS 121: Trigamma Function.
 * Applied Statistics, Volume 27, Number 1, pages 97-99, 1978.
 *
 * The main modification with respect to the published version is the
 * addition of three extra terms to the asymptotic expansion for x >= B.
 */

double trigamma(double x)
{
    double ret = 0;
    double A = 0.000001; /* threshold for "small" argument */
    double B = 5.0;      /* threshold for "large" argument */

    /* the Bernoulli numbers */
    double b2 =  1.0/6;
    double b4 = -1.0/30;
    double b6 =  1.0/42;
    double b8 = -1.0/30;
    double b10 = 5.0/66;
    double b12 = -691.0/2730;
    double b14 = 7.0/6;

    if (x <= 0) {
        ret = NADBL;
    } else if (x < A) {
        ret = 1.0 / (x * x);
    } else {
        double y, a1, a2, z = x;

        while (z < B) {
            ret += 1.0 / (z * z);
            z += 1.0;
        }
        /* Apply asymptotic formula for argument >= B */
        y = 1.0 / (z * z);
        a1 = 0.5 * y;
        a2 = (1 + y *
              (b2 + y *
               (b4 + y *
                (b6 + y *
                 (b8 + y *
                  (b10 + y *
                   (b12 + y * b14))))))) / z;
        ret += a1 + a2;
    }

    return ret;
}

/**
 * hypergeo:
 * @a: argument.
 * @b: argument.
 * @c: argument.
 * @x: absolute value must be less than 1.0.
 *
 * Returns: the Gauss hypergeometric function 2F1 of
 * the given arguments.
 */

double hypergeo(double a, double b, double c, double x)
{
    double ret = hyp2f1(a, b, c, x);

    if (get_cephes_errno()) {
        ret = NADBL;
    }

    return ret;
}
/**
 * beta_cdf:
 *
 * Returns the probability that a B(a,b) random variable is
 * between 0 and z, or #NADBL on failure.
 */

double beta_cdf(double a, double b, double z)
{
    double x = NADBL;

    if (a > 0 && b > 0 && z >= 0 && z <= 1) {
        x = btdtr(a, b, z);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/**
 * binomial_cdf:
 * @p: probability of success on each trial.
 * @n: number of trials.
 * @k: maximum number of successes.
 *
 * Returns: the probability of @k or less successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf(double p, int n, int k)
{
    double x = NADBL;

    if (p >= 0 && n >= 0 && k >= 0) {
        x = bdtr(k, n, p);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/**
 * binomial_cdf_comp:
 * @p: probability of success on each trial.
 * @n: number of trials.
 * @k: maximum number of successes.
 *
 * Returns: the probability of @k + 1 or more successes on
 * @n trials given binomial probability @p, or
 * #NADBL on failure.
 */

double binomial_cdf_comp(double p, int n, int k)
{
    double x = NADBL;

    if (p >= 0 && n >= 0 && k >= 0) {
        x = bdtrc(k, n, p);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/* Binomial inverse CDF via binary search */

private int bininv_binary(double p, int n, double u)
{
    int klo = 0;
    int khi = n;
    int kmid, ke, ke1;
    double f, f1;
    double fe, fe1;
    double mu, pdmax;
    int evals = 0;
    int nk, emax;
    int ret = -1;

    if (u > 0.5) {
        klo = cast(int) floor(n*p) - 1;
    } else if (u < 0.5) {
        khi = cast(int) floor(n*p) + 1;
    }

    kmid = (klo + khi) / 2;
    nk = khi - klo + 1;
    emax = cast(int) (2 * log2(cast(double) nk));

    /* Find the max step-up of the CDF, which occurs
       from expected value minus 1 to expected value.
    */
    mu = n*p;
    ke = cast(int) (u > 0.5 ? ceil(mu) : floor(mu));
    if (ke > 0) {
        ke1 = ke - 1;
        fe = binomial_cdf(p, n, ke);
        fe1 = binomial_cdf(p, n, ke1);
        pdmax = fe - fe1;
    } else {
        pdmax = 0;
        ke = ke1 = -1;
        fe = fe1 = 0;
    }

    while (evals < emax) {
        if (kmid == ke || kmid == ke1) {
            f = (kmid == ke)? fe : fe1;
        } else {
            f = binomial_cdf(p, n, kmid);
            evals++;
        }
        if (f < u) {
            /* we need to look on the right */
            klo = kmid;
            kmid = (klo + khi) / 2;
            if (kmid == klo) {
                ret = kmid;
                break;
            }
        } else if (f == u) {
            /* on the nail (unlikely) */
            ret = kmid;
            break;
        } else {
            /* f > u */
            if (f - u <= pdmax && kmid > 0) {
                /* check out k - 1 */
                int k1 = kmid - 1;

                if (k1 == ke || k1 == ke1) {
                    f1 = (k1 == ke)? fe : fe1;
                } else {
                    f1 = binomial_cdf(p, n, k1);
                    evals++;
                }
                if (f1 < u) {
                    ret = kmid;
                    break;
                }
            }
            /* we need to look on the left */
            khi = kmid;
            kmid = (klo + khi) / 2;
            if (kmid == klo) {
                ret = kmid;
                break;
            }
        }
    }

    if (ret < 0) {
        fprintf(stderr, "bininv_binary, bad: evals=%d, p=%g, n=%d, u=%g\n",
                evals, p, n, u);
    }

    return ret;
}

/* Binomial inverse CDF via bottom-up or top-down summation */

private int bininv_sum(double p, int n, double u)
{
    double f, a, s;
    int k = 1;

    f = binomial_cdf(p, n, n/2);

    if (u <= f) {
        /* bottom-up */
        a = pow(1 - p, n);
        s = a - u;
        while (s < 0) {
            a = (a*p/(1 - p)) * (n - k + 1)/k;
            s += a;
            k++;
        }
        return k - 1;
    } else {
        /* top-down */
        a = pow(p, n);
        s = 1 - u - a;
        while (s >= 0) {
            a = (a*(1 - p)/p) * (n - k + 1)/k;
            s -= a;
            k++;
        }
        return n - k + 1;
    }
}

/* Returns the smallest x such that the probability of obtaining
   x successes on @n trials with success probability @p is at least
   the given cumulative probability @u. FIXME efficiency?
*/

private double binomial_cdf_inverse(double p, int n, double u)
{
    double x = NADBL;
    int k;

    if (u <= 0 || u > 1 || n <= 0 || p <= 0 || p >= 1) {
        ;
    } else if (u == 1.0) {
        x = n;
    } else if (n < 500) {
        k = bininv_sum(p, n, u);
        if (k >= 0) {
            x = k;
        }
    } else {
        k = bininv_binary(p, n, u);
        if (k >= 0) {
            x = k;
        }
    }

    return x;
}

/* The following is no doubt horribly inefficient */

private double binomial_critval(double p, int n, double a)
{
    double pk, ac = 1 - a;
    int k;

    if (n <= 0 || p <= 0 || p >= 1 || a <= 0 || a >= 1) {
        return NADBL;
    }

    for (k=n; k>0; k--) {
        pk = binomial_cdf(p, n, k);
        if (pk < ac) {
            break;
        }
    }

    return cast(double) (k + 1);
}

/**
 * binomial_pmf:
 * @p: success probability.
 * @n: number of trials.
 * @k: number of successes.
 *
 * Returns: the probability mass for @k successes in @n
 * binomial trials with success probability @p.
 */

double binomial_pmf(double p, int n, int k)
{
    double pm = binomial_cdf(p, n, k);

    if (k > 0 && !na(pm)) {
        pm -= binomial_cdf(p, n, k - 1);
    }

    return pm;
}

private int binomial_pmf_array(double p, int n, double * x, int T)
{
    double pm;
    int i, k;

    for (i=0; i<T; i++) {
        k = cast(int) x[i];
        pm = binomial_cdf(p, n, k);
        if (k > 0 && !na(pm)) {
            pm -= binomial_cdf(p, n, k - 1);
        }
        x[i] = pm;
    }

    return 0;
}

/**
 * x_factorial:
 * @x: input value.
 *
 * Returns: the factorial of int(x), cast to a double, or
 * NADBL on failure.
 */

double x_factorial(double x)
{
    double fact;
    int n = cast(int) x;

    if (x < 0) {
        fact = NADBL;
    } else if (x > 12) {
        fact = cephes_gamma(1 + x);
        if (get_cephes_errno()) {
            fact = NADBL;
        }
    } else if (n == 0) {
        fact = 1;
    } else {
        fact = n;
        while (--n > 1) {
            fact *= n;
        }
    }

    return fact;
}

/**
 * log_x_factorial:
 * @x: input value.
 *
 * Returns: the log of the factorial of int(x), cast to a double, or
 * NADBL on failure.
 */

double log_x_factorial(double x)
{
    double lfact;
    int n = cast(int) x;

    if (x < 0) {
        lfact = NADBL;
    } else if (x > 12) {
        lfact = cephes_lgamma(1 + x);
        if (get_cephes_errno()) {
            lfact = NADBL;
        }
    } else if (n == 0) {
        lfact = 0;
    } else {
        lfact = n;
        while (--n > 1) {
            lfact *= n;
        }
        lfact = log(lfact);
    }

    return lfact;
}

/**
 * tcrit95:
 * @df: degrees of freedom.
 *
 * Returns: the two-sided 95 percent critical value for the t
 * distribution with @df degrees of freedom, or #NADBL on
 * failure.
 */

double tcrit95(int df)
{
    double x = NADBL;

    if (df > 0) {
        x = stdtri(df, 0.975);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/**
 * rhocrit:
 * @n: sample size.
 * @alpha: significance level as decimal fraction.
 *
 * Computes the two-sided 100 * @alpha critical value of the sample
 * correlation coefficient for a sample of size @n. This is based
 * on the inverse of the function which maps from the correlation
 * coefficient, r, to a student t statistic, namely
 *
 *            t = r / sqrt[(1—r^2) / (n-2)]
 *
 * The inverse is r = sqrt(t^2 / (t^2 + n - 2)).
 *
 * Returns: the critical value, or #NADBL on failure.
 */

double rhocrit(int n, double alpha)
{
    double rc = NADBL;

    if (n - 2 > 0) {
        double tmax = 1.0 - 0.5 * alpha;
        double tc = stdtri(n - 2, tmax);

        if (get_cephes_errno() == 0) {
            double tc2 = tc * tc;

            rc = sqrt(tc2 / (tc2 + n - 2));
        }
    }

    return rc;
}

/**
 * normal_pvalue_2:
 * @x: double-precision value.
 *
 * Calculates the two-sided p-value for @x in relation to the
 * standard normal distribution.
 *
 * Returns: 2 times (1 minus the value of the standard normal
 * CDF evaluated at abs(@x)), or 0 on underflow.
 */

double normal_pvalue_2(double x)
{
    double p = (x < 0)? ndtr(x) : ndtr(-x);

    return 2 * p;
}

/**
 * normal_pvalue_1:
 * @x: double-precision value.
 *
 * Calculates the one-sided p-value for @x in relation to the
 * standard normal distribution (that is, the probability that a
 * random variable distributed as N(0, 1) is greater than @x).
 *
 * Returns: 1 minus the value of the standard normal CDF
 * evaluated at @x.
 */

double normal_pvalue_1(double x)
{
    double p = ndtr(x);

    if (get_cephes_errno()) {
        p = NADBL;
    } else {
        p = 1 - p;
    }

    return p;
}

/**
 * student_cdf:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from minus infinity to @x of
 * the Student's t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

double student_cdf(double df, double x)
{
    double p = NADBL;

    if (df > 0) {
        p = stdtr(df, x);
        if (get_cephes_errno()) {
            p = NADBL;
        }
    }

    return p;
}

/**
 * student_cdf_comp:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from @x to infinity of
 * the t distribution with @df degrees of freedom, or
 * #NADBL on failure.
 */

private double student_cdf_comp(double df, double x)
{
    double p = NADBL;

    if (df > 0) {
        if (x > 0) {
            p = stdtr(df, -x);
            if (get_cephes_errno()) {
                p = NADBL;
            }
        } else {
            p = stdtr(df, x);
            if (get_cephes_errno()) {
                p = NADBL;
            } else {
                p = 1 - p;
            }
        }
    }

    return p;
}

/**
 * normal_cdf_comp:
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from @x to infinity of the standard
 * normal distribution, or #NADBL on failure.
 */

double normal_cdf_comp(double x)
{
    double p;

    if (x > 0) {
        p = ndtr(-x);
        if (get_cephes_errno()) {
            p = NADBL;
        }
    } else {
        p = ndtr(x);
        if (get_cephes_errno()) {
            p = NADBL;
        } else {
            p = 1 - p;
        }
    }

    return p;
}

/**
 * student_pvalue_2:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the probability that t(@df) is greater than @x
 * (two-sided, using the absolute value of @x), or
 * #NADBL on failure.
 */

double student_pvalue_2(double df, double x)
{
    double p = NADBL;

    if (df > 0) {
        if (x < 0.0) {
            p = stdtr(df, x);
        } else {
            p = stdtr(df, -x);
        }
        if (get_cephes_errno()) {
            p = NADBL;
        } else {
            p *= 2;
        }
    }

    return p;
}

/**
 * student_pvalue_1:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the probability that t(@df) is greater than @x,
 * or #NADBL on failure.
 */

double student_pvalue_1(double df, double x)
{
    double p = NADBL;

    if (df > 0) {
        p = stdtr(df, x);
        if (get_cephes_errno()) {
            p = NADBL;
        } else {
            p = (p >= 1.0)? 0.0 : 1 - p;
        }
    }

    return p;
}

/**
 * student_critval:
 * @df: degrees of freedom.
 * @a: right-tail probability.
 *
 * Returns: the argument x such that the integral from x to
 * infinity of the t(@df) density is equal to the given
 * probability @a, or #NADBL on failure.
 */

double student_critval(double df, double a)
{
    double x;

    if (df < 0) {
        return NADBL;
    }

    if (a > .10) {
        x = stdtri(df, 1 - a);
    } else {
        x = -stdtri(df, a);
    }

    if (get_cephes_errno()) {
        x = NADBL;
    }

    return x;
}

/**
 * student_cdf_inverse:
 * @df: degrees of freedom.
 * @a: probability.
 *
 * Returns: the argument x such that the integral from
 * minus infinity to @x of the t(@df) density is equal to
 * the given probability @a, or #NADBL on failure.
 */

double student_cdf_inverse(double df, double a)
{
    double x;

    if (df < 0) {
        return NADBL;
    }

    x = stdtri(df, a);

    if (get_cephes_errno()) {
        x = NADBL;
    }

    return x;
}

/**
 * chisq_cdf:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from 0 to @x of the chi-square
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf(double df, double x)
{
    double p = chdtr(df, x);

    /* 2016-07-22: should we return NA instead of zero
       in case of CEPHES_UNDERFLOW? */
    if (get_cephes_errno() == CephesError.domain) {
        p = NADBL;
    }

    return p;
}

/**
 * chisq_cdf_comp:
 * @df: degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral from @x to infinity of the chi-square
 * distribution with @df degrees of freedom, or #NADBL
 * on failure.
 */

double chisq_cdf_comp(double df, double x)
{
    double p = chdtrc(df, x);

    if (get_cephes_errno() == CephesError.domain) {
        p = NADBL;
    }

    return p;
}

/**
 * chisq_critval:
 * @df: degrees of freedom.
 * @a: right-tail probability.
 *
 * Returns: the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given probability @a, or #NADBL on failure.
 */

private double chisq_critval(int df, double a)
{
    double x = NADBL;

    if (a == 1.0) {
        return 0.0;
    } else if (a == 0.0) {
        return 1.0/0.0;
    }

    x = chdtri(df, a);
    if (get_cephes_errno() == CephesError.domain) {
        x = NADBL;
    }

    return x;
}

private double chisq_cdf_inverse(int df, double a)
{
    double x = NADBL;

    if (a == 1.0) {
        return 1.0/0.0;
    } else if (a == 0.0) {
        return 0.0;
    }

    if (df >= 1 && a >= 0.0) {
        x = chdtri(df, 1 - a);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/**
 * snedecor_cdf:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from 0 to @x, or #NADBL on failure.
 */

double snedecor_cdf(double dfn, double dfd, double x)
{
    double p = NADBL;

    if (dfn > 0 && dfd > 0 && x >= 0) {
        p = fdtr(dfn, dfd, x);
        if (get_cephes_errno()) {
            p = NADBL;
        }
    }

    return p;
}

/**
 * snedecor_cdf_comp:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @x: the cutoff point in the distribution.
 *
 * Returns: the integral of the F distribution with @dfn and
 * @dfd degrees of freedom, from @x to infinity, or #NADBL
 * on failure.
 */

double snedecor_cdf_comp(double dfn, double dfd, double x)
{
    double p = NADBL;

    if (dfn > 0 && dfd > 0 && x >= 0) {
        p = fdtrc(dfn, dfd, x);
        if (get_cephes_errno()) {
            p = NADBL;
        }
    }

    return p;
}

/**
 * snedecor_pvalue_2:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @x: observed F-value.
 *
 * Returns: the two-sided p-value for the null hypothesis that
 * the population value of F equals 1.0, as in a test for
 * equality of variances, or #NADBL on failure.
 */

double snedecor_pvalue_2(double dfn, double dfd, double x)
{
    double p = NADBL;

    if (dfn > 0 && dfd > 0 && x >= 0) {
        p = 2.0 * fdtrc(dfn, dfd, x);
        if (x < 1) {
            p = 2.0 - p;
        }
        if (get_cephes_errno()) {
            p = NADBL;
        }
    }

    return p;
}

/**
 * snedecor_critval:
 * @dfn: numerator degrees of freedom.
 * @dfd: denominator degrees of freedom.
 * @a: right-tail probability.
 *
 * Returns: the F argument x such that the integral
 * from x to infinity of the F density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double snedecor_critval(double dfn, double dfd, double a)
{
    double x = NADBL;

    if (dfn > 0 && dfd > 0 && a >= 0) {
        x = fdtri(dfn, dfd, a);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

private double snedecor_cdf_inverse(double dfn, double dfd, double a)
{
    double x = NADBL;

    if (dfn > 0 && dfd > 0 && a >= 0) {
        x = fdtri(dfn, dfd, 1 - a);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/* PDFs */

private double Binv(double p, double q)
{
    double f = NADBL;

    errno = 0;

    if (p > 0 && q > 0) {
        double x1 = lngamma(p + q);
        double x2 = lngamma(p);
        double x3 = lngamma(q);

        if (!na(x1) && !na(x2) && !na(x3)) {
            f = exp(x1 - x2 - x3);
            if (errno) {
                f = NADBL;
            }
        }
    }

    return f;
}

private int snedecor_pdf_array(double v1, double v2, double * x, int n)
{
    int i, err = 0;

    errno = 0;

    if (v1 > 0 && v2 > 0) {
        double xm = v1, xn = v2;
        double vr = xm / xn;
        double x1 = Binv(0.5*xm, 0.5*xn);
        double x2 = pow(vr, 0.5*xm);
        double x3, x4;

        if (errno) {
            err = GretlError.nan;
        } else {
            for (i=0; i<n; i++) {
                if (!na(x[i]) && x[i] > 0) {
                    errno = 0;
                    x3 = pow(x[i], 0.5*xm - 1.0);
                    x4 = pow(1.0 + vr * x[i], 0.5 * (xm+xn));
                    x[i] = x1 * x2 * x3 / x4;
                    if (errno) {
                        x[i] = NADBL;
                    }
                } else {
                    x[i] = NADBL;
                }
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    errno = 0;

    return err;
}

private double snedecor_pdf(double m, double n, double x)
{
    snedecor_pdf_array(m, n, &x, 1);

    return x;
}

private int chisq_pdf_array(double m, double * x, int n)
{
    int i, err = 0;

    errno = 0;

    if (m > 0) {
        double m2 = m / 2.0;
        double x1 = pow(.5, m2);
        double x2 = gammafun(m2);
        double x3, x4;

        if (errno) {
            err = GretlError.nan;
        } else {
            for (i=0; i<n; i++) {
                if (!na(x[i]) && x[i] >= 0) {
                    errno = 0;
                    x3 = pow(x[i], m2 - 1.0);
                    x4 = exp(-x[i] / 2.0);
                    x[i] = (x1/x2) * x3 * x4;
                    if (errno) {
                        x[i] = NADBL;
                    }
                } else {
                    x[i] = NADBL;
                }
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    errno = 0;

    return err;
}

private double chisq_pdf(double m, double x)
{
    chisq_pdf_array(m, &x, 1);

    return x;
}

private int student_pdf_array(double m, double * x, int n)
{
    int i, err = 0;

    errno = 0;

    if (m > 0 && !na(m)) {
        double x1 = Binv(0.5 * m, 0.5) / sqrt(m);
        double x3 = 0.5 * (m + 1.0);
        double x2;

        if (errno) {
            err = GretlError.nan;
        } else {
            for (i=0; i<n; i++) {
                if (!na(x[i])) {
                    errno = 0;
                    x2 = m / (m + x[i] * x[i]);
                    x[i] = x1 * pow(x2, x3);
                    if (errno) {
                        x[i] = NADBL;
                    }
                }
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    errno = 0;

    return err;
}

private double student_pdf(double m, double x)
{
    student_pdf_array(m, &x, 1);

    return x;
}

private int weibull_pdf_array(double k, double l,
                              double * x, int n)
{
    int i, err = 0;

    errno = 0;

    if (!na(k) && k > 0 && !na(l) && l > 0) {
        double x1 = k / l;
        double x2, x3, x4;

        if (errno) {
            err = GretlError.nan;
        } else {
            for (i=0; i<n; i++) {
                if (!na(x[i]) && x[i] >= 0) {
                    errno = 0;
                    x3 = x[i] / l;
                    x2 = pow(x3, k - 1.0);
                    x3 *= x2;
                    x4 = exp(-x3);
                    x[i] = x1 * x2 * x4;
                    if (errno) {
                        x[i] = NADBL;
                    }
                } else {
                    x[i] = NADBL;
                }
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    errno = 0;

    return err;
}

private int exponential_pdf_array(double mu, double * x, int n)
{
    int i, err = 0;

    if (na(mu) || mu <= 0) {
        err = GretlError.invarg;
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    } else {
        double b = 1.0 / mu;

        for (i=0; i<n; i++) {
            if (x[i] < 0) {
                x[i] = 0;
            } else if (!na(x[i])) {
                x[i] = b * exp(-b * x[i]);
            } else {
                x[i] = NADBL;
            }
        }
    }

    return err;
}

private int beta_pdf_array(double a, double b, double * x, int n)
{
    int i, err = 0;

    if (na(a) || a <= 0 || na(b) || b <= 0) {
        err = GretlError.invarg;
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    } else {
        double k = exp(lngamma(a + b) - lngamma(a) - lngamma(b));

        for (i=0; i<n; i++) {
            x[i] = k * pow(x[i], a-1) * pow(1.0 - x[i], b-1);
        }
    }

    return err;
}

private double weibull_pdf(double k, double l, double x)
{
    weibull_pdf_array(k, l, &x, 1);

    return x;
}

private double exponential_pdf(double mu, double x)
{
    exponential_pdf_array(mu, &x, 1);

    return x;
}

private double beta_pdf(double a, double b, double x)
{
    beta_pdf_array(a, b, &x, 1);

    return x;
}

private int tnormal_pdf_array(double lo, double hi, double * x, int n)
{
    double P0 = na(lo) ? 0.0 : normal_cdf(lo);
    double P1 = na(hi) ? 1.0 : normal_cdf(hi);
    double den = P1 - P0;
    int i, err = 0;

    for (i=0; i<n; i++) {
        if (x[i] < lo || x[i] > hi) {
            x[i] = 0.0;
        } else {
            x[i] = normal_pdf(x[i]) / den;
        }
    }

    return err;
}

private double tnormal_pdf(double lo, double hi, double x)
{
    tnormal_pdf_array(lo, hi, &x, 1);

    return x;
}

/**
 * normal_cdf:
 * @x: double-precision value.
 *
 * Returns: the value of the standard normal CDF evaluated
 * at @x, or #NADBL on failure.
 */

double normal_cdf(double x)
{
    double y = ndtr(x);

    if (get_cephes_errno()) {
        y = NADBL;
    } else if (y == 1.0) {
        /* do we want to do this? */
        y = NORM_CDF_MAX;
    }

    return y;
}

/**
 * normal_cdf_inverse:
 * @x: double-precision value.
 *
 * Returns: the argument, y, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to y) is equal to x, or #NADBL on failure.
 */

double normal_cdf_inverse(double x)
{
    double y = ndtri(x);

    if (get_cephes_errno()) {
        y = NADBL;
    }

    return y;
}

/**
 * normal_critval:
 * @a: right-tail probability.
 *
 * Returns: the argument z such that the integral from z to
 * infinity of the standard normal density is equal
 * to the given probability @a, or #NADBL on failure.
 */

double normal_critval(double a)
{
    double z;

    if (a > 0.10) {
        z = ndtri(1.0 - a);
    } else {
        z = -ndtri(a);
    }

    if (get_cephes_errno()) {
        z = NADBL;
    }

    return z;
}

private int normal_pdf_array(double * x, int n)
{
    double s = 1.0 / SQRT_2_PI;
    int i;

    for (i=0; i<n; i++) {
        if (!na(x[i])) {
            errno = 0;
            x[i] = s * exp(-0.5 * x[i] * x[i]);
            if (errno) {
                x[i] = NADBL;
            }
        }
    }

    errno = 0;

    return 0;
}

/**
 * normal_pdf:
 * @x: double-precision value.
 *
 * Returns: the value of the standard normal PDF evaluated
 * at @x.
 */

double normal_pdf(double x)
{
    return exp(-0.5 * x * x) / SQRT_2_PI;
}

/**
 * log_normal_pdf:
 * @x: double-precision value.
 *
 * Returns: the value of the log-normal PDF evaluated
 * at @x.
 */

double log_normal_pdf(double x)
{
    return -0.5 * x * x - LN_SQRT_2_PI;
}

/**
 * tnormal_cdf:
 * @lo: left truncation point.
 * @hi: right truncation point.
 * @x: reference value.
 *
 * Calculates the value of the CDF of the truncated normal
 * distribution at @x, with truncation limits @lo and @hi

 * Returns: the calculated probability, or #NADBL on failure.
 */

double tnormal_cdf(double lo, double hi, double x)
{
    double P0 = na(lo) ? 0.0 : normal_cdf(lo);
    double P1 = na(hi) ? 1.0 : normal_cdf(hi);
    double c = P1 - P0;
    double p;

    if (x < lo) {
        p = 0;
    } else if (x > hi) {
        p = 1;
    } else {
        p = (normal_cdf(x) - P0) / c;
    }

    return p;
}

/**
 * gamma_cdf:
 * @s1: first parameter.
 * @s2: second parameter.
 * @x: reference value.
 * @control: see below.
 *
 * Calculates the value of the CDF of the gamma distribution
 * at @x.  If @control equals 1, then it is assumed that the
 * parameters @s1 and @s2 represent the shape and scale,
 * respectively, otherwise it is assumed they give mean and
 * variance.

 * Returns: the calculated probability, or #NADBL on failure.
 */

double gamma_cdf(double s1, double s2, double x, int control)
{
    double scale, shape, p;

    if (control == 1) {
        shape = s1;
        scale = s2;
    } else {
        scale = s2 / s1;
        shape = s1 / scale;
    }

    /* for the cephes functions, the parameterization is
       inverse-scale (or "rate"), shape
    */

    p = gdtr(1.0 / scale, shape, x);
    if (get_cephes_errno()) {
        p = NADBL;
    }

    return p;
}

/**
 * gamma_cdf_comp:
 * @s1: first parameter.
 * @s2: second parameter.
 * @x: reference value.
 * @control: see below.
 *
 * Calculates the complement of the CDF of the gamma distribution
 * at @x.  If @control equals 1, then it is assumed that the
 * parameters @s1 and @s2 represent the shape and scale,
 * respectively, otherwise it is assumed they give mean and
 * variance.

 * Returns: the calculated probability, or #NADBL on failure.
 */

double gamma_cdf_comp(double s1, double s2, double x, int control)
{
    double shape, scale, p;

    if (control == 1) {
        shape = s1;
        scale = s2;
    } else {
        scale = s2 / s1;    /* variance / mean */
        shape = s1 / scale; /* mean / scale */
    }

    p = gdtrc(1.0 / scale, shape, x);
    if (get_cephes_errno()) {
        p = NADBL;
    }

    return p;
}

/**
 * gamma_cdf_inverse:
 * @shape: shape.
 * @scale: scale.
 * @p: probability.
 *
 * Returns: the argument x such that the integral from zero to @x of
 * the gamma density with given scale and shape parameters is equal to
 * the given probability @p, or #NADBL on failure. Note that the
 * alternate parametrization (mean, variance) is not supported.
 */

double gamma_cdf_inverse(double shape, double scale, double p)
{
    double x = NADBL;

    if (p==0) {
        return 0;
    }

    if (shape > 0 && scale > 0 && p > 0 && p < 1) {
        x = igami(shape, 1-p) * scale;
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

private int gamma_pdf_array(double shape, double scale,
                            double * x, int n)
{
    int i, err = 0;

    errno = 0;

    if (!na(shape) && shape > 0 && !na(scale) && scale > 0) {
        double x1, x2;
        double x3 = pow(scale, shape);
        double x4 = gammafun(shape);

        if (errno || na(x4)) {
            err = GretlError.nan;
        } else {
            for (i=0; i<n; i++) {
                if (!na(x[i]) && x[i] > 0) {
                    errno = 0;
                    x1 = pow(x[i], shape - 1.0);
                    x2 = exp(-x[i]/scale);
                    x[i] = x1 * x2 / (x3 * x4);
                    if (errno) {
                        x[i] = NADBL;
                    }
                } else {
                    x[i] = NADBL;
                }
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    errno = 0;

    return err;
}

/**
 * gamma_pdf:
 * @shape: shape parameter.
 * @scale: scale parameter.
 * @x: reference value.
 *
 * Returns: the value of the gamma pdf at @x, or #NADBL on failure.
 */

double gamma_pdf(double shape, double scale, double x)
{
    gamma_pdf_array(shape, scale, &x, 1);

    return x;
}

/**
 * poisson_cdf:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability of X <= @k, for X an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

private double poisson_cdf(double lambda, int k)
{
    double x = NADBL;

    if (lambda >= 0 && k >= 0) {
        x = pdtr(k, lambda);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

/**
 * poisson_cdf_comp:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability of X > @k, for X an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

private double poisson_cdf_comp(double lambda, int k)
{
    double x = NADBL;

    if (lambda >= 0 && k >= 0) {
        x = pdtrc(k, lambda);
        if (get_cephes_errno()) {
            x = NADBL;
        }
    }

    return x;
}

private int poisson_pmf_array(double lambda, double * x, int n)
{
    double den, l0, p;
    int i, j, k;

    if (lambda <= 0) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
        return 0;
    }

    l0 = exp(-lambda);

    for (i=0; i<n; i++) {
        k = cast(int) x[i];
        den = x_factorial(cast(double) k);
        if (na(den)) {
            p = NADBL;
        } else {
            p = l0 * pow(lambda, cast(double) k) / den;
        }
        if (na(p)) {
            p = l0;
            for (j=1; j<=k; j++) {
                p *= lambda / j;
            }
        }
        x[i] = p;
    }

    return 0;
}

/**
 * poisson_pmf:
 * @lambda: mean (also variance).
 * @k: test value.
 *
 * Returns: the probability mass at @k, for an r.v. that follows
 * the Poisson distribution with parameter @lambda.
 */

double poisson_pmf(double lambda, int k)
{
    double den, l0, p;

    if (lambda <= 0 || k < 0) {
        return NADBL;
    }

    den = x_factorial(cast(double) k);
    l0 = exp(-lambda);

    if (na(den)) {
        p = NADBL;
    } else {
        p = l0 * pow(lambda, cast(double) k) / den;
    }

    if (na(p)) {
        int i;

        p = l0;
        for (i=1; i<=k; i++) {
            p *= lambda / i;
        }
    }

    return p;
}

/* The following is probably horribly inefficient --
   investigate the possibility of using Chebyshev,
   or perhaps binary search.
*/

private double poisson_critval(double mu, double a)
{
    double pk = 0.0;
    double ac = 1 - a;
    int k, k0 = 0;

    if (mu <= 0 || a <= 0 || a >= 1) {
        return NADBL;
    }

    if (mu >= 10 && a < 0.5) {
        k0 = cast(int) (mu - 1);
        pk = poisson_cdf(mu, k0++);
    }

    for (k=k0; ; k++) {
        pk = poisson_cdf(mu, k);
        if (pk >= ac) {
            break;
        }
    }

    return cast(double) k;
}

/**
 * poisson_cdf_inverse:
 * @lambda: Poisson parameter (mean = variance).
 * @p: cumulative probability.
 *
 * Returns: the Poisson variable @x such that the integral
 * from 0 to @x of the Poisson density is equal to the
 * given probability @p.
 */

/* The C source includes "poisson.h" here ("use Mike Giles's code");
   its contents are inlined below. */

/*
   This is Mike Giles's Poisson inverse code, under GPL. See
   https://people.maths.ox.ac.uk/gilesm/codes/poissinv/

   The method is explicated in Giles's paper "Algorithm 955:
   Approximation of the Inverse Poisson Cumulative Distribution
   Function", ACM Transactions on Mathematical Software, Volume 42
   Issue 1, February 2016.
*/

/*
  This double precision function computes the inverse
  of the Poisson CDF

  u = CDF value in range (0,1)
  lam = Poisson rate

  For lam < 1e15,  max |error| no more than 1
  ave |error| < 1e-16*max(4,lam) for lam < 1e9
              < 1e-6             for lam < 1e15

  For lam > 1e15, the errors will be about 1 ulp.
*/

private double poissinv_core(double U, double V, double Lam)
{
    double Xi, W, T, Del, R, R2, S, S2, Eta, B0, B1;
    double Lami = 1/Lam;
    double X = 0.0;
    int i;

    if (U == 0.0) return 0.0;
    if (V == 0.0) return cast(double) FP_INFINITE; /* sic: the integer constant, as in C */
    if (!(U > 0.0 && V > 0.0)) return cast(double) FP_NAN; /* sic: the integer constant, as in C */

    if (Lam > 4.0) {
        W = normal_cdf_inverse(fmin(U, V));
        if (U > V) W = -W;

        /* use polynomial approximations in central region */
        if (fabs(W) < 3.0) {
            double Lam_root = sqrt(Lam);

            S = Lam_root*W + (1.0/3.0 + (1.0/6.0)*W*W)*(1.0 - W/(12.0*Lam_root));

            Del = (1.0 /160.0);
            Del = (1.0 / 80.0) + Del*(W*W);
            Del = (1.0 / 40.0) + Del*(W*W);
            Del = Del * Lami;

            S = Lam + (S + Del);
        } else {
            /* otherwise use Newton iteration */
            S = W / sqrt(Lam);
            R = 1.0 + S;
            if (R < 0.1) R = 0.1;

            do {
                T = log(R);
                R2 = R;
                S2 = sqrt(2.0*((1.0-R) + R*T));
                if (R < 1.0) S2 = -S2;
                R = R2 - (S2-S)*S2/T;
                if (R < 0.1*R2) R = 0.1*R2;
            } while (fabs(R-R2) > 1e-8);

            T = log(R);
            S = Lam*R + log(sqrt(2.0*R*((1.0-R)+R*T))/fabs(R-1.0)) / T;
            S = S - 0.0218/(S+0.065*Lam);
            Del = 0.01/S;
            S = S + Del;
        }

        /* if x > 10, round down to nearest integer, and check accuracy */
        X = floor(S);

        if (S > 10.0 && S < X + 2.0*Del) {
            /* correction procedure based on Temme approximation */
            if (X > 0.5*Lam && X < 2.0*Lam) {
                Xi = 1.0 / X;
                Eta = X * Lami;
                Eta = sqrt(2.0*(1.0-Eta+Eta*log(Eta))/Eta);
                if (X > Lam) Eta = -Eta;

                B1 =  8.0995211567045583e-16;              S = B1;
                B0 = -1.9752288294349411e-15;              S = B0 + S*Eta;
                B1 = -5.1391118342426808e-16 + 25.0*B1*Xi; S = B1 + S*Eta;
                B0 =  2.8534893807047458e-14 + 24.0*B0*Xi; S = B0 + S*Eta;
                B1 = -1.3923887224181616e-13 + 23.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.3717632624009806e-13 + 22.0*B0*Xi; S = B0 + S*Eta;
                B1 =  1.1004392031956284e-13 + 21.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.0276692801141763e-12 + 20.0*B0*Xi; S = B0 + S*Eta;
                B1 =  2.4361948020667402e-11 + 19.0*B1*Xi; S = B1 + S*Eta;
                B0 = -5.8307721325504166e-11 + 18.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.5514193994946487e-11 + 17.0*B1*Xi; S = B1 + S*Eta;
                B0 =  9.1476995822367933e-10 + 16.0*B0*Xi; S = B0 + S*Eta;
                B1 = -4.3820360184533521e-09 + 15.0*B1*Xi; S = B1 + S*Eta;
                B0 =  1.0261809784240299e-08 + 14.0*B0*Xi; S = B0 + S*Eta;
                B1 =  6.7078535434015332e-09 + 13.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.7665952736826086e-07 + 12.0*B0*Xi; S = B0 + S*Eta;
                B1 =  8.2967113409530833e-07 + 11.0*B1*Xi; S = B1 + S*Eta;
                B0 = -1.8540622107151585e-06 + 10.0*B0*Xi; S = B0 + S*Eta;
                B1 = -2.1854485106799979e-06 +  9.0*B1*Xi; S = B1 + S*Eta;
                B0 =  3.9192631785224383e-05 +  8.0*B0*Xi; S = B0 + S*Eta;
                B1 = -0.00017875514403292177 +  7.0*B1*Xi; S = B1 + S*Eta;
                B0 =  0.00035273368606701921 +  6.0*B0*Xi; S = B0 + S*Eta;
                B1 =   0.0011574074074074078 +  5.0*B1*Xi; S = B1 + S*Eta;
                B0 =   -0.014814814814814815 +  4.0*B0*Xi; S = B0 + S*Eta;
                B1 =    0.083333333333333329 +  3.0*B1*Xi; S = B1 + S*Eta;
                B0 =    -0.33333333333333331 +  2.0*B0*Xi; S = B0 + S*Eta;
                S  = S / (1.0 + B1*Xi);

                S = S*exp(-0.5*X*Eta*Eta)/sqrt(2.0*3.141592653589793*X);
                if (X < Lam) {
                    S += 0.5 * erfc(Eta * sqrt(0.5*X));
                    if (S > U) X -= 1.0;
                } else {
                    S -= 0.5 * erfc(-Eta * sqrt(0.5*X));
                    if (S > -V) X -= 1.0;
                }
            } else {
                /* sum downwards or upwards */
                Xi = 1.0 / X;
                S = -(691.0/360360.0);
                S =  (1.0/1188.0) + S*Xi*Xi;
                S = -(1.0/1680.0) + S*Xi*Xi;
                S =  (1.0/1260.0) + S*Xi*Xi;
                S = -(1.0/360.0)  + S*Xi*Xi;
                S =  (1.0/12.0)   + S*Xi*Xi;
                S = S*Xi;
                S = (X - Lam) - X * log(X*Lami) - S;

                if (X < Lam) {
                    T = exp(-0.5*S);
                    S = 1.0 - T*(U*T) * sqrt(2.0*3.141592653589793*Xi) * Lam;
                    T = 1.0;
                    Xi = X;
                    for (i=1; i<50; i++) {
                        Xi -= 1.0;
                        T *= Xi*Lami;
                        S += T;
                    }
                    if (S > 0.0) X -= 1.0;
                } else {
                    T = exp(-0.5*S);
                    S = 1.0 - T*(V*T) * sqrt(2.0*3.141592653589793*X);
                    Xi = X;
                    for (i=0; i<50; i++) {
                        Xi += 1.0;
                        S = S*Xi*Lami + 1.0;
                    }
                    if (S < 0.0) {
                        X -= 1.0;
                    }
                }
            }
        }
    }

    /* bottom-up summation */
    if (X < 10.0) {
        X = 0.0;
        T = exp(0.5*Lam);
        Del = 0.0;
        if (U > 0.5) {
            Del = T*(1e-13*T);
        }
        S = 1.0 - T*(U*T) + Del;

        while (S < 0.0) {
            X += 1.0;
            T = X*Lami;
            Del = T*Del;
            S = T*S + 1.0;
        }

        /* top-down summation if needed */
        if (S < 2.0*Del) {
            Del = 1e13*Del;
            T = 1e17*Del;
            Del = V*Del;

            while (Del < T) {
                X += 1.0;
                Del *= X*Lami;
            }

            S = Del;
            T = 1.0;
            while (S > 0.0) {
                T *= X*Lami;
                S -= T;
                X -= 1.0;
            }
        }
    }

    return X;
}

private double poissinv(double U, double Lam)
{
    return poissinv_core(U, 1.0 - U, Lam);
}

private double poisson_cdf_inverse(double lambda, double p)
{
    return poissinv(p, lambda);
}

private double weibull_critval(double shape, double scale,
                               double rtail)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0 && rtail > 0 && rtail < 1) {
        ret = scale * pow(-log(rtail), 1.0 / shape);
    }

    return ret;
}

/**
 * weibull_cdf:
 * @shape: shape parameter > 0.
 * @scale: scale parameter > 0.
 * @x: test value.
 *
 * Returns: the probability of X <= @x, for X an r.v. that follows
 * the Weibull distribution with parameters @shape and @scale.
 */

double weibull_cdf(double shape, double scale, double x)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0 && !na(x)) {
        if (x == 0.0) {
            ret = 0.0;
        } else if (x > 0.0) {
            ret = 1.0 - exp(-pow(x/scale, shape));
        }
    }

    return ret;
}

/**
 * exponential_cdf:
 * @mu: scale parameter > 0.
 * @x: test value.
 *
 * Returns: the probability of X <= @x, for X an r.v. that follows
 * the exponential distribution with scale parameter @mu.
 */

double exponential_cdf(double mu, double x)
{
    double ret = NADBL;

    if (mu > 0 && !na(x)) {
        if (x < 0.0) {
            ret = 0.0;
        } else {
            ret = 1.0 - exp(-x / mu);
        }
    }

    return ret;
}

private double weibull_cdf_comp(double shape, double scale, double x)
{
    double ret = NADBL;

    if (shape > 0 && scale > 0 && !na(x)) {
        if (x == 0.0) {
            ret = 1.0;
        } else if (x > 0.0) {
            ret = exp(-pow(x/scale, shape));
        }
    }

    return ret;
}

/**
 * GED_pdf:
 * @nu: shape parameter.
 * @x: reference value.
 *
 * Returns: the density function of the Generalized Error distribution
 * with shape parameter @nu at @x, or #NADBL on failure.
 */

double GED_pdf(double nu, double x)
{
    if (nu > 0) {
        double lg1 = lngamma(1/nu);
        double lg3 = lngamma(3/nu);
        double lC  = 0.5*(lg3 - 3*lg1);
        double k   = pow(0.5, 1/nu) * exp(0.5*(lg1 - lg3));
        double znu = pow(fabs(x/k), nu);

        return (0.5 * nu) * exp(lC - 0.5 * znu);
    } else {
        return NADBL;
    }
}

private int GED_pdf_array(double nu, double * x, int n)
{
    int i, err = 0;

    if (nu > 0) {
        double lg1 = lngamma(1/nu);
        double lg3 = lngamma(3/nu);
        double lC  = 0.5*(lg3 - 3*lg1);
        double k   = pow(0.5, 1/nu) * exp(0.5*(lg1 - lg3));
        double znu;

        for (i=0; i<n; i++) {
            if (!na(x[i])) {
                znu = pow(fabs(x[i]/k), nu);
                x[i] = (0.5 * nu) * exp(lC - 0.5 * znu);
            } else {
                x[i] = NADBL;
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    return err;
}

/**
 * GED_cdf:
 * @nu: shape parameter.
 * @x: reference value.
 *
 * Calculates the value of the CDF of the Generalized Error
 * distribution with parameter @nu at @x. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double GED_cdf(double nu, double x)
{
    if (nu > 0) {
        int sgn    = (x > 0)? 1 : -1;
        double p   = 1/nu;
        double lg1 = lngamma(p);
        double lg3 = lngamma(3*p);
        double k   = pow(0.5, p) * exp(0.5*(lg1 - lg3));
        double znu = pow(fabs(x/k), nu);
        double P   = gamma_cdf(p, 2, znu, 1);

        P = 0.5 * (1 + sgn*P);
        return P;
    } else {
        return NADBL;
    }
}

/**
 * GED_cdf_comp:
 * @nu: shape parameter.
 * @x: reference value.
 *
 * Calculates the complement of the CDF of the Generalized Error
 * distribution with parameter @nu at @x. We exploit the fact that
 * if x ~ GED(n), then |x/k|^n is a Gamma rv.
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double GED_cdf_comp(double nu, double x)
{
    if (nu > 0) {
        int sgn    = (x > 0)? 1 : -1;
        double p   = 1/nu;
        double lg1 = lngamma(p);
        double lg3 = lngamma(3*p);
        double k   = pow(0.5, p) * exp(0.5*(lg1 - lg3));
        double znu = pow(fabs(x/k), nu);
        double P   = gamma_cdf_comp(p, 2, znu, 1);

        P = (sgn == 1) ? 0.5 * P : 1 - 0.5 * P;
        return P;
    } else {
        return NADBL;
    }
}

/**
 * GED_cdf_inverse:
 * @nu: shape parameter.
 * @a: probability.
 *
 * Returns: the argument x such that the integral from minus infinity
 * to @x of the standardized GED density with shape parameter @nu is
 * equal to the given probability @a, or #NADBL on failure. We exploit
 * the well-known relationship between the standardized GED and the
 * Gamma variates.
 */

double GED_cdf_inverse(double nu, double a)
{
    if (nu > 0 && a < 1 && a > 0) {
        double a2, p, lg1, lg3, sd, x;
        int sgn;

        if (a > 0.5) {
            a2 = 2*a - 1;
            sgn = 1;
        } else {
            a2 = 1 - 2*a;
            sgn = -1;
        }

        p   = 1/nu;
        lg1 = lngamma(p);
        lg3 = lngamma(3*p);
        sd  = pow(2.0, p) * exp(0.5*(lg3 - lg1));
        x   = gamma_cdf_inverse(p, 2.0, a2);

        return sgn * pow(x, p) / sd;
    } else {
        return NADBL;
    }
}

/**
 * laplace_pdf:
 * @mu: mean.
 * @b: scale (greater than 0).
 * @x: reference value.
 *
 * Returns: the density function of the Laplace distribution
 * with mean @mu and scale @b evaluated at @x, or #NADBL on failure.
 */

double laplace_pdf(double mu, double b, double x)
{
    if (b > 0) {
        return exp(-fabs(x - mu)/b) / (2*b);
    } else {
        return NADBL;
    }
}

private int laplace_pdf_array(double mu, double b,
                              double * x, int n)
{
    int i, err = 0;

    if (b > 0) {
        for (i=0; i<n; i++) {
            if (!na(x[i])) {
                x[i] = exp(-fabs(x[i] - mu)/b) / (2*b);
            } else {
                x[i] = NADBL;
            }
        }
    } else {
        err = GretlError.data;
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    return err;
}

/**
 * laplace_cdf:
 * @mu: mean.
 * @b: scale (greater than 0).
 * @x: reference value.
 *
 * Returns: the CDF of the Laplace distribution
 * with mean @mu and scale @b evaluated at @x, or #NADBL on failure.
 */

double laplace_cdf(double mu, double b, double x)
{
    if (b > 0) {
        if (x < mu) {
            return 0.5 * exp((x-mu)/b);
        } else {
            return 1 - 0.5 * exp(-(x-mu)/b);
        }
    } else {
        return NADBL;
    }
}

/**
 * laplace_cdf_comp:
 * @mu: mean.
 * @b: scale (greater than 0).
 * @x: reference value.
 *
 * Returns: the complement of the CDF of the Laplace distribution
 * with mean @mu and scale @b evaluated at @x, or #NADBL on failure.
 */

double laplace_cdf_comp(double mu, double b, double x)
{
    if (b > 0) {
        if (x < mu) {
            return 1 - 0.5 * exp((x-mu)/b);
        } else {
            return 0.5 * exp(-(x-mu)/b);
        }
    } else {
        return NADBL;
    }
}

/**
 * laplace_cdf_inverse:
 * @mu: mean.
 * @b: scale (greater than 0).
 * @a: probability.
 *
 * Returns: the argument x such that the integral from minus infinity
 * to @x of the Laplace density with mean @mu and scale @b is
 * equal to the given probability @a, or #NADBL on failure.
 */

double laplace_cdf_inverse(double mu, double b, double a)
{
    if (b <= 0 || a < 0 || a > 1) {
        return NADBL;
    } else if (a == 0.0) {
        return -1.0 / 0.0;
    } else if (a == 1.0) {
        return 1.0 / 0.0;
    } else {
        int sgn = a - 0.5 < 0 ? -1 : 1;

        return mu - b*sgn * log(1.0 - 2*fabs(a - 0.5));
    }
}

/**
 * johansen_trace_pval:
 * @N: the number of potentially cointegrated variables
 * minus the cointegration rank under H0.
 * @det: index of deterministic setup of model (0 to 4).
 * @T: the sample size, or 0 for an asymptotic result.
 * @tr: the trace statistic.
 *
 * Returns: the p-value of the trace statistic, computed
 * via Doornik's gamma approximation, or #NADBL on failure.
 */

double johansen_trace_pval(int N, int det, int T, double tr)
{
    TracePvalueFn pvfunc;
    double pv = NADBL;

    pvfunc = cast(TracePvalueFn) get_plugin_function("trace_pvalue");

    if (pvfunc != null) {
        pv = (*pvfunc)(tr, N, det, T);
    }

    return pv;
}

/* below: non-central distributions: chi-square, F and
   Student's t
*/

private bool qsmall(double sum, double x)
{
    return (sum < 1.0e-30 || (x) < 1.0e-8*sum);
}

/**
 * nc_chisq_cdf:
 * @df: degrees of freedom.
 * @delta: noncentrality parameter.
 * @x: reference value.
 *
 * Calculates the value at @x of the CDF of the noncentral chi^2
 * distribution with @df dof and noncentrality parameter equal to
 * @delta.
 *
 * This is a version of cumchn() from dcdflib, de-spaghettized by
 * Jack Lucchetti (2015-06-21). The original algorithm uses formula
 * 26.4.25 from Abramowitz and Stegun, Handbook of Mathematical
 * Functions, US NBS (1966).
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double nc_chisq_cdf(double df, double delta, double x)
{
    double adj, centaj, centwt, chid2, dfd2, lcntaj, lcntwt;
    double lfact, pcent, pterm, sum, sumadj, term, wt, xnonc;
    double T1, T2, T3;
    int i, icent, iterb, iterf;
    int itermax = 1000;

    if (x < 0.0) {
        return 1.0;
    }

    if (df <= 0.0) {
        return NADBL;
    }

    if (delta <= 1.0e-10) {
        /* When non-centrality parameter is (essentially) zero,
          use ordinary chi-square distribution
        */
        return chisq_cdf(df, x);
    }

    xnonc = delta / 2.0;

    /* The following code calculates the weight, chi-square, and
       adjustment term for the central term in the infinite series.  The
       central term is the one in which the poisson weight is greatest.
       The adjustment term is the amount that must be subtracted from
       the chi-square to move up two degrees of freedom.
    */
    icent = (xnonc < 1.0) ? 1 : cast(int) trunc(xnonc);
    chid2 = x / 2.0;

    /* Calculate central weight term */
    T1 = cast(double) (icent + 1);
    lfact = lngamma(T1);
    lcntwt = -xnonc + cast(double) icent*log(xnonc) - lfact;
    centwt = exp(lcntwt);

    /* Calculate central chi-square */
    T2 = df + 2.0 * cast(double) icent;
    pcent = chisq_cdf(T2, x);

    /* Calculate central adjustment term */
    dfd2 = df / 2.0 + icent;
    T3 = dfd2 + 1.0;
    lfact = lngamma(T3);
    lcntaj = dfd2 * log(chid2) - chid2 - lfact;
    centaj = exp(lcntaj);
    sum = centwt * pcent;

    /* Sum backwards from the central term towards zero.
       Quit whenever either
       (1) the zero term is reached, or
       (2) the term gets small relative to the sum, or
       (3) More than NTIRED terms are totaled.
    */

    iterb = 0;
    sumadj = 0.0;
    adj = centaj;
    wt = centwt;
    i = icent;

    do {
        dfd2 = df/2.0 + i;
        /* Adjust chi-square for two fewer degrees of freedom.
           The adjusted value ends up in PTERM.
        */
        adj = adj * dfd2/chid2;
        sumadj += adj;
        pterm = pcent + sumadj;
        /* Adjust poisson weight for J decreased by one */
        wt *= (cast(double) i/xnonc);
        term = wt * pterm;
        sum += term;
        i--;
        iterb++;
    } while (iterb <= itermax && !qsmall(sum, term) && i > 0);

    /* Sum forward from the central term towards infinity.
       Quit when either
       (1) the term gets small relative to the sum, or
       (2) More than NTIRED terms are totaled.
    */
    iterf = 0;
    sumadj = adj = centaj;
    wt = centwt;
    i = icent;

    do {
        /* Update weights for next higher J */
        wt *= (xnonc/cast(double)(i+1));
        /* Calculate PTERM and add term to sum */
        pterm = pcent - sumadj;
        term = wt * pterm;
        sum += term;
        /* Update DF adjustment term for next iteration */
        i++;
        dfd2 = df/2.0 + i;
        adj = adj * chid2/dfd2;
        sumadj += adj;
        iterf++;
    } while (iterf <= itermax && !qsmall(sum, term));

    return sum;
}

/**
 * nc_chisq_pdf_array:
 * @p: degrees of freedom.
 * @c: noncentrality parameter.
 * @x: array of arguments (overwritten on exit).
 * @n: no. of elements in x.
 *
 * Calculates the value at @x of the CDF of the noncentral chi^2
 * distribution with @p dof and noncentrality parameter equal
 * to @c.
 *
 * Returns: an error code, as appropriate.
 */

private int nc_chisq_pdf_array(double p, double c, double * x, int n)
{
    int i, err = 0;
    double k, a, b;

    if (fabs(c) < 1.0e-10) {
        return chisq_pdf_array(cast(int) floor(p), x, n);
    }

    if (p <= 0 || c < 0) {
        return GretlError.data;
    }

    k = p/2.0 - 1;

    for (i=0; i<n; i++) {
        if (na(x[i]) || x[i] < 0) {
            x[i] = NADBL;
        } else {
            a = exp(-0.5*(x[i]+c) + k/2.0 * log(x[i]/c)) / 2.0;
            b = gretl_bessel('I', k, sqrt(x[i]*c), &err);
            if (err) {
                break;
            }
            x[i] = a*b;
        }
    }

    return err;
}

double nc_chisq_pdf(double p, double c, double x)
{
    nc_chisq_pdf_array(p, c, &x, 1);
    return x;
}

/**
 * nc_chisq_cdf_inverse:
 * @p: degrees of freedom
 * @c: noncentrality parameter.
 * @q: probability.
 *
 * Calculates the @q-th quantile of the noncentral chi^2
 * distribution with @n1, @n2 dof and noncentrality parameter equal to
 * @c via a rough and not particularly clever root-finding
 * algorithm. Maybe this can be more efficient by using logs. Some
 * experimentation needed.
 *
 * Returns: the calculated quantile, or #NADBL on failure.
 */

private double nc_chisq_cdf_inverse(double p, double c, double q)
{
    double x, d0, d1;
    int iter, subiter, retry;
    double F, f, dir;

    if (p < 0 || c < 0 || q <= 0 || q >= 1) {
        return NADBL;
    }

    if (fabs(c) < 1.0e-10) {
        /* don't bother for infinitesimal c */
        return chisq_cdf_inverse(cast(int) p, q);
    }

    /* start from the mean (safe bet) */
    x = p + c;
    d0 = 1.0e7;
    iter = 0;

    while (fabs(d0) > 1.0e-10 && iter < 1000) {
        F = nc_chisq_cdf(p, c, x);
        f = nc_chisq_pdf(p, c, x);
        d0 = F - q;
        dir = d0/f;
        d1 = 1.0e7;
        retry = 1;
        subiter = 0;

        while (retry && subiter < 100) {
            if ((x-dir) > 0) {
                d1 = F - nc_chisq_cdf(p, c, x - dir);
            }
            dir /= 2.0;
            retry = (x-dir) < 0 || fabs(d1) > fabs(d0);
            subiter++;
        }

        if (subiter >= 100) {
            x = NADBL;
            break;
        } else {
            x -= dir*2;
            d0 = d1;
            iter++;
        }
    }

    if (iter >= 1000) {
        x = NADBL;
    }

    return x;
}

/**
 * nc_snedecor_cdf:
 * @dfn: degrees of freedom (num).
 * @dfd: degrees of freedom (den).
 * @delta: noncentrality parameter.
 * @x: reference value.
 *
 * Calculates the value at @x of the CDF of the noncentral F
 * distribution with @dfn, @dfd dof and noncentrality parameter equal
 * to @delta.
 *
 * This is a version of cumfnc() from dcdflib, de-spaghettized by
 * Jack Lucchetti (2015-06-21). The original algorithm uses formula
 * 26.6.18 from Abramowitz and Stegun, Handbook of Mathematical
 * Functions, US NBS (1966).
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double nc_snedecor_cdf(double dfn, double dfd, double delta, double x)
{
    double dsum, prod, xx, yy, adn, aup, b, betdn, betup;
    double centwt, dnterm, sum, upterm, xmult, xnonc;
    double T1, T2, T3, T4, T5, T6;
    int i, icent;

    if (x < 0.0) {
        return 1.0;
    }

    if (dfn <= 0.0 || dfd <= 0.0) {
        return NADBL;
    }

    if (delta <= 1.0e-10) {
        /* When non-centrality parameter is (essentially) zero, use
           ordinary F distribution
        */
        return snedecor_cdf(dfn, dfd, x);
    } else {
        xnonc = delta / 2.0;
    }

    /* Calculate the central term of the poisson weighting factor */
    icent = (xnonc < 1.0) ? 1 : cast(int) trunc(xnonc);

    /* Compute central weight term */
    T1 = cast(double) (icent + 1);
    centwt = exp(-xnonc + cast(double) icent * log(xnonc) - lngamma(T1));

    /* Compute central incomplete beta term. Ensure that minimum of arg
       to beta and 1 - arg is computed accurately.
    */
    prod = dfn * x;
    dsum = dfd + prod;
    yy = dfd / dsum;
    if (yy > 0.5) {
        xx = prod / dsum;
        yy = 1.0 - xx;
    } else {
        xx = 1.0 - yy;
    }

    T2 = dfn / 2.0 + cast(double) icent;
    T3 = dfd / 2.0;
    betdn = incbet(T2, T3, xx);
    adn = dfn / 2.0 + cast(double) icent;
    aup = adn;
    b = dfd / 2.0;
    betup = betdn;
    sum = centwt * betdn;

    /* Sum terms backward from icent until convergence or all done */

    xmult = centwt;
    i = icent;
    T4 = adn + b;
    T5 = adn + 1.0;
    dnterm = exp(lngamma(T4)-lngamma(T5)-lngamma(b) + adn*log(xx) + b*log(yy));

    while (!qsmall(sum, xmult*betdn) && i > 0) {
        xmult *= cast(double) i /xnonc;
        i--;
        adn -= 1.0;
        dnterm *= (adn + 1.0)/((adn+b)*xx);
        betdn += dnterm;
        sum += xmult * betdn;
    }

    i = icent+1;

    /* Sum forwards until convergence */

    xmult = centwt;
    if (aup-1.0+b == 0) {
        upterm = exp(-lngamma(aup) - lngamma(b) +
                     (aup-1.0)*log(xx) + b*log(yy));
    } else {
        T6 = aup - 1.0 + b;
        upterm = exp(lngamma(T6) - lngamma(aup) - lngamma(b) +
                     (aup-1.0)*log(xx) + b*log(yy));
    }

    while (!qsmall(sum, xmult*betup)) {
        xmult *= (xnonc/cast(double) i);
        i++;
        aup += 1.0;
        upterm *= (aup + b - 2.0) * xx/(aup-1.0);
        betup -= upterm;
        sum += (xmult*betup);
    }

    return sum;
}

/**
 * nc_snedecor_pdf_array:
 * @dfn: degrees of freedom (num).
 * @dfd: degrees of freedom (den).
 * @c: noncentrality parameter.
 * @x: array of arguments (overwritten on exit).
 * @n: no. of elements in x.
 *
 * Calculates the value at @x of the CDF of the noncentral F
 * distribution with @dfn, @dfd dof and noncentrality parameter equal
 * to @c.
 *
 * Source: S. Kay, Fundamentals of Statistical Signal Processing:
 * Detection Theory, (New Jersey: Prentice Hall, 1998),
 * <TeX>
 * p(x) = \sum\limits_{k=0}^\infty
 *   \frac{e^{-c/2}(c/2)^k}{k!}  % Poisson weights
 *   \frac{1}{B\left(\frac{\nu_2}{2},\frac{\nu_1}{2}+k\right)} % Beta function
 *   \left(\frac{\nu_1}{\nu_2}\right)^{\frac{\nu_1}{2}+k}
 *   \left(\frac{\nu_2}{\nu_2+\nu_1f}\right)^{\frac{\nu_1+\nu_2}{2}+k}
 *   x^{\nu_1/2-1+k}
 * </TeX>
 * coded in C by Jack Lucchetti (2015-06-27).
 *
 * Returns: an error code, as appropriate.
 */

private int nc_snedecor_pdf_array(double dfn, double dfd, double c,
                                  double * x, int n)
{
    double ch, k1, k2, k;
    double a, b, l, pw, beta;
    double pwi, betai;
    double errtol = 1.0e-16;
    double maxit = 256;
    double * vx, vz;
    int i, t, start, iter;
    int err = 0;

    if (dfd <= 0.0 || dfn <= 0.0 || c < 0.0) {
        return GretlError.data;
    }

    if (fabs(c) <= 1.0e-10) {
        /* When non-centrality parameter is (essentially) zero,
           use ordinary F distribution
        */
        return snedecor_pdf_array(dfn, dfd, x, n);
    }

    vx  = cast(double *) malloc(n * double.sizeof);
    vz  = cast(double *) malloc(n * double.sizeof);

    if (vx == null || vz == null) {
        err = GretlError.alloc;
        goto bailout;
    }

    /* fill up auxiliary vectors */
    ch = c/2.0;
    k1 = dfn/2;
    k2 = (dfn + dfd)/2;
    k = log(dfn) - log(dfd);

    for(t=0; t<n; t++) {
        if (na(x[t]) || x[t] < 0) {
            vx[t] = vz[t] = NADBL;
        } else {
            vx[t] = log(x[t]);
            vz[t] = log(dfd) - log(dfd + dfn * x[t]);
        }
    }

    /* start from central Poisson weight */

    start = cast(int) floor(ch);

    /* Poisson weight */
    pw = exp(-ch - lngamma(start+1) + start * log(ch));
    /* Beta */
    beta = exp(lngamma(k1+start) + lngamma(dfd/2.0) - lngamma(start + k2));
    a = pw / beta;
    b = (k1 + start) * k;

    for (t=0; t<n; t++) {
        if (na(x[t]) || x[t] < 0) {
            x[t] = NADBL;
        } else {
            l = b + (start + k1 - 1) * vx[t] + (start + k2) * vz[t];
            x[t] = a * exp(l);
        }
    }

    /* Go back from start to 0 */

    pwi = pw;
    betai = beta;
    iter = 0;

    for (i = start-1; i>=0 && pwi>errtol && iter < maxit; i--) {
        iter++;
        pwi *= (i + 1.0)/ch;
        betai *= (k2 + i)/(k1 + i);
        a = pwi / betai;
        b = (k1 + i) * k;
        for (t=0; t<n; t++) {
            if (!na(x[t])) {
                l = b + (i + k1 - 1) * vx[t] + (i + k2) * vz[t];
                x[t] += a * exp(l);
            }
        }
    }

    /* Go from start all the way up as necessary */

    iter = 0;
    pwi = pw;
    betai = beta;

    for (i = start+1; pwi>errtol && iter<maxit; i++) {
        iter++;
        pwi *= ch/i;
        betai *= (k1 + i - 1.0)/(k2 + i - 1.0);
        a = pwi / betai;
        b = (k1 + i) * k;
        for (t=0; t<n; t++) {
            if (!na(x[t])) {
                l = b + (i + k1 - 1) * vx[t] + (i + k2) * vz[t];
                x[t] += a * exp(l);
            }
        }
    }

 bailout:

    free(vx);
    free(vz);

    return err;
}

double nc_snedecor_pdf(double dfn, double dfd, double c, double x)
{
    nc_snedecor_pdf_array(dfn, dfd, c, &x, 1);

    return x;
}

/**
 * nc_snedecor_cdf_inverse:
 * @n1: degrees of freedom (numerator).
 * @n2: degrees of freedom (denominator).
 * @c: noncentrality parameter.
 * @q: probability.
 *
 * Calculates the @q-th quantile of the noncentral F distribution with
 * @n1, @n2 dof and noncentrality parameter equal to @c via a rough
 * and not particularly clever root-finding algorithm. Maybe this can
 * be more efficient by using logs. Some experimentation needed.
 *
 * Returns: the calculated quantile, or #NADBL on failure.
 */

private double nc_snedecor_cdf_inverse(double n1, double n2,
                                       double c, double q)
{
    double x, d0, d1;
    int iter, subiter;
    double F, f, dir;

    if (n2 < 1 || n1 < 1 || c < 0 || q <= 0 || q >= 1) {
        return NADBL;
    }

    x = 0.5;
    d0 = 1.0e7;
    iter = 0;

    while (fabs(d0) > 1.0e-10 && iter < 1000) {
        F = nc_snedecor_cdf(n1, n2, c, x);
        f = nc_snedecor_pdf(n1, n2, c, x);
        d0 = F - q;
        dir = d0/f;
        d1 = 1.0e7;
        subiter = 0;

        while (fabs(d1) > fabs(d0) && subiter < 100) {
            d1 = F - nc_snedecor_cdf(n1, n2, c, x - dir);
            dir /= 2.0;
            subiter++;
        }

        if (subiter >= 100) {
            x = NADBL;
            break;
        } else {
            x -= dir*2;
            d0 = d1;
            iter++;
        }
    }

    if (iter >= 1000) {
        x = NADBL;
    }

    return x;
}

private enum ISQRT_2 = .707106781186547524401;

/**
 * nc_student_cdf:
 * @df: degrees of freedom.
 * @delta: noncentrality parameter.
 * @x: reference value.
 *
 * Calculates the value at @x of the CDF of the noncentral Student t
 * distribution with @df dof and noncentrality parameter equal to
 * @delta. The algorithm is by Benson-Krishnamoorthy (2003) CSDA 43,
 * with minimal changes.
 *
 * Returns: the calculated probability, or #NADBL on failure.
 */

double nc_student_cdf(double df, double delta, double x)
{
    double errtol = 1.0e-16;
    double maxit = 512;
    double ax, y, del, dels, k, a, b, c;
    double pkf, pkb, qkf, qkb, pbetaf, pbetab, qbetaf, qbetab;
    double pgamf, pgamb, qgamf, qgamb, tmp, ret;
    double rempois, sum, ptermf, qtermf, ptermb, qtermb;
    int i;

    if (df <= 0.0) {
        return NADBL;
    }

    if (fabs(delta) <= 1.0e-10) {
        /* When non-centrality parameter is (essentially) zero,
           use ordinary t distribution
        */
        return student_cdf(df, x);
    }

    ax = fabs(x);

    del = x > 0 ? delta : -delta;
    ret = normal_cdf(-del);

    if (ax < 1.0e-12) {
        return 1.0 - ret;
    }

    /* settings */

    y = ax*ax / (df + ax*ax);
    dels = del * del / 2.0;
    k = cast(int) floor(dels);
    a = k + 0.5;
    c = k + 1;
    b = df * 0.5;

    /* Initialization to compute the P_k and Q_k terms
       and the respective incomplete beta functions
    */

    tmp = -dels + k * log(dels);
    pkf = pkb = exp(tmp - lngamma(k + 1));
    qkf = qkb = exp(tmp - lngamma(k + 1.5));
    pbetaf = pbetab = incbet(a, b, y);
    qbetaf = qbetab = incbet(c, b, y);

    /* Initialization to compute the incomplete beta functions
       associated with the P_i and the Q_i recursively:
    */

    tmp = b * log(1-y) - lngamma(b);
    pgamf = exp(lngamma(a+b-1) - lngamma(a) + (a-1) * log(y) + tmp);
    pgamb = pgamf * y * (a + b - 1)/a;

    qgamf = exp(lngamma(c+b-1) - lngamma(c) + (c-1) * log(y) + tmp);
    qgamb = qgamf * y * (c + b - 1)/c;

    /* Compute the remainder of the Poisson weights */

    rempois = 1.0 - pkf;
    sum = pkf * pbetaf + del * qkf * qbetaf * ISQRT_2;

    for (i = 1; i<=k && rempois>errtol; i++) {
        /* first block, backwards */
        pgamb *= (a-i+1)/(y * (a+b-i));
        pbetab += pgamb;
        pkb *= (k-i+1)/dels;
        ptermb = pkb * pbetab;

        /* second block, backwards */
        qgamb *= (c-i+1)/(y * (c+b-i));
        qbetab += qgamb;
        qkb *= (k-i+1.5)/dels;
        qtermb = qkb * qbetab;

        /* accumulate */
        sum += ptermb + del * qtermb * ISQRT_2;
        rempois -= pkb;
    }

    for (i = 1; i<maxit && rempois > errtol; i++) {
        /* first block, forwards */
        pgamf *= y * (a+b-2+i)/(a+i-1);
        pbetaf -= pgamf;
        pkf *= dels/(k+i);
        ptermf = pkf * pbetaf;

        /* second block, forwards */
        qgamf *= y * (c+b-2+i)/(c+i-1);
        qbetaf -= qgamf;
        qkf *= dels/(k+i+0.5);
        qtermf = qkf * qbetaf;

        /* accumulate */
        sum += ptermf + del * qtermf * ISQRT_2;
        rempois -= pkf;
    }

    ret += sum/2.0;

    return x < 0 ? (1.0 - ret) : ret;
}

/**
 * nc_student_pdf:
 * @df: degrees of freedom.
 * @delta: noncentrality parameter.
 * @x: reference value.
 *
 * Calculates the value at @x of the PDF of the noncentral Student t
 * distribution with @df dof and noncentrality parameter equal to
 * @delta. The algorithm is from Wikipedia, apparently used in R too.
 *
 * Returns: the calculated density, or #NADBL on failure.
 */

double nc_student_pdf(double df, double delta, double x)
{
    double ret, tmp;

    if (df <= 0.0) {
        return NADBL;
    }

    if (fabs(delta) <= 1.0e-10) {
        /* When non-centrality parameter is (essentially) zero,
           use ordinary t distribution
        */
        return student_pdf(df, x);
    }

    if (fabs(x) < 1.0e-12) {
        tmp = lngamma((df+1)/2) - lngamma(df/2);
        ret = exp(tmp - 0.5 * delta*delta) / (sqrt(M_PI * df));
    } else {
        tmp = nc_student_cdf(df+2, delta, x * sqrt(1 + 2.0/df)) -
            nc_student_cdf(df, delta, x);
        ret = tmp * (df / x);
    }

    return ret;
}


private int nc_student_pdf_array(double df, double delta,
                                 double * x, int n)
{
    int i, err = 0;

    if (df > 0) {
        for (i=0; i<n; i++) {
            if (!na(x[i])) {
                x[i] = nc_student_pdf(df, delta, x[i]);
            } else {
                x[i] = NADBL;
            }
        }
    } else {
        err = GretlError.data;
    }

    if (err) {
        for (i=0; i<n; i++) {
            x[i] = NADBL;
        }
    }

    return err;
}

/**
 * nc_student_cdf_inverse:
 * @p: degrees of freedom.
 * @c: noncentrality parameter.
 * @q: probability.
 *
 * Calculates the @q-th quantile of the noncentral Student t
 * distribution with @c dof and noncentrality parameter equal to @c
 * via a rough and not particularly clever root-finding
 * algorithm. Maybe this can be more efficient by using logs. Some
 * experimentation needed.
 *
 * Returns: the calculated quantile, or #NADBL on failure.
 */

private double nc_student_cdf_inverse(double p, double c, double q)
{
    double x, d0, d1;
    int iter, subiter;
    double F, f, dir;

    if (p < 1 || c < 0 || q <= 0 || q >= 1) {
        return NADBL;
    }

    if (fabs(c) < 1.0e-10) {
        /* don't bother for infinitesimal c */
        return student_cdf_inverse(p, q);
    }

    x = c + student_cdf_inverse(p, q) / sqrt(p - 0.5);
    d0 = 1.0e7;
    iter = 0;

    while (fabs(d0) > 1.0e-10 && iter < 1000) {
        F = nc_student_cdf(p, c, x);
        f = nc_student_pdf(p, c, x);
        d0 = F - q;
        dir = d0/f;
        d1 = 1.0e7;
        subiter = 0;

        while (fabs(d1) > fabs(d0) && subiter < 100) {
            d1 = F - nc_student_cdf(p, c, x - dir);
            dir /= 2.0;
            subiter++;
        }

        if (subiter >= 100) {
            x = NADBL;
            break;
        } else {
            x -= dir*2;
            d0 = d1;
            iter++;
        }
    }

    if (iter >= 1000) {
        x = NADBL;
    }

    return x;
}

private struct Distmap {
    int code;
    const(char) * s;
}

int dist_code_from_string(const(char) * s)
{
    Distmap[26] dmap = [
        Distmap(DistCode.uniform,   "u"),
        Distmap(DistCode.udiscrt,   "i"),
        Distmap(DistCode.normal,    "z"),
        Distmap(DistCode.student,   "t"),
        Distmap(DistCode.chisq,     "x"),
        Distmap(DistCode.snedecor,  "f"),
        Distmap(DistCode.binomial,  "b"),
        Distmap(DistCode.poisson,   "p"),
        Distmap(DistCode.expon,     "exp"),
        Distmap(DistCode.weibull,   "w"),
        Distmap(DistCode.gamma,     "g"),
        Distmap(DistCode.ged,       "e"),
        Distmap(DistCode.laplace,   "l"),
        Distmap(DistCode.beta,      "beta"),
        Distmap(DistCode.dw,        "d"),
        Distmap(DistCode.binorm,    "D"),
        Distmap(DistCode.johansen,  "J"),
        Distmap(DistCode.betabin,   "bb"),
        Distmap(DistCode.nc_chisq,  "ncx"),
        Distmap(DistCode.nc_f,      "ncf"),
        Distmap(DistCode.nc_t,      "nct"),
        Distmap(DistCode.logistic,  "s"),
        Distmap(DistCode.dirichlet, "dir"),
        Distmap(DistCode.discrete,  "disc"),
        Distmap(DistCode.tnormal,   "tn"),
        Distmap(DistCode.none,      null)
    ];
    char[8] test;
    int i;

    if (!strcmp(s, "D")) {
        /* special: case counts for bivariate normal */
        return DistCode.binorm;
    }

    /* otherwise we'll ignore case */
    for (i=0; i<8 && s[i]; i++) {
        test[i] = cast(char) tolower(s[i]);
    }
    test[i] = '\0';

    for (i=0; dmap[i].code > 0; i++) {
        if (!strcmp(test.ptr, dmap[i].s)) {
            return dmap[i].code;
        }
    }

    /* backward compatibility */
    if (!strcmp(test.ptr, "n")) {
        return DistCode.normal;
    } else if (!strcmp(test.ptr, "c")) {
        return DistCode.chisq;
    } else if (!strcmp(test.ptr, "lgt")) {
        return DistCode.logistic;
    }

    return DistCode.none;
}

/**
 * print_critval:
 * @dist: distribution code.
 * @parm: array holding 0 to 2 parameter values.
 * @a: alpha.
 * @c: the critical value.
 * @prn: gretl printer.
 *
 * Prints the critical value information in a consistent manner.
 */

void print_critval(int dist, const(double) * parm, double a, double c, Prn * prn)
{
    switch (dist) {
    case DistCode.normal:
        pprintf(prn, "%s", _("Standard normal distribution").ptr);
        break;
    case DistCode.student:
        pprintf(prn, "t(%g)", parm[0]);
        break;
    case DistCode.chisq:
        pprintf(prn, "%s(%g)", _("Chi-square").ptr, parm[0]);
        break;
    case DistCode.snedecor:
        pprintf(prn, "F(%g, %g)", parm[0], parm[1]);
        break;
    case DistCode.binomial:
        pprintf(prn, _("Binomial (P = %g, %g trials)").ptr, parm[0], parm[1]);
        break;
    case DistCode.poisson:
        pprintf(prn, _("Poisson (mean = %g)").ptr, parm[0]);
        break;
    case DistCode.expon:
        pprintf(prn, _("Exponential (scale = %g)").ptr, parm[0]);
        break;
    case DistCode.weibull:
        pprintf(prn, _("Weibull (shape = %g, scale = %g)").ptr, parm[0], parm[1]);
        break;
    default:
        break;
    }

    pputs(prn, "\n ");
    pprintf(prn, _("right-tail probability = %g").ptr, a);
    pputs(prn, "\n ");
    pprintf(prn, _("complementary probability = %g").ptr, 1.0 - a);
    if (a < 0.5 && (dist == DistCode.normal || dist == DistCode.student)) {
        pputs(prn, "\n ");
        pprintf(prn, _("two-tailed probability = %g").ptr, 2.0 * a);
    }
    pputs(prn, "\n\n ");
    pprintf(prn, _("Critical value = %g").ptr, c);
    pputc(prn, '\n');
}

private int pdist_check_input(int dist, const(double) * parm,
                              double x)
{
    int i, np = 1; /* default */

    if (na(x)) {
        return GretlError.missdata;
    }

    if (dist == DistCode.normal) {
        np = 0;
    } else if (dist == DistCode.snedecor || dist == DistCode.gamma ||
               dist == DistCode.binomial || dist == DistCode.weibull ||
               dist == DistCode.nc_chisq || dist == DistCode.nc_t ||
               dist == DistCode.laplace  || dist == DistCode.beta) {
        np = 2;
    } else if (dist == DistCode.johansen || dist == DistCode.betabin ||
               dist == DistCode.nc_f) {
        np = 3;
    }

    for (i=0; i<np; i++) {
        if (na(parm[i])) {
            return GretlError.missdata;
        }
    }

    return 0;
}

/**
 * gretl_get_cdf_inverse:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @a: probability value.
 *
 * Returns: the argument, y, for which the area under the PDF
 * specified by @dist and @parm, integrated from its minimum to y,
 * is equal to @a, or #NADBL on failure.
 */

double gretl_get_cdf_inverse(int dist, const(double) * parm,
                             double a)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, a) == GretlError.missdata) {
        return y;
    }

    if (dist == DistCode.normal) {
        y = normal_cdf_inverse(a);
    } else if (dist == DistCode.student) {
        y = student_cdf_inverse(parm[0], a);
    } else if (dist == DistCode.chisq) {
        y = chisq_cdf_inverse(cast(int) parm[0], a);
    } else if (dist == DistCode.gamma) {
        y = gamma_cdf_inverse(parm[0], parm[1], a);
    } else if (dist == DistCode.snedecor) {
        y = snedecor_cdf_inverse(parm[0], parm[1], a);
    } else if (dist == DistCode.binomial) {
        y = binomial_cdf_inverse(parm[0], cast(int) parm[1], a);
    } else if (dist == DistCode.poisson) {
        y = poisson_cdf_inverse(parm[0], a);
    } else if (dist == DistCode.ged) {
        y = GED_cdf_inverse(parm[0], a);
    } else if (dist == DistCode.laplace) {
        y = laplace_cdf_inverse(parm[0], parm[1], a);
    } else if (dist == DistCode.nc_f) {
        y = nc_snedecor_cdf_inverse(parm[0], parm[1], parm[2], a);
    } else if (dist == DistCode.nc_chisq) {
        y = nc_chisq_cdf_inverse(parm[0], parm[1], a);
    } else if (dist == DistCode.nc_t) {
        y = nc_student_cdf_inverse(parm[0], parm[1], a);
    }

    return y;
}

/**
 * gretl_get_critval:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @a: right-tail probability.
 *
 * Returns: the abcsissa value x for the distribution specified
 * by @dist and @parm, such that P(X >= x) = @a, or #NADBL on
 * failure.
 */

double gretl_get_critval(int dist, const(double) * parm, double a)
{
    double x = NADBL;

    if (pdist_check_input(dist, parm, a) == GretlError.missdata) {
        return x;
    }

    if (dist == DistCode.normal) {
        x = normal_critval(a);
    } else if (dist == DistCode.student) {
        x = student_critval(parm[0], a);
    } else if (dist == DistCode.chisq) {
        x = chisq_critval(cast(int) parm[0], a);
    } else if (dist == DistCode.snedecor) {
        x = snedecor_critval(cast(int) parm[0], cast(int) parm[1], a);
    } else if (dist == DistCode.binomial) {
        x = binomial_critval(parm[0], cast(int) parm[1], a);
    } else if (dist == DistCode.poisson) {
        x = poisson_critval(parm[0], a);
    } else if (dist == DistCode.weibull) {
        x = weibull_critval(parm[0], parm[1], a);
    } else if (dist == DistCode.gamma) {
        x = gamma_cdf_inverse(parm[0], parm[1], 1-a);
    } else if (dist == DistCode.expon) {
        /* special case of Weibull */
        x = weibull_critval(1.0, parm[0], a);
    } else if (dist == DistCode.ged) {
        x = GED_cdf_inverse(parm[0], 1-a);
    } else if (dist == DistCode.laplace) {
        x = laplace_cdf_inverse(parm[0], parm[1], 1-a);
    }

    return x;
}

/**
 * gretl_get_cdf:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the CDF for the distribution specified by
 * @dist and @parm applicable at @x.
 *
 * Returns: the CDF value, or #NADBL on error.
 */

double gretl_get_cdf(int dist, const(double) * parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == GretlError.missdata) {
        return y;
    }

    if (dist == DistCode.normal) {
        y = normal_cdf(x);
    } else if (dist == DistCode.student) {
        y = student_cdf(parm[0], x);
    } else if (dist == DistCode.chisq) {
        y = chisq_cdf(cast(int) parm[0], x);
    } else if (dist == DistCode.snedecor) {
        y = snedecor_cdf(cast(int) parm[0], cast(int) parm[1], x);
    } else if (dist == DistCode.gamma) {
        y = gamma_cdf(parm[0], parm[1], x, 1);
    } else if (dist == DistCode.binomial) {
        y = binomial_cdf(parm[0], cast(int) parm[1], cast(int) x);
    } else if (dist == DistCode.poisson) {
        y = poisson_cdf(parm[0], cast(int) x);
    } else if (dist == DistCode.expon) {
        y = exponential_cdf(parm[0], x);
    } else if (dist == DistCode.weibull) {
        y = weibull_cdf(parm[0], parm[1], x);
    } else if (dist == DistCode.ged) {
        y = GED_cdf(parm[0], x);
    } else if (dist == DistCode.laplace) {
        y = laplace_cdf(parm[0], parm[1], x);
    } else if (dist == DistCode.nc_chisq) {
        y = nc_chisq_cdf(parm[0], parm[1], x);
    } else if (dist == DistCode.nc_f) {
        y = nc_snedecor_cdf(parm[0], parm[1], parm[2], x);
    } else if (dist == DistCode.nc_t) {
        y = nc_student_cdf(parm[0], parm[1], x);
    } else if (dist == DistCode.logistic) {
        y = logistic_cdf(x);
    } else if (dist == DistCode.beta) {
        y = beta_cdf(parm[0], parm[1], x);
    } else if (dist == DistCode.tnormal) {
        y = tnormal_cdf(parm[0], parm[1], x);
    }

    return y;
}

/**
 * gretl_get_pdf:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Evaluates the PDF for the distribution specified by
 * @dist and @parm at @x.
 *
 * Returns: the PDF value, or #NADBL on error.
 */

double gretl_get_pdf(int dist, const(double) * parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == GretlError.missdata) {
        return y;
    }

    if (dist == DistCode.normal) {
        y = normal_pdf(x);
    } else if (dist == DistCode.student) {
        y = student_pdf(parm[0], x);
    } else if (dist == DistCode.chisq) {
        y = chisq_pdf(parm[0], x);
    } else if (dist == DistCode.snedecor) {
        y = snedecor_pdf(cast(int) parm[0], cast(int) parm[1], x);
    } else if (dist == DistCode.gamma) {
        y = gamma_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.binomial) {
        y = binomial_pmf(parm[0], cast(int) parm[1], cast(int) x);
    } else if (dist == DistCode.poisson) {
        y = poisson_pmf(parm[0], cast(int) x);
    } else if (dist == DistCode.expon) {
        y = exponential_pdf(parm[0], x);
    } else if (dist == DistCode.weibull) {
        y = weibull_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.ged) {
        y = GED_pdf(parm[0], x);
    } else if (dist == DistCode.laplace) {
        y = laplace_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.nc_f) {
        y = nc_snedecor_pdf(parm[0], parm[1], parm[2], x);
    } else if (dist == DistCode.nc_t) {
        y = nc_student_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.nc_chisq) {
        y = nc_chisq_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.beta) {
        y = beta_pdf(parm[0], parm[1], x);
    } else if (dist == DistCode.tnormal) {
        y = tnormal_pdf(parm[0], parm[1], x);
    }

    return y;
}

/**
 * gretl_fill_pdf_array:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: see below.
 * @n: number of elements in @x.
 *
 * On input, @x contains an array of abscissae at which the
 * PDF specified by @dist and @parm should be evaluated. On
 * output it contains the corresponding PDF values.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_fill_pdf_array(int dist, const(double) * parm,
                         double * x, int n)
{
    int err = GretlError.data;

    if (pdist_check_input(dist, parm, 0) == GretlError.missdata) {
        return GretlError.missdata;
    }

    if (dist == DistCode.normal) {
        err = normal_pdf_array(x, n);
    } else if (dist == DistCode.student) {
        err = student_pdf_array(parm[0], x, n);
    } else if (dist == DistCode.chisq) {
        err = chisq_pdf_array(parm[0], x, n);
    } else if (dist == DistCode.snedecor) {
        err = snedecor_pdf_array(cast(int) parm[0], cast(int) parm[1], x, n);
    } else if (dist == DistCode.gamma) {
        err = gamma_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.binomial) {
        err = binomial_pmf_array(parm[0], cast(int) parm[1], x, n);
    } else if (dist == DistCode.poisson) {
        err = poisson_pmf_array(parm[0], x, n);
    } else if (dist == DistCode.expon) {
        err = exponential_pdf_array(parm[0], x, n);
    } else if (dist == DistCode.weibull) {
        err = weibull_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.ged) {
        err = GED_pdf_array(parm[0], x, n);
    } else if (dist == DistCode.laplace) {
        err = laplace_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.nc_f) {
        err = nc_snedecor_pdf_array(parm[0], parm[1], parm[2], x, n);
    } else if (dist == DistCode.nc_t) {
        err = nc_student_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.nc_chisq) {
        err = nc_chisq_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.beta) {
        err = beta_pdf_array(parm[0], parm[1], x, n);
    } else if (dist == DistCode.tnormal) {
        err = tnormal_pdf_array(parm[0], parm[1], x, n);
    }

    return err;
}

/* This apparatus remembers the parameters (p) and argument (x) from the
   last internal p-value assessment, via gretl_get_pvalue or
   batch_pvalue().
*/

private __gshared double[3] pvargs = 0;

private void remember_pvalue_args(const(double) * p, double x)
{
    pvargs[0] = p[0];
    pvargs[1] = p[1];
    pvargs[2] = x;
}

/**
 * gretl_get_pvalue:
 * @dist: distribution code.
 * @parm: array holding from zero to two parameter values,
 * depending on the distribution.
 * @x: abscissa value.
 *
 * Returns: the integral of the PDF specified by @dist and
 * @parm from @x to infinity, or #NADBL on error.
 */

double gretl_get_pvalue(int dist, const(double) * parm, double x)
{
    double y = NADBL;

    if (pdist_check_input(dist, parm, x) == GretlError.missdata) {
        return y;
    }

    if (dist == DistCode.normal) {
        y = normal_cdf_comp(x);
    } else if (dist == DistCode.student) {
        y = student_cdf_comp(parm[0], x);
    } else if (dist == DistCode.chisq) {
        y = chisq_cdf_comp(cast(int) parm[0], x);
    } else if (dist == DistCode.snedecor) {
        y = snedecor_cdf_comp(parm[0], parm[1], x);
    } else if (dist == DistCode.gamma) {
        y = gamma_cdf_comp(parm[0], parm[1], x, 1);
    } else if (dist == DistCode.binomial) {
        y = binomial_cdf_comp(parm[0], cast(int) parm[1], cast(int) x);
    } else if (dist == DistCode.poisson) {
        y = poisson_cdf_comp(parm[0], cast(int) x);
    } else if (dist == DistCode.expon) {
        y = weibull_cdf_comp(1.0, parm[0], x);
    } else if (dist == DistCode.weibull) {
        y = weibull_cdf_comp(parm[0], parm[1], x);
    } else if (dist == DistCode.ged) {
        y = GED_cdf_comp(parm[0], x);
    } else if (dist == DistCode.laplace) {
        y = laplace_cdf_comp(parm[0], parm[1], x);
    } else if (dist == DistCode.johansen) {
        y = johansen_trace_pval(cast(int) parm[0], cast(int) parm[1],
                                cast(int) parm[2], x);
    }

    remember_pvalue_args(parm, x);

    return y;
}

private int gretl_fill_random_array(double * x, int t1, int t2,
                                    int dist, const(double) * parm,
                                    const(double) * vecp1,
                                    const(double) * vecp2)
{
    int t, err = 0;

    if (dist == DistCode.uniform) {
        /* uniform, continuous */
        double min = parm[0], max = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) min = vecp1[t];
                if (vecp2 != null) max = vecp2[t];
                err = gretl_rand_uniform_minmax(x, t, t, min, max);
            }
        } else {
            err = gretl_rand_uniform_minmax(x, t1, t2, min, max);
        }
    } else if (dist == DistCode.udiscrt) {
        /* uniform, discrete */
        int min = cast(int) parm[0], max = cast(int) parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) min = cast(int) vecp1[t];
                if (vecp2 != null) max = cast(int) vecp2[t];
                err = gretl_rand_uniform_int_minmax(x, t, t, min, max,
                                                    Gretlopt.OPT_NONE);
            }
        } else {
            err = gretl_rand_uniform_int_minmax(x, t1, t2, min, max,
                                                Gretlopt.OPT_NONE);
        }
    } else if (dist == DistCode.normal) {
        double mu = parm[0], sd = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) mu = vecp1[t];
                if (vecp2 != null) sd = vecp2[t];
                err = gretl_rand_normal_full(x, t, t, mu, sd);
            }
        } else {
            err = gretl_rand_normal_full(x, t1, t2, mu, sd);
        }
    } else if (dist == DistCode.student) {
        /* Student's t */
        double v = parm[0];

        if (vecp1 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                v = vecp1[t];
                err = gretl_rand_student(x, t, t, v);
            }
        } else {
            err = gretl_rand_student(x, t1, t2, v);
        }
    } else if (dist == DistCode.chisq) {
        /* chi-square */
        int v = cast(int) parm[0];

        if (vecp1 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                v = cast(int) vecp1[t];
                err = gretl_rand_chisq(x, t, t, v);
            }
        } else {
            err = gretl_rand_chisq(x, t1, t2, v);
        }
    } else if (dist == DistCode.snedecor) {
        int v1 = cast(int) parm[0], v2 = cast(int) parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) v1 = cast(int) vecp1[t];
                if (vecp2 != null) v2 = cast(int) vecp2[t];
                err = gretl_rand_F(x, t, t, v1, v2);
            }
        } else {
            err = gretl_rand_F(x, t1, t2, v1, v2);
        }
    } else if (dist == DistCode.gamma) {
        double shape = parm[0], scale = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) shape = vecp1[t];
                if (vecp2 != null) scale = vecp2[t];
                err = gretl_rand_gamma(x, t, t, shape, scale);
            }
        } else {
            err = gretl_rand_gamma(x, t1, t2, shape, scale);
        }
    } else if (dist == DistCode.binomial) {
        double pr = parm[0];
        int n = cast(int) parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) pr = vecp1[t];
                if (vecp2 != null) n = cast(int) vecp2[1];
                err = gretl_rand_binomial(x, t, t, n, pr);
            }
        } else {
            err = gretl_rand_binomial(x, t1, t2, n, pr);
        }
    } else if (dist == DistCode.poisson) {
        double m = parm[0];

        if (vecp1 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                m = vecp1[t];
                err = gretl_rand_poisson(x, t, t, &m, 0);
            }
        } else {
            err = gretl_rand_poisson(x, t1, t2, &m, 0);
        }
    } else if (dist == DistCode.expon) {
        double scale = parm[0];

        if (vecp1 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                scale = vecp1[t];
                err = gretl_rand_exponential(x, t, t, scale);
            }
        } else {
            err = gretl_rand_exponential(x, t1, t2, scale);
        }
    } else if (dist == DistCode.weibull) {
        double shape = parm[0], scale = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) shape = vecp1[t];
                if (vecp2 != null) scale = vecp2[t];
                err = gretl_rand_weibull(x, t, t, shape, scale);
            }
        } else {
            err = gretl_rand_weibull(x, t1, t2, shape, scale);
        }
    } else if (dist == DistCode.ged) {
        double nu = parm[0];

        if (vecp1 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) nu = vecp1[t];
                err = gretl_rand_GED(x, t, t, nu);
            }
        } else {
            err = gretl_rand_GED(x, t1, t2, nu);
        }
    } else if (dist == DistCode.laplace) {
        double mu = parm[0], b = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) mu = vecp1[t];
                if (vecp2 != null) b = vecp2[t];
                err = gretl_rand_laplace(x, t, t, mu, b);
            }
        } else {
            err = gretl_rand_laplace(x, t1, t2, mu, b);
        }
    } else if (dist == DistCode.beta) {
        double shape1 = parm[0], shape2 = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) shape1 = vecp1[t];
                if (vecp2 != null) shape2 = vecp2[t];
                err = gretl_rand_beta(x, t, t, shape1, shape2);
            }
        } else {
            err = gretl_rand_beta(x, t1, t2, shape1, shape2);
        }
    } else if (dist == DistCode.betabin) {
        double shape1 = parm[1], shape2 = parm[2];
        int n = cast(int) parm[0];

        err = gretl_rand_beta_binomial(x, t1, t2, n, shape1, shape2);
    } else if (dist == DistCode.logistic) {
        double loc = parm[0], shape = parm[1];

        err = gretl_rand_logistic(x, t1, t2, loc, shape);
    } else if (dist == DistCode.tnormal) {
        double lo = parm[0], hi = parm[1];

        if (vecp1 != null || vecp2 != null) {
            for (t=t1; t<=t2 && !err; t++) {
                if (vecp1 != null) lo = vecp1[t];
                if (vecp2 != null) hi = vecp2[t];
                err = gretl_rand_tnormal(x, t, t, lo, hi);
            }
        } else {
            err = gretl_rand_tnormal(x, t1, t2, lo, hi);
        }
    }

    return err;
}

/**
 * gretl_fill_random_series:
 * @x: series to fill (must be of length dset->n).
 * @dist: distribution code.
 * @parm: array holding either one or two scalar
 * parameter values, depending on the distribution.
 * @vecp1: series containing values for first param,
 * or %NULL.
 * @vecp2: series containing values for second param,
 * or %NULL.
 * @dset: dataset information.
 *
 * Fills @x with random values conforming to the distribution
 * given by @dist, which may require specification of either
 * one or two parameters. These parameters are either
 * given as (scalar) elements of @parm, or they may vary
 * by observation, in which case they are given as the
 * elements of @vecp1 or @vecp2.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_fill_random_series(double * x, int dist,
                             const(double) * parm,
                             const(double) * vecp1,
                             const(double) * vecp2,
                             const(DatasetRaw) * dset)
{
    return gretl_fill_random_array(x, dset.t1, dset.t2,
                                   dist, parm, vecp1, vecp2);
}

GretlMatrixRaw * gretl_get_random_matrix(int dist,
                                      const(double) * parm,
                                      const(double) * vecp1,
                                      const(double) * vecp2,
                                      int rows, int cols,
                                      int * err)
{
    GretlMatrixRaw * m = null;
    int n = rows * cols;

    if (n <= 0) {
        *err = GretlError.invarg;
    } else {
        m = gretl_matrix_alloc(rows, cols);
        if (m == null) {
            *err = GretlError.alloc;
            return null;
        } else {
            *err = gretl_fill_random_array(m.val, 0, n-1, dist,
                                           parm, vecp1, vecp2);
        }
    }

    return m;
}

double gretl_get_random_scalar(int dist, const(double) * parm,
                               int * err)
{
    double x;

    *err = gretl_fill_random_array(&x, 0, 0, dist,
                                   parm, null, null);

    return x;
}

private int print_pv_string(double x, double p, Prn * prn)
{
    char[32] numstr;

    if (na(p)) {
        pprintf(prn, _("area to the right of %g: NA\n").ptr, x);
        return 1;
    }

    sprintf(numstr.ptr, "%g", p);

    if (!strcmp(numstr.ptr, "1") || !strcmp(numstr.ptr, "0")) {
        pprintf(prn, _("area to the right of %g =~ %g\n").ptr, x, p);
    } else {
        pprintf(prn, _("area to the right of %g = %g\n").ptr, x, p);
    }

    return 0;
}

/**
 * print_pvalue:
 * @dist: distribution code.
 * @parm: array holding 1 or 2 parameter values.
 * @x: the value in the distribution.
 * @pv: the p-value.
 * @prn: gretl printer.
 *
 * Prints the p-value information in a consistent manner.
 */

void print_pvalue(int dist, const(double) * parm, double x,
                  double pv, Prn * prn)
{
    double pc;
    int err;

    switch (dist) {

    case DistCode.normal:
        pprintf(prn, "%s: ", _("Standard normal").ptr);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        if (pv < 0.5) {
            pprintf(prn, _("(two-tailed value = %g; complement = %g)\n").ptr,
                    2 * pv, 1 - 2 * pv);
        } else {
            pc = normal_cdf(x);
            pprintf(prn, _("(to the left: %g)\n").ptr, pc);
            pprintf(prn, _("(two-tailed value = %g; complement = %g)\n").ptr,
                    2 * pc, 1 - 2 * pc);
        }
        break;

    case DistCode.student:
        pprintf(prn, "t(%d): ", cast(int) parm[0]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        if (pv < 0.5) {
            pprintf(prn, _("(two-tailed value = %g; complement = %g)\n").ptr,
                    2 * pv, 1 - 2 * pv);
        } else {
            pc = student_cdf(parm[0], x);
            pprintf(prn, _("(to the left: %g)\n").ptr, pc);
            pprintf(prn, _("(two-tailed value = %g; complement = %g)\n").ptr,
                    2 * pc, 1 - 2 * pc);
        }
        break;

    case DistCode.chisq:
        pprintf(prn, "%s(%d): ", _("Chi-square").ptr, cast(int) parm[0]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = chisq_cdf(parm[0], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    case DistCode.snedecor:
        pprintf(prn, "F(%d, %d): ", cast(int) parm[0], cast(int) parm[1]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = snedecor_cdf(cast(int) parm[0], cast(int) parm[1], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    case DistCode.gamma:
        pprintf(prn, _("Gamma (shape %g, scale %g, mean %g, variance %g):"
                       ~ "\n area to the right of %g = %g\n").ptr,
                parm[0], parm[1], parm[0] * parm[1],
                parm[0] * parm[1] * parm[1],
                x, pv);
        break;

    case DistCode.binomial:
        pprintf(prn, _("Binomial (p = %g, n = %d):"
                       ~ "\n Prob(x > %d) = %g\n").ptr,
                parm[0], cast(int) parm[1], cast(int) x, pv);
        pc = binomial_cdf(parm[0], cast(int) parm[1], cast(int) x);
        if (x > 0) {
            double pm1 = binomial_cdf(parm[0], cast(int) parm[1], cast(int) (x - 1));

            pprintf(prn, _(" Prob(x <= %d) = %g\n").ptr, cast(int) x, pc);
            pprintf(prn, _(" Prob(x >= %d) = %g\n").ptr, cast(int) x, 1.0 - pm1);
            pprintf(prn, _(" Prob(x = %d) = %g\n").ptr, cast(int) x, pc - pm1);
        } else {
            pprintf(prn, _(" Prob(x = %d) = %g\n").ptr, cast(int) x, pc);
        }
        break;

    case DistCode.poisson:
        pprintf(prn, _("Poisson (mean = %g): ").ptr, parm[0]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = poisson_cdf(parm[0], cast(int) x);
        if (x > 0) {
            pprintf(prn, _(" Prob(x <= %d) = %g\n").ptr, cast(int) x, pc);
            pprintf(prn, _(" Prob(x = %d) = %g\n").ptr, cast(int) x,
                    poisson_pmf(parm[0], cast(int) x));
        } else {
            pprintf(prn, _(" Prob(x = %d) = %g\n").ptr, cast(int) x, pc);
        }
        break;

    case DistCode.expon:
        pprintf(prn, _("Exponential (scale = %g): ").ptr, parm[0]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = exponential_cdf(parm[0], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    case DistCode.weibull:
        pprintf(prn, _("Weibull (shape = %g, scale = %g): ").ptr,
                parm[0], parm[1]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = weibull_cdf(parm[0], parm[1], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    case DistCode.ged:
        pprintf(prn, _("GED (shape = %g): ").ptr, parm[0]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = GED_cdf(parm[0], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    case DistCode.laplace:
        pprintf(prn, _("Laplace (mean = %g, scale = %g): ").ptr, parm[0], parm[1]);
        err = print_pv_string(x, pv, prn);
        if (err) return;
        pc = laplace_cdf(parm[0], parm[1], x);
        pprintf(prn, _("(to the left: %g)\n").ptr, pc);
        break;

    default:
        break;
    }
}

/**
 * batch_pvalue:
 * @str: the command line, which should be of one of the following forms:
 * pvalue z x (Normal distribution);
 * pvalue t df x (t-distribution);
 * pvalue X df x (Chi-square);
 * pvalue F dfn dfd x (F-distribution); or
 * pvalue G mean variance x (Gamma distribution).
 * pvalue B prob n x (Binomial distribution).
 * pvalue P mean k (Poisson distribution).
 * pvalue W shape scale x (Weibull distribution).
 * @dset: dataset struct.
 * @prn: gretl printing struct.
 *
 * Calculates and prints the probability that a random variable
 * distributed as specified in the command line @str exceeds the
 * value indicated in @str.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int batch_pvalue(const(char) * str, DatasetRaw * dset, Prn * prn)
{
    double pv = NADBL;
    char[MAXLEN] line;
    char ** S;
    int dist;
    int i, n, m;
    int err = 0;

    if (str == null || *str == '\0') {
        return GretlError.args;
    }

    if (!strncmp(str, "pvalue ", 7)) {
        str += 7;
    }

    S = gretl_string_split(str, &n, null);
    if (S == null) {
        return GretlError.alloc;
    }

    dist = dist_code_from_string(S[0]);

    if (dist == DistCode.none) {
        err = GretlError.invarg;
    } else {
        strcpy(line.ptr, "pvalue(");
        m = 8;
        for (i=0; i<n && !err; i++) {
            m += cast(int) strlen(S[i]) + 1;
            if (m > MAXLEN) {
                err = GretlError.data;
            } else {
                strcat(line.ptr, S[i]);
                strcat(line.ptr, (i == n - 1)? ")" : ",");
            }
        }
    }

    strings_array_free(S, n);

    if (!err) {
        pv = generate_scalar(line.ptr, dset, &err);
    }

    if (!err) {
        print_pvalue(dist, pvargs.ptr, pvargs[2], pv, prn);
    }

    return err;
}

/**
 * gretl_get_DW:
 * @n: number of observations
 * @k: number of regressors excluding the constant.
 * @err: location to receive error code.
 *
 * Consults a table of Durbin-Watson critical values and
 * returns the results in a gretl_matrix.
 *
 * Returns: on success, a 4-vector containing the lower
 * and upper Durbin-Watson values, dl and du, along with
 * the values actually used for @n and @k (which may differ
 * from those given on input if the exact values are not
 * found in the table and have to be approximated).
 * On error, returns %NULL.
 */

GretlMatrixRaw * gretl_get_DW(int n, int k, int * err)
{
    DwLookupFn dw_lookup;
    GretlMatrixRaw * m = null;

    dw_lookup = cast(DwLookupFn) get_plugin_function("dw_lookup");

    if (dw_lookup == null) {
        *err = GretlError.fopen;
        return null;
    }

    *err = (*dw_lookup)(n, k, &m);

    return m;
}
