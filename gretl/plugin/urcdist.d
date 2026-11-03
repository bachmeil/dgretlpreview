/* Port of plugin/urcdist.c: driver for James MacKinnon's "urcval" function, to calculate p-values for unit root tests.

   See James G. MacKinnon, "Numerical Distribution Functions for Unit Root and Cointegration Tests", Journal of Applied
   Econometrics, Vol. 11, No. 6, 1996, pp. 601-618, and also http://qed.econ.queensu.ca/pub/faculty/mackinnon/numdist/

   The calculation code here is Copyright (c) James G. MacKinnon, 1996 (corrected 2003-5-5). This "wrapper" written by
   Allin Cottrell, 2004; revised to load MacKinnon's tables in binary format, 2015.

   In the C this is a gretl plugin loaded with get_plugin_function("mackinnon_pvalue") (and "dfgls_pvalue"), which reads
   its tables from binary files in the "data" directory of the plugin path. Here it is an ordinary module that plugins.d
   returns for those names (both functions are extern (C) so that their addresses have the plugin-function types used by
   adf_kpss.d), and the tables are D arrays compiled in (gretl.plugin.urcdata, generated from the .bin files by
   gretl/data/make_urcdata_d.py): the file open, seek and read of the C, and the byte swapping for big-endian machines,
   are replaced by copies out of those arrays. Nothing else is changed. */

module gretl.plugin.urcdist;

import core.stdc.math : fabs, pow, sqrt;
import core.stdc.stdio : fprintf, stderr;
import core.stdc.stdlib : free, malloc;
import core.stdc.string : memcpy, strcmp;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.missingh : NADBL;
import gretl.plugin.urcdata : EmbeddedMatrix, embedded_matrices, urcdata;
import gretl.pvalues : normal_cdf, normal_cdf_inverse;

enum NIVMAX = 8;
enum URCLEN = 221;
enum BIGLEN = 884;

/* Based on Fortran code copyright (c) James G. MacKinnon,
   1995. Routine to evaluate response surface for specified betas and
   sample size.
*/

private void eval_all_crit(double tau, double * b, int nb, int T, double * crit, int * imin)
{
    double d1 = 0, d2 = 0, d3 = 0;
    double diff, diffm = 1000;
    int i;

    if (T > 0) {
        d1 = 1.0 / T;
        d2 = d1 * d1;
        d3 = d1 * d2;
    }

    for (i=0; i<URCLEN; i++) {
        if (T == 0) {
            crit[i] = b[0];
        } else if (nb == 3) {
            crit[i] = b[0] + b[1]*d1 + b[2]*d2;
        } else if (nb == 4) {
            crit[i] = b[0] + b[1]*d1 + b[2]*d2 + b[3]*d3;
        }
        diff = fabs(tau - crit[i]);
        if (diff < diffm) {
            diffm = diff;
            *imin = i + 1; /* has to be 1-based */
        }
        b += nb;
    }
}

/* Copyright (c) James G. MacKinnon, 1993.  This routine uses the
   Cholesky decomposition to invert a real symmetric matrix.
*/

private int cholx(double * a, int m, int n)
{
    int i, j, k, kl;
    double t, ooa = 0.0;
    int err = 0;

    /* Parameter adjustment */
    a -= 1 + m;

    for (i = 1; i <= n; ++i) {
        kl = i - 1;
        for (j = i; j <= n; ++j) {
            if (i > 1) {
                for (k = 1; k <= kl; ++k) {
                    a[i + j * m] -= a[k + i * m] * a[k + j * m];
                }
            } else if (a[i + i * m] <= 0.0) {
                /* error: get out */
                fprintf(stderr, "cholx: failed at i = %d\n", i);
                err = GretlError.notpd;
                goto cholx_exit;
            }
            if (i == j) {
                a[i + i * m] = sqrt(a[i + i * m]);
            } else {
                if (j == i + 1) {
                    ooa = 1. / a[i + i * m];
                }
                a[i + j * m] *= ooa;
            }
        }
    }

    for (i = 1; i <= n; ++i) {
        for (j = i; j <= n; ++j) {
            ooa = 1.0 / a[j + j * m];
            if (i >= j) {
                t = 1.0;
            } else {
                kl = j - 1;
                t = 0.0;
                for (k = i; k <= kl; ++k) {
                    t -= a[i + k * m] * a[k + j * m];
                }
            }
            a[i + j * m] = t * ooa;
        }
    }

    for (i = 1; i <= n; ++i) {
        for (j = i; j <= n; ++j) {
            t = 0.0;
            for (k = j; k <= n; ++k) {
                t += a[i + k * m] * a[j + k * m];
            }
            a[i + j * m] = t;
            a[j + i * m] = t;
        }
    }

 cholx_exit:

    return err;
}

/* Copyright (c) James G. MacKinnon, 1995.  Subroutine to do GLS
   estimation the obvious way.  Use only when sample size is small
   (nobs <= 50). 1995-1-3
*/

private int gls(double * xmat, double * yvec, double * omega,
                double * beta, double * xomx, double * fits,
                double * resid, double * ssr, double * ssrt,
                int T, int ivrt)
{
    int nomax = 20;
    int nvmax = 4;
    int nvar = 4 - ivrt;
    int omega_offset = 1 + nomax;
    int xomx_offset = 1 + nvmax;
    int xmat_offset = 1 + nomax;
    int i, j, k, l;
    double[50] xomy;
    int err = 0;

    /* xomx is covariance matrix of parameter estimates if omega is
       truly known. First, invert omega matrix if ivrt=0. Original one
       gets replaced.
    */

    /* parameter adjustments */
    omega -= omega_offset;
    xomx -= xomx_offset;
    xmat -= xmat_offset;
    --resid;
    --fits;
    --yvec;
    --beta;

    if (ivrt == 0) {
        err = cholx(&omega[omega_offset], nomax, T);
        if (err) {
            return err;
        }
    }

    /* form xomx matrix and xomy vector */
    for (j = 1; j <= nvar; ++j) {
        xomy[j - 1] = 0.;
        for (l = j; l <= nvar; ++l) {
            xomx[j + l * nvmax] = 0.;
        }
    }
    for (i = 1; i <= T; ++i) {
        for (k = 1; k <= T; ++k) {
            for (j = 1; j <= nvar; ++j) {
                xomy[j - 1] += xmat[i + j * nomax] *
                    omega[k + i * nomax] * yvec[k];
                for (l = j; l <= nvar; ++l) {
                    xomx[j + l * nvmax] += xmat[i + j * nomax] *
                        omega[k + i * nomax] *
                        xmat[k + l * nomax];
                }
            }
        }
    }

    for (j = 1; j <= nvar; ++j) {
        for (l = j; l <= nvar; ++l) {
            xomx[l + j * nvmax] = xomx[j + l * nvmax];
        }
    }

    /* invert xomx matrix */
    err = cholx(&xomx[xomx_offset], nvmax, nvar);
    if (err) {
        return err;
    }

    /* form estimates of beta */
    for (i = 1; i <= nvar; ++i) {
        beta[i] = 0.0;
        for (j = 1; j <= nvar; ++j) {
            beta[i] += xomx[i + j * nvmax] * xomy[j - 1];
        }
    }

    /* find ssr, fitted values, and residuals */
    *ssr = 0.0;
    for (i = 1; i <= T; ++i) {
        fits[i] = 0.0;
        for (j = 1; j <= nvar; ++j) {
            fits[i] += xmat[i + j * nomax] * beta[j];
        }
        resid[i] = yvec[i] - fits[i];
        *ssr += resid[i] * resid[i];
    }

    /* find ssr from transformed regression */
    *ssrt = 0.0;
    for (i = 1; i <= T; ++i) {
        for (k = 1; k <= T; ++k) {
            *ssrt += resid[i] * omega[k + i * nomax] * resid[k];
        }
    }

    return err;
}

/* Based on Fortran code copyright (c) James G. MacKinnon, 1995.
   Routine to find P-value for any specified test statistic.
*/

private double fpval(double * beta, int nbeta, double * wght,
                     double * prob, double * cnorm,
                     double tau, int T, int * err)
{
    double d1, precrt = 2.0;
    int i, j, ic, jc, imin = 0;
    int np1, nph, nptop, np = 9;
    double bot, top, ssr, ssrt;
    double se3, ttest, crfit;
    double * yvec, fits, resid;
    double * xmat, xomx, gamma, omega;
    double[URCLEN] crits;
    double pval = 0.0;

    /* sizes of arrays:
       yvec[20], fits[20], resid[20]
       xmat[80], xomx[16], gamma[4], omega[400]
    */
    int wlen = 3*20 + 20 + 400;
    double * wspace = cast(double *) malloc(wlen * double.sizeof);

    if (wspace == null) {
        *err = GretlError.alloc;
        return NADBL;
    }

    /* set offsets into workspace */
    yvec = wspace;
    fits = yvec + 20;
    resid = fits + 20;
    xmat = resid + 20;
    xomx = xmat + 80;
    gamma = xomx + 16;
    omega = gamma + 4;

    /* first compute all the estimated critical values,
       and find the one closest to the test statistic,
       indexed by @imin.
    */
    eval_all_crit(tau, beta, nbeta, T, crits.ptr, &imin);

    nph = np / 2;
    nptop = URCLEN - nph;

    if (imin > nph && imin < nptop) {
        /* imin is not too close to the end.
           Use np points around tau.
        */
        for (i=1; i<=np; i++) {
            ic = imin - nph - 1 + i;
            yvec[i - 1] = cnorm[ic];
            xmat[i - 1] = 1.0;
            xmat[i + 19] = crits[ic - 1];
            xmat[i + 39] = xmat[i + 19] * crits[ic - 1];
            xmat[i + 59] = xmat[i + 39] * crits[ic - 1];
        }

        /* form omega matrix */
        for (i=1; i<=np; i++) {
            for (j=i; j<=np; j++) {
                ic = imin - nph - 1 + i;
                jc = imin - nph - 1 + j;
                top = prob[ic] * (1 - prob[jc]);
                bot = prob[jc] * (1 - prob[ic]);
                omega[i + j * 20 - 21] = wght[ic] * wght[jc] * sqrt(top / bot);
            }
        }
        for (i=1; i<=np; i++) {
            for (j=i; j<=np; j++) {
                omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
            }
        }

        *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
                   &ssr, &ssrt, np, 0);
        if (*err) {
            goto bailout;
        }

        /* check: is a cubic term actually needed? */
        se3 = sqrt(ssrt / (np - 4) * xomx[15]);
        ttest = fabs(gamma[3]) / se3;
        d1 = tau;

        if (ttest > precrt) {
            crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1) +
                gamma[3] * (d1 * d1 * d1);
        } else {
            *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
                       &ssr, &ssrt, np, 1);
            if (*err) {
                goto bailout;
            }
            crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1);
        }
        pval = normal_cdf(crfit);
    } else {
        /* imin is close to one of the ends. Use points from
           imin +/- nph to end.
        */
        if (imin < np) {
            np1 = imin + nph;
            if (np1 < 5) {
                np1 = 5;
            }
            for (i = 1; i <= np1; ++i) {
                yvec[i - 1] = cnorm[i];
                xmat[i - 1] = 1.0;
                xmat[i + 19] = crits[i - 1];
                xmat[i + 39] = xmat[i + 19] * crits[i-1];
                xmat[i + 59] = xmat[i + 39] * crits[i-1];
            }
        } else {
            np1 = (URCLEN + 1) - imin + nph;
            if (np1 < 5) {
                np1 = 5;
            }
            for (i = 1; i <= np1; ++i) {
                ic = (URCLEN + 1) - i;
                yvec[i - 1] = cnorm[ic];
                xmat[i - 1] = 1.;
                xmat[i + 19] = crits[ic - 1];
                xmat[i + 39] = xmat[i + 19] * crits[ic-1];
                xmat[i + 59] = xmat[i + 39] * crits[ic-1];
            }
        }

        /* form omega matrix */
        for (i = 1; i <= np1; ++i) {
            for (j = i; j <= np1; ++j) {
                if (imin < np) {
                    top = prob[i] * (1.0 - prob[j]);
                    bot = prob[j] * (1.0 - prob[i]);
                    omega[i + j * 20 - 21] = wght[i] * wght[j] * sqrt(top / bot);
                } else {
                    /* avoid numerical singularities at the upper end */
                    omega[i + j * 20 - 21] = 0.;
                    if (i == j) {
                        omega[i + i * 20 - 21] = 1.;
                    }
                }
            }
        }

        for (i = 1; i <= np1; ++i) {
            for (j = i; j <= np1; ++j) {
                omega[j + i * 20 - 21] = omega[i + j * 20 - 21];
            }
        }

        *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
                   &ssr, &ssrt, np1, 0);
        if (*err) {
            goto bailout;
        }

        /* is gamma[3] needed? */
        se3 = sqrt(ssrt / (np1 - 4) * xomx[15]);
        ttest = fabs(gamma[3]) / se3;
        d1 = tau;

        if (ttest > precrt) {
            crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1) +
                gamma[3] * (d1 * d1 * d1);
        } else {
            *err = gls(xmat, yvec, omega, gamma, xomx, fits, resid,
                       &ssr, &ssrt, np1, 1);
            if (*err) {
                goto bailout;
            }
            crfit = gamma[0] + gamma[1] * d1 + gamma[2] * (d1 * d1);
        }
        pval = normal_cdf(crfit);

        /* check that nothing crazy has happened at the ends */
        if (imin == 1 && pval > prob[1]) {
            pval = prob[1];
        }
        if (imin == URCLEN && pval < prob[URCLEN]) {
            pval = prob[URCLEN];
        }
    }

 bailout:

    free(wspace);
    if (*err) {
        pval = NADBL;
    }

    return pval;
}

private struct UrcInfo {
    int nz;
    int nreg;
    int model;
    int Tmin;
    int pos;
}

/*
   niv = # of integrated variables
   itv = appropriate ur_code for nc, c, ct, ctt models
   T = sample size (0 for asymptotic)
   tau = test statistic
   err = location to receive error code
*/

private double urcval(int niv, int itv, int T, double tau, int * err)
{
    double[BIGLEN] beta;
    double[URCLEN+1] wght;
    double[URCLEN+1] prob;
    double[URCLEN+1] cnorm;
    static immutable UrcInfo[34] uis = [
        UrcInfo(0,  1, 2, 20,      0), /* dfnc: table 1 */
        UrcInfo(0,  2, 2, 20,   7072), /* dfc */
        UrcInfo(0,  3, 3, 20,  14144), /* dfct */
        UrcInfo(0,  4, 3, 20,  22984), /* dfctt */
        UrcInfo(1,  2, 2, 20,  31824), /* conc: table 2 */
        UrcInfo(1,  3, 2, 20,  38896), /* coc */
        UrcInfo(1,  4, 3, 25,  45968), /* coct */
        UrcInfo(1,  5, 3, 20,  54808), /* coctt */
        UrcInfo(2,  3, 2, 25,  63648), /* conc: table 3 */
        UrcInfo(2,  4, 2, 20,  70720), /* coc */
        UrcInfo(2,  5, 2, 20,  77792), /* coct */
        UrcInfo(2,  6, 3, 20,  84864), /* coctt */
        UrcInfo(3,  4, 3, 20,  93704), /* conc: table 4 */
        UrcInfo(3,  5, 2, 25, 102544), /* coc */
        UrcInfo(3,  6, 3, 20, 109616), /* coct */
        UrcInfo(3,  7, 2, 30, 118456), /* coctt */
        UrcInfo(4,  5, 2, 25, 125528), /* conc: table 5 */
        UrcInfo(4,  6, 3, 20, 132600), /* coc */
        UrcInfo(4,  7, 3, 20, 141440), /* coct */
        UrcInfo(4,  8, 3, 20, 150280), /* coctt */
        UrcInfo(5,  6, 2, 30, 159120), /* conc: table 6 */
        UrcInfo(5,  7, 2, 30, 166192), /* coc */
        UrcInfo(5,  8, 2, 30, 173264), /* coct */
        UrcInfo(5,  9, 3, 25, 180336), /* coctt */
        UrcInfo(6,  7, 3, 25, 189176), /* conc: table 7 */
        UrcInfo(6,  8, 3, 25, 198016), /* coc */
        UrcInfo(6,  9, 3, 30, 206856), /* coct */
        UrcInfo(6, 10, 3, 30, 215696), /* coctt */
        UrcInfo(7,  8, 2, 40, 224536), /* conc: table 8 */
        UrcInfo(7,  9, 2, 35, 231608), /* coc */
        UrcInfo(7, 10, 2, 40, 238680), /* coct */
        UrcInfo(7, 11, 2, 40, 245752), /* coctt */
        UrcInfo(0,  0, 0,  0, 252824), /* prob */
        UrcInfo(0,  0, 0, -1, 254592)  /* cnorm */
    ];
    const(UrcInfo) * ui;
    int i, nbeta;
    double pval = NADBL;

    /* Check that parameters are valid */
    if (niv < 1 || niv > NIVMAX) {
        *err = GretlError.data;
        return pval;
    }

    if (itv < 1 || itv > 4) {
        /* these limits correspond to UR_NO_CONST and UR_QUAD_TREND
           in lib/src/adf_kpss.c */
        *err = GretlError.data;
        return pval;
    }

    /* "skip to appropriate location in data file": here, the index into the table of doubles
       (the offsets are in bytes of urcdata.bin) */
    i = (niv-1) * 4 + (itv - 1);
    ui = &uis[i];
    const size_t pos = ui.pos / double.sizeof;

    /* the number of coefficients in the critical values
       equations */
    nbeta = ui.model == 2 ? 3 : 4;

    /* these arrays are padded by one for fortran */
    wght[0] = prob[0] = cnorm[0] = 0.0;

    /* read coefficients and weights */
    memcpy(beta.ptr, &urcdata[pos], nbeta * URCLEN * double.sizeof);
    memcpy(wght.ptr + 1, &urcdata[pos + nbeta * URCLEN], URCLEN * double.sizeof);

    if (!*err) {
        /* read from embedded "probs.tab" data */
        const size_t ppos = uis[32].pos / double.sizeof;

        memcpy(prob.ptr + 1, &urcdata[ppos], URCLEN * double.sizeof);
        memcpy(cnorm.ptr + 1, &urcdata[ppos + URCLEN], URCLEN * double.sizeof);
    }

    if (!*err && T > 0 && T < ui.Tmin) {
        /* error, or warning? */
        fprintf(stderr, "Warning, too few observations!\n");
        /* *err = E_TOOFEW; */
    }

    if (!*err) {
        pval = fpval(beta.ptr, nbeta, wght.ptr, prob.ptr, cnorm.ptr,
                     tau, T, err);
    }

    return pval;
}

/*
   tau = test statistic
   T = sample size (or 0 for asymptotic)
   niv = # of integrated variables
   itv = 1, 2, 3, or 4 for nc, c, ct, ctt models.

   returns: the computed P-value
*/

extern (C) double mackinnon_pvalue(double tau, int T, int niv, int itv)
{
    double pval = NADBL;
    int err = 0;

    pval = urcval(niv, itv, T, tau, &err);

    return pval;
}

/* Extra: code to obtain approximate finite-sample p-values for
   DF-GLS tests a la Elliott-Rothenberg-Stock, using Sephton's
   response surfaces for the test-down cases, otherwise Cottrell
   and Komashko surfaces.
*/

enum N_ALPHA_C = 25;
enum N_ALPHA_T = 26;

private immutable double[N_ALPHA_C] alpha_c = [
    0.001, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.2, 0.3,
    0.4, 0.5, 0.6, 0.7, 0.9, 0.99
];

private immutable double[N_ALPHA_T] alpha_t = [
    0.001, 0.0025, 0.005, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.15, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.9, 0.99
];

private immutable string[2] bnames = [
    "dfgls-beta-c.bin", "dfgls-beta-t.bin"
];

private immutable string[4] snames = [
    "npc.bin", "npt.bin", "pqc.bin", "pqt.bin"
];

/* In the C this reads a gretl binary matrix file from the data directory of the plugin path; here the matrix is
   made from the embedded copy of that file. */

private GretlMatrixRaw * get_data_matrix(string base, int * err)
{
    foreach (ref em; embedded_matrices) {
        if (em.name == base) {
            GretlMatrixRaw * m = gretl_matrix_alloc(em.rows, em.cols);

            if (m == null) {
                *err = GretlError.alloc;
                return null;
            }
            memcpy(m.val, em.val.ptr, em.rows * em.cols * double.sizeof);
            return m;
        }
    }

    fprintf(stderr, "Couldn't open %s\n", base.ptr);
    *err = GretlError.fopen;

    return null;
}

extern (C) double dfgls_pvalue(double tau, int T, int trend,
                               int kmax, int PQ, int * err)
{
    string fbase;
    GretlMatrixBlock * B;
    GretlMatrixRaw * alpha = null;
    GretlMatrixRaw * beta = null;
    GretlMatrixRaw * C, V, g;
    GretlMatrixRaw * r, y, X;
    const(double) * a;
    const(double) * b;
    double d, dmin = 1.0e6;
    double ci, s2, pval = NADBL;
    int i1, i2, npoints, np2;
    int sephton = (kmax > 0);
    int nreg, ncoeff, nalpha;
    int i, j, k, imin = 0;

    fbase = sephton ? snames[trend + 2*PQ] : bnames[trend];
    beta = get_data_matrix(fbase, err);

    if (!*err && sephton) {
        alpha = get_data_matrix("s_alpha.bin", err);
        if (!*err) {
            a = alpha.val;
            nalpha = alpha.rows;
            npoints = 13;
        }
    } else if (!*err) {
        a = trend ? alpha_t.ptr : alpha_c.ptr;
        nalpha = trend ? N_ALPHA_T : N_ALPHA_C;
        npoints = 5;
    }

    if (*err) {
        gretl_matrix_free(beta);
        return pval;
    }

    ncoeff = beta.rows;
    b = beta.val;
    np2 = npoints / 2;

    /* max # of second-stage regressors */
    nreg = 4;

    /* allocate all storage */
    B = gretl_matrix_block_new(&C, nalpha, 1,
                               &y, npoints, 1,
                               &X, npoints, nreg,
                               &g, nreg, 1,
                               &V, nreg, nreg,
                               &r, ncoeff, 1, null);
    if (B == null) {
        *err = GretlError.alloc;
        goto bailout;
    }

    if (T == 0) {
        /* asymptotic value */
        r.val[0] = 1.0;
        for (j=1; j<ncoeff; j++) {
            r.val[j] = 0;
        }
    } else {
        /* finite sample */
        int npow = ncoeff - 2 - sephton;
        double Tr = 1.0 / T;

        r.val[0] = 1.0;
        r.val[1] = Tr;
        for (j=0; j<npow; j++) {
            r.val[j+2] = pow(Tr, j+2);
        }
        if (sephton) {
            r.val[ncoeff-1] = kmax * Tr;
        }
    }

    /* compute all @nalpha critical values and determine
       which is closest to @tau.
    */
    k = 0;
    for (i=0; i<nalpha; i++) {
        ci = 0;
        for (j=0; j<ncoeff; j++) {
            ci += b[k+j] * r.val[j];
        }
        gretl_matrix_set(C, i, 0, ci);
        if (i == 0) {
            if (tau < 1.25 * C.val[0]) {
                /* extrapolation unlikely to work well for big
                   negative tau: employ conservative fudge
                */
                pval = 1.0e-5;
                goto bailout;
            }
        }
        d = fabs(ci - tau);
        if (d < dmin) {
            dmin = d;
            imin = i;
        }
        k += ncoeff;
    }

    /* select starting point @i1 for critvals range */
    i1 = imin - np2; i2 = imin + np2;
    if (i1 < 0) {
        i1 = 0;
    } else if (i2 >= nalpha) {
        i2 = nalpha - 1;
        i1 = i2 - npoints + 1;
    }

    /* fill out y */
    for (i=0; i<npoints; i++) {
        y.val[i] = normal_cdf_inverse(a[i1 + i]);
    }

    /* fill out X */
    for (j=0; j<nreg; j++) {
        for (i=0; i<npoints; i++) {
            if (j == 0) {
                gretl_matrix_set(X, i, j, 1.0);
            } else {
                ci = C.val[i1 + i];
                if (j == 1) {
                    gretl_matrix_set(X, i, j, ci);
                } else {
                    gretl_matrix_set(X, i, j, pow(ci, j));
                }
            }
        }
    }

    /* we use SVD here just to avoid spamming stderr with
       reports of near-perfect collinearity
    */
    *err = gretl_matrix_SVD_ols(y, X, g, V, null, &s2);

    if (!*err) {
        double t3 = g.val[3] / sqrt(gretl_matrix_get(V, 3, 3));

        if (fabs(t3) < 2) {
            /* drop the cubic term */
            X.cols = 3;
            g.rows = 3;
            *err = gretl_matrix_SVD_ols(y, X, g, null, null, null);
        }
        if (!*err) {
            ci = g.val[0] + g.val[1] * tau;
            for (i=2; i<g.rows; i++) {
                ci += pow(tau, i) * g.val[i];
            }
            pval = normal_cdf(ci);
        }
    }

    if (!*err && (pval < 0 || pval > 1)) {
        /* allow for a little inaccuracy at the extremes */
        pval = pval < 0 ? 0 : 1;
    }

 bailout:

    gretl_matrix_free(beta);
    gretl_matrix_free(alpha);
    gretl_matrix_block_destroy(B);

    return pval;
}
