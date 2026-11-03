/* Tests for plugin/stats_tables.d: the Durbin-Watson table look-up dw_lookup (the only part of stats_tables.c ported so far).
   The compiled-in table is compared with the gzip file it was generated from (gretl/data/dwdata.gz); every (n, k) entry that
   exists for n = 6 .. 200 is compared with an independent parse of that file; the published Savin-White values, the choice of
   the nearest tabulated n and k, and the error cases are checked; gretl_get_DW (pvalues.d) finds the function.
   Build and run with testing/run_tests.sh. */

import core.stdc.math : fabs;
import core.stdc.stdio : _IOLBF, printf, setvbuf, sscanf, stdout;

import std.array : split;
import std.conv : to;
import std.file : read;
import std.string : splitLines, strip;
import std.zlib : HeaderFormat, UnCompress;

import gretl.gretl_matrix : gretl_matrix_free;
import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.plugin.dwdata : dwdata;
import gretl.plugin.stats_tables : dw_lookup, get_IPS_critvals, IPS_tbar_moments, IPS_tbar_rho_moments, qlr_asy_pvalue, qlr_critval_15_05, rank_sum_lookup, stock_yogo_lookup;
import gretl.plugins : get_plugin_function;
import gretl.pvalues : gretl_get_DW;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

bool close(double a, double b, double tol = 1e-12)
{
    return fabs(a - b) <= tol;
}

/* look up, returning the error code and the 4-vector (dl, du, effective n, effective k) */
int look(int n, int k, out double dl, out double du, out int dn, out int dk)
{
    GretlMatrixRaw * m = null;
    const err = dw_lookup(n, k, &m);

    if (err == 0) {
        dl = m.val[0];
        du = m.val[1];
        dn = cast(int) m.val[2];
        dk = cast(int) m.val[3];
        gretl_matrix_free(m);
    }
    return err;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* the compiled-in table is the file */
    auto gz = cast(ubyte[]) read("../gretl/data/dwdata.gz");
    auto un = new UnCompress(HeaderFormat.gzip);
    string text = cast(string) (un.uncompress(gz) ~ un.flush()).idup;

    check(text.length == 71400 && dwdata == text, "the compiled-in table is the decompressed content of dwdata.gz");

    /* every entry for n = 6 .. 200, against an independent parse of the file (row n - 6, entry k) */
    auto rows = text.splitLines();
    bool all = rows.length == 255;
    int compared;
    foreach (n; 6 .. 201) {
        auto fields = rows[n - 6].strip.split;

        foreach (k; 1 .. 21) {
            const el = to!double(fields[2 * (k - 1)]), eu = to!double(fields[2 * (k - 1) + 1]);
            double dl, du;
            int dn, dk;
            const err = look(n, k, dl, du, dn, dk);

            if (el == 0 || eu == 0) {
                all = all && err != 0;
            } else {
                all = all && err == 0 && close(dl, el) && close(du, eu) && dn == n && dk == k;
                compared++;
            }
        }
    }
    printf("  %d entries compared\n", compared);
    check(all && compared > 3000, "all entries for n = 6 .. 200 agree with the table (and the empty ones are errors)");

    /* published 5 percent values (Savin and White) */
    double dl, du;
    int dn, dk;
    check(look(15, 1, dl, du, dn, dk) == 0 && close(dl, 1.0770) && close(du, 1.3605), "n = 15, k = 1: 1.077, 1.3605 (Savin-White 1.08, 1.36)");
    check(look(20, 2, dl, du, dn, dk) == 0 && close(dl, 1.1004, 5e-4) && close(du, 1.5367, 5e-4), "n = 20, k = 2: about 1.100, 1.537");
    check(look(100, 1, dl, du, dn, dk) == 0 && close(dl, 1.6540) && close(du, 1.6944), "n = 100, k = 1: 1.654, 1.6944 (Savin-White 1.65, 1.69)");
    check(look(6, 1, dl, du, dn, dk) == 0 && close(dl, 0.6102) && close(du, 1.4002), "the first entry, n = 6, k = 1");

    /* the nearest tabulated n and k */
    check(look(203, 3, dl, du, dn, dk) == 0 && dn == 200 && dk == 3, "n = 203 is looked up as n = 200");
    check(look(207, 3, dl, du, dn, dk) == 0 && dn == 210, "n = 207 is looked up as n = 210");
    check(look(520, 3, dl, du, dn, dk) == 0 && dn == 500 && look(530, 3, dl, du, dn, dk) == 0 && dn == 550, "n = 520 -> 500 and 530 -> 550 (steps of 50)");
    check(look(2000, 2, dl, du, dn, dk) == 0 && dn == 2000 && look(2500, 2, dl, du, dn, dk) == 0 && dn == 2000, "n of 2000 or more is looked up as 2000");
    check(look(100, 25, dl, du, dn, dk) == 0 && dk == 20, "k above 20 is looked up as 20");
    double d1, u1, d2, u2;
    int n1, k1, n2, k2;
    look(200, 2, d1, u1, n1, k1);
    look(500, 2, d2, u2, n2, k2);
    check(d2 > d1 && (u2 - d2) < (u1 - d1) && look(1000, 2, dl, du, dn, dk) == 0 && dl > d2 && (du - dl) < (u2 - d2), "the bounds rise towards 2 and the gap between them shrinks as n grows");

    /* errors */
    GretlMatrixRaw * m = null;
    check(dw_lookup(5, 1, &m) != 0 && m is null, "n below 6 is an error");
    check(look(6, 5, dl, du, dn, dk) != 0, "an entry that is not in the table (n = 6, k = 5) is an error");
    check(look(10, 20, dl, du, dn, dk) != 0, "another empty entry (n = 10, k = 20) is an error");
    check((look(50, 0, dl, du, dn, dk) || true) && (look(50, -3, dl, du, dn, dk) || true), "k below 1 does not crash (bugs.md #113: the C does not check it)");

    /* the library finds it */
    int err = 0;
    auto v = gretl_get_DW(15, 1, &err);
    check(get_plugin_function("dw_lookup") !is null && err == 0 && v !is null && v.rows == 1 && v.cols == 4 && close(v.val[0], 1.0770) && close(v.val[1], 1.3605), "gretl_get_DW (pvalues.d) finds dw_lookup and returns the 4-vector");
    gretl_matrix_free(v);

    /* ---- rank_sum_lookup: critical values of the Wilcoxon rank-sum test for small samples ---- */
    check(get_plugin_function("rank_sum_lookup") !is null, "the library finds rank_sum_lookup");
    static string lookupText(int na, int nb)
    {
        import gretl.gretl_prn : gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
        import gretl.gretl_prnh : PrnType;
        import std.string : fromStringz;

        int perr;
        auto prn = gretl_print_new(PrnType.buffer, &perr);

        rank_sum_lookup(na, nb, prn);
        auto text = fromStringz(gretl_print_get_buffer(prn)).idup;

        gretl_print_destroy(prn);
        return text;
    }
    // the C prints the upper tail of the first row from columns 1 and 2 of the table instead of 0 and 1: "10% 25, 5% -1" for 23 and 25 (bugs.md #128); the port keeps that
    check(lookupText(4, 4) == "\nCritical values:\n  lower tail:  5% 11, 10% 13\n  upper tail: 10% 25,  5% -1\n", "(4, 4): only the 5% and 10% values exist (the upper tail is printed from the wrong columns, as in the C)");
    check(lookupText(4, 5) == "\nCritical values:\n  lower tail:  1% 10,  5% 12, 10% 14\n  upper tail: 10% 26,  5% 28,  1% 30\n", "(4, 5): all three levels");
    check(lookupText(9, 12) == "\nCritical values:\n  lower tail:  1% 66,  5% 75, 10% 80\n  upper tail: 10% 118,  5% 123,  1% 132\n", "(9, 12): the last row");
    check(lookupText(3, 5).length == 0 && lookupText(5, 4).length == 0 && lookupText(10, 10).length == 0 && lookupText(4, 13).length == 0, "sizes outside the table print nothing");
    // each critical value pair is symmetric: lower + upper = nA (nA + nB + 1), since the distribution of the rank sum is symmetric about nA (nA + nB + 1) / 2
    bool symmetric = true;
    int tableRows = 0;
    foreach (na; 4 .. 10) {
        foreach (nb; na .. 13) {
            const t = lookupText(na, nb);
            int l1 = -1, l5, l10, u10, u5, uu1 = -1;
            const total = na * (na + nb + 1);

            if (na == 4 && nb == 4) {
                sscanf(t.ptr, "\nCritical values:\n  lower tail:  5%% %d, 10%% %d\n", &l5, &l10);
                symmetric = symmetric && l5 == 11 && l10 == 13 && total - l5 == 25 && total - l10 == 23;      // the upper values by symmetry: 23 and 25
            } else {
                sscanf(t.ptr, "\nCritical values:\n  lower tail:  1%% %d,  5%% %d, 10%% %d\n  upper tail: 10%% %d,  5%% %d,  1%% %d\n", &l1, &l5, &l10, &u10, &u5, &uu1);
                if (l5 + u5 != total || l10 + u10 != total || l1 + uu1 != total) {
                    printf("  asymmetric: nA %d nB %d: %d+%d %d+%d %d+%d vs %d\n", na, nb, l1, uu1, l5, u5, l10, u10, total);
                    symmetric = false;
                }
            }
            tableRows++;
        }
    }
    check(tableRows == 39 && symmetric, "all 39 rows: lower + upper critical values add to nA (nA + nB + 1) at each level");

    /* ---- stock_yogo_lookup: Stock and Yogo (2003) critical values for weak instruments ---- */
    check(get_plugin_function("stock_yogo_lookup") !is null, "the library finds stock_yogo_lookup");
    static bool row(int n, int K2, int which, double[4] crit, double[4] vals)
    {
        auto m = stock_yogo_lookup(n, K2, which);

        if (m is null || m.rows != 2 || m.cols != 4) return false;
        bool ok = true;

        foreach (i; 0 .. 4) {
            ok = ok && close(m.val[2 * i], crit[i], 1e-12) && close(m.val[2 * i + 1], vals[i], 1e-12);      // column-major: (0, i), (1, i)
        }
        gretl_matrix_free(m);
        return ok;
    }
    immutable double[4] bvals = [0.05, 0.10, 0.20, 0.30], rvals = [0.10, 0.15, 0.20, 0.25];
    // Table 1, relative bias of TSLS: the published first and last rows of each column block
    check(row(1, 3, 1, bvals, [13.91, 9.08, 6.46, 5.39]) && row(1, 30, 1, bvals, [21.42, 11.32, 6.09, 4.29]), "Table 1 (TSLS bias), one endogenous regressor: K2 = 3 and 30");
    check(row(2, 4, 1, bvals, [11.04, 7.56, 5.57, 4.73]) && row(2, 30, 1, bvals, [20.86, 11.05, 5.99, 4.23]), "Table 1, two endogenous regressors: K2 = 4 and 30");
    check(row(3, 5, 1, bvals, [9.53, 6.61, 4.99, 4.30]) && row(3, 30, 1, bvals, [20.27, 10.77, 5.87, 4.17]), "Table 1, three endogenous regressors: K2 = 5 and 30");
    // Table 2, size of TSLS tests at a nominal 5%
    check(row(1, 1, 2, rvals, [16.38, 8.96, 6.66, 5.53]) && row(1, 30, 2, rvals, [86.17, 44.78, 30.72, 23.65]) && row(2, 2, 2, rvals, [7.03, 4.58, 3.95, 3.63]) && row(2, 30, 2, rvals, [63.51, 33.61, 23.51, 18.35]),
          "Table 2 (TSLS size): K2 = 1, 2 and 30");
    // Table 3, size of LIML tests
    check(row(1, 1, 3, rvals, [16.38, 8.96, 6.66, 5.53]) && row(1, 2, 3, rvals, [8.68, 5.33, 4.42, 3.92]) && row(1, 30, 3, rvals, [3.88, 2.18, 1.89, 1.75]) && row(2, 30, 3, rvals, [4.12, 2.39, 1.95, 1.75]),
          "Table 3 (LIML size): K2 = 1, 2 and 30");
    // what is out of range gives NULL
    check(stock_yogo_lookup(0, 5, 1) is null && stock_yogo_lookup(4, 10, 1) is null && stock_yogo_lookup(3, 10, 2) is null && stock_yogo_lookup(3, 10, 3) is null, "n below 1 or above 3 (bias) / 2 (size) gives NULL");
    check(stock_yogo_lookup(1, 2, 1) is null && stock_yogo_lookup(1, 31, 1) is null && stock_yogo_lookup(1, 0, 2) is null && stock_yogo_lookup(2, 1, 3) is null, "K2 below the minimum (3 for bias, 1 for size), above 30, or below n gives NULL");
    // structure: in every valid case the critical values fall as the tolerated distortion rises, and TSLS size values rise with K2 (n = 1)
    static void structure(ref bool ok, ref int cases)
    {
        foreach (which; 1 .. 4) {
            foreach (n; 1 .. (which == 1 ? 4 : 3)) {
                foreach (K2; (which == 1 ? 3 : 1) .. 31) {
                    auto mm = stock_yogo_lookup(n, K2, which);

                    if (mm is null) continue;
                    double[4] cc, vv;

                    foreach (i; 0 .. 4) { cc[i] = mm.val[2 * i]; vv[i] = mm.val[2 * i + 1]; }
                    gretl_matrix_free(mm);
                    if (vv[0] == 0) continue;                    // an empty cell of Table 1 (see below)
                    cases++;
                    ok = ok && vv[0] > vv[1] && vv[1] > vv[2] && vv[2] > vv[3] && cc[0] < cc[1] && cc[1] < cc[2] && cc[2] < cc[3];
                }
            }
        }
    }
    bool structureOk = true;
    int cases = 0;

    structure(structureOk, cases);
    // 28 + 27 + 26 rows of Table 1, 30 + 29 of Table 2, 30 + 29 of Table 3
    check(structureOk && cases == 28 + 27 + 26 + 30 + 29 + 30 + 29, "all 199 valid rows: the critical values fall as the tolerated bias or size rises");
    double prev = 0;
    bool rising = true;
    foreach (K2; 1 .. 31) {
        auto mm = stock_yogo_lookup(1, K2, 2);

        rising = rising && mm.val[1] > prev;
        prev = mm.val[1];
        gretl_matrix_free(mm);
    }
    check(rising, "TSLS size critical values for one regressor rise with K2");
    // the cells of Table 1 that do not exist (two regressors with K2 = 3, three with K2 = 3 or 4) give a matrix of zeros, not NULL (bugs.md #129)
    auto zz = stock_yogo_lookup(2, 3, 1);
    check(zz !is null && zz.val[1] == 0 && zz.val[3] == 0 && zz.val[5] == 0 && zz.val[7] == 0 && stock_yogo_lookup(3, 4, 1) !is null, "empty cells of Table 1 give zeros (as in the C)");
    gretl_matrix_free(zz);

    {
        import gretl.libgretlh : Gretlopt;
        import gretl.modelprint : get_stock_yogo_critvals;
        GretlMatrixRaw * p1 = null;
        GretlMatrixRaw * p2 = null;

        get_stock_yogo_critvals(1, 4, Gretlopt.OPT_NONE, &p1, &p2);
        check(p1 !is null && p2 !is null && close(p1.val[1], 16.85) && close(p2.val[1], 24.58), "get_stock_yogo_critvals (modelprint.d) finds the plugin function: TSLS bias and size rows");
        gretl_matrix_free(p1);
        gretl_matrix_free(p2);
        p1 = p2 = null;
        get_stock_yogo_critvals(1, 4, Gretlopt.OPT_L, &p1, &p2);
        check(p1 is null && p2 !is null && close(p2.val[1], 5.44) && close(p2.val[7], 2.98), "with the LIML option only the LIML size row is looked up");
        gretl_matrix_free(p2);
    }

    /* ---- qlr_asy_pvalue: Hansen's approximate asymptotic p-values of the sup-Wald (QLR) test ---- */
    {
        import gretl.pvalues : chisq_cdf_comp;
        import core.stdc.math : sqrt;

        check(get_plugin_function("qlr_asy_pvalue") !is null, "the library finds qlr_asy_pvalue");
        // lambda = 1 is tau = 0.5, the symmetric case where the statistic is the ordinary chi-square
        bool chi = true;
        foreach (df; [1, 2, 5, 12, 40]) foreach (x; [0.5, 3.0, 9.0, 20.0]) chi = chi && close(qlr_asy_pvalue(x, df, 1.0), chisq_cdf_comp(df, x), 1e-12);
        check(chi, "lambda = 1 (tau = 0.5): the chi-square(df) tail probability");
        // lambda below 1 is tau itself; 15% trimming is lambda = (0.85 / 0.15)^2, tau = 1 / (1 + sqrt(lambda)) = 0.15
        const lam15 = (0.85 * 0.85) / (0.15 * 0.15);
        check(close(qlr_asy_pvalue(9.0, 3, lam15), qlr_asy_pvalue(9.0, 3, 0.15), 1e-12) && close(qlr_asy_pvalue(9.0, 3, 4.0), qlr_asy_pvalue(9.0, 3, 1.0 / 3), 1e-12), "lambda above 1 is converted to tau = 1 / (1 + sqrt(lambda))");
        // Andrews (1993) 15% trimming: the 5% critical value is about 8.7 to 8.9 for one restriction, the 1% about 12.2 to 12.4, the 10% about 7.1 to 7.2
        const p5 = qlr_asy_pvalue(8.8, 1, lam15), p1 = qlr_asy_pvalue(12.3, 1, lam15), p10 = qlr_asy_pvalue(7.15, 1, lam15);
        check(p5 > 0.04 && p5 < 0.06 && p1 > 0.005 && p1 < 0.015 && p10 > 0.08 && p10 < 0.12, "15% trimming, one restriction: p(8.8) near 5%, p(12.3) near 1%, p(7.15) near 10%");
        check(qlr_asy_pvalue(11.8, 2, lam15) > 0.04 && qlr_asy_pvalue(11.8, 2, lam15) < 0.06 && qlr_asy_pvalue(13.9, 3, lam15) > 0.04 && qlr_asy_pvalue(13.9, 3, lam15) < 0.06, "15% trimming: p(11.8) near 5% for two restrictions, p(13.9) for three");
        // the sup statistic has a larger tail probability than the chi-square, and the more so the more it is trimmed less (smaller tau)
        bool ordered = true;
        foreach (df; [1, 3, 8, 20]) {
            foreach (x; [6.0, 10.0, 15.0, 25.0]) {
                double last = chisq_cdf_comp(df, x);

                foreach (tau; [0.49, 0.45, 0.35, 0.25, 0.15, 0.05, 0.01]) {
                    const pv = qlr_asy_pvalue(x, df, tau);

                    ordered = ordered && pv >= last * 0.9 && pv >= 0 && pv <= 1;       // the approximation is not exactly monotone far in the tail (df = 1, x = 25: 2.6e-5 at tau = 0.01 against 2.8e-5 at 0.05)
                    last = pv;
                }
            }
        }
        check(ordered, "p-values are probabilities and rise as tau falls from 0.49 to 0.01 (to within 10% in the far tail), starting from the chi-square value");
        // decreasing in the statistic
        bool falling = true;
        foreach (tau; [0.4, 0.15, 0.02]) {
            double last = 2;

            foreach (x; [1.0, 3.0, 6.0, 9.0, 12.0, 16.0, 20.0, 30.0]) {
                const pv = qlr_asy_pvalue(x, 4, tau);

                falling = falling && pv <= last && (last >= 1 || pv < last);       // 1 for small statistics (x near 0), strictly falling after
                last = pv;
            }
        }
        check(falling, "the p-value falls as the statistic rises (it is 1 for very small statistics)");
        // linear interpolation between the grid points tau = 0.49, 0.47, ...: the midpoint is the mean; the ends are continuous
        const mid = qlr_asy_pvalue(10.0, 2, 0.48), hi = qlr_asy_pvalue(10.0, 2, 0.49), lo = qlr_asy_pvalue(10.0, 2, 0.47);
        check(close(mid, (hi + lo) / 2, 1e-9), "between grid points the p-value is interpolated linearly in tau");
        check(close(qlr_asy_pvalue(10.0, 2, 0.4999999), chisq_cdf_comp(2, 10.0), 1e-5) && close(qlr_asy_pvalue(10.0, 2, 0.4900001), hi, 1e-5) && close(qlr_asy_pvalue(10.0, 2, 0.0100001), qlr_asy_pvalue(10.0, 2, 0.01), 1e-4), "continuous at tau = 0.5, 0.49 and near 0.01");
        check(close(qlr_asy_pvalue(10.0, 41, lam15), qlr_asy_pvalue(10.0, 40, lam15), 1e-12) && close(qlr_asy_pvalue(10.0, 100, 0.3), qlr_asy_pvalue(10.0, 40, 0.3), 1e-12), "more than 40 degrees of freedom are treated as 40");
    }

    /* ---- qlr_critval_15_05: the 5% critical value of the QLR test with 15% trimming ---- */
    {
        import std.math : isNaN;
        const lam15 = (0.85 * 0.85) / (0.15 * 0.15);

        check(get_plugin_function("qlr_critval_15_05") !is null, "the library finds qlr_critval_15_05");
        bool found = true, close5 = true, growing = true;
        double lastc = 0;

        foreach (df; 1 .. 41) {
            const c = qlr_critval_15_05(df);

            found = found && !isNaN(c);
            if (isNaN(c)) continue;
            const pv = qlr_asy_pvalue(c, df, lam15);

            close5 = close5 && pv >= 0.049 && pv <= 0.051;            // the stopping rule of the bisection
            growing = growing && c > lastc;
            lastc = c;
        }
        check(found && close5, "df = 1 .. 40: a critical value is found and its p-value is within 0.049 and 0.051");
        check(growing, "the critical value rises with the degrees of freedom");
        // Andrews (1993, Table 1, 15% trimming, 5%): 8.68, 11.79, 13.89, 16.17, 18.11 for one to five restrictions
        immutable double[5] andrews = [8.68, 11.79, 13.89, 16.17, 18.11];
        bool near = true;

        foreach (i; 0 .. 5) near = near && fabs(qlr_critval_15_05(i + 1) - andrews[i]) < 0.35;
        check(near, "within 0.35 of Andrews' published 5% critical values for one to five restrictions");
        const c1 = qlr_critval_15_05(1);
        check(qlr_asy_pvalue(c1 - 0.5, 1, lam15) > 0.05 && qlr_asy_pvalue(c1 + 0.5, 1, lam15) < 0.05, "the value separates p above 5% (below it) from p below 5% (above it)");
        check(qlr_critval_15_05(41) == qlr_critval_15_05(40) && qlr_critval_15_05(100) == qlr_critval_15_05(40), "more than 40 degrees of freedom are treated as 40");
    }

    /* ---- get_IPS_critvals: critical values of the Im-Pesaran-Shin tbar statistic ---- */
    {
        import gretl.gretl_errors : gretl_errmsg_get, gretl_error_clear;
        import std.string : fromStringz;

        check(get_plugin_function("get_IPS_critvals") !is null, "the library finds get_IPS_critvals");
        static double[3] ips(int N, int T, int trend, ref int err)
        {
            double[3] c = double.nan;

            err = get_IPS_critvals(N, T, trend, c.ptr);
            return c;
        }
        static bool same(double[3] c, double a, double b, double d)
        {
            return fabs(c[0] - a) < 1e-12 && fabs(c[1] - b) < 1e-12 && fabs(c[2] - d) < 1e-12;      // 10%, 5%, 1%
        }
        int e;

        // the values on the grid are the table entries (IPS 2003, Table 2)
        check(same(ips(5, 5, 0, e), -2.38, -2.76, -3.79) && e == 0 && same(ips(10, 20, 0, e), -1.89, -1.99, -2.21) && same(ips(25, 30, 0, e), -1.75, -1.82, -1.94), "grid points without trend: (5, 5), (10, 20), (25, 30)");
        check(same(ips(100, 100, 0, e), -1.64, -1.67, -1.73) && same(ips(100, 100, 1, e), -2.28, -2.31, -2.36) && same(ips(15, 40, 1, e), -2.44, -2.52, -2.67), "grid points: (100, 100) with and without trend, (15, 40) with trend");
        // interpolation: one dimension, weights 1 / distance
        double[3] c = ips(12, 10, 0, e);
        const w0 = 1.0 / 2, w1 = 1.0 / 3;
        check(fabs(c[1] - (w0 * -2.06 + w1 * -1.95) / (w0 + w1)) < 1e-12, "N between grid values, T on the grid: inverse-distance weights of the two neighbours (N = 12 between 10 and 15)");
        c = ips(10, 12, 0, e);
        check(fabs(c[1] - (w0 * -2.06 + w1 * -2.02) / (w0 + w1)) < 1e-12, "T between grid values, N on the grid (T = 12 between 10 and 15)");
        // both between: four corners with Euclidean-distance weights
        import core.stdc.math : sqrt;
        c = ips(12, 12, 0, e);
        const q0 = 1.0 / sqrt(8.0), q1 = 1.0 / sqrt(13.0), q2 = 1.0 / sqrt(18.0), q3 = 1.0 / sqrt(13.0);
        check(fabs(c[1] - (q0 * -2.06 + q1 * -2.02 + q2 * -1.92 + q3 * -1.95) / (q0 + q1 + q2 + q3)) < 1e-12, "N and T both between grid values: four corners, Euclidean-distance weights");
        check(c[1] > -2.06 && c[1] < -1.92 && c[2] < c[1] && c[1] < c[0], "the interpolated value lies among the corners and the levels are ordered");
        // beyond the published sizes
        check(same(ips(500, 500, 0, e), -1.64, -1.67, -1.73) && same(ips(100, 300, 1, e), -2.28, -2.31, -2.36) && same(ips(300, 100, 0, e), -1.64, -1.67, -1.73), "N and T both at or above 100: the (100, 100) values");
        check(fabs(ips(150, 20, 0, e)[1] - -1.67) < 1e-12 && fabs(ips(150, 20, 0, e)[0] - -1.64) < 1e-12 && e == 0, "N above 100 with a short T (20): the values of the N = 100 row");
        check(fabs(ips(5, 500, 0, e)[1] - -2.15) < 1e-12 && e == 0, "T above 100 with N = 5: the T = 100 column");
        // errors
        double[3] keep = [1.0, 2.0, 3.0];
        gretl_error_clear();
        e = get_IPS_critvals(4, 20, 0, keep.ptr);
        check(e != 0 && keep[0] == 1.0 && fromStringz(gretl_errmsg_get()) == "IPS test: too few groups", "fewer than 5 groups is an error and leaves the output alone");
        gretl_error_clear();
        e = get_IPS_critvals(10, 4, 0, keep.ptr);
        check(e != 0 && keep[1] == 2.0 && fromStringz(gretl_errmsg_get()) == "IPS test: time series too short", "T below 5 is an error");
        // structure of the whole table: at each point 1% < 5% < 10% < 0, a trend makes the values more negative, and they rise with N and T
        bool ordered = true, trendLower = true;
        foreach (N; [5, 6, 7, 9, 10, 12, 15, 17, 20, 22, 25, 30, 50, 75, 100, 130]) {
            foreach (T; [5, 7, 10, 12, 15, 17, 20, 22, 25, 28, 30, 35, 40, 45, 50, 55, 60, 65, 70, 85, 100, 140]) {
                auto a = ips(N, T, 0, e), b = ips(N, T, 1, e);

                ordered = ordered && e == 0 && a[2] < a[1] && a[1] < a[0] && a[0] < 0 && b[2] < b[1] && b[1] < b[0] && b[0] < 0;
                trendLower = trendLower && b[0] < a[0] && b[1] < a[1] && b[2] < a[2];
            }
        }
        check(ordered && trendLower, "over a grid of N and T (interpolated too): 1% < 5% < 10% < 0, and a trend gives lower values");
        immutable int[8] gridN = [5, 7, 10, 15, 20, 25, 50, 100];
        immutable int[11] gridT = [5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 100];
        bool rises = true;
        int anomalies = 0;

        foreach (trend; 0 .. 2) {
            foreach (lev; 0 .. 3) {
                foreach (i; 0 .. 8) {
                    foreach (j; 0 .. 11) {
                        const here = ips(gridN[i], gridT[j], trend, e)[lev];

                        if (j < 10) {
                            const next = ips(gridN[i], gridT[j + 1], trend, e)[lev];

                            if (next < here - 0.011) { anomalies++; rises = false; }
                        }
                        if (i < 7) {
                            const up = ips(gridN[i + 1], gridT[j], trend, e)[lev];

                            if (up < here - 0.011) { anomalies++; rises = false; }
                        }
                    }
                }
            }
        }
        // the one entry that was out of line in the C (the 1% value with trend for N = 7, T = 50: -2.99 after -2.90 and before -2.88) is -2.89 in the port (bugs.md #131)
        check(rises && anomalies == 0 && fabs(ips(7, 50, 1, e)[2] - -2.89) < 1e-12 && fabs(ips(7, 40, 1, e)[2] - -2.90) < 1e-12 && fabs(ips(7, 60, 1, e)[2] - -2.88) < 1e-12,
              "the table values rise with N and T to within rounding everywhere, N = 7, T = 50 (1%, trend) being -2.89 as corrected (bugs.md #131)");
        // and the interpolation near it uses the corrected value: N = 7, T = 45 is between -2.89 (T = 50) and -2.90 (T = 40)
        const mid = ips(7, 45, 1, e)[2];
        check(mid > -2.90 && mid < -2.89 + 1e-12 && fabs(mid - (-2.895)) < 1e-9, "interpolation at N = 7, T = 45 (trend, 1%) uses the corrected entry: -2.895");
    }

    /* ---- IPS_tbar_moments: moments of the t-bar statistic (IPS 2003, Table 1) ---- */
    {
        import std.math : isNaN;

        check(get_plugin_function("IPS_tbar_moments") !is null, "the library finds IPS_tbar_moments");
        static bool mom(int T, double E, double V)
        {
            double e = 0, v = 0;
            const err = IPS_tbar_moments(T, &e, &v);

            return err == 0 && fabs(e - E) < 1e-12 && fabs(v - V) < 1e-12;
        }
        // the values on the grid are the table entries
        check(mom(6, -1.520, 1.745) && mom(7, -1.514, 1.414) && mom(10, -1.504, 1.069) && mom(15, -1.514, 0.923) && mom(25, -1.520, 0.809) && mom(50, -1.527, 0.760), "grid values T = 6, 7, 10, 15, 25, 50");
        check(mom(100, -1.532, 0.735) && mom(500, -1.531, 0.715) && mom(1000, -1.529, 0.707) && mom(5000, -1.529, 0.707) && mom(100000, -1.529, 0.707), "grid values T = 100, 500, 1000, and 1000 or more is the T = 1000 row");
        // interpolation with inverse-distance weights
        const a1 = 1.0 / 2, a2 = 1.0 / 3;                           // T = 12 between 10 and 15
        check(mom(12, (a1 * -1.504 + a2 * -1.514) / (a1 + a2), (a1 * 1.069 + a2 * 0.923) / (a1 + a2)), "T = 12 between 10 and 15: inverse-distance weights 1/2 and 1/3");
        const b1 = 1.0 / 200, b2 = 1.0 / 300;                       // T = 700 between 500 and 1000
        check(mom(700, (b1 * -1.531 + b2 * -1.529) / (b1 + b2), (b1 * 0.715 + b2 * 0.707) / (b1 + b2)), "T = 700 between 500 and 1000");
        check(mom(8, -1.501, 1.228) && mom(9, -1.501, 1.132), "T = 8 and 9 on the grid");
        // too short a series
        double me = 1, mv = 1;
        check(IPS_tbar_moments(5, &me, &mv) != 0 && isNaN(me) && isNaN(mv) && IPS_tbar_moments(0, &me, &mv) != 0 && IPS_tbar_moments(-3, &me, &mv) != 0, "T below 6 is an error with NaN moments");
        // structure over T = 6 .. 1200: the variance falls with T to about 0.7, the mean stays near -1.52, interpolated values lie between their neighbours
        bool varFalls = true, meanRange = true;
        double lastV = 1e9;

        foreach (T; 6 .. 1201) {
            double ee, vv;

            IPS_tbar_moments(T, &ee, &vv);
            varFalls = varFalls && vv <= lastV + 1e-12 && vv > 0.7 - 1e-12 && vv < 1.75 + 1e-12;
            meanRange = meanRange && ee > -1.54 && ee < -1.49;
            lastV = vv;
        }
        check(varFalls && meanRange, "T = 6 .. 1200: the variance falls monotonically from 1.745 to 0.707, the mean stays between -1.54 and -1.49");
        double e1, v1, e2, v2, em, vm;
        IPS_tbar_moments(40, &e1, &v1);
        IPS_tbar_moments(50, &e2, &v2);
        IPS_tbar_moments(43, &em, &vm);
        check(vm < v1 && vm > v2 && em >= (e1 < e2 ? e1 : e2) && em <= (e1 > e2 ? e1 : e2) && vm - v2 > v1 - vm - 0.1, "an interpolated value (T = 43) lies between those at 40 and 50, nearer the one at 40");
    }

    /* ---- IPS_tbar_rho_moments: moments of t_T(rho, 0) for lag order p (IPS 2003, Table 3) ---- */
    {
        import std.math : isNaN;

        check(get_plugin_function("IPS_tbar_rho_moments") !is null, "the library finds IPS_tbar_rho_moments");
        static int rho(int p, int T, int trend, ref double E, ref double V)
        {
            return IPS_tbar_rho_moments(p, T, trend, &E, &V);
        }
        static bool cell(int p, int T, int trend, double E, double V)
        {
            double e, v;

            return rho(p, T, trend, e, v) == 0 && fabs(e - E) < 1e-12 && fabs(v - V) < 1e-12;
        }
        // table entries: without trend
        check(cell(0, 10, 0, -1.504, 1.069) && cell(4, 10, 0, -1.171, 2.080) && cell(8, 25, 0, -1.212, 1.217) && cell(5, 20, 0, -1.313, 1.171) && cell(3, 60, 0, -1.486, 0.802) && cell(8, 100, 0, -1.456, 0.818), "table entries without trend");
        // with trend
        check(cell(0, 10, 1, -2.166, 1.132) && cell(4, 10, 1, -1.750, 3.947) && cell(8, 25, 1, -1.761, 1.208) && cell(2, 40, 1, -2.120, 0.702) && cell(8, 100, 1, -2.088, 0.670), "table entries with trend");
        check(cell(2, 150, 0, -1.514, 0.754) && cell(6, 5000, 1, -2.113, 0.650), "T of 100 or more is the T = 100 row");
        // for p = 0 the no-trend column is Table 1 of IPS_tbar_moments, an independent table
        bool agree = true;

        foreach (T; [10, 15, 20, 25, 30, 40, 50, 100]) {
            double e1, v1, e2, v2;

            IPS_tbar_moments(T, &e1, &v1);
            rho(0, T, 0, e2, v2);
            agree = agree && fabs(e1 - e2) < 1e-12 && fabs(v1 - v2) < 1e-12;
        }
        check(agree, "the p = 0 column without trend equals IPS_tbar_moments at T = 10, 15, 20, 25, 30, 40, 50, 100");
        // interpolation by inverse distances
        double ie, iv;
        rho(2, 12, 0, ie, iv);
        const a1 = 1.0 / 2, a2 = 1.0 / 3;                                // T = 12 between 10 and 15
        check(fabs(ie - (a1 * -1.319 + a2 * -1.387) / (a1 + a2)) < 1e-12 && fabs(iv - (a1 * 1.421 + a2 * 1.078) / (a1 + a2)) < 1e-12, "T = 12, p = 2: inverse-distance weights 1/2 and 1/3");
        rho(7, 85, 1, ie, iv);
        const b1 = 1.0 / 15, b2 = 1.0 / 15;                              // T = 85 between 70 and 100: equal distances
        check(fabs(ie - (b1 * -2.079 + b2 * -2.112) / (b1 + b2)) < 1e-12 && fabs(iv - (b1 * 0.713 + b2 * 0.661) / (b1 + b2)) < 1e-12, "T = 85, p = 7, with trend: the mean of the T = 70 and T = 100 entries");
        // unhandled cases: the zero cells, T below 10, p above 8
        double me = 1, mv = 1;
        check(rho(5, 10, 0, me, mv) != 0 && isNaN(me) && isNaN(mv) && rho(8, 20, 0, me, mv) != 0 && rho(6, 20, 1, me, mv) != 0 && rho(5, 15, 1, me, mv) != 0, "the zero cells of the table are errors with NaN moments");
        check(rho(5, 12, 0, me, mv) != 0 && rho(6, 22, 0, me, mv) != 0 && rho(5, 22, 0, me, mv) == 0, "interpolating from a zero cell is an error; from the T = 20 row at p = 5 it is not");
        check(rho(0, 9, 0, me, mv) != 0 && isNaN(me) && rho(0, 0, 1, me, mv) != 0 && rho(9, 50, 0, me, mv) != 0 && isNaN(mv) && rho(100, 50, 1, me, mv) != 0, "T below 10 or p above 8 is an error with NaN moments");
        // structure over all valid cells: the variance falls with T, rises with p, a trend lowers the mean
        bool varT = true, varP = true, trendLow = true;
        immutable int[10] gridT = [10, 15, 20, 25, 30, 40, 50, 60, 70, 100];

        foreach (trend; 0 .. 2) {
            foreach (p; 0 .. 9) {
                foreach (i; 0 .. 10) {
                    double e0, v0, e1, v1;
                    const ok0 = rho(p, gridT[i], trend, e0, v0) == 0;

                    if (!ok0) continue;
                    if (i < 9 && rho(p, gridT[i + 1], trend, e1, v1) == 0) varT = varT && v1 < v0;
                    if (p < 8 && rho(p + 1, gridT[i], trend, e1, v1) == 0) varP = varP && v1 > v0;
                    double en, vn;
                    if (trend) {
                        rho(p, gridT[i], 0, en, vn);
                        trendLow = trendLow && e0 < en;
                    }
                }
            }
        }
        check(varT && varP && trendLow, "over all valid cells: the variance falls with T and rises with p, and a trend gives a lower mean");
    }

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
