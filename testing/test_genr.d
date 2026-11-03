/* End-to-end tests of 'genr' through generate(): the lexer (genlex.d), the parser (gensyntax.d), the evaluator (geneval0..5.d) and the saving of
   results (genmain.d), on scalars, series, matrices and strings. Every result is compared with an independent computation in D.

   This is the first time the ported evaluator is driven by real expressions (before the parser existed it could only be reached directly), so the tests
   are deliberately broad: arithmetic and precedence, comparison and logical operators, the ternary operator, built-in functions, series with lags, leads,
   differences and cumulation, series statistics, matrix literals and operators, slicing, concatenation, determinants and inverses, string operations,
   and the errors for malformed expressions.

   gretl_get_object_by_name (objstack.c, not ported yet) is given a stand-in that knows no saved objects.

   Build and run with testing/run_tests.sh. */

import core.stdc.math : exp, fabs, log, sin, sqrt, floor, ceil, pow;
import core.stdc.stdio : printf, setvbuf, stdout, _IOLBF;
import core.stdc.string : strcmp, strcpy;

import std.math : PI, isNaN;
import std.string : toStringz;

import gretl.dataset : create_new_dataset, dataset_set_time_series;
import gretl.genmain : current_series_index, generate;
import gretl.gretl_errors : gretl_error_clear;
import gretl.gretl_matrixh : GretlMatrixRaw, gretl_matrix_get;
import gretl.libgretlh : DatasetRaw, GretlType, Gretlopt;
import gretl.usermat : get_matrix_by_name;
import gretl.uservar : get_string_by_name, gretl_scalar_get_value;


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
    if (isNaN(a) && isNaN(b)) return true;
    return fabs(a - b) <= tol * (1 + fabs(b));
}

enum T = 12;
DatasetRaw * dset;
double[T] yv, xv;

int run(string src, GretlType t)
{
    gretl_error_clear();
    return generate(src.toStringz, dset, t, Gretlopt.OPT_NONE, null);
}

/* the value of a scalar expression (NaN on error) */
double sc(string expr)
{
    if (run("tmp_s = " ~ expr, GretlType.double_) != 0) return double.nan;
    int e;

    return gretl_scalar_get_value("tmp_s", &e);
}

/* the values of a series expression (empty on error) */
double[] se(string expr)
{
    if (run("tmp_x = " ~ expr, GretlType.series) != 0) return [];
    const v = current_series_index(dset, "tmp_x");

    return dset.Z[v][0 .. T].dup;
}

/* the values of a matrix expression, as rows (empty on error) */
double[][] ma(string expr)
{
    if (run("tmp_m = " ~ expr, GretlType.matrix) != 0) return [];
    auto m = get_matrix_by_name("tmp_m");
    double[][] r;

    foreach (i; 0 .. m.rows) {
        double[] row;

        foreach (j; 0 .. m.cols) row ~= gretl_matrix_get(m, i, j);
        r ~= row;
    }
    return r;
}

bool sameSeries(double[] a, double[] b, double tol = 1e-12)
{
    if (a.length != b.length) return false;
    foreach (i; 0 .. a.length) if (!close(a[i], b[i], tol)) return false;
    return true;
}

bool sameMatrix(double[][] a, double[][] b)
{
    if (a.length != b.length) return false;
    foreach (i; 0 .. a.length) if (!sameSeries(a[i], b[i])) return false;
    return true;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    dset = create_new_dataset(3, T, 0);
    strcpy(dset.varname[1], "y");
    strcpy(dset.varname[2], "x");
    foreach (t; 0 .. T) {
        yv[t] = 3.0 + 0.5 * t + ((t % 3) - 1) * 0.7;
        xv[t] = (t % 4) + 1.0;
        dset.Z[1][t] = yv[t];
        dset.Z[2][t] = xv[t];
    }
    dataset_set_time_series(dset, 1, 1900, 1);

    /* ---- scalars ---- */
    check(sc("2 + 3 * 4") == 14 && sc("(2 + 3) * 4") == 20 && sc("10 - 4 - 3") == 3 && sc("2 * 3 ^ 2") == 18, "arithmetic and precedence");
    check(sc("7 % 3") == 1 && sc("7.5 / 2.5") == 3, "modulus and division");
    check(close(sc("sqrt(2) * sqrt(2)"), 2.0000000000000004) && close(sc("exp(1)"), exp(1.0)) && close(sc("log(10)"), log(10.0)) && close(sc("sin(1)"), sin(1.0)), "sqrt, exp, log, sin");
    check(sc("abs(-3.5)") == 3.5 && sc("floor(2.7)") == 2 && sc("ceil(2.2)") == 3 && sc("round(2.5)") == 3 && sc("int(-2.7)") == -2, "abs, floor, ceil, round, int");
    check(sc("3 > 2") == 1 && sc("3 < 2") == 0 && sc("2 >= 2") == 1 && sc("2 != 2") == 0 && sc("2 == 2") == 1, "comparisons");
    check(sc("(1 > 0) && (2 > 1)") == 1 && sc("(1 > 0) && (2 < 1)") == 0 && sc("(0 > 1) || (2 > 1)") == 1 && sc("!(1 > 0)") == 0, "logical operators");
    check(sc("(3 > 2) ? 10 : 20") == 10 && sc("(3 < 2) ? 10 : 20") == 20, "the ternary operator");
    check(close(sc("$pi"), PI) && isNaN(sc("NA + 1")), "$pi and NA propagation");
    check(close(sc("gammafun(5)"), 24.0, 1e-12) && close(sc("cnorm(0)"), 0.5) && close(sc("dnorm(0)"), 1.0 / sqrt(2 * PI)) && close(sc("qnorm(0.975)"), 1.959963984540054, 1e-9), "gammafun, cnorm, dnorm, qnorm");
    check(sc("atan2(1, 1)") > 0.785398 && sc("atan2(1, 1)") < 0.785399, "atan2");

    /* ---- series ---- */
    {
        double[] lag(double[T] a, int k) { auto r = new double[T]; foreach (t; 0 .. T) r[t] = (t - k >= 0 && t - k < T) ? a[t - k] : double.nan; return r; }
        double[] sum2 = new double[T], sum3 = new double[T];

        foreach (t; 0 .. T) { sum2[t] = yv[t] + xv[t]; sum3[t] = yv[t] * 2 - xv[t] / 2; }
        check(sameSeries(se("y + x"), sum2) && sameSeries(se("y * 2 - x / 2"), sum3), "series arithmetic");
        check(sameSeries(se("y(-1)"), lag(yv, 1)) && sameSeries(se("y(-2)"), lag(yv, 2)), "lags y(-1), y(-2)");
        check(sameSeries(se("y(+1)"), lag(yv, -1)), "a lead y(+1)");
        auto dy = new double[T];

        foreach (t; 0 .. T) dy[t] = (t >= 1) ? yv[t] - yv[t - 1] : double.nan;
        check(sameSeries(se("diff(y)"), dy), "diff");
        auto ly = new double[T];

        foreach (t; 0 .. T) ly[t] = (t >= 1) ? log(yv[t]) - log(yv[t - 1]) : double.nan;
        check(sameSeries(se("ldiff(y)"), ly), "ldiff");
        auto cu = new double[T];
        double acc = 0;

        foreach (t; 0 .. T) { acc += yv[t]; cu[t] = acc; }
        check(sameSeries(se("cum(y)"), cu), "cum");
        auto srt = yv.dup;

        srt.length = T;
        foreach (i; 0 .. T) foreach (j; i + 1 .. T) if (srt[j] < srt[i]) { auto tmp = srt[i]; srt[i] = srt[j]; srt[j] = tmp; }
        check(sameSeries(se("sort(y)"), srt), "sort");
        auto cnd = new double[T];

        foreach (t; 0 .. T) cnd[t] = (xv[t] > 2) ? yv[t] : -1.0;
        check(sameSeries(se("(x > 2) ? y : -1"), cnd), "a conditional series");
        auto tm = new double[T];

        foreach (t; 0 .. T) tm[t] = t + 1;
        check(sameSeries(se("sqrt(y) ^ 2"), yv.dup), "sqrt then square of a series");
        double sy = 0, sy2 = 0, mn = yv[0], mx = yv[0];

        foreach (t; 0 .. T) { sy += yv[t]; sy2 += yv[t] * yv[t]; if (yv[t] < mn) mn = yv[t]; if (yv[t] > mx) mx = yv[t]; }
        const meany = sy / T, sdy = sqrt((sy2 - T * meany * meany) / (T - 1));

        check(close(sc("sum(y)"), sy) && close(sc("mean(y)"), meany) && close(sc("sd(y)"), sdy) && close(sc("min(y)"), mn) && close(sc("max(y)"), mx) && sc("sum(x > 2)") == 6, "sum, mean, sd, min, max of a series; counting a condition");
        double cxy = 0, mx_ = 0;

        foreach (t; 0 .. T) mx_ += xv[t] / T;
        foreach (t; 0 .. T) cxy += (yv[t] - meany) * (xv[t] - mx_);
        double sx2 = 0;

        foreach (t; 0 .. T) sx2 += (xv[t] - mx_) ^^ 2;
        const sdx = sqrt(sx2 / (T - 1));

        check(close(sc("corr(y, x)"), cxy / ((T - 1) * sdy * sdx), 1e-10), "corr");
    }

    /* ---- matrices ---- */
    {
        check(sameMatrix(ma("{1, 2; 3, 4}"), [[1, 2], [3, 4]]), "a matrix literal");
        check(sameMatrix(ma("{1, 2; 3, 4} * {5, 6; 7, 8}"), [[19, 22], [43, 50]]), "matrix product");
        check(sameMatrix(ma("{1, 2; 3, 4}'"), [[1, 3], [2, 4]]), "transpose");
        check(sameMatrix(ma("{1, 2; 3, 4} .* {5, 6; 7, 8}"), [[5, 12], [21, 32]]) && sameMatrix(ma("{1, 2; 3, 4} + 10"), [[11, 12], [13, 14]]), "elementwise product and scalar addition");
            check(sameMatrix(ma("inv({2, 0; 0, 4})"), [[0.5, 0], [0, 0.25]]) && close(sc("det({1, 2; 3, 4})"), -2) && close(sc("tr({1, 2; 3, 4})"), 5), "inverse, determinant, trace");
        check(sameMatrix(ma("{1, 2; 3, 4} ~ {5; 6}"), [[1, 2, 5], [3, 4, 6]]) && sameMatrix(ma("{1, 2} | {3, 4}"), [[1, 2], [3, 4]]), "horizontal (~) and vertical (|) concatenation");
            check(sc("rows({1, 2; 3, 4; 5, 6})") == 3 && sc("cols({1, 2; 3, 4; 5, 6})") == 2, "rows and cols");
        check(sameMatrix(ma("sumc({1, 2; 3, 4})"), [[4, 6]]) && sameMatrix(ma("meanc({1, 2; 3, 4})"), [[2, 3]]) && sameMatrix(ma("sumr({1, 2; 3, 4})"), [[3], [7]]), "sumc, meanc, sumr");
        check(sc("{1, 2; 3, 4}[2,1]") == 3 && sameMatrix(ma("{1, 2, 3; 4, 5, 6}[2,]"), [[4, 5, 6]]) && sameMatrix(ma("{1, 2, 3; 4, 5, 6}[,2:3]"), [[2, 3], [5, 6]]), "element, row and column-range slices");
        check(sameMatrix(ma("I(2)"), [[1, 0], [0, 1]]) && sameMatrix(ma("ones(1,3)"), [[1, 1, 1]]) && sameMatrix(ma("zeros(2,1)"), [[0], [0]]) && sameMatrix(ma("seq(1,4)"), [[1, 2, 3, 4]]), "I, ones, zeros, seq");
        check(sc("sum(({1, 2; 3, 4} * {1; 1}) .= {3; 7})") == 2, "elementwise comparison of matrices");
        // a linear system solved by the \ operator: x = A^-1 b
        auto sol = ma("{4, 1; 1, 3} \\ {1; 2}");

        check(sol.length == 2 && close(sol[0][0], 1.0 / 11) && close(sol[1][0], 7.0 / 11), "matrix left division solves a linear system");
        auto ev = ma("eigensym({2, 1; 1, 2})'");

        check(ev.length == 1 && ev[0].length == 2 && close(ev[0][0], 1) && close(ev[0][1], 3), "eigenvalues of a symmetric matrix");
        check(close(sc("rcond({1, 0; 0, 1})"), 1.0), "rcond of the identity");
        auto cm = ma("cholesky({4, 2; 2, 3})");

        check(cm.length == 2 && close(cm[0][0], 2) && close(cm[1][0], 1) && close(cm[1][1], sqrt(2.0)) && cm[0][1] == 0, "Cholesky factor");
    }

    /* ---- strings ---- */
    {
        check(run("str1 = \"abc\" ~ \"def\"", GretlType.string_) == 0 && strcmp(get_string_by_name("str1"), "abcdef") == 0, "string concatenation with ~");
        check(sc("strlen(\"hello\")") == 5, "strlen");
        check(run("str2 = strsub(\"hello world\", \"o\", \"0\")", GretlType.string_) == 0 && strcmp(get_string_by_name("str2"), "hell0 w0rld") == 0, "strsub");
        check(run("str3 = sprintf(\"%d-%s-%.2f\", 7, \"a\", 3.14159)", GretlType.string_) == 0 && strcmp(get_string_by_name("str3"), "7-a-3.14") == 0, "sprintf");
        check(run("str4 = toupper(\"abc\")", GretlType.string_) == 0 && strcmp(get_string_by_name("str4"), "ABC") == 0, "toupper");
    }

    /* ---- errors ---- */
    check(run("e1 = 2 +", GretlType.double_) != 0, "an incomplete expression is an error");
    check(run("e2 = (2 + 3", GretlType.double_) != 0, "an unmatched parenthesis is an error");
    check(run("e3 = nosuchvar + 1", GretlType.double_) != 0, "an undefined name is an error");
    check(run("e4 = nosuchfunc(2)", GretlType.double_) != 0, "an unknown function is an error");
    check(run("e5 = {1, 2} * {3, 4}", GretlType.matrix) != 0, "non-conformable matrices are an error");
    check(run("e6 = y + ", GretlType.series) != 0, "an incomplete series expression is an error");

    printf("\n%d tests, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
