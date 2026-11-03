# ljungBox

## Signature

```d
PortmanteauTest ljungBox(X)(X x, int lags, int fitted = 0);
TestResult ljungBox(Dataset d, string name, int lags, int fitted = 0);

struct TestResult {
    double statistic;
    double pvalue;
    int df;
}
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `lags`: the number of autocorrelations used, from 1 to `n - 1`.
- `fitted`: the number of ARMA parameters that were estimated, when `x` holds the residuals of a fitted model. It is subtracted from the degrees of freedom and must be less than `lags`. The default is 0.

## Return value

A `TestResult` (also called `PortmanteauTest`) for the Ljung-Box test of the null hypothesis of no autocorrelation up to lag `lags`: `statistic` is `Q`, `df` is `lags - fitted` and `pvalue` is the upper tail of the chi-square distribution with `df` degrees of freedom. A small p-value is evidence of autocorrelation. In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

## Example

```d
import std.stdio;
import dgretl;

double[] arSeries(size_t n, double phi) {     // an AR(1) series from the library's random numbers
    auto e = rnorm(n);
    auto y = new double[n];

    y[0] = e[0];
    foreach (t; 1 .. n) {
        y[t] = phi * y[t - 1] + e[t];
    }
    return y;
}

void main() {
    setSeed(15);
    auto noise = rnorm(200);
    auto t = ljungBox(noise, 10);
    writefln("white noise: Q = %.3f, df = %d, p = %.4f", t.statistic, t.df, t.pvalue);

    auto ar = arSeries(200, 0.5);
    auto u = ljungBox(ar, 10);
    writefln("AR(1): Q = %.3f, p = %.6f", u.statistic, u.pvalue);

    auto r = ljungBox(noise, 10, 2);                 // residuals of a model with 2 estimated parameters
    writeln(r.df);
}
```

Output:

```
white noise: Q = 4.628, df = 10, p = 0.9146
AR(1): Q = 41.829, p = 0.000008
8
```
