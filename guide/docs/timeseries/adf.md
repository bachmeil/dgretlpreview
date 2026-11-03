# adf

## Signature

```d
UnitRootTest adf(X)(X x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off, bool gls = false,
                    bool seasonals = false, bool difference = false, int frequency = 1, string start = "1");
UnitRootTest adf(Dataset d, string name, int lags = 0, ...);   // the same options; frequency and start come from the Dataset

enum Trend { none, constant, trend, quadratic }
enum TestDown { off, tstat, aic, bic }

struct UnitRootTest {
    string test;
    double statistic;
    double pvalue;
    string report;
}
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`. More than `lags + 10` observations are needed.
- `lags`: the number of lagged differences in the augmented regression. With `testDown` it is the largest lag order tried. The default is 0 (the plain Dickey-Fuller test).
- `trend`: the deterministic terms of the test regression: `Trend.none` (no constant), `Trend.constant` (the default), `Trend.trend` (constant and linear trend) or `Trend.quadratic` (constant, linear and quadratic trend).
- `testDown`: choose the lag order from 0 to `lags` by a criterion: `TestDown.tstat` (start from `lags` and drop the last lag while its t-ratio is not significant at the 10% level), `TestDown.aic` or `TestDown.bic` (minimize the criterion). The default `TestDown.off` uses `lags` as given.
- `gls`: if `true`, use Elliott-Rothenberg-Stock GLS detrending (the DF-GLS test); `trend` must be `constant` or `trend`. The default is `false`.
- `seasonals`: if `true`, add seasonal dummies (the data need a frequency above 1). The default is `false`.
- `difference`: if `true`, test the first difference of `x`. The default is `false`.
- `frequency`, `start`: the calendar of the series, needed for the seasonal dummies: the observations per year (4, 12, ...) and the first observation ("1990:1" or "1990Q1"). The default is annual data from 1.

## Return value

A `UnitRootTest` struct with `test` ("ADF", or "DF-GLS" when `gls` is set), the `statistic` (the t-ratio of the lagged level), the `pvalue` and `report`, the text the library prints (it shows the regression and, with `testDown`, the lag order chosen). The null hypothesis is that the series has a unit root; a small p-value is evidence against it. The p-value is NaN for a case in which the library has no p-value.

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
    setSeed(17);
    auto walk = arSeries(200, 1.0);                  // a random walk: it has a unit root
    auto stationary = arSeries(200, 0.5);

    auto a = adf(walk);
    writefln("random walk: statistic %.3f, p = %.3f", a.statistic, a.pvalue);

    auto b = adf(stationary);
    writefln("AR(0.5):     statistic %.3f, p = %.3f", b.statistic, b.pvalue);

    auto c = adf(walk, 4, Trend.trend, TestDown.aic);    // with a trend, choosing the lags (up to 4) by AIC
    writefln("with trend and AIC lag selection: statistic %.3f, p = %.3f", c.statistic, c.pvalue);

    auto g = adf(walk, 2, Trend.constant, TestDown.off, true);   // DF-GLS
    writefln("%s: statistic %.3f", g.test, g.statistic);
}
```

Output:

```
random walk: statistic -1.847, p = 0.357
AR(0.5):     statistic -7.175, p = 0.000
with trend and AIC lag selection: statistic -0.673, p = 0.974
DF-GLS: statistic 1.149
```
