# engleGranger

## Signature

```d
UnitRootTest engleGranger(Y, X)(Y y, X x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off);
UnitRootTest engleGranger(Dataset d, string y, string[] x, int lags = 0, Trend trend = Trend.constant, TestDown testDown = TestDown.off);

enum Trend { none, constant, trend, quadratic }
enum TestDown { off, tstat, aic, bic }
```

## Parameters

- `y`: the dependent series (an array or a one-column matrix).
- `x`: the other series, as the columns of a matrix (or a single array). It must have as many rows as `y`. In the `Dataset` form, `y` and `x` are the names of the series.
- `lags`: the number of lagged differences in the ADF regression on the residuals. The default is 0. More than `lags + 10 + ncol(x)` observations are needed.
- `trend`: the deterministic terms of the cointegrating regression: `Trend.none`, `Trend.constant` (the default), `Trend.trend` or `Trend.quadratic`.
- `testDown`: choose the lag order from 0 to `lags` by `TestDown.tstat`, `aic` or `bic`; the default `off` uses `lags` as given.

## Return value

A `UnitRootTest` struct (see `adf`) for the Engle-Granger test of the null hypothesis that `y` and the columns of `x` are not cointegrated. The test regresses `y` on `x` and applies an augmented Dickey-Fuller test to the residuals. `statistic` is its t-ratio, `pvalue` the library's (MacKinnon) p-value, and `report` shows the cointegrating regression and the unit-root tests of the individual series. A small p-value is evidence of cointegration.

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
    setSeed(19);
    auto x = arSeries(200, 1.0);                     // a random walk
    auto e = arSeries(200, 0.3);                     // a stationary disturbance
    auto y = new double[200];
    foreach (t; 0 .. 200) {
        y[t] = 1.0 + 2.0 * x[t] + e[t];              // y and x share a stochastic trend
    }

    auto t = engleGranger(y, x);
    writefln("cointegrated: statistic %.3f, p = %.4f", t.statistic, t.pvalue);

    auto u = arSeries(200, 1.0);                     // an unrelated random walk
    auto v = engleGranger(u, x);
    writefln("unrelated:    statistic %.3f, p = %.4f", v.statistic, v.pvalue);
}
```

Output:

```
cointegrated: statistic -10.902, p = 0.0000
unrelated:    statistic -2.789, p = 0.1739
```
