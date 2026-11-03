# kpss

## Signature

```d
UnitRootTest kpss(X)(X x, int lags, Trend trend = Trend.constant, bool difference = false);
UnitRootTest kpss(Dataset d, string name, int lags, Trend trend = Trend.constant, bool difference = false);

enum Trend { none, constant, trend, quadratic }
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`. It needs more than `lags + 2` observations.
- `lags`: the number of lags in the Bartlett window used for the long-run variance, from 0 to `n - 3`. A common choice is `4 * (n / 100)^0.25`.
- `trend`: `Trend.constant` (the default) tests stationarity around a constant, `Trend.trend` around a linear trend. No other value is meaningful, and others throw a `GretlException`.
- `difference`: if `true`, test the first difference of `x`. The default is `false`.

## Return value

A `UnitRootTest` struct (see `adf`) with `test` "KPSS", the `statistic`, the `pvalue` (NaN if the library gives none: the KPSS critical values are in the `report`) and `report`. The null hypothesis is the opposite of the Dickey-Fuller one: that the series is stationary. A large statistic (a small p-value) is evidence against stationarity.

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
    setSeed(18);
    auto walk = arSeries(200, 1.0);
    auto stationary = arSeries(200, 0.5);

    auto a = kpss(stationary, 4);
    writefln("stationary: statistic %.3f", a.statistic);
    auto b = kpss(walk, 4);
    writefln("random walk: statistic %.3f", b.statistic);          // much larger: stationarity is rejected
    auto c = kpss(walk, 4, Trend.trend);
    writefln("random walk, around a trend: statistic %.3f", c.statistic);
}
```

Output:

```
stationary: statistic 0.463
random walk: statistic 1.062
random walk, around a trend: statistic 0.386
```
