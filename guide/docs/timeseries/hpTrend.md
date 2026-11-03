# hpTrend

## Signature

```d
GretlMatrix hpTrend(X)(X x, double lambda = 1600);
GretlMatrix hpTrend(Dataset d, string name, double lambda = 1600);
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `lambda`: the smoothing parameter, which must be positive. The usual choices are 1600 for quarterly data (the default), 100 for annual and 14400 for monthly data. A larger value gives a smoother trend. At least 4 observations are needed.

## Return value

An n x 1 `GretlMatrix` with the trend component of the Hodrick-Prescott filter, the smooth series that minimizes the sum of squared deviations of `x` from it plus `lambda` times the sum of squared second differences of the trend. The cycle is `x - trend` (`hpCycle` returns it directly). In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

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
    setSeed(11);
    auto x = arSeries(120, 0.95);                    // a persistent series
    auto trend = hpTrend(x);
    auto cycle = hpCycle(x);
    foreach (t; 0 .. 4) {
        writefln("%8.4f %8.4f %8.4f", x[t], trend[t, 0], cycle[t, 0]);   // x = trend + cycle
    }

    auto smoother = hpTrend(x, 14400);               // a larger lambda: a smoother trend
    writefln("%.4f %.4f", trend[60, 0], smoother[60, 0]);
}
```

Output:

```
  0.3788  -0.0989   0.4777
  0.3970   0.0310   0.3659
 -0.7334   0.1612  -0.8947
  0.6600   0.2923   0.3678
6.1045 4.8638
```
