# hpCycle

## Signature

```d
GretlMatrix hpCycle(X)(X x, double lambda = 1600);
GretlMatrix hpCycle(Dataset d, string name, double lambda = 1600);
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `lambda`: the smoothing parameter, which must be positive: 1600 for quarterly data (the default), 100 for annual, 14400 for monthly. At least 4 observations are needed.

## Return value

An n x 1 `GretlMatrix` with the cyclical component of the Hodrick-Prescott filter, `x - hpTrend(x, lambda)`. In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

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
    auto x = arSeries(120, 0.95);
    auto cycle = hpCycle(x);

    double s = 0;
    foreach (t; 0 .. 120) {
        s += cycle[t, 0];
    }
    writefln("the cycle averages %.6f", s / 120);       // about 0
    writefln("standard deviation of x %.3f, of the cycle %.3f", stddev(x), stddev(Col(cycle, 0)));
}
```

Output:

```
the cycle averages 0.000000
standard deviation of x 3.088, of the cycle 1.224
```
