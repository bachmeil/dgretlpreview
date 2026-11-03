# correlogram

## Signature

```d
Correlogram correlogram(X)(X x, int maxLag);
Correlogram correlogram(Dataset d, string name, int maxLag);

struct Correlogram {
    GretlMatrix acf;
    GretlMatrix pacf;
    GretlMatrix ljungBox;
}
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `maxLag`: the number of lags, from 1 to `n - 1`.

## Return value

A `Correlogram` struct with three `maxLag` x 1 matrices, element `k - 1` being for lag `k`:

- `acf`: the sample autocorrelation function.
- `pacf`: the sample partial autocorrelation function.
- `ljungBox`: the Ljung-Box statistic `Q(k)` for each `k` from 1 to `maxLag`.

For a white-noise series the autocorrelations are within about `2 / sqrt(n)` of zero. An AR(p) process has a partial autocorrelation function that is near zero after lag `p`. In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

## Example

```d
import std.stdio;
import std.math: sqrt;
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
    setSeed(14);
    auto x = arSeries(300, 0.7);
    auto c = correlogram(x, 6);
    foreach (k; 0 .. 6) {
        writefln("lag %d: acf %7.4f  pacf %7.4f  Q %8.3f", k + 1, c.acf[k, 0], c.pacf[k, 0], c.ljungBox[k, 0]);
    }
    writefln("two standard errors: %.3f", 2.0 / sqrt(300.0));
}
```

Output:

```
lag 1: acf  0.6911  pacf  0.6911  Q  144.732
lag 2: acf  0.4344  pacf -0.0829  Q  202.092
lag 3: acf  0.2643  pacf -0.0068  Q  223.394
lag 4: acf  0.1679  pacf  0.0170  Q  232.026
lag 5: acf  0.0905  pacf -0.0351  Q  234.538
lag 6: acf  0.0224  pacf -0.0405  Q  234.694
two standard errors: 0.115
```
