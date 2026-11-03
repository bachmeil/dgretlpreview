# fractionalIntegration

## Signature

```d
FractionalIntegration fractionalIntegration(X)(X x, int bandwidth = 0);

struct FractionalIntegration {
    int bandwidth;
    double lweD, lweSe, lweZ, lwePvalue;
    double gphD, gphSe, gphT, gphPvalue;
}
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`. At least 12 observations are needed.
- `bandwidth`: the number of Fourier frequencies used, at least 2 and at most `n / 2`. The default 0 takes the library's choice, about `n^0.6`.

## Return value

A `FractionalIntegration` struct with estimates of the degree of fractional integration `d` of the series (0 for a stationary series, 1 for a unit root) by two methods, using the frequencies up to `bandwidth` (the number used is in `bandwidth`):

- the local Whittle estimator of Robinson: `lweD`, its standard error `lweSe` (`1 / (2 sqrt(m))`), `lweZ = lweD / lweSe` and the two-sided normal p-value `lwePvalue` for the null `d = 0`.
- the Geweke and Porter-Hudak log-periodogram regression estimator: `gphD`, `gphSe`, `gphT = gphD / gphSe` and the two-sided p-value `gphPvalue` from the t distribution with `bandwidth - 2` degrees of freedom.

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
    setSeed(21);
    auto noise = rnorm(500);
    auto f = fractionalIntegration(noise);
    writefln("white noise: local Whittle d = %.3f (se %.3f), GPH d = %.3f", f.lweD, f.lweSe, f.gphD);

    auto walk = arSeries(500, 1.0);
    auto g = fractionalIntegration(walk);
    writefln("random walk: local Whittle d = %.3f, GPH d = %.3f", g.lweD, g.gphD);
}
```

Output:

```
white noise: local Whittle d = -0.049 (se 0.079), GPH d = 0.008
random walk: local Whittle d = 0.936, GPH d = 0.937
```
