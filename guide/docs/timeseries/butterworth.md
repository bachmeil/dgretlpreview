# butterworth

## Signature

```d
GretlMatrix butterworth(X)(X x, int order = 2, double cutoff = 67.5);
GretlMatrix butterworth(Dataset d, string name, int order = 2, double cutoff = 67.5);
```

## Parameters

- `x`: the series: an array of numbers, or a `GretlMatrix` or `Submatrix` with one column. Missing values (NaN) are not allowed and throw a `GretlException`.
- `order`: the order of the filter, at least 1. The default is 2.
- `cutoff`: the cutoff frequency in degrees, strictly between 0 and 180. The default is 67.5. A lower cutoff keeps only slower movements.

## Return value

An n x 1 `GretlMatrix` with the low-pass Butterworth filtered series: the slowly varying part of `x`. In the `Dataset` form, `name` is the name of a series of the `Dataset` instead of the data themselves (the series must have no missing values).

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
    setSeed(13);
    auto x = arSeries(100, 0.3);
    auto smooth = butterworth(x);
    auto smoother = butterworth(x, 2, 30);           // a lower cutoff: smoother still

    writefln("standard deviation of x %.3f, filtered %.3f, filtered with a lower cutoff %.3f", stddev(x), stddev(Col(smooth, 0)), stddev(Col(smoother, 0)));
}
```

Output:

```
standard deviation of x 0.985, filtered 0.670, filtered with a lower cutoff 0.437
```
