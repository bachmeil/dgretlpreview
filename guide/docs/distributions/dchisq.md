# dchisq

## Signature

```d
double dchisq(double x, double df, double ncp = 0);
double[] dchisq(const(double)[] x, double df, double ncp = 0);
GretlMatrix dchisq(GretlMatrix x, double df, double ncp = 0);      // also Submatrix
double[] dchisq(Dataset ds, string series, double df, double ncp = 0);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dchisq(ds, "series", ...)`.
- `df`: the degrees of freedom. It must be positive and finite; unlike the F distribution it need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary chi-square distribution. The default is `0`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dchisq(x, ncp: values)` is allowed.

## Return value

The density of the chi-square distribution at each point: `x^(df/2 - 1) exp(-x/2) / (2^(df/2) gamma(df/2))` when `ncp = 0` (the density of the non-central chi-square otherwise). The support is x >= 0; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dchisq(3.0, 3));
    auto xs = [1.0, 3.0, 7.815];
    writefln("%(%.6f %)", dchisq(xs, 3));         // an array: one value per element
    writefln("%(%.6f %)", dchisq(xs, [2.0, 3, 4]));  // the parameter varies by element
    writefln("%.6f", dchisq(3.0, 3, ncp: 2));         // non-central
}
```

Output:

```
0.154180
0.241971 0.154180 0.022406
0.303265 0.154180 0.039252
0.134100
```
