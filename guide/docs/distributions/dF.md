# dF

## Signature

```d
double dF(double x, double df1, double df2, double ncp = 0);
double[] dF(const(double)[] x, double df1, double df2, double ncp = 0);
GretlMatrix dF(GretlMatrix x, double df1, double df2, double ncp = 0);      // also Submatrix
double[] dF(Dataset ds, string series, double df1, double df2, double ncp = 0);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dF(ds, "series", ...)`.
- `df1`: the numerator degrees of freedom, a whole number of at least 1. It has no default.
- `df2`: the denominator degrees of freedom, a whole number of at least 1. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary F distribution. The default is `0`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dF(x, ncp: values)` is allowed.

## Return value

The density of the F distribution at each point: the F density with `df1` and `df2` degrees of freedom (the non-central density when `ncp` is not 0). The support is x >= 0; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

The library accepts only whole-number degrees of freedom for the F distribution, so a non-integer `df1` or `df2` throws.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dF(2.0, 4, 20));
    auto xs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", dF(xs, 4, 20));         // an array: one value per element
    writefln("%(%.6f %)", dF(xs, 4, [10, 20, 30]));  // the parameter varies by element
    writefln("%.6f", dF(2.0, 4, 20, ncp: 1));      // non-central
}
```

Output:

```
0.155219
0.700988 0.493489 0.046896
0.669796 0.493489 0.041979
0.204464
```
