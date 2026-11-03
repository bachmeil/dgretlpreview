# pF

## Signature

```d
double pF(double q, double df1, double df2, double ncp = 0, bool lower = true);
double[] pF(const(double)[] q, double df1, double df2, double ncp = 0, bool lower = true);
GretlMatrix pF(GretlMatrix q, double df1, double df2, double ncp = 0, bool lower = true);      // also Submatrix
double[] pF(Dataset ds, string series, double df1, double df2, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pF(ds, "series", ...)`.
- `df1`: the numerator degrees of freedom, a whole number of at least 1. It has no default.
- `df2`: the denominator degrees of freedom, a whole number of at least 1. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary F distribution. The default is `0`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `pF(q, ncp: values)` is allowed.

## Return value

The distribution function of the F distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is x >= 0, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

The library accepts only whole-number degrees of freedom for the F distribution, so a non-integer `df1` or `df2` throws.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pF(2.0, 4, 20));
    writefln("%.6f", pF(2.0, 4, 20, lower: false));    // the upper tail, P(X > q)
    auto qs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", pF(qs, 4, 20));         // an array
    writefln("%(%.6f %)", pF(qs, 4, [10, 20, 30]));  // the parameter varies by element
    auto m = cbind([0.5, 1.0, 3.0], [0.5, 1.0, 3.0]);
    auto pm = pF(m, 4, 20);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
    writefln("%.6f", pF(2.0, 4, 20, ncp: 1));      // non-central
}
```

Output:

```
0.866652
0.133348
0.263963 0.569318 0.956799
0.263224 0.569318 0.966023
3 x 2, element [1, 1] = 0.569318
0.786959
```
