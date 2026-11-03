# pchisq

## Signature

```d
double pchisq(double q, double df, double ncp = 0, bool lower = true);
double[] pchisq(const(double)[] q, double df, double ncp = 0, bool lower = true);
GretlMatrix pchisq(GretlMatrix q, double df, double ncp = 0, bool lower = true);      // also Submatrix
double[] pchisq(Dataset ds, string series, double df, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pchisq(ds, "series", ...)`.
- `df`: the degrees of freedom. It must be positive and finite; unlike the F distribution it need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary chi-square distribution. The default is `0`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `pchisq(q, ncp: values)` is allowed.

## Return value

The distribution function of the chi-square distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is x >= 0, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pchisq(3.0, 3));
    writefln("%.6f", pchisq(3.0, 3, lower: false));    // the upper tail, P(X > q)
    auto qs = [1.0, 3.0, 7.815];
    writefln("%(%.6f %)", pchisq(qs, 3));         // an array
    writefln("%(%.6f %)", pchisq(qs, [2.0, 3, 4]));  // the parameter varies by element
    auto m = cbind([1.0, 3.0, 7.815], [1.0, 3.0, 7.815]);
    auto pm = pchisq(m, 3);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
    writefln("%.6f", pchisq(3.0, 3, ncp: 2));         // non-central
}
```

Output:

```
0.608375
0.391625
0.198748 0.608375 0.950006
0.393469 0.608375 0.901405
3 x 2, element [1, 1] = 0.608375
0.357668
```
