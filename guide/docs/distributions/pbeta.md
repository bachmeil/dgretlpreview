# pbeta

## Signature

```d
double pbeta(double q, double shape1, double shape2, bool lower = true);
double[] pbeta(const(double)[] q, double shape1, double shape2, bool lower = true);
GretlMatrix pbeta(GretlMatrix q, double shape1, double shape2, bool lower = true);      // also Submatrix
double[] pbeta(Dataset ds, string series, double shape1, double shape2, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pbeta(ds, "series", ...)`.
- `shape1`: the first shape parameter. It must be positive and finite. It has no default.
- `shape2`: the second shape parameter. It must be positive and finite. It has no default.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The distribution function of the beta distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is 0 <= x <= 1, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pbeta(0.25, 2, 6));
    writefln("%.6f", pbeta(0.25, 2, 6, lower: false));    // the upper tail, P(X > q)
    auto qs = [0.1, 0.25, 0.5];
    writefln("%(%.6f %)", pbeta(qs, 2, 6));         // an array
    writefln("%(%.6f %)", pbeta(qs, 2, [4.0, 6, 8]));  // the parameter varies by element
    auto m = cbind([0.1, 0.25, 0.5], [0.1, 0.25, 0.5]);
    auto pm = pbeta(m, 2, 6);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.555054
0.444946
0.149694 0.555054 0.937500
0.081460 0.555054 0.980469
3 x 2, element [1, 1] = 0.555054
```
