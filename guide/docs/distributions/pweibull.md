# pweibull

## Signature

```d
double pweibull(double q, double shape, double scale = 1, bool lower = true);
double[] pweibull(const(double)[] q, double shape, double scale = 1, bool lower = true);
GretlMatrix pweibull(GretlMatrix q, double shape, double scale = 1, bool lower = true);      // also Submatrix
double[] pweibull(Dataset ds, string series, double shape, double scale = 1, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pweibull(ds, "series", ...)`.
- `shape`: the shape parameter. It must be positive and finite. It has no default.
- `scale`: the scale parameter. It must be positive and finite. The default is `1`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `pweibull(q, scale: values)` is allowed.

## Return value

The distribution function of the Weibull distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is x >= 0, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pweibull(1.0, 1.5, 2));
    writefln("%.6f", pweibull(1.0, 1.5, 2, lower: false));    // the upper tail, P(X > q)
    auto qs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", pweibull(qs, 1.5, 2));         // an array
    writefln("%(%.6f %)", pweibull(qs, 1.5, [1.0, 2, 3]));  // the parameter varies by element
    auto m = cbind([0.5, 1.0, 3.0], [0.5, 1.0, 3.0]);
    auto pm = pweibull(m, 1.5, 2);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.297811
0.702189
0.117503 0.297811 0.840724
0.297811 0.297811 0.632121
3 x 2, element [1, 1] = 0.297811
```
