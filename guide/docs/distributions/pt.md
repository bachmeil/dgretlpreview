# pt

## Signature

```d
double pt(double q, double df, double ncp = 0, bool lower = true);
double[] pt(const(double)[] q, double df, double ncp = 0, bool lower = true);
GretlMatrix pt(GretlMatrix q, double df, double ncp = 0, bool lower = true);      // also Submatrix
double[] pt(Dataset ds, string series, double df, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pt(ds, "series", ...)`.
- `df`: the degrees of freedom. It must be positive and finite, and need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite. The default 0 is the ordinary (central) t distribution; any other value gives the non-central t distribution. The default is `0`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `pt(q, ncp: values)` is allowed.

## Return value

The distribution function of the Student's t distribution, P(X <= q) (or P(X > q) with `lower: false`); the tails reach 0 and 1 only at minus and plus infinity.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pt(2.0, 10));
    writefln("%.6f", pt(2.0, 10, lower: false));    // the upper tail, P(X > q)
    auto qs = [-2.0, 0.0, 2.0];
    writefln("%(%.6f %)", pt(qs, 10));         // an array
    writefln("%(%.6f %)", pt(qs, [5.0, 10, 20]));  // the parameter varies by element
    auto m = cbind([-2.0, 0.0, 2.0], [-2.0, 0.0, 2.0]);
    auto pm = pt(m, 10);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
    writefln("%.6f", pt(2.0, 10, ncp: 1.5));       // non-central
}
```

Output:

```
0.963306
0.036694
0.036694 0.500000 0.963306
0.050970 0.500000 0.970367
3 x 2, element [1, 1] = 0.500000
0.659154
```
