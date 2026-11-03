# pged

## Signature

```d
double pged(double q, double shape, bool lower = true);
double[] pged(const(double)[] q, double shape, bool lower = true);
GretlMatrix pged(GretlMatrix q, double shape, bool lower = true);      // also Submatrix
double[] pged(Dataset ds, string series, double shape, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pged(ds, "series", ...)`.
- `shape`: the shape parameter, in the library's parameterization. It must be positive and finite. A shape of 2 is the standard normal distribution; smaller values give heavier tails. It has no default.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The distribution function of the generalized error distribution, P(X <= q) (or P(X > q) with `lower: false`); the tails reach 0 and 1 only at minus and plus infinity.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pged(0.5, 1.5));
    writefln("%.6f", pged(0.5, 1.5, lower: false));    // the upper tail, P(X > q)
    auto qs = [-1.0, 0.0, 1.0];
    writefln("%(%.6f %)", pged(qs, 1.5));         // an array
    writefln("%(%.6f %)", pged(qs, [1.0, 1.5, 2.0]));  // the parameter varies by element
    auto m = cbind([-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]);
    auto pm = pged(m, 1.5);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.713379
0.286621
0.144229 0.500000 0.855771
0.121558 0.500000 0.841345
3 x 2, element [1, 1] = 0.500000
```
