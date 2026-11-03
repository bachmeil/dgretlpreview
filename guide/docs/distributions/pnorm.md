# pnorm

## Signature

```d
double pnorm(double q, double mean = 0, double sd = 1, bool lower = true);
double[] pnorm(const(double)[] q, double mean = 0, double sd = 1, bool lower = true);
GretlMatrix pnorm(GretlMatrix q, double mean = 0, double sd = 1, bool lower = true);      // also Submatrix
double[] pnorm(Dataset ds, string series, double mean = 0, double sd = 1, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pnorm(ds, "series", ...)`.
- `mean`: the mean. It must be finite. The default is `0`.
- `sd`: the standard deviation. It must be positive and finite. The default is `1`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `pnorm(q, sd: values)` is allowed.

## Return value

The distribution function of the normal distribution, P(X <= q) (or P(X > q) with `lower: false`); the tails reach 0 and 1 only at minus and plus infinity.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

`pnorm(1.96)` is 0.975 to the accuracy shown; `qnorm(0.975)` returns about 1.96.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pnorm(0.5, 1, 2));
    writefln("%.6f", pnorm(0.5, 1, 2, lower: false));    // the upper tail, P(X > q)
    auto qs = [-1.0, 0.0, 1.0];
    writefln("%(%.6f %)", pnorm(qs, 1, 2));         // an array
    writefln("%(%.6f %)", pnorm(qs, sd: [1.0, 2, 3]));  // the parameter varies by element
    auto m = cbind([-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]);
    auto pm = pnorm(m, 1, 2);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.401294
0.598706
0.158655 0.308538 0.500000
0.158655 0.500000 0.630559
3 x 2, element [1, 1] = 0.308538
```
