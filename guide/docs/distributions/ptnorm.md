# ptnorm

## Signature

```d
double ptnorm(double q, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
double[] ptnorm(const(double)[] q, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
GretlMatrix ptnorm(GretlMatrix q, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);      // also Submatrix
double[] ptnorm(Dataset ds, string series, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `ptnorm(ds, "series", ...)`.
- `mean`: the mean of the normal distribution before truncation. It must be finite. The default is `0`.
- `sd`: the standard deviation of the normal distribution before truncation. It must be positive and finite. The default is `1`.
- `lo`: the lower bound of the truncation. The default means no lower bound. The default is `-double.infinity`.
- `hi`: the upper bound of the truncation. It must be greater than `lo`. The default means no upper bound. The default is `double.infinity`.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `ptnorm(q, hi: values)` is allowed.

## Return value

The distribution function of the truncated normal distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is lo <= x <= hi, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

`mean` and `sd` describe the normal distribution before truncation, not the mean and standard deviation of the truncated distribution. Parameters can be skipped and given by name, as in `ptnorm(x, hi: 2)` (a normal distribution truncated above at 2).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", ptnorm(1.0, 1, 2, 0, 4));
    writefln("%.6f", ptnorm(1.0, 1, 2, 0, 4, lower: false));    // the upper tail, P(X > q)
    auto qs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", ptnorm(qs, 1, 2, 0, 4));         // an array
    writefln("%(%.6f %)", ptnorm(qs, mean: [0.0, 1, 2]));  // the parameter varies by element
    auto m = cbind([0.5, 1.0, 3.0], [0.5, 1.0, 3.0]);
    auto pm = ptnorm(m, 1, 2, 0, 4);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.306509
0.693491
0.148492 0.306509 0.852962
0.691462 0.500000 0.841345
3 x 2, element [1, 1] = 0.306509
```
