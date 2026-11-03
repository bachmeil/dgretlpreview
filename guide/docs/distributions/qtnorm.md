# qtnorm

## Signature

```d
double qtnorm(double p, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
double[] qtnorm(const(double)[] p, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
GretlMatrix qtnorm(GretlMatrix p, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);      // also Submatrix
double[] qtnorm(Dataset ds, string series, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qtnorm(ds, "series", ...)`. It should lie in [0, 1].
- `mean`: the mean of the normal distribution before truncation. It must be finite. The default is `0`.
- `sd`: the standard deviation of the normal distribution before truncation. It must be positive and finite. The default is `1`.
- `lo`: the lower bound of the truncation. The default means no lower bound. The default is `-double.infinity`.
- `hi`: the upper bound of the truncation. It must be greater than `lo`. The default means no upper bound. The default is `double.infinity`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qtnorm(p, hi: values)` is allowed.

## Return value

The quantile function of the truncated normal distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `ptnorm`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

`mean` and `sd` describe the normal distribution before truncation, not the mean and standard deviation of the truncated distribution. Parameters can be skipped and given by name, as in `ptnorm(x, hi: 2)` (a normal distribution truncated above at 2).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qtnorm(0.95, 1, 2, 0, 4);
    writefln("%.6f", x);
    writefln("%.6f", ptnorm(x, 1, 2, 0, 4));                 // P(X <= x) at that point
    writefln("%.6f", qtnorm(0.05, 1, 2, 0, 4, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qtnorm(ps, 1, 2, 0, 4));         // an array
    writefln("%(%.6f %)", qtnorm(ps, mean: [0.0, 1, 2]));  // the parameter varies by element
}
```

Output:

```
3.585602
0.950000
3.585602
0.341604 1.954066 3.259671
-1.281552 1.253347 3.281552
```
