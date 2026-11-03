# qgamma

## Signature

```d
double qgamma(double p, double shape, double scale = 1, bool lower = true);
double[] qgamma(const(double)[] p, double shape, double scale = 1, bool lower = true);
GretlMatrix qgamma(GretlMatrix p, double shape, double scale = 1, bool lower = true);      // also Submatrix
double[] qgamma(Dataset ds, string series, double shape, double scale = 1, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qgamma(ds, "series", ...)`. It should lie in [0, 1].
- `shape`: the shape parameter. It must be positive and finite. It has no default.
- `scale`: the scale parameter (not the rate). It must be positive and finite. The mean of the distribution is `shape * scale`. The default is `1`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qgamma(p, scale: values)` is allowed.

## Return value

The quantile function of the gamma distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pgamma`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qgamma(0.95, 3, 2);
    writefln("%.6f", x);
    writefln("%.6f", pgamma(x, 3, 2));                 // P(X <= x) at that point
    writefln("%.6f", qgamma(0.05, 3, 2, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qgamma(ps, 3, 2));         // an array
    writefln("%(%.6f %)", qgamma(ps, 3, [1.0, 2, 3]));  // the parameter varies by element
}
```

Output:

```
12.591587
0.950000
12.591587
2.204131 6.210757 10.644641
1.102065 6.210757 15.966961
```
