# qged

## Signature

```d
double qged(double p, double shape, bool lower = true);
double[] qged(const(double)[] p, double shape, bool lower = true);
GretlMatrix qged(GretlMatrix p, double shape, bool lower = true);      // also Submatrix
double[] qged(Dataset ds, string series, double shape, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qged(ds, "series", ...)`. It should lie in [0, 1].
- `shape`: the shape parameter, in the library's parameterization. It must be positive and finite. A shape of 2 is the standard normal distribution; smaller values give heavier tails. It has no default.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The quantile function of the generalized error distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pged`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qged(0.95, 1.5);
    writefln("%.6f", x);
    writefln("%.6f", pged(x, 1.5));                 // P(X <= x) at that point
    writefln("%.6f", qged(0.05, 1.5, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qged(ps, 1.5));         // an array
    writefln("%(%.6f %)", qged(ps, [1.0, 1.5, 2.0]));  // the parameter varies by element
}
```

Output:

```
1.652739
0.950000
1.652739
-1.238020 0.216909 1.238020
-1.138044 0.216909 1.281552
```
