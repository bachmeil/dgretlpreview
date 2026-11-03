# qt

## Signature

```d
double qt(double p, double df, double ncp = 0, bool lower = true);
double[] qt(const(double)[] p, double df, double ncp = 0, bool lower = true);
GretlMatrix qt(GretlMatrix p, double df, double ncp = 0, bool lower = true);      // also Submatrix
double[] qt(Dataset ds, string series, double df, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qt(ds, "series", ...)`. It should lie in [0, 1].
- `df`: the degrees of freedom. It must be positive and finite, and need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite. The default 0 is the ordinary (central) t distribution; any other value gives the non-central t distribution. The default is `0`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qt(p, ncp: values)` is allowed.

## Return value

The quantile function of the Student's t distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pt`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qt(0.975, 10);
    writefln("%.6f", x);
    writefln("%.6f", pt(x, 10));                 // P(X <= x) at that point
    writefln("%.6f", qt(0.05, 10, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qt(ps, 10));         // an array
    writefln("%(%.6f %)", qt(ps, [5.0, 10, 20]));  // the parameter varies by element
    writefln("%.6f", qt(0.975, 10, ncp: 1.5));       // non-central
}
```

Output:

```
2.228139
0.975000
1.812461
-1.372184 0.260185 1.372184
-1.475884 0.260185 1.325341
4.237300
```
