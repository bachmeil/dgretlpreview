# qchisq

## Signature

```d
double qchisq(double p, double df, double ncp = 0, bool lower = true);
double[] qchisq(const(double)[] p, double df, double ncp = 0, bool lower = true);
GretlMatrix qchisq(GretlMatrix p, double df, double ncp = 0, bool lower = true);      // also Submatrix
double[] qchisq(Dataset ds, string series, double df, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qchisq(ds, "series", ...)`. It should lie in [0, 1].
- `df`: the degrees of freedom. It must be positive and finite; unlike the F distribution it need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary chi-square distribution. The default is `0`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qchisq(p, ncp: values)` is allowed.

## Return value

The quantile function of the chi-square distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pchisq`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

`qchisq(0.95, 3)` is 7.815, the usual 5% critical value.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qchisq(0.95, 3);
    writefln("%.6f", x);
    writefln("%.6f", pchisq(x, 3));                 // P(X <= x) at that point
    writefln("%.6f", qchisq(0.05, 3, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qchisq(ps, 3));         // an array
    writefln("%(%.6f %)", qchisq(ps, [2.0, 3, 4]));  // the parameter varies by element
    writefln("%.6f", qchisq(0.95, 3, ncp: 2));         // non-central
}
```

Output:

```
7.814728
0.950000
7.814728
0.584374 2.946166 6.251389
0.210721 2.946166 7.779440
12.283754
```
