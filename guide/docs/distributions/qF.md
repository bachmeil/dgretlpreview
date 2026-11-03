# qF

## Signature

```d
double qF(double p, double df1, double df2, double ncp = 0, bool lower = true);
double[] qF(const(double)[] p, double df1, double df2, double ncp = 0, bool lower = true);
GretlMatrix qF(GretlMatrix p, double df1, double df2, double ncp = 0, bool lower = true);      // also Submatrix
double[] qF(Dataset ds, string series, double df1, double df2, double ncp = 0, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qF(ds, "series", ...)`. It should lie in [0, 1].
- `df1`: the numerator degrees of freedom, a whole number of at least 1. It has no default.
- `df2`: the denominator degrees of freedom, a whole number of at least 1. It has no default.
- `ncp`: the non-centrality parameter. It must be finite and not negative. The default 0 is the ordinary F distribution. The default is `0`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qF(p, ncp: values)` is allowed.

## Return value

The quantile function of the F distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pF`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

The library accepts only whole-number degrees of freedom for the F distribution, so a non-integer `df1` or `df2` throws.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qF(0.95, 4, 20);
    writefln("%.6f", x);
    writefln("%.6f", pF(x, 4, 20));                 // P(X <= x) at that point
    writefln("%.6f", qF(0.05, 4, 20, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qF(ps, 4, 20));         // an array
    writefln("%(%.6f %)", qF(ps, 4, [10, 20, 30]));  // the parameter varies by element
    writefln("%.6f", qF(0.95, 4, 20, ncp: 1));      // non-central
}
```

Output:

```
2.866081
0.950000
2.866081
0.260123 1.064218 2.248934
0.255110 1.064218 2.142235
3.546698
```
