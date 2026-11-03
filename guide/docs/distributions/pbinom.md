# pbinom

## Signature

```d
double pbinom(double q, double size, double prob, bool lower = true);
double[] pbinom(const(double)[] q, double size, double prob, bool lower = true);
GretlMatrix pbinom(GretlMatrix q, double size, double prob, bool lower = true);      // also Submatrix
double[] pbinom(Dataset ds, string series, double size, double prob, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `q`: the point (or points) at which to evaluate the distribution function. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `pbinom(ds, "series", ...)`.
- `size`: the number of trials, a whole number of at least 0. It has no default.
- `prob`: the probability of success on each trial, from 0 to 1. It has no default.
- `lower`: if `true` (the default) the function returns the lower tail, P(X <= q); if `false` it returns the upper tail, P(X > q), which keeps its accuracy far out in the tail where `1 - p` would round to 0.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The distribution function of the binomial distribution, P(X <= q) (or P(X > q) with `lower: false`); the support is the whole numbers 0 to `size`, so below it the lower tail is 0 (the upper tail 1) and above it the lower tail is 1 (the upper tail 0). For a discrete distribution a non-integer `q` is rounded down: P(X <= q) is P(X <= floor(q)).

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", pbinom(3.0, 10, 0.3));
    writefln("%.6f", pbinom(3.0, 10, 0.3, lower: false));    // the upper tail, P(X > q)
    auto qs = [2.0, 3.0, 4.0];
    writefln("%(%.6f %)", pbinom(qs, 10, 0.3));         // an array
    writefln("%(%.6f %)", pbinom(qs, 10, [0.2, 0.3, 0.4]));  // the parameter varies by element
    auto m = cbind([2.0, 3.0, 4.0], [2.0, 3.0, 4.0]);
    auto pm = pbinom(m, 10, 0.3);                          // a matrix gives a matrix
    writefln("%d x %d, element [1, 1] = %.6f", pm.nrow, pm.ncol, pm[1, 1]);
}
```

Output:

```
0.649611
0.350389
0.382783 0.649611 0.849732
0.677800 0.649611 0.633103
3 x 2, element [1, 1] = 0.649611
```
