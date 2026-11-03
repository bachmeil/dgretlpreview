# qbinom

## Signature

```d
double qbinom(double p, double size, double prob, bool lower = true);
double[] qbinom(const(double)[] p, double size, double prob, bool lower = true);
GretlMatrix qbinom(GretlMatrix p, double size, double prob, bool lower = true);      // also Submatrix
double[] qbinom(Dataset ds, string series, double size, double prob, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qbinom(ds, "series", ...)`. It should lie in [0, 1].
- `size`: the number of trials, a whole number of at least 0. It has no default.
- `prob`: the probability of success on each trial, from 0 to 1. It has no default.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The quantile function of the binomial distribution: the smallest whole number `x` with P(X <= x) >= p (with `lower: false`, the smallest `x` with P(X > x) <= p). It is the inverse of `pbinom` in that sense, and a value of `p` strictly between two cumulative probabilities gives the whole number above. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qbinom(0.95, 10, 0.3);
    writefln("%.6f", x);
    writefln("%.6f", pbinom(x, 10, 0.3));                 // P(X <= x) at that point
    writefln("%.6f", qbinom(0.05, 10, 0.3, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qbinom(ps, 10, 0.3));         // an array
    writefln("%(%.6f %)", qbinom(ps, 10, [0.2, 0.3, 0.4]));  // the parameter varies by element
}
```

Output:

```
5.000000
0.952651
5.000000
1.000000 3.000000 5.000000
0.000000 3.000000 6.000000
```
