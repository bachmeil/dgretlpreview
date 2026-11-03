# qexp

## Signature

```d
double qexp(double p, double rate = 1, bool lower = true);
double[] qexp(const(double)[] p, double rate = 1, bool lower = true);
GretlMatrix qexp(GretlMatrix p, double rate = 1, bool lower = true);      // also Submatrix
double[] qexp(Dataset ds, string series, double rate = 1, bool lower = true);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `p`: the probability (or probabilities). A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `qexp(ds, "series", ...)`. It should lie in [0, 1].
- `rate`: the rate. It must be positive and finite. The mean of the distribution is `1 / rate`. (The random generator `rexp` takes the mean instead.) The default is `1`.
- `lower`: if `true` (the default) `p` is a lower-tail probability, P(X <= x) = p; if `false` it is an upper-tail probability, P(X > x) = p.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `qexp(p, rate: values)` is allowed.

## Return value

The quantile function of the exponential distribution: the `x` with P(X <= x) = p (with `lower: false`, the `x` with P(X > x) = p), which inverts `pexp`. A `p` outside [0, 1] gives NaN. `p = 0` gives the lower end of the support and `p = 1` the upper end (with `lower: false` the other way round), so an unbounded support gives `-double.infinity` or `double.infinity`.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    const x = qexp(0.95, 0.5);
    writefln("%.6f", x);
    writefln("%.6f", pexp(x, 0.5));                 // P(X <= x) at that point
    writefln("%.6f", qexp(0.05, 0.5, lower: false));  // the point with 5% in the upper tail
    auto ps = [0.1, 0.6, 0.9];
    writefln("%(%.6f %)", qexp(ps, 0.5));         // an array
    writefln("%(%.6f %)", qexp(ps, [0.5, 1.0, 2.0]));  // the parameter varies by element
}
```

Output:

```
5.991465
0.950000
5.991465
0.210721 1.832581 4.605170
0.210721 0.916291 1.151293
```
