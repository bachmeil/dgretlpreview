# dbinom

## Signature

```d
double dbinom(double x, double size, double prob);
double[] dbinom(const(double)[] x, double size, double prob);
GretlMatrix dbinom(GretlMatrix x, double size, double prob);      // also Submatrix
double[] dbinom(Dataset ds, string series, double size, double prob);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the probability mass. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dbinom(ds, "series", ...)`.
- `size`: the number of trials, a whole number of at least 0. It has no default.
- `prob`: the probability of success on each trial, from 0 to 1. It has no default.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The probability mass of the binomial distribution at each point: the probability mass `choose(size, x) prob^x (1 - prob)^(size - x)` at a whole number `x` from 0 to `size`. The support is the whole numbers 0 to `size`; the mass is 0 outside it and at any `x` that is not a whole number.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dbinom(3.0, 10, 0.3));
    auto xs = [2.0, 3.0, 4.0];
    writefln("%(%.6f %)", dbinom(xs, 10, 0.3));         // an array: one value per element
    writefln("%(%.6f %)", dbinom(xs, 10, [0.2, 0.3, 0.4]));  // the parameter varies by element
    writefln("%.6f", dbinom(3.0 + 0.5, 10, 0.3));         // not a whole number: probability 0
}
```

Output:

```
0.266828
0.233474 0.266828 0.200121
0.301990 0.266828 0.250823
0.000000
```
