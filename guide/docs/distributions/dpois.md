# dpois

## Signature

```d
double dpois(double x, double lambda);
double[] dpois(const(double)[] x, double lambda);
GretlMatrix dpois(GretlMatrix x, double lambda);      // also Submatrix
double[] dpois(Dataset ds, string series, double lambda);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the probability mass. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dpois(ds, "series", ...)`.
- `lambda`: the mean. It must be positive and finite. It has no default.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The probability mass of the Poisson distribution at each point: the probability mass `exp(-lambda) lambda^x / x!` at a whole number `x`. The support is the whole numbers 0, 1, 2, ...; the mass is 0 outside it and at any `x` that is not a whole number.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dpois(4.0, 4.5));
    auto xs = [2.0, 4.0, 6.0];
    writefln("%(%.6f %)", dpois(xs, 4.5));         // an array: one value per element
    writefln("%(%.6f %)", dpois(xs, [2.0, 4.5, 8.0]));  // the parameter varies by element
    writefln("%.6f", dpois(4.0 + 0.5, 4.5));         // not a whole number: probability 0
}
```

Output:

```
0.189808
0.112479 0.189808 0.128120
0.270671 0.189808 0.122138
0.000000
```
