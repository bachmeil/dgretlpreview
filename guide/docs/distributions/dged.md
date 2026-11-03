# dged

## Signature

```d
double dged(double x, double shape);
double[] dged(const(double)[] x, double shape);
GretlMatrix dged(GretlMatrix x, double shape);      // also Submatrix
double[] dged(Dataset ds, string series, double shape);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dged(ds, "series", ...)`.
- `shape`: the shape parameter, in the library's parameterization. It must be positive and finite. A shape of 2 is the standard normal distribution; smaller values give heavier tails. It has no default.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The density of the generalized error distribution at each point: the density of the generalized error distribution, centred at 0, as the library defines it (the standard normal density when `shape = 2`). It is positive everywhere.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dged(0.5, 1.5));
    auto xs = [-1.0, 0.0, 1.0];
    writefln("%(%.6f %)", dged(xs, 1.5));         // an array: one value per element
    writefln("%(%.6f %)", dged(xs, [1.0, 1.5, 2.0]));  // the parameter varies by element
}
```

Output:

```
0.359134
0.214587 0.475967 0.214587
0.171909 0.475967 0.241971
```
