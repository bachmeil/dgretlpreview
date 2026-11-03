# dbeta

## Signature

```d
double dbeta(double x, double shape1, double shape2);
double[] dbeta(const(double)[] x, double shape1, double shape2);
GretlMatrix dbeta(GretlMatrix x, double shape1, double shape2);      // also Submatrix
double[] dbeta(Dataset ds, string series, double shape1, double shape2);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dbeta(ds, "series", ...)`.
- `shape1`: the first shape parameter. It must be positive and finite. It has no default.
- `shape2`: the second shape parameter. It must be positive and finite. It has no default.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name.

## Return value

The density of the beta distribution at each point: `x^(shape1 - 1) (1 - x)^(shape2 - 1) / B(shape1, shape2)`. The support is 0 <= x <= 1; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dbeta(0.25, 2, 6));
    auto xs = [0.1, 0.25, 0.5];
    writefln("%(%.6f %)", dbeta(xs, 2, 6));         // an array: one value per element
    writefln("%(%.6f %)", dbeta(xs, 2, [4.0, 6, 8]));  // the parameter varies by element
}
```

Output:

```
2.491699
2.480058 2.491699 0.656250
1.458000 2.491699 0.281250
```
