# dlaplace

## Signature

```d
double dlaplace(double x, double mean = 0, double scale = 1);
double[] dlaplace(const(double)[] x, double mean = 0, double scale = 1);
GretlMatrix dlaplace(GretlMatrix x, double mean = 0, double scale = 1);      // also Submatrix
double[] dlaplace(Dataset ds, string series, double mean = 0, double scale = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dlaplace(ds, "series", ...)`.
- `mean`: the mean and median. It must be finite. The default is `0`.
- `scale`: the scale parameter. It must be positive and finite. The variance is `2 scale^2`. The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dlaplace(x, scale: values)` is allowed.

## Return value

The density of the Laplace distribution at each point: `exp(-|x - mean| / scale) / (2 scale)`. It is positive everywhere.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dlaplace(2.0, 1, 2));
    auto xs = [-1.0, 1.0, 3.0];
    writefln("%(%.6f %)", dlaplace(xs, 1, 2));         // an array: one value per element
    writefln("%(%.6f %)", dlaplace(xs, scale: [1.0, 2, 3]));  // the parameter varies by element
}
```

Output:

```
0.151633
0.091970 0.250000 0.091970
0.183940 0.151633 0.061313
```
