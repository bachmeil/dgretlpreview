# dlogis

## Signature

```d
double dlogis(double x, double location = 0, double scale = 1);
double[] dlogis(const(double)[] x, double location = 0, double scale = 1);
GretlMatrix dlogis(GretlMatrix x, double location = 0, double scale = 1);      // also Submatrix
double[] dlogis(Dataset ds, string series, double location = 0, double scale = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dlogis(ds, "series", ...)`.
- `location`: the mean and median. It must be finite. The default is `0`.
- `scale`: the scale parameter. It must be positive and finite. The standard deviation is `scale * pi / sqrt(3)`. The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dlogis(x, scale: values)` is allowed.

## Return value

The density of the logistic distribution at each point: `exp(-z) / (scale (1 + exp(-z))^2)` with `z = (x - location) / scale`. It is positive everywhere.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dlogis(2.0, 1, 2));
    auto xs = [-1.0, 1.0, 3.0];
    writefln("%(%.6f %)", dlogis(xs, 1, 2));         // an array: one value per element
    writefln("%(%.6f %)", dlogis(xs, scale: [1.0, 2, 3]));  // the parameter varies by element
}
```

Output:

```
0.117502
0.098306 0.125000 0.098306
0.196612 0.117502 0.065537
```
