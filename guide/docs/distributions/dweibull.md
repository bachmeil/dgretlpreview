# dweibull

## Signature

```d
double dweibull(double x, double shape, double scale = 1);
double[] dweibull(const(double)[] x, double shape, double scale = 1);
GretlMatrix dweibull(GretlMatrix x, double shape, double scale = 1);      // also Submatrix
double[] dweibull(Dataset ds, string series, double shape, double scale = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dweibull(ds, "series", ...)`.
- `shape`: the shape parameter. It must be positive and finite. It has no default.
- `scale`: the scale parameter. It must be positive and finite. The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dweibull(x, scale: values)` is allowed.

## Return value

The density of the Weibull distribution at each point: `(shape / scale) (x / scale)^(shape - 1) exp(-(x / scale)^shape)`. The support is x >= 0; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dweibull(1.0, 1.5, 2));
    auto xs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", dweibull(xs, 1.5, 2));         // an array: one value per element
    writefln("%(%.6f %)", dweibull(xs, 1.5, [1.0, 2, 3]));  // the parameter varies by element
}
```

Output:

```
0.372392
0.330936 0.372392 0.146304
0.744783 0.372392 0.183940
```
