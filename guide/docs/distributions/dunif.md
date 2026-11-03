# dunif

## Signature

```d
double dunif(double x, double min = 0, double max = 1);
double[] dunif(const(double)[] x, double min = 0, double max = 1);
GretlMatrix dunif(GretlMatrix x, double min = 0, double max = 1);      // also Submatrix
double[] dunif(Dataset ds, string series, double min = 0, double max = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dunif(ds, "series", ...)`.
- `min`: the lower end of the interval. It must be finite. The default is `0`.
- `max`: the upper end of the interval. It must be finite and greater than `min`. The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dunif(x, max: values)` is allowed.

## Return value

The density of the uniform distribution at each point: `1 / (max - min)` on the interval. The support is min <= x <= max; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dunif(3.0, 2, 5));
    auto xs = [1.0, 3.0, 4.5];
    writefln("%(%.6f %)", dunif(xs, 2, 5));         // an array: one value per element
    writefln("%(%.6f %)", dunif(xs, max: [4.0, 5, 6]));  // the parameter varies by element
}
```

Output:

```
0.333333
0.000000 0.333333 0.333333
0.250000 0.200000 0.166667
```
