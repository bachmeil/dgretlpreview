# dexp

## Signature

```d
double dexp(double x, double rate = 1);
double[] dexp(const(double)[] x, double rate = 1);
GretlMatrix dexp(GretlMatrix x, double rate = 1);      // also Submatrix
double[] dexp(Dataset ds, string series, double rate = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dexp(ds, "series", ...)`.
- `rate`: the rate. It must be positive and finite. The mean of the distribution is `1 / rate`. (The random generator `rexp` takes the mean instead.) The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dexp(x, rate: values)` is allowed.

## Return value

The density of the exponential distribution at each point: `rate exp(-rate x)`. The support is x >= 0; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dexp(1.0, 0.5));
    auto xs = [0.5, 1.0, 2.0];
    writefln("%(%.6f %)", dexp(xs, 0.5));         // an array: one value per element
    writefln("%(%.6f %)", dexp(xs, [0.5, 1.0, 2.0]));  // the parameter varies by element
}
```

Output:

```
0.303265
0.389400 0.303265 0.183940
0.389400 0.367879 0.036631
```
