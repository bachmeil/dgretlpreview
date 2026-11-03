# dnorm

## Signature

```d
double dnorm(double x, double mean = 0, double sd = 1);
double[] dnorm(const(double)[] x, double mean = 0, double sd = 1);
GretlMatrix dnorm(GretlMatrix x, double mean = 0, double sd = 1);      // also Submatrix
double[] dnorm(Dataset ds, string series, double mean = 0, double sd = 1);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dnorm(ds, "series", ...)`.
- `mean`: the mean. It must be finite. The default is `0`.
- `sd`: the standard deviation. It must be positive and finite. The default is `1`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dnorm(x, sd: values)` is allowed.

## Return value

The density of the normal distribution at each point: `exp(-(x - mean)^2 / (2 sd^2)) / (sd sqrt(2 pi))`. It is positive everywhere.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dnorm(0.5, 1, 2));
    auto xs = [-1.0, 0.0, 1.0];
    writefln("%(%.6f %)", dnorm(xs, 1, 2));         // an array: one value per element
    writefln("%(%.6f %)", dnorm(xs, sd: [1.0, 2, 3]));  // the parameter varies by element
}
```

Output:

```
0.193334
0.120985 0.176033 0.199471
0.241971 0.199471 0.125794
```
