# dt

## Signature

```d
double dt(double x, double df, double ncp = 0);
double[] dt(const(double)[] x, double df, double ncp = 0);
GretlMatrix dt(GretlMatrix x, double df, double ncp = 0);      // also Submatrix
double[] dt(Dataset ds, string series, double df, double ncp = 0);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dt(ds, "series", ...)`.
- `df`: the degrees of freedom. It must be positive and finite, and need not be an integer. It has no default.
- `ncp`: the non-centrality parameter. It must be finite. The default 0 is the ordinary (central) t distribution; any other value gives the non-central t distribution. The default is `0`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dt(x, ncp: values)` is allowed.

## Return value

The density of the Student's t distribution at each point: `gamma((df + 1) / 2) / (sqrt(df pi) gamma(df / 2)) * (1 + x^2 / df)^(-(df + 1) / 2)` when `ncp = 0` (the density of the non-central t otherwise). It is positive everywhere.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dt(2.0, 10));
    auto xs = [-2.0, 0.0, 2.0];
    writefln("%(%.6f %)", dt(xs, 10));         // an array: one value per element
    writefln("%(%.6f %)", dt(xs, [5.0, 10, 20]));  // the parameter varies by element
    writefln("%.6f", dt(2.0, 10, ncp: 1.5));       // non-central
}
```

Output:

```
0.061146
0.061146 0.389108 0.061146
0.065090 0.389108 0.058087
0.314606
```
