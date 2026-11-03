# dtnorm

## Signature

```d
double dtnorm(double x, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity);
double[] dtnorm(const(double)[] x, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity);
GretlMatrix dtnorm(GretlMatrix x, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity);      // also Submatrix
double[] dtnorm(Dataset ds, string series, double mean = 0, double sd = 1, double lo = -double.infinity, double hi = double.infinity);
```

The `Dataset` form takes the name of one of its series; the others take the data themselves. The first three are templates, so any numeric array or `Submatrix` also works.

## Parameters

- `x`: where to evaluate the density. A number, a numeric array, a `GretlMatrix` or `Submatrix`, or (in the `Dataset` form) the name of one of its series, given with the `Dataset` as `dtnorm(ds, "series", ...)`.
- `mean`: the mean of the normal distribution before truncation. It must be finite. The default is `0`.
- `sd`: the standard deviation of the normal distribution before truncation. It must be positive and finite. The default is `1`.
- `lo`: the lower bound of the truncation. The default means no lower bound. The default is `-double.infinity`.
- `hi`: the upper bound of the truncation. It must be greater than `lo`. The default means no upper bound. The default is `double.infinity`.
- Any parameter of the distribution may also be an array, a matrix, a submatrix or (with a `Dataset`) a series name instead of a number, with one value per element of the first argument; a number is repeated. Parameters can be given by name, and those that are not given take their defaults, so `dtnorm(x, hi: values)` is allowed.

## Return value

The density of the truncated normal distribution at each point: the normal density with the given `mean` and `sd`, divided by the probability that a normal draw falls in [`lo`, `hi`], on that interval. The support is lo <= x <= hi; the density is 0 outside it.

The result has the shape of the first argument: a `double` for a number, a `double[]` for an array, a `GretlMatrix` of the same dimensions for a `GretlMatrix` or `Submatrix`, and a `double[]` with one value per observation for a series of a `Dataset`. If the first argument is a number but a parameter is an array or matrix, the result is a `double[]` (or `GretlMatrix`) with one value per element of the parameter. A missing value (NaN) in the first argument gives NaN, and so does a NaN in an array parameter (for that element). A `GretlException` is thrown if a parameter is invalid (the message says which), if a single-number parameter is NaN, or if arrays or matrices of different sizes are combined.

`mean` and `sd` describe the normal distribution before truncation, not the mean and standard deviation of the truncated distribution. Parameters can be skipped and given by name, as in `ptnorm(x, hi: 2)` (a normal distribution truncated above at 2).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    writefln("%.6f", dtnorm(1.0, 1, 2, 0, 4));
    auto xs = [0.5, 1.0, 3.0];
    writefln("%(%.6f %)", dtnorm(xs, 1, 2, 0, 4));         // an array: one value per element
    writefln("%(%.6f %)", dtnorm(xs, mean: [0.0, 1, 2]));  // the parameter varies by element
}
```

Output:

```
0.319330
0.309505 0.319330 0.193683
0.352065 0.398942 0.241971
```
