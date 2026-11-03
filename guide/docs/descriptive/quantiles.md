# quantiles

## Signature

```d
double[] quantiles(S)(S x, const(double)[] ps);
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.
- `ps`: the probabilities, each from 0 to 1 (a value outside that range throws a `GretlException`).

## Return value

A `double[]` with one quantile for each element of `ps`, in the same order, each as `quantile(x, p)` with the default rule (`QuantileType.q6`). A quantile the sample is too small for is NaN. Only series are accepted, not matrices.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto x = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    auto q = quantiles(x, [0.2, 0.25, 0.5, 0.75, 0.8]);
    writefln("%(%.3f %)", q);
}
```

Output:

```
3.600 4.000 4.500 6.500 7.400
```
