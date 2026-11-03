# fiveNumbers

## Signature

```d
FiveNumbers fiveNumbers(S)(S x);

struct FiveNumbers {
    double min, q1, median, q3, max;
}
```

## Parameters

- `x`: a series: an array of numbers (`double[]`, `int[]`, ...) or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). Missing values (NaN) are skipped.

## Return value

A `FiveNumbers` struct with the five numbers: `min`, the first quartile `q1` (the 0.25 quantile), `median`, the third quartile `q3` (the 0.75 quantile), and `max`, with the default quantile rule. A quartile that the sample is too small for is NaN. Only series are accepted, not matrices.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto x = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    auto f = fiveNumbers(x);
    writefln("min %.2f, q1 %.2f, median %.2f, q3 %.2f, max %.2f", f.min, f.q1, f.median, f.q3, f.max);
}
```

Output:

```
min 2.00, q1 4.00, median 4.50, q3 6.50, max 9.00
```
