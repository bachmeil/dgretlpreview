# summary

## Signature

```d
SummaryStats summary(X)(X x, string[] names = null, const(double)[] weights = null);
SummaryStats summary(Dataset d, string[] names, const(double)[] weights = null);

struct SummaryStats {
    string[] names;
    GretlMatrix valid, missing, mean, median, sd, min, max, cv, skewness, excessKurtosis, p5, p95, iqr;   // each 1 x k
    string toString();
}
```

## Parameters

- `x`: a series (an array, or a row or column of a matrix) or a matrix whose columns are the series.
- `names`: the names of the series, one for each column. The default is `x` for a single series and `x1`, `x2`, ... for several. A wrong number of names throws a `GretlException`.
- `weights`: one weight for each observation, used for every column. The statistics are then weighted, and observations whose weight is 0 or NaN are left out. The default is no weights.
- `d`: a `Dataset`, with `names` the series to summarize (required in this form).

## Return value

A `SummaryStats` struct. Each statistic is a 1 x k `GretlMatrix`, with one element for each series in the order given; a series with no valid observation has NaN in every field:

- `valid`: the number of valid observations (with weights, those with a non-zero weight).
- `missing`: the number of missing observations.
- `mean`, `median`, `min`, `max`.
- `sd`: the standard deviation (divisor `n - 1`; weighted, the library's weighted formula).
- `cv`: the coefficient of variation `|sd / mean|` (NaN if the mean is 0).
- `skewness` and `excessKurtosis`.
- `p5` and `p95`: the 5th and 95th percentiles.
- `iqr`: the interquartile range, the 0.75 quantile minus the 0.25 quantile.

`names` holds the series names, and `toString()` formats the table, so `writeln(s)` prints it. Less than one series or one observation throws a `GretlException`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto x = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
    auto y = [1.0, 3.0, 2.0, 5.0, 4.0, 6.0, double.nan, 8.0];

    writeln(summary(x));                                    // a single series is called x

    auto s = summary(cbind(x, y), names: ["x", "y"]);
    writeln(s);
    writefln("y has %.0f valid values and %.0f missing", s.valid[0, 1], s.missing[0, 1]);
    writefln("the mean of x is %.2f", s.mean[0, 0]);

    auto w = summary(x, weights: [1.0, 1, 1, 1, 2, 2, 0, 1]);    // the observation with weight 0 is left out
    writeln(w);
}
```

Output:

```
                n       mean     median        min        max       s.d.       skew    ex.kurt
x               8          5        4.5          2          9      2.138     0.6562    -0.2188

                n       mean     median        min        max       s.d.       skew    ex.kurt
x               8          5        4.5          2          9      2.138     0.6562    -0.2188
y               7      4.143          4          1          8       2.41     0.2834    -0.9598

y has 7 valid values and 1 missing
the mean of x is 5.00
                n       mean     median        min        max       s.d.       skew    ex.kurt
x               7      4.778          4          2          9       1.89      1.088      1.486
```
