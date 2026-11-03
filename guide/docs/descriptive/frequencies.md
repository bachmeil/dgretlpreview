# frequencies

## Signature

```d
FrequencyTable frequencies(S)(S x, int bins = 0, double min = double.nan, double width = double.nan, bool normalTest = false);

struct FrequencyTable {
    bool discrete;
    GretlMatrix midpoint, lower, upper, count, fraction, cumulative;   // each nbins x 1
    int nobs;
    double mean, sd;
    double normalTest, normalPvalue;
}
```

## Parameters

- `x`: a series: an array of numbers or one row or column of a matrix. Missing values are left out. At least 8 valid values are needed, and the series must not be constant, otherwise a `GretlException` is thrown.
- `bins`: the number of bins, at least 2. The default 0 chooses it from the sample size (an odd number from 5 to 29).
- `min` and `width`: the left end of the first bin and the bin width, which must be given together (and not together with `bins`); the number of bins follows.
- `normalTest`: if `true`, also compute the Doornik-Hansen test of normality (needs at least 8 values). The default is `false`.

## Return value

A `FrequencyTable` struct with one row per bin; each matrix field is nbins x 1:

- `midpoint`: the middle of each bin.
- `lower`, `upper`: the ends of each bin. The first and last bins also hold anything below or above the outer ends.
- `count`: the number of observations in each bin.
- `fraction`: `count / nobs`.
- `cumulative`: the cumulative fractions.
- `nobs`: the number of observations used, and `mean` and `sd`: the mean and standard deviation of the series.
- `normalTest`, `normalPvalue`: the chi-square(2) statistic and p-value of the Doornik-Hansen test, if asked for; otherwise NaN.
- `discrete`: `true` when the series takes few distinct values, in which case there is one row for each value (`midpoint`, `lower` and `upper` are then the value). `valueCounts` forces this form.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(1);
    auto x = rnorm(200);
    auto f = frequencies(x);
    writefln("%d bins, %d observations, mean %.3f", f.count.nrow, f.nobs, f.mean);
    foreach (i; 0 .. f.count.nrow) {
        writefln("[%6.2f, %6.2f)  %3.0f  %.3f", f.lower[i, 0], f.upper[i, 0], f.count[i, 0], f.fraction[i, 0]);
    }

    auto g = frequencies(x, bins: 5, normalTest: true);
    writefln("5 bins; normality statistic %.3f, p-value %.3f", g.normalTest, g.normalPvalue);

    auto h = frequencies(x, min: -3, width: 1);          // bins of width 1 starting at -3
    writefln("%d bins of width 1", h.count.nrow);
}
```

Output:

```
15 bins, 200 observations, mean 0.038
[ -2.22,  -1.88)    4  0.020
[ -1.88,  -1.53)    5  0.025
[ -1.53,  -1.19)   12  0.060
[ -1.19,  -0.84)   16  0.080
[ -0.84,  -0.50)   22  0.110
[ -0.50,  -0.15)   28  0.140
[ -0.15,   0.19)   24  0.120
[  0.19,   0.54)   29  0.145
[  0.54,   0.88)   24  0.120
[  0.88,   1.23)   15  0.075
[  1.23,   1.57)   11  0.055
[  1.57,   1.91)    5  0.025
[  1.91,   2.26)    4  0.020
[  2.26,   2.60)    0  0.000
[  2.60,   2.95)    1  0.005
5 bins; normality statistic 0.785, p-value 0.675
6 bins of width 1
```
