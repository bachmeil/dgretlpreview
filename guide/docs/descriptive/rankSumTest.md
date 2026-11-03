# rankSumTest

## Signature

```d
RankSumTest rankSumTest(S, T)(S x, T y);
RankSumTest rankSumTest(Dataset d, string first, string second);

struct RankSumTest {
    size_t nx, ny;
    double rankSumX, rankSumY, rankSumSmaller;
    double statistic, pvalue, pvalueOneSided;
    RankSumCritical critical;
    string report;
}

struct RankSumCritical {
    bool available;
    double lower1, lower5, lower10, upper10, upper5, upper1;
}
```

## Parameters

- `x`, `y`: two independent samples (arrays or a row or column of a matrix), or in the `Dataset` form the names of two series. They may have different lengths; missing values are left out and each sample needs at least one valid value.

## Return value

A `RankSumTest` struct, the Wilcoxon rank-sum (Mann-Whitney) test of the null hypothesis that the two samples come from populations with the same median:

- `nx`, `ny`: the numbers of valid observations.
- `rankSumX`, `rankSumY`: the sums of the (mid)ranks of each sample in the pooled ranking; they add to `(nx + ny)(nx + ny + 1) / 2`.
- `rankSumSmaller`: the rank sum of the smaller sample (of `x` if the sizes are equal).
- `statistic`: the z statistic for the smaller sample's rank sum; `pvalue`: the two-sided p-value; `pvalueOneSided`: the tail beyond z (the one the library records). These are NaN unless both samples have at least 10 observations, when the normal approximation is used.
- `critical`: for small samples (the smaller of 4 to 9 observations and the larger at most 12), the critical values of the library's table for `rankSumSmaller`: reject at the 1%, 5% or 10% level if it is at most `lower1`, `lower5`, `lower10` or at least `upper1`, `upper5`, `upper10` (`available` is `true` if the table has an entry). Below that the samples are too small for a test.
- `report`: the text the library prints.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(7);
    auto x = rnorm(25, 0.0);
    auto y = rnorm(30, 0.8);

    auto t = rankSumTest(x, y);
    writefln("rank sums %.1f and %.1f (they add to %.0f)", t.rankSumX, t.rankSumY, t.rankSumX + t.rankSumY);
    writefln("z = %.3f, two-sided p = %.5f", t.statistic, t.pvalue);

    auto small = rankSumTest([1.1, 2.3, 0.7, 3.1, 2.8], [4.2, 3.9, 5.1, 4.8, 3.3, 6.0]);    // too small for a p-value
    writefln("small samples: p = %s, critical values available: %s", small.pvalue, small.critical.available);
    writefln("rank sum of the smaller sample %.0f, 5%% limits %.0f and %.0f", small.rankSumSmaller, small.critical.lower5, small.critical.upper5);
}
```

Output:

```
rank sums 571.0 and 969.0 (they add to 1540)
z = -2.180, two-sided p = 0.02922
small samples: p = nan, critical values available: true
rank sum of the smaller sample 15, 5% limits 20 and 40
```
