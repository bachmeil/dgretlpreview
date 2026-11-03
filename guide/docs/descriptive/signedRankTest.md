# signedRankTest

## Signature

```d
SignedRankTest signedRankTest(S, T)(S x, T y);
SignedRankTest signedRankTest(Dataset d, string first, string second);

struct SignedRankTest {
    size_t pairs, n, zeroDifferences, nonZeroTies;
    double wPlus, wMinus, expected, variance;
    double statistic, pvalue, pvalueOneSided;
    string report;
}
```

## Parameters

- `x`, `y`: two paired series of the same length (arrays or a row or column of a matrix), or in the `Dataset` form the names of two series. Different lengths throw a `GretlException`. Pairs where either value is missing are left out.

## Return value

A `SignedRankTest` struct, the Wilcoxon signed-rank test of the null hypothesis that the median of the differences `x - y` is zero, using Cureton's normal approximation with a continuity correction, with ties and zero differences handled (zero differences take the lowest ranks).

- `pairs`: the pairs with valid values; `n`: those with a non-zero difference; `zeroDifferences`: the rest; `nonZeroTies`: the groups of equal absolute differences.
- `wPlus`, `wMinus`: the sums of the ranks of the positive and of the negative differences.
- `expected`, `variance`: of `wPlus` under the null.
- `statistic`: the z statistic, `pvalue` the two-sided p-value and `pvalueOneSided` the tail beyond z (the one the library records).
- `report`: the text the library prints.

When `n` is 8 or less the normal approximation is not used: `expected`, `variance`, `statistic` and the p-values are NaN, and `report` gives 5% critical values for `n` from 6 to 8. If every difference is zero, `n` is 0 and the report says that the null is not rejected.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto x = [8.2, 7.9, 9.1, 8.8, 7.5, 9.4, 8.0, 8.9, 9.6, 7.7, 8.4, 9.0];
    auto y = [7.9, 7.8, 8.5, 8.9, 7.1, 8.8, 8.1, 8.2, 9.0, 7.2, 8.4, 8.6];

    auto t = signedRankTest(x, y);
    writefln("%d pairs, %d non-zero differences, %d zero", t.pairs, t.n, t.zeroDifferences);
    writefln("W+ = %.1f, W- = %.1f", t.wPlus, t.wMinus);
    writefln("z = %.3f, two-sided p = %.4f", t.statistic, t.pvalue);
}
```

Output:

```
12 pairs, 11 non-zero differences, 1 zero
W+ = 72.0, W- = 5.0
z = 2.592, two-sided p = 0.0095
```
