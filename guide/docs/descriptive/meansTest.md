# meansTest

## Signature

```d
TwoSampleTest meansTest(S, T)(S x, T y, MeansTest kind = MeansTest.equalVariances);
TwoSampleTest meansTest(Dataset d, string first, string second, MeansTest kind = MeansTest.equalVariances);

enum MeansTest { equalVariances, unequalVariances, paired }

struct TwoSampleTest {
    string test;
    double statistic, pvalue, df, df2;
    double difference, se;
    size_t n1, n2;
}
```

## Parameters

- `x`, `y`: the two samples (arrays or a row or column of a matrix), or in the `Dataset` form the names of two series. Missing values are left out, and each sample needs at least two valid values. The samples may have different lengths, except for `MeansTest.paired`, where they must have the same length and a pair with a missing value is left out.
- `kind`: `MeansTest.equalVariances` (the default: the pooled-variance t test), `MeansTest.unequalVariances` (Welch's t test, with the Satterthwaite degrees of freedom, which the library rounds down to an integer) or `MeansTest.paired` (the t test of the mean of the differences of paired observations).

## Return value

A `TwoSampleTest` struct, the test of the null hypothesis that the two populations have the same mean:

- `test`: the name of the test.
- `statistic`: the t statistic.
- `pvalue`: the two-sided p-value.
- `df`: the degrees of freedom (`n1 + n2 - 2` for pooled variances, `n1 - 1` for paired).
- `difference` and `se`: the difference of the sample means (for a paired test, the mean difference) and its standard error.
- `n1`, `n2`: the numbers of valid observations in each sample (for a paired test, the number of pairs in `n1` and 0 in `n2`).
- `df2`: NaN (it is used only by the F test, see `variancesTest`).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(5);
    auto x = rnorm(40, 0.0);
    auto y = rnorm(50, 0.5);                       // a larger mean

    auto t = meansTest(x, y);
    writefln("%s: t = %.3f, df = %.0f, p = %.4f", t.test, t.statistic, t.df, t.pvalue);
    writefln("difference %.3f, standard error %.3f", t.difference, t.se);

    auto w = meansTest(x, y, MeansTest.unequalVariances);
    writefln("%s: t = %.3f, df = %.0f, p = %.4f", w.test, w.statistic, w.df, w.pvalue);

    auto before = rnorm(30);
    auto after = before.dup;
    foreach (i, ref v; after) {
        v += 0.4 + 0.1 * (cast(int) (i % 3) - 1);
    }
    auto p = meansTest(before, after, MeansTest.paired);
    writefln("%s: t = %.3f, p = %.5f", p.test, p.statistic, p.pvalue);
}
```

Output:

```
Equality of means, equal variances: t = -3.287, df = 88, p = 0.0015
difference -0.557, standard error 0.169
Equality of means, unequal variances: t = -3.301, df = 84, p = 0.0014
Paired differences: t = -26.382, p = 0.00000
```
