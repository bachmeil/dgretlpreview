# variancesTest

## Signature

```d
VarianceTest variancesTest(S, T)(S x, T y);
VarianceTest variancesTest(Dataset d, string first, string second);

struct VarianceTest {
    TwoSampleTest test;
    double pvalueUpper;
}
```

## Parameters

- `x`, `y`: the two samples (arrays or a row or column of a matrix), or in the `Dataset` form the names of two series. Missing values are left out and each sample needs at least two valid values. The samples may have different lengths.

## Return value

A `VarianceTest` struct, the F test of the null hypothesis that the two populations have the same variance: `F = var(x) / var(y)` with `n1 - 1` and `n2 - 1` degrees of freedom.

- `test.statistic`: F. `test.df` and `test.df2`: the numerator and denominator degrees of freedom. `test.n1`, `test.n2`: the sample sizes. `test.test`: the name of the test.
- `test.pvalue`: the two-sided p-value, `2 min(P(F' > F), P(F' < F))`.
- `pvalueUpper`: the one-sided p-value `P(F' > F)`, which is the one the library itself records (see `bugs.md`, item 124).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(6);
    auto x = rnorm(60, 0.0, 1.0);
    auto y = rnorm(50, 0.0, 2.0);                  // a larger variance

    auto v = variancesTest(x, y);
    writefln("F = %.4f with (%.0f, %.0f) df", v.test.statistic, v.test.df, v.test.df2);
    writefln("two-sided p = %.6f, upper p = %.6f", v.test.pvalue, v.pvalueUpper);
}
```

Output:

```
F = 0.3182 with (59, 49) df
two-sided p = 0.000034, upper p = 0.999983
```
