# normalityTests

## Signature

```d
NormalityTests normalityTests(S)(S x);

struct NormalityTests {
    StatTest doornikHansen, shapiroWilk, jarqueBera, lilliefors;
}
struct StatTest {
    double statistic, pvalue;
}
```

## Parameters

- `x`: a series: an array of numbers or one row or column of a matrix. Missing values are left out. At least 8 valid values are needed and the series must not be constant, otherwise a `GretlException` is thrown.

## Return value

A `NormalityTests` struct with the four tests of the null hypothesis that the series is normally distributed, each a `StatTest` with a `statistic` and a `pvalue` (NaN for a test that cannot be computed for the sample). A small p-value is evidence against normality.

- `doornikHansen`: the Doornik-Hansen test, chi-square(2), from the skewness and kurtosis.
- `shapiroWilk`: the Shapiro-Wilk W (the p-value is approximate).
- `jarqueBera`: the Jarque-Bera test, chi-square(2).
- `lilliefors`: the Kolmogorov-Smirnov test with estimated mean and variance (the p-value is approximate).

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(2);
    auto normal = rnorm(100);
    auto n = normalityTests(normal);
    writefln("Doornik-Hansen  %.3f  p = %.3f", n.doornikHansen.statistic, n.doornikHansen.pvalue);
    writefln("Shapiro-Wilk    %.3f  p = %.3f", n.shapiroWilk.statistic, n.shapiroWilk.pvalue);
    writefln("Jarque-Bera     %.3f  p = %.3f", n.jarqueBera.statistic, n.jarqueBera.pvalue);
    writefln("Lilliefors      %.3f  p = %.3f", n.lilliefors.statistic, n.lilliefors.pvalue);

    auto skewed = rexp(100);                       // exponential draws are far from normal
    writefln("Doornik-Hansen for exponential data: p = %.6f", normalityTests(skewed).doornikHansen.pvalue);
}
```

Output:

```
Doornik-Hansen  3.985  p = 0.136
Shapiro-Wilk    0.970  p = 0.021
Jarque-Bera     3.595  p = 0.166
Lilliefors      0.098  p = 0.020
Doornik-Hansen for exponential data: p = 0.000000
```
