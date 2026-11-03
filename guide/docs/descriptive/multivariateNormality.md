# multivariateNormality

## Signature

```d
StatTest multivariateNormality(E)(E e, GretlMatrix sigma = GretlMatrix.init);

struct StatTest {
    double statistic, pvalue;
}
```

## Parameters

- `e`: a `GretlMatrix` or `Submatrix` with n rows and p columns, for example the residuals of a system. It needs at least two columns and more than `p + 8` rows, and no missing values.
- `sigma`: the p x p covariance matrix of the columns. The default (an empty matrix) uses `e'e / n`.

## Return value

A `StatTest` with the Doornik-Hansen statistic of the null hypothesis that the rows of `e` are drawn from a multivariate normal distribution, which is chi-square with `2 p` degrees of freedom, and its p-value. A small p-value is evidence against multivariate normality. Too few rows or columns, or a wrong-sized `sigma`, throws a `GretlException`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(4);
    auto e = randomMatrix!rnorm(200, 3);          // 200 observations of three independent normals
    auto t = multivariateNormality(e);
    writefln("statistic %.3f, p-value %.3f", t.statistic, t.pvalue);
}
```

Output:

```
statistic 10.068, p-value 0.122
```
