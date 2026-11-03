# crosstab

## Signature

```d
CrossTab crosstab(S, T)(S x, T y);
CrossTab crosstab(Dataset d, string rows, string cols);

struct CrossTab {
    GretlMatrix counts, rowValues, colValues, rowTotals, colTotals;
    int nobs, df;
    double pearson, pvalue;
}
```

## Parameters

- `x`, `y`: two discrete series of the same length (integer valued, or taking few distinct values), as arrays or a row or column of a matrix; in the `Dataset` form, the names of the series for the rows and columns. Different lengths, or fewer than 2 observations, throw a `GretlException`. Observations with a missing value in either series are left out.

## Return value

A `CrossTab` struct, the table of counts of pairs of values with Pearson's chi-square test of independence:

- `counts`: a rows x cols `GretlMatrix`; `counts[i, j]` is the number of observations with the i-th value of `x` and the j-th value of `y`.
- `rowValues` (rows x 1) and `colValues` (cols x 1): the distinct values of `x` and of `y`, increasing.
- `rowTotals` (rows x 1) and `colTotals` (cols x 1): the marginal counts.
- `nobs`: the number of observations counted.
- `pearson`: the statistic of independence, the sum of `(O - E)^2 / E`; `df`: `(rows - 1)(cols - 1)`; `pvalue`: the chi-square(`df`) p-value (NaN if `df` is 0).

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%9.4f", m[i, j]);
        }
        writeln();
    }
}

void main() {
    auto smoker = [0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0];
    auto illness = [0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0];

    auto t = crosstab(smoker, illness);
    show(t.counts);
    writefln("row totals %(%.0f %), column totals %(%.0f %)", Col(t.rowTotals, 0).dup, Col(t.colTotals, 0).dup);
    writefln("Pearson chi-square %.3f with %d df, p = %.5f", t.pearson, t.df, t.pvalue);
}
```

Output:

```
  10.0000   1.0000
   1.0000   8.0000
row totals 11 9, column totals 11 9
Pearson chi-square 12.735 with 1 df, p = 0.00036
```
