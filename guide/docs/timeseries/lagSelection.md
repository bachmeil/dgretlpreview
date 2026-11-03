# lagSelection

## Signature

```d
LagSelection lagSelection(Y)(Y y, int maxLag, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool trend = false, bool seasonals = false,
                             int frequency = 1, string start = "1", string fitStart = null, string fitEnd = null);
LagSelection lagSelection(Dataset d, string[] endog, int maxLag, string[] exog = null, bool constant = true, bool trend = false, bool seasonals = false,
                          int frequency = 0, string start = null, string fitStart = null, string fitEnd = null);

struct LagSelection {
    GretlMatrix table;
    string[] columns;
    int bestAic, bestBic, bestHqc;
}
```

## Parameters

- `y`: the data of the endogenous series as the columns of a `GretlMatrix` or `Submatrix` (n rows, one column for each series), or an array for a single series. Missing values are not allowed.
- `maxLag`: the largest lag order to try, at least 2.
- `exog`, `constant`, `trend`, `seasonals`: the model, as for `var`.
- `frequency`, `start`: the calendar of the rows: the observations per year (1, 4, 12, ...) and the first observation ("1990:1" or "1990Q1"). The default is annual data from 1. They matter for seasonal dummies, the labels of the sample and the date arguments `fitStart` and `fitEnd`.
- `fitStart`, `fitEnd`: the estimation sample as observations of that calendar (for example "1960:2"). The default is all the rows, of which the first `order` are used as presample values.

## Return value

A `LagSelection` struct. VARs of order 1 to `maxLag` are estimated on the same sample (the one of the longest), and `table` has one row for each lag order with the columns named in `columns`: `lags`, `loglik` (the log-likelihood), `p(LR)` (the p-value of the likelihood-ratio test of each lag order against the one before), `AIC`, `BIC` and `HQC` (the information criteria; smaller is better). `bestAic`, `bestBic` and `bestHqc` are the lag orders chosen by each criterion. In the `Dataset` forms the series are named: `endog` names the endogenous series, `exog` the exogenous ones, and the calendar of the `Dataset` is used if it has one (`frequency` 0 and `start` null take it from the data).

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%10.4f", m[i, j]);
        }
        writeln();
    }
}

void main() {
    setSeed(31);
    auto n = 250;
    auto e1 = rnorm(n);
    auto e2 = rnorm(n);
    auto y1 = new double[n];
    auto y2 = new double[n];
    y1[0] = e1[0];
    y2[0] = e2[0];
    foreach (t; 1 .. n) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + e1[t];
        y2[t] = 0.3 * y1[t - 1] + 0.4 * y2[t - 1] + e2[t];
    }

    auto ls = lagSelection(cbind(y1, y2), 4);
    writeln(ls.columns);
    show(ls.table);
    writefln("best order: AIC %d, BIC %d, HQC %d", ls.bestAic, ls.bestBic, ls.bestHqc);
}
```

Output:

```
["lags", "loglik", "p(LR)", "AIC", "BIC", "HQC"]
    1.0000 -503.6965    0.0000    4.1439    4.2294    4.1783
    2.0000 -502.2393    0.5722    4.1645    4.3070    4.2219
    3.0000 -501.4757    0.8218    4.1909    4.3903    4.2712
    4.0000 -500.3403    0.6861    4.2141    4.4706    4.3174
best order: AIC 1, BIC 1, HQC 1
```
