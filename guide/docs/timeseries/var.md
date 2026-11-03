# var

## Signature

```d
VarFit var(Y)(Y y, int order, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool trend = false, bool seasonals = false,
              int[] lags = null, bool robust = false, string[] names = null, string[] exogNames = null, int frequency = 1,
              string start = "1", string fitStart = null, string fitEnd = null);
VarFit var(Dataset d, string[] endog, int order, string[] exog = null, bool constant = true, bool trend = false, bool seasonals = false,
           int[] lags = null, bool robust = false, int frequency = 0, string start = null, string fitStart = null, string fitEnd = null);
```

## Parameters

- `y`: the data of the endogenous series as the columns of a `GretlMatrix` or `Submatrix` (n rows, one column for each series), or an array for a single series. Missing values are not allowed.
- `order`: the lag order, at least 1: the VAR has lags 1 to `order` of every series.
- `exog`: a `GretlMatrix` of exogenous regressors with the same number of rows as `y` (in the `Dataset` form, `exog` is an array of series names). The default is none.
- `constant`: include a constant. The default is `true`.
- `trend`: add a linear trend. The default is `false`.
- `seasonals`: add seasonal dummies; this needs a frequency above 1. The default is `false`.
- `lags`: a gappy list of lags such as `[1, 4]`, in which case `order` must be the largest. The default (null) uses every lag from 1 to `order`.
- `robust`: use heteroskedasticity-consistent standard errors. The default is `false`.
- `names`, `exogNames`: the names of the columns of `y` and of `exog`. The default is `y1`, `y2`, ... and `x1`, `x2`, ....
- `frequency`, `start`: the calendar of the rows: the observations per year (1, 4, 12, ...) and the first observation ("1990:1" or "1990Q1"). The default is annual data from 1. They matter for seasonal dummies, the labels of the sample and the date arguments `fitStart` and `fitEnd`.
- `fitStart`, `fitEnd`: the estimation sample as observations of that calendar (for example "1960:2"). The default is all the rows, of which the first `order` are used as presample values.

## Return value

A `VarFit` struct, the estimated vector autoregression: the coefficients, standard errors, residual covariance matrix, residuals, roots of the companion matrix, Granger-causality F tests and information criteria (see the page of `VarFit` for the fields). It also keeps the library's estimate, so that `irf`, `fevd` and `forecast` can be applied to it, and `equation(i)` gives the regression of equation `i` as an `Estimate`. Copies of a `VarFit` share that data, which is freed when the last copy goes away. In the `Dataset` forms the series are named: `endog` names the endogenous series, `exog` the exogenous ones, and the calendar of the `Dataset` is used if it has one (`frequency` 0 and `start` null take it from the data). The series must have no missing values in the sample.

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
    auto n = 200;
    auto e1 = rnorm(n);
    auto e2 = rnorm(n);
    auto y1 = new double[n];                         // a bivariate VAR(1)
    auto y2 = new double[n];
    y1[0] = e1[0];
    y2[0] = e2[0];
    foreach (t; 1 .. n) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + e1[t];
        y2[t] = 0.3 * y1[t - 1] + 0.4 * y2[t - 1] + e2[t];
    }

    auto fit = var(cbind(y1, y2), 1);
    writeln(fit.coefNames);
    show(fit.coef);                                  // one column for each equation
    writefln("%d observations, sample %s to %s", fit.nobs, fit.sampleStart, fit.sampleEnd);

    auto named = var(cbind(y1, y2), 2, names: ["a", "b"]);     // a VAR(2) with named series
    writeln(named.names);
    writefln("AIC of the VAR(1) %.4f, of the VAR(2) %.4f", fit.aic, named.aic);
}
```

Output:

```
["const", "y1_1", "y2_1"]
   -0.1326    0.0543
    0.5427    0.2739
    0.1903    0.3882
199 observations, sample 2 to 200
["a", "b"]
AIC of the VAR(1) 4.1781, of the VAR(2) 4.1919
```
