# ols

## Signature

```d
Estimate ols(Y, X)(Y y, X x, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
                   bool prewhiten = false, const(double)[] cluster = null);
Estimate ols(Dataset d, string y, string[] x, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
             bool prewhiten = false, string cluster = null);

enum Vcov { classical, hc0, hc1, hc2, hc3, hac }
enum HACKernel { bartlett, parzen, qs }
```

## Parameters

- `y`: the dependent variable, an array of numbers or a `GretlMatrix` or `Submatrix` with one column (n values).
- `x`: the regressors: a `GretlMatrix` or `Submatrix` with n rows and one column for each regressor, or an array for a single regressor. In the `Dataset` form, `y` and `x` are the names of the series (`x` an array of names), observations with a missing value in any series used are left out, and the coefficients carry the names of the series.
- `constant`: include a constant, named `const`, as the first coefficient. The default is `true`.
- `vcov`: how the covariance matrix of the coefficients, and so the standard errors, is computed: `Vcov.classical` (the default, `sigma^2 (X'X)^-1`), `Vcov.hc0` (White's heteroskedasticity-consistent estimator), `Vcov.hc1` (HC0 times `n / (n - k)`: the usual "robust" standard errors), `Vcov.hc2` (squared residuals divided by `1 - h`, with `h` the leverage), `Vcov.hc3` (divided by `(1 - h)^2`, the jackknife-like version) or `Vcov.hac` (heteroskedasticity and autocorrelation consistent, Newey-West and relatives, with the data taken to be in time order).
- `lag`: for `Vcov.hac`, the number of lags; the default 0 lets the library choose from the sample size.
- `kernel`: for `Vcov.hac`, `HACKernel.bartlett` (the default), `HACKernel.parzen` or `HACKernel.qs`.
- `prewhiten`: for `Vcov.hac`, pre-whiten the data. The default is `false`.
- `cluster`: a series (an array with one value per observation; in the `Dataset` form, the name of a series) that identifies clusters. The standard errors are then made robust to correlation within clusters (cluster-robust), whatever `vcov` is. The default (null) is no clustering.

## Return value

An `Estimate` (see its page). The library settings that the options change (the HC version, the HAC kernel and lag, pre-whitening) are restored afterwards. A `GretlException` is thrown for missing or non-finite values, a wrong number of rows, too few observations, or a failure of the estimation. The `Estimate` holds the coefficients `coef` (named `const`, then `x1`, `x2`, ... or the series names), `se`, `tstat`, `pvalue` and `vcov`, the `resid` and `fitted` values, `rsq`, `adjrsq`, `fstat`, `sigma`, `ess`, `tss`, `lnL` and the information criteria (the Durbin-Watson statistic `dw` is NaN here, because data given as a matrix are not declared to have a time order); `method` is "OLS" and `vcovType` says which standard errors were used. A regressor that is perfectly collinear with the others is dropped by the library: the result then has fewer coefficients than columns, and `names` shows which remain.

The classical p-values use the t distribution with `df = nobs - ncoef` degrees of freedom.

## Example

```d
import std.stdio;
import std.math: fabs;
import dgretl;

void main() {
    setSeed(71);
    auto n = 200;
    auto x1 = rnorm(n);
    auto x2 = rnorm(n);
    auto e = rnorm(n);
    auto y = new double[n];
    foreach (t; 0 .. n) {
        y[t] = 1.0 + 2.0 * x1[t] - 1.0 * x2[t] + (1 + 0.8 * fabs(x1[t])) * e[t];     // the variance grows with |x1|
    }
    auto x = cbind(x1, x2);

    auto fit = ols(y, x);
    writeln(fit);                                     // a table of the coefficients
    writefln("coefficients: %(%.3f %)", Col(fit.coef, 0).dup);

    auto robust = ols(y, x, true, Vcov.hc1);          // heteroskedasticity-robust standard errors
    writefln("classical se of x1 %.4f, robust %.4f (%s)", fit.se[1, 0], robust.se[1, 0], robust.vcovType);

    auto noConstant = ols(y, x, false);               // without a constant
    writeln(noConstant.names);

    auto hac = ols(y, x, true, Vcov.hac, 4);          // HAC standard errors with 4 lags
    writefln("HAC se of x1 %.4f (%s)", hac.se[1, 0], hac.vcovType);
}
```

Output:

```
OLS, n = 200, df = 197, standard errors: classical
                     coef      std.err    t-ratio   p-value
const             1.05554     0.116413     9.0672    0.0000
x1                 1.9559     0.116069    16.8512    0.0000
x2               -1.03208     0.117985    -8.7475    0.0000
R-squared 0.685845, adjusted 0.682656, sigma 1.63992, F 215
log-likelihood -381.205891, AIC 768.4118, BIC 778.3067

coefficients: 1.056 1.956 -1.032
classical se of x1 0.1161, robust 0.1448 (HC1)
["x1", "x2"]
HAC se of x1 0.1431 (HAC)
```
