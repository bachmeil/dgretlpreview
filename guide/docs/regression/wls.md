# wls

## Signature

```d
Estimate wls(Y, X)(Y y, X x, const(double)[] w, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0, HACKernel kernel = HACKernel.bartlett,
                    bool prewhiten = false, const(double)[] cluster = null);
Estimate wls(Dataset d, string y, string[] x, string weights, bool constant = true, Vcov vcov = Vcov.classical, int lag = 0,
             HACKernel kernel = HACKernel.bartlett, bool prewhiten = false, string cluster = null);
```

## Parameters

- `y`: the dependent variable, an array of numbers or a `GretlMatrix` or `Submatrix` with one column (n values).
- `x`: the regressors: a `GretlMatrix` or `Submatrix` with n rows and one column for each regressor, or an array for a single regressor. In the `Dataset` form, `y` and `x` are the names of the series (`x` an array of names), observations with a missing value in any series used are left out, and the coefficients carry the names of the series.
- `w`: the weights, one for each observation (in the `Dataset` form, `weights` is the name of a series). A weight is the inverse of the error variance of the observation: the sum of `w * residual^2` is minimized, so an observation with weight 4 counts as four observations of weight 1. Weights must be positive and finite, except that a series of only 0s and 1s is also accepted and leaves out the observations with weight 0. Any other zero weight, or a negative weight, throws a `GretlException`.
- `constant`: include a constant, named `const`, as the first coefficient. The default is `true`.
- `vcov`: how the covariance matrix of the coefficients, and so the standard errors, is computed: `Vcov.classical` (the default, `sigma^2 (X'X)^-1`), `Vcov.hc0` (White's heteroskedasticity-consistent estimator), `Vcov.hc1` (HC0 times `n / (n - k)`: the usual "robust" standard errors), `Vcov.hc2` (squared residuals divided by `1 - h`, with `h` the leverage), `Vcov.hc3` (divided by `(1 - h)^2`, the jackknife-like version) or `Vcov.hac` (heteroskedasticity and autocorrelation consistent, Newey-West and relatives, with the data taken to be in time order).
- `lag`, `kernel`, `prewhiten`: for `Vcov.hac`, as for `ols`.
- `cluster`: a series (an array with one value per observation; in the `Dataset` form, the name of a series) that identifies clusters. The standard errors are then made robust to correlation within clusters (cluster-robust), whatever `vcov` is. The default (null) is no clustering.

## Return value

An `Estimate` (see its page). The library settings that the options change (the HC version, the HAC kernel and lag, pre-whitening) are restored afterwards. A `GretlException` is thrown for missing or non-finite values, a wrong number of rows, too few observations, or a failure of the estimation. The `Estimate` is as for `ols`, with `method` "WLS". The `resid` and `fitted` values are the unweighted ones, `y - x * coef`; `sigma`, `ess`, `rsq` and the other summary statistics are those of the weighted regression.

## Example

```d
import std.stdio;
import std.math: fabs;
import dgretl;

void main() {
    setSeed(72);
    auto n = 200;
    auto x = rnorm(n);
    auto sd = new double[n];                          // the error standard deviation grows with x
    auto y = new double[n];
    auto w = new double[n];
    auto e = rnorm(n);
    foreach (t; 0 .. n) {
        sd[t] = 0.5 + fabs(x[t]);
        y[t] = 1.0 + 2.0 * x[t] + sd[t] * e[t];
        w[t] = 1.0 / (sd[t] * sd[t]);                 // the weight is the inverse of the variance
    }

    auto plain = ols(y, x);
    auto weighted = wls(y, x, w);
    writefln("OLS slope %.3f (se %.4f)", plain.coef[1, 0], plain.se[1, 0]);
    writefln("WLS slope %.3f (se %.4f), method %s", weighted.coef[1, 0], weighted.se[1, 0], weighted.method);

    auto sample = new double[n];                      // a series of 0s and 1s leaves out observations
    foreach (t; 0 .. n) {
        sample[t] = (t < 150) ? 1 : 0;
    }
    auto subset = wls(y, x, sample);
    writefln("using the first %d observations: slope %.3f", subset.nobs, subset.coef[1, 0]);
}
```

Output:

```
OLS slope 2.104 (se 0.1086)
WLS slope 2.083 (se 0.1336), method WLS
using the first 150 observations: slope 2.002
```
