# tsls

## Signature

```d
Estimate tsls(Y, E, Z)(Y y, E endog, Z instruments, GretlMatrix exog = GretlMatrix.init, bool constant = true, Vcov vcov = Vcov.classical,
                       int lag = 0, HACKernel kernel = HACKernel.bartlett, bool prewhiten = false);
Estimate tsls(Dataset d, string y, string[] endog, string[] instruments, string[] exog = null, bool constant = true, Vcov vcov = Vcov.classical,
              int lag = 0, HACKernel kernel = HACKernel.bartlett, bool prewhiten = false);
```

## Parameters

- `y`: the dependent variable, an array of numbers or a `GretlMatrix` or `Submatrix` with one column (n values).
- `endog`: the endogenous regressors, a matrix with n rows and `k1` columns (or an array for one). In the `Dataset` form, an array of series names.
- `instruments`: the EXCLUDED instruments, a matrix with n rows and `m` columns (or an array for one): the variables that are correlated with the endogenous regressors but not with the error. The included exogenous regressors and the constant are used as instruments as well. The model is identified if `m >= k1`; fewer instruments throw a `GretlException`. In the `Dataset` form, an array of series names.
- `exog`: the included exogenous regressors, a `GretlMatrix` with n rows. The default is none. In the `Dataset` form, an array of series names.
- `constant`: include a constant, named `const`, as the first coefficient. The default is `true`.
- `vcov`: how the covariance matrix of the coefficients, and so the standard errors, is computed: `Vcov.classical` (the default, `sigma^2 (X'X)^-1`), `Vcov.hc0` (White's heteroskedasticity-consistent estimator), `Vcov.hc1` (HC0 times `n / (n - k)`: the usual "robust" standard errors), `Vcov.hc2` (squared residuals divided by `1 - h`, with `h` the leverage), `Vcov.hc3` (divided by `(1 - h)^2`, the jackknife-like version) or `Vcov.hac` (heteroskedasticity and autocorrelation consistent, Newey-West and relatives, with the data taken to be in time order).
- `lag`, `kernel`, `prewhiten`: for `Vcov.hac`, as for `ols`. (A cluster-robust form is not available for two-stage least squares.)

## Return value

An `Estimate` (see its page). The library settings that the options change (the HC version, the HAC kernel and lag, pre-whitening) are restored afterwards. A `GretlException` is thrown for missing or non-finite values, a wrong number of rows, too few observations, or a failure of the estimation. The `Estimate` is as for `ols`, with `method` "2SLS". The coefficients come in the order: the constant, the exogenous regressors (named `w1`, `w2`, ...), then the endogenous regressors (`x1`, `x2`, ...); in the `Dataset` form they carry the names of the series. The library also computes the Sargan overidentification test for the model it builds, which is not copied into the `Estimate`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(73);
    auto n = 500;
    auto z = rnorm(n);                                // the instrument
    auto u = rnorm(n);                                // the unobserved factor that makes x endogenous
    auto v = rnorm(n);
    auto x = new double[n];
    auto y = new double[n];
    foreach (t; 0 .. n) {
        x[t] = 0.8 * z[t] + u[t] + 0.5 * v[t];
        y[t] = 1.0 + 2.0 * x[t] + 1.5 * u[t];         // the true slope is 2, but x and the error are correlated
    }

    auto biased = ols(y, x);
    auto iv = tsls(y, x, z);
    writefln("OLS slope %.3f (biased)", biased.coef[1, 0]);
    writefln("2SLS slope %.3f (se %.4f), method %s", iv.coef[1, 0], iv.se[1, 0], iv.method);
    writeln(iv.names);

    auto w = rnorm(n);                                // an included exogenous regressor
    auto withExog = tsls(y, x, z, cbind(w));
    writeln(withExog.names);
}
```

Output:

```
OLS slope 2.846 (biased)
2SLS slope 1.960 (se 0.0960), method 2SLS
["const", "x1"]
["const", "w1", "x1"]
```
