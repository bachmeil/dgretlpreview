# probit

## Signature

```d
Estimate probit(Y, X)(Y y, X x, bool constant = true, bool robust = false);
Estimate probit(Dataset d, string y, string[] x, bool constant = true, bool robust = false);
```

## Parameters

- `y`: the dependent variable, an array or a one-column matrix holding only 0s and 1s. Any other value throws a `GretlException`.
- `x`: the regressors: a `GretlMatrix` or `Submatrix` with n rows and one column for each regressor, or an array for a single regressor. In the `Dataset` form, `y` and `x` are the names of the series (`x` an array of names), observations with a missing value in any series used are left out, and the coefficients carry the names of the series.
- `constant`: include a constant, named `const`, as the first coefficient. The default is `true`.
- `robust`: use the quasi-maximum likelihood (sandwich) covariance matrix instead of the one from the Hessian. The default is `false`.

## Return value

An `Estimate` (see its page). The library settings that the options change (the HC version, the HAC kernel and lag, pre-whitening) are restored afterwards. A `GretlException` is thrown for missing or non-finite values, a wrong number of rows, too few observations, or a failure of the estimation. The `Estimate` is for the model `P(y = 1) = F(x b)` with `F` the standard normal distribution function, estimated by maximum likelihood. `method` is "Probit" and `vcovType` is "ML (Hessian)" or "QML (sandwich)". `tstat` holds the z-ratios and `pvalue` the normal p-values (`df` is 0). `lnL` is the log-likelihood and `aic`, `bic` and `hqc` the information criteria; `rsq` is McFadden's R squared (`adjrsq`, `fstat` and `sigma` are whatever the library records for the model and should not be relied on). `fitted` is the estimated probability that `y = 1` and `resid` the generalized residuals of the library. `logit` is the same model with the logistic distribution function.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(75);
    auto n = 400;
    auto x = rnorm(n);
    auto e = rnorm(n);
    auto y = new double[n];
    foreach (t; 0 .. n) {
        y[t] = (0.3 + 1.0 * x[t] + e[t] > 0) ? 1 : 0;       // a latent-variable probit with slope 1
    }

    auto fit = probit(y, x);
    writeln(fit);
    writefln("log-likelihood %.3f, AIC %.3f", fit.lnL, fit.aic);

    auto l = logit(y, x);                                   // the logit slope is about 1.7 times larger
    writefln("probit slope %.3f, logit slope %.3f", fit.coef[1, 0], l.coef[1, 0]);
}
```

Output:

```
Probit, n = 400, df = 0, standard errors: ML (Hessian)
                     coef      std.err    t-ratio   p-value
const            0.223421    0.0753418     2.9654    0.0030
x1                1.21431     0.113899    10.6613    0.0000
R-squared 0.322710, adjusted 0.315427, sigma 0.401141, F 215.2
log-likelihood -185.989209, AIC 375.9784, BIC 383.9613

log-likelihood -185.989, AIC 375.978
probit slope 1.214, logit slope 2.051
```
