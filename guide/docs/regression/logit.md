# logit

## Signature

```d
Estimate logit(Y, X)(Y y, X x, bool constant = true, bool robust = false);
Estimate logit(Dataset d, string y, string[] x, bool constant = true, bool robust = false);
```

## Parameters

- `y`: the dependent variable, an array or a one-column matrix holding only 0s and 1s. Any other value throws a `GretlException`.
- `x`: the regressors: a `GretlMatrix` or `Submatrix` with n rows and one column for each regressor, or an array for a single regressor. In the `Dataset` form, `y` and `x` are the names of the series (`x` an array of names), observations with a missing value in any series used are left out, and the coefficients carry the names of the series.
- `constant`: include a constant, named `const`, as the first coefficient. The default is `true`.
- `robust`: use the quasi-maximum likelihood (sandwich) covariance matrix instead of the one from the Hessian. The default is `false`.

## Return value

An `Estimate` (see its page). The library settings that the options change (the HC version, the HAC kernel and lag, pre-whitening) are restored afterwards. A `GretlException` is thrown for missing or non-finite values, a wrong number of rows, too few observations, or a failure of the estimation. The `Estimate` is for the model `P(y = 1) = F(x b)` with `F` the logistic distribution function, estimated by maximum likelihood. `method` is "Logit" and `vcovType` is "ML (Hessian)" or "QML (sandwich)". `tstat` holds the z-ratios and `pvalue` the normal p-values (`df` is 0). `lnL` is the log-likelihood and `aic`, `bic` and `hqc` the information criteria; `rsq` is McFadden's R squared (`adjrsq`, `fstat` and `sigma` are whatever the library records for the model and should not be relied on). `fitted` is the estimated probability that `y = 1` and `resid` the generalized residuals of the library. `probit` is the same model with the normal distribution function.

## Example

```d
import std.stdio;
import std.math: exp;
import dgretl;

void main() {
    setSeed(74);
    auto n = 400;
    auto x = rnorm(n);
    auto u = runif(n);
    auto y = new double[n];
    foreach (t; 0 .. n) {
        const p = 1.0 / (1.0 + exp(-(0.5 + 1.5 * x[t])));     // the true probability
        y[t] = (u[t] < p) ? 1 : 0;
    }

    auto fit = logit(y, x);
    writeln(fit);
    writefln("McFadden R squared %.3f", fit.rsq);
    writefln("probability of y = 1 for the first observation: %.3f", fit.fitted[0, 0]);

    auto robust = logit(y, x, true, true);
    writefln("sandwich se of the slope %.4f (%s)", robust.se[1, 0], robust.vcovType);
}
```

Output:

```
Logit, n = 400, df = 0, standard errors: ML (Hessian)
                     coef      std.err    t-ratio   p-value
const             0.56924     0.124968     4.5551    0.0000
x1                1.48416     0.163014     9.1045    0.0000
R-squared 0.243192, adjusted 0.235805, sigma 0.417088, F 158.2
log-likelihood -204.900804, AIC 413.8016, BIC 421.7845

McFadden R squared 0.243
probability of y = 1 for the first observation: 0.711
sandwich se of the slope 0.1642 (QML (sandwich))
```
