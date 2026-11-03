# garch

## Signature

```d
GarchFit garch(Y)(Y y, int p, int q, GretlMatrix exog = GretlMatrix.init, bool constant = true, bool robust = false, bool fcp = false,
                  bool standardized = false, bool armaInit = false, string[] exogNames = null);
GarchFit garch(Dataset data, string y, int p, int q, string[] exog = null, bool constant = true, bool robust = false, bool fcp = false,
               bool standardized = false, bool armaInit = false);

struct GarchFit {
    Estimate fit;   // alias this: coef, se, names, resid, ...
    int p, q;
    GretlMatrix variance;
    double lrStat;
    int lrDf;
}
```

## Parameters

- `y`: the series (an array or a one-column matrix), or in the `Dataset` form the name of a series. No missing values, and more than `p + q + ncol(exog) + 20` observations.
- `p`: the number of GARCH lags (of the conditional variance `h`), at least 0.
- `q`: the number of ARCH lags (of the squared shock), at least 1. As in the library, `garch(y, 1, 1)` is the usual GARCH(1,1); `p + q` must not exceed 7.
- `exog`: regressors of the mean equation, a matrix with the same rows as `y` (in the `Dataset` form, an array of series names). The default is none.
- `constant`: include a constant in the mean equation. The default is `true`.
- `robust`: use the quasi-maximum likelihood (robust) covariance matrix. The default is `false`.
- `fcp`: use the Fiorentini-Calzolari-Panattoni algorithm instead of the library's BFGS-based estimator. The default is `false`.
- `standardized`: return the residuals divided by the conditional standard deviation. The default is `false`.
- `armaInit`: start the variance parameters from an ARMA fit of the squared residuals. The default is `false`.
- `exogNames`: the names of the columns of `exog`.

## Return value

A `GarchFit` struct for the model `y(t) = const + x(t) b + e(t)` with `e(t)` normal given the past with variance `h(t) = alpha_0 + alpha_1 e(t-1)^2 + ... + alpha_q e(t-q)^2 + beta_1 h(t-1) + ... + beta_p h(t-p)`. It has everything of an `Estimate` (through `alias this`): the mean equation and the variance parameters together in `coef`, `se`, `tstat`, `pvalue` and `vcov`, with `names` `const`, the regressors, `alpha(0)`, `alpha(1)`, ..., `beta(1)`, ..., and `resid`, `lnL`, `aic`, `bic`. In addition:

- `p`, `q`: the orders.
- `variance`: an `nobs` x 1 matrix with the conditional variance `h(t)`, in the units of the data.
- `lrStat`, `lrDf`: the likelihood-ratio statistic of the GARCH terms against the plain regression and its degrees of freedom. Its p-value is not the usual chi-square one, because the null is on the boundary of the parameter space.

The fit keeps the library's model, so `forecast` and `forecastRange` can be applied to it.

## Example

```d
import std.stdio;
import std.math: sqrt;
import dgretl;

void main() {
    setSeed(61);
    auto n = 1000;
    auto z = rnorm(n);
    auto y = new double[n];                          // a GARCH(1,1) process
    auto h = new double[n];
    h[0] = 0.2 / (1 - 0.15 - 0.8);
    y[0] = sqrt(h[0]) * z[0];
    foreach (t; 1 .. n) {
        h[t] = 0.2 + 0.15 * y[t - 1] ^^ 2 + 0.8 * h[t - 1];
        y[t] = sqrt(h[t]) * z[t];
    }

    auto fit = garch(y, 1, 1);
    foreach (i, name; fit.names) {
        writefln("%-8s %8.4f (se %.4f)", name, fit.coef[i, 0], fit.se[i, 0]);
    }
    writefln("conditional variance: first %.3f, last %.3f", fit.variance[0, 0], fit.variance[fit.variance.nrow - 1, 0]);
    writefln("log-likelihood %.3f", fit.lnL);
}
```

Output:

```
const     -0.0754 (se 0.0578)
alpha(0)   0.3776 (se 0.1094)
alpha(1)   0.1763 (se 0.0333)
beta(1)    0.7456 (se 0.0442)
conditional variance: first 4.901, last 9.521
log-likelihood -2123.899
```
