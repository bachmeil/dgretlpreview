# Estimate

## Signature

```d
struct Estimate {
    string method;            // "OLS", "WLS", "2SLS", "Logit", "Probit", ...
    string vcovType;          // "classical", "HC1", "HAC", ...
    string[] names;
    GretlMatrix coef, se, tstat, pvalue;   // each k x 1
    GretlMatrix vcov;                      // k x k
    GretlMatrix resid, fitted;             // each n x 1
    size_t nobs, ncoef;
    int df;
    double sigma, ess, tss, rsq, adjrsq, fstat, lnL, aic, bic, hqc, dw, rho;
    GretlMatrix arcoef, arse;
    int[] arlags;
    string toString();
}
```

## Parameters

- `method`: the name of the estimator ("OLS", "WLS", "2SLS", "Logit", "Probit", "Cochrane-Orcutt", "AR", "ARCH-WLS", "GARCH", an ARIMA label, ...).
- `vcovType`: what the standard errors are ("classical", "HC1", "HAC", "cluster-robust", "ML (Hessian)", ...).
- `names`: the names of the coefficients, in order: `const`, then the names of the regressors (`x1`, `x2`, ... for matrix data, the series names for a `Dataset`).
- `coef`, `se`, `tstat`, `pvalue`: k x 1 matrices of the coefficients, their standard errors, the ratios `coef / se` (z-ratios for the maximum likelihood estimators) and the two-sided p-values (Student t with `df` degrees of freedom, or normal).
- `vcov`: the k x k covariance matrix of the coefficients.
- `resid`, `fitted`: n x 1 matrices of the residuals (`y - fitted` for the least squares estimators) and the fitted values.
- `nobs`: the number of observations used. `ncoef`: the number of coefficients (a perfectly collinear regressor is dropped by the library). `df`: the residual degrees of freedom, `nobs - ncoef` for least squares (0 for the maximum likelihood binary-choice estimators).
- `sigma`: the standard error of the regression; `ess`, `tss`: the error and the total sum of squares; `rsq`, `adjrsq`: R squared and adjusted R squared; `fstat`: the overall F statistic. Fields that do not apply to an estimator are NaN; for `logit` and `probit` only `rsq` (McFadden's) is meaningful among them.
- `lnL`: the log-likelihood; `aic`, `bic`, `hqc`: the information criteria (the library's definitions).
- `dw`: the Durbin-Watson statistic and `rho`: the first-order autocorrelation of the residuals; NaN where the library does not compute them, as for the regressions on matrices.
- `arcoef`, `arse`, `arlags`: for models with autoregressive errors (`ar1`, `ar`), the coefficients of the errors, their standard errors and their lags; empty otherwise.
- `toString()`: a table of the coefficients with the usual summary lines, so `writeln(estimate)` prints it.

## Return value

An `Estimate` is what the estimators return (`ols`, `wls`, `tsls`, `logit`, `probit`, and `ar1`, `ar`, `arch` of the time-series section; `ArmaFit` and `GarchFit` contain one and expose its fields through `alias this`). It is a plain struct: the matrices are `GretlMatrix` values and the rest are numbers and strings, so nothing has to be freed by hand. The coefficient statistics are k x 1 matrices in the order of `names`, and `resid` and `fitted` are n x 1.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(71);
    auto n = 100;
    auto x = rnorm(n);
    auto e = rnorm(n);
    auto y = new double[n];
    foreach (t; 0 .. n) {
        y[t] = 1.0 + 2.0 * x[t] + e[t];
    }
    auto fit = ols(y, x);

    writeln(fit.method, ", ", fit.vcovType, ", ", fit.names);
    writefln("%d observations, %d coefficients, %d degrees of freedom", fit.nobs, fit.ncoef, fit.df);
    foreach (i, name; fit.names) {
        writefln("%-6s coef %7.4f  se %.4f  t %8.3f  p %.4f", name, fit.coef[i, 0], fit.se[i, 0], fit.tstat[i, 0], fit.pvalue[i, 0]);
    }
    writefln("R squared %.4f, adjusted %.4f, F %.2f, sigma %.4f", fit.rsq, fit.adjrsq, fit.fstat, fit.sigma);
    writefln("log-likelihood %.3f, AIC %.3f, BIC %.3f, Durbin-Watson %s", fit.lnL, fit.aic, fit.bic, fit.dw);
    writefln("residuals: %d x %d, the first is %.4f; fitted %.4f", fit.resid.nrow, fit.resid.ncol, fit.resid[0, 0], fit.fitted[0, 0]);
    writefln("covariance matrix: %d x %d, variance of the slope %.5f", fit.vcov.nrow, fit.vcov.ncol, fit.vcov[1, 1]);
}
```

Output:

```
OLS, classical, ["const", "x1"]
100 observations, 2 coefficients, 98 degrees of freedom
const  coef  1.0664  se 0.1091  t    9.774  p 0.0000
x1     coef  1.9833  se 0.1141  t   17.382  p 0.0000
R squared 0.7551, adjusted 0.7526, F 302.15, sigma 1.0874
log-likelihood -149.267, AIC 302.533, BIC 307.744, Durbin-Watson nan
residuals: 100 x 1, the first is -0.2327; fitted 2.4039
covariance matrix: 2 x 2, variance of the slope 0.01302
```
