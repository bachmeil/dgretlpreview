# arch

## Signature

```d
Estimate arch(Y, X)(Y y, X x, int order, bool constant = true);
Estimate arch(Dataset d, string y, string[] x, int order, bool constant = true);
```

## Parameters

- `y`: the dependent series (an array or a one-column matrix).
- `x`: the regressors, as the columns of a matrix (or a single array), with as many rows as `y`. In the `Dataset` form, `y` and `x` are the names of the series. There must be no missing values.
- `order`: the number of ARCH lags, at least 1, with more than `order + ncol(x) + 3` observations.
- `constant`: include a constant. The default is `true`.

## Return value

An `Estimate` from the two-step weighted least squares estimate of a regression with ARCH(`order`) errors: OLS of `y` on `x`, a regression of the squared residuals on their own `order` lags, then WLS with the fitted variances as weights. The result (`method` "ARCH-WLS") is the final WLS regression, with `coef`, `se`, `names` and the other fields. For maximum likelihood estimation of the variance equation see `garch`.

## Example

```d
import std.stdio;
import std.math: sqrt;
import dgretl;

void main() {
    setSeed(25);
    auto x = rnorm(500);
    auto z = rnorm(500);
    auto u = new double[500];                        // ARCH(1) errors
    u[0] = z[0];
    foreach (t; 1 .. 500) {
        u[t] = sqrt(0.3 + 0.6 * u[t - 1] ^^ 2) * z[t];
    }
    auto y = new double[500];
    foreach (t; 0 .. 500) {
        y[t] = 1.0 + 2.0 * x[t] + u[t];
    }

    auto e = arch(y, x, 1);
    writeln(e);
}
```

Output:

```
ARCH-WLS, n = 499, df = 497, standard errors: classical
                     coef      std.err    t-ratio   p-value
const             1.04276    0.0315904    33.0086    0.0000
x1                2.03141    0.0326289    62.2578    0.0000
R-squared 0.886349, adjusted 0.886120, sigma 0.957827, F 3876
log-likelihood -685.547395, AIC 1375.0948, BIC 1383.5200
```
