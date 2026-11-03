# ar1

## Signature

```d
Estimate ar1(Y, X)(Y y, X x, bool constant = true, Ar1Method method = Ar1Method.cochraneOrcutt);
Estimate ar1(Dataset d, string y, string[] x, bool constant = true, Ar1Method method = Ar1Method.cochraneOrcutt);

enum Ar1Method { cochraneOrcutt, praisWinsten, hildrethLu }
```

## Parameters

- `y`: the dependent series (an array or a one-column matrix).
- `x`: the regressors, as the columns of a matrix (or a single array), with as many rows as `y`. In the `Dataset` form, `y` and `x` are the names of the series. More than `ncol(x) + 5` observations are needed and there must be no missing values.
- `constant`: include a constant. The default is `true`.
- `method`: the estimator: `Ar1Method.cochraneOrcutt` (iterated Cochrane-Orcutt, the default), `praisWinsten` (iterated Prais-Winsten, which keeps the first observation) or `hildrethLu` (the Hildreth-Lu grid search).

## Return value

An `Estimate` for the regression `y = x b + e` with first-order autoregressive errors `e(t) = rho e(t - 1) + u(t)`, estimated on the transformed data. `coef`, `se`, `tstat`, `pvalue`, `names` ("const", `x1`, `x2`, ... or the series names), `resid`, `rsq` and the other fields are those of that regression. The error coefficient `rho` is in `arcoef` (a 1 x 1 matrix) with `arlags = [1]`; the library gives no standard error for it, so `arse` is NaN. `method` holds the name of the estimator. `writeln(estimate)` prints a table of the coefficients.

## Example

```d
import std.stdio;
import dgretl;

double[] arSeries(size_t n, double phi) {     // an AR(1) series from the library's random numbers
    auto e = rnorm(n);
    auto y = new double[n];

    y[0] = e[0];
    foreach (t; 1 .. n) {
        y[t] = phi * y[t - 1] + e[t];
    }
    return y;
}

void main() {
    setSeed(23);
    auto x = rnorm(200);
    auto u = arSeries(200, 0.6);                     // AR(1) errors with rho = 0.6
    auto y = new double[200];
    foreach (t; 0 .. 200) {
        y[t] = 1.0 + 2.0 * x[t] + u[t];
    }

    auto e = ar1(y, x);
    writeln(e);
    writefln("rho = %.3f", e.arcoef[0, 0]);

    auto p = ar1(y, x, true, Ar1Method.praisWinsten);
    writefln("Prais-Winsten: slope %.3f, rho %.3f", p.coef[1, 0], p.arcoef[0, 0]);
}
```

Output:

```
Cochrane-Orcutt, n = 199, df = 197, standard errors: classical
                     coef      std.err    t-ratio   p-value
const            0.975216     0.174362     5.5930    0.0000
x1                1.98081    0.0608607    32.5467    0.0000
R-squared 0.796945, adjusted 0.795914, sigma 1.02691, F 1059

rho = 0.582
Prais-Winsten: slope 1.981, rho 0.582
```
