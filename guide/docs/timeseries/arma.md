# arma

## Signature

```d
ArmaFit arma(Y)(Y y, int p, int q, int d = 0, int[] seasonal = null, GretlMatrix exog = GretlMatrix.init, bool constant = true,
                ArmaMethod method = ArmaMethod.exact, bool opg = false, string[] exogNames = null, int frequency = 1, string start = "1",
                string fitStart = null, string fitEnd = null);
ArmaFit arma(Dataset data, string y, int p, int q, int d = 0, int[] seasonal = null, string[] exog = null, bool constant = true,
             ArmaMethod method = ArmaMethod.exact, bool opg = false, int frequency = 0, string start = null, string fitStart = null, string fitEnd = null);

enum ArmaMethod { exact, kalman, conditional }

struct ArmaFit {
    Estimate fit;   // alias this: coef, se, names, resid, ...
    int p, d, q, sP, sD, sQ;
    GretlMatrix arRoots, maRoots;
}
```

## Parameters

- `y`: the series (an array or a one-column matrix), or in the `Dataset` form the name of a series. No missing values.
- `p`, `q`: the orders of the autoregressive and moving-average parts. At least one of `p`, `q` and the seasonal orders must be positive.
- `d`: the order of differencing for an ARIMA model. The default is 0.
- `seasonal`: `[P, D, Q]`, the seasonal AR order, the number of seasonal differences and the seasonal MA order (needs a frequency above 1). The default (null) is none.
- `exog`: regressors for an ARMAX model, a matrix with the same rows as `y` (in the `Dataset` form, an array of series names). The default is none.
- `constant`: include a constant. The default is `true`; for an integrated model there is none to estimate and the regressors are differenced too.
- `method`: `ArmaMethod.exact` (exact Gaussian maximum likelihood by the library's default routines), `kalman` (exact maximum likelihood by the Kalman filter) or `conditional` (conditional maximum likelihood, with the pre-sample innovations set to zero).
- `opg`: estimate the covariance matrix from the outer product of the gradient instead of the Hessian. The default is `false`.
- `exogNames`, `frequency`, `start`, `fitStart`, `fitEnd`: as for `var`. In the `Dataset` form the calendar of the `Dataset` is used.

## Return value

An `ArmaFit` struct for the model `y(t) = const + phi_1 y(t - 1) + ... + phi_p y(t - p) + e(t) + theta_1 e(t - 1) + ... + theta_q e(t - q)` (for `d > 0` or a seasonal difference, the ARIMA model for the differenced series). It has everything of an `Estimate` (through `alias this`): `coef`, `se`, `tstat`, `pvalue`, `vcov`, `names` (`const`, `phi_1`, ..., `theta_1`, ..., `Phi_1`, ..., `Theta_1`, ..., then the regressors), `resid` and `fitted` (the one-step predictions of the differenced model), `lnL`, `aic`, `bic` and `hqc`. In addition:

- `p`, `d`, `q`, `sP`, `sD`, `sQ`: the orders.
- `arRoots`, `maRoots`: the roots of the AR and MA polynomials (non-seasonal, then seasonal), one row each, with the real part, the imaginary part and the modulus. Moduli above 1 mean stationarity (AR) and invertibility (MA).

The fit keeps the library's model, so `forecast` and `forecastRange` can be applied to it. A GretlException is thrown for missing values, too few observations, or a failure to estimate.

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
    setSeed(51);
    auto n = 400;
    auto e = rnorm(n + 1);
    auto y = new double[n];                          // y(t) = 0.6 y(t-1) + e(t) + 0.4 e(t-1)
    y[0] = e[1];
    foreach (t; 1 .. n) {
        y[t] = 0.6 * y[t - 1] + e[t + 1] + 0.4 * e[t];
    }

    auto fit = arma(y, 1, 1);
    writeln(fit.names);
    writefln("phi %.3f (se %.3f), theta %.3f (se %.3f)", fit.coef[1, 0], fit.se[1, 0], fit.coef[2, 0], fit.se[2, 0]);
    writefln("AR root modulus %.3f, MA root modulus %.3f", fit.arRoots[0, 2], fit.maRoots[0, 2]);

    auto cond = arma(y, 1, 1, 0, null, GretlMatrix.init, true, ArmaMethod.conditional);
    writefln("conditional ML: phi %.3f, theta %.3f", cond.coef[1, 0], cond.coef[2, 0]);

    auto integrated = arma(arSeries(n, 1.0), 1, 0, 1);         // an ARIMA(1,1,0) for a random walk
    writefln("ARIMA(1,1,0): %s", integrated.method);
}
```

Output:

```
["const", "phi_1", "theta_1"]
phi 0.555 (se 0.051), theta 0.450 (se 0.054)
AR root modulus 1.803, MA root modulus 2.221
conditional ML: phi 0.559, theta 0.446
ARIMA(1,1,0): ARIMA(1,1,0)
```
