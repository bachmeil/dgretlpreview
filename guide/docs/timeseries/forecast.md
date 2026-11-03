# forecast

## Signature

```d
VarForecast forecast(F)(ref F fit, int steps, GretlMatrix futureExog = GretlMatrix.init);   // a VarFit or VecmFit
SeriesForecast forecast(ref ArmaFit fit, int steps, GretlMatrix futureExog = GretlMatrix.init);
SeriesForecast forecast(ref GarchFit fit, int steps, GretlMatrix futureExog = GretlMatrix.init);

struct VarForecast {
    int first;
    GretlMatrix mean, se;      // rows x neqns
}
struct SeriesForecast {
    int first;
    GretlMatrix mean, se;      // rows x 1
}
```

## Parameters

- `fit`: a `VarFit` or `VecmFit` (from `var`, `vecm`), an `ArmaFit` (from `arma`) or a `GarchFit` (from `garch`).
- `steps`: the number of observations to forecast after the end of the estimation sample, at least 1.
- `futureExog`: the future values of the exogenous regressors, a `steps` x m matrix in the order of the model's `exog`, when the model has any. Seasonal dummies and the trend continue by themselves. The default is none.

## Return value

The forecasts for the `steps` observations after the estimation sample, made dynamically (each forecast uses the earlier ones).

- For a VAR or VECM, a `VarForecast`: `mean` and `se` are `steps` x `neqns` matrices, the point forecast and the standard error of the forecast error for each variable at each step.
- For an ARMA or GARCH model, a `SeriesForecast` with `steps` x 1 matrices. For an ARIMA model the forecasts are of the level of the series. For a GARCH model `mean` is the forecast of the mean equation and `se` the forecast of the conditional standard deviation `sqrt(h)`.

`first` is the 0-based observation number of row 0. The standard errors come from the library and ignore the uncertainty of the coefficients. Too few future regressor rows throw a `GretlException`. `forecastRange` forecasts any range of observations, including ones inside the sample.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(31);
    auto n = 200;
    auto e1 = rnorm(n);
    auto e2 = rnorm(n);
    auto y1 = new double[n];
    auto y2 = new double[n];
    y1[0] = e1[0];
    y2[0] = e2[0];
    foreach (t; 1 .. n) {
        y1[t] = 0.5 * y1[t - 1] + 0.2 * y2[t - 1] + e1[t];
        y2[t] = 0.3 * y1[t - 1] + 0.4 * y2[t - 1] + e2[t];
    }
    auto fit = var(cbind(y1, y2), 1);

    auto f = forecast(fit, 4);
    writefln("the forecasts start at observation %d", f.first);
    foreach (h; 0 .. 4) {
        writefln("step %d: y1 %7.4f (se %.3f)  y2 %7.4f (se %.3f)", h + 1, f.mean[h, 0], f.se[h, 0], f.mean[h, 1], f.se[h, 1]);
    }

    auto ar = arma(y1, 1, 0);                        // the same function for an ARMA model
    auto g = forecast(ar, 3);
    writefln("ARMA forecasts: %(%.4f %)", Col(g.mean, 0).dup);
}
```

Output:

```
the forecasts start at observation 200
step 1: y1 -0.9275 (se 1.009)  y2 -0.6151 (se 1.052)
step 2: y1 -0.7530 (se 1.162)  y2 -0.4386 (se 1.158)
step 3: y1 -0.6247 (se 1.226)  y2 -0.3223 (se 1.203)
step 4: y1 -0.5330 (se 1.256)  y2 -0.2420 (se 1.226)
ARMA forecasts: -0.7891 -0.5790 -0.4542
```
