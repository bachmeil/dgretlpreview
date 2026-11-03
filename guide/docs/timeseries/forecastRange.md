# forecastRange

## Signature

```d
VarForecast forecastRange(F)(ref F fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init);   // a VarFit or VecmFit
SeriesForecast forecastRange(ref ArmaFit fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init);
SeriesForecast forecastRange(ref GarchFit fit, int first, int last, ForecastMethod method = ForecastMethod.automatic, GretlMatrix futureExog = GretlMatrix.init);

enum ForecastMethod { automatic, oneStep, dynamic }
```

## Parameters

- `fit`: a `VarFit` or `VecmFit`, an `ArmaFit` or a `GarchFit`.
- `first`, `last`: the range of observations to forecast, as 0-based rows counting from the beginning of the data the model was estimated on, with `0 <= first <= last`. The range may lie inside the estimation sample, after it, or across its end.
- `method`: `ForecastMethod.automatic` (the default) makes one-step-ahead forecasts from the actual data inside the estimation sample and dynamic forecasts after it; `oneStep` and `dynamic` force one kind over the whole range (a dynamic forecast started inside the sample uses the model's own forecasts from `first` on). A dynamic forecast is not applicable to a GARCH model and throws a `GretlException`.
- `futureExog`: the regressors for observations after the end of the data, if the model has any: row 0 is the first observation after the data.

## Return value

A `VarForecast` (for a VAR or VECM) or a `SeriesForecast` (for ARMA and GARCH): `first` is the observation of row 0 and `mean` and `se` hold the forecasts and the standard errors of the forecast errors for observations `first` to `last` (so `last - first + 1` rows; see `forecast` for the shapes). The standard errors of one-step forecasts of a VAR are those of the residuals; one-step forecasts of an ARMA model have no standard error (NaN); for a GARCH model it is the conditional standard deviation. The first `order` observations of a VAR have no forecast (NaN).

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

    auto inside = forecastRange(fit, 100, 102);      // one-step forecasts inside the sample
    foreach (h; 0 .. 3) {
        writefln("observation %d: forecast %7.4f, actual %7.4f", inside.first + h, inside.mean[h, 0], y1[100 + h]);
    }

    auto across = forecastRange(fit, 198, 201);      // two inside the sample, two after its end
    writefln("%d rows from observation %d", across.mean.nrow, across.first);

    auto dyn = forecastRange(fit, 100, 103, ForecastMethod.dynamic);   // dynamic from observation 100
    writefln("dynamic forecast for y1: %(%.4f %)", Col(dyn.mean, 0).dup);
}
```

Output:

```
observation 100: forecast  0.7244, actual  0.3617
observation 101: forecast  0.0785, actual -0.0288
observation 102: forecast -0.1561, actual -0.1338
4 rows from observation 198
dynamic forecast for y1: 0.7244 0.3872 0.1747 0.0305
```
