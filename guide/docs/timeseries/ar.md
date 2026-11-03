# ar

## Signature

```d
Estimate ar(Y, X)(Y y, X x, int[] lags, bool constant = true);
Estimate ar(Y, X)(Y y, X x, int p, bool constant = true);
Estimate ar(Dataset d, string y, string[] x, int[] lags, bool constant = true);
Estimate ar(Dataset d, string y, string[] x, int p, bool constant = true);
```

## Parameters

- `y`: the dependent series (an array or a one-column matrix).
- `x`: the regressors, as the columns of a matrix (or a single array), with as many rows as `y`. In the `Dataset` forms, `y` and `x` are the names of the series. There must be no missing values.
- `lags`: the lags of the autoregressive error, as an array such as `[1, 2]` or `[1, 4]` (at least one, each at least 1); or `p`, a single integer for lags 1 to `p`.
- `constant`: include a constant. The default is `true`.

## Return value

An `Estimate` for the regression `y = x b + e` with autoregressive errors `e(t) = rho_1 e(t - lag_1) + rho_2 e(t - lag_2) + ... + u(t)` at the given lags, estimated by the generalized Cochrane-Orcutt iteration. The regression coefficients are in `coef`, `se` and the other fields, with `names`; the error coefficients are in `arcoef` (one row for each lag), their standard errors in `arse`, and the lags they belong to in `arlags`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    setSeed(24);
    auto x = rnorm(300);
    auto e = rnorm(300);
    auto u = new double[300];                        // errors with AR(2) structure
    u[0] = e[0];
    u[1] = e[1];
    foreach (t; 2 .. 300) {
        u[t] = 0.5 * u[t - 1] + 0.3 * u[t - 2] + e[t];
    }
    auto y = new double[300];
    foreach (t; 0 .. 300) {
        y[t] = 1.0 + 2.0 * x[t] + u[t];
    }

    auto m = ar(y, x, 2);                            // errors with lags 1 and 2
    writefln("slope %.3f", m.coef[1, 0]);
    writefln("error coefficients %.3f and %.3f at lags %s", m.arcoef[0, 0], m.arcoef[1, 0], m.arlags);

    auto s = ar(y, x, [1, 4]);                       // lags 1 and 4 only
    writeln(s.arlags);
}
```

Output:

```
slope 1.942
error coefficients 0.415 and 0.331 at lags [1, 2]
[1, 4]
```
