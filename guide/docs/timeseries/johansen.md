# johansen

## Signature

```d
JohansenTest johansen(Y)(Y y, int order, Deterministic deterministic = Deterministic.unrestrictedConstant, GretlMatrix exog = GretlMatrix.init,
                         bool seasonals = false, bool asymptotic = false, string[] names = null, string[] exogNames = null, int frequency = 1,
                         string start = "1", string fitStart = null, string fitEnd = null);
JohansenTest johansen(Dataset d, string[] endog, int order, Deterministic deterministic = Deterministic.unrestrictedConstant, string[] exog = null,
                      bool seasonals = false, bool asymptotic = false, int frequency = 0, string start = null, string fitStart = null, string fitEnd = null);

enum Deterministic { none, restrictedConstant, unrestrictedConstant, restrictedTrend, unrestrictedTrend }

struct JohansenTest {
    size_t neqns, nobs;
    int order;
    GretlMatrix eigenvalues, trace, lmax, tracePvalue, lmaxPvalue, beta, alpha;
    double lnL;
    string report;
}
```

## Parameters

- `y`: the data of the series as the columns of a matrix, at least two. No missing values.
- `order`: the lag order of the VAR in levels, at least 1 (so `order - 1` lagged differences in the error-correction form).
- `deterministic`: which of Johansen's five cases: `none` (no constant), `restrictedConstant` (a constant restricted to the cointegrating space), `unrestrictedConstant` (the default: drift in the levels), `restrictedTrend` (an unrestricted constant and a trend restricted to the cointegrating space) or `unrestrictedTrend` (an unrestricted constant and trend).
- `exog`: unrestricted exogenous regressors (a matrix with the same rows as `y`). The default is none.
- `seasonals`: add centered seasonal dummies (needs a frequency above 1).
- `asymptotic`: use asymptotic p-values for the trace test instead of the ones corrected for the sample size. The default is `false`.
- `names`, `exogNames`, `frequency`, `start`, `fitStart`, `fitEnd`: as for `var`.

## Return value

A `JohansenTest` struct with the trace and maximum-eigenvalue tests of the cointegration rank. Each of the `n` x 1 matrices has row `r` for the null hypothesis that the rank is at most `r`:

- `eigenvalues`: the eigenvalues of the Johansen problem, decreasing.
- `trace`, `lmax`: the trace and maximum-eigenvalue statistics, with `tracePvalue` and `lmaxPvalue` (the trace p-values are corrected for the sample size unless `asymptotic`; those for `lmax` are asymptotic).
- `beta`: the eigenvectors, normalized so that `beta' S11 beta = I` (restricted deterministic terms add rows) and `alpha`: the matching adjustment coefficients.
- `neqns`, `nobs`, `order`, `lnL` (the log-likelihood) and `report`, the text of the library's table.

The first rank `r` for which the null is not rejected is the estimate of the rank. `vecm` estimates the model for a chosen rank.

## Example

```d
import std.stdio;
import dgretl;

void show(M)(M m) {
    foreach (i; 0 .. m.nrow) {
        foreach (j; 0 .. m.ncol) {
            writef("%10.4f", m[i, j]);
        }
        writeln();
    }
}

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
    setSeed(41);
    auto n = 250;
    auto walk = arSeries(n, 1.0);                    // a common stochastic trend
    auto e = arSeries(n, 0.3);
    auto y1 = walk.dup;
    auto y2 = new double[n];
    foreach (t; 0 .. n) {
        y2[t] = walk[t] + e[t];                      // y1 and y2 are cointegrated
    }

    auto j = johansen(cbind(y1, y2), 2);
    foreach (r; 0 .. 2) {
        writefln("rank <= %d: eigenvalue %.4f, trace %8.3f (p = %.4f), lmax %8.3f (p = %.4f)", r, j.eigenvalues[r, 0], j.trace[r, 0], j.tracePvalue[r, 0], j.lmax[r, 0], j.lmaxPvalue[r, 0]);
    }
    show(j.beta);                                    // the first column is the cointegrating vector
}
```

Output:

```
rank <= 0: eigenvalue 0.3251, trace  101.148 (p = 0.0000), lmax   97.496 (p = 0.0000)
rank <= 1: eigenvalue 0.0146, trace    3.652 (p = 0.0571), lmax    3.652 (p = 0.0560)
   -1.3260   -0.2666
    1.3269    0.0047
```
