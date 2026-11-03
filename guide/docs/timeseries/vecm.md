# vecm

## Signature

```d
VecmFit vecm(Y)(Y y, int order, int rank, Deterministic deterministic = Deterministic.unrestrictedConstant, GretlMatrix exog = GretlMatrix.init,
                bool seasonals = false, BetaRestriction restriction = BetaRestriction.init, AlphaRestriction alphaRestriction = AlphaRestriction.init,
                string[] names = null, string[] exogNames = null, int frequency = 1, string start = "1", string fitStart = null, string fitEnd = null);
VecmFit vecm(Dataset d, string[] endog, int order, int rank, Deterministic deterministic = Deterministic.unrestrictedConstant, string[] exog = null,
             bool seasonals = false, BetaRestriction restriction = BetaRestriction.init, AlphaRestriction alphaRestriction = AlphaRestriction.init,
             int frequency = 0, string start = null, string fitStart = null, string fitEnd = null);

struct BetaRestriction { GretlMatrix R, q; }
struct AlphaRestriction { GretlMatrix R, q; }

struct VecmFit {
    string[] names, coefNames, betaNames;
    size_t neqns, rank, order, nobs;
    int df;
    GretlMatrix beta, alpha, betaSe, alphaSe, pi;
    GretlMatrix coef, se, sigma, resid, eigenvalues;
    double lnL, aic, bic, hqc;
    double restrictionLR;
    int restrictionDf;
    string sampleStart, sampleEnd;
    Estimate equation(size_t i);
}
```

## Parameters

- `y`: the data of the series as the columns of a matrix, at least two. No missing values.
- `order`: the lag order of the VAR in levels, at least 1.
- `rank`: the cointegration rank, from 1 to the number of series minus 1.
- `deterministic`, `exog`, `seasonals`, `names`, `exogNames`, `frequency`, `start`, `fitStart`, `fitEnd`: as for `johansen`.
- `restriction`: a `BetaRestriction` on the cointegrating vectors, a homogeneous linear restriction `R beta = 0` that applies to every cointegrating vector (every column of `beta`): `R = [1 1]` says that each is proportional to `(1, -1)`. `beta` is then estimated subject to it and `restrictionLR` tests it. A `q` that is not zero is refused. The default is none.
- `alphaRestriction`: an `AlphaRestriction` on the adjustment coefficients, `R alpha = 0` for every column of `alpha` (for example `R = [1 0]`: the first variable does not adjust, i.e. is weakly exogenous). Not allowed together with a restriction on `beta`. The default is none.

## Return value

A `VecmFit` struct, the vector error-correction model `dy(t) = alpha beta' y(t - 1) + Gamma_1 dy(t - 1) + ... + deterministic terms + e(t)` estimated by Johansen's maximum likelihood. Its fields:

- `names`, `coefNames` (the rows of `coef`: `const`, `d_y1_1`, ..., the exogenous variables, `EC1`, ...) and `betaNames` (the rows of `beta`: the variables, then `const` or `trend` if restricted to the cointegrating space).
- `beta`: the cointegrating vectors (normalized unless restricted) and `alpha`: the adjustment coefficients (`neqns` x `rank`), with standard errors `betaSe` (0 for the normalizing elements) and `alphaSe`, and `pi = alpha beta'`, the long-run matrix.
- `coef`, `se`, `sigma`, `resid`: the coefficients and standard errors of the equations given `beta` (`ncoef` x `neqns`), the residual covariance matrix and the residuals.
- `eigenvalues` (empty when `beta` is restricted), `lnL`, `aic`, `bic`, `hqc`, `restrictionLR` and `restrictionDf` (the likelihood-ratio test of a restriction against the unrestricted model), `sampleStart`, `sampleEnd`.
- `equation(i)`: equation `i` of the error-correction form as an `Estimate` (the dependent variable is the difference; `vcov` is empty).

Like a `VarFit` it keeps the library's estimate, so `irf`, `fevd` and `forecast` work on it. `testBeta` and `testAlpha` test restrictions on a fit.

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
    auto walk = arSeries(n, 1.0);
    auto e = arSeries(n, 0.3);
    auto y1 = walk.dup;
    auto y2 = new double[n];
    foreach (t; 0 .. n) {
        y2[t] = walk[t] + e[t];
    }

    auto fit = vecm(cbind(y1, y2), 2, 1);
    writeln(fit.betaNames);
    writefln("beta: (%.3f, %.3f)", fit.beta[0, 0], fit.beta[1, 0]);        // proportional to (1, -1)
    writefln("alpha: (%.3f, %.3f)", fit.alpha[0, 0], fit.alpha[1, 0]);
    show(fit.pi);

    auto r = BetaRestriction(cbind([1.0], [1.0]));  // R = [1 1]: beta proportional to (1, -1)
    auto restricted = vecm(cbind(y1, y2), 2, 1, restriction: r);
    writefln("LR test of the restriction: %.3f with %d df", restricted.restrictionLR, restricted.restrictionDf);
}
```

Output:

```
["y1", "y2"]
beta: (1.000, -1.001)
alpha: (-0.046, 0.810)
   -0.0460    0.0460
    0.8099   -0.8105
LR test of the restriction: 0.001 with 1 df
```
