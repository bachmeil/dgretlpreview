# VarFit

## Signature

```d
struct VarFit {
    string[] names;                // the endogenous variables
    string[] coefNames;            // the rows of coef
    int[] lags;                    // the lags used
    size_t order, nobs, neqns;
    int df;
    GretlMatrix coef, se;          // ncoef x neqns
    GretlMatrix sigma;             // neqns x neqns
    GretlMatrix resid;             // nobs x neqns
    GretlMatrix roots;             // neqns * order x 3
    GretlMatrix grangerF, grangerP;  // neqns x neqns
    double lnL, aic, bic, hqc, ldet;
    string sampleStart, sampleEnd;

    Estimate equation(size_t i);
    string toString();
}
```

## Parameters

- `names`: the names of the endogenous variables, in the order of the columns of `y`.
- `coefNames`: the names of the rows of `coef` and `se`: `const`, then for each variable its lags (`y1_1`, `y1_2`, `y2_1`, ...), then the exogenous variables, the seasonal dummies (`S1`, ...) and `time` if present.
- `lags`: the lags used (1 to `order` unless a gappy list was given); `order`: the largest lag; `nobs`: the observations used by each equation; `neqns`: the number of equations; `df`: the residual degrees of freedom of each equation.
- `coef`, `se`: `ncoef` x `neqns` matrices of the coefficients and standard errors (robust, if requested), one column for each equation.
- `sigma`: the `neqns` x `neqns` covariance matrix of the residuals; `resid`: the `nobs` x `neqns` residuals.
- `roots`: a `neqns * order` x 3 matrix of the eigenvalues of the companion matrix: real part, imaginary part and modulus. The VAR is stable if all the moduli are less than 1.
- `grangerF`, `grangerP`: `neqns` x `neqns` matrices; element `[i, j]` is the F statistic (and its p-value) of the test that all the lags of variable `j` are zero in equation `i`. Its diagonal is NaN when the test is not available.
- `lnL`, `aic`, `bic`, `hqc`: the log-likelihood and the information criteria of the library (per observation); `ldet`: the log determinant of `sigma`.
- `sampleStart`, `sampleEnd`: labels of the first and last observation used.
- `equation(i)`: equation `i` (counting from 0) as an `Estimate`, with the coefficients, standard errors, t-ratios and R squared of that regression.
- `toString()`: a printable summary (`writeln(fit)`).

## Return value

A `VarFit` is what `var` returns. It holds the results as plain matrices and numbers, and also keeps the library's own estimate (reference counted), which is what `irf`, `fevd` and `forecast` use; a `VecmFit` has the same handle, so those functions work for it too. Copies share the data, which is freed when the last copy goes away. See `var` for how to estimate one and `lagSelection` to choose the order.

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

    writeln(fit.names, " ", fit.coefNames);
    writefln("modulus of the roots: %.3f and %.3f (below 1: stable)", fit.roots[0, 2], fit.roots[1, 2]);
    writefln("Granger: do the lags of y2 help to predict y1? F = %.3f, p = %.4f", fit.grangerF[0, 1], fit.grangerP[0, 1]);
    writefln("standard errors of the first equation: %.3f %.3f %.3f", fit.se[0, 0], fit.se[1, 0], fit.se[2, 0]);
    show(fit.sigma);

    auto eq = fit.equation(0);                       // the first equation as a regression
    writefln("%s: R squared %.3f", eq.method, eq.rsq);
}
```

Output:

```
["y1", "y2"] ["const", "y1_1", "y2_1"]
modulus of the roots: 0.706 and 0.224 (below 1: stable)
Granger: do the lags of y2 help to predict y1? F = 10.053, p = 0.0018
standard errors of the first equation: 0.074 0.058 0.060
    1.0186   -0.0414
   -0.0414    1.1062
VAR equation for y1: R squared 0.384
```
