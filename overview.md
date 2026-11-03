# Overview of dgretl

This file contains the primary functions of interest when learning the dgretl library.

Everything below comes from one import, `import dgretl;`. There is nothing to start or shut down: the library is
initialised before `main` runs and cleaned up after it returns, data live in reference-counted structs that free
themselves, and failures are `GretlException`s that carry the library's message.

## Loading Data

You can load .csv files.

`Dataset importCsv(string fname, string[] options...)`  
`GretlMatrix importCsvAsMatrix(string fname)`

Example usage:

```d
import dgretl;
import std.stdio: writeln;

void main() {
    auto d = importCsv("example.csv");
    writeln(d.nobs, " observations: ", d.names);
    writeln("gdp: ", d["gdp"][0 .. 2]);

    auto m = importCsvAsMatrix("example.csv");
    writeln(m.nrow, " x ", m.ncol, ", m[2, 1] = ", m[2, 1]);
}
```

For a file `example.csv` that contains

```
year,gdp,unemp
2020,1.5,4.1
2021,2.5,3.9
2022,3.1,3.6
2023,2.2,3.7
```

the program prints

```
4 observations: ["gdp", "unemp"]
gdp: [1.5, 2.5]
4 x 3, m[2, 1] = 3.1
```

Things to notice:

- `importCsv` returns a `Dataset`: one series per column, one observation per row. `d.nobs` and `d.nseries` are its
  size, `d.names` the series names, and `d["gdp"]` (by name) or `d[0]` (by number, from 0) a `const(double)[]` with a
  series' values, valid as long as the `Dataset` or a copy of it is alive. Copies share the data, which is freed
  when the last copy goes away.
- The options are optional strings after the file name, in any number: `importCsv("x.csv", "noHeader", "allColumns")`.
  `"noHeader"` is for a file without a header row (the series are then called `v1`, `v2`, ...), `"allColumns"` keeps
  a first column of labels or years as a series, `"codedStrings"` turns columns of text into coded values, and
  `"verbose"` shows the importer's progress messages, which are otherwise kept back. A name that is not one of these
  stops the program with a message (an `assert`, so only in builds that keep assertions).
- When the options are known when you write the program, you can give them at compile time:
  `importCsv!("noHeader", "allColumns")("x.csv")`. The compiler works out the flags, a misspelt option name is then
  a compile-time error (`"noHeaders" is not an option for function importCsv.`), and nothing is looked up when the
  call runs. For `importCsv` the difference is negligible next to reading a file; the form is there for the cheap
  functions that follow the same pattern, which may be called thousands of times in a bootstrap or a simulation.
- A first column such as `year` that looks like observation labels is taken as the observation index (here it
  also identifies the data as annual time series) and does not become a series. That is why `d` has two series.
- `importCsvAsMatrix` returns a `GretlMatrix`. It keeps every column, so the `year` column is column 0. Indices
  are 0-based: `m[i, j]` is row `i`, column `j`, and `m[i, j] = x` assigns. Elements are stored column by column
  underneath.
- `m.nrow` and `m.ncol` are the sizes. `m[a .. b, c .. d]` makes a **copy** of that block, so passing it to a function
  that changes it cannot change `m`. To work on part of a matrix in place, ask for a reference, a `Submatrix`:
  `m.sub[a .. b, c .. d]`, `m.rows[a .. b]` (those rows) or `m.cols[a .. b]` (those columns). It reads and writes the
  elements of `m`, keeps `m`'s data alive, and has its own `sub`, `rows` and `cols`.
- `m.byRow` and `m.byColumn` iterate over the rows or columns (like `byKey` for an associative array), giving a
  `Submatrix` reference for each: `foreach (row; m.byRow) row[0, 1] = 0;` changes `m`. Both also work on a `Submatrix`.
- `Row(m, i)` and `Col(m, j)` give one row or column as a range of doubles that behaves like a `double[]`:
  `foreach (val; Row(m, 1))`, `foreach (j, val; Col(m, 2))`, and `foreach (ref val; Row(m, 1)) val *= 2;` to change `m`.
  They can be indexed, have `length`, work with `std.algorithm`, and `.dup` gives a `double[]` copy.
- `m.el` iterates over all the elements, in no particular order: `foreach (val; m.el)`, `foreach (ref val; m.el) val *= 2;`
  (changes `m`), and `foreach (el, val; m.el)` where `el.row` and `el.col` are the `size_t` position of `val`.
- Operators give a new `GretlMatrix`; `m` and `mm` may each be a `GretlMatrix` or a `Submatrix`, `a` a scalar:
  `m + mm`, `m - mm`, `m + a`, `a - m` and so on; `m * a`, `a * m`, `m / a`, and `a / m` (a divided by each element);
  `m * mm` is matrix multiplication. Element-by-element products and quotients are written with `.el` on either or both
  sides: `m.el * mm`, `m * mm.el`, `m.el * mm.el`, and likewise for `/`. `m + mm.el` and the like are the same as without `.el`. `m / mm` is Matlab's matrix right division, the `x`
  with `x * mm = m` (`m * inv(mm)` for a square `mm`, the QR least-squares solution when `mm` has more columns than rows), so
  for a row `y` (1 x n) and `X'` (p x n), `y / X'` is the row of OLS coefficients. Sizes that do not
  fit throw `GretlException`.
- `solve(X, y)` is Matlab's `X \ y` for data in columns: `X` is n x p, `y` is n x 1 (or n x m), and the result is the p x 1 (p x m)
  matrix of least-squares (OLS) coefficients, computed by QR. It needs n >= p and full-rank `X`, else `GretlException`.
- `t(m)` is the transpose, `inv(m)` the inverse (a singular or non-square matrix throws `GretlException`), and `kron(a, b)`
  (or `a.kron(b)`) the Kronecker product. Each returns a new `GretlMatrix` and accepts a `GretlMatrix` or a `Submatrix`.
- Powers: `m.el ^^ 2` squares each element (`a ^^ m.el` and `m.el ^^ mm` work too), and `power(m, 2)` is `m * m`;
  `power(m, n)` is the matrix power for an integer `n` (0 gives the identity, a negative `n` uses `inv(m)`). `m ^^ 2` on its
  own does not compile, because it could mean either.
- `chol(m)` (or `m.chol`) is the Cholesky factor of a symmetric positive definite matrix. It is the **upper** triangle `U`
  with `t(U) * U == m` (as in R and Matlab); `chol(m, Triangle.lower)` is the **lower** triangle `L` with `L * t(L) == m`.
  Only the lower triangle of `m` is read, and a matrix that is not positive definite throws `GretlException`.
- `diag(m)` is a reference to the diagonal, like `Row` and `Col`, and assigning to it writes into the matrix:
  `diag(m) = [1.1, 3.4, 1.7]`, `diag(m) = Row(m, 1)`, `diag(m) = 0`. (`Row(m, i) = [..]` and `Col(m, j) = [..]` work the
  same way.) The number of values must match, otherwise `GretlException`.
- `eigen(m)` (or `m.eigen`) is the eigen decomposition of a real symmetric matrix. It returns an `Eigen` with `values`
  (an n x 1 matrix, ascending) and `vectors` (the orthonormal eigenvectors as columns), so `m * vectors` equals
  `vectors * diag(values)`; `eigen(m, false)` leaves out the vectors. A non-symmetric matrix throws `GretlException`
  (the general case with complex eigenvalues is not available yet), and the signs of the vectors are LAPACK's.
- `eigenGeneral(m)` (or `m.eigenGeneral`) does the same for any real square matrix, whose eigenvalues can be complex. A
  `GretlMatrix` holds only real numbers, so the result has `values` as an n x 2 matrix (real parts, imaginary parts) and the
  eigenvectors as `vectorsRe` and `vectorsIm` (n x n each); `isReal` says whether all eigenvalues are real. The order is
  LAPACK's (not sorted, conjugate pairs together). Use `eigen` for symmetric matrices.
- `svd(m)` (or `m.svd`) is the singular value decomposition `m = u * diag(s) * t(v)`. It returns an `SVD` with `s` (a
  k x 1 matrix, k = min(nrow, ncol), largest first), `u` (nrow x k) and `v` (ncol x k, V itself and not its transpose).
  `svd(m, false)` gives only `s`, and `svd(m, true, true)` the full-size `u` and `v` for a non-square `m`.
- `qr(m)` (or `m.qr`) is the thin QR decomposition `m = q * r`: a `QR` with `q` (the shape of `m`, orthonormal columns) and
  `r` (ncol x ncol, upper triangular). `m` needs at least as many rows as columns, otherwise `GretlException`.
- `lu(m)` (or `m.lu`) is the LU decomposition with partial pivoting, `p * m = l * u`, for a matrix of any shape. The `LU`
  result has `l` (unit lower triangular), `u` (upper triangular), the permutation matrix `p`, the same permutation as 0-based
  row indices `perm` (row `i` of `p * m` is row `perm[i]` of `m`), and `singular`, true if `u` has an exact zero on its diagonal.
- `schur(m)` (or `m.schur`) is the real Schur decomposition `m = z * t * t(z)` of a square matrix. The `Schur` result has `t`
  (upper triangular except for 2 x 2 blocks on the diagonal, one per complex conjugate pair of eigenvalues), the orthogonal
  `z`, and `values` (n x 2, real and imaginary parts of the eigenvalues in the order of the diagonal of `t`).
  `schur(m, false)` leaves out `z`.
- `pinv(m)` (or `m.pinv`) is the Moore-Penrose pseudoinverse (ncol x nrow), computed from the SVD. Singular values not larger
  than `tol * max(nrow, ncol) * (largest)` count as zero (`pinv(m, tol)`; the default is machine epsilon, as in numpy and
  Matlab). It equals `inv(m)` for a nonsingular matrix and `inv(t(m) * m) * t(m)` for a tall matrix of full column rank.
- `det(m)`, `logdet(m)` and `trace(m)` (or `m.det` and so on) take a square matrix and return a `double`. `logdet` is the log of
  the determinant from the LU factors, so it stays finite when `det` would overflow; a zero determinant gives `-infinity` and a
  negative one `NaN`, unless you ask for the log of the absolute value with `logdet(m, true)`. A non-square matrix throws
  `GretlException`.
- `expm(m)` (or `m.expm`) is the matrix exponential `I + m + m^2/2! + ...` of a square matrix, as a new `GretlMatrix`. (It is
  `expm` and not `exp` so that it does not clash with `std.math.exp`.) `det(expm(m))` equals `exp(trace(m))`.
- `rank(m)` (or `m.rank`) is the numerical rank: the number of singular values larger than `tol * max(nrow, ncol) * (largest)`,
  with `rank(m, tol)` and machine epsilon as the default (as in numpy and Matlab). It returns a `size_t`.
- `eye(n)` is the n x n identity matrix, and `eye(nr, nc)` has ones on the diagonal of an nr x nc matrix.
- `norm(m, kind)` (or `m.norm(kind)`) is a matrix norm as a `double`, with `kind` one of `Norm.frobenius` (the default),
  `Norm.one` (largest column sum of absolute values), `Norm.infinity` (largest row sum), `Norm.spectral` (the 2-norm, the largest
  singular value, which is what Matlab's `norm(A)` is) and `Norm.max` (largest absolute element). For a vector, `frobenius` and
  `spectral` are the Euclidean length.
- `cond(m)` (or `m.cond`) is the condition number: by default (`Norm.spectral`) the largest singular value over the smallest, as in
  Matlab, which also works for a non-square matrix. `cond(m, Norm.one)`, `Norm.infinity` and `Norm.frobenius` are
  `norm(m, kind) * norm(inv(m), kind)` for a square matrix. A singular matrix gives `infinity`.
- `sqrtm(m)` (or `m.sqrtm`) is the principal matrix square root, the `s` with `s * s == m`. (It is `sqrtm` and not `sqrt`, which
  would suggest the elements: for those write `m.el ^^ 0.5`.) A symmetric `m` must be positive semidefinite and gives the
  symmetric root through `eigen`. Any other square `m` is done by the Denman-Beavers iteration and needs to be nonsingular and
  to have no eigenvalue on the negative real axis; if there is no real root (`[1 2; 3 4]` has a negative eigenvalue) it throws
  `GretlException`. (D reads `2 * m.el * mm` as `(2 * m.el) * mm`, so parenthesise as needed.)
- `logm(m)` (or `m.logm`) is the principal matrix logarithm, the inverse of `expm`: `expm(logm(m))` is `m`. A symmetric `m` must be
  positive definite; any other `m` must be nonsingular with no eigenvalue on the negative real axis (it is computed by repeated square
  roots and a Padé-type approximation), otherwise `GretlException`. There is no function for the elementwise logarithm.
- `order(m, col, ...)` and `sort(m, col, ...)` (or `m.order(1, 0)`, `m.sort(1, 0)`) put the rows in ascending order of column `col`;
  where rows tie on it the next column given decides, and so on, and rows that tie on all of them keep their order. `order`
  returns the row indices as a `size_t[]` (as in R) and `sort` a new matrix with the rows in that order. NaN sorts last. The
  sorting is done by `std.algorithm.sorting`, and `import std.algorithm;` together with `import dgretl;` is fine.
- `cumsum(m)` and `cumprod(m)` (or `m.cumsum`) give the cumulative sums and products as a new matrix of the same size, down each column
  by default. `cumsum(m, Along.row)` accumulates across each row instead. A NaN carries along the rest of its column or row.
- `cbind(a, b, ...)` and `rbind(a, b, ...)` join matrices side by side or one below the other, as in R. Each argument is a
  matrix or view, a number (repeated to fill a column or row: `cbind(1, m)` adds a column of ones), or an array of numbers
  (a column for `cbind`, a row for `rbind`; a series from a `Dataset` works). Sizes that do not fit throw `GretlException`. A
  `GretlMatrix.init` with no data is left out, so `acc = rbind(acc, row)` can start from an empty `acc`.
- `m ~ x` adds the columns of `x` to `m` (it is `cbind(m, x)`), and `m ~= x` does it to `m` itself: the storage grows in place
  when nothing else refers to the data of `m`, and otherwise `m` becomes a new matrix, so other copies and views are never disturbed.
- If a file cannot be read, both functions throw `GretlException`, whose message has the file name and the
  library's explanation, for example `importCsv("nofile.csv") failed: ...  Couldn't open nofile.csv`. An uncaught
  one ends the program with that message.
- Your program may have types of its own called `Dataset` or `GretlMatrix`. Import selectively
  (`import dgretl: importCsv;`) or use the qualified names (`dgretl.Dataset`).
- To build a program like this one, see `guide/build.sh` and `guide/loading-and-plotting.md`.

### Where things are

The source tree has three layers:

- `gretl/` (and `gretl/plugin/`): the direct ports of the files of libgretl, one for one next to the original C
  files, as the modules `gretl.<file>`, for example `gretl.csvdata` or `gretl.gretl_matrix`. They keep the C names
  and the C-style interface.
- `cephes/`: the port of the cephes probability and special-function library (`cephes.<file>`).
- `dgretl/`: the D-friendly front end. `import dgretl;` brings in `dgretl.dataset` (`Dataset`), `dgretl.matrix`
  (`GretlMatrix`), `dgretl.importers` (the three import functions), `dgretl.exception` (`GretlException`) and
  `dgretl.startup` (the automatic start-up and shut-down). The two type modules also import the related modules of
  the direct port publicly (`gretl.dataset` and `gretl.dataseth`, `gretl.gretl_matrix` and `gretl.gretl_matrixh`).

### Lower level

`Dataset` and `GretlMatrix` wrap the C-layout structs `DatasetRaw` (gretl.libgretlh) and `GretlMatrixRaw`
(gretl.gretl_matrixh) that the rest of the library works with. The `ptr` property of each wrapper gives the raw
pointer for calling the lower-level functions, and `Dataset(raw)` hands a raw dataset you made yourself to a
`Dataset` that frees it (`Dataset(nvar, nobs)` makes a new one). The functions behind the two above are

`int import_csv(const(char) * fname, DatasetRaw * dset, Gretlopt opt, Prn * prn)`  (gretl.csvdata)  
`GretlMatrixRaw * import_csv_as_matrix(const(char) * fname, int * err)`  (gretl.csvdata)

and they work as in the C library: you create the dataset and the `Prn` yourself, free them yourself, and call
`libgretl_init()` / `libgretl_cleanup()` if `dgretl` is not imported.

## Estimation

`dgretl.estimate` fits models on matrices (or arrays) and returns an `Estimate`, a plain struct with the results as `GretlMatrix`
fields and numbers; there is nothing to free. A constant is added unless you say `constant: false`.

```d
auto fit = ols(y, X);                                        // classical standard errors
auto rob = ols(y, X, vcov: Vcov.hc1);                        // robust: hc0 .. hc3, or Vcov.hac (with lag:, kernel:, prewhiten:)
auto cl  = ols(y, X, cluster: groupId);                      // cluster-robust
writeln(fit);                                                // a table of coefficients
fit.coef, fit.se, fit.tstat, fit.pvalue, fit.vcov, fit.resid, fit.fitted, fit.rsq, fit.sigma, fit.aic, ...
```

- `wls(y, X, w)` weights by `w`, the inverse of each observation's error variance; a series of only 0s and 1s leaves out the 0s.
- `tsls(y, endog, instruments, exog: W)` is two-stage least squares: `instruments` are the *excluded* instruments, and the
  constant and `exog` are used as instruments too. Coefficients come out as constant, `exog`, `endog`.
- `logit(y, X)` and `probit(y, X)` (y is 0 or 1) are estimated by maximum likelihood; `robust: true` gives sandwich standard errors.
  The z-ratios have normal p-values and `fitted` holds the probabilities.
- All of them take a `GretlMatrix`, a `Submatrix` or an array for `y` and `X`, throw `GretlException` for missing values or sizes that
  do not fit, and leave the library settings as they were. The call costs about 10 microseconds for 100 observations, so it is fine inside
  a bootstrap or Monte Carlo loop.

With a `Dataset` you name the series instead:

```d
auto d = importCsv("macro.csv");
auto fit = ols(d, "gdp", ["unemp", "rate"]);                       // coefficients are called const, unemp, rate
auto rob = ols(d, "gdp", ["unemp", "rate"], vcov: Vcov.hc1, cluster: "state");
auto iv  = tsls(d, "ye", ["xe"], ["zin"], exog: ["w"]);
auto lg  = logit(d, "default", ["income", "age"]);
```

`wls(d, y, x, "weights")`, `probit` and the options are the same as for matrices; `cluster` is the name of a series here. Observations with a missing value
in any series used are left out, and `nobs` says how many were used.

## Time series

`dgretl.timeseries` (also part of `import dgretl;`) takes a series as an array, a `GretlMatrix` or one column of a matrix, and returns plain structs; missing
values throw `GretlException`. Functions that need a calendar take `frequency:` (observations per year: 1, 4, 12, ...) and `start:` (`"1990:1"`).

```d
auto c  = correlogram(x, 12);                       // c.acf, c.pacf, c.ljungBox: 12 x 1
auto lb = ljungBox(resid, 8, fitted: 2);            // .statistic, .pvalue, .df
auto t  = adf(x, lags: 2, trend: Trend.constant);   // .statistic, .pvalue, .report (also testDown:, gls:, seasonals:); kpss(x, 4); engleGranger(y, X, 1)
auto lm = archTest(resid, 4);                       // Engle's LM test
auto trend = hpTrend(x, 1600);  auto cycle = hpCycle(x);  bkFilter(x, 6, 32, 12);  butterworth(x, 2, 67.5);  fracDiff(x, 0.4);
auto l = lag(x, 1);  auto d = diff(x);              // NaN-padded, same length
auto p = periodogram(x);  periodogram(x, 8);        // p.frequency, p.period, p.density; a positive window smooths with the Bartlett window
auto fi = fractionalIntegration(x);                 // local Whittle and GPH estimates of d with standard errors and p-values
auto a = ar1(y, X);  auto g = ar(y, X, [1, 4]);  auto h = arch(y, X, 2);   // Estimate; a.arcoef, a.arse hold the error coefficients
```

A VAR keeps the library's structure inside a counted reference, so the analysis calls can use it:

```d
auto fit = var(Y, 2, frequency: 4, start: "1990:1");     // also exog:, constant:, trend:, seasonals:, lags: [1, 4], robust:, fitStart:, fitEnd:
fit.coef, fit.se, fit.sigma, fit.resid, fit.roots, fit.grangerF, fit.grangerP, fit.aic, fit.bic, fit.lnL;   // coef is ncoef x neqns
fit.equation(0);                                          // one equation as an Estimate
auto ls = lagSelection(Y, 8);                             // ls.table, ls.bestBic, ...
auto r  = irf(fit, 12, alpha: 0.1, seed: 42);             // r.at(target, shock): 12 x 1, or 12 x 3 (point, lower, upper) with bootstrap bands
auto v  = fevd(fit, 12);                                  // v.at(target, shock): shares of the forecast error variance
auto f  = forecast(fit, 4);                               // f.mean, f.se: 4 x neqns (futureExog: for models with exogenous variables)
auto g  = forecastRange(fit, 80, 103);                    // observations 80 .. 103: one-step inside the sample, dynamic after it (also method: ForecastMethod.dynamic); f.first = 80
```

The calendar and sample arguments take dates in that calendar (`fitStart: "1992:1"`); `var(d, ["y1", "y2"], 2)` does the same for the series of a `Dataset` (using its calendar).
ARMA and GARCH return an `Estimate` with extras (the fields of `Estimate` are reached directly):

```d
auto m = arma(y, 1, 1);                                   // ARMA(1,1); arma(y, 1, 0, 1) is ARIMA(1,1,0); seasonal: [1, 0, 0], frequency: 4
m.names, m.coef, m.se, m.lnL, m.aic, m.resid, m.arRoots, m.maRoots;   // also exog:, constant:, method: ArmaMethod.kalman / conditional, opg:
auto fc = forecast(m, 8);                                 // fc.mean, fc.se: 8 x 1 (ARIMA: levels; ARMAX: futureExog:)
auto g = garch(r, 1, 1, exog: X);                         // GARCH(1,1) with regressors; g.variance is the conditional variance, g.lrStat the LR test
auto gf = forecast(g, 5, futureX);                        // gf.mean, gf.se = the conditional standard deviation forecast
arma(d, "price", 1, 1);  garch(d, "ret", 1, 1, ["x"]);    // the same on the series of a Dataset
```

Cointegration:

```d
auto jt = johansen(Y, 2);                                  // jt.trace, jt.lmax, jt.tracePvalue, jt.lmaxPvalue, jt.eigenvalues, jt.report; Deterministic.restrictedConstant etc.
auto v  = vecm(Y, 2, 1);                                   // lag order 2 in levels, rank 1: v.beta, v.alpha, v.betaSe, v.pi, v.coef, v.sigma, v.lnL, v.equation(0)
auto r  = irf(v, 12, alpha: 0.1);  fevd(v, 12);  forecast(v, 4);        // the same analysis functions as for a VAR
auto t  = testBeta(v, BetaRestriction(R));                 // LR test of R beta = 0 (every cointegrating vector); vecm(Y, 2, 1, restriction: BetaRestriction(R)) estimates under it
auto ta = testAlpha(v, AlphaRestriction(Ra));             // R alpha = 0, e.g. Ra = [1 0]: the first variable is weakly exogenous; vecm(..., alphaRestriction: ...)
```

Every test, filter and model above also has a form for a `Dataset` that names the series: `adf(d, "gdp", 2)`, `hpTrend(d, "gdp")`, `ar1(d, "y", ["x1", "x2"])`, `lagSelection(d, ["y1", "y2"], 8)`, `engleGranger(d, "y", ["x"])`, `lag(d, ["a", "b"])` (a missing value in a series throws).

Not wrapped yet: restrictions with a right-hand side q different from 0, and restrictions on beta and alpha together.

## Descriptive statistics

`dgretl.descriptive` (also part of `import dgretl;`) works on an array, a `GretlMatrix`/`Submatrix` or one row or column of a matrix (`Col(M, 1)`). For a series a statistic is a `double`; for a matrix it is a 1 x k
`GretlMatrix`, one value per column. NaN means missing and is skipped.

```d
double m = mean(x);  median(x);  variance(x);  stddev(x, population: true);  skewness(x);  kurtosis(x);   // kurtosis is the EXCESS kurtosis
quantile(x, 0.9);  quantile(x, 0.9, QuantileType.q7);  quantiles(x, [0.25, 0.5, 0.75]);  fiveNumbers(x);
total(x);  minimum(x);  maximum(x);  validCount(x);  gini(x);  longRunVariance(x, 4);
auto means = mean(M);                                   // 1 x k, the columns' means
covariance(x, y);  correlation(x, y);  correlation(M);  covariance(M);    // k x k for a matrix, pairwise deletion
auto s = summary(M, names: ["a", "b"], weights: w);     // s.mean, s.median, s.sd, s.min, s.max, s.cv, s.skewness, s.excessKurtosis, s.p5, s.p95, s.iqr, s.valid, s.missing
auto f = frequencies(x, bins: 9);                       // f.midpoint, f.lower, f.upper, f.count, f.fraction, f.cumulative; normalTest: true; valueCounts(x) for discrete data
auto n = normalityTests(x);                             // n.doornikHansen, n.shapiroWilk, n.jarqueBera, n.lilliefors: .statistic, .pvalue;  multivariateNormality(E)
auto t = meansTest(x, y, MeansTest.unequalVariances);   // .statistic, .pvalue, .df, .difference, .se;  MeansTest.paired;  variancesTest(x, y): .test.pvalue (two-sided), .pvalueUpper
auto d = rowDistances(X, Y, DistanceMetric.manhattan);  rowDistanceMatrix(X);  mahalanobis(X);        // between rows
auto c = crosstab(a, b);                                // c.counts, c.rowValues, c.colValues, c.pearson, c.pvalue
auto st = signTest(x, y);  auto sr = signedRankTest(x, y);  auto rs = rankSumTest(x, y);   // .statistic, .pvalue (two-sided), .pvalueOneSided; the rank sums, W+ and W-, ...
auto cc = crossCorrelogram(x, y, 8);                    // (timeseries) cc.lags, cc.ccf
```

For a `Dataset` the series go straight in (`mean(d["gdp"])`), and `summary`, `covariance`, `correlation`, `mahalanobis`, `crosstab`, `meansTest`, `variancesTest` also take the names.

## Probability distributions

`dgretl.distributions` (part of `import dgretl;`) gives the density (or probability mass), the distribution function and the quantile function of the library's distributions, with the names and argument order of R: `dnorm`, `pnorm`, `qnorm`, and the same for `t`, `chisq`, `F`, `gamma`, `beta`, `exp`, `weibull`, `logis`, `unif`, `laplace`, `ged`, `tnorm` (truncated normal), `binom` and `pois`.

```d
double p = pnorm(1.96);                      // 0.975
double c = qchisq(0.95, 3);                  // 7.81
double u = pt(2.5, 10, lower: false);        // upper tail P(T > 2.5), accurate far out
auto d = dnorm(x, mean: 1, sd: 2);           // x: double[] or a matrix; the same shape comes back
double w = pt(2.0, 8, ncp: 1.5);             // non-central t (also chi-square and F)
auto pit = pnorm(d, "y", "mu", "sigma");  // parameters can vary by element: arrays, matrices or, for a Dataset, series names
auto z = pnorm(x, mu, 2.0);                  // arrays and numbers mix; mu: double[] with one value per element of x
auto s = pnorm(d, "z", mean: 1, sd: 2);     // a series of a Dataset, by name: a double[] with one value per observation
```

A quantile of a discrete distribution is the smallest x with P(X <= x) >= p. Missing values and probabilities outside [0, 1] give NaN; invalid parameters throw `GretlException`. The F distribution takes integer degrees of freedom (the library's restriction); the chi-square distribution accepts any positive df.

## Random numbers

`dgretl.random` (also part of `import dgretl;`) gives the library's generators as D functions, and parallel streams that drive all of them.

```d
setSeed(42);                                          // the library's generator: sequential, reproducible
auto x = rnorm(1000, mean: 1, sd: 2);                 // double[]; also runif, rgamma, rchisq, rt, rF, rbinom, rpois, rweibull, rexp, rlogis, rged, rlaplace, rbeta, rbetabinom, rtnorm
auto k = runifInt(10, 1, 100, withoutReplacement: true);   rdiscrete(100, [0.2, 0.5, 0.3]);       // int[] / size_t[]
auto d = rdirichlet([1, 2, 3], 500);  auto w = rinvwishart(S, 10);                                  // GretlMatrix
auto z = randomMatrix!rgamma(100, 3, 2.0);            // a 100 x 3 matrix of gamma(2) draws

auto streams = makeStreams(8, seed: 2024);            // 8 non-overlapping streams of L'Ecuyer's MRG31k3p
auto draws = parallelDraws(streams, (size_t i) => rgamma(10_000, 2.0));      // one task per stream, in parallel; the functions draw from that task's stream
auto y = streams[3].rnorm(100);                       // every distribution is also a method of a stream
```

A result depends only on the stream, not on the number of processors or the scheduling, so parallel simulations reproduce from the seed. A stream is used by one thread at a time.
