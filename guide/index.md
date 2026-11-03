# dgretl: what version 1 contains

dgretl is a D front end to the port of libgretl. One import, `import dgretl;`, gives everything below. There is nothing to start or shut down (the library is
initialised before `main` and cleaned up after it), data live in reference-counted structs that free themselves, and every failure is a `GretlException` carrying the
library's message. Results are plain structs of `GretlMatrix` fields and numbers.

Each bullet is a topic for a guide of its own; link the guides from here as they are written.

## Getting started

- [Installing, building and running a first program](installing-and-building.md) (`ldmd2`/`dmd`, the C libraries needed)
- [Conventions used throughout](conventions.md): arrays, `GretlMatrix`, `Submatrix` and rows or columns of a matrix are accepted wherever data go in; missing values are NaN; named arguments for options
- Errors: [`GretlException`](docs/errors/GretlException.md)
- [Loading a spreadsheet and plotting a series](loading-and-plotting.md) (a worked example)

## Data

- [`Dataset`](docs/data/Dataset.md): a reference-counted set of series (`nobs`, `nseries`, `names`, a series by name or number)
- Reading CSV files: [`importCsv`](docs/data/importCsv.md) (header or none, label columns, coded strings, options chosen at compile time), [`importCsvAsMatrix`](docs/data/importCsvAsMatrix.md)
- Reading spreadsheets: [`importSpreadsheet`](docs/data/importSpreadsheet.md) (Excel `.xlsx`/`.xls`, OpenDocument `.ods`, Gnumeric; choice of sheet and offsets)

## Matrices

- [`GretlMatrix`](docs/matrix/GretlMatrix.md): reference-counted, 0-based, column-major; `nrow`, `ncol`, element access and assignment
- Copies and views: slicing copies; [`sub`](docs/matrix/sub.md), [`rows`](docs/matrix/rows.md) and [`cols`](docs/matrix/cols.md) give writable [`Submatrix`](docs/matrix/Submatrix.md) views; [`Row`](docs/matrix/Row.md), [`Col`](docs/matrix/Col.md) and [`diag`](docs/matrix/diag.md) give lines that can be assigned to
- Iteration: [`byRow`](docs/matrix/byRow.md), [`byColumn`](docs/matrix/byColumn.md), [`el`](docs/matrix/el.md) (the elements)
- Arithmetic: [`operators`](docs/matrix/operators.md) (`+`, `-`, `*` matrix product, `/` right division, scalars, element-wise operations with `.el`)
- Building: [`eye`](docs/matrix/eye.md), [`cbind`](docs/matrix/cbind.md), [`rbind`](docs/matrix/rbind.md), `~` and `~=`
- Shape and structure: [`t`](docs/matrix/t.md) (transpose), [`diag`](docs/matrix/diag.md), [`sort`](docs/matrix/sort.md), [`order`](docs/matrix/order.md), [`cumsum`](docs/matrix/cumsum.md), [`cumprod`](docs/matrix/cumprod.md)
- Solving and inverting: [`solve`](docs/matrix/solve.md), [`inv`](docs/matrix/inv.md), [`pinv`](docs/matrix/pinv.md), [`power`](docs/matrix/power.md), [`kron`](docs/matrix/kron.md)
- Scalars of a matrix: [`det`](docs/matrix/det.md), [`logdet`](docs/matrix/logdet.md), [`trace`](docs/matrix/trace.md), [`rank`](docs/matrix/rank.md), [`norm`](docs/matrix/norm.md) (several norms), [`cond`](docs/matrix/cond.md)
- Decompositions: [`chol`](docs/matrix/chol.md), [`eigen`](docs/matrix/eigen.md) (symmetric), [`eigenGeneral`](docs/matrix/eigenGeneral.md) (complex eigenvalues), [`svd`](docs/matrix/svd.md), [`qr`](docs/matrix/qr.md), [`lu`](docs/matrix/lu.md), [`schur`](docs/matrix/schur.md)
- Matrix functions: [`expm`](docs/matrix/expm.md), [`sqrtm`](docs/matrix/sqrtm.md), [`logm`](docs/matrix/logm.md)

## Regression

- Least squares: [`ols`](docs/regression/ols.md), [`wls`](docs/regression/wls.md) (weights), [`tsls`](docs/regression/tsls.md) (two-stage least squares with excluded instruments)
- Binary choice: [`logit`](docs/regression/logit.md), [`probit`](docs/regression/probit.md)
- Standard errors: classical, heteroskedasticity-robust (HC0 to HC3), HAC (Newey-West kernels, lag, pre-whitening), cluster-robust (the `vcov` options of [`ols`](docs/regression/ols.md))
- The [`Estimate`](docs/regression/Estimate.md) result: coefficients, standard errors, t-ratios, p-values, covariance matrix, residuals, fitted values, R squared, F, log-likelihood, information criteria, Durbin-Watson
- Regressions on matrices or on named series of a `Dataset`

## Descriptive statistics

- Location, spread and shape: [`mean`](docs/descriptive/mean.md), [`median`](docs/descriptive/median.md), [`total`](docs/descriptive/total.md), [`minimum`](docs/descriptive/minimum.md), [`maximum`](docs/descriptive/maximum.md), [`sst`](docs/descriptive/sst.md), [`variance`](docs/descriptive/variance.md), [`stddev`](docs/descriptive/stddev.md), [`skewness`](docs/descriptive/skewness.md), [`kurtosis`](docs/descriptive/kurtosis.md), [`gini`](docs/descriptive/gini.md), [`longRunVariance`](docs/descriptive/longRunVariance.md), [`validCount`](docs/descriptive/validCount.md)
- Quantiles: [`quantile`](docs/descriptive/quantile.md) (three definitions), [`quantiles`](docs/descriptive/quantiles.md), [`fiveNumbers`](docs/descriptive/fiveNumbers.md)
- Series or whole matrices: the same functions give one value per column of a matrix
- Covariance and correlation: [`covariance`](docs/descriptive/covariance.md) and [`correlation`](docs/descriptive/correlation.md), for two series and as matrices (pairwise deletion of missing values)
- Summary table: [`summary`](docs/descriptive/summary.md) (with weights)
- Frequency distributions: [`frequencies`](docs/descriptive/frequencies.md), [`valueCounts`](docs/descriptive/valueCounts.md), cross-tabulation [`crosstab`](docs/descriptive/crosstab.md) with Pearson's chi-square
- Normality: [`normalityTests`](docs/descriptive/normalityTests.md) (Doornik-Hansen, Shapiro-Wilk, Jarque-Bera, Lilliefors), [`multivariateNormality`](docs/descriptive/multivariateNormality.md)
- Two-sample tests: [`meansTest`](docs/descriptive/meansTest.md) (pooled, Welch, paired), [`variancesTest`](docs/descriptive/variancesTest.md)
- Nonparametric tests: [`signTest`](docs/descriptive/signTest.md), [`signedRankTest`](docs/descriptive/signedRankTest.md), [`rankSumTest`](docs/descriptive/rankSumTest.md) (with small-sample critical values)
- Distances: [`rowDistances`](docs/descriptive/rowDistances.md), [`rowDistanceMatrix`](docs/descriptive/rowDistanceMatrix.md) (six metrics), [`mahalanobis`](docs/descriptive/mahalanobis.md)

## Time series: transformations and diagnostics

- Lags and differences: [`lag`](docs/timeseries/lag.md), [`diff`](docs/timeseries/diff.md) (including seasonal differences)
- Filters: [`hpTrend`](docs/timeseries/hpTrend.md), [`hpCycle`](docs/timeseries/hpCycle.md) (Hodrick-Prescott), [`bkFilter`](docs/timeseries/bkFilter.md) (Baxter-King), [`butterworth`](docs/timeseries/butterworth.md), [`fracDiff`](docs/timeseries/fracDiff.md)
- Correlograms: [`correlogram`](docs/timeseries/correlogram.md) (ACF, PACF, Ljung-Box), [`ljungBox`](docs/timeseries/ljungBox.md), [`crossCorrelogram`](docs/timeseries/crossCorrelogram.md)
- Unit roots and stationarity: [`adf`](docs/timeseries/adf.md) (trends, lag selection, GLS detrending, seasonals), [`kpss`](docs/timeseries/kpss.md), [`engleGranger`](docs/timeseries/engleGranger.md)
- Spectrum and long memory: [`periodogram`](docs/timeseries/periodogram.md), [`fractionalIntegration`](docs/timeseries/fractionalIntegration.md) (local Whittle and GPH)
- Heteroskedasticity: [`archTest`](docs/timeseries/archTest.md)
- Single-equation models with autocorrelated or ARCH errors: [`ar1`](docs/timeseries/ar1.md) (Cochrane-Orcutt, Prais-Winsten, Hildreth-Lu), [`ar`](docs/timeseries/ar.md), [`arch`](docs/timeseries/arch.md)
- Calendars: giving a series a frequency and starting date, and choosing an estimation sample by date
- All of the above also exist for the named series of a `Dataset`

## Time series: models

- Vector autoregression: [`var`](docs/timeseries/var.md) (exogenous variables, trend, seasonals, gappy lags, robust errors), [`VarFit`](docs/timeseries/VarFit.md) (coefficients, standard errors, residual covariance, roots, Granger tests, information criteria, single equations), [`lagSelection`](docs/timeseries/lagSelection.md)
- Analysis of a VAR or VECM: [`irf`](docs/timeseries/irf.md) (with bootstrap bands), [`fevd`](docs/timeseries/fevd.md), [`forecast`](docs/timeseries/forecast.md), [`forecastRange`](docs/timeseries/forecastRange.md) (in-sample, out-of-sample or mixed; static or dynamic)
- Cointegration: [`johansen`](docs/timeseries/johansen.md) (five deterministic cases), [`vecm`](docs/timeseries/vecm.md), restrictions on alpha and beta (`AlphaRestriction`, `BetaRestriction`, [`testAlpha`](docs/timeseries/testAlpha.md), [`testBeta`](docs/timeseries/testBeta.md))
- ARMA and ARIMA: [`arma`](docs/timeseries/arma.md) (seasonal terms, regressors, exact, Kalman or conditional estimation), roots, forecasts
- GARCH: [`garch`](docs/timeseries/garch.md) (robust errors, regressors, two estimators), conditional variance, forecasts

## Probability distributions

- Density (or mass), distribution function and quantile function with R's names (`d`, `p` and `q` followed by the name):
  - normal: [`dnorm`](docs/distributions/dnorm.md), [`pnorm`](docs/distributions/pnorm.md), [`qnorm`](docs/distributions/qnorm.md)
  - Student's t: [`dt`](docs/distributions/dt.md), [`pt`](docs/distributions/pt.md), [`qt`](docs/distributions/qt.md)
  - chi-square: [`dchisq`](docs/distributions/dchisq.md), [`pchisq`](docs/distributions/pchisq.md), [`qchisq`](docs/distributions/qchisq.md)
  - F: [`dF`](docs/distributions/dF.md), [`pF`](docs/distributions/pF.md), [`qF`](docs/distributions/qF.md)
  - gamma: [`dgamma`](docs/distributions/dgamma.md), [`pgamma`](docs/distributions/pgamma.md), [`qgamma`](docs/distributions/qgamma.md)
  - beta: [`dbeta`](docs/distributions/dbeta.md), [`pbeta`](docs/distributions/pbeta.md), [`qbeta`](docs/distributions/qbeta.md)
  - exponential: [`dexp`](docs/distributions/dexp.md), [`pexp`](docs/distributions/pexp.md), [`qexp`](docs/distributions/qexp.md)
  - Weibull: [`dweibull`](docs/distributions/dweibull.md), [`pweibull`](docs/distributions/pweibull.md), [`qweibull`](docs/distributions/qweibull.md)
  - logistic: [`dlogis`](docs/distributions/dlogis.md), [`plogis`](docs/distributions/plogis.md), [`qlogis`](docs/distributions/qlogis.md)
  - uniform: [`dunif`](docs/distributions/dunif.md), [`punif`](docs/distributions/punif.md), [`qunif`](docs/distributions/qunif.md)
  - Laplace: [`dlaplace`](docs/distributions/dlaplace.md), [`plaplace`](docs/distributions/plaplace.md), [`qlaplace`](docs/distributions/qlaplace.md)
  - generalized error: [`dged`](docs/distributions/dged.md), [`pged`](docs/distributions/pged.md), [`qged`](docs/distributions/qged.md)
  - truncated normal: [`dtnorm`](docs/distributions/dtnorm.md), [`ptnorm`](docs/distributions/ptnorm.md), [`qtnorm`](docs/distributions/qtnorm.md)
  - binomial: [`dbinom`](docs/distributions/dbinom.md), [`pbinom`](docs/distributions/pbinom.md), [`qbinom`](docs/distributions/qbinom.md)
  - Poisson: [`dpois`](docs/distributions/dpois.md), [`ppois`](docs/distributions/ppois.md), [`qpois`](docs/distributions/qpois.md)
- Non-central t, chi-square and F (`ncp`); upper tails (`lower: false`) that stay accurate in the far tail
- Scalars, arrays and matrices in, the same shape out; a series of a `Dataset` by name; parameters that vary by observation (arrays, matrices, series names); invalid parameters throw

## Random numbers

- Seeding the library's generator: [`setSeed`](docs/random/setSeed.md), [`getSeed`](docs/random/getSeed.md)
- Distributions: [`runif`](docs/random/runif.md), [`runifInt`](docs/random/runifInt.md), [`rnorm`](docs/random/rnorm.md), [`rtnorm`](docs/random/rtnorm.md), [`rgamma`](docs/random/rgamma.md), [`rchisq`](docs/random/rchisq.md), [`rt`](docs/random/rt.md), [`rF`](docs/random/rF.md), [`rbinom`](docs/random/rbinom.md), [`rbetabinom`](docs/random/rbetabinom.md), [`rpois`](docs/random/rpois.md), [`rweibull`](docs/random/rweibull.md), [`rexp`](docs/random/rexp.md), [`rlogis`](docs/random/rlogis.md), [`rged`](docs/random/rged.md), [`rlaplace`](docs/random/rlaplace.md), [`rbeta`](docs/random/rbeta.md), [`rdiscrete`](docs/random/rdiscrete.md), [`rdirichlet`](docs/random/rdirichlet.md), [`rinvwishart`](docs/random/rinvwishart.md)
- Random matrices: [`randomMatrix`](docs/random/randomMatrix.md)
- Parallel streams: [`makeStreams`](docs/random/makeStreams.md) (L'Ecuyer's MRG31k3p; non-overlapping streams), [`parallelDraws`](docs/random/parallelDraws.md), [`RandomStream`](docs/random/RandomStream.md) (every distribution as a method of a stream)
- Reproducible parallel simulation whatever the number of processors

## Reference

- The direct ports in `gretl/`: one D module for each libgretl file, with the library's own names, for when the front end is not enough
- Tests: what is checked and how (`testing/README.md`)
- Library quirks found while porting: `bugs.md` in the development repository, <https://github.com/bachmeil/gretlport>

## Not in version 1

- Running gretl scripts and user-defined functions (the command interpreter, loops and saved objects are not ported)
- GMM, quantile regression, kernel density estimation, count-data and limited-dependent-variable models other than logit and probit, dynamic panels, stepwise selection, leverage and VIF, principal components, k-means
- Importing EViews, SPSS, SAS, JMulti, JSON, ODBC and binary data files; database access
- TeX output of results
