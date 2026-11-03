# Smoke tests for the D port

`./run_tests.sh [test_name ...]` builds and runs every `test_*.d` in this
directory against all the ported modules in the parent directory. These are
smoke tests: they check that each function compiles, links and produces output,
not that the numbers are right.

## Layout

- `test_adf_kpss.d`: `get_urc_pvalue`, `kpss_critvals`, `adf_test` (several option
  combinations), `kpss_test`. Output of every call is printed so it can be eyeballed.
  `engle_granger_test`; `levin_lin_test` needs a plugin and is not tested.
- `test_qr_estimate.d`: QR/Cholesky regression, robust (HC, HAC, clustered) covariance
  matrices, the HAC helpers. Unlike the other smoke tests it also checks that the QR and
  Cholesky estimators agree.
- `test_boxplots.d`: quartiles, whiskers, outliers, the median interval and summary printing (boxplots.d); full plots are not covered.
- `test_genfuncs.d`: differences, cumulation, standardization, sort-by, ranks, Box-Cox, fractional differencing, EMA, AR filter, HP filter, Bessel/logistic, vector stats, value substitution (genfuncs.d); bandpass/Butterworth, panel, loess, aggregation and interpolation are not covered.
- `test_gretl_color.d`: channel decomposition and colormix_array() string results (gretl_color.d); the gnuplot display of the mixes is not covered.
- `test_nonparam.d`: Spearman, Kendall, sign/signed-rank/rank-sum tests, sort_pairs_by_x, loess, kernel bandwidth (nonparam.d); printing forms, runs and Locke/VGE tests are not covered.
- `test_gretl_typemap.d`: type relations, names, parsing, gen-type keywords and ordering (gretl_typemap.d).
- `test_uservar.d`: named scalars, matrices, strings and lists, lookups, replace and delete (uservar.d); function-level scoping and XML serialisation are not covered.
- `test_gretl_bundle.d`: bundle members of each kind, copy/union/append/equality and printing (gretl_bundle.d); serialisation and Kalman/virtual members are not covered.
- `test_gretl_array.d`: arrays of strings and matrices (gretl_array.d); bundles/lists/arrays of arrays, user-function sorting and serialisation are not covered.
- `test_gretl_plot.d`: the "plot" block state machine (gretl_plot.d); finalizing into gnuplot is not covered.
- `test_gretl_gridplot.d`: grid layout arithmetic and small helpers (gretl_gridplot.d); collecting and finalizing plots are not covered.
- `test_graphing.d`: colours, PNG scale/size and naming helpers (graphing.d); plot construction and gnuplot launching are not covered.
- `test_plotbands.d`: band styles and the gnuplot text for band lines (plotbands.d); plot_with_band and band specifications are not covered.
- `test_plotspec.d`: plot specification management and label/arrow output (plotspec.d); full gnuplot output, fits and bars files are not covered.
- `test_system.d`: building a two-equation system, its accessors and estimation through the sysest plugin module (system.d, sysest.d): OLS and SUR; 3SLS, TSLS, restrictions and FIML/LIML are not covered.
- `test_nls.d`: nonlinear-model specification, parameter lists and the Gauss-Newton regression (nls.d); estimation proper needs the genr parser and is not covered.
- `test_gretl_restrict.d`: Wald/F test of a hand-built restriction set (gretl_restrict.d); parsing, systems, VECMs and nonlinear restrictions are not covered.
- `test_compare.d`: omit/add/RESET/non-linearity/autocorrelation/Chow tests (compare.d); plugin- and gnuplot-based tests are not covered.
- `test_tsls.d`: two-stage least squares against the two-step OLS procedure (tsls.d); LIML, GMM and Hausman tests are not covered.
- `test_subsample.d`: sample restriction by mask, dummy, missing values and at random, restore, set_sample (subsample.d); expressions, panels and permanent samples are not covered.
- `test_irfboot.d`: bootstrap bands for VAR impulse responses (irfboot.d); VECMs are not covered.
- `test_tsjoin.d`: the time-series join helpers (tsjoin.d).
- `test_varprint.d`: VAR text printout, covariance, impulse responses and variance decompositions (varprint.d); VECM and TeX/RTF are not covered.
- `test_vartest.d`: VAR diagnostics and omission tests (vartest.d); only their presence, not their values, is checked.
- `test_var.d`: a simulated VAR(1) through gretl_VAR (var.d); VECMs, lag selection, tests and serialization are not covered.
- `test_printout.d`: printout.d (bufgets family, digits, layout helpers, list_series, printdata); distributions, forecasts and RTF/TeX are not covered.
- `test_printscan.d`: sprintf/printf with literal numeric arguments (printscan.d); arguments needing the genr parser are not covered.
- `test_modelprint.d`: model printing (modelprint.d) of an OLS model in plain text, RTF and CSV; TeX output is not covered.
- `test_geneval.d`: the genr evaluator (geneval1.d - geneval5.d). Syntax trees are built by hand, as the parser
  would build them (the lexer/parser, gensyntax.c and genlex.c, are not ported), and passed to `eval()`: scalar,
  series and matrix arithmetic, comparisons and short-circuit logic, `? :`, one-argument functions on every
  type, series/matrix/string functions, lags, `seq`, `epochday` and a few helpers.
- `stubs.d`: small working versions of libgretl functions that are not ported yet but
  are needed at run time (`series_index`, `ntolabel`, `printmodel`, ...). Delete each
  one when the real module is ported.
- `stubs_unported.d`: generated, plain D (compiles with `dmd` or `ldmd2`), one aborting `extern (C)` stub per still-unported function
  that the ported modules reference. If a test reaches one it aborts with
  `*** call to unported libgretl function NAME()`. Regenerate with
  `./make_unported_stubs.sh` after porting more files.
- `build/`: build output (ignored by git), including `libceph.a`, built from the cephes
  sources of the gretl tree (set `CEPHES_DIR`, default `~/Downloads/gretl-2026c/cephes`).

## Adding a test

Create `test_<name>.d` with a `main()`. It can use the helpers in `test_adf_kpss.d` as a
model (a `check()` function, a deterministic random generator, a dataset builder).
If it hits an aborting stub, either port the missing file or add a minimal stand-in to
`stubs.d`.
- `test_gretl_commands.d`: command word table round trip, aliases, help index and completion (gretl_commands.d).
- `test_usermat.d`: matrix lookup by name, range/element submatrices, replacement, scalar assignment, vec/vech/unvech, in-place ops, matrix OLS, eigenvalues, names (usermat.d); matrix/dummy selection, complex matrices, QR/SVD, RLS and GHK are not covered.
- `test_interact_extra.d`: check_stringvar_name and the plot-produced flag (interact_extra.d).
- `test_random.d`: the generator core, Ziggurat normal and gamma draws checked against reference values from the original C; distribution moments, Halton, Dirichlet, discrete, inverse Wishart and hex strings (random.d). The boxplot and irfboot tests now seed the generator explicitly.
- `test_cephes_const.d`: the cephes constants (machine epsilon, largest number, logs, pi multiples, minus zero), consistency, and a C function (ndtri) reading them through their C names (cephes/constants.d).
- `test_cephes_gamma.d`: gamma and log-gamma values, branches, special values and signs, agreement with libm over grids, and the still-C `igam` calling the D `lgam` (cephes/gamma.d).
- `test_cephes_kv.d`: K_v(x) against tabulated values, half-integer closed forms, scaling, symmetry, recurrence, and the error cases (cephes/kv.d).
- `test_cephes_jv.d`: J_v(x) for integer orders (including orders above 500) against libm jn, half-integer closed forms, the three-term recurrence, reflections, zero arguments, Hankel expansion, and the domain error (cephes/jv.d).
- `test_cephes_hyp2f1.d`: 2F1 against closed forms (log, arcsine, binomial, polynomial cases), Gauss's theorem, the x < -1 and x near 1 transformations, symmetry, Euler's transformation and the divergent case (cephes/hyp2f1.d).
- `test_cephes_airy.d`: Ai, Ai', Bi, Bi' at tabulated points, the Wronskian across all branches and the out-of-range return (cephes/airy.d).
- `test_cephes_hyperg.d`: 1F1 against closed forms (exponential, erf, polynomial), Kummer's transformation, small and large |x|, the singular case, hyp2f0, and the recorded upstream inaccuracy at 1F1(1;2;-3) (cephes/hyperg.d).
- `test_cephes_incbet.d`: the incomplete beta integral against closed forms, reflection, large-argument symmetry, monotonicity, the binomial tail relation, domain cases, and the still-C stdtr calling it (cephes/incbet.d).
- `test_cephes_ndtr.d`: the normal distribution function at tabulated points, against libm erfc over the tails, symmetry, extreme arguments, special values, and the round trip with the still-C ndtri (cephes/ndtr.d).
- `test_cephes_j0.d`: J0 and Y0 at tabulated points, against libm over all regimes, symmetry, the first zero, the tiny-argument branch and the Y0 domain error (cephes/j0.d).
- `test_cephes_fdtr.d`: the F distribution (fdtr, fdtrc, fdtri) against closed forms, tabulated critical values, complement and inversion round trips, and the domain errors (cephes/fdtr.d).
- `test_cephes_incbi.d`: the inverse incomplete beta integral against closed forms, round trips through incbet, reflection, boundaries, and the still-C bdtri calling it (cephes/incbi.d).
- `test_cephes_bdtr.d`: the binomial distribution (bdtr, bdtrc, bdtri) against direct summation, closed forms, complement, inversion round trips and the domain errors (cephes/bdtr.d).
- `test_cephes_j1.d`: J1 and Y1 at tabulated points, against libm, oddness for |x| <= 5, the first zero, the recorded NaN of J1 below -5 and the Y1 domain error (cephes/j1.d).
- `test_cephes_ndtri.d`: the inverse normal distribution function at tabulated quantiles, symmetry, round trips with ndtr from 1e-300 to 1/2, and the domain errors (cephes/ndtri.d).
- `test_cephes_kn.d`: K_n(x) at tabulated points, against the netlib routine (cephes.kv), recurrence, symmetry and the error cases (cephes/kn.d).
- `test_cephes_stdtr.d`: Student's t distribution and inverse against closed forms (k = 1, 2), tabulated critical values, symmetry, large and non-integer degrees of freedom, inversion round trips and the domain errors (cephes/stdtr.d).
- `test_cephes_polrt.d`: polynomial roots (real, complex pair, repeated, zero, degree 8), the return codes, and the use through gretl_matrix_polroots (cephes/polrt.d).
- `test_cephes_struve.d`: the Struve function against numerical integration, Y_v (half-integer closed forms, integer orders), the negative half-integer case, onef2 and threef0, and the recorded failure of H_{1/2}(40) (cephes/struve.d).
- `test_cephes_igam.d`: the incomplete gamma integral (igam, igamc) against closed forms, the Poisson sum, the complement across branches, large arguments, underflow, the domain cases and the round trip with the still-C igami (cephes/igam.d).
- `test_cephes_i1.d`: I1 and i1e against the power series, tabulated values, oddness, the scaling relation and the large-argument form (cephes/i1.d).
- `test_cephes_i0.d`: I0 and i0e against the power series, tabulated values, evenness, the scaling relation, the derivative I0' = I1 and the large-argument form (cephes/i0.d).
- `test_cephes_k0.d`: K0 and k0e against tabulated values, the netlib routine, the scaling relation, the Wronskian with I0, I1 and K1, and the domain errors (cephes/k0.d).
- `test_cephes_k1.d`: K1 and k1e against tabulated values, the netlib routine, the scaling relation, the Wronskian with I0, I1 and K0, and the domain errors (cephes/k1.d).
- `test_cephes_nbdtr.d`: the negative binomial distribution (nbdtr, nbdtrc, nbdtri) against direct summation, the geometric case, complement, inversion round trip and the domain errors (cephes/nbdtr.d).
- `test_cephes_chdtr.d`: the chi-square distribution (chdtr, chdtrc, chdtri) against closed forms, tabulated critical values, complement, inversion round trip and the domain errors (cephes/chdtr.d).
- `test_cephes_mtherr.d`: the cephes error messages (captured from stderr), the recorded code and its clearing, out-of-range codes and the hush switch (cephes/mtherr.d).
- `test_cephes_psi.d`: the digamma function at tabulated values, the recurrence and reflection formulas, large arguments and the poles (cephes/psi.d).
- `test_cephes_pdtr.d`: the Poisson distribution (pdtr, pdtrc, pdtri) against direct summation, the k = 0 form, complement, inversion round trip and the domain errors (cephes/pdtr.d).
- `test_cephes_cbrt.d`: the cephes cube root: cubes, signs, special values, accuracy against the correctly rounded cbrtl and cubing back (cephes/cbrt.d).
- `test_cephes_unity.d`: log(1+x), exp(x)-1 and cos(x)-1 near zero against libm and an identity, continuity across the switch points and special values (cephes/unity.d).
- `test_cephes_iv.d`: I_v(x) against I0/I1, the closed forms for orders 3/2 and 5/2, the recurrence, symmetry rules, special values and the recorded inaccuracy at order 1/2 (cephes/iv.d).
- `test_cephes_jn.d`: J_n(x) against libm jn over x <= 30, parity, the recurrence, the tiny-argument case and the recorded degradation for large x (cephes/jn.d).
- `test_cephes_gdtr.d`: the gamma distribution (gdtr, gdtrc) against closed forms, the chi-square relation, complement and the domain error (cephes/gdtr.d).
- `test_cephes_yn.d`: Y_n(x) against libm yn, tabulated values, parity, the Wronskian with libm jn, the recurrence and the domain error (cephes/yn.d).
- `test_cephes_chbevl.d`: the Chebyshev series evaluation (at x/2) against hand cases and direct summation of Chebyshev polynomials (cephes/chbevl.d).
- `test_cephes_polevl.d`: polevl and p1evl on hand-worked polynomials and against direct power sums (cephes/polevl.d).
- `test_cephes_expx2.d`: exp(x^2) and exp(-x^2) against expl, the gain over the naive formula for large x, evenness, overflow and underflow (cephes/expx2.d).
- `test_cephes_btdtr.d`: btdtr equals incbet, closed forms and the domain case (cephes/btdtr.d).
- `test_cephes_igami.d`: the inverse of the complemented incomplete gamma integral against closed forms, round trips, chi-square critical values, monotonicity and the end points (cephes/igami.d).
- `test_dgretl.d`: the D front end (dgretl/): automatic start-up, `importCsv` and its string options (run time and compile time), `Dataset`/`GretlMatrix` copies, indexing and constructors, `GretlException` messages, spreadsheet errors, idempotent init and a heap check over 300 imports.
- `test_dgretl_matrix.d`: `GretlMatrix` `nrow`/`ncol`, `byRow`/`byColumn` iteration, `Row`/`Col` as ranges of doubles, `el` over all elements, arithmetic operators (including `.el` element-wise forms), `solve`, `t`, `inv`, `kron`, `.el ^^`, `power`, `chol`, `diag` and assignment to rows, columns and the diagonal, `eigen`, `eigenGeneral`, `svd`, `qr`, `lu`, `schur`, `pinv`, `det`, `logdet`, `trace`, `expm`, `rank`, `eye`, `norm`, `cond`, `sqrtm`, `logm`, `order` and `sort`, `cumsum` and `cumprod`, `cbind`, `rbind`, `~` and `~=`, copying `m[a..b, c..d]`, the reference views `sub`/`rows`/`cols`, strides, views of views, lifetime of the data and bounds errors (dgretl/matrix.d).
- `test_dgretl_estimate.d`: `dgretl.estimate`: `ols` (classical, HC0-HC3, HAC, cluster-robust), `wls`, `tsls`, `logit`, `probit` (on matrices and on a `Dataset` by series name) checked against by-hand formulas with the matrix operators, option restoration, errors, collinearity, and a heap check.
- `test_urcdist.d`: MacKinnon and DF-GLS p-values (plugin/urcdist.d): the compiled-in tables against the .bin files bit for bit, reference values from the original C code, published critical values, DF-GLS surfaces, plugin look-up.
- `test_stats_tables.d`: plugin/stats_tables.d: the Durbin-Watson table look-up `dw_lookup` (the compiled-in table against dwdata.gz, all entries against an independent parse, published values, nearest n and k, errors, `gretl_get_DW`) `rank_sum_lookup` (the printed critical values of the rank-sum test, and all 78 table rows against the symmetry lower + upper = nA (nA + nB + 1)) `stock_yogo_lookup` (the published Stock and Yogo values, the out-of-range cases, the structure of all 199 rows, `get_stock_yogo_critvals`) `qlr_asy_pvalue` (the chi-square case, the lambda conversion, Andrews' critical values, ordering, continuity and interpolation) `qlr_critval_15_05` (the bisection's stopping interval, growth with df, Andrews' values) and `get_IPS_critvals` (grid values, the three interpolation rules, beyond the table, errors, ordering and monotonicity of all six tables, the corrected N = 7, T = 50 entry) `IPS_tbar_moments` (grid values, inverse-distance interpolation, T >= 1000, the error for T < 6, monotone variance) and `IPS_tbar_rho_moments` (all four tables, interpolation, the zero cells as errors, agreement with Table 1 at p = 0, structure)..
- `test_trace_pvalue.d`: Johansen trace-test p-values (plugin/johansen.d `trace_pvalue`): tables against the C source text, 80 values against an independent computation, published critical values, plugin look-up.
- `test_johansen.d`: the Johansen cointegration test and VECM estimation (plugin/johansen.d) in the five deterministic cases against an independent computation with the matrix operators; VECM alpha, beta, standard errors and equations; a simple restriction on beta (LR test and restricted estimation); bootstrap impulse responses.
- `test_jalpha.d`: restrictions on alpha in a VECM (plugin/jalpha.d): the LR test and the restricted estimation against an independent computation of the conditional moment matrices and eigenproblem.
- `test_jrestrict.d`: general restrictions on beta and alpha in a VECM (plugin/jrestrict.d): identification, switching algorithm, scale removal, alpha restrictions, r = 2, estimation from scratch and bootstrap; all against independent computations.
- `test_lbfgsb.d`: the L-BFGS-B routines (gretl/minpack/lbfgsb.d): five optimization problems identical bit for bit to the original C code, known solutions, an error message, and `BFGS_max` with the LBFGS option.
- `test_kalman.d`: the Kalman filter, smoothers and simulator (gretl/kalman.d) through the Kalman bundle interface: filter outputs and log-likelihood against an independent filter (one and two states, correlated observables, sequential and de Jong code paths, Q/R form), exact diffuse initialization, state and disturbance smoothing against the RTS and Durbin-Koopman recursions, simulation, copy, serialization, errors. The state and disturbance smoothers are run for all three filter code paths (the standard-filter ones come from kalman_legacy.c).
- `test_kalman_arma.d`: exact ML for ARMA with the Kalman filter (plugin/kalman_arma.d): ARMA, MA, AR, seasonal and ARIMA-by-levels models (with missing observations) on simulated data against an independent exact likelihood; the needed functions of arma.c are defined in the test.
- `test_arma.d`: the ARMA/ARIMA driver and initialization (plugin/arma.d, arma_init.d): arma_model with exact ML by the Kalman filter (own and user initial values; ARMA, ARIMA, seasonal, ARMAX), conditional AR by OLS, exact ML by AS 197 (default) and AS 154 (also with missing observations) against the Kalman route, OPG covariance, conditional ML by BHHH (ARMA, MA, seasonal MA), order selection, against independent computations; error returns, plugin table. The seasonal and ARMAX default initialization (NLS through genr and MINPACK) is covered.
- `test_garch.d`: the GARCH driver (plugin/garch.d): list errors, scaling, model statistics, parameter names, LR test, standardized residuals, ARMA-based initialization, verbose output, plugin table.
- `test_fcp.d`: GARCH by the FCP algorithm (plugin/fcp.d, option --fcp): GARCH(1,1), ARCH(2), GARCH(1,2) against an independent likelihood, score, covariance matrix.
- `test_mod_garch.d`: GARCH by BFGS / Newton (plugin/mod_garch.d, the default estimator): GARCH(1,1), ARCH(2), GARCH(1,2) against an independent likelihood, score, covariance matrix and its variants, agreement with FCP.
- `test_gretl_func.d`: user-defined functions and function packages (gretl_func0..4.d): indentation/normalization of hansl code, function definition and introspection, call structures, packages (properties, writing and reloading a .gfn); body execution waits for the interpreter.
- `test_genlex.d`: the genr lexer (genlex.d): token streams, symbol tables and lookups, function pointers, getsymb, function overrides.
- `test_genr.d`: genr end to end through generate() (genlex, gensyntax, geneval0-5, genmain): scalar, series, matrix and string expressions against independent computations, and error cases.
- `test_minpack.d`: MINPACK (enorm, qrfac, qrsolv, lmpar, fdjac2, lmder, lmdif) against normal equations, analytic Jacobians, OLS and known solutions.

- `test_dgretl_timeseries.d`: dgretl/timeseries.d (correlogram, Ljung-Box, ARCH test, ADF/KPSS/Engle-Granger, filters, lag/diff, AR/ARCH models, VAR with irf/fevd/forecast/lag selection, ARMA/ARIMA, GARCH, Johansen test and VECM (with alpha and beta restrictions), ARMA and GARCH forecasts, forecasts over a range, periodogram and fractional integration, the Dataset forms) against independent computations.
- `test_dgretl_descriptive.d`: dgretl/descriptive.d (statistics, quantiles, correlations, summary, frequency tables, normality tests, two-sample tests, distances, Mahalanobis, cross-tabulation, the sign, signed-rank and rank-sum tests) against direct computations.
- `test_dgretl_random.d`: dgretl/random.d (the distributions against theory, the source hook of gretl/random.d, prng.d's MRG31k3p streams against an independent implementation and the 2^134 jump, streams driving the library's generators, parallelDraws equal to the sequential run, no leaks).
- `test_dgretl_distributions.d`: dgretl/distributions.d (densities, distribution functions and quantiles of all 15 distributions against closed forms, std.mathspecial, direct sums, numerical integration and the defining property of a quantile; non-central t, chi-square and F; both tails and their far ends; shapes, the Dataset forms, NaN, edges, parameter errors, leaks).
