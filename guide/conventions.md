# Conventions used throughout dgretl

The functions of dgretl follow the same rules, so what you learn from one carries over to the others. This page
collects them. The code in the examples was compiled and run, and the output shown is real.

## One import

`import dgretl;` brings in everything: the data types, the statistics, the models, the distributions and the random
numbers. If a program has a type of its own called `Dataset` or `GretlMatrix`, use a selective import
(`import dgretl: ols;`) or the qualified name (`dgretl.matrix.GretlMatrix`).

## Data go in as arrays, matrices or their parts

Wherever a function takes data, it takes an array of numbers (`double[]`, `int[]`, ...), a `GretlMatrix`, a
`Submatrix` (a view of part of a matrix), or one row or column of a matrix (`Row(m, i)`, `Col(m, j)`). A function
that needs a series accepts any of them; a function that needs a matrix accepts a `GretlMatrix` or a `Submatrix`.

```d
auto a = [1.0, 2.0, 3.0, 4.0];
auto m = cbind(a, [10.0, 20.0, 40.0, 80.0]);        // a 4 x 2 matrix: the arguments are its columns

writeln(mean(a));                                   // 2.5: an array
writeln(mean(Col(m, 1)));                           // 37.5: one column of a matrix
writeln(mean(m.cols[0 .. 1])[0, 0]);                // 2.5: part of a matrix
show(mean(m));                                      // for a matrix: one value per column
```

```
2.5
37.5
2.5
   2.500  37.500
```

(`show` here prints a matrix; it is a three-line helper, not part of dgretl.) The statistics of a matrix are those
of each column, as a 1 x k `GretlMatrix`. Many functions also have a form that takes a `Dataset` and the *names* of
series, such as `ols(d, "y", ["x1", "x2"])`; see [`Dataset`](docs/data/Dataset.md).

## Results are plain values and structs

A function returns what it computed and nothing else has to be done with it: a `double` for one number, a
`double[]` for a list of numbers, a `GretlMatrix` for a matrix, or a struct with named fields for a result that has
several parts (`Estimate`, `VarFit`, `SummaryStats`, `UnitRootTest`, ...). The fields are documented on the page of
the function that returns the struct. Nothing is freed by hand: matrices and datasets are reference counted, and
free themselves when the last copy goes away.

```d
auto fit = ols(m.cols[1 .. 2], m.cols[0 .. 1]);
writeln(typeof(fit).stringof, " with ", fit.names);
```

```
Estimate with ["const", "x1"]
```

## Copies, views and `dup`

Copying a `GretlMatrix` (assignment, passing it to a function) shares the data, so a change through one copy is seen
through the others. `m.dup` makes an independent copy. Slicing with `m[a .. b, c .. d]` makes a copy;
`m.sub[a .. b, c .. d]`, `m.rows[a .. b]` and `m.cols[a .. b]` give a `Submatrix`, a reference into `m` (see
[`GretlMatrix`](docs/matrix/GretlMatrix.md) and [`Submatrix`](docs/matrix/Submatrix.md)). The series of a `Dataset`
are read-only views; use `.dup` to keep one after the `Dataset` is gone.

```d
auto alias_ = m;                    // shares the data
alias_[0, 0] = 100;
writeln(m[0, 0]);                   // 100
auto own = m.dup;                   // an independent copy
own[0, 0] = -1;
writeln(m[0, 0]);                   // still 100
```

## Indexing starts at 0, and matrices are stored by column

Rows, columns, lags and observation numbers start at 0: `m[3, 1]` is row 3 and column 1 of the matrix, the
fourth row and the second column. The arguments of `cbind` become columns and those of `rbind` become rows. Inside
the library the elements are stored column by column; that matters only when a parameter that is a matrix is read
element by element (the distribution functions read it column by column).

## Missing values are NaN

A missing value is `double.nan`, and the statistics skip it, as the library does:

```d
writeln(mean([1.0, double.nan, 3.0]));              // 2
writeln(validCount([1.0, double.nan, 3.0]));        // 2
```

A statistic of no valid observation is NaN. The functions that cannot cope with a missing value (the time-series
tests, the models) throw a `GretlException` instead. `lag` and `diff` put NaN where nothing can be computed, and the
importers read empty fields and `NA` as NaN.

## Names, defaults and options

The function names follow R where there is an R function: `dnorm`, `pnorm`, `qnorm`, `rnorm`, `quantile`. Where an R or standard-library name would clash with D's, the name is changed:
`total` is the sum (not `std.algorithm.sum`), `validCount` the number of valid values, `rowDistances` the
library's `distance`, `expm`, `sqrtm` and `logm` the matrix functions that are not the elementwise `exp`, `sqrt`
and `log`, and `kurtosis` is the *excess* kurtosis (zero for a normal distribution).

Parameters after the required ones have defaults, and any of them can be given by name, in any order:

```d
writeln(pnorm(1.96));                               // 0.975002: the standard normal
writeln(pnorm(1.96, mean: 0, sd: 2));
writeln(pnorm(1.96, sd: 2));                        // 0.836457: the same, the mean left at its default
```

Options that are a choice among a few are enums, with the enum's name spelled out so that a call reads clearly:
`adf(x, 0, Trend.trend)`, `ols(y, x, true, Vcov.hc1)`, `norm(m, Norm.spectral)`, `quantile(x, 0.9, QuantileType.q7)`.
`bool` options are `bool` parameters with a default.

## Parameters that vary, in the distribution functions

The functions of the probability distributions take the first argument as one number, an array, a matrix or (with a
`Dataset`) the name of a series, and each parameter of the distribution as a number, an array, a matrix or a
series name, one value per element of the first argument. A number is repeated. See the pages in
[`docs/distributions/`](docs/distributions/pnorm.md).

## Random numbers

The functions of the random numbers draw from the library's generator, which is seeded from the system. Call
`setSeed(seed)` first to make a program reproducible. For simulation in parallel, `makeStreams` gives
non-overlapping streams, one for each task, and the result does not depend on the number of processors
([`makeStreams`](docs/random/makeStreams.md), [`parallelDraws`](docs/random/parallelDraws.md)).

## Time information

A series given as an array or a matrix has no dates. Where a function needs them (seasonal dummies in `adf`, the
calendar of a VAR) the optional arguments `frequency` (observations per year: 1, 4, 12, ...) and `start` (the first
observation, as `"1990:1"` or `"1990Q1"`) give them; without them the series is annual from 1. The estimation sample
can be chosen with `fitStart` and `fitEnd`, observations of that calendar. A `Dataset` with a time structure supplies
them itself.

## Errors

Every failure is a `GretlException` with a message that says what went wrong, and a function that throws leaves
nothing half done ([`GretlException`](docs/errors/GretlException.md)):

```d
try {
    rnorm(3, 0, -1);
} catch (GretlException e) {
    writeln(e.msg);                 // rnorm: need a finite mean and a positive sd
}
```
