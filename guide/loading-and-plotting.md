# Loading a spreadsheet and plotting a series

This chapter shows how to use the D port of libgretl from your own program: read the state GDP spreadsheet that
ships in this directory (`gdp-by-state.xlsx`), pick out one state, and plot it. Everything below was run with
`ldmd2` (the LDC front end that accepts `dmd`'s command line), and `kansas.d` in this directory is the complete,
working program.

![Kansas real GDP growth, annualized quarterly rates, 2005Q2 to 2026Q2](kansas.png)

## The data

`gdp-by-state.xlsx` is the BEA table SQGDP1: real GDP, percent change from the preceding period, annualized.

| Rows 1-4 | Row 5 | Row 6 | Rows 7-66 | Rows 67-69 |
|---|---|---|---|---|
| title, units, source, "State or DC" | blank | header: `GeoFIPS`, `GeoName`, `2005:Q2` ... `2026:Q2` | the United States, the 50 states and DC, and 8 BEA regions | footnotes |

Columns A and B identify the geography, and columns C to CI hold the 85 quarters from 2005Q2 to 2026Q2. So the
spreadsheet is the transpose of what gretl normally expects (one series per column, one observation per row):
here each *row* is a series and each *column* is a date. The program therefore imports the sheet as it is and
then builds a proper quarterly time series from the Kansas row.

## What you need

- `ldmd2` (LDC 1.39 was used)
- the development packages for libxml2, LAPACK, BLAS and FFTW3
- `gnuplot` on your `PATH` (only for drawing the plot)

## Building

```bash
guide/build.sh kansas
```

This script (`guide/build.sh`) does two things, and you can reproduce them by hand if you prefer:

1. **Compiles `testing/stubs_unported.d`** on its own (`ldmd2 -c -of=build/stubs_unported.o ../testing/stubs_unported.d`).
2. **Runs `ldmd2`** on your program together with every `.d` file of the port (`gretl/` holds the direct ports of the libgretl files, `gretl/plugin/` the plugins, `cephes/` the
probability and special-function library, and `dgretl/` the D-friendly front end that `import dgretl;` brings in) and
`testing/stubs.d`:

   ```bash
   ldmd2 -I.. -of=kansas kansas.d ../gretl/*.d ../gretl/plugin/*.d ../gretl/minpack/*.d ../cephes/*.d ../dgretl/*.d \
         ../testing/stubs.d \
         build/stubs_unported.o \
         -L-lm -L-lxml2 -L-lz -L-llapack -L-lblas -L-lfftw3
   ```

   Parts of libgretl have not been ported yet. `testing/stubs_unported.d` gives every function the ported code
   could call in them a stand-in that prints `*** call to unported libgretl function NAME()` and aborts, so a
   missing piece fails loudly instead of returning garbage (`testing/make_unported_stubs.sh` regenerates the file
   when more of the library is ported). `testing/stubs.d` holds working stand-ins for the unported pieces that the library calls in normal use, for
   example "no user functions are defined" or "single-threaded". Until those modules are ported, any program
   that uses the library links it in, not only the test suite.

There is no separate library to install: the port is compiled together with your program, which takes a few
seconds.

## The program, step by step

The complete program is `kansas.d`. It starts with a single import:

```d
import dgretl;
```

### Starting the library

There is nothing to do. `dgretl` starts libgretl before `main` runs (settings, random number generator, XML
library, a per-process scratch directory, the gnuplot program and a default PNG font for plots) and shuts it down,
removing the scratch directory, after `main` returns. A program that does not import `dgretl` has to call
`libgretl_init()` and `libgretl_cleanup()` itself and set those paths by hand.

### Read the worksheet

```d
auto raw = importSpreadsheet("gdp-by-state.xlsx", rowOffset: 5);
```

`importSpreadsheet` reads one worksheet of an .xlsx, .xls, .ods or Gnumeric file into a new `Dataset`. The sheet
is chosen by number (`sheet:`, the first by default) and `columnOffset:` and `rowOffset:` skip that many columns and
rows before the header row. A row offset of 5 skips the title block, so row 6 is read as the header. The importer's
progress messages are collected and only shown with `verbose: true`; if the import fails, `importSpreadsheet` throws a
`GretlException` whose message contains them.

A `Dataset` is a reference-counted struct: copies share the data, and the data are freed when the last copy goes
away. `raw` now has:

| | |
|---|---|
| series 0 `GeoFIPS` | the FIPS codes, as strings |
| series 1 `GeoName` | the names, as strings |
| series 2 to 86 | the 85 quarters, as numbers |
| `raw.nobs` = 63 | the 60 geographies of rows 7-66, plus the 3 footnote rows, which have no numbers |

The series of a `Dataset` are numbered from 0; the constant that libgretl keeps as series 0 internally is not
visible. The header labels `2005:Q2`, `2005:Q3` and so on are not valid names, so the importer reduces each to
`Q2`, `Q3`, ... and disambiguates repeats (`Q2_aau`), with a warning. Never look a quarter up by name: **use the
position**, as the program does (`raw[2 + q]`).

### Find Kansas

```d
foreach (t; 0 .. raw.nobs) {
    const(char) * name = series_get_string_for_obs(raw.ptr, 2, cast(int) t);
    if (name !is null && strcmp(name, "Kansas") == 0) { row = t; break; }
}
```

String-valued series keep their text in a table, which the wrapper does not expose yet; `raw.ptr` is the underlying
`DatasetRaw *` (the C-layout struct all the lower-level functions take), and
`series_get_string_for_obs(dset, series, obs)` returns the string for one observation. At this low level the
constant counts as series 0, so `GeoName` is series 2. Observations are 0-based.

### Build a quarterly time series

```d
auto ks = Dataset(2, nquarters);                                 // the constant plus one series, 85 observations
dataset_set_time_series(ks.ptr, 4, 2005, 2);                     // quarterly, starting 2005Q2
foreach (q; 0 .. nquarters)
    ks.ptr.Z[1][q] = raw[2 + q][row];
```

This is the transposition. `Dataset(nvars, nobs)` makes a new dataset (the constant counts in `nvars`, as in the
low-level `create_new_dataset`); `Dataset(raw)` instead takes over a raw dataset you made yourself. Either way the
`Dataset` frees it when it goes away.
`dataset_set_time_series(dset, pd, year, minor)` declares the frequency (4 for quarterly) and the start; it is
this that gives the plot a date axis. (The series name is set directly in `ks.ptr.varname[1]`.)

### Plot

```d
set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "kansas.png");   // --output=kansas.png
int[2] list = [1, 1];                                                      // one series: number 1
int err = gnuplot(list.ptr, null, ks.ptr, cast(Gretlopt) (Gretlopt.OPT_T | Gretlopt.OPT_O));
```

- `gnuplot(list, literal, dset, opt)` (in `graphing.d`, low level: it takes the raw dataset) is the `gnuplot`
  command. `list` holds the series to plot, with the constant counted as series 0, so the Kansas series is number
  1. `literal` can hold extra gnuplot commands (`null` for none).
- The options are gretl's command options, as `Gretlopt` flags. `OPT_T` is `--time-series`: with a single
  series in the list it is plotted against time. `OPT_O` is `--with-lines`. Without `OPT_T`, a one-series list
  is an error, because gretl then expects `y` and `x`.
- The output file is chosen the way `--output=` would on the command line, through an option value attached
  to the command. The file type follows from the extension (`.png`, `.pdf`, `.eps`, `.fig` and `.emf` are recognized). The special
  value `"display"` opens an interactive gnuplot window instead.

The program writes `kansas.png` into the current directory. There is no cleanup to do: both datasets are freed when
`main` returns and the library shuts itself down.

## Running it

```bash
cd guide
./kansas
```

```
Kansas: 85 quarters, first -0.2, last -0.6
wrote kansas.png
```

(A warning about duplicated variable names, described above, is printed first.)

The first and last values match cells C24 and CI24 of the spreadsheet. The series is noisy, as annualized
quarterly state growth is, and the 2020Q2 pandemic collapse (about -28%) and the sharp rebound after it stand
out.

## Things to know

- **Indexing.** A `Dataset` numbers its series from 0 and hides the constant; the low-level functions count the
  constant as series 0 (as in the plot list above). Observations are 0-based everywhere.
- **Memory.** `Dataset` and `GretlMatrix` free themselves. The raw structs (`DatasetRaw`, `GretlMatrixRaw`) behind
  them, and everything the lower-level functions return, are still allocated with `malloc` and have to be freed by
  hand unless they are handed to the `Dataset(raw)` constructor.
- **Errors.** The functions of `dgretl` throw `GretlException` with the library's message. The lower-level
  functions return 0 on success and a `GretlError` code otherwise; `gretl_errmsg_get()` has the message, if there is one.
- **Not available yet.** Anything that reaches an unported piece (the command interpreter, user functions,
  loops, MIDAS, maps, the GUI) aborts with the "call to unported libgretl function" message. Plotting from
  D, as here, does not need those.
- **`ldmd2` and `dmd`.** The same command works with `dmd`. LDC is stricter in one respect: it reports two
  `extern (C)` functions with the same name instead of silently picking one, which found and fixed a clash in
  the port.
