# Dataset

## Signature

```d
struct Dataset {
    this(T1, T2)(T1 nvar, T2 nobs, bool markers = false);
    this(DatasetRaw* raw);
    size_t nobs();
    size_t nseries();
    string[] names();
    const(double)[] opIndex(size_t i);
    const(double)[] opIndex(string name);
    DatasetRaw* ptr();
}
```

## Parameters

- `nobs()`: the number of observations (rows).
- `nseries()`: the number of series (columns), not counting the constant that the library keeps as series 0.
- `names()`: the names of the series, in order, as a `string[]`.
- `d[i]`: the values of series number `i`, counting from 0, as a `const(double)[]`. An index outside the range throws a `GretlException`.
- `d["name"]`: the values of the series with that name, as a `const(double)[]`. A name that does not exist throws a `GretlException`.
- `this(nvar, nobs, markers)`: a new dataset with `nvar` series (counting the constant, as the library's `create_new_dataset` does) and `nobs` observations; `markers` asks for observation labels. This is for programs that build a dataset with the low-level functions.
- `this(DatasetRaw* raw)`: take over a raw dataset made with the low-level functions. From then on the `Dataset` owns it and frees it with `destroy_dataset`, so the raw pointer must not be freed or used after the last copy has gone.
- `ptr()`: the raw pointer (`DatasetRaw*`), to call the low-level functions of the port with; null if there is no data.

## Return value

A `Dataset` is a set of series of equal length, as read from a data file by `importCsv` or `importSpreadsheet`. It is reference counted: copies share the same data, which is freed when the last copy goes away, so nothing has to be freed by hand. The series returned by `d[i]` and `d["name"]` are views into that data (read-only) and are valid as long as the `Dataset` or a copy of it lives; use `.dup` to keep a copy that outlives it. The series pass straight into the statistics, tests and estimators of dgretl (for example `mean(d["math"])` or `ols(d["y"], d["x"])`), and many of those functions also take a `Dataset` and the names of series directly.

A program that has a type of its own called `Dataset` should use a selective import or the qualified name `dgretl.dataset.Dataset`.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    import std.file: write, remove;

    write("scores.csv", "obs,math,reading\nann,81,90\nbob,72,65\ncy,95,88\ndee,60,71\n");     // a small file to read

    auto d = importCsv("scores.csv");                 // the first column, obs, holds labels: it is not a series
    writefln("%d observations, %d series: %s", d.nobs, d.nseries, d.names);

    auto math = d["math"];                            // a series by name
    writeln(math);
    writeln(d[1]);                                    // a series by number: the second one, reading
    writefln("mean of math %.2f", mean(d["math"]));

    auto own = d["reading"].dup;                      // a copy that outlives the Dataset
    own[0] = 0;
    writeln(d["reading"][0]);                         // the data are unchanged

    try {
        auto missing = d["history"];
    } catch (GretlException e) {
        writeln("error: ", e.msg);
    }
    remove("scores.csv");
}
```

Output:

```
4 observations, 2 series: ["math", "reading"]
[81, 72, 95, 60]
[90, 65, 88, 71]
mean of math 77.00
90
error: Dataset: there is no series called "history"
```
