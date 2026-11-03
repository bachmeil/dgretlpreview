# importCsv

## Signature

```d
Dataset importCsv(string fname, string[] options...);
Dataset importCsv(Options...)(string fname);   // the options as compile-time strings: importCsv!("noHeader")("x.csv")
```

## Parameters

- `fname`: the name of the comma-separated-values file.
- `options`: any of these strings:
    - `"noHeader"`: the file has no header row, so the series are named `v1`, `v2`, ... instead of from the first row.
    - `"allColumns"`: keep a first column of labels or dates as a series. By default such a column becomes the observation index (the observation labels), not a series.
    - `"codedStrings"`: turn columns of text into coded numeric values (1, 2, ... for each distinct string) instead of failing.
    - `"verbose"`: print the importer's messages.
  An unknown option is a failed `assert` in the run-time form and a compile-time error in the compile-time form.

## Return value

A new `Dataset` with one series for each column of the file and one observation for each row after the header. The first row holds the series names (unless `"noHeader"` is given), and a first column of labels, dates or years (header empty, `obs`, `date`, `year`) becomes the observation index instead of a series (unless `"allColumns"`). A column of text becomes a series of coded values, numbered 1, 2, ... in order of first appearance. A value that is not a number in a numeric column (such as `n/a`), an empty field and `NA` are read as missing values (NaN). The importer's progress messages go to an internal buffer and are shown with `"verbose"`; if the import fails (a missing file, an unreadable format, text in a column without `"codedStrings"`) a `GretlException` is thrown whose message includes the library's error text and those messages.

The options can also be given at compile time, `importCsv!("noHeader", "allColumns")("x.csv")`: the compiler then works out the flags and nothing is looked up when the call runs, which matters when the call is repeated many times. `importCsvAsMatrix` reads the numbers into a matrix instead.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    import std.file: write, remove;

    write("scores.csv", "obs,math,reading\nann,81,90\nbob,72,65\ncy,95,88\ndee,60,71\n");     // a small file to read

    auto d = importCsv("scores.csv");                 // the first column, obs, holds labels: not a series
    writeln(d.names);

    auto all = importCsv("scores.csv", "allColumns"); // ... unless asked to keep it (renamed obs_)
    writeln(all.names);
    writeln(all["obs_"]);                             // text is coded: 1, 2, ... in order of first appearance

    write("grades.csv", "math,grade
    81,A
    72,B
    95,A
    60,C
    ");
    auto g = importCsv("grades.csv", "codedStrings");  // a column of text becomes coded values
    writeln(g["grade"]);

    write("raw.csv", "1.5,2.5\n3.5,4.5\n5.5,6.5\n");
    auto raw = importCsv("raw.csv", "noHeader");       // no header row: the series are v1, v2
    writeln(raw.names);
    auto fast = importCsv!("noHeader")("raw.csv");      // the same, with the option fixed at compile time
    writeln(fast["v2"]);

    try {
        importCsv("no-such-file.csv");
    } catch (GretlException e) {
        writeln("error: ", e.msg);
    }
    remove("scores.csv");
    remove("grades.csv");
    remove("raw.csv");
}
```

Output:

```
["math", "reading"]
["obs_", "math", "reading"]
[1, 2, 3, 4]
[1, 2, 1, 3]
["v1", "v2"]
[2.5, 4.5, 6.5]
error: importCsv("no-such-file.csv") failed: gzopen: No such file or directory
Couldn't open no-such-file.csv
```
