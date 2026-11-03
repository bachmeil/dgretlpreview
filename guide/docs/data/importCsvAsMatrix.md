# importCsvAsMatrix

## Signature

```d
GretlMatrix importCsvAsMatrix(string fname);
```

## Parameters

- `fname`: the name of the comma-separated-values file.

## Return value

A new `GretlMatrix` with the numbers of the file, every column included: unlike `importCsv`, a first column of years or labels is a matrix column too. A failure (a missing file, text that is not a number) throws a `GretlException` with the library's error text.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    import std.file: write, remove;

    write("scores.csv", "year,math,reading\n2021,81,90\n2022,72,65\n2023,95,88\n2024,60,71\n");     // a small file to read

    auto m = importCsvAsMatrix("scores.csv");
    writefln("%d x %d", m.nrow, m.ncol);
    writeln(m[0, 0]);                                 // the year column is a column too
    writeln(Col(m, 1).dup);                           // the math column

    try {
        importCsvAsMatrix("no-such-file.csv");
    } catch (GretlException e) {
        writeln("error: ", e.msg);
    }
    remove("scores.csv");
}
```

Output:

```
4 x 3
2021
[81, 72, 95, 60]
error: importCsvAsMatrix("no-such-file.csv") failed: gzopen: No such file or directory
```
