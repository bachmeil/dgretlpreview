# importSpreadsheet

## Signature

```d
Dataset importSpreadsheet(string fname, int sheet = 1, int columnOffset = 0, int rowOffset = 0, bool verbose = false);
```

## Parameters

- `fname`: the name of the file. Its extension chooses the format: `.xlsx` or `.xls` (Excel), `.ods` (OpenDocument) or `.gnumeric` (Gnumeric); any other extension throws a `GretlException`.
- `sheet`: the worksheet to read, counting from 1 (the first). The default is 1.
- `columnOffset`: the number of columns to skip at the left, counting from 0, before the header row. The default is 0.
- `rowOffset`: the number of rows to skip at the top, counting from 0, so that the next row is the header row with the series names (as in the gretl `open` command). The default is 0.
- `verbose`: print the importer's messages. The default is `false`.

## Return value

A new `Dataset` with one series for each column of the worksheet after the offsets and one observation for each row below the header. The importer's progress messages are collected and shown only with `verbose`; if the import fails (a missing file, no such sheet, an unreadable workbook) a `GretlException` is thrown whose message includes the library's error text and those messages. Column headers that are not valid series names are changed by the library: a heading such as `2005:Q2` becomes `Q2`, and a repeated name gets a suffix. A worksheet of the kind the BEA publishes, with title rows above the header row, is read by skipping those rows with `rowOffset`.

If the worksheet has a series per row rather than per column, as the state GDP table does, the dataset has to be turned around by hand; the chapter [Loading a spreadsheet and plotting a series](../../loading-and-plotting.md) works through that example.

## Example

```d
import std.stdio;
import dgretl;

void main() {
    auto raw = importSpreadsheet("gdp-by-state.xlsx", rowOffset: 5);     // skip the five title rows
    writefln("%d observations, %d series", raw.nobs, raw.nseries);
    writeln(raw.names[0 .. 4]);                       // GeoFIPS, GeoName, and the first two quarters
    writeln(raw.names[$ - 1]);

    try {
        importSpreadsheet("data.csv");                // not a spreadsheet extension
    } catch (GretlException e) {
        writeln("error: ", e.msg);
    }
}
```

Output:

```
63 observations, 87 series
["GeoFIPS", "GeoName", "Q2", "Q3"]
Q2_aau
error: importSpreadsheet("data.csv"): the file name must end in .xlsx, .xls, .ods or .gnumeric
```
