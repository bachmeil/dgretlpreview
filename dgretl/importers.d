/* Functions that read data files into the dgretl types: importCsv, importCsvAsMatrix and importSpreadsheet. The
   importer's progress messages go to an internal buffer; they are shown with the "verbose" option and are part of
   the exception message (GretlException) when an import fails. */

module dgretl.importers;

import std.exception: enforce;
import std.path: extension;
import std.stdio: write;
import std.string: fromStringz, toStringz;
import std.uni: toLower;

import gretl.csvdata: import_csv, import_csv_as_matrix;
import gretl.dataio: import_spreadsheet;
import gretl.dataioh: GretlFileType;
import gretl.dataset: datainfo_new;
import gretl.gretl_errors: gretl_errmsg_get;
import gretl.gretl_matrix: gretl_matrix_free;
import gretl.gretl_matrixh: GretlMatrixRaw;
import gretl.gretl_prn: Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh: PrnType;
import gretl.libgretlh: DatasetRaw, Gretlopt;

import dgretl.dataset: Dataset;
import dgretl.exception: GretlException;
import dgretl.matrix: GretlMatrix;

private string errorText() {
    const(char) * s = gretl_errmsg_get();
    return (s is null) ? "" : fromStringz(s).idup;
}

/* The options of importCsv, given as strings: "codedStrings", "noHeader", "allColumns" and "verbose". One function
   does the lookup, so that it serves both the run-time call importCsv("x.csv", "noHeader") and the compile-time call
   importCsv!("noHeader")("x.csv"); it allocates nothing and can be evaluated by the compiler. */
private bool csvOption(string name, ref Gretlopt opt, ref bool verbose) pure nothrow @nogc @safe {
    switch (name) {
        case "codedStrings": opt |= Gretlopt.OPT_N; return true;
        case "noHeader": opt |= Gretlopt.OPT_H; return true;
        case "allColumns": opt |= Gretlopt.OPT_A; return true;
        case "verbose": verbose = true; return true;
        default: return false;
    }
}

/* the options as the compiler sees them: the Gretlopt flags and the verbose switch of a list of strings, with a compile
   time error for a name that is not an option */
private struct CsvChoice {
    Gretlopt opt;
    bool verbose;
}

private template csvChoice(Options...) {
    private CsvChoice compute() {
        CsvChoice c;
        static foreach (name; Options) {
            static assert(is(typeof(name) == string), "importCsv: the options must be strings, not " ~ typeof(name).stringof);
            static assert(csvOptionKnown(name), "\"" ~ name ~ "\" is not an option for function importCsv.");
            csvOption(name, c.opt, c.verbose);
        }
        return c;
    }

    enum csvChoice = compute();
}

private bool csvOptionKnown(string name) pure nothrow @nogc @safe {
    Gretlopt o;
    bool v;
    return csvOption(name, o, v);
}

/* Read a comma-separated-values file into a new Dataset. The first row holds the series names (unless "noHeader" is
   given) and a first column of labels or dates becomes the observation index, not a series (unless "allColumns" is
   given). With "codedStrings", columns of text are turned into coded values instead of being an error. "verbose" shows
   the importer's messages. Example: importCsv("x.csv", "noHeader", "allColumns").

   The options can also be given at compile time, importCsv!("noHeader", "allColumns")("x.csv"): the flags are then
   worked out by the compiler (a name that is not an option is a compile-time error) and nothing is built or looked up
   when the call runs, which matters when the call is repeated many times. */

Dataset importCsv(string fname, string[] options...) @trusted {
    Gretlopt opt = Gretlopt.OPT_NONE;
    bool verbose = false;

    foreach (option; options) {
        assert(csvOption(option, opt, verbose), option ~ " is not an option for function importCsv.");
    }

    return importCsvWith(fname, opt, verbose);
}

/// ditto
Dataset importCsv(Options...)(string fname) @trusted if (Options.length > 0) {
    enum choice = csvChoice!Options;

    return importCsvWith(fname, choice.opt, choice.verbose);
}

private Dataset importCsvWith(string fname, Gretlopt opt, bool verbose) @trusted {
    DatasetRaw * raw = datainfo_new();
    enforce!GretlException(raw !is null, "importCsv: out of memory");

    auto ds = Dataset(raw);    /* owns the raw dataset from here on, also when an exception is thrown */
    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    enforce!GretlException(prn !is null, "importCsv: could not create a message buffer");
    scope (exit) gretl_print_destroy(prn);

    err = import_csv(fname.toStringz, raw, opt, prn);

    const(char) * buf = gretl_print_get_buffer(prn);
    string messages = (buf is null) ? "" : fromStringz(buf).idup;

    if (err != 0) {
        string text = errorText();
        throw new GretlException("importCsv(\"" ~ fname ~ "\") failed" ~ (text.length ? ": " ~ text : "") ~ (messages.length ? "\n" ~ messages : ""));
    }

    if (verbose) {
        write(messages);
    }

    return ds;
}

/* Read the numbers of a comma-separated-values file into a matrix, every column included (a first column of years or
   labels is a matrix column too). */

GretlMatrix importCsvAsMatrix(string fname) @trusted {
    int err = 0;
    GretlMatrixRaw * m = import_csv_as_matrix(fname.toStringz, &err);

    if (err != 0 || m is null) {
        if (m !is null) {
            gretl_matrix_free(m);
        }

        string text = errorText();
        throw new GretlException("importCsvAsMatrix(\"" ~ fname ~ "\") failed" ~ (text.length ? ": " ~ text : ""));
    }

    return GretlMatrix(m);
}

/* Read one worksheet of an Excel (.xlsx, .xls), OpenDocument (.ods) or Gnumeric file into a new Dataset. The sheet is
   chosen by number (1 is the first) and columnOffset / rowOffset skip that many columns and rows before the header
   row (the offsets are 0-based counts, as in the gretl "open" command). */

Dataset importSpreadsheet(string fname, int sheet = 1, int columnOffset = 0, int rowOffset = 0, bool verbose = false) @trusted {
    GretlFileType ftype;
    string ext = extension(fname).toLower;

    if (ext == ".xlsx") ftype = GretlFileType.xlsx;
    else if (ext == ".xls") ftype = GretlFileType.xls;
    else if (ext == ".ods") ftype = GretlFileType.ods;
    else if (ext == ".gnumeric") ftype = GretlFileType.gnumeric;
    else throw new GretlException("importSpreadsheet(\"" ~ fname ~ "\"): the file name must end in .xlsx, .xls, .ods or .gnumeric");

    DatasetRaw * raw = datainfo_new();
    enforce!GretlException(raw !is null, "importSpreadsheet: out of memory");

    auto ds = Dataset(raw);
    int err = 0;
    Prn * prn = gretl_print_new(PrnType.buffer, &err);

    enforce!GretlException(prn !is null, "importSpreadsheet: could not create a message buffer");
    scope (exit) gretl_print_destroy(prn);

    int[4] params = [3, sheet, columnOffset, rowOffset];

    err = import_spreadsheet(fname.toStringz, ftype, params.ptr, null, raw, Gretlopt.OPT_NONE, prn);

    const(char) * buf = gretl_print_get_buffer(prn);
    string messages = (buf is null) ? "" : fromStringz(buf).idup;

    if (err != 0) {
        string text = errorText();
        throw new GretlException("importSpreadsheet(\"" ~ fname ~ "\") failed" ~ (text.length ? ": " ~ text : "") ~ (messages.length ? "\n" ~ messages : ""));
    }

    if (verbose) {
        write(messages);
    }

    return ds;
}
