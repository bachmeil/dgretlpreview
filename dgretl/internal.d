/* Helpers shared by the front-end modules (estimate, timeseries): temporary datasets and matrices made from D arrays, library settings, error
   texts. Internal: dgretl.package does not import this module. */

module dgretl.internal;

import core.stdc.string: strcpy;
import std.conv: to;
import std.exception: enforce;
import std.math: isFinite;
import std.string: fromStringz, toStringz;

import dgretl.dataset: Dataset;
import dgretl.exception;
import dgretl.matrix;
import gretl.gretl_errors: gretl_errmsg_get;
import gretl.gretl_list: gretl_list_new;
import gretl.libgretlh: Gretlopt, Model;
import gretl.gretl_prn: Prn, gretl_print_destroy, gretl_print_get_buffer, gretl_print_new;
import gretl.gretl_prnh: PrnType;
import gretl.libset: execute_set;

GretlMatrix fresh(size_t r, size_t c) @trusted {
    auto raw = gretl_matrix_alloc(cast(int) r, cast(int) c);

    enforce!GretlException(raw !is null || r * c == 0, "out of memory");
    return (raw is null) ? GretlMatrix.init : GretlMatrix(raw);
}

void requireFinite(GretlMatrix m, string what) @trusted {
    foreach (j; 0 .. m.ncol) {
        foreach (i; 0 .. m.nrow) {
            enforce!GretlException(isFinite(m[i, j]), what ~ ": the value in row " ~ to!string(i) ~ ", column " ~ to!string(j) ~ " is missing or not finite");
        }
    }
}

void setLibrary(string name, string value) @trusted {
    const err = execute_set(name.toStringz, value.toStringz, null, Gretlopt.OPT_NONE, null);

    enforce!GretlException(err == 0, "could not set " ~ name ~ " = " ~ value);
}

string errorText(const Model m) @trusted {
    auto msg = fromStringz(gretl_errmsg_get()).idup;

    return "error code " ~ to!string(m.errcode) ~ (msg.length ? ": " ~ msg : "");
}

/* the first `n` entries of an array of the library, as an n x 1 matrix */
GretlMatrix columnOf(const(double) * a, size_t first, size_t n) @trusted {
    auto m = fresh(n, 1);

    foreach (i; 0 .. n) m.ptr.val[i] = a[first + i];
    return m;
}

/* the temporary dataset: series 0 is the constant, then one series for each column in `cols`, called names[1], names[2], ... */
Dataset makeDataset(const(double)[][] cols, string[] names) @trusted {
    const n = cols[0].length;
    auto ds = Dataset(1 + cols.length, n);
    auto raw = ds.ptr;

    foreach (v, nm; names) strcpy(raw.varname[v], (nm ~ "\0").ptr);
    foreach (j, col; cols) {
        foreach (i; 0 .. n) raw.Z[1 + j][i] = col[i];
    }
    return ds;
}

/* the columns of a (contiguous) GretlMatrix as slices of its storage, for makeDataset */
const(double)[][] columnsOf(GretlMatrix m) @trusted {
    const(double)[][] cols;

    foreach (j; 0 .. m.ncol) cols ~= m.ptr.val[j * m.nrow .. (j + 1) * m.nrow];
    return cols;
}

string[] numbered(string stem, size_t count) {
    string[] r;

    foreach (j; 0 .. count) r ~= stem ~ to!string(j + 1);
    return r;
}

int * makeList(int[] v) @trusted {
    auto l = gretl_list_new(cast(int) v.length);

    enforce!GretlException(l !is null, "out of memory");
    foreach (i, x; v) l[i + 1] = x;
    return l;
}

/* the chosen series of a Dataset, with the observations that have a missing (NaN) value in any of them left out (listwise deletion),
   so that datasets read from files with gaps can be used directly; `kept` says how many observations are left */
struct Selection {
    double[][] cols;            // one array per series, all of the same length
    size_t kept;
}

Selection selectSeries(Dataset d, string[][] groups) @trusted {
    const(double)[][] series;

    foreach (g; groups) {
        foreach (nm; g) series ~= d[nm];
    }
    Selection s;
    const n = d.nobs;
    bool[] keep = new bool[n];

    foreach (i; 0 .. n) {
        keep[i] = true;
        foreach (col; series) {
            if (!isFinite(col[i])) keep[i] = false;
        }
    }
    foreach (col; series) {
        double[] c;

        foreach (i; 0 .. n) {
            if (keep[i]) c ~= col[i];
        }
        s.cols ~= c;
    }
    s.kept = s.cols.length ? s.cols[0].length : 0;
    return s;
}

GretlMatrix matrixFrom(double[][] cols) @trusted {
    const n = cols.length ? cols[0].length : 0;
    auto m = fresh(n, cols.length);

    foreach (j, c; cols) {
        foreach (i; 0 .. n) m.ptr.val[j * n + i] = c[i];
    }
    return m;
}


/* the text the library prints into a buffer printer */
struct Printer {
    Prn * prn;

    this(int dummy) @trusted {
        int err;

        prn = gretl_print_new(PrnType.buffer, &err);
        enforce!GretlException(prn !is null && err == 0, "could not create a print buffer");
    }

    ~this() @trusted {
        if (prn !is null) gretl_print_destroy(prn);
    }

    @disable this(this);

    string text() @trusted {
        return fromStringz(gretl_print_get_buffer(prn)).idup;
    }
}

