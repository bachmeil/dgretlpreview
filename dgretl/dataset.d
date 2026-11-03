/* Dataset: a reference-counted dataset (series of equal length, as read from a data file). Copies share the same
   data, which is freed when the last copy goes away. Series number 0 inside the library is the constant; here the
   series are numbered from 0 and the constant is not visible.

   The related modules of the port are imported publicly, so that `import dgretl.dataset;` also brings in the
   functions that work on the raw dataset (gretl.dataset, gretl.dataseth); `ptr` gives the raw pointer to call them
   with. If a program has a type of its own called Dataset, use a selective import or the qualified name
   (`dgretl.dataset.Dataset`). */

module dgretl.dataset;

public import gretl.dataset;
public import gretl.dataseth;
public import gretl.libgretlh: DatasetRaw;
public import dgretl.exception;

import std.conv: to;
import std.exception: enforce;
import std.string: fromStringz;
import std.typecons: RefCountedAutoInitialize, SafeRefCounted;

import gretl.libgretlh: Gretlopt;

struct Dataset {
    private struct Payload {
        DatasetRaw * raw;

        this(DatasetRaw * p) {
            raw = p;
        }

        ~this() {
            if (raw !is null) {
                destroy_dataset(raw);
                raw = null;
            }
        }

        @disable this(this);
    }

    private SafeRefCounted!(Payload, RefCountedAutoInitialize.no) _p;

    /* Take over a raw dataset made with the low-level functions (create_new_dataset and so on). From here on the
       Dataset owns it and frees it with destroy_dataset: the raw pointer must not be freed or used after the last
       copy has gone. */
    this(DatasetRaw * raw) {
        enforce!GretlException(raw !is null, "Dataset: null raw dataset");
        _p = SafeRefCounted!(Payload, RefCountedAutoInitialize.no)(raw);
    }

    /* A new dataset with nvar series (counting the constant, as create_new_dataset does) and nobs observations;
       markers asks for observation labels. */
    this(T1, T2)(T1 nvar, T2 nobs, bool markers = false) {
        this(create_new_dataset(nvar.to!int, nobs.to!int, markers ? 1 : 0));
    }

    /* the raw pointer, for the low-level functions of the port (null if there is no data) */
    DatasetRaw * ptr() @trusted {
        return _p.refCountedStore.isInitialized ? _p.refCountedPayload.raw : null;
    }

    /* number of observations (rows) */
    size_t nobs() @trusted {
        auto r = ptr;
        return (r is null) ? 0 : r.n;
    }

    /* number of series, not counting the constant */
    size_t nseries() @trusted {
        auto r = ptr;
        return (r is null || r.v < 1) ? 0 : r.v - 1;
    }

    /* the series names, in order */
    string[] names() @trusted {
        string[] ret;
        auto r = ptr;

        if (r !is null) {
            foreach (v; 1 .. r.v) {
                ret ~= fromStringz(r.varname[v]).idup;
            }
        }

        return ret;
    }

    /* the values of series number i (0-based); valid as long as this Dataset (or a copy) lives */
    const(double)[] opIndex(size_t i) @trusted {
        auto r = ptr;

        enforce!GretlException(r !is null && i < nseries, "Dataset: series number " ~ to!string(i) ~ " is out of range");
        return r.Z[i + 1][0 .. r.n];
    }

    /* the values of the series with the given name */
    const(double)[] opIndex(string name) @trusted {
        auto r = ptr;

        if (r !is null) {
            foreach (v; 1 .. r.v) {
                if (fromStringz(r.varname[v]) == name) {
                    return r.Z[v][0 .. r.n];
                }
            }
        }

        throw new GretlException("Dataset: there is no series called \"" ~ name ~ "\"");
    }
}
