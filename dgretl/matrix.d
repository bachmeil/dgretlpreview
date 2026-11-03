/* GretlMatrix: a reference-counted matrix of doubles. Copies share the same data, which is freed when the last copy
   goes away. Indices are 0-based; the storage underneath is column by column. `nrow` and `ncol` are the sizes.

   m[a .. b, c .. d] makes a COPY (a new GretlMatrix), so a function that changes it cannot change the original.
   A reference is asked for explicitly, with a Submatrix:
       m.sub[a .. b, c .. d]   the elements in rows a .. b-1 and columns c .. d-1
       m.rows[a .. b]          those rows, all columns
       m.cols[a .. b]          those columns, all rows
   A Submatrix reads and writes the elements of the matrix it was cut from, keeps that matrix's data alive, and has
   `sub`, `rows` and `cols` of its own. `.dup` makes a GretlMatrix copy of either type. A Submatrix has a stride (the
   distance between its columns in memory), because a block is not contiguous in the column-by-column storage.

   The related modules of the port are imported publicly (gretl.gretl_matrix, gretl.gretl_matrixh); `ptr` gives the raw
   pointer to call them with. If a program has a type of its own called GretlMatrix, use a selective import or the
   qualified name (`dgretl.matrix.GretlMatrix`). */

module dgretl.matrix;

public import gretl.gretl_matrix;
public import gretl.gretl_matrixh;
public import dgretl.exception;

import std.conv: to;
import std.exception: enforce;
import std.meta: allSatisfy;
import std.range: ElementType;
import std.traits: isArray, isNumeric;
import std.typecons: RefCountedAutoInitialize, SafeRefCounted;

struct GretlMatrix {
    private struct Payload {
        GretlMatrixRaw * raw;

        this(GretlMatrixRaw * p) {
            raw = p;
        }

        ~this() {
            if (raw !is null) {
                gretl_matrix_free(raw);
                raw = null;
            }
        }

        @disable this(this);
    }

    private SafeRefCounted!(Payload, RefCountedAutoInitialize.no) _p;

    mixin MatrixOps;

    /* Take over a raw matrix made with the low-level functions; it is freed with gretl_matrix_free when the last copy
       goes away. */
    this(GretlMatrixRaw * raw) {
        enforce!GretlException(raw !is null, "GretlMatrix: null raw matrix");
        _p = SafeRefCounted!(Payload, RefCountedAutoInitialize.no)(raw);
    }

    /* the raw pointer, for the low-level functions of the port (null if there is no matrix) */
    GretlMatrixRaw * ptr() @trusted {
        return _p.refCountedStore.isInitialized ? _p.refCountedPayload.raw : null;
    }

    size_t nrow() @trusted {
        auto r = ptr;
        return (r is null) ? 0 : r.rows;
    }

    size_t ncol() @trusted {
        auto r = ptr;
        return (r is null) ? 0 : r.cols;
    }

    /* element in row i, column j */
    double opIndex(size_t i, size_t j) @trusted {
        auto r = ptr;

        enforce!GretlException(r !is null && i < r.rows && j < r.cols, "GretlMatrix: index [" ~ to!string(i) ~ ", " ~ to!string(j) ~ "] is out of range");
        return gretl_matrix_get(r, cast(int) i, cast(int) j);
    }

    double opIndexAssign(double x, size_t i, size_t j) @trusted {
        auto r = ptr;

        enforce!GretlException(r !is null && i < r.rows && j < r.cols, "GretlMatrix: index [" ~ to!string(i) ~ ", " ~ to!string(j) ~ "] is out of range");
        gretl_matrix_set(r, cast(int) i, cast(int) j, x);
        return x;
    }

    /* m ~= x adds the columns of x to the right of m (x a matrix, a scalar or an array, as in cbind). If nothing else refers to the data of m
       (no other GretlMatrix, Submatrix, Row, Col or iterator made from it) the storage is grown in place; otherwise m is rebound to a new matrix,
       so that no other handle ever sees a change of size or a moved data pointer. */
    void opOpAssign(string op, T)(T rhs) @trusted if (op == "~" && isCbindOperand!T) {
        auto raw = ptr;

        if (raw is null || _p.refCountedStore.refCount != 1) {
            this = (raw is null) ? cbind(rhs) : cbind(this, rhs);
            return;
        }
        auto piece = makePiece(rhs);
        const nr = cast(size_t) raw.rows, oldcols = cast(size_t) raw.cols;
        const extra = piece.cols;

        cbindRows([piece], nr, true);               // checks that a matrix or an array fits the number of rows

        enforce!GretlException(oldcols + extra <= int.max, "~=: the matrix would be too large");
        if (extra == 0) return;
        enforce!GretlException(gretl_matrix_realloc(raw, cast(int) nr, cast(int) (oldcols + extra)) == 0, "~=: out of memory");
        fillColumns(raw.val + oldcols * nr, nr, piece);
    }

    /* m[a .. b, c .. d]: a copy of the block */
    GretlMatrix opIndex(MatrixRange r, MatrixRange c) @trusted {
        return wholeView(this)[r, c];
    }

    size_t opDollar(size_t dim)() {
        static if (dim == 0) return nrow; else return ncol;
    }

    MatrixRange opSlice(size_t dim)(size_t lo, size_t hi) {
        return MatrixRange(lo, hi);
    }

    GretlMatrix dup() @trusted {
        return wholeView(this).dup;
    }
}

/* a half-open range of rows or columns, from a .. b */
struct MatrixRange {
    size_t lo, hi;
}

/* A reference to a block of a matrix, with the stride between columns */
struct Submatrix {
    private SafeRefCounted!(GretlMatrix.Payload, RefCountedAutoInitialize.no) _p;   // keeps the data alive
    private double * _data;
    private size_t _nrow, _ncol, _stride;

    mixin MatrixOps;

    size_t nrow() { return _nrow; }
    size_t ncol() { return _ncol; }
    size_t stride() { return _stride; }

    /* true if the elements are one run of memory (all of them, column by column), so that `data` can be passed on */
    bool isContiguous() { return _stride == _nrow || _ncol <= 1; }

    /* pointer to element [0, 0]; element [i, j] is at data[j * stride + i] */
    double * data() @system { return _data; }

    double opIndex(size_t i, size_t j) @trusted {
        enforce!GretlException(i < _nrow && j < _ncol, "Submatrix: index [" ~ to!string(i) ~ ", " ~ to!string(j) ~ "] is out of range");
        return _data[j * _stride + i];
    }

    double opIndexAssign(double x, size_t i, size_t j) @trusted {
        enforce!GretlException(i < _nrow && j < _ncol, "Submatrix: index [" ~ to!string(i) ~ ", " ~ to!string(j) ~ "] is out of range");
        _data[j * _stride + i] = x;
        return x;
    }

    /* s[a .. b, c .. d]: a copy of the block, as for GretlMatrix */
    GretlMatrix opIndex(MatrixRange r, MatrixRange c) @trusted {
        return block(r, c).dup;
    }

    size_t opDollar(size_t dim)() {
        static if (dim == 0) return _nrow; else return _ncol;
    }

    MatrixRange opSlice(size_t dim)(size_t lo, size_t hi) {
        return MatrixRange(lo, hi);
    }

    /* a new GretlMatrix with a copy of the elements */
    GretlMatrix dup() @trusted {
        enforce!GretlException(_nrow <= int.max && _ncol <= int.max, "Submatrix: too large to copy");
        auto raw = gretl_matrix_alloc(cast(int) _nrow, cast(int) _ncol);

        enforce!GretlException(raw !is null || _nrow * _ncol == 0, "Submatrix: out of memory");
        if (raw is null) {
            return GretlMatrix.init;
        }
        foreach (j; 0 .. _ncol) {
            foreach (i; 0 .. _nrow) {
                raw.val[j * _nrow + i] = _data[j * _stride + i];
            }
        }
        return GretlMatrix(raw);
    }

    private Submatrix block(MatrixRange r, MatrixRange c) @trusted {
        enforce!GretlException(r.lo <= r.hi && r.hi <= _nrow, "Submatrix: row range " ~ to!string(r.lo) ~ " .. " ~ to!string(r.hi) ~ " is out of range");
        enforce!GretlException(c.lo <= c.hi && c.hi <= _ncol, "Submatrix: column range " ~ to!string(c.lo) ~ " .. " ~ to!string(c.hi) ~ " is out of range");
        Submatrix s;

        s._p = _p;
        s._data = _data + c.lo * _stride + r.lo;
        s._nrow = r.hi - r.lo;
        s._ncol = c.hi - c.lo;
        s._stride = _stride;
        return s;
    }
}

private Submatrix wholeView(GretlMatrix m) @trusted {
    Submatrix s;
    auto r = m.ptr;

    if (r !is null) {
        s._p = m._p;
        s._data = r.val;
        s._nrow = r.rows;
        s._ncol = r.cols;
        s._stride = r.rows;
    }
    return s;
}

private Submatrix wholeView(Submatrix s) {
    return s;
}

/* The explicit reference views. `sub`, `rows` and `cols` take no arguments and return something to index, so they
   read as m.sub[1 .. 8, 2 .. 4], m.rows[2 .. 6] and m.cols[2 .. 6]; they work on a GretlMatrix or a Submatrix. */
struct SubIndexer {
    private Submatrix s;

    Submatrix opIndex(MatrixRange r, MatrixRange c) {
        return s.block(r, c);
    }

    size_t opDollar(size_t dim)() {
        static if (dim == 0) return s._nrow; else return s._ncol;
    }

    MatrixRange opSlice(size_t dim)(size_t lo, size_t hi) {
        return MatrixRange(lo, hi);
    }
}

struct RowsIndexer {
    private Submatrix s;

    Submatrix opIndex(MatrixRange r) {
        return s.block(r, MatrixRange(0, s._ncol));
    }

    size_t opDollar(size_t dim)() {
        return s._nrow;
    }

    MatrixRange opSlice(size_t dim)(size_t lo, size_t hi) {
        return MatrixRange(lo, hi);
    }
}

struct ColsIndexer {
    private Submatrix s;

    Submatrix opIndex(MatrixRange c) {
        return s.block(MatrixRange(0, s._nrow), c);
    }

    size_t opDollar(size_t dim)() {
        return s._ncol;
    }

    MatrixRange opSlice(size_t dim)(size_t lo, size_t hi) {
        return MatrixRange(lo, hi);
    }
}

SubIndexer sub(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    return SubIndexer(wholeView(m));
}

RowsIndexer rows(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    return RowsIndexer(wholeView(m));
}

ColsIndexer cols(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    return ColsIndexer(wholeView(m));
}

/* Iteration over the rows or the columns, as byKey does for an associative array: m.byRow and m.byColumn (a GretlMatrix
   or a Submatrix) are random-access ranges whose elements are Submatrix references, a 1 x ncol view of each row or an
   nrow x 1 view of each column, so `foreach (row; m.byRow) row[0, 2] = 0;` changes m. They keep the data alive. */
struct ByRow {
    private Submatrix s;
    private size_t lo, hi;

    bool empty() { return lo >= hi; }
    size_t length() { return hi - lo; }
    Submatrix front() { return s.block(MatrixRange(lo, lo + 1), MatrixRange(0, s._ncol)); }
    Submatrix back() { return s.block(MatrixRange(hi - 1, hi), MatrixRange(0, s._ncol)); }
    void popFront() { lo++; }
    void popBack() { hi--; }
    ByRow save() { return this; }

    Submatrix opIndex(size_t i) {
        enforce!GretlException(i < length, "ByRow: index " ~ to!string(i) ~ " is out of range");
        return s.block(MatrixRange(lo + i, lo + i + 1), MatrixRange(0, s._ncol));
    }
}

struct ByColumn {
    private Submatrix s;
    private size_t lo, hi;

    bool empty() { return lo >= hi; }
    size_t length() { return hi - lo; }
    Submatrix front() { return s.block(MatrixRange(0, s._nrow), MatrixRange(lo, lo + 1)); }
    Submatrix back() { return s.block(MatrixRange(0, s._nrow), MatrixRange(hi - 1, hi)); }
    void popFront() { lo++; }
    void popBack() { hi--; }
    ByColumn save() { return this; }

    Submatrix opIndex(size_t i) {
        enforce!GretlException(i < length, "ByColumn: index " ~ to!string(i) ~ " is out of range");
        return s.block(MatrixRange(0, s._nrow), MatrixRange(lo + i, lo + i + 1));
    }
}

ByRow byRow(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    auto v = wholeView(m);
    return ByRow(v, 0, v._nrow);
}

ByColumn byColumn(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    auto v = wholeView(m);
    return ByColumn(v, 0, v._ncol);
}

/* One row or column as a range of doubles, like a double[]: foreach (val; Row(m, 1)), foreach (j, val; Col(m, 2)),
   `ref` to change the matrix (foreach (ref val; Row(m, 1)) val *= 2). It is a reference into the matrix (a column is
   one run of memory, a row is spaced by the stride), keeps the data alive, and is also a random-access range with
   length, indexing, front/back and `dup` (a double[] copy). */
struct MatrixLine {
    private Submatrix s;        // a 1 x n or n x 1 block, for the data pointer and to keep the data alive
    private size_t n, step;     // number of elements and distance between them in memory
    private size_t lo, hi;      // the part still in the range

    size_t length() { return hi - lo; }
    bool empty() { return lo >= hi; }
    size_t opDollar() { return hi - lo; }

    private double * at(size_t i) @trusted { return s._data + (lo + i) * step; }

    ref double front() @trusted { return *at(0); }
    ref double back() @trusted { return *at(length - 1); }
    void popFront() { lo++; }
    void popBack() { hi--; }
    MatrixLine save() { return this; }

    ref double opIndex(size_t i) @trusted {
        enforce!GretlException(i < length, "Row/Col: index " ~ to!string(i) ~ " is out of range");
        return *at(i);
    }

    double opIndexAssign(double x, size_t i) @trusted {
        enforce!GretlException(i < length, "Row/Col: index " ~ to!string(i) ~ " is out of range");
        *at(i) = x;
        return x;
    }

    int opApply(scope int delegate(ref double) @system dg) @system {
        foreach (i; 0 .. length) {
            if (auto r = dg(*at(i))) return r;
        }
        return 0;
    }

    int opApply(scope int delegate(size_t, ref double) @system dg) @system {
        foreach (i; 0 .. length) {
            if (auto r = dg(i, *at(i))) return r;
        }
        return 0;
    }

    int opApply(scope int delegate(ref double) @safe dg) @trusted {
        foreach (i; 0 .. length) {
            if (auto r = dg(*at(i))) return r;
        }
        return 0;
    }

    int opApply(scope int delegate(size_t, ref double) @safe dg) @trusted {
        foreach (i; 0 .. length) {
            if (auto r = dg(i, *at(i))) return r;
        }
        return 0;
    }

    /* Assigning to a line writes into the matrix (it does not re-point the line), so diag(m) = [1.1, 3.4, 1.7] and
       diag(m) = Row(m, 1) work: from a double[] or another line of the same length, or one scalar for every element. */
    void opAssign(double[] src) @trusted {
        enforce!GretlException(src.length == length, "assignment: " ~ to!string(src.length) ~ " values for " ~ to!string(length) ~ " elements");
        foreach (i; 0 .. length) *at(i) = src[i];
    }

    void opAssign(MatrixLine src) @trusted {
        opAssign(src.dup);          // a copy first, in case the two overlap
    }

    void opAssign(double x) @trusted {
        foreach (i; 0 .. length) *at(i) = x;
    }

    double[] dup() @trusted {
        auto a = new double[length];

        foreach (i; 0 .. length) a[i] = *at(i);
        return a;
    }
}

/* row i of a GretlMatrix or Submatrix */
MatrixLine Row(M)(M m, size_t i) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    auto v = wholeView(m);
    auto b = v.block(MatrixRange(i, i + 1), MatrixRange(0, v._ncol));
    return MatrixLine(b, v._ncol, v._stride, 0, v._ncol);
}

/* column j of a GretlMatrix or Submatrix */
MatrixLine Col(M)(M m, size_t j) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    auto v = wholeView(m);
    auto b = v.block(MatrixRange(0, v._nrow), MatrixRange(j, j + 1));
    return MatrixLine(b, v._nrow, 1, 0, v._nrow);
}

/* All the elements of a matrix: m.el (a GretlMatrix or a Submatrix) is a reference to be iterated over, in memory order
   (column by column), so the order is not part of the contract.
       foreach (val; m.el) ...            foreach (ref val; m.el) val *= 2;   // changes m
       foreach (el, val; m.el) ...        el.row and el.col are the (size_t) row and column of val */
struct ElementIndex {
    size_t row, col;
}

struct Elements {
    private Submatrix s;

    size_t length() { return s._nrow * s._ncol; }

    /* element by element: m.el * mm, m.el * mm.el, m.el * a, and the same for /; + and - are element by element
       anyway, so m.el + mm is m + mm */
    GretlMatrix opBinary(string op, T)(T rhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && isNumeric!T) {
        return scalarRight!op(s, cast(double) rhs);
    }

    GretlMatrix opBinary(string op, T)(T rhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && (isMatrixLike!T || is(T == Elements))) {
        return elementwise!op(s, asView(rhs));
    }

    /* element powers: m.el ^^ 2 squares each element, m.el ^^ mm, a ^^ m.el */
    GretlMatrix opBinary(string op, T)(T rhs) if (op == "^^" && isNumeric!T) {
        return scalarRight!op(s, cast(double) rhs);
    }

    GretlMatrix opBinary(string op, T)(T rhs) if (op == "^^" && (isMatrixLike!T || is(T == Elements))) {
        return elementwise!op(s, asView(rhs));
    }

    GretlMatrix opBinaryRight(string op, T)(T lhs) if (op == "^^" && isNumeric!T) {
        return scalarLeft!op(cast(double) lhs, s);
    }

    GretlMatrix opBinaryRight(string op, T)(T lhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && isNumeric!T) {
        return scalarLeft!op(cast(double) lhs, s);
    }

    int opApply(scope int delegate(ref double) @system dg) @system {
        foreach (j; 0 .. s._ncol) {
            double * p = s._data + j * s._stride;

            foreach (i; 0 .. s._nrow) {
                if (auto r = dg(p[i])) return r;
            }
        }
        return 0;
    }

    int opApply(scope int delegate(ElementIndex, ref double) @system dg) @system {
        foreach (j; 0 .. s._ncol) {
            double * p = s._data + j * s._stride;

            foreach (i; 0 .. s._nrow) {
                if (auto r = dg(ElementIndex(i, j), p[i])) return r;
            }
        }
        return 0;
    }

    int opApply(scope int delegate(ref double) @safe dg) @trusted {
        foreach (j; 0 .. s._ncol) {
            double * p = s._data + j * s._stride;

            foreach (i; 0 .. s._nrow) {
                if (auto r = dg(p[i])) return r;
            }
        }
        return 0;
    }

    int opApply(scope int delegate(ElementIndex, ref double) @safe dg) @trusted {
        foreach (j; 0 .. s._ncol) {
            double * p = s._data + j * s._stride;

            foreach (i; 0 .. s._nrow) {
                if (auto r = dg(ElementIndex(i, j), p[i])) return r;
            }
        }
        return 0;
    }
}

Elements el(M)(M m) if (is(M == GretlMatrix) || is(M == Submatrix)) {
    return Elements(wholeView(m));
}

/* Operators. Every result is a new GretlMatrix. In the table m and mm are a GretlMatrix or a Submatrix, a is a scalar:
       m + mm  m - mm           element by element, same size
       m + a   a + m   m - a   a - m
       m * a   a * m   m / a   a / m     (a / m divides a by each element)
       m * mm                   matrix multiplication (ncol of m equals nrow of mm)
       m / mm                   matrix right division as in Matlab: the x that solves x * mm = m, that is m * inv(mm) for a
                                square mm and the least-squares solution (QR) when mm has more columns than rows
       m.el * mm   m * mm.el   m.el * mm.el   m.el * a        element by element
       m.el / mm   m / mm.el   m.el / mm.el   m.el / a  a / m.el
       m.el ^^ a   a ^^ m.el   m.el ^^ mm   m.el ^^ mm.el        element powers; m ^^ n alone is refused (see power)
       m.el + mm, m + mm.el, ...  the same as without .el (they are element by element anyway)
   Sizes that do not fit throw GretlException. */
enum isMatrixLike(T) = is(T == GretlMatrix) || is(T == Submatrix);

mixin template MatrixOps() {
    GretlMatrix opBinary(string op, T)(T rhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && isNumeric!T) {
        return scalarRight!op(wholeView(this), cast(double) rhs);
    }

    GretlMatrix opBinary(string op, T)(T rhs) if ((op == "+" || op == "-") && isMatrixLike!T) {
        return elementwise!op(wholeView(this), asView(rhs));
    }

    GretlMatrix opBinary(string op, T)(T rhs) if (op == "*" && isMatrixLike!T) {
        static if (is(typeof(this) == GretlMatrix) && is(T == GretlMatrix)) {
            return matmulRaw(this, rhs);
        } else {
            return matmul(wholeView(this), asView(rhs));
        }
    }

    GretlMatrix opBinary(string op, T)(T rhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && is(T == Elements)) {
        return elementwise!op(wholeView(this), asView(rhs));
    }

    GretlMatrix opBinary(string op, T)(T rhs) if (op == "/" && isMatrixLike!T) {
        return rdivide(wholeView(this), asView(rhs));
    }

    GretlMatrix opBinaryRight(string op, T)(T lhs) if ((op == "+" || op == "-" || op == "*" || op == "/") && isNumeric!T) {
        return scalarLeft!op(cast(double) lhs, wholeView(this));
    }

    /* m ~ x: the columns of x added to the right of m (cbind); x is a matrix, a scalar or an array, see cbind */
    GretlMatrix opBinary(string op, T)(T rhs) if (op == "~" && isCbindOperand!T) {
        return cbind(this, rhs);
    }

    GretlMatrix opBinaryRight(string op, T)(T lhs) if (op == "~" && (isNumeric!T || isNumericArray!T)) {
        return cbind(lhs, this);
    }

    /* ^^ on a matrix itself is refused, because it could mean either of two things */
    void opBinary(string op, T)(T rhs) if (op == "^^") {
        static assert(false, "m ^^ n is not defined: use m.el ^^ n for the power of each element, or power(m, n) for the matrix power");
    }

    void opBinaryRight(string op, T)(T lhs) if (op == "^^") {
        static assert(false, "a ^^ m is not defined: use a ^^ m.el for each element as the exponent");
    }
}

private Submatrix asView(T)(T x) {
    static if (is(T == GretlMatrix)) return wholeView(x);
    else static if (is(T == Submatrix)) return x;
    else static if (is(T == Elements)) return x.s;
    else static assert(false, "not a matrix type");
}

/* a new nr x nc matrix with unspecified contents (an empty GretlMatrix if there are no elements) */
private GretlMatrix newMatrix(size_t nr, size_t nc) @trusted {
    enforce!GretlException(nr <= int.max && nc <= int.max, "GretlMatrix: too large");
    auto raw = gretl_matrix_alloc(cast(int) nr, cast(int) nc);

    if (raw is null) {
        enforce!GretlException(nr * nc == 0, "GretlMatrix: out of memory");
        return GretlMatrix.init;
    }
    return GretlMatrix(raw);
}

private string sizeMsg(string op, Submatrix a, Submatrix b) {
    return "operator " ~ op ~ ": sizes " ~ to!string(a._nrow) ~ " x " ~ to!string(a._ncol) ~ " and "
        ~ to!string(b._nrow) ~ " x " ~ to!string(b._ncol) ~ " do not fit";
}

private GretlMatrix elementwise(string op)(Submatrix a, Submatrix b) @trusted {
    enforce!GretlException(a._nrow == b._nrow && a._ncol == b._ncol, sizeMsg(op, a, b));
    auto c = newMatrix(a._nrow, a._ncol);

    if (a._nrow * a._ncol == 0) return c;
    double * z = c.ptr.val;

    foreach (j; 0 .. a._ncol) {
        const double * pa = a._data + j * a._stride;
        const double * pb = b._data + j * b._stride;

        foreach (i; 0 .. a._nrow) {
            z[j * a._nrow + i] = mixin("pa[i] " ~ op ~ " pb[i]");
        }
    }
    return c;
}

private GretlMatrix scalarRight(string op)(Submatrix a, double x) @trusted {
    auto c = newMatrix(a._nrow, a._ncol);

    if (a._nrow * a._ncol == 0) return c;
    double * z = c.ptr.val;

    foreach (j; 0 .. a._ncol) {
        const double * pa = a._data + j * a._stride;

        foreach (i; 0 .. a._nrow) {
            z[j * a._nrow + i] = mixin("pa[i] " ~ op ~ " x");
        }
    }
    return c;
}

private GretlMatrix scalarLeft(string op)(double x, Submatrix a) @trusted {
    auto c = newMatrix(a._nrow, a._ncol);

    if (a._nrow * a._ncol == 0) return c;
    double * z = c.ptr.val;

    foreach (j; 0 .. a._ncol) {
        const double * pa = a._data + j * a._stride;

        foreach (i; 0 .. a._nrow) {
            z[j * a._nrow + i] = mixin("x " ~ op ~ " pa[i]");
        }
    }
    return c;
}

/* matrix product of two views (for views and for sizes that need no BLAS) */
private GretlMatrix matmul(Submatrix a, Submatrix b) @trusted {
    enforce!GretlException(a._ncol == b._nrow, sizeMsg("*", a, b));
    auto c = newMatrix(a._nrow, b._ncol);

    if (a._nrow * b._ncol == 0) return c;
    double * z = c.ptr.val;

    z[0 .. a._nrow * b._ncol] = 0.0;
    foreach (j; 0 .. b._ncol) {
        foreach (k; 0 .. a._ncol) {
            const double bkj = b._data[j * b._stride + k];
            const double * pa = a._data + k * a._stride;

            foreach (i; 0 .. a._nrow) {
                z[j * a._nrow + i] += pa[i] * bkj;
            }
        }
    }
    return c;
}

/* product of two whole GretlMatrixes, with the library's multiplication (BLAS) */
private GretlMatrix matmulRaw(GretlMatrix a, GretlMatrix b) @trusted {
    auto pa = a.ptr, pb = b.ptr;

    if (pa is null || pb is null || pa.rows == 0 || pa.cols == 0 || pb.cols == 0) {
        return matmul(wholeView(a), wholeView(b));
    }
    enforce!GretlException(pa.cols == pb.rows, sizeMsg("*", wholeView(a), wholeView(b)));
    auto c = newMatrix(pa.rows, pb.cols);
    int err = gretl_matrix_multiply(pa, pb, c.ptr);

    enforce!GretlException(err == 0, "operator *: gretl_matrix_multiply failed");
    return c;
}

private GretlMatrix transposeCopy(Submatrix a) @trusted {
    auto t = newMatrix(a._ncol, a._nrow);

    if (a._nrow * a._ncol == 0) return t;
    double * z = t.ptr.val;

    foreach (j; 0 .. a._ncol) {
        foreach (i; 0 .. a._nrow) {
            z[i * a._ncol + j] = a._data[j * a._stride + i];
        }
    }
    return t;
}

/* a / b, Matlab's mrdivide: x with x * b = a, as (b' \ a')'. With a (m x n) and b (p x n), the regression of the n rows of
   a' on the p columns of b' by QR, so for a row vector y (1 x n) and b = X' (p x n), y / X' is the OLS coefficient row.
   Needs n >= p and b' of full column rank. */
private GretlMatrix rdivide(Submatrix a, Submatrix b) @trusted {
    enforce!GretlException(a._ncol == b._ncol, sizeMsg("/", a, b));
    enforce!GretlException(b._ncol >= b._nrow, "operator /: the divisor has more rows (" ~ to!string(b._nrow) ~ ") than columns ("
                           ~ to!string(b._ncol) ~ "), the system is underdetermined");
    if (a._nrow * b._nrow == 0) return newMatrix(a._nrow, b._nrow);

    auto yt = transposeCopy(a);          // n x m
    auto xt = transposeCopy(b);          // n x p
    auto coef = newMatrix(b._nrow, a._nrow);    // p x m
    int err = gretl_matrix_QR_ols(yt.ptr, xt.ptr, coef.ptr, null, null, null);

    enforce!GretlException(err == 0, "operator /: the least-squares solve failed (is the divisor rank deficient?), error code " ~ to!string(err));
    return transposeCopy(wholeView(coef));
}

/* solve(X, y): Matlab's X \ y for data in columns, X (n x p) and y (n x m); the p x m matrix b that minimizes the sum of
   squares of y - X b, by QR. For a vector y (n x 1) these are the OLS coefficients, a p x 1 matrix. Each argument is a
   GretlMatrix or a Submatrix. Needs n >= p and X of full column rank, else GretlException. */
GretlMatrix solve(MX, MY)(MX X, MY y) @trusted if (isMatrixLike!MX && isMatrixLike!MY) {
    auto xv = wholeView(X), yv = wholeView(y);

    enforce!GretlException(xv._nrow == yv._nrow, sizeMsg("solve", xv, yv));
    enforce!GretlException(xv._nrow >= xv._ncol, "solve: X has more columns (" ~ to!string(xv._ncol) ~ ") than rows ("
                           ~ to!string(xv._nrow) ~ "), the system is underdetermined");
    if (xv._ncol * yv._ncol == 0) return newMatrix(xv._ncol, yv._ncol);

    static if (is(MX == GretlMatrix)) { auto xm = X; } else { auto xm = xv.dup; }
    static if (is(MY == GretlMatrix)) { auto ym = y; } else { auto ym = yv.dup; }
    auto b = newMatrix(xv._ncol, yv._ncol);
    int err = gretl_matrix_QR_ols(ym.ptr, xm.ptr, b.ptr, null, null, null);

    enforce!GretlException(err == 0, "solve: the least-squares solve failed (is X rank deficient?), error code " ~ to!string(err));
    return b;
}

/* t(m): the transpose, as a new GretlMatrix. */
GretlMatrix t(M)(M m) if (isMatrixLike!M) {
    return transposeCopy(wholeView(m));
}

/* inv(m): the inverse of a square matrix, as a new GretlMatrix (gretl_invert_matrix: Cholesky for a symmetric positive
   definite matrix, LU otherwise). A singular matrix throws GretlException. */
GretlMatrix inv(M)(M m) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "inv: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    auto c = v.dup;

    if (v._nrow == 0) return c;
    int err = gretl_invert_matrix(c.ptr);

    enforce!GretlException(err == 0, "inv: the matrix could not be inverted (singular?), error code " ~ to!string(err));
    return c;
}

/* kron(a, b), also a.kron(b): the Kronecker product, (a.nrow * b.nrow) x (a.ncol * b.ncol), a new GretlMatrix. */
GretlMatrix kron(MA, MB)(MA a, MB b) @trusted if (isMatrixLike!MA && isMatrixLike!MB) {
    auto av = wholeView(a), bv = wholeView(b);
    auto k = newMatrix(av._nrow * bv._nrow, av._ncol * bv._ncol);

    if (k.nrow * k.ncol == 0) return k;
    double * z = k.ptr.val;
    const size_t kr = k.nrow;

    foreach (j; 0 .. av._ncol) {
        foreach (l; 0 .. bv._ncol) {
            foreach (i; 0 .. av._nrow) {
                const double aij = av._data[j * av._stride + i];

                foreach (r; 0 .. bv._nrow) {
                    z[(j * bv._ncol + l) * kr + i * bv._nrow + r] = aij * bv._data[l * bv._stride + r];
                }
            }
        }
    }
    return k;
}

/* power(m, n): the matrix power m * m * ... * m (n factors) of a square matrix, as a new GretlMatrix, by repeated squaring;
   power(m, 2) is m * m, power(m, 0) the identity, and a negative n is the power of inv(m). */
GretlMatrix power(M)(M m, long n) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "power: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (n < 0) {
        enforce!GretlException(n != long.min, "power: exponent out of range");
        return power(inv(m), -n);
    }
    auto result = newMatrix(v._nrow, v._ncol);

    if (v._nrow == 0) return result;
    result.ptr.val[0 .. v._nrow * v._nrow] = 0.0;
    foreach (i; 0 .. v._nrow) result.ptr.val[i * v._nrow + i] = 1.0;
    if (n == 0) return result;
    auto base = v.dup;
    bool first = true;

    while (n > 0) {
        if (n & 1) {
            result = first ? base.dup : result * base;
            first = false;
        }
        n >>= 1;
        if (n > 0) base = base * base;
    }
    return result;
}

/* Which triangle a factor is stored in. */
enum Triangle { upper, lower }

/* chol(m), m.chol: the Cholesky factor of a symmetric positive definite matrix, as a new GretlMatrix. The default is
   the UPPER triangle U with U'U = m (as in R and Matlab); chol(m, Triangle.lower) is the LOWER triangle L with L L' = m
   (as gretl's own cholesky() returns it), which is the transpose of the upper one. Only the lower triangle of m is read,
   symmetry is not checked. A matrix that is not positive definite throws GretlException. */
GretlMatrix chol(M)(M m, Triangle tri = Triangle.upper) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "chol: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    auto c = v.dup;

    if (v._nrow == 0) return c;
    int err = gretl_matrix_cholesky_decomp(c.ptr);

    enforce!GretlException(err == 0, "chol: the matrix is not positive definite, error code " ~ to!string(err));
    return (tri == Triangle.lower) ? c : transposeCopy(wholeView(c));
}

/* diag(m): the diagonal of a GretlMatrix or a Submatrix, as a reference (a MatrixLine with min(nrow, ncol) elements), so that
   diag(m) = [1.1, 3.4, 1.7], diag(m) = Row(m, 1), diag(m) = 0 and foreach (ref val; diag(m)) change the matrix. */
MatrixLine diag(M)(M m) if (isMatrixLike!M) {
    auto v = wholeView(m);
    const size_t n = (v._nrow < v._ncol) ? v._nrow : v._ncol;

    return MatrixLine(v, n, v._stride + 1, 0, n);
}

/* The result of eigen: `values` is an n x 1 matrix of eigenvalues in ascending order, `vectors` the n x n matrix whose
   columns are the matching orthonormal eigenvectors (an empty GretlMatrix if they were not asked for). */
struct Eigen {
    GretlMatrix values;
    GretlMatrix vectors;
}

/* eigen(m), m.eigen: the eigen decomposition of a real SYMMETRIC matrix, so that m * v = v * diag(values) with
   v = vectors and t(v) * v the identity. The eigenvalues are real and come in ascending order; the sign of each eigenvector
   is whatever LAPACK (dsyev) returns, so only compare them up to sign. `eigen(m, false)` skips the eigenvectors, which
   is cheaper. A matrix that is not square or not symmetric (the same to 1e-12 relative to its largest element) throws
   GretlException; the general (non-symmetric) case, with complex eigenvalues, is not available yet. */
Eigen eigen(M)(M m, bool vectors = true) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "eigen: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    double big = 0;

    foreach (j; 0 .. n) {
        foreach (i; 0 .. n) {
            const double x = v._data[j * v._stride + i];

            if (x > big) big = x; else if (-x > big) big = -x;
        }
    }
    foreach (j; 0 .. n) {
        foreach (i; j + 1 .. n) {
            const double d = v._data[j * v._stride + i] - v._data[i * v._stride + j];

            enforce!GretlException((d < 0 ? -d : d) <= 1e-12 * big, "eigen: the matrix is not symmetric (elements [" ~ to!string(i) ~ ", "
                                   ~ to!string(j) ~ "] and [" ~ to!string(j) ~ ", " ~ to!string(i) ~ "] differ)");
        }
    }
    Eigen result;

    if (n == 0) return result;
    auto work = v.dup;
    int err = 0;
    auto vals = gretl_symmetric_matrix_eigenvals(work.ptr, vectors ? 1 : 0, &err);

    enforce!GretlException(err == 0 && vals !is null, "eigen: the decomposition failed, error code " ~ to!string(err));
    result.values = GretlMatrix(vals);
    if (vectors) result.vectors = work;
    return result;
}

/* The result of eigenGeneral. The eigenvalues of a real non-symmetric matrix can be complex, and GretlMatrix holds only real
   numbers, so every complex quantity comes as a pair of real matrices:
       values        n x 2: the real parts in column 0, the imaginary parts in column 1 (complex ones come in conjugate pairs)
       vectorsRe,
       vectorsIm     n x n: the real and the imaginary parts of the right eigenvectors, column j for eigenvalue j
   `isReal` is true if every eigenvalue is real (then `vectorsIm` is all zero). */
struct EigenGeneral {
    GretlMatrix values;
    GretlMatrix vectorsRe;
    GretlMatrix vectorsIm;

    bool isReal() @trusted {
        foreach (i; 0 .. values.nrow) {
            if (values[i, 1] != 0.0) return false;
        }
        return true;
    }
}

/* eigenGeneral(m), m.eigenGeneral: the eigen decomposition of a general real square matrix (LAPACK dgeev, through
   gretl_dgeev): m v = lambda v for each right eigenvector v. The eigenvalues are NOT sorted (LAPACK's order, conjugate pairs
   adjacent, positive imaginary part first); each eigenvector has unit Euclidean norm and its largest component real. For
   a symmetric matrix use eigen, which returns real, sorted values and orthogonal vectors. `eigenGeneral(m, false)` skips
   the vectors (`vectorsRe` and `vectorsIm` are then empty). A non-square matrix throws GretlException. */
EigenGeneral eigenGeneral(M)(M m, bool vectors = true) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "eigenGeneral: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    EigenGeneral result;

    if (n == 0) return result;
    result.values = newMatrix(n, 2);
    if (n == 1) {
        result.values[0, 0] = v._data[0];
        result.values[0, 1] = 0.0;
        if (vectors) {
            result.vectorsRe = newMatrix(1, 1);
            result.vectorsIm = newMatrix(1, 1);
            result.vectorsRe[0, 0] = 1.0;
            result.vectorsIm[0, 0] = 0.0;
        }
        return result;
    }
    auto a = v.dup;
    GretlMatrixRaw * vr = vectors ? gretl_matrix_alloc(1, 1) : null;   // becomes the complex n x n eigenvector matrix
    int err = 0;
    GretlMatrixRaw * lam = gretl_dgeev(a.ptr, vr, null, &err);

    scope (exit) {
        gretl_matrix_free(lam);
        gretl_matrix_free(vr);
    }
    enforce!GretlException(err == 0 && lam !is null, "eigenGeneral: the decomposition failed, error code " ~ to!string(err));
    enforce!GretlException(lam.is_complex && lam.rows == n, "eigenGeneral: unexpected result from gretl_dgeev");
    foreach (i; 0 .. n) {
        const z = gretl_cmatrix_get(lam, cast(int) i, 0);

        result.values[i, 0] = z.r;
        result.values[i, 1] = z.i;
    }
    if (vectors) {
        enforce!GretlException(vr.is_complex && vr.rows == n && vr.cols == n, "eigenGeneral: unexpected eigenvector result from gretl_dgeev");
        result.vectorsRe = newMatrix(n, n);
        result.vectorsIm = newMatrix(n, n);
        foreach (j; 0 .. n) {
            foreach (i; 0 .. n) {
                const z = gretl_cmatrix_get(vr, cast(int) i, cast(int) j);

                result.vectorsRe[i, j] = z.r;
                result.vectorsIm[i, j] = z.i;
            }
        }
    }
    return result;
}

/* The result of svd: m = u * diag(s) * t(v) (with k = min(nrow, ncol) singular values, or all of u and v if `full`).
       s   k x 1, the singular values, largest first
       u   nrow x k   (nrow x nrow if full)       orthonormal columns: the left singular vectors
       v   ncol x k   (ncol x ncol if full)       orthonormal columns: the right singular vectors (V, not V')
   u and v are empty if only the singular values were asked for. */
struct SVD {
    GretlMatrix u;
    GretlMatrix s;
    GretlMatrix v;
}

/* svd(m), m.svd: the singular value decomposition, m = u * diag(s) * t(v), with the singular values largest first.
   svd(m, false) gives only the singular values (cheaper); svd(m, true, true) gives the full u (nrow x nrow) and v (ncol x
   ncol) when m is not square. The signs of the singular vectors are LAPACK's.

   For a matrix with more rows than columns the library's own gretl_matrix_SVD would use a fast method (eigenvalues of
   m'm) that squares the condition number, so small singular values lose their accuracy. To avoid that this function
   decomposes t(m) in that case and swaps the roles of u and v: the same result, computed by LAPACK's real SVD routines. */
SVD svd(M)(M m, bool vectors = true, bool full = false) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);
    SVD result;

    if (v._nrow * v._ncol == 0) return result;
    const bool tall = v._nrow > v._ncol;
    auto x = tall ? transposeCopy(v) : v.dup;       // the wide or square matrix that is handed to the library
    GretlMatrixRaw * pu = null, ps = null, pvt = null;
    int err = gretl_matrix_SVD(x.ptr, vectors ? &pu : null, &ps, vectors ? &pvt : null, full ? 1 : 0);

    scope (exit) {
        gretl_matrix_free(pu);
        gretl_matrix_free(ps);
        gretl_matrix_free(pvt);
    }
    enforce!GretlException(err == 0 && ps !is null, "svd: the decomposition failed, error code " ~ to!string(err));
    const k = cast(size_t) ps.rows * ps.cols;

    result.s = newMatrix(k, 1);
    foreach (i; 0 .. k) result.s.ptr.val[i] = ps.val[i];
    if (vectors) {
        enforce!GretlException(pu !is null && pvt !is null, "svd: the singular vectors were not returned");
        // x = pu * diag(s) * pvt; for the tall case m = t(x), so u = t(pvt) and v = pu
        auto left = GretlMatrix(gretl_matrix_copy(pu));
        auto right = transposeCopy(rawView(pvt));
        if (tall) {
            result.u = right;
            result.v = left;
        } else {
            result.u = left;
            result.v = right;
        }
    }
    return result;
}

/* a view of a raw matrix (no ownership) */
private Submatrix rawView(GretlMatrixRaw * raw) @trusted {
    Submatrix s;

    s._data = raw.val;
    s._nrow = raw.rows;
    s._ncol = raw.cols;
    s._stride = raw.rows;
    return s;
}

/* The result of qr: m = q * r, with q nrow x ncol having orthonormal columns and r ncol x ncol upper triangular. */
struct QR {
    GretlMatrix q;
    GretlMatrix r;
}

/* qr(m), m.qr: the (thin) QR decomposition m = q * r by Householder reflections (LAPACK dgeqrf and dorgqr, through
   gretl_matrix_QR_decomp). q has the shape of m and orthonormal columns, r is the ncol x ncol upper triangle; the
   diagonal of r can have either sign. Needs nrow >= ncol (the library does not decompose wide matrices), else
   GretlException. No pivoting. */
QR qr(M)(M m) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow >= v._ncol, "qr: the matrix has more columns (" ~ to!string(v._ncol) ~ ") than rows ("
                           ~ to!string(v._nrow) ~ "), it must have rows >= columns");
    QR result;

    if (v._nrow * v._ncol == 0) {
        result.q = newMatrix(v._nrow, v._ncol);
        result.r = newMatrix(v._ncol, v._ncol);
        return result;
    }
    auto q = v.dup;                                 // overwritten by Q
    auto r = newMatrix(v._ncol, v._ncol);
    int err = gretl_matrix_QR_decomp(q.ptr, r.ptr);

    enforce!GretlException(err == 0, "qr: the decomposition failed, error code " ~ to!string(err));
    result.q = q;
    result.r = r;
    return result;
}

/* The result of lu: p * m = l * u, with k = min(nrow, ncol)
       l      nrow x k   unit lower triangular (ones on the diagonal, zeros above)
       u      k x ncol   upper triangular (upper trapezoidal if not square)
       p      nrow x nrow  permutation matrix of the row interchanges
       perm   the same permutation as indices: row i of p * m is row perm[i] of m (0-based)
       singular  true if an exact zero turned up on the diagonal of u: the factors are still valid but u cannot be used to solve */
struct LU {
    GretlMatrix l;
    GretlMatrix u;
    GretlMatrix p;
    size_t[] perm;
    bool singular;
}

/* lu(m), m.lu: the LU decomposition with partial pivoting, p * m = l * u (LAPACK dgetrf, called directly because the library has
   no function that returns the factors). Works for any shape, square or not. The input is copied, not changed. */
LU lu(M)(M m) @trusted if (isMatrixLike!M) {
    import gretl.clapack_doubleh: dgetrf_;
    import gretl.gretl_f2ch: integer;

    auto v = wholeView(m);
    const nr = v._nrow, nc = v._ncol, k = nr < nc ? nr : nc;
    LU result;

    if (nr * nc == 0) {
        result.l = newMatrix(nr, k);
        result.u = newMatrix(k, nc);
        result.p = newMatrix(nr, nr);
        foreach (i; 0 .. nr) result.p[i, i] = 1.0;
        foreach (i; 0 .. nr) result.perm ~= i;
        return result;
    }
    enforce!GretlException(nr <= int.max && nc <= int.max, "lu: the matrix is too large");
    auto a = v.dup;
    auto ipiv = new integer[k];
    integer mm = cast(integer) nr, nn = cast(integer) nc, lda = cast(integer) nr, info = 0;

    dgetrf_(&mm, &nn, a.ptr.val, &lda, ipiv.ptr, &info);
    enforce!GretlException(info >= 0, "lu: dgetrf failed, info = " ~ to!string(info));
    result.singular = info > 0;
    result.l = newMatrix(nr, k);
    result.u = newMatrix(k, nc);
    foreach (j; 0 .. k) {
        foreach (i; 0 .. nr) {
            result.l.ptr.val[j * nr + i] = (i == j) ? 1.0 : (i > j ? a.ptr.val[j * nr + i] : 0.0);
        }
    }
    foreach (j; 0 .. nc) {
        foreach (i; 0 .. k) {
            result.u.ptr.val[j * k + i] = (i <= j) ? a.ptr.val[j * nr + i] : 0.0;
        }
    }
    result.perm = new size_t[nr];
    foreach (i; 0 .. nr) result.perm[i] = i;
    foreach (i; 0 .. k) {
        const size_t r = cast(size_t) ipiv[i] - 1;      // dgetrf interchanged rows i and r, in this order

        if (r != i) {
            const tmp = result.perm[i];

            result.perm[i] = result.perm[r];
            result.perm[r] = tmp;
        }
    }
    result.p = newMatrix(nr, nr);
    result.p.ptr.val[0 .. nr * nr] = 0.0;
    foreach (i; 0 .. nr) result.p[i, result.perm[i]] = 1.0;
    return result;
}

/* The result of schur: m = z * t * t(z) (the real Schur decomposition)
       t       n x n   quasi upper triangular: upper triangular except for 2 x 2 blocks on the diagonal, one for each
               complex conjugate pair of eigenvalues
       z       n x n   orthogonal: the Schur vectors (empty if they were not asked for)
       values  n x 2   the eigenvalues read off t, real parts and imaginary parts, in the order of the diagonal of t */
struct Schur {
    GretlMatrix t;
    GretlMatrix z;
    GretlMatrix values;
}

/* schur(m), m.schur: the real Schur decomposition of a real square matrix, m = z * t * t(z) (LAPACK dgees, called directly
   because the library only has the complex version, gretl_zgees). The eigenvalues are not sorted. schur(m, false) skips z.
   A non-square matrix throws GretlException, as does a failure of the QR iteration to converge. */
Schur schur(M)(M m, bool vectors = true) @trusted if (isMatrixLike!M) {
    import gretl.clapack_doubleh: dgees_;
    import gretl.gretl_f2ch: integer, logical;

    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "schur: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    Schur result;

    if (n == 0) return result;
    enforce!GretlException(n <= int.max, "schur: the matrix is too large");
    auto a = v.dup;                                     // becomes T
    auto z = vectors ? newMatrix(n, n) : GretlMatrix.init;
    auto vals = newMatrix(n, 2);
    char jobvs = vectors ? 'V' : 'N';
    char sort = 'N';
    integer nn = cast(integer) n, lda = cast(integer) n, sdim = 0, ldvs = vectors ? cast(integer) n : 1, lwork = -1, info = 0;
    double * wr = vals.ptr.val, wi = wr + n;
    double * vs = vectors ? z.ptr.val : null;
    double[1] query;

    dgees_(&jobvs, &sort, null, &nn, a.ptr.val, &lda, &sdim, wr, wi, vs, &ldvs, query.ptr, &lwork, null, &info);
    enforce!GretlException(info == 0, "schur: the workspace query failed, info = " ~ to!string(info));
    lwork = cast(integer) query[0];
    auto work = new double[lwork];

    dgees_(&jobvs, &sort, null, &nn, a.ptr.val, &lda, &sdim, wr, wi, vs, &ldvs, work.ptr, &lwork, null, &info);
    enforce!GretlException(info == 0, "schur: the decomposition failed (the QR iteration did not converge), info = " ~ to!string(info));
    result.t = a;
    result.z = z;
    result.values = vals;
    return result;
}

/* pinv(m), m.pinv: the Moore-Penrose pseudoinverse, an ncol x nrow matrix, from the SVD (gretl_matrix_moore_penrose): the
   singular values that are not larger than tol * max(nrow, ncol) * (largest singular value) count as zero. The default
   tol = 0 means machine epsilon (2.2e-16), as in numpy and Matlab. For a square nonsingular matrix it is inv(m); for a
   tall matrix of full column rank it is inv(t(m) * m) * t(m); for the zero matrix it is the zero matrix. */
GretlMatrix pinv(M)(M m, double tol = 0.0) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    if (v._nrow * v._ncol == 0) return newMatrix(v._ncol, v._nrow);
    bool allZero = true;

    foreach (j; 0 .. v._ncol) {
        foreach (i; 0 .. v._nrow) {
            if (v._data[j * v._stride + i] != 0.0) allZero = false;
        }
    }
    if (allZero) {                  // the library cannot handle a matrix without nonzero singular values
        auto z = newMatrix(v._ncol, v._nrow);

        z.ptr.val[0 .. v._nrow * v._ncol] = 0.0;
        return z;
    }
    auto c = v.dup;                 // overwritten by the pseudoinverse (and turned into ncol x nrow)
    int err = gretl_matrix_moore_penrose(c.ptr, tol);

    enforce!GretlException(err == 0, "pinv: the computation failed, error code " ~ to!string(err));
    return c;
}

/* det(m), m.det: the determinant of a square matrix (closed form for 2 x 2 and 3 x 3, otherwise LU, as the library does it).
   The determinant of a 0 x 0 matrix is 1. A non-square matrix throws GretlException. */
double det(M)(M m) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "det: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (v._nrow == 0) return 1.0;
    auto c = v.dup;                 // the library overwrites the matrix with its LU factors
    int err = 0;
    const d = gretl_matrix_determinant(c.ptr, &err);

    enforce!GretlException(err == 0, "det: the computation failed, error code " ~ to!string(err));
    return d;
}

/* logdet(m), m.logdet: the natural log of the determinant of a square matrix, from one LU factorization (LAPACK dgetrf, summing the
   logs of the pivots, so it stays finite when the determinant itself would overflow or underflow). A determinant of zero gives
   -infinity and a negative determinant gives NaN (there is no real log); `logdet(m, true)` returns the log of the ABSOLUTE value of
   the determinant instead, which is finite for any nonsingular matrix. (The library's own log-determinant returns its NA value
   with a message on stderr in both of those cases; this one reports them as plain doubles.) The log-determinant of a 0 x 0 matrix is 0. */
double logdet(M)(M m, bool absolute = false) @trusted if (isMatrixLike!M) {
    import std.math: log, fabs;
    import gretl.clapack_doubleh: dgetrf_;
    import gretl.gretl_f2ch: integer;

    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "logdet: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (n == 0) return 0.0;
    enforce!GretlException(n <= int.max, "logdet: the matrix is too large");
    auto a = v.dup;
    auto ipiv = new integer[n];
    integer nn = cast(integer) n, lda = cast(integer) n, info = 0;

    dgetrf_(&nn, &nn, a.ptr.val, &lda, ipiv.ptr, &info);
    enforce!GretlException(info >= 0, "logdet: dgetrf failed, info = " ~ to!string(info));
    if (info > 0) return -double.infinity;          // an exactly zero pivot: singular
    double sum = 0.0;
    bool negative = false;

    foreach (i; 0 .. n) {
        const di = a.ptr.val[i * n + i];

        if (ipiv[i] != cast(integer) i + 1) negative = !negative;
        if (di < 0) negative = !negative;
        sum += log(fabs(di));
    }
    if (negative && !absolute) return double.nan;
    return sum;
}

/* trace(m), m.trace: the sum of the diagonal of a square matrix; a non-square matrix throws GretlException. */
double trace(M)(M m) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "trace: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    double sum = 0.0;

    foreach (i; 0 .. v._nrow) sum += v._data[i * v._stride + i];
    return sum;
}

/* expm(m), m.expm: the matrix exponential, exp(m) = I + m + m^2/2! + ..., of a square matrix, as a new GretlMatrix (gretl_matrix_exp:
   scaling and squaring with a Pade approximation, Golub and Van Loan 11.3.1). It is called expm and not exp, to stay clear of
   std.math.exp, which takes a number; for a diagonal matrix it is the exponential of each diagonal element, and for any m
   det(expm(m)) = exp(trace(m)). A non-square matrix throws GretlException. */
GretlMatrix expm(M)(M m) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(v._nrow == v._ncol, "expm: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (v._nrow == 0) return newMatrix(0, 0);
    auto c = v.dup;                 // gretl_matrix_exp takes a gretl matrix, so a view has to be copied first
    int err = 0;
    auto raw = gretl_matrix_exp(c.ptr, &err);

    enforce!GretlException(err == 0 && raw !is null, "expm: the computation failed, error code " ~ to!string(err));
    return GretlMatrix(raw);
}

/* rank(m), m.rank: the numerical rank, the number of singular values larger than tol * max(nrow, ncol) * (largest singular
   value); the default tol = 0 means machine epsilon (2.2e-16), as in numpy and Matlab. Counted from svd(m, false), not with
   the library's gretl_matrix_rank, which for a matrix with more than 4 times as many rows as columns (or the reverse) takes
   the SVD of m'm and so applies the tolerance to squared singular values (bugs.md #111). */
size_t rank(M)(M m, double tol = 0.0) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    if (v._nrow * v._ncol == 0) return 0;
    auto s = svd(m, false).s;
    const dmax = (v._nrow > v._ncol) ? v._nrow : v._ncol;
    const eps = (tol > 0.0) ? tol : 2.20e-16;
    const smin = dmax * eps * s[0, 0];
    size_t r = 0;

    foreach (i; 0 .. s.nrow) {
        if (s[i, 0] > smin) r++;
    }
    return r;
}

/* eye(n): the n x n identity matrix; eye(nr, nc): an nr x nc matrix with ones on the diagonal and zeros elsewhere. */
GretlMatrix eye(size_t n) {
    return eye(n, n);
}

GretlMatrix eye(size_t nr, size_t nc) @trusted {
    auto e = newMatrix(nr, nc);

    if (nr * nc == 0) return e;
    e.ptr.val[0 .. nr * nc] = 0.0;
    const k = (nr < nc) ? nr : nc;

    foreach (i; 0 .. k) e.ptr.val[i * nr + i] = 1.0;
    return e;
}

/* Which matrix norm. */
enum Norm {
    frobenius,      /* the square root of the sum of squares of all elements */
    one,            /* the largest column sum of absolute values */
    infinity,       /* the largest row sum of absolute values */
    spectral,       /* the 2-norm: the largest singular value */
    max             /* the largest absolute value of an element (not a submultiplicative norm) */
}

/* norm(m, kind), m.norm(kind): a norm of a matrix (or of a vector, an n x 1 or 1 x n matrix, for which frobenius and spectral are
   the usual Euclidean length). The default is Norm.frobenius, which is cheap; Norm.spectral is what Matlab's norm(A) means but needs an
   SVD. frobenius scales by the largest element so that it does not overflow or underflow in the squares. An empty matrix has norm 0,
   a matrix with a NaN has norm NaN. */
double norm(M)(M m, Norm kind = Norm.frobenius) @trusted if (isMatrixLike!M) {
    import std.math: fabs, isNaN, sqrt;

    auto v = wholeView(m);

    if (v._nrow * v._ncol == 0) return 0.0;
    if (kind == Norm.spectral) {
        foreach (j; 0 .. v._ncol) {
            foreach (i; 0 .. v._nrow) {
                if (isNaN(v._data[j * v._stride + i])) return double.nan;
            }
        }
        return svd(m, false).s[0, 0];
    }
    double result = 0.0;

    final switch (kind) {
        case Norm.one:
            foreach (j; 0 .. v._ncol) {
                double colsum = 0.0;

                foreach (i; 0 .. v._nrow) colsum += fabs(v._data[j * v._stride + i]);
                if (isNaN(colsum)) return double.nan;
                if (colsum > result) result = colsum;
            }
            break;
        case Norm.infinity:
            auto rowsum = new double[v._nrow];

            rowsum[] = 0.0;
            foreach (j; 0 .. v._ncol) {
                foreach (i; 0 .. v._nrow) rowsum[i] += fabs(v._data[j * v._stride + i]);
            }
            foreach (x; rowsum) {
                if (isNaN(x)) return double.nan;
                if (x > result) result = x;
            }
            break;
        case Norm.max, Norm.frobenius:
            foreach (j; 0 .. v._ncol) {
                foreach (i; 0 .. v._nrow) {
                    const a = fabs(v._data[j * v._stride + i]);

                    if (isNaN(a)) return double.nan;
                    if (a > result) result = a;
                }
            }
            if (kind == Norm.frobenius && result > 0.0 && result < double.infinity) {
                double sumsq = 0.0;

                foreach (j; 0 .. v._ncol) {
                    foreach (i; 0 .. v._nrow) {
                        const q = v._data[j * v._stride + i] / result;

                        sumsq += q * q;
                    }
                }
                result *= sqrt(sumsq);
            }
            break;
        case Norm.spectral:
            break;      // handled above
    }
    return result;
}

/* cond(m, kind), m.cond(kind): the condition number of a matrix in a norm, norm(m) * norm(inv(m)); the default is Norm.spectral, the ratio of
   the largest to the smallest singular value (so, as in Matlab, it works for a non-square matrix too, with min(nrow, ncol) singular values).
   Norm.one, Norm.infinity and Norm.frobenius are computed as norm(m, kind) * norm(inv(m), kind) and need a square matrix. A singular
   matrix has condition number infinity; an empty one 0. Norm.max is not a condition number and throws GretlException, as does a
   non-square matrix for the other norms. A value near 1 / 2.2e-16 means m is singular to machine precision. */
double cond(M)(M m, Norm kind = Norm.spectral) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);

    enforce!GretlException(kind != Norm.max, "cond: Norm.max does not define a condition number");
    if (v._nrow * v._ncol == 0) return 0.0;
    if (kind == Norm.spectral) {
        auto s = svd(m, false).s;
        const smin = s[s.nrow - 1, 0];

        return (smin == 0.0) ? double.infinity : s[0, 0] / smin;
    }
    enforce!GretlException(v._nrow == v._ncol, "cond: the matrix is " ~ to!string(v._nrow) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (lu(m).singular) return double.infinity;
    return norm(m, kind) * norm(inv(m), kind);
}

private bool isSymmetricView(Submatrix v) @trusted {
    import std.math: fabs;

    double big = 0;

    foreach (j; 0 .. v._ncol) {
        foreach (i; 0 .. v._nrow) big = (fabs(v._data[j * v._stride + i]) > big) ? fabs(v._data[j * v._stride + i]) : big;
    }
    foreach (j; 0 .. v._ncol) {
        foreach (i; j + 1 .. v._nrow) {
            if (fabs(v._data[j * v._stride + i] - v._data[i * v._stride + j]) > 1e-12 * big) return false;
        }
    }
    return true;
}

/* sqrtm(m), m.sqrtm: the principal matrix square root, the matrix s with s * s = m whose eigenvalues all have a non-negative real part,
   as a new GretlMatrix. Written sqrtm, like expm, so that it cannot be mistaken for std.math.sqrt of the elements (that is `m.el ^^ 0.5`).
     - A symmetric matrix (checked to 1e-12 relative) must be positive semidefinite, and the root is the symmetric positive semidefinite
       one, v * diag(sqrt(values)) * t(v) from eigen(m); rounding errors that make an eigenvalue negative by no more than n * 2.2e-16 times
       the largest one are treated as zero. A singular symmetric matrix is fine.
     - Any other square matrix uses the Denman-Beavers iteration (Y <- (Y + inv(Z)) / 2, Z <- (Z + inv(Y)) / 2, which needs inverses, so
       m must be nonsingular), converging quadratically when m has no eigenvalues on the closed negative real axis, in which case the root is
       real. Otherwise there is no real principal root, the iteration does not converge, and GretlException is thrown.
   A non-square matrix throws GretlException; a 0 x 0 matrix gives a 0 x 0 matrix. */
GretlMatrix sqrtm(M)(M m) @trusted if (isMatrixLike!M) {
    import std.math: fabs, sqrt;

    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "sqrtm: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (n == 0) return newMatrix(0, 0);
    if (isSymmetricView(v)) {
        auto e = eigen(m);
        double amax = 0;

        foreach (i; 0 .. n) amax = (fabs(e.values[i, 0]) > amax) ? fabs(e.values[i, 0]) : amax;
        const slack = n * 2.2e-16 * amax;

        enforce!GretlException(e.values[0, 0] >= -slack, "sqrtm: the symmetric matrix is not positive semidefinite (smallest eigenvalue "
                               ~ to!string(e.values[0, 0]) ~ ")");
        auto w = e.vectors;                         // becomes v * diag(lambda^(1/4)), and then root = w * w'
        foreach (j; 0 .. n) {
            const lam = e.values[j, 0];
            const q = (lam > 0) ? sqrt(sqrt(lam)) : 0.0;

            foreach (i; 0 .. n) w[i, j] = w[i, j] * q;
        }
        return w * t(w);
    }
    auto y = v.dup;
    auto z = eye(n);
    const scale = norm(y);

    foreach (iter; 0 .. 100) {
        GretlMatrix yi, zi;

        try {
            yi = inv(y);
            zi = inv(z);
        } catch (GretlException ex) {
            throw new GretlException("sqrtm: the iteration needs inverses and met a singular matrix (is m singular?)");
        }
        auto ynew = (y + zi) * 0.5;
        auto znew = (z + yi) * 0.5;
        const change = norm(ynew - y) / (scale > 0 ? scale : 1.0);

        y = ynew;
        z = znew;
        if (!(change == change) || change > 1e100) break;           // NaN or runaway: it is diverging
        if (change < 1e-15) {
            auto check = y * y - v.dup;
            enforce!GretlException(norm(check) <= 1e-8 * scale, "sqrtm: no real principal square root (the iteration converged to a wrong matrix)");
            return y;
        }
    }
    throw new GretlException("sqrtm: the iteration did not converge, there is no real principal square root (does m have a negative real eigenvalue?)");
}

/* logm(m), m.logm: the principal matrix logarithm, the matrix l with expm(l) = m whose eigenvalues have imaginary parts between -pi and pi,
   as a new GretlMatrix. Written logm, like expm, so that it is not mistaken for the elementwise log.
     - A symmetric matrix (checked to 1e-12 relative) must be positive definite, and the logarithm is v * diag(log(values)) * t(v) from eigen(m).
     - Any other square matrix uses inverse scaling and squaring: take square roots (sqrtm) until m is within 0.25 of the identity in the one
       norm, evaluate log(I + x) there with an 8-point Gauss-Legendre quadrature of the integral of x * inv(I + t x) over t in [0, 1] (which is
       the diagonal Pade approximant of order 8), and multiply the result by 2^(number of roots). This needs m nonsingular with no eigenvalue on
       the closed negative real axis, which is what a real principal logarithm requires; otherwise GretlException is thrown.
   A non-square matrix throws GretlException; a 0 x 0 matrix gives a 0 x 0 matrix. */
GretlMatrix logm(M)(M m) @trusted if (isMatrixLike!M) {
    import std.math: log;

    auto v = wholeView(m);
    const n = v._nrow;

    enforce!GretlException(n == v._ncol, "logm: the matrix is " ~ to!string(n) ~ " x " ~ to!string(v._ncol) ~ ", not square");
    if (n == 0) return newMatrix(0, 0);
    if (isSymmetricView(v)) {
        auto e = eigen(m);

        enforce!GretlException(e.values[0, 0] > 0.0, "logm: the symmetric matrix is not positive definite (smallest eigenvalue "
                               ~ to!string(e.values[0, 0]) ~ "), so there is no real logarithm");
        auto w = e.vectors;                         // v * diag(log(values)) * t(v), scaling a copy of the columns
        auto vl = e.vectors.dup;

        foreach (j; 0 .. n) {
            const lg = log(e.values[j, 0]);

            foreach (i; 0 .. n) vl[i, j] = vl[i, j] * lg;
        }
        return vl * t(w);
    }
    /* Gauss-Legendre nodes and weights on [0, 1] (8 points: nodes (1 + xi) / 2, weights w / 2 for the nodes xi, w on [-1, 1]) */
    static immutable double[4] xi = [0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363];
    static immutable double[4] wi = [0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763];

    auto a = v.dup;
    auto id = eye(n);
    int roots = 0;

    while (norm(a - id, Norm.one) > 0.25) {
        enforce!GretlException(roots < 60, "logm: too many square roots, the matrix is too far from the identity");
        try {
            a = sqrtm(a);
        } catch (GretlException ex) {
            throw new GretlException("logm: no real principal logarithm (" ~ ex.msg ~ ")");
        }
        roots++;
    }
    auto x = a - id;
    auto sum = x * 0.0;

    foreach (j; 0 .. 4) {
        foreach (sign; [-1.0, 1.0]) {
            const node = 0.5 * (1.0 + sign * xi[j]);
            const weight = 0.5 * wi[j];

            sum = sum + (inv(id + x * node) * x) * weight;
        }
    }
    return sum * (2.0 ^^ roots);
}

/* order(m, col, ...), m.order(col, ...): the row indices (0-based) that put the rows of m in ascending order of column `col`; where
   rows tie on it the next column given decides, and so on. Rows that tie on every column given keep their original order (the sort is
   stable, std.algorithm.sorting.sort with SwapStrategy.stable). NaN sorts after every number. It is like R's order(): m.order(0)[0] is the
   row with the smallest first column, and sort(m, ...) is m with its rows taken in this order. A column outside the matrix throws
   GretlException. */
size_t[] order(M)(M m, size_t first, size_t[] rest...) @trusted if (isMatrixLike!M) {
    import std.algorithm.mutation: SwapStrategy;
    import std.algorithm.sorting: stdSort = sort;
    import std.array: array;
    import std.math: isNaN;
    import std.range: iota;

    auto v = wholeView(m);
    auto cols = [first] ~ rest;

    foreach (c; cols) {
        enforce!GretlException(c < v._ncol, "order: column " ~ to!string(c) ~ " is out of range for a matrix with " ~ to!string(v._ncol) ~ " columns");
    }
    auto idx = iota(v._nrow).array;
    static bool lessThan(double a, double b) {
        return isNaN(a) ? false : (isNaN(b) ? true : a < b);
    }
    bool before(size_t i, size_t j) {
        foreach (c; cols) {
            const x = v._data[c * v._stride + i], y = v._data[c * v._stride + j];

            if (lessThan(x, y)) return true;
            if (lessThan(y, x)) return false;
        }
        return false;
    }
    stdSort!(before, SwapStrategy.stable)(idx);
    return idx;
}

/* sort(m, col, ...), m.sort(col, ...): a new GretlMatrix with the rows of m in ascending order of column `col`, ties broken by the next
   column given, and so on (see order; stable, NaN last). The matrix itself is not changed. For example m.sort(1, 0) sorts by the second
   column and, where its values are equal, by the first. */
GretlMatrix sort(M)(M m, size_t first, size_t[] rest...) @trusted if (isMatrixLike!M) {
    auto v = wholeView(m);
    auto idx = order(m, first, rest);
    auto result = newMatrix(v._nrow, v._ncol);

    if (v._nrow * v._ncol == 0) return result;
    foreach (j; 0 .. v._ncol) {
        foreach (r; 0 .. v._nrow) {
            result.ptr.val[j * v._nrow + r] = v._data[j * v._stride + idx[r]];
        }
    }
    return result;
}

/* The direction of a cumulative operation. */
enum Along {
    column,     /* down each column: element [i, j] accumulates [0, j] .. [i, j] (Matlab's dimension 1, gretl's cum()) */
    row         /* across each row:  element [i, j] accumulates [i, 0] .. [i, j] (Matlab's dimension 2) */
}

private GretlMatrix cumulate(string op)(Submatrix v, Along dir) @trusted {
    auto c = newMatrix(v._nrow, v._ncol);

    if (v._nrow * v._ncol == 0) return c;
    double * z = c.ptr.val;
    const nr = v._nrow;

    if (dir == Along.column) {
        foreach (j; 0 .. v._ncol) {
            double acc = v._data[j * v._stride];

            z[j * nr] = acc;
            foreach (i; 1 .. nr) {
                acc = mixin("acc " ~ op ~ " v._data[j * v._stride + i]");
                z[j * nr + i] = acc;
            }
        }
    } else {
        foreach (i; 0 .. nr) z[i] = v._data[i];
        foreach (j; 1 .. v._ncol) {
            foreach (i; 0 .. nr) {
                z[j * nr + i] = mixin("z[(j - 1) * nr + i] " ~ op ~ " v._data[j * v._stride + i]");
            }
        }
    }
    return c;
}

/* cumsum(m, dir), m.cumsum(dir): the cumulative sums, a new matrix of the same size. The default Along.column accumulates down each column,
   so a column vector is its running total; Along.row accumulates across each row. (An n x 1 vector is the same in either direction
   as a 1 x n one is the other way round.) A NaN makes the rest of its column or row NaN. An empty matrix gives an empty matrix. */
GretlMatrix cumsum(M)(M m, Along dir = Along.column) if (isMatrixLike!M) {
    return cumulate!"+"(wholeView(m), dir);
}

/* cumprod(m, dir), m.cumprod(dir): the cumulative products, as cumsum but multiplying. */
GretlMatrix cumprod(M)(M m, Along dir = Along.column) if (isMatrixLike!M) {
    return cumulate!"*"(wholeView(m), dir);
}

/* ---- rbind, cbind, ~ and ~= ---- */

enum isNumericArray(T) = isArray!T && !is(T == string) && isNumeric!(ElementType!T);
enum isOperand(T) = isMatrixLike!T || isNumeric!T || isNumericArray!T;
enum isCbindOperand(T) = isOperand!T;

/* what is stitched together: a matrix (view), a vector (an array), or a scalar */
private struct Piece {
    enum Kind { matrix, vec, scalar }
    Kind kind;
    Submatrix view;
    const(double)[] vec;
    double x;

    /* a 0 x 0 matrix (such as GretlMatrix.init) is left out, so that  acc = rbind(acc, row)  can start from an empty one */
    bool skipped() { return kind == Kind.matrix && view._nrow == 0 && view._ncol == 0; }
    size_t cols() { return skipped ? 0 : (kind == Kind.matrix ? view._ncol : 1); }     // as a piece of cbind
    size_t rows() { return skipped ? 0 : (kind == Kind.matrix ? view._nrow : 1); }     // as a piece of rbind
}

private Piece makePiece(T)(T a) {
    Piece p;

    static if (isMatrixLike!T) {
        p.kind = Piece.Kind.matrix;
        p.view = wholeView(a);
    } else static if (isNumeric!T) {
        p.kind = Piece.Kind.scalar;
        p.x = cast(double) a;
    } else static if (isNumericArray!T) {
        p.kind = Piece.Kind.vec;
        static if (is(T : const(double)[])) {
            p.vec = a;
        } else {
            double[] conv;

            foreach (x; a) conv ~= cast(double) x;
            p.vec = conv;
        }
    } else {
        static assert(false, "cbind and rbind take matrices, scalars and arrays of numbers");
    }
    return p;
}

/* the common number of rows of the pieces of a cbind (the length of vectors, the nrow of matrices; scalars fit anything); `hint` is a row
   count that is already fixed if `hinted` */
private size_t cbindRows(Piece[] ps, size_t hint, bool hinted) {
    size_t nr = hint;
    bool known = hinted;

    foreach (p; ps) {
        if (p.skipped || p.kind == Piece.Kind.scalar) continue;
        const len = (p.kind == Piece.Kind.matrix) ? p.view._nrow : p.vec.length;

        enforce!GretlException(!known || len == nr, "cbind: " ~ (p.kind == Piece.Kind.matrix ? "a matrix with " : "a vector of length ")
                               ~ to!string(len) ~ (p.kind == Piece.Kind.matrix ? " rows" : "") ~ " does not fit " ~ to!string(nr) ~ " rows");
        nr = len;
        known = true;
    }
    return known ? nr : 1;
}

private size_t rbindCols(Piece[] ps) {
    size_t nc = 0;
    bool known = false;

    foreach (p; ps) {
        if (p.skipped || p.kind == Piece.Kind.scalar) continue;
        const len = (p.kind == Piece.Kind.matrix) ? p.view._ncol : p.vec.length;

        enforce!GretlException(!known || len == nc, "rbind: " ~ (p.kind == Piece.Kind.matrix ? "a matrix with " : "a vector of length ")
                               ~ to!string(len) ~ (p.kind == Piece.Kind.matrix ? " columns" : "") ~ " does not fit " ~ to!string(nc) ~ " columns");
        nc = len;
        known = true;
    }
    return known ? nc : 1;
}

/* write the columns of a piece (nr rows each) one after another at z */
private void fillColumns(double * z, size_t nr, Piece p) @trusted {
    final switch (p.kind) {
        case Piece.Kind.matrix:
            foreach (j; 0 .. p.view._ncol) {
                foreach (i; 0 .. nr) z[j * nr + i] = p.view._data[j * p.view._stride + i];
            }
            break;
        case Piece.Kind.vec:
            foreach (i; 0 .. nr) z[i] = p.vec[i];
            break;
        case Piece.Kind.scalar:
            foreach (i; 0 .. nr) z[i] = p.x;
            break;
    }
}

private GretlMatrix cbindPieces(Piece[] ps) @trusted {
    size_t total = 0;

    bool allSkipped = true;

    foreach (p; ps) {
        total += p.cols;
        if (!p.skipped) allSkipped = false;
    }
    if (allSkipped) return newMatrix(0, 0);
    const nr = cbindRows(ps, 0, false);
    auto result = newMatrix(nr, total);
    size_t col = 0;

    if (nr * total == 0) return result;
    foreach (p; ps) {
        if (p.skipped) continue;
        fillColumns(result.ptr.val + col * nr, nr, p);
        col += p.cols;
    }
    return result;
}

private GretlMatrix rbindPieces(Piece[] ps) @trusted {
    size_t total = 0;

    bool allSkipped = true;

    foreach (p; ps) {
        total += p.rows;
        if (!p.skipped) allSkipped = false;
    }
    if (allSkipped) return newMatrix(0, 0);
    const nc = rbindCols(ps);
    auto result = newMatrix(total, nc);
    size_t row = 0;

    if (total * nc == 0) return result;
    foreach (p; ps) {
        if (p.skipped) continue;
        const pr = p.rows;

        foreach (j; 0 .. nc) {
            foreach (i; 0 .. pr) {
                double x;

                final switch (p.kind) {
                    case Piece.Kind.matrix: x = p.view._data[j * p.view._stride + i]; break;
                    case Piece.Kind.vec: x = p.vec[j]; break;
                    case Piece.Kind.scalar: x = p.x; break;
                }
                result.ptr.val[j * total + row + i] = x;
            }
        }
        row += pr;
    }
    return result;
}

/* cbind(a, b, ...): the arguments side by side, as in R: a new GretlMatrix with the columns of a, then those of b, and so on. Each argument is a
   GretlMatrix or Submatrix (all with the same number of rows), a number (recycled to fill a column of that height), or an array of numbers
   (one column, whose length must equal the number of rows). cbind(1, m) adds a column of ones. With only numbers the result is 1 row. A 0 x 0 matrix
   is left out, so a result can be started from GretlMatrix.init. m ~ x is cbind(m, x). */
GretlMatrix cbind(Args...)(Args args) if (Args.length >= 1 && allSatisfy!(isOperand, Args)) {
    Piece[] ps;

    static foreach (a; args) ps ~= makePiece(a);
    return cbindPieces(ps);
}

/* rbind(a, b, ...): the arguments one below the other, as in R: a new GretlMatrix with the rows of a, then those of b, and so on. Each argument is a
   GretlMatrix or Submatrix (all with the same number of columns), a number (recycled to fill a row of that width), or an array of numbers (one row,
   whose length must equal the number of columns). A 0 x 0 matrix is left out. */
GretlMatrix rbind(Args...)(Args args) if (Args.length >= 1 && allSatisfy!(isOperand, Args)) {
    Piece[] ps;

    static foreach (a; args) ps ~= makePiece(a);
    return rbindPieces(ps);
}
