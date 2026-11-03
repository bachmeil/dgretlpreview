/* Port of lib/src/usermat.c: user-matrix lookup by name, matrix sub-specifications (selection by index, range,
   matrix, "diag" and so on: checking, offsets, slices, assignment of scalars and matrices to submatrices,
   extraction of submatrices), row and column names of user matrices, in-place matrix operations (invert,
   Cholesky, transpose, X'X), vec/vech/unvech, QR and SVD decompositions, OLS/RLS on matrices, the GHK
   simulator, symmetric eigenproblems and gretl_matrix_set_part. The header usermat.h is usermath.d.

   The mspec convenience macros of the C are the small functions at the top of this module; they are templates
   over the pointer type of their argument because the C uses them on both specs and selectors. */

module gretl.usermat;


import core.stdc.math : floor, fabs, isnan, sqrt;
import core.stdc.stdio : FILE, fprintf, fputs, printf, sprintf, stderr;
import core.stdc.stdlib : calloc, free, malloc, realloc;
import core.stdc.string : memcpy, memset, strcmp, strlen;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.usermath;

/* mspec convenience macros */

private int mspec_get_offset(S)(S m) { return m.lsel.range[0]; }
private int mspec_get_n_elem(S)(S m) { return m.lsel.range[1]; }

private void mspec_set_offset(S)(S m, int o) { m.lsel.range[0] = o; }
private void mspec_set_n_elem(S)(S m, int n) { m.lsel.range[1] = n; }

private void mspec_set_element(S)(S m, int i) { m.lsel.range[0] = i; }

private bool spec_is_single(S)(S s) { return s.ltype == SelType.single || s.rtype == SelType.single; }
private int spec_single_val(S)(S s) { return s.ltype == SelType.single ? s.lsel.range[0] : s.rsel.range[0]; }

private bool singleton_left_range(S)(S s)
{
    return s.ltype == SelType.range && s.lsel.range[0] == s.lsel.range[1];
}

private bool singleton_right_range(S)(S s)
{
    return s.rtype == SelType.range && s.rsel.range[0] == s.rsel.range[1];
}

private bool all_or_null(SelType t) { return t == SelType.all || t == SelType.null_; }

private bool lhs_is_scalar(S, M)(S s, M m)
{
    return s.ltype == SelType.element || (m.rows == 1 && s.ltype == SelType.all) || singleton_left_range(s);
}

private bool rhs_is_scalar(S, M)(S s, M m)
{
    return s.rtype == SelType.element || (m.cols == 1 && all_or_null(s.rtype)) || singleton_right_range(s);
}

private int rowmax(S, M)(S s, M m) { return s.range[1] == MSEL_MAX ? m.rows : s.range[1]; }
private int colmax(S, M)(S s, M m) { return s.range[1] == MSEL_MAX ? m.cols : s.range[1]; }

/* end mspec convenience macros */

import gretl.gretl_array : gretl_array_get_length, gretl_array_get_strings;
import gretl.gretl_cmatrix : gretl_cmatrix_SVD, gretl_cmatrix_extract, gretl_cmatrix_from_scalar, gretl_cmatrix_set_diagonal, gretl_cmatrix_set_part, gretl_cmatrix_set_triangle, gretl_zheev;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_intl : gretl_is_ascii;
import gretl.gretl_list : gretl_list_new, gretl_null_list;
import gretl.gretl_normal : gretl_GHK2;
import gretl.gretl_prn : Prn, pprintf;
import gretl.matrix_extra : gretl_matrix_from_scalar, gretl_matrix_mp_ols;
import gretl.strutils : gretl_strdup, gretl_string_split, gretl_strndup, string_is_blank, strings_array_dup, strings_array_free, strings_array_new;
import gretl.uservar : get_user_var_of_type_by_name, user_var_get_value, user_var_steal_value;
import gretl.uservarh : UserVar;





/**
 * get_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

GretlMatrixRaw * get_matrix_by_name(const(char) * name)
{
    GretlMatrixRaw * ret = null;

    if (name != null && *name != '\0') {
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.matrix);
    if (u != null) {
        ret = cast(GretlMatrixRaw *) user_var_get_value(u);
    }
    }

    return ret;
}

/**
 * steal_matrix_by_name:
 * @name: name of the matrix.
 *
 * Looks up a user-defined matrix by name and if found,
 * grabs the matrix, leaving the matrix pointer on the
 * named matrix as %NULL.
 *
 * Returns: pointer to matrix, or %NULL if not found.
 */

GretlMatrixRaw * steal_matrix_by_name(const(char) * name)
{
    GretlMatrixRaw * ret = null;

    if (name != null && *name != '\0') {
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.matrix);

    if (u != null) {
        ret = cast(GretlMatrixRaw *) user_var_steal_value(u);
    }
    }

    return ret;
}

/**
 * get_matrix_copy_by_name:
 * @name: name of the matrix.
 * @err: location to receive error code;
 *
 * Looks up a user-defined matrix by name.
 *
 * Returns: a copy of the named matrix, or %NULL if not found.
 */

GretlMatrixRaw * get_matrix_copy_by_name(const(char) * name, int * err)
{
    GretlMatrixRaw * m = get_matrix_by_name(name);
    GretlMatrixRaw * ret = null;

    if (m == null) {
    *err = GretlError.unkvar;
    } else {
    ret = gretl_matrix_copy(m);
    if (ret == null) {
        *err = GretlError.alloc;
    }
    }

    return ret;
}

/* @s is the user-supplied selection vector.
   @n is the max possible value of row or col.
   @pslice is location to receive positive selection
   list.
*/

int handle_vector_exclusion(const(GretlVector) * s, int n, int ** pslice)
{
    int slen = gretl_vector_get_length(s);
    int i, j, k, nsel = n;
    int err = 0;

    /* do we need this check here? */
    for (i=0; i<slen; i++) {
    k = cast(int) fabs(s.val[i]);
    if (k > n) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, k);
        return GretlError.bounds;
    }
    }

    /* how many of the @n values are selected? */
    for (i=1; i<=n; i++) {
    for (j=0; j<slen; j++) {
        if (s.val[j] == -i) {
        /* element i is dropped */
        nsel--;
        break;
        }
    }
    }

    if (nsel == 0) {
    *pslice = gretl_null_list();
    } else {
    /* compose positive selection list */
    int * slice = gretl_list_new(nsel);
    int excluded, kk = 1;

    if (slice != null) {
        for (i=1; i<=n; i++) {
        excluded = 0;
        for (j=0; j<slen; j++) {
            if (s.val[j] == -i) {
            excluded = 1;
            break;
            }
        }
        if (!excluded) {
            slice[kk++] = i;
        }
        }
        *pslice = slice;
    }
    }

    if (*pslice == null) {
    err = GretlError.alloc;
    }

    return err;
}

int bad_sel_vector(const(GretlVector) * v, int n)
{
    int i, k, len = gretl_vector_get_length(v);
    int exclude = v.val[0] < 0;

    for (i=0; i<len; i++) {
    k = cast(int) (exclude ? -v.val[i] : v.val[i]);
    if (k < 1 || k > n) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, k);
        return GretlError.bounds;
    }
    }

    return 0;
}

int bad_sel_range(int * range, int n)
{
    int i, k, err = 0;

    for (i=0; i<2; i++) {
    k = range[i];
    if (k != MSEL_MAX && (k < 1 || k > n)) {
        err = GretlError.bounds;
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, k);
        break;
    }
    }

    return err;
}

int bad_sel_single(int * pk, int n)
{
    int err = 0;

    if (*pk != MSEL_MAX && (*pk < 1 || *pk > n)) {
    err = GretlError.bounds;
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, *pk);
    }

    return err;
}

/* convert a matrix subspec component into list of rows
   or columns: n is the maximal dimension of selectable items
*/

int * mspec_make_list(int type, Msel * sel, int n, int * err)
{
    int * slice = null;
    int single_exclude = 0;
    int i, ns = 0;

    if (type == SelType.all || type == SelType.null_) {
    return null;
    }

    if (type == SelType.matrix) {
    if (sel.m == null) {
        gretl_errmsg_set(_("Range is non-positive!").ptr);
        *err = GretlError.data;
    } else if (gretl_vector_get_length(sel.m) == 0) {
        ns = 0;
    } else if (sel.m.val[0] < 0) {
        *err = handle_vector_exclusion(sel.m, n, &slice);
        return slice; /* we're done */
    } else {
        ns = gretl_vector_get_length(sel.m);
    }
    } else if (type == SelType.str) {
        *err = GretlError.types;
        return null;
    } else {
    /* range or single exclusion */
    int sr0 = sel.range[0];
    int sr1 = sel.range[1];

    if (sr1 == MSEL_MAX) {
        sr1 = sel.range[1] = n;
    }
    if (sr0 < 0 && sr1 == sr0) {
        /* excluding a single row or column? */
        sr0 = -sr0;
        if (sr0 > n) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr,
                     sr0);
        *err = GretlError.bounds;
        } else {
        ns = n - 1;
        single_exclude = sr0;
        }
    } else {
        *err = bad_sel_range(sel.range.ptr, n);
        if (!*err) {
        ns = sel.range[1] - sel.range[0] + 1;
        if (ns <= 0) {
            gretl_errmsg_sprintf(_("Range %d to %d is non-positive!").ptr,
                     sel.range[0], sel.range[1]);
            *err = GretlError.invarg;
        }
        }
    }
    }

    if (!*err) {
    slice = gretl_list_new(ns);
    if (slice == null) {
        *err = GretlError.alloc;
    }
    }

    if (!*err) {
    /* compose and check slice */
    if (single_exclude > 0) {
        int k = 1;

        for (i=1; i<=slice[0]; i++) {
        if (i == single_exclude) {
            k++;
        }
        slice[i] = k++;
        }
    } else {
        for (i=0; i<slice[0]; i++) {
        if (type == SelType.matrix) {
            slice[i+1] = cast(int) sel.m.val[i];
        } else {
            slice[i+1] = sel.range[0] + i;
        }
        }
    }

    for (i=1; i<=slice[0] && !*err; i++) {
        if (slice[i] < 1 || slice[i] > n) {
        gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr,
                     slice[i]);
        *err = GretlError.bounds;
        }
    }
    }

    if (*err) {
    free(slice);
    slice = null;
    }

    return slice;
}

int set_element_index(MatrixSubspec * spec, const(GretlMatrixRaw) * m)
{
    int i = spec.lsel.range[0];
    int j = spec.rsel.range[0];
    int k;

    if (spec.ltype == SelType.all && spec.rtype == SelType.all) {
    k = 0;
    } else if (spec.ltype == SelType.element ||
    (spec.ltype == SelType.range && spec.rtype == SelType.range)) {
    k = (j-1) * m.rows + (i-1);
    } else {
    if (spec.ltype == SelType.contig) {
        fprintf(stderr, "CONTIG on left\n".ptr);
    }
    if (spec.rtype == SelType.contig) {
        fprintf(stderr, "CONTIG on right\n".ptr);
    }
    k = i > j ? (i-1) : (j-1);
    }


    /* record single (zero-based) index */
    mspec_set_offset(spec, k);
    mspec_set_n_elem(spec, 1);

    return 0;
}

int spec_check_dimensions(MatrixSubspec * spec, const(GretlMatrixRaw) * m)
{
    int lt = spec.ltype;
    int rt = spec.rtype;
    int err = 0;

    if (m == null) {
    return GretlError.data;
    }

    /* check spec->lsel against rows */
    if (lt == SelType.matrix) {
    err = bad_sel_vector(spec.lsel.m, m.rows);
    } else if (lt == SelType.range || lt == SelType.element) {
    err = bad_sel_range(spec.lsel.range.ptr, m.rows);
    }
    if (!err) {
    /* check spec->rsel against cols */
    if (rt == SelType.matrix) {
        err = bad_sel_vector(spec.rsel.m, m.cols);
    } else if (rt == SelType.range || rt == SelType.element) {
        err = bad_sel_range(spec.rsel.range.ptr, m.cols);
    }
    }

    return err;
}


/* When we come here we've already screened out DIAG
   selection and single-element selection.
*/

int submatrix_contig(MatrixSubspec * spec, const(GretlMatrixRaw) * m)
{
    int contig = 0;

    if (spec.ltype == SelType.matrix || spec.rtype == SelType.matrix ||
    spec.ltype == SelType.excl || spec.rtype == SelType.excl) {
    ; /* cannot treat as contig */
    } else if (m.rows == 1 || m.cols == 1) {
    /* must be contig */
    contig = 1;
    } else if (singleton_right_range(spec)) {
    /* single column selected */
    contig = 1;
    } else if ((spec.rtype == SelType.all || spec.rtype == SelType.range) &&
           spec.ltype == SelType.all) {
    /* range of columns, all rows */
    contig = 1;
    }

    return contig;
}

/* Determine the offset from the start of m->val, and the
   number of elements, in a matrix slice that has been found
   to consist of contiguous data.
*/

int get_offset_nelem(MatrixSubspec * spec, const(GretlMatrixRaw) * m, int * poff, int * pn)
{
    SelType stype = spec.ltype;
    Msel *sel = &spec.lsel;
    int offset = 0, nelem = 0;

    if (m.cols == 1) {
    /* column vector */
    if (stype == SelType.all) {
        nelem = m.rows;
    } else {
        int rmax = rowmax(sel, m);

        offset = sel.range[0] - 1;
        nelem = rmax - sel.range[0] + 1;
    }
    } else if (m.rows == 1) {
    /* row vector */
    if (spec.rtype != SelType.null_) {
        stype = spec.rtype;
        sel = &spec.rsel;
    }
    if (stype == SelType.all) {
        nelem = m.cols;
    } else {
        int cmax = colmax(sel, m);

        offset = sel.range[0] - 1;
        nelem = cmax - sel.range[0] + 1;
    }
    } else if (singleton_right_range(spec)) {
    /* a single matrix column */
    if (stype == SelType.all) {
        nelem = m.rows;
    } else {
        int rmax = rowmax(sel, m);

        offset = sel.range[0] - 1;
        nelem = rmax - sel.range[0] + 1;
    }
    offset += m.rows * (spec.rsel.range[0] - 1);
    } else {
    /* multiple columns, all rows */
    if (spec.rtype == SelType.all) {
        nelem = m.rows * m.cols;
    } else {
        int cmax, nc;

        sel = &spec.rsel;
        cmax = colmax(sel, m);
        nc = cmax - sel.range[0] + 1;
        offset = m.rows * (sel.range[0] - 1);
        nelem = m.rows * nc;
    }
    }

    *poff = offset;
    *pn = nelem;

    return 0;
}

/* to make clear that spec->lsel pertains to columns, move
   it to the right position */

void commute_selectors(MatrixSubspec * spec)
{
    spec.rtype = spec.ltype;
    spec.rsel = spec.lsel;
    memset(&spec.lsel, 0, spec.lsel.sizeof);
    spec.ltype = SelType.all;
}

/* Catch the case of an implicit column or row specification for
   a sub-matrix of an (n x 1) or (1 x m) matrix; also catch the
   error of giving just one row/col spec for a matrix that has
   more than one row and more than one column.
*/

int check_matrix_subspec(MatrixSubspec * spec, const(GretlMatrixRaw) * m)
{
    int veclen, err = 0;

    if (spec.checked) {
    return 0;
    } else if (spec.ltype == SelType.real_ && !m.is_complex) {
    spec.ltype = spec.rtype = SelType.all;
    goto check_done;
    } else if (is_sel_dummy(spec.ltype)) {
    /* SEL_DIAG, SEL_UPPER, SEL_LOWER, SEL_IMAG:
       nothing to be done here?
    */
    goto check_done;
    }


    if (spec.rtype == SelType.null_ && m.rows == 1 && m.cols > 1) {
    /* row vector: transfer spec to column dimension */
    commute_selectors(spec);
    }

    veclen = gretl_vector_get_length(m);

    if (spec.rtype == SelType.null_ && m.rows != 1 && m.cols != 1) {
    gretl_errmsg_set(_("Ambiguous matrix index").ptr);
    return GretlError.data;
    }

    if (spec_is_single(spec)) {
    int k = spec_single_val(spec);

    err = bad_sel_single(&k, veclen);
    if (!err) {
        spec.ltype = spec.rtype = SelType.element;
        mspec_set_offset(spec, k - 1);
        mspec_set_n_elem(spec, 1);
    }
    /* nothing more to do */
    goto check_done;
    }

    err = spec_check_dimensions(spec, m);
    if (err) {
    return err;
    }

    if (m.rows == 0 || m.cols == 0) {
    /* nothing more to do */
    goto check_done;
    }

    if (lhs_is_scalar(spec, m) && rhs_is_scalar(spec, m)) {
    /* we're looking at just one element */
    set_element_index(spec, m);
    spec.ltype = spec.rtype = SelType.element;
    goto check_done;
    }

    if (submatrix_contig(spec, m)) {
    int offset, nelem;

    get_offset_nelem(spec, m, &offset, &nelem);
    spec.ltype = SelType.contig;
    mspec_set_offset(spec, offset);
    mspec_set_n_elem(spec, nelem);
    if (offset < 0 || nelem <= 0) {
        fprintf(stderr, "*** offset = %d, nelem = %d (m: %dx%d) ***\n".ptr,
            offset, nelem, m.rows, m.cols);
            gretl_errmsg_set(_("Invalid submatrix specification").ptr);
        return GretlError.data;
    }
    }

 check_done:

    if (!err) {
    spec.checked = 1;
    }

    return err;
}

const(char) * mspec_get_string(MatrixSubspec * spec, int i)
{
    if (spec == null || i < 0 || i > 1) {
    return null;
    } else if (i == 0) {
    return spec.ltype != SelType.str ? null : spec.lsel.str;
    } else {
    return spec.rtype != SelType.str ? null : spec.rsel.str;
    }
}

int get_slices(MatrixSubspec * spec, const(GretlMatrixRaw) * M)
{
    int err = 0;

    spec.rslice = mspec_make_list(spec.ltype, &spec.lsel,
                   M.rows, &err);
    if (!err) {
    spec.cslice = mspec_make_list(spec.rtype, &spec.rsel,
                       M.cols, &err);
    }


    return err;
}

int assign_scalar_to_submatrix(GretlMatrixRaw * M, const(GretlMatrixRaw) * S, double x, MatrixSubspec * spec)
{
    Cmplx z = Cmplx(NADBL, 0);
    int mr = M.rows;
    int mc = M.cols;
    int i, err = 0;

    if (spec == null) {
    fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n".ptr);
    return GretlError.data;
    }

    if (M.is_complex) {
    z = (S != null)? S.z[0] : Cmplx(x, 0);
    }

    if (spec.ltype == SelType.contig) {
    int ini = mspec_get_offset(spec);
    int fin = ini + mspec_get_n_elem(spec);

    for (i=ini; i<fin; i++) {
        if (M.is_complex) {
        M.z[i] = z;
        } else {
        M.val[i] = x;
        }
    }
    return 0;
    }

    if (is_sel_dummy(spec.ltype)) {
    if (M.is_complex) {
        err = gretl_matrix_set_part(M, S, x, spec.ltype);
    } else {
        err = gretl_matrix_set_part(M, null, x, spec.ltype);
    }
    return err; /* we're done */
    }

    if (spec.rslice == null && spec.cslice == null) {
    /* parse @spec into lists of affected rows and columns */
    err = get_slices(spec, M);
    }

    if (!err) {
    int sr = (spec.rslice == null)? mr : spec.rslice[0];
    int sc = (spec.cslice == null)? mc : spec.cslice[0];
    int j, l, k = 0;
    int mi, mj;

    for (i=0; i<sr; i++) {
        mi = (spec.rslice == null)? k++ : spec.rslice[i+1] - 1;
        l = 0;
        for (j=0; j<sc; j++) {
        mj = (spec.cslice == null)? l++ : spec.cslice[j+1] - 1;
        if (M.is_complex) {
            gretl_cmatrix_set(M, mi, mj, z);
        } else {
            gretl_matrix_set(M, mi, mj, x);
        }
        }
    }
    }

    return err;
}

MatrixSubspec * matrix_subspec_new()
{
    MatrixSubspec * spec = cast(MatrixSubspec *) calloc(1, MatrixSubspec.sizeof);

    if (spec != null) {
    spec.rslice = spec.cslice = null;
    }

    return spec;
}

int contig_cols(MatrixSubspec * spec, const(GretlMatrixRaw) * m)
{
    if (spec.rtype == SelType.all) {
    return m.cols;
    } else if (spec.rtype == SelType.range) {
    int r1 = spec.rsel.range[1];
    int cmax = r1 == MSEL_MAX ? m.cols : r1;

    return cmax - spec.rsel.range[0] + 1;
    } else {
    return 1;
    }
}

/* for 'real' matrices with doubles as elements */

void transcribe_cols_8(GretlMatrixRaw * M, const(GretlMatrixRaw) * S, const(int) * cslice, int sscalar)
{
    int i, j, mcol, nr = M.rows;

    if (sscalar) {
    /* write RHS scalar into all rows of each selected col */
    for (j=1; j<=cslice[0]; j++) {
        mcol = cslice[j] - 1;
        for (i=0; i<nr; i++) {
        gretl_matrix_set(M, i, mcol, S.val[0]);
        }
    }
    } else {
    /* zap from (cols of) @S into selected cols of @M */
    double * xtarg;
    const(double) * xsrc = S.val;
    size_t colsize = nr * (*M.val).sizeof;

    for (j=1; j<=cslice[0]; j++) {
        mcol = cslice[j] - 1;
        xtarg = M.val + mcol * nr;
        memcpy(xtarg, xsrc, colsize);
        xsrc += nr;
    }
    }
}

/* for complex matrices */

void transcribe_cols_16(GretlMatrixRaw * M, const(GretlMatrixRaw) * S, const(int) * cslice, int sscalar)
{
    int i, j, mcol, nr = M.rows;

    if (sscalar) {
    /* write RHS scalar into all rows of each selected col */
    Cmplx z = S.is_complex ? S.z[0] : Cmplx(S.val[0], 0);

    for (j=1; j<=cslice[0]; j++) {
        mcol = cslice[j] - 1;
        for (i=0; i<nr; i++) {
        gretl_cmatrix_set(M, i, mcol, z);
        }
    }
    } else {
    /* zap from (cols of) S into selected cols of M */
    Cmplx * ztarg;
    const(Cmplx) * zsrc = S.z;
    const(double) * xsrc = S.val;
    size_t colsize = nr * (*M.z).sizeof;

    for (j=1; j<=cslice[0]; j++) {
        mcol = cslice[j] - 1;
        ztarg = M.z + mcol * nr;
        if (!S.is_complex) {
        for (i=0; i<nr; i++) {
            ztarg[i] = Cmplx(xsrc[i], 0);
        }
        xsrc += nr;
        } else {
        memcpy(ztarg, zsrc, colsize);
        zsrc += nr;
        }
    }
    }
}

/* @M is the target for partial replacement, @S is the source to
   substitute, and @spec tells how/where to make the
   substitution.
*/

int matrix_replace_submatrix(GretlMatrixRaw * M, const(GretlMatrixRaw) * S, MatrixSubspec * spec)
{
    int mr = gretl_matrix_rows(M);
    int mc = gretl_matrix_cols(M);
    int sr = gretl_matrix_rows(S);
    int sc = gretl_matrix_cols(S);
    int sscalar = 0;
    int i, j;
    int err = 0;

    if (!M.is_complex && S.is_complex) {
    fputs("matrix_replace_submatrix: M is real but S is complex\n", stderr);
    return GretlError.mixed;
    }


    if (spec == null) {
    fprintf(stderr, "matrix_replace_submatrix: spec is NULL!\n".ptr);
    return GretlError.data;
    }

    if (is_sel_dummy(spec.ltype)) {
    return gretl_matrix_set_part(M, S, 0, spec.ltype);
    }

    if (spec.ltype == SelType.contig) {
    int ini = mspec_get_offset(spec);
    int n = mspec_get_n_elem(spec);
    int ccols = contig_cols(spec, M);

    if (M.is_complex && S.rows == 1 && S.cols == 1) {
        for (i=0; i<n; i++) {
        M.z[ini + i] = S.z[0];
        }
    } else if (S.rows * S.cols != n) {
        err = GretlError.nonconf;
    } else if (ccols > 1 && M.rows != S.rows) {
        err = GretlError.nonconf;
    } else if (M.is_complex && !S.is_complex) {
        /* can't use memcpy here! */
        for (i=0; i<n; i++) {
        M.z[ini + i] = Cmplx(S.val[i], 0);
        }
    } else if (M.is_complex) {
        memcpy(M.z + ini, S.z, n * (*M.z).sizeof);
    } else {
        memcpy(M.val + ini, S.val, n * (*M.val).sizeof);
    }
    return err;
    }

    if (sr > mr || sc > mc) {
    /* the replacement matrix won't fit into M */
    fprintf(stderr, "matrix_replace_submatrix: target is %d x %d but replacement part is %d x %d\n".ptr, mr, mc, sr, sc);
    return GretlError.nonconf;
    }

    if (spec.rslice == null && spec.cslice == null) {
    /* parse mspec into lists of affected rows and columns */
    err = get_slices(spec, M);
    if (err) {
        return err;
    }
    }


    if (sr == 1 && sc == 1) {
    /* the replacement is a scalar */
    sscalar = 1;
    sr = (spec.rslice == null)? mr : spec.rslice[0];
    sc = (spec.cslice == null)? mc : spec.cslice[0];
    } else if (spec.rslice != null && spec.rslice[0] != sr) {
    fprintf(stderr, "mspec has %d rows but substitute matrix has %d\n".ptr,
        spec.rslice[0], sr);
    err = GretlError.nonconf;
    } else if (spec.cslice != null && spec.cslice[0] != sc) {
    fprintf(stderr, "mspec has %d cols but substitute matrix has %d\n".ptr,
        spec.cslice[0], sc);
    err = GretlError.nonconf;
    }

    if (!err && spec.rslice == null && spec.cslice != null) {
    /* the target is just specified by column(s) */
    if (M.is_complex) {
        transcribe_cols_16(M, S, spec.cslice, sscalar);
    } else {
        transcribe_cols_8(M, S, spec.cslice, sscalar);
    }
    } else if (!err) {
    /* the general case, no special shortcuts */
    Cmplx z = Cmplx(0, 0);
    double x = 0;
    int l, k = 0;
    int mi, mj;

    if (sscalar && S.is_complex) {
        z = S.z[0];
    } else if (sscalar) {
        x = S.val[0];
        z = Cmplx(x, 0); /* in case M is complex */
    }

    for (j=0; j<sc; j++) {
        mj = (spec.cslice == null)? k++ : spec.cslice[j+1] - 1;
        l = 0;
        for (i=0; i<sr; i++) {
        mi = (spec.rslice == null)? l++ : spec.rslice[i+1] - 1;
        if (!sscalar) {
            if (S.is_complex) {
            z = gretl_cmatrix_get(S, i, j);
            } else {
            x = gretl_matrix_get(S, i, j);
            }
        }
        if (M.is_complex) {
            gretl_cmatrix_set(M, mi, mj, z);
        } else {
            gretl_matrix_set(M, mi, mj, x);
        }
        }
    }
    }

    return err;
}

GretlMatrixRaw * matrix_get_submatrix(const(GretlMatrixRaw) * M, MatrixSubspec * spec, int prechecked, int * err)
{
    GretlMatrixRaw * S = null;
    int r = -1;
    int c = -1;

    if (M == null || spec == null) {
    *err = GretlError.data;
    return null;
    }


    if (!prechecked) {
    *err = check_matrix_subspec(spec, M);
    if (*err) {
        return null;
    }
    }

    if (spec.ltype == SelType.diag) {
    return gretl_matrix_get_diagonal(M, err);
    } else if (spec.ltype == SelType.upper || spec.ltype == SelType.lower) {
    int upper = (spec.ltype == SelType.upper);

    return gretl_matrix_get_triangle(M, upper, err);
    } else if (spec.ltype == SelType.real_ || spec.ltype == SelType.imag) {
    int im = (spec.ltype == SelType.imag);

    return gretl_cmatrix_extract(M, im, err);
    } else if (spec.ltype == SelType.contig) {
    return matrix_get_chunk(M, spec, err);
    } else if (spec.ltype == SelType.element) {
    int i = mspec_get_element(spec);

    if (M.is_complex) {
        S = gretl_cmatrix_from_scalar(M.z[i], err);
    } else {
        S = gretl_matrix_from_scalar(M.val[i]);
    }
    if (S == null) {
        *err = GretlError.alloc;
    }
    return S;
    }

    if (spec.rslice == null && spec.cslice == null) {
    *err = get_slices(spec, M);
    if (*err) {
        return null;
    }
    }


    if (M.rows == 1 && M.cols == 1) {
        if (spec.ltype == SelType.excl && spec.rtype == SelType.null_) {
            r = c = 0;
        } else if (spec.ltype == SelType.excl && spec.rtype == SelType.excl) {
            r = c = 0;
        } else if (spec.ltype == SelType.excl && spec.rtype == SelType.all) {
            r = 0; c = 1;
        } else if (spec.ltype == SelType.all && spec.rtype == SelType.excl) {
            r = 1; c = 0;
        }
    }
    if (r < 0 && c < 0) {
        /* dimensions not yet determined */
        r = (spec.rslice == null)? M.rows : spec.rslice[0];
        c = (spec.cslice == null)? M.cols : spec.cslice[0];
    }

    S = gretl_matching_matrix_new(r, c, M);

    if (S == null) {
    *err = GretlError.alloc;
    } else if (r * c == 0) {
        /* empty matrix: nothing more to be done */
        return S;
    } else {
    int j, mj;

    if (spec.cslice != null && spec.rslice == null) {
        /* copying entire columns */
        if (M.is_complex) {
        Cmplx * dest = S.z;
        size_t csize = r * (*dest).sizeof;

        for (j=0; j<c; j++) {
            mj = spec.cslice[j+1] - 1;
            memcpy(dest, M.z + mj * r, csize);
            dest += r;
        }
        } else {
        double * dest = S.val;
        size_t csize = r * (*dest).sizeof;

        for (j=0; j<c; j++) {
            mj = spec.cslice[j+1] - 1;
            memcpy(dest, M.val + mj * r, csize);
            dest += r;
        }
        }
    } else {
        int i, mi, l, k = 0;
        Cmplx z;
        double x;

        for (j=0; j<c; j++) {
        mj = (spec.cslice == null)? k++ : spec.cslice[j+1] - 1;
        l = 0;
        for (i=0; i<r; i++) {
            mi = (spec.rslice == null)? l++ : spec.rslice[i+1] - 1;
            if (M.is_complex) {
            z = gretl_cmatrix_get(M, mi, mj);
            gretl_cmatrix_set(S, i, j, z);
            } else {
            x = gretl_matrix_get(M, mi, mj);
            gretl_matrix_set(S, i, j, x);
            }
        }
        }
    }
    }

    if (S != null) {
    /* try transcribing metadata on @M if applicable */
    if (S.rows == M.rows && gretl_matrix_is_dated(M)) {
        gretl_matrix_transcribe_obs_info(S, M);
    }
    if (S.cols == M.cols) {
        const(char) ** cnames = gretl_matrix_get_colnames(M);

        if (cnames != null) {
        char ** cpy = strings_array_dup(cast(char **) cnames, M.cols);

        if (cpy != null) {
            gretl_matrix_set_colnames(S, cpy);
        }
        }
    }
    }

    return S;
}

/* Return the element of @M that occupies 0-based position @i in
   its internal vec representation.
*/

double matrix_get_element(const(GretlMatrixRaw) * M, int i, int * err)
{
    double x = NADBL;

    if (M == null) {
    *err = GretlError.data;
    } else if (i < 0 || i >= M.rows * M.cols) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    *err = GretlError.bounds;
    } else {
    x = M.val[i];
    }

    return x;
}

/* Copy a contiguous chunk of data out of @M, the offset and
   extent of which are given by @spec.
*/

GretlMatrixRaw * matrix_get_chunk(const(GretlMatrixRaw) * M, MatrixSubspec * spec, int * err)
{
    int offset = spec.lsel.range[0];
    int nelem = spec.lsel.range[1];
    GretlMatrixRaw * ret;
    int rows, cols;

    if (offset < 0) {
    fprintf(stderr, "matrix_get_chunk: offset = %d\n".ptr, offset);
    *err = GretlError.data;
    return null;
    }

    if (M.cols > 1 && M.rows > 1) {
    cols = contig_cols(spec, M);
    rows = nelem / cols;
    } else if (M.rows == 1) {
    rows = 1;
    cols = nelem;
    } else {
    rows = nelem;
    cols = 1;
    }

    ret = gretl_matching_matrix_new(rows, cols, M);

    if (ret == null) {
    *err = GretlError.alloc;
    } else {
    size_t sz;

    if (M.is_complex) {
        sz = nelem * (*M.z).sizeof;
        memcpy(ret.z, M.z + offset, sz);
    } else {
        sz = nelem * (*M.val).sizeof;
        memcpy(ret.val, M.val + offset, sz);
    }
    if (M.rows > 1 && rows == M.rows && offset % rows == 0 &&
        gretl_matrix_is_dated(M)) {
        gretl_matrix_transcribe_obs_info(ret, M);
    }
    }

    return ret;
}

/* Handle the case where we got a single string as argument
   to colnames() or rownames(), for a matrix with more than
   one column or row: construct specific names by appending
   a column or row index.
*/

char ** expand_names(const(char) * s, int n, int * err)
{
    char ** S = null;

    if (!gretl_is_ascii(s)) {
    *err = GretlError.invarg;
    return null;
    }

    S = strings_array_new(n);

    if (S == null) {
    *err = GretlError.alloc;
    } else {
    char[12] tmp;
    int i, m;

    for (i=0; i<n && !*err; i++) {
        sprintf(tmp.ptr, "%d".ptr, i+1);
        m = cast(int) strlen(tmp.ptr);
        sprintf(tmp.ptr, "%.*s%d".ptr, 9-m, s, i+1);
        S[i] = gretl_strdup(tmp.ptr);
        if (S[i] == null) {
        *err = GretlError.alloc;
        }
    }
    if (*err) {
       strings_array_free(S, n);
       S = null;
    }
    }

    return S;
}

int real_umatrix_set_names(GretlMatrixRaw * M, char ** S, int byrow)
{
    int n = byrow ? M.rows : M.cols;
    int set_done = 0;
    int err = 0;

    if (S != null) {
    int i, n_empty = 0;

    for (i=0; i<n && !err; i++) {
        if (S[i] == null) {
        gretl_errmsg_sprintf("Missing string in %s".ptr,
                     byrow? "rnameset".ptr : "cnameset".ptr);
        err = GretlError.invarg;
        } else if (S[i][0] == '\0') {
                n_empty++;
            }
    }
    if (!err && n_empty < n) {
        if (byrow) {
        gretl_matrix_set_rownames(M, S);
        } else {
        gretl_matrix_set_colnames(M, S);
        }
            set_done = 1;
    }
    }

    /* !set_done implies that either (a) S = NULL or (b) n_empty = n
       or (c) an error was flagged
    */
    if (!set_done) {
        strings_array_free(S, n);
        if (!err && byrow) {
            gretl_matrix_set_rownames(M, null);
        } else if (!err) {
            gretl_matrix_set_colnames(M, null);
        }
    }

    return err;
}

int n_names_check(int n, int ns, int byrow, Prn * prn)
{
    if (ns != n) {
     pprintf(prn, "*** Warning: %s: got %d names but matrix has %d %s ***\n".ptr,
                byrow ? "rnameset".ptr : "cnameset".ptr,
                ns, n, byrow ? "row(s)".ptr : "column(s)".ptr);
    return GretlError.invarg;
    } else {
    return 0;
    }
}

int umatrix_set_names_from_string(GretlMatrixRaw * M, const(char) * s, int byrow, Prn * prn)
{
    int n = byrow ? M.rows : M.cols;
    int err = 0;

    if (s == null || *s == '\0') {
        /* OK: delete any existing names */
        err = real_umatrix_set_names(M, null, byrow);
    } else if (string_is_blank(s)) {
        err = GretlError.invarg;
    } else {
        int ns = 0;
    char ** S = gretl_string_split(s, &ns, " \n\t");

    if (S == null) {
        err = GretlError.alloc;
    } else if (ns == 1 && n > 1) {
        strings_array_free(S, ns);
        S = expand_names(s, n, &err);
    } else {
        err = n_names_check(n, ns, byrow, prn);
        if (err) {
        strings_array_free(S, ns);
        }
    }
        if (!err) {
            err = real_umatrix_set_names(M, S, byrow);
        }
    }

    return err;
}

int umatrix_set_names_from_array(GretlMatrixRaw * M, void * data, int byrow, Prn * prn)
{
    GretlArray * A = cast(GretlArray *) data;
    char ** S = null;
    int n, ns = 0;
    int err = 0;

    n = byrow ? M.rows : M.cols;

    if (A != null && gretl_array_get_length(A) > 0) {
    char ** AS = gretl_array_get_strings(A, &ns);

    err = n_names_check(n, ns, byrow, prn);
    if (!err) {
        S = strings_array_dup(AS, n);
        if (S == null) {
        err = GretlError.alloc;
        }
    }
    }

    if (!err) {
    err = real_umatrix_set_names(M, S, byrow);
    }

    return err;
}

int umatrix_set_names_from_list(GretlMatrixRaw * M, const(int) * list, const(DatasetRaw) * dset, int byrow, Prn * prn)
{
    char ** S = null;
    int n, err = 0;

    n = byrow ? M.rows : M.cols;

    if (list != null && list[0] > 0) {
    int i;

    err = n_names_check(n, list[0], byrow, prn);
    if (!err) {
        S = strings_array_new(n);
        if (S == null) {
        err = GretlError.alloc;
        }
    }
    for (i=0; i<n && !err; i++) {
        S[i] = gretl_strndup(dset.varname[list[i+1]], 12);
        if (S[i] == null) {
        err = GretlError.alloc;
        }
    }
    if (err) {
        strings_array_free(S, n);
    }
    }

    if (!err) {
    err = real_umatrix_set_names(M, S, byrow);
    }

    return err;
}

char * user_matrix_get_column_name(const(GretlMatrixRaw) * M, int col, int * err)
{
    char * ret = null;

    if (M == null || col < 1 || col > M.cols) {
    *err = GretlError.data;
    } else {
    const(char) ** S = gretl_matrix_get_colnames(M);

    if (S == null) {
        ret = gretl_strdup("");
    } else {
        ret = gretl_strdup(S[col-1]);
    }
    if (ret == null) {
        *err = GretlError.alloc;
    }
    }

    return ret;
}

char * user_matrix_get_row_name(const(GretlMatrixRaw) * M, int row, int * err)
{
    char * ret = null;

    if (M == null || row < 1 || row > M.rows) {
    *err = GretlError.data;
    } else {
    const(char) ** S = gretl_matrix_get_rownames(M);

    if (S == null) {
        ret = gretl_strdup("");
    } else {
        ret = gretl_strdup(S[row-1]);
    }
    if (ret == null) {
        *err = GretlError.alloc;
    }
    }

    return ret;
}

double user_matrix_get_determinant(GretlMatrixRaw * m, int tmpmat, int ldet, int * err)
{
    GretlMatrixRaw * R = null;
    double d = NADBL;

    if (gretl_is_null_matrix(m)) {
    return d;
    } else if (tmpmat) {
    /* it's OK to overwrite @m */
    R = m;
    } else {
    /* @m should not be over-written! */
    R = gretl_matrix_copy(m);
    }

    if (R != null) {
    if (ldet) {
        d = gretl_matrix_log_determinant(R, err);
    } else {
        d = gretl_matrix_determinant(R, err);
    }
    if (R != m) {
        gretl_matrix_free(R);
    }
    }

    return d;
}

void maybe_replace_content(GretlMatrixRaw * targ, GretlMatrixRaw * tmp, int err)
{
    if (!err) {
    gretl_matrix_replace_content(targ, tmp);
    }
    gretl_matrix_free(tmp);
}

int matrix_invert_in_place(GretlMatrixRaw * m)
{
    GretlMatrixRaw * tmp = gretl_matrix_copy(m);
    int err = 0;

    if (tmp == null) {
    err = GretlError.alloc;
    } else {
    err = gretl_invert_matrix(tmp);
    maybe_replace_content(m, tmp, err);
    }

    return err;
}

int matrix_cholesky_in_place(GretlMatrixRaw * m)
{
    GretlMatrixRaw * tmp = gretl_matrix_copy(m);
    int err = 0;

    if (tmp == null) {
    err = GretlError.alloc;
    } else {
    err = gretl_matrix_cholesky_decomp(tmp);
    maybe_replace_content(m, tmp, err);
    }

    return err;
}

int matrix_transpose_in_place(GretlMatrixRaw * m)
{
    GretlMatrixRaw * tmp = gretl_matrix_copy_transpose(m);
    int err = 0;

    if (tmp == null) {
    err = GretlError.alloc;
    } else {
    maybe_replace_content(m, tmp, 0);
    }

    return err;
}

int matrix_XTX_in_place(GretlMatrixRaw * m)
{
    GretlMatrixRaw * tmp = gretl_matrix_alloc(m.cols, m.cols);
    int err;

    if (tmp == null) {
    err = GretlError.alloc;
    } else {
    err = gretl_matrix_multiply_mod(m, GretlMatrixMod.transpose,
                    m, GretlMatrixMod.none,
                    tmp, GretlMatrixMod.none);
    maybe_replace_content(m, tmp, err);
    }

    return err;
}

GretlMatrixRaw * user_matrix_vec(const(GretlMatrixRaw) * m, int * err)
{
    GretlMatrixRaw * R = null;

    if (gretl_is_null_matrix(m)) {
    R = gretl_null_matrix_new();
    } else {
    R = gretl_matching_matrix_new(m.rows * m.cols, 1, m);
    if (R != null) {
        gretl_matrix_vectorize(R, m);
    }
    }

    if (R == null) {
    *err = GretlError.alloc;
    }

    return R;
}

GretlMatrixRaw * user_matrix_vech(const(GretlMatrixRaw) * m, int omit_diag, int * err)
{
    GretlMatrixRaw * R = null;

    if (gretl_is_null_matrix(m)) {
    R = gretl_null_matrix_new();
    } else if (m.rows != m.cols) {
    *err = GretlError.nonconf;
    } else {
    int k, n = m.rows;

    if (omit_diag) {
        k = n * (n - 1) / 2;
    } else {
        k = n * (n + 1) / 2;
    }
    R = gretl_matching_matrix_new(k, 1, m);
    if (R != null) {
        if (omit_diag) {
        *err = gretl_matrix_vectorize_h_skip(R, m);
        } else {
        *err = gretl_matrix_vectorize_h(R, m);
        }
    }
    }

    if (R == null && !*err) {
    *err = GretlError.alloc;
    }

    return R;
}

GretlMatrixRaw * user_matrix_unvech(const(GretlMatrixRaw) * m, double diag, int * err)
{
    GretlMatrixRaw * R = null;
    int k;

    if (gretl_is_null_matrix(m)) {
    R = gretl_null_matrix_new();
    } else if ((k = gretl_vector_get_length(m)) == 0) {
    *err = GretlError.nonconf;
    } else {
    int n;

    if (na(diag)) {
        n = cast(int) ((sqrt(1.0 + 8.0 * m.rows) - 1.0) / 2.0);
    } else {
        n = cast(int) ((sqrt(1.0 + 8.0 * m.rows) + 1.0) / 2.0);
    }
    R = gretl_matching_matrix_new(n, n, m);
    if (R != null) {
        if (na(diag)) {
        *err = gretl_matrix_unvectorize_h(R, m);
        } else {
        *err = gretl_matrix_unvectorize_h_diag(R, m, diag);
        }
    }
    }

    if (R == null && !*err) {
    *err = GretlError.alloc;
    }

    return R;
}

int real_user_matrix_QR_decomp(const(GretlMatrixRaw) * m, GretlMatrixRaw ** Q, GretlMatrixRaw ** R, GretlMatrixRaw ** P)
{
    int mc = gretl_matrix_cols(m);
    int err = 0;

    *Q = gretl_matrix_copy(m);

    if (*Q == null) {
    err = GretlError.alloc;
    } else {
    if (R != null) {
        *R = gretl_matrix_alloc(mc, mc);
        if (*R == null) {
        err = GretlError.alloc;
        }
    }
    if (!err && P != null) {
        *P = gretl_matrix_alloc(1, mc);
        if (*P == null) {
        err = GretlError.alloc;
        }
    }
    }

    if (!err) {
    if (P != null) {
        err = gretl_matrix_QR_pivot_decomp(*Q, (R == null)? null : *R, *P);
    } else {
        err = gretl_matrix_QR_decomp(*Q, (R == null)? null : *R);
    }
    }

    if (err) {
    gretl_errmsg_set(_("Matrix decomposition failed").ptr);
    gretl_matrix_free(*Q);
    *Q = null;
    if (R != null) {
        gretl_matrix_free(*R);
        *R = null;
    }
    if (P != null) {
        gretl_matrix_free(*P);
        *P = null;
    }
    }

    return err;
}

GretlMatrixRaw * user_matrix_QR_decomp(const(GretlMatrixRaw) * m, GretlMatrixRaw * R, GretlMatrixRaw * P, int * err)
{
    GretlMatrixRaw * Q = null;
    GretlMatrixRaw * Rtmp = null;
    GretlMatrixRaw * Ptmp = null;
    GretlMatrixRaw ** pR;
    GretlMatrixRaw ** pP;

    if (gretl_is_null_matrix(m)) {
    *err = GretlError.data;
    return null;
    }

    if (m.cols > m.rows) {
    gretl_errmsg_set(_("qrdecomp: the input must have rows >= columns").ptr);
    *err = GretlError.nonconf;
    return null;
    }

    pR = (R != null)? &Rtmp : null;
    pP = (P != null)? &Ptmp : null;

    *err = real_user_matrix_QR_decomp(m, &Q, pR, pP);

    if (Rtmp != null) {
    maybe_replace_content(R, Rtmp, *err);
    }
    if (Ptmp != null) {
    maybe_replace_content(P, Ptmp, *err);
    }

    return Q;
}

GretlMatrixRaw * user_matrix_SVD(const(GretlMatrixRaw) * m, GretlMatrixRaw * U, GretlMatrixRaw * V, int * err)
{
    GretlMatrixRaw * S = null;
    GretlMatrixRaw * Utmp = null;
    GretlMatrixRaw * Vtmp = null;
    GretlMatrixRaw ** pU, pV;

    if (gretl_is_null_matrix(m)) {
    *err = GretlError.data;
    return null;
    }

    pU = (U != null)? &Utmp : null;
    pV = (V != null)? &Vtmp : null;

    if (m.is_complex) {
    *err = gretl_cmatrix_SVD(m, pU, &S, pV, 0);
    } else {
    *err = gretl_matrix_SVD(m, pU, &S, pV, 0);
    }
    if (!*err) {
    if (Utmp != null) {
        maybe_replace_content(U, Utmp, *err);
    }
    if (Vtmp != null) {
        maybe_replace_content(V, Vtmp, *err);
    }
    }

    return S;
}

GretlMatrixRaw * null_OLS(const(GretlMatrixRaw) * Y, GretlMatrixRaw * U, GretlMatrixRaw * V, int * err)
{
    GretlMatrixRaw * B = null;

    if (U != null) {
    /* residuals = Y */
    GretlMatrixRaw * Ycpy;

    Ycpy = gretl_matrix_copy(Y);
    if (Ycpy == null) {
        *err = GretlError.alloc;
    } else {
        gretl_matrix_replace_content(U, Ycpy);
        gretl_matrix_free(Ycpy);
    }
    }

    if (!*err && V != null) {
    /* variance is empty */
    GretlMatrixRaw * nullV;

    nullV = gretl_null_matrix_new();
    gretl_matrix_replace_content(V, nullV);
    gretl_matrix_free(nullV);
    }

    if (!*err) {
    B = gretl_matrix_alloc(0, Y.cols);
    }

    return B;
}

GretlMatrixRaw * user_matrix_ols(const(GretlMatrixRaw) * Y, const(GretlMatrixRaw) * X, const(GretlMatrixRaw) * F, GretlMatrixRaw * U, GretlMatrixRaw * V, GretlMatrixRaw * A, Gretlopt opt, int * err)
{
    GretlMatrixRaw * B = null;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || X == null) {
    *err = GretlError.data;
    } else if (opt & Gretlopt.OPT_F) {
    /* factorized */
    if (gretl_vector_get_length(F) != Y.rows) {
        *err = GretlError.invarg;
    }
    }

    if (*err) {
    return null;
    }

    if (gretl_is_complex(Y) || gretl_is_complex(X) ||
    gretl_is_complex(U) || gretl_is_complex(V) ||
    gretl_is_complex(F)) {
    fprintf(stderr, "E_CMPLX in user_matrix_ols\n".ptr);
    *err = GretlError.cmplx;
    } else {
    T = Y.rows;
    k = X.cols;
    g = Y.cols;
    if (X.rows != T && X.cols > 0) {
        *err = GretlError.nonconf;
    }
    }

    if (*err) {
    return null;
    }

    if (X.cols == 0 && !(opt & Gretlopt.OPT_F)) {
    /* handle the null model case */
    return null_OLS(Y, U, V, err);
    }

    if (g > 1 && (opt & Gretlopt.OPT_M)) {
    /* multiple precision: we accept only one y var */
    *err = GretlError.invarg;
    return null;
    }

    if (U != null) {
    *err = gretl_matrix_realloc(U, T, g);
    if (*err) {
        return null;
    }
    }

    if (V != null && g == 1) {
    /* single regressand case */
    *err = gretl_matrix_realloc(V, k, k);
    }

    if (!*err) {
    B = gretl_matrix_alloc(k, g);
    if (B == null) {
        *err = GretlError.alloc;
    }
    }

    if (!*err) {
    if (opt & Gretlopt.OPT_F) {
        /* factorized: Y can be single or multiple */
        *err = gretl_matrix_factorized_ols(Y, X, F, B, V, U, A);
    } else if (g == 1) {
        /* single regressand */
        double s2 = 0;
        double * ps2 = (V != null)? &s2 : null;

        if (opt & Gretlopt.OPT_M) {
        /* use multiple precision */
        *err = gretl_matrix_mp_ols(Y, X, B, V, U, ps2);
        } else {
        /* regular C arithmetic */
        *err = gretl_matrix_ols(Y, X, B, V, U, ps2);
        }
    } else {
        /* multiple regressands: @V will actually be (X'X)^{-1} */
        GretlMatrixRaw * Vtmp = null;
        GretlMatrixRaw ** Vp = (V != null)? &Vtmp : null;

        *err = gretl_matrix_multi_ols(Y, X, B, U, Vp);
        if (Vtmp != null) {
        maybe_replace_content(V, Vtmp, *err);
        }
    }
    }

    if (*err) {
    gretl_matrix_free(B);
    B = null;
    }

    return B;
}

GretlMatrixRaw * user_matrix_rls(const(GretlMatrixRaw) * Y, const(GretlMatrixRaw) * X, const(GretlMatrixRaw) * R, const(GretlMatrixRaw) * Q, GretlMatrixRaw * U, GretlMatrixRaw * V, int * err)
{
    GretlMatrixRaw * B = null;
    int g, k, T;

    if (gretl_is_null_matrix(Y) || gretl_is_null_matrix(X)) {
    *err = GretlError.data;
    return null;
    }

    if (gretl_is_complex(Y) || gretl_is_complex(X) ||
    gretl_is_complex(R) || gretl_is_complex(Q)) {
    fprintf(stderr, "E_CMPLX in user_matrix_rls\n".ptr);
    *err = GretlError.cmplx;
    return null;
    }

    T = Y.rows;
    k = X.cols;
    g = Y.cols;

    if (X.rows != T) {
    *err = GretlError.nonconf;
    return null;
    }

    if (U != null) {
    *err = gretl_matrix_realloc(U, T, g);
    if (*err) {
        return null;
    }
    }

    if (!*err) {
    B = gretl_matrix_alloc(k, g);
    if (B == null) {
        *err = GretlError.alloc;
    }
    }

    if (!*err) {
    /* note: &V will actually be M (X'X)^{-1} M' */
    GretlMatrixRaw * Vtmp = null;
    GretlMatrixRaw ** Vp = (V != null)? &Vtmp : null;

    *err = gretl_matrix_restricted_multi_ols(Y, X, R, Q, B,
                         U, Vp);
    if (Vtmp != null) {
        maybe_replace_content(V, Vtmp, *err);
    }
    }

    if (*err) {
    gretl_matrix_free(B);
    B = null;
    }

    return B;
}

GretlMatrixRaw * user_matrix_GHK(const(GretlMatrixRaw) * C, const(GretlMatrixRaw) * A, const(GretlMatrixRaw) * B, const(GretlMatrixRaw) * U, GretlMatrixRaw * dP, int * err)
{
    GretlMatrixRaw * P = null;
    int n, m, npar;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(C)) {
    *err = GretlError.data;
    return null;
    }

    n = A.rows;
    m = C.rows;
    npar = m + m + m*(m+1)/2;

    if (dP != null) {
    *err = gretl_matrix_realloc(dP, n, npar);
    }

    if (!*err) {
    P = gretl_GHK2(C, A, B, U, dP, err);
    }

    return P;
}

GretlMatrixRaw * user_matrix_eigensym(const(GretlMatrixRaw) * m, GretlMatrixRaw * R, int * err)
{
    GretlMatrixRaw * C = null;
    GretlMatrixRaw * E = null;
    int vecs = 0;

    if (gretl_is_null_matrix(m)) {
    *err = GretlError.data;
    return null;
    } else if (gretl_matrix_is_nonfinite(m)) {
    *err = GretlError.nonfin;
    return null;
    }

    if (R != null) {
    /* computing (right) eigenvectors */
    vecs = 1;
    }

    /* @m would be destroyed by this operation */
    C = gretl_matrix_copy(m);
    if (C == null) {
    *err = GretlError.alloc;
    }

    if (!*err) {
    if (m.is_complex) {
        E = gretl_zheev(C, vecs, err);
    } else {
        E = gretl_symmetric_matrix_eigenvals(C, vecs, err);
    }
    }

    if (!*err && vecs) {
    maybe_replace_content(R, C, 0);
    }

    if (!vecs) {
    gretl_matrix_free(C);
    }

    return E;
}

GretlMatrixRaw * user_gensymm_eigenvals(const(GretlMatrixRaw) * A, const(GretlMatrixRaw) * B, GretlMatrixRaw * V, int * err)
{
    GretlMatrixRaw * E = null;

    if (gretl_is_null_matrix(A) || gretl_is_null_matrix(B)) {
    *err = GretlError.data;
    return null;
    }

    if (gretl_matrix_is_nonfinite(A) || gretl_matrix_is_nonfinite(B)) {
    *err = GretlError.nonfin;
    return null;
    }

    if (V != null) {
    *err = gretl_matrix_realloc(V, B.cols, A.rows);
    }

    if (!*err) {
    E = gretl_gensymm_eigenvals(A, B, V, err);
    }

    return E;
}

int gretl_matrix_set_part(GretlMatrixRaw * targ, const(GretlMatrixRaw) * src, double x, SelType sel)
{
    int err = 0;

    if (sel == SelType.diag) {
    if (targ.is_complex) {
        err = gretl_cmatrix_set_diagonal(targ, src, x);
    } else {
        err = gretl_matrix_set_diagonal(targ, src, x);
    }
    } else if (sel == SelType.lower || sel == SelType.upper) {
    int upper = (sel == SelType.upper);

    if (targ.is_complex) {
        err = gretl_cmatrix_set_triangle(targ, src, x, upper);
    } else {
        err = gretl_matrix_set_triangle(targ, src, x, upper);
    }
    } else if (sel == SelType.real_ || sel == SelType.imag) {
    if (targ.is_complex) {
        int im = (sel == SelType.imag);

        err = gretl_cmatrix_set_part(targ, src, x, im);
    } else {
        err = GretlError.types;
    }
    }

    return err;
}

