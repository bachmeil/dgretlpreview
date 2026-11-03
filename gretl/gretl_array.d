/* Port of lib/src/gretl_array.c: gretl arrays (of strings, matrices, bundles, lists and arrays): creation,
   element access, copying, joining, sorting and (de)serialisation. The type GretlArray (struct gretl_array_ in
   the C, opaque in gretl_array.h) is declared in libgretlh.d. The header gretl_array.h has only prototypes. */

module gretl.gretl_array;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : round, nearbyint, atan2, cos, sin, ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fread, fwrite, FILE, fclose, fopen, fgets, fputc, fputs, rewind, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atof, atoi, calloc, free, getenv, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, strtok, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;




import core.stdc.config : c_long;

/* the C library's strtol takes a plain char ** end pointer */
import gretl.dataset : is_string_valued, series_get_string_for_obs;
import gretl.geneval0 : Fncall, Ufunc;
import gretl.geneval1 : set_user_qsorting;
import gretl.genmain : generate_int;
import gretl.gretl_cmatrix : real_to_complex_fill;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_list_cmp, gretl_list_copy, gretl_list_print, gretl_null_list;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_utils : gretl_compare_strings;
import gretl.gretl_xml : gretl_matrix_serialize, gretl_xml_get_list, gretl_xml_get_matrix, gretl_xml_get_prop_as_bool, gretl_xml_get_prop_as_int, gretl_xml_get_string, gretl_xml_get_type_property, gretl_xml_put_tagged_list, gretl_xml_put_tagged_string;
import gretl.gstrfuncs : g_strcmp0;
import gretl.strutils : gretl_strdup, strings_array_dup, strings_array_free, strings_array_new;
import gretl.uservarh : UserVar;
import gretl.gretl_bundle : gretl_bundle_copy, gretl_bundle_deserialize, gretl_bundle_destroy, gretl_bundle_new, gretl_bundle_serialize, gretl_bundles_are_equal;
import gretl.uservar : get_user_var_of_type_by_name, user_var_add, user_var_get_value, user_var_unstack_value;
import gretl.gretl_typemap : gretl_type_from_string, gretl_type_get_name, gretl_type_get_singular;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void fncall_destroy(void * ptr);
    GretlType fncall_get_return_type(Fncall * fc);
    Fncall * fncall_new(Ufunc * fun, int preserve);
    int gretl_function_exec(Fncall * call, int rtype, DatasetRaw * dset, void * ret, Prn * prn);
    int push_anon_function_arg(Fncall * fc, GretlType type, void * value);
    int set_anon_function_arg(Fncall * fc, int i, GretlType type, void * value);
    Ufunc * get_user_function_by_name(const(char) * name);
}

/**
 * gretl_array:
 *
 * An opaque type; use the relevant accessor functions.
 */

void gretl_array_destroy_data(GretlArray * A)
{
    int i;

    if (A.data != null) {
    if (A.type == GretlType.strings ||
        A.type == GretlType.lists) {
        /* a simple "free" will do */
        for (i=0; i<A.n; i++) {
        free(A.data[i]);
        }
    } else if (A.type == GretlType.matrices) {
        if (A.mdata != null) {
        free(A.mdata);
        for (i=0; i<A.n; i++) {
            free(A.data[i]);
        }
        } else {
        for (i=0; i<A.n; i++) {
            gretl_matrix_free(cast(GretlMatrixRaw *) A.data[i]);
        }
        }
    } else if (A.type == GretlType.bundles) {
        for (i=0; i<A.n; i++) {
        gretl_bundle_destroy(cast(GretlBundle *) A.data[i]);
        }
    } else if (A.type == GretlType.arrays) {
        for (i=0; i<A.n; i++) {
        gretl_array_destroy(cast(GretlArray *) A.data[i]);
        }
    }
    free(A.data);
    A.data = null;
    }
}

/* Destroy the whole array, including freeing its
   contents.
*/

void gretl_array_destroy(GretlArray * A)
{
    if (A != null) {
    gretl_array_destroy_data(A);
    free(A);
    }
}

/* Reduce the array to empty status, freeing the
   contents but not the array structure itself.
*/

void gretl_array_void_content(GretlArray * A)
{
    if (A != null) {
    gretl_array_destroy_data(A);
    A.n = 0;
    }
}

/* Reduce the array to empty status, setting the
   entire data array to NULL without freeing anything.
   Makes sense only if the data array is actually
   "owned" elsewhere.
*/

void gretl_array_nullify_content(GretlArray * A)
{
    if (A != null) {
    A.data = null;
    A.n = 0;
    }
}

/* Set the array's element pointers to NULL, without
   freeing them. Makes sense only if the element
   pointers are actually "owned" by something else
   (have not been copied into the array).
*/

void gretl_array_nullify_elements(GretlArray * A)
{
    if (A != null && A.data != null) {
    int i;

    for (i=0; i<A.n; i++) {
        A.data[i] = null;
    }
    }
}

int array_allocate_storage(GretlArray * A)
{
    int i, err = 0;

    A.data = cast(typeof(A.data)) malloc(A.n * (*A.data).sizeof);

    if (A.data == null) {
    err = GretlError.alloc;
    } else {
    for (i=0; i<A.n; i++) {
        A.data[i] = null;
    }
    }

    return err;
}

int array_extend_content(GretlArray * A, int plus)
{
    if (plus == 0) {
    return 0; /* no-op */
    } else if (plus < 0) {
    return GretlError.data;
    } else {
    void ** data;
    int n = A.n + plus;
    int i, err = 0;

    data = cast(typeof(data)) realloc(A.data, n * (*A.data).sizeof);

    if (data == null) {
        err = GretlError.alloc;
    } else {
        for (i=A.n; i<n; i++) {
        data[i] = null;
        }
        A.data = data;
        A.n = n;
    }

    return err;
    }
}

/* Create a new array of type @type. The idea is that it's OK to have
   n = 0: in that case the array starts empty but can be extended by
   the += operator. Or it can be sized in advance.
*/

GretlArray * gretl_array_new(GretlType type, int n, int * err)
{
    GretlArray * A;

    if (type != GretlType.strings &&
    type != GretlType.matrices &&
    type != GretlType.bundles &&
    type != GretlType.lists &&
    type != GretlType.arrays &&
    type != GretlType.any) {
    *err = GretlError.types;
    return null;
    } else if (n < 0) {
    *err = GretlError.data;
    return null;
    }

    A = cast(typeof(A)) malloc((*A).sizeof);

    if (A == null) {
    *err = GretlError.alloc;
    } else {
    A.type = type;
    A.n = n;
    A.data = null;
    A.mdata = null;
    if (n > 0) {
        *err = array_allocate_storage(A);
        if (*err) {
        gretl_array_destroy(A);
        A = null;
        }
    }
    }

    return A;
}

GretlArray * gretl_singleton_array(void * ptr, GretlType atype, int copy, int * err)
{
    GretlArray * A = gretl_array_new(atype, 1, err);

    if (A != null) {
    GretlType t = gretl_type_get_singular(atype);

    *err = gretl_array_set_element(A, 0, ptr, t, copy);
    if (*err) {
        free(A);
        A = null;
    }
    }

    return A;
}

GretlArray * gretl_array_from_strings(char ** S, int n, int copy, int * err)
{
    GretlArray * A;

    A = gretl_array_new(GretlType.strings, 0, err);

    if (A != null && n > 0) {
    if (copy) {
        A.data = cast(void **) strings_array_dup(S, n);
        if (A.data == null) {
        *err = GretlError.alloc;
        }
    } else {
        A.data = cast(void **) S;
    }
    if (!*err) {
        A.n = n;
    }
    }

    return A;
}

GretlArray * gretl_matrix_array_sized(int n, int r, int c, int * err)
{
    GretlArray * A;

    A = gretl_array_new(GretlType.matrices, n, err);

    if (A != null && n > 0) {
    size_t msize = n * r * c * double.sizeof;
    double * ai_val;
    GretlMatrixRaw * ai;
    int i, rc = r * c;

    /* common block for matrix data */
    ai_val = A.mdata = cast(typeof(A.mdata)) malloc(msize);
    if (A.mdata == null) {
        *err = GretlError.alloc;
    }

    for (i=0; i<n && !*err; i++) {
        ai = gretl_null_matrix_new();
        if (ai == null) {
        *err = GretlError.alloc;
        break;
        }
        ai.val = ai_val;
        ai.rows = r;
        ai.cols = c;
        A.data[i] = ai;
        ai_val += rc;
    }
    }

    if (*err && A != null) {
    gretl_array_destroy(A);
    A = null;
    }

    return A;
}

/* When we're returning an array of strings, ensure
   that any NULL elements are converted to empty
   strings.
*/

int strings_array_null_check(GretlArray * A)
{
    int i;

    for (i=0; i<A.n; i++) {
    if (A.data[i] == null) {
        A.data[i] = gretl_strdup("");
        if (A.data[i] == null) {
        return GretlError.alloc;
        }
    }
    }

    return 0;
}

/* note: don't modify the returned value */

char ** gretl_array_get_strings(GretlArray * A, int * ns)
{
    char ** AS = null;

    *ns = 0;

    if (A != null && A.type == GretlType.strings) {
    int err = strings_array_null_check(A);

    if (!err) {
        *ns = A.n;
        AS = cast(char **) A.data;
    }
    }

    return AS;
}

/* note: take ownership of the returned value */

char ** gretl_array_steal_strings(GretlArray * A, int * ns)
{
    char ** AS = null;

    *ns = 0;

    if (A != null && A.type == GretlType.strings) {
    int err = strings_array_null_check(A);

    if (!err) {
        *ns = A.n;
        AS = cast(char **) A.data;
        A.n = 0;
        A.data = null;
    }
    }

    return AS;
}

char ** gretl_array_get_stringify_strings(GretlArray * A, int nreq, int * pns, int * err)
{
    char ** S = null;

    *pns = 0;

    if (A == null) {
    *err = GretlError.data;
    } else if (A.type != GretlType.strings) {
    *err = GretlError.types;
    } else if (A.n < nreq) {
    gretl_errmsg_sprintf("Too few strings: %d given but %d needed".ptr,
                 A.n, nreq);
    *err = GretlError.data;
    }

    if (!*err) {
    char ** AS = cast(char **) A.data;

    S = strings_array_new(A.n);
    if (S == null) {
        *err = GretlError.alloc;
    } else {
        int myerr = 0;
        int ndone = 0;
        int i, j;

        for (i=0; i<A.n && !myerr; i++) {
        if (AS[i] == null || AS[i][0] == '\0') {
            myerr = GretlError.data;
        } else {
            S[i] = gretl_strdup(AS[i]);
            if (i > 0) {
            for (j=0; j<i; j++) {
                if (!strcmp(AS[j], S[i])) {
                gretl_errmsg_sprintf("Duplicated string '%s'".ptr, S[i]);
                myerr = GretlError.data;
                }
            }
            }
        }
        if (!myerr) {
            ndone++;
        }
        }
        if (myerr && ndone < nreq) {
        *err = myerr;
        strings_array_free(S, A.n);
        S = null;
        } else {
        *pns = ndone;
        }
    }
    }

    return S;
}

/* note: the return value is newly allocated, and owned by the caller */

char * gretl_strings_array_flatten(GretlArray * A, const(char) * sep, int * err)
{
    char * s = null;

    if (A == null) {
    *err = GretlError.data;
    } else if (A.type != GretlType.strings) {
    *err = GretlError.types;
    } else {
        int ns = cast(int) (cast(int) strlen(sep));
    int i, len = 1; /* for terminating NUL */

    for (i=0; i<A.n; i++) {
        if (A.data[i] == null) {
        *err = GretlError.missdata;
        break;
        } else {
        len += strlen(cast(const(char) *) A.data[i]);
                if (i < A.n - 1) {
                    len += ns;
                }
        }
    }

    if (!*err) {
        s = cast(typeof(s)) calloc(len, 1);
        if (s == null) {
        *err = GretlError.alloc;
        } else {
        for (i=0; i<A.n; i++) {
            strcat(s, cast(const(char) *) A.data[i]);
            if (i < A.n - 1) {
            strcat(s, sep);
            }
        }
        }
    }
    }

    return s;
}

/* Return a column vector holding the position(s) in
   the strings array @A at which the string @s is
   matched -- or an empty matrix in case of no matches.
*/

GretlMatrixRaw * gretl_strings_array_pos(GretlArray * A, const(char) * s, int * err)
{
    GretlMatrixRaw * ret = null;
    const(char) * si;
    int i, np = 0;

    for (i=0; i<A.n; i++) {
    si = (A.data[i] == null)? "".ptr : cast(const(char) *) A.data[i];
    if (strcmp(si, s) == 0) {
        np++;
    }
    }

    if (np == 0) {
    ret = gretl_null_matrix_new();
    } else {
    ret = gretl_matrix_alloc(np, 1);
    }

    if (ret == null) {
    *err = GretlError.alloc;
    } else if (np > 0) {
    int j = 0;

    for (i=0; i<A.n; i++) {
        si = (A.data[i] == null)? "".ptr : cast(const(char) *) A.data[i];
        if (strcmp(si, s) == 0) {
        ret.val[j++] = i+1;
        }
    }
    }

    return ret;
}

/* Return 1 if @A includes @s, otherwise 0 */

int gretl_strings_array_includes(GretlArray * A, const(char) * s)
{
    int ret = 0;
    const(char) * si;
    int i;

    for (i=0; i<A.n; i++) {
    si = (A.data[i] == null)? "".ptr : cast(const(char) *) A.data[i];
    if (strcmp(si, s) == 0) {
        ret = 1;
        break;
    }
    }

    return ret;
}

void * gretl_array_get_element(GretlArray * A, int i, GretlType * type, int * err)
{
    void * ret = null;

    /* Note that the data returned here are not "deep copied", we just
       pass the pointer. It's up to the caller to decide if a copy has
       to be made, given that the pointer from here should not be
       modified.
    */

    if (A == null) {
    *err = GretlError.unkvar;
    } else if (i < 0 || i >= A.n) {
    *err = GretlError.invarg;
    } else if (A.type == GretlType.arrays && A.data[i] == null) {
    *err = GretlError.invarg;
    }
    if (*err) {
    return null;
    }

    if (type != null) {
    *type = gretl_type_get_singular(A.type);
    }
    if (A.type == GretlType.strings) {
    if (A.data[i] == null) {
        A.data[i] = gretl_strdup("");
    }
    } else if (A.type == GretlType.matrices) {
    if (A.data[i] == null) {
        A.data[i] = gretl_null_matrix_new();
    }
    } else if (A.type == GretlType.bundles) {
    if (A.data[i] == null) {
        A.data[i] = gretl_bundle_new();
    }
    } else if (A.type == GretlType.lists) {
    if (A.data[i] == null) {
        A.data[i] = gretl_null_list();
    }
    }

    ret = A.data[i];
    if (ret == null) {
    *err = GretlError.alloc;
    }

    return ret;
}

/* Note: no type-checking, the caller is supposed to
   be on the ball. Use with caution.
*/

void * gretl_array_get_data(GretlArray * A, int i)
{
    if (A == null || i < 0 || i >= A.n) {
    return null;
    } else {
    return A.data[i];
    }
}

int gretl_array_set_data(GretlArray * A, int i, void * ptr)
{
    if (A == null || i < 0 || i >= A.n) {
    return GretlError.data;
    } else {
    A.data[i] = ptr;
    return 0;
    }
}

void * gretl_array_get_all_data(GretlArray * A)
{
    if (A == null) {
    return null;
    } else {
    return A.data;
    }
}

int check_list_bounds(int * list, int arrdim)
{
    int i;

    for (i=1; i<=list[0]; i++) {
    if (list[i] <= 0 || list[i] > arrdim) {
        return 0;
    }
    }

    return 1;
}

GretlArray * gretl_array_copy_subspec(GretlArray * A, int * list, int * err)
{
    GretlArray * C = null;

    if (A == null) {
    *err = GretlError.data;
    } else if (!check_list_bounds(list, A.n)) {
    *err = GretlError.invarg;
    } else {
    int i, k, m = list[0];

    C = gretl_array_new(A.type, m, err);

    for (i=0; i<m && !*err; i++) {
        k = list[i+1] - 1;
        if (A.data[k] != null) {
        if (A.type == GretlType.strings) {
            C.data[i] = gretl_strdup(cast(const(char) *) A.data[k]);
        } else if (A.type == GretlType.matrices) {
            C.data[i] = gretl_matrix_copy(cast(const(GretlMatrixRaw) *) A.data[k]);
        } else if (A.type == GretlType.bundles) {
            C.data[i] = gretl_bundle_copy(cast(GretlBundle *) A.data[k], err);
        } else if (A.type == GretlType.arrays) {
            C.data[i] = gretl_array_copy(cast(const(GretlArray) *) A.data[k], err);
        } else if (A.type == GretlType.lists) {
            C.data[i] = gretl_list_copy(cast(const(int) *) A.data[k]);
        }
        if (!*err && C.data[i] == null) {
            *err = GretlError.alloc;
        }
        }
    }
    }

    if (*err && C != null) {
    gretl_array_destroy(C);
    C = null;
    }

    return C;
}

GretlBundle * gretl_array_get_bundle(GretlArray * A, int i)
{
    GretlBundle * b = null;

    /* The bundle returned here is not "deep copied",
       we just pass the pointer. It's up to the caller to
       decide if a copy has to be made, given that the
       pointer from here should not be modified.
    */

    if (A != null && i >= 0 && i < A.n &&
    A.type == GretlType.bundles) {
    b = cast(GretlBundle *) A.data[i];
    }

    return b;
}

/* "flatten" an array of matrices, yielding a single matrix.
   Concatenation may be horizontal (@mode = 0, default), vertical
   (@mode = 1) or (for lack of a better term) "vec-wise" (@mode =
   2). An error is flagged if the matrices are not conformable for the
   operation.
*/

GretlMatrixRaw * gretl_matrix_array_flatten(GretlArray * A, int mode, int * err)
{
    enum {HORIZ, VERT, VECWISE};
    GretlMatrixRaw * ret = null;
    GretlMatrixRaw * m;
    int sum_r = 0;
    int sum_c = 0;
    int r = 0, c = 0;
    int nvals = 0;
    int cmplx = 0;
    int real_ = 0;
    int nmats = 0;
    int i;

    if (A.type != GretlType.matrices) {
    *err = GretlError.types;
    return null;
    }

    for (i=0; i<A.n && !*err; i++) {
    m = cast(GretlMatrixRaw *) A.data[i];
    if (gretl_is_null_matrix(m)) {
            continue;
        }
        nmats++;
        if (m.is_complex) {
            cmplx = 1;
        } else {
            real_ = 1;
        }
        if (c == 0) {
            /* initialize */
            r = m.rows;
            c = m.cols;
            nvals = r * c;
        } else {
            /* check conformability */
            if (mode == HORIZ && m.rows != r) {
                *err = GretlError.nonconf;
            } else if (mode == VERT && m.cols != c) {
                *err = GretlError.nonconf;
            } else if (mode == VECWISE && m.rows * m.cols != nvals) {
                *err = GretlError.nonconf;
            }
        }
        sum_r += m.rows;
        sum_c += m.cols;
    }

    if (*err) {
        return null;
    }

    if (mode == HORIZ) {
        c = sum_c;
    } else if (mode == VERT) {
        r = sum_r;
    }

    if (mode == VECWISE) {
    if (cmplx) {
        ret = gretl_cmatrix_new(nvals, nmats);
    } else {
        ret = gretl_matrix_alloc(nvals, nmats);
    }
    } else if (cmplx) {
        ret = gretl_cmatrix_new(r, c);
    } else {
        ret = gretl_matrix_alloc(r, c);
    }

    if (ret == null) {
        *err = GretlError.alloc;
        return ret;
    } else if (gretl_is_null_matrix(ret)) {
        return ret;
    }

    if (mode == HORIZ) {
        if (cmplx && real_) {
            /* horizontal: mixed matrices */
            Cmplx * dest = ret.z;
            int n, k = 0;

            for (i=0; i<A.n; i++) {
                m = cast(GretlMatrixRaw *) A.data[i];
                if (!gretl_is_null_matrix(m)) {
                    n = m.rows * m.cols;
                    if (m.is_complex) {
                        memcpy(dest, m.z, n * (*dest).sizeof);
                    } else {
                        real_to_complex_fill(ret, m, 0, k);
                    }
                    dest += n;
                    k += m.cols;
                }
            }
        } else {
            /* horizontal: the easy case */
            double * dest = ret.val;
            int n, p = cmplx ? 2 : 1;

            for (i=0; i<A.n; i++) {
                m = cast(GretlMatrixRaw *) A.data[i];
                if (!gretl_is_null_matrix(m)) {
                    n = p * m.rows * m.cols;
                    memcpy(dest, m.val, n * (*dest).sizeof);
                    dest += n;
                }
            }
        }
    } else if (mode == VERT) {
    /* vertical concatenation */
    int ii, j, k, rpos = 0;
    Cmplx z;
    double x;

    for (k=0; k<A.n; k++) {
        m = cast(GretlMatrixRaw *) A.data[k];
        if (!gretl_is_null_matrix(m)) {
        for (j=0; j<c; j++) {
            for (i=0, ii=rpos; i<m.rows; i++, ii++) {
            if (cmplx && !m.is_complex) {
                z = Cmplx(gretl_matrix_get(m, i, j), 0.0);
                gretl_cmatrix_set(ret, ii, j, z);
            } else if (cmplx) {
                z = gretl_cmatrix_get(m, i, j);
                gretl_cmatrix_set(ret, ii, j, z);
            } else {
                x = gretl_matrix_get(m, i, j);
                gretl_matrix_set(ret, ii, j, x);
            }
            }
        }
        rpos += m.rows;
        }
    }
    } else if (mode == VECWISE) {
    if (cmplx && real_) {
        /* vecwise: mixed matrices */
        Cmplx * dest = ret.z;
        Cmplx z;
            int k = 0;

        for (i=0; i<A.n; i++) {
        m = cast(GretlMatrixRaw *) A.data[i];
        if (!gretl_is_null_matrix(m)) {
            if (m.is_complex) {
            memcpy(dest, m.z, nvals * (*dest).sizeof);
            } else {
            for (k=0; k<nvals; k++) {
                z = Cmplx(m.val[k], 0.0);
                gretl_cmatrix_set(ret, k, i, z);
            }
            }
            dest += nvals;
        }
        }
    } else {
        /* vecwise: the easy case */
        double * dest = ret.val;
        int n, p = cmplx ? 2 : 1;

        for (i=0; i<A.n; i++) {
        m = cast(GretlMatrixRaw *) A.data[i];
        if (!gretl_is_null_matrix(m)) {
            n = p * nvals;
            memcpy(dest, m.val, n * (*dest).sizeof);
            dest += n;
        }
        }
    }
    }

    return ret;
}

int split_matrix_by_chunks(GretlArray * A, int n, const(GretlMatrixRaw) * X, int chunk, int colwise)
{
    int rows = colwise ? X.rows : chunk;
    int cols = colwise ? chunk : X.cols;
    GretlMatrixRaw * tmp;
    int i, err = 0;

    /* allocate the matrices */
    for (i=0; i<n && !err; i++) {
    tmp = gretl_matrix_alloc(rows, cols);
    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        gretl_array_set_element(A, i, tmp, GretlType.matrix, 0);
    }
    }

    /* fill the matrices */
    if (!err && colwise) {
    int nelem = rows * cols;
    size_t csize = nelem * (*X.val).sizeof;
    const(double) * src = X.val;

    for (i=0; i<n; i++) {
        tmp = cast(GretlMatrixRaw *) A.data[i];
        memcpy(tmp.val, src, csize);
        src += nelem;
    }
    } else if (!err) {
    int k, j, ii = 0;
    double x;

    for (k=0; k<n; k++) {
        tmp = cast(GretlMatrixRaw *) A.data[k];
        for (j=0; j<cols; j++) {
        for (i=0; i<rows; i++) {
            x = gretl_matrix_get(X, ii+i, j);
            gretl_matrix_set(tmp, i, j, x);
        }
        }
        ii += rows;
    }
    }

    return err;
}

int split_matrix_by_vector(GretlArray * A, int n, const(GretlMatrixRaw) * X, int dim, int colwise, GretlMatrixRaw * vals, const(double) * sel)
{
    int * nj = null;
    GretlMatrixRaw * tmp;
    int rows, cols;
    int i, j, k;
    int err = 0;

    /* How many rows or columns will each matrix need? */
    nj = cast(typeof(nj)) calloc(n, (*nj).sizeof);
    if (nj == null) {
    err = GretlError.alloc;
    } else {
    for (i=0; i<dim; i++) {
        j = cast(int) (sel[i] - 1);
        nj[j] += 1;
    }
    }
    for (i=0; i<n && !err; i++) {
    if (nj[i] > 0) {
        rows = colwise ? X.rows : nj[i];
        cols = colwise ? nj[i] : X.cols;
        tmp = gretl_matrix_alloc(rows, cols);
        if (tmp == null) {
        err = GretlError.alloc;
        } else {
        gretl_array_set_element(A, i, tmp, GretlType.matrix, 0);
        /* @nj will serve as an array of counters below */
        nj[i] = 0;
        }
    }
    }

    /* fill the individual matrices */
    if (!err && colwise) {
    size_t csize = X.rows * (*X.val).sizeof;
    const(double) * src = X.val;
    double * targ;

    for (i=0; i<X.cols; i++) {
        k = cast(int) (sel[i] - 1);
        tmp = cast(GretlMatrixRaw *) A.data[k];
        targ = tmp.val + nj[k];
        memcpy(targ, src, csize);
        /* advance the read and write positions */
        src += X.rows;
        nj[k] += X.rows;
    }
    } else if (!err) {
    double x;
    int jj;

    for (i=0; i<X.rows; i++) {
        k = cast(int) (sel[i] - 1);
        tmp = cast(GretlMatrixRaw *) A.data[k];
        jj = 0;
        for (j=0; j<X.cols; j++) {
        x = gretl_matrix_get(X, i, j);
        gretl_matrix_set(tmp, nj[k], jj++, x);
        }
        /* advance the write position */
        nj[k] += 1;
    }
    }

    free(nj);

    return err;
}

GretlArray * gretl_matrix_split_by(const(GretlMatrixRaw) * X, const(GretlMatrixRaw) * v, int colwise, int chunks, int * err)
{
    GretlArray * ret = null;
    int i, dim, nm = 0;
    int chunksize = 0;
    GretlMatrixRaw * vals = null;
    const(double) * sel;
    double x;

    dim = colwise ? X.cols : X.rows;

    if (chunks) {
    /* interpret single element of @v as chunk size */
    chunksize = cast(int) v.val[0];
    if (chunksize <= 0 || chunksize > dim) {
        *err = GretlError.invarg;
    } else if (dim % chunksize != 0) {
        *err = GretlError.nonconf;
    }
    } else if (gretl_vector_get_length(v) == dim) {
    /* vector of indices: only positive integers allowed */
    for (i=0; i<dim; i++) {
        x = v.val[i];
        if (x != floor(x) || x <= 0 || x >= int.max) {
        *err = GretlError.invarg;
        break;
        }
    }
    } else {
    *err = GretlError.invarg;
    }

    if (*err) {
    return null;
    }

    /* How many matrices do we need ? */
    if (chunksize > 0) {
    nm = dim / chunksize;
    } else {
    sel = v.val;
    vals = gretl_matrix_values(sel, dim, Gretlopt.OPT_NONE, err);
    if (!*err) {
        /* get the maximum index value */
        for (i=0; i<vals.rows; i++) {
        nm = cast(int) (vals.val[i] > nm ? vals.val[i] : nm);
        }
    }
    }

    if (!*err) {
    ret = gretl_array_new(GretlType.matrices, nm, err);
    }

    if (!*err && chunksize > 0) {
    /* the easier case */
    *err = split_matrix_by_chunks(ret, nm, X, chunksize, colwise);
    } else if (!*err) {
    /* the fiddly but general one */
    *err = split_matrix_by_vector(ret, nm, X, dim, colwise, vals, sel);
    }

    gretl_matrix_free(vals);

    if (*err && ret != null) {
    gretl_array_destroy(ret);
    ret = null;
    }

    return ret;
}

int gretl_array_set_type(GretlArray * A, GretlType type)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (type != GretlType.strings &&
           type != GretlType.matrices &&
           type != GretlType.bundles &&
           type != GretlType.lists &&
           type != GretlType.arrays) {
    err = GretlError.types;
    } else if (type == A.type) {
    /* no-op */
    return 0;
    } else if (A.n > 0) {
    /* we can (re-)set the type only if no data have
       been entered already
    */
    int i;

    for (i=0; i<A.n; i++) {
        if (A.data[i] != null) {
        err = GretlError.data;
        break;
        }
    }
    }

    if (!err) {
    A.type = type;
    }

    return err;
}

GretlType gretl_array_get_type(GretlArray * A)
{
    return (A != null)? A.type : GretlType.none;
}

GretlType gretl_array_get_content_type(GretlArray * A)
{
    return (A != null)? gretl_type_get_singular(A.type) : GretlType.none;
}

int gretl_array_get_length(GretlArray * A)
{
    return (A != null)? A.n : 0;
}

/* Return the 0-based index of the first empty slot
   in @A, or -1 on failure (if @A is NULL or all its
   slots are already filled).
*/

int gretl_array_get_next_index(GretlArray * A)
{
    int ret = -1;

    if (A != null) {
    int i;

    for (i=0; i<A.n; i++) {
        if (A.data[i] == null) {
        ret = i;
        break;
        }
    }
    }

    return ret;
}

int set_string(GretlArray * A, int i, char * s, int copy)
{
    int err = 0;

    if (copy) {
    A.data[i] = gretl_strdup(s);
    if (A.data[i] == null) {
        err = GretlError.alloc;
    }
    } else {
    A.data[i] = s;
    }

    return err;
}

/* The static set_*() functions assume that error-checking
   has already been done.
*/

int set_matrix(GretlArray * A, int i, GretlMatrixRaw * m, int copy)
{
    int err = 0;

    if (copy) {
    A.data[i] = gretl_matrix_copy(m);
    if (A.data[i] == null) {
        err = GretlError.alloc;
    }
    } else {
    A.data[i] = m;
    }

    return err;
}

int set_bundle(GretlArray * A, int i, GretlBundle * b, int copy)
{
    int err = 0;

    if (copy) {
    A.data[i] = gretl_bundle_copy(b, &err);
    } else {
    A.data[i] = b;
    }

    return err;
}

int set_array(GretlArray * A, int i, GretlArray * a, int copy)
{
    int err = 0;

    if (copy) {
    A.data[i] = gretl_array_copy(a, &err);
    } else {
    A.data[i] = a;
    }

    return err;
}

int set_list(GretlArray * A, int i, int * L, int copy)
{
    int err = 0;

    if (copy) {
    A.data[i] = gretl_list_copy(L);
    if (A.data[i] == null) {
        err = GretlError.alloc;
    }
    } else {
    A.data[i] = L;
    }

    return err;
}

int set_type_error(GretlArray * A, GretlType required)
{
    if (A.type == GretlType.any) {
    A.type = required;
    return 0;
    } else if (A.type != required) {
    GretlType reqt = gretl_type_get_singular(required);

    gretl_errmsg_sprintf("Cannot add %s to array of %s".ptr,
                 gretl_type_get_name(reqt),
                 gretl_type_get_name(A.type));
    return 1;
    } else {
    return 0;
    }
}

/* In the functions below we assume the @copy parameter
   will be set appropriately by "genr": if the incoming
   object is a named variable in its own right it should
   be copied, but if it's on on-the-fly thing there's
   no need to copy it, just "donate" it to the array.
*/

/* respond to A[i] = s */

int gretl_array_set_string(GretlArray * A, int i, char * s, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.strings)) {
    err = GretlError.types;
    } else if (i < 0 || i >= A.n) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    err = GretlError.data;
    } else if (s != A.data[i]) {
    free(A.data[i]);
    err = set_string(A, i, s, copy);
    }

    return err;
}

/* respond to A += s */

int gretl_array_append_string(GretlArray * A, char * s, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.strings)) {
    err = GretlError.types;
    } else {
    err = array_extend_content(A, 1);
    if (!err) {
        err = set_string(A, A.n - 1, s, copy);
    }
    }

    return err;
}

/* respond to A[i] = m */

int gretl_array_set_matrix(GretlArray * A, int i, GretlMatrixRaw * m, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.matrices)) {
    err = GretlError.types;
    } else if (i < 0 || i >= A.n) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    err = GretlError.data;
    } else if (m != A.data[i]) {
    gretl_matrix_free(cast(GretlMatrixRaw *) A.data[i]);
    err = set_matrix(A, i, m, copy);
    }

    return err;
}

/* respond to A += m */

int gretl_array_append_matrix(GretlArray * A, GretlMatrixRaw * m, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.matrices)) {
    err = GretlError.types;
    } else {
    err = array_extend_content(A, 1);
    if (!err) {
        err = set_matrix(A, A.n - 1, m, copy);
    }
    }

    return err;
}

/* respond to A[i] = b */

int gretl_array_set_bundle(GretlArray * A, int i, GretlBundle * b, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.bundles)) {
    err = GretlError.types;
    } else if (i < 0 || i >= A.n) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    err = GretlError.data;
    } else if (b != A.data[i]) {
    gretl_bundle_destroy(cast(GretlBundle *) A.data[i]);
    err = set_bundle(A, i, b, copy);
    }

    return err;
}

/* respond to A[i] = a */

int gretl_array_set_array(GretlArray * A, int i, GretlArray * a, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.arrays)) {
    err = GretlError.types;
    } else if (i < 0 || i >= A.n) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    err = GretlError.data;
    } else if (a != A.data[i]) {
    gretl_array_destroy(cast(GretlArray *) A.data[i]);
    err = set_array(A, i, a, copy);
    }

    return err;
}

/* respond to A += b */

int gretl_array_append_bundle(GretlArray * A, GretlBundle * b, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.bundles)) {
    err = GretlError.types;
    } else {
    err = array_extend_content(A, 1);
    if (!err) {
        err = set_bundle(A, A.n - 1, b, copy);
    }
    }

    return err;
}

/* respond to A += a */

int gretl_array_append_array(GretlArray * A, GretlArray * a, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.arrays)) {
    err = GretlError.types;
    } else {
    err = array_extend_content(A, 1);
    if (!err) {
        err = set_array(A, A.n - 1, a, copy);
    }
    }

    return err;
}

/* respond to A[i] = L */

int gretl_array_set_list(GretlArray * A, int i, int * L, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.lists)) {
    err = GretlError.types;
    } else if (i < 0 || i >= A.n) {
    gretl_errmsg_sprintf(_("Index value %d is out of bounds").ptr, i+1);
    err = GretlError.data;
    } else if (L != A.data[i]) {
    free(A.data[i]);
    err = set_list(A, i, L, copy);
    }

    return err;
}

/* respond to A += L */

int gretl_array_append_list(GretlArray * A, int * L, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (set_type_error(A, GretlType.lists)) {
    err = GretlError.types;

    } else {
    err = array_extend_content(A, 1);
    if (!err) {
        err = set_list(A, A.n - 1, L, copy);
    }
    }

    return err;
}

int gretl_array_set_element(GretlArray * A, int i, void * ptr, GretlType type, int copy)
{
    int err = 0;

    if (type == GretlType.matrix) {
    err = gretl_array_set_matrix(A, i, cast(GretlMatrixRaw *) ptr, copy);
    } else if (type == GretlType.string_) {
    err = gretl_array_set_string(A, i, cast(char *) ptr, copy);
    } else if (type == GretlType.bundle) {
    err = gretl_array_set_bundle(A, i, cast(GretlBundle *) ptr, copy);
    } else if (type == GretlType.list) {
    err = gretl_array_set_list(A, i, cast(int *) ptr, copy);
    } else if (type == GretlType.array) {
    err = gretl_array_set_array(A, i, cast(GretlArray *) ptr, copy);
    }

    return err;
}

void free_array_element(GretlArray * A, int i)
{
    if (A.type == GretlType.strings ||
    A.type == GretlType.lists) {
    free(A.data[i]);
    } else if (A.type == GretlType.matrices) {
    gretl_matrix_free(cast(GretlMatrixRaw *) A.data[i]);
    } else if (A.type == GretlType.bundles) {
    gretl_bundle_destroy(cast(GretlBundle *) A.data[i]);
    } else if (A.type == GretlType.arrays) {
    gretl_array_destroy(cast(GretlArray *) A.data[i]);
    }
}

int gretl_array_delete_element(GretlArray * A, int i)
{
    if (A == null) {
    return GretlError.data;
    } else if (i < 0 || i >= A.n) {
    return GretlError.invarg;
    } else {
    int j;

    if (A.data[i] != null) {
        free_array_element(A, i);
    }
    /* shift the higher-numbered elements down */
    for (j=i; j<A.n-1; j++) {
        A.data[j] = A.data[j+1];
    }
    /* decrement the element count */
    A.n -= 1;

    return 0;
    }
}

/* Drop any/all occurrences of string @s from array @A */

int gretl_array_drop_string(GretlArray * A, const(char) * s)
{
    if (A.type != GretlType.strings) {
    return GretlError.types;
    } else if (A.data == null || A.n == 0 || s == null) {
    gretl_errmsg_set("Array is empty");
    return GretlError.invarg;
    } else {
    int i, j, rem = A.n;
    int n_orig = A.n;
    size_t sz;

    for (i=0; ; ) {
        if (A.data[i] != null && !strcmp(cast(const(char) *) A.data[i], s)) {
        free(A.data[i]);
        j = i + 1;
        rem = A.n - j;
        sz = rem * (*A.data).sizeof;
        memmove(A.data + i, A.data + j, sz);
        A.n -= 1;
        } else {
        i++;
        }
        if (i == A.n || rem == 0) {
        break;
        }
    }
    if (A.n == 0) {
        free(A.data);
        A.data = null;
    } else if (A.n < n_orig) {
        A.data = cast(typeof(A.data)) realloc(A.data, A.n * (*A.data).sizeof);
    }
    }

    return 0;
}

int gretl_array_drop_null(GretlArray * A)
{
    int i, j, rem = A.n;
    int n_orig = A.n;
    size_t sz;

    if (A.data == null || rem == 0) {
    gretl_errmsg_set("Array is empty");
    return GretlError.invarg;
    }

    for (i=0; ; ) {
    if (A.data[i] == null) {
        j = i + 1;
        rem = A.n - j;
        sz = rem * (*A.data).sizeof;
        memmove(A.data + i, A.data + j, sz);
        A.n -= 1;
    } else {
        i++;
    }
    if (i == A.n || rem == 0) {
        break;
    }
    }
    if (A.n == 0) {
    free(A.data);
    A.data = null;
    } else if (A.n < n_orig) {
    A.data = cast(typeof(A.data)) realloc(A.data, A.n * (*A.data).sizeof);
    }

    return 0;
}

/* @ptr must be pre-checked as matching the array type */

int gretl_array_append_object(GretlArray * A, void * ptr, int copy)
{
    int err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (A.type == GretlType.matrices) {
    gretl_array_append_matrix(A, cast(GretlMatrixRaw *) ptr, copy);
    } else if (A.type == GretlType.strings) {
    gretl_array_append_string(A, cast(char *) ptr, copy);
    } else if (A.type == GretlType.bundles) {
    gretl_array_append_bundle(A, cast(GretlBundle *) ptr, copy);
    } else if (A.type == GretlType.lists) {
    gretl_array_append_list(A, cast(int *) ptr, copy);
    } else if (A.type == GretlType.arrays) {
    gretl_array_append_array(A, cast(GretlArray *) ptr, copy);
    }

    return err;
}

int gretl_array_copy_content(GretlArray * Acpy, const(GretlArray) * A, int write_offset)
{
    int i, j, err = 0;

    for (i=0; i<A.n && !err; i++) {
    if (A.data[i] != null) {
        j = i + write_offset;
        if (A.type == GretlType.strings) {
        Acpy.data[j] = gretl_strdup(cast(const(char) *) A.data[i]);
        } else if (A.type == GretlType.matrices) {
        Acpy.data[j] = gretl_matrix_copy(cast(const(GretlMatrixRaw) *) A.data[i]);
        } else if (A.type == GretlType.bundles) {
        Acpy.data[j] = gretl_bundle_copy(cast(GretlBundle *) A.data[i], &err);
        } else if (A.type == GretlType.arrays) {
        Acpy.data[j] = gretl_array_copy(cast(const(GretlArray) *) A.data[i], &err);
        } else if (A.type == GretlType.lists) {
        Acpy.data[j] = gretl_list_copy(cast(const(int) *) A.data[i]);
        } else {
        err = GretlError.types;
        }
        if (!err && Acpy.data[j] == null) {
        err = GretlError.alloc;
        }
    }
    }

    return err;
}

GretlArray * gretl_array_copy(const(GretlArray) * A, int * err)
{
    GretlArray * Acpy = null;

    if (A != null) {
    Acpy = gretl_array_new(A.type, A.n, err);
    if (!*err) {
        *err = gretl_array_copy_content(Acpy, A, 0);
    }
    }

    return Acpy;
}

extern (C) int inverse_compare_strings(const(void) * a, const(void) * b)
{
    const(char) ** sa = cast(const(char) **) a;
    const(char) ** sb = cast(const(char) **) b;

    return -g_strcmp0(*sa, *sb);
}

GretlArray * gretl_strings_sort(const(GretlArray) * A, int descending, int * err)
{
    GretlArray * Asrt = null;

    if (A != null) {
    if (A.type != GretlType.strings) {
        *err = GretlError.types;
    } else {
        Asrt = gretl_array_new(A.type, A.n, err);
    }
    if (!*err) {
        *err = gretl_array_copy_content(Asrt, A, 0);
    }
    if (!*err) {
        qsort(Asrt.data, Asrt.n, (*Asrt.data).sizeof,
          descending ? &inverse_compare_strings : &gretl_compare_strings);
    }
    }

    return Asrt;
}

/* respond to A1 += A2 */

int gretl_array_copy_into(GretlArray * A1, const(GretlArray) * A2)
{
    GretlArray * Tmp = null;
    int old_n = 0, err = 0;

    if (A1 == null || A2 == null) {
    err = GretlError.data;
    } else if (A1.type != A2.type) {
    err = GretlError.types;
    } else if (A1 == A2) {
        Tmp = gretl_array_copy(A1, &err);
        if (!err) {
            A2 = Tmp;
        }
    }

    if (!err) {
    old_n = A1.n;
    err = array_extend_content(A1, A2.n);
    }
    if (!err) {
    err = gretl_array_copy_content(A1, A2, old_n);
    }
    if (Tmp != null) {
        gretl_array_destroy(Tmp);
    }

    return err;
}

/* respond to C = A + B */

GretlArray * gretl_arrays_join(GretlArray * A, GretlArray * B, int * err)
{
    GretlArray * C = null;

    if (A == null || B == null) {
    *err = GretlError.data;
    } else if (A.type != B.type) {
    *err = GretlError.types;
    } else {
    int n = A.n + B.n;

    C = gretl_array_new(A.type, n, err);
    }

    if (!*err) {
    *err = gretl_array_copy_content(C, A, 0);
    }

    if (!*err) {
    *err = gretl_array_copy_content(C, B, A.n);
    }

    if (*err && C != null) {
    gretl_array_destroy(C);
    C = null;
    }

    return C;
}

/* respond to C = A || B */

GretlArray * gretl_arrays_union(GretlArray * A, GretlArray * B, int * err)
{
    GretlArray * C = null;
    const(char) * sa, sb;
    char * copy = null;
    int i, j, n = 0;

    if (A == null || B == null) {
    *err = GretlError.data;
    } else if (A.type != GretlType.strings ||
           B.type != GretlType.strings) {
    *err = GretlError.types;
    } else {
    if (B.n > 0) {
        copy = cast(typeof(copy)) calloc(1, B.n);
    }
    n = A.n;
    for (j=0; j<B.n; j++) {
        sb = cast(const(char) *) B.data[j];
        if (sb == null || *sb == '\0') {
        continue;
        }
        copy[j] = 1;
        for (i=0; i<A.n; i++) {
        sa = cast(const(char) *) A.data[i];
        if (sa == null || *sa == '\0') {
            continue;
        }
        if (strcmp(sa, sb) == 0) {
            copy[j] = 0;
            break;
        }
        }
        n += copy[j];
    }
    C = gretl_array_new(A.type, n, err);
    }

    if (!*err) {
    *err = gretl_array_copy_content(C, A, 0);
    }

    if (!*err && n > A.n) {
    i = A.n;
    for (j=0; j<B.n; j++) {
        if (copy[j]) {
        C.data[i++] = gretl_strdup(cast(const(char) *) B.data[j]);
        }
    }
    }
    free(copy);

    if (*err && C != null) {
    gretl_array_destroy(C);
    C = null;
    }

    return C;
}

/* respond to C = A && B */

GretlArray * gretl_arrays_intersection(GretlArray * A, GretlArray * B, int * err)
{
    GretlArray * C = null;
    const(char) * sa, sb;
    char * copy = null;
    int i, j, n = 0;

    if (A == null || B == null) {
    *err = GretlError.data;
    } else if (A.type != GretlType.strings ||
           B.type != GretlType.strings) {
    *err = GretlError.types;
    } else {
    if (A.n > 0) {
        copy = cast(typeof(copy)) calloc(1, A.n);
    }
    for (i=0; i<A.n; i++) {
        sa = cast(const(char) *) A.data[i];
        if (sa == null || *sa == '\0') {
        continue;
        }
        for (j=0; j<B.n; j++) {
        sb = cast(const(char) *) B.data[j];
        if (sb == null || *sb == '\0') {
            continue;
        }
        if (strcmp(sa, sb) == 0) {
            copy[i] = 1;
            break;
        }
        }
        n += copy[i];
    }
    C = gretl_array_new(A.type, n, err);
    }

    if (!*err && n > 0) {
    j = 0;
    for (i=0; i<A.n; i++) {
        if (copy[i]) {
        C.data[j++] = gretl_strdup(cast(const(char) *) A.data[i]);
        }
    }
    }
    free(copy);

    if (*err && C != null) {
    gretl_array_destroy(C);
    C = null;
    }

    return C;
}

/* respond to C = A ~ B, for strings arrays only so far */

GretlArray * gretl_arrays_concat(const(GretlArray) * A, const(GretlArray) * B, int * err)
{
    GretlArray * C = null;

    if (A == null || A.type != GretlType.strings ||
        B == null || B.type != GretlType.strings ||
        A.n != B.n) {
        *err = GretlError.types;
        return null;
    }

    C = gretl_array_new(GretlType.strings, A.n, err);

    if (C != null) {
        char * si;
        const(char) * sa;
        const(char) * sb;
        int i;

        for (i=0; i<A.n; i++) {
            sa = (A.data[i] == null)? "".ptr : cast(const(char) *) A.data[i];
            sb = (B.data[i] == null)? "".ptr : cast(const(char) *) B.data[i];
            si = cast(typeof(si)) calloc(strlen(sa) + strlen(sb) + 1, 1);
            sprintf(si, "%s%s".ptr, sa, sb);
            C.data[i] = si;
        }
    }

    return C;
}

/**
 * gretl_array_copy_as:
 * @name: name of source array.
 * @copyname: name for copy.
 * @copytpe: the type specified for the copied array, or 0.
 *
 * Look for a saved array specified by @name, and if found,
 * make a full copy and save it under @copyname. This is
 * called from geneval.c on completion of assignment to a
 * array named @copyname, where the returned value on the
 * right-hand side is a pre-existing array.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_array_copy_as(const(char) * name, const(char) * copyname, GretlType copytype)
{
    GretlArray * A0, A1 = null;
    UserVar * u;
    int err = 0;

    u = get_user_var_of_type_by_name(name, GretlType.array);
    if (u == null) {
    return GretlError.data;
    } else {
    A0 = cast(GretlArray *) user_var_get_value(u);
    }

    if (copytype > 0 && A0.type != copytype) {
    return GretlError.types;
    }

    /* is there a pre-existing array named @copyname? */
    u = get_user_var_of_type_by_name(copyname, GretlType.array);
    if (u != null) {
    A1 = cast(GretlArray *) user_var_get_value(u);
    }

    if (A1 != null) {
    if (A1.type != A0.type) {
        err = GretlError.types;
    } else {
        gretl_array_void_content(A1);
        A1.n = A0.n;
        err = array_allocate_storage(A1);
        if (!err) {
        err = gretl_array_copy_content(A1, A0, 0);
        }
    }
    } else {
    A1 = gretl_array_copy(A0, &err);
    if (!err) {
        err = user_var_add(copyname, A1.type, A1);
    }
    }

    return err;
}

/**
 * get_array_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved array, if found, else NULL.
 */

GretlArray * get_array_by_name(const(char) * name)
{
    GretlArray * a = null;

    if (name != null && *name != '\0') {
    UserVar * u = get_user_var_of_type_by_name(name, GretlType.array);

    if (u != null) {
        a = cast(GretlArray *) user_var_get_value(u);
    }
    }

    return a;
}

/**
 * get_strings_array_by_name:
 * @name: the name to look up.
 *
 * Returns: pointer to a saved array of strings, if found, else NULL.
 */

GretlArray * get_strings_array_by_name(const(char) * name)
{
    GretlArray * ret = null;

    if (name != null && *name != '\0') {
    UserVar * u = get_user_var_of_type_by_name(name, GretlType.array);

    if (u != null) {
        ret = cast(GretlArray *) user_var_get_value(u);
        if (ret.type != GretlType.strings) {
        ret = null;
        }
    }
    }

    return ret;
}

/**
 * get_strings_array_from_series:
 * @dset: gretl dataset.
 * @v: series ID number.
 * @err: location to receive error code.
 *
 * Returns: a newly allocated array of strings holding the string
 * values of series @v over the current sample range, or NULL on
 * failure.
 */

GretlArray * get_strings_array_from_series(DatasetRaw * dset, int v, int * err)
{
    GretlArray * ret = null;
    const(char) * st;
    int i, t, n;

    if (!is_string_valued(dset, v)) {
    *err = GretlError.types;
    return null;
    }

    n = dset.t2 - dset.t1 + 1;
    ret = gretl_array_new(GretlType.strings, n, err);
    if (ret == null) {
    return null;
    }

    for (t=dset.t1, i=0; t<=dset.t2; t++, i++) {
    st = series_get_string_for_obs(dset, v, t);
    ret.data[i] = gretl_strdup(st);
    }

    return ret;
}

GretlArray * gretl_array_pull_from_stack(const(char) * name, int * err)
{
    GretlArray * a = null;
    UserVar * u;

    u = get_user_var_of_type_by_name(name, GretlType.array);

    if (u != null) {
    a = cast(GretlArray *) user_var_unstack_value(u);
    }

    if (a == null) {
    *err = GretlError.data;
    }

    return a;
}

void print_array_string(const(char) * s, Prn * prn)
{
    int n = cast(int) strcspn(s, "\r\n");
    int m = cast(int) (cast(int) strlen(s));

    if (n > 72) {
    pprintf(prn, "\"%.69s...\"\n".ptr, s);
    } else if (n < m) {
    pprintf(prn, "\"%.*s...\"\n".ptr, n, s);
    } else {
    pprintf(prn, "\"%s\"\n".ptr, s);
    }
}

void print_array_elements(GretlArray * A, int imin, int imax, int range, Prn * prn)
{
    int i, lim = min(A.n, imax);

    for (i=imin; i<lim; i++) {
    pprintf(prn, "[%d] ".ptr, i+1);
    if (A.data[i] == null) {
        pputs(prn, "null\n");
    } else if (A.type == GretlType.strings) {
        const(char) * s = cast(const(char) *) A.data[i];

        print_array_string(s, prn);
    } else if (A.type == GretlType.matrices) {
        const(GretlMatrixRaw) * m = cast(const(GretlMatrixRaw) *) A.data[i];

        pprintf(prn, "%d x %d\n".ptr, m.rows, m.cols);
    } else if (A.type == GretlType.lists) {
        const(int) * list = cast(const(int) *) A.data[i];

        gretl_list_print(list, null, prn);
    }
    }

    if (!range && A.n > lim) {
    pputs(prn, "...\n\n");
    } else {
    pputc(prn, '\n');
    }
}

int gretl_array_print(GretlArray * A, Prn * prn)
{
    if (A != null) {
    const(char) * s = gretl_type_get_name(A.type);
    int nmax = 10;

    pprintf(prn, _("Array of %s, length %d\n").ptr, s, A.n);

    if (A.n > 0 &&
        A.type != GretlType.bundles && A.type != GretlType.arrays) {
        print_array_elements(A, 0, nmax, 0, prn);
    }
    }

    return 0;
}

int gretl_array_print_range(GretlArray * A, int imin, int imax, Prn * prn)
{
    if (A != null) {
    const(char) * s = gretl_type_get_name(A.type);

    pprintf(prn, _("Array of %s, length %d\n").ptr, s, A.n);

    if (A.type != GretlType.bundles && A.type != GretlType.arrays) {
        print_array_elements(A, imin, imax, 1, prn);
    }
    }

    return 0;
}

/* Called from gretl_bundle.c when serializing a bundle
   which contains one or more arrays.
*/

int gretl_array_serialize(GretlArray * A, Prn * prn)
{
    GretlType type;
    const(char) * subname;
    void * ptr;
    int i;
    int err = 0;

    type = gretl_type_get_singular(A.type);
    subname = gretl_type_get_name(type);

    pprintf(prn, "<gretl-array type=\"%s\" length=\"%d\">\n".ptr,
        gretl_type_get_name(A.type), A.n);

    for (i=0; i<A.n && !err; i++) {
    ptr = A.data[i];
    if (ptr == null) {
        pprintf(prn, "<%s placeholder=\"true\"/>\n".ptr, subname);
    } else if (type == GretlType.string_) {
        gretl_xml_put_tagged_string("string", cast(const(char) *) ptr, prn);
    } else if (type == GretlType.matrix) {
        err = gretl_matrix_serialize(cast(const(GretlMatrixRaw) *) ptr, null, prn);
    } else if (type == GretlType.bundle) {
        err = gretl_bundle_serialize(cast(GretlBundle *) ptr, null, prn);
    } else if (type == GretlType.array) {
        err = gretl_array_serialize(cast(GretlArray *) ptr, prn);
    } else if (type == GretlType.list) {
        gretl_xml_put_tagged_list("list", cast(const(int) *) ptr, prn);
    }
    }

    pputs(prn, "</gretl-array>\n");

    return err;
}

int name_matches_array_type(char * s, GretlType type)
{
    if (!strncmp(s, "gretl-", 6)) {
    /* we're expecting "short-form" type strings */
    s += 6;
    }

    return (gretl_type_from_string(s) == type);
}

int deserialize_array_elements(GretlArray * A, xmlNodePtr cur, xmlDocPtr doc)
{
    GretlType type = gretl_type_get_singular(A.type);
    int i = 0, err = 0;

    while (cur != null && !err && i < A.n) {
    if (!name_matches_array_type(cast(char *) cur.name, type)) {
        fprintf(stderr, "deserialize array: mismatched element '%s'\n".ptr,
            cast(char *) cur.name);
        err = GretlError.data;
    } else if (gretl_xml_get_prop_as_bool(cur, "placeholder")) {
        ; /* null array element: no-op */
    } else if (A.type == GretlType.strings) {
        A.data[i] = gretl_xml_get_string(cur, doc);
    } else if (A.type == GretlType.matrices) {
        A.data[i] = gretl_xml_get_matrix(cur, doc, &err);
    } else if (A.type == GretlType.bundles) {
        A.data[i] = gretl_bundle_deserialize(cur, doc, &err);
    } else if (A.type == GretlType.arrays) {
        A.data[i] = gretl_array_deserialize(cur, doc, &err);
    } else if (A.type == GretlType.lists) {
        A.data[i] = gretl_xml_get_list(cur, doc, &err);
    }
    /* note: arrays of scalars not handled */
    i++;
    cur = cur.next;
    }

    if (!err && i != A.n) {
    fprintf(stderr, "deserialize array: array is corrupted\n".ptr);
    err = GretlError.data;
    }

    return err;
}

/* For internal use only: @p1 should be of type xmlNodePtr and @p2
   should be an xmlDocPtr. We suppress the actual pointer types in the
   prototype so that it's possible for a module to include
   gretl_array.h without including the full libxml headers.
*/

GretlArray * gretl_array_deserialize(void * p1, void * p2, int * err)
{
    xmlNodePtr node = cast(xmlNode *) p1;
    xmlDocPtr doc = cast(xmlDoc *) p2;
    GretlType type = cast(GretlType) 0;
    int n = 0;
    GretlArray * A = null;

    if (xmlStrcmp(node.name, cast(XUC) "gretl-array")) {
    fprintf(stderr, "deserialize array: node is not gretl-array!\n".ptr);
    *err = GretlError.data;
    } else {
    type = gretl_xml_get_type_property(node);
    if (type == 0) {
        fprintf(stderr, "deserialize array: couldn't get array type\n".ptr);
        *err = GretlError.data;
    }
    }

    if (!*err && !gretl_xml_get_prop_as_int(node, "length", &n)) {
    fprintf(stderr, "deserialize array: couldn't get length\n".ptr);
    *err = GretlError.data;
    }

    if (!*err) {
    A = gretl_array_new(type, n, err);
    }

    if (A != null && n > 0) {
    *err = deserialize_array_elements(A, node.children, doc);
    if (*err) {
        gretl_array_destroy(A);
        A = null;
    }
    }

    return A;
}

/* In case a matrix is too wide for comfortable printing, split it by
   column. If there are one or more "leading" columns that should be
   displayed with each chunk of the matrix this should be signalled
   via a non-zero value for @leadcols. The maximum number of columns
   per chunk (including the leading columns, if any) is set by the
   @maxcols argument.

   The return value is an array of suitably sized matrices. If @m
   has column names attached these are distributed to the matrices
   in the array.
*/

GretlArray * gretl_matrix_col_split(const(GretlMatrixRaw) * m, int leadcols, int maxcols, int * err)
{
    GretlArray * a = null;
    int maincols, nm, rem;

    if (gretl_is_null_matrix(m)) {
    *err = GretlError.invarg;
    } else {
    maincols = m.cols - leadcols;
    if (maincols <= 0) {
        *err = GretlError.invarg;
    }
    }

    if (!*err) {
    /* how many matrices will we need? */
    nm = maincols / (maxcols - leadcols);
    rem = maincols % (maxcols - leadcols);
    nm += rem > 0;
    if (nm == 1) {
        /* nothing to be done */
        *err = GretlError.invarg;
    }
    }

    if (!*err) {
    a = gretl_array_new(GretlType.matrices, nm, err);
    }

    if (!*err) {
    const(char) ** S0 = gretl_matrix_get_colnames(m);
    char ** Si = null;
    size_t rsize = m.rows * double.sizeof;
    const(double) * src;
    double * targ;
    GretlMatrixRaw * ai;
    int i, j, cols, spos;

    /* initial read position for "main" columns */
    src = m.val + leadcols * m.rows;
    spos = leadcols;

    for (i=0; i<nm && !*err; i++) {
        cols = (i == nm-1 && rem > 0)? rem + leadcols : maxcols;
        ai = gretl_zero_matrix_new(m.rows, cols);
        if (ai == null) {
        *err = GretlError.alloc;
        } else {
        Si = S0 == null ? null : strings_array_new(cols);
        /* initial write position */
        targ = ai.val;
        if (leadcols > 0) {
            memcpy(targ, m.val, leadcols * rsize);
            targ += leadcols * m.rows;
            if (Si != null) {
            /* transcribe column names */
            for (j=0; j<leadcols; j++) {
                Si[j] = gretl_strdup(S0[j]);
            }
            }
        }
        /* now handle the "main" columns */
        cols -= leadcols;
        memcpy(targ, src, cols * rsize);
        /* advance the read position */
        src += cols * m.rows;
        if (Si != null) {
            for (j=0; j<cols; j++) {
            Si[leadcols+j] = gretl_strdup(S0[spos++]);
            }
            gretl_matrix_set_colnames(ai, Si);
        }
        /* stick matrix @i into the array */
        gretl_array_set_matrix(a, i, ai, 0);
        }
    }
    }

    if (*err && a != null) {
    gretl_array_destroy(a);
    a = null;
    }

    return a;
}

/* is_strings_array_element() tests @str for representation of an
   element of an array of strings. This should work if @str takes the
   form

   <arrayname>[<index>]

   where @arrayname identifies a strings array and <index> represents
   a valid (1-based) index into this array. In that case @arrayname is
   written into @aname (which should be of length 32 bytes or more),
   @index is written into @idx (should be 8 bytes or more), and 1 is
   returned. On failure, 0 is returned.
*/

int is_strings_array_element(const(char) * str, char * aname, char * idx)
{
    int ret = 0;

    if (strchr(str, '[') != null) {
    GretlArray * A = null;
    int i, err = 0;

    if (sscanf(str, "%31[^[][%7[^]]".ptr, aname, idx) == 2) {
        A = get_strings_array_by_name(aname);
    }
    if (A != null) {
            if (!strcmp(idx, "end")) {
                i = A.n;
                sprintf(idx, "%d".ptr, i);
            } else {
                i = generate_int(idx, null, &err);
            }
        if (!err && i > 0 && i <= A.n) {
        ret = 1;
        }
    }
    }

    return ret;
}

/* start apparatus for sorting a gretl array via qsort */

struct Asorter {
    Fncall * call;
    GretlType type;
    DatasetRaw * dset;
    Prn * prn;
    int err;
}

__gshared Asorter asort;

void asort_cleanup()
{
    fncall_destroy(asort.call);
    asort.call = null;
    asort.type = cast(GretlType) 0;
    asort.dset = null;
    asort.prn = null;
    asort.err = 0;
}

int asort_setup(const(char) * fname, GretlType type, DatasetRaw * dset, Prn * prn)
{
    Fncall * call = null;
    GretlType rtype;
    Ufunc * func;
    int i, err = 0;

    func = get_user_function_by_name(fname);
    if (func == null) {
    return GretlError.data;
    }

    call = fncall_new(func, 1);
    if (call == null) {
    return GretlError.alloc;
    }

    rtype = fncall_get_return_type(call);
    if (rtype != GretlType.double_) {
    fncall_destroy(call);
    return GretlError.types;
    }

    for (i=0; i<2 && !err; i++) {
    err = push_anon_function_arg(call, type, null);
    }

    if (err) {
    fncall_destroy(call);
    } else {
    asort.call = call;
    asort.type = type;
    asort.dset = dset;
    asort.prn = prn;
    asort.err = 0;
    }

    return err;
}

extern (C) int compare_array_elements(const(void) * a, const(void) * b)
{
    void * va = *cast(void **) a;
    void * vb = *cast(void **) b;
    double d = 0;
    double * pd = &d;

    if (asort.err) {
    /* get out asap */
    return 0;
    }

    set_anon_function_arg(asort.call, 0, asort.type, va);
    set_anon_function_arg(asort.call, 1, asort.type, vb);
    asort.err = gretl_function_exec(asort.call, GretlType.double_,
                    asort.dset, &pd, asort.prn);
    if (asort.err) {
    fprintf(stderr, " compare_array_elements: err=%d, d=%g\n".ptr, asort.err, *pd);
        return 0;
    } else {
        d = *pd;
        return d > 0 ? 1 : d < 0 ? -1 : 0;
    }
}

int gretl_array_qsort(GretlArray * a, const(char) * fname, DatasetRaw * dset, Prn * prn)
{
    int n = gretl_array_get_length(a);
    GretlType atype;
    void * aptr;
    int err;

    if (n < 2) {
    /* nothing to be done */
    return 0;
    }

    atype = gretl_array_get_content_type(a);
    err = asort_setup(fname, atype, dset, prn);
    if (err) {
    return err;
    }

    aptr = gretl_array_get_all_data(a);

    set_user_qsorting(1);
    qsort(aptr, n, (void *).sizeof, &compare_array_elements);
    set_user_qsorting(0);

    err = asort.err;
    asort_cleanup();

    return err;
}

/* end apparatus for sorting a gretl array via qsort */

/* Check for equality between two arrays. Returns 1 if the two arrays
   are of the same type and length, and all their elements compare
   equal, else returns 0.
*/

int gretl_arrays_are_equal(const(GretlArray) * a, const(GretlArray) * b)
{
    if (a == b) {
    return 1;
    } else if (a.type != b.type || a.n != b.n) {
    return 0;
    } else {
    int i, eq, nulls;
    int err = 0;

    for (i=0; i<a.n; i++) {
        if (a.type == GretlType.strings) {
        nulls = (a.data[i] == null) + (b.data[i] == null);
        if (nulls == 1) {
            return 0;
        } else if (nulls == 0 && strcmp(cast(const(char) *) a.data[i], cast(const(char) *) b.data[i])) {
            return 0;
        }
        } else if (a.type == GretlType.matrices) {
        eq = gretl_matrices_are_equal(cast(const(GretlMatrixRaw) *) a.data[i], cast(const(GretlMatrixRaw) *) b.data[i], 0, &err);
        if (eq != 1) {
            /* could be 0 or non-conformable (-1) */
            return eq; /* or return 0? */
        }
        } else if (a.type == GretlType.bundles) {
        if (!gretl_bundles_are_equal(cast(GretlBundle *) a.data[i], cast(GretlBundle *) b.data[i])) {
            return 0;
        }
        } else if (a.type == GretlType.arrays) {
        if (!gretl_arrays_are_equal(cast(const(GretlArray) *) a.data[i], cast(const(GretlArray) *) b.data[i])) {
            return 0;
        }
        } else if (a.type == GretlType.lists) {
        if (gretl_list_cmp(cast(const(int) *) a.data[i], cast(const(int) *) b.data[i])) {
            return 0;
        }
        }
    }
    return 1;
    }
}

/* arglist_validate() is used only in geneval.c but it's defined here
   because (a) geneval.c is already overburdened and (b) the job can
   be done more efficiently with access to the internals of
   gretl_array. We're checking that @args is an array of strings
   holding exactly the strings in @keys minus "arglist" itself.
*/

int arglist_validate(GretlArray * keys, GretlArray * args)
{
    const(char) * sk, sa;
    int n_keys = keys.n - 1;
    int i, j, found = 0;

    if (args.type != GretlType.strings || args.n != n_keys) {
    return 0;
    }

    for (i=0; i<keys.n; i++) {
    sk = cast(const(char) *) keys.data[i];
    for (j=0; j<args.n; j++) {
        sa = cast(const(char) *) args.data[j];
        if (sa == null || !strcmp(sa, "arglist")) {
        return 0;
        } else if (!strcmp(sa, sk)) {
        found++;
        break;
        }
    }
    }

    return found == n_keys;
}

