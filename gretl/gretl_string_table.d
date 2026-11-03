/* Port of gretl_string_table.c (gretl_string_table.h has only prototypes: no header module)

   The hash tables of series tables are the AA-backed ghashtable.d (not glib); USE_CURL (web resources)
   is not defined, so retrieve_file_content() reports "Internet access not supported" for URLs.
   The two structs defined here (SeriesTable, GretlStringTable) were opaque declarations in
   dataseth.d and libgretlh.d, which now import them from this module. */

module gretl.gretl_string_table;


import core.stdc.errno : errno;
import core.stdc.math : lrint;
import core.stdc.stdio : fprintf, stderr;
import core.stdc.stdlib : atof, atoi, free, getenv, malloc, realloc, strtod;
import core.stdc.string : strchr, strcmp, strcpy, strlen, strncat, strncmp, strspn;
import etc.c.zlib : gzFile, gzclose, gzgetc, gzread, gzrewind;

import gretl.dataioh : GretlFileType;
import gretl.dataseth : VarFlags;
import gretl.dataset : series_attach_string_table, series_get_string_table, series_set_flag, series_unset_flag;
import gretl.dataio : is_gzipped, ntolabel;
import gretl.gerror : GError, g_error_free;
import gretl.ghashtable : GHashTable, g_hash_table_destroy, g_hash_table_insert, g_hash_table_lookup,
    g_hash_table_new, g_str_equal, g_str_hash;
import gretl.glib : gchar, gpointer, gsize, gssize;
import gretl.gfileutils : g_file_get_contents;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_new,
    in_gretl_list;
import gretl.gretl_matrixh : GretlMatrixRaw, gretl_vector_get_length;
import gretl.gretl_paths : gretl_addpath, gretl_gzopen;
import gretl.gretl_prn : Prn, gretl_print_steal_buffer, pprintf, pputc;
import gretl.gstrfuncs : g_convert, g_get_charset, g_locale_to_utf8, g_strdup, g_utf8_validate;
import gretl.gretl_utils : gretl_int_from_double, gretl_isdummy;
import gretl.libgretlh;
import gretl.libset : gretl_in_gui_mode, libset_get_bool;
import gretl.libseth : SetKey;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_namechar_spn, gretl_strdup, gretl_strndup, integer_string, strings_array_add,
    strings_array_dup, strings_array_free, strings_array_realloc_with_length;
import gretl.strutilsh : SLASH;
import gretl.gretl_array : gretl_array_destroy, gretl_array_new, gretl_array_set_data;
import gretl.uservar : get_string_by_name;

private enum StFlags {
    quoted  = 1 << 0,
    allints = 1 << 1,
    alldbls = 1 << 2
}

struct SeriesTable {
    int n_strs;         /* number of strings in table */
    char ** strs;       /* saved strings */
    GHashTable * ht;    /* hash table for quick lookup */
    int flags;          /* status flags (above) */
}

struct GretlStringTable {
    int * cols_list;        /* list of included columns */
    SeriesTable ** cols;    /* per-column tables (see above) */
    char * extra;           /* extra information, if any */
}

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void gretl_shell_grab(const(char) * arg, char ** sout);
}

private bool st_quoted(const(SeriesTable) * t) { return (t.flags & StFlags.quoted) != 0; }
private bool all_ints(const(SeriesTable) * t)  { return (t.flags & StFlags.allints) != 0; }
private bool all_dbls(const(SeriesTable) * t)  { return (t.flags & StFlags.alldbls) != 0; }
private bool all_num(const(SeriesTable) * t)   { return (t.flags & (StFlags.allints | StFlags.alldbls)) != 0; }

private gpointer index_to_ptr(int n) { return cast(gpointer) cast(size_t) n; }
private int ptr_to_index(gpointer p) { return cast(int) cast(size_t) p; }

private SeriesTable * series_table_alloc()
{
    SeriesTable * st = cast(SeriesTable *) malloc(SeriesTable.sizeof);

    if (st != null) {
        st.strs = null;
        st.n_strs = 0;
        st.ht = g_hash_table_new(&g_str_hash, &g_str_equal);
        st.flags = 0;
    }

    return st;
}

private GretlStringTable * gretl_string_table_alloc()
{
    GretlStringTable * gst = cast(GretlStringTable *) malloc(GretlStringTable.sizeof);

    if (gst != null) {
        gst.cols_list = null;
        gst.cols = null;
        gst.extra = null;
    }

    return gst;
}

/**
 * gretl_string_table_new:
 * @list: list of series IDs whose values are to be
 * given a string representation, or NULL.
 *
 * These values in @list should correspond to the 0-based indices
 * of the series in question within the dataset.  For example,
 * if strings are to be recorded for variables 2, 5 and 10 the
 * @list argument would be {3, 2, 5, 10}. If NULL is passed for
 * @list the return value is an initialized, empty string table.
 *
 * Returns: pointer to a newly allocated string table or NULL
 * on failure.
 */

GretlStringTable * gretl_string_table_new(const(int) * list)
{
    GretlStringTable * gst;
    int ncols = 0;
    int err = 0;

    gst = gretl_string_table_alloc();
    if (gst == null) {
        return null;
    }

    if (list != null && list[0] > 0) {
        gst.cols_list = gretl_list_copy(list);
        if (gst.cols_list == null) {
            err = GretlError.alloc;
        } else {
            ncols = list[0];
        }
    }

    if (ncols > 0) {
        gst.cols = cast(SeriesTable **) malloc(ncols * (*gst.cols).sizeof);
        if (gst.cols == null) {
            err = GretlError.alloc;
        } else {
            int i, j;

            for (i=0; i<ncols && !err; i++) {
                gst.cols[i] = series_table_alloc();
                if (gst.cols[i] == null) {
                    err = GretlError.alloc;
                    for (j=0; j<i; j++) {
                        free(gst.cols[j]);
                    }
                    free(gst.cols);
                }
            }
        }
    }

    if (err) {
        free(gst.cols_list);
        free(gst);
        gst = null;
    }

    return gst;
}

private int series_table_get_index(const(SeriesTable) * st,
                                   const(char) * s)
{
    int ret = 0;

    if (st != null && st.ht != null && s != null) {
        gpointer p = g_hash_table_lookup(cast(GHashTable *) st.ht, s);

        if (p != null) {
            ret = ptr_to_index(p);
        }
    }

    return ret;
}

/**
 * series_table_get_value:
 * @st: a gretl series table.
 * @s: the string to look up.
 *
 * Returns: the numerical value associated with @s in the
 * given series table, or #NADBL in case there is no match.
 */

double series_table_get_value(const(SeriesTable) * st, const(char) * s)
{
    int k = series_table_get_index(st, s);

    return (k > 0)? cast(double) k : NADBL;
}

/**
 * series_table_get_string:
 * @st: a gretl series table.
 * @val: the numerical value to look up.
 *
 * Returns: the string associated with @val in the
 * given series table, or NULL in case there is no match.
 */

const(char) * series_table_get_string(const(SeriesTable) * st, double val)
{
    const(char) * ret = null;

    if (!na(val)) {
        int k = cast(int) lrint(val);

        if (k > 0 && k <= st.n_strs) {
            ret = st.strs[k-1];
        }
    }

    return ret;
}

/**
 * series_table_map:
 * @st_from: gretl series table.
 * @st_to: gretl series table.
 *
 * Constructs a mapping from the integer codes in @st_from
 * to those in @st_to. For example, if the string "foo"
 * has code 3 in @st_from and code 12 in @st_to, then
 * element 3 in the mapping array will have value 12.
 * For any strings in @st_from that are not matched
 * in @st_to, the associated element of the map is set
 * to -1.
 *
 * Element 0 of the map holds the number of following
 * elements, which is the same as the number of strings in
 * @st_from.
 *
 * Returns: allocated array of int or NULL in case of failure.
 */

int * series_table_map(const(SeriesTable) * st_from, const(SeriesTable) * st_to)
{
    int * map = null;
    int n = st_from.n_strs;

    map = gretl_list_new(n);

    if (map != null) {
        const(char) * s1;
        int i, i2;

        for (i=0; i<n; i++) {
            s1 = st_from.strs[i];
            i2 = series_table_get_index(st_to, s1);
            map[i+1] = i2 == 0 ? -1 : i2;
        }
    }

    return map;
}

/**
 * series_table_get_strings:
 * @st: a gretl series table.
 * @n_strs: location to receive the number of strings, or NULL.
 *
 * Returns: the array of strings associated with @st. These
 * should not be modified in any way.
 */

char ** series_table_get_strings(const(SeriesTable) * st, int * n_strs)
{
    if (st != null) {
        if (n_strs != null) {
            *n_strs = st.n_strs;
        }
        return cast(char **) st.strs;
    } else {
        return null;
    }
}

int series_table_get_n_strings(const(SeriesTable) * st)
{
    if (st != null) {
        return st.n_strs;
    } else {
        return 0;
    }
}

private char * get_unquoted(const(char) * s)
{
    char * tmp = null;
    int n = cast(int) strlen(s);

    if (s[n-1] == s[0]) {
        tmp = gretl_strndup(s+1, n-2);
    }

    return tmp;
}

/**
 * series_table_add_string:
 * @st: a gretl series table.
 * @s: new string to add.
 *
 * Returns: the index of the new string within the table, or
 * -1 on failure.
 */

int series_table_add_string(SeriesTable * st, const(char) * s)
{
    char * tmp = null;
    int n, err;

    if (s == null) {
        return -1;
    }

    if (*s == '"' || *s == '\'') {
        tmp = get_unquoted(s);
    }

    if (tmp != null) {
        st.flags |= StFlags.quoted;
        err = strings_array_add(&st.strs, &st.n_strs, tmp);
        free(tmp);
    } else {
        err = strings_array_add(&st.strs, &st.n_strs, s);
    }

    if (err) {
        n = -1;
    } else {
        n = st.n_strs;
        g_hash_table_insert(st.ht, cast(gpointer) st.strs[n-1],
                            index_to_ptr(n));
    }

    return n;
}

/**
 * series_table_add_strings:
 * @st: a gretl series table.
 * @S: array of new strings to add.
 * @ns: number of elements of @S.
 *
 * Appends the strings in @S to @st.
 *
 * Returns: 0 on successful completion, or error code on error.
 */

int series_table_add_strings(SeriesTable * st,
                             const(char) ** S,
                             int ns)
{
    int oldn, newn;
    int i, j, err = 0;

    if (S == null || ns <= 0) {
        return 0;
    }

    oldn = st.n_strs;
    newn = oldn + ns;
    st.strs = strings_array_realloc_with_length(&st.strs,
                                                oldn, newn,
                                                0);
    if (st.strs == null) {
        err = GretlError.alloc;
    } else {
        st.n_strs = newn;
        for (i=0, j=oldn; i<ns; i++, j++) {
            st.strs[j] = gretl_strdup(S[i]);
            g_hash_table_insert(st.ht, cast(gpointer) st.strs[j],
                                index_to_ptr(j+1));
        }
    }

    return err;
}

SeriesTable * series_table_new(char ** strs, int n_strs, int * err)
{
    SeriesTable * st = series_table_alloc();
    int i;

    if (st == null || strs == null) {
        *err = GretlError.alloc;
    } else {
        st.n_strs = n_strs;
        st.strs = strs;
        for (i=0; i<n_strs; i++) {
            if (st.strs[i] == null) {
                fprintf(stderr, "series_table_new: str %d is NULL\n", i);
                *err = GretlError.data;
            } else {
                g_hash_table_insert(st.ht, cast(gpointer) st.strs[i],
                                    index_to_ptr(i+1));
            }
        }
    }

    return st;
}

SeriesTable * series_table_copy(const(SeriesTable) * st)
{
    SeriesTable * ret = null;

    if (st != null) {
        ret = series_table_alloc();
    }

    if (ret != null) {
        char ** S = strings_array_dup(cast(const(char *) *) st.strs, st.n_strs);
        int i;

        if (S == null) {
            series_table_destroy(ret);
            ret = null;
        } else {
            ret.n_strs = st.n_strs;
            ret.strs = S;
            for (i=0; i<ret.n_strs; i++) {
                g_hash_table_insert(ret.ht, cast(gpointer) ret.strs[i],
                                    index_to_ptr(i+1));
            }
        }
    }

    return ret;
}

private SeriesTable * gretl_string_table_add_column(GretlStringTable * gst, int colnum)
{
    SeriesTable ** cols;
    int * newlist;
    int n, err = 0;

    newlist = gretl_list_append_term(&gst.cols_list, colnum);
    if (newlist == null) {
        return null;
    }

    n = gst.cols_list[0];

    cols = cast(SeriesTable **) realloc(gst.cols, n * (*cols).sizeof);
    if (cols == null) {
        err = GretlError.alloc;
    } else {
        gst.cols = cols;
        cols[n-1] = series_table_alloc();
        if (cols[n-1] == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        gst.cols_list[0] -= 1;
        return null;
    } else {
        return cols[n-1];
    }
}

/**
 * gretl_string_table_index:
 * @gst: a gretl string table.
 * @s: the string to look up or add.
 * @col: index of the column to be accessed or created.
 * @addcol: non-zero to indicate that a column should be
 * added if it's not already present.
 * @prn: gretl printer (or %NULL).
 *
 * This function has two main uses: for lookup in the context of
 * a completed string table, or for constructing such a
 * table (with @addcol non-zero).  The returned index reflects
 * any additions to the table that may be required (if column
 * @col does not already exist, or if string @s is not already
 * stored for column @col).
 *
 * Returns: the 1-based index of @s within the column of
 * @st that has index @col, if available, otherwise 0.
 */

int gretl_string_table_index(GretlStringTable * gst, const(char) * s,
                             int col, int addcol, Prn * prn)
{
    SeriesTable * st = null;
    char * tmp = null;
    int i, idx = 0;

    if (gst == null) {
        return idx;
    }

    if (gst.cols_list != null) {
        for (i=1; i<=gst.cols_list[0]; i++) {
            if (gst.cols_list[i] == col) {
                st = gst.cols[i-1];
                break;
            }
        }
    }

    if (*s == '"') {
        tmp = get_unquoted(s);
    }

    if (st != null) {
        /* there's a table for this column already */
        idx = series_table_get_index(st, tmp != null ? tmp : s);
    } else if (addcol) {
        /* no table for this column yet: start one now */
        st = gretl_string_table_add_column(gst, col);
        if (st != null) {
            pprintf(prn, (_("variable %d: translating from strings to ")
                          ~ _("code numbers\n")).ptr, col);
        }
    }

    if (idx == 0 && st != null) {
        idx = series_table_add_string(st, s);
    }

    free(tmp);

    return idx;
}

/* Used in the context of deletion of "empty" variables from
   an imported dataset: the index of a given "column" in
   a string table is adjusted to match the new position of
   the variable in question.
*/

int gretl_string_table_reset_column_id(GretlStringTable * gst,
                                       int oldid, int newid)
{
    if (gst != null) {
        int i;

        for (i=1; i<=gst.cols_list[0]; i++) {
            if (gst.cols_list[i] == oldid) {
                gst.cols_list[i] = newid;
                return 0;
            }
        }
    }

    return GretlError.data;
}

SeriesTable * gretl_string_table_detach_col(GretlStringTable * gst,
                                            int col)
{
    SeriesTable * st = null;

    if (gst != null) {
        int pos = in_gretl_list(gst.cols_list, col);
        int i, n = gst.cols_list[0];

        if (pos > 0) {
            st = gst.cols[pos-1];
            for (i=pos-1; i<n-1; i++) {
                gst.cols[i] = gst.cols[i+1];
            }
            gst.cols[n-1] = null;
            gretl_list_delete_at_pos(gst.cols_list, pos);
        }
    }

    return st;
}

int in_string_table(GretlStringTable * gst, int id)
{
    if (gst != null) {
        return in_gretl_list(gst.cols_list, id);
    } else {
        return 0;
    }
}

int * string_table_copy_list(GretlStringTable * gst)
{
    if (gst != null) {
        return gretl_list_copy(gst.cols_list);
    } else {
        return null;
    }
}

int string_table_replace_list(GretlStringTable * gst,
                              int * newlist)
{
    if (gst != null) {
        /* FIXME pruning? */
        free(gst.cols_list);
        gst.cols_list = newlist;
    }

    return GretlError.data;
}

/**
 * series_table_destroy:
 * @st: series string table.
 *
 * Frees all resources associated with @st.
 */

void series_table_destroy(SeriesTable * st)
{
    if (st != null) {
        strings_array_free(st.strs, st.n_strs);
        if (st.ht != null) {
            g_hash_table_destroy(st.ht);
        }
        free(st);
    }
}

/**
 * series_table_free_shallow:
 * @st: series string table.
 *
 * Does a "shallow"free on @st, without freeing the
 * array of strings to which the table has a reference.
 */

void series_table_free_shallow(SeriesTable * st)
{
    if (st != null) {
        if (st.ht != null) {
            g_hash_table_destroy(st.ht);
        }
        free(st);
    }
}

/**
 * gretl_string_table_destroy:
 * @gst: gretl string table.
 *
 * Frees all resources associated with @gst.
 */

void gretl_string_table_destroy(GretlStringTable * gst)
{
    int i, ncols;

    if (gst == null) return;

    ncols = (gst.cols_list != null)? gst.cols_list[0] : 0;

    for (i=0; i<ncols; i++) {
        series_table_destroy(gst.cols[i]);
    }
    free(gst.cols);

    free(gst.cols_list);

    if (gst.extra != null) {
        free(gst.extra);
    }

    free(gst);
}

/* Given a series_table in which all the strings are just
   representations of integers, write the integer values
   into the series and destroy the table, while marking
   the series as "coded". Or, if all the strings are valid
   doubles, convert to numeric.
*/

private void series_commute_string_table(DatasetRaw * dset, int i,
                                         SeriesTable * st)
{
    if (dset != null && i > 0 && i < dset.v) {
        const(char) * s;
        double val;
        int t;

        for (t=0; t<dset.n; t++) {
            val = dset.Z[i][t];
            if (!na(val)) {
                s = series_table_get_string(st, val);
                if (all_ints(st)) {
                    dset.Z[i][t] = cast(double) atoi(s);
                } else {
                    dset.Z[i][t] = atof(s);
                }
            }
        }
        if (all_ints(st)) {
            series_set_flag(dset, i, VarFlags.discrete);
            if (!gretl_isdummy(0, dset.n - 1, dset.Z[i])) {
                series_set_flag(dset, i, VarFlags.coded);
            }
        } else {
            series_unset_flag(dset, i, VarFlags.discrete);
        }
        series_table_destroy(st);
    }
}

void series_table_print(DatasetRaw * dset, int i, Prn * prn)
{
    SeriesTable * st = series_get_string_table(dset, i);
    int j;

    if (st != null) {
        pprintf(prn, _("String code table for variable %d (%s):\n").ptr,
                i, dset.varname[i]);
        pputc(prn, '\n');
        for (j=0; j<st.n_strs; j++) {
            pprintf(prn, "%5d \"%s\"\n", j+1, st.strs[j]);
        }
        pputc(prn, '\n');
    }
}

/**
 * gretl_string_table_finalize:
 * @gst: gretl string table.
 * @dset: dataset.
 *
 * Attaches the string-value information in @gst to the dataset
 * @dset. However, if one or more of the series referenced by @gst
 * are deemed to be integer codes or misclassified numeric data,
 * their "series tables" are commuted into numeric form.
 *
 * Returns: the number of series tables actually attached.
 */

int gretl_string_table_finalize(GretlStringTable * gst, DatasetRaw * dset)
{
    SeriesTable * st;
    int i, vi, ret = 0;

    if (gst == null || gst.cols_list == null ||
        dset == null || dset.varinfo == null) {
        return 0;
    }

    for (i=0; i<gst.cols_list[0]; i++) {
        st = gst.cols[i];
        if (st != null) {
            vi = gst.cols_list[i+1];
            if (all_num(st)) {
                series_commute_string_table(dset, vi, st);
            } else {
                series_attach_string_table(dset, vi, st);
                ret++;
            }
            gst.cols[i] = null;
        }
    }

    return ret;
}

private int string_is_double(const(char) * s)
{
    const(char) * test;

    errno = 0;
    strtod(s, &test);
    return errno == 0 && *test == '\0';
}

/**
 * gretl_string_table_validate:
 * @gst: gretl string table.
 * @opt: may include OPT_S to indicate that the string
 * table was constructed from a spreadsheet file, in which
 * case certain cells may be explicitly marked as holding
 * strings, which information may be not be reliable.
 *
 * Checks that the "string values" in @gst are not in fact
 * undigested quasi-numerical values. We run this on
 * imported CSV data to ensure we don't produce
 * misleading results. In the spreadsheet case we're on
 * the lookout for purely numerical data wrongly identified
 * as strings; we don't treat that as fatal but mark the
 * series table(s) for conversion to numeric.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_string_table_validate(GretlStringTable * gst,
                                Gretlopt opt)
{
    const(char) * test = "0123456789.,";
    int ssheet = (opt & Gretlopt.OPT_S) != 0;
    int i, ncols = 0;
    int err = 0;

    if (gst != null && gst.cols_list != null) {
        ncols = gst.cols_list[0];
    }

    for (i=0; i<ncols; i++) {
        SeriesTable * st = gst.cols[i];
        const(char) * s;
        int nint = 0;
        int ndbl = 0;
        int j, myerr = GretlError.data;

        if (st_quoted(st)) {
            myerr = 0;
        }

        for (j=0; j<st.n_strs; j++) {
            s = st.strs[j];
            if (st_quoted(st)) {
                if (integer_string(s)) {
                    nint++;
                }
            } else {
                /* not quoted */
                if (ssheet && (*s == '\0' || string_is_double(s))) {
                    /* could really be numeric? (2020-12-25) */
                    ndbl++;
                    continue;
                } else if (*s == '-' || *s == '+') {
                    s++;
                }
                if (strspn(s, test) < strlen(s)) {
                    /* not quasi-numeric */
                    myerr = 0;
                    break;
                }
            }
        }

        if (nint == st.n_strs) {
            /* treat as integer codes */
            st.flags |= StFlags.allints;
        } else if (ndbl == st.n_strs) {
            /* all really numeric */
            st.flags |= StFlags.alldbls;
            myerr = 0;
        }

        if (myerr) {
            err = myerr;
            break;
        }
    }

    return err;
}

/**
 * gretl_string_table_save:
 * @gst: gretl string table.
 * @dset: dataset information (for names of variables).
 *
 * Attaches the content of @gst to @dset.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_string_table_save(GretlStringTable * gst, DatasetRaw * dset)
{
    SeriesTable * st;
    int i, vi, ncols = 0;

    if (gst == null || dset.varinfo == null) {
        return GretlError.data;
    }

    ncols = (gst.cols_list != null)? gst.cols_list[0] : 0;

    for (i=0; i<ncols; i++) {
        st = gst.cols[i];
        if (st != null) {
            vi = gst.cols_list[i+1];
            st = gst.cols[i];
            series_attach_string_table(dset, vi, st);
            gst.cols[i] = null;
        }
    }

    return 0;
}

/**
 * gretl_string_table_add_extra:
 * @gst: gretl string table.
 * @prn: gretl printer.
 *
 * Steals the printing buffer from @prn and adds it to @gst.
 * The buffer will be appended when @gst is printed via
 * gretl_string_table_finalize().
 */

void gretl_string_table_add_extra(GretlStringTable * gst, Prn * prn)
{
    if (gst != null && prn != null) {
        gst.extra = gretl_print_steal_buffer(prn);
    }
}

/* apparatus for built-in strings: the names, and (parallel array) the current values;
   the C has an array of structs with a char[VNAMELEN] name each */

private immutable string[12] built_in_names = [
    "gretldir", "dotdir", "workdir", "gnuplot", "x12a", "x12adir",
    "tramo", "tramodir", "seats", "pkgdir", "lang", "logfile"
];

private __gshared gchar *[12] built_in_vals;

void builtin_strings_cleanup()
{
    int i, n = cast(int) built_in_names.length;

    for (i=0; i<n; i++) {
        g_free(built_in_vals[i]);
    }
}

/**
 * gretl_insert_builtin_string:
 * @name: the name of the string to be added or replaced.
 * @s: the value for this string variable.
 *
 * Inserts value @s for string @name in gretl's table
 * of built-in string variables.
 */

void gretl_insert_builtin_string(const(char) * name, const(char) * s)
{
    int i, n = cast(int) built_in_names.length;
    int m, gui = gretl_in_gui_mode();

    for (i=0; i<n; i++) {
        if (!strcmp(name, built_in_names[i].ptr)) {
            g_free(built_in_vals[i]);
            if (s == null) {
                built_in_vals[i] = null;
            } else if (gui && !g_utf8_validate(s, -1, null)) {
                /* handle non-ASCII Windows paths */
                gsize bytes;
                gchar * u;

                u = g_locale_to_utf8(s, -1, null, &bytes, null);
                if (u != null) {
                    m = cast(int) strlen(u);
                    if (m > 0 && u[m-1] == SLASH) {
                        u[m-1] = '\0';
                    }
                }
                built_in_vals[i] = u;
            } else {
                m = cast(int) strlen(s);
                if (m > 0 && s[m-1] == SLASH) {
                    /* drop trailing dir separator for paths */
                    built_in_vals[i] = gretl_strndup(s, m - 1);
                } else {
                    built_in_vals[i] = g_strdup(s);
                }
            }
            return;
        }
    }
}

char * get_built_in_string_by_name(const(char) * name)
{
    int i, n = cast(int) built_in_names.length;

    for (i=0; i<n; i++) {
        if (!strcmp(name, built_in_names[i].ptr)) {
            char * s = built_in_vals[i];

            return s != null ? s : cast(char *) "".ptr;
        }
    }

    return null;
}

/* Try to recode the content of a local file or web resource
   to UTF-8. Can be tricky since we don't know the original
   encoding of the content.
*/

private gchar * recode_content(gchar * orig, const(char) * codeset,
                               int * err)
{
    const(gchar) * charset = null;
    GError * gerr = null;
    gsize wrote = 0;
    gchar * tr;

    if (codeset != null) {
        /* the user specified the source encoding */
        tr = g_convert(orig, -1, "UTF-8", codeset,
                       null, &wrote, &gerr);
    } else if (g_get_charset(&charset)) {
        /* we're in a UTF-8 locale, so we know that
           g_locale_to_utf8 won't do the job; so guess
           the content is iso-8859-something?
        */
        tr = g_convert(orig, -1, "UTF-8", "ISO-8859-15",
                       null, &wrote, &gerr);
    } else {
        /* try assuming the material is in the locale
           encoding */
        tr = g_locale_to_utf8(orig, -1, null, &wrote, &gerr);
        if (gerr != null) {
            /* failed: try iso-8859-15? */
            g_error_free(gerr);
            gerr = null;
            tr = g_convert(orig, -1, "UTF-8", "ISO-8859-15",
                           null, &wrote, &gerr);
        }
    }

    if (gerr != null) {
        gretl_errmsg_set(gerr.message);
        *err = GretlError.data;
        g_error_free(gerr);
    }

    g_free(orig);

    return tr;
}

private int shell_grab(const(char) * arg, char ** sout)
{
    int err = 0;

    if (arg == null || *arg == '\0') {
        return GretlError.parse;
    }

    if (!libset_get_bool(SetKey.shell_ok)) {
        gretl_errmsg_set(_("The shell command is not activated.").ptr);
        return 1;
    }

    gretl_shell_grab(arg, sout);

    if (sout != null && *sout != null) {
        char * content = *sout;

        if (!g_utf8_validate(content, -1, null)) {
            content = recode_content(content, null, &err);
            *sout = content;
        }

        if (content != null) {
            /* trim trailing newline */
            int n = cast(int) strlen(content);

            if (content[n-1] == '\n') {
                content[n-1] = '\0';
            }
        }
    }

    return err;
}

char * gretl_backtick(const(char) * arg, int * err)
{
    char * val = null;

    *err = shell_grab(arg, &val);

    if (!*err && val == null) {
        val = gretl_strdup("");
        if (val == null) {
            *err = GretlError.alloc;
        }
    }

    return val;
}

char * gretl_getenv(const(char) * key, int * defined, int * err)
{
    char * test = getenv(key);
    char * val = null;

    if (test == null) {
        *defined = 0;
        val = gretl_strdup("");
    } else {
        *defined = 1;
        val = gretl_strdup(test);
    }

    if (val == null) {
        *err = GretlError.alloc;
    }

    return val;
}

char * retrieve_date_string(int t, const(DatasetRaw) * dset, int * err)
{
    char * ret = null;

    if (t <= 0 || t > dset.n) {
        *err = GretlError.data;
    } else if (dset.S != null) {
        ret = gretl_strdup(dset.S[t-1]);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    } else {
        char[OBSLEN] datestr = 0;

        ntolabel(datestr.ptr, t - 1, dset);
        ret = gretl_strdup(datestr.ptr);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }

    return ret;
}

/* returns a gretl_array of strings on success */

GretlArray * retrieve_date_strings(const(GretlMatrixRaw) * v,
                                   const(DatasetRaw) * dset,
                                   int * err)
{
    GretlArray * ret = null;
    char * s = null;
    int i, t, n;

    n = gretl_vector_get_length(v);
    if (n == 0) {
        *err = GretlError.invarg;
    } else {
        ret = gretl_array_new(GretlType.strings, n, err);
    }

    for (i=0; i<n && !*err; i++) {
        t = gretl_int_from_double(v.val[i], err);
        if (!*err) {
            s = retrieve_date_string(t, dset, err);
        }
        if (!*err) {
            gretl_array_set_data(ret, i, s);
        }
    }

    if (*err && ret != null) {
        gretl_array_destroy(ret);
        ret = null;
    }

    return ret;
}

private int is_web_resource(const(char) * s)
{
    if (!strncmp(s, "http://", 7) ||
        !strncmp(s, "https://", 8) ||
        !strncmp(s, "ftp://", 6)) {
        return 1;
    } else {
        return 0;
    }
}

private gchar * gzipped_file_get_content(const(char) * fname,
                                         int * err)
{
    gzFile fz = gretl_gzopen(fname, "rb");
    gchar * ret = null;

    if (fz == null) {
        *err = GretlError.fopen;
    } else {
        size_t len = 0;
        int chk;

        while (gzgetc(fz) > 0) {
            len++;
        }
        if (len > 0) {
            gzrewind(fz);
            ret = cast(gchar *) malloc(len + 1);
            if (ret == null) {
                *err = GretlError.alloc;
            } else {
                chk = gzread(fz, ret, cast(uint) len);
                if (chk <= 0) {
                    *err = GretlError.data;
                }
                ret[len] = '\0';
            }
        } else {
            ret = g_strdup("");
        }
        gzclose(fz);
    }

    return ret;
}

private gchar * regular_file_get_content(const(char) * fname,
                                         int * err)
{
    GError * gerr = null;
    gchar * ret = null;
    size_t len = 0;

    g_file_get_contents(fname, &ret, &len, &gerr);

    if (gerr != null) {
        gretl_errmsg_set(gerr.message);
        *err = GretlError.fopen;
        g_error_free(gerr);
    }

    return ret;
}

char * retrieve_file_content(const(char) * fname, const(char) * codeset,
                             int * err)
{
    char * ret = null;
    gchar * content = null;
    size_t len = 0;
    gssize sz;

    if (fname == null || *fname == '\0') {
        *err = GretlError.invarg;
    } else {
        len = strlen(fname);
        if (len >= MAXLEN) {
            gretl_errmsg_sprintf(_("filename too long (%d bytes)").ptr, cast(int) len);
            *err = GretlError.invarg;
        }
        len = 0;
    }

    if (*err) {
        return null;
    }

    if (is_web_resource(fname)) {
        /* C: USE_CURL undefined here */
        gretl_errmsg_set(_("Internet access not supported").ptr);
        *err = GretlError.data;
    } else {
        char[FILENAME_MAX] fullname = 0;

        strcpy(fullname.ptr, fname);
        gretl_addpath(fullname.ptr, 0);
        if (is_gzipped(fullname.ptr)) {
            content = gzipped_file_get_content(fullname.ptr, err);
        } else {
            content = regular_file_get_content(fullname.ptr, err);
        }
    }

    sz = (len > 0)? cast(gssize) len : -1;
    if (content != null && !g_utf8_validate(content, sz, null)) {
        content = recode_content(content, codeset, err);
    }

    if (content != null) {
        if (*err == 0) {
            ret = gretl_strdup(content);
        }
        g_free(content);
    }

    return ret;
}

/* inserting string into format portion of (s)printf command:
   double any backslashes to avoid breakage of Windows paths
*/

private char * mod_strdup(const(char) * s)
{
    char * ret = null;
    int i, n = cast(int) strlen(s);
    int bs = 0;

    for (i=0; i<n; i++) {
        if (s[i] == '\\' && (i == n - 1 || s[i+1] != '\\')) {
            bs++;
        }
    }

    ret = cast(char *) malloc(n + 1 + bs);
    if (ret == null) {
        return null;
    }

    if (bs == 0) {
        strcpy(ret, s);
    } else {
        int j = 0;

        for (i=0; i<n; i++) {
            if (s[i] == '\\' && (i == n - 1 || s[i+1] != '\\')) {
                ret[j++] = '\\';
            }
            ret[j++] = s[i];
        }
        ret[j] = '\0';
    }

    return ret;
}

private char * maybe_get_subst(char * name, int * n, int quoted,
                               int * freeit)
{
    char * s = null;
    int k = *n - 1;
    char * ret = null;

    while (k >= 0) {
        s = cast(char *) get_string_by_name(name);
        if (s != null) {
            *n = k + 1;
            ret = s;
            break;
        }
        name[k--] = '\0';
    }

    if (ret != null) {
        if (quoted && strchr(ret, '\\')) {
            ret = mod_strdup(ret);
            *freeit = 1;
        }
    }

    return ret;
}

private void too_long()
{
    gretl_errmsg_sprintf((_("Maximum length of command line ")
                          ~ _("(%d bytes) exceeded\n")).ptr, MAXLINE);
}

int substitute_named_strings(char * line, int * subst)
{
    char[VNAMELEN] sname = 0;
    int len = cast(int) strlen(line);
    char * sub;
    char * tmp;
    char * s = line;
    int bs = 0, in_format = 0;
    int freeit;
    int n, m, err = 0;

    *subst = 0;

    if (*s == '#' || strchr(s, '@') == null) {
        return 0;
    }

    if (!strncmp(line, "printf", 6) || !strncmp(line, "sprintf", 7)) {
        s = strchr(s, '"');
        if (s == null) {
            /* no format string */
            return GretlError.parse;
        }
        s++;
        in_format = 1;
    }

    while (*s && !err) {
        if (in_format) {
            if (*s == '"' && (bs % 2 == 0)) {
                /* reached end of (s)printf format string */
                in_format = 0;
            }
            if (*s == '\\') {
                bs++;
            } else {
                bs = 0;
            }
        }
        if (*s == '@') {
            n = gretl_namechar_spn(s + 1);
            if (n > 0) {
                if (n >= VNAMELEN) {
                    n = VNAMELEN - 1;
                }
                sname[0] = '\0';
                strncat(sname.ptr, s + 1, n);
                freeit = 0;
                sub = maybe_get_subst(sname.ptr, &n, in_format, &freeit);
                if (sub != null) {
                    m = cast(int) strlen(sub);
                    if (len + m + 2 >= MAXLINE) {
                        too_long();
                        err = 1;
                        break;
                    }
                    tmp = gretl_strdup(s + n + 1);
                    if (tmp == null) {
                        err = GretlError.alloc;
                    } else {
                        strcpy(s, sub);
                        strcpy(s + m, tmp);
                        free(tmp);
                        len += m - (n + 1);
                        s += m - 1;
                        *subst = 1;
                    }
                    if (freeit) {
                        free(sub);
                    }
                }
            }
        }
        s++;
    }

    return err;
}
