/* Port of gretl_list.c */

module gretl.gretl_list;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.errno : errno;
import core.stdc.math : ceil, floor, log, log10, fabs;
import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : abs, calloc, free, malloc, qsort, realloc, strtol;
import core.stdc.string : memcpy, strcat, strcmp, strcpy, strcspn, strlen, strncat, strspn;

import gretl.dataset : dataset_add_NA_series, dataset_add_series, dataset_drop_last_variables,
    series_get_mtime, series_get_midas_period, series_get_stack_level,
    series_is_discrete, series_is_hidden, series_is_midas_anchor, series_record_label,
    series_set_midas_anchor, series_set_midas_freq, series_set_midas_period;
import gretl.dataseth : sample_size;
import gretl.gmem : g_free;
import gretl.gpatternspec;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_matrix_free, gretl_matrix_get_colnames, gretl_matrix_get_t1,
    gretl_matrix_get_t2, gretl_matrix_values, gretl_null_matrix_new;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gstrfuncs : g_strchomp, g_strchug, g_strdup_printf;
import gretl.libgretlh;
import gretl.missingh : na, NADBL;
import gretl.strutils : count_fields, gretl_charsub, gretl_strdup, gretl_string_split, strings_array_free, strings_array_new, tailstrip;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_utils : gretl_compare_ints, gretl_int_from_double, positive_int_from_string, true_const;
import gretl.dataio : check_varname;
import gretl.genmain : current_series_index, gretl_is_series, series_index;
import gretl.uservar : get_list_by_name, gretl_is_user_var;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    int gretl_function_depth();
    int get_midas_frequency(const(DatasetRaw) * dset, int m);
    int is_valid_midas_frequency_ratio(const(DatasetRaw) * dset, int m);
}

/**
 * SECTION:gretl_list
 * @short_description: handling of lists of integers
 * @title: Lists
 * @include: libgretl.h
 *
 * Lists of integers are used in many contexts in libgretl, e.g.
 * for holding the ID numbers of variables in a regression
 * specification. A gretl "list" is simply an array of ints
 * following a definite convention: the value at position 0
 * gives the number of elements that follow. The total number
 * of ints in the list foo is therefore foo[0] + 1, and reading
 * the substantive members of foo involves looping from
 * position 1 to position foo[0].
 */

/**
 * gretl_is_midas_list:
 * @list: the list array.
 * @dset: pointer to dataset.
 *
 * Returns: 1 if @list has been set as a MIDAS list in an
 * approved manner, 0 otherwise.
 */

int gretl_is_midas_list(const(int) * list, const(DatasetRaw) * dset)
{
    int ret = 0;

    if (list != null && list[0] > 2) {
        int i, m, p;

        if (!series_is_midas_anchor(dset, list[1])) {
            return 0;
        }

        m = series_get_midas_period(dset, list[1]);
        if (!is_valid_midas_frequency_ratio(dset, m) || list[0] != m) {
            return 0;
        }

        ret = 1;
        for (i=2; i<=list[0] && ret; i++) {
            p = series_get_midas_period(dset, list[i]);
            if (p != m - 1) {
                ret = 0;
            }
            m = p;
        }
    }

    return ret;
}

/**
 * gretl_list_set_midas:
 * @list: the list array.
 * @dset: pointer to dataset.
 *
 * Attempts to set the MIDAS flag on the members of @list.
 *
 * 0 on success, non-zero code on failure.
 */

int gretl_list_set_midas(const(int) * list, DatasetRaw * dset)
{
    int err = 0;

    if (list != null) {
        int i, m = list[0];

        if (!is_valid_midas_frequency_ratio(dset, m)) {
            err = GretlError.invarg;
        }

        for (i=1; i<=list[0] && !err; i++) {
            if (list[i] < 1 || list[i] >= dset.v) {
                err = GretlError.invarg;
            }
        }

        if (!err) {
            int freq = get_midas_frequency(dset, m);

            series_set_midas_anchor(dset, list[1]);
            series_set_midas_period(dset, list[1], m);
            series_set_midas_freq(dset, list[1], freq);

            for (i=2; i<=list[0]; i++) {
                series_set_midas_period(dset, list[i], m - i + 1);
                series_set_midas_freq(dset, list[i], freq);
            }
        }
    }

    return err;
}

/**
 * gretl_list_new:
 * @nterms: the maximum number of elements to be stored in the list.
 *
 * Creates a newly allocated list with space for @nterms elements,
 * besides the leading element, which in a gretl list always
 * holds a count of the number of elements that follow.  This
 * leading element is initialized appropriately.  For example, if
 * @nterms = 4, space for 5 integers is allocated and the first
 * element of the array is set to 4.  The other elements of
 * the list are initialized to 0.
 *
 * Returns: the newly allocated list, or NULL on failure.
 */

int * gretl_list_new(int nterms)
{
    int * list = null;
    int i;

    if (nterms < 0) {
        return null;
    }

    /* Guard against creation of a list with more elements than
       can be indexed.
    */
    if (cast(size_t) nterms > int.max - 1) {
        return null;
    }

    list = cast(int *) malloc((cast(size_t) nterms + 1) * (*list).sizeof);

    if (list != null) {
        list[0] = nterms;
        for (i=1; i<=nterms; i++) {
            list[i] = 0;
        }
    }

    return list;
}

/* experimental wrapper for use with gretl4py */

void gretl_list_free(int * list)
{
    free(list);
}

/**
 * gretl_list_array_new:
 * @nlists: the number of lists to create.
 * @nterms: the maximum number of elements to be stored in each list.
 *
 * Creates an array of newly allocated lists, each of which as described
 * in connection with gretl_list_new().
 *
 * Returns: the newly allocated lists, or NULL on failure.
 */

int ** gretl_list_array_new(int nlists, int nterms)
{
    int ** lists = null;
    int i;

    if (nlists < 0) {
        return null;
    }

    lists = cast(int **) malloc(nlists * (*lists).sizeof);

    if (lists != null) {
        for (i=0; i<nlists; i++) {
            lists[i] = gretl_list_new(nterms);
        }
    }

    return lists;
}

/**
 * gretl_list_array_free:
 * @lists: array of gretl lists.
 * @nlists: the number of lists in @lists.
 *
 * Frees all the lists in @lists and also the top-level pointer.
 */

void gretl_list_array_free(int ** lists, int nlists)
{
    if (lists != null) {
        int i;

        for (i=0; i<nlists; i++) {
            free(lists[i]);
        }
        free(lists);
    }
}

/**
 * gretl_consecutive_list_new:
 * @lmin: starting value for consecutive list elements.
 * @lmax: ending value.
 *
 * Creates a newly allocated list whose elements run from
 * @lmin to @lmax consecutively.
 *
 * Returns: the newly allocated list, or NULL on failure.
 */

int * gretl_consecutive_list_new(int lmin, int lmax)
{
    int * list = null;
    int i, n;

    n = lmax - lmin + 1;
    if (n <= 0) {
        return null;
    }

    list = gretl_list_new(n);

    if (list != null) {
        for (i=0; i<n; i++) {
            list[i+1] = lmin + i;
        }
    }

    return list;
}

/**
 * gretl_list_resize:
 * @oldlist: pointer to list to be resized.
 * @nterms: the new maximum number of elements for the list.
 *
 * Resizes the content of @oldlist to hold @nterms, and adjusts
 * the first element to reflect the new size.  If the new
 * list is longer than the old, the extra elements are initialized
 * to zero.
 *
 * Returns: the resized list, or NULL on failure.
 */

int * gretl_list_resize(int ** oldlist, int nterms)
{
    int * list = null;
    int i, oldn = 0;

    if (nterms < 0 || oldlist == null) {
        return null;
    }

    if (*oldlist != null) {
        oldn = (*oldlist)[0];
        if (nterms == oldn) {
            /* no-op */
            return *oldlist;
        }
    }

    list = cast(int *) realloc(*oldlist, (nterms + 1) * (*list).sizeof);

    if (list != null) {
        list[0] = nterms;
        *oldlist = list;
        for (i=oldn+1; i<=list[0]; i++) {
            list[i] = 0;
        }
    } else {
        free(*oldlist);
        *oldlist = null;
    }

    return list;
}

/**
 * gretl_list_append_term:
 * @plist: pointer to list to be augmented.
 * @v: the term to be added.
 *
 * Resizes (or allocates from scratch) the content of @plist,
 * so that it can hold one extra element, and sets the last
 * element to @v.
 *
 * Returns: the augmented list, or NULL on failure.
 */

int * gretl_list_append_term(int ** plist, int v)
{
    int * list = null;

    if (*plist == null) {
        list = gretl_list_new(1);
        if (list != null) {
            list[1] = v;
        }
    } else {
        int oldn = (*plist)[0];

        list = cast(int *) realloc(*plist, (oldn + 2) * (*list).sizeof);
        if (list != null) {
            list[0] += 1;
            list[list[0]] = v;
        } else {
            free(*plist);
        }
    }

    *plist = list;

    return list;
}

/**
 * gretl_list_sort:
 * @list: list to be sorted.
 *
 * Sorts the elements from position 1 to the end of @list
 * in ascending order.
 *
 * Returns: the sorted list.
 */

int * gretl_list_sort(int * list)
{
    int i, sorted = 1;

    for (i=1; i<list[0]; i++) {
        if (list[i] > list[i+1]) {
            sorted = 0;
            break;
        }
    }

    if (!sorted) {
        qsort(list + 1, list[0], (*list).sizeof, &gretl_compare_ints);
    }

    return list;
}

/**
 * gretl_list_cmp:
 * @list1: gretl list.
 * @list2: gretl list.
 *
 * Returns: 0 if @list1 and @list2 have identical content,
 * otherwise 1.
 */

int gretl_list_cmp(const(int) * list1, const(int) * list2)
{
    int i;

    if (list1 == null && list2 != null) {
        return 1;
    } else if (list1 != null && list2 == null) {
        return 1;
    } else if (list1 == null && list2 == null) {
        return 0;
    }

    for (i=0; i<=list1[0]; i++) {
        if (list2[i] != list1[i]) {
            return 1;
        }
    }

    return 0;
}

/**
 * gretl_null_list:
 *
 * Creates a newly allocated "list" with only one member,
 * which is set to zero.
 *
 * Returns: the newly allocated list, or NULL on failure.
 */

int * gretl_null_list()
{
    int * list = cast(int *) malloc(int.sizeof);

    if (list != null) {
        list[0] = 0;
    }

    return list;
}

/**
 * gretl_is_null_list:
 * @list: list to check.
 *
 * Returns: 1 if @list is NULL or empty.
 */

int gretl_is_null_list(const(int) * list)
{
    return (list == null || list[0] == 0);
}

/**
 * gretl_list_copy:
 * @src: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Returns: an allocated copy @src (or NULL if @src is NULL).
 */

int * gretl_list_copy(const(int) * src)
{
    int * targ = null;

    if (src != null) {
        int n = src[0] + 1;

        targ = cast(int *) malloc(n * (*targ).sizeof);
        if (targ != null) {
            memcpy(targ, src, n * (*targ).sizeof);
        }
    }

    return targ;
}

/**
 * gretl_list_copy_from_pos:
 * @src: an array of integers, the first element of which holds
 * a count of the number of elements following.
 *
 * Returns: an allocated copy @src from position @pos onward
 * (or NULL on failure).
 */

int * gretl_list_copy_from_pos(const(int) * src, int pos)
{
    int * targ = null;
    int i, n;

    if (src != null && ((n = src[0] - pos + 1) > 0)) {
        targ = gretl_list_new(n);
        if (targ != null) {
            for (i=1; i<=n; i++) {
                targ[i] = src[i+pos-1];
            }
        }
    }

    return targ;
}

/**
 * gretl_list_from_string:
 * @str: string representation of list of integers.
 * @err: location to receive error code.
 *
 * Reads a string containing a list of integers, separated by
 * spaces and/or commas and possibly wrapped in parentheses,
 * and constructs an array of these integers.  The first
 * element is the number of integers that follow.
 * This function supports an abbreviation for consecutive
 * (increasing) integers in the list, using the notation, e.g.,
 * "1-4" as shorthand for "1 2 3 4".
 *
 * Returns: the allocated array, or NULL on failure.
 */

int * gretl_list_from_string(const(char) * str, int * err)
{
    char * p, q, s, next;
    int i, r1, r2, rg;
    int * list;
    int n = 0;

    if (str == null || *str == '\0') {
        *err = GretlError.invarg;
        return null;
    }

    /* 'p' marks the memory to be freed */
    p = s = gretl_strdup(str);
    if (s == null) {
        *err = GretlError.alloc;
        return null;
    }

    *err = 0;

    /* strip white space at both ends */
    while (isspace(*s)) s++;
    tailstrip(s);

    /* strip parentheses or braces, if present */
    if (*s == '(' || *s == '{') {
        char close = (*s == '(')? ')' : '}';

        n = cast(int) strlen(s);
        if (s[n-1] != close) {
            /* got opening grouping character but no close */
            *err = GretlError.parse;
            return null;
        }
        s[n-1] = '\0';
        s++;
        while (isspace(*s)) s++;
        tailstrip(s);
    }

    q = s; /* copy relevant starting point */

    gretl_charsub(s, ',', ' ');

    errno = 0;

    /* first pass: figure out the number of values
       in the list, checking for errors as we go
    */

    n = 0; /* value counter */

    while (*s && !*err) {
        s += strspn(s, " ");
        if (n > 0 && *s == ';') {
            /* list separator */
            n++;
            s++;
            continue;
        }
        r1 = cast(int) strtol(s, &next, 10);
        if (errno || next == s) {
            fprintf(stderr, "gretl_list_from_string: '%s'\n", s);
            *err = GretlError.parse;
        } else {
            s = next;
            if (*s == '-') {
                /* hyphen indicating range? */
                s++;
                r2 = cast(int) strtol(s, &next, 10);
                if (errno || next == s) {
                    *err = GretlError.parse;
                } else if (r2 < r1) {
                    *err = GretlError.parse;
                } else {
                    n += r2 - r1 + 1;
                }
                s = next;
            } else {
                /* single numerical value */
                n++;
            }
        }
    }

    if (*err || n == 0) {
        free(p);
        return null;
    }

    list = gretl_list_new(n);
    if (list == null) {
        *err = GretlError.alloc;
        free(p);
        return null;
    }

    /* second pass: fill out the list (no error
       checking should be needed at this stage)
    */

    s = q; /* back to start of string */
    n = 1; /* list position indicator */

    while (*s) {
        s += strspn(s, " ");
        if (*s == ';') {
            list[n++] = LISTSEP;
            s++;
            continue;
        }
        r1 = cast(int) strtol(s, &s, 10);
        if (*s == '-') {
            s++;
            r2 = cast(int) strtol(s, &s, 10);
            rg = r2 - r1 + 1;
            for (i=0; i<rg; i++) {
                list[n++] = r1 + i;
            }
        } else {
            list[n++] = r1;
        }
    }

    free(p);

    return list;
}

/**
 * gretl_list_from_varnames:
 * @str: string holding space-separated series names.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Returns: an allocated gretl list holding the ID numbers of
 * the series named in @str, or NULL on failure.
 */

int * gretl_list_from_varnames(const(char) * str,
                               const(DatasetRaw) * dset,
                               int * err)
{
    int * list = null;
    char ** S;
    int n = 0;

    if (str == null || count_fields(str, null) < 1) {
        *err = GretlError.data;
        return null;
    }

    S = gretl_string_split(str, &n, null);
    if (S == null) {
        *err = GretlError.alloc;
        return null;
    }

    list = gretl_list_new(n);

    if (list == null) {
        *err = GretlError.alloc;
    } else {
        int i, vi;

        for (i=0; i<n; i++) {
            if (!strcmp(S[i], "0")) {
                vi = 0;
            } else {
                vi = current_series_index(dset, S[i]);
            }
            if (vi < 0) {
                *err = GretlError.unkvar;
                break;
            } else {
                list[i+1] = vi;
            }
        }
    }

    strings_array_free(S, n);

    return list;
}

private enum LN_10 = 2.30258509299404590;

private int integer_length(int k)
{
    double x = k;
    int len = 0;

    if (x < 0) {
        x = fabs(x);
        len = 1;
    }

    if (x < 10) {
        len += 1;
    } else {
        len += cast(int) ceil(log(x)/LN_10);
        len += (k % 10 == 0);
    }

    return len;
}

/**
 * gretl_list_to_numeric_string:
 * @list: array of integers.
 *
 * Prints the given @list of integers into a newly
 * allocated string, separated by single spaces and with
 * one leading space.
 *
 * Returns: The string representation of the list on success,
 * or NULL on failure.
 */

char * gretl_list_to_numeric_string(const(int) * list)
{
    char * buf;
    int i, len = 1;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            len += 2;
        } else {
            len += integer_length(list[i]) + 1;
        }
    }

    if (len > MAXLINE - 32) {
        /* string would be too long for command line */
        return null;
    }

    buf = cast(char *) malloc(len);

    if (buf != null) {
        char[16] numstr;

        *buf = '\0';
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                strcat(buf, " ;");
            } else {
                sprintf(numstr.ptr, " %d", list[i]);
                strcat(buf, numstr.ptr);
            }
        }
    }

    return buf;
}

/**
 * gretl_list_to_string:
 * @list: array of integers.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Returns: allocated string representation of @list, with ID
 * numbers cashed out as series names (and with one leading
 * space), or NULL on failure. The list separator #LISTSEP,
 * is accepted in the incoming @list, otherwise all terms
 * must be integers in the range 0 to the greatest current
 * series ID within @dset.
 */

char * gretl_list_to_string(const(int) * list,
                            const(DatasetRaw) * dset,
                            int * err)
{
    char * buf = null;
    int len = 1;
    int i, vi;

    if (list == null) {
        *err = GretlError.data;
        return null;
    }

    if (list[0] == 0) {
        return gretl_strdup("");
    }

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi == LISTSEP) {
            len += 2;
        } else if (vi >= 0 && vi < dset.v) {
            len += strlen(dset.varname[vi]) + 1;
        } else {
            *err = GretlError.data;
            return null;
        }
    }

    buf = cast(char *) calloc(len, 1);

    if (buf == null) {
        *err = GretlError.alloc;
    } else {
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            if (vi == LISTSEP) {
                strcat(buf, " ;");
            } else {
                strcat(buf, " ");
                strcat(buf, dset.varname[vi]);
            }
        }
    }

    return buf;
}

private int get_excluded_id(const(int) * list, int n)
{
    int i;

    for (i=1; i<n; i++) {
        if (!in_gretl_list(list, i)) {
            return i;
        }
    }

    return 0;
}

/**
 * gretl_list_to_compact_string:
 * @list: array of integers.
 * @dset: pointer to dataset.
 * @argstyle: boolean for function-argument variant.
 * @err: location to receive error code.
 *
 * Designed to produce a reasonably compact string representation
 * of the @list argument. If the list is short enough the task is
 * handed off to gretl_list_to_string(). Otherwise we try to use
 * either the "dataset" macro, or the ".." notation for a list of
 * contiguous series IDs. The fallback for a large list, if these
 * tactics do not seem applicable, is to compose a string holding
 * the series ID numbers.
 *
 * If @argstyle is non-zero the returned string will suitable for
 * giving as an argument to the userspace deflist function.
 *
 * Returns: allocated string representation of @list, or NULL
 * on error.
 */

char * gretl_list_to_compact_string(const(int) * list,
                                    const(DatasetRaw) * dset,
                                    int argstyle,
                                    int * err)
{
    char * buf = null;
    int * mylist = null;
    char[8] cstr;
    char sep;
    int buflen;
    int nv = 0;      /* number of 'real' series */
    int v0 = -1;     /* first non-zero series ID */
    int vn = -1;     /* last non-zero series ID */
    int contig = 1;  /* contiguity flag */
    int vprev = -1;  /* ID of previous series */
    int idmax = 0;   /* greatest series ID */
    int cpos = 0;    /* position of constant */
    int i, vi;

    if (argstyle) {
        strcpy(cstr.ptr, "const,");
        sep = ',';
    } else {
        strcpy(cstr.ptr, "const ");
        sep = ' ';
    }

    if (list == null) {
        *err = GretlError.data;
        return null;
    } else if (list[0] == 0) {
        /* empty list */
        return gretl_strdup("");
    } else if (list[0] <= 6) {
        /* short list: don't bother with the fancy stuff */
        buf = gretl_list_to_string(list, dset, err);
        if (*err) {
            return buf;
        } else if (argstyle) {
            return gretl_charsub(g_strchug(buf), ' ', ',');
        } else {
            return g_strchug(buf);
        }
    }

    /* analyse the incoming list */
    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi == LISTSEP) {
            gretl_errmsg_set("Invalid list element");
            *err = GretlError.invarg;
            return null;
        } else if (vi == 0) {
            /* record the const position, if present */
            cpos = i;
        } else {
            if (v0 < 0) {
                /* record the first "real" series */
                v0 = vi;
            }
            if (vi > idmax) {
                /* highest-numbered series in list */
                idmax = vi;
            }
            if (contig && vprev > 0 && vi != vprev + 1) {
                /* list is not contiguous */
                contig = 0;
            }
            vprev = vi;
        }
    }

    if (v0 < 0) {
        /* degenerate case: no "real" (non-zero) series ID was found
           in the list, e.g. it consists entirely of duplicated
           zeros (the constant). The "almost complete", ".." and
           numeric-ID branches below all assume that at least one
           real series ID is present, so fall back to the plain,
           uncompacted rendering instead.
        */
        buf = gretl_list_to_string(list, dset, err);
        if (*err) {
            return buf;
        } else if (argstyle) {
            return gretl_charsub(g_strchug(buf), ' ', ',');
        } else {
            return g_strchug(buf);
        }
    }

    if (cpos > 0) {
        /* delete the constant */
        mylist = gretl_list_copy(list);
        gretl_list_delete_at_pos(mylist, cpos);
    } else {
        mylist = cast(int *) list;
    }

    nv = mylist[0];
    vn = mylist[nv];

    if (nv >= dset.v - 2) {
        /* an almost-complete list: use the "dataset" macro,
           maybe with const prepended and/or one exclusion
        */
        if (nv == dset.v - 1) {
            /* no "real" series excluded */
            buflen = 6 * (cpos > 0) + 9;
            buf = cast(char *) calloc(buflen, 1);
            sprintf(buf, "%sdataset", cpos > 0 ? cstr.ptr : "");
        } else {
            /* one "real" series excluded */
            int vx = get_excluded_id(mylist, dset.v);
            const(char) * sx = dset.varname[vx];

            buflen = cast(int) (12 + 6 * (cpos > 0) + strlen(sx));
            buf = cast(char *) calloc(buflen, 1);
            sprintf(buf, "%sdataset-%s", cpos > 0 ? cstr.ptr : "", sx);
        }
    } else if (contig) {
        /* use ".." notation for contiguous series IDs from
           @v0 to @vn, possibly with const prepended
        */
        const(char) * s0 = dset.varname[v0];
        const(char) * sn = dset.varname[vn];

        buflen = cast(int) (4 + 6 * (cpos > 0) + strlen(s0) + strlen(sn));
        buf = cast(char *) calloc(buflen, 1);
        sprintf(buf, "%s%s..%s", cpos > 0 ? cstr.ptr : "", s0, sn);
    } else {
        /* use straight numeric IDs */
        int numlen = cast(int) floor(log10(cast(double) idmax));
        char[16] tmp;

        buflen = 2 + 2 * (cpos > 0) + nv * (2 + numlen);
        buf = cast(char *) calloc(buflen, 1);
        sprintf(buf, "%s%c", cpos > 0 ? "0".ptr : "".ptr, cpos > 0 ? cast(int) sep : 0);
        for (i=1; i<=mylist[0]; i++) {
            sprintf(tmp.ptr, "%d%c", mylist[i], cast(int) sep);
            strcat(buf, tmp.ptr);
        }
        g_strchomp(buf);
    }

    if (mylist != list) {
        free(mylist);
    }

    return buf;
}

/**
 * gretl_list_to_vector:
 * @list: array of integers.
 * @err: location to receive error code.
 *
 * Returns: allocated representation of @list as a row vector
 * or NULL on failure.
 */

GretlMatrixRaw * gretl_list_to_vector(const(int) * list, int * err)
{
    GretlVector * v = null;

    if (list == null) {
        *err = GretlError.data;
    } else {
        int i, n = list[0];

        if (n == 0) {
            v = gretl_null_matrix_new();
            if (v == null) {
                *err = GretlError.alloc;
            }
        } else if (n > 0) {
            v = gretl_vector_alloc(n);
            if (v == null) {
                *err = GretlError.alloc;
            } else {
                for (i=0; i<n; i++) {
                    v.val[i] = list[i+1];
                }
            }
        } else {
            *err = GretlError.data;
        }
    }

    return v;
}

/**
 * gretl_list_get_names:
 * @list: array of integers.
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Prints the names of the members of @list of integers into
 * a newly allocated string, separated by commas.
 *
 * Returns: allocated string on success or NULL on failure.
 */

char * gretl_list_get_names(const(int) * list, const(DatasetRaw) * dset,
                            int * err)
{
    char * buf = null;
    int len = 0;
    int i, vi;

    if (list == null) {
        *err = GretlError.data;
        return null;
    }

    if (list[0] == 0) {
        return gretl_strdup("");
    }

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi < 0 || vi >= dset.v) {
            len += cast(int) strlen("unknown") + 1;
        } else {
            len += cast(int) strlen(dset.varname[vi]) + 1;
        }
    }

    buf = cast(char *) malloc(len);
    if (buf == null) {
        *err = GretlError.alloc;
        return null;
    }

    *buf = '\0';

    for (i=1; i<=list[0]; i++) {
        vi = list[i];
        if (vi < 0 || vi >= dset.v) {
            strcat(buf, "unknown");
        } else {
            strcat(buf, dset.varname[vi]);
        }
        if (i < list[0]) {
            strcat(buf, ",");
        }
    }

    return buf;
}

/**
 * gretl_list_get_names_array:
 * @list: array of integers.
 * @dset: dataset information.
 * @err: location to receive error code, or #NULL.
 *
 * Returns: An array of strings holding the names of the
 * members of @list, or NULL on failure.
 */

char ** gretl_list_get_names_array(const(int) * list,
                                   const(DatasetRaw) * dset,
                                   int * err)
{
    char ** S = null;
    int i, vi, n;

    if (list == null || dset == null) {
        if (err != null) {
            *err = GretlError.data;
        }
        return null;
    }

    if (list[0] == 0) {
        return null;
    }

    n = list[0];

    S = strings_array_new(n);
    if (S == null) {
        if (err != null) {
            *err = GretlError.alloc;
        }
        return null;
    }

    for (i=0; i<n; i++) {
        vi = list[i+1];
        if (vi < 0 || vi >= dset.v) {
            S[i] = gretl_strdup("unknown");
        } else {
            S[i] = gretl_strdup(dset.varname[vi]);
        }
        if (S[i] == null) {
            if (err != null) {
                *err = GretlError.alloc;
            }
            strings_array_free(S, n);
            S = null;
            break;
        }
    }

    return S;
}

/**
 * gretl_list_to_lags_string:
 * @list: array of integers.
 * @err: location to receive error code.
 *
 * Prints the given @list of integers into a newly
 * allocated string, enclosed by braces and separated by commas.
 * Will fail if @list contains any numbers greater than 998.
 *
 * Returns: The string representation of the list on success,
 * or NULL on failure.
 */

char * gretl_list_to_lags_string(const(int) * list, int * err)
{
    char * buf;
    char[16] numstr;
    int len, i;

    len = 4 * (list[0] + 1) + 2;

    if (len > MAXLINE - 32) {
        *err = GretlError.data;
        return null;
    }

    buf = cast(char *) calloc(len, 1);
    if (buf == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=1; i<=list[0]; i++) {
        if (abs(list[i] >= 999)) {
            *err = GretlError.data;
            break;
        } else {
            if (i == 1) {
                sprintf(numstr.ptr, "{%d", list[i]);
            } else {
                sprintf(numstr.ptr, ",%d", list[i]);
            }
            strcat(buf, numstr.ptr);
        }
    }
    strcat(buf, "}");

    if (*err) {
        free(buf);
        buf = null;
    }

    return buf;
}

/**
 * in_gretl_list:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @k: integer to test.
 *
 * Checks whether @k is present among the members of @list,
 * in position 1 or higher.
 *
 * Returns: the position of @k in @list, or 0 if @k is not
 * present.
 */

int in_gretl_list(const(int) * list, int k)
{
    int i;

    if (list != null) {
        for (i=1; i<=list[0]; i++) {
            if (list[i] == k) {
                return i;
            }
        }
    }

    return 0;
}

private void reglist_move_const(int * list, int k)
{
    int i, cnum = list[k];

    for (i=k; i>2; i--) {
        list[i] = list[i-1];
    }

    list[2] = cnum;
}

/**
 * reglist_check_for_const:
 * @list: regression list suitable for use with a gretl
 * model (should not contain #LISTSEP).
 * @dset: dataset struct.
 *
 * Checks @list for an intercept term (a variable all of
 * whose valid values in sample are 1).  If such a variable
 * is present, it is moved to position 2 in the list.
 *
 * Returns: 1 if the list contains an intercept, else 0.
 */

int reglist_check_for_const(int * list, const(DatasetRaw) * dset)
{
    int cpos = gretl_list_const_pos(list, 2, dset);
    int ret = 0;

    if (cpos > 1) {
        ret = 1;
    }

    if (cpos > 2) {
        reglist_move_const(list, cpos);
    }

    return ret;
}

/**
 * gretl_list_delete_at_pos:
 * @list: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @pos: position at which to delete list element.
 *
 * Deletes the element at position @pos from @list and moves any
 * remaining elements forward.  Decrements the value of the first,
 * counter, element of @list.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_list_delete_at_pos(int * list, int pos)
{
    int i, err = 0;

    if (pos < 1 || pos > list[0]) {
        err = GretlError.invarg;
    } else {
        for (i=pos; i<list[0]; i++) {
            list[i] = list[i + 1];
        }
        list[list[0]] = 0;
        list[0] -= 1;
    }

    return err;
}

/**
 * gretl_list_purge_const:
 * @list: list of variable ID numbers.
 * @dset: dataset struct.
 *
 * Checks @list from position 1 onward for the presence of a
 * variable whose valid values in sample all equal 1.0.  If
 * such a variable is found, it is deleted from @list (that is,
 * any following elements are moved forward by one and list[0]
 * is decremented by 1).
 *
 * Returns: 1 if a constant was found and deleted, else 0.
 */

int gretl_list_purge_const(int * list, const(DatasetRaw) * dset)
{
    int i, gotc = 0;
    int l0 = list[0];

    /* Handle the case where the constant comes last; if it's
       the only element behind the list separator, remove both
       the constant and the separator.
    */
    if (list[l0] == 0 || true_const(list[l0], dset)) {
        gotc = 1;
        list[0] -= 1;
        if (list[l0 - 1] == LISTSEP) {
            list[l0 - 1] = 0;
            list[0] -= 1;
        }
    } else {
        for (i=1; i<l0; i++) {
            if (list[i] == 0 || true_const(list[i], dset)) {
                for ( ; i<l0; i++) {
                    list[i] = list[i+1];
                }
                list[l0] = 0;
                list[0] -= 1;
                gotc = 1;
                break;
            }
        }
    }

    return gotc;
}

/**
 * gretl_list_add:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @add: list of variables to be added.
 * @err: location to receive error code.
 *
 * Creates a new list containing the elements of @orig plus
 * those of @add. It is an error -- flagged by %E_ADDDUP --
 * if one or more elements of @add were already present in
 * @orig.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_add(const(int) * orig, const(int) * add, int * err)
{
    int n_orig = orig[0];
    int n_add = add[0];
    int i, j;
    int * ret;

    if (n_add == 0) {
        *err = GretlError.noadd;
        return null;
    }

    *err = 0;

    ret = gretl_list_new(n_orig + n_add);
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    for (i=0; i<=n_orig; i++) {
        ret[i] = orig[i];
    }

    j = n_orig + 1;
    for (i=1; i<=n_add && !*err; i++) {
        if (in_gretl_list(ret, add[i])) {
            /* a supposedly "new" var was already present */
            *err = GretlError.adddup;
        } else {
            ret[0] += 1;
            ret[j++] = add[i];
        }
    }

    if (*err) {
        free(ret);
        ret = null;
    }

    return ret;
}

/**
 * gretl_list_plus:
 * @l1: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @l2: list of variables to be added.
 * @err: location to receive error code.
 *
 * Creates a list containing all elements of @l1 followed
 * by all elements of @l2. This differs from gretl_list_union()
 * in that some elements may end up being repeated in the
 * returned list.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_plus(const(int) * l1, const(int) * l2, int * err)
{
    int n1 = l1[0];
    int n2 = l2[0];
    int i, j;
    int * ret;

    ret = gretl_list_new(n1 + n2);
    if (ret == null) {
        *err = GretlError.alloc;
        return null;
    }

    j = 1;

    for (i=1; i<=n1; i++) {
        ret[j++] = l1[i];
    }

    for (i=1; i<=n2; i++) {
        ret[j++] = l2[i];
    }

    return ret;
}

private int list_mark_duplicates(const(int) * l1, int * l2)
{
    int n_dups = 0;
    int i, j;

    for (i=1; i<=l2[0]; i++) {
        if (l2[i] == -1) {
            continue;
        }
        if (in_gretl_list(l1, l2[i])) {
            /* element already present in l1 */
            n_dups++;
            l2[i] = -1;
        } else {
            /* not present in @l1, but check for duplicates of
               this element in @l2
            */
            for (j=1; j<=l2[0]; j++) {
                if (j != i && l2[j] == l2[i]) {
                    n_dups++;
                    l2[j] = -1;
                }
            }
        }
    }

    return n_dups;
}

/**
 * gretl_list_union:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the union of @l1 and @l2.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_union(const(int) * l1, const(int) * l2, int * err)
{
    int * ret, lcopy;
    int n_orig = l1[0];
    int n_add = l2[0];

    *err = 0;

    lcopy = gretl_list_copy(l2);

    if (lcopy == null) {
        *err = GretlError.alloc;
        return null;
    } else {
        n_add -= list_mark_duplicates(l1, lcopy);
    }

    if (n_add == 0) {
        ret = gretl_list_copy(l1);
    } else {
        ret = gretl_list_new(n_orig + n_add);
    }

    if (ret == null) {
        *err = GretlError.alloc;
    } else if (n_add > 0) {
        int i, j;

        for (i=1; i<=n_orig; i++) {
            ret[i] = l1[i];
        }
        j = n_orig + 1;
        for (i=1; i<=lcopy[0]; i++) {
            if (lcopy[i] != -1) {
                ret[j++] = lcopy[i];
            }
        }
    }

    free(lcopy);

    return ret;
}

/* return a copy of @list from which any duplicated elements
   have been deleted
*/

private int * gretl_list_uniq(const(int) * list)
{
    int * ret = gretl_list_copy(list);
    int i, j;

    for (j=list[0]; j>0; j--) {
        for (i=1; i<j; i++) {
            if (list[j] == list[i]) {
                gretl_list_delete_at_pos(ret, j);
            }
        }
    }

    return ret;
}

/**
 * gretl_list_merge_list:
 * @targ: pointer to original list.
 * @src: list to merge.
 *
 * Appends to *targ any unique elements of @src that are not
 * already present in *targ. @targ must not be NULL, but it's
 * OK if the list to which it points is NULL.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_merge_list(int ** targ, const(int) * src)
{
    int * tmp = null;
    int * scpy = null;
    int n1 = 0, n2 = 0;
    int i, j;
    int err = 0;

    if (targ == null) {
        return GretlError.invarg;
    } else if (src == null || src[0] == 0) {
        /* no-op */
        return 0;
    }

    if (src[0] > 1) {
        src = scpy = gretl_list_uniq(src);
    }

    if (*targ == null) {
        /* special case */
        *targ = gretl_list_copy(src);
        if (*targ == null) {
            err = GretlError.alloc;
        }
        free(scpy);
        return err;
    }

    n1 = (*targ)[0];
    for (i=1; i<=src[0]; i++) {
        if (!in_gretl_list(*targ, src[i])) {
            n2++;
        }
    }
    if (n2 == 0) {
        /* nothing to be appended */
        free(scpy);
        return 0;
    }

    tmp = cast(int *) realloc(*targ, (n1 + n2 + 1) * (*tmp).sizeof);
    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        j = n1 + 1;
        for (i=1; i<=src[0]; i++) {
            if (!in_gretl_list(tmp, src[i])) {
                tmp[j++] = src[i];
            }
        }
        tmp[0] += n2;
        *targ = tmp;
    }

    free(scpy);

    return err;
}

/**
 * gretl_list_intersection:
 * @l1: list of integers.
 * @l2: list of integers.
 * @err: location to receive error code.
 *
 * Creates a list holding the intersection of @l1 and @l2.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_intersection(const(int) * l1, const(int) * l2, int * err)
{
    int * ret = null;
    int i, j, n = 0;

    if (l1 == null || l2 == null) {
        *err = GretlError.invarg;
        return null;
    }

    for (i=1; i<=l1[0]; i++) {
        if (in_gretl_list(l2, l1[i])) {
            n++;
        }
    }

    if (n == 0) {
        ret = gretl_null_list();
    } else {
        ret = gretl_list_new(n);
        if (ret != null) {
            j = 1;
            for (i=1; i<=l1[0]; i++) {
                if (in_gretl_list(l2, l1[i])) {
                    ret[j++] = l1[i];
                }
            }
        }
    }

    if (ret == null) {
        *err = GretlError.alloc;
    }

    return ret;
}

private void name_xprod_term(char * vname, int vi, int vj,
                             int di, const(DatasetRaw) * dset)
{
    const(char) * si = dset.varname[vi];
    const(char) * sj = dset.varname[vj];
    int ilen = cast(int) strlen(si);
    int jlen = cast(int) strlen(sj);
    int totlen = ilen + jlen + 2;
    char[16] numstr;

    sprintf(numstr.ptr, "%d", di);
    totlen += cast(int) strlen(numstr.ptr);

    if (totlen >= VNAMELEN) {
        int decr = 1 + totlen - VNAMELEN;

        while (decr > 0) {
            if (ilen > jlen) {
                ilen--;
            } else {
                jlen--;
            }
            decr--;
        }
    }

    sprintf(vname, "%.*s_%.*s_%s", ilen, si, jlen, sj, numstr.ptr);
}

private void set_xprod_label(int v, int vi, int vj,
                             double val, DatasetRaw * dset)
{
    const(char) * si = dset.varname[vi];
    const(char) * sj = dset.varname[vj];
    char[MAXLABEL] label;

    sprintf(label.ptr, "interaction of %s and (%s == %g)", si, sj, val);
    series_record_label(dset, v, label.ptr);
}

private int nonneg_integer_series(const(DatasetRaw) * dset, int v)
{
    double xt;
    int t;

    for (t=dset.t1; t<=dset.t2; t++) {
        xt = dset.Z[v][t];
        if (!na(xt) && (xt != floor(xt) || xt < 0)) {
            return 0;
        }
    }

    return 1;
}

/**
 * gretl_list_product:
 * @X: list of integers (representing discrete variables).
 * @Y: list of integers.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Creates a list holding the Cartesian product of @X and @Y.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_product(const(int) * X, const(int) * Y,
                         DatasetRaw * dset, int * err)
{
    int * ret = null;
    GretlMatrixRaw * xvals;
    char[VNAMELEN] vname;
    const(double) * x;
    const(double) * y;
    int * x_is_int = null;
    int newv, n_old = 0;
    int n, vi, vj;
    int i, j, k, t;

    if (X == null || Y == null) {
        *err = GretlError.data;
        return null;
    }

    if (X[0] == 0 || Y[0] == 0) {
        ret = gretl_null_list();
        if (ret == null) {
            *err = GretlError.alloc;
        }
        return ret;
    }

    x_is_int = gretl_list_new(X[0]);
    if (x_is_int == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* check the X list for discreteness */

    for (j=1; j<=X[0] && !*err; j++) {
        vj = X[j];
        if (nonneg_integer_series(dset, vj)) {
            x_is_int[j] = 1;
        } else if (!series_is_discrete(dset, vj)) {
            gretl_errmsg_sprintf(_("The variable '%s' is not discrete").ptr,
                                 dset.varname[vj]);
            *err = GretlError.data;
        }
    }

    if (*err) {
        free(x_is_int);
        return null;
    }

    n = sample_size(dset);
    newv = dset.v;

    for (j=1; j<=X[0] && !*err; j++) {
        vj = X[j];
        x = dset.Z[vj];
        xvals = gretl_matrix_values(x + dset.t1, n, Gretlopt.OPT_S, err);
        if (!*err) {
            *err = dataset_add_series(dset, Y[0] * xvals.rows);
            if (!*err) {
                for (i=1; i<=Y[0] && !*err; i++) {
                    vi = Y[i];
                    y = dset.Z[vi];
                    for (k=0; k<xvals.rows && !*err; k++) {
                        int v, oldv, iik;
                        double xik;

                        xik = gretl_vector_get(xvals, k);
                        iik = x_is_int[j] ? cast(int) xik : (k + 1);
                        name_xprod_term(vname.ptr, vi, vj, iik, dset);
                        oldv = current_series_index(dset, vname.ptr);
                        if (oldv > 0) {
                            /* reuse existing series of the same name */
                            v = oldv;
                            n_old++;
                        } else {
                            /* make a new series */
                            v = newv++;
                        }
                        for (t=dset.t1; t<=dset.t2; t++) {
                            if (na(x[t]) || na(xik)) {
                                dset.Z[v][t] = NADBL;
                            } else {
                                dset.Z[v][t] = (x[t] == xik)? y[t] : 0;
                            }
                        }
                        gretl_list_append_term(&ret, v);
                        if (ret == null) {
                            *err = GretlError.alloc;
                        } else {
                            if (v != oldv) {
                                strcpy(dset.varname[v], vname.ptr);
                            }
                            set_xprod_label(v, vi, vj, xik, dset);
                        }
                    }
                }
            }
            gretl_matrix_free(xvals);
        }
    }

    free(x_is_int);

    if (n_old > 0) {
        /* we added more series than were actually needed */
        dataset_drop_last_variables(dset, n_old);
    }

    return ret;
}

/**
 * gretl_list_omit_last:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @err: location to receive error code.
 *
 * Creates a list containing all but the last element of @orig,
 * which must not contain #LISTSEP and must contain at least
 * two members.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_omit_last(const(int) * orig, int * err)
{
    int * list = null;
    int n = orig[0];
    int i;

    *err = 0;

    if (n < 2) {
        *err = GretlError.novars;
    } else {
        for (i=1; i<=n; i++) {
            if (orig[i] == LISTSEP) {
                /* can't handle compound lists */
                *err = 1;
                break;
            }
        }
    }

    if (!*err) {
        list = gretl_list_new(n - 1);
        if (list == null) {
            *err = GretlError.alloc;
        } else {
            for (i=1; i<n; i++) {
                list[i] = orig[i];
            }
        }
    }

    return list;
}

private int list_count(const(int) * list)
{
    int i, k = 0;

    if (list != null) {
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                break;
            } else {
                k++;
            }
        }
    }

    return k;
}

/**
 * gretl_list_omit:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @omit: list of variables to drop.
 * @minpos: minimum position to check. This should be 2 for a regular
 * regression list, to skip the dependent var in position 1; but in
 * other contexts it may be 1 to start from the first element of @orig.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @omit.  It is an error if the @omit list contains
 * members that are not present in @orig.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_omit(const(int) * orig, const(int) * omit,
                      int minpos, int * err)
{
    int n_omit = omit[0];
    int n_orig = list_count(orig);
    int n_rem = n_orig - n_omit;
    int * ret = null;
    int i, j;

    if (n_rem < 0) {
        *err = GretlError.data;
        return null;
    }

    *err = 0;

    /* check for spurious "omissions" */
    for (i=1; i<=omit[0]; i++) {
        j = in_gretl_list(orig, omit[i]);
        if (j < minpos) {
            gretl_errmsg_sprintf(_("Variable %d was not in the original list").ptr,
                                 omit[i]);
            *err = 1;
            return null;
        }
    }

    ret = gretl_list_new(n_rem);

    if (ret == null) {
        *err = GretlError.alloc;
    } else if (n_rem > 0) {
        j = 1;
        for (i=1; i<=n_orig; i++) {
            if (i < minpos) {
                ret[j++] = orig[i];
            } else if (!in_gretl_list(omit, orig[i])) {
                ret[j++] = orig[i];
            }
        }
    }

    return ret;
}

/**
 * gretl_list_drop:
 * @orig: an array of integers, the first element of which holds
 * a count of the number of elements following.
 * @drop: list of variables to drop.
 * @err: pointer to receive error code.
 *
 * Creates a list containing the elements of @orig that are not
 * present in @drop.  Unlike gretl_list_omit(), processing always
 * starts from position 1 in @orig, and it is not an error if
 * some members of @drop are not present in @orig.
 *
 * Returns: new list on success, NULL on error.
 */

int * gretl_list_drop(const(int) * orig, const(int) * drop, int * err)
{
    int * lcopy = null;
    int * ret = null;
    int n_omit = 0;
    int i, j;

    *err = 0;

    lcopy = gretl_list_copy(orig);

    if (lcopy == null) {
        *err = GretlError.alloc;
        return null;
    } else {
        /* determine how many terms we're omitting */
        for (i=1; i<=drop[0]; i++) {
            j = in_gretl_list(lcopy, drop[i]);
            if (j > 0) {
                n_omit++;
                lcopy[j] = -1;
            }
        }
    }

    if (n_omit == 0) {
        ret = lcopy;
    } else {
        int n_rem = orig[0] - n_omit;

        ret = gretl_list_new(n_rem);
        if (ret == null) {
            *err = GretlError.alloc;
        } else if (n_rem > 0) {
            j = 1;
            for (i=1; i<=orig[0]; i++) {
                if (lcopy[i] >= 0) {
                    ret[j++] = orig[i];
                }
            }
        }
        free(lcopy);
    }

    return ret;
}

/**
 * gretl_list_diff_new:
 * @biglist: inclusive list.
 * @sublist: subset of biglist.
 * @minpos: position in lists at which to start.
 *
 * Returns: a newly allocated list including the elements of @biglist,
 * from position @minpos onwards, that are not present in @sublist,
 * again from @minpos onwards, or NULL on failure.  Note that
 * comparison stops whenever a list separator is found; i.e. only
 * the pre-separator portions of the lists are compared.
 */

int * gretl_list_diff_new(const(int) * biglist, const(int) * sublist,
                          int minpos)
{
    int * targ = null;
    int i, j, bi;
    int match;

    if (biglist == null || sublist == null) {
        return null;
    }

    targ = gretl_null_list();
    if (targ == null) {
        return null;
    }

    for (i=minpos; i<=biglist[0]; i++) {
        bi = biglist[i];
        if (bi == LISTSEP) {
            break;
        }
        match = 0;
        for (j=minpos; j<=sublist[0]; j++) {
            if (sublist[j] == LISTSEP) {
                break;
            } else if (sublist[j] == bi) {
                match = 1;
                break;
            }
        }
        if (!match) {
            /* but is this var already accounted for? */
            for (j=1; j<=targ[0]; j++) {
                if (targ[j] == bi) {
                    match = 1;
                    break;
                }
            }
        }
        if (!match) {
            targ = gretl_list_append_term(&targ, biglist[i]);
            if (targ == null) {
                break;
            }
        }
    }

    return targ;
}

/**
 * gretl_list_append_list:
 * @targ: location of list to which @src should be added.
 * @src: list to be added to @targ.
 *
 * Adds @src onto the end of @targ.  The length of @targ becomes the
 * sum of the lengths of the two original lists. Note that the result
 * may contain duplicated elements.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_append_list(int ** targ, const(int) * src)
{
    int * tmp;
    int i, n1, n2;
    int err = 0;

    if (targ == null || *targ == null) {
        return GretlError.invarg;
    } else if (src == null || src[0] == 0) {
        /* no-op */
        return 0;
    }

    n1 = (*targ)[0];
    n2 = src[0];

    tmp = cast(int *) realloc(*targ, (n1 + n2 + 1) * (*tmp).sizeof);

    if (tmp == null) {
        err = GretlError.alloc;
    } else {
        tmp[0] = n1 + n2;
        for (i=1; i<=src[0]; i++) {
            tmp[n1 + i] = src[i];
        }
        *targ = tmp;
    }

    return err;
}

/**
 * gretl_list_insert_list:
 * @targ: location of list into which @src should be inserted.
 * @src: list to be inserted.
 * @pos: zero-based position at which @src should be inserted.
 *
 * Inserts @src into @targ at @pos.  The length of @targ becomes the
 * sum of the lengths of the two original lists.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_insert_list(int ** targ, const(int) * src, int pos)
{
    int * big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int bign = n1 + n2;
    int i, err = 0;

    if (pos > n1 + 1) {
        return 1;
    }

    big = cast(int *) realloc(*targ, (bign + 1) * (*big).sizeof);

    if (big == null) {
        err = GretlError.alloc;
    } else {
        big[0] = bign;
        for (i=bign; i>=pos+n2; i--) {
            big[i] = big[i-n2];
        }
        for (i=1; i<=src[0]; i++) {
            big[pos+i-1] = src[i];
        }
        *targ = big;
    }

    return err;
}

/**
 * gretl_list_insert_list_minus:
 * @targ: location of list into which @src should be inserted.
 * @src: list to be inserted.
 * @pos: zero-based position at which @src should be inserted.
 *
 * Inserts @src into @targ at @pos.  The length of @targ becomes the
 * sum of the lengths of the two original lists minus one.  This
 * can be useful if we were expecting to insert a single variable
 * but found we had to insert a list instead.  Insertion of @src
 * overwrites any entries in @targ beyond @pos (the expectation is
 * that this function will be called in the process of assembling
 * @targ, in left-to-right mode).
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_list_insert_list_minus(int ** targ, const(int) * src, int pos)
{
    int * big;
    int n1 = (*targ)[0];
    int n2 = src[0];
    int bign = n1 - 1 + n2;
    int i, err = 0;

    if (pos > n1 + 1) {
        return 1;
    }

    big = cast(int *) realloc(*targ, (bign + 1) * (*big).sizeof);
    if (big == null) {
        err = GretlError.alloc;
    } else {
        big[0] = bign;
        for (i=1; i<=src[0]; i++) {
            big[pos+i-1] = src[i];
        }
        *targ = big;
    }

    return err;
}

/**
 * gretl_list_sublist:
 * @list: the source list.
 * @pos0: the first position.
 * @pos1: the last position, or -1 to continue to the end.
 *
 * Returns: a newly allocated sublist containing elements @pos0
 * to @pos1 of the source.
 */

int * gretl_list_sublist(const(int) * list, int pos0, int pos1)
{
    int * ret = null;
    int n;

    if (list == null) {
        return null;
    }

    if (pos1 < 0) {
        pos1 = list[0];
    }
    n = pos1 - pos0 + 1;

    if (n >= 0) {
        ret = gretl_list_new(n);
        if (n > 0 && ret != null) {
            int i, j = 1;

            for (i=pos0; i<=pos1; i++) {
                ret[j++] = list[i];
            }
        }
    }

    return ret;
}

/**
 * gretl_list_select:
 * @list: the source list.
 * @sel: the selection list.
 *
 * Returns: a newly allocated list containing the elements
 * of @list that are selected by @sel.
 */

int * gretl_list_select(const(int) * list, const(int) * sel, int * err)
{
    int * ret = gretl_list_new(sel[0]);
    int idxmax = list[0];
    int i, idx;

    for (i=1; i<=sel[0]; i++) {
        idx = sel[i];
        if (idx > 0 && idx <= idxmax) {
            ret[i] = list[idx];
        } else {
            *err = GretlError.invarg;
            break;
        }
    }

    if (*err) {
        free(ret);
        ret = null;
    }

    return ret;
}

/**
 * list_members_replaced:
 * @pmod: the model whose list is to be tested.
 * @dset: dataset information.
 *
 * Checks whether any variable used in @pmod has been redefined
 * since the model in question was estimated.
 *
 * Returns: non-zero if any variables have been replaced, 0 otherwise.
 */

int list_members_replaced(const(Model) * pmod, const(DatasetRaw) * dset)
{
    string errmsg = N_("Can't do this: some vars in original "
                       ~ "model have been redefined");
    int i, vi;

    if (pmod.list == null) {
        return 0;
    }

    for (i=1; i<=pmod.list[0]; i++) {
        vi = pmod.list[i];
        if (vi == LISTSEP) {
            continue;
        }
        if (vi >= dset.v) {
            gretl_errmsg_set(_(errmsg).ptr);
            return GretlError.data;
        }
        if (series_get_mtime(dset, vi) > pmod.esttime) {
            gretl_errmsg_set(_(errmsg).ptr);
            return GretlError.data;
        }
    }

    return 0;
}

/**
 * gretl_list_const_pos:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @minpos: position in @list at which to start the search (>= 1).
 * @dset: dataset struct.
 *
 * Checks @list for the presence, in position @minpos or higher, of
 * a variable whose valid values in sample all equal 1.  This usually
 * amounts to checking whether a list of regressors includes
 * an intercept term.
 *
 * Returns: The list position of the const, or 0 if none is
 * found.
 */

int gretl_list_const_pos(const(int) * list, int minpos,
                         const(DatasetRaw) * dset)
{
    int i;

    if (minpos < 1) {
        return 0;
    }

    /* we give preference to the "official" const... */
    for (i=minpos; i<=list[0]; i++) {
        if (list[i] == 0) {
            return i;
        }
    }

    /* ... but if it's not found */
    for (i=minpos; i<=list[0]; i++) {
        if (true_const(list[i], dset)) {
            return i;
        }
    }

    return 0;
}

/**
 * gretl_list_separator_position:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: if @list contains the separator for compound
 * lists, #LISTSEP, the position in @list at which this is found,
 * else 0.  The search begins at position 1.
 */

int gretl_list_separator_position(const(int) * list)
{
    int i;

    if (list != null) {
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                return i;
            }
        }
    }

    return 0;
}

/**
 * gretl_list_has_separator:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 *
 * Returns: 1 if @list contains the separator for compound
 * lists, #LISTSEP, else 0.  The search begins at position 1.
 */

int gretl_list_has_separator(const(int) * list)
{
    return gretl_list_separator_position(list) > 0;
}

/**
 * gretl_list_split_on_separator:
 * @list: source list.
 * @plist1: pointer to accept first sub-list, or NULL.
 * @plist2: pointer to accept second sub-list, or NULL.
 *
 * If @list contains the list separator, #LISTSEP, creates two
 * sub-lists, one containing the elements of @list preceding
 * the separator and one containing the elements following
 * the separator.  The sub-lists are newly allocated, and assigned
 * as the content of @plist1 and @plist2 respectively. Note, however,
 * that one or other of the sublists can be discarded by passing
 * NULL as the second or third argument.
 *
 * Returns: 0 on success, %E_ALLOC is memory allocation fails,
 * or %E_DATA if @list does not contain a separator.
 */

int gretl_list_split_on_separator(const(int) * list,
                                  int ** plist1,
                                  int ** plist2)
{
    int * list1 = null;
    int * list2 = null;
    int i, n = 0;

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            n = i;
            break;
        }
    }

    if (n == 0) {
        return GretlError.parse;
    }

    if (plist1 != null) {
        if (n > 1) {
            list1 = gretl_list_new(n - 1);
            if (list1 == null) {
                return GretlError.alloc;
            }
            for (i=1; i<n; i++) {
                list1[i] = list[i];
            }
        }
        *plist1 = list1;
    }

    if (plist2 != null) {
        if (n < list[0]) {
            list2 = gretl_list_new(list[0] - n);
            if (list2 == null) {
                free(list1);
                return GretlError.alloc;
            }
            for (i=1; i<=list2[0]; i++) {
                list2[i] = list[i + n];
            }
        }
        *plist2 = list2;
    }

    return 0;
}

/**
 * gretl_lists_join_with_separator:
 * @list1: first sub-list.
 * @list2: second sub-list.
 *
 * Concatenates the content of @list2 onto @list1, after first
 * appending the list separator.  It is acceptable that @list1
 * be NULL, in which case the returned list is just @list2
 * with the separator prepended.  But it is not acceptable that
 * @list2 be null; in that this function returns NULL.
 *
 * Returns: allocated list on success or NULL on failure.
 */

int * gretl_lists_join_with_separator(const(int) * list1, const(int) * list2)
{
    int * biglist;
    int i, j, n;

    if (list2 == null) {
        return null;
    }

    n = (list1 != null)? list1[0] : 0;


    n += list2[0] + 1;
    biglist = gretl_list_new(n);

    if (biglist == null) {
        return null;
    }

    j = 1;

    if (list1 != null) {
        for (i=1; i<=list1[0]; i++) {
            biglist[j++] = list1[i];
        }
    }

    biglist[j++] = LISTSEP;

    for (i=1; i<=list2[0]; i++) {
        biglist[j++] = list2[i];
    }

    return biglist;
}

private int real_list_dup(const(int) * list, int start, int stop)
{
    int i, j, ret = -1;

    for (i=start; i<stop && ret<0; i++) {
        for (j=i+1; j<=stop && ret<0; j++) {
            if (list[i] == list[j]) {
                ret = list[i];
            }
        }
    }

    return ret;
}

/**
 * gretl_list_duplicates:
 * @list: an array of integer variable ID numbers, the first element
 * of which holds a count of the number of elements following.
 * @ci: index of gretl command (for context).
 *
 * Checks whether or not a gretl list contains duplicated elements.
 * Exactly what counts as duplication depends on the context of the
 * command in which @list will be used, which is given by @ci.
 *
 * Returns: the ID number of the first duplicated variable found,
 * or -1 in case of no duplication.
 */

int gretl_list_duplicates(const(int) * list, GretlCmdIndex ci)
{
    int done = 0;
    int start = 2; /* perhaps should be 1? */
    int i, ret = -1;

    if (ci == GretlCmdIndex.coint || ci == GretlCmdIndex.anova ||
        ci == GretlCmdIndex.deleet) {
        start = 1;
    } else if (ci == GretlCmdIndex.arch) {
        start = 3;
    } else if (ci == GretlCmdIndex.arma) {
        for (i=list[0]-1; i>2; i--) {
            if (list[i] == LISTSEP) {
                start = i+1;
                break;
            }
        }
    } else if (ci == GretlCmdIndex.lags && list[0] > 1 && list[2] == LISTSEP) {
        start = 3;
    } else if (ci == GretlCmdIndex.ar || ci == GretlCmdIndex.scatters ||
               ci == GretlCmdIndex.mpols || ci == GretlCmdIndex.garch) {
        for (i=2; i<list[0]; i++) {
            if (list[i] == LISTSEP) {
                start = i+1;
                break;
            }
        }
    } else if (ci == GretlCmdIndex.ivreg || ci == GretlCmdIndex.heckit ||
               ci == GretlCmdIndex.equation) {
        for (i=2; i<list[0]; i++) {
            if (list[i] == LISTSEP) {
                start = i+1;
                break;
            }
        }
        ret = real_list_dup(list, start, list[0]);
        if (ret == -1) {
            ret = real_list_dup(list, 2, start - 2);
        }
        done = 1;
    } else if (ci == GretlCmdIndex.dpanel) {
        int stop = 0;

        for (i=2; i<list[0]; i++) {
            if (list[i] == LISTSEP) {
                start = i;
                break;
            }
        }
        for (i=list[0]-1; i>=2; i--) {
            if (list[i] == LISTSEP) {
                stop = i;
                break;
            }
        }
        if (stop == start) {
            ret = real_list_dup(list, start + 1, list[0]);
        } else {
            ret = real_list_dup(list, start + 1, stop - 1);
            if (ret == -1) {
                ret = real_list_dup(list, stop + 1, list[0]);
            }
        }
        done = 1;
    } else if (ci == GretlCmdIndex.biprobit) {
        if (list[1] == list[2]) {
            ret = 1;
        }
        if (ret == -1) {
            for (i=3; i<list[0]; i++) {
                if (list[i] == LISTSEP) {
                    start = i+1;
                    break;
                }
            }
            ret = real_list_dup(list, start, list[0]);
            if (ret == -1) {
                ret = real_list_dup(list, 3, start - 2);
            }
        }
        done = 1;
    }

    if (!done) {
        ret = real_list_dup(list, start, list[0]);
    }

    return ret;
}

/**
 * gretl_lists_share_members:
 * @list1:
 * @list2:
 *
 * Returns: the number of elements that are in common between
 * @list1 and @list2.
 */

int gretl_lists_share_members(const(int) * list1, const(int) * list2)
{
    int i, n = 0;

    if (list1 != null && list2 != null) {
        for (i=1; i<=list1[0]; i++) {
            if (in_gretl_list(list2, list1[i])) {
                n++;
            }
        }
    }

    return n;
}

/**
 * gretl_list_n_distinct_members:
 * @list: list to test.
 *
 * Returns: the count of distinct elements in list from position
 * 1 onward, not counting #LISTSEP if present.
 */

int gretl_list_n_distinct_members(const(int) * list)
{
    int i, j, n = list[0];

    for (i=1; i<=list[0]; i++) {
        if (list[i] == LISTSEP) {
            n--;
        } else {
            for (j=2; j<i; j++) {
                if (list[i] == list[j]) {
                    n--;
                    break;
                }
            }
        }
    }

    return n;
}

/**
 * full_var_list:
 * @dset: dataset information.
 * @nvars: location for return of number of elements in full list.
 *
 * Creates a newly allocated list including all series in the
 * dataset that are not hidden variables, and are accessible
 * at the current level of function execution.
 * The return value is NULL in case either (a) allocation of
 * memory failed, or (b) the resulting list would be empty.
 * The caller can distinguish between these possibilities by
 * examining the value returned in @nvars, which will be zero if
 * and only if the resulting list would be empty.  If this is
 * not of interest to the caller, @nvars may be given as NULL.
 *
 * Returns: the allocated list, or NULL.
 */

int * full_var_list(const(DatasetRaw) * dset, int * nvars)
{
    int fsd = gretl_function_depth();
    int i, j, nv = 0;
    int * list = null;

    if (dset == null) {
        if (nvars != null) {
            *nvars = 0;
        }
        return null;
    }

    for (i=1; i<dset.v; i++) {
        if (!series_is_hidden(dset, i) &&
            series_get_stack_level(dset, i) == fsd) {
            nv++;
        }
    }

    if (nvars != null) {
        *nvars = nv;
    }

    if (nv > 0) {
        list = gretl_list_new(nv);
    }

    if (list != null) {
        j = 1;
        for (i=1; i<dset.v; i++) {
            if (!series_is_hidden(dset, i) &&
                series_get_stack_level(dset, i) == fsd) {
                list[j++] = i;
            }
        }
    }

    return list;
}

/**
 * gretl_list_is_consecutive:
 * @list: list to check.
 *
 * Returns: 1 if the elements of @list, from position 1 onward,
 * are consecutive integer values, else 0.
 */

int gretl_list_is_consecutive(const(int) * list)
{
    int i, ret = 1;

    for (i=2; i<=list[0]; i++) {
        if (list[i] != list[i-1] + 1) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_list_build:
 * @s: string list specification.
 * @dset: dataset information.
 * @err: location to receive error code
 *
 * Builds a list based on the specification in @s, which may include
 * the ID numbers of variables, the names of variables, and/or the
 * names of previously defined lists (all separated by spaces).
 *
 * Returns: the constructed list, or NULL on failure.
 */

int * gretl_list_build(const(char) * s, const(DatasetRaw) * dset, int * err)
{
    char[32] test;
    int * list = null;
    int * nlist;
    int i, v, len, nf;

    list = gretl_null_list();
    if (list == null) {
        *err = GretlError.alloc;
        return null;
    }

    nf = count_fields(s, null);

    for (i=0; i<nf && !*err; i++) {
        s += strspn(s, " ");
        len = cast(int) strcspn(s, " ");
        if (len > 31) {
            *err = GretlError.parse;
        } else {
            test[0] = 0;
            strncat(test.ptr, s, len);

            /* valid elements: integers, varnames, named lists */

            if (isdigit(test[0])) {
                v = positive_int_from_string(test.ptr);
                if (v >= 0) {
                    list = gretl_list_append_term(&list, v);
                } else {
                    *err = GretlError.parse;
                }
            } else {
                v = series_index(dset, test.ptr);
                if (v < dset.v) {
                    list = gretl_list_append_term(&list, v);
                } else {
                    nlist = get_list_by_name(test.ptr);
                    if (nlist != null) {
                        *err = gretl_list_append_list(&list, nlist);
                    } else {
                        *err = GretlError.unkvar;
                    }
                }
            }

            if (list == null) {
                *err = GretlError.alloc;
            }
        }
        s += len;
    }

    if (*err) {
        free(list);
        list = null;
    }

    return list;
}

/**
 * gretl_list_print:
 * @list: list to print.
 * @dset: dataset information.
 * @prn: gretl printing struct.
 *
 * Prints to @prn the given @list of variables, by name
 * if @dset is non-NULL otherwise by ID number.
 */

void gretl_list_print(const(int) * list, const(DatasetRaw) * dset,
                      Prn * prn)
{
    int testlen = 62;
    int i, li, len = 0;

    if (list == null) {
        pputs(prn, "null\n");
    } else if (list[0] == 0) {
        pputs(prn, "empty\n");
    } else {
        for (i=1; i<=list[0]; i++) {
            li = list[i];
            if (li == LISTSEP) {
                len += pputs(prn, "; ");
            } else if (dset == null) {
                len += pprintf(prn, "%d ", li);
            } else if (li < 0 || li >= dset.v) {
                len += pputs(prn, "?? ");
            } else {
                len += pprintf(prn, "%s ", dset.varname[li]);
                if (len > testlen && i < list[0]) {
                    pputs(prn, "\\\n ");
                    len = 1;
                }
            }
        }
        pputc(prn, '\n');
    }
}

/**
 * varname_match_list:
 * @dset: pointer to dataset information.
 * @pattern: pattern to be matched.
 * @err: location to receive error code.
 *
 * Returns: a list of ID numbers of variables whose names
 * match @pattern, or NULL if there are no matches.
 */

int * varname_match_list(const(DatasetRaw) * dset, const(char) * pattern,
                         int * err)
{
    GPatternSpec * pspec;
    int * list = null;
    int i, fd, n = 0;

    if (dset == null || dset.v == 0) {
        return null;
    }

    fd = gretl_function_depth();

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<dset.v; i++) {
        if (fd == 0 || fd == series_get_stack_level(dset, i)) {
            if (g_pattern_match_string(pspec, dset.varname[i])) {
                n++;
            }
        }
    }

    if (n > 0) {
        list = gretl_list_new(n);
        if (list == null) {
            *err = GretlError.alloc;
        } else {
            int j = 1;

            for (i=1; i<dset.v; i++) {
                if (fd == 0 || fd == series_get_stack_level(dset, i)) {
                    if (g_pattern_match_string(pspec, dset.varname[i])) {
                        list[j++] = i;
                    }
                }
            }
        }
    }

    g_pattern_spec_free(pspec);

    return list;
}

/**
 * ellipsis_list:
 * @dset: pointer to dataset information.
 * @v1: index of first variable.
 * @v2: index of last variable.
 * @err: location to receive error code.
 *
 * Returns: a list of ID numbers of variables running
 * from @v1 to @v2.
 */

int * ellipsis_list(const(DatasetRaw) * dset, int v1, int v2, int * err)
{
    int * list = null;
    int i, fd, n = 0;

    if (dset == null || dset.v == 0) {
        return null;
    }

    fd = gretl_function_depth();

    for (i=v1; i<=v2; i++) {
        if (fd == 0 || fd == series_get_stack_level(dset, i)) {
            n++;
        }
    }

    if (n > 0) {
        list = gretl_list_new(n);
        if (list == null) {
            *err = GretlError.alloc;
        } else {
            int j = 1;

            for (i=v1; i<=v2; i++) {
                if (fd == 0 || fd == series_get_stack_level(dset, i)) {
                    list[j++] = i;
                }
            }
        }
    }

    return list;
}

/**
 * gretl_list_from_vector:
 * @v: matrix (must be a vector).
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Tries to interpret the matrix @v as a list of ID
 * numbers of series. This can work only if @v is a
 * vector, and all its elements have integer values
 * k satisfying 0 <= k < v, where v is the number
 * of series in @dset. In the special case where @v
 * is a null matrix, an empty list is returned.
 *
 * Returns: a gretl list, or NULL on failure.
 */

int * gretl_list_from_vector(const(GretlMatrixRaw) * v,
                             const(DatasetRaw) * dset,
                             int * err)
{
    int * list = null;

    if (gretl_is_null_matrix(v)) {
        list = gretl_null_list();
        if (list == null) {
            *err = GretlError.alloc;
        }
    } else {
        int i, vi, k = gretl_vector_get_length(v);

        if (k == 0) {
            *err = GretlError.types;
        } else {
            for (i=0; i<k && !*err; i++) {
                vi = gretl_int_from_double(v.val[i], err);
                if (!*err && (vi >= dset.v || (vi < 0 && vi != LISTSEP))) {
                    gretl_errmsg_sprintf("list from vector: series ID %d "
                                         ~ "is out of bounds", vi);
                    *err = GretlError.unkvar;
                }
            }
            if (!*err) {
                list = gretl_list_new(k);
                if (list == null) {
                    *err = GretlError.alloc;
                } else {
                    for (i=0; i<k; i++) {
                        list[i+1] = cast(int) v.val[i];
                    }
                }
            }
        }
    }

    return list;
}

private int colnames_ok_for_series(const(char) ** S, int n)
{
    int i, j, err = 0;

    for (i=0; i<n; i++) {
        err = check_varname(S[i]);
        if (!err && gretl_is_user_var(S[i])) {
            gretl_errmsg_sprintf("'%s': name conflicts with different type", S[i]);
            err = GretlError.types;
        } else {
            for (j=0; j<n; j++) {
                if (j != i && strcmp(S[j], S[i]) == 0) {
                    gretl_errmsg_sprintf("'%s': name is not unique", S[i]);
                    err = 1;
                }
            }
        }
        if (err) {
            return 0;
        }
    }

    return 1;
}

private int try_list_vname(char * chkname,
                           const(char) * pfx,
                           int j, int n)
{
    char * tmp = g_strdup_printf("%s%0*d", pfx, n, j);
    int err = check_varname(tmp);

    if (!err && gretl_is_user_var(tmp)) {
        gretl_errmsg_sprintf("'%s': name conflicts with different type", tmp);
        err = GretlError.types;
    }
    if (!err) {
        strcpy(chkname, tmp);
    }
    g_free(tmp);

    return err;
}

/**
 * gretl_list_from_matrix:
 * @X: matrix.
 * @prefix: prefix for series names, or NULL.
 * @dset: pointer to dataset.
 * @err: location to receive error code.
 *
 * Tries to add the columns of @X to @dset as series, and
 * if successful constructs a list holding the added series.
 * numbers of series. This is possible only if the rows of
 * @X can be matched to observations in @dset.
 *
 * Returns: allocated list, or NULL on failure.
 */

int * gretl_list_from_matrix(const(GretlMatrixRaw) * X,
                             const(char) * prefix,
                             DatasetRaw * dset,
                             int * err)
{
    int * list = null;
    const(char) ** S = null;
    int orig_v = dset.v;

    if (gretl_is_null_matrix(X)) {
        list = gretl_null_list();
        if (list == null) {
            *err = GretlError.alloc;
        }
        return list;
    }

    S = gretl_matrix_get_colnames(X);
    if (S == null && prefix == null) {
        *err = GretlError.invarg;
    } else if (S != null) {
        /* check the colnames */
        if (!colnames_ok_for_series(S, X.cols)) {
            *err = GretlError.invarg;
        }
    } else {
        /* check the prefix */
        if (strlen(prefix) > VNAMELEN - 3) {
            *err = GretlError.invarg;
        }
    }

    if (!*err) {
        char[VNAMELEN] chkname;
        int mt1 = gretl_matrix_get_t1(X);
        int mt2 = gretl_matrix_get_t2(X);
        int s1 = 0, s2 = 0;
        int n_add = X.cols;
        int j, slen = 0;
        int ok = 0;

        if (mt2 > 0) {
            if (X.rows != mt2 - mt1 + 1) {
                ; /* time indices not usable! */
            } else if (mt2 < dset.t1 || mt1 > dset.t2) {
                ; /* no overlap of ranges */
            } else {
                s1 = dset.t1 - mt1;
                s2 = X.rows - s1;
                ok = 1;
            }
        } else if (X.rows == sample_size(dset)) {
            /* length matches current sample */
            s1 = 0;
            s2 = X.rows - 1;
            ok = 1;
        } else if (X.rows == dset.n) {
            /* length matches full series length */
            s1 = dset.t1;
            s2 = dset.t2;
            ok = 1;
        }
        if (!ok) {
            gretl_errmsg_set("matrix to list: data ranges could not be matched");
            *err = GretlError.data;
        } else {
            list = gretl_list_new(X.cols);
            if (list == null) {
                *err = GretlError.alloc;
            }
        }

        if (S == null) {
            slen = cast(int) floor(log10(cast(double) X.cols)) + 1;
            if (slen > 6) {
                *err = GretlError.data;
            }
        }

        /* first pass, check the putative series names */
        for (j=0; j<X.cols && !*err; j++) {
            if (S != null) {
                strcpy(chkname.ptr, S[j]);
            } else {
                *err = try_list_vname(chkname.ptr, prefix, j+1, slen);
            }
            if (!*err && gretl_is_series(chkname.ptr, dset)) {
                /* an existing series, decrement the count
                   of series to be added to @dset
                */
                n_add--;
            }
        }

        if (!*err && n_add > 0) {
            *err = dataset_add_NA_series(dset, n_add);
        }
        if (!*err) {
            int vnew = orig_v;
            int vj, t, s;

            for (j=0; j<X.cols && !*err; j++) {
                if (S != null) {
                    strcpy(chkname.ptr, S[j]);
                } else if (prefix != null) {
                    sprintf(chkname.ptr, "%s%0*d", prefix, slen, j+1);
                }
                vj = current_series_index(dset, chkname.ptr);
                if (vj < 0) {
                    vj = vnew++;
                    strcpy(dset.varname[vj], chkname.ptr);
                }
                for (t=dset.t1, s=s1; t<=dset.t2 && s<=s2; t++, s++) {
                    if (s < 0) {
                        dset.Z[vj][t] = NADBL;
                    } else {
                        dset.Z[vj][t] = gretl_matrix_get(X, s, j);
                    }
                }
                list[j+1] = vj;
            }
        }
    }

    if (*err) {
        dataset_drop_last_variables(dset, dset.v > orig_v);
        free(list);
        list = null;
    }

    return list;
}

/**
 * gretl_auxlist_from_vector:
 * @v: source vector.
 * @err: location to receive error code.
 *
 * Unlike gretl_list_from_vector() this function does not
 * require that the elements of @v are valid series IDs. They
 * may represent a set of orders.
 *
 * Returns: a newly allocated gretl list containing the values
 * in @v, or NULL on failure. Note that it is an error if
 * @v is NULL, or is not a vector.
 */

int * gretl_auxlist_from_vector(const(GretlVector) * v, int * err)
{
    int i, n = gretl_vector_get_length(v);
    int * list = null;

    if (n == 0) {
        *err = GretlError.data;
    } else {
        list = gretl_list_new(n);
        if (list == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0; i<n; i++) {
                list[i+1] = cast(int) v.val[i];
            }
        }
    }

    return list;
}

/**
 * varname_match_any:
 * @dset: pointer to dataset information.
 * @pattern: pattern to be matched.
 *
 * Returns: 1 if at least one variable in the dataset has a
 * name that matches @pattern, otherwise 0.
 */

int varname_match_any(const(DatasetRaw) * dset, const(char) * pattern)
{
    GPatternSpec * pspec;
    int i, fd, ret = 0;

    fd = gretl_function_depth();

    pspec = g_pattern_spec_new(pattern);

    for (i=1; i<dset.v; i++) {
        if (fd == 0 || fd == series_get_stack_level(dset, i)) {
            if (g_pattern_match_string(pspec, dset.varname[i])) {
                ret = 1;
                break;
            }
        }
    }

    g_pattern_spec_free(pspec);

    return ret;
}
