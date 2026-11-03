/* Port of gretl_utils.c */

module gretl.gretl_utils;


import core.cpuid : avx, avx2, avx512f, processor, vendor;
import core.stdc.ctype : isalpha, isdigit;
import core.stdc.errno : EDOM, ERANGE, errno;
import core.stdc.math : ceil, fabs, floor, isnan, log, lrint, nearbyint;
import core.stdc.stdio : FILE, fclose, fputc, fputs, fprintf, fread, fscanf, fwrite, snprintf,
    sprintf, sscanf, stderr;
import core.stdc.stdlib : atoi, calloc, free, getenv, malloc, qsort, realloc, strtod, strtol;
import core.stdc.string : memcpy, strchr, strcmp, strcpy, strlen, strncat, strncmp, strrchr,
    strspn, strstr, strtok;
version (Posix) import core.sys.posix.dirent : DIR, closedir, dirent, opendir, readdir, rewinddir; else import gretl.winport : DIR, closedir, dirent, opendir, readdir, rewinddir;
version (Posix) import core.sys.posix.dlfcn : RTLD_NOW, dlopen, dlsym; else import gretl.winport : RTLD_NOW, dlopen, dlsym;
version (Posix) import core.sys.posix.sys.resource : RLIMIT_STACK, getrlimit, rlimit; else import gretl.winport : RLIMIT_STACK, getrlimit, rlimit;
version (Posix) import core.sys.posix.sys.stat : S_IFLNK, S_IFREG, S_IXGRP, S_IXOTH, S_IXUSR, stat_t; else import gretl.winport : S_IFLNK, S_IFREG, S_IXGRP, S_IXOTH, S_IXUSR, stat_t;
version (Posix) import core.sys.posix.unistd : getgid, getuid, gid_t, uid_t; else import gretl.winport : getgid, getuid, gid_t, uid_t;
import etc.c.zlib : gzFile, gzclose, gzread, gzwrite;

import gretl.dataset : dataset_destroy_obs_markers, series_is_discrete;
import gretl.dataseth : dataset_is_panel;
import gretl.gmem : g_free;
import gretl.glib : GDateWeekday, gint64;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_matrix : gretl_matrix_copy, gretl_matrix_free, gretl_zero_matrix_new, lapack_mem_free;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputs;
import gretl.gstrfuncs : g_strcmp0, g_strdup_printf, g_strstrip;
import gretl.transforms : gretl_transforms_cleanup;
import gretl.gerror : GError, g_error_free;
import gretl.gspawn : G_SPAWN_SEARCH_PATH, g_spawn_command_line_sync, g_spawn_sync;
import gretl.gtimer : g_get_monotonic_time;
import gretl.libgretlh;
import gretl.libset : gretl_function_depth, gretl_in_gui_mode, gretl_in_tool_mode, libset_cleanup, libset_init;
import gretl.missingh : NADBL, na;
import gretl.options : option_printing_cleanup, setopt_cleanup, stored_options_cleanup;
import gretl.strutils : dot_atof, gretl_strdup, gretl_strdup_printf, gretl_strndup, has_suffix;
import gretl.gretl_utilsh : BlasVariant, SessionClear, floatneq, ok_int;
import gretl.versionh : GRETL_VERSION, LIBGRETL_CURRENT;
import gretl.dataio : ntolabel;
import gretl.gretl_xml : gretl_xml_cleanup, gretl_xml_init;
import gretl.gretl_paths : gretl_chdir, gretl_deltree, gretl_dotdir, gretl_fopen, gretl_gnuplot_path, gretl_gzopen, gretl_isdir, gretl_remove, gretl_script_dirs_cleanup, gretl_stat, gretl_workdir;
import gretl.gretl_string_table : builtin_strings_cleanup;
import gretl.plugins : plugins_cleanup;
import gretl.calendar : calendar_date_string, get_epoch_day, weekday_from_epoch_day, ymd_extended_from_epoch_day;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.describe : last_result_cleanup;
import gretl.genmain : generate, gretl_is_series;
import gretl.gretl_panel : set_panel_group_strings, set_panel_structure_from_varnames, switch_panel_orientation;
import gretl.forecast : forecast_matrix_cleanup;
import gretl.printout : bufgets_cleanup;
import gretl.graphing : gnuplot_cleanup;
import gretl.gretl_bundle : gretl_bundle_cleanup;
import gretl.uservar : destroy_user_vars, get_scalar_value_by_name, gretl_is_user_var, gretl_lists_cleanup, user_var_delete_by_name, user_var_get_type_by_name;
import gretl.gretl_typemap : gretl_typemap_cleanup;
import gretl.genfuncs : gretl_plotx;
import gretl.gretl_commands : gretl_command_hash_cleanup;
import gretl.random : gretl_rand_init;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    void * get_function_package_by_name(const(char) * pkgname);
    const(char) * get_function_package_path_by_name(const(char) * pkgname);
    void function_package_unload_full_by_filename(const(char) * fname);
    int object_is_function_arg(const(char) * name);
    void num_threads_init(int blas_type);
    void gretl_saved_objects_cleanup();
    void gretl_functions_cleanup();
    void gretl_function_hash_cleanup();
}

private enum FN_NAMELEN = 32;     /* gretl_func.h */

/* The installation prefix is a configure-time value in the C build
   (GRETL_PREFIX); it is used only to guess the path of libgretl when
   probing for the BLAS variant via ldd. */

private enum GRETL_PREFIX = "/usr/local";

/**
 * date_as_double:
 * @t: observation number (zero-based).
 * @pd: data periodicity or frequency.
 * @sd0: floating point representation of starting date.
 *
 * Returns: the date corresponding to @t, as a double-precision number.
 */

double date_as_double(int t, int pd, double sd0)
{
    int ysd, yy, pp, yp;
    int p10 = 10;

    ysd = cast(int) sd0;

    if (pd == 1) {
        return cast(double) (ysd + t);
    }

    pp = pd;
    while ((pp = pp / 10) != 0) {
        p10 *= 10;
    }

    pp = cast(int) (t % pd + p10 * (sd0 - ysd) + .5);
    if (pp != pd)  {
        yy = cast(int) (ysd + t/pd  + pp/pd + .5);
        yp = pp % pd;
    }  else {
        yy = cast(int) (ysd + t/pd + .5);
        yp = pp;
    }

    return yy + cast(double) yp / p10;
}

/**
 * ijton:
 * @i: row number (0-based)
 * @j: column number (0-based)
 * @nrows: number of rows (and columns) in symmetric matrix.
 *
 * Given a (row, column) reference into a symmetric matrix A, finds
 * the index into a 1-dimensional array x composed of the
 * non-redundant (lower) elements of A.
 *
 * E.g. for the 3 x 3 case with 6 non-redundant elements, 0 to 5,
 *
 *    A(0,0) = x[0]  A(0,1) = x[1]  A(0,2) = x[2]
 *    A(1,0) = x[1]  A(1,1) = x[3]  A(1,2) = x[4]
 *    A(2,0) = x[2]  A(2,1) = x[4]  A(2,1) = x[5]
 *
 * Returns: 0-based index into flat array.
 */

int ijton(int i, int j, int nrows)
{
    if (i > j) {
        int tmp = i;

        i = j;
        j = tmp;
    }

    return nrows * i + j - i - ((i - 1) * i / 2);
}

/**
 * transcribe_array:
 * @targ: array to which to write.
 * @src: array from which to read.
 * @dset: data information struct.
 *
 * Copy from @src to @targ, skipping any missing values,
 * over the sample range defined in @dset.
 *
 * Returns: the number of valid observations put into @targ.
 */

int transcribe_array(double * targ, const(double) * src,
                     const(DatasetRaw) * dset)
{
    int t, n = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (!na(src[t])) {
            targ[n++] = src[t];
        }
    }

    return n;
}

/**
 * gretl_isdummy:
 * @t1: starting observation.
 * @t2: ending observation.
 * @x: data series to examine.
 *
 * Check whether series @x has only 0 or 1 values over the
 * given sample range (or possibly missing values).
 *
 * Returns: 0 if the variable is not a 0/1 dummy, otherwise the
 * number of 1s in the series.
 */

int gretl_isdummy(int t1, int t2, const(double) * x)
{
    int t, m = 0, goodobs = 0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (x[t] != 0.0 && x[t] != 1.0) {
            return 0;
        }
        if (x[t] == 1.0) {
            m++;
        }
        goodobs++;
    }

    if (m < goodobs) {
        return m;
    }

    return 0;
}

/**
 * gretl_iszero:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x has only zero (or missing) values
 * over the given sample range.
 *
 * Returns: 0 if the series contains any valid non-zero values,
 * otherwise 1.
 */

int gretl_iszero(int t1, int t2, const(double) * x)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != 0) {
            return 0;
        }
    }

    return 1;
}

/**
 * gretl_isconst:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x is constant over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if the variable is constant, otherwise 0.
 */

int gretl_isconst(int t1, int t2, const(double) * x)
{
    int t, ret = 1;

    while (na(x[t1]) && t1 <= t2) {
        t1++;
    }

    if (t1 >= t2) {
        return 0;
    }

    for (t=t1+1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (floatneq(x[t], x[t1])) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_isstoch:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x is stochastic, and contains a
 * contiguous set of valid values within the given sample span.
 * The simple and fallible heuristic for a stochastic series is
 * that the differences between successive values are all the
 * same.
 *
 * Returns: 1 if the variable appears to be stochastic on
 * the specified criterion, otherwise 0.
 */

int gretl_isstoch(int t1, int t2, const(double) * x)
{
    double dx0;
    int na_count = 0;
    int multidiff = 0;
    int s1 = -1, s2 = -1;
    int t;

    if (t1 >= t2) {
        return 0;
    }

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            s1 = t;
            break;
        }
    }
    for (t=t2; t>=s1; t--) {
        if (!na(x[t])) {
            s2 = t;
            break;
        }
    }
    if (s1 < 0 || s2 < 0 || s1 >= s2) {
        return 0;
    }

    dx0 = x[s1+1] - x[s1];

    for (t=s1; t<=s2; t++) {
        if (na(x[t])) {
            na_count++;
            break;
        }
        if (t > s1 && !multidiff) {
            if (x[t] - x[t-1] != dx0) {
                multidiff = 1;
            }
        }
    }

    return (na_count == 0 && multidiff);
}

/**
 * gretl_isunits:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x equals 1 over the
 * given sample range (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_isunits(int t1, int t2, const(double) * x)
{
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != 1.0) {
            ret = 0;
            break;
        }
    }

    return ret;
}

/**
 * gretl_iscount:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x contains nothing but non-negative
 * integer values (some of which are > 1) over the
 * given sample range.
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_iscount(int t1, int t2, const(double) * x)
{
    int t, xi;
    int g1 = 0;
    int ret = 1;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        if (x[t] < 0.0) {
            ret = 0;
            break;
        }
        xi = cast(int) x[t];
        if (x[t] != cast(double) xi) {
            ret = 0;
            break;
        }
        if (x[t] > 1.0) {
            g1 = 1;
        }
    }

    if (g1 == 0) {
        ret = 0;
    }

    return ret;
}

/**
 * gretl_isint:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Check whether series @x contains nothing but integer
 * values (disregarding NAs) over the given sample range.
 *
 * Returns: 1 if so, otherwise 0.
 */

int gretl_isint(int t1, int t2, const(double) * x)
{
    int t;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != floor(x[t])) {
            return 0;
        }
    }

    return 1;
}

private enum FEW_VALS = 32;
private enum FEWER_VALS = 8;

private int few_vals(int t1, int t2, const(double) * x, double * ratio)
{
    double[FEW_VALS] test;
    int match;
    int i, t, n = 0, nv = 0;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t])) {
            match = 0;
            for (i=0; i<nv; i++) {
                if (x[t] == test[i]) {
                    /* we've already seen this value */
                    match = 1;
                    break;
                }
            }
            if (!match) {
                /* x[t] is a "new" value */
                if (nv == FEW_VALS) {
                    /* hit the limit */
                    nv++;
                    break;
                }
                test[nv++] = x[t];
            }
            n++;
        }
    }

    /* ratio of distinct values to observations */
    *ratio = nv / cast(double) n;

    return nv;
}

/**
 * gretl_isdiscrete:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Checks the variable @x over the range @t1 to @t2 for discreteness.
 * This is a heuristic whose components are (a) whether the values
 * are "fairly round" (multiples of 0.25) or not, and, if test (a) is
 * passed, (b) whether the variable takes on only "few" distinct
 * values.
 *
 * Returns: 0 if test (a) is not passed or the number of distinct values
 * is > 32; else 1 if the number of distinct values is <= 32; else 2 if
 * the number of distinct values is <= 8.  A return of 1 is supposed
 * to indicate that it's "reasonable" to treat @x as discrete, while
 * a return of 2 indicates that it's probably unreasonable _not_ to
 * treat @x as discrete for the purpose of drawing up a frequency
 * distribution.
 */

int gretl_isdiscrete(int t1, int t2, const(double) * x)
{
    int t, n = 0, disc = 1;
    int allints = 1;
    int isconst = 1;
    double x0 = NADBL;
    double r = 0;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        } else if (na(x0)) {
            x0 = x[t];
        } else if (x[t] != x0) {
            isconst = 0;
        }
        n++;
        if (!ok_int(x[t])) {
            allints = 0;
        }
        r = x[t] - floor(x[t]);
        if (allints && r != 0) {
            allints = 0;
        }
        if (r != 0.0 && r != 0.25 && r != 0.5 && r != 0.75) {
            disc = 0;
            break;
        }
    }

    if (n == 0 || isconst) {
        /* we'll not treat a single-valued series as "discrete" */
        disc = 0;
    }

    if (disc) {
        n = few_vals(t1, t2, x, &r);
        if (allints) {
            if (n <= FEW_VALS && r < 0.2) {
                disc = 2;
            } else {
                disc = (n <= FEWER_VALS)? 2 : 1;
            }
        } else if (n > FEW_VALS) {
            disc = 0;
        } else if (r > 0.9 && n > 30) {
            /* somewhat arbitrary: but if r (= ratio of distinct
               values to number of cases) is "too high", and the
               number of observations is not tiny, perhaps we should
               not automatically take the var as discrete
            */
            disc = 0;
        } else if (n <= FEWER_VALS) {
            disc = 2;
        }
    }

    return disc;
}

int accept_as_discrete(const(DatasetRaw) * dset, int v, int strict)
{
    if (series_is_discrete(dset, v)) {
        /* the series has been explicitly marked as discrete */
        return 1;
    } else {
        /* check for plausibility of discreteness */
        int d = gretl_isdiscrete(dset.t1, dset.t2, dset.Z[v]);

        return strict ? d > 1 : d > 0;
    }
}

/**
 * gretl_ispositive:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 * @strict: boolean, strict inequality.
 *
 * Returns: 1 if all non-missing values in @x, over the range @t1 to
 * @t2, are positive (if strict==1) or non-negative (if strict==0),
 * otherwise 0.
 */

int gretl_ispositive(int t1, int t2, const(double) * x, int strict)
{
    int t;

    if (strict) {
        for (t=t1; t<=t2; t++) {
            if (x[t] <= 0) {
                return 0;
            }
        }
    } else {
        for (t=t1; t<=t2; t++) {
            if (x[t] < 0) {
                return 0;
            }
        }
    }

    return 1;
}

/**
 * gretl_is_oprobit_ok:
 * @x: data series to examine.
 * @t1: starting observation.
 * @t2: ending observation.
 *
 * Checks the variable @x over the range @t1 to @t2 for its
 * suitability as the dependent variable in an ordered probit
 * analysis.  The criterion used is that the variable has
 * only non-negative integer values.
 *
 * Returns: 1 if the test succeeds, otherwise 0.
 */

int gretl_is_oprobit_ok(int t1, int t2, const(double) * x)
{
    int t, n = 0, d = 1;

    for (t=t1; t<=t2; t++) {
        if (na(x[t])) {
            continue;
        }
        n++;
        if (x[t] != floor(x[t]) || x[t] < 0) {
            d = 0;
            break;
        }
    }

    return (d > 0 && n > 0);
}

/**
 * true_const:
 * @v: index number of variable to test.
 * @dset: dataset struct.
 *
 * Check whether variable Z[v] equals 1 over the sample
 * range given in @dset (aside from any missing values).
 *
 * Returns: 1 if so, otherwise 0.
 */

int true_const(int v, const(DatasetRaw) * dset)
{
    if (v < 0 || v >= dset.v) {
        return 0;
    }

    return gretl_isunits(dset.t1, dset.t2, dset.Z[v]);
}

/**
 * gretl_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * ascending order.
 *
 * Returns: appropriate value for qsort.
 */

extern (C) int gretl_compare_doubles(const(void) * a, const(void) * b)
{
    const(double) * da = cast(const(double) *) a;
    const(double) * db = cast(const(double) *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return -1;
        } else if (!isnan(*db)) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return (*da > *db) - (*da < *db);
    }
}

/**
 * gretl_inverse_compare_doubles:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts doubles in
 * descending order.
 *
 * Returns: appropriate value for qsort.
 */

extern (C) int gretl_inverse_compare_doubles(const(void) * a, const(void) * b)
{
    const(double) * da = cast(const(double) *) a;
    const(double) * db = cast(const(double) *) b;

    if (isnan(*da) || isnan(*db)) {
        if (!isnan(*da)) {
            return 1;
        } else if (!isnan(*db)) {
            return -1;
        } else {
            return 0;
        }
    } else {
        return (*da < *db) - (*da > *db);
    }
}

/**
 * gretl_compare_ints:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts integers in
 * ascending order.
 *
 * Returns: appropriate value for qsort.
 */

extern (C) int gretl_compare_ints(const(void) * a, const(void) * b)
{
    const(int) * ia = cast(const(int) *) a;
    const(int) * ib = cast(const(int) *) b;

    return *ia - *ib;
}

/**
 * gretl_compare_strings:
 * @a: pointer to first element to compare.
 * @b: pointer to second element to compare.
 *
 * Comparison function for use with qsort.  Sorts strings in
 * alphabetical order.
 *
 * Returns: appropriate value for qsort.
 */

extern (C) int gretl_compare_strings(const(void) * a, const(void) * b)
{
    const(char *) * sa = cast(const(char *) *) a;
    const(char *) * sb = cast(const(char *) *) b;

    return g_strcmp0(*sa, *sb);
}

/**
 * count_distinct_values:
 * @x: sorted array of doubles.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_values(const(double) * x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
        if (x[i] != x[i-1]) {
            c++;
        }
    }

    return c;
}

/**
 * count_distinct_int_values:
 * @x: sorted array of ints.
 * @n: number of elements in array.
 *
 * Returns: the number of distinct values in array @x,
 * provided that @x is already sorted.
 */

int count_distinct_int_values(const(int) * x, int n)
{
    int i, c = 1;

    for (i=1; i<n; i++) {
        if (x[i] != x[i-1]) c++;
    }

    return c;
}

/**
 * rearrange_id_array:
 * @x: sorted array of doubles.
 * @m: number of distinct values in array.
 * @n: number of elements in array.
 *
 * Rearranges the sorted array @x such that the first @m
 * elements contain the @m distinct values in sorted order.
 *
 * Returns: 0 on success, 1 on error (in case @m is greater
 * than @n).
 */

int rearrange_id_array(double * x, int m, int n)
{
    int i, k = 1;

    if (m >= n || m == 1) {
        return 1;
    }

    for (i=1; i<n && k<m; i++) {
        if (x[i] != x[i-1]) {
            x[k++] = x[i];
        }
    }

    return 0;
}

/**
 * printlist:
 * @list: array of integers.
 * @msg: message to print along with @list (or NULL).
 *
 * Prints to stderr the given @list of integers along with a message.
 */

void printlist(const(int) * list, const(char) * msg)
{
    int i;

    if (msg) {
        fprintf(stderr, "%s:\n", msg);
    } else {
        fprintf(stderr, "list: ");
    }

    if (list == null) {
        fputs("list is NULL", stderr);
    } else {
        fprintf(stderr, "%d : ", list[0]);
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                fputs("; ", stderr);
            } else {
                fprintf(stderr, "%d ", list[i]);
            }
        }
    }

    fputc('\n', stderr);
}

/**
 * gretl_calculate_criteria:
 * @ess: error sum of squares.
 * @n: number of observations.
 * @k: number of parameters estimated.
 * @ll: pointer to recieve loglikelihood.
 * @aic: pointer to recieve Akaike criterion.
 * @bic: pointer to recieve Schwarz Bayesian criterion.
 * @hqc: pointer to recieve Hannan-Quinn criterion.
 *
 * Calculates model selection criteria based on @ess, @nobs and
 * @k, for a model estimated via least squares.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_calculate_criteria(double ess, int n, int k,
                             double * ll, double * aic, double * bic,
                             double * hqc)
{
    double lnl = 0;
    double[3] c;
    int err = 0;

    if (na(ess) || ess <= 0.0 || n <= k) {
        err = 1;
    } else {
        errno = 0;
        lnl = -.5 * n * log(ess);
        if (errno == EDOM || errno == ERANGE) {
            err = 1;
        } else {
            lnl += -.5 * n * (1 + LN_2_PI - log(n));
            c[0] = -2.0 * lnl + 2 * k;
            c[1] = -2.0 * lnl + k * log(n);
            c[2] = -2.0 * lnl + 2 * k * log(log(n));
        }
    }

    if (err) {
        *ll = *aic = *bic = *hqc = NADBL;
    } else {
        *ll = lnl;
        *aic = c[0];
        *bic = c[1];
        *hqc = c[2];
    }

    return err;
}

/**
 * ls_criteria:
 * @pmod: pointer to gretl model structure.
 *
 * Fills out the model selection criteria members of @pmod, using
 * gretl_calculate_criteria().
 *
 * Returns: 0 on success, non-zero on error.
 */

int ls_criteria(Model * pmod)
{
    double ll, aic, bic, hqc;
    int err;

    err = gretl_calculate_criteria(pmod.ess, pmod.nobs, pmod.ncoeff,
                                   &ll, &aic, &bic, &hqc);

    pmod.lnL = ll;
    pmod.criterion[InfoCriterion.aic] = aic;
    pmod.criterion[InfoCriterion.bic] = bic;
    pmod.criterion[InfoCriterion.hqc] = hqc;

    return err;
}

private char * real_format_obs(char * obs, int maj, int min, int pd, char sep)
{
    int n;

    if (pd >= 10) {
        int pdp = pd / 10, minlen = 2;
        char[18] fmt;

        while ((pdp = pdp / 10) != 0) minlen++;
        sprintf(fmt.ptr, "%%d%c%%0%dd", cast(int) sep, minlen);
        n = snprintf(obs, OBSLEN, fmt.ptr, maj, min);
    } else {
        n = snprintf(obs, OBSLEN, "%d%c%d", maj, cast(int) sep, min);
    }

    if (n < 0 || n >= OBSLEN) {
        /* The formatted observation string does not fit into the
           #OBSLEN-sized target buffer: signal failure rather than
           silently truncating or overflowing it.
        */
        return null;
    }

    return obs;
}

/**
 * format_obs:
 * @obs: target string (should be of length #OBSLEN).
 * @maj: major period (e.g. year).
 * @min: minor period (e.g. quarter, month).
 * @pd: data frequency.
 *
 * Prints to @obs the gretl-type date string representing
 * the observation given by @maj, @min and @pd.
 *
 * Returns: @obs on success or NULL on invalid input.
 */

char * format_obs(char * obs, int maj, int min, int pd)
{
    return real_format_obs(obs, maj, min, pd, ':');
}

private int get_stobs_maj_min(char * stobs, int structure,
                              int * maj, int * min)
{
    int dotc = 0;
    char * p = stobs;
    int err = 0;

    if (stobs == null || *stobs == '\0') {
        /* Empty starting-obs string is never valid: reject it here,
           before any structure- or pd-dependent branch, and before
           any length-based indexing below.
        */
        return 1;
    }

    while (*p) {
        if (*p == ':') {
            *p = '.';
            dotc++;
        } else if (*p == '.') {
            dotc++;
        } else if (!isdigit(cast(ubyte) *p)) {
            err = 1;
            break;
        }
        p++;
    }

    if (!err) {
        if (dotc > 1 || *stobs == '.' ||
            stobs[strlen(stobs) - 1] == '.') {
            err = 1;
        }
    }

    if (!err) {
        if (dotc > 0) {
            sscanf(stobs, "%d.%d", maj, min);
            if (*maj <= 0 || *min <= 0) {
                err = 1;
            }
        } else {
            sscanf(stobs, "%d", maj);
            if (*maj <= 0 && structure != TsCodes.special_time_series) {
                err = 1;
            }
        }
    }

    return err;
}

private int catch_setobs_errors(const(char) * stobs, int pd, int min, int structure)
{
    int panel = structure == TsCodes.stacked_time_series ||
        structure == TsCodes.stacked_cross_section;
    int err = 0;

    if (pd == 1) {
        if (min > 0) {
            gretl_errmsg_set(_("no ':' allowed in starting obs with "
                               ~ "frequency 1").ptr);
            err = 1;
        } else if (panel) {
            gretl_errmsg_set(_("panel data must have frequency > 1").ptr);
            err = 1;
        }
    } else {
        if (min == 0) {
            gretl_errmsg_set(_("starting obs must contain a ':' with "
                               ~ "frequency > 1").ptr);
            err = 1;
        } else if (min > pd) {
            gretl_errmsg_sprintf(_("starting obs '%s' is incompatible with frequency").ptr,
                                 stobs);
            err = 1;
        } else if (structure == TsCodes.cross_section) {
            gretl_errmsg_set(_("cross-sectional data: frequency must be 1").ptr);
            err = 1;
        }
    }

    return err;
}

private int invalid_stobs(const(char) * s)
{
    gretl_errmsg_sprintf(_("starting obs '%s' is invalid").ptr, s);
    return GretlError.data;
}

private void maybe_fix_daily_start(uint * ed, int pd)
{
    int wday = weekday_from_epoch_day(*ed);
    int fix = 0;

    if (wday == GDateWeekday.sunday) {
        /* 5- or 6-day data: sunday not valid */
        fix = 1;
    } else if (pd == 5 && wday == GDateWeekday.saturday) {
        /* 5-day data: saturday not valid */
        fix = 2;
    }

    if (fix) {
        char * fixed;
        char * msg;

        *ed += fix;
        fixed = ymd_extended_from_epoch_day(*ed, 0, null);
        msg = gretl_strdup_printf(_("the starting date was corrected to Monday %s").ptr,
                                  fixed);
        gretl_warnmsg_set(msg);
        free(msg);
        free(fixed);
    }
}

private bool likely_calendar_obs_string(const(char) * s)
{
    return strchr(s, '-') || strchr(s, '/');
}

private bool recognized_ts_frequency(int f)
{
    return f == 4 || f == 12 || f == 24;
}

private bool calendar_pd(int pd)
{
    return pd == 5 || pd == 6 || pd == 7 || pd == 52;
}

private int process_starting_obs(const(char) * stobs_in, int pd,
                                 int * pstructure, double * psd0,
                                 uint * ped0, Gretlopt opt)
{
    char[OBSLEN] stobs;
    int structure = *pstructure;
    double sd0 = 0.0;
    int maybe_tseries = 1;
    int dated = 0;
    int err = 0;

    if (structure == TsCodes.cross_section ||
        structure == TsCodes.stacked_time_series ||
        structure == TsCodes.stacked_cross_section) {
        maybe_tseries = 0;
    }

    stobs[0] = '\0';
    strncat(stobs.ptr, stobs_in, OBSLEN - 1);

    /* check for possible calendar date */
    if (likely_calendar_obs_string(stobs.ptr)) {
        if (maybe_tseries && calendar_pd(pd)) {
            dated = 1;
        } else {
            return invalid_stobs(stobs.ptr);
        }
    }

    if (dated) {
        /* calendar-dated data, daily or weekly */
        uint ed0 = get_epoch_day(stobs.ptr);

        if (ed0 <= 0) {
            return invalid_stobs(stobs.ptr);
        } else {
            if (pd < 7) {
                maybe_fix_daily_start(&ed0, pd);
            }
            sd0 = ed0;
            *ped0 = ed0;
            structure = TsCodes.time_series;
        }
    } else if (structure == TsCodes.time_series && pd == 10) {
        /* decennial data */
        sd0 = cast(double) atoi(stobs.ptr);
    } else {
        int maj = 0, min = 0;

        if (get_stobs_maj_min(stobs.ptr, structure, &maj, &min)) {
            return invalid_stobs(stobs.ptr);
        }

        if (calendar_pd(pd) && min == 0 && maybe_tseries) {
            /* catch undated daily or weekly data */
            structure = TsCodes.time_series;
        } else {
            if (catch_setobs_errors(stobs.ptr, pd, min, structure)) {
                return GretlError.data;
            } else if (pd == 1) {
                sprintf(stobs.ptr, "%d", maj);
                if (structure == TsCodes.structure_unknown) {
                    if (maj > 1) {
                        structure = TsCodes.time_series; /* annual? */
                    } else {
                        structure = TsCodes.cross_section;
                    }
                }
            } else {
                if (structure == TsCodes.time_series && min > 0 &&
                    !recognized_ts_frequency(pd)) {
                    if (opt & Gretlopt.OPT_I) {
                        /* --panel-time */
                        return invalid_stobs(stobs.ptr);
                    } else {
                        structure = TsCodes.special_time_series;
                    }
                }
                if (real_format_obs(stobs.ptr, maj, min, pd, '.') == null) {
                    /* The starting obs, as formatted for this frequency,
                       is too long for a valid gretl observation string.
                    */
                    return invalid_stobs(stobs_in);
                }
                if (structure == TsCodes.structure_unknown &&
                    recognized_ts_frequency(pd)) {
                    structure = TsCodes.time_series;
                }
            }
        }

        /* for non-calendar data */
        sd0 = dot_atof(stobs.ptr);
    }

    if (!err) {
        *pstructure = structure;
        *psd0 = sd0;
    }

    return err;
}

/**
 * set_obs:
 * @parm1: first parameter (periodicity: integer, as string).
 * @parm2: second parameter (starting observation, string).
 * @dset: dataset struct.
 * @opt: %OPT_S for stacked time-series, %OPT_C for stacked cross-section,
 * %OPT_T for time series, %OPT_X for cross section, %OPT_P to set
 * panel structure via two variables representing unit and period
 * respectively. For data already set as panel, %OPT_G to set panel
 * group names or %OPT_I to set panel time-dimension information.
 *
 * Set the frequency and initial observation for a dataset, or
 * in the case of a panel dataset, extra group or time information.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int set_obs(const(char) * parm1, const(char) * parm2,
            DatasetRaw * dset, Gretlopt opt)
{
    const(char) * stobs = null;
    int structure = TsCodes.structure_unknown;
    double sd0 = dset.sd0;
    uint ed0 = 0;
    int pd, panel = 0;
    int err = 0;

    if (dset == null) {
        return GretlError.nodata;
    }

    if ((opt & (Gretlopt.OPT_R | Gretlopt.OPT_P)) && dset.Z == null) {
        return GretlError.nodata;
    }

    if ((opt & (Gretlopt.OPT_G | Gretlopt.OPT_I)) && !dataset_is_panel(dset)) {
        return GretlError.data;
    }

    /* prevent substantive reorganization of the dataset
       within a function */
    if ((opt & (Gretlopt.OPT_P | Gretlopt.OPT_C)) && gretl_function_depth() > 0) {
        gretl_errmsg_set(_("You cannot do this within a function").ptr);
        return GretlError.data;
    }

    gretl_error_clear();

    if (opt & Gretlopt.OPT_R) {
        /* restructure panel: "hidden" option */
        return switch_panel_orientation(dset);
    }

    if (opt == Gretlopt.OPT_NONE && parm1 == null && parm2 == null) {
        /* we'll just print current obs info */
        return 0;
    }

    if (parm1 == null || (!(opt & Gretlopt.OPT_G) && parm2 == null)) {
        /* all cases need at least one param, most need two */
        return GretlError.args;
    }

    if (opt & Gretlopt.OPT_P) {
        return set_panel_structure_from_varnames(parm1, parm2, dset);
    } else if (opt & Gretlopt.OPT_G) {
        /* --panel-groups */
        return set_panel_group_strings(parm1, parm2, dset);
    }

    /* now we get down to business */

    pd = gretl_int_from_string(parm1, &err);
    if (!err && pd < 1) {
        gretl_errmsg_sprintf(_("frequency (%d) does not seem to make sense").ptr, pd);
        return GretlError.data;
    }

    stobs = parm2;

    /* if an explicit structure option was passed in, respect it */
    if (opt == Gretlopt.OPT_X) {
        structure = TsCodes.cross_section;
    } else if (opt == Gretlopt.OPT_T) {
        structure = TsCodes.time_series;
    } else if (opt == Gretlopt.OPT_S) {
        structure = TsCodes.stacked_time_series;
        panel = 1;
    } else if (opt == Gretlopt.OPT_C) {
        structure = TsCodes.stacked_cross_section;
        panel = 1;
    } else if (opt == Gretlopt.OPT_N) {
        structure = TsCodes.special_time_series;
    } else if (opt == Gretlopt.OPT_I) {
        /* --panel-time */
        structure = TsCodes.time_series;
    }

    if (panel && dset.n > 0 && pd > dset.n) {
        gretl_errmsg_sprintf(_("frequency (%d) does not seem to make sense").ptr, pd);
        return 1;
    }

    err = process_starting_obs(stobs, pd, &structure, &sd0, &ed0, opt);

    if (err) {
        return err;
    }

    if (opt == Gretlopt.OPT_I) {
        /* --panel-time */
        dset.panel_pd = pd;
        dset.panel_sd0 = sd0;
        return 0;
    }

    if (panel && dset.n % pd != 0) {
        int sts = structure == TsCodes.stacked_time_series;

        gretl_errmsg_sprintf(_("Panel datasets must be balanced.\n"
                               ~ "The number of observations (%d) is not a multiple\n"
                               ~ "of the number of %s (%d).").ptr,
                             dset.n, sts ? _("periods").ptr : _("units").ptr, pd);
        return GretlError.data;
    }

    if (ed0 > 0) {
        /* replace any existing markers with date strings */
        dataset_destroy_obs_markers(dset);
    } else if (structure == TsCodes.time_series && (pd == 1 || pd == 4 || pd == 12)) {
        /* force use of regular time-series obs labels */
        dataset_destroy_obs_markers(dset);
    }

    dset.pd = pd;
    dset.structure = cast(typeof(dset.structure)) structure;
    dset.sd0 = sd0;

    if (ed0 > 0) {
        calendar_date_string(dset.stobs.ptr, 0, dset);
        calendar_date_string(dset.endobs.ptr, dset.n - 1, dset);
    } else {
        ntolabel(dset.stobs.ptr, 0, dset);
        ntolabel(dset.endobs.ptr, dset.n - 1, dset);
    }

    /* pre-process stacked cross-sectional panels: put into canonical
       stacked time series form
    */
    if (dset.structure == TsCodes.stacked_cross_section) {
        err = switch_panel_orientation(dset);
    }
    if (dataset_is_panel(dset) && dset.panel_pd == 0) {
        dset.panel_pd = 1;
    }

    return err;
}

/**
 * simple_set_obs:
 * @dset: pointer to data information struct.
 * @pd: data frequency.
 * @stobs: string representation of starting observation.
 * @opt: options flags.
 *
 * See the "setobs" command.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int simple_set_obs(DatasetRaw * dset, int pd, const(char) * stobs,
                   Gretlopt opt)
{
    int structure = TsCodes.structure_unknown;
    double sd0 = dset.sd0;
    uint ed0 = 0;
    int panel = 0;
    int err = 0;

    if (opt == Gretlopt.OPT_X) {
        structure = TsCodes.cross_section;
    } else if (opt == Gretlopt.OPT_T) {
        structure = TsCodes.time_series;
    } else if (opt == Gretlopt.OPT_S) {
        structure = TsCodes.stacked_time_series;
        panel = 1;
    } else if (opt == Gretlopt.OPT_C) {
        structure = TsCodes.stacked_cross_section;
        panel = 1;
    } else if (opt == Gretlopt.OPT_N) {
        structure = TsCodes.special_time_series;
    } else if (opt == Gretlopt.OPT_I) {
        /* --panel-time */
        structure = TsCodes.time_series;
    }

    err = process_starting_obs(stobs, pd, &structure, &sd0, &ed0, Gretlopt.OPT_NONE);

    if (err) {
        return err;
    }

    if (panel && dset.n % pd != 0) {
        int sts = structure == TsCodes.stacked_time_series;

        gretl_errmsg_sprintf(_("Panel datasets must be balanced.\n"
                               ~ "The number of observations (%d) is not a multiple\n"
                               ~ "of the number of %s (%d).").ptr,
                             dset.n, sts ? _("periods").ptr : _("units").ptr, pd);
        return GretlError.data;
    }

    if (ed0 > 0) {
        /* replace any existing markers with date strings */
        dataset_destroy_obs_markers(dset);
    } else if (structure == TsCodes.time_series && (pd == 1 || pd == 4 || pd == 12)) {
        /* force use of regular time-series obs labels */
        dataset_destroy_obs_markers(dset);
    }

    dset.pd = pd;
    dset.structure = cast(typeof(dset.structure)) structure;
    dset.sd0 = sd0;

    if (ed0 > 0) {
        calendar_date_string(dset.stobs.ptr, 0, dset);
        calendar_date_string(dset.endobs.ptr, dset.n - 1, dset);
    } else {
        ntolabel(dset.stobs.ptr, 0, dset);
        ntolabel(dset.endobs.ptr, dset.n - 1, dset);
    }

    if (dset.structure == TsCodes.stacked_cross_section) {
        err = switch_panel_orientation(dset);
    }

    return err;
}

/**
 * gretl_double_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of a double,
 * return its value, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * otherwise set the content of @err to a non-zero value.
 *
 * Returns: double value.
 */

double gretl_double_from_string(const(char) * s, int * err)
{
    const(char) * test;
    double x;

    if (s == null || *s == '\0') {
        *err = GretlError.data;
        return NADBL;
    }

    if (isalpha(*s)) {
        return get_scalar_value_by_name(s, err);
    }

    gretl_push_c_numeric_locale();
    errno = 0;
    x = strtod(s, &test);
    gretl_pop_c_numeric_locale();

    if (*test != '\0' || errno == ERANGE) {
        *err = GretlError.data;
        errno = 0;
        return NADBL;
    }

    return x;
}

/**
 * gretl_int_from_string:
 * @s: string to examine.
 * @err: location to receive error code.
 *
 * If @s is a valid string representation of an integer,
 * return that integer, otherwise if @s is the name of a
 * scalar variable, return the value of that variable,
 * provided it can be converted to an integer, otherwise
 * set the content of @err to a non-zero value.
 *
 * Returns: integer value.
 */

int gretl_int_from_string(const(char) * s, int * err)
{
    const(char) * test;
    int n = 0;

    if (s == null || *s == '\0') {
        *err = GretlError.data;
        return 0;
    }

    if (isalpha(*s)) {
        double x = get_scalar_value_by_name(s, err);

        if (!*err) {
            n = gretl_int_from_double(x, err);
        }
        return n;
    }

    errno = 0;
    n = cast(int) strtol(s, &test, 10);

    if (*test != '\0' || errno == ERANGE) {
        *err = GretlError.data;
        errno = 0;
        return 0;
    }

    return n;
}

/**
 * positive_int_from_string:
 * @s: string to examine.
 *
 * If @s is a valid string representation of a positive integer
 * that can be represented as a 32-bit signed value return that
 * integer, otherwise return -1.
 *
 * Returns: integer value.
 */

int positive_int_from_string(const(char) * s)
{
    int ret = -1;

    if (s != null && *s != '\0') {
        const(char) * test;
        long il;

        errno = 0;

        il = strtol(s, &test, 10);
        if (*test != '\0' || strcmp(s, test) == 0 ||
            errno == ERANGE || il <= 0 || il > int.max) {
            ret = -1;
        } else {
            ret = cast(int) il;
        }
    }

    return ret;
}

private int letter_to_int(char c)
{
    const(char) * s = "abcdefghij";
    int i = 0;

    while (*s) {
        if (c == *s) {
            return i;
        }
        s++;
        i++;
    }

    return 0;
}

/**
 * libgretl_version:
 *
 * Returns: the integer version of libgretl.
 */

int libgretl_version()
{
    return LIBGRETL_CURRENT;
}

/**
 * gretl_version_number:
 * @version: gretl program version in string form.
 *
 * Returns: the integer gretl version number.
 */

int gretl_version_number(const(char) * version_)
{
    int vnum = 0;

    if (!strcmp(version_, "@VERSION@")) {
        version_ = GRETL_VERSION.ptr;
    }

    if (atoi(version_) >= 2015) {
        /* as in "2015d" and subsequent releases */
        int Y;
        char c;

        sscanf(version_, "%d%c", &Y, &c);
        vnum = 10 * Y + letter_to_int(c);
    } else {
        /* old style: "1.9.14" or whatever */
        int x, y, z;

        sscanf(version_, "%d.%d.%d", &x, &y, &z);
        vnum = 10000 * x + 100 * y + z;
    }

    return vnum;
}

/* Mapping from new-style version numbers, based on
   year, to old style based on maj.min.pl. This
   covers gretl 1.9.4 to 1.10.2.
*/

private int old_style_gretl_version_number(int v)
{
    static immutable int[2][18] vtrans = [
        [10904, 20110],
        [10905, 20111],
        [10906, 20112],
        [10907, 20113],
        [10908, 20120],
        [10909, 20121],
        [10910, 20122],
        [10911, 20123],
        [10912, 20130],
        [10913, 20131],
        [10914, 20132],
        [10990, 20140],
        [10991, 20141],
        [10992, 20142],
        [11000, 20150],
        [11001, 20151],
        [11002, 20152],
        [11003, 20153]
    ];
    int i;

    for (i=0; i<17; i++) {
        if (v == vtrans[i][1] || v < vtrans[i+1][1]) {
            return vtrans[i][0];
        }
    }

    /* default to 1.9.4 */

    return vtrans[0][0];
}

/**
 * gretl_version_string:
 * @targ: string into which to write (9 bytes minimum).
 * @vnum: integer program version number.
 *
 * Returns: the string representation of @vnum, in @targ.
 */

char * gretl_version_string(char * targ, int vnum)
{
    if (vnum >= 20153) {
        /* "2105d" or higher */
        const(char) * s = "abcdefghij";
        int y, idx;
        char c;

        y = vnum / 10;
        idx = vnum - 10 * y;

        if (idx >= 0 && idx < 10) {
            c = s[idx];
        } else {
            c = 'a';
        }

        sprintf(targ, "%d%c", y, cast(int) c);
    } else {
        int x, y, z;

        if (vnum > 20000) {
            /* translate back to old-style */
            vnum = old_style_gretl_version_number(vnum);
        }

        x = vnum / 10000;
        y = (vnum - x * 10000) / 100;
        z = vnum % 100;

        sprintf(targ, "%d.%d.%d", x, y, z);
    }

    return targ;
}

/**
 * varnum_from_string:
 * @str: string representation of an integer ID number.
 * @dset: dataset information.
 *
 * Returns: integer ID number, or -1 on failure.
 */

int varnum_from_string(const(char) * str, DatasetRaw * dset)
{
    int v = positive_int_from_string(str);

    if (v <= 0 || v >= dset.v) {
        v = -1;
    }

    return v;
}

/**
 * gretl_int_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an integer, if
 * possible. Otherwise returns -1 with @err set to a
 * non-zero value. Note that it is considered an
 * error if @x is "too far" from the nearest integer;
 * it must be "almost integral", with tolerance 0.001.
 */

int gretl_int_from_double(double x, int * err)
{
    int k = -1;

    if (na(x) || fabs(x) > int.max || fabs(x - nearbyint(x)) > 0.001) {
        *err = GretlError.invarg;
    } else {
        k = cast(int) lrint(x);
    }

    return k;
}

/**
 * gretl_uint32_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an unsigned
 * 32-bit integer, if possible. Otherwise returns -1
 * with @err set to a non-zero value. Note that it is
 * considered an error if @x is "too far" from the
 * nearest integer; it must be "almost integral", with
 * tolerance 1.0e-6.
 */

uint gretl_uint32_from_double(double x, int * err)
{
    uint u = 0;

    if (na(x) || x < 0 || fabs(x) > uint.max) {
        *err = GretlError.invarg;
    } else {
        double f = floor(x);
        double c = ceil(x);

        if (x - f < 1e-6) {
            u = cast(uint) f;
        } else if (c - x < 1e-6) {
            u = cast(uint) c;
        } else {
            *err = GretlError.invarg;
        }
    }

    return u;
}

/**
 * gretl_uint64_from_double:
 * @x: double-precision floating point value
 * @err: location to receive error code.
 *
 * Returns: the value of @x converted to an unsigned 64-bit
 * integer, if possible. Otherwise returns -1 with @err set
 * to a non-zero value. It is considered an error if @x is
 * not "almost integral", with tolerance 1.0e-6.
 */

ulong gretl_uint64_from_double(double x, int * err)
{
    ulong k = 0;

    if (na(x) || x > 9007199254740991 /* 2^53 - 1 */) {
        *err = GretlError.invarg;
    } else {
        double f = floor(x);
        double c = ceil(x);

        if (x - f < 1e-6) {
            k = cast(ulong) f;
        } else if (c - x < 1e-6) {
            k = cast(ulong) c;
        } else {
            *err = GretlError.invarg;
        }
    }

    return k;
}

/**
 * gretl_type_from_name:
 * @s: the name to check.
 * @dset: dataset information.
 *
 * Returns: non-zero if @s is the name of an existing
 * series, matrix, scalar, list, string or bundle,
 * otherwise 0.
 */

GretlType gretl_type_from_name(const(char) * s, const(DatasetRaw) * dset)
{
    if (gretl_is_series(s, dset)) {
        return GretlType.series;
    } else {
        return user_var_get_type_by_name(s);
    }
}

/**
 * copyvec:
 * @src: array of doubles.
 * @n: number of elements to copy.
 *
 * Returns: an allocated copy of the first @n elements of
 * array @src, or %NULL on failure.
 */

double * copyvec(const(double) * src, int n)
{
    double * targ = null;
    size_t sz = n * (*targ).sizeof;

    if (n > 0 && src != null) {
        targ = cast(double *) malloc(sz);
    }

    if (targ != null) {
        memcpy(targ, src, sz);
    }

    return targ;
}

/**
 * doubles_array_free:
 * @X: 2-dimensional array of doubles.
 * @m: number of sub-arrays.
 *
 * Frees a 2-dimensional array of doubles, first freeing
 * each sub-array.
 */

void doubles_array_free(double ** X, int m)
{
    if (X != null) {
        int i;

        for (i=0; i<m; i++) {
            free(X[i]);
        }
        free(X);
    }
}

/**
 * doubles_array_new:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Allocates a 2-dimensional array of doubles, that is,
 * @m arrays each containing @n elements.  If @n is
 * zero the sub-arrays are just set to %NULL.
 *
 * Returns: the allocated array, or %NULL on failure.
 */

double ** doubles_array_new(int m, int n)
{
    double ** X;
    int i;

    if (m == 0) {
        return null;
    }

    X = cast(double **) malloc(m * (*X).sizeof);
    if (X == null) {
        return X;
    }

    for (i=0; i<m; i++) {
        X[i] = null;
    }

    if (n > 0) {
        for (i=0; i<m; i++) {
            X[i] = cast(double *) malloc(n * (**X).sizeof);
            if (X[i] == null) {
                doubles_array_free(X, m);
                X = null;
                break;
            }
        }
    }

    return X;
}

/**
 * doubles_array_new0:
 * @m: number of sub-arrays.
 * @n: length of each sub-array.
 *
 * Works just as doubles_array_new(), except that on
 * successful allocation all values in the arrays are
 * set to zero.
 *
 * Returns: the allocated array, or %NULL on failure.
 */

double ** doubles_array_new0(int m, int n)
{
    double ** X = doubles_array_new(m, n);

    if (X != null && n > 0) {
        int i, j;

        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                X[i][j] = 0.0;
            }
        }
    }

    return X;
}

/**
 * doubles_array_adjust_length:
 * @X: original two-dimensional array.
 * @m: number of sub-arrays in @X.
 * @new_n: new length of each sub-array.
 *
 * For each of the @m sub-arrays in @X, reallocate to
 * a length of @new_n.
 *
 * Returns: 0 on success, non-zero on error.
 */

int doubles_array_adjust_length(double ** X, int m, int new_n)
{
    int err = 0;

    if (X == null || m == 0) {
        ; /* no-op */
    } else {
        double * xi;
        int i;

        for (i=0; i<m && !err; i++) {
            if (new_n == 0) {
                free(X[i]);
                X[i] = null;
            } else {
                xi = cast(double *) realloc(X[i], new_n * (*xi).sizeof);
                if (xi == null) {
                    err = GretlError.alloc;
                } else {
                    X[i] = xi;
                }
            }
        }
    }

    return err;
}

/**
 * data_array_from_model:
 * @pmod: reference model.
 * @Z: main data array.
 * @missv: should equal 1 if there are missing values to be
 * skipped, else 0.
 *
 * Constructs a dataset containing all the variables referenced in
 * @pmod.  The arrays start at the correct sample offset for @pmod,
 * and are contiguous.  If @missvals equals 0, this is done by creating
 * a set of pointers into the main dataset, but if there are missing
 * values to be handled, the sub-arrays are newly allocated and purged
 * of NAs.
 *
 * Returns: two-dimensional array, or %NULL on failure.
 */

double ** data_array_from_model(const(Model) * pmod, double ** Z, int missv)
{
    double ** X;
    int nv = pmod.list[0];
    int offset = pmod.t1;
    int v, i;

    if (missv) {
        X = doubles_array_new(nv, pmod.nobs);
    } else {
        X = cast(double **) malloc(nv * (*X).sizeof);
    }

    if (X == null) return null;

    if (missv) {
        int t, s;

        for (t=0; t<pmod.nobs; t++) {
            X[0][t] = 1.0;
        }

        for (i=1; i<nv; i++) {
            v = (i == 1)? pmod.list[1] : pmod.list[i + 1];
            s = 0;
            for (t=pmod.t1; t<=pmod.t2; t++) {
                if (!na(pmod.uhat[t])) {
                    X[i][s++] = Z[v][t];
                }
            }
        }
    } else {
        /* constant in slot 0 */
        X[0] = Z[0] + offset;

        /* dependent var in slot 1 */
        X[1] = Z[pmod.list[1]] + offset;

        /* independent vars in slots 2, 3, ... */
        for (i=2; i<nv; i++) {
            v = pmod.list[i + 1];
            X[i] = Z[v] + offset;
        }
    }

    return X;
}

/* (the C code's gretl_spawn() etc. use glib's g_spawn_*; see gspawn.d) */

int gretl_spawn(char * cmdline)
{
    GError * error = null;
    char * errout = null;
    char * sout = null;
    int ok, status;
    int ret = 0;

    gretl_error_clear();

    ok = g_spawn_command_line_sync(cmdline,
                                   &sout,   /* standard output */
                                   &errout, /* standard error */
                                   &status, /* exit status */
                                   &error);

    if (!ok) {
        gretl_errmsg_set(error.message);
        fprintf(stderr, "gretl_spawn: '%s'\n", error.message);
        g_error_free(error);
        if (errout != null) {
            fprintf(stderr, " stderr = '%s'\n", errout);
        }
        ret = 1;
    } else if (status != 0) {
        if (errout != null && *errout) {
            gretl_errmsg_set(errout);
            fprintf(stderr, "gretl_spawn: status = %d: '%s'\n", status, errout);
        } else if (sout != null && *sout) {
            gretl_errmsg_set(sout);
            fprintf(stderr, "gretl_spawn: status = %d: '%s'\n", status, sout);
        } else {
            gretl_errmsg_set(_("Command failed").ptr);
            fprintf(stderr, "gretl_spawn: status = %d\n", status);
        }
        ret = 1;
    }

    if (errout != null) g_free(errout);
    if (sout != null) g_free(sout);

    if (ret) {
        fprintf(stderr, "Failed command: '%s'\n", cmdline);
    }

    return ret;
}

int gnuplot_make_image(const(char) * input_fname)
{
    char * cmdline;
    int err;

    cmdline = g_strdup_printf("\"%s\" \"%s\"",
                              gretl_gnuplot_path(),
                              input_fname);
    err = gretl_spawn(cmdline);
    g_free(cmdline);

    return err;
}

int gretl_pipe_output(char ** argv, char ** envp,
                      const(char) * currdir, Prn * prn,
                      char ** errp)
{
    GError * error = null;
    char * errout = null;
    char * sout = null;
    int ok, status;
    int err = 0;

    gretl_error_clear();

    ok = g_spawn_sync(currdir,
                      argv,
                      envp,    /* may be NULL to inherit env */
                      G_SPAWN_SEARCH_PATH, /* ? */
                      null,    /* child_setup */
                      null,    /* user_data */
                      &sout,   /* standard output */
                      &errout, /* standard error */
                      &status, /* exit status */
                      &error);

    if (!ok) {
        gretl_errmsg_set(error.message);
        fprintf(stderr, "gretl_pipe_output: '%s'\n", error.message);
        g_error_free(error);
        if (errout != null) {
            fprintf(stderr, " stderr = '%s'\n", errout);
        }
        err = 1;
    } else if (status != 0) {
        if (errout != null && *errout) {
            gretl_errmsg_set(errout);
            fprintf(stderr, "gretl_pipe_output: status = %d: '%s'\n", status, errout);
        } else {
            gretl_errmsg_set(_("Command failed").ptr);
            fprintf(stderr, "gretl_pipe_output: status = %d\n", status);
        }
        err = 1;
    }

    if (sout != null) {
        if (*sout != '\0') {
            pputs(prn, sout);
        }
        g_free(sout);
    }
    if (errout != null) {
        if (*errout != '\0') {
            if (errp != null) {
                *errp = errout;
                errout = null;
            } else {
                fputs(errout, stderr);
            }
        }
        g_free(errout);
    }

    return err;
}

/* file copying */

int gretl_copy_file(const(char) * src, const(char) * dest)
{
    FILE * srcfd;
    FILE * destfd;
    char[8192] buf;
    size_t n;

    if (!strcmp(src, dest)) {
        /* @src and @dest are the same: no-op */
        return 0;
    }

    if ((srcfd = gretl_fopen(src, "rb")) == null) {
        gretl_errmsg_sprintf(_("Couldn't open %s").ptr, src);
        return GretlError.fopen;
    }

    if ((destfd = gretl_fopen(dest, "wb")) == null) {
        gretl_errmsg_sprintf(_("Couldn't write to %s").ptr, dest);
        fclose(srcfd);
        return GretlError.fopen;
    }

    while ((n = fread(buf.ptr, 1, buf.sizeof, srcfd)) > 0) {
        fwrite(buf.ptr, 1, n, destfd);
    }

    fclose(srcfd);
    fclose(destfd);

    return 0;
}

private int maybe_unload_function_package(const(char) * s,
                                          Prn * prn)
{
    char[FN_NAMELEN+4] pkgname;
    char * p;
    int done = 0;

    pkgname[0] = '\0';
    strncat(pkgname.ptr, s, FN_NAMELEN+3);
    p = strrchr(pkgname.ptr, '.');

    /* Note: @p will be NULL if @s was long enough that truncation
       to FN_NAMELEN+3 bytes cut off the ".gfn" suffix; in that
       case there's no valid package name to look up, so just
       fall through to "not loaded" below.
    */
    if (p != null) {
        void * pkg = null;
        const(char) * path;

        *p = '\0';
        pkg = get_function_package_by_name(pkgname.ptr);
        if (pkg != null) {
            path = get_function_package_path_by_name(pkgname.ptr);
            if (path != null) {
                function_package_unload_full_by_filename(path);
                done = 1;
            }
        }
    }

    if (done) {
        pprintf(prn, _("Unloaded package %s\n").ptr, pkgname.ptr);
    } else {
        pprintf(prn, _("Package %s was not loaded\n").ptr, pkgname.ptr);
    }

    return 0;
}

/* note: this is not used for series */

int gretl_delete_var_by_name(const(char) * s, Prn * prn)
{
    int err = 0;

    if (s == null || *s == '\0') {
        return GretlError.parse;
    }

    if (object_is_function_arg(s)) {
        gretl_errmsg_sprintf(_("The variable %s is read-only").ptr, s);
        return GretlError.data;
    }

    if (has_suffix(s, ".gfn")) {
        err = maybe_unload_function_package(s, prn);
    } else if (gretl_is_user_var(s)) {
        err = user_var_delete_by_name(s, prn);
    } else {
        /* try for a bundle member? */
        char * genstr = g_strdup_printf("%s=null", s);

        err = generate(genstr, null, GretlType.any, Gretlopt.OPT_P, prn);
        g_free(genstr);
    }

    return err;
}

/* apparatus to support gretl's stopwatch */

/* (The C code has two timer variants: "xtimers", based on OpenMP's
   omp_get_wtime() or MPI's timer, and "itimers", based on
   g_get_monotonic_time(). This port drops OpenMP and MPI, so only the
   itimer variant is kept, with a D dynamic array in place of the
   GPtrArray.) */

private struct ITimer {
    int level;
    gint64 t0;
}

private __gshared ITimer[] itimers;

private ITimer * new_itimer(int level)
{
    itimers ~= ITimer(level, g_get_monotonic_time());
    return &itimers[$ - 1];
}

private void itimer_init()
{
    if (itimers.length == 0) {
        new_itimer(gretl_function_depth());
    }
}

private ITimer * get_itimer()
{
    int i, lev = gretl_function_depth();

    for (i=0; i<itimers.length; i++) {
        if (itimers[i].level == lev) {
            return &itimers[i];
        }
    }

    return new_itimer(lev);
}

private double itimer_stopwatch()
{
    ITimer * it = get_itimer();
    gint64 t1 = g_get_monotonic_time();
    double ret = (t1 - it.t0) * 1.0e-6;

    it.t0 = t1;

    return ret;
}

private void gretl_stopwatch_cleanup()
{
    itimers = null;
}

double gretl_stopwatch()
{
    return itimer_stopwatch();
}

private void gretl_stopwatch_init()
{
    itimer_init();
}

/* BLAS detection and analysis */

private __gshared char[64] blas_core = 0;
private __gshared char[32] blas_parallel = 0;
private __gshared char[32] blas_version = 0;

private __gshared int blas_variant;

/* @ldd_line should look something like this:

   libopenblaso.so.0 => /usr/lib64/libopenblaso.so.0 ...
*/

private void retry_for_blas_details(char * ldd_line,
                                    int blastype)
{
    char * s = strchr(ldd_line, '/');
    char * libpath = null;

    if (s != null) {
        char * p;

        s += strspn(s, " ");
        p = strchr(s, ' ');
        if (p != null) {
            /* copy the library path */
            libpath = gretl_strndup(s, p-s);
        }
    }

    if (libpath != null) {
        void * ptr = dlopen(libpath, RTLD_NOW);

        if (ptr != null) {
            if (blastype == BlasVariant.openblas) {
                register_openblas_details(ptr);
            } else if (blastype == BlasVariant.mkl) {
                register_mkl_details(ptr);
            } else if (blastype == BlasVariant.blis) {
                register_blis_details(ptr);
            }
        }
        free(libpath);
    } else {
        fprintf(stderr, "retry_for_blas_details: couldn't get @libpath\n");
    }
}

private int parse_ldd_output(const(char) * s)
{
    char[512] line;
    int i = 0;
    int ret = BlasVariant.unknown;

    line[0] = '\0';

    while (*s && ret == BlasVariant.unknown) {
        if (*s == '\n') {
            /* got to the end of a line: check it */
            line[i] = '\0';
            if (strstr(line.ptr, "libopenblas")) {
                ret = BlasVariant.openblas;
            } else if (strstr(line.ptr, "libblis")) {
                ret = BlasVariant.blis;
            } else if (strstr(line.ptr, "Accelerate.frame")) {
                ret = BlasVariant.veclib;
            } else if (strstr(line.ptr, "libmkl")) {
                ret = BlasVariant.mkl;
            } else if (strstr(line.ptr, "libatlas")) {
                ret = BlasVariant.atlas;
            } else if (strstr(line.ptr, "libblas")) {
                ret = BlasVariant.netlib;
            } else {
                /* prepare for next line */
                line[0] = '\0';
                i = 0;
            }
        } else {
            /* cumulate this line */
            line[i++] = *s;
        }
        s++;
    }

    if (ret == BlasVariant.openblas || ret == BlasVariant.blis || ret == BlasVariant.mkl) {
        retry_for_blas_details(line.ptr, ret);
    } else if (ret == BlasVariant.unknown) {
        fputs("detect blas: found no relevant libs!\n", stderr);
    }

    return ret;
}

private int detect_blas_via_ldd()
{
    const(char) * targ;
    char * sout = null;
    char * errout = null;
    int status = 0;
    GError * gerr = null;
    int variant = BlasVariant.unknown;
    char *[3] argv;

    version (Windows) {
        return BlasVariant.unknown;    /* no ldd */
    }

    targ = getenv("GRETL_LIB_NAME");
    if (targ == null || *targ == '\0') {
        targ = (GRETL_PREFIX ~ "/lib/libgretl-1.0.so").ptr;
    }
    argv[0] = cast(char *) "ldd".ptr;
    argv[1] = cast(char *) targ;
    argv[2] = null;

    g_spawn_sync(null, argv.ptr, null, G_SPAWN_SEARCH_PATH,
                 null, null, &sout, &errout,
                 &status, &gerr);

    if (gerr != null) {
        fprintf(stderr, "%s\n", gerr.message);
        g_error_free(gerr);
    } else if (status != 0) {
        fprintf(stderr, "%s exited with status %d\n", argv[0], status);
    } else if (sout != null) {
        variant = parse_ldd_output(sout);
    } else {
        fprintf(stderr, "%s: %s\n", argv[0], "Got no output".ptr);
    }

    g_free(sout);
    g_free(errout);

    return variant;
}

private extern (C) {
    alias ObSetNumThreadsFn = void function(int);
    alias ObGetNumThreadsFn = int function();
    alias BlisSetNumThreadsFn = void function(int);
    alias BlisGetNumThreadsFn = int function();
    alias VoidFn = void function();
    alias MklDomainSetNumThreadsFn = void function(int, int);
    alias MklDomainGetMaxThreadsFn = int function(int);
    alias IntFn = int function();
}

private __gshared ObSetNumThreadsFn OB_set_num_threads;
private __gshared ObGetNumThreadsFn OB_get_num_threads;
private __gshared BlisSetNumThreadsFn BLIS_set_num_threads;
private __gshared BlisGetNumThreadsFn BLIS_get_num_threads;
private __gshared VoidFn BLIS_init;
private __gshared VoidFn BLIS_finalize;
private __gshared VoidFn MKL_finalize;
private __gshared MklDomainSetNumThreadsFn MKL_domain_set_num_threads;
private __gshared MklDomainGetMaxThreadsFn MKL_domain_get_max_threads;
private __gshared VoidFn FLAME_init;
private __gshared IntFn FLAME_initialized;
private __gshared VoidFn FLAME_finalize;

private void register_openblas_details(void * handle)
{
    alias GetCorenameFn = extern (C) char * function();
    alias GetConfigFn = extern (C) char * function();

    GetCorenameFn OB_get_corename;
    IntFn OB_get_parallel;
    GetConfigFn OB_get_config;

    OB_get_corename = cast(GetCorenameFn) dlsym(handle, "openblas_get_corename");
    OB_get_parallel = cast(IntFn) dlsym(handle, "openblas_get_parallel");
    OB_get_config = cast(GetConfigFn) dlsym(handle, "openblas_get_config");

    if (OB_get_corename != null) {
        char * s = OB_get_corename();

        if (s != null) {
            blas_core[0] = '\0';
            strncat(blas_core.ptr, s, 31);
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_corename()\n");
    }

    if (OB_get_parallel != null) {
        int p = OB_get_parallel();

        if (p == 0) {
            strcpy(blas_parallel.ptr, "none");
        } else if (p == 1) {
            strcpy(blas_parallel.ptr, "pthreads");
        } else if (p == 2) {
            strcpy(blas_parallel.ptr, "OpenMP");
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_parallel()\n");
    }

    if (OB_get_config != null) {
        char * ver = OB_get_config();

        if (ver != null) {
            blas_version[0] = '\0';
            sscanf(ver, "OpenBLAS %16[^ ]", blas_version.ptr);
        }
    } else {
        fprintf(stderr, "Couldn't find openblas_get_config()\n");
    }
}

private void register_blis_details(void * handle)
{
    /* The last element on arch_t enum in libblis =>
       => must be updated whenever new architecture/cpu model appears*/
    /* Shouldn't this come from a header? (Allin) */
    /* -> Well, this enum is defined in blis.h which we do not include.
       And we can't retrieve it via dlopen(), can we? (Marcin) */
    const int BLIS_NUM_ARCHS = 26;
    const(char) * buf = null;
    int id;

    alias GetVersionStrFn = extern (C) char * function();
    alias ArchStringFn = extern (C) char * function(int);

    /* Functions from libblis we need. */
    GetVersionStrFn BLIS_get_version_str;
    IntFn BLIS_get_enable_threading;
    IntFn BLIS_get_enable_openmp;
    IntFn BLIS_get_enable_pthreads;
    ArchStringFn BLIS_arch_string;
    IntFn BLIS_arch_query_id;

    BLIS_get_version_str = cast(GetVersionStrFn) dlsym(handle, "bli_info_get_version_str");
    BLIS_get_enable_threading = cast(IntFn) dlsym(handle, "bli_info_get_enable_threading");
    BLIS_get_enable_openmp = cast(IntFn) dlsym(handle, "bli_info_get_enable_openmp");
    BLIS_get_enable_pthreads = cast(IntFn) dlsym(handle, "bli_info_get_enable_pthreads");
    BLIS_arch_string = cast(ArchStringFn) dlsym(handle, "bli_arch_string");
    BLIS_arch_query_id = cast(IntFn) dlsym(handle, "bli_arch_query_id");

    if (BLIS_get_version_str != null) {
        buf = BLIS_get_version_str();
        if (buf != null) {
            blas_version[0] = '\0';
            strncat(blas_version.ptr, buf, 31);
            buf = null;
        }
    } else {
        fprintf(stderr, "Couldn't find bli_info_get_version_str()\n");
    }

    /* Model we have: threaded or sequential */
    if (!BLIS_get_enable_threading()) {
        buf = "sequential (non-threading)";
    } else {
        if (BLIS_get_enable_openmp()) {
            buf = "OpenMP";
        } else if (BLIS_get_enable_pthreads()) {
            buf = "pthreads";
        } else {
            buf = "unrecognized";
        }
    }
    if (buf != null) {
        blas_parallel[0] = '\0';
        strncat(blas_parallel.ptr, buf, 31);
        buf = null;
    }

    /* BLIS core in use */
    if (BLIS_arch_query_id != null) {
        id = BLIS_arch_query_id();
    } else {
        fprintf(stderr, "Couldn't find bli_arch_query_id()\n");
        id = BLIS_NUM_ARCHS-1;
    }
    if (BLIS_arch_string != null) {
        buf = BLIS_arch_string(id);
    } else {
        fprintf(stderr, "Couldn't find bli_arch_string()\n");
        buf = "unrecognized";
    }
    if (buf != null) {
        blas_core[0] = '\0';
        strncat(blas_core.ptr, buf, 31);
    }
}

private void register_mkl_details(void * handle)
{
    static struct MKL_version {
        int MajorVersion;
        int MinorVersion;
        int UpdateVersion;
        char * ProductStatus;
        char * Build;
        char * Processor;
        char * Platform;
    }

    alias GetVersionFn = extern (C) void function(MKL_version *);

    MKL_version * version_;
    const(char) * buf = null;
    int nt = 0, blas_nt = 0, id = 0;

    /* Functions from libmkl_intel_lp64.so we need. */
    GetVersionFn MKL_get_version;
    IntFn MKL_get_max_threads;
    IntFn MKL_cbwr_get_auto_branch;

    MKL_get_version = cast(GetVersionFn) dlsym(handle, "MKL_Get_Version");
    MKL_get_max_threads = cast(IntFn) dlsym(handle, "MKL_Get_Max_Threads");
    MKL_cbwr_get_auto_branch = cast(IntFn) dlsym(handle, "MKL_CBWR_Get_Auto_Branch");

    /* MKL version */
    if (MKL_get_version != null) {
        version_ = cast(MKL_version *) malloc((*version_).sizeof);
        if (version_ == null) {
            return;
        } else {
            MKL_get_version(version_);
            if (version_ != null) {
                blas_version[0] = '\0';
                snprintf(blas_version.ptr, 31, "%d.%d build %s",
                         version_.MajorVersion,
                         version_.UpdateVersion,
                         version_.Build);
            }
        }
        free(version_);
    } else {
        fprintf(stderr, "Couldn't find MKL_Get_Version()\n");
    }

    /* Model we have: threaded or sequential */
    if (MKL_get_max_threads != null && MKL_domain_get_max_threads != null) {
        nt = MKL_get_max_threads();
        blas_nt = MKL_domain_get_max_threads(1); /* MKL_DOMAIN_BLAS */
        if (nt == 1)  {
            buf = "sequential (non-threading)";
        } else {
            if (blas_nt == 1) {
                buf = "TBB";
            } else if (blas_nt > 1) {
                buf = "OpenMP";
            } else {
                buf = "unrecognized";
            }
        }
        strncat(blas_parallel.ptr, buf, 31);
        buf = null;
    } else {
        fprintf(stderr, "Couldn't find MKL_get_max_threads() and/or MKL_domain_get_max_threads()\n");
    }

    /* MKL core in use */
    if (MKL_cbwr_get_auto_branch != null) {
        id = MKL_cbwr_get_auto_branch();
        switch (id) {
        case 3: /* MKL_CBWR_COMPATIBLE */
            buf = "SSE2 without rcpps/rsqrtps instructions";
            break;
        case 4: /* MKL_CBWR_SSE2 */
            buf = "SSE2";
            break;
        case 5: /* MKL_CBWR_SSE3 */
            buf = "SSE2 (deprecated SSE3)";
            break;
        case 6: /* MKL_CBWR_SSSE3 */
            buf = "SSSE3";
            break;
        case 7: /* MKL_CBWR_SSE4_1 */
            buf = "SSE4-1";
            break;
        case 8: /* MKL_CBWR_SSE4_2 */
            buf = "SSE4-2";
            break;
        case 9: /* MKL_CBWR_AVX */
            buf = "AVX";
            break;
        case 10: /* MKL_CBWR_AVX2 */
            buf = "AVX2";
            break;
        case 11: /* MKL_CBWR_AVX512_MIC */
            buf = "AVX2 (deprecated AVX-512MIC)";
            break;
        case 12: /* MKL_CBWR_AVX512 */
            buf = "AVX-512";
            break;
        case 13: /* MKL_CBWR_AVX512_MIC_E1 */
            buf = "AVX2 (deprecated AVX-512MICE1)";
            break;
        case 14: /* MKL_CBWR_AVX512_E1 */
            buf = "AVX-512 with support of Vector Neural Network Instructions";
            break;
        case 15: /* MKL_CBWR_AVX10 */
            buf = "AVX-10.2 with support for Advanced Matrix Extensions";
            break;
        default:
            buf = "unknown";
        }
        strncat(blas_core.ptr, buf, 63);
    }
}

private void libflame_utils(void * handle)
{
    /* Functions form libflame we need */
    FLAME_initialized = cast(IntFn) dlsym(handle, "FLA_Initialized");
    FLAME_finalize = cast(VoidFn) dlsym(handle, "FLA_Finalize");
}

/* below: called in creating $sysinfo bundle */

int get_blas_details(char ** s1, char ** s2, char ** s3)
{
    if (blas_core[0] == '\0' || blas_parallel[0] == '\0') {
        return 0;
    } else {
        *s1 = blas_core.ptr;
        *s2 = blas_parallel.ptr;
        if (s3 != null && blas_version[0] != '\0') {
            *s3 = blas_version.ptr;
        }
        return 1;
    }
}

/* for $sysinfo bundle */

int get_stack_size()
{
    rlimit rl;

    if (getrlimit(RLIMIT_STACK, &rl) == 0) {
        return cast(int) rl.rlim_cur;
    }

    return 0;
}

const(char) * blas_variant_string()
{
    if (blas_variant == BlasVariant.netlib) {
        return "netlib";
    } else if (blas_variant == BlasVariant.atlas) {
        return "atlas";
    } else if (blas_variant == BlasVariant.openblas) {
        return "openblas";
    } else if (blas_variant == BlasVariant.mkl) {
        return "mkl";
    } else if (blas_variant == BlasVariant.veclib) {
        return "veclib";
    } else if (blas_variant == BlasVariant.blis) {
        return "blis";
    } else {
        return "unknown";
    }
}

int blas_is_openblas()
{
    return blas_variant == BlasVariant.openblas;
}

/* for gretl_matrix.c, gretl_mt.c */

int blas_is_threaded()
{
    return blas_variant == BlasVariant.openblas ||
        blas_variant == BlasVariant.blis ||
        blas_variant == BlasVariant.mkl;
}

void blas_set_num_threads(int nt)
{
    if (OB_set_num_threads != null) {
        OB_set_num_threads(nt);
    } else if (BLIS_set_num_threads != null) {
        BLIS_set_num_threads(nt);
    } else if (MKL_domain_set_num_threads != null) {
        MKL_domain_set_num_threads(nt, 1); /* MKL_DOMAIN_BLAS */
    }
}

int blas_get_num_threads()
{
    if (OB_get_num_threads != null) {
        return OB_get_num_threads();
    } else if (BLIS_get_num_threads != null) {
        return BLIS_get_num_threads();
    } else if (MKL_domain_get_max_threads != null) {
        return MKL_domain_get_max_threads(1); /* MKL_DOMAIN_BLAS */
    } else {
        return 0;
    }
}

private void blas_init()
{
    void * ptr = null;

    blas_variant = BlasVariant.unknown;
    ptr = dlopen(null, RTLD_NOW);
    if (ptr == null) {
        goto try_ldd;
    }

    /* try for OpenBLAS */
    OB_set_num_threads = cast(ObSetNumThreadsFn) dlsym(ptr, "openblas_set_num_threads");
    OB_get_num_threads = cast(ObGetNumThreadsFn) dlsym(ptr, "openblas_get_num_threads");
    if (OB_set_num_threads != null) {
        blas_variant = BlasVariant.openblas;
        register_openblas_details(ptr);
        return;
    }

    /* try for Intel's MKL */
    MKL_finalize = cast(VoidFn) dlsym(ptr, "mkl_finalize");
    MKL_domain_get_max_threads = cast(MklDomainGetMaxThreadsFn) dlsym(ptr, "MKL_Domain_Get_Max_Threads");
    MKL_domain_set_num_threads = cast(MklDomainSetNumThreadsFn) dlsym(ptr, "MKL_Domain_Set_Num_Threads");
    if (MKL_domain_set_num_threads != null) {
        blas_variant = BlasVariant.mkl;
        register_mkl_details(ptr);
        return;
    }

    BLIS_init = cast(VoidFn) dlsym(ptr, "bli_init");
    BLIS_finalize = cast(VoidFn) dlsym(ptr, "bli_finalize");
    BLIS_set_num_threads = cast(BlisSetNumThreadsFn) dlsym(ptr, "bli_thread_set_num_threads");
    BLIS_get_num_threads = cast(BlisGetNumThreadsFn) dlsym(ptr, "bli_thread_get_num_threads");
    if (BLIS_init != null) {
        blas_variant = BlasVariant.blis;
        BLIS_init(); /* This is only to be sure that BLIS was initialized */
        register_blis_details(ptr);
        /* don't return in case libflame is used */
    }

    if (FLAME_init == null) {
        FLAME_init = cast(VoidFn) dlsym(ptr, "FLA_Init");
        if (FLAME_init != null) {
            FLAME_init();
            libflame_utils(ptr);
        }
    }

 try_ldd:
    if (blas_variant == BlasVariant.unknown) {
        blas_variant = detect_blas_via_ldd();
    }
}

void blas_cleanup()
{
    if (blas_variant == BlasVariant.blis && BLIS_finalize != null) {
        BLIS_finalize();
    } else if (blas_variant == BlasVariant.mkl && MKL_finalize != null) {
        MKL_finalize();
    }

    if (FLAME_initialized != null && FLAME_finalize != null) {
        if (FLAME_initialized()) {
            FLAME_finalize();
        }
    }
}

/* CPUID detection and analysis (for $sysinfo bundle) */

/* (the C code reads the CPUID instruction directly, via <cpuid.h>; here
   druntime's core.cpuid supplies the vendor and brand strings) */

char * get_cpu_details()
{
    char * ret = null;
    char[16] vendor_buf = 0;
    char[64] brand_buf = 0;

    strncat(vendor_buf.ptr, vendor.ptr, 15);
    strncat(brand_buf.ptr, processor.ptr, 63);

    if (brand_buf[0]) {
        ret = gretl_strdup(g_strstrip(brand_buf.ptr));
    } else if (vendor_buf[0]) {
        ret = gretl_strdup(g_strstrip(vendor_buf.ptr));
    } else {
        ret = gretl_strdup("unknown");
    }

    return ret;
}

/* AVX support detection */

int avx_support()
{
    if (avx512f) {
        return 512;
    } else if (avx2) {
        return 2;
    } else if (avx) {
        return 1;
    } else {
        return 0;
    }
}

/* library init and cleanup functions */

/**
 * libgretl_init:
 *
 * In a program that uses libgretl, this function should be
 * called once, before any other libgretl functions are
 * used. See also libgretl_cleanup(), and libgretl_mpi_init().
 **/

private __gshared bool library_started = false;

void libgretl_init()
{
    if (library_started) {
        return;     /* idempotent: dgretl.d calls this at program start-up */
    }
    library_started = true;

    libset_init();
    gretl_rand_init();
    gretl_xml_init();
    gretl_stopwatch_init();
    if (!gretl_in_tool_mode()) {
        blas_init();
        num_threads_init(blas_variant);
    }
}

int gretl_mpi_initialized()
{
    /* stub for non-MPI build */
    return 0;
}

private int other_gretl_running(const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != null) {
        /* Each entry in gretl.pid will contain two
           integers: PID and sequence number. So if
           we can read more than two values that must
           mean there's more than one gretl process
           running.
        */
        long[4] lv;
        int np;

        np = fscanf(fp, "%lld %lld %lld %lld", &lv[0],
                    &lv[1], &lv[2], &lv[3]);
        if (np > 2) {
            ret = 1;
        }
        fclose(fp);
    }

    return ret;
}

private void dotdir_cleanup()
{
    const(char) * dotdir = gretl_dotdir();
    const(char) * workdir = gretl_workdir();
    int err;

    if (!strcmp(dotdir, workdir)) {
        return;
    }

    err = gretl_chdir(dotdir);

    if (!err) {
        /* (the C code uses glib's GDir for the directory listing) */
        DIR * dir = opendir(".");
        dirent * ent;
        const(char) * fname;
        int skipit = 0;

        if (dir != null) {
            while ((ent = readdir(dir)) != null) {
                fname = ent.d_name.ptr;
                if (!strcmp(fname, ".") || !strcmp(fname, "..")) {
                    continue;   /* g_dir_read_name() never returns these */
                }
                if (!strcmp(fname, "gretl.pid")) {
                    skipit = other_gretl_running(fname);
                    break;
                }
            }
            if (!skipit) {
                rewinddir(dir);
                while ((ent = readdir(dir)) != null) {
                    fname = ent.d_name.ptr;
                    if (!strcmp(fname, ".") || !strcmp(fname, "..")) {
                        continue;
                    }
                    if (gretl_isdir(fname)) {
                        if (*fname == '.' && strlen(fname) > 3) {
                            /* failed auto-dot directory? */
                            gretl_deltree(fname);
                        }
                    } else if (!strncmp(fname, "prntmp.", 7) &&
                               !gretl_in_gui_mode()) {
                        ; /* don't let gretlcli delete a GUI file */
                    } else if (strcmp(fname, "..") &&
                               strcmp(fname, ".") &&
                               strcmp(fname, ".gretl2rc") &&
                               strcmp(fname, "gretl.pid") &&
                               strcmp(fname, "addons.idx") &&
                               strcmp(fname, "mail.dat")) {
                        gretl_remove(fname);
                    }
                }
            }
            closedir(dir);
        }
    }
}

void libgretl_session_cleanup(int mode)
{
    /* "last model" for accessors */
    gretl_saved_objects_cleanup();

    /* trash dataset-related items */
    gretl_transforms_cleanup();
    gretl_lists_cleanup();
    gretl_tests_cleanup();
    gretl_plotx(null, Gretlopt.OPT_NONE);

    /* scrub options set via "setopt" */
    setopt_cleanup();

    if (mode != SessionClear.dataset) {
        destroy_user_vars();
    }
}

/**
 * libgretl_cleanup:
 *
 * In a program that uses libgretl, this function may be
 * called to free various chunks of memory after the program
 * is finished with libgretl. Do not attempt to call any
 * other libgretl functions after invoking this cleanup.
 *
 * See also libgretl_init().
 **/

void libgretl_cleanup()
{
    if (!library_started) {
        return;     /* nothing to clean up, or already done */
    }
    library_started = false;

    libgretl_session_cleanup(SessionClear.all);

    gretl_functions_cleanup();
    libset_cleanup();
    gretl_command_hash_cleanup();
    gretl_function_hash_cleanup();
    lapack_mem_free();
    forecast_matrix_cleanup();
    stored_options_cleanup();
    option_printing_cleanup();
    gnuplot_cleanup();
    bufgets_cleanup();
    plugins_cleanup();
    gretl_bundle_cleanup();
    gretl_typemap_cleanup();
    gretl_stopwatch_cleanup();
    builtin_strings_cleanup();
    last_result_cleanup();

    dotdir_cleanup();

    gretl_script_dirs_cleanup();
    gretl_xml_cleanup();
    blas_cleanup();
}

/* record and retrieve hypothesis test results */

private enum TestRecordCode {
    set_test_stat,
    get_test_stat,
    get_test_pval,
    get_test_lnl,
    get_test_brk,
    tests_cleanup
}

private bool getcode(int c)
{
    return c == TestRecordCode.get_test_stat || c == TestRecordCode.get_test_pval ||
        c == TestRecordCode.get_test_lnl || c == TestRecordCode.get_test_brk;
}

private __gshared GretlType last_test_type;

private __gshared double rtr_val = NADBL;
private __gshared double rtr_pv = NADBL;
private __gshared double rtr_ll = NADBL;
private __gshared double rtr_brk = NADBL;

private double record_or_get_test_result(double teststat,
                                         double pval,
                                         double lnl,
                                         double inbrk,
                                         int code)
{
    double ret = NADBL;

    if (code == TestRecordCode.set_test_stat) {
        last_test_type = GretlType.double_;
        rtr_val = teststat;
        rtr_pv = pval;
        rtr_ll = lnl;
        rtr_brk = inbrk;
    } else if (getcode(code) && last_test_type == GretlType.double_) {
        if (code == TestRecordCode.get_test_stat) {
            ret = rtr_val;
        } else if (code == TestRecordCode.get_test_pval) {
            ret = rtr_pv;
        } else if (code == TestRecordCode.get_test_lnl) {
            ret = rtr_ll;
        } else if (code == TestRecordCode.get_test_brk) {
            ret = rtr_brk;
        }
    }

    return ret;
}

private __gshared GretlMatrixRaw * rtm_vals;
private __gshared GretlMatrixRaw * rtm_pvs;

private GretlMatrixRaw * record_or_get_test_matrix(GretlMatrixRaw * tests,
                                                GretlMatrixRaw * pvals,
                                                int code, int * err)
{
    GretlMatrixRaw * ret = null;

    if (code == TestRecordCode.tests_cleanup) {
        gretl_matrix_free(rtm_vals);
        gretl_matrix_free(rtm_pvs);
        rtm_vals = rtm_pvs = null;
        last_test_type = GretlType.none;
        return null;
    }

    if (code == TestRecordCode.set_test_stat) {
        last_test_type = GretlType.matrix;
        gretl_matrix_free(rtm_vals);
        rtm_vals = tests;
        gretl_matrix_free(rtm_pvs);
        rtm_pvs = pvals;
    } else if (getcode(code)) {
        GretlMatrixRaw * src = (code == TestRecordCode.get_test_stat)? rtm_vals : rtm_pvs;

        if (src != null) {
            ret = gretl_matrix_copy(src);
            if (ret == null) {
                *err = GretlError.alloc;
            }
        } else {
            *err = GretlError.badstat;
        }
    }

    return ret;
}

private void gretl_tests_cleanup()
{
    record_or_get_test_matrix(null, null, TestRecordCode.tests_cleanup, null);
}

/* Returns the type (scalar, matrix or none) of the last recorded
   test statistic/p-value pair.
*/

int get_last_test_type()
{
    return last_test_type;
}

void record_test_result(double teststat, double pval)
{
    record_or_get_test_result(teststat, pval, NADBL, NADBL, TestRecordCode.set_test_stat);
}

void record_LR_test_result(double teststat, double pval, double lnl)
{
    record_or_get_test_result(teststat, pval, lnl, NADBL, TestRecordCode.set_test_stat);
}

void record_QLR_test_result(double teststat, double pval, double brk)
{
    record_or_get_test_result(teststat, pval, NADBL, brk, TestRecordCode.set_test_stat);
}

/* Note: the hypothesis-test recorder "takes ownership" of the
   two input matrices, @tests and @pvals, which should therefore
   NOT be freed by the caller.
*/

void record_matrix_test_result(GretlMatrixRaw * tests,
                               GretlMatrixRaw * pvals)
{
    record_or_get_test_matrix(tests, pvals, TestRecordCode.set_test_stat, null);
}

double get_last_test_statistic()
{
    return record_or_get_test_result(0, 0, 0, 0, TestRecordCode.get_test_stat);
}

double get_last_pvalue()
{
    return record_or_get_test_result(0, 0, 0, 0, TestRecordCode.get_test_pval);
}

double get_last_lnl()
{
    return record_or_get_test_result(0, 0, 0, 0, TestRecordCode.get_test_lnl);
}

double get_last_break()
{
    return record_or_get_test_result(0, 0, 0, 0, TestRecordCode.get_test_brk);
}

GretlMatrixRaw * get_last_test_matrix(int * err)
{
    return record_or_get_test_matrix(null, null, TestRecordCode.get_test_stat, err);
}

GretlMatrixRaw * get_last_pvals_matrix(int * err)
{
    return record_or_get_test_matrix(null, null, TestRecordCode.get_test_pval, err);
}

/* (the C code's NEED_ALIGNED_MALLOC block, which is disabled, is not ported) */

private int is_executable(const(char) * s, uid_t myid, gid_t mygrp)
{
    stat_t buf;
    int ok = 0;

    if (gretl_stat(s, &buf) != 0) {
        return 0;
    }

    if (buf.st_mode & (S_IFREG|S_IFLNK)) {
        if (buf.st_uid == myid && (buf.st_mode & S_IXUSR)) {
            ok = 1;
        } else if (buf.st_gid == mygrp && (buf.st_mode & S_IXGRP)) {
            ok = 1;
        } else if (buf.st_uid != myid && buf.st_gid != mygrp &&
                   (buf.st_mode & S_IXOTH)) {
            ok = 1;
        }
    }

    return ok;
}

/**
 * check_for_program:
 * @prog: name of program.
 *
 * Returns: 1 if @prog is found in the PATH, otherwise 0.
 */

int check_for_program(const(char) * prog)
{
    uid_t myid = getuid();
    gid_t mygrp = getgid();
    char * path;
    char * pathcpy;
    char ** dirs;
    char * fullpath;
    char * p;
    int max_dlen = 0;
    int found = 0;
    int i, ndirs;

    if (prog == null || *prog == '\0') {
        return 0;
    }

    if (*prog == '/') {
        return is_executable(prog, myid, mygrp);
    }

    path = getenv("PATH");
    if (path == null || *path == '\0') {
        return 0;
    }

    pathcpy = gretl_strdup(path);
    if (pathcpy == null) {
        return 0;
    }

    ndirs = 1;
    p = pathcpy;
    while (*p) {
        if (*p == ':') ndirs++;
        p++;
    }

    dirs = cast(char **) malloc(ndirs * (*dirs).sizeof);
    if (dirs == null) {
        free(pathcpy);
        return 0;
    }

    if (ndirs == 1) {
        dirs[0] = pathcpy;
        max_dlen = cast(int) strlen(pathcpy);
    } else {
        for (i=0; i<ndirs; i++) {
            int dlen;

            dirs[i] = strtok((i == 0)? pathcpy : null, ":");
            if (dirs[i] == null) {
                ndirs = i;
                break;
            }
            dlen = cast(int) strlen(dirs[i]);
            if (dlen > max_dlen) {
                max_dlen = dlen;
            }
        }
    }

    if (ndirs == 0 ||
        (fullpath = cast(char *) malloc(max_dlen + strlen(prog) + 2)) == null) {
        free(dirs);
        free(pathcpy);
        return 0;
    }

    for (i=0; i<ndirs && !found; i++) {
        sprintf(fullpath, "%s/%s", dirs[i], prog);
        found = is_executable(fullpath, myid, mygrp);
    }

    free(dirs);
    free(pathcpy);
    free(fullpath);

    return found;
}

gint64 gretl_monotonic_time()
{
    return g_get_monotonic_time();
}

/* implementations of gzip and gunzip for unitary files */

private enum GRETL_ZBUFSIZ = 131072; /* 128k bytes */

int gretl_gzip(char * fname, char * zname)
{
    char[GRETL_ZBUFSIZ] buf;
    size_t len;
    FILE * fp;
    gzFile fz;

    fp = gretl_fopen(fname, "rb");
    if (fp == null) {
        return GretlError.fopen;
    }

    fz = gretl_gzopen(zname, "wb");
    if (fz == null) {
        fclose(fp);
        return GretlError.fopen;
    }

    while ((len = fread(buf.ptr, 1, GRETL_ZBUFSIZ, fp)) > 0) {
        gzwrite(fz, buf.ptr, cast(uint) len);
    }

    fclose(fp);
    gzclose(fz);

    return 0;
}

int gretl_gunzip(char * zname, char * fname)
{
    char[GRETL_ZBUFSIZ] buf;
    size_t len;
    FILE * fp;
    gzFile fz;

    fz = gretl_gzopen(zname, "rb");
    if (fz == null) {
        return GretlError.fopen;
    }

    fp = gretl_fopen(fname, "wb");
    if (fp == null) {
        gzclose(fz);
        return GretlError.fopen;
    }

    while ((len = gzread(fz, buf.ptr, GRETL_ZBUFSIZ)) > 0) {
        fwrite(buf.ptr, 1, len, fp);
    }

    gzclose(fz);
    fclose(fp);

    return 0;
}

/* Check if we're able to launch an automatic MPI routine, on the
   local machine, from within some time-consuming function.  This
   requires that MPI is available but not already running, and
   the local machine has at least 2 processors.
*/

int auto_mpi_ok()
{
    int ret = 0;

    /* (no MPI in this port) */

    return ret;
}

GretlMatrixRaw * dec2bin(double x, const(GretlMatrixRaw) * v, int * err)
{
    GretlMatrixRaw * ret = null;
    const(double) * val;
    uint ui;
    int n = 1;
    int i, j;

    if (v != null) {
        n = gretl_vector_get_length(v);
        if (n == 0) {
            *err = GretlError.invarg;
            return null;
        }
        val = v.val;
    } else {
        val = &x;
    }

    ret = gretl_zero_matrix_new(n, 32);
    if (ret == null) {
        *err = GretlError.alloc;
    }

    for (i=0; i<n; i++) {
        ui = gretl_uint32_from_double(val[i], err);
        if (*err) {
            break;
        }
        j = 0;
        while (ui > 0) {
            gretl_matrix_set(ret, i, j++, ui % 2);
            ui = ui >> 1;
        }
    }

    if (*err) {
        gretl_matrix_free(ret);
        ret = null;
    }

    return ret;
}

GretlMatrixRaw * bin2dec(const(GretlMatrixRaw) * m, int * err)
{
    int r = gretl_matrix_rows(m); /* must be 1 or more */
    int c = gretl_matrix_cols(m); /* must be >= 1 && <= 32 */
    GretlMatrixRaw * ret = null;

    if (r == 0 || c == 0 || c > 32) {
        *err = GretlError.invarg;
    } else {
        /* allocate column vector for return */
        ret = gretl_zero_matrix_new(r, 1);
        if (ret == null) {
            *err = GretlError.alloc;
        }
    }

    if (!*err) {
        uint ui, k;
        double mij;
        int i, j;

        for (i=0; i<r && !*err; i++) {
            ui = 0;
            k = 0x01;
            for (j=0; j<c && !*err; j++) {
                mij = gretl_matrix_get(m, i, j);
                if (mij == 1) {
                    ui += k;
                } else if (mij != 0) {
                    *err = GretlError.invarg;
                    break;
                }
                k = k << 1;
            }
            ret.val[i] = cast(double) ui;
        }
    }

    return ret;
}

/* apparatus for combining categorical indices */

private struct TriSorter {
    int a;
    int b;
    int pos;
}

private extern (C) int compare_ts(const(void) * a, const(void) * b)
{
    const(TriSorter) * ts1 = cast(const(TriSorter) *) a;
    const(TriSorter) * ts2 = cast(const(TriSorter) *) b;
    int ret = ts1.a - ts2.a;

    if (ret == 0) {
        ret = ts1.b - ts2.b;
    }

    return ret;
}

/* Given discrete series @v0 and @v1, fill series @v2 with
   indices of the combination of values of the first two
   series.
*/

int combine_categories(DatasetRaw * dset, int v0, int v1, int v2)
{
    TriSorter * tsr;
    int i, k;

    tsr = cast(TriSorter *) malloc(dset.n * (*tsr).sizeof);
    if (tsr == null) {
        return GretlError.alloc;
    }

    for (i=0; i<dset.n; i++) {
        tsr[i].a = cast(int) dset.Z[v0][i];
        tsr[i].b = cast(int) dset.Z[v1][i];
        tsr[i].pos = i;
    }

    qsort(tsr, dset.n, (*tsr).sizeof, &compare_ts);

    dset.Z[v2][tsr[0].pos] = k = 1;
    for (i=1; i<dset.n; i++) {
        if (tsr[i].b != tsr[i-1].b || tsr[i].a != tsr[i-1].a) {
            k++;
        }
        dset.Z[v2][tsr[i].pos] = k;
    }

    free(tsr);

    return 0;
}

/* Generate the permutations of @v from element i to n-1,
   with a shout-out to Daniel A. Jiménez, now at Texas A&M. See
   https://www.cs.utexas.edu/~djimenez/utsa/cs3343/lecture25.html
*/

private void perm(double * v, int n, int i, double * p, int * pk)
{
    int j, k = *pk;

    if (i == n) {
        for (j=0; j<n; j++) {
            p[k] = v[j];
            k++;
        }
        *pk = k;
    } else {
        /* recurse over the permutations from index i to n-1 */
        int tmp;

        for (j=i; j<n; j++) {
            /* swap elements i and j */
            tmp = cast(int) v[i];
            v[i] = v[j];
            v[j] = tmp;
            perm(v, n, i+1, p, pk);
            /* put them back as they were */
            tmp = cast(int) v[i];
            v[i] = v[j];
            v[j] = tmp;
        }
    }
}

/**
 * gretl_permute:
 * @v: m-vector (1 < m <= 8) whose elements should be permuted.
 * @err: location to receive error code.
 *
 * Returns: on successful completion, an m x n matrix, where
 * n = m! is the number of permutations of m items and each column
 * holds a distinct permutation. The underlying function is
 * recursive.
 */

GretlMatrixRaw * gretl_permute(GretlVector * v, int * err)
{
    GretlMatrixRaw * P;
    int m = gretl_vector_get_length(v);
    int n, i, k;

    if (m < 2 || m > 8) {
        /* maybe 8 is too cautious? */
        *err = GretlError.invarg;
        return null;
    }

    /* determine @m factorial */
    n = 1;
    for (i=2; i<=m; i++) {
        n *= i;
    }

    /* allocate storage */
    P = gretl_zero_matrix_new(m, n);
    if (P == null) {
        *err = GretlError.alloc;
        return null;
    }

    /* fill @P */
    k = 0;
    perm(v.val, m, 0, P.val, &k);

    return P;
}
