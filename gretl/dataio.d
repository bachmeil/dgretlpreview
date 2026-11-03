/* Port of dataio.c */

module gretl.dataio;


import core.stdc.ctype : isalpha, isdigit, isprint, isspace;
import core.stdc.errno : ERANGE, errno;
import core.stdc.math : ceil, fabs, floor, log10, lrint;
import core.stdc.stdio : EOF, FILE, fclose, fgetc, fgets, fprintf, fputc, fputs, fread, fscanf,
    printf, rewind, sprintf, sscanf, stderr, ungetc;
import core.stdc.config : c_long;
import core.stdc.stdlib : atoi, free, malloc, realloc, strtod, strtol;
import gretl.glib : gchar, gsize;
import core.stdc.string : memcpy, strcat, strchr, strcmp, strcpy, strlen, strncat, strrchr, strspn,
    strstr;
import etc.c.zlib : gzFile, gzclose, gzprintf, gzputc, gzputs;

import gretl.dataseth : SeriesTable, calendar_data, dataset_has_markers, dataset_is_cross_section, dataset_is_daily,
    dataset_is_decennial, dataset_is_hourly, dataset_is_panel, dataset_is_seasonal,
    dataset_is_time_series, dataset_is_weekly, dated_daily_data, dated_weekly_data,
    DatasetMarkerType, quarterly_or_monthly;
import gretl.matrix_extrah : MMissingCode;
import gretl.dataset : allocate_Z, clear_datainfo, copy_varinfo, create_new_dataset, datainfo_new, dataset_add_series, dataset_allocate_varnames, dataset_destroy_obs_markers, dataset_obs_info_default, destroy_dataset, free_Z, get_current_map, is_string_valued, merge_string_tables, series_attach_string_table, series_get_label, series_get_string_for_obs, series_get_string_table, series_is_coded, series_set_label, set_current_dataset, start_new_Z;
import gretl.dataioh : DataClearCode, GretlFileType, DBNA;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_set,
    gretl_warnmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : full_var_list;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs, printing_to_standard_stream;
import gretl.gerror : GError, g_error_free;
import gretl.gfileutils : g_file_get_contents, g_get_current_dir;
import gretl.gstrfuncs : g_ascii_strcasecmp, g_locale_to_utf8, g_strdup, g_strdup_printf, g_strstrip,
    g_utf8_validate;
import gretl.gretl_utils : date_as_double, positive_int_from_string;
import gretl.libgretlh;
import gretl.libset : get_csv_na_write_string, get_data_export_decpoint, get_data_export_delimiter, gretl_function_depth, gretl_messages_on, libset_get_int;
import gretl.libseth : SetKey;
import gretl.missingh : NADBL, na;
import gretl.options : get_optval_string, incompatible_options;
import gretl.strutils : colonize_obs, gretl_charsub, gretl_strdup, gretl_strndup, gretl_utf8_truncate_b,
    has_native_data_suffix, has_suffix, integer_string, is_greek_letter, modify_date_for_csv,
    obs_str_to_double, print_time, string_is_blank, strings_array_free,
    strings_array_new_with_length, strrslash;
import gretl.strutilsh : CTRLZ;
import gretl.dataioh;
import gretl.gretl_xml : gretl_read_gdt, gretl_write_gdt;
import gretl.gretl_data_io : gretl_read_foreign_data;
import gretl.csvdata : import_na_init;
import gretl.gretl_paths : gretl_addpath, gretl_chdir, gretl_fopen, gretl_gzopen, gretl_is_xml_file, gretl_maybe_switch_dir, gretl_test_fopen;
import gretl.plugins : get_plugin_function;
import gretl.calendar : calendar_date_string, calendar_obs_number, get_epoch_day, weekday_from_date, weekday_from_epoch_day;
import gretl.gretl_intl : get_local_decpoint, gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.genmain : current_series_index, gretl_is_series, gretl_reserved_word, series_index;
import gretl.matrix_extra : gretl_matrix_data_subset;
import gretl.gretl_array : get_array_by_name, get_strings_array_by_name, gretl_array_destroy, gretl_array_get_data, gretl_array_get_length, gretl_array_get_type, gretl_array_new, gretl_array_set_element;
import gretl.gretl_bundle : gretl_bundle_destroy, gretl_bundle_write_to_file;
import gretl.uservar : gretl_is_user_var, user_var_add_or_replace;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    int function_lookup(const(char) * s);
    void * get_user_function_by_name(const(char) * name);
    int is_function_alias(const(char) * name);
    int write_db_data(const(char) * fname, const(int) * list, Gretlopt opt, const(DatasetRaw) * dset);
}

private enum GretlDataFormat {
    gdt,       /* standard gretl XML data */
    binary,    /* native gretl binary data */
    octave,    /* data in Gnu Octave format */
    csv,       /* data in Comma Separated Values format */
    r,         /* data in Gnu R format */
    dat,       /* data in PcGive format */
    db,        /* gretl native database format */
    jm,        /* JMulti ascii data */
    dta,       /* Stata .dta format */
    json       /* geojson (maps) */
}

private bool IS_DATE_SEP(char c)
{
    return c == '.' || c == ':' || c == ',';
}

/**
 * SECTION:dataio
 * @short_description: data handling (internal)
 * @title: Data support
 * @include: gretl/libgretl.h
 *
 * The following data handling functions are basically internal to
 * gretl and not in a state where they can be readily
 * documented as public APIs.
 *
 */

/**
 * get_date_x:
 * @pd: frequency of data.
 * @obs: observation string.
 *
 * Returns: the floating-point representation of @obs.
 */

double get_date_x(int pd, const(char) * obs)
{
    double x = 1.0;

    if ((pd == 5 || pd == 6 || pd == 7 || pd == 52) && strlen(obs) > 4) {
        /* calendar data */
        uint ed = get_epoch_day(obs);

        if (ed > 0) {
            x = ed;
        }
    } else {
        x = obs_str_to_double(obs);
    }

    return x;
}

private int real_check_varname(const(char) * vname,
                               int is_series)
{
    int testchar = 'a';
    int firstbad = 0;
    int err = 0;

    gretl_error_clear();

    if (vname == null || *vname == '\0') {
        gretl_errmsg_set(_("Expected an identifier").ptr);
        return GretlError.parse;
    }

    if (strlen(vname) >= VNAMELEN) {
        gretl_errmsg_set(_("Varname exceeds the maximum of 31 characters").ptr);
        err = GretlError.data;
    } else if (strcmp(vname, "return") && gretl_reserved_word(vname)) {
        err = GretlError.data;
    } else if (!(isalpha(cast(ubyte) *vname))) {
        firstbad = 1;
        testchar = *vname;
        err = GretlError.data;
    } else if (is_series && (function_lookup(vname) ||
                             get_user_function_by_name(vname) != null ||
                             is_function_alias(vname))) {
        gretl_warnmsg_sprintf(_("'%s' shadows a function of the same name").ptr,
                              vname);
    } else {
        const(char) * p = vname;

        while (*p && testchar == 'a') {
            if (!(isalpha(cast(ubyte) *p))
                && !(isdigit(cast(ubyte) *p))
                && *p != '_') {
                testchar = *p;
                err = GretlError.data;
            }
            p++;
        }
    }

    if (err && !is_series && strlen(vname) == 2 && is_greek_letter(vname)) {
        return 0;
    }

    if (testchar != 'a') {
        if (isprint(cast(ubyte) testchar)) {
            if (firstbad) {
                gretl_errmsg_sprintf(_("First char of varname '%s' is bad\n"
                                       ~ "(first must be alphabetical)").ptr,
                                     vname);
            } else {
                gretl_errmsg_sprintf(_("Varname '%s' contains illegal character '%c'\n"
                                       ~ "Use only letters, digits and underscore").ptr,
                                     vname, cast(int) cast(ubyte) testchar);
            }
        } else {
            if (firstbad) {
                gretl_errmsg_sprintf(_("First char of varname (0x%x) is bad\n"
                                       ~ "(first must be alphabetical)").ptr,
                                     cast(uint) testchar);
            } else {
                gretl_errmsg_sprintf(_("Varname contains illegal character 0x%x\n"
                                       ~ "Use only letters, digits and underscore").ptr,
                                     cast(uint) testchar);
            }
        }
    }

    return err;
}

/**
 * check_varname:
 * @varname: putative name for variable (or object).
 *
 * Check a variable/object name for legality: the name
 * must start with a letter, and be composed of letters,
 * numbers or the underscore character, and nothing else.
 *
 * Returns: 0 if name is OK, non-zero if not.
 */

int check_varname(const(char) * varname)
{
    return real_check_varname(varname, 1);
}

int check_identifier(const(char) * varname)
{
    /* FIXME series case? */
    return real_check_varname(varname, 0);
}

private int bad_date_string(const(char) * s)
{
    int err = 0;

    gretl_error_clear();

    while (*s && !err) {
        if (!isdigit(cast(ubyte) *s) && !IS_DATE_SEP(*s)) {
            if (isprint(cast(ubyte) *s)) {
                gretl_errmsg_sprintf(_("Bad character '%c' in date string").ptr, cast(int) *s);
            } else {
                gretl_errmsg_sprintf(_("Bad character %d in date string").ptr, cast(int) *s);
            }
            err = 1;
        }
        s++;
    }

    return err;
}

/* Copy (and if applicable, unquote) @src into @targ, a buffer of size
   @targlen. The content is only ever written if it is known to fit,
   including the terminating NUL: if the label content (after stripping
   a matching pair of quotes, when present) is too long for @targlen,
   nothing is written to @targ and 1 is returned so the caller can treat
   this as "no match" rather than risking a buffer overflow. Returns 0
   on success (when @targ was filled in).
*/

private int maybe_unquote_label(char * targ, const(char) * src, int targlen)
{
    int quoted = (*src == '"' || *src == '\'');
    const(char) * content = quoted ? src + 1 : src;
    int n = cast(int) strlen(content);

    if (quoted && n > 0 && (content[n-1] == '"' || content[n-1] == '\'')) {
        /* the closing quote is not part of the label content */
        n--;
    }
    if (n >= targlen) {
        /* content (minus quotes) won't fit in @targ: bail out */
        return 1;
    }

    memcpy(targ, content, n);
    targ[n] = '\0';

    return 0;
}

private int get_dot_pos(const(char) * s)
{
    int i, pos = 0;

    for (i=0; *s != '\0'; i++, s++) {
        if (IS_DATE_SEP(*s)) {
            pos = i;
            break;
        }
    }

    return pos;
}

private int match_obs_marker(const(char) * s, const(DatasetRaw) * dset)
{
    char[OBSLEN] test;
    int t;

    if (maybe_unquote_label(test.ptr, s, cast(int) test.sizeof)) {
        /* the (unquoted) label content can't fit in test[], so it
           cannot match any marker string of length < OBSLEN
        */
        return -1;
    }

    for (t=0; t<dset.n; t++) {
        if (!strcmp(test.ptr, dset.S[t])) {
            /* handled */
            return t;
        }
    }

    return -1;
}

private int caldate_precheck(const(char) * s, int * y2, int * slashed)
{
    int n = cast(int) strlen(s);
    int ok = 0;

    if (n < 8 || !isdigit(cast(ubyte) s[0]) || !isdigit(cast(ubyte) s[1])) {
        /* can't be date */
        return 0;
    }

    if (n == 10) {
        if (s[4] == '-' && s[7] == '-') {
            ok = 1; /* could be YYYY-MM-DD */
        } else if (s[4] == '/' && s[7] == '/') {
            *slashed = 1; /* could be YYYY/MM/DD */
            ok = 1;
        }
    } else if (n == 8) {
        if (s[2] == '-' && s[5] == '-') {
            ok = 1; /* could be YY-MM-DD */
        } else if (s[2] == '/' && s[5] == '/') {
            *slashed = 1; /* could be YY/MM/DD */
            ok = 1;
        }
        if (ok) *y2 = 1;
    }

    return ok;
}

private int datecmp(const(char) * s1, int y21, int slash1,
                    const(char) * s2, int y22, int slash2)
{
    if (y21 && !y22) {
        s2 += 2;
    } else if (!y21 && y22) {
        s1 += 2;
    }
    if (slash1 && !slash2) {
        char[OBSLEN] revs1;

        strcpy(revs1.ptr, s1);
        gretl_charsub(revs1.ptr, '/', '-');
        return strcmp(revs1.ptr, s2);
    } else if (slash2 && !slash1) {
        char[OBSLEN] revs2;

        strcpy(revs2.ptr, s2);
        gretl_charsub(revs2.ptr, '/', '-');
        return strcmp(s1, revs2.ptr);
    } else {
        return strcmp(s1, s2);
    }
}

private int real_dateton(const(char) * date, const(DatasetRaw) * dset, int nolimit)
{
    int handled = 0;
    int t, n = -1;

    /* first check if this is calendar data and if so,
       treat accordingly */

    if (calendar_data(dset)) {
        if (dataset_has_markers(dset)) {
            /* "hard-wired" calendar dates as strings */
            int tryit, y21 = 0, slash1 = 0;
            int y22 = 0, slash2 = 0;

            tryit = caldate_precheck(date, &y21, &slash1);
            if (!tryit) {
                return -1;
            }
            tryit = caldate_precheck(dset.S[0], &y22, &slash2);
            if (!tryit) {
                return -1;
            }
            for (t=0; t<dset.n; t++) {
                if (!datecmp(date, y21, slash1, dset.S[t], y22, slash2)) {
                    /* handled */
                    return t;
                }
            }
            /* out of options: abort */
            return -1;
        } else {
            /* automatic calendar dates */
            n = calendar_obs_number(date, dset, nolimit);
            handled = 1;
        }
    } else if (dataset_is_daily(dset) ||
               dataset_is_weekly(dset)) {
        t = positive_int_from_string(date);
        if (t > 0) {
            n = t - 1;
            handled = 1;
        }
    } else if (dataset_is_decennial(dset)) {
        t = positive_int_from_string(date);
        if (t > 0) {
            n = cast(int) ((t - dset.sd0) / 10);
            handled = 1;
        }
    } else if (dataset_has_markers(dset)) {
        t = match_obs_marker(date, dset);
        if (t >= 0) {
            return t;
        }
        /* else maybe just a straight obs number */
        t = positive_int_from_string(date);
        if (t > 0) {
            n = t - 1;
            handled = 1;
        }
    }

    if (!handled) {
        int pos1, pos2;

        if (bad_date_string(date)) {
            return -1;
        }

        pos1 = get_dot_pos(date);
        pos2 = get_dot_pos(dset.stobs.ptr);

        if ((pos1 && !pos2) || (pos2 && !pos1)) {
            gretl_errmsg_sprintf(_("'%s': invalid observation index").ptr,
                                 date);
        } else if (!pos1 && !pos2) {
            n = atoi(date) - atoi(dset.stobs.ptr);
        } else if (pos1 > OBSLEN - 2) {
            gretl_errmsg_sprintf(_("'%s': invalid observation index").ptr,
                                 date);
        } else {
            char[OBSLEN] tmp;
            int maj, min;

            tmp[0] = '\0';
            strncat(tmp.ptr, date, OBSLEN-1);
            tmp[pos1] = '\0';
            maj = positive_int_from_string(tmp.ptr);
            min = positive_int_from_string(tmp.ptr + pos1 + 1);

            if (maj <= 0 || min <= 0 || min > dset.pd) {
                gretl_errmsg_sprintf(_("'%s': invalid observation index").ptr,
                                     date);
                n = -1;
            } else {
                int maj0, min0;

                tmp[0] = '\0';
                strncat(tmp.ptr, dset.stobs.ptr, OBSLEN-1);
                tmp[pos2] = '\0';
                maj0 = atoi(tmp.ptr);
                min0 = atoi(tmp.ptr + pos2 + 1);

                n = dset.pd * (maj - maj0) + (min - min0);
            }
        }
    }

    if (!nolimit && dset.n > 0 && n >= dset.n) {
        fprintf(stderr, "n = %d, dset->n = %d: out of bounds\n", n, dset.n);
        gretl_errmsg_set(_("Observation number out of bounds").ptr);
        n = -1;
    }

    return n;
}

/**
 * dateton:
 * @date: string representation of date for processing.
 * @dset: pointer to data information struct.
 *
 * Determines the observation number corresponding to @date,
 * relative to @dset. It is an error if @date represents an
 * observation that lies outside of the full data range
 * specified in @dset.
 *
 * Returns: zero-based observation number, or -1 on error.
 */

int dateton(const(char) * date, const(DatasetRaw) * dset)
{
    return real_dateton(date, dset, 0);
}

/**
 * merge_dateton:
 * @date: string representation of date for processing.
 * @dset: pointer to data information struct.
 *
 * Works just as dateton(), except that for this function it
 * is not an error if @date represents an observation that
 * lies beyond the data range specified in @dset. This is
 * intended for use when merging data, or when creating a new
 * dataset.
 *
 * Returns: zero-based observation number, or -1 on error.
 */

int merge_dateton(const(char) * date, const(DatasetRaw) * dset)
{
    return real_dateton(date, dset, 1);
}

private char * panel_obs(char * s, int t, const(DatasetRaw) * dset)
{
    int i = t / dset.pd + 1;
    int j = (t + 1) % dset.pd;
    int d = cast(int) (1 + floor(log10(cast(double) dset.pd)));

    if (j == 0) {
        j = dset.pd;
    }

    sprintf(s, "%d:%0*d", i, d, j);

    return s;
}

/**
 * ntolabel:
 * @datestr: char array to which date is to be printed.
 * @t: zero-based observation number.
 * @dset: data information struct.
 *
 * Prints to @datestr (which must be at least #OBSLEN bytes)
 * the calendar representation of observation number @t.
 *
 * Returns: the observation string.
 */

char * ntolabel(char * datestr, int t, const(DatasetRaw) * dset)
{
    double x;

    if (calendar_data(dset)) {
        /* handles both daily and dated weekly data */
        if (dataset_has_markers(dset)) {
            strcpy(datestr, dset.S[t]);
            if (strchr(datestr, '/')) {
                gretl_charsub(datestr, '/', '-');
            }
        } else {
            calendar_date_string(datestr, t, dset);
        }
        return datestr;
    } else if (dataset_is_daily(dset) ||
               dataset_is_weekly(dset)) {
        /* undated time series */
        x = date_as_double(t, 1, dset.sd0);
        sprintf(datestr, "%d", cast(int) x);
        return datestr;
    } else if (dataset_is_decennial(dset)) {
        x = dset.sd0 + 10 * t;
        sprintf(datestr, "%d", cast(int) x);
        return datestr;
    } else if (dataset_is_panel(dset)) {
        panel_obs(datestr, t, dset);
        return datestr;
    } else if (dataset_is_hourly(dset) && dataset_has_markers(dset)) {
        strcpy(datestr, dset.S[t]);
        return datestr;
    }

    x = date_as_double(t, dset.pd, dset.sd0);

    if (dset.pd == 1) {
        sprintf(datestr, "%d", cast(int) x);
    } else {
        int pdp = dset.pd;
        short len = 1;
        char[10] fmt;

        while ((pdp = pdp / 10) != 0) len++;
        sprintf(fmt.ptr, "%%.%df", cast(int) len);
        sprintf(datestr, fmt.ptr, x);
        colonize_obs(datestr);
    }

    return datestr;
}

/* print observation date in ISO 8601 extended format */

char * ntolabel_8601(char * datestr, int t, const(DatasetRaw) * dset)
{
    *datestr = '\0';

    if (calendar_data(dset)) {
        if (dataset_has_markers(dset)) {
            strcpy(datestr, dset.S[t]);
        } else {
            calendar_date_string(datestr, t, dset);
        }
    } else if (dataset_is_decennial(dset)) {
        double x = dset.sd0 + 10 * t;
        int yr = cast(int) lrint(x);

        sprintf(datestr, "%d-01-01", yr);
    } else {
        double x = date_as_double(t, dset.pd, dset.sd0);
        int maj = cast(int) lrint(floor(x));

        if (dset.pd == 1) {
            sprintf(datestr, "%d-01-01", maj);
        } else if (dset.pd == 12) {
            int min = cast(int) lrint(100 * (x - floor(x)));

            sprintf(datestr, "%d-%02d-01", maj, min);
        } else if (dset.pd == 4) {
            int min = cast(int) lrint(10 * (x - floor(x)));
            int mo = min==2 ? 4 : min==3? 7 : min==4? 10 : min;

            sprintf(datestr, "%d-%02d-01", maj, mo);
        }
    }

    return datestr;
}

private double xround(double x)
{
    return ((x-floor(x))>.5)? ceil(x) : floor(x);
}

/**
 * get_subperiod:
 * @t: zero-based observation number.
 * @dset: data information struct.
 * @err: location to receive error code, or NULL.
 *
 * For "seasonal" time series data (in a broad sense),
 * determines the sub-period at observation @t. The "sub-period"
 * might be a quarter, month, hour or whatever.  The value
 * returned is zero-based (e.g. first quarter = 0).
 * If the data are not "seasonal", 0 is returned and if
 * @err is non-NULL it receives a non-zero error code.
 *
 * Returns: the sub-period.
 */

int get_subperiod(int t, const(DatasetRaw) * dset, int * err)
{
    int ret = 0;

    if (!dataset_is_seasonal(dset)) {
        if (err != null) {
            *err = GretlError.pdwrong;
        }
        return 0;
    }

    if (dataset_is_weekly(dset)) {
        /* bodge -- what else to do? */
        ret = t % dset.pd;
    } else if (dated_daily_data(dset)) {
        char[OBSLEN] obs;

        ntolabel(obs.ptr, t, dset);
        ret = weekday_from_date(obs.ptr);
        /* 2023-08-06: for now, ensure backward compatibility,
           though I'm not sure that's actually right.
        */
        ret = ret == 7 ? 0 : ret;
    } else if (dataset_is_daily(dset)) {
        /* bodge, again */
        ret = t % dset.pd;
    } else {
        /* quarterly, monthly, hourly... */
        double x = date_as_double(t, dset.pd, dset.sd0);
        int i, d = cast(int) ceil(log10(cast(double) dset.pd));

        x -= floor(x);
        for (i=0; i<d; i++) {
            x *= 10;
        }
        ret = cast(int) xround(x) - 1;
    }

    return ret;
}

/**
 * get_precision:
 * @x: data array.
 * @n: length of @x.
 * @placemax: the maximum number of decimal places to try.
 *
 * Find the number of decimal places required to represent a given
 * data series uniformly and accurately, if possible.
 *
 * Returns: the required number of decimal places or
 * #PMAX_NOT_AVAILABLE if it can't be done.
 */

int get_precision(const(double) * x, int n, int placemax)
{
    int t, p, pmax = 0;
    char * s;
    char[64] numstr;
    int len, n_ok = 0;
    double z;

    for (t=0; t<n; t++) {
        if (!na(x[t])) {
            z = fabs(x[t]);
            /* escape clause: numbers are too big or too small for
               this treatment */
            if (z > 0 && (z < 1.0e-6 || z > 1.0e+8)) {
                return PMAX_NOT_AVAILABLE;
            }
            n_ok++;
        }
    }

    if (n_ok == 0) {
        return PMAX_NOT_AVAILABLE;
    }

    for (t=0; t<n; t++) {
        if (!na(x[t])) {
            p = placemax;
            sprintf(numstr.ptr, "%.*f", p, fabs(x[t]));
            /* go to the end and drop trailing zeros */
            len = cast(int) strlen(numstr.ptr);
            s = numstr.ptr + len - 1;
            while (*s-- == '0') {
                p--;
                len--;
            }
            if (len > 10) {
                /* this is going to be too big */
                return PMAX_NOT_AVAILABLE;
            }
            if (p > pmax) {
                pmax = p;
            }
        }
    }

    return pmax;
}

private struct ExtMap {
    GretlFileType ftype;
    string ext;
}

private immutable ExtMap[12] data_ftype_map = [
    ExtMap(GretlFileType.xml_data,     ".gdt"),
    ExtMap(GretlFileType.binary_data,  ".gdtb"),
    ExtMap(GretlFileType.csv,          ".csv"),
    ExtMap(GretlFileType.octave,       ".m"),
    ExtMap(GretlFileType.gnumeric,     ".gnumeric"),
    ExtMap(GretlFileType.xls,          ".xls"),
    ExtMap(GretlFileType.xlsx,         ".xlsx"),
    ExtMap(GretlFileType.ods,          ".ods"),
    ExtMap(GretlFileType.wf1,          ".wf1"),
    ExtMap(GretlFileType.dta,          ".dta"),
    ExtMap(GretlFileType.sav,          ".sav"),
    ExtMap(GretlFileType.sas,          ".xpt")
];

private immutable string[4] map_suffixes = [
    ".json",
    ".geojson",
    ".dbf",
    ".shp"
];

private const(char) * get_filename_extension(const(char) * fname)
{
    const(char) * ext = strrchr(fname, '.');

    if (ext != null && strchr(ext, '/')) {
        /* the rightmost dot is not in the basename */
        ext = null;
    }

    return ext;
}

private GretlFileType dat_filetype(const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "rb");
    GretlFileType ret = GretlFileType.csv;

    if (fp != null) {
        char[3] test = 0;

        if (fread(test.ptr, 1, 2, fp) == 2 && !strcmp(test.ptr, "/*")) {
            ret = GretlFileType.jmulti;
        }
        fclose(fp);
    }

    return ret;
}

private GretlFileType data_file_type_from_extension(const(char) * ext,
                                                    const(char) * fname)
{
    int i, n = cast(int) data_ftype_map.length;

    for (i=0; i<n; i++) {
        if (!g_ascii_strcasecmp(ext, data_ftype_map[i].ext.ptr)) {
            return data_ftype_map[i].ftype;
        }
    }

    /* try to handle JMulTi */
    if (!g_ascii_strcasecmp(ext, ".dat")) {
        return dat_filetype(fname);
    }

    /* a few extras */
    if (!g_ascii_strcasecmp(ext, ".txt") ||
        !g_ascii_strcasecmp(ext, ".asc") ||
        !g_ascii_strcasecmp(ext, ".gz")) {
        return GretlFileType.csv;
    }

    /* map metadata */
    n = cast(int) map_suffixes.length;
    for (i=0; i<n; i++) {
        if (!g_ascii_strcasecmp(ext, map_suffixes[i].ptr)) {
            return GretlFileType.map;
        }
    }

    return GretlFileType.unrecognized;
}

GretlFileType data_file_type_from_name(const(char) * fname)
{
    GretlFileType ret = GretlFileType.unrecognized;
    const(char) * ext = strrchr(fname, '.');

    if (ext != null && strchr(ext, '/')) {
        /* the rightmost dot is not in the basename */
        ext = null;
    }
    if (ext != null) {
        ret = data_file_type_from_extension(ext, fname);
    }

    return ret;
}

private bool non_native(Gretlopt o)
{
    return (o & (Gretlopt.OPT_M | Gretlopt.OPT_R | Gretlopt.OPT_D |
                 Gretlopt.OPT_G | Gretlopt.OPT_J)) != 0;
}

private GretlDataFormat format_from_opt_or_name(Gretlopt opt, const(char) * fname,
                                                char * delim, int * add_ext,
                                                int * gzip, int * err)
{
    GretlDataFormat fmt = GretlDataFormat.gdt;

    if (has_suffix(fname, ".gdt")) {
        if (non_native(opt)) {
            *err = GretlError.badopt;
        }
        return GretlDataFormat.gdt;
    } else if (has_suffix(fname, ".gdtb")) {
        if (non_native(opt)) {
            *err = GretlError.badopt;
        }
        return GretlDataFormat.binary;
    } else if (has_suffix(fname, ".geojson")) {
        return GretlDataFormat.json;
    }

    if (opt & Gretlopt.OPT_M) {
        fmt = GretlDataFormat.octave;
    } else if (opt & Gretlopt.OPT_R) {
        fmt = GretlDataFormat.r;
    } else if (opt & Gretlopt.OPT_D) {
        fmt = GretlDataFormat.db;
    } else if (opt & Gretlopt.OPT_G) {
        fmt = GretlDataFormat.dat;
    } else if (opt & Gretlopt.OPT_J) {
        fmt = GretlDataFormat.jm;
    }

    if (fmt == GretlDataFormat.gdt) {
        if (has_suffix(fname, ".R")) {
            fmt = GretlDataFormat.r;
        } else if (has_suffix(fname, ".csv")) {
            fmt = GretlDataFormat.csv;
        } else if (has_suffix(fname, ".csv.gz")) {
            fmt = GretlDataFormat.csv;
            *gzip = 1;
        } else if (has_suffix(fname, ".m")) {
            fmt = GretlDataFormat.octave;
        } else if (has_suffix(fname, ".txt") ||
                   has_suffix(fname, ".asc")) {
            fmt = GretlDataFormat.csv;
            *delim = ' ';
        } else if (has_suffix(fname, ".dta")) {
            fmt = GretlDataFormat.dta;
        } else if (has_suffix(fname, ".bin")) {
            fmt = GretlDataFormat.db;
        }
    }

    if (fmt == GretlDataFormat.gdt) {
        *add_ext = 1;
    }

    return fmt;
}

void date_maj_min(int t, const(DatasetRaw) * dset, int * maj, int * min)
{
    char[OBSLEN] obs;

    ntolabel(obs.ptr, t, dset);

    if (maj != null) {
        *maj = atoi(obs.ptr);
    }

    if (min != null) {
        char * s;
        char sep = ':';

        if (strchr(obs.ptr, sep) == null) {
            if (dset.pd == 4 && strchr(obs.ptr, 'Q')) {
                sep = 'Q';
            } else if (dset.pd == 12 && strchr(obs.ptr, 'M')) {
                sep = 'M';
            }
        }

        s = strchr(obs.ptr, sep);
        if (s != null && strlen(s) > 1) {
            *min = atoi(s + 1);
        } else {
            *min = 1;
        }
    }
}

private enum TMPLEN = 64;

private void csv_data_out(const(DatasetRaw) * dset, const(int) * list,
                          int print_obs, int digits, char decpoint,
                          char delim, FILE * fp, gzFile fz)
{
    const(char) * NA = get_csv_na_write_string();
    char[TMPLEN] tmp;
    double xt;
    int popit = 0, dotsub = 0;
    int t, i, vi;

    if (decpoint == '.' && get_local_decpoint() == ',') {
        gretl_push_c_numeric_locale();
        popit = 1;
    } else if (decpoint == ',' && get_local_decpoint() == '.') {
        dotsub = 1;
    }

    for (t=dset.t1; t<=dset.t2; t++) {
        if (print_obs) {
            if (dset.S != null) {
                if (fz != null) {
                    gzprintf(fz, "\"%s\"%c", dset.S[t], cast(int) delim);
                } else {
                    fprintf(fp, "\"%s\"%c", dset.S[t], cast(int) delim);
                }
            } else {
                ntolabel(tmp.ptr, t, dset);
                if (quarterly_or_monthly(dset)) {
                    modify_date_for_csv(tmp.ptr, dset.pd);
                }
                if (fz != null) {
                    gzprintf(fz, "%s%c", tmp.ptr, cast(int) delim);
                } else {
                    fprintf(fp, "%s%c", tmp.ptr, cast(int) delim);
                }
            }
        }

        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            xt = dset.Z[vi][t];
            if (na(xt)) {
                fputs(NA, fp);
            } else {
                if (is_string_valued(dset, vi)) {
                    const(char) * st;

                    st = series_get_string_for_obs(dset, vi, t);
                    if (st != null) {
                        tmp[0] = '\0';
                        strcat(tmp.ptr, "\"");
                        strncat(tmp.ptr, st, TMPLEN - 3);
                        strcat(tmp.ptr, "\"");
                    } else {
                        fprintf(stderr, "missing string at t=%d, vi=%d, xt=%g\n",
                                t, vi, xt);
                        strcpy(tmp.ptr, "\"NA\"");
                    }
                } else if (series_is_coded(dset, vi)) {
                    sprintf(tmp.ptr, "\"%d\"", cast(int) xt);
                } else {
                    sprintf(tmp.ptr, "%.*g", digits, xt);
                }
                if (dotsub) {
                    gretl_charsub(tmp.ptr, '.', ',');
                }
                if (fz != null) {
                    gzputs(fz, tmp.ptr);
                } else {
                    fputs(tmp.ptr, fp);
                }
            }
            if (fz != null) {
                gzputc(fz, i < list[0] ? delim : '\n');
            } else {
                fputc(i < list[0] ? delim : '\n', fp);
            }
        }
    }

    if (popit) {
        gretl_pop_c_numeric_locale();
    }
}

private int markers_are_unique(const(DatasetRaw) * dset)
{
    int t, s;

    for (t=dset.t1; t<dset.t2; t++) {
        for (s=t+1; s<=dset.t2; s++) {
            if (strcmp(dset.S[t], dset.S[s]) == 0) {
                return 0;
            }
        }
    }

    return 1;
}

private void R_data_out(const(DatasetRaw) * dset, const(int) * list,
                        int digits, FILE * fp)
{
    int print_markers = 0;
    double xt;
    int t, i, vi;

    if (dset.S != null) {
        print_markers = markers_are_unique(dset);
    }

    for (t=dset.t1; t<=dset.t2; t++) {
        if (print_markers) {
            fprintf(fp, "\"%s\" ", dset.S[t]);
        }
        for (i=1; i<=list[0]; i++) {
            vi = list[i];
            xt = dset.Z[vi][t];
            if (na(xt)) {
                fputs("NA", fp);
            } else if (is_string_valued(dset, vi)) {
                fprintf(fp, "\"%s\"", series_get_string_for_obs(dset, vi, t));
            } else if (series_is_coded(dset, vi)) {
                fprintf(fp, "\"%d\"", cast(int) xt);
            } else {
                fprintf(fp, "%.*g", digits, xt);
            }
            fputc(i < list[0] ? ' ' : '\n', fp);
        }
    }
}

private int write_dta_data(const(char) * fname, const(int) * list,
                           Gretlopt opt, const(DatasetRaw) * dset)
{
    alias ExporterFn = extern (C) int function(const(char) *, const(int) *, Gretlopt,
                                               const(DatasetRaw) *);
    ExporterFn exporter;
    int err = 0;

    exporter = cast(ExporterFn) get_plugin_function("stata_export");

    if (exporter == null) {
        err = 1;
    } else {
        err = (*exporter)(fname, list, opt, dset);
    }

    return err;
}

private int write_map_data(const(char) * fname,
                           const(int) * list,
                           const(DatasetRaw) * dset)
{
    GretlBundle * b = null;
    int err = 0;

    b = get_current_map(dset, list, &err);

    if (!err) {
        err = gretl_bundle_write_to_file(b, fname, 0);
    }

    gretl_bundle_destroy(b);

    return err;
}

private enum DEFAULT_CSV_DIGITS = 15;

private int real_write_data(const(char) * fname, int * list,
                            const(DatasetRaw) * dset, Gretlopt opt,
                            int progress, Prn * prn)
{
    int i, t, v, l0;
    GretlDataFormat fmt;
    char[MAXLEN] datfile;
    int n = dset.n;
    int pop_locale = 0;
    char delim = 0;
    FILE * fp = null;
    gzFile fz = null;
    int freelist = 0;
    int csv_digits = 0;
    int add_ext = 0;
    int gzip = 0;
    double xx;
    int err = 0;

    gretl_error_clear();

    if (list != null && list[0] == 0) {
        return GretlError.args;
    }

    fmt = format_from_opt_or_name(opt, fname, &delim, &add_ext,
                                  &gzip, &err);
    if (err) {
        return err;
    }

    if (list == null) {
        list = full_var_list(dset, &l0);
        if (l0 == 0) {
            return GretlError.args;
        } else if (list == null) {
            return GretlError.alloc;
        } else {
            freelist = 1;
        }
    }

    l0 = list[0];
    fname = gretl_maybe_switch_dir(fname);

    if (fmt == GretlDataFormat.gdt || fmt == GretlDataFormat.binary) {
        /* write native data file (.gdt or .gdtb) */
        err = gretl_write_gdt(fname, list, dset, opt, progress);
        goto write_exit;
    }

    if (fmt == GretlDataFormat.db) {
        /* native type database file */
        err = write_db_data(fname, list, opt, dset);
        goto write_exit;
    }

    if (fmt == GretlDataFormat.dta) {
        /* Stata */
        err = write_dta_data(fname, list, opt, dset);
        goto write_exit;
    }

    if (fmt == GretlDataFormat.json) {
        /* writing map as geojson */
        err = write_map_data(fname, list, dset);
        goto write_exit;
    }

    if (strlen(fname) >= MAXLEN) {
        gretl_errmsg_sprintf(_("filename too long (%d bytes)").ptr,
                             cast(int) strlen(fname));
        err = GretlError.data;
        goto write_exit;
    }
    strcpy(datfile.ptr, fname);

    /* open file for output */
    if (gzip) {
        fz = gretl_gzopen(datfile.ptr, "wb1");
    } else {
        fp = gretl_fopen(datfile.ptr, "wb");
    }
    if (fp == null && fz == null) {
        err = GretlError.fopen;
        goto write_exit;
    }

    csv_digits = libset_get_int(SetKey.csv_digits);

    if (csv_digits <= 0) {
        csv_digits = DEFAULT_CSV_DIGITS;
    }

    if (fmt != GretlDataFormat.csv) {
        /* ensure C locale for data output */
        gretl_push_c_numeric_locale();
        pop_locale = 1;
    }

    if (fmt == GretlDataFormat.csv) {
        const(char) * msg = get_optval_string(GretlCmdIndex.store, Gretlopt.OPT_E);
        char decpoint = get_data_export_decpoint();
        int print_obs = 0;

        if (opt & Gretlopt.OPT_I) {
            /* the CSV --decimal-comma option */
            decpoint = ',';
            delim = ';';
        } else if (delim == 0) {
            delim = get_data_export_delimiter();
        }

        if (msg != null && *msg != '\0') {
            if (gzip) {
                gzprintf(fz, "# %s\n", msg);
            } else {
                fprintf(fp, "# %s\n", msg);
            }
        }

        if (!(opt & Gretlopt.OPT_X)) {
            /* OPT_X prohibits printing of observation strings */
            print_obs = dataset_is_time_series(dset) || dset.S != null;
        }

        if (!(opt & Gretlopt.OPT_N)) {
            /* header: variable names */
            if (print_obs && (dset.S != null || dset.structure != TsCodes.cross_section)) {
                if (gzip) {
                    gzprintf(fz, "obs%c", cast(int) delim);
                } else {
                    fprintf(fp, "obs%c", cast(int) delim);
                }
            }
            for (i=1; i<l0; i++) {
                if (gzip) {
                    gzprintf(fz, "%s%c", dset.varname[list[i]], cast(int) delim);
                } else {
                    fprintf(fp, "%s%c", dset.varname[list[i]], cast(int) delim);
                }
            }
            if (gzip) {
                gzprintf(fz, "%s\n", dset.varname[list[l0]]);
            } else {
                fprintf(fp, "%s\n", dset.varname[list[l0]]);
            }
        }

        csv_data_out(dset, list, print_obs, csv_digits,
                     decpoint, delim, fp, fz);
    } else if (fmt == GretlDataFormat.r) {
        /* friendly to GNU R */
        if (dataset_is_time_series(dset)) {
            char[OBSLEN] datestr;

            ntolabel(datestr.ptr, dset.t1, dset);
            fprintf(fp, "# time-series data: start = %s, frequency = %d\n",
                    datestr.ptr, dset.pd);
        }

        for (i=1; i<l0; i++) {
            fprintf(fp, "%s ", dset.varname[list[i]]);
        }
        fprintf(fp, "%s\n", dset.varname[list[l0]]);

        R_data_out(dset, list, csv_digits, fp);
    } else if (fmt == GretlDataFormat.octave) {
        /* GNU Octave: write out data as several matrices (one per
           series) in the same file */

        for (i=1; i<=list[0]; i++) {
            v = list[i];
            fprintf(fp, "# name: %s\n# type: matrix\n# rows: %d\n# columns: 1\n",
                    dset.varname[v], n);
            for (t=dset.t1; t<=dset.t2; t++) {
                xx = dset.Z[v][t];
                if (na(xx)) {
                    fputs("NaN ", fp);
                } else {
                    fprintf(fp, "%.*g\n", csv_digits, xx);
                }
            }
        }
    } else if (fmt == GretlDataFormat.dat) {
        /* PcGive: data file with load info */
        int pd = dset.pd;

        for (i=1; i<=list[0]; i++) {
            fprintf(fp, ">%s ", dset.varname[list[i]]);
            if (dset.structure == TsCodes.time_series &&
                (pd == 1 || pd == 4 || pd == 12)) {
                int maj, min;

                date_maj_min(dset.t1, dset, &maj, &min);
                fprintf(fp, "%d %d ", maj, min);
                date_maj_min(dset.t2, dset, &maj, &min);
                fprintf(fp, "%d %d %d", maj, min, pd);
            } else {
                fprintf(fp, "%d 1 %d 1 1", dset.t1, dset.t2);
            }

            fputc('\n', fp);

            for (t=dset.t1; t<=dset.t2; t++) {
                v = list[i];
                xx = dset.Z[v][t];
                if (na(xx)) {
                    fprintf(fp, "-9999.99");
                } else {
                    fprintf(fp, "%.*g", csv_digits, xx);
                }
                fputc('\n', fp);
            }
            fputc('\n', fp);
        }
    } else if (fmt == GretlDataFormat.jm) {
        /* JMulti: ascii with comments and date info */
        const(char) * vlabel;
        int maj, min;

        fputs("/*\n", fp);
        for (i=1; i<=list[0]; i++) {
            v = list[i];
            vlabel = series_get_label(dset, v);
            fprintf(fp, " %s: %s\n", dset.varname[v],
                    vlabel == null ? "".ptr : vlabel);
        }
        fputs("*/\n", fp);
        date_maj_min(dset.t1, dset, &maj, &min);
        if (dset.pd == 4 || dset.pd == 12) {
            fprintf(fp, "<%d %c%d>\n", maj, (dset.pd == 4)? cast(int) 'Q' : cast(int) 'M', min);
        } else if (dset.pd == 1) {
            fprintf(fp, "<%d>\n", maj);
        } else {
            fputs("<1>\n", fp);
        }
        for (i=1; i<=list[0]; i++) {
            v = list[i];
            fprintf(fp, " %s", dset.varname[v]);
        }
        fputc('\n', fp);
        for (t=dset.t1; t<=dset.t2; t++) {
            for (i=1; i<=list[0]; i++) {
                v = list[i];
                if (na(dset.Z[v][t])) {
                    fputs("NaN ", fp);
                } else {
                    fprintf(fp, "%.*g ", csv_digits, dset.Z[v][t]);
                }
            }
            fputc('\n', fp);
        }
    }

    if (pop_locale) {
        gretl_pop_c_numeric_locale();
    }

    if (fz != null) {
        gzclose(fz);
    } else if (fp != null) {
        fclose(fp);
    }

 write_exit:

    if (!err && prn != null && gretl_messages_on()) {
        if (add_ext) {
            pprintf(prn, _("wrote %s.gdt\n").ptr, fname);
        } else {
            pprintf(prn, _("wrote %s\n").ptr, fname);
        }
    }

    if (freelist) {
        free(list);
    }

    return err;
}

/**
 * write_data:
 * @fname: name of file to write.
 * @list: list of variables to write (or %NULL to write all series).
 * @dset: dataset struct.
 * @opt: option flag indicating format in which to write the data.
 * @prn: gretl printer or NULL.
 *
 * Write out a data file containing the values of the given set
 * of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int write_data(const(char) * fname, int * list, const(DatasetRaw) * dset,
               Gretlopt opt, Prn * prn)
{
    return real_write_data(fname, list, dset, opt, 0, prn);
}

int gui_write_data(const(char) * fname, int * list, const(DatasetRaw) * dset,
                   Gretlopt opt)
{
    return real_write_data(fname, list, dset, opt, 1, null);
}

/**
 * is_gzipped:
 * @fname: filename to examine.
 *
 * Determine if the given file is gzipped.
 *
 * Returns: 1 in case of a gzipped file, 0 if not gzipped or
 * inaccessible.
 *
 */

int is_gzipped(const(char) * fname)
{
    FILE * fp;
    int gz = 0;

    if (fname == null || *fname == '\0') {
        return 0;
    }

    fp = gretl_fopen(fname, "rb");
    if (fp == null) {
        return 0;
    }

    if (fgetc(fp) == 0x1f && fgetc(fp) == 0x8b) {   /* octal 037, 0213 */
        gz = 1;
    }

    fclose(fp);

    return gz;
}

/**
 * gretl_seek_data:
 * @fname: name of file to try.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: where messages should be written.
 *
 * Reads data from file (either "native" or "CSV") into
 * gretl's workspace, allocating storage as required.
 *
 * Path-searching is performed on @fname, and the function
 * will try adding the .gdt or .gdtb extension to @fname
 * if this is not given.
 *
 * The only applicable option is that @opt may contain
 * OPT_T when appending data to a panel dataset: in
 * that case we try to interpret the new data as time
 * series, in common across all panel units. In most
 * cases, just give OPT_NONE.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_seek_data(char * fname, DatasetRaw * dset,
                    Gretlopt opt, Prn * prn)
{
    Gretlopt myopt = Gretlopt.OPT_NONE;
    char * test;
    int err = 0;

    gretl_error_clear();
    test = gretl_addpath(fname, 0);
    if (test == null) {
        return GretlError.fopen;
    }

    if (opt & Gretlopt.OPT_T) {
        myopt = Gretlopt.OPT_T;
    }

    if (has_native_data_suffix(fname)) {
        /* specific processing for gretl datafiles */
        err = gretl_read_gdt(fname, dset, myopt, prn);
    } else {
        /* all other file types */
        GretlFileType ft = data_file_type_from_name(fname);

        if (ft == GretlFileType.unrecognized) {
            gretl_errmsg_set(_("Unknown data import type").ptr);
            err = GretlError.data;
        } else {
            err = gretl_read_foreign_data(fname, ft, dset, prn);
        }
    }

    return err;
}

/**
 * gretl_get_data:
 * @fname: name of file to try.
 * @dset: dataset struct.
 * @opt: option flags.
 * @prn: where any messages should be written.
 *
 * Wrapper for gretl_seek_data() designed for use by gretl4py,
 * respecting the const modifier for @fname.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int gretl_get_data(const(char) * fname, DatasetRaw * dset,
                   Gretlopt opt, Prn * prn)
{
    char[FILENAME_MAX] tmp = 0;

    if (fname == null || strlen(fname) >= tmp.sizeof) {
        return GretlError.data;
    }

    strcpy(tmp.ptr, fname);
    return gretl_seek_data(tmp.ptr, dset, opt, prn);
}

/**
 * open_nulldata:
 * @dset: dataset struct.
 * @data_status: indicator for whether a data file is currently open
 * in gretl's work space (1) or not (0).
 * @length: desired length of data series.
 * @opt: may contain OPT_N to suppress addition of an index series.
 * @prn: gretl printing struct.
 *
 * Create an empty "dummy" data set, suitable for simulations.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 *
 */

int open_nulldata(DatasetRaw * dset, int data_status, int length,
                  Gretlopt opt, Prn * prn)
{
    int t;

    /* clear any existing data info */
    if (data_status) {
        clear_datainfo(dset, DataClearCode.full);
    }

    /* dummy up the data info */
    dset.n = length;
    dset.v = (opt & Gretlopt.OPT_N)? 1 : 2;
    dataset_obs_info_default(dset);

    if (dataset_allocate_varnames(dset)) {
        return GretlError.alloc;
    }

    /* allocate data storage */
    if (allocate_Z(dset, Gretlopt.OPT_NONE)) {
        return GretlError.alloc;
    }

    if (dset.v > 1) {
        /* add an index var */
        strcpy(dset.varname[1], "index");
        series_set_label(dset, 1, _("index variable").ptr);
        for (t=0; t<dset.n; t++) {
            dset.Z[1][t] = cast(double) (t + 1);
        }
    }

    if (prn != null && gretl_messages_on()) {
        /* print basic info */
        pprintf(prn, _("periodicity: %d, maxobs: %d\n"
                       ~ "observations range: %s to %s\n").ptr,
                dset.pd, dset.n, dset.stobs.ptr, dset.endobs.ptr);
    }

    /* Set sample range to entire length of data-set by default */
    dset.t1 = 0;
    dset.t2 = dset.n - 1;

    set_current_dataset(dset);

    return 0;
}

private int extend_markers(DatasetRaw * dset, int old_n, int new_n)
{
    char ** S = cast(char **) realloc(dset.S, new_n * (char *).sizeof);
    int t, err = 0;

    if (S == null) {
        err = 1;
    } else {
        dset.S = S;
        for (t=old_n; t<new_n && !err; t++) {
            S[t] = cast(char *) malloc(OBSLEN);
            if (S[t] == null) {
                err = 1;
            }
        }
    }

    return err;
}

private void merge_error(const(char) * msg, Prn * prn)
{
    pputs(prn, msg);
    if (!printing_to_standard_stream(prn)) {
        /* hack to avoid revising translations */
        char * s = gretl_strndup(msg, strlen(msg) - 1);

        gretl_errmsg_set(s);
        g_free(s);
    }
}

private void merge_warning(const(char) * msg, Prn * prn)
{
    pprintf(prn, "%s\n", msg);
    if (!printing_to_standard_stream(prn)) {
        gretl_warnmsg_set(msg);
    }
}

private void merge_name_error(const(char) * objname, Prn * prn)
{
    char * msg;

    msg = g_strdup_printf(_("Can't replace %s with a series").ptr, objname);
    pprintf(prn, "%s\n", msg);
    if (!printing_to_standard_stream(prn)) {
        gretl_errmsg_set(msg);
    }
    g_free(msg);
}

private int count_new_vars(const(DatasetRaw) * d1, const(DatasetRaw) * d2,
                           Prn * prn)
{
    const(char) * vname;
    int addvars = d2.v - 1;
    int i, j;

    /* We start by assuming that all the series in @d2 are new,
       then subtract those we find to be already present. We also
       check for collision between the names of series to be added and
       the names of existing objects other than series.
    */

    for (i=1; i<d2.v && addvars >= 0; i++) {
        vname = d2.varname[i];
        if (gretl_is_user_var(vname)) {
            merge_name_error(vname, prn);
            addvars = -1;
        } else if (gretl_function_depth() > 0) {
            if (current_series_index(d1, vname) > 0) {
                addvars--;
            }
        } else {
            for (j=1; j<d1.v; j++) {
                if (!strcmp(vname, d1.varname[j])) {
                    addvars--;
                    break;
                }
            }
        }
    }

    return addvars;
}

private enum TimeSpecial {
    day = 1,
    yr = 2
}

private int year_special_markers(const(DatasetRaw) * ldset,
                                 const(DatasetRaw) * rdset)
{
    const(char) * test = null;
    int overlap = 0;
    int i, t, err = 0;

    /* See if we can match obs markers in @addset against years in
       @dset: we'll try this if all the markers in addset are integer
       strings, at least some of them are within the obs range of @dset,
       and none of them are outside of the "sanity" range of 1 to 2500.
    */

    if (!dataset_is_time_series(ldset) || ldset.pd != 1) {
        return 0;
    }

    if (ldset.markers || !rdset.markers) {
        return 0;
    }

    errno = 0;

    for (i=0; i<rdset.n; i++) {
        t = cast(int) strtol(rdset.S[i], &test, 10);
        if (*test || errno) {
            errno = 0;
            err = 1;
            break;
        }
        if (t < 1 || t > 2500) {
            err = 1;
            break;
        }
        if (!overlap) {
            t = dateton(rdset.S[i], ldset);
            if (t >= 0 && t < ldset.n) {
                overlap = 1;
            }
        }
    }

    return (!err && overlap) ? TimeSpecial.yr : 0;
}

private int compare_ranges(const(DatasetRaw) * ldset,
                           const(DatasetRaw) * rdset,
                           int newvars,
                           int * offset,
                           int * tspecial,
                           Prn * prn,
                           int * err)
{
    int ed0 = dateton(ldset.endobs.ptr, ldset);
    int sd1, ed1, addobs = -1;
    int range_err = 0;

    if (dataset_is_cross_section(ldset) &&
        dataset_is_cross_section(rdset) &&
        !(ldset.markers && rdset.markers)) {
        if (newvars == 0) {
            if (rdset.markers) {
                /* pass the problem on to just_append_rows */
                return 0;
            } else {
                /* assume the new data should be appended length-wise */
                *offset = ed0 + 1;
                return rdset.n;
            }
        } else {
            /* we've already determined that the series length in
               @rdset doesn't match either the full series length or
               the current sample range in @ldset; we therefore have
               no information with which to match rows for new
               series
            */
            pputs(prn, _("append: don't know how to align the new series!").ptr);
            pputc(prn, '\n');
            *err = GretlError.data;
            return -1;
        }
    }

    sd1 = merge_dateton(rdset.stobs.ptr, ldset);
    ed1 = merge_dateton(rdset.endobs.ptr, ldset);

    if (sd1 < 0) {
        /* case: new data start earlier than old */
        if (ed1 < 0) {
            range_err = 1;
        } else if (ed1 > ed0) {
            range_err = 2;
        } else {
            *offset = sd1;
            addobs = 0;
        }
    } else if (sd1 == 0 && ed1 == ed0) {
        /* case: exact match of ranges */
        *offset = 0;
        addobs = 0;
    } else if (sd1 == 0) {
        /* case: starting obs the same */
        *offset = 0;
        if (ed1 > ed0) {
            addobs = ed1 - ed0;
        } else {
            addobs = 0;
        }
    } else if (sd1 == ed0 + 1) {
        /* case: new data start right after end of old */
        *offset = sd1;
        addobs = rdset.n;
    } else if (sd1 > 0) {
        /* case: new data start later than old */
        if (sd1 <= ed0) {
            /* but there's some overlap */
            *offset = sd1;
            if (ed1 > ed0) {
                addobs = ed1 - ed0;
            } else {
                addobs = 0;
            }
        }
    }

    if (range_err) {
        /* try another approach? */
        *tspecial = year_special_markers(ldset, rdset);
        if (*tspecial) {
            addobs = 0;
        }
    }

    if (addobs < 0) {
        if (range_err == 1) {
            pputs(prn, "compare_ranges: no overlap, can't merge\n");
        } else if (range_err == 2) {
            pputs(prn, "compare ranges: new data start earlier, end later\n");
        } else {
            pputs(prn, "compare_ranges: unclassified error\n");
        }
    }

    return addobs;
}

/* Determine whether there's any overlap between the calendar
   in @addset and that in @dset. Return 0 on success (there is
   an overlap), non-zero otherwise.
*/

private int check_for_overlap(const(DatasetRaw) * dset,
                              const(DatasetRaw) * addset,
                              int * offset)
{
    int at1 = merge_dateton(addset.stobs.ptr, dset);
    int at2 = merge_dateton(addset.endobs.ptr, dset);

    if (!(at1 >= dset.n) && !(at2 < 0)) {
        /* OK, there must be some overlap */
        *offset = at1;
        return 0;
    } else {
        /* either the "add" data start after the original data end,
           or they end before the originals start, no there's no
           overlap
        */
        gretl_errmsg_set(_("No overlap in data ranges").ptr);
        return GretlError.data;
    }
}

/* When appending data to a current panel dataset, and the length of
   the series in the new data is less than the full panel size (N *
   T), try to determine if it's OK to expand the incoming data to
   match.

   We'll say it's OK if the new series length equals the panel T: in
   that case we'll take the new data to be time-series, which should
   be replicated for each panel unit -- unless the dataset represented
   by @addinfo is itself a panel.

   A second possibility arises if the length of the new series equals
   the panel N: in that case we could treat it as a time-invariant
   characteristic of the panel unit, which should be replicated for
   each time period.  But note that if OPT_T is given, this second
   expansion is forbidden: the user has stipulated that the new data
   are time-varying.
*/

private int panel_expand_ok(DatasetRaw * dset, DatasetRaw * addinfo,
                            Gretlopt opt)
{
    int N = dset.n / dset.pd;
    int T = dset.pd;
    int ok = 0;

    if (addinfo.n == T && !dataset_is_panel(addinfo)) {
        ok = 1;
    } else if (!(opt & Gretlopt.OPT_T) &&
               addinfo.n == N &&
               addinfo.pd == 1) {
        ok = 1;
    }

    return ok;
}

private int panel_append_special(DatasetRaw * dset,
                                 DatasetRaw * addset,
                                 int addvars,
                                 Gretlopt opt,
                                 Prn * prn)
{
    int n = dset.n / dset.pd;
    int T = dset.pd;
    int k = dset.v;
    int tsdata;
    int i, j, s, p, t;
    int err = 0;

    if (addvars > 0 && dataset_add_series(dset, addvars)) {
        merge_error(_("Out of memory!\n").ptr, prn);
        err = GretlError.alloc;
    }

    tsdata = ((opt & Gretlopt.OPT_T) || addset.n != n);

    for (i=1; i<addset.v && !err; i++) {
        int v = series_index(dset, addset.varname[i]);

        if (v >= k) {
            /* a new variable */
            v = k++;
            strcpy(dset.varname[v], addset.varname[i]);
            copy_varinfo(dset.varinfo[v], addset.varinfo[i]);
        }

        s = 0;
        for (j=0; j<n; j++) {
            /* loop across units */
            for (t=0; t<T; t++) {
                /* loop across periods */
                p = (tsdata)? t : j;
                dset.Z[v][s++] = addset.Z[i][p];
            }
        }
    }

    return err;
}

private int markers_compatible(const(DatasetRaw) * d1, DatasetRaw * d2,
                               int * offset)
{
    int ret = 0;

    if (d1.markers == 0 && d2.markers == 0) {
        *offset = d1.n;
        ret = 1;
    } else if (d1.markers == 0) {
        /* markers "on the right only": are they consecutive
           integers starting between 1 and d1->n + 1?
        */
        if (integer_string(d2.S[0])) {
            int k0 = atoi(d2.S[0]);

            if (k0 >= 1 && k0 <= d1.n + 1) {
                int i, k1;

                ret = 1;
                for (i=1; i<d2.n && ret; i++) {
                    if (!integer_string(d2.S[i])) {
                        ret = 0;
                    } else if ((k1 = atoi(d2.S[i])) != k0 + 1) {
                        ret = 0;
                    } else {
                        k0 = k1;
                    }
                }
            }
            if (ret) {
                *offset = atoi(d2.S[0]) - 1;
                /* the @d2 markers have done their job -- yielding
                   an @offset value -- and they can now be trashed
                */
                dataset_destroy_obs_markers(d2);
            }
        }
    } else {
        /* markers on both sides: are they totally distinct? */
        int i, j;

        ret = 1;
        for (i=0; i<d2.n && ret; i++) {
            for (j=0; j<d1.n && ret; j++) {
                if (!strcmp(d2.S[i], d1.S[j])) {
                    /* no, not totally distinct */
                    ret = 0;
                }
            }
        }
    }

    return ret;
}

private int just_append_rows(const(DatasetRaw) * targ, DatasetRaw * src, int * offset)
{
    int ret = 0;

    if (targ.structure == TsCodes.cross_section &&
        src.structure == TsCodes.cross_section &&
        targ.sd0 == 1 && src.sd0 == 1) {
        int ok, test_offset = -1;

        ok = markers_compatible(targ, src, &test_offset);
        if (ok) {
            /* note: we do this only if we're not adding any new
               series: we'll append to existing series lengthwise
               (or perhaps write data into existing rows)
            */
            *offset = test_offset;
            ret = src.n - (targ.n - *offset);
            if (ret < 0) {
                ret = 0;
            }
        }
    }

    return ret;
}

private int simple_range_match(const(DatasetRaw) * targ, const(DatasetRaw) * src,
                               int * offset)
{
    int ret = 0;

    if (src.pd == 1 && src.structure == TsCodes.cross_section) {
        if (src.n == targ.n) {
            ret = 1;
        } else if (src.n == targ.t2 - targ.t1 + 1) {
            ret = 1;
            *offset = targ.t1;
        }
    }

    return ret;
}

private int merge_series_data(DatasetRaw * ldset,
                              DatasetRaw * rdset,
                              int orig_v,
                              int orig_n,
                              int addvars,
                              int addobs,
                              int offset,
                              int update_overlap,
                              int tspecial)
{
    int k = orig_v;
    int newvar;
    int t, tmin;
    int i, v;
    int err = 0;

    for (i=1; i<rdset.v && !err; i++) {
        v = series_index(ldset, rdset.varname[i]);
        newvar = v >= k;
        if (!newvar && !update_overlap) {
            tmin = orig_n;
        } else {
            tmin = 0;
        }
        if (newvar) {
            v = k++;
            strcpy(ldset.varname[v], rdset.varname[i]);
            copy_varinfo(ldset.varinfo[v], rdset.varinfo[i]);
            if (is_string_valued(rdset, i) &&
                rdset.n == ldset.n && offset == 0 &&
                addobs == 0) {
                /* attach the string table to the target
                   series and detach it from @rdset
                */
                SeriesTable * st;

                st = series_get_string_table(rdset, i);
                series_attach_string_table(ldset, v, st);
                series_attach_string_table(rdset, i, null);
            }
        } else {
            /* not a new series */
            int lsval = is_string_valued(ldset, v);
            int rsval = is_string_valued(rdset, i);

            if (lsval + rsval == 1) {
                gretl_errmsg_set(_("Can't concatenate string-valued and numeric series").ptr);
                err = GretlError.data;
            } else if (lsval) {
                err = merge_string_tables(ldset, v, rdset, i);
            }
            if (err) {
                break;
            }
        }

        if (tspecial == TimeSpecial.day) {
            char[OBSLEN] obs;
            int s;

            for (t=tmin; t<ldset.n; t++) {
                ntolabel(obs.ptr, t, ldset);
                s = dateton(obs.ptr, rdset);
                if (s >= 0 && s < rdset.n) {
                    ldset.Z[v][t] = rdset.Z[i][s];
                } else {
                    ldset.Z[v][t] = NADBL;
                }
            }
        } else if (tspecial == TimeSpecial.yr) {
            int s;

            if (newvar) {
                for (t=0; t<ldset.n; t++) {
                    ldset.Z[v][t] = NADBL;
                }
            }
            for (s=0; s<rdset.n; s++) {
                t = dateton(rdset.S[s], ldset);
                if (t >= tmin && t < ldset.n) {
                    ldset.Z[v][t] = rdset.Z[i][s];
                }
            }
        } else {
            for (t=tmin; t<ldset.n; t++) {
                if (t >= offset && t - offset < rdset.n) {
                    ldset.Z[v][t] = rdset.Z[i][t - offset];
                } else if (newvar) {
                    ldset.Z[v][t] = NADBL;
                }
            }
        }
    }

    return err;
}

private int merge_lengthen_series(DatasetRaw * dset,
                                  const(DatasetRaw) * addset,
                                  int addobs,
                                  int offset)
{
    int i, t, new_n = dset.n + addobs;
    int err = 0;

    if (dset.markers) {
        err = extend_markers(dset, dset.n, new_n);
        if (!err) {
            for (t=dset.n; t<new_n; t++) {
                strcpy(dset.S[t], addset.S[t-offset]);
            }
        }
    }

    for (i=0; i<dset.v && !err; i++) {
        double * x;

        x = cast(double *) realloc(dset.Z[i], new_n * (*x).sizeof);
        if (x == null) {
            err = GretlError.alloc;
            break;
        }

        for (t=dset.n; t<new_n; t++) {
            if (i == 0) {
                x[t] = 1.0;
            } else {
                x[t] = NADBL;
            }
        }
        dset.Z[i] = x;
    }

    if (!err) {
        dset.n = new_n;
        ntolabel(dset.endobs.ptr, new_n - 1, dset);
        dset.t2 = dset.n - 1;
    }

    return err;
}

/* When we're merging dated weekly data, the day-of-week on which
   weeks are supposed to start should be consistent across the
   two datasets.
*/

private int check_week_start_dates(const(DatasetRaw) * ds0,
                                   const(DatasetRaw) * ds1)
{
    int wd0 = weekday_from_epoch_day(cast(uint) ds0.sd0);
    int wd1 = weekday_from_epoch_day(cast(uint) ds1.sd0);

    if (wd1 != wd0) {
        gretl_errmsg_set("Error appending data: inconsistent dating of weeks");
        return GretlError.data;
    } else {
        return 0;
    }
}

private int two_panels_check(DatasetRaw * dset, DatasetRaw * addset)
{
    if (dset.pd == addset.pd && dset.n == addset.n) {
        return 0;
    } else {
        gretl_errmsg_set(_("Panel data not conformable for appending").ptr);
        return GretlError.data;
    }
}

/**
 * merge_data:
 * @dset: dataset struct.
 * @addset: dataset to be merged in.
 * @opt: may include OPT_T to force a time-series interpretation
 * when appending to a panel dataset; may include OPT_U to update
 * values of overlapping observations.
 * @prn: print struct to accept messages.
 *
 * Attempt to merge the content of a newly opened data file into
 * gretl's current working data set.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

private int merge_data(DatasetRaw * dset, DatasetRaw * addset,
                       Gretlopt opt, Prn * prn)
{
    int update_overlap = (opt & Gretlopt.OPT_U);
    int orig_n = dset.n;
    int tspecial = 0;
    int fixsample = 0;
    int addsimple = 0;
    int addpanel = 0;
    int addvars = 0;
    int addobs = 0;
    int offset = 0;
    int err = 0;

    if (dataset_is_panel(dset) && dataset_is_panel(addset)) {
        err = two_panels_check(dset, addset);
        if (err) {
            return err;
        }
    }

    /* first see how many new vars we have */
    addvars = count_new_vars(dset, addset, prn);
    if (addvars < 0) {
        return 1;
    }

    if (dated_daily_data(dset) && dated_daily_data(addset)) {
        tspecial = TimeSpecial.day;
    } else if (dated_weekly_data(dset) && dated_weekly_data(addset)) {
        err = check_week_start_dates(dset, addset);
        if (err) {
            return err;
        }
    }

    if (opt & Gretlopt.OPT_X) {
        fixsample = 1;
    } else if (simple_range_match(dset, addset, &offset)) {
        /* we'll allow undated data to be merged with the existing
           dateset, sideways, provided the number of observations
           matches OK
        */
        addsimple = 1;
    } else if (dataset_is_panel(dset) &&
               panel_expand_ok(dset, addset, opt)) {
        /* allow appending to panel when the number of obs matches
           either the cross-section size or the time-series length
        */
        addpanel = 1;
    } else if (dset.pd != addset.pd) {
        merge_error(_("Data frequency does not match\n").ptr, prn);
        err = 1;
    }

    if (!err && fixsample) {
        err = check_for_overlap(dset, addset, &offset);
    } else if (!err && gretl_function_depth() > 0) {
        /* we won't add observations within a function, but
           we should still check for an error from compare_ranges()
        */
        if (!addsimple && !addpanel) {
            addobs = compare_ranges(dset, addset, addvars, &offset,
                                    &tspecial, prn, &err);
            if (!err && addobs > 0) {
                addobs = 0;
            }
        }
    } else if (!err) {
        if (!addsimple && !addpanel) {
            addobs = compare_ranges(dset, addset, addvars, &offset,
                                    &tspecial, prn, &err);
        }
        if (!err && addobs <= 0 && addvars == 0) {
            addobs = just_append_rows(dset, addset, &offset);
        }
    }

    if (!err && (addobs < 0 || addvars < 0)) {
        merge_error(_("New data not conformable for appending\n").ptr, prn);
        err = GretlError.data;
    } else if (!err && addobs == 0 && addvars == 0) {
        merge_warning(_("No data were found for appending").ptr, prn);
        return 0;
    }

    if (!err && !addpanel && dset.markers != addset.markers) {
        if (addobs == 0 && addvars == 0) {
            if (update_overlap) {
                ; /* might be OK? */
            } else {
                gretl_errmsg_set(_("Found no data conformable for appending").ptr);
                err = GretlError.data;
            }
        } else if (addset.n != dset.n && !tspecial) {
            merge_error(_("Inconsistency in observation markers\n").ptr, prn);
            err = GretlError.data;
        } else if (addset.markers && !dset.markers && !tspecial) {
            dataset_destroy_obs_markers(addset);
        }
    }

    /* if checks are passed, try merging the data */

    if (!err && addobs > 0) {
        err = merge_lengthen_series(dset, addset, addobs, offset);
        if (err) {
            merge_error(_("Out of memory!\n").ptr, prn);
        }
    }

    if (!err && addpanel) {
        err = panel_append_special(dset, addset, addvars, opt, prn);
    } else if (!err) {
        int orig_v = dset.v;

        if (addvars > 0 && dataset_add_series(dset, addvars)) {
            merge_error(_("Out of memory!\n").ptr, prn);
            err = GretlError.alloc;
        } else {
            err = merge_series_data(dset, addset, orig_v, orig_n,
                                    addvars, addobs, offset,
                                    update_overlap, tspecial);
        }
    }

    if (!err && (addvars || addobs) && gretl_messages_on()) {
        pputs(prn, _("Data appended OK\n").ptr);
    }

    return err;
}

/* We want to ensure that calendar dates are recorded as per
   ISO 8601 -- that is, YYYY-MM-DD; here we remedy dates
   recorded in the form YYYY/MM/DD.
*/

private void maybe_fix_calendar_dates(DatasetRaw * dset)
{
    if (strchr(dset.stobs.ptr, '/') != null) {
        gretl_charsub(dset.stobs.ptr, '/', '-');
        gretl_charsub(dset.endobs.ptr, '/', '-');
        if (dset.S != null && dset.markers == DatasetMarkerType.daily_date_strings) {
            int t;

            for (t=0; t<dset.n; t++) {
                gretl_charsub(dset.S[t], '/', '-');
            }
        }
    }
}

/**
 * get_merge_opts:
 * @opt: gretl options flags.
 *
 * Returns: just those components of @opt (if any) that
 * can be passed to merge_or_replace_data(); may be
 * useful when calling that function in the context
 * of a command only some of whose options should be
 * forwarded.
 */

Gretlopt get_merge_opts(Gretlopt opt)
{
    Gretlopt merge_opt = Gretlopt.OPT_NONE;

    if (opt & Gretlopt.OPT_T) {
        /* panel, common time-series */
        merge_opt |= Gretlopt.OPT_T;
    }
    if (opt & Gretlopt.OPT_U) {
        /* update overlapping observations */
        merge_opt |= Gretlopt.OPT_U;
    }
    if (opt & Gretlopt.OPT_X) {
        /* fixed sample range */
        merge_opt |= Gretlopt.OPT_X;
    }

    return merge_opt;
}

/* Apparatus for converting a dataset read from file
   into a gretl matrix, as opposed to replacing an
   existing dataset or merging with it.
*/

private __gshared GretlMatrixRaw ** dset_matrix;

void set_dset_matrix_target(GretlMatrixRaw ** pm)
{
    dset_matrix = pm;
}

void * get_dset_matrix_target()
{
    return dset_matrix;
}

/**
 * merge_or_replace_data:
 * @dset0: original dataset struct.
 * @pdset1: new dataset struct.
 * @opt: zero or more option flags (OPT_K presrves @pdset1,
 * otherwise it is destroyed).
 * @prn: print struct to accept messages.
 *
 * Given a newly-created dataset, pointed to by @pdset1, either
 * attempt to merge it with @dset0, if the original data array
 * is non-NULL, or replace the content of the original pointer
 * with the new dataset.
 *
 * In case merging is not successful, the new dataset is
 * destroyed.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int merge_or_replace_data(DatasetRaw * dset0, DatasetRaw ** pdset1,
                          Gretlopt opt, Prn * prn)
{
    int keep = (opt & Gretlopt.OPT_K);
    int err = 0;

    if (dset_matrix != null) {
        /* Convert the new dataset to matrix; don't touch
           the existing dataset, if any.
        */
        const(DatasetRaw) * dset = *pdset1;

        *dset_matrix = gretl_matrix_data_subset(null, dset,
                                                0, dset.n - 1,
                                                MMissingCode.ok,
                                                &err);
        destroy_dataset(*pdset1);
        *pdset1 = null;
        return err;
    }

    if (dset0.Z != null) {
        /* we have an existing dataset into which the new data
           should be merged */
        Gretlopt merge_opt = Gretlopt.OPT_NONE;

        if (opt & Gretlopt.OPT_T) {
            /* panel, common time-series */
            merge_opt |= Gretlopt.OPT_T;
        }
        if (opt & Gretlopt.OPT_U) {
            /* update overlapping observations */
            merge_opt |= Gretlopt.OPT_U;
        }
        if (opt & Gretlopt.OPT_X) {
            /* fixed sample range */
            merge_opt |= Gretlopt.OPT_X;
        }
        err = merge_data(dset0, *pdset1, merge_opt, prn);
        if (!keep) {
            destroy_dataset(*pdset1);
        }
    } else {
        /* starting from scratch */
        *dset0 = **pdset1;
        free(*pdset1);
        if (calendar_data(dset0)) {
            maybe_fix_calendar_dates(dset0);
        }
    }

    if (!keep) {
        *pdset1 = null;
    }

    set_current_dataset(dset0);

    return err;
}

private int check_imported_string(char * src, int i, size_t len)
{
    int err = 0;

    if (!g_utf8_validate(src, -1, null)) {
        char * trstr = null;
        size_t bytes;

        trstr = g_locale_to_utf8(src, -1, null, &bytes, null);

        if (trstr == null) {
            gretl_errmsg_sprintf(_("Invalid characters in imported string, line %d").ptr, i);
            err = GretlError.data;
        } else {
            *src = '\0';
            strncat(src, trstr, len - 1);
            g_free(trstr);
        }
    }

    return err;
}

private int count_markers(FILE * fp, char * line, int linelen,
                          char * marker)
{
    int n = 0;

    while (fgets(line, linelen, fp)) {
        if (sscanf(line, "%31[^\n\r]", marker) == 1) {
            g_strstrip(marker);
            if (*marker != '\0') {
                n++;
            }
        }
    }

    rewind(fp);

    return n;
}

private void finalize_add_markers(DatasetRaw * dset, char ** S, int err)
{
    if (err) {
        strings_array_free(S, dset.n);
    } else {
        if (dset.S != null) {
            strings_array_free(dset.S, dset.n);
        }
        dset.markers = cast(char) DatasetMarkerType.regular_markers;
        dset.S = S;
    }
}

/**
 * add_obs_markers_from_file:
 * @dset: data information struct.
 * @fname: name of file containing case markers.
 *
 * Read case markers (strings of %OBSLEN - 1 characters or less that identify
 * the observations) from a file, and associate them with the
 * current data set.  The file should contain one marker per line,
 * with a number of lines equal to the number of observations in
 * the current data set.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int add_obs_markers_from_file(DatasetRaw * dset, const(char) * fname)
{
    char ** S = null;
    FILE * fp;
    char[128] line;
    char[32] marker;
    int ns, panel_units = 0;
    int t, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == null) {
        return GretlError.fopen;
    }

    ns = count_markers(fp, line.ptr, cast(int) line.sizeof, marker.ptr);
    if (dataset_is_panel(dset) && ns == dset.n / dset.pd) {
        panel_units = 1; /* OK */
    } else if (ns < dset.n) {
        return GretlError.invarg;
    }

    S = strings_array_new_with_length(dset.n, OBSLEN);
    if (S == null) {
        fclose(fp);
        return GretlError.alloc;
    }

    if (panel_units) {
        int T = dset.pd;
        int t2, i = 0;

        while (fgets(line.ptr, cast(int) line.sizeof, fp) && !err) {
            marker[0] = '\0';
            if (sscanf(line.ptr, "%31[^\n\r]", marker.ptr) == 1) {
                g_strstrip(marker.ptr);
                strncat(S[i], marker.ptr, OBSLEN - 1);
                err = check_imported_string(S[i], i+1, OBSLEN);
                if (!err) {
                    /* copy to remaining observations */
                    for (t2=1; t2<T; t2++) {
                        strcpy(S[i+t2], S[i]);
                    }
                }
                i += T;
            }
        }
    } else {
        for (t=0; t<dset.n && !err; t++) {
            if (fgets(line.ptr, cast(int) line.sizeof, fp) == null) {
                gretl_errmsg_sprintf(_("Expected %d markers; found %d\n").ptr,
                                     dset.n, t);
                err = GretlError.data;
            } else if (sscanf(line.ptr, "%31[^\n\r]", marker.ptr) != 1) {
                gretl_errmsg_sprintf(_("Couldn't read marker on line %d").ptr, t+1);
                err = GretlError.data;
            } else {
                g_strstrip(marker.ptr);
                strncat(S[t], marker.ptr, OBSLEN - 1);
                err = check_imported_string(S[t], t+1, OBSLEN);
            }
        }
    }

    finalize_add_markers(dset, S, err);

    return err;
}

private int add_obs_markers_from_array(DatasetRaw * dset, const(char) * aname)
{
    GretlArray * a;
    char ** S = null;
    const(char) * s;
    const(char) * white = " \t\r\n";
    char[32] marker;
    int ns, panel_units = 0;
    int t, err = 0;

    a = get_strings_array_by_name(aname);
    if (a == null) {
        return GretlError.invarg;
    }

    ns = gretl_array_get_length(a);
    if (dataset_is_panel(dset) && ns == dset.n / dset.pd) {
        panel_units = 1; /* OK */
    } else if (ns < dset.n) {
        return GretlError.invarg;
    }

    S = strings_array_new_with_length(dset.n, OBSLEN);
    if (S == null) {
        return GretlError.alloc;
    }

    if (panel_units) {
        int N = dset.n / dset.pd;
        int T = dset.pd;
        int i, j = 0;

        for (i=0; i<N; i++) {
            s = cast(const(char) *) gretl_array_get_data(a, i);
            s += strspn(s, white);
            marker[0] = '\0';
            strncat(marker.ptr, s, 31);
            g_strstrip(marker.ptr);
            gretl_utf8_truncate_b(marker.ptr, OBSLEN-1);
            strcpy(S[j++], marker.ptr);
            for (t=1; t<T; t++) {
                strcpy(S[j++], marker.ptr);
            }
        }
    } else {
        for (t=0; t<dset.n; t++) {
            s = cast(const(char) *) gretl_array_get_data(a, t);
            s += strspn(s, white);
            marker[0] = '\0';
            strncat(marker.ptr, s, 31);
            g_strstrip(marker.ptr);
            gretl_utf8_truncate_b(marker.ptr, OBSLEN-1);
            strcpy(S[t], marker.ptr);
        }
    }

    finalize_add_markers(dset, S, err);

    return err;
}

private int add_obs_markers_from_series(DatasetRaw * dset, const(char) * vname)
{
    char ** S = null;
    const(char) * s;
    const(char) * white = " \t\r\n";
    char[32] marker;
    int v, t, err = 0;

    v = current_series_index(dset, vname);
    if (v <= 0) {
        return GretlError.unkvar;
    } else if (!is_string_valued(dset, v)) {
        return GretlError.invarg;
    }

    S = strings_array_new_with_length(dset.n, OBSLEN);
    if (S == null) {
        return GretlError.alloc;
    }

    for (t=0; t<dset.n; t++) {
        s = series_get_string_for_obs(dset, v, t);
        if (s == null) {
            gretl_errmsg_sprintf(_("Missing marker at observation %d").ptr, t+1);
            err = GretlError.data;
            break;
        } else {
            s += strspn(s, white);
            marker[0] = '\0';
            strncat(marker.ptr, s, 31);
            g_strstrip(marker.ptr);
            gretl_utf8_truncate_b(marker.ptr, OBSLEN-1);
            strcpy(S[t], marker.ptr);
        }
    }

    finalize_add_markers(dset, S, err);

    return err;
}

/**
 * dataset_has_var_labels:
 * @dset: data information struct.
 *
 * Returns: 1 if at least one variable in the current dataset
 * has a descriptive label, otherwise 0.
 */

int dataset_has_var_labels(const(DatasetRaw) * dset)
{
    const(char) * vlabel;
    int i, imin = 1;

    if (dset.v > 1) {
        if (!strcmp(dset.varname[1], "index")) {
            vlabel = series_get_label(dset, 1);
            if (vlabel != null && !strcmp(vlabel, _("index variable").ptr)) {
                imin = 2;
            }
        }
    }

    for (i=imin; i<dset.v; i++) {
        vlabel = series_get_label(dset, i);
        if (vlabel != null && *vlabel != '\0') {
            return 1;
        }
    }

    return 0;
}

/**
 * save_var_labels_to_file:
 * @dset: data information struct.
 * @fname: name of file containing labels.
 *
 * Writes to @fname the descriptive labels for the series in
 * the current dataset.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int save_var_labels_to_file(const(DatasetRaw) * dset,
                            const(char) * fname)
{
    const(char) * vlabel;
    FILE * fp;
    int i, err = 0;

    fp = gretl_fopen(fname, "w");

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        for (i=1; i<dset.v; i++) {
            vlabel = series_get_label(dset, i);
            fprintf(fp, "%s\n", vlabel == null ? "".ptr : vlabel);
        }
        fclose(fp);
    }

    return err;
}

private int save_dset_description_to_file(const(DatasetRaw) * dset,
                                          const(char) * fname)
{
    FILE * fp;
    int err = 0;

    fp = gretl_fopen(fname, "w");

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        fputs(dset.descrip, fp);
        fclose(fp);
    }

    return err;
}

private int save_var_labels_to_array(const(DatasetRaw) * dset,
                                     const(char) * aname)
{
    GretlArray * a = null;
    int err = 0;

    if (gretl_is_series(aname, dset)) {
        err = GretlError.types;
    } else {
        err = check_identifier(aname);
    }

    if (!err) {
        a = gretl_array_new(GretlType.strings, dset.v - 1, &err);
    }

    if (!err) {
        err = user_var_add_or_replace(aname, GretlType.strings, a);
    }

    if (!err) {
        char * vlabel;
        int i;

        for (i=1; i<dset.v; i++) {
            vlabel = cast(char *) series_get_label(dset, i);
            gretl_array_set_element(a, i-1, vlabel != null ? vlabel : cast(char *) "".ptr,
                                    GretlType.string_, 1);
        }
    }

    if (err && a != null) {
        gretl_array_destroy(a);
        a = null;
    }

    return err;
}

private int save_obs_markers_to_array(const(DatasetRaw) * dset,
                                      const(char) * aname)
{
    GretlArray * a = null;
    int err = 0;

    if (gretl_is_series(aname, dset)) {
        err = GretlError.types;
    } else {
        err = check_identifier(aname);
    }

    if (!err) {
        a = gretl_array_new(GretlType.strings, dset.n, &err);
    }

    if (!err) {
        err = user_var_add_or_replace(aname, GretlType.strings, a);
    }

    if (!err) {
        char * marker;
        int i;

        for (i=0; i<dset.n; i++) {
            marker = cast(char *) dset.S[i];
            gretl_array_set_element(a, i, marker != null ? marker : cast(char *) "".ptr,
                                    GretlType.string_, 1);
        }
    }

    if (err && a != null) {
        gretl_array_destroy(a);
        a = null;
    }

    return err;
}

/**
 * add_var_labels_from_file:
 * @dset: data information struct.
 * @fname: name of file containing labels.
 *
 * Read descriptive variables for labels (strings of %MAXLABEL - 1
 * characters or less) from a file, and associate them with the
 * current data set.  The file should contain one label per line,
 * with a number of lines equal to the number of variables in
 * the current data set, excluding the constant.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int add_var_labels_from_file(DatasetRaw * dset, const(char) * fname)
{
    FILE * fp;
    char[1024] line;
    char * label;
    int nlabels = 0;
    int i, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == null) {
        return GretlError.fopen;
    }

    for (i=1; i<dset.v && !err; i++) {
        if (fgets(line.ptr, cast(int) line.sizeof, fp) == null) {
            break;
        } else {
            label = g_strstrip(g_strdup(line.ptr));
            if (strlen(label) > 0) {
                if (!g_utf8_validate(label, -1, null)) {
                    char * trstr = null;
                    size_t bytes;

                    trstr = g_locale_to_utf8(label, -1, null,
                                             &bytes, null);
                    if (trstr != null) {
                        series_set_label(dset, i, trstr);
                        nlabels++;
                        g_free(trstr);
                    }
                } else {
                    series_set_label(dset, i, label);
                    nlabels++;
                }
            }
            g_free(label);
        }
    }

    if (!err && nlabels == 0) {
        gretl_errmsg_set(_("No labels found").ptr);
        err = GretlError.data;
    }

    return err;
}

private int add_dset_description_from_file(DatasetRaw * dset,
                                           const(char) * fname)
{
    char * buf = null;
    GError * gerr = null;
    size_t size;
    int err = 0;

    g_file_get_contents(fname, &buf, &size, &gerr);

    if (gerr != null) {
        gretl_errmsg_set(gerr.message);
        g_error_free(gerr);
        return GretlError.data;
    }

    if (size == 0) {
        gretl_errmsg_set("No content found");
        err = GretlError.data;
    } else if (g_utf8_validate(buf, -1, null)) {
        free(dset.descrip);
        dset.descrip = gretl_strdup(buf);
    } else {
        char * trstr = null;

        trstr = g_locale_to_utf8(buf, -1, null,
                                 &size, &gerr);
        if (gerr != null) {
            gretl_errmsg_set(gerr.message);
            g_error_free(gerr);
            err = GretlError.data;
        } else {
            free(dset.descrip);
            dset.descrip = gretl_strdup(trstr);
            g_free(trstr);
        }
    }

    g_free(buf);

    return err;
}

private int add_var_labels_from_array(DatasetRaw * dset, const(char) * aname)
{
    GretlArray * a = get_array_by_name(aname);
    int i, err = 0;

    if (a == null) {
        gretl_errmsg_sprintf(_("%s: no such array").ptr, aname);
        err = GretlError.data;
    } else if (gretl_array_get_type(a) != GretlType.strings ||
               gretl_array_get_length(a) < dset.v - 1) {
        err = GretlError.types;
    }

    for (i=1; i<dset.v && !err; i++) {
        const(char) * s = cast(const(char) *) gretl_array_get_data(a, i-1);

        series_set_label(dset, i, s);
    }

    return err;
}

int read_or_write_var_labels(Gretlopt opt, DatasetRaw * dset, Prn * prn)
{
    const(char) * lname = null;
    int err;

    err = incompatible_options(opt, Gretlopt.OPT_D | Gretlopt.OPT_T | Gretlopt.OPT_F |
                               Gretlopt.OPT_A | Gretlopt.OPT_R);
    if (err) {
        return err;
    }

    if (opt & (Gretlopt.OPT_T | Gretlopt.OPT_F | Gretlopt.OPT_A | Gretlopt.OPT_R)) {
        lname = get_optval_string(GretlCmdIndex.labels, opt);
        if (lname == null) {
            return GretlError.badopt;
        } else if (opt & (Gretlopt.OPT_T | Gretlopt.OPT_F)) {
            gretl_maybe_switch_dir(lname);
        }
    }

    if (opt & Gretlopt.OPT_D) {
        /* delete */
        int i;

        for (i=1; i<dset.v; i++) {
            series_set_label(dset, i, "");
        }
    } else if (opt & (Gretlopt.OPT_T | Gretlopt.OPT_R)) {
        /* to-file, to-array */
        if (!dataset_has_var_labels(dset)) {
            pprintf(prn, _("No labels are available for writing\n").ptr);
            err = GretlError.data;
        } else {
            if (opt & Gretlopt.OPT_T) {
                err = save_var_labels_to_file(dset, lname);
            } else {
                err = save_var_labels_to_array(dset, lname);
            }
            if (!err && gretl_messages_on()) {
                pprintf(prn, _("Labels written OK\n").ptr);
            }
        }
    } else if (opt & (Gretlopt.OPT_F | Gretlopt.OPT_A)) {
        /* from-file, from-array */
        if (opt & Gretlopt.OPT_F) {
            err = add_var_labels_from_file(dset, lname);
        } else {
            err = add_var_labels_from_array(dset, lname);
        }
        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Labels loaded OK\n").ptr);
        }
    }

    return err;
}

int read_or_write_dset_description(Gretlopt opt, DatasetRaw * dset, Prn * prn)
{
    const(char) * fname = null;
    int err;

    if (dset == null) {
        return GretlError.nodata;
    }

    err = incompatible_options(opt, Gretlopt.OPT_F | Gretlopt.OPT_T);
    if (err) {
        return err;
    }

    fname = get_optval_string(GretlCmdIndex.info, opt);
    if (fname == null) {
        err = GretlError.badopt;
    }

    gretl_maybe_switch_dir(fname);

    if (opt & Gretlopt.OPT_T) {
        /* to-file */
        if (dset.descrip == null) {
            pputs(prn, _("No description is available\n").ptr);
            err = GretlError.data;
        } else {
            err = save_dset_description_to_file(dset, fname);
            if (!err && gretl_messages_on()) {
                pputs(prn, _("dataset description written OK\n").ptr);
            }
        }
    } else if (opt & Gretlopt.OPT_F) {
        /* from-file */
        int have_desc = dset.descrip != null;

        err = add_dset_description_from_file(dset, fname);
        if (!err && gretl_messages_on()) {
            if (have_desc) {
                pputs(prn, _("description replaced OK\n").ptr);
            } else {
                pputs(prn, _("description loaded OK\n").ptr);
            }
        }
    }

    return err;
}

private int save_obs_markers_to_file(DatasetRaw * dset, const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "w");
    int err = 0;

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        int i;

        for (i=0; i<dset.n; i++) {
            fprintf(fp, "%s\n", dset.S[i]);
        }
        fclose(fp);
    }

    return err;
}

int read_or_write_obs_markers(Gretlopt opt, DatasetRaw * dset, Prn * prn)
{
    const(char) * fname = null;
    int err;

    err = incompatible_options(opt, Gretlopt.OPT_D | Gretlopt.OPT_T | Gretlopt.OPT_F |
                               Gretlopt.OPT_A | Gretlopt.OPT_R);
    if (err) {
        return err;
    }

    if (opt & (Gretlopt.OPT_T | Gretlopt.OPT_F)) {
        fname = get_optval_string(GretlCmdIndex.markers, opt);
        if (fname == null) {
            return GretlError.badopt;
        } else {
            fname = gretl_maybe_switch_dir(fname);
        }
    }

    if (opt & (Gretlopt.OPT_A | Gretlopt.OPT_T)) {
        /* writing to file or array */
        if (dset.S == null) {
            gretl_errmsg_set(_("No markers are available for writing").ptr);
            return GretlError.data;
        }
    }

    if (opt & Gretlopt.OPT_D) {
        /* delete */
        dataset_destroy_obs_markers(dset);
    } else if (opt & Gretlopt.OPT_T) {
        /* to-file */
        err = save_obs_markers_to_file(dset, fname);
        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Markers written OK\n").ptr);
        }
    } else if (opt & Gretlopt.OPT_F) {
        /* from-file */
        err = add_obs_markers_from_file(dset, fname);
        if (!err && gretl_messages_on()) {
            pprintf(prn, _("Markers loaded OK\n").ptr);
        }
    } else if (opt & Gretlopt.OPT_A) {
        /* to-array */
        const(char) * aname = get_optval_string(GretlCmdIndex.markers, Gretlopt.OPT_A);

        err = save_obs_markers_to_array(dset, aname);
    } else if (opt & Gretlopt.OPT_R) {
        /* from-array */
        const(char) * aname = get_optval_string(GretlCmdIndex.markers, Gretlopt.OPT_R);

        err = add_obs_markers_from_array(dset, aname);
    } else if (opt & Gretlopt.OPT_S) {
        /* from-series */
        const(char) * vname = get_optval_string(GretlCmdIndex.markers, Gretlopt.OPT_S);

        err = add_obs_markers_from_series(dset, vname);
    }

    return err;
}

private void octave_varname(char * name, const(char) * s, int nnum, int v)
{
    char[12] nstr = 0;
    int len, tr;

    if (nnum == 0) {
        strcpy(name, s);
    } else {
        sprintf(nstr.ptr, "%d", nnum);
        len = cast(int) strlen(nstr.ptr);
        tr = VNAMELEN - len - 1; /* leave room for the NUL terminator */

        if (tr > 0) {
            strncat(name, s, tr);
            strcat(name, nstr.ptr);
        } else {
            sprintf(name, "v%d", v);
        }
    }
}


private int get_max_line_length(FILE * fp, Prn * prn)
{
    int c, c1, cc = 0;
    int maxlen = 0;

    while ((c = fgetc(fp)) != EOF) {
        if (c == 0x0d) {
            /* CR */
            c1 = fgetc(fp);
            if (c1 == EOF) {
                break;
            } else if (c1 == 0x0a) {
                /* CR + LF -> LF */
                c = c1;
            } else {
                /* Mac-style: CR not followed by LF */
                c = 0x0a;
                ungetc(c1, fp);
            }
        }
        if (c == 0x0a) {
            if (cc > maxlen) {
                maxlen = cc;
            }
            cc = 0;
            continue;
        }
        if (!isspace(cast(ubyte) c) && !isprint(cast(ubyte) c) &&
            !(c == CTRLZ)) {
            pprintf(prn, _("Binary data (%d) encountered: this is not a valid "
                           ~ "text file\n").ptr, c);
            return -1;
        }
        cc++;
    }

    if (maxlen == 0) {
        pprintf(prn, _("Data file is empty\n").ptr);
    }

    if (maxlen > 0) {
        /* allow for newline and null terminator */
        maxlen += 3;
    }

    return maxlen;
}

private int import_octave(const(char) * fname, DatasetRaw * dset,
                          Gretlopt opt, Prn * prn)
{
    DatasetRaw * octset = null;
    FILE * fp = null;
    char * line = null;
    char[8] tmp = 0;
    char[16] fmt = 0;
    char[32] name = 0;
    int nrows = 0, ncols = 0, nblocks = 0;
    int brows = 0, bcols = 0, oldbcols = 0;
    int maxlen, got_type = 0, got_name = 0;
    int i, t, err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == null) {
        return GretlError.fopen;
    }

    pprintf(prn, "%s %s...\n", _("parsing").ptr, fname);

    maxlen = get_max_line_length(fp, prn);
    if (maxlen <= 0) {
        err = GretlError.data;
        goto oct_bailout;
    }

    line = cast(char *) malloc(maxlen);
    if (line == null) {
        err = GretlError.alloc;
        goto oct_bailout;
    }

    pprintf(prn, _("   longest line: %d characters\n").ptr, maxlen - 1);

    rewind(fp);

    while (fgets(line, maxlen, fp) && !err) {
        if (*line == '#') {
            if (!got_name) {
                if (sscanf(line, "# name: %31s", name.ptr) == 1) {
                    got_name = 1;
                    nblocks++;
                    continue;
                }
            }
            if (!got_type) {
                if (sscanf(line, "# type: %7s", tmp.ptr) == 1) {
                    if (!got_name || strcmp(tmp.ptr, "matrix")) {
                        err = 1;
                    } else {
                        got_type = 1;
                    }
                    continue;
                }
            }
            if (brows == 0) {
                if (sscanf(line, "# rows: %d", &brows) == 1) {
                    if (!got_name || !got_type || brows <= 0) {
                        err = 1;
                    } else if (nrows > 0 && brows != nrows) {
                        err = 1;
                    } else {
                        nrows = brows;
                    }
                    continue;
                }
            }
            if (bcols == 0) {
                if (sscanf(line, "# columns: %d", &bcols) == 1) {
                    if (!got_name || !got_type || bcols <= 0) {
                        err = 1;
                    } else {
                        ncols += bcols;
                        pprintf(prn, _("   Found matrix '%s' with "
                                       ~ "%d rows, %d columns\n").ptr, name.ptr, brows, bcols);
                    }
                    continue;
                }
            }
        } else if (string_is_blank(line)) {
            continue;
        } else {
            got_name = 0;
            got_type = 0;
            brows = 0;
            bcols = 0;
        }
    }

    if (err || nrows == 0 || ncols == 0) {
        pputs(prn, _("Invalid data file\n").ptr);
        err = GretlError.data;
        goto oct_bailout;
    }

    /* initialize datainfo and Z */

    octset = datainfo_new();
    if (octset == null) {
        pputs(prn, _("Out of memory!\n").ptr);
        err = GretlError.alloc;
        goto oct_bailout;
    }

    octset.n = nrows;
    octset.v = ncols + 1;

    if (start_new_Z(octset, cast(Gretlopt) 0)) {
        pputs(prn, _("Out of memory!\n").ptr);
        err = GretlError.alloc;
        goto oct_bailout;
    }

    rewind(fp);

    pprintf(prn, _("   number of variables: %d\n").ptr, ncols);
    pprintf(prn, _("   number of observations: %d\n").ptr, nrows);
    pprintf(prn, _("   number of data blocks: %d\n").ptr, nblocks);

    i = 1;
    t = 0;

    sprintf(fmt.ptr, "# name: %%%ds", VNAMELEN - 1);

    while (fgets(line, maxlen, fp) && !err) {
        char * s = line;
        int j;

        if (*s == '#') {
            if (sscanf(line, fmt.ptr, name.ptr) == 1) {
                ;
            } else if (sscanf(line, "# rows: %d", &brows) == 1) {
                t = 0;
            } else if (sscanf(line, "# columns: %d", &bcols) == 1) {
                i += oldbcols;
                oldbcols = bcols;
            }
        }

        if (*s == '#' || string_is_blank(s)) {
            continue;
        }

        if (t >= octset.n) {
            err = 1;
        }

        for (j=0; j<bcols && !err; j++) {
            double x;
            int v = i + j;

            if (t == 0) {
                int nnum = (bcols > 1)? j + 1 : 0;

                octave_varname(octset.varname[i+j], name.ptr, nnum, v);
            }

            while (isspace(*s)) s++;
            if (sscanf(s, "%lf", &x) != 1) {
                fprintf(stderr, "error: '%s', didn't get double\n", s);
                err = 1;
            } else {
                octset.Z[v][t] = x;
                while (!isspace(*s)) s++;
            }
        }
        t++;
    }

    if (err) {
        pputs(prn, _("Invalid data file\n").ptr);
        err = GretlError.data;
    } else {
        int merge = dset.Z != null;
        Gretlopt merge_opt = cast(Gretlopt) 0;

        if (merge && (opt & Gretlopt.OPT_T)) {
            merge_opt = Gretlopt.OPT_T;
        }
        err = merge_or_replace_data(dset, &octset, merge_opt, prn);
    }

 oct_bailout:

    if (fp != null) {
        fclose(fp);
    }

    if (line != null) {
        free(line);
    }

    if (octset != null) {
        clear_datainfo(octset, DataClearCode.full);
    }

    return err;
}

/**
 * import_other:
 * @fname: name of file.
 * @ftype: type of data file.
 * @dset: pointer to dataset struct.
 * @opt: option flag; see gretl_get_data().
 * @prn: gretl printing struct.
 *
 * Open a data file of a type that requires a special plugin.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_other(const(char) * fname, GretlFileType ftype,
                 DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    extern (C) int function(const(char) *, DatasetRaw *,
                            Gretlopt, Prn *) importer;
    int err = 0;

    if (gretl_test_fopen(fname, "r") != 0) {
        pprintf(prn, _("Couldn't open %s\n").ptr, fname);
        return GretlError.fopen;
    }

    if (ftype == GretlFileType.octave) {
        /* plugin not needed */
        return import_octave(fname, dset, opt, prn);
    }

    if (ftype == GretlFileType.wf1) {
        importer = cast(typeof(importer)) get_plugin_function("wf1_get_data");
    } else if (ftype == GretlFileType.dta) {
        importer = cast(typeof(importer)) get_plugin_function("dta_get_data");
    } else if (ftype == GretlFileType.sav) {
        importer = cast(typeof(importer)) get_plugin_function("sav_get_data");
    } else if (ftype == GretlFileType.sas) {
        importer = cast(typeof(importer)) get_plugin_function("xport_get_data");
    } else if (ftype == GretlFileType.jmulti) {
        importer = cast(typeof(importer)) get_plugin_function("jmulti_get_data");
    } else if (ftype == GretlFileType.map) {
        importer = cast(typeof(importer)) get_plugin_function("map_get_data");
    } else {
        pprintf(prn, _("Unrecognized data type").ptr);
        pputc(prn, '\n');
        return GretlError.data;
    }

    if (importer == null) {
        err = 1;
    } else {
        err = importer(fname, dset, opt, prn);
    }

    return err;
}

/**
 * import_spreadsheet:
 * @fname: name of file.
 * @ftype: type of data file.
 * @list: list of parameters for spreadsheet import, or NULL.
 * @sheetname: name of specific worksheet, or NULL.
 * @dset: dataset struct.
 * @opt: option flag; see gretl_get_data().
 * @prn: gretl printing struct.
 *
 * Open a data file of a type that requires a special plugin.
 * Acceptable values for @ftype are %GRETL_GNUMERIC,
 * %GRETL_XLS, %GRETL_XLSX and %GRETL_ODS.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_spreadsheet(const(char) * fname, GretlFileType ftype,
                       int * list, char * sheetname,
                       DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    extern (C) int function(const(char) *, int *, char *,
                            DatasetRaw *, Gretlopt, Prn *,
                            void *) importer;
    int err = 0;

    if (gretl_test_fopen(fname, "r") != 0) {
        pprintf(prn, _("Couldn't open %s\n").ptr, fname);
        return GretlError.fopen;
    }

    import_na_init();

    if (ftype == GretlFileType.gnumeric) {
        importer = cast(typeof(importer)) get_plugin_function("gnumeric_get_data");
    } else if (ftype == GretlFileType.xls) {
        importer = cast(typeof(importer)) get_plugin_function("xls_get_data");
    } else if (ftype == GretlFileType.xlsx) {
        importer = cast(typeof(importer)) get_plugin_function("xlsx_get_data");
    } else if (ftype == GretlFileType.ods) {
        importer = cast(typeof(importer)) get_plugin_function("ods_get_data");
    } else {
        pprintf(prn, _("Unrecognized data type").ptr);
        pputc(prn, '\n');
        return GretlError.data;
    }

    if (importer == null) {
        fprintf(stderr, "import_spreadsheet: importer not found\n");
        err = 1;
    } else {
        gchar * thisdir = g_get_current_dir();

        err = importer(fname, list, sheetname, dset, opt, prn, null);
        if (thisdir != null) {
            /* come back out of dotdir? */
            gretl_chdir(thisdir);
            g_free(thisdir);
        }
    }

    return err;
}

private int is_jmulti_datafile(const(char) * fname)
{
    FILE * fp;
    int ret = 0;

    fp = gretl_fopen(fname, "r");

    if (fp != null) {
        char[128] test = 0;
        int gotobs = 0;
        int gotcomm = 0;
        int incomm = 0;

        /* look for characteristic C-style comment and
           <obs stuff> field, outside of comment */

        while (fgets(test.ptr, cast(int) test.sizeof, fp)) {
            if (!incomm && strstr(test.ptr, "/*")) {
                gotcomm = 1;
                incomm = 1;
            }
            if (incomm && strstr(test.ptr, "*/")) {
                incomm = 0;
            }
            if (!incomm && test[0] == '<' && strchr(test.ptr, '>')) {
                gotobs = 1;
            }
            if (gotcomm && gotobs) {
                ret = 1;
                break;
            }
        }
        fclose(fp);
    }

    return ret;
}

/**
 * gretl_is_pkzip_file:
 * @fname: name of file to examine.
 *
 * Returns: 1 if @fname is readable and is a PKZIP file,
 * else 0.
 */

int gretl_is_pkzip_file(const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "rb");
    int ret = 0;

    if (fp != null) {
        char[3] test = 0;

        if (fread(test.ptr, 1, 2, fp) == 2 && !strcmp(test.ptr, "PK")) {
            ret = 1;
        }
        fclose(fp);
    }

    return ret;
}

/**
 * detect_filetype:
 * @fname: the name of the file to test.
 * @opt: OPT_P may be included to permit path-searching if @fname
 * is not an absolute path; in that case the @fname argument
 * may be modified, otherwise it will be left unchanged.
 *
 * Attempts to determine the type of a file to be opened in gretl:
 * data file (of various formats), or command script. If OPT_P
 * is given, the @fname argument must be an array of length
 * at least FILENAME_MAX (libgretl.h): a path may be prepended
 * and in some cases an extension may be appended.
 *
 * Returns: integer code indicating the type of file.
 */

GretlFileType detect_filetype(char * fname, Gretlopt opt)
{
    const(char) * ext = get_filename_extension(fname);
    GretlFileType ftype = GretlFileType.unrecognized;

    if (ext != null) {
        /* First try judging the type by extension */
        if (!strcmp(ext, ".inp")) {
            ftype = GretlFileType.script;
        } else if (!strcmp(ext, ".gretl")) {
            if (gretl_is_pkzip_file(fname)) {
                ftype = GretlFileType.session;
            } else {
                ftype = GretlFileType.script;
            }
        } else {
            ftype = data_file_type_from_extension(ext, fname);
            if (ftype == GretlFileType.unrecognized) {
                /* check for database types */
                if (!strcmp(ext, ".bin")) {
                    ftype = GretlFileType.native_db;
                } else if (!strcmp(ext, ".rat")) {
                    ftype = GretlFileType.rats_db;
                } else if (!strcmp(ext, ".bn7")) {
                    ftype = GretlFileType.pcgive_db;
                }
            }
        }
        if (ftype != GretlFileType.unrecognized) {
            /* We got a type from the extension, but can we find
               the file "as is"? If so, we're done.
            */
            if (gretl_test_fopen(fname, "r") == 0) {
                return ftype;
            }
        }
    }

    if ((opt & Gretlopt.OPT_P) && gretl_addpath(fname, 0) != null) {
        ext = get_filename_extension(fname);
        if (ext != null) {
            /* check again for known data file types */
            ftype = data_file_type_from_extension(ext, fname);
        }
    }

    if (ftype == GretlFileType.unrecognized) {
        /* last gasp */
        if (gretl_is_xml_file(fname)) {
            ftype = GretlFileType.xml_data;
        } else if (has_suffix(fname, ".dat") && is_jmulti_datafile(fname)) {
            ftype = GretlFileType.jmulti;
        } else {
            /* default to assuming plain text data */
            ftype = GretlFileType.csv;
        }
    }

    return ftype;
}

/**
 * check_atof:
 * @numstr: string to check.
 *
 * Returns: 0 if @numstr is blank, or is a valid string representation
 * of a floating point number, else non-zero. What counts as valid
 * input depends on the current locale decimal character.
 */

int check_atof(const(char) * numstr)
{
    int err = 0;
    const(char) * test;

    /* accept blank entries */
    if (*numstr == '\0') {
        return 0;
    }

    errno = 0;
    strtod(numstr, &test);

    if (*test != '\0') {
        if (strcmp(numstr, test) == 0) {
            gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!").ptr, numstr);
        } else if (isprint(*test)) {
            gretl_errmsg_sprintf(_("Extraneous character '%c' in data").ptr, *test);
        } else {
            gretl_errmsg_sprintf(_("Extraneous character (0x%x) in data").ptr, *test);
        }
        err = 1;
    } else if (errno == ERANGE) {
        gretl_errmsg_sprintf(_("'%s' -- number out of range!").ptr, numstr);
        err = 1;
    }

    return err;
}

/**
 * check_atoi:
 * @numstr: string to check.
 *
 * Returns: 0 if @numstr is blank, or is a valid string representation
 * of a 32-bit signed integer, else non-zero.
 */

int check_atoi(const(char) * numstr)
{
    int err = 0;
    c_long val;
    const(char) * test;

    /* accept blank entries */
    if (*numstr == '\0') {
        return 0;
    }

    errno = 0;
    val = strtol(numstr, &test, 10);

    if (*test != '\0') {
        if (strcmp(numstr, test) == 0) {
            gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!").ptr, numstr);
        } else if (isprint(*test)) {
            gretl_errmsg_sprintf(_("Extraneous character '%c' in data").ptr, *test);
        } else {
            gretl_errmsg_sprintf(_("Extraneous character (0x%x) in data").ptr, *test);
        }
        err = 1;
    } else if (errno == ERANGE || val < int.min || val > int.max) {
        gretl_errmsg_sprintf(_("'%s' -- number out of range!").ptr, numstr);
        err = 1;
    }

    return err;
}

private int transpose_varname_used(const(char) * vname,
                                   DatasetRaw * dinfo,
                                   int imax)
{
    int i;

    for (i=0; i<imax; i++) {
        if (!strcmp(vname, dinfo.varname[i])) {
            return 1;
        }
    }

    return 0;
}

/**
 * transpose_data:
 * @dset: pointer to dataset information struct.
 *
 * Attempts to transpose the current dataset, so that each
 * variable becomes interpreted as an observation and each
 * observation as a variable.
 *
 * Returns: 0 on success, non-zero error code on error.
 */

int transpose_data(DatasetRaw * dset)
{
    DatasetRaw * tset;
    int k = dset.n + 1;
    int T = dset.v - 1;
    int i, t;

    tset = create_new_dataset(k, T, 0);
    if (tset == null) {
        return GretlError.alloc;
    }

    for (i=1; i<dset.v; i++) {
        for (t=0; t<dset.n; t++) {
            tset.Z[t+1][i-1] = dset.Z[i][t];
        }
    }

    for (t=0; t<dset.n; t++) {
        int kt = t + 1; /* C: inner "k", shadows the outer one */
        char * targ = tset.varname[kt];

        if (dset.S != null && dset.S[t][0] != '\0') {
            int err;

            *targ = '\0';
            strncat(targ, dset.S[t], VNAMELEN - 1);
            gretl_charsub(targ, ' ', '_');
            err = check_varname(targ);
            if (err) {
                sprintf(targ, "v%d", kt);
                gretl_error_clear();
            } else if (transpose_varname_used(targ, tset, kt)) {
                sprintf(targ, "v%d", kt);
            }
        } else {
            sprintf(targ, "v%d", kt);
        }
    }

    free_Z(dset);
    dset.Z = tset.Z;

    clear_datainfo(dset, DataClearCode.full);

    dset.v = k;
    dset.n = T;
    dset.t1 = 0;
    dset.t2 = dset.n - 1;

    dset.varname = tset.varname;
    dset.varinfo = tset.varinfo;

    dataset_obs_info_default(dset);

    free(tset);

    return 0;
}

void dataset_set_regular_markers(DatasetRaw * dset)
{
    dset.markers = cast(char) DatasetMarkerType.regular_markers;
}

private struct FiletypeInfo {
    GretlFileType type;
    const(char) * src;
}

/**
 * dataset_add_import_info:
 * @dset: pointer to dataset information struct.
 * @fname: the name of a file from which data have been imported.
 * @type: code representing the type of the file identified by
 * @fname.
 *
 * On successful import of data from some "foreign" format,
 * add a note to the "descrip" member of the new dataset
 * saying where it came from and when.
 */

void dataset_add_import_info(DatasetRaw * dset, const(char) * fname,
                             GretlFileType type)
{
    static immutable FiletypeInfo[10] ftypes = [
        FiletypeInfo(GretlFileType.csv,      "CSV"),
        FiletypeInfo(GretlFileType.gnumeric, "Gnumeric"),
        FiletypeInfo(GretlFileType.xls,      "Excel"),
        FiletypeInfo(GretlFileType.xlsx,     "Excel"),
        FiletypeInfo(GretlFileType.ods,      "Open Document"),
        FiletypeInfo(GretlFileType.wf1,      "Eviews"),
        FiletypeInfo(GretlFileType.dta,      "Stata"),
        FiletypeInfo(GretlFileType.sav,      "SPSS"),
        FiletypeInfo(GretlFileType.sas,      "SAS"),
        FiletypeInfo(GretlFileType.jmulti,   "JMulTi")
    ];
    int i, nt = cast(int) ftypes.length;
    const(char) * src = null;
    gchar * note = null;
    char[128] tstr = 0;

    for (i=0; i<nt; i++) {
        if (type == ftypes[i].type) {
            src = ftypes[i].src;
            break;
        }
    }

    if (src == null) {
        return;
    }

    print_time(tstr.ptr);

    if (g_utf8_validate(fname, -1, null)) {
        const(char) * p = strrslash(cast(char *) fname);

        if (p != null) {
            fname = p + 1;
        }
        note = g_strdup_printf(_("Data imported from %s file '%s', %s\n").ptr,
                               src, fname, tstr.ptr);
    } else {
        note = g_strdup_printf(_("Data imported from %s, %s\n").ptr,
                               src, tstr.ptr);
    }

    if (note != null) {
        if (dset.descrip == null) {
            dset.descrip = gretl_strdup(note);
        } else {
            int dlen = cast(int) strlen(dset.descrip);
            int nlen = cast(int) strlen(note);
            char * tmp = cast(char *) realloc(dset.descrip, dlen + nlen + 5);

            if (tmp != null) {
                dset.descrip = tmp;
                strcat(dset.descrip, "\n\n");
                strncat(dset.descrip, note, nlen);
            }
        }
        g_free(note);
    }
}

private int is_weekend(int t, int pd, int sat0, int sun0)
{
    int sat = 0, sun = 0;

    /* this is intended to identify weekend days for
       both 7- and 6-day data */

    if (sat0 >= 0) {
        sat = (t - sat0) % pd == 0;
    }

    if (!sat && sun0 >= 0) {
        sun = (t - sun0) % pd == 0;
    }

    return sat || sun;
}

/* Scan imported daily data for missing values, so as to
   be able to offer the user some options.

   Return values:

   0 : no missing values
   1 : all weekend data are missing (or possibly just all
       Saturdays, or just all Sundays), but no weekday
       data missing
   2 : as with 1, but also some weekdays missing
   3 : scattering of weekend and/or weekday data missing
*/

int analyse_daily_import(const(DatasetRaw) * dset, Prn * prn)
{
    int all_weekends_blank = 0;
    int blank_weekends = 0;
    int blank_weekdays = 0;
    int n_weekdays = 0;
    int n_weekend_days = 0;
    int sat0 = -1, sun0 = -1;
    int i, t, pd = dset.pd;
    int all_missing, weekend;
    int ret = 0;

    if (pd > 5) {
        char[OBSLEN] datestr = 0;
        int wkday;

        /* start by finding first Sat and/or Sun */
        for (t=0; t<dset.n; t++) {
            ntolabel(datestr.ptr, t, dset);
            wkday = weekday_from_date(datestr.ptr);
            if (wkday == 6 && sat0 < 0) {
                sat0 = t;
            } else if (wkday == 7 && sun0 < 0) {
                sun0 = t;
            }
            if (sat0 >= 0 && sun0 >= 0) {
                break;
            } else if (dset.pd == 6 && (sat0 >= 0 || sun0 >= 0)) {
                break;
            }
        }
        all_weekends_blank = 1; /* may be revised below */
    } else {
        /* there are no weekend days in 5-day data */
        weekend = 0;
    }

    for (t=0; t<dset.n; t++) {
        if (pd > 5) {
            weekend = is_weekend(t, pd, sat0, sun0);
        }
        all_missing = 1;
        for (i=1; i<dset.v; i++) {
            if (!na(dset.Z[i][t])) {
                all_missing = 0;
                break;
            }
        }
        if (weekend) {
            n_weekend_days++;
            if (!all_missing) {
                /* not all weekend data are missing */
                all_weekends_blank = 0;
            }
        } else {
            n_weekdays++;
        }
        if (all_missing) {
            if (weekend) {
                blank_weekends++;
            } else {
                blank_weekdays++;
            }
        }
    }

    if (all_weekends_blank) {
        double misspc = 100.0 * blank_weekdays / cast(double) n_weekdays;

        if (pd == 7) {
            pputs(prn, _("This dataset is on 7-day calendar, but weekends are blank.").ptr);
        } else {
            pprintf(prn, _("This dataset is on 6-day calendar, but %s are blank.").ptr,
                    sat0 >= 0 ? _("Sundays").ptr : _("Saturdays").ptr);
        }
        ret = 1;
        if (misspc > 0.0) {
            pputc(prn, '\n');
            pprintf(prn, "%s ", _("In addition,").ptr);
            if (misspc >= 0.01) {
                pprintf(prn, _("%.2f percent of weekday observations are missing.").ptr,
                        misspc);
            } else {
                pprintf(prn, _("%g percent of weekday observations are missing.").ptr,
                        misspc);
            }
            if (misspc < 10.0) {
                ret = 2;
            }
        }
    } else if (blank_weekdays || blank_weekends) {
        int ndays = n_weekdays + n_weekend_days;
        int nmiss = blank_weekdays + blank_weekends;
        double misspc = 100.0 * nmiss / cast(double) ndays;

        if (misspc >= 0.01) {
            pprintf(prn, _("%.2f percent of daily observations are missing.").ptr,
                    misspc);
        } else {
            pprintf(prn, _("%g percent of daily observations are missing.").ptr,
                    misspc);
        }
        if (misspc < 10) {
            ret = 3;
        }
    }

    return ret;
}
