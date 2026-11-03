/* Port of csvdata.c (csvdata.h is in csvdatah.d) */

module gretl.csvdata;


import core.stdc.config : c_long;
import core.stdc.ctype : isalnum, isalpha, isdigit, isprint, ispunct, isspace, isupper, tolower, toupper;
import core.stdc.errno : ERANGE, errno;
import core.stdc.stdio : EOF, FILE, fclose, fgetc, fgets, fprintf, fputs, fseek, ftell, getc, printf,
    rewind, sprintf, snprintf, sscanf, stderr, ungetc, SEEK_SET;
import core.stdc.stdlib : atof, atoi, calloc, free, malloc, realloc, strtod, strtol;
import core.stdc.string : memcpy, memmove, memset, strcat, strchr, strcmp, strcpy, strlen, strncat, strncmp,
    strpbrk, strrchr, strspn, strstr;

import core.stdc.math : ceil, fabs;
import core.stdc.stdio : perror;
import core.stdc.stdlib : abs;
import core.stdc.time : tm;
version (Posix) import core.sys.posix.time : strptime; else import gretl.winport : strptime;
import etc.c.zlib : gzFile, gzclose, gzgetc, gzgets, gzread, gzrewind, gzungetc;
import gretl.csvdatah : NON_NUMERIC;
import gretl.dataio : check_varname, dataset_add_import_info, get_merge_opts, merge_or_replace_data;
import gretl.dataseth : SeriesTable;
import gretl.gerror : GError, g_error_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.matrix_extrah : MMissingCode;
import gretl.options : get_optval_string, incompatible_options;
import gretl.libset : get_csv_na_read_string, get_data_export_decpoint, get_data_export_delimiter, gretl_messages_on;
import gretl.strutilsh : CTRLZ;
import gretl.gstrfuncs : g_convert, g_locale_to_utf8, g_utf8_get_char, g_utf8_next_char, g_utf8_strlen;
import gretl.glib : gchar, gsize;
import gretl.gretl_prn : gretl_print_destroy, gretl_print_reset_buffer;
import gretl.dataio : ntolabel;
import gretl.dataseth;
import gretl.dataset;
import gretl.dataioh : DataClearCode, GretlFileType;
import gretl.gmem : g_free;
import gretl.gretl_errors : errmsg, gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear, gretl_warnmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_append_term, gretl_list_from_string, gretl_list_new, in_gretl_list;
import gretl.gretl_matrixh;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gstrfuncs : g_strdup_printf, g_utf8_validate;
import gretl.join_privh : JoinSpec;
import gretl.libgretlh;
import gretl.libset : libset_get_bool, libset_get_int;
import gretl.libseth : SetKey;
import gretl.missingh : NADBL, na;
import gretl.strutils;
import gretl.gretl_paths : gretl_dotdir, gretl_gzopen, gretl_maybe_prepend_dir, gretl_remove;
import gretl.gretl_string_table : gretl_string_table_destroy, gretl_string_table_finalize, gretl_string_table_index, gretl_string_table_new, gretl_string_table_save, gretl_string_table_validate;
import gretl.calendar : calendar_obs_number, get_epoch_day, guess_daily_pd, weekday_from_date;
import gretl.gretl_intl : get_local_decpoint, gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale, gretl_recode_file, iso_to_ascii;
import gretl.genmain : gretl_reserved_word;
import gretl.matrix_extra : gretl_matrix_data_subset;
import gretl.usermat : get_matrix_by_name;

private enum CSVSTRLEN = 128;

private enum CsvFlags {
    havedata = 1 << 0,
    gotdelim = 1 << 1,
    gottab   = 1 << 2,
    gotsemi  = 1 << 3,
    gotpipe  = 1 << 4,
    blank1   = 1 << 5,
    obs1     = 1 << 6,
    trail    = 1 << 7,
    autoname = 1 << 8,
    reversed = 1 << 9,
    dotsub   = 1 << 10,
    allcols  = 1 << 11,
    bom      = 1 << 12,
    verbose  = 1 << 13,
    thousep  = 1 << 14,
    noheader = 1 << 15,
    quotes   = 1 << 16,
    as_mat   = 1 << 17
}

struct CsvProbe {
    DatasetRaw * dset; /* more info might be wanted */
}

struct CsvData {
    CsvFlags flags;
    char delim;
    char decpoint;
    byte thousep;       /* C char is signed: -1 is used as a flag */
    char qchar;
    int markerpd;
    int maxlinelen;
    int real_n;
    char * line;
    DatasetRaw * dset;
    int ncols, nrows;
    c_long datapos;
    char[CSVSTRLEN] str;
    char[8] skipstr;
    int * codelist;
    char * descrip;
    const(char) * user_na;
    GretlStringTable * st;
    int * cols_list;
    int * width_list;
    const(GretlMatrixRaw) * rowmask;
    int masklen;
    JoinSpec * jspec;   /* info used for "join" command */
    CsvProbe * probe;   /* used in connection with "join" */
}

/* the C accessor macros, as private functions */

private int csv_has_trailing_comma(const(CsvData) * c) { return (c.flags & CsvFlags.trail) != 0; }
private int csv_has_obs_column(const(CsvData) * c)     { return (c.flags & CsvFlags.obs1) != 0; }
private int csv_has_blank_column(const(CsvData) * c)   { return (c.flags & CsvFlags.blank1) != 0; }
private int csv_got_tab(const(CsvData) * c)            { return (c.flags & CsvFlags.gottab) != 0; }
private int csv_got_semi(const(CsvData) * c)           { return (c.flags & CsvFlags.gotsemi) != 0; }
private int csv_got_pipe(const(CsvData) * c)           { return (c.flags & CsvFlags.gotpipe) != 0; }
private int csv_got_delim(const(CsvData) * c)          { return (c.flags & CsvFlags.gotdelim) != 0; }
private int csv_autoname(const(CsvData) * c)           { return (c.flags & CsvFlags.autoname) != 0; }
private int csv_skip_col_1(const(CsvData) * c)         { return (c.flags & (CsvFlags.obs1 | CsvFlags.blank1)) != 0; }
private int csv_have_data(const(CsvData) * c)          { return (c.flags & CsvFlags.havedata) != 0; }
private int csv_data_reversed(const(CsvData) * c)      { return (c.flags & CsvFlags.reversed) != 0; }
private int csv_do_dotsub(const(CsvData) * c)          { return (c.flags & CsvFlags.dotsub) != 0; }
private int csv_all_cols(const(CsvData) * c)           { return (c.flags & CsvFlags.allcols) != 0; }
private int csv_has_bom(const(CsvData) * c)            { return (c.flags & CsvFlags.bom) != 0; }
private int csv_is_verbose(const(CsvData) * c)         { return (c.flags & CsvFlags.verbose) != 0; }
private int csv_scrub_thousep(const(CsvData) * c)      { return (c.flags & CsvFlags.thousep) != 0; }
private int csv_no_header(const(CsvData) * c)          { return (c.flags & CsvFlags.noheader) != 0; }
private int csv_keep_quotes(const(CsvData) * c)        { return (c.flags & CsvFlags.quotes) != 0; }
private int csv_as_matrix(const(CsvData) * c)          { return (c.flags & CsvFlags.as_mat) != 0; }

private void csv_set_trailing_comma(CsvData * c)   { c.flags |= CsvFlags.trail; }
private void csv_unset_trailing_comma(CsvData * c) { c.flags &= ~CsvFlags.trail; }
private void csv_set_obs_column(CsvData * c)       { c.flags |= CsvFlags.obs1; }
private void csv_set_blank_column(CsvData * c)     { c.flags |= CsvFlags.blank1; }
private void csv_set_got_tab(CsvData * c)          { c.flags |= CsvFlags.gottab; }
private void csv_set_got_semi(CsvData * c)         { c.flags |= CsvFlags.gotsemi; }
private void csv_set_got_pipe(CsvData * c)         { c.flags |= CsvFlags.gotpipe; }
private void csv_set_got_delim(CsvData * c)        { c.flags |= CsvFlags.gotdelim; }
private void csv_set_autoname(CsvData * c)         { c.flags |= CsvFlags.autoname; }
private void csv_set_data_reversed(CsvData * c)    { c.flags |= CsvFlags.reversed; }
private void csv_set_dotsub(CsvData * c)           { c.flags |= CsvFlags.dotsub; }
private void csv_set_all_cols(CsvData * c)         { c.flags |= CsvFlags.allcols; }
private void csv_set_has_bom(CsvData * c)          { c.flags |= CsvFlags.bom; }
private void csv_set_verbose(CsvData * c)          { c.flags |= CsvFlags.verbose; }
private void csv_set_scrub_thousep(CsvData * c)    { c.flags |= CsvFlags.thousep; }
private void csv_set_no_header(CsvData * c)        { c.flags |= CsvFlags.noheader; }
private void csv_unset_keep_quotes(CsvData * c)    { c.flags &= ~CsvFlags.quotes; }
private void csv_set_as_matrix(CsvData * c)        { c.flags |= CsvFlags.as_mat; }

private int csv_skip_bad(const(CsvData) * c)        { return c.skipstr[0] != '\0'; }
private int csv_has_non_numeric(const(CsvData) * c) { return c.st != null; }

private int fixed_format(const(CsvData) * c) { return c.cols_list != null && c.width_list != null; }
private int cols_subset(const(CsvData) * c)  { return c.cols_list != null && c.width_list == null; }
private int rows_subset(const(CsvData) * c)  { return c.rowmask != null; }

private int joining(const(CsvData) * c) { return c.jspec != null; }
private int probing(const(CsvData) * c) { return c.probe != null; }

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    c_long gzseek(gzFile file, c_long offset, int whence);   /* z_off_t is a long here */
    c_long gztell(gzFile file);
    int timecol_get_format(const(DatasetRaw) * dset, int v, char ** pfmt, int * q);
    int try_http(const(char) * s, char * fname, int * err);
}

/* file-scope global */
private __gshared char[8] import_na = 0;

/* for use in gretl_join.c */

DatasetRaw * csvdata_get_dataset(CsvData * c)
{
    return c.dset;
}

/* shared with gretl_join.c */

void csvdata_free(CsvData * c)
{
    if (c == null) {
        return;
    }

    if (c.descrip != null) {
        free(c.descrip);
    }

    if (c.st != null) {
        gretl_string_table_destroy(c.st);
    }

    if (c.codelist != null) {
        free(c.codelist);
    }

    if (c.line != null) {
        free(c.line);
    }

    if (c.cols_list != null) {
        free(c.cols_list);
        free(c.width_list);
    }

    destroy_dataset(c.dset);

    free(c);
}

private CsvData * csvdata_new(DatasetRaw * dset)
{
    CsvData * c = cast(CsvData *) malloc(CsvData.sizeof);

    if (c == null) {
        return null;
    }

    c.flags = CsvFlags.quotes;
    c.delim = '\t';
    c.thousep = 0;
    c.qchar = 0;
    c.markerpd = -1;
    c.maxlinelen = 0;
    c.real_n = 0;
    c.line = null;
    c.dset = null;
    c.ncols = 0;
    c.nrows = 0;
    c.datapos = 0;
    c.str[0] = '\0';
    c.skipstr[0] = '\0';
    c.codelist = null;
    c.descrip = null;
    c.user_na = null;
    c.st = null;
    c.cols_list = null;
    c.width_list = null;
    c.rowmask = null;
    c.masklen = 0;

    if (strcmp(import_na.ptr, "default")) {
        c.user_na = import_na.ptr;
    }

    c.jspec = null;
    c.probe = null;

    c.dset = datainfo_new();

    if (c.dset == null) {
        free(c);
        c = null;
    } else {
        c.delim = get_data_export_delimiter();
        c.decpoint = get_data_export_decpoint();
        if (dset != null && dset.Z != null) {
            c.flags |= CsvFlags.havedata;
        }
    }

    return c;
}

private int * cols_list_from_matrix(const(char) * s, int * err)
{
    GretlMatrixRaw * m = get_matrix_by_name(s);
    int i, n = gretl_vector_get_length(m);
    int * list = null;

    if (n == 0) {
        *err = GretlError.data;
    } else {
        list = gretl_list_new(n);
        if (list == null) {
            *err = GretlError.alloc;
        } else {
            for (i=0; i<n; i++) {
                list[i+1] = cast(int) gretl_vector_get(m, i);
            }
        }
    }

    return list;
}

/* The interpretation of the "cols" specification depends on
   @opt: if this includes OPT_L then it should provide a 1-based
   list of columns to be read; but if @opt includes OPT_F it
   should provide a fixed-format spec, consisting of pairs
   (start column, width).
*/

private int csvdata_add_cols_list(CsvData * c, const(char) * s,
                                  Gretlopt opt)
{
    int delimited = (opt & Gretlopt.OPT_L) != 0;
    int * list = null;
    int * clist = null;
    int * wlist = null;
    int i, n, m = 0;
    int err = 0;

    if (get_matrix_by_name(s)) {
        list = cols_list_from_matrix(s, &err);
    } else {
        list = gretl_list_from_string(s, &err);
    }

    if (!err) {
        n = list[0];
        if (n == 0) {
            err = GretlError.data;
        } else if (delimited) {
            m = n;
            clist = list;
        } else {
            /* fixed format: we need two lists */
            if (n % 2 != 0) {
                err = GretlError.data;
            } else {
                m = n / 2;
                clist = gretl_list_new(m);
                wlist = gretl_list_new(m);
                if (clist == null || wlist == null) {
                    err = GretlError.alloc;
                } else {
                    int j = 1;

                    for (i=1; i<=n; i+=2, j++) {
                        clist[j] = list[i];
                        wlist[j] = list[i+1];
                    }
                }
            }
        }
    }

    /* clist = column (start) list: must be a set of increasing
       positive integers; and wlist = respective column widths,
       must all be positive, if present
    */

    for (i=1; i<=m && !err; i++) {
        if (clist[i] <= 0 || (i > 1 && clist[i] <= clist[i-1])) {
            err = GretlError.data;
        } else if (wlist != null && wlist[i] <= 0) {
            err = GretlError.data;
        } else if (wlist != null && wlist[i] >= CSVSTRLEN) {
            fprintf(stderr, "Warning: field %d too wide (%d), truncating\n",
                    i, wlist[i]);
            wlist[i] = CSVSTRLEN - 1;
        }
    }

    if (list != clist) {
        free(list);
    }

    if (!err) {
        c.cols_list = clist;
        c.width_list = wlist;
    } else {
        free(clist);
        free(wlist);
        if (err == GretlError.data) {
            gretl_errmsg_set(_("Invalid column specification").ptr);
        }
    }

    return err;
}

private int csvdata_add_row_mask(CsvData * c, const(char) * s)
{
    int err = 0;

    c.rowmask = get_matrix_by_name(s);
    if (c.rowmask == null) {
        gretl_errmsg_sprintf(_("'%s': no such matrix").ptr, s);
        err = GretlError.data;
    } else {
        c.masklen = gretl_vector_get_length(c.rowmask);
        if (c.masklen == 0) {
            err = GretlError.nonconf;
        }
    }

    return err;
}

private int n_from_row_mask(CsvData * c)
{
    int i, n = 0;

    for (i=0; i<c.masklen && i<=c.nrows; i++) {
        if (gretl_vector_get(c.rowmask, i) != 0) {
            n++;
        }
    }

    return n;
}

private int add_obs_marker(DatasetRaw * dset, int n)
{
    char ** S = cast(char **) realloc(dset.S, n * (char *).sizeof);
    int err = 0;

    if (S == null) {
        err = GretlError.alloc;
    } else {
        dset.S = S;
        dset.S[n-1] = cast(char *) malloc(OBSLEN);
        if (dset.S[n-1] == null) {
            err = GretlError.alloc;
        } else {
            strcpy(dset.S[n-1], "NA");
        }
    }

    return err;
}

private int add_single_obs(DatasetRaw * dset)
{
    double * x;
    int i, err = 0;

    for (i=0; i<dset.v && !err; i++) {
        x = cast(double *) realloc(dset.Z[i], (dset.n + 1) * double.sizeof);
        if (x != null) {
            dset.Z[i] = x;
        } else {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        dset.n += 1;
        dset.Z[0][dset.n - 1] = 1.0;
        for (i=1; i<dset.v; i++) {
            dset.Z[i][dset.n - 1] = NADBL;
        }
        if (dset.S != null) {
            err = add_obs_marker(dset, dset.n);
        }
    }

    return err;
}

private int pad_weekly_data(DatasetRaw * dset, int add)
{
    int oldn = dset.n;
    int ttarg, offset = 0, skip = 0;
    int i, s, t, tc, err;

    err = dataset_add_observations(dset, add, Gretlopt.OPT_A);

    if (!err) {
        for (t=0; t<oldn; t++) {
            tc = calendar_obs_number(dset.S[t], dset, 0) - offset;
            if (tc != t) {
                skip = tc - t;
                fprintf(stderr, "Gap of size %d at original t = %d\n", skip, t);
                offset += skip;
                ttarg = oldn - 1 + offset;
                for (s=0; s<oldn-t+skip; s++) {
                    for (i=1; i<dset.v; i++) {
                        if (s < oldn - t) {
                            if (s == 0 || s == oldn-t-1) {
                                fprintf(stderr, "shifting obs %d to obs %d\n",
                                        ttarg-skip, ttarg);
                            }
                            dset.Z[i][ttarg] = dset.Z[i][ttarg - skip];
                        } else {
                            fprintf(stderr, "inserting NA at obs %d\n", ttarg);
                            dset.Z[i][ttarg] = NADBL;
                        }
                    }
                    ttarg--;
                }
            }
        }
    }

    return err;
}

/* FIXME the following needs to be made more flexible? */

private int csv_weekly_data(DatasetRaw * dset)
{
    char * lbl2 = dset.S[dset.n - 1];
    int ret = 1;
    int misscount = 0;
    int t, tc;

    for (t=0; t<dset.n; t++) {
        tc = calendar_obs_number(dset.S[t], dset, 0) - misscount;
        if (tc != t) {
            misscount += tc - t;
        }
    }

    if (misscount > 0) {
        double missfrac = cast(double) misscount / dset.n;

        fprintf(stderr, "nobs = %d, misscount = %d (%.2f%%)\n",
                dset.n, misscount, 100.0 * missfrac);
        if (missfrac > 0.05) {
            ret = 0;
        } else {
            int Tc = calendar_obs_number(lbl2, dset, 0) + 1;
            int altmiss = Tc - dset.n;

            fprintf(stderr, "check: Tc = %d, missing = %d\n", Tc, altmiss);
            if (altmiss != misscount) {
                ret = 0;
            } else if (dset.Z != null) {
                int err;

                fprintf(stderr, "OK, consistent\n");
                err = pad_weekly_data(dset, misscount);
                if (err) ret = 0;
            }
        }
    }

    return ret;
}

/* C: #define DAY_DEBUG 1, so the diagnostics below are always on */

private int check_daily_dates(DatasetRaw * dset, int * pd,
                              int * reversed, Prn * prn)
{
    int T = dset.n;
    char * lbl1 = dset.S[0];
    char * lbl2 = dset.S[T - 1];
    int fulln = 0, n, t, nbak;
    int alt_pd = 0;
    int oldpd = dset.pd;
    double oldsd0 = dset.sd0;
    uint ed1, ed2;
    int nmiss = 0, err = 0;

    *pd = 0;

    ed1 = get_epoch_day(lbl1);
    ed2 = get_epoch_day(lbl2);
    if (ed1 <= 0 || ed2 <= 0) {
        err = 1;
    }

    fprintf(stderr, "check_daily_dates: '%s' -> %d, '%s' -> %d\n",
            lbl1, cast(int) ed1, lbl2, cast(int) ed2);

    dset.pd = guess_daily_pd(dset);
    dset.structure = TsCodes.time_series;

    fprintf(stderr, "guessed at daily pd = %d\n", dset.pd);

    if (!err) {
        if (ed2 < ed1) {
            fprintf(stderr, "check_daily_dates: data are reversed?\n");
            dset.sd0 = ed2;
            *reversed = 1;
        } else {
            dset.sd0 = ed1;
        }
    }

 recompute:

    alt_pd = 0;
    nbak = 0;

    if (!err) {
        uint n1 = (*reversed)? ed2 : ed1;
        uint n2 = (*reversed)? ed1 : ed2;

        fulln = n2 - n1 + 1;

        if (T > fulln) {
            err = 1;
        } else {
            nmiss = fulln - T;
            pprintf(prn, _("Observations: %d; days in sample: %d\n").ptr,
                    T, fulln);
            if (nmiss > 300 * T) {
                pprintf(prn, _("Probably annual data\n").ptr);
                *pd = 1;
            } else if (nmiss > 50 * T) {
                pprintf(prn, _("Probably quarterly data\n").ptr);
                *pd = 4;
            } else if (nmiss > 20 * T) {
                pprintf(prn, _("Probably monthly data\n").ptr);
                *pd = 12;
            } else if (nmiss > 3 * T) {
                pprintf(prn, _("Probably weekly data\n").ptr);
                *pd = dset.pd = 52;
            } else {
                pprintf(prn, _("Missing daily rows: %d\n").ptr, nmiss);
            }
        }
    }

    nbak = 0;

    for (t=0; t<dset.n && !err; t++) {
        int wd, s = (*reversed)? (dset.n - 1 - t) : t;

        wd = weekday_from_date(dset.S[s]);

        if (dset.pd == 5 && (wd == 6 || wd == 7)) {
            /* Got Sat or Sun, can't be 5-day daily? */
            alt_pd = wd;
            pprintf(prn, _("Found a Saturday (%s): re-trying with pd = %d\n").ptr,
                    dset.S[s], alt_pd);
            break;
        } else if (dset.pd == 6 && wd == 7) {
            /* Got Sun, can't be 6-day daily? */
            alt_pd = wd;
            pprintf(prn, _("Found a Sunday (%s): re-trying with pd = %d\n").ptr,
                    dset.S[s], alt_pd);
            break;
        }

        n = calendar_obs_number(dset.S[s], dset, 0);
        if (n < t) {
            pprintf(prn, (_("Daily dates error at t = %d:\n")
                    ~ _("  calendar_obs_number() for '%s' = %d but t = %d\n")).ptr,
                    t, dset.S[s], n, t);
            err = 1;
        } else if (n > fulln - 1) {
            pprintf(prn, _("Error: date '%s' out of bounds\n").ptr, dset.S[s]);
            err = 1;
        } else if (nbak > 0 && n == nbak) {
            pprintf(prn, _("Error: date '%s' is repeated\n").ptr, dset.S[s]);
            err = 1;
        }
        nbak = n;
    }

    if (alt_pd > 0) {
        dset.pd = alt_pd;
        goto recompute;
    }

    if (err) {
        dset.pd = oldpd;
        dset.sd0 = oldsd0;
        dset.structure = TsCodes.cross_section;
    } else {
        strcpy(dset.stobs.ptr, (*reversed)? lbl2 : lbl1);
        strcpy(dset.endobs.ptr, (*reversed)? lbl1 : lbl2);
        dset.t2 = dset.n - 1;
        if (nmiss > 0 && *pd == 0) {
            dset.markers = cast(char) DatasetMarkerType.daily_date_strings;
        }
    }

    fprintf(stderr, "check_daily_dates: daily pd = %d, reversed = %d, err = %d\n",
            dset.pd, *reversed, err);

    return (err)? -1 : dset.pd;
}

/* convert from daily date label to a lower frequency --
   annual, monthly or quarterly -- if @pd indicates this
   is required
*/

private void convert_daily_label(char * targ, const(char) * src,
                                 int pd)
{
    int y, m, d;

    sscanf(src, YMD_READ_FMT.ptr, &y, &m, &d);

    if (pd == 1) {
        sprintf(targ, "%d", y);
    } else if (pd == 12) {
        sprintf(targ, "%d:%02d", y, m);
    } else if (pd == 4) {
        sprintf(targ, "%d:%d", y, m / 3 + (m % 3 != 0));
    }
}

/* There's a special case (ugh!) where observation strings are
   given as in monthly data, but the frequency is in fact
   quarterly, as in:

   1947.06
   1947.09
   1947.12
   1948.03

   we'll make a brave attempt to handle this.
*/

private bool fakequarter(int m)
{
    return m==3 || m==6 || m==9 || m==12;
}

private int consistent_qm_labels(DatasetRaw * dset, int reversed,
                                 int convert_pd, char * skipstr,
                                 int * ppd, const(char) * fmt,
                                 int * extra_zero, Prn * prn)
{
    char[16] bad = 0;
    char[8] skip = 0;
    char[OBSLEN] label = 0;
    int Ey; /* expected year */
    int Ep; /* expected sub-period */
    int t, s, yr, per;
    int pmin = 1;
    int pd, pd0;
    int ret = 1;

    pd = pd0 = *ppd;

 restart:

    s = reversed ? (dset.n - 1) : 0;

    if (convert_pd) {
        convert_daily_label(label.ptr, dset.S[s], pd);
    } else {
        strcpy(label.ptr, dset.S[s]);
    }

    if (sscanf(label.ptr, fmt, &yr, &per) != 2) {
        return 0;
    }

    for (t=1; t<dset.n; t++) {
        s = (reversed)? (dset.n - 1 - t) : t;
        Ey = (per == pd)? yr + 1 : yr;
        Ep = (per == pd)? pmin : per + pmin;

        if (convert_pd) {
            convert_daily_label(label.ptr, dset.S[s], pd);
        } else {
            strcpy(label.ptr, dset.S[s]);
        }

        if (sscanf(label.ptr, fmt, &yr, &per) != 2) {
            ret = 0;
        } else if (Ep == 1 && pd == pd0 && per == pd + 1
                   && skipstr != null) {
            skip[0] = bad[0] = '\0';
            strncat(skip.ptr, label.ptr + 4, 7);
            strncat(bad.ptr, label.ptr, OBSLEN-1);
            pd = pd0 + 1;
            goto restart;
        } else if (per == Ep + 2 && pmin == 1 && fakequarter(per)) {
            bad[0] = '\0';
            strncat(bad.ptr, label.ptr, OBSLEN-1);
            pmin = 3;
            goto restart;
        } else if (pd == 12 && Ep == 5 && per == 1 && yr == Ey + 1) {
            /* apparently monthly but really quarterly? */
            pprintf(prn, _("   \"%s\": quarterly date with spurious zero?\n").ptr, label.ptr);
            *extra_zero = 1;
            *ppd = pd0 = pd = 4;
            goto restart;
        } else if (yr != Ey || per != Ep) {
            ret = 0;
        }

        if (!ret) {
            pprintf(prn, "   %s: %s\n", label.ptr, _("not a consistent date").ptr);
            break;
        }
    }

    if (ret) {
        if (pmin == 3) {
            pprintf(prn, "   \"%s\": %s\n", bad.ptr, _("quarterly data pretending to be monthly?").ptr);
            *ppd = 4;
        } else if (pd == pd0 + 1) {
            pprintf(prn, "   \"%s\": %s\n", bad.ptr, _("BLS-type nonsense? Trying again").ptr);
            strcpy(skipstr, skip.ptr);
        }
    }

    return ret;
}

private int consistent_year_labels(const(DatasetRaw) * dset,
                                   int reversed,
                                   int convert_pd)
{
    char[OBSLEN] label = 0;
    int s, t, yr, yprev;
    int ret = 1;

    s = (reversed)? (dset.n - 1) : 0;
    yprev = atoi(dset.S[s]);

    for (t=1; t<dset.n; t++) {
        s = reversed ? (dset.n - 1 - t) : t;
        if (convert_pd) {
            convert_daily_label(label.ptr, dset.S[s], 1);
            yr = atoi(label.ptr);
        } else {
            yr = atoi(dset.S[s]);
        }
        if (yr != yprev + 1) {
            ret = 0;
            break;
        }
        yprev = yr;
    }

    return ret;
}

/* check for all 1s in first column of dates: this may
   indicate start-of-period dates, day first */

private int all_day_ones(DatasetRaw * dset)
{
    int t;

    for (t=1; t<dset.n; t++) {
        if (atoi(dset.S[t]) != 1) {
            return 0;
        } else if (t > 31) {
            /* "1" can't mean January */
            return 1;
        }
    }

    return 0;
}

private enum DateOrder {
    yyyymmdd = 1,
    mmddyyyy,
    ddmmyyyy
}

private int get_date_order(int f0, int fn, DatasetRaw * dset)
{
    if (f0 > 31 || fn > 31) {
        /* first field must be year */
        return DateOrder.yyyymmdd;
    } else if (f0 > 12 || fn > 12) {
        /* first field must be day */
        return DateOrder.ddmmyyyy;
    } else if (f0 == 1 && fn == 1 && all_day_ones(dset)) {
        /* start-of-period dates, day first? */
        return DateOrder.ddmmyyyy;
    } else {
        /* could be wrong here */
        return DateOrder.mmddyyyy;
    }
}

private void retransform_daily_dates(DatasetRaw * dset)
{
    int t, y, m, d;

    /* we apparently guessed wrongly at MMDDYYYY, so
       put the dates back as they were for another try,
       at DDMMYYYY.
    */

    for (t=0; t<dset.n; t++) {
        sscanf(dset.S[t], YMD_READ_FMT.ptr, &y, &d, &m);
        sprintf(dset.S[t], YMD_WRITE_FMT.ptr, d, m, y);
    }
}

private int transform_daily_dates(DatasetRaw * dset, int dorder,
                                  char sep)
{
    char * label;
    char[16] fmt = 0;
    int t, yr, mon, day;
    int n, err = 0;

    if (sep > 0) {
        sprintf(fmt.ptr, "%%d%c%%d%c%%d", sep, sep);
    } else {
        strcpy(fmt.ptr, "%4d%2d%2d");
    }

    for (t=0; t<dset.n && !err; t++) {
        label = dset.S[t];
        if (dorder == DateOrder.yyyymmdd) {
            n = sscanf(label, fmt.ptr, &yr, &mon, &day);
        } else if (dorder == DateOrder.ddmmyyyy) {
            n = sscanf(label, fmt.ptr, &day, &mon, &yr);
        } else {
            n = sscanf(label, fmt.ptr, &mon, &day, &yr);
        }
        /* Guard against oversized field values: @fmt places no
           width limit on the sscanf conversions above (the
           separator-derived form has none, and even the fallback
           "%4d%2d%2d" only bounds the width of the *input* field,
           not the numeric value produced), so a crafted label can
           yield a mon/day/yr combination that would not fit back
           into the fixed OBSLEN-byte label buffer once reformatted
           by YMD_WRITE_Y2_FMT below. Reject any row whose parsed
           values fall outside sane calendar-date bounds before
           attempting to write them back.
        */
        if (n == 3 && mon >= 1 && mon <= 12 && day >= 1 && day <= 31 &&
            yr >= 0 && yr <= 9999) {
            snprintf(label, OBSLEN, YMD_WRITE_Y2_FMT.ptr, yr, mon, day);
        } else {
            err = 1;
        }
    }

    return err;
}

void reverse_data(DatasetRaw * dset, Prn * prn)
{
    char[OBSLEN] tmp = 0;
    double x;
    int T = dset.n / 2;
    int i, t, s;

    pprintf(prn, _("reversing the data!\n").ptr);

    for (t=0; t<T; t++) {
        s = dset.n - 1 - t;
        for (i=1; i<dset.v; i++) {
            x = dset.Z[i][t];
            dset.Z[i][t] = dset.Z[i][s];
            dset.Z[i][s] = x;
        }
        if (dset.S != null) {
            strcpy(tmp.ptr, dset.S[t]);
            strcpy(dset.S[t], dset.S[s]);
            strcpy(dset.S[s], tmp.ptr);
        }
    }
}

private int csv_daily_date_check(DatasetRaw * dset, int * reversed,
                                 char * skipstr, Prn * prn)
{
    int[3] d1 = 0;
    int[3] d2 = 0;
    char s1 = 0, s2 = 0;
    char * lbl1 = dset.S[0];
    char * lbl2 = dset.S[dset.n - 1];
    int dorder = 0;

    if ((sscanf(lbl1, "%d%c%d%c%d", &d1[0], &s1, &d1[1], &s2, &d1[2]) == 5 &&
         sscanf(lbl2, "%d%c%d%c%d", &d2[0], &s1, &d2[1], &s2, &d2[2]) == 5 &&
         s1 == s2 && ispunct(s1)) ||
        (sscanf(lbl1, "%4d%2d%2d", &d1[0], &d1[1], &d1[2]) == 3 &&
         sscanf(lbl2, "%4d%2d%2d", &d2[0], &d2[1], &d2[2]) == 3)) {
        int mon1, day1;
        int mon2, day2;
        int pd, ret = 0;

        dorder = get_date_order(d1[0], d2[0], dset);

    tryagain:

        if (dorder == DateOrder.yyyymmdd) {
            pputs(prn, _("Trying date order YYYYMMDD\n").ptr);
            mon1 = d1[1];
            day1 = d1[2];
            mon2 = d2[1];
            day2 = d2[2];
        } else if (dorder == DateOrder.ddmmyyyy) {
            pputs(prn, _("Trying date order DDMMYYYY\n").ptr);
            day1 = d1[0];
            mon1 = d1[1];
            day2 = d2[0];
            mon2 = d2[1];
        } else {
            pputs(prn, _("Trying date order MMDDYYYY\n").ptr);
            mon1 = d1[0];
            day1 = d1[1];
            mon2 = d2[0];
            day2 = d2[1];
        }

        if (mon1 > 0 && mon1 < 13 &&
            mon2 > 0 && mon2 < 13 &&
            day1 > 0 && day1 < 32 &&
            day2 > 0 && day2 < 32) {
            /* looks promising for calendar dates, but check
               further if we don't have the canonical order
               or separator
            */
            if (dorder != DateOrder.yyyymmdd || s1 != '-') {
                if (transform_daily_dates(dset, dorder, s1)) {
                    return -1;
                }
                s1 = '-';
            }
            pprintf(prn, _("Could be %s - %s\n").ptr, lbl1, lbl2);
            ret = check_daily_dates(dset, &pd, reversed, prn);
            if (ret >= 0 && pd > 0) {
                if (pd == 52) {
                    if (csv_weekly_data(dset)) {
                        ret = 52;
                    } else if (dorder == DateOrder.mmddyyyy) {
                        /* maybe we guessed wrong */
                        retransform_daily_dates(dset);
                        dorder = DateOrder.ddmmyyyy;
                        goto tryagain;
                    } else {
                        ret = -1;
                    }
                } else {
                    int convert_pd = 0;

                    if (pd == 1 || pd == 4 || pd == 12) {
                        convert_pd = pd;
                    }
                    ret = time_series_label_check(dset,
                                                  *reversed,
                                                  skipstr,
                                                  convert_pd,
                                                  prn);
                    if (ret < 0 && dorder == DateOrder.mmddyyyy) {
                        retransform_daily_dates(dset);
                        dorder = DateOrder.ddmmyyyy;
                        goto tryagain;
                    }
                }
            }
            return ret;
        }
    } else {
        pprintf(prn, _("'%s' and '%s': couldn't get dates\n").ptr, lbl1, lbl2);
    }

    return -1;
}

private int pd_from_date_label(const(char) * lbl, char * year, char * subp,
                               char * format, Prn * prn)
{
    const(char) * subchars = ".:QqMmPp-";
    int len = cast(int) strlen(lbl);
    int tryval, pd = -1; /* C: "try" is a D keyword */

    strncat(year, lbl, 4);
    tryval = atoi(year);

    if (tryval > 0 && tryval < 3000) {
        pprintf(prn, _("   %s: probably a year... ").ptr, year);
    } else {
        pprintf(prn, _("   %s: probably not a year\n").ptr, year);
    }

    if (len == 5) {
        pputs(prn, _("   but I can't make sense of the extra bit\n").ptr);
    } else if (len == 4) {
        pputs(prn, _("and just a year\n").ptr);
        pd = 1;
    } else {
        char sep = lbl[4];
        char[3] sub = 0;
        char * s = null;
        int dashQ = 0;
        int p;

        if (strchr(subchars, sep)) {
            sub[0] = '\0';
            strncat(sub.ptr, lbl + 5, 2);
            s = sub.ptr;
            if (len == 6 || (len == 7 && (sep == 'q' || sep == 'Q'))) {
                if (len == 7) s++;
                p = atoi(s);
                if (p > 0 && p < 5) {
                    pprintf(prn, _("quarter %s?\n").ptr, s);
                    pd = 4;
                } else {
                    pprintf(prn, _("quarter %d: not possible\n").ptr, p);
                }
            } else if (len == 7) {
                if (*s == 'Q') {
                    /* YYYY-Qn? This is supported by SDMX */
                    dashQ = 1;
                    s++;
                }
                p = atoi(s);
                if (dashQ) {
                    if (p > 0 && p < 5) {
                        pprintf(prn, _("quarter %d?\n").ptr, p);
                        pd = 4;
                    } else {
                        pprintf(prn, _("quarter %d: not possible\n").ptr, p);
                    }
                } else {
                    if (p > 0 && p < 13) {
                        pprintf(prn, _("month %s?\n").ptr, s);
                        pd = 12;
                    } else {
                        pprintf(prn, _("month %d: not possible\n").ptr, p);
                    }
                }
            }
            strcpy(subp, s);
            if (format != null && (pd == 4 || pd == 12)) {
                if (dashQ) {
                    sprintf(format, "%%d%cQ%%d", sep);
                } else {
                    sprintf(format, "%%d%c%%d", sep);
                }
            }
        }
    }

    return pd;
}

private int time_series_label_check(DatasetRaw * dset, int reversed,
                                    char * skipstr, int convert_pd,
                                    Prn * prn)
{
    char[5] year = 0;
    char[3] sub = 0;
    char[8] format = 0;
    char * lbl1 = dset.S[0];
    char * lbl2 = dset.S[dset.n - 1];
    char * label;
    int pd = -1;

    year[0] = sub[0] = '\0';
    label = reversed ? lbl2 : lbl1;

    if (convert_pd) {
        char[OBSLEN] altobs = 0;

        convert_daily_label(altobs.ptr, label, convert_pd);
        pd = pd_from_date_label(altobs.ptr, year.ptr, sub.ptr, format.ptr, prn);
    } else {
        pd = pd_from_date_label(label, year.ptr, sub.ptr, format.ptr, prn);
    }

    if (pd == 1) {
        if (consistent_year_labels(dset, reversed, convert_pd)) {
            dset.pd = pd;
            strcpy(dset.stobs.ptr, year.ptr);
            dset.sd0 = atof(dset.stobs.ptr);
            strcpy(dset.endobs.ptr, lbl2);
            dset.structure = TsCodes.time_series;
        } else {
            pputs(prn, _("   but the dates are not complete and consistent\n").ptr);
            pd = -1;
        }
    } else if (pd == 4 || pd == 12) {
        int savepd = pd;
        int extra_zero = 0;

        if (consistent_qm_labels(dset, reversed, convert_pd,
                                 skipstr, &pd, format.ptr,
                                 &extra_zero, prn)) {
            dset.pd = pd;
            if (savepd == 12 && pd == 4) {
                /* we switched the interpretation from
                   monthly to quarterly */
                int s;

                if (extra_zero) {
                    /* e.g. 1960Q1 written as 1960:01 */
                    s = atoi(sub.ptr + 1);
                } else {
                    /* e.g. 1960Q1 written as 1960:03 */
                    s = atoi(sub.ptr) / 3;
                }
                sprintf(dset.stobs.ptr, "%s:%d", year.ptr, s);
            } else {
                sprintf(dset.stobs.ptr, "%s:%s", year.ptr, sub.ptr);
            }
            dset.sd0 = obs_str_to_double(dset.stobs.ptr);
            ntolabel(dset.endobs.ptr, dset.n - 1, dset);
        } else {
            pputs(prn, _("   but the dates are not complete and consistent\n").ptr);
            pd = -1;
        }
    }

    return pd;
}

private int dates_maybe_reversed(const(char) * s1,
                                 const(char) * s2,
                                 Prn * prn)
{
    char[5] d1 = 0;
    char[5] d2 = 0;
    int ret = 0;

    d1[0] = d2[0] = '\0';

    strncat(d1.ptr, s1, 4);
    strncat(d2.ptr, s2, 4);

    ret = atoi(d1.ptr) > atoi(d2.ptr);

    if (ret) {
        pputs(prn, _("   dates are reversed?\n").ptr);
    }

    return ret;
}

/* e.g. "M1 1957", "M12 2009" */

private int fix_IFS_data_labels(DatasetRaw * dset)
{
    char * s1 = dset.S[0];
    char * s2 = dset.S[dset.n - 1];
    int ret = 0;

    if ((*s1 == 'M' || *s1 == 'Q') && *s2 == *s1) {
        int n1 = cast(int) strlen(s1);
        int n2 = cast(int) strlen(s2);

        if ((n1 == 7 || n1 == 8) && (n2 == 7 || n2 == 8) &&
            isdigit(s1[1]) && isdigit(s2[1])) {
            int pmax = (*s1 == 'M')? 12 : 4;
            char c;
            char[8] tmp = 0;
            char * s;
            int y, p, pbak = 0;
            int i, n, doit = 1;

            for (i=0; i<dset.n; i++) {
                s = dset.S[i];
                n = cast(int) strlen(s);
                if (n != 7 && n != 8) {
                    doit = 0;
                    break;
                }
                n = sscanf(s, "%c%d %d", &c, &p, &y);
                if (n != 3 || c != *s1) {
                    doit = 0;
                    break;
                }
                if (y < 1800 || y > 2500 || p <= 0 || p > pmax) {
                    doit = 0;
                    break;
                }
                if (i > 0 && p != pbak + 1 && p != 1) {
                    doit = 0;
                    break;
                }
                pbak = p;
            }

            if (doit) {
                for (i=0; i<dset.n; i++) {
                    s = dset.S[i];
                    sscanf(s, "%c%d %d", &c, &p, &y);
                    if (pmax == 12) {
                        sprintf(tmp.ptr, "%d:%02d", y, p);
                    } else {
                        sprintf(tmp.ptr, "%d:%d", y, p);
                    }
                    if (strlen(tmp.ptr) > strlen(s)) {
                        free(s);
                        dset.S[i] = gretl_strdup(tmp.ptr);
                    } else {
                        strcpy(s, tmp.ptr);
                    }
                }
                ret = 1;
            }
        }
    }

    return ret;
}

private int month_number(char * s)
{
    static immutable string[12] mo = [
        "jan", "feb", "mar", "apr",
        "may", "jun", "jul", "aug",
        "sep", "oct", "nov", "dec"
    ];
    int i;

    gretl_lower(s);

    for (i=0; i<12; i++) {
        if (!strcmp(s, mo[i].ptr)) {
            return i+1;
        }
    }

    return 0;
}

/* e.g. "Jan-1980", for monthly or quarterly data */

private int fix_mon_year_labels(DatasetRaw * dset)
{
    char * s1 = dset.S[0];
    char * s2 = dset.S[dset.n - 1];
    char[4] m1 = 0;
    char[4] m2 = 0;
    int yr1 = 0, yr2 = 0;
    int ret = 0;

    if (strlen(s1) == 8 && strlen(s2) == 8 &&
        s1[3] == '-' && s2[3] == '-') {
        yr1 = atoi(s1 + 4);
        yr2 = atoi(s2 + 4);
        strncat(m1.ptr, s1, 3);
        strncat(m2.ptr, s2, 3);
    }

    if (yr1 > 999 && yr1 < 3000 && yr2 > 999 && yr2 < 3000 &&
        month_number(m1.ptr) && month_number(m2.ptr)) {
        int i, p, pbak = 0;
        int dt, pd = 0;
        char * s;

        for (i=0; i<dset.n; i++) {
            s = dset.S[i];
            if (strlen(s) != 8 || s[3] != '-') {
                pd = 0;
                break;
            }
            yr1 = atoi(s + 4);
            m1[0] = '\0';
            strncat(m1.ptr, s, 3);
            if (yr1 < 1000 || yr1 >= 3000 ||
                (p = month_number(m1.ptr)) < 1) {
                pd = 0;
                break;
            }
            if (i > 0) {
                dt = p - pbak;
                if (dt != 1 && dt != 3 && p != 1) {
                    pd = 0;
                    break;
                }
                if (pd == 0 && dt > 0) {
                    pd = (dt == 1)? 12 : 4;
                }
            }
            pbak = p;
        }

        if (pd > 0) {
            for (i=0; i<dset.n; i++) {
                s = dset.S[i];
                yr1 = atoi(s + 4);
                m1[0] = '\0';
                strncat(m1.ptr, s, 3);
                p = month_number(m1.ptr);
                if (pd == 12) {
                    sprintf(dset.S[i], "%d:%02d", yr1, p);
                } else {
                    sprintf(dset.S[i], "%d:%g", yr1, ceil((3+p)/4.0));
                }
            }
            ret = 1;
        }
    }

    return ret;
}

/* Attempt to parse CSV row labels as dates.  Return -1 if this
   doesn't work out, or 0 if the labels seem to be just integer
   observation numbers, else return the inferred data frequency.
*/

int test_markers_for_dates(DatasetRaw * dset, int * reversed,
                           char * skipstr, Prn * prn)
{
    char[OBSLEN] endobs = 0;
    char * lbl1;
    char * lbl2;
    int len1, len2;
    int pd = -1;
    int n = dset.n;

    if (dset.S == null) {
        return pd;
    }

    lbl1 = dset.S[0];
    lbl2 = dset.S[n - 1];
    len1 = cast(int) strlen(lbl1);
    len2 = cast(int) strlen(lbl2);

    if (skipstr != null && *skipstr != '\0') {
        return time_series_label_check(dset, *reversed, skipstr, 0, prn);
    }

    pprintf(prn, _("   first row label \"%s\", last label \"%s\"\n").ptr,
            lbl1, lbl2);

    /* are the labels (probably) just 1, 2, 3 etc.? */
    sprintf(endobs.ptr, "%d", n);
    if (!strcmp(lbl1, "1") && !strcmp(lbl2, endobs.ptr)) {
        return 0;
    }

    if (fix_IFS_data_labels(dset) || fix_mon_year_labels(dset)) {
        lbl1 = dset.S[0];
        lbl2 = dset.S[n - 1];
        len1 = cast(int) strlen(lbl1);
    }

    /* labels are of different lengths? */
    if (len1 != len2) {
        if (abs(len1 - len2) > 1) {
            return -1;
        } else if (len2 > len1) {
            len1 = len2;
        }
    }

    pputs(prn, _("trying to parse row labels as dates...\n").ptr);

    if (len1 == 8 || len1 == 10) {
        /* daily data? */
        pd = csv_daily_date_check(dset, reversed, skipstr, prn);
    } else if (len1 >= 4) {
        /* annual, quarterly, monthly? */
        if (isdigit(cast(ubyte) lbl1[0]) &&
            isdigit(cast(ubyte) lbl1[1]) &&
            isdigit(cast(ubyte) lbl1[2]) &&
            isdigit(cast(ubyte) lbl1[3])) {
            *reversed = dates_maybe_reversed(lbl1, lbl2, prn);
            pd = time_series_label_check(dset, *reversed, skipstr, 0, prn);
        } else {
            pputs(prn, _("   definitely not a four-digit year\n").ptr);
        }
    }

    if (pd <= 0 && *reversed) {
        /* give up the "reversed" notion if we didn't get
           a workable time-series interpretation */
        *reversed = 0;
    }

    return pd;
}

private int utf8_ok(gzFile fp, int pos)
{
    c_long mark = gztell(fp);
    int len = pos + 9;
    char * test = cast(char *) malloc(len + 1);
    int i, ret = 0;

    gzseek(fp, mark - pos - 1, SEEK_SET);

    for (i=0; i<len; i++) {
        test[i] = cast(char) gzgetc(fp);
    }
    test[i] = '\0';

    if (g_utf8_validate(test, -1, null)) {
        ret = 1;
    } else {
        GError * gerr = null;
        gsize wrote = 0;
        gchar * tr;

        /* try for iso-8859? */
        tr = g_convert(test, -1, "UTF-8", "ISO-8859-15",
                       null, &wrote, &gerr);
        if (gerr != null) {
            g_error_free(gerr);
        } else {
            g_free(tr);
            ret = 1;
        }
    }

    free(test);

    gzseek(fp, mark, SEEK_SET);

    return ret;
}

private enum UnicodeBOM {
    utf_8 = 1,
    utf_16,
    utf_32
}

/* If we got a UTF-16 or UTF-32 BOM, try recoding to UTF-8 before
   parsing data. We write the recoded text to a temporary file in the
   user's "dotdir" (and then delete that file once we're done).
*/

private int csv_recode_input(gzFile * fpp,
                             const(char) * fname,
                             gchar ** pfname,
                             int ucode,
                             Prn * prn)
{
    const(gchar) * from_set =
        (ucode == UnicodeBOM.utf_32)? "UTF-32" : "UTF-16";
    gchar * altname = null;
    int err = 0;

    /* the current stream is not useable as is,
       so shut it down
    */
    gzclose(*fpp);
    *fpp = null;

    /* we'll recode to a temp file in dotdir */
    altname = g_strdup_printf("%srecode_tmp.u8", gretl_dotdir());

    err = gretl_recode_file(fname, altname,
                            from_set, "UTF-8",
                            prn);

    if (!err) {
        /* try reattaching the stream */
        *fpp = gretl_gzopen(altname, "rb");
        if (*fpp == null) {
            gretl_remove(altname);
            err = GretlError.fopen;
        } else {
            pputs(prn, "switched to recoded input\n");
            *pfname = altname;
            altname = null;
        }
    }

    g_free(altname);

    return err;
}

/* Check the first 4 bytes of "CSV" input for a Byte Order Mark. If we
   find the UTF-8 BOM (typically written by Microsoft tools), simply
   record the fact so that we can skip it on reading. But if we find a
   BOM indicating a 16-bit or 32-bit unicode encoding, flag this by
   returning a non-zero @ucode value; in that case we'll attempt a full
   recording of the input (via GLib) before we start reading data.
*/

private int csv_unicode_check(gzFile fp, CsvData * c, Prn * prn)
{
    ubyte[4] b = 0;
    int n = gzread(fp, b.ptr, 4);
    int ucode = 0;

    if (n == 4) {
        if (b[0] == 0xEF && b[1] == 0xBB && b[2] == 0xBF) {
            pputs(prn, "got UTF-8 BOM\n");
            ucode = UnicodeBOM.utf_8;
        } else if (b[0] == 0xFE && b[1] == 0xFF) {
            pputs(prn, "got UTF-16BE, will try recoding\n");
            ucode = UnicodeBOM.utf_16;
        } else if (b[0] == 0xFF && b[1] == 0xFE) {
            if (b[2] == 0 && b[3] == 0) {
                pputs(prn, "got UTF-32LE, will try recoding\n");
                ucode = UnicodeBOM.utf_32;
            } else {
                pputs(prn, "got UTF-16LE, will try recoding\n");
                ucode = UnicodeBOM.utf_16;
            }
        } else if (b[0] == 0 && b[1] == 0 &&
                   b[0] == 0xFE && b[1] == 0xFF) {
            pputs(prn, "got UTF-32BE, will try recoding\n");
            ucode = UnicodeBOM.utf_32;
        }
    }

    if (ucode == UnicodeBOM.utf_8) {
        csv_set_has_bom(c);
        gzseek(fp, 3, SEEK_SET);
        ucode = 0;
    } else {
        gzrewind(fp);
    }

    return ucode;
}

private int report_quotation_info(int * qtotal, int * nqmax,
                                  char qchar, Prn * prn)
{
    int err = 0;

    if (qtotal[0] > 0) {
        pprintf(prn, _("Found %d double-quotes, max %d per line\n").ptr,
                qtotal[0], nqmax[0]);
    }
    if (qtotal[1] > 0) {
        pprintf(prn, _("Found %d single-quotes, max %d per line\n").ptr,
                qtotal[1], nqmax[1]);
    }
    if (qchar == '"') {
        pputs(prn, (_("Assuming double-quote is the relevant ")
                    ~ _("quotation character\n")).ptr);
    } else if (qchar == '\'') {
        pputs(prn, (_("Assuming single-quote is the relevant ")
                    ~ _("quotation character\n")).ptr);
    } else {
        gretl_warnmsg_set(_("unmatched quotation marks in CSV file").ptr);
        pputs(prn, _("Quotation broken: marks not paired\n").ptr);
    }

    return err;
}

private int handle_CR(int * crlf, gzFile fp)
{
    int c, c1 = gzgetc(fp);

    if (c1 == EOF) {
        return c1;
    } else if (c1 == 0x0a) {
        /* CR + LF -> LF */
        *crlf = 1;
        c = c1;
    } else {
        /* old Mac-style: CR not followed by LF */
        c = 0x0a;
        gzungetc(c1, fp);
    }

    return c;
}

/* The function below checks for the maximum line length in the given
   file.  It also checks for extraneous binary data (the file is
   supposed to be plain text), and checks whether the 'delim'
   character is present in the file, on a non-comment line (where
   a comment line is one that starts with '#').

   In addition, we check whether the file has a trailing comma on every
   line, and for the numbers of double- and single-quote characters
   to try to determine which, if either, is used to indicate quoted
   fields in the input.

   We return the maximum line length, or -1 on error.
*/

private int csv_max_line_length(gzFile fp, CsvData * cdata, Prn * prn)
{
    int c, c1, cbak = 0, cc = 0;
    int comment = 0, maxlinelen = 0;
    int[2] qcand = [1, 1];
    int[2] nq = 0;
    int[2] nqmax = 0;
    int[2] qtotal = 0;
    int dquoted = 0;
    int crlf = 0, lines = 0;
    int i, err = 0;

    csv_set_trailing_comma(cdata); /* just provisionally */

    while ((c = gzgetc(fp)) != EOF) {
        if (c == 0x0d) {
            c = handle_CR(&crlf, fp);
            if (c == EOF) {
                break;
            }
        }
        if (c == 0x0a) {
            /* we reached the end of a line */
            if (cc > maxlinelen) {
                maxlinelen = cc;
            }
            cc = 0;
            if (cbak != 0 && cbak != ',') {
                csv_unset_trailing_comma(cdata);
            }
            lines++;
            for (i=0; i<2; i++) {
                /* quotation--mark info, double and single */
                if (nq[i] % 2) {
                    /* unpaired: not a valid candidate */
                    qcand[i] = 0;
                }
                if (nq[i] > nqmax[i]) {
                    nqmax[i] = nq[i];
                }
                qtotal[i] += nq[i];
                nq[i] = 0; /* reset count */
            }
            dquoted = 0;
            continue;
        }
        cbak = c;
        if (!isspace(cast(ubyte) c) && !isprint(cast(ubyte) c) &&
            !(c == CTRLZ) && !utf8_ok(fp, cc)) {
            pprintf(prn, (_("Binary data (%d) encountered (line %d:%d): ")
                          ~ _("this is not a valid text file\n")).ptr,
                    c, lines + 1, cc + 1);
            return -1;
        }
        if (cc == 0) {
            /* line starts with hash mark */
            comment = (c == '#');
        }
        if (!comment) {
            if (c == '\t') {
                /* let's ignore trailing tabs in this heuristic */
                c1 = gzgetc(fp);
                if (c1 != 0x0d && c1 != 0x0a) {
                    csv_set_got_tab(cdata);
                }
                gzungetc(c1, fp);
            }
            if (c == ';') {
                csv_set_got_semi(cdata);
            } else if (c == '|') {
                csv_set_got_pipe(cdata);
            }
            if (c == cdata.delim && !dquoted) {
                csv_set_got_delim(cdata);
            } else if (c == '"') {
                dquoted = !dquoted;
                nq[0] += 1;
            } else if (c == '\'') {
                nq[1] += 1;
            }
        }
        cc++;
    }

    if (maxlinelen == 0) {
        pputs(prn, _("Data file is empty\n").ptr);
    } else if (csv_has_trailing_comma(cdata)) {
        pputs(prn, _("Data file has trailing commas\n").ptr);
    }

    for (i=0; i<2; i++) {
        if (qcand[i] && qtotal[i] == 0) {
            qcand[i] = 0;
        }
    }

    if (qtotal[0] || qtotal[1]) {
        if (qcand[0] && qcand[1]) {
            /* quotes: double and single are both paired,
               maybe prefer the more numerous as @qchar?
            */
            cdata.qchar = qtotal[1] > qtotal[0] ? '\'' : '"';
        } else if (qcand[0]) {
            cdata.qchar = '"';
        } else if (qcand[1]) {
            cdata.qchar = '\'';
        }
        err = report_quotation_info(qtotal.ptr, nqmax.ptr, cdata.qchar, prn);
    }

    if (!err && maxlinelen > 0) {
        /* allow for newline and null terminator */
        maxlinelen += 2 + crlf;
    }

    return err ? -1 : maxlinelen;
}

private bool nonspace_delim(char d)
{
    return d != ',' && d != ';' && d != '\t' && d != '|';
}

private int count_csv_fields(CsvData * c)
{
    const(char) * s = c.line;
    int inquote = 0;
    int cbak, nf = 0;

    if (*s == c.delim && *s == ' ') {
        s++;
    }

    while (*s) {
        if (csv_keep_quotes(c) && *s == c.qchar) {
            inquote = !inquote;
        } else if (!inquote && *s == c.delim) {
            nf++;
        }
        cbak = *s;
        s++;
        /* Problem: (when) should a trailing delimiter be read as an
           implicit NA?  For now we'll so treat it if the delimiter
           is not plain space.
        */
        if (*s == '\0' && cbak == c.delim && nonspace_delim(c.delim)) {
            nf--;
        }
    }

    return nf + 1;
}

private void purge_quoted_delimiters(CsvData * c)
{
    char * s = c.line;
    int inquote = 0;

    while (*s) {
        if (*s == '"') {
            inquote = !inquote;
        } else if (inquote && *s == c.delim) {
            *s = (c.delim == ' ')? '_' : ' ';
        }
        s++;
    }
}

private void purge_unquoted_spaces(char * s)
{
    int inquote = 0;

    while (*s) {
        if (*s == '"') {
            inquote = !inquote;
        } else if (!inquote && *s == ' ') {
            shift_string_left(s, 1);
        }
        s++;
    }
}

private void compress_csv_line(CsvData * c, int nospace)
{
    int n = cast(int) strlen(c.line);
    char * p = c.line + n - 1;

    if (*p == 0x0a) {
        *p = '\0';
        p--;
    }

    if (*p == 0x0d) {
        *p = '\0';
    }

    if (!csv_keep_quotes(c)) {
        purge_quoted_delimiters(c);
    }

    if (c.delim != ' ') {
        if (nospace) {
            purge_unquoted_spaces(c.line);
        }
    } else {
        compress_spaces(c.line);
    }

    if (!csv_keep_quotes(c)) {
        gretl_delchar('"', c.line);
    }

    if (csv_has_trailing_comma(c)) {
        /* chop trailing comma */
        n = cast(int) strlen(c.line);
        if (n > 0) {
            c.line[n-1] = '\0';
        }
    }
}

int import_obs_label(const(char) * s)
{
    char[VNAMELEN] tmp = 0;

    if (s == null) {
        return 1;
    }

    if (!strcmp(s, "\"\"") || !strcmp(s, "''")) {
        return 1;
    }

    if (*s == '"' || *s == '\'') s++;

    if (*s == '\0') {
        return 1;
    }

    if (strlen(s) > VNAMELEN - 1) {
        return 0;
    }

    tmp[0] = '\0';
    strncat(tmp.ptr, s, VNAMELEN - 1);
    gretl_lower(tmp.ptr);

    return (!strcmp(tmp.ptr, "obs") ||
            !strcmp(tmp.ptr, "date") ||
            !strcmp(tmp.ptr, "year") ||
            !strcmp(tmp.ptr, "period") ||
            !strcmp(tmp.ptr, "observation") ||
            !strcmp(tmp.ptr, "observation_date"));
}

private int join_wants_col_zero(CsvData * c, const(char) * s)
{
    const(char) * colname;
    int i;

    if (*s == '\0') {
        return 0;
    }

    for (i=0; i<c.jspec.ncols; i++) {
        colname = c.jspec.colnames[i];
        if (colname != null && !strcmp(s, colname)) {
            return 1;
        }
    }

    return 0;
}

private void check_first_field(const(char) * line, CsvData * c, Prn * prn)
{
    const(char) * s;

 tryagain:
    s = line;

    if (c.delim != ' ' && *s == c.delim) {
        csv_set_blank_column(c);
    } else {
        char[VNAMELEN] field1 = 0;
        int quoted = 0;
        int i = 0;

        if (c.delim == ' ' && *s == ' ') {
            s++;
        }

        if (*s == c.qchar) {
            quoted = 1;
        }

        while (*s && i < field1.sizeof - 1) {
            if (!quoted && *s == c.delim) {
                break;
            } else if (!quoted && *s == '\t') {
                /* presence of a tab must indicate tab-separation? */
                c.delim = '\t';
                goto tryagain;
            } else if (i > 0 && *s == c.qchar) {
                quoted = 0;
            }
            field1[i++] = *s++;
        }

        field1[i] = '\0';
        iso_to_ascii(field1.ptr);

        if (joining(c) && join_wants_col_zero(c, field1.ptr)) {
            return;
        } else if (csv_all_cols(c)) {
            /* open/append wants all columns as data */
            return;
        }

        pprintf(prn, _("   first field: '%s'\n").ptr, field1.ptr);

        if (import_obs_label(field1.ptr)) {
            pputs(prn, _("   seems to be observation label\n").ptr);
            csv_set_obs_column(c);
        }
    }
}

void import_na_init()
{
    const(char) * s = get_csv_na_read_string();

    import_na[0] = '\0';
    strncat(import_na.ptr, s, import_na.sizeof - 1);
}

/* Returns 1 if the string @s should be counted as representing a
   missing value, 0 otherwise. If there is a user-set "csv_read_na"
   value we consult it, otherwise we consult a set of default values.
*/

int import_na_string(const(char) * s)
{
    static immutable string[14] defaults = [
        "NA",
        "N.A.",
        "n.a.",
        "na",
        "n/a",
        "N/A",
        "#N/A",
        "NaN",
        ".NaN",
        ".",
        "..",
        "-999",
        "-9999",
        "-"
    ];
    int i;

    if (import_na[0] != '\0' && strcmp(import_na.ptr, "default")) {
        /* the user has set a specific "NA" string */
        return !strcmp(s, import_na.ptr);
    }

    /* Should we continue in the following way if "csv_read_na"
       is user-specified? */
    for (i=0; i<defaults.length; i++) {
        if (!strcmp(s, defaults[i].ptr)) {
            return 1;
        }
    }

    return 0;
}

private int csv_missval(const(char) * str, int i, int t,
                        int * miss_shown, Prn * prn)
{
    int miss = 0;

    if (*str == '\0' || !strcmp(str, "\"\"")) {
        /* 2021-03-03: let '""' indicate missing */
        if (miss_shown != null) {
            if (t < 80 || *miss_shown < i) {
                pprintf(prn, (_("   the cell for variable %d, obs %d ")
                              ~ _("is empty: treating as missing value\n")).ptr,
                        i, t);
                *miss_shown += 1;
            }
        }
        miss = 1;
    }

    if (import_na_string(str)) {
        if (miss_shown != null) {
            if (t < 80 || *miss_shown < i) {
                pprintf(prn, (_("   warning: missing value for variable ")
                              ~ _("%d, obs %d\n")).ptr, i, t);
                *miss_shown += 1;
            }
        }
        miss = 1;
    }

    return miss;
}

/* In the case where we think we've found thousands separators in
   numerical input, provisionally mark all "non-numeric" values as NAs;
   we do this prior to a second pass through the data.
*/

private void revise_non_numeric_values(CsvData * c)
{
    int i, t;

    for (i=1; i<c.dset.v; i++) {
        for (t=0; t<c.dset.n; t++) {
            if (c.dset.Z[i][t] == NON_NUMERIC) {
                c.dset.Z[i][t] = NADBL;
            }
        }
    }
}

int non_numeric_check(DatasetRaw * dset, int ** plist,
                      GretlStringTable ** pst,
                      Prn * prn)
{
    int * list = null;
    int i, j, t, nn = 0;
    int err = 0;

    if (pst == null) {
        /* not interested in string-valued series/columns */
        for (i=1; i<dset.v; i++) {
            for (t=0; t<dset.n; t++) {
                if (dset.Z[i][t] == NON_NUMERIC) {
                    dset.Z[i][t] = NADBL;
                }
            }
        }
        return 0;
    }

    for (i=1; i<dset.v; i++) {
        for (t=0; t<dset.n; t++) {
            if (dset.Z[i][t] == NON_NUMERIC) {
                nn++;
                break;
            }
        }
    }

    if (nn == 0) {
        return 0; /* nothing to be done */
    }

    list = gretl_list_new(nn);
    if (list == null) {
        return GretlError.alloc;
    }

    j = 1;
    for (i=1; i<dset.v; i++) {
        for (t=0; t<dset.n; t++) {
            if (dset.Z[i][t] == NON_NUMERIC) {
                list[j++] = i;
                break;
            }
        }
    }

    for (i=1; i<=list[0]; i++) {
        /* check each member of @list */
        double nnfrac;
        int nnon = 0;
        int tnon = -1;
        int nok = 0;
        int v = list[i];

        series_set_flag(dset, v, VarFlags.discrete);

        for (t=0; t<dset.n; t++) {
            if (dset.Z[v][t] == NON_NUMERIC) {
                if (tnon < 0) {
                    /* record the first non-numeric obs */
                    tnon = t + 1;
                }
                nnon++;
            } else if (!na(dset.Z[v][t])) {
                nok++;
            }
        }

        nnfrac = (nok == 0)? 1.0 : nnon / cast(double) (nnon + nok);
        pprintf(prn, (_("variable %d (%s): non-numeric values = %d ")
                      ~ _("(%.2f percent)\n")).ptr, v, dset.varname[v],
                nnon, 100 * nnfrac);
        if ((nnon < 2 && dset.n > 2) || nnfrac < 0.05) {
            /* if we got just a few non-numeric values, we'll assume
               that the data file is broken
            */
            pprintf(prn, (_("ERROR: variable %d (%s), observation %d, ")
                          ~ _("expected numeric value\n")).ptr,
                    v, dset.varname[v], tnon);
            err = GretlError.data;
            break;
        }
    }

    if (!err) {
        pputs(prn, _("allocating string table\n").ptr);
        *pst = gretl_string_table_new(list);
        if (*pst == null) {
            err = GretlError.alloc;
        }
    }

    if (err) {
        free(list);
    } else {
        *plist = list;
    }

    return err;
}

private int csv_non_numeric_check(CsvData * c, Prn * prn)
{
    GretlStringTable * st = null;
    int * nlist = null;
    int err = 0;

    if (csv_as_matrix(c)) {
        err = non_numeric_check(c.dset, &nlist, null, prn);
    } else {
        err = non_numeric_check(c.dset, &nlist, &st, prn);
    }

    if (!err) {
        c.codelist = nlist;
        c.st = st;
    }

    return err;
}

/* Handle the case in "join" where the user specified some time
   columns for conversion to numeric and also gave a specific format
   for the conversion.
*/

private double special_time_val(const(char) * s, const(char) * fmt,
                                int m_means_q)
{
    tm t;
    char * test;

    /* C: struct tm t = {0,0,0,1,0,0,0,0,-1} */
    t.tm_mday = 1;
    t.tm_isdst = -1;

    test = strptime(s, fmt, &t);

    if (test == null || *test != '\0') {
        /* conversion didn't work right */
        return NADBL;
    } else {
        int y, m, d;

        y = t.tm_year + 1900;
        m = t.tm_mon + 1;
        d = t.tm_mday;

        if (m_means_q) {
            /* convert to 1st month of quarter */
            if (m == 2) m = 4;
            else if (m == 3) m = 7;
            else if (m == 4) m = 10;
            else if (m != 1) {
                return NADBL;
            }
        }

        if (d == 0) d = 1;

        return 10000*y + 100*m + d;
    }
}

private int char_count(char c, const(char) * s)
{
    int n = 0;

    while (*s) {
        if (*s == c) n++;
        s++;
    }

    return n;
}

/* Follow-up check for the case where we think we might
   have found a thousands separator: each occurrence of
   the putative separator must be followed by exactly 3
   digits: we set c->thousep to an invalid value if this
   is not the case.
*/

private void validate_thousep(CsvData * c, const(char) * s)
{
    int nd;

    while (*s) {
        if (*s == c.thousep) {
            nd = 0;
            s++;
            while (*s) {
                if (isdigit(*s)) {
                    nd++;
                    s++;
                } else {
                    break;
                }
            }
            if (nd != 3) {
                /* nope! */
                c.thousep = -1;
                break;
            }
        } else {
            s++;
        }
    }
}

/* Initial heuristic for detecting a thousands separator,
   where the string @s has been determined to contain
   nothing but digits, dot and comma (allowing for a leading
   minus).

   1) If the string contains both comma and dot, whichever
   character appears to the left cannot be the decimal
   separator and may be a thousands separator.

   2) If more than one comma appears in the string, comma
   cannot be the decimal character and might be a thousands
   separator; mutatis mutandis for dot.
*/

private void test_for_thousands_sep(CsvData * c, const(char) * s)
{
    const(char) * p1 = strrchr(s, '.');
    const(char) * p2 = strrchr(s, ',');
    char thousep = 0;

    if (p1 != null && p2 != null) {
        thousep = (p2 - p1 > 0)? '.' : ',';
    } else if (p1 != null && char_count('.', s) > 0) {
        thousep = '.';
    } else if (p2 != null && char_count(',', s) > 0) {
        thousep = ',';
    }

    if (c.thousep > 0) {
        if (thousep != 0 && thousep != c.thousep) {
            /* no consistent interpretation exists */
            c.thousep = -1; /* invalid */
        }
    } else if (thousep != 0) {
        /* we have a candidate for testing */
        const(char) * test;
        char[CSVSTRLEN] tmp = 0;

        strcpy(tmp.ptr, s);
        gretl_delchar(thousep, tmp.ptr);
        if (thousep == '.' && get_local_decpoint() == '.') {
            gretl_charsub(tmp.ptr, ',', '.');
        }
        errno = 0;
        strtod(tmp.ptr, &test);
        if (*test == '\0' && errno == 0) {
            c.thousep = cast(byte) thousep;
        }
    }

    if (c.thousep && thousep != 0) {
        validate_thousep(c, s);
    }
}

private int all_digits_and_seps(const(char) * s)
{
    const(char) * test = "0123456789.,";

    if (*s == '-') s++;

    return strspn(s, test) == strlen(s);
}

private double eval_non_numeric(CsvData * c, int i, const(char) * s)
{
    double x = NON_NUMERIC;

    if (series_get_flags(c.dset, i) & VarFlags.timecol) {
        char * fmt = null;
        int mq = 0;

        if (timecol_get_format(c.dset, i, &fmt, &mq)) {
            /* the user gave a specific format for this */
            x = special_time_val(s, fmt, mq);
        } else {
            /* default: ISO 8601 extended */
            int y, m, d, n;

            n = sscanf(s, "%d-%d-%d", &y, &m, &d);
            if (n == 3) {
                x = 10000*y + 100*m + d;
            } else {
                x = NADBL;
            }
        }
    } else if (c.thousep >= 0 && !csv_scrub_thousep(c)) {
        /* Here we consider the possibility although @s does not
           validate as numeric according to the C library, it is by
           intent numeric but includes one or more thousands
           separators.

           The condition c->thousep >= 0 requires that we haven't
           already ruled out this interpretation due to inconsistency,
           and !csv_scrub_thousep(c) requires that we're not on a
           second pass through the data.
        */
        if (all_digits_and_seps(s)) {
            test_for_thousands_sep(c, s);
        }
    }

    return x;
}

private int converted_ok(const(char) * s, const(char) * test, double x)
{
    if (*test != '\0') {
        if (errno) perror(s);
        return 0; /* definitely not OK */
    } else if (errno == ERANGE && fabs(x) > 0 && fabs(x) < 0.001) {
        return 1; /* subnormal, but we'll let that pass */
    } else if (errno) {
        perror(s);
        return 0;
    } else {
        return 1;
    }
}

private char * csv_unquote(char * s)
{
    if (s[0] == '"') {
        int i, n = cast(int) strlen(s);

        if (n > 1 && s[n-1] == '"') {
            for (i=0; i<n-2; i++) {
                s[i] = s[i+1];
            }
            s[i] = '\0';
        }
    }
    return s;
}

private double csv_atof(CsvData * c, int i)
{
    char[CSVSTRLEN] tmp = 0;
    char[CSVSTRLEN] clean = 0;
    double x = NON_NUMERIC;
    const(char) * s = c.str.ptr;
    const(char) * test;

    if (csv_scrub_thousep(c) && strchr(s, c.thousep) &&
        all_digits_and_seps(s)) {
        /* second pass through the data: pre-process fields
           that we reckon include thousands separators
        */
        strcpy(clean.ptr, s);
        gretl_delchar(c.thousep, clean.ptr);
        s = clean.ptr;
    }

    if (c.decpoint == '.' || !csv_do_dotsub(c) || strchr(s, ',') == null) {
        /* either we're currently set to the correct locale,
           or there's no problematic decimal point in @s
        */
        errno = 0;
        x = strtod(s, &test);
        if (converted_ok(s, test, x)) {
            return x; /* handled */
        }
    } else if (csv_do_dotsub(c)) {
        /* in C numeric locale: substitute dot for comma */
        strcpy(tmp.ptr, s);
        gretl_charsub(tmp.ptr, ',', '.');
        errno = 0;
        x = strtod(tmp.ptr, &test);
        if (converted_ok(s, test, x)) {
            return x; /* handled */
        }
    }

    if (c.decpoint == '.' && strchr(s, ',') != null) {
        /* try remediation for decimal comma? */
        strcpy(tmp.ptr, s);
        gretl_charsub(tmp.ptr, ',', '.');
        errno = 0;
        x = strtod(tmp.ptr, &test);
        if (converted_ok(s, test, x)) {
            return x; /* handled */
        }
    }

    /* fallback */
    /* revised 2020-02-13 to use csv_unquote */
    return eval_non_numeric(c, i, csv_unquote(c.str.ptr));
}

private int process_csv_obs(CsvData * c, int i, int t, int * miss_shown,
                            Prn * prn)
{
    int err = 0;

    if (c.st != null) {
        /* second round, handling string-valued variables */
        if (in_gretl_list(c.codelist, i)) {
            double zit = c.dset.Z[i][t];
            int ix;

            if (na(zit) && c.str[0] != '\0' && c.user_na == null) {
                /* by default (no user_na) only blanks count as NAs */
                zit = NON_NUMERIC;
            }
            if (!na(zit)) {
                ix = gretl_string_table_index(c.st, c.str.ptr, i, 0, prn);
                if (ix > 0) {
                    c.dset.Z[i][t] = cast(double) ix;
                } else {
                    err = GretlError.data;
                }
            }
        }
    } else if (csv_missval(c.str.ptr, i, t+1, miss_shown, prn)) {
        c.dset.Z[i][t] = NADBL;
    } else {
        gretl_strstrip(c.str.ptr);
        c.dset.Z[i][t] = csv_atof(c, i);
    }

    return err;
}

/* Emulation of fgets(), designed to handle any sort of line
   termination (unix, DOS, Mac or even an unholy mixture).
   Line-endings are converted to LF (0x0a).
*/

private char * csv_fgets(CsvData * cdata, gzFile fp)
{
    char * s = cdata.line;
    int n = cdata.maxlinelen;
    int i, c1, c = 0;

    for (i=0; i<n-1 && c!=0x0a; i++) {
        c = gzgetc(fp);
        if (c == EOF) {
            if (i == 0) {
                /* signal end of read */
                return null;
            } else {
                break;
            }
        } else if (c == 0x0d) {
            /* CR: convert to LF and peek at next char: if it's
               LF swallow it, otherwise put it back */
            c = 0x0a;
            c1 = gzgetc(fp);
            if (c1 != 0x0a) {
                gzungetc(c1, fp);
            }
        }
        s[i] = cast(char) c;
    }

    s[i] = '\0';

    return s;
}

/* pick up any comments following the data block in a CSV file */

private char * get_csv_descrip(CsvData * c, gzFile fp)
{
    char * line = c.line;
    char * desc = null;
    size_t llen, totlen;

    while (csv_fgets(c, fp)) {
        tailstrip(line);
        llen = strlen(line);
        if (desc == null) {
            totlen = llen + 4;
            desc = cast(char *) malloc(totlen);
            if (desc == null) {
                return null;
            }
            sprintf(desc, "%s\n", line);
        } else {
            char * tmp;

            totlen = strlen(desc) + llen + 4;
            tmp = cast(char *) realloc(desc, totlen);
            if (tmp == null) {
                free(desc);
                return null;
            }
            desc = tmp;
            strcat(desc, line);
            strcat(desc, "\n");
        }
    }

    if (desc != null && string_is_blank(desc)) {
        free(desc);
        desc = null;
    }

    return desc;
}

private immutable string csv_msg = N_("\nPlease note:\n"
             ~ "- The first row of the CSV file should contain the "
             ~ "names of the variables.\n"
             ~ "- The first column may optionally contain date "
             ~ "strings or other 'markers':\n  in that case its row 1 entry "
             ~ "should be blank, or should say 'obs' or 'date'.\n"
             ~ "- The remainder of the file must be a rectangular "
             ~ "array of data.\n");

/* Here we check whether we get a consistent reading on the number of
   fields per line in the CSV file.
*/

private int csv_fields_check(gzFile fp, CsvData * c, Prn * prn)
{
    int gotdata = 0;
    int chkcols = 0;
    int err = 0;

    c.ncols = c.nrows = 0;

    if (csv_has_bom(c)) {
        gzseek(fp, 3, SEEK_SET);
    }

    while (csv_fgets(c, fp) && !err) {
        /* skip comment lines */
        if (*c.line == '#') {
            continue;
        }
        /* skip blank lines -- but finish if the blank comes after data */
        if (string_is_blank(c.line)) {
            if (gotdata) {
                if (!csv_have_data(c)) {
                    c.descrip = get_csv_descrip(c, fp);
                }
                break;
            } else {
                continue;
            }
        }

        c.nrows += 1;

        if (fixed_format(c)) {
            tailstrip(c.line);
            gotdata = 1;
            chkcols = cast(int) strlen(c.line);
            if (chkcols < c.cols_list[c.cols_list[0]]) {
                gretl_errmsg_set(_("Invalid column specification").ptr);
                err = GretlError.data;
                break;
            } else {
                continue;
            }
        }

        compress_csv_line(c, 1);

        if (!gotdata) {
            /* scrutinize the first "real" line */
            check_first_field(c.line, c, prn);
            gotdata = 1;
        }

        chkcols = count_csv_fields(c);
        if (c.ncols == 0) {
            c.ncols = chkcols;
            pprintf(prn, _("   number of columns = %d\n").ptr, c.ncols);
        } else if (chkcols > c.ncols) {
            pprintf(prn, _("   ...but row %d has %d fields: aborting\n").ptr,
                    c.nrows, chkcols);
            err = GretlError.data;
        } else if (chkcols < c.ncols) {
            pprintf(prn, _("Warning: row %d has %d fields: %d entries set to missing\n").ptr,
                    c.nrows, chkcols, c.ncols - chkcols);
        } else if (cols_subset(c)) {
            int datacols = csv_skip_col_1(c) ? (c.ncols - 1) : c.ncols;

            if (c.cols_list[c.cols_list[0]] > datacols) {
                gretl_errmsg_set(_("Invalid column specification").ptr);
                err = GretlError.data;
            }
        }
    }

    if (!err && fixed_format(c)) {
        c.ncols = c.cols_list[0];
    }

    return err;
}

private void strip_illegals(char * s)
{
    char[VNAMELEN] name = 0;
    int i, j = 0;

    for (i=0; s[i] != '\0'; i++) {
        if (isalnum(s[i]) || s[i] == '_') {
            name[j++] = s[i];
        }
    }

    name[j] = '\0';
    strcpy(s, name.ptr);
}

private int intercept_nan_as_name(const(char) * s)
{
    if (strlen(s) == 3) {
        char[4] screen = 0;

        strcpy(screen.ptr, s);
        gretl_lower(screen.ptr);
        if (!strcmp(screen.ptr, "nan")) {
            return 1;
        }
    }

    return 0;
}

private int csv_is_numeric(const(char) * s, CsvData * c)
{
    int ret = 0;

    if (!strcmp(s, ".")) {
        /* treat as NA? */
        ret = 1;
    } else if (c.decpoint == '.') {
        ret = numeric_string(s);
    } else {
        /* decimal comma in force */
        char * tmp = gretl_strdup(s);

        gretl_charsub(tmp, ',', '.');
        ret = numeric_string(tmp);
        free(tmp);
    }

    return ret;
}

private int process_csv_varname(CsvData * c, int j, int * numcount,
                                Prn * prn)
{
    char * vname = c.dset.varname[j];
    char * src = c.str.ptr;
    int err = 0;

    *vname = '\0';

    if (intercept_nan_as_name(src)) {
        gretl_errmsg_sprintf((_("If '%s' is intended as the name of a variable, ")
                              ~ _("please change it --\nstrings of this sort usually ")
                              ~ _("mean 'not a number'.")).ptr, src);
        err = GretlError.data;
    } else if (*src == '\0') {
        fprintf(stderr, "variable name %d is missing\n", j);
        sprintf(vname, "v%d", j);
    } else if (csv_is_numeric(src, c)) {
        *numcount += 1;
    } else {
        const(char) * s = src;

        while (*s && !isalpha(*s)) s++;
        if (*s == '\0') {
            fprintf(stderr, "variable name %d (%s) is garbage\n", j, src);
            sprintf(vname, "v%d", j);
        } else {
            strncat(vname, s, VNAMELEN - 1);
        }
        iso_to_ascii(vname);
        strip_illegals(vname);
        if (gretl_reserved_word(vname)) {
            /* try a fix for this */
            int n = cast(int) strlen(vname);

            if (n < VNAMELEN-1) {
                strcat(vname, "_");
            } else {
                vname[n-1] = '_';
            }
        }
        if (check_varname(vname)) {
            errmsg(1, prn);
            err = GretlError.data;
        }
    }

    return err;
}

private int csv_reconfigure_for_markers(DatasetRaw * dset)
{
    int err = dataset_allocate_obs_markers(dset);

    if (!err) {
        err = dataset_drop_last_variables(dset, 1);
    }

    return err;
}

private int skip_data_column(CsvData * c, int k)
{
    int col = csv_skip_col_1(c) ? k : k + 1;

    if (!in_gretl_list(c.cols_list, col)) {
        return 1;
    } else {
        return 0;
    }
}

private int update_join_cols_list(CsvData * c, int k)
{
    int * test;
    int err = 0;

    test = gretl_list_append_term(&c.cols_list, k);
    if (test == null) {
        err = GretlError.alloc;
    }

    return err;
}

/* handle_join_varname: the index @k contains the column number
   relative to the entire CSV file, while @pj points to j, the column
   number relative to the reduced dataset that will be constructed by
   selection of columns from the file.

   Here we're examining a column heading read from file (c->str) to
   see whether it matches any of the column-names required for an
   ongoing join operation (held in c->jspec->colnames). If so, we
   write the index j into the appropriate slot in c->jspec->colnums
   (which starts off filled with zeros), so the joiner will know where
   to find the required data. (The j value is bound to be at least 1
   since column 0 is reserved to the constant.)

   In some cases a given named column may perform more than one role in
   a join operation -- for example, it may serve as an element in a
   filter and also as the auxiliary variable in an "aggregation"
   method. To allow for this we don't stop scanning at the first match
   of c->str with a required column name.

   The call to update_join_cols_list() uses the index @k to record the
   overall column position of "wanted data", for use by the CSV
   reader.
*/

private int handle_join_varname(CsvData * c, int k, int * pj)
{
    const(char) * colname;
    char[VNAMELEN] okname = 0;
    int matched = 0;
    int i, j = *pj;

    if (!csv_skip_col_1(c)) {
        k++;
    }

    if (csv_no_header(c)) {
        sprintf(okname.ptr, "col%d", k);
    } else {
        /* convert to valid gretl identifier */
        gretl_normalize_varname(okname.ptr, c.str.ptr, 0, k);
    }

    for (i=0; i<c.jspec.ncols; i++) {
        /* find "wanted name" i */
        colname = c.jspec.colnames[i];
        if (colname == null || c.jspec.colnums[i] > 0) {
            /* name not wanted, or already found */
            continue;
        }
        if (!strcmp(okname.ptr, colname)) {
            c.jspec.colnums[i] = j;
            if (!matched) {
                matched = 1;
                strcpy(c.dset.varname[j], okname.ptr);
                update_join_cols_list(c, k);
                *pj += 1;
                if (in_gretl_list(c.jspec.timecols, i)) {
                    series_set_flag(c.dset, j, VarFlags.timecol);
                }
            }
        }
    }

    return 0;
}

private bool starts_number(int c)
{
    return isdigit(cast(ubyte) c) || c == '-' || c == '+' || c == '.';
}

private bool obs_labels_no_varnames(int o, const(DatasetRaw) * c, int n)
{
    return !o && c.v > 3 && n == c.v - 2;
}

private int csv_varname_scan(CsvData * c, gzFile fp, Prn * prn, Prn * mprn)
{
    char * p;
    int obscol = csv_has_obs_column(c);
    int i, j, k, numcount;
    int err = 0;

    if (!csv_no_header(c)) {
        pputs(mprn, _("scanning for variable names...\n").ptr);
    }

    if (csv_has_bom(c)) {
        gzseek(fp, 3, SEEK_SET);
    }

    while (csv_fgets(c, fp)) {
        if (*c.line == '#' || string_is_blank(c.line)) {
            continue;
        } else {
            break;
        }
    }

    c.datapos = gztell(fp);

    compress_csv_line(c, 1);

    p = c.line;
    if (c.delim == ' ' && *p == ' ') p++;
    iso_to_ascii(p);

    if (strlen(p) > 118) {
        pprintf(mprn, _("   line: %.115s...\n").ptr, p);
    } else {
        pprintf(mprn, _("   line: %s\n").ptr, p);
    }

    numcount = 0;
    j = 1; /* for the constant */

    for (k=0; k<c.ncols && !err; k++) {
        int quoted = (*p == c.qchar);

        i = 0;
        while (*p && (*p != c.delim || quoted)) {
            if (i < CSVSTRLEN - 1) {
                c.str[i++] = *p;
            }
            if (i > 1 && *p == c.qchar) {
                quoted = 0;
            }
            p++;
        }
        c.str[i] = '\0';
        if (*p == c.delim) p++;

        if (k == 0 && csv_skip_col_1(c)) {
            ; /* no-op */
        } else if (!joining(c) && cols_subset(c) && skip_data_column(c, k)) {
            ; /* no-op */
        } else {
            if (joining(c)) {
                handle_join_varname(c, k, &j);
            } else if (probing(c) && csv_no_header(c)) {
                sprintf(c.dset.varname[j], "col%d", j);
                j++;
            } else {
                err = process_csv_varname(c, j, &numcount, prn);
                j++;
            }
        }
        if (j == c.dset.v) {
            break;
        }
    }

    if (!err && joining(c) && c.cols_list == null) {
        /* no relevant columns were found */
        gretl_errmsg_set(_("No relevant columns were found").ptr);
        err = GretlError.unkvar;
    }

    if (err) {
        return err;
    }

    if (csv_no_header(c) || numcount == c.dset.v - 1 ||
        obs_labels_no_varnames(obscol, c.dset, numcount)) {
        if (!csv_no_header(c)) {
            pputs(prn, _("it seems there are no variable names\n").ptr);
            /* then we undercounted the observations by one? */
            if (!rows_subset(c)) {
                err = add_single_obs(c.dset);
            }
        }
        if (!err) {
            /* set up to handle the "no varnames" case */
            csv_set_autoname(c);
            c.datapos = csv_has_bom(c) ? 3 : 0;
            if (!csv_all_cols(c)) {
                if (obs_labels_no_varnames(obscol, c.dset, numcount)) {
                    err = csv_reconfigure_for_markers(c.dset);
                    if (!err) {
                        csv_set_obs_column(c);
                    }
                }
            }
        }
    } else if (numcount > 0) {
        for (i=1; i<c.dset.v; i++) {
            if (check_varname(c.dset.varname[i])) {
                errmsg(1, prn);
                break;
            }
        }
        fprintf(stderr, "numcount = %d\n", numcount);
        err = GretlError.data;
    }

    return err;
}

private int row_not_wanted(CsvData * c, int t)
{
    if (c.rowmask != null) {
        if (t >= c.masklen) {
            return 1;
        } else if (gretl_vector_get(c.rowmask, t) == 0) {
            return 1;
        }
    }

    return 0;
}

/* read numerical data when we've been given a fixed column-reading
   specification */

private int fixed_format_read(CsvData * c, gzFile fp, Prn * prn)
{
    char * p;
    int miss_shown = 0;
    int * missp = null;
    int t = 0, s = 0;
    int i, k, n, m;
    int err = 0;

    c.real_n = c.dset.n;

    if (csv_has_bom(c)) {
        gzseek(fp, 3, SEEK_SET);
    }

    if (csv_is_verbose(c)) {
        missp = &miss_shown;
    }

    while (csv_fgets(c, fp) && !err) {
        tailstrip(c.line);
        if (*c.line == '#' || string_is_blank(c.line)) {
            continue;
        }
        if (row_not_wanted(c, s)) {
            s++;
            continue;
        }
        m = cast(int) strlen(c.line);
        for (i=1; i<=c.ncols && !err; i++) {
            k = c.cols_list[i];
            n = c.width_list[i];
            if (k + n - 1 > m) {
                /* attempting to read out of bounds */
                fprintf(stderr, "row %d, column %d: start=%d, width=%d, "
                        ~ "but line length = %d\n", t+1, i, k, n, m);
                err = GretlError.data;
                break;
            }
            p = c.line + k - 1;
            c.str[0] = '\0';
            strncat(c.str.ptr, p, n);
            /* Added 2016-11-16: allow trailing blanks in a field
               of specified width. This is required for handling
               US CPS data.
            */
            tailstrip(c.str.ptr);
            if (csv_missval(c.str.ptr, i, t+1, missp, prn)) {
                c.dset.Z[i][t] = NADBL;
            } else {
                c.dset.Z[i][t] = csv_atof(c, i);
                if (c.dset.Z[i][t] == NON_NUMERIC) {
                    gretl_errmsg_sprintf(_("At row %d, column %d:\n").ptr, t+1, k);
                    gretl_errmsg_sprintf(_("'%s' -- no numeric conversion performed!").ptr,
                                         c.str.ptr);
                    err = GretlError.data;
                }
            }
        }
        s++;
        if (++t == c.dset.n) {
            break;
        }
    }

    if (err == GretlError.data) {
        gretl_errmsg_set(_("Invalid column specification").ptr);
    }

    return err;
}

private bool XML1_OK(uint u)
{
    return (u>=0x0020 && u<=0xD7FF) || (u>=0xE000 && u<=0xFFFD);
}

/* Check that an observation label contains only
   valid UTF-8, and moreover that every character
   is valid in XML 1.0. If not, try recoding from
   ISO 8859.
*/

private int maybe_fix_csv_string(gchar * s)
{
    int err = 0;

    if (!g_utf8_validate(s, -1, null)) {
        GError * gerr = null;
        gsize wrote = 0;
        gchar * tr;

        /* try for iso-8859? */
        tr = g_convert(s, -1, "UTF-8", "ISO-8859-15",
                       null, &wrote, &gerr);
        if (gerr != null) {
            gretl_errmsg_set(gerr.message);
            g_error_free(gerr);
            err = GretlError.data;
        } else {
            *s = '\0';
            gretl_utf8_strncat(s, tr, CSVSTRLEN-1);
            g_free(tr);
        }
    }

    if (!err) {
        int i, n = cast(int) g_utf8_strlen(s, -1);
        uint u;

        for (i=0; i<n; i++) {
            u = g_utf8_get_char(s);
            if (!XML1_OK(u)) {
                return 0;
            }
            s = g_utf8_next_char(s);
        }
    }

    return err;
}

private void transcribe_obs_label(CsvData * c, int t)
{
    char * s = c.str.ptr;
    char c0 = *s;
    int n = cast(int) strlen(s);

    /* skip a leading quote, and unquote fully
       if a matching trailing quote is found
    */

    if (c0 == '"' || c0 == '\'') {
        if (s[n-1] == c0) {
            s[n-1] = '\0';
            n--;
        }
        s++;
        n--;
        /* and once more, with feeling... */
        if (s[0] == '\'') {
            s++;
            n--;
        }
    }

    if (n > OBSLEN - 1) {
        n = OBSLEN - 1;
    }

    c.dset.S[t][0] = '\0';
    gretl_utf8_strncat(c.dset.S[t], s, n);
}

private int real_read_labels_and_data(CsvData * c, gzFile fp, Prn * prn)
{
    char * p;
    int miss_shown = 0;
    int * missp = null;
    int truncated = 0;
    int t = 0, s = 0;
    int i, j, k;
    int err = 0;

    if (csv_is_verbose(c)) {
        missp = &miss_shown;
    }

    c.real_n = c.dset.n;

    while (csv_fgets(c, fp) && !err) {
        int inquote = 0;

        if (*c.line == '#' || string_is_blank(c.line)) {
            continue;
        } else if (c.skipstr[0] != '\0' && strstr(c.line, c.skipstr.ptr)) {
            c.real_n -= 1;
            continue;
        } else if (row_not_wanted(c, s)) {
            s++;
            continue;
        }

        compress_csv_line(c, 0);
        p = c.line;

        if (c.delim == ' ') {
            if (*p == ' ') p++;
        } else {
            p += strspn(p, " ");
        }

        j = 1;
        for (k=0; k<c.ncols && !err; k++) {
            i = 0;
            while (*p) {
                if (csv_keep_quotes(c) && *p == c.qchar) {
                    inquote = !inquote;
                } else if (!inquote && *p == c.delim) {
                    break;
                }
                if (i < CSVSTRLEN - 1) {
                    c.str[i++] = *p;
                } else {
                    truncated++;
                }
                p++;
            }
            c.str[i] = '\0';
            err = maybe_fix_csv_string(c.str.ptr);
            if (!err) {
                if (k == 0 && csv_skip_col_1(c) && c.dset.S != null) {
                    transcribe_obs_label(c, t);
                } else if (cols_subset(c) && skip_data_column(c, k)) {
                    ; /* no-op */
                } else {
                    err = process_csv_obs(c, j++, t, missp, prn);
                }
            }
            if (!err) {
                /* prep for next column */
                if (*p == c.delim) {
                    p++;
                }
                if (c.delim != ' ') {
                    p += strspn(p, " ");
                }
            }
        }

        s++;
        if (++t == c.dset.n) {
            break;
        }
    }

    if (truncated) {
        pprintf(prn, _("warning: %d strings were truncated.\n").ptr, truncated);
    }

    if (!err && c.real_n < c.dset.n) {
        int drop = c.dset.n - c.real_n;

        err = dataset_drop_observations(c.dset, drop);
    }

    return err;
}

/* When reading a CSV file, should we attempt to parse observation
   strings as dates (and impose time-series structure on the data
   if this is successful)? In general, yes, but maybe not if we're
   reading the data in the context of a "join" operation, since
   in this case automatic detection may collide with time-key
   information supplied by the user. Current status: we'll skip
   the auto-dating stuff when joining unless (a) it's a MIDAS
   join (mixed frequencies) and the user has _not_ supplied any
   time key specification.
*/

private int csv_skip_dates(CsvData * c)
{
    if (c.jspec != null) {
        /* with --aggr=spread (MIDAS) we'll need dates info,
           unless the user have a time key spec
        */
        return c.jspec.auto_midas == 0;
    } else {
        return 0;
    }
}

private int csv_read_data(CsvData * c, gzFile fp, Prn * prn, Prn * mprn)
{
    int reversed = csv_data_reversed(c);
    int err;

    if (mprn != null) {
        if (csv_all_cols(c)) {
            pputs(mprn, _("scanning for data...\n").ptr);
        } else {
            pputs(mprn, _("scanning for row labels and data...\n").ptr);
        }
    }

    gzseek(fp, c.datapos, SEEK_SET);

    err = real_read_labels_and_data(c, fp, prn);

    if (!err && csv_skip_col_1(c) && !rows_subset(c) && !csv_skip_dates(c)) {
        c.markerpd = test_markers_for_dates(c.dset, &reversed,
                                            c.skipstr.ptr, prn);
        if (reversed) {
            csv_set_data_reversed(c);
        }
    }

    return err;
}

private void print_csv_parsing_header(const(char) * fname, Prn * prn)
{
    if (!g_utf8_validate(fname, -1, null)) {
        gchar * trfname = g_locale_to_utf8(fname, -1, null, null, null);

        pprintf(prn, "%s %s...\n", _("parsing").ptr, trfname);
        g_free(trfname);
    } else {
        pprintf(prn, "%s %s...\n", _("parsing").ptr, fname);
    }
}

private int join_unique_columns(CsvData * c)
{
    const(char) ** cnames = c.jspec.colnames;
    char * counted;
    int i, j, ncols = 0;

    counted = cast(char *) calloc(c.jspec.ncols, 1);

    for (i=0; i<c.jspec.ncols; i++) {
        if (cnames[i] != null && counted[i] == 0) {
            counted[i] = 1;
            /* mark any duplicates as counted too */
            for (j=i+1; j<c.jspec.ncols; j++) {
                if (cnames[j] != null && !strcmp(cnames[i], cnames[j])) {
                    counted[j] = 1;
                }
            }
            ncols++;
        }
    }

    free(counted);

    return ncols;
}

private int csv_set_dataset_dimensions(CsvData * c)
{
    int err = 0;

    c.dset.v = 0;

    if (rows_subset(c)) {
        c.dset.n = n_from_row_mask(c);
    }

    if (fixed_format(c)) {
        if (c.dset.n == 0) {
            c.dset.n = c.nrows;
        }
        c.dset.v = c.ncols + 1;
    } else {
        int cols_wanted, cols_present;

        if (c.dset.n == 0) {
            if (csv_no_header(c)) {
                c.dset.n = c.nrows;
            } else if (c.nrows == 1) {
                csv_set_no_header(c);
                c.dset.n = 1;
            } else {
                /* allow for varnames row */
                c.dset.n = c.nrows - 1;
            }
        }

        cols_present = csv_skip_col_1(c) ? (c.ncols - 1) : c.ncols;

        if (joining(c)) {
            cols_wanted = join_unique_columns(c);
        } else if (cols_subset(c)) {
            cols_wanted = c.cols_list[0];
        } else {
            cols_wanted = cols_present;
        }

        if (cols_wanted > cols_present) {
            gretl_errmsg_set(_("Invalid column specification").ptr);
            err = GretlError.data;
        } else {
            /* allow for the constant */
            c.dset.v = cols_wanted + 1;
        }
    }

    if (probing(c)) {
        /* don't allocate tons of space for data that
           we won't read right now */
        c.dset.n = 1;
    }

    return err;
}

/*
 * real_import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @cols: column specification.
 * @rows: row specification.
 * @join: specification pertaining to "join" command.
 * @probe: also pertains to "join" (via GUI).
 * @pm: location of matrix to accept the data or NULL.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; use OPT_H
 * to indicate absence of a header row; use OPT_A to indicate that
 * all columns should be read as data series (i.e. do not try to
 * interpret the first column as observation labels); for use of
 * OPT_T see the help text for the "append" command.
 * @prn: gretl printing struct (or NULL).
 *
 * Open a Comma-Separated Values data file and read the data into
 * the current work space. Shared with gretl_join.c.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int real_import_csv(const(char) * fname,
                    DatasetRaw * dset,
                    const(char) * cols,
                    const(char) * rows,
                    JoinSpec * join,
                    void * probe,
                    GretlMatrixRaw ** pm,
                    Gretlopt opt,
                    Prn * prn)
{
    CsvData * c = null;
    gzFile fp = null;
    Prn * mprn = null;
    gchar * altname = null;
    int recode = 0;
    int popit = 0;
    int i, err = 0;

    import_na_init();

    if (gretl_messages_on()) {
        mprn = prn;
    }

    fp = gretl_gzopen(fname, "rb");
    if (fp == null) {
        pprintf(prn, _("Couldn't open %s\n").ptr, fname);
        err = GretlError.fopen;
        goto csv_bailout;
    }

    c = csvdata_new(dset);
    if (c == null) {
        err = GretlError.alloc;
        goto csv_bailout;
    }

    recode = csv_unicode_check(fp, c, prn);
    if (recode) {
        err = csv_recode_input(&fp, fname, &altname, recode, prn);
        if (err) {
            goto csv_bailout;
        }
    }

    if (cols != null) {
        err = csvdata_add_cols_list(c, cols, opt);
        if (err) {
            goto csv_bailout;
        } else if (fixed_format(c)) {
            pprintf(mprn, _("using fixed column format\n").ptr);
        }
    }

    if (rows != null) {
        err = csvdata_add_row_mask(c, rows);
        if (err) {
            goto csv_bailout;
        }
    }

    if (opt & Gretlopt.OPT_H) {
        csv_set_no_header(c);
    }

    if (join != null) {
        c.jspec = join;
        c.flags |= CsvFlags.havedata;
    } else if (probe != null) {
        c.probe = cast(CsvProbe *) probe;
        c.flags |= CsvFlags.havedata;
    } else {
        if (pm != null) {
            csv_set_as_matrix(c);
        }
        if (opt & Gretlopt.OPT_A) {
            csv_set_all_cols(c);
        }
        if (opt & Gretlopt.OPT_V) {
            csv_set_verbose(c);
        }
    }

    if (opt & Gretlopt.OPT_I) {
        csv_unset_keep_quotes(c);
    }

    if (mprn != null) {
        print_csv_parsing_header(fname, mprn);
    }

    /* get line length, also check for binary data, etc. */
    c.maxlinelen = csv_max_line_length(fp, c, prn);
    if (c.maxlinelen <= 0) {
        err = GretlError.data;
        goto csv_bailout;
    }

    if (csv_as_matrix(c) && csv_got_semi(c)) {
        if (c.delim == ',' && csv_got_delim(c)) {
            c.decpoint = ',';
        }
        c.delim = ';';
    } else if (!fixed_format(c) && !csv_got_delim(c)) {
        /* set default delimiter */
        if (csv_got_tab(c)) {
            c.delim = '\t';
        } else if (csv_got_semi(c)) {
            c.delim = ';';
        } else if (csv_got_pipe(c)) {
            c.delim = '|';
        } else {
            c.delim = ' ';
        }
    }

    /* buffer to hold lines */
    c.line = cast(char *) malloc(c.maxlinelen);
    if (c.line == null) {
        err = GretlError.alloc;
        goto csv_bailout;
    }

 alt_delim:

    if (mprn != null) {
        if (!fixed_format(c)) {
            pprintf(mprn, _("using delimiter '%c'\n").ptr, c.delim);
        }
        pprintf(mprn, _("   longest line: %d characters\n").ptr, c.maxlinelen - 1);
    }

    if (csv_has_trailing_comma(c) && c.delim != ',') {
        csv_unset_trailing_comma(c);
    }

    gzrewind(fp);

    /* read lines, check for consistency in number of fields */
    err = csv_fields_check(fp, c, mprn);
    if (err && !fixed_format(c)) {
        if (c.delim != ';' && csv_got_semi(c)) {
            c.delim = ';';
            err = 0;
            goto alt_delim;
        } else if (c.delim != '|' && csv_got_pipe(c)) {
            /* C: the delimiter is not changed here (see bugs.md #57) */
            err = 0;
            goto alt_delim;
        }
        pputs(prn, _(csv_msg).ptr);
        goto csv_bailout;
    }

    err = csv_set_dataset_dimensions(c);
    if (err) {
        err = GretlError.data;
        goto csv_bailout;
    }

    pprintf(mprn, _("   number of variables: %d\n").ptr, c.dset.v - 1);
    pprintf(mprn, _("   number of non-blank lines: %d\n").ptr, c.nrows);

    if (c.dset.v == 1 || c.dset.n == 0) {
        pputs(prn, _("Invalid data file\n").ptr);
        err = GretlError.data;
        goto csv_bailout;
    }

    /* initialize CSV dataset */
    err = start_new_Z(c.dset, cast(Gretlopt) 0);
    if (!err && csv_skip_col_1(c)) {
        err = dataset_allocate_obs_markers(c.dset);
    }

    if (err) {
        goto csv_bailout;
    }

    /* second pass */

    gzrewind(fp);

    if (fixed_format(c)) {
        err = fixed_format_read(c, fp, prn);
        if (err) {
            goto csv_bailout;
        } else {
            csv_set_autoname(c);
            goto csv_continue;
        }
    }

    err = csv_varname_scan(c, fp, prn, mprn);
    if (err || probing(c)) {
        goto csv_bailout;
    }

    if (c.decpoint == '.' && get_local_decpoint() == ',') {
        /* we're in a locale that uses decimal comma:
           switch to the C locale */
        gretl_push_c_numeric_locale();
        popit = 1;
    } else if (c.decpoint == ',' && get_local_decpoint() == '.') {
        /* dotsub: define this if we're in a '.' locale and
           we've figured that the decimal character is ',' in
           the file we're reading
        */
        csv_set_dotsub(c);
    }

    err = csv_read_data(c, fp, prn, mprn);

    if (!err) {
        /* try again, under certain conditions */
        if (csv_skip_bad(c)) {
            err = csv_read_data(c, fp, prn, null);
        } else if (c.thousep > 0) {
            pprintf(mprn, (_("WARNING: it seems '%c' is being used ")
                           ~ _("as thousands separator\n")).ptr, c.thousep);
            c.decpoint = (c.thousep == '.')? ',' : '.';
            if (c.decpoint == ',') {
                if (get_local_decpoint() == '.') {
                    csv_set_dotsub(c);
                } else if (popit) {
                    gretl_pop_c_numeric_locale();
                    popit = 0;
                }
            }
            revise_non_numeric_values(c);
            csv_set_scrub_thousep(c);
            err = csv_read_data(c, fp, prn, null);
        }
    }

    if (!err && !probing(c)) {
        err = csv_non_numeric_check(c, prn);
        if (!err && csv_has_non_numeric(c)) {
            /* try once more */
            err = csv_read_data(c, fp, prn, null);
        }
    }

    if (popit) {
        gretl_pop_c_numeric_locale();
    }

    if (err) {
        goto csv_bailout;
    }

    if (csv_data_reversed(c)) {
        reverse_data(c.dset, mprn);
    }

 csv_continue:

    c.dset.t1 = 0;
    c.dset.t2 = c.dset.n - 1;

    if (c.markerpd > 0) {
        pputs(mprn, _("taking date information from row labels\n\n").ptr);
        if (csv_skip_bad(c)) {
            pprintf(prn, (_("WARNING: Check your data! gretl has stripped out ")
                    ~ _("what appear to be\nextraneous lines in a %s dataset: ")
                    ~ _("this may not be right.\n\n")).ptr,
                    (c.dset.pd == 4)? "quarterly".ptr : "monthly".ptr);
        }
    } else {
        pputs(mprn, _("treating these as undated data\n\n").ptr);
        dataset_obs_info_default(c.dset);
    }

    if (c.dset.pd != 1 || strcmp(c.dset.stobs.ptr, "1")) {
        c.dset.structure = TsCodes.time_series;
    }

    if (c.st != null) {
        err = gretl_string_table_validate(c.st, Gretlopt.OPT_NONE);
        if (err) {
            pputs(prn, _("Failed to interpret the data as numeric\n").ptr);
            goto csv_bailout;
        } else if (joining(c)) {
            gretl_string_table_save(c.st, c.dset);
        } else {
            gretl_string_table_finalize(c.st, c.dset);
        }
    }

    if (csv_as_matrix(c)) {
        /* FIXME placement of this */
        if (csv_autoname(c)) {
            strings_array_free(c.dset.varname, c.dset.v);
            c.dset.varname = null;
        }
        *pm = gretl_matrix_data_subset(null, c.dset, -1, -1,
                                       MMissingCode.ok, &err);
        goto csv_bailout;
    }

    /* If there were observation labels and they were not interpretable
       as dates, and they weren't simply "1, 2, 3, ...", then they
       should probably be preserved; otherwise discard them.
    */
    if (c.dset.S != null && c.markerpd >= 0 &&
        c.dset.markers != DatasetMarkerType.daily_date_strings) {
        dataset_destroy_obs_markers(c.dset);
    }

    if (csv_autoname(c)) {
        /* no variable names were found */
        for (i=1; i<c.dset.v; i++) {
            sprintf(c.dset.varname[i], "v%d", i);
        }
    } else {
        if (fix_varname_duplicates(c.dset)) {
            pputs(prn, _("warning: some variable names were duplicated\n").ptr);
        }
    }

    if (!joining(c) && !probing(c)) {
        int newdata = (dset.Z == null);

        /* not doing a special "join" operation */
        err = merge_or_replace_data(dset, &c.dset, get_merge_opts(opt), prn);

        if (!err && newdata && c.descrip != null) {
            dset.descrip = c.descrip;
            c.descrip = null;
        }
        if (!err && newdata) {
            dataset_add_import_info(dset, fname, GretlFileType.csv);
        }
    }

 csv_bailout:

    if (fp != null) {
        gzclose(fp);
    }

    if (!err && c.jspec != null) {
        c.jspec.c = c;
    } else if (!err && c.probe != null) {
        c.probe.dset = c.dset;
        c.dset = null;
        csvdata_free(c);
    } else {
        csvdata_free(c);
    }

    if (altname != null) {
        gretl_remove(altname);
        g_free(altname);
    }

    if (err == GretlError.alloc) {
        pputs(prn, _("Out of memory\n").ptr);
    }

    return err;
}

/**
 * import_csv:
 * @fname: name of CSV file.
 * @dset: dataset struct.
 * @opt: use OPT_N to force interpretation of data colums containing
 * strings as coded (non-numeric) values and not errors; for use of
 * OPT_T see the help for "append".
 * @prn: gretl printing struct (or NULL).
 *
 * Open a Comma-Separated Values data file and read the data into
 * the current work space.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int import_csv(const(char) * fname, DatasetRaw * dset, Gretlopt opt, Prn * prn) {
    const(char) * cols = null;
    const(char) * rows = null;
    GretlCmdIndex ci;
    int err;

    err = incompatible_options(opt, Gretlopt.OPT_F | Gretlopt.OPT_L);
    if (err) {
        /* --cols and --fixed-cols */
        return err;
    }

    ci = (dset != null && dset.v > 0)? GretlCmdIndex.append : GretlCmdIndex.open;

    if (opt & Gretlopt.OPT_F) {
        /* we should have a "--fixed-cols=XXX" specification */
        cols = get_optval_string(ci, Gretlopt.OPT_F);
        if (cols == null || *cols == '\0') {
            return GretlError.parse;
        }
    } else if (opt & Gretlopt.OPT_L) {
        /* should have a "--cols=XXX" specification */
        cols = get_optval_string(ci, Gretlopt.OPT_L);
        if (cols == null || *cols == '\0') {
            return GretlError.parse;
        }
    }

    if (opt & Gretlopt.OPT_M) {
        /* we should have a "--rowmask=XXX" specification */
        rows = get_optval_string(ci, Gretlopt.OPT_M);
        if (rows == null || *rows == '\0') {
            return GretlError.parse;
        }
    }

    return real_import_csv(fname, dset, cols, rows, null, null, null, opt, prn);
}

GretlMatrixRaw * import_csv_as_matrix(const(char) * fname, int * err) {
    Prn * prn = null;
    GretlMatrixRaw * m = null;
    char[MAXLEN] csvname = 0;
    Gretlopt opt = Gretlopt.OPT_A; /* --all-cols */
    int http = 0;

    *err = try_http(fname, csvname.ptr, &http);

    if (!*err && http) {
        *err = real_import_csv(csvname.ptr, null, null, null,
                               null, null, &m, opt, prn);
    } else if (!*err) {
        char[FILENAME_MAX] fullname = 0;

        if (strlen(fname) >= fullname.sizeof) {
            gretl_errmsg_set(_("Filename is too long").ptr);
            *err = GretlError.data;
        } else {
            strcpy(fullname.ptr, fname);
            gretl_maybe_prepend_dir(fullname.ptr);
            *err = real_import_csv(fullname.ptr, null, null, null,
                                   null, null, &m, opt, prn);
        }
    }

    gretl_print_destroy(prn);

    return m;
}

private int probe_varnames_check(DatasetRaw * dset, Gretlopt opt,
                                 int * rerun)
{
    int missnames = 0;
    int i, err = 0;

    for (i=1; i<dset.v; i++) {
        if (dset.varname[i][0] == '\0') {
            missnames = 1;
            break;
        }
    }

    if (missnames) {
        if (opt & Gretlopt.OPT_H) {
            gretl_errmsg_set(_("Couldn't find all variable names").ptr);
            err = GretlError.data;
        } else {
            *rerun = 1;
        }
    }

    return err;
}

/**
 * probe_csv:
 * @fname: name of CSV file.
 * @varnames: location to receive variable names.
 * @nvars: location to receive number of variables (columns).
 * @opt: on input, may contain any extra options to pass to
 * real_import_csv(); on return, OPT_H (indicating that the
 * CSV file has no header) may be added if it seems to be
 * required (no header).
 *
 * Open a Comma-Separated Values data file and read enough to
 * determine the variable names.
 *
 * Returns: 0 on successful completion, non-zero otherwise.
 */

int probe_csv(const(char) * fname, char *** varnames,
              int * nvars, Gretlopt * opt)
{
    CsvProbe probe;
    int err;

    probe.dset = null;

    err = real_import_csv(fname, null, null, null, null,
                          &probe, null, *opt, null);

    if (!err) {
        int rerun = 0;

        err = probe_varnames_check(probe.dset, *opt, &rerun);

        if (err || rerun) {
            destroy_dataset(probe.dset);
            probe.dset = null;
        }

        if (!err && rerun) {
            /* try again with --no-header flag */
            *opt |= Gretlopt.OPT_H;
            err = real_import_csv(fname, null, null, null, null,
                                  &probe, null, *opt, null);
        }

        if (!err) {
            /* steal the varname array */
            *varnames = probe.dset.varname;
            *nvars = probe.dset.v;
            probe.dset.varname = null;
        }

        destroy_dataset(probe.dset);
    }

    return err;
}

int csv_open_needs_matrix(Gretlopt opt)
{
    int ret = 0;

    if (opt & Gretlopt.OPT_M) {
        /* --rowmask=matrix */
        ret = 1;
    } else if (opt & Gretlopt.OPT_F) {
        /* --fixed-cols=whatever */
        const(char) * s = get_optval_string(GretlCmdIndex.open, Gretlopt.OPT_F);

        ret = get_matrix_by_name(s) != null;
    }

    return ret;
}

/* Read the first @n_lines lines of @fname and write the
   content to @prn.
*/

int peek_at_csv(const(char) * fname, int n_lines, Prn * prn)
{
    CsvData * c = null;
    gzFile fp = null;
    gchar * altname = null;
    int recode;
    int err = 0;

    fp = gretl_gzopen(fname, "rb");
    if (fp == null) {
        return GretlError.fopen;
    }

    c = csvdata_new(null);
    if (c == null) {
        err = GretlError.alloc;
    }

    if (!err) {
        recode = csv_unicode_check(fp, c, prn);
        if (recode) {
            err = csv_recode_input(&fp, fname, &altname, recode, prn);
        }
    }

    if (!err) {
        c.maxlinelen = csv_max_line_length(fp, c, prn);
        if (c.maxlinelen <= 0) {
            err = GretlError.data;
        }
    }

    if (!err) {
        c.line = cast(char *) malloc(c.maxlinelen);
        if (c.line == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        const(char) * s;
        int i;

        gzrewind(fp);
        gretl_print_reset_buffer(prn);

        for (i=0; i<n_lines; i++) {
            s = gzgets(fp, c.line, c.maxlinelen);
            if (s == null) {
                break;
            } else {
                pputs(prn, s);
            }
        }
    }

    gzclose(fp);
    csvdata_free(c);

    if (altname != null) {
        gretl_remove(altname);
        g_free(altname);
    }

    return err;
}
