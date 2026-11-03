/* Port of plugin/stata_export.c: writer for Stata .dta files (format 113, Stata 8).

   In the C this is a gretl plugin (loaded by get_plugin_function("stata_export")); here it is an ordinary
   module that plugins.d returns for that name. stata_export() keeps its C signature.

   As in the C, the file is written in the host byte order (the header's byteorder byte says which). The
   C's use of u8_to_ascii_convert() (gretl_intl.c) and gretl_minmax() (describe.c), both not ported yet, is
   through per-file extern (C) declarations. The C has "use_strtypes = 1" ("experiment"), which makes the
   value-label branches (dta_value_labels_write, labeled > 0) dead code; they are ported anyway. */

module gretl.plugin.stata_export;


import core.stdc.ctype : tolower;
import core.stdc.locale : LC_TIME, setlocale;
import core.stdc.math : floor;
import core.stdc.stdio : FILE, fclose, fprintf, fwrite, sprintf, sscanf, stderr;
import core.stdc.stdlib : atoi, calloc, free, malloc;
import core.stdc.string : memset, strcat, strcpy, strlen, strspn;
import core.stdc.time : localtime, strftime, time, time_t, tm;

import gretl.calendar : epoch_day_from_t, get_epoch_day;
import gretl.dataio : ntolabel;
import gretl.dataseth : calendar_data, dataset_has_markers, dataset_is_time_series;
import gretl.dataset : is_string_valued, series_get_label, series_get_string_for_obs, series_get_string_table,
    series_get_string_vals;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : in_gretl_list;
import gretl.gretl_paths : gretl_fopen;
import gretl.gretl_string_table : SeriesTable, series_table_get_strings;
import gretl.gretl_utils : gretl_isdummy;
import gretl.libgretlh;
import gretl.missingh : na;
import gretl.versionh : GRETL_VERSION;
import gretl.gretl_intl : u8_to_ascii_convert;
import gretl.describe : gretl_minmax;
import gretl.genmain : current_series_index;

/* unported dependencies: per-file declarations from the 2026c headers */
private enum StataType {
    byte_ = 251,
    int_,
    long_,
    float_,
    double_
}

private enum TimeKind {
    auto_ = 1,
    cal
}

/* default missing-value codes (".") */
private enum STATA_DOUBLE_NA = 0x1.0000000000000p+1023;
private enum STATA_LONG_NA   = 2147483621;
private enum STATA_INT_NA    = 32741;
private enum STATA_BYTE_NA   = 101;

/* subtract this from epoch day to get Stata date */
private enum STATA_DAY_OFFSET = 715523;

private bool dta_113_string_type(int t) { return t > 0 && t < 245; }

private struct Dta113Hdr {
    ubyte ds_format;         /* contains 113 = 0x71 */
    ubyte byteorder;         /* 0x01 -> HILO, 0x02 -> LOHI */
    ubyte filetype;          /* 0x01 */
    ubyte unused;            /* 0x01 */
    short nvar;              /* number of vars */
    int nobs;                /* number of obs */
    char[81] data_label;     /* dataset label */
    char[18] time_stamp;     /* date/time saved */
}

private struct Dta113ValueLabelTable {
    int n;          /* number of entries */
    int txtlen;     /* length of txt[] */
    int off;        /* txt[] offset table (array) */
    int val;        /* sorted value table (array) */
    char * txt;     /* text table */
}

private struct Dta113ValueLabel {
    int len;                    /* length of vltable */
    char[33] labname;           /* name of table */
    char[3] padding;
    Dta113ValueLabelTable table;    /* the table itself */
}

private void dta_value_labels_write(FILE * fp, const(DatasetRaw) * dset,
                                    int v, long * w)
{
    Dta113ValueLabel lbl;
    Dta113ValueLabelTable tab;
    char ** S;
    char * tmp;
    int len, maxlen = 0;
    int i, ns;

    S = series_get_string_vals(dset, v, &ns, 0);

    for (i=0; i<ns; i++) {
        len = cast(int) strlen(S[i]);
        if (len > maxlen) {
            maxlen = len;
        }
    }

    tmp = cast(char *) malloc(maxlen + 1);
    if (tmp == null) {
        fprintf(stderr, "dta_value_labels_write: out of memory\n");
        return;
    }

    len = 0;
    for (i=0; i<ns; i++) {
        u8_to_ascii_convert(tmp, S[i], maxlen, '?');
        len += cast(int) strlen(tmp) + 1;
    }

    lbl.len = len + (ns + 1) * 8;
    memset(lbl.labname.ptr, 0, 33);
    sprintf(lbl.labname.ptr, "S%d", v);
    memset(lbl.padding.ptr, 0, 3);

    *w += fwrite(&lbl.len, 4, 1, fp);
    *w += fwrite(lbl.labname.ptr, 1, 33, fp);
    *w += fwrite(lbl.padding.ptr, 1, 3, fp);

    tab = lbl.table;
    tab.n = ns;
    tab.txtlen = len;

    *w += fwrite(&tab.n, 4, 1, fp);
    *w += fwrite(&tab.txtlen, 4, 1, fp);

    tab.off = 0;
    for (i=0; i<tab.n; i++) {
        *w += fwrite(&tab.off, 4, 1, fp);
        tab.off += cast(int) strlen(S[i]) + 1;
    }

    tab.val = 1;
    for (i=0; i<tab.n; i++) {
        *w += fwrite(&tab.val, 4, 1, fp);
        tab.val++;
    }

    for (i=0; i<tab.n; i++) {
        u8_to_ascii_convert(tmp, S[i], maxlen, '?');
        *w += fwrite(tmp, 1, strlen(tmp) + 1, fp);
    }

    free(tmp);
}

private void dta_hdr_write(Dta113Hdr * hdr, FILE * fp, long * w)
{
    *w += fwrite(&hdr.ds_format, 1, 1, fp);
    *w += fwrite(&hdr.byteorder, 1, 1, fp);
    *w += fwrite(&hdr.filetype, 1, 1, fp);
    *w += fwrite(&hdr.unused, 1, 1, fp);
    *w += fwrite(&hdr.nvar, 2, 1, fp);
    *w += fwrite(&hdr.nobs, 4, 1, fp);
    *w += fwrite(hdr.data_label.ptr, 81, 1, fp);
    *w += fwrite(hdr.time_stamp.ptr, 18, 1, fp);
}

private void asciify_to_length(char * targ, const(char) * src, int len)
{
    if (src != null && *src != '\0') {
        /* skip leading white space */
        src += strspn(src, " \t\n");
        /* and do our best to replace UTF-8 with ASCII */
        u8_to_ascii_convert(targ, src, len, '?');
    } else {
        *targ = '\0';
    }
}

private int write_stata_strval(const(DatasetRaw) * dset, int v, int t,
                               ubyte len, FILE * fp)
{
    /* @len can be as large as 244 (see get_string_type()), and
       asciify_to_length()/u8_to_ascii_convert() always append a
       terminating NUL after writing up to @len bytes, so the
       buffer must hold @len + 1 bytes: size it for the full
       range of guint8 to be safe regardless of @len
    */
    char[256] buf = 0;
    const(char) * s;

    s = series_get_string_for_obs(dset, v, t);
    if (s != null) {
        asciify_to_length(buf.ptr, s, len);
    }

    return cast(int) fwrite(buf.ptr, 1, len, fp);
}

/* Note: see https://www.stata.com/help.cgi?dta_113 for detail
   on the dta_113 format. Also see
   https://www.loc.gov/preservation/digital/formats/fdd/fdd000471.shtml
   according to which the 113 format was introduced with Stata 8
   in January 2003.
*/

private void dta_hdr_init(Dta113Hdr * hdr, const(DatasetRaw) * dset,
                          int nvars)
{
    /* ENABLE_NLS is assumed defined, as in the stock build */
    char * localt = null;
    long now = time(null);
    tm * local;

    hdr.ds_format = 113;
    version (BigEndian) {
        hdr.byteorder = 0x01;
    } else {
        hdr.byteorder = 0x02;
    }
    hdr.filetype = 0x01;
    hdr.unused = 0x01;
    hdr.nvar = cast(short) nvars;
    hdr.nobs = dset.t2 - dset.t1 + 1;

    memset(hdr.data_label.ptr, 0, 81);
    sprintf(hdr.data_label.ptr, "Written by gretl %s", GRETL_VERSION.ptr);

    memset(hdr.time_stamp.ptr, 0, 18);
    localt = setlocale(LC_TIME, null);
    setlocale(LC_TIME, "C");
    local = localtime(cast(time_t *) &now);
    /* dd Mon yyyy hh:mm (in C locale) */
    strftime(hdr.time_stamp.ptr, 18, "%d %b %Y %I:%M", local);
    setlocale(LC_TIME, localt);
}

private int include_var(const(int) * list, int i)
{
    if (list == null) {
        return 1;
    } else {
        return in_gretl_list(list, i);
    }
}

private int x_isint(int t1, int t2, const(double) * x)
{
    int t, ret = 1;

    for (t=t1; t<=t2; t++) {
        if (!na(x[t]) && x[t] != floor(x[t])) {
            ret = 0;
            break;
        }
    }

    return ret;
}

private int get_string_type(const(DatasetRaw) * dset, int v)
{
    const(SeriesTable) * st = series_get_string_table(dset, v);
    int ret = 0;

    if (st != null) {
        int i, ni, ns;
        char ** S = series_table_get_strings(st, &ns);

        for (i=0; i<ns; i++) {
            ni = cast(int) strlen(S[i]);
            if (ni > ret) {
                ret = ni;
            }
        }
        if (ret > 244) {
            ret = 244;
        } else if (ret < 244) {
            ret++;
        }
    }

    return ret;
}

private ubyte * make_types_array(const(DatasetRaw) * dset,
                                 const(int) * list,
                                 int use_strtypes,
                                 int * nvars)
{
    ubyte st;
    ubyte * t;

    if (list != null) {
        *nvars = list[0];
    } else {
        *nvars = dset.v - 1;
    }

    t = cast(ubyte *) calloc(*nvars, 1);

    if (t != null) {
        const(double) * x;
        double xmin, xmax;
        int i, n_ok, j = 0;

        for (i=1; i<dset.v; i++) {
            if (include_var(list, i)) {
                if (use_strtypes && is_string_valued(dset, i)) {
                    st = cast(ubyte) get_string_type(dset, i);
                    if (st > 0) {
                        t[j] = st;
                    }
                }
                if (t[j] == 0) {
                    x = dset.Z[i];
                    n_ok = gretl_minmax(dset.t1, dset.t2, x, &xmin, &xmax);
                    if (n_ok == 0) {
                        /* all missing, hmm */
                        t[j] = StataType.byte_;
                    } else if (gretl_isdummy(dset.t1, dset.t2, x)) {
                        t[j] = StataType.byte_;
                    } else if (x_isint(dset.t1, dset.t2, x)) {
                        if (xmin > -10000 && xmax < 10000) {
                            t[j] = StataType.int_;
                        } else if (xmin > -200000000 && xmax < 200000000) {
                            t[j] = StataType.long_;
                        } else {
                            t[j] = StataType.double_;
                        }
                    } else {
                        t[j] = StataType.double_;
                    }
                }
                j++;
            }
        }
    }

    return t;
}

/* For time series, get the starting observation in
   the numerical form wanted by Stata; also get a
   suitable name for the added variable.
*/

private int get_stata_t0(const(DatasetRaw) * dset,
                         char * timevar,
                         int * add_time)
{
    int t0 = 0;
    char[OBSLEN] obs = 0;

    *timevar = '\0';
    ntolabel(obs.ptr, dset.t1, dset);

    if (dset.pd == 1 && dset.sd0 > 999) {
        /* we just want the year */
        t0 = atoi(obs.ptr);
        strcpy(timevar, "year");
        *add_time = TimeKind.auto_;
    } else if ((dset.pd == 4 || dset.pd == 12) && dset.sd0 > 999) {
        /* quarters or months since the start of 1960 */
        int y, p;
        char c;

        sscanf(obs.ptr, "%d%c%d", &y, &c, &p);
        t0 = dset.pd * (y - 1960) + p - 1;
        strcpy(timevar, dset.pd == 4 ? "quarter" : "month");
        *add_time = TimeKind.auto_;
    } else if (calendar_data(dset)) {
        strcpy(timevar, "date");
        *add_time = TimeKind.cal;
    } else {
        t0 = atoi(obs.ptr);
        strcpy(timevar, "generic_t");
        *add_time = TimeKind.auto_;
    }

    if (*timevar != '\0') {
        /* don't collide with regular series */
        if (current_series_index(dset, timevar) > 0) {
            strcat(timevar, "_t");
        }
    }

    return t0;
}

private void make_timevar_label(char * buf, const(char) * tvar)
{
    if (*tvar == 'y') {
        strcpy(buf, "years since 0 CE");
    } else if (*tvar == 'q') {
        strcpy(buf, "quarters since 1960q1");
    } else if (*tvar == 'm') {
        strcpy(buf, "months since 1960m1");
    } else if (*tvar == 'g') {
        strcpy(buf, "period of observation");
    } else if (*tvar == 'd') {
        strcpy(buf, "date of observation");
    }
}

private int stata_date(int t, const(DatasetRaw) * dset)
{
    long ed = -1;

    if (dataset_has_markers(dset)) {
        ed = get_epoch_day(dset.S[t]);
    } else {
        ed = epoch_day_from_t(t, dset);
    }

    return cast(int) (ed - STATA_DAY_OFFSET);
}

private void write_lower_case(char * targ, const(char) * src)
{
    int i, n = cast(int) strlen(src);

    for (i=0; i<n; i++) {
        targ[i] = cast(char) tolower(src[i]);
    }
    targ[i] = '\0';
}

extern (C) int stata_export(const(char) * fname,
                            const(int) * list,
                            Gretlopt opt,
                            const(DatasetRaw) * dset)
{
    FILE * fp;
    Dta113Hdr hdr;
    char[96] buf = 0;
    char[16] timevar = 0;
    long w = 0;
    ubyte * types;
    int use_strtypes;
    double xit;
    int i32, t32 = 0;
    short i16;
    ubyte i8;
    int missing;
    int add_time = 0;
    int i, j, t;
    int labeled = 0;
    int nv = 0;
    int err = 0;

    /* experiment */
    use_strtypes = 1;

    fp = gretl_fopen(fname, "wb");
    if (fp == null) {
        return GretlError.fopen;
    }

    types = make_types_array(dset, list, use_strtypes, &nv);
    if (types == null) {
        fclose(fp);
        return GretlError.alloc;
    }

    timevar[0] = '\0';

    if (dataset_is_time_series(dset)) {
        t32 = get_stata_t0(dset, timevar.ptr, &add_time);
    }

    dta_hdr_init(&hdr, dset, nv + (add_time > 0));
    dta_hdr_write(&hdr, fp, &w);

    /*
        typlist        nvar    byte array
        varlist     33*nvar    char array
        srtlist  2*(nvar+1)    int array
        fmtlist     12*nvar    char array
        lbllist     33*nvar    char array
    */

    /* typlist */
    if (add_time) {
        i8 = StataType.long_;
        w += fwrite(&i8, i8.sizeof, 1, fp);
    }
    for (j=0; j<nv; j++) {
        w += fwrite(&types[j], 1, 1, fp);
    }

    /* varlist (names) */
    if (add_time) {
        memset(buf.ptr, 0, 33);
        strcat(buf.ptr, timevar.ptr);
        w += fwrite(buf.ptr, 1, 33, fp);
    }
    for (i=1; i<dset.v; i++) {
        if (include_var(list, i)) {
            memset(buf.ptr, 0, 33);
            if (opt & Gretlopt.OPT_L) {
                /* --lcnames */
                write_lower_case(buf.ptr, dset.varname[i]);
            } else {
                strcat(buf.ptr, dset.varname[i]);
            }
            w += fwrite(buf.ptr, 1, 33, fp);
        }
    }

    /* srtlist */
    i16 = 0;
    if (add_time) {
        w += fwrite(&i16, 2, 1, fp);
    }
    for (j=0; j<=nv; j++) {
        w += fwrite(&i16, 2, 1, fp);
    }

    /* fmtlist */
    if (add_time) {
        memset(buf.ptr, 0, 12);
        sprintf(buf.ptr, "%%t%c", cast(int) timevar[0]);
        w += fwrite(buf.ptr, 1, 12, fp);
    }
    for (j=0; j<nv; j++) {
        memset(buf.ptr, 0, 12);
        if (dta_113_string_type(types[j])) {
            /* the C has a commented-out sprintf(buf, "%%%ds", types[j]) here */
        } else if (types[j] == StataType.byte_ || types[j] == StataType.int_) {
            strcpy(buf.ptr, "%8.0g");
        } else if (types[j] == StataType.long_) {
            strcpy(buf.ptr, "%12.0g");
        } else {
            strcpy(buf.ptr, "%9.0g");
        }
        w += fwrite(buf.ptr, 1, 12, fp);
        fprintf(stderr, "wrote format %d, '%s'\n", j+1, buf.ptr);
    }

    /* lbllist */
    if (add_time) {
        memset(buf.ptr, 0, 33);
        w += fwrite(buf.ptr, 1, 33, fp);
    }
    for (i=1; i<dset.v; i++) {
        if (include_var(list, i)) {
            memset(buf.ptr, 0, 33);
            if (!use_strtypes && is_string_valued(dset, i)) {
                sprintf(buf.ptr, "S%d", i);
                labeled++;
            }
            w += fwrite(buf.ptr, 1, 33, fp);
        }
    }

    /* Variable labels */
    if (add_time) {
        memset(buf.ptr, 0, 81);
        make_timevar_label(buf.ptr, timevar.ptr);
        w += fwrite(buf.ptr, 1, 81, fp);
    }
    for (i=1; i<dset.v; i++) {
        if (include_var(list, i)) {
            const(char) * vlabel = series_get_label(dset, i);

            memset(buf.ptr, 0, 81);
            if (vlabel != null && *vlabel != '\0') {
                asciify_to_length(buf.ptr, vlabel, 80);
            }
            w += fwrite(buf.ptr, 1, 81, fp);
        }
    }

    /* "Expansion fields" (a Stata-private mystery) */
    i8 = 0;
    for (i=0; i<5; i++) {
        w += fwrite(&i8, 1, 1, fp);
    }

    /* Data */
    for (t=dset.t1; t<=dset.t2; t++) {
        if (add_time == TimeKind.auto_) {
            w += fwrite(&t32, 4, 1, fp);
            t32++;
        } else if (add_time == TimeKind.cal) {
            t32 = stata_date(t, dset);
            w += fwrite(&t32, 4, 1, fp);
        }
        j = 0;
        for (i=1; i<dset.v; i++) {
            if (include_var(list, i)) {
                xit = dset.Z[i][t];
                missing = na(xit);
                if (dta_113_string_type(types[j])) {
                    w += write_stata_strval(dset, i, t, types[j], fp);
                } else if (types[j] == StataType.byte_) {
                    i8 = missing ? STATA_BYTE_NA : cast(ubyte) xit;
                    w += fwrite(&i8, i8.sizeof, 1, fp);
                } else if (types[j] == StataType.int_) {
                    i16 = missing ? STATA_INT_NA : cast(short) xit;
                    w += fwrite(&i16, i16.sizeof, 1, fp);
                } else if (types[j] == StataType.long_) {
                    i32 = missing ? STATA_LONG_NA : cast(int) xit;
                    w += fwrite(&i32, i32.sizeof, 1, fp);
                } else {
                    xit = missing ? STATA_DOUBLE_NA : xit;
                    w += fwrite(&xit, xit.sizeof, 1, fp);
                }
                j++;
            }
        }
    }

    if (labeled > 0) {
        /* write string values as value labels */
        for (i=1; i<dset.v; i++) {
            if (include_var(list, i) && is_string_valued(dset, i)) {
                dta_value_labels_write(fp, dset, i, &w);
            }
        }
    }

    fclose(fp);
    free(types);

    fprintf(stderr, "stata_export: wrote %d elements\n", cast(int) w);

    return err;
}
