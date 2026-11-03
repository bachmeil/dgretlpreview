/* Port of plugin/stata_import.c: reader for Stata .dta files, versions 5.0 to 14 (formats 105 to 118).

   Originally based on stataread.c from the GNU R "foreign" package (Thomas Lumley, Saikat DebRoy) and
   modified for gretl by Allin Cottrell (2005-2015).

   In the C this is a gretl plugin (loaded by get_plugin_function("dta_get_data")); here it is an
   ordinary module that plugins.d returns for that name. dta_get_data() keeps its C signature (note
   that, unlike the spreadsheet importers, it takes no @list, @sheetname or parent arguments).

   Differences from the C:
   - the file-scope globals of the C (stata_version, stata_OLD, ...) are private __gshared variables,
     re-initialised on each call as in the C;
   - the anonymous enums of the C are StataSeType (STATA_SE_*) and Stata13Type (STATA_13_*);
   - the function-like macros (stata_type_float(), STATA_FLOAT_NA(), ...) are private functions;
   - only the 64-bit-offset POSIX ftello/fseeko are used (the Windows variants are dropped);
   - the STRDEBUG and HDR_DEBUG diagnostics (compiled out in the C) are dropped;
   - the GINTn_FROM_BE/LE macros are the small helper functions below. */

module gretl.plugin.stata_import;


import core.bitop : bswap;
import core.stdc.limits : INT_MAX;
import core.stdc.stdio : FILE, SEEK_CUR, SEEK_SET, fclose, feof, fprintf, fputc, fputs, fread, printf,
    rewind, sprintf, sscanf, stderr;
import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memcpy, strcat, strcmp, strcpy, strlen, strncat;
version (Posix) import core.sys.posix.stdio : fseeko, ftello; else import gretl.winport : fseeko, ftello;

import gretl.calendar : ymd_bits_from_epoch_day, ymd_extended_from_epoch_day;
import gretl.dataio : check_varname, dataset_add_import_info, get_date_x, get_merge_opts, merge_or_replace_data,
    ntolabel;
import gretl.dataioh : GretlFileType, free_datainfo;
import gretl.dataseth : DatasetMarkerType;
import gretl.dataset : datainfo_new, dataset_allocate_obs_markers, dataset_destroy_obs_markers,
    dataset_obs_info_default, destroy_dataset, series_set_discrete, series_set_label, start_new_Z;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_append_term, gretl_list_new;
import gretl.gretl_paths : gretl_fopen;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new, pprintf, pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_string_table : GretlStringTable, gretl_string_table_add_extra, gretl_string_table_destroy,
    gretl_string_table_finalize, gretl_string_table_index, gretl_string_table_new;
import gretl.glib : gchar, gsize;
import gretl.gmem : g_free;
import gretl.gstrfuncs : g_convert, g_get_charset, g_locale_to_utf8, g_strdup, g_utf8_validate;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.strutils : fix_varname_duplicates, gretl_namechar_spn, gretl_utf8_truncate_b, strings_array_add,
    strings_array_free;
import gretl.swap_bytesh : reverse_double, reverse_float, reverse_int;
import gretl.gretl_intl : asciify_utf8_varname, iso_to_ascii;

/* unported dependencies: per-file declarations from the 2026c headers */
/* extra error code for new dta version not yet supported */
private enum E_NOTYET = GretlError.max + 1;

/* see http://www.stata.com/help.cgi?dta */

/* Stata versions */
private enum VERSION_5   = 105;
private enum VERSION_6   = 108;
private enum VERSION_7   = 110;
private enum VERSION_7SE = 111;
private enum VERSION_8   = 113;
private enum VERSION_10  = 114;
private enum VERSION_12  = 115;
private enum VERSION_13  = 117;
private enum VERSION_14  = 118;

/* old Stata format constants */
private enum STATA_STRINGOFFSET = 0x7f;
private enum STATA_FLOAT  = 'f';
private enum STATA_DOUBLE = 'd';
private enum STATA_LONG   = 'l';
private enum STATA_INT    = 'i';
private enum STATA_BYTE   = 'b';

/* Stata SE format constants */
private enum STATA_SE_STRINGOFFSET = 0;

private enum StataSeType {
    byte_ = 251,
    int_,
    long_,
    float_,
    double_
}

private enum Stata13Type {
    double_ = 65526,
    float_,
    long_,
    int_,
    byte_
}

/* Stata 13+ format constants */
private enum STATA_13_STRL = 32768;
private enum STATA_STRF_MAX = 2045;

/* size of the stack buffer used to hold the dataset label in
   read_dta_117_data(): the length read from the file must not
   be allowed to exceed this
*/
private enum DTA_LABEL_BUFSIZE = 322;

private enum STATA_FLOAT_MAX  = 1.701e+38;
private enum STATA_DOUBLE_MAX = 8.988e+307;
private enum STATA_LONG_MAX   = 2147483620;
private enum STATA_INT_MAX    = 32740;

private enum STATA_INT_CUT = 32767;

/* Stata missing value codes */
private bool STATA_FLOAT_NA(double x)  { return x > STATA_FLOAT_MAX; }
private bool STATA_DOUBLE_NA(double x) { return x > STATA_DOUBLE_MAX; }
private bool STATA_LONG_NA(int i)      { return i > STATA_LONG_MAX; }
private bool STATA_INT_NA(int i)       { return i > STATA_INT_MAX; }
private bool STATA_BYTE_NA(int b, int v) { return (v<8 && b==127) || b>=101; }

/* Add this to Stata daily date values to get epoch days */
private enum STATA_DAY_OFFSET = 715523;

private enum NA_INT = -999;

/* glib's G_LITTLE_ENDIAN / G_BIG_ENDIAN */
private enum Endian {
    little = 1234,
    big = 4321
}

version (BigEndian) {
    private enum HOST_ENDIAN = Endian.big;
} else {
    private enum HOST_ENDIAN = Endian.little;
}

private ulong u64_from_be(ulong u) { version (LittleEndian) return bswap(u); else return u; }
private ulong u64_from_le(ulong u) { version (LittleEndian) return u; else return bswap(u); }
private uint u32_from_be(uint u) { version (LittleEndian) return bswap(u); else return u; }
private uint u32_from_le(uint u) { version (LittleEndian) return u; else return bswap(u); }
private ushort u16_from_be(ushort u)
{
    version (LittleEndian) return cast(ushort) (((u >> 8) & 0xff) | ((u & 0xff) << 8)); else return u;
}
private ushort u16_from_le(ushort u)
{
    version (LittleEndian) return u; else return cast(ushort) (((u >> 8) & 0xff) | ((u & 0xff) << 8));
}

/* it's convenient to have these as file-scope globals */
__gshared private int stata_version;
__gshared private int stata_OLD;
__gshared private int stata_SE;
__gshared private int stata_13;
__gshared private int stata_endian;
__gshared private int swapends;

private struct DtaTable {
    int version_;       /* 117 or higher */
    int nvar;           /* number of variables */
    int nobs;           /* number of observations */
    int dlabelpos;      /* position of dataset label */
    int dlabellen;      /* dataset label length, bytes */
    int timepos;        /* position of time-stamp */
    long vtype_pos;     /* position of variable type info */
    long vname_pos;     /* position of variable names */
    long vfmt_pos;      /* position of variable formats */
    long vallblnam_pos; /* position of value label names */
    long varlabel_pos;  /* position of variable labels */
    long data_pos;      /* position of actual data */
    long strl_pos;      /* position of StrL variables */
    long vallabel_pos;  /* position of value labels */
}

private void stata_globals_init()
{
    stata_version = 0;
    stata_OLD = 0;
    stata_SE = 0;
    stata_13 = 0;
    stata_endian = 0;
    swapends = 0;
}

private DtaTable * dta_table_new(int * err)
{
    DtaTable * dtab = cast(DtaTable *) calloc(1, DtaTable.sizeof);

    if (dtab == null) {
        *err = GretlError.alloc;
    }

    return dtab;
}

private void bin_error(int * err, const(char) * func)
{
    fprintf(stderr, "binary read error in %s()\n", func);
    *err = 1;
}

private long stata_read_int64(FILE * fp, int * err)
{
    long i;

    if (fread(&i, i.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_int64");
        return NA_INT;
    }

    if (swapends) {
        if (stata_endian == Endian.big) {
            i = cast(long) u64_from_be(cast(ulong) i);
        } else {
            i = cast(long) u64_from_le(cast(ulong) i);
        }
    }

    return i;
}

private ulong stata_read_uint64(FILE * fp, int * err)
{
    ulong u;

    if (fread(&u, u.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_uint64");
        return 0;
    }

    if (swapends) {
        if (stata_endian == Endian.big) {
            u = u64_from_be(u);
        } else {
            u = u64_from_le(u);
        }
    }

    return u;
}

private uint stata_read_uint32(FILE * fp, int * err)
{
    uint u;

    if (fread(&u, u.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_uint32");
        return 0;
    }

    if (swapends) {
        if (stata_endian == Endian.big) {
            u = u32_from_be(u);
        } else {
            u = u32_from_le(u);
        }
    }

    return u;
}

/* 4-byte signed int */

private int stata_read_int32(FILE * fp, int naok, int * err)
{
    int i;

    if (fread(&i, i.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_int32");
        return NA_INT;
    }

    if (swapends) {
        reverse_int(i);
    }

    return (STATA_LONG_NA(i) & !naok)? NA_INT : i;
}

private int stata_read_signed_byte(FILE * fp, int naok, int * err)
{
    byte b;
    int ret;

    if (fread(&b, 1, 1, fp) != 1) {
        bin_error(err, "stata_read_signed_byte");
        ret = NA_INT;
    } else {
        ret = cast(int) b;

        if (!naok && STATA_BYTE_NA(b, stata_version)) {
            ret = NA_INT;
        }
    }

    return ret;
}

private int stata_read_byte(FILE * fp, int * err)
{
    ubyte u;

    if (fread(&u, 1, 1, fp) != 1) {
        bin_error(err, "stata_read_byte");
        return NA_INT;
    }

    return cast(int) u;
}

/* 2-byte signed int */

private int stata_read_short(FILE * fp, int naok, int * err)
{
    uint first, second;
    int s;

    first = stata_read_byte(fp, err);
    second = stata_read_byte(fp, err);

    if (stata_endian == Endian.big) {
        s = (first << 8) | second;
    } else {
        s = (second << 8) | first;
    }

    if (s > STATA_INT_CUT) {
        /* ?? */
        s -= 65536;
    }

    return (STATA_INT_NA(s) && !naok)? NA_INT : s;
}

/* 2-byte unsigned int */

private ushort stata_read_uint16(FILE * fp, int * err)
{
    ushort u;

    if (fread(&u, u.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_uint16");
        return cast(ushort) NA_INT;
    }

    if (swapends) {
        if (stata_endian == Endian.big) {
            u = u16_from_be(u);
        } else {
            u = u16_from_le(u);
        }
    }

    return u;
}

private double stata_read_double(FILE * fp, int * err)
{
    double d = 0;

    if (fread(&d, d.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_double");
    }

    if (swapends) {
        reverse_double(d);
    }

    return (STATA_DOUBLE_NA(d))? NADBL : d;
}

private double stata_read_float(FILE * fp, int * err)
{
    float f = 0;

    if (fread(&f, f.sizeof, 1, fp) != 1) {
        bin_error(err, "stata_read_float");
    }

    if (swapends) {
        reverse_float(f);
    }

    return (STATA_FLOAT_NA(f))? NADBL : cast(double) f;
}

/* In Stata files some strings are NUL-terminated but some may not be.
   This function simply reads @nc bytes from @fp. If the string may
   not be properly terminated the caller should ensure a NUL is
   appended to @buf.
*/

private void stata_read_string(FILE * fp, int nc, char * buf, int * err)
{
    if (fread(buf, 1, nc, fp) != nc) {
        bin_error(err, "stata_read_string");
    }
}

private int stata_get_version_and_namelen(ubyte u, int * vnamelen)
{
    int err = 0;

    switch (u) {
    case VERSION_5:
        stata_version = 5;
        *vnamelen = 8;
        break;
    case VERSION_6:
        stata_version = 6;
        *vnamelen = 8;
        break;
    case VERSION_7:
        stata_version = 7;
        *vnamelen = 32;
        break;
    case VERSION_7SE:
        stata_version = 7;
        stata_SE = 1;
        *vnamelen = 32;
        break;
    case VERSION_8:
        stata_version = 8; /* versions >= 8 automatically use 'SE' format */
        stata_SE = 1;
        *vnamelen = 32;
        break;
    case VERSION_10:
        stata_version = 10;
        stata_SE = 1;
        *vnamelen = 32;
        break;
    case VERSION_12:
        stata_version = 12;
        stata_SE = 1;
        *vnamelen = 32;
        break;
    default:
        err = 1;
    }

    stata_OLD = !stata_SE;

    return err;
}

private int stata_get_endianness(FILE * fp, int * err)
{
    int i = cast(int) stata_read_byte(fp, err);

    return (i == 0x01)? Endian.big : Endian.little;
}

/* the following are _not_ to be used with dta >= 117 */

private bool stata_type_float(int t)
{
    return (stata_SE && t == StataSeType.float_) || (stata_OLD && t == STATA_FLOAT);
}

private bool stata_type_double(int t)
{
    return (stata_SE && t == StataSeType.double_) || (stata_OLD && t == STATA_DOUBLE);
}

private bool stata_type_long(int t)
{
    return (stata_SE && t == StataSeType.long_) || (stata_OLD && t == STATA_LONG);
}

private bool stata_type_int(int t)
{
    return (stata_SE && t == StataSeType.int_) || (stata_OLD && t == STATA_INT);
}

private bool stata_type_byte(int t)
{
    return (stata_SE && t == StataSeType.byte_) || (stata_OLD && t == STATA_BYTE);
}

private bool stata_type_string(int t)
{
    return (stata_SE && t <= 244) || (stata_OLD && t >= STATA_STRINGOFFSET);
}

private int check_old_variable_types(FILE * fp, int * types,
                                     int nvar, int * nsv,
                                     Prn * prn)
{
    int i, err = 0;

    *nsv = 0;

    for (i=0; i<nvar && !err; i++) {
        ubyte u = cast(ubyte) stata_read_byte(fp, &err);

        types[i] = u;
        if (stata_type_float(u) || stata_type_double(u)) {
            pprintf(prn, "variable %d: float type\n", i+1);
        } else if (stata_type_long(u)) {
            pprintf(prn, "variable %d: long type\n", i+1);
        } else if (stata_type_int(u)) {
            pprintf(prn, "variable %d: int type\n", i+1);
        } else if (stata_type_byte(u)) {
            pprintf(prn, "variable %d: byte type\n", i+1);
        } else if (stata_type_string(u)) {
            pprintf(prn, "variable %d: string type (%d)\n", i+1, cast(int) u);
            *nsv += 1;
        } else {
            pputs(prn, _("unknown data type").ptr);
            pputc(prn, '\n');
            err = 1;
        }
    }

    return err;
}

private int check_new_variable_types(FILE * fp, int * types,
                                     int nvar, int * nsv,
                                     Prn * prn)
{
    int i, err = 0;

    *nsv = 0;

    for (i=0; i<nvar && !err; i++) {
        ushort u = stata_read_uint16(fp, &err);

        types[i] = u;
        if (u == Stata13Type.float_ || u == Stata13Type.double_) {
            pprintf(prn, "variable %d: float type\n", i+1);
        } else if (u == Stata13Type.long_) {
            pprintf(prn, "variable %d: long type\n", i+1);
        } else if (u == Stata13Type.int_) {
            pprintf(prn, "variable %d: int type\n", i+1);
        } else if (u == Stata13Type.byte_) {
            pprintf(prn, "variable %d: byte type\n", i+1);
        } else if (u <= STATA_STRF_MAX) {
            pprintf(prn, "variable %d: string type (fixed length %d)\n",
                    i+1, cast(int) u);
            *nsv += 1;
        } else if (u == STATA_13_STRL) {
            pprintf(prn, "variable %d: long string type\n", i+1);
            *nsv += 1;
        } else {
            pputs(prn, _("unknown data type").ptr);
            pputc(prn, '\n');
            err = 1;
        }
    }

    return err;
}

/* mechanism for handling (coding) non-numeric variables */

private GretlStringTable * dta_make_string_table(int * types, int nvar, int ncols)
{
    GretlStringTable * st;
    int * list;
    int i, j;

    list = gretl_list_new(ncols);
    if (list == null) {
        return null;
    }

    j = 1;

    if (stata_13) {
        for (i=0; i<nvar && j<=list[0]; i++) {
            if (types[i] <= STATA_STRF_MAX ||
                types[i] == STATA_13_STRL) {
                list[j++] = i + 1;
            }
        }
    } else {
        for (i=0; i<nvar && j<=list[0]; i++) {
            if (!stata_type_float(types[i]) &&
                !stata_type_double(types[i]) &&
                !stata_type_long(types[i]) &&
                !stata_type_int(types[i]) &&
                !stata_type_byte(types[i])) {
                list[j++] = i + 1;
            }
        }
    }

    st = gretl_string_table_new(list);

    free(list);

    return st;
}

/* for use with dta <= 117 */

private gchar * recode_stata_string(const(char) * s)
{
    const(gchar) * cset;
    gchar * tr;
    gsize sz;

    /* most likely to be Windows CP 1252 */
    tr = g_convert(s, -1, "UTF-8", "CP1252", null, &sz, null);

    if (tr == null && !g_get_charset(&cset)) {
        /* fallback: try recoding from current locale, if
           the locale is not UTF-8
        */
        tr = g_locale_to_utf8(s, -1, null, &sz, null);
    }

    return tr;
}

private void save_dataset_info(DatasetRaw * dinfo, const(char) * label, const(char) * stamp)
{
    int dlen = cast(int) strlen(stamp);
    gchar * tr = null;

    if (*label != '\0') {
        if (g_utf8_validate(label, -1, null)) {
            tr = g_strdup(label);
        } else {
            tr = recode_stata_string(label);
        }
    }

    if (tr != null) {
        dlen += cast(int) strlen(tr);
    }

    if (dlen > 0) {
        dinfo.descrip = cast(char *) malloc(dlen + 2);
    }

    if (dinfo.descrip != null) {
        *dinfo.descrip = '\0';
        if (tr != null) {
            strcat(dinfo.descrip, tr);
            strcat(dinfo.descrip, "\n");
        }
        strcat(dinfo.descrip, stamp);
    }

    g_free(tr);
}

private int try_fix_varname(char * name)
{
    char[VNAMELEN] test = 0;
    int err = 0;

    test[0] = '\0';

    if (*name == '_') {
        strcat(test.ptr, "x");
        strncat(test.ptr, name, VNAMELEN - 2);
    } else {
        strncat(test.ptr, name, VNAMELEN - 2);
        strcat(test.ptr, "1");
    }

    err = check_varname(test.ptr);
    if (!err) {
        fprintf(stderr, "Warning: illegal name '%s' changed to '%s'\n",
                name, test.ptr);
        strcpy(name, test.ptr);
    } else {
        /* get the right error message in place */
        check_varname(name);
    }

    return err;
}

/* Check a Stata series that supposedly contains daily dates
   in their 1960-based format. See if we can figure out an
   appropriate daily frequency (5, 6 or 7) and in addition
   register whether the time series is complete (for the given
   frequency) or has some implicit missing values.
*/

private int stata_daily_pd(DatasetRaw * dset, int tv, int * complete)
{
    double * d = dset.Z[tv];
    int T = dset.n;
    int[4] dfreq = 0;
    int wdelta = 0;
    int mdelta = 0;
    int qdelta = 0;
    int t, delta;
    int pd = 0;

    /* find frequency of daily deltas from 1 to 4 days;
       but also check for possible weekly, monthly
       or quarterly daily deltas
    */

    for (t=1; t<T; t++) {
        delta = cast(int) d[t] - cast(int) d[t-1];
        if (delta <= 0) {
            /* not right! */
            pd = -1;
            break;
        } else if (delta <= 4) {
            dfreq[delta-1] += 1;
        } else if (delta == 7) {
            wdelta++;
        } else if (delta >= 28 && delta <= 31) {
            mdelta++;
        } else if (delta >= 90 && delta <= 92) {
            qdelta++;
        }
    }

    if (pd == 0) {
        if (wdelta == T - 1) {
            pd = 52;
            *complete = 1;
        } else if (mdelta == T - 1) {
            pd = 12;
        } else if (qdelta == T - 1) {
            pd = 4;
        } else if (dfreq[0] == T - 1) {
            /* all days consecutive */
            pd = 7;
            *complete = 1;
        } else {
            /* heuristic for 5- or 6-day daily data: "most" of
               the deltas should be 1 day
            */
            if (dfreq[0] / cast(double) T > 0.6) {
                if (dfreq[1] > dfreq[2]) {
                    /* skipping one day per week, in general? */
                    *complete = (dfreq[0] + dfreq[1] == T - 1);
                    pd = 6;
                } else if (dfreq[2] > dfreq[1]) {
                    /* skipping two days per week, in general? */
                    *complete = (dfreq[0] + dfreq[2] == T - 1);
                    pd = 5;
                }
            }
        }
    }

    if (pd == 12 || pd == 4) {
        long ed0 = cast(int) d[0] + STATA_DAY_OFFSET;
        int yr = 0, mo = 0, day = 0;

        ymd_bits_from_epoch_day(cast(uint) ed0, &yr, &mo, &day);
        if (day != 1) {
            /* insist that months and quarters start on
               the first day? */
            pd = 0;
        } else if (pd == 12) {
            sprintf(dset.stobs.ptr, "%d:%02d", yr, mo);
        } else if (pd == 4) {
            int q = mo < 4 ? 1 : mo < 7 ? 2 : mo < 10 ? 3 : 4;

            sprintf(dset.stobs.ptr, "%d:%d", yr, q);
        }
    }

    return pd;
}

/* In case we get what appear to be valid daily dates, but
   they are not complete, add observation markers to the
   dataset and write in the specific daily dates.
*/

private int add_daily_labels(DatasetRaw * dset, int tv)
{
    int err = dataset_allocate_obs_markers(dset);

    if (!err) {
        int t, y, m, d;
        long ed;

        for (t=0; t<dset.n && !err; t++) {
            ed = cast(long) (dset.Z[tv][t] + STATA_DAY_OFFSET);
            err = ymd_bits_from_epoch_day(cast(uint) ed, &y, &m, &d);
            if (!err) {
                sprintf(dset.S[t], "%04d-%02d-%02d", y, m, d);
            }
        }
        if (err) {
            dataset_destroy_obs_markers(dset);
        } else {
            dset.markers = cast(char) DatasetMarkerType.daily_date_strings;
        }
    }

    return err;
}

/* Try using Stata's "date formats" to reconstruct time series
   information. (Stata dates are all zero at the start of 1960.)
*/

private int set_time_info(DatasetRaw * dset, int tv, int pd)
{
    int t1 = cast(int) dset.Z[tv][0];

    dset.stobs[0] = '\0';

    if (pd == 12) {
        int y = (t1 / 12) + 1960;
        int m = t1 % 12 + 1;

        sprintf(dset.stobs.ptr, "%d:%02d", y, m);
    } else if (pd == 4) {
        int y = (t1 / 4) + 1960;
        int q = t1 % 4 + 1;

        sprintf(dset.stobs.ptr, "%d:%d", y, q);
    } else if (pd == 5) {
        /* daily dates present, but may not really be daily */
        int complete = 0;
        int err = 0;

        pd = stata_daily_pd(dset, tv, &complete);

        if (pd == 12 || pd == 4) {
            /* got excess precision: handled above */
        } else if (pd > 0) {
            long ed0 = t1 + STATA_DAY_OFFSET;
            char * ymd = ymd_extended_from_epoch_day(cast(uint) ed0, 0, &err);

            if (!err) {
                strcpy(dset.stobs.ptr, ymd);
                free(ymd);
                if (!complete) {
                    /* let's add observation labels */
                    err = add_daily_labels(dset, tv);
                    fprintf(stderr, "add_daily_labels: err = %d\n", err);
                }
            }
            if (err) {
                /* scrub it */
                dset.stobs[0] = '\0';
                pd = 1;
            }
        }
    } else {
        int y = t1 + 1960;

        if (y > 2050) {
            /* Can't really be a starting year? */
        } else {
            sprintf(dset.stobs.ptr, "%d", y);
        }
    }

    if (dset.stobs[0] != '\0') {
        dset.pd = pd;
        printf("starting obs seems to be %s\n", dset.stobs.ptr);
        dset.structure = TsCodes.time_series;
        dset.sd0 = get_date_x(dset.pd, dset.stobs.ptr);
        ntolabel(dset.endobs.ptr, dset.n - 1, dset);
    } else {
        dset.pd = 1;
    }

    return 0;
}

private int dta_value_labels_setup(Prn ** pprn, GretlStringTable ** pst)
{
    int err = 0;

    *pprn = gretl_print_new(PrnType.buffer, &err);

    if (*pprn != null) {
        if (*pst == null) {
            *pst = gretl_string_table_new(null);
            if (*pst == null) {
                gretl_print_destroy(*pprn);
                *pprn = null;
            }
        }
    }

    return err;
}

private int push_label_info(int ** pv, char *** pS, int v, const(char) * lname)
{
    int n = (*pv == null)? 0 : (*pv)[0];
    int err = 0;

    if (gretl_list_append_term(pv, v) == null) {
        err = GretlError.alloc;
    } else {
        err = strings_array_add(pS, &n, lname);
    }

    return err;
}

private int label_array_header(const(int) * list, char ** names,
                               const(char) * lname, const(DatasetRaw) * dset,
                               Prn * prn)
{
    int i, v = 0, n = 0;

    for (i=1; i<=list[0]; i++) {
        if (!strcmp(names[i-1], lname)) {
            v = list[i];
            n++;
        }
    }

    if (n == 0) {
        return 0;
    }

    if (n == 1) {
        pprintf(prn, _("\nValue -> label mappings for variable %d (%s)\n").ptr,
                v, dset.varname[v]);
    } else {
        pprintf(prn, _("\nValue -> label mappings for the following %d variables\n").ptr, n);
        for (i=1; i<=list[0]; i++) {
            if (!strcmp(names[i-1], lname)) {
                v = list[i];
                pprintf(prn, " %3d (%s)\n", v, dset.varname[v]);
            }
        }
    }

    return 1;
}

private int stata_seek(FILE * fp, long offset, int whence)
{
    if (fseeko(fp, offset, whence) < 0) {
        fputs("binary read error in stata_seek()\n", stderr);
        return GretlError.data;
    } else {
        return 0;
    }
}

/* The following pertains to numeric variables that have value
   labels attached, as opposed to variables of string type.
*/

private int process_value_labels(FILE * fp, DatasetRaw * dset, int j,
                                 int * lvars, char ** lnames, int namelen,
                                 GretlStringTable ** pst,
                                 Prn ** pst_prn, Prn * vprn)
{
    Prn * st_prn;
    double * level = null;
    char * txt = null;
    int * off = null;
    char[130] buf = 0;
    int nlabels, totlen;
    int i, err = 0;

    if (stata_13) {
        int llen = stata_read_int32(fp, 0, &err);

        pprintf(vprn, "labels %d: value_label_table = %d bytes\n", j, llen);
    }

    stata_read_string(fp, namelen + 1, buf.ptr, &err);
    if (!err) {
        /* the field is namelen+1 bytes; don't trust that the file
           put a NUL in the last one
        */
        buf[namelen] = '\0';
    }
    pprintf(vprn, "labels %d: (namelen=%d) name = '%s'\n", j, namelen, buf.ptr);

    /* padding */
    err = stata_seek(fp, 3, SEEK_CUR);

    if (!err) {
        nlabels = stata_read_int32(fp, 1, &err);
        totlen = stata_read_int32(fp, 1, &err);
        if (nlabels <= 0 || totlen <= 0) {
            return 0;
        }
    }

    if (!err && *pst_prn == null) {
        err = dta_value_labels_setup(pst_prn, pst);
    }

    if (err) {
        return err;
    }

    st_prn = *pst_prn;

    off = cast(int *) malloc(nlabels * (*off).sizeof);
    if (off == null) {
        return GretlError.alloc;
    }

    level = cast(double *) malloc(nlabels * (*level).sizeof);
    if (level == null) {
        free(off);
        return GretlError.alloc;
    }

    label_array_header(lvars, lnames, buf.ptr, dset, st_prn);

    for (i=0; i<nlabels && !err; i++) {
        off[i] = stata_read_int32(fp, 1, &err);
    }

    for (i=0; i<nlabels && !err; i++) {
        level[i] = cast(double) stata_read_int32(fp, 0, &err);
        pprintf(vprn, " level %d = %g\n", i, level[i]);
    }

    if (!err) {
        txt = cast(char *) calloc(totlen, 1);
        if (txt == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        stata_read_string(fp, totlen, txt, &err);
        if (!err) {
            /* guard against a missing NUL at the very end of the
               text block: force one so a valid-looking offset can
               never run off the end of the buffer
            */
            txt[totlen - 1] = '\0';
        }
    }

    for (i=0; i<nlabels && !err; i++) {
        const(char) * vlabel;

        if (off[i] < 0 || off[i] >= totlen) {
            /* @off[i] is an unvalidated on-disk value: without this
               check it can point anywhere relative to @txt
            */
            fprintf(stderr, "value label: implausible text offset "
                    ~ "%d (totlen=%d)\n", off[i], totlen);
            err = GretlError.data;
            break;
        }

        vlabel = txt + off[i];

        pprintf(vprn, " label %d = '%s'\n", i, vlabel);
        if (g_utf8_validate(vlabel, -1, null)) {
            pprintf(st_prn, "%10g -> '%s'\n", level[i], vlabel);
        } else {
            gchar * tr = recode_stata_string(vlabel);

            if (tr != null) {
                pprintf(st_prn, "%10g -> '%s'\n", level[i], tr);
                g_free(tr);
            } else {
                pprintf(st_prn, "%10g -> 'unknown'\n", level[i]);
            }
        }
    }

    free(off);
    free(level);
    free(txt);

    return err;
}

private int process_stata_varname(FILE * fp, char * buf, int namelen,
                                  DatasetRaw * dset, int v, Prn * vprn)
{
    int err = 0;

    stata_read_string(fp, namelen + 1, buf, &err);

    if (!err) {
        /* the field is namelen+1 bytes; don't trust that the file
           put a NUL in the last one -- @buf is sized namelen+1 by
           every caller
        */
        buf[namelen] = '\0';
    }

    if (!err) {
        if (gretl_namechar_spn(buf) == strlen(buf)) {
            pprintf(vprn, "variable %d: name = '%s'\n", v, buf);
        } else {
            pprintf(vprn, "variable %d: raw name = '%s'\n", v, buf);
        }
        if (namelen == 32) {
            /* dta 117: try to fix possible bad encoding */
            iso_to_ascii(buf);
        } else {
            /* dta > 117: UTF-8 varnames of up to 128 bytes */
            asciify_utf8_varname(buf);
        }
        if (*buf == '\0') {
            sprintf(buf, "v%d", v);
        } else {
            err = check_varname(buf);
            if (err) {
                err = try_fix_varname(buf);
            }
        }
    }

    if (!err) {
        dset.varname[v][0] = '\0';
        strncat(dset.varname[v], buf, VNAMELEN - 1);
    }

    return err;
}

/* If we get a format that seems to represent date/time,
   write the associated periodicity to @pd and record
   the series ID in @tvar for later use.
*/

private void process_stata_format(char * buf, int v,
                                  int * pd, int * tvar,
                                  Prn * vprn)
{
    if (*buf != '\0' && buf[strlen(buf)-1] != 'g') {
        pprintf(vprn, "variable %d: format = '%s'\n", v, buf);
        if (!strcmp(buf, "%tm")) {
            *pd = 12;
            *tvar = v;
        } else if (!strcmp(buf, "%tq")) {
            *pd = 4;
            *tvar = v;
        } else if (!strcmp(buf, "%ty")) {
            *pd = 1;
            *tvar = v;
        } else if (!strcmp(buf, "%td")) {
            *pd = 5; /* may be revised later */
            *tvar = v;
        }
    }
}

private int process_stata_varlabel(char * label, DatasetRaw * dset,
                                   int v, Prn * vprn)
{
    int err = 0;

    pprintf(vprn, "variable %d: label = '%s'\n", v, label);

    if (g_utf8_validate(label, -1, null)) {
        if (stata_version > 13) {
            gretl_utf8_truncate_b(label, MAXLABEL-1);
        }
        series_set_label(dset, v, label);
    } else {
        gchar * tr = recode_stata_string(label);

        if (tr != null) {
            series_set_label(dset, v, tr);
            g_free(tr);
        } else {
            err = GretlError.data;
        }
    }

    return err;
}

/* @buf must be of length @bufsize; if @len > @bufsize - 1 the read is
   truncated to @bufsize - 1 bytes (and nul-termination is handled), but
   regardless we move the read position to the end of the full-length
   region of the stream as given by @len.
*/

private int stata_read_buffer(char * buf, int bufsize,
                              int len, FILE * fp)
{
    int n = bufsize - 1;
    int err = 0;

    *buf = '\0';

    if (len < 0) {
        /* the on-disk length field is unsigned (guint32); a value
           greater than INT_MAX will have come through as negative --
           reject it rather than passing it to fread()/array-indexing
           as a huge size_t
        */
        fputs("stata_read_buffer: implausible (negative) length\n", stderr);
        return GretlError.data;
    }

    if (len > n) {
        stata_read_string(fp, n, buf, &err);
        buf[n] = '\0';
        if (stata_version > 13) {
            /* beware broken UTF-8! */
            int pos = n - 1;

            while (!g_utf8_validate(buf, -1, null)) {
                buf[pos--] = '\0';
            }
        }
        err = stata_seek(fp, len - n, SEEK_CUR);
    } else {
        stata_read_string(fp, len, buf, &err);
        buf[len] = '\0';
    }

    return err;
}

/* The following pertains to a non-numeric variable whose
   values are actually strings (as opposed to a numeric
   variable with value labels attached). According to the
   Stata 117 dta spec such string values (@buf) should be
   ASCII but in practice this may not be the case. It
   seems that in dta 118 they may be UTF-8.
*/

private void process_string_value(char * buf, GretlStringTable * st,
                                  DatasetRaw * dset, int v, int t,
                                  Prn * prn)
{
    if (st != null && strcmp(buf, ".")) {
        int ix = 0;

        if (g_utf8_validate(buf, -1, null)) {
            ix = gretl_string_table_index(st, buf, v, 0, prn);
        } else {
            gchar * tr = recode_stata_string(buf);

            if (tr != null) {
                ix = gretl_string_table_index(st, tr, v, 0, prn);
                g_free(tr);
            }
        }

        if (ix > 0) {
            dset.Z[v][t] = ix;
            if (t == 0) {
                series_set_discrete(dset, v, 1);
            }
        }
    }
}

private void maybe_fix_varlabel_pos(FILE * fp, DtaTable * dtab)
{
    /* we're at the end of value_label_names: maybe we have
       variable labels up next?
    */
    char[32] test = 0;
    int err = 0;

    stata_read_string(fp, 20, test.ptr, &err);
    test[20] = '\0';

    if (!strcmp(test.ptr, "</value_label_names>")) {
        /* so far so good */
        stata_read_string(fp, 17, test.ptr, &err);
        test[17] = '\0';
        if (!strcmp(test.ptr, "<variable_labels>")) {
            /* got it */
            dtab.varlabel_pos = ftello(fp);
        }
    }
}

private int find_strl_value(char * buf, int bufsize,
                            int vnum, int obs, FILE * fp)
{
    char[4] test = 0;
    uint v, len;
    ulong o;
    int err = 0;

    stata_read_string(fp, 3, test.ptr, &err);
    test[3] = '\0';
    if (!err && strcmp(test.ptr, "GSO")) {
        fprintf(stderr, "empty or bad strls block\n");
        err = GretlError.data;
    }

    while (!err) {
        v = stata_read_uint32(fp, &err);
        o = stata_read_uint64(fp, &err);
        /* next is 't', 'binary vs ASCII' -- not needed? */
        stata_read_byte(fp, &err);
        len = stata_read_uint32(fp, &err);
        if (!err) {
            if (v == vnum && o == obs) {
                /* found what we're looking for */
                err = stata_read_buffer(buf, bufsize, len, fp);
                break;
            } else {
                /* try for the next GSO */
                err = stata_seek(fp, len, SEEK_CUR);
                stata_read_string(fp, 3, test.ptr, &err);
                test[3] = '\0';
                if (strcmp(test.ptr, "GSO")) {
                    if (strcmp(test.ptr, "</s")) {
                        fprintf(stderr, "broken strls zone?\n");
                    }
                    break;
                }
            }
        }
    }

    return err;
}

private int get_strl_data(char * buf, int bufsize,
                          DtaTable * dtab, FILE * fp)
{
    char[8] vobytes = 0;
    int err = 0;

    *buf = '\0';
    vobytes[0] = '\0';
    stata_read_string(fp, 8, vobytes.ptr, &err);

    if (!err) {
        long o = 0;
        ushort v;
        version (BigEndian) {
            char * targ = cast(char *) &o + 2;
        } else {
            char * targ = cast(char *) &o;
        }

        memcpy(&v, vobytes.ptr, 2);
        memcpy(targ, vobytes.ptr + 2, 6);

        if (o > INT_MAX) {
            fprintf(stderr, "strL: obs reference too big!\n");
            err = GretlError.data;
        } else if (v == 0 && o == 0) {
            /* empty string */
        } else {
            long thispos = ftello(fp);

            err = stata_seek(fp, dtab.strl_pos, SEEK_SET);
            if (!err) {
                err = find_strl_value(buf, bufsize, cast(int) v, cast(int) o, fp);
            }
            if (!err) {
                /* get back in position */
                err = stata_seek(fp, thispos, SEEK_SET);
            }
        }
    }

    return err;
}

private void set_discreteness_from_types(DatasetRaw * dset,
                                         const(int) * types,
                                         int newtypes)
{
    int i;

    for (i=1; i<dset.v; i++) {
        int discrete = 0;
        int t = types[i-1];

        if (newtypes) {
            if (t != Stata13Type.float_ && t != Stata13Type.double_ &&
                t != Stata13Type.long_) {
                discrete = 1;
            }
        } else {
            if (!stata_type_float(t) && !stata_type_double(t) &&
                !stata_type_long(t)) {
                discrete = 1;
            }
        }
        if (discrete) {
            series_set_discrete(dset, i, 1);
        }
    }
}

/* main reader for dta format 117+ */

private int read_dta_117_data(FILE * fp, DatasetRaw * dset,
                              GretlStringTable ** pst,
                              DtaTable * dtab, Prn * prn,
                              Prn * vprn)
{
    int i, j, t;
    int namelen = 32;
    int fmtlen = 49;
    int vlabellen = 81;
    int nvar = dset.v - 1;
    int nsv = 0;
    int pd = 0, tvar = -1;
    char[DTA_LABEL_BUFSIZE] label = 0; /* dataset label */
    char[130] aname = 0; /* variable names */
    char[60] c60 = 0;    /* misc strings */
    int * types = null;
    int * lvars = null;
    char ** lnames = null;
    char[256] buf = 0;
    int st_err = 0;
    int err = 0;

    if (dtab.version_ > 117) {
        /* dta > 117 supports UTF-8 strings */
        fmtlen = 57;
        namelen = 128;
        vlabellen = 321;
    }

    label[0] = c60[0] = '\0';

    if (dtab.dlabellen > 0) {
        err = stata_seek(fp, dtab.dlabelpos, SEEK_SET);
        if (!err) {
            stata_read_string(fp, dtab.dlabellen, label.ptr, &err);
        }
        if (!err) {
            label[dtab.dlabellen] = '\0';
            pprintf(vprn, "dataset label: '%s'\n", label.ptr);
        }
    }

    if (!err && dtab.timepos > 0) {
        err = stata_seek(fp, dtab.timepos, SEEK_SET);
        if (!err) {
            stata_read_string(fp, 17, c60.ptr, &err);
        }
        if (!err) {
            c60[17] = '\0';
            pprintf(vprn, "timestamp: '%s'\n", c60.ptr);
        }
    }

    if (err) {
        return err;
    }

    if (label[0] != '\0') {
        save_dataset_info(dset, label.ptr, c60.ptr);
    }

    /** read variable descriptors **/

    /* types */

    types = cast(int *) malloc(nvar * (*types).sizeof);
    if (types == null) {
        return GretlError.alloc;
    }

    err = stata_seek(fp, dtab.vtype_pos, SEEK_SET);

    if (!err) {
        err = check_new_variable_types(fp, types, nvar, &nsv, vprn);
    }

    if (err) {
        free(types);
        return err;
    }

    if (nsv > 0) {
        /* we have 1 or more non-numeric variables */
        *pst = dta_make_string_table(types, nvar, nsv);
    }

    err = stata_seek(fp, dtab.vname_pos, SEEK_SET);

    /* variable names */
    for (i=0; i<nvar && !err; i++) {
        err = process_stata_varname(fp, aname.ptr, namelen,
                                    dset, i+1, vprn);
    }

    if (!err && dtab.vfmt_pos > 0) {
        /* format list (use it to extract time-series info?) */
        err = stata_seek(fp, dtab.vfmt_pos, SEEK_SET);
        for (i=0; i<nvar && !err; i++) {
            stata_read_string(fp, fmtlen, c60.ptr, &err);
            if (!err) {
                /* @fmtlen never exceeds sizeof c60 - 1: the file is
                   not guaranteed to have NUL-terminated the field
                */
                c60[fmtlen] = '\0';
            }
            if (!err && types[i] >= Stata13Type.double_) {
                process_stata_format(c60.ptr, i+1, &pd, &tvar, vprn);
            }
        }
    }

    if (err) {
        free(types);
        return err;
    }

    err = stata_seek(fp, dtab.vallblnam_pos, SEEK_SET);

    /* value-label names: the _names of label formats,
       which are themselves stored later in the file
    */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, namelen + 1, aname.ptr, &err);
        if (!err) {
            /* the field is namelen+1 bytes, its own last byte is
               meant to be the NUL -- don't trust that the file
               actually put one there
            */
            aname[namelen] = '\0';
        }
        if (!err && aname[0] != '\0' && !st_err) {
            pprintf(vprn, "variable %d: value-label name = '%s'\n", i+1, aname.ptr);
            st_err = push_label_info(&lvars, &lnames, i+1, aname.ptr);
        }
    }

    if (!err && dtab.varlabel_pos == 0) {
        /* work around bug in some Stata files */
        maybe_fix_varlabel_pos(fp, dtab);
    }

    if (!err && dtab.varlabel_pos > 0) {
        /* variable descriptive labels */
        err = stata_seek(fp, dtab.varlabel_pos, SEEK_SET);
        for (i=0; i<nvar && !err; i++) {
            stata_read_string(fp, vlabellen, label.ptr, &err);
            if (!err) {
                /* @vlabellen never exceeds sizeof label - 1 */
                label[vlabellen] = '\0';
            }
            if (label[0] != '\0') {
                process_stata_varlabel(label.ptr, dset, i+1, vprn);
            }
        }
    }

    if (err) {
        free(types);
        return err;
    }

    err = stata_seek(fp, dtab.data_pos, SEEK_SET);

    /* actual data values */
    for (t=0; t<dset.n && !err; t++) {
        for (i=0; i<nvar && !err; i++) {
            int ix, v = i + 1;
            int typ = types[i];

            dset.Z[v][t] = NADBL;

            if (typ == Stata13Type.float_) {
                dset.Z[v][t] = cast(double) stata_read_float(fp, &err);
            } else if (typ == Stata13Type.double_) {
                dset.Z[v][t] = stata_read_double(fp, &err);
            } else if (typ == Stata13Type.long_) {
                ix = stata_read_int32(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : cast(double) ix;
            } else if (typ == Stata13Type.int_) {
                ix = stata_read_short(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : cast(double) ix;
            } else if (typ == Stata13Type.byte_) {
                ix = stata_read_signed_byte(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : cast(double) ix;
            } else if (typ <= STATA_STRF_MAX) {
                /* fixed length string */
                err = stata_read_buffer(buf.ptr, cast(int) buf.length, typ, fp);
                if (!err && buf[0] != '\0') {
                    process_string_value(buf.ptr, *pst, dset, v, t, prn);
                }
            } else if (typ == STATA_13_STRL) {
                /* arbitrarily long string! */
                err = get_strl_data(buf.ptr, cast(int) buf.length, dtab, fp);
                if (!err && buf[0] != '\0') {
                    process_string_value(buf.ptr, *pst, dset, v, t, prn);
                }
            }
        }
    }

    if (!err && tvar > 0) {
        set_time_info(dset, tvar, pd);
    }

    if (!err) {
        err = stata_seek(fp, dtab.vallabel_pos, SEEK_SET);
    }

    /* value labels */

    if (!err && !st_err && lvars != null) {
        Prn * st_prn = null;

        for (j=0; j<nvar && !st_err; j++) {
            /* check for closing "</value_labels>" */
            char[16] test = 0;

            if (fread(test.ptr, 1, 15, fp) < 15 || !strcmp(test.ptr, "</value_labels>")) {
                pprintf(vprn, _("breaking on end of value labels\n").ptr);
                break;
            }
            /* otherwise it's supposed to be "<lbl>" */
            st_err = stata_seek(fp, -10, SEEK_CUR);
            if (!st_err) {
                st_err = process_value_labels(fp, dset, j, lvars, lnames,
                                              namelen, pst, &st_prn, vprn);
            }
            if (!st_err) {
                /* skip pseudo-XML for next read */
                st_err = stata_seek(fp, strlen("</lbl>"), SEEK_CUR);
            }
        }

        if (st_prn != null) {
            if (!st_err) {
                gretl_string_table_add_extra(*pst, st_prn);
            }
            gretl_print_destroy(st_prn);
        }
    }

    if (!err) {
        set_discreteness_from_types(dset, types, 1);
    }

    free(types);

    if (lvars != null) {
        strings_array_free(lnames, lvars[0]);
        free(lvars);
    }

    return err;
}

/* main reader for dta format <= 115 */

private int read_old_dta_data(FILE * fp, DatasetRaw * dset,
                              GretlStringTable ** pst,
                              int namelen, Prn * prn,
                              Prn * vprn)
{
    int i, j, t, clen;
    int labellen, fmtlen;
    int nvar = dset.v - 1;
    int nsv = 0;
    int soffset, pd = 0, tvar = -1;
    char[81] label = 0;
    char[50] c50 = 0;
    char[33] aname = 0;
    int * types = null;
    int * lvars = null;
    char ** lnames = null;
    char[256] strbuf = 0;
    int st_err = 0;
    int err = 0;

    labellen = (stata_version == 5)? 32 : 81;
    fmtlen = (stata_version < 10)? 12 : 49;
    soffset = (stata_SE)? STATA_SE_STRINGOFFSET : STATA_STRINGOFFSET;

    pprintf(vprn, "Max length of labels = %d\n", labellen);

    /* dataset label: fixed length, nominally NUL-terminated, but
       don't trust a hostile file to have actually done so
    */
    stata_read_string(fp, labellen, label.ptr, &err);
    if (!err) {
        label[labellen-1] = '\0';
    }
    pprintf(vprn, "dataset label: '%s'\n", label.ptr);

    /* timestamp: fixed length, likewise force termination */
    stata_read_string(fp, 18, c50.ptr, &err);
    if (!err) {
        c50[18] = '\0';
    }
    pprintf(vprn, "timestamp: '%s'\n", c50.ptr);

    if (label[0] != '\0' || c50[0] != '\0') {
        save_dataset_info(dset, label.ptr, c50.ptr);
    }

    /** read variable descriptors **/

    /* types */

    types = cast(int *) malloc(nvar * (*types).sizeof);
    if (types == null) {
        return GretlError.alloc;
    }

    err = check_old_variable_types(fp, types, nvar, &nsv, vprn);
    if (err) {
        free(types);
        return err;
    }

    if (nsv > 0) {
        /* we have 1 or more non-numeric variables */
        *pst = dta_make_string_table(types, nvar, nsv);
    }

    /* variable names */
    for (i=0; i<nvar && !err; i++) {
        err = process_stata_varname(fp, aname.ptr, namelen,
                                    dset, i+1, vprn);
    }

    /* sortlist -- not relevant */
    for (i=0; i<2*(nvar+1) && !err; i++) {
        stata_read_byte(fp, &err);
    }

    /* format list (use it to extract time-series info?) */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, fmtlen, c50.ptr, &err);
        if (!err) {
            /* @fmtlen never exceeds sizeof c50 - 1 */
            c50[fmtlen] = '\0';
            process_stata_format(c50.ptr, i+1, &pd, &tvar, vprn);
        }
    }

    /* value-label names: the names of label formats,
       which are themselves stored later in the file
    */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, namelen + 1, aname.ptr, &err);
        if (!err) {
            aname[namelen] = '\0';
        }
        if (aname[0] != '\0' && !st_err) {
            pprintf(vprn, "variable %d: value-label name = '%s'\n", i+1, aname.ptr);
            st_err = push_label_info(&lvars, &lnames, i+1, aname.ptr);
        }
    }

    /* variable descriptive labels */
    for (i=0; i<nvar && !err; i++) {
        stata_read_string(fp, labellen, label.ptr, &err);
        if (!err) {
            label[labellen-1] = '\0';
        }
        if (label[0] != '\0') {
            process_stata_varlabel(label.ptr, dset, i+1, vprn);
        }
    }

    /* variable 'characteristics': not handled */
    if (!err) {
        /* this block ends with a zero byte */
        while (stata_read_byte(fp, &err)) {
            if (stata_version >= 7) {
                /* manual is wrong here */
                clen = stata_read_int32(fp, 1, &err);
            } else {
                clen = stata_read_short(fp, 1, &err);
            }
            for (i=0; i<clen; i++) {
                stata_read_signed_byte(fp, 1, &err);
            }
        }
        if (stata_version >= 7) {
            clen = stata_read_int32(fp, 1, &err);
        } else {
            clen = stata_read_short(fp, 1, &err);
        }
        if (clen != 0) {
            fputs(_("something strange in the file\n"
                    ~ "(Type 0 characteristic of nonzero length)").ptr, stderr);
            fputc('\n', stderr);
        }
    }

    /* actual data values */
    for (t=0; t<dset.n && !err; t++) {
        for (i=0; i<nvar && !err; i++) {
            int ix, v = i + 1;

            dset.Z[v][t] = NADBL;

            if (stata_type_float(types[i])) {
                dset.Z[v][t] = stata_read_float(fp, &err);
            } else if (stata_type_double(types[i])) {
                dset.Z[v][t] = stata_read_double(fp, &err);
            } else if (stata_type_long(types[i])) {
                ix = stata_read_int32(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : ix;
            } else if (stata_type_int(types[i])) {
                ix = stata_read_short(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : ix;
            } else if (stata_type_byte(types[i])) {
                ix = stata_read_signed_byte(fp, 0, &err);
                dset.Z[v][t] = (ix == NA_INT)? NADBL : ix;
            } else {
                /* a string variable */
                clen = types[i] - soffset;
                if (clen > 255) {
                    /* release <= 115: the max is supposed to be 244 */
                    clen = 255;
                }
                stata_read_string(fp, clen, strbuf.ptr, &err);
                strbuf[clen] = '\0';
                if (!err && strbuf[0] != '\0') {
                    process_string_value(strbuf.ptr, *pst, dset, v, t, prn);
                }
            }
        }
    }

    if (!err && tvar > 0) {
        set_time_info(dset, tvar, pd);
    }

    /* value labels */

    if (!err && !st_err && lvars != null && stata_version > 5) {
        Prn * st_prn = null;

        for (j=0; j<nvar; j++) {
            /* first int not needed, use fread directly to trigger EOF */
            int first;
            size_t k = fread(&first, int.sizeof, 1, fp);

            if (k == 0 || feof(fp)) {
                pprintf(vprn, "breaking on feof\n");
                break;
            }
            st_err = process_value_labels(fp, dset, j, lvars, lnames,
                                          namelen, pst, &st_prn, vprn);
        }
        if (st_prn != null) {
            if (!st_err) {
                gretl_string_table_add_extra(*pst, st_prn);
            }
            gretl_print_destroy(st_prn);
        }
    }

    if (!err) {
        set_discreteness_from_types(dset, types, 0);
    }

    free(types);

    if (lvars != null) {
        strings_array_free(lnames, lvars[0]);
        free(lvars);
    }

    return err;
}

/* Save the offsets read from <map>...</map> for later use.
   In general it should be fatal if any of these offsets are
   invalid, however many Stata 13 files have 0 in place of
   the offset to the variable labels block, and we're prepared
   to work around that.
*/

private int dtab_save_offset(DtaTable * dtab, int i, long offset)
{
    if (offset <= 0) {
        int err = 0;

        if (i == 5 || i == 7 || i != 10) {
            /* maybe we can struggle on? */
        } else {
            err = GretlError.data;
        }

        if (i == 7) {
            /* semi-expected */
        } else {
            fprintf(stderr, "map: bad offset (<= 0) for element %d\n", i);
        }

        return err;
    }

    if (i == 2) {
        dtab.vtype_pos = offset + strlen("<variable_types>");
    } else if (i == 3) {
        dtab.vname_pos = offset + strlen("<varnames>");
    } else if (i == 5) {
        dtab.vfmt_pos = offset + strlen("<formats>");
    } else if (i == 6) {
        dtab.vallblnam_pos = offset + strlen("<value_label_names>");
    } else if (i == 7) {
        dtab.varlabel_pos = offset + strlen("<variable_labels>");
    } else if (i == 9) {
        dtab.data_pos = offset + strlen("<data>");
    } else if (i == 10) {
        dtab.strl_pos = offset + strlen("<strls>");
    } else if (i == 11) {
        dtab.vallabel_pos = offset + strlen("<value_labels>");
    }

    return 0;
}

private int stata_get_nobs(FILE * fp, int * err)
{
    ulong N;
    int n = 0;

    N = stata_read_uint64(fp, err);

    /* current gretl can't handle more than INT_MAX
       observations */

    if (!*err && N > INT_MAX) {
        *err = GretlError.data;
    }

    if (!*err) {
        n = cast(int) N;
    }

    return n;
}

/* dta header (format 117, Stata 13):

  file format id     <release>...</release>
  byteorder          <byteorder>...</byteorder>
  # of variables     <K>...</K>
  # of observations  <N>...</N>
  dataset label      <label>...</label>
  datetime stamp     <timestamp>...</timestamp>
*/

private int parse_dta_117_header(FILE * fp, DtaTable * dtab,
                                 Prn * prn, Prn * vprn)
{
    int clen = 0;
    char[8] order = 0;
    char[128] buf = 0;
    size_t b;
    int i, err = 0;

    b = fread(buf.ptr, 1, 96, fp);
    buf[b] = '\0';

    if (sscanf(buf.ptr, "<header><release>%d</release>"
               ~ "<byteorder>%3[^<]</byteorder>",
               &dtab.version_, order.ptr) != 2) {
        err = 1;
    } else {
        if (dtab.version_ != 117 && dtab.version_ != 118) {
            err = 1;
        } else if (!strcmp(order.ptr, "LSF")) {
            stata_endian = Endian.little;
        } else if (!strcmp(order.ptr, "MSF")) {
            stata_endian = Endian.big;
        } else {
            err = 1;
        }
    }

    if (!err) {
        pprintf(prn, "Stata dta version %d, byte-order %s\n",
                dtab.version_, order.ptr);
        swapends = stata_endian != HOST_ENDIAN;
        err = stata_seek(fp, 70, SEEK_SET);
        if (!err) {
            /* read K = number of variables */
            dtab.nvar = stata_read_short(fp, 1, &err);
        }
    }

    if (dtab.version_ > 118) {
        pprintf(prn, _("This dta version not yet supported\n").ptr);
        return E_NOTYET;
    }

    if (!err) {
        /* skip "</K><N>" */
        err = stata_seek(fp, 7, SEEK_CUR);
        if (!err) {
            /* read N = number of observations */
            if (dtab.version_ == 118) {
                dtab.nobs = stata_get_nobs(fp, &err);
            } else {
                dtab.nobs = stata_read_int32(fp, 1, &err);
            }
        }
    }

    if (!err) {
        /* skip "</N><label>" */
        err = stata_seek(fp, 11, SEEK_CUR);
        if (!err) {
            if (dtab.version_ == 118) {
                clen = stata_read_uint16(fp, &err);
            } else {
                clen = stata_read_byte(fp, &err);
            }
            if (!err && clen > 0) {
                if (clen > DTA_LABEL_BUFSIZE - 1) {
                    fprintf(stderr, "dta: implausible dataset label "
                            ~ "length %d\n", clen);
                    err = 1;
                } else {
                    dtab.dlabellen = clen;
                    dtab.dlabelpos = cast(int) ftello(fp);
                    err = stata_seek(fp, clen, SEEK_CUR);
                }
            }
        }
    }

    if (!err) {
        /* skip "</label><timestamp>" */
        err = stata_seek(fp, 19, SEEK_CUR);
        if (!err) {
            clen = stata_read_byte(fp, &err);
            if (!err) {
                if (clen > 0) {
                    dtab.timepos = cast(int) ftello(fp);
                    err = stata_seek(fp, clen, SEEK_CUR);
                }
            }
        }
    }

    if (!err) {
        /* skip "</timestamp></header>" */
        err = stata_seek(fp, 21, SEEK_CUR);
        /* and read 5 bytes */
        if (!err && fread(buf.ptr, 1, 5, fp) != 5) {
            err = 1;
        }
    }

    if (!err) {
        buf[5] = '\0';
        if (strcmp(buf.ptr, "<map>")) {
            err = 1;
        } else {
            long offset;

            for (i=0; i<14 && !err; i++) {
                offset = stata_read_int64(fp, &err);
                if (i > 0 && !err) {
                    err = dtab_save_offset(dtab, i, offset);
                }
            }
        }
    }

    if (!err && (dtab.nvar <= 0 || dtab.nobs < 0)) {
        err = 1;
    }

    if (!err) {
        stata_version = (dtab.version_ == 117)? 13 : 14;
        stata_13 = 1;
    }

    return err;
}

/* for dta versions < 117 */

private int parse_old_dta_header(FILE * fp, int * namelen,
                                 int * nvar, int * nobs,
                                 Prn * prn, Prn * vprn)
{
    ubyte u;
    int err = 0;

    u = cast(ubyte) stata_read_byte(fp, &err); /* release version */

    if (!err) {
        err = stata_get_version_and_namelen(u, namelen);
    }

    if (err) {
        fputs("not a Stata version 5-12 .dta file\n", stderr);
        return err;
    }

    pprintf(prn, _("Stata file version %d\n").ptr, stata_version);

    /* these are file-scope globals */
    stata_endian = stata_get_endianness(fp, &err);
    swapends = stata_endian != HOST_ENDIAN;

    stata_read_byte(fp, &err);              /* filetype -- junk */
    stata_read_byte(fp, &err);              /* padding */
    *nvar = stata_read_short(fp, 1, &err);  /* number of variables */
    *nobs = stata_read_int32(fp, 1, &err);  /* number of observations */

    if (!err && (*nvar <= 0 || *nobs <= 0)) {
        err = 1;
    }

    return err;
}

private int new_style_dta_file(FILE * fp)
{
    char[12] test = 0;

    if (fread(test.ptr, 1, 11, fp) != 11) {
        rewind(fp);
        return 0;
    } else if (!strcmp(test.ptr, "<stata_dta>")) {
        return 1;
    } else {
        rewind(fp);
        return 0;
    }
}

extern (C) int dta_get_data(const(char) * fname, DatasetRaw * dset,
                            Gretlopt opt, Prn * prn)
{
    DtaTable * dtab = null;
    int namelen = 32;
    int nvar = 0, nobs = 0;
    FILE * fp;
    DatasetRaw * newset = null;
    GretlStringTable * st = null;
    Prn * vprn = prn;
    int err = 0;

    /* (re-)initialize file-scope globals */
    stata_globals_init();

    fp = gretl_fopen(fname, "rb");
    if (fp == null) {
        return GretlError.fopen;
    }

    if (opt & Gretlopt.OPT_Q) {
        vprn = null;
    }

    if (new_style_dta_file(fp)) {
        dtab = dta_table_new(&err);
        if (!err) {
            err = parse_dta_117_header(fp, dtab, prn, vprn);
        }
        if (!err) {
            nvar = dtab.nvar;
            nobs = dtab.nobs;
        }
    } else {
        err = parse_old_dta_header(fp, &namelen, &nvar, &nobs, prn, vprn);
    }

    if (err) {
        if (err != GretlError.alloc && err != E_NOTYET) {
            pputs(prn, _("This file does not seem to be a valid Stata data file").ptr);
            pputc(prn, '\n');
        }
        fclose(fp);
        return GretlError.data;
    }

    if (vprn != null) {
        pprintf(vprn, "endianness: %s\n", (stata_endian == Endian.big)?
                "big".ptr : "little".ptr);
        pprintf(vprn, "number of variables = %d\n", nvar);
        pprintf(vprn, "number of observations = %d\n", nobs);
        pprintf(vprn, "length of varnames = %d\n", namelen);
    }

    newset = datainfo_new();
    if (newset == null) {
        pputs(prn, _("Out of memory\n").ptr);
        fclose(fp);
        return GretlError.alloc;
    }

    newset.v = nvar + 1;
    newset.n = nobs;
    dataset_obs_info_default(newset);

    err = start_new_Z(newset, cast(Gretlopt) 0);
    if (err) {
        pputs(prn, _("Out of memory\n").ptr);
        free_datainfo(newset);
        fclose(fp);
        return GretlError.alloc;
    }

    if (stata_13) {
        if (dtab == null) {
            fprintf(stderr, "Got stata_13 but @dtab is NULL!\n");
            err = GretlError.data;
        } else {
            err = read_dta_117_data(fp, newset, &st, dtab, prn, vprn);
        }
    } else {
        err = read_old_dta_data(fp, newset, &st, namelen, prn, vprn);
    }

    if (err) {
        destroy_dataset(newset);
        if (st != null) {
            gretl_string_table_destroy(st);
        }
    } else {
        int merge = (dset.Z != null);

        if (fix_varname_duplicates(newset)) {
            pputs(prn, _("warning: some variable names were duplicated\n").ptr);
        }

        if (st != null) {
            gretl_string_table_finalize(st, newset);
            gretl_string_table_destroy(st);
        }

        err = merge_or_replace_data(dset, &newset, get_merge_opts(opt), prn);

        if (!err && !merge) {
            dataset_add_import_info(dset, fname, GretlFileType.dta);
        }
    }

    fclose(fp);
    free(dtab);

    return err;
}
