/* Port of plugin/excel_import.c: importer for Excel (.xls, BIFF) workbooks.

   In the C this is a gretl plugin (loaded by get_plugin_function("xls_get_data")); here it is an
   ordinary module that plugins.d returns for that name. The GTK sheet-selection dialog is dropped:
   xls_get_data() keeps its C signature (including the unused trailing GtkWidget * argument, now
   void *), and the OPT_G ("GUI") branch behaves like the command-line one: the sheet and the column
   and row offsets come from @sheetname and @list only. Originally based on the Gnumeric excel
   plugin by Michael Meeks.

   Differences from the C:
   - glib's g_malloc/g_strdup/g_free are malloc/g_strdup/free;
   - g_utf16_to_utf8() (used for codepage-1200 strings) is the local utf16le_to_utf8();
   - only little-endian hosts are supported by get_le_double() (the G_BIG_ENDIAN branch is dropped);
   - the C "#undef FORMAT_INFO" branches (COLINFO/XF/FORMAT diagnostics) are dropped;
   - BUILD_DATE (from the generated build.h) is the compile date of this module. */

module gretl.plugin.excel_import;


import core.stdc.ctype : isalnum, isalpha, isdigit, ispunct;
import core.stdc.math : isnan;
import core.stdc.stdarg;
import core.stdc.stdio : fflush, fprintf, fputs, printf, puts, sprintf, stderr, vfprintf;
import core.stdc.stdlib : atof, free, getenv, malloc, realloc;
import core.stdc.string : memcpy, memset, strchr, strcmp, strlen, strncat, strspn;

import gretl.plugin.biffh;
import gretl.csvdata : import_na_string, import_obs_label, reverse_data;
import gretl.dataio : check_varname, dataset_add_import_info, get_date_x, get_merge_opts, merge_or_replace_data;
import gretl.dataioh : GretlFileType;
import gretl.dataset : dataset_allocate_obs_markers, dataset_obs_info_default, datainfo_new, destroy_dataset;
import gretl.gretl_errors : gretl_errmsg_get, gretl_errmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : gretl_list_append_term, in_gretl_list;
import gretl.gretl_prn : Prn, gretl_print_get_buffer, pprintf, pputc, pputs;
import gretl.gretl_string_table : GretlStringTable, gretl_string_table_destroy, gretl_string_table_finalize,
    gretl_string_table_index, gretl_string_table_new, gretl_string_table_validate;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.plugin.import_common;
import gretl.plugin.importerh;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.plugin.ms_ole : ms_ole_destroy, ms_ole_open, ms_ole_stream_close, ms_ole_stream_open_workbook,
    ms_ole_stream_position;
import gretl.plugin.ms_oleh;
import core.stdc.config : c_ulong;
import gretl.gretl_utils : printlist;
import gretl.versionh : GRETL_VERSION;
import gretl.strutils : fix_varname_duplicates, gretl_charsub, gretl_utf8_strncat_trim, numeric_string, tailstrip;
import gretl.plugin.workbook : excel_book_get_info, ms_biff_query_destroy, ms_biff_query_new, ms_biff_query_next;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale, iso_to_ascii;

/* unported dependencies: per-file declarations from the 2026c headers */
private enum BUILD_DATE = __DATE__;

private struct SheetRow {
    int last, end;
    char ** cells;
}

private struct XlsInfo {
    int codepage;
    char ** sst;
    int sstsize;
    int sstnext;
    int datacols;
    int totcols;
    int nrows;
    SheetRow * rows;
    char * blank_col;
    int * codelist;
    GretlStringTable * st;
}

__gshared int debug_print;

private bool cell_record(int r)
{
    return (r == BiffRecord.label ||
            r == BiffRecord.string ||
            r == BiffRecord.number ||
            r == BiffRecord.rk ||
            r == BiffRecord.mulrk ||
            r == BiffRecord.formula ||
            r == BiffRecord.labelsst);
}

enum VarnameError {
    ok,
    null_,
    notstr,
    invalid,
    none
}

__gshared const(char) * adjust_rc = "Perhaps you need to adjust the starting column or row?";

private extern (C) int dbprintf(const(char) * format, ...)
{
    va_list args;
    int len = 0;

    if (debug_print) {
        va_start(args, format);
        len = vfprintf(stderr, format, args);
        va_end(args);
        fflush(stderr);
    }

    return len;
}

private void print_version()
{
    dbprintf("gretl, version %s, %s %s\n", GRETL_VERSION.ptr,
             _("build date").ptr, BUILD_DATE.ptr);
}

private double get_le_double(const(ubyte) * rec)
{
    double d;

    memcpy(&d, rec, 8);

    return d;
}

private double biff_get_rk(const(ubyte) * ptr)
{
    enum EType {
        ieee = 0,
        ieee_x100,
        int_,
        int_x100
    }
    int number;
    int type;

    number = cast(int) MS_OLE_GET_GUINT32(ptr);
    type = (number & 0x3);

    switch (type) {
    case EType.ieee:
    case EType.ieee_x100:
        {
            ubyte[8] tmp;
            double answer;
            int lp;

            for (lp = 0; lp < 4; lp++) {
                tmp[lp + 4] = (lp > 0) ? ptr[lp]: (ptr[lp] & 0xfc);
                tmp[lp] = 0;
            }
            answer = get_le_double(tmp.ptr);
            return (type == EType.ieee_x100)? answer / 100 : answer;
        }
    case EType.int_:
        return cast(double) (number >> 2);
    case EType.int_x100:
        number >>= 2;
        if ((number % 100) == 0) {
            return cast(double) (number/100);
        } else {
            return cast(double) (number/100.0);
        }
    default:
        break;
    }

    return NADBL;
}

private char * convert8to7(const(char) * s, int count)
{
    char * dest;
    int n;

    /* we'll skip any leading space */
    n = cast(int) strspn(s, " \t");
    count -= n;

    if (count <= 0) {
        dest = g_strdup("");
    } else {
        if (count > VNAMELEN - 1) {
            count = VNAMELEN - 1;
        }
        dest = cast(char *) malloc(VNAMELEN);
        *dest = '\0';
        s += n;
        strncat(dest, s, count);
        iso_to_ascii(dest);
        tailstrip(dest);
    }

    dbprintf("convert8to7: returning '%s'\n", dest);

    return dest;
}

private char * convert16to7(const(ubyte) * s, int count)
{
    char * dest;
    int i, u, j = 0;

    dest = cast(char *) malloc(VNAMELEN);
    if (dest == null) {
        return null;
    }

    memset(dest, 0, VNAMELEN);

    for (i=0; i<count && j<VNAMELEN-1; i++) {
        u = MS_OLE_GET_GUINT16(s);
        s += 2;
        if ((isalnum(u) || ispunct(u)) && u < 128) {
            dest[j++] = cast(char) u;
        }
    }

    dbprintf("convert16to7: returning '%s'\n", dest);

    return dest;
}

/* Replacement for g_utf16_to_utf8() on @len UTF-16LE code units at @s: returns a malloc'd
   string, or null on an unpaired surrogate (as with glib, the offending input is an error).
   *got is the number of units consumed and *wrote the number of bytes written. */

private char * utf16le_to_utf8(const(ubyte) * s, int len, int * got, int * wrote)
{
    char * ret = cast(char *) malloc(cast(size_t) len * 3 + 1);
    int i = 0, n = 0;

    while (i < len) {
        uint c = MS_OLE_GET_GUINT16(s + 2*i);

        i++;
        if (c >= 0xd800 && c < 0xdc00) {
            uint c2 = (i < len) ? MS_OLE_GET_GUINT16(s + 2*i) : 0;

            if (c2 >= 0xdc00 && c2 < 0xe000) {
                c = 0x10000 + ((c - 0xd800) << 10) + (c2 - 0xdc00);
                i++;
            } else {
                *got = i - 1;
                *wrote = n;
                free(ret);
                return null;
            }
        } else if (c >= 0xdc00 && c < 0xe000) {
            *got = i - 1;
            *wrote = n;
            free(ret);
            return null;
        }
        if (c < 0x80) {
            ret[n++] = cast(char) c;
        } else if (c < 0x800) {
            ret[n++] = cast(char) (0xc0 | (c >> 6));
            ret[n++] = cast(char) (0x80 | (c & 0x3f));
        } else if (c < 0x10000) {
            ret[n++] = cast(char) (0xe0 | (c >> 12));
            ret[n++] = cast(char) (0x80 | ((c >> 6) & 0x3f));
            ret[n++] = cast(char) (0x80 | (c & 0x3f));
        } else {
            ret[n++] = cast(char) (0xf0 | (c >> 18));
            ret[n++] = cast(char) (0x80 | ((c >> 12) & 0x3f));
            ret[n++] = cast(char) (0x80 | ((c >> 6) & 0x3f));
            ret[n++] = cast(char) (0x80 | (c & 0x3f));
        }
    }

    ret[n] = '\0';
    *got = i;
    *wrote = n;

    return ret;
}

private char * copy_unicode_string(XlsInfo * xi, ubyte * src, int remlen,
                                   int * skip, int * slop)
{
    int count = MS_OLE_GET_GUINT16(src);
    ubyte flags = *(src + 2);
    int this_skip = 3, skip_to_next = 3;
    int csize = (flags & 0x01)? 2 : 1;
    char * ret = null;

    dbprintf("copy_unicode_string: count = %d, csize = %d\n",
             count, csize);

    if (flags & 0x08) {
        dbprintf(" contains Rich-Text info\n");
    }
    if (flags & 0x04) {
        dbprintf(" contains Far-East info\n");
    }

    skip_to_next += count * csize;

    if (flags & 0x08) {
        ushort rich_text_info_len = 0;

        rich_text_info_len = cast(ushort) (4 * MS_OLE_GET_GUINT16(src + 3));
        this_skip += 2;
        skip_to_next += 2 + rich_text_info_len;
    }

    if (flags & 0x04) {
        uint far_east_info_len = 0;
        int far_east_offset = 3;

        if (flags & 0x08) {
            far_east_offset = 5;
        }
        far_east_info_len = MS_OLE_GET_GUINT32(src + far_east_offset);
        this_skip += 4;
        skip_to_next += 4 + far_east_info_len;
    }

    /* skip for the next read */
    if (skip != null) {
        *skip = skip_to_next;
    }

    /* size check */
    if (slop != null) {
        if (remlen > 0 && this_skip + count > remlen) {
            *slop = this_skip + count - remlen;
        } else {
            *slop = 0;
        }
    }

    if (count > 64) {
        /* let's not mess with excessive strings */
        ret = g_strdup("bigstr");
    } else if (csize == 1) {
        char[68] show = 0;

        strncat(show.ptr, cast(char *) src + this_skip, count);
        dbprintf("original string = '%s'\n", show.ptr);
        ret = convert8to7(cast(char *) src + this_skip, count);
    } else {
        if (xi.codepage == 1200) {
            int got, wrote;

            ret = utf16le_to_utf8(src + this_skip, count, &got, &wrote);
            dbprintf("utf16_to_utf8: got=%d, wrote=%d\n", got, wrote);
            if (ret == null) {
                fprintf(stderr, "Invalid sequence in conversion input\n");
            }
        }

        if (ret == null) {
            /* fallback */
            ret = convert16to7(src + this_skip, count);
        }
    }

    return ret;
}

private char * make_string(char * str)
{
    char * ret = null;

    if (str != null) {
        ret = g_strdup_printf("\"%s", str);
        free(str);
    } else {
        ret = g_strdup("\"");
    }

    return ret;
}

__gshared private int row_col_err_prevrow = -1, row_col_err_prevcol = -1;

private int row_col_err(int row, int col, Prn * prn)
{
    int err = 0;

    if (row < 0 || col < 0) {
        fprintf(stderr, "Error: got row=%d, col=%d\n", row, col);
        err = 1;
    } else if (row == row_col_err_prevrow && col == row_col_err_prevcol) {
        pprintf(prn, "Error: found a second cell entry for cell (%d, %d)\n",
                row_col_err_prevrow, row_col_err_prevcol);
        err = 1;
    }

    row_col_err_prevrow = row;
    row_col_err_prevcol = col;

    return err;
}

/* This function is called on LABELSST records: we check for possible
   NAs and also for numeric values that may have strayed into the
   string table, and that shouldn't really be treated as quoted
   strings.  We don't just use strtod at the outset, because some
   XLS files (put out by agencies that should know better!) contain
   numerical values that are "nicely formatted" as text using spaces
   or commas for thousands separators -- we try stripping this
   junk out first before doing the numeric-value test.
*/

__gshared private int check_copy_string_warned = 0;

private int check_copy_string(SheetRow * prow, int row, int col,
                              int idx, const(char) * s)
{
    dbprintf("inspecting sst[%d] = '%s'\n", idx, s);

    if (row > 0 && col > 0) {
        const(char) * numok = "0123456789 -,.";
        int i, len = cast(int) strlen(s);
        int commas = 0, digits = 0;

        if (len == 0) {
            dbprintf(" converting to NA\n");
            prow.cells[col] = g_strdup("-999");
            return 0;
        }

        for (i=0; i<len; i++) {
            if (strchr(numok, s[i]) == null) {
                /* does not look promising for numerical value */
                len = 0;
                break;
            }
            if (isdigit(s[i])) {
                digits++;
            } else if (s[i] == ',') {
                commas++;
            }
        }

        if (len > 0 && digits > 0) {
            /* may be numerical? */
            char * p, q = cast(char *) malloc(len + 1);

            if (q == null) return 1;

            p = q;
            for (i=0; i<len; i++) {
                if (s[i] != ' ' && s[i] != ',') {
                    *p++ = s[i];
                }
                if (commas == 1 && s[i] == ',') {
                    /* single comma could be for 000s, or decimal */
                    if (!check_copy_string_warned) {
                        fprintf(stderr, "Warning: found ambiguous comma in '%s'\n", s);
                        check_copy_string_warned = 1;
                    }
                    if (len - i != 4) {
                        /* comma is probably decimal separator? */
                        *p++ = '.';
                    }
                }
            }
            *p = '\0';

            /* If we don't do a rigorous check on q, as below, we
               may end up with zeros where there should be NAs.
            */
            if (numeric_string(q)) {
                dbprintf(" taking '%s' to be numeric string: %s\n", s, q);
                prow.cells[col] = q;
                return 0;
            } else {
                free(q);
            }
        }
    }

    dbprintf(" copying '%s' into place as string\n", s);
    prow.cells[col] = g_strdup_printf("\"%s", s);

    return 0;
}

private int is_date_format(int fmt)
{
    int ret = 0;

    fprintf(stderr, "is_date_format? fmt=%d\n", fmt);

    if (fmt >= 14 && fmt <= 22) {
        ret = 1;
    } else if (fmt >= 45 && fmt <= 47) {
        ret = 1;
    } else if (fmt >= 50 && fmt <= 58) {
        ret = 1;
    } else if (fmt == 164) {
        /* user-defined: FRED uses this */
        ret = 1;
    }

    return ret;
}

private int wbook_find_format(WBook * book, int xfref)
{
    int fmt = -1;

    if (book.xf_list != null && xfref < book.xf_list[0]) {
        fmt = book.xf_list[xfref + 1];
    }

    return fmt;
}

private int func_is_date(ubyte * data, int version_)
{
    /* check for built-in DATE function */
    if (version_ < MsBiffVersion.v4) {
        return MS_OLE_GET_GUINT8(data) == 65;
    } else {
        return MS_OLE_GET_GUINT16(data) == 65;
    }
}

private int t_func_size(int v) { return (v < MsBiffVersion.v4)? 2 : 3; }
private int t_ref_size(int v)  { return (v < MsBiffVersion.v8)? 4 : 5; }

/* Could be a date formula?  If so, it should have 3 cell reference fields
   and a trailing function ID == 65 */

private void check_for_date_formula(BiffQuery * q, WBook * book)
{
    int version_ = book.version_;
    int offset = (version_ < MsBiffVersion.v5)? 16 : 20;
    ubyte * fdata = q.data + offset;
    ushort sz, targ;
    ubyte u1;
    int i;

    targ = cast(ushort) (3 * t_ref_size(version_) + t_func_size(version_));

    if (version_ < MsBiffVersion.v3) {
        sz = MS_OLE_GET_GUINT8(fdata);
        fdata += 1;
    } else {
        sz = MS_OLE_GET_GUINT16(fdata);
        fdata += 2;
    }

    /* There's a one-byte ambiguity over the size of the
       function ID field in the OpenOffice.org doc for
       BIFF, so we'll allow sz to be one byte bigger than
       targ.
    */
    if (sz != targ && sz != targ + 1) {
        return;
    }

    for (i=0; i<3; i++) {
        /* token ID */
        u1 = MS_OLE_GET_GUINT8(fdata);
        if (u1 != 0x44) { /* 0x44 = tRef */
            return;
        }
        fdata += t_ref_size(version_);
    }

    u1 = MS_OLE_GET_GUINT8(fdata);

    if (u1 == 0x41 && func_is_date(fdata + 1, version_)) { /* 0x41 = tFunc */
        fprintf(stderr, "Got DATE formula in first column\n");
        book_set_numeric_dates(book);
    }
}

/* Excel's NA() formula stores a "result" of 0.0, so if
   we get a formula result of zero we need to check for
   NA() and react accordingly. Use of NA() can be found
   in XLS files downloaded from FRED.
*/

private int is_na_formula(ubyte * ptr, WBook * book)
{
    int version_ = book.version_;
    int offset = (version_ < MsBiffVersion.v5)? 16 : 20;
    ushort sz;

    ptr += offset;
    sz = MS_OLE_GET_GUINT16(ptr);

    /* 0x42: indicates a built-in function
       10 is the index of the BIFF NA() function
    */
    if (sz == 4 && ptr[2] == 0x42 && ptr[4] == 10) {
        return 1;
    } else {
        return 0;
    }
}

__gshared private char ** process_item_string_targ;
__gshared private int process_item_slop; /* SST overslop */

private int process_item(BiffQuery * q, WBook * book, XlsInfo * xi,
                         Prn * prn)
{
    SheetRow * prow = null;
    ubyte * ptr = null;
    int i = 0, j = 0;
    double val;

    if (cell_record(q.ls_op)) {
        i = ex_getrow(q);
        j = ex_getcol(q);
        if (row_col_err(i, j, prn)) {
            return 1;
        }
        if (q.ls_op == BiffRecord.number || q.ls_op == BiffRecord.rk || q.ls_op == BiffRecord.mulrk) {
            ushort xfref = cast(ushort) ex_getxf(q);
            int fmt = wbook_find_format(book, xfref);

            if (i == book.row_offset + 1 &&
                j == book.col_offset &&
                is_date_format(fmt)) {
                fprintf(stderr, "Testing first obs cell (%d, %d): date format %d\n",
                        i, j, fmt);
                book_set_numeric_dates(book);
            }
        }
    }

    switch (cast(int) q.ls_op) {

    case BiffRecord.sst: {
        int k, skip, remlen, oldsz = xi.sstsize;
        ushort sz;

        if (xi.sst != null) {
            fprintf(stderr, "Got a second string table: this is nonsense\n");
            return 1;
        }

        sz = MS_OLE_GET_GUINT16(q.data + 4);
        xi.sstsize += sz;
        xi.sst = cast(char **) realloc(xi.sst, xi.sstsize * (*xi.sst).sizeof);
        if (xi.sst == null) {
            return 1;
        }

        dbprintf("Got SST: allocated for %d strings (%d bytes), %p\n",
                 xi.sstsize, cast(int) (xi.sstsize * (*xi.sst).sizeof), cast(void *) xi.sst);

        for (k=oldsz; k<xi.sstsize; k++) {
            /* careful: initialize all pointers to NULL */
            xi.sst[k] = null;
        }

        ptr = q.data + 8;

        for (k=oldsz; k<xi.sstsize; k++) {
            remlen = cast(int) (q.length - (ptr - q.data));
            dbprintf("Working on sst[%d], data offset=%d, remlen=%d\n",
                     k, cast(int) (ptr - q.data), remlen);
            if (remlen <= 0) {
                break;
            }
            xi.sst[k] = copy_unicode_string(xi, ptr, remlen, &skip, &process_item_slop);
            ptr += skip;
        }

        if (k < xi.sstsize) {
            xi.sstnext = k;
        }

        break;
    }

    case BiffRecord.continue_:
        dbprintf("Got CONTINUE, xi->sstnext = %d, len = %d\n",
                 xi.sstnext, cast(int) q.length);
        if (xi.sstnext > 0) {
            int k, skip, remlen;

            ptr = q.data;
            if (process_item_slop > 0) {
                ubyte flags = *ptr;
                int csize = (flags & 0x01)? 2 : 1;

                dbprintf("BIFF_CONTINUE: slop = %d, csize = %d\n", process_item_slop,
                         csize);
                ptr += 1 + csize * process_item_slop;
            }
            for (k=xi.sstnext; k<xi.sstsize; k++) {
                remlen = cast(int) (q.length - (ptr - q.data));
                if (remlen <= 0) {
                    break;
                }
                dbprintf("Working on sst[%d], remlen = %d\n", k, remlen);
                if (xi.sst[k] != null) {
                    free(xi.sst[k]);
                }
                xi.sst[k] = copy_unicode_string(xi, ptr, remlen, &skip, &process_item_slop);
                ptr += skip;
            }
            if (k < xi.sstsize) {
                xi.sstnext = k;
            }
        }
        break;

    case BiffRecord.label:
        dbprintf("Got LABEL, row=%d, col=%d\n", i, j);
        if (allocate_row_col(i, j, book, xi)) {
            return GretlError.alloc;
        } else {
            uint len = MS_OLE_GET_GUINT16(q.data + 6);

            prow = xi.rows + i;
            ptr = q.data + 8;
            fprintf(stderr, "BIFF_LABEL: calling convert8to7\n");
            prow.cells[j] = make_string(convert8to7(cast(char *) ptr, len));
        }
        break;

    case BiffRecord.labelsst:
        dbprintf("Got LABELSST, row=%d, col=%d\n", i, j);
        if (allocate_row_col(i, j, book, xi)) {
            return GretlError.alloc;
        } else {
            uint sidx = MS_OLE_GET_GUINT16(q.data + 6);

            prow = xi.rows + i;
            if (sidx >= xi.sstsize) {
                pprintf(prn, _("String index too large").ptr);
                pputc(prn, '\n');
            } else if (xi.sst[sidx] != null) {
                check_copy_string(prow, i, j, sidx, xi.sst[sidx]);
            } else {
                dbprintf("sst[%d] seems to be NULL, leaving string blank\n", cast(int) sidx);
                prow.cells[j] = cast(char *) malloc(2);
                if (prow.cells[j] != null) {
                    prow.cells[j][0] = '\0';
                }
            }
        }
        break;

    case BiffRecord.number:
        if (allocate_row_col(i, j, book, xi)) {
            return GretlError.alloc;
        } else {
            val = get_le_double(q.data + 6);
            prow = xi.rows + i;
            prow.cells[j] = g_strdup_printf("%.15g", val);
            dbprintf("Got NUMBER (%g), row=%d, col=%d\n", val, i, j);
        }
        break;

    case BiffRecord.rk:
        if (allocate_row_col(i, j, book, xi)) {
            return GretlError.alloc;
        } else {
            val = biff_get_rk(q.data + 6);
            prow = xi.rows + i;
            prow.cells[j] = g_strdup_printf("%.15g", val);
            dbprintf("Got RK (%g), row=%d, col=%d\n", val, i, j);
        }
        break;

    case BiffRecord.mulrk: {
        int k, ncols = (q.length - 6) / 6;

        dbprintf("Got MULRK, row=%d, first_col=%d, ncols=%d\n", i, j, ncols);
        for (k=0; k<ncols; k++) {
            if (allocate_row_col(i, j, book, xi)) {
                return GretlError.alloc;
            }
            val = biff_get_rk(q.data + 6 + 6 * k);
            prow = xi.rows + i; /* might have moved */
            prow.cells[j] = g_strdup_printf("%.15g", val);
            dbprintf(" MULRK[col=%d] = %g\n", j, val);
            j++;
        }
        break;
    }

    case BiffRecord.formula:
        dbprintf("Got FORMULA, row=%d, col=%d\n", i, j);
        if (allocate_row_col(i, j, book, xi)) {
            return GretlError.alloc;
        } else {
            /* the result of the formula should be at offset 6 */
            ptr = q.data + 6;
            prow = xi.rows + i;
            if (ptr[6] == 0xff && ptr[7] == 0xff) {
                /* string, boolean or error result */
                ubyte fcode = ptr[0];

                if (fcode == 0x0) {
                    /* string formula: record the target for following
                       STRING record */
                    process_item_string_targ = prow.cells + j;
                } else if (fcode == 0x1) {
                    /* boolean value */
                    prow.cells[j] = g_strdup((ptr[2])? "1" : "0");
                } else if (fcode == 0x2 || fcode == 0x3) {
                    /* error code or empty */
                    prow.cells[j] = g_strdup("-999");
                } else {
                    fprintf(stderr, "Bad formula code 0x%u\n",
                            cast(uint) fcode);
                    prow.cells[j] = g_strdup("-999");
                }
            } else {
                /* should have a floating-point result */
                val = get_le_double(ptr);
                if (val == 0.0 && is_na_formula(q.data, book)) {
                    dbprintf("floating-point value = na()\n");
                    prow.cells[j] = g_strdup("-999");
                } else if (isnan(val)) {
                    dbprintf("floating-point value is NaN\n");
                    prow.cells[j] = g_strdup("-999");
                } else {
                    dbprintf(" floating-point value = %g\n", val);
                    prow.cells[j] = g_strdup_printf("%.15g", val);
                    if (i == book.row_offset + 1 && j == book.col_offset) {
                        /* could be a date formula? */
                        check_for_date_formula(q, book);
                    }
                }
            }
        }
        break;

    case BiffRecord.string:
        if (process_item_string_targ == null) {
            dbprintf("String record without preceding string formula\n");
        } else {
            char * tmp = copy_unicode_string(xi, q.data, 0, null, null);

            *process_item_string_targ = make_string(tmp);
            dbprintf("Filled out string formula with '%s'\n", *process_item_string_targ);
            process_item_string_targ = null; /* handled */
        }
        break;

    case BiffRecord.bof:
        if (xi.rows != null) {
            fprintf(stderr, "BOF when current sheet is not flushed\n");
            return 1;
        }
        {
            uint version_, boftype;

            version_ = MS_OLE_GET_GUINT16(q.data + 0);
            boftype = MS_OLE_GET_GUINT16(q.data + 2);
            dbprintf("Got BOF: version=%x, type=%x\n", version_, boftype);
        }
        break;

    case BiffRecord.format: {
        int idx = MS_OLE_GET_GUINT16(q.data + 0);

        if (idx >= 14 && idx <= 17) {
            fprintf(stderr, "Got date format: index %d\n", idx);
        }
        break;
    }

    default:
        break;
    }

    return 0;
}

private int handled_record(BiffQuery * q)
{
    if (q.opcode == BiffRecord.sst ||
        q.opcode == BiffRecord.continue_ ||
        q.opcode == BiffRecord.labelsst ||
        q.opcode == BiffRecord.mulrk ||
        q.opcode == BiffRecord.formula) {
        return 1;
    }

    if (q.ms_op == 0x02) {
        if (q.ls_op == BiffRecord.label ||
            q.ls_op == BiffRecord.number ||
            q.ls_op == BiffRecord.rk ||
            q.ls_op == BiffRecord.string) {
            return 1;
        }
    }

    if (q.ms_op == 0x08) {
        if (q.ls_op == BiffRecord.bof) return 1;
    }

    return 0;
}

__gshared private int process_sheet_skipped;

private int process_sheet(const(char) * filename, WBook * book, XlsInfo * xi,
                          Prn * prn)
{
    int err = 0, gotbof = 0, eofcount = 0;
    long offset = book.byte_offsets[book.selected];
    MsOleStream * stream;
    MsOleErr result;
    BiffQuery * q;
    MsOle * file;

    if (ms_ole_open(&file, filename)) {
        return 1;
    }

    result = ms_ole_stream_open_workbook(&stream, file);

    if (result != MsOleErr.ok) {
        ms_ole_destroy(&file);
        return 1;
    }

    fputs("Reading file...\n", stderr);
    q = ms_biff_query_new(stream);

    while (!gotbof && ms_biff_query_next(q)) {
        if (q.ls_op == BiffRecord.bof) {
            gotbof = 1;
            break;
        }
    }

    if (!gotbof) {
        pprintf(prn, _("%s: No BOF record found").ptr, filename);
        return 1;
    }

    while (!err && ms_biff_query_next(q)) {
        dbprintf("At %lu: q->opcode=0x%02x\n", cast(c_ulong) q.streamPos, cast(uint) q.opcode);
        if (q.opcode == BiffRecord.eof) {
            dbprintf("got MSEOF at %lu\n", cast(c_ulong) ms_ole_stream_position(stream));
            eofcount++;

            if (eofcount == 1) {
                if (ms_ole_stream_position(stream) < offset) {
                    /* skip to the worksheet we want? */
                    while (q.streamPos < offset && ms_biff_query_next(q)) {
                        /* keep reading */
                    }
                    fprintf(stderr, "skipped forward to %lu\n",
                            cast(c_ulong) q.streamPos);
                } else {
                    fprintf(stderr, "reading worksheet at %lu\n",
                            cast(c_ulong) ms_ole_stream_position(stream));
                }
            }

            if (eofcount == 2) {
                break;
            } else {
                continue;
            }
        }

        if (handled_record(q)) {
            err = process_item(q, book, xi, prn);
        } else if (q.ms_op == 0x02 && q.ls_op == BiffRecord.row) {
            dbprintf("Got BIFF_ROW\n");
        } else if (q.opcode == BiffRecord.dbcell) {
            dbprintf("Got BIFF_DBCELL\n");
        } else if (q.opcode == 0x42) {
            if (q.length == 2 && q.data != null) {
                int cp = MS_OLE_GET_GUINT16(q.data);

                fprintf(stderr, "CODEPAGE: got %d\n", cp);
                xi.codepage = cp;
            }
        } else {
            if (q.opcode != process_sheet_skipped) {
                dbprintf("skipping unhandled opcode 0x%02x\n", cast(uint) q.opcode);
            }
            process_sheet_skipped = q.opcode;
        }
    }

    ms_biff_query_destroy(q);
    ms_ole_stream_close(&stream);
    ms_ole_destroy(&file);

    return err;
}

private void row_init(SheetRow * row)
{
    row.last = 0;
    row.end = 0;
    row.cells = null;
}

__gshared private int allocate_row_col_started;

private int allocate_row_col(int i, int j, WBook * book,
                             XlsInfo * xi)
{
    int k;

    if (!allocate_row_col_started && i > book.row_offset) {
        book.row_offset = i;
        fprintf(stderr, "Missing rows: trying an offset of %d\n", i);
    }

    allocate_row_col_started = 1;

    dbprintf("allocate: row=%d, col=%d, nrows=%d\n", i, j, xi.nrows);

    if (i >= xi.nrows) {
        int new_nrows = (i / 16 + 1) * 16;
        SheetRow * myrows;

        myrows = cast(SheetRow *) realloc(xi.rows, new_nrows * (*myrows).sizeof);
        if (myrows == null) {
            return 1;
        }

        xi.rows = myrows;

        for (k=xi.nrows; k<new_nrows; k++) {
            dbprintf("allocate: initing rows[%d]\n", k);
            row_init(&xi.rows[k]);
            dbprintf("rows[%d].end=%d\n", i, xi.rows[k].end);
        }
        xi.nrows = new_nrows;
    }

    dbprintf("allocate: col=%d and rows[%d].end = %d\n", j, i, xi.rows[i].end);

    if (j >= xi.rows[i].end) {
        int newcol = (j / 16 + 1) * 16;
        char ** cells;

        dbprintf("allocate: reallocing rows[%d].cells to size %d\n", i, newcol);
        cells = cast(char **) realloc(xi.rows[i].cells, newcol * (*cells).sizeof);

        if (cells == null) {
            return 1;
        }

        xi.rows[i].cells = cells;

        for (k=xi.rows[i].end; k<newcol; k++) {
            xi.rows[i].cells[k] = null;
        }
        xi.rows[i].end = newcol;
    }

    if (j > xi.rows[i].last) {
        xi.rows[i].last = j;
    }

    return 0;
}

private void xls_info_init(XlsInfo * xi)
{
    xi.codepage = 0;
    xi.sst = null;
    xi.sstsize = 0;
    xi.datacols = 0;
    xi.totcols = 0;
    xi.nrows = 0;
    xi.rows = null;
    xi.blank_col = null;
    xi.codelist = null;
    xi.st = null;
}

private void free_xls_info(XlsInfo * xi)
{
    int i, j;

    dbprintf("free_xls_info(), nrows=%d\n", xi.nrows);

    /* free shared string table */
    if (xi.sst != null) {
        for (i=0; i<xi.sstsize; i++) {
            free(xi.sst[i]);
        }
        free(xi.sst);
        xi.sst = null;
    }

    /* free cells */
    if (xi.rows != null) {
        for (i=0; i<xi.nrows; i++) {
            if (xi.rows[i].cells == null) {
                dbprintf("rows[%d].cells = NULL, skipping free\n", i);
                continue;
            }
            for (j=0; j<xi.rows[i].end; j++) {
                if (xi.rows[i].cells[j] != null) {
                    dbprintf("Freeing rows[%d].cells[%d] at %p\n",
                             i, j, cast(void *) xi.rows[i].cells[j]);
                    free(xi.rows[i].cells[j]);
                }
            }
            dbprintf("Freeing rows[%d].cells at %p\n", i, cast(void *) xi.rows[i].cells);
            free(xi.rows[i].cells);
        }
        free(xi.rows);
    }

    free(xi.blank_col);
    free(xi.codelist);
    if (xi.st != null) {
        gretl_string_table_destroy(xi.st);
    }
}

private bool IS_STRING(const(char) * v) { return v[0] == '"'; }

/* check for full set of strings in first column to be read (which may
   be at an offset into the worksheet)
 */

private int first_col_strings(WBook * book, XlsInfo * xi)
{
    int i, j = book.col_offset;
    int startrow = book.row_offset + 1;
    int ret = 1;

    dbprintf("checking for first column strings...\n");

    for (i=startrow; i<xi.nrows; i++) {
        dbprintf("book->row_offset=%d, i=%d\n", book.row_offset, i);
        dbprintf("rows = %p\n", cast(void *) xi.rows);
        if (xi.rows == null || xi.rows[i].cells == null ||
            xi.rows[i].cells[j] == null ||
            !IS_STRING(xi.rows[i].cells[j])) {
            dbprintf("no: not a string at row %d\n", i);
            ret = 0;
            break;
        }
        dbprintf("first_col_strings: rows[%d].cells[%d]: '%s'\n", i, j,
                 xi.rows[i].cells[j]);
    }

    if (ret) {
        book_set_obs_labels(book);
    }

    return ret;
}

private bool obs_string(const(char) * s)
{
    return !strcmp(s, "obs") || !strcmp(s, "id");
}

private int check_all_varnames(WBook * book, XlsInfo * xi, Prn * prn)
{
    int j, i = book.row_offset;
    int startcol = book.col_offset;
    int realcols = 0;
    int vnames = 0;
    int ret = VarnameError.none;

    if (book_obs_labels(book) || book_numeric_dates(book)) {
        startcol++;
    }

    if (xi.rows[i].cells == null) {
        fprintf(stderr, "Row %d is empty, trying lower...\n", i);
        while (i < xi.nrows - 1 && xi.rows[i].cells == null) {
            book.row_offset += 1;
            i++;
        }
    }

    for (j=startcol; j<xi.totcols; j++) {
        if (xi.blank_col[j]) {
            continue;
        }
        if (xi.rows[i].cells[j] == null) {
            dbprintf("got_varnames: rows[%d].cells[%d] is NULL\n", i, j);
            break;
        }

        dbprintf("got_varnames: rows[%d].cells[%d] is '%s'\n", i, j,
                 xi.rows[i].cells[j]);

        if (IS_STRING(xi.rows[i].cells[j])) {
            /* skip beyond the quote */
            char * test = xi.rows[i].cells[j] + 1;

            /* "obs" or "id" is OK in the first col of the selection,
               but not thereafter */
            if (j == startcol && obs_string(test)) {
                /* pass along */
            } else {
                int verr = check_imported_varname(test, 0, i, j, prn);

                if (verr) {
                    return verr;
                }
            }
            vnames++;
        }
        realcols++;
    }

    if (vnames == realcols) {
        ret = VarnameError.ok;
    } else if (vnames > 0) {
        ret = VarnameError.notstr;
    }

    return ret;
}

private int missval_string(const(char) * s)
{
    s++;

    return (*s == '\0' || import_na_string(s));
}

private struct StringErr {
    int row;
    int column;
    char * str;
}

private void clear_string_err(StringErr * strerr)
{
    strerr.row = 0;
    strerr.column = 0;
    free(strerr.str);
    strerr.str = null;
}

private char * xls_cell(XlsInfo * x, int i, int j)
{
    return x.rows[i].cells[j];
}

/* check for invalid data in the selected data block */

private int check_data_block(WBook * book, XlsInfo * xi, int * missvals,
                             StringErr * strerr)
{
    int * codelist = null;
    int startcol = book.col_offset;
    int startrow = book.row_offset + 1;
    int j, i, err = 0;

    if (book_obs_labels(book) || book_numeric_dates(book)) {
        startcol++;
    }

    strerr.row = 0;
    strerr.column = 0;
    strerr.str = null;

    for (j=startcol; j<xi.totcols && !err; j++) {
        int strvals = 0;

        dbprintf("data_block: col=%d\n", j);
        if (xi.blank_col[j]) {
            continue;
        }
        for (i=startrow; i<xi.nrows; i++) {
            dbprintf(" rows[%d], end = %d\n", i, xi.rows[i].end);
            if (xi.rows[i].cells  == null) {
                dbprintf("  rows[%d].cells = NULL\n", i);
                *missvals = 1;
            } else if (j >= xi.rows[i].end) {
                dbprintf("  short row, fell off the end\n");
                *missvals = 1;
            } else if (xls_cell(xi, i, j) == null) {
                dbprintf("  rows[%d].cells[%d] = NULL\n", i, j);
                xi.rows[i].cells[j] = g_strdup("-999");
                *missvals = 1;
            } else if (IS_STRING(xls_cell(xi, i, j))) {
                if (missval_string(xls_cell(xi, i, j))) {
                    dbprintf("  rows[%d].cells[%d] = missval\n", i, j);
                    free(xi.rows[i].cells[j]);
                    xi.rows[i].cells[j] = g_strdup("-999");
                    *missvals = 1;
                } else {
                    dbprintf("  rows[%d].cells[%d]: %s (string)\n",
                             i, j, xls_cell(xi, i, j));
                    strvals++;
                    if (strerr.row == 0) {
                        strerr.row = i + 1;
                        strerr.column = j + 1;
                        strerr.str = g_strdup(xls_cell(xi, i, j));
                    }
                }
            } else {
                dbprintf("  rows[%d].cells[%d]: %s (numeric?)\n",
                         i, j, xls_cell(xi, i, j));
            }
        }
        if (strvals > 0) {
            dbprintf(" col %d: %d string values\n", j, strvals);
            if (strvals == xi.nrows - startrow) {
                int k = j - startcol + 1;

                fprintf(stderr, "col %d: all strings -> accept\n", j);
                codelist = gretl_list_append_term(&codelist, k);
                clear_string_err(strerr);
            } else {
                err = GretlError.data;
            }
        }
    }

    if (codelist != null) {
        printlist(codelist, "codelist");
        if (err) {
            free(codelist);
        } else {
            xi.codelist = codelist;
        }
    }

    return err;
}

/* determine the number of actual data columns, starting from a
   given offset into the worksheet, and allowing for the
   possibility that the first selected column contains string
   labels
*/

private int n_vars_from_col(WBook * book, int totcols, char * blank_col)
{
    int offset = book.col_offset;
    int i, nv = 1;

    if (book_time_series(book) || book_obs_labels(book)) {
        /* got a first non-data column */
        offset++;
    }

    for (i=offset; i<totcols; i++) {
        if (!blank_col[i]) nv++;
    }

    dbprintf("n_vars_from_col: totcols=%d, nv=%d\n", totcols, nv);

    return nv;
}

private int transcribe_data(WBook * book, XlsInfo * xi, DatasetRaw * dset,
                            Prn * prn)
{
    int startcol = book.col_offset;
    int roff = book.row_offset;
    int i, t, j = 1;
    int err = 0;

    if (book_obs_labels(book) || book_time_series(book)) {
        startcol++;
    }

    if (xi.codelist != null) {
        xi.st = gretl_string_table_new(xi.codelist);
        if (xi.st == null) {
            return GretlError.alloc;
        }
    }

    for (i=startcol; i<xi.totcols && !err; i++) {
        const(char) * val = null;
        int ts, strvals = 0;

        if (xi.blank_col[i]) {
            continue;
        }

        if (j >= dset.v) {
            break;
        }

        dset.varname[j][0] = 0;
        if (book_auto_varnames(book)) {
            sprintf(dset.varname[j], "v%d", j);
        } else if (xi.rows[roff].cells[i] == null) {
            sprintf(dset.varname[j], "v%d", j);
        } else if (i >= xi.rows[roff].end) {
            sprintf(dset.varname[j], "v%d", j);
        } else {
            strncat(dset.varname[j], xi.rows[roff].cells[i] + 1,
                    VNAMELEN - 1);
            dbprintf("accessing rows[%d].cells[%d] at %p\n",
                     roff, i, cast(void *) xi.rows[roff].cells[i]);
        }

        /* remedial: replace space with underscore */
        gretl_charsub(dset.varname[j], ' ', '_');

        err = check_varname(dset.varname[j]);
        if (err) {
            pprintf(prn, "%s\n", gretl_errmsg_get());
            break;
        }

        dbprintf("set varname[%d] = '%s'\n", j, dset.varname[j]);

        if (in_gretl_list(xi.codelist, j)) {
            strvals = 1;
        }

        for (t=0; t<dset.n && !err; t++) {
            ts = t + 1 + roff;
            if (xi.rows[ts].cells == null || i >= xi.rows[ts].end ||
                xi.rows[ts].cells[i] == null) {
                continue;
            }

            val = xi.rows[ts].cells[i];
            if (val != null && *val == '"') {
                val++;
            }

            dbprintf("accessing rows[%d].cells[%d] at %p\n", ts, i,
                     cast(void *) val);
            dbprintf("setting Z[%d][%d] = rows[%d].cells[%d] "
                     ~ "= '%s'\n", j, t, i, ts, val);

            if (strvals) {
                int xjt = gretl_string_table_index(xi.st, val, j, 0, prn);

                if (xjt > 0) {
                    dset.Z[j][t] = xjt;
                } else {
                    err = GretlError.data;
                }
            } else {
                dset.Z[j][t] = atof(val);
                if (dset.Z[j][t] == -999 || dset.Z[j][t] == -9999) {
                    dset.Z[j][t] = NADBL;
                }
            }
        }

        j++;
    }

    return err;
}

private int get_sheet_dimensions(WBook * book, XlsInfo * xi, Prn * prn)
{
    char * blanks = null;
    int i, j;

    /* trim any trailing blank rows */
    for (i=xi.nrows-1; i>=0; i--) {
        if (xi.rows[i].cells == null) {
            xi.nrows -= 1;
        } else {
            break;
        }
    }

    for (i=0; i<xi.nrows; i++) {
        if (xi.rows[i].cells != null) {
            if (xi.rows[i].last + 1 > xi.totcols) {
                xi.totcols = xi.rows[i].last + 1;
            }
        }
    }

    if (xi.totcols <= 0 || xi.nrows < 1) {
        pputs(prn, _("No data found.\n").ptr);
        pputs(prn, adjust_rc);
        return 1;
    }

    blanks = cast(char *) malloc(xi.totcols);
    if (blanks == null) {
        return GretlError.alloc;
    }

    memset(blanks, 1, xi.totcols);

    for (i=0; i<xi.nrows; i++) {
        if (xi.rows[i].cells == null) {
            continue;
        }
        for (j=0; j<=xi.rows[i].last; j++) {
            if (xi.rows[i].cells[j] != null) {
                if (blanks[j]) {
                    blanks[j] = 0;
                }
            }
        }
    }

    for (i=0; i<xi.totcols; i++) {
        if (!blanks[i]) {
            xi.datacols += 1;
        }
    }

    if (book_numeric_dates(book)) {
        xi.datacols -= 1;
    }

    printf("rows=%d, total cols=%d, data cols=%d\n", xi.nrows,
           xi.totcols, xi.datacols);

    if (xi.datacols < 1) {
        pputs(prn, _("No data found.\n").ptr);
        pputs(prn, adjust_rc);
        free(blanks);
        return 1;
    }

    xi.blank_col = blanks;

    return 0;
}

private int col0_is_numeric(XlsInfo * xi, int row_offset, int col_offset)
{
    int t, tstart = 1 + row_offset;
    int nx = 0;
    char * test;

    fprintf(stderr, "testing for all numerical values in col %d\n",
            col_offset);

    for (t=tstart; t<xi.nrows; t++) {
        test = xi.rows[t].cells[col_offset];
        if (!numeric_string(test)) {
            fprintf(stderr, " no: non-numeric cell at row %d\n", t + 1);
            return 0;
        } else if (test != null && *test != '\0') {
            nx++;
        }
    }

    return nx > 0;
}

private int alpha_cell(const(char) * s)
{
    if (s != null) {
        if (*s == '"' || *s == '\'') s++;
        return isalpha(*s);
    }

    return 0;
}

private void book_time_series_setup(WBook * book, DatasetRaw * newinfo, int pd)
{
    newinfo.pd = pd;
    newinfo.structure = TsCodes.time_series;

    fprintf(stderr, "stobs='%s'\n", newinfo.stobs.ptr);
    newinfo.sd0 = get_date_x(newinfo.pd, newinfo.stobs.ptr);
    fprintf(stderr, "sd0=%g\n", newinfo.sd0);

    book_set_time_series(book);
    book_unset_obs_labels(book);
}

/* Make a contiguous array of observation labels for the
   purpose of checking for dated observations. All we need
   here is a "shell"; the actual strings are already
   allocated.
*/

private char ** labels_array(XlsInfo * xi, int row_offset, int j,
                             DatasetRaw * newset)
{
    char * s;
    char ** labels = null;
    int i, t, ok = 1;

    for (t=0; t<newset.n; t++) {
        i = t + row_offset;
        s = xi.rows[i].cells != null ? xi.rows[i].cells[j] : null;
        if (s == null || *s == '\0') {
            ok = 0;
            break;
        }
    }

    if (ok) {
        labels = cast(char **) malloc(newset.n * (*labels).sizeof);
        if (labels != null) {
            for (t=0; t<newset.n; t++) {
                i = t + row_offset;
                labels[t] = xi.rows[i].cells[j];
            }
        }
    }

    return labels;
}

private void maybe_revise_xls_codelist(XlsInfo * xi)
{
    if (xi.codelist != null) {
        int i;

        for (i=1; i<=xi.codelist[0]; i++) {
            xi.codelist[i] += 1;
        }
    }
}

extern (C) int xls_get_data(const(char) * fname, int * list, char * sheetname,
                            DatasetRaw * dset, Gretlopt opt, Prn * prn,
                            void * parent)
{
    int gui = (opt & Gretlopt.OPT_G) != 0;
    WBook xbook;
    WBook * book = &xbook;
    XlsInfo xlsi;
    XlsInfo * xi = &xlsi;
    DatasetRaw * newset;
    StringErr strerr;
    int ts_markers = 0;
    int missvals = 0;
    char ** ts_S = null;
    int r0, c0;
    int merge = (dset.Z != null);
    int i, t, pd = 0;
    int err = 0;

    newset = datainfo_new();
    if (newset == null) {
        pputs(prn, _("Out of memory\n").ptr);
        return 1;
    }

    if (sheetname != null) {
        fprintf(stderr, "xls_get_data: sheetname='%s'\n", sheetname);
    }

    gretl_push_c_numeric_locale();

    wbook_init(book, list, sheetname);
    xls_info_init(xi);

    if (excel_book_get_info(fname, book)) {
        pputs(prn, _("Failed to get workbook info").ptr);
        err = 1;
    } else if (book.nsheets == 0) {
        pputs(prn, _("No worksheets found").ptr);
        err = 1;
    } else {
        wbook_print_info(book);
    }

    if (!err) {
        if (gui) {
            /* the C runs the sheet-selection dialog here; with it dropped we only
               honour the debugging flag */
            if (book_debugging(book)) {
                debug_print = 1;
                print_version();
            }
        } else {
            if (getenv("IMPORT_DEBUG") != null) {
                debug_print = 1;
            }
            err = wbook_check_params(book);
            if (err) {
                gretl_errmsg_set(_("Invalid argument for worksheet import").ptr);
            }
        }
    }

    dbprintf("sheet selected=%d; import offsets: col=%d, row=%d\n",
             book.selected, book.col_offset, book.row_offset);

    if (book.selected == -1) {
        /* canceled */
        err = -1;
    }

    if (err) goto getout;

    /* processing for specific worksheet */
    err = process_sheet(fname, book, xi, prn);

    if (err) {
        const(char) * buf = gretl_print_get_buffer(prn);

        if (*buf == 0) {
            pputs(prn, _("Failed to process Excel file").ptr);
            buf = gretl_print_get_buffer(prn);
        }
        fprintf(stderr, "%s\n", buf);
        goto getout;
    }

    /* get sizes and locate any blank columns */
    err = get_sheet_dimensions(book, xi, prn);

    if (err) goto getout;

    /* check feasibility of offsets */
    if (book.row_offset >= xi.nrows) {
        pputs(prn, _("Starting row is out of bounds.\n").ptr);
        err = 1;
    } else if (book.col_offset >= xi.totcols) {
        pputs(prn, _("Starting column is out of bounds.\n").ptr);
        err = 1;
    }

    if (err) goto getout;

    if (first_col_strings(book, xi)) {
        puts("found label strings in first imported column");
    } else if (book_numeric_dates(book)) {
        puts("found calendar dates in first imported column");
    } else {
        puts("check for label strings in first imported column: not found");
    }

    /* any bad or missing variable names? */
    err = check_all_varnames(book, xi, prn);

    if (err == VarnameError.null_ || err == VarnameError.notstr) {
        pputs(prn, _("One or more variable names are missing.\n").ptr);
        pputs(prn, adjust_rc);
    } else if (err == VarnameError.none) {
        pputs(prn, _("it seems there are no variable names\n").ptr);
        book_set_auto_varnames(book);
        book.row_offset -= 1;
        err = 0;
    }

    if (err) goto getout;

    /* any bad data? */
    err = check_data_block(book, xi, &missvals, &strerr);

    if (err) {
        pprintf(prn, _("Expected numeric data, found string:\n"
                       ~ "%s\" at row %d, column %d\n").ptr,
                strerr.str, strerr.row, strerr.column);
        free(strerr.str);
        pputs(prn, adjust_rc);
        goto getout;
    } else if (missvals) {
        pputs(prn, _("Warning: there were missing values\n").ptr);
    }

    r0 = book.row_offset;
    c0 = book.col_offset;
    newset.n = xi.nrows - 1 - r0;

    if (book_numeric_dates(book) ||
        (!book_auto_varnames(book) && import_obs_label(xls_cell(xi, r0, c0)))) {
        char ** labels = labels_array(xi, r0 + 1, c0, newset);

        if (labels != null) {
            pd = importer_dates_check(labels, &book.flags, newset, prn, &err);
            free(labels);
            if (pd < 0) {
                book_unset_numeric_dates(book);
            }
        }

        if (pd > 0) {
            /* got time-series info from dates/labels */
            book_time_series_setup(book, newset, pd);
            ts_markers = newset.markers;
            ts_S = newset.S;
        } else if (!book_numeric_dates(book) &&
                   alpha_cell(xls_cell(xi, r0, c0)) &&
                   col0_is_numeric(xi, r0, c0)) {
            book_unset_obs_labels(book);
            maybe_revise_xls_codelist(xi);
        }
    }

    /* dimensions of the dataset */
    newset.v = n_vars_from_col(book, xi.totcols, xi.blank_col);
    fprintf(stderr, "newset->v = %d, newset->n = %d\n",
            newset.v, newset.n);

    /* create import dataset */
    err = worksheet_start_dataset(newset);
    if (err) {
        goto getout;
    }

    if (book_time_series(book)) {
        newset.markers = cast(char) ts_markers;
        newset.S = ts_S;
    } else {
        dataset_obs_info_default(newset);
    }

    /* OK: actually populate the dataset */
    err = transcribe_data(book, xi, newset, prn);
    if (err) {
        goto getout;
    }

    if (fix_varname_duplicates(newset)) {
        pputs(prn, _("warning: some variable names were duplicated\n").ptr);
    }

    if (book_obs_labels(book)) {
        dataset_allocate_obs_markers(newset);
        if (newset.S != null) {
            i = book.col_offset;
            for (t=0; t<newset.n; t++) {
                int ts = t + 1 + book.row_offset;
                char * src = xls_cell(xi, ts, i);

                if (src != null) {
                    gretl_utf8_strncat_trim(newset.S[t], src + 1, OBSLEN - 1);
                }
            }
        }
    }

    if (book.flags & BookFlag.data_reversed) {
        reverse_data(newset, prn);
    }

    if (!err && xi.st != null) {
        err = gretl_string_table_validate(xi.st, Gretlopt.OPT_S);
        if (err) {
            pputs(prn, _("Failed to interpret the data as numeric\n").ptr);
        } else {
            gretl_string_table_finalize(xi.st, newset);
        }
    }

    err = merge_or_replace_data(dset, &newset, get_merge_opts(opt), prn);

    if (!err && !merge) {
        dataset_add_import_info(dset, fname, GretlFileType.xls);
    }

    /* the C has "if (!err && gui) wbook_record_params(book, list)" here: only meaningful
       with the (dropped) GUI dialog, but kept for fidelity */
    if (!err && gui) {
        wbook_record_params(book, list);
    }

 getout:

    free_xls_info(xi);
    wbook_free(book);
    gretl_pop_c_numeric_locale();

    if (newset != null) {
        destroy_dataset(newset);
    }

    return err;
}
