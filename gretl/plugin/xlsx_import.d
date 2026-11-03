/* Port of plugin/xlsx_import.c: importer for OOXML (.xlsx) workbooks.

   In the C this is a gretl plugin (loaded by get_plugin_function("xlsx_get_data")); here it is an
   ordinary module that plugins.d returns for that name. The GTK sheet-selection dialog is dropped:
   xlsx_get_data() keeps its C signature (including the unused trailing GtkWidget * argument, now
   void *), and the OPT_G ("GUI") branch behaves like the command-line one: the sheet and the column
   and row offsets come from @sheetname and @list only. */

module gretl.plugin.xlsx_import;


import core.stdc.ctype : isalpha;
import core.stdc.stdio : fprintf, snprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : atof, atoi, calloc, free;
import core.stdc.string : strcmp, strcpy, strncat, strncmp;

import gretl.calendar : MS_excel_date_string;
import gretl.csvdatah : NON_NUMERIC;
import gretl.csvdata : import_na_string, import_obs_label, non_numeric_check;
import gretl.dataio : check_atof, dataset_add_import_info, get_merge_opts, merge_or_replace_data;
import gretl.dataioh : GretlFileType;
import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_errors : gretl_errmsg_get, gretl_errmsg_set, gretl_errmsg_sprintf, maybe_save_gretl_errmsg;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : in_gretl_list;
import gretl.gretl_prn : Prn, gretl_print_flush_stream, pprintf, pputc, pputs;
import gretl.gretl_string_table : GretlStringTable, gretl_string_table_destroy, gretl_string_table_finalize,
    gretl_string_table_index, gretl_string_table_validate;
import gretl.gretl_xml : gretl_xml_open_doc_root;
import gretl.glib : gchar;
import gretl.gmem : g_free;
import gretl.gstrfuncs : g_strdup_printf, g_strstrip;
import gretl.plugin.import_common;
import gretl.plugin.importerh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.strutils : fix_varname_duplicates, gretl_strdup, gretl_utf8_strncat_trim, integer_string,
    strings_array_add, strings_array_free, strings_array_new;
import gretl.strutilsh : SLASH;
import gretl.xmlh;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;

private alias XUC = const(xmlChar) *;

/* unported dependencies: per-file declarations from the 2026c headers */
/* possibly useful reference:
   http://officeopenxml.com/SScontentOverview.php
*/

private struct XlsxInfo {
    BookFlag flags;
    Gretlopt opt;
    int trydates;
    int maxrow;
    int maxcol;
    int namerow;
    int obscol;
    int xoffset;
    int yoffset;
    char[FILENAME_MAX] sheetfile;
    char[FILENAME_MAX] stringsfile;
    int n_sheets;
    char ** sheetnames;
    char ** filenames;
    int selsheet;
    int n_strings;
    char ** strings;
    DatasetRaw * dset;
    int * codelist;
    GretlStringTable * st;
}

private enum CellType {
    none,
    number,
    stringref,
    formula,
    istring
}

private void xlsx_info_init(XlsxInfo * xinfo, Gretlopt opt)
{
    xinfo.flags = BookFlag.top_left_empty;
    xinfo.opt = opt;
    xinfo.trydates = 0;
    xinfo.maxrow = 0;
    xinfo.maxcol = 0;
    xinfo.namerow = -1;
    xinfo.obscol = -1;
    xinfo.xoffset = 0;
    xinfo.yoffset = 0;
    xinfo.sheetfile[0] = '\0';
    xinfo.stringsfile[0] = '\0';
    xinfo.n_sheets = 0;
    xinfo.sheetnames = null;
    xinfo.filenames = null;
    xinfo.selsheet = 0;
    xinfo.n_strings = 0;
    xinfo.strings = null;
    xinfo.dset = null;
    xinfo.codelist = null;
    xinfo.st = null;
}

private void xlsx_info_reset(XlsxInfo * xinfo, Gretlopt opt)
{
    xinfo.flags = BookFlag.top_left_empty;
    xinfo.opt = opt;
    xinfo.trydates = 0;
    xinfo.maxrow = 0;
    xinfo.maxcol = 0;
    xinfo.namerow = -1;
    xinfo.obscol = -1;
    strings_array_free(xinfo.strings, xinfo.n_strings);
    xinfo.strings = null;
    xinfo.n_strings = 0;
    destroy_dataset(xinfo.dset);
    xinfo.dset = null;
    free(xinfo.codelist);
    xinfo.codelist = null;
    if (xinfo.st != null) {
        gretl_string_table_destroy(xinfo.st);
        xinfo.st = null;
    }
}

private void xlsx_info_free(XlsxInfo * xinfo)
{
    if (xinfo != null) {
        strings_array_free(xinfo.sheetnames, xinfo.n_sheets);
        strings_array_free(xinfo.filenames, xinfo.n_sheets);
        strings_array_free(xinfo.strings, xinfo.n_strings);
        destroy_dataset(xinfo.dset);
        free(xinfo.codelist);
        if (xinfo.st != null) {
            gretl_string_table_destroy(xinfo.st);
        }
    }
}

/* Parse what we need from sharedStrings.xml */

private int xlsx_read_shared_strings(XlsxInfo * xinfo, Prn * prn)
{
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    xmlNodePtr val;
    char * tmp;
    int i, n = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(xinfo.stringsfile.ptr, "sst",
                                  &doc, &cur);

    if (err) {
        pprintf(prn, _("Couldn't find shared strings table\n").ptr);
        pprintf(prn, "%s", gretl_errmsg_get());
        return err;
    }

    tmp = cast(char *) xmlGetProp(cur, cast(XUC) "uniqueCount");
    if (tmp == null) {
        tmp = cast(char *) xmlGetProp(cur, cast(XUC) "count");
    }

    if (tmp == null) {
        pprintf(prn, _("didn't get sst count\n").ptr);
        err = GretlError.data;
    } else {
        n = atoi(tmp);
        if (n <= 0) {
            pprintf(prn, _("didn't get valid sst count\n").ptr);
            err = GretlError.data;
        }
        free(tmp);
    }

    if (!err) {
        xinfo.strings = strings_array_new(n);
        if (xinfo.strings == null) {
            err = GretlError.alloc;
        }
    }

    cur = cur.children;

    /* The strings in an <sst> are mostly set up as

       <si><t>XXX</t></si>
       <si><t>YYY</t></si> ...

       But there are also weird cases where junk is interposed
       and the structure becomes

       <si><r>...<t>XXX</t></r><r>...<t>YYY</t></r></si> ...

       That is, an <si> element may contain more than one <r>
       element, which embeds a <t> along with formatting crap.
    */

    i = 0;
    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "si")) {
            int gotstr = 0;

            val = cur.children;
            while (val != null && !err && !gotstr) {
                if (!xmlStrcmp(val.name, cast(XUC) "t")) {
                    /* got a regular <t> element */
                    tmp = cast(char *) xmlNodeGetContent(val);
                    if (tmp == null) {
                        pprintf(prn, _("failed reading string %d\n").ptr, i);
                        err = GretlError.data;
                    } else {
                        xinfo.strings[i++] = g_strstrip(tmp);
                        gotstr = 1;
                    }
                } else if (!xmlStrcmp(val.name, cast(XUC) "r")) {
                    /* hunt for <t> inside an <r> element */
                    xmlNodePtr sub = val.children;

                    while (sub != null && !err && i < n) {
                        if (!xmlStrcmp(sub.name, cast(XUC) "t")) {
                            tmp = cast(char *) xmlNodeGetContent(sub);
                            if (tmp == null) {
                                pprintf(prn, _("failed reading string %d\n").ptr, i);
                                err = GretlError.data;
                            } else {
                                xinfo.strings[i++] = g_strstrip(tmp);
                                gotstr = 1;
                            }
                        }
                        sub = sub.next;
                    }
                }
                val = val.next;
            }
        }
        if (i == n) {
            break;
        }
        cur = cur.next;
    }

    if (!err && i < n) {
        pprintf(prn, _("expected %d shared strings but only found %d\n").ptr,
                n, i);
        err = GretlError.data;
    }

    if (!err) {
        xinfo.n_strings = i;
    } else if (xinfo.strings != null) {
        strings_array_free(xinfo.strings, n);
        xinfo.strings = null;
    }

    xmlFreeDoc(doc);

    return err;
}

/* Given the string representation of the shared-string index for
   a cell with a string value, look up the target string. If the
   shared strings XML file has not yet been read, its reading is
   triggered.
*/

private const(char) * xlsx_string_value(const(char) * idx, XlsxInfo * xinfo,
                                        Prn * prn)
{
    const(char) * ret = null;
    int err = 0;

    if (xinfo.n_strings == 0) {
        err = xlsx_read_shared_strings(xinfo, prn);
    }

    if (!err) {
        int i = atoi(idx);

        if (i >= 0 && i < xinfo.n_strings) {
            ret = xinfo.strings[i];
        }
    }

    return ret;
}

private int letter_val(char c)
{
    static immutable string UC = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int i;

    for (i=0; i<26; i++) {
        if (c == UC[i]) {
            return i+1;
        }
    }

    return 0;
}

/* Given an Excel cell reference such as "AB65" split it out
   into 1-based row and column indices.
*/

private int xlsx_cell_get_coordinates(const(char) * s,
                                      int * row,
                                      int * col)
{
    char[8] colref = 0;
    int i, k, v;
    int err = 0;

    for (i=0; i<7; i++) {
        if (isalpha(s[i])) {
            colref[i] = s[i];
        } else {
            break;
        }
    }

    colref[i] = '\0';

    *row = atoi(s + i--);
    *col = 0;

    k = 1;
    while (i >= 0 && !err) {
        v = letter_val(colref[i]);
        if (v == 0) {
            err = GretlError.data;
        } else {
            *col += k * v;
            k *= 26;
        }
        i--;
    }

    return err;
}

private int xlsx_set_varname(XlsxInfo * xinfo, int i, const(char) * s,
                             int row, int col, Prn * prn)
{
    int err = 0;

    if (i == -1) {
        return 0;
    }

    if (i < 1 || i >= xinfo.dset.v) {
        fprintf(stderr, "error in xlsx_set_varname: i = %d\n", i);
        err = GretlError.data;
    } else {
        *xinfo.dset.varname[i] = '\0';
        if (s != null) {
            strncat(xinfo.dset.varname[i], s, VNAMELEN - 1);
        }
        err = check_imported_varname(xinfo.dset.varname[i],
                                     i, row, col, prn);
    }

    return err;
}

private void xlsx_real_set_obs_string(XlsxInfo * xinfo, int t, const(char) * s)
{
    *xinfo.dset.S[t] = '\0';
    if (s != null) {
        gretl_utf8_strncat_trim(xinfo.dset.S[t], s, OBSLEN - 1);
    }
}

private int xlsx_set_obs_string(XlsxInfo * xinfo, int row, int col,
                                int t, const(char) * s, Prn * prn)
{
    int err = 0;

    if (t == -1) {
        return 0;
    }

    if (xinfo.dset.S == null) {
        /* "can't happen" */
        fprintf(stderr, "error in xlsx_set_obs_string: S not allocated!\n");
        pprintf(prn, (_("Expected numeric data, found string:\n")
                      ~ _("'%s' at row %d, column %d\n")).ptr, s, row, col);
        err = GretlError.data;
    } else if (t < 0 || t >= xinfo.dset.n) {
        fprintf(stderr, "error in xlsx_set_obs_string: t = %d\n", t);
        err = GretlError.data;
    } else {
        xlsx_real_set_obs_string(xinfo, t, s);
    }

    return err;
}

private int xlsx_set_value(XlsxInfo * xinfo, int i, int t, double x)
{
    if (i == -1 && t >= 0 && t < xinfo.dset.n) {
        /* maybe this should really be an obs label? */
        if (xinfo.flags & BookFlag.obs_labels) {
            gchar * tmp = g_strdup_printf("%g", x);

            xlsx_real_set_obs_string(xinfo, t, tmp);
            g_free(tmp);
            xinfo.trydates = 1;
            return 0;
        }
    }

    if (i == -1 || t == -1) {
        return 0;
    }

    if (i < 1 || i >= xinfo.dset.v ||
        t < 0 || t >= xinfo.dset.n) {
        fprintf(stderr, "error in xlsx_set_value: i = %d, t = %d\n", i, t);
        return GretlError.data;
    } else {
        xinfo.dset.Z[i][t] = x;
        return 0;
    }
}

/* for use on third xlsx pass, if applicable */

private int xlsx_handle_stringval3(XlsxInfo * xinfo,
                                   int i, int t,
                                   const(char) * s,
                                   Prn * prn)
{
    const(char) * strval;
    gchar * tmp = null;
    int ix, err = 0;

    if (xinfo.dset.Z[i][t] == NON_NUMERIC) {
        strval = s;
    } else if (!na(xinfo.dset.Z[i][t])) {
        tmp = g_strdup_printf("%g", xinfo.dset.Z[i][t]);
        strval = tmp;
    } else {
        return 0;
    }

    ix = gretl_string_table_index(xinfo.st, strval, i, 0, prn);
    if (ix > 0) {
        xinfo.dset.Z[i][t] = cast(double) ix;
    } else {
        err = GretlError.data;
    }

    g_free(tmp);

    return err;
}

/* for use on second xlsx pass */

private int xlsx_handle_stringval(XlsxInfo * xinfo, int i, int t,
                                  const(char) * s, int r, int c,
                                  Prn * prn)
{
    double xit = NADBL;
    int err = 0;

    if (import_na_string(s)) {
        ; /* OK */
    } else if (*s == '\0') {
        ; /* well, OK */
    } else if (!strncmp(s, "$$ER:", 5)) {
        ; /* failed accessing external data? */
    } else {
        xit = NON_NUMERIC;
    }

    if (i < 1 || i >= xinfo.dset.v ||
        t < 0 || t >= xinfo.dset.n) {
        fprintf(stderr, "error in xlsx_handle_stringval: i=%d, t=%d, s='%s'\n",
                i, t, s);
        err = GretlError.data;
    } else {
        xinfo.dset.Z[i][t] = xit;
    }

    return err;
}

private int xlsx_var_index(XlsxInfo * xinfo, int col)
{
    int i = col - xinfo.xoffset - 1;

    if (i == 0 && (xinfo.flags & BookFlag.obs_labels)) {
        i = -1; /* the first column holds labels */
    } else if (i >= 0 && !(xinfo.flags & BookFlag.obs_labels)) {
        i++; /* skip the constant in position 0 */
    }

    return i;
}

private int xlsx_obs_index(XlsxInfo * xinfo, int row)
{
    int t = row - xinfo.yoffset - 1;

    if (!(xinfo.flags & BookFlag.auto_varnames)) {
        t--; /* the first row holds varnames */
    }

    return t;
}

private int may_have_obs_label(XlsxInfo * xinfo)
{
    if (xinfo.opt & Gretlopt.OPT_A) {
        /* --all-cols */
        return 0;
    } else {
        return 1;
    }
}

private void xlsx_check_top_left(XlsxInfo * xinfo, int r, int c,
                                 const(char) * s, double x,
                                 Prn * myprn)
{
    if (r == xinfo.yoffset + 1 && c == xinfo.xoffset + 1) {
        /* We're in the top left cell of the reading area:
           this could be blank, or could hold the first
           varname, could hold "obs" or similar, or could
           be the first numerical value.
        */
        if (!na(x)) {
            /* Got a valid numerical value: that means we don't
               have variable names on the top row.
            */
            xinfo.flags |= BookFlag.auto_varnames;
        } else if (may_have_obs_label(xinfo) && import_obs_label(s)) {
            /* blank or "obs" or similar */
            xinfo.flags |= BookFlag.obs_labels;
            xinfo.obscol = c;
        }
        if (!na(x) || s != null) {
            /* record the fact that the top-left corner is not empty */
            xinfo.flags &= ~BookFlag.top_left_empty;
        }
    } else if (r == xinfo.yoffset + 1 && c == xinfo.xoffset + 2) {
        /* first row, second column */
        if (!na(x)) {
            /* got a number, not a varname */
            xinfo.flags |= BookFlag.auto_varnames;
        } else {
            pprintf(myprn, "set namerow = %d\n", r);
            xinfo.namerow = r;
        }
    }
}

/* Minimal handling of formulae: we're really just looking for the
   case where dates in the observations column are defined as
   previous value plus constant, or minor variations on that
   notion.
*/

private void xlsx_maybe_handle_formula(XlsxInfo * xinfo,
                                       const(char) * src,
                                       int i, int t)
{
    if (i == -1 && t > 0 && (xinfo.flags & BookFlag.obs_labels)) {
        char c;
        char[8] cref = 0;
        int k;

        if (sscanf(src, "%7[^+-]%c%d", cref.ptr, &c, &k) == 3) {
            int err, st = 0, col = 0;

            if (c == '-') {
                k = -k;
            } else if (c != '+') {
                /* we'll only only handle +/- */
                return;
            }
            err = xlsx_cell_get_coordinates(cref.ptr, &st, &col);
            if (err || col != xinfo.xoffset + 1) {
                return;
            }
            st = st - xinfo.yoffset - 1;
            if (!(xinfo.flags & BookFlag.auto_varnames)) {
                st--;
            }
            if (st >= 0 && st < t) {
                const(char) * s = xinfo.dset.S[st];

                if (integer_string(s)) {
                    gchar * tmp = g_strdup_printf("%d", atoi(s) + k);

                    xlsx_real_set_obs_string(xinfo, t, tmp);
                    g_free(tmp);
                }
            }
        }
    }
}

private void xlsx_set_dims(XlsxInfo * xinfo, int r, int c)
{
    if (r > xinfo.maxrow) {
        xinfo.maxrow = r;
    }
    if (c > xinfo.maxcol) {
        xinfo.maxcol = c;
    }
}

private char * node_get_inline_string(xmlNodePtr val)
{
    xmlNodePtr ist = val.children;
    char * ret = null;

    if (ist != null && !xmlStrcmp(ist.name, cast(XUC) "t")) {
        char * tmp = cast(char *) xmlNodeGetContent(ist);

        if (tmp != null) {
            ret = gretl_strdup(g_strstrip(tmp));
            free(tmp);
        }
    }

    return ret;
}

private int has_formula_child(xmlNodePtr val)
{
    int ret = 0;

    while (val) {
        if (!xmlStrcmp(val.name, cast(XUC) "f")) {
            ret = 1;
            break;
        }
        val = val.next;
    }

    return ret;
}

private int get_cell_basics(xmlNodePtr cur,
                            char ** cref, int * celltype,
                            int * row, int * col,
                            Prn * prn)
{
    char * ctype = null;
    int err = 0;

    pputs(prn, " cell");

    *cref = cast(char *) xmlGetProp(cur, cast(XUC) "r");
    if (*cref == null) {
        pprintf(prn, ": %s\n", _("couldn't find 'r' property").ptr);
        return GretlError.data;
    }

    err = xlsx_cell_get_coordinates(*cref, row, col);
    if (err) {
        pprintf(prn, ": %s\n", _("couldn't find coordinates").ptr);
        return GretlError.data;
    }

    ctype = cast(char *) xmlGetProp(cur, cast(XUC) "t");
    if (ctype == null) {
        /* default to numeric? */
        ctype = gretl_strdup("n");
    }

    pprintf(prn, "(%d, %d; %s %s)", *row, *col, _("type").ptr, ctype);

    if (!strcmp(ctype, "n") || !strcmp(ctype, "b")) {
        /* numeric or boolean (0/1) */
        *celltype = CellType.number;
    } else if (!strcmp(ctype, "d")) {
        /* date: try numeric? */
        *celltype = CellType.number;
    } else if (!strcmp(ctype, "s")) {
        /* reference to entry in string table */
        *celltype = CellType.stringref;
    } else if (!strcmp(ctype, "str")) {
        /* string representing a formula */
        *celltype = CellType.formula;
    } else if (!strcmp(ctype, "inlineStr")) {
        /* in-cell string literal (or maybe formula!) */
        if (has_formula_child(cur.children)) {
            *celltype = CellType.formula;
        } else {
            *celltype = CellType.istring;
        }
    }

    free(ctype);

    return err;
}

/* Read the cells in a given row: the basic info we want from
   each cell is its reference ("r"), e.g. "A2"; its type
   ("t"), e.g. "s", if present; and its value, which is
   not a property but a sub-element "<v>...</v>".
*/

private int xlsx_read_row(xmlNodePtr cur, XlsxInfo * xinfo,
                          int i, Prn * prn)
{
    Prn * myprn = null;
    xmlNodePtr val;
    char * tmp;
    int row = -1, col = -1;
    int pass, empty = 1;
    int colmax = 0;
    int err = 0;

    if (xinfo.codelist != null) {
        pass = 3;
    } else if (xinfo.dset != null) {
        /* if we never found anything in the (possibly notional)
           top-left cell, set the observations column number
        */
        if (xinfo.flags & BookFlag.top_left_empty) {
            xinfo.obscol = xinfo.xoffset + 1;
        }
        pass = 2;
    } else {
        pass = 1;
    }

    cur = cur.children;

    /* loop across cells in row */

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "c")) {
            /* we got a cell in the given row */
            char * cref = null;
            char * formula = null;
            const(char) * strval = null;
            double xval = NADBL;
            int celltype = CellType.none;
            int gotv = 0;
            int gotf = 0;

            err = get_cell_basics(cur, &cref, &celltype, &row, &col, myprn);
            if (err) {
                break;
            }

            if (pass >= 2 && (row > xinfo.maxrow || col > xinfo.maxcol)) {
                goto skipit;
            }

            val = cur.children;

            /* try to find a value in the current row/cell */

            while (val && !err && !gotv) {
                if (!xmlStrcmp(val.name, cast(XUC) "v")) {
                    /* value element */
                    int bogus = 0;

                    tmp = cast(char *) xmlNodeGetContent(val);
                    if (tmp != null) {
                        if (celltype == CellType.number) {
                            pprintf(myprn, " %s = %s\n", _("value").ptr, tmp);
                            if (*tmp != '\0' && check_atof(tmp) == 0) {
                                xval = atof(tmp);
                            }
                        } else if (celltype == CellType.stringref) {
                            /* look up string table */
                            strval = xlsx_string_value(tmp, xinfo, prn);
                            if (strval == null) {
                                pprintf(myprn, " %s = NULL\n", _("value").ptr);
                                err = GretlError.data;
                            } else if (*strval == '\0') {
                                pprintf(myprn, " %s = empty (CELL_STRINGREF)\n",
                                        _("value").ptr, strval);
                                bogus = 1; /* 2025-12-15: is this OK? */
                            } else {
                                pprintf(myprn, " %s = '%s'\n", _("value").ptr, strval);
                            }
                        }
                        free(tmp);
                        if (!bogus) {
                            gotv = 1;
                        }
                    }
                } else if (celltype == CellType.formula &&
                           !xmlStrcmp(val.name, cast(XUC) "f")) {
                    /* formula element */
                    formula = cast(char *) xmlNodeGetContent(val);
                    if (formula != null) {
                        gotf = 1;
                    }
                } else if (celltype == CellType.istring &&
                           !xmlStrcmp(val.name, cast(XUC) "is")) {
                    /* string element */
                    strval = node_get_inline_string(val);
                    if (strval != null) {
                        gotv = 1;
                    }
                }
                /* proceed to next child of given cell */
                val = val.next;
            }

            if (!err && xinfo.dset == null) {
                /* on the first pass, check for obs column, varname status */
                xlsx_check_top_left(xinfo, row, col, strval, xval, myprn);
            }

            if (err) {
                pprintf(myprn, ": (%s) %s", cref, _("error").ptr);
            } else if (!gotv) {
                pprintf(myprn, ": (%s) %s", cref, _("no data value").ptr);
                if (gotf) {
                    pprintf(myprn, ": %s = '%s'\n", _("formula").ptr, formula);
                } else {
                    pputc(myprn, '\n');
                }
            }

            if (gotv || gotf) {
                colmax = col;
            }

            if (!err && xinfo.dset != null &&
                col > xinfo.xoffset &&
                row > xinfo.yoffset) {
                int vi = xlsx_var_index(xinfo, col); /* C: "i", shadowing the parameter */
                int t = xlsx_obs_index(xinfo, row);
                int strcell = celltype == CellType.stringref || celltype == CellType.istring;

                /* here we're on the second or possibly third pass,
                   with a dataset allocated */

                if (pass == 3) {
                    if (vi > 0 && t >= 0) {
                        if (in_gretl_list(xinfo.codelist, vi)) {
                            xlsx_handle_stringval3(xinfo, vi, t, strval, prn);
                        }
                    }
                } else if (strcell) {
                    if (row == xinfo.namerow) {
                        err = xlsx_set_varname(xinfo, vi, strval, row, col, prn);
                    } else if (col == xinfo.obscol) {
                        err = xlsx_set_obs_string(xinfo, row, col, t, strval, prn);
                    } else if (strval != null) {
                        err = xlsx_handle_stringval(xinfo, vi, t, strval,
                                                    row, col, prn);
                    }
                    if (celltype == CellType.istring) {
                        /* finished with copy of string literal */
                        free(cast(char *) strval);
                    }
                } else if (row == xinfo.namerow) {
                    err = xlsx_set_varname(xinfo, vi, strval, row, col, prn);
                } else if (gotv) {
                    err = xlsx_set_value(xinfo, vi, t, xval);
                } else if (gotf) {
                    xlsx_maybe_handle_formula(xinfo, formula, vi, t);
                }
            } else if (celltype == CellType.istring && strval != null) {
                /* finished with copy of string literal */
                free(cast(char *) strval);
            }

            if (gotv || gotf) {
                empty = 0;
            }

        skipit:

            free(cref);
            free(formula);
        }

        /* move onto next cell in row */
        cur = cur.next;
    } /* end loop across cells in row */

    if (!err) {
        if (empty) {
            pputs(myprn, " xlsx_read_row: empty row!\n");
        } else if (pass == 1) {
            xlsx_set_dims(xinfo, row, colmax);
        }
    }

    if (err) {
        fprintf(stderr, "xlsx_read_row: returning %d\n", err);
    }

    return err;
}

private int xlsx_check_dimensions(XlsxInfo * xinfo, Prn * prn)
{
    int v = xinfo.maxcol - xinfo.xoffset;
    int n = xinfo.maxrow - xinfo.yoffset;
    int err = 0;

    if (xinfo.flags & BookFlag.top_left_empty) {
        xinfo.flags |= BookFlag.obs_labels;
    }

    if (xinfo.flags & BookFlag.obs_labels) {
        /* subtract a column for obs labels */
        v--;
    }

    if (!(xinfo.flags & BookFlag.auto_varnames)) {
        /* subtract a row for varnames */
        n--;
    }

    if (v < 0) v = 0;
    if (n < 0) n = 0;

    pprintf(prn, _("Found %d variables and %d observations\n").ptr, v, n);

    if (v == 0 || n == 0) {
        pputs(prn, _("File contains no data\n").ptr);
        err = GretlError.data;
    } else {
        int labels = (xinfo.flags & BookFlag.obs_labels) != 0;

        xinfo.dset = create_new_dataset(v + 1, n, labels);
        if (xinfo.dset == null) {
            err = GretlError.alloc;
        } else if (xinfo.flags & BookFlag.auto_varnames) {
            /* write fallback variable names */
            int i;

            for (i=1; i<=v; i++) {
                sprintf(xinfo.dset.varname[i], "v%d", i);
            }
        }
    }

    return err;
}

private int xlsx_non_numeric_check(XlsxInfo * xinfo, Prn * prn)
{
    GretlStringTable * st = null;
    int * nlist = null;
    int err = 0;

    err = non_numeric_check(xinfo.dset, &nlist, &st, prn);

    if (!err) {
        xinfo.codelist = nlist;
        xinfo.st = st;
    }

    return err;
}

private int xlsx_read_worksheet(XlsxInfo * xinfo,
                                const(char) * fname,
                                Prn * prn)
{
    xmlDocPtr doc = null;
    xmlNodePtr data_node = null;
    xmlNodePtr cur = null;
    xmlNodePtr c1;
    int i, gotdata = 0;
    int err = 0;

    if (snprintf(xinfo.sheetfile.ptr, FILENAME_MAX, "xl%c%s", SLASH,
                 xinfo.filenames[xinfo.selsheet]) >= FILENAME_MAX) {
        fprintf(stderr, "xlsx_read_worksheet: internal sheet path too long\n");
        return GretlError.data;
    }

    sprintf(xinfo.stringsfile.ptr, "xl%csharedStrings.xml", SLASH);

    err = gretl_xml_open_doc_root(xinfo.sheetfile.ptr, "worksheet",
                                  &doc, &cur);

    if (err) {
        pprintf(prn, _("didn't get worksheet\n").ptr);
        pprintf(prn, "%s", gretl_errmsg_get());
        return err;
    }

    /* walk the tree, first pass */
    i = 0;
    cur = cur.children;
    while (cur != null && !err && !gotdata) {
        if (!xmlStrcmp(cur.name, cast(XUC) "sheetData")) {
            data_node = c1 = cur.children;
            while (c1 != null && !err) {
                if (!xmlStrcmp(c1.name, cast(XUC) "row")) {
                    err = xlsx_read_row(c1, xinfo, i++, prn);
                }
                c1 = c1.next;
            }
            gotdata = 1;
        }
        cur = cur.next;
    }

    if (!err && xinfo.dset == null) {
        err = xlsx_check_dimensions(xinfo, prn);
        if (!err) {
            /* second pass: get actual data */
            i = 0;
            gretl_push_c_numeric_locale();
            c1 = data_node;
            while (c1 != null && !err) {
                if (!xmlStrcmp(c1.name, cast(XUC) "row")) {
                    err = xlsx_read_row(c1, xinfo, i++, prn);
                }
                c1 = c1.next;
            }
            gretl_pop_c_numeric_locale();
        }
    }

    if (!err) {
        err = xlsx_non_numeric_check(xinfo, prn);
        if (!err && xinfo.codelist != null) {
            /* third pass, if needed: get string-valued vars */
            c1 = data_node;
            i = 0;
            while (c1 != null && !err) {
                if (!xmlStrcmp(c1.name, cast(XUC) "row")) {
                    err = xlsx_read_row(c1, xinfo, i++, prn);
                }
                c1 = c1.next;
            }
            if (!err) {
                err = gretl_string_table_validate(xinfo.st, Gretlopt.OPT_S);
                if (err) {
                    pputs(prn, _("Failed to interpret the data as numeric\n").ptr);
                } else {
                    gretl_string_table_finalize(xinfo.st, xinfo.dset);
                }
            }
        }
    }

    xmlFreeDoc(doc);

    return err;
}

/* A quick check to see if a worksheet XML file contains
   any actual data. */

private int xlsx_sheet_has_data(const(char) * fname)
{
    xmlDocPtr doc = null;
    xmlNodePtr c1, cur = null;
    gchar * fullname;
    int err, ret = 0;

    fullname = g_strdup_printf("xl%c%s", SLASH, fname);

    err = gretl_xml_open_doc_root(fullname, "worksheet",
                                  &doc, &cur);
    if (!err) {
        cur = cur.children;
        while (cur != null && ret == 0) {
            if (!xmlStrcmp(cur.name, cast(XUC) "sheetData")) {
                c1 = cur.children;
                while (c1 != null && ret == 0) {
                    if (!xmlStrcmp(c1.name, cast(XUC) "row")) {
                        ret = 1;
                    }
                    c1 = c1.next;
                }
            }
            cur = cur.next;
        }
        xmlFreeDoc(doc);
    }

    if (!ret) {
        fprintf(stderr, "%s: contains no data\n", fname);
    }

    g_free(fullname);

    return ret;
}

private int xlsx_workbook_get_sheetnames(XlsxInfo * xinfo,
                                         const(char) * fname)
{
    xmlDocPtr doc = null;
    xmlNodePtr c1, cur = null;
    char * ID;
    char * sheetname;
    int ns = 0, found = 0;
    int err;

    err = gretl_xml_open_doc_root(fname, "workbook",
                                  &doc, &cur);
    if (!err) {
        cur = cur.children;
        while (cur != null && !found) {
            if (!xmlStrcmp(cur.name, cast(XUC) "sheets")) {
                c1 = cur.children;
                while (c1 != null) {
                    if (!xmlStrcmp(c1.name, cast(XUC) "sheet")) {
                        ID = cast(char *) xmlGetProp(c1, cast(XUC) "id");
                        sheetname = cast(char *) xmlGetProp(c1, cast(XUC) "name");
                        if (ID != null && sheetname != null) {
                            strings_array_add(&xinfo.sheetnames,
                                              &xinfo.n_sheets,
                                              sheetname);
                            strings_array_add(&xinfo.filenames,
                                              &ns, ID);
                        }
                        free(ID);
                        free(sheetname);
                    }
                    c1 = c1.next;
                }
                found = 1;
            }
            cur = cur.next;
        }
        xmlFreeDoc(doc);
    }

    return err;
}

private void xlsx_expunge_sheet(XlsxInfo * xinfo, int i)
{
    int j;

    free(xinfo.sheetnames[i]);
    free(xinfo.filenames[i]);

    for (j=i; j<xinfo.n_sheets-1; j++) {
        xinfo.sheetnames[j] = xinfo.sheetnames[j+1];
        xinfo.filenames[j] = xinfo.filenames[j+1];
    }

    xinfo.n_sheets -= 1;
}

private int xlsx_match_sheet_id(XlsxInfo * xinfo,
                                const(char) * ID)
{
    int i;

    if (ID == null) {
        return -1;
    }

    for (i=0; i<xinfo.n_sheets; i++) {
        if (!strcmp(ID, xinfo.filenames[i])) {
            return i;
        }
    }

    return -1;
}

/* For the sheet names we got from the main workbook file:
   check that (1) we can track down the corresponding
   worksheet filename and (2) the worksheet file contains
   some data.
*/

private int xlsx_verify_sheets(XlsxInfo * xinfo, Prn * prn)
{
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    char * checker;
    char * ID;
    char * fname;
    int i, err;

    checker = cast(char *) calloc(xinfo.n_sheets, 1);
    if (checker == null) {
        return GretlError.alloc;
    }

    err = gretl_xml_open_doc_root("xl/_rels/workbook.xml.rels",
                                  "Relationships",
                                  &doc, &cur);
    if (!err) {
        cur = cur.children;
        while (cur != null) {
            if (!xmlStrcmp(cur.name, cast(XUC) "Relationship")) {
                ID = cast(char *) xmlGetProp(cur, cast(XUC) "Id");
                if ((i = xlsx_match_sheet_id(xinfo, ID)) >= 0) {
                    fname = cast(char *) xmlGetProp(cur, cast(XUC) "Target");
                    if (fname != null) {
                        if (xlsx_sheet_has_data(fname)) {
                            checker[i] = 1;
                            free(xinfo.filenames[i]);
                            xinfo.filenames[i] = fname;
                        } else {
                            free(fname);
                        }
                    }
                }
                free(ID);
            }
            cur = cur.next;
        }
        xmlFreeDoc(doc);
    }

    if (!err) {
        int j = 0;

        for (i=0; i<xinfo.n_sheets; i++) {
            if (checker[j++] == 0) {
                fprintf(stderr, "dropping sheet '%s'\n",
                        xinfo.sheetnames[i]);
                xlsx_expunge_sheet(xinfo, i--);
            }
        }
    }

    free(checker);

    return err;
}

private int xlsx_verify_specific_sheet(XlsxInfo * xinfo,
                                       int idx, Prn * prn)
{
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    char * ID;
    char * fname;
    int err;

    err = gretl_xml_open_doc_root("xl/_rels/workbook.xml.rels",
                                  "Relationships",
                                  &doc, &cur);
    if (!err) {
        const(char) * sname = xinfo.sheetnames[idx];
        int found = 0;

        cur = cur.children;
        while (cur != null && !found) {
            if (!xmlStrcmp(cur.name, cast(XUC) "Relationship")) {
                ID = cast(char *) xmlGetProp(cur, cast(XUC) "Id");
                if (xlsx_match_sheet_id(xinfo, ID) == idx) {
                    found = 1;
                    fname = cast(char *) xmlGetProp(cur, cast(XUC) "Target");
                    if (fname == null) {
                        pprintf(prn, _("'%s': couldn't find filename\n").ptr, sname);
                        err = GretlError.data;
                    } else if (xlsx_sheet_has_data(fname)) {
                        pprintf(prn, "'%s' -> %s\n", sname, fname);
                        free(xinfo.filenames[idx]);
                        xinfo.filenames[idx] = fname;
                    } else {
                        pprintf(prn, _("'%s': contains no data\n").ptr, sname);
                        err = GretlError.data;
                        free(fname);
                    }
                }
                free(ID);
            }
            cur = cur.next;
        }

        xmlFreeDoc(doc);

        if (!found) {
            pprintf(prn, _("'%s': couldn't find file Id\n").ptr, sname);
            err = GretlError.data;
        }

        if (!err) {
            /* record the pre-checked sheet selection */
            xinfo.selsheet = idx;
        }
    }

    return err;
}

private int xlsx_gather_sheet_names(XlsxInfo * xinfo,
                                    const(char) * insheet,
                                    int * list,
                                    Prn * prn)
{
    gchar * wb_name = g_strdup_printf("xl%cworkbook.xml", SLASH);
    int i, err;

    err = xlsx_workbook_get_sheetnames(xinfo, wb_name);
    g_free(wb_name);

    if (!err && xinfo.n_sheets == 0) {
        pputs(prn, _("Found no sheets\n").ptr);
        err = GretlError.data;
    }

    if (!err) {
        const(char) * sname;
        int have_insheet_name = 0;
        int target_sheet_number = 0;
        int insheet_idx = 0;

        if (insheet != null && *insheet != '\0') {
            /* try looking for a worksheet specified by name */
            have_insheet_name = 1;
            insheet_idx = -1; /* invalidate index */
        } else if (list != null && list[0] >= 1 && list[1] > 0) {
            /* or sheet specified by (1-based) sequence number */
            target_sheet_number = list[1];
            insheet_idx = -1; /* invalidate index */
        }

        fprintf(stderr, "Found these worksheets:\n");
        for (i=0; i<xinfo.n_sheets; i++) {
            sname = xinfo.sheetnames[i];
            fprintf(stderr, "%d: '%s' (%s)\n", i+1, sname, xinfo.filenames[i]);
            if (have_insheet_name) {
                if (insheet_idx < 0 && !strcmp(insheet, sname)) {
                    /* found the name we were looking for */
                    insheet_idx = i;
                }
            } else if (target_sheet_number == i + 1) {
                /* found the number we were looking for */
                insheet_idx = i;
            }
        }

        if (insheet_idx < 0) {
            /* a sheet was pre-specified but was not found */
            if (have_insheet_name && integer_string(insheet)) {
                /* try interpreting as plain integer index? */
                insheet_idx = atoi(insheet) - 1;
                if (insheet_idx < 0 || insheet_idx >= xinfo.n_sheets) {
                    /* no, it doesn't work */
                    insheet_idx = -1;
                }
            }
            if (insheet_idx < 0) {
                if (have_insheet_name) {
                    gretl_errmsg_sprintf(_("'%s': no such worksheet").ptr, insheet);
                } else {
                    gretl_errmsg_set(_("Invalid argument for worksheet import").ptr);
                }
                err = GretlError.data;
            }
        }

        if (!err) {
            if (have_insheet_name || target_sheet_number > 0) {
                /* avoid the expense of checking all sheets */
                err = xlsx_verify_specific_sheet(xinfo, insheet_idx, prn);
            } else {
                err = xlsx_verify_sheets(xinfo, prn);
                if (xinfo.n_sheets == 0) {
                    pputs(prn, _("\nFound no valid sheets\n").ptr);
                    err = GretlError.data;
                } else {
                    pprintf(prn, _("\nFound %d valid sheet(s)\n").ptr, xinfo.n_sheets);
                }
            }
        }
    }

    return err;
}

private void record_xlsx_params(XlsxInfo * xinfo, int * list)
{
    if (list != null && list[0] == 3) {
        list[1] = xinfo.selsheet + 1;
        list[2] = xinfo.xoffset;
        list[3] = xinfo.yoffset;
    }
}

/* When we're called from the command-line, the @list
   argument may contain a sheet number followed by
   x and y offsets. If a sheet number was given, it is
   already handled by this point: here we check and
   record the offsets, if present.
*/

private int set_xlsx_offsets_from_cli(XlsxInfo * xinfo,
                                      const(int) * list)
{
    int err = 0;

    if (list != null && list[0] == 3) {
        int xoff = list[2];
        int yoff = list[3];

        if (xoff < 0 || yoff < 0) {
            gretl_errmsg_set(_("Invalid argument for worksheet import").ptr);
            err = GretlError.data;
        } else {
            xinfo.xoffset = xoff;
            xinfo.yoffset = yoff;
        }
    }

    return err;
}

private void xlsx_dates_check(DatasetRaw * dset)
{
    int t, maybe_dates = 1;
    int date_min = 0, date_max = 0;
    int d, delta_min = 0, delta_max = 0;
    int dzero_count = 0;

    /* We're dealing here with the case where our prior heuristics
       suggest we got an "observations" column, yet the values we
       found there were numeric (and we converted them to strings).
       Here we see if it might be reasonable to interpret the
       labels as representing MS dates (days since Dec 31, 1899).

       For this purpose we'll require that all the obs labels are
       integer strings, and that the gap between successive values
       should be constant, or variable to a degree that's consistent
       with a sane time-series frequency.

       We should bear in mind, however, that the numeric values that
       we started with could be plain years rather than MS dates.
    */

    for (t=0; t<dset.n && maybe_dates; t++) {
        if (!integer_string(dset.S[t])) {
            maybe_dates = 0;
        } else if (t == 0) {
            if (!strcmp(dset.S[0], "1")) {
                maybe_dates = 0;
            } else {
                date_min = date_max = atoi(dset.S[t]);
            }
        } else {
            d = atoi(dset.S[t]);
            if (d < date_min) {
                date_min = d;
            }
            if (d > date_max) {
                date_max = d;
            }
            d = atoi(dset.S[t]) - atoi(dset.S[t-1]);
            if (d == 0) {
                if (++dzero_count < 20) {
                    fprintf(stderr, "dates? delta = 0 at obs %d\n", t+1);
                }
            }
            if (t == 1) {
                delta_min = delta_max = d;
            } else if (d < delta_min) {
                delta_min = d;
            } else if (d > delta_max) {
                delta_max = d;
            }
        }
    }

    if (maybe_dates && delta_max < 0) {
        /* allow for the possibility that time runs backwards */
        int tmp = delta_min;

        delta_min = -delta_max;
        delta_max = -tmp;
        fprintf(stderr, " xlsx_dates_check: diffmin=%d, diffmax=%d\n",
                delta_min, delta_max);
    }

    if (maybe_dates) {
        /* are these things in fact more plausibly years? */
        if (delta_min == 1 && delta_max == 1 &&
            date_min > 1749 && date_max < 2050) {
            maybe_dates = 0;
        }
    }

    if (maybe_dates) {
        if (delta_min >= 364 && delta_max <= 365) {
            ; /* annual? */
        } else if (delta_min >= 90 && delta_max <= 92) {
            ; /* quarterly? */
        } else if (delta_min >= 28 && delta_max <= 31) {
            ; /* monthly? */
        } else if (delta_min == 7 && delta_max == 7) {
            ; /* weekly? */
        } else if (delta_min == 1 && delta_max <= 7) {
            /* daily? delta_max of 7 is high, but could be
               valid for (e.g.) US financial market data that
               includes "9/11" (2001-09-11)
            */
            ;
        } else {
            /* unsupported frequency or nonsensical */
            maybe_dates = 0;
        }
    }

    if (maybe_dates) {
        char[OBSLEN] datestr = 0;

        for (t=0; t<dset.n; t++) {
            /* FIXME detect use of 1904-based dates? */
            MS_excel_date_string(datestr.ptr, atoi(dset.S[t]), 0, 0);
            strcpy(dset.S[t], datestr.ptr);
        }
    }
}

private int finalize_xlsx_import(DatasetRaw * dset,
                                 XlsxInfo * xinfo,
                                 const(char) * fname,
                                 Gretlopt opt,
                                 Prn * prn)
{
    int merge = (dset.Z != null);
    int err;

    err = import_prune_columns(xinfo.dset);

    if (!err) {
        int i;

        for (i=1; i<xinfo.dset.v && !err; i++) {
            if (*xinfo.dset.varname[i] == '\0') {
                pprintf(prn, _("Name missing for variable %d\n").ptr, i);
                err = GretlError.data;
            }
        }
    }

    if (!err && fix_varname_duplicates(xinfo.dset)) {
        pputs(prn, _("warning: some variable names were duplicated\n").ptr);
    }

    if (!err) {
        if (xinfo.trydates) {
            xlsx_dates_check(xinfo.dset);
        }
        if (xinfo.dset.S != null) {
            int mpd = import_ts_check(xinfo.dset);

            if (mpd < 0) {
                err = E_OBSCOL;
            }
        }
    }
    if (!err) {
        err = merge_or_replace_data(dset, &xinfo.dset,
                                    get_merge_opts(opt), prn);
    }
    if (!err && !merge) {
        dataset_add_import_info(dset, fname, GretlFileType.xlsx);
    }

    return err;
}

/* Public driver function for retrieving data from OOXML workbook
   file.

   If we're coming from the command line @sheetname may be non-NULL
   and/or @list may contain specification of sheet number (1-based),
   row offset and/or column offset.

   If we're coming from the GUI @sheetname will be NULL but @list will
   be non-NULL, and is used to record choices made here via the
   (here dropped) sheet dialog.
*/

extern (C) int xlsx_get_data(const(char) * fname, int * list, char * sheetname,
                             DatasetRaw * dset, Gretlopt opt, Prn * prn,
                             void * parent)
{
    int gui = (opt & Gretlopt.OPT_G) != 0;
    char * save_errmsg = null;
    XlsxInfo xinfo;
    char[32] dname = 0;
    int err;

    err = open_import_zipfile(fname, dname.ptr, prn);
    if (err) {
        return err;
    }

    xlsx_info_init(&xinfo, opt);
    err = xlsx_gather_sheet_names(&xinfo, sheetname, list, prn);

    if (!err) {
        /* C: if (gui) the sheet dialog was run here; without it, both cases use the offsets given
           in @list */
        err = set_xlsx_offsets_from_cli(&xinfo, list);
    }

 try_again:

    if (!err) {
        err = xlsx_read_worksheet(&xinfo, fname, prn);
    }

    if (!err && xinfo.dset != null) {
        err = finalize_xlsx_import(dset, &xinfo, fname, opt, prn);
        if (err == E_OBSCOL) {
            opt |= Gretlopt.OPT_A;
            xlsx_info_reset(&xinfo, opt);
            err = 0;
            goto try_again;
        } else if (!err && gui) {
            record_xlsx_params(&xinfo, list);
        }
    }

    gretl_print_flush_stream(prn);
    xlsx_info_free(&xinfo);

    if (err > 0) {
        /* preserve any error message so it's not overwritten
           by remove_temp_dir();
        */
        save_errmsg = maybe_save_gretl_errmsg(err);
    }

    remove_temp_dir(dname.ptr);

    if (save_errmsg != null) {
        gretl_errmsg_set(save_errmsg);
        free(save_errmsg);
    } else if (err > 0) {
        pputs(prn, _("Error reading xlsx data\n").ptr);
    }

    return err;
}
