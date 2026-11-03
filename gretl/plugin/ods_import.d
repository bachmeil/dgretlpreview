/* Port of plugin/ods_import.c: importer for OpenDocument spreadsheets (.ods).

   In the C this is a gretl plugin; here it is an ordinary module that plugins.d returns for
   "ods_get_data". The GTK sheet dialog (ods_sheet_dialog, ods_book_init, ods_min_offset) is dropped:
   ods_get_data() keeps its C signature (the trailing GtkWidget * is void *) and the OPT_G branch uses
   the sheet and offsets given in @sheetname and @list, like the command-line branch. The unconditional
   diagnostics that the C writes to stderr (table dimensions, top-left analysis) are kept. */

module gretl.plugin.ods_import;


import core.stdc.stdio : FILE, fclose, fopen, fputc, fprintf, fread, sprintf, stderr;
import core.stdc.stdlib : atof, atoi, free, malloc, realloc;
import core.stdc.string : strcmp, strcpy, strncat;

import gretl.csvdatah : NON_NUMERIC;
import gretl.csvdata : import_na_string, import_obs_label, non_numeric_check;
import gretl.dataio : dataset_add_import_info, get_merge_opts, merge_or_replace_data;
import gretl.dataioh : GretlFileType;
import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_errors : gretl_errmsg_get, gretl_errmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : in_gretl_list;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_string_table : GretlStringTable, gretl_string_table_destroy, gretl_string_table_finalize,
    gretl_string_table_index, gretl_string_table_validate;
import gretl.gretl_utils : printlist;
import gretl.gretl_xml : gretl_xml_open_doc_root;
import gretl.plugin.import_common;
import gretl.plugin.importerh;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.strutils : gretl_utf8_strncat_trim, integer_string;
import gretl.xmlh;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;

private alias XUC = const(xmlChar) *;

/* unported dependencies: per-file declarations from the 2026c headers */
private enum OdsType {
    none,
    numeric,
    date,
    time,
    bool_,
    string_
}

/* Is this (always) right? */
private bool ods_string_type(int t)
{
    return t == OdsType.string_ || t == OdsType.date || t == OdsType.time;
}

private enum XOFF_UNDEF = 999999;

private struct OdsTable {
    char * name;
    xmlNodePtr node;
    int rows;        /* total rows defined */
    int cols;        /* total columns */
    int xoffset;     /* offset to first non-blank column */
    int yoffset;     /* offset to first non-blank row */
    int empty;       /* has any content (0) or not (1) */
}

private struct OdsSheet {
    BookFlag flags;
    xmlDocPtr doc;          /* document pointer */
    int n_tables;           /* number of tables */
    OdsTable ** tables;     /* pointers to table info */
    int seltab;             /* number of selected table */
    int xoffset;            /* col offset chosen by user */
    int yoffset;            /* row offset chosen by user */
    DatasetRaw * dset;         /* dataset struct */
    int * codelist;         /* list of string-valued variables */
    GretlStringTable * st;  /* table for the above */
}

private OdsTable * ods_table_new(xmlNodePtr node, int * err)
{
    OdsTable * tab = null;
    char * name;

    name = cast(char *) xmlGetProp(node, cast(XUC) "name");
    if (name == null) {
        *err = GretlError.data;
        return null;
    }

    tab = cast(OdsTable *) malloc(OdsTable.sizeof);

    if (tab != null) {
        tab.name = name;
        tab.node = node;
        tab.rows = tab.cols = 0;
        tab.xoffset = XOFF_UNDEF;
        tab.yoffset = 0;
        tab.empty = 1;
    } else {
        *err = GretlError.alloc;
        free(name);
    }

    return tab;
}

private int ods_sheet_add_table(OdsSheet * sheet, OdsTable * tab)
{
    int n = sheet.n_tables;
    OdsTable ** tabs = null;

    tabs = cast(OdsTable **) realloc(sheet.tables, (n+1) * (*sheet.tables).sizeof);
    if (tabs == null) {
        return GretlError.alloc;
    }

    sheet.tables = tabs;
    sheet.tables[n] = tab;
    sheet.n_tables += 1;

    return 0;
}

private OdsSheet * ods_sheet_new(xmlDocPtr doc, int * err)
{
    OdsSheet * sheet;

    sheet = cast(OdsSheet *) malloc(OdsSheet.sizeof);

    if (sheet != null) {
        sheet.flags = cast(BookFlag) 0;
        sheet.doc = doc;
        sheet.n_tables = 0;
        sheet.tables = null;
        sheet.seltab = -1;
        sheet.xoffset = 0;
        sheet.yoffset = 0;
        sheet.dset = null;
        sheet.codelist = null;
        sheet.st = null;
    } else {
        *err = GretlError.alloc;
    }

    return sheet;
}

private void ods_table_free(OdsTable * tab)
{
    free(tab.name);
    free(tab);
}

private void ods_sheet_free(OdsSheet * sheet)
{
    if (sheet != null) {
        int i;

        for (i=0; i<sheet.n_tables; i++) {
            ods_table_free(sheet.tables[i]);
        }
        free(sheet.tables);

        if (sheet.doc != null) {
            xmlFreeDoc(sheet.doc);
        }

        destroy_dataset(sheet.dset);
        free(sheet.codelist);
        if (sheet.st != null) {
            gretl_string_table_destroy(sheet.st);
        }

        free(sheet);
    }
}

private void ods_table_print(OdsTable * tab)
{
    fprintf(stderr, "Table \"%s\": ", tab.name);

    if (tab.empty) {
        fprintf(stderr, "(empty)\n");
    } else {
        fprintf(stderr, "%d x %d, xoff = %d, yoff = %d, data area %d x %d\n",
                tab.rows, tab.cols, tab.xoffset, tab.yoffset,
                tab.rows - tab.yoffset, tab.cols - tab.xoffset);
    }
}

private void ods_sheet_print(OdsSheet * sheet)
{
    if (sheet != null) {
        int i;

        fprintf(stderr, "Sheet: %d tables\n", sheet.n_tables);

        for (i=0; i<sheet.n_tables; i++) {
            ods_table_print(sheet.tables[i]);
        }
    }
}

private int ods_sheet_prune(OdsSheet * sheet, Prn * prn)
{
    int err = 0;

    if (sheet == null) {
        err = GretlError.data;
    } else {
        int i, j;

        for (i=0; i<sheet.n_tables; i++) {
            if (sheet.tables[i].empty) {
                ods_table_free(sheet.tables[i]);
                sheet.n_tables -= 1;
                for (j=i; j<sheet.n_tables; j++) {
                    sheet.tables[j] = sheet.tables[j+1];
                }
                i--;
            }
        }

        if (sheet.n_tables == 0) {
            pputs(prn, _("File contains no data").ptr);
            err = GretlError.data;
        }
    }

    return err;
}

private int get_ods_value_type(xmlNodePtr node)
{
    int ret = OdsType.none;
    char * s;

    s = cast(char *) xmlGetProp(node, cast(XUC) "value-type");
    if (s == null) {
        return ret;
    }

    if (!strcmp(s, "float") ||
        !strcmp(s, "percentage") ||
        !strcmp(s, "currency")) {
        ret = OdsType.numeric;
    } else if (!strcmp(s, "date")) {
        ret = OdsType.date;
    } else if (!strcmp(s, "time")) {
        ret = OdsType.time;
    } else if (!strcmp(s, "boolean")) {
        ret = OdsType.bool_;
    } else if (!strcmp(s, "string")) {
        ret = OdsType.string_;
    }

    free(s);

    return ret;
}

private int ods_non_numeric_check(OdsSheet * sheet, Prn * prn)
{
    GretlStringTable * st = null;
    int * nlist = null;
    int err = 0;

    err = non_numeric_check(sheet.dset, &nlist, &st, prn);

    if (!err) {
        sheet.codelist = nlist;
        sheet.st = st;
    }

    return err;
}

private int p_content_NA(xmlNodePtr cur)
{
    char * s = null;
    int ret = 0;

    cur = cur.children;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "p")) {
            s = cast(char *) xmlNodeGetContent(cur);
            break;
        }
        cur = cur.next;
    }

    if (s != null) {
        if (!strcmp(s, "#N/A")) {
            ret = 1;
        }
        free(s);
    }

    return ret;
}

private char * get_ods_string_value(xmlNodePtr cur)
{
    char * sval;

    sval = cast(char *) xmlGetProp(cur, cast(XUC) "string-value");
    if (sval != null) {
        return sval;
    }

    cur = cur.children;

    while (cur != null) {
        if (!xmlStrcmp(cur.name, cast(XUC) "p")) {
            sval = cast(char *) xmlNodeGetContent(cur);
            break;
        }
        cur = cur.next;
    }

    return sval;
}

private int get_ods_bool_value(xmlNodePtr cur)
{
    char * tmp;
    int ret = 0;

    tmp = cast(char *) xmlGetProp(cur, cast(XUC) "boolean-value");
    if (tmp != null) {
        ret = (strcmp(tmp, "true") == 0);
        free(tmp);
    }

    return ret;
}

private double get_ods_numeric_value(xmlNodePtr cur)
{
    char * tmp;
    double ret = NADBL;

    tmp = cast(char *) xmlGetProp(cur, cast(XUC) "value");
    if (tmp != null) {
        if (!strcmp(tmp, "0") && p_content_NA(cur)) {
            ret = NADBL;
        } else {
            ret = atof(tmp);
        }
        free(tmp);
    }

    return ret;
}

private int ods_cell_has_content(xmlNodePtr node)
{
    return (get_ods_value_type(node) != OdsType.none);
}

private int ods_row_height(xmlNodePtr p)
{
    char * s;
    int h = 1;

    s = cast(char *) xmlGetProp(p, cast(XUC) "number-rows-repeated");
    if (s != null) {
        if (*s != '\0') {
            h = atoi(s);
        }
        free(s);
    }

    return h;
}

private int ods_cell_width(xmlNodePtr p)
{
    char * s;
    int w = 1;

    s = cast(char *) xmlGetProp(p, cast(XUC) "number-columns-repeated");
    if (s != null) {
        if (*s != '\0') {
            w = atoi(s);
        }
        free(s);
    }

    return w;
}

private const(char) * ods_name(int t)
{
    if (t == OdsType.none)
        return "blank";
    if (t == OdsType.numeric)
        return "numeric";
    if (t == OdsType.date)
        return "date";
    if (t == OdsType.time)
        return "time";
    if (t == OdsType.bool_)
        return "boolean";
    if (t == OdsType.string_)
        return "string";

    return "blank";
}

private int ods_error(OdsSheet * sheet,
                      int i, int j,
                      int etype, int vtype,
                      Prn * prn)
{
    int si = i + sheet.yoffset + 1;
    int sj = j + sheet.xoffset + 1;

    pprintf(prn, _("Sheet row %d, column %d").ptr, si, sj);

    if ((sheet.flags & BookFlag.auto_varnames) || i == 0) {
        pputs(prn, ":\n");
    } else {
        int v = (sheet.flags & BookFlag.obs_labels)? j : j + 1;

        if (v > 0 && v < sheet.dset.v) {
            pprintf(prn, " (\"%s\"):\n", sheet.dset.varname[v]);
        } else {
            pputs(prn, ":\n");
        }
    }

    pprintf(prn, _("expected %s but found %s").ptr,
            ods_name(etype), ods_name(vtype));
    pputc(prn, '\n');

    return GretlError.data;
}

private int ods_handle_stringval(OdsSheet * sheet,
                                 int i, int t,
                                 const(char) * s,
                                 int nr,
                                 Prn * prn)
{
    int err = 0;

    if (sheet.dset.Z[i][t] == NON_NUMERIC) {
        int j, vj, ix;

        for (j=0, vj=i; j<nr && vj<sheet.dset.v; j++, vj++) {
            ix = gretl_string_table_index(sheet.st, s, vj, 0, prn);
            if (ix > 0) {
                sheet.dset.Z[vj][t] = cast(double) ix;
            } else {
                err = GretlError.data;
                break;
            }
        }
    }

    return err;
}

private int real_read_cell(xmlNodePtr cur,
                           OdsSheet * sheet, int pass,
                           int iread, int * preadcol,
                           Prn * prn)
{
    int verbose = 0;
    char * val = null;
    int jread = *preadcol;
    int obscol = (sheet.flags & BookFlag.obs_labels)? 1 : 0;
    int blank0 = (sheet.flags & BookFlag.obs_blank)? 1 : 0;
    int vnames = (sheet.flags & BookFlag.auto_varnames)? 0 : 1;
    int nr, j, v, vj, t, vtype;
    double x = NADBL;
    int err = 0;

    v = jread + 1 - obscol;
    t = iread - vnames;

    if (v >= sheet.dset.v || t >= sheet.dset.n) {
        fprintf(stderr, "v = %d, t = %d: out of bounds?\n", v, t);
        return GretlError.data;
    }

    vtype = get_ods_value_type(cur);
    nr = ods_cell_width(cur);

    *preadcol += nr;

    if (pass == 2) {
        /* just going after string-valued variables on second pass */
        if (iread == 0 && vnames) {
            return 0;
        } else if (jread == 0 && obscol) {
            return 0;
        } else if (ods_string_type(vtype) && in_gretl_list(sheet.codelist, v)) {
            val = get_ods_string_value(cur);
            err = ods_handle_stringval(sheet, v, t, val, nr, prn);
            free(val);
        }
        return err;
    }

    if (iread == 0 && vnames) {
        /* reading a variable name */
        jread += blank0;
        v += blank0;
        if (jread == 0 && obscol) {
            return 0;
        }
        if (vtype == OdsType.string_) {
            val = get_ods_string_value(cur);
            if (val != null) {
                *sheet.dset.varname[v] = '\0';
                strncat(sheet.dset.varname[v], val, VNAMELEN - 1);
                err = check_imported_varname(sheet.dset.varname[v],
                                             v, iread, jread, prn);
                free(val);
            } else {
                err = ods_error(sheet, iread, jread, OdsType.string_,
                                OdsType.none, prn);
            }
        } else if (vtype != OdsType.none) {
            err = ods_error(sheet, iread, jread, OdsType.string_,
                            vtype, prn);
        }
        return err;
    }

    if (jread == 0 && obscol) {
        /* reading an observation label */
        if (vtype == OdsType.string_) {
            val = get_ods_string_value(cur);
            if (val != null) {
                if (verbose) {
                    fprintf(stderr, " obs string: '%s'\n", val);
                }
            } else {
                err = ods_error(sheet, iread, jread, OdsType.string_,
                                OdsType.none, prn);
            }
        } else if (vtype == OdsType.date) {
            val = cast(char *) xmlGetProp(cur, cast(XUC) "date-value");
            if (val != null) {
                if (verbose) {
                    fprintf(stderr, " date: '%s'\n", val);
                }
            } else {
                err = ods_error(sheet, iread, jread, OdsType.date,
                                OdsType.none, prn);
            }
        } else if (vtype == OdsType.numeric) {
            val = cast(char *) xmlGetProp(cur, cast(XUC) "value");
            if (val != null) {
                if (verbose) {
                    fprintf(stderr, " numeric obs: '%s'\n", val);
                }
            } else {
                err = ods_error(sheet, iread, jread, OdsType.numeric,
                                OdsType.none, prn);
            }
        } else {
            err = ods_error(sheet, iread, jread, OdsType.date, vtype, prn);
        }
        if (!err) {
            gretl_utf8_strncat_trim(sheet.dset.S[t], val, OBSLEN - 1);
        }

        free(val);

        return err;
    }

    /* reading actual data */

    if (vtype == OdsType.numeric) {
        x = get_ods_numeric_value(cur);
    } else if (vtype == OdsType.bool_) {
        x = get_ods_bool_value(cur);
    } else if (vtype == OdsType.none) {
        ;
    } else if (ods_string_type(vtype)) {
        val = get_ods_string_value(cur);
        if (val != null && import_na_string(val)) {
            ;
        } else if (val != null && *val != '\0') {
            x = NON_NUMERIC;
            sheet.flags |= BookFlag.non_numeric;
        }
        free(val);
    } else {
        fprintf(stderr, " unrecognized value-type %d\n", vtype);
        err = GretlError.data;
    }

    if (err) {
        ods_error(sheet, iread, jread, OdsType.numeric, vtype, prn);
    } else {
        for (j=0, vj=v; j<nr && vj<sheet.dset.v; j++, vj++) {
            sheet.dset.Z[vj][t] = x;
        }
    }

    return err;
}

private int read_data_row(xmlNodePtr cur,
                          OdsTable * tab,
                          OdsSheet * sheet,
                          int pass,
                          int readrow,
                          Prn * prn)
{
    int readmax = tab.cols - sheet.xoffset;
    int tabcol = 0, readcol = 0;
    int err = 0;

    cur = cur.children;

    while (cur != null && !err && readcol < readmax) {
        if (!xmlStrcmp(cur.name, cast(XUC) "table-cell")) {
            if (tabcol >= sheet.xoffset) {
                err = real_read_cell(cur, sheet, pass,
                                     readrow, &readcol,
                                     prn);
            }
            tabcol += ods_cell_width(cur);
        }
        cur = cur.next;
    }

    return err;
}

private int sheet_allocate_data(OdsSheet * sheet,
                                OdsTable * tab)
{
    int n = tab.rows - sheet.yoffset - 1;
    int v = tab.cols - sheet.xoffset + 1;
    int i, labels = 0;

    if (sheet.flags & BookFlag.auto_varnames) {
        n++;
    }

    if (sheet.flags & BookFlag.obs_labels) {
        labels = 1;
        v--;
    }

    if (n <= 0 || v <= 1) {
        return GretlError.data;
    }

    fprintf(stderr, "sheet_allocate_data: n=%d, v=%d\n",
            n, v);

    sheet.dset = create_new_dataset(v, n, labels);
    if (sheet.dset == null) {
        return GretlError.alloc;
    }

    /* write fallback variable names */
    for (i=1; i<v; i++) {
        sprintf(sheet.dset.varname[i], "v%d", i);
    }

    return 0;
}

/* Look at the cells in the top left-hand corner of the reading area
   of the table: try to determine (a) if we have an observations
   column and (b) if we have a varnames row.
*/

private int analyse_top_left(OdsSheet * sheet, OdsTable * tab)
{
    xmlNodePtr colp;
    xmlNodePtr rowp = tab.node.children;
    xmlNodePtr p00 = null;
    xmlNodePtr p01 = null;
    xmlNodePtr p10 = null;
    int readcols = tab.cols - sheet.xoffset;
    int nr, tabcol, tabrow = 0;
    int done = 0;
    int err = 0;

    fprintf(stderr, "analyse_top_left: sheet->xoffset = %d, readcols = %d\n",
            sheet.xoffset, readcols);

    while (!err && rowp != null && !done) {
        if (!xmlStrcmp(rowp.name, cast(XUC) "table-row")) {
            nr = ods_row_height(rowp);
            if (tabrow == sheet.yoffset) {
                colp = rowp.children;
                tabcol = 0;
                while (colp != null && !err) {
                    if (!xmlStrcmp(colp.name, cast(XUC) "table-cell")) {
                        if (tabcol == sheet.xoffset) {
                            p00 = colp;
                        } else if (tabcol == sheet.xoffset + 1) {
                            p01 = colp;
                        }
                        tabcol += ods_cell_width(colp);
                    }
                    colp = colp.next;
                }
            } else if (tabrow == sheet.yoffset + 1) {
                colp = rowp.children;
                tabcol = 0;
                while (colp != null && !err) {
                    if (!xmlStrcmp(colp.name, cast(XUC) "table-cell")) {
                        if (tabcol == sheet.xoffset) {
                            p10 = colp;
                        }
                        tabcol += ods_cell_width(colp);
                    }
                    colp = colp.next;
                }
            }
            tabrow += nr;
        }
        if (readcols > 1) {
            done = (p01 != null && p10 != null);
        } else {
            done = (p10 != null);
        }
        rowp = rowp.next;
    }

    if (!done) {
        fprintf(stderr, "analyse_top_left: failed\n");
        err = GretlError.data;
    } else {
        int vt10, vt00 = OdsType.none, vt01 = OdsType.none;

        if (p00 != null) {
            vt00 = get_ods_value_type(p00);
            fprintf(stderr, "cell(0,0): type = %s\n", ods_name(vt00));
        } else {
            fprintf(stderr, "cell(0,0): blank\n");
            sheet.flags |= BookFlag.obs_blank;
        }

        vt10 = get_ods_value_type(p10);
        fprintf(stderr, "cell(1,0): type = %s\n", ods_name(vt10));

        if (p01 == null) {
            /* single column */
            if (vt00 != OdsType.string_) {
                sheet.flags |= BookFlag.auto_varnames;
            }
        } else {
            vt01 = get_ods_value_type(p01);
            fprintf(stderr, "cell(0,1): type = %s\n", ods_name(vt01));
            if (vt01 != OdsType.string_) {
                sheet.flags |= BookFlag.auto_varnames;
            }
            if (vt00 == OdsType.none) {
                sheet.flags |= BookFlag.obs_labels;
            } else if (vt00 == OdsType.string_) {
                char * val = get_ods_string_value(p00);

                fprintf(stderr, "cell(0,0): val = '%s'\n", val);
                if (import_obs_label(val)) {
                    fprintf(stderr, "looks like obs label\n");
                    sheet.flags |= BookFlag.obs_labels;
                }
                free(val);
            }
        }
    }

    fprintf(stderr, "analyse_top_left: vnames=%d, obscol=%d, returning %d\n",
            (sheet.flags & BookFlag.auto_varnames)? 0 : 1,
            (sheet.flags & BookFlag.obs_labels)? 1 : 0, err);

    return err;
}

private int check_repeat_data_row(OdsSheet * sheet,
                                  int iread,
                                  Prn * prn)
{
    int vnames = (sheet.flags & BookFlag.auto_varnames)? 0 : 1;
    int i, t = iread - vnames;

    if (t < 1 || t >= sheet.dset.n) {
        pprintf(prn, _("Found a repeated row in the wrong place\n").ptr);
        return GretlError.data;
    }

    for (i=1; i<sheet.dset.v; i++) {
        sheet.dset.Z[i][t] = sheet.dset.Z[i][t-1];
    }

    if (sheet.dset.S != null) {
        strcpy(sheet.dset.S[t], sheet.dset.S[t-1]);
    }

    return 0;
}

private int read_table_content(OdsSheet * sheet, Prn * prn)
{
    OdsTable * tab;
    xmlNodePtr top;
    xmlNodePtr cur;
    int i, nr, maxrow;
    int tabrow = 0, readrow = 0;
    int pass = 1;
    int err = 0;

    if (sheet.seltab < 0 || sheet.seltab >= sheet.n_tables) {
        return GretlError.data;
    }

    tab = sheet.tables[sheet.seltab];
    err = analyse_top_left(sheet, tab);

    if (!err) {
        err = sheet_allocate_data(sheet, tab);
    }

    if (err) {
        return err;
    }

    maxrow = tab.rows - sheet.yoffset;
    top = cur = tab.node.children;

    gretl_push_c_numeric_locale();

 tryagain:

    while (cur != null && !err && readrow < maxrow) {
        if (!xmlStrcmp(cur.name, cast(XUC) "table-row")) {
            nr = ods_row_height(cur);
            if (tabrow >= sheet.yoffset) {
                err = read_data_row(cur, tab, sheet, pass, readrow++, prn);
                for (i=1; i<nr && !err; i++) {
                    err = check_repeat_data_row(sheet, readrow++, prn);
                }
            }
            tabrow += nr;
        }
        cur = cur.next;
    }

    if (pass == 1 && (sheet.flags & BookFlag.non_numeric)) {
        err = ods_non_numeric_check(sheet, prn);
        if (sheet.codelist != null) {
            tabrow = readrow = 0;
            cur = top;
            pass = 2;
            goto tryagain;
        }
    }

    gretl_pop_c_numeric_locale();

    return err;
}

private int get_table_dimensions(xmlNodePtr cur, OdsSheet * sheet)
{
    OdsTable * tab = null;
    xmlNodePtr rowp;
    int hascont, nr, nc, row_empty;
    int cols, xoffset, xtrail;
    int rows, rchk;
    int err = 0;

    tab = ods_table_new(cur, &err);
    if (tab == null) {
        return err;
    }

    err = ods_sheet_add_table(sheet, tab);
    if (err) {
        return err;
    }

    cur = cur.children;
    rows = rchk = 0;

    while (cur != null && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "table-row")) {
            nr = ods_row_height(cur);
            cols = xoffset = xtrail = 0;
            row_empty = 1;
            rowp = cur.children;
            while (rowp != null && !err) {
                if (!xmlStrcmp(rowp.name, cast(XUC) "table-cell")) {
                    hascont = ods_cell_has_content(rowp);
                    nc = ods_cell_width(rowp);
                    if (hascont) {
                        row_empty = 0;
                        tab.empty = 0;
                        xtrail = 0;
                    } else if (row_empty) {
                        xoffset += nc;
                    }
                    if (rowp.next == null) {
                        /* last cell(s) in row: ignore if blank */
                        cols += (hascont)? nc : 0;
                    } else {
                        cols += nc;
                        if (!hascont) {
                            xtrail += nc;
                        }
                    }
                }
                rowp = rowp.next;
            }
            cols -= xtrail;
            if (!err) {
                rows += nr;
                if (!row_empty) {
                    rchk = rows;
                }
                if (cols > tab.cols) {
                    tab.cols = cols;
                }
                if (tab.empty) {
                    tab.yoffset += nr;
                }
                if (!row_empty && xoffset < tab.xoffset) {
                    tab.xoffset = xoffset;
                }
            }
        }
        cur = cur.next;
    }

    tab.rows = rchk;

    return err;
}

private OdsSheet * ods_read_content(Prn * prn, int * err)
{
    OdsSheet * sheet = null;
    xmlDocPtr doc = null;
    xmlNodePtr cur = null;
    xmlNodePtr c1;
    xmlNodePtr c2;

    *err = gretl_xml_open_doc_root("content.xml",
                                   "document-content",
                                   &doc, &cur);

    if (*err) {
        pprintf(prn, _("didn't get office:document-content\n").ptr);
        pprintf(prn, "%s", gretl_errmsg_get());
        return null;
    }

    sheet = ods_sheet_new(doc, err);
    if (sheet == null) {
        xmlFreeDoc(doc);
        return null;
    }

    /* Now walk the tree */
    cur = cur.children;
    while (cur != null && !*err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "body")) {
            c1 = cur.children;
            while (c1 != null && !*err) {
                if (!xmlStrcmp(c1.name, cast(XUC) "spreadsheet")) {
                    c2 = c1.children;
                    while (c2 != null && !*err) {
                        if (!xmlStrcmp(c2.name, cast(XUC) "table")) {
                            *err = get_table_dimensions(c2, sheet);
                        }
                        c2 = c2.next;
                    }
                }
                c1 = c1.next;
            }
        }
        cur = cur.next;
    }

    return sheet;
}

private int check_mimetype(Prn * prn)
{
    const(char) * odsmime =
        "application/vnd.oasis.opendocument.spreadsheet";
    char[48] buf = 0;
    FILE * fp;
    int err = 0;

    fp = fopen("mimetype", "r");
    if (fp == null) {
        pprintf(prn, _("Couldn't find mimetype\n").ptr);
        err = GretlError.fopen;
    } else {
        if (fread(buf.ptr, 1, 46, fp) != 46 ||
            strcmp(buf.ptr, odsmime)) {
            pprintf(prn, _("Wrong or missing mime type,\n should be '%s'\n").ptr,
                    odsmime);
            err = GretlError.data;
        }
        fclose(fp);
    }

    return err;
}

private void record_ods_params(OdsSheet * sheet, int * list)
{
    if (list != null && list[0] == 3) {
        list[1] = sheet.seltab + 1;
        list[2] = sheet.xoffset;
        list[3] = sheet.yoffset;
    }
}

private int set_ods_params_from_cli(OdsSheet * sheet,
                                    const(int) * list,
                                    char * sheetname)
{
    int gotname = (sheetname != null && *sheetname != '\0');
    int gotlist = (list != null && list[0] == 3);
    int i;

    sheet.seltab = 0; /* default to first */

    if (!gotname && !gotlist) {
        /* no user specs */
        sheet.xoffset = sheet.tables[0].xoffset;
        sheet.yoffset = sheet.tables[0].yoffset;
        return 0;
    }

    if (gotname || (gotlist && list[1] > 0)) {
        /* invalidate this pro tem */
        sheet.seltab = -1;
    }

    if (gotname) {
        for (i=0; i<sheet.n_tables; i++) {
            if (!strcmp(sheetname, sheet.tables[i].name)) {
                sheet.seltab = i;
                break;
            }
        }
        if (sheet.seltab < 0 && integer_string(sheetname)) {
            i = atoi(sheetname);
            if (i >= 1 && i <= sheet.n_tables) {
                sheet.seltab = i - 1;
            }
        }
    }

    if (gotlist) {
        if (!gotname && list[1] > 0) {
            /* convert to zero-based */
            sheet.seltab = list[1] - 1;
        }
        sheet.xoffset = list[2];
        sheet.yoffset = list[3];
    }

    if (sheet.seltab < 0 || sheet.seltab >= sheet.n_tables ||
        sheet.xoffset < 0 || sheet.yoffset < 0) {
        gretl_errmsg_set(_("Invalid argument for worksheet import").ptr);
        fprintf(stderr, "seltab=%d, xoffset=%d, yoffset=%d\n",
                sheet.seltab, sheet.xoffset, sheet.yoffset);
        return GretlError.data;
    }

    return 0;
}

private int finalize_ods_import(DatasetRaw * dset,
                                OdsSheet * sheet,
                                const(char) * fname,
                                Gretlopt opt,
                                Prn * prn)
{
    int err = import_prune_columns(sheet.dset);
    int merge = (dset.Z != null);

    if (!err && sheet.dset.S != null) {
        import_ts_check(sheet.dset);
    }

    if (!err) {
        err = merge_or_replace_data(dset, &sheet.dset,
                                    get_merge_opts(opt), prn);
    }

    if (!err && !merge) {
        dataset_add_import_info(dset, fname, GretlFileType.ods);
    }

    return err;
}

extern (C) int ods_get_data(const(char) * fname, int * list, char * sheetname,
                            DatasetRaw * dset, Gretlopt opt, Prn * prn,
                            void * parent)
{
    int gui = (opt & Gretlopt.OPT_G) != 0;
    OdsSheet * sheet = null;
    char[32] dname = 0;
    int err;

    err = open_import_zipfile(fname, dname.ptr, prn);
    if (err) {
        return err;
    }

    if (!err) {
        err = check_mimetype(prn);
    }

    if (!err) {
        sheet = ods_read_content(prn, &err);
    }

    remove_temp_dir(dname.ptr);

    ods_sheet_print(sheet);

    if (!err) {
        err = ods_sheet_prune(sheet, prn);
    }

    printlist(list, "ods list");
    fprintf(stderr, "sheetname='%s'\n", sheetname);

    if (!err) {
        /* C: if (gui) the sheet dialog was run here; without it both cases use @list and @sheetname */
        err = set_ods_params_from_cli(sheet, list, sheetname);
    }

    if (!err) {
        err = read_table_content(sheet, prn);
    }

    if (!err && sheet.st != null) {
        err = gretl_string_table_validate(sheet.st, Gretlopt.OPT_S);
        if (err) {
            pputs(prn, _("Failed to interpret the data as numeric\n").ptr);
        } else {
            gretl_string_table_finalize(sheet.st, sheet.dset);
        }
    }

    if (!err) {
        err = finalize_ods_import(dset, sheet, fname, opt, prn);
        if (!err && gui) {
            record_ods_params(sheet, list);
        }
    }

    ods_sheet_free(sheet);

    return err;
}
