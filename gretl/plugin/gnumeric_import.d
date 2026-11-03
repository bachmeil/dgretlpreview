/* Port of plugin/gnumeric_import.c: importer for Gnumeric workbooks (.gnumeric: gzipped XML,
   which libxml2 reads directly).

   In the C this is a gretl plugin; here it is an ordinary module that plugins.d returns for
   "gnumeric_get_data". As with xlsx_import.d the GTK sheet dialog is dropped: gnumeric_get_data() keeps its
   C signature (the trailing GtkWidget * is void *), and the OPT_G branch uses the sheet and offsets
   given in @sheetname and @list, like the command-line branch. */

module gretl.plugin.gnumeric_import;


import core.stdc.stdio : fprintf, sprintf, stderr;
import core.stdc.stdlib : atof, atoi, free, realloc;
import core.stdc.string : strchr, strcmp, strcpy, strncat, strstr;

import gretl.calendar : MS_excel_date_string;
import gretl.csvdatah : NON_NUMERIC;
import gretl.csvdata : import_na_string, import_obs_label, non_numeric_check;
import gretl.dataio : get_merge_opts, merge_or_replace_data, dataset_add_import_info;
import gretl.dataioh : GretlFileType;
import gretl.dataset : create_new_dataset, destroy_dataset;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_list : in_gretl_list;
import gretl.gretl_prn : Prn, pprintf, pputs;
import gretl.gretl_string_table : GretlStringTable, gretl_string_table_destroy, gretl_string_table_finalize,
    gretl_string_table_index, gretl_string_table_validate;
import gretl.gretl_xml : gretl_xml_open_doc_root;
import gretl.plugin.import_common;
import gretl.plugin.importerh;
import gretl.libgretlh;
import gretl.missingh : NADBL;
import gretl.strutils : gretl_strdup, gretl_utf8_strncat_trim, numeric_string, string_is_blank, tailstrip;
import gretl.xmlh;
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;

private alias XUC = const(xmlChar) *;

/* unported dependencies: per-file declarations from the 2026c headers */
/* from gnumeric's value.h */
private enum ValueType {
    empty     = 10,
    boolean   = 20,
    integer   = 30,
    float_    = 40,
    error     = 50,
    string_   = 60,
    cellrange = 70,
    array     = 80
}

private struct WSheet {
    BookFlag flags;
    Gretlopt opt;
    int maxcol, maxrow;
    int col_offset;
    int row_offset;
    int colheads;
    int ID;
    char * name;
    DatasetRaw * dset;
    int * codelist;
    GretlStringTable * st;
}

private void wsheet_init(WSheet * sheet, Gretlopt opt)
{
    sheet.flags = cast(BookFlag) 0;
    sheet.opt = opt;
    sheet.col_offset = sheet.row_offset = 0;
    sheet.maxcol = sheet.maxrow = 0;
    sheet.colheads = 0;
    sheet.ID = 0;
    sheet.name = null;
    sheet.dset = null;
    sheet.codelist = null;
    sheet.st = null;
}

private void wsheet_free(WSheet * sheet)
{
    destroy_dataset(sheet.dset);
    free(sheet.name);
    free(sheet.codelist);
    if (sheet.st != null) {
        gretl_string_table_destroy(sheet.st);
    }
    wsheet_init(sheet, Gretlopt.OPT_NONE);
}

private bool vtype_is_numeric(int v)
{
    return v == ValueType.boolean || v == ValueType.integer || v == ValueType.float_;
}

private int wsheet_allocate(WSheet * sheet, int cols, int rows)
{
    int labels = (sheet.flags & BookFlag.obs_labels)? 1 : 0;
    int nvars = cols + 1 - labels;
    int nobs = rows - sheet.colheads;
    int err = 0;

    sheet.dset = create_new_dataset(nvars, nobs, labels);

    if (sheet.dset == null) {
        err = GretlError.alloc;
    } else if (!sheet.colheads) {
        /* write fallback variable names */
        int i;

        for (i=1; i<nvars; i++) {
            sprintf(sheet.dset.varname[i], "v%d", i);
        }
    }

    return err;
}

private int node_get_vtype_and_content(xmlNodePtr p, int * vtype,
                                       char ** content)
{
    char * tmp;
    int err = 0;

    tmp = cast(char *) xmlGetProp(p, cast(XUC) "ValueType");

    if (tmp != null) {
        *vtype = atoi(tmp);
        free(tmp);
        if (content != null) {
            *content = cast(char *) xmlNodeGetContent(p);
        }
    } else {
        err = GretlError.data;
    }

    return err;
}

private int inspect_top_left(xmlNodePtr p, WSheet * sheet)
{
    char * content = null;
    int err, vtype = 0;

    err = node_get_vtype_and_content(p, &vtype, &content);

    if (!err) {
        if (vtype == ValueType.empty) {
            sheet.flags |= BookFlag.obs_labels;
        } else if (vtype == ValueType.string_) {
            if (sheet.opt & Gretlopt.OPT_A) {
                if (content != null && *content != '\0') {
                    sheet.colheads = 1;
                }
            } else if (import_obs_label(content)) {
                sheet.flags |= BookFlag.obs_labels;
            }
        }
    }

    free(content);

    return err;
}

private int check_for_vname(xmlNodePtr p, WSheet * sheet)
{
    char * content = null;
    int err, vtype = 0;

    err = node_get_vtype_and_content(p, &vtype, &content);

    if (!err && vtype == ValueType.string_ &&
        content != null && *content != '\0') {
        sheet.colheads = 1;
    }

    free(content);

    return err;
}

private int cell_get_index(xmlNodePtr p, const(char) * which)
{
    char * s = cast(char *) xmlGetProp(p, cast(XUC) which);
    int ret = -1;

    if (s != null) {
        ret = atoi(s);
        free(s);
    }

    return ret;
}

/* Crawl over all the cells and determine the maximum row and column
   indices. While we're at it, inspect the top left cell.
*/

private int wsheet_get_real_size_etc(xmlNodePtr node, WSheet * sheet)
{
    xmlNodePtr p = node.children;
    int topleft_found = 0;
    int i, j;
    int err = 0;

    sheet.maxrow = 0;
    sheet.maxcol = 0;

    while (p != null && !err) {
        if (!xmlStrcmp(p.name, cast(XUC) "Cell")) {
            i = cell_get_index(p, "Row");
            if (i > sheet.maxrow) {
                sheet.maxrow = i;
            }
            j = cell_get_index(p, "Col");
            if (j > sheet.maxcol) {
                /* watch out for bogus cols! */
                sheet.maxcol = j;
            }
            if (i == sheet.row_offset) {
                if (j == sheet.col_offset) {
                    topleft_found = 1;
                    err = inspect_top_left(p, sheet);
                } else if (j > sheet.col_offset) {
                    /* do we have varnames? */
                    err = check_for_vname(p, sheet);
                }
            }
        }
        p = p.next;
    }

    if (!err) {
        if (!topleft_found) {
            /* No top-left cell at all: a sign that we have
               an observations column?
            */
            sheet.flags |= BookFlag.obs_labels;
        }
        fprintf(stderr, "wsheet_get_real_size: maxrow=%d, maxcol=%d\n",
                sheet.maxrow, sheet.maxcol);
    }

    return err;
}

private int gnumeric_non_numeric_check(WSheet * sheet, Prn * prn)
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

private int cell_get_data2(WSheet * sheet,
                           int i, int t,
                           const(char) * s,
                           Prn * prn)
{
    int err = 0;

    if (i > 0 && t >= 0 && sheet.dset.Z[i][t] == NON_NUMERIC) {
        int ix = gretl_string_table_index(sheet.st, s, i, 0, prn);

        if (ix > 0) {
            sheet.dset.Z[i][t] = cast(double) ix;
        } else {
            err = GretlError.data;
        }
    }

    return err;
}

private double cell_get_data(int vtype, const(char) * s,
                             int * any_nonnum, int * err)
{
    double x = NADBL;

    if (vtype_is_numeric(vtype)) {
        x = atof(s);
    } else if (vtype == ValueType.string_) {
        if (string_is_blank(s)) {
            ; /* OK, NA */
        } else if (import_na_string(s)) {
            ; /* OK, NA */
        } else if (numeric_string(s)) {
            x = atof(s);
        } else {
            *any_nonnum = 1;
            x = NON_NUMERIC;
        }
    } else if (vtype != ValueType.empty) {
        *err = GretlError.data;
    }

    return x;
}

private int gnumeric_set_obs_label(WSheet * sheet, int real_t,
                                   xmlNodePtr p, int vtype,
                                   const(char) * s)
{
    int err = 0;

    if (vtype_is_numeric(vtype)) {
        char * fmt = cast(char *) xmlGetProp(p, cast(XUC) "ValueFormat");
        int done = 0;

        /* can we read the numeric value as a Microsoft-style
           date? */

        if (fmt != null) {
            if (strchr(fmt, '/') ||
                (strstr(fmt, "mm") && !(strchr(fmt, ':'))) ||
                strstr(fmt, "yy")) {
                char[12] targ = 0;
                int mst = atoi(s);

                MS_excel_date_string(targ.ptr, mst, 0, 0);
                strcpy(sheet.dset.S[real_t], targ.ptr);
                done = 1;
            }
            free(fmt);
        }
        if (done) {
            return 0;
        }
    }

    if (vtype_is_numeric(vtype) || vtype == ValueType.string_) {
        gretl_utf8_strncat_trim(sheet.dset.S[real_t], s, OBSLEN - 1);
    } else {
        err = GretlError.data;
    }

    return err;
}

/* Note that below we're being agnostic regarding the presence/absence
   of observation labels in the first column. We're writing first
   column values into sheet->labels if they're of string type and
   entering them into sheet->Z if they're numeric. Once we're finished
   we can decide what to do with the first column and the labels.
*/

private int wsheet_parse_cells(xmlNodePtr node, WSheet * sheet,
                               Prn * prn)
{
    xmlNodePtr p;
    xmlNodePtr top = node.children;
    char * tmp;
    int vtype = 0;
    int any_nonnum = 0;
    int have_labels;
    int cols, rows;
    int i, t, r, c;
    int pass = 1;
    int err = 0;

    cols = sheet.maxcol + 1 - sheet.col_offset;
    rows = sheet.maxrow + 1 - sheet.row_offset;

    if (rows < 1) {
        pputs(prn, _("Starting row is out of bounds.\n").ptr);
        return 1;
    }

    if (cols < 1) {
        pputs(prn, _("Starting column is out of bounds.\n").ptr);
        return 1;
    }

    if (wsheet_allocate(sheet, cols, rows)) {
        return 1;
    }

    if (sheet.opt & Gretlopt.OPT_A) {
        /* --all-cols */
        have_labels = 0;
    } else {
        have_labels = (sheet.flags & BookFlag.obs_labels)? 1 : 0;
    }

 tryagain:

    p = top;

    while (p != null && !err) {
        if (!xmlStrcmp(p.name, cast(XUC) "Cell")) {
            int dset_i, dset_t;

            /* what column are we in? */
            c = cell_get_index(p, "Col");
            i = c - sheet.col_offset;

            /* what row are we on? */
            r = cell_get_index(p, "Row");
            t = r - sheet.row_offset;

            if (i < 0 || t < 0) {
                /* we're not in the user-specified reading area */
                p = p.next;
                continue;
            }

            /* get cell type and content */
            err = node_get_vtype_and_content(p, &vtype, &tmp);
            if (err) {
                /* a formula perhaps? */
                pprintf(prn, (_("Couldn't get value for col %d, row %d.\n")
                              ~ _("Maybe there's a formula in the sheet?")).ptr,
                        c+1, r+1);
                break;
            }

            dset_i = i + 1 - have_labels;
            dset_t = t - sheet.colheads;

            if (i == 0 && t == 0 && have_labels) {
                ; /* obs labels heading, ignore */
            } else if (pass == 1 && i == 0 && have_labels) {
                /* should be obs label */
                err = gnumeric_set_obs_label(sheet, dset_t, p, vtype, tmp);
            } else if (pass == 1 && t == 0 && sheet.colheads) {
                /* should be varname */
                if (vtype == ValueType.string_) {
                    strncat(sheet.dset.varname[dset_i], tmp, VNAMELEN - 1);
                    err = check_imported_varname(sheet.dset.varname[dset_i],
                                                 i, r, c, prn);
                } else {
                    err = GretlError.data;
                }
            } else if (pass == 1) {
                /* should be actual data */
                double x = cell_get_data(vtype, tmp, &any_nonnum, &err);

                if (!err) {
                    sheet.dset.Z[dset_i][dset_t] = x;
                }
            } else if (in_gretl_list(sheet.codelist, dset_i)) {
                /* second pass for strings */
                err = cell_get_data2(sheet, dset_i, dset_t, tmp, prn);
            }

            free(tmp);
        }
        p = p.next;
    }

    if (pass == 1 && !err && any_nonnum) {
        pass = 2;
        err = gnumeric_non_numeric_check(sheet, prn);
        if (!err && sheet.st != null) {
            goto tryagain;
        }
    }

    return err;
}

private int wsheet_get_data(const(char) * fname, WSheet * sheet,
                            Prn * prn)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlNodePtr sub;
    char * tmp = null;
    int got_sheet = 0;
    int err;

    err = gretl_xml_open_doc_root(fname, "Workbook", &doc, &cur);
    if (err) {
        return err;
    }

    cur = cur.children;

    /* Now walk the tree */

    while (!err && cur != null && !got_sheet) {
        if (!xmlStrcmp(cur.name, cast(XUC) "Sheets")) {
            int sheetcount = 0;

            sub = cur.children;

            while (sub != null && !got_sheet && !err) {
                if (!xmlStrcmp(sub.name, cast(XUC) "Sheet")) {
                    xmlNodePtr snode = sub.children;

                    while (snode != null && !err) {
                        if (!xmlStrcmp(snode.name, cast(XUC) "Name")) {
                            sheetcount++;
                            tmp = cast(char *) xmlNodeGetContent(snode);
                            if (tmp) {
                                tailstrip(tmp);
                                if (!strcmp(tmp, sheet.name) &&
                                    sheetcount == sheet.ID + 1) {
                                    got_sheet = 1;
                                }
                                free(tmp);
                            }
                        } else if (got_sheet && !xmlStrcmp(snode.name, cast(XUC) "Cells")) {
                            err = wsheet_get_real_size_etc(snode, sheet);
                            if (!err) {
                                err = wsheet_parse_cells(snode, sheet, prn);
                            }
                        }
                        snode = snode.next;
                    }
                }
                sub = sub.next;
            }
        }
        cur = cur.next;
    }

    xmlFreeDoc(doc);

    if (!got_sheet) {
        fprintf(stderr, "'%s': couldn't find the requested sheet\n", sheet.name);
        err = 1;
    }

    return err;
}

private int wbook_record_name(char * name, WBook * book)
{
    char ** sheetnames;
    int ns = book.nsheets + 1;

    sheetnames = cast(char **) realloc(book.sheetnames, ns * (*sheetnames).sizeof);
    if (sheetnames == null) {
        return 1;
    }

    book.sheetnames = sheetnames;
    book.nsheets = ns;
    book.sheetnames[ns - 1] = name;
    tailstrip(name);

    return 0;
}

private int wbook_get_info(const(char) * fname, const(int) * list,
                           char * sheetname, WBook * book,
                           Prn * prn)
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlNodePtr sub;
    char * tmp = null;
    int got_index = 0;
    int err = 0;

    err = gretl_xml_open_doc_root(fname, "Workbook",
                                  &doc, &cur);
    if (err) {
        return err;
    }

    wbook_init(book, list, sheetname);

    /* Now walk the tree */
    cur = cur.children;
    while (cur != null && !got_index && !err) {
        if (!xmlStrcmp(cur.name, cast(XUC) "SheetNameIndex")) {
            got_index = 1;
            sub = cur.children;
            while (sub != null && !err) {
                if (!xmlStrcmp(sub.name, cast(XUC) "SheetName")) {
                    tmp = cast(char *) xmlNodeGetContent(sub);
                    if (tmp != null) {
                        if (wbook_record_name(tmp, book)) {
                            err = 1;
                            free(tmp);
                        }
                    }
                }
                sub = sub.next;
            }
        }
        cur = cur.next;
    }

    xmlFreeDoc(doc);

    return err;
}

private int wsheet_setup(WSheet * sheet, WBook * book, int n)
{
    int err = 0;

    sheet.name = gretl_strdup(book.sheetnames[n]);

    if (sheet.name == null) {
        err = GretlError.alloc;
    } else {
        sheet.ID = n;
        sheet.col_offset = book.col_offset;
        sheet.row_offset = book.row_offset;
    }

    return err;
}

private int finalize_gnumeric_import(DatasetRaw * dset,
                                     WSheet * sheet,
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
        dataset_add_import_info(dset, fname, GretlFileType.gnumeric);
    }

    return err;
}

extern (C) int gnumeric_get_data(const(char) * fname, int * list, char * sheetname,
                                 DatasetRaw * dset, Gretlopt opt, Prn * prn,
                                 void * parent)
{
    int gui = (opt & Gretlopt.OPT_G) != 0;
    WBook gbook;
    WBook * book = &gbook;
    WSheet gsheet;
    WSheet * sheet = &gsheet;
    int sheetnum = -1;
    int err = 0;

    /* C: gbook is uninitialised here, and wbook_free(book) below runs even when wbook_get_info() failed
       before calling wbook_init() (see bugs.md #67); D's default-initialised struct makes that safe */

    wsheet_init(sheet, opt);

    gretl_push_c_numeric_locale();

    if (wbook_get_info(fname, list, sheetname, book, prn)) {
        pputs(prn, _("Failed to get workbook info").ptr);
        err = 1;
        goto getout;
    }

    wbook_print_info(book);

    if (book.nsheets == 0) {
        pputs(prn, _("No worksheets found").ptr);
        err = 1;
        goto getout;
    }

    /* C: if (gui) the sheet dialog was run here; without it both cases use the parameters in @list
       and @sheetname */
    err = wbook_check_params(book);
    if (err) {
        gretl_errmsg_set(_("Invalid argument for worksheet import").ptr);
    } else if (book.selected >= 0) {
        sheetnum = book.selected;
    } else {
        sheetnum = 0;
    }

    if (book.selected == -1) {
        /* canceled */
        err = -1;
    }

    if (!err && sheetnum >= 0) {
        fprintf(stderr, "Getting data...\n");
        err = wsheet_setup(sheet, book, sheetnum);
        if (!err) {
            err = wsheet_get_data(fname, sheet, prn);
            if (err) {
                fprintf(stderr, "wsheet_get_data returned %d\n", err);
            } else {
                book.flags |= sheet.flags;
            }
        }
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
        err = finalize_gnumeric_import(dset, sheet, fname, opt, prn);
        if (!err && gui) {
            wbook_record_params(book, list);
        }
    }

 getout:

    wbook_free(book);
    wsheet_free(sheet);

    gretl_pop_c_numeric_locale();

    return err;
}
