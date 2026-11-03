/* Port of plugin/import_common.c: helpers shared by the spreadsheet importers.

   In the C file this is #included into each plugin, with macros (XLSX_IMPORTER, EXCEL_IMPORTER, ...)
   selecting parts; here it is one module with all the non-GUI parts. The GTK dialog code of that
   file (the "Start import at" / sheet-selection dialog and its callbacks) is dropped: there is no
   GUI in the port, so the sheet, column offset and row offset come only from the command-line style
   arguments. */

module gretl.plugin.import_common;


import core.stdc.ctype : isalnum, isalpha, isdigit;
import core.stdc.errno : errno;
import core.stdc.stdio : fprintf, sprintf, sscanf, stderr;
import core.stdc.stdlib : atoi, free;
import core.stdc.string : strcmp, strcpy, strlen, strncat;
version (Posix) import core.sys.posix.stdlib : mkdtemp; else import gretl.winport : mkdtemp;

import gretl.csvdata : reverse_data, test_markers_for_dates;
import gretl.dataio : check_varname;
import gretl.dataseth : DatasetMarkerType;
import gretl.dataset : dataset_allocate_obs_markers, dataset_destroy_obs_markers, start_new_Z;
import gretl.glib : gchar;
import gretl.gfileutils : g_build_filename, g_get_current_dir, g_path_is_absolute;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_errmsg_get, gretl_errmsg_set, gretl_errmsg_set_from_errno, gretl_errmsg_sprintf,
    gretl_warnmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_paths : gretl_chdir, gretl_deltree, gretl_dotdir, gretl_mkdir, gretl_remove, gretl_test_fopen;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new, pprintf, pputc, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.plugin.importerh;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.strutils : gretl_strdup, numeric_string, tailstrip;
import gretl.gretl_zip : gretl_unzip;
import gretl.calendar : MS_excel_date_string;
import gretl.gretl_intl : iso_to_ascii;

enum E_OBSCOL = 99;

/* unported dependencies: per-file declarations from the 2026c headers */
void strip_vname_illegals(char * s)
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

private int missing_varname()
{
    gretl_errmsg_set(_("Variable name is missing").ptr);
    return GretlError.data;
}

/* The following is modeled on process_csv_varname() in csvdata.c,
   adapted slightly for the specifics of the spreadsheet
   importers.
*/

int check_imported_varname(char * vname, int vnum,
                           int row, int col,
                           Prn * prn)
{
    int err = 0;

    if (*vname == '\0') {
        if (vnum > 0) {
            fprintf(stderr, "variable name %d is missing\n", vnum);
            sprintf(vname, "v%d", vnum);
        } else {
            err = missing_varname();
        }
    } else if (numeric_string(vname)) {
        err = check_varname(vname);
    } else {
        char * s;
        char[VNAMELEN] tmp = 0;

        tmp[0] = '\0';
        strncat(tmp.ptr, vname, VNAMELEN - 1);
        s = tmp.ptr;
        *vname = '\0';

        while (*s && !isalpha(*s)) s++;
        if (*s == '\0') {
            if (vnum > 0) {
                fprintf(stderr, "variable name %d is garbage\n", vnum);
                sprintf(vname, "v%d", vnum);
            } else {
                err = missing_varname();
            }
        } else {
            strncat(vname, s, VNAMELEN - 1);
        }
        iso_to_ascii(vname);
        strip_vname_illegals(vname);
        err = check_varname(vname);
    }

    if (err) {
        err = GretlError.data;
        if (row >= 0 && col >= 0) {
            pputc(prn, '\n');
            pprintf(prn, _("At row %d, column %d:\n").ptr, row+1, col+1);
        }
        pputs(prn, gretl_errmsg_get());
    }

    return err;
}

/* C: "#ifndef EXCEL_IMPORTER" */

int import_ts_check(DatasetRaw * dset)
{
    Prn * prn = gretl_print_new(PrnType.stderr, null);
    int reversed = 0;
    int mpd;

    mpd = test_markers_for_dates(dset, &reversed, null, prn);

    if (mpd > 0) {
        pputs(prn, _("taking date information from row labels\n\n").ptr);
        if (dset.markers != DatasetMarkerType.daily_date_strings) {
            dataset_destroy_obs_markers(dset);
        }
        if (reversed) {
            reverse_data(dset, prn);
        }
    }

    if (dset.pd != 1 || strcmp(dset.stobs.ptr, "1")) {
        dset.structure = TsCodes.time_series;
    }

    gretl_print_destroy(prn);

    return mpd;
}

/* check for empty columns at the right of the sheet */

int import_prune_columns(DatasetRaw * dset)
{
    int allmiss = 1, n_empty = 0;
    int i, t, err = 0;

    for (i=dset.v-1; i>0 && allmiss; i--) {
        for (t=0; t<dset.n; t++) {
            if (!na(dset.Z[i][t])) {
                allmiss = 0;
                break;
            }
        }
        if (allmiss) n_empty++;
    }

    if (n_empty == dset.v - 1) {
        gretl_errmsg_set(_("No numeric data were found").ptr);
        err = GretlError.data;
    } else if (n_empty > 0) {
        /* changed 2025-05-21 */
        gretl_warnmsg_sprintf(_("Sheet has %d trailing empty variables").ptr, n_empty);
    }

    return err;
}

/* we want this for unzipping purposes */

private gchar * get_absolute_path(const(char) * fname)
{
    gchar * cwd = g_get_current_dir();
    gchar * ret = null;

    if (cwd != null) {
        ret = g_build_filename(cwd, fname, null);
    }

    return ret;
}

void remove_temp_dir(char * dname)
{
    if (gretl_chdir(gretl_dotdir()) == 0) {
        gretl_deltree(dname);
    }
}

private int gretl_make_tempdir(char * dname)
{
    char * s;
    int err = 0;

    strcpy(dname, ".gretl-ssheet-XXXXXX");
    s = mkdtemp(dname);

    if (s == null) {
        gretl_errmsg_set_from_errno("gretl_make_tempdir", errno);
        err = GretlError.fopen;
    }

    return err;
}

/* For ODS and XLSX: unzip the target file in the user's "dotdir". On
   successful completion @dname holds the name of the temporary
   subdirectory, in the dotdir, holding the contents of the zipfile.
*/

int open_import_zipfile(const(char) * fname, char * dname,
                        Prn * prn)
{
    const(char) * real_fname = fname;
    gchar * abspath = null;
    int err = 0;

    errno = 0;
    *dname = '\0';

    if (gretl_test_fopen(real_fname, "r") != 0) {
        return GretlError.fopen;
    }

    /* by doing chdir, we may lose track of the file if
       its path is relative */
    if (!g_path_is_absolute(real_fname)) {
        abspath = get_absolute_path(real_fname);
        if (abspath != null) {
            real_fname = abspath;
        }
    }

    /* cd to dotdir and make temporary dir */
    if (gretl_chdir(gretl_dotdir()) != 0) {
        err = GretlError.fopen;
    } else {
        err = gretl_make_tempdir(dname);
        if (!err) {
            err = gretl_chdir(dname);
            if (err) {
                gretl_remove(dname);
            }
        }
    }

    if (!err) {
        /* if all has gone OK, we're now in the temporary
           directory under dotdir, and @real_fname is the
           absolute path to the file to be unzipped.
        */
        err = gretl_unzip(real_fname);
        if (err) {
            pprintf(prn, "gretl_unzip: %s\n", gretl_errmsg_get());
        }
    }

    g_free(abspath);

    return err;
}

/* the parts for the non-ODS importers (Excel, Gnumeric, ...) */

int worksheet_start_dataset(DatasetRaw * newinfo)
{
    if (newinfo.v == 1) {
        /* only the constant is present! */
        gretl_errmsg_set(_("No numeric data were found").ptr);
        return GretlError.data;
    } else {
        /* create import dataset */
        return start_new_Z(newinfo, cast(Gretlopt) 0);
    }
}

int importer_dates_check(char ** labels, BookFlag * pflags,
                         DatasetRaw * newset, Prn * prn,
                         int * err)
{
    int d, t;
    char[12] dstr = 0;
    char * s;
    int ret = 0;

    for (t=0; t<newset.n; t++) {
        s = labels[t];
        if (s == null || *s == '\0') {
            fprintf(stderr, "importer_dates_check: got blank label\n");
            return 0;
        }
    }

    *err = dataset_allocate_obs_markers(newset);
    if (*err) {
        return 0;
    }

    for (t=0; t<newset.n && !*err; t++) {
        s = labels[t];
        if (*s == '"' || *s == '\'') s++;
        if (*pflags & BookFlag.numeric_dates) {
            if (sscanf(s, "%d", &d)) {
                MS_excel_date_string(dstr.ptr, d, 0, (*pflags & BookFlag.date_base_1904) != 0);
                s = dstr.ptr;
            } else {
                pprintf(prn, _("Bad date on row %d: '%s'\n").ptr, t+1, s);
                *err = GretlError.data;
            }
        }
        strncat(newset.S[t], s, OBSLEN - 1);
    }

    if (!*err) {
        int reversed = 0;

        ret = test_markers_for_dates(newset, &reversed, null, prn);
        if (reversed) {
            *pflags |= BookFlag.data_reversed;
        } else if (ret < 0) {
            /* not really time series */
            fprintf(stderr, "importer_dates_check: scrubbing time series\n");
        }
    }

    if (newset.markers != DatasetMarkerType.daily_date_strings) {
        dataset_destroy_obs_markers(newset);
    }

    return ret;
}

void wbook_print_info(WBook * book)
{
    int i;

    fprintf(stderr, "Found %d sheet%s\n", book.nsheets,
            (book.nsheets > 1)? "s".ptr : "".ptr);

    for (i=0; i<book.nsheets; i++) {
        if (book.byte_offsets != null) {
            fprintf(stderr, "%d: '%s' at offset %u\n", i,
                    book.sheetnames[i], book.byte_offsets[i]);
        } else {
            fprintf(stderr, "%d: '%s'\n", i, book.sheetnames[i]);
        }
    }
}

void wbook_free(WBook * book)
{
    int i;

    for (i=0; i<book.nsheets; i++) {
        free(book.sheetnames[i]);
    }
    free(book.sheetnames);
    free(book.targname);
    free(book.byte_offsets);
    free(book.xf_list);
}

int wbook_check_params(WBook * book)
{
    if (book.targname != null) {
        int i;

        book.selected = -1;
        for (i=0; i<book.nsheets; i++) {
            if (!strcmp(book.targname, book.sheetnames[i])) {
                book.selected = i;
            }
        }
        if (book.selected < 0) {
            gretl_errmsg_sprintf(_("\"%s\": no such sheet").ptr, book.targname);
            return GretlError.data;
        }
    }

    if (book.selected < 0 || book.selected >= book.nsheets) {
        return GretlError.data;
    } else if (book.col_offset < 0 || book.row_offset < 0) {
        return GretlError.data;
    } else {
        return 0;
    }
}

void wbook_record_params(WBook * book, int * list)
{
    if (list != null && list[0] == 3) {
        list[1] = book.selected + 1;
        list[2] = book.col_offset;
        list[3] = book.row_offset;
    }
}

/* @list may contain sheet number, row and/or column offset;
   @sheetname may contain the name of a specific sheet; but
   both may be NULL
*/

void wbook_init(WBook * book, const(int) * list, char * sheetname)
{
    book.version_ = 0;
    book.nsheets = 0;
    book.col_offset = book.row_offset = 0;
    book.targname = null;
    book.sheetnames = null;
    book.byte_offsets = null;
    book.selected = 0;
    book.flags = cast(BookFlag) 0;
    book.xf_list = null;
    book.get_min_offset = null;
    book.data = null;

    if (sheetname != null && *sheetname != '\0') {
        book.targname = gretl_strdup(sheetname);
        tailstrip(book.targname);
    }

    if (list != null && list[0] == 3) {
        if (book.targname == null && list[1] > 0) {
            book.selected = list[1] - 1;
        }
        book.col_offset = list[2];
        book.row_offset = list[3];
    }
}

private __gshared char[5] label_buf = 0;

const(char) * column_label(int col)
{
    static immutable string abc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    char a1, a2;

    if (col < 26) {
        a2 = abc[col];
        sprintf(label_buf.ptr, "(%c)", a2);
    } else {
        a1 = abc[col / 26 - 1];
        a2 = abc[col % 26];
        sprintf(label_buf.ptr, "(%c%c)", a1, a2);
    }

    return label_buf.ptr;
}

int book_get_min_offset(WBook * book, int k)
{
    if (book.get_min_offset != null) {
        return book.get_min_offset(book, k);
    } else {
        return 1;
    }
}
