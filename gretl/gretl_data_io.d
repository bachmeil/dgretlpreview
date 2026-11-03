/* Port of gretl_data_io.c: wrappers for functions in dataio.c and
   gretl_xml.c, presenting a consistent, simplified API for reading
   and writing data files. */

module gretl.gretl_data_io;


import core.stdc.stdio : fprintf, stderr;

import gretl.csvdata : import_csv;
import gretl.dataio : detect_filetype, import_other, import_spreadsheet;
import gretl.dataioh : GretlFileType, OTHER_IMPORT, SPREADSHEET_IMPORT;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_prn : Prn;
import gretl.gretl_xml : gretl_read_gdt, gretl_write_gdt;
import gretl.libgretlh;

/**
 * gretl_read_native_data:
 * @fname: path to a native gretl (.gdt) data file.
 * @dset: dataset struct.
 *
 * Read data from file into gretl's work space, allocating memory
 * as required.
 *
 * The argument @dset represents a pointer-to-DATASET. It should
 * either be given as the address of a #DATASET struct that exists
 * at the caller level, or it can be a pointer obtained via the
 * libgretl function datainfo_new().
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_read_native_data(const(char) * fname, DatasetRaw * dset)
{
    int err = 0;

    if (fname == null || dset == null) {
        err = GretlError.invarg;
    } else if (dset.Z != null) {
        fprintf(stderr, "gretl_read_native_data: Z must be NULL on entry\n");
        err = GretlError.invarg;
    } else {
        err = gretl_read_gdt(fname, dset, Gretlopt.OPT_NONE, null);
    }

    return err;
}

/**
 * gretl_write_native_data:
 * @fname: name of file to write.
 * @list: list of ID numbers of series to write (or NULL to write all).
 * @dset: dataset struct.
 *
 * Write out in native gretl (.gdt) format a data file containing
 * the values of the given set of variables.
 *
 * Returns: 0 on successful completion, non-zero on error.
 */

int gretl_write_native_data(const(char) * fname, const(int) * list,
                            const(DatasetRaw) * dset)
{
    if (fname == null || dset == null || dset.Z == null) {
        return GretlError.invarg;
    }

    return gretl_write_gdt(fname, list, dset, Gretlopt.OPT_Z, 0);
}

/**
 * gretl_read_foreign_data:
 * @fname: path to a data file of a type that gretl can handle.
 * @file_type: code representing the format of the data file.
 * @dset: dataset struct.
 * @prn: printer for diagnostic info, or NULL.
 *
 * Read data from a "foreign" format data file into gretl's work space,
 * allocating memory as required. For comments on the arguments @pZ and @dset,
 * see gretl_read_native_data().
 *
 * @file_type must be one of GRETL_CSV, GRETL_OCTAVE,
 * GRETL_GNUMERIC, GRETL_XLS, GRETL_XLSX, GRETL_ODS, GRETL_WF1,
 * GRETL_DTA, GRETL_SAV, GRETL_SAS or GRETL_JMULTI. If you are
 * unsure of the type of the file you may call gretl_detect_filetype()
 * first. Note that the GRETL_CSV type is quite "permissive",
 * including plain text data files with the data values separated
 * in various ways, not just comma-separated values in the strict
 * sense. (The separation must, however, be consistent within the
 * given file.)
 *
 * If @prn is non-NULL, diagnostic information will be printed.
 * This can be useful to gauge how successful the import was.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gretl_read_foreign_data(const(char) * fname, GretlFileType file_type,
                            DatasetRaw * dset, Prn * prn)
{
    int err = 0;

    if (fname == null || dset == null) {
        err = GretlError.invarg;
    } else if (dset.Z != null) {
        fprintf(stderr, "gretl_read_foreign_data: Z must be NULL on entry\n");
        err = GretlError.invarg;
    } /* C: no "else" before the next "if" (see bugs.md #55) */
    if (file_type == GretlFileType.csv) {
        err = import_csv(fname, dset, Gretlopt.OPT_NONE, prn);
    } else if (SPREADSHEET_IMPORT(file_type)) {
        err = import_spreadsheet(fname, file_type, null, null,
                                 dset, Gretlopt.OPT_NONE, prn);
    } else if (OTHER_IMPORT(file_type)) {
        err = import_other(fname, file_type, dset, Gretlopt.OPT_NONE, prn);
    } else {
        gretl_errmsg_set(_("Unknown data import type").ptr);
        err = GretlError.invarg;
    }

    return err;
}

/**
 * gretl_detect_filetype:
 * @fname: name of the file to be examined.
 *
 * Attempts to determine if the named file is of a type from
 * which gretl can read data.
 *
 * Returns: code representing the type of the data file, or
 * %GRETL_UNRECOGNIZED if the file doesn't seem to be something
 * gretl can work with.
 */

GretlFileType gretl_detect_filetype(const(char) * fname)
{
    if (fname == null) {
        return GretlFileType.unrecognized;
    }

    /* Note: OPT_NONE means that @fname will not be
       modified, so the cast to (char *) is alright.
    */
    return detect_filetype(cast(char *) fname, Gretlopt.OPT_NONE);
}
