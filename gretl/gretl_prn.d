/* Port of gretl_prn.c */

module gretl.gretl_prn;


/*
 * gretl printing struct
 *
 * Most libgretl functions that print output call for a
 * pointer-to-Prn as one of their arguments. The Prn type is
 * an opaque structure, constructed and manipulated by the
 * functions listed here. It is used as a means of generalizing
 * the printing operation, which may be to a regular file, a
 * temporary file or a buffer.
 *
 * To get hold of a Prn use gretl_print_new() or one of the more
 * specific constructors, and to free it use gretl_print_destroy().
 * If you want to use a Prn directly for printing in your own code,
 * use the functions pprintf(), pputs() and pputc(). These are
 * counterparts to the standard C functions fprintf, fputs and
 * fputc, but note that with pputs and pputc the Prn argument
 * must be given first (unlike fputs and fputc in which the
 * FILE argument goes last).
 *
 * Note that whenever a Prn appears as a function parameter
 * in libgretl it is OK to give a null argument: in that case
 * pprintf(), pputs() and pputc() are no-ops.
 */

import core.stdc.stdarg : va_end, va_list, va_start;
import core.stdc.stdio : fclose, fflush, fputc, fputs, fread, fseek, ftell,
    FILE, fprintf, rewind, SEEK_END, SEEK_SET, sprintf, stderr, stdout,
    vfprintf, vsnprintf;
import core.stdc.stdlib : calloc, free, malloc, realloc;
import core.stdc.string : memset, strcpy, strlen, strncmp, strstr;
import etc.c.zlib : gzclose, gzFile, gzflush, gzputs, gzsetparams,
    Z_DEFAULT_STRATEGY, Z_SYNC_FLUSH;

import gretl.garray;
import gretl.glib;
import gretl.gmem;
import gretl.gretl_commandsh;
import gretl.gretl_errorsh;
import gretl.gretl_prnh;
import gretl.gstrfuncs;
import gretl.libgretlh;
import gretl.strutils : gretl_strdup;
import gretl.gretl_errors : gretl_errmsg_get;
import gretl.uservarh : UserVar;
import gretl.gretl_paths : gretl_dotdir, gretl_fopen, gretl_gzopen, gretl_mktemp, gretl_remove, gretl_rename;
import gretl.printout : bufgets, bufgets_finalize, bufgets_init;
import gretl.uservar : get_user_var_of_type_by_name, user_var_replace_value;

/* ---- declarations of not-yet-ported items used below ----
   Replace each with an import when its defining file is ported. */

extern (C) {
    /* gretl utilities */
    int gretl_function_depth();

    /* user variables */
}

/* ---- end of declarations ---- */

/* An opaque structure accessed only via gretl_print functions. */
struct Prn {
    FILE * fp;          /* file to which to print, or null */
    gzFile fz;          /* gzipped file target, or null */
    char * buf;         /* buffer to which to print, or null */
    size_t bufsize;     /* allocated size of buffer */
    size_t blen;        /* string length of buffer */
    int savepos;        /* saved position in stream or buffer */
    GArray * fplist;    /* stack for use with output redirection */
    PrnFormat format;   /* plain, TeX, RTF */
    gint8 fixed;        /* non-zero for fixed-size buffer */
    gint8 gbuf;         /* non-zero for buffer obtained via GLib */
    guint8 nlcount;     /* count of trailing newlines */
    char delim;         /* CSV field delimiter */
    char * fname;       /* temp file name, or null */
}

private struct FpInfo {
    FILE * fp;       /* stream to which we're printing */
    int level;       /* level of redirection */
    gchar * fname;   /* name of file or null */
    gchar * strvar;  /* associated string variable or null */
}

/* equivalent of g_array_index(a, FpInfo, i) */
private FpInfo * fplist_index(GArray * a, int i)
{
    return (cast(FpInfo *) a.data) + i;
}

private void prn_destroy_fp_list(Prn * prn)
{
    int i, n = prn.fplist.len;
    FpInfo * fi;

    for (i=n-1; i>=0; i--) {
        fi = fplist_index(prn.fplist, i);
        if (fi != null) {
            if (fi.fp != null && fi.fp != prn.fp &&
                fi.fp != stdout && fi.fp != stderr) {
                fclose(fi.fp);
            }
            if (fi.fname != null) {
                g_free(fi.fname);
            }
            if (fi.strvar != null) {
                g_free(fi.strvar);
            }
            g_array_remove_index(prn.fplist, i);
        }
    }

    g_array_free(prn.fplist, 1);
    prn.fplist = null;
}

private int prn_fp_list_active(Prn * prn)
{
    return prn.fplist != null && prn.fplist.len > 0;
}

/*
 * gretl_print_destroy:
 * @prn: pointer to gretl printing struct.
 *
 * Close a gretl printing struct and free any associated resources.
 */

void gretl_print_destroy(Prn * prn)
{
    if (prn == null) {
        return;
    }

    if (prn.fplist != null) {
        prn_destroy_fp_list(prn);
    }

    if (prn.fp != null) {
        if (prn.fp == stdout) {
            fflush(stdout);
        } else if (prn.fp != stderr) {
            fclose(prn.fp);
        }
    } else if (prn.fz != null) {
        gzclose(prn.fz);
    }

    if (prn.fname != null) {
        /* prn had a tempfile attached */
        gretl_remove(prn.fname);
        free(prn.fname);
    }

    if (prn.buf != null) {
        if (prn.gbuf) {
            /* the buffer was obtained from GLib */
            g_free(prn.buf);
        } else {
            free(prn.buf);
        }
    }

    free(prn);
}

private int prn_add_tempfile(Prn * prn)
{
    const(char) * dotdir = gretl_dotdir();
    int n = cast(int) strlen(dotdir) + 16;
    int err = 0;

    prn.fname = cast(char *) malloc(n);
    if (prn.fname == null) {
        return GretlError.alloc;
    }

    sprintf(prn.fname, "%sprntmp.XXXXXX", dotdir);
    prn.fp = gretl_mktemp(prn.fname, "w+b");

    if (prn.fp == null) {
        free(prn.fname);
        prn.fname = null;
        err = GretlError.fopen;
    }

    return err;
}

private Prn * real_gretl_print_new(PrnType ptype,
                                   const(char) * fname,
                                   char * buf,
                                   FILE * fp,
                                   int zcomp,
                                   int * perr)
{
    Prn * prn = cast(Prn *) malloc(Prn.sizeof);
    int err = 0;

    if (prn == null) {
        if (perr != null) {
            *perr = GretlError.alloc;
        }
        return null;
    }

    prn.fp = null;
    prn.fz = null;
    prn.fplist = null;
    prn.buf = null;
    prn.bufsize = 0;
    prn.blen = 0;
    prn.savepos = -1;
    prn.format = PrnFormat.txt;
    prn.fixed = 0;
    prn.gbuf = 0;
    prn.nlcount = 0;
    prn.delim = ',';
    prn.fname = null;

    if (ptype == PrnType.stream) {
        prn.fp = fp;
    } else if (ptype == PrnType.file) {
        prn.fp = gretl_fopen(fname, "wb");
        if (prn.fp == null) {
            err = GretlError.fopen;
            free(prn);
            prn = null;
        }
    } else if (ptype == PrnType.gzfile) {
        prn.fz = gretl_gzopen(fname, "wb");
        if (prn.fz == null) {
            err = GretlError.fopen;
            free(prn);
            prn = null;
        } else if (zcomp >= 0) {
            gzsetparams(prn.fz, zcomp, Z_DEFAULT_STRATEGY);
        }
    } else if (ptype == PrnType.tempfile) {
        err = prn_add_tempfile(prn);
        if (err) {
            free(prn);
            prn = null;
        } else {
            prn.savepos = 0;
        }
    } else if (ptype == PrnType.stdout) {
        prn.fp = stdout;
    } else if (ptype == PrnType.stderr) {
        prn.fp = stderr;
    } else if (ptype == PrnType.buffer) {
        if (buf != null) {
            prn.buf = buf;
            prn.fixed = 1;
        } else {
            int p = pprintf(prn, "@init");

            if (p < 0) {
                err = GretlError.alloc;
                free(prn);
                prn = null;
            }
        }
    }

    if (perr != null) {
        *perr = err;
    }

    return prn;
}

/*
 * gretl_print_new:
 * @ptype: code indicating the desired printing mode.
 * @err: location to receive error code, or null.
 *
 * Create and initialize a gretl printing struct. If @ptype
 * is PrnType.buffer, output will go to an automatically
 * resized buffer; if @ptype is PrnType.stdout or
 * PrnType.stderr output goes to stdout or stderr
 * respectively.
 *
 * If you want a named file associated with the struct, use
 * gretl_print_new_with_filename instead; if you want to
 * attach a fixed, pre-allocated text buffer, use
 * gretl_print_new_with_buffer.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new(PrnType ptype, int * err)
{
    if (ptype == PrnType.file) {
        fprintf(stderr, "gretl_print_new: needs a filename\n");
        return null;
    }

    return real_gretl_print_new(ptype, null, null, null, -1, err);
}

/*
 * gretl_print_new_with_filename:
 * @fname: name of the file to be opened for writing.
 * @err: location to receive error code.
 *
 * Create and initialize a gretl printing struct, with
 * output directed to the named file.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new_with_filename(const(char) * fname, int * err)
{
    if (fname == null) {
        fprintf(stderr, "gretl_prn_new: Must supply a filename\n");
        return null;
    }

    return real_gretl_print_new(PrnType.file, fname,
                                null, null, -1, err);
}

/*
 * gretl_gzip_print_new:
 * @fname: name of the compressed file to be opened for writing.
 * @comp_level: -1 for default, or integer 0 to 9.
 * @err: location to receive error code.
 *
 * Create and initialize a gretl printing struct, with
 * output directed to the named compressed file.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_gzip_print_new(const(char) * fname, int comp_level,
                           int * err)
{
    if (fname == null) {
        fprintf(stderr, "gretl_gzip_print_new: must supply a filename\n");
        return null;
    } else if (comp_level < -1 || comp_level > 9) {
        fprintf(stderr, "gretl_gzip_print_new: invalid compression level\n");
        *err = GretlError.invarg;
        return null;
    }

    return real_gretl_print_new(PrnType.gzfile, fname,
                                null, null, comp_level, err);
}

/*
 * gretl_print_new_with_tempfile:
 * @err: location to receive error code.
 *
 * Create and initialize a gretl printing struct, with
 * output directed to a temporary file, which is deleted
 * when the printing struct is destroyed.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new_with_tempfile(int * err)
{
    return real_gretl_print_new(PrnType.tempfile, null,
                                null, null, -1, err);
}

/*
 * gretl_print_has_tempfile:
 * @prn: printing struct to test.
 *
 * Returns: 1 if @prn has a tempfile attached, else 0.
 */

int gretl_print_has_tempfile(Prn * prn)
{
    if (prn != null && prn.fname != null && prn.fp != null) {
        return strstr(prn.fname, "prntmp.") != null;
    } else {
        return 0;
    }
}

/*
 * gretl_print_get_tempfile_name:
 * @prn: printing struct to test.
 *
 * Returns: if @prn has a tempfile attached, return the name
 * of that file, otherwise return null.
 */

const(char) * gretl_print_get_tempfile_name(Prn * prn)
{
    if (prn != null) {
        return prn.fname;
    } else {
        return null;
    }
}

/*
 * gretl_print_new_with_buffer:
 * @buf: pre-allocated text buffer.
 *
 * Creates and initializes a gretl printing struct, with
 * fixed text buffer @buf.  Fails if @buf is null.  This is a
 * convenience function: you can't use pprintf, pputs or
 * pputc with a printing struct obtained in this way.
 *
 * Note that @buf will be freed if and when gretl_print_destroy
 * is called on the Prn pointer obtained.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new_with_buffer(char * buf)
{
    if (buf == null) {
        return null;
    } else {
        return real_gretl_print_new(PrnType.buffer, null,
                                    buf, null, -1, null);
    }
}

/*
 * gretl_print_new_with_gchar_buffer:
 * @buf: pre-allocated text buffer.
 *
 * Just as gretl_print_new_with_buffer() except that the buffer is
 * of pointer-to-gchar type, as obtained from one or other GLib
 * function. This means that when the Prn is detroyed the
 * buffer will be freed using GLib's g_free function rather than
 * the standard C library's free function.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new_with_gchar_buffer(gchar * buf)
{
    Prn * prn = null;

    if (buf != null) {
        prn = real_gretl_print_new(PrnType.buffer, null,
                                   buf, null, -1, null);
        if (prn != null) {
            prn.gbuf = 1;
        }
    }

    return prn;
}

/*
 * gretl_print_new_with_stream:
 * @fp: pre-opened stream.
 *
 * Creates and initializes a gretl printing struct, with
 * printing to @fp.
 *
 * Note that @fp will be closed if and when gretl_print_destroy
 * is called on the Prn pointer obtained.
 *
 * Returns: pointer to newly created struct, or null on failure.
 */

Prn * gretl_print_new_with_stream(FILE * fp)
{
    if (fp == null) {
        return null;
    } else {
        return real_gretl_print_new(PrnType.stream, null,
                                    null, fp, -1, null);
    }
}

/*
 * gretl_print_detach_stream
 * @prn: printing struct to operate on.
 *
 * Sets the stream member of @prn to null so that @prn can
 * be destroyed without closing the associated stream.  May be
 * used in conjunction with gretl_print_new_with_stream().
 */

void gretl_print_detach_stream(Prn * prn)
{
    prn.fp = null;
}

/*
 * gretl_print_rename_file:
 * @prn: printing struct to operate on.
 * @oldpath: name of current file (or null if @prn was
 * set up using gretl_print_new_with_tempfile).
 * @newpath: new name for file.
 *
 * If @prn is printing to a FILE pointer, rename the
 * file to which it is printing and switch the stream
 * to the new file.
 *
 * Returns: 0 on success, 1 on error.
 */

int gretl_print_rename_file(Prn * prn, const(char) * oldpath,
                            const(char) * newpath)
{
    int err = 0;

    if (prn == null) {
        fprintf(stderr, "gretl_print_rename_file: prn is NULL\n");
        return GretlError.data;
    }

    if (prn.fp == null || prn_fp_list_active(prn)) {
        return GretlError.data;
    }

    fclose(prn.fp);
    prn.fp = null;

    if (oldpath == null && prn.fname != null) {
        /* renaming from tempfile */
        err = gretl_rename(prn.fname, newpath);
    } else {
        err = gretl_rename(oldpath, newpath);
    }

    if (err) {
        fprintf(stderr, "%s\n", gretl_errmsg_get());
    } else {
        /* re-open the stream under its new name */
        prn.fp = gretl_fopen(newpath, "a");
        if (prn.fname != null) {
            /* @prn originally used a tempfile: the record of
               the temporary filename should be deleted
            */
            free(prn.fname);
            prn.fname = null;
        }
    }

    return err;
}

/*
 * gretl_print_reset_buffer:
 * @prn: printing struct to operate on.
 *
 * If @prn has an attached buffer, write a NUL byte to
 * the start of the buffer and reset the count of bytes
 * printed to zero.  The next call to @pprintf or
 * similar will then overwite rather than cumulating
 * the printed content.
 *
 * Returns: 0 on success, 1 if @prn has no buffer.
 */

int gretl_print_reset_buffer(Prn * prn)
{
    int err = 0;

    if (prn != null && prn.buf != null) {
        *prn.buf = '\0';
        prn.blen = 0;
    } else {
        err = 1;
    }

    return err;
}

/*
 * gretl_print_get_buffer:
 * @prn: printing struct.
 *
 * Obtain a pointer to the buffer associated
 * with @prn, if any.  This pointer must not be
 * modified in any way.
 *
 * Returns: the buffer, or null on failure.
 */

const(char) * gretl_print_get_buffer(Prn * prn)
{
    return (prn != null) ? prn.buf : null;
}

/*
 * gretl_print_get_trimmed_buffer:
 * @prn: printing struct.
 *
 * Obtain a pointer to the buffer associated
 * with @prn, if any.  This pointer must not be
 * modified in any way.  An opening newline
 * is skipped and any trailing white space is
 * substituted by NULs.
 *
 * Returns: the buffer, or null on failure.
 */

const(char) * gretl_print_get_trimmed_buffer(Prn * prn)
{
    char * buf = (prn != null) ? prn.buf : null;

    if (buf != null) {
        int i, n;

        if (*buf == '\n') {
            buf++;
        }
        n = cast(int) strlen(buf);
        for (i=n-1; i>0; i--) {
            if (buf[i] == '\n' && buf[i-1] == '\n') {
                buf[i] = '\0';
            } else {
                break;
            }
        }
    }

    return buf;
}

/*
 * gretl_print_get_size:
 * @prn: printing struct.
 * @width: location to receive width, or null.
 * @height: location to receive height, or null.
 *
 * If @prn has a non-null buffer attached, provides
 * the width and/or height of the buffer, the width in
 * characters and the height in lines. This function is
 * intended for use with text designed for printing in
 * a GUI window, and with "reasonable" line lengths for
 * that context; if the lines are too long (more than
 * 120 UTF-8 characters) the values written to @width
 * and/or @height will be zero.
 */

void gretl_print_get_size(Prn * prn, int * width, int * height)
{
    int w = 0, h = 0;

    if (prn != null && prn.buf != null) {
        char[128] line;
        int lw;

        bufgets_init(prn.buf);
        while (bufgets(line.ptr, line.sizeof, prn.buf)) {
            lw = cast(int) g_utf8_strlen(line.ptr, -1) - 1;
            if (lw > 120) {
                w = h = 0;
                break;
            }
            if (lw > w) {
                w = lw;
            }
            h++;
        }
        bufgets_finalize(prn.buf);
    }

    if (width != null) {
        *width = w;
    }

    if (height != null) {
        *height = h;
    }
}

/*
 * gretl_print_steal_buffer:
 * @prn: printing struct.
 *
 * Obtain a pointer to the buffer associated with @prn,
 * if any.  The pointer on @prn itself is set to null
 * and the caller takes responsibility for freeing the
 * buffer.
 *
 * Returns: the buffer, or null on failure.
 */

char * gretl_print_steal_buffer(Prn * prn)
{
    char * buf = null;

    if (prn != null) {
        buf = prn.buf;
        prn.buf = null;
    }

    return buf;
}

/*
 * gretl_print_replace_buffer:
 * @prn: printing struct.
 * @buf: malloced replacement buffer.
 *
 * If @prn currently has a printing buffer in place,
 * destroy the original and replace it with @buf. Note
 * @prn "takes ownership" of @buf, which will be freed
 * when @prn is destroyed.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_print_replace_buffer(Prn * prn, char * buf)
{
    int err = 0;

    if (prn == null || prn.buf == null || buf == null) {
        err = GretlError.data;
    } else {
        free(prn.buf);
        prn.buf = buf;
        prn.blen = strlen(buf);
        prn.bufsize = prn.blen + 1;
        prn.savepos = -1;
    }

    return err;
}

/*
 * gretl_print_set_save_position:
 * @prn: printing struct.
 *
 * Sets the current buffer offset as the position from
 * which a chunk of the current buffer may be saved,
 * using gretl_print_get_chunk().
 *
 * Returns: 0 on success, error code if @prn is not
 * connected to a buffer.
 */

int gretl_print_set_save_position(Prn * prn)
{
    if (prn == null || prn.buf == null) {
        return GretlError.data;
    } else {
        prn.savepos = cast(int) prn.blen;
        return 0;
    }
}

/*
 * gretl_print_unset_save_position:
 * @prn: printing struct.
 *
 * Erases the "save position" offset as set by
 * gretl_print_set_save_position().
 */

void gretl_print_unset_save_position(Prn * prn)
{
    if (prn != null) {
        prn.savepos = -1;
    }
}

/*
 * gretl_print_get_chunk:
 * @prn: printing struct.
 *
 * Retrieves a copy of the buffer associated with @prn,
 * starting at the offset from the start of the buffer
 * as set by gretl_print_set_save_position().
 *
 * Returns: allocated buffer on success, or null on error.
 */

char * gretl_print_get_chunk(Prn * prn)
{
    /* NB: as in the C source, @prn is dereferenced here
       before the null check below */
    int maxpos = cast(int) prn.blen;
    char * ret;

    if (prn == null || prn.buf == null ||
        prn.savepos < 0 || prn.savepos > maxpos) {
        return null;
    }

    ret = gretl_strdup(prn.buf + prn.savepos);
    prn.savepos = -1;

    return ret;
}

/*
 * gretl_print_get_chunk_at:
 * @prn: printing struct.
 * @pos: the byte position at which to start.
 *
 * Retrieves a copy of the buffer associated with @prn,
 * starting at the offset from the start of the buffer
 * specified by @pos.
 *
 * Returns: allocated buffer on success, or null on error.
 */

char * gretl_print_get_chunk_at(Prn * prn, int pos)
{
    /* NB: as in the C source, @prn is dereferenced here
       before the null check below */
    int maxpos = cast(int) prn.blen;
    char * ret;

    if (prn == null || prn.buf == null ||
        pos < 0 || pos > maxpos) {
        return null;
    }

    ret = gretl_strdup(prn.buf + pos);

    return ret;
}

/*
 * gretl_print_tell:
 * @prn: printing struct.
 *
 * Returns: the current write position, if @prn
 * has a buffer attached, otherwise -1.
 */

int gretl_print_tell(Prn * prn)
{
    if (prn == null || prn.buf == null) {
        return -1;
    } else {
        return cast(int) prn.blen;
    }
}

private void set_default_csv_delim(Prn * prn)
{
    char[4] test;

    sprintf(test.ptr, "%.1f", 1.0);

    if (test[1] == ',') {
        prn.delim = ';';
    } else {
        prn.delim = ',';
    }
}

/*
 * gretl_print_set_format:
 * @prn: printing struct.
 * @format: desired format code.
 *
 * Sets the printing format on @prn.
 */

void gretl_print_set_format(Prn * prn, PrnFormat format)
{
    if (prn != null) {
        if (format == PrnFormat.rtf) {
            format = PrnFormat.rtf | PrnFormat.doc;
        }
        if (format == PrnFormat.csv) {
            set_default_csv_delim(prn);
        }
        prn.format = format;
    }
}

/*
 * gretl_print_set_delim:
 * @prn: printing struct.
 * @delim: desired CSV field-delimiter.
 *
 * Sets the CSV delimiter on @prn.
 */

void gretl_print_set_delim(Prn * prn, char delim)
{
    if (prn != null) {
        prn.delim = delim;
    }
}

/*
 * gretl_print_toggle_doc_flag:
 * @prn: printing struct.
 *
 * Toggles the PrnFormat.doc flag on @prn.
 */

void gretl_print_toggle_doc_flag(Prn * prn)
{
    if (prn != null) {
        prn.format ^= PrnFormat.doc;
    }
}

/*
 * gretl_print_set_has_minus:
 * @prn: printing struct.
 *
 * Sets the PrnFormat.has_minus flag on @prn,
 * indicating that the Unicode minus sign is
 * supported.
 */

void gretl_print_set_has_minus(Prn * prn)
{
    if (prn != null) {
        prn.format |= PrnFormat.has_minus;
    }
}

/*
 * gretl_print_has_minus:
 * @prn: printing struct.
 *
 * Returns: 1 if @prn supports Unicode minus sign, else 0.
 */

int gretl_print_has_minus(Prn * prn)
{
    if (prn != null) {
        return (prn.format & PrnFormat.has_minus) ? 1 : 0;
    } else {
        return 0;
    }
}

private int realloc_prn_buffer(Prn * prn, size_t newlen)
{
    char * tmp;
    int err = 0;

    if (newlen <= 0) {
        newlen = prn.bufsize * 2;
    }

    tmp = cast(char *) realloc(prn.buf, newlen);

    if (tmp == null) {
        err = 1;
    } else {
        prn.bufsize = newlen;
        prn.buf = tmp;
    }

    memset(prn.buf + prn.blen, 0, 1);

    return err;
}

private enum INIT_SIZE = 2048;
private enum MINREM = 1024;

private int pprintf_init(Prn * prn)
{
    int ret = 0;

    prn.bufsize = INIT_SIZE;
    prn.buf = cast(char *) malloc(prn.bufsize);

    if (prn.buf == null) {
        ret = -1;
    } else {
        memset(prn.buf, 0, 1);
        prn.blen = 0;
    }

    return ret;
}

private int get_nl_count(const(char) * s)
{
    int i, n = cast(int) (strlen(s) - 1);
    int nnl = 0;

    for (i=n; i>=0; i--) {
        if (s[i] == '\n') nnl++;
        else break;
    }

    return nnl;
}

/*
 * pprintf:
 * @prn: gretl printing struct.
 * @format: as in the C library's printf().
 * @Varargs: arguments to be printed.
 *
 * Multi-purpose printing function: can output to stream, to buffer
 * or to nowhere (silently discarding the output), depending on
 * how @prn was initialized.  Note that it's preferable to use
 * pputs() for large chunks of fixed text.
 *
 * Declared extern (C) because it is a C-style variadic: callers
 * pass C strings (use .ptr on D strings).
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

extern (C) int pprintf(Prn * prn, const(char) * format, ...)
{
    va_list args;
    int rem, plen = 0;

    if (prn == null || prn.fixed) {
        return 0;
    }

    if (format != null && *format != '\0') {
        prn.nlcount = cast(guint8) get_nl_count(format);
    }

    if (prn.fp != null) {
        /* printing to stream: straightforward */
        va_start(args, format);
        plen = vfprintf(prn.fp, format, args);
        va_end(args);
        return plen;
    } else if (prn.fz != null) {
        gchar * tmp = null;

        va_start(args, format);
        plen = g_vasprintf(&tmp, format, args);
        va_end(args);
        gzputs(prn.fz, tmp);
        g_free(tmp);
        return plen;
    }

    if (strncmp(format, "@init", 5) == 0) {
        return pprintf_init(prn);
    }

    if (prn.buf == null) {
        return 0;
    }

    rem = cast(int) (prn.bufsize - prn.blen);
    if (rem < MINREM) {
        if (realloc_prn_buffer(prn, 0)) {
            return -1;
        }
    }

    /* printing to buffer: be careful not to overrun */
    rem = cast(int) (prn.bufsize - prn.blen - 1);
    va_start(args, format);
    plen = vsnprintf(prn.buf + prn.blen, rem, format, args);
    va_end(args);

    if (plen >= rem) {
        /* buffer not big enough: try again */
        size_t newsize = prn.bufsize + plen + 1024;

        if (realloc_prn_buffer(prn, newsize)) {
            return -1;
        }
        rem = cast(int) (prn.bufsize - prn.blen - 1);
        va_start(args, format);
        plen = vsnprintf(prn.buf + prn.blen, rem, format, args);
        va_end(args);
    }

    if (plen > 0) {
        prn.blen += plen;
    }

    return plen;
}

/* void variant of pprintf() for use in certain callbacks */

extern (C) void pprintf2(Prn * prn, const(char) * format, ...)
{
    va_list args;
    int rem, plen = 0;

    if (prn == null || prn.fixed) {
        return;
    }

    if (prn.fp != null) {
        /* printing to stream: straightforward */
        va_start(args, format);
        vfprintf(prn.fp, format, args);
        va_end(args);
        return;
    }

    if (strncmp(format, "@init", 5) == 0) {
        pprintf_init(prn);
        return;
    }

    if (prn.buf == null) {
        return;
    }

    rem = cast(int) (prn.bufsize - prn.blen);
    if (rem < MINREM) {
        if (realloc_prn_buffer(prn, 0)) {
            return;
        }
    }

    rem = cast(int) (prn.bufsize - prn.blen - 1);
    va_start(args, format);
    plen = vsnprintf(prn.buf + prn.blen, rem, format, args);
    va_end(args);

    if (plen >= rem) {
        size_t newsize = prn.bufsize + plen + 1024;

        if (realloc_prn_buffer(prn, newsize)) {
            return;
        }
        rem = cast(int) (prn.bufsize - prn.blen - 1);
        va_start(args, format);
        plen = vsnprintf(prn.buf + prn.blen, rem, format, args);
        va_end(args);
    }

    if (plen > 0) {
        prn.blen += plen;
    }
}

/*
 * pputs:
 * @prn: gretl printing struct.
 * @s: constant string to print.
 *
 * Takes a C string (use .ptr on a D string literal).
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

extern (C) int pputs(Prn * prn, const(char) * s)
{
    int slen, rem;

    if (prn == null || prn.fixed) {
        return 0;
    }

    slen = cast(int) strlen(s);
    if (slen > 0) {
        prn.nlcount = cast(guint8) get_nl_count(s);
    }

    if (prn.fp != null) {
        fputs(s, prn.fp);
        return slen;
    } else if (prn.fz != null) {
        return gzputs(prn.fz, s);
    }

    if (prn.buf == null) {
        return 0;
    }

    rem = cast(int) (prn.bufsize - prn.blen);

    while (rem < MINREM || rem <= slen) {
        if (realloc_prn_buffer(prn, 0)) {
            return -1;
        }
        rem = cast(int) (prn.bufsize - prn.blen);
    }

    strcpy(prn.buf + prn.blen, s);
    prn.blen += slen;

    return slen;
}

/*
 * pputc:
 * @prn: gretl printing struct.
 * @c: character to print
 *
 * Returns: the number of bytes printed, or -1 on memory allocation
 * failure.
 */

extern (C) int pputc(Prn * prn, int c)
{
    if (prn == null || prn.fixed) {
        return 0;
    }

    if (c == '\n') prn.nlcount += 1;
    else if (c != '\0') prn.nlcount = 0;

    if (prn.fp != null) {
        fputc(c, prn.fp);
        return 1;
    } else if (prn.fz != null) {
        char[2] s = [cast(char) c, '\0'];

        return gzputs(prn.fz, s.ptr);
    }

    if (prn.buf == null) {
        return 0;
    }

    if (prn.bufsize - prn.blen < MINREM) {
        if (realloc_prn_buffer(prn, 0)) {
            return -1;
        }
    }

    prn.buf[prn.blen] = cast(char) c;
    prn.buf[prn.blen + 1] = '\0';
    prn.blen += 1;

    return 1;
}

/*
 * gretl_prn_newline:
 * @prn: gretl printing struct.
 *
 * Print a line break, in the mode appropriate to the
 * format of @prn (plain text, TeX or RTF).
 */

void gretl_prn_newline(Prn * prn)
{
    if (tex_format(prn)) {
        pputs(prn, "\\\\\n");
    } else if (rtf_format(prn)) {
        pputs(prn, "\\par\n");
    } else {
        pputc(prn, '\n');
    }
}

/*
 * gretl_print_ensure_vspace:
 * @prn: gretl printing struct.
 *
 */

void gretl_print_ensure_vspace(Prn * prn)
{
    if (prn != null && prn.nlcount < 2) {
        pputc(prn, '\n');
    }
}

/*
 * gretl_print_flush_stream:
 * @prn: gretl printing struct.
 *
 * If the output of @prn is directed to a stream, flush
 * that stream.
 */

void gretl_print_flush_stream(Prn * prn)
{
    if (prn != null) {
        if (prn.fp != null) {
            fflush(prn.fp);
        } else if (prn.fz != null) {
            gzflush(prn.fz, Z_SYNC_FLUSH);
        }
    }
}

/*
 * gretl_print_close_file:
 * @prn: gretl printing struct.
 *
 * If the output of @prn is directed to a file, close
 * the file. Also frees and sets to null the filename
 * associated with the stream (but does not remove the
 * file).
 */

void gretl_print_close_stream(Prn * prn)
{
    if (prn == null) return;

    if (prn.fp != null) {
        fclose(prn.fp);
        prn.fp = null;
        free(prn.fname);
        prn.fname = null;
    } else if (prn.fz != null) {
        gzclose(prn.fz);
        prn.fz = null;
        free(prn.fname);
        prn.fname = null;
    }
}

/*
 * printing_to_standard_stream:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the output of @prn is directed to one of the
 * standard streams, stdout or stderr, else 0.
 */

int printing_to_standard_stream(Prn * prn)
{
    int ret = 0;

    if (prn != null && (prn.fp == stdout || prn.fp == stderr)) {
        ret = 1;
    }

    return ret;
}

private void prn_push_stream(Prn * prn, FILE * fp, const(char) * fname,
                             const(char) * strvar)
{
    FpInfo fi = FpInfo(prn.fp, 0, null, null);

    if (prn.fplist == null) {
        prn.fplist = g_array_new(0, 0, cast(guint) FpInfo.sizeof);
    }
    fi.level = gretl_function_depth();
    if (fname != null) {
        fi.fname = g_strdup(fname);
    }
    if (strvar != null) {
        fi.strvar = g_strdup(strvar);
    }
    g_array_append_vals(prn.fplist, &fi, 1);

    prn.fp = fp;
}

private int handle_outbuf_content(FILE * fp, FpInfo * fi)
{
    UserVar * uv;
    int err = 0;

    /* get hold of the string variable specified in the
       initial call to "outfile"
    */
    uv = get_user_var_of_type_by_name(fi.strvar,
                                      GretlType.string_);
    if (uv == null) {
        err = GretlError.data;
    } else {
        /* grab the content of the tempfile */
        char * buf = null;
        long pos;
        size_t len;

        fflush(fp);
        fseek(fp, 0, SEEK_END);
        pos = ftell(fp);

        if (pos > 0) {
            buf = cast(char *) calloc(pos + 1, 1);
            if (buf == null) {
                err = GretlError.alloc;
            } else {
                fseek(fp, 0, SEEK_SET);
                len = fread(buf, 1, pos, fp);
                if (len < pos) {
                    err = GretlError.data;
                    free(buf);
                }
            }
        } else {
            buf = gretl_strdup("");
        }

        if (!err) {
            /* stick the content into the string variable */
            err = user_var_replace_value(uv, buf, GretlType.string_);
        }
    }

    return err;
}

private int prn_pop_stream(Prn * prn)
{
    FILE * prev = null;
    int err = 0;

    if (prn.fplist != null) {
        int n = prn.fplist.len;
        FpInfo * fi;

        if (n > 0) {
            fi = fplist_index(prn.fplist, n-1);
            prev = fi.fp;
            if (fi.strvar != null) {
                /* the --buffer case */
                err = handle_outbuf_content(prn.fp, fi);
                fclose(prn.fp);
                prn.fp = null;
                gretl_remove(fi.fname);
                g_free(fi.strvar);
                fi.strvar = null;
            }
            if (fi.fname != null) {
                g_free(fi.fname);
                fi.fname = null;
            }
            g_array_remove_index(prn.fplist, n-1);
        }
    }

    if (prn.fp != null && prn.fp != stdout && prn.fp != stderr) {
        fclose(prn.fp);
    }

    prn.fp = prev;

    return err;
}

/*
 * print_redirection_level:
 * @prn: gretl printing struct.
 *
 * Returns: 0 if the output of @prn has not been redirected
 * relative to its original setting, else the level of
 * (possibly nested) redirection.
 */

int print_redirection_level(Prn * prn)
{
    int ret = 0;

    if (prn != null) {
        if (prn.fp != null && prn.buf != null) {
            ret = 1;
        } else if (prn.fixed) {
            ret = 1;
        }
        if (prn_fp_list_active(prn)) {
            ret += prn.fplist.len; /* is this right? */
        }
    }

    return ret;
}

/*
 * print_redirection_filename:
 * @prn: gretl printing struct.
 *
 * Returns: the name of the file to which output is currently
 * redirected, if applicable, otherwise null.
 */

const(char) * print_redirection_filename(Prn * prn)
{
    if (prn != null && prn_fp_list_active(prn)) {
        int i = prn.fplist.len - 1;
        FpInfo * fi;

        fi = fplist_index(prn.fplist, i);
        if (fi != null && fi.fname != null) {
            return fi.fname;
        }
    }

    return null;
}

int print_redirected_at_level(Prn * prn, int level)
{
    if (prn.fplist != null) {
        FpInfo * fi;
        int i;

        for (i=0; i<prn.fplist.len; i++) {
            fi = fplist_index(prn.fplist, i);
            if (fi.level == level) {
                return 1;
            }
        }
    }

    return 0;
}

/*
 * print_start_redirection:
 * @prn: gretl printing struct.
 * @fp: stream to which output should be redirected.
 * @fname: name of the file or null.
 * @strvar: name of associated string variable, or null.
 *
 * Redirects output of @prn to @fp.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_start_redirection(Prn * prn, FILE * fp,
                            const(char) * fname,
                            const(char) * strvar)
{
    if (prn == null || prn.fixed) {
        return 1;
    }

    if (prn.fp != null) {
        /* flush the current stream */
        fflush(prn.fp);
    }
    if (fp == null) {
        /* "redirection" means just disable printing */
        prn.fixed = 1;
    } else {
        /* record current stream in prn->fplist, and
           hook output to specified stream
        */
        prn_push_stream(prn, fp, fname, strvar);
    }

    return 0;
}

/*
 * print_end_redirection:
 * @prn: gretl printing struct.
 *
 * If the output of @prn has been diverted relative to
 * its original setting, terminate the redirection.
 *
 * Returns: 0 on success, 1 on error.
 */

int print_end_redirection(Prn * prn)
{
    int err = 0;

    if (prn != null) {
        if (prn.fixed) {
            prn.fixed = 0;
        } else if (prn.fp != null) {
            err = prn_pop_stream(prn);
        }
    } else {
        err = GretlError.data;
    }

    return err;
}

/*
 * plain_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is plain text, else 0.
 */

int plain_format(Prn * prn)
{
    return (prn != null && (prn.format & PrnFormat.txt));
}

/*
 * rtf_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is RTF, else 0.
 */

int rtf_format(Prn * prn)
{
    if (prn != null) {
        if (prn.format & (PrnFormat.rtf | PrnFormat.rtf_txt)) {
            return 1;
        }
    }

    return 0;
}

/*
 * rtf_doc_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is RTF, and the
 * RTF preamble should be printed, else 0.
 */

int rtf_doc_format(Prn * prn)
{
    return (prn != null &&
            (prn.format & PrnFormat.rtf) &&
            (prn.format & PrnFormat.doc));
}

/*
 * tex_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is TeX, else 0.
 */

int tex_format(Prn * prn)
{
    return (prn != null && (prn.format & PrnFormat.tex));
}

/*
 * tex_doc_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is TeX and a complete
 * TeX document is being produced, else 0.
 */

int tex_doc_format(Prn * prn)
{
    return (prn != null &&
            (prn.format & PrnFormat.tex) &&
            (prn.format & PrnFormat.doc));
}

/*
 * tex_eqn_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is TeX and a model
 * is being represented as an equation, rather than in
 * tabular form, else 0.
 */

int tex_eqn_format(Prn * prn)
{
    return (prn != null &&
            (prn.format & PrnFormat.tex) &&
            (prn.format & PrnFormat.eqn));
}

/*
 * csv_format:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if the format of @prn is CSV, else 0.
 */

int csv_format(Prn * prn)
{
    return (prn != null && (prn.format & PrnFormat.csv));
}

/*
 * prn_format:
 * @prn: gretl printing struct.
 *
 * Returns: The formatting flags for @prn.
 */

int prn_format(Prn * prn)
{
    return (prn != null) ? prn.format : 0;
}

/*
 * prn_delim:
 * @prn: gretl printing struct.
 *
 * Returns: The character to be used as delimiter for CSV.
 */

char prn_delim(Prn * prn)
{
    return (prn != null) ? prn.delim : ',';
}

/*
 * gretl_print_has_buffer:
 * @prn: gretl printing struct.
 *
 * Returns: 1 if @prn is using a buffer for printing,
 * otherwise 0 (e.g., if a file is being used).
 */

int gretl_print_has_buffer(Prn * prn)
{
    return (prn != null && prn.buf != null);
}

/*
 * gretl_print_alloc:
 * @prn: gretl printing struct.
 * @s: desired size in bytes.
 *
 * If @prn is using a buffer for printing, try to
 * ensure that at least @s bytes are available for
 * further printing.  May speed things up if we
 * know there is a large amount of text to be
 * printed.
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_print_alloc(Prn * prn, size_t s)
{
    size_t rem;
    int err = 0;

    if (prn == null || prn.buf == null || prn.fixed) {
        return GretlError.data;
    }

    rem = prn.bufsize - prn.blen;

    if (rem <= s) {
        size_t newlen = prn.blen + s + 1;
        char * tmp;

        tmp = cast(char *) realloc(prn.buf, newlen);
        if (tmp == null) {
            err = GretlError.alloc;
        } else {
            prn.buf = tmp;
            prn.bufsize = newlen;
            memset(prn.buf + prn.blen, 0, 1);
        }
    }

    return err;
}

/*
 * gretl_print_read_tempfile:
 * @prn: gretl printing struct.
 * @err: location to receive error code.
 *
 * Retrieves a copy of the content of @prn, if it takes the
 * form of a file that has been opened in "w+" or "a+" mode.
 *
 * Returns: allocated buffer, or null on failure.
 */

char * gretl_print_read_tempfile(Prn * prn, int * err)
{
    size_t chk, bsize;
    import core.stdc.config : c_long;
    c_long pos0;
    char * buf = null;

    if (prn == null || prn.fp == null) {
        fprintf(stderr, "gretl_print_read_tempfile: no temp file to read!\n");
        *err = GretlError.data;
        return null;
    }

    fflush(prn.fp);
    pos0 = ftell(prn.fp);
    bsize = pos0;
    buf = cast(char *) calloc(bsize + 1, 1);

    if (buf == null) {
        *err = GretlError.alloc;
    } else {
        rewind(prn.fp);
        chk = fread(buf, 1, bsize, prn.fp);
        if (chk != bsize) {
            fprintf(stderr, "gretl_print_read_tempfile: bsize=%d but chk=%d\n",
                    cast(int) bsize, cast(int) chk);
            *err = GretlError.data;
            free(buf);
            buf = null;
        }
        /* get back to where we were */
        fseek(prn.fp, pos0, SEEK_SET);
    }

    return buf;
}
