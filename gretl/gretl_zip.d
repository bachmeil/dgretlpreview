/* Port of gretl_zip.c (gretl_zip.h has only prototypes: no header module)

   The C has two interchangeable back ends: libgsf (USE_GSF) or the native "gretlzip" plugin
   (plugin/gretlzip.c plus plugin/zipunzip). Neither is ported; instead the four primitives
   (unzip, unzip a session file, make a zipfile from a tree, zip the data.xml/data.bin pair)
   are implemented directly with Phobos' std.zip, which is the one back end of this port:
   - zip entries whose names could escape the extraction directory (absolute paths, ".." components,
     drive-style names) are rejected, as in the libgsf back end;
   - the deflate level of gretl_zip_datafile() is not adjustable in std.zip: level 0 stores the
     files, any other level uses zlib's default;
   - modification times are not preserved on extraction;
   - gretl_make_zipfile() adds a directory entry for every directory, with @path as the top directory
     of the archive.
   The package-building functions (package_make_zipfile and helpers) are ported as in the C. */

module gretl.gretl_zip;


import core.stdc.stdio : fprintf, stderr;
import core.stdc.string : strcmp, strlen, strrchr, strcpy;
version (Posix) import core.sys.posix.sys.stat : S_IFDIR, stat, stat_t; else import gretl.winport : S_IFDIR, stat, stat_t;

import std.algorithm.sorting : sort;
import std.array : array, replace;
import std.file : DirEntry, SpanMode, dirEntries, exists, isDir, isFile, mkdirRecurse, read, write;
import std.string : fromStringz, toStringz;
import std.zip : ArchiveMember, CompressionMethod, ZipArchive;

import gretl.glib : gchar;
import gretl.gfileutils : g_build_filename, g_get_current_dir, g_path_is_absolute;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_paths : gretl_chdir, gretl_deltree, gretl_dotdir, gretl_mkdir, gretl_remove;
import gretl.gretl_prn : Prn, pprintf, pputs;
import gretl.gretl_utils : gretl_copy_file;
import gretl.gstrfuncs : g_strdup, g_strdup_printf;
import gretl.libgretlh;
import gretl.strutils : gretl_strndup, has_suffix, strrslash;
import gretl.strutilsh : SLASH;

/* the error reporting of the C's handle_zip_error(), for a D exception message */

private int zip_error(const(char) * fname, string msg)
{
    fprintf(stderr, "handle_zip_error: '%.*s'\n", cast(int) msg.length, msg.ptr);
    gretl_errmsg_sprintf("%s: %.*s".ptr, fname, cast(int) msg.length, msg.ptr);
    return 1;
}

/* Reject zip/tree entry names that could be used to escape the
   intended extraction directory ("zip-slip"): absolute paths,
   any ".." path component, or a Windows drive-relative path such
   as "C:foo" (no separator required after the colon).
*/

private bool entry_name_is_unsafe(string name)
{
    if (name.length == 0) {
        return true;
    }

    /* absolute path, Unix- or Windows-style */
    if (name[0] == '/' || name[0] == '\\') {
        return true;
    }

    /* Windows drive-relative, e.g. "C:foo" or "C:\foo" */
    if (name.length > 1 && name[1] == ':') {
        return true;
    }

    /* any ".." path component */
    foreach (i; 0 .. name.length) {
        if (name[i] == '.' && i + 1 < name.length && name[i+1] == '.' &&
            (i == 0 || name[i-1] == '/' || name[i-1] == '\\') &&
            (i + 2 == name.length || name[i+2] == '/' || name[i+2] == '\\')) {
            return true;
        }
    }

    return false;
}

/* the name of the top-level directory of a gretl session file: the directory that holds
   "session.xml", falling back on the first path component of the first entry */

private string zip_topdir(ZipArchive zip)
{
    string first;

    foreach (name, am; zip.directory) {
        if (first.length == 0) {
            first = name;
        }
        if (name.length > 13 && name[$ - 11 .. $] == "session.xml") {
            string top = name[0 .. $ - 11];

            if (top.length && (top[$ - 1] == '/' || top[$ - 1] == '\\')) {
                top = top[0 .. $ - 1];
            }
            return top;
        }
    }

    foreach (i; 0 .. first.length) {
        if (first[i] == '/') {
            return first[0 .. i];
        }
    }
    return first;
}

/* @fname: the name of the file to be unzipped.
   @path: if non-NULL, the location into which the unzipping
   should be done (if NULL, use pwd).
   @zdirname: if non-NULL, set its content to the name of
   the top-level directory within the zipfile (via g_strdup).
*/

private int real_unzip(const(char) * fname, const(char) * path, gchar ** zdirname)
{
    try {
        auto zip = new ZipArchive(read(fromStringz(fname)));
        string base = path != null ? cast(string) fromStringz(path).idup : "";

        if (zdirname != null) {
            *zdirname = g_strdup(toStringz(zip_topdir(zip)));
        }
        if (base.length) {
            mkdirRecurse(base);
        }

        foreach (name, am; zip.directory) {
            if (entry_name_is_unsafe(name)) {
                gretl_errmsg_sprintf(_("Invalid archive entry name '%s'").ptr, toStringz(name));
                return 1;
            }
            string target = base.length ? base ~ "/" ~ name : name;

            if (name[$ - 1] == '/') {
                mkdirRecurse(target);
                continue;
            }
            auto slash = target.length;

            while (slash > 0 && target[slash - 1] != '/') {
                slash--;
            }
            if (slash > 1) {
                mkdirRecurse(target[0 .. slash - 1]);
            }
            zip.expand(am);
            write(target, am.expandedData);
        }
    } catch (Exception e) {
        return zip_error(fname, e.msg);
    }

    return 0;
}

/* a tree of files to a zipfile: entries are named "<path>/<relative name>" */

private int real_make_zipfile(const(char) * fname, const(char) * path)
{
    try {
        auto zip = new ZipArchive();
        string root = cast(string) fromStringz(path).idup;

        auto top = new ArchiveMember();

        top.name = root ~ "/";
        top.compressionMethod = CompressionMethod.none;
        zip.addMember(top);

        auto entries = dirEntries(root, SpanMode.breadth).array;

        entries.sort!((a, b) => a.name < b.name);
        foreach (e; entries) {
            auto m = new ArchiveMember();

            /* zip entry names always use '/' */
            string ename = e.name.replace("\\", "/");

            if (e.isDir) {
                m.name = ename ~ "/";
                m.compressionMethod = CompressionMethod.none;
            } else {
                m.name = ename;
                m.expandedData = cast(ubyte[]) read(e.name);
                m.compressionMethod = CompressionMethod.deflate;
            }
            zip.addMember(m);
        }
        write(fromStringz(fname), zip.build());
    } catch (Exception e) {
        return zip_error(fname, e.msg);
    }

    return 0;
}

/* @fname: zipfile to create; @path: directory holding data.xml and data.bin (both required, as with
   the libgsf back end); the entries are named plain "data.xml" and "data.bin" */

private int real_zip_datafile(const(char) * fname, const(char) * path, int level)
{
    try {
        auto zip = new ZipArchive();
        string dir = cast(string) fromStringz(path).idup;

        foreach (nm; ["data.xml", "data.bin"]) {
            auto m = new ArchiveMember();

            m.name = nm;
            m.expandedData = cast(ubyte[]) read(dir ~ "/" ~ nm);
            m.compressionMethod = level == 0 ? CompressionMethod.none : CompressionMethod.deflate;
            zip.addMember(m);
        }
        write(fromStringz(fname), zip.build());
    } catch (Exception e) {
        return zip_error(fname, e.msg);
    }

    return 0;
}

/**
 * gretl_unzip:
 * @fname: name of the file to unzip.
 *
 * Unzips @fname in the current directory, preserving any
 * internal directory structure.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip(const(char) * fname)
{
    return real_unzip(fname, null, null);
}

/**
 * gretl_unzip_into:
 * @fname: name of the file to unzip.
 * @dirname: the name of the directory in which unzipping
 * should take place.
 *
 * Unzips @fname in the specified directory, preserving any
 * internal directory structure.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip_into(const(char) * fname, const(char) * dirname)
{
    return real_unzip(fname, dirname, null);
}

/**
 * gretl_make zipfile:
 * @fname: name of the zip file to create.
 * @path: the path to the content which should be zipped.
 *
 * Creates a zip file of the specified name, whose content
 * is determined (recursively) by @path.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_make_zipfile(const(char) * fname, const(char) * path)
{
    return real_make_zipfile(fname, path);
}

/**
 * gretl_unzip_session_file:
 * @fname: name of the file to unzip.
 * @zdirname: location to receive the name of the top-level
 * directory within the zip archive.
 *
 * Specialized (slightly) unzipper for gretl session files.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_unzip_session_file(const(char) * fname, gchar ** zdirname)
{
    return real_unzip(fname, null, zdirname);
}

/**
 * gretl_zip_datafile:
 * @fname: name of the zip file to create.
 * @path: the path to the content which should be zipped.
 * @level: the zlib compression level to apply.
 *
 * Creates a zip file of the specified name, with content
 * given by @path, using the specified compression level.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_zip_datafile(const(char) * fname, const(char) * path,
                       int level)
{
    return real_zip_datafile(fname, path, level);
}

/* below: apparatus for making a zipfile for a function
   package
*/

private void zip_report(int err, int nf, Gretlopt opt, Prn * prn)
{
    if (opt & Gretlopt.OPT_G) {
        /* GUI use */
        if (err && nf) {
            pprintf(prn, "<@fail> (%s)\n", _("not found").ptr);
        } else if (err) {
            pputs(prn, "<@fail>\n");
        } else {
            pputs(prn, "<@ok>\n");
        }
    } else {
        if (err && nf) {
            pprintf(prn, _("failed (%s)\n").ptr, _("not found").ptr);
        } else if (err) {
            pputs(prn, _("failed\n").ptr);
        } else {
            pputs(prn, "OK\n");
        }
    }
}

private int pkg_zipfile_add(const(char) * fname,
                            const(char) * dotpath,
                            Gretlopt opt,
                            Prn * prn)
{
    gchar * dest = null;
    stat_t sbuf;
    int nf = 0;
    int err = 0;

    pprintf(prn, _("Copying %s... ").ptr, fname);

    if (stat(fname, &sbuf) != 0) {
        nf = err = GretlError.data;
    } else if (sbuf.st_mode & S_IFDIR) {
        /* aha, we've got a subdir */
        gchar * ziptmp;

        ziptmp = g_strdup_printf("%s%cpkgtmp.zip", dotpath, SLASH);
        err = gretl_make_zipfile(ziptmp, fname);
        if (!err) {
            err = gretl_unzip_into(ziptmp, dotpath);
            gretl_remove(ziptmp);
        }
        g_free(ziptmp);
    } else {
        /* a regular file, we hope; but if this is the PDF file
           for a package we'll allow that it may be found in a
           subdirectory (probably "doc")
        */
        const(char) * p = strrslash(cast(char *) fname);

        if (p != null && has_suffix(fname, ".pdf")) {
            dest = g_strdup_printf("%s%c%s", dotpath, SLASH, p + 1);
        } else {
            dest = g_strdup_printf("%s%c%s", dotpath, SLASH, fname);
        }
        err = gretl_copy_file(fname, dest);
    }

    zip_report(err, nf, opt, prn);
    g_free(dest);

    return err;
}

/**
 * package_make_zipfile:
 * @gfnname: name of the gfn file to be zip-packaged.
 * @pdfdoc: package has PDF documentation? (0/1).
 * @datafiles: names of any data files to include (or NULL).
 * @n_datafiles: the number of strings in @datafiles.
 * @pzipname: location to receive "dotdir" zipname, or NULL.
 * @dest: specific path for output zipfile, or NULL.
 * @opt: use OPT_G for GUI use (only);
 * @prn: gretl printer for progress.
 *
 * Collects the specified gfn file plus any additional files
 * it references (PDF doc and/or data files) and makes a
 * zip archive, using the user's "dotdir" as workspace.
 * If @pzipname is non-NULL this is taken as signal to
 * leave the zipfile where it has been created, and to
 * "return" its full path via this pointer. Otherwise, if
 * @dest is non-NULL it is taken to stipulate a path to
 * which the zipfile should be moved/renamed.
 *
 * Errors are flagged if both @pzipname and @dest are NULL,
 * if @gfnname does not have the ".gfn" extension, or if
 * the basename of the gfn file (minus extension) is over
 * 32 bytes long.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int package_make_zipfile(const(char) * gfnname,
                         int pdfdoc,
                         char ** datafiles,
                         int n_datafiles,
                         gchar ** pzipname,
                         const(char) * dest,
                         Gretlopt opt,
                         Prn * prn)
{
    char[FILENAME_MAX] pkgbase = 0;
    gchar * pkgname = null;
    gchar * origdir = null;
    gchar * tmp;
    gchar * dotpath = null;
    int len, dir_made = 0;
    int err = 0;

    if (pzipname == null && dest == null) {
        /* we need one or the other of these */
        return GretlError.data;
    }

    if (!has_suffix(gfnname, ".gfn")) {
        gretl_errmsg_set(_("Input must have extension \".gfn\"").ptr);
        return GretlError.data;
    }

    if (strlen(gfnname) >= pkgbase.sizeof) {
        gretl_errmsg_set(_("Path name is too long").ptr);
        return GretlError.data;
    }

    /* determine the common path to files for packaging,
       and also the name of the package
    */
    strcpy(pkgbase.ptr, gfnname);
    tmp = strrslash(pkgbase.ptr);
    if (tmp != null) {
        /* got some path component */
        *tmp = '\0';
        tmp++;
    } else {
        /* using current working directory */
        pkgbase[0] = '\0';
        tmp = cast(gchar *) gfnname;
    }

    len = cast(int) strlen(tmp) - 4; /* minus 4 for ".gfn" */
    if (len < 0 || len > 31) {
        gretl_errmsg_set(_("Invalid package name (31 bytes max)").ptr);
        return GretlError.data;
    }

    pkgname = gretl_strndup(tmp, len);

    /* record where we are now */
    origdir = g_get_current_dir();

    if (pkgbase[0] != '\0') {
        /* get into place for copying */
        pputs(prn, _("Getting in place... ").ptr);
        err = gretl_chdir(pkgbase.ptr);
        zip_report(err, 0, opt, prn);
    }

    if (!err) {
        /* path to temporary dir for zipping */
        dotpath = g_strdup_printf("%s%s", gretl_dotdir(), pkgname);
        pputs(prn, _("Making temporary directory... ").ptr);
        err = gretl_mkdir(dotpath);
        zip_report(err, 0, opt, prn);
    }

    if (!err) {
        dir_made = 1;
    }

    if (!err) {
        /* copy the gfn file into place */
        tmp = g_strdup_printf("%s.gfn", pkgname);
        err = pkg_zipfile_add(tmp, dotpath, opt, prn);
        g_free(tmp);
    }

    if (!err && pdfdoc) {
        /* copy PDF file into place */
        stat_t sbuf;

        tmp = g_strdup_printf("%s.pdf", pkgname);
        if (stat(tmp, &sbuf) != 0) {
            g_free(tmp);
            tmp = g_strdup_printf("doc/%s.pdf", pkgname);
        }
        err = pkg_zipfile_add(tmp, dotpath, opt, prn);
        g_free(tmp);
    }

    if (!err && datafiles != null) {
        /* copy data files into place, if any */
        int i;

        for (i=0; i<n_datafiles && !err; i++) {
            err = pkg_zipfile_add(datafiles[i],
                                  dotpath, opt, prn);
        }
    }

    if (!err) {
        /* get into place for making zipfile */
        err = gretl_chdir(gretl_dotdir());
        if (!err) {
            tmp = g_strdup_printf("%s.zip", pkgname);
            pprintf(prn, _("Making %s... ").ptr, tmp);
            err = gretl_make_zipfile(tmp, pkgname);
            zip_report(err, 0, opt, prn);
            if (!err) {
                if (pzipname != null) {
                    /* Signals that we should leave the zipfile in
                       the user's dotdir and return its name via
                       this pointer.
                    */
                    *pzipname = g_strdup_printf("%s%s.zip", gretl_dotdir(),
                                                pkgname);
                } else {
                    /* If we get here dest must be non-NULL,
                       but it may be a relative path, in which
                       case we need to prepend @origdir for the
                       zipfile copying operation.
                    */
                    const(char) * realdest = dest;
                    gchar * zipname = null;

                    pprintf(prn, _("Copying %s... ").ptr, tmp);

                    if (origdir != null && !g_path_is_absolute(dest)) {
                        zipname = g_build_filename(origdir, dest, null);
                        realdest = zipname;
                    }
                    err = gretl_copy_file(tmp, realdest);
                    zip_report(err, 0, opt, prn);
                    if (strcmp(tmp, realdest)) {
                        gretl_remove(tmp);
                    }
                    g_free(zipname);
                }
            }
            g_free(tmp);
        }
    }

    if (origdir != null) {
        /* get back to where we came from */
        gretl_chdir(origdir);
        g_free(origdir);
    }

    if (dir_made) {
        /* delete the temporary zip directory */
        gretl_deltree(dotpath);
    }

    g_free(dotpath);
    g_free(pkgname);

    return err;
}
