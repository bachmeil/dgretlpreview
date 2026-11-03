/* Port of gretl_paths.c (gretl_paths.h is in gretl_pathsh.d)

   Build-time choices of the C configuration, which the port fixes as follows (the Windows, macOS,
   MPI, R-library and libcurl branches are dropped, as everywhere in the port):
   - PLUGIN_SFX is "gretl-gtk3" (USE_GTK3);
   - GRETL_PREFIX is "/usr/local" (the configure default); LIBDIR and RLIBPATH are not defined;
   - HAVE_X12A and HAVE_TRAMO are not defined, so no x12arima/tramo working directories are made
     and their default program names are empty.
*/

module gretl.gretl_paths;


import core.stdc.errno : EACCES, errno;
import core.stdc.stdarg : va_arg, va_end, va_list, va_start;
import core.stdc.stdio : FILE, fclose, fgets, fopen, fprintf, perror, printf, remove, rename,
    sprintf, sscanf, stderr;
import core.stdc.stdlib : calloc, free, getenv, malloc, realloc;
version (Posix) import core.sys.posix.stdio : fdopen; else import gretl.winport : fdopen;
version (Posix) import core.sys.posix.stdlib : setenv; else import gretl.winport : setenv;
import core.stdc.string : memcpy, memset, strcat, strchr, strcmp, strlen, strncat, strncmp, strrchr, strstr,
    strtok, strcpy, strspn;
version (Posix) import core.sys.posix.fcntl : O_CREAT, open; else import gretl.winport : O_CREAT, open;
version (Posix) import core.sys.posix.sys.stat : S_IFDIR, S_IFREG, mode_t, stat, stat_t; else import gretl.winport : S_IFDIR, S_IFREG, mode_t, stat, stat_t;
version (Posix) import core.sys.posix.unistd : W_OK, access, chdir, getpid, readlink; else import gretl.winport : W_OK, access, chdir, getpid, readlink;
import etc.c.zlib : gzFile, gzclose, gzopen, gzread;
import std.conv : octal;
import std.string : fromStringz, toStringz;

import gretl.gerror : GError, g_error_free;
import gretl.gfileutils : GDir, g_build_filename, g_chdir, g_dir_close, g_dir_open, g_dir_read_name,
    g_dir_rewind, g_file_test, g_find_program_in_path, g_get_current_dir, g_mkdir_with_parents, g_mkstemp,
    g_path_get_dirname, g_path_is_absolute, g_remove, GFileTest;
import gretl.glib : gchar;
import gretl.gmem : g_free;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_set, gretl_errmsg_set_from_errno, gretl_errmsg_sprintf,
    gretl_strerror, set_gretl_alarm;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_pathsh : ConfigPaths, HelpPaths, PkgType, SearchType;
import gretl.gretl_prn : Prn, pprintf;
import gretl.gstrfuncs : g_strdup, g_strdup_printf, g_strlcat, g_strlcpy, g_strstrip, g_utf8_validate;
import gretl.libgretlh;
import gretl.libset : gretl_in_tool_mode, libset_set_bool, set_garch_alt_vcv, set_panel_hccme, set_tseries_hccme,
    set_xsect_hccme;
import gretl.libseth : SetKey;
import gretl.options : get_optval_string;
import gretl.strutils : gretl_strdup, gretl_strstrip, has_suffix, strings_array_add, strings_array_free, strrslash;
import gretl.strutilsh : SLASH, SLASHSTR;
import gretl.gretl_string_table : gretl_insert_builtin_string;
import gretl.gretl_intl : chinese_locale, set_lcnumeric;
import gretl.graphing : set_graph_color_from_string, set_plotstyle;
import gretl.gretl_bundle : gretl_bundle_new, gretl_bundle_set_int;

private enum PLUGIN_SFX = "gretl-gtk3";
private enum GRETL_PREFIX = "/usr/local";

/* local alias for ease of typing */
private enum PLEN = FILENAME_MAX;

private struct InternalPaths {
    char[PLEN] gretldir = 0;
    char[PLEN] dotdir = 0;
    char[PLEN] workdir = 0;
    char[PLEN] gnuplot = 0;
    char[PLEN] plotfile = 0;
    char[PLEN] plugpath = 0;
    char[PLEN+4] binbase = 0;
    char[PLEN] x12a = 0;
    char[PLEN] x12adir = 0;
    char[PLEN] tramo = 0;
    char[PLEN] tramodir = 0;
    char[PLEN] rbinpath = 0;
    char[PLEN] rlibpath = 0;
    char[PLEN] oxlpath = 0;
    char[PLEN] octpath = 0;
    char[PLEN] statapath = 0;
    char[PLEN] pypath = 0;
    char[PLEN] jlpath = 0;
    char[PLEN] lppath = 0;
    char[PLEN] mpiexec = 0;
    char[PLEN] mpi_hosts = 0;
    char[128] pngfont = 0;
}

private __gshared InternalPaths paths;

/* recorder for directories from which scripts were loaded: a GList in the C, here a plain array
   (appended at the end, searched from the end back to the start) */
private __gshared char *[] script_dirs;

private __gshared int force_en_cmdref = 0;
private __gshared int force_en_fnref = 0;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void set_gretl_tex_preamble();
    int gretl_is_updated(const(char) * s);
    int is_gretl_addon(const(char) * name);
    char * gretl_addon_get_path(const(char) * name);
    void update_addons_index(const(char) ** pkgs);
    const(char) * get_function_package_path_by_name(const(char) * name);
}

private enum LANG_AUTO = 0; /* first member of the language-id enum in gretl_intl.h */

private void set_helpfile_option(Gretlopt opt)
{
    if (opt & Gretlopt.OPT_N) {
        force_en_cmdref = 1;
        force_en_fnref = 1;
    }
}

private immutable string[5] helpfiles = [
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_cli_cmdref.pt -- otherwise leave it untranslated
    */
    N_("gretl_cli_cmdref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_gui_cmdref.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_cmdref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_commands_en.xml has been translated for your language,
       as in gretl_gui_help.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_help.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_functions_en.xml has been translated for your language,
       as in gretl_cli_fnref.pt -- otherwise leave it untranslated
    */
    N_("gretl_cli_fnref.en"),
    /* TRANSLATORS: you may change the two-letter language code
       if gretl_functions_en.xml has been translated for your language,
       as in gretl_gui_fnref.pt -- otherwise leave it untranslated
    */
    N_("gretl_gui_fnref.en")
];

private __gshared char[PLEN+19] hpath = 0;

const(char) * helpfile_path(int id, int cli, int en)
{
    int i = -1;

    hpath[0] = '\0';

    if ((id == HelpPaths.cmdref && force_en_cmdref) ||
        (id == HelpPaths.funcref && force_en_fnref)) {
        en = 1;
    }

    if (cli) {
        /* Command-line program */
        if (id == HelpPaths.cmdref) {
            i = 0;
        } else if (id == HelpPaths.funcref) {
            i = 3;
        }
    } else {
        /* GUI program */
        if (id == HelpPaths.cmdref) {
            i = 1;
        } else if (id == HelpPaths.gui_help) {
            i = 2;
        } else if (id == HelpPaths.funcref) {
            i = 4;
        }
    }

    if (i >= 0) {
        if (en || (_(helpfiles[i]).length != helpfiles[i].length)) {
            sprintf(hpath.ptr, "%s%s", paths.gretldir.ptr, helpfiles[i].ptr);
        } else {
            sprintf(hpath.ptr, "%s%s", paths.gretldir.ptr, _(helpfiles[i]).ptr);
        }
    }

    return hpath.ptr;
}

int using_translated_helpfile(int id)
{
    int ret = 0;
    int i = 0;

    if (id == HelpPaths.cmdref) {
        if (force_en_cmdref) return 0;
        i = 1;
    } else if (id == HelpPaths.funcref) {
        if (force_en_fnref) return 0;
        i = 4;
    } else {
        return 0;
    }

    /* If we're not forcing English help, the criterion is
       that the relevant help file has a "translated" filename
       and the translation can actually be opened.
    */

    if (helpfiles[i] != _(helpfiles[i])) {
        gchar * test;
        int err;

        test = g_strdup_printf("%s%s", paths.gretldir.ptr, _(helpfiles[i]).ptr);
        err = gretl_test_fopen(test, "r");
        if (err) {
            if (id == HelpPaths.cmdref) {
                force_en_cmdref = 1;
            } else {
                force_en_fnref = 1;
            }
        } else {
            ret = 1;
        }
        g_free(test);
    }

    return ret;
}

/* If @fname does not already have suffix @sfx, add it.
   With the qualification that if the @fname bears either of
   the standard gretl data-file suffixes, ".gdt" or ".gdtb",
   we won't stick the other one onto the end.
*/

private int maybe_add_suffix(char * fname, const(char) * sfx)
{
    if (has_suffix(fname, ".gdtb") && !strcmp(sfx, ".gdt")) {
        return 0;
    } else if (has_suffix(fname, ".gdt") && !strcmp(sfx, ".gdtb")) {
        return 0;
    } else if (!has_suffix(fname, sfx)) {
        strcat(fname, sfx);
        return 1;
    }

    return 0;
}

/**
 * utf8_encoded:
 * @s: the string to examine.
 *
 * The primary use of this function is to determine
 * whether a filename can be passed to regular C-library
 * functions on MS Windows: if it's in UTF-8 the answer
 * is No -- unless it's in the ASCII subset of UTF-8.
 *
 * Returns: non-zero if @s validates as UTF-8 and
 * contains bytes that are not printable ASCII,
 * otherwise zero.
 */

int utf8_encoded(const(char) * s)
{
    int ret = 0;

    if (g_utf8_validate(s, -1, null)) {
        const(ubyte) * p = cast(const(ubyte) *) s;

        while (*p) {
            if (*p < 32 || *p > 126) {
                /* not printable ASCII */
                ret = 1;
                break;
            }
            p++;
        }
    }

    return ret;
}

/**
 * gretl_fopen:
 * @fname: name of file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for the C library's fopen(), adding some error handling.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE * gretl_fopen(const(char) * fname, const(char) * mode)
{
    FILE * fp = null;

    gretl_error_clear();

    fp = fopen(fname, mode);

    if (errno != 0) {
        gretl_errmsg_set_from_errno(fname, errno);
    }

    return fp;
}

/**
 * gretl_mktemp:
 * @pattern: template for filename; must include "XXXXXX".
 * @mode: e.g. "w" for text use or "wb" for binary mode.
 *
 * A wrapper for the combination of mkstemp() and fdopen().
 * On successful exit @pattern holds the name of the newly
 * created file.
 *
 * Returns: file pointer, or %NULL on failure.
 */

FILE * gretl_mktemp(char * pattern, const(char) * mode)
{
    FILE * fp = null;
    int fd;

    gretl_error_clear();

    fd = g_mkstemp(pattern);

    if (errno != 0) {
        gretl_errmsg_set_from_errno(null, errno);
    } else if (fd != -1) {
        fp = fdopen(fd, mode);
    }

    return fp;
}

/**
 * gretl_test_fopen:
 * @fname: name of file to be opened.
 * @mode: mode as used with fopen().
 *
 * Attempts to open @fname in the given mode; if the opening
 * is successful the stream is then closed.
 *
 * Returns: 0 on success, -1 on filename encoding
 * failure, or the system errno on failure of fopen().
 */

int gretl_test_fopen(const(char) * fname, const(char) * mode)
{
    FILE * fp = null;
    int err = 0;

    gretl_error_clear();

    fp = fopen(fname, mode);

    if (fp == null) {
        err = errno;
    } else {
        fclose(fp);
        if (*mode == 'w') {
            remove(fname);
        }
    }

    return err;
}

/**
 * gretl_open:
 * @pathname: name of file to be opened.
 * @flags: flags to pass to the system open().
 * @mode: ignored unless @flags contains O_CREAT
 * or O_TMPFILE.
 *
 * A wrapper for the C library's open(), adding some error handling.
 *
 * Returns: new file descriptor, or -1 on error.
 */

int gretl_open(const(char) * pathname, int flags, int mode)
{
    mode_t m = 0;
    int fd = -1;

    gretl_error_clear();

    if (flags & O_CREAT) {
        m = cast(mode_t) mode;
    }

    fd = open(pathname, flags, m);

    if (errno != 0) {
        gretl_errmsg_set_from_errno(pathname, errno);
    }

    return fd;
}

/**
 * gretl_stat:
 * @fname: name of file to be examined.
 * @buf: pointer to a C struct stat (or NULL is OK if
 * the caller just wants the return value from the stat
 * call).
 *
 * A wrapper for the C library's stat().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_stat(const(char) * fname, stat_t * buf)
{
    stat_t tmp;

    gretl_error_clear();

    return stat(fname, buf == null ? &tmp : buf);
}

/**
 * gretl_file_exists:
 * @fname: name of file to be examined.
 *
 * Uses the C library's stat() function.
 *
 * Returns: 1 if @fname is the name of an existing file,
 * otherwise 0.
 */

int gretl_file_exists(const(char) * fname)
{
    return gretl_stat(fname, null) == 0;
}

/**
 * gretl_rename:
 * @oldpath: name of file to be opened.
 * @newpath: new name to give the file.
 *
 * A wrapper for the C library's rename().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_rename(const(char) * oldpath, const(char) * newpath)
{
    int err = 0;

    if (!strcmp(oldpath, newpath)) {
        /* check for no-op */
        return 0;
    }

    gretl_error_clear();

    err = rename(oldpath, newpath);

    if (errno != 0) {
        perror("rename");
        gretl_errmsg_set_from_errno("gretl_rename", errno);
    }

    return err;
}

/**
 * gretl_remove:
 * @path: name of file or directory to remove.
 *
 * A wrapper for remove().
 *
 * Returns: 0 on sucess, non-zero on failure.
 */

int gretl_remove(const(char) * path)
{
    return remove(path);
}

/**
 * gretl_gzopen:
 * @fname: name of gzipped file to be opened.
 * @mode: mode in which to open the file.
 *
 * A wrapper for zlib's gzopen().
 *
 * Returns: pointer to gzip stream, or %NULL on failure.
 */

gzFile gretl_gzopen(const(char) * fname, const(char) * mode)
{
    gzFile fz = null;

    gretl_error_clear();

    fz = gzopen(fname, mode);

    if (errno != 0) {
        gretl_errmsg_set_from_errno("gzopen", errno);
    }

    return fz;
}

/**
 * gretl_chdir:
 * @path: name of directory.
 *
 * A wrapper for POSIX chdir().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_chdir(const(char) * path)
{
    int err = 0;

    gretl_error_clear();

    err = chdir(path);

    if (errno != 0) {
        gretl_errmsg_set_from_errno("chdir", errno);
    }

    return err;
}

/**
 * gretl_isdir:
 * @path: path to check.
 *
 * A test for whether or not @path is the name of a directory.
 *
 * Returns: 1 if @path is the name of a directory, else 0.
 */

int gretl_isdir(const(char) * path)
{
    return g_file_test(path, GFileTest.is_dir);
}

/**
 * gretl_mkdir:
 * @path: name of directory to be created.
 *
 * Calls the underlying library function to create the
 * specified directory with mode 0755.  If the directory in
 * question already exists, this does not count as an error.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_mkdir(const(char) * path)
{
    int err;

    errno = 0;
    err = g_mkdir_with_parents(path, octal!755);

    if (err) {
        fprintf(stderr, "%s: %s\n", path, gretl_strerror(errno));
        err = 1;
    }

    return err;
}

private int real_delete_recursive(const(char) * path)
{
    GDir * dir;
    int err = 0;

    errno = 0;
    dir = g_dir_open(path, 0, null);

    if (dir == null) {
        err = 1;
    } else {
        const(gchar) * fname;

        err = g_chdir(path);
        while ((fname = g_dir_read_name(dir)) != null && !err) {
            /* recursively delete dir's contents */
            if (strcmp(fname, ".") && strcmp(fname, "..")) {
                if (g_file_test(fname, GFileTest.is_dir)) {
                    err = real_delete_recursive(fname);
                } else {
                    err = g_remove(fname);
                }
            }
        }
        if (!err) {
            g_dir_close(dir);
            /* delete the directory itself */
            if (g_chdir("..") == 0) {
                err = g_remove(path);
            }
        }
    }

    if (err) {
        gretl_errmsg_set_from_errno(path, errno);
        err = GretlError.fopen;
    }

    return err;
}

/**
 * gretl_deltree:
 * @path: name of directory to be deleted.
 *
 * Carries out recursive deletion of the specified directory.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_deltree(const(char) * path)
{
    gchar * savedir = null;
    int err;

    savedir = g_get_current_dir();
    err = real_delete_recursive(path);

    if (savedir != null) {
        g_chdir(savedir);
        g_free(savedir);
    }

    return err;
}

GDir * gretl_opendir(const(char) * name)
{
    GError * error = null;
    GDir * dir;

    dir = g_dir_open(name, 0, &error);

    if (error != null) {
        gretl_errmsg_set(error.message);
        g_error_free(error);
    }

    return dir;
}

/**
 * gretl_setenv:
 * @name: name of variable to be set.
 * @value: value to set.
 *
 * Wrapper for setenv().
 *
 * Returns: 0 on success, non-zero on failure.
 */

int gretl_setenv(const(char) * name, const(char) * value)
{
    return setenv(name, value, 1);
}

/**
 * gretl_write_access:
 * @fname: name of file to test.
 *
 * Returns: 0 on success (meaning that the current user has
 * write access to @fname), non-zero on failure.
 */

int gretl_write_access(char * fname)
{
    int err;

    gretl_error_clear();

    err = access(fname, W_OK);
    if (errno != 0) {
        gretl_errmsg_set_from_errno(fname, errno);
    }

    return err;
}

/**
 * gretl_is_xml_file:
 * @fname: name of file to test.
 *
 * Returns: 1 if @fname appears to be a (possibly gzipped) XML file,
 * otherwise 0.
 */

int gretl_is_xml_file(const(char) * fname)
{
    gzFile fz;
    char[6] test = 0;
    int ret = 0;

    fz = gretl_gzopen(fname, "rb");
    if (fz != null) {
        if (gzread(fz, test.ptr, 5)) {
            test[5] = '\0';
            if (!strcmp(test.ptr, "<?xml")) ret = 1;
        }
        gzclose(fz);
    }

    gretl_error_clear();

    return ret;
}

/**
 * gretl_path_prepend:
 * @file: target filename.
 * @path: path to prepend.
 *
 * Creates a path string by prepending @path, plus an appropriate
 * separator if needed, to @file.  The result is written back into
 * @file: this variable is assumed to have storage for at least
 * #FILENAME_MAX (libgretl.h) characters.
 *
 * Returns: 0 on success, or 1 if the final path string would
 * exceed #FILENAME_MAX characters (including nul-termination).
 */

int gretl_path_prepend(char * file, const(char) * path)
{
    int n = cast(int) (strlen(file) + strlen(path) + 1);
    char[PLEN] temp = 0;

    if (n > PLEN) {
        return 1;
    }

    strcpy(temp.ptr, path);
    n = cast(int) strlen(temp.ptr);

    if (temp[n-1] != SLASH && n < PLEN - 1) {
        strcat(temp.ptr, SLASHSTR.ptr);
    }

    strcat(temp.ptr, file);
    strcpy(file, temp.ptr);

    return 0;
}

private enum PathFlags {
    add_gdt = 1 << 0,
    add_gfn = 1 << 1,
    subdirs = 1 << 2
}

private int try_open_file(char * targ, const(char) * finddir,
                          const(gchar) * dname, int flags)
{
    char[PLEN] tmp = 0;
    int err, found = 0;

    g_strlcpy(tmp.ptr, finddir, tmp.sizeof);
    g_strlcat(tmp.ptr, dname, tmp.sizeof);
    g_strlcat(tmp.ptr, SLASHSTR.ptr, tmp.sizeof);
    g_strlcat(tmp.ptr, targ, tmp.sizeof);

    err = gretl_test_fopen(tmp.ptr, "r");

    /* Only attempt to append a suffix (".gdt", or ".gdtb" via the
       extra "b") if there is guaranteed room left in @tmp: this
       guards against maybe_add_suffix()'s own unchecked strcat()
       running past the end of the fixed-size buffer when @tmp is
       already (close to) full.
    */
    if (err && (flags & PathFlags.add_gdt) && strlen(tmp.ptr) + 5 < tmp.sizeof) {
        if (maybe_add_suffix(tmp.ptr, ".gdt")) {
            err = gretl_test_fopen(tmp.ptr, "r");
            if (err) {
                /* try .gdtb also */
                g_strlcat(tmp.ptr, "b", tmp.sizeof);
                err = gretl_test_fopen(tmp.ptr, "r");
            }
        }
    }

    if (!err) {
        strcpy(targ, tmp.ptr);
        found = 1;
    }

    return found;
}

private void make_finddir(char * targ, const(char) * src)
{
    int n;

    g_strlcpy(targ, src, PLEN);
    n = cast(int) strlen(targ);

    if (n > 0 && targ[n-1] != SLASH) {
        g_strlcat(targ, SLASHSTR.ptr, PLEN);
    }
}

private int got_subdir(const(char) * topdir, const(gchar) * dname)
{
    int ret = 0;

    if (strcmp(dname, ".") && strcmp(dname, "..")) {
        char[PLEN] tmp = 0;
        GDir * sub;

        g_strlcpy(tmp.ptr, topdir, tmp.sizeof);
        g_strlcat(tmp.ptr, dname, tmp.sizeof);
        sub = gretl_opendir(tmp.ptr);
        if (sub != null) {
            g_dir_close(sub);
            ret = 1;
        }
    }

    return ret;
}

/* Try to find @fname in a first-level subdirectory of @topdir.
   Return 1 if found, otherwise 0.
*/

private int find_in_subdir(const(char) * topdir, char * fname, int flags)
{
    GDir * dir;
    const(gchar) * dname;
    char[PLEN] finddir = 0;
    int found = 0;

    /* make find target */
    make_finddir(finddir.ptr, topdir);

    dir = gretl_opendir(finddir.ptr);
    if (dir != null) {
        while (!found && (dname = g_dir_read_name(dir)) != null) {
            if (got_subdir(finddir.ptr, dname)) {
                found = try_open_file(fname, finddir.ptr, dname, flags);
            }
        }
        g_dir_close(dir);
    }

    return found;
}

private char * search_dir(char * fname, const(char) * topdir, int flags)
{
    char[PLEN] orig = 0;
    int err;

    strcpy(orig.ptr, fname);

    if (gretl_path_prepend(fname, topdir) == 0) {
        err = gretl_test_fopen(fname, "r");
        if (!err) {
            return fname;
        }
        if (flags & PathFlags.add_gdt) {
            if (maybe_add_suffix(fname, ".gdt")) {
                err = gretl_test_fopen(fname, "r");
                if (!err) {
                    return fname;
                }
            } else if (maybe_add_suffix(fname, ".gdtb")) {
                err = gretl_test_fopen(fname, "r");
                if (!err) {
                    return fname;
                }
            }
        } else if (flags & PathFlags.add_gfn) {
            if (maybe_add_suffix(fname, ".gfn")) {
                err = gretl_test_fopen(fname, "r");
                if (!err) {
                    return fname;
                }
            }
        }
        strcpy(fname, orig.ptr);
        if (flags & PathFlags.subdirs) {
            if (find_in_subdir(topdir, fname, flags)) {
                return fname;
            }
        }
    }

    return null;
}

/**
 * get_plausible_search_dirs:
 * @stype: DATA_SEARCH for data file packages, DB_SEARCH for
 * gretl databases, FUNCS_SEARCH for function packages, or
 * SCRIPT_SEARCH for hansl scripts.
 * @n_dirs: location to receive the number of directories.
 *
 * Returns: an array of plausible search paths, depending on the
 * @type of search. The array should be freed when you are done
 * with it, using strings_array_free().
 */

char ** get_plausible_search_dirs(SearchType stype, int * n_dirs)
{
    char ** dirs = null;
    const(char) * subdir;
    const(char) * wd;
    char[PLEN] dirname = 0;

    *n_dirs = 0;

    if (stype == SearchType.funcs_search) {
        /* for testing of gfns */
        char * forcepath = getenv("FORCE_GFN_PATH");

        if (forcepath != null) {
            strings_array_add(&dirs, n_dirs, forcepath);
        }
    }

    if (stype == SearchType.data_search) {
        subdir = "data";
    } else if (stype == SearchType.db_search) {
        subdir = "db";
    } else if (stype == SearchType.funcs_search) {
        subdir = "functions";
    } else if (stype == SearchType.script_search) {
        subdir = "scripts";
    } else {
        fprintf(stderr, "get_plausible_search_dir: no type specified\n");
        return null;
    }

    /* the user's dotdir */
    gretl_build_path(dirname.ptr, gretl_dotdir(), subdir, null);
    strings_array_add(&dirs, n_dirs, dirname.ptr);

    /* add system dir */
    gretl_build_path(dirname.ptr, gretl_home(), subdir, null);
    strings_array_add(&dirs, n_dirs, dirname.ptr);

    /* the user's working dir */
    gretl_build_path(dirname.ptr, gretl_workdir(), subdir, null);
    strings_array_add(&dirs, n_dirs, dirname.ptr);

    if (stype != SearchType.funcs_search) {
        /* working dir, no subdir */
        strcpy(dirname.ptr, gretl_workdir());
        strings_array_add(&dirs, n_dirs, dirname.ptr);
    }

    /* a legacy thing: some files may have been written to
       the "default workdir" in the past
    */
    wd = maybe_get_default_workdir();
    if (wd != null) {
        gretl_build_path(dirname.ptr, wd, subdir, null);
        strings_array_add(&dirs, n_dirs, dirname.ptr);
        if (stype != SearchType.funcs_search) {
            strcpy(dirname.ptr, wd);
            strings_array_add(&dirs, n_dirs, dirname.ptr);
        }
    }

    return dirs;
}

/* it's a dirent thing */
private enum NAME_MAX = 255;

/**
 * gretl_function_package_get_path:
 * @name: the name of the package to find, without the .gfn extension.
 * @type: %PKG_SUBDIR for a package that lives in its own subdirectory,
 * %PKG_TOPLEV for a package file not in a subdirectory, or %PKG_ALL
 * for a package that may be at either level.
 *
 * Searches a list of directories in which we might expect to find
 * function packages, and, if the package in question is found,
 * returns a newly allocated string holding the full path to
 * the package's gfn file. Public (system) directories are
 * searched first, then directories in the user's filespace.
 *
 * Returns: allocated path on success, otherwise NULL.
 */

char * gretl_function_package_get_path(const(char) * name,
                                       PkgType type)
{
    char * ret = null;
    char[FILENAME_MAX] path = 0;
    char ** dirs;
    int err, found = 0;
    int i, n_dirs;

    if (type == PkgType.all || type == PkgType.subdir) {
        if (is_gretl_addon(name)) {
            return gretl_addon_get_path(name);
        }
    }

    path[0] = '\0';
    dirs = get_plausible_search_dirs(SearchType.funcs_search, &n_dirs);

    for (i=0; i<n_dirs && !found; i++) {
        const(char) * fndir = dirs[i];
        const(char) * dname;
        char * p;
        char[NAME_MAX+1] test = 0;
        GDir * dir;

        if ((dir = gretl_opendir(fndir)) == null) {
            continue;
        }
        if (type != PkgType.toplev) {
            /* look preferentially for .gfn files in their own
               subdirectories */
            while ((dname = g_dir_read_name(dir)) != null && !found) {
                if (!strcmp(dname, name)) {
                    sprintf(path.ptr, "%s%c%s%c%s.gfn", fndir, SLASH, dname, SLASH, dname);
                    err = gretl_test_fopen(path.ptr, "r");
                    if (!err) {
                        found = 1;
                    } else {
                        path[0] = '\0';
                    }
                }
            }
        }
        if (!found && type != PkgType.subdir) {
            /* look for .gfn files in the top-level functions
               directory */
            g_dir_rewind(dir);
            while ((dname = g_dir_read_name(dir)) != null && !found) {
                if (has_suffix(dname, ".gfn")) {
                    strcpy(test.ptr, dname);
                    p = strrchr(test.ptr, '.');
                    *p = '\0';
                    if (!strcmp(test.ptr, name)) {
                        sprintf(path.ptr, "%s%c%s", fndir, SLASH, dname);
                        found = 1;
                    }
                }
            }
        }
        g_dir_close(dir);
    }

    strings_array_free(dirs, n_dirs);

    if (path[0] != '\0') {
        ret = gretl_strdup(path.ptr);
    }

    return ret;
}

/* Search for file with basename @fname in directory
   @dirname, or in any subdirectory of same up to
   depth @maxdepth. If the file is found, write its
   path into @fullname and return 1, otherwise
   return 0 (= not found).
*/

private int find_file_in_dir(const(char) * fname,
                             const(char) * dirname,
                             char * fullname,
                             int maxdepth,
                             int depth)
{
    char[FILENAME_MAX] tmp = 0;
    const(gchar) * dname;
    stat_t sbuf;
    GDir * dir;
    int found = 0;

    dir = gretl_opendir(dirname);

    if (dir == null) {
        return 0;
    }

    /* look for top-level plain file first */
    while (!found && (dname = g_dir_read_name(dir)) != null) {
        if (!strcmp(dname, ".") ||
            !strcmp(dname, "..")) {
            continue;
        }
        sprintf(tmp.ptr, "%s%c%s", dirname, SLASH, dname);
        if (gretl_stat(tmp.ptr, &sbuf) < 0) {
            continue;
        } else if ((sbuf.st_mode & S_IFREG) &&
                   !strcmp(dname, fname)) {
            strcpy(fullname, tmp.ptr);
            found = 1;
        }
    }

    if (!found && depth < maxdepth) {
        /* then look in subdirs */
        g_dir_rewind(dir);
        depth++;
        while (!found && (dname = g_dir_read_name(dir)) != null) {
            if (!strcmp(dname, ".") ||
                !strcmp(dname, "..")) {
                continue;
            }
            sprintf(tmp.ptr, "%s%c%s", dirname, SLASH, dname);
            if (gretl_stat(tmp.ptr, &sbuf) < 0) {
                continue;
            } else if (sbuf.st_mode & S_IFDIR) {
                found = find_file_in_dir(fname, tmp.ptr, fullname,
                                         maxdepth, depth);
            }
        }
    }

    g_dir_close(dir);

    return found;
}

/**
 * get_package_data_path:
 * @ci: index of the active command.
 * @fname: the basename of the file whose full path is wanted.
 * @fullname: location to which the full path should be written
 * (should be at least FILENAME_MAX bytes).
 *
 * Looks for @fname in association with the name of a function
 * package, which must have been set previously using the
 * --frompkg option with the "open" command.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int get_package_data_path(int ci, const(char) * fname, char * fullname)
{
    const(char) * pkgname = null;
    int err = 0;

    *fullname = '\0';
    if (ci == GretlCmdIndex.open || ci == GretlCmdIndex.append) {
        pkgname = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_K);
    } else if (ci == GretlCmdIndex.join) {
        pkgname = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_R);
    }

    if (pkgname == null) {
        err = GretlError.data;
    } else {
        const(char) * ppath;
        char * gfnpath;

        ppath = get_function_package_path_by_name(pkgname);

        if (ppath != null) {
            gfnpath = gretl_strdup(ppath);
        } else {
            gfnpath = gretl_addon_get_path(pkgname);
        }

        if (gfnpath == null) {
            gretl_errmsg_sprintf(_("Couldn't find package %s").ptr,
                                 pkgname);
            err = GretlError.data;
        } else {
            char * p = strrslash(gfnpath);
            const(char) * needle;

            if (p != null) {
                *p = '\0';
            }

            /* trim path from @fname if present */
            needle = strrslash(cast(char *) fname);
            if (needle != null) {
                needle++;
            } else {
                needle = fname;
            }

            if (!find_file_in_dir(needle, gfnpath, fullname, 2, 0)) {
                gretl_errmsg_sprintf(_("Couldn't find file %s for package %s").ptr,
                                     needle, pkgname);
                *fullname = '\0';
                err = GretlError.fopen;
            }
            free(gfnpath);
        }
    }

    return err;
}

/**
 * gretl_addpath:
 * @fname: on input, the initially given file name; on output
 * a path may be prepended and/or a suffix may be appended.
 * This variable must be of size at least #FILENAME_MAX bytes
 * (as defined in libgretl.h) to allow for possible additions.
 * @script: if non-zero, assume the file we're looking for
 * is a hansl script.
 *
 * Elementary path-searching: try adding various paths to the given
 * @fname and see if it can be opened. Usually called by get_full_filename().
 * If @fname does not already have a dot-extension we may also try adding
 * an appropriate gretl extension in case no file is found.
 *
 * Returns: the path to the file that was found (in @fname), or
 * NULL if no file could be found even allowing for prepending
 * a path and/or adding a suffix.
 */

char * gretl_addpath(char * fname, int script)
{
    char[FILENAME_MAX] orig = 0;
    char * test;
    int found = 0;
    int err;

    if (fname == null || *fname == '\0' ||
        strlen(fname) >= orig.sizeof) {
        return null;
    }

    /* keep a backup of the original input */
    strcpy(orig.ptr, fname);

    if (g_path_is_absolute(fname)) {
        err = gretl_test_fopen(fname, "r");
        if (err && !script && maybe_add_suffix(fname, ".gdt")) {
            err = gretl_test_fopen(fname, "r");
            if (err) {
                strcpy(fname, orig.ptr);
            }
        }
        return err ? null : fname;
    }

    /* try workdir first */
    gretl_build_path(fname, paths.workdir.ptr, orig.ptr, null);
    err = gretl_test_fopen(fname, "r");
    if (!err) {
        /* got it */
        return fname;
    }

    if (script_dirs.length > 0) {
        int flags = script ? 0 : PathFlags.add_gdt;
        const(char) * gpath;

        /* the directories are tried from the most recently recorded back to the first */
        for (ptrdiff_t k = cast(ptrdiff_t) script_dirs.length - 1; k >= 0 && !found; k--) {
            strcpy(fname, orig.ptr);
            gpath = script_dirs[k];
            test = search_dir(fname, gpath, flags);
            if (test != null) {
                found = 1;
            }
        }
        if (found) {
            return fname;
        }
    }

    if (!found) {
        const(char) * gpath = gretl_home();
        char[PLEN] trydir = 0;

        strcpy(fname, orig.ptr);

        if (*gpath != '\0') {
            /* try searching some standard gretl paths */
            if (script) {
                sprintf(trydir.ptr, "%sscripts", gpath);
                test = search_dir(fname, trydir.ptr, PathFlags.subdirs);
                if (test != null) {
                    return fname;
                } else {
                    strcpy(fname, orig.ptr);
                    sprintf(trydir.ptr, "%sfunctions", gpath);
                    test = search_dir(fname, trydir.ptr, PathFlags.add_gfn | PathFlags.subdirs);
                    if (test != null) {
                        return fname;
                    }
                }
            } else if (has_suffix(fname, ".bin")) {
                /* database? */
                sprintf(trydir.ptr, "%sdb", gpath);
                test = search_dir(fname, trydir.ptr, 0);
                if (test != null) {
                    return fname;
                }
            } else {
                /* data file */
                sprintf(trydir.ptr, "%sdata", gpath);
                test = search_dir(fname, trydir.ptr, PathFlags.add_gdt | PathFlags.subdirs);
                if (test != null) {
                    return fname;
                }
            }
        }

        strcpy(fname, orig.ptr);

        /* now try user's personal filespace */
        gpath = gretl_dotdir();

        if (*gpath != '\0') {
            /* try looking in ~/Library or dotdir */
            if (script) {
                sprintf(trydir.ptr, "%sscripts", gpath);
                test = search_dir(fname, trydir.ptr, PathFlags.subdirs);
                if (test != null) {
                    return fname;
                } else {
                    strcpy(fname, orig.ptr);
                    sprintf(trydir.ptr, "%sfunctions", gpath);
                    test = search_dir(fname, trydir.ptr, PathFlags.add_gfn | PathFlags.subdirs);
                    if (test != null) {
                        return fname;
                    }
                }
            } else if (has_suffix(fname, ".bin")) {
                /* database? */
                sprintf(trydir.ptr, "%sdb", gpath);
                test = search_dir(fname, trydir.ptr, 0);
                if (test != null) {
                    return fname;
                }
            } else {
                /* data file? */
                sprintf(trydir.ptr, "%sdata", gpath);
                test = search_dir(fname, trydir.ptr, PathFlags.add_gdt | PathFlags.subdirs);
                if (test != null) {
                    return fname;
                }
            }
        }

        strcpy(fname, orig.ptr);
        gpath = gretl_workdir();

        if (*gpath != '\0') {
            /* try looking in user's dir (and subdirs) */
            test = search_dir(fname, gpath, PathFlags.subdirs);
            if (test != null) {
                return fname;
            }
        }

        strcpy(fname, orig.ptr);
        gpath = maybe_get_default_workdir();

        if (gpath != null && *gpath != '\0') {
            /* try looking in default workdir? */
            test = search_dir(fname, gpath, PathFlags.subdirs);
            if (test != null) {
                return fname;
            }
        }
    }

    strcpy(fname, orig.ptr);

    gretl_error_clear();

    return null;
}

/* It is assumed here that @fname starts with "~/" */

char * gretl_prepend_homedir(const(char) * fname, int * err)
{
    char * homedir = getenv("HOME");
    char * ret = null;

    if (homedir != null) {
        ret = cast(char *) malloc(strlen(homedir) + strlen(fname));
        if (ret == null) {
            *err = GretlError.alloc;
        } else {
            strcpy(ret, homedir);
            strcat(ret, fname + 1);
        }
    } else {
        *err = GretlError.data;
    }

    return ret;
}

private int substitute_homedir(char * fname)
{
    char * homedir = getenv("HOME");
    int err = 0;

    if (homedir != null) {
        int len = cast(int) strlen(fname);
        int homelen = cast(int) strlen(homedir);

        if (len + homelen > PLEN) {
            err = 1;
        } else {
            char[PLEN] tmp = 0;

            strcpy(tmp.ptr, homedir);
            strcat(tmp.ptr, fname + 1);
            strcpy(fname, tmp.ptr);
        }
    }

    return err;
}

private int get_gfn_special(char * fname)
{
    int ok = 0;

    if (!strchr(fname, '/') && !strchr(fname, '\\')) {
        /* no extra path elements */
        char * p;
        char[64] pkgname = 0;
        char * pkgpath;
        size_t len;

        /* Locate the ".gfn" suffix in the full, untruncated @fname
           first: searching for it only after copying into the
           fixed-size @pkgname buffer risks chopping the suffix off
           for names longer than the buffer, which would leave us
           with nothing to find and a NULL to dereference below.
        */
        p = strstr(fname, ".gfn");
        if (p != null && (len = p - fname) < pkgname.sizeof) {
            memcpy(pkgname.ptr, fname, len);
            pkgname[len] = '\0';
            pkgpath = gretl_function_package_get_path(pkgname.ptr, PkgType.all);
            if (pkgpath != null) {
                strcpy(fname, pkgpath);
                free(pkgpath);
                ok = 1;
            }
        }
    }

    return ok;
}

/**
 * get_full_filename:
 * @fname: input filename.
 * @fullname: filename to be filled out: must be at least #FILENAME_MAX bytes.
 * @opt: if OPT_S, treat as a script; if OPT_I we're responding
 * to the "include" command; if OPT_W, pass @fname through as is.
 *
 * Includes elementary path-searching: try adding various paths to the
 * given @fname, if appropriate, and see if it can be opened. For
 * internal gretl use.
 *
 * Returns: 0 on success, non-zero on error.
 */

int get_full_filename(const(char) * fname, char * fullname, Gretlopt opt)
{
    int script = (opt & (Gretlopt.OPT_S | Gretlopt.OPT_I))? 1 : 0;
    char * test = null;
    int err = 0;

    *fullname = '\0';

    if (fname == null || *fname == '\0') {
        return GretlError.data;
    }

    strncat(fullname, fname, PLEN - 1);

    if (opt & Gretlopt.OPT_W) {
        /* remote database: just use original name */
        return 0;
    }

    if (fullname[0] == '~' && fullname[1] == '/') {
        /* handle tilde == HOME */
        substitute_homedir(fullname);
    }

    if (g_path_is_absolute(fullname)) {
        goto test_open;
    }

    if (opt & Gretlopt.OPT_I) {
        /* respect special "include" setting if present */
        char * ipath = getenv("GRETL_INCLUDE_PATH");

        if (ipath != null && *ipath != '\0') {
            gretl_build_path(fullname, ipath, fname, null);
            if (gretl_test_fopen(fullname, "r") == 0) {
                return 0;
            } else {
                /* on failure, try ignoring @ipath? */
                *fullname = '\0';
                strncat(fullname, fname, PLEN - 1);
            }
        }
    }

    if (has_suffix(fullname, ".gfn") && get_gfn_special(fullname)) {
        /* checked for existence */
        return 0;
    }

    /* try a basic path search on this filename */
    test = gretl_addpath(fullname, script);

    gretl_normalize_path(fullname);

 test_open:

    if (!err && test == null) {
        err = gretl_test_fopen(fullname, "r");
        if (err) {
            /* ensure we return a gretl error code */
            err = GretlError.fopen;
        }
    }

    if (test != null && (opt & Gretlopt.OPT_S)) {
        /* If @test is non-NULL that means we actually found
           the file somewhere, so if it's a script we'll
           record the directory in which it was found.
        */
        gretl_set_script_dir(fullname);
    }

    return err;
}

int has_system_prefix(const(char) * fname, SearchType stype)
{
    const(char) * gretldir = gretl_home();
    int n = cast(int) strlen(gretldir);
    int ret = 0;

    if (strlen(fname) < n) {
        return 0;
    }

    if (!strncmp(fname, gretldir, n)) {
        if (stype == SearchType.data_search &&
            !strncmp(fname + n, "data", 4)) {
            ret = 1;
        } else if (stype == SearchType.script_search &&
                   !strncmp(fname + n, "scripts", 7)) {
            ret = 1;
        }
    }

    return ret;
}

private void set_gretl_plugpath(const(char) * path)
{
    paths.plugpath[0] = '\0';

    /* C: the branch without LIBDIR, see the note at the top of the file */
    const(char) * p = strstr(path, "/share");

    if (p != null) {
        /* back up a level */
        size_t len = p - path;

        strncat(paths.plugpath.ptr, path, len);
    } else {
        strcpy(paths.plugpath.ptr, path);
    }
    slash_terminate(paths.plugpath.ptr);
    strcat(paths.plugpath.ptr, "lib");
    slash_terminate(paths.plugpath.ptr);
    strcat(paths.plugpath.ptr, PLUGIN_SFX.ptr);

    slash_terminate(paths.plugpath.ptr);
}

private void set_gretl_binbase(const(char) * path)
{
    sprintf(paths.binbase.ptr, "%sdb", path);
}

/* This should be called after we're fairly confident that we
   have the dotdir setting right */

private int set_extra_dot_paths()
{
    char[PLEN+128] dirname = 0;
    int err = 0;

    /* the personal function package directory */
    dirname[0] = '\0';
    gretl_build_path(dirname.ptr, paths.dotdir.ptr, "functions".ptr, null);
    gretl_mkdir(dirname.ptr);

    paths.tramodir[0] = '\0';
    paths.x12adir[0] = '\0';

    /* C: "#if !defined(HAVE_TRAMO) && !defined(HAVE_X12A) return 0;" --
       neither is defined here, so the x12arima and tramo directories
       are not made */
    return err;
}

private void set_builtin_path_strings(int update)
{
    gretl_insert_builtin_string("gretldir", paths.gretldir.ptr);
    gretl_insert_builtin_string("gnuplot",  paths.gnuplot.ptr);
    gretl_insert_builtin_string("x12a",     paths.x12a.ptr);
    gretl_insert_builtin_string("tramo",    paths.tramo.ptr);

    if (!update) {
        /* these only have to be done once */
        gretl_insert_builtin_string("dotdir",   paths.dotdir.ptr);
        gretl_insert_builtin_string("workdir",  paths.workdir.ptr);
        gretl_insert_builtin_string("x12adir",  paths.x12adir.ptr);
        gretl_insert_builtin_string("tramodir", paths.tramodir.ptr);
    }

    if (paths.tramo[0] != '\0') {
        char[PLEN] s = 0;
        int n;

        s[0] = '\0';
        strncat(s.ptr, paths.tramo.ptr, PLEN - 1);
        n = cast(int) strlen(s.ptr);
        if (n >= 5 && !strcmp(s.ptr + n - 5, "tramo")) {
            strcpy(s.ptr + n - 5, "seats");
            gretl_insert_builtin_string("seats", s.ptr);
        }
    }
}

const(char) * gretl_home()
{
    if (paths.gretldir[0] == '\0') {
        /* try to remedy this! */
        char * chk = getenv("GRETL_HOME");

        if (chk != null) {
            strncat(paths.gretldir.ptr, chk, PLEN-2);
        }
    }

    if (paths.gretldir[0] != '\0') {
        slash_terminate(paths.gretldir.ptr);
    }

    return paths.gretldir.ptr;
}

private __gshared char[PLEN] bindir = 0;

const(char) * gretl_bindir()
{
    if (bindir[0] == '\0') {
        char * p;

        strcpy(bindir.ptr, paths.gretldir.ptr);
        p = strstr(bindir.ptr, "share/gretl");
        if (p != null) {
            *p = '\0';
            strcat(p, "bin/");
        }
    }

    return bindir.ptr;
}

private __gshared int plugin_path_set = 0;

const(char) * gretl_plugin_path()
{
    if (!plugin_path_set) {
        char * epath = getenv("GRETL_PLUGIN_PATH");

        if (epath != null) {
            paths.plugpath[0] = '\0';
            strncat(paths.plugpath.ptr, epath, PLEN - 2);
            slash_terminate(paths.plugpath.ptr);
        }

        /* if blank, try drawing on compiled-in values */
        if (paths.plugpath[0] == '\0') {
            strcat(paths.plugpath.ptr, GRETL_PREFIX.ptr);
            slash_terminate(paths.plugpath.ptr);
            strcat(paths.plugpath.ptr, "lib");
            slash_terminate(paths.plugpath.ptr);
            strcat(paths.plugpath.ptr, PLUGIN_SFX.ptr);
            slash_terminate(paths.plugpath.ptr);
        }
        plugin_path_set = 1;
    }

    return paths.plugpath.ptr;
}

const(char) * gretl_dotdir()
{
    return paths.dotdir.ptr;
}

gchar * gretl_make_dotpath(const(char) * basename)
{
    return g_build_filename(paths.dotdir.ptr, basename, null);
}

const(char) * gretl_workdir()
{
    return paths.workdir.ptr;
}

private __gshared char[PLEN] default_workdir = 0;

private const(char) * regular_default_workdir()
{
    char * home = getenv("HOME");
    const(char) * retval = null;

    if (home != null) {
        sprintf(default_workdir.ptr, "%s/gretl/", home);
        if (strcmp(default_workdir.ptr, paths.workdir.ptr)) {
            GDir * dir = gretl_opendir(default_workdir.ptr);

            if (dir != null) {
                g_dir_close(dir);
                retval = default_workdir.ptr;
            }
        }
    }

    return retval;
}

/**
 * maybe_get_default_workdir:
 *
 * Figures the full path to the default value of the
 * user's gretl working directory; call this @defdir.
 *
 * If this @defdir turns out to be the same as the
 * current gretl working directory, as would be returned
 * by gretl_workdir(), this function returns NULL,
 * otherwise it returns the @defdir value.
 *
 * Returns: a path, or NULL.
 */

const(char) * maybe_get_default_workdir()
{
    return regular_default_workdir();
}

private int validate_writedir(const(char) * dirname)
{
    int err = 0;

    if (*dirname == '\0') {
        gretl_errmsg_set(_("User directory is not set").ptr);
        return GretlError.data;
    }

    err = gretl_mkdir(dirname);
    if (err) {
        gretl_errmsg_sprintf(_("Couldn't create directory '%s'").ptr, dirname);
    }

    if (!err) {
        /* ensure the directory is writable */
        long mypid = cast(long) getpid();
        char[FILENAME_MAX] testname = 0;
        gchar * chkfile;
        FILE * fp;

        chkfile = g_strdup_printf("write.%ld.chk", mypid);
        gretl_build_path(testname.ptr, dirname, chkfile, null);
        fp = gretl_fopen(testname.ptr, "w");
        if (fp == null) {
            gretl_errmsg_sprintf((_("Couldn't write to '%s': ")
                                  ~ _("gretl will not work properly!")).ptr,
                                 dirname);
            err = GretlError.fopen;
        } else {
            fclose(fp);
            gretl_remove(testname.ptr);
        }
        g_free(chkfile);
    }

    if (err) {
        set_gretl_alarm(1);
    }

    return err;
}

private int set_gretl_workdir(const(char) * path)
{
    GDir * test;
    int err = 0;

    errno = 0;

    test = gretl_opendir(path);

    if (test == null) {
        gretl_errmsg_set_from_errno(path, errno);
        fprintf(stderr, "set_gretl_work_dir: '%s': failed\n", path);
        err = GretlError.fopen;
    } else {
        g_dir_close(test);
        strcpy(paths.workdir.ptr, path);
        slash_terminate(paths.workdir.ptr);
        gretl_insert_builtin_string("workdir", paths.workdir.ptr);
    }

    return err;
}

const(char) * gretl_gnuplot_path()
{
    return paths.gnuplot.ptr;
}

const(char) * gretl_plotfile()
{
    return paths.plotfile.ptr;
}

void report_plot_written(Prn * prn)
{
    if (prn != null) {
        pprintf(prn, _("wrote %s\n").ptr, paths.plotfile.ptr);
    }
}

const(char) * gretl_binbase()
{
    return paths.binbase.ptr;
}

const(char) * gretl_tramo()
{
    return paths.tramo.ptr;
}

const(char) * gretl_tramo_dir()
{
    return paths.tramodir.ptr;
}

const(char) * gretl_x12_arima()
{
    return paths.x12a.ptr;
}

const(char) * gretl_x12_arima_dir()
{
    return paths.x12adir.ptr;
}

int gretl_x12_is_x13()
{
    return strstr(paths.x12a.ptr, "x13") != null;
}

const(char) * gretl_rbin_path()
{
    return paths.rbinpath.ptr;
}

const(char) * gretl_rlib_path()
{
    return paths.rlibpath.ptr;
}

const(char) * gretl_oxl_path()
{
    return paths.oxlpath.ptr;
}

const(char) * gretl_octave_path()
{
    return paths.octpath.ptr;
}

const(char) * gretl_stata_path()
{
    return paths.statapath.ptr;
}

const(char) * gretl_python_path()
{
    return paths.pypath.ptr;
}

const(char) * gretl_julia_path()
{
    return paths.jlpath.ptr;
}

const(char) * gretl_lpsolve_path()
{
    return paths.lppath.ptr;
}

const(char) * gretl_mpi_hosts()
{
    return paths.mpi_hosts.ptr;
}

const(char) * gretl_mpiexec()
{
    return paths.mpiexec.ptr;
}

void gretl_set_script_dir(const(char) * s)
{
    gchar * add = g_path_get_dirname(s);
    ptrdiff_t k, found = -1;

    /* first matching record, as with g_list_find_custom() */
    for (k=0; k<cast(ptrdiff_t) script_dirs.length; k++) {
        if (strcmp(script_dirs[k], add) == 0) {
            found = k;
            break;
        }
    }

    if (found >= 0) {
        /* this directory is already in the list */
        if (found + 1 < cast(ptrdiff_t) script_dirs.length) {
            /* delete intervening record */
            g_free(script_dirs[found+1]);
            script_dirs = script_dirs[0 .. found+1] ~ script_dirs[found+2 .. $];
        }
        g_free(add);
    } else {
        script_dirs ~= add;
    }
}

void gretl_script_dirs_cleanup()
{
    foreach (d; script_dirs) {
        g_free(d);
    }
    script_dirs = null;
}

const(char) * gretl_png_font()
{
    return paths.pngfont.ptr;
}

void set_gretl_png_font(const(char) * s)
{
    strcpy(paths.pngfont.ptr, s);
}

void show_paths()
{
    printf(_("gretl: using these basic search paths:\n").ptr);
    printf("gretldir: %s\n", paths.gretldir.ptr);
    printf("workdir: %s\n", paths.workdir.ptr);
    printf("dotdir: %s\n", paths.dotdir.ptr);
    printf("gnuplot: %s\n", paths.gnuplot.ptr);
}

/* We have paths.gretldir in place: now test it by seeing if we can
   open the the GPL file "COPYING", which definitely should be in that
   directory.  If that doesn't work, try some remedial measures.
   Note, @config_path is the path garnered from the config file,
   which we may or may not have used to write paths.gretldir (and
   which may indeed be an empty string).
*/

private void check_gretldir(char * config_path)
{
    char[FILENAME_MAX] testname = 0;
    FILE * fp;
    int gotit = 0;

    gretl_build_path(testname.ptr, paths.gretldir.ptr, "COPYING".ptr, null);
    fp = gretl_fopen(testname.ptr, "r");

    if (fp != null) {
        /* should be fine as is */
        fclose(fp);
        gotit = 1;
    } else if (*config_path != '\0') {
        slash_terminate(config_path);
        if (strcmp(config_path, paths.gretldir.ptr)) {
            /* we weren't using the config-file version: try it now */
            sprintf(testname.ptr, "%sCOPYING", config_path);
            fp = gretl_fopen(testname.ptr, "r");
            if (fp != null) {
                strcpy(paths.gretldir.ptr, config_path);
                fclose(fp);
                gotit = 1;
            }
        }
    }

    if (!gotit && !gretl_in_tool_mode()) {
        /* we're messed up; try to recover */
        int pid = getpid();
        gchar * proc_exe;
        const(char) * s;
        ptrdiff_t nr;

        proc_exe = g_strdup_printf("/proc/%d/exe", pid);
        nr = readlink(proc_exe, testname.ptr, FILENAME_MAX - 1);

        if (nr > 0) {
            testname[nr] = '\0';
            fprintf(stderr, "gretl is process %d, '%s'\n", cast(int) pid, testname.ptr);
            /* should be something like /foo/bar/bin/gretl; we
               want the /foo/bar bit to append to
            */
            s = strstr(testname.ptr, "bin/gretl");
            if (s != null) {
                paths.gretldir[0] = '\0';
                strncat(paths.gretldir.ptr, testname.ptr, s - testname.ptr);
                strcat(paths.gretldir.ptr, "share/gretl/");
                fprintf(stderr, "gretldir is maybe '%s'?\n",
                        paths.gretldir.ptr);
            }
        }

        g_free(proc_exe);
    }
}

/* Called at start-up only: the @dirname argument is the value taken
   from the config file or registry.  In case we end up using a value
   other than the incoming one, sync back to @dirname.
*/

private void initialize_gretldir(char * dirname, Gretlopt opt)
{
    char * ghome = getenv("GRETL_HOME");
    int err = 0;

    if (ghome != null) {
        /* environment setting, if any, takes precedence */
        strcpy(paths.gretldir.ptr, ghome);
        slash_terminate(paths.gretldir.ptr);
    } else if (dirname != null && *dirname != '\0' &&
               paths.gretldir[0] == '\0') {
        /* use value from config/registry, unless we already got
           a value somehow */
        strcpy(paths.gretldir.ptr, dirname);
        slash_terminate(paths.gretldir.ptr);
    }

    if (paths.gretldir[0] == '\0') {
        /* use the compile-time value */
        strcpy(paths.gretldir.ptr, GRETL_PREFIX.ptr);
        strcat(paths.gretldir.ptr, "/share/gretl/");
    }

    check_gretldir(dirname);

    if (!err) {
        set_helpfile_option(opt);
        set_gretl_plugpath(paths.gretldir.ptr);
        set_gretl_binbase(paths.gretldir.ptr);
    }

    strcpy(dirname, paths.gretldir.ptr);
}

/**
 * set_gretl_plugin_path:
 * @path: path to the gretl plugins.
 *
 * For use by third-party code: the purpose of this function
 * is to ensure that libgretl can find its plugins.
 *
 * @prefix, if given, should be the path under which the plugins
 * are installed. On a unix-type system this might be, for example,
 * /usr/local/lib/gretl-gtk2; on MS Windows it might be
 * c:\program files\gretl\plugins.
 **/

void set_gretl_plugin_path(const(char) * path)
{
    if (path != null) {
        paths.plugpath[0] = '\0';
        strncat(paths.plugpath.ptr, path, PLEN - 2);
        slash_terminate(paths.plugpath.ptr);
    }
}

int gretl_set_path_by_name(const(char) * name, const(char) * path)
{
    char * targ = null;
    int builtin = 0;

    if (name == null || path == null) {
        return 1;
    } else if (!strcmp(name, "workdir")) {
        return set_gretl_workdir(path);
    } else if (!strcmp(name, "gretldir")) {
        targ = paths.gretldir.ptr;
        builtin = 1;
    } else if (!strcmp(name, "dotdir")) {
        targ = paths.dotdir.ptr;
        builtin = 1;
    } else if (!strcmp(name, "gnuplot")) {
        targ = paths.gnuplot.ptr;
        builtin = 1;
    } else if (!strcmp(name, "plotfile")) {
        targ = paths.plotfile.ptr;
    } else if (!strcmp(name, "rlibpath")) {
        targ = paths.rlibpath.ptr;
    } else if (!strcmp(name, "tramo")) {
        targ = paths.tramo.ptr;
        builtin = 1;
    } else if (!strcmp(name, "x12a")) {
        targ = paths.x12a.ptr;
        builtin = 1;
    } else {
        fprintf(stderr, "gretl_set_path_by_name: target '%s' not recognized\n",
                name);
        return 1;
    }

    if (targ != null) {
        *targ = '\0';
        strncat(targ, path, PLEN - 2);
        if (builtin) {
            gretl_insert_builtin_string(name, targ);
        }
    }

    return 0;
}

/* Called at start-up only: set the "hidden" working dir,
   which is not user-configurable.
*/

private int initialize_dotdir()
{
    char * dirname;
    int err = 0;

    paths.dotdir[0] = '\0';

    dirname = getenv("HOME");
    if (dirname != null) {
        sprintf(paths.dotdir.ptr, "%s/.gretl/", dirname);
    }

    err = validate_writedir(paths.dotdir.ptr);

    if (err) {
        paths.x12adir[0] = '\0';
        paths.tramodir[0] = '\0';
    } else {
        /* these paths depend on dotdir */
        err = set_extra_dot_paths();
    }

    return err;
}

private enum PathFlag {
    needs_slash = 1 << 0,
    blank_ok    = 1 << 1
}

/* Updating a gretl paths element: transcribe the new value unless it
   is unchanged; if it's a directory string that needs to be
   slash-terminated, check that; return 1 if any change was made to
   the internally recorded value, @targ, otherwise return 0.  Note
   that we ignore empty @src unless the PATH_BLANK_OK flag is given.
*/

private int maybe_transcribe_path(char * targ, char * src, int flags)
{
    int ret = 0;

    if (*targ != '\0' && g_path_is_absolute(targ) &&
        *src != '\0' && g_path_is_absolute(src)) {
        int st = gretl_stat(targ, null);
        int ss = gretl_stat(src, null);

        if (st == 0 && ss != 0) {
            /* don't replace OK path with broken? */
            return 0;
        }
    }

    if (*src == '\0' && (flags & PathFlag.blank_ok)) {
        if (*targ != '\0') {
            *targ = '\0';
            ret = 1;
        }
    } else if (*src != '\0') {
        if (flags & PathFlag.needs_slash) {
            slash_terminate(src);
        }
        if (strcmp(src, targ)) {
            strcpy(targ, src);
            ret = 1;
        }
    } else {
        /* back-sync */
        strcpy(src, targ);
    }

    return ret;
}

/* gretl_update_paths is called from the GUI preferences dialog. The
   internal path elements that can be set in this way are:

   gretldir
   gnuplot (but not for MS Windows package)
   tramo, x12a, rbinpath, rlibpath, oxlpath, octpath, statapath,
     pypath, jlpath, lppath

   * paths.workdir is updated via the separate working directory
     dialog

   * paths.pngfont is updated separately via the plot editing
     dialog

   The @opt argument can include OPT_N to force use of the English-
   language help file where this would not be the default.
*/

int gretl_update_paths(ConfigPaths * cpaths, Gretlopt opt)
{
    int ndelta = 0;

    if (maybe_transcribe_path(paths.gretldir.ptr, cpaths.gretldir.ptr,
                              PathFlag.needs_slash)) {
        set_helpfile_option(opt);
        set_gretl_plugpath(paths.gretldir.ptr);
        ndelta++;
    }

    /* gnuplot path */
    ndelta += maybe_transcribe_path(paths.gnuplot.ptr, cpaths.gnuplot.ptr, 0);

    /* other external programs */
    ndelta += maybe_transcribe_path(paths.x12a.ptr, cpaths.x12a.ptr, 0);
    ndelta += maybe_transcribe_path(paths.tramo.ptr, cpaths.tramo.ptr, 0);
    ndelta += maybe_transcribe_path(paths.rbinpath.ptr, cpaths.rbinpath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.oxlpath.ptr, cpaths.oxlpath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.octpath.ptr, cpaths.octpath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.statapath.ptr, cpaths.statapath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.pypath.ptr, cpaths.pypath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.jlpath.ptr, cpaths.jlpath.ptr, 0);
    ndelta += maybe_transcribe_path(paths.lppath.ptr, cpaths.lppath.ptr, 0);

    if (ndelta > 0) {
        /* we changed at least one thing that should be
           recorded in the builtin path strings */
        set_builtin_path_strings(1);
    }

    return 0;
}

/* unix-type variants of defaults for any paths that we need
   that were not found in the gretl config file.
*/

private void load_default_workdir(char * targ)
{
    char * home = getenv("HOME");

    if (home != null) {
        sprintf(targ, "%s/gretl/", home);
    } else {
        home = getenv("GRETL_WORKDIR");
        if (home != null) {
            strcpy(targ, home);
        } else {
            gretl_path_compose(targ, PLEN, paths.gretldir.ptr, "user/");
        }
    }
}

private void load_default_path(char * targ)
{
    static immutable string[4] app_paths = [
        "R",
        "octave",
        "oxl",
        "stata"
    ];

    if (targ == paths.workdir.ptr) {
        load_default_workdir(targ);
    } else if (targ == paths.gnuplot.ptr) {
        strcpy(targ, "gnuplot");
    } else if (targ == paths.x12a.ptr) {
        /* C: HAVE_X12A undefined here */
        *targ = '\0';
    } else if (targ == paths.tramo.ptr) {
        /* C: HAVE_TRAMO undefined here */
        *targ = '\0';
    } else if (targ == paths.rlibpath.ptr) {
        /* C: RLIBPATH undefined here */
        paths.rlibpath[0] = '\0';
    } else if (targ == paths.rbinpath.ptr) {
        strcpy(paths.rbinpath.ptr, app_paths[0].ptr);
    } else if (targ == paths.octpath.ptr) {
        strcpy(paths.octpath.ptr, app_paths[1].ptr);
    } else if (targ == paths.oxlpath.ptr) {
        strcpy(paths.oxlpath.ptr, app_paths[2].ptr);
    } else if (targ == paths.statapath.ptr) {
        strcpy(paths.statapath.ptr, app_paths[3].ptr);
    } else if (targ == paths.pypath.ptr) {
        strcpy(paths.pypath.ptr, "python");
    } else if (targ == paths.jlpath.ptr) {
        strcpy(paths.jlpath.ptr, "julia");
    } else if (targ == paths.lppath.ptr) {
        strcpy(paths.lppath.ptr, "liblpsolve55.so");
    } else if (targ == paths.mpiexec.ptr) {
        strcpy(paths.mpiexec.ptr, "mpiexec");
    } else if (targ == paths.mpi_hosts.ptr) {
        paths.mpi_hosts[0] = '\0';
    } else if (targ == paths.pngfont.ptr) {
        if (chinese_locale()) {
            strcpy(targ, "NSimSun 10");
        } else {
            strcpy(targ, "Vera 9");
        }
    }
}

int add_slash(char * s)
{
    if (s[strlen(s)-1] != SLASH) {
        strcat(s, SLASHSTR.ptr);
        return 1;
    }

    return 0;
}

private void path_init(char * targ, char * src, int needs_slash)
{
    if (*src) {
        strcpy(targ, src);
        if (needs_slash && slash_terminate(targ)) {
            strcpy(src, targ);
        }
    } else {
        load_default_path(targ);
        strcpy(src, targ);
    }
}

/* Set paths, falling back to defaults if no value has been supplied.
   We do this only at startup.  If the path that we record differs
   from that given in @cpaths, sync the value back to @cpaths
   (via path_init, above).
*/

private void copy_paths_with_fallback(ConfigPaths * cpaths)
{
    /* working directory */
    path_init(paths.workdir.ptr, cpaths.workdir.ptr, 1);

    /* gnuplot */
    path_init(paths.gnuplot.ptr, cpaths.gnuplot.ptr, 0);

    /* other external programs */
    path_init(paths.x12a.ptr, cpaths.x12a.ptr, 0);
    path_init(paths.tramo.ptr, cpaths.tramo.ptr, 0);
    path_init(paths.rbinpath.ptr, cpaths.rbinpath.ptr, 0);
    path_init(paths.rlibpath.ptr, cpaths.rlibpath.ptr, 0);
    path_init(paths.oxlpath.ptr, cpaths.oxlpath.ptr, 0);
    path_init(paths.octpath.ptr, cpaths.octpath.ptr, 0);
    path_init(paths.statapath.ptr, cpaths.statapath.ptr, 0);
    path_init(paths.pypath.ptr, cpaths.pypath.ptr, 0);
    path_init(paths.jlpath.ptr, cpaths.jlpath.ptr, 0);
    path_init(paths.lppath.ptr, cpaths.lppath.ptr, 0);
    path_init(paths.mpiexec.ptr, cpaths.mpiexec.ptr, 0);
    path_init(paths.mpi_hosts.ptr, cpaths.mpi_hosts.ptr, 0);

    /* graphing font */
    path_init(paths.pngfont.ptr, cpaths.pngfont.ptr, 0);
}

/* This is called after reading the gretl config file at startup
   (and only then).  Subsequent updates to paths via the GUI (if any)
   are handled by the function gretl_update_paths().

   The no_dotdir member of cpaths is used only when gretlcli is
   operating in "slave" mode (e.g. under a webserver). It forces gretl
   to use paths.workdir as the "dotdir" rather than using a directory
   under the executing user's HOME.  See
   http://gretl.sourceforge.net/slave/
*/

int gretl_set_paths(ConfigPaths * cpaths)
{
    int err0 = 0, err1 = 0;
    int retval = 0;

    paths.workdir[0] = '\0';
    paths.plotfile[0] = '\0';

    initialize_gretldir(cpaths.gretldir.ptr, Gretlopt.OPT_NONE);

    if (!cpaths.no_dotdir) {
        err0 = initialize_dotdir();
    }

    copy_paths_with_fallback(cpaths);

    if (cpaths.no_dotdir) {
        strcpy(paths.dotdir.ptr, paths.workdir.ptr);
    }

    if (strcmp(paths.dotdir.ptr, paths.workdir.ptr)) {
        err1 = validate_writedir(paths.workdir.ptr);
        if (err1) {
            /* try falling back on the default working dir */
            const(char) * defpath = maybe_get_default_workdir();

            if (defpath != null && *defpath != '\0' &&
                strcmp(defpath, paths.workdir.ptr)) {
                err1 = validate_writedir(defpath);
                if (err1 == 0) {
                    strcpy(paths.workdir.ptr, defpath);
                }
            }
        }
    }

    set_builtin_path_strings(0);
    set_gretl_tex_preamble();

    retval = (err0)? err0 : err1;

    return retval;
}

/* For writing a file, name given by user: if the path is not
   absolute, switch to the gretl "workdir" unless @fname begins
   with '~' in which case we switch to the user's HOME.
*/

const(char) * gretl_maybe_switch_dir(const(char) * fname)
{
    if (fname[0] == '~' && fname[1] == '/') {
        char * home = getenv("HOME");

        if (home != null && gretl_chdir(home) == 0) {
            fname += 2;
        }
    } else if (!g_path_is_absolute(fname)) {
        gretl_chdir(paths.workdir.ptr);
    }

    return fname;
}

/**
 * gretl_maybe_prepend_dir:
 * @fname: the original filename, which should be in a
 * location of length FILENAME_MAX.
 *
 * If @fname starts with the construction "~/" to indicate
 * the user's HOME, replace this with the full path to that
 * directory.  Otherwise, if @fname is not already an
 * absolute path, prepend the user's gretl working directory.
 * Otherwise do nothing.
 *
 * Returns: the possibly modified filename.
 */

char * gretl_maybe_prepend_dir(char * fname)
{
    char[FILENAME_MAX] tmp = 0;

    tmp[0] = '\0';

    if (fname[0] == '~' && fname[1] == '/') {
        char * home = getenv("HOME");

        if (home != null) {
            gretl_build_path(tmp.ptr, home, fname + 2, null);
        }
    } else if (!g_path_is_absolute(fname)) {
        gretl_build_path(tmp.ptr, paths.workdir.ptr, fname, null);
    }

    if (tmp[0] != '\0') {
        strcpy(fname, tmp.ptr);
    }

    return fname;
}

/**
 * gretl_read_user_file:
 * @fname: name of file to open.
 *
 * Attempts to open @fname in read-only mode.  If the file
 * is not found when the name is used "as is", we use
 * gretl_maybe_prepend_dir() to prepend the user's gretl
 * working directory and try again.
 *
 * Returns: file pointer, or NULL on failure.
 */

FILE * gretl_read_user_file(const(char) * fname)
{
    FILE * fp = gretl_fopen(fname, "r");

    if (fp == null) {
        char[FILENAME_MAX] fullname = 0;

        strcpy(fullname.ptr, fname);
        gretl_maybe_prepend_dir(fullname.ptr);
        if (fullname[0] != '\0') {
            fp = gretl_fopen(fullname.ptr, "r");
        }
    }

    return fp;
}

/* remove '.' and '..' from @path */

int gretl_normalize_path(char * path)
{
    char[FILENAME_MAX] tmp = 0;
    char[3] split = "/\0\0";
    char[2] slash = "/\0";
    char * pcpy, pbit, s = path;
    char ** S;
    char ** P = null;
    int i, n;
    int err = 0;

    if (*path == '\0') {
        return 0;
    }

    if (strstr(path, slash.ptr) == null) {
        return 0;
    }

    if (*path == '.') {
        /* absolutize the path first, if necessary */
        gchar * cwd = g_get_current_dir();

        if (cwd != null) {
            char * tmp2 = gretl_strdup(path + 1); /* C: inner "tmp" shadows the buffer */

            gretl_build_path(path, cwd, tmp2, null);
            free(tmp2);
            g_free(cwd);
        }
    }

    pcpy = gretl_strdup(path);
    if (pcpy == null) {
        return GretlError.alloc;
    }

    tmp[0] = '\0';
    s = pcpy;

    /* split string @s on the path separator and cumulate
       the pieces in array P, skipping any pieces which
       are just "." */

    n = 0;
    while ((pbit = strtok(s, split.ptr)) != null && !err) {
        if (strcmp(pbit, ".")) {
            S = cast(char **) realloc(P, (n+1) * (*P).sizeof);
            if (S == null) {
                err = GretlError.alloc;
            } else {
                P = S;
                P[n++] = pbit;
            }
        }
        s = null; /* for subsequent strtok calls */
    }

    if (!err) {
        int j;

        /* let each ".." annihilate the preceding path chunk */
        for (i=n-1; i>0; i--) {
            if (P[i] != null && !strcmp(P[i], "..")) {
                for (j=i-1; j>0; j--) {
                    if (P[j] != null && strcmp(P[j], "..")) {
                        P[j] = null;
                        break;
                    }
                }
            }
        }
        /* re-assemble the path */
        for (i=0; i<n; i++) {
            if (P[i] != null && strcmp(P[i], "..")) {
                strcat(tmp.ptr, slash.ptr);
                strcat(tmp.ptr, P[i]);
            }
        }
        strcpy(path, tmp.ptr);
    }

    free(P);
    free(pcpy);

    return err;
}

/**
 * slash_terminate:
 * @path: path string.
 *
 * Check whether @path is already slash-terminated, and if
 * not, append a #SLASH; @path should be a large enough
 * array to accept an extra byte.
 *
 * Returns: 1 if a slash was appended, otherwise 0.
 */

int slash_terminate(char * path)
{
    if (path != null && *path != '\0') {
        if (path[strlen(path) - 1] != '/') {
            strcat(path, "/");
            return 1;
        }
    }

    return 0;
}

private void rc_set_gp_extra_colors(const(char) * s)
{
    char[12][2] cstr = 0;

    cstr[0][0] = cstr[1][0] = '\0';

    if (sscanf(s, "%10s %10s", cstr[0].ptr, cstr[1].ptr) == 2) {
        set_graph_color_from_string(0, cstr[0].ptr);
        set_graph_color_from_string(1, cstr[1].ptr);
    }
}

private int rc_bool(const(char) * s)
{
    if (!strcmp(s, "true") || !strcmp(s, "1")) {
        return 1;
    } else {
        return 0;
    }
}

private void handle_use_cwd(int use_cwd, ConfigPaths * cpaths)
{
    libset_set_bool(SetKey.use_cwd, use_cwd);

    if (use_cwd) {
        gchar * cwd = g_get_current_dir();

        if (cwd != null) {
            cpaths.workdir[0] = '\0';
            strncat(cpaths.workdir.ptr, cwd, PLEN - 2);
            slash_terminate(cpaths.workdir.ptr);
            g_free(cwd);
        }
    }
}

private enum PROXLEN = 64;
private enum GRETLCLI_USE_CWD = 1;

/* called only on behalf of gretlcli (for all platforms) */

void get_gretl_config_from_file(FILE * fp, ConfigPaths * cpaths,
                                char * dbproxy, int * use_proxy,
                                int * updated, gchar ** gptheme)
{
    char[PLEN] line = 0;
    char[32] key = 0;
    char[PLEN] val = 0;

    while (fgets(line.ptr, cast(int) line.sizeof, fp) != null) {
        if (line[0] == '#') {
            continue;
        }
        if (!strncmp(line.ptr, "recent", 6)) {
            /* reached the "recent files" section */
            break;
        }
        if (sscanf(line.ptr, "%s", key.ptr) != 1) {
            continue;
        }
        val[0] = '\0';
        /* get the string that follows " = " */
        strncat(val.ptr, line.ptr + strlen(key.ptr) + 3, PLEN - 1);
        gretl_strstrip(val.ptr);
        if (!strcmp(key.ptr, "gretldir")) {
            cpaths.gretldir[0] = '\0';
            strncat(cpaths.gretldir.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "gnuplot")) {
            cpaths.gnuplot[0] = '\0';
            strncat(cpaths.gnuplot.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "workdir") || !strcmp(key.ptr, "userdir")) {
            /* "userdir" is a legacy thing */
            cpaths.workdir[0] = '\0';
            strncat(cpaths.workdir.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "no_dotdir")) {
            cpaths.no_dotdir = rc_bool(val.ptr);
        } else if (!strcmp(key.ptr, "shellok")) {
            libset_set_bool(SetKey.shell_ok, rc_bool(val.ptr));
        } else if (!strcmp(key.ptr, "usecwd")) {
            ; /* handled later (GRETLCLI_USE_CWD) */
        } else if (!strcmp(key.ptr, "lcnumeric")) {
            set_lcnumeric(LANG_AUTO, rc_bool(val.ptr));
        } else if (!strcmp(key.ptr, "dbproxy")) {
            *dbproxy = '\0';
            strncat(dbproxy, val.ptr, PROXLEN - 1);
        } else if (!strcmp(key.ptr, "useproxy")) {
            *use_proxy = rc_bool(val.ptr);
        } else if (!strcmp(key.ptr, "x12a")) {
            cpaths.x12a[0] = '\0';
            strncat(cpaths.x12a.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "tramo")) {
            cpaths.tramo[0] = '\0';
            strncat(cpaths.tramo.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "Rbin")) {
            cpaths.rbinpath[0] = '\0';
            strncat(cpaths.rbinpath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "Rlib")) {
            cpaths.rlibpath[0] = '\0';
            strncat(cpaths.rlibpath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "ox")) {
            cpaths.oxlpath[0] = '\0';
            strncat(cpaths.oxlpath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "octave")) {
            cpaths.octpath[0] = '\0';
            strncat(cpaths.octpath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "stata")) {
            cpaths.statapath[0] = '\0';
            strncat(cpaths.statapath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "python")) {
            cpaths.pypath[0] = '\0';
            strncat(cpaths.pypath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "julia")) {
            cpaths.jlpath[0] = '\0';
            strncat(cpaths.jlpath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "lpsolve")) {
            cpaths.lppath[0] = '\0';
            strncat(cpaths.lppath.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "mpiexec")) {
            cpaths.mpiexec[0] = '\0';
            strncat(cpaths.mpiexec.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "mpi_hosts")) {
            cpaths.mpi_hosts[0] = '\0';
            strncat(cpaths.mpi_hosts.ptr, val.ptr, PLEN - 1);
        } else if (!strcmp(key.ptr, "mpi_pref")) {
            ; /* C: set_mpi_variant() only with HAVE_MPI */
        } else if (!strcmp(key.ptr, "Png_font")) {
            cpaths.pngfont[0] = '\0';
            strncat(cpaths.pngfont.ptr, val.ptr, 128 - 1);
        } else if (!strcmp(key.ptr, "Gp_extra_colors")) {
            rc_set_gp_extra_colors(val.ptr);
        } else if (!strcmp(key.ptr, "HC_xsect")) {
            set_xsect_hccme(val.ptr);
        } else if (!strcmp(key.ptr, "HC_tseri")) {
            set_tseries_hccme(val.ptr);
        } else if (!strcmp(key.ptr, "HC_panel")) {
            set_panel_hccme(val.ptr);
        } else if (!strcmp(key.ptr, "HC_garch")) {
            set_garch_alt_vcv(val.ptr);
        } else if (!strcmp(key.ptr, "graph_theme")) {
            *gptheme = g_strdup(val.ptr);
        } else if (!strcmp(key.ptr, "build_date")) {
            *updated = gretl_is_updated(val.ptr);
        }
    }

    /* "workdir" is always the current directory */
    handle_use_cwd(1, cpaths);
}

void get_gretl_rc_path(char * rcfile)
{
    char * path = getenv("GRETL_CONFIG_FILE");

    if (path != null) {
        *rcfile = '\0';
        strncat(rcfile, path, FILENAME_MAX - 1);
    } else {
        path = getenv("HOME");
        if (path != null) {
            sprintf(rcfile, "%s/.gretl2rc", path);
        } else {
            strcpy(rcfile, ".gretl2rc");
        }
    }
}

/* non-Windows read of the gretl configuration file on behalf
   of the CLI program, gretlcli; the Windows variant of this,
   win32_cli_read_rc(), is in gretl_win32.c
*/

int cli_read_rc()
{
    ConfigPaths cpaths;
    char[FILENAME_MAX] rcfile = 0;
    char[PROXLEN] dbproxy = 0;
    gchar * gptheme = null;
    int use_proxy = 0;
    int updated = 0;
    FILE * fp;
    int err = 0;

    /* C: ConfigPaths cpaths = {0}: D's default for the char arrays is 0xFF, so clear them */
    memset(&cpaths, 0, ConfigPaths.sizeof);

    get_gretl_rc_path(rcfile.ptr);
    fp = gretl_fopen(rcfile.ptr, "r");

    if (fp == null) {
        err = GretlError.fopen;
    } else {
        get_gretl_config_from_file(fp, &cpaths, dbproxy.ptr,
                                   &use_proxy, &updated,
                                   &gptheme);
        fclose(fp);
    }

    if (err) {
        gretl_set_paths(&cpaths);
    } else {
        err = gretl_set_paths(&cpaths);
    }

    if (gptheme != null) {
        set_plotstyle(gptheme);
        g_free(gptheme);
    }

    if (updated) {
        update_addons_index(null);
    }

    return err;
}

private int dir_is_writable(const(char) * dirname)
{
    int ok = 0;

    if (gretl_mkdir(dirname) == 0) {
        gchar * test = g_strdup_printf("%s%c%s", dirname, SLASH, "wtest".ptr);

        ok = (gretl_test_fopen(test, "w") == 0);
        g_free(test);
    }

    return ok;
}

private int get_user_install_path(char * path, const(char) * subdir)
{
    const(char) * dirname = gretl_dotdir();
    int err = 0;

    if (dirname == null || *dirname == '\0') {
        err = GretlError.fopen;
    } else {
        sprintf(path, "%s%s", dirname, subdir);
        err = (dir_is_writable(path) == 0);
    }

    return err;
}

private int get_system_install_path(char * path, const(char) * subdir)
{
    sprintf(path, "%s%s", gretl_home(), subdir);

    if (dir_is_writable(path)) {
        return 0;
    } else {
        return GretlError.fopen;
    }
}

private __gshared char[FILENAME_MAX] install_path = 0;

/* get a path that's suitable for writing a function
   package or scripts package on installation
*/

const(char) * gretl_package_install_path(const(char) * payload)
{
    if (install_path[0] == '\0') {
        int sys_first = 1;
        int err = 0;

        if (sys_first) {
            err = get_system_install_path(install_path.ptr, payload);
        }
        if (err || !sys_first) {
            err = get_user_install_path(install_path.ptr, payload);
        }
        if (err) {
            install_path[0] = '\0';
        } else {
            slash_terminate(install_path.ptr);
        }
    }

    return install_path.ptr;
}

int gretl_path_compose(char * targ, int len,
                       const(char) * s1,
                       const(char) * s2)
{
    targ[0] = '\0';
    if (strlen(s1) + strlen(s2) >= len) {
        gretl_errmsg_set("filename is too long");
        return GretlError.data;
    } else {
        strcpy(targ, s1);
        strcat(targ, s2);
        return 0;
    }
}

private int strcat_bounded(char * dest, const(char) * src)
{
    size_t len = strlen(dest) + strlen(src);

    if (len >= FILENAME_MAX) {
        fprintf(stderr, "strcat_bounded: FILENAME_MAX would be exceeded\n");
        return 1;
    } else {
        strcat(dest, src);
        return 0;
    }
}

private int strncat_bounded(char * dest, const(char) * src, size_t n)
{
    size_t len = strlen(dest) + n;

    if (len >= FILENAME_MAX) {
        fprintf(stderr, "strncat_bounded: FILENAME_MAX would be exceeded\n");
        return 1;
    } else {
        strncat(dest, src, n);
        return 0;
    }
}

/* Code borrowed from GLib (gfileutils.c) and adapted to write to an
   input char * (@targ) instead of returning a newly allocated
   string. The code is also somewhat simplified by the assumption that
   if the platform is not MS Windows the directory separator is always
   '/': this is a safe assumption for the platforms supported by gretl.
*/

private void real_build_path(char * targ,
                             const(gchar) * first_element,
                             ref va_list args)
{
    bool is_first = true;
    bool have_leading = false;
    const(gchar) * single_element = null;
    const(gchar) * next_element;
    const(gchar) * last_trailing = null;
    int err = 0;

    next_element = first_element;

    while (!err) {
        const(gchar) * element;
        const(gchar) * start;
        const(gchar) * end;

        if (next_element) {
            element = next_element;
            next_element = va_arg!(const(gchar) *)(args);
        } else {
            break;
        }

        /* ignore empty elements */
        if (*element == '\0') {
            continue;
        }

        start = element;
        while (*start == '/') {
            start++;
        }

        end = start + strlen(start);
        while (end >= start + 1 && end[-1] == '/') {
            end--;
        }

        last_trailing = end;
        while (last_trailing >= element + 1 && last_trailing[-1] == '/') {
            last_trailing--;
        }

        if (!have_leading) {
            /* If the leading and trailing separator strings are in the
               same element and overlap, the result is exactly that
               element.
            */
            if (last_trailing <= start) {
                single_element = element;
            }
            err = strncat_bounded(targ, element, start - element);
            have_leading = true;
        } else {
            single_element = null;
        }

        if (err) {
            break;
        } else if (end == start) {
            continue;
        }

        if (!is_first) {
            err = strcat_bounded(targ, "/");
        }
        err = strncat_bounded(targ, start, end - start);
        is_first = false;
    }

    if (single_element) {
        *targ = '\0';
        strcat_bounded(targ, single_element);
    } else if (last_trailing) {
        strcat_bounded(targ, last_trailing);
    }
}

/**
 * gretl_build_path:
 * @targ: target string to write to. Must be pre-allocated on
 * the stack or heap with size >= FILENAME_MAX as defined
 * in libgretl.h.
 * @first_element: first component of path.
 *
 * Writes to @targ a path composed of @first_element
 * plus any additional string arguments supplied before the
 * required terminating NULL. An appropriate separator is inserted
 * between the components of the path.
 *
 * Returns: the target string, @targ.
 */

extern (C) char * gretl_build_path(char * targ, const(gchar) * first_element, ...)
{
    va_list args;

    *targ = '\0';

    va_start(args, first_element);
    real_build_path(targ, first_element, args);
    va_end(args);

    return targ;
}

private struct ForeignPath {
    const(char) * id;
    const(char) * path;
}

GretlBundle * foreign_info()
{
    /* C: a file-scope static table; here built on entry because it points into @paths */
    ForeignPath[8] fpaths = [
        ForeignPath("Rbin",    paths.rbinpath.ptr),
        ForeignPath("Rlib",    paths.rlibpath.ptr),
        ForeignPath("ox",      paths.oxlpath.ptr),
        ForeignPath("octave",  paths.octpath.ptr),
        ForeignPath("stata",   paths.statapath.ptr),
        ForeignPath("python",  paths.pypath.ptr),
        ForeignPath("julia",   paths.jlpath.ptr),
        ForeignPath(null, null)
];

    GretlBundle * b = gretl_bundle_new();
    gchar * fullpath;
    int found, i;
    int dbg = 0;

    char * s;
    if ((s = getenv("FOREIGN_DEBUG")) != null) {
        dbg = 1;
    }

    for (i=0; fpaths[i].id != null; i++) {
        if (dbg) {
            fprintf(stderr, "'%s' -> '%s'\n", fpaths[i].id, fpaths[i].path);
        }
        if (fpaths[i].path[0] == '\0') {
            gretl_bundle_set_int(b, fpaths[i].id, 0);
        } else if (g_path_is_absolute(fpaths[i].path)) {
            found = gretl_stat(fpaths[i].path, null) == 0;
            gretl_bundle_set_int(b, fpaths[i].id, found);
        } else {
            fullpath = g_find_program_in_path(fpaths[i].path);
            if (fullpath == null) {
                gretl_bundle_set_int(b, fpaths[i].id, 0);
            } else {
                gretl_bundle_set_int(b, fpaths[i].id, 1);
                g_free(fullpath);
            }
        }
    }

    return b;
}

/* next: determining the path for gretl package downloads
   (the GString of the C is a D string, NUL-terminated where passed to C) */

private int get_target_in_home(ref string gs, const(char) * dlname)
{
    const(char) * savedir = gretl_dotdir();
    string dl = cast(string) fromStringz(dlname).idup;
    int subdir = 1;
    int err = 0;

    if (savedir == null || *savedir == '\0') {
        return GretlError.fopen;
    }

    gs = cast(string) fromStringz(savedir).idup;

    if (strstr(dlname, ".gfn") || strstr(dlname, ".zip")) {
        gs ~= "functions";
    } else if (strstr(dlname, ".ggz")) {
        gs ~= "db";
    } else if (strstr(dlname, "addons")) {
        gs ~= "functions";
    } else if (strstr(dlname, ".tar.gz")) {
        gs ~= "data";
    } else {
        gs ~= dl;
        subdir = 0;
    }

    if (subdir) {
        err = gretl_mkdir(toStringz(gs));
        if (!err) {
            gs ~= SLASH;
            gs ~= dl;
        }
    }

    return err;
}

private void get_system_target(ref string gs, const(char) * dlname)
{
    gs = cast(string) fromStringz(gretl_home()).idup;

    if (strstr(dlname, ".ggz")) {
        gs ~= "db";
    } else if (strstr(dlname, "addons")) {
        gs ~= "functions";
    } else if (strstr(dlname, ".tar.gz")) {
        gs ~= "data";
    } else {
        gs ~= "functions";
    }

    gs ~= SLASH;
    gs ~= cast(string) fromStringz(dlname).idup;
}

/* Given the name of a file to be downloaded, @dlname, construct a
   suitable path (for which the user has write permission) into which
   the downloaded content should be written.
*/

gchar * get_download_path(const(char) * dlname, int * err)
{
    string gs;
    gchar * targ = null;
    int done_home = 0;

    get_system_target(gs, dlname);

    if (!*err) {
        *err = gretl_test_fopen(toStringz(gs), "w");
        if (*err == EACCES && !done_home) {
            /* permissions problem: write to home dir instead */
            *err = get_target_in_home(gs, dlname);
        }
    }

    if (!*err) {
        targ = g_strdup(toStringz(gs));
    }

    return targ;
}
