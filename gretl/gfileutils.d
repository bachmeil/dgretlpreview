/*
   Replacements for the glib file helpers used by the C code:
   g_file_get_contents(), g_file_set_contents() and g_get_current_dir(), written with std.file and
   getcwd so that the port does not depend on glib for them.

   - g_file_get_contents() reads a whole file into a malloc'd, NUL-terminated
     buffer (so it also works for text), returns its length in *length, and
     reports failure (missing file, permissions, ...) through *error, with a
     message. Release the buffer with g_free()/free().
   - g_get_current_dir() returns the working directory as a malloc'd string,
     or null on failure.

   - g_file_test(), g_mkdir_with_parents(), g_mkstemp(), g_remove(), g_chdir(),
     g_path_is_absolute(), g_path_get_dirname(), g_build_filename() (variadic, as in glib) and
     g_find_program_in_path() are small POSIX-only versions (the port drops the Windows branches).
   - GDir and g_dir_open/g_dir_read_name/g_dir_rewind/g_dir_close wrap opendir/readdir; as in glib,
     g_dir_read_name() skips "." and "..", and the returned name is only valid until the next read.

   Used by: dataio.d, gretl_paths.d.
*/

module gretl.gfileutils;


import core.stdc.errno : EEXIST, ENOTDIR, errno;
import core.stdc.stdarg : va_arg, va_end, va_list, va_start;
import core.stdc.stdio : remove;
import core.stdc.stdlib : getenv, malloc;
version (Posix) import core.sys.posix.stdlib : mkstemp; else import gretl.winport : mkstemp;
import core.stdc.string : memcpy, strcmp, strerror, strlen, strrchr;
version (Posix) import core.sys.posix.dirent : DIR, closedir, dirent, opendir, readdir, rewinddir; else import gretl.winport : DIR, closedir, dirent, opendir, readdir, rewinddir;
version (Posix) import core.sys.posix.sys.stat : S_ISDIR, S_ISLNK, S_ISREG, lstat, mkdir, stat, stat_t; else import gretl.winport : S_ISDIR, S_ISLNK, S_ISREG, lstat, mkdir, stat, stat_t;
version (Posix) import core.sys.posix.unistd : X_OK, access, chdir, getcwd; else import gretl.winport : X_OK, access, chdir, getcwd;

import std.file : read, write, FileException;

import gretl.gerror : GError, g_set_error_literal_msg;
import gretl.glib : gboolean, gchar, gint, gsize;
import gretl.gstrfuncs : g_strdup, g_strndup;

import std.algorithm.iteration : splitter;
import std.string : fromStringz;

gboolean g_file_get_contents(const(gchar) * filename, gchar ** contents, gsize * length,
                             GError ** error)
{
    try {
        ubyte[] data = cast(ubyte[]) read(fromStringz(filename));
        gchar * buf = cast(gchar *) malloc(data.length + 1);

        if (data.length > 0) {
            memcpy(buf, data.ptr, data.length);
        }
        buf[data.length] = 0;
        *contents = buf;
        if (length != null) {
            *length = data.length;
        }
        return true;
    } catch (FileException e) {
        g_set_error_literal_msg(error, e.msg);
    } catch (Exception e) {
        g_set_error_literal_msg(error, e.msg);
    }
    *contents = null;
    if (length != null) {
        *length = 0;
    }
    return false;
}

/* glib writes to a temporary file and renames it over @filename; this version writes @filename
   directly (so a failure part-way can leave a partial file). @length < 0 means NUL-terminated. */
gboolean g_file_set_contents(const(gchar) * filename, const(gchar) * contents, ptrdiff_t length,
                             GError ** error)
{
    try {
        size_t n = (length < 0) ? strlen(contents) : cast(size_t) length;

        write(fromStringz(filename), (cast(const(ubyte) *) contents)[0 .. n]);
        return true;
    } catch (FileException e) {
        g_set_error_literal_msg(error, e.msg);
    } catch (Exception e) {
        g_set_error_literal_msg(error, e.msg);
    }
    return false;
}

gchar * g_get_current_dir()
{
    char[4096] buf;

    if (getcwd(buf.ptr, buf.length) == null) {
        return null;
    }
    size_t n = 0;

    while (buf[n] != 0) {
        n++;
    }
    gchar * ret = cast(gchar *) malloc(n + 1);

    memcpy(ret, buf.ptr, n + 1);
    return ret;
}

enum GFileTest {
    is_regular    = 1 << 0,
    is_symlink    = 1 << 1,
    is_dir        = 1 << 2,
    is_executable = 1 << 3,
    exists        = 1 << 4
}

gboolean g_file_test(const(gchar) * filename, GFileTest test)
{
    stat_t st;

    if ((test & GFileTest.exists) && access(filename, 0) == 0) {
        return true;
    }
    if ((test & GFileTest.is_executable) && access(filename, X_OK) == 0) {
        /* glib: running as root, X_OK is true for any file: also require a non-directory */
        if (stat(filename, &st) == 0 && !S_ISDIR(st.st_mode)) {
            return true;
        }
    }
    if ((test & GFileTest.is_symlink) && lstat(filename, &st) == 0 && S_ISLNK(st.st_mode)) {
        return true;
    }
    if (test & (GFileTest.is_regular | GFileTest.is_dir)) {
        if (stat(filename, &st) == 0) {
            if ((test & GFileTest.is_regular) && S_ISREG(st.st_mode)) {
                return true;
            }
            if ((test & GFileTest.is_dir) && S_ISDIR(st.st_mode)) {
                return true;
            }
        }
    }
    return false;
}

gint g_mkdir_with_parents(const(gchar) * pathname, gint mode)
{
    char * fn = g_strdup(pathname);
    char * p = fn;

    if (pathname == null || *pathname == '\0') {
        errno = 22; /* EINVAL */
        return -1;
    }

    while (*p == '/') {
        p++;
    }

    gint ret = 0;

    do {
        while (*p && *p != '/') {
            p++;
        }
        if (!*p) {
            p = null;
        } else {
            *p = '\0';
        }
        if (!g_file_test(fn, GFileTest.exists)) {
            if (mkdir(fn, mode) == -1 && errno != EEXIST) {
                ret = -1;
                break;
            }
        } else if (!g_file_test(fn, GFileTest.is_dir)) {
            errno = ENOTDIR;
            ret = -1;
            break;
        }
        if (p) {
            *p++ = '/';
            while (*p && *p == '/') {
                p++;
            }
            if (!*p) {
                break;
            }
        }
    } while (p);

    import core.stdc.stdlib : free;
    free(fn);
    return ret;
}

gint g_mkstemp(gchar * tmpl)
{
    return mkstemp(tmpl);
}

gint g_remove(const(gchar) * filename)
{
    version (Windows) {
        /* the C runtime's remove() cannot delete a directory */
        import gretl.winport : rmdir;

        stat_t sb;

        if (stat(filename, &sb) == 0 && S_ISDIR(sb.st_mode)) {
            return rmdir(filename);
        }
    }
    return remove(filename);
}

gint g_chdir(const(gchar) * path)
{
    return chdir(path);
}

gboolean g_path_is_absolute(const(gchar) * file_name)
{
    version (Windows) {
        return file_name[0] == '/' || file_name[0] == '\\' ||
            (file_name[0] != 0 && file_name[1] == ':' && (file_name[2] == '/' || file_name[2] == '\\'));
    } else {
        return file_name[0] == '/';
    }
}

gchar * g_path_get_dirname(const(gchar) * file_name)
{
    const(gchar) * base = strrchr(file_name, '/');

    if (base == null) {
        return g_strdup(".");
    }
    while (base > file_name && *base == '/') {
        base--;
    }
    size_t len = base + 1 - file_name;
    gchar * ret = cast(gchar *) malloc(len + 1);

    memcpy(ret, file_name, len);
    ret[len] = 0;
    return ret;
}

/* g_build_filename: join the NULL-terminated list of path elements with "/" (collapsing repeated
   separators at the joins, ignoring empty elements), as glib does on POSIX systems. */

extern (C) gchar * g_build_filename(const(gchar) * first_element, ...)
{
    va_list args;
    const(gchar) *[32] elems;
    size_t n = 0, total = 1;

    va_start(args, first_element);
    const(gchar) * e = first_element;

    while (e != null && n < elems.length) {
        elems[n++] = e;
        total += strlen(e) + 1;
        e = va_arg!(const(gchar) *)(args);
    }
    va_end(args);

    gchar * buf = cast(gchar *) malloc(total + 1);
    size_t len = 0;
    bool is_first = true, have_leading = false;
    const(gchar) * single_element = null, last_trailing = null;

    buf[0] = 0;
    foreach (i; 0 .. n) {
        const(gchar) * element = elems[i];

        if (*element == 0) {
            continue;
        }
        const(gchar) * start = element;

        while (*start == '/') {
            start++;
        }
        const(gchar) * end = start + strlen(start);

        while (end >= start + 1 && end[-1] == '/') {
            end--;
        }
        last_trailing = end;
        while (last_trailing >= element + 1 && last_trailing[-1] == '/') {
            last_trailing--;
        }
        if (!have_leading) {
            if (last_trailing <= start) {
                single_element = element;
            }
            memcpy(buf + len, element, start - element);
            len += start - element;
            have_leading = true;
        } else {
            single_element = null;
        }
        if (end == start) {
            continue;
        }
        if (!is_first) {
            buf[len++] = '/';
        }
        memcpy(buf + len, start, end - start);
        len += end - start;
        is_first = false;
    }
    buf[len] = 0;

    if (single_element) {
        len = strlen(single_element);
        memcpy(buf, single_element, len + 1);
    } else if (last_trailing) {
        size_t tl = strlen(last_trailing);

        memcpy(buf + len, last_trailing, tl + 1);
    }
    return buf;
}

/* g_find_program_in_path: absolute names are tested directly; otherwise each directory in PATH
   (default "/bin:/usr/bin:." as in glib) is tried. Returns a malloc'd path or null. */

gchar * g_find_program_in_path(const(gchar) * program)
{
    if (g_path_is_absolute(program)) {
        return g_file_test(program, GFileTest.is_executable) &&
            !g_file_test(program, GFileTest.is_dir) ? g_strdup(program) : null;
    }
    const(char) * path = getenv("PATH");
    string p = path != null ? cast(string) fromStringz(path).idup : "/bin:/usr/bin:.";
    string prog = cast(string) fromStringz(program).idup;

    foreach (dir; p.splitter(':')) {
        string d = dir.length == 0 ? "." : dir;
        string cand = d ~ "/" ~ prog ~ "\0";

        if (g_file_test(cand.ptr, GFileTest.is_executable) && !g_file_test(cand.ptr, GFileTest.is_dir)) {
            return g_strdup(cand.ptr);
        }
    }
    return null;
}

struct GDir {
    DIR * dir;
}

GDir * g_dir_open(const(gchar) * path, uint flags, GError ** error)
{
    DIR * d = opendir(path);

    if (d == null) {
        g_set_error_literal_msg(error, cast(string) fromStringz(strerror(errno)).idup);
        return null;
    }
    GDir * ret = cast(GDir *) malloc(GDir.sizeof);

    ret.dir = d;
    return ret;
}

const(gchar) * g_dir_read_name(GDir * dir)
{
    dirent * entry;

    while ((entry = readdir(dir.dir)) != null) {
        if (strcmp(entry.d_name.ptr, ".") != 0 && strcmp(entry.d_name.ptr, "..") != 0) {
            return entry.d_name.ptr;
        }
    }
    return null;
}

void g_dir_rewind(GDir * dir)
{
    rewinddir(dir.dir);
}

void g_dir_close(GDir * dir)
{
    import core.stdc.stdlib : free;

    if (dir != null) {
        closedir(dir.dir);
        free(dir);
    }
}

unittest {
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp, strlen;
    import gretl.gerror : g_error_free;
    import std.file : remove, write;

    write("/tmp/gfileutils_test.txt", "line1\nline2\n");
    gchar * buf;
    gsize n;
    GError * err;

    assert(g_file_get_contents("/tmp/gfileutils_test.txt", &buf, &n, &err));
    assert(n == 12 && strcmp(buf, "line1\nline2\n") == 0);
    free(buf);
    write("/tmp/gfileutils_test.txt", "");
    assert(g_file_get_contents("/tmp/gfileutils_test.txt", &buf, &n, &err) && n == 0 && strlen(buf) == 0);
    free(buf);
    remove("/tmp/gfileutils_test.txt");
    assert(!g_file_get_contents("/tmp/gfileutils_no_such_file", &buf, &n, &err));
    assert(err !is null && strlen(err.message) > 0 && buf is null);
    g_error_free(err);

    gchar * cwd = g_get_current_dir();

    assert(cwd !is null && cwd[0] == '/');
    free(cwd);
}

unittest {
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp;
    import std.conv : octal;
    import std.file : exists, isDir, rmdirRecurse;

    /* path helpers */
    assert(g_path_is_absolute("/a/b") && !g_path_is_absolute("a/b"));
    gchar * d = g_path_get_dirname("/usr/lib/x.so");
    assert(strcmp(d, "/usr/lib") == 0);
    free(d);
    d = g_path_get_dirname("x.so");
    assert(strcmp(d, ".") == 0);
    free(d);
    d = g_path_get_dirname("/x");
    assert(strcmp(d, "/") == 0);
    free(d);

    gchar * f = g_build_filename("/a/".ptr, "/b/".ptr, "c.txt".ptr, null);
    assert(strcmp(f, "/a/b/c.txt") == 0);
    free(f);
    f = g_build_filename("a".ptr, "".ptr, "b".ptr, null);
    assert(strcmp(f, "a/b") == 0);
    free(f);

    /* directories */
    string top = "/tmp/gfileutils_dirtest";

    if (exists(top)) {
        rmdirRecurse(top);
    }
    assert(g_mkdir_with_parents("/tmp/gfileutils_dirtest/a/b/c", octal!755) == 0);
    assert(g_mkdir_with_parents("/tmp/gfileutils_dirtest/a/b/c", octal!755) == 0);
    assert(g_file_test("/tmp/gfileutils_dirtest/a/b", GFileTest.is_dir));
    assert(!g_file_test("/tmp/gfileutils_dirtest/a/b", GFileTest.is_regular));

    char[64] tmpl = 0;
    tmpl[0 .. 36] = "/tmp/gfileutils_dirtest/t.XXXXXX\0\0\0\0\0"[0 .. 36];
    int fd = g_mkstemp(tmpl.ptr);
    assert(fd >= 0);
    version (Posix) import core.sys.posix.unistd : close; else import gretl.winport : close;
    close(fd);
    assert(g_file_test(tmpl.ptr, GFileTest.is_regular));

    GDir * dir = g_dir_open("/tmp/gfileutils_dirtest", 0, null);
    assert(dir !is null);
    int n = 0;
    while (g_dir_read_name(dir) !is null) {
        n++;
    }
    assert(n == 2); /* "a" and the temp file; "." and ".." are skipped */
    g_dir_rewind(dir);
    assert(g_dir_read_name(dir) !is null);
    g_dir_close(dir);
    assert(g_dir_open("/tmp/gfileutils_no_such_dir", 0, null) is null);

    gchar * sh = g_find_program_in_path("sh");
    assert(sh !is null);
    free(sh);
    assert(g_find_program_in_path("no_such_program_xyz") is null);

    rmdirRecurse(top);
}

/* g_path_get_basename(): the last element of @file_name, ignoring trailing slashes; a newly allocated string (free with g_free) */
gchar * g_path_get_basename(const(gchar) * file_name)
{
    import core.stdc.string : strlen;

    size_t n = strlen(file_name);

    if (n == 0) {
        return g_strdup(".");
    }
    while (n > 1 && file_name[n - 1] == '/') {
        n--;
    }
    if (n == 1 && file_name[0] == '/') {
        return g_strdup("/");
    }
    size_t i = n;

    while (i > 0 && file_name[i - 1] != '/') {
        i--;
    }
    return g_strndup(file_name + i, n - i);
}
