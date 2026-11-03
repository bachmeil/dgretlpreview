/*
   Replacements for the glib process-spawning functions used by the C code:
   g_spawn_sync(), g_spawn_command_line_sync(), g_shell_parse_argv(),
   together with the GError/GSpawnFlags bits that go with them, written on
   top of D's std.process so that the port does not depend on glib.

   Semantics (following glib's documentation):

   - g_spawn_sync() runs argv[0] with the given arguments, waits for it to
     finish, and returns its standard output and standard error as
     malloc'd, NUL-terminated strings (free with g_free()/free()), and its
     "wait status" in *exit_status (exit code * 256 for a normal exit, or
     the signal number for a process killed by a signal; 0 means success).
     It returns false, with a message in *error, if the program could not
     be started (e.g. not found). argv is a NULL-terminated C array.
   - working_directory may be null (inherit); envp may be null (inherit
     the environment) or a NULL-terminated array of "NAME=value" strings,
     in which case it REPLACES the environment (as in glib).
   - The search of PATH is always done (as with G_SPAWN_SEARCH_PATH, the only
     flag the callers use); the flags argument is accepted but ignored.
     child_setup and user_data are not supported and must be null.
   - The child's standard input is /dev/null. The two outputs are captured
     through temporary files (in the system temporary directory), so no pipe can fill up and block the child.
   - g_spawn_command_line_sync() splits a command line with
     g_shell_parse_argv() (no shell is involved).
   - g_shell_parse_argv() implements the usual quoting rules: whitespace
     separates words; single quotes preserve everything; inside double
     quotes a backslash only escapes ", \, $ and `, and is otherwise
     literal; outside quotes a backslash escapes the next character.
     Unmatched quotes, an empty command line or a trailing backslash are
     errors. No expansion (variables, globs, ~) is done.

   - g_spawn_async_with_pipes() starts the program with fork/exec (always searching PATH) and, for each of
     the stdin/stdout/stderr descriptor pointers that is not null, connects a pipe and returns the parent's end;
     a null stdout/stderr pointer gives /dev/null with the matching G_SPAWN_*_TO_DEV_NULL flag and otherwise
     is inherited. The child is never reaped here (the caller uses waitpid, as with G_SPAWN_DO_NOT_REAP_CHILD).
     A failed exec is reported through a close-on-exec pipe, as glib does. envp, child_setup and the other
     flags are not supported.

   Differences from glib: a GError here is a plain struct {code, message}
   allocated with malloc (free with g_error_free()); the code is always 1.
   The GC is used by std.process.

   Used by: gretl_utils.d, graphing.d.
*/

module gretl.gspawn;


import core.stdc.stdlib : free, malloc;
import core.stdc.string : memcpy, strlen;

import std.process : ProcessException, Config, spawnProcess, wait;
import std.stdio : File;
import std.string : fromStringz;

import gretl.gerror : GError, g_error_free, g_set_error_literal_msg;
import gretl.glib : gboolean, gchar, gint;

public import gretl.gerror : GError, g_error_free;
version (Windows) public import gretl.gspawn_win : g_spawn_async_with_pipes;

enum GSpawnFlags {
    default_ = 0,
    leave_descriptors_open = 1 << 0,
    do_not_reap_child = 1 << 1,
    search_path = 1 << 2,
    stdout_to_dev_null = 1 << 3,
    stderr_to_dev_null = 1 << 4,
    child_inherits_stdin = 1 << 5
}

alias G_SPAWN_SEARCH_PATH = GSpawnFlags.search_path;
alias G_SPAWN_DO_NOT_REAP_CHILD = GSpawnFlags.do_not_reap_child;
alias G_SPAWN_STDOUT_TO_DEV_NULL = GSpawnFlags.stdout_to_dev_null;
alias G_SPAWN_STDERR_TO_DEV_NULL = GSpawnFlags.stderr_to_dev_null;

private alias set_error = g_set_error_literal_msg;

version (Windows) private enum devnull = "NUL"; else private enum devnull = "/dev/null";

/* returns false (and sets @msg) on a parse error */

bool shell_parse(string s, ref string[] argv, ref string msg)
{
    char[] cur;
    bool have = false;
    size_t i = 0;

    argv = null;

    while (i < s.length) {
        char c = s[i];

        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            if (have) {
                argv ~= cur.idup;
                cur = null;
                have = false;
            }
            i++;
        } else if (c == '\'') {
            size_t j = i + 1;

            while (j < s.length && s[j] != '\'') {
                cur ~= s[j++];
            }
            if (j >= s.length) {
                msg = "Text ended before matching quote was found for '";
                return false;
            }
            have = true;
            i = j + 1;
        } else if (c == '"') {
            size_t j = i + 1;

            while (j < s.length && s[j] != '"') {
                if (s[j] == '\\' && j + 1 < s.length &&
                    (s[j+1] == '"' || s[j+1] == '\\' || s[j+1] == '$' || s[j+1] == '`')) {
                    cur ~= s[j+1];
                    j += 2;
                } else if (s[j] == '\\' && j + 1 < s.length && s[j+1] == '\n') {
                    j += 2;
                } else {
                    cur ~= s[j++];
                }
            }
            if (j >= s.length) {
                msg = "Text ended before matching quote was found for \"";
                return false;
            }
            have = true;
            i = j + 1;
        } else if (c == '\\') {
            if (i + 1 >= s.length) {
                msg = "Text ended just after a '\\' character";
                return false;
            }
            if (s[i+1] != '\n') {
                cur ~= s[i+1];
                have = true;
            }
            i += 2;
        } else {
            cur ~= c;
            have = true;
            i++;
        }
    }
    if (have) {
        argv ~= cur.idup;
    }
    if (argv.length == 0) {
        msg = "Text was empty (or contained only whitespace)";
        return false;
    }
    return true;
}

/* a unique temporary file name (the file is created empty and closed) */

private string temp_path()
{
    version (Posix) {
        import core.sys.posix.stdlib : mkstemp;
        import core.sys.posix.unistd : close;
    } else {
        import gretl.winport : mkstemp, close;
    }
    import std.file : tempDir;

    char[] tmpl = (tempDir() ~ "/gspawnXXXXXX\0").dup;
    int fd = mkstemp(tmpl.ptr);

    if (fd < 0) {
        throw new Exception("could not create a temporary file");
    }
    close(fd);
    return tmpl[0 .. $ - 1].idup;
}

private gchar * dup_file_contents(string path)
{
    import std.file : read;

    ubyte[] data = cast(ubyte[]) read(path);

    gchar * ret = cast(gchar *) malloc(data.length + 1);

    if (data.length > 0) {
        memcpy(ret, data.ptr, data.length);
    }
    ret[data.length] = 0;
    return ret;
}

private gboolean spawn_sync_d(string workdir, string[] args, string[string] env, bool replace_env,
                              gchar ** standard_output, gchar ** standard_error,
                              gint * exit_status, GError ** error)
{
    string outp, errp;

    scope(exit) {
        import std.file : exists, remove;

        if (outp.length && exists(outp)) remove(outp);
        if (errp.length && exists(errp)) remove(errp);
    }

    try {
        outp = temp_path();
        errp = temp_path();
        Config cfg = replace_env ? Config.newEnv : Config.none;
        auto pid = spawnProcess(args, File(devnull, "r"), File(outp, "w"), File(errp, "w"), env, cfg,
                                workdir.length ? workdir : null);
        int code = wait(pid);

        if (standard_output != null) {
            *standard_output = dup_file_contents(outp);
        }
        if (standard_error != null) {
            *standard_error = dup_file_contents(errp);
        }
        if (exit_status != null) {
            /* wait status as returned by waitpid() */
            *exit_status = code >= 0 ? code << 8 : -code;
        }
        return true;
    } catch (ProcessException e) {
        set_error(error, e.msg);
    } catch (Exception e) {
        set_error(error, e.msg);
    }
    if (standard_output != null) {
        *standard_output = null;
    }
    if (standard_error != null) {
        *standard_error = null;
    }
    return false;
}

gboolean g_spawn_sync(const(gchar) * working_directory, gchar ** argv, gchar ** envp,
                      GSpawnFlags flags, void * child_setup, void * user_data,
                      gchar ** standard_output, gchar ** standard_error,
                      gint * exit_status, GError ** error)
{
    assert(child_setup == null && user_data == null, "gspawn.d: child_setup is not supported");

    string[] args;
    string[string] env;

    for (size_t i = 0; argv[i] != null; i++) {
        args ~= fromStringz(argv[i]).idup;
    }
    if (envp != null) {
        for (size_t i = 0; envp[i] != null; i++) {
            string kv = fromStringz(envp[i]).idup;
            size_t eq = 0;

            while (eq < kv.length && kv[eq] != '=') {
                eq++;
            }
            if (eq < kv.length) {
                env[kv[0 .. eq]] = kv[eq + 1 .. $];
            }
        }
    }
    string wd = working_directory != null ? fromStringz(working_directory).idup : "";

    return spawn_sync_d(wd, args, env, envp != null, standard_output, standard_error,
                        exit_status, error);
}

gboolean g_spawn_command_line_sync(const(gchar) * command_line,
                                   gchar ** standard_output, gchar ** standard_error,
                                   gint * exit_status, GError ** error)
{
    string[] args;
    string msg;

    if (!shell_parse(fromStringz(command_line).idup, args, msg)) {
        set_error(error, msg);
        return false;
    }

    return spawn_sync_d("", args, null, false, standard_output, standard_error,
                        exit_status, error);
}

unittest {
    import core.stdc.string : strcmp;

    string[] a;
    string m;

    assert(shell_parse("a b  c", a, m) && a == ["a", "b", "c"]);
    assert(shell_parse(`say "hello world" 'x y' z\ w`, a, m) && a == ["say", "hello world", "x y", "z w"]);
    assert(shell_parse(`"a\"b" 'c\d'`, a, m) && a == [`a"b`, `c\d`]);
    assert(shell_parse(`"" x`, a, m) && a == ["", "x"]);
    assert(!shell_parse("", a, m) && !shell_parse("   ", a, m));
    assert(!shell_parse(`"abc`, a, m) && !shell_parse(`it's`, a, m) && !shell_parse(`abc\`, a, m));

    gchar * so, se;
    gint st;
    GError * err;

    { int ok = g_spawn_command_line_sync("echo 'hi there'", &so, &se, &st, &err); if (!ok) { import core.stdc.stdio : puts; puts(err.message); } assert(ok); }
    assert(strcmp(so, "hi there\n") == 0 && strcmp(se, "") == 0 && st == 0);
    free(so);
    free(se);

    assert(g_spawn_command_line_sync("sh -c 'echo out; echo err >&2; exit 3'", &so, &se, &st, &err));
    assert(strcmp(so, "out\n") == 0 && strcmp(se, "err\n") == 0 && st == 3 << 8);
    free(so);
    free(se);

    assert(!g_spawn_command_line_sync("no-such-program-xyz", &so, &se, &st, &err));
    assert(err !is null && strlen(err.message) > 0 && so is null);
    g_error_free(err);
    err = null;

    gchar *[4] av = [cast(gchar *) "sh".ptr, cast(gchar *) "-c".ptr, cast(gchar *) "echo $FOO; pwd".ptr, null];
    gchar *[2] ev = [cast(gchar *) "FOO=bar".ptr, null];

    assert(g_spawn_sync("/tmp", av.ptr, ev.ptr, G_SPAWN_SEARCH_PATH, null, null, &so, &se, &st, &err));
    assert(strcmp(so, "bar\n/tmp\n") == 0);
    free(so);
    free(se);

    /* large stderr must not block */
    assert(g_spawn_command_line_sync("sh -c 'head -c 300000 /dev/zero | tr \"\\0\" x >&2'", &so, &se, &st, &err));
    assert(strlen(se) == 300000);
    free(so);
    free(se);
}


import core.stdc.errno : errno;

version (Posix):

import core.sys.posix.fcntl : F_SETFD, FD_CLOEXEC, O_RDONLY, O_WRONLY, fcntl, open;
import core.sys.posix.sys.wait : waitpid;
import core.sys.posix.unistd : _exit, chdir, close, dup2, execvp, fork, pipe, read;

private bool make_pipe(ref int[2] fds)
{
    return pipe(fds) == 0;
}

gboolean g_spawn_async_with_pipes(const(gchar) * working_directory, gchar ** argv, gchar ** envp,
                                  GSpawnFlags flags, void * child_setup, void * user_data,
                                  int * child_pid, gint * standard_input, gint * standard_output,
                                  gint * standard_error, GError ** error)
{
    assert(child_setup == null && user_data == null && envp == null,
           "gspawn.d: child_setup, user_data and envp are not supported");

    int[2] pin = [-1, -1], pout = [-1, -1], perr = [-1, -1], pex = [-1, -1];

    if ((standard_input != null && !make_pipe(pin)) || (standard_output != null && !make_pipe(pout))
        || (standard_error != null && !make_pipe(perr)) || !make_pipe(pex)) {
        set_error(error, "Failed to create pipes");
        return false;
    }
    fcntl(pex[1], F_SETFD, FD_CLOEXEC);

    int pid = fork();

    if (pid < 0) {
        set_error(error, "Failed to fork");
        return false;
    }

    if (pid == 0) {
        /* child */
        close(pex[0]);
        if (standard_input != null) {
            close(pin[1]);
            dup2(pin[0], 0);
            close(pin[0]);
        }
        if (standard_output != null) {
            close(pout[0]);
            dup2(pout[1], 1);
            close(pout[1]);
        } else if (flags & GSpawnFlags.stdout_to_dev_null) {
            int dn = open("/dev/null", O_WRONLY);

            dup2(dn, 1);
            close(dn);
        }
        if (standard_error != null) {
            close(perr[0]);
            dup2(perr[1], 2);
            close(perr[1]);
        } else if (flags & GSpawnFlags.stderr_to_dev_null) {
            int dn = open("/dev/null", O_WRONLY);

            dup2(dn, 2);
            close(dn);
        }
        if (working_directory != null) {
            chdir(working_directory);
        }
        execvp(argv[0], argv);
        int e = errno;

        .write(pex[1], &e, e.sizeof);
        _exit(127);
    }

    /* parent */
    close(pex[1]);
    if (standard_input != null) {
        close(pin[0]);
    }
    if (standard_output != null) {
        close(pout[1]);
    }
    if (standard_error != null) {
        close(perr[1]);
    }

    int exec_errno = 0;
    ptrdiff_t n = read(pex[0], &exec_errno, exec_errno.sizeof);

    close(pex[0]);
    if (n > 0) {
        int status;

        waitpid(pid, &status, 0);
        foreach (fd; [pin[1], pout[0], perr[0]]) {
            if (fd >= 0) {
                close(fd);
            }
        }
        set_error(error, "Failed to execute child process");
        return false;
    }

    if (child_pid != null) {
        *child_pid = pid;
    }
    if (standard_input != null) {
        *standard_input = pin[1];
    }
    if (standard_output != null) {
        *standard_output = pout[0];
    }
    if (standard_error != null) {
        *standard_error = perr[0];
    }
    return true;
}

import core.sys.posix.unistd : write;
