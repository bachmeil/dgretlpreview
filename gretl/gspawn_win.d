/* g_spawn_async_with_pipes() for Windows, over CreateProcess and anonymous pipes. The three pipe ends handed to the
   caller are C runtime file descriptors, so the read(), write() and close() of the callers work on them. The
   limitations of gspawn.d apply (no child_setup, user_data or envp). */

module gretl.gspawn_win;

version (Windows):

import core.sys.windows.windows;
import core.stdc.string : strlen;

import gretl.gerror : GError, g_set_error_literal_msg;
import gretl.glib : gboolean, gchar, gint;
import gretl.gspawn : GSpawnFlags;
import gretl.winport : registerChild;

private extern (C) int _open_osfhandle(ptrdiff_t, int);

private string quote_arg(const(char) * a)
{
    import std.string : fromStringz;

    string s = fromStringz(a).idup;
    bool plain = s.length > 0;

    foreach (c; s) {
        if (c == ' ' || c == '\t' || c == '"') {
            plain = false;
        }
    }
    if (plain) {
        return s;
    }

    string r = "\"";
    size_t bs = 0;

    foreach (c; s) {
        if (c == '\\') {
            bs++;
        } else if (c == '"') {
            foreach (i; 0 .. 2 * bs + 1) r ~= '\\';
            r ~= '"';
            bs = 0;
        } else {
            foreach (i; 0 .. bs) r ~= '\\';
            r ~= c;
            bs = 0;
        }
    }
    foreach (i; 0 .. 2 * bs) r ~= '\\';
    return r ~ "\"";
}

/* create a pipe whose end @child_end is inheritable and whose other end is not */

private bool make_pipe(HANDLE * parent_end, HANDLE * child_end, bool child_reads)
{
    SECURITY_ATTRIBUTES sa;
    HANDLE r, w;

    sa.nLength = SECURITY_ATTRIBUTES.sizeof;
    sa.bInheritHandle = TRUE;
    if (!CreatePipe(&r, &w, &sa, 0)) {
        return false;
    }
    *child_end = child_reads ? r : w;
    *parent_end = child_reads ? w : r;
    SetHandleInformation(*parent_end, HANDLE_FLAG_INHERIT, 0);
    return true;
}

private HANDLE open_nul()
{
    SECURITY_ATTRIBUTES sa;

    sa.nLength = SECURITY_ATTRIBUTES.sizeof;
    sa.bInheritHandle = TRUE;
    return CreateFileA("NUL", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa,
                       OPEN_EXISTING, 0, null);
}

gboolean g_spawn_async_with_pipes(const(gchar) * working_directory, gchar ** argv, gchar ** envp,
                                  GSpawnFlags flags, void * child_setup, void * user_data,
                                  int * child_pid, gint * standard_input, gint * standard_output,
                                  gint * standard_error, GError ** error)
{
    assert(child_setup == null && user_data == null && envp == null,
           "gspawn.d: child_setup, user_data and envp are not supported");

    HANDLE[3] parent = INVALID_HANDLE_VALUE, child = INVALID_HANDLE_VALUE;
    bool[3] own = false;    /* child[i] is ours to close */
    gint *[3] want = [standard_input, standard_output, standard_error];

    foreach (i; 0 .. 3) {
        if (want[i] != null) {
            if (!make_pipe(&parent[i], &child[i], i == 0)) {
                g_set_error_literal_msg(error, "Failed to create pipes");
                return false;
            }
            own[i] = true;
        } else if ((i == 1 && (flags & GSpawnFlags.stdout_to_dev_null)) ||
                   (i == 2 && (flags & GSpawnFlags.stderr_to_dev_null))) {
            child[i] = open_nul();
            own[i] = true;
        } else {
            child[i] = GetStdHandle(i == 0 ? STD_INPUT_HANDLE : i == 1 ? STD_OUTPUT_HANDLE : STD_ERROR_HANDLE);
        }
    }

    string cmd;

    for (size_t i = 0; argv[i] != null; i++) {
        cmd ~= (i ? " " : "") ~ quote_arg(argv[i]);
    }

    STARTUPINFOA si;
    PROCESS_INFORMATION pi;

    si.cb = STARTUPINFOA.sizeof;
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = child[0];
    si.hStdOutput = child[1];
    si.hStdError = child[2];

    char[] cl = (cmd ~ '\0').dup;
    BOOL ok = CreateProcessA(null, cl.ptr, null, null, TRUE, CREATE_NO_WINDOW, null,
                             working_directory, &si, &pi);

    foreach (i; 0 .. 3) {
        if (own[i]) {
            CloseHandle(child[i]);
        }
    }
    if (!ok) {
        foreach (i; 0 .. 3) {
            if (parent[i] != INVALID_HANDLE_VALUE) {
                CloseHandle(parent[i]);
            }
        }
        g_set_error_literal_msg(error, "Failed to execute child process");
        return false;
    }

    CloseHandle(pi.hThread);
    registerChild(cast(int) pi.dwProcessId, pi.hProcess);
    if (child_pid != null) {
        *child_pid = cast(int) pi.dwProcessId;
    }
    foreach (i; 0 .. 3) {
        if (want[i] != null) {
            *want[i] = _open_osfhandle(cast(ptrdiff_t) parent[i], i == 0 ? 0 : 0x8000 /* _O_RDONLY | _O_BINARY */);
        }
    }
    return true;
}
