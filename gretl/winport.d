/* Windows replacements for the POSIX functions that the rest of the port uses. On Posix this module is empty
   (the files import core.sys.posix directly); on Windows it supplies the same names, built on the C runtime
   and Win32. Nothing here is in the original C library: see log.md, "Windows port". */

module gretl.winport;

version (Windows):

import core.stdc.stdio : FILE;
import core.stdc.time : time_t, tm;
import core.sys.windows.windows;

alias mode_t = int;
alias uid_t = int;
alias gid_t = int;
alias pid_t = int;
alias ssize_t = ptrdiff_t;

enum S_IFMT = 0xF000;
enum S_IFDIR = 0x4000;
enum S_IFREG = 0x8000;
enum S_IFLNK = 0xA000;      // never set by the Windows C runtime
enum S_IXUSR = 0x40;
enum S_IXGRP = 0x8;
enum S_IXOTH = 0x1;
enum O_RDONLY = 0;
enum O_WRONLY = 1;
enum O_CREAT = 0x100;
enum O_BINARY = 0x8000;
enum X_OK = 0;              // the C runtime has no execute test: existence is the nearest thing
enum W_OK = 2;
enum RTLD_NOW = 2;
enum RTLD_LAZY = 1;

bool S_ISDIR(int m) { return (m & S_IFMT) == S_IFDIR; }
bool S_ISREG(int m) { return (m & S_IFMT) == S_IFREG; }
bool S_ISLNK(int m) { return false; }

/* the C runtime's struct _stat64 */
struct stat_t {
    uint st_dev;
    ushort st_ino;
    ushort st_mode;
    short st_nlink;
    short st_uid;
    short st_gid;
    uint st_rdev;
    long st_size;
    long st_atime;
    long st_mtime;
    long st_ctime;
}

extern (C) nothrow @nogc {
    private int _stat64(const(char) *, stat_t *);
    private int _fstat64(int, stat_t *);
    private int _access(const(char) *, int);
    private int _mkdir(const(char) *);
    private int _rmdir(const(char) *);
    private int _chdir(const(char) *);
    private char * _getcwd(char *, int);
    private int _getpid();
    private FILE * _fdopen(int, const(char) *);
    private int _putenv_s(const(char) *, const(char) *);
    private int _open(const(char) *, int, ...);
    private int _close(int);
    private int _read(int, void *, uint);
    private int _write(int, const(void) *, uint);
    private int _fseeki64(FILE *, long, int);
    private long _ftelli64(FILE *);
    private int _mktemp_s(char *, size_t);
    pragma(mangle, "_localtime64_s") private int localtime_s(tm *, const(long) *);
    pragma(mangle, "_gmtime64_s") private int gmtime_s(tm *, const(long) *);
    private long _mkgmtime(tm *);
}

int stat(const(char) * path, stat_t * buf) { return _stat64(path, buf); }
int lstat(const(char) * path, stat_t * buf) { return _stat64(path, buf); }
int fstat(int fd, stat_t * buf) { return _fstat64(fd, buf); }
int access(const(char) * path, int mode) { return _access(path, mode); }
int mkdir(const(char) * path, int mode = 0) { return _mkdir(path); }
int rmdir(const(char) * path) { return _rmdir(path); }
int chdir(const(char) * path) { return _chdir(path); }
char * getcwd(char * buf, size_t n) { return _getcwd(buf, cast(int) n); }
int getpid() { return _getpid(); }
FILE * fdopen(int fd, const(char) * mode) { return _fdopen(fd, mode); }
int setenv(const(char) * name, const(char) * value, int overwrite) { return _putenv_s(name, value); }
int open(const(char) * path, int flags, int mode = 0) { return _open(path, flags | O_BINARY, mode); }
int close(int fd) { return _close(fd); }
ssize_t read(int fd, void * buf, size_t n) { return _read(fd, buf, cast(uint) n); }
ssize_t write(int fd, const(void) * buf, size_t n) { return _write(fd, buf, cast(uint) n); }
int fseeko(FILE * fp, long off, int whence) { return _fseeki64(fp, off, whence); }
long ftello(FILE * fp) { return _ftelli64(fp); }

/* There are no users or groups: every file counts as owned by the caller. */
uid_t getuid() { return 0; }
gid_t getgid() { return 0; }

int gethostname(char * name, size_t n)
{
    DWORD len = cast(DWORD) n;

    return GetComputerNameA(name, &len) ? 0 : -1;
}

/* mkstemp: @tmpl ends in XXXXXX; the file is created in binary read-write mode. */
int mkstemp(char * tmpl)
{
    import core.stdc.string : strlen;

    if (_mktemp_s(tmpl, strlen(tmpl) + 1) != 0) {
        return -1;
    }
    return _open(tmpl, O_CREAT | 2 | O_BINARY | 0x0200 /* _O_EXCL */, 0x180 /* _S_IREAD | _S_IWRITE */);
}

/* mkdtemp: the same, creating a directory */
char * mkdtemp(char * tmpl)
{
    import core.stdc.string : strlen;

    if (_mktemp_s(tmpl, strlen(tmpl) + 1) != 0 || _mkdir(tmpl) != 0) {
        return null;
    }
    return tmpl;
}

/* directory reading, with the interface of opendir() and friends over FindFirstFile */

struct dirent {
    char[MAX_PATH] d_name;
}

struct DIR {
    HANDLE h;
    WIN32_FIND_DATAA data;
    bool first;
    char[MAX_PATH + 3] pattern;
    dirent ent;
}

DIR * opendir(const(char) * path)
{
    import core.stdc.stdlib : calloc, free;
    import core.stdc.string : strlen, strcpy;
    import core.stdc.errno : errno, ENOENT;

    DIR * d = cast(DIR *) calloc(1, DIR.sizeof);
    size_t n = strlen(path);

    if (d == null || n > MAX_PATH - 1) {
        free(d);
        return null;
    }
    strcpy(d.pattern.ptr, path);
    if (n > 0 && path[n-1] != '/' && path[n-1] != '\\') {
        d.pattern[n++] = '\\';
    }
    d.pattern[n] = '*';
    d.pattern[n+1] = 0;
    d.h = FindFirstFileA(d.pattern.ptr, &d.data);
    if (d.h == INVALID_HANDLE_VALUE) {
        free(d);
        errno = ENOENT;
        return null;
    }
    d.first = true;
    return d;
}

dirent * readdir(DIR * d)
{
    import core.stdc.string : strcpy;

    if (d.first) {
        d.first = false;
    } else if (!FindNextFileA(d.h, &d.data)) {
        return null;
    }
    strcpy(d.ent.d_name.ptr, d.data.cFileName.ptr);
    return &d.ent;
}

void rewinddir(DIR * d)
{
    FindClose(d.h);
    d.h = FindFirstFileA(d.pattern.ptr, &d.data);
    d.first = true;
}

int closedir(DIR * d)
{
    import core.stdc.stdlib : free;

    FindClose(d.h);
    free(d);
    return 0;
}

/* dynamic loading */

void * dlopen(const(char) * path, int flags)
{
    return path == null ? cast(void *) GetModuleHandleA(null)
        : cast(void *) LoadLibraryExA(path, null, LOAD_WITH_ALTERED_SEARCH_PATH);
}

void * dlsym(void * handle, const(char) * name) { return cast(void *) GetProcAddress(cast(HMODULE) handle, name); }
int dlclose(void * handle) { return FreeLibrary(cast(HMODULE) handle) ? 0 : -1; }

const(char) * dlerror()
{
    return "could not load the library or find the symbol";
}

/* The stack limit: the default of 1 MB for a Windows thread. */

enum RLIMIT_STACK = 3;

struct rlimit {
    ulong rlim_cur;
    ulong rlim_max;
}

int getrlimit(int resource, rlimit * rl)
{
    rl.rlim_cur = rl.rlim_max = 1024 * 1024;
    return 0;
}

/* time */

/* time_t is 32 bits in the D runtime on Windows but the C runtime's time functions use 64-bit values: callers
   keep their times in a long */
tm * localtime_r(const(long) * t, tm * out_) { return localtime_s(out_, t) == 0 ? out_ : null; }
tm * gmtime_r(const(long) * t, tm * out_) { return gmtime_s(out_, t) == 0 ? out_ : null; }
long timegm(tm * t) { return _mkgmtime(t); }

/* strptime for the C locale: the conversions %a %A %b %B %h %c-less %d %e %H %I %j %m %M %n %t %p %S %y %Y %T %D %F %R
   %% and white space. Returns a pointer to the first unparsed character, or null. */

private immutable string[12] monthNames = ["january", "february", "march", "april", "may", "june", "july",
                                           "august", "september", "october", "november", "december"];
private immutable string[7] dayNames = ["sunday", "monday", "tuesday", "wednesday", "thursday", "friday", "saturday"];

private bool matchName(ref const(char) * s, const(string[]) names, out int idx)
{
    import core.stdc.ctype : tolower;

    foreach (i, nm; names) {
        foreach (len; [nm.length, 3]) {
            size_t k = 0;

            while (k < len && s[k] != 0 && tolower(s[k]) == nm[k]) {
                k++;
            }
            if (k == len) {
                idx = cast(int) i;
                s += len;
                return true;
            }
        }
    }
    return false;
}

private bool getNum(ref const(char) * s, int minv, int maxv, int maxdigits, out int val)
{
    int n = 0, v = 0;

    while (*s == ' ') {
        s++;
    }
    while (n < maxdigits && *s >= '0' && *s <= '9') {
        v = v * 10 + (*s - '0');
        s++;
        n++;
    }
    val = v;
    return n > 0 && v >= minv && v <= maxv;
}

char * strptime(const(char) * s, const(char) * fmt, tm * t)
{
    import core.stdc.ctype : isspace;

    int v;
    bool pm = false, have_pm = false;

    while (*fmt) {
        char c = *fmt++;

        if (isspace(c)) {
            while (isspace(*s)) {
                s++;
            }
            continue;
        }
        if (c != '%') {
            if (*s != c) {
                return null;
            }
            s++;
            continue;
        }
        c = *fmt++;
        if (c == 'E' || c == 'O') {
            c = *fmt++;
        }
        switch (c) {
        case '%':
            if (*s++ != '%') {
                return null;
            }
            break;
        case 'n': case 't':
            while (isspace(*s)) {
                s++;
            }
            break;
        case 'Y':
            if (!getNum(s, 0, 9999, 4, v)) {
                return null;
            }
            t.tm_year = v - 1900;
            break;
        case 'y':
            if (!getNum(s, 0, 99, 2, v)) {
                return null;
            }
            t.tm_year = v < 69 ? v + 100 : v;
            break;
        case 'C':
            if (!getNum(s, 0, 99, 2, v)) {
                return null;
            }
            break;
        case 'm':
            if (!getNum(s, 1, 12, 2, v)) {
                return null;
            }
            t.tm_mon = v - 1;
            break;
        case 'd': case 'e':
            if (!getNum(s, 1, 31, 2, v)) {
                return null;
            }
            t.tm_mday = v;
            break;
        case 'H': case 'k':
            if (!getNum(s, 0, 23, 2, v)) {
                return null;
            }
            t.tm_hour = v;
            break;
        case 'I': case 'l':
            if (!getNum(s, 1, 12, 2, v)) {
                return null;
            }
            t.tm_hour = v % 12;
            break;
        case 'M':
            if (!getNum(s, 0, 59, 2, v)) {
                return null;
            }
            t.tm_min = v;
            break;
        case 'S':
            if (!getNum(s, 0, 61, 2, v)) {
                return null;
            }
            t.tm_sec = v;
            break;
        case 'j':
            if (!getNum(s, 1, 366, 3, v)) {
                return null;
            }
            t.tm_yday = v - 1;
            break;
        case 'p': case 'P':
            if ((s[0] == 'P' || s[0] == 'p') && (s[1] == 'M' || s[1] == 'm')) {
                pm = true;
            } else if (!((s[0] == 'A' || s[0] == 'a') && (s[1] == 'M' || s[1] == 'm'))) {
                return null;
            }
            have_pm = true;
            s += 2;
            break;
        case 'b': case 'B': case 'h':
            if (!matchName(s, monthNames[], v)) {
                return null;
            }
            t.tm_mon = v;
            break;
        case 'a': case 'A':
            if (!matchName(s, dayNames[], v)) {
                return null;
            }
            t.tm_wday = v;
            break;
        case 'T': case 'D': case 'F': case 'R': {
            string sub = c == 'T' ? "%H:%M:%S" : c == 'D' ? "%m/%d/%y" : c == 'F' ? "%Y-%m-%d" : "%H:%M";
            char[16] buf = 0;

            buf[0 .. sub.length] = sub[];
            s = strptime(s, buf.ptr, t);
            if (s == null) {
                return null;
            }
            break;
        }
        default:
            return null;
        }
    }
    if (have_pm && pm) {
        t.tm_hour += 12;
    }
    return cast(char *) s;
}

ptrdiff_t readlink(const(char) * path, char * buf, size_t n) { return -1; }

/* Child processes. gspawn.d registers the process handle of a child it starts, under its process id; waitpid()
   waits on it. The status has the layout of a POSIX wait status: the exit code in bits 8-15. */

private __gshared HANDLE[int] children;

void registerChild(int pid, HANDLE h)
{
    synchronized {
        children[pid] = h;
    }
}

int waitpid(int pid, int * status, int options)
{
    HANDLE h;

    synchronized {
        auto p = pid in children;

        if (p is null) {
            return -1;
        }
        h = *p;
        children.remove(pid);
    }
    WaitForSingleObject(h, INFINITE);

    DWORD code = 0;

    GetExitCodeProcess(h, &code);
    CloseHandle(h);
    if (status != null) {
        *status = cast(int) ((code & 0xff) << 8);
    }
    return pid;
}

bool WIFEXITED(int status) { return true; }
int WEXITSTATUS(int status) { return (status >> 8) & 0xff; }

/* read-only file mapping, enough for mmap(null, len, PROT_READ, MAP_SHARED, fd, 0) */

enum PROT_READ = 1;
enum MAP_SHARED = 1;
enum MAP_FAILED = cast(void *) -1;

private extern (C) ptrdiff_t _get_osfhandle(int);

void * mmap(void * addr, size_t length, int prot, int flags, int fd, long offset)
{
    HANDLE fh = cast(HANDLE) _get_osfhandle(fd);
    HANDLE mh = CreateFileMappingA(fh, null, PAGE_READONLY, 0, 0, null);

    if (mh == null) {
        return MAP_FAILED;
    }

    void * p = MapViewOfFile(mh, FILE_MAP_READ, cast(DWORD) (offset >> 32), cast(DWORD) offset, length);

    CloseHandle(mh);
    return p == null ? MAP_FAILED : p;
}

int munmap(void * addr, size_t length) { return UnmapViewOfFile(addr) ? 0 : -1; }

/* nl_langinfo(CODESET): the ANSI code page as an iconv name, such as "CP1252" */

char * nl_langinfo(int item)
{
    import core.stdc.stdio : snprintf;

    __gshared char[16] cp;

    snprintf(cp.ptr, cp.length, "CP%u", GetACP());
    return cp.ptr;
}

private extern (C) nothrow @nogc {
    int _dup(int);
    int _dup2(int, int);
}

int dup(int fd) { return _dup(fd); }
int dup2(int fd, int fd2) { return _dup2(fd, fd2); }

/* The C runtime opens files in text mode unless told otherwise, translating "\n" to "\r\n" on output; the port (like
   the C code on Unix) reads and writes the bytes it asks for, so files are binary by default. */

private extern (C) int _set_fmode(int);

shared static this()
{
    _set_fmode(0x8000);     // _O_BINARY
}
