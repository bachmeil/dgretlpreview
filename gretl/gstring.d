/*
   Replacement for the glib GString (a growable string), written with malloc'd storage.

   GString is a struct with the string in `str` (NUL-terminated), its length in `len` and the allocated size in
   `allocated_len`, as in glib. Provided: g_string_new, g_string_sized_new, g_string_append, g_string_append_c,
   g_string_append_printf (an extern (C) variadic, as in the C), g_string_assign, g_string_truncate and
   g_string_free (which returns the character data if free_segment is false, else null). Add others here as ported
   files need them.

   Used by: gretl_bundle.d.
*/

module gretl.gstring;


import core.stdc.stdarg;
import core.stdc.stdio : vsnprintf;
import core.stdc.stdlib : free, malloc, realloc;
import core.stdc.string : memcpy, strlen;

import gretl.glib : gchar, gboolean, gsize, gssize;

struct GString {
    gchar * str;
    gsize len;
    gsize allocated_len;
}

private void ensure(GString * s, gsize extra)
{
    gsize need = s.len + extra + 1;

    if (need > s.allocated_len) {
        gsize n = s.allocated_len < 16 ? 16 : s.allocated_len;

        while (n < need) {
            n *= 2;
        }
        s.str = cast(gchar *) realloc(s.str, n);
        s.allocated_len = n;
    }
}

GString * g_string_sized_new(gsize dfl_size)
{
    GString * s = cast(GString *) malloc(GString.sizeof);

    s.str = null;
    s.len = 0;
    s.allocated_len = 0;
    ensure(s, dfl_size);
    s.str[0] = 0;
    return s;
}

GString * g_string_new(const(gchar) * init)
{
    GString * s = g_string_sized_new(init != null ? strlen(init) : 2);

    if (init != null) {
        g_string_append(s, init);
    }
    return s;
}

GString * g_string_append(GString * s, const(gchar) * val)
{
    gsize n = strlen(val);

    ensure(s, n);
    memcpy(s.str + s.len, val, n + 1);
    s.len += n;
    return s;
}

GString * g_string_append_c(GString * s, gchar c)
{
    ensure(s, 1);
    s.str[s.len++] = c;
    s.str[s.len] = 0;
    return s;
}

GString * g_string_assign(GString * s, const(gchar) * rval)
{
    s.len = 0;
    s.str[0] = 0;
    return g_string_append(s, rval);
}

GString * g_string_truncate(GString * s, gsize len)
{
    if (len < s.len) {
        s.len = len;
        s.str[len] = 0;
    }
    return s;
}

extern (C) void g_string_append_printf(GString * s, const(gchar) * format, ...)
{
    va_list args, args2;

    va_start(args, format);
    va_copy(args2, args);
    int n = vsnprintf(null, 0, format, args2);

    va_end(args2);
    if (n > 0) {
        ensure(s, n);
        vsnprintf(s.str + s.len, n + 1, format, args);
        s.len += n;
    }
    va_end(args);
}

/* returns the string (to be freed by the caller) if free_segment is false, else frees it and returns null */

gchar * g_string_free(GString * s, gboolean free_segment)
{
    gchar * ret = free_segment ? null : s.str;

    if (free_segment) {
        free(s.str);
    }
    free(s);
    return ret;
}

unittest
{
    GString * s = g_string_new("ab");

    g_string_append_c(s, 'c');
    g_string_append_printf(s, "%d-%s", 12, "x".ptr);
    assert(s.len == 7);
    gchar * r = g_string_free(s, false);

    assert(r[0 .. 7] == "abc12-x");
    assert(r[7] == 0);
    free(r);
}
