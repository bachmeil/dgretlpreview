/*
   Replacements for the glib string functions g_strchug, g_strchomp,
   g_strstrip, g_ascii_strcasecmp, g_strdup, g_vasprintf, g_strdup_printf,
   g_strcmp0, g_strlcpy, g_strlcat, g_utf8_strlen and g_utf8_validate, g_strsplit, g_strfreev,
   g_utf8_next_char, g_utf8_prev_char, g_utf8_find_next_char, g_utf8_offset_to_pointer,
   g_utf8_substring and g_locale_to_utf8, written in plain D so that the port does not depend on glib for them.

   Behaviour follows glib's documentation:

   - Whitespace is ASCII whitespace only (g_ascii_isspace): space, \t, \n,
     \v, \f, \r. Locale is ignored.
   - g_strchug and g_strchomp modify the string in place and return it.
     g_strchug moves the remaining characters (including the terminating
     NUL) to the start of the buffer.
   - g_strstrip is a macro in glib; here it is a function.
   - g_ascii_strcasecmp compares ASCII-case-insensitively and returns the
     difference of the first pair of differing (lower-cased) bytes, 0 if
     equal.

   - g_strdup, g_vasprintf and g_strdup_printf return malloc'd memory,
     which may be released with g_free (gmem.d) or free(). g_strdup(null)
     is null. g_vasprintf returns the length, or -1 on failure.
     g_strdup_printf is an extern (C) variadic, as in the C.
   - g_utf8_strlen counts UTF-8 characters by skipping over each lead
     byte's sequence (no validation), up to max bytes, or to the NUL if
     max < 0.

   - g_strlcpy copies at most dest_size-1 bytes, always NUL-terminates (if
     dest_size > 0) and returns the length of src.
   - g_strcmp0 is strcmp that accepts null (null sorts before any string).
   - g_utf8_validate accepts only well-formed UTF-8 (no overlong forms,
     surrogates, noncharacters or values above U+10FFFF), up to max_len
     bytes or to the NUL if max_len < 0; *end, if given, receives the
     position at which validation stopped.

   - g_strsplit splits at each occurrence of a (non-empty) delimiter string, keeping empty
     tokens, up to max_tokens pieces (max_tokens < 1 means no limit); the result is a
     NULL-terminated malloc'd array of malloc'd strings, released with g_strfreev.
   - g_utf8_next_char is a macro in glib; here a function that skips by the lead byte's
     sequence length (like g_utf8_strlen). g_utf8_prev_char steps back over continuation
     bytes. g_utf8_find_next_char returns the start of the next character, or null if it would
     reach @end (or, with @end null, the position of the terminating NUL). g_utf8_offset_to_pointer
     advances by a non-negative number of characters. g_utf8_substring returns a malloc'd copy of
     characters [start_pos, end_pos), clamped to the string length.

   - g_locale_to_utf8 converts from the current locale's character set (as reported by
     nl_langinfo(CODESET), after the program has called setlocale) to UTF-8 with iconv, returning
     a malloc'd string, or null with a message in *error if the input is not valid in the locale's
     charset (so in a UTF-8 locale, invalid UTF-8 input is an error, as in glib). @len < 0 means
     NUL-terminated; *bytes_read and *bytes_written are set if non-null.

   Only these functions are provided so far; add others here (g_ascii_strncasecmp,
   ...) as ported files need them.

   Used by: libset.d, gretl_prn.d, adf_kpss.d, dataset.d, gretl_model.d, gretl_list.d, strutils.d.
*/

module gretl.gstrfuncs;


import core.stdc.stdarg : va_copy, va_end, va_list, va_start;
import core.stdc.stdio : vsnprintf;
import core.stdc.stdlib : malloc;
import core.stdc.string : memcpy, memmove, strlen, strstr;

import core.stdc.string : strcmp;

import core.stdc.stdlib : free;

import gretl.gerror : GError, g_set_error_literal_msg;
import gretl.glib : gboolean, gchar, gint, glong, gsize, gssize;

private bool g_ascii_isspace(char c)
{
    return c == ' ' || (c >= '\t' && c <= '\r');
}

private int g_ascii_tolower(int c)
{
    return (c >= 'A' && c <= 'Z') ? c + ('a' - 'A') : c;
}

gchar * g_strchug(gchar * string)
{
    gchar * start;

    for (start = string; *start != 0 && g_ascii_isspace(*start); start++) {
        ;
    }
    memmove(string, start, strlen(start) + 1);
    return string;
}

gchar * g_strchomp(gchar * string)
{
    size_t len = strlen(string);

    while (len > 0 && g_ascii_isspace(string[len - 1])) {
        string[--len] = 0;
    }
    return string;
}

gchar * g_strstrip(gchar * string)
{
    return g_strchomp(g_strchug(string));
}

gint g_ascii_strcasecmp(const(gchar) * s1, const(gchar) * s2)
{
    int c1, c2;

    while (*s1 != 0 && *s2 != 0) {
        c1 = g_ascii_tolower(cast(ubyte) *s1);
        c2 = g_ascii_tolower(cast(ubyte) *s2);
        if (c1 != c2) {
            return c1 - c2;
        }
        s1++;
        s2++;
    }

    return g_ascii_tolower(cast(ubyte) *s1) - g_ascii_tolower(cast(ubyte) *s2);
}

gchar * g_strdup(const(gchar) * str)
{
    size_t n;
    gchar * ret;

    if (str == null) {
        return null;
    }
    n = strlen(str) + 1;
    ret = cast(gchar *) malloc(n);
    if (ret != null) {
        memcpy(ret, str, n);
    }
    return ret;
}

gint g_vasprintf(gchar ** string, const(gchar) * format, va_list args)
{
    va_list args2;
    int len;

    va_copy(args2, args);
    len = vsnprintf(null, 0, format, args2);
    va_end(args2);

    if (len < 0) {
        *string = null;
        return -1;
    }

    *string = cast(gchar *) malloc(len + 1);
    if (*string == null) {
        return -1;
    }
    vsnprintf(*string, len + 1, format, args);
    return len;
}

extern (C) gchar * g_strdup_printf(const(gchar) * format, ...)
{
    va_list args;
    gchar * ret = null;

    va_start(args, format);
    g_vasprintf(&ret, format, args);
    va_end(args);
    return ret;
}

gint g_strcmp0(const(gchar) * str1, const(gchar) * str2)
{
    if (str1 == null) {
        return -(str1 != str2);
    }
    if (str2 == null) {
        return str1 != str2;
    }
    return strcmp(str1, str2);
}

private bool unicode_valid(uint c)
{
    return c < 0x110000 && (c & 0xFFFFF800) != 0xD800 &&
        (c < 0xFDD0 || c > 0xFDEF) && (c & 0xFFFE) != 0xFFFE;
}

gboolean g_utf8_validate(const(gchar) * str, gssize max_len, const(gchar) ** end)
{
    const(ubyte) * p = cast(const(ubyte) *) str;
    const(ubyte) * stop = max_len < 0 ? null : p + max_len;
    gboolean ok = 1;

    while ((stop == null || p < stop) && *p != 0) {
        uint c = *p;
        int len, i;
        uint min;

        if (c < 0x80) {
            p++;
            continue;
        } else if (c >= 0xC2 && c < 0xE0) {
            len = 2; c &= 0x1F; min = 0x80;
        } else if (c >= 0xE0 && c < 0xF0) {
            len = 3; c &= 0x0F; min = 0x800;
        } else if (c >= 0xF0 && c < 0xF5) {
            len = 4; c &= 0x07; min = 0x10000;
        } else {
            ok = 0;
            break;
        }
        if (stop != null && p + len > stop) {
            ok = 0;
            break;
        }
        for (i=1; i<len; i++) {
            if ((p[i] & 0xC0) != 0x80) {
                ok = 0;
                break;
            }
            c = (c << 6) | (p[i] & 0x3F);
        }
        if (!ok || c < min || !unicode_valid(c)) {
            ok = 0;
            break;
        }
        p += len;
    }

    if (end != null) {
        *end = cast(const(gchar) *) p;
    }
    return ok && (stop != null ? p == stop : *p == 0);
}

gsize g_strlcpy(gchar * dest, const(gchar) * src, gsize dest_size)
{
    gsize srclen = strlen(src);

    if (dest_size > 0) {
        gsize n = srclen < dest_size - 1 ? srclen : dest_size - 1;

        memcpy(dest, src, n);
        dest[n] = 0;
    }
    return srclen;
}

/* g_strlcat: append @src to the NUL-terminated @dest (whose total size is @dest_size), truncating
   if needed; returns the length the full result would have had (glib semantics). */

gsize g_strlcat(gchar * dest, const(gchar) * src, gsize dest_size)
{
    gsize dlen = 0;

    while (dlen < dest_size && dest[dlen] != 0) {
        dlen++;
    }
    if (dlen == dest_size) {
        return dest_size + strlen(src);
    }
    return dlen + g_strlcpy(dest + dlen, src, dest_size - dlen);
}

glong g_utf8_strlen(const(gchar) * p, gssize max)
{
    const(ubyte) * s = cast(const(ubyte) *) p;
    const(ubyte) * start = s;
    glong len = 0;

    while ((max < 0 || s - start < max) && *s != 0) {
        ubyte c = *s;
        size_t skip = c < 0xC0 ? 1 : c < 0xE0 ? 2 : c < 0xF0 ? 3 :
            c < 0xF8 ? 4 : c < 0xFC ? 5 : c < 0xFE ? 6 : 1;
        s += skip;
        len++;
    }
    return len;
}

private size_t utf8_skip(ubyte c)
{
    return c < 0xC0 ? 1 : c < 0xE0 ? 2 : c < 0xF0 ? 3 :
        c < 0xF8 ? 4 : c < 0xFC ? 5 : c < 0xFE ? 6 : 1;
}

gchar * g_utf8_next_char(const(gchar) * p)
{
    return cast(gchar *) p + utf8_skip(cast(ubyte) *p);
}

/* g_utf8_get_char: decode the character at @p (which must point to valid UTF-8, as in glib) */
uint g_utf8_get_char(const(gchar) * p)
{
    ubyte c = cast(ubyte) p[0];
    size_t len = utf8_skip(c);
    uint u;

    if (len == 1) {
        return c;
    }
    u = c & (0x7F >> len);
    foreach (i; 1 .. len) {
        u = (u << 6) | (cast(ubyte) p[i] & 0x3F);
    }
    return u;
}

gchar * g_utf8_prev_char(const(gchar) * p)
{
    do {
        p--;
    } while ((cast(ubyte) *p & 0xC0) == 0x80);
    return cast(gchar *) p;
}

gchar * g_utf8_find_next_char(const(gchar) * p, const(gchar) * end)
{
    if (*p) {
        if (end != null) {
            for (++p; p < end && (cast(ubyte) *p & 0xC0) == 0x80; ++p) {}
        } else {
            for (++p; (cast(ubyte) *p & 0xC0) == 0x80; ++p) {}
        }
    }
    return (p == end) ? null : cast(gchar *) p;
}

gchar * g_utf8_offset_to_pointer(const(gchar) * str, glong offset)
{
    const(gchar) * p = str;

    while (offset > 0) {
        p = g_utf8_next_char(p);
        offset--;
    }
    return cast(gchar *) p;
}

gchar * g_utf8_substring(const(gchar) * str, glong start_pos, glong end_pos)
{
    const(gchar) * start = g_utf8_offset_to_pointer(str, start_pos);
    const(gchar) * end = start;
    glong i;
    size_t n;
    gchar * ret;

    for (i = start_pos; i < end_pos && *end != 0; i++) {
        end = g_utf8_next_char(end);
    }
    n = end - start;
    ret = cast(gchar *) malloc(n + 1);
    if (ret != null) {
        memcpy(ret, start, n);
        ret[n] = 0;
    }
    return ret;
}

gchar ** g_strsplit(const(gchar) * string, const(gchar) * delimiter, gint max_tokens)
{
    import core.stdc.stdlib : free, realloc;

    size_t dlen = strlen(delimiter);
    gchar ** v = null;
    size_t n = 0;
    const(gchar) * s = string;
    const(gchar) * p;

    if (dlen == 0) {
        return null;
    }
    if (max_tokens < 1) {
        max_tokens = int.max;
    }

    for (;;) {
        size_t len;

        p = (n + 1 < cast(size_t) max_tokens) ? strstr(s, delimiter) : null;
        len = (p != null) ? cast(size_t) (p - s) : strlen(s);
        v = cast(gchar **) realloc(v, (n + 2) * (gchar *).sizeof);
        v[n] = cast(gchar *) malloc(len + 1);
        memcpy(v[n], s, len);
        v[n][len] = 0;
        n++;
        if (p == null) {
            break;
        }
        s = p + dlen;
    }
    v[n] = null;
    return v;
}

void g_strfreev(gchar ** str_array)
{
    import core.stdc.stdlib : free;

    if (str_array != null) {
        for (size_t i = 0; str_array[i] != null; i++) {
            free(str_array[i]);
        }
        free(str_array);
    }
}

private alias iconv_t = void *;

version (Windows) {
    /* GNU libiconv (vcpkg's iconv port) exports its functions with a "lib" prefix; nl_langinfo is in winport.d */
    import gretl.winport : nl_langinfo;

    private extern (C) {
        pragma(mangle, "libiconv_open") iconv_t iconv_open(const(char) * tocode, const(char) * fromcode);
        pragma(mangle, "libiconv") size_t iconv(iconv_t cd, char ** inbuf, size_t * inbytesleft, char ** outbuf, size_t * outbytesleft);
        pragma(mangle, "libiconv_close") int iconv_close(iconv_t cd);
    }
} else {
    private extern (C) {
        char * nl_langinfo(int item);
        iconv_t iconv_open(const(char) * tocode, const(char) * fromcode);
        size_t iconv(iconv_t cd, char ** inbuf, size_t * inbytesleft, char ** outbuf, size_t * outbytesleft);
        int iconv_close(iconv_t cd);
    }
}

/* g_convert: convert @str from character set @from_codeset to @to_codeset with iconv. As with
   g_locale_to_utf8, this is a plain converter: no fallback characters, and an invalid or
   incomplete input sequence is an error. */

gchar * g_convert(const(gchar) * opstring, gssize len, const(gchar) * to_codeset,
                  const(gchar) * from_codeset, gsize * bytes_read,
                  gsize * bytes_written, GError ** error)
{
    size_t inlen = len < 0 ? strlen(opstring) : cast(size_t) len;
    iconv_t cd = iconv_open(to_codeset, from_codeset);

    if (cd == cast(iconv_t) -1) {
        g_set_error_literal_msg(error, "Conversion from character set is not supported");
        return null;
    }

    size_t outsize = inlen * 4 + 1;
    char * out_ = cast(char *) malloc(outsize);
    char * inp = cast(char *) opstring;
    char * outp = out_;
    size_t inleft = inlen;
    size_t outleft = outsize - 1;
    size_t r = iconv(cd, &inp, &inleft, &outp, &outleft);

    iconv_close(cd);
    if (r == cast(size_t) -1) {
        free(out_);
        g_set_error_literal_msg(error, "Invalid byte sequence in conversion input");
        return null;
    }
    *outp = 0;
    if (bytes_read != null) {
        *bytes_read = inp - opstring;
    }
    if (bytes_written != null) {
        *bytes_written = outp - out_;
    }
    return out_;
}

/* g_get_charset: the current locale's character set; returns true if it is UTF-8 (glib semantics,
   with the codeset taken from nl_langinfo as in g_locale_to_utf8). */

gboolean g_get_charset(const(char) ** charset)
{
    enum CODESET = 14;   /* glibc value of CODESET in <langinfo.h> (Linux) */
    const(char) * codeset = nl_langinfo(CODESET);

    if (charset != null) {
        *charset = codeset;
    }
    return strcmp(codeset, "UTF-8") == 0 || strcmp(codeset, "utf8") == 0;
}

gchar * g_locale_to_utf8(const(gchar) * opstring, gssize len, gsize * bytes_read,
                         gsize * bytes_written, GError ** error)
{
    enum CODESET = 14;   /* glibc value of CODESET in <langinfo.h> (Linux) */
    const(char) * codeset = nl_langinfo(CODESET);

    return g_convert(opstring, len, "UTF-8", codeset, bytes_read, bytes_written, error);
}

unittest
{
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp;

    gchar * s = g_strdup("hello");
    assert(strcmp(s, "hello") == 0);
    free(s);
    assert(g_strdup(null) is null);

    s = g_strdup_printf("%s %.3f %d", "pv".ptr, 0.12345, 7);
    assert(strcmp(s, "pv 0.123 7") == 0);
    free(s);

    {
        char[8] cb = 0;

        assert(g_strlcpy(cb.ptr, "abc".ptr, cb.length) == 3);
        assert(g_strlcat(cb.ptr, "defgh".ptr, cb.length) == 8 && strcmp(cb.ptr, "abcdefg") == 0);
    }
    assert(g_utf8_get_char("A".ptr) == 'A');
    assert(g_utf8_get_char("\u00e9".ptr) == 0xe9);
    assert(g_utf8_get_char("\u20ac".ptr) == 0x20ac);
    assert(g_utf8_strlen("abc".ptr, -1) == 3);
    assert(g_utf8_strlen("héllo€".ptr, -1) == 6);
    assert(g_utf8_strlen("abcdef".ptr, 3) == 3);

    assert(g_strcmp0(null, null) == 0);
    assert(g_strcmp0(null, "a") < 0);
    assert(g_strcmp0("a", null) > 0);
    assert(g_strcmp0("a", "b") < 0);
    assert(g_utf8_validate("h\u00e9llo\u20ac".ptr, -1, null));
    assert(!g_utf8_validate("\xff".ptr, -1, null));
    assert(!g_utf8_validate("\xc0\x80".ptr, -1, null));
    assert(!g_utf8_validate("\xed\xa0\x80".ptr, -1, null));
    assert(g_utf8_validate("abcdef".ptr, 3, null));
    assert(!g_utf8_validate("\xe2\x82".ptr, -1, null));
}

unittest
{
    char[16] buf;

    buf[0 .. 8] = "  ab c \0"[0 .. 8];
    g_strstrip(buf.ptr);
    assert(buf[0 .. 5] == "ab c\0");

    buf[0 .. 4] = "   \0"[0 .. 4];
    g_strstrip(buf.ptr);
    assert(buf[0] == 0);

    assert(g_ascii_strcasecmp("Newton", "nEWTON") == 0);
    assert(g_ascii_strcasecmp("abc", "abd") < 0);
    assert(g_ascii_strcasecmp("abd", "abc") > 0);
    assert(g_ascii_strcasecmp("ab", "abc") < 0);
    assert(g_ascii_strcasecmp("", "") == 0);
}

unittest
{
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp;

    gchar ** v = g_strsplit("a,,b,c", ",", -1);
    assert(strcmp(v[0], "a") == 0 && strcmp(v[1], "") == 0 && strcmp(v[2], "b") == 0
           && strcmp(v[3], "c") == 0 && v[4] == null);
    g_strfreev(v);
    v = g_strsplit("a::b::c", "::", 2);
    assert(strcmp(v[0], "a") == 0 && strcmp(v[1], "b::c") == 0 && v[2] == null);
    g_strfreev(v);
    v = g_strsplit("", ",", -1);
    assert(strcmp(v[0], "") == 0 && v[1] == null);
    g_strfreev(v);
    v = g_strsplit("abc", ",", -1);
    assert(strcmp(v[0], "abc") == 0 && v[1] == null);
    g_strfreev(v);

    const(gchar) * u = "h\u00e9llo\u20ac";
    assert(g_utf8_next_char(u + 1) == u + 3);
    assert(g_utf8_prev_char(u + 3) == u + 1);
    assert(g_utf8_offset_to_pointer(u, 2) == u + 3);
    assert(g_utf8_find_next_char(u + 1, null) == u + 3);
    assert(g_utf8_find_next_char("a", null) !is null);
    gchar * sub = g_utf8_substring(u, 1, 3);
    assert(strcmp(sub, "\u00e9l") == 0);
    free(sub);
    sub = g_utf8_substring(u, 4, 99);
    assert(strcmp(sub, "o\u20ac") == 0);
    free(sub);
}

unittest
{
    import core.stdc.locale : LC_ALL, setlocale;
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp;
    import gretl.gerror : g_error_free;

    GError * err;
    gsize nw;

    setlocale(LC_ALL, "C");
    gchar * u = g_locale_to_utf8("plain ascii", -1, null, &nw, &err);
    assert(u !is null && strcmp(u, "plain ascii") == 0 && nw == 11);
    free(u);
    assert(g_locale_to_utf8("caf\xe9", -1, null, null, &err) is null && err !is null);
    g_error_free(err);
    err = null;

    if (setlocale(LC_ALL, "C.UTF-8") !is null) {
        u = g_locale_to_utf8("caf\u00e9", -1, null, &nw, &err);
        assert(u !is null && nw == 5);
        free(u);
        assert(g_locale_to_utf8("caf\xe9", -1, null, null, &err) is null);
        g_error_free(err);
    }
    setlocale(LC_ALL, "C");
}

/* g_utf8_strdown / g_utf8_strup: full Unicode case conversion, done with
   D's std.uni rather than glib; the result is malloc'd (as g_malloc is
   malloc in this port) and NUL-terminated. @len < 0 means @str is
   NUL-terminated. */

private gchar * utf8_case_convert(const(gchar) * str, gssize len, bool up)
{
    import std.uni : toLower, toUpper;
    import std.utf : UTFException;

    if (str is null) {
        return null;
    }
    const(char)[] src = (len < 0)? str[0 .. strlen(str)] : str[0 .. len];
    string conv;
    try {
        conv = up ? toUpper(src.idup) : toLower(src.idup);
    } catch (UTFException e) {
        conv = src.idup;
    }
    gchar * ret = cast(gchar *) malloc(conv.length + 1);
    if (ret !is null) {
        memcpy(ret, conv.ptr, conv.length);
        ret[conv.length] = '\0';
    }
    return ret;
}

gchar * g_utf8_strdown(const(gchar) * str, gssize len)
{
    return utf8_case_convert(str, len, false);
}

gchar * g_utf8_strup(const(gchar) * str, gssize len)
{
    return utf8_case_convert(str, len, true);
}

unittest {
    gchar * d = g_utf8_strdown("ÄBc€".ptr, -1);
    assert(strcmp(d, "äbc€") == 0);
    free(d);
    gchar * u = g_utf8_strup("straße x".ptr, -1);
    assert(strcmp(u, "STRASSE X") == 0 || strcmp(u, "STRAßE X") == 0);
    free(u);
}

/* g_strndup: copy of at most @n bytes of @str (stopping at a NUL), NUL-terminated; malloc'd */

gchar * g_strndup(const(gchar) * str, gsize n)
{
    if (str is null) {
        return null;
    }
    gsize len = 0;
    while (len < n && str[len] != '\0') {
        len++;
    }
    gchar * ret = cast(gchar *) malloc(len + 1);
    if (ret !is null) {
        memcpy(ret, str, len);
        ret[len] = '\0';
    }
    return ret;
}

unittest {
    gchar * t = g_strndup("hello".ptr, 3);
    assert(strcmp(t, "hel") == 0);
    free(t);
    t = g_strndup("hi".ptr, 10);
    assert(strcmp(t, "hi") == 0);
    free(t);
}

/* g_unichar_iswide: true if @c is a "wide" (East Asian full-width) character.
   glib consults a generated Unicode table; this uses the principal wide
   ranges of the same Unicode East Asian Width property (W and F), which is
   adequate for column-width calculations. */
bool g_unichar_iswide(uint c)
{
    static immutable uint[2][] ranges = [
        [0x1100, 0x115F], [0x2329, 0x232A], [0x2E80, 0x303E], [0x3041, 0x33FF],
        [0x3400, 0x4DBF], [0x4E00, 0xA4CF], [0xA960, 0xA97F], [0xAC00, 0xD7A3],
        [0xF900, 0xFAFF], [0xFE10, 0xFE19], [0xFE30, 0xFE6F], [0xFF00, 0xFF60],
        [0xFFE0, 0xFFE6], [0x1B000, 0x1B001], [0x1F200, 0x1F251], [0x20000, 0x3FFFD]
    ];

    foreach (r; ranges) {
        if (c >= r[0] && c <= r[1]) {
            return true;
        }
    }
    return false;
}
