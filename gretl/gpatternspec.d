/*
   Replacement for the glib GPatternSpec (g_pattern_spec_new,
   g_pattern_match_string, g_pattern_spec_free), in plain D.

   Why this file exists: the C code matches series names against glob
   patterns using glib's GPatternSpec. This module provides the same
   function names and signatures so the ported call sites read like the C,
   without depending on glib.

   Semantics (as in glib): '*' matches any sequence of characters
   (including none), '?' matches any single (UTF-8) character, and every
   other character matches itself. There is no escape character and no
   character-class syntax. Matching is case sensitive and must cover the
   whole string.

   Differences from glib: glib pre-compiles the pattern into special cases
   (exact, prefix, suffix, ...). This version keeps a copy of the pattern
   and runs an iterative wildcard match on every call, which has the same
   results.

   Used by: gretl_list.d.
*/

module gretl.gpatternspec;


import core.stdc.stdlib : free, malloc;
import core.stdc.string : memcpy, strlen;

import gretl.glib : gboolean;

struct GPatternSpec {
    char * pattern;   /* private NUL-terminated copy */
    size_t len;
}

GPatternSpec * g_pattern_spec_new(const(char) * pattern)
{
    if (pattern == null) {
        return null;
    }
    GPatternSpec * ps = cast(GPatternSpec *) malloc(GPatternSpec.sizeof);
    if (ps == null) {
        return null;
    }
    ps.len = strlen(pattern);
    ps.pattern = cast(char *) malloc(ps.len + 1);
    if (ps.pattern == null) {
        free(ps);
        return null;
    }
    memcpy(ps.pattern, pattern, ps.len + 1);
    return ps;
}

void g_pattern_spec_free(GPatternSpec * pspec)
{
    if (pspec != null) {
        free(pspec.pattern);
        free(pspec);
    }
}

/* length in bytes of the UTF-8 sequence starting with lead byte c
   (1 for stray continuation or invalid bytes) */

private size_t utf8_step(const(char) * s)
{
    ubyte c = cast(ubyte) *s;
    size_t n = c < 0xc0 ? 1 : c < 0xe0 ? 2 : c < 0xf0 ? 3 : 4;
    size_t i;

    for (i = 1; i < n; i++) {
        if ((cast(ubyte) s[i] & 0xc0) != 0x80) {
            return i;   /* truncated sequence: stop at the bad byte */
        }
    }
    return n;
}

gboolean g_pattern_match_string(GPatternSpec * pspec, const(char) * str)
{
    const(char) * p = pspec.pattern;
    const(char) * s = str;
    const(char) * star = null;     /* position of last '*' in pattern */
    const(char) * mark = null;     /* position in @str to retry from */

    while (*s != '\0') {
        if (*p == '*') {
            star = p++;
            mark = s;
        } else if (*p == '?') {
            p++;
            s += utf8_step(s);
        } else if (*p != '\0' && *p == *s) {
            p++;
            s++;
        } else if (star != null) {
            p = star + 1;
            mark += utf8_step(mark);
            s = mark;
        } else {
            return 0;
        }
    }
    while (*p == '*') {
        p++;
    }
    return *p == '\0';
}

unittest {
    import std.string : toStringz;

    bool m(string pat, string s)
    {
        auto ps = g_pattern_spec_new(pat.toStringz);
        scope(exit) g_pattern_spec_free(ps);
        return g_pattern_match_string(ps, s.toStringz) != 0;
    }

    assert(m("abc", "abc"));
    assert(!m("abc", "abcd"));
    assert(!m("abc", "ab"));
    assert(m("*", ""));
    assert(m("*", "anything"));
    assert(m("x*", "x"));
    assert(m("x*", "xyz"));
    assert(!m("x*", "yx"));
    assert(m("*z", "xyz"));
    assert(!m("*z", "xyzq"));
    assert(m("a*c", "abbbc"));
    assert(m("a*b*c", "aXbYbZc"));
    assert(!m("a*b*c", "aXbYbZ"));
    assert(m("x?", "x1"));
    assert(!m("x?", "x"));
    assert(!m("x?", "x12"));
    assert(m("?", "é"));          /* one UTF-8 char */
    assert(!m("?", "éé"));
    assert(m("*?", "ab"));
    assert(m("**a", "ba"));
    assert(m("", ""));
    assert(!m("", "a"));
}
