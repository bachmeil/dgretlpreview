/*
   Replacement for the one use the port makes of glib's GRegex:
   g_regex_new() + g_regex_replace() + g_regex_unref(), as used by
   gretl_regexp_replace() in strutils.c.

   Why this file exists: the C code calls glib's GRegex (PCRE). To stay
   independent of glib this module implements "replace every match of a
   regular expression in a string" on top of D's std.regex.

   Interface: one function, gregex_replace(), returning a malloc'd C string
   (free with free()/g_free()), or null on error with the message in
   *errmsg.

   Differences from glib/PCRE that callers should know about:

   - Pattern syntax is std.regex's (Perl-like): the common syntax (classes,
     groups, alternation, quantifiers, anchors, \d \w \s, non-greedy
     quantifiers, lookahead) agrees with PCRE, but some PCRE features
     (e.g. possessive quantifiers, recursion, \K) are not available. An
     unsupported or invalid pattern is reported as an error.
   - Replacement syntax follows glib's: \0 (whole match) .. \9, \g<n> and
     \g<name> for groups, and the escapes \\ \n \t \r \f \v \a \e. A group
     that did not participate in the match expands to the empty string.
     glib's case-conversion escapes (\L \U \E \l \u) and \x{...} are NOT
     supported and are reported as an error.
   - std.regex is GC-backed, so this module is a GC dependency (as is
     ghashtable.d). If a betterc build is ever wanted, replace it.

   Used by: strutils.d.
*/

module gretl.gregex;


import core.stdc.stdlib : malloc;
import core.stdc.string : memcpy;
import std.regex : Captures, matchAll, regex, RegexException;
import std.string : fromStringz;

private bool expand_replacement(ref char[] out_, string repl, Captures!string m,
                                ref string errmsg)
{
    size_t i = 0;

    while (i < repl.length) {
        char c = repl[i];

        if (c != '\\') {
            out_ ~= c;
            i++;
            continue;
        }
        i++;
        if (i >= repl.length) {
            errmsg = "regular expression replacement ends with a backslash";
            return false;
        }
        c = repl[i];
        if (c >= '0' && c <= '9') {
            size_t k = c - '0';

            if (k < m.length) {
                out_ ~= m[k];
            }
            i++;
        } else if (c == 'g') {
            size_t j = i + 1;

            if (j >= repl.length || repl[j] != '<') {
                errmsg = "invalid group reference in replacement";
                return false;
            }
            size_t close = j + 1;
            while (close < repl.length && repl[close] != '>') {
                close++;
            }
            if (close >= repl.length) {
                errmsg = "unterminated group reference in replacement";
                return false;
            }
            string name = repl[j + 1 .. close];
            bool numeric = name.length > 0;

            foreach (ch; name) {
                if (ch < '0' || ch > '9') {
                    numeric = false;
                }
            }
            if (numeric) {
                size_t k = 0;

                foreach (ch; name) {
                    k = k * 10 + (ch - '0');
                }
                if (k < m.length) {
                    out_ ~= m[k];
                }
            } else {
                try {
                    out_ ~= m[name];
                } catch (Exception e) {
                    errmsg = "unknown group name '" ~ name ~ "' in replacement";
                    return false;
                }
            }
            i = close + 1;
        } else if (c == '\\') {
            out_ ~= '\\';
            i++;
        } else if (c == 'n') {
            out_ ~= '\n';
            i++;
        } else if (c == 't') {
            out_ ~= '\t';
            i++;
        } else if (c == 'r') {
            out_ ~= '\r';
            i++;
        } else if (c == 'f') {
            out_ ~= '\f';
            i++;
        } else if (c == 'v') {
            out_ ~= '\v';
            i++;
        } else if (c == 'a') {
            out_ ~= '\a';
            i++;
        } else if (c == 'e') {
            out_ ~= '\x1b';
            i++;
        } else {
            errmsg = "unsupported escape in regular expression replacement";
            return false;
        }
    }

    return true;
}

char * gregex_replace(const(char) * orig, const(char) * pattern,
                      const(char) * repl, string * errmsg)
{
    string s = fromStringz(orig).idup;
    string r = fromStringz(repl).idup;
    char[] result;
    size_t pos = 0;
    string err;

    try {
        auto re = regex(fromStringz(pattern).idup);

        foreach (m; matchAll(s, re)) {
            size_t b = m.hit.ptr - s.ptr;

            result ~= s[pos .. b];
            if (!expand_replacement(result, r, m, err)) {
                *errmsg = err;
                return null;
            }
            pos = b + m.hit.length;
        }
    } catch (RegexException e) {
        *errmsg = e.msg;
        return null;
    }

    result ~= s[pos .. $];

    char * ret = cast(char *) malloc(result.length + 1);
    if (ret != null) {
        memcpy(ret, result.ptr, result.length);
        ret[result.length] = '\0';
    }
    return ret;
}

unittest {
    import core.stdc.stdlib : free;
    import core.stdc.string : strcmp;

    string e;
    char * p = gregex_replace("a1b22c", "[0-9]+", "#", &e);
    assert(strcmp(p, "a#b#c") == 0);
    free(p);
    p = gregex_replace("john smith", "(\\w+) (\\w+)", "\\2, \\1", &e);
    assert(strcmp(p, "smith, john") == 0);
    free(p);
    p = gregex_replace("john smith", "(?P<f>\\w+) (?P<l>\\w+)", "\\g<l>-\\g<f>-\\g<2>", &e);
    assert(strcmp(p, "smith-john-smith") == 0);
    free(p);
    p = gregex_replace("abc", "x", "y", &e);
    assert(strcmp(p, "abc") == 0);
    free(p);
    p = gregex_replace("a.b", "\\.", "\\\\", &e);
    assert(strcmp(p, "a\\b") == 0);
    free(p);
    p = gregex_replace("abc", "b", "\\U", &e);
    assert(p is null && e.length > 0);
    p = gregex_replace("abc", "(", "x", &e);
    assert(p is null);
    p = gregex_replace("aaa", "a*", "-", &e);
    assert(p !is null);
    free(p);
}
