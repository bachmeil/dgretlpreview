/* Port of lib/src/gretl_intl.c: locale handling and some character-set utilities.

   Differences from the C:
   - The Windows branches (win32 locale names, g_win32_getlocale, textdomain) are dropped.
   - ENABLE_NLS: the C compiles two sets of the locale functions. Here the NLS-enabled set is the default
     (so gretl_push_c_numeric_locale() and the others really call setlocale/localeconv), and the stub set for
     a build without NLS is selected with the D version identifier NoNLS. Translation of messages is not part
     of the port in either case: _() (libgretlh.d) is the identity, so doing_nls() is always 0.
   - The tables of the C (langs[]) are immutable arrays of LangInfo.
   - u8_to_ascii_convert() is ported exactly, including the places where the C writes a character without
     advancing the output pointer (see bugs.md #78).
   - g_file_get/set_contents() are the replacements in gfileutils.d (g_file_set_contents writes the file
     directly rather than through a temporary file). */

module gretl.gretl_intl;


import core.stdc.locale : LC_ALL, LC_COLLATE, LC_NUMERIC, lconv, localeconv, setlocale;
import core.stdc.stdio : fputs, sprintf, stderr;
import core.stdc.stdlib : free, getenv, malloc;
import core.stdc.string : strcmp, strcpy, strlen, strncmp, strrchr;

import gretl.gerror : GError, g_error_free;
import gretl.gfileutils : g_file_get_contents, g_file_set_contents;
import gretl.glib : gchar, gsize;
import gretl.gmem : g_free;
import gretl.gretl_errors : gretl_error_clear, gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_errorsh : GretlError;
import gretl.gretl_intlh : GretlLangCode;
import gretl.gretl_paths : gretl_setenv;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_new, gretl_print_steal_buffer, pprintf, pputc;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_string_table : get_built_in_string_by_name, gretl_insert_builtin_string;
import gretl.gstrfuncs : g_convert, g_strdup, g_strdup_printf, g_utf8_get_char, g_utf8_next_char, g_utf8_strlen;
import gretl.libgretlh : _;
import gretl.strutils : gretl_strdup, set_atof_point;

version (NoNLS) {
    /* stubs for NLS-disabled case */

    void set_lcnumeric(int langid, int lcnumeric)
    {
        return;
    }

    int prefer_locale_decimal()
    {
        return 0;
    }

    int test_locale(const(char) * langstr)
    {
        return 1;
    }

    int force_language(int langid)
    {
        return 1;
    }

    void gretl_push_c_numeric_locale()
    {
        return;
    }

    void gretl_pop_c_numeric_locale()
    {
        return;
    }

    int doing_nls()
    {
        return 0;
    }

    int reset_local_decpoint()
    {
        return '.';
    }

    int get_local_decpoint()
    {
        return '.';
    }
} else {
    __gshared private int numeric_c_locale_depth = 0;
    __gshared private char * numeric_locale = null;
    __gshared private int native_dot = -1;

    /**
     * gretl_push_c_numeric_locale:
     *
     * Description: Saves the current %LC_NUMERIC locale and sets it to "C".
     * This way you can safely read/write floating point numbers all in the
     * same format, using '.' as the decimal character.  You should make sure
     * that code between gretl_push_c_numeric_locale() and gretl_pop_c_numeric_locale()
     * doesn't do any setlocale calls, or locale may end up in a strange setting.
     * Also make sure to always pop the C numeric locale after you've pushed it.
     * The calls can be nested.
     **/

    void gretl_push_c_numeric_locale()
    {
        if (native_dot == -1) {
            /* not yet determined */
            lconv * lc = localeconv();

            native_dot = (*lc.decimal_point == '.');
        }

        if (native_dot == 1) {
            /* no-op */
            return;
        }

        if (numeric_c_locale_depth == 0) {
            free(numeric_locale);
            numeric_locale = gretl_strdup(setlocale(LC_NUMERIC, null));
            setlocale(LC_NUMERIC, "C");
        }

        numeric_c_locale_depth++;
    }

    /**
     * gretl_pop_c_numeric_locale:
     *
     * Description:  Restores the LC_NUMERIC locale to what it was
     * before the matching gretl_push_c_numeric_locale(). If these calls
     * were nested, then this is a no-op until we get to the most outermost
     * layer. Code in between these should not do any setlocale calls
     * to change the %LC_NUMERIC locale or things may come out very strange.
     **/

    void gretl_pop_c_numeric_locale()
    {
        if (numeric_c_locale_depth == 0) {
            /* no-op */
            return;
        }

        numeric_c_locale_depth--;

        if (numeric_c_locale_depth == 0 && numeric_locale != null) {
            setlocale(LC_NUMERIC, numeric_locale);
            free(numeric_locale);
            numeric_locale = null;
        }
    }

    /**
     * doing_nls:
     *
     * Returns: 1 if NLS translation is in effect, 0 otherwise.
     */

    __gshared private int doing_nls_called, doing_nls_nls;

    int doing_nls()
    {
        if (!doing_nls_called) {
            doing_nls_nls = (strcmp("_Open data", _("_Open data").ptr) ||
                             strcmp("Test statistic", _("Test statistic").ptr) ||
                             strcmp("annual", _("annual").ptr));
            doing_nls_called = 1;
        }

        return doing_nls_nls;
    }

    __gshared private int decpoint;

    /**
     * reset_local_decpoint:
     *
     * Uses localeconv() to determine the representation of the decimal
     * point in the current locale.
     *
     * Returns: the decimal character for the current locale.
     */

    int reset_local_decpoint()
    {
        lconv * lc = localeconv();

        if (lc == null) {
            fputs("localeconv() gave NULL!\n", stderr);
            decpoint = '.';
        } else if (lc.decimal_point == null) {
            fputs("lc->decimal_point is NULL!\n", stderr);
            decpoint = '.';
        } else {
            decpoint = *lc.decimal_point;
        }

        set_atof_point(cast(char) decpoint);

        return decpoint;
    }

    /**
     * get_local_decpoint:
     *
     * Returns: the decimal character for the current locale.
     */

    int get_local_decpoint()
    {
        if (decpoint == 0) {
            decpoint = reset_local_decpoint();
        }
        return decpoint;
    }
}

int chinese_locale()
{
    int ret = 0;

    version (NoNLS) {
    } else {
        char * loc = setlocale(LC_ALL, null);

        ret = (loc != null && !strncmp(loc, "zh", 2));
    }

    return ret;
}

int japanese_locale()
{
    int ret = 0;

    version (NoNLS) {
    } else {
        char * loc = setlocale(LC_ALL, null);

        ret = (loc != null && !strncmp(loc, "ja", 2));
    }

    return ret;
}

int east_asian_locale()
{
    int ret = 0;

    version (NoNLS) {
    } else {
        char * loc = setlocale(LC_ALL, null);

        ret = (loc != null && (!strncmp(loc, "zh", 2) ||
                               !strncmp(loc, "ja", 2)));
    }

    return ret;
}

private struct LangInfo {
    int id;
    string name;
    string code;
}

private immutable LangInfo[] langs = [
    LangInfo(GretlLangCode.auto_,  "Automatic",            null),
    LangInfo(GretlLangCode.c,      "English",              "en_US"), /* was "C" */
    LangInfo(GretlLangCode.sq,     "Albanian",             "sq_AL"),
    LangInfo(GretlLangCode.eu,     "Basque",               "eu_ES"),
    LangInfo(GretlLangCode.bg,     "Bulgarian",            "bg_BG"),
    LangInfo(GretlLangCode.ca,     "Catalan",              "ca_ES"),
    LangInfo(GretlLangCode.zh_tw,  "Chinese (Taiwan)",     "zh_TW"),
    LangInfo(GretlLangCode.zh_cn,  "Chinese (simplified)", "zh_CN"),
    LangInfo(GretlLangCode.cs,     "Czech",                "cs_CZ"),
    LangInfo(GretlLangCode.fr,     "French",               "fr_FR"),
    LangInfo(GretlLangCode.gl,     "Galician",             "gl_ES"),
    LangInfo(GretlLangCode.de,     "German",               "de_DE"),
    LangInfo(GretlLangCode.el,     "Greek",                "el_GR"),
    LangInfo(GretlLangCode.it,     "Italian",              "it_IT"),
    LangInfo(GretlLangCode.ja,     "Japanese",             "ja_JP"),
    LangInfo(GretlLangCode.pl,     "Polish",               "pl_PL"),
    LangInfo(GretlLangCode.pt,     "Portuguese",           "pt_PT"),
    LangInfo(GretlLangCode.pt_br,  "Portuguese (Brazil)",  "pt_BR"),
    LangInfo(GretlLangCode.ro,     "Romanian",             "ro_RO"),
    LangInfo(GretlLangCode.ru,     "Russian",              "ru_RU"),
    LangInfo(GretlLangCode.es,     "Spanish",              "es_ES"),
    LangInfo(GretlLangCode.tr,     "Turkish",              "tr_TR"),
    LangInfo(GretlLangCode.uk,     "Ukrainian",            "uk_UA"),
    LangInfo(GretlLangCode.max,    null,                   null)
];

/* called from settings.c */

const(char) * gretl_lang_string_from_id(int langid)
{
    int i;

    for (i=0; i<GretlLangCode.max; i++) {
        if (langid == langs[i].id) {
            return langs[i].name.ptr;
        }
    }

    return null;
}

int gretl_lang_id_from_name(const(char) * s)
{
    int i;

    if (s == null) {
        s = get_built_in_string_by_name("lang");
        for (i=GretlLangCode.c; i<GretlLangCode.max; i++) {
            if (!strcmp(s, langs[i].code.ptr)) {
                return langs[i].id;
            }
        }
    }

    if (s != null && *s != '\0') {
        for (i=0; i<GretlLangCode.max; i++) {
            if (!strcmp(s, langs[i].name.ptr)) {
                return langs[i].id;
            }
        }
    }

    return 0;
}

const(char) * gretl_lang_code_from_id(int langid)
{
    int i;

    for (i=0; i<GretlLangCode.max; i++) {
        if (langid == langs[i].id) {
            return langs[i].code.ptr;
        }
    }

    return null;
}

version (NoNLS) {
} else {
    private char * other_set_numeric(const(char) * lang)
    {
        char * set = setlocale(LC_NUMERIC, lang);

        if (set == null) {
            char[32] lfix = 0;

            sprintf(lfix.ptr, "%s.UTF-8", lang);
            set = setlocale(LC_NUMERIC, lfix.ptr);
        }

        return set;
    }

    __gshared private int decimal_use_locale;

    int prefer_locale_decimal()
    {
        return decimal_use_locale;
    }

    /* more functions conditional on NLS enabled */

    void set_lcnumeric(int langid, int lcnumeric)
    {
        /* record the user's preference */
        decimal_use_locale = lcnumeric;

        if (!lcnumeric || langid == GretlLangCode.c) {
            setlocale(LC_NUMERIC, "C");
            gretl_setenv("LC_NUMERIC", "C");
        } else {
            /* lcnumeric is selected and we're not in LANG_C */
            const(char) * lang;
            char * set = null;

            if (langid == GretlLangCode.auto_) {
                /* respect the system LANG setting */
                lang = getenv("LANG");
            } else {
                /* fake it from user preference */
                lang = gretl_lang_code_from_id(langid);
            }

            if (lang != null) {
                set = other_set_numeric(lang);
            }
            if (set == null) {
                setlocale(LC_NUMERIC, "");
                gretl_setenv("LC_NUMERIC", "");
            }
        }

        reset_local_decpoint();
    }

    private int set_locale_with_workaround(int langid, const(char) * lcode,
                                           char ** locp)
    {
        gchar * orig = g_strdup(setlocale(LC_COLLATE, null));
        gchar * lcopy = null;
        char * test;
        int err = 0;

        if (lcode != null) {
            lcopy = g_strdup_printf("%s.UTF-8", lcode);
            test = setlocale(LC_ALL, lcopy);
        } else {
            test = setlocale(LC_ALL, null);
        }

        if (test == null) {
            /* setlocale() failed */
            err = GretlError.data;
        } else if (lcopy != null) {
            /* belt and braces */
            gretl_setenv("LANG", lcopy);
        }

        g_free(orig);
        g_free(lcopy);

        if (!err && locp != null && test != null) {
            *locp = gretl_strdup(test);
        }

        return err;
    }

    private const(char) * get_setlocale_string(int i)
    {
        return gretl_lang_code_from_id(i);
    }

    /* Here we're testing the validity, on the current system, of a request
       to use a particular locale as specified by @langstr, the English name
       of the selected language as displayed in the GUI Preferences dialog
       (e.g. "German", "French").
    */

    int test_locale(const(char) * langstr)
    {
        const(char) * lcode;
        gchar * ocpy = null;
        char * orig;
        int langid, err = 0;

        langid = gretl_lang_id_from_name(langstr);
        lcode = get_setlocale_string(langid);
        orig = setlocale(LC_ALL, null);
        if (orig != null) {
            ocpy = g_strdup(orig);
        }

        gretl_error_clear();
        err = set_locale_with_workaround(langid, lcode, null);

        if (err) {
            gretl_errmsg_sprintf(_("%s: locale is not supported "
                                   ~ "on this system").ptr, lcode);
        }
        if (ocpy != null) {
            /* restore the original locale */
            setlocale(LC_ALL, ocpy);
            g_free(ocpy);
        }

        return err;
    }

    private void record_locale(char * locale)
    {
        int done = 0;
        char * lang = getenv("LANG");

        if (lang != null) {
            /* prefer using LANG */
            if (strrchr(lang, '.') == null) {
                gretl_insert_builtin_string("lang", lang);
            } else {
                char * tmp = gretl_strdup(lang);
                char * p = strrchr(tmp, '.');

                *p = '\0';
                gretl_insert_builtin_string("lang", tmp);
                free(tmp);
            }
            done = 1;
        } else if (locale != null) {
            /* use locale as fallback */
            if (strrchr(locale, '.') == null) {
                gretl_insert_builtin_string("lang", locale);
            } else {
                char * p = strrchr(locale, '.');

                *p = '\0';
            }
            gretl_insert_builtin_string("lang", locale);
            done = 1;
        }

        if (!done) {
            gretl_insert_builtin_string("lang", "unknown");
        }
    }

    int force_language(int langid)
    {
        const(char) * lcode = null;
        char * locale = null;
        int err = 0;

        if (langid == GretlLangCode.auto_) {
            /* note: avoid getting long spew from Windows */
            locale = gretl_strdup(setlocale(LC_COLLATE, null));
            goto record;
        }

        if (langid == GretlLangCode.c) {
            gretl_setenv("LANGUAGE", "english");
            gretl_setenv("LANG", "C");
            setlocale(LC_ALL, "C");
        } else {
            /* setting a specific language other than English */
            lcode = get_setlocale_string(langid);
            if (lcode != null) {
                err = set_locale_with_workaround(langid, lcode, &locale);
            }
        }

        if (langid != GretlLangCode.c) {
            lcode = gretl_lang_code_from_id(langid);
            if (lcode != null) {
                gretl_setenv("LANGUAGE", lcode);
                gretl_setenv("LANG", lcode);
            }
        }

     record:

        record_locale(locale);
        free(locale);

        return err;
    }
}

private void iso_to_ascii_translate(char * targ, const(char) * src, int latin)
{
    char * p;
    const(char) * q;

    p = targ;
    q = src;

    if (latin == 1) {
        while (*q) {
            ubyte c = *q;

            if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
                *p++ = c;
            } else if (c >= 192 && c <= 198) {
                *p++ = 'A';
            } else if (c == 199) {
                *p++ = 'C';
            } else if (c >= 200 && c <= 203) {
                *p++ = 'E';
            } else if (c >= 204 && c <= 207) {
                *p++ = 'I';
            } else if (c == 208) {
                *p++ = 'D';
            } else if (c == 209) {
                *p++ = 'N';
            } else if (c >= 210 && c <= 214) {
                *p++ = 'O';
            } else if (c == 216) {
                *p++ = 'O';
            } else if (c >= 217 && c <= 220) {
                *p++ = 'U';
            } else if (c == 221) {
                *p++ = 'Y';
            } else if (c >= 224 && c <= 230) {
                *p++ = 'a';
            } else if (c == 231) {
                *p++ = 'c';
            } else if (c >= 232 && c <= 235) {
                *p++ = 'e';
            } else if (c >= 236 && c <= 239) {
                *p++ = 'i';
            } else if (c == 240) {
                *p++ = 'd';
            } else if (c == 241) {
                *p++ = 'n';
            } else if (c >= 242 && c <= 246) {
                *p++ = 'o';
            } else if (c == 248) {
                *p++ = 'o';
            } else if (c >= 249 && c <= 252) {
                *p++ = 'u';
            } else if (c == 253) {
                *p++ = 'y';
            }
            q++;
        }
    } else if (latin == 2) {
        while (*q) {
            ubyte c = *q;

            if (c == '\t' || c == '\n' || (c >= 32 && c <= 126)) {
                *p++ = c;
            }

            if (c==161 || c==193 || c==194 || c==195 || c==196) {
                *p++ = 'A';
            }
            else if (c==198 || c==199 || c==200) {
                *p++ = 'C';
            }
            else if (c==207 || c==208) {
                *p++ = 'D';
            }
            else if (c==201 || c==202 || c==203 || c==204) {
                *p++ = 'E';
            }
            else if (c==205 || c==206) {
                *p++ = 'I';
            }
            else if (c==163 || c==165 || c==197) {
                *p++ = 'L';
            }
            else if (c==209 || c==210) {
                *p++ = 'N';
            }
            else if (c==211 || c==212 || c==213 || c==214) {
                *p++ = 'O';
            }
            else if (c==192 || c==216) {
                *p++ = 'R';
            }
            else if (c==166 || c==169 || c==170) {
                *p++ = 'S';
            }
            else if (c==171 || c==222) {
                *p++ = 'T';
            }
            else if (c==217 || c==218 || c==219 || c==220) {
                *p++ = 'U';
            }
            else if (c==221) {
                *p++ = 'Y';
            }
            else if (c==172 || c==174 || c==175) {
                *p++ = 'Z';
            }
            else if (c==177 || c==225 || c==226 || c==227 || c==228) {
                *p++ = 'a';
            }
            else if (c==230 || c==231 || c==232) {
                *p++ = 'c';
            }
            else if (c==239 || c==240) {
                *p++ = 'd';
            }
            else if (c==233 || c==234 || c==235 || c==236) {
                *p++ = 'e';
            }
            else if (c==237 || c==238) {
                *p++ = 'i';
            }
            else if (c==179 || c==181 || c==229) {
                *p++ = 'l';
            }
            else if (c==241 || c==242) {
                *p++ = 'n';
            }
            else if (c==243 || c==244 || c==245 || c==246) {
                *p++ = 'o';
            }
            else if (c==224 || c==248) {
                *p++ = 'r';
            }
            else if (c==182 || c==185 || c==186 || c==223) {
                *p++ = 's';
            }
            else if (c==187 || c==254) {
                *p++ = 't';
            }
            else if (c==249 || c==250 || c==251 || c==252) {
                *p++ = 'u';
            }
            else if (c==253) {
                *p++ = 'y';
            }
            else if (c==188 || c==190 || c==191) {
                *p++ = 'z';
            }
            q++;
        }
    }

    *p = '\0';
}

/* If @maxlen > 0 we limit the write to @targ to at most
   @maxlen bytes (excluding the terminating nul byte).
   If @sub > 0 we write this byte to @targ in place of
   UTF-8 characters that we can't represent in ASCII,
   provided they are lower than 0x0180.
*/

extern (C) char * u8_to_ascii_convert(char * targ, const(char) * src,
                                      int maxlen, char sub)
{
    int prevspace = 0;
    const(char) * q = src;
    char * p = targ;
    uint u;
    int c, skip;
    int len = 0;

    *p = '\0';

    /* If sub == 0 we assume we're doing varnames and
       so we skip all characters that are not valid in
       a gretl varname. But if sub > 0 we pass through
       all printable ASCII characters.
    */

    while (q && *q) {
        skip = 0;
        c = cast(byte) *q;      /* C char is signed */
        if (sub > 0 && ((c >= 32 && c <= 126) || c == 9 || c == 10)) {
            /* ASCII printables */
            *p++ = cast(char) c;
            q++;
        } else if (c >= 0x0030 && c <= 0x0039) {
            /* digits 0-9 */
            *p++ = cast(char) c;
            q++;
        } else if (c >= 0x0041 && c <= 0x005A) {
            /* upper-case ASCII letters */
            *p++ = cast(char) c;
            q++;
        } else if (c >= 0x0061 && c <= 0x007A) {
            /* lower-case ASCII letters */
            *p++ = cast(char) c;
            q++;
        } else if (c == 0x005F) {
            /* underscore */
            *p++ = cast(char) c;
            q++;
        } else if (c == 0x0020) {
            if (!prevspace) {
                prevspace = 1;
                *p++ = '_';
            } else {
                skip = 1;
            }
            q++;
        } else {
            /* handle Latin-1 and Latin-2, only */
            u = g_utf8_get_char(q);
            if (u >= 0x0180) {
                skip = 1; /* can't handle */
            } else if ((u >= 0x00C0 && u <= 0x00C6) || u == 0x0102 || u == 0x0104) {
                *p++ = 'A';
            } else if (u == 0x00C7 || u == 0x0106 || u == 0x010C) {
                *p++ = 'C';
            } else if ((u >= 0x00C8 && u <= 0x00CB) || u == 0x0118 || u == 0x011A) {
                *p++ = 'E';
            } else if (u >= 0x00CC && u <= 0x00CF) {
                *p++ = 'I';
            } else if (u == 0x00D0 || u == 0x010E || u == 0x0110 || u == 0x010E) {
                *p++ = 'D';
            } else if (u == 0x00D1 || u == 0x0143 || u == 0x0147) {
                *p++ = 'N';
            } else if (u == 0x00D8 || (u >= 0x00D2 && u <= 0x00D6) || u == 0x0150) {
                *p++ = 'O';
            } else if ((u >= 0x00D9 && u <= 0x00DC) || u == 0x016E || u == 0x0170) {
                *p++ = 'U';
            } else if (u == 0x00DD) {
                *p++ = 'Y';
            } else if (u == 0x00DE || u == 0x0164) {
                *p++ = 'T';
            } else if (u == 0x00DF) {
                *p++ = 's';
            } else if ((u >= 0x00E0 && u <= 0x00E6) || u == 0x0103) {
                *p++ = 'a';
            } else if (u == 0x00E7 || u == 0x0107) {
                *p++ = 'c';
            } else if ((u >= 0x00E8 && u <= 0x00EB) || u == 0x0119 || u == 0x011B) {
                *p++ = 'e';
            } else if (u >= 0x00EC && u <= 0x00EF) {
                *p++ = 'i';
            } else if (u == 0x00F0 || u == 0x0111 || u == 0x010F) {
                *p++ = 'd';
            } else if (u == 0x00F1 || u == 0x0144 || u == 0x0148) {
                *p++ = 'n';
            } else if (u == 0x00F8 || u == 0x0151 || (u >= 0x00F2 && u <= 0x00F6)) {
                *p++ = 'o';
            } else if ((u >= 0x00F9 && u <= 0x00FC) || u == 0x016F || u == 0x0171) {
                *p++ = 'u';
            } else if (u == 0x00FD || u == 0x00FF) {
                *p++ = 'y';
            } else if (u == 0x00FE || u == 0x0163) {
                *p++ = 't';
            } else if (u == 0x0141 || u == 0x013D || u == 0x0139) {
                *p++ = 'L';
            } else if (u == 0x0142 || u == 0x013E || u == 0x013A) {
                *p++ = 'l';
            } else if (u == 0x0154 || u == 0x0158) {
                *p++ = 'R';
            } else if (u == 0x0155 || u == 0x0159) {
                *p++ = 'r';
            } else if (u == 0x0160 || u == 0x015E) {
                *p++ = 'S';
            } else if (u == 0x0161 || u == 0x015F) {
                *p = 's';       /* no p++ (as in the C) */
            } else if (u == 0x0179 || u == 0x017D || u == 0x0178) {
                *p = 'Z';       /* no p++ (as in the C) */
            } else if (u == 0x017A || u == 0x017E || u == 0x017C) {
                *p = 'z';       /* no p++ (as in the C) */
            } else if (sub > 0) {
                *p = sub;       /* no p++ (as in the C) */
            } else {
                skip = 1;
            }
            q = g_utf8_next_char(q);
        }
        if (c != 0x0020) {
            prevspace = 0;
        }
        if (!skip) len++;
        if (maxlen > 0 && len == maxlen) {
            break;
        }
    }

    *p = '\0';

    return targ;
}

private char * real_iso_to_ascii(char * s, int latin)
{
    char * tmp;

    tmp = cast(char *) malloc(strlen(s) + 1);
    if (tmp == null) {
        return null;
    }

    if (latin != 1 && latin != 2) {
        /* fallback?? */
        latin = 1;
    }

    iso_to_ascii_translate(tmp, s, latin);

    strcpy(s, tmp);
    free(tmp);

    return s;
}

extern (C) char * iso_to_ascii(char * s)
{
    return real_iso_to_ascii(s, 1);
}

char * sprint_l2_to_ascii(char * targ, const(char) * s, size_t len)
{
    iso_to_ascii_translate(targ, s, 2);

    return targ;
}

extern (C) char * asciify_utf8_varname(char * s)
{
    char * tmp = cast(char *) malloc(32);

    if (tmp != null) {
        u8_to_ascii_convert(tmp, s, 31, 0);
        strcpy(s, tmp);
        free(tmp);
    }

    return s;
}

/* Convert from UTF-8 text in @s to a form suitable for
   inclusion in RTF, where non-ASCII characters are
   recoded to escaped Unicode numbering. Return the
   converted text in a newly allocated string.
*/

char * utf8_to_rtf(const(char) * s)
{
    const(char) * nextp;
    const(char) * p = s;
    short k;
    Prn * prn;
    char * ret = null;
    int err = 0;

    prn = gretl_print_new(PrnType.buffer, &err);
    if (prn == null) {
        return null;
    }

    while (*p) {
        nextp = g_utf8_next_char(p);
        if (nextp - p > 1) {
            k = cast(short) g_utf8_get_char(p);
            pprintf(prn, "\\u%d?", cast(int) k);
        } else {
            pputc(prn, *p);
        }
        p = nextp;
    }

    ret = gretl_print_steal_buffer(prn);
    gretl_print_destroy(prn);

    return ret;
}

private enum CTRLZ = 26;

private bool ascii_ctrl(int a)
{
    return a == '\t' || a == '\n' || a == '\r' || a == CTRLZ;
}

int gretl_is_ascii(const(char) * buf)
{
    int a;

    while (*buf) {
        a = cast(byte) *buf;    /* C char is signed */
        if (a > 126 || (a < 32 && !(ascii_ctrl(a)))) {
            return 0;
        }
        buf++;
    }

    return 1;
}

/* We want to print @str in a field of @width (visible) characters,
   but @str may contain multi-byte characters. In that case, determine
   the adjustment to @width that is needed to avoid underrun and
   return the adjusted value.
*/

int get_utf_width(const(char) * str, int width)
{
    /* the number of "invisible" bytes */
    int invis = cast(int) (strlen(str) - g_utf8_strlen(str, -1));

    return width + invis;
}

int get_translated_width(const(char) * str)
{
    int w = cast(int) strlen(str);

    w += w - cast(int) g_utf8_strlen(str, -1);

    return w;
}

/* utility functionality: recoding of an entire file:
   we start with a couple of static "helpers"
*/

private gchar * file_get_content(const(char) * fname,
                                 gsize * bytes,
                                 Prn * prn,
                                 int * err)
{
    GError * gerr = null;
    gchar * buf = null;
    int ok;

    ok = g_file_get_contents(fname, &buf, bytes, &gerr);

    if (ok) {
        pprintf(prn, "got content, %zu bytes\n", *bytes);
    } else {
        *err = GretlError.fopen;
        if (gerr != null) {
            gretl_errmsg_set(gerr.message);
            g_error_free(gerr);
        }
    }

    return buf;
}

private int file_set_content(const(char) * fname,
                             const(gchar) * buf,
                             gsize buflen)
{
    GError * gerr = null;
    int ok, err = 0;

    ok = g_file_set_contents(fname, buf, buflen, &gerr);

    if (!ok) {
        err = GretlError.fopen;
        if (gerr != null) {
            gretl_errmsg_set(gerr.message);
            g_error_free(gerr);
        }
    }

    return err;
}

private gchar * glib_recode_buffer(const(char) * buf,
                                   const(char) * from_set,
                                   const(char) * to_set,
                                   gsize bytes,
                                   gsize * written,
                                   int * err)
{
    gchar * trbuf = null;
    GError * gerr = null;
    gsize got = 0;

    trbuf = g_convert(buf, bytes, to_set, from_set,
                      &got, written, &gerr);

    if (gerr != null) {
        *err = GretlError.data;
        gretl_errmsg_set(gerr.message);
        g_error_free(gerr);
    }

    return trbuf;
}

/**
 * gretl_recode_file:
 * @path1: path to original file.
 * @path2: path to file to be written.
 * @from_set: the codeset of the original file.
 * @to_set: the codeset for the recoded file.
 * @prn: gretl printer (for a few comments) or NULL.
 *
 * Returns: 0 on success or non-zero code on error.
 */

int gretl_recode_file(const(char) * path1, const(char) * path2,
                      const(char) * from_set, const(char) * to_set,
                      Prn * prn)
{
    gchar * buf = null;
    gsize bytes = 0;
    int err = 0;

    /* get entire content of original file */
    buf = file_get_content(path1, &bytes, prn, &err);

    if (!err) {
        gsize written = 0;
        gchar * trbuf = glib_recode_buffer(buf, from_set, to_set,
                                           bytes, &written, &err);

        if (!err) {
            /* write recoded text to file */
            pprintf(prn, "recoded: %zu bytes\n", written);
            err = file_set_content(path2, trbuf, written);
        }
        g_free(trbuf);
    }

    g_free(buf);

    return err;
}
