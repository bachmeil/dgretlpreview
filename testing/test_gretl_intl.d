/* Smoke tests for gretl_intl.d: the ASCII-conversion helpers (Latin-1 bytes, UTF-8 with and without a
   substitution character, maxlen, the varname mode), gretl_is_ascii, get_utf_width/get_translated_width,
   utf8_to_rtf, the language tables, the numeric-locale functions in the C locale, and gretl_recode_file
   (Latin-1 -> UTF-8). No locale with a decimal comma is installed on the build machine, so the push/pop
   functions are only checked as no-ops, and set_lcnumeric/force_language are not exercised beyond "C".

   Build and run with testing/run_tests.sh. */

import core.stdc.stdio : _IOLBF, printf, setvbuf, stdout;
import core.stdc.stdlib : free;
import core.stdc.string : strcmp, strlen;

import std.file : exists, read, remove, write;

import gretl.gretl_intl;
import gretl.gretl_intlh : GretlLangCode;

int n_tests, n_failed;

void check(bool ok, string what)
{
    n_tests++;
    if (!ok) {
        n_failed++;
        printf("FAIL: %.*s\n", cast(int) what.length, what.ptr);
    } else {
        printf("ok:   %.*s\n", cast(int) what.length, what.ptr);
    }
}

string u8(string s, int maxlen, char sub)
{
    char[128] buf = 0;

    u8_to_ascii_convert(buf.ptr, (s ~ "\0").ptr, maxlen, sub);
    return buf[0 .. strlen(buf.ptr)].idup;
}

int main()
{
    setvbuf(stdout, null, _IOLBF, 4096);

    /* Latin-1 bytes -> ASCII */
    char[32] l1 = 0;
    string src = "caf\xe9 \xc0\xd1\xdc \xe7a\xf1o\n";

    l1[0 .. src.length] = src[];
    iso_to_ascii(l1.ptr);
    check(strcmp(l1.ptr, "cafe ANU cano\n") == 0, "iso_to_ascii transliterates Latin-1 and keeps ASCII");

    /* UTF-8 -> ASCII, varname mode (sub == 0): spaces become one underscore, others dropped */
    check(u8("Ab 1_x", 0, 0) == "Ab_1_x", "varname mode: letters, digits, underscore kept");
    check(u8("a   b", 0, 0) == "a_b", "varname mode: a run of spaces gives one underscore");
    check(u8("a-b+c", 0, 0) == "abc", "varname mode: other ASCII punctuation is dropped");
    check(u8("café ÀÑÜ", 0, 0) == "cafe_ANU", "varname mode: Latin-1 letters transliterated");
    check(u8("\u0142\u00f3d\u0159", 0, 0) == "lodr", "varname mode: Latin-2 letters transliterated");
    check(u8("\u0105\u015b\u017b", 0, 0) == "", "bug #79: a-ogonek, s-acute and Z-dot are not in the table");
    check(u8("中文x", 0, 0) == "x", "varname mode: CJK is skipped");

    /* sub > 0: printable ASCII passes through */
    check(u8("a-b c.d", 0, '?') == "a-b c.d", "sub mode: printable ASCII kept");
    check(u8("naïve", 0, '?') == "naive", "sub mode: Latin-1 transliterated");
    check(u8("abcdef", 3, '?') == "abc", "maxlen limits the output");

    /* the C does not advance the output pointer after writing 's', 'Z', 'z' or the substitute (bugs.md #78) */
    check(u8("šx", 0, 0) == "x", "bug #78: s-caron is written without advancing, then overwritten");
    check(u8("a©b", 0, '?') == "ab", "bug #78: the substitute character is lost");

    char[16] vn = 0;
    string vsrc = "naïve näme";

    vn[0 .. vsrc.length] = vsrc[];
    asciify_utf8_varname(vn.ptr);
    check(strcmp(vn.ptr, "naive_name") == 0, "asciify_utf8_varname");

    /* ASCII test */
    check(gretl_is_ascii("plain text\r\n\t") == 1, "gretl_is_ascii: plain text with control characters");
    check(gretl_is_ascii("café") == 0, "gretl_is_ascii: UTF-8 is not ASCII");
    check(gretl_is_ascii("bell\x07") == 0, "gretl_is_ascii: other control characters are not");

    /* widths */
    check(get_utf_width("café", 10) == 11, "get_utf_width adds the invisible bytes");
    check(get_translated_width("café") == 6, "get_translated_width");

    /* RTF */
    char * rtf = utf8_to_rtf("aéb");

    check(rtf !is null && strcmp(rtf, "a\\u233?b") == 0, "utf8_to_rtf escapes non-ASCII");
    free(rtf);

    /* language tables */
    check(strcmp(gretl_lang_string_from_id(GretlLangCode.de), "German") == 0, "gretl_lang_string_from_id");
    check(strcmp(gretl_lang_code_from_id(GretlLangCode.pt_br), "pt_BR") == 0, "gretl_lang_code_from_id");
    check(gretl_lang_code_from_id(GretlLangCode.auto_) is null, "automatic has no code");
    check(gretl_lang_id_from_name("French") == GretlLangCode.fr, "gretl_lang_id_from_name");
    check(gretl_lang_id_from_name("Klingon") == 0, "unknown language name -> 0");
    check(gretl_lang_string_from_id(99) is null, "bad language id -> NULL");

    /* locale functions in the C locale */
    check(doing_nls() == 0, "doing_nls: no translation in the port");
    check(reset_local_decpoint() == '.' && get_local_decpoint() == '.', "decimal point is '.' in the C locale");
    gretl_push_c_numeric_locale();
    gretl_push_c_numeric_locale();
    gretl_pop_c_numeric_locale();
    gretl_pop_c_numeric_locale();
    gretl_pop_c_numeric_locale();
    check(get_local_decpoint() == '.', "push/pop (nested, and an extra pop) leave the decimal point alone");
    set_lcnumeric(GretlLangCode.c, 1);
    check(prefer_locale_decimal() == 1 && get_local_decpoint() == '.', "set_lcnumeric records the preference");
    set_lcnumeric(GretlLangCode.c, 0);
    check(prefer_locale_decimal() == 0, "set_lcnumeric(.., 0)");
    check(chinese_locale() == 0 && japanese_locale() == 0 && east_asian_locale() == 0, "not an East Asian locale");

    /* recoding a file: Latin-1 -> UTF-8 */
    write("build/intl_in.txt", cast(ubyte[]) "caf\xe9\n".dup);
    int err = gretl_recode_file("build/intl_in.txt", "build/intl_out.txt", "ISO-8859-1", "UTF-8", null);

    check(err == 0, "gretl_recode_file returns 0");
    check(exists("build/intl_out.txt") && cast(string) read("build/intl_out.txt") == "café\n",
          "recoded content is UTF-8");
    err = gretl_recode_file("build/no_such_file.txt", "build/intl_out.txt", "ISO-8859-1", "UTF-8", null);
    check(err != 0, "missing input file is an error");
    remove("build/intl_in.txt");
    remove("build/intl_out.txt");

    printf("\n%d checks, %d failed\n", n_tests, n_failed);
    return n_failed ? 1 : 0;
}
