/* Port of gretl_intl.h */

module gretl.gretl_intlh;


/* the following enumeration is organized by alphabetical order of
   English name of language: Albanian, Basque, Bulgarian, ... */
enum GretlLangCode {
    auto_ = 0,
    c,
    sq,
    eu,
    bg,
    ca,
    zh_tw,
    zh_cn,
    cs,
    fr,
    gl,
    de,
    el,
    it,
    ja,
    pl,
    pt,
    pt_br,
    ro,
    ru,
    es,
    tr,
    uk,
    max
}

/* The !ENABLE_NLS versions of these macros */
string gettext(string s) { return s; }
int UTF_WIDTH(string s, int w) { return w; }
int UTF_WIDTH(const(char) * s, int w) { return w; }
size_t TRANSLATED_WIDTH(string s) { return s.length; }
