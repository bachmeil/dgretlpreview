/* Port of gretl_prn.h */

module gretl.gretl_prnh;


enum PrnType {
    stdout,
    stderr,
    file,
    buffer,
    tempfile,
    stream,
    gzfile
}

enum PrnFormat {
    txt       = 1 << 0,
    tex       = 1 << 1,
    doc       = 1 << 2,
    rtf       = 1 << 3,
    rtf_txt   = 1 << 4,
    eqn       = 1 << 5,
    selection = 1 << 6,
    csv       = 1 << 7,
    tab       = 1 << 8,
    modeltab  = 1 << 9,
    landscape = 1 << 10,
    has_minus = 1 << 11,
    xml       = 1 << 12
}
