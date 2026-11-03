/* Port of plugin/biff.h: MS Excel BIFF record types and the query structure
   (by Michael Meeks, from the Gnumeric excel plugin). The BIFF_* record-type #defines are the
   members of BiffRecord (prefix dropped, lower case; BIFF_1904 is date1904). A few values appear
   twice in the C (codepage, paramqry/sxext) and are kept as in the header. The BiffQuery functions
   (ms_biff_query_*) are in workbook.d, where the C defines them. */

module gretl.plugin.biffh;


import core.stdc.stdint : uint16_t, uint32_t;

import gretl.plugin.ms_oleh;

struct BiffQuery {
    ubyte ms_op;
    ubyte ls_op;
    ushort opcode;

    ubyte * data;
    int data_malloced;
    uint length;

    uint streamPos;
    MsOleStream * pos;
}

enum MsBiffVersion {
    unknown = 0,
    v2 = 2,
    v3 = 3,
    v4 = 4,
    v5 = 5,   /* Excel 5.0 */
    v7 = 7,   /* Excel 95 */
    v8 = 8    /* Excel 97 */
}

/* Pass these a BiffQuery * */
int ex_getrow(const(BiffQuery) * p)    { return MS_OLE_GET_GUINT16(p.data + 0); }
int ex_getcol(const(BiffQuery) * p)    { return MS_OLE_GET_GUINT16(p.data + 2); }
int ex_getxf(const(BiffQuery) * p)     { return MS_OLE_GET_GUINT16(p.data + 4); }
int ex_getstrlen(const(BiffQuery) * p) { return MS_OLE_GET_GUINT16(p.data + 6); }

/* Version info types as found in various Biff records */
enum MsBiffFileType {
    workbook,
    vbmodule,
    worksheet,
    chart,
    macrosheet,
    workspace,
    unknown
}

/* Cell / XF types */
enum MsBiffHidden {
    visible,
    hidden,
    very_hidden
}

enum MsBiffXfType {
    style,
    cell
}

enum MsBiffFormat {
    ms,
    lotus
}

enum MsBiffEastern {
    context,
    left_to_right,
    right_to_left
}

struct MsBiffBofData {
    MsBiffVersion version_;     /* C: "version" is a D keyword */
    MsBiffFileType type;
}

enum MsBiffMaxRows {
    v7 = 16384,
    v8 = 65536
}

/* the BIFF record types ("A long and dull list of BIFF types", by Michael Meeks) */
enum BiffRecord {
    dimensions = 0x00,  /* 2, NOT 1,10 */
    blank = 0x01,  /* 2, NOT 10 */
    number = 0x03,  /* 2, NOT 1,10 */
    label = 0x04,  /* 2 */
    boolerr = 0x05,  /* 2, NOT 10 */
    formula = 0x06,  /* 4, NOT 10 */
    string = 0x07,  /* 2 */
    row = 0x08,  /* 2 */
    bof = 0x09,  /* 8, NOT 10 */
    eof = 0x0a,  /* 0, NOT 10 */
    index = 0x0b,  /* 2, NOT 10 */
    calccount = 0x0c,  /* 0, NOT 10 */
    calcmode = 0x0d,  /* 0, NOT 10 */
    precision = 0x0e,  /* 0 */
    refmode = 0x0f,  /* 0 */
    delta = 0x10,  /* 0 */
    iteration = 0x11,  /* 0 */
    protect = 0x12,  /* 0 */
    password = 0x13,  /* 0 */
    header = 0x14,  /* 0, NOT 10 */
    footer = 0x15,  /* 0, NOT 10 */
    externcount = 0x16,  /* 0, NOT 10 */
    externsheet = 0x17,  /* 0, NOT 10 */
    name = 0x18,  /* 2, NOT 10 */
    windowprotect = 0x19,  /* 0, NOT 10 */
    verticalpagebreaks = 0x1a,  /* 0, NOT 10 */
    horizontalpagebreaks = 0x1b,  /* 0, NOT 10 */
    note = 0x1c,  /* 0, NOT 10 */
    selection = 0x1d,  /* 0, NOT 10 */
    format = 0x1e,  /* 4, NOT 10 */
    array = 0x21,  /* 2, NOT 10 */
    date1904 = 0x22,  /* 0, NOT 1,10 */
    externname = 0x23,  /* 2 */
    defaultrowheight = 0x25,  /* 2, NOT 10 */
    left_margin = 0x26,  /* 0, NOT 10 */
    right_margin = 0x27,  /* 0, NOT 10 */
    top_margin = 0x28,  /* 0 */
    bottom_margin = 0x29,  /* 0 */
    printheaders = 0x2a,  /* 0 */
    printgridlines = 0x2b,  /* 0 */
    filepass = 0x2f,  /* 0 */
    font = 0x31,  /* 2 */
    printsize = 0x33,  /* 0, Undocumented */
    table = 0x36,  /* 2 */
    continue_ = 0x3c,  /* 0, NOT 10 */
    window1 = 0x3d,  /* 0, NOT 1,10 */
    window2 = 0x3e,  /* 2, NOT 10 */
    backup = 0x40,  /* 0, NOT 10 */
    pane = 0x41,  /* 0, NOT 10 */
    codepage = 0x42,  /* DUPLICATE 42 */
    xf_old = 0x43,  /* What is this ?, NOT 10 */
    pls = 0x4d,  /* 0 */
    dcon = 0x50,  /* 0, NOT 10 */
    dconref = 0x51,  /* 0, NOT 10 */
    dconname = 0x52,  /* 0 */
    defcolwidth = 0x55,  /* 0 */
    xct = 0x59,  /* 0 */
    crn = 0x5a,  /* 0 */
    filesharing = 0x5b,  /* 0, NOT 10 */
    writeaccess = 0x5c,  /* 0 */
    obj = 0x5d,  /* 0, NOT 10 */
    uncalced = 0x5e,  /* 0 */
    saverecalc = 0x5f,  /* 0 */
    template_ = 0x60,  /* 0, NOT 1,10 */
    tab_color = 0x62,  /* 8 */
    objprotect = 0x63,  /* 0, NOT 10 */
    colinfo = 0x7d,  /* 0 */
    rk = 0x7e,  /* 0 */
    imdata = 0x7f,  /* 0 */
    guts = 0x80,  /* 0 */
    wsbool = 0x81,  /* 0 */
    gridset = 0x82,  /* 0 */
    hcenter = 0x83,  /* 0 */
    vcenter = 0x84,  /* 0 */
    boundsheet = 0x85,  /* 0 */
    writeprot = 0x86,  /* 0 */
    addin = 0x87,  /* 0 */
    edg = 0x88,  /* 0 */
    pub = 0x89,  /* 0 */
    country = 0x8c,  /* 0 */
    hideobj = 0x8d,  /* 0 */
    sort = 0x90,  /* 0 */
    sub = 0x91,  /* 0 */
    palette = 0x92,  /* 0 */
    style = 0x93,  /* 2 */
    lhrecord = 0x94,  /* 0 */
    lhngraph = 0x95,  /* 0 */
    sound = 0x96,  /* 0 */
    lpr = 0x98,  /* 0 */
    standardwidth = 0x99,  /* 0 */
    fngroupname = 0x9a,  /* 0 */
    filtermode = 0x9b,  /* 0 */
    fngroupcount = 0x9c,  /* 0 */
    autofilterinfo = 0x9d,  /* 0 */
    autofilter = 0x9e,  /* 0 */
    scl = 0xa0,  /* 0 */
    setup = 0xa1,  /* 0 */
    coordlist = 0xa9,  /* 0 NOT 1 */
    gcw = 0xab,  /* 0 NOT 1 */
    scenman = 0xae,  /* 0 NOT 1 */
    scenario = 0xaf,  /* 0 NOT 1 */
    sxview = 0xb0,  /* 0 */
    sxvd = 0xb1,  /* 0 NOT 1 */
    sxvi = 0xb2,  /* 0 NOT 1 */
    sxivd = 0xb4,  /* 0 */
    sxli = 0xb5,  /* 0 NOT 1 */
    sxpi = 0xb6,  /* 0 NOT 1 */
    docroute = 0xb8,  /* 0 NOT 1 */
    recipname = 0xb9,  /* 0 */
    shrfmla = 0xbc,  /* 0 NOT 1 */
    mulrk = 0xbd,  /* 0 */
    mulblank = 0xbe,  /* 0 NOT 1 */
    mms = 0xc1,  /* 0 */
    addmenu = 0xc2,  /* 0 */
    delmenu = 0xc3,  /* 0 */
    sxdi = 0xc5,  /* 0 */
    sxdb = 0xc6,  /* 0 */
    sxstring = 0xcd,  /* 0 */
    sxtbl = 0xd0,  /* 0 */
    sxtbrgiitm = 0xd1,  /* 0 */
    sxtbpg = 0xd2,  /* 0 */
    obproj = 0xd3,  /* 0 */
    sxidstm = 0xd5,  /* 0 */
    rstring = 0xd6,  /* 0 */
    dbcell = 0xd7,  /* 0 */
    bookbool = 0xda,  /* 0 */
    paramqry = 0xdc,  /* DUPLICATE dc */
    sxext = 0xdc,  /* DUPLICATE dc */
    scenprotect = 0xdd,  /* 0 */
    olesize = 0xde,  /* 0 */
    uddesc = 0xdf,  /* 0 */
    xf = 0xe0,  /* 0 */
    interfacehdr = 0xe1,  /* 0 */
    interfaceend = 0xe2,  /* 0 */
    sxvx = 0xe3,  /* 0 */
    mergecells = 0xe5,  /* Undocumented */
    bg_pic = 0xe9,  /* Undocumented */
    tabidconf = 0xea,  /* 0 */
    ms_o_drawing_group = 0xeb,  /* 0 */
    ms_o_drawing = 0xec,  /* 0 */
    ms_o_drawing_selection = 0xed,  /* 0 */
    phonetic = 0xef,  /* semi-Undocumented */
    sxrule = 0xf0,  /* 0 */
    sxex = 0xf1,  /* 0 */
    sxfilt = 0xf2,  /* 0 */
    sxname = 0xf6,  /* 0 */
    sxselect = 0xf7,  /* 0 */
    sxpair = 0xf8,  /* 0 */
    sxfmla = 0xf9,  /* 0 */
    sxformat = 0xfb,  /* 0 */
    sst = 0xfc,  /* 0 */
    labelsst = 0xfd,  /* 0 */
    extsst = 0xff,  /* 0 */
    sxvdex = 0x100,  /* ONLY 1 */
    sxformula = 0x103,  /* ONLY 1 */
    sxdbex = 0x122,  /* ONLY 1 */
    tabid = 0x13d,  /* ONLY 1 */
    useselfs = 0x160,  /* ONLY 1 */
    dsf = 0x161,  /* ONLY 1 */
    xl5modify = 0x162,  /* ONLY 1 */
    filesharing2 = 0x1a5,
    userdbview = 0x1a9,  /* ONLY 1 */
    usersviewbegin = 0x1aa,
    usersviewend = 0x1ab,  /* ONLY 1 */
    qsi = 0x1ad,
    supbook = 0x1ae,  /* ONLY 1 */
    prot4rev = 0x1af,  /* ONLY 1 */
    condfmt = 0x1b0,  /* ONLY 1 */
    cf = 0x1b1,  /* ONLY 1 */
    dval = 0x1b2,  /* ONLY 1 */
    dconbin = 0x1b5,  /* ONLY 1 */
    txo = 0x1b6,  /* ONLY 1 */
    refreshall = 0x1b7,  /* ONLY 1 */
    hlink = 0x1b8,  /* ONLY 1 */
    codename = 0x1ba,  /* ONLY 1, TYPO in MS Docs */
    sxfdbtype = 0x1bb,
    prot4revpass = 0x1bc,  /* ONLY 1 */
    dv = 0x1be,  /* ONLY 1 */
    xl9file = 0x1c0,  /* ONLY 1 */
    recalcid = 0x1c1,  /* ONLY 1 */
    unknown_1 = 0x810,  /* what this is */
}
