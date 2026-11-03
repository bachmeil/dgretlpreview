/* Port of plugin/libole2/ms-ole.h: MS Office OLE2 (compound document) reading, for the Excel importer
   (from Gnumeric, by Michael Meeks and Arturo Tena; adapted for gretl by Allin Cottrell).

   The C header declares MsOle opaquely and defines MsOleStream; here MsOle (private in ms-ole.c) is
   defined too because MsOleStream points to it. glib's GArray/GList are D dynamic arrays:
   the FAT (bb), the small-block table (sb), the small-block file (sbf) and a stream's block list are
   BLP[]; the property-storage tree keeps its children in a PPS*[] array. The MS_OLE_GET_GUINT*
   macros are functions. */

module gretl.plugin.ms_oleh;


import core.stdc.stdint : uint32_t;

enum MsOleErr {
    ok,
    exist,
    invalid,
    format,
    perm,
    mem,
    space,
    notempty,
    badarg
}

enum MsOleSeek {
    set,
    cur,
    end
}

enum MsOleType {
    storage = 1,
    stream  = 2,
    root    = 5
}

alias MsOlePos = uint;
alias MsOleSPos = int;
alias BLP = uint;       /* Block pointer */

struct MsOleStat {
    MsOleType type;
    MsOlePos size;
}

struct PPS {
    int sig;
    char * name;
    PPS *[] children;
    PPS * parent;
    uint size;
    BLP start;
    MsOleType type;
    uint idx;           /* only used on write */
}

/* Structure describing an OLE file */
struct MsOle {
    int ref_count;
    bool ole_mmap;
    ubyte * mem;
    uint length;
    char mode;
    int file_des;
    BLP[] bb;           /* Big blocks status (the FAT) */
    BLP[] sb;           /* Small block status */
    BLP[] sbf;          /* The small block file */
    uint num_pps;       /* Count of number of property sets */
    PPS * pps;          /* Property Storage: the root entry, or null */
}

enum MsOleBlockKind {
    small_block,
    large_block
}

struct MsOleStream {
    MsOlePos size;
    int function(MsOleStream *, ubyte *, MsOlePos) read_copy;
    ubyte * function(MsOleStream *, MsOlePos) read_ptr;
    MsOleSPos function(MsOleStream *, MsOleSPos, MsOleSeek) lseek;
    MsOlePos function(MsOleStream *) tell;
    MsOlePos function(MsOleStream *, ubyte *, MsOlePos) write;

    /* Private. */
    MsOleBlockKind type;
    MsOle * file;
    void * pps;         /* Straight PPS */
    BLP[] blocks;       /* A list of the blocks in the file */
    bool has_blocks;    /* C: blocks != NULL (no file when false) */
    MsOlePos position;  /* Current offset into file: points to the next byte to read */
}

ubyte MS_OLE_GET_GUINT8(const(void) * p)
{
    return *(cast(const(ubyte) *) p);
}

ushort MS_OLE_GET_GUINT16(const(void) * p)
{
    const(ubyte) * b = cast(const(ubyte) *) p;

    return cast(ushort) (b[0] | (b[1] << 8));
}

uint MS_OLE_GET_GUINT32(const(void) * p)
{
    const(ubyte) * b = cast(const(ubyte) *) p;

    return cast(uint) (b[0] | (b[1] << 8) | (b[2] << 16) | (cast(uint) b[3] << 24));
}

ulong MS_OLE_GET_GUINT64(const(void) * p)
{
    const(ubyte) * b = cast(const(ubyte) *) p;

    return MS_OLE_GET_GUINT32(b) | ((cast(ulong) MS_OLE_GET_GUINT32(b + 4)) << 32);
}
