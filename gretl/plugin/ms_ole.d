/* Port of plugin/libole2/ms-ole.c: reading MS Office OLE2 files (as used for Excel .xls workbooks).
   Originally from Gnumeric (Michael Meeks, Arturo Tena), "adapted and stripped down for gretl by
   Allin Cottrell": read-only.

   Differences from the C:
   - glib's GArray/GList are D dynamic arrays (see ms_oleh.d); the children of a storage are kept in
     the order found rather than sorted by name, which only matters for duplicate names (lookup is by name);
   - the file is always memory-mapped (mmap, PROT_READ); the C's fallback of reading blocks through an LRU
     cache when mmap fails (and its BBBlkAttr machinery) is dropped, so a failing mmap is an error;
   - g_warning() calls are fprintf(stderr, ...);
   - the open file is read with gretl_open() and the descriptor kept, as in the C. */

module gretl.plugin.ms_ole;


import core.stdc.stdio : fprintf, fputs, stderr;
import core.stdc.stdlib : calloc, free, malloc;
import core.stdc.string : memcpy, strcmp, strdup;
version (Posix) {
import core.sys.posix.fcntl : O_RDONLY;
import core.sys.posix.sys.mman : MAP_FAILED, MAP_SHARED, PROT_READ, mmap, munmap;
import core.sys.posix.sys.stat : S_ISREG, fstat, stat_t;
import core.sys.posix.unistd : close;
} else import gretl.winport : O_RDONLY, MAP_FAILED, MAP_SHARED, PROT_READ, mmap, munmap, S_ISREG, fstat, stat_t, close;

import gretl.gretl_paths : gretl_open;
import gretl.plugin.ms_oleh;

private enum ADD_BBD_LIST_BLOCK = 0xfffffffc;       /* -4 */
private enum SPECIAL_BLOCK = 0xfffffffd;            /* -3 (BBD_LIST BLOCK) */
private enum END_OF_CHAIN = 0xfffffffe;             /* -2 */
private enum UNUSED_BLOCK = 0xffffffff;             /* -1 */

/* FIXME tenix laola reads this from the header */
private enum BB_BLOCK_SIZE = 512;
private enum SB_BLOCK_SIZE = 64;

private uint BLOCK_COUNT(const(MsOle) * f)
{
    return (f.length + BB_BLOCK_SIZE - 1) / BB_BLOCK_SIZE;
}

private int isregfile(int fd)
{
    stat_t st;

    if (fstat(fd, &st)) {
        return 0;
    }

    return S_ISREG(st.st_mode);
}

private int getfilesize(int fd, uint * size)
{
    stat_t st;

    if (fstat(fd, &st)) {
        return -1;
    }

    *size = cast(uint) st.st_size;

    return 0;
}

private alias PPS_IDX = uint;

private enum BB_THRESHOLD = 0x1000;

private enum PPS_ROOT_INDEX = 0;
private enum PPS_BLOCK_SIZE = 0x80;
private enum PPS_END_OF_CHAIN = 0xffffffff;

private enum PPS_SIG = 0x13579753;

private bool IS_PPS(const(PPS) * p)
{
    return p.sig == PPS_SIG;
}

/* with the file mapped, a big block is at mem + (b+1) * BB_BLOCK_SIZE */
private ubyte * BB_R_PTR(MsOle * f, BLP b)
{
    return f.mem + (cast(size_t) b + 1) * BB_BLOCK_SIZE;
}

private ubyte * GET_SB_R_PTR(MsOle * f, BLP b)
{
    return BB_R_PTR(f, f.sbf[b / (BB_BLOCK_SIZE/SB_BLOCK_SIZE)])
        + ((b % (BB_BLOCK_SIZE/SB_BLOCK_SIZE)) * SB_BLOCK_SIZE);
}

/* This is a list of big blocks which contain a flat description of all blocks
   in the file. Effectively inside these blocks is a FAT of chains of other BBs,
   so the theoretical max size = 128 BB Fat blocks, thus = 128*512*512/4 blocks
   ~= 8.4MBytes */
/* The number of Big Block Descriptor (fat) Blocks */
private uint GET_NUM_BBD_BLOCKS(const(MsOle) * f) { return MS_OLE_GET_GUINT32(f.mem + 0x2c); }
/* The block locations of the Big Block Descriptor Blocks */
private enum MAX_SIZE_BBD_LIST = 109;
private uint GET_BBD_LIST(const(MsOle) * f, uint i) { return MS_OLE_GET_GUINT32(f.mem + 0x4c + i*4); }
/* Additional Big Block Descriptor (fat) Blocks */
private enum MAX_SIZE_ADD_BBD_LIST = 127;
private uint GET_NUM_ADD_BBD_LISTS(const(MsOle) * f) { return MS_OLE_GET_GUINT32(f.mem + 0x48); }
private uint GET_FIRST_ADD_BBD_LIST(const(MsOle) * f) { return MS_OLE_GET_GUINT32(f.mem + 0x44); }

/* Get the start block of the root directory ( PPS ) chain */
private uint GET_ROOT_STARTBLOCK(const(MsOle) * f) { return MS_OLE_GET_GUINT32(f.mem + 0x30); }
/* Get the start block of the SBD chain */
private uint GET_SBD_STARTBLOCK(const(MsOle) * f) { return MS_OLE_GET_GUINT32(f.mem + 0x3c); }

/* NB it is misleading to assume that Microsofts linked lists link correctly.
   It is not the case that pps_next(f, pps_prev(f, n)) = n ! For the final list
   item there are no valid links. Cretins. */
private uint PPS_GET_NAME_LEN(const(ubyte) * p) { return MS_OLE_GET_GUINT16(p + 0x40); }
private PPS_IDX PPS_GET_PREV(const(ubyte) * p) { return MS_OLE_GET_GUINT32(p + 0x44); }
private PPS_IDX PPS_GET_NEXT(const(ubyte) * p) { return MS_OLE_GET_GUINT32(p + 0x48); }
private PPS_IDX PPS_GET_DIR(const(ubyte) * p) { return MS_OLE_GET_GUINT32(p + 0x4c); }
/* These get other interesting stuff from the PPS record */
private uint PPS_GET_STARTBLOCK(const(ubyte) * p) { return MS_OLE_GET_GUINT32(p + 0x74); }
private uint PPS_GET_SIZE(const(ubyte) * p) { return MS_OLE_GET_GUINT32(p + 0x78); }
private MsOleType PPS_GET_TYPE(const(ubyte) * p) { return cast(MsOleType) MS_OLE_GET_GUINT8(p + 0x42); }

/* Try to mark the Big Block "b" as as unused if it is marked as "c", in the
   FAT "f". */
private void TRY_MARK_UNUSED_BLOCK(BLP[] f, uint block, uint mark)
{
    if (f[block] != mark) {
        fprintf(stderr, "Tried to mark as unused the block %d which has %d\n", block, f[block]);
    } else {
        f[block] = UNUSED_BLOCK;
    }
}

/* FIXME: This needs proper unicode support ! current support is a guess */
/* Length is in bytes == 1/2 the final text length */
/* NB. Different from biff_get_text, looks like a bug ! */
private char * pps_get_text(ubyte * ptr, int length)
{
    int i;
    char * ans;
    ushort c;
    ubyte * inb;

    length = (length+1)/2;

    if (length <= 0 ||
        length > (PPS_BLOCK_SIZE/4)) {
        return null;
    }

    ans = cast(char *) calloc(length + 1, 1);

    inb = ptr;
    for (i = 0; i < length; i++) {
        c = MS_OLE_GET_GUINT16(inb);
        if (c > 0 && c < 32) {
            i--;
        } else {
            ans[i] = cast(char) c;
        }
        inb += 2;
    }

    return ans;
}

/*
 * get_next_block:
 * @f:   the file handle
 * @blk: an index into the big block fat
 *
 * Return value: the block index of the BBD block.
 */
private BLP get_next_block(MsOle * f, BLP blk, bool * err)
{
    BLP bbd = GET_BBD_LIST(f, blk / (BB_BLOCK_SIZE / 4));

    if (bbd > BLOCK_COUNT(f)) {
        *err = true;
        return 0;
    } else {
        *err = false;
    }

    return MS_OLE_GET_GUINT32(BB_R_PTR(f, bbd) +
                              4 * (blk % (BB_BLOCK_SIZE / 4)));
}

/* Builds the FAT */
private int read_bb(MsOle * f)
{
    uint numbbd;
    BLP lp;
    uint num_add_bbd_lists;
    BLP missing_lps = 0;
    BLP missing_bbds;
    uint visited_add_bbd_list;
    BLP tmp;
    BLP bbd;

    if (f == null || f.mem == null) {
        return 0;
    }

    f.bb = null;
    numbbd = GET_NUM_BBD_BLOCKS(f);

    /* Add BBD's that live in the BBD list */
    for (lp = 0; (lp < BLOCK_COUNT(f) - 1) &&
             (lp < MAX_SIZE_BBD_LIST * BB_BLOCK_SIZE / 4); lp++) {
        bool err;

        tmp = get_next_block(f, lp, &err);
        if (err)
            return 0;

        f.bb ~= tmp;
    }

    /* Add BBD's that live in the additional BBD lists */
    num_add_bbd_lists = GET_NUM_ADD_BBD_LISTS(f);
    if (num_add_bbd_lists > 0) {
        if (lp != MAX_SIZE_BBD_LIST * BB_BLOCK_SIZE / 4)
            return 0;

        visited_add_bbd_list = GET_FIRST_ADD_BBD_LIST(f);
        missing_lps = BLOCK_COUNT(f) - 1
            - MAX_SIZE_BBD_LIST*BB_BLOCK_SIZE/4;
        for (lp = 0; lp < missing_lps; lp++) {
            if ((lp != 0) && !(lp%(MAX_SIZE_ADD_BBD_LIST*
                                   (BB_BLOCK_SIZE/4)))) {
                /* This lp lives in the next add bbd list */
                visited_add_bbd_list = MS_OLE_GET_GUINT32(
                    BB_R_PTR(f, visited_add_bbd_list)
                    + 4*MAX_SIZE_ADD_BBD_LIST);
                if (visited_add_bbd_list == END_OF_CHAIN) {
                    if (lp + 1 != missing_lps) {
                        /* FIXME tenix error */
                    }
                }
            }

            /* tmp here means the number of one block that
               belongs to the fat */
            bbd = MS_OLE_GET_GUINT32(BB_R_PTR(f, visited_add_bbd_list) + 4*((lp/(BB_BLOCK_SIZE/4))%MAX_SIZE_ADD_BBD_LIST));
            tmp = MS_OLE_GET_GUINT32(BB_R_PTR(f, bbd) +
                                     4 * (lp % (BB_BLOCK_SIZE / 4)));
            f.bb ~= tmp;
        }
        /* FIXME tenix do we check if we have visited all lp's but
           there are more additional lists? */
    }

    /* Mark the bbd list blocks as unused */
    for (lp = 0; lp < (numbbd < MAX_SIZE_BBD_LIST ? numbbd : MAX_SIZE_BBD_LIST); lp++) {
        TRY_MARK_UNUSED_BLOCK(f.bb, GET_BBD_LIST(f, lp),
                              SPECIAL_BLOCK);
    }
    if (num_add_bbd_lists > 0) {
        visited_add_bbd_list = GET_FIRST_ADD_BBD_LIST(f);
        TRY_MARK_UNUSED_BLOCK(f.bb, visited_add_bbd_list,
                              ADD_BBD_LIST_BLOCK);
        missing_bbds = numbbd - MAX_SIZE_BBD_LIST;
        for (lp = 0; lp < missing_bbds; lp++) {
            if ((lp != 0) && !(lp % (MAX_SIZE_ADD_BBD_LIST))) {
                /* This lp lives in the next add bbd list */
                visited_add_bbd_list = MS_OLE_GET_GUINT32(
                    BB_R_PTR(f, visited_add_bbd_list)
                    + 4*MAX_SIZE_ADD_BBD_LIST);
                if (visited_add_bbd_list == END_OF_CHAIN) {
                    if (lp + 1 != missing_lps) {
                        /* FIXME tenix error */
                    }
                }
                TRY_MARK_UNUSED_BLOCK(f.bb,
                                      visited_add_bbd_list,
                                      ADD_BBD_LIST_BLOCK);
            }

            bbd = MS_OLE_GET_GUINT32(BB_R_PTR(f, visited_add_bbd_list) + 4*(lp%MAX_SIZE_ADD_BBD_LIST));
            TRY_MARK_UNUSED_BLOCK(f.bb, bbd, SPECIAL_BLOCK);
        }
    }

    assert(f.bb.length < BLOCK_COUNT(f));

    return 1;
}

private ubyte * get_pps_ptr(MsOle * f, PPS_IDX i)
{
    BLP blk = GET_ROOT_STARTBLOCK(f);
    int lp = i/(BB_BLOCK_SIZE/PPS_BLOCK_SIZE);

    while (lp && blk != END_OF_CHAIN) {
        if (blk == SPECIAL_BLOCK ||
            blk == UNUSED_BLOCK) {
            fprintf(stderr, "Duff block in root chain\n");
            return null;
        }
        lp--;
        blk = f.bb[blk];
    }

    if (blk == END_OF_CHAIN) {
        fprintf(stderr, "Serious error finding pps %d\n", i);
        return null;
    }

    return BB_R_PTR(f, blk) + (i%(BB_BLOCK_SIZE/PPS_BLOCK_SIZE))*PPS_BLOCK_SIZE;
}

private void pps_decode_tree(MsOle * f, PPS_IDX p, PPS * parent)
{
    PPS * pps;
    ubyte * mem;

    if (p == PPS_END_OF_CHAIN)
        return;

    pps = cast(PPS *) calloc(1, PPS.sizeof);
    pps.sig = PPS_SIG;
    mem = get_pps_ptr(f, p);
    if (!mem) {
        fprintf(stderr, "Serious directory error %d\n", p);
        f.pps = null;
        free(pps);
        return;
    }

    pps.name = pps_get_text(mem, PPS_GET_NAME_LEN(mem));
    pps.type = PPS_GET_TYPE(mem);
    pps.size = PPS_GET_SIZE(mem);
    pps.children = null;
    pps.parent = parent;
    pps.idx = 0;

    if (!pps.name) { /* Make safe */
        fputs("how odd: blank named file in directory\n", stderr);
        free(pps);
        return;
    }

    f.num_pps++;

    if (parent) {
        parent.children ~= pps;
    } else {
        f.pps = pps;
    }

    if (PPS_GET_NEXT(mem) != PPS_END_OF_CHAIN)
        pps_decode_tree(f, PPS_GET_NEXT(mem), parent);

    if (PPS_GET_PREV(mem) != PPS_END_OF_CHAIN)
        pps_decode_tree(f, PPS_GET_PREV(mem), parent);

    if (PPS_GET_DIR(mem) != PPS_END_OF_CHAIN)
        pps_decode_tree(f, PPS_GET_DIR(mem), pps);

    pps.start = PPS_GET_STARTBLOCK(mem);
}

private int read_pps(MsOle * f)
{
    PPS * pps;

    if (f == null) {
        return 0;
    }

    f.num_pps = 0;
    pps_decode_tree(f, PPS_ROOT_INDEX, null);

    if (!f.pps) {
        fprintf(stderr, "Invalid root chain\n");
        return 0;
    }

    /* Fiddle root, perhaps our get_text is broken */
    /* perhaps it is just an MS oddity in coding */
    pps = f.pps;
    if (pps.name)
        free(pps.name);
    pps.name = strdup("Root Entry");

    { /* Free up the root chain */
        BLP blk, last;

        last = blk = GET_ROOT_STARTBLOCK(f);
        while (blk != END_OF_CHAIN) {
            last = blk;
            blk = f.bb[blk];
            f.bb[last] = UNUSED_BLOCK;
        }
    }

    return 1;
}

private void destroy_pps(PPS * pps)
{
    if (pps == null) {
        return;
    }

    if (pps.name) {
        free(pps.name);
        pps.name = null;
    }
    foreach (ch; pps.children) {
        destroy_pps(ch);
    }
    pps.children = null;
    free(pps);
}

private int read_sb(MsOle * f)
{
    BLP ptr;
    int lastidx, idx;
    PPS * root;

    if (f == null || f.pps == null) {
        return 0;
    }

    root = f.pps;

    f.sbf = null;
    f.sb = null;

    /* List of big blocks in SB file */
    ptr = root.start;
    while (ptr != END_OF_CHAIN) {
        if (ptr == UNUSED_BLOCK ||
            ptr == SPECIAL_BLOCK) {
            fprintf(stderr, "Corrupt small block file: serious error, "
                    ~ "invalid block in chain\n");
            f.sbf = null;
            return 0;
        }
        f.sbf ~= ptr;
        ptr = f.bb[ptr];
    }

    /* Description of small blocks */
    lastidx = -1;
    idx = 0;
    ptr = GET_SBD_STARTBLOCK(f);

    if (f.sbf.length == 0 && ptr != END_OF_CHAIN) {
        ptr = END_OF_CHAIN;
    }

    while (ptr != END_OF_CHAIN) {
        uint lp;

        if (ptr == UNUSED_BLOCK || ptr == SPECIAL_BLOCK) {
            fprintf(stderr, "Corrupt file descriptor: serious error, "
                    ~ "invalid block in chain\n");
            f.sb = null;
            return 0;
        }
        for (lp = 0;lp<BB_BLOCK_SIZE/4;lp++) {
            BLP p = MS_OLE_GET_GUINT32(BB_R_PTR(f, ptr) + lp*4);

            f.sb ~= p;
            if (p != UNUSED_BLOCK) {
                lastidx = idx;
            }
            idx++;
        }
        ptr = f.bb[ptr];
    }

    if (lastidx > 0) {
        f.sb.length = lastidx + 1;
    }

    if (f.sbf.length * BB_BLOCK_SIZE < f.sb.length*SB_BLOCK_SIZE) {
        fprintf(stderr, "Not enough small block file for descriptors\n"
                ~ "sbf->len == %d, sb->len == %d\n", cast(int) f.sbf.length,
                cast(int) f.sb.length);
        return 0;
    }

    return 1;
}

private int ms_ole_setup(MsOle * f)
{
    if (read_bb(f) && read_pps(f) && read_sb(f)) {
        return 1;
    }

    return 0;
}

private MsOle * ms_ole_new()
{
    MsOle * f = cast(MsOle *) calloc(1, MsOle.sizeof);

    f.mem = null;
    f.length = 0;
    f.mode = 'r';

    return f;
}

/**
 * ms_ole_ref:
 * @fs: filesystem object.
 *
 * Increment by one the count of references to the filesystem.
 **/
private void ms_ole_ref(MsOle * fs)
{
    if (fs != null) {
        fs.ref_count++;
    }
}

/**
 * ms_ole_unref:
 * @fs: filesystem object.
 *
 * Decrement by one the count of references to the filesystem.
 **/

private void ms_ole_unref(MsOle * fs)
{
    if (fs != null) {
        fs.ref_count--;
    }
}

/**
 * ms_ole_open:
 * @fs: filesystem object.
 * @path: path to the filesystem-in-the file on the actual filesystem.
 *
 * Opens the filesystem-in-the-file @path and creates the filesystem object @fs.
 *
 * Return value: a #MsOleErr code.
 **/

MsOleErr ms_ole_open(MsOle ** fs, const(char) * name)
{
    MsOle * f;
    int file;

    if (!fs) {
        return MsOleErr.badarg;
    }

    f = *fs = ms_ole_new();

    /* Allin Cottrell modification: we only need to read */

    f.file_des = file = gretl_open(name, O_RDONLY, 0);
    f.ref_count = 0;

    if ((file == -1) || !isregfile(file)) {
        fprintf(stderr, "No such file '%s'\n", name);
        free(f);
        *fs = null;
        return MsOleErr.exist;
    }

    if (getfilesize(file, &(f.length))) {
        fprintf(stderr, "Couldn't get the size of file '%s'\n", name);
        close(file);
        free(f);
        *fs = null;
        return MsOleErr.exist;
    }

    if (f.length <= 0x4c) { /* Bad show */
        close(file);
        free(f);
        *fs = null;
        return MsOleErr.format;
    }

    f.mem = null;

    void * m = mmap(null, f.length, PROT_READ, MAP_SHARED, file, 0);

    if (m == null || m == MAP_FAILED) {
        fprintf(stderr, "I can't mmap that file '%s'\n", name);
        close(file);
        free(f);
        *fs = null;
        return MsOleErr.exist;
    }
    f.mem = cast(ubyte *) m;
    f.ole_mmap = true;

    if (MS_OLE_GET_GUINT32(f.mem    ) != 0xe011cfd0 ||
        MS_OLE_GET_GUINT32(f.mem + 4) != 0xe11ab1a1) {
        fprintf(stderr, "Failed OLE2 magic number %x %x\n",
                MS_OLE_GET_GUINT32(f.mem), MS_OLE_GET_GUINT32(f.mem+4));
        ms_ole_destroy(fs);
        return MsOleErr.format;
    }

    {
        ushort bbs = MS_OLE_GET_GUINT16(f.mem + 0x1e);
        ushort sbs = MS_OLE_GET_GUINT16(f.mem + 0x20);

        if ((1 << bbs) != BB_BLOCK_SIZE)
            fprintf(stderr, "Big-block-size mismatch [%d] -- expect trouble.\n", bbs);

        if ((1 << sbs) != SB_BLOCK_SIZE)
            fprintf(stderr, "Small-block-size mismatch [%d] -- expect trouble.\n", sbs);
    }

    if (!ms_ole_setup(f)) {
        fprintf(stderr, "'%s': duff file!\n", name);
        ms_ole_destroy(fs);
        return MsOleErr.format;
    }

    assert(f.bb.length < BLOCK_COUNT(f));

    return MsOleErr.ok;
}

/**
 * ms_ole_destroy:
 * @fs: filesystem object.
 *
 * Closes the filesystem @fs and truncates any free blocks.
 **/

void ms_ole_destroy(MsOle ** ptr)
{
    MsOle * f = *ptr;

    if (f) {
        if (f.ref_count != 0)
            fprintf(stderr, "Unclosed files exist on this OLE stream\n");

        if (f.mem != null && f.ole_mmap) {
            munmap(f.mem, f.length);
        }
        f.mem = null;

        destroy_pps(f.pps);
        f.pps = null;

        f.bb = null;
        f.sb = null;
        f.sbf = null;

        close(f.file_des);
        free(f);
    }

    *ptr = null;
}

private MsOlePos tell_pos(MsOleStream * s)
{
    return s.position;
}

/**
 * ms_ole_lseek:
 * @s: stream object.
 * @bytes: number of bytes to set the stream pointer.
 * @type: relative from where the stream pointer will be set.
 *
 * Set the stream pointer for @s as many as @bytes bytes according to @type.
 *
 * Return value: the new position of the stream pointer.
 **/

private MsOleSPos ms_ole_lseek(MsOleStream * s, MsOleSPos bytes, MsOleSeek type)
{
    MsOleSPos newpos;

    if (s == null) {
        return -1;
    }

    if (type == MsOleSeek.set)
        newpos = bytes;
    else if (type == MsOleSeek.cur)
        newpos = s.position + bytes;
    else
        newpos = s.size + bytes;

    if (newpos > s.size || newpos < 0) {
        fprintf(stderr, "Invalid seek\n");
        return -1;
    }

    s.position = newpos;

    return newpos;
}

private ubyte * ms_ole_read_ptr_bb(MsOleStream * s, MsOlePos length)
{
    int blklen;
    ubyte * ans;
    uint len = length;
    int blockidx = s.position / BB_BLOCK_SIZE;

    if (s == null) {
        return null;
    }

    if (!s.has_blocks || blockidx >= s.blocks.length) {
        fprintf(stderr, "Reading from NULL file\n");
        return null;
    }

    blklen = BB_BLOCK_SIZE - s.position % BB_BLOCK_SIZE;

    if (len > blklen && !s.file.ole_mmap) {
        return null;
    }

    while (len > blklen) {
        len -= blklen;
        blklen = BB_BLOCK_SIZE;
        if (blockidx >= (cast(int) s.blocks.length - 1)
            || (s.blocks[blockidx] != blockidx + 1))
            return null;
        blockidx++;
    }

    /* Straight map, simply return a pointer */
    ans = BB_R_PTR(s.file, s.blocks[s.position / BB_BLOCK_SIZE])
        + s.position % BB_BLOCK_SIZE;
    ms_ole_lseek(s, length, MsOleSeek.cur);

    return ans;
}

private ubyte * ms_ole_read_ptr_sb(MsOleStream * s, MsOlePos length)
{
    int blklen;
    ubyte * ans;
    uint len = length;
    int blockidx = s.position / SB_BLOCK_SIZE;

    if (s == null) {
        return null;
    }

    if (!s.has_blocks || blockidx >= s.blocks.length) {
        fprintf(stderr, "Reading from NULL file\n");
        return null;
    }

    blklen = SB_BLOCK_SIZE - s.position % SB_BLOCK_SIZE;

    if (len > blklen && !s.file.ole_mmap) {
        return null;
    }

    while (len > blklen) {
        len -= blklen;
        blklen = SB_BLOCK_SIZE;
        if (blockidx >= (cast(int) s.blocks.length - 1)
            || (s.blocks[blockidx] != blockidx + 1))
            return null;
        blockidx++;
    }

    /* Straight map, simply return a pointer */
    ans = GET_SB_R_PTR(s.file, s.blocks[s.position / SB_BLOCK_SIZE])
        + s.position%SB_BLOCK_SIZE;
    ms_ole_lseek(s, length, MsOleSeek.cur);

    return ans;
}

/*
 *  Returns:
 *  zero    - on error
 *  no zero - on success
 */

private int ms_ole_read_copy_bb(MsOleStream * s, ubyte * ptr, MsOlePos length)
{
    ubyte * src;
    int offset = s.position % BB_BLOCK_SIZE;
    int blkidx = s.position / BB_BLOCK_SIZE;

    if (s == null || ptr == null) {
        return 0;
    }

    if (!s.has_blocks) {
        fprintf(stderr, "Reading from NULL file\n");
        return 0;
    }

    while (length > 0) {
        BLP block;
        int cpylen = BB_BLOCK_SIZE - offset;

        if (cpylen > length)
            cpylen = length;

        if (s.position + cpylen > s.size || blkidx == s.blocks.length) {
            return 0;
        }
        assert(blkidx < s.blocks.length);
        block = s.blocks[blkidx];
        src = BB_R_PTR(s.file, block) + offset;

        memcpy(ptr, src, cpylen);
        ptr += cpylen;
        length -= cpylen;

        offset = 0;

        blkidx++;
        s.position += cpylen;
    }

    return 1;
}

/*
 *  Returns:
 *  zero    - on error
 *  no zero - on success
 */

private int ms_ole_read_copy_sb(MsOleStream * s, ubyte * ptr, MsOlePos length)
{
    int offset = s.position % SB_BLOCK_SIZE;
    int blkidx = s.position / SB_BLOCK_SIZE;
    ubyte * src;

    if (s == null || ptr == null) {
        return 0;
    }

    if (!s.has_blocks) {
        fprintf(stderr, "Reading from NULL file\n");
        return 0;
    }

    while (length > 0) {
        int cpylen = SB_BLOCK_SIZE - offset;
        BLP block;

        if (cpylen>length)
            cpylen = length;
        if (s.position + cpylen > s.size || blkidx == s.blocks.length) {
            return 0;
        }
        assert(blkidx < s.blocks.length);
        block = s.blocks[blkidx];
        src = GET_SB_R_PTR(s.file, block) + offset;

        memcpy(ptr, src, cpylen);
        ptr += cpylen;
        length -= cpylen;

        offset = 0;

        blkidx++;
        s.position += cpylen;
    }

    return 1;
}

/**
 * find_in_pps:
 * @l: the parent storage.
 *
 * Find the right Stream ... John 4:13-14 ...
 * in a storage
 *
 * Return value: %NULL if not found or pointer to the child
 **/

private PPS * find_in_pps(PPS * pps, const(char) * name)
{
    if (!IS_PPS(pps)) {
        return null;
    }

    foreach (cur; pps.children) {
        if (!IS_PPS(cur)) {
            return null;
        }

        if (cur.name && !strcmp(cur.name, name)) {
            return cur;
        }
    }

    return null;
}

private MsOleErr path_to_workbook(PPS ** pps, MsOle * f)
{
    PPS * cur = f.pps;

    if (!cur)
        return MsOleErr.exist;

    cur = find_in_pps(f.pps, "Workbook");

    if (!cur) {
        cur = find_in_pps(f.pps, "Book");
    }

    if (!cur) {
        return MsOleErr.exist;
    }

    *pps = cur;
    if (!IS_PPS(cur)) {
        return MsOleErr.invalid;
    }
    return MsOleErr.ok;
}

/**
 * real_ms_ole_open_workbook:
 * @stream: stream object.
 * @fs: filesystem object.
 *
 * Opens the stream in "/" with the name "Workbook" (or failing
 * that, "Book" and creates the stream object @stream.  The open
 * is read-only.
 *
 * Return value: a #MsOleErr code.
 **/

private MsOleErr real_ms_ole_open_workbook(MsOleStream ** stream, MsOle * f)
{
    PPS * p;
    MsOleStream * s;
    int lp, panic = 0;
    MsOleErr result;

    if (!stream)
        return MsOleErr.badarg;

    *stream = null;

    if ((result = path_to_workbook(&p, f)) != MsOleErr.ok)
        return result;

    if (p.type != MsOleType.stream)
        return MsOleErr.invalid;

    s = cast(MsOleStream *) calloc(1, MsOleStream.sizeof);
    s.file = f;
    s.pps = p;
    s.position = 0;
    s.size = p.size;
    s.blocks = null;
    s.has_blocks = false;

    if (s.size >= BB_THRESHOLD) {
        BLP b = p.start;

        s.read_copy = &ms_ole_read_copy_bb;
        s.read_ptr = &ms_ole_read_ptr_bb;
        s.lseek = &ms_ole_lseek;
        s.tell = &tell_pos;
        s.write = null;

        s.has_blocks = true;
        s.type = MsOleBlockKind.large_block;
        for (lp = 0; !panic && (lp < (s.size + BB_BLOCK_SIZE - 1) / BB_BLOCK_SIZE); lp++) {
            s.blocks ~= b;
            if (b == END_OF_CHAIN ||
                b == SPECIAL_BLOCK ||
                b == UNUSED_BLOCK) {

                fprintf(stderr, "Panic: broken stream, truncating to block %d\n", lp);
                s.size = (lp-1)*BB_BLOCK_SIZE;
                panic = 1;
            } else {
                b = f.bb[b];
            }
        }
        if (b != END_OF_CHAIN) {
            BLP next;

            fprintf(stderr, "Panic: extra unused blocks on end of '%s', %x wiping it\n",
                    p.name, b);
            while (b != END_OF_CHAIN &&
                   b != UNUSED_BLOCK &&
                   b != SPECIAL_BLOCK &&
                   b < f.bb.length) {
                next = f.bb[b];
                f.bb[b] = END_OF_CHAIN;
                b = next;
            }
        }
    } else {
        BLP b = p.start;

        s.read_copy = &ms_ole_read_copy_sb;
        s.read_ptr = &ms_ole_read_ptr_sb;
        s.lseek = &ms_ole_lseek;
        s.tell = &tell_pos;
        s.write = null;

        s.has_blocks = s.size > 0;

        s.type = MsOleBlockKind.small_block;

        for (lp = 0; !panic && (lp < (s.size + SB_BLOCK_SIZE - 1) / SB_BLOCK_SIZE); lp++) {
            s.blocks ~= b;
            if (b == END_OF_CHAIN ||
                b == SPECIAL_BLOCK ||
                b == UNUSED_BLOCK) {

                fprintf(stderr, "Panic: broken stream, truncating to block %d\n", lp);
                s.size = (lp-1) * SB_BLOCK_SIZE;
                panic = 1;
            } else {
                b = f.sb[b];
            }
        }
        if (b != END_OF_CHAIN) {
            BLP next;

            fprintf(stderr, "Panic: extra unused blocks on end of '%s', wiping it\n",
                    p.name);
            while (b != END_OF_CHAIN &&
                   b != UNUSED_BLOCK &&
                   b != SPECIAL_BLOCK &&
                   b < f.sb.length) {
                next = f.sb[b];
                f.sb[b] = END_OF_CHAIN;
                b = next;
            }
            if (b != END_OF_CHAIN)
                fprintf(stderr, "Panic: even more serious block error\n");
        }
    }

    *stream = s;
    ms_ole_ref(s.file);

    return MsOleErr.ok;
}

MsOleErr ms_ole_stream_open_workbook(MsOleStream ** stream, MsOle * f)
{
    MsOleErr result;

    result = real_ms_ole_open_workbook(stream, f);

    if (result != MsOleErr.ok) {
        ms_ole_stream_close(stream);
        fputs("No Workbook or Book streams found\n", stderr);
    }

    return result;
}

/**
 * ms_ole_stream_close:
 * @stream: stream object to be closed.
 *
 * Closes the @stream.
 *
 * Return value: a #MsOleErr code.
 **/

MsOleErr ms_ole_stream_close(MsOleStream ** s)
{
    if (*s) {
        (*s).blocks = null;

        ms_ole_unref((*s).file);

        free(*s);
        *s = null;

        return MsOleErr.ok;
    }

    return MsOleErr.badarg;
}

MsOlePos ms_ole_stream_position(const(MsOleStream) * s)
{
    return s.position;
}
