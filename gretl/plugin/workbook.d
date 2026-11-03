/* Port of plugin/workbook.c: reading the BIFF record stream of an Excel workbook to get the sheet
   names and offsets (and the 1904-date flag and XF format list). Mostly borrowed in the C from
   Gnumeric's Excel importer, written by Michael Meeks.

   In the C this is compiled into the excel_import plugin; the BiffQuery structure and the
   record-type constants are in biffh.d. g_malloc/g_new0/g_free are malloc/calloc/free. */

module gretl.plugin.workbook;


import core.stdc.stdio : fprintf, fputs, printf, stderr;
import core.stdc.stdlib : calloc, free, malloc, realloc;
import core.stdc.string : strcmp;

import gretl.plugin.biffh;
import gretl.gstrfuncs : g_strdup;
import gretl.plugin.importerh;
import gretl.libgretlh : _;
import gretl.plugin.ms_ole : ms_ole_destroy, ms_ole_open, ms_ole_stream_close, ms_ole_stream_open_workbook;
import gretl.plugin.ms_oleh;
import gretl.strutils : tailstrip;

private struct BiffBoundsheetData {
    ushort index;
    uint streamStartPos;
    MsBiffFileType type;
    char * name;
}

BiffQuery * ms_biff_query_new(MsOleStream * ptr)
{
    BiffQuery * bq;

    if (ptr == null) return null;

    bq = cast(BiffQuery *) calloc(1, BiffQuery.sizeof);
    bq.opcode        = 0;
    bq.length        = 0;
    bq.data_malloced = 0;
    bq.pos           = ptr;

    return bq;
}

/**
 * Returns 0 if has hit end
 **/
int ms_biff_query_next(BiffQuery * bq)
{
    ubyte[4] tmp;
    int ans = 1;

    if (bq == null || bq.pos.position >= bq.pos.size) {
        return 0;
    }

    if (bq.data_malloced) {
        free(bq.data);
        bq.data = null;
        bq.data_malloced = 0;
    }

    bq.streamPos = bq.pos.position;
    if (!bq.pos.read_copy(bq.pos, tmp.ptr, 4)) {
        return 0;
    }

    bq.opcode = MS_OLE_GET_GUINT16(tmp.ptr);
    bq.length = MS_OLE_GET_GUINT16(tmp.ptr + 2);
    bq.ms_op  = cast(ubyte) (bq.opcode >> 8);
    bq.ls_op  = (bq.opcode & 0xff);

    if (bq.length > 0) {
        bq.data = bq.pos.read_ptr(bq.pos, bq.length);
    }

    if (bq.length > 0 && bq.data == null) {
        bq.data = cast(ubyte *) calloc(bq.length, 1);
        if (!bq.pos.read_copy(bq.pos, bq.data, bq.length)) {
            ans = 0;
            free(bq.data);
            bq.data = null;
            bq.length = 0;
        } else {
            bq.data_malloced = 1;
        }
    }

    if (bq.length == 0) {
        bq.data = null;
        return 1;
    }

    return ans;
}

void ms_biff_query_destroy(BiffQuery * bq)
{
    if (bq != null) {
        if (bq.data_malloced) {
            free(bq.data);
            bq.data = null;
            bq.data_malloced = 0;
        }
        free(bq);
    }
}

private bool biff_string_get_flags(const(ubyte) * ptr, bool * word_chars,
                                   bool * extended, bool * rich)
{
    ubyte header;

    header = MS_OLE_GET_GUINT8(ptr);

    if (((header & 0xf2) == 0)) {
        *word_chars = (header & 0x1) != 0;
        *extended   = (header & 0x4) != 0;
        *rich       = (header & 0x8) != 0;
        return true;
    } else {
        *word_chars = false;
        *extended   = false;
        *rich       = false;
        return false;
    }
}

private char * biff_get_text(const(ubyte) * ptr, uint len)
{
    char * ans = null;
    bool header;
    bool high_byte;
    bool ext_str;
    bool rich_str;
    uint i;

    if (len == 0) {
        return null;
    }

    ans = cast(char *) malloc(len + 2);

    header = biff_string_get_flags(ptr,
                                   &high_byte,
                                   &ext_str,
                                   &rich_str);

    if (header)   ptr += 1;
    if (rich_str) ptr += 2;
    if (ext_str)  ptr += 4;

    for (i=0; i<len; i++) {
        ushort c;
        ubyte hdr;

        if (((hdr = MS_OLE_GET_GUINT8(ptr)) & 0xf2) == 0) {
            high_byte  = (hdr & 0x1) != 0;
            ext_str    = (hdr & 0x4) != 0;
            rich_str   = (hdr & 0x8) != 0;
            ptr += 1;
            i--;
        } else if (high_byte) {
            c = MS_OLE_GET_GUINT16(ptr);
            ptr += 2;
            ans[i] = cast(char) c;
        } else {
            c = MS_OLE_GET_GUINT8(ptr);
            ptr += 1;
            ans[i] = cast(char) c;
        }
    }

    if (i > 0) {
        ans[i] = '\0';
    } else {
        fprintf(stderr, "Warning unterminated string floating\n");
    }

    return ans;
}

private BiffBoundsheetData * biff_boundsheet_data_new(BiffQuery * q, MsBiffVersion ver)
{
    BiffBoundsheetData * ans = null;
    uint startpos;

    if (ver != MsBiffVersion.v5 && ver != MsBiffVersion.v7 && ver != MsBiffVersion.v8) {
        printf("Unknown BIFF Boundsheet spec. Assuming same as Biff7\n");
        ver = MsBiffVersion.v7;
    }

    startpos = MS_OLE_GET_GUINT32(q.data);
    if (!(startpos == MS_OLE_GET_GUINT32(q.data))) {
        return null;
    }

    if (MS_OLE_GET_GUINT8(q.data + 4) == 0 &&
        ((MS_OLE_GET_GUINT8(q.data + 5)) & 0x3) == 0) {
        /* worksheet, visible */
        ans = cast(BiffBoundsheetData *) malloc(BiffBoundsheetData.sizeof);
        ans.streamStartPos = startpos;
        ans.name = biff_get_text(q.data + 7,
                                 MS_OLE_GET_GUINT8(q.data + 6));
    }

    return ans;
}

private MsBiffBofData * ms_biff_bof_data_new(BiffQuery * q)
{
    MsBiffBofData * ans = cast(MsBiffBofData *) malloc(MsBiffBofData.sizeof);

    if ((q.opcode & 0xff) == BiffRecord.bof && (q.length >= 4)) {

        switch (q.opcode >> 8) {
        case 0:
            ans.version_ = MsBiffVersion.v2;
            break;
        case 2:
            ans.version_ = MsBiffVersion.v3;
            break;
        case 4:
            ans.version_ = MsBiffVersion.v4;
            break;
        case 8: {
            switch (MS_OLE_GET_GUINT16(q.data)) {
            case 0x0600:
                ans.version_ = MsBiffVersion.v8;
                break;
            case 0x500:
                ans.version_ = MsBiffVersion.v7;
                break;
            default:
                printf("Unknown BIFF sub-number in BOF %x\n", q.opcode);
                ans.version_ = MsBiffVersion.unknown;
            }
            break;
        }

        default:
            printf("Unknown BIFF number in BOF %x\n", q.opcode);
            ans.version_ = MsBiffVersion.unknown;
            printf("Biff version %d\n", cast(int) ans.version_);
        }

        switch (MS_OLE_GET_GUINT16(q.data + 2)) {
        case 0x0005: ans.type = MsBiffFileType.workbook; break;
        case 0x0006: ans.type = MsBiffFileType.vbmodule; break;
        case 0x0010: ans.type = MsBiffFileType.worksheet; break;
        case 0x0020: ans.type = MsBiffFileType.chart; break;
        case 0x0040: ans.type = MsBiffFileType.macrosheet; break;
        case 0x0100: ans.type = MsBiffFileType.workspace; break;
        default:
            ans.type = MsBiffFileType.unknown;
            printf("Unknown BIFF type in BOF %x\n",
                   MS_OLE_GET_GUINT16(q.data + 2));
            break;
        }
    } else {
        printf("Not a BOF!\n");
        ans.version_ = MsBiffVersion.unknown;
        ans.type = MsBiffFileType.unknown;
    }

    return ans;
}

private void ms_biff_bof_data_destroy(MsBiffBofData * data)
{
    free(data);
}

private void ms_excel_read_bof(BiffQuery * q, MsBiffBofData ** version_)
{
    /* The first BOF seems to be OK, the rest lie? */
    MsBiffVersion vv = MsBiffVersion.unknown;
    MsBiffBofData * ver = *version_;

    if (ver) {
        vv = ver.version_;
        ms_biff_bof_data_destroy(ver);
    }

    *version_ = ver = ms_biff_bof_data_new(q);
    if (vv != MsBiffVersion.unknown) {
        ver.version_ = vv;
    }

    if (ver.type == MsBiffFileType.workbook) {
        if (ver.version_ >= MsBiffVersion.v8) {
            uint ver2 = MS_OLE_GET_GUINT32(q.data + 4);
            if (ver2 == 0x4107cd18) {
                printf("Excel 2000 ?\n");
            } else {
                printf("Excel 97 +\n");
            }
        } else if (ver.version_ >= MsBiffVersion.v7)
            printf("Excel 95\n");
        else if (ver.version_ >= MsBiffVersion.v5)
            printf("Excel 5.x\n");
        else if (ver.version_ >= MsBiffVersion.v4)
            printf("Excel 4.x\n");
        else if (ver.version_ >= MsBiffVersion.v3)
            printf("Excel 3.x\n");
        else if (ver.version_ >= MsBiffVersion.v2)
            printf("Excel 2.x\n");
    } else if (ver.type == MsBiffFileType.worksheet || ver.type == MsBiffFileType.chart) {
        /* nothing to do */
    } else if (ver.type == MsBiffFileType.vbmodule ||
               ver.type == MsBiffFileType.macrosheet) {
        /* printf ("VB Module or Macrosheet.\n"); */
        while (ms_biff_query_next(q) && q.opcode != BiffRecord.eof) {
            /* skip to the end of the record stream */
        }
        if (q.opcode != BiffRecord.eof)
            fprintf(stderr, "XLS file format error. Missing BIFF_EOF\n");
    } else {
        printf("Unknown BOF (%x)\n", cast(uint) ver.type);
    }
}

private int extend_xf_list(WBook * book, int fmt)
{
    int * xf_list = null;
    int nxf = 0, ret = 0;

    if (book.xf_list == null) {
        xf_list = cast(int *) malloc(2 * int.sizeof);
    } else {
        nxf = book.xf_list[0];
        xf_list = cast(int *) realloc(book.xf_list, (nxf + 2) * int.sizeof);
    }

    if (xf_list != null) {
        nxf++;
        xf_list[0] = nxf;
        xf_list[nxf] = fmt;
        book.xf_list = xf_list;
    } else {
        ret = 1;
    }

    return ret;
}

private int ms_excel_read_workbook(MsOle * file, BiffBoundsheetData *** bounds,
                                   WBook * book)
{
    MsOleStream * stream;
    MsOleErr result;
    BiffQuery * q;
    MsBiffBofData * ver = null;
    char * problem_loading = null;
    int excel_version = MsBiffVersion.unknown;

    result = ms_ole_stream_open_workbook(&stream, file);
    if (result != MsOleErr.ok) {
        return excel_version;
    }

    q = ms_biff_query_new(stream);

    while (problem_loading == null && ms_biff_query_next(q)) {

        if (0x1 == q.ms_op) {
            /* BIFF_DSF, XL9FILE, REFRESHALL, USESELFS, TABID, PROT4REV, PROT4REVPASS,
               CODENAME, SUPBOOK are recognized in the C, but nothing is done for them
               (or for unexpected tokens) */
            continue;
        }

        switch (cast(int) q.ls_op) {
        case BiffRecord.bof:
            ms_excel_read_bof(q, &ver);
            break;

        case BiffRecord.boundsheet: {
            BiffBoundsheetData * ans;

            ans = biff_boundsheet_data_new(q, ver.version_);
            if (ans != null) {
                *bounds = cast(BiffBoundsheetData **) realloc(*bounds, (book.nsheets + 1) * (**bounds).sizeof);
                (*bounds)[book.nsheets] = ans;
                book.nsheets += 1;
            }
            break;
        }

        case BiffRecord.xf:
            if (q.data != null) {
                ushort val = cast(ushort) *(q.data + 2);

                extend_xf_list(book, cast(int) val);
            }
            break;

        case BiffRecord.filepass:
            /* All records after this are encrypted */
            problem_loading = g_strdup(_("Password protected workbooks "
                                         ~ "are not supported yet.").ptr);
            break;

        case BiffRecord.date1904:
            fprintf(stderr, "Got BIFF_1904: value ");
            if (q.data != null) {
                ushort val = cast(ushort) *q.data;

                fprintf(stderr, "%d\n", cast(int) val);
                if (val) {
                    book_set_base_1904(book);
                }
            } else {
                fputs("unknown\n", stderr);
            }
            break;

        default:
            /* the other record types listed in the C (EOF, PALETTE, FONT, XF_OLD, SST,
               EXTSST, EXTERNSHEET, FORMAT, CODEPAGE, ..., DIMENSIONS, OBJ, ...) need no
               processing here */
            break;
        }
    }

    ms_biff_query_destroy(q);
    if (ver) {
        excel_version = ver.version_;
        ms_biff_bof_data_destroy(ver);
    }
    ms_ole_stream_close(&stream);

    return excel_version;
}

/* public interface */

int excel_book_get_info(const(char) * fname, WBook * book)
{
    MsOleErr ole_error;
    MsOle * f;
    BiffBoundsheetData ** bounds = null;
    int i;

    ole_error = ms_ole_open(&f, fname);

    if (ole_error != MsOleErr.ok) {
        string msg = (ole_error == MsOleErr.invalid ||
                      ole_error == MsOleErr.format) ?
            _("This file is not an 'OLE' file -- it may be too "
              ~ "old for gretl to read\n")
            : _("Unexpected error reading the file\n");
        ms_ole_destroy(&f);
        fputs(msg.ptr, stderr);
        return 1;
    }

    book.version_ = ms_excel_read_workbook(f, &bounds, book);
    ms_ole_destroy(&f);

    if (book.nsheets == 0 || bounds == null) {
        return 1;
    }

    book.sheetnames = cast(char **) malloc(book.nsheets * (char *).sizeof);
    if (book.sheetnames == null) return 1;

    book.byte_offsets = cast(uint *) malloc(book.nsheets * uint.sizeof);
    if (book.byte_offsets == null) return 1;

    for (i=0; i<book.nsheets; i++) {
        if (bounds[i].name == null) {
            book.sheetnames[i] = g_strdup("");
        } else {
            book.sheetnames[i] = bounds[i].name;
            tailstrip(book.sheetnames[i]);
        }
        book.byte_offsets[i] = bounds[i].streamStartPos;
        free(bounds[i]);
    }

    free(bounds);

    return 0;
}
