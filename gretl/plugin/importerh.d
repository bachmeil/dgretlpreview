/* Port of plugin/importer.h: shared declarations for the spreadsheet importers
   (the BookFlag bits and the workbook descriptor). The macros book_numeric_dates()
   etc. are functions here. */

module gretl.plugin.importerh;


import core.stdc.stdint : uint32_t;

import gretl.libgretlh;

enum BookFlag {
    numeric_dates   = 1 << 0,
    date_base_1904  = 1 << 1,
    auto_varnames   = 1 << 2,
    time_series     = 1 << 3,
    obs_labels      = 1 << 4,
    obs_blank       = 1 << 5,
    debug_          = 1 << 6,   /* C: BOOK_DEBUG ("debug" is a D keyword) */
    data_reversed   = 1 << 7,
    top_left_empty  = 1 << 8,
    non_numeric     = 1 << 9
}

enum OffsetKind {
    col_offset,
    row_offset
}

struct WBook {
    BookFlag flags;
    int version_;       /* C: "version" is a D keyword */
    int nsheets;
    int selected;
    int col_offset, row_offset;
    char * targname;
    char ** sheetnames;
    uint32_t * byte_offsets;
    void * colspin, rowspin;    /* GUI widgets in the C: unused here */
    int * xf_list;
    int function(WBook *, int) get_min_offset;
    void * data;
}

int book_numeric_dates(const(WBook) * b) { return (b.flags & BookFlag.numeric_dates) != 0; }
int book_base_1904(const(WBook) * b)     { return (b.flags & BookFlag.date_base_1904) != 0; }
int book_auto_varnames(const(WBook) * b) { return (b.flags & BookFlag.auto_varnames) != 0; }
int book_time_series(const(WBook) * b)   { return (b.flags & BookFlag.time_series) != 0; }
int book_obs_labels(const(WBook) * b)    { return (b.flags & BookFlag.obs_labels) != 0; }
int book_debugging(const(WBook) * b)     { return (b.flags & BookFlag.debug_) != 0; }

void book_set_numeric_dates(WBook * b)   { b.flags |= BookFlag.numeric_dates; }
void book_unset_numeric_dates(WBook * b) { b.flags &= ~BookFlag.numeric_dates; }
void book_set_base_1904(WBook * b)       { b.flags |= BookFlag.date_base_1904; }
void book_set_auto_varnames(WBook * b)   { b.flags |= BookFlag.auto_varnames; }
void book_set_time_series(WBook * b)     { b.flags |= BookFlag.time_series; }
void book_set_obs_labels(WBook * b)      { b.flags |= BookFlag.obs_labels; }
void book_set_debug(WBook * b)           { b.flags |= BookFlag.debug_; }

void book_unset_obs_labels(WBook * b)    { b.flags &= ~BookFlag.obs_labels; }
void book_unset_time_series(WBook * b)   { b.flags &= ~BookFlag.time_series; }
