/* Port of csvdata.h: the NON_NUMERIC marker value. The prototypes are in
   csvdata.d, which defines them. (The header also declares
   normalize_join_colname(), which has no definition in csvdata.c; it is
   not ported here.) */

module gretl.csvdatah;


enum NON_NUMERIC = 1.0e99;
