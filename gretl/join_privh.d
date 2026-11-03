/* Port of join_priv.h: stuff that needs to be shared between csvdata.c
   and gretl_join.c. The function prototypes are in the modules that
   define them (csvdata.d), where CsvData is also defined. */

module gretl.join_privh;


import gretl.csvdata : CsvData;
import gretl.libgretlh;


struct JoinSpec {
    int ncols;
    const(char) ** colnames;
    const(char) * mdsbase;
    int * colnums;
    int * timecols;
    CsvData * c;
    DatasetRaw * tmpset;
    int wildcard;
    int auto_midas;
    int midas_pd;
    int user_tkey;
    char ** wildnames;
    char ** mdsnames;
    char ** tmpnames;
    int n_tmp;
}
