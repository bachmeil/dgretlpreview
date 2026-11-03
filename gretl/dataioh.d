/* Port of dataio.h */

module gretl.dataioh;


import core.stdc.stdlib : free;
import gretl.dataset : clear_datainfo;
import gretl.libgretlh;

enum GretlFileType {
    xml_data,       /* gretl XML data file (.gdt) */
    binary_data,    /* gretl_binary data file (.gdtb) */
    csv,            /* comma-separated or other plain text data */
    octave,         /* GNU octave ascii data file */
    gnumeric,       /* gnumeric workbook data */
    xls,            /* MS Excel spreadsheet data */
    xlsx,           /* MS Office Open XML spreadsheet data */
    ods,            /* Open Document Spreadsheet data */
    wf1,            /* Eviews workfile data */
    dta,            /* Stata .dta data */
    sav,            /* SPSS .sav data */
    sas,            /* SAS xport data file */
    jmulti,         /* JMulTi data file */
    data_max,       /* -- place marker -- */
    script,         /* file containing gretl commands */
    session,        /* zipped session file */
    native_db,      /* gretl database */
    native_db_www,  /* gretl database, accessed via internet */
    rats_db,        /* RATS 4.0 database */
    pcgive_db,      /* PcGive bn7/in7 pair */
    odbc,           /* Open DataBase Connectivity */
    dbnomics,       /* DB.NOMICS access */
    map,            /* shapefile or GeoJSON */
    unrecognized    /* none of the above */
}

enum DataClearCode {
    full,           /* fully clear the dataset */
    subsample       /* dataset is sub-sampled: clear partially */
}

bool SPREADSHEET_IMPORT(GretlFileType f)
{
    return (f == GretlFileType.xls ||
            f == GretlFileType.xlsx ||
            f == GretlFileType.gnumeric ||
            f == GretlFileType.ods);
}

bool OTHER_IMPORT(GretlFileType f)
{
    return (f == GretlFileType.dta ||
            f == GretlFileType.sav ||
            f == GretlFileType.sas ||
            f == GretlFileType.jmulti ||
            f == GretlFileType.octave ||
            f == GretlFileType.wf1 ||
            f == GretlFileType.map);
}

void free_datainfo(DatasetRaw * p)
{
    if (p != null) {
        clear_datainfo(p, 0);
        free(p);
    }
}

enum DBNA = -999.0; /* missing value code for gretl databases */

enum GRETL_SCALAR_DIGITS = 12;
