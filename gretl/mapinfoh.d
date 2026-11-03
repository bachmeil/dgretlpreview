/* Port of lib/src/mapinfo.h: the "geoplot" map state shared by mapinfo.c, plugin/geoplot.c and graphing.c. The
   enumerations lose their prefixes (MapProj keeps the projection names in lower case; MapFlags drops MAP_ and
   NaAction drops NA_) and the struct `mapinfo_` is `MapInfo`. The prototypes belong to mapinfo.c, geoplot.c
   and graphing.c. */

module gretl.mapinfoh;


import gretl.gretl_matrixh : GretlMatrixRaw;
import gretl.libgretlh;

enum MapProj {
    prj0,      /* geoplot default: quasi-Mercator */
    wgs84,     /* null projection */
    epsg3857,  /* "proper" Mercator */
    epsg2163,  /* US National Atlas Equal Area */
    epsg3035   /* Europe Equal Area */
}

enum MapFlags {
    non_std  = 1 << 0,
    is_image = 1 << 1,
    display  = 1 << 2,
    dummy    = 1 << 3,
    multi    = 1 << 4
}

enum NaAction {
    skip,     /* exclude regions with missing payload */
    outline,  /* show outlines of such areas */
    fill      /* give such regions a specific fill color */
}

struct MapInfo {
    const(char) * mapfile;   /* name of file containing polygons */
    GretlBundle * map;       /* bundle containing polygons */
    GretlMatrixRaw * zvec;      /* "payload" data to plot, or NULL */
    char * zname;            /* name of payload series */
    char * plotfile;         /* plot filename specified by caller */
    char * datfile;          /* auxiliary file: gnuplot data */
    GretlMatrixRaw * bbox;      /* plot bounding box */
    GretlMatrixRaw * zrange;    /* range of the payload data */
    GretlMatrixRaw * zvals;     /* ordered unique values (discrete data) */
    char ** zlabels;         /* labels for @zvals (discrete data) */
    int n_codes;             /* number of coded payload values */
    GretlBundle * opts;      /* options specified by caller, or NULL */
    MapFlags flags;          /* state flags */
    NaAction na_action;      /* method for handling missing values */
    MapProj proj;            /* choice of projection */
}
