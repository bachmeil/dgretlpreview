/* Port of gretl_paths.h */

module gretl.gretl_pathsh;


import gretl.libgretlh;

enum SearchType {
    data_search,
    db_search,
    funcs_search,
    script_search,
    user_search
}

enum HelpPaths {
    cmdref,
    funcref,
    gui_help
}

enum PkgType {
    all,
    subdir,
    toplev
}

/* these are all the gretl paths which are recorded in
   the gretl config file or Windows registry entries */
struct ConfigPaths {
    char[FILENAME_MAX] gretldir;
    char[FILENAME_MAX] workdir;
    char[FILENAME_MAX] gnuplot;
    char[FILENAME_MAX] x12a;
    char[FILENAME_MAX] tramo;
    char[FILENAME_MAX] rbinpath;
    char[FILENAME_MAX] rlibpath;
    char[FILENAME_MAX] oxlpath;
    char[FILENAME_MAX] octpath;
    char[FILENAME_MAX] statapath;
    char[FILENAME_MAX] pypath;
    char[FILENAME_MAX] jlpath;
    char[FILENAME_MAX] lppath;
    char[FILENAME_MAX] mpiexec;
    char[FILENAME_MAX] mpi_hosts;
    char[128] pngfont;
    int no_dotdir;
}
