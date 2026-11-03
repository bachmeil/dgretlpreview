/* Port of lib/src/gretl_gridplot.c: the "gridplot" command: collecting several plots (gnuplot buffers) into one
   multi-panel figure, via a block of commands or an array of plot specifications. The header
   gretl_gridplot.h has only prototypes. */

module gretl.gretl_gridplot;


import core.stdc.ctype : isdigit, isspace;
import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : round, nearbyint, atan2, cos, sin, ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : fread, fwrite, FILE, fclose, fopen, fgets, fputc, fputs, rewind, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atof, atoi, calloc, free, getenv, malloc, qsort, realloc;
import core.stdc.string : memcpy, memset, strcmp, strcpy, strlen, strchr, strstr, strcat, strncat, strncmp, strspn, strcspn, strtok, memmove, strrchr;

import std.algorithm : max, min;

import gretl.gretl_errorsh : GretlError;
import gretl.gretl_matrix;
import gretl.gretl_matrixh;
import gretl.johansenh : JohansenCode, JohansenInfo, jcode, jrank, auto_restr;
import gretl.libgretlh;
import gretl.missingh : NADBL, na;
import gretl.varh;

import gretl.xmlh;


private alias XUC = const(xmlChar) *;




import core.stdc.config : c_long;

/* the C library's strtol takes a plain char ** end pointer */
extern (C) c_long strtol(const(char) * nptr, char ** endptr, int base);
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gstrfuncs : g_strdup;
import gretl.plot_privh : BPMode, GPNA, GnuplotInfo;
import gretl.mapinfoh : MapFlags, MapInfo, MapProj, NaAction;
version (Posix) import core.sys.posix.sys.wait : WEXITSTATUS, WIFEXITED, waitpid; else import gretl.winport : WEXITSTATUS, WIFEXITED, waitpid;
version (Posix) import core.sys.posix.unistd : close, read, write; else import gretl.winport : close, read, write;
enum PI = 3.14159265358979323846;
import core.stdc.config : c_ulong;
extern (C) c_ulong strtoul(const(char) * nptr, char ** endptr, int base);
import gretl.dataio : check_identifier, detect_filetype;
import gretl.dataioh : GretlFileType;
import gretl.genmain : current_series_index;
import gretl.gfileutils : g_file_get_contents;
import gretl.glib : gsize;
import gretl.gmem : g_free;
import gretl.graphing : finalize_plot_input_file, open_plot_input_file, set_special_font_size, set_special_plot_size;
import gretl.graphingh : GptFlags, PlotType;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.options : clear_stored_options_for_command, get_optval_int, get_optval_string, incompatible_options;
import gretl.printout : bufgets_finalize, bufgets_init, safe_bufgets;
import gretl.strutils : gretl_strdup, has_suffix;
import gretl.gretl_array : get_array_by_name, get_strings_array_by_name, gretl_array_append_string, gretl_array_destroy, gretl_array_get_data, gretl_array_get_length, gretl_array_new, gretl_array_set_data;
import gretl.uservar : gretl_is_user_var, user_var_add_or_replace;
import gretl.gretl_commands : gretl_command_word;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}

struct MplotOption {
    Gretlopt flag;
    int * target;
    int min;
    int max;
    int def;
}

__gshared GretlArray * mp_array = null;
__gshared char * array_name = null;
__gshared int mp_fontsize = 10;
__gshared int mp_width = 800;
__gshared int mp_height = 600;
__gshared int mp_rows = 0;
__gshared int mp_cols = 0;
__gshared int mp_collecting = 0;

/* called in graphing.c */

int gretl_gridplot_collecting()
{
    return mp_array != null && mp_collecting;
}

__gshared MplotOption[] mp_options = [

    MplotOption(Gretlopt.OPT_F, &mp_fontsize, 4, 24, 10),
    MplotOption(Gretlopt.OPT_W, &mp_width, 200, 2048, 800),
    MplotOption(Gretlopt.OPT_H, &mp_height, 200, 2048, 600),
    MplotOption(Gretlopt.OPT_R, &mp_rows, 1, 12, 0),
    MplotOption(Gretlopt.OPT_C, &mp_cols, 1, 12, 0)
];

enum n_mp_options = 5;   /* the number of entries of mp_options */

int set_multiplot_sizes(Gretlopt opt)
{
    MplotOption * mpo;
    int i, k, err = 0;

    for (i=0; i<n_mp_options && !err; i++) {
    mpo = &mp_options[i];
    if (opt & mpo.flag) {
        k = get_optval_int(GretlCmdIndex.gridplot, mpo.flag, &err);
        if (!err && (k < mpo.min || k > mpo.max)) {
        gretl_errmsg_set("gridplot: out-of-bounds option value");
        err = GretlError.invarg;
        }
        if (!err) {
        *mpo.target = k;
        }
    }
    }

    return err;
}

void set_multiplot_defaults()
{
    int i;

    for (i=0; i<n_mp_options; i++) {
    *(mp_options[i].target) = mp_options[i].def;
    }
}

/* Set up the strings array in which gpbuild will cumulate
   sub-plots, with checks that @param is not the name of
   an existing object of incompatible type, and is a
   valid gretl identifier.
*/

int initialize_mp_array(const(char) * param, DatasetRaw * dset)
{
    int err = 0;

    mp_array = get_strings_array_by_name(param);

    if (mp_array != null) {
    ; /* OK, there's already a strings array of this name */
    } else if (gretl_is_user_var(param) ||
           current_series_index(dset, param) >= 0) {
    /* there's already an object of this name, but of the wrong type */
    err = GretlError.types;
    } else {
    /* check validity of @param as identifier */
    err = check_identifier(param);
    }

    if (!err) {
    mp_array = gretl_array_new(GretlType.strings, 0, &err);
    }

    if (mp_array != null) {
    array_name = g_strdup(param);
    mp_collecting = 1;
    }

    return err;
}

void gretl_gridplot_clear(int err)
{
    if (mp_array != null) {
    if (err) {
        gretl_array_destroy(mp_array);
    } else {
        err = user_var_add_or_replace(array_name,
                      GretlType.strings,
                      mp_array);
    }
        mp_array = null;
    }

    g_free(array_name);
    array_name = null;
    mp_collecting = 0;
}

/* called from interact.c: process_command_error() */

void gretl_gridplot_destroy()
{
    gretl_gridplot_clear(1);
}

/* This responds to the starting command for a "gpbuild" block */

int gretl_gridplot_start(const(char) * param, Gretlopt opt, DatasetRaw * dset)
{
    int err = 0;

    if (mp_array == null) {
    err = initialize_mp_array(param, dset);
    } else {
        gretl_errmsg_set(_("gpbuild: cannot be nested").ptr);
        err = GretlError.data;
    }

    return err;
}

int invalid_mp_error(int ci)
{
    gretl_errmsg_sprintf(_("%s: invalid (multiplot) plot specification").ptr,
             gretl_command_word(ci));
    return GretlError.invarg;
}

/* Append a plot specification to the @multiplot array. This is
   called in graphing.c, but only inside a gpbuild block.
*/

int gretl_gridplot_add_plot(char * buf)
{
    int err = 0;

    if (mp_array != null && buf != null) {
        if (strstr(buf, "set multiplot")) {
            err = invalid_mp_error(GretlCmdIndex.gpbuild);
        } else {
        gretl_array_append_string(mp_array, buf, 1);
        }
    } else {
    gretl_errmsg_set("gretl_gridplot_add_plot: failed");
        err = GretlError.data;
    }

    return err;
}

int set_mp_grid(int n_plots)
{
    int err = 0;

    if (mp_rows > n_plots) {
    gretl_errmsg_sprintf("invalid --rows specification for %d plots".ptr, n_plots);
    err = GretlError.invarg;
    } else if (mp_cols > n_plots) {
    gretl_errmsg_sprintf("invalid --cols specification for %d plots".ptr, n_plots);
    err = GretlError.invarg;
    }

    if (!err) {
    if (mp_rows > 0) {
        /* automatic cols value */
        mp_cols = cast(int) ceil(n_plots / cast(double) mp_rows);
    } else if (mp_cols > 0) {
        /* automatic rows value */
        mp_rows = cast(int) ceil(n_plots / cast(double) mp_cols);
    } else {
        /* fully automatic */
        mp_rows = cast(int) ceil(sqrt(cast(double) n_plots));
        mp_cols = cast(int) ceil(n_plots / cast(double) mp_rows);
    }
    }

    return err;
}

/* Write subplot buffer @buf to file, stripping out any
   "set term..." statements.
*/

void filter_subplot(const(char) * buf, FILE * fp)
{
    size_t sz = 2048;
    char * line = cast(char *) calloc(sz, 1);

    bufgets_init(buf);
    while (safe_bufgets(&line, &sz, buf)) {
    if (strncmp(line, "set term", 8)) {
        fputs(line, fp);
    }
    }
    bufgets_finalize(buf);
    free(line);
}

int get_subplot_index(GretlMatrixRaw * m, int i, int j)
{
    return cast(int) gretl_matrix_get(m, i, j) - 1;
}

/* Write a multiplot specification to file, drawing on the
   strings array @a.

   @np indicates the number of included plots and @maxp the
   total number of subplots available (the length of the
   relevant array).
*/

int output_multiplot_script(GretlArray * a, GretlMatrixRaw * m, int np, int maxp, const(char) * title)
{
    const(char) * buf;
    int i, j, k, p;
    int err = 0;
    FILE * fp;

    fp = open_plot_input_file(PlotType.gridplot, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    fputs("# literal lines = 1\n", fp);
    fprintf(fp, "# grid_params: plots=%d, fontsize=%d, width=%d, height=%d, ".ptr,
        np, mp_fontsize, mp_width, mp_height);
    fprintf(fp, "rows=%d, cols=%d\n".ptr, mp_rows, mp_cols);
    fprintf(fp, "set multiplot layout %d,%d rowsfirst".ptr, mp_rows, mp_cols);
    if (title != null) {
    fprintf(fp, " title '%s'\n".ptr, title);
    } else {
    fputc('\n', fp);
    }
    gretl_push_c_numeric_locale();

    k = -1;
    p = 0;
    for (i=0; i<mp_rows; i++) {
    for (j=0; j<mp_cols; j++) {
        if (m != null) {
        k = get_subplot_index(m, i, j);
        } else {
        k++;
        }
        if (k < 0) {
        fputs("set multiplot next\n", fp);
        } else {
        buf = null;
        if (k < maxp) {
            buf = cast(char *) gretl_array_get_data(a, k);
        }
        if (buf != null) {
            if (i + j > 0) {
            fputs("reset\n", fp);
            }
            fprintf(fp, "# subplot %d\n".ptr, ++p);
            if (strstr(buf, "set term")) {
            filter_subplot(buf, fp);
            } else {
            fputs(buf, fp);
            }
        }
        }
    }
    }

    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);
    err = finalize_plot_input_file(fp);

    return err;
}

/* set multiplot layout using a matrix argument */

int set_mp_layout(GretlMatrixRaw ** pm, int * np)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gridplot, Gretlopt.OPT_L);
    GretlMatrixRaw * m = null;
    int err = 0;

    if (s != null) {
    m = get_matrix_by_name(s);
    }
    if (gretl_is_null_matrix(m)) {
    err = GretlError.invarg;
    } else {
    int i, n = m.rows * m.cols;
        int nonzero = 0;
    double mi;

    for (i=0; i<n; i++) {
        mi = m.val[i];
        if (na(mi) || mi != floor(mi) ||
        mi < 0 || mi > *np) {
        gretl_errmsg_set(_("Invalid layout specification").ptr);
        err = GretlError.invarg;
        break;
        } else if (mi != 0) {
                nonzero++;
            }
    }
    if (!err) {
        *pm = m;
        *np = nonzero;
        mp_rows = m.rows;
        mp_cols = m.cols;
    }
    }

    return err;
}

/* respond to "end gpbuild" */

int gretl_gridplot_finalize(Gretlopt opt)
{
    int err = 0;

    if (mp_array == null) {
    gretl_errmsg_set("end gpbuild: building not started");
    err = GretlError.data;
    } else {
    gretl_gridplot_clear(0);
    }

    return err;
}

int is_graphic_filename(const(char) * s)
{
    static const(char) *[7] exts = [".png".ptr, ".pdf".ptr, ".svg".ptr, ".eps".ptr, ".fig".ptr, ".emf".ptr, ".html".ptr];
    int i;

    for (i=0; i < cast(int) exts.length; i++) {
        if (has_suffix(s, exts[i])) {
            return 1;
        }
    }

    return 0;
}

/* If we find that @s is not really a gnuplot commands buffer, try
   treating it as the name of file that may contain gnuplot
   commands. This requires ruling out the possibility that it's a
   graphic file (e.g. a PNG).  If it looks plausible, return the
   contents of the file.
*/

char * try_for_plot_filename(const(char) * s)
{
    char * ret = null;

    if (strlen(s) >= FILENAME_MAX) {
    return null;
    }

    if (!is_graphic_filename(s)) {
        char[FILENAME_MAX] fname;
        GretlFileType ft;

        strcpy(fname.ptr, s);
        ft = detect_filetype(fname.ptr, Gretlopt.OPT_P);
        if (ft == GretlFileType.csv) {
            /* generic text data: might be gnuplot commands */
            char * tmp = null;
            gsize sz = 0;

            g_file_get_contents(fname.ptr, &tmp, &sz, null);
            if (tmp != null) {
                ret = gretl_strdup(tmp);
                g_free(tmp);
            }
        }
    }

    return ret;
}

/* Find and check the strings array identified by @argname.
   Called when gridplot is used independently of gpbuild.
*/

int retrieve_plots_array(const(char) * argname, GretlArray ** pa, int * pnp)
{
    GretlArray * a = null;
    int msg_set = 0;
    int err = 0;

    a = get_array_by_name(argname);
    if (a == null) {
    err = GretlError.invarg;
    } else {
    int i, n = gretl_array_get_length(a);
        char * content = null;
    char * buf;

    for (i=0; i<n; i++) {
        buf = cast(char *) gretl_array_get_data(a, i);
            if (buf != null && strchr(buf, '\n') == null) {
                /* @buf can't be an actual plot buffer */
                content = try_for_plot_filename(buf);
                if (content != null) {
                    free(buf);
                    gretl_array_set_data(a, i, content);
                }
                buf = content;
            }
            /* check for embedded multiplots */
            if (buf == null || strstr(buf, "set multiplot")) {
        err = invalid_mp_error(GretlCmdIndex.gridplot);
        msg_set = 1;
        break;
        }
    }
    if (!err) {
        *pa = a;
        *pnp = n;
    }
    }

    if (err && !msg_set) {
    gretl_errmsg_set("Didn't get a valid array of plot strings");
    }

    return err;
}

/* This supports production of a multiplot from the array of
   plot-specification strings identified by @param. Such an
   array may be produced by "gpbuild", or may be assembled
   manually by the user.
*/

int gretl_gridplot_from_array(const(char) * param, Gretlopt opt)
{
    GretlArray * a = null;
    GretlMatrixRaw * m = null;
    const(char) * title = null;
    int maxp, np = 0;
    int err = 0;

    if (mp_collecting) {
        gretl_errmsg_set("gridplot: a block is in progress");
        return GretlError.data;
    }

    err = retrieve_plots_array(param, &a, &np);
    if (err) {
        return err;
    }

    maxp = np;
    set_multiplot_defaults();

    if (opt) {
    err = set_multiplot_sizes(opt);
    if (!err) {
        if (opt & Gretlopt.OPT_L) {
        err = set_mp_layout(&m, &np);
        } else {
        err = set_mp_grid(np);
        }
        if (opt & Gretlopt.OPT_T) {
        title = get_optval_string(GretlCmdIndex.gridplot, Gretlopt.OPT_T);
        }
    }
    } else {
    /* no options supplied, figure the default grid */
    err = set_mp_grid(np);
    }

    if (!err) {
    set_special_plot_size(mp_width, mp_height);
    set_special_font_size(mp_fontsize);
        err = output_multiplot_script(a, m, np, maxp, title);
    }

    clear_stored_options_for_command(GretlCmdIndex.gridplot);

    return err;
}

int check_gridplot_options(Gretlopt opt)
{
    int err;

    /* can't have both --output and --outbuf */
    err = incompatible_options(opt, Gretlopt.OPT_U | Gretlopt.OPT_B);
    if (!err) {
    /* can't have more than one of --rows, --cols, --layout */
    err = incompatible_options(opt, Gretlopt.OPT_R | Gretlopt.OPT_C | Gretlopt.OPT_L);
    }

    return err;
}

