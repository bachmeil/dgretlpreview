/* Port of lib/src/graphing.c: construction and launching of gnuplot graphs: plot-file creation and the gnuplot
   command line, time-series, scatter, factorized, frequency, correlogram, forecast, impulse-response, boxplot
   and other special plots, the plot-type and terminal bookkeeping and the user-visible "gnuplot" command.
   The private header is plot_privh.d and the public header graphingh.d.

   Differences from the C: the Windows branches and the GUI-only code are not ported, and HAVE_SYS_WAIT_H is
   taken as true. */

module gretl.graphing;


import core.stdc.time : tm, strftime;
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
import gretl.gretl_intl : get_local_decpoint, gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.gstrfuncs : g_strdup;
import gretl.plot_privh : BPMode, GPNA, GnuplotInfo;
import gretl.mapinfoh : MapFlags, MapInfo, MapProj, NaAction;
version (Posix) import core.sys.posix.sys.wait : WEXITSTATUS, WIFEXITED, waitpid; else import gretl.winport : WEXITSTATUS, WIFEXITED, waitpid;
version (Posix) import core.sys.posix.unistd : close, read, write; else import gretl.winport : close, read, write;
enum PI = 3.14159265358979323846;
import core.stdc.config : c_ulong;
extern (C) c_ulong strtoul(const(char) * nptr, char ** endptr, int base);
import gretl.dataio : date_maj_min, ntolabel;
import gretl.dataset : create_auxiliary_dataset, dataset_allocate_obs_markers, dataset_period_label, destroy_dataset, get_panel_group_name, get_panel_group_table, panel_group_names_ok, panel_group_names_varname, panel_sample_size, series_get_display_name, series_get_label, series_get_string_table, series_record_display_name, series_set_display_name;
import gretl.dataseth : calendar_data, dataset_is_panel, dataset_is_time_series, dated_daily_data, dated_weekly_data, quarterly_or_monthly, sample_size;
import gretl.describe : gretl_max, gretl_mean, gretl_minmax, gretl_sorted_series, gretl_stddev;
import gretl.describeh : FreqDist;
import gretl.estimate : lsq;
import gretl.forecasth : FitResid;
import gretl.genmain : current_series_index, generate_int;
import gretl.gerror : GError, g_error_free;
import gretl.gfileutils : g_build_filename, g_file_get_contents;
import gretl.glib : gsize, guint32, guint8;
import gretl.gmem : g_free;
import gretl.graphingh : FitType, GP_HEIGHT, GP_LB_HEIGHT, GP_LB_WIDTH, GP_SQ_SIZE, GP_WIDTH, GP_XL_HEIGHT, GP_XL_WIDTH, GP_XW_WIDTH, GP_XXL_HEIGHT, GP_XXL_WIDTH, GptFlags, N_GP_LINETYPES, PlotType, TermType, frequency_plot_code, gretlRGB, maybe_big_multiplot;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_error_clear, invalid_field_error;
import gretl.gretl_list : gretl_consecutive_list_new, gretl_list_append_term, gretl_list_copy, gretl_list_delete_at_pos, gretl_list_has_separator, gretl_list_new, gretl_list_separator_position, gretl_list_split_on_separator, in_gretl_list;
import gretl.gretl_listh : LISTSEP;
import gretl.gretl_model : arma_spectrum_plot_data, clear_model, gretl_model_get_data, gretl_model_get_double;
import gretl.gretl_panel : plausible_panel_time_var, series_is_group_invariant, time_series_from_panel;
import gretl.gretl_paths : gretl_dotdir, gretl_fopen, gretl_gnuplot_path, gretl_home, gretl_make_dotpath, gretl_maybe_prepend_dir, gretl_mktemp, gretl_plotfile, gretl_png_font, gretl_remove, gretl_set_path_by_name;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_detach_stream, gretl_print_get_buffer, gretl_print_new, gretl_print_new_with_stream, pprintf, pputs;
import gretl.gretl_prnh : PrnType;
import gretl.gretl_string_table : SeriesTable, series_table_get_string;
import gretl.gretl_utils : accept_as_discrete, gretl_copy_file, gretl_isdummy, gretl_spawn, ijton, set_obs;
import gretl.gspawn : GSpawnFlags, g_spawn_async_with_pipes;
import gretl.gstrfuncs : g_locale_to_utf8, g_strchomp, g_strdup_printf, g_strfreev, g_strsplit, g_strstrip, g_utf8_validate;
import gretl.libset : get_tex_plot_opts, gretl_in_batch_mode, gretl_in_gui_mode, gretl_no_plots;
import gretl.matrix_extra : gretl_dataset_from_matrix, gretl_plotfit_matrices;
import gretl.missing : list_adjust_sample;
import gretl.options : get_optval_double, get_optval_string, incompatible_options, set_optval_string, transcribe_option_flags;
import gretl.plotbands : plot_with_band;
import gretl.plotspec : print_auto_fit_string, print_keypos_string, set_plotfit_line;
import gretl.plotspech : GpKeySpec, GptLine, GptSpec;
import gretl.plugins : get_plugin_function;
import gretl.printout : bufgets, bufgets_finalize, bufgets_init;
import gretl.pvalues : chisq_cdf_comp, normal_critval, normal_pvalue_2, snedecor_cdf_comp, student_cdf_inverse, student_critval, student_pvalue_2;
import gretl.pvaluesh : DistCode;
import gretl.strutils : gretl_charsub, gretl_scan_varname, gretl_strdup, gretl_utf8_truncate, has_suffix, integer_string, string_is_blank;
import gretl.strutilsh : SLASH;
import gretl.system : system_get_depvar;
import gretl.systemh : EquationSystem;
import gretl.var : VECM_get_EC_matrix, gretl_VAR_get_FEVD_matrix, gretl_VAR_get_impulse_response, gretl_VAR_get_residual_matrix, gretl_VAR_get_roots, gretl_VAR_get_variable_number;
import gretl.gretl_gridplot : gretl_gridplot_add_plot, gretl_gridplot_collecting;
import gretl.gretl_plot : plot_command_active;
import gretl.gretl_array : get_strings_array_by_name, gretl_array_destroy, gretl_array_get_data, gretl_array_new, gretl_array_set_element, gretl_array_set_string, is_strings_array_element;
import gretl.gretl_bundle : gretl_bundle_get_int, gretl_bundle_get_scalar, gretl_bundle_get_string, gretl_bundle_has_key, gretl_bundle_set_matrix, gretl_bundle_set_string;
import gretl.uservar : get_list_by_name, get_string_by_name, user_var_add_or_replace;
import gretl.boxplots : boxplots;
import gretl.genfuncs : gretl_plotx;
import gretl.gretl_color : colormix_array;
import gretl.nonparam : loess_fit, sort_pairs_by_x;
import gretl.interact_extra : check_stringvar_name, set_plot_produced;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
    void manufacture_gui_callback(int ci);
    const(char) * midas_pdstr(const(DatasetRaw) * dset, int cfac);
    int print_map_palette(MapInfo * mi, double gpver, FILE * fp);
    void set_map_plot_limits(MapInfo * mi, double * xlim, double * ylim, double margin);
}


/* length of buffer for "set term ..." */
enum TERMLEN = 256;

__gshared char[MAXLEN] gnuplot_path = 0;
__gshared int gp_small_font_size = 0;
__gshared double default_png_scale = 1.0;
__gshared int xwide = 0;

__gshared char[64] ad_hoc_font = 0;
__gshared char[32] plot_buffer_name = 0;
__gshared char[8] plot_buffer_idx = 0;
__gshared double ad_hoc_alpha = 0;

enum WithStyle {
    points,
    lines,
    impulses,
    lp,
    boxes,
    steps
}

enum MAX_LETTERBOX_LINES = 8;

bool ts_plot(const(GptSpec) * g)
{
    return (g.flags & GptFlags.ts) != 0;
}

struct PlotTypeInfo {
    PlotType ptype;
    const(char) * pstr;
}

__gshared PlotTypeInfo[] ptinfo = [

    PlotTypeInfo(PlotType.regular, null),
    PlotTypeInfo(PlotType.correlogram, "correlogram".ptr),
    PlotTypeInfo(PlotType.cusum, "CUSUM test".ptr),
    PlotTypeInfo(PlotType.forecast, "forecast plot".ptr),
    PlotTypeInfo(PlotType.freq_simple, "frequency plot (simple)".ptr),
    PlotTypeInfo(PlotType.freq_normal, "frequency plot (against normal)".ptr),
    PlotTypeInfo(PlotType.freq_gamma, "frequency plot (against gamma)".ptr),
    PlotTypeInfo(PlotType.freq_discrete, "frequency plot (discrete)".ptr),
    { PlotType.garch,          "GARCH residual plot" },
    { PlotType.hurst,          "rescaled range plot" },
    PlotTypeInfo(PlotType.irfboot, "impulse response plot with quantiles".ptr),
    { PlotType.kernel,         "kernel density plot" },
    { PlotType.leverage,       "leverage/influence plot" },
    PlotTypeInfo(PlotType.multi_basic, "multiple small plots".ptr),
    PlotTypeInfo(PlotType.periodogram, "periodogram".ptr),
    { PlotType.range_mean, "range-mean plot" },
    PlotTypeInfo(PlotType.h_test, "sampling distribution".ptr),
    PlotTypeInfo(PlotType.prob_dist, "probability distribution".ptr),
    PlotTypeInfo(PlotType.tri_graph, "TRAMO / X12A tri-graph".ptr),
    PlotTypeInfo(PlotType.roots, "roots plot".ptr),
    PlotTypeInfo(PlotType.ellipse, "confidence ellipse plot".ptr),
    PlotTypeInfo(PlotType.multi_irf, "multiple impulse responses".ptr),
    { PlotType.panel,          "multiple panel plots" },
    PlotTypeInfo(PlotType.bi_graph, "double time-series plot".ptr),
    PlotTypeInfo(PlotType.many_ts, "multiple timeseries".ptr),
    PlotTypeInfo(PlotType.rq_tau, "tau sequence plot".ptr),
    PlotTypeInfo(PlotType.factorized, "factorized scatter".ptr),
    PlotTypeInfo(PlotType.boxplots, "boxplots".ptr),
    PlotTypeInfo(PlotType.curve, "curve".ptr),
    PlotTypeInfo(PlotType.qq, "QQ plot".ptr),
    { PlotType.user,           "user-defined plot" },
    PlotTypeInfo(PlotType.xcorrelogram, "cross-correlogram".ptr),
    PlotTypeInfo(PlotType.bar, "bars".ptr),
    PlotTypeInfo(PlotType.stacked_bar, "stacked-bars".ptr),
    PlotTypeInfo(PlotType.plot_3d, "3-D plot".ptr),
    PlotTypeInfo(PlotType.band, "band plot".ptr),
    PlotTypeInfo(PlotType.heatmap, "heatmap".ptr),
    PlotTypeInfo(PlotType.geomap, "geoplot".ptr),
    { PlotType.gridplot,       "user-multi" },
    PlotTypeInfo(PlotType.type_max, null)
];

/**
 * gnuplot_test_command:
 * @cmd: gnuplot command string.
 *
 * See if the installed version of gnuplot will accept a given
 * command.
 *
 * Returns: 0 if gnuplot successfully handles the given command,
 * 1 on error.
 */

int gnuplot_test_command(const(char) * cmd)
{
    int ok, ret = 1;
    int child_pid = 0, sinp = 0, serr = 0;
    GError *error = null;
    char *[2] argv = [null, null];

    if (*gnuplot_path.ptr == '\0') {
        strcpy(gnuplot_path.ptr, gretl_gnuplot_path());
    }

    argv[0] = gnuplot_path.ptr;

    ok = g_spawn_async_with_pipes (null,
                                   argv.ptr,
                                   null,
                                   GSpawnFlags.search_path |
                                   GSpawnFlags.stdout_to_dev_null |
                                   GSpawnFlags.do_not_reap_child,
                                   null,
                                   null,
                                   &child_pid,
                                   &sinp,
                                   null,
                                   &serr,
                                   &error);

    if (ok) {
        char[128] errbuf;
        int test, status;
        int errbytes = 0;

        errbytes += write(sinp, cmd, strlen(cmd));
        errbytes += write(sinp, "\n".ptr, 1);
        close(sinp);
        test = waitpid(child_pid, &status, 0);
        if (test == child_pid && WIFEXITED(status)) {
            ret = WEXITSTATUS(status);
        }
        errbytes = cast(int) read(serr, errbuf.ptr, errbuf.sizeof - 1);
        if (errbytes > 0) {
            errbuf[errbytes] = '\0';
            if (strstr(errbuf.ptr, "not find/open font")) {
                if (strstr(cmd, "font") != null) {
                    ret = 1;
                }
            }
        }
        close(serr);
    } else {
        fprintf(stderr, "error: '%s'\n".ptr, error.message);
        g_error_free(error);
    }

    return ret;
}

/* Apparatus for getting the gnuplot version, either
   in numerical form or as a string
*/

__gshared char[16] gpver_string = 0;
__gshared double gpver_number = 0;

void set_gpver_number(const(char) * s)
{
    const(char) * fmt = "%d.%d.%d";
    int maj, min, plev;

    if (sscanf(s, fmt, &maj, &min, &plev) == 3) {
        gpver_number = maj + min / 10.0 + plev / 100.0;
    }
}

int get_gp_version_info()
{
    char * qname = null;
    char * fname = null;
    char * buf = null;
    FILE * fp;
    int err = 0;

    if (*gnuplot_path.ptr == '\0') {
        strcpy(gnuplot_path.ptr, gretl_gnuplot_path());
    }

    if (*gnuplot_path.ptr == '\0') {
        return 0;
    }

    qname = gretl_make_dotpath("gpver_query");
    fname = gretl_make_dotpath("gpver.txt");

    fp = gretl_fopen(qname, "w");
    if (fp == null) {
        err = GretlError.fopen;
    } else {
        /* write a little query script */
        fputs("gpver = sprintf(\"%.1f.%s\", GPVAL_VERSION, GPVAL_PATCHLEVEL)\n".ptr, fp);
        fprintf(fp, "set print \"%s\"\n".ptr, fname);
        fputs("print gpver\n", fp);
        fclose(fp);
    }

    if (!err) {
        /* send script to gnuplot for execution */
        buf = g_strdup_printf("\"%s\" \"%s\"".ptr, gnuplot_path.ptr, qname);
        err = gretl_spawn(buf);
    }

    if (!err) {
        /* read version info from response */
        fp = gretl_fopen(fname, "r");
        if (fp != null) {
            if (fgets(gpver_string.ptr, gpver_string.sizeof, fp)) {
                g_strstrip(gpver_string.ptr);
                set_gpver_number(gpver_string.ptr);
            }
            fclose(fp);
        }
    }

    gretl_remove(qname);
    gretl_remove(fname);
    g_free(qname);
    g_free(fname);
    g_free(buf);

    return err;
}

double gnuplot_version()
{
    if (gpver_number == 0) {
        get_gp_version_info();
    }

    return gpver_number;
}

const(char) * gnuplot_version_string()
{
    if (gpver_number == 0) {
        get_gp_version_info();
    }

    return gpver_number == 0 ? "unknown" : gpver_string.ptr;
}

/* end gnuplot version apparatus */

int gp_list_pos(const(char) * s, const(int) * list, const(DatasetRaw) * dset)
{
    int k;

    if (integer_string(s)) {
        k = atoi(s);
    } else {
        k = current_series_index(dset, s);
    }

    return in_gretl_list(list, k);
}

__gshared int plot_ci = GretlCmdIndex.gnuplot;

void set_effective_plot_ci(int ci)
{
    plot_ci = ci;
}

void reset_effective_plot_ci()
{
    plot_ci = GretlCmdIndex.gnuplot;
}

int get_effective_plot_ci()
{
    return plot_ci;
}

/* When we get from the user something like

   --with-lines=foo,bar

   or

   --with-lines=<listname>

   this indicates that the "with lines" format should be
   applied to selected y-axis variables, not all.
*/

int gp_set_non_point_info(GnuplotInfo * gi, const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    const(char) * s = get_optval_string(cast(GretlCmdIndex) plot_ci, opt);
    int withval = WithStyle.points;
    int i, imax = gi.withlist[0];
    int * slist = null;

    if (opt == Gretlopt.OPT_O) {
        withval = WithStyle.lines;
    } else if (opt == Gretlopt.OPT_M) {
        withval = WithStyle.impulses;
    } else if (opt == Gretlopt.OPT_P) {
        withval = WithStyle.lp;
    } else if (opt == Gretlopt.OPT_B) {
        withval = WithStyle.boxes;
    } else if (opt & Gretlopt.OPT_Q) {
        withval = WithStyle.steps;
    }

    if (s == null) {
        /* spec applies to all members of list */
        for (i=1; i<=imax; i++) {
            if (gi.withlist[i] == WithStyle.points) {
                gi.withlist[i] = withval;
            }
        }
    } else if (strchr(s, ',') != null) {
        /* spec has multiple components */
        char ** strs = g_strsplit(s, ",", 0);
        int j;

        for (j=0; strs[j]!=null; j++) {
            i = gp_list_pos(strs[j], list, dset);
            if (i > 0 && i <= imax) {
                gi.withlist[i] = withval;
            }
        }
        g_strfreev(strs);
    } else if ((slist = get_list_by_name(s)) != null) {
        /* got a named list spec */
        int j;

        for (j=1; j<=slist[0]; j++) {
            i = in_gretl_list(list, slist[j]);
            if (i > 0 && i <= imax) {
                gi.withlist[i] = withval;
            }
        }
    } else {
        /* just one component */
        i = gp_list_pos(s, list, dset);
        if (i > 0 && i <= imax) {
            gi.withlist[i] = withval;
        }
    }

    return 0;
}

/* note: @plot_ci is a file-scope global */

int plain_lines_spec(Gretlopt opt)
{
    if ((opt & Gretlopt.OPT_O) && !(opt & (Gretlopt.OPT_M | Gretlopt.OPT_B | Gretlopt.OPT_P | Gretlopt.OPT_Q))) {
        return get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_O) == null;
    } else {
        return 0;
    }
}

int plain_impulses_spec(Gretlopt opt)
{
    if ((opt & Gretlopt.OPT_M) && !(opt & (Gretlopt.OPT_O | Gretlopt.OPT_B | Gretlopt.OPT_P | Gretlopt.OPT_Q))) {
        return get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_M) == null;
    } else {
        return 0;
    }
}

int plain_steps_spec(Gretlopt opt)
{
    if ((opt & Gretlopt.OPT_Q) && !(opt & (Gretlopt.OPT_O | Gretlopt.OPT_M | Gretlopt.OPT_B | Gretlopt.OPT_P))) {
        return get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_Q) == null;
    } else {
        return 0;
    }
}

int get_fit_type(GnuplotInfo * gi)
{
    const(char) * ftype = get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_F);
    int err = 0;

    if (ftype == null || *ftype == '\0') {
        err = GretlError.data;
    } else if (!strcmp(ftype, "none")) {
        gi.flags |= GptFlags.fit_omit;
    } else if (!strcmp(ftype, "linear")) {
        gi.fit = FitType.ols;
    } else if (!strcmp(ftype, "quadratic")) {
        gi.fit = FitType.quadratic;
    } else if (!strcmp(ftype, "cubic")) {
        gi.fit = FitType.cubic;
    } else if (!strcmp(ftype, "inverse")) {
        gi.fit = FitType.inverse;
    } else if (!strcmp(ftype, "loess")) {
        gi.fit = FitType.loess;
    } else if (!strcmp(ftype, "semilog")) {
        gi.fit = FitType.loglin;
    } else if (!strcmp(ftype, "linlog")) {
        gi.fit = FitType.linlog;
    } else {
        err = invalid_field_error(ftype);
    }

    return err;
}

void maybe_record_font_choice(Gretlopt opt)
{
    const(char) * s = get_optval_string(cast(GretlCmdIndex) plot_ci, opt);

    if (s != null) {
        ad_hoc_font[0] = '\0';
        strcat(ad_hoc_font.ptr, s);
        gretl_charsub(ad_hoc_font.ptr, ',', ' ');
    }
}

void maybe_record_alpha(Gretlopt opt)
{
    double alpha;
    int err = 0;

    alpha = get_optval_double(cast(GretlCmdIndex) plot_ci, opt, &err);
    if (!na(alpha) && alpha > 0 && alpha <= 1) {
        ad_hoc_alpha = alpha;
    } else {
        ad_hoc_alpha = 0;
    }
}

/* We come here in response to the --buffer option in plot
   commands, and if applicable we set the static variables
   @plot_buffer_name and @plot_buffer_idx.
*/

int set_plot_buffer_name(const(char) * bname)
{
    int err = 0;

    if (bname == null || *bname == '\0') {
        *plot_buffer_name.ptr = '\0';
        *plot_buffer_idx.ptr = '\0';
    } else if (is_strings_array_element(bname, plot_buffer_name.ptr, plot_buffer_idx.ptr)) {
        ; /* handled */
    } else {
        *plot_buffer_idx.ptr = '\0';
        err = check_stringvar_name(bname, 1, null);
        if (!err) {
            strcpy(plot_buffer_name.ptr, bname);
        }
    }

    return err;
}

int plot_output_to_buffer()
{
    return *plot_buffer_name.ptr != '\0' ||
        gretl_gridplot_collecting();
}

int make_plot_commands_buffer(const(char) * fname)
{
    char * contents = null;
    GError *gerr = null;
    gsize size = 0;
    int err = 0;

    g_file_get_contents(fname, &contents, &size, &gerr);

    if (gerr != null) {
        gretl_errmsg_set(gerr.message);
        g_error_free(gerr);
        err = GretlError.fopen;
    } else if (gretl_gridplot_collecting()) {
        err = gretl_gridplot_add_plot(contents);
    } else if (*plot_buffer_idx.ptr != '\0') {
        GretlArray * a = get_strings_array_by_name(plot_buffer_name.ptr);
        int i = generate_int(plot_buffer_idx.ptr, null, &err);

        if (a != null && !err) {
            err = gretl_array_set_string(a, i-1, contents, 1);
        }
    } else {
        char * buf = gretl_strdup(contents);

        err = user_var_add_or_replace(plot_buffer_name.ptr,
                                      GretlType.string_,
                                      buf);
    }

    g_free(contents);
    gretl_remove(fname);

    return err;
}

int get_gp_flags(GnuplotInfo * gi, Gretlopt opt, const(int) * list, const(DatasetRaw) * dset)
{
    int n_yvars = list[0] - 1;
    int err = 0;

    gi.flags = cast(GptFlags) 0;

    /* can't have --single-yaxis and --y2axis */
    err = incompatible_options(opt, Gretlopt.OPT_Y | Gretlopt.OPT_D);

    if (opt & (Gretlopt.OPT_N | Gretlopt.OPT_A)) {
        /* --band or --bands */
        if (opt & Gretlopt.OPT_T) {
            /* --time-series */
            gi.flags |= (GptFlags.ts | GptFlags.idx);
            /* there's no xvar in @list */
            n_yvars++;
        }
        gi.flags |= GptFlags.fit_omit;
        gi.band = 1;
        goto linespec;
    }

    if (opt & Gretlopt.OPT_W) {
        /* --font=<fontspec> */
        maybe_record_font_choice(Gretlopt.OPT_W);
    }

     if (opt & Gretlopt.OPT_E) {
        /* --alpha=<val> */
        maybe_record_alpha(Gretlopt.OPT_E);
    }

    if (opt & Gretlopt.OPT_L) {
        /* log y axis */
        const(char) * sbase = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_L);

        gi.flags |= GptFlags.logy;
        gi.ybase = 10;
        if (sbase != null) {
            gi.ybase = atof(sbase);
            if (gi.ybase <= 0) {
                gi.ybase = 10;
            }
        }
    }

    if (opt & Gretlopt.OPT_S) {
        /* the old --suppress-fitted option may still be used
           internally, for some plot types */
        gi.flags |= GptFlags.fit_omit;
    }

    if (opt & Gretlopt.OPT_R) {
        /* internal option for residual plot */
        gi.flags |= GptFlags.resids;
    } else if (opt & Gretlopt.OPT_A) {
        /* internal option for fitted-actual plot */
        gi.flags |= GptFlags.fa;
    }

    if (opt & Gretlopt.OPT_Z) {
        /* --factorized */
        gi.flags |= (GptFlags.factor | GptFlags.fit_omit);
    } else if (opt & Gretlopt.OPT_C) {
        /* --control */
        gi.flags |= GptFlags.xyz;
    } else {
        if (opt & Gretlopt.OPT_T) {
            /* --time-series */
            gi.flags |= GptFlags.idx;
            /* there's no xvar in @list */
            n_yvars++;
        }
    }

 linespec:

    if (plain_lines_spec(opt)) {
        /* just using lines */
        gi.flags |= GptFlags.lines;
    } else if (plain_impulses_spec(opt)) {
        /* just using impulses */
        gi.flags |= GptFlags.impulses;
    } else if (plain_steps_spec(opt)) {
        /* just using steps */
        gi.flags |= GptFlags.steps;
    } else if (opt & (Gretlopt.OPT_M | Gretlopt.OPT_O | Gretlopt.OPT_P | Gretlopt.OPT_B)) {
        /* for handling per-variable "plot with" options */
        gi.withlist = gretl_list_new(n_yvars);
    }

    if (gi.withlist != null) {
        if (opt & Gretlopt.OPT_M) {
            /* --with-impulses */
            gp_set_non_point_info(gi, list, dset, Gretlopt.OPT_M);
        }
        if (opt & Gretlopt.OPT_O) {
            /* --with-lines */
            gp_set_non_point_info(gi, list, dset, Gretlopt.OPT_O);
        }
        if (opt & Gretlopt.OPT_P) {
            /* --with-lp */
            gp_set_non_point_info(gi, list, dset, Gretlopt.OPT_P);
        }
        if (opt & Gretlopt.OPT_B) {
            /* --with-boxes */
            gp_set_non_point_info(gi, list, dset, Gretlopt.OPT_B);
        }
    }

    if (opt & Gretlopt.OPT_G) {
    /* internal option, saving as icon */
    gi.flags |= GptFlags.icon;
    }

    gi.fit = FitType.none;

    if (!(gi.flags & GptFlags.fit_omit) && n_yvars == 1) {
        if (opt & Gretlopt.OPT_F) {
            /* the --fit=fitspec option */
            err = get_fit_type(gi);
        }
    }

    if (xwide) {
        /* access file-scope global */
        gi.flags |= GptFlags.xw;
        xwide = 0;
    }

    return err;
}

void write_gp_dataval(double x, FILE * fp, int final_)
{
    if (final_) {
        if (na(x)) {
            fprintf(fp, "%s\n".ptr, GPNA.ptr);
        } else {
            fprintf(fp, "%.10g\n".ptr, x);
        }
    } else {
        if (na(x)) {
            fprintf(fp, "%s ".ptr, GPNA.ptr);
        } else {
            fprintf(fp, "%.10g ".ptr, x);
        }
    }
}

void printvars(FILE * fp, int t, const(int) * list, const(DatasetRaw) * dset, GnuplotInfo * gi, const(char) * label, double offset)
{
    const(double) * x = (gi != null)? gi.x : null;
    double xt;
    int i;

    if (x != null) {
        xt = x[t] + offset;
        if (gi.flags & GptFlags.timefmt) {
            fprintf(fp, "%.0f ".ptr, xt);
        } else {
            fprintf(fp, "%.10g ".ptr, xt);
        }
    }

    for (i=1; i<=list[0]; i++) {
        xt = dset.Z[list[i]][t];
        if (!na(xt) && x == null && i == 1) {
            /* the x variable */
            xt += offset;
        }
        write_gp_dataval(xt, fp, 0);
    }

    if (label != null) {
        fprintf(fp, "# %s".ptr, label);
    }

    fputc('\n', fp);
}

int y_var_is_factorized(GnuplotInfo * gi, int vj)
{
    if (gi.flist == null) {
        return 1;
    } else if (in_gretl_list(gi.flist, vj)) {
        return 1;
    } else {
        return 0;
    }
}

/* Here we see if the user added a list or series parameter to the
   --factorized plot option, to limit factorization to a subset of the
   y-variables. If so, we check the specified series: it/they must be
   present in the "y" portion of the plot list or else we flag an
   error. If there's no error we attach this information to @gi under
   the name "flist".
*/

int maybe_add_factorization_list(GnuplotInfo * gi, const(DatasetRaw) * dset)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_Z);
    int err = 0;

    if (s != null) {
        int lmax = gi.list[gi.list[0]];
        int * flist = null;
        int vf = 0;
        int i, k;

        if ((flist = get_list_by_name(s)) != null) {
            for (i=1; i<=flist[0]; i++) {
                k = in_gretl_list(gi.list, flist[i]);
                if (k == 0 || k == lmax) {
                    err = GretlError.invarg;
                    break;
                }
            }
            if (!err) {
                gi.flist = gretl_list_copy(flist);
            }
        } else if (strchr(s, ',') != null) {
            /* spec has multiple comma-separated components */
            char ** strs = g_strsplit(s, ",", 0);
            int ns = 0;

            for (i=0; strs[i] != null; i++) {
                vf = current_series_index(dset, strs[i]);
                k = in_gretl_list(gi.list, vf);
                if (k == 0 || k == lmax) {
                    err = GretlError.invarg;
                    break;
                } else {
                    ns++;
                }
            }
            if (!err) {
                gi.flist = gretl_list_new(ns);
                for (i=0; strs[i]!=null; i++) {
                    vf = current_series_index(dset, strs[i]);
                    gi.flist[i+1] = vf;
                }
            }
            g_strfreev(strs);
        } else if ((vf = current_series_index(dset, s)) > 0) {
            k = in_gretl_list(gi.list, vf);
            if (k == 0 || k == lmax) {
                err = GretlError.invarg;
            } else {
                gi.flist = gretl_list_new(1);
                gi.flist[1] = vf;
            }
        } else {
            err = GretlError.invarg;
        }
    }

    return err;
}

int factor_check(GnuplotInfo * gi, const(DatasetRaw) * dset)
{
    int err = 0;

    if (gi.list[0] < 2) {
        err = GretlError.data;
    } else if (!accept_as_discrete(dset, gi.fid, 0)) {
        err = GretlError.data;
    }

    if (err) {
        gretl_errmsg_set(_("You must supply at least three variables, the last of which is discrete").ptr);
    } else {
        const(double) * d = dset.Z[gi.fid] + gi.t1;
        int T = gi.t2 - gi.t1 + 1;

        gi.fvals = gretl_matrix_values(d, T, Gretlopt.OPT_S, &err);
        if (!err) {
            gi.n_fvals = gretl_vector_get_length(gi.fvals);
            err = maybe_add_factorization_list(gi, dset);
        }
    }

    return err;
}

int gnuplot_has_wxt()
{
    __gshared int err = -1;

    if (err == -1) {
        err = gnuplot_test_command("set term wxt");
    }

    return !err;
}

int gnuplot_has_x11()
{
    __gshared int err = -1;

    if (err == -1) {
        err = gnuplot_test_command("set term x11");
    }

    return !err;
}

int gnuplot_has_qt()
{
    __gshared int err = -1;

    if (err == -1) {
        err = gnuplot_test_command("set term qt");
    }

    return !err;
}

int gnuplot_has_tikz()
{
    __gshared int err = -1;

    if (err == -1) {
        err = gnuplot_test_command("set term tikz");
    }

    return !err;
}

__gshared gretlRGB[N_GP_LINETYPES] user_color = [0xff0000, 0x0000ff, 0x00cc00, 0xbf25b2, 0x8faab3, 0xffa500, 0xe51e10, 0x000000];

/* By default, returns an array of strings holding the string-valued
   specifications of the gnuplot line colors currently selected by the
   user.  However, if @f is non-missing, non-zero and in (-1.0, 1.0)
   each color returned is shifted lighter (> 0, towards white) or darker
   (< 0, towards black) relative to the original.
*/

GretlArray * get_user_colors(double f, int * err)
{
    GretlArray * ret;

    ret = gretl_array_new(GretlType.strings, N_GP_LINETYPES, err);

    if (ret != null) {
    guint32 c0, c1 = 0;
    char[10] s;
    int i;

    if (f != 0) {
        if (f < 0.0) {
        /* move towards white */
        c1 = 0xFFFFFF;
        }
        f = fabs(f);
    }
    for (i=0; i<N_GP_LINETYPES; i++) {
        c0 = user_color[i];
        if (f != 0.0) {
        GretlArray * mix = colormix_array(c0, c1, &f, 1, 0, err);

        sprintf(s.ptr, "%s".ptr, cast(char *) gretl_array_get_data(mix, 0));
        gretl_array_destroy(mix);
        } else {
        sprintf(s.ptr, "#%06x".ptr, c0);
        }
        gretl_array_set_element(ret, i, s.ptr, GretlType.string_, 1);
    }
    }

    return ret;
}

/* apparatus for handling plot "extra" colors */

__gshared gretlRGB[3] extra_color = [0x5f6b84, 0xdddddd, 0x66aaaaaa];

__gshared gretlRGB[3] user_extra_color = [0x5f6b84, 0xdddddd, 0x66aaaaaa];

gretlRGB get_boxcolor()
{
    return user_extra_color[0];
}

void set_boxcolor (gretlRGB color) {
    user_extra_color[0] = color;
}

gretlRGB get_shadecolor()
{
    return user_extra_color[1];
}

gretlRGB get_alpha_shadecolor()
{
    return user_extra_color[2];
}

void set_shadecolor (gretlRGB color) {
    user_extra_color[1] = color;
}

void print_rgb_hash(char * s, gretlRGB color)
{
    if (color > 0xffffff) {
        sprintf(s, "#%08X".ptr, color);
    } else {
        sprintf(s, "#%06X".ptr, color);
    }
}

gretlRGB gretl_rgb_get(const(char) * s)
{
    gretlRGB x = 0;

    if (sscanf(s, "#%x".ptr, &x) != 1) {
        x = 0;
    }

    return x;
}

/* called only by print_palette_string() below */

void print_argb_hex(char * s, gretlRGB color)
{
    if (color > 0xffffff) {
        sprintf(s, "x%08X".ptr, color);
    } else {
        sprintf(s, "x%06X".ptr, color);
    }
}

/* We use this format only when writing to the gretl rc file */

void print_palette_string(char * s)
{
    char[10][2] cstr;

    print_argb_hex(cstr[0].ptr, user_extra_color[0]);
    print_argb_hex(cstr[1].ptr, user_extra_color[1]);
    sprintf(s, "%s %s".ptr, cstr[0].ptr, cstr[1].ptr);
}

gretlRGB get_graph_color(int i)
{
    return (i >= 0 && i < N_GP_LINETYPES)? user_color[i] : 0;
}

void set_graph_color_from_string(int i, const(char) * s)
{
    int err = 0;

    if (i >= 0 && i < 2) {
        gretlRGB u;

        if (sscanf(s + 1, "%x".ptr, &u) == 1) {
            user_extra_color[i] = u;
        } else {
            err = 1;
        }
    } else {
        err = 1;
    }

    if (err) {
        fprintf(stderr, "Error in set_graph_palette_from_string(%d, '%s')\n".ptr,
                i, s);
    }
}

void graph_palette_reset(int i)
{
    if (i >= 0 && i < 2) {
        user_extra_color[i] = extra_color[i];
    }
}

/* Given a string @s such as "Sans 8" or "Bodoni MT 12", write
   the name part into @name and the point-size part into @psz.
   Return 2 if we got both a name and a size, 1 if we just
   got a name, 0 if we got nothing.
*/

int split_graph_fontspec(const(char) * s, char * name, int * psz)
{
    int i, k = 0;
    int n, nf = 0;

    *name = '\0';
    s += cast(int) strspn(s, " ");

    if (isdigit(*s)) {
        /* we only got a size? */
        *psz = atoi(s);
        strcpy(name, "sans");
        return 2;
    }

    /* look for digits at the end */
    n = cast(int) (strlen(s));
    for (i=n-1; i>0; i--) {
        if (isdigit(s[i])) k++;
        else break;
    }

    if (k > 0) {
        /* got a size */
        char[8] ptstr;

        *ptstr.ptr = *name = '\0';
        strncat(ptstr.ptr, s + n - k, k);
        *psz = atoi(ptstr.ptr);
        strncat(name, s, n - k - 1);
        nf = 2;
    } else if (*s != '\0') {
        nf = 1;
        strcpy(name, s);
    }

    return nf;
}

char * gretl_png_font_string()
{
    const(char) * s = gretl_png_font();
    char[256] fstr;
    char[128] name;
    int nf, ptsize = 0;

    fstr[0] = '\0';
    nf = split_graph_fontspec(s, name.ptr, &ptsize);
    if (nf == 2) {
        sprintf(fstr.ptr, " font \"%s,%d\"".ptr, name.ptr, ptsize);
    } else if (nf == 1) {
        sprintf(fstr.ptr, " font \"%s\"".ptr, name.ptr);
    }

    return gretl_strdup(fstr.ptr);
}

void maybe_set_small_font(int nplots)
{
    gp_small_font_size = (nplots > 4)? 6 : 0;
}

/* apparatus for plots at custom sizes (e.g. maps) */

__gshared float special_width = 0;
__gshared float special_height = 0;
__gshared int special_fontsize = 0;

void set_special_plot_size(float width, float height)
{
    special_width = width;
    special_height = height;
}

void set_special_font_size(int fsize)
{
    special_fontsize = fsize;
}

void clear_special_plot_size()
{
    special_width = special_height = 0;
}

void clear_special_font_size()
{
    special_fontsize = 0;
}

int special_plot_size_is_set()
{
    return special_width > 0 && special_height > 0;
}

int special_font_size_is_set()
{
    return special_fontsize > 0;
}

/* end special size apparatus */

void write_png_font_string(char * fstr, char * ad_hoc_fontspec, PlotType ptype, const(char) * specfont, double scale)
{
    const(char) * thisfont;

    if (specfont != null) {
        thisfont = specfont;
    } else if (*ad_hoc_font.ptr != '\0') {
        thisfont = ad_hoc_font.ptr;
    } else {
        thisfont = gretl_png_font();
    }

    if (*thisfont == '\0') {
        /* fallback */
        thisfont = getenv("GRETL_PNG_GRAPH_FONT");
    }

    if (*thisfont == '\0') {
        *fstr = '\0';
        return;
    } else {
        char[128] fname;
        int nf, fsize = 0;

        nf = split_graph_fontspec(thisfont, fname.ptr, &fsize);
        if (special_font_size_is_set()) {
            fsize = special_fontsize;
        } else if (nf == 2) {
            if (maybe_big_multiplot(ptype) && gp_small_font_size > 0) {
                fsize = gp_small_font_size;
            }
            if (scale > 1.0) {
                fsize = cast(int) round(scale * fsize);
            }
        }
        if (fsize > 0 && fsize < 100) {
            sprintf(fstr, " font \"%s,%d\"".ptr, fname.ptr, fsize);
        } else if (nf == 1) {
            sprintf(fstr, " font \"%s\"".ptr, fname.ptr);
        }
        if (thisfont == ad_hoc_font.ptr && ad_hoc_fontspec != null) {
            strcpy(ad_hoc_fontspec, thisfont);
        }
        /* ensure these settings don't outstay their welcome */
        ad_hoc_font[0] = '\0';
        clear_special_font_size();
    }
}

/* for gnuplot pdfcairo, epscairo output */

char * write_other_font_string(int stdsize)
{
    char * fstr = null;

    if (ad_hoc_font[0] != '\0') {
        char[128] fname;
        int nf, fsize = 0;

        nf = split_graph_fontspec(ad_hoc_font.ptr, fname.ptr, &fsize);
        if (nf == 2) {
            fstr = g_strdup_printf("%s,%d".ptr, fname.ptr, fsize);
        } else if (nf == 1) {
            fstr = g_strdup_printf("%s,%d".ptr, fname.ptr, stdsize);
        }
        ad_hoc_font[0] = '\0';
    } else if (special_font_size_is_set()) {
        fstr = g_strdup_printf("sans,%d".ptr, special_fontsize);
    } else {
        fstr = g_strdup_printf("sans,%d".ptr, stdsize);
    }

    return fstr;
}

/* gnuplot styles apparatus */

__gshared char[32] gp_style = 0; /* basename of style file */
__gshared char * alt_sty = null;    /* content of style file */

__gshared const(char) * classic_sty =
    "# gpstyle classic\nset linetype 1 pt 1 lc rgb \"#FF0000\"\nset linetype 2 pt 2 lc rgb \"#0000FF\"\nset linetype 3 pt 3 lc rgb \"#00CC00\"\nset linetype 4 pt 4 lc rgb \"#BF25B2\"\nset linetype 5 pt 5 lc rgb \"#8FAAB3\"\nset linetype 6 pt 6 lc rgb \"#FFA500\"\nset linetype 7 pt 7 lc rgb \"#E51E10\"\nset linetype 8 pt 8 lc rgb \"#000000\"\n"; /* black */

/* Read 8 line colors out of @s and transcribe to the
   @user_color array. Fail if we can't get all 8.
*/

int transcribe_style(const(char) * s)
{
    const(char) * p = s;
    gretlRGB rgb;
    int i, got;

    for (i=0; i<N_GP_LINETYPES; i++) {
        p = strstr(p, " lc rgb ");
        got = 0;
        if (p != null) {
            p += 8;
            p += cast(int) strspn(p, " ");
            if (sscanf(p+1, "#%x".ptr, &rgb) == 1) {
                user_color[i] = rgb;
                got = 1;
                p += 8;
            }
        }
        if (!got) {
            break;
        }
    }

    return i < N_GP_LINETYPES ? GretlError.data : 0;
}

/* Given a style-name (provisionally copied to @gp_style),
   try to find its file and check it for conformance.
*/

int try_set_alt_sty()
{
    GError *gerr = null;
    gsize len = 0;
    char * try_ = null;
    char * fname;
    int err = 0;

    fname = g_build_filename(gretl_home(), "data".ptr,
                             "gnuplot".ptr, gp_style.ptr, null);

    if (g_file_get_contents(fname, &try_, &len, &gerr)) {
        err = transcribe_style(try_);
        if (err) {
            /* failed to conform to spec */
            fprintf(stderr, "%s failed spec check\n".ptr, gp_style.ptr);
            set_plotstyle("classic");
        } else {
            /* OK, put the style in place */
            g_free(alt_sty);
            alt_sty = try_;
        }
    } else {
        /* couldn't find the file */
        if (gerr != null) {
            fprintf(stderr, "%s\n".ptr, gerr.message);
            g_error_free(gerr);
        }
        err = GretlError.fopen;
    }

    g_free(fname);

    return err;
}

/* callback from libset.c, in response to "set plot_style <name>" */

int set_plotstyle(const(char) * style)
{
    int to_classic = 0;

    if (!strcmp(style, "classic") || !strcmp(style, "default")) {
        /* "default" is just for backward compat */
        to_classic = 1;
    }

    if (!strcmp(style, gp_style.ptr)) {
        return 0; /* no-op */
    } else if (to_classic && *gp_style.ptr == '\0') {
        return 0; /* no-op */
    } else if (to_classic) {
        /* replace alt with classic */
        g_free(alt_sty);
        alt_sty = null;
        gp_style[0] = '\0';
        transcribe_style(classic_sty);
        return 0;
    } else {
        /* try replacing current with what's requested */
        sprintf(gp_style.ptr, "%s.gpsty".ptr, style);
        return try_set_alt_sty();
    }
}

/* called from libset.c */

const(char) * get_plotstyle()
{
    __gshared char[32] pstyle = 0;
    char * p;

    strcpy(pstyle.ptr, gp_style.ptr);
    p = strrchr(pstyle.ptr, '.');
    if (p != null) {
        *p = '\0';
    }

    if (pstyle[0] == '\0') {
        return "classic";
    } else {
        return pstyle.ptr;
    }
}

/* e.g. line = "set linetype 1 pt 1 lc rgb \"#FF0000\"\n" */

void put_line_with_added_alpha(const(char) * line, guint8 a, FILE * fp)
{
    char * s = g_strdup(line);
    char * p = strstr(s, "rgb ");
    guint32 u;

    *p = '\0';
    fputs(s, fp);
    sscanf(p+5, "#%x".ptr, &u);
    /* put 'a' = alpha into the high bits */
    fprintf(fp, "rgb \"#%08X\"\n".ptr, (a << 24) | u);
    g_free(s);
}

/* Write the content of either the default, or an alternative,
   gnuplot style into @fp. The @offset argument allows for
   skipping one or more leading linetype definitions. The
   @linetypes argument allows for skipping the linetypes
   altogether.
*/

void inject_gp_style(int offset, int linetypes, FILE * fp)
{
    const(char) * sty = alt_sty != null ? alt_sty : classic_sty;

    if (linetypes == 0) {
        /* just print what follows the linetypes, if anything */
        const(char) * sub = strstr(sty, "set bord");

        if (sub != null) {
            fputs(sub, fp);
        }
    } else {
        const(char) * src = sty;
        guint8 a = 0;

        if (!na(ad_hoc_alpha) && ad_hoc_alpha > 0) {
            a = cast(ubyte) nearbyint(ad_hoc_alpha * 255);
            ad_hoc_alpha = 0; /* don't persist */
        }

        if (offset > 0) {
            /* start at a specified linetype */
            char[32] targ;

            sprintf(targ.ptr, "set linetype %d".ptr, offset + 1);
            src = strstr(sty, targ.ptr);
            if (src == null) {
                src = sty;
            }
        }
        if (a == 0) {
            fputs(src, fp);
        } else {
            char[128] line;

            bufgets_init(src);
            while (bufgets(line.ptr, line.sizeof, src)) {
                if (strstr(line.ptr, "linetype") && strstr(line.ptr, "rgb ")) {
                    put_line_with_added_alpha(line.ptr, a, fp);
                } else {
                    fputs(line.ptr, fp);
                }
            }
            bufgets_finalize(src);
        }
    }
}

void write_plot_line_styles(int ptype, FILE * fp)
{
    char[12] cstr;
    int i;

    if (ptype == PlotType.plot_3d) {
        for (i=0; i<2; i++) {
            print_rgb_hash(cstr.ptr, user_color[i]);
            fprintf(fp, "set linetype %d lc rgb \"%s\"\n".ptr, i+1, cstr.ptr);
        }
    } else if (ptype == PlotType.boxplots) {
        for (i=0; i<2; i++) {
            print_rgb_hash(cstr.ptr, user_color[i+1]);
            fprintf(fp, "set linetype %d lc rgb \"%s\"\n".ptr, i+1, cstr.ptr);
        }
        inject_gp_style(0, 0, fp);
    } else if (frequency_plot_code(cast(PlotType) ptype)) {
        print_rgb_hash(cstr.ptr, get_boxcolor());
        fprintf(fp, "set linetype 1 lc rgb \"%s\"\n".ptr, cstr.ptr);
        fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
        inject_gp_style(0, 0, fp);
    } else if (ptype == PlotType.rq_tau) {
        fputs("set linetype 1 lc rgb \"#000000\"\n", fp);
        fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
        fputs("set linetype 3 lc rgb \"#0000FF\"\n", fp);
        fputs("set linetype 4 lc rgb \"#0000FF\"\n", fp);
        fputs("set linetype 5 lc rgb \"#0000FF\"\n", fp);
        fputs("set linetype 6 lc rgb \"#FFA500\"\n", fp);
        fputs("set linetype 7 lc rgb \"#E51E10\"\n", fp);
        fputs("set linetype 8 lc rgb \"#000000\"\n", fp);
        inject_gp_style(0, 0, fp);
    } else if (ptype == PlotType.heatmap || ptype == PlotType.geomap) {
        ; /* these are handled specially */
    } else {
        /* the primary default case */
        inject_gp_style(0, 1, fp);
    }
}

/* end gnuplot styles apparatus */

/* Get gnuplot to print the dimensions of a PNG plot, in terms
   of both pixels and data bounds (gnuplot >= 4.4.0).
*/

int write_plot_bounding_box_request(FILE * fp)
{
    char * fname = gretl_make_dotpath("gretltmp.png.bounds");

    fprintf(fp, "set print \"%s\"\n".ptr, fname);
    g_free(fname);

    fputs("print \"pixel_bounds: \", GPVAL_TERM_XMIN, GPVAL_TERM_XMAX, GPVAL_TERM_YMIN, GPVAL_TERM_YMAX\n", fp);
    fputs("print \"data_bounds: \", GPVAL_X_MIN, GPVAL_X_MAX, GPVAL_Y_MIN, GPVAL_Y_MAX\n", fp);
    fputs("print \"term_size: \", GPVAL_TERM_XSIZE, GPVAL_TERM_YSIZE, GPVAL_TERM_SCALE\n", fp);

    return 0;
}

int do_plot_bounding_box()
{
    FILE * fp = gretl_fopen(gretl_plotfile(), "ab");
    int err = 0;

    if (fp != null) {
        err = write_plot_bounding_box_request(fp);
        fclose(fp);
    } else {
        err = GretlError.fopen;
    }

    return err;
}

void maybe_set_eps_pdf_dims(char * s, PlotType ptype, GptFlags flags)
{
    double w = 0, h = 0;

    if (special_plot_size_is_set()) {
        w = (5.0 * special_width) / GP_WIDTH;
        h = (3.5 * special_height) / GP_HEIGHT;
        clear_special_plot_size();
    } else if (flags & GptFlags.letterbox) {
        /* for time series */
        w = (5.0 * GP_LB_WIDTH) / GP_WIDTH;
        h = (3.5 * GP_LB_HEIGHT) / GP_HEIGHT;
    } else if (flags & GptFlags.xl) {
        /* large */
        w = (5.0 * GP_XL_WIDTH) / GP_WIDTH;
        h = (3.5 * GP_XL_HEIGHT) / GP_HEIGHT;
    } else if (flags & GptFlags.xxl) {
        /* extra large */
        w = h = (5.0 * GP_XXL_WIDTH) / GP_WIDTH;
    } else if (flags & GptFlags.xw) {
        /* extra wide */
        w = (5.0 * GP_XW_WIDTH) / GP_WIDTH;
        h = 3.5;
    } else if (ptype == PlotType.roots || ptype == PlotType.qq) {
        /* square plots */
        w = h = 3.5;
    }

    if (w > 0 && h > 0) {
        char[32] size_str;

        gretl_push_c_numeric_locale();
        sprintf(size_str.ptr, " size %.2f,%.2f".ptr, w, h);
        gretl_pop_c_numeric_locale();
        strcat(s, size_str.ptr);
    }
}

void append_gp_encoding(char * s)
{
    strcat(s, "\nset encoding utf8");
}

/* In pdf and eps term lines: should "dashed" be appended when
   "mono" is specified? Try experimenting?
*/

char * eps_pdf_term_line(char * term_line, PlotType ptype, GptFlags flags, TermType ttype)
{
    char * font_string;
    const(char) * tname;
    int ptsize;

    ptsize = (ptype == PlotType.multi_basic)? 6 : 12;
    tname = (ttype == TermType.eps)? "epscairo" : "pdfcairo";

    font_string = write_other_font_string(ptsize);
    sprintf(term_line, "set term %s noenhanced font \"%s\"".ptr,
            tname, font_string);
    g_free(font_string);

    maybe_set_eps_pdf_dims(term_line, ptype, flags);
    append_gp_encoding(term_line);

    return term_line;
}

void maybe_append_tex_opts(char * s)
{
    const(char) * opts = get_tex_plot_opts();

    if (opts != null) {
        strcat(s, " ");
        strcat(s, opts);
    }
}

char * tex_term_line(char * term_line, PlotType ptype, GptFlags flags)
{
    if (gpver_number < 5.4) {
        strcpy(term_line, "set term cairolatex");
    } else {
        strcpy(term_line, "set term pict2e");
    }
    maybe_append_tex_opts(term_line);
    append_gp_encoding(term_line);

    return term_line;
}

char * tikz_term_line(char * term_line, PlotType ptype, GptFlags flags)
{
    strcpy(term_line, "set term tikz");
    maybe_append_tex_opts(term_line);
    append_gp_encoding(term_line);

    return term_line;
}

void plot_get_scaled_dimensions(int * width, int * height, double scale)
{
    *width = cast(int) (*width * scale);
    *height = cast(int) (*height * scale);

    /* PNG: round up to an even number of pixels if need be */
    if (*width % 2) *width += 1;
    if (*height % 2) *height += 1;
}

void write_png_size_string(char * s, PlotType ptype, GptFlags flags, double scale)
{
    int w = GP_WIDTH, h = GP_HEIGHT;

    if (special_plot_size_is_set()) {
        w = cast(int) special_width;
        h = cast(int) special_height;
        clear_special_plot_size();
    } else if (flags & GptFlags.letterbox) {
        /* time series plots */
        w = GP_LB_WIDTH;
        h = GP_LB_HEIGHT;
    } else if (flags & GptFlags.xl) {
        /* large */
        w = GP_XL_WIDTH;
        h = GP_XL_HEIGHT;
    } else if (flags & GptFlags.xxl) {
        /* extra large */
        w = GP_XXL_WIDTH;
        h = GP_XXL_HEIGHT;
    } else if (flags & GptFlags.xw) {
        /* extra wide */
        w = GP_XW_WIDTH;
        h = GP_HEIGHT;
    } else if (ptype == PlotType.roots || ptype == PlotType.qq) {
        /* square plots */
        w = h = GP_SQ_SIZE;
    }

    if (scale != 1.0 && ptype != PlotType.gridplot) {
        plot_get_scaled_dimensions(&w, &h, scale);
    }

    *s = '\0';
    sprintf(s, " size %d,%d".ptr, w, h);
}

/* platform-specific on-screen term string */

char * var_term_line(char * term_line, int ptype, GptFlags flags)
{
    char[140] font_string;
    char[16] size_string;
    const(char) * varterm = null;

    const(char) * pref = getenv("GNUTERM");
    const(char) *[] strs = ["wxt", "qt", "x11"];
    guint8[3] opts;
    int i;

    opts[0] = cast(ubyte) gnuplot_has_wxt();
    opts[1] = cast(ubyte) gnuplot_has_qt();
    opts[2] = cast(ubyte) gnuplot_has_x11();

    if (pref != null) {
        for (i=0; i<3; i++) {
            if (!strcmp(pref, strs[i]) && opts[i]) {
                varterm = strs[i];
            }
        }
    } else {
        for (i=0; i<3; i++) {
            if (opts[i]) {
                varterm = strs[i];
            }
        }
    }
    if (varterm == null) {
        varterm = "x11"; /* ?? */
    }

    *font_string.ptr = *size_string.ptr = '\0';
    write_png_font_string(font_string.ptr, null, cast(PlotType) ptype, null, 1.0);
    write_png_size_string(size_string.ptr, cast(PlotType) ptype, flags, 1.0);

    sprintf(term_line, "set term %s%s%s noenhanced".ptr,
            varterm, font_string.ptr, size_string.ptr);
    append_gp_encoding(term_line);

    return term_line;
}

char * real_png_term_line(char * term_line, PlotType ptype, GptFlags flags, const(char) * specfont, double scale)
{
    char[128] ad_hoc_fontspec;
    char[140] font_string;
    char[16] size_string;

    *font_string.ptr = *size_string.ptr = *ad_hoc_fontspec.ptr = '\0';

    /* Note: @specfont will be non-NULL only if we're coming
       from the GUI with a GptFlags.spec struct in hand.
    */

    write_png_font_string(font_string.ptr, ad_hoc_fontspec.ptr,
                          ptype, specfont, scale);
    write_png_size_string(size_string.ptr, cast(PlotType) ptype, flags, scale);

    sprintf(term_line, "set term pngcairo%s%s noenhanced".ptr,
            font_string.ptr, size_string.ptr);
    append_gp_encoding(term_line);

    if (*ad_hoc_fontspec.ptr != '\0') {
        strcat(term_line, "\n# fontspec: ");
        strcat(term_line, ad_hoc_fontspec.ptr);
    }

    return term_line;
}

char * gretl_png_term_line(char * term_line, PlotType ptype, GptFlags flags)
{
    if (ptype == PlotType.geomap) {
        return real_png_term_line(term_line, ptype, flags, null, 1.0);
    } else {
        double s = default_png_scale;

        return real_png_term_line(term_line, ptype, flags, null, s);
    }
}

/**
 * gretl_gnuplot_term_line:
 * @ttype: code for the gnuplot "terminal" type.
 * @ptype: indication of the sort of plot to be made, which
 * may make a difference to the color palette chosen.
 * @flags: plot option flags.
 * @font: if non-NULL, try to respect a specified font.
 *
 * Constructs a suitable line for sending to gnuplot to invoke
 * the specified "terminal".
 *
 * Returns: a static char * pointer.
 */

const(char) * gretl_gnuplot_term_line(TermType ttype, PlotType ptype, GptFlags flags, const(char) * font)
{
    __gshared char[TERMLEN] term_line = 0;

    *term_line.ptr = '\0';

    if (font != null && ad_hoc_font[0] == '\0') {
        strcpy(ad_hoc_font.ptr, font);
    }

    if (ttype == TermType.png) {
        double s = default_png_scale;

        real_png_term_line(term_line.ptr, ptype, flags, null, s);
    } else if (ttype == TermType.eps || ttype == TermType.pdf) {
        eps_pdf_term_line(term_line.ptr, cast(PlotType) ptype, flags, cast(TermType) ttype);
    } else if (ttype == TermType.emf) {
        gretl_emf_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.svg) {
        strcpy(term_line.ptr, "set term svg noenhanced");
        append_gp_encoding(term_line.ptr);
    } else if (ttype == TermType.htm) {
        strcpy(term_line.ptr, "set term canvas noenhanced");
        append_gp_encoding(term_line.ptr);
    } else if (ttype == TermType.fig) {
        strcpy(term_line.ptr, "set term fig");
        append_gp_encoding(term_line.ptr);
    } else if (ttype == TermType.tex) {
        tex_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.tkz) {
        tikz_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.var) {
        var_term_line(term_line.ptr, ptype, flags);
    }

    return term_line.ptr;
}

const(char) * get_png_line_for_plotspec(GptSpec * spec)
{
    __gshared char[TERMLEN] term_line = 0;

    real_png_term_line(term_line.ptr, spec.code, spec.flags,
                       spec.fontstr, spec.scale);
    return term_line.ptr;
}

void set_default_png_scale(double s)
{
    if (s >= 0.5 && s <= 2.0) {
        default_png_scale = s;
    }
}

double get_default_png_scale()
{
    return default_png_scale;
}

void write_emf_font_string(char * fstr)
{
    const(char) * src = null;
    int stdsize = 16;
    int adhoc = 0;

    if (ad_hoc_font[0] != '\0') {
        src = ad_hoc_font.ptr;
        adhoc = 1;
    } else {
        src = gretl_png_font();
    }

    if (src != null) {
        char[128] fname;
        int nf, fsize = 0;

        nf = split_graph_fontspec(src, fname.ptr, &fsize);
        if (nf == 2) {
            if (adhoc) {
                /* go with what the user specified */
                sprintf(fstr, "font \"%s,%d\"".ptr, fname.ptr, fsize);
            } else {
                /* adjust size to avoid tiny text? */
                fsize = (fsize <= 8)? 12 : 16;
                sprintf(fstr, "font \"%s,%d\"".ptr, fname.ptr, fsize);
            }
        } else if (nf == 1) {
            sprintf(fstr, "font \"%s,%d\"".ptr, fname.ptr, stdsize);
        }
        ad_hoc_font[0] = '\0';
    } else {
        sprintf(fstr, "font \"sans,%d\"".ptr, stdsize);
    }
}

char * gretl_emf_term_line(char * term_line, PlotType ptype, GptFlags flags)
{
    char * size_string = null;
    char[140] font_string;

    *font_string.ptr = '\0';
    write_emf_font_string(font_string.ptr);

    if (special_plot_size_is_set()) {
        size_string = g_strdup_printf("size %d,%d ".ptr,
                                      cast(int) special_width,
                                      cast(int) special_height);
        clear_special_plot_size();
    }

    if (flags & GptFlags.mono) {
        strcat(term_line, "set term emf dash noenhanced ");
    } else {
        strcat(term_line, "set term emf color noenhanced ");
    }

    if (size_string != null) {
        strcat(term_line, size_string);
        g_free(size_string);
    }

    if (*font_string.ptr != '\0') {
        strcat(term_line, font_string.ptr);
    }

    append_gp_encoding(term_line);

    return term_line;
}

/**
 * plot_type_from_string:
 * @str: initial comment line from plot file.
 *
 * Returns: the special plot code corresponding to the initial
 * comment string in the plot file, or %PlotType.regular if no special
 * comment is recognized.
 */

PlotType plot_type_from_string(const(char) * str)
{
    int i, len;
    PlotType ret = PlotType.regular;

    for (i=1; ptinfo[i].pstr != null; i++) {
        len = cast(int) (strlen(ptinfo[i].pstr));
        if (!strncmp(str + 2, ptinfo[i].pstr, len)) {
            ret = ptinfo[i].ptype;
            break;
        }
    }

    return ret;
}

int write_plot_type_string(PlotType ptype, GptFlags flags, FILE * fp)
{
    int i, ret = 0;

    if (ptype == PlotType.geomap) {
        /* handled specially */
        return 0;
    }

    for (i=1; i<PlotType.type_max; i++) {
        if (ptype == ptinfo[i].ptype) {
            if (flags & GptFlags.xl) {
                fprintf(fp, "# %s (large)\n".ptr, ptinfo[i].pstr);
            } else if (flags & GptFlags.xxl) {
                fprintf(fp, "# %s (extra-large)\n".ptr, ptinfo[i].pstr);
            } else if (flags & GptFlags.xw) {
                fprintf(fp, "# %s (extra-wide)\n".ptr, ptinfo[i].pstr);
            } else {
                fprintf(fp, "# %s\n".ptr, ptinfo[i].pstr);
            }
            ret = 1;
            break;
        }
    }

    if (ret == 0 && (flags & GptFlags.xw)) {
        fputs("# extra-wide\n", fp);
    }

    if (get_local_decpoint() == ',') {
        /* is this right? */
        fputs("set decimalsign ','\n", fp);
    }

    return ret;
}

void print_term_string(int ttype, PlotType ptype, GptFlags flags, FILE * fp)
{
    char[TERMLEN] term_line;

    *term_line.ptr = '\0';

    if (ttype == TermType.eps || ttype == TermType.pdf) {
        eps_pdf_term_line(term_line.ptr, cast(PlotType) ptype, flags, cast(TermType) ttype);
    } else if (ttype == TermType.png) {
        gretl_png_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.emf) {
        gretl_emf_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.fig) {
        strcpy(term_line.ptr, "set term fig\nset encoding utf8");
    } else if (ttype == TermType.svg) {
        strcpy(term_line.ptr, "set term svg noenhanced\nset encoding utf8");
    } else if (ttype == TermType.htm) {
        strcpy(term_line.ptr, "set term canvas noenhanced\nset encoding utf8");
    } else if (ttype == TermType.tex) {
        tex_term_line(term_line.ptr, ptype, flags);
    } else if (ttype == TermType.tkz) {
        tikz_term_line(term_line.ptr, ptype, flags);
    }

    if (*term_line.ptr != '\0') {
        fprintf(fp, "%s\n".ptr, term_line.ptr);
        if (flags & GptFlags.mono) {
            fputs("set mono\n", fp);
        } else if (ptype != PlotType.gridplot) {
            write_plot_line_styles(ptype, fp);
        }
    }
}

__gshared int gretl_plot_count = 0;
__gshared int this_term_type = 0;

const(char) * plotname(const(DatasetRaw) * dset, int v, int full)
{
    __gshared char[32] texname = 0;
    const(char) * ret = dset.varname[v];

    if (full) {
        const(char) * alt = series_get_display_name(dset, v);

        if (alt != null && *alt != '\0') {
            ret = alt;
        }
    }

    if ((this_term_type == TermType.tex || this_term_type == TermType.tkz) &&
         strchr(ret, '_') != null) {
        int i = 0;

        while (*ret) {
            if (*ret == '_') {
                texname[i++] = '\\';
                texname[i++] = '_';
            } else {
                texname[i++] = *ret;
            }
            ret++;
        }
        texname[i] = '\0';
        return texname.ptr;
    }

    return ret;
}

/* The trick used in plotname() -- with a static char object to hold a
   TeX-corrected name -- won't work if the caller wants to define two
   or more plotname strings before making use of them, since the first
   will be modified before it's used. So we have to make a copy, as
   below.
*/

char * make_plotname(const(DatasetRaw) * dset, int v)
{
    const(char) * s = plotname(dset, v, 1);

    return gretl_strdup(s);
}

/* recorder for filename given via --output=foo */
__gshared char[FILENAME_MAX] gnuplot_outname = 0;

int set_term_type_from_fname(const(char) * fname, int * err)
{
    if (has_suffix(fname, ".eps")) {
        this_term_type = TermType.eps;
    } else if (has_suffix(fname, ".ps")) {
        this_term_type = TermType.eps;
    } else if (has_suffix(fname, ".pdf")) {
        this_term_type = TermType.pdf;
    } else if (has_suffix(fname, ".png")) {
        this_term_type = TermType.png;
    } else if (has_suffix(fname, ".fig")) {
        this_term_type = TermType.fig;
    } else if (has_suffix(fname, ".emf")) {
        this_term_type = TermType.emf;
    } else if (has_suffix(fname, ".svg")) {
        this_term_type = TermType.svg;
    } else if (has_suffix(fname, ".html") ||
               has_suffix(fname, ".htm")) {
        this_term_type = TermType.htm;
    } else if (has_suffix(fname, ".tex")) {
        this_term_type = TermType.tex;
    } else if (has_suffix(fname, ".tikz")) {
        if (gnuplot_has_tikz()) {
            this_term_type = TermType.tkz;
        } else {
            gretl_errmsg_set("gnuplot does not support lua/tikz");
            this_term_type = TermType.none;
            *err = 1;
        }
    } else if (!strcmp(fname, "gnuplot")) {
        this_term_type = TermType.var;
    }

    return this_term_type;
}

int specified_gp_output_format()
{
    return this_term_type;
}

void reset_plot_count()
{
    gretl_plot_count = 0;
}

/* We're writing the 'set output...' line for a gnuplot script.
   If @path is non-NULL we use it, otherwise we make a path
   using dotdir and "gretltmp.png".
*/

int write_plot_output_line(const(char) * path, FILE * fp)
{
    const(char) * fname;
    char * tmp = null;

    if (path == null) {
        fname = tmp = gretl_make_dotpath("gretltmp.png");
    } else {
        fname = path;
    }

    fprintf(fp, "set output \"%s\"\n".ptr, fname);
    g_free(tmp);

    return 0;
}

FILE * gp_set_up_batch(char * fname, PlotType ptype, GptFlags flags, const(char) * outspec, int * err)
{
    int fmt = TermType.none;
    FILE * fp = null;

    if (plot_output_to_buffer()) {
        this_term_type = TermType.plt;
        sprintf(fname, "%sgpttmp.XXXXXX".ptr, gretl_dotdir());
        fp = gretl_mktemp(fname, "w");
    } else if (outspec != null) {
        /* user gave --output=<filename> */
        fmt = set_term_type_from_fname(outspec, err);
        if (*err) {
            return null;
        } else if (fmt) {
            /* input needs processing */
            strcpy(gnuplot_outname.ptr, outspec);
            gretl_maybe_prepend_dir(gnuplot_outname.ptr);
            sprintf(fname, "%sgpttmp.XXXXXX".ptr, gretl_dotdir());
            fp = gretl_mktemp(fname, "w");
        } else {
            /* just passing gnuplot commands through */
            this_term_type = TermType.plt;
            strcpy(fname, outspec);
            gretl_maybe_prepend_dir(fname);
            fp = gretl_fopen(fname, "w");
        }
    } else {
        /* auto-constructed gnuplot commands filename */
        this_term_type = TermType.plt;
        sprintf(fname, "gpttmp%02d.plt".ptr, ++gretl_plot_count);
        gretl_maybe_prepend_dir(fname);
        fp = gretl_fopen(fname, "w");
    }

    if (fp == null) {
        *err = GretlError.fopen;
    } else {
        gretl_set_path_by_name("plotfile", fname);
        if (*gnuplot_outname.ptr != '\0') {
            /* write terminal/style/output lines */
            print_term_string(fmt, ptype, flags, fp);
            write_plot_output_line(gnuplot_outname.ptr, fp);
        } else if (ptype != PlotType.gridplot) {
            /* just write style lines */
            write_plot_line_styles(ptype, fp);
        }
        if (get_local_decpoint() == ',') {
            fputs("set decimalsign ','\n", fp);
        }
    }

    return fp;
}

__gshared char * iact_gpfile = null;

/* Set-up for an "interactive" plot: we open a file in the user's
   dotdir into which gnuplot commands will be written.  If we're
   running the GUI program this command file will eventually be used
   to create a PNG file for display in a gretl window; otherwise
   (gretlcli) the commands will eventually be sent to gnuplot for
   "direct" display (e.g. using the wxt or windows terminal).

   In this function we just open the file for writing; if in GUI
   mode insert a suitable PNG terminal line and output spec line;
   and write some header-type material including our line style
   specifications.
*/

FILE * gp_set_up_interactive(char * fname, PlotType ptype, GptFlags flags, int * err)
{
    int gui = gretl_in_gui_mode();
    FILE * fp = null;

    if (iact_gpfile != null) {
        fname = iact_gpfile;
        fp = gretl_fopen(fname, "wb");
        iact_gpfile = null;
    } else if (gui) {
        /* the filename should be unique */
        sprintf(fname, "%sgpttmp.XXXXXX".ptr, gretl_dotdir());
        fp = gretl_mktemp(fname, "wb");
    } else {
        /* gretlcli: no need for uniqueness */
        sprintf(fname, "%sgpttmp.plt".ptr, gretl_dotdir());
        fp = gretl_fopen(fname, "wb");
    }

    if (fp == null) {
        *err = GretlError.fopen;
    } else {
        gretl_set_path_by_name("plotfile", fname);
        if (gui) {
            /* set up for PNG output */
            fprintf(fp, "%s\n".ptr, gretl_gnuplot_term_line(TermType.png, ptype,
                                                        flags, null));
            if (default_png_scale != 1.0) {
                gretl_push_c_numeric_locale();
                fprintf(fp, "# scale = %.1f\n".ptr, default_png_scale);
                fputs("# auto linewidth\n", fp);
                gretl_pop_c_numeric_locale();
            }
            write_plot_output_line(null, fp);
        } else if (ptype == PlotType.geomap) {
            fprintf(fp, "%s\n".ptr, gretl_gnuplot_term_line(TermType.var, ptype,
                                                        flags, null));
        } else {
            fprintf(fp, "%s\n".ptr, gretl_gnuplot_term_line(TermType.var, ptype,
                                                        flags, null));
        }
        write_plot_type_string(ptype, flags, fp);
        if (ptype != PlotType.gridplot) {
            write_plot_line_styles(ptype, fp);
        }
    }

    return fp;
}

int gnuplot_too_old()
{
    __gshared double gpv = 0;
    int msg_done = 0;
    int ret = 0;

    if (gpv == 0.0) {
        gpv = gnuplot_version();
    }

    if (gpv == 0.0) {
        if (!msg_done) {
            if (gretl_gnuplot_path()[0] == '\0') {
                gretl_errmsg_set("Path to gnuplot is not available");
            } else {
                gretl_errmsg_sprintf("'%s': gnuplot is broken or too old: must be >= version 5.0".ptr,
                                     gretl_gnuplot_path());
            }
        }
        ret = 1;
    } else if (gpv < 5.0) {
        gretl_errmsg_sprintf("Gnuplot %g is too old: must be >= version 5.0".ptr, gpv);
        ret = 1;
    }

    return ret;
}

int got_display_option(const(char) * s)
{
    return s != null && !strcmp(s, "display");
}

int got_none_option(const(char) * s)
{
    return s != null && !strcmp(s, "none");
}

int gpbuild_output_conflict(const(char) * s)
{
    const(char) * msg = "incompatible output specification".ptr;

    if (s == null) {
        gretl_errmsg_sprintf("gpbuild: %s".ptr, _(msg));
    } else {
        gretl_errmsg_sprintf("gpbuild: %s '%s'".ptr, _(msg), s);
    }

    return GretlError.badopt;
}

const(char) * plot_output_option(PlotType p, int * pci, int * err)
{
    int grid_mode = gretl_gridplot_collecting();
    int ci = get_effective_plot_ci();
    const(char) * s;

    /* clear any previous setting */
    set_plot_buffer_name(null);

    /* set a more specific command index based on
       the plot type, if applicable */

    if (p == PlotType.multi_basic) {
        ci = GretlCmdIndex.scatters;
    } else if (p == PlotType.boxplots) {
        ci = GretlCmdIndex.bxplot;
    } else if (p == PlotType.forecast) {
        ci = GretlCmdIndex.fcast;
    } else if (p == PlotType.correlogram) {
        ci = GretlCmdIndex.corrgm;
    } else if (p == PlotType.xcorrelogram) {
        ci = GretlCmdIndex.xcorrgm;
    } else if (p == PlotType.periodogram) {
        ci = GretlCmdIndex.pergm;
    } else if (p == PlotType.hurst) {
        ci = GretlCmdIndex.hurst;
    } else if (p == PlotType.qq) {
        ci = GretlCmdIndex.qqplot;
    } else if (p == PlotType.range_mean) {
        ci = GretlCmdIndex.rmplot;
    } else if (p == PlotType.leverage) {
        ci = GretlCmdIndex.leverage;
    } else if (p == PlotType.freq_simple ||
               p == PlotType.freq_normal ||
               p == PlotType.freq_gamma ||
               p == PlotType.freq_discrete) {
        ci = GretlCmdIndex.freq;
    } else if (p == PlotType.heatmap) {
        ci = GretlCmdIndex.corr;
    } else if (p == PlotType.cusum) {
        ci = GretlCmdIndex.cusum;
    } else if (p == PlotType.kernel) {
        ci = GretlCmdIndex.kdplot;
    } else if (p == PlotType.gridplot) {
        ci = GretlCmdIndex.gridplot;
    }

    /* --output or --plot, depending on the command */
    if (plot_command_active()) {
        s = get_optval_string(GretlCmdIndex.plot, Gretlopt.OPT_U);
    } else {
        s = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_U);
    }

    if (grid_mode && s != null) {
        if (!strcmp(s, "display")) {
            /* let this pass: hansl functions that do plots
               generally seem to default to "display"
            */
            s = null;
    } else if (!strcmp(s, "none")) {
        /* we should exclude this plot from the
           gridplot-collecting routine
        */
        return s;
        } else if (p == PlotType.geomap) {
            s = null;
        } else {
            *err = gpbuild_output_conflict(s);
            return null;
        }
    }

    if (s != null && *s == '\0') {
        /* cancel an empty @s */
        s = null;
    } else if (s == null) {
        /* try for --outbuf=<strname> */
        if (plot_command_active()) {
            s = get_optval_string(GretlCmdIndex.plot, Gretlopt.OPT_B);
        } else {
            s = get_optval_string(cast(GretlCmdIndex) ci, Gretlopt.OPT_B);
        }
        if (s != null && grid_mode) {
            char * tmp = g_strdup_printf("--outbuf=%s".ptr, s);

            *err = gpbuild_output_conflict(tmp);
            g_free(tmp);
            return null;
        } else if (s != null && *s == '\0') {
            s = null;
        } else {
            *err = set_plot_buffer_name(s);
        }
    }

    if (pci != null) {
        *pci = ci;
    }

    return s;
}

/**
 * open_plot_input_file:
 * @ptype: indication of the sort of plot to be made.
 * @flags: may inflect some characteristics of plot.
 * @err: location to receive error code.
 *
 * Opens a file into which gnuplot commands will be written.
 * Depending on the prospective use of the stream, some
 * header-type material may be written into it (the primary
 * case being when we're going to produce PNG output
 * for display in the gretl GUI). The prospective use is
 * figured out based on the program state, @ptype and
 * @flags.
 *
 * Returns: writable stream on success, %NULL on failure.
 */

FILE * open_plot_input_file(PlotType ptype, GptFlags flags, int * err)
{
    char[FILENAME_MAX] fname = 0;
    const(char) * outspec = null;
    int ci, interactive = 0;
    FILE * fp = null;

    /* ensure we have 'gnuplot_path' in place (file-scope static var) */
    if (*gnuplot_path.ptr == '\0') {
        strcpy(gnuplot_path.ptr, gretl_gnuplot_path());
    }

    if (gnuplot_too_old()) {
        *err = GretlError.external;
        return null;
    }

    /* initialize */
    this_term_type = TermType.none;
    *gnuplot_outname.ptr = '\0';

    /* check for --output=whatever or --outbuf=whatever */
    outspec = plot_output_option(ptype, &ci, err);

    if (!*err) {
    if (gretl_gridplot_collecting()) {
        interactive = 0;
    } else if (got_display_option(outspec)) {
        /* --output=display specified */
        interactive = 1;
    } else if (outspec != null) {
        /* --output=filename or --buffer=strvar specified */
        interactive = 0;
    } else if (flags & GptFlags.icon) {
        interactive = 1;
    } else {
        /* default */
        interactive = !gretl_in_batch_mode();
    }
    if (interactive) {
        fp = gp_set_up_interactive(fname.ptr, ptype, flags, err);
    } else {
        fp = gp_set_up_batch(fname.ptr, ptype, flags, outspec, err);
    }
    }

    return fp;
}

/* For plot-types that are generated by commands that also produce
   textual output: based on the program state and command option,
   figure out if a plot is actually wanted or not. In interactive mode
   the answer will be Yes unless --plot=none has been given; in batch
   mode the answer will be No unless --plot=display or --plot=fname
   has been given, or we're inside a "gpbuild" block.

   The @err location will receive an error code only in the gpbuild
   case, if the command in question a --plot specification that's
   incompatible with gpbuild. But err = NULL is accepted, in which
   case the error code is (obviously) not passed on, but the "wanted"
   return value is zero. This is a dispensation for commands where the
   plot is incidental.
*/

int gnuplot_graph_wanted(PlotType ptype, Gretlopt opt, int * err)
{
    const(char) * optname = null;
    int ret = 0;

    if (gretl_no_plots()) {
        return 0;
    }

    if (opt & (Gretlopt.OPT_U | Gretlopt.OPT_B)) {
        /* check for --plot or --outbuf option */
        int griderr = 0;

        optname = plot_output_option(ptype, null, &griderr);
        if (griderr) {
            if (err != null) {
                *err = griderr;
            }
            return 0;
        }
    }

    if (got_none_option(optname)) {
        /* --plot=none specified */
        ret = 0;
    } else if (optname != null) {
        /* --plot=display or --plot=fname specified */
        ret = 1;
    } else if (gretl_gridplot_collecting()) {
        ret = 1;
    } else {
        /* defaults */
        ret = !gretl_in_batch_mode();
    }

    return ret;
}

/**
 * gnuplot_cleanup:
 *
 * Removes any temporary gnuplot input file written in
 * the user's dot directory.
 */

void gnuplot_cleanup()
{
    const(char) * p, fname = gretl_plotfile();

    p = strstr(fname, "gpttmp");

    if (p != null) {
        int pnum;

        if (sscanf(p, "gpttmp%d.plt".ptr, &pnum) == 0) {
            gretl_remove(fname);
        }
    }

    if (alt_sty != null) {
        /* clean up alternative gnuplot stylesheet */
        g_free(alt_sty);
        alt_sty = null;
    }
}

__gshared int graph_file_written = 0;

int graph_written_to_file()
{
    return graph_file_written;
}

void set_graph_file_written(int w)
{
    graph_file_written = w;
}

__gshared int graph_file_shown = 0;

int graph_displayed()
{
    return graph_file_shown;
}

void remove_old_png()
{
    char * tmp = gretl_make_dotpath("gretltmp.png");

    gretl_remove(tmp);
    g_free(tmp);
}

/*
 * gnuplot_make_graph:
 *
 * Executes gnuplot to produce a graph: in the gretl GUI in interactive
 * mode this will be a PNG file. In the CLI program in interactive mode
 * there will be a direct call to gnuplot to display the graph. In batch
 * mode the type of file written depends on the options selected by the
 * user.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gnuplot_make_graph()
{
    char[MAXLEN] buf;
    const(char) * fname = gretl_plotfile();
    int gui = gretl_in_gui_mode();
    int fmt, err = 0;

    set_graph_file_written(0);
    graph_file_shown = 0;
    fmt = specified_gp_output_format();

    if (fmt == TermType.plt) {
        /* just the gnuplot commands are wanted */
        if (plot_output_to_buffer()) {
            err = make_plot_commands_buffer(fname);
        } else {
            set_graph_file_written(1);
        }
        return err;
    } else if (fmt == TermType.none && gui) {
        do_plot_bounding_box();
        /* ensure we don't get stale output */
        remove_old_png();
    }

    if (gui || fmt) {
        sprintf(buf.ptr, "%s \"%s\"".ptr, gretl_gnuplot_path(), fname);
    } else {
        /* gretlcli, interactive */
        sprintf(buf.ptr, "%s -persist \"%s\"".ptr, gretl_gnuplot_path(), fname);
    }
    err = gretl_spawn(buf.ptr);

    if (fmt) {
        /* got a user-specified output format */
        if (err) {
            /* leave the bad file for diagnostic purposes */
            fprintf(stderr, "err = %d: bad file is '%s'\n".ptr, err, fname);
        } else {
            /* remove the temporary input file */
            gretl_remove(fname);
            if (fmt == TermType.var) {
                graph_file_shown = 1;
            } else {
                gretl_set_path_by_name("plotfile", gnuplot_outname.ptr);
                set_graph_file_written(1);
            }
        }
    }

    return err;
}

/**
 * finalize_plot_input_file:
 * @fp: stream to which gnuplot commands have been written.
 *
 * Closes @fp and attempts to "make" the graph that it specifies.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int finalize_plot_input_file(FILE * fp)
{
    int err;

    if (fp != null) {
        fclose(fp);
        err = gnuplot_make_graph();
        if (!err) {
            /* for the benefit of gretl_cmd_exec() in interact.c,
               indicate that we actually produced a plot
            */
            set_plot_produced();
        }
    } else {
        err = 1;
    }

    return err;
}

/**
 * finalize_3d_plot_input_file:
 * @fp: stream to which gnuplot commands have been written.
 *
 * Closes @fp and alerts libgretl to the fact that an interactive
 * 3-D plot is wanted.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int finalize_3d_plot_input_file(FILE * fp)
{
    int err = 0;

    if (fp != null) {
        fclose(fp);
        set_graph_file_written(1);
    } else {
        err = 1;
    }

    return err;
}

enum GraphTitle {
    vls,
    resid,
    af,
    afv,
    factor
}

void make_gtitle(GnuplotInfo * gi, int code, const(char) * s1, const(char) * s2, FILE * fp)
{
    char[VNAMELEN] depvar;
    char * title = null;

    switch (code) {
    case GraphTitle.vls:
        if (gi.fit == FitType.ols) {
            title = g_strdup_printf(_("%s versus %s (with least squares fit)").ptr,
                                    s1, s2);
        } else if (gi.fit == FitType.inverse) {
            title = g_strdup_printf(_("%s versus %s (with inverse fit)").ptr,
                                    s1, s2);
        } else if (gi.fit == FitType.quadratic) {
            title = g_strdup_printf(_("%s versus %s (with quadratic fit)").ptr,
                                    s1, s2);
        } else if (gi.fit == FitType.cubic) {
            title = g_strdup_printf(_("%s versus %s (with cubic fit)").ptr,
                                    s1, s2);
        }
        break;
    case GraphTitle.resid:
        if (strncmp(s1, "residual for ", 13) == 0 &&
            gretl_scan_varname(s1 + 13, depvar.ptr) == 1) {
            title = g_strdup_printf(_("Regression residuals (= observed - fitted %s)").ptr,
                                    depvar.ptr);
        }
        break;
    case GraphTitle.af:
        title = g_strdup_printf(_("Actual and fitted %s").ptr, s1);
        break;
    case GraphTitle.afv:
        if (s2 == null || (gi.flags & GptFlags.ts)) {
            title = g_strdup_printf(_("Actual and fitted %s").ptr, s1);
        } else {
            title = g_strdup_printf(_("Actual and fitted %s versus %s").ptr, s1, s2);
        }
        break;
    case GraphTitle.factor:
        title = g_strdup_printf(_("%s versus %s with separation by %s").ptr,
                                s1, s2, gi.facname);
        break;
    default:
        break;
    }

    if (title != null) {
        fprintf(fp, "set title '%s'\n".ptr, title);
        g_free(title);
    }
}

void print_axis_label(char axis, const(char) * s, FILE * fp)
{
    if (strchr(s, '\'')) {
        fprintf(fp, "set %clabel \"%s\"\n".ptr, axis, s);
    } else {
        fprintf(fp, "set %clabel '%s'\n".ptr, axis, s);
    }
}

int literal_line_out(const(char) * s, int len, FILE * fp)
{
    char * q;
    char * p = cast(char *) malloc(len + 1);
    int n, warn = 0;

    if (p != null) {
        *p = '\0';
        strncat(p, s, len);
        q = p + cast(int) strspn(p, " \t");
        n = cast(int) (strlen(q));
        if (n > 0) {
            /* note: allow "set termoption ..." */
            if (!strncmp(q, "set term ", 9)) {
                warn = 1;
            } else {
                fputs(q, fp);
                if (q[n-1] != '\n') {
                    fputc('\n', fp);
                }
            }
        }
        free(p);
    }

    return warn;
}

int gnuplot_literal_from_string(const(char) * s, FILE * fp)
{
    const(char) * p;
    int braces = 1;
    int wi, warn = 0;

    s += cast(int) strspn(s, " \t{");
    p = s;

    fputs("# start literal lines\n", fp);

    while (*s) {
        if (*s == '{') {
            braces++;
        } else if (*s == '}') {
            braces--;
        }
        if (braces == 0) {
            break;
        }
        if (*s == ';') {
            wi = literal_line_out(p, cast(int) (s - p), fp);
            if (wi && !warn) {
                warn = 1;
            }
            p = s + 1;
        }
        s++;
    }

    fputs("# end literal lines\n", fp);

    return warn;
}

void print_gnuplot_literal_lines(const(char) * s, FILE * fp)
{
    if (s != null && *s != '\0') {
        gnuplot_literal_from_string(s, fp);
    }
}

void print_extra_literal_lines(char ** S, int ns, FILE * fp)
{
    int i, n;

    for (i=0; i<ns; i++) {
        if (S[i] != null) {
            n = cast(int) (strlen(S[i]));
            if (n > 0) {
                fputs(S[i], fp);
                if (S[i][n-1] != '\n') {
                    fputc('\n', fp);
                }
            }
        }
    }
}

int loess_plot(GnuplotInfo * gi, const(char) * literal, const(DatasetRaw) * dset)
{
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * x = null;
    GretlMatrixRaw * yh = null;
    int xno, yno = gi.list[1];
    const(double) * yvar = dset.Z[yno];
    const(double) * xvar;
    FILE * fp = null;
    char * title = null;
    int t, T, d = 1;
    double q = 0.5;
    int err = 0;

    if (gi.x != null) {
        xno = 0;
        xvar = gi.x;
    } else {
        xno = gi.list[2];
        xvar = dset.Z[xno];
    }

    err = graph_list_adjust_sample(gi, dset, 2);
    if (!err && gi.list[0] > 2) {
        err = GretlError.data;
    }
    if (err) {
        return err;
    }

    fp = open_plot_input_file(PlotType.regular, gi.flags, &err);
    if (err) {
        return GretlError.fopen;
    }

    err = gretl_plotfit_matrices(yvar, xvar, FitType.loess,
                                 gi.t1, gi.t2, &y, &x);

    if (!err) {
        err = sort_pairs_by_x(x, y, null, null); /* markers! */
    }

    if (!err) {
        yh = loess_fit(x, y, d, q, Gretlopt.OPT_R, &err);
    }

    if (err) {
        fclose(fp);
        goto bailout;
    }

    if (xno > 0) {
        char * s1 = make_plotname(dset, yno);
        char * s2 = make_plotname(dset, xno);

        title = g_strdup_printf(_("%s versus %s (with loess fit)").ptr, s1, s2);
        print_keypos_string(GpKeySpec.left | GpKeySpec.top, fp);
        fprintf(fp, "set title '%s'\n".ptr, title);
        g_free(title);
        print_axis_label('y', s1, fp);
        print_axis_label('x', s2, fp);
        free(s1);
        free(s2);
    } else {
        print_keypos_string(GpKeySpec.left | GpKeySpec.top, fp);
        print_axis_label('y', plotname(dset, yno, 1), fp);
    }

    print_auto_fit_string(FitType.loess, fp);

    if (literal != null) {
    print_gnuplot_literal_lines(literal, fp);
    }

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    title = g_strdup_printf(_("loess fit, d = %d, q = %g").ptr, d, q);
    fprintf(fp, " '-' using 1:2 title \"%s\" w lines\n".ptr, title);
    g_free(title);

    T = gretl_vector_get_length(yh);

    gretl_push_c_numeric_locale();

    for (t=0; t<T; t++) {
        fprintf(fp, "%.10g %.10g\n".ptr, x.val[t], y.val[t]);
    }
    fputs("e\n", fp);

    for (t=0; t<T; t++) {
        fprintf(fp, "%.10g %.10g\n".ptr, x.val[t], yh.val[t]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(x);
    gretl_matrix_free(yh);
    clear_gpinfo(gi);

    return err;
}

int get_fitted_line(GnuplotInfo * gi, const(DatasetRaw) * dset, char ** targ)
{
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * V = null;
    const(double) * yvar, xvar = null;
    double x0;
    double s2 = 0;
    double * ps2 = null;
    int allow_err = 0;
    int err = 0;

    if (gi.x != null && (dset.pd == 1 || dset.pd == 4 || dset.pd == 12)) {
        /* starting value of time index */
        x0 = gi.x[gi.t1];
    } else {
        xvar = dset.Z[gi.list[2]];
        x0 = NADBL;
    }

    yvar = dset.Z[gi.list[1]];

    if (gi.fit == FitType.none) {
        /* Doing first-time automatic OLS: we want to check for
           statistical significance of the slope coefficient
           to see if it's worth drawing the fitted line, so
           we have to allocate the variance matrix; otherwise
           we only need the coefficients.
        */
        V = gretl_matrix_alloc(2, 2);
        if (V == null) {
            return GretlError.alloc;
        }
        ps2 = &s2;
        /* if the fit attempt is automatic we'll allow it to
           fail without aborting the plot */
        allow_err = 1;
    }

    err = gretl_plotfit_matrices(yvar, xvar, gi.fit,
                                 dset.t1, dset.t2,
                                 &y, &X);

    if (!err) {
        int k = 2;

        if (gi.fit == FitType.cubic) {
            k = 4;
        } else if (gi.fit == FitType.quadratic) {
            k = 3;
        }
        b = gretl_column_vector_alloc(k);
        if (b == null) {
            err = GretlError.alloc;
        } else if (gi.fit == FitType.loglin) {
            ps2 = &s2;
        }
    }

    if (!err) {
        err = gretl_matrix_ols(y, X, b, V, null, ps2);
    }

    if (!err && gi.fit == FitType.none) {
        /* the "automatic" case */
        double pv, v = gretl_matrix_get(V, 1, 1);
        int T = gretl_vector_get_length(y);

        pv = student_pvalue_2(T - 2, b.val[1] / sqrt(v));
        /* show the line if the two-tailed p-value for the slope coeff
           is less than 0.1, otherwise discard it */
        if (pv < 0.10) {
            gi.fit = FitType.ols;
        }
    }

    if (!err && gi.fit != FitType.none) {
        GptLine line;
    memset(&line, 0, line.sizeof);
        double pd = dset.pd;

        if (gi.fit == FitType.loglin) {
            b.val[0] += s2 / 2;
        }
        set_plotfit_line(&line, gi.fit, b.val, x0, pd);
        *targ = g_strdup_printf("%s title \"%s\" w lines\n".ptr,
                                line.formula, line.title);
        g_free(line.formula);
        g_free(line.title);
        gi.flags |= GptFlags.auto_fit;
    }

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(V);

    if (err && allow_err) {
        err = 0;
    }

    return err;
}

/* support the "fit" options for a single time-series plot */

int time_fit_plot(GnuplotInfo * gi, const(char) * literal, const(DatasetRaw) * dset)
{
    int yno = gi.list[1];
    const(double) * yvar = dset.Z[yno];
    char * fitline = null;
    FILE * fp = null;
    Prn * prn;
    int t, err = 0;

    if (gi.x == null) {
        return GretlError.data;
    }

    err = graph_list_adjust_sample(gi, dset, 1);
    if (err) {
        return err;
    }

    err = get_fitted_line(gi, dset, &fitline);
    if (err) {
        return err;
    }

    gi.flags |= GptFlags.letterbox;

    fp = open_plot_input_file(PlotType.regular, gi.flags, &err);
    if (err) {
        g_free(fitline);
        return err;
    }

    prn = gretl_print_new_with_stream(fp);

    if (prn != null) {
        make_time_tics(gi, dset, 0, null, prn);
        gretl_print_detach_stream(prn);
        gretl_print_destroy(prn);
    }

    print_keypos_string(GpKeySpec.left | GpKeySpec.top, fp);
    print_axis_label('y', plotname(dset, yno, 1), fp);
    print_auto_fit_string(gi.fit, fp);
    print_gnuplot_literal_lines(literal, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w lines, \\\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, " %s".ptr, fitline);
    g_free(fitline);

    for (t=gi.t1; t<=gi.t2; t++) {
        if (gi.flags & GptFlags.timefmt) {
            fprintf(fp, "%.0f %.10g\n".ptr, gi.x[t], yvar[t]);
        } else {
            fprintf(fp, "%.10g %.10g\n".ptr, gi.x[t], yvar[t]);
        }
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    clear_gpinfo(gi);

    return err;
}

int check_tic_labels(double vmin, double vmax, GnuplotInfo * gi, char axis)
{
    char[32] s1;
    char[32] s2;
    int d, err = 0;

    for (d=6; d<12; d++) {
        sprintf(s1.ptr, "%.*g".ptr, d, vmin);
        sprintf(s2.ptr, "%.*g".ptr, d, vmax);
        if (strcmp(s1.ptr, s2.ptr)) {
            break;
        }
    }

    if (d > 6) {
        if (axis == 'x') {
            sprintf(gi.xfmt.ptr, "%% .%dg".ptr, d+1);
            if (vmax > vmin) {
                sprintf(gi.xtics.ptr, "%.*g %#.6g".ptr, d+1, vmin,
                        (vmax - vmin)/ 4.0);
            }
        } else {
            sprintf(gi.yfmt.ptr, "%% .%dg".ptr, d+1);
        }
    }

    return err;
}

void check_y_tics(GnuplotInfo * gi, const(double) ** Z, FILE * fp)
{
    double ymin, ymax;

    *gi.yfmt.ptr = '\0';

    gretl_minmax(gi.t1, gi.t2, Z[gi.list[1]], &ymin, &ymax);
    check_tic_labels(ymin, ymax, gi, 'y');

    if (*gi.yfmt.ptr != '\0') {
        fprintf(fp, "set format y \"%s\"\n".ptr, gi.yfmt.ptr);
    }
}

/* Find the minimum and maximum x-axis values and construct the gnuplot
   x-range.  We have to be a bit careful here to include only values
   that will actually display on the plot, i.e. x-values that are
   accompanied by at least one non-missing y-axis value.

   In the case of a "factorized" plot the factor variable must also
   be non-missing in order to include a given data point.

   We also have to avoid creating an "empty x range" that will choke
   gnuplot.
*/

void print_x_range_from_list(GnuplotInfo * gi, const(DatasetRaw) * dset, FILE * fp)
{
    const(double) * x, d = null;
    int k = gi.list[0];

    if (gi.flags & GptFlags.factor) {
        d = dset.Z[gi.fid];
    }

    x = dset.Z[gi.list[k]];

    if (gretl_isdummy(gi.t1, gi.t2, x)) {
        fputs("set xrange [-1:2]\n", fp);
        fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
        gi.xrange = 3;
    } else {
        double xmin, xmin0 = NADBL;
        double xmax, xmax0 = NADBL;
        int t, i, vy, obs_ok;

        for (t=gi.t1; t<=gi.t2; t++) {
            obs_ok = 0;
            if (!na(x[t]) && (d == null || !na(d[t]))) {
                for (i=1; i<k; i++) {
                    vy = gi.list[i];
                    if (!na(dset.Z[vy][t])) {
                        /* got x obs and at least one y obs */
                        obs_ok = 1;
                        break;
                    }
                }
            }
            if (obs_ok) {
                if (na(xmin0) || x[t] < xmin0) {
                    xmin0 = x[t];
                }
                if (na(xmax0) || x[t] > xmax0) {
                    xmax0 = x[t];
                }
            }
        }

        gi.xrange = xmax0 - xmin0;

        if (gi.xrange == 0.0) {
            /* construct a non-empty range */
            xmin = xmin0 - 0.5;
            xmax = xmin0 + 0.5;
        } else {
            xmin = xmin0 - gi.xrange * .025;
            if (xmin0 >= 0.0 && xmin < 0.0) {
                xmin = 0.0;
            }
            xmax = xmax0 + gi.xrange * .025;
        }

        fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);
        gi.xrange = xmax - xmin;
        check_tic_labels(xmin0, xmax0, gi, 'x');
    }
}

void print_x_range(GnuplotInfo * gi, FILE * fp)
{
    if (gretl_isdummy(gi.t1, gi.t2, gi.x)) {
        fputs("set xrange [-1:2]\n", fp);
        fputs("set xtics (\"0\" 0, \"1\" 1)\n", fp);
        gi.xrange = 3;
    } else {
        double xmin0, xmin, xmax0, xmax;

        gretl_minmax(gi.t1, gi.t2, gi.x, &xmin0, &xmax0);
        gi.xrange = xmax0 - xmin0;
        xmin = xmin0 - gi.xrange * .025;
        if (xmin0 >= 0.0 && xmin < 0.0) {
            xmin = 0.0;
        }
        xmax = xmax0 + gi.xrange * .025;
        if (xmax == xmin) {
            xmax += 0.1;
            xmin -= 0.1;
        }
        fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);
        gi.xrange = xmax - xmin;
    }
}

/* Two or more y vars plotted against some x: test to see if we want
   to use two y axes. We come here only if the user has not already
   specified one or more variables to place on the second y axis.
*/

void check_for_yscale(GnuplotInfo * gi, const(double) ** Z)
{
    double[6] ymin;
    double[6] ymax;
    double ratio;
    int lmax = gi.list[0];
    int i, j, oddcount;
    int oddman = 0;

    if (gi.flags & GptFlags.idx) {
        /* do this only if we haven't added a 0 at the end of
           the list */
        if (gi.list[lmax] != 0) {
            lmax++;
        }
    }

    /* find minima, maxima of the y-axis vars */
    for (i=1; i<lmax; i++) {
        gretl_minmax(gi.t1, gi.t2, Z[gi.list[i]],
                     &ymin[i-1], &ymax[i-1]);
    }

    gi.flags &= ~GptFlags.y2axis;

    for (i=lmax-1; i>0; i--) {
        oddcount = 0;
        for (j=1; j<lmax; j++) {
            if (j != i) {
                ratio = ymax[i-1] / ymax[j-1];
                if (ratio > 5.0 || ratio < 0.2) {
                    gi.flags |= GptFlags.y2axis;
                    oddcount++;
                }
            }
        }
        if (oddcount == lmax - 2) {
            /* the series at list position i differs considerably
               in scale from all the others in the list
            */
            oddman = i;
            break;
        }
    }

    if (oddman == 0) {
        gi.flags &= ~GptFlags.y2axis;
    } else {
        gi.y2vars = gretl_list_new(1);
        gi.y2vars[1] = gi.list[oddman]; /* ? */
        fprintf(stderr, "HERE oddman = %d, gi->list[%d] = %d\n".ptr,
                oddman, oddman, gi.list[oddman]);
    }
}

int get_y2_vars(GnuplotInfo * gi, const(DatasetRaw) * dset)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_D);
    const(int) * plist;
    int i, pos;
    int err = 0;

    if ((plist = get_list_by_name(s)) != null) {
        /* got a names list */
        for (i=1; i<=plist[0]; i++) {
            if (!in_gretl_list(gi.list, plist[i])) {
                err = GretlError.invarg;
                break;
            }
        }
        if (!err) {
            gi.y2vars = gretl_list_copy(plist);
        }
    } else if (strchr(s, ',')) {
        /* comma-separated terms */
        char ** strs = g_strsplit(s, ",", 0);

        for (i=0; strs[i] != null; i++) {
            pos = gp_list_pos(strs[i], gi.list, dset);
            if (pos <= 0) {
                err = GretlError.invarg;
            } else {
                gretl_list_append_term(&gi.y2vars, gi.list[pos]);
            }
        }
        g_strfreev(strs);
    } else {
        /* a single term */
        gi.y2vars = gretl_list_new(1);
        pos = gp_list_pos(s, gi.list, dset);
        if (pos <= 0) {
            err = GretlError.invarg;
        } else {
            gi.y2vars[1] = gi.list[pos];
        }
    }

    return err;
}

int is_y2_var(GnuplotInfo * gi, int i)
{
    if (gi.y2vars != null) {
        return in_gretl_list(gi.y2vars, gi.list[i]);
    } else {
        return 0;
    }
}

int print_gp_factorized_data(GnuplotInfo * gi, const(DatasetRaw) * dset, FILE * fp)
{
    const(double) * d = null;
    const(double) * y = null;
    const(double) * x = null;
    double yjt;
    int len;
    int i, j, vj, t;

    len = gi.list[0];
    d = dset.Z[gi.fid];
    x = dset.Z[gi.list[len]];

    for (j=1; j<len; j++) {
        vj = gi.list[j];
        y = dset.Z[vj];
        if (!y_var_is_factorized(gi, vj)) {
            for (t=gi.t1; t<=gi.t2; t++) {
                if (na(x[t]) || na(d[t]) || na(y[t])) {
                    continue;
                }
                fprintf(fp, "%.10g %.10g".ptr, x[t], y[t]);
                if (dset.markers) {
                    fprintf(fp, " # %s".ptr, dset.S[t]);
                } else if (dataset_is_time_series(dset)) {
                    char[OBSLEN] obs;

                    ntolabel(obs.ptr, t, dset);
                    fprintf(fp, " # %s".ptr, obs.ptr);
                }
                fputc('\n', fp);
            }
            fputs("e\n", fp);
            continue;
        }
        for (i=0; i<gi.n_fvals; i++) {
            for (t=gi.t1; t<=gi.t2; t++) {
                yjt = (d[t] == gi.fvals.val[i])? y[t] : NADBL;
                if (na(x[t]) || na(d[t]) || na(yjt)) {
                    continue;
                }
                fprintf(fp, "%.10g %.10g".ptr, x[t], y[t]);
                if (dset.markers) {
                    fprintf(fp, " # %s".ptr, dset.S[t]);
                } else if (dataset_is_time_series(dset)) {
                    char[OBSLEN] obs;

                    ntolabel(obs.ptr, t, dset);
                    fprintf(fp, " # %s".ptr, obs.ptr);
                }
                fputc('\n', fp);
            }
            fputs("e\n", fp);
        }
    }

    return 0;
}

/* for printing panel time-series graph: insert a discontinuity
   between the panel units */

void maybe_print_panel_jot(int t, const(DatasetRaw) * dset, FILE * fp)
{
    char[OBSLEN] obs;
    int maj, min;

    ntolabel(obs.ptr, t, dset);
    sscanf(obs.ptr, "%d:%d".ptr, &maj, &min);
    if (maj > 1 && min == 1) {
        fprintf(fp, "%g %s\n".ptr, t + 0.5, GPNA.ptr);
    }
}

/* sanity check for totally empty graph */

int all_graph_data_missing(const(int) * list, int t, const(double) ** Z)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        if (!na(Z[list[i]][t])) {
            return 0;
        }
    }

    return 1;
}

int use_impulses(GnuplotInfo * gi)
{
    if (gi.withlist != null) {
        int i;

        for (i=1; i<=gi.withlist[0]; i++) {
            if (gi.withlist[i] == WithStyle.impulses) {
                return 1;
            }
        }
    }

    return 0;
}

int use_lines(GnuplotInfo * gi)
{
    if (gi.withlist != null) {
        int i;

        for (i=1; i<=gi.withlist[0]; i++) {
            if (gi.withlist[i] == WithStyle.lines) {
                return 1;
            }
        }
    }

    return 0;
}

/* list of series IDs for which we should skip observations
   with NAs when printing the plot data */
__gshared int * na_skiplist = null;

void print_gp_data(GnuplotInfo * gi, const(DatasetRaw) * dset, FILE * fp)
{
    int n = gi.t2 - gi.t1 + 1;
    double offset = 0.0;
    int[3] datlist;
    int lmax, ynum = 2;
    int nomarkers = 0;
    int i, t;

    /* multi impulse plot? calculate offset for lines */
    if (use_impulses(gi) && gi.list[0] > 2) {
        offset = 0.10 * gi.xrange / n;
    }

    if (gi.x != null) {
        lmax = gi.list[0] - 1;
        datlist[0] = 1;
        ynum = 1;
    } else {
        lmax = gi.list[0] - 1;
        datlist[0] = 2;
        datlist[1] = gi.list[gi.list[0]];
    }

    if (use_impulses(gi) || use_lines(gi)) {
        nomarkers = 1;
    }

    /* loop across the variables, printing x then y[i] for each i */

    for (i=1; i<=lmax; i++) {
        double xoff = offset * (i - 1);

        datlist[ynum] = gi.list[i];

        for (t=gi.t1; t<=gi.t2; t++) {
            const(char) * label = null;
            char[OBSLEN] obs;

            if (in_gretl_list(na_skiplist, datlist[ynum]) &&
                na(dset.Z[datlist[ynum]][t])) {
                continue;
            } else if (gi.x == null &&
                       all_graph_data_missing(gi.list, t, cast(const(double) **) dset.Z)) {
                continue;
            }
            if (!(gi.flags & GptFlags.ts) && i == 1) {
                if (dset.markers) {
                    label = dset.S[t];
                } else if (!nomarkers && dataset_is_time_series(dset)) {
                    ntolabel(obs.ptr, t, dset);
                    label = obs.ptr;
                }
            }
            if ((gi.flags & GptFlags.ts) && dset.structure == TsCodes.stacked_time_series) {
                maybe_print_panel_jot(t, dset, fp);
            }
            printvars(fp, t, datlist.ptr, dset, gi, label, xoff);
        }

        fputs("e\n", fp);
    }
}

int gpinfo_init(GnuplotInfo * gi, Gretlopt opt, const(int) * list, const(char) * literal, const(DatasetRaw) * dset)
{
    int l0 = list[0];
    int err = 0;

    gi.withlist = null;
    gi.yformula = null;
    gi.facname = null;
    gi.x = null;
    gi.list = null;
    gi.fvals = null;
    gi.flist = null;
    gi.y2vars = null;
    gi.n_fvals = 0;
    gi.fid = 0;
    gi.band = 0;

    err = get_gp_flags(gi, opt, list, dset);
    if (err) {
        return err;
    }

    if (gi.band) {
        /* force single y-axis in "band" case */
        opt |= Gretlopt.OPT_Y;
    }

    if (gi.fit == FitType.none) {
        gi.flags |= GptFlags.ts; /* may be renounced later */
    }

    if (dset.t2 - dset.t1 + 1 <= 0) {
        /* null sample range */
        return GretlError.data;
    }

    gi.t1 = dset.t1;
    gi.t2 = dset.t2;
    gi.xrange = 0.0;
    gi.xtics[0] = '\0';
    gi.xfmt[0] = '\0';
    gi.yfmt[0] = '\0';

    if (l0 < 2 && !(gi.flags & GptFlags.idx)) {
        return GretlError.args;
    }

    if ((gi.flags & GptFlags.factor) && (gi.flags & GptFlags.idx)) {
        return GretlError.badopt;
    }

    gi.list = gretl_list_copy(list);
    if (gi.list == null) {
        return GretlError.alloc;
    }
    if (gi.flags & GptFlags.factor) {
        /* record the factor var and cut it off the plot list */
        gi.fid = gi.list[gi.list[0]];
        gi.facname = dset.varname[gi.fid];
        gi.list[0] -= 1;
    }

    if (opt & Gretlopt.OPT_D) {
        /* got explicit --y2axis spec */
        gi.flags |= GptFlags.y2axis;
    } else if ((l0 > 2 || (l0 > 1 && (gi.flags & GptFlags.idx))) &&
        l0 < 7 && !(gi.flags & GptFlags.resids) && !(gi.flags & GptFlags.fa)
        && !(gi.flags & GptFlags.factor) && !(opt & Gretlopt.OPT_Y)) {
        /* FIXME GptFlags.xyz ? */
        /* allow probe for using two y axes */
        gi.flags |= GptFlags.y2axis;
    }

    if ((gi.flags & GptFlags.fa) && literal != null &&
        !strncmp(literal, "yformula: ", 10)) {
        /* fitted vs actual plot with fitted given by formula */
        gi.yformula = literal + 10;
    }

    if (literal != null && strstr(literal, "set style data")) {
        gi.flags |= GptFlags.data_style;
    }

    return 0;
}

void clear_gpinfo(GnuplotInfo * gi)
{
    free(gi.list);
    free(gi.flist);
    gretl_matrix_free(gi.fvals);
    free(gi.withlist);
    free(gi.y2vars);
}

void set_lwstr(char * s)
{
    if (default_png_scale > 1.0) {
        strcpy(s, " lw 2");
    } else {
        *s = '\0';
    }
}

void set_plot_withstr(GnuplotInfo * gi, int i, char * str)
{
    if (gi.flags & GptFlags.data_style) {
        *str = '\0';
    } else if (gi.withlist != null) {
        int withval = WithStyle.points;

        if (i > 0 && i <= gi.withlist[0]) {
            withval = gi.withlist[i];
        }

        if (withval == WithStyle.lines) {
            strcpy(str, "w lines");
        } else if (withval == WithStyle.impulses) {
            strcpy(str, "w impulses");
        } else if (withval == WithStyle.lp) {
            strcpy(str, "w linespoints");
        } else if (withval == WithStyle.boxes) {
            strcpy(str, "w boxes");
        } else if (withval == WithStyle.steps) {
            strcpy(str, "w steps");
        } else {
            strcpy(str, "w points");
        }
    } else if (gi.flags & GptFlags.lines) {
        strcpy(str, "w lines");
    } else if (gi.flags & GptFlags.impulses) {
        strcpy(str, "w impulses");
    } else if (gi.flags & GptFlags.steps) {
        strcpy(str, "w steps");
    } else {
        strcpy(str, "w points");
    }
}

int graph_list_adjust_sample(GnuplotInfo * gi, const(DatasetRaw) * dset, int listmin)
{
    int t1min = gi.t1;
    int t2max = gi.t2;
    int i, t, vi, t_ok;
    int err = 0;

    for (t=t1min; t<=t2max; t++) {
        t_ok = 0;
        for (i=1; i<=gi.list[0]; i++) {
            vi = gi.list[i];
            if (vi > 0 && !na(dset.Z[vi][t])) {
                t_ok = 1;
                break;
            }
        }
        if (t_ok) {
            break;
        }
        t1min++;
    }

    for (t=t2max; t>t1min; t--) {
        t_ok = 0;
        for (i=1; i<=gi.list[0]; i++) {
            vi = gi.list[i];
            if (vi > 0 && !na(dset.Z[vi][t])) {
                t_ok = 1;
                break;
            }
        }
        if (t_ok) {
            break;
        }
        t2max--;
    }

    if (t2max > t1min) {
        for (i=1; i<=gi.list[0]; i++) {
            int all_missing = 1;

            vi = gi.list[i];
            for (t=t1min; t<=t2max; t++) {
                if (!na(dset.Z[vi][t])) {
                    all_missing = 0;
                    break;
                }
            }
            if (all_missing) {
                gretl_list_delete_at_pos(gi.list, i);
                i--;
            }
        }
    }

    gi.t1 = t1min;
    gi.t2 = t2max;

    if (gi.t1 > gi.t2 || gi.list[0] < listmin) {
        err = GretlError.data;
        gretl_errmsg_set(_("No data are available for plotting").ptr);
    }

    return err;
}

int timefmt_useable(const(DatasetRaw) * dset)
{
    return dated_daily_data(dset) || dated_weekly_data(dset);
}

int maybe_add_plotx(GnuplotInfo * gi, int time_fit, const(DatasetRaw) * dset)
{
    Gretlopt xopt = Gretlopt.OPT_NONE;
    int k = gi.list[0];
    int add0 = 0;

    /* are we really doing a time-series plot? */
    if (k > 1 && !strcmp(dset.varname[gi.list[k]], "time")) {
        ; /* yes */
    } else if (gi.flags & GptFlags.idx) {
        add0 = 1; /* yes */
    } else {
        /* no: get out */
        gi.flags &= ~GptFlags.ts;
        return 0;
    }

    if (!time_fit && timefmt_useable(dset)) {
        /* experimental */
        gi.flags |= GptFlags.timefmt;
        xopt = Gretlopt.OPT_T;
    }

    gi.x = gretl_plotx(dset, xopt);
    if (gi.x == null) {
        return GretlError.alloc;
    }

    /* a bit nasty, but add a dummy list entry for
       the 'virtual' plot variable
    */
    if (add0) {
        gretl_list_append_term(&gi.list, 0);
        if (gi.list == null) {
            return GretlError.alloc;
        }
    }

    return 0;
}

void gnuplot_missval_string(FILE * fp)
{
    fputs("set datafile missing \"?\"\n", fp);
}

void graph_month_name(char * mname, int m)
{
    tm mt;

    mt.tm_sec = 0;
    mt.tm_min = 0;
    mt.tm_hour = 0;
    mt.tm_mday = 1;
    mt.tm_mon = m - 1;
    mt.tm_year = 100;

    strftime(mname, 7, "%b", &mt);

    if (!g_utf8_validate(mname, -1, null)) {
        /* we might be in a non-UTF-8 locale */
        char * tmp;
        gsize bytes;

        tmp = g_locale_to_utf8(mname, -1, null, &bytes, null);
        if (tmp != null) {
            strcpy(mname, tmp);
            g_free(tmp);
        }
    }
}

/* for short daily time-series plots: write month names
   into the xtics */

void make_named_month_tics(const(GnuplotInfo) * gi, double yrs, Prn * prn)
{
    double t0 = gi.x[gi.t1];
    double t1 = gi.x[gi.t2];
    double x, tw = 1.0/12;
    int i, m, n = 0;
    char[16] mname;
    int notfirst = 0;
    int scale = cast(int) (yrs * 1.5);

    if (scale == 0) {
        scale = 1;
    }

    t0 += (1.0 - (t0 - floor(t0)) * 12.0) / 12.0;
    for (x=t0; x<t1; x+=tw) n++;

    x = (t0 - floor(t0)) * 12;
    m = cast(int) (1 + ((x - floor(x) > .8)? ceil(x) : floor(x)));
    if (m > 12) m -= 12;

    pputs(prn, "set xtics (");
    x = t0;

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        if (m == 1) {
            if (notfirst) {
                pputs(prn, ", ");
            }
            pprintf(prn, "\"%4.0f\" %.10g".ptr, x, x);
            notfirst = 1;
        } else if ((scale == 1) || (m % scale == 1)) {
            graph_month_name(mname.ptr, m);
            if (notfirst) {
                pputs(prn, ", ");
            }
            pprintf(prn, "\"%s\" %.10g".ptr, mname.ptr, x);
            notfirst = 1;
        }
        m++;
        x += tw;
        if (m > 12) m -= 12;
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

int continues_unit(const(DatasetRaw) * dset, int t)
{
    return t / dset.pd == (t-1) / dset.pd;
}

/* Below: we're making a combined time series plot for panel data.
   That is, time series for unit 1, followed by time series for unit
   2, etc.  We'd like to show tic marks to represent the start of each
   unit's time series, but we have to watch out for the case where
   there are "too many" units -- we don't want a dense fudge of marks
   on the x-axis.  In that case we put a tic mark only for every k'th
   unit.
*/

void make_panel_unit_tics(const(DatasetRaw) * dset, GnuplotInfo * gi, Prn * prn)
{
    int maxtics, ticskip;
    double ntics;
    int printed;
    int u, t, n;

    pputs(prn, "set xtics (");

    gretl_push_c_numeric_locale();

    /* how many panel units are included in the plot? */
    maxtics = gi.t2 / dset.pd - gi.t1 / dset.pd + 1;

    ntics = maxtics;
    while (ntics > 20) {
        ntics /= 1.5;
    }

    ticskip = cast(int) (maxtics / ceil(ntics));

    if (ticskip == 1 && ntics < maxtics) {
        /* otherwise we'll get an incomplete scale */
        ntics = maxtics;
    }

    n = printed = 0;
    u = gi.t1 / dset.pd;

    for (t=gi.t1; t<=gi.t2 && printed<ntics; t++) {
        if (t == gi.t1 || !continues_unit(dset, t)) {
            u++;
            if (n % ticskip == 0) {
                pprintf(prn, "\"%d\" %.10g".ptr, u, gi.x[t]);
                if (++printed < ntics) {
                    pputs(prn, ", ");
                }
            }
            n++;
        }
    }

    gretl_pop_c_numeric_locale();

    pputs(prn, ")\n");
}

void make_calendar_tics(const(DatasetRaw) * dset, const(GnuplotInfo) * gi, Prn * prn)
{
    int T = gi.t2 - gi.t1 + 1;
    double yrs;

    if (dset.pd == 52) {
        yrs = T / 52.0;
    } else {
        yrs = T / (dset.pd * 52.0);
    }

    if (yrs <= 3) {
        make_named_month_tics(gi, yrs, prn);
    } else if (yrs < 6) {
        /* don't show ugly "fractions of years" */
        pputs(prn, "set xtics 1\n");
        if (yrs < 3) {
            /* put monthly minor tics */
            pputs(prn, "set mxtics 12\n");
        } else if (yrs < 5) {
            /* quarterly minor tics */
            pputs(prn, "set mxtics 4\n");
        }
    }
}

int multiple_groups(const(DatasetRaw) * dset, int t1, int t2)
{
    int ret = 0;

    if (dataset_is_panel(dset)) {
        ret = (t2 / dset.pd > t1 / dset.pd);
    }

    return ret;
}

int single_year_sample(const(DatasetRaw) * dset, int t1, int t2)
{
    char[OBSLEN] obs;
    int y1, y2;

    ntolabel(obs.ptr, t1, dset);
    y1 = atoi(obs.ptr);
    ntolabel(obs.ptr, t2, dset);
    y2 = atoi(obs.ptr);

    return y2 == y1;
}

/* Use year-and-month major tics when plotting fewer than ??
   monthly observations.
*/

/* stands in for g_date_time_new_local() + g_date_time_format(dt, "%b %Y") */

char * month_year_string(int y, int m)
{
    char[32] buf;
    tm mt;

    memset(&mt, 0, mt.sizeof);
    mt.tm_year = y - 1900;
    mt.tm_mon = m - 1;
    mt.tm_mday = 1;
    strftime(buf.ptr, buf.length, "%b %Y", &mt);
    return g_strdup(buf.ptr);
}

void short_monthly_tics(GnuplotInfo * gi, int T, const(DatasetRaw) * dset, Prn * prn)
{
    int dy = 0, dm = 0;
    char * tstr = null;
    int use_names = 1;
    int ticskip = 1;
    int y, m, t, j;
    double x;

    if (T > 12) {
        ticskip += (T > 12) + (T > 24);
    }

    date_maj_min(gi.t1, dset, &y, &m);
    if (use_names) {
        dy = y;
        dm = m;
    }

    pprintf(prn, "# xtics lines = %d\n".ptr, T);
    pprintf(prn, "set xtics (".ptr);
    gretl_push_c_numeric_locale();

    for (t=gi.t1, j=0; t<=gi.t2; t++, j++) {
        x = y + (m - 1) / 12.0;
        if (j % ticskip == 0) {
            /* major (labeled) tic */
            if (use_names) {
                tstr = month_year_string(dy, dm);
                pprintf(prn, "\"%s\" %g 0".ptr, tstr, x);
                dm += ticskip;
                while (dm > 12) {
                    dm -= 12;
                    dy++;
                }
                g_free(tstr);
            } else {
                pprintf(prn, "\"%d-%02d\" %g 0".ptr, y, m, x);
            }
        } else {
            /* minor (unlabeled) tic */
            pprintf(prn, "\"\" %g 1".ptr, x);
        }
        pputs(prn, t < gi.t2 ? ", \\\n" : ")\n");
        if (m == 12) {
            y++;
            m = 1;
        } else {
            m++;
        }
    }

    gretl_pop_c_numeric_locale();
    if (T > 8) {
        pputs(prn, "set xtics rotate by -45\n");
    }
}

/* special tics for time series plots */

void make_time_tics(GnuplotInfo * gi, const(DatasetRaw) * dset, int many, char * xlabel, Prn * prn)
{
    int T = gi.t2 - gi.t1 + 1;
    int few_years = 8;

    if (many) {
        pprintf(prn, "# multiple timeseries %d\n".ptr, dset.pd);
    } else {
        pprintf(prn, "# timeseries %d".ptr, dset.pd);
        gi.flags |= GptFlags.letterbox;
        pputs(prn, " (letterbox)\n");
    }

    if (gi.flags & GptFlags.timefmt) {
        double yr100 = -59011441438;

        pputs(prn, "set xdata time\n");
        pputs(prn, "set timefmt \"%s\"\n");
        if (single_year_sample(dset, gi.t1, gi.t2)) {
            strcpy(gi.xfmt.ptr, "%m-%d");
        } else if (gi.x != null && gi.x[gi.t2] < yr100) {
            /* century not known */
            strcpy(gi.xfmt.ptr, "%y-%m-%d");
        } else {
            strcpy(gi.xfmt.ptr, "%Y-%m-%d");
        }
        pprintf(prn, "set format x \"%s\"\n".ptr, gi.xfmt.ptr);
        pputs(prn, "set xtics rotate by -45\n");
        pputs(prn, "set xtics nomirror\n");
        return;
    }

    if (dset.pd == 1 && T < few_years) {
        pputs(prn, "set xtics nomirror 1\n");
    } else if (dset.pd == 4 && T / 4 < few_years) {
        pputs(prn, "set xtics nomirror 0,1\n");
        pputs(prn, "set mxtics 4\n");
    } else if (dset.pd == 12) {
        if (T < 36) {
            short_monthly_tics(gi, T, dset, prn);
        } else if (T / 12 < few_years) {
            pputs(prn, "set xtics nomirror 0,1\n");
            pputs(prn, "set mxtics 12\n");
        }
    } else if (dated_daily_data(dset) || dated_weekly_data(dset)) {
        make_calendar_tics(dset, gi, prn);
    } else if (multiple_groups(dset, gi.t1, gi.t2)) {
        make_panel_unit_tics(dset, gi, prn);
        if (xlabel != null) {
            strcpy(xlabel, _("time series by group").ptr);
        }
    }
}

/* Handle the use of a matrix in the context of the "plot" command
   block, or create a plot directly from a matrix and a plot list.
*/

int matrix_plot(GretlMatrixRaw * m, const(int) * list, const(char) * literal, Gretlopt opt)
{
    DatasetRaw * dset = null;
    int * plotlist = null;
    int pmax, err = 0;

    if (gretl_is_null_matrix(m)) {
        err = GretlError.data;
    } else if (list != null && list[0] == 0) {
        /* let an empty list be equivalent to NULL */
        dset = gretl_dataset_from_matrix(m, null, Gretlopt.OPT_B, &err);
    } else {
        dset = gretl_dataset_from_matrix(m, list, Gretlopt.OPT_B, &err);
    }
    if (err) {
        return err;
    }

    pmax = dset.v - 1;
    if (pmax <= 0) {
        err = GretlError.data;
    } else {
        plotlist = gretl_consecutive_list_new(1, pmax);
        if (plotlist == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        if (opt & (Gretlopt.OPT_N | Gretlopt.OPT_A)) {
            /* --band or --bands */
            GnuplotInfo gi;

            err = gpinfo_init(&gi, opt, plotlist, literal, dset);
            if (!err) {
                err = maybe_add_plotx(&gi, 0, dset);
            }
            if (!err) {
                err = plot_with_band(BPMode.blockmat, &gi, literal,
                                     dset, opt);
            }
        } else {
            err = gnuplot(plotlist, literal, dset, opt);
        }
    }

    destroy_dataset(dset);
    free(plotlist);

    return err;
}

int plotlist_is_group_invariant(const(int) * list, const(DatasetRaw) * dset)
{
    int i;

    for (i=1; i<=list[0]; i++) {
        if (!series_is_group_invariant(dset, list[i])) {
            return 0;
        }
    }

    return 1;
}

int panel_group_invariant_plot(const(int) * plotlist, const(char) * literal, DatasetRaw * dset, Gretlopt opt)
{
    DatasetRaw orig = *dset;
    int err;

    /* limit sample to first group */
    dset.t1 = 0;
    dset.t2 = dset.pd - 1;

    /* and mark as time-series data */
    if (dset.panel_pd > 0) {
        dset.pd = dset.panel_pd;
        dset.sd0 = dset.panel_sd0;
        dset.structure = TsCodes.time_series;
    } else {
        dset.structure = TsCodes.special_time_series;
        dset.pd = 1;
    }

    err = gnuplot(plotlist, literal, dset, opt);

    /* put everything back as it was */
    *dset = orig;

    return err;
}

int time_fit_wanted(Gretlopt * popt)
{
    if ((*popt & Gretlopt.OPT_T) && (*popt & Gretlopt.OPT_F)) {
        const(char) * s = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_F);

        if (s == null || !strcmp(s, "none")) {
            *popt = cast(Gretlopt) (*popt & ~Gretlopt.OPT_F);
        } else {
            return 1;
        }
    }

    return 0;
}

/**
 * gnuplot:
 * @plotlist: list of variables to plot, by ID number.
 * @literal: commands to be passed to gnuplot.
 * @dset: dataset struct.
 * @opt: option flags.
 *
 * Writes a gnuplot plot file to display the values of the
 * variables in @list and calls gnuplot to make the graph.
 * The option flags available in hansl are set out in
 * options.c (see GNUPLOT). Three additional flags are
 * supported for use in C code: OPT_S is equivalent to
 * the hansl option --fit=none; OPT_A and OPT_S signal
 * a residual plot and a fitted/actual plot, respectively.
 *
 * Returns: 0 on successful completion, non-zero code on error.
 */

int gnuplot(const(int) * plotlist, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    Prn * prn = null;
    FILE * fp = null;
    char[MAXDISP] s1 = 0;
    char[MAXDISP] s2 = 0;
    char[MAXDISP] xlabel = 0;
    char[16] withstr = 0;
    char[8] lwstr = 0;
    char[48] keystr = 0;
    char * fitline = null;
    int time_fit = 0;
    int many = 0;
    int set_xrange = 1;
    PlotType ptype;
    GnuplotInfo gi;
    int i, err = 0;

    gretl_error_clear();

    if (time_fit_wanted(&opt)) {
        if (plotlist[0] > 1 || !dataset_is_time_series(dset)) {
            return GretlError.badopt;
        } else {
            time_fit = 1;
        }
    }

    if (literal != null && strstr(literal, "set xdata time")) {
        set_xrange = 0;
    }

    if (dataset_is_panel(dset) &&
        plotlist_is_group_invariant(plotlist, dset)) {
        return panel_group_invariant_plot(plotlist, literal,
                                          cast(DatasetRaw *) dset, opt);
    }

    err = gpinfo_init(&gi, opt, plotlist, literal, dset);
    if (err) {
        goto bailout;
    }

    if ((opt & Gretlopt.OPT_D) && gi.list[0] > 1) {
        err = get_y2_vars(&gi, dset);
        if (err) {
            goto bailout;
        }
    }

    err = maybe_add_plotx(&gi, time_fit, dset);
    if (err) {
        goto bailout;
    }

    /* hive off some special cases */
    if (gi.fit == FitType.loess) {
        return loess_plot(&gi, literal, dset);
    } else if (time_fit) {
        return time_fit_plot(&gi, literal, dset);
    } else if (gi.band) {
        return plot_with_band(BPMode.regular, &gi, literal,
                              cast(DatasetRaw *) dset, opt);
    }

    if (gi.list[0] > MAX_LETTERBOX_LINES + 1) {
        many = 1;
    }

    prn = gretl_print_new(PrnType.buffer, &err);
    if (err) {
        goto bailout;
    }

    /* adjust sample range, and reject if it's empty */
    err = graph_list_adjust_sample(&gi, dset, 2);
    if (err) {
        goto bailout;
    }

    /* add a regression line if appropriate */
    if (!use_impulses(&gi) && !(gi.flags & GptFlags.fit_omit) && gi.list[0] == 2 &&
        !(gi.flags & GptFlags.ts) && !(gi.flags & GptFlags.resids)) {
        err = get_fitted_line(&gi, dset, &fitline);
        if (err) {
            goto bailout;
        } else {
            const(char) * xname = dset.varname[gi.list[2]];
            const(char) * yname = dset.varname[gi.list[1]];

            if (*xname != '\0' && *yname != '\0') {
                pprintf(prn, "# X = '%s' (%d)\n".ptr, xname, gi.list[2]);
                pprintf(prn, "# Y = '%s' (%d)\n".ptr, yname, gi.list[1]);
            }
        }
    }

    ptype = PlotType.regular;

    /* separation by factor: create special vars */
    if (gi.flags & GptFlags.factor) {
        err = factor_check(&gi, dset);
        if (err) {
            goto bailout;
        }
        ptype = PlotType.factorized;
    }

    /* special tics for time series plots */
    if (gi.flags & GptFlags.ts) {
        make_time_tics(&gi, dset, many, xlabel.ptr, prn);
    }

    /* open file and, if that goes OK, dump the prn into it
       after writing the header
    */
    fp = open_plot_input_file(cast(PlotType) ptype, gi.flags, &err);
    if (err) {
        gretl_print_destroy(prn);
        goto bailout;
    }

    fputs(gretl_print_get_buffer(prn), fp);
    gretl_print_destroy(prn);

    /* set x-axis label for non-time series plots */
    if (!(gi.flags & GptFlags.ts)) {
        int vx = gi.list[gi.list[0]];

        strcpy(xlabel.ptr, plotname(dset, vx, 1));
    }

    print_axis_label('x', xlabel.ptr, fp);
    fputs("set xzeroaxis\n", fp);
    gnuplot_missval_string(fp);

    if (gi.flags & GptFlags.logy) {
        fprintf(fp, "set logscale y %g\n".ptr, gi.ybase);
    }

    /* key: default to left top */
    strcpy(keystr.ptr, "set key left top\n");

    if (gi.list[0] == 1) {
        /* only one variable (time series) */
        print_axis_label('y', plotname(dset, gi.list[1], 1), fp);
        strcpy(keystr.ptr, "set nokey\n");
    } else if (gi.list[0] == 2) {
        /* plotting two variables */
        int no_key = 1;

        if (gi.flags & GptFlags.auto_fit) {
            print_auto_fit_string(gi.fit, fp);
            if (gi.flags & GptFlags.fa) {
                char * tmp1 = make_plotname(dset, gi.list[1]);
                char * tmp2 = make_plotname(dset, gi.list[2]);

                make_gtitle(&gi, GraphTitle.afv, tmp1, tmp2, fp);
                free(tmp1);
                free(tmp2);
            } else {
                make_gtitle(&gi, GraphTitle.vls, plotname(dset, gi.list[1], 1),
                            xlabel.ptr, fp);
            }
            no_key = 0;
        } else if (gi.flags & GptFlags.factor) {
            char * tmp1 = make_plotname(dset, gi.list[1]);
            char * tmp2 = make_plotname(dset, gi.list[2]);

            make_gtitle(&gi, GraphTitle.factor, tmp1, tmp2, fp);
            no_key = 0;
        }
        if (gi.flags & GptFlags.resids) {
            const(char) * vlabel = series_get_label(dset, gi.list[1]);

            make_gtitle(&gi, GraphTitle.resid, vlabel == null ? "residual" : vlabel,
                        null, fp);
            fprintf(fp, "set ylabel '%s'\n".ptr, _("residual").ptr);
        } else {
            print_axis_label('y', plotname(dset, gi.list[1], 1), fp);
        }
        if (no_key) {
            strcpy(keystr.ptr, "set nokey\n");
        }
    } else if ((gi.flags & GptFlags.resids) && (gi.flags & GptFlags.factor)) {
        const(char) * vlabel = series_get_label(dset, gi.list[1]);

        make_gtitle(&gi, GraphTitle.resid, vlabel == null ? "residual" : vlabel,
                    null, fp);
        fprintf(fp, "set ylabel '%s'\n".ptr, _("residual").ptr);
    } else if (gi.flags & GptFlags.fa) {
        if (gi.list[3] == dset.v - 1) {
            /* x var is just time or index: is this always right? */
            make_gtitle(&gi, GraphTitle.af, plotname(dset, gi.list[2], 1),
                        null, fp);
        } else {
            char * s2p = make_plotname(dset, gi.list[2]);
            char * s3 = make_plotname(dset, gi.list[3]);

            make_gtitle(&gi, GraphTitle.afv, s2p, s3, fp);
            free(s2p);
            free(s3);
        }
        print_axis_label('y', plotname(dset, gi.list[2], 1), fp);
    } else if (gi.flags & GptFlags.factor) {
        if (gi.list[0] == 2) {
            print_axis_label('y', plotname(dset, gi.list[1], 1), fp);
        }
    }

    if (many) {
        strcpy(keystr.ptr, "set key outside\n");
    }

    fputs(keystr.ptr, fp);

    gretl_push_c_numeric_locale();

    if (set_xrange) {
        if (gi.x != null) {
            print_x_range(&gi, fp);
        } else {
            print_x_range_from_list(&gi, dset, fp);
        }
    }

    if (!(gi.flags & GptFlags.timefmt) && *gi.xfmt.ptr != '\0' && *gi.xtics.ptr != '\0') {
        /* remedial handling of broken tics */
        fprintf(fp, "set format x \"%s\"\n".ptr, gi.xfmt.ptr);
        fprintf(fp, "set xtics %s\n".ptr, gi.xtics.ptr);
    }

    if (gi.flags & GptFlags.y2axis) {
        if (gi.y2vars == null) {
            check_for_yscale(&gi, cast(const(double) **) dset.Z);
        }
        if (gi.flags & GptFlags.y2axis) {
            fputs("set ytics nomirror\n", fp);
            fputs("set y2tics\n", fp);
        }
    } else if (gi.yformula == null && gi.list[0] == 2) {
        check_y_tics(&gi, cast(const(double) **) dset.Z, fp);
    }

    if (gi.yformula != null) {
        /* cut out the "dummy" yvar that is in fact represented
           by a formula rather than raw data */
        gi.list[1] = gi.list[2];
        gi.list[2] = gi.list[3];
        gi.list[0] = 2;
    } else {
        print_gnuplot_literal_lines(literal, fp);
    }

    /* now print the 'plot' lines */
    fputs("plot \\\n", fp);
    if (gi.flags & GptFlags.y2axis) {
        /* using two y axes */
        int lmax = gi.list[0];

        if ((gi.flags & GptFlags.idx) && gi.list[lmax] != 0) {
            lmax++;
        }
        for (i=1; i<lmax; i++) {
            set_lwstr(lwstr.ptr);
            set_plot_withstr(&gi, i, withstr.ptr);
            fprintf(fp, " '-' using 1:2 axes %s title \"%s (%s)\" %s%s%s".ptr,
                    is_y2_var(&gi, i) ? "x1y2".ptr : "x1y1".ptr,
                    plotname(dset, gi.list[i], 1),
                    is_y2_var(&gi, i) ? _("right").ptr : _("left").ptr,
                    withstr.ptr,
                    lwstr.ptr,
                    (i == lmax - 1)? "\n".ptr : ", \\\n".ptr);
        }
    } else if (gi.flags & GptFlags.factor) {
        /* plot shows separation by a discrete factor variable */
        int lmax = gi.list[0];
        SeriesTable * st;
        double di;
        int j, vj;

        strcpy(s2.ptr, plotname(dset, gi.fid, 1));
        st = series_get_string_table(dset, gi.fid);

        for (j=1; j<lmax; j++) {
            vj = gi.list[j];
            set_plot_withstr(&gi, j, withstr.ptr);
            strcpy(s1.ptr, plotname(dset, vj, 1));
            if (!y_var_is_factorized(&gi, vj)) {
                fprintf(fp, " '-' using 1:2 title \"%s\" %s".ptr, s1.ptr, withstr.ptr);
                if (j < lmax - 1) {
                    fputs(", \\\n", fp);
                } else {
                    fputc('\n', fp);
                }
                continue;
            }
            for (i=0; i<gi.n_fvals; i++) {
                di = gretl_vector_get(gi.fvals, i);
                if (st != null) {
                    fprintf(fp, " '-' using 1:2 title \"%s\" %s".ptr,
                            series_table_get_string(st, di), withstr.ptr);
                } else {
                    fprintf(fp, " '-' using 1:2 title \"%s,%s=%g\" %s".ptr,
                            s1.ptr, s2.ptr, di, withstr.ptr);
                }
                if (i < gi.n_fvals - 1 || j < lmax - 1) {
                    fputs(", \\\n", fp);
                } else {
                    fputc('\n', fp);
                }
            }
        }
    } else if (gi.yformula != null) {
        /* we have a formula to plot, not just data */
        fprintf(fp, " '-' using 1:2 title \"%s\" w points, \\\n".ptr, _("actual").ptr);
        fprintf(fp, "%s title '%s' w lines\n".ptr, gi.yformula, _("fitted").ptr);
    } else if (gi.flags & GptFlags.fa) {
        /* this is a fitted vs actual plot */
        int tmp = gi.list[1];

        gi.list[1] = gi.list[2];
        gi.list[2] = tmp;
        set_plot_withstr(&gi, 1, withstr.ptr);
        fprintf(fp, " '-' using 1:2 title \"%s\" %s, \\\n".ptr, _("actual").ptr, withstr.ptr);
        fprintf(fp, " '-' using 1:2 title \"%s\" %s\n".ptr, _("fitted").ptr, withstr.ptr);
    } else {
        /* all other cases */
        int lmax = gi.list[0] - 1;

        for (i=1; i<=lmax; i++)  {
            set_lwstr(lwstr.ptr);
            if (gi.list[0] == 2 && !(gi.flags & GptFlags.timefmt)) {
                *s1.ptr = '\0';
            } else {
                strcpy(s1.ptr, plotname(dset, gi.list[i], 1));
            }
            set_plot_withstr(&gi, i, withstr.ptr);
            if (string_is_blank(s1.ptr)) {
                fprintf(fp, " '-' using 1:2 notitle %s%s".ptr, withstr.ptr, lwstr.ptr);
            } else {
                fprintf(fp, " '-' using 1:2 title \"%s\" %s%s".ptr, s1.ptr, withstr.ptr, lwstr.ptr);
            }
            if (i < lmax || (gi.flags & GptFlags.auto_fit)) {
                fputs(", \\\n", fp);
            } else {
                fputc('\n', fp);
            }
        }
    }

    if (fitline != null) {
        fputs(fitline, fp);
    }

    /* print the data to be graphed */
    if (gi.flags & GptFlags.factor) {
        print_gp_factorized_data(&gi, dset, fp);
    } else {
        print_gp_data(&gi, dset, fp);
    }

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    g_free(fitline);
    clear_gpinfo(&gi);

    return err;
}

int theil_forecast_plot(const(int) * plotlist, const(DatasetRaw) * dset)
{
    FILE * fp = null;
    GnuplotInfo gi;
    int vx, vy;
    int err = 0;

    gretl_error_clear();

    if (plotlist[0] != 2) {
        return GretlError.data;
    }

    /* OPT_S: don't show auto-fitted line */
    err = gpinfo_init(&gi, Gretlopt.OPT_S, plotlist, null, dset);
    if (err) {
        goto bailout;
    }

    /* ensure the time-series flag is unset */
    gi.flags &= ~GptFlags.ts;

    err = graph_list_adjust_sample(&gi, dset, 1);
    if (err) {
        goto bailout;
    }

    fp = open_plot_input_file(PlotType.regular, gi.flags, &err);
    if (err) {
        goto bailout;
    }

    vx = gi.list[2];
    vy = gi.list[1];

    print_axis_label('x', plotname(dset, vx, 1), fp);
    print_axis_label('y', plotname(dset, vy, 1), fp);

    fputs("set xzeroaxis\n", fp);
    gnuplot_missval_string(fp);
    fputs("set key left top\n", fp);

    gretl_push_c_numeric_locale();

    print_x_range_from_list(&gi, dset, fp);

    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    fprintf(fp, " x title \"%s\" w lines\n".ptr, _("actual = predicted").ptr);

    print_gp_data(&gi, dset, fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);

 bailout:

    clear_gpinfo(&gi);

    return err;
}

void multi_time_tics(const(double) * obs, const(DatasetRaw) * dset, FILE * fp)
{
    double startdate = obs[dset.t1];
    double enddate = obs[dset.t2];
    double obsrange = enddate - startdate;
    int k1 = cast(int) ceil(startdate);
    int k2 = cast(int) floor(enddate);
    int nmaj = k2 - k1 + 1;

    fputs("set xtics nomirror\n", fp);

    if (obsrange > 8) {
        double incr = obsrange / 4;

        fprintf(fp, "set xrange [%g:%g]\n".ptr, floor(startdate), ceil(enddate));
        fprintf(fp, "set xtics %g,%g\n".ptr, ceil(startdate), floor(incr));
    } else if (nmaj == 0) {
        fputs("set format x ''\n", fp);
    } else {
        /* integer major tics plus minor */
        int T = dset.t2 - dset.t1 + 1;

        fprintf(fp, "set xrange [%g:%g]\n".ptr, startdate, enddate);
        fprintf(fp, "set xtics %g,1\n".ptr, floor(startdate));
        if (T < 55) {
            fprintf(fp, "set mxtics %d\n".ptr, dset.pd);
        }
    }
}

void multi_set_timefmt(const(DatasetRaw) * dset, const(double) * obs, FILE * fp)
{
    double T = obs[dset.t2] - obs[dset.t1];
    int ntics = 6;

    fputs("set xdata time\n", fp);
    fputs("set timefmt \"%s\"\n", fp);
    if (single_year_sample(dset, dset.t1, dset.t2)) {
        fputs("set format x \"%m-%d\"\n", fp);
    } else {
        fputs("set format x \"%Y-%m-%d\"\n", fp);
    }
    fputs("set xtics rotate by -45\n", fp);
    fprintf(fp, "set xtics %g\n".ptr, round(T/ntics));
}

/**
 * multi_plots:
 * @list: list of variables to plot, by ID number.
 * @dset: dataset struct.
 * @opt: can include %OPT_O to use lines, %PT_T to
 * produce time-series plots, %OPT_U to direct output
 * to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the variables in @list, and calls gnuplot to make
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int multi_plots(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    GptFlags flags = cast(GptFlags) 0;
    int xvar = 0, yvar = 0;
    const(double) * x = null;
    const(double) * y = null;
    const(double) * obs = null;
    int rows, cols, tseries = 0;
    int use_timefmt = 0;
    int tsdata = 0;
    int * plotlist = null;
    int pos, nplots = 0;
    FILE * fp = null;
    int i, t, err = 0;

    pos = gretl_list_separator_position(list);
    tsdata = dataset_is_time_series(dset);

    if (opt & Gretlopt.OPT_O) {
        flags |= GptFlags.lines;
    }

    if (opt & Gretlopt.OPT_T) {
        /* the "tsplots" command */
        if (pos > 0) {
            /* separator not accepted */
            err = GretlError.invarg;
        } else if (!tsdata) {
            err = GretlError.pdwrong;
        } else {
            flags |= GptFlags.lines;
        }
    }

    if (pos == 0) {
        /* plot against time or index */
        plotlist = gretl_list_copy(list);
        flags |= GptFlags.lines;
        if (tsdata) {
            tseries = 1;
            if (calendar_data(dset)) {
                use_timefmt = 1;
            }
        }
        obs = gretl_plotx(dset, use_timefmt ? Gretlopt.OPT_T : Gretlopt.OPT_S);
        if (obs == null) {
            return GretlError.alloc;
        }
    } else if (pos > 2) {
        /* plot several yvars against one xvar */
        plotlist = gretl_list_new(pos - 1);
        xvar = list[list[0]];
        x = dset.Z[xvar];
    } else {
        /* plot one yvar against several xvars */
        plotlist = gretl_list_new(list[0] - pos);
        yvar = list[1];
        y = dset.Z[yvar];
    }

    if (plotlist == null) {
        return GretlError.alloc;
    }

    if (yvar) {
        for (i=1; i<=plotlist[0]; i++) {
            plotlist[i] = list[i + pos];
        }
    } else if (xvar) {
        for (i=1; i<pos; i++) {
            plotlist[i] = list[i];
        }
    }

    /* max 16 plots */
    if (plotlist[0] > 16) {
        plotlist[0] = 16;
    }

    nplots = plotlist[0];

    if (nplots > 1) {
        get_multiplot_layout(nplots, tseries, &rows, &cols);
        if (use_timefmt) {
            gp_small_font_size = nplots > 2 ? 6 : 0;
        } else {
            maybe_set_small_font(nplots);
        }
        if (nplots > 12) {
            flags |= GptFlags.xxl;
        } else if (nplots > 9) {
            flags |= GptFlags.xl;
        }
    }

    fp = open_plot_input_file(PlotType.multi_basic, flags, &err);
    if (err) {
        return err;
    }

    if (nplots > 1) {
        fprintf(fp, "set multiplot layout %d,%d\n".ptr, rows, cols);
    }
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (use_timefmt) {
        fprintf(fp, "set xrange [%.0f:%.0f]\n".ptr, obs[dset.t1], obs[dset.t2]);
        multi_set_timefmt(dset, obs, fp);
    } else if (obs != null) {
        multi_time_tics(obs, dset, fp);
    } else {
        /* avoid having points sticking to the axes */
        fputs("set offsets graph 0.02, graph 0.02, graph 0.02, graph 0.02\n", fp);
        fputs("set noxtics\nset noytics\n", fp);
    }

    fputs("set xzeroaxis\n", fp);

    for (i=0; i<nplots; i++) {
        int j = plotlist[i+1];

        if (obs != null) {
            fputs("set noxlabel\n", fp);
            fputs("set noylabel\n", fp);
            fprintf(fp, "set title '%s'\n".ptr, plotname(dset, j, 1));
        } else {
            fprintf(fp, "set xlabel '%s'\n".ptr,
                    (yvar)? dset.varname[j] :
                    dset.varname[xvar]);
            fprintf(fp, "set ylabel '%s'\n".ptr,
                    (yvar)? dset.varname[yvar] :
                    dset.varname[j]);
        }
        fputs("plot '-' using 1:2", fp);
        if (flags & GptFlags.lines) {
            fputs(" with lines", fp);
        }
        fputc('\n', fp);

        for (t=dset.t1; t<=dset.t2; t++) {
            double xt = yvar ? dset.Z[j][t] : xvar ? x[t] : obs[t];
            double yt = yvar ? y[t] : dset.Z[j][t];

            write_gp_dataval(xt, fp, 0);
            write_gp_dataval(yt, fp, 1);
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    if (nplots > 1) {
        fputs("unset multiplot\n", fp);
    }

    free(plotlist);

    return finalize_plot_input_file(fp);
}

int matrix_plotx_ok(const(GretlMatrixRaw) * m, const(DatasetRaw) * dset, int * pt1, int * pt2, int * ppd)
{
    if (dset == null) {
        return 0;
    } else if (m.rows == dset.n) {
        return 1;
    } else {
        int t1 = gretl_matrix_get_t1(m);
        int t2 = gretl_matrix_get_t2(m);

        if (t2 > t1 && t2 < dset.n) {
            *pt1 = t1;
            *pt2 = t2;
            *ppd = dset.pd;
            return 1;
        }
    }

    return 0;
}

const(double) * matrix_col(const(GretlMatrixRaw) * m, int j)
{
    const(double) * x = m.val;

    return x + (j-1) * m.rows;
}

void plot_colname(char * s, const(char) ** colnames, int j)
{
    if (colnames != null) {
        const(char) * name = colnames[j-1];

        *s = '\0';
        if (strlen(name) >= 16) {
            strncat(s, name, 14);
            strcat(s, "~");
        } else {
            strncat(s, name, 15);
        }
    } else {
        sprintf(s, "col %d".ptr, j);
    }
}

double get_obsx(const(double) * obs, int t, int s)
{
    return (obs != null)? obs[t] : s;
}

/**
 * matrix_multi_plots:
 * @m: matrix containing data to plot.
 * @list: list of columns to plot, or NULL.
 * @dset: dataset pointer, or NULL.
 * @opt: can include %OPT_O to use lines, %OPT_U to
 * direct output to a named file.
 *
 * Writes a gnuplot plot file to display up to 16 small graphs
 * based on the data in @m, and calls gnuplot to make
 * the graph.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int matrix_multi_plots(const(GretlMatrixRaw) * m, const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    GptFlags flags = cast(GptFlags) 0;
    const(double) * x = null;
    const(double) * y = null;
    const(double) * obs = null;
    const(char) ** colnames = null;
    FILE * fp = null;
    int * plotlist = null;
    int need_list = 0;
    int rows, cols;
    int xcol = 0, ycol = 0;
    int t1 = 0, t2 = 0, pd = 1;
    int pos = 0, nplots = 0;
    int simple_obs = 0;
    int i, t, err = 0;

    if (gretl_is_null_matrix(m)) {
        return GretlError.data;
    }

    if (opt & Gretlopt.OPT_O) {
        flags |= GptFlags.lines;
    }

    if (list != null) {
        for (i=1; i<=list[0]; i++) {
            if (list[i] == LISTSEP) {
                pos = i;
            } else if (list[i] < 1 || list[i] > m.cols) {
                err = GretlError.invarg;
                break;
            }
        }
    }

    if (err) {
        return err;
    }

    t1 = 0;
    t2 = m.rows - 1;

    if (pos == 0) {
        /* plot against time or index */
        if (matrix_plotx_ok(m, dset, &t1, &t2, &pd)) {
            obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
            if (obs == null) {
                return GretlError.alloc;
            }
        } else {
            simple_obs = 1;
        }
        if (list != null && list[0] > 0) {
            need_list = 1;
            plotlist = gretl_list_copy(list);
        }
        flags |= GptFlags.lines;
    } else if (pos > 2) {
        /* plot several yvars against one xvar */
        need_list = 1;
        plotlist = gretl_list_new(pos - 1);
        xcol = list[list[0]];
        x = matrix_col(m, xcol);
    } else {
        /* plot one yvar against several xvars */
        need_list = 1;
        plotlist = gretl_list_new(list[0] - pos);
        ycol = list[1];
        y = matrix_col(m, ycol);
    }

    if (need_list && plotlist == null) {
        return GretlError.alloc;
    }

    if (plotlist != null) {
        if (y != null) {
            for (i=1; i<=plotlist[0]; i++) {
                plotlist[i] = list[i + pos];
            }
        } else if (x != null) {
            for (i=1; i<pos; i++) {
                plotlist[i] = list[i];
            }
        }
        /* max 16 plots */
        if (plotlist[0] > 16) {
            plotlist[0] = 16;
        }
        nplots = plotlist[0];
    } else {
        nplots = (m.cols > 16)? 16 : m.cols;
    }

    get_multiplot_layout(nplots, 0, &rows, &cols);
    maybe_set_small_font(nplots);

    if (nplots > 12) {
        flags |= GptFlags.xxl;
    } else if (nplots > 9) {
        flags |= GptFlags.xl;
    }

    fp = open_plot_input_file(PlotType.multi_basic, flags, &err);
    if (err) {
        return err;
    }

    colnames = gretl_matrix_get_colnames(m);

    fprintf(fp, "set multiplot layout %d,%d\n".ptr, rows, cols);
    fputs("set xzeroaxis\n", fp);
    fputs("set nokey\n", fp);

    gretl_push_c_numeric_locale();

    if (obs != null) {
        double startdate = obs[t1];
        double enddate = obs[t2];
        int incr, T = t2 - t1 + 1;

        fprintf(fp, "set xrange [%g:%g]\n".ptr, floor(startdate), ceil(enddate));

        incr = (pd == 1)? (T / 6) : (T / (4 * pd));
        if (incr > 0) {
            fprintf(fp, "set xtics %g, %d\n".ptr, ceil(startdate), incr);
        }
    } else if (simple_obs) {
        int incr = m.rows / 6;

        fprintf(fp, "set xrange [0:%d]\n".ptr, m.rows - 1);
        if (incr > 0) {
            fprintf(fp, "set xtics 0, %d\n".ptr, incr);
        }
    } else {
        fputs("set noxtics\nset noytics\n", fp);
    }

    if (obs != null || simple_obs) {
        fputs("set noxlabel\n", fp);
        fputs("set noylabel\n", fp);
    }

    for (i=0; i<nplots; i++) {
        int j = (plotlist == null)? (i+1) : plotlist[i+1];
        const(double) * zj = matrix_col(m, j);
        char[16] label;

        if (obs != null || simple_obs) {
            plot_colname(label.ptr, colnames, j);
            fprintf(fp, "set title '%s'\n".ptr, label.ptr);
        } else {
            plot_colname(label.ptr, colnames, (y != null)? j : xcol);
            fprintf(fp, "set xlabel '%s'\n".ptr, label.ptr);
            plot_colname(label.ptr, colnames, (y != null)? ycol : j);
            fprintf(fp, "set ylabel '%s'\n".ptr, label.ptr);
        }

        fputs("plot '-' using 1:2", fp);
        if (flags & GptFlags.lines) {
            fputs(" with lines", fp);
        }
        fputc('\n', fp);

        for (t=t1; t<=t2; t++) {
            int s = t - t1;
            double xt = ycol ? zj[s] : xcol ? x[s] : get_obsx(obs, t, s);
            double yt = (y != null)? y[s] : zj[s];

            write_gp_dataval(xt, fp, 0);
            write_gp_dataval(yt, fp, 1);
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

    free(plotlist);

    return finalize_plot_input_file(fp);
}

FILE * get_3d_input_file(int * err)
{
    FILE * fp = null;
    char[MAXLEN] fname;

    sprintf(fname.ptr, "%sgpttmp.plt".ptr, gretl_dotdir());
    fp = gretl_fopen(fname.ptr, "w");

    if (fp == null) {
        *err = GretlError.fopen;
    } else {
        gretl_set_path_by_name("plotfile", fname.ptr);
    }

    return fp;
}

char * maybe_get_surface(const(int) * list, DatasetRaw * dset, int force_show)
{
    Model smod;
    int[5] olslist;
    int show = force_show;
    char * ret = null;

    olslist[0] = 4;
    olslist[1] = list[3];
    olslist[2] = 0;
    olslist[3] = list[2];
    olslist[4] = list[1];

    /* get OLS estimates */
    smod = lsq(olslist.ptr, dset, GretlCmdIndex.ols, Gretlopt.OPT_A);

    if (smod.errcode) {
        show = 0;
    } else if (!show) {
        /* show fit if "good enough" */
        show = !na(smod.fstt) &&
            snedecor_cdf_comp(smod.dfn, smod.dfd, smod.fstt) < 0.10;
    }

    if (show) {
        /* formulate parametric plot line */
        double adj, umin, umax, vmin, vmax;
        double * b = smod.coeff;

        gretl_minmax(dset.t1, dset.t2, dset.Z[list[2]], &umin, &umax);
        gretl_minmax(dset.t1, dset.t2, dset.Z[list[1]], &vmin, &vmax);
        adj = (umax - umin) * 0.02;
        umin -= adj; umax += adj;
        adj = (vmax - vmin) * 0.02;
        vmin -= adj; vmax += adj;

        ret = g_strdup_printf("[u=%g:%g] [v=%g:%g] %g+(%g)*u+(%g)*v notitle".ptr,
                              umin, umax, vmin, vmax,
                              b[0], b[1], b[2]);
    }

    clear_model(&smod);

    return ret;
}

/**
 * gnuplot_3d:
 * @list: list of variables to plot, by ID number: Y, X, Z
 * @literal: literal command(s) to pass to gnuplot (or NULL)
 * @dset: pointer to dataset.
 * @show_surface: show surface unconditionally.
 * @interactive: on input, whether an interactive 3-D plot
 * is wanted; on output, whether it's actually implemented.
 *
 * Writes a gnuplot plot file to display a 3-D plot (Z on
 * the vertical axis, X and Y on base plane).
 *
 * Returns: 0 on successful completion, error code on error.
 */

int gnuplot_3d(int * list, const(char) * literal, DatasetRaw * dset, int show_surface, int * interactive)
{
    FILE * fp = null;
    int t1 = dset.t1;
    int t2 = dset.t2;
    int save_t1, save_t2;
    int lo = list[0];
    int[4] datlist;
    const(char) * term = null;
    char * surface = null;
    int t, err = 0;

    if (lo != 3) {
        fprintf(stderr, "gnuplot_3d needs three variables (only)\n".ptr);
        return GretlError.data;
    }

    list_adjust_sample(list, &t1, &t2, dset, null);

    /* if resulting sample range is empty, complain */
    if (t1 >= t2) {
        return GretlError.missdata;
    }

    if (*interactive) {
        /* On Windows we let the gnuplot terminal default to
           "win"; on other systems we need a suitable
           terminal for interactive 3-D display.
        */
        if (gnuplot_has_wxt()) {
            term = "wxt size 640,420 noenhanced";
        } else if (gnuplot_has_x11()) {
            term = "x11";
        } else if (gnuplot_has_qt()) {
            term = "qt";
        } else {
            *interactive = 0;
        }
    }

    if (*interactive) {
        fp = get_3d_input_file(&err);
    } else {
        fp = open_plot_input_file(PlotType.plot_3d, cast(GptFlags) 0, &err);
    }

    if (err) {
        return err;
    }

    /* record and (potentially) modify sample range */
    save_t1 = dset.t1;
    save_t2 = dset.t2;
    dset.t1 = t1;
    dset.t2 = t2;

    if (*interactive) {
        if (term != null) {
            fprintf(fp, "set term %s\n".ptr, term);
        }
        write_plot_line_styles(PlotType.plot_3d, fp);
    }

    gretl_push_c_numeric_locale();

    print_axis_label('x', plotname(dset, list[2], 1), fp);
    print_axis_label('y', plotname(dset, list[1], 1), fp);
    print_axis_label('z', plotname(dset, list[3], 1), fp);

    gnuplot_missval_string(fp);

    if (literal != null) {
        print_gnuplot_literal_lines(literal, fp);
    }

    surface = maybe_get_surface(list, dset, show_surface);
    if (surface != null) {
        fprintf(fp, "splot %s, \\\n'-' notitle w p\n".ptr, surface);
        g_free(surface);
    } else {
        fputs("splot '-' notitle w p\n", fp);
    }

    datlist[0] = 3;
    datlist[1] = list[2];
    datlist[2] = list[1];
    datlist[3] = list[3];

    for (t=t1; t<=t2; t++) {
        const(char) * label = null;

        if (dset.markers) {
            label = dset.S[t];
        }
        printvars(fp, t, datlist.ptr, dset, null, label, 0.0);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    dset.t1 = save_t1;
    dset.t2 = save_t2;

    if (interactive) {
        fputs("pause mouse close\n", fp);
        fclose(fp);
    } else {
        err = finalize_plot_input_file(fp);
    }

    return err;
}

/**
 * open_3d_plot_input_file:
 * @iact: on input, non-zero if an interactive plot is
 * preferred, 0 otherwise; on output, non-zero if interactive
 * status can be supported, 0 otherwise.
 *
 * Writes a gnuplot plot file to display a 3D plot
 * (interactive if requested and feasible).
 *
 * Returns: FILE pointer on success, NULL on error.
 */

FILE * open_3d_plot_input_file(int * iact)
{
    const(char) * term = null;
    FILE * fp = null;
    int err = 0;

    if (*iact != 0) {
        /* On Windows we let the gnuplot terminal default to
           "win"; on other operating systems we need a suitable
           terminal for interactive 3-D display.
        */
        if (gnuplot_has_wxt()) {
            term = "wxt size 640,420 noenhanced";
        } else if (gnuplot_has_x11()) {
            term = "x11";
        } else if (gnuplot_has_qt()) {
            term = "qt";
        } else {
            /* can't do it? */
            *iact = 0;
        }
    }

    if (*iact != 0) {
        fp = get_3d_input_file(&err);
    } else {
        fp = open_plot_input_file(PlotType.plot_3d, cast(GptFlags) 0, &err);
    }

    if (*iact) {
        if (term != null) {
            fprintf(fp, "set term %s\n".ptr, term);
        }
        write_plot_line_styles(PlotType.plot_3d, fp);
    }

    return fp;
}

char * make_freq_test_label(int teststat, double v, double pv)
{
    char * s = null;

    gretl_pop_c_numeric_locale();
    if (teststat == TestStats.z) {
        s = g_strdup_printf("z = %.3f [%.4f]".ptr, v, pv);
    } else if (teststat == TestStats.normal_chisq) {
        s = g_strdup_printf("%s(2) = %.3f [%.4f]".ptr, _("Chi-square").ptr, v, pv);
    }
    gretl_push_c_numeric_locale();

    return s;
}

char * make_freq_dist_label(int dist, double x, double y)
{
    char * s = null;
    char c;
    char[10] test;

    gretl_pop_c_numeric_locale();
    sprintf(test.ptr, "%g".ptr, 0.5);
    c = strchr(test.ptr, ',') ? ' ' : ',';

    if (dist == DistCode.normal) {
        s = g_strdup_printf("N(%.5g%c%.5g)".ptr, x, c, y);
    } else if (dist == DistCode.gamma) {
        s = g_strdup_printf("gamma(%.5g%c%.5g)".ptr, x, c, y);
    }
    gretl_push_c_numeric_locale();

    return s;
}

/* Below: a fix for the case where the y-range is by default
   degenerate, in which case gnuplot produces a graph OK, but
   issues a warning and returns non-zero.
*/

void maybe_set_yrange(FreqDist * freq, double lambda, FILE * fp)
{
    double ymin = 1.0e+200;
    double ymax = -1.0e+200;
    int i;

    for (i=0; i<freq.numbins; i++) {
        if (freq.f[i] > ymax) {
            ymax = freq.f[i];
        }
        if (freq.f[i] < ymin) {
            ymin = freq.f[i];
        }
    }

    if (ymax == ymin) {
        fprintf(fp, "set yrange [%.10g:%.10g]\n".ptr, ymax * lambda * 0.99,
                ymax * lambda * 1.01);
    } else {
        fprintf(fp, "set yrange [0.0:%.10g]\n".ptr, ymax * lambda * 1.1);
    }
}

double discrete_minskip(FreqDist * freq)
{
    double s, ms = freq.midpt[1] - freq.midpt[0];
    int i;

    for (i=2; i<freq.numbins; i++) {
        s = freq.midpt[i] - freq.midpt[i-1];
        if (s < ms) {
            ms = s;
        }
    }

    return ms;
}

/**
 * plot_freq:
 * @freq: pointer to frequency distribution struct.
 * @dist: probability distribution code.
 *
 * Plot the actual frequency distribution for a variable versus a
 * theoretical distribution: Gaussian, gamma or none.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_freq(FreqDist * freq, DistCode dist)
{
    double alpha = 0.0, beta = 0.0, lambda = 1.0;
    FILE * fp = null;
    int i, K = freq.numbins;
    char[32] withstr = 0;
    char * label = null;
    double plotmin = 0.0, plotmax = 0.0;
    double barwidth;
    const(double) * endpt;
    int plottype, use_boxes = 1;
    char ** S = null;
    int ns = 0, real_ns = 0;
    int err = 0;

    if (K == 0) {
        return GretlError.data;
    } else if (K == 1) {
        gretl_errmsg_sprintf(_("'%s' is a constant").ptr, freq.varname.ptr);
        return GretlError.data;
    }

    if (freq.strvals) {
        dist = DistCode.none; /* just to be safe */
    }

    if (dist == DistCode.normal) {
        plottype = PlotType.freq_normal;
    } else if (dist == DistCode.gamma) {
        plottype = PlotType.freq_gamma;
    } else if (freq.discrete) {
        plottype = PlotType.freq_discrete;
    } else {
        plottype = PlotType.freq_simple;
    }

    fp = open_plot_input_file(cast(PlotType) plottype, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    if (freq.strvals) {
        endpt = null;
        barwidth = 1;
        use_boxes = 0;
    } else if (freq.discrete) {
        endpt = freq.midpt;
        barwidth = discrete_minskip(freq);
        use_boxes = 0;
    } else {
        /* equally sized bins, width to be determined */
        endpt = freq.endpt;
        barwidth = freq.endpt[K-1] - freq.endpt[K-2];
    }

    gretl_push_c_numeric_locale();

    if (dist) {
        int nlit = 2 + 2 * (!na(freq.test)) + real_ns;

        lambda = 1.0 / (freq.n * barwidth);

        if (dist == DistCode.normal) {
            fprintf(fp, "# literal lines = %d\n".ptr, nlit);
            fprintf(fp, "sigma = %g\n".ptr, freq.sdx);
            fprintf(fp, "mu = %g\n".ptr, freq.xbar);

            plotmin = endpt[0] - barwidth;
            if (plotmin > freq.xbar - 3.3 * freq.sdx) {
                plotmin = freq.xbar - 3.3 * freq.sdx;
            }
            plotmax = endpt[K-1] + barwidth;
            if (plotmax < freq.xbar + 3.3 * freq.sdx) {
                plotmax = freq.xbar + 3.3 * freq.sdx;
            }
            if (!na(freq.test)) {
                fprintf(fp, "set label \"%s:\" at graph .03, graph .97 front\n".ptr,
                        _("Test statistic for normality").ptr);
                label = make_freq_test_label(TestStats.normal_chisq, freq.test,
                                             chisq_cdf_comp(2, freq.test));
                fprintf(fp, "set label '%s' at graph .03, graph .93 front\n".ptr,
                        label);
                g_free(label);
            }
            if (real_ns > 0) {
                print_extra_literal_lines(S, ns, fp);
            }
        } else if (dist == DistCode.gamma) {
            double var = freq.sdx * freq.sdx;

            /* scale param = variance/mean */
            beta = var / freq.xbar;
            /* shape param = mean/scale */
            alpha = freq.xbar / beta;

            fprintf(fp, "# literal lines = %d\n".ptr, nlit);
            fprintf(fp, "beta = %g\n".ptr, beta);
            fprintf(fp, "alpha = %g\n".ptr, alpha);
            plotmin = 0.0;
            plotmax = freq.xbar + 4.0 * freq.sdx;

            if (!na(freq.test)) {
                fprintf(fp, "set label '%s:' at graph .03, graph .97 front\n".ptr,
                        _("Test statistic for gamma").ptr);
                label = make_freq_test_label(TestStats.z, freq.test,
                                             normal_pvalue_2(freq.test));
                fprintf(fp, "set label '%s' at graph .03, graph .93 front\n".ptr,
                        label);
                g_free(label);
            }
            if (real_ns > 0) {
                print_extra_literal_lines(S, ns, fp);
            }
        }

        /* adjust min, max if needed */
        if (freq.midpt[0] < plotmin) {
            plotmin = freq.midpt[0];
        }
        if (freq.midpt[K-1] > plotmax) {
            plotmax = freq.midpt[K-1];
        }

        fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, plotmin, plotmax);
        fputs("set key right top\n", fp);
    } else {
        /* plain frequency plot (no theoretical distribution shown) */
        lambda = 1.0 / freq.n;
        if (freq.strvals) {
            plotmin = 0.5;
            plotmax = K + 0.5;
        } else {
            plotmin = freq.midpt[0] - barwidth;
            plotmax = freq.midpt[K-1] + barwidth;
        }
        fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, plotmin, plotmax);
        maybe_set_yrange(freq, lambda, fp);
        fputs("set nokey\n", fp);

        if (real_ns > 0) {
            fprintf(fp, "# literal lines = %d\n".ptr, real_ns);
            print_extra_literal_lines(S, ns, fp);
        }
    }

    if (isnan(lambda)) {
        if (fp != null) {
            fclose(fp);
        }
        return 1;
    }

    if (freq.strvals) {
        fprintf(fp, "set title '%s: relative frequencies'\n".ptr,
                freq.varname.ptr);
    } else {
        fprintf(fp, "set xlabel '%s'\n".ptr, freq.gname.ptr);
        if (dist) {
            fprintf(fp, "set ylabel '%s'\n".ptr, _("Density").ptr);
        } else {
            fprintf(fp, "set ylabel '%s'\n".ptr, _("Relative frequency").ptr);
        }
    }

    if (freq.strvals) {
        fputs("set xtics rotate by -45\n", fp);
        fputs("set xtics (", fp);
        for (i=0; i<K; i++) {
            label = g_strdup(freq.S[i]);
            gretl_utf8_truncate(label, 6);
            fprintf(fp, "\"%s\" %d".ptr, label, i+1);
            if (i < K-1) {
                fputs(", ", fp);
            }
            g_free(label);
        }
        fputs(")\n", fp);
    } else if (freq.discrete > 1 && K < 10 && fabs(freq.midpt[K-1]) < 11) {
        /* few values, all integers, none too big: force integer tic marks */
        fprintf(fp, "set xtics %.0f, 1, %.0f\n".ptr, freq.midpt[0],
                freq.midpt[K-1]);
    }

    /* plot instructions */
    if (use_boxes) {
        fputs("set style fill solid 0.6\n", fp);
        strcpy(withstr.ptr, "w boxes");
    } else {
        strcpy(withstr.ptr, "w impulses lw 3");
    }

    if (!dist) {
        fprintf(fp, "plot '-' using 1:2 %s\n".ptr, withstr.ptr);
    } else if (dist == DistCode.normal) {
        label = make_freq_dist_label(dist, freq.xbar, freq.sdx);
        fputs("plot \\\n", fp);
        fprintf(fp, "'-' using 1:2 title \"%s\" %s, \\\n1.0/(sqrt(2.0*pi)*sigma)*exp(-.5*((x-mu)/sigma)**2) title \"%s\" w lines\n".ptr,
                _("relative frequency").ptr, withstr.ptr, label);
        g_free(label);
    } else if (dist == DistCode.gamma) {
        label = make_freq_dist_label(dist, alpha, beta);
        fputs("plot \\\n", fp);
        fprintf(fp, "'-' using 1:2 title '%s' %s, \\\nx**(alpha-1.0)*exp(-x/beta)/(exp(lgamma(alpha))*(beta**alpha)) title \"%s\" w lines\n".ptr,
                _("relative frequency").ptr, withstr.ptr, label);
        g_free(label);
    }

    for (i=0; i<K; i++) {
        if (freq.midpt == null) {
            fprintf(fp, "%d %.10g\n".ptr, i + 1, lambda * freq.f[i]);
        } else {
            fprintf(fp, "%.10g %.10g\n".ptr, freq.midpt[i], lambda * freq.f[i]);
        }
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* Return a string which doubles any singleton %-signs in @s,
   so that gnuplot doesn't choke when @s is used in a tics
   expression.
*/

char * fix_percent(const(char) * s)
{
    char * ret = null;
    const(char) * p = s;
    int len = cast(int) (cast(int) strlen(s));
    int n = 0;

    while (*p) {
    if (*p == '%' && *(p+1) != '%') {
        n++;
    }
    p++;
    }

    if (n > 0) {
    int i = 0;

    ret = cast(typeof(ret)) calloc(len + n + 1, 1);
    while (*s) {
        ret[i++] = *s;
        if (*s == '%' && *(s+1) != '%') {
        ret[i++] = '%';
        }
        s++;
    }
    } else {
    ret = gretl_strdup(s);
    }

    return ret;
}

/**
 * plot_corrmat:
 * @corr: pointer to correlation matrix struct.
 * @opt: can use OPT_T for triangular representation.
 *
 * Produces a heatmap plot based on a correlation matrix.
 *
 * Returns: 0 on successful completion, error code on error.
 */

int plot_corrmat(VMatrix * corr, const(DatasetRaw) * dset, Gretlopt opt)
{
    FILE * fp;
    const(char) * pname;
    char * tmp;
    double rcrit = 0.0;
    int i, j, df, n, idx;
    int allpos = 1;
    int err = 0;

    fp = open_plot_input_file(PlotType.heatmap, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    n = corr.dim;

    /* are all the correlations non-negative? */
    for (i=0; i<n; i++) {
        for (j=i+1; j<n; j++) {
            idx = ijton(i, j, n);
            if (corr.vec[idx] < 0) {
                allpos = 0;
                break;
            }
        }
    }

    df = corr.ncrit - 2;
    if (df > 1) {
        /* determine 20% critical value */
        double tc = student_critval(df, 0.10);
        double t2 = tc * tc;

        rcrit = sqrt(t2 / (t2 + df));
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set title '%s'\n".ptr, _("Correlation matrix").ptr);
    fputs("set nokey\n", fp);
    fputs("set tics nomirror\n", fp);

    if (allpos) {
        fputs("set cbrange [0:1]\n", fp);
        if (rcrit > 0) {
            fprintf(fp, "set palette defined (0 'white', %.4f 'white', 1 'red')\n".ptr,
                    rcrit);
        } else {
            fputs("set palette defined (0 'white', 1 'red')\n", fp);
        }
    } else {
        fputs("set cbrange [-1:1]\n", fp);
        if (rcrit > 0) {
            fprintf(fp, "set palette defined (-1 'blue', %.4f 'white', %.4f 'white', 1 'red')\n".ptr,
                    -rcrit, rcrit);
        } else {
            fputs("set palette defined (-1 'blue', 0 'white', 1 'red')\n", fp);
        }
    }

    if (opt & Gretlopt.OPT_T) {
        fputs("set border 3\n", fp);
    }

    /* for grid lines */
    fputs("set x2tics 1 format '' scale 0,0.001\n", fp);
    fputs("set y2tics 1 format '' scale 0,0.001\n", fp);
    fputs("set mx2tics 2\n", fp);
    fputs("set my2tics 2\n", fp);

    /* y-axis tics */
    fputs("set ytics (", fp);
    for (i=0; i<n; i++) {
    pname = plotname(dset, corr.list[i+1], 1);
    if (strchr(pname, '%')) {
        tmp = fix_percent(pname);
        fprintf(fp, "\"%s\" %d".ptr, tmp, n-i-1);
        free(tmp);
    } else {
        fprintf(fp, "\"%s\" %d".ptr, pname, n-i-1);
    }
        if (i < n - 1) {
            fputs(", ", fp);
        }
    }
    fputs(") out\n", fp);

    /* x-axis tics */
    fputs("set xtics (", fp);
    for (i=0; i<n; i++) {
    pname = plotname(dset, corr.list[i+1], 1);
    if (strchr(pname, '%')) {
        tmp = fix_percent(pname);
        fprintf(fp, "\"%s\" %d".ptr, tmp, i);
        free(tmp);
    } else {
        fprintf(fp, "\"%s\" %d".ptr, pname, i);
    }
        if (i < n - 1) {
            fputs(", ", fp);
        }
    }
    fputs(") out\n", fp);
    fputs("set xtics rotate by 45 right\n", fp);

    /* note: "set link" requires gnuplot 5 */
    fputs("set autoscale fix\n", fp);
    fputs("set link x\n", fp);
    fputs("set link y\n", fp);
    fputs("set grid front mx2tics my2tics lw 2 lt -1 lc rgb 'white'\n", fp);

    gnuplot_missval_string(fp);
    fprintf(fp, "printcorr = %d\n".ptr, n <= 16 ? 1 : 0);

    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if ((opt & Gretlopt.OPT_T) && j > n-i-1) {
                write_gp_dataval(NADBL, fp, 0);
            } else {
                idx = ijton(n-i-1, j, n);
                fprintf(fp, "%.4f ".ptr, corr.vec[idx]);
            }
        }
        fputc('\n', fp);
    }
    fputs("EOD\n", fp);
    fputs("# end inline data\n", fp);
    fputs("if (printcorr) {\n", fp);
    fputs("plot $data matrix with image, $data matrix using 1:2:", fp);
    if (opt & Gretlopt.OPT_T) {
        fputs("($3!=$3 ? \"\" : sprintf(\"%.1f\",$3)):($3 < -0.6 ? 0xffffff : 0x000000) with labels tc rgb var\n".ptr, fp);
    } else {
        fputs("(sprintf(\"%.1f\",$3)):($3 < -0.6 ? 0xffffff : 0x000000) with labels tc rgb var\n".ptr, fp);
    }
    fputs("} else {\n", fp);
    fputs("plot $data matrix with image\n", fp);
    fputs("}\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

/* Print the y-axis data for a forecast with errors plot. */

void fcast_print_y_data(const(double) * x, const(double) * y, int t0, int t1, int t2, const(DatasetRaw) * dset, FILE * fp)
{
    int i, t, n = t2 - t0 + 1;
    double yt;

    for (i=0; i<n; i++) {
        t = t0 + i;
        yt = t < t1 ? NADBL : y[t];
        fprintf(fp, "%.10g ".ptr, x[t]);
        write_gp_dataval(yt, fp, 1);
        if (dataset_is_panel(dset) && ((t+1) % dset.pd == 0)) {
            /* hack, 2024-02-28 */
            fprintf(fp, "%g %s\n".ptr, x[t] + 0.5, GPNA.ptr);
        }
    }

    fputs("e\n", fp);
}

/* public because called from plotbands.c */

void print_user_y_data(const(double) * x, const(double) * y, int t1, int t2, FILE * fp)
{
    int t;

    for (t=t1; t<=t2; t++) {
        fprintf(fp, "%.10g ".ptr, x[t]);
        write_gp_dataval(y[t], fp, 1);
    }

    fputs("e\n", fp);
}

enum ConfBand {
    bars,
    fill,
    low,
    high
}

void print_confband_data(const(double) * x, const(FitResid) * fr, const(double) * me, int t0, int t1, int t2, const(DatasetRaw) * dset, int mode, FILE * fp)
{
    const(double) * yhat = fr.fitted;
    int fix_ci = (fr.opt & Gretlopt.OPT_X);
    int i, t, n = t2 - t0 + 1;
    double lo, hi, lyh;

    /* In the "fix_ci" case the @yhat values are the result of
       exponentiating log forecasts. To obtain the limits of the
       confidence intervals we need to take the log of @yhat, subtract
       and add the maximum error @me (which is in log terms), then
       exponentiate the result.
    */

    for (i=0; i<n; i++) {
        t = t0 + i;
        if (t < t1 || na(yhat[t]) || na(me[t])) {
            if (mode == ConfBand.low || mode == ConfBand.high) {
                fprintf(fp, "%.10g %s\n".ptr, x[t], GPNA.ptr);
            } else {
                fprintf(fp, "%.10g %s %s\n".ptr, x[t], GPNA.ptr, GPNA.ptr);
            }
        } else {
        if (fix_ci) {
        lyh = log(yhat[t]);
        lo = exp(lyh - me[t]);
        hi = exp(lyh + me[t]);
        } else {
        lo = yhat[t] - me[t];
        hi = yhat[t] + me[t];
        }
        if (mode == ConfBand.fill) {
        fprintf(fp, "%.10g %.10g %.10g\n".ptr, x[t], lo, hi);
        } else if (mode == ConfBand.low) {
        fprintf(fp, "%.10g %.10g\n".ptr, x[t], lo);
        } else if (mode == ConfBand.high) {
        fprintf(fp, "%.10g %.10g\n".ptr, x[t], hi);
        } else {
        /* ConfBand.bars */
        if (fix_ci) {
            fprintf(fp, "%.10g %.10g %.10g %.10g\n".ptr, x[t], yhat[t], lo, hi);
        } else {
            fprintf(fp, "%.10g %.10g %.10g\n".ptr, x[t], yhat[t], me[t]);
        }
        }
    }
        if (dataset_is_panel(dset) && ((t+1) % dset.pd == 0)) {
            /* hack, 2024-02-28 */
            if (mode == ConfBand.low || mode == ConfBand.high) {
                fprintf(fp, "%g %s\n".ptr, x[t] + 0.5, GPNA.ptr);
            } else {
                fprintf(fp, "%g %s %s\n".ptr, x[t] + 0.5, GPNA.ptr, GPNA.ptr);
            }
        }
    }

    fputs("e\n", fp);
}

void print_x_confband_data(const(double) * x, const(double) * y, const(double) * se, const(int) * order, double tval, int t1, int t2, int mode, FILE * fp)
{
    int i, t, n = t2 - t1 + 1;
    double et;

    for (i=0; i<n; i++) {
        t = order[i];
        if (na(y[t]) || na(se[t])) {
            if (!na(x[t])) {
                fprintf(fp, "%.10g %s\n".ptr, x[t], GPNA.ptr);
            }
        } else {
            et = tval * se[t];
            if (mode == ConfBand.low) {
                fprintf(fp, "%.10g %.10g\n".ptr, x[t], y[t] - et);
            } else if (mode == ConfBand.high) {
                fprintf(fp, "%.10g %.10g\n".ptr, x[t], y[t] + et);
            }
        }
    }

    fputs("e\n", fp);
}

struct Fsorter {
    int obs;
    double y;
}

extern (C) int compare_fs(const(void) * a, const(void) * b)
{
    const(Fsorter) * fa = cast(const(Fsorter) *) a;
    const(Fsorter) * fb = cast(const(Fsorter) *) b;

    return (fa.y > fb.y) - (fa.y < fb.y);
}

void print_filledcurve_line(const(char) * title, const(char) * rgb, int use_alpha, FILE * fp)
{
    char[10] cstr;

    if (rgb != null && *rgb != '\0') {
        *cstr.ptr = '\0';
        strncat(cstr.ptr, rgb, 9);
    } else if (use_alpha) {
        print_rgb_hash(cstr.ptr, get_alpha_shadecolor());
        /* mark this line as continuation */
        fputs(", \\\n", fp);
    } else {
        print_rgb_hash(cstr.ptr, get_shadecolor());
    }

    if (title == null) {
        fprintf(fp, "'-' using 1:2:3 notitle lc rgb \"%s\" w filledcurve".ptr,
                cstr.ptr);
    } else {
        fprintf(fp, "'-' using 1:2:3 title '%s' lc rgb \"%s\" w filledcurve".ptr,
                title, cstr.ptr);
    }
    fputs(use_alpha ? "\n" : ", \\\n", fp);
}

/* Note: if @opt includes OPT_H, that says to show fitted values for the
   pre-forecast range.
*/

int plot_fcast_errs(const(FitResid) * fr, const(double) * maxerr, const(DatasetRaw) * dset, Gretlopt opt)
{
    FILE * fp = null;
    const(double) * obs = null;
    GptFlags flags = cast(GptFlags) 0;
    PlotType ptype = PlotType.forecast;
    double xmin, xmax, xrange;
    int depvar_present = 0;
    int use_fill = 0;
    int use_lines = 0;
    int use_alpha = 0;
    int do_errs = (maxerr != null);
    int expon = (fr.opt & Gretlopt.OPT_X);
    char * cistr = null;
    int t2 = fr.t2;
    int t1, yhmin;
    int t, n, err = 0;

    /* @yhmin is the first obs at which to start plotting y-hat */
    if (do_errs) {
        t1 = fr.t0;
        yhmin = (opt & Gretlopt.OPT_H)? fr.t0 : fr.t1;
    } else {
        t1 = (fr.t0 >= 0)? fr.t0 : 0;
        yhmin = (opt & Gretlopt.OPT_H)? t1 : fr.t1;
    }

    /* don't graph empty trailing portion of forecast */
    for (t=fr.t2; t>=t1; t--) {
        if (na(fr.actual[t]) && na(fr.fitted[t])) {
            t2--;
        } else {
            break;
        }
    }

    n = t2 - t1 + 1;
    if (n < 2) {
        /* we won't draw a graph for fewer than 2 data points */
        return 1;
    }

    if (do_errs) {
        /* handle style options */
        if (opt & Gretlopt.OPT_F) {
            use_fill = 1;
        } else if (opt & Gretlopt.OPT_L) {
            use_lines = 1;
        } else if (n > 150) {
            use_fill = 1;
        }
    }

    if (use_fill) {
        ptype = PlotType.band;
    set_effective_plot_ci(GretlCmdIndex.fcast);
    }

    obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
    if (obs == null) {
    reset_effective_plot_ci();
        return GretlError.alloc;
    }

    if (dataset_is_time_series(dset)) {
        flags |= GptFlags.letterbox;
    }

    fp = open_plot_input_file(cast(PlotType) ptype, flags, &err);
    if (err) {
    reset_effective_plot_ci();
        return err;
    }

    /* check that we have any values for the actual var */
    for (t=t1; t<=t2; t++) {
        if (!na(fr.actual[t])) {
            depvar_present = 1;
            break;
        }
    }

    gretl_minmax(t1, t2, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();
    fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);
    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    if (dataset_is_time_series(dset)) {
        fprintf(fp, "# timeseries %d (letterbox)\n".ptr, dset.pd);
    }

    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("plot \\\n", fp);

    if (do_errs) {
        cistr = g_strdup_printf(_("%g percent interval").ptr, 100 * (1 - fr.alpha));
    }

    if (use_fill && !use_alpha) {
        /* plot the confidence band first so the other lines
           come out on top */
        if (do_errs) {
            print_filledcurve_line(cistr, null, 0, fp);
        }
        if (depvar_present) {
            fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1, \\\n".ptr,
                    fr.ylabel.ptr);
        }
        fprintf(fp, "'-' using 1:2 title '%s' w lines lt 2\n".ptr, _("forecast").ptr);
    } else {
        /* plot confidence bands last */
        if (depvar_present) {
            const(char) * wstr = dataset_is_time_series(dset) ? "lines" : "points";

            fprintf(fp, "'-' using 1:2 title '%s' w %s, \\\n".ptr,
                    fr.ylabel.ptr, wstr);
        }
        fprintf(fp, "'-' using 1:2 title '%s' w lines".ptr, _("forecast").ptr);
        if (do_errs) {
            if (use_fill && use_alpha) {
                print_filledcurve_line(cistr, null, 1, fp);
            } else if (use_lines) {
                fprintf(fp, ", \\\n'-' using 1:2 title '%s' w lines, \\\n".ptr,
                        cistr);
                fputs("'-' using 1:2 notitle '%s' w lines lt 3\n", fp);
            } else if (expon) {
                fprintf(fp, ", \\\n'-' using 1:2:3:4 title '%s' w errorbars\n".ptr,
                        cistr);
        } else {
                fprintf(fp, ", \\\n'-' using 1:2:3 title '%s' w errorbars\n".ptr,
                        cistr);
            }
        } else {
            fputc('\n', fp);
        }
    }

    g_free(cistr);

    gretl_push_c_numeric_locale();

    /* write out the inline data, the order depending on whether
       or not we're using fill style for the confidence bands
    */

    if (use_fill && !use_alpha) {
        if (do_errs) {
            print_confband_data(obs, fr, maxerr, t1, yhmin, t2,
                dset, ConfBand.fill, fp);
        }
        if (depvar_present) {
            fcast_print_y_data(obs, fr.actual, t1, t1, t2, dset, fp);
        }
        fcast_print_y_data(obs, fr.fitted, t1, yhmin, t2, dset, fp);
    } else {
        if (depvar_present) {
            fcast_print_y_data(obs, fr.actual, t1, t1, t2, dset, fp);
        }
        fcast_print_y_data(obs, fr.fitted, t1, yhmin, t2, dset, fp);
        if (do_errs) {
            if (use_fill && use_alpha) {
                print_confband_data(obs, fr, maxerr, t1, yhmin, t2,
                    dset, ConfBand.fill, fp);
            } else if (use_lines) {
                print_confband_data(obs, fr, maxerr, t1, yhmin, t2,
                    dset, ConfBand.low, fp);
                print_confband_data(obs, fr, maxerr, t1, yhmin, t2,
                    dset, ConfBand.high, fp);
            } else {
                print_confband_data(obs, fr, maxerr, t1, yhmin, t2,
                    dset, ConfBand.bars, fp);
            }
        }
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

gretlRGB numeric_color_from_string(const(char) * s, int * err)
{
    gretlRGB ret = 0;
    char * test = null;
    int len = cast(int) (cast(int) strlen(s));
    int offset = 0;

    if (s[0] == '#') {
        offset = 1;
    } else if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
        offset = 2;
    }

    if (offset > 0) {
        /* should be numeric [AA]RRGGBB */
        const(char) * p = s + offset;

        len -= offset;
        if (len != 6 && len != 8) {
             *err = invalid_field_error(s);
        } else {
            ret = cast(gretlRGB) strtoul(p, &test, 16);
            if (*test != '\0') {
                *err = invalid_field_error(s);
            }
        }
        return ret;
    }

    /* otherwise we should have a gnuplot color name */
    if (len < 3 || len > 17) {
        *err = invalid_field_error(s);
    } else {
        char[FILENAME_MAX] fname;
        FILE * fp;

        sprintf(fname.ptr, "%sdata%cgnuplot%cgpcolors.txt".ptr,
                gretl_home(), SLASH, SLASH);
        fp = gretl_fopen(fname.ptr, "r");

        if (fp == null) {
            *err = GretlError.fopen;
        } else {
            char[32] line;
            char[18] cname;
            int found = 0;
            guint32 u;

            while (fgets(line.ptr, line.sizeof, fp)) {
                if (sscanf(line.ptr, "%s %x".ptr, cname.ptr, &u) == 2 &&
                    strcmp(s, cname.ptr) == 0) {
                    ret = u;
                    found = 1;
                    break;
                }
            }
            fclose(fp);
            if (!found) {
                *err = invalid_field_error(s);
                ret = 0;
            }
        }
    }

    return ret;
}

int * get_x_sorted_order(const(FitResid) * fr, const(double) * x, int t1, int * pt2)
{
    int * order = null;
    Fsorter * fs;
    int nmiss = 0;
    int t2 = *pt2;
    int n = t2 - t1 + 1;
    int i, t;

    for (t=t1; t<=t2; t++) {
        if (na(fr.actual[t])) {
            nmiss++;
        }
    }

    fs = cast(typeof(fs)) malloc(n * (*fs).sizeof);
    if (fs == null) {
        return null;
    }

    order = cast(typeof(order)) malloc(n * (*order).sizeof);
    if (order == null) {
        free(fs);
        return null;
    }

    for (i=0, t=t1; t<=t2; t++, i++) {
        fs[i].obs = t;
        fs[i].y = x[t];
    }

    qsort(fs, n, (*fs).sizeof, &compare_fs);

    for (i=0; i<n; i++) {
        order[i] = fs[i].obs;
    }

    free(fs);

    if (nmiss > 0) {
        /* chop off trailing NAs */
        *pt2 = (n - nmiss) + t1 - 1;
    }

    return order;
}

void print_x_ordered_data(const(double) * x, const(double) * y, const(int) * order, int t1, int t2, FILE * fp)
{
    int i, t, n = t2 - t1 + 1;

    for (i=0; i<n; i++) {
        t = order[i];
        if (na(x[t])) {
            continue;
        } else if (na(y[t])) {
            fprintf(fp, "%.10g %s\n".ptr, x[t], GPNA.ptr);
        } else {
            fprintf(fp, "%.10g %.10g\n".ptr, x[t], y[t]);
        }
    }

    fputs("e\n", fp);
}

/* Plotting routine for a simple regression where we want to show
   actual and forecast y, plus confidence bands, against x.
   In this case We use lines to indicate the confidence bands.
*/

int plot_simple_fcast_bands(const(Model) * pmod, const(FitResid) * fr, const(DatasetRaw) * dset, Gretlopt opt)
{
    FILE * fp = null;
    const(double) * x = null;
    int * order = null;
    GptFlags flags = cast(GptFlags) 0;
    double a, xmin, xmax, xrange, tval;
    char * cistr;
    int t2 = fr.t2;
    int xv, t1, yhmin;
    int t, n, err = 0;

    /* @yhmin is the first obs at which to start plotting y-hat */
    t1 = fr.t0;
    yhmin = (opt & Gretlopt.OPT_H)? fr.t0 : fr.t1;

    /* don't graph empty trailing portion of forecast */
    for (t=fr.t2; t>=t1; t--) {
        if (na(fr.actual[t]) && na(fr.fitted[t])) {
            t2--;
        } else {
            break;
        }
    }

    n = t2 - t1 + 1;

    if (n < 2) {
        /* won't draw a graph for a single data point */
        return 1;
    }

    xv = pmod.list[pmod.list[0]];
    x = dset.Z[xv];

    order = get_x_sorted_order(fr, x, t1, &t2);
    if (order == null) {
        return GretlError.alloc;
    }

    fp = open_plot_input_file(PlotType.forecast, flags, &err);
    if (err) {
        return err;
    }

    gretl_minmax(t1, t2, x, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;

    gretl_push_c_numeric_locale();
    fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);
    gretl_pop_c_numeric_locale();

    gnuplot_missval_string(fp);

    fprintf(fp, "set xlabel '%s'\n".ptr, dset.varname[xv]);
    fprintf(fp, "set ylabel '%s'\n".ptr, fr.ylabel.ptr);

    fputs("set key left top\n", fp);
    fputs("plot \\\n", fp);

    a = 100 * (1 - fr.alpha);
    tval = student_critval(fr.df, fr.alpha / 2);

    if (opt & Gretlopt.OPT_M) {
        cistr = g_strdup_printf(_("%g%% interval for mean").ptr, a);
    } else {
        cistr = g_strdup_printf(_("%g percent interval").ptr, a);
    }

    fprintf(fp, "'-' using 1:2 title '%s' w points, \\\n".ptr, _("actual").ptr);
    fprintf(fp, "'-' using 1:2 title '%s' w lp, \\\n".ptr, _("fitted").ptr);
    fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n".ptr, cistr);
    fputs("'-' using 1:2 notitle w lines lt 3\n", fp);
    g_free(cistr);

    gretl_push_c_numeric_locale();

    print_x_ordered_data(x, fr.actual, order, t1, t2, fp);
    print_x_ordered_data(x, fr.fitted, order, yhmin, t2, fp);
    print_x_confband_data(x, fr.fitted, fr.sderr, order,
                          tval, yhmin, t2, ConfBand.low, fp);
    print_x_confband_data(x, fr.fitted, fr.sderr, order,
                          tval, yhmin, t2, ConfBand.high, fp);

    gretl_pop_c_numeric_locale();

    free(order);

    return finalize_plot_input_file(fp);
}

int plot_tau_sequence(const(Model) * pmod, const(DatasetRaw) * dset, int k)
{
    FILE * fp;
    GretlMatrixRaw * tau = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "rq_tauvec");
    GretlMatrixRaw * B = cast(GretlMatrixRaw *) gretl_model_get_data(pmod, "rq_sequence");
    double tau_i, bi, se, blo, bhi;
    double alpha, cval, tcrit, olsband;
    double[2] ymin;
    double[2] ymax;
    char * tmp;
    int ntau, bcols;
    int i, j, err = 0;

    if (tau == null || B == null) {
        return GretlError.data;
    }

    ntau = gretl_vector_get_length(tau);
    if (ntau == 0) {
        return GretlError.data;
    }

    fp = open_plot_input_file(PlotType.rq_tau, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    bcols = gretl_matrix_cols(B);

    alpha = gretl_model_get_double(pmod, "rq_alpha");
    if (na(alpha)) {
        alpha = .05;
    }

    cval = 100 * (1 - alpha);
    tcrit = student_cdf_inverse(pmod.dfd, 1 - alpha/2);
    olsband = tcrit * pmod.sderr[k];

    /* Try to figure best placement of key */

    j = k * ntau;
    if (bcols == 3) {
        blo = gretl_matrix_get(B, j, 1);
        bhi = gretl_matrix_get(B, j, 2);
    } else {
        bi = gretl_matrix_get(B, j, 0);
        se = gretl_matrix_get(B, j, 1);
        blo = bi - tcrit * se;
        bhi = bi + tcrit * se;
    }
    ymin[0] = min(blo, pmod.coeff[k] - olsband);
    ymax[0] = max(bhi, pmod.coeff[k] + olsband);

    j += ntau - 1;
    if (bcols == 3) {
        blo = gretl_matrix_get(B, j, 1);
        bhi = gretl_matrix_get(B, j, 2);
    } else {
        bi = gretl_matrix_get(B, j, 0);
        se = gretl_matrix_get(B, j, 1);
        blo = bi - tcrit * se;
        bhi = bi + tcrit * se;
    }
    ymin[1] = min(blo, pmod.coeff[k] - olsband);
    ymax[1] = max(bhi, pmod.coeff[k] + olsband);

    fputs("set xrange [0.0:1.0]\n", fp);
    fputs("set xlabel 'tau'\n", fp);

    tmp = g_strdup_printf(_("Coefficient on %s").ptr,
                          plotname(dset, pmod.list[k+2], 1));
    fprintf(fp, "set title '%s'\n".ptr, tmp);
    g_free(tmp);

    fputs("set style fill solid 0.5\n", fp);

    if (ymax[0] < .88 * ymax[1]) {
        fputs("set key left top\n", fp);
    } else if (ymax[1] < .88 * ymax[0]) {
        fputs("set key right top\n", fp);
    } else if (ymin[0] < .88 * ymin[1]) {
        fputs("set key right bottom\n", fp);
    } else {
        fputs("set key left bottom\n", fp);
    }

    fputs("plot \\\n", fp);

    /* plot the rq confidence band first so the other lines
       come out on top */
    print_filledcurve_line(null, null, 0, fp);

    /* rq estimates */
    tmp = g_strdup_printf(_("Quantile estimates with %g%% band").ptr, cval);
    fprintf(fp, "'-' using 1:2 title '%s' w lp, \\\n".ptr, tmp);
    g_free(tmp);

    /* numeric output coming up! */
    gretl_push_c_numeric_locale();

    /* ols estimate plus (1 - alpha) band */
    tmp = g_strdup_printf(_("OLS estimate with %g%% band").ptr, cval);
    fprintf(fp, "%g title '%s' w l, \\\n".ptr, pmod.coeff[k], tmp);
    g_free(tmp);
    fprintf(fp, "%g notitle w l dt 2, \\\n".ptr, pmod.coeff[k] + olsband);
    fprintf(fp, "%g notitle w l dt 2\n".ptr, pmod.coeff[k] - olsband);

    /* write out the interval values */

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
        tau_i = gretl_vector_get(tau, i);
        if (bcols == 3) {
            blo = gretl_matrix_get(B, j, 1);
            bhi = gretl_matrix_get(B, j, 2);
        } else {
            bi = gretl_matrix_get(B, j, 0);
            se = gretl_matrix_get(B, j, 1);
            blo = bi - tcrit * se;
            bhi = bi + tcrit * se;
        }
        fprintf(fp, "%.10g %.10g %.10g\n".ptr, tau_i, blo, bhi);
    }
    fputs("e\n", fp);

    for (i=0, j=k*ntau; i<ntau; i++, j++) {
        tau_i = gretl_vector_get(tau, i);
        bi = gretl_matrix_get(B, j, 0);
        fprintf(fp, "%.10g %.10g\n".ptr, tau_i, bi);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int garch_resid_plot(const(Model) * pmod, const(DatasetRaw) * dset)
{
    FILE * fp = null;
    const(double) * obs;
    const(double) * h;
    double sd2;
    int t, err = 0;

    h = cast(typeof(h)) gretl_model_get_data(pmod, "garch_h");
    if (h == null) {
        return GretlError.data;
    }

    obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
    if (obs == null) {
        return GretlError.alloc;
    }

    fp = open_plot_input_file(PlotType.garch, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    fputs("set key left top\n", fp);

    fprintf(fp, "plot \\\n'-' using 1:2 title '%s' w lines, \\\n'-' using 1:2 title '%s' w lines lt 2, \\\n'-' using 1:2 notitle w lines lt 2\n".ptr,
            _("residual").ptr, _("+- sqrt(h(t))").ptr);

    gretl_push_c_numeric_locale();

    for (t=pmod.t1; t<=pmod.t2; t++) {
        fprintf(fp, "%.10g %.10g\n".ptr, obs[t], pmod.uhat[t]);
    }
    fputs("e\n", fp);

    for (t=pmod.t1; t<=pmod.t2; t++) {
        sd2 = -sqrt(h[t]);
        fprintf(fp, "%.10g %.10g\n".ptr, obs[t], sd2);
    }
    fputs("e\n", fp);

    for (t=pmod.t1; t<=pmod.t2; t++) {
        sd2 = sqrt(h[t]);
        fprintf(fp, "%.10g %.10g\n".ptr, obs[t], sd2);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int rmplot(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    alias RangeMeanGraphFn = extern (C) int function(int, const(DatasetRaw) *, Gretlopt, Prn *);
    RangeMeanGraphFn range_mean_graph;

    range_mean_graph = cast(RangeMeanGraphFn) get_plugin_function("range_mean_graph");
    if (range_mean_graph == null) {
        return 1;
    }

    return range_mean_graph(list[1], dset, opt, prn);
}

int hurstplot(const(int) * list, DatasetRaw * dset, Gretlopt opt, Prn * prn)
{
    alias HurstExponentFn = extern (C) int function(int, const(DatasetRaw) *, Gretlopt, Prn *);
    HurstExponentFn hurst_exponent;

    hurst_exponent = cast(HurstExponentFn) get_plugin_function("hurst_exponent");
    if (hurst_exponent == null) {
        return 1;
    }

    return hurst_exponent(list[1], dset, opt, prn);
}

void get_multiplot_layout(int n, int tseries, int * rows, int * cols)
{
    if (n < 3) {
        if (tseries) {
            *cols = 1;
            *rows = 2;
        } else {
            *cols = 2;
            *rows = 1;
        }
    } else if (n < 5) {
        *cols = *rows = 2;
    } else if (n < 7) {
        *cols = 3;
        *rows = 2;
    } else if (n < 10) {
        *cols = *rows = 3;
    } else if (n < 13) {
        *cols = 4;
        *rows = 3;
    } else if (n < 17) {
        *cols = *rows = 4;
    } else {
        *cols = *rows = 0;
    }
}

int panel_ytic_width(double ymin, double ymax)
{
    char[16] s1;
    char[16] s2;
    int n1, n2;

    if (ymin < 0 && ymax > 0) {
        sprintf(s1.ptr, "% g".ptr, ymin);
        sprintf(s2.ptr, "% g".ptr, ymax);
    } else {
        sprintf(s1.ptr, "%g".ptr, ymin);
        sprintf(s2.ptr, "%g".ptr, ymax);
    }

    n1 = cast(int) (strlen(s1.ptr));
    n2 = cast(int) (strlen(s2.ptr));

    return (n1 > n2)? n1 : n2;
}

/* Panel: produce a single time-series plot showing the mean across
   groups of the series in question. Note that this function calls
   the main gnuplot() function after some panel-specific setup.
   The incoming @opt will be forwarded to gnuplot() and must not
   contain any "foreign matter".
*/

int panel_means_ts_plot(const(int) vnum, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    DatasetRaw * gset;
    int nunits, T = dset.pd;
    int[2] list = [1, 1];
    char * my_literal = null;
    char * title = null;
    int i, t, s, s0;
    int err = 0;

    gset = create_auxiliary_dataset(2, T, Gretlopt.OPT_NONE);
    if (gset == null) {
        return GretlError.alloc;
    }

    strcpy(gset.varname[1], dset.varname[vnum]);
    series_set_display_name(gset, 1, series_get_display_name(dset, vnum));

    nunits = panel_sample_size(dset);
    time_series_from_panel(gset, dset);

    s0 = dset.t1;

    for (t=0; t<T; t++) {
        double xit, xsum = 0.0;
        int n = 0;

        for (i=0; i<nunits; i++) {
            s = s0 + i * T + t;
            xit = dset.Z[vnum][s];
            if (!na(xit)) {
                xsum += xit;
                n++;
            }
        }
        gset.Z[1][t] = (n == 0)? NADBL : xsum / n;
    }

    /* ensure use of lines, time series */
    opt |= (Gretlopt.OPT_O | Gretlopt.OPT_T);

    title = g_strdup_printf(_("mean %s").ptr, plotname(dset, vnum, 1));
    my_literal = g_strdup_printf("set ylabel '%s' ; set xlabel ;".ptr,
                                 title);
    if (literal != null) {
        char * all_lit;

        all_lit = g_strdup_printf("%s %s".ptr, my_literal, literal);
        err = gnuplot(list.ptr, all_lit, gset, opt);
        g_free(all_lit);
    } else {
        err = gnuplot(list.ptr, my_literal, gset, opt);
    }

    g_free(title);
    g_free(my_literal);
    destroy_dataset(gset);

    return err;
}

int panel_means_XY_scatter(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    DatasetRaw * gset;
    int N, T = dset.pd;
    int[3] glist = [2, 1, 2];
    char * my_literal = null;
    int grpnames = 0;
    int yvar, xvar;
    int i, t, s;
    int err = 0;

    if (list == null || list[0] != 2) {
        return GretlError.data;
    }

    N = panel_sample_size(dset);

    gset = create_auxiliary_dataset(3, N, Gretlopt.OPT_NONE);
    if (gset == null) {
        return GretlError.alloc;
    }

    /* If we have valid panel group names, use them
       as obs markers here */
    grpnames = panel_group_names_ok(dset);
    if (grpnames) {
        dataset_allocate_obs_markers(gset);
    }

    yvar = list[1];
    xvar = list[2];

    strcpy(gset.varname[1], dset.varname[yvar]);
    series_set_display_name(gset, 1, series_get_display_name(dset, yvar));

    strcpy(gset.varname[2], dset.varname[xvar]);
    series_set_display_name(gset, 2, series_get_display_name(dset, xvar));

    s = dset.t1;

    for (i=0; i<N; i++) {
        double yit, ysum = 0.0;
        double xit, xsum = 0.0;
        int ny = 0, nx = 0;
        int s0 = s;

        for (t=0; t<T; t++) {
            yit = dset.Z[yvar][s];
            xit = dset.Z[xvar][s];
            if (!na(yit)) {
                ysum += yit;
                ny++;
            }
            if (!na(xit)) {
                xsum += xit;
                nx++;
            }
            s++;
        }
        gset.Z[1][i] = ny == 0 ? NADBL : ysum / ny;
        gset.Z[2][i] = nx == 0 ? NADBL : xsum / nx;
        if (gset.S != null) {
            strcpy(gset.S[i], get_panel_group_name(dset, s0));
        }
    }

    my_literal = g_strdup_printf("set title '%s';".ptr, _("Group means").ptr);

    if (literal != null) {
        char * all_lit = g_strdup_printf("%s %s".ptr, my_literal, literal);

        err = gnuplot(glist.ptr, all_lit, gset, opt);
        g_free(all_lit);
    } else {
        err = gnuplot(glist.ptr, my_literal, gset, opt);
    }

    g_free(my_literal);
    destroy_dataset(gset);

    return err;
}

void copy_string_stripped(char * targ, const(char) * src, int strip)
{
    char * tmp = gretl_strdup(src);
    size_t len = strlen(tmp) - strip;

    strcpy(targ, gretl_utf8_truncate(tmp, len));
    free(tmp);
}

/* Here we're trying to find out if the observation labels
   for a panel dataset are such that they uniquely identify
   the units/individuals (e.g. country or city names,
   repeated for each time-series observation on the given
   unit).
*/

int dataset_has_panel_labels(const(DatasetRaw) * dset, int maxlen, int * use, int * strip)
{
    int t, u, ubak = -1;
    int len, lmax = 0, fail = 0;
    int ret = 0;

    if (dset.S == null) {
        return 0;
    }

    for (t=dset.t1; t<=dset.t2 && !fail; t++) {
        u = t / dset.pd;
        if (u == ubak && strcmp(dset.S[t], dset.S[t-1])) {
            /* same unit, different label: no */
            fail = 1;
        } else if (ubak >= 0 && u != ubak &&
                   !strcmp(dset.S[t], dset.S[t-1])) {
            /* different unit, same label: no */
            fail = 2;
        }
        if (!fail) {
            len = cast(int) (strlen(dset.S[t]));
            if (len > lmax) {
                lmax = len;
            }
        }
        ubak = u;
    }

    if (!fail) {
        /* the full obs labels satisfy the criterion,
           but are they perhaps too long? */
        if (maxlen > 0 && lmax > maxlen) {
            ret = 0;
        } else {
            ret = 1;
        }
    } else if (fail == 1) {
        /* Try for a leading portion of the obs labels: for
           example we might have AUS1990, AUS1991, ...
           followed by USA1990, USA1991, ... or some such.
           We try to identify a trailing portion of the obs
           string that varies by time, and which should be
           omitted in forming "panel labels".
        */
        const(char) * s;
        int i, n, len2t, len2 = 0;
        int obslen = 0;

        fail = 0;
        for (t=dset.t1; t<=dset.t2 && !fail; t++) {
            s = dset.S[t];
            n = cast(int) (strlen(s));
            len2t = 0;
            for (i=n-1; i>0; i--) {
                if (isdigit(s[i]) || s[i] == ':' ||
                    s[i] == '-' || s[i] == '_') {
                    len2t++;
                } else {
                    break;
                }
            }
            if (len2t == 0) {
                /* no "tail" string (e.g. year) */
                fail = 1;
            } else if (len2 == 0) {
                /* starting */
                len2 = len2t;
                obslen = n;
            } else if (len2t != len2) {
                /* the "tails" don't have a common length */
                fail = 1;
            } else if (n != obslen) {
                /* the obs strings are of differing lengths */
                obslen = 0;
            }
        }

        if (!fail) {
            char[OBSLEN] s0;
            char[OBSLEN] s1;

            if (obslen > 0) {
                *use = obslen - len2;
            } else {
                *strip = len2;
            }
            /* now check that the leading portion really
               is in common for each unit/individual
            */
            *s0.ptr = '\0';
            ubak = -1;
            for (t=dset.t1; t<=dset.t2 && !fail; t++) {
                u = t / dset.pd;
                *s1.ptr = '\0';
                if (*use > 0) {
                    strncat(s1.ptr, dset.S[t], *use);
                } else {
                    copy_string_stripped(s1.ptr, dset.S[t], *strip);
                }
                if (u == ubak && strcmp(s1.ptr, s0.ptr)) {
                    /* same unit, different label: no */
                    fail = 1;
                } else if (ubak >= 0 && u != ubak && !strcmp(s1.ptr, s0.ptr)) {
                    /* different unit, same label: no */
                    fail = 2;
                }
                if (!fail && maxlen > 0 && strlen(s1.ptr) > maxlen) {
                    fail = 3;
                }
                strcpy(s0.ptr, s1.ptr);
                ubak = u;
            }
            if (fail) {
                *use = *strip = 0;
            } else {
                ret = 1;
            }
        }
    }

    /* There's a loophole above: unit m might have the same
       label as some other unit, although we've checked that
       it doesn't have the same label as unit m - 1. But
       that seems ike a corner case and I can't be bothered
       checking for it right now.
    */

    return ret;
}

/* Panel: plot one series using separate lines for each
   cross-sectional unit. The individuals' series are overlaid, in the
   same manner as a plot of several distinct time series. To do this
   we construct on the fly a notional time-series dataset.  But note:
   if it turns out the series in question is invariant across groups,
   just show a single line.

   This function calls gnuplot() after performing the required setup.
   and forwards the incoming @opt flags.
*/

int panel_overlay_ts_plot(const(int) vnum, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    DatasetRaw * gset;
    int u0, nunits, T = dset.pd;
    int * list = null;
    char * my_literal = null;
    char * title = null;
    SeriesTable * gst = null;
    const(char) * sval;
    int vg = 0;
    int panel_labels = 0;
    int maxlen = 0;
    int single_series;
    int use = 0, strip = 0;
    int i, s, s0;
    int err = 0;

    single_series = series_is_group_invariant(dset, vnum);

    if (single_series) {
        nunits = 1;
    } else {
        nunits = panel_sample_size(dset);
    }

    u0 = dset.t1 / T;

    gset = create_auxiliary_dataset(nunits + 1, T, Gretlopt.OPT_B);
    if (gset == null) {
        return GretlError.alloc;
    }

    /* add time series info to @gset */
    time_series_from_panel(gset, dset);

    list = gretl_consecutive_list_new(1, nunits);
    if (list == null) {
        destroy_dataset(gset);
        return GretlError.alloc;
    }

    if (nunits > 80) {
        /* FIXME calibrate this properly */
        maxlen = 3;
    }

    if (!single_series) {
        gst = get_panel_group_table(dset, maxlen, &vg);
        if (gst == null && dset.S != null) {
            /* maybe we have obs markers that are usable */
            panel_labels =
                dataset_has_panel_labels(dset, maxlen, &use, &strip);
        }
    }

    s0 = dset.t1;

    /* compose varnames and transcribe data */
    for (i=0; i<nunits; i++) {
        s = s0 + i * T;
        if (single_series) {
            strcpy(gset.varname[i+1], dset.varname[vnum]);
        } else if (gst != null) {
            /* look up the string for this unit/group */
            sval = series_table_get_string(gst, dset.Z[vg][s]);
            if (sval != null) {
                strncat(gset.varname[i+1], sval, VNAMELEN-1);
            } else {
                sprintf(gset.varname[i+1], "%d".ptr, u0+i+1);
            }
        } else if (panel_labels) {
            if (use > 0) {
                strncat(gset.varname[i+1], dset.S[s], use);
            } else if (strip > 0) {
                copy_string_stripped(gset.varname[i+1], dset.S[s], strip);
            } else {
                strcpy(gset.varname[i+1], dset.S[s]);
            }
        } else {
            sprintf(gset.varname[i+1], "%d".ptr, u0+i+1);
        }
        /* borrow the appropriate chunk of data */
        gset.Z[i+1] = cast(double *) dset.Z[vnum] + s;
    }

    if (nunits > 9 && T < 50) {
        opt |= Gretlopt.OPT_P; /* lines/points */
    } else {
        opt |= Gretlopt.OPT_O; /* use lines */
    }

    if (!single_series) {
        const(char) * gname = plotname(dset, vnum, 1);
        const(char) * vname = panel_group_names_varname(dset);

        if (vname != null) {
            title = g_strdup_printf(_("%s by %s").ptr, gname, vname);
        } else {
            title = g_strdup_printf(_("%s by group").ptr, gname);
        }
        my_literal = g_strdup_printf("set title '%s' ; set xlabel ;".ptr, title);
    }

    if (nunits > 80) {
        /* set file-scope global */
        xwide = 1;
    }

    if (my_literal != null && literal != null) {
        char * all_lit = g_strdup_printf("%s %s".ptr, my_literal, literal);

        err = gnuplot(list, all_lit, gset, opt | Gretlopt.OPT_T);
        g_free(all_lit);
    } else if (literal != null) {
        err = gnuplot(list, literal, gset, opt | Gretlopt.OPT_T);
    } else {
        err = gnuplot(list, my_literal, gset, opt | Gretlopt.OPT_T);
    }

    xwide = 0;
    g_free(title);
    g_free(my_literal);
    destroy_dataset(gset);
    free(list);

    return err;
}

/* Panel: plot one variable as time series, with separate plots for
   each cross-sectional unit.  By default we arrange the plots in a
   grid, but if @vstack is non-zero we make each plot full width and
   stack the plots vertically on the "page". This option is limited
   to a relatively small number of units.

   This function employs gnuplot's multiplot mechanism, and
   does not call libgretl's gnuplot() function.
*/

int panel_grid_ts_plot(int vnum, const(DatasetRaw) * dset, int vstack)
{
    FILE * fp = null;
    int w, rows, cols;
    const(double) * y, x = null;
    const(char) * vname;
    char[OBSLEN] uname;
    double xt, yt, ymin, ymax, incr;
    int u0, nunits, T = dset.pd;
    int n_ok_units = 0;
    int panel_labels = 0;
    int use = 0, strip = 0;
    int * badlist = null;
    int i, s, t, t0;
    int err = 0;

    n_ok_units = nunits = panel_sample_size(dset);
    u0 = dset.t1 / dset.pd;
    y = dset.Z[vnum];

    /* check for "blank" units */
    t0 = dset.t1;
    for (i=0; i<nunits; i++) {
        int ok = 0;

        for (t=0; t<T; t++) {
            if (!na(y[t+t0])) {
                ok = 1;
                break;
            }
        }
        if (!ok) {
            badlist = gretl_list_append_term(&badlist, i);
            n_ok_units--;
        }
        t0 += T;
    }

    if (n_ok_units < 2) {
        free(badlist);
        return GretlError.missdata;
    }

    if (vstack) {
        int xvar = plausible_panel_time_var(dset);

        if (xvar > 0) {
            x = dset.Z[xvar];
        }
        cols = 1;
        rows = n_ok_units;
    } else {
        get_multiplot_layout(n_ok_units, 0, &rows, &cols);
    }

    if (rows == 0 || cols == 0) {
        return GretlError.data;
    }

    maybe_set_small_font(nunits);

    fp = open_plot_input_file(PlotType.panel, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    if (dset.S != null) {
        panel_labels = dataset_has_panel_labels(dset, 0, &use, &strip);
    }

    vname = dset.varname[vnum];
    gretl_minmax(dset.t1, dset.t2, y, &ymin, &ymax);
    w = panel_ytic_width(ymin, ymax);

    fputs("set key left top\n", fp);
    gnuplot_missval_string(fp);
    fputs("set xtics nomirror\n", fp);
    fputs("set ytics nomirror\n", fp);
    fprintf(fp, "set format y \"%%%dg\"\n".ptr, w);
    fprintf(fp, "set multiplot layout %d,%d\n".ptr, rows, cols);

    if (vstack) {
        fputs("set noxlabel\n", fp);
    } else {
        fprintf(fp, "set xlabel '%s'\n".ptr, _("time").ptr);
    }

    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();

    t0 = dset.t1;
    for (i=0; i<nunits; i++) {
        if (in_gretl_list(badlist, i)) {
            t0 += T;
            continue;
        }
        if (panel_labels) {
            *uname.ptr = '\0';
            s = (u0 + i) * dset.pd;
            if (use > 0) {
                strncat(uname.ptr, dset.S[s], use);
            } else if (strip > 0) {
                copy_string_stripped(uname.ptr, dset.S[s], strip);
            } else {
                strcpy(uname.ptr, dset.S[s]);
            }
        } else {
            sprintf(uname.ptr, "%d".ptr, u0+i+1);
        }
        if (vstack) {
            gretl_minmax(t0, t0 + T - 1, y, &ymin, &ymax);
            incr = (ymax - ymin) / 2.0;
            fprintf(fp, "set ytics %g\n".ptr, incr);
            fprintf(fp, "set ylabel '%s (%s)'\n".ptr, vname, uname.ptr);
        } else {
            fprintf(fp, "set title '%s (%s)'\n".ptr, vname, uname.ptr);
        }

        fputs("plot \\\n'-' using 1:2 notitle w lines\n", fp);

        for (t=0; t<T; t++) {
            if (x != null) {
                xt = x[t+t0];
            } else {
                xt = t + 1;
            }
            yt = y[t+t0];
            if (na(yt)) {
                fprintf(fp, "%g %s\n".ptr, xt, GPNA.ptr);
            } else {
                fprintf(fp, "%g %.10g\n".ptr, xt, yt);
            }
        }
        fputs("e\n", fp);
        t0 += T;
    }

    gretl_pop_c_numeric_locale();

    fputs("unset multiplot\n", fp);

    free(badlist);

    return finalize_plot_input_file(fp);
}

/* This is called from the GUI only: see do_panel_plot() in
   gui/library.c. The incoming @opt flags here are as described
   below for cli_panel_plot(); they are not gnuplot options.
*/

int gretl_panel_ts_plot(int vnum, DatasetRaw * dset, Gretlopt opt)
{
    if (opt & Gretlopt.OPT_D) {
        /* small multiples in grid */
        return panel_grid_ts_plot(vnum, dset, 0);
    } else if (opt & Gretlopt.OPT_A) {
        /* small multiples in vertical stack */
        return panel_grid_ts_plot(vnum, dset, 1);
    } else if (opt & Gretlopt.OPT_M) {
        /* time-series of group means */
        return panel_means_ts_plot(vnum, null, dset, Gretlopt.OPT_NONE);
    } else {
        /* default: overlaid time series (OPT_V) */
        return panel_overlay_ts_plot(vnum, null, dset, Gretlopt.OPT_NONE);
    }
}

/* The following implements the "panplot" command.
   Option flags specific to this command are:

   OPT_M --means
   OPT_V --overlay
   OPT_S --sequence
   OPT_D --grid
   OPT_A --stack
   OPT_B --boxplots
   OPT_C --boxplot

   In addition there are some gnuplot-compatible options:

   OPT_Y --single-yaxis
   OPT_W --font
   OPT_U --output=...
   Gretlopt.OPT_B --outbuf=...
*/

int cli_panel_plot(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    Gretlopt gp_opt = Gretlopt.OPT_NONE;
    int N, vnum = list[1];
    int err;

    /* condition on multi_unit_panel_sample() ? */

    if (!dataset_is_panel(dset)) {
        gretl_errmsg_set(_("This command needs panel data").ptr);
        err = GretlError.data;
    } else {
        err = incompatible_options(opt, Gretlopt.OPT_M | Gretlopt.OPT_V | Gretlopt.OPT_S |
                                   Gretlopt.OPT_D | Gretlopt.OPT_A | Gretlopt.OPT_B | Gretlopt.OPT_C);
    }
    if (err) {
        return err;
    }

    N = panel_sample_size(dset);

    /* check for too many groups */
    if ((opt & (Gretlopt.OPT_V | Gretlopt.OPT_S)) && N > 130) {
        err = GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_B) && N > 150) {
        err = GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_D) && N > 16) {
        err = GretlError.badopt;
    } else if ((opt & Gretlopt.OPT_A) && N > 6) {
        err = GretlError.badopt;
    }
    if (err) {
        gretl_errmsg_set("Too many groups for the specified plot");
        return err;
    }

    /* select a default if no panplot-specific option given */
    if (!(opt & (Gretlopt.OPT_M | Gretlopt.OPT_V | Gretlopt.OPT_S | Gretlopt.OPT_D |
                 Gretlopt.OPT_A | Gretlopt.OPT_B | Gretlopt.OPT_C))) {
        if (N <= 130) {
            opt |= Gretlopt.OPT_V; /* --overlay */
        } else if (N <= 150) {
            opt |= Gretlopt.OPT_B; /* --boxplots */
        } else {
            opt |= Gretlopt.OPT_M; /* --means */
        }
    }

    /* transcribe any gnuplot-compatible flags from the
       incoming @opt to @gp_opt
    */
    transcribe_option_flags(&gp_opt, opt, Gretlopt.OPT_U | Gretlopt.OPT_B | Gretlopt.OPT_Y | Gretlopt.OPT_W);

    if (opt & (Gretlopt.OPT_U | Gretlopt.OPT_B)) {
        /* handle output or outbuf spec */
        Gretlopt outopt = (opt & Gretlopt.OPT_U) ? Gretlopt.OPT_U : Gretlopt.OPT_B;
        const(char) * s = get_optval_string(GretlCmdIndex.panplot, outopt);
        int pci = (opt & (Gretlopt.OPT_B | Gretlopt.OPT_C)) ? GretlCmdIndex.bxplot : GretlCmdIndex.gnuplot;

        if (s != null) {
            set_optval_string(cast(GretlCmdIndex) pci, outopt, s);
        }
    }

    set_effective_plot_ci(GretlCmdIndex.panplot);

    if (opt & Gretlopt.OPT_M) {
        /* --means */
        err = panel_means_ts_plot(vnum, literal, dset, gp_opt);
    } else if (opt & Gretlopt.OPT_V) {
        /* --overlay */
        err = panel_overlay_ts_plot(vnum, literal, dset, gp_opt);
    } else if (opt & Gretlopt.OPT_S) {
        /* --sequence */
        err = gnuplot(list, literal, dset, gp_opt | Gretlopt.OPT_O | Gretlopt.OPT_T);
    } else if (opt & Gretlopt.OPT_D) {
        /* --grid (uses multiplot) */
        err = panel_grid_ts_plot(vnum, dset, 0);
    } else if (opt & Gretlopt.OPT_A) {
        /* --stack (uses multiplot) */
        err = panel_grid_ts_plot(vnum, dset, 1);
    } else if (opt & Gretlopt.OPT_B) {
        /* --boxplots */
        gp_opt = cast(Gretlopt)(gp_opt & ~Gretlopt.OPT_Y);
        err = boxplots(list, literal, dset, gp_opt | Gretlopt.OPT_P);
    } else if (opt & Gretlopt.OPT_C) {
        /* --boxplot */
        gp_opt = cast(Gretlopt)(gp_opt & ~Gretlopt.OPT_Y);
        err = boxplots(list, literal, dset, gp_opt);
    }

    reset_effective_plot_ci();

    return err;
}

int data_straddle_zero(const(GretlMatrixRaw) * m)
{
    int t, lt0 = 0, gt0 = 0;

    for (t=0; t<m.rows; t++) {
        if (gretl_matrix_get(m, t, 1) < 0) {
            lt0 = 1;
        }
        if (gretl_matrix_get(m, t, 2) > 0) {
            gt0 = 1;
        }
        if (lt0 && gt0) {
            return 1;
        }
    }

    return 0;
}

void real_irf_print_plot(const(GretlMatrixRaw) * resp, const(char) * targname, const(char) * shockname, const(char) * perlabel, double alpha, int confint, int use_fill, FILE * fp)
{
    int periods = gretl_matrix_rows(resp);
    char * title = null;
    int t;

    if (!confint) {
        fputs("# impulse response plot\n", fp);
    }

    if (confint) {
        fputs("set key left top\n", fp);
        title = g_strdup_printf(_("response of %s to a shock in %s, with bootstrap confidence interval").ptr,
                                targname, shockname);
    } else {
        fputs("set nokey\n", fp);
        title = g_strdup_printf(_("response of %s to a shock in %s").ptr,
                                targname, shockname);
    }

    fprintf(fp, "set xlabel '%s'\n".ptr, perlabel);
    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [-1:%d]\n".ptr, periods);
    fprintf(fp, "set title '%s'\n".ptr, title);
    g_free(title);

    if (confint) {
        double ql = alpha / 2;
        double qh = 1.0 - ql;

        fputs("plot \\\n", fp);
        if (use_fill) {
            title = g_strdup_printf(_("%g percent confidence band").ptr, 100 * (1 - alpha));
            print_filledcurve_line(title, null, 0, fp);
            g_free(title);
            if (data_straddle_zero(resp)) {
                fputs("0 notitle w lines lt 0, \\\n", fp);
            }
            fprintf(fp, "'-' using 1:2 title '%s' w lines lt 1\n".ptr, _("point estimate").ptr);
        } else {
            fprintf(fp, "'-' using 1:2 title '%s' w lines, \\\n".ptr,
                    _("point estimate").ptr);
            title = g_strdup_printf(_("%g and %g quantiles").ptr, ql, qh);
            fprintf(fp, "'-' using 1:2:3:4 title '%s' w errorbars\n".ptr, title);
            g_free(title);
        }
    } else {
        fputs("plot \\\n'-' using 1:2 w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    if (confint && use_fill) {
        for (t=0; t<periods; t++) {
            fprintf(fp, "%d %.10g %.10g\n".ptr, t,
                    gretl_matrix_get(resp, t, 1),
                    gretl_matrix_get(resp, t, 2));
        }
        fputs("e\n", fp);
    }

    for (t=0; t<periods; t++) {
        fprintf(fp, "%d %.10g\n".ptr, t, gretl_matrix_get(resp, t, 0));
    }
    fputs("e\n", fp);

    if (confint && !use_fill) {
        for (t=0; t<periods; t++) {
            fprintf(fp, "%d %.10g %.10g %.10g\n".ptr, t,
                    gretl_matrix_get(resp, t, 0),
                    gretl_matrix_get(resp, t, 1),
                    gretl_matrix_get(resp, t, 2));
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();
}

int gretl_VAR_plot_impulse_response(GretlVar * var, int targ, int shock, int periods, double alpha, const(DatasetRaw) * dset, int error_bars)
{
    int use_fill = !error_bars;
    GretlMatrixRaw * resp;
    int err = 0;

    if (alpha != 0 && (alpha < 0.01 || alpha > 0.5)) {
        return GretlError.data;
    }

    resp = gretl_VAR_get_impulse_response(var, targ, shock, periods,
                                          alpha, dset, &err);

    if (!err) {
        int vtarg = gretl_VAR_get_variable_number(var, targ);
        int vshock = gretl_VAR_get_variable_number(var, shock);
        int confint = (resp.cols > 1);
        FILE * fp;

        fp = open_plot_input_file((confint)? PlotType.irfboot : PlotType.regular, cast(GptFlags) 0, &err);
        if (!err) {
            real_irf_print_plot(resp, dset.varname[vtarg],
                                dset.varname[vshock],
                                dataset_period_label(dset),
                                alpha, confint, use_fill,
                                fp);
            err = finalize_plot_input_file(fp);
        }
        gretl_matrix_free(resp);
    }

    return err;
}

/* Show as a histogram if requested, otherwise do a regular plot */

int gretl_VAR_plot_FEVD(GretlVar * var, int targ, int periods, const(DatasetRaw) * dset, int histogram)
{
    FILE * fp = null;
    GretlMatrixRaw * V;
    char * title;
    int i, t, v;
    PlotType ptype;
    int err = 0;

    V = gretl_VAR_get_FEVD_matrix(var, targ, -1, periods, dset, &err);
    if (V == null) {
        return GretlError.alloc;
    }

    ptype = histogram ? PlotType.stacked_bar : PlotType.regular;

    fp = open_plot_input_file(cast(PlotType) ptype, cast(GptFlags) 0, &err);
    if (err) {
        gretl_matrix_free(V);
        return err;
    }

    v = gretl_VAR_get_variable_number(var, targ);

    fprintf(fp, "set xlabel '%s'\n".ptr, dataset_period_label(dset));
    title = g_strdup_printf(_("forecast variance decomposition for %s").ptr,
                            dset.varname[v]);
    fprintf(fp, "set title '%s'\n".ptr, title);
    g_free(title);

    if (histogram) {
        fputs("set key outside\n", fp);
        fputs("# literal lines = 3\n", fp);
        fputs("set style fill solid 0.35\n", fp);
        fputs("set style histogram rowstacked\n", fp);
        fputs("set style data histogram\n", fp);
        fprintf(fp, "set xrange [-1:%d]\n".ptr, periods);
    } else {
        fputs("set key left top\n", fp);
        fputs("set xzeroaxis\n", fp);
    }

    fputs("set yrange [0:100]\n", fp);
    fputs("plot \\\n", fp);

    for (i=0; i<var.neqns; i++) {
        v = gretl_VAR_get_variable_number(var, i);
        if (histogram) {
            fprintf(fp, "'-' using 2 title \"%s\"".ptr, dset.varname[v]);
        } else {
            fprintf(fp, "'-' using 1:2 title \"%s\" w lines".ptr, dset.varname[v]);
        }
        if (i < var.neqns - 1) {
            fputs(", \\\n", fp);
        } else {
            fputc('\n', fp);
        }
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<var.neqns; i++) {
        for (t=0; t<periods; t++) {
            fprintf(fp, "%d %.4f\n".ptr, t, 100 * gretl_matrix_get(V, t, i));
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    gretl_matrix_free(V);

    return finalize_plot_input_file(fp);
}

int gretl_VAR_plot_multiple_irf(GretlVar * var, int periods, double alpha, const(DatasetRaw) * dset, int error_bars)
{
    FILE * fp = null;
    GptFlags flags = cast(GptFlags) 0;
    int confint = 0;
    int use_fill = !error_bars;
    char * title = null;
    int n = var.neqns;
    int nplots = n * n;
    int vtarg, vshock;
    GretlMatrixRaw * R = null;
    int Rcol, Rstep;
    int t, i, j;
    int err = 0;

    maybe_set_small_font(nplots);

    if (nplots > 12) {
        flags |= GptFlags.xxl;
    } else if (nplots > 9) {
        flags |= GptFlags.xl;
    }

    fp = open_plot_input_file(PlotType.multi_irf, flags, &err);
    if (err) {
        return err;
    }

    fprintf(fp, "set multiplot layout %d,%d\n".ptr, n, n);

    if (n < 4) {
        fprintf(fp, "set xlabel '%s'\n".ptr, dataset_period_label(dset));
    } else {
        fputs("set noxlabel\n", fp);
    }

    fputs("set xzeroaxis\n", fp);
    fprintf(fp, "set xrange [-1:%d]\n".ptr, periods);

    gretl_push_c_numeric_locale();

    /* Use facility to get all impulse responses
       via one call
    */
    R = gretl_VAR_get_impulse_response(var, -1, -1, periods,
                                       alpha, dset, &err);
    if (!err && R.cols > nplots) {
        confint = 1;
    }
    Rcol = 0;
    Rstep = confint ? 3 : 1;

    for (i=0; i<n && !err; i++) {
        vtarg = gretl_VAR_get_variable_number(var, i);

        for (j=0; j<n; j++) {
            vshock = gretl_VAR_get_variable_number(var, j);

            if (i == 0 && j == 0) {
                /* the first plot */
                if (confint) {
                    fputs("set key left top\n", fp);
                } else {
                    fputs("set nokey\n", fp);
                }
            }
            title = g_strdup_printf("%s -> %s".ptr, dset.varname[vshock],
                                    dset.varname[vtarg]);
            fprintf(fp, "set title '%s'\n".ptr, title);
            g_free(title);

            fputs("plot \\\n", fp);
            if (confint && use_fill) {
                print_filledcurve_line(null, null, 0, fp);
                fputs("'-' using 1:2 notitle w lines lt 1\n", fp);
            } else if (confint) {
                fputs("'-' using 1:2 notitle w lines, \\\n", fp);
                fputs("'-' using 1:2:3:4 notitle w errorbars\n", fp);
            } else {
                fputs("'-' using 1:2 notitle w lines\n", fp);
            }

            if (confint && use_fill) {
                for (t=0; t<periods; t++) {
                    fprintf(fp, "%d %.10g %.10g\n".ptr, t,
                            gretl_matrix_get(R, t, Rcol+1),
                            gretl_matrix_get(R, t, Rcol+2));
                }
                fputs("e\n", fp);
            }

            for (t=0; t<periods; t++) {
                fprintf(fp, "%d %.10g\n".ptr, t, gretl_matrix_get(R, t, Rcol));
            }
            fputs("e\n", fp);

            if (confint && !use_fill) {
                for (t=0; t<periods; t++) {
                    fprintf(fp, "%d %.10g %.10g %.10g\n".ptr, t,
                            gretl_matrix_get(R, t, Rcol),
                            gretl_matrix_get(R, t, Rcol+1),
                            gretl_matrix_get(R, t, Rcol+2));
                }
                fputs("e\n", fp);
            }
            Rcol += Rstep;
        }
    }
    gretl_matrix_free(R);

    gretl_pop_c_numeric_locale();

    if (err) {
        fclose(fp);
        return err;
    }

    fputs("unset multiplot\n", fp);

    return finalize_plot_input_file(fp);
}

void write_xrange(FILE * fp, const(double) * obs, int t1, int nobs)
{
    /* for use in gretl_system_residual_plot() and
       gretl_VECM_combined_EC_plot()
    */
    double x0, x1, z;

    x0 = obs[t1];
    z = floor(x0);
    x0 = (x0 - z) > 0.333 ? z : z - 0.5;

    x1 = obs[t1+nobs-1];
    z = ceil(x1);
    x1 = (z - x1) > 0.333 ? z : z + 0.5;

    fprintf(fp, "set xrange [%g:%g]\n".ptr, x0, x1);
}

int gretl_system_residual_plot(void * p, int ci, int eqn, const(DatasetRaw) * dset)
{
    GretlVar * var = null;
    EquationSystem * sys = null;
    const(GretlMatrixRaw) * E = null;
    FILE * fp = null;
    const(double) * obs;
    char[8] lwstr;
    int single = 0;
    int nvars, nobs;
    int i, v, t, t1;
    int imin, imax;
    int err = 0;

    if (ci == GretlCmdIndex.var || ci == GretlCmdIndex.vecm) {
        var = cast(GretlVar *) p;
        E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == GretlCmdIndex.system) {
        sys = cast(EquationSystem *) p;
        E = sys.E;
    }

    if (E == null) {
        return GretlError.data;
    }

    nvars = gretl_matrix_cols(E);
    nobs = gretl_matrix_rows(E);
    t1 = gretl_matrix_get_t1(E);

    if (eqn > 0 && eqn <= nvars) {
        imin = eqn - 1;
        imax = imin + 1;
        single = 1;
    } else {
        imin = 0;
        imax = nvars;
        single = (nvars == 1);
    }

    fp = open_plot_input_file(PlotType.regular, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    obs = gretl_plotx(dset, Gretlopt.OPT_NONE);

    if (quarterly_or_monthly(dset)) {
        fprintf(fp, "# timeseries %d\n".ptr, dset.pd);
    }

    if (!single) {
        fputs("# system residual plot\n", fp);
    }

    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();
    write_xrange(fp, obs, t1, nobs);

    if (ci == GretlCmdIndex.var) {
        fprintf(fp, "set title '%s'\n".ptr, _("VAR residuals").ptr);
    } else {
        fprintf(fp, "set title '%s'\n".ptr, _("System residuals").ptr);
    }

    set_lwstr(lwstr.ptr);

    if (single) {
        fputs("plot ", fp);
    } else {
        fputs("plot \\\n", fp);
    }

    for (i=imin; i<imax; i++) {
        if (var != null) {
            v = gretl_VAR_get_variable_number(var, i);
        } else {
            v = system_get_depvar(sys, i);
        }
        fprintf(fp, "'-' using 1:2 title '%s' w lines%s".ptr,
                dset.varname[v], lwstr.ptr);
        if (i == imax - 1) {
            fputc('\n', fp);
        } else {
            fputs(", \\\n", fp);
        }
    }

    for (i=imin; i<imax; i++) {
        for (t=0; t<nobs; t++) {
            double eti = gretl_matrix_get(E, t, i);

            if (obs != null) {
                fprintf(fp, "%g %.10g\n".ptr, obs[t+t1], eti);
            } else {
                fprintf(fp, "%d %.10g\n".ptr, t+1, eti);
            }
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int gretl_VECM_combined_EC_plot(GretlVar * var, const(DatasetRaw) * dset)
{
    const(GretlMatrixRaw) * EC = null;
    FILE * fp = null;
    const(double) * obs;
    int nvars, nobs;
    int i, t, t1;
    int err = 0;

    EC = VECM_get_EC_matrix(var, dset, &err);
    if (err) {
        return err;
    }

    t1 = gretl_matrix_get_t1(EC);

    fp = open_plot_input_file(PlotType.regular, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    obs = gretl_plotx(dset, Gretlopt.OPT_NONE);

    nvars = gretl_matrix_cols(EC);
    nobs = gretl_matrix_rows(EC);

    fputs("# VECM EC plot\n", fp);
    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);

    gretl_push_c_numeric_locale();
    write_xrange(fp, obs, t1, nobs);

    if (nvars > 1) {
        fprintf(fp, "set title '%s'\n".ptr, _("EC terms").ptr);
    } else {
        fprintf(fp, "set title '%s'\n".ptr, _("EC term").ptr);
    }

    fputs("plot \\\n", fp);
    for (i=0; i<nvars; i++) {
        if (nvars > 1) {
            fprintf(fp, "'-' using 1:2 title 'EC %d' w lines".ptr, i + 1);
        } else {
            fprintf(fp, "'-' using 1:2 notitle w lines".ptr);
        }
        if (i == nvars - 1) {
            fputc('\n', fp);
        } else {
            fputs(", \\\n", fp);
        }
    }

    for (i=0; i<nvars; i++) {
        for (t=0; t<nobs; t++) {
            double eti = gretl_matrix_get(EC, t, i);

            if (obs != null) {
                fprintf(fp, "%g %.10g\n".ptr, obs[t+t1], eti);
            } else {
                fprintf(fp, "%d %.10g\n".ptr, t+1, eti);
            }
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int gretl_system_residual_mplot(void * p, int ci, const(DatasetRaw) * dset)
{
    const(GretlMatrixRaw) * E = null;
    GretlVar * var = null;
    EquationSystem * sys = null;
    FILE * fp = null;
    const(double) * obs;
    double startdate;
    double xmin, xmax, xrange;
    int nvars, nobs, incr;
    int rows, cols;
    int i, v, t, t1;
    int err = 0;

    if (ci == GretlCmdIndex.var || ci == GretlCmdIndex.vecm) {
        var = cast(GretlVar *) p;
        E = gretl_VAR_get_residual_matrix(var);
    } else if (ci == GretlCmdIndex.system) {
        sys = cast(EquationSystem *) p;
        E = sys.E;
    }

    if (E == null) {
        return GretlError.data;
    }

    nvars = gretl_matrix_cols(E);
    if (nvars > 9) {
        return 1;
    }

    obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
    if (obs == null) {
        return GretlError.alloc;
    }

    nobs = gretl_matrix_rows(E);
    t1 = gretl_matrix_get_t1(E);

    fp = open_plot_input_file(PlotType.multi_basic, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    if (nvars <= 3) {
        rows = nvars;
        cols = 1;
    } else if (nvars == 4) {
        rows = cols = 2;
    } else if (nvars == 4) {
        rows = cols = 2;
    } else {
        rows = 3;
        cols = cast(int) (1 + floor((nvars - 0.1) / 3));
    }

    fprintf(fp, "set multiplot layout %d,%d\n".ptr, rows, cols);
    fputs("set nokey\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set noxlabel\n", fp);
    fputs("set noylabel\n", fp);

    gretl_push_c_numeric_locale();

    startdate = obs[t1];
    incr = nobs / (2 * dset.pd);
    if (incr > 0) {
        fprintf(fp, "set xtics %g, %d\n".ptr, ceil(startdate), incr);
    }

    gretl_minmax(t1, t1 + nobs - 1, obs, &xmin, &xmax);
    xrange = xmax - xmin;
    xmin -= xrange * .025;
    xmax += xrange * .025;
    fprintf(fp, "set xrange [%.10g:%.10g]\n".ptr, xmin, xmax);

    for (i=0; i<nvars; i++) {
        if (var != null) {
            v = gretl_VAR_get_variable_number(var, i);
        } else {
            v = system_get_depvar(sys, i);
        }

        fprintf(fp, "set title '%s'\n".ptr, dset.varname[v]);
        fputs("plot '-' using 1:2 with lines\n", fp);

        for (t=0; t<nobs; t++) {
            double eti;

            fprintf(fp, "%.10g\t".ptr, obs[t+t1]);
            eti = gretl_matrix_get(E, t, i);
            write_gp_dataval(eti, fp, 1);
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();
    fputs("unset multiplot\n", fp);

    return finalize_plot_input_file(fp);
}

int gretl_VAR_roots_plot(GretlVar * var)
{
    const(GretlMatrixRaw) * lam;
    FILE * fp = null;
    double x, y;
    double px, py;
    int i, n, err = 0;

    lam = gretl_VAR_get_roots(var, &err);
    if (err) {
        return err;
    }

    fp = open_plot_input_file(PlotType.roots, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    n = gretl_matrix_rows(lam);

    fprintf(fp, "set title '%s'\n".ptr,
            _("VAR inverse roots in relation to the unit circle").ptr);
    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1 notitle w lines, \\\n'-' using 1:2 w points pt 7\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(lam, i, 0);
        y = gretl_matrix_get(lam, i, 1);
        /* in polar form */
        px = atan2(y, x);
        py = sqrt(x * x + y * y);
        fprintf(fp, "%.8f %.8f # %.4f,%.4f\n".ptr, px, py, x, y);
    }

    gretl_pop_c_numeric_locale();

    fputs("e\n", fp);

    return finalize_plot_input_file(fp);
}

/**
 * confidence_ellipse_plot:
 * @V: 2x2 covariance matrix.
 * @b: 2-vector containing point estimates
 * @tcrit: critical t-value for 1 - alpha confidence.
 * @Fcrit: critical F-value for 1 - alpha confidence.
 * @alpha: nominal non-coverage, as decimal.
 * @iname: name of first parameter.
 * @jname: name of second parameter.
 *
 * Plots a 95% confidence ellipse for the parameter estimates
 * in @b with covariance @V.
 *
 * Returns: 0 on success, non-zero on error.
 */

int confidence_ellipse_plot(GretlMatrixRaw * V, double * b, double tcrit, double Fcrit, double alpha, const(char) * iname, const(char) * jname)
{
    FILE * fp = null;
    double[2] maxerr;
    double[2] xcoeff;
    double[2] ycoeff;
    double cval = 100 * (1 - alpha);
    GretlMatrixRaw * e = null;
    char * title;
    int i, err = 0;

    maxerr[0] = tcrit * sqrt(gretl_matrix_get(V, 0, 0));
    maxerr[1] = tcrit * sqrt(gretl_matrix_get(V, 1, 1));

    err = gretl_invert_symmetric_matrix(V);
    if (err) {
        return err;
    }

    e = gretl_symmetric_matrix_eigenvals(V, 1, &err);
    if (err) {
        return err;
    }

    for (i=0; i<2; i++) {
        e.val[i] = sqrt(1.0 / e.val[i] * Fcrit);
        xcoeff[i] = e.val[i] * gretl_matrix_get(V, 0, i);
        ycoeff[i] = e.val[i] * gretl_matrix_get(V, 1, i);
    }

    gretl_matrix_free(e);

    fp = open_plot_input_file(PlotType.ellipse, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    title = g_strdup_printf(_("%g%% confidence ellipse and %g%% marginal intervals").ptr,
                            cval, cval);
    fprintf(fp, "set title '%s'\n".ptr, title);
    g_free(title);

    fputs("# literal lines = 9\n", fp);
    fputs("set parametric\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);

    fprintf(fp, "set xlabel '%s'\n".ptr, iname);
    fprintf(fp, "set ylabel '%s'\n".ptr, jname);
    fprintf(fp, "set label '%.3g, %.3g' at ".ptr, b[0], b[1]);

    gretl_push_c_numeric_locale();

    fprintf(fp, "%g,%g point lt 2 pt 1 offset 3,3\n".ptr, b[0], b[1]);

    fprintf(fp, "x(t) = %g*cos(t)%+g*sin(t)%+g\n".ptr, xcoeff[0], xcoeff[1], b[0]);
    fprintf(fp, "y(t) = %g*cos(t)%+g*sin(t)%+g\n".ptr, ycoeff[0], ycoeff[1], b[1]);

    fputs("plot x(t), y(t) notitle, \\\n", fp);
    fprintf(fp, "%g, y(t) notitle w lines lt 2, \\\n".ptr, b[0] - maxerr[0]);
    fprintf(fp, "%g, y(t) notitle w lines lt 2, \\\n".ptr, b[0] + maxerr[0]);
    fprintf(fp, "x(t), %g notitle w lines lt 2, \\\n".ptr, b[1] - maxerr[1]);
    fprintf(fp, "x(t), %g notitle w lines lt 2\n".ptr, b[1] + maxerr[1]);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

void corrgm_min_max(const(double) * acf, const(double) * pacf, const(GretlMatrixRaw) * PM, int m, double pm, double * pymin, double * pymax)
{
    double ymax = 0;
    double ymin = 0;
    double pmk;
    int PMcol = 0;
    int k;

    if (PM != null && PM.cols == 3) {
        PMcol = 1;
    }

    /* Basic "pm" is relevant if PM (Bartlett) is not present
       or pacf is present.
    */
    if (PM == null || pacf != null) {
        ymax = pm;
        ymin = -pm;
    }

    /* adjust the above based on max and min of ACF, PM and PACF */
    for (k=0; k<m; k++) {
        if (acf[k] > ymax) {
            ymax = acf[k];
        } else if (acf[k] < ymin) {
            ymin = acf[k];
        }
        if (PM != null) {
            pmk = gretl_matrix_get(PM, k, PMcol);
            if (pmk > ymax) {
                ymax = pmk;
            } else if (-pmk < ymin) {
                ymin = -pmk;
            }
        }
        if (pacf != null) {
            if (pacf[k] > ymax) {
                ymax = pacf[k];
            } else if (pacf[k] < ymin) {
                ymin = pacf[k];
            }
        }
    }

    if (pacf == null) {
        /* just create a little "breathing space" */
        ymax *= 1.2;
        ymin *= 1.2;
    } else {
        if (ymax > 0.6) {
            ymax = 1;
        } else {
            ymax *= 1.2;
        }
        if (ymin < -0.6) {
            ymin = -1;
        } else {
            ymin *= 1.2;
        }
        /* make the range symmetrical */
        if (fabs(ymin) > ymax) {
            ymax = -ymin;
        } else if (ymax > fabs(ymin)) {
            ymin = -ymax;
        }
    }

    if (ymax > 1.0) ymax = 1;
    if (ymin < -1.0) ymin = -1;

    *pymax = ymax;
    *pymin = ymin;
}

int real_correlogram_print_plot(const(char) * vname, const(double) * acf, const(double) * pacf, const(GretlMatrixRaw) * PM, int m, double pm, Gretlopt opt, FILE * fp)
{
    char[16] pm_title;
    double ymin, ymax;
    int multi = pacf != null;
    int PMcol = 0;
    int k;

    if (PM != null) {
        PMcol = 1;
    }
    if (multi || PM == null) {
        sprintf(pm_title.ptr, "%.2f/T^%.1f".ptr, 1.96, 0.5);
    }

    corrgm_min_max(acf, pacf, PM, m, pm, &ymin, &ymax);

    gretl_push_c_numeric_locale();

    /* create two separate plots, if both are OK */
    if (pacf != null) {
        fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    }
    fputs("set xzeroaxis\n", fp);
    print_keypos_string(GpKeySpec.right | GpKeySpec.top, fp);
    fprintf(fp, "set xlabel '%s'\n".ptr, _("lag").ptr);

    fprintf(fp, "set yrange [%.2f:%.2f]\n".ptr, ymin, ymax);

    /* upper (or only) plot: ACF */
    if (multi) {
        fputs("set origin 0.0,0.50\n", fp);
    }
    if (opt & Gretlopt.OPT_R) {
        fprintf(fp, "set title '%s'\n".ptr, _("Residual ACF").ptr);
    } else {
        fprintf(fp, "set title '%s %s'\n".ptr, _("ACF for").ptr, vname);
    }
    fprintf(fp, "set xrange [0:%d]\n".ptr, m + 1);

    if (PM != null) {
        /* xgettext:no-c-format */
        const(char) * PM_title = "95% Bartlett interval".ptr;

        fprintf(fp, "plot \\\n'-' using 1:2 notitle w impulses lw 5, \\\n'-' title '%s' w lines lt 2, \\\n'-' notitle w lines lt 2\n".ptr, _(PM_title));
    } else {
        fprintf(fp, "plot \\\n'-' using 1:2 notitle w impulses lw 5, \\\n%g title '+- %s' lt 2, \\\n%g notitle lt 2\n".ptr, pm, pm_title.ptr, -pm);
    }

    /* data for ACF plot */
    for (k=0; k<m; k++) {
        fprintf(fp, "%d %g\n".ptr, k + 1, acf[k]);
    }
    fputs("e\n", fp);

    if (PM != null) {
        /* Bartlett-type confidence band data */
        for (k=0; k<m; k++) {
            fprintf(fp, "%d %g\n".ptr, k + 1, gretl_matrix_get(PM, k, PMcol));
        }
        fputs("e\n", fp);
        for (k=0; k<m; k++) {
            fprintf(fp, "%d -%g\n".ptr, k + 1, gretl_matrix_get(PM, k, PMcol));
        }
        fputs("e\n", fp);
    }

    if (multi) {
        /* lower plot: PACF */
        fputs("set origin 0.0,0.0\n", fp);
        if (opt & Gretlopt.OPT_R) {
            fprintf(fp, "set title '%s'\n".ptr, _("Residual PACF").ptr);
        } else {
            fprintf(fp, "set title '%s %s'\n".ptr, _("PACF for").ptr, vname);
        }
        fprintf(fp, "set xrange [0:%d]\n".ptr, m + 1);
        fprintf(fp, "plot \\\n'-' using 1:2 notitle w impulses lw 5, \\\n%g title '+- %s' lt 2, \\\n%g notitle lt 2\n".ptr, pm, pm_title.ptr, -pm);
        /* PACF data */
        for (k=0; k<m; k++) {
            fprintf(fp, "%d %g\n".ptr, k + 1, pacf[k]);
        }
        fputs("e\n", fp);
    }

    if (multi) {
        fputs("unset multiplot\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

int correlogram_plot(const(char) * vname, const(double) * acf, const(double) * pacf, const(GretlMatrixRaw) * PM, int m, double pm, Gretlopt opt)
{
    FILE * fp;
    int err = 0;

    fp = open_plot_input_file(PlotType.correlogram, cast(GptFlags) 0, &err);

    if (!err) {
        real_correlogram_print_plot(vname, acf, pacf,
                                    PM, m, pm, opt, fp);
        err = finalize_plot_input_file(fp);
    }

    return err;
}

int roundup_mod(int i, double x)
{
    return cast(int) ceil(cast(double) x * i);
}

/* option flags:
   OPT_R --radians: use radians as unit
   OPT_D --degrees: use degrees as unit
   OPT_L --log: use log scale
*/

int real_pergm_plot(const(char) * vname, int T, int L, const(double) * dens, Gretlopt opt, FILE * fp)
{
    char[80] s;
    double ft;
    int T2 = T / 2;
    int k, t, err = 0;

    fputs("set xtics nomirror\n", fp);

    fprintf(fp, "set x2label '%s'\n".ptr, _("periods").ptr);
    fprintf(fp, "set x2range [0:%d]\n".ptr, roundup_mod(T, 2.0));

    fputs("set x2tics (", fp);
    k = T2 / 6;
    for (t = 1; t <= T2; t += k) {
        fprintf(fp, "\"%.1f\" %d, ".ptr, cast(double) T / t, 4 * t);
    }
    fprintf(fp, "\"\" %d)\n".ptr, 2 * T);

    if (opt & Gretlopt.OPT_R) {
        fprintf(fp, "set xlabel '%s'\n".ptr, _("radians").ptr);
    } else if (opt & Gretlopt.OPT_D) {
        fprintf(fp, "set xlabel '%s'\n".ptr, _("degrees").ptr);
    } else {
        fprintf(fp, "set xlabel '%s'\n".ptr, _("scaled frequency").ptr);
    }

    fputs("set xzeroaxis\n", fp);
    fputs("set nokey\n", fp);

    /* open gnuplot title string */
    fputs("set title '", fp);

    if (vname == null) {
        fputs(_("Residual spectrum").ptr, fp);
    } else {
        sprintf(s.ptr, _("Spectrum of %s").ptr, vname);
        fputs(s.ptr, fp);
    }

    if (opt & Gretlopt.OPT_O) {
        fputs(" (", fp);
        fprintf(fp, _("Bartlett window, length %d").ptr, L);
        fputc(')', fp);
    }

    if (opt & Gretlopt.OPT_L) {
        fputs(" (", fp);
        fputs(_("log scale").ptr, fp);
        fputc(')', fp);
    }

    /* close gnuplot title string */
    fputs("'\n", fp);

    gretl_push_c_numeric_locale();

    if (opt & Gretlopt.OPT_R) {
        /* frequency scale in radians */
        fputs("set xrange [0:3.1416]\n", fp);
    } else if (opt & Gretlopt.OPT_D) {
        /* frequency scale in degrees */
        fputs("set xrange [0:180]\n", fp);
    } else {
        /* data-scaled frequency */
        fprintf(fp, "set xrange [0:%d]\n".ptr, roundup_mod(T, 0.5));
    }

    if (!(opt & Gretlopt.OPT_L)) {
        fprintf(fp, "set yrange [0:%g]\n".ptr, 1.2 * gretl_max(0, T/2 - 1, dens));
    }

    if (opt & Gretlopt.OPT_R) {
        fputs("set xtics (\"0\" 0, \"π/4\" pi/4, \"π/2\" pi/2, \"3π/4\" 3*pi/4, \"π\" pi)\n", fp);
    }

    fputs("plot '-' using 1:2 w lines\n", fp);

    for (t=1; t<=T2; t++) {
        if (opt & Gretlopt.OPT_R) {
            ft = PI * cast(double) t / T2;
        } else if (opt & Gretlopt.OPT_D) {
            ft = 180 * cast(double) t / T2;
        } else {
            ft = t;
        }
        fprintf(fp, "%g %g\n".ptr, ft, (opt & Gretlopt.OPT_L)? log(dens[t-1]) : dens[t-1]);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    return err;
}

int periodogram_plot(const(char) * vname, int T, int L, const(double) * x, Gretlopt opt)
{
    FILE * fp;
    int err = 0;

    fp = open_plot_input_file(PlotType.periodogram, cast(GptFlags) 0, &err);

    if (!err) {
        real_pergm_plot(vname, T, L, x, opt, fp);
        err = finalize_plot_input_file(fp);
    }

    return err;
}

int arma_spectrum_plot(Model * pmod, const(DatasetRaw) * dset, Gretlopt opt)
{
    GretlMatrixRaw * pdata = null;
    FILE * fp;
    int err = 0;

    pdata = arma_spectrum_plot_data(pmod, dset, &err);
    if (err) {
        return err;
    }

    fp = open_plot_input_file(PlotType.periodogram, cast(GptFlags) 0, &err);

    if (!err) {
        double px, pRe, pIm, scale = pmod.nobs * M_2PI;
        int i, grid = pdata.rows;

        gretl_push_c_numeric_locale();

        fprintf(fp, "set xrange [0:%g]\n".ptr, PI);
        switch (dset.pd) {
        case 12:
            fputs("set xtics (\"0\" 0, \"π/6\" pi/6, \"π/3\" pi/3, \"π/2\" pi/2, \"2π/3\" 2*pi/3, \"5π/6\" 5*pi/6, \"π\" pi)\n", fp);
            break;
        case 6:
            fputs("set xtics (\"0\" 0, \"π/3\" pi/3, \"2π/3\" 2*pi/3, \"π\" pi)\n", fp);
            break;
        case 5:
            fputs("set xtics (\"0\" 0, \"π/5\" pi/5, \"2π/5\" 2*pi/5, \"3π/5\" 3*pi/5, \"4π/5\" 4*pi/5, \"π\" pi)\n", fp);
            break;
        default:
            fputs("set xtics (\"0\" 0, \"π/4\" pi/4, \"π/2\" pi/2, \"3π/4\" 3*pi/4, \"π\" pi)\n", fp);
        }
        fprintf(fp, "set title '%s (%s)'\n".ptr, _("Sample periodogram vs ARMA Spectrum").ptr,
                _("log scale").ptr);
        fprintf(fp, "plot '-' using 1:2 w lines title '%s' lw 2, \\\n".ptr, _("spectrum").ptr);
        fprintf(fp, "'-' using 1:2 w lines title '%s' lw 0.5\n".ptr, _("periodogram").ptr);

        for (i=0; i<grid; i++) {
            fprintf(fp, "%7.5f %12.7f\n".ptr, gretl_matrix_get(pdata, i, 0),
                    log(gretl_matrix_get(pdata, i, 1)));
        }
        fputs("e\n", fp);

        for (i=0; i<grid; i++) {
            pRe = gretl_matrix_get(pdata, i, 2);
            pIm = gretl_matrix_get(pdata, i, 3);
            px = (pRe * pRe + pIm * pIm) / scale;
            fprintf(fp, "%7.5f %12.7f\n".ptr, gretl_matrix_get(pdata, i, 0), log(px));
        }
        fputs("e\n", fp);

        gretl_pop_c_numeric_locale();
        err = finalize_plot_input_file(fp);
    }

    return err;
}

/* Probability of non-exceedance of the kth value in a set of n
   rank-ordered values.  See L. Makkonen, 'Bringing Closure to the
   Plotting Position Controversy', Communications in Statistics -
   Theory and Methods, vol 37, January 2008, for an argument in favor
   of using k/(n + 1); but also see many uses of (k - 1/2)/n in
   the literature.
*/

double plotpos(int k, int n)
{
    return (k - 0.5) / n;
}

double quantile_interp(const(double) * y, int n, double ftarg)
{
    double f, ret = NADBL;
    int i;

    for (i=0; i<n; i++) {
        f = plotpos(i+1, n);
        if (f >= ftarg) {
            if (f > ftarg && i > 0) {
                double f0 = plotpos(i, n);
                double d = (ftarg - f0) / (f - f0);

                ret = (1-d) * y[i-1] + d * y[i];
            } else {
                ret = y[i];
            }
            break;
        }
    }

    return ret;
}

int qq_plot_two_series(const(int) * list, const(DatasetRaw) * dset)
{
    double * x = null;
    double * y = null;
    double f, qx, qy;
    FILE * fp = null;
    int vx = list[1];
    int vy = list[2];
    int nx = 10, ny = 10;
    int i, n, err = 0;

    x = gretl_sorted_series(vx, dset, Gretlopt.OPT_NONE, &nx, &err);

    if (!err) {
        y = gretl_sorted_series(vy, dset, Gretlopt.OPT_NONE, &ny, &err);
        if (err) {
            free(x);
            x = null;
        }
    }

    if (!err) {
        /* take the smaller sample as basis */
        n = (nx > ny)? ny : nx;
    }

    if (!err) {
        fp = open_plot_input_file(PlotType.qq, cast(GptFlags) 0, &err);
    }

    if (err) {
        free(x);
        free(y);
        return err;
    }

    fprintf(fp, "set title '%s'\n".ptr, _("Q-Q plot").ptr);
    gnuplot_missval_string(fp);
    fputs("set key top left\n", fp);
    fprintf(fp, "set xlabel '%s'\n".ptr, plotname(dset, vx, 1));
    fprintf(fp, "set ylabel '%s'\n".ptr, plotname(dset, vy, 1));
    fputs("plot \\\n", fp);
    fputs(" '-' using 1:2 notitle w points, \\\n", fp);
    fputs(" x notitle w lines\n", fp);

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        f = plotpos(i+1, n);

        if (nx == ny) {
            qx = x[i];
            qy = y[i];
        } else if (nx == n) {
            qx = x[i];
            qy = quantile_interp(y, ny, f);
        } else {
            qx = quantile_interp(x, nx, f);
            qy = y[i];
        }

        if (!na(qx) && !na(qy)) {
            fprintf(fp, "%.12g %.12g\n".ptr, qx, qy);
        }
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    free(x);
    free(y);

    return finalize_plot_input_file(fp);
}

/* Option flags: OPT_R --raw, OPT_Z --z-scores */

int normal_qq_plot(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    GptFlags flags = cast(GptFlags) 0;
    char * title = null;
    int zscores = 0;
    double ym = 0, ys = 1;
    double p, qx, qy;
    double * y = null;
    FILE * fp = null;
    int v = list[1];
    int i, n = 20;
    int err = 0;

    y = gretl_sorted_series(v, dset, Gretlopt.OPT_NONE, &n, &err);

    if (!err && y[0] == y[n-1]) {
        gretl_errmsg_sprintf(_("%s is a constant").ptr, dset.varname[v]);
        err = GretlError.data;
    }

    if (err) {
        return err;
    }

    if (opt & Gretlopt.OPT_Z) {
        /* standardize the data */
        zscores = 1;
    }

    if (opt & Gretlopt.OPT_G) {
    flags = GptFlags.icon;
    }

    if (!(opt & Gretlopt.OPT_R)) {
        /* without the --raw option */
        ym = gretl_mean(0, n-1, y);
        ys = gretl_stddev(0, n-1, y);

        if (zscores) {
            /* standardize y */
            for (i=0; i<n; i++) {
                y[i] = (y[i] - ym) / ys;
            }
        }
    }

    fp = open_plot_input_file(PlotType.qq, flags, &err);
    if (err) {
        free(y);
        return err;
    }

    title = g_strdup_printf(_("Q-Q plot for %s").ptr, plotname(dset, v, 1));
    fprintf(fp, "set title '%s'\n".ptr, title);
    g_free(title);
    gnuplot_missval_string(fp);
    fprintf(fp, "set xlabel '%s'\n".ptr, _("Normal quantiles").ptr);

    if (opt & Gretlopt.OPT_R) {
        fputs("set nokey\n", fp);
        fputs("plot \\\n", fp);
        fputs(" '-' using 1:2 notitle w points\n", fp);
    } else {
        fputs("set key top left\n", fp);
        fputs("plot \\\n", fp);
        fputs(" '-' using 1:2 notitle w points, \\\n", fp);
        fputs(" x title 'y = x' w lines\n", fp);
    }

    gretl_push_c_numeric_locale();

    for (i=0; i<n; i++) {
        p = plotpos(i+1, n);
        /* empirical quantile */
        qy = y[i];
        /* normal quantile */
        qx = normal_critval(1 - p);
        if (!na(qx) && !zscores && !(opt & Gretlopt.OPT_R)) {
            qx = ys * qx + ym;
        }
        if (!na(qx) && !na(qy)) {
            fprintf(fp, "%.12g %.12g\n".ptr, qx, qy);
        }
    }

    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    free(y);

    return finalize_plot_input_file(fp);
}

int qq_plot(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    int err;

    if (list[0] == 1) {
        /* one series against normal */
        err = normal_qq_plot(list, dset, opt);
    } else if (list[0] == 2) {
        /* two empirical series */
        err = qq_plot_two_series(list, dset);
    } else {
        err = GretlError.data;
    }

    return err;
}

int kd_plot(const(int) * list, const(DatasetRaw) * dset, Gretlopt opt)
{
    alias KdfuncFn = extern (C) int function(const(double) *, int, double, const(char) *, Gretlopt);
    KdfuncFn kdfunc;
    int err = 0;

    kdfunc = cast(KdfuncFn) get_plugin_function("kernel_density");

    if (kdfunc == null) {
        err = GretlError.fopen;
    } else {
        int v = list[1];
        const(double) * y = dset.Z[v] + dset.t1;
        const(char) * label = dset.varname[v];
        int n = sample_size(dset);
        double bws = 1.0;

        if (opt & Gretlopt.OPT_S) {
            /* bandwidth scale */
            bws = get_optval_double(GretlCmdIndex.kdplot, Gretlopt.OPT_S, &err);
        }
        if (!err) {
            err = kdfunc(y, n, bws, label, opt);
        }
    }

    return err;
}

int pd_from_compfac(const(DatasetRaw) * dset, int compfac, char * stobs)
{
    int pd = -1;

    if (dset.pd == 1 && (compfac == 12 || compfac == 4)) {
        /* annual from monthly or quarterly */
        pd = compfac;
    } else if (dset.pd == 4 && compfac == 3) {
        /* quarterly from monthly */
        pd = 12;
    } else if (dset.pd == 4) {
        /* maybe quarterly from daily? */
        if (compfac >= 60 && compfac <= 69) {
            pd = 5;
        } else if (compfac >= 71 && compfac <= 81) {
            pd = 6;
        } else if (compfac >= 82 && compfac <= 93) {
            return 7;
        }
    } else if (dset.pd == 12) {
        /* maybe monthly from daily? */
        if (compfac >= 20 && compfac <= 23) {
            pd = 5;
        } else if (compfac >= 24 && compfac <= 27) {
            pd = 6;
        } else if (compfac >= 28 && compfac <= 31) {
            pd = 7;
        }
    }

    if (pd > 0) {
        char * p;
        char[OBSLEN] tmp;
        int y, q, m;

        ntolabel(tmp.ptr, dset.t1, dset);
        y = atoi(tmp.ptr);
        p = strchr(tmp.ptr, ':');

        if ((dset.pd == 4 || dset.pd == 12) && p == null) {
            return -1;
        }

        if (dset.pd == 1) {
            if (pd == 4) {
                sprintf(stobs, "%d:1".ptr, y);
            } else if (pd == 12) {
                sprintf(stobs, "%d:01".ptr, y);
            }
        } else if (dset.pd == 4) {
            q = atoi(p + 1);
            m = (q==1)? 1 : (q==2)? 4 : (q==3)? 7 : 10;
            if (pd == 12) {
                sprintf(stobs, "%d:%02d".ptr, y, m);
            } else {
                /* daily */
                sprintf(stobs, "%d-%02d-01".ptr, y, m);
            }
        } else if (dset.pd == 12) {
            /* daily */
            m = atoi(p + 1);
            sprintf(stobs, "%d-%02d-01".ptr, y, m);
        }
    }

    return pd;
}

void transcribe_graph_name(DatasetRaw * targ, int i, const(DatasetRaw) * src, int j)
{
    const(char) * s = series_get_display_name(src, j);

    if (s != null && *s != '\0') {
        series_record_display_name(targ, i, s);
    }
}

/* high-frequency plot for MIDAS */

int hf_plot(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    DatasetRaw * hset;
    double xit;
    char[OBSLEN] stobs;
    int * gplist = null;
    int * hflist = null;
    int * lflist = null;
    char * mylit = null;
    char * p;
    Gretlopt plotopt = Gretlopt.OPT_T;
    int plotpd = 0;
    int nv, nlf = 0;
    int cfac;
    int i, s, t, T;
    int err;

    if (list == null || list[0] < 3) {
        return GretlError.invarg;
    } else if (!dataset_is_time_series(dset)) {
        return GretlError.pdwrong;
    }

    if (gretl_list_has_separator(list)) {
        err = gretl_list_split_on_separator(list, &hflist, &lflist);
        if (err) {
            return err;
        } else {
            cfac = hflist[0];
            nlf = lflist[0];
            nv = 2 + nlf;
        }
    } else {
        cfac = list[0];
        nv = 2;
    }

    T = sample_size(dset) * cfac;

    hset = create_auxiliary_dataset(nv, T, Gretlopt.OPT_NONE);
    if (hset == null) {
        return GretlError.alloc;
    }

    /* set the hf series name */
    strcpy(hset.varname[1], dset.varname[list[1]]);
    p = strrchr(hset.varname[1], '_');
    if (p != null) {
        *p = '\0';
    }
    transcribe_graph_name(hset, 1, dset, list[1]);

    s = 0;
    /* transcribe high-frequency data */
    for (t=dset.t1; t<=dset.t2; t++) {
        for (i=cfac; i>0; i--) {
            xit = dset.Z[list[i]][t];
            hset.Z[1][s++] = xit;
        }
    }

    if (lflist != null) {
        /* add low-frequency term(s), if any */
        for (i=1; i<=nlf; i++) {
            int vi = lflist[i];

            strcpy(hset.varname[i+1], dset.varname[vi]);
            transcribe_graph_name(hset, i+1, dset, vi);
            for (s=0; s<hset.n; s++) {
                hset.Z[i+1][s] = NADBL;
            }
            s = 0;
            for (t=dset.t1; t<=dset.t2; t++) {
                xit = dset.Z[vi][t];
                hset.Z[i+1][s] = xit;
                s += cfac;
            }
        }
    }

    gplist = gretl_consecutive_list_new(1, nv - 1);

    if (lflist != null) {
        free(lflist);
        lflist = gretl_consecutive_list_new(2, nv - 1);
    }

    /* try to set a suitable time-series interpretation
       on the data to be plotted
    */
    plotpd = pd_from_compfac(dset, cfac, stobs.ptr);
    if (plotpd > 0) {
        char[12] numstr;

        sprintf(numstr.ptr, "%d".ptr, plotpd);
        set_obs(numstr.ptr, stobs.ptr, hset, Gretlopt.OPT_T);
    }

    if (opt & Gretlopt.OPT_O) {
        plotopt |= Gretlopt.OPT_O;
    }
    if (opt & Gretlopt.OPT_U) {
        plotopt |= Gretlopt.OPT_U;
    }

    if (literal == null) {
        const(char) * pdstr = midas_pdstr(dset, cfac);
        char * title;

        title = g_strdup_printf("%s (%s)".ptr, hset.varname[1], _(pdstr));
        mylit = g_strdup_printf("{ set ylabel ''; set title '%s'; }".ptr, title);
        g_free(title);
    }

    set_effective_plot_ci(GretlCmdIndex.hfplot);
    na_skiplist = lflist; /* file-scope global */
    err = gnuplot(gplist, literal != null ? literal : mylit,
                  hset, plotopt);
    na_skiplist = null;
    reset_effective_plot_ci();

    free(gplist);
    free(hflist);
    free(lflist);
    g_free(mylit);

    destroy_dataset(hset);

    return err;
}

/**
 * xy_plot_with_control:
 * @list: list of variables by ID number: Y, X, control(s).
 * @literal: extra gnuplot commands or %NULL.
 * @dset: dataset struct.
 * @opt: can add "gnuplot" options.
 *
 * Constructs a scatterplot of modified Y against modified X,
 * where the modification consists in taking the residuals from
 * OLS regression of the variable in question on the control
 * variables, a la Frisch-Waugh-Lovell.
 *
 * Returns: 0 on success, non-zero on error.
 */

int xy_plot_with_control(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    int t1 = dset.t1, t2 = dset.t2;
    int * mlist = null;
    char[MAXDISP] dname;
    Model mod;
    DatasetRaw * gset = null;
    int l0, vy, vx;
    int i, s, t, T;
    int missvals = 0;
    int skip = 0;
    int err = 0;

    if (list == null || list[0] < 3) {
        return GretlError.args;
    }

    l0 = list[0];
    vy = list[1];
    vx = list[2];

    list_adjust_sample(list, &t1, &t2, dset, &missvals);

    /* maximum usable observations */
    T = t2 - t1 + 1 - missvals;

    if (T <= l0) {
        return GretlError.df;
    }

    /* create temporary dataset */
    gset = create_auxiliary_dataset(l0 + 1, T, Gretlopt.OPT_NONE);
    if (gset == null) {
        return GretlError.alloc;
    }

    sprintf(dname.ptr, _("adjusted %s").ptr, dset.varname[vy]);
    series_set_display_name(gset, 1, dname.ptr);
    sprintf(dname.ptr, _("adjusted %s").ptr, dset.varname[vx]);
    series_set_display_name(gset, 2, dname.ptr);

    s = 0;
    for (t=t1; t<=t2; t++) {
        skip = 0;
        for (i=1; i<=list[0]; i++) {
            if (na(dset.Z[list[i]][t])) {
                skip = 1;
                break;
            }
        }
        if (!skip) {
            for (i=1; i<=l0; i++) {
                gset.Z[i][s] = dset.Z[list[i]][t];
            }
            s++;
        }
    }

    /* construct regression list */
    mlist = gretl_list_new(l0);
    for (i=3; i<=l0; i++) {
        mlist[i] = i;
    }

    /* regress Y on Z and save the residuals in series 1 */
    mlist[1] = 1;
    mod = lsq(mlist, gset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = mod.errcode;
    if (err) {
        clear_model(&mod);
        goto bailout;
    } else {
        for (t=0; t<mod.nobs; t++) {
            gset.Z[1][t] = mod.uhat[t];
        }
        clear_model(&mod);
    }

    /* regress X on Z and save the residuals in series 2 */
    mlist[1] = 2;
    mod = lsq(mlist, gset, GretlCmdIndex.ols, Gretlopt.OPT_A);
    err = mod.errcode;
    if (err) {
        clear_model(&mod);
        goto bailout;
    } else {
        for (t=0; t<mod.nobs; t++) {
            gset.Z[2][t] = mod.uhat[t];
        }
        clear_model(&mod);
    }

    /* call for scatter of Y-residuals against X-residuals */
    mlist[0] = 2;
    mlist[1] = 1;
    mlist[2] = 2;
    err = gnuplot(mlist, literal, gset, opt | Gretlopt.OPT_C);

 bailout:

    /* trash the temporary dataset */
    destroy_dataset(gset);
    free(mlist);

    return err;
}

int is_auto_fit_string(const(char) * s)
{
    /* FIXME? */
    if (strstr(s, "automatic fit")) return 1;
    if (strstr(s, _("with least squares fit").ptr)) return 1;
    return 0;
}

void inbuf_inject_literal(const(char) * buf, const(char) * literal, FILE * fp)
{
    const(char) * p = strstr(buf, "\nplot ");

    if (p != null) {
        int n = cast(int) (p - buf);

        fwrite(buf, 1, n, fp);
        print_gnuplot_literal_lines(literal, fp);
        fputs(p + 1, fp);
    } else {
        fputs(buf, fp);
    }
}

void infile_inject_literal(FILE * fin, const(char) * literal, FILE * fout)
{
    char[1024] line;

    while (fgets(line.ptr, line.sizeof, fin)) {
        if (!strncmp(line.ptr, "plot ", 5)) {
            print_gnuplot_literal_lines(literal, fout);
        }
        fputs(line.ptr, fout);
    }
}

/**
 * gnuplot_process_input:
 * @literal: literal "extra" gnuplot commands.
 * @opt: should be OPT_I (input file) or Gretlopt.OPT_I (input buffer)
 * @prn: gretl printing struct.
 *
 * Respond to the "gnuplot" command with %OPT_I, to specify
 * that input should be taken from a user-created gnuplot
 * command file, or %Gretlopt.OPT_I, to take input from a named
 * string variable. If @literal is non-NULL we'll inject
 * the literal gnuplot commands immediately preceding the
 * "plot" line in the file or buffer input.
 *
 * Returns: 0 on success, error code on error.
 */

int gnuplot_process_input(const(char) * literal, Gretlopt opt, Prn * prn)
{
    const(char) * iname = null;
    const(char) * buf = null;
    FILE * fq, fp = null;
    int err = 0;

    if (opt & Gretlopt.OPT_I) {
        iname = get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_I);
    } else if (opt & Gretlopt.OPT_I) {
        iname = get_optval_string(cast(GretlCmdIndex) plot_ci, Gretlopt.OPT_I);
    }

    if (iname != null && *iname != '\0') {
        if (opt & Gretlopt.OPT_I) {
            /* open the input file for reading */
            fp = gretl_fopen(iname, "r");
        } else {
            /* find the input buffer */
            buf = get_string_by_name(iname);
        }
    }

    if (fp == null && buf == null) {
        gretl_errmsg_set("Couldn't find the specified input");
        return GretlError.invarg;
    }

    if (opt & Gretlopt.OPT_W) {
        /* --font=<fontspec> */
        maybe_record_font_choice(Gretlopt.OPT_W);
    }
    if (opt & Gretlopt.OPT_E) {
        /* --alpha=<val> */
        maybe_record_alpha(Gretlopt.OPT_E);
    }

    /* open our own file for writing */
    fq = open_plot_input_file(PlotType.user, cast(GptFlags) 0, &err);

    if (!err) {
        if (literal != null && *literal != '\0') {
            if (fp != null) {
                infile_inject_literal(fp, literal, fq);
            } else {
                inbuf_inject_literal(buf, literal, fq);
            }
        } else if (fp != null) {
            char[1024] line;

            while (fgets(line.ptr, line.sizeof, fp)) {
                fputs(line.ptr, fq);
            }
        } else {
            fputs(buf, fq);
        }
        err = finalize_plot_input_file(fq);
    }

    if (fp != null) {
        fclose(fp);
    }

    return err;
}

/* geoplot-specific functions */

int inline_map_data(const(char) * datfile, FILE * fp)
{
    char[8192] buf;
    FILE * fsrc;
    size_t n;

    fsrc = gretl_fopen(datfile, "rb");
    if (fsrc == null) {
        gretl_errmsg_sprintf(_("Couldn't open %s").ptr, datfile);
        return GretlError.fopen;
    }

    fputs("$MapData << EOD\n", fp);
    while ((n = cast(int) fread(buf.ptr, 1, buf.sizeof, fsrc)) > 0) {
        fwrite(buf.ptr, 1, n, fp);
    }
    fputs("EOD\n", fp);

    fclose(fsrc);
    gretl_remove(datfile); /* not needed any more */

    return 0;
}

GretlMatrixRaw * geoplot_dimensions(double * xlim, double * ylim, int height, int have_payload, int * non_standard)
{
    GretlMatrixRaw * ret = gretl_matrix_alloc(1, 2);
    double xyr = (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]);
    int width;

    if (*non_standard == 0) {
        if (fabs(ylim[0]) > 180 || fabs(ylim[1]) > 180 ||
            fabs(xlim[0]) > 360 || fabs(xlim[1]) > 360) {
            /* Quick and dirty check in the absence of prior
               information that the X, Y data are not in degrees
               of latitude and longitude.
            */
            fprintf(stderr, "alt coordinates: X %g to %g, Y %g to %g\n".ptr,
                    xlim[0], xlim[1], ylim[0], ylim[1]);
            *non_standard = 1;
        }
    }

    if (*non_standard == 0) {
        /* We'll calculate a width-to-height ratio which varies
           inversely with the (mid-point of) latitude so that we don't
           don't get severe stretching in the x dimension when showing
           an area far from the equator. In effect this is a cheap and
           cheerful version of Mercator.
        */
        double ymid = (ylim[0] + ylim[1]) / 2;

        xyr *= cos(ymid * PI/180);
    }

    if (have_payload) {
        /* 1.05 is to compensate for the colorbox */
        width = cast(int) floor(xyr * height * 1.05);
    } else {
        width = cast(int) floor(xyr * height);
    }

    set_special_plot_size(width, height);

    ret.val[0] = width;
    ret.val[1] = height;

    return ret;
}

void fputs_literal(const(char) * s, FILE * fp)
{
    char * tmp = g_strdup(s);

    fputs(g_strchomp(tmp), fp);
    fputc('\n', fp);
    g_free(tmp);
}

const(char) * map_linecolor(const(char) * optlc, int have_payload, NaAction action)
{
    if (optlc != null) {
        /* respect the user's choice */
        return optlc;
    } else if (have_payload) {
        return (action == NaAction.outline)? "gray" : "white";
    } else {
        /* outlines only */
        return "black";
    }
}

void output_map_plot_lines(MapInfo * mi, const(char) * datasrc, const(char) * optlc, int have_payload, int do_key, double linewidth, FILE * fp)
{
    const(char) * lc = map_linecolor(optlc, have_payload, mi.na_action);
    char * bline = g_strdup_printf("lc '%s' lw %g".ptr, lc, linewidth);

    if (have_payload) {
        int do_lines = linewidth > 0;
        const(char) * with_ = "with filledcurves fc palette";
        const(char) * cont = ", \\\n";
        int i, nv = mi.n_codes;

        if (do_key) {
            const(char) * kp = gretl_bundle_get_string(mi.opts, "keypos", null);

            /* ensure colorbox is omitted and key boxes are filled */
            fputs("unset colorbox\n", fp); /* should be handled already? */
            fputs("set style fill solid\n", fp);
            if (kp != null && *kp != '\0') {
                fprintf(fp, "set key %s\n".ptr, kp);
            } else {
                /* this seems to work better than the default on average? */
                fputs("set key bottom right\n", fp);
            }
        }

        /* polygons */
        fprintf(fp, "plot for [i=0:*] %s index i %s notitle%s".ptr, datasrc, with_,
                (do_lines || do_key)? cont : "\n".ptr);
        if (linewidth > 0) {
            /* plus feature outlines */
            fprintf(fp, "  %s using 1:2 with lines %s notitle%s".ptr, datasrc, bline,
                    do_key ? cont : "\n".ptr);
        }
        if (do_key) {
            /* plus key for discrete payload */
            int v;
    int v0 = cast(int) mi.zvals.val[0];

            for (i=0; i<nv; i++) {
                v = i + v0;
                if (mi.zlabels != null) {
                    if (strcmp(mi.zlabels[i], "empty string")) {
                        fprintf(fp, "keyentry with boxes fc palette cb %d title \"%s\"%s".ptr,
                                v, mi.zlabels[i], (i < nv - 1)? cont : "\n".ptr);
                    } else if (i == nv - 1) {
                        fputc('\n', fp);
                    }
                } else {
                    fprintf(fp, "keyentry with boxes fc palette cb %d title \"%s=%d\"%s".ptr,
                            v, mi.zname, v, (i < nv - 1)? cont : "\n".ptr);
                }
            }
        }
    } else {
        /* just show feature outlines */
        fprintf(fp, "plot %s using 1:2 with lines %s\n".ptr, datasrc, bline);
    }

    g_free(bline);
}

int map_do_inlining(MapInfo * mi)
{
    if (gretl_bundle_get_int(mi.opts, "inlined", null)) {
        return 1;
    } else if (mi.flags & MapFlags.multi) {
        return 1;
    } else {
        return 0;
    }
}

/* called from the geoplot plugin to finalize a map */

int write_map_gp_file(void * ptr)
{
    MapInfo * mi = cast(MapInfo *) ptr;
    GretlBundle * opts = mi.opts;
    double gpver = gnuplot_version();
    double[2] xlim;
    double[2] ylim;
    GretlMatrixRaw * dims = null;
    const(char) * optlc = null;
    const(char) * sval;
    FILE * fp = null;
    char * datasrc = null;
    double linewidth = 1.0;
    double margin = 0.02;
    int display, have_payload = 0;
    int non_standard;
    int use_arg0 = 0;
    int do_key = 0;
    int height = 600;
    int border = 1;
    int notics = 1;
    int err = 0;

    if (mi.zrange != null) {
        have_payload = 1;
        if (mi.n_codes > 0 && gpver >= 5.4) {
            if (mi.zlabels != null || mi.zname != null) {
                do_key = 1;
            }
        }
    } else if (opts == null) {
        /* the simple outlines case */
        border = 1;
        notics = 0;
    }

    display = (mi.flags & MapFlags.display) != 0;
    non_standard = (mi.flags & MapFlags.non_std) != 0;

    set_map_plot_limits(mi, xlim.ptr, ylim.ptr, margin);

    if (gretl_bundle_has_key(opts, "height")) {
        height = cast(int) gretl_bundle_get_scalar(opts, "height", &err);
        if (display && height <= 0) {
            height = 600;
        }
    }

    gretl_push_c_numeric_locale();

    if (height > 0) {
        dims = geoplot_dimensions(xlim.ptr, ylim.ptr, height, have_payload,
                                  &non_standard);
    }
    if (display) {
        set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "display");
        if (mi.plotfile != null) {
            iact_gpfile = cast(char *) mi.plotfile;
        }
    } else if (mi.flags & MapFlags.is_image) {
        set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, mi.plotfile);
    }

    fp = open_plot_input_file(PlotType.geomap, cast(GptFlags) 0, &err);
    if (err) {
        return err;
    }

    fprintf(fp, "# geoplot %g %g\n".ptr, dims.val[0], dims.val[1]);

    if (!do_key) {
        fputs("unset key\n", fp);
    }

    if (have_payload) {
        err = print_map_palette(mi, gpver, fp);
        if (err) {
            fclose(fp);
            gretl_remove(gretl_plotfile());
            return err;
        }
    }

    fprintf(fp, "set xrange [%g:%g]\n".ptr, xlim[0], xlim[1]);
    fprintf(fp, "set yrange [%g:%g]\n".ptr, ylim[0], ylim[1]);

    if (gretl_bundle_has_key(opts, "title")) {
        sval = gretl_bundle_get_string(opts, "title", null);
        if (sval != null) {
            fprintf(fp, "set title '%s'\n".ptr, sval);
        }
    }

    if (gretl_bundle_get_int(opts, "tics", null)) {
        notics = 0;
    }
    if (notics) {
        fputs("set noxtics\n", fp);
        fputs("set noytics\n", fp);
    }

    if (gretl_bundle_get_int(opts, "logscale", null)) {
        fputs("set logscale cb\n", fp);
    }

    if (gretl_bundle_has_key(opts, "border")) {
        /* allow override of default */
        border = gretl_bundle_get_int(opts, "border", null);
    }

    if (border == 0) {
        fputs("unset border\n", fp);
    }

    if ((sval = gretl_bundle_get_string(opts, "literal", null)) != null) {
        fputs_literal(sval, fp);
    }

    if (gretl_bundle_has_key(opts, "linewidth")) {
        double lw = gretl_bundle_get_scalar(opts, "linewidth", &err);

        if (!err && lw >= 0) {
            if (have_payload) {
                linewidth = lw;
            } else if (lw >= 0.1) {
                linewidth = lw;
            }
        }
    }
    if (gretl_bundle_has_key(opts, "linecolor")) {
        sval = gretl_bundle_get_string(opts, "linecolor", &err);
        if (!err) {
            optlc = sval;
        }
    }

    gnuplot_missval_string(fp);

    if (map_do_inlining(mi)) {
        err = inline_map_data(mi.datfile, fp);
        if (!err) {
            datasrc = g_strdup("$MapData");
        }
    } else if (mi.flags & MapFlags.is_image) {
        /* @plotfile and @datfile are both disposable, no need
           to bother about name alignment */
        datasrc = g_strdup_printf("\"%s\"".ptr, mi.datfile);
    } else if (mi.plotfile != null) {
        /* the names of @plotfile and @datfile will already be
           correctly aligned */
        use_arg0 = 1;
    } else {
        /* rename @datfile to match the auto-named plot file */
        char * tmp = g_strdup_printf("%s.dat".ptr, gretl_plotfile());

        gretl_copy_file(mi.datfile, tmp);
        gretl_remove(mi.datfile);
        g_free(tmp);
        use_arg0 = 1;
    }

    if (use_arg0) {
        fputs("datafile = sprintf(\"%s.dat\", ARG0)\n".ptr, fp);
        datasrc = g_strdup("datafile");
    }

    if (!err) {
        output_map_plot_lines(mi, datasrc, optlc, have_payload,
                              do_key, linewidth, fp);
    }

    err = finalize_plot_input_file(fp);
    if (!err) {
        if (display && gretl_in_gui_mode()) {
            if (gretl_bundle_get_int(opts, "gui_auto", null)) {
                gretl_bundle_set_string(opts, "plotfile", gretl_plotfile());
                gretl_bundle_set_matrix(opts, "dims", dims);
            } else {
                manufacture_gui_callback(GretlCmdIndex.gnuplot);
            }
        }
    }

    gretl_pop_c_numeric_locale();

    gretl_matrix_free(dims);
    g_free(datasrc);

    return err;
}

/* Transcribe geoplot map file from @src to @dest, allowing for
   the possibility that data contained in a separate datafile
   have to be inlined. If @datname is NULL it is assumed that
   the datafile will be named as @src, with ".dat" appended.
*/

int transcribe_geoplot_file(const(char) * src, const(char) * dest, const(char) * datname)
{
    FILE * f1 = null, f2 = null;
    const(char) * mapdata = "$MapData";
    char[8196] buf;
    int integrate = -1;
    int n, err = 0;

    f1 = gretl_fopen(src, "rb");
    f2 = gretl_fopen(dest, "wb");

    if (f1 == null || f2 == null) {
        err = GretlError.fopen;
        goto bailout;
    }

    while (integrate < 0 && fgets(buf.ptr, buf.sizeof, f1)) {
        if (strstr(buf.ptr, mapdata)) {
            integrate = 0;
        } else if (!strncmp(buf.ptr, "datafile =", 10)) {
            integrate = 1;
        } else {
            fputs(buf.ptr, f2);
        }
    }

    if (integrate == 1) {
        /* open the datafile and inject its content */
        char * s, dattmp = null;
        FILE * fdat = null;
        int i;

        if (datname != null) {
            fdat = gretl_fopen(datname, "rb");
        } else {
            dattmp = g_strdup_printf("%s.dat".ptr, src);
            fdat = gretl_fopen(dattmp, "rb");
            g_free(dattmp);
        }

        if (fdat == null) {
            err = GretlError.fopen;
        } else {
            /* inject data */
            fprintf(f2, "%s << EOD\n".ptr, mapdata);
            while ((n = cast(int) fread(buf.ptr, 1, buf.sizeof, fdat)) > 0) {
                fwrite(buf.ptr, 1, n, f2);
            }
            fputs("EOD\n", f2);
            fclose(fdat);
            buf[0] = '\0';

            /* and transcribe the remainder of @src */
            while (fgets(buf.ptr, buf.sizeof, f1)) {
                if ((s = strstr(buf.ptr, "datafile")) != null) {
                    for (i=0; i<8; i++) {
                        s[i] = mapdata[i];
                    }
                }
                fputs(buf.ptr, f2);
            }
        }
    } else if (integrate == 0) {
        /* integration not needed */
        while ((n = cast(int) fread(buf.ptr, 1, buf.sizeof, f1)) > 0) {
            fwrite(buf.ptr, 1, n, f2);
        }
    } else {
        /* ?? */
        err = GretlError.data;
    }

 bailout:

    if (f1 != null) fclose(f1);
    if (f2 != null) fclose(f2);

    return err;
}

/* called from the interpolate plugin */

int write_tdisagg_plot(const(GretlMatrixRaw) * YY, int mult, const(char) * title, DatasetRaw * dset)
{
    const(double) * obs = null;
    char[16] mstr = 0;
    int t, T = YY.rows;
    double y0t;
    FILE * fp;
    int err = 0;

    set_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_U, "display");
    fp = open_plot_input_file(PlotType.regular, GptFlags.letterbox, &err);
    if (err) {
        return err;
    }

    if (dset != null) {
        fprintf(fp, "# timeseries %d (letterbox)\n".ptr, dset.pd);
        obs = gretl_plotx(dset, Gretlopt.OPT_NONE);
    } else {
        fputs("# timeseries 1 (letterbox)\n", fp);
    }
    fputs("set key left top\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (title != null) {
        fprintf(fp, "set title '%s'\n".ptr, title);
    }

    gretl_push_c_numeric_locale();

    if (obs != null) {
        double d1 = obs[dset.t1];
        double d2 = obs[dset.t2];

        fprintf(fp, "set xrange [%g:%g]\n".ptr, floor(d1), ceil(d2));
    }

    gnuplot_missval_string(fp);
    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (t=0; t<T; t++) {
        if (obs != null) {
            fprintf(fp, "%g ".ptr, obs[t+dset.t1]);
        } else {
            fprintf(fp, "%d ".ptr, t + 1);
        }
        y0t = gretl_matrix_get(YY, t, 0);
        if (na(y0t)) {
            fputs("? ", fp);
        } else {
            fprintf(fp, "%.10g ".ptr, y0t);
        }
        fprintf(fp, "%.10g\n".ptr, gretl_matrix_get(YY, t, 1));
    }
    fputs("EOD\n", fp);

    if (mult > 1) {
        sprintf(mstr.ptr, " * %d".ptr, mult);
    }

    fprintf(fp, "plot $data using 1:2 title \"%s\" w steps, \\\n".ptr,
            _("original data").ptr);
    fprintf(fp, " $data using 1:3 title \"%s%s\" w lines\n".ptr,
            _("final_ series").ptr, mstr.ptr);

    err = finalize_plot_input_file(fp);

    if (!err && gretl_in_gui_mode()) {
        manufacture_gui_callback(GretlCmdIndex.gnuplot);
    }

    gretl_pop_c_numeric_locale();

    return err;
}

