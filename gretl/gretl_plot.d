/* Port of lib/src/gretl_plot.c: the "plot" command block: collecting literal gnuplot lines and a data source
   (series, list or matrix) in a plot-in-progress, then finalizing it into a gnuplot script. The header
   gretl_plot.h has only prototypes. */

module gretl.gretl_plot;


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
import gretl.gretl_intl : gretl_pop_c_numeric_locale, gretl_push_c_numeric_locale;
import gretl.dataseth : dataset_is_time_series;
import gretl.genmain : current_series_index, generate_string;
import gretl.gmem : g_free;
import gretl.graphing : finalize_plot_input_file, gnuplot, matrix_plot, open_plot_input_file;
import gretl.graphingh : GptFlags, PlotType;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_list_new;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.libset : gretl_no_plots;
import gretl.options : push_option_param, valid_long_opt;
import gretl.optionsh : OptStatus;
import gretl.strutils : gretl_strdup, gretl_strndup, shift_string_left, strings_array_add, strings_array_free;
import gretl.gretl_array : get_array_by_name, gretl_array_get_strings, gretl_array_get_type;
import gretl.uservar : get_list_by_name, user_var_get_value_and_type;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


struct GretlPlot {
    int in_progress;
    Gretlopt opt;
    char * datasource;
    GretlType datatype;
    char ** lines;
    int nlines;
}

__gshared GretlPlot plot;

void clear_plot()
{
    plot.in_progress = 0;
    plot.opt = Gretlopt.OPT_NONE;
    free(plot.datasource);
    plot.datasource = null;
    plot.datatype = cast(GretlType) 0;
    strings_array_free(plot.lines, plot.nlines);
    plot.lines = null;
    plot.nlines = 0;
}

int plot_command_active()
{
    return plot.in_progress;
}

int no_data_plot(Gretlopt opt)
{
    FILE * fp = null;
    int i, np = 0;
    int err = 0;

    /* check that there's something to plot */
    for (i=0; i<plot.nlines; i++) {
    if (!strncmp(plot.lines[i], "plot ", 5) ||
        !strncmp(plot.lines[i], "splot ", 6)) {
        np++;
    }
    }

    if (np == 0) {
    gretl_errmsg_set("plot: nothing to plot");
    return GretlError.args;
    }

    fp = open_plot_input_file(PlotType.user, cast(GptFlags) 0, &err);

    if (!err) {
    for (i=0; i<plot.nlines; i++) {
        fputs(plot.lines[i], fp);
        fputc('\n', fp);
    }
    err = finalize_plot_input_file(fp);
    }

    return err;
}

/* In the following function we turn the array of "literal"
   lines from "plot" into the form

   { literal 1; literal 2; ...}

   which is the form wanted by the gnuplot() function at
   present. For future reference, it would be more efficient to
   revise gnuplot() so that it accepts an array of strings as
   input.
*/

char * construct_literal_arg(int * err)
{
    char * literal = null;
    size_t litlen = 0;
    int i;

    for (i=0; i<plot.nlines; i++) {
    litlen += strlen(plot.lines[i]);
    }

    if (litlen > 0) {
    litlen += 2 * plot.nlines + 4;
    literal = cast(typeof(literal)) calloc(litlen, 1);
    if (literal == null) {
        *err = GretlError.alloc;
    } else {
        strcpy(literal, "{ ");
        for (i=0; i<plot.nlines; i++) {
        strcat(literal, plot.lines[i]);
        strcat(literal, "; ");
        }
        strcat(literal, "}");
    }
    }

    return literal;
}

int execute_plot(const(DatasetRaw) * dset, Gretlopt opt)
{
    int * list = null;
    char * literal = null;
    int free_list = 0;
    int err = 0;

    plot.opt |= opt;

    if (plot.datatype == GretlType.list) {
    list = get_list_by_name(plot.datasource);
    } else if (plot.datatype == GretlType.series) {
    /* convert to singleton list */
    list = gretl_list_new(1);
    list[1] = current_series_index(dset, plot.datasource);
    free_list = 1;
    }

    if (plot.datatype != 0 && plot.nlines > 0) {
    literal = construct_literal_arg(&err);
    }

    if (!err && list != null && list[0] == 1) {
    /* default to time-series plot */
    if (dataset_is_time_series(dset)) {
        plot.opt |= Gretlopt.OPT_T;
    }
    }

    if (!err) {
    if (plot.datatype == 0) {
        err = no_data_plot(plot.opt);
    } else if (plot.datatype == GretlType.matrix) {
        GretlMatrixRaw * m = get_matrix_by_name(plot.datasource);

        err = matrix_plot(m, null, literal, plot.opt);
    } else {
        err = gnuplot(list, literal, dset, plot.opt);
    }
    }

    free(literal);

    if (free_list) {
    free(list);
    }

    return err;
}

/* we'll accept a single series, list or matrix as the
   "datasource" for a plot block */

int check_plot_data_source(const(char) * s, const(DatasetRaw) * dset)
{
    int err = 0;

    if (plot.datatype != 0) {
    /* duplicated data spec */
    err = GretlError.data;
    } else if (s == null) {
    /* no param given */
    return 0;
    } else {
    int id = current_series_index(dset, s);

    if (id >= 0) {
        plot.datatype = GretlType.series;
    } else {
        GretlType type = cast(GretlType) 0;

        user_var_get_value_and_type(s, &type);
        if (type == GretlType.matrix ||
        type == GretlType.list) {
        plot.datatype = type;
        }
    }
    }

    if (!err) {
    if (plot.datatype == 0) {
        err = GretlError.data;
    } else {
        plot.datasource = gretl_strdup(s);
    }
    }

    return err;
}

void maybe_unquote_param(char * s)
{
    if (*s == '"') {
    int n;

    shift_string_left(s, 1);
    n = cast(int) (strlen(s));
    if (s[n-1] == '"') {
        s[n-1] = '\0';
    }
    }
}

int check_plot_option(const(char) * s)
{
    const(char) * p = null;
    char * param = null;
    OptStatus status;
    Gretlopt opt;
    int err = 0;

    if (!strncmp(s, "--", 2)) {
    /* tolerate (non-required) option dashes */
    s += 2;
    }

    if ((p = strchr(s, '=')) != null) {
    char * flag = gretl_strndup(s, p - s);

    param = gretl_strdup(p + 1);
    opt = valid_long_opt(GretlCmdIndex.gnuplot, flag, &status);
    free(flag);
    } else {
    opt = valid_long_opt(GretlCmdIndex.gnuplot, s, &status);
    }

    if (opt == Gretlopt.OPT_NONE) {
    fprintf(stderr, "plot option: got OPT_NONE for '%s'\n".ptr, s);
    err = GretlError.badopt;
    } else if (status == OptStatus.no_parm && param != null) {
    fprintf(stderr, "plot option: got spurious param: '%s'\n".ptr, s);
    err = GretlError.badopt;
    } else if (status == OptStatus.needs_parm && param == null) {
    fprintf(stderr, "plot option: missing param: '%s'\n".ptr, s);
    err = GretlError.args;
    } else {
    plot.opt |= opt;
    if (param != null) {
        if (opt == Gretlopt.OPT_W) {
        /* --font=whatever */
        maybe_unquote_param(param);
        }
        err = push_option_param(GretlCmdIndex.gnuplot, opt, param);
        if (err) {
        fprintf(stderr, "plot option: error pushing param: '%s'\n".ptr, s);
        } else {
        param = null;
        }
    }
    }

    free(param);

    return err;
}

int plot_printf(const(char) * s, const(DatasetRaw) * dset)
{
    char * genline;
    char * genout;
    int err = 0;

    genline = g_strdup_printf("sprintf(%s)".ptr, s);
    gretl_push_c_numeric_locale();
    genout = generate_string(genline, cast(DatasetRaw *) dset, &err);
    gretl_pop_c_numeric_locale();

    if (err) {
    fprintf(stderr, "plot_printf: error %d: genline='%s'\n".ptr,
        err, genline);
    } else {
    strings_array_add(&plot.lines, &plot.nlines, genout);
    free(genout);
    }

    g_free(genline);

    return err;
}

/* This could go into strutils.c? */

int strings_array_splice(char *** pS1, int * pn1, char ** S2, int n2)
{
    char ** S = null;
    int i, j, n;
    int err = 0;

    if (S2 == null || n2 == 0) {
    return 0;
    }

    n = *pn1 + n2;

    S = cast(typeof(S)) realloc(*pS1, n * (*S).sizeof);
    if (S == null) {
    return GretlError.alloc;
    }

    for (i=0; i<n2; i++) {
    S[*pn1 + i] = null;
    }

    for (i=0; i<n2 && !err; i++) {
    j = *pn1 + i;
    if (S2[i] != null) {
        S[j] = gretl_strdup(S2[i]);
        if (S[j] == null) {
        err = GretlError.alloc;
        }
    }
    }

    *pS1 = S;
    *pn1 = n;

    return err;
}

int plot_append_strings(char *** pS, int * pns, const(char) * s)
{
    GretlArray * A = get_array_by_name(s);
    int ns = 0, err = 0;

    if (A == null) {
    err = GretlError.data;
    } else if (gretl_array_get_type(A) != GretlType.strings) {
    err = GretlError.types;
    } else {
    char ** S = gretl_array_get_strings(A, &ns);

    if (S != null) {
        err = strings_array_splice(pS, pns, S, ns);
    }
    }

    return err;
}

enum PlotField {
    last_field = 1,
    empty_ok
}

const(char) * get_plot_field_and_advance(const(char) * s, char * field, size_t maxlen, int flag, int * err)
{
    const(char) * p = s;
    int quoted = 0;
    size_t i = 0;

    while (isspace(*s)) s++;

    *field = '\0';

    while (*s) {
    if (*s == '"') {
        quoted = !quoted;
    } else if (!quoted && isspace(*s)) {
        break;
    }
    if (i < maxlen) {
        field[i++] = *s;
    } else {
        fprintf(stderr, "plot field: overflow!\n".ptr);
        *field = '\0';
        return p;
    }
    s++;
    }

    field[i] = '\0';
    s += cast(int) strspn(s, " \t\r\n");

    if (*field == '\0' && flag != PlotField.empty_ok) {
    *err = GretlError.args;
    } else if (flag == PlotField.last_field && *s != '\0') {
    gretl_errmsg_sprintf(_("Parse error at unexpected token '%s'").ptr, s);
    *err = GretlError.parse;
    }

    return s;
}

/**
 * gretl_plot_append_line:
 * @s: the line to append.
 * @dset: pointer to dataset.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_append_line(const(char) * s, const(DatasetRaw) * dset)
{
    char[256] field;
    int err = 0;

    if (!plot.in_progress) {
    return GretlError.data;
    }

    s = get_plot_field_and_advance(s, field.ptr, 16, 0, &err);

    if (!strcmp(field.ptr, "option")) {
    s = get_plot_field_and_advance(s, field.ptr, field.sizeof, PlotField.last_field, &err);
    if (!err) {
        err = check_plot_option(field.ptr);
        if (err) {
        fprintf(stderr, "Invalid plot option '%s'\n".ptr, field.ptr);
        }
    }
    } else if (!strcmp(field.ptr, "options")) {
    int flag = 0;

    while (1) {
        s = get_plot_field_and_advance(s, field.ptr, field.sizeof, flag, &err);
        if (err || *field.ptr == '\0') {
        break;
        } else {
        err = check_plot_option(field.ptr);
        }
        flag = PlotField.empty_ok;
    }
    } else if (!strcmp(field.ptr, "literal")) {
    err = strings_array_add(&plot.lines, &plot.nlines, s);
    } else if (!strcmp(field.ptr, "printf")) {
    err = plot_printf(s, dset);
    } else if (!strcmp(field.ptr, "strings")) {
    err = plot_append_strings(&plot.lines, &plot.nlines, s);
    } else {
    fprintf(stderr, "plot: invalid field '%s'\n".ptr, field.ptr);
    err = GretlError.parse;
    }

    if (err) {
    clear_plot();
    }

    return err;
}

/**
 * gretl_plot_start:
 * @param: name of the variable supplying the data.
 * @dset: pointer to dataset.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_start(const(char) * param, const(DatasetRaw) * dset)
{
    int err = 0;

    if (plot.in_progress) {
    clear_plot();
    err = GretlError.data;
    } else {
    err = check_plot_data_source(param, dset);
    }

    if (!err) {
    plot.in_progress = 1;
    }

    return err;
}

/**
 * gretl_plot_finalize:
 * @s: unused.
 * @dset: pointer to dataset.
 * @opt: may contain OPT_U (output spec).
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_plot_finalize(const(char) * s, const(DatasetRaw) * dset, Gretlopt opt)
{
    int err = 0;

    if (gretl_no_plots()) {
    ; /* skip it */
    } else {
    err = execute_plot(dset, opt);
    }
    clear_plot();

    return err;
}

