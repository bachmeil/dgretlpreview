/* Port of lib/src/plotspec.c: the plot specification structure (GPT_SPEC) and its management: creating, copying
   and freeing specs, lines, labels and arrows, writing the gnuplot input file (data, styles, bars, fit lines)
   and reading it back. The header plotspec.h is plotspech.d. */

module gretl.plotspec;


import core.stdc.ctype : isspace;
import core.stdc.float_ : DBL_EPSILON;
import core.stdc.math : ceil, exp, fabs, floor, fmod, isinf, isnan, log, log10, pow, sqrt;
import core.stdc.stdio : FILE, fclose, fopen, fgets, fputc, fputs, rewind, fprintf, printf, puts, sprintf, snprintf, sscanf, stderr;
import core.stdc.stdlib : abs, atoi, calloc, free, getenv, malloc, qsort, realloc;
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
import gretl.calendar : gretl_strptime;
import gretl.gmem : g_free;
import gretl.graphingh : FitType, GptFlags, N_GP_LINETYPES, PlotType, TermType, frequency_plot_code, get_png_output, gretlRGB;
import gretl.gretl_errors : gretl_errmsg_set;
import gretl.gretl_list : gretl_list_copy;
import gretl.gretl_paths : gretl_fopen;
import gretl.gretl_utils : doubles_array_free, doubles_array_new;
import gretl.gspawn : GSpawnFlags;
import gretl.gstrfuncs : g_strdup_printf;
import gretl.plotspech : GpArrowFlags, GpJustCodes, GpKeySpec, GpLineFlags, GpLineStyle, GpRange, GpStyleSpec, GptArrow, GptLabel, GptLine, GptSpec, LT_AUTO, PlotBars;
import gretl.strutils : gretl_strdup, string_is_blank, strings_array_free;
import gretl.graphing : get_boxcolor, get_graph_color, get_shadecolor, gnuplot_missval_string, print_rgb_hash, write_gp_dataval, write_plot_bounding_box_request, write_plot_type_string;
import gretl.nonparam : loess_fit, sort_pairs_by_x;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


/* for handling "recession bars" or similar */

GptSpec * plotspec_new()
{
    GptSpec * spec;
    int i;

    spec = cast(typeof(spec)) malloc((*spec).sizeof);
    if (spec == null) {
    return null;
    }

    spec.n_datacols = 0;
    spec.heredata = 0;
    spec.heredata_cols = null;

    spec.lines = null;
    spec.n_lines = 0;

    spec.labels = null;
    spec.n_labels = 0;

    spec.arrows = null;
    spec.n_arrows = 0;

    spec.literal = null;
    spec.n_literal = 0;

    spec.multi_xtics = null;
    spec.n_xtics = 0;

    for (i=0; i<5; i++) {
    spec.titles[i] = null;
    }

    *spec.xvarname.ptr = '\0';
    *spec.yvarname.ptr = '\0';

    spec.xticstr = null;
    spec.x2ticstr = null;

    spec.xfmt[0] = 0;
    spec.xtics[0] = 0;
    spec.mxtics[0] = 0;
    spec.yfmt[0] = 0;
    spec.ytics[0] = 0;
    spec.fname[0] = 0;
    spec.keyspec = GpKeySpec.left | GpKeySpec.top;

    for (i=0; i<6; i++) {
    spec.range[i][0] = NADBL;
    spec.range[i][1] = NADBL;
    if (i < 3) {
        spec.logbase[i] = 0.0;
    }
    }

    spec.b_ols = null;
    spec.b_quad = null;
    spec.b_cub = null;
    spec.b_inv = null;
    spec.b_log = null;
    spec.b_linlog = null;

    spec.code = PlotType.regular;
    spec.flags = cast(GptFlags) 0;
    spec.fit = FitType.none;
    spec.fp = null;
    spec.data = null;
    spec.auxdata = null;
    spec.markers = null;
    spec.n_markers = 0;
    spec.scale = 1.0;
    spec.labeled = null;
    spec.ptr = null;
    spec.bars = null;
    spec.fontstr = null;
    spec.reglist = null;
    spec.nobs = 0;
    spec.okobs = 0;
    spec.pd = 0;
    spec.nbars = 0;
    spec.boxwidth = 0;
    spec.fillfrac = 0;
    spec.samples = 0;
    spec.border = GSpawnFlags.default_;
    spec.border_lc[0] = '\0';
    spec.bmargin = 0;

    spec.termtype = TermType.none;

    return spec;
}

void plotbars_free(PlotBars * bars)
{
    if (bars != null) {
    doubles_array_free(bars.dx, bars.n);
    free(bars);
    }
}

void free_plotspec_lines(GptSpec * spec)
{
    int i;

    for (i=0; i<spec.n_lines; i++) {
    free(spec.lines[i].ustr);
    free(spec.lines[i].mcols);
    }

    free(spec.lines);
}

void plotspec_destroy(GptSpec * spec)
{
    if (spec == null) {
    return;
    }

    if (spec.lines != null) {
    free_plotspec_lines(spec);
    }
    if (spec.heredata_cols != null) {
    free(spec.heredata_cols);
    }
    if (spec.labels != null) {
    free(spec.labels);
    }
    if (spec.arrows != null) {
    free(spec.arrows);
    }
    if (spec.data != null) {
    gretl_matrix_free(spec.data);
    }
    if (spec.auxdata != null) {
    gretl_matrix_free(spec.auxdata);
    }
    if (spec.reglist != null) {
    free(spec.reglist);
    }
    if (spec.literal != null) {
    strings_array_free(spec.literal, spec.n_literal);
    }
    if (spec.multi_xtics != null) {
    strings_array_free(spec.multi_xtics, spec.n_xtics);
    }
    if (spec.markers != null) {
    strings_array_free(spec.markers, spec.n_markers);
    }
    if (spec.labeled != null) {
    free(spec.labeled);
    }
    if (spec.bars != null) {
    plotbars_free(spec.bars);
    }
    if (spec.fontstr != null) {
    free(spec.fontstr);
    }
    if (spec.xticstr != null) {
    free(spec.xticstr);
    }
    if (spec.x2ticstr != null) {
    free(spec.x2ticstr);
    }
    gretl_matrix_free(spec.b_ols);
    gretl_matrix_free(spec.b_quad);
    gretl_matrix_free(spec.b_cub);
    gretl_matrix_free(spec.b_inv);
    gretl_matrix_free(spec.b_log);
    gretl_matrix_free(spec.b_linlog);

    free(spec);
}

void plotspec_get_xt1_xt2(GptSpec * spec, double * xt1, double * xt2)
{
    if (spec.data != null) {
    *xt1 = spec.data.val[0];
    *xt2 = spec.data.val[spec.nobs - 1];
    }
}

__gshared GpStyleSpec[11] style_specs = [
    GpStyleSpec(GpLineStyle.lines, "lines".ptr, "lines".ptr),
    GpStyleSpec(GpLineStyle.points, "points".ptr, "points".ptr),
    GpStyleSpec(GpLineStyle.linespoints, "linespoints".ptr, "lines/points".ptr),
    GpStyleSpec(GpLineStyle.impulses, "impulses".ptr, "impulses".ptr),
    GpStyleSpec(GpLineStyle.dots, "dots".ptr, "dots".ptr),
    GpStyleSpec(GpLineStyle.steps, "steps".ptr, "steps".ptr),
    GpStyleSpec(GpLineStyle.boxes, "boxes".ptr, "boxes".ptr),
    GpStyleSpec(GpLineStyle.errorbars, "errorbars".ptr, "error bars".ptr),
    GpStyleSpec(GpLineStyle.filledcurve, "filledcurve".ptr, "filled curve".ptr),
    GpStyleSpec(GpLineStyle.candlesticks, "candlesticks".ptr, "candlesticks".ptr),
    GpStyleSpec(0, null, null)
];

const(char) * gp_line_style_display_name(int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
    if (t == style_specs[i].id) {
        return style_specs[i].trname;
    }
    }

    return "lines".ptr;
}

const(char) * gp_line_style_name(int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
    if (t == style_specs[i].id) {
        return style_specs[i].name;
    }
    }

    return "lines";
}

int gp_style_index_from_name(const(char) * s)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
    if (!strcmp(s, style_specs[i].name)) {
        return style_specs[i].id;
    }
    }

    /* allowed abbreviations */

    if (!strcmp(s, "l")) {
    return GpLineStyle.lines;
    } else if (!strcmp(s, "p")) {
    return GpLineStyle.points;
    } else if (!strcmp(s, "lp")) {
    return GpLineStyle.linespoints;
    } else if (!strcmp(s, "i")) {
    return GpLineStyle.impulses;
    } else if (!strcmp(s, "b")) {
    return GpLineStyle.boxes;
    } else if (!strncmp(s, "can", 3)) {
    return GpLineStyle.candlesticks;
    }

    /* fallback */
    return GpLineStyle.auto_;
}

int gp_style_index_from_display_name(const(char) * s)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
    if (!strcmp(s, _(style_specs[i].trname))) {
        return style_specs[i].id;
    }
    }

    return GpLineStyle.lines;
}

GpStyleSpec * get_style_spec(int t)
{
    int i;

    for (i=0; style_specs[i].id != 0; i++) {
    if (t == style_specs[i].id) {
        return &style_specs[i];
    }
    }

    return null;
}

char * gp_keyspec_string(GpKeySpec t)
{
    char[64] buf;

    buf[0] = '\0';

    if (t & GpKeySpec.left) {
        strcat(buf.ptr, " left");
    } else if (t & GpKeySpec.right) {
        strcat(buf.ptr, " right");
    }
    if (t & GpKeySpec.top) {
        strcat(buf.ptr, " top");
    } else if (t & GpKeySpec.bottom) {
        strcat(buf.ptr, " bottom");
    } else if (t & GpKeySpec.outside) {
        strcat(buf.ptr, " outside");
    }
    if (t & GpKeySpec.plain) {
        strcat(buf.ptr, " noenhanced");
    }

    return g_strdup(buf.ptr);   /* the C builds this in a GString */
}

void print_keypos_string(GpKeySpec t, FILE * fp)
{
    if (t == GpKeySpec.none) {
        fputs("set nokey\n", fp);
    } else {
        char * s = gp_keyspec_string(t);

        fprintf(fp, "set key %s\n".ptr, s);
        g_free(s);
    }
}

GpKeySpec gp_keyspec_from_string(const(char) * s)
{
    GpKeySpec ret = GpKeySpec.none;

    if (strstr(s, "left")) {
        ret |= GpKeySpec.left;
    } else if (strstr(s, "right")) {
        ret |= GpKeySpec.right;
    }
    if (strstr(s, "top")) {
        ret |= GpKeySpec.top;
    } else if (strstr(s, "bottom")) {
        ret |= GpKeySpec.bottom;
    } else if (strstr(s, "outside")) {
        ret |= GpKeySpec.outside;
    }
    if (strstr(s, "noenhanced")) {
        ret |= GpKeySpec.plain;
    }

    return ret;
}

int plotspec_add_line(GptSpec * spec)
{
    GptLine *lines;
    int n = spec.n_lines;

    lines = cast(typeof(lines)) realloc(spec.lines, (n + 1) * (*lines).sizeof);
    if (lines == null) {
    return GretlError.alloc;
    }

    spec.lines = lines;
    spec.n_lines += 1;

    lines[n].varnum = 0;
    lines[n].style = GpLineStyle.none;
    lines[n].pscale = 1.0;
    lines[n].rgb[0] = '\0';
    lines[n].yaxis = 1;
    lines[n].type = LT_AUTO;
    lines[n].ptype = 0;
    lines[n].dtype = 0;
    lines[n].width = 1.0;
    lines[n].ncols = 0;
    lines[n].whiskwidth = 0;
    lines[n].flags = 0;

    lines[n].title = null;
    lines[n].formula = null;
    lines[n].ustr = null;
    lines[n].mcols = null;

    return 0;
}

void copy_line_content(GptLine * targ, GptLine * src)
{
    targ.varnum = src.varnum;
    targ.style = src.style;
    targ.pscale = src.pscale;
    strcpy(targ.rgb.ptr, src.rgb.ptr);
    targ.yaxis = src.yaxis;
    targ.type = src.type;
    targ.ptype = src.ptype;
    targ.dtype = src.dtype;
    targ.width = src.width;
    targ.ncols = src.ncols;
    targ.whiskwidth = src.whiskwidth;
    targ.flags = src.flags;

    if (src.title != null) {
    g_free(targ.title);
    targ.title = g_strdup(src.title);
    }
    if (src.formula != null) {
    g_free(targ.formula);
    targ.formula = g_strdup(src.formula);
    }
    if (src.ustr != null) {
    free(targ.ustr);
    targ.ustr = gretl_strdup(src.ustr);
    }
    if (src.mcols != null) {
    free(targ.mcols);
    targ.mcols = gretl_list_copy(src.mcols);
    }
}

int plotspec_delete_line(GptSpec * spec, int i)
{
    GptLine *lines = spec.lines;
    int j, n = spec.n_lines;

    if (i < 0 || i >= n) {
    return GretlError.data;
    }

    g_free(lines[i].title);
    g_free(lines[i].formula);
    free(lines[i].ustr);
    free(lines[i].mcols);

    /* move the remaining line contents down */
    for (j=i; j<n-1; j++) {
    copy_line_content(&lines[j], &lines[j+1]);
    }

    spec.n_lines -= 1;

    lines = cast(typeof(lines)) realloc(spec.lines, (n - 1) * (*lines).sizeof);
    if (lines == null) {
    return GretlError.alloc;
    }

    spec.lines = lines;

    return 0;
}

GptLine * plotspec_clone_lines(GptSpec * spec, int * err)
{
    GptLine *lines = null;
    int i;

    if (spec.n_lines == 0) {
    /* no-op */
    return null;
    }

    lines = cast(typeof(lines)) calloc(spec.n_lines, (*lines).sizeof);

    if (lines == null) {
    *err = GretlError.alloc;
    } else {
    for (i=0; i<spec.n_lines; i++) {
        lines[i].ustr = null;
        lines[i].mcols = null;
        copy_line_content(&lines[i], &spec.lines[i]);
    }
    }

    return lines;
}

int plotspec_max_line_width(GptSpec * spec)
{
    int i, lw = 0;

    for (i=0; i<spec.n_lines; i++) {
    if (spec.lines[i].width > lw) {
        lw = cast(int) spec.lines[i].width;
    }
    }

    return lw;
}

int plotspec_add_label(GptSpec * spec)
{
    GptLabel *labels;
    int n = spec.n_labels;

    labels = cast(typeof(labels)) realloc(spec.labels, (n + 1) * (*labels).sizeof);
    if (labels == null) {
    return GretlError.alloc;
    }

    spec.labels = labels;
    spec.n_labels += 1;

    labels[n].text[0] = '\0';
    labels[n].pos[0] = NADBL;
    labels[n].pos[1] = NADBL;
    labels[n].just = GpJustCodes.left;

    return 0;
}

void copy_label_content(GptLabel * targ, GptLabel * src)
{
    strcpy(targ.text.ptr, src.text.ptr);
    targ.pos[0] = src.pos[0];
    targ.pos[1] = src.pos[1];
    targ.just = src.just;
}

int plotspec_delete_label(GptSpec * spec, int i)
{
    GptLabel *labels = spec.labels;
    int j, n = spec.n_labels;
    int err = 0;

    if (i < 0 || i >= n) {
    return GretlError.data;
    }

    for (j=i; j<n-1; j++) {
    copy_label_content(&labels[j], &labels[j+1]);
    }

    spec.n_labels -= 1;

    if (spec.n_labels == 0) {
    free(spec.labels);
    spec.labels = null;
    } else {
    labels = cast(typeof(labels)) realloc(spec.labels, (n - 1) * (*labels).sizeof);
    if (labels == null) {
        err = GretlError.alloc;
    } else {
        spec.labels = labels;
    }
    }

    return err;
}

GptLabel * plotspec_clone_labels(GptSpec * spec, int * err)
{
    GptLabel *labels = null;
    int i;

    if (spec.n_labels == 0) {
    /* no-op */
    return null;
    }

    labels = cast(typeof(labels)) malloc(spec.n_labels * (*labels).sizeof);
    if (labels == null) {
    *err = GretlError.alloc;
    } else {
    for (i=0; i<spec.n_labels; i++) {
        copy_label_content(&labels[i], &spec.labels[i]);
    }
    }

    return labels;
}

int plotspec_add_arrow(GptSpec * spec)
{
    GptArrow *arrows;
    int n = spec.n_arrows;

    arrows = cast(typeof(arrows)) realloc(spec.arrows, (n + 1) * (*arrows).sizeof);
    if (arrows == null) {
    return GretlError.alloc;
    }

    spec.arrows = arrows;
    spec.n_arrows += 1;
    arrows[n].x0 = 0;
    arrows[n].y0 = 0;
    arrows[n].x1 = 0;
    arrows[n].y1 = 0;
    arrows[n].flags = 0;

    return 0;
}

void copy_arrow_content(GptArrow * targ, GptArrow * src)
{
    targ.x0 = src.x0;
    targ.y0 = src.y0;
    targ.x1 = src.x1;
    targ.y1 = src.y1;
    targ.flags = src.flags;
}

int plotspec_delete_arrow(GptSpec * spec, int i)
{
    GptArrow *arrows = spec.arrows;
    int j, n = spec.n_arrows;
    int err = 0;

    if (i < 0 || i >= n) {
    return GretlError.data;
    }

    for (j=i; j<n-1; j++) {
    copy_arrow_content(&arrows[j], &arrows[j+1]);
    }

    spec.n_arrows -= 1;

    if (spec.n_arrows == 0) {
    free(spec.arrows);
    spec.arrows = null;
    } else {
    arrows = cast(typeof(arrows)) realloc(spec.arrows, (n - 1) * (*arrows).sizeof);
    if (arrows == null) {
        err = GretlError.alloc;
    } else {
        spec.arrows = arrows;
    }
    }

    return err;
}

GptArrow * plotspec_clone_arrows(GptSpec * spec, int * err)
{
    GptArrow *arrows = null;
    int i;

    if (spec.n_arrows == 0) {
    /* no-op */
    return null;
    }

    arrows = cast(typeof(arrows)) malloc(spec.n_arrows * (*arrows).sizeof);
    if (arrows == null) {
    *err = GretlError.alloc;
    } else {
    for (i=0; i<spec.n_arrows; i++) {
        copy_arrow_content(&arrows[i], &spec.arrows[i]);
    }
    }

    return arrows;
}

char * escape_quotes(const(char) * s)
{
    if (strchr(s, '"') == null) {
    return null;
    } else {
    int qcount = 0;
    char * ret, r;
    const(char) * p = s;

    while (*p) {
        if (*p == '"') qcount++;
        p++;
    }

    ret = cast(typeof(ret)) malloc(strlen(s) + 1 + qcount);
    if (ret == null) {
        return null;
    }

    r = ret;
    while (*s) {
        if (*s == '"') {
        *r++ = '\\';
        *r++ = '"';
        } else {
        *r++ = *s;
        }
        s++;
    }
    *r = 0;

    return ret;
    }
}

const(char) * gp_justification_string(int j)
{
    if (j == GpJustCodes.left) {
    return "left".ptr;
    } else if (j == GpJustCodes.center) {
    return "center".ptr;
    } else if (j == GpJustCodes.right) {
    return "right".ptr;
    } else {
    return "left".ptr;
    }
}

void print_plot_labelspec(GptLabel * lbl, FILE * fp)
{
    char * label = escape_quotes(lbl.text.ptr);

    fprintf(fp, "set label \"%s\" ".ptr, (label != null)?
        label : lbl.text.ptr);

    gretl_push_c_numeric_locale();

    fprintf(fp, "at %g,%g %s front\n".ptr,
        lbl.pos[0], lbl.pos[1],
        gp_justification_string(lbl.just));

    gretl_pop_c_numeric_locale();

    if (label != null) {
    free(label);
    }
}

void print_plot_arrow(const(GptArrow) * arrow, FILE * fp)
{
    gretl_push_c_numeric_locale();
    fprintf(fp, "set arrow from %g,%g to %g,%g %s".ptr,
        arrow.x0, arrow.y0, arrow.x1, arrow.y1,
        (arrow.flags & GpArrowFlags.head)? "head".ptr : "nohead".ptr);
    if (arrow.flags & GpArrowFlags.dots) {
    fputs(" lt 0", fp);
    }
    fputc('\n', fp);
    gretl_pop_c_numeric_locale();
}

int print_polar_labels(GptSpec * spec, FILE * fp)
{
    double x, y;
    int t;

    fputs("# printing data labels\n", fp);

    gretl_push_c_numeric_locale();

    for (t=0; t<spec.nobs; t++) {
    if (spec.labeled != null && !spec.labeled[t]) {
        /* printing only specified labels */
        continue;
    }
    if (sscanf(spec.markers[t], "%lf,%lf".ptr, &x, &y) == 2) {
        fprintf(fp, "set label \"%.3f,%.3f\" at %.3f,%.3f\n".ptr,
            x, y, x, y + .04);
    }
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

int usable_obs(const(double) * x, const(double) * y0, const(double) * y1, int t, const(double) ** py)
{
    *py = y0;

    if (!na(x[t]) && !na(y0[t])) {
    return 1;
    } else if (!na(x[t]) && y1 != null && !na(y1[t])) {
    *py = y1;
    return 1;
    } else {
    return 0;
    }
}

int plotspec_line_is_formula(GptSpec * spec, int i)
{
    int ret = 0;

    if (i >= 0 && i < spec.n_lines) {
    GptLine *line = &spec.lines[i];

    if (!string_is_blank(line.formula)) {
        ret = 1;
    } else if (line.flags & GpLineFlags.user) {
        ret = 1;
    } else if (i == 1 && (spec.flags & GptFlags.auto_fit)) {
        ret = 1;
    }
    }

    return ret;
}

int print_data_labels(GptSpec * spec, FILE * fp)
{
    const(double) * x, y, y0;
    const(double) * y1 = null;
    double xrange, yrange;
    double xmin, xmax;
    double ymin, ymax;
    double yoff;
    int nlines;
    int i, t;

    if (spec.markers == null) {
    /* "can't happen" */
    return 1;
    }

    if (spec.flags & GptFlags.polar) {
    return print_polar_labels(spec, fp);
    }

    nlines = 0;
    for (i=0; i<spec.n_lines; i++) {
    if (!plotspec_line_is_formula(spec, i)) {
        nlines++;
    }
    }
    if (nlines > 2 || spec.lines[0].ncols != 2) {
    return 1;
    }

    x = spec.data.val;
    y0 = x + spec.nobs;

    if (spec.code == PlotType.factorized) {
    y1 = y0 + spec.nobs;
    }

    xmin = spec.range[GpRange.x_range][0];
    xmax = spec.range[GpRange.x_range][1];
    xrange = (na(xmin) || na(xmax)) ? 0.0 : xmax - xmin;

    ymin = spec.range[GpRange.y_range][0];
    ymax = spec.range[GpRange.y_range][1];
    yrange = (na(ymin) || na(ymax)) ? 0.0 : ymax - ymin;

    if (xrange == 0.0 || yrange == 0.0) {
    ymin = 1.0e+50; ymax = -1.0e+50;
    xmin = 1.0e+50; xmax = -1.0e+50;

    for (t=0; t<spec.nobs; t++) {
        if (usable_obs(x, y0, y1, t, &y)) {
        if (yrange == 0.0) {
            if (y[t] < ymin) {
            ymin = y[t];
            }
            if (y[t] > ymax) {
            ymax = y[t];
            }
        }
        if (xrange == 0.0) {
            if (x[t] < xmin) {
            xmin = x[t];
            }
            if (x[t] > xmax) {
            xmax = x[t];
            }
        }
        }
    }
    if (yrange == 0.0) {
        yrange = ymax - ymin;
    }
    if (xrange == 0.0) {
        xrange = xmax - xmin;
    }
    }

    yoff = 0.03 * yrange;

    fputs("# printing data labels\n", fp);

    gretl_push_c_numeric_locale();

    for (t=0; t<spec.nobs; t++) {
    double xoff = 0.0;

    if (spec.labeled != null && !spec.labeled[t]) {
        /* printing only specified labels */
        continue;
    }

    if (usable_obs(x, y0, y1, t, &y)) {
        if (x[t] > .90 * xrange) {
        xoff = -.02 * xrange;
        }
        fprintf(fp, "set label \"%s\" at %.8g,%.8g\n".ptr, spec.markers[t],
            x[t] + xoff, y[t] + yoff);
    }
    }

    gretl_pop_c_numeric_locale();

    return 0;
}

void print_auto_fit_string(FitType fit, FILE * fp)
{
    if (fit == FitType.ols) {
    fputs("# plot includes automatic fit: OLS\n", fp);
    } else if (fit == FitType.quadratic) {
    fputs("# plot includes automatic fit: quadratic\n", fp);
    } else if (fit == FitType.cubic) {
    fputs("# plot includes automatic fit: cubic\n", fp);
    } else if (fit == FitType.inverse) {
    fputs("# plot includes automatic fit: inverse\n", fp);
    } else if (fit == FitType.loess) {
    fputs("# plot includes automatic fit: loess\n", fp);
    } else if (fit == FitType.loglin) {
    fputs("# plot includes automatic fit: semilog\n", fp);
    } else if (fit == FitType.linlog) {
    fputs("# plot includes automatic fit: linlog\n", fp);
    }
}

int range_is_unset(GptSpec * spec, int i)
{
    const(double) * r = spec.range[i].ptr;

    return r[0] == 0 && r[1] == 0;
}

void print_plot_ranges_etc(GptSpec * spec, FILE * fp)
{
    static const(char) *[6] rstrs = ["x".ptr, "y".ptr, "y2".ptr, "t".ptr, "x2".ptr, "cb".ptr];
    int i;

    gretl_push_c_numeric_locale();

    for (i=0; i<6; i++) {
    if (i < 3 && spec.logbase[i] > 0.0) {
        fprintf(fp, "set logscale %s %g\n".ptr, rstrs[i], spec.logbase[i]);
    }
    if (na(spec.range[i][0]) || na(spec.range[i][1]) ||
        spec.range[i][0] == spec.range[i][1]) {
        continue;
    } else if ((i == GpRange.y2_range && !(spec.flags & GptFlags.y2axis)) ||
        (i == GpRange.t_range && !(spec.flags & GptFlags.parametric))) {
        continue;
    } else if (i == GpRange.cb_range && range_is_unset(spec, i)) {
            continue;
        }

    fprintf(fp, "set %srange [%.10g:%.10g]\n".ptr, rstrs[i],
        spec.range[i][0], spec.range[i][1]);

    if (i == 4 && spec.code == PlotType.prob_dist && spec.samples == 0) {
        int ns = cast(int) (spec.range[i][1] - spec.range[i][0] + 1);

        fprintf(fp, "set samples %d\n".ptr, ns);
    }
    }

    if (spec.boxwidth > 0 && spec.boxwidth < 1) {
    fprintf(fp, "set boxwidth %g relative\n".ptr, cast(double) spec.boxwidth);
    } else if (spec.boxwidth < 0 && spec.boxwidth > -1) {
    fprintf(fp, "set boxwidth %g absolute\n".ptr, cast(double) -spec.boxwidth);
    }

    gretl_pop_c_numeric_locale();
}

int blank_user_line(GptLine * line)
{
    return ((line.flags & GpLineFlags.user) &&
        string_is_blank(line.formula));
}

void print_user_lines_info(GptSpec * spec, FILE * fp)
{
    int i, n = 0;

    for (i=0; i<spec.n_lines; i++) {
    if (spec.lines[i].flags & GpLineFlags.user) {
        n++;
    }
    }

    if (n > 0) {
    fprintf(fp, "# %d user-defined lines: ".ptr, n);
    for (i=0; i<spec.n_lines; i++) {
        if (spec.lines[i].flags & GpLineFlags.user) {
        fprintf(fp, "%d ".ptr, i);
        }
    }
    fputc('\n', fp);
    }
}

int more_lines(GptSpec * spec, int i, int skipline)
{
    if (i == spec.n_lines - 1) {
    /* at last line */
    return 0;
    } else {
    int j;

    for (j=i+1; j<spec.n_lines; j++) {
        if (j != skipline && !blank_user_line(&spec.lines[j])) {
        return 1;
        }
    }
    return 0;
    }
}

void print_linestyle(GptSpec * spec, int targ, int src, FILE * fp)
{
    GptLine *line;
    int done = 0;

    fprintf(fp, "set linetype %d ".ptr, targ + 1);

    if (src < spec.n_lines) {
    line = &spec.lines[src];
    if (*line.rgb.ptr != '\0' && line.style != GpLineStyle.filledcurve) {
        fprintf(fp, "lc rgb \"%s\"\n".ptr, line.rgb.ptr);
        done = 1;
    }
    }

    if (!done) {
    /* hmm, should we be doing this? */
    gretlRGB color = get_graph_color(targ);
    char[10] cstr;

    print_rgb_hash(cstr.ptr, color);
    fprintf(fp, "lc rgb \"%s\"\n".ptr, cstr.ptr);
    }
}

void write_styles_from_plotspec(GptSpec * spec, FILE * fp)
{
    int i;

    if (frequency_plot_code(spec.code)) {
    gretlRGB color = get_boxcolor();
    char[10] cstr;

    print_rgb_hash(cstr.ptr, color);
    fprintf(fp, "set linetype 1 lc rgb \"%s\"\n".ptr, cstr.ptr);
    fputs("set linetype 2 lc rgb \"#000000\"\n", fp);
    } else if (spec.code == PlotType.rq_tau) {
    fputs("set linetype 1 lc rgb \"#000000\"\n", fp);
    for (i=1; i<N_GP_LINETYPES; i++) {
        print_linestyle(spec, i, i, fp);
    }
    } else if (spec.code == PlotType.boxplots) {
    /* the first two elements use the same color */
    print_linestyle(spec, 0, 0, fp);
    print_linestyle(spec, 1, 2, fp);
    } else {
    for (i=0; i<N_GP_LINETYPES; i++) {
        print_linestyle(spec, i, i, fp);
    }
    }
}

void maybe_print_point_info(GptLine * line, FILE * fp)
{
    if (line.style == GpLineStyle.points ||
    line.style == GpLineStyle.linespoints) {
    if (line.ptype != 0) {
        fprintf(fp, " pt %d".ptr, line.ptype);
    }
    if (line.pscale != 1.0) {
        fprintf(fp, " ps %g".ptr, cast(double) line.pscale);
    }
    }

    if (line.style == GpLineStyle.lines ||
    line.style == GpLineStyle.linespoints) {
    if (line.dtype > 1) {
        fprintf(fp, " dt %d".ptr, line.dtype);
    }
    }
}

bool show_fit(GptSpec * s)
{
    return (s.fit == FitType.ols || s.fit == FitType.quadratic || s.fit == FitType.cubic || s.fit == FitType.inverse || s.fit == FitType.loess || s.fit == FitType.loglin || s.fit == FitType.linlog);
}

int gp_line_data_columns(GptSpec * spec, int i)
{
    if (spec.lines[i].flags & GpLineFlags.auxdata) {
    return 0;
    } else {
    return spec.lines[i].ncols;
    }
}

void plotspec_print_data(GptSpec * spec, int skipline, int * miss, FILE * fp)
{
    double *[5] x;
    int started_data_lines = 0;
    int i, j, t;

    for (i=0; i<spec.n_lines; i++) {
    int ncols = gp_line_data_columns(spec, i);

    if (ncols == 0 || i == skipline) {
        /* no (regular) data to print */
        continue;
    }

    if (!started_data_lines) {
        x[0] = spec.data.val;
        x[1] = x[0] + spec.nobs;
        started_data_lines = 1;
    }

    x[2] = x[1] + spec.nobs;
    x[3] = x[2] + spec.nobs;
    x[4] = x[3] + spec.nobs;

    for (t=0; t<spec.nobs; t++) {
        /* print x-axis value */
        if (na(x[0][t])) {
        write_gp_dataval(NADBL, fp, 0);
        *miss = 1;
        } else if (spec.flags & GptFlags.timefmt) {
        fprintf(fp, "%.0f ".ptr, x[0][t]);
        } else {
        fprintf(fp, "%.10g ".ptr, x[0][t]);
        }
        /* print y-axis value(s) */
        for (j=1; j<ncols; j++) {
        if (na(x[j][t])) {
            *miss = 1;
        }
        write_gp_dataval(x[j][t], fp, 0);
        }
        /* append markers if applicable */
        if (spec.markers != null && i == 0) {
        fprintf(fp, " # %s".ptr, spec.markers[t]);
        }
        fputc('\n', fp);
    }

    fputs("e\n", fp);

    if (spec.code == PlotType.factorized) {
        x[0] += ncols * spec.nobs;
        x[1] += ncols * spec.nobs;
    } else {
        x[1] += (ncols - 1) * spec.nobs;
    }
    }

    if (spec.auxdata != null) {
    int rows = gretl_matrix_rows(spec.auxdata);
    int cols = gretl_matrix_cols(spec.auxdata);

    fprintf(fp, "# auxdata %d %d\n".ptr, rows, cols);

    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
        fprintf(fp, "%.10g ".ptr, gretl_matrix_get(spec.auxdata, i, j));
        }
        fputc('\n', fp);
    }
    fputs("e\n", fp);
    }
}

void handle_line_width(GptLine * line, GptSpec * spec, FILE * fp)
{
    if (line.width == 1.0) {
        if (spec.scale >= 1.4) {
            line.width = 2.0;
        } else if (spec.scale > 1.0) {
            line.width = 1.5;
        }
    }
    fprintf(fp, " lw %g".ptr, cast(double) line.width);
}

void plotspec_print_heredata(GptSpec * spec, int * miss, FILE * fp)
{
    GretlMatrixRaw * m = null;
    double mti;
    int i, t;

    if (spec.auxdata != null) {
    int r = spec.auxdata.rows;
    int c = spec.auxdata.cols;
    int j;

    fprintf(fp, "# auxdata %d %d\n".ptr, r, c);
    fputs("$aux << EOA\n", fp);

    for (i=0; i<r; i++) {
        for (j=0; j<c; j++) {
        fprintf(fp, "%.10g ".ptr, gretl_matrix_get(spec.auxdata, i, j));
        }
        fputc('\n', fp);
    }
    fputs("EOA\n", fp);
    }

    if (spec.data == null || spec.n_datacols == 0) {
    return;
    }

    m = spec.data;

    fputs("$data << EOD\n", fp);

    for (t=0; t<spec.nobs; t++) {
    for (i=0; i<m.cols; i++) {
        mti = gretl_matrix_get(m, t, i);
        if (na(mti)) {
        write_gp_dataval(NADBL, fp, 0);
        *miss = 1;
        } else if (i == 0 && (spec.flags & GptFlags.timefmt)) {
        fprintf(fp, "%.0f ".ptr, mti);
        } else {
        fprintf(fp, "%.10g ".ptr, mti);
        }
    }
    if (spec.markers != null) {
        fprintf(fp, " # %s".ptr, spec.markers[t]);
    }
    fputc('\n', fp);
    }

    fputs("EOD\n", fp);
}

int plotspec_print(GptSpec * spec, FILE * fp)
{
    int png = get_png_output(spec);
    int mono = (spec.flags & GptFlags.mono);
    double xt1 = 0, xt2 = 0;
    char[16] src = "'-'";
    int skipline = -1;
    int ycol, any_y2 = 0;
    int anydata = 0;
    int i, k, miss = 0;

    if (spec.pd > 0) {
    fprintf(fp, "# timeseries %d".ptr, spec.pd);
    if (spec.flags & GptFlags.letterbox) {
        fputs(" (letterbox)\n", fp);
    } else {
        fputc('\n', fp);
    }
    }

    if (spec.scale != 1.0) {
    gretl_push_c_numeric_locale();
    fprintf(fp, "# scale = %.1f\n".ptr, spec.scale);
    gretl_pop_c_numeric_locale();
    }

    if (mono) {
    fputs("set mono\n", fp);
    } else {
    write_styles_from_plotspec(spec, fp);
    }

    if (!string_is_blank(spec.titles[0])) {
    fprintf(fp, "set title \"%s\"\n".ptr, spec.titles[0]);
    }

    if (!string_is_blank(spec.titles[1])) {
    fprintf(fp, "set xlabel \"%s\"\n".ptr, spec.titles[1]);
    }

    if (!string_is_blank(spec.titles[2])) {
    fprintf(fp, "set ylabel \"%s\"\n".ptr, spec.titles[2]);
    }

    if ((spec.flags & GptFlags.y2axis) && !string_is_blank(spec.titles[3])) {
    fprintf(fp, "set y2label \"%s\"\n".ptr, spec.titles[3]);
    }

    if (!string_is_blank(spec.titles[4])) {
    fprintf(fp, "set x2label \"%s\"\n".ptr, spec.titles[4]);
    }

    for (i=0; i<spec.n_labels; i++) {
    if (!string_is_blank(spec.labels[i].text.ptr)) {
        print_plot_labelspec(&spec.labels[i], fp);
    }
    }

    for (i=0; i<spec.n_arrows; i++) {
    print_plot_arrow(&spec.arrows[i], fp);
    }

    if (spec.flags & GptFlags.xzeroaxis) {
    fputs("set xzeroaxis\n", fp);
    }

    if (spec.flags & GptFlags.yzeroaxis) {
    fputs("set yzeroaxis\n", fp);
    }

    if (spec.flags & (GptFlags.grid_y | GptFlags.grid_x)) {
    if (!(spec.flags & GptFlags.grid_x)) {
        fputs("set grid ytics\n", fp);
    } else if (!(spec.flags & GptFlags.grid_y)) {
        fputs("set grid xtics\n", fp);
    } else {
        /* default, both */
        fputs("set grid\n", fp);
    }
    }

    if (spec.flags & GptFlags.parametric) {
    fputs("set parametric\n", fp);
    if (spec.samples > 0) {
        fprintf(fp, "set samples %d\n".ptr, spec.samples);
    }
    }

    if (spec.flags & GptFlags.polar) {
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    }

    gnuplot_missval_string(fp);

    print_keypos_string(spec.keyspec, fp);

    print_plot_ranges_etc(spec, fp);

    /* using time format for x-axis? */
    if (spec.flags & GptFlags.timefmt) {
    fputs("set xdata time\n", fp);
    fputs("set timefmt \"%s\"\n", fp);
    }

    /* special x and/or y format? */
    if (*spec.xfmt.ptr != '\0') {
    fprintf(fp, "set format x \"%s\"\n".ptr, spec.xfmt.ptr);
    }
    if (*spec.yfmt.ptr != '\0') {
    fprintf(fp, "set format y \"%s\"\n".ptr, spec.yfmt.ptr);
    }

    /* customized xtics? */
    if (!strcmp(spec.xtics.ptr, "none")) {
    fputs("set noxtics\n", fp);
    } else {
    if (!string_is_blank(spec.xtics.ptr)) {
        /* may contain "nomirror" */
        fprintf(fp, "set xtics %s\n".ptr, spec.xtics.ptr);
    }
    if (spec.xticstr != null) {
        fprintf(fp, "set xtics %s\n".ptr, spec.xticstr);
    }
    if (!string_is_blank(spec.mxtics.ptr)) {
        fprintf(fp, "set mxtics %s\n".ptr, spec.mxtics.ptr);
    }
    }

    /* using x2tics? */
    if (spec.x2ticstr != null) {
    fprintf(fp, "set x2tics %s\n".ptr, spec.x2ticstr);
    }

    /* customized ytics? */
    if (!strcmp(spec.ytics.ptr, "none")) {
    fputs("set noytics\n", fp);
    } else {
    if (!string_is_blank(spec.ytics.ptr)) {
        fprintf(fp, "set ytics %s\n".ptr, spec.ytics.ptr);
    }
    }

    if (spec.flags & GptFlags.y2axis) {
    /* using two y axes */
    fputs("set ytics nomirror\n", fp);
    fputs("set y2tics\n", fp);
    } else if (spec.border != GSpawnFlags.default_) {
    /* suppressing all or part of border */
    if (spec.border == 0) {
        fputs("unset border\n", fp);
    } else if (spec.border_lc[0] != '\0') {
        fprintf(fp, "set border %d lc rgb \"%s\"\n".ptr,
            spec.border, spec.border_lc.ptr);
    } else {
        fprintf(fp, "set border %d\n".ptr, spec.border);
    }
    if (string_is_blank(spec.xtics.ptr)) {
        fputs("set xtics nomirror\n", fp);
    }
    if (string_is_blank(spec.ytics.ptr)) {
        fputs("set ytics nomirror\n", fp);
    }
    }

    if (spec.bmargin > 0) {
    fprintf(fp, "set bmargin %d\n".ptr, spec.bmargin);
    }

    /* In case of plots that are editable (see gui client), it is
       important to write out the comment string that identifies the
       sort of graph, so that it will be recognized by type when
       it is redisplayed.
    */
    write_plot_type_string(spec.code, spec.flags, fp);

    if (spec.n_literal > 0) {
    fprintf(fp, "# literal lines = %d\n".ptr, spec.n_literal);
    for (i=0; i<spec.n_literal; i++) {
        if (spec.literal[i] != null && *spec.literal[i] != '\0') {
        fprintf(fp, "%s\n".ptr, spec.literal[i]);
        } else {
        fputs("# empty line!\n", fp);
        }
    }
    }

    if (spec.n_xtics > 0) {
    fprintf(fp, "# xtics lines = %d\n".ptr, spec.n_xtics);
    for (i=0; i<spec.n_xtics; i++) {
        if (spec.multi_xtics[i] != null && *spec.multi_xtics[i] != '\0') {
        fprintf(fp, "%s\n".ptr, spec.multi_xtics[i]);
        }
    }
    }

    if (show_fit(spec)) {
    print_auto_fit_string(spec.fit, fp);
    }

    if (spec.bars != null) {
    plotspec_get_xt1_xt2(spec, &xt1, &xt2);
    spec.nbars = n_bars_shown(xt1, xt2, spec.bars);
    if (spec.nbars > 0) {
        fprintf(fp, "# n_bars = %d\n".ptr, spec.nbars);
    }
    }

    print_user_lines_info(spec, fp);

    if ((spec.code == PlotType.freq_simple ||
     spec.code == PlotType.freq_normal ||
     spec.code == PlotType.freq_gamma)) {
    if (mono) {
        fputs("set style fill solid 0.3\n", fp);
    } else {
        fputs("set style fill solid 0.6\n", fp);
    }
    } else if (spec.fillfrac > 0) {
    gretl_push_c_numeric_locale();
    fprintf(fp, "set style fill solid %g\n".ptr, cast(double) spec.fillfrac);
    gretl_pop_c_numeric_locale();
    }

    if (spec.flags & GptFlags.print_markers) {
    print_data_labels(spec, fp);
    }

    if (spec.flags & GptFlags.fit_hidden) {
    skipline = 1;
    }

    for (i=0; i<spec.n_lines; i++) {
    if (i == skipline || blank_user_line(&spec.lines[i])) {
        continue;
    }
    if ((spec.flags & GptFlags.y2axis) && spec.lines[i].yaxis != 1) {
        any_y2 = 1;
        break;
    }
    }

    if ((spec.bars != null && spec.nbars > 0) ||
    spec.data != null || spec.auxdata != null) {
    anydata = 1;
    }

    if (spec.heredata && anydata) {
    fputs("# start inline data\n", fp);
    if (spec.bars != null && spec.nbars > 0) {
        print_plotbars(spec, fp);
    }
    if (spec.data != null || spec.auxdata != null) {
        plotspec_print_heredata(spec, &miss, fp);
    }
    fputs("# end inline data\n", fp);
    strcpy(src.ptr, "$data");
    }

    fputs("plot \\\n", fp);

    gretl_push_c_numeric_locale();

    if (spec.nbars > 0) {
    print_bars_header(spec, fp);
    }

    ycol = 2;

    for (i=0; i<spec.n_lines; i++) {
    GptLine *line = &spec.lines[i];

    if (i == skipline || blank_user_line(line)) {
        if (i < spec.n_lines - 1) {
        continue;
        } else {
        break;
        }
    }
    if (!string_is_blank(line.formula)) {
        fprintf(fp, "%s ".ptr, line.formula);
    } else if (line.ustr != null) {
        fprintf(fp, "%s using %s ".ptr, src.ptr, line.ustr);
    } else if (line.mcols != null) {
        int j;

        fprintf(fp, "%s using ".ptr, src.ptr);
        for (j=1; j<=line.mcols[0]; j++) {
        fprintf(fp, "%d".ptr, line.mcols[j]);
        fputc(j < line.mcols[0] ? ':' : ' ', fp);
        }
    } else {
        fprintf(fp, "%s using 1".ptr, src.ptr);
        if (line.ncols == 5) {
        /* Note: boxplot candlesticks, hard-wired! */
        fprintf(fp, ":%d:%d:%d:%d".ptr, ycol+1, ycol, ycol+3, ycol+2);
        } else if (line.ncols == 2) {
        if (line.flags & GpLineFlags.boxdata) {
            /* boxplot median, hard-wired */
            fprintf(fp, ":%d:%d:%d:%d".ptr, ycol, ycol, ycol, ycol);
        } else {
            fprintf(fp, ":%d".ptr, ycol);
        }
        } else {
        for (k=2; k<=line.ncols; k++) {
            fprintf(fp, ":%d".ptr, ycol + k - 2);
        }
        }
        fputc(' ', fp);
    }

    if ((spec.flags & GptFlags.y2axis) && line.yaxis != 1) {
        fprintf(fp, "axes x1y%d ".ptr, line.yaxis);
    }

    if (string_is_blank(line.title)) {
        fputs("notitle ", fp);
    } else {
        fprintf(fp, "title \"%s".ptr, line.title);
        if (any_y2 && strchr(line.title, '(') == null) {
        if (line.yaxis == 1) {
            fprintf(fp, " (%s)\" ".ptr, _("left").ptr);
        } else {
            fprintf(fp, " (%s)\" ".ptr, _("right").ptr);
        }
        } else {
        fputs("\" ", fp);
        }
    }

    if (line.style == GpLineStyle.filledcurve && line.type == LT_AUTO) {
        if (*line.rgb.ptr != '\0') {
        fprintf(fp, "lc rgb \"%s\" ".ptr, line.rgb.ptr);
        } else {
        print_filledcurve_color(fp);
        }
    }

    if (line.style != GpLineStyle.auto_) {
        fprintf(fp, "w %s".ptr, gp_line_style_name(line.style));
    }

    /* FIXME handling of linetype for PlotType.band? */
    if (spec.code != PlotType.band && line.type != LT_AUTO) {
        fprintf(fp, " lt %d".ptr, line.type);
    } else if (spec.nbars > 0 && line.style != GpLineStyle.filledcurve) {
        fprintf(fp, " lt %d".ptr, i + 1);
    }

    maybe_print_point_info(line, fp);

    if (line.style != GpLineStyle.filledcurve) {
            handle_line_width(line, spec, fp);
        }
    if (line.whiskwidth > 0) {
        fprintf(fp, " whiskerbars %g".ptr, cast(double) line.whiskwidth);
    }
    if (spec.heredata) {
        ycol += line.ncols - 1;
    }
    if (more_lines(spec, i, skipline)) {
        fputs(", \\", fp);
    }
    fputc('\n', fp);
    }

    /* supply the data to gnuplot inline, old style */

    if (!spec.heredata && anydata) {
    if (spec.bars != null && spec.nbars > 0) {
        print_plotbars(spec, fp);
    }
    if (spec.data != null || spec.auxdata != null) {
        plotspec_print_data(spec, skipline, &miss, fp);
    }
    }

    gretl_pop_c_numeric_locale();

    if (png) {
    write_plot_bounding_box_request(fp);
    }

    return miss;
}

int set_loess_fit(GptSpec * spec, int d, double q, GretlMatrixRaw * x, GretlMatrixRaw * y, GretlMatrixRaw * yh)
{
    int i, t, T = gretl_vector_get_length(y);
    const(double) * src;
    GretlMatrixRaw * m;

    m = gretl_matrix_alloc(T, 3);
    if (m == null) {
    return GretlError.data;
    }

    gretl_matrix_free(spec.data);
    spec.data = m;

    for (i=0; i<3; i++) {
    src = i == 0 ? x.val : i == 1 ? y.val : yh.val;
    for (t=0; t<T; t++) {
        gretl_matrix_set(m, t, i, src[t]);
    }
    }

    spec.data = m;
    spec.n_datacols = 3;
    spec.nobs = spec.okobs = T;

    g_free(spec.lines[1].title);
    spec.lines[1].title = g_strdup_printf(_("loess fit, d = %d, q = %g").ptr, d, q);
    spec.lines[1].pscale = 1.0;
    spec.lines[1].style = GpLineStyle.lines;
    spec.lines[1].ncols = 2;

    spec.fit = FitType.loess;

    return 0;
}

void set_plotfit_formula(GptLine * line, FitType f, const(double) * b, double x0, double pd)
{
    char * s = null;
    char[32] x;

    g_free(line.formula);
    line.formula = null;
    gretl_push_c_numeric_locale();

    /* If x0 is not NA that indicates that we've calculated a fit
       against time (t = 1,2,3,...,T), and we have to transform the
       equation so that it produces correct fitted values given
       x.ptr-data in the form year[.fraction].
    */
    if (!na(x0)) {
    if (pd > 1) {
        sprintf(x.ptr, "(%.1f*(x-%.8g)+1)".ptr, pd, x0);
    } else {
        sprintf(x.ptr, "(x-%.8g+1)".ptr, x0);
    }
    } else {
    strcpy(x.ptr, "x");
    }

    if (f == FitType.ols) {
    s = g_strdup_printf("%.8g + %.8g*%s".ptr, b[0], b[1], x.ptr);
    } else if (f == FitType.quadratic) {
    s = g_strdup_printf("%.8g + %.8g*%s + %.8g*%s**2".ptr, b[0], b[1], x.ptr, b[2], x.ptr);
    } else if (f == FitType.cubic) {
    s = g_strdup_printf("%.8g + %.8g*%s + %.8g*%s**2 + %.8g*%s**3".ptr,
                b[0], b[1], x.ptr, b[2], x.ptr, b[3], x.ptr);
    } else if (f == FitType.inverse) {
    s = g_strdup_printf("%.8g + %.8g/%s".ptr, b[0], b[1], x.ptr);
    } else if (f == FitType.loglin) {
    s = g_strdup_printf("exp(%.8g + %.8g*%s)".ptr, b[0], b[1], x.ptr);
    } else if (f == FitType.linlog) {
    if (!na(x0)) {
        s = g_strdup_printf("%.8g + %.8g*log%s".ptr, b[0], b[1], x.ptr);
    } else {
        s = g_strdup_printf("%.8g + %.8g*log(x)".ptr, b[0], b[1]);
    }
    }

    gretl_pop_c_numeric_locale();
    line.formula = s;
}

void set_plotfit_line(GptLine * line, FitType f, const(double) * b, double x0, double pd)
{
    char xc = (na(x0)) ? 'X' : 't';
    char * s = null;

    g_free(line.title);
    line.title = null;

    /* first compose the key string for the fitted line */
    if (f == FitType.ols) {
    s = g_strdup_printf("Y = %#.3g %c %#.3g%c".ptr, b[0],
                (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    } else if (f == FitType.quadratic) {
    s = g_strdup_printf("Y = %#.3g %c %#.3g%c %c %#.3g%c^2".ptr, b[0],
                (b[1] > 0)? '+' : '-', fabs(b[1]), xc,
                (b[2] > 0)? '+' : '-', fabs(b[2]), xc);
    } else if (f == FitType.cubic) {
    s = g_strdup_printf("Y = %#.3g %c %#.3g%c %c %#.3g%c^2 %c %#.3g%c^3".ptr,
                b[0], (b[1] > 0)? '+' : '-', fabs(b[1]), xc,
                (b[2] > 0)? '+' : '-', fabs(b[2]), xc,
                (b[3] > 0)? '+' : '-', fabs(b[3]), xc);
    } else if (f == FitType.inverse) {
    s = g_strdup_printf("Y = %#.3g %c %#.3g(1/%c)".ptr, b[0],
                (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    } else if (f == FitType.loglin) {
    s = g_strdup_printf("logY = %#.3g %c %#.3g%c".ptr, b[0],
                (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    if (xc == 't' && (pd == 1 || pd == 4 || pd == 12)) {
        /* display annual growth rate in title */
        double g = 100 * (pow(exp(b[1]), pd) - 1);
        char * gstr, tmp;

        gstr = g_strdup_printf("\\n(%s %.2f%%)".ptr, _("annual growth").ptr, g);
        tmp = g_strdup_printf("%s%s".ptr, s, gstr);
        g_free(s);
        s = tmp;
        g_free(gstr);
    }
    } else if (f == FitType.linlog) {
    s = g_strdup_printf("Y = %#.3g %c %#.3glog(%c)".ptr, b[0],
                (b[1] > 0)? '+' : '-', fabs(b[1]), xc);
    }

    line.title = s;

    /* then set the formula itself */
    set_plotfit_formula(line, f, b, x0, pd);
}

void plotspec_set_fitted_line(GptSpec * spec, FitType f, double x0)
{
    double pd = spec.pd;
    const(double) * b;

    if (f == FitType.ols) {
    b = spec.b_ols.val;
    } else if (f == FitType.quadratic) {
    b = spec.b_quad.val;
    } else if (f == FitType.cubic) {
    b = spec.b_cub.val;
    } else if (f == FitType.inverse) {
    b = spec.b_inv.val;
    } else if (f == FitType.loglin) {
    b = spec.b_log.val;
    } else if (f == FitType.linlog) {
    b = spec.b_linlog.val;
    } else {
    return;
    }

    set_plotfit_line(&spec.lines[1], f, b, x0, pd);

    spec.fit = f;
    spec.lines[1].pscale = 1.0;
    spec.lines[1].style = GpLineStyle.lines;
    spec.lines[1].ncols = 0;
}

int plotspec_add_fit(GptSpec * spec, FitType f)
{
    GretlMatrixRaw * y = null;
    GretlMatrixRaw * X = null;
    GretlMatrixRaw * b = null;
    GretlMatrixRaw * yh = null;
    const(double) * px = spec.data.val;
    const(double) * py;
    int T = spec.okobs;
    double x0 = NADBL;
    double xt, q = 0.5;
    int d = 1;
    int i, t, k;
    int err = 0;

    if ((spec.flags & GptFlags.ts) && f != FitType.loess) {
    if (spec.pd == 1 || spec.pd == 4 || spec.pd == 12) {
        /* Annual, quarterly or monthly data: we'll estimate
           the fit against t = 1,2,3,...,T to give a more
           comprehensible displayed equation. So we have
           to record the initial x-value (e.g. 2008.25) to
           get the fitted line right on the plot.
        */
        x0 = px[0];
    }
    }

    if ((f == FitType.ols && spec.b_ols != null) ||
    (f == FitType.quadratic && spec.b_quad != null) ||
    (f == FitType.cubic && spec.b_cub != null) ||
    (f == FitType.inverse && spec.b_inv != null) ||
    (f == FitType.loglin && spec.b_log != null) ||
    (f == FitType.linlog && spec.b_linlog != null)) {
    /* just activate existing setup */
    plotspec_set_fitted_line(spec, f, x0);
    return 0;
    }

    if (f == FitType.ols || f == FitType.inverse ||
    f == FitType.loglin || f == FitType.linlog) {
    k = 2;
    } else if (f == FitType.quadratic) {
    k = 3;
    } else if (f == FitType.cubic) {
    k = 4;
    } else if (f == FitType.loess) {
    k = 1;
    } else {
    return GretlError.data;
    }

    y = gretl_column_vector_alloc(T);
    X = gretl_matrix_alloc(T, k);
    if (y == null || X == null) {
    err = GretlError.alloc;
    goto bailout;
    }

    if (f != FitType.loess) {
    b = gretl_column_vector_alloc(k);
    if (b == null) {
        err = GretlError.alloc;
        goto bailout;
    }
    }

    px = spec.data.val;
    py = px + spec.nobs;

    i = 0;
    for (t=0; t<spec.nobs; t++) {
    if (!na(x0)) {
        /* use a 1-based time index */
        xt = t + 1;
    } else {
        xt = px[t];
    }
    if (!na(py[t]) && !na(xt)) {
        if (f == FitType.loglin) {
        /* y must be positive */
        if (py[t] > 0.0) {
            y.val[i] = log(py[t]);
        } else {
            gretl_errmsg_set(_("Non-positive values encountered").ptr);
            err = GretlError.data;
            goto bailout;
        }
        } else {
        y.val[i] = py[t];
        }
        if (f == FitType.loess) {
        gretl_matrix_set(X, i, 0, xt);
        } else {
        gretl_matrix_set(X, i, 0, 1.0);
        if (f == FitType.linlog) {
            /* x must be positive */
            if (xt > 0.0) {
            gretl_matrix_set(X, i, 1, log(xt));
            } else {
            gretl_errmsg_set(_("Non-positive values encountered").ptr);
            err = GretlError.data;
            goto bailout;
            }
        } else if (f == FitType.inverse) {
            gretl_matrix_set(X, i, 1, 1.0 / xt);
        } else {
            gretl_matrix_set(X, i, 1, xt);
        }
        if (f == FitType.quadratic || f == FitType.cubic) {
            gretl_matrix_set(X, i, 2, xt * xt);
        }
        if (f == FitType.cubic) {
            gretl_matrix_set(X, i, 3, xt * xt * xt);
        }
        }
        i++;
    }
    }

    if (f == FitType.loess) {
    err = sort_pairs_by_x(X, y, null, spec.markers);
    if (!err) {
        yh = loess_fit(X, y, d, q, Gretlopt.OPT_R, &err);
    }
    } else if (f == FitType.loglin) {
    double s2;

    err = gretl_matrix_ols(y, X, b, null, null, &s2);
    if (!err) {
        b.val[0] += s2 / 2;
    }
    } else {
    err = gretl_matrix_ols(y, X, b, null, null, null);
    }

    if (!err && spec.n_lines == 1) {
    err = plotspec_add_line(spec);
    }

    if (!err) {
    if (f == FitType.ols) {
        spec.b_ols = b;
        b = null;
    } else if (f == FitType.quadratic) {
        spec.b_quad = b;
        b = null;
    } else if (f == FitType.cubic) {
        spec.b_cub = b;
        b = null;
    } else if (f == FitType.inverse) {
        spec.b_inv = b;
        b = null;
    } else if (f == FitType.loglin) {
        spec.b_log = b;
        b = null;
    } else if (f == FitType.linlog) {
        spec.b_linlog = b;
        b = null;
    }
    }

    if (!err) {
    if (f == FitType.loess) {
        set_loess_fit(spec, d, q, X, y, yh);
        g_free(spec.lines[1].formula);
        spec.lines[1].formula = null;
    } else {
        plotspec_set_fitted_line(spec, f, x0);
    }
    }

 bailout:

    gretl_matrix_free(y);
    gretl_matrix_free(X);
    gretl_matrix_free(b);
    gretl_matrix_free(yh);

    return err;
}

/* next: apparatus for producing shaded bars in time-series plots */

PlotBars * plotbars_new(int n)
{
    PlotBars * bars = cast(PlotBars *) malloc(PlotBars.sizeof);

    if (bars != null) {
    bars.dx = doubles_array_new(n, 2);
    if (bars.dx == null) {
        free(bars);
        bars = null;
    } else {
        bars.n = n;
        bars.t1 = bars.t2 = 0;
        bars.ymin = bars.ymax = 0;
    }
    }

    return bars;
}

/* Heuristic for detecting the case where the x-axis of a time-series
   plot has unix time as its unit.
*/

int use_unix_time(double xmin, double xmax)
{
    if (xmax > 1.0e7 && xmax == floor(xmax) && xmin == floor(xmin)) {
    return 1;
    } else {
    return 0;
    }
}

/* Read and check a plain text "plotbars" file.  Such a file may contain
   comment lines starting with '#'; other than comments, each line
   should contain a pair of space-separated date strings in gretl
   monthly date format (YYYY:MM).  These strings represent,
   respectively, the start and end of an episode of some kind (paradigm
   case: peak and trough of business cycle).  Each pair may be used to
   define a vertical bar on a time-series plot. The file can contain any
   number of such pairs.  An example is provided in
   <prefix>/share/data/plotbars/nber.txt.

   The dates are converted internally into the format, year plus decimal
   fraction of year, and may be used when plotting annual, quarterly or
   monthly time series.
*/

PlotBars * parse_bars_file(const(char) * fname, double xmin, double xmax, int * err)
{
    FILE * fp;
    PlotBars * bars = null;
    char[128] line;
    int ncolon = 0;
    int ndash = 0;
    int nyrs = 0;
    int nother = 0;
    int ntypes = 0;
    int y1, y2;
    int m1, m2;

    fp = gretl_fopen(fname, "r");
    if (fp == null) {
    *err = GretlError.fopen;
    return null;
    }

    /* first count the data lines */
    while (fgets(line.ptr, line.sizeof, fp)) {
    if (*line.ptr == '#' || string_is_blank(line.ptr)) {
        continue;
    } else if (sscanf(line.ptr, "%d:%d %d:%d".ptr, &y1, &m1, &y2, &m2) == 4) {
        ncolon++;
    } else if (sscanf(line.ptr, "%d-%d %d-%d".ptr, &y1, &m1, &y2, &m2) == 4) {
        ndash++;
    } else if (sscanf(line.ptr, "%d %d".ptr, &y1, &y2) == 2) {
        nyrs++;
    } else {
        nother++;
    }
    }

    ntypes = (ncolon > 0) + (ndash > 0) + (nyrs > 0);

    /* initial check and allocation */
    if (nother > 0) {
    /* non-comment, non-blank, non-parseable cruft found */
    *err = GretlError.data;
    } else if (ntypes == 0 || ntypes > 1) {
    /* no parseable data found, or a mash-up */
    *err = GretlError.data;
    } else {
    bars = plotbars_new(ncolon + ndash + nyrs);
    if (bars == null) {
        *err = GretlError.alloc;
    }
    }

    if (*err == 0) {
    double x0, x1;
    int unix_time;
    int i = 0;

    rewind(fp);
    unix_time = use_unix_time(xmin, xmax);

    /* now read, check, convert and record the date pairs */
    while (fgets(line.ptr, line.sizeof, fp) && !*err) {
        if (*line.ptr == '#' || string_is_blank(line.ptr)) {
        continue;
        }
        if (nyrs) {
        sscanf(line.ptr, "%d %d".ptr, &y1, &y2);
        m1 = 1;
        m2 = 12;
        } else if (ncolon) {
        sscanf(line.ptr, "%d:%d %d:%d".ptr, &y1, &m1, &y2, &m2);
        } else {
        sscanf(line.ptr, "%d-%d %d-%d".ptr, &y1, &m1, &y2, &m2);
        }
        if (unix_time) {
        const(char) * prtfmt = "%d-%02d-%02d";
        const(char) * argfmt = "%Y-%m-%d";
        char[12] buf;

        sprintf(buf.ptr, prtfmt, y1, m1, 1);
        gretl_strptime(buf.ptr, argfmt, &x0);
        sprintf(buf.ptr, prtfmt, y2, m2, 28);
        gretl_strptime(buf.ptr, argfmt, &x1);
        } else {
        x0 = y1 + (m1 - 1.0) / 12;
        x1 = y2 + (nyrs ? m2 / 12.0 : (m2 - 1.0) / 12);
        }
        if (x1 < x0) {
        *err = GretlError.data;
        } else {
        bars.dx[i][0] = x0;
        bars.dx[i][1] = x1;
        i++;
        }
    }
    }

    fclose(fp);

    if (*err == GretlError.data) {
    gretl_errmsg_set(_("Dates file does not conform to the specification").ptr);
    }

    if (*err && bars != null) {
    plotbars_free(bars);
    bars = null;
    }

    return bars;
}

/* Output data representing vertical shaded bars: each bar is
   represented by a pair of rows.
*/

void print_plotbars(GptSpec * spec, FILE * fp)
{
    PlotBars * bars = spec.bars;
    double y0 = bars.ymin;
    double y1 = bars.ymax;
    int i, started, stopped;

    if (spec.heredata) {
    fputs("$bars << EOB\n", fp);
    }

    for (i=0; i<bars.n; i++) {
    if (bars.dx[i][1] < bars.t1) {
        continue;
    }
    if (bars.dx[i][0] > bars.t2) {
        break;
    }
    started = stopped = 0;
    if (bars.dx[i][0] >= bars.t1) {
        /* start is in range */
        fprintf(fp, "%.3f %.10g %.10g\n".ptr, bars.dx[i][0], y0, y1);
        started = 1;
    }
    if (bars.dx[i][1] <= bars.t2) {
        /* stop is in range */
        if (!started) {
        /* but the start was not */
        fprintf(fp, "%.3f %.10g %.10g\n".ptr, bars.t1, y0, y1);
        }
        fprintf(fp, "%.3f %.10g %.10g\n".ptr, bars.dx[i][1], y0, y1);
        if (!spec.heredata) {
        fputs("e\n", fp);
        }
        stopped = 1;
    }
    if (started && !stopped) {
        /* truncated */
        fprintf(fp, "%.3f %.10g %.10g\n".ptr, bars.t2, y0, y1);
        if (!spec.heredata) {
        fputs("e\n", fp);
        }
    }
    }

    if (spec.heredata) {
    fputs("EOB\n", fp);
    }
}

/* output special plot lines for creating vertical shaded bars
   in a time-series graph: the corresponding data will take
   the form x:ymin:ymax */

void print_bars_header(GptSpec * spec, FILE * fp)
{
    char[8] cstr;
    int i;

    print_rgb_hash(cstr.ptr, get_shadecolor());

    if (spec.heredata) {
    int j = 0;

    for (i=0; i<spec.nbars; i++) {
        fprintf(fp, "'$bars' every ::%d::%d using 1:2:3 notitle lc rgb \"%s\" w filledcurve, \\\n".ptr, j, j+1, cstr.ptr);
        j += 2;
    }
    } else {
    for (i=0; i<spec.nbars; i++) {
        fprintf(fp, "'-' using 1:2:3 notitle lc rgb \"%s\" w filledcurve, \\\n".ptr,
            cstr.ptr);
    }
    }
}

void print_filledcurve_color(FILE * fp)
{
    char[8] cstr;

    print_rgb_hash(cstr.ptr, get_shadecolor());
    fprintf(fp, "lc rgb \"%s\" ".ptr, cstr.ptr);
}

/* Given the info in @bars, calculate how many of its start-stop pairs
   fall (at least partially) within the current x-axis range.
*/

int n_bars_shown(double xmin, double xmax, PlotBars * bars)
{
    double ** dx = bars.dx;
    int i, n = 0;

    for (i=0; i<bars.n; i++) {
    if (dx[i][1] < xmin) {
        continue;
    } else if (dx[i][0] > xmax) {
        break;
    } else if (dx[i][0] >= xmin || dx[i][1] <= xmax) {
        n++;
    }
    }

    bars.t1 = xmin;
    bars.t2 = xmax;

    return n;
}

/* Given the limits of the data area of a plot, @xmin et al, and the
   name of a plain text file containing "bars" information in the form
   of a set of start-stop pairs, try attaching this information to
   @spec. We fail if there's something wrong with the info in the file,
   or if none of the bars defined therein would be visible in the
   current range, @xmin to @xmax.
*/

int plotspec_add_bars_info(GptSpec * spec, double xmin, double xmax, double ymin, double ymax, const(char) * fname)
{
    PlotBars * bars = null;
    int err = 0;

    if (spec.bars != null) {
    /* clear out any stale info */
    plotbars_free(bars);
    spec.bars = null;
    spec.nbars = 0;
    }

    bars = parse_bars_file(fname, xmin, xmax, &err);

    if (!err) {
    int n = n_bars_shown(xmin, xmax, bars);

    if (n > 0) {
        spec.bars = bars;
        bars.ymin = ymin;
        bars.ymax = ymax;
        spec.nbars = n;
    } else {
        gretl_errmsg_set(_("No bars fall within the plot range").ptr);
        plotbars_free(bars);
        err = 1;
    }
    }

    return err;
}

/* The following is a public interface because it's also used when
   reconstituting @spec from a gnuplot command file.
*/

int plotspec_allocate_bars(GptSpec * spec)
{
    int err = 0;

    spec.bars = plotbars_new(spec.nbars);
    if (spec.bars == null) {
    err = GretlError.alloc;
    }

    return err;
}

/* for use when reconstituting @spec from a gnuplot
   command file: set the start/stop values for a
   particular bar @i.
*/

int plotspec_set_bar_info(GptSpec * spec, int i, double t1, double t2)
{
    if (i < spec.bars.n) {
    spec.bars.dx[i][0] = t1;
    spec.bars.dx[i][1] = t2;
    return 0;
    } else {
    return GretlError.data;
    }
}

void plotspec_set_bars_limits(GptSpec * spec, double t1, double t2, double ymin, double ymax)
{
    if (spec.bars != null) {
    spec.bars.t1 = t1;
    spec.bars.t2 = t2;
    spec.bars.ymin = ymin;
    spec.bars.ymax = ymax;
    }
}

void plotspec_remove_bars(GptSpec * spec)
{
    plotbars_free(spec.bars);
    spec.bars = null;
    spec.nbars = 0;
}

