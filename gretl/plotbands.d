/* Port of lib/src/plotbands.c (with the legacy band syntax of pb_legacy.c, which it #includes when
   SUPPORT_LEGACY is set): support for multiple "bands" (confidence bands, ranges, bars) in a gnuplot plot. The
   private header plot_priv.h is plot_privh.d. */

module gretl.plotbands;


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
import gretl.plot_privh : BPMode, GPNA, GnuplotInfo;
import gretl.dataset : dataset_add_series, dataset_drop_last_variables, matrix_dataset_expand_Z;
import gretl.genmain : current_series_index;
import gretl.glib : guint32;
import gretl.gmem : g_free;
import gretl.graphingh : GptFlags, PlotType, gretlRGB;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf, gretl_warnmsg_set, invalid_field_error;
import gretl.gretl_list : gretl_list_append_term, gretl_list_delete_at_pos, in_gretl_list;
import gretl.gretl_paths : gretl_plotfile, gretl_remove;
import gretl.gretl_prn : Prn, gretl_print_destroy, gretl_print_detach_stream, gretl_print_new_with_stream;
import gretl.gretl_utils : gretl_isdummy;
import gretl.gstrfuncs : g_strfreev, g_strndup, g_strsplit;
import gretl.options : get_optval_string;
import gretl.strutils : dot_atof, integer_string, numeric_string, strings_array_free;
import gretl.graphing : check_for_yscale, clear_gpinfo, finalize_plot_input_file, get_default_png_scale, get_shadecolor, graph_list_adjust_sample, is_y2_var, make_time_tics, numeric_color_from_string, open_plot_input_file, plotname, print_gnuplot_literal_lines, print_rgb_hash, print_user_y_data, print_x_range, set_plot_withstr;
import gretl.gretl_array : get_array_by_name, gretl_array_get_data, gretl_array_get_length, gretl_array_get_type;
import gretl.gretl_bundle : get_bundle_by_name, gretl_bundle_get_data, gretl_bundle_get_keys_raw, gretl_bundle_get_member_type, gretl_bundle_get_scalar, gretl_bundle_get_string, gretl_bundle_get_uint32, gretl_bundle_has_key;
import gretl.uservar : gretl_is_scalar, gretl_scalar_get_value;
import gretl.gretl_typemap : gretl_is_scalar_type, gretl_type_get_name;
import gretl.usermat : get_matrix_by_name;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}


enum BandStyle {
    line,
    fill,
    dash,
    bars,
    step
}

static const(char) *[5] style_specs = ["line".ptr, "fill".ptr, "dash".ptr, "bars".ptr, "step".ptr];

struct BandInfo {
    int center; /* series ID of center-of-band data */
    int width; /* series ID of width-of-band data */
    double factor; /* factor by which to multiply width */
    int bdummy; /* flag for drawing rectangles */
    BandStyle style; /* see enumeration above */
    char[10] rgb; /* specific color, if wanted */
    const(char) * title;
}

BandInfo * band_info_new()
{
    BandInfo * bi = cast(BandInfo *) malloc(BandInfo.sizeof);

    if (bi != null) {
        bi.center = -1;
        bi.width = -1;
        bi.factor = 1.0;
        bi.bdummy = 0;
        bi.style = BandStyle.line;
        bi.rgb[0] = '\0';
    bi.title = null;
    }

    return bi;
}

void free_bbi(BandInfo ** bbi, int n)
{
    if (bbi != null) {
    int i;

    for (i=0; i<n; i++) {
        free(bbi[i]);
    }
    free(bbi);
    }
}

BandStyle style_from_string(const(char) * s, int * err)
{
    int i, n = cast(int) style_specs.length;

    for (i=0; i<n; i++) {
        if (!strcmp(s, style_specs[i])) {
            return cast(BandStyle) i;
        }
    }

    *err = invalid_field_error(s);
    return cast(BandStyle) 0;
}

/* handle the band-as-matrix case */

int process_band_matrix(const(GretlMatrixRaw) * m, BandInfo * bi, GnuplotInfo * gi, DatasetRaw * dset)
{
    int err = 0;

    if (m == null || m.cols != 2 || m.rows != dset.n) {
    /* matrix missing or non-conformable */
    err = GretlError.nonconf;
    } else {
        /* enlarge the dset->Z array and stick bi->center and
           bi->width onto the end of it
        */
    err = matrix_dataset_expand_Z(dset, 2);
        if (!err) {
            bi.center = dset.v - 2;
            bi.width = dset.v - 1;
            dset.Z[bi.center] = cast(double *) (m.val);
            dset.Z[bi.width] = cast(double *) (m.val + m.rows);
        }
    }

    if (!err) {
        gretl_list_append_term(&gi.list, bi.center);
    bi.center = gi.list[0] + 1;
        gretl_list_append_term(&gi.list, bi.width);
    bi.width = gi.list[0] + 1;
    }

    return err;
}

/* check the validity of a series that's supposed to be a dummy
   variable (to produce something resembling recession bars)
*/

int check_assign_bdummy(BandInfo * bi, int vnum, const(DatasetRaw) * dset)
{
    int err = 0;

    if (gretl_isdummy(dset.t1, dset.t2, dset.Z[vnum])) {
    bi.bdummy = vnum;
    } else {
    err = GretlError.invarg;
    fprintf(stderr, "%s: not a dummy variable\n".ptr,
        dset.varname[vnum]);
    }

    return err;
}

void do_center_or_width(BandInfo * bi, GnuplotInfo * gi, int i, int v)
{
    int pos = in_gretl_list(gi.list, v);

    if (pos == 0) {
    /* not already in list */
    gretl_list_append_term(&gi.list, v);
    pos = gi.list[0];
    }

    if (i == 0) {
    bi.center = pos + 1;
    } else {
    bi.width = pos + 1;
    }
}

int has_bogus_key(GretlBundle * b)
{
    static const(char) *[8] known_keys = ["center".ptr, "width".ptr, "dummy".ptr, "style".ptr, "color".ptr, "factor".ptr, "bandmat".ptr, "title".ptr];
    char ** S = null;
    int i, j, ns = 0;
    int ok = 0;

    S = gretl_bundle_get_keys_raw(b, &ns);

    for (i=0; i<ns; i++) {
    ok = 0;
    for (j=0; known_keys[j] != null; j++) {
        if (!strcmp(S[i], known_keys[j])) {
        ok = 1;
        break;
        }
    }
    if (!ok) {
        gretl_errmsg_sprintf("Unknown key '%s'".ptr, S[i]);
        break;
    }
    }

    strings_array_free(S, ns);

    return !ok;
}

int get_input_color(GretlBundle * b, char * rgb)
{
    GretlType t = gretl_bundle_get_member_type(b, "color", null);
    guint32 u;
    int err = 0;

    if (t == GretlType.string_) {
    const(char) * s = gretl_bundle_get_string(b, "color", &err);

    if (!err) {
        if (!strcmp(s, "default")) {
        /* ensure @rgb is blank */
        rgb[0] = '\0';
        return 0;
        } else {
        u = numeric_color_from_string(s, &err);
        }
    }
    } else if (gretl_is_scalar_type(t)) {
    u = gretl_bundle_get_uint32(b, "color", &err);
    } else {
    err = GretlError.types;
    }

    if (!err) {
    sprintf(rgb, "#%x".ptr, u);
    }

    return err;
}

const(GretlMatrixRaw) * retrieve_bandmat(GretlBundle * b, int * err)
{
    GretlMatrixRaw * m = null;
    GretlType t;
    void * ptr;

    ptr = gretl_bundle_get_data(b, "bandmat", &t, err);

    if (!*err) {
        if (t == GretlType.matrix) {
            m = cast(GretlMatrixRaw *) ptr;
        } else if (t == GretlType.string_) {
            m = get_matrix_by_name(cast(const(char) *) ptr);
        } else {
            *err = GretlError.types;
        }
    }

    return m;
}

struct BandData {
    GretlType t;
    char * s;
    double c;
};

int band_type_error(BandData * bdata, const(char) * key)
{
    gretl_errmsg_sprintf("%s: invalid type %s in band bundle".ptr,
             key, gretl_type_get_name(bdata.t));
    return GretlError.types;
}

int band_data_error(const(char) * key)
{
    gretl_errmsg_sprintf("%s: failed to get a series".ptr, key);
    return GretlError.data;
}

int make_const_series(double c, DatasetRaw * dset, int * v_next)
{
    if (c == 1.0) {
    return 0; /* use "const" */
    } else {
    int t, v = *v_next;

    for (t=dset.t1; t<=dset.t2; t++) {
        dset.Z[v][t] = c;
    }
    *v_next += 1;
    return v;
    }
}

BandInfo * band_info_from_bundle(int matrix_mode, GretlBundle * b, GnuplotInfo * gi, DatasetRaw * dset, int * err)
{
    BandInfo *bi = band_info_new();
    BandData[3] bdata;
    memset(bdata.ptr, 0, bdata.sizeof);
    const(char) *[3] keys = ["center".ptr, "width".ptr, "dummy".ptr];
    int v_orig = dset.v;

    if (has_bogus_key(b)) {
    *err = GretlError.data;
    return null;
    }

    /* Allow for the "matrix_mode" case where the band (center, width)
       is represented by a two-column matrix. This arises when
       "gnuplot" or "plot" is fed data in matrix form.
    */
    if (matrix_mode) {
        const(GretlMatrixRaw) * m = retrieve_bandmat(b, err);

    if (!*err) {
        *err = process_band_matrix(m, bi, gi, dset);
    }
    if (*err) {
        return null;
    }
    } else {
    int[3] got = 0;
    int v_next = v_orig;
    int v_add = 0;
    int i, v = 0;
    int imin = 0;
    int imax = 3;
    void * data;

    /* first pass */
    for (i=imin; i<imax && !*err; i++) {
        data = gretl_bundle_get_data(b, keys[i], &(bdata[i].t), null);
        if (data != null) {
        if (bdata[i].t == GretlType.string_) {
            bdata[i].s = cast(char *) data;
        } else if (bdata[i].t == GretlType.double_) {
            if (i < 2) {
            bdata[i].c = *cast(double *) data;
            v_add += (bdata[i].c != 1.0);
            } else {
            *err = band_type_error(&bdata[i], keys[i]);
            }
        } else {
            *err = band_type_error(&bdata[i], keys[i]);
        }
        got[i] = 1;
        }
    }
    /* we should have either center and width, or dummy alone */
    if (got[0] && got[1] && !got[2]) {
        imax = 2;
    } else if (got[2] && !got[0] && !got[1]) {
        imin = 2;
    } else {
        gretl_errmsg_set(_("Invalid band specification").ptr);
        *err = GretlError.invarg;
    }
    if (!*err && v_add > 0) {
        *err = dataset_add_series(dset, v_add);
    }
    /* second pass */
    for (i=imin; i<imax && !*err; i++) {
        v = -1;
        if (bdata[i].t == GretlType.string_) {
        v = current_series_index(dset, bdata[i].s);
        } else if (bdata[i].t == GretlType.double_) {
        v = make_const_series(bdata[i].c, dset, &v_next);
        }
        if (v < 0) {
        *err = band_data_error(keys[i]);
        }
        if (!*err) {
        if (i < 2) {
            do_center_or_width(bi, gi, i, v);
        } else {
            /* should be a dummy */
            *err = check_assign_bdummy(bi, v, dset);
        }
        }
        }
    }

    if (*err) {
    goto bailout;
    }

    /* optional bundle members */
    if (!*err && gretl_bundle_has_key(b, "color")) {
    *err = get_input_color(b, bi.rgb.ptr);
    }
    if (!*err && gretl_bundle_has_key(b, "factor")) {
        bi.factor = gretl_bundle_get_scalar(b, "factor", err);
        if (!*err && (bi.factor <= 0 || na(bi.factor))) {
            *err = GretlError.invarg;
        }
    }
    if (!*err && gretl_bundle_has_key(b, "style")) {
    const(char) * s = gretl_bundle_get_string(b, "style", err);
    if (s != null) {
        bi.style = style_from_string(s, err);
    }
    }
    if (!*err && gretl_bundle_has_key(b, "title")) {
        bi.title = gretl_bundle_get_string(b, "title", err);
    }

 bailout:

    if (*err) {
    if (dset.v > v_orig) {
        dataset_drop_last_variables(dset, dset.v - v_orig);
    }
        free(bi);
        bi = null;
    }

    return bi;
}

int band_matrices_differ(const(GretlMatrixRaw) * a, const(GretlMatrixRaw) * b)
{
    int n = a.rows * a.cols;

    if (n != b.rows * b.cols) {
        return 1;
    } else {
        int i;

        for (i=0; i<n; i++) {
        if (na(a.val[i]) && na(b.val[i])) {
        continue;
        } else if (a.val[i] != b.val[i]) {
                return 1;
            }
        }
    }

    return 0;
}

/* Here we're trying to determine a plausible z-order for two bands
   that both employ the "fill" style, to avoid having a narrower band
   obscured by a wider one. The idea is that if a multiplicative
   "factor" is present in both cases, and the two bands have the same
   center and width, they can be ordered by the factor. Otherwise we
   treat them as unordered.
*/

int swap_bands_order(GretlArray * a, int * err)
{
    GretlBundle * b0 = cast(GretlBundle *) gretl_array_get_data(a, 0);
    GretlBundle * b1 = cast(GretlBundle *) gretl_array_get_data(a, 1);
    const(char) * s0;
    const(char) * s1;
    double f0, f1;

    f0 = gretl_bundle_get_scalar(b0, "factor", null);
    f1 = gretl_bundle_get_scalar(b1, "factor", null);
    if (na(f0) || na(f1)) {
        /* we won't try to order the two */
        return 0;
    }

    s0 = gretl_bundle_get_string(b0, "style", null);
    s1 = gretl_bundle_get_string(b1, "style", null);
    if (s0 == null || s1 == null) {
        return 0;
    } else if (strcmp(s0, "fill") || strcmp(s1, "fill")) {
        /* condition ordering on "fill" for both bands? */
        return 0;
    }

    if (gretl_bundle_has_key(b0, "bandmat") &&
        gretl_bundle_has_key(b1, "bandmat")) {
        const(GretlMatrixRaw) * m0 = retrieve_bandmat(b0, err);
        const(GretlMatrixRaw) * m1 = retrieve_bandmat(b1, err);

        if (m0 == null || m1 == null || band_matrices_differ(m0, m1)) {
            return 0;
        }
    } else if (gretl_bundle_has_key(b0, "center") &&
               gretl_bundle_has_key(b1, "center") &&
               gretl_bundle_has_key(b0, "width") &&
               gretl_bundle_has_key(b1, "width")) {
        s0 = gretl_bundle_get_string(b0, "center", null);
        s1 = gretl_bundle_get_string(b1, "center", null);
        if (s0 == null || s1 == null || strcmp(s0, s1)) {
            /* we don't have the same centers */
            return 0;
        }
        s0 = gretl_bundle_get_string(b0, "width", null);
        s1 = gretl_bundle_get_string(b1, "width", null);
        if (s0 == null || s1 == null || strcmp(s0, s1)) {
            /* we don't have the same widths */
            return 0;
        }
    }

    return f1 > f0;
}

/* Here we're supposing we got --bands=@aname. We try to retrieve an
   array of bundles, and if that succeeds we build and return an array
   of BandInfo structs.
*/

BandInfo ** get_band_info_array(int matrix_mode, int * n_bands, GnuplotInfo * gi, DatasetRaw * dset, int * err)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_A);
    GretlArray * a = get_array_by_name(s);
    BandInfo **pbi = null;
    int nb = 0;
    int swap = 0;

    if ((nb = gretl_array_get_length(a)) < 1) {
    fprintf(stderr, "get_band_info_array: array has no content\n".ptr);
        *err = GretlError.data;
        return null;
    } else if (gretl_array_get_type(a) != GretlType.bundles) {
    fprintf(stderr, "get_band_info_array: array has wrong type\n".ptr);
        *err = GretlError.types;
        return null;
    }

    if (nb == 2) {
        swap = swap_bands_order(a, err);
    }

    pbi = cast(typeof(pbi)) calloc(nb, (*pbi).sizeof);
    if (pbi == null) {
        *err = GretlError.alloc;
    } else {
    GretlBundle * b;
        int i, j;

        for (i=0; i<nb && !*err; i++) {
            b = cast(GretlBundle *) gretl_array_get_data(a, swap ? 1 - i : i);
            pbi[i] = band_info_from_bundle(matrix_mode, b, gi, dset, err);
        }
        if (*err) {
            for (j=0; j<=i; j++) {
                free(pbi[j]);
            }
            free(pbi);
            pbi = null;
        } else {
            *n_bands = nb;
        }
    }

    return pbi;
}

const(char) * get_lw_string()
{
    if (get_default_png_scale() > 1.0) {
    return " lw 2";
    } else {
    return "";
    }
}

void gp_newline(int contd, FILE * fp)
{
    if (contd) {
    fputs(", \\\n", fp);
    } else {
    fputc('\n', fp);
    }
}

void band_title(BandInfo * bi, FILE * fp)
{
    if (bi.title != null) {
    fprintf(fp, "title '%s' ".ptr, bi.title);
    } else {
    fputs("notitle ", fp);
    }
}

void print_pm_filledcurve(BandInfo * bi, int contd, FILE * fp)
{
    double f = bi.factor;
    int c = bi.center;
    int w = bi.width;
    char[10] cstr;

    if (bi.rgb[0] != '\0') {
        strcpy(cstr.ptr, bi.rgb.ptr);
    } else {
        print_rgb_hash(cstr.ptr, get_shadecolor());
    }

    if (bi.factor == 1.0) {
    fprintf(fp, "$data using 1:($%d-$%d):($%d+$%d) ".ptr, c, w, c, w);
    } else {
    fprintf(fp, "$data using 1:($%d-%g*$%d):($%d+%g*$%d) ".ptr,
        c, f, w, c, f, w);
    }
    band_title(bi, fp);
    fprintf(fp, "lc rgb \"%s\" w filledcurve".ptr, cstr.ptr);
    gp_newline(contd, fp);
}

/* for plain lines, dashed lines, or steps */

void print_pm_lines(BandInfo * bi, int n_yvars, const(char) * lw, int contd, FILE * fp)
{
    const(char) * wstr = (bi.style == BandStyle.step)? "steps".ptr : "lines".ptr;
    char[24] lspec;
    char[8] dspec = 0;
    double f = bi.factor;
    int c = bi.center;
    int w = bi.width;

    if (bi.rgb[0] != '\0') {
    sprintf(lspec.ptr, "lc rgb \"%s\"".ptr, bi.rgb.ptr);
    } else {
    /* default to gray */
    strcpy(lspec.ptr, "lc rgb \"#bbbbbb\"");
    }
    if (bi.style == BandStyle.dash) {
    strcpy(dspec.ptr, " dt 2");
    }
    /* lower line */
    fprintf(fp, "$data using 1:($%d-%g*$%d) ".ptr, c, f, w);
    band_title(bi, fp);
    fprintf(fp, "w %s %s%s%s, \\\n".ptr, wstr, lspec.ptr, dspec.ptr, lw);
    /* upper line */
    fprintf(fp, "$data using 1:($%d+%g*$%d) notitle w %s %s%s%s".ptr,
        c, f, w, wstr, lspec.ptr, dspec.ptr, lw);
    gp_newline(contd, fp);
}

/* for errorbars only */

void print_pm_bars(BandInfo * bi, int n_yvars, int contd, FILE * fp)
{
    char[24] lspec;

    /* note: "pt 7" is a solid circle */

    if (bi.rgb[0] != '\0') {
    sprintf(lspec.ptr, "lc rgb \"%s\" pt 7".ptr, bi.rgb.ptr);
    } else {
    sprintf(lspec.ptr, "lt %d pt 7".ptr, n_yvars + 1);
    }
    fprintf(fp, "$data using 1:%d:(%g*$%d) ".ptr, bi.center, bi.factor, bi.width);
    band_title(bi, fp);
    fprintf(fp, "w errorbars %s".ptr, lspec.ptr);
    gp_newline(contd, fp);
}

/* Print the data in the form of a gnuplot "heredoc": this should be
   relatively efficient, avoiding duplication of the x-axis variable
   and the band center.
*/

void print_data_block(GnuplotInfo * gi, const(double) * x, const(DatasetRaw) * dset, int * show_zero, FILE * fp)
{
    double yt;
    int gt_zero = 0;
    int lt_zero = 0;
    int n = gi.list[0];
    int t, i;

    fputs("# start inline data\n", fp);
    fputs("$data << EOD\n", fp);
    for (t=gi.t1; t<=gi.t2; t++) {
    fprintf(fp, "%g".ptr, x[t]);
    for (i=0; i<n; i++) {
        yt = dset.Z[gi.list[i+1]][t];
        if (na(yt)) {
        fprintf(fp, " %s".ptr, GPNA.ptr);
        } else {
        fprintf(fp, " %.10g".ptr, yt);
        }
        if (yt > 0) {
        gt_zero++;
        } else if (yt < 0) {
        lt_zero++;
        }
    }
    fputc('\n', fp);
    }
    fputs("EOD\n", fp);
    fputs("# end inline data\n", fp);

    if (gt_zero && lt_zero) {
    *show_zero = 1;
    }
}

BandInfo ** get_single_band_info(int matrix_mode, GnuplotInfo * gi, DatasetRaw * dset, Gretlopt opt, int * err)
{
    const(char) * spec = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_N);
    GretlBundle * b = null;
    BandInfo **bi = null;

    if (spec == null || *spec == '\0') {
        *err = GretlError.invarg;
    return null;
    }

    b = get_bundle_by_name(spec);
    if (b != null) {
        bi = cast(typeof(bi)) malloc((*bi).sizeof);
        bi[0] = band_info_from_bundle(matrix_mode, b, gi, dset, err);
        if (*err) {
            free(bi);
            bi = null;
        }
        return bi; /* done */
    }

    bi = legacy_get_band_info(spec, matrix_mode, gi, dset, opt, err);
    gretl_warnmsg_set(_("Deprecated --band syntax: please use the current syntax\nas described in the help for the 'gnuplot' command").ptr);

    return bi;
}

/* Write full-height bars as gnuplot rectangle objects, using
   the dummy variable bi->bdummy for on/off information.
*/

int write_rectangles(BandInfo * bi, GnuplotInfo * gi, const(double) * x, DatasetRaw * dset, FILE * fp)
{
    double * d = dset.Z[bi.bdummy];
    char[16] stobs;
    char[16] endobs;
    int bar_on = 0, obj = 1;
    int t, err = 0;

    if (bi.rgb[0] == '\0') {
        strcpy(bi.rgb.ptr, "#dddddd");
    }

    *stobs.ptr = *endobs.ptr = '\0';

    for (t=gi.t1; t<=gi.t2; t++) {
        if (na(d[t])) {
            err = GretlError.missdata;
            break;
        }
        if (bar_on && d[t] == 0) {
            /* finalize a bar */
            sprintf(endobs.ptr, "%g".ptr, x[t]);
            fprintf(fp, "set object %d rectangle from %s, graph 0 to %s, graph 1 back fillstyle solid 0.5 noborder fc rgb \"%s\"\n".ptr,
                    obj++, stobs.ptr, endobs.ptr, bi.rgb.ptr);
            bar_on = 0;
        } else if (!bar_on && d[t] != 0) {
            /* start a bar */
            sprintf(stobs.ptr, "%g".ptr, x[t]);
            bar_on = 1;
        }
    }

    if (bar_on) {
        /* terminate an unfinished bar */
        sprintf(endobs.ptr, "%g".ptr, x[gi.t2]);
        fprintf(fp, "set object rectangle from %s, graph 0 to %s, graph 1 back fillstyle solid 0.5 noborder fc rgb \"%s\"\n".ptr,
                stobs.ptr, endobs.ptr, bi.rgb.ptr);
    }

    return err;
}

/* Below: the bars in question don't have to indicate recessions,
   but they're the same sort of thing: full-height bars indicating
   the presence of absence of some binary time-varying attribute.
   This function handles the case where the only "band" instance
   is a case of such bars.
*/

void recession_bars_plot(BandInfo * bi, GnuplotInfo * gi, int n_yvars, const(double) * x, DatasetRaw * dset, Gretlopt opt, FILE * fp)
{
    const(double) * y;
    char[16] wspec = 0;
    int i;

    /* write out the rectangles as objects */
    write_rectangles(bi, gi, x, dset, fp);

    if (!(opt & Gretlopt.OPT_Y)) {
    check_for_yscale(gi, cast(const(double) **) dset.Z);
    if (gi.flags & GptFlags.y2axis) {
        fputs("set ytics nomirror\n", fp);
        fputs("set y2tics\n", fp);
    }
    }

    fputs("plot \\\n", fp);

    /* plot the actual data */
    for (i=1; i<=n_yvars; i++) {
    const(char) * iname = plotname(dset, gi.list[i], 1);

    set_plot_withstr(gi, i, wspec.ptr);
    if (gi.flags & GptFlags.y2axis) {
        fprintf(fp, "'-' using 1:2 axes %s title \"%s (%s)\" %s lt %d".ptr,
            is_y2_var(gi, i) ? "x1y2".ptr : "x1y1".ptr, iname,
            is_y2_var(gi, i) ? _("right").ptr : _("left").ptr,
            wspec.ptr, i);
    } else {
        fprintf(fp, "'-' using 1:2 title \"%s\" %s lt %d".ptr, iname, wspec.ptr, i);
    }
    gp_newline(i < n_yvars, fp);
    }

    /* and write the data block */
    for (i=0; i<n_yvars; i++) {
    y = dset.Z[gi.list[i+1]];
    print_user_y_data(x, y, gi.t1, gi.t2, fp);
    }
}

/* The last entry in gi->list pertains to the x-axis variable:
   either the index of a "genuine" x-var or a dummy placeholder
   of 0. Either way we want to drop this term from the list,
   and return a pointer to the actual x-axis data.
*/

const(double) * gi_get_xdata(GnuplotInfo * gi, char * xname, DatasetRaw * dset)
{
    int xpos = gi.list[0];
    const(double) * x;

    if (gi.x != null) {
    /* the "plotx" (dummy) case */
    x = gi.x;
    *xname = '\0';
        if (gi.flags & GptFlags.ts) {
            gi.flags |= GptFlags.letterbox;
        }
    } else {
    int xno = gi.list[xpos];

    x = dset.Z[xno];
    strcpy(xname, plotname(dset, xno, 1));
    }

    gretl_list_delete_at_pos(gi.list, xpos);

    return x;
}

/* If we have a single data series to plot and no band is
   titled, omit the key and put the name of the series on
   the y axis.
*/

int maybe_suppress_key(BandInfo ** bb, int n_bands, char * yname, const(char) * src, FILE * fp)
{
    int i, have_titles = 0;

    for (i=0; i<n_bands; i++) {
    if (bb[i].title != null) {
        have_titles = 1;
        break;
    }
    }
    if (!have_titles) {
    strcpy(yname, src);
    fprintf(fp, "set ylabel \"%s\"\n".ptr, yname);
    return 1;
    } else {
    return 0;
    }
}

int plot_with_band(BPMode mode, GnuplotInfo * gi, const(char) * literal, DatasetRaw * dset, Gretlopt opt)
{
    BandInfo **bbi = null;
    BandInfo *bi = null;
    FILE * fp = null;
    const(double) * x = null;
    const(char) * lwstr;
    char[MAXDISP] yname;
    char[MAXDISP] xname;
    char[16] wspec = 0;
    int show_zero = 0;
    int no_key = 0;
    int i, j, n_yvars = 0;
    int matrix_mode;
    int n_bands = 1;
    int pushed = 0;
    int v_orig;
    int err = 0;

    matrix_mode = (opt & Gretlopt.OPT_X) || mode == BPMode.blockmat;
    v_orig = dset.v;

    /* subtract 1 for x */
    n_yvars = gi.list[0] - 1;

    x = gi_get_xdata(gi, xname.ptr, dset);

    if (opt & Gretlopt.OPT_A) {
    /* --bands=@array */
    bbi = get_band_info_array(matrix_mode, &n_bands, gi, dset, &err);
    } else {
    /* --band=<whatever> */
    bbi = get_single_band_info(matrix_mode, gi, dset, opt, &err);
    }
    if (err) {
    return err;
    }

    err = graph_list_adjust_sample(gi, dset, 1);
    if (err) {
    goto finish;
    }

    fp = open_plot_input_file(PlotType.band, gi.flags, &err);
    if (err) {
    goto finish;
    }

    if (gi.flags & GptFlags.ts) {
        Prn * prn = gretl_print_new_with_stream(fp);

        make_time_tics(gi, dset, 0, null, prn);
        gretl_print_detach_stream(prn);
        gretl_print_destroy(prn);
    }

    if (n_yvars == 1) {
    const(char) * s = plotname(dset, gi.list[1], 1);

    no_key = maybe_suppress_key(bbi, n_bands, yname.ptr, s, fp);
    }

    if (no_key) {
    fputs("set nokey\n", fp);
    } else {
    fputs("set key left top\n", fp);
    }

    if (*xname.ptr != '\0') {
        fprintf(fp, "set xlabel \"%s\"\n".ptr, xname.ptr);
    }

    /* for adjusting line width if wanted */
    lwstr = get_lw_string();

    gretl_push_c_numeric_locale();
    pushed = 1;

    if (gi.x != null) {
        print_x_range(gi, fp);
    }

    if (literal != null) {
    print_gnuplot_literal_lines(literal, fp);
    }

    if (show_zero && bi.style != BandStyle.fill) {
    /* in the "fill" case xzeroaxis won't be visible */
        fputs("set xzeroaxis\n", fp);
    }

    if (bbi[0].bdummy) {
    /* special case */
    if (n_bands == 1) {
        recession_bars_plot(bbi[0], gi, n_yvars,
                x, dset, opt, fp);
        goto finish;
    } else {
        err = write_rectangles(bbi[0], gi, x, dset, fp);
    }
    }

    print_data_block(gi, x, dset, &show_zero, fp);

    fputs("plot \\\n", fp);

    /* bands first */
    for (j=0; j<n_bands; j++) {
    bi = bbi[j];
    if (bi.bdummy) {
        continue;
    } else if (bi.style == BandStyle.fill) {
        print_pm_filledcurve(bi, 1, fp);
        if (show_zero) {
        /* custom xzeroaxis on top of the fill */
        fputs("0 notitle w lines lc rgb \"#888888\" dt 2, \\\n", fp);
        }
    } else if (bi.style == BandStyle.bars) {
        print_pm_bars(bi, n_yvars, 1, fp);
    } else {
        print_pm_lines(bi, n_yvars, lwstr, 1, fp);
    }
    }
    /* then the non-band data */
    for (i=1; i<=n_yvars; i++) {
    const(char) * iname = plotname(dset, gi.list[i], 1);

    set_plot_withstr(gi, i, wspec.ptr);
    fprintf(fp, "$data using 1:%d title '%s' %s lt %d%s".ptr, i+1,
        iname, wspec.ptr, i, lwstr);
    gp_newline(i < n_yvars, fp);
    }

 finish:

    if (pushed) {
    gretl_pop_c_numeric_locale();
    }

    if (!matrix_mode && dset.v > v_orig) {
    dataset_drop_last_variables(dset, dset.v - v_orig);
    }

    if (err) {
    if (fp != null) {
        fclose(fp);
        gretl_remove(gretl_plotfile());
    }
    } else {
    err = finalize_plot_input_file(fp);
    }

    clear_gpinfo(gi);
    free_bbi(bbi, n_bands);

    if (matrix_mode) {
        /* hide the two extra dataset columns
           representing the band */
        dset.v -= 2;
    }

    return err;
}

/* Here: functions supporting the legacy plot band syntax that was in
   force prior to the gretl 2023c release.
*/

/* Handle the case where we get to the band-plot code from a command
   in which the data to be plotted (and hence also the band
   specification) are given in matrix form. By this point the
   plot-data have been converted to (temporary) DATASET form; here we
   retrieve the band-spec matrix, check it for conformability, and
   stick the two extra columns onto the dataset (borrowing pointers
   into the matrix content).
*/

int parse_band_matrix_option(BandInfo * bi, const(char) * spec, GnuplotInfo * gi, DatasetRaw * dset)
{
    const(char) * mname = null;
    char ** S;
    int i = 0;
    int err = 0;

    S = g_strsplit(spec, ",", -1);

    while (S != null && S[i] != null && !err) {
        if (i == 0) {
            mname = S[i];
        } else if (i == 1) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi.factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi.factor = gretl_scalar_get_value(S[i], &err);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    if (!err && (bi.factor <= 0 || na(bi.factor))) {
        err = GretlError.invarg;
    }
    if (!err) {
        GretlMatrixRaw * m = get_matrix_by_name(mname);

        if (m == null) {
            err = GretlError.invarg;
        } else {
            err = process_band_matrix(m, bi, gi, dset);
        }
    }

    g_strfreev(S);

    return err;
}

int handle_recession_bars(BandInfo * bi, const(char) * s, const(DatasetRaw) * dset)
{
    int v = current_series_index(dset, s);
    int err = 0;

    if (v >= 0 && v < dset.v) {
        err = check_assign_bdummy(bi, v, dset);
    } else {
        err = GretlError.invarg;
    }

    return err;
}

/* Handle the band plus-minus option for all cases apart from the one
   handled just above. Here we require two comma-separated series
   identifiers for center and width.
*/

int parse_band_pm_option(BandInfo * bi, const(char) * spec, GnuplotInfo * gi, const(DatasetRaw) * dset, Gretlopt opt)
{
    char ** S;
    int v, i = 0;
    int err = 0;

    if (strchr(spec, ',') == null) {
        return handle_recession_bars(bi, spec, dset);
    }

    S = g_strsplit(spec, ",", -1);

    while (S != null && S[i] != null && !err) {
        if (i < 2) {
            /* specs for the "center" and "width" series: required */
            if (integer_string(S[i])) {
                /* var ID number? */
                v = atoi(S[i]);
            } else {
                /* varname? */
                v = current_series_index(dset, S[i]);
            }
            if (v >= 0 && v < dset.v) {
                do_center_or_width(bi, gi, i, v);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else if (i == 2) {
            /* spec for width multiplier: optional */
            if (numeric_string(S[i])) {
                bi.factor = dot_atof(S[i]);
            } else if (gretl_is_scalar(S[i])) {
                bi.factor = gretl_scalar_get_value(S[i], &err);
            } else {
                err = invalid_field_error(S[i]);
            }
        } else {
            /* we got too many comma-separated terms */
            err = invalid_field_error(S[i]);
        }
        i++;
    }

    g_strfreev(S);

    if (!err) {
        if (bi.center < 0 || bi.width < 0 ||
            bi.factor < 0 || na(bi.factor)) {
            err = GretlError.invarg;
        }
    }

    return err;
}

/* We're looking here for any one of three patterns:

   <style>
   <style>,<color>
   <color>

   where <style> should be "fill", "dash" or "line" (the default) and
   <color> should be a hex string such as "#00ff00" or "0x00ff00".
*/

int parse_band_style_option(BandInfo * bi)
{
    const(char) * s = get_optval_string(GretlCmdIndex.gnuplot, Gretlopt.OPT_J);
    int err = 0;

    if (s != null) {
        const(char) * p = strchr(s, ',');
        int do_color = 1;
        guint32 u = 0;

        if (bi.bdummy && *s != ',') {
            /* must be just a color */
            u = numeric_color_from_string(s, &err);
        } else if (*s == ',') {
            /* skipping field 1, going straight to color */
            u = numeric_color_from_string(s + 1, &err);
        } else if (p == null) {
            /* just got field 1, style spec */
            bi.style = style_from_string(s, &err);
            do_color = 0;
        } else {
            /* embedded comma: style + color */
            if (strlen(s) >= 8 && s[4] == ',') {
                char * tmp = g_strndup(s, 4);

                bi.style = style_from_string(tmp, &err);
                g_free(tmp);
            } else {
                err = invalid_field_error(s);
            }
            if (!err) {
                u = numeric_color_from_string(s + 5, &err);
            }
            if (!err && do_color) {
                sprintf(bi.rgb.ptr, "#%x".ptr, u);
            }
        }
    }

    return err;
}

BandInfo ** legacy_get_band_info(const(char) * spec, int matrix_mode, GnuplotInfo * gi, DatasetRaw * dset, Gretlopt opt, int * err)
{
    BandInfo **bi = null;

    bi = cast(typeof(bi)) malloc((*bi).sizeof);
    bi[0] = band_info_new();
    if (bi[0] == null) {
        *err = GretlError.alloc;
        return null;
    }

    if (matrix_mode) {
        /* In this case the band should be given in the form of a
           named matrix with two columns holding center and width,
           respectively.
        */
        *err = parse_band_matrix_option(bi[0], spec, gi, dset);
    } else {
        *err = parse_band_pm_option(bi[0], spec, gi, dset, opt);
    }
    if (!*err && (opt & Gretlopt.OPT_J)) {
        *err = parse_band_style_option(bi[0]);
    }

    if (*err) {
        free_bbi(bi, 1);
        bi = null;
    }

    return bi;
}

