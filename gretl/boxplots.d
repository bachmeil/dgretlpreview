/* Port of lib/src/boxplots.c: boxplots (single, grouped and factorized) written as gnuplot scripts: quartiles,
   whiskers and outliers, bootstrap confidence intervals for the median, the numerical summary and the reading
   back of a boxplot file. The header boxplots.h has only prototypes. */

module gretl.boxplots;


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
import gretl.dataset : create_auxiliary_dataset, destroy_dataset, multi_unit_panel_sample, panel_sample_size, series_get_label, series_get_string_table, series_set_label;
import gretl.dataseth : sample_size;
import gretl.describe : gretl_array_quantile, gretl_array_quantiles, gretl_mean;
import gretl.graphing : finalize_plot_input_file, open_plot_input_file, plotname, print_gnuplot_literal_lines;
import gretl.graphingh : GptFlags, PlotType;
import gretl.gretl_commandsh : GretlCmdIndex;
import gretl.gretl_errors : gretl_errmsg_set, gretl_errmsg_sprintf;
import gretl.gretl_list : gretl_consecutive_list_new;
import gretl.gretl_paths : gretl_fopen;
import gretl.gretl_prn : Prn, pprintf, pputc, pputs;
import gretl.gretl_string_table : SeriesTable, series_table_get_string;
import gretl.gretl_utils : accept_as_discrete, gretl_compare_doubles, gretl_isconst, transcribe_array;
import gretl.options : get_optval_double;
import gretl.strutils : gretl_namechar_spn, gretl_quoted_string_strdup, gretl_scan_varname, gretl_strdup, gretl_strdup_printf, gretl_strndup, numeric_string, strings_array_free, strings_array_new;
import gretl.random : gretl_rand_int_max;

/* unported dependencies: per-file declarations from the 2026c headers */
extern (C) {
}

struct Boxplot {
    double mean;
    double median;
    double[2] conf; /* bounds of confidence interval */
    double uq, lq; /* upper and lower quartiles */
    double max, min; /* data max and min */
    double wmax, wmin; /* whisker max and min */
    int n;
    char[VNAMELEN] varname;
    char * boolstr; /* FIXME */
    GretlMatrixRaw * outliers;
}

enum BoxFlags {
    show_mean = 1 << 0,
    intervals = 1 << 1,
    whiskbars = 1 << 2,
    xtics     = 1 << 3,
    factor    = 1 << 4
}

struct PlotGroup {
    int flags;
    int nplots;
    const(int) * list;
    char * title;
    int n_bools;
    int pan_min;
    Boxplot * plots;
    double gmin, gmax;
    double limit;
    double * x;
    GretlMatrixRaw * dvals;
    char ** strvals;
    char[VNAMELEN] xlabel;
    char[VNAMELEN] ylabel;
    const(char) * literal;
}

bool show_mean(const(PlotGroup) * g)
{
    return (g.flags & BoxFlags.show_mean) != 0;
}
bool do_intervals(const(PlotGroup) * g)
{
    return (g.flags & BoxFlags.intervals) != 0;
}
bool whiskerbars(const(PlotGroup) * g)
{
    return (g.flags & BoxFlags.whiskbars) != 0;
}
bool use_xtics(const(PlotGroup) * g)
{
    return (g.flags & BoxFlags.xtics) != 0;
}
bool factorized(const(PlotGroup) * g)
{
    return (g.flags & BoxFlags.factor) != 0;
}

void set_show_mean(PlotGroup * g)
{
    g.flags |= BoxFlags.show_mean;
}
void set_do_intervals(PlotGroup * g)
{
    g.flags |= BoxFlags.intervals;
}
void set_use_xtics(PlotGroup * g)
{
    g.flags |= BoxFlags.xtics;
}
void set_box_factor(PlotGroup * g)
{
    g.flags |= BoxFlags.factor;
}

void unset_show_mean(PlotGroup * g)
{
    g.flags &= ~BoxFlags.show_mean;
}
void unset_do_intervals(PlotGroup * g)
{
    g.flags &= ~BoxFlags.intervals;
}

void quartiles_etc(double * x, int n, Boxplot * plt, double limit)
{
    double[3] p = [0.25, 0.5, 0.75];

    if (n <= 3) {
        /* special treatment for tiny samples */
        int i;

        if (n == 1) {
            plt.wmin = plt.min = x[0];
            plt.wmax = plt.max = x[0];
            plt.mean = plt.median = NADBL;
            plt.lq = plt.uq = NADBL;
        } else if (n == 2) {
            plt.wmin = plt.min = x[0];
            plt.wmax = plt.max = x[1];
            plt.mean = (x[0] + x[1]) / 2.0;
            plt.median = plt.lq = plt.uq = NADBL;
        } else {
            /* n = 3 */
            plt.wmin = plt.min = x[0];
            plt.wmax = plt.max = x[2];
            plt.median = x[1];
            plt.lq = plt.uq = NADBL;
            plt.mean = (x[0] + x[1] + x[2]) / 3.0;
        }

        plt.outliers = gretl_vector_alloc(n);

        if (plt.outliers != null) {
            for (i=0; i<n; i++) {
                gretl_vector_set(plt.outliers, i, x[i]);
            }
        }

        return;
    }

    plt.wmin = plt.min = x[0];
    plt.wmax = plt.max = x[n-1];

    gretl_array_quantiles(x, n, p.ptr, 3);
    plt.lq = p[0];
    plt.median = p[1];
    plt.uq = p[2];

    if (na(plt.lq) || na(plt.uq)) {
        /* ensure that Q1 and Q3 are either both present
           or both missing
        */
        plt.lq = plt.uq = NADBL;
    }

    plt.mean = gretl_mean(0, n-1, x);

    if (limit > 0 && !na(plt.lq)) {
        /* set suitable multiple of IQR */
        double d = limit * (plt.uq - plt.lq);
        /* set limits beyond which points are outliers */
        double xlo = plt.lq - d;
        double xhi = plt.uq + d;
        int i, nlo = 0, nhi = 0;

        /* initialize whiskers to zero length */
        plt.wmin = plt.lq;
        plt.wmax = plt.uq;

        /* work on lower region... */
        for (i=0; x[i] < plt.lq; i++) {
            if (x[i] < xlo) {
                /* outlier */
                nlo++;
            } else {
                /* let the lower whisker extend to the
                   least non-outlying value below Q1
                */
                plt.wmin = x[i];
                break;
            }
        }

        /* and upper region */
        for (i=n-1; x[i] > plt.uq; i--) {
            if (x[i] > xhi) {
                /* outlier */
                nhi++;
            } else {
                /* let the upper whisker extend to the
                   greatest non-outlying value above Q3
                */
                plt.wmax = x[i];
                break;
            }
        }

        if (nlo + nhi > 0) {
            /* build vector of outliers */
            int j = 0;

            plt.outliers = gretl_vector_alloc(nlo + nhi);
            if (plt.outliers != null) {
                for (i=0; i<nlo; i++) {
                    gretl_vector_set(plt.outliers, j++, x[i]);
                }
                for (i=n-nhi; i<n; i++) {
                    gretl_vector_set(plt.outliers, j++, x[i]);
                }
            }
        }
    }
}

/* backward compatibility, for when we were writing NAs into
   gretl boxplot files using "%g"
*/
bool bpna(double x)
{
    return na(x) || x > 1.0e308;
}

void six_numbers(Boxplot * p, int do_mean, Prn * prn)
{
    double[5] vals = [p.min, p.lq, p.median, p.uq, p.max];
    int i;

    if (do_mean) {
        if (bpna(p.mean)) {
            pprintf(prn, "%11s".ptr, "NA".ptr);
        } else {
            pprintf(prn, "%#11.5g".ptr, p.mean);
        }
    }

    for (i=0; i<5; i++) {
        if (bpna(vals[i])) {
            pprintf(prn, "%11s".ptr, "NA".ptr);
        } else {
            pprintf(prn, "%#11.5g".ptr, vals[i]);
        }
    }

    if (p.n > 0) {
        pprintf(prn, "  (n=%d)\n".ptr, p.n);
    } else {
        pputc(prn, '\n');
    }
}

void five_numbers_with_interval(Boxplot * plt, Prn * prn)
{
    char[32] tmp;

    sprintf(tmp.ptr, "%#.5g - %#.5g".ptr, plt.conf[0], plt.conf[1]);
    pprintf(prn, "%#8.5g%#10.5g%#10.5g%17s%#10.5g%#10.5g\n".ptr,
            plt.min, plt.lq, plt.median,
            tmp.ptr, plt.uq, plt.max);
}

void box_stats_leader(PlotGroup * grp, int i, int offset, Prn * prn)
{
    Boxplot * plt = &grp.plots[i];

    if (factorized(grp)) {
        char[32] numstr = "";

    if (grp.strvals != null) {
        snprintf(numstr.ptr, 32, "%s".ptr, grp.strvals[i]);
    } else if (grp.dvals != null) {
            sprintf(numstr.ptr, "%g".ptr, gretl_vector_get(grp.dvals, i));
        } else if (grp.pan_min > 0) {
            sprintf(numstr.ptr, "%d".ptr, grp.pan_min + i);
        }
        pprintf(prn, "%-*s".ptr, offset, numstr.ptr);
    } else if (plt.boolstr != null) {
        pprintf(prn, "%s\n %-*s".ptr, plt.varname.ptr, offset - 1, plt.boolstr);
    } else {
        pprintf(prn, "%-*s".ptr, offset, plt.varname.ptr);
    }
}

int has_mean(PlotGroup * grp)
{
    Boxplot * plt;
    int i;

    for (i=0; i<grp.nplots; i++) {
        plt = &grp.plots[i];
        if (na(plt.mean)) {
            if (plt.n == 1) {
                plt.mean = plt.min;
            }
        }
        if (!na(plt.mean)) {
            return 1;
        }
    }

    return 0;
}

int get_format_offset(PlotGroup * grp)
{
    int L = 6;
    int i, n;

    for (i=0; i<grp.nplots; i++) {
        if (grp.plots[i].boolstr != null) {
            n = cast(int) (strlen(grp.plots[i].boolstr));
        } else {
            n = cast(int) (strlen(grp.plots[i].varname.ptr));
        }
        if (n > L) {
            L = n;
        }
    }

    return L + 2;
}

int get_factorized_offset(PlotGroup * grp)
{
    int base = cast(int) (cast(int) strlen(grp.xlabel.ptr));

    if (grp.strvals != null) {
    int i, len;

    for (i=0; i<grp.nplots; i++) {
        len = cast(int) (strlen(grp.strvals[i]));
        if (len > base) {
        base = len;
        }
    }
    }

    return base + 2;
}

int boxplot_print_stats(PlotGroup * grp, Prn * prn)
{
    int do_mean = 0;
    int offset, pad, i;

    if (factorized(grp)) {
    offset = get_factorized_offset(grp);
    } else {
        offset = get_format_offset(grp);
    }

    if (do_intervals(grp)) {
        pprintf(prn, "%s\n\n".ptr, _("Numerical summary with bootstrapped confidence interval for median").ptr);
        pad = offset + 8;
    } else {
    if (*grp.ylabel.ptr != '\0' && *grp.xlabel.ptr != '\0') {
        pprintf(prn, _("Numerical summary for %s by %s").ptr, grp.ylabel.ptr, grp.xlabel.ptr);
        pputs(prn, "\n\n");
        } else if (*grp.ylabel.ptr != '\0') {
            pprintf(prn, _("Numerical summary for %s").ptr, grp.ylabel.ptr);
            pputs(prn, "\n\n");
        } else {
            pprintf(prn, "%s\n\n".ptr, _("Numerical summary").ptr);
        }
        do_mean = has_mean(grp);
        pad = (do_mean)? (offset + 9) : (offset + 10);
    }

    if (factorized(grp)) {
        pputs(prn, grp.xlabel.ptr);
        pad -= strlen(grp.xlabel.ptr);
    }

    if (do_intervals(grp)) {
        pprintf(prn, "%*s%10s%10s%17s%10s%10s\n".ptr, pad,
                "min".ptr, "Q1".ptr, _("median").ptr,
                /* xgettext:no-c-format */
                _("(90% interval)").ptr,
                "Q3".ptr, "max".ptr);
        for (i=0; i<grp.nplots; i++) {
            box_stats_leader(grp, i, offset, prn);
            five_numbers_with_interval(&grp.plots[i], prn);
        }
    } else {
        if (do_mean) {
            pprintf(prn, "%*s%11s%11s%11s%11s%11s\n".ptr, pad, _("mean").ptr,
                    "min".ptr, "Q1".ptr, _("median").ptr, "Q3".ptr, "max".ptr);
        } else {
            pprintf(prn, "%*s%11s%11s%11s%11s\n".ptr, pad,
                    "min".ptr, "Q1".ptr, _("median").ptr, "Q3".ptr, "max".ptr);
        }
        for (i=0; i<grp.nplots; i++) {
            box_stats_leader(grp, i, offset, prn);
            six_numbers(&grp.plots[i], do_mean, prn);
        }
    }

    return 0;
}

enum ITERS = 560;
enum CONFIDENCE = 0.90;

/* obtain bootstrap estimate of 90% confidence interval
   for the sample median of data series x; return low and
   high values in 'low' and 'high' */

int median_interval(double * x, int n, double * low, double * high)
{
    double * medians, samp;
    double[2] p;
    int i, j, t;
    int err;

    medians = cast(typeof(medians)) malloc(ITERS * (*medians).sizeof);
    if (medians == null) return 1;

    samp = cast(typeof(samp)) malloc(n * (*samp).sizeof);
    if (samp == null) {
        free(medians);
        return 1;
    }

    for (i=0; i<ITERS; i++) {
        /* sample with replacement from x */
        for (j=0; j<n; j++) {
            t = gretl_rand_int_max(n);
            samp[j] = x[t];
        }
        /* find the median of the sample */
        medians[i] = gretl_array_quantile(samp, n, 0.5);
    }

    p[0] = (1.0 - CONFIDENCE) / 2.0;
    p[1] = 1.0 - p[0];

    err = gretl_array_quantiles(medians, ITERS, p.ptr, 2);

    *low = p[0];
    *high = p[1];

    free(samp);
    free(medians);

    return err;
}

void plotgroup_destroy(PlotGroup * grp)
{
    int i;

    if (grp == null) {
        return;
    }

    for (i=0; i<grp.nplots; i++) {
        gretl_matrix_free(grp.plots[i].outliers);
        free(grp.plots[i].boolstr);
    }

    free(grp.title);
    free(grp.x);
    free(grp.plots);
    gretl_matrix_free(grp.dvals);

    if (grp.strvals != null) {
    strings_array_free(grp.strvals, grp.nplots);
    }

    free(grp);
}

void boxplot_init(Boxplot * bp)
{
    bp.mean = NADBL;
    bp.median = NADBL;
    bp.conf[0] = bp.conf[1] = NADBL;
    bp.uq = bp.lq = NADBL;
    bp.max = bp.min = NADBL;
    bp.wmax = bp.wmin = NADBL;
    bp.n = 0;
    bp.varname[0] = '\0';
    bp.boolstr = null;
    bp.outliers = null;
}

int plotgroup_factor_setup(PlotGroup * grp, const(DatasetRaw) * dset)
{
    int T = sample_size(dset);
    int dv = grp.list[2];
    int err = 0;

    grp.dvals = gretl_matrix_values(dset.Z[dv] + dset.t1, T,
                                     Gretlopt.OPT_S, &err);
    if (!err) {
        grp.nplots = gretl_vector_get_length(grp.dvals);
    }

    return err;
}

PlotGroup * plotgroup_new(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, double limit, int pan_min, Gretlopt opt)
{
    PlotGroup * grp;
    int i, n, err = 0;

    grp = cast(typeof(grp)) malloc((*grp).sizeof);
    if (grp == null) {
        return null;
    }

    grp.nplots = 0;
    grp.list = null;
    grp.plots = null;
    grp.x = null;
    grp.dvals = null;
    grp.strvals = null;
    grp.title = null;

    grp.literal = literal;

    if (na(limit)) {
        grp.limit = 1.5;
    } else {
        grp.limit = limit;
    }

    /* for identifying panel groups */
    grp.pan_min = pan_min;

    *grp.xlabel.ptr = '\0';
    *grp.ylabel.ptr = '\0';

    if (dset != null) {
        grp.list = list;
        n = sample_size(dset);
        grp.x = cast(typeof(grp.x)) malloc(n * (*grp.x).sizeof);
        if (grp.x == null) {
            err = GretlError.alloc;
        }
    }

    if (!err) {
        if (opt & Gretlopt.OPT_Z) {
            /* --factorized */
            if (dset != null) {
                strcpy(grp.ylabel.ptr, dset.varname[list[1]]);
                strcpy(grp.xlabel.ptr, dset.varname[list[2]]);
                err = plotgroup_factor_setup(grp, dset);
            }
        } else {
            grp.nplots = list[0];
        }
        if ((opt & Gretlopt.OPT_P) && dset != null) {
            const(char) * vlabel = series_get_label(dset, 1);

            strcpy(grp.ylabel.ptr, vlabel == null ? "" : vlabel);
            strcpy(grp.xlabel.ptr, _("group").ptr);
        }
    }

    if (!err) {
        grp.plots = cast(typeof(grp.plots)) malloc(grp.nplots * (*grp.plots).sizeof);
        if (grp.plots == null) {
            err = GretlError.alloc;
            grp.nplots = 0;
        }
    }

    if (!err) {
        for (i=0; i<grp.nplots; i++) {
            boxplot_init(&grp.plots[i]);
        }
        grp.flags = 0;
        grp.n_bools = 0;
        grp.gmax = grp.gmin = 0.0;
        if (opt & Gretlopt.OPT_B) {
            grp.flags |= BoxFlags.whiskbars;
        }
    }

    if (err) {
        plotgroup_destroy(grp);
        grp = null;
    }

    return grp;
}

int test_for_outliers(PlotGroup * grp)
{
    int i, n = 0;

    for (i=0; i<grp.nplots; i++) {
        if (grp.plots[i].outliers != null) {
            n += gretl_vector_get_length(grp.plots[i].outliers);
        }
    }

    return n;
}

int test_for_bool(PlotGroup * grp)
{
    int i;

    for (i=0; i<grp.nplots; i++) {
        if (grp.plots[i].boolstr != null) {
            return 1;
        }
    }

    return 0;
}

void get_box_y_range(PlotGroup * grp, double gyrange, double * ymin, double * ymax)
{
    if (grp.gmin >= 0 && grp.gmin < 0.05 &&
        grp.gmax <= 1 && grp.gmax > 0.95) {
        /* proportions? */
        *ymin = 0.0;
        *ymax = 1.0;
    } else {
        *ymin = grp.gmin - 0.04 * gyrange;
        *ymax = grp.gmax + 0.04 * gyrange;
    }
}

void do_boxplot_xtics(PlotGroup * grp, Gretlopt opt, FILE * fp)
{
    Boxplot * bp;
    const(char) * tstr;
    int maxtics = 20;
    int ticskip = 1;
    int ticlen = 1;
    int ntics, ilen;
    int i, j;

    if (grp.nplots > maxtics) {
        ticskip = cast(int) ceil(grp.nplots / cast(double) maxtics);
        ntics = cast(int) ceil(grp.nplots / cast(double) ticskip);
    } else {
        ntics = grp.nplots;
    }

    fputs("set xtics (", fp);

    for (i=0, j=0; i<grp.nplots; i+=ticskip, j++) {
        bp = &grp.plots[i];
        tstr = (opt & Gretlopt.OPT_Z)? bp.boolstr : bp.varname.ptr;
        ilen = cast(int) (strlen(tstr));
        if (ilen > ticlen) {
            ticlen = ilen;
        }
        fprintf(fp, "\"%s\" %d".ptr, tstr, i+1);
        if (j < ntics - 1) {
            fputs(", ", fp);
        } else {
            fputs(") scale 0.0\n", fp);
        }
    }

    if (ntics * ticlen > 80) {
        /* try to avoid xtic mash-up */
        fputs("set xtics rotate by -45\n", fp);
    }

    fputs("set xtics nomirror\n", fp);
}

void do_boxplot_labels(PlotGroup * grp, int anybool, double gyrange, FILE * fp)
{
    Boxplot * bp;
    double h = gyrange / 12.5; /* was 20 */
    double lpos = grp.gmin - h;
    int i;

    fputs("set noxtics\n", fp);
    fprintf(fp, "set bmargin %d\n".ptr, anybool ? 4 : 3);

    for (i=0; i<grp.nplots; i++) {
        bp = &grp.plots[i];
        fprintf(fp, "set label \"%s\" at %d,%g center\n".ptr,
                bp.varname.ptr, i+1, lpos);
        if (bp.boolstr != null) {
            /* FIXME positioning? */
            fprintf(fp, "set label \"%s\" at %d,%g center\n".ptr,
                    bp.boolstr, i+1, lpos - .8 * h);
        }
    }
}

int write_gnuplot_boxplot(PlotGroup * grp, Gretlopt opt)
{
    FILE * fp = null;
    Boxplot * bp;
    int np = grp.nplots;
    double gyrange;
    double ymin, ymax;
    int anybool, n_outliers;
    int lwidth;
    int i, err = 0;

    fp = open_plot_input_file(PlotType.boxplots, cast(GptFlags) 0, &err);

    if (err) {
        return err;
    }

    anybool = test_for_bool(grp);
    n_outliers = test_for_outliers(grp);
    lwidth = np > 30 ? 1 : 2;

    gyrange = grp.gmax - grp.gmin;
    get_box_y_range(grp, gyrange, &ymin, &ymax);

    if (np > 1 && (opt & (Gretlopt.OPT_Z | Gretlopt.OPT_P))) {
        set_use_xtics(grp);
    }

    if (grp.pan_min > 0) {
        fprintf(fp, "# panel_min %d\n".ptr, grp.pan_min);
    }

    if (grp.title != null) {
        fprintf(fp, "set title \"%s\"\n".ptr, grp.title);
    }

    fputs("set nokey\n", fp);
    fputs("set datafile missing \"?\"\n", fp);

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [0:%d]\n".ptr, np + 1);
    fprintf(fp, "set yrange [%g:%g]\n".ptr, ymin, ymax);

    if (*grp.ylabel.ptr != '\0') {
        fprintf(fp, "set ylabel \"%s\"\n".ptr, grp.ylabel.ptr);
    }
    if (*grp.xlabel.ptr != '\0') {
        fprintf(fp, "set xlabel \"%s\"\n".ptr, grp.xlabel.ptr);
    }

    if (use_xtics(grp)) {
        do_boxplot_xtics(grp, opt, fp);
    } else if (grp.nplots > 1 || anybool) {
        do_boxplot_labels(grp, anybool, gyrange, fp);
    } else {
        fputs("set noxtics\n", fp);
    }

    fputs("set boxwidth 0.4 absolute\n", fp);

    if (grp.literal != null) {
        print_gnuplot_literal_lines(grp.literal, fp);
    }

    fputs("plot \\\n", fp);
    /* quartiles and whiskers */
    fprintf(fp, "'-' using 1:3:2:5:4 w candlesticks lw %d notitle%s, \\\n".ptr, lwidth, (whiskerbars(grp))?
            " whiskerbars 0.5".ptr : "".ptr);
    /* median */
    fputs("'-' using 1:2:2:2:2 w candlesticks lt 1 notitle", fp);

    if (show_mean(grp) || do_intervals(grp) || n_outliers > 0) {
        /* continue plot lines array */
        fputs(", \\\n", fp);
    } else {
        fputc('\n', fp);
    }

    if (show_mean(grp)) {
        /* plot the mean as point */
        fputs("'-' using 1:2 w points lt 2 pt 1 notitle", fp);
    } else if (do_intervals(grp)) {
        /* upper and lower bounds of median c.i. */
        fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle, \\\n", fp);
        fputs("'-' using 1:2:2:2:2 w candlesticks lt 0 notitle", fp);
    }

    if (n_outliers > 0) {
        fputs(", \\\n", fp);
        fprintf(fp, "'-' using 1:2 w points lt 0 pt 7 ps 0.25 notitle\n".ptr);
    } else {
        fputc('\n', fp);
    }

    /* data for quartiles and whiskers, plus plot->n */
    for (i=0; i<np; i++) {
        bp = &grp.plots[i];
        if (na(bp.wmin)) bp.wmin = bp.min;
        if (na(bp.wmax)) bp.wmax = bp.max;
        fprintf(fp, "%d %g".ptr, i+1, bp.wmin);
        if (na(bp.lq)) {
            fputs(" ?", fp);
        } else {
            fprintf(fp, " %g".ptr, bp.lq);
        }
        if (na(bp.uq)) {
            fputs(" ?", fp);
        } else {
            fprintf(fp, " %g".ptr, bp.uq);
        }
        fprintf(fp, " %g %d\n".ptr, bp.wmax, bp.n);
    }
    fputs("e\n", fp);

    /* data for median */
    for (i=0; i<np; i++) {
        bp = &grp.plots[i];
        if (na(bp.median)) {
            fprintf(fp, "%d ?\n".ptr, i+1);
        } else {
            fprintf(fp, "%d %g\n".ptr, i+1, bp.median);
        }
    }
    fputs("e\n", fp);

    if (show_mean(grp)) {
        /* data for mean */
        for (i=0; i<np; i++) {
            bp = &grp.plots[i];
            if (na(bp.mean)) {
                fprintf(fp, "%d ?\n".ptr, i+1);
            } else {
                fprintf(fp, "%d %g\n".ptr, i+1, bp.mean);
            }
        }
        fputs("e\n", fp);
    } else if (do_intervals(grp)) {
        /* data for lower median bound */
        for (i=0; i<np; i++) {
            bp = &grp.plots[i];
            fprintf(fp, "%d %g\n".ptr, i+1, bp.conf[0]);
        }
        fputs("e\n", fp);
        /* data for upper median bound */
        for (i=0; i<np; i++) {
            bp = &grp.plots[i];
            fprintf(fp, "%d %g\n".ptr, i+1, bp.conf[1]);
        }
        fputs("e\n", fp);
    }

    if (n_outliers > 0) {
        int j, nout;

        fprintf(fp, "# auxdata %d 2\n".ptr, n_outliers);
        for (i=0; i<np; i++) {
            bp = &grp.plots[i];
            if (bp.outliers != null) {
                nout = gretl_vector_get_length(bp.outliers);
                for (j=0; j<nout; j++) {
                    fprintf(fp, "%d %.10g\n".ptr, i+1,
                            gretl_vector_get(bp.outliers, j));
                }
            }
        }
        fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

int transcribe_array_factorized(PlotGroup * grp, int i, const(DatasetRaw) * dset)
{
    const(double) * y = dset.Z[grp.list[1]];
    const(double) * d = dset.Z[grp.list[2]];
    int t, n = 0;

    for (t=dset.t1; t<=dset.t2; t++) {
        if (d[t] == grp.dvals.val[i] && !na(y[t])) {
            grp.x[n++] = y[t];
        }
    }

    return n;
}

int factorized_boxplot_check(const(int) * list, const(DatasetRaw) * dset)
{
    int err = 0;

    if (list[0] != 2) {
        /* exactly two variables needed */
        err = GretlError.data;
    } else {
        /* the second of which must be discrete */
        int v2 = list[2];

        if (!accept_as_discrete(dset, v2, 0)) {
            gretl_errmsg_set(_("You must supply two variables, the second of which is discrete").ptr);
            err = GretlError.data;
        }
    }

    return err;
}

/* build the boxplot group data structure, then pass it to
   write_gnuplot_boxplot() to write the command file */

int real_boxplots(const(int) * list, const(DatasetRaw) * dset, const(char) * literal, int pan_min, Gretlopt opt)
{
    PlotGroup * grp;
    SeriesTable * st = null;
    double lim = NADBL;
    int i, k, np;
    int err = 0;

    if (opt & Gretlopt.OPT_L) {
        /* we should have an explicit outlier limit */
        lim = get_optval_double(GretlCmdIndex.bxplot, Gretlopt.OPT_L, &err);
        if (err) {
            return err;
        } else if (na(lim) || lim < 0) {
            return GretlError.invarg;
        }
    }

    if (!(opt & Gretlopt.OPT_P)) {
        /* OPT_P enforces a by-group panel plot */
        for (i=1; i<=list[0]; i++) {
            int v = list[i];

            if (gretl_isconst(dset.t1, dset.t2, dset.Z[v])) {
                gretl_errmsg_sprintf(_("%s is a constant").ptr, dset.varname[v]);
                return GretlError.data;
            }
        }
    }

    if (opt & Gretlopt.OPT_Z) {
        /* --factorized */
        err = factorized_boxplot_check(list, dset);
        if (err) {
            return err;
        } else {
            st = series_get_string_table(dset, list[2]);
        }
    }

    grp = plotgroup_new(list, literal, dset, lim, pan_min, opt);
    if (grp == null) {
        return GretlError.alloc;
    }

    if (*grp.ylabel.ptr != '\0' && *grp.xlabel.ptr != '\0') {
        grp.title = gretl_strdup_printf(_("Distribution of %s by %s").ptr,
                                         grp.ylabel.ptr, grp.xlabel.ptr);
    } else if (list[0] == 1) {
        /* doing a single plot */
        const(char) * s = series_get_label(dset, list[1]);
        char[VNAMELEN] yname;

        if (s != null && strncmp(s, "residual for ", 13) == 0 &&
            gretl_scan_varname(s + 13, yname.ptr) == 1) {
            grp.title =
                gretl_strdup_printf(_("Regression residuals (= observed - fitted %s)").ptr,
                                    yname.ptr);
        } else {
            grp.title = gretl_strdup(plotname(dset,list[1], 1));
        }
    }

    if (opt & Gretlopt.OPT_O) {
        set_do_intervals(grp);
    }

    np = grp.nplots;
    k = 0;

    for (i=0; i<np && !err; i++) {
        Boxplot * plot;
        const(char) * vname;
        int n;

        if (opt & Gretlopt.OPT_Z) {
            n = transcribe_array_factorized(grp, i, dset);
            vname = dset.varname[list[1]];
        } else {
            n = transcribe_array(grp.x, dset.Z[list[i+1]], dset);
            vname = dset.varname[list[i+1]];
        }

        if (n < 2) {
            /* less than two observations */
            if (grp.nplots == 1) {
                /* we'll be out of data */
                err = GretlError.data;
                break;
            } else if (n == 0) {
                /* discard this particular plot */
                grp.nplots -= 1;
                continue;
            }
        }

        plot = &grp.plots[k]; /* note: k rather than i */
        qsort(grp.x, n, (*grp.x).sizeof, &gretl_compare_doubles);
        quartiles_etc(grp.x, n, plot, grp.limit);

        plot.n = n;

        if (k == 0 || plot.min < grp.gmin) {
            grp.gmin = plot.min;
        }
        if (k == 0 || plot.max > grp.gmax) {
            grp.gmax = plot.max;
        }
        if (do_intervals(grp)) {
            if (median_interval(grp.x, n, &plot.conf[0], &plot.conf[1])) {
                plot.conf[0] = plot.conf[1] = NADBL;
                unset_do_intervals(grp);
            }
        }

        strcpy(plot.varname.ptr, vname);

        if (opt & Gretlopt.OPT_Z) {
            double di = grp.dvals.val[i];

            if (st != null) {
                plot.boolstr = gretl_strdup(series_table_get_string(st, di));
            } else {
                plot.boolstr = gretl_strdup_printf("%g", di);
            }
            grp.n_bools += 1;
        }

        k++; /* increment the actually used plot index */
    }

    if (!err) {
        if (!do_intervals(grp)) {
            set_show_mean(grp);
        }
        err = write_gnuplot_boxplot(grp, opt);
    }

    plotgroup_destroy(grp);

    return err;
}

/* boxplots using a single selected series, by panel group */

int panel_group_boxplots(int vnum, const(DatasetRaw) * dset, const(char) * literal, Gretlopt opt)
{
    DatasetRaw * gdset;
    int nunits, T = dset.pd;
    int * list = null;
    int nv, u1, i, t, s;
    int err = 0;

    nunits = panel_sample_size(dset);
    nv = nunits + 1;
    u1 = 1 + dset.t1 / T;

    gdset = create_auxiliary_dataset(nv, T, Gretlopt.OPT_NONE);
    if (gdset == null) {
        return GretlError.alloc;
    }

    list = gretl_consecutive_list_new(1, nv - 1);
    if (list == null) {
        destroy_dataset(gdset);
        return GretlError.alloc;
    }

    /* record the name of the original variable */
    series_set_label(gdset, 1, dset.varname[vnum]);

    for (i=0; i<nunits; i++) {
        sprintf(gdset.varname[i+1], "%d".ptr, u1 + i);
        s = dset.t1 + i * T;
        for (t=0; t<T; t++) {
            gdset.Z[i+1][t] = dset.Z[vnum][s++];
        }
    }

    err = real_boxplots(list, gdset, literal, u1, opt);

    destroy_dataset(gdset);
    free(list);

    return err;
}

/*
 * boxplots:
 * @list: list of variables to plot.
 * @dset: dataset struct.
 * @opt: may include OPT_P for a plot by panel group (in which
 * case @list should have just one member), or OPT_Z for a
 * "factorized" plot (@list should have exactly two members),
 * and/or OPT_L to set the limit for outliers.
 *
 * Constructs a boxplot command file for use with gnuplot.
 *
 * Returns: 0 on success, error code on error.
 */

int boxplots(const(int) * list, const(char) * literal, const(DatasetRaw) * dset, Gretlopt opt)
{
    int err;

    if (opt & Gretlopt.OPT_P) {
        if (!multi_unit_panel_sample(dset) ||
            list[0] > 1 || (opt & Gretlopt.OPT_Z)) {
            err = GretlError.badopt;
        } else {
            err = panel_group_boxplots(list[1], dset, literal, opt);
        }
    } else {
        err = real_boxplots(list, dset, literal, 0, opt);
    }

    return err;
}

int get_strvals_from_tics(PlotGroup * grp, const(char) * s)
{
    const(char) * p;
    int i, len;

    grp.strvals = strings_array_new(grp.nplots);
    if (grp.strvals == null) {
        return GretlError.alloc;
    }

    for (i=0; i<grp.nplots; i++) {
        s = strchr(s, '"');
    s++; /* skip the double quote */
    p = s + cast(int) strcspn(s, "\"");
    len = cast(int) (p - s);
    grp.strvals[i] = gretl_strndup(s, len);
    s = p + 1;
    }

    return 0;
}

/* Try parsing numeric factor values out of the xtics line in a gretl
   boxplot file, such as

   set xtics ("0" 1, "1" 2, "4", 3)

   We're looking for numeric values in quotes, but we should be prepared
   to find string values instead (in case the factor variable happens to
   be string-valued).
*/

int get_vals_from_tics(PlotGroup * grp, const(char) * s, int * numvals)
{
    const(char) * p;
    char * tmp;
    double x;
    int len;
    int i, err = 0;

    grp.dvals = gretl_column_vector_alloc(grp.nplots);
    if (grp.dvals == null) {
        return GretlError.alloc;
    }

    for (i=0; i<grp.nplots; i++) {
        s = strchr(s, '"');
        if (s == null) {
        gretl_errmsg_set("Couldn't find factor values");
            err = GretlError.data;
        break;
    } else {
        s++; /* skip the double quote */
        p = s + cast(int) strcspn(s, "\"");
        len = cast(int) (p - s);
        if (len == 0) {
        err = GretlError.data;
        break;
        } else {
        tmp = gretl_strndup(s, len);
        if (numeric_string(tmp)) {
            sscanf(s, "%lf".ptr, &x);
            gretl_vector_set(grp.dvals, i, x);
            *numvals += 1;
        }
        free(tmp);
        }
        s = p + 1;
    }
    }

    if (err || *numvals < grp.nplots) {
        gretl_matrix_free(grp.dvals);
        grp.dvals = null;
    }

    return err;
}

/* read outliers data out of a boxplots file and revise the
   min and max values for the individual plots */

int read_plot_outliers(PlotGroup * grp, char * line, size_t lsize, FILE * fp)
{
    Boxplot * plot;
    double x;
    int i, k, n;
    int err = 0;

    if (sscanf(line + 9, "%d".ptr, &n) != 1) {
        return GretlError.data;
    }

    for (i=0; i<n && !err; i++) {
        if (fgets(line, cast(int) lsize, fp) == null) {
            err = GretlError.data;
        } else if (sscanf(line, "%d %lf".ptr, &k, &x) != 2) {
            err = GretlError.data;
        } else if (k < 1 || k > grp.nplots) {
            err = GretlError.data;
        } else {
            plot = &grp.plots[k-1];
            if (x < plot.min) {
                plot.min = x;
            } else if (x > plot.max) {
                plot.max = x;
            }
        }
    }

    return err;
}

bool is_misscode(const(char) * s)
{
    return strcmp(s, "?") == 0 || strcmp(s, "NaN") == 0;
}

int scan_five(Boxplot * plot, int * j, const(char) * s)
{
    char[16] lq;
    char[16] uq;
    int nf, nobs = 0;

    nf = sscanf(s, "%d %lf %15s %15s %lf %d".ptr, j, &plot.min,
                lq.ptr, uq.ptr, &plot.max, &nobs);

    if (nf >= 5) {
        if (nf == 6) {
            plot.n = nobs;
        }
        if (is_misscode(lq.ptr)) {
            plot.lq = NADBL;
        } else {
            nf--;
            nf += sscanf(lq.ptr, "%lf".ptr, &plot.lq);
        }
        if (is_misscode(uq.ptr)) {
            plot.uq = NADBL;
        } else {
            nf--;
            nf += sscanf(uq.ptr, "%lf".ptr, &plot.uq);
        }
    }

    return nf;
}

int scan_two(int * j, double * x, const(char) * s)
{
    char[16] test;
    int nf;

    nf = sscanf(s, "%d %15s".ptr, j, test.ptr);

    if (nf == 2) {
        if (is_misscode(test.ptr)) {
            *x = NADBL;
        } else {
            nf--;
            nf += sscanf(test.ptr, "%lf".ptr, x);
        }
    }

    return nf;
}

/* read data out of a boxplots file in order to reconstruct
   numerical summary */

int read_plot_data_block(PlotGroup * grp, char * line, size_t lsize, int * block, FILE * fp)
{
    Boxplot * plot;
    int i, j, nf, targ;
    int err = 0;

    for (i=0; i<grp.nplots && !err; i++) {
        plot = &grp.plots[i];
        targ = 2; /* the number of values we should read */
        nf = 0;   /* the number of values actually read */

        if (*block == 0) {
            /* min, Q1, Q3, max */
            targ = 5;
            nf = scan_five(plot, &j, line);
        } else if (*block == 1) {
            /* median */
            nf = scan_two(&j, &plot.median, line);
        } else if (*block == 2) {
            /* mean, or lower bound of c.i. */
            if (do_intervals(grp)) {
                nf = scan_two(&j, &plot.conf[0], line);
            } else {
                nf = scan_two(&j, &plot.mean, line);
            }
        } else if (*block == 3) {
            /* upper bound of c.i. */
            nf = scan_two(&j, &plot.conf[1], line);
        }

        if (nf < targ || j != i + 1) {
            /* data read screwed up */
            err = GretlError.data;
        } else if (fgets(line, cast(int) lsize, fp) == null) {
            /* can't get data line for next plot */
            err = GretlError.data;
        }
    }

    *block += 1;

    return err;
}

void title_to_varname(PlotGroup * grp, const(char) * title)
{
    int len = cast(int) (cast(int) strlen(title));

    if (len < 16 && len == gretl_namechar_spn(title)) {
        strcpy(grp.plots[0].varname.ptr, title);
    }
}

/* Given a gnuplot boxplot file created by gretl, parse out the
   "numerical summary" information.  Note that this requires a
   strictly regimented plot file, so in case of any changes to the
   way such files are printed this function should be checked for
   breakage.
*/

int boxplot_numerical_summary(const(char) * fname, Prn * prn)
{
    PlotGroup * grp = null;
    int factorized = 0;
    int plotlines = 0;
    char[VNAMELEN] yvar;
    char[VNAMELEN] xvar;
    char[128] title;
    char[16] fmt;
    char[16] tfmt;
    char[512] line;
    int parsing = 1;
    int pan_min = 0;
    FILE * fp;
    int n = 0;
    int err = 0;

    fp = gretl_fopen(fname, "r");
    if (fp == null) {
        return GretlError.fopen;
    }

    *yvar.ptr = *xvar.ptr = *title.ptr = '\0';
    sprintf(fmt.ptr, "%%%d[^\\\"]".ptr, VNAMELEN-1);
    sprintf(tfmt.ptr, "%%%d[^\\\"]".ptr, 127);

    /* first pass: count the plots, using labels as heuristic */

    while (fgets(line.ptr, line.sizeof, fp)) {
        if (parsing && !strncmp(line.ptr, "# start literal", 15)) {
            parsing = 0;
        } else if (!parsing && !strncmp(line.ptr, "# end literal", 13)) {
            parsing = 1;
        }
        if (!parsing) {
            continue;
        }
        if (!strncmp(line.ptr, "# panel_min", 11)) {
            pan_min = atoi(line.ptr + 12);
        } else if (!strncmp(line.ptr, "1 ", 2)) {
            /* reached first data block */
            break;
        } else if (!strncmp(line.ptr, "set ylabel ", 11)) {
            sscanf(line.ptr + 12, fmt.ptr, yvar.ptr);
        } else if (!strncmp(line.ptr, "set xlabel ", 11)) {
            sscanf(line.ptr + 12, fmt.ptr, xvar.ptr);
        } else if (!strncmp(line.ptr, "set title ", 10)) {
            sscanf(line.ptr + 11, tfmt.ptr, title.ptr);
        } else if (!strncmp(line.ptr, "set label ", 10)) {
            if (!strncmp(line.ptr + 10, "\"(", 2)) {
                ; /* boolean specifier */
            } else {
                /* should correspond to a plot */
                n++;
            }
        } else if (!strncmp(line.ptr, "'-'", 3)) {
            plotlines++;
        }
    }

    if (n == 0 && *yvar.ptr != '\0' && *xvar.ptr != '\0') {
        /* try getting a count of plots via the first data block */
        while (fgets(line.ptr, line.sizeof, fp)) {
            sscanf(line.ptr, "%d".ptr, &n);
            if (*line.ptr == 'e') {
                /* end of block */
                break;
            }
        }
        if (n > 0) {
            factorized = 1;
        }
    } else if (n == 0 && *title.ptr != '\0') {
        /* plotting a single variable, no labels */
        n = 1;
    }

    if (n == 0) {
        err = GretlError.data;
    } else {
        int[2] list = [n, 0];

        grp = plotgroup_new(list.ptr, null, null, NADBL, pan_min, Gretlopt.OPT_NONE);
        if (grp == null) {
            err = GretlError.alloc;
        }
    }

    if (!err && plotlines > 3) {
        /* note: this is provisional */
        set_do_intervals(grp);
    }

    if (!err) {
        /* second pass: read the data */
        int i = 0, block = 0, outliers = 0;
    int numvals = 0;
        char * label;

        rewind(fp);
        gretl_push_c_numeric_locale();
        parsing = 1;

        while (fgets(line.ptr, line.sizeof, fp) && !err) {
            if (parsing && !strncmp(line.ptr, "# start literal", 15)) {
                parsing = 0;
            } else if (!parsing && !strncmp(line.ptr, "# end literal", 13)) {
                parsing = 1;
            }
            if (!parsing) {
                continue;
            }
            if (factorized && !strncmp(line.ptr, "set xtics (", 11)) {
                err = get_vals_from_tics(grp, line.ptr + 10, &numvals);
        if (grp.dvals == null && !err) {
            err = get_strvals_from_tics(grp, line.ptr + 10);
        }
            } else if (!factorized && !strncmp(line.ptr, "set label ", 10)) {
                label = gretl_quoted_string_strdup(line.ptr + 10, null);
                if (label == null) {
                    err = GretlError.data;
                } else if (*label == '(') {
                    /* should be boolean specifier, following varname */
                    if (i > 0) {
                        grp.plots[i-1].boolstr = label;
                    } else {
                        free(label);
                        err = GretlError.data;
                    }
                } else {
                    strncat(grp.plots[i++].varname.ptr, label, VNAMELEN-1);
                    free(label);
                }
            } else if (!strncmp(line.ptr, "# auxdata", 9)) {
                outliers = 1;
                err = read_plot_outliers(grp, line.ptr, line.sizeof, fp);
            } else if (!strncmp(line.ptr, "1 ", 2)) {
                /* start of a data block */
                err = read_plot_data_block(grp, line.ptr, line.sizeof, &block, fp);
            }
        }

        gretl_pop_c_numeric_locale();

        if (!err && outliers && plotlines == 4) {
            /* we didn't really get median intervals after all */
            for (i=0; i<grp.nplots; i++) {
                grp.plots[i].mean = grp.plots[i].conf[0];
                grp.plots[i].conf[0] = NADBL;
            }
            unset_do_intervals(grp);
        }
    }

    fclose(fp);

    if (!err) {
        if (n == 1 && grp.plots[0].varname[0] == '\0' &&
            *title.ptr != '\0') {
            title_to_varname(grp, title.ptr);
        }
        strcpy(grp.ylabel.ptr, yvar.ptr);
        strcpy(grp.xlabel.ptr, xvar.ptr);
        if (factorized) {
            set_box_factor(grp);
        }
        if (!na(grp.plots[0].conf[0])) {
            set_do_intervals(grp);
        }
        boxplot_print_stats(grp, prn);
    }

    plotgroup_destroy(grp);

    return err;
}

