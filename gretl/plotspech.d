/* Port of plotspec.h */

module gretl.plotspech;


import core.stdc.stdio : FILE;
import gretl.glib;
import gretl.graphingh;
import gretl.gretl_matrixh;
import gretl.libgretlh;

enum GP_BORDER_DEFAULT = -1;
enum LT_AUTO = -2;

enum GpLineStyle {
    none,
    lines,
    points,
    linespoints,
    impulses,
    dots,
    steps,
    boxes,
    errorbars,
    filledcurve,
    candlesticks,
    auto_
}

enum GpKeySpec {
    none    = 0,
    left    = 1 << 0,
    right   = 1 << 1,
    top     = 1 << 2,
    bottom  = 1 << 3,
    outside = 1 << 4,
    plain   = 1 << 5
}

enum GpRange {
    x_range,
    y_range,
    y2_range,
    t_range,
    x2_range,
    cb_range
}

struct GpStyleSpec {
    int id;
    const(char) * name;
    const(char) * trname;
}

enum GpLineFlags {
    user    = 1 << 0,
    boxdata = 1 << 1,
    auxdata = 1 << 2
}

/* information about a line within a gnuplot graph */
struct GptLine {
    int varnum;          /* ID number of variable to plot */
    GpLineStyle style;   /* lines, points, etc. */
    gchar * title;       /* key or legend title */
    gchar * formula;     /* expression to plot (rather than data) */
    char * ustr;         /* custom 'using' string */
    int * mcols;         /* data-matrix columns used */
    float pscale;        /* scale factor for points */
    char[20] rgb;        /* rgb color specification */
    char yaxis;          /* 1 for left, 2 for right */
    int type;            /* 1, 2, ... (style reference) */
    int ptype;           /* point type */
    int dtype;           /* dash type */
    float width;         /* line width, default 1.0 */
    char ncols;          /* number of data columns (0 for formula) */
    float whiskwidth;    /* whiskerbar width (boxplots) */
    char flags;          /* additional options */
}

enum PLOT_LABEL_TEXT_LEN = 31;

enum GpJustCodes {
    left,
    center,
    right
}

struct GptLabel {
    char[PLOT_LABEL_TEXT_LEN + 1] text;
    double[2] pos; /* x,y position on plot */
    guint8 just;   /* justification of text, see GpJustCodes */
}

enum GpArrowFlags {
    head = 1 << 0,
    dots = 1 << 1
}

struct GptArrow {
    double x0, y0, x1, y1;
    int flags;
}

/* "global" information concerning a gnuplot graph specification */
/* struct plotbars_ is private to plotspec.c in the C (opaque in plotspec.h); it is defined here because
   plotspec.d needs the layout */
struct PlotBars {
    int n;         /* number of available start-stop "bar" pairs */
    double t1;     /* current min. value on time axis */
    double t2;     /* current max. value on time axis */
    double ymin;   /* current min. value on y axis */
    double ymax;   /* current max. value on y axis */
    double ** dx;  /* start-stop time pairs */
}

struct GptSpec {
    FILE * fp;
    char[FILENAME_MAX] fname;  /* for gui purposes */
    PlotType code;             /* to deal with FREQ, FCAST... */
    GptFlags flags;            /* bitwise OR of options */
    int heredata;              /* plot file uses "heredoc" for data */
    FitType fit;               /* type of fitted line shown */
    int nobs;                  /* number of observations */
    int okobs;                 /* number of fully valid observations */
    int n_datacols;            /* number of data columns */
    int * heredata_cols;       /* lists of "heredata" columns */
    int pd;                    /* frequency (time series data) */
    int nbars;                 /* number of time-series shaded bars */
    char[MAXDISP] xvarname;    /* name of x variable */
    char[MAXDISP] yvarname;    /* name of y variable */
    gchar *[5] titles;         /* main, x, y, y2, x2 */
    double[2][6] range;        /* axis range specifiers: x, y, y2, t, x2, cb */
    double[3] logbase;         /* axis log-scales base (0 for linear) */
    GpKeySpec keyspec;         /* specification of key (or none) */
    char[16] xfmt;             /* x-axis tic format */
    char[64] xtics;            /* x-axis tic marks */
    char[4] mxtics;            /* minor tics */
    char[16] yfmt;             /* y-axis tic format */
    char[64] ytics;            /* y-axis tic marks */
    char * xticstr;            /* special xtics */
    char * x2ticstr;           /* special x2tics */
    int termtype;              /* gnuplot "terminal" code */
    int n_lines;               /* number of lines */
    int samples;               /* number of samples for parametric plots */
    int border;                /* gnuplot border code */
    char[10] border_lc;        /* border line color */
    int bmargin;               /* bottom margin */
    float boxwidth;            /* when using box style for frequency plots */
    float fillfrac;            /* fraction for "set style fill solid" */
    GptLine * lines;           /* details on individual lines */
    char ** literal;           /* additional commands */
    int n_literal;             /* number of the above */
    char ** multi_xtics;       /* array of xtics strings and positions */
    int n_xtics;               /* number of the above */
    GretlMatrixRaw * data;        /* data to plot */
    char ** markers;           /* data-point markers (not always present) */
    int n_markers;             /* number of such markers */
    GptLabel * labels;         /* labels written onto graph */
    int n_labels;              /* number of the above */
    GptArrow * arrows;         /* arrows drawn onto graph */
    int n_arrows;              /* number of the above */
    int * reglist;             /* regression list for X-Y plot with fitted line */
    GretlMatrixRaw * b_ols;       /* coeffs for linear fit */
    GretlMatrixRaw * b_quad;      /* coeffs for quadratic fit */
    GretlMatrixRaw * b_cub;       /* coeffs for cubic fit */
    GretlMatrixRaw * b_inv;       /* coeffs for inverse fit */
    GretlMatrixRaw * b_log;       /* coeffs for semi-log fit */
    GretlMatrixRaw * b_linlog;    /* coeffs for linear-log fit */
    double scale;              /* for GUI use */
    char * labeled;            /* for GUI use */
    void * ptr;                /* for GUI use */
    PlotBars * bars;           /* for GUI use */
    char * fontstr;            /* for GUI use */
    GretlMatrixRaw * auxdata;     /* for GUI use */
}
