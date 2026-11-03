/* Port of graphing.h */

module gretl.graphingh;


import gretl.gretl_commandsh;
import gretl.glib;
import gretl.plotspech;

enum GptFlags {
    impulses      = 1 << 0, /* use impulses for plotting */
    lines         = 1 << 1, /* force use of lines for plotting */
    resids        = 1 << 2, /* doing residual plot */
    fa            = 1 << 3, /* doing fitted/actual plot */
    factor        = 1 << 4, /* using a factor for separation */
    xyz           = 1 << 5, /* X-Y, controlling for Z */
    fit_omit      = 1 << 6, /* user said don't draw fitted line on graph */
    data_style    = 1 << 7, /* data style is set by user */
    idx           = 1 << 8, /* plot against time or obs index */
    ts            = 1 << 9, /* doing time series plot */
    y2axis        = 1 << 10, /* plot has second y-axis */
    auto_fit      = 1 << 11, /* automatic (OLS) fitted line was added */
    fit_hidden    = 1 << 12, /* autofit line calculated, but suppressed */
    png_output    = 1 << 13, /* output is to PNG file */
    print_markers = 1 << 14, /* print observation markers */
    letterbox     = 1 << 15, /* special format for time series graphs */
    parametric    = 1 << 16, /* gnuplot should be in parametric mode */
    xzeroaxis     = 1 << 17, /* show x = 0 line */
    yzeroaxis     = 1 << 18, /* show y = 0 line */
    mono          = 1 << 19, /* monochrome output */
    grid_y        = 1 << 20, /* display horizontal grid lines */
    grid_x        = 1 << 21, /* display vertical grid lines */
    polar         = 1 << 22, /* plot is in polar mode */
    xl            = 1 << 23, /* large */
    xxl           = 1 << 24, /* extra-large */
    xw            = 1 << 25, /* extra-wide */
    timefmt       = 1 << 26, /* using gnuplot "timefmt" */
    steps         = 1 << 27, /* force steps for plot */
    logy          = 1 << 28, /* log y axis */
    icon          = 1 << 29  /* saving plot "as icon" */
}

/* an extra "command" for use with GUI callback */
enum GP_ASYNC = GretlCmdIndex.nc + 1;

alias gretlRGB = guint32;

enum N_GP_LINETYPES = 8;

enum GP_WIDTH      = 640;
enum GP_HEIGHT     = 480;
enum GP_LB_WIDTH   = 680;
enum GP_LB_HEIGHT  = 400;
enum GP_XL_WIDTH   = 680;
enum GP_XL_HEIGHT  = 510;
enum GP_XXL_WIDTH  = 680;
enum GP_XXL_HEIGHT = 680;
enum GP_XW_WIDTH   = 800;
enum GP_SQ_SIZE    = 480;

enum PlotType {
    regular = 0,
    h_test,
    prob_dist,
    forecast,
    garch,
    freq_simple,
    freq_normal,
    freq_gamma,
    freq_discrete,
    periodogram,
    correlogram,
    cusum,
    multi_basic,
    tri_graph,
    range_mean,
    hurst,
    leverage,
    irfboot,
    kernel,
    roots,
    ellipse,
    multi_irf,
    panel,
    bi_graph,
    many_ts,
    rq_tau,
    factorized,
    boxplots,
    curve,
    qq,
    user,
    xcorrelogram,
    bar,
    stacked_bar,
    plot_3d,
    band,
    heatmap,
    geomap,
    gridplot,
    type_max
}

enum FitType {
    none,
    ols,
    quadratic,
    cubic,
    inverse,
    loess,
    loglin,
    linlog,
    na  /* fit option not applicable */
}

enum TermType {
    none,
    png,
    eps,
    pdf,
    fig,
    tex,
    tkz,
    emf,
    svg,
    htm,
    plt,
    var
}

bool maybe_big_multiplot(PlotType c)
{
    return (c == PlotType.multi_irf ||
            c == PlotType.multi_basic ||
            c == PlotType.panel);
}

bool frequency_plot_code(PlotType c)
{
    return (c == PlotType.freq_simple ||
            c == PlotType.freq_normal ||
            c == PlotType.freq_gamma ||
            c == PlotType.freq_discrete);
}

void set_png_output(GptSpec * p)
{
    p.flags |= GptFlags.png_output;
}

bool get_png_output(const(GptSpec) * p)
{
    return (p.flags & GptFlags.png_output) != 0;
}

void unset_png_output(GptSpec * p)
{
    p.flags &= ~GptFlags.png_output;
}
