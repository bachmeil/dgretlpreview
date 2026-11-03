/* Port of gretl_func.h: the public definitions of the user-function module (gretl.gretl_func0 / gretl_func1 .. gretl_func4, which
   the C has as one file, gretl_func.c).

   The function prototypes of the header are not ported here: the functions are defined in the modules above. The enums are named in
   the D style (members lower case, redundant prefix dropped): DataReq (FN_NEEDS_*), UfunRole (UFUN_*), UfunAttrs (UFUN_PRIVATE ...).
   The string constants NEEDS_TS and so on are enums, and the function-like macro ok_function_return_type is a function. The typedefs
   ufunc, fncall and fnpkg are the structs Ufunc, Fncall and Fnpkg of gretl_func0.d. */

module gretl.gretl_funch;

import gretl.libgretlh : GretlType;

enum DataReq {
    needs_data,  /* needs some (any) sort of dataset */
    needs_ts,    /* function requires time-series data */
    needs_qm,    /* function requires quarterly or monthly data */
    needs_panel, /* function requires panel data */
    nodata_ok    /* function does not require a dataset */
}

enum UfunRole {
    none,
    bundle_print,
    bundle_plot,
    bundle_test,
    bundle_fcast,
    bundle_extra,
    gui_main,
    gui_precheck,
    plot_precheck,
    list_maker,
    r_setup,
    ui_maker,
    max
}

enum UfunAttrs {
    private_  = 1 << 0, /* is private to a package */
    noprint   = 1 << 1, /* offers no printed output */
    menu_only = 1 << 2, /* is GUI-only */
    uses_set  = 1 << 3, /* includes the "set" command */
    has_flow  = 1 << 4, /* includes flow-control (ifs, loops) */
    recurses  = 1 << 5, /* includes a call to itself */
    assign    = 1 << 6  /* GUI usage: return must be assigned */
}

enum NEEDS_TS = "needs-time-series-data";
enum NEEDS_QM = "needs-qm-data";
enum NEEDS_PANEL = "needs-panel-data";
enum NO_DATA_OK = "no-data-ok";

enum FN_NAMELEN = 32;

bool ok_function_return_type(GretlType r)
{
    return r == GretlType.double_ ||
           r == GretlType.series ||
           r == GretlType.matrix ||
           r == GretlType.list ||
           r == GretlType.string_ ||
           r == GretlType.bundle ||
           r == GretlType.strings ||
           r == GretlType.matrices ||
           r == GretlType.bundles ||
           r == GretlType.lists ||
           r == GretlType.arrays ||
           r == GretlType.void_ ||
           r == GretlType.numeric;
}
