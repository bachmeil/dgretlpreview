/* Port of gretl_bundle.h */

module gretl.gretl_bundleh;


import gretl.glib;
import gretl.libgretlh;

enum BUNDLE_PRINT = "bundle-print";
enum BUNDLE_PLOT = "bundle-plot";
enum BUNDLE_TEST = "bundle-test";
enum BUNDLE_FCAST = "bundle-fcast";
enum BUNDLE_EXTRA = "bundle-extra";
enum GUI_MAIN = "gui-main";
enum GUI_PRECHECK = "gui-precheck";
enum PLOT_PRECHECK = "plot-precheck";
enum LIST_MAKER = "list-maker";
enum R_SETUP = "R-setup";
enum UI_MAKER = "ui-maker";

enum BundleType {
    plain,
    kalman
}

struct BundledItem {
    GretlType type;
    int is_virtual;
    gpointer data;
    char * note;
    char * key;
}
