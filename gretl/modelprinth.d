/* Port of modelprint.h */

module gretl.modelprinth;


enum MC_NAMELEN = 32;

struct ModelCoeff {
    double b;
    double se;
    double tval;
    double pval;
    double slope;
    double lo;
    double hi;
    int show_tval;
    int show_pval;
    int df_pval;
    int multi;
    char[MC_NAMELEN] name;
}
