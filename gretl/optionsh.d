/* Port of options.h */

module gretl.optionsh;


enum OptStatus {
    no_parm = 0,
    accepts_parm,
    needs_parm,
    ambiguous
}
