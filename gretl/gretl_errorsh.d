/* Port of gretl_errors.h */

module gretl.gretl_errorsh;


enum GretlError {
    data = 2,
    singular,  /* 3 */
    df,  /* 4 */
    zero,  /* 5 */
    tss,  /* 6 */
    ess,  /* 7 */
    notimp,  /* 8 */
    unspec,  /* 9 */
    pdwrong,  /* 10 */
    fopen,  /* 11 */
    alloc,  /* 12 */
    eqn,  /* 13 */
    unkvar,  /* 14 */
    args,  /* 15 */
    olsonly,  /* 16 */
    invarg,  /* 17 */
    parse,  /* 18 */
    novars,  /* 19 */
    noomit,  /* 20 */
    noadd,  /* 21 */
    adddup,  /* 22 */
    logs,  /* 23 */
    squares,  /* 24 */
    lags,  /* 25 */
    sqrt,  /* 26 */
    high,  /* 27 */
    obs,  /* 28 */
    noconst,  /* 29 */
    badstat,  /* 30 */
    nomerge,  /* 31 */
    noconv,  /* 32 */
    cancel,  /* 33 */
    missdata,  /* 34 */
    nan,  /* 35 */
    nonconf,  /* 36 */
    types,  /* 37 */
    badopt,  /* 38 */
    noident,  /* 39 */
    external,  /* 40 */
    toolong,  /* 41 */
    nodata,  /* 42 */
    notpd,  /* 43 */
    jacobian,  /* 44 */
    toofew,  /* 45 */
    fnest,  /* 46 */
    bounds,  /* 47 */
    funcerr,  /* 48 : error set by function writer */
    stop,  /* 49 : user aborted execution */
    badcatch,  /* 50 : "catch" used where it's not valid */
    cmplx,  /* 51 : complex arguments/operands not supported */
    mixed,  /* 52 : mixed complex/real operands not supported */
    depends,  /* 53 : gfn dependencies not met */
    nonfin,  /* 54 : non-finite values found in matrix argument */
    tsonly,  /* 55 : command needs time-series data */
    db_dup,  /* 56 : duplicate vars found when saving to database */
    ok,  /* 57 : not really an error */
    max  /* 58 */
}

enum GretlWarningCodes {
    gradient = 1,
    genmiss,  /* 2 */
    gennan,  /* 3 */
    max  /* 4 */
}
