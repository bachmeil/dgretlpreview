/* Port of the enum in lib/src/gretl_midas.h (the remaining declarations
   there are function prototypes, which get no D module). */

module gretl.gretl_midash;


enum MidasType {
    u,        /* "unrestricted MIDAS" */
    nealmon,  /* normalized exponential Almon */
    beta0,    /* normalized beta, last lag 0 */
    betan,    /* normalized beta, non-zero last lag */
    almonp,   /* plain Almon polynomial */
    beta1,    /* as beta0, but 1st coeff clamped */
    max       /* sentinel */
}
