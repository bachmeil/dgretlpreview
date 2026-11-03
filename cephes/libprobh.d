/* Port of cephes/libprob.h: the cephes error conditions. The rest of the header is function prototypes
   (bdtr, bdtrc, bdtri, btdtr, chdtr, chdtrc, chdtri, fdtr, fdtrc, fdtri, stdtr, stdtri, ndtr, ndtri, pdtr,
   pdtrc, pdtri, gdtr, gdtrc, igami, cephes_gamma, cephes_lgamma, get_cephes_sgngam, incbet, psi, polrt,
   hyp2f1, get_cephes_errno, set_cephes_hush), which belong to the modules that port the C files defining
   them. */

module cephes.libprobh;

/* cephes error conditions */
enum CephesError {
    domain = 1,  /* argument domain error */
    sing,        /* argument singularity */
    overflow,    /* overflow range error */
    underflow,   /* underflow range error */
    tloss,       /* total loss of precision */
    ploss,       /* partial loss of precision */
    unknown      /* unspecified error */
}
