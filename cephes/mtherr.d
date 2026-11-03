/* Port of cephes/mtherr.c: the error reporting of the cephes routines. mtherr(name, code) and
   mtherr_with_arg(name, code, arg) print "name <kind> error" (with the offending argument for the second form)
   to stderr and record the code; get_cephes_errno() returns the last recorded code and clears it (the p-value
   code in pvalues.d uses it to detect that a cephes function reported an error); set_cephes_hush(1) silences the
   messages (and the recording) until set_cephes_hush(0).

   Cephes Math Library Release 2.0, April 1987, Copyright 1984, 1987 by Stephen L. Moshier.

   The functions have D linkage. The state is the two __gshared variables of the C (cephes_errno, and the
   function-local static hush, which is zero-initialised). */

module cephes.mtherr;

import core.stdc.stdio : fprintf, stderr;

import cephes.libprobh : CephesError;

private __gshared int cephes_errno = 0;

/* Notice: the order of appearance of the following messages is bound to the error codes defined in
   cephes.libprobh (CephesError) */

private __gshared const(char) *[8] ermsg = [
    "no error".ptr,
    "domain".ptr,       /* error code 1 */
    "singularity".ptr,  /* et seq.      */
    "overflow".ptr,
    "underflow".ptr,
    "total loss of precision".ptr,
    "partial loss of precision".ptr,
    "unknown".ptr
];

/* @name is supposed to be the name of the function in
   which the error occurred; @code is an index into
   the array of error messages above; @arg is the offending
   argument, if @have_arg is non-zero, otherwise it is
   ignored.
*/

private int real_mtherr(const(char) * name, int code, double arg, int have_arg, int quiet)
{
    static __gshared int hush = 0;

    if (name == null) {
        hush = quiet;
        return 0;
    }

    if (hush) {
        return 0;
    }

    fprintf(stderr, "%s ", name);

    if (code <= 0 || code > CephesError.unknown) {
        code = CephesError.unknown;
    }

    cephes_errno = code;

    if (have_arg) {
        fprintf(stderr, "%s error (arg = %g)\n", ermsg[code], arg);
    } else {
        fprintf(stderr, "%s error\n", ermsg[code]);
    }

    return 0;
}

int mtherr_with_arg(const(char) * name, int code, double arg)
{
    return real_mtherr(name, code, arg, 1, 0);
}

int mtherr(const(char) * name, int code)
{
    return real_mtherr(name, code, 0, 0, 0);
}

int get_cephes_errno()
{
    int ret = cephes_errno;

    cephes_errno = 0; /* clear the code */

    return ret;
}

void set_cephes_hush(int s)
{
    real_mtherr(null, 0, 0, 0, s);
}
