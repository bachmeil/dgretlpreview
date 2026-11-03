/* Port of gretl_f2c.h */

module gretl.gretl_f2ch;


/* Minimal support for f2c'd source files: we assume that
   most fortran-isms will be purged from such files,
   including real -> float and doublereal -> double.
*/

alias integer = int;     /* int32_t */
alias logical = int;
alias real_ = float;     /* C "real"; renamed because "real" is a D type */
alias doublereal = double;

/* I/O stuff */

alias flag = int;
alias ftnlen = int;

/* C declares this as `logical (*L_fp)()`, a function pointer with an
   unspecified argument list; D has no such type, so callbacks of this
   kind need a cast when passed. */
alias L_fp = extern (C) logical function();
