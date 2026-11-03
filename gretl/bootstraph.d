/* Port of the type of lib/src/bootstrap.h (the prototypes belong to bootstrap.c, not ported yet). The BOOT_METHOD_
   prefix is dropped. */

module gretl.bootstraph;


enum BootMethod {
    residuals = 1,
    pairs,
    wild,
    parametric
}
