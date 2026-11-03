/* Port of lib/src/bhhh_max.h: the callback type of bhhh_max() (the prototype of bhhh_max() itself is in
   bhhh_max.d). */

module gretl.bhhh_maxh;


import gretl.gretl_matrixh : GretlMatrixRaw;

alias BHHH_FUNC = extern (C) double function(double *, GretlMatrixRaw *, void *, int, int *);
