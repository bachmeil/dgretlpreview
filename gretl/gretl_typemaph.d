/* Port of gretl_typemap.h: only the NUMERIC_TYPE macro (a function here); the prototypes are not ported. */

module gretl.gretl_typemaph;


import gretl.libgretlh : GretlType;

bool NUMERIC_TYPE(GretlType t)
{
    return (t == GretlType.int_ ||
            t == GretlType.double_ ||
            t == GretlType.series ||
            t == GretlType.useries ||
            t == GretlType.list ||
            t == GretlType.matrix);
}
