/* Port of gretl_utils.h */

module gretl.gretl_utilsh;


import std.math : fabs;

bool floateq(double x, double y)
{
    return fabs(x - y) < double.epsilon;
}

bool floatneq(double x, double y)
{
    return fabs(x - y) > double.epsilon;
}

bool floatgt(double x, double y)
{
    return x - y > double.epsilon;
}

bool floatlt(double x, double y)
{
    return y - x > double.epsilon;
}

bool ok_int(double x)
{
    return x <= cast(double) int.max && x >= cast(double) int.min;
}

enum SessionClear {
    all,
    dataset
}

enum BlasVariant {
    unknown,
    netlib,
    atlas,
    openblas,
    mkl,
    veclib,
    blis
}
