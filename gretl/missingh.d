/* Port of missing.h */

module gretl.missingh;


import std.math : isInfinity, isNaN;

enum NADBL = double.nan;

bool na(double x)
{
    return isNaN(x) || isInfinity(x);
}
