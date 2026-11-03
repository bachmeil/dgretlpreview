/* Port of missing_private.h */

module gretl.missing_privateh;


bool masked(const(char) * m, int t)
{
    return m != null && m[t] == '1';
}
