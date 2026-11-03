/* Port of gretl_cmatrix.h: the header holds only function prototypes (they get no D module) and the macro
   complex_scalar, which is a function here. */

module gretl.gretl_cmatrixh;


import gretl.gretl_matrixh : GretlMatrixRaw;

bool complex_scalar(const(GretlMatrixRaw) * m)
{
    return m.is_complex && m.rows == 1 && m.cols == 1;
}
