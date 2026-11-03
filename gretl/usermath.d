/* Port of usermat.h: matrix sub-specification types (the function prototypes are not ported here; they are
   declared extern (C) where used, until usermat.c is ported). */

module gretl.usermath;


import core.stdc.limits : INT_MIN;

import gretl.gretl_matrixh : GretlMatrixRaw;

enum MSEL_MAX = INT_MIN;

enum SelType {
    null_,      /* nothing supplied */
    range,      /* integer range p:q provided */
    matrix,     /* selection matrix provided */
    all,        /* comma-separated blank */
    diag,       /* the "diag" dummy constant */
    upper,      /* the "upper" dummy constant */
    lower,      /* the "lower" dummy constant */
    real_,      /* the "real" dummy constant */
    imag,       /* the "imag" dummy constant */
    element,    /* derived: selection is a single element */
    contig,     /* derived: selection is contiguous */
    excl,       /* single exclusion (negative index) */
    single,     /* derived: degenerate range + null */
    str         /* for use with bundles only */
}

bool is_sel_dummy(SelType s)
{
    return s >= SelType.diag && s <= SelType.imag;
}

/* Note SEL_EXCL is flagged only in the case of a single negative
   index. SEL_MATRIX can also do exclusion, if all the elements
   of the vector are negative.
*/

union Msel {
    int[2] range;
    GretlMatrixRaw * m;
    char * str;
}

struct MatrixSubspec {
    int checked;
    SelType ltype, rtype;
    Msel lsel, rsel;
    int * rslice;
    int * cslice;
}

int mspec_get_row_index(const(MatrixSubspec) * m) { return m.lsel.range[0]; }
int mspec_get_col_index(const(MatrixSubspec) * m) { return m.rsel.range[0]; }

void mspec_set_row_index(MatrixSubspec * m, int i) { m.lsel.range[0] = m.lsel.range[1] = i; }
void mspec_set_col_index(MatrixSubspec * m, int j) { m.rsel.range[0] = m.rsel.range[1] = j; }

int mspec_get_element(const(MatrixSubspec) * m) { return m.lsel.range[0]; }
