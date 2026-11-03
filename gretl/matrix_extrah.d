/* Port of matrix_extra.h */

module gretl.matrix_extrah;


enum MMissingCode {
    ok,
    error,
    skip,
    trim
}

enum QuadMethod {
    ghermite = 1,
    legendre,
    laguerre,
    invalid
}
