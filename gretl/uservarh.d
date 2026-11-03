/* Port of uservar.h: the user-variable types (the function prototypes are not ported here).

   Only the enums, the callback type and UserVar are needed so far;
   the C functions are declared extern (C) where they are used, until uservar.c is ported. */

module gretl.uservarh;


import gretl.libgretlh : GretlType, VNAMELEN;

enum UvarAction {
    add = 1,
    delete_
}

enum UVFlags {
    private_ = 1 << 0,
    shell    = 1 << 1,
    main     = 1 << 2,
    nodecl   = 1 << 3,
    norepl   = 1 << 4
}

alias USER_VAR_FUNC = extern (C) int function(const(char) *, GretlType, int);

/* "user_var": opaque outside uservar.c in the C, whose private header uservar_priv.h defines the struct:
   ported here (as the only module for the type) because geneval.c includes uservar_priv.h */
struct UserVar {
    GretlType type;
    int level;
    UVFlags flags;
    char[VNAMELEN] name;
    void * ptr;
}

void uvar_set_scalar_fast(UserVar * u, double x) { *cast(double *) u.ptr = x; }

double uvar_get_scalar_value(const(UserVar) * u) { return *cast(double *) u.ptr; }
