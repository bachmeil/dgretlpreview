/* Facade for the port of gretl_func.c: the C file (11,869 lines) is ported in the five modules gretl_func0 (shared head: structs, flags, state),
   gretl_func1 .. gretl_func4 (the functions in their original order); see gretl_func0.d. This module publicly imports them, so that
   `import gretl.gretl_func : name;` works as it would for the C header gretl_func.h. The types and enums of the header are in gretl_funch.d. */

module gretl.gretl_func;

public import gretl.gretl_func0;
public import gretl.gretl_func1;
public import gretl.gretl_func2;
public import gretl.gretl_func3;
public import gretl.gretl_func4;
public import gretl.gretl_funch;
