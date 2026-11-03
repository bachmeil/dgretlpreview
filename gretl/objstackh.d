/* Port of objstack.h */

module gretl.objstackh;


/* Note: the following enumeration may be added to, but the
   existing entries up to GRETL_OBJ_SCALARS must not be
   moved.
*/

enum GretlObjType {
    null_,
    eqn,
    sys,
    var,
    dset,
    info,
    stats,
    corr,
    script,
    notes,
    modtab,
    gpage,
    graph,
    plot,
    text,
    matrix,
    scalars,
    bundle,
    session,
    any,
    max
}

enum SavedObjectFlags {
    in_gui_session = 1 << 0,
    in_named_stack = 1 << 1,
    in_model_table = 1 << 2,
    is_last_model  = 1 << 3,
    in_graph_page  = 1 << 4
}

/* the anonymous enum of object actions */
enum ObjAction {
    none,
    invalid,
    null_,
    show,
    free
}
