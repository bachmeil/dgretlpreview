/*
   Replacement for the glib GList (a doubly linked list), written with malloc'd nodes so that it behaves like
   glib's: the list pointer is the first node, an empty list is null, and the functions return the (possibly new)
   head. The functions have the glib names and signatures where the ported code needs them.

   Provided: g_list_append, g_list_prepend, g_list_free, g_list_length, g_list_first, g_list_last, g_list_nth_data,
   g_list_foreach, g_list_sort (stable merge sort, as glib's), g_list_remove, g_list_next (a function standing in
   for the glib macro) and g_list_reverse. Add others here as ported files need them.

   Callbacks use D linkage (GFunc, GCompareFunc): the callers are D functions.

   Used by: gretl_bundle.d.
*/

module gretl.glist;


import core.stdc.stdlib : free, malloc;

import gretl.glib : gconstpointer, gpointer, guint;

struct GList {
    gpointer data;
    GList * next;
    GList * prev;
}

alias GFunc = void function(gpointer data, gpointer user_data);
alias GCompareFunc = int function(gconstpointer a, gconstpointer b);

private GList * new_node(gpointer data)
{
    GList * n = cast(GList *) malloc(GList.sizeof);

    n.data = data;
    n.next = null;
    n.prev = null;
    return n;
}

GList * g_list_next(GList * l)
{
    return l != null ? l.next : null;
}

GList * g_list_last(GList * list)
{
    if (list != null) {
        while (list.next != null) {
            list = list.next;
        }
    }
    return list;
}

GList * g_list_first(GList * list)
{
    if (list != null) {
        while (list.prev != null) {
            list = list.prev;
        }
    }
    return list;
}

GList * g_list_append(GList * list, gpointer data)
{
    GList * n = new_node(data);

    if (list == null) {
        return n;
    }
    GList * last = g_list_last(list);

    last.next = n;
    n.prev = last;
    return list;
}

GList * g_list_prepend(GList * list, gpointer data)
{
    GList * n = new_node(data);

    n.next = list;
    if (list != null) {
        n.prev = list.prev;
        if (list.prev != null) {
            list.prev.next = n;
        }
        list.prev = n;
    }
    return n;
}

void g_list_free(GList * list)
{
    while (list != null) {
        GList * nx = list.next;

        free(list);
        list = nx;
    }
}

guint g_list_length(GList * list)
{
    guint n = 0;

    for (; list != null; list = list.next) {
        n++;
    }
    return n;
}

gpointer g_list_nth_data(GList * list, guint n)
{
    while (n > 0 && list != null) {
        list = list.next;
        n--;
    }
    return list != null ? list.data : null;
}

void g_list_foreach(GList * list, GFunc func, gpointer user_data)
{
    while (list != null) {
        GList * nx = list.next;   /* as glib: the callback may free the node's data */

        func(list.data, user_data);
        list = nx;
    }
}

GList * g_list_remove(GList * list, gconstpointer data)
{
    for (GList * l = list; l != null; l = l.next) {
        if (l.data == data) {
            if (l.prev != null) {
                l.prev.next = l.next;
            } else {
                list = l.next;
            }
            if (l.next != null) {
                l.next.prev = l.prev;
            }
            free(l);
            break;
        }
    }
    return list;
}

GList * g_list_reverse(GList * list)
{
    GList * last = null;

    while (list != null) {
        last = list;
        list = last.next;
        last.next = last.prev;
        last.prev = list;
    }
    return last;
}

/* stable sort: insertion into a new list is O(n^2) but the lists are short (bundle keys); the order of equal
   elements is preserved, as with glib's merge sort */

GList * g_list_sort(GList * list, GCompareFunc compare_func)
{
    GList * sorted = null;

    for (GList * l = list; l != null; l = l.next) {
        GList * n = new_node(l.data);

        if (sorted == null || compare_func(n.data, sorted.data) < 0) {
            n.next = sorted;
            if (sorted != null) {
                sorted.prev = n;
            }
            sorted = n;
        } else {
            GList * p = sorted;

            while (p.next != null && compare_func(n.data, p.next.data) >= 0) {
                p = p.next;
            }
            n.next = p.next;
            n.prev = p;
            if (p.next != null) {
                p.next.prev = n;
            }
            p.next = n;
        }
    }
    g_list_free(list);
    return sorted;
}

unittest
{
    GList * l = null;
    int a = 3, b = 1, c = 2;

    l = g_list_append(l, &a);
    l = g_list_append(l, &b);
    l = g_list_append(l, &c);
    assert(g_list_length(l) == 3 && g_list_nth_data(l, 1) == &b);
    l = g_list_sort(l, function int(gconstpointer x, gconstpointer y) {
        return *cast(const(int) *) x - *cast(const(int) *) y;
    });
    assert(*cast(int *) l.data == 1 && *cast(int *) g_list_last(l).data == 3);
    l = g_list_remove(l, &c);
    assert(g_list_length(l) == 2);
    g_list_free(l);
}

/* stands in for the glib macro */
GList * g_list_previous(GList * list)
{
    return (list != null) ? list.prev : null;
}

alias GDestroyNotifyC = extern (C) void function(gpointer);

void g_list_free_full(GList * list, GDestroyNotifyC free_func)
{
    for (GList * l = list; l != null; l = l.next) {
        if (free_func != null) {
            free_func(l.data);
        }
    }
    g_list_free(list);
}
