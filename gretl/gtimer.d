/*
   Replacement for glib's g_get_monotonic_time(), using druntime's
   MonoTime, so that the port does not depend on glib for it.

   g_get_monotonic_time() returns the system's monotonic time in
   microseconds, as a 64-bit signed integer. Only the difference between
   two readings is meaningful.

   Used by: gretl_utils.d.
*/

module gretl.gtimer;


import core.time : MonoTime;

import gretl.glib : gint64;

gint64 g_get_monotonic_time()
{
    long ticks = MonoTime.currTime.ticks;
    long tps = MonoTime.ticksPerSecond;

    /* avoid overflow of ticks * 1_000_000 */
    return (ticks / tps) * 1_000_000 + ((ticks % tps) * 1_000_000) / tps;
}

unittest {
    import core.thread : Thread;
    import core.time : msecs;

    gint64 a = g_get_monotonic_time();
    Thread.sleep(20.msecs);
    gint64 b = g_get_monotonic_time();

    assert(b - a >= 15_000 && b - a < 2_000_000);
}
