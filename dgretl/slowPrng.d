/* slowPrng.d: a minimal example of "false sharing" between small GC-allocated objects, and of how much it costs.

   Build and run (no dependencies):

       dmd -O -release -inline slowPrng.d && ./slowPrng        (or: ldc2 -O3 -release slowPrng.d && ./slowPrng)

   What it shows. A random number generator is a small mutable object: here a class with six ints of state that every call to next() reads and writes. A parallel simulation
   gives each task its own generator, so nothing is shared and no locking is needed. The generators are made one after another with `new`, and the allocator puts small objects
   next to each other: several generators land on the same 64-byte cache line. Every write by one thread to its own generator then invalidates that cache line in the cores of
   the threads whose generators are on the same line ("false sharing"), and the threads slow each other down although no data is shared. The work per task is the same for every
   number of threads, so with no sharing the time would stay flat as the number of threads grows.

   Three variants of the same generator are timed with 1, 2, 4, ... tasks, all doing the same fixed amount of work in parallel (std.parallelism):

     Prng        the class written naturally: six ints of state
     PaddedPrng  the same with 64 bytes of padding before and after the state, so that no two objects can share a cache line (this is the workaround)
     local       the state in local variables of the loop (the best case: nothing in memory)

   The program prints the addresses of the first objects, then a table of times in milliseconds (the best of three runs). On an Intel Core i9-13900 (32 hardware threads;
   Linux 6.8) it prints, with dmd -O -release -inline (ldc2 -O3 -release gives the same shape):

       tasks      Prng  PaddedPrng     local  Prng / Padded
           1        37          37        21           1.0x
           2        37          42        22           0.9x
           4       156          60        29           2.6x
           8       218          65        41           3.4x
          16       282          75        44           3.8x
          32       398         129        78           3.1x

   With one task nothing is shared and the classes are the same speed. With 4 or more tasks the natural class is 2.6 to 3.8 times slower than the padded one, up to 10 times slower
   than its own single-task time, although every task does the same work on its own object. In the program this was taken from (parallel streams of a random number generator,
   each a small class, in a simulation) the padding made 16 streams about 3 times faster.

   Two notes. (1) `pragma(inline, false)` on next() stands for what real code looks like (the generator is called from elsewhere, so its state is in memory between calls); if the
   compiler can inline next() into the loop and keep the fields in registers, as ldc2 does without the pragma, there is nothing to share and no slowdown. (2) The object sizes
   (printed by the program) are why neighbours end up on one line: the allocator rounds 40 bytes up to a 48-byte size class.

   The question for the compiler and runtime developers: the padding does not change the program's meaning and costs nothing when nothing is shared, but nobody writes it until
   they have been bitten, and the loss is a large multiple of the run time. Is there a way for the language, the compiler or the GC to avoid the problem automatically? Some
   ideas: allocating class instances (or small objects) cache-line aligned or rounded up in size on request, an `align(64)` that works for class instances (the GC only guarantees
   16 bytes), a GC allocation function that takes an alignment, or thread-local allocation buffers, so that the objects that a thread allocates are kept apart from those of other
   threads.

   (Not shown here: with fresh `core.thread.Thread`s instead of the task pool I measured large slowdowns for BOTH classes in this loop; I do not understand that and have left it
   out so as not to mix two things.)
*/

module slowPrng;

import std.algorithm: min;
import std.datetime.stopwatch: StopWatch, AutoStart;
import std.parallelism: TaskPool, totalCPUs;
import std.range: iota;
import std.stdio: writefln, writeln;

enum callsPerTask = 30_000_000;

// A small generator whose state (six ints) is changed by every call, like the state of a multiple recursive generator.
final class Prng {
    int a, b, c, d, e, f;

    this(int seed) {
        a = seed; b = seed + 1; c = seed + 2; d = seed + 3; e = seed + 4; f = seed + 5;
    }

    pragma(inline, false)               // as in real code, where the generator is called from elsewhere: the state has to live in memory between calls
    int next() {
        const y = a ^ (b << 3) ^ (c >>> 5);
        a = b; b = c; c = y;
        const z = d + e * 7 + f;
        d = e; e = f; f = z;
        return y - z;
    }
}

// The same generator with 64 bytes of padding on both sides of the state: two objects can no longer share a cache line.
final class PaddedPrng {
    ubyte[64] padBefore;
    int a, b, c, d, e, f;
    ubyte[64] padAfter;

    this(int seed) {
        a = seed; b = seed + 1; c = seed + 2; d = seed + 3; e = seed + 4; f = seed + 5;
    }

    pragma(inline, false)               // as in real code, where the generator is called from elsewhere: the state has to live in memory between calls
    int next() {
        const y = a ^ (b << 3) ^ (c >>> 5);
        a = b; b = c; c = y;
        const z = d + e * 7 + f;
        d = e; e = f; f = z;
        return y - z;
    }
}

// One task per generator, run in parallel on a pool of n threads (the calling thread is one of them); returns milliseconds.
long timeParallel(G)(size_t n) {
    G[] gens;

    foreach (i; 0 .. n) gens ~= new G(cast(int) i + 1);          // made one after another, in the main thread
    auto result = new long[n * 16];                              // one slot per task, a cache line apart
    auto pool = new TaskPool(n > 1 ? n - 1 : 0);

    scope (exit) pool.finish(true);
    auto sw = StopWatch(AutoStart.yes);

    foreach (i; pool.parallel(iota(n), 1)) {
        auto g = gens[i];
        long sum = 0;

        foreach (_; 0 .. callsPerTask) sum += g.next();
        result[i * 16] = sum;
    }
    return sw.peek.total!"msecs";
}

// The same loop with the state in local variables.
long timeLocal(size_t n) {
    auto result = new long[n * 16];
    auto pool = new TaskPool(n > 1 ? n - 1 : 0);

    scope (exit) pool.finish(true);
    auto sw = StopWatch(AutoStart.yes);

    foreach (i; pool.parallel(iota(n), 1)) {
        int a = cast(int) i + 1, b = a + 1, c = a + 2, d = a + 3, e = a + 4, f = a + 5;
        long sum = 0;

        foreach (_; 0 .. callsPerTask) {
            const y = a ^ (b << 3) ^ (c >>> 5);
            a = b; b = c; c = y;
            const z = d + e * 7 + f;
            d = e; e = f; f = z;
            sum += y - z;
        }
        result[i * 16] = sum;
    }
    return sw.peek.total!"msecs";
}

long best(alias f, Args...)(Args args) {
    long b = long.max;

    foreach (_; 0 .. 3) b = min(b, f(args));
    return b;
}

void main() {
    // where do consecutively allocated objects go?
    Prng[] p;

    foreach (i; 0 .. 6) p ~= new Prng(cast(int) i + 1);
    writeln("addresses of six Prng objects made one after another (a cache line is 64 bytes):");
    foreach (i, o; p) {
        const addr = cast(size_t) cast(void *) o;
        string note = "";

        if (i > 0) note = (addr / 64 == (cast(size_t) cast(void *) p[i - 1]) / 64) ? "   <- same cache line as the previous object" : "";
        writefln("  object %s at %#x (line %s)%s", i, addr, addr / 64, note);
    }
    writefln("  one instance is %s bytes, of which the state is 24; a PaddedPrng is %s bytes\n", __traits(classInstanceSize, Prng), __traits(classInstanceSize, PaddedPrng));

    writefln("%s hardware threads; every task makes %s calls of next(); best of three runs, milliseconds\n", totalCPUs, callsPerTask);
    writefln("%7s %9s %11s %9s %14s", "tasks", "Prng", "PaddedPrng", "local", "Prng / Padded");
    size_t[] counts = [1, 2, 4, 8, 16];

    if (totalCPUs >= 32) counts ~= 32;
    foreach (n; counts) {
        if (n > totalCPUs) break;
        const a = best!(timeParallel!Prng)(n);
        const b = best!(timeParallel!PaddedPrng)(n);
        const c = best!timeLocal(n);

        writefln("%7s %9s %11s %9s %13.1fx", n, a, b, c, cast(double) a / b);
    }
}
