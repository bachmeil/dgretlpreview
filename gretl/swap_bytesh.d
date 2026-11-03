/* Port of swap_bytes.h: byte-order reversal helpers.

   The C header is a set of macros (a GNU bswap variant and a portable
   fallback); here they are functions built on core.bitop.bswap, taking
   the target by reference so that the reverse_*() forms (which swap in
   place) are the same call. The macro pairs swap_bytes_X(from, to) become
   `to = swap_bytes_X(from)`.
*/

module gretl.swap_bytesh;


import core.bitop : bswap;

ushort swap_bytes_16(ushort from)
{
    return cast(ushort) (((from >> 8) & 0xff) | ((from & 0xff) << 8));
}

uint swap_bytes_32(uint from)
{
    return bswap(from);
}

double swap_bytes_double(double from)
{
    union U { ulong u64; double d; }
    U u;

    u.d = from;
    u.u64 = bswap(u.u64);
    return u.d;
}

float swap_bytes_float(float from)
{
    union U { uint u32; float f; }
    U u;

    u.f = from;
    u.u32 = bswap(u.u32);
    return u.f;
}

short swap_bytes_short(short from)
{
    return cast(short) swap_bytes_16(cast(ushort) from);
}

int swap_bytes_int(int from)
{
    return cast(int) swap_bytes_32(cast(uint) from);
}

alias swap_bytes_ushort = swap_bytes_16;
alias swap_bytes_uint = swap_bytes_32;

void reverse_ushort(ref ushort x) { x = swap_bytes_16(x); }
void reverse_short(ref short x) { x = swap_bytes_short(x); }
void reverse_uint(ref uint x) { x = swap_bytes_32(x); }
void reverse_int(ref int x) { x = swap_bytes_int(x); }
void reverse_float(ref float x) { x = swap_bytes_float(x); }
void reverse_double(ref double x) { x = swap_bytes_double(x); }
