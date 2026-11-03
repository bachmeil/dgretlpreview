#!/bin/bash
# Faster version of run_tests.sh for LDC: the library is compiled once (in one run, so generated template code is
# shared) and every test is linked against the objects, instead of recompiling ~250 modules per test.
# Usage: testing/run_tests_fast.sh [test_name ...]   (default: every test_*.d here)
# Needs: ldc2, libxml2, lapack, blas and fftw3. Objects are rebuilt when `rm -rf build/fast` is done, or when
# a library source is newer than the stamp.
set -e
cd "$(dirname "$0")"
OUT=build/fast
mkdir -p $OUT/lib $OUT/dg $OUT/exe

LIBSRC="$(ls ../gretl/*.d ../gretl/plugin/*.d ../gretl/minpack/*.d ../cephes/*.d | tr '\n' ' ')"
DGSRC="$(ls ../dgretl/*.d | grep -v -e /testgslprng.d -e /slowPrng.d | tr '\n' ' ')"

newer() { [ ! -e "$1" ] || [ -n "$(find $2 stubs.d -newer "$1" -print -quit)" ]; }

if newer $OUT/lib.stamp "$LIBSRC"; then
    echo "=== compiling the library"
    rm -f $OUT/lib/*.o
    ldc2 -c -g -oq -od=$OUT/lib -I.. $LIBSRC stubs.d
    touch $OUT/lib.stamp
fi
if newer $OUT/dg.stamp "$DGSRC $LIBSRC"; then
    echo "=== compiling dgretl"
    rm -f $OUT/dg/*.o
    ldc2 -c -g -oq -od=$OUT/dg -I.. $DGSRC
    touch $OUT/dg.stamp
fi
ldc2 -c -of=$OUT/stubs_unported.o stubs_unported.d

tests=${@:-$(ls test_*.d | sed 's/\.d$//')}
status=0
failed=""
for t in $tests; do
    echo "=== building $t"
    objs="$OUT/lib/*.o"
    if grep -q '^import dgretl' $t.d; then objs="$objs $OUT/dg/*.o"; fi
    if ! ldc2 -g -I.. -of=$OUT/exe/$t $t.d $objs $OUT/stubs_unported.o \
        -L-lm -L-lxml2 -L-lz -L-llapack -L-lblas -L-lfftw3; then
        echo "BUILD FAILED $t"; failed="$failed $t(build)"; status=1; continue
    fi
    echo "=== running $t"
    $OUT/exe/$t || { failed="$failed $t"; status=1; }
done
echo "=== done; failed:${failed:- none}"
exit $status
