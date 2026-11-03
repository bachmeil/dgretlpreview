#!/bin/bash
# Build and run the smoke tests.
# Usage: testing/run_tests.sh [test_name ...]   (default: every test_*.d here)
#
# Needs: dmd, libxml2, lapack, blas and fftw3. (The cephes library is part of the D port now: ../cephes/*.d.)
set -e
cd "$(dirname "$0")"
OUT=build
mkdir -p $OUT

# stubs for libgretl functions that are not ported yet (see make_unported_stubs.sh); a separate compilation, because the
# ported modules declare these functions with their real signatures and the stubs do not
dmd -c -w -of=$OUT/stubs_unported.o stubs_unported.d

tests=${@:-$(ls test_*.d | sed 's/\.d$//')}
status=0
for t in $tests; do
    echo "=== building $t"
    # the direct ports (gretl/, gretl/plugin/, gretl/minpack/, cephes/); the dgretl modules start the library before main and clean it up
    # at exit, so they are only compiled into the tests that import them
    srcs="../gretl/*.d ../gretl/plugin/*.d ../gretl/minpack/*.d ../cephes/*.d"
    # testgslprng.d (a test of the GSL version of the parallel generator: it has a main and needs GSL bindings that are not in this tree) and slowPrng.d (a stand-alone demonstration
    # with a main) are not compiled
    if grep -q '^import dgretl' $t.d; then srcs="$srcs $(ls ../dgretl/*.d | grep -v -e /testgslprng.d -e /slowPrng.d | tr '\n' ' ')"; fi
    dmd -g -w -I.. -of=$OUT/$t $t.d $srcs stubs.d $OUT/stubs_unported.o \
        -L-lm -L-lxml2 -L-lz -L-llapack -L-lblas -L-lfftw3
    echo "=== running $t"
    $OUT/$t || status=1
done
exit $status
