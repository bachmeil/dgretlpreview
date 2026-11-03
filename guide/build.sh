#!/bin/bash
# Build the guide's example programs with ldmd2 (the LDC front end that accepts dmd's command line).
# Usage: guide/build.sh kansas            -> guide/kansas
# Needs: ldmd2 (or dmd) and the libxml2, lapack, blas and fftw3 development libraries.
set -e
cd "$(dirname "$0")"
PORT=..
OUT=build
mkdir -p $OUT

# 1. the stand-ins for the libgretl functions that are not ported yet (they abort if they are ever called): compiled on
#    their own, because the ported modules declare these functions with their real signatures and the stubs do not
ldmd2 -c -of=$OUT/stubs_unported.o $PORT/testing/stubs_unported.d

# 2. the program, the whole port, and the small runtime stand-ins in testing/stubs.d
for prog in "$@"; do
    ldmd2 -I$PORT -of=$prog $prog.d $PORT/gretl/*.d $PORT/gretl/plugin/*.d $PORT/gretl/minpack/*.d $PORT/cephes/*.d $(ls $PORT/dgretl/*.d | grep -v -e /testgslprng.d -e /slowPrng.d) $PORT/testing/stubs.d \
        $OUT/stubs_unported.o \
        -L-lm -L-lxml2 -L-lz -L-llapack -L-lblas -L-lfftw3
done
