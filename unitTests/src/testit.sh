#! /bin/bash

: "${abs_top_srcdir:=$(cd ../.. && pwd)}"
: "${abs_top_builddir:=$(cd ../.. && pwd)}"

cd "${abs_top_builddir}/unitTests/src" || exit 1
./unitTest_twoERIs > outfile 2>&1
diff -b -B outfile "${abs_top_srcdir}/unitTests/src/OUTPUT/out"
