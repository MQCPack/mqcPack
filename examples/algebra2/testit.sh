#! /bin/bash

: "${abs_top_srcdir:=$(cd ../.. && pwd)}"
: "${abs_top_builddir:=$(cd ../.. && pwd)}"

cd "${abs_top_builddir}/examples/algebra2"
./algebra_fun &> outfile
diff -b -B outfile "${abs_top_srcdir}/examples/algebra2/OUTPUT/out"

