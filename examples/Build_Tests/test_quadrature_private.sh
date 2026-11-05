#!/bin/sh
# Verify a successful public import before requiring privacy diagnostics.
set -eu
work=$(mktemp -d "${TMPDIR:-/tmp}/mqc-quadrature-private.XXXXXX")
trap 'rm -rf "$work"' EXIT HUP INT TERM
cd "$work"
cat >public.f03 <<'EOF'
program public_grid
use MQC_Integrals, only: quadrature_grid,gridLebedev
use iso_fortran_env, only: int64
implicit none
type(quadrature_grid)::grid
call gridLebedev(grid,6_int64)
end program
EOF
${FC} ${MQC_TEST_FCFLAGS} -I"$abs_top_builddir/src" -c public.f03
cat >private.f03 <<'EOF'
program private_grid
use MQC_Integrals, only: quadrature_grid
implicit none
type(quadrature_grid)::grid
print *,grid%weights
end program
EOF
if ${FC} ${MQC_TEST_FCFLAGS} -I"$abs_top_builddir/src" -c private.f03 >private.log 2>&1; then
  echo 'Private quadrature weights were accessible.'
  exit 1
fi
cat >helper.f03 <<'EOF'
program private_helper
use MQC_Integrals, only: checkRules
implicit none
end program
EOF
if ${FC} ${MQC_TEST_FCFLAGS} -I"$abs_top_builddir/src" -c helper.f03 >helper.log 2>&1; then
  echo 'Private quadrature helper was accessible.'
  exit 1
fi
echo 'Private quadrature storage and helper imports rejected.'
