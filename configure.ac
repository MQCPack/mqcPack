AC_INIT([MQC_Pack], [26.7.1], [hhratchian@ucmerced.edu])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/mqc_general.F03])

AC_PROG_CPP
AC_PROG_CC([nvc pgcc gcc icc])
AC_PROG_FC([pgfortran nvfortran gfortran ifort])
AC_PROG_RANLIB
AC_LANG([Fortran])

AC_ARG_WITH([gauopen],
  [AS_HELP_STRING([--with-gauopen=DIR],
    [enable MatrixFile support using GauOpen sources in DIR])],
  [AS_CASE([$withval],
    [no], [mqc_with_gauopen=no; MQC_GAUOPEN_DIR=],
    [yes], [AC_MSG_ERROR([--with-gauopen requires a directory])],
    [mqc_with_gauopen=yes; MQC_GAUOPEN_DIR=$withval])],
  [mqc_with_gauopen=no; MQC_GAUOPEN_DIR=])

AC_ARG_WITH([gauopen-integer-bytes],
  [AS_HELP_STRING([--with-gauopen-integer-bytes=BYTES],
    [select the GauOpen integer ABI; currently only 8 is supported @<:@default=8@:>@])],
  [MQC_GAUOPEN_INTEGER_BYTES=$withval],
  [MQC_GAUOPEN_INTEGER_BYTES=8])

AS_CASE([$MQC_GAUOPEN_INTEGER_BYTES],
  [8], [],
  [AC_MSG_ERROR([only the 8-byte GauOpen integer ABI is currently supported])])

AS_IF([test "x$mqc_with_gauopen" = xyes], [
  AS_IF([test ! -d "$MQC_GAUOPEN_DIR"],
    [AC_MSG_ERROR([GauOpen directory not found: $MQC_GAUOPEN_DIR])])
  AS_IF([test ! -f "$MQC_GAUOPEN_DIR/qcmatrixio.F" ||
         test ! -f "$MQC_GAUOPEN_DIR/qcmatrix.F"],
    [AC_MSG_ERROR([GauOpen directory must contain qcmatrixio.F and qcmatrix.F])])
])

AM_CONDITIONAL([WITH_GAUOPEN],
  [test "x$mqc_with_gauopen" = xyes])

MQC_SPECIFIC_FCFLAGS=
MQC_CPPFLAGS=
MQC_LDFLAGS=

AS_CASE([$FC],
  [*pgfortran*|*nvfortran*], [
    MQC_SPECIFIC_FCFLAGS="-r8 -i8 -Mallocatable=03 -module $PWD/src"
    MQC_CPPFLAGS="-DPGI"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -mp"])
  ],
  [*gfortran*], [
    MQC_SPECIFIC_FCFLAGS="-std=f2008 -fdefault-real-8 -fdefault-integer-8 -J$PWD/src"
    MQC_CPPFLAGS="-DGNU"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -fopenmp"])
  ],
  [*ifort*], [
    MQC_SPECIFIC_FCFLAGS="-stand f03 -real-size 64 -integer-size 64 -module $PWD/src"
    MQC_CPPFLAGS="-DINTEL"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -openmp"])
  ],
  [AC_MSG_ERROR([unsupported Fortran compiler: $FC])])

AS_IF([test "x$mqc_with_gauopen" = xyes], [
  MQC_CPPFLAGS="$MQC_CPPFLAGS -DMATRIXFILE"
  AC_CHECK_LIB([blas], [dswap],
    [AC_MSG_NOTICE([using the detected BLAS/LAPACK link environment])],
    [AS_CASE([$FC],
      [*gfortran*],
        [MQC_LDFLAGS="-lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl"],
      [MQC_LDFLAGS="-lmkl_rt -lpthread -lm -ldl"])])
], [
  MQC_CPPFLAGS="$MQC_CPPFLAGS -DFORMCHK"
])

AC_SUBST([MQC_SPECIFIC_FCFLAGS])
AC_SUBST([MQC_CPPFLAGS])
AC_SUBST([MQC_LDFLAGS])
AC_SUBST([MQC_GAUOPEN_DIR])
AC_SUBST([MQC_GAUOPEN_INTEGER_BYTES])

AC_MSG_NOTICE([GauOpen support: $mqc_with_gauopen])
AS_IF([test "x$mqc_with_gauopen" = xyes], [
  AC_MSG_NOTICE([GauOpen source directory: $MQC_GAUOPEN_DIR])
  AC_MSG_NOTICE([GauOpen integer ABI: $MQC_GAUOPEN_INTEGER_BYTES bytes])
])

AC_CONFIG_FILES([Makefile src/Makefile examples/Build_Tests/Makefile])
AC_OUTPUT
