AC_INIT([MQC_Pack], [26.8.0], [hhratchian@ucmerced.edu])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/mqc_general.F03])

AC_PROG_CPP
AC_PROG_CC([nvc pgcc gcc icc])
AC_PROG_FC([pgfortran nvfortran gfortran ifort])
AC_PROG_RANLIB
AC_PROG_GREP
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

AC_ARG_WITH([gauopen-api],
  [AS_HELP_STRING([--with-gauopen-api=API],
    [select the GauOpen API family: auto, public, or frontier @<:@default=auto@:>@])],
  [mqc_gauopen_api_request=$withval],
  [mqc_gauopen_api_request=auto])

AS_CASE([$mqc_gauopen_api_request],
  [auto|public|frontier], [],
  [AC_MSG_ERROR([--with-gauopen-api must be auto, public, or frontier])])

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

MQC_GAUOPEN_API=none
AS_IF([test "x$mqc_with_gauopen" = xyes], [
  mqc_gauopen_detected_api=unknown
  mqc_gauopen_has_close_faf=no
  mqc_gauopen_has_close_matf=no
  mqc_gauopen_has_rd_chbuf=no
  mqc_gauopen_has_wr_chbuf=no
  mqc_gauopen_has_frontier_rind=no
  mqc_gauopen_has_public_rind=no

  AS_IF([$GREP -F -i 'Subroutine Close_FAF' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_close_faf=yes])
  AS_IF([$GREP -F -i 'Subroutine Close_MatF' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_close_matf=yes])
  AS_IF([$GREP -F -i 'Subroutine Rd_ChBuf' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_rd_chbuf=yes])
  AS_IF([$GREP -F -i 'Subroutine Wr_ChBuf' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_wr_chbuf=yes])
  AS_IF([$GREP -F -i 'Rd_RInd(IU,NR,LR,NRLR' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_frontier_rind=yes])
  AS_IF([$GREP -F -i 'Rd_RInd(IU,NR,LR,NTot,LenBuf,LNZ' \
      "$MQC_GAUOPEN_DIR/qcmatrixio.F" >/dev/null 2>&1],
    [mqc_gauopen_has_public_rind=yes])

  AS_IF([test "x$mqc_gauopen_has_close_faf" = xyes &&
         test "x$mqc_gauopen_has_close_matf" = xno &&
         test "x$mqc_gauopen_has_rd_chbuf" = xyes &&
         test "x$mqc_gauopen_has_wr_chbuf" = xyes &&
         test "x$mqc_gauopen_has_frontier_rind" = xyes],
    [mqc_gauopen_detected_api=frontier],
    [AS_IF([test "x$mqc_gauopen_has_close_matf" = xyes &&
            test "x$mqc_gauopen_has_close_faf" = xno &&
            test "x$mqc_gauopen_has_rd_chbuf" = xno &&
            test "x$mqc_gauopen_has_wr_chbuf" = xno &&
            test "x$mqc_gauopen_has_public_rind" = xyes],
       [mqc_gauopen_detected_api=public])])

  AS_CASE([$mqc_gauopen_api_request],
    [auto], [
      AS_IF([test "x$mqc_gauopen_detected_api" = xunknown],
        [AC_MSG_ERROR([unable to identify the GauOpen API family in $MQC_GAUOPEN_DIR; use --with-gauopen-api only after verifying its procedure signatures])])
      MQC_GAUOPEN_API=$mqc_gauopen_detected_api
    ],
    [public|frontier], [
      AS_IF([test "x$mqc_gauopen_detected_api" != xunknown &&
             test "x$mqc_gauopen_detected_api" != "x$mqc_gauopen_api_request"],
        [AC_MSG_ERROR([requested GauOpen API $mqc_gauopen_api_request conflicts with detected API $mqc_gauopen_detected_api])])
      AS_IF([test "x$mqc_gauopen_detected_api" = xunknown],
        [AC_MSG_WARN([forcing unrecognized GauOpen sources to use the $mqc_gauopen_api_request API profile])])
      MQC_GAUOPEN_API=$mqc_gauopen_api_request
    ])
], [
  AS_IF([test "x$mqc_gauopen_api_request" != xauto],
    [AC_MSG_ERROR([--with-gauopen-api requires --with-gauopen=DIR])])
])

AM_CONDITIONAL([WITH_GAUOPEN],
  [test "x$mqc_with_gauopen" = xyes])
AM_CONDITIONAL([GAUOPEN_PUBLIC],
  [test "x$MQC_GAUOPEN_API" = xpublic])
AM_CONDITIONAL([GAUOPEN_FRONTIER],
  [test "x$MQC_GAUOPEN_API" = xfrontier])

MQC_SPECIFIC_FCFLAGS=
MQC_GAUOPEN_FCFLAGS=
MQC_CPPFLAGS=
MQC_LDFLAGS=

AS_CASE([$FC],
  [*pgfortran*|*nvfortran*], [
    MQC_SPECIFIC_FCFLAGS="-r8 -i8 -Mallocatable=03 -module $PWD/src"
    MQC_GAUOPEN_FCFLAGS="-r8 -i8 -DUSE_I8"
    MQC_CPPFLAGS="-DPGI"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -mp"])
  ],
  [*gfortran*], [
    MQC_SPECIFIC_FCFLAGS="-std=f2008 -fdefault-real-8 -fdefault-integer-8 -J$PWD/src"
    MQC_GAUOPEN_FCFLAGS="-std=legacy -fallow-argument-mismatch -fdefault-real-8 -fdefault-integer-8 -DUSE_I8"
    MQC_CPPFLAGS="-DGNU"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -fopenmp"])
  ],
  [*ifort*], [
    MQC_SPECIFIC_FCFLAGS="-stand f03 -real-size 64 -integer-size 64 -module $PWD/src"
    MQC_GAUOPEN_FCFLAGS="-real-size 64 -integer-size 64 -DUSE_I8"
    MQC_CPPFLAGS="-DINTEL"
    AS_IF([test "x$mqc_with_gauopen" = xyes],
      [MQC_SPECIFIC_FCFLAGS="$MQC_SPECIFIC_FCFLAGS -openmp"])
  ],
  [AC_MSG_ERROR([unsupported Fortran compiler: $FC])])

AS_IF([test "x$mqc_with_gauopen" = xyes], [
  MQC_CPPFLAGS="$MQC_CPPFLAGS -DMATRIXFILE"
  AS_CASE([$MQC_GAUOPEN_API],
    [public],
      [MQC_CPPFLAGS="$MQC_CPPFLAGS -DMQC_GAUOPEN_API_PUBLIC"],
    [frontier],
      [MQC_CPPFLAGS="$MQC_CPPFLAGS -DMQC_GAUOPEN_API_FRONTIER"],
    [AC_MSG_ERROR([internal error: no GauOpen API profile selected])])
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
AC_SUBST([MQC_GAUOPEN_FCFLAGS])
AC_SUBST([MQC_CPPFLAGS])
AC_SUBST([MQC_LDFLAGS])
AC_SUBST([MQC_GAUOPEN_DIR])
AC_SUBST([MQC_GAUOPEN_INTEGER_BYTES])
AC_SUBST([MQC_GAUOPEN_API])

AC_MSG_NOTICE([GauOpen support: $mqc_with_gauopen])
AS_IF([test "x$mqc_with_gauopen" = xyes], [
  AC_MSG_NOTICE([GauOpen source directory: $MQC_GAUOPEN_DIR])
  AC_MSG_NOTICE([GauOpen API family: $MQC_GAUOPEN_API])
  AC_MSG_NOTICE([GauOpen integer ABI: $MQC_GAUOPEN_INTEGER_BYTES bytes])
])

AC_CONFIG_FILES([Makefile src/Makefile examples/Build_Tests/Makefile])
AC_OUTPUT
