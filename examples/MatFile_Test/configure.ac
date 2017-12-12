AC_INIT([MQC_Pack], [0.1.1], [hhratchian@ucmerced.edu])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([../algebra2/algebra_fun.f03])
AC_PROG_CPP
AC_PROG_CC( [pgcc gcc icc] )
AC_PROG_FC( [pgfortran gfortran ifort] )
AC_PROG_F77
MQC_SPECIFIC_FCFLAGS=;
if `echo $FC 2>&1 | grep 'pgfortran' >/dev/null 2>&1` ; then
  MQC_SPECIFIC_FCFLAGS="-r8 -Mallocatable=03 -module `pwd`/../../src"
else 
if `echo $FC 2>&1 | grep 'gfortran' >/dev/null 2>&1` ; then
  MQC_SPECIFIC_FCFLAGS="-std=f2008 -fdefault-real-8 -J`pwd`/../../src"
else 
  MQC_SPECIFIC_FCFLAGS="ERROR need to add instance of the special flags for this compiler into configure.ac" 
fi
fi
AC_PROG_RANLIB
AC_SUBST(MQC_SPECIFIC_FCFLAGS, $MQC_SPECIFIC_FCFLAGS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
