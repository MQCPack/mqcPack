AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
TESTS = ../setup/set_permissions ../algebra2/testit.sh
check_PROGRAMS = ../setup/set_permissions ../algebra2/algebra_fun
MQC_INSTALLDIR=/usr/bin
SOURCES = ../setup/set_permissions.c ../algebra2/algebra_fun.f03
___setup_set_permissions_SOURCES = ../setup/set_permissions.c
___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.f03
LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}
build: ../setup/set_permissions.$(OBJEXT) ../algebra2/algebra_fun.$(OBJEXT)

