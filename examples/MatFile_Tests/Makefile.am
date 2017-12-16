AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
lib_LIBRARIES = ../../src/libmqc.a
______src_libmqc_a_SOURCES = ../../src/mqc_gaussian.f03
______src_libmqc_a_LIBADD = ../../src/mqc_general.$(OBJEXT) ../../src/mqc_algebra2.$(OBJEXT) ../../src/mqc_algebra.$(OBJEXT) ../../src/mqc_datastructures.$(OBJEXT) ../../src/mqc_files.$(OBJEXT) ../../src/mqc_molecule.$(OBJEXT) ../../src/mqc_est.$(OBJEXT) ../../src/mqc_util.$(OBJEXT) ../../.other_libs/gauopen/qcmatrixio.$(OBJEXT) ../../.other_libs/gauopen/qcmatrix.$(OBJEXT)

TESTS = ../setup/set_permissions ../algebra2/testit.sh ../tools/testit.sh
check_PROGRAMS = ../setup/set_permissions ../algebra2/algebra_fun ../tools/CheckGauInput
SOURCES = ../setup/set_permissions.c ../algebra2/algebra_fun.f03 ../tools/CheckGauInput.c
___setup_set_permissions_SOURCES = ../setup/set_permissions.c
___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.f03
___tools_CheckGauInput_SOURCES = ../tools/Check.c ../tools/main_CheckGauInput.f03 ../tools/C_To_F_interface.f03

LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}

build: ../../src/mqc_gaussian.$(OBJEXT) ../setup/set_permissions.$(OBJEXT) ../algebra2/algebra_fun.$(OBJEXT) ../tools/CheckGauInput.$(OBJEXT)


