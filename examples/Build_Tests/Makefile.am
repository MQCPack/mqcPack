AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
AM_LDFLAGS = $(MQC_LDFLAGS)
AUTOMAKE_OPTIONS = subdir-objects
TESTS_ENVIRONMENT = abs_top_srcdir='$(abs_top_srcdir)' \
  abs_top_builddir='$(abs_top_builddir)'
CLEANFILES = ../algebra2/outfile ../fullci/outfile ../hartreefock/outfile \
  unitTest07-numerical.faf

TESTS = ../algebra2/testit.sh unitTest01 unitTest02 \
  ../../unitTests/src/testit.sh
check_PROGRAMS = ../algebra2/algebra_fun unitTest01 unitTest02 \
  ../../unitTests/src/unitTest_twoERIs

if WITH_GAUOPEN
TESTS += ../fullci/testit.sh ../hartreefock/testit.sh \
  unitTest07
check_PROGRAMS += ../fullci/fullci ../hartreefock/hartreefock unitTest07
if GAUOPEN_FRONTIER
TESTS += unitTest03 unitTest04
check_PROGRAMS += unitTest03 unitTest04
else
TESTS += test_public_gauopen_character.sh
check_PROGRAMS += unitTest03
endif
else
TESTS += unitTest05 test_no_gauopen_matrixfile.sh
check_PROGRAMS += unitTest05 unitTest06
endif

___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.F03
___fullci_fullci_SOURCES = ../fullci/fullci.F03
___hartreefock_hartreefock_SOURCES = ../hartreefock/hartreefock.F03
unitTest01_SOURCES = ../../unitTests/src/unitTest01.f03
unitTest02_SOURCES = ../../unitTests/src/unitTest02.f03
unitTest03_SOURCES = ../../unitTests/src/unitTest03.f03
unitTest04_SOURCES = ../../unitTests/src/unitTest04.f03
unitTest05_SOURCES = ../../unitTests/src/unitTest05.f03
unitTest06_SOURCES = ../../unitTests/src/unitTest06.f03
unitTest07_SOURCES = ../../unitTests/src/unitTest07.f03
______unitTests_src_unitTest_twoERIs_SOURCES = \
  ../../unitTests/src/unitTest_twoERIs.F03
EXTRA_DIST = test_no_gauopen_matrixfile.sh \
  test_public_gauopen_character.sh ../algebra2/testit.sh \
  ../algebra2/OUTPUT/out ../fullci/testit.sh ../fullci/OUTPUT/out \
  ../hartreefock/testit.sh ../hartreefock/OUTPUT/out \
  ../data/MatrixFile/rhf_h2-sto3g.mat ../../unitTests/src/testit.sh \
  ../../unitTests/src/OUTPUT/out

LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}

BUILD_PROGRAM_OBJECTS = ../algebra2/algebra_fun.$(OBJEXT)
if WITH_GAUOPEN
BUILD_PROGRAM_OBJECTS += ../fullci/fullci.$(OBJEXT) \
  ../hartreefock/hartreefock.$(OBJEXT) unitTest07.$(OBJEXT)
endif

build: $(BUILD_PROGRAM_OBJECTS)
