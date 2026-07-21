AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
AM_LDFLAGS = $(MQC_LDFLAGS)
AUTOMAKE_OPTIONS = subdir-objects

TESTS = ../algebra2/testit.sh unitTest01 unitTest02
check_PROGRAMS = ../algebra2/algebra_fun unitTest01 unitTest02

if WITH_GAUOPEN
TESTS += ../fullci/testit.sh ../hartreefock/testit.sh \
  ../WriteTest/testit.sh unitTest03 unitTest04
check_PROGRAMS += ../fullci/fullci ../hartreefock/hartreefock \
  ../WriteTest/ReadWrite unitTest03 unitTest04
else
TESTS += unitTest05 test_no_gauopen_matrixfile.sh
check_PROGRAMS += unitTest05 unitTest06
endif

___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.F03
___fullci_fullci_SOURCES = ../fullci/fullci.F03
___hartreefock_hartreefock_SOURCES = ../hartreefock/hartreefock.F03
___WriteTest_ReadWrite_SOURCES = ../WriteTest/ReadWrite.F03
unitTest01_SOURCES = ../../unitTests/src/unitTest01.f03
unitTest02_SOURCES = ../../unitTests/src/unitTest02.f03
unitTest03_SOURCES = ../../unitTests/src/unitTest03.f03
unitTest04_SOURCES = ../../unitTests/src/unitTest04.f03
unitTest05_SOURCES = ../../unitTests/src/unitTest05.f03
unitTest06_SOURCES = ../../unitTests/src/unitTest06.f03

EXTRA_DIST = test_no_gauopen_matrixfile.sh

LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}

BUILD_PROGRAM_OBJECTS = ../algebra2/algebra_fun.$(OBJEXT)
if WITH_GAUOPEN
BUILD_PROGRAM_OBJECTS += ../fullci/fullci.$(OBJEXT) \
  ../hartreefock/hartreefock.$(OBJEXT) \
  ../WriteTest/ReadWrite.$(OBJEXT)
endif

build: $(BUILD_PROGRAM_OBJECTS)
