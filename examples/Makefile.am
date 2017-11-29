AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
TESTS = algebra2/testit.sh readTest/testit.sh
check_PROGRAMS = algebra2/algebra_fun readTest/algebra_fun1

SOURCES = algebra2/algebra_fun.f03 readTest/algebra_fun1.f03
algebra2_algebra_fun_SOURCES = algebra2/algebra_fun.f03
readTest_algebra_fun1_SOURCES = readTest/algebra_fun1.f03
LDADD = ../src/libmqc.a
build: algebra2/algebra_fun.$(OBJEXT) readTest/algebra_fun1.$(OBJEXT)
