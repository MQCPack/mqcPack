AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
TESTS = ../algebra2/testit.sh ../CreateDataList/testit.sh ../tools/testit.sh
check_PROGRAMS = ../algebra2/algebra_fun ../CreateDataList/CheckInput ../tools/FullWavefunctionRW
SOURCES = ../algebra2/algebra_fun.F03 ../CreateDataList/CheckInput.F03 ../tools/FullWavefunctionRW.F03
___algebra2_algebra_fun_SOURCES = ../algebra2/algebra_fun.F03
___CreateDataList_CheckInput_SOURCES = ../CreateDataList/CheckInput.F03
___tools_FullWavefunctionRW_SOURCES = ../tools/FullWavefunctionRW.F03
LDADD = ../../src/libmqc.a ${MQC_LAPACK} ${MQC_BLAS}
build: ../algebra2/algebra_fun.$(OBJEXT) ../CreateDataList/CheckInput.$(OBJEXT) ../tools/FullWavefunctionRW.$(OBJEXT)


