AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
TESTS = setup/set_permissions tools/testit.sh algebra2/testit.sh
check_PROGRAMS = setup/set_permissions tools/CheckGauInput algebra2/algebra_fun

SOURCES = setup/set_permissions.c tools/CheckGauInput.c algebra2/algebra_fun.f03
setup_set_permissions_SOURCES = setup/set_permissions.c
tools_CheckGauInput_SOURCES = tools/CheckGauInput.c
algebra2_algebra_fun_SOURCES = algebra2/algebra_fun.f03
LDADD = ../src/libmqc.a
build: setup/set_permissions.$(OBJEXT) tools/CheckGauInput.$(OBJEXT) algebra2/algebra_fun.$(OBJEXT)

