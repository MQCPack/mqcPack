AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
mqcmoddir = $(prefix)/mod
MQC_MODULE_FILES = mqc_general.mod mqc_binary.mod \
  mqc_datastructures.mod mqc_files.mod mqc_integrals.mod \
  mqc_algebra.mod mqc_algebra2.mod mqc_molecule.mod mqc_est.mod \
  mqc_nef.mod mqc_matwrapper.mod mqc_gaussian.mod \
  mqc_fullwavefunction.mod
nodist_mqcmod_DATA = $(MQC_MODULE_FILES)
EXTRA_DIST = mqc_general_lapack.F03
lib_LIBRARIES = libmqc.a
CLEANFILES = gauopen_qcmatrix.$(OBJEXT) gauopen_qcmatrixio.$(OBJEXT) \
  $(MQC_MODULE_FILES)
if WITH_GAUOPEN
libmqc_a_LIBADD = gauopen_qcmatrix.$(OBJEXT) gauopen_qcmatrixio.$(OBJEXT)
endif
libmqc_a_SOURCES = mqc_general.F03 mqc_binary.F03 mqc_datastructures.F03 mqc_files.F03 mqc_integrals.F03 mqc_algebra.F03 mqc_algebra2.F03 mqc_molecule.F03 mqc_est.F03 mqc_nef.F03 mqc_util.c
if WITH_GAUOPEN
libmqc_a_SOURCES += mqc_matwrapper_direct.F03
else
libmqc_a_SOURCES += mqc_matwrapper_stub.F03
endif
libmqc_a_SOURCES += mqc_gaussian.F03 mqc_FullWavefunction.F03 mqc_interface.F03

build: $(libmqc_a_OBJECTS)

mqc_binary.$(OBJEXT): mqc_general.$(OBJEXT)
mqc_datastructures.$(OBJEXT): mqc_general.$(OBJEXT) \
  mqc_binary.$(OBJEXT)
mqc_files.$(OBJEXT): mqc_general.$(OBJEXT)
mqc_integrals.$(OBJEXT): mqc_general.$(OBJEXT) \
  mqc_datastructures.$(OBJEXT)
mqc_algebra.$(OBJEXT): mqc_general.$(OBJEXT)
mqc_algebra2.$(OBJEXT): mqc_general.$(OBJEXT)
mqc_molecule.$(OBJEXT): mqc_general.$(OBJEXT) mqc_algebra.$(OBJEXT) \
  mqc_datastructures.$(OBJEXT) mqc_files.$(OBJEXT)
mqc_est.$(OBJEXT): mqc_general.$(OBJEXT) mqc_binary.$(OBJEXT) \
  mqc_algebra.$(OBJEXT) mqc_datastructures.$(OBJEXT) \
  mqc_files.$(OBJEXT)
mqc_nef.$(OBJEXT): mqc_est.$(OBJEXT)
mqc_matwrapper_direct.$(OBJEXT): mqc_general.$(OBJEXT)
mqc_matwrapper_stub.$(OBJEXT): mqc_general.$(OBJEXT)
if WITH_GAUOPEN
MQC_MATWRAPPER_OBJECT = mqc_matwrapper_direct.$(OBJEXT)
else
MQC_MATWRAPPER_OBJECT = mqc_matwrapper_stub.$(OBJEXT)
endif
mqc_gaussian.$(OBJEXT): mqc_general.$(OBJEXT) mqc_files.$(OBJEXT) \
  mqc_integrals.$(OBJEXT) mqc_algebra.$(OBJEXT) \
  mqc_algebra2.$(OBJEXT) mqc_est.$(OBJEXT) mqc_molecule.$(OBJEXT) \
  $(MQC_MATWRAPPER_OBJECT)
mqc_FullWavefunction.$(OBJEXT): mqc_algebra.$(OBJEXT) \
  mqc_datastructures.$(OBJEXT) mqc_files.$(OBJEXT) \
  mqc_est.$(OBJEXT) mqc_gaussian.$(OBJEXT)
mqc_interface.$(OBJEXT): mqc_general.$(OBJEXT)

$(MQC_MODULE_FILES): libmqc.a
	@test -f $@ || { \
	  echo "Expected Fortran module file $@ was not generated." >&2; \
	  exit 1; \
	}

if WITH_GAUOPEN
gauopen_qcmatrix.$(OBJEXT): $(MQC_GAUOPEN_DIR)/qcmatrix.F
	$(FC) $(MQC_GAUOPEN_FCFLAGS) -c -o $@ "$<"

gauopen_qcmatrixio.$(OBJEXT): $(MQC_GAUOPEN_DIR)/qcmatrixio.F
	$(FC) $(MQC_GAUOPEN_FCFLAGS) -c -o $@ "$<"
endif
