AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
lib_LIBRARIES = libmqc.a
libmqc_a_LIBADD = GauOpen/qcmatrix_4.o GauOpen/qcmatrix_8.o GauOpen/qcmatrixio_4.o GauOpen/qcmatrixio_8.o
libmqc_a_SOURCES = mqc_general.F03 mqc_binary.F03 mqc_datastructures.F03 mqc_files.F03 mqc_integrals.F03 mqc_algebra.F03 mqc_algebra2.F03 mqc_molecule.F03 mqc_est.F03 mqc_util.c mqc_matwrapper.F03 mqc_gaussian.F03 mqc_FullWavefunction.F03 mqc_interface.F03
build: mqc_general.$(OBJEXT) mqc_binary.$(OBJEXT) mqc_datastructures.$(OBJEXT) mqc_files.$(OBJEXT) mqc_integrals.$(OBJEXT) mqc_algebra.$(OBJEXT) mqc_algebra2.$(OBJEXT) mqc_molecule.$(OBJEXT) mqc_est.$(OBJEXT) mqc_util.$(OBJEXT) mqc_matwrapper.$(OBJEXT) mqc_gaussian.$(OBJEXT) mqc_FullWavefunction.$(OBJEXT) mqc_interface.$(OBJEXT)


