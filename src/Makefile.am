AM_FCFLAGS = $(MQC_SPECIFIC_FCFLAGS)
AM_CPPFLAGS = $(MQC_CPPFLAGS)
lib_LIBRARIES = libmqc.a
if WITH_GAUOPEN
libmqc_a_LIBADD = GauOpen/qcmatrix_4.o GauOpen/qcmatrix_8.o GauOpen/qcmatrixio_4.o GauOpen/qcmatrixio_8.o
endif
libmqc_a_SOURCES = mqc_general.F03 mqc_binary.F03 mqc_datastructures.F03 mqc_files.F03 mqc_integrals.F03 mqc_algebra.F03 mqc_algebra2.F03 mqc_molecule.F03 mqc_est.F03 mqc_util.c
if WITH_GAUOPEN
libmqc_a_SOURCES += mqc_matwrapper.F03
else
libmqc_a_SOURCES += mqc_matwrapper_stub.F03
endif
libmqc_a_SOURCES += mqc_gaussian.F03 mqc_FullWavefunction.F03 mqc_interface.F03

build: $(libmqc_a_OBJECTS)
