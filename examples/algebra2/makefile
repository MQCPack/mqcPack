#
# Set the variables for the working directory, module directory, compile
# options, and library locations
#
MQCDir     = $(mqcroot)
GAUOPENDIR = $(gauopenroot)
#LIBDIR    = $(mqcroot)/lapack
LIBS       = -llapack -lblas -L$(MQCDir)/bin
MODS       = $(MQCDir)/bin
#F03Flags   = -std=f2008
#RunF      = gfortran -fdefault-integer-8 -fdefault-real-8 
F03Flags   = 
RunF       = pgfortran -i8 -r8 -Mallocatable=03
#Prof      = -g -O0
#
# List all of the modules.
#
MQCMODULESF =                                              \
	$(MQCDir)/src/mqc_general.f03                      \
	$(MQCDir)/src/mqc_algebra.f03                      \
	$(MQCDir)/src/mqc_datastructures.f03               \
	$(MQCDir)/src/mqc_files.f03                        \
	$(MQCDir)/src/mqc_est.f03                          \
	$(MQCDir)/src/mqc_molecule.f03                     \
	$(MQCDir)/src/mqc_gaussian.f03                     

MQCMODULESMOD =                                            \
	$(MQCDir)/src/mqc_general.mod                      \
	$(MQCDir)/src/mqc_algebra.mod                      \
	$(MQCDir)/src/mqc_datastructures.mod               \
	$(MQCDir)/src/mqc_files.mod                        \
	$(MQCDir)/src/mqc_est.mod                          \
	$(MQCDir)/src/mqc_molecule.mod                     \
	$(MQCDir)/src/mqc_gaussian.mod                     

#
# The 'all' rule.
#
all: algebra_fun.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
%.mod: %.f03
	$(RunF) $(F03Flags) -c $*.f03

%.o: %.f03
	$(RunF) $(F03Flags) -c $*.f03

%.o: %.F
	$(RunF) -DUSE_I8 -c $*.F

#
# Generic rule for building general executable program (*.exe) from a standard
# Fortran 2003 source (*.f03) file.
#
%.exe: %.f03 $(MQCDir)/bin/mqc.a
	$(RunF) $(LIBS) $(Prof) -I$(MQCDir)/bin -o $*.exe $*.f03 $(MQCDir)/bin/mqc.a

#
# Rules for building GauOpen object files. These rules take into account the
# compiler chosen for the rest of MQC building in this makefile and build
# qcmatrix.o and qcmatrixio.o from GAUOPENDIR. These object files are then
# moved to the bin directory here so other build versions of the GauOpen object
# files can be built for other purposes.
#
gauopen:
	echo Hrant mqcdir at 1 is $(MQCDir);               \
	(cd $(GAUOPENDIR);                                 \
	make -f $(MQCDir)/makefile qcmatrix.o;             \
	echo Hrant mqcdir at 2 is $(MQCDir);               \
	ar rlv $(MQCDir)/bin/mqc.a qcmatrix.o;             \
	make -f $(MQCDir)/makefile qcmatrixio.o;           \
	ar rlv $(MQCDir)/bin/mqc.a qcmatrixio.o;           \
	)

#
# Rules for building the MQC modules and filling up the bin directory.
#
mqcmods: $(MQCMODULESF)
	(cd $(MQCDir)/src;                                 \
	make -f $(MQCDir)/makefile mqc_general.mod;        \
	make -f $(MQCDir)/makefile mqc_algebra.mod;        \
	make -f $(MQCDir)/makefile mqc_datastructures.mod; \
	make -f $(MQCDir)/makefile mqc_files.mod;          \
	make -f $(MQCDir)/makefile mqc_est.mod;            \
	make -f $(MQCDir)/makefile mqc_molecule.mod;       \
	make -f $(MQCDir)/makefile mqc_gaussian.mod;       \
	)	

mqcbin: cleanmqcbin mqcmods
	(cd $(MQCDir)/bin;                                 \
	ln -s $(MQCDir)/src/mqc_general.mod .;             \
	ln -s $(MQCDir)/src/mqc_general.o   .;             \
	ln -s $(MQCDir)/src/mqc_algebra.mod .;             \
	ln -s $(MQCDir)/src/mqc_algebra.o   .;             \
	ln -s $(MQCDir)/src/mqc_datastructures.mod .;      \
	ln -s $(MQCDir)/src/mqc_datastructures.o   .;      \
	ln -s $(MQCDir)/src/mqc_files.mod .;               \
	ln -s $(MQCDir)/src/mqc_files.o   .;               \
	ln -s $(MQCDir)/src/mqc_est.mod .;                 \
	ln -s $(MQCDir)/src/mqc_est.o   .;                 \
	ln -s $(MQCDir)/src/mqc_molecule.mod .;            \
	ln -s $(MQCDir)/src/mqc_molecule.o   .;            \
	ln -s $(MQCDir)/src/mqc_gaussian.mod .;            \
	ln -s $(MQCDir)/src/mqc_gaussian.o   .;            \
	ar rlv mqc.a $(MQCDir)/src/*.o;                    \
	)

#
# Rules for building test programs.
#
test1: test1.f03 $(MQCDir)/bin/mqc.a
	$(RunF) $(LIBS) -I$(MQCDir)/bin -o test1.exe test1.f03 $(MQCDir)/bin/mqc.a
#
# Various clean rules.
#
clean: cleanmqcbin cleanmqcsrc

cleanmqcbin:
	(cd $(MQCDir)/bin ; rm -f *.o *.mod *.a)

cleanmqcsrc:
	(cd $(MQCDir)/src ; rm -f *.o *.mod)

cleantest:
	(rm -f *.o *.mod *.a *.exe)

