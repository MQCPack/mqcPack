#
# Set the variables for the working directory, module directory, compile
# options, and library locations
#
ModDir = ../../src
LIBDIR = ../../lapack
LIBS   = -llapack -lrefblas
#RunF = gfortran -fdefault-integer-8 -fdefault-real-8 -std=f2003 -L$(LIBDIR)
#RunF = pgfortran -i8 -r8 -Mallocatable=03 
#RunF = pgfortran -i8 -r8 -Mallocatable=03 -L/opt/pgi/osx86-64/14.3/lib/ -L/Volumes/Gaussian/plugin/Modules-v1.1-20150506/lapack/lapack-3.5.0/
RunF = pgfortran -i8 -r8 -Mallocatable=03 -L$(mqcroot)/lapack/ -O2 -g
#
# List all of the modules.
#
ALLMODULES = $(ModDir)/mqc_files.f03 $(ModDir)/mqc_algebra.f03 $(ModDir)/mqc_gaussian.f03 $(ModDir)/mqc_datastructures.f03
#
# Specific rules.
#
all: extendedap 

extendedap: extendedap.f03 $(ALLMODULES)
	$(RunF) -I$(mqcroot)/bin -o extendedap ../../bin/*.o subroutines.f03 ../fullci/mqc_pscf.f03 extendedap.f03 $(LIBS)
