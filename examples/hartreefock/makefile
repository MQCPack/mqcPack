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
RunF = pgfortran -i8 -r8 -Mallocatable=03 -g -L$(mqcroot)/lapack/ 
#
# List all of the modules.
#
ALLMODULES = $(ModDir)/mqc_files.f03 $(ModDir)/mqc_algebra.f03 $(ModDir)/mqc_gaussian.f03 $(ModDir)/mqc_datastructures.f03
#
# Specific rules.
#
all: hartreefock 

hartreefock: hartreefock.f03 $(ALLMODULES)
	$(RunF) -I$(mqcroot)/bin -o hartreefock ../../bin/*.o hartreefock.f03 $(LIBS)
