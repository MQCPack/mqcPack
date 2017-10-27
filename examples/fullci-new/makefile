#
# Set the variables for the working directory, module directory, compile
# options, and library locations
#
#LIBDIR = $(mqcroot)/lapack
LIBS   = -llapack -lblas
RunF = gfortran -fdefault-integer-8 -fdefault-real-8 -std=f2003 #-L$(LIBDIR)
#RunF = pgfortran -i8 -r8 -Mallocatable=03 #-L$(mqcroot)/lapack/
#Prof = -g -O0
#
# Specific rules.
#
all: fullci-new 

clean:
	(rm -f *.o *.exe fullci-new)

fullci-new: fullci-new.f03
	$(RunF) $(Prof) -o fullci-new fullci-new.f03 -I$(mqcroot)/bin $(mqcroot)/bin/mqc.a $(LIBS)


