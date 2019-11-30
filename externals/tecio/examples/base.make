# Set to appropriate C++ compilers
           CPP=g++
        MPICPP=mpic++

# Set to appropriate Fortran compilers
            FC=gfortran
         MPIFC=mpif90
        FFLAGS=-fcray-pointer

#
# Determine if this is the Tecio source distribution or a Tecplot installation:
#
ISTECIOSRC := $(shell test -d ../../teciosrc; echo $$?)

ifeq ($(ISTECIOSRC),0)
	TECIOMPILIB=../../teciompisrc/libteciompi.a
	TECIOLIB=../../teciosrc/libtecio.a
    EXTRALIBS=-lstdc++
    EXTRAINCLUDES=-I../../teciosrc
else
    ifeq ($(shell uname -s),Darwin)
	    TECIOLIB=../../../../MacOS/libtecio.dylib
        TECIOMPILIB=../../../../MacOS/libteciompi.dylib
    else
	    TECIOMPILIB=../../../../bin/libteciompi.so
	    TECIOLIB=../../../../bin/libtecio.so
    endif
    FOUND_INSTALLED_LIBSTDCXX_S := $(shell test -f ../../../../bin/sys/libstdc++.so.6 && echo found || echo missing)
    #
    # Note:
    #     On Linux the examples must link against the libstdc++.so.6 that libtecio.so was built
    #     with. For customer installations this is located up and over in the bin directory of the
    #     installation. For internal developer builds, libstdc++.so.6 isn't in the CMake binary
    #     directory and is located with the compiler.
    #
    ifeq ($(FOUND_INSTALLED_LIBSTDCXX_S),found)
        EXTRALIBS=../../../../bin/sys/libstdc++.so.6
    else
        EXTRALIBS=-lstdc++
    endif
    EXTRAINCLUDES=-I../../../../include
endif


#
# (If not needed reset to empty in actual Makefile)
# link libraries
      LINKLIBS=-lpthread -lm

   PATHTOEXECUTABLE=.

all: $(TARGETS)

clean:
	rm -f $(PATHTOEXECUTABLE)/$(EXECUTABLE)

cppbuild:
	$(CPP) $(CPPFILES) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)

mpicppbuild:
	$(MPICPP) -DTECIOMPI $(CPPFILES) $(EXTRAINCLUDES) $(TECIOMPILIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-mpi

fbuild:
	$(FC) $(FFILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-f

f90build:
	$(FC) $(F90FILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOLIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-f90

mpif90build:
	$(MPIFC) -DTECIOMPI $(F90FILES) $(FFLAGS) $(EXTRAINCLUDES) $(TECIOMPILIB) $(EXTRALIBS) $(LINKLIBS) -o $(PATHTOEXECUTABLE)/$(EXECUTABLE)-mpif90
