# Include local definitions
# the full path to the root TACS directory
SU2_DIR= /home/dthomas/InstalledSoftware/SU2/SU2-feature_pyWrapper
CXX = mpicxx
RM = rm -f

SO_EXT=so

SWIG_SRC = SU2_API.i

SWIG_WRAP_REAL = SU2_APIPYTHON_wrap

SWIG_SO_REAL = _SU2Solver.${SO_EXT}

PYTHON_INCLUDE = -I/usr/include/python2.7 -I/usr/include/python2.7

NUMPY_DIR = /usr/lib/python2.7/dist-packages/numpy/core/include

MPI4PY_DIR = /usr/lib/python2.7/dist-packages/mpi4py/include

SWIG = swig

pySU2_INCLUDE = -I${SU2_DIR}/Common/include \
	-I${SU2_DIR}/SU2_CFD/include \
	-I${SU2_DIR}/SU2_DEF/include \
	-I${SU2_DIR}/SU2_DOT/include \
	-I${SU2_DIR}/SU2_DOT/include \
	-I${SU2_DIR}/SU2_GEO/include \
	-I${SU2_DIR}/SU2_MSH/include \
	-I${SU2_DIR}/SU2_SOL/include \
	-I${SU2_DIR}/SU2_API/include \
	-I${SU2_DIR}/SU2_API/obj


PY_INCLUDE = ${PYTHON_INCLUDE} -I${NUMPY_DIR} -I${MPI4PY_DIR}
SWIG_INCLUDE = -I${MPI4PY_DIR} -I${SU2_DIR}/SU2_API/include

SUBDIR_EXEC = ${SU2_DIR}/SU2_API/_SU2_API

default:
	@echo "Building SU2 Interface"; \
	${MAKE} real; 
	
# List the dependencies
real: ${SWIG_SO_REAL}

EXTRA_CC_FLAGS = -fPIC -O3

METIS_INCLUDE = -I${SU2_DIR}/externals/metis/include
METIS_LIB = ${SU2_DIR}/externals/metis/libmetis.a

PARMETIS_INCLUDE = -I${SU2_DIR}/externals/parmetis/include
PARMETIS_LIB = ${SU2_DIR}/externals/parmetis/libparmetis.a

TECIO_INCLUDE = -I${SU2_DIR}/externals/tecio/include
TECIO_LIB = -I${SU2_DIR}/externals/tecio/libtecio.a

SU2SOLVER_LIB = ${SU2_DIR}/SU2_API/obj/libSU2Solver.a

# Set the linking flags to use
pySU2_LD_FLAGS = -L${SU2_DIR}/SU2_API/obj/ -Wl,-rpath,${SU2_DIR}/SU2_API/obj -lSU2Solver ${PARMETIS_LIB} ${METIS_LIB}

# Set the command line flags to use for compilation
pySU2_CC_FLAGS = ${EXTRA_CC_FLAGS} ${METIS_INCLUDE} ${PARMETIS_INCLUDE} ${TECIO_INCLUDE} ${pySU2_INCLUDE}

SO_LINK_FLAGS=-fPIC -shared 

# Default rule for creating the _wrap.cxx file from the .i file
${SWIG_WRAP_REAL}.cxx: ${SWIG_SRC}
	${SWIG} ${TACS_DEF} -Wall ${SWIG_INCLUDE} -outdir ./ -o ${SWIG_WRAP_REAL}.cxx -c++ -python ${SWIG_SRC} 

#Default rule for compiling the .o files
${SWIG_WRAP_REAL}.o: ${SWIG_WRAP_REAL}.cxx
	${CXX} ${pySU2_CC_FLAGS} ${PY_INCLUDE} -c ${SWIG_WRAP_REAL}.cxx -o ${SWIG_WRAP_REAL}.o

#_%.${SO_EXT}: %.o
#	${CXX} ${SO_LINK_FLAGS} -o _$*.${SO_EXT} $*_wrap.o ${pySU2_LD_FLAGS}

${SWIG_SO_REAL}: ${SWIG_WRAP_REAL}.o
	${CXX} ${SO_LINK_FLAGS} -o ${SWIG_SO_REAL} ${SWIG_WRAP_REAL}.o ${pySU2_LD_FLAGS}

.SECONDARY:

clean: 
	rm -f *.cxx
	rm -f *.${SO_EXT}
	rm -f *.o
	rm -f *.h
	rm -f *.py

