################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 Python wrapper.
# \author D. Thomas
# \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################


CURRENT_DIR = ${abs_top_builddir}/SU2_PY/pySU2

SO_EXT=so

SWIG_SRC = ${abs_top_srcdir}/SU2_PY/pySU2/pySU2.i
SWIG_WRAP_REAL = SU2_APIPYTHON_wrap
SWIG_SO_REAL = _pysu2.${SO_EXT}

PYTHON_SITE_PACKAGES=$(shell python -c "import site; print(site.getsitepackages()[0])")
#NUMPY_INCLUDE = /usr/lib/python2.7/dist-packages/numpy/core/include
MPI4PY_INCLUDE = ${HOME}/.local/lib/python2.7/site-packages/mpi4py/include \
                 -I${PYTHON_SITE_PACKAGES}/mpi4py/include \
                 -I/Library/Python/2.7/site-packages/mpi4py/include

SWIG = swig

pySU2_INCLUDE = -I${abs_top_srcdir}/Common/include \
	-I${abs_top_srcdir}/SU2_CFD/include

PY_INCLUDE = ${PYTHON_INCLUDE} -I${MPI4PY_INCLUDE}

PY_LIB = ${PYTHON_LIBS} \
         -L${PYTHON_EXEC_PREFIX}/lib \
         -L/usr/lib \
         -L/usr/lib/x86_64-linux-gnu
         

SWIG_INCLUDE = ${PY_INCLUDE}

SUBDIR_EXEC = ${bindir}

all:	
	${MAKE} real; 
	
# List the dependencies
real: ${SWIG_SO_REAL}

EXTRA_CC_FLAGS = -fPIC -O3

#SU2Core lib
SU2SOLVER_LIB = ${abs_top_builddir}/SU2_CFD/obj/libSU2Core.a \
		${abs_top_builddir}/Common/lib/libSU2.a

# Set the linking flags to use
pySU2_LD_FLAGS = ${SU2SOLVER_LIB} @su2_externals_LIBS@ ${PY_LIB}

# Set the command line flags to use for compilation
pySU2_CC_FLAGS = ${CPPFLAGS} ${EXTRA_CC_FLAGS} @su2_externals_INCLUDES@ ${pySU2_INCLUDE}

SO_LINK_FLAGS= -fPIC -shared 

# Default rule for creating the _wrap.cxx file from the .i file
${SWIG_WRAP_REAL}.cxx: ${SU2SOLVER_LIB} ${SWIG_SRC}
	${SWIG} ${CPPFLAGS} -Wall ${SWIG_INCLUDE} -outdir ./ -o ${SWIG_WRAP_REAL}.cxx -c++ -python ${SWIG_SRC} 

#Default rule for compiling the .o files
${SWIG_WRAP_REAL}.o: ${SWIG_WRAP_REAL}.cxx
	${CXX} ${pySU2_CC_FLAGS} ${PY_INCLUDE} -c ${SWIG_WRAP_REAL}.cxx -o ${SWIG_WRAP_REAL}.o

#_%.${SO_EXT}: %.o
#	${CXX} ${SO_LINK_FLAGS} -o _$*.${SO_EXT} $*_wrap.o ${pySU2_LD_FLAGS}

${SWIG_SO_REAL}: ${SWIG_WRAP_REAL}.o
	${CXX} ${SO_LINK_FLAGS} -o ${SWIG_SO_REAL} ${SWIG_WRAP_REAL}.o ${pySU2_LD_FLAGS}

.SECONDARY:

clean: 
	rm -f ${CURRENT_DIR}/*.cxx
	rm -f ${CURRENT_DIR}/*.${SO_EXT}
	rm -f ${CURRENT_DIR}/*.o
	rm -f ${CURRENT_DIR}/*.h
	rm -f ${CURRENT_DIR}/*.py

install: all
	${INSTALL} pysu2.py ${bindir}
	${INSTALL} _pysu2.so ${bindir}

