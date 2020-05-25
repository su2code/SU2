################################################################################
#
# \file Makefile.am
# \brief Makefile for the SU2 Python wrapper.
# \author D. Thomas
# \version 7.0.4 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

if BUILD_NORMAL
  SWIG_SRC = ${abs_top_srcdir}/SU2_PY/pySU2/pySU2.i
  SWIG_WRAP_REAL = SU2_APIPYTHON_wrap
  SWIG_SO_REAL = _pysu2.${SO_EXT}
endif

if BUILD_REVERSE
  SWIG_SRC = ${abs_top_srcdir}/SU2_PY/pySU2/pySU2ad.i
  SWIG_WRAP_REAL = SU2_APIPYTHON_wrap_AD
  SWIG_SO_REAL = _pysu2ad.${SO_EXT}
endif

PYTHON_SITE_PACKAGES=$(shell python -c "import site; print(site.getsitepackages()[0])")
#NUMPY_INCLUDE = /usr/lib/python2.7/dist-packages/numpy/core/include
MPI4PY_INCLUDE = ${HOME}/.local/lib/python2.7/site-packages/mpi4py/include \
                 -I${PYTHON_SITE_PACKAGES}/mpi4py/include \
                 -I/Library/Python/2.7/site-packages/mpi4py/include

SWIG = swig

pySU2_INCLUDE = -I${abs_top_builddir}/Common/include \
	-I${abs_top_builddir}/SU2_CFD/include

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

if BUILD_NORMAL
  EXTRA_CC_FLAGS = -fPIC -O3 -std=c++11
endif

if BUILD_REVERSE
  EXTRA_CC_FLAGS = -fPIC -O3 -std=c++11 -DCODI_REVERSE_TYPE -I$(top_srcdir)/externals/codi/include -I$(top_srcdir)/externals/medi/include
#  EXTRA_CC_FLAGS = -fPIC -g3 -std=c++11 -DCODI_REVERSE_TYPE -I$(top_srcdir)/externals/codi/include -I$(top_srcdir)/externals/medi/include
endif

#SU2_BASE lib
if BUILD_NORMAL
  SU2SOLVER_LIB = ${abs_top_builddir}/SU2_CFD/obj/libSU2Core.a \
    ${abs_top_builddir}/Common/lib/libSU2.a
  PYLIB = pysu2.py
  SOLIB = _pysu2.so
endif

#SU2_AD lib
if BUILD_REVERSE
  SU2SOLVER_LIB = ${abs_top_builddir}/SU2_CFD/obj/libSU2Core_AD.a \
    ${abs_top_builddir}/Common/lib/libSU2_AD.a
  PYLIB = pysu2ad.py
  SOLIB = _pysu2ad.so
endif

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
	${INSTALL} ${PYLIB} ${bindir}
	${INSTALL} ${SOLIB} ${bindir}

