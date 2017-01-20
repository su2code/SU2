/*
################################################################################
#
# \file pySU2.i
# \brief Configuration file for the Swig compilation of the Python wrapper.
# \author D. Thomas
# \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
*/

%feature("autodoc","1");

%module(docstring=
"'pysu2' module",
directors="1",
threads="1"
) pysu2
%{

#include "../../SU2_CFD/include/driver_structure.hpp"

%}

// ----------- USED MODULES ------------
%import "../../Common/include/datatypes/primitive_structure.hpp"
%import "../../Common/include/mpi_structure.hpp"
%include "std_string.i"
%include "typemaps.i"
//%include "numpy.i"
#ifdef HAVE_MPI			//Need mpi4py only for a parallel build of the wrapper.
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
#endif
// ----------- API CLASSES ----------------
%include "../../SU2_CFD/include/driver_structure.hpp"
