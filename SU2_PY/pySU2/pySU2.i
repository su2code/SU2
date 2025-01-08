/*
################################################################################
#
# \file pySU2.i
# \brief Configuration file for the Swig compilation of the Python wrapper.
# \author D. Thomas
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
#include "../../Common/include/containers/CPyWrapperMatrixView.hpp"
#include "../../SU2_CFD/include/drivers/CDiscAdjSinglezoneDriver.hpp"
#include "../../SU2_CFD/include/drivers/CDriver.hpp"
#include "../../SU2_CFD/include/drivers/CDriverBase.hpp"
#include "../../SU2_CFD/include/drivers/CMultizoneDriver.hpp"
#include "../../SU2_CFD/include/drivers/CSinglezoneDriver.hpp"
#include "../../SU2_DEF/include/drivers/CDeformationDriver.hpp"
%}

// ----------- USED MODULES ------------
%import "../../Common/include/code_config.hpp"
%import "../../Common/include/basic_types/datatype_structure.hpp"
%import "../../Common/include/parallelization/mpi_structure.hpp"

%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_map.i"
%include "typemaps.i"
//%include "numpy.i"
#ifdef HAVE_MPI                    //Need mpi4py only for a parallel build of the wrapper.
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
#endif

namespace std {
   %template() vector<bool>;
   %template() vector<unsigned short>;
   %template() vector<unsigned long>;
   %template() vector<double>;
   %template() vector<string>;
   %template() map<string, unsigned short>;
   %template() map<string, string>;
   %template() pair<unsigned long, unsigned long>;
}

// ----------- API CLASSES ----------------

//Constants definitions
/*!
 * \brief different software components of SU2
 */
enum class SU2_COMPONENT {
  SU2_CFD, /*!< \brief Running the SU2_CFD software. */
  SU2_DEF, /*!< \brief Running the SU2_DEF software. */
  SU2_DOT, /*!< \brief Running the SU2_DOT software. */
  SU2_GEO, /*!< \brief Running the SU2_GEO software. */
  SU2_SOL  /*!< \brief Running the SU2_SOL software. */
};

const unsigned int MESH_0 = 0; /*!< \brief Definition of the finest grid level. */
const unsigned int MESH_1 = 1; /*!< \brief Definition of the finest grid level. */
const unsigned int ZONE_0 = 0; /*!< \brief Definition of the first grid domain. */
const unsigned int ZONE_1 = 1; /*!< \brief Definition of the first grid domain. */

%include "../../Common/include/containers/CPyWrapperMatrixView.hpp"
%include "../../SU2_CFD/include/drivers/CDriverBase.hpp"
%include "../../SU2_CFD/include/drivers/CDriver.hpp"
%include "../../SU2_CFD/include/drivers/CSinglezoneDriver.hpp"
%include "../../SU2_CFD/include/drivers/CMultizoneDriver.hpp"
%include "../../SU2_CFD/include/drivers/CDiscAdjSinglezoneDriver.hpp"
%include "../../SU2_DEF/include/drivers/CDeformationDriver.hpp"
