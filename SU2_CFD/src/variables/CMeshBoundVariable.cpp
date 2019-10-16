/*!
 * \file CMeshBoundVariable.cpp
 * \brief Definition of the boundary variables for mesh motion using a pseudo-elastic approach.
 * \author Ruben Sanchez
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/variables/CMeshBoundVariable.hpp"

CMeshBoundVariable::CMeshBoundVariable(unsigned long npoint, unsigned long ndim, CConfig *config) :
  CMeshVariable(npoint, ndim, config) {

  VertexMap.Reset(nPoint);
}

void CMeshBoundVariable::AllocateBoundaryVariables(CConfig *config) {

  if (VertexMap.GetIsValid()) return; // nothing to do

  /*--- Count number of vertices and build map ---*/

  unsigned long nBoundPt = VertexMap.Build();

  /*--- Allocate ---*/

  Boundary_Displacement.resize(nBoundPt,nDim) = su2double(0.0);
}

void CMeshBoundVariable::Register_BoundDisp(bool input) {
  if (input) {  
    for (unsigned long iVertex = 0; iVertex < Boundary_Displacement.rows(); iVertex++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterInput(Boundary_Displacement(iVertex,iVar));
  }
  else {
    for (unsigned long iVertex = 0; iVertex < Boundary_Displacement.rows(); iVertex++)
      for (unsigned long iVar = 0; iVar < nVar; iVar++)
        AD::RegisterOutput(Boundary_Displacement(iVertex,iVar));
  }
}
