/*!
 * \file CFEM_DG_SolverBase.cpp
 * \brief Main subroutines for the base class of the finite element flow solvers.
 * \author J. Alonso, E. van der Weide, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CFEM_DG_SolverBase.hpp"

CFEM_DG_SolverBase::CFEM_DG_SolverBase(void)
  : CSolver() {
}

CFEM_DG_SolverBase::CFEM_DG_SolverBase(CGeometry      *geometry,
                                       CConfig        *config,
                                       unsigned short iMesh)
  : CSolver() {

  /*--- Store the multigrid level and retrieve the number of dimensions
        and the number of markers. ---*/
  MGLevel = iMesh;
  nDim    = geometry->GetnDim();
  nMarker = config->GetnMarker_All();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);
  if( !DGGeometry) SU2_MPI::Error(string("Dynamic cast failed"), CURRENT_FUNCTION);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  ownedElemAdjLowTimeLevel = DGGeometry->GetOwnedElemAdjLowTimeLevel();
  haloElemAdjLowTimeLevel  = DGGeometry->GetHaloElemAdjLowTimeLevel();

  nMeshPoints = DGGeometry->GetNMeshPoints();
  meshPoints  = DGGeometry->GetMeshPoints();

  nMatchingInternalFacesWithHaloElem = DGGeometry->GetNMatchingFacesWithHaloElem();
  nMatchingInternalFacesLocalElem    = DGGeometry->GetNMatchingFacesInternal();
  matchingInternalFaces              = DGGeometry->GetMatchingFaces();

  boundaries = DGGeometry->GetBoundaries();

  timeCoefADER_DG                        = DGGeometry->GetTimeCoefADER_DG();
  timeInterpolDOFToIntegrationADER_DG    = DGGeometry->GetTimeInterpolDOFToIntegrationADER_DG();
  timeInterpolAdjDOFToIntegrationADER_DG = DGGeometry->GetTimeInterpolAdjDOFToIntegrationADER_DG();

  /*--- Determine the maximum number of integration points and DOFs
        used on this rank. ---*/
  const unsigned short nIntegrationMax = DGGeometry->DetermineMaxNIntegration();
  const unsigned short nDOFsMax        = DGGeometry->DetermineMaxNDOFs();

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

CFEM_DG_SolverBase::~CFEM_DG_SolverBase(void) {

  delete FluidModel;
}
