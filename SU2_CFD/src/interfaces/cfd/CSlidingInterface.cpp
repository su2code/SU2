/*!
 * \file CSlidingInterface.cpp
 * \brief Declaration and inlines of the class to transfer conservative variables
 *        from a generic zone into another
 * \author G. Gori Politecnico di Milano
 * \version 8.0.1 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/interfaces/cfd/CSlidingInterface.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"																	

CSlidingInterface::CSlidingInterface(unsigned short val_nVar, unsigned short val_nConst) : CInterface() {

  Physical_Constants = new su2double[val_nConst] ();
  Donor_Variable     = new su2double[val_nVar] ();
  Target_Variable    = new su2double[val_nVar+1] ();

  valAggregated      = false;

  nVar = val_nVar;

}

void CSlidingInterface::GetDonor_Velocity_RotatingFrame(const CConfig *donor_config, CGeometry *donor_geometry, 
                                 const CConfig *target_config, CGeometry *target_geometry){
	
	/*---rotate the velocity for rotating frame---*/
	if (donor_config->GetRotating_Frame()==YES){
		
		unsigned short nDim = donor_geometry->GetnDim();
		
		su2double Theta, Phi, Psi;
		su2double oriVel[3] = {0.0, 0.0, 0.0};	
		su2double rotVel[3] = {0.0, 0.0, 0.0};	
		su2double rotMatrix[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};	
		su2double Omega[3] = {0.0, 0.0, 0.0};
		const su2double zeros[3] = {0.0};
		su2double dt = donor_config->GetDelta_UnstTimeND();
		unsigned long TimeIter = donor_config->GetTimeIter();
		for (unsigned short iDim=0; iDim<3; iDim++){
		  Omega[iDim] = donor_config->GetRotation_Rate(iDim)/donor_config->GetOmega_Ref();
		}
		
		/*--- Compute the rotation matrix. Note that the implicit
		 ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
		Theta    = Omega[0]*dt*TimeIter;      Phi = Omega[1]*dt*TimeIter;     Psi = Omega[2]*dt*TimeIter;
		GeometryToolbox::RotationMatrix(Theta, Phi, Psi, rotMatrix);
		
		/*--- Velocities before rotating ---*/
		for (unsigned short iDim = 0; iDim < nDim; ++iDim)
			oriVel[iDim] = Donor_Variable[iDim+1];
		
		/*--- Compute transformed velocities. ---*/
		GeometryToolbox::Rotate(rotMatrix, zeros, oriVel, rotVel);
		
		/*--- set the rotated velocity ---*/
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		  Donor_Variable[iDim+1] = rotVel[iDim];
		}
		
	}
	
	if (target_config->GetRotating_Frame()==YES){
		
		unsigned short nDim = target_geometry->GetnDim();
		
		su2double Theta, Phi, Psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
		su2double oriVel[3] = {0.0, 0.0, 0.0};	
		su2double rotVel[3] = {0.0, 0.0, 0.0};	
		su2double rotMatrix[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};	
		su2double Omega[3] = {0.0, 0.0, 0.0};
		const su2double zeros[3] = {0.0};
		su2double dt = target_config->GetDelta_UnstTimeND();
		unsigned long TimeIter = donor_config->GetTimeIter();
		for (unsigned short iDim=0; iDim<3; iDim++){
		  Omega[iDim] = -target_config->GetRotation_Rate(iDim)/target_config->GetOmega_Ref();
		}

		/*--- Compute the rotation matrix. Note that the implicit
		 ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
		Theta    = Omega[0]*dt*TimeIter;      Phi = Omega[1]*dt*TimeIter;     Psi = Omega[2]*dt*TimeIter;
		GeometryToolbox::RotationMatrix(Theta, Phi, Psi, rotMatrix);
		
		/*--- velocities before rotating ---*/
		for (unsigned short iDim = 0; iDim < nDim; ++iDim)
			oriVel[iDim] = Donor_Variable[iDim+1];
		
		/*--- Compute transformed velocities. ---*/
		GeometryToolbox::Rotate(rotMatrix, zeros, oriVel, rotVel);
		
		/*--- set the rotated velocity ---*/
		for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		  Donor_Variable[iDim+1] = rotVel[iDim];
		}
		
	}
	
	

}
void CSlidingInterface::GetDonor_Variable(CSolver *donor_solution, CGeometry *donor_geometry,
                                          const CConfig *donor_config, unsigned long Marker_Donor,
                                          unsigned long Vertex_Donor, unsigned long Point_Donor) {

  const auto nDonorVar = donor_solution->GetnPrimVar();
  /// TODO: Replace with approach compatible with any number of variables (e.g. encapsulate in a "solver info" object).
  const bool turbulent = nDonorVar <= 2;
  if (turbulent){
    /*---  For turbulent solver retrieve solution variables and set then as the donor variables. ---*/
    for (unsigned short iVar = 0; iVar < nDonorVar; iVar++) {
      Donor_Variable[iVar] = donor_solution->GetNodes()->GetSolution(Point_Donor, iVar);
    }
  } else {
    /*---  For flow solver retrieve primitive variables and set them as the donor variables. ---*/
    for (unsigned short iVar = 0; iVar < nDonorVar; iVar++) {
      Donor_Variable[iVar] = donor_solution->GetNodes()->GetPrimitive(Point_Donor, iVar);
    }
  }
}

void CSlidingInterface::InitializeTarget_Variable(CSolver *target_solution, unsigned long Marker_Target,
                                                  unsigned long Vertex_Target, unsigned short nDonorPoints) {

  target_solution->SetnSlidingStates(Marker_Target, Vertex_Target, nDonorPoints); // This is to allocate
  target_solution->SetSlidingStateStructure(Marker_Target, Vertex_Target);
  target_solution->SetnSlidingStates(Marker_Target, Vertex_Target, 0); // Reset counter to 0

}

void CSlidingInterface::SetTarget_Variable(CSolver *target_solution, CGeometry *target_geometry,
                                           const CConfig *target_config, unsigned long Marker_Target,
                                           unsigned long Vertex_Target, unsigned long Point_Target) {

  unsigned short iVar, iDonorVertex, nTargetVar;
  nTargetVar = target_solution->GetnPrimVar();
  /*--- Set the Sliding solution with the value of the Target Variable ---*/

  iDonorVertex = target_solution->GetnSlidingStates(Marker_Target, Vertex_Target);

  for (iVar = 0; iVar < nTargetVar+1; iVar++)
    target_solution->SetSlidingState(Marker_Target, Vertex_Target, iVar, iDonorVertex, Target_Variable[iVar]);

  target_solution->SetnSlidingStates( Marker_Target, Vertex_Target, iDonorVertex + 1 );
}
