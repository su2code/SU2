/*!
 * \file CTurbVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 7.0.3 "Blackbird"
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


#include "../../include/variables/CTurbVariable.hpp"


CTurbVariable::CTurbVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CVariable(npoint, ndim, nvar, config), Gradient_Reconstruction(config->GetReconstructionGradientRequired() ? Gradient_Aux : Gradient) {

  /*--- Allocate space for the harmonic balance source terms ---*/

  if (config->GetTime_Marching() == HARMONIC_BALANCE) {
    HB_Source.resize(nPoint,nVar) = su2double(0.0);
  }

  /*--- Gradient related fields ---*/

  Gradient.resize(nPoint,nVar,nDim,0.0);

  if (config->GetReconstructionGradientRequired()) {
    Gradient_Aux.resize(nPoint,nVar,nDim,0.0);
  }

  if (config->GetLeastSquaresRequired()) {
    Rmatrix.resize(nPoint,nDim,nDim,0.0);
  }

  /*--- Always allocate the slope limiter, and the auxiliar
   variables (check the logic - JST with 2nd order Turb model) ---*/

  Limiter.resize(nPoint,nVar) = su2double(0.0);
  Solution_Max.resize(nPoint,nVar) = su2double(0.0);
  Solution_Min.resize(nPoint,nVar) = su2double(0.0);

  Delta_Time.resize(nPoint) = su2double(0.0);

  /* Under-relaxation parameter. */
  UnderRelaxation.resize(nPoint) = su2double(1.0);
  LocalCFL.resize(nPoint) = su2double(0.0);

}

void CTurbVariable::SetVortex_Tilting(unsigned long iPoint, const su2double* const* PrimGrad_Flow,
                                        const su2double* Vorticity, su2double LaminarViscosity) {

  su2double Strain[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}}, Omega, StrainDotVort[3], numVecVort[3];
  su2double numerator, trace0, trace1, denominator;

  AD::StartPreacc();
  AD::SetPreaccIn(PrimGrad_Flow, nDim+1, nDim);
  AD::SetPreaccIn(Vorticity, 3);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(muT(iPoint));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(LaminarViscosity);

  Strain[0][0] = PrimGrad_Flow[1][0];
  Strain[1][0] = 0.5*(PrimGrad_Flow[2][0] + PrimGrad_Flow[1][1]);
  Strain[0][1] = 0.5*(PrimGrad_Flow[1][1] + PrimGrad_Flow[2][0]);
  Strain[1][1] = PrimGrad_Flow[2][1];
  if (nDim == 3){
    Strain[0][2] = 0.5*(PrimGrad_Flow[3][0] + PrimGrad_Flow[1][2]);
    Strain[1][2] = 0.5*(PrimGrad_Flow[3][1] + PrimGrad_Flow[2][2]);
    Strain[2][0] = 0.5*(PrimGrad_Flow[1][2] + PrimGrad_Flow[3][0]);
    Strain[2][1] = 0.5*(PrimGrad_Flow[2][2] + PrimGrad_Flow[3][1]);
    Strain[2][2] = PrimGrad_Flow[3][2];
  }

  Omega = sqrt(Vorticity[0]*Vorticity[0] + Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);

  StrainDotVort[0] = Strain[0][0]*Vorticity[0]+Strain[0][1]*Vorticity[1]+Strain[0][2]*Vorticity[2];
  StrainDotVort[1] = Strain[1][0]*Vorticity[0]+Strain[1][1]*Vorticity[1]+Strain[1][2]*Vorticity[2];
  StrainDotVort[2] = Strain[2][0]*Vorticity[0]+Strain[2][1]*Vorticity[1]+Strain[2][2]*Vorticity[2];

  numVecVort[0] = StrainDotVort[1]*Vorticity[2] - StrainDotVort[2]*Vorticity[1];
  numVecVort[1] = StrainDotVort[2]*Vorticity[0] - StrainDotVort[0]*Vorticity[2];
  numVecVort[2] = StrainDotVort[0]*Vorticity[1] - StrainDotVort[1]*Vorticity[0];

  numerator = sqrt(6.0) * sqrt(numVecVort[0]*numVecVort[0] + numVecVort[1]*numVecVort[1] + numVecVort[2]*numVecVort[2]);
  trace0 = 3.0*(pow(Strain[0][0],2.0) + pow(Strain[1][1],2.0) + pow(Strain[2][2],2.0));
  trace1 = pow(Strain[0][0] + Strain[1][1] + Strain[2][2],2.0);
  denominator = pow(Omega, 2.0) * sqrt(trace0-trace1);

  Vortex_Tilting(iPoint) = (numerator/denominator) * max(1.0,0.2*LaminarViscosity/muT(iPoint));
  
  AD::SetPreaccOut(Vortex_Tilting(iPoint));
  AD::EndPreacc();
}
