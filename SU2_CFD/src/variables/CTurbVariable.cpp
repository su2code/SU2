/*!
 * \file CTurbVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
  : CScalarVariable(npoint, ndim, nvar, config) {

    turb_index.resize(nPoint) = su2double(1.0);
    intermittency.resize(nPoint) = su2double(1.0);

    DES_LengthScale.resize(nPoint) = su2double(0.0);
    
    Vortex_Tilting.resize(nPoint);

   }


void CTurbVariable::SetVortex_Tilting(unsigned long iPoint, CMatrixView<const su2double> PrimGrad_Flow,
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
  /*--- 3 tr(S^2) - tr(S)^2, where tr(S^2) is the sum of the squares of all the components of S. ---*/
  trace0 = 0.0;
  for (auto iDim = 0u; iDim < 3; iDim++)
    for (auto jDim = 0u; jDim < 3; jDim++) trace0 += 3.0 * pow(Strain[iDim][jDim], 2);
  trace1 = pow(Strain[0][0] + Strain[1][1] + Strain[2][2],2.0);
  denominator = pow(Omega, 2.0) * sqrt(max(trace0-trace1, 0.0));

  /*--- Without vorticity or strain anisotropy the measure is undefined, make it neutral (F_KH = 1). ---*/
  if (denominator < 1e-20) {
    Vortex_Tilting(iPoint) = 1.0;
  } else {
    Vortex_Tilting(iPoint) = (numerator/denominator) * max(1.0,0.2*LaminarViscosity/max(muT(iPoint), 1e-20));
  }

  AD::SetPreaccOut(Vortex_Tilting(iPoint));
  AD::EndPreacc();
}

void CTurbVariable::RegisterEddyViscosity(bool input) {
  RegisterContainer(input, muT);
}