/*!
 * \file CTurbVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
 * \version 8.4.0 "Harrier"
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
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"


CTurbVariable::CTurbVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config)
  : CScalarVariable(npoint, ndim, nvar, config) {

    turb_index.resize(nPoint) = su2double(1.0);
    intermittency.resize(nPoint) = su2double(1.0);

    DES_LengthScale.resize(nPoint) = su2double(0.0);
    
    Vortex_Tilting.resize(nPoint);

   }


void CTurbVariable::SetVortex_Tilting(unsigned long iPoint, su2double **Strain,
                                      const su2double* Vorticity, su2double LaminarViscosity) {

  su2double Omega, StrainDotVort[3], numVecVort[3];
  su2double numerator, trace0, trace1, denominator;

  AD::StartPreacc();
  AD::SetPreaccIn(Strain, nDim, nDim);
  AD::SetPreaccIn(Vorticity, 3);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(muT(iPoint));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(LaminarViscosity);

  Omega = GeometryToolbox::Norm(3, Vorticity);

  StrainDotVort[0] = GeometryToolbox::DotProduct(3, Strain[0], Vorticity);
  StrainDotVort[1] = GeometryToolbox::DotProduct(3, Strain[1], Vorticity);
  StrainDotVort[2] = GeometryToolbox::DotProduct(3, Strain[2], Vorticity);

  GeometryToolbox::CrossProduct(StrainDotVort, Vorticity, numVecVort);

  numerator = sqrt(6.0) * GeometryToolbox::Norm(3, numVecVort);
  trace0 = 3.0*(pow(Strain[0][0],2.0) + pow(Strain[1][1],2.0) + pow(Strain[2][2],2.0));
  trace1 = pow(Strain[0][0] + Strain[1][1] + Strain[2][2],2.0);
  denominator = pow(Omega, 2.0) * sqrt(trace0-trace1);

  Vortex_Tilting(iPoint) = (numerator/denominator) * max(1.0,0.2*LaminarViscosity/muT(iPoint));

  AD::SetPreaccOut(Vortex_Tilting(iPoint));
  AD::EndPreacc();
}

void CTurbVariable::RegisterEddyViscosity(bool input) {
  RegisterContainer(input, muT);
}