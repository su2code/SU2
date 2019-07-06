/*!
 * \file CTurbSAVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../include/variables/CTurbSAVariable.hpp"


CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim,
                                 unsigned short val_nvar, CConfig *config) : CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of S-A variables ---*/
  Solution[0] = val_nu_tilde;    Solution_Old[0] = val_nu_tilde;

  /*--- Initialization of the eddy viscosity ---*/
  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_nu_tilde;
    Solution_time_n1[0] = val_nu_tilde;
  }

  DES_LengthScale = 0.0;

}

void CTurbSAVariable::SetVortex_Tilting(su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity){

  su2double Strain[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}}, Omega, StrainDotVort[3], numVecVort[3];
  su2double numerator, trace0, trace1, denominator;

  AD::StartPreacc();
  AD::SetPreaccIn(PrimGrad_Flow, nDim+1, nDim);
  AD::SetPreaccIn(Vorticity, 3);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(muT);
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

  Vortex_Tilting = (numerator/denominator) * max(1.0,0.2*LaminarViscosity/muT);

  AD::SetPreaccOut(Vortex_Tilting);
  AD::EndPreacc();
}

CTurbSAVariable::~CTurbSAVariable(void) {}
