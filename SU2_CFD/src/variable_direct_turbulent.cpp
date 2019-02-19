/*!
 * \file variable_direct_turbulent.cpp
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

#include "../include/variable_structure.hpp"

CTurbVariable::CTurbVariable(void) : CVariable() {
  
  /*--- Array initialization ---*/
  HB_Source = NULL;
  
}

CTurbVariable::CTurbVariable(unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CVariable(val_nDim, val_nvar, config) {
  
  unsigned short iVar;

  /*--- Array initialization ---*/
  
  HB_Source = NULL;
  
  /*--- Allocate space for the harmonic balance source terms ---*/
  
  if (config->GetUnsteady_Simulation() == HARMONIC_BALANCE) {
    HB_Source = new su2double[nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      HB_Source[iVar] = 0.0;
  }
  
  /*--- Always allocate the slope limiter,
   and the auxiliar variables (check the logic - JST with 2nd order Turb model - ) ---*/

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;
  
  Solution_Max = new su2double [nVar];
  Solution_Min = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }
  
}

CTurbVariable::~CTurbVariable(void) { }

su2double CTurbVariable::GetmuT() { return muT; }

void CTurbVariable::SetmuT(su2double val_muT) { muT = val_muT; }

CTurbSAVariable::CTurbSAVariable(void) : CTurbVariable() { }

CTurbSAVariable::CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {
  
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

CTurbSAVariable::~CTurbSAVariable(void) {
  
  if (HB_Source != NULL) delete [] HB_Source;
  
}

CTurbSSTVariable::CTurbSSTVariable(void) : CTurbVariable() { }

CTurbSSTVariable::CTurbSSTVariable(su2double val_kine, su2double val_omega, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar,
                                   su2double *constants, CConfig *config)
: CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- Initialization of variables ---*/
  
  Solution[0] = val_kine;     Solution_Old[0] = val_kine;
  Solution[1] = val_omega;  Solution_Old[1] = val_omega;
  
  sigma_om2 = constants[3];
  beta_star = constants[6];
  
  F1   = 1.0;
  F2   = 0.0;
  CDkw = 0.0;
  
  /*--- Initialization of eddy viscosity ---*/
  
  muT = val_muT;
  
  /*--- Allocate and initialize solution for the dual time strategy ---*/
  
  if (dual_time) {
    Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_omega;
    Solution_time_n1[0]  = val_kine; Solution_time_n1[1]  = val_omega;
  }
    
}

CTurbSSTVariable::~CTurbSSTVariable(void) {

  if (HB_Source != NULL) delete [] HB_Source;
  
}

void CTurbSSTVariable::SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density) {
  unsigned short iDim;
  su2double arg2, arg2A, arg2B, arg1;

  AD::StartPreacc();
  AD::SetPreaccIn(val_viscosity);  AD::SetPreaccIn(val_dist);
  AD::SetPreaccIn(val_density);
  AD::SetPreaccIn(Solution, nVar);
  AD::SetPreaccIn(Gradient, nVar, nDim);
  
  /*--- Cross diffusion ---*/
  
  CDkw = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    CDkw += Gradient[0][iDim]*Gradient[1][iDim];
  CDkw *= 2.0*val_density*sigma_om2/Solution[1];
  CDkw = max(CDkw, pow(10.0, -20.0));
  
  /*--- F1 ---*/
  
  arg2A = sqrt(Solution[0])/(beta_star*Solution[1]*val_dist+EPS*EPS);
  arg2B = 500.0*val_viscosity / (val_density*val_dist*val_dist*Solution[1]+EPS*EPS);
  arg2 = max(arg2A, arg2B);
  arg1 = min(arg2, 4.0*val_density*sigma_om2*Solution[0] / (CDkw*val_dist*val_dist+EPS*EPS));
  F1 = tanh(pow(arg1, 4.0));
  
  /*--- F2 ---*/
  
  arg2 = max(2.0*arg2A, arg2B);
  F2 = tanh(pow(arg2, 2.0));

  AD::SetPreaccOut(F1); AD::SetPreaccOut(F2); AD::SetPreaccOut(CDkw);
  AD::EndPreacc();
  
}
