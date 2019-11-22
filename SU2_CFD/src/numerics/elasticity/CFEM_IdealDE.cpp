/*!
 * \file CFEM_IdealDE.cpp
 * \brief Definition of ideal dielectric elastomer.
 * \author R. Sanchez
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/elasticity/CFEM_IdealDE.hpp"


CFEM_IdealDE::CFEM_IdealDE(unsigned short val_nDim, unsigned short val_nVar,
                           CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {

  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Zhao, X. and Suo, Z., Method to analyze programmable deformation of
   * -- dielectric elastomer layers, Applied Physics Letters 93, 251902 (2008).
   * --
   * -- http://imechanica.org/node/4234
   */

  trbbar = 0.0;
  Eg     = 0.0;
  Eg23   = 0.0;
  Ek     = 0.0;
  Pr     = 0.0;

}

void CFEM_IdealDE::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.", CURRENT_FUNCTION);

}

void CFEM_IdealDE::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  if (nDim == 2){

    D_Mat[0][0] = Eg23*(b_Mat_Iso[0][0]+trbbar)+Ek;
    D_Mat[1][1] = Eg23*(b_Mat_Iso[1][1]+trbbar)+Ek;

    D_Mat[0][1] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;
    D_Mat[1][0] = -Eg23*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1]-trbbar)+Ek;

    D_Mat[0][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][0] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[1][2] = Eg23*b_Mat_Iso[0][1]/2.0;
    D_Mat[2][1] = Eg23*b_Mat_Iso[0][1]/2.0;

    D_Mat[2][2] = Eg*(b_Mat_Iso[0][0]+b_Mat_Iso[1][1])/2.0;

  }
  else {
    SU2_MPI::Error("This material model cannot be used for 3D yet.", CURRENT_FUNCTION);
  }

}

void CFEM_IdealDE::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  /* -- Zhao, X. and Suo, Z. (2008) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;
  su2double dij = 0.0;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  cout.precision(15);
  // Stress terms

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  Eg = Mu / J_F;
  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);
  Eg23 = 2.0 * Eg / 3.0;

  // Stress tensor

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      if (iVar == jVar) dij = 1.0;
      else if (iVar != jVar) dij = 0.0;
      Stress_Tensor[iVar][jVar] = Eg * ( b_Mat_Iso[iVar][jVar] - dij * trbbar) + dij * Pr ;
    }
  }

}

