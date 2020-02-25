/*!
 * \file CFEM_Knowles_NearInc.cpp
 * \brief FE numerics for nearly incompressible Knowles material model.
 * \author R. Sanchez
 * \version 7.0.1 "Blackbird"
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

#include "../../../include/numerics/elasticity/CFEM_Knowles_NearInc.hpp"


CFEM_Knowles_NearInc::CFEM_Knowles_NearInc(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CFEANonlinearElasticity(val_nDim, val_nVar, config) {

  /* -- The formulation adopted for this material model has been described by:
   * --
   * -- Suchocki, C., A Finite Element Implementation of Knowles stored-energy function:
   * -- theory, coding and applications, Archive of Mechanical Engineering, Vol. 58, pp. 319-346 (2011).
   * --
   * -- DOI: 10.2478/v10180-011-0021-7
   */

  Bk = config->GetKnowles_B();
  Nk = config->GetKnowles_N();

  trbbar = 0.0;
  Ek     = 0.0;
  Pr     = 0.0;

}

void CFEM_Knowles_NearInc::Compute_Plane_Stress_Term(CElement *element, CConfig *config) {

  SU2_MPI::Error("This material model cannot (yet) be used for plane stress.",CURRENT_FUNCTION);

}

void CFEM_Knowles_NearInc::Compute_Constitutive_Matrix(CElement *element, CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  /*--- Computation of the tensor cijkl ---*/

  unsigned short iVar, jVar, kVar, lVar;

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      for (kVar = 0; kVar < 3; kVar++){
        for (lVar = 0; lVar < 3; lVar++){
          cijkl[iVar][jVar][kVar][lVar] =
            term1 * ((1.0/2.0)*( deltaij(iVar,kVar)*b_Mat_Iso[jVar][lVar]
                      +deltaij(jVar,lVar)*b_Mat_Iso[iVar][kVar]
                      +deltaij(iVar,lVar)*b_Mat_Iso[jVar][kVar]
                      +deltaij(jVar,kVar)*b_Mat_Iso[iVar][lVar])
                 +(2.0/3.0)*(trbbar*deltaij(iVar,jVar)*deltaij(kVar,lVar)
                      -b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                      -deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar]))
             +term2 * ( b_Mat_Iso[iVar][jVar]*b_Mat_Iso[kVar][lVar]
                - trbbar*(b_Mat_Iso[iVar][jVar]*deltaij(kVar,lVar)
                     +deltaij(iVar,jVar)*b_Mat_Iso[kVar][lVar])
                + pow(trbbar,2.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar))
             +Kappa * (2.0 * J_F - 1.0) * deltaij(iVar,jVar) * deltaij(kVar,lVar);

        }
      }
    }
  }

  /*--- Reorganizing the tensor into the matrix D ---*/

  Assign_cijkl_D_Mat();


}

void CFEM_Knowles_NearInc::Compute_Stress_Tensor(CElement *element, CConfig *config) {

  /* -- Suchocki (2011) (full reference in class constructor). ---*/

  unsigned short iVar, jVar;

  /*--- Compute the isochoric deformation gradient Fbar and left Cauchy-Green tensor bbar ---*/
  Compute_Isochoric_F_b();

  trbbar = (b_Mat_Iso[0][0] + b_Mat_Iso[1][1] + b_Mat_Iso[2][2]) / 3.0;
  term1 = (Mu / J_F) * pow((1 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-1.0));
  term2 = 2.0 * (Mu / J_F) * (Bk * (Nk - 1.0) / Nk) *
      pow((1.0 + (Bk / Nk) * (3.0 * trbbar - 3.0)), (Nk-2.0));

  Ek = Kappa * (2.0 * J_F - 1.0);
  Pr = Kappa * (J_F - 1.0);

  for (iVar = 0; iVar < 3; iVar++){
    for (jVar = 0; jVar < 3; jVar++){
      Stress_Tensor[iVar][jVar] = term1 * (b_Mat_Iso[iVar][jVar] - deltaij(iVar,jVar)*trbbar) +
                                  deltaij(iVar,jVar) * Pr;
    }
  }

}

