/*!
 * \file numerics_viscous_ad.cpp
 * \brief This file contains all the viscous term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.4
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"

void CAvgGrad_AdjDiscTurbSA::SetDirectResidual_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CAvgGrad_TurbSA__SetResidual

    int ii2;
    int ii1;
    Density_id = U_id[0];
    Density_i = U_i[0];
    Density_jd = U_jd[0];
    Density_j = U_j[0];
    /*--- Compute mean effective viscosity ---*/
    nu_id = (Laminar_Viscosity_id*Density_i-Laminar_Viscosity_i*Density_id)/(
        Density_i*Density_i);
    nu_i = Laminar_Viscosity_i/Density_i;
    nu_jd = (Laminar_Viscosity_jd*Density_j-Laminar_Viscosity_j*Density_jd)/(
        Density_j*Density_j);
    nu_j = Laminar_Viscosity_j/Density_j;
    nu_ed = 0.5*(nu_id+nu_jd+TurbVar_id[0]+TurbVar_jd[0]);
    nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
    /*--- Compute vector going from iPoint to jPoint ---*/
    dist_ij_2 = 0;
    proj_vector_ij = 0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Edge_Vectord[iDim] = 0.0;
        Edge_Vector[iDim] = Coord_j[iDim] - Coord_i[iDim];
        dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
        proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
    }
    proj_vector_ij = proj_vector_ij/dist_ij_2;
    for (ii1 = 0; ii1 < nVar; ++ii1)
        for (ii2 = 0; ii2 < nDim; ++ii2)
            Mean_GradTurbVard[ii1][ii2] = 0.0;
    for (ii1 = 0; ii1 < nVar; ++ii1)
        Proj_Mean_GradTurbVar_Kappad[ii1] = 0.0;
    /*--- Mean gradient approximation ---*/
    // to normalize vectors
    for (iVar = 0; iVar < nVar; ++iVar) {
        Proj_Mean_GradTurbVar_Kappad[iVar] = 0.0;
        Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
        Proj_Mean_GradTurbVar_Edged[iVar] = 0.0;
        Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            Mean_GradTurbVard[iVar][iDim] = 0.5*(TurbVar_Grad_id[iVar][iDim]+
                TurbVar_Grad_jd[iVar][iDim]);
            Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim]+
                TurbVar_Grad_j[iVar][iDim]);
            Proj_Mean_GradTurbVar_Kappad[iVar] = Proj_Mean_GradTurbVar_Kappad[
                iVar] + Normal[iDim]*Mean_GradTurbVard[iVar][iDim];
            Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*
            Normal[iDim];
        }
    }
    val_residuald[0] = (nu_ed*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*
        Proj_Mean_GradTurbVar_Kappad[0])/sigma;
    val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;


//SU2_DIFF END CAvgGrad_TurbSA__SetResidual

}
