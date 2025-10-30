/*!
 * \file species_diffusion.hpp
 * \brief Declarations of numerics classes for discretization of
 *        viscous fluxes in species problems.
 * \author T. Kattmann
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include "../scalar/scalar_diffusion.hpp"
#include "../CNumerics.hpp"

/*!
 * \class CAvgGrad_Species
 * \brief Class for computing viscous term using average of gradients (species transport model).
 * \ingroup ViscDiscr
 */
template <class FlowIndices>
class CAvgGrad_Species final : public CAvgGrad_Scalar<FlowIndices> {
 private:
  using Base = CAvgGrad_Scalar<FlowIndices>;
  using Base::nVar;
  using Base::Eddy_Viscosity_i;
  using Base::Eddy_Viscosity_j;
  using Base::Diffusion_Coeff_i;
  using Base::Diffusion_Coeff_j;
  using Base::Density_i;
  using Base::Density_j;
  using Base::ScalarVar_i;
  using Base::ScalarVar_j;
  using Base::Proj_Mean_GradScalarVar;
  using Base::proj_vector_ij;
  using Base::Flux;
  using Base::Jacobian_i;
  using Base::Jacobian_j;
  using Base::Gamma;
  using Base::Gamma_Minus_One;
  using Base::Gas_Constant;
  using Base::Prandtl_Lam;
  using Base::Prandtl_Turb;
  using Base::Laminar_Viscosity_i;  
  using Base::Laminar_Viscosity_j;
  using Base::nDim;
  using Base::V_i;       // primitive vars at i
  using Base::V_j;
  using Base::PrimVar_Grad_i;
  using Base::PrimVar_Grad_j;
  using Base::Coord_i;
  using Base::Coord_j;
  using Base::Normal;
  using Base::correct_gradient;
  using Base::Thermal_Conductivity_i;
  using Base::Thermal_Conductivity_j;
  using Base::Cp_i;
  using Base::Cp_j;


  

  const bool turbulence;

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override {
    AD::SetPreaccIn(Diffusion_Coeff_i, nVar);
    AD::SetPreaccIn(Diffusion_Coeff_j, nVar);
    AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
    AD::SetPreaccIn(PrimVar_Grad_j, nDim+1, nDim);
  }

void CorrectGradient(su2double* GradPrimVar,
                                    const su2double val_PrimVar_i,
                                    const su2double val_PrimVar_j,
                                    const su2double* val_edge_vector,
                                    const su2double val_dist_ij_2,
                                    const unsigned short val_nPrimVar) {
    su2double Proj_Mean_GradPrimVar_Edge = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Proj_Mean_GradPrimVar_Edge += GradPrimVar[iDim]*val_edge_vector[iDim];
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      GradPrimVar[iDim] -= (Proj_Mean_GradPrimVar_Edge -
                                 (val_PrimVar_j-val_PrimVar_i))*val_edge_vector[iDim] / val_dist_ij_2;
    }
}



  /*!
   * \brief Species transport specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override {
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      
      switch(iVar)   {
      case 0:  {
      const su2double Diffusivity_Lam = 0.5 * (Density_i * Diffusion_Coeff_i[iVar] + Density_j * Diffusion_Coeff_j[iVar]);

      su2double Diffusivity_Turb = 0.0;

      if (turbulence) {
        const su2double Sc_t = config->GetSchmidt_Number_Turbulent();
        Diffusivity_Turb = 0.5 * (Eddy_Viscosity_i / Sc_t + Eddy_Viscosity_j / Sc_t);
      }

      const su2double Diffusivity = Diffusivity_Lam + Diffusivity_Turb;

      Flux[iVar] = Diffusivity * Proj_Mean_GradScalarVar[iVar];

      /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

      /*--- Off-diagonal entries are all zero. ---*/
      const su2double proj_on_rhoi = proj_vector_ij / Density_i;
      Jacobian_i[iVar][iVar] = -Diffusivity * proj_on_rhoi;

      const su2double proj_on_rhoj = proj_vector_ij / Density_j;
      Jacobian_j[iVar][iVar] = Diffusivity * proj_on_rhoj;
       break;
       }

       case 1: {

         const su2double Conductivity_lam = 0.5 * (Laminar_Viscosity_i + Laminar_Viscosity_j)/Prandtl_Lam;
         su2double Conductivity_Turb = 0.0;

         if (turbulence) {
           Conductivity_Turb = 0.5 * (Eddy_Viscosity_i + Eddy_Viscosity_j) / Prandtl_Turb;
         }
         su2double Conductivity = Conductivity_lam + Conductivity_Turb;

         Jacobian_i[iVar][iVar] = -Conductivity * proj_vector_ij;
         Jacobian_j[iVar][iVar] =  Conductivity * proj_vector_ij;

        Flux[iVar] = Conductivity * Proj_Mean_GradScalarVar[iVar];
        break;
       }
     // iVar
     }
  }
}

 public:



  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Whether to correct gradient for skewness.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Species(unsigned short val_nDim, unsigned short val_nVar, bool correct_grad, const CConfig* config)
    : CAvgGrad_Scalar<FlowIndices>(val_nDim, val_nVar, correct_grad, config),
      turbulence(config->GetKind_Turb_Model() != TURB_MODEL::NONE) {}
};
