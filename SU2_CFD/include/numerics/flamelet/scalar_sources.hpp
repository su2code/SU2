/*!
 * \file scalar_sources.hpp
 * \brief Delarations of numerics classes for integration of source
 *        terms in scalar problems.
 * \author T. Economon, N. Beishuizen
 * \version 7.1.1 "Blackbird"
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

#pragma once

#include "../CNumerics.hpp"
#include "../../solvers/CScalarSolver.hpp"

/*!
 * \class CSourcePieceWise_transportedScalar_general
 * \brief Class for integrating the source terms of the transported scalar turbulence model equations.
 * \ingroup SourceDiscr
 * \author N. Beishuizen
 */
class CSourcePieceWise_transportedScalar_general final : public CNumerics {

private:
  su2double *Residual = nullptr;
  su2double **Jacobian_i = nullptr;
  bool flame;
  su2double *scalar_sources = nullptr;

  su2double source_pv;

  su2double Sc_t; 

  bool incompressible;
  bool viscous;
  bool axisymmetric;
  bool implicit;
  bool inc_rans;


  /*!
   * \brief Add contribution due to axisymmetric formulation to 2D residual
   */
  inline void ResidualAxisymmetric(void){
    su2double yinv,Density_i,Velocity_i[3];

    if (Coord_i[1] > EPS) {

      AD::SetPreaccIn(Coord_i[1]);

      yinv = 1.0/Coord_i[1];

      /*--- the incompressible density. Note that this is different for compressible flows ---*/

      Density_i = V_i[nDim+2];

      /*--- Set primitive variables at points iPoint. ---*/
    
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Velocity_i[iDim] = V_i[iDim+1];

      /*--- Inviscid component of the source term. ---*/
    

      for (auto iVar=0u; iVar < nVar; iVar++)
        Residual[iVar] -= yinv*Volume*Density_i*ScalarVar_i[iVar]*Velocity_i[1]; 

      if (implicit) {
        for (auto iVar=0u; iVar < nVar; iVar++) {
          for (auto jVar=0u; jVar < nVar; jVar++) {
            if (iVar == jVar) Jacobian_i[iVar][jVar] = -Velocity_i[1];
            Jacobian_i[iVar][jVar] *= yinv*Volume*Density_i; 
          }
        }
      }

      /*--- Add the viscous terms if necessary. ---*/
    
      if (viscous) {

        Laminar_Viscosity_i    = V_i[nDim+4];
        Eddy_Viscosity_i       = V_i[nDim+5];
        Thermal_Conductivity_i = V_i[nDim+6];

        su2double Mass_Diffusivity_Tur = 0.0;
        if (inc_rans) 
          Mass_Diffusivity_Tur = Eddy_Viscosity_i/Sc_t;

      // In case of flamelet, diffusion_coeff is actually rho*D
      // remember (rho*cp*D) = k when Le=1
      // and in case of transported scalar, diffusion_coeff is binary diffusion coefficient
        for (auto iVar=0u; iVar < nVar; iVar++){
          if (flame) 
            Residual[iVar] += yinv*Volume*(Diffusion_Coeff_i[iVar]+Mass_Diffusivity_Tur)*ScalarVar_Grad_i[iVar][1];
          else
            Residual[iVar] += yinv*Volume*Density_i*(Diffusion_Coeff_i[iVar]+Mass_Diffusivity_Tur)*ScalarVar_Grad_i[iVar][1];
        } 
      }

    } else {
    
      for (auto iVar=0u; iVar < nVar; iVar++)
        Residual[iVar] = 0.0;
  
      if (implicit) {
        for (auto iVar=0u; iVar < nVar; iVar++) {
          for (auto jVar=0u; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }
    
    }

  }

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_transportedScalar_general(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;

  /*!
   * \brief Set the value of the progress variable source term for the flamelet model.
   * \param[in] val_sourcepv_i - Value of the source term at point i.
   * \param[in] val_sourcepv_j - Value of the source term at point j.
   */
  inline void SetScalarSources(su2double *val_scalar_sources) override {
    for (auto i_var=0u; i_var < nVar; i_var++)
      scalar_sources[i_var] = val_scalar_sources[i_var];
  }

  inline void SetSourcePV(su2double val_sourcepv) override {
    source_pv = val_sourcepv;
  }

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_transportedScalar_general(void) override;


};


