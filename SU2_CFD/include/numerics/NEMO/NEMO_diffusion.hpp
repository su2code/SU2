/*!
 * \file NEMO_diffusion.hpp
 * \brief Declarations of numerics classes for viscous flux computation.
 * \author S.R. Copeland, W. Maier, C. Garbacz.
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "CNEMONumerics.hpp"

/*!
 * \class CAvgGrad_NEMO
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author S.R. Copeland, W. Maier, C. Garbacz
 * \version 8.0.0 "Harrier"
 */
class CAvgGrad_NEMO : public CNEMONumerics {
private:

  su2double *Mean_PrimVar,    /*!< \brief Mean primitive variables. */
  *Mean_U,
  **Mean_GU,
  *Mean_dTdU,
  *Mean_dTvedU,
  *Mean_dPdU,
  *Mean_Eve,
  *Mean_Cvve,
  *PrimVar_i, *PrimVar_j,       /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,           /*!< \brief Mean value of the gradient. */
  *Mean_Diffusion_Coeff,        /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity,       /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,          /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,    /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */
  *ProjFlux,                    /*!< \brief Projection of the viscous fluxes. */
  dist_ij;                      /*!< \brief Length of the edge and face. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_NEMO(unsigned short val_nDim,
                unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_NEMO(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;
};

/*!
 * \class CAvgGradCorrected_NEMO
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author C. Garbacz, W. Maier, S.R. Copeland.
 * \version 8.0.0 "Harrier"
 */
class CAvgGradCorrected_NEMO : public CNEMONumerics {
private:
  unsigned short
  nPrimVar, nPrimVarGrad;       /*!< \brief Iterators in dimension an variable. */
  su2double
  *Mean_PrimVar,                /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,       /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,           /*!< \brief Mean value of the gradient. */
  *Mean_Eve,                    /*!< \brief Mean value of eve. */
  *Mean_Cvve,                   /*!< \brief Mean value of cvve. */
  Edge_Vector[MAXNDIM]={0.0},   /*!< \brief Vector from point i to point j. */
  *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Inner product of the Mean gradient and the edge vector. */
  *Mean_Diffusion_Coeff,        /*!< \brief Mean value of the species diffusion coefficient. */
  Mean_Laminar_Viscosity,       /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,          /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,    /*!< \brief Mean value of the thermal conductivity. */
  Mean_Thermal_Conductivity_ve, /*!< \brief Mean value of the vib-el. thermal conductivity. */
  *ProjFlux,                    /*!< \brief Projection of the viscous fluxes. */
  dist_ij;                      /*!< \brief Length of the edge and face. */
  bool implicit;                /*!< \brief Implicit calculus. */

  su2double* Flux = nullptr;    /*!< \brief The flux / residual across the edge. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive variables of the problem.
   * \param[in] val_nPrimVarGrad - Number of variables in the primitive variable gradient.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_NEMO(unsigned short val_nDim,
                unsigned short val_nVar,
                unsigned short val_nPrimVar,
                unsigned short val_nPrimVarGrad,
                CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_NEMO(void);

    /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) final;

};
