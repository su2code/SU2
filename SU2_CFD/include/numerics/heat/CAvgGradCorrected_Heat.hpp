/*!
 * \file CAvgGradCorrected_Heat.hpp
 * \brief Delaration of numerics class CAvgGradCorrected_Heat, the
 *        implementation is in the CAvgGradCorrected_Heat.cpp file.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "../CNumerics.hpp"

/*!
 * \class CAvgGradCorrected_Heat
 * \brief Class for computing viscous term using average of gradients with correction (heat equation).
 * \ingroup ViscDiscr
 * \author O. Burghardt.
 * \version 7.0.0 "Blackbird"
 */
class CAvgGradCorrected_Heat : public CNumerics {
private:
  su2double **Mean_GradHeatVar;
  su2double *Proj_Mean_GradHeatVar_Kappa, *Proj_Mean_GradHeatVar_Edge, *Proj_Mean_GradHeatVar_Corrected;
  su2double *Edge_Vector;
  bool implicit;
  su2double dist_ij_2, proj_vector_ij, Thermal_Diffusivity_Mean;
  unsigned short iVar, iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_Heat(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_Heat(void);

  /*!
   * \brief Compute the viscous heat residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config);
};
