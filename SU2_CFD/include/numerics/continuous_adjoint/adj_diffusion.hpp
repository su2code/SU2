/*!
 * \file adj_diffusion.hpp
 * \brief Declarations of numerics classes for continuous adjoint
 *        diffusion discretization. Implemented in adj_diffusion.cpp.
 * \author F. Palacios, T. Economon
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

#include "../CNumerics.hpp"

/*!
 * \class CAvgGrad_AdjFlow
 * \brief Class for computing the adjoint viscous terms.
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CAvgGrad_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;  /*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;  /*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double *Mean_GradPsiE;  /*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;  /*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;    /*!< \brief Vector going from node i to node j. */
  bool implicit;             /*!< \brief Implicit calculus. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjFlow(void) override;

  /*!
   * \brief Residual computation.
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total residual at point j.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) override;
};

/*!
 * \class CAvgGradCorrected_AdjFlow
 * \brief Class for computing the adjoint viscous terms, including correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGradCorrected_AdjFlow : public CNumerics {
private:
  su2double *Velocity_i;  /*!< \brief Auxiliary vector for storing the velocity of point i. */
  su2double *Velocity_j;  /*!< \brief Auxiliary vector for storing the velocity of point j. */
  su2double *Mean_Velocity;
  su2double **Mean_GradPsiVar;           /*!< \brief Counter for dimensions of the problem. */
  su2double *Edge_Vector;                /*!< \brief Vector going from node i to node j. */
  su2double *Proj_Mean_GradPsiVar_Edge;  /*!< \brief Projection of Mean_GradPsiVar onto Edge_Vector. */
  su2double *Mean_GradPsiE;              /*!< \brief Counter for dimensions of the problem. */
  su2double **Mean_GradPhi;              /*!< \brief Counter for dimensions of the problem. */
  bool implicit;                         /*!< \brief Boolean controlling Jacobian calculations. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjFlow(void) override;

  /*!
   * \brief Compute the adjoint flow viscous residual in a non-conservative way using an average of gradients and derivative correction.
   * \param[out] val_residual_i - Pointer to the viscous residual at point i.
   * \param[out] val_residual_j - Pointer to the viscous residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) override;
};

/*!
 * \class CAvgGrad_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author F. Palacios
 */
class CAvgGrad_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_AdjTurb(void) override;

  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */

  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;

  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) override;
};

/*!
 * \class CAvgGradCorrected_AdjTurb
 * \brief Class for adjoint turbulent using average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGradCorrected_AdjTurb : public CNumerics {
private:
  su2double **Mean_GradTurbPsi;
  su2double *Proj_Mean_GradTurbPsi_Kappa, *Proj_Mean_GradTurbPsi_Edge, *Proj_Mean_GradTurbPsi_Corrected;
  su2double *Edge_Vector;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_AdjTurb(void) override;

  /*!
   * \brief Compute the adjoint turbulent residual using average of gradients and a derivative correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */

  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;

  /*!
   * \overload
   * \param[out] val_residual_i - Pointer to the total residual at point i.
   * \param[out] val_residual_j - Pointer to the total viscosity residual at point j.
   * \param[out] val_Jacobian_ii - Jacobian of the numerical method at node i (implicit computation) from node i.
   * \param[out] val_Jacobian_ij - Jacobian of the numerical method at node i (implicit computation) from node j.
   * \param[out] val_Jacobian_ji - Jacobian of the numerical method at node j (implicit computation) from node i.
   * \param[out] val_Jacobian_jj - Jacobian of the numerical method at node j (implicit computation) from node j.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual_i, su2double *val_residual_j,
                       su2double **val_Jacobian_ii, su2double **val_Jacobian_ij,
                       su2double **val_Jacobian_ji, su2double **val_Jacobian_jj, CConfig *config) override;
};
