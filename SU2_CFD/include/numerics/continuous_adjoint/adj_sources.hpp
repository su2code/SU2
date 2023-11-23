/*!
 * \file adj_sources.hpp
 * \brief Declarations of numerics classes for continuous adjoint
 *        source term integration. Implemented in adj_sources.cpp.
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
 * \class CSourceAxisymmetric_AdjFlow
 * \brief Class for source term for solving axisymmetric problems.
 * \ingroup SourceDiscr
 * \author F. Palacios
 */
class CSourceAxisymmetric_AdjFlow : public CNumerics {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_AdjFlow(void) override;

  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config) override;
};

class CSourceConservative_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *Residual_i, *Residual_j, *Mean_Residual;
  su2double **Mean_PrimVar_Grad;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjFlow(void) override;

  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config) override;
};

/*!
 * \class CSourceRotatingFrame_AdjFlow
 * \brief Source term class for rotating frame adjoint.
 * \ingroup SourceDiscr
 * \author T. Economon.
 */
class CSourceRotatingFrame_AdjFlow : public CNumerics {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceRotatingFrame_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceRotatingFrame_AdjFlow(void) override;

  /*!
   * \brief Residual of the adjoint rotating frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, CConfig *config) override;
};

/*!
 * \class CSourceViscous_AdjFlow
 * \brief Class for source term integration in adjoint problem.
 * \ingroup SourceDiscr
 * \author F. Palacios
 */
class CSourceViscous_AdjFlow : public CNumerics {
private:
  su2double *Velocity, *GradDensity, *GradInvDensity, *dPoDensity2, *alpha, *beta, *Sigma_5_vec;
  su2double **GradVel_o_Rho, **sigma, **Sigma_phi, **Sigma_5_Tensor, **Sigma;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceViscous_AdjFlow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceViscous_AdjFlow(void) override;

  /*!
   * \brief Source term integration of the flow adjoint equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, CConfig *config) override;

};

/*!
 * \class CSourceConservative_AdjTurb
 * \brief Class for source term integration in adjoint turbulent problem using a conservative scheme.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 */
class CSourceConservative_AdjTurb : public CNumerics {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceConservative_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourceConservative_AdjTurb(void) override;

  /*!
   * \brief Source term integration using a conservative scheme.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;
};

/*!
 * \class CSourcePieceWise_AdjTurb
 * \brief Class for source term integration of the adjoint turbulent equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 */
class CSourcePieceWise_AdjTurb : public CNumerics {
private:
  su2double **tau, *Velocity;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_AdjTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_AdjTurb(void) override;

  /*!
   * \brief Source term integration of the adjoint turbulence equation.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;
};
