/*!
 * \file transition.hpp
 * \brief Delarations of numerics classes for transition problems.
 * \author F. Palacios, T. Economon
 * \version 7.0.2 "Blackbird"
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

#include "CNumerics.hpp"
#include "turbulent/turb_convection.hpp"
#include "turbulent/turb_diffusion.hpp"

/*!
 * \class CUpwSca_TransLM
 * \brief Class for doing a scalar upwind solver for the LM transition model.
 * \ingroup ConvDiscr
 * \author A. Aranake, E. van der Weide.
 */
class CUpwSca_TransLM : public CUpwScalar {
private:
  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn() override;

  /*!
   * \brief LM specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig *config) override;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
};


/*!
 * \class CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradient with correction (LM transition model).
 * \ingroup ViscDiscr
 * \author A. Bueno, E. van der Weide.
 */
class CAvgGrad_TransLM : public CAvgGrad_Scalar {
private:
  su2double sigma_intermittency; /*!< \brief Constant for the viscous term of the intermittency equation. */
  su2double sigma_Re_theta;      /*!< \brief Constant for the viscous term of the Re_theta equation. */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void) override;

  /*!
   * \brief LM specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig *config) override;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim     - Number of dimensions of the problem.
   * \param[in] val_nVar     - Number of variables of the problem.
   * \param[in] constants    - Array containing the constants used in the LM model.
   * \param[in] correct_grad - Whether or not the gradients must be corrected.
   * \param[in] config       - Definition of the particular problem.
   */
  CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                   const su2double* constants, bool correct_grad, CConfig *config);
};

/*!
 * \class CSourcePieceWise_TransLM
 * \brief Class for integrating the source terms of the LM transition model equation.
 * \ingroup SourceDiscr
 * \author E. van der Weide
 */
class CSourcePieceWise_TransLM : public CNumerics {
private:

  su2double ca1, ca2, ce1, ce2, cthetat;        /*!< \brief Constants in the source term of the LM model. */
  su2double C_crossflow;                        /*!< \brief Constants in the source term of the LM model
                                                            related to cross flow instabilities. */
  su2double Flength_CF, C_Fonset1_CF, CHe_max;  /*!< \brief Constants in the source term of the LM model
                                                            related to cross flow instabilities. */

  bool incompressible;  /*!< \brief Whether or not an incompressible simulation is carried out. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim  - Number of dimensions of the problem.
   * \param[in] val_nVar  - Number of variables of the problem.
   * \param[in] constants - Constants used in the LM transition model.
   * \param[in] config    - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TransLM(void);

  /*!
   * \brief Static member function to compute the Re_theta from turbulence intensity.
   * \param[in] var_tu - Turbulence intensity.
   * \return Value of the Re_theta.
   */
  static su2double GetREth(const su2double var_tu);

  /*!
   * \brief Static member function to compute the critical Reynolds_theta. This is the
            Reynolds number where the intermittency starts to increase in
            the boundary layer.
   * \param[in] var_Re_theta - Value of Reynolds theta. The critical Reynolds number
                               is a correlation based on this value.
   * \return Value of the critical Re_theta.
   */
  static su2double GetREth_crit(const su2double var_Re_theta);

  /*!
   * \brief Static member function to compute the empirical correlation for Flength,
            which controls the length of the transition region.
   * \param[in] var_Re_theta - Value of Reynolds theta. Flength is a
                               correlation based on this value.
   * \return Value of the Flength.
   */
  static su2double GetFlength(const su2double var_Re_theta);

  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};
