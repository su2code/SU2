/*!
 * \file transition.hpp
 * \brief Delarations of numerics classes for transition problems.
 * \author F. Palacios, T. Economon
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class CUpwLin_TransLM
 * \brief Class for performing a linear upwind solver for the Spalart-Allmaras turbulence model equations with transition
 * \ingroup ConvDiscr
 * \author A. Aranake
 */
class CUpwLin_TransLM : public CNumerics {
private:
  su2double *Velocity_i;
  su2double *Velocity_j;
  bool implicit, incompressible;
  su2double Density_i, Density_j, q_ij, a0, a1;
  unsigned short iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwLin_TransLM(void) override;

  /*!
   * \brief Compute the upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;
};

/*!
 * \class CUpwSca_TransLM
 * \brief Class for doing a scalar upwind solver for the Spalart-Allmaras turbulence model equations with transition.
 * \ingroup ConvDiscr
 * \author A. Aranake.
 */
class CUpwSca_TransLM : public CNumerics {
private:
  su2double *Velocity_i, *Velocity_j;
  bool implicit;
  su2double q_ij, a0, a1;
  unsigned short iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TransLM(void) override;

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) override;
};

/*!
 * \class CAvgGrad_TransLM
 * \brief Class for computing viscous term using average of gradients (Spalart-Allmaras Turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGrad_TransLM : public CNumerics {
private:
  su2double **Mean_GradTransVar;
  su2double *Proj_Mean_GradTransVar_Kappa, *Proj_Mean_GradTransVar_Edge;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;
  //su2double dist_ij_2;
  //su2double proj_vector_ij;
  //unsigned short iVar, iDim;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TransLM(void) override;

  /*!
   * \brief Compute the viscous turbulence terms residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) override;
};

/*!
 * \class CAvgGradCorrected_TransLM
 * \brief Class for computing viscous term using average of gradients with correction (Spalart-Allmaras turbulence model).
 * \ingroup ViscDiscr
 * \author A. Bueno.
 */
class CAvgGradCorrected_TransLM : public CNumerics {
private:
  su2double **Mean_GradTurbVar;
  su2double *Proj_Mean_GradTurbVar_Kappa, *Proj_Mean_GradTurbVar_Edge, *Proj_Mean_GradTurbVar_Corrected;
  su2double *Edge_Vector;
  bool implicit, incompressible;
  su2double sigma;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_TransLM(void) override;

  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) override;
};

/*!
 * \class CSourcePieceWise_TransLM
 * \brief Class for integrating the source terms of the Spalart-Allmaras turbulence model equation.
 * \ingroup SourceDiscr
 * \author A. Bueno.
 */
class CSourcePieceWise_TransLM : public CNumerics {
private:

  /*-- SA model constants --*/
  su2double cv1_3;
  su2double k2;
  su2double cb1;
  su2double cw2;
  su2double cw3_6;
  su2double sigma;
  su2double cb2;
  su2double cw1;

  /*-- gamma-theta model constants --*/
  su2double c_e1;
  su2double c_a1;
  su2double c_e2;
  su2double c_a2;
  su2double sigmaf;
  su2double s1;
  su2double c_theta;
  su2double sigmat;

  /*-- Correlation constants --*/
  su2double flen_global;
  su2double alpha_global;
  su2double Vorticity;

  bool implicit;

public:
  bool debugme; // For debugging only, remove this. -AA

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TransLM(void) override;

  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual_TransLM(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config, su2double &gamma_sep) override;

  void CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(const su2double *TransVar_i, const su2double *TransVar_id, su2double *val_residual, su2double *val_residuald, CConfig *config);
};
