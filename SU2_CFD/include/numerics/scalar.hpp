/*!
 * \file scalar.hpp
 * \brief Delarations of numerics classes for scalar transport problems.
 * \author F. Palacios, T. Economon, D. Mayer
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

/*!
 * \class CUpwScalar
 * \brief Template class for scalar upwind fluxes between nodes i and j.
 * \details This class serves as a template for the scalar upwinding residual
 *   classes.  The general structure of a scalar upwinding calculation is the
 *   same for many different  models, which leads to a lot of repeated code.
 *   By using the template design pattern, these sections of repeated code are
 *   moved to this shared base class, and the specifics of each model
 *   are implemented by derived classes.  In order to add a new residual
 *   calculation for a convection residual, extend this class and implement
 *   the pure virtual functions with model-specific behavior.
 * \ingroup ConvDiscr
 * \author C. Pederson, A. Bueno., and A. Campos.
 */
class CUpwtransportedScalar : public CNumerics {
protected:
  su2double
  q_ij = 0.0,                  /*!< \brief Projected velocity at the face. */
  a0 = 0.0,                    /*!< \brief The maximum of the face-normal velocity and 0 */
  a1 = 0.0,                    /*!< \brief The minimum of the face-normal velocity and 0 */
  *Flux = nullptr,             /*!< \brief Final result, diffusive flux/residual. */
  **Jacobian_i = nullptr,      /*!< \brief Flux Jacobian w.r.t. node i. */
  **Jacobian_j = nullptr;      /*!< \brief Flux Jacobian w.r.t. node j. */

  const bool implicit = false, incompressible = false, dynamic_grid = false;

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method, derived classes
   *        compute the Flux and its Jacobians via this method.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(const CConfig* config) = 0;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwtransportedScalar(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwtransportedScalar(void);

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config);

};

/*!
 * \class CUpwtransportedScalar_General
 * \brief Class for a general scalar upwind solver for an arbitrary number of scalar transport eqns.
 * \ingroup ConvDiscr
 * \author T. Economon
 */
class CUpwtransportedScalar_General : public CUpwtransportedScalar {
private:

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn();

  /*!
   * \brief SST specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                          su2double **Jacobian_j, CConfig *config);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(const CConfig* config) override;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwtransportedScalar_General(unsigned short val_nDim,
                     unsigned short val_nVar,
                     CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwtransportedScalar_General(void);
};

/*!
 * \class CAvgGradtransportedScalar
 * \brief Template class for computing viscous residual of scalar values
 * \ingroup ViscDiscr
 * \author C. Pederson, T. Economon
 */
class CAvgGradtransportedScalar : public CNumerics {
private:

  /*!
   * \brief A pure virtual function; Adds any extra variables to AD
   */
  virtual void ExtraADPreaccIn() = 0;

  /*!
   * \brief Model-specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  virtual void FinishResidualCalc(su2double *val_residual,
                                  su2double **Jacobian_i,
                                  su2double **Jacobian_j,
                                  CConfig *config) = 0;

protected:
  bool implicit, incompressible;
  bool correct_gradient;
  unsigned short iVar, iDim;
  su2double **Mean_GradScalarVar;   /*!< \brief Average of scalar var gradients at cell face */
  su2double *Edge_Vector,           /*!< \brief Vector from node i to node j. */
  *Proj_Mean_GradScalarVar_Normal,  /*!< \brief Mean_GradScalarVar DOT normal */
  *Proj_Mean_GradScalarVar_Edge,    /*!< \brief Mean_GradScalarVar DOT Edge_Vector */
  *Proj_Mean_GradScalarVar;         /*!< \brief Mean_GradScalarVar DOT normal, corrected if required */
  su2double dist_ij_2,              /*!< \brief |Edge_Vector|^2 */
  proj_vector_ij;                   /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradtransportedScalar(unsigned short val_nDim, unsigned short val_nVar,
                  bool correct_gradient, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradtransportedScalar(void);

  /*!
   * \brief Compute the viscous residual using an average of gradients without correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i,
                       su2double **Jacobian_j, CConfig *config);
};

  /*!
 * \class CAvgGradtransportedScalar_General
 * \brief Class for computing viscous term using average of gradients for a passive scalar.
 * \ingroup ViscDiscr
 * \author T. Economon
 */
class CAvgGradtransportedScalar_General : public CAvgGradtransportedScalar {
protected:
  su2double *Mean_Diffusivity;   /*!< \brief Average of mass diffusivities at cell face */

private:

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief SA specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                          su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradtransportedScalar_General(unsigned short val_nDim, unsigned short val_nVar,
                         bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradtransportedScalar_General(void);
};


/*!
 * \class CSourcePieceWise_Scalar
 * \brief Class for integrating the source terms of scalar transport equations.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourcePieceWise_Scalar : public CNumerics {
private:
  bool incompressible;  /*!< \brief Flag defining compressibility. */
  bool implicit;        /*!< \brief Flag defining implicit scheme. */

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_Scalar(unsigned short val_nDim,
                          unsigned short val_nVar,
                          CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_Scalar(void);
  
  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i,
                       su2double **val_Jacobian_j,
                       CConfig *config);

};

/*!
 * \class CSourceAxisymmetric_Scalar
 * \brief Class for source term for solving scalar axisymmetric problems.
 * \ingroup SourceDiscr
 * \author T. Economon
 */
class CSourceAxisymmetric_Scalar : public CNumerics {
  bool implicit, /*!< \brief Implicit calculation. */
  viscous, /*!< \brief Viscous incompressible flows. */
  energy; /*!< \brief computation with the energy equation. */

  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourceAxisymmetric_Scalar(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSourceAxisymmetric_Scalar(void);
  
  /*!
   * \brief Residual of the rotational frame source term.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **Jacobian_i, CConfig *config);
  
};