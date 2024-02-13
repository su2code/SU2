/*!
 * \file pbflow.hpp
 * \brief Delarations of numerics classes for heat transfer problems.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
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
 * \class CUpwPB_Flow
 * \brief Class for solving an approximate Riemann solver of Roe for the pressure based incompressible flow equations.
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CUpwPB_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  bool gravity;
  su2double Froude, Upw_i, Upw_j;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity, *Velocity_upw;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, MeanDensity, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
  ProjVelocity, FaceVel, Face_Flux;
  unsigned short iDim, iVar, jVar, kVar;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;
  su2double **Jacobian_upw = nullptr;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwPB_Flow(void);
  
  /*!
   * \brief Compute the Roe's flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
  
  /*!
   * \brief Set the value of face velocity. This is used as a proxy for massflux at a face.
   * \param[in] val_FaceVel.
   */
  inline void SetFaceVel(su2double val_FaceVel){ FaceVel = val_FaceVel; }
};


/*!
 * \class CCentLaxArtComp_Flow
 * \brief Class for computing the Lax-Friedrich centered scheme (artificial compressibility).
 * \ingroup ConvDiscr
 * \author F. Palacios
 */
class CCentPB_Flow : public CNumerics {
private:
  bool implicit, dynamic_grid;
  bool gravity;
  su2double Froude, Upw_i, Upw_j;
  su2double *Diff_U;
  su2double *Velocity_i, *Velocity_j, *MeanVelocity, *Velocity_upw;
  su2double *ProjFlux_i, *ProjFlux_j;
  su2double *Lambda, *Epsilon;
  su2double **P_Tensor, **invP_Tensor,**val_Jacobian_upw;
  su2double Proj_ModJac_Tensor_ij, Pressure_i,
  Pressure_j, MeanDensity, MeanSoundSpeed, MeanPressure, MeanBetaInc2,
  ProjVelocity, FaceVel;
  unsigned short iDim, iVar, jVar, kVar;
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;
  su2double **Jacobian_upw = nullptr;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CCentPB_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CCentPB_Flow(void);
  
  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
};

/*!
 * \class CAvgGradInc_Flow
 * \brief Class for computing viscous term using an average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, F. Palacios, T. Economon
 */
class CAvgGradPBInc_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,           /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j;           /*!< \brief Primitives variables at point i and j. */
  su2double **Mean_GradPrimVar,          /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, /*!< \brief Mean value of the viscosity. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij,              /*!< \brief Length of the edge and face. */
  proj_vector_ij;                  /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  bool implicit;        /*!< \brief Implicit calculus. */
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;

public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradPBInc_Flow(void);
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
};


/*!
 * \class CAvgGradCorrectedPBInc_Flow
 * \brief Class for computing viscous term using an average of gradients with correction (incompressible).
 * \ingroup ViscDiscr
 */
class CAvgGradCorrectedPBInc_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;  /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar;           /*!< \brief Mean primitive variables. */
  su2double *PrimVar_i, *PrimVar_j,      /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,      /*!< \brief Mean value of the viscosity. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2,          /*!< \brief Length of the edge and face. */
  proj_vector_ij;     /*!< \brief (Edge_Vector DOT normal)/|Edge_Vector|^2 */
  bool implicit;      /*!< \brief Implicit calculus. */
  
  su2double *Flux = nullptr;
  su2double **Jacobian_i = nullptr;
  su2double **Jacobian_j = nullptr;
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrectedPBInc_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrectedPBInc_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override;
};
