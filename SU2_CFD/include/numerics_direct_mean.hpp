/*!
 * \file numerics_direct_mean.hpp
 * \brief
 * \author C. Pederson
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "numerics_structure.hpp"

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 */
class CAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;     /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,           /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,           /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,             /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij,           /*!< \brief Length of the edge and face. */
  Mean_TauWall,     /*!< \brief Mean wall shear stress (wall functions). */
  TauWall_i, TauWall_j;  /*!< \brief Wall shear stress at point i and j (wall functions). */
  bool implicit; /*!< \brief Implicit calculus. */
  su2double* heat_flux_vector;

  void GetTau(const su2double *val_primvar,
              su2double **val_gradprimvar,
              const su2double val_turb_ke,
              const su2double val_laminar_viscosity,
              const su2double val_eddy_viscosity);

  void AddQCR(su2double **val_gradprimvar);

  void AddTauWall(const su2double *val_normal,
                  const su2double val_tau_wall);

  void GetHeatFluxVector(su2double **val_gradprimvar,
                         const su2double val_laminar_viscosity,
                         const su2double val_eddy_viscosity);

  void GetViscousProjFlux(const su2double *val_primvar,
                          const su2double *val_normal);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);


  /*!
   * \brief Set the value of the wall shear stress at point i and j (wall functions).
   * \param[in] val_tauwall_i - Value of the wall shear stress at point i.
   * \param[in] val_tauwall_j - Value of the wall shear stress at point j.
   */
  void SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j);
};

/*!
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 */

class CGeneralAvgGrad_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;     /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,           /*!< \brief Mean primitive variables. */
  *Mean_SecVar,                   /*!< \brief Mean secondary variables. */
  *PrimVar_i, *PrimVar_j,           /*!< \brief Primitives variables at point i and 1. */
  **Mean_GradPrimVar,             /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,             /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                               /*!< \brief Mean value of the Cp. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;            /*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

/*!
 * \class CAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 */
class CAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;    /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,          /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,        /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                  /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,         /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2,           /*!< \brief Length of the edge and face. */
  TauWall_i, TauWall_j,  /*!< \brief Wall shear stress at point i and j (wall functions). */
  Mean_TauWall;     /*!< \brief Mean wall shear stress (wall functions). */
  bool implicit;      /*!< \brief Implicit calculus. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGradCorrected_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);


  /*!
   * \brief Set the value of the wall shear stress at point i and j (wall functions).
   * \param[in] val_tauwall_i - Value of the wall shear stress at point i.
   * \param[in] val_tauwall_j - Value of the wall shear stress at point j.
   */
  void SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j);
};


/*!
 * \class CGeneralAvgGradCorrected_Flow
 * \brief Class for computing viscous term using the average of gradients with a correction.
 * \ingroup ViscDiscr
 * \author M. Pini, S. Vitale
 */
class CGeneralAvgGradCorrected_Flow : public CNumerics {
private:
  unsigned short iDim, iVar, jVar;    /*!< \brief Iterators in dimension an variable. */
  su2double *Mean_PrimVar,          /*!< \brief Mean primitive variables. */
  *Mean_SecVar,                  /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,            /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                  /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,      /*!< \brief Mean value of the laminar viscosity. */
  Mean_Eddy_Viscosity,         /*!< \brief Mean value of the eddy viscosity. */
  Mean_Thermal_Conductivity,   /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp,                     /*!< \brief Mean value of the specific heat. */
  Mean_turb_ke,         /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij_2;           /*!< \brief Length of the edge and face. */
  bool implicit;      /*!< \brief Implicit calculus. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGradCorrected_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGradCorrected_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients with correction.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);
};

// TODO: Move these to an inline file after refactoring at the end
inline void CAvgGrad_Flow::SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) {
  TauWall_i = val_tauwall_i;
  TauWall_j = val_tauwall_j;
}

inline void CAvgGradCorrected_Flow::SetTauWall(su2double val_tauwall_i, su2double val_tauwall_j) {
  TauWall_i = val_tauwall_i;
  TauWall_j = val_tauwall_j;
}
