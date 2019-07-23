/*!
 * \file CDiscAdjFEAVariable.hpp
 * \brief Main class for defining the variables of the adjoint solver.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "CVariable.hpp"

/*!
 * \class CDiscAdjFEAVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring, R. Sanchez.
 * \version 6.2.0 "Falcon"
 */
class CDiscAdjFEAVariable : public CVariable {
private:
  su2double* Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  su2double* Solution_Direct;

  su2double* Dynamic_Derivative;
  su2double* Dynamic_Derivative_n;
  su2double* Dynamic_Derivative_Vel;
  su2double* Dynamic_Derivative_Vel_n;
  su2double* Dynamic_Derivative_Accel;
  su2double* Dynamic_Derivative_Accel_n;

  su2double* Solution_Vel;
  su2double* Solution_Accel;

  su2double* Solution_Vel_time_n;
  su2double* Solution_Accel_time_n;

  su2double* Solution_Old_Vel;
  su2double* Solution_Old_Accel;

  su2double* Solution_Direct_Vel;
  su2double* Solution_Direct_Accel;

  su2double* Cross_Term_Derivative;
  su2double* Geometry_CrossTerm_Derivative;

  su2double* Solution_BGS;
  su2double* Solution_BGS_k;

public:
  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjFEAVariable(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEAVariable(void);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the adjoint value (initialization value).
   * \param[in] val_ndim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEAVariable(su2double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution       - Pointer to the adjoint value (initialization value).
   * \param[in] val_solution_accel - Pointer to the adjoint value (initialization value).
   * \param[in] val_solution_vel   - Pointer to the adjoint value (initialization value).
   * \param[in] val_ndim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEAVariable(su2double *val_solution, su2double *val_solution_accel, su2double *val_solution_vel, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned short iDim, su2double val) {Sensitivity[iDim] = val;}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned short iDim) {return Sensitivity[iDim];}

  inline void SetDynamic_Derivative(unsigned short iVar, su2double der) {Dynamic_Derivative[iVar] = der; }

  inline void SetDynamic_Derivative_n(unsigned short iVar, su2double der) {Dynamic_Derivative_n[iVar] = der; }

  inline su2double GetDynamic_Derivative(unsigned short iVar) {return Dynamic_Derivative[iVar]; }

  inline su2double GetDynamic_Derivative_n(unsigned short iVar) {return Dynamic_Derivative_n[iVar]; }

  inline void SetDynamic_Derivative_Vel(unsigned short iVar, su2double der) {Dynamic_Derivative_Vel[iVar] = der; }

  inline void SetDynamic_Derivative_Vel_n(unsigned short iVar, su2double der) {Dynamic_Derivative_Vel_n[iVar] = der; }

  inline su2double GetDynamic_Derivative_Vel(unsigned short iVar) {return Dynamic_Derivative_Vel[iVar]; }

  inline su2double GetDynamic_Derivative_Vel_n(unsigned short iVar) {return Dynamic_Derivative_Vel_n[iVar]; }

  inline void SetDynamic_Derivative_Accel(unsigned short iVar, su2double der) {Dynamic_Derivative_Accel[iVar] = der; }

  inline void SetDynamic_Derivative_Accel_n(unsigned short iVar, su2double der) {Dynamic_Derivative_Accel_n[iVar] = der; }

  inline su2double GetDynamic_Derivative_Accel(unsigned short iVar) {return Dynamic_Derivative_Accel[iVar]; }

  inline su2double GetDynamic_Derivative_Accel_n(unsigned short iVar) {return Dynamic_Derivative_Accel_n[iVar]; }

  inline void SetSolution_Direct(su2double *val_solution_direct) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct[iVar] = val_solution_direct[iVar];
  }

  inline void SetSolution_Vel_Direct(su2double *val_solution_direct) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct_Vel[iVar] = val_solution_direct[iVar];
  }

  inline void SetSolution_Accel_Direct(su2double *val_solution_direct) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Direct_Accel[iVar] = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct() {return Solution_Direct; }

  inline su2double* GetSolution_Vel_Direct() {return Solution_Direct_Vel; }

  inline su2double* GetSolution_Accel_Direct() {return Solution_Direct_Accel; }

  inline su2double GetSolution_Old_Vel(unsigned short iVar) {return Solution_Old_Vel[iVar]; }

  inline su2double GetSolution_Old_Accel(unsigned short iVar) {return Solution_Old_Accel[iVar]; }

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_Accel(unsigned short val_var) {return Solution_Accel[val_var]; }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Accel_time_n(unsigned short val_var) {return Solution_Accel_time_n[val_var]; }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_Vel(unsigned short val_var) {return Solution_Vel[val_var]; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Vel_time_n(unsigned short val_var) {return Solution_Vel_time_n[val_var]; }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_time_n[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set the value of the acceleration (Structural Analysis - adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline void SetSolution_Accel(su2double *val_solution_accel) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_Accel[iVar] = val_solution_accel[iVar];
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis - adjoint).
   * \param[in] val_solution - Solution of the problem (velocity).
   */
  inline void SetSolution_Vel(su2double *val_solution_vel) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Vel[iVar] = val_solution_vel[iVar];
  }

  /*!
   * \brief Set the value of the adjoint acceleration (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Accel_time_n[iVar] = val_solution_accel_time_n[iVar];
  }

  /*!
   * \brief Set the value of the adjoint velocity (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Vel_time_n[iVar] = val_solution_vel_time_n[iVar];
  }

  /*!
   * \brief Set the value of the old acceleration (Structural Analysis - adjoint).
   * \param[in] val_solution - Old solution of the problem (acceleration).
   */
  inline void Set_OldSolution_Accel(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Old_Accel[iVar] = Solution_Accel[iVar];
  }

  /*!
   * \brief Set the value of the old velocity (Structural Analysis - adjoint).
   * \param[in] val_solution - Old solution of the problem (velocity).
   */
  inline void Set_OldSolution_Vel(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Old_Vel[iVar] = Solution_Vel[iVar];
  }

  /*!
   * \brief Set the contribution of crossed terms into the derivative.
   */
  inline void SetCross_Term_Derivative(unsigned short iVar, su2double der) {Cross_Term_Derivative[iVar] = der; }

  /*!
   * \brief Get the contribution of crossed terms into the derivative.
   */
  inline su2double GetCross_Term_Derivative(unsigned short iVar) {return Cross_Term_Derivative[iVar]; }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetGeometry_CrossTerm_Derivative(unsigned short val_var) {return Geometry_CrossTerm_Derivative[val_var];}

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] der - cross term derivative.
   */
  inline void SetGeometry_CrossTerm_Derivative(unsigned short iDim, su2double der) {Geometry_CrossTerm_Derivative[iDim] = der;}

  /*!
   * \brief Set the value of the adjoint solution in the current BGS subiteration.
   */
  inline void Set_BGSSolution(unsigned short iDim, su2double val_solution) {Solution_BGS[iDim] = val_solution;}

  /*!
   * \brief Set the value of the adjoint solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_k(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_BGS_k[iDim] = Solution_BGS[iDim];
  }

  /*!
   * \brief Get the value of the adjoint solution in the previous BGS subiteration.
   * \param[out] val_solution - adjoint solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution(unsigned short iDim) {return Solution_BGS[iDim];}

  /*!
   * \brief Get the value of the adjoint solution in the previous BGS subiteration.
   * \param[out] val_solution - adjoint solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution_k(unsigned short iDim) {return Solution_BGS_k[iDim];}

};
