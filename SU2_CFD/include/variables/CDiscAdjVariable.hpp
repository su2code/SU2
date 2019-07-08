/*!
 * \file CDiscAdjVariable.hpp
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
 * \class CDiscAdjVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring.
 */
class CDiscAdjVariable : public CVariable {
private:
  su2double* Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  su2double* Solution_Direct;
  su2double* DualTime_Derivative;
  su2double* DualTime_Derivative_n;

  su2double* Cross_Term_Derivative;
  su2double* Geometry_CrossTerm_Derivative;
  su2double* Geometry_CrossTerm_Derivative_Flow;

  su2double* Solution_Geometry;
  su2double* Solution_Geometry_Old;
  su2double* Geometry_Direct;

  su2double* Solution_BGS;
  su2double* Solution_BGS_k;
  su2double* Solution_Geometry_BGS_k;

public:
  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjVariable(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjVariable(void);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the adjoint value (initialization value).
   * \param[in] val_ndim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjVariable(su2double *val_solution, unsigned short val_ndim, unsigned short val_nvar, CConfig *config);

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

  inline void SetDual_Time_Derivative(unsigned short iVar, su2double der) {DualTime_Derivative[iVar] = der;}

  inline void SetDual_Time_Derivative_n(unsigned short iVar, su2double der) {DualTime_Derivative_n[iVar] = der;}

  inline su2double GetDual_Time_Derivative(unsigned short iVar) {return DualTime_Derivative[iVar];}

  inline su2double GetDual_Time_Derivative_n(unsigned short iVar) {return DualTime_Derivative_n[iVar];}

  inline void SetSolution_Direct(su2double *val_solution_direct) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_Direct[iVar] = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct() {return Solution_Direct; }

  /*!
   * \brief Set the restart geometry (coordinate of the converged solution)
   * \param[in] val_geometry_direct - Value of the restart coordinate.
   */
  inline void SetGeometry_Direct(su2double *val_geometry_direct) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Geometry_Direct[iDim] = val_geometry_direct[iDim];
  }

  /*!
   * \brief Get the restart geometry (coordinate of the converged solution).
   * \return Pointer to the restart coordinate vector.
   */
  inline su2double *GetGeometry_Direct(void) {return Geometry_Direct;}

  /*!
   * \brief Get the restart geometry (coordinate of the converged solution).
   * \return Coordinate val_dim of the geometry_direct vector.
   */
  inline su2double GetGeometry_Direct(unsigned short val_dim) {return Geometry_Direct[val_dim]; }

  /*!
   * \brief Get the geometry solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_Geometry(unsigned short val_var) {return Solution_Geometry[val_var];}

  /*!
   * \brief Set the value of the mesh solution (adjoint).
   * \param[in] val_solution_geometry - Solution of the problem (acceleration).
   */
  inline void SetSolution_Geometry(su2double *val_solution_geometry) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Geometry[iDim] = val_solution_geometry[iDim];
  }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution_geometry - Solution of the problem (acceleration).
   */
  inline void SetSolution_Geometry(unsigned short val_var, su2double val_solution_geometry) {
    Solution_Geometry[val_var] = val_solution_geometry;
  }

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
   * \brief Get the mesh cross term derivative from the flow solution.
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetGeometry_CrossTerm_Derivative_Flow(unsigned short val_var) {return Geometry_CrossTerm_Derivative_Flow[val_var];}

  /*!
   * \brief Set the value of the mesh cross term derivative from the flow solution (adjoint).
   * \param[in] der - cross term derivative.
   */
  inline void SetGeometry_CrossTerm_Derivative_Flow(unsigned short iDim, su2double der) {Geometry_CrossTerm_Derivative_Flow[iDim] = der;}

  /*!
   * \brief Set the value of the mesh solution (adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline void Set_OldSolution_Geometry(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Geometry_Old[iDim] = Solution_Geometry[iDim];
  }

  /*!
   * \brief Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline su2double Get_OldSolution_Geometry(unsigned short iDim) {return Solution_Geometry_Old[iDim];}

  /*!
   * \brief Set the value of the adjoint solution in the current BGS subiteration.
   */
  inline void Set_BGSSolution(unsigned short iDim, su2double val_solution) {Solution_BGS[iDim] = val_solution;}

  /*!
   * \brief Set the value of the adjoint solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_k(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_BGS_k[iVar] = Solution_BGS[iVar];
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

  /*!
   * \brief Set the value of the adjoint geometry solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_Geometry(void) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Geometry_BGS_k[iDim] = Solution_Geometry[iDim];
  }

  /*!
   * \brief Get the value of the adjoint geometry solution in the previous BGS subiteration.
   * \param[out] val_solution - geometrical adjoint solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution_Geometry(unsigned short iDim) {return Solution_Geometry_BGS_k[iDim];}

  /*!
   * \brief Set the contribution of crossed terms into the derivative.
   */
  inline void SetCross_Term_Derivative(unsigned short iVar, su2double der) {Cross_Term_Derivative[iVar] = der; }

  /*!
   * \brief Get the contribution of crossed terms into the derivative.
   */
  inline su2double GetCross_Term_Derivative(unsigned short iVar) {return Cross_Term_Derivative[iVar]; }

};
