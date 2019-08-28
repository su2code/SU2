/*!
 * \file CFEAVariable.hpp
 * \brief Class for defining the variables of the FEM structural problem.
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
 * \class CFEAVariable
 * \brief Class for defining the variables of the FEM structural problem.
 * \ingroup Structural Finite Element Analysis Variables
 * \author F. Palacios, R. Sanchez.
 * \version 6.2.0 "Falcon"
 */
class CFEAVariable : public CVariable {
protected:

  su2double *Stress;              /*!< \brief Stress tensor. */

  su2double *Residual_Ext_Body;   /*!< \brief Term of the residual due to body forces */

  su2double VonMises_Stress;      /*!< \brief Von Mises stress. */

  su2double *Solution_Vel,        /*!< \brief Velocity of the nodes. */
  *Solution_Vel_time_n;           /*!< \brief Velocity of the nodes at time n. */

  su2double *Solution_Accel,      /*!< \brief Acceleration of the nodes. */
  *Solution_Accel_time_n;         /*!< \brief Acceleration of the nodes at time n. */

  su2double *Solution_Pred,       /*!< \brief Predictor of the solution for FSI purposes */
  *Solution_Pred_Old;             /*!< \brief Predictor of the solution at time n for FSI purposes */

  su2double *Reference_Geometry;  /*!< \brief Reference solution for optimization problems */

  su2double *Prestretch;          /*!< \brief Prestretch geometry */

  su2double* Solution_BGS_k;      /*!< \brief Old solution container for BGS iterations ---*/


public:

  /*!
   * \brief Constructor of the class.
   */
  CFEAVariable(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAVariable(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEAVariable(void);

  /*!
   * \brief Get the value of the stress.
   * \return Value of the stress.
   */
  inline su2double *GetStress_FEM(void) {return Stress; }

  /*!
   * \brief Set the value of the stress at the node
   * \param[in] iVar - index of the stress term
   * \param[in] val_stress - value of the stress
   */
  inline void SetStress_FEM(unsigned short iVar, su2double val_stress) {Stress[iVar] = val_stress; }

  /*!
   * \brief Add a certain value to the value of the stress at the node
   * \param[in] iVar - index of the stress term
   * \param[in] val_stress - value of the stress
   */
  inline void AddStress_FEM(unsigned short iVar, su2double val_stress) {Stress[iVar] += val_stress; }

  /*!
   * \brief Add body forces to the residual term.
   */
  inline void Add_BodyForces_Res(su2double *val_bodyForce) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Residual_Ext_Body[iVar] += val_bodyForce[iVar];
  }

  /*!
   * \brief Clear the surface load residual
   */
  inline void Clear_BodyForces_Res(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_Ext_Body[iVar] = 0.0;
  }

  /*!
   * \brief Get the body forces.
   */
  inline su2double Get_BodyForces_Res(unsigned short iVar) {return Residual_Ext_Body[iVar];}

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_time_n[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_time_n(su2double *val_solution_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_time_n[iVar] = val_solution_time_n[iVar];
  }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_time_n(unsigned short val_var, su2double val_solution) {
    Solution_time_n[val_var] = val_solution;
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis).
   * \param[in] val_solution - Solution of the problem (velocity).
   */
  void SetSolution_Vel(su2double *val_solution_vel) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Vel[iVar] = val_solution_vel[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Vel(unsigned short val_var, su2double val_solution_vel) {Solution_Vel[val_var] = val_solution_vel; }

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  void SetSolution_Vel_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Vel_time_n[iVar] = Solution_Vel[iVar];
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  void SetSolution_Vel_time_n(su2double *val_solution_vel_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Vel_time_n[iVar] = val_solution_vel_time_n[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_old - Value of the old solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Vel_time_n(unsigned short val_var, su2double val_solution_vel_time_n) {Solution_Vel_time_n[val_var] = val_solution_vel_time_n; }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_Vel(unsigned short val_var) {return Solution_Vel[val_var]; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline su2double *GetSolution_Vel(void) {return Solution_Vel; }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Vel_time_n(unsigned short val_var) {return Solution_Vel_time_n[val_var]; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_Vel_time_n(void) {return Solution_Vel_time_n; }

  /*!
   * \brief Set the value of the acceleration (Structural Analysis).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline void SetSolution_Accel(su2double *val_solution_accel) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Accel[iVar] = val_solution_accel[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Accel(unsigned short val_var, su2double val_solution_accel) {Solution_Accel[val_var] = val_solution_accel;}

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Accel_time_n(su2double *val_solution_accel_time_n) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Accel_time_n[iVar] = val_solution_accel_time_n[iVar];
  }

  /*!
   * \brief Set the value of the acceleration (Structural Analysis) at time n.
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline void SetSolution_Accel_time_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Accel_time_n[iVar] = Solution_Accel[iVar];
  }

  /*!
   * \overload
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution_old - Value of the old solution for the index <i>val_var</i>.
   */
  inline void SetSolution_Accel_time_n(unsigned short val_var, su2double val_solution_accel_time_n) {Solution_Accel_time_n[val_var] = val_solution_accel_time_n; }

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] val_var - Index of the variable.
   * \return Value of the solution for the index <i>val_var</i>.
   */
  inline su2double GetSolution_Accel(unsigned short val_var) {return Solution_Accel[val_var]; }

  /*!
   * \brief Get the solution of the problem.
   * \return Pointer to the solution vector.
   */
  inline su2double *GetSolution_Accel(void) {return Solution_Accel; }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Accel_time_n(unsigned short val_var) {return Solution_Accel_time_n[val_var]; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_Accel_time_n(void) {return Solution_Accel_time_n; }

  /*!
   * \brief Set the value of the solution predictor.
   */
  inline void SetSolution_Pred(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Pred[iVar] = Solution[iVar];
  }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Pred(su2double *val_solution_pred) {Solution_Pred = val_solution_pred;  }

  /*!
   * \brief  Set the value of the predicted solution.
   * \param[in] val_var - Index of the variable
   * \param[in] val_solution_pred - Value of the predicted solution.
   */
  inline void SetSolution_Pred(unsigned short val_var, su2double val_solution_pred) {Solution_Pred[val_var] = val_solution_pred;  }

  /*!
   * \brief Get the value of the solution predictor.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Pred(unsigned short val_var) {return Solution_Pred[val_var]; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_Pred(void) {return Solution_Pred; }

  /*!
   * \brief Set the value of the solution predictor.
   */
  inline void SetSolution_Pred_Old(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution_Pred_Old[iVar] = Solution_Pred[iVar];
  }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Pred_Old(su2double *val_solution_pred_Old) {Solution_Pred_Old = val_solution_pred_Old;  }

  /*!
   * \brief  A virtual member. Set the value of the old solution predicted.
   * \param[in] val_var - Index of the variable
   * \param[in] val_solution_pred_old - Value of the old predicted solution.
   */
  inline void SetSolution_Pred_Old(unsigned short val_var, su2double val_solution_pred_old) {Solution_Pred_Old[val_var] = val_solution_pred_old;  }

  /*!
   * \brief Get the value of the solution predictor.
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Pred_Old(unsigned short val_var) {return Solution_Pred_Old[val_var]; }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline su2double *GetSolution_Pred_Old(void) {return Solution_Pred_Old; }

  /*!
   * \brief A virtual member.
   */
  inline void SetPrestretch(unsigned short iVar, su2double val_prestretch) {Prestretch[iVar] = val_prestretch;}

  /*!
   * \brief A virtual member.
   */
  inline su2double *GetPrestretch(void) {return Prestretch; }

  /*!
   * \brief A virtual member.
   */
  inline su2double GetPrestretch(unsigned short iVar) {return Prestretch[iVar]; }

  /*!
   * \brief Set the value of the Von Mises stress.
   * \param[in] val_stress - Value of the Von Mises stress.
   */
  inline void SetVonMises_Stress(su2double val_stress) {VonMises_Stress = val_stress; }

  /*!
   * \brief Get the value of the Von Mises stress.
   * \return Value of the Von Mises stress.
   */
  inline su2double GetVonMises_Stress(void) {return VonMises_Stress; }

  /*!
   * \brief Set the reference geometry.
   * \return Pointer to the solution (at time n) vector.
   */
  inline void SetReference_Geometry(unsigned short iVar, su2double ref_geometry) {Reference_Geometry[iVar] = ref_geometry;}

  /*!
   * \brief Get the pointer to the reference geometry
   */
  inline su2double *GetReference_Geometry(void) {return Reference_Geometry; }

  /*!
   * \brief Get the value of the reference geometry for the coordinate iVar
   */
  inline su2double GetReference_Geometry(unsigned short iVar) {return Reference_Geometry[iVar]; }

  /*!
   * \brief Register the variables in the solution time_n array as input/output variable.
   * \param[in] input - input or output variables.
   */
  inline void Register_femSolution_time_n(void) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_time_n[iVar]);
  }

  /*!
   * \brief Register the variables in the velocity array as input/output variable.
   * \param[in] input - input or output variables.
   */
  inline void RegisterSolution_Vel(bool input) {
	  if (input) {
	    for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterInput(Solution_Vel[iVar]);
	  }
	  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterOutput(Solution_Vel[iVar]);}
  }

  /*!
   * \brief Register the variables in the velocity time_n array as input/output variable.
   */
  inline void RegisterSolution_Vel_time_n(void) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_Vel_time_n[iVar]);
  }

  /*!
   * \brief Register the variables in the acceleration array as input/output variable.
   * \param[in] input - input or output variables.
   */
  inline void RegisterSolution_Accel(bool input) {
	  if (input) {
	    for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterInput(Solution_Accel[iVar]);
	  }
	  else { for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      AD::RegisterOutput(Solution_Accel[iVar]);}
  }

  /*!
   * \brief Register the variables in the acceleration time_n array as input/output variable.
   */
  inline void RegisterSolution_Accel_time_n(void){
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	    AD::RegisterInput(Solution_Accel_time_n[iVar]);
  }

  /*!
   * \brief Set the velocity adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_Vel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_Vel[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the velocity adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_Vel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel[iVar]);
  }

  /*!
   * \brief Set the velocity adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  void SetAdjointSolution_Vel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_Vel_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the velocity adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_Vel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Vel_time_n[iVar]);
  }

  /*!
   * \brief Set the acceleration adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void SetAdjointSolution_Accel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_Accel[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the acceleration adjoint values of the solution.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_Accel(su2double *adj_sol) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel[iVar]);
  }

  /*!
   * \brief Set the acceleration adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  void SetAdjointSolution_Accel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
      SU2_TYPE::SetDerivative(Solution_Accel_time_n[iVar], SU2_TYPE::GetValue(adj_sol[iVar]));
  }

  /*!
   * \brief Get the acceleration adjoint values of the solution at time n.
   * \param[in] adj_sol - The adjoint values of the solution.
   */
  inline void GetAdjointSolution_Accel_time_n(su2double *adj_sol) {
	  for (unsigned short iVar = 0; iVar < nVar; iVar++)
	      adj_sol[iVar] = SU2_TYPE::GetDerivative(Solution_Accel_time_n[iVar]);
  }

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_k(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_BGS_k[iVar] = Solution[iVar];
  }

  /*!
   * \brief Get the value of the solution in the previous BGS subiteration.
   * \param[out] val_solution - solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution_k(unsigned short iDim) {return Solution_BGS_k[iDim];}

};
