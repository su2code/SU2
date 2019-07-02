/*!
 * \file CIncEulerVariable.hpp
 * \brief Class for defining the variables of the incompressible Euler solver.
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
 * \class CIncEulerVariable
 * \brief Class for defining the variables of the incompressible Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncEulerVariable : public CVariable {
protected:
  su2double Velocity2;      /*!< \brief Square of the velocity vector. */

  /*--- Primitive variable definition ---*/

  su2double *Primitive;  /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  su2double **Gradient_Primitive;  /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */

  /*--- Old solution container for BGS iterations ---*/

  su2double* Solution_BGS_k;

  /*--- Old density for variable density turbulent flows (SST). ---*/

  su2double Density_Old;

public:

  /*!
   * \brief Constructor of the class.
   */
  CIncEulerVariable(void);

  /*!
   * \overload
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_temperature - Value of the temperature (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncEulerVariable(su2double val_pressure, su2double *val_velocity, su2double val_temperature, unsigned short val_nDim,
                    unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIncEulerVariable(void);

  /*!
   * \brief Set to zero the gradient of the primitive variables.
   */
  void SetGradient_PrimitiveZero(unsigned short val_primvar);

  /*!
   * \brief Add <i>val_value</i> to the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to add to the gradient of the primitive variables.
   */
  inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Primitive[val_var][val_dim] += val_value; }

  /*!
   * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
   */
  inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Primitive[val_var][val_dim] -= val_value; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) {return Gradient_Primitive[val_var][val_dim]; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned short val_var) {return Limiter_Primitive[val_var]; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Primitive[val_var][val_dim] = val_value; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetLimiter_Primitive(unsigned short val_var, su2double val_value) {Limiter_Primitive[val_var] = val_value; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double **GetGradient_Primitive(void) {return Gradient_Primitive; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double *GetLimiter_Primitive(void) {return Limiter_Primitive; }

  /*!
   * \brief Set the value of the pressure.
   */
  inline void SetPressure(void) {Primitive[0] = Solution[0];}

  /*!
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetPrimitive(unsigned short val_var) {return Primitive[val_var]; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitive(unsigned short val_var, su2double val_prim) {Primitive[val_var] = val_prim; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitive(su2double *val_prim) {
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
      Primitive[iVar] = val_prim[iVar];
  }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(void) {return Primitive; }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   */
  inline bool SetDensity(su2double val_density) {
    Primitive[nDim+2] = val_density;
    if (Primitive[nDim+2] > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   */
  inline void SetVelocity(void) {
    Velocity2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Primitive[iDim+1] = Solution[iDim+1];
      Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
    }
  }

  /*!
   * \brief Set the value of the temperature for incompressible flows with energy equation.
   */
  inline bool SetTemperature(su2double val_temperature) {
    Primitive[nDim+1] = val_temperature;
    if (Primitive[nDim+1] > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the beta coeffient for incompressible flows.
   */
  inline void SetBetaInc2(su2double val_betainc2) {Primitive[nDim+3] = val_betainc2; }

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(void) {return Velocity2; }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(void) {return Primitive[0]; }

  /*!
   * \brief Get the value of beta squared for the incompressible flow
   * \return Value of beta squared.
   */
  inline su2double GetBetaInc2(void) {return Primitive[nDim+3]; }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(void) {return Primitive[nDim+2]; }

  /*!
   * \brief Get the density of the flow from the previous iteration.
   * \return Old value of the density of the flow.
   */
  inline su2double GetDensity_Old(void) {return Density_Old; }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(void) {return Primitive[nDim+1]; }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline su2double GetVelocity(unsigned short val_dim) {return Primitive[val_dim+1]; }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  su2double GetProjVel(su2double *val_vector);

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(su2double *val_velocity) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Old[iDim+1] = val_velocity[iDim];
  }

  /*!
   * \brief Set all the primitive variables for incompressible flows.
   */
  bool SetPrimVar(CFluidModel *FluidModel);

  /*!
   * \brief Set the specific heat Cp.
   */
  inline void SetSpecificHeatCp(su2double val_Cp) {Primitive[nDim+7] = val_Cp;}

  /*!
   * \brief Set the specific heat Cv.
   */
  inline void SetSpecificHeatCv(su2double val_Cv) {Primitive[nDim+8] = val_Cv;}

  /*!
   * \brief Get the specific heat at constant P of the flow.
   * \return Value of the specific heat at constant P of the flow.
   */
  inline su2double GetSpecificHeatCp(void) {return Primitive[nDim+7]; }

  /*!
   * \brief Get the specific heat at constant V of the flow.
   * \return Value of the specific heat at constant V of the flow.
   */
  inline su2double GetSpecificHeatCv(void) {return Primitive[nDim+8]; }

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
