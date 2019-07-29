/*!
 * \file CEulerVariable.hpp
 * \brief Class for defining the variables of the compressible Euler solver.
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
 * \class CEulerVariable
 * \brief Class for defining the variables of the compressible Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon
 */
class CEulerVariable : public CVariable {
protected:
  su2double  Velocity2;      /*!< \brief Square of the velocity vector. */
  su2double *HB_Source;     /*!< \brief harmonic balance source term. */
  su2double  Precond_Beta;  /*!< \brief Low Mach number preconditioner value, Beta. */
  su2double *WindGust;      /*! < \brief Wind gust value */
  su2double *WindGustDer;   /*! < \brief Wind gust derivatives value */

  /*--- Primitive variable definition ---*/

  su2double *Primitive;  /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  su2double **Gradient_Primitive;  /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */

  /*--- Secondary variable definition ---*/

  su2double *Secondary;            /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  su2double **Gradient_Secondary;  /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double *Limiter_Secondary;   /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */

  /*--- New solution container for Classical RK4 ---*/

  su2double *Solution_New;

  /*--- Old solution container for BGS iterations ---*/
  su2double* Solution_BGS_k;

public:

  /*!
   * \brief Constructor of the class.
   */
  CEulerVariable(void);

  /*!
   * \overload
   * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_energy - Value of the flow energy (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CEulerVariable(su2double val_density, su2double *val_velocity, su2double val_energy, unsigned short val_nDim,
                 unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CEulerVariable(void);

  /*!
   * \brief Get the new solution of the problem (Classical RK4).
   * \param[in] val_var - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_New(unsigned short val_var) {return Solution_New[val_var]; }

  /*!
   * \brief Set the new solution container for Classical RK4.
   */
  inline void SetSolution_New(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_New[iVar] = Solution[iVar];
  }

  /*!
   * \brief Add a value to the new solution container for Classical RK4.
   * \param[in] val_var - Number of the variable.
   * \param[in] val_solution - Value that we want to add to the solution.
   */
  inline void AddSolution_New(unsigned short val_var, su2double val_solution) {Solution_New[val_var] += val_solution;}

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
   * \brief Set to zero the gradient of the primitive variables.
   */
  void SetGradient_SecondaryZero(unsigned short val_secondaryvar);

  /*!
   * \brief Add <i>val_value</i> to the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to add to the gradient of the primitive variables.
   */
  inline void AddGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Secondary[val_var][val_dim] += val_value; }

  /*!
   * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
   */
  inline void SubtractGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Secondary[val_var][val_dim] -= val_value; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Secondary(unsigned short val_var, unsigned short val_dim) {return Gradient_Secondary[val_var][val_dim]; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Secondary(unsigned short val_var) {return Limiter_Secondary[val_var]; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGradient_Secondary(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Secondary[val_var][val_dim] = val_value; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetLimiter_Secondary(unsigned short val_var, su2double val_value) {Limiter_Secondary[val_var] = val_value; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double **GetGradient_Secondary(void) {return Gradient_Secondary; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double *GetLimiter_Secondary(void) {return Limiter_Secondary; }

  /*!
   * \brief A virtual member.
   */
  inline void SetdPdrho_e(su2double dPdrho_e) {Secondary[0] = dPdrho_e;}

  /*!
   * \brief A virtual member.
   */
  inline void SetdPde_rho(su2double dPde_rho) {Secondary[1] = dPde_rho;}

  /*!
   * \brief Set the value of the pressure.
   */
  inline bool SetPressure(su2double pressure) {
    Primitive[nDim+1] = pressure;
    if (Primitive[nDim+1] > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the speed of the sound.
   * \param[in] soundspeed2 - Value of soundspeed^2.
   */
  bool SetSoundSpeed(su2double soundspeed2) {
    su2double radical = soundspeed2;
    if (radical < 0.0) return true;
    else {
      Primitive[nDim+4] = sqrt(radical);
      return false;
    }
  }

  /*!
   * \brief Set the value of the enthalpy.
   */
  inline void SetEnthalpy(void) {Primitive[nDim+3] = (Solution[nVar-1] + Primitive[nDim+1]) / Solution[0]; }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(CFluidModel *FluidModel);

  /*!
   * \brief A virtual member.
   */
  void SetSecondaryVar(CFluidModel *FluidModel);

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
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetSecondary(unsigned short val_var) {return Secondary[val_var]; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetSecondary(unsigned short val_var, su2double val_secondary) {Secondary[val_var] = val_secondary; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetSecondary(su2double *val_secondary) {
    for (unsigned short iVar = 0; iVar < nSecondaryVar; iVar++)
      Secondary[iVar] = val_secondary[iVar];
  }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetSecondary(void) {return Secondary; }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   */
  inline bool SetDensity(void) {
    Primitive[nDim+2] = Solution[0];
    if (Primitive[nDim+2] > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the temperature.
   * \param[in] temperature - how agitated the particles are :)
   */
  inline bool SetTemperature(su2double temperature) {
    Primitive[0] = temperature;
    if (Primitive[0] > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(void) {return Velocity2; }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(void) {return Primitive[nDim+1]; }

  /*!
   * \brief Get the speed of the sound.
   * \return Value of speed of the sound.
   */
  inline su2double GetSoundSpeed(void) {return Primitive[nDim+4]; }

  /*!
   * \brief Get the enthalpy of the flow.
   * \return Value of the enthalpy of the flow.
   */
  inline su2double GetEnthalpy(void) {return Primitive[nDim+3]; }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(void) {return Solution[0]; }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(void) {return Solution[nVar-1]/Solution[0]; };

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(void) {return Primitive[0]; }

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
   * \brief Set the velocity vector from the solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity(void) {
    Velocity2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Primitive[iDim+1] = Solution[iDim+1] / Solution[0];
      Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
    }
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(su2double *val_velocity) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Solution_Old[iDim+1] = val_velocity[iDim]*Solution[0];
  }

  /*!
   * \brief Set the harmonic balance source term.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_solution - Value of the harmonic balance source term. for the index <i>val_var</i>.
   */
  inline void SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) {HB_Source[val_var] = val_source; }

  /*!
   * \brief Get the harmonic balance source term.
   * \param[in] val_var - Index of the variable.
   * \return Value of the harmonic balance source term for the index <i>val_var</i>.
   */
  inline su2double GetHarmonicBalance_Source(unsigned short val_var) {return HB_Source[val_var]; }

  /*!
   * \brief Get the value of the preconditioner Beta.
   * \return Value of the low Mach preconditioner variable Beta
   */
  inline su2double GetPreconditioner_Beta() {return Precond_Beta; }

  /*!
   * \brief Set the value of the preconditioner Beta.
   * \param[in] Value of the low Mach preconditioner variable Beta
   */
  inline void SetPreconditioner_Beta(su2double val_Beta) {Precond_Beta = val_Beta; }

  /*!
   * \brief Get the value of the wind gust
   * \return Value of the wind gust
   */
  inline su2double* GetWindGust() {return WindGust;}

  /*!
   * \brief Set the value of the wind gust
   * \param[in] Value of the wind gust
   */
  inline void SetWindGust(su2double* val_WindGust) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      WindGust[iDim] = val_WindGust[iDim];
  }

  /*!
   * \brief Get the value of the derivatives of the wind gust
   * \return Value of the derivatives of the wind gust
   */
  inline su2double* GetWindGustDer() {return WindGustDer;}

  /*!
   * \brief Set the value of the derivatives of the wind gust
   * \param[in] Value of the derivatives of the wind gust
   */
  inline void SetWindGustDer(su2double* val_WindGustDer) {
    for (unsigned short iDim = 0; iDim < nDim+1; iDim++)
      WindGustDer[iDim] = val_WindGustDer[iDim];
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
