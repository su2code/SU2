/*!
 * \file CEulerVariable.hpp
 * \brief Class for defining the variables of the compressible Euler solver.
 * \author F. Palacios, T. Economon, W. Maier, S.R. Copeland
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
 * \class CTNE2EulerVariable
 * \brief Main class for defining the variables of the TNE2 Euler's solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 6.2.0
 */
class CTNE2EulerVariable : public CVariable {
protected:
  bool ionization;          /*!< \brief Presence of charged species in gas mixture. */
  unsigned short nSpecies;  /*!< \brief Number of species in the gas mixture. */
  su2double Velocity2;		  /*!< \brief Square of the velocity vector. */
  su2double Precond_Beta;	  /*!< \brief Low Mach number preconditioner value, Beta. */

  /*--- Primitive variable definition ---*/

  su2double *Primitive;	            /*!< \brief Primitive variables (T,vx,vy,vz,P,rho,h,c) in compressible flows. */
  su2double **Gradient_Primitive;	  /*!< \brief Gradient of the primitive variables (T,vx,vy,vz,P,rho). */
  su2double *Limiter_Primitive;     /*!< \brief Limiter of the primitive variables (T,vx,vy,vz,P,rho). */

  /*--- Secondary variable definition ---*/

  su2double *Secondary;             /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  su2double **Gradient_Secondary;   /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double *Limiter_Secondary;     /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */

  /*--- New solution container for Classical RK4 ---*/

  su2double *Solution_New;

  /*--- Other Necessary Variable Definition ---*/

  su2double *dPdU;   /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  su2double *dTdU;   /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */
  su2double *dTvedU; /*!< \brief Partial derivative of vib.-el. temperature w.r.t. conserved variables. */
  su2double *eves;   /*!< \brief energy of vib-el mode w.r.t. species. */
  su2double *Cvves;  /*!< \brief Specific heat of vib-el mode w.r.t. species. */

  /*--- Index Definition ---*/

  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;

public:

  /*!
   * \brief Constructor of the class.
   */
  CTNE2EulerVariable(void);

  /*!
   * \overload
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive vars of the problem
   * \param[in] val_nPrimVarGrad - Number of prim vars gradients of the problem
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2EulerVariable(unsigned short val_nDim, unsigned short val_nVar,
                     unsigned short val_nPrimVar,
                     unsigned short val_nPrimVarGrad,
                     CConfig *config);

  /*!
   * \overload
   * \param[in] val_pressure - Value of the flow pressure (initialization value).
   * \param[in] val_massfrac - Value of the mass fraction (initialization value).
   * \param[in] val_mach - Value of the Mach number (initialization value).
   * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_temperature_ve - Value of the flow temperature_ve (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nVarPrim - Number of primitive variables.
   * \param[in] val_nVarPrimGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2EulerVariable(su2double val_pressure, su2double *val_massfrac,
                     su2double *val_mach, su2double val_temperature,
                     su2double val_temperature_ve, unsigned short val_nDim,
                     unsigned short val_nVar, unsigned short val_nVarPrim,
                     unsigned short val_nVarPrimGrad, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_nPrimVar - Number of primitive vars of the problem
   * \param[in] val_nPrimVarGrad - Number of prim vars gradients of the problem
   * \param[in] config - Definition of the particular problem.
   */
  CTNE2EulerVariable(su2double *val_solution, unsigned short val_nDim,
                     unsigned short val_nVar, unsigned short val_nVarPrim,
                     unsigned short val_nVarPrimGrad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTNE2EulerVariable(void);

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
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) {Gradient_Primitive[val_var][val_dim] = val_value; }

  /*!
    * \brief Get the value of the primitive variables gradient.
    * \return Value of the primitive variables gradient.
    */
  inline su2double **GetGradient_Primitive(void) {return Gradient_Primitive; }

  /*!
   * \brief Set the value of the velocity*velocity.
   */
  void SetVelocity2(void);

  /*!
   * \brief Set the value of the mixture density.
   */
  bool SetDensity(void);

  /*!
   * \brief Set the value of the pressure.  Requires T&Tve calculation.
   */
  bool SetPressure(CConfig *config);

  /*!
   * \brief Set the value of the speed of the sound.
   * \param[in] config
   */
  bool SetSoundSpeed(CConfig *config);

  /*!
   * \brief Set the value of the enthalpy.
   */
  inline void SetEnthalpy(void) { Primitive[H_INDEX] = (Solution[nSpecies+nDim] + Primitive[P_INDEX]) / Primitive[RHO_INDEX]; }

  /*!
   * \brief Sets gas mixture quantities (\f$\rho C^{trans-rot}_v\f$ & \f$\rho C^{vib-el}_v\f$)
   */
  void SetGasProperties(CConfig *config);

  /*!
   * \brief Calculates vib.-el. energy per mass, \f$e^{vib-el}_s\f$, for input species (not including KE)
   */
  su2double CalcEve(CConfig *config, su2double val_Tve, unsigned short val_Species);

  /*!
   * \brief Returns the stored value of Eve at the specified node
   */
  inline su2double *GetEve(void) { return eves; }

  /*!
   * \brief Calculates enthalpy per mass, \f$h^{vib-el}_s\f$, for input species (not including KE)
   */
  su2double CalcHs(CConfig *config, su2double val_T, su2double val_eves,
                   unsigned short val_Species);

  /*!
   * \brief Calculates enthalpy per mass, \f$C^{vib-el}_{v_s}\f$, for input species (not including KE)
   */
  su2double CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species);

  /*!
   * \brief Returns the value of Cvve at the specified node
   */
  su2double *GetCvve(void) { return Cvves; }

  /*!
   * \brief Calculates partial derivative of pressure w.r.t. conserved variables \f$\frac{\partial P}{\partial U}\f$
   * \param[in] config - Configuration settings
   * \param[in] dPdU - Passed-by-reference array to assign the derivatives
   */
  void CalcdPdU(su2double *V, su2double *val_eves, CConfig *config, su2double *dPdU);

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void CalcdTdU(su2double *V, CConfig *config, su2double *dTdU);

  /*!
   * \brief Set partial derivative of vib.-el. temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void CalcdTvedU(su2double *V, su2double *val_eves, CConfig *config, su2double *dTvedU);

  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void SetdTdrhos(CConfig *config);

  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  void SetdTvedrhos(CConfig *config);

  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  inline su2double *GetdPdU(void) { return dPdU; }

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
   */
  inline su2double *GetdTdU(void) { return dTdU; }

  /*!
   * \brief Set partial derivative of vib.-el. temperature w.r.t. density \f$\frac{\partial T^{V-E}}{\partial \rho_s}\f$
   */
  inline su2double *GetdTvedU(void) { return dTvedU; }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar_Compressible(CConfig *config);

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  void SetPrimVar_Gradient(CConfig *config);

  /*!
   * \brief Set all the conserved variables.
   */
  bool Cons2PrimVar(CConfig *config, su2double *U, su2double *V, su2double *dPdU,
                    su2double *dTdU, su2double *dTvedU, su2double *val_eves,
                    su2double *val_Cvves);

  /*!
   * \brief Set Gradient of the primitive variables from
   */
  bool GradCons2GradPrimVar(CConfig *config, su2double *U, su2double *V,
                            su2double **GradU, su2double **GradV);

  /*!
   * \brief Set all the conserved variables.
   */
  void Prim2ConsVar(CConfig *config, su2double *V, su2double *U);

  /*!
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetPrimVar(unsigned short val_var) {return Primitive[val_var]; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimVar(unsigned short val_var, su2double val_prim) {Primitive[val_var] = val_prim; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
   inline void SetPrimVar(su2double *val_prim) {
     for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
       Primitive[iVar] = val_prim[iVar];
   };

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimVar(void) { return Primitive; }

  /*!
   * \brief A virtual member.
   * \param[in] config - Configuration parameters.
   */
  bool SetTemperature(CConfig *config);

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(void) {return Velocity2; }
  
  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(void) { return Primitive[P_INDEX]; }

  /*!
   * \brief Get the speed of the sound.
   * \return Value of speed of the sound.
   */
  inline su2double GetSoundSpeed(void) { return Primitive[A_INDEX]; }

  /*!
   * \brief Get the enthalpy of the flow.
   * \return Value of the enthalpy of the flow.
   */
  inline su2double GetEnthalpy(void) { return Primitive[H_INDEX]; }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(void) { return Primitive[RHO_INDEX]; }

  /*!
   * \brief Get the mass fraction \f$\rho_s / \rho \f$ of species s.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline su2double GetMassFraction(unsigned short val_Species)  {
    return Primitive[RHOS_INDEX+val_Species] / Primitive[RHO_INDEX];
  }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(void) { return Solution[nSpecies+nDim]/Primitive[RHO_INDEX]; }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(void) { return Primitive[T_INDEX]; }

  /*!
   * \brief Sets the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline bool SetTemperature(su2double val_T) { Primitive[T_INDEX] = val_T; return false; }

  /*!
   * \brief A virtual member.
   * \return Value of the vibrational-electronic temperature.
   */
  inline su2double GetTemperature_ve(void) { return Primitive[TVE_INDEX]; }

  /*!
   * \brief Sets the vibrational electronic temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline bool SetTemperature_ve(su2double val_Tve) { Primitive[TVE_INDEX] = val_Tve; return false; }

  /*!
   * \brief Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  inline su2double GetRhoCv_tr(void) { return Primitive[RHOCVTR_INDEX]; }

  /*!
   * \brief Get the mixture specific heat at constant volume (vib.-el.).
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  inline su2double GetRhoCv_ve(void) { return Primitive[RHOCVVE_INDEX]; }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline su2double GetVelocity(unsigned short val_dim) { return Primitive[VEL_INDEX+val_dim]; }

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
     Solution_Old[nSpecies+iDim] = val_velocity[iDim]*Primitive[RHO_INDEX];
  }

  /*!
   * \brief Get the value of the limiter.
   */
  inline su2double *GetLimiter_Primitive(void) {return Limiter_Primitive; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned short val_var) {return Limiter_Primitive[val_var]; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetLimiter_Primitive(unsigned short val_var, su2double val_value) {Limiter_Primitive[val_var] = val_value; }

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
   * \brief Retrieves the value of the species density in the primitive variable vector.
   */
  inline unsigned short GetRhosIndex(void) { return RHOS_INDEX; }

  /*!
   * \brief Retrieves the value of the total density in the primitive variable vector.
   */
  inline unsigned short GetRhoIndex(void) { return RHO_INDEX; }

  /*!
   * \brief Retrieves the value of the pressure in the primitive variable vector.
   */
  inline unsigned short GetPIndex(void) { return P_INDEX; }

  /*!
   * \brief Retrieves the value of the in temperature the primitive variable vector.
   */
  inline unsigned short GetTIndex(void) { return T_INDEX; }

  /*!
   * \brief Retrieves the value of the vibe-elec temperature in the primitive variable vector.
   */
  inline unsigned short GetTveIndex(void) { return TVE_INDEX; }

  /*!
   * \brief Retrieves the value of the velocity  in the primitive variable vector.
   */
  inline unsigned short GetVelIndex(void) { return VEL_INDEX; }

  /*!
   * \brief Retrieves the value of the enthalpy in the primitive variable vector.
   */
  inline unsigned short GetHIndex(void) { return H_INDEX; }

  /*!
   * \brief Retrieves the value of the soundspeed in the primitive variable vector.
   */
  inline unsigned short GetAIndex(void) { return A_INDEX; }

  /*!
   * \brief Retrieves the value of the RhoCvtr in the primitive variable vector.
   */
  inline unsigned short GetRhoCvtrIndex(void) { return RHOCVTR_INDEX; }

  /*!
   * \brief Retrieves the value of the RhoCvve in the primitive variable vector.
   */
  inline unsigned short GetRhoCvveIndex(void) { return RHOCVVE_INDEX; }

};
