/*!
 * \file CNEMOEulerVariable.hpp
 * \brief Class for defining the variables of the compressible NEMO Euler solver.
 * \author C. Garbacz, W. Maier, S.R. Copeland
 * \version 7.1.0 "Blackbird"
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

#include "CVariable.hpp"
#include "../fluid/CNEMOGas.hpp"

/*!
 * \class CNEMOEulerVariable
 * \brief Main class for defining the variables of the NEMO Euler's solver.
 * \ingroup Euler_Equations
 * \author S. R. Copeland, F. Palacios, W. Maier, C. Garbacz
 * \version 7.1.0
 */
class CNEMOEulerVariable : public CVariable {
public:
  static constexpr size_t MAXNVAR = 25;

protected:

  bool ionization;          /*!< \brief Presence of charged species in gas mixture. */
  bool monoatomic = false;  /*!< \brief Presence of single species gas. */

  VectorType Velocity2;     /*!< \brief Square of the velocity vector. */
  MatrixType Precond_Beta;  /*!< \brief Low Mach number preconditioner value, Beta. */

  CVectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the conservative variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix  Gradient_Aux;             /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

  /*--- Primitive variable definition ---*/
  MatrixType Primitive;                /*!< \brief Primitive variables (rhos_s, T, Tve, ...) in compressible flows. */
  MatrixType Primitive_Aux;            /*!< \brief Primitive auxiliary variables (Y_s, T, Tve, ...) in compressible flows. */
  CVectorOfMatrix Gradient_Primitive;  /*!< \brief Gradient of the primitive variables (rhos_s, T, Tve, ...). */
  MatrixType Limiter_Primitive;        /*!< \brief Limiter of the primitive variables  (rhos_s, T, Tve, ...). */

  /*--- Secondary variable definition ---*/
  MatrixType Secondary;                /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  CVectorOfMatrix Gradient_Secondary;  /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */

  /*--- New solution container for Classical RK4 ---*/
  MatrixType Solution_New;  /*!< \brief New solution container for Classical RK4. */

  /*--- Other Necessary Variable Definition ---*/
  MatrixType dPdU;   /*!< \brief Partial derivative of pressure w.r.t. conserved variables. */
  MatrixType dTdU;   /*!< \brief Partial derivative of temperature w.r.t. conserved variables. */
  MatrixType dTvedU; /*!< \brief Partial derivative of vib.-el. temperature w.r.t. conserved variables. */
  MatrixType eves;   /*!< \brief energy of vib-el mode w.r.t. species. */
  MatrixType Cvves;  /*!< \brief Specific heat of vib-el mode w.r.t. species. */
  VectorType Gamma;  /*!< \brief Ratio of specific heats. */

  CNEMOGas *fluidmodel;

  /*!< \brief Index definition for NEMO pritimive variables. */
  unsigned long RHOS_INDEX, T_INDEX, TVE_INDEX, VEL_INDEX, P_INDEX,
  RHO_INDEX, H_INDEX, A_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX,
  LAM_VISC_INDEX, EDDY_VISC_INDEX, nSpecies;

  su2double Tve_Freestream; /*!< \brief Freestream vib-el temperature. */

  VectorType symmetry;      /*!< \brief Nodes in symmetry planes. */

  VectorType x_new_freestream;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_pressure - Value of the flow pressure (initialization value).
   * \param[in] val_massfrac - Value of the mass fraction (initialization value).
   * \param[in] val_mach - Value of the Mach number (initialization value).
   * \param[in] val_temperature - Value of the flow temperature (initialization value).
   * \param[in] val_temperature_ve - Value of the flow temperature_ve (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of conserved variables.
   * \param[in] val_nVarPrim - Number of primitive variables.
   * \param[in] val_nVarPrimGrad - Number of primitive gradient variables.
   * \param[in] config - Definition of the particular problem.
   */
  CNEMOEulerVariable(su2double val_pressure, const su2double *val_massfrac,
                     su2double *val_mach, su2double val_temperature,
                     su2double val_temperature_ve, unsigned long npoint,
                     unsigned long ndim,
                     unsigned long nvar, unsigned long nvalprim,
                     unsigned long nvarprimgrad, CConfig *config, CNEMOGas *fluidmodel);

  /*!
   * \brief Destructor of the class.
   */
  ~CNEMOEulerVariable() override = default;

  /*---------------------------------------*/
  /*---          U,V,S Routines         ---*/
  /*---------------------------------------*/

  /*!
   * \brief Get the new solution of the problem (Classical RK4).
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_New(unsigned long iPoint, unsigned long iVar) const final { return Solution_New(iPoint,iVar); }

  /*!
   * \brief Set the new solution container for Classical RK4.
   */
  void SetSolution_New() final;

  /*!
   * \brief Add a value to the new solution container for Classical RK4.
   * \param[in] iVar - Number of the variable.
   * \param[in] val_solution - Value that we want to add to the solution.
   */
  inline void AddSolution_New(unsigned long iPoint, unsigned long iVar, su2double val_solution) final {
    Solution_New(iPoint,iVar) += val_solution;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) override { Primitive(iPoint,iVar) = val_prim; }

inline bool Set_x_new_freestream(unsigned long iPoint, su2double val_x_new_freestream) final
                               { x_new_freestream(iPoint) = val_x_new_freestream; return false; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) final {
    for (unsigned long iVar = 0; iVar < nPrimVar; iVar++)
      Primitive(iPoint,iVar) = val_prim[iVar];
  }

  /*!
   * \brief Get the primitive variables limiter.
   * \return Primitive variables limiter for the entire domain.
   */
  inline MatrixType& GetLimiter_Primitive(void) {return Limiter_Primitive; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const final {return Limiter_Primitive(iPoint,iVar); }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double *GetLimiter_Primitive(unsigned long iPoint) final { return Limiter_Primitive[iPoint]; }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] value - Value of the gradient.
   */
  inline void SetLimiter_Primitive(unsigned long iPoint, unsigned long iVar, su2double value) final {
    Limiter_Primitive(iPoint,iVar) = value;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetSecondary(unsigned long iPoint, unsigned long iVar, su2double val_secondary) final {Secondary(iPoint,iVar) = val_secondary; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetSecondary(unsigned long iPoint, const su2double *val_secondary) final {
    for (unsigned long iVar = 0; iVar < nSecondaryVar; iVar++)
      Secondary(iPoint,iVar) = val_secondary[iVar];
  }

  /*!
   * \brief Set the value of the primitive auxiliary variables - with mass fractions.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive_Aux(unsigned long iPoint, unsigned long iVar, su2double val_prim) { Primitive_Aux(iPoint,iVar) = val_prim; }


  /*!
   * \brief Get the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final { return Primitive(iPoint,iVar); }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(unsigned long iPoint) final {return Primitive[iPoint]; }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive(void) const { return Primitive; }

   /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive_Aux(void) const { return Primitive_Aux; }


  /*!
   * \brief Get the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetSecondary(unsigned long iPoint, unsigned long iVar) const final {return Secondary(iPoint,iVar); }

  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetSecondary(unsigned long iPoint) final { return Secondary[iPoint]; }

  /*---------------------------------------*/
  /*---  Gradient Routines  ---*/
  /*---------------------------------------*/

 /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction(void) final { return Gradient_Reconstruction; }

  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \return Value of the reconstruction variables gradient at a node.
   */
  inline su2double GetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Reconstruction(iPoint,iVar,iDim);
  }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline su2double **GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \param[in] value  - Value of the reconstruction gradient component.
   */
  inline void SetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) override {
    Gradient_Reconstruction(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Set to zero the gradient of the primitive variables.
   */
  void SetGradient_PrimitiveZero();

  /*!
   * \brief Add <i>value</i> to the gradient of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value to add to the gradient of the primitive variables.
   */
  inline void AddGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) += value;
  }

  /*!
   * \brief Subtract <i>value</i> to the gradient of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value to subtract to the gradient of the primitive variables.
   */
  inline void SubtractGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) {
    Gradient_Primitive(iPoint,iVar,iDim) -= value;
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Primitive(iPoint,iVar,iDim);
  }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the gradient.
   */
  inline void SetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double **GetGradient_Primitive(unsigned long iPoint) final { return Gradient_Primitive[iPoint]; }

  /*!
   * \brief Get the primitive variable gradients for all points.
   * \return Reference to primitive variable gradient.
   */
  inline CVectorOfMatrix& GetGradient_Primitive(void) { return Gradient_Primitive; }

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) override;

   /*!
  * \brief Set all the primitive and secondary variables from the conserved vector.
  */
  bool Cons2PrimVar(su2double *U, su2double *V, su2double *dPdU,
                    su2double *dTdU, su2double *dTvedU, su2double *val_eves,
                    su2double *val_Cvves);

  /*---------------------------------------*/
  /*---   Specific variable routines    ---*/
  /*---------------------------------------*/

   /*!
   * \brief Set the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  void SetVelocity2(unsigned long iPoint) override;

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(unsigned long iPoint) const final { return Velocity2(iPoint); }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint,P_INDEX); }

  /*!
   * \brief Get the speed of the sound.
   * \return Value of speed of the sound.
   */
  inline su2double GetSoundSpeed(unsigned long iPoint) const final { return Primitive(iPoint,A_INDEX); }

  /*!
   * \brief Get the enthalpy of the flow.
   * \return Value of the enthalpy of the flow.
   */
  inline su2double GetEnthalpy(unsigned long iPoint) const final { return Primitive(iPoint,H_INDEX); }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint) const final { return Primitive(iPoint,RHO_INDEX); }

  /*!
   * \brief Get the specie density of the flow.
   * \return Value of the specie density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint, unsigned long val_Species) const final { return Primitive(iPoint,RHOS_INDEX+val_Species); }

  /*!
   * \brief Get the energy of the flow.
   * \return Value of the energy of the flow.
   */
  inline su2double GetEnergy(unsigned long iPoint) const final { return Solution(iPoint,nSpecies+nDim)/Primitive(iPoint,RHO_INDEX); }

  /*!
   * \brief Get the temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline su2double GetTemperature(unsigned long iPoint) const final { return Primitive(iPoint,T_INDEX); }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final { return Primitive(iPoint,VEL_INDEX+iDim); }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
    su2double ProjVel = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      ProjVel += Primitive(iPoint,VEL_INDEX+iDim)*val_vector[iDim];
    return ProjVel;
  }

  /*!
   * \brief Set the velocity vector from the solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      Primitive(iPoint,VEL_INDEX+iDim) = Solution(iPoint,nSpecies+iDim) / Primitive(iPoint,RHO_INDEX);
      Velocity2(iPoint) += pow(Primitive(iPoint,VEL_INDEX+iDim),2);
    }
  }

  /*!
   * \brief Set the Energy vector from the old solution.
   * \param[in] val_energy - Pointer to energy.
   */
  inline void SetEnergy_Old(unsigned long iPoint, const vector<su2double>& val_energy)  {
    for (unsigned long i = 0; i < 2; i++){
      Solution_Old(iPoint,nSpecies+nDim+i) = val_energy[i]*Primitive(iPoint,RHO_INDEX);
     }
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++){
      Solution_Old(iPoint,nSpecies+iDim) = val_velocity[iDim]*Primitive(iPoint,RHO_INDEX);
     }
  }

  /*!
   * \brief A virtual member.
   * \return Value of the vibrational-electronic temperature.
   */
  inline su2double GetTemperature_ve(unsigned long iPoint) const override
                                    { return Primitive(iPoint,TVE_INDEX); }

  /*!
   * \brief Sets the vibrational electronic temperature of the flow.
   * \return Value of the temperature of the flow.
   */
  inline bool SetTemperature_ve(unsigned long iPoint, su2double val_Tve) override
                               { Primitive(iPoint,TVE_INDEX) = val_Tve; return false; }

  /*!
   * \brief Get the mixture specific heat at constant volume (trans.-rot.).
   * \return \f$\rho C^{t-r}_{v} \f$
   */
  inline su2double GetRhoCv_tr(unsigned long iPoint) const override
                              { return Primitive(iPoint,RHOCVTR_INDEX); }

  /*!
   * \brief Get the mixture specific heat at constant volume (vib.-el.).
   * \return \f$\rho C^{v-e}_{v} \f$
   */
  inline su2double GetRhoCv_ve(unsigned long iPoint) const override
                              { return Primitive(iPoint,RHOCVVE_INDEX); }

  /*!
   * \brief Returns the stored value of Eve at the specified node
   */
  inline su2double *GetEve(unsigned long iPoint) { return eves[iPoint]; }

  /*!
   * \brief Returns the value of Cvve at the specified node
   */
  su2double *GetCvve(unsigned long iPoint) { return Cvves[iPoint]; }

  /*!
   * \brief Set partial derivative of pressure w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   */
  inline su2double *GetdPdU(unsigned long iPoint) override { return dPdU[iPoint]; }

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial T}{\partial \rho_s}\f$
   */
  inline su2double *GetdTdU(unsigned long iPoint) override { return dTdU[iPoint]; }

  /*!
   * \brief Set partial derivative of vib.-el. temperature w.r.t. density \f$\frac{\partial T^{V-E}}{\partial \rho_s}\f$
   */
  inline su2double *GetdTvedU(unsigned long iPoint) override { return dTvedU[iPoint]; }

  /*!
   * \brief Get the mass fraction \f$\rho_s / \rho \f$ of species s.
   * \param[in] val_Species - Index of species s.
   * \return Value of the mass fraction of species s.
   */
  inline su2double GetMassFraction(unsigned long iPoint, unsigned short val_Species) const {
    return Primitive(iPoint,RHOS_INDEX+val_Species) / Primitive(iPoint,RHO_INDEX);
  }

  /*!
   * \brief Returns the stored value of Gamma at the specified node
   */
  inline su2double GetGamma(unsigned long iPoint) { return Gamma(iPoint); }

  /*!
   * \brief Retrieves the number of symmetry planes at the specified node.
   */
  inline su2double GetSymmetry(unsigned long iPoint) { return symmetry[iPoint]; }

  /*!
   * \brief Increases the number of symmetry planes at the specified node by one.
   */
  inline void SetSymmetry(unsigned long iPoint) {symmetry[iPoint] += 1.0;}

  /*---------------------------------------*/
  /*---           NEMO indices          ---*/
  /*---------------------------------------*/

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
