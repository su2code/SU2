/*!
 * \file CPBIncEulerVariable.hpp
 * \brief Class for defining the variables of the pressure based incompressible Euler solver.
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
 * \brief Main class for defining the variables of the incompressible Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CPBIncEulerVariable : public CVariable {
protected:
  su2double Velocity2;      /*!< \brief Square of the velocity vector. */
  su2double *VelocityOld;  /*!< \brief Velocity at prev time step. */
  su2double *WindGust;           /*! < \brief Wind gust value */
  su2double *WindGustDer;        /*! < \brief Wind gust derivatives value */
  
  /*--- Primitive variable definition ---*/
  
  su2double *Primitive;  /*!< \brief Primitive variables (T, vx, vy, vz, P, rho, h, c) in compressible flows. */
  su2double **Gradient_Primitive;  /*!< \brief Gradient of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double *Limiter_Primitive;    /*!< \brief Limiter of the primitive variables (T, vx, vy, vz, P, rho). */
  su2double MassFlux;   /*!< \brief Massflux associated with each CV */
  su2double *Mom_Coeff;
  su2double *Mom_Coeff_nb;
  su2double *PrimitiveMGCorr;
  bool strong_bc;

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CPBIncEulerVariable(void);
  
  /*!
   * \overload
   * \param[in] val_pressure - value of the pressure.
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBIncEulerVariable(su2double val_pressure, su2double *val_velocity, unsigned short val_nDim,
                    unsigned short val_nvar, CConfig *config);
  
  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBIncEulerVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPBIncEulerVariable(void);
  
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
  inline void AddGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] += val_value; }
  
  /*!
   * \brief Subtract <i>val_value</i> to the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value to subtract to the gradient of the primitive variables.
   */
  inline void SubtractGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] -= val_value; }
  
  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Primitive(unsigned short val_var, unsigned short val_dim) { return Gradient_Primitive[val_var][val_dim]; }
  
  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned short val_var) { return Limiter_Primitive[val_var]; }
  
  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_dim - Index of the dimension.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetGradient_Primitive(unsigned short val_var, unsigned short val_dim, su2double val_value) { Gradient_Primitive[val_var][val_dim] = val_value; }
  
  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_value - Value of the gradient.
   */
  inline void SetLimiter_Primitive(unsigned short val_var, su2double val_value) { Limiter_Primitive[val_var] = val_value; }
  
  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double **GetGradient_Primitive(void) { return Gradient_Primitive; }
  
  /*!
   * \brief Get the value of the primitive variables gradient.
   * \return Value of the primitive variables gradient.
   */
  inline su2double *GetLimiter_Primitive(void) { return Limiter_Primitive; }
  
  /*!
   * \brief Set the value of the pressure.
   */
  inline void SetPressure_val(su2double val_pressure) { Primitive[0] = val_pressure; }
  
  /*!
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetPrimitive(unsigned short val_var) { return Primitive[val_var]; }
  
  /*!
   * \brief Get the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \return Value of the primitive variable for the index <i>val_var</i>.
   */
  inline su2double GetPrimitiveMGCorr(unsigned short val_var) { return PrimitiveMGCorr[val_var]; }
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitive(unsigned short val_var, su2double val_prim) { Primitive[val_var] = val_prim; }
  
  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_var - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitiveMGCorr(unsigned short val_var, su2double val_prim) { PrimitiveMGCorr[val_var] = val_prim; }
  
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
   * \brief Set the value of the primitive variables.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>val_var</i>.
   */
  inline void SetPrimitiveMGCorr(su2double *val_prim_new) {
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++)
    PrimitiveMGCorr[iVar] = val_prim_new[iVar];
  }
  
  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(void) { return Primitive; }
  
  /*!
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitiveMGCorr(void) { return PrimitiveMGCorr; }
  
  /*!
   * \brief Set the value of the density for the incompressible flows.
   */
  inline bool SetDensity(su2double val_density) { Primitive[nDim+1] = val_density; return true;}
  
  /*!
   * \brief Set the value of the square of velocity for the incompressible flows.
   */
  inline void SetVelocity(void){
	  Velocity2 = 0.0;
	  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
		  Primitive[iDim+1] = Solution[iDim] ;
		  Velocity2 += Primitive[iDim+1]*Primitive[iDim+1];
	  }
  }
    
  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
  inline su2double GetVelocity2(void) { return Velocity2; }
  
  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(void) { return Primitive[0]; }
  
  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(void) { return Primitive[nDim+1]; }
  
  /*!
   * \brief Get the velocity of the flow.
   * \param[in] val_dim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>val_dim</i>.
   */
  inline su2double GetVelocity(unsigned short val_dim) { return Primitive[val_dim+1]; }
  
  /*!
   * \brief Get the projected velocity in a unitary vector direction 
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  su2double GetProjVel(su2double *val_vector);
  
  /*!
   * \brief Get the mass flux 
   * \param[in] val_density - Direction of projection.
   * \return Value of the projected velocity.
   */  
  inline void SetMassFluxZero() { MassFlux = 0.0 ; }

  /*!
   * \brief Get the mass flux 
   * \param[in] val_density - Direction of projection.
   * \return Value of the projected velocity.
   */  
  inline void AddMassFlux(su2double val_MassFlux) { MassFlux += val_MassFlux ; }
  
  /*!
   * \brief Get the mass flux 
   * \param[in] val_density - Direction of projection.
   * \return Value of the projected velocity.
   */  
  inline void SetMassFlux(su2double val_MassFlux) { MassFlux = val_MassFlux ; }
  
  /*!
   * \brief Get the mass flux 
   * \param[in] val_density - Direction of projection.
   * \return Value of the projected velocity.
   */  
  inline void SubtractMassFlux(su2double val_MassFlux) { MassFlux -= val_MassFlux ; }
  
  inline su2double GetMassFlux() {return MassFlux;}

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(su2double *val_velocity) {
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	     Solution_Old[iDim] = val_velocity[iDim];
  }
  
  /*!
   * \brief Set all the primitive variables for incompressible flows.
   */
  bool SetPrimVar(su2double Density_Inf, CConfig *config);

  inline void Set_Mom_CoeffZero() {
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff[iDim] = 0.0;
  }

  inline su2double Get_Mom_Coeff(unsigned short val_Var) { return Mom_Coeff[val_Var];}

  inline void Set_Mom_Coeff(su2double *val_Mom_Coeff) { 
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff[iDim] = val_Mom_Coeff[iDim]; 
  }
    
  inline void Set_Mom_Coeff(unsigned short val_Var, su2double val_Mom_Coeff) { Mom_Coeff[val_Var] = val_Mom_Coeff; }
    
  inline su2double Get_Mom_Coeff_nb(unsigned short val_Var) { return Mom_Coeff_nb[val_Var];}

  inline void Set_Mom_Coeff_nb(su2double *val_Mom_Coeff) { 
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff_nb[iDim] = val_Mom_Coeff[iDim]; 
  }
  
  inline void Set_Mom_Coeff_nb(unsigned short val_Var, su2double val_Mom_Coeff) { Mom_Coeff_nb[val_Var] = val_Mom_Coeff; }
  
  inline void Set_Mom_Coeff_nbZero() {
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff_nb[iDim] = 0.0;
  }
  
  inline void Add_Mom_Coeff_nb(su2double val_coeff_nb, unsigned short val_Var) { Mom_Coeff_nb[val_Var] += val_coeff_nb;}
  
  inline void SetStrongBC() { strong_bc = true; }
  
  inline void ResetStrongBC() { strong_bc = false; }
  
  inline bool GetStrongBC() { return strong_bc; }

};


/*!
 * \class CPoissonVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 */
class CPoissonVariable : public CVariable {
  su2double SourceTerm;
  su2double Poisson_Coeff;
  bool strong_bc;
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CPoissonVariable(void);
  
  /*!
   * \overload
   * \param[in] val_potential - Value of the potential solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPoissonVariable(su2double val_SourceTerm, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPoissonVariable(void);
  
  inline su2double GetPoisson_Coeff() { return Poisson_Coeff;}
  
  inline void SetPoisson_Coeff(su2double val_Poisson_Coeff) { Poisson_Coeff = val_Poisson_Coeff ; }
  
  inline void SetSourceTerm(su2double val_SourceTerm) { SourceTerm = val_SourceTerm ; }
  
  inline su2double GetSourceTerm() { return SourceTerm;}
  
  inline void SetStrongBC() { strong_bc = true; }
  
  inline bool GetStrongBC() { return strong_bc; }
  
  inline void ResetStrongBC() { strong_bc = false; }
  
};
