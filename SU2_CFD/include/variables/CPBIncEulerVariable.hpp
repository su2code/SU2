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
#include "CIncEulerVariable.hpp"

/*!
 * \class CPBIncEulerVariable
 * \brief Class for defining the variables of the pressure based incompressible Euler solver.
 * \ingroup Euler_Equations
 * \author A. Koodly
 */
class CPBIncEulerVariable : public CFlowVariable {

public:
mutable su2vector<int8_t> NonPhysicalEdgeCounter;  /*!< \brief Non-physical reconstruction counter for each edge. */
   static constexpr size_t MAXNVAR = 7;
  template <class IndexType>
  struct CIndices {
    const IndexType nDim;
    CIndices(IndexType ndim, IndexType) : nDim(ndim) {}
    inline IndexType NDim() const { return nDim; }
    // inline IndexType NSpecies() const { return 0; }
    inline IndexType Pressure() const { return 0; }
    inline IndexType Velocity() const { return 1; }
    // inline IndexType Temperature() const { return nDim+1; }
    inline IndexType Density() const { return nDim+1; }
    // inline IndexType Beta() const { return nDim+3; }
    // inline IndexType SoundSpeed() const { return Beta(); }
    inline IndexType LaminarViscosity() const { return nDim+2; }
    inline IndexType EddyViscosity() const { return nDim+3; }
    inline IndexType Temperature() const { return nDim+4; }
    // inline IndexType ThermalConductivity() const { return nDim+4; }
    // inline IndexType CpTotal() const { return nDim+7; }
    // inline IndexType CvTotal() const { return nDim+8; }

    /*--- For compatible interface with NEMO. ---*/
    // inline IndexType SpeciesDensities() const { return std::numeric_limits<IndexType>::max(); }
    inline IndexType Temperature_ve() const { return std::numeric_limits<IndexType>::max(); }
    // inline IndexType Enthalpy() const { return std::numeric_limits<IndexType>::max(); }
    
  };
  protected:
    const CIndices<unsigned long> indices;
//   VectorType Velocity2;                    /*!< \brief Square of the velocity vector. */
//   MatrixType Primitive;                    /*!< \brief Primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) in incompressible flows. */
//   CVectorOfMatrix Gradient_Primitive;       /*!< \brief Gradient of the primitive variables (P, vx, vy, vz, T, rho, beta). */
//   CVectorOfMatrix& Gradient_Reconstruction; /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix Gradient_Aux;             /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */
  MatrixType Limiter_Primitive;            /*!< \brief Limiter of the primitive variables (P, vx, vy, vz, T, rho, beta). */
  VectorType Density_Old;                  /*!< \brief Old density for variable density turbulent flows (SST). */
  
  VectorType MassFlux;                     /*!< \brief Massflux associated with each CV */
  MatrixType Mom_Coeff;
  MatrixType Mom_Coeff_nb;
  using BoolVectorType = C2DContainer<unsigned long, bool, StorageType::ColumnMajor, 64, DynamicSize, 1>;
  BoolVectorType strong_bc;
  // VectorType StrainMag; /*!< \brief Magnitude of rate of strain tensor. */


  public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBIncEulerVariable(su2double density, su2double pressure, const su2double *velocity, 
                    unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPBIncEulerVariable() = default;
  
  /*!
   * \brief Get the primitive variable gradients for all points.
   * \return Reference to primitive variable gradient.
   */
//   inline CVectorOfMatrix& GetGradient_Primitive(void) { return Gradient_Primitive; }
  
   /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
//   inline CVectorOfMatrix& GetGradient_Reconstruction(void) final { return Gradient_Reconstruction; }

  /*!
   * \brief Add <i>value</i> to the gradient of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value to add to the gradient of the primitive variables.
   */
  inline void AddGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) += value;
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
//   inline su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
//     return Gradient_Primitive(iPoint,iVar,iDim);
//   }

  /*!
   * \brief Get the primitive variables limiter.
   * \return Primitive variables limiter for the entire domain.
   */
//   inline MatrixType& GetLimiter_Primitive(void) {return Limiter_Primitive; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
//   inline su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const final {
//     return Limiter_Primitive(iPoint,iVar);
//   }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the gradient.
   */
  inline void SetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Primitive(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Set the gradient of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] value - Value of the gradient.
   */
  inline void SetLimiter_Primitive(unsigned long iPoint, unsigned long iVar, su2double value) final {
    Limiter_Primitive(iPoint,iVar) = value;
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the primitive variables gradient.
   */
//   inline CMatrixView<su2double> GetGradient_Primitive(unsigned long iPoint, unsigned long iVar=0) final { return Gradient_Primitive[iPoint]; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the primitive variables gradient.
   */
//   inline su2double *GetLimiter_Primitive(unsigned long iPoint) final { return Limiter_Primitive[iPoint]; }
  
  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \return Value of the reconstruction variables gradient at a node.
   */
//   inline su2double GetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
//     return Gradient_Reconstruction(iPoint,iVar,iDim);
//   }
  
  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \param[in] value  - Value of the reconstruction gradient component.
   */
//   inline void SetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
//     Gradient_Reconstruction(iPoint,iVar,iDim) = value;
//   }
  
  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
//   inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
//   inline const MatrixType& GetPrimitive(void) const { return Primitive; }
  
  /*!
   * \brief Get the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
//   inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final { return Primitive(iPoint,iVar); }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
//   inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final { Primitive(iPoint,iVar) = val_prim; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] val_prim - Primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
//   inline void SetPrimitive(unsigned long iPoint, const su2double *val_prim) final {
//     for (unsigned long iVar = 0; iVar < nPrimVar; iVar++) Primitive(iPoint,iVar) = val_prim[iVar];
//   }

  /*!
   * \brief Get the primitive variables of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the primitive variable vector.
   */
//   inline su2double *GetPrimitive(unsigned long iPoint) final { return Primitive[iPoint]; }

  /*!
   * \brief Set the value of the pressure.
   */
//   inline bool SetPressure(unsigned long iPoint, su2double val_pressure) override { 
//     Primitive(iPoint,0) = val_pressure; 
//     return true;
//     }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline bool SetDensity(unsigned long iPoint, su2double val_density) final {
    Primitive(iPoint,nDim+1) = val_density;
    if (Primitive(iPoint,nDim+1) > 0.0) return false;
    else return true;
  }

  /*!
   * \brief Set the value of the density for the incompressible flows.
   * \param[in] iPoint - Point index.
   */
  inline void SetVelocity(unsigned long iPoint) final {
    Velocity2(iPoint) = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      Primitive(iPoint,iDim+1) = Solution(iPoint,iDim)/Primitive(iPoint, nDim+1);
      Velocity2(iPoint) += pow(Primitive(iPoint,iDim+1),2);
    }
  }

  /*!
   * \brief Get the norm 2 of the velocity.
   * \return Norm 2 of the velocity vector.
   */
//   inline su2double GetVelocity2(unsigned long iPoint) const final { return Velocity2(iPoint); }

  /*!
   * \brief Get the flow pressure.
   * \return Value of the flow pressure.
   */
  inline su2double GetPressure(unsigned long iPoint) const final { return Primitive(iPoint,0); }

  /*!
   * \brief Get the density of the flow.
   * \return Value of the density of the flow.
   */
  inline su2double GetDensity(unsigned long iPoint) const final { return Primitive(iPoint,nDim+1); }

  /*!
   * \brief Get the density of the flow from the previous iteration.
   * \return Old value of the density of the flow.
   */
  inline su2double GetDensity_Old(unsigned long iPoint) const final { return Density_Old(iPoint); }

  /*!
   * \brief Get the velocity of the flow.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the velocity for the dimension <i>iDim</i>.
   */
  inline su2double GetVelocity(unsigned long iPoint, unsigned long iDim) const final { return Primitive(iPoint,iDim+1); }

  /*!
   * \brief Get the projected velocity in a unitary vector direction (compressible solver).
   * \param[in] val_vector - Direction of projection.
   * \return Value of the projected velocity.
   */
  inline su2double GetProjVel(unsigned long iPoint, const su2double *val_vector) const final {
    su2double ProjVel = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      ProjVel += Primitive(iPoint,iDim+1)*val_vector[iDim];
    return ProjVel;
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim) = Primitive(iPoint,nDim+1)*val_velocity[iDim];
  }

  /*!
   * \brief Set all the primitive variables for incompressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, su2double Density_Inf,  CConfig *config);

  inline void Set_Mom_CoeffZero(unsigned long iPoint) {
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff(iPoint,iDim) = 0.0;
  }

  inline su2double Get_Mom_Coeff(unsigned long iPoint, unsigned short val_Var) { return Mom_Coeff(iPoint,val_Var);}

  inline void Set_Mom_Coeff(unsigned long iPoint, const su2double *val_Mom_Coeff) { 
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff(iPoint,iDim) = val_Mom_Coeff[iDim]; 
  }
    
  inline void Set_Mom_Coeff(unsigned long iPoint, unsigned short val_Var, su2double val_Mom_Coeff) { Mom_Coeff(iPoint,val_Var) = val_Mom_Coeff; }
    
  inline su2double Get_Mom_Coeff_nb(unsigned long iPoint, unsigned short val_Var) { return Mom_Coeff_nb(iPoint,val_Var);}

  inline void Set_Mom_Coeff_nb(unsigned long iPoint, su2double *val_Mom_Coeff) { 
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff_nb(iPoint,iDim) = val_Mom_Coeff[iDim]; 
  }
  
  inline void Set_Mom_Coeff_nb(unsigned long iPoint, unsigned short val_Var, su2double val_Mom_Coeff) { Mom_Coeff_nb(iPoint,val_Var) = val_Mom_Coeff; }
  
  inline void Set_Mom_Coeff_nbZero(unsigned long iPoint) {
	  for (unsigned short iDim = 0; iDim < nDim; iDim++)
	        Mom_Coeff_nb(iPoint,iDim) = 0.0;
  }
  
  inline void Add_Mom_Coeff_nb(unsigned long iPoint, su2double val_coeff_nb, unsigned short val_Var) { Mom_Coeff_nb(iPoint,val_Var) += val_coeff_nb;}

  inline void Add_Mom_Coeff(unsigned long iPoint, su2double val_coeff, unsigned short val_Var) { Mom_Coeff(iPoint,val_Var) += val_coeff;}
  
  inline void SetStrongBC(unsigned long iPoint) { strong_bc(iPoint) = true; }
  
  inline void ResetStrongBC(unsigned long iPoint) { strong_bc(iPoint) = false; }
  
  inline bool GetStrongBC(unsigned long iPoint) { return strong_bc(iPoint); }
  
  inline void SetMassFluxZero(unsigned iPoint) {MassFlux(iPoint) = 0.0; }
  
  inline void SetMassFlux(unsigned long iPoint, su2double val_MassFlux) { MassFlux(iPoint) = val_MassFlux; }
  
  inline void AddMassFlux(unsigned long iPoint, su2double val_MassFlux) {MassFlux(iPoint) += val_MassFlux; }
  
  inline void SubtractMassFlux(unsigned long iPoint, su2double val_MassFlux) {MassFlux(iPoint) -= val_MassFlux; }
  
  inline su2double GetMassFlux(unsigned long iPoint) const { return MassFlux(iPoint); }

    /*!
   * \brief Get the entire vector of the rate of strain magnitude.
   * \return Vector of magnitudes.
   */
//   inline su2activevector& GetStrainMag() { return StrainMag; }
  
  inline CMatrixView<const su2double> GetVelocityGradient(unsigned long iPoint) const final {
    return Gradient_Primitive(iPoint, indices.Velocity());
  }
};




/*!
 * \class CPoissonVariable
 * \brief Main class for defining the variables of the potential solver.
 * \ingroup Potential_Flow_Equation
 * \author F. Palacios
 */
class CPoissonVariable : public CVariable {
  VectorType SourceTerm;
  VectorType Poisson_Coeff;
  using BoolVectorType = C2DContainer<unsigned long, bool, StorageType::ColumnMajor, 64, DynamicSize, 1>;
  BoolVectorType strong_bc;
public:
  /*!
   * \overload
   * \param[in] val_potential - Value of the potential solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPoissonVariable(su2double val_SourceTerm, unsigned long npoint, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPoissonVariable() = default;
  
  inline su2double GetPoisson_Coeff(unsigned long iPoint) { return Poisson_Coeff(iPoint);}
  
  inline void SetPoisson_Coeff(unsigned long iPoint, su2double val_Poisson_Coeff) { Poisson_Coeff(iPoint) = val_Poisson_Coeff ; }
  
  inline void SetSourceTerm(unsigned long iPoint, su2double val_SourceTerm) { SourceTerm(iPoint) = val_SourceTerm ; }
  
  inline su2double GetSourceTerm(unsigned long iPoint) { return SourceTerm(iPoint);}
  
  inline void SetStrongBC(unsigned long iPoint) { strong_bc(iPoint) = true; }
  
  inline bool GetStrongBC(unsigned long iPoint) { return strong_bc(iPoint); }
  
  inline void ResetStrongBC(unsigned long iPoint) { strong_bc(iPoint) = false; }
  
};
