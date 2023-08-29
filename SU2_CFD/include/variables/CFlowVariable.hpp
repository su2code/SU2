/*!
 * \file CFlowVariable.hpp
 * \brief Class for defining the common variables of flow solvers.
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CFlowVariable
 * \brief Class for defining the common variables of flow solvers.
 */
class CFlowVariable : public CVariable {
 protected:
  /*--- Primitive variable definition. ---*/
  MatrixType Primitive;                     /*!< \brief Primitive variables. */
  CVectorOfMatrix Gradient_Primitive;       /*!< \brief Gradient of the primitive variables. */
  CVectorOfMatrix& Gradient_Reconstruction; /*!< \brief Reference to the gradient of the primitive variables for MUSCL
                                               reconstruction for the convective term */
  CVectorOfMatrix
      Gradient_Aux; /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */
  MatrixType Limiter_Primitive; /*!< \brief Limiter of the primitive variables. */
  VectorType Velocity2;         /*!< \brief Squared norm of velocity. */

  MatrixType Solution_New; /*!< \brief New solution container for Classical RK4. */
  MatrixType HB_Source;    /*!< \brief harmonic balance source term. */

  /*--- NS Variables declared here to make it easier to re-use code between compressible and incompressible solvers.
   * ---*/
  MatrixType Vorticity; /*!< \brief Vorticity of the flow field. */
  VectorType StrainMag; /*!< \brief Magnitude of rate of strain tensor. */

  /*!
   * \brief Constructor of the class.
   * \note This class is not meant to be instantiated directly, it is only a building block.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions.
   * \param[in] nvar - Number of variables.
   * \param[in] nprimvar - Number of primitive variables.
   * \param[in] nprimvargrad - Number of primitive variables for which to compute gradients.
   * \param[in] config - Definition of the particular problem.
   */
  CFlowVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, unsigned long nprimvar,
                unsigned long nprimvargrad, const CConfig* config);

 public:
  mutable su2vector<int8_t> NonPhysicalEdgeCounter;  /*!< \brief Non-physical reconstruction counter for each edge. */
  /*!
   * \brief Updates the non-physical counter of an edge.
   * \param[in] iEdge - Edge index.
   * \param[in] isNonPhys - Should be 1 (true) if a non-physical reconstruction was occurred.
   * \return Whether the reconstruction should be limited to first order, based on the counter.
   */
  template <class T>
  inline T UpdateNonPhysicalEdgeCounter(unsigned long iEdge, const T& isNonPhys) const {
    if (isNonPhys != 0) {
      /*--- Force 1st order for this edge for at least 20 iterations. ---*/
      NonPhysicalEdgeCounter[iEdge] = 21;
    }
    NonPhysicalEdgeCounter[iEdge] = std::max<int8_t>(0, NonPhysicalEdgeCounter[iEdge] - 1);
    return static_cast<T>(NonPhysicalEdgeCounter[iEdge] > 0);
  }

  /*!
   * \brief Get a primitive variable.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final {
    return Primitive(iPoint, iVar);
  }

  /*!
   * \brief Get all the primitive variables of the problem.
   * \param[in] iPoint - Point index.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double* GetPrimitive(unsigned long iPoint) final { return Primitive[iPoint]; }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive() const final { return Primitive; }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_prim - Value of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final {
    Primitive(iPoint, iVar) = val_prim;
  }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iPoint - Point index.
   * \param[in] val_prim - Values of primitive variables.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, const su2double* val_prim) final {
    for (unsigned long iVar = 0; iVar < nPrimVar; iVar++) Primitive(iPoint, iVar) = val_prim[iVar];
  }

  /*!
   * \brief Get the squared norm of the velocity.
   * \param[in] iPoint - Point index.
   * \return Squared norm of the velocity vector.
   */
  inline su2double GetVelocity2(unsigned long iPoint) const final { return Velocity2(iPoint); }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetGradient_Primitive(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Primitive(iPoint, iVar, iDim);
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVarStart - Offset from which to start reading.
   * \return Value of the primitive variables gradient.
   */
  inline CMatrixView<su2double> GetGradient_Primitive(unsigned long iPoint, unsigned long iVarStart = 0) final {
    return Gradient_Primitive(iPoint, iVarStart);
  }

  /*!
   * \brief Get the primitive variable gradients for all points.
   * \return Reference to primitive variable gradient.
   */
  inline CVectorOfMatrix& GetGradient_Primitive() final { return Gradient_Primitive; }
  inline const CVectorOfMatrix& GetGradient_Primitive() const final { return Gradient_Primitive; }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Point index.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final {
    return Gradient_Reconstruction[iPoint];
  }

  /*!
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction() final { return Gradient_Reconstruction; }
  inline const CVectorOfMatrix& GetGradient_Reconstruction() const final { return Gradient_Reconstruction; }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variables gradient.
   */
  inline su2double GetLimiter_Primitive(unsigned long iPoint, unsigned long iVar) const final {
    return Limiter_Primitive(iPoint, iVar);
  }

  /*!
   * \brief Get the value of the primitive variables gradient.
   * \param[in] iPoint - Point index.
   * \return Value of the primitive variables gradient.
   */
  inline su2double* GetLimiter_Primitive(unsigned long iPoint) final { return Limiter_Primitive[iPoint]; }

  /*!
   * \brief Get the primitive variables limiter.
   * \return Primitive variables limiter for the entire domain.
   */
  inline MatrixType& GetLimiter_Primitive() final { return Limiter_Primitive; }
  inline const MatrixType& GetLimiter_Primitive() const final { return Limiter_Primitive; }

  /*!
   * \brief Get the new solution of the problem (Classical RK4).
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_New(unsigned long iPoint, unsigned long iVar) const final {
    return Solution_New(iPoint, iVar);
  }

  /*!
   * \brief Add a value to the new solution container for Classical RK4.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_solution - Value that we want to add to the new solution.
   */
  inline void AddSolution_New(unsigned long iPoint, unsigned long iVar, su2double val_solution) final {
    Solution_New(iPoint, iVar) += val_solution;
  }

  /*!
   * \brief Set the new solution container for Classical RK4.
   */
  void SetSolution_New() final;

  /*!
   * \brief Get the harmonic balance source term.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \return Value of the harmonic balance source term for the index <i>iVar</i>.
   */
  inline su2double GetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar) const final {
    return HB_Source(iPoint, iVar);
  }

  /*!
   * \brief Set the harmonic balance source term.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_solution - Value of the harmonic balance source term. for the index <i>iVar</i>.
   */
  inline void SetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar, su2double val_source) final {
    HB_Source(iPoint, iVar) = val_source;
  }

  /*!
   * \brief Get the values of the vorticity (3 values also in 2D).
   * \param[in] iPoint - Point index.
   * \return Vorticity array.
   */
  inline su2double* GetVorticity(unsigned long iPoint) final { return Vorticity[iPoint]; }
  inline const su2double* GetVorticity(unsigned long iPoint) const final { return Vorticity[iPoint]; }

  /*!
   * \brief Get the magnitude of rate of strain.
   * \param[in] iPoint - Point index.
   * \return Value of magnitude.
   */
  inline su2double GetStrainMag(unsigned long iPoint) const final { return StrainMag(iPoint); }

  /*!
   * \brief Get the entire vector of the rate of strain magnitude.
   * \return Vector of magnitudes.
   */
  inline su2activevector& GetStrainMag() { return StrainMag; }
};
