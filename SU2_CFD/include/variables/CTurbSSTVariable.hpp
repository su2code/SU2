/*!
 * \file CTurbSSTVariable.hpp
 * \brief Declaration of the variables of the SST turbulence model.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "CTurbVariable.hpp"

/*!
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSSTVariable final : public CTurbVariable {
protected:
  su2double sigma_om2;
  su2double beta_star;
  VectorType F1;
  VectorType F2;    /*!< \brief Menter blending function for blending of k-w and k-eps. */
  VectorType CDkw;  /*!< \brief Cross-diffusion. */

  MatrixType Primitive; /*!< \brief Primitive form of the solution. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] kine - Turbulence kinetic energy (k) (initialization value).
   * \param[in] omega - Turbulent variable value (initialization value).
   * \param[in] mut - Eddy viscosity (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants -
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSSTVariable(su2double kine, su2double omega, su2double mut, unsigned long npoint,
                   unsigned long ndim, unsigned long nvar, const su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSSTVariable() = default;

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   */
  void SetBlendingFunc(unsigned long iPoint, su2double val_viscosity, su2double val_dist, su2double val_density) override;

  /*!
   * \brief Get the first blending function.
   */
  inline su2double GetF1blending(unsigned long iPoint) const override { return F1(iPoint); }

  /*!
   * \brief Get the second blending function.
   */
  inline su2double GetF2blending(unsigned long iPoint) const override { return F2(iPoint); }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline su2double GetCrossDiff(unsigned long iPoint) const override { return CDkw(iPoint); }

  /*!
   * \brief Get the primitive variables for all points.
   * \return Reference to primitives.
   */
  inline const MatrixType& GetPrimitive(void) const { return Primitive; }

  /*!
   * \brief Get the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \return Value of the primitive variable for the index <i>iVar</i>.
   */
  inline su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const final { return Primitive(iPoint,iVar); }

  /*!
   * \brief Set the value of the primitive variables.
   * \param[in] iVar - Index of the variable.
   * \param[in] iVar - Index of the variable.
   * \return Set the value of the primitive variable for the index <i>iVar</i>.
   */
  inline void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) final { Primitive(iPoint,iVar) = val_prim; }

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
   * \brief Get the primitive variables of the problem.
   * \return Pointer to the primitive variable vector.
   */
  inline su2double *GetPrimitive(unsigned long iPoint) final {return Primitive[iPoint]; }

  /*!
   * \brief Register the variables in the solution array as input/output variable.
   * \param[in] input - input or output variables.
   * \param[in] push_index - boolean whether we want to push the index or save it in a member variable.
   */
  void RegisterConservativeSolution(bool input, bool push_index = true);

  /*!
   * \brief Register the variables in the solution_time_n array as input/output variable.
   */
  void RegisterConservativeSolution_time_n();

  /*!
   * \brief Register the variables in the solution_time_n1 array as input/output variable.
   */
  void RegisterConservativeSolution_time_n1();

  /*!
   * \brief Add a value to the solution.
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Number of the variable.
   * \param[in] solution - Value that we want to add to the solution.
   */
  inline void AddConservative(unsigned long iPoint, unsigned long iVar, su2double solution,
                              su2double lowerlimit, su2double upperlimit) { 
    // su2double cons_new = Conservative(iPoint, iVar) + solution;
    // Conservative(iPoint,iVar) = min(max(cons_new, lowerlimit), upperlimit);
  }
};
