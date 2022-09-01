/*!
 * \file CSpeciesFlameletVariable.hpp
 * \brief Definition of the variable fields for the flamelet class.
 * \author D. Mayer, T. Economon
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

#include "CSpeciesVariable.hpp"

/*!
 * \class CSpeciesFlameletVariable
 * \brief Main class for defining the variables of the flamelet model.
 * \author D. Mayer, T. Economon
 */
class CSpeciesFlameletVariable final : public CSpeciesVariable {
protected:
  //MatrixType source_prog;            /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  MatrixType source_scalar;          /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  MatrixType lookup_scalar;          /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  //MatrixType source_scalar_rescaled; /*!< \brief Vector of the source terms from the lookup table for each scalar equation */
  //MatrixType scalar_table_paraview;  /*!< \brief Vector of the values of the scalar from the lookup table for visualization */
  
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_scalar_inf - Scalar variable value (initialization value).
   * \param[in] npoint         - Number of points/nodes/vertices in the domain.
   * \param[in] ndim           - Number of dimensions of the problem.
   * \param[in] nvar           - Number of variables of the problem.
   * \param[in] config         - Definition of the particular problem.
   */
  CSpeciesFlameletVariable(su2double     *val_scalar_inf,
                    unsigned long npoint,
                    unsigned long ndim,
                    unsigned long nvar,
                    CConfig       *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CSpeciesFlameletVariable() = default;
  
    
  /*!
   * \brief Set the value of the transported scalar source term
   * \param[in] val_- the .
   * \param[in] val_ivar        - eqn. index to the .
   */
  inline void SetSourceScalar(unsigned long  iPoint,
                              su2double      val_source_scalar,
                              unsigned short val_ivar) override {
    source_scalar(iPoint,val_ivar) = val_source_scalar;
  }
  
  /*!
   * \brief Set the value of the transported scalar source term
   * \param[in] val_- the .
   * \param[in] val_ivar        - eqn. index to the .
   */
  inline void SetLookupScalar(unsigned long  iPoint,
                            su2double      val_lookup_scalar,
                            unsigned short val_ivar) override {
    lookup_scalar(iPoint,val_ivar) = val_lookup_scalar;
  }

  /*!
   * \brief Get the value of the transported scalar source term
   * \param[in] val_ivar - eqn. index to the transported scalar source term
   * \return Value of the progress variable source term
   */
  inline su2double GetScalarSources(unsigned long  iPoint,
                                 unsigned short val_ivar) override {
    return source_scalar(iPoint,val_ivar);
  }

   /*!
   * \brief Get the value of the looked up scalar field
   * \param[in] val_ivar - eqn. index to the looked up scalar field
   * \return Value of the looked up scalar field
   */
  inline su2double GetScalarLookups(unsigned long  iPoint,
                                 unsigned short val_ivar) override {
    return lookup_scalar(iPoint,val_ivar);
  }
  

  /*!
   * \brief Get the value of the transported scalars source term
   * \return Pointer to the transported scalars source term
   */
  inline su2double *GetScalarSources(unsigned long iPoint) override {
    return source_scalar[iPoint];
  }

  /*!
   * \brief Get the value of the looked up table based on the transported scalar
   * \return Pointer to the transported scalars source term
   */
  inline su2double *GetScalarLookups(unsigned long iPoint) override {
    return lookup_scalar[iPoint];
  }
  
  
};
