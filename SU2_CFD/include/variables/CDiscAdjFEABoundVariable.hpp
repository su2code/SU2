/*!
 * \file CDiscAdjFEABoundVariable.hpp
 * \brief Main class for defining the variables of the adjoint FEA solver at the boundary.
 * \author Ruben Sanchez
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

#include "CDiscAdjFEAVariable.hpp"

/*!
 * \class CDiscAdjFEABoundVariable
 * \brief Main class for defining the variables on the FEA boundaries for adjoint applications.
 * \ingroup Discrete_Adjoint
 * \author R. Sanchez.
 * \version 6.2.0 "Falcon"
 */
class CDiscAdjFEABoundVariable : public CDiscAdjFEAVariable {
protected:

  su2double *FlowTraction_Sens;         /*!< \brief Adjoint of the flow tractions. */
  su2double *SourceTerm_DispAdjoint;    /*!< \brief Source term applied into the displacement adjoint
                                                    coming from external solvers. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CDiscAdjFEABoundVariable(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEABoundVariable(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEABoundVariable(void);

  /*!
   * \brief Set the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetFlowTractionSensitivity(unsigned short iDim, su2double val) { FlowTraction_Sens[iDim] = val; }

  /*!
   * \brief Get the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetFlowTractionSensitivity(unsigned short iDim) { return FlowTraction_Sens[iDim]; }

  /*!
   * \brief Set the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \param[in] val - value of the source term
   */
  inline void SetSourceTerm_DispAdjoint(unsigned short iDim, su2double val){
    SourceTerm_DispAdjoint[iDim] = val;
  }

  /*!
   * \brief Get the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \return value of the source term
   */
  inline su2double GetSourceTerm_DispAdjoint(unsigned short iDim){
    return SourceTerm_DispAdjoint[iDim];
  }

  /*!
   * \brief Get whether this node is on the boundary
   */
  inline bool Get_isVertex(void) { return true; }

};
