/*!
 * \file CFEABoundVariable.hpp
 * \brief Class for defining the variables on the FEA boundaries for FSI applications.
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

#include "CFEAVariable.hpp"

/*!
 * \class CFEABoundVariable
 * \brief Class for defining the variables on the FEA boundaries for FSI applications.
 * \ingroup Structural Finite Element Analysis Variables
 * \author R. Sanchez.
 * \version 6.2.0 "Falcon"
 */
class CFEABoundVariable : public CFEAVariable {
protected:

  su2double *FlowTraction;        /*!< \brief Traction from the fluid field. */
  su2double *FlowTraction_n;      /*!< \brief Traction from the fluid field at time n. */

  su2double *Residual_Ext_Surf;   /*!< \brief Term of the residual due to external forces */
  su2double *Residual_Ext_Surf_n; /*!< \brief Term of the residual due to external forces at time n */

public:

  /*!
   * \brief Constructor of the class.
   */
  CFEABoundVariable(void);

  /*!
   * \overload
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEABoundVariable(su2double *val_fea, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEABoundVariable(void);

  /*!
   * \brief Add surface load to the residual term
   */
  inline void Add_SurfaceLoad_Res(su2double *val_surfForce) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_Ext_Surf[iVar] += val_surfForce[iVar];
  }

  /*!
   * \brief Set surface load of the residual term (for dampers - deletes all the other loads)
   */
  inline void Set_SurfaceLoad_Res(unsigned short iVar, su2double val_surfForce) {Residual_Ext_Surf[iVar] = val_surfForce;}

  /*!
   * \brief Get the residual term due to surface load
   */
  inline su2double Get_SurfaceLoad_Res(unsigned short iVar) {return Residual_Ext_Surf[iVar];}

  /*!
   * \brief Clear the surface load residual
   */
  inline void Clear_SurfaceLoad_Res(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) Residual_Ext_Surf[iVar] = 0.0;
  }

  /*!
   * \brief Store the surface load as the load for the previous time step.
   */
  inline void Set_SurfaceLoad_Res_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)  Residual_Ext_Surf_n[iVar] = Residual_Ext_Surf[iVar];
  }

  /*!
   * \brief Get the surface load from the previous time step.
   */
  inline su2double Get_SurfaceLoad_Res_n(unsigned short iVar) {return Residual_Ext_Surf_n[iVar]; }

  /*!
   * \brief Set the flow traction at a node on the structural side
   */
  inline void Set_FlowTraction(su2double *val_flowTraction) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) FlowTraction[iVar] = val_flowTraction[iVar];
  }

  /*!
   * \brief Add a value to the flow traction at a node on the structural side
   */
  inline void Add_FlowTraction(su2double *val_flowTraction) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) FlowTraction[iVar] += val_flowTraction[iVar];
  }

  /*!
   * \brief Get the residual term due to the flow traction
   */
  inline su2double Get_FlowTraction(unsigned short iVar) {return FlowTraction[iVar]; }

  /*!
   * \brief Set the value of the flow traction at the previous time step.
   */
  void Set_FlowTraction_n(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) FlowTraction_n[iVar] = FlowTraction[iVar];
  }

  /*!
   * \brief Retrieve the value of the flow traction from the previous time step.
   */
  inline su2double Get_FlowTraction_n(unsigned short iVar) {return FlowTraction_n[iVar]; }

  /*!
   * \brief Clear the flow traction residual
   */
  inline void Clear_FlowTraction(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) FlowTraction[iVar] = 0.0;
  }

  /*!
   * \brief Get whether this node is on the boundary
   */
  inline bool Get_isVertex(void) {return true; }

};
