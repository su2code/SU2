/*!
 * \file CFEAMeshNumerics.hpp
 * \brief Declaration and inlines of the class to compute
 *        the stiffness matrix of a linear, pseudo-elastic mesh problem.
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

#include "../numerics_structure.hpp"

class CFEAMeshElasticity : public CFEALinearElasticity {

  bool element_based;
  bool stiffness_set;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAMeshElasticity(unsigned short val_nDim, unsigned short val_nVar, unsigned long val_nElem, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEAMeshElasticity(void);

  inline void SetElement_Properties(CElement *element, CConfig *config){
    if(element_based){E = E_i[element->Get_iProp()];  Compute_Lame_Parameters();}
  }

  /*!
   * \brief Set the element-based local properties in mesh problems
   * \param[in] element_container - Element structure for the particular element integrated.
   */
  inline void SetMeshElasticProperties(unsigned long iElem, su2double val_E){
    if (element_based){ E_i[iElem]  = val_E;}
  }

};
