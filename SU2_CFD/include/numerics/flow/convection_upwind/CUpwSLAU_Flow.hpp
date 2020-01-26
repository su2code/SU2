/*!
 * \file CUpwSLAU_Flow.hpp
 * \brief Delaration of numerics class CUpwSLAU_Flow, the
 *        implementation is in the CUpwSLAU_Flow.cpp file.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
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

#include "CUpwAUSMPLUS_SLAU_Base_Flow.hpp"

/*!
 * \class CUpwSLAU_Flow
 * \brief Class for solving the Low-Dissipation AUSM.
 * \ingroup ConvDiscr
 * \author E. Molina
 */
class CUpwSLAU_Flow : public CUpwAUSMPLUS_SLAU_Base_Flow {
protected:
  bool slau_low_diss;
  bool slau2;
  
  /*!
   * \brief Mass flux and pressure for the SLAU and SLAU2 schemes.
   * \param[in] config - Definition of the particular problem.
   * \param[out] mdot - The mass flux.
   * \param[out] pressure - The pressure at the control volume face.
   */
  void ComputeMassAndPressureFluxes(CConfig *config, su2double &mdot, su2double &pressure);
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSLAU_Flow(void);

};

/*!
 * \class CUpwSLAU2_Flow
 * \brief Class for solving the Simple Low-Dissipation AUSM 2.
 * \ingroup ConvDiscr
 * \author E. Molina
 */
class CUpwSLAU2_Flow : public CUpwSLAU_Flow {
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSLAU2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config, bool val_low_dissipation);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSLAU2_Flow(void);

};
