/*!
 * \file CUpwAUSMPLUSUP2_Flow.hpp
 * \brief Delaration of numerics class CUpwAUSMPLUSUP2_Flow, the
 *        implementation is in the CUpwAUSMPLUSUP2_Flow.cpp file.
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
 * \class CUpwAUSMPLUSUP2_Flow
 * \brief Class for solving an approximate Riemann AUSM+ -up.
 * \ingroup ConvDiscr
 * \author Amit Sachdeva
 */
class CUpwAUSMPLUSUP2_Flow : public CUpwAUSMPLUS_SLAU_Base_Flow {
private:
  su2double Kp, sigma;

  /*!
   * \brief Mass flux and pressure for the AUSM+up2 scheme.
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
  CUpwAUSMPLUSUP2_Flow(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwAUSMPLUSUP2_Flow(void);

};
