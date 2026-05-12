/*!
 * \file CPBIncNSVariable.hpp
 * \brief Class for defining the variables of the pressure based incompressible
          Navier-Stokes solver.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "CPBIncEulerVariable.hpp"

/*!
 * \class CPBIncNSVariable
 * \brief Class for defining the variables of the incompressible Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CPBIncNSVariable final : public CPBIncEulerVariable {
private:
  const bool Energy;          /*!< \brief Flag for Energy equation in incompressible flows. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CPBIncNSVariable(su2double pressure, const su2double *velocity, su2double temperature,
                 unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);


};
