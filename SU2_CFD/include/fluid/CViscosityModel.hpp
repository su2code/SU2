/*!
 * \file CViscosityModel.hpp
 * \brief Interface class for defining laminar viscosity models.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
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

#include "../../../Common/include/basic_types/datatype_structure.hpp"

using namespace std;

/*!
 * \class CViscosityModel
 * \brief Interface class for defining the laminar viscosity model.
 * \author S.Vitale, M.Pini
 */
class CViscosityModel {
 public:
  CViscosityModel() = default;
  CViscosityModel(const CViscosityModel&) = delete;
  void operator=(const CViscosityModel&) = delete;
  virtual ~CViscosityModel() {}

  /*!
   * \brief return viscosity value.
   */
  inline su2double GetViscosity() const { return mu_; }

  /*!
   * \brief return viscosity partial derivative value.
   */
  inline su2double Getdmudrho_T() const { return dmudrho_t_; }

  /*!
   * \brief return viscosity partial derivative value.
   */
  inline su2double GetdmudT_rho() const { return dmudt_rho_; }

  /*!
   * \brief Set Viscosity.
   */
  virtual void SetViscosity(su2double t, su2double rho) = 0;

 protected:
  su2double mu_{0.0};        /*!< \brief Dynamic viscosity. */
  su2double dmudrho_t_{0.0}; /*!< \brief DmuDrho_T. */
  su2double dmudt_rho_{0.0}; /*!< \brief DmuDT_rho. */
};
