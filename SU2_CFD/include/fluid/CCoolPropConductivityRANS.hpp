/*!
* \file CCoolPropConductivityRANS.hpp
* \brief Defines conductivity model for RANS problems from CoolProp.
* \author P.YAn
* \version 7.4.0 "Blackbird"
*
* SU2 Project Website: https://su2code.github.io
*
* The SU2 Project is maintained by the SU2 Foundation
* (http://su2foundation.org)
*
* Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#include "CConductivityModel.hpp"
#include "CCoolProp.hpp"

#ifdef USE_COOLPROP
#include "CoolProp.h"
#include "AbstractState.h"

/*!
* \class CCoolPropConductivityRANS
* \brief Defines laminar thermal conductivity from CoolProp along
*        with a turbulent Prandtl number for including effects of
*        turbulent heat transfer. Returns the effective conductivity.
* \author P.Yan
*/
class CCoolPropConductivityRANS final : public CConductivityModel {
public:
 /*!
  * \brief Constructor of the class.
  */
 CCoolPropConductivityRANS(su2double kt_lam_const, su2double pr_turb)
     : kt_lam_const_(kt_lam_const), pr_turb_(pr_turb) {}

 /*!
  * \brief Set thermal conductivity.
  */
 void SetConductivity(su2double, su2double, su2double, su2double mu_turb, su2double cp,
                      su2double, su2double) override {
   kt_ = kt_lam_const_ + cp * mu_turb / pr_turb_;
 }

private:
 const su2double kt_lam_const_{0.0}; /*!< \brief Constant laminar conductivity. */
 const su2double pr_turb_{0.0};      /*!< \brief Turbulent Prandtl number. */
};
#endif