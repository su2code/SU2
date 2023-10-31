/*!
 * \file CCoolPropConductivity.hpp
 * \brief Defines laminar thermal conductivity model from CoolProp.
 * \author P.YAn, G. Gori, A. Guardone
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

#include "CConductivityModel.hpp"
#include "CCoolProp.hpp"

#ifdef USE_COOLPROP
#include "AbstractState.h"
#include "CoolProp.h"
#endif
/*!
 * \class CCoolPropConductivity
 * \brief Defines conductivity model from CoolProp.
 * \author P.Yan
 */
class CCoolPropConductivity final : public CConductivityModel {
 private:
#ifdef USE_COOLPROP
  std::unique_ptr<CoolProp::AbstractState> fluid_entity; /*!< \brief fluid entity */
#endif
  /*!< \brief list of fluids whose conductivity model is available in CoolProp */
  const vector<string> fluidNameList = {
      "Air",      "Ammonia",  "Argon",      "Benzene",   "CarbonDioxide", "Cyclopentane", "EthylBenzene",
      "Ethane",   "Ethanol",  "HeavyWater", "Helium",    "Hydrogen",      "IsoButane",    "Isopentane",
      "Methane",  "Methanol", "Nitrogen",   "Oxygen",    "Propylene",     "ParaHydrogen", "R11",
      "R116",     "R12",      "R123",       "R1234yf",   "R1234ze(E)",    "R124",         "R125",
      "R13",      "R134a",    "R14",        "R141b",     "R142b",         "R143a",        "R152A",
      "R218",     "R22",      "R227EA",     "R23",       "R236EA",        "R236FA",       "R245fa",
      "R32",      "R404A",    "R407C",      "R410A",     "R507A",         "RC318",        "SulfurHexafluoride",
      "Toluene",  "Water",    "m-Xylene",   "n-Butane",  "n-Decane",      "n-Dodecane",   "n-Heptane",
      "n-Hexane", "n-Nonane", "n-Octane",   "n-Pentane", "n-Propane",     "o-Xylene",     "p-Xylene"};

 public:
  /*!
   * \brief Constructor of the class.
   */
  CCoolPropConductivity(const string& fluidname) {
#ifdef USE_COOLPROP
    fluid_entity = std::unique_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS", fluidname));
    if (std::find(fluidNameList.begin(), fluidNameList.end(), fluidname) == fluidNameList.end()) {
      SU2_MPI::Error("Conductivity model not available for this fluid in CoolProp library.", CURRENT_FUNCTION);
    }
#else
    SU2_MPI::Error(
        "SU2 was not compiled with CoolProp (-Denable-coolprop=true). Note that CoolProp cannot be used with "
        "directdiff or autodiff",
        CURRENT_FUNCTION);
#endif
  }

  /*!
   * \brief Set thermal conductivity.
   */
  void SetConductivity(su2double t, su2double rho, su2double, su2double, su2double, su2double, su2double) override {
#ifdef USE_COOLPROP
    fluid_entity->update(CoolProp::DmassT_INPUTS, rho, t);
    kt_ = fluid_entity->conductivity();
#endif
  }
};
