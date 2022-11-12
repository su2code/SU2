/*!
* \file CCoolPropViscosity.hpp
* \brief Defines CoolPropviscosity model.
* \author P.Yan, G. Gori, A. Guardone,
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

#include "CViscosityModel.hpp"
#include "CCoolProp.hpp"
#include "CoolProp.h"
#include "AbstractState.h"
/*!
* \class CCoolPropViscosity
* \brief Defines CoolProp viscosity model.
* \author P.Yan
*/
class CCoolPropViscosity final : public CViscosityModel {
 private:
  std::unique_ptr<CoolProp::AbstractState> fluid_entity;   /*!< \brief fluid entity */
public:
 /*!
  * \brief Constructor of the class.
  */
 CCoolPropViscosity (string fluidname) {
   fluid_entity = std::unique_ptr<CoolProp::AbstractState>(CoolProp::AbstractState::factory("HEOS",fluidname));
 }


 /*!
  * \brief Set Viscosity.
  */
 void SetViscosity(su2double t, su2double rho) override {
   fluid_entity->update(CoolProp::DmassT_INPUTS, rho, t);
   if(fluid_entity->name()=="Air" || fluid_entity->name() == "Ammonia" || fluid_entity->name() == "Argon" || fluid_entity->name() == "Benzene" || fluid_entity->name() == "CarbonDioxide"
       || fluid_entity->name() == "CycloHexane" || fluid_entity->name() == "Cyclopentane" || fluid_entity->name() == "DimethylEther" || fluid_entity->name() == "Ethane"
       || fluid_entity->name() == "Ethanol" || fluid_entity->name() == "HeavyWater" || fluid_entity->name() == "Helium" || fluid_entity->name() == "Hydrogen"
       || fluid_entity->name() == "HydrogenSulfide" || fluid_entity->name() == "IsoButane" || fluid_entity->name() == "Isopentane" || fluid_entity->name() == "Methane"
       || fluid_entity->name() == "Methanol" || fluid_entity->name() == "Nitrogen" || fluid_entity->name() == "Oxygen" || fluid_entity->name() == "Propylene"
       || fluid_entity->name() == "R11" || fluid_entity->name() == "R116" || fluid_entity->name() == "R12" || fluid_entity->name() == "R123"
       || fluid_entity->name() == "R1233zd(E)" || fluid_entity->name() == "R1234yf" || fluid_entity->name() == "R1234ze(E)" || fluid_entity->name() == "R124"
       || fluid_entity->name() == "R125" || fluid_entity->name() == "R13" || fluid_entity->name() == "R134a" || fluid_entity->name() == "R14"
       || fluid_entity->name() == "R141b" || fluid_entity->name() == "R142b" || fluid_entity->name() == "R143a" || fluid_entity->name() == "R152A"
       || fluid_entity->name() == "R218" || fluid_entity->name() == "R22" || fluid_entity->name() == "R227EA" || fluid_entity->name() == "R23"
       || fluid_entity->name() == "R236EA" || fluid_entity->name() == "R236FA" || fluid_entity->name() == "R245fa" || fluid_entity->name() == "R32"
       || fluid_entity->name() == "R404A" || fluid_entity->name() == "R407C" || fluid_entity->name() == "R410A" || fluid_entity->name() == "R507A"
       || fluid_entity->name() == "RC318" || fluid_entity->name() == "SulfurHexafluoride" || fluid_entity->name() == "Toluene" || fluid_entity->name() == "Water"
       || fluid_entity->name() == "m-Xylene" || fluid_entity->name() == "n-Butane" || fluid_entity->name() == "n-Decane" || fluid_entity->name() == "n-Dodecane"
       || fluid_entity->name() == "n-Heptane" || fluid_entity->name() == "n-Hexane" || fluid_entity->name() == "n-Nonane" || fluid_entity->name() == "n-Octane"
       || fluid_entity->name() == "n-Pentane" || fluid_entity->name() == "n-Propane" || fluid_entity->name() == "o-Xylene" || fluid_entity->name() == "p-Xylene") {
     mu_ = fluid_entity->viscosity();
   }
   else
   {
     SU2_MPI::Error("Viscosity model not available for this fluid in CoolProp library.", CURRENT_FUNCTION);
   }
 }
};
