/*!
 * \file CFluidModel.cpp
 * \brief Source of the fluid model base class containing thermo-physical subroutines.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna, T. Economon
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

#include "../../include/fluid/CFluidModel.hpp"

#include <utility>

#include "../../include/fluid/CConstantConductivity.hpp"
#include "../../include/fluid/CConstantConductivityRANS.hpp"
#include "../../include/fluid/CConstantPrandtl.hpp"
#include "../../include/fluid/CConstantPrandtlRANS.hpp"
#include "../../include/fluid/CConstantViscosity.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CPolynomialConductivity.hpp"
#include "../../include/fluid/CPolynomialConductivityRANS.hpp"
#include "../../include/fluid/CPolynomialViscosity.hpp"
#include "../../include/fluid/CSutherland.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/fluid/CConstantDiffusivity.hpp"

unique_ptr<CViscosityModel> CFluidModel::MakeLaminarViscosityModel(const CConfig* config, unsigned short iSpecies) {
  switch (config->GetKind_ViscosityModel()) {
    case VISCOSITYMODEL::CONSTANT:
      return unique_ptr<CConstantViscosity>(new CConstantViscosity(config->GetMu_ConstantND(iSpecies)));
    case VISCOSITYMODEL::SUTHERLAND:
      return unique_ptr<CSutherland>(new CSutherland(config->GetMu_RefND(iSpecies),
                                                     config->GetMu_Temperature_RefND(iSpecies),
                                                     config->GetMu_SND(iSpecies)));
    case VISCOSITYMODEL::POLYNOMIAL:
      return unique_ptr<CPolynomialViscosity<N_POLY_COEFFS>>(
          new CPolynomialViscosity<N_POLY_COEFFS>(config->GetMu_PolyCoeffND()));
    case VISCOSITYMODEL::FLAMELET:
      /*--- Viscosity is obtained from the LUT ---*/
      break;
    default:
      SU2_MPI::Error("Viscosity model not available.", CURRENT_FUNCTION);
      return nullptr;
  }
}

void CFluidModel::SetLaminarViscosityModel(const CConfig* config) {
  LaminarViscosity = MakeLaminarViscosityModel(config, 0);
}

unique_ptr<CConductivityModel> CFluidModel::MakeThermalConductivityModel(const CConfig* config,
                                                                         unsigned short iSpecies) {
  switch (config->GetKind_ConductivityModel()) {
    case CONDUCTIVITYMODEL::CONSTANT:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CConstantConductivityRANS>(
            new CConstantConductivityRANS(config->GetThermal_Conductivity_ConstantND(iSpecies),
                                          config->GetPrandtl_Turb(iSpecies)));
      } else {
        return unique_ptr<CConstantConductivity>(
            new CConstantConductivity(config->GetThermal_Conductivity_ConstantND(iSpecies)));
      }
      break;
    case CONDUCTIVITYMODEL::CONSTANT_PRANDTL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CConstantPrandtlRANS>(
            new CConstantPrandtlRANS(config->GetPrandtl_Lam(iSpecies), config->GetPrandtl_Turb(iSpecies)));
      } else {
        return unique_ptr<CConstantPrandtl>(new CConstantPrandtl(config->GetPrandtl_Lam(iSpecies)));
      }
      break;
    case CONDUCTIVITYMODEL::POLYNOMIAL:
      if (config->GetKind_ConductivityModel_Turb() == CONDUCTIVITYMODEL_TURB::CONSTANT_PRANDTL) {
        return unique_ptr<CPolynomialConductivityRANS<N_POLY_COEFFS>>(
            new CPolynomialConductivityRANS<N_POLY_COEFFS>(config->GetKt_PolyCoeffND(), config->GetPrandtl_Turb()));
      } else {
        return unique_ptr<CPolynomialConductivity<N_POLY_COEFFS>>(
            new CPolynomialConductivity<N_POLY_COEFFS>(config->GetKt_PolyCoeffND()));
      }
      break;
    case CONDUCTIVITYMODEL::FLAMELET:
      /*--- conductivity is obtained from the LUT ---*/
      break;
    default:
      SU2_MPI::Error("Conductivity model not available.", CURRENT_FUNCTION);
      return nullptr;
      break;
  }
}

void CFluidModel::SetThermalConductivityModel(const CConfig* config) {
  ThermalConductivity = MakeThermalConductivityModel(config, 0);
}


void CFluidModel::SetMassDiffusivityModel (const CConfig* config) {
  switch (config->GetKind_DiffusivityModel()) {
    // nijso TODO 
    case DIFFUSIVITYMODEL::CONSTANT_DIFFUSIVITY:
      MassDiffusivity = unique_ptr<CConstantDiffusivity>(new CConstantDiffusivity(config->GetDiffusivity_ConstantND()));
      break;
    //case DIFFUSIVITYMODEL::CONSTANT_SCHMIDT:
    //  if ((config->GetKind_Solver() == RANS) || (config->GetKind_Solver() == DISC_ADJ_RANS)) {
    //    MassDiffusivity = unique_ptr<CConstantSchmidtRANS>(new CConstantSchmidtRANS(config->GetSchmidt_Lam(),config->GetSchmidt_Turb()));
    //  } else {
    //    MassDiffusivity = unique_ptr<CConstantSchmidt>(new CConstantSchmidt(config->GetSchmidt_Lam()));
    //  }
    //  break;
    case DIFFUSIVITYMODEL::FLAMELET:
      /* do nothing. Diffusivity is obtained from the table and set in setTDState_T */
      break;
    default:
      SU2_MPI::Error("Diffusivity model not available.", CURRENT_FUNCTION);
      break;
  }
}
