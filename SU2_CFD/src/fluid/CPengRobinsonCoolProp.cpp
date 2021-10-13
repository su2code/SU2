/*!
 * \file CPengRobinsonCoolProp.cpp
 * \brief Source of the Cubic Peng-Robinson Equation of State model based on CoolProp library:
 * \ http://coolprop.org/
 * \ https://github.com/CoolProp/CoolProp
 *
 * \author Dmitry Brezgin
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/fluid/CPengRobinsonCoolProp.hpp"
#include "../../include/fluid/CPengRobinsonCoolPropAncillaries.hpp" 



CPengRobinsonCoolProp::CPengRobinsonCoolProp(std::string fluidName) : CFluidModel() {

  /*transform to uppper case letters*/
  for (auto& c : fluidName) c = (char)toupper(c);

  std::string errstr, el_name;
  bool isFluidFound = false;

  /* First we validate the json string against the schema */
  cpjson::schema_validation_code val_code = cpjson::validate_schema(all_cubics_extended_JSON, all_cubics_extended_JSON, errstr);
  /* Then we check the validation code */
  if (val_code == cpjson::SCHEMA_VALIDATION_OK) {
    rapidjson::Document dd;

    dd.Parse<0>(all_cubics_extended_JSON.c_str());
    if (dd.HasParseError()) {
      SU2_MPI::Error("Cubics JSON is not valid JSON", CURRENT_FUNCTION);
    } else {            

      for (rapidjson::Value::ValueIterator itr = dd.Begin(); itr != dd.End(); ++itr) {

        el_name = cpjson::get_string(*itr, "name");
        if (el_name.compare(fluidName) != 0) continue;

        isFluidFound = true;

        acentric = cpjson::get_double(*itr, "acentric");
        Molar_mass = cpjson::get_double(*itr, "molemass");
        Tcrit = cpjson::get_double(*itr, "Tc");
        Pcrit = cpjson::get_double(*itr, "pc");
        Rhocrit_molar = cpjson::get_double(*itr, "rhomolarc");
        Rgas = Ru / Molar_mass;

        if (itr->HasMember("pS") && (*itr)["pS"].IsObject()) {
          rapidjson::Value& pS = (*itr)["pS"];
          std::vector<su2double> n = cpjson::get_double_array(pS, "n");
          std::vector<su2double> t = cpjson::get_double_array(pS, "t");
          su2double pS_reducing = cpjson::get_double(pS, "reducing_value");
          su2double Tr = cpjson::get_double(pS, "T_r");
          su2double Tmax = cpjson::get_double(pS, "Tmax");
          su2double Tmin = cpjson::get_double(pS, "Tmin");
          bool using_tau_r = cpjson::get_bool(pS, "using_tau_r");

          PsContainer = SaturationLinelPs(n, t, Tr, Tmax, Tmin, pS_reducing, using_tau_r);
        }

        if (itr->HasMember("triple_vapor") && (*itr)["triple_vapor"].IsObject()) {
          rapidjson::Value& triple_data = (*itr)["triple_vapor"];
          Ttriple = cpjson::get_double(triple_data, "T");
          Rhotriple = cpjson::get_double(triple_data, "rhomolar") * Molar_mass;
        }

        if (itr->HasMember("alpha0") && (*itr)["alpha0"].IsArray()) {
          for (rapidjson::Value::ConstValueIterator itrA = (*itr)["alpha0"].Begin(); itrA != (*itr)["alpha0"].End();
               ++itrA) {
            /*A reference for code cleanness*/
            const rapidjson::Value& contribution = *itrA;

            /* Get the type (required!)*/
            std::string type = contribution["type"].GetString();
            if (!type.compare("IdealGasHelmholtzLead")) {
              su2double a1 = cpjson::get_double(contribution, "a1");
              su2double a2 = cpjson::get_double(contribution, "a2");

              Lead = IdealHelmholtzLead(a1, a2);
            } else if (!type.compare("IdealGasHelmholtzLogTau")) {
              su2double a = cpjson::get_double(contribution, "a");

              LogTau = IdealHelmholtzLogTau(a);
            } else if (!type.compare("IdealGasHelmholtzPlanckEinsteinGeneralized")) {
              /*Retrieve the values*/
              std::vector<su2double> n = cpjson::get_long_double_array(contribution["n"]);
              std::vector<su2double> t = cpjson::get_long_double_array(contribution["t"]);

              std::vector<su2double> c = cpjson::get_long_double_array(contribution["c"]);
              std::vector<su2double> d = cpjson::get_long_double_array(contribution["d"]);

              if (PlanckEinstein.is_enabled() == true) {
                PlanckEinstein.extend(n, t, c, d);
              } else {
                PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
              }
            } else if (!type.compare("IdealGasHelmholtzPlanckEinstein")) {
              /*Retrieve the values*/
              std::vector<su2double> n = cpjson::get_long_double_array(contribution["n"]);
              std::vector<su2double> t = cpjson::get_long_double_array(contribution["t"]);
              /* Flip the sign of theta */
              for (std::size_t i = 0; i < t.size(); ++i) {
                t[i] *= -1;
              }
              std::vector<su2double> c(n.size(), 1);
              std::vector<su2double> d(c.size(), -1);
              if (PlanckEinstein.is_enabled() == true) {
                PlanckEinstein.extend(n, t, c, d);
              } else {
                PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
              }
            } else if (!type.compare("IdealGasHelmholtzEnthalpyEntropyOffset")) {
              su2double a1 = cpjson::get_double(contribution, "a1");
              su2double a2 = cpjson::get_double(contribution, "a2");
              std::string reference = cpjson::get_string(contribution, "reference");

              EnthEntrOffsetCore = IdealHelmholtzEnthalpyEntropyOffset(a1, a2, reference);
            } else if (!type.compare("IdealGasHelmholtzPower")) {
              std::vector<su2double> n = cpjson::get_long_double_array(contribution["n"]);
              std::vector<su2double> t = cpjson::get_long_double_array(contribution["t"]);

              Power = IdealHelmholtzPower(n, t);
            } else if (!type.compare("IdealGasHelmholtzCP0PolyT")) {
              std::vector<su2double> c = cpjson::get_long_double_array(contribution["c"]);
              std::vector<su2double> t = cpjson::get_long_double_array(contribution["t"]);
              su2double Tc = cpjson::get_double(contribution, "Tc");
              su2double T0 = cpjson::get_double(contribution, "T0");

              CP0PolyT = IdealHelmholtzCP0PolyT(c, t, Tc, T0);
            } else if (!type.compare("IdealGasHelmholtzCP0Constant")) {
              su2double cp_over_R = cpjson::get_double(contribution, "cp_over_R");
              su2double Tc = cpjson::get_double(contribution, "Tc");
              su2double T0 = cpjson::get_double(contribution, "T0");

              CP0Constant = IdealHelmholtzCP0Constant(cp_over_R, Tc, T0);
            } else if (!type.compare("IdealGasHelmholtzCP0AlyLee")) {
              std::vector<su2double> constants = cpjson::get_long_double_array(contribution["c"]);
              su2double Tc = cpjson::get_double(contribution, "Tc");
              su2double T0 = cpjson::get_double(contribution, "T0");

              /* Take the constant term if nonzero and set it as a polyT term */
              if (abs(constants[0]) > 1e-14) {
                std::vector<su2double> c(1, constants[0]), t(1, 0);
                if (CP0PolyT.is_enabled() == true) {
                  CP0PolyT.extend(c, t);
                } else {
                  CP0PolyT = IdealHelmholtzCP0PolyT(c, t, Tc, T0);
                }
              }
              std::vector<su2double> n, c, d, t;
              if (abs(constants[1]) > 1e-14) {
                n.push_back(constants[1]);
                t.push_back(-2 * constants[2] / Tc);
                c.push_back(1);
                d.push_back(-1);
              }
              if (abs(constants[3]) > 1e-14) {
                n.push_back(-constants[3]);
                t.push_back(-2 * constants[4] / Tc);
                c.push_back(1);
                d.push_back(1);
              }
              if (PlanckEinstein.is_enabled() == true) {
                PlanckEinstein.extend(n, t, c, d);
              } else {
                PlanckEinstein = IdealHelmholtzPlanckEinsteinGeneralized(n, t, c, d);
              }
            } else {
              SU2_MPI::Error(std::string("Unsupported ideal-gas Helmholtz type: %s\n") + type.c_str(), CURRENT_FUNCTION);
            }
          }
        }

        Zed = 1.0;

        /* !!!!!!  Peng-Robinson supplementary functions  */
        a = 0.45724 * Rgas * Rgas * Tcrit * Tcrit / Pcrit; /*it uses Specific gas constant*/
        b = 0.0778 * Rgas * Tcrit / Pcrit;
        k = 0.37464 + 1.54226 * acentric - 0.26992 * acentric * acentric;
        a0 = 0.45724 * Ru * Ru * Tcrit * Tcrit / Pcrit; /*it uses Universal gas constant*/
        /* !!!!!! */

        Tr_over_Tci = 1 / Tcrit;
        sqrt_Tr_Tci = sqrt(Tr_over_Tci);
      }
      if (!isFluidFound) SU2_MPI::Error("Couldn't find data for fluid: " + fluidName, CURRENT_FUNCTION);
    }
  } else {
    SU2_MPI::Error("Unable to validate cubics library against schema with error: ", CURRENT_FUNCTION);
  }
}

CPengRobinsonCoolProp::~CPengRobinsonCoolProp() {}

void CPengRobinsonCoolProp::SetTDState_rhoe(su2double rho, su2double e) {
  StaticEnergy = e;
  su2double temperature = GetTemperatureByRhoUmass(rho, StaticEnergy);
  SetInnerTDState_rhoT(rho, temperature, true, true, false);
}

su2double CPengRobinsonCoolProp::Estimate_Ps_SatLine(su2double T) {
  su2double Ps;
  Ps = PsContainer.Evaluate(T);
  return Ps;
}

su2double CPengRobinsonCoolProp::Get_Saturation_Temperature(su2double P) {
  su2double f0, f0Prime, f0left, f0right, temperature, deltax, Ts_est, theta;

  deltax = 1E-6;

  /*  initial guess of saturation temperature */
  theta = -log10(P / Pcrit) * (1.0 / 0.7 - 1.0) / (acentric + 1.0);
  Ts_est = Tcrit / (theta + 1.0);

  su2double z0 = Ts_est;
  su2double relError = 1.0;
  su2double z = z0;
  /*  iteration counter */
  unsigned short iter = 0;

  while (abs(relError) > 1E-3) {
    /* function */
    f0 = P - Estimate_Ps_SatLine(z);
    /* function to the left */
    f0left = P - Estimate_Ps_SatLine(z - deltax * 0.5);
    /*  function to the right */
    f0right = P - Estimate_Ps_SatLine(z + deltax * 0.5);
    /* derivative at the point */
    f0Prime = (f0right - f0left) / deltax;

    z = z0 - f0 / f0Prime;
    relError = (z - z0) / z;
    z0 = z;

    iter++;
  }
  temperature = z;
  return temperature;
}

su2double CPengRobinsonCoolProp::GetTemperatureByPRho(su2double P, su2double rho) {
  su2double A, B, C, T, vb1, vb2;

  AD::StartPreacc();
  AD::SetPreaccIn(P);
  AD::SetPreaccIn(rho);

  vb1 = (1 / rho - b);
  vb2 = (1 / rho / rho + 2 * b / rho - b * b);

  A = Rgas / vb1 - a * k * k / Tcrit / vb2;
  B = 2 * a * k * (k + 1) / sqrt(Tcrit) / vb2;
  C = -P - a * (1 + k) * (1 + k) / vb2;

  T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
  T *= T;
  return T;
}

su2double CPengRobinsonCoolProp::GetTemperatureByRhoUmass(su2double rho, su2double umass) {
  su2double f0, f0Prime, f0left, f0right, temperature, deltax;
  
  deltax = 0.000001;

  /* initial guess of temperature */
  su2double z0 = Ttriple;
  su2double relError = 1.0;
  su2double z = z0;
  unsigned short iter = 0;

  while (abs(relError) > 1E-3) {
    /* function*/
    f0 = umass - GetUmassByRhoT(rho, z);
    /* function to the left */
    f0left = umass - GetUmassByRhoT(rho, z - deltax * 0.5);
    /* function to the right */
    f0right = umass - GetUmassByRhoT(rho, z + deltax * 0.5);
    /* derivative at the point */
    f0Prime = (f0right - f0left) / deltax;

    z = z0 - f0 / f0Prime;
    relError = (z - z0) / z;
    z0 = z;

    iter++;
  }
  temperature = z;
  return temperature;
}

void CPengRobinsonCoolProp::SetTDState_PT(su2double P, su2double T) {
  Pressure = P;
 
  /*su2double rhoV0(0.0), rhoV1(0.0), rhoV2(0.0);
  su2double theta(0.0), Ts_est(0.0), RhoV_est(0.0);
  int NsolnV = 0;*/
  
  su2double rho0(0.0), rho1(0.0), rho2(0.0), rho(-1.0);
  int Nsoln = 0;

  Rho_Tp_cubic(T, P, Nsoln, rho0, rho1, rho2); /*  Densities are sorted in increasing order */

  std::string msg;
  std::stringstream Ts, Ps;
  
  if (Nsoln == 1) {
    rho = rho0;
  } else if (Nsoln == 3) {
    /* Use imposed phase to select root */
    if (rho0 > 0) {
      rho = rho0;
    } else if (rho1 > 0) {
      rho = rho1;
    } else if (rho2 > 0) {
      rho = rho2;
    } else {
      Ts << T;
      Ps << P;
      SU2_MPI::Error(std::string("Unable to find gaseous density for T: ") + Ts.str() + " K, P: " + Ps.str() + " Pa",  CURRENT_FUNCTION);
    }

    /* If phase is not imposed */
    /* if (P < Pcrit) {

        Ts_est = Get_Saturation_Temperature(P);
        if (T > Ts_est)
        {
            // There are three density solutions, we know the highest is the liquid, the lowest is the vapor (rhoV0)
            Rho_Tp_cubic(Ts_est, P, NsolnV, rhoV0, rhoV1, rhoV2);
            // Gas
            if (rho0 > 0 && rho0 < rhoV0) {
                rho = rho0;
            }
            else if (rho1 > 0 && rho1 < rhoV0) {
                rho = rho1;
            }
            else {
                Ts << T; Ps << P;
                msg = "Unable to find gaseous density for T: '" + Ts.str() + "', P: " + Ps.str() + "'";
                SU2_MPI::Error(msg, CURRENT_FUNCTION);
            }
        }
        else {
            // Liquid
            //rho = rho2;
            SU2_MPI::Error(string("Liquid phase is detected (fluid Temperature is below the saturated one)\n") +
            string("Check the boundary or initial conditions!"), CURRENT_FUNCTION);
        }
    }
     else
        SU2_MPI::Error("Cubic has three roots, but phase not imposed and guess density not provided", CURRENT_FUNCTION); */
  }

  SetInnerTDState_rhoT(rho, T, false, true, true);
}

void CPengRobinsonCoolProp::SetInnerTDState_rhoT(su2double rho, su2double T, bool isSetPressure, bool isSetEntropy,
                                               bool isSetEnergy) {
  su2double rho_molar, taustar, deltastar, tau;
  su2double psi_plus0, psi_plus1, tau_times_a0, tau_times_a1;

  su2double dOf_dT, dOf_drho, dPdRho_Tconst, dPdT_Rhoconst, dUdRho_Tconst, DpDd_T, DpDT_d, DeDd_T;
  su2double umolar(0.0), umass(0.0), smolar(0.0), smass(0.0), cpmolar(0.0), cpmass(0.0), cvmolar(0.0), cvmass(0.0);
  su2double s(0.0), da0_dTau(0.0), dar_dTau(0.0), dar_dDelta(0.0), d2a0_dTau2(0.0), d2ar_dTau2(0.0),
      d2ar_dDelta_dTau(0.0), d2ar_dDelta2(0.0), alphar(0.0), alpha0(0.0);

  Density = rho;
  Temperature = T;

  AD::StartPreacc();
  AD::SetPreaccIn(Density);
  AD::SetPreaccIn(Temperature);

  tau = 1.0 / Temperature;
  rho_molar = 1.0 / Molar_mass * Density;
  taustar = Tcrit / Temperature;
  deltastar = 1.0 / Rhocrit_molar * rho_molar;

  psi_plus0 = psi_plus(rho_molar, 0);
  psi_plus1 = psi_plus(rho_molar, 1);
  tau_times_a0 = tau_times_a(tau, 0);
  tau_times_a1 = tau_times_a(tau, 1);

  HelmholtzDerivatives ders = CalcAllHelmholtzIdeal(taustar, deltastar);
  da0_dTau = ders.dalpha0_dtau;
  da0_dTau /= 1.0 / Tcrit;

  dar_dDelta = psi_minus(rho_molar, 0, 1) - tau_times_a0 / Ru * psi_plus1;
  if (isSetPressure) {
    Pressure = rho_molar * Ru * Temperature * (1 + rho_molar * dar_dDelta);
    AD::SetPreaccOut(Pressure);
  }

  dar_dTau = psi_minus(rho_molar, 1, 0) - tau_times_a1 / Ru * psi_plus0;

  if (isSetEnergy) {
    umolar = Ru * T * tau * (da0_dTau + dar_dTau);
    /*specific internal energy*/
    umass = umolar / Molar_mass;
    StaticEnergy = umass;
    AD::SetPreaccOut(StaticEnergy);
  }

  if (isSetEntropy) {
    alpha0 = ders.alpha0;
    alphar = psi_minus(rho_molar, 0, 0) - tau_times_a0 / Ru * psi_plus0;

    /* molar entropy */
    smolar = Ru * (tau * (da0_dTau + dar_dTau) - alpha0 - alphar);
    smass = smolar / Molar_mass;
    Entropy = smass;
    AD::SetPreaccOut(Entropy);
  }

  /* Cp*/
  d2a0_dTau2 = ders.d2alpha0_dtau2;
  d2a0_dTau2 /= pow(1.0 / Tcrit, 2.0);
  d2ar_dTau2 = psi_minus(rho_molar, 2, 0) - tau_times_a(tau, 2) / Ru * psi_plus0;
  d2ar_dDelta_dTau = psi_minus(rho_molar, 1, 1) - tau_times_a1 / Ru * psi_plus1;
  d2ar_dDelta2 = psi_minus(rho_molar, 0, 2) - tau_times_a0 / Ru * psi_plus(rho_molar, 2);
  cpmolar = Ru * (-pow(tau, 2.0) * (d2ar_dTau2 + d2a0_dTau2) +
                   pow(1.0 + rho_molar * dar_dDelta - rho_molar * tau * d2ar_dDelta_dTau, 2.0) /
                       (1.0 + 2.0 * rho_molar * dar_dDelta + pow(rho_molar, 2.0) * d2ar_dDelta2));
  cpmass = cpmolar / Molar_mass;

  /* Cv*/
  cvmolar = -Ru * pow(tau, 2.0) * (d2ar_dTau2 + d2a0_dTau2);
  cvmass = cvmolar / Molar_mass;

  /* DERIVATIVES*/
  dOf_drho = Ru * Temperature * (1.0 + 2.0 * rho_molar * dar_dDelta + pow(rho_molar, 2.0) * d2ar_dDelta2);
  dOf_dT = rho_molar * Ru * (1.0 + rho_molar * dar_dDelta - tau * rho_molar * d2ar_dDelta_dTau);
  dPdRho_Tconst = dOf_drho / Molar_mass;
  dPdT_Rhoconst = dOf_dT;
  dOf_drho = (Temperature * Ru / rho_molar * (tau * rho_molar * d2ar_dDelta_dTau)) / Molar_mass;
  dUdRho_Tconst = dOf_drho / Molar_mass;

  DpDd_T = dPdRho_Tconst;
  DpDT_d = dPdT_Rhoconst;
  DeDd_T = dUdRho_Tconst;

  Cp = cpmass;
  Cv = cvmass;

  dPde_rho = DpDT_d / Cv;
  dPdrho_e = DpDd_T - dPde_rho * DeDd_T;
  SoundSpeed2 = dPdrho_e + Pressure / (Density * Density) * dPde_rho;
  dTde_rho = 1 / Cv;
  Zed = Pressure / (Rgas * Temperature * Density);

  AD::SetPreaccOut(SoundSpeed2);
  AD::SetPreaccOut(dPde_rho);
  AD::SetPreaccOut(dPdrho_e);
  AD::SetPreaccOut(Zed);
  AD::SetPreaccOut(dTde_rho);

  AD::EndPreacc();
}

void CPengRobinsonCoolProp::SetTDState_Prho(su2double P, su2double rho) {
  Pressure = P;
  su2double T = GetTemperatureByPRho(P, rho);
  SetInnerTDState_rhoT(rho, T, false, true, true);
}

void CPengRobinsonCoolProp::SetEnergy_Prho(su2double P, su2double rho) {

  su2double T;
  su2double rho_molar, taustar, tau;
  su2double da0_dTau(0.0), dar_dTau(0.0), s(0.0), umolar(0.0), umass(0.0);

  AD::StartPreacc();
  AD::SetPreaccIn(P);
  AD::SetPreaccIn(rho);

  // vb1 = (1 / rho - b);
  // vb2 = (1 / rho / rho + 2 * b / rho - b * b);
  // A = Rgas / vb1 - a * k * k / Tcrit / vb2;
  // B = 2 * a * k * (k + 1) / sqrt(Tcrit) / vb2;
  // C = -P - a * (1 + k) * (1 + k) / vb2;
  // T = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
  // T *= T;
  // su2double T_matlab_symbolic_solve = (2 * Tcrit * a * k * vb1 * (k + 1) * (a * k * k * vb1 - sqrt(vb1 * vb2 * (Rgas
  // * Tcrit * a + Rgas * Tcrit * a * k * k
  //    - P * a * k * k * vb1 + P * Rgas * Tcrit * vb2 + 2 * Rgas * Tcrit * a * k)) + a * k * vb1)) / pow(a * k * k *
  //    vb1 - Rgas * Tcrit * vb2, 2)
  //    - (Tcrit * vb1 * (a * k * k + 2 * a * k + a + P * vb2)) / (a * vb1 * k * k - Rgas * Tcrit * vb2);

  T = GetTemperatureByPRho(P, rho);

  rho_molar = 1.0 / Molar_mass * rho;
  taustar = Tcrit / T;
  tau = 1.0 / T;

  da0_dTau = Get_da0_dTau(taustar);
  da0_dTau /= 1.0 / Tcrit;
  dar_dTau = -tau_times_a(tau, 1) / Ru * psi_plus(rho_molar, 1);

  umolar = Ru * T * tau * (da0_dTau + dar_dTau);
  /*specific internal energy*/
  umass = umolar / Molar_mass;

  StaticEnergy = umass;

  AD::SetPreaccOut(StaticEnergy);
  AD::EndPreacc();
}

void CPengRobinsonCoolProp::SetTDState_hs(su2double h, su2double s) {

  /* solves system of two nonlinear equations (density-temperature) by means of 2D Newton-Raphson method*/

  su2double cX(0.0), nX(0.0), cY(0.0), nY(0.0), f1(0.0), f2(0.0), jacA1(0.0), jacA2(0.0), jacB1(0.0), jacB2(0.0), detJac(0.0), detA1(0.0), detA2(0.0);
  su2double f1LeftX, f1RightX, f2LeftX, f2RightX, f1LeftY, f1RightY, f2LeftY, f2RightY;

  su2double tolerance = 1E-3;
  su2double delta = 1E-10;
  std::stringstream hs, ss;
  std::string msg;

  /* Initial guess of density.
  *  !!!
     In general, suggesting the vapor density at triple point is suitable.
     However, some fluids posses large triple point densities (e.g., carbon dioxide) and the possible workaround is to decrease the initial guess (e.g., 0.001)
  */
  cX = Rhotriple;
   /*initial guess of temperature (triple point)*/
  cY = Ttriple;
  su2double relError = 1.0;
   /*iteration counter*/
  unsigned short iter = 0;

  while (relError > tolerance) {
    if (cX < 0.0 || cY < 0.0) {
      hs << h;
      ss << s;
      SU2_MPI::Error("Check the boundary conditions!\n Enthalpy is '" + hs.str() + "', Entropy is '" + ss.str() + "'\n", CURRENT_FUNCTION);
    }

    f1 = h - GetHmassByRhoT(cX, cY);
    f2 = s - GetSmassByRhoT(cX, cY);

     /*function by tiny shift along X to the left*/
    f1LeftX = h - GetHmassByRhoT(cX - delta * 0.5, cY);
     /*function by tiny shift along X to the right*/
    f1RightX = h - GetHmassByRhoT(cX + delta * 0.5, cY);
    /*dh/dRho_T derivative at the point*/
    jacA1 = (f1RightX - f1LeftX) / delta;

     /*function by tiny shift along Y to the top*/
    f1LeftY = h - GetHmassByRhoT(cX, cY - delta * 0.5);
    /*function by tiny shift along Y to the bottom*/
    f1RightY = h - GetHmassByRhoT(cX, cY + delta * 0.5);
     /*dh/dT_Rho derivative at the point*/
    jacA2 = (f1RightY - f1LeftY) / delta;

    /* function by tiny shift along X to the left*/
    f2LeftX = s - GetSmassByRhoT(cX - delta * 0.5, cY);
    /* function by tiny shift along X to the right*/
    f2RightX = s - GetSmassByRhoT(cX + delta * 0.5, cY);
    /* ds/dRho_T derivative at the point*/
    jacB1 = (f2RightX - f2LeftX) / delta;

     /*function by tiny shift along Y to the top*/
    f2LeftY = s - GetSmassByRhoT(cX, cY - delta * 0.5);
     /*function by tiny shift along Y to the bottom*/
    f2RightY = s - GetSmassByRhoT(cX, cY + delta * 0.5);
     /*ds/dT_Rho derivative at the point*/
    jacB2 = (f2RightY - f2LeftY) / delta;

    detJac = jacA1 * jacB2 - jacB1 * jacA2;

    detA1 = f1 * jacB2 - jacA2 * f2;
    detA2 = jacA1 * f2 - jacB1 * f1;
    ;

    nX = cX - detA1 / detJac;
    nY = cY - detA2 / detJac;

    relError = max(abs((nX - cX) / cX), abs((nY - cY) / cY));
    if (isnan(relError)) {
      hs << h;
      ss << s;
      msg = "Check the initial guess of Density-Temperature pair: Enthalpy is '" + hs.str() + "', Entropy is '" +  ss.str() + "'";
      SU2_MPI::Error(msg, CURRENT_FUNCTION);
    }
    cX = nX;
    cY = nY;
    iter++;
    if (iter > 25) {
      hs << h;
      ss << s;
      msg = "Newton-Rapson solver is diverged (MAX number of iterations is exceeded)! Check the Enthalpy-Entropy pair: Enthalpy is '" + hs.str() + "', Entropy is '" + ss.str() + "'";
      SU2_MPI::Error(msg, CURRENT_FUNCTION);
    }
  }

  Entropy = s - f2;

  SetInnerTDState_rhoT(cX, cY, true, false, true);
}

su2double CPengRobinsonCoolProp::GetPressureByRhoT(su2double rho, su2double T) {
  su2double pressure, tau, rho_molar, taustar, dar_dDelta(0.0);

  tau = 1.0 / T;
  rho_molar = 1.0 / Molar_mass * rho;
  taustar = Tcrit / T;

  dar_dDelta = psi_minus(rho_molar, 0, 1) - tau_times_a(tau, 0) / Ru * psi_plus(rho_molar, 1);
  pressure = rho_molar * Ru * T * (1.0 + rho_molar * dar_dDelta);
  return pressure;
}

su2double CPengRobinsonCoolProp::GetUmassByRhoT(su2double rho, su2double T) {
  su2double tau(0.0), taustar(0.0), da0_dTau(0.0), dar_dTau(0.0), umolar(0.0), umass(0.0), rho_molar(0.0);

  rho_molar = 1.0 / Molar_mass * rho;
  taustar = Tcrit / T;
  tau = 1.0 / T;

  da0_dTau = Get_da0_dTau(taustar);
  da0_dTau /= 1.0 / Tcrit;
  dar_dTau = psi_minus(rho_molar, 1, 0) - tau_times_a(tau, 1) / Ru * psi_plus(rho_molar, 0);

  umolar = Ru * T * tau * (da0_dTau + dar_dTau);
  umass = umolar / Molar_mass;
  return umass;
}

su2double CPengRobinsonCoolProp::GetHmassByRhoT(su2double rho, su2double T) {
  su2double tau, rho_molar, taustar;
  su2double dar_dDelta(0.0), da0_dTau(0.0), dar_dTau(0.0), hmolar, hmass;

  tau = 1.0 / T;
  rho_molar = 1.0 / Molar_mass * rho;
  taustar = Tcrit / T;

  dar_dDelta = psi_minus(rho_molar, 0, 1) - tau_times_a(tau, 0) / Ru * psi_plus(rho_molar, 1);

  da0_dTau = Get_da0_dTau(taustar);
  da0_dTau /= 1.0 / Tcrit;

  dar_dTau = psi_minus(rho_molar, 1, 0) - tau_times_a(tau, 1) / Ru * psi_plus(rho_molar, 0);
  /* molar enthalpy */
  hmolar = Ru * T * (1 + tau * (da0_dTau + dar_dTau) + rho_molar * dar_dDelta);
  /* mass enthalpy */
  hmass = hmolar / Molar_mass;
  return hmass;
}

su2double CPengRobinsonCoolProp::GetSmassByRhoT(su2double rho, su2double T) {
  su2double tau, rho_molar, taustar, deltastar;
  su2double da0_dTau(0.0), dar_dTau(0.0), alpha0(0.0), alphar(0.0), smolar(0.0), smass(0.0);
  su2double psi_plus0;

  tau = 1.0 / T;
  rho_molar = 1.0 / Molar_mass * rho;
  taustar = Tcrit / T;
  deltastar = 1.0 / Rhocrit_molar * rho_molar;

  psi_plus0 = psi_plus(rho_molar, 0);

  HelmholtzDerivatives ders = CalcAllHelmholtzIdeal(taustar, deltastar);
  alpha0 = ders.alpha0;
  da0_dTau = ders.dalpha0_dtau;
  da0_dTau /= 1.0 / Tcrit;

  alphar = psi_minus(rho_molar, 0, 0) - tau_times_a(tau, 0) / Ru * psi_plus0;
  dar_dTau = psi_minus(rho_molar, 1, 0) - tau_times_a(tau, 1) / Ru * psi_plus0;

  /* molar entropy */
  smolar = Ru * (tau * (da0_dTau + dar_dTau) - alpha0 - alphar);
  smass = smolar / Molar_mass;

  return smass;
}

su2double CPengRobinsonCoolProp::psi_minus(su2double delta, std::size_t itau, std::size_t idelta) {
  if (itau > 0) return 0.0;
  su2double cm_term(0.0);
  su2double bmc = bm_term() - cm_term; /* appears only in the form(b - c) in the equations */
  su2double bracket = 1 - bmc * delta;

  switch (idelta) {
    case 0:
      return -log(bracket);
    case 1:
      return bmc / bracket;
    case 2:
      return pow(bmc / bracket, 2.0);
    case 3:
      return 2.0 * pow(bmc / bracket, 3.0);
    case 4:
      return 6.0 * pow(bmc / bracket, 4.0);
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}
su2double CPengRobinsonCoolProp::bm_term() {
  su2double summer = 0;

  summer = 0.07780 * Ru * Tcrit / Pcrit;
  return summer;
}
su2double CPengRobinsonCoolProp::A_term(su2double delta) {
  su2double bm = bm_term();
  su2double cm = 0.0;
  return log((delta * (Delta_1 * bm + cm) + 1.0) / (delta * (Delta_2 * bm + cm) + 1.0));
}
su2double CPengRobinsonCoolProp::psi_plus(su2double delta, std::size_t idelta) {
  su2double c_term = 1.0 / bm_term();

  switch (idelta) {
    case 0:
      return A_term(delta) * c_term / (Delta_1 - Delta_2);
    case 1:
      return 1.0 / PI_12(delta, 0);
    case 2:
      return -1.0 / pow(PI_12(delta, 0), 2.0) * PI_12(delta, 1);
    case 3:
      return (-PI_12(delta, 0) * PI_12(delta, 2) + 2.0 * pow(PI_12(delta, 1), 2.0)) / pow(PI_12(delta, 0), 3.0);
    case 4:
      /* Term -PI_12(delta,x,0)*PI_12(delta,x,3) in the numerator is zero (and removed) since PI_12(delta,x,3) = 0*/
      return (6.0 * PI_12(delta, 0) * PI_12(delta, 1) * PI_12(delta, 2.0) - 6.0 * pow(PI_12(delta, 1), 3.0)) /
             pow(PI_12(delta, 0), 4.0);
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}
su2double CPengRobinsonCoolProp::PI_12(su2double delta, std::size_t idelta) {
  su2double bm = bm_term();
  su2double cm = 0.0;

  switch (idelta) {
    case 0:
      return (1.0 + (Delta_1 * bm + cm) * delta) * (1.0 + (Delta_2 * bm + cm) * delta);
    case 1:
      return (2.0 * (bm * Delta_1 + cm) * (bm * Delta_2 + cm) * delta + (Delta_1 + Delta_2) * bm + 2.0 * cm);
    case 2:
      return 2.0 * (Delta_1 * bm + cm) * (Delta_2 * bm + cm);
    case 3:
      return 0.0;
    case 4:
      return 0.0;
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}
su2double CPengRobinsonCoolProp::tau_times_a(su2double tau, std::size_t itau) {
  su2double res(0.0);
  if (itau == 0) {
    res = tau * am_term(tau, 0);
  } else {
    res = tau * am_term(tau, itau) + itau * am_term(tau, itau - 1);
  }
  return res;
}
su2double CPengRobinsonCoolProp::am_term(su2double tau, std::size_t itau) {
  su2double summer = aij_term(tau, 0, 0, itau);
  return summer;
}
su2double CPengRobinsonCoolProp::aij_term(su2double& tau, std::size_t i, std::size_t j, std::size_t& itau) {
  su2double u = u_term(tau, i, j, 0);
  switch (itau) {
    case 0:
      return sqrt(u);
    case 1:
      return 1.0 / (2.0 * sqrt(u)) * u_term(tau, i, j, 1);
    case 2:
      return 1.0 / (4.0 * pow(u, 3.0 / 2.0)) * (2.0 * u * u_term(tau, i, j, 2) - pow(u_term(tau, i, j, 1), 2.0));
    case 3:
      return 1.0 / (8.0 * pow(u, 5.0 / 2.0)) *
             (4.0 * pow(u, 2.0) * u_term(tau, i, j, 3) - 6.0 * u * u_term(tau, i, j, 1) * u_term(tau, i, j, 2) +
              3.0 * pow(u_term(tau, i, j, 1), 3.0));
    case 4:
      return 1.0 / (16.0 * pow(u, 7.0 / 2.0)) *
             (-4.0 * pow(u, 2.0) *
                  (4.0 * u_term(tau, i, j, 1) * u_term(tau, i, j, 3) + 3.0 * pow(u_term(tau, i, j, 2), 2.0)) +
              8.0 * pow(u, 3.0) * u_term(tau, i, j, 4) +
              36.0 * u * pow(u_term(tau, i, j, 1), 2.0) * u_term(tau, i, j, 2) - 15.0 * pow(u_term(tau, i, j, 1), 4.0));
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}
su2double CPengRobinsonCoolProp::u_term(su2double& tau, std::size_t& i, std::size_t& j, const std::size_t& itau) {
  su2double aii = aii_term(tau, i, 0), ajj = aii_term(tau, j, 0);
  switch (itau) {
    case 0:
      return aii * ajj;
    case 1:
      return aii * aii_term(tau, j, 1) + ajj * aii_term(tau, i, 1);
    case 2:
      return (aii * aii_term(tau, j, 2) + 2.0 * aii_term(tau, i, 1) * aii_term(tau, j, 1) + ajj * aii_term(tau, i, 2));
    case 3:
      return (aii * aii_term(tau, j, 3) + 3.0 * aii_term(tau, i, 1) * aii_term(tau, j, 2) +
              3.0 * aii_term(tau, i, 2) * aii_term(tau, j, 1) + ajj * aii_term(tau, i, 3));
    case 4:
      return (aii * aii_term(tau, j, 4) + 4.0 * aii_term(tau, i, 1) * aii_term(tau, j, 3) +
              6.0 * aii_term(tau, i, 2) * aii_term(tau, j, 2) + 4.0 * aii_term(tau, i, 3) * aii_term(tau, j, 1) +
              ajj * aii_term(tau, i, 4));
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}
su2double CPengRobinsonCoolProp::aii_term(su2double& tau, std::size_t& i, const std::size_t& itau) {
  su2double m = k;
  su2double B = 1.0 + m * (1.0 - sqrt_Tr_Tci * sqrt(1.0 / tau));
  switch (itau) {
    case 0:
      return a0 * B * B;
    case 1:
      return a0 * m * B / pow(tau, 3.0 / 2.0) * sqrt_Tr_Tci;
    case 2:
      return a0 * m / 2.0 * (m / pow(tau, 3.0) * Tr_over_Tci - 3.0 * B / pow(tau, 5.0 / 2.0) * sqrt_Tr_Tci);
    case 3:
      return (3.0 / 4.0) * a0 * m *
             (-3.0 * m / pow(tau, 4.0) * Tr_over_Tci + 5.0 * B / pow(tau, 7.0 / 2.0) * sqrt_Tr_Tci);
    case 4:
      return (3.0 / 8.0) * a0 * m *
             (29.0 * m / pow(tau, 5.0) * Tr_over_Tci - 35.0 * B / pow(tau, 9.0 / 2.0) * sqrt_Tr_Tci);
	default:
      SU2_MPI::Error("case argument exception", CURRENT_FUNCTION);
  }
}

su2double CPengRobinsonCoolProp::Get_alpha0(const su2double& taustar, const su2double& deltastar) {
  su2double alpha0{0.0};

  alpha0 += Lead.Get_alpha0(taustar, deltastar);
  alpha0 += EnthEntrOffsetCore.Get_alpha0(taustar);
  alpha0 += LogTau.Get_alpha0(taustar);
  alpha0 += Power.Get_alpha0(taustar);
  alpha0 += PlanckEinstein.Get_alpha0(taustar);
  alpha0 += CP0PolyT.Get_alpha0(taustar);
  alpha0 += CP0Constant.Get_alpha0(taustar);
  return alpha0;
}

su2double CPengRobinsonCoolProp::Get_da0_dTau(const su2double& taustar) {
  su2double da0_dTau{0.0};

  da0_dTau += Lead.Get_dalpha0_dtau();
  da0_dTau += EnthEntrOffsetCore.Get_dalpha0_dtau();
  da0_dTau += LogTau.Get_dalpha0_dtau(taustar);
  da0_dTau += Power.Get_dalpha0_dtau(taustar);
  da0_dTau += PlanckEinstein.Get_dalpha0_dtau(taustar);
  da0_dTau += CP0PolyT.Get_dalpha0_dtau(taustar);
  da0_dTau += CP0Constant.Get_dalpha0_dtau(taustar);
  return da0_dTau;
}

su2double CPengRobinsonCoolProp::Get_d2a0_dTau2(const su2double& taustar) {
  su2double d2a0_dTau2{0.0};

  d2a0_dTau2 += LogTau.Get_d2alpha0_dtau2(taustar);
  d2a0_dTau2 += Power.Get_d2alpha0_dtau2(taustar);
  d2a0_dTau2 += PlanckEinstein.Get_d2alpha0_dtau2(taustar);
  d2a0_dTau2 += CP0PolyT.Get_d2alpha0_dtau2(taustar);
  d2a0_dTau2 += CP0Constant.Get_d2alpha0_dtau2(taustar);
  return d2a0_dTau2;
}

HelmholtzDerivatives CPengRobinsonCoolProp::CalcAllHelmholtzIdeal(const su2double& taustar,
                                                                  const su2double& deltastar) {
  HelmholtzDerivatives derivs;  // zeros out the elements

  Lead.SolveAll(taustar, deltastar, derivs);
  EnthEntrOffsetCore.SolveAll(taustar, deltastar, derivs);
  LogTau.SolveAll(taustar, deltastar, derivs);
  Power.SolveAll(taustar, deltastar, derivs);
  PlanckEinstein.SolveAll(taustar, deltastar, derivs);
  CP0Constant.SolveAll(taustar, deltastar, derivs);
  CP0PolyT.SolveAll(taustar, deltastar, derivs);

  return derivs;
}

void CPengRobinsonCoolProp::Rho_Tp_cubic(su2double T, su2double P, int& Nsolns, su2double& rho0, su2double& rho1, su2double& rho2) {
	
  su2double am(0.0), bm(0.0), cm(0.0), d1(0.0), d2(0.0), d3(0.0), crho0(0.0),crho1(0.0),crho2(0.0),crho3(0.0),DELTA(0.0), p(0.0), q(0.0);
  su2double t0(0.0), t1(0.0), t2(0.0);
  
  am = am_term(1 / T, 0);
  bm = bm_term();
  cm = 0.0;

   /*Introducing new variables to simplify the equation:*/
  d1 = cm - bm;
  d2 = cm + Delta_1 * bm;
  d3 = cm + Delta_2 * bm;

   /*Cubic coefficients:*/
  crho0 = -P;
  crho1 = Ru * T - P * (d1 + d2 + d3);
  crho2 = Ru * T * (d2 + d3) - P * (d1 * (d2 + d3) + d2 * d3) - am;
  crho3 = Ru * T * d2 * d3 - P * d1 * d2 * d3 - d1 * am;

  /*Discriminant*/
  DELTA = 18 * crho3 * crho2 * crho1 * crho0 - 4 * crho2 * crho2 * crho2 * crho0 +
                    crho2 * crho2 * crho1 * crho1 - 4 * crho3 * crho1 * crho1 * crho1 -
                    27 * crho3 * crho3 * crho0 * crho0;
   /*Coefficients for the depressed cubic t^3+p*t+q = 0*/
  p = (3 * crho3 * crho1 - crho2 * crho2) / (3 * crho3 * crho3);
  q = (2 * crho2 * crho2 * crho2 - 9 * crho3 * crho2 * crho1 + 27 * crho3 * crho3 * crho0) /
                (27 * crho3 * crho3 * crho3);

  if (DELTA < 0) {
    /* One real root*/
    if (4 * p * p * p + 27 * q * q > 0 && p < 0) {
      t0 = -2.0 * abs(q) / q * sqrt(-p / 3.0) * cosh(1.0 / 3.0 * acosh(-3.0 * abs(q) / (2.0 * p) * sqrt(-3.0 / p)));
    } else {
      t0 = -2.0 * sqrt(p / 3.0) * sinh(1.0 / 3.0 * asinh(3.0 * q / (2.0 * p) * sqrt(3.0 / p)));
    }
    Nsolns = 1;
    rho0 = t0 - crho2 / (3 * crho3);
    rho1 = t0 - crho2 / (3 * crho3);
    rho2 = t0 - crho2 / (3 * crho3);
  } else  /*(DELTA>0)*/
  {
     /*Three real roots*/
    t0 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 0 * 2.0 * PI_NUMBER / 3.0);
    t1 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 1 * 2.0 * PI_NUMBER / 3.0);
    t2 = 2.0 * sqrt(-p / 3.0) * cos(1.0 / 3.0 * acos(3.0 * q / (2.0 * p) * sqrt(-3.0 / p)) - 2 * 2.0 * PI_NUMBER / 3.0);

    Nsolns = 3;
    rho0 = t0 - crho2 / (3 * crho3);
    rho1 = t1 - crho2 / (3 * crho3);
    rho2 = t2 - crho2 / (3 * crho3);
  }

  rho0 *= Molar_mass;
  rho1 *= Molar_mass;
  rho2 *= Molar_mass;
  sort3(rho0, rho1, rho2);
  return;
}