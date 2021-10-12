/*!
 * \file CPengRobinsonCoolProp.hpp
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

#pragma once

#include "CFluidModel.hpp"

#define LIST_OF_DERIVATIVE_VARIABLES \
  X(alpha0)                          \
  X(dalpha0_dtau)                    \
  X(d2alpha0_dtau2)

struct HelmholtzDerivatives {
#define X(name) su2double name;
  LIST_OF_DERIVATIVE_VARIABLES
#undef X

  void reset(su2double v) {
#define X(name) name = v;
    LIST_OF_DERIVATIVE_VARIABLES
#undef X
  }
  HelmholtzDerivatives() { reset(0.0); };
};
#undef LIST_OF_DERIVATIVE_VARIABLES


class IdealHelmholtzLead  {
 private:
  su2double _a1, _a2;
  bool enabled;

 public:
  IdealHelmholtzLead() : _a1(DBL_MAX), _a2(DBL_MAX), enabled(false) {}

  IdealHelmholtzLead(su2double a1, su2double a2) : _a1(a1), _a2(a2), enabled(true) {}

  bool is_enabled() const { return enabled; }

  inline su2double Get_dalpha0_dtau() const {
    if (!enabled) {
      return 0.0;
    }
    return _a2;
  }
  inline su2double Get_alpha0(const su2double& taustar, const su2double& deltastar) const {
    if (!enabled) {
      return 0.0;
    }
    return log(deltastar) + _a1 + _a2 * taustar;
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    derivs.alpha0 += log(deltastar) + _a1 + _a2 * taustar;
    derivs.dalpha0_dtau += _a2;
  }
};

class IdealHelmholtzEnthalpyEntropyOffset {
 private:
  su2double _a1, _a2;  
  std::string reference;
  bool enabled;

 public:
  IdealHelmholtzEnthalpyEntropyOffset() : _a1(DBL_MAX), _a2(DBL_MAX), enabled(false) {}

  IdealHelmholtzEnthalpyEntropyOffset(su2double a1, su2double a2, const std::string& ref)
      : _a1(a1), _a2(a2), reference(ref), enabled(true) {}

  bool is_enabled() const { return enabled; };

  inline su2double Get_dalpha0_dtau() const {
    if (!enabled) {
      return 0.0;
    }
    return _a2;
  }
  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return _a1 + _a2 * taustar;
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    derivs.alpha0 += _a1 + _a2 * taustar;
    derivs.dalpha0_dtau += _a2;
  }
};

class IdealHelmholtzLogTau {
 private:
  su2double _a1;
  bool enabled;

 public:
  
  IdealHelmholtzLogTau() : _a1(DBL_MAX), enabled(false) {}


  IdealHelmholtzLogTau(su2double a1) : _a1(a1), enabled(true) {}

  bool is_enabled() const { return enabled; };

  inline su2double Get_dalpha0_dtau(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return _a1 / taustar;
  }
  inline su2double Get_d2alpha0_dtau2(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return -_a1 / taustar / taustar;
  }
  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return _a1 * log(taustar);
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    derivs.alpha0 += _a1 * log(taustar);
    derivs.dalpha0_dtau += _a1 / taustar;
    derivs.d2alpha0_dtau2 += -_a1 / taustar / taustar;
  }
};

class IdealHelmholtzPower {
 private:
  std::vector<su2double> n, t;  
  std::size_t N;
  bool enabled;

 public:
  IdealHelmholtzPower() : N(0), enabled(false){};
  // Constructor
  IdealHelmholtzPower(const std::vector<su2double>& n, const std::vector<su2double>& t)
      : n(n), t(t), N(n.size()), enabled(true){};

  bool is_enabled() const { return enabled; };

  inline su2double Get_dalpha0_dtau(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * t[i] * pow(taustar, t[i] - 1);
    }
    return s;
  }
  inline su2double Get_d2alpha0_dtau2(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * t[i] * (t[i] - 1) * pow(taustar, t[i] - 2);
    }
    return s;
  }
  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * pow(taustar, t[i]);
    }
    return s;
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    {
      su2double s = 0;
      for (std::size_t i = 0; i < N; ++i) {
        s += n[i] * pow(taustar, t[i]);
      }
      derivs.alpha0 += s;
    }
    {
      su2double s = 0;
      for (std::size_t i = 0; i < N; ++i) {
        s += n[i] * t[i] * pow(taustar, t[i] - 1);
      }
      derivs.dalpha0_dtau += s;
    }
    {
      su2double s = 0;
      for (std::size_t i = 0; i < N; ++i) {
        s += n[i] * t[i] * (t[i] - 1) * pow(taustar, t[i] - 2);
      }
      derivs.d2alpha0_dtau2 += s;
    }
  }
};

class IdealHelmholtzPlanckEinsteinGeneralized {
 private:
  std::vector<su2double> n, theta, c, d; 
  std::size_t N;
  bool enabled;

 public:
  IdealHelmholtzPlanckEinsteinGeneralized() : N(0), enabled(false) {}
  
  IdealHelmholtzPlanckEinsteinGeneralized(const std::vector<su2double>& n, const std::vector<su2double>& theta,
                                          const std::vector<su2double>& c, const std::vector<su2double>& d)
      : n(n), theta(theta), c(c), d(d), N(n.size()), enabled(true) {}

  bool is_enabled() const { return enabled; };

   /*Extend the vectors to allow for multiple instances feeding values to this function*/
  void extend(const vector<su2double>& n, const vector<su2double>& theta, const vector<su2double>& c,
              const vector<su2double>& d) {
    this->n.insert(this->n.end(), n.begin(), n.end());
    this->theta.insert(this->theta.end(), theta.begin(), theta.end());
    this->c.insert(this->c.end(), c.begin(), c.end());
    this->d.insert(this->d.end(), d.begin(), d.end());
    N += n.size();
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    {
      std::vector<su2double> expthetatau(N);
      for (size_t i = 0; i < N; ++i) {
        expthetatau[i] = exp(theta[i] * taustar);
      }
      {
        su2double s = 0;
        for (size_t i = 0; i < N; ++i) {
          s += n[i] * log(c[i] + d[i] * expthetatau[i]);
        }
        derivs.alpha0 += s;
      }
      {
        su2double s = 0;
        for (size_t i = 0; i < N; ++i) {
          s += n[i] * theta[i] * d[i] * expthetatau[i] / (c[i] + d[i] * expthetatau[i]);
        }
        derivs.dalpha0_dtau += s;
      }
      {
        su2double s = 0;
        for (size_t i = 0; i < N; ++i) {
          s += n[i] * pow(theta[i], 2) * c[i] * d[i] * expthetatau[i] / pow(c[i] + d[i] * expthetatau[i], 2);
        }
        derivs.d2alpha0_dtau2 += s;
      }
    }
  }

  inline su2double Get_dalpha0_dtau(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    std::vector<su2double> expthetatau(N);
    for (size_t i = 0; i < N; ++i) {
      expthetatau[i] = exp(theta[i] * taustar);
    }
    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * theta[i] * d[i] * expthetatau[i] / (c[i] + d[i] * expthetatau[i]);
    }
    return s;
  }

  inline su2double Get_d2alpha0_dtau2(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    std::vector<su2double> expthetatau(N);
    for (size_t i = 0; i < N; ++i) {
      expthetatau[i] = exp(theta[i] * taustar);
    }
    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * pow(theta[i], 2.0) * c[i] * d[i] * expthetatau[i] / pow(c[i] + d[i] * expthetatau[i], 2);
    }
    return s;
  }

  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    std::vector<su2double> expthetatau(N);
    for (size_t i = 0; i < N; ++i) {
      expthetatau[i] = exp(theta[i] * taustar);
    }
    su2double s = 0;
    for (std::size_t i = 0; i < N; ++i) {
      s += n[i] * log(c[i] + d[i] * expthetatau[i]);
    }
    return s;
  }
};

class IdealHelmholtzCP0PolyT {
 private:
  std::vector<su2double> c, t;
  su2double Tc, T0, tau0;  
  std::size_t N;
  bool enabled;

 public:
  IdealHelmholtzCP0PolyT() : Tc(DBL_MAX), T0(DBL_MAX), tau0(DBL_MAX), N(0), enabled(false) {}

  IdealHelmholtzCP0PolyT(const std::vector<su2double>& c, const std::vector<su2double>& t, su2double Tc, su2double T0)
      : c(c), t(t), Tc(Tc), T0(T0), tau0(Tc / T0), N(c.size()), enabled(true) {}

  bool is_enabled() const { return enabled; };

  void extend(const vector<su2double>& c, const vector<su2double>& t) {
    this->c.insert(this->c.end(), c.begin(), c.end());
    this->t.insert(this->t.end(), t.begin(), t.end());
    N += c.size();
  }

  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
      if (abs(t[i]) < 10 * EPS) {
        sum += c[i] - c[i] * taustar / tau0 + c[i] * log(taustar / tau0);
      } else if (abs(t[i] + 1) < 10 * EPS) {
        sum += c[i] * taustar / Tc * log(tau0 / taustar) + c[i] / Tc * (taustar - tau0);
      } else {
        sum += -c[i] * pow(Tc, t[i]) * pow(taustar, -t[i]) / (t[i] * (t[i] + 1)) -
               c[i] * pow(T0, t[i] + 1) * taustar / (Tc * (t[i] + 1)) + c[i] * pow(T0, t[i]) / t[i];
      }
    }
    return sum;
  }

  inline su2double Get_dalpha0_dtau(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
      if (abs(t[i]) < 10 * EPS) {
        sum += c[i] / taustar - c[i] / tau0;
      } else if (abs(t[i] + 1) < 10 * EPS) {
        sum += c[i] / Tc * log(tau0 / taustar);
      } else {
        sum += c[i] * pow(Tc, t[i]) * pow(taustar, -t[i] - 1) / (t[i] + 1) -
               c[i] * pow(Tc, t[i]) / (pow(tau0, t[i] + 1) * (t[i] + 1));
      }
    }
    return sum;
  }

  inline su2double Get_d2alpha0_dtau2(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }

    su2double sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
      if (abs(t[i]) < 10 * EPS) {
        sum += -c[i] / (taustar * taustar);
      } else if (abs(t[i] + 1) < 10 * EPS) {
        sum += -c[i] / (taustar * Tc);
      } else {
        sum += -c[i] * pow(Tc / taustar, t[i]) / (taustar * taustar);
      }
    }
    return sum;
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    {
      su2double sum = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (abs(t[i]) < 10 * EPS) {
          sum += c[i] - c[i] * taustar / tau0 + c[i] * log(taustar / tau0);
        } else if (abs(t[i] + 1) < 10 * EPS) {
          sum += c[i] * taustar / Tc * log(tau0 / taustar) + c[i] / Tc * (taustar - tau0);
        } else {
          sum += -c[i] * pow(Tc, t[i]) * pow(taustar, -t[i]) / (t[i] * (t[i] + 1)) -
                 c[i] * pow(T0, t[i] + 1) * taustar / (Tc * (t[i] + 1)) + c[i] * pow(T0, t[i]) / t[i];
        }
      }
      derivs.alpha0 += sum;
    }
    {
      su2double sum = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (abs(t[i]) < 10 * EPS) {
          sum += c[i] / taustar - c[i] / tau0;
        } else if (abs(t[i] + 1) < 10 * EPS) {
          sum += c[i] / Tc * log(tau0 / taustar);
        } else {
          sum += c[i] * pow(Tc, t[i]) * pow(taustar, -t[i] - 1) / (t[i] + 1) -
                 c[i] * pow(Tc, t[i]) / (pow(tau0, t[i] + 1) * (t[i] + 1));
        }
      }
      derivs.dalpha0_dtau += sum;
    }
    {
      su2double sum = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (abs(t[i]) < 10 * EPS) {
          sum += -c[i] / (taustar * taustar);
        } else if (abs(t[i] + 1) < 10 * EPS) {
          sum += -c[i] / (taustar * Tc);
        } else {
          sum += -c[i] * pow(Tc / taustar, t[i]) / (taustar * taustar);
        }
      }
      derivs.d2alpha0_dtau2 += sum;
    }
  }
};

class IdealHelmholtzCP0Constant {
 private:
  su2double cp_over_R, Tc, T0, tau0;  /* Use these variables internally */
  bool enabled;

 public:
  /*  Default constructor */
  IdealHelmholtzCP0Constant() : cp_over_R(DBL_MAX), Tc(DBL_MAX), T0(DBL_MAX), tau0(DBL_MAX), enabled(false){};

  /* Constructor with three double values */
  IdealHelmholtzCP0Constant(su2double cp_over_R, su2double Tc, su2double T0)
      : cp_over_R(cp_over_R), Tc(Tc), T0(T0), enabled(true) {
    tau0 = Tc / T0;
  };

  ~IdealHelmholtzCP0Constant(){};

  bool is_enabled() const { return enabled; };

  inline su2double Get_dalpha0_dtau(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return cp_over_R / taustar - cp_over_R / tau0;
  }
  inline su2double Get_d2alpha0_dtau2(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return -cp_over_R / (taustar * taustar);
  }
  inline su2double Get_alpha0(const su2double& taustar) const {
    if (!enabled) {
      return 0.0;
    }
    return cp_over_R - cp_over_R * taustar / tau0 + cp_over_R * log(taustar / tau0);
  }

  inline void SolveAll(const su2double& taustar, const su2double& deltastar, HelmholtzDerivatives& derivs) {
    if (!enabled) {
      return;
    }
    derivs.alpha0 += cp_over_R - cp_over_R * taustar / tau0 + cp_over_R * log(taustar / tau0);
    derivs.dalpha0_dtau += cp_over_R / taustar - cp_over_R / tau0;
    derivs.d2alpha0_dtau2 += -cp_over_R / (taustar * taustar);
  }
};

class SaturationLinelPs {
 private:
  std::vector<su2double> n, t;
  su2double reducing_value, T_r, Tmax, Tmin;
  size_t N;
  bool using_tau_r;

 public:
  SaturationLinelPs() : T_r(DBL_MAX), Tmax(DBL_MAX), Tmin(DBL_MAX), reducing_value(DBL_MAX), N(0), using_tau_r(false) {}

  /*  Constructor with std::vectors */
  SaturationLinelPs(const std::vector<su2double>& n, const std::vector<su2double>& t, su2double T_r, su2double Tmax,
                    su2double Tmin, su2double reducing_value, bool using_tau_r)
      : n(n),
        t(t),
        T_r(T_r),
        Tmax(Tmax),
        Tmin(Tmin),
        reducing_value(reducing_value),
        N(n.size()),
        using_tau_r(using_tau_r) {}

  inline su2double Evaluate(su2double T) {
    su2double theta, summ{0.0};

    theta = 1 - T / T_r;
    for (std::size_t i = 0; i < N; ++i) {
      summ += n[i] * pow(theta, t[i]);
    }

    su2double tau_r_value;
    if (using_tau_r)
      tau_r_value = T_r / T;
    else
      tau_r_value = 1.0;

    return reducing_value * exp(tau_r_value * summ);
  }
};


class CPengRobinsonCoolProp final : public CFluidModel {
 private:
  su2double a{0.0};   /*!< \brief Peng-Robinson model parameter. */
  su2double a0{0.0};  /*!< \brief Peng-Robinson model parameter. */
  su2double b{0.0};   /*!< \brief Peng-Robinson model parameter. */
  su2double k{0.0};   /*!< \brief Peng-Robinson model parameter (computed with acentric factor). */
  su2double Zed{0.0}; /*!< \brief compressibility factor. */
  su2double acentric{0.0};

  su2double Tr_over_Tci{0.0}; /*!< \brief Specific CoolProp cubic model parameter.*/
  su2double sqrt_Tr_Tci{0.0}; /*!< \brief Specific CoolProp cubic model parameter.*/

  SaturationLinelPs PsContainer;                            /*!< \brief Saturation pressure solver object */
  IdealHelmholtzLead Lead;                                  /*!< \brief Lead term of the Helmholtz ideal part derivatives */
  IdealHelmholtzEnthalpyEntropyOffset EnthEntrOffsetCore;   /*!< \brief Enthalpy/Entropy offset term of the Helmholtz ideal part derivatives */
  IdealHelmholtzLogTau LogTau;                              /*!< \brief LogTau term of the Helmholtz ideal part derivatives */
  IdealHelmholtzPower Power;                                /*!< \brief Power term of the Helmholtz ideal part derivatives */
  IdealHelmholtzPlanckEinsteinGeneralized PlanckEinstein;   /*!< \brief Planck-Einstein term of the Helmholtz ideal part derivatives */
  IdealHelmholtzCP0Constant CP0Constant;                    /*!< \brief CP0Constant term of the Helmholtz ideal part derivatives */
  IdealHelmholtzCP0PolyT CP0PolyT;                          /*!< \brief CP0PolyT term of the Helmholtz ideal part derivatives */

  su2double Tcrit{0.0};             /*!< \brief Critical temperature. */
  su2double Ttriple{0.0};           /*!< \brief Temperature of the vapor at the triple point */
  su2double Rhotriple{0.0};         /*!< \brief Density of the vapor at the triple point */
  su2double Pcrit{0.0};             /*!< \brief Critical pressure. */
  su2double Rhocrit_molar{0.0};     /*!< \brief Critical molar density. */
  su2double Rgas{0.0};              /*!< \brief Specific gas constant */
  const su2double Ru{UNIVERSAL_GAS_CONSTANT}; /*!< \brief Universal gas constant */
  su2double Molar_mass{0.0};        /*!< \brief Molar mass */

  const su2double Delta_1 = 1 + sqrt(2.0);
  const su2double Delta_2 = 1 - sqrt(2.0);

  /*
  Set of supplementary functions utilized in calculation of Helmholtz real part derivatives.
  */
  su2double bm_term();
  su2double A_term(su2double delta);
  su2double psi_plus(su2double delta, std::size_t idelta);
  su2double tau_times_a(su2double tau, std::size_t itau);
  su2double am_term(su2double tau, std::size_t itau);
  su2double aij_term(su2double& tau, std::size_t i, std::size_t j, std::size_t& itau);
  su2double u_term(su2double& tau, std::size_t& i, std::size_t& j, const std::size_t& itau);
  su2double aii_term(su2double& tau, std::size_t& i, const std::size_t& itau);
  su2double psi_minus(su2double delta, std::size_t itau, std::size_t idelta);
  su2double PI_12(su2double delta, std::size_t idelta);

  /*!
   * \brief returns Temperature by using Density and internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] umass - second thermodynamic variable.
   */
  su2double GetTemperatureByRhoUmass(su2double rho, su2double umass);

  /*!
   * \brief returns Temperature by using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  su2double GetTemperatureByPRho(su2double P, su2double rho);
  /*!
   * \brief returns speicific internal Energy by using Density and Temperature
   * \param[in] rho - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  su2double GetUmassByRhoT(su2double rho, su2double T);
  /*!
   * \brief returns Pressure by using Density and Temperature
   * \param[in] rho - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  su2double GetPressureByRhoT(su2double rho, su2double T);

  /*!
   * \brief returns specific Enthalpy by using Density and Temperature
   * \param[in] rho - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  su2double GetHmassByRhoT(su2double rho, su2double T);

  /*!
   * \brief returns specific Entropy by using Density and Temperature
   * \param[in] rho - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  su2double GetSmassByRhoT(su2double rho, su2double T);

  /*!
   * \brief set thermodynamic state by using Density and Temperature
   * \param[in] rho - first thermodynamic variable (rho).
   * \param[in] T - second thermodynamic variable (T).
   * \param[in] isSetPressure - boolean flag whether to set Pressure
   * \param[in] isSetEntropy - boolean flag whether to set specific Entropy
   * \param[in] isSetEnergy - boolean flag whether to set internal Energy
   *
   */
  void SetInnerTDState_rhoT(su2double rho, su2double T, bool isSetPressure, bool isSetEntropy, bool isSetEnergy);

   /*!
   * \brief returns saturation Pressure by using saturation Temperature
   * \param[in] T - Saturation temperature.
   */
  su2double Estimate_Ps_SatLine(su2double T);

   /*!
   * \brief returns saturation Temperature by using saturation Pressure
   * \param[in] P - Saturation Pressure.
   */
  su2double Get_Saturation_Temperature(su2double P);

  /*
  * !
  * \brief solves cubic eqution for Density by using Temperature and Pressure  
  * \param[in] T - first thermodynamic variable (Temperature).
  * \param[in] P - second thermodynamic variable (Pressure).
  * \param[out] Nsolns - number of real roots
  * \param[out] rho0 - first root
  * \param[out] rho1 - second root
  * \param[out] rho2 - third root
  */
  void Rho_Tp_cubic(su2double T, su2double P, int& Nsolns, su2double& rho0, su2double& rho1, su2double& rho2);

  inline su2double acosh(su2double x) { return log(x + sqrt(x * x - 1.0)); }
  inline su2double asinh(su2double value) {
    if (value > 0) {
      return log(value + sqrt(value * value + 1));
    } else {
      return -log(-value + sqrt(value * value + 1));
    }
  }

  template <typename T>
  inline void sort3(T& a, T& b, T& c) {
    if (a > b) {
      std::swap(a, b);
    }
    if (a > c) {
      std::swap(a, c);
    }
    if (b > c) {
      std::swap(b, c);
    }
  }

  su2double Get_da0_dTau(const su2double& taustar);
  su2double Get_d2a0_dTau2(const su2double& taustar);
  su2double Get_alpha0(const su2double& taustar, const su2double& deltastar);

  /*
  * \brief Solves for all Helmholtz ideal part derivatives in one-shot
  */
  HelmholtzDerivatives CalcAllHelmholtzIdeal(const su2double& taustar, const su2double& deltastar);

 public:
  /*!
   * \brief Constructor of the class
   */
  CPengRobinsonCoolProp(std::string fluidName);
  /*
  * \brief Destructor of the class
  */
  ~CPengRobinsonCoolProp();

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe(su2double rho, su2double e) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT(su2double P, su2double T) override;

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho(su2double P, su2double rho) override;

  /*!
   * \brief Set the Dimensionless Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho(su2double P, su2double rho) override;

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("hs").
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   *
   */
  void SetTDState_hs(su2double h, su2double s) override;
};