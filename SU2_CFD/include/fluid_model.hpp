/*!
 * \file fluid_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

using namespace std;

#include "../include/transport_model.hpp"
#include "../include/liquid_phase_model.hpp"
#include "../include/heat_capacity.hpp"
#include "../../Common/include/config_structure.hpp"


/*!
 * \class CFluidModel
 * \brief Main class for defining the Thermo-Physical Model
 * a child class for each particular Model (Ideal-Gas, Van der Waals, etc.)
 * \author: S.Vitale, G.Gori, M.Pini
 * \version 5.0.0 "Raven"
 */
class CFluidModel {
protected:
su2double      StaticEnergy,      /*!< \brief Internal Energy. */
       Entropy,          /*!< \brief Entropy. */
       Density,          /*!< \brief Density. */
       Pressure,         /*!< \brief Pressure. */
       SoundSpeed2,       /*!< \brief SpeedSound. */
       Temperature,      /*!< \brief Temperature. */
       dPdrho_e,         /*!< \brief DpDd_e. */
       dPde_rho,         /*!< \brief DpDe_d. */
       dTdrho_e,         /*!< \brief DTDd_e. */
       dTde_rho,         /*!< \brief DTDe_d. */
             Cp,                    /*!< \brief Specific Heat Capacity at constant pressure. */
			 Cv,
			 Cv0,
       Mu,          /*!< \brief Specific Heat Capacity at constant pressure. */
         dmudrho_T,       /*!< \brief Specific Heat Capacity at constant pressure. */
         dmudT_rho,        /*!< \brief Specific Heat Capacity at constant pressure. */
         Kt,          /*!< \brief Specific Heat Capacity at constant pressure. */
         dktdrho_T,       /*!< \brief Specific Heat Capacity at constant pressure. */
         dktdT_rho;        /*!< \brief Specific Heat Capacity at constant pressure. */

CViscosityModel *LaminarViscosity;            /*!< \brief Laminar Viscosity Model */
CConductivityModel *ThermalConductivity;    /*!< \brief Thermal Conductivity Model */

CLiquidModel *Liquid_Prop;

CHeatCapacity *HeatCapacity;    /*!< \brief Heat Capacityy Model */


public:

  /*!
     * \brief Constructor of the class.
     */
    CFluidModel(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CFluidModel(void);

    /*!
     * \brief Get fluid pressure.
     */
    su2double GetPressure ();

    /*!
     * \brief Get fluid temperature.
     */
    su2double GetTemperature ();

    /*!
     * \brief Get fluid entropy.
     */
    su2double GetEntropy ();

    /*!
     * \brief Get fluid internal energy.
     */
    su2double GetStaticEnergy ();

    /*!
     * \brief Get fluid density.
     */
    su2double GetDensity ();

    /*!
     * \brief Get fluid speed of sound.
     */
    su2double GetSoundSpeed ();

    /*!
     * \brief Get fluid speed of sound squared.
     */
    su2double GetSoundSpeed2 ();

    /*!
     * \brief Get fluid specific heat at constant pressure.
     */
    su2double GetCp ();

    /*!
     * \brief Get fluid dynamic viscosity
     */

    su2double GetLaminarViscosity ();

    /*!
     * \brief Get fluid thermal conductivity
     */

    su2double GetThermalConductivity ();

    /*!
     * \brief Get fluid pressure partial derivative.
     */
    su2double GetdPdrho_e ();

    /*!
     * \brief Get fluid pressure partial derivative.
     */
    su2double GetdPde_rho ();

    /*!
     * \brief Get fluid temperature partial derivative.
     */
    su2double GetdTdrho_e ();

    /*!
     * \brief Get fluid temperature partial derivative.
     */
    su2double GetdTde_rho ();

    /*!
     * \brief Get fluid dynamic viscosity partial derivative.
     */
    su2double Getdmudrho_T ();

    /*!
     * \brief Get fluid dynamic viscosity partial derivative.
     */
    su2double GetdmudT_rho ();

    /*!
     * \brief Get fluid thermal conductivity partial derivative.
     */
    su2double Getdktdrho_T ();

    /*!
     * \brief Get fluid thermal conductivity partial derivative.
     */
    su2double GetdktdT_rho ();

    /*!
     * \brief Set viscosity model.
     */
    void SetLaminarViscosityModel (CConfig *config);

    /*!
     * \brief Set thermal conductivity model.
     */
    void SetThermalConductivityModel (CConfig *config);

    /*!

     * \brief Set liquid model.
     */
    void SetLiquidPhaseModel (CConfig *config);

     /* \brief Set heat capacity model.
     */
    void SetHeatCapacityModel (CConfig *config);


    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("e, rho").
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */

    virtual void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("PT").
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (T).
     */

    virtual void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("Pv").
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (v).
     */

    virtual void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("Pv").
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (v).
     *
     */

    virtual void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("hs").
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    virtual void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    virtual void SetTDState_rhoT (su2double rho, su2double T );


    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("Pv").
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    virtual void SetTDState_Ps (su2double P, su2double s );


    void SetLiquidProp (su2double P, su2double T, su2double rho, su2double h_v, su2double Rcritical, su2double R, su2double mom3);

    su2double GetLiquidDensity ();
    su2double GetMixtureDensity ();
    su2double GetLiquidTemperature ();
    su2double GetLiquidEnthalpy ();
    su2double GetSurfaceTension ();
    su2double GetTsat ();
    su2double GetPsat ();
    su2double GetCriticalRadius ();
    su2double GetRadius ();


    virtual void SetGamma_Trho ();

    virtual su2double GetGamma ();

    void SetHeatCapacityModel_Dimensionless (CConfig *config);

    void SetHeatCapacityModel_Dimensional   (CConfig *config);

    virtual void Set_Cv(su2double T, su2double v);
};


/*!
 * \class CIdealGas
 * \brief Child class for defining ideal gas model.
 * \author: S.Vitale, M.Pini.
 * \version 5.0.0 "Raven"
 */
class CIdealGas : public CFluidModel {

protected:
  su2double Gamma,             /*!< \brief Heat Capacity Ratio. */
          Gamma_Minus_One,       /*!< \brief Heat Capacity Ratio Minus One. */
          Gas_Constant;        /*!< \brief Gas Constant. */


public:

     /*!
     * \brief Constructor of the class.
     */
    CIdealGas(void);

    /*!
     * \brief Constructor of the class.
     */
    CIdealGas(su2double gamma, su2double R);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CIdealGas(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */

    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure  and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */

    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Internal Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless State using Enthalpy and Entropy
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief Set the Dimensionless State using Density and Temperature
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

    /*!
     * \brief Return gamma value.
     */

    su2double GetGamma ();


};


class CIdealGas_Generic : public CIdealGas {

protected:
  su2double Gamma,             /*!< \brief Heat Capacity Ratio. */
          Gamma_Minus_One,       /*!< \brief Heat Capacity Ratio Minus One. */
          Gas_Constant;        /*!< \brief Gas Constant. */


public:

     /*!
     * \brief Constructor of the class.
     */
    CIdealGas_Generic(void);

    /*!
     * \brief Constructor of the class.
     */
    CIdealGas_Generic(su2double gamma, su2double R);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CIdealGas_Generic(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */

    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure  and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */

    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Internal Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless State using Enthalpy and Entropy
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief Set the Dimensionless State using Density and Temperature
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

    /*!
     * \brief Return gamma value.
     */

    void SetGamma_Trho();
};

/*!
 * derived class CVanDerWaalsGas
 * \brief Child class for defining the Van der Waals model.
 * \author: S.Vitale, M.Pini
 * \version 5.0.0 "Raven"
 */
class CVanDerWaalsGas : public CIdealGas {

protected:
  su2double
      a, b, Zed;             /*!< \brief Parameters for the Dimensionless Equation. */


public:

     /*!
     * \brief Constructor of the class.
     */
    CVanDerWaalsGas(void);

    /*!
     * \brief Constructor of the class.
     */
    CVanDerWaalsGas(su2double gamma, su2double R, su2double Pstar, su2double Tstar);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CVanDerWaalsGas(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */
    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */
    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Internal Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */


    void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless state using Enthalpy and Entropy
     * \param[in] h - first thermodynamic variable (h).
     * \param[in] s - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief Set the Dimensionless state using Density and Temperature
     * \param[in] rho - first thermodynamic variable (rho).
     * \param[in] T - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] P - first thermodynamic variable (P).
     * \param[in] s - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );


};



class CVanDerWaalsGas_Generic : public CVanDerWaalsGas {

protected:
  su2double
      a, b, Zed;             /*!< \brief Parameters for the Dimensionless Equation. */
  su2double Cv, Cv0, Cp, Cp0; /* brief auxiliary variables for gamma evaluation*/
  su2double TstarCrit;

public:


     /*!
     * \brief Constructor of the class.
     */
    CVanDerWaalsGas_Generic(void);

    /*!
     * \brief Destructor of the class.
     */
    ~CVanDerWaalsGas_Generic(void);

    CVanDerWaalsGas_Generic (su2double gamma, su2double R, su2double Pstar, su2double Tstar);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */
    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */
    void SetTDState_PT (su2double P, su2double T );

    void SetTDState_Prho (su2double P, su2double rho );
    /*!
     * \brief Set the Dimensionless Internal Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless state using Enthalpy and Entropy
     * \param[in] h - first thermodynamic variable (h).
     * \param[in] s - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief Set the Dimensionless state using Density and Temperature
     * \param[in] rho - first thermodynamic variable (rho).
     * \param[in] T - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] P - first thermodynamic variable (P).
     * \param[in] s - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

    void SetGamma_Trho ();

};


/*!
 * \derived class CPengRobinson
 * \brief Child class for defining the Peng-Robinson model.
 * \author: S.Vitale, G. Gori
 * \version 5.0.0 "Raven"
 */
class CPengRobinson : public CIdealGas {

protected:
  su2double  a,             /*!< \brief model parameter. */
        b,             /*!< \brief model parameter. */
        k,             /*!< \brief model parameter (computed with acentric factor). */
        Zed,             /*!< \brief compressibility factor. */
        TstarCrit;        /*!< \brief Critical temperature. */

private:

       /*!
      * \brief Internal model parameter.
      */
      su2double  alpha2 (su2double T);


     /*!
    * \brief Internal function for the implicit call hs.
    */
    su2double  T_v_h (su2double v, su2double h);
    /*!
    * \brief Internal function for the implicit call Ps.
    */
    su2double T_P_rho(su2double P, su2double rho);



public:

      /*!
     * \brief Constructor of the class.
     */
    CPengRobinson(void);

    /*!
     * \brief Constructor of the class.
     */
    CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CPengRobinson(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */
    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */
    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetEnergy_Prho (su2double P, su2double rho );
    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("hs").
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

};



class CPengRobinson_Generic : public CPengRobinson {

protected:
  su2double  a,             /*!< \brief model parameter. */
        b,             /*!< \brief model parameter. */
        k,             /*!< \brief model parameter (computed with acentric factor). */
        Zed,             /*!< \brief compressibility factor. */
        TstarCrit,        /*!< \brief Critical temperature. */
        Cv0;

private:

  /*!
      * \brief Internal model parameter.
      */
      su2double  alpha2 (su2double T);

      su2double  dalphadT (su2double T);

      su2double  dalpha2dT2 (su2double T);

      su2double  dalpha3dT3 (su2double T);

      su2double  dCvdT (su2double T, su2double rho);

      su2double  dCvdrho (su2double T, su2double rho);

      su2double  datanh(su2double rho);

     /*!
    * \brief Internal function for the implicit call hs.
    */
    su2double  T_v_h (su2double v, su2double h);
    /*!
    * \brief Internal function for the implicit call Ps.
    */
    su2double T_P_rho(su2double P, su2double rho);



public:

      /*!
     * \brief Constructor of the class.
     */
    CPengRobinson_Generic(void);

    /*!
     * \brief Constructor of the class.
     */
    CPengRobinson_Generic(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CPengRobinson_Generic(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */
    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */
    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetEnergy_Prho (su2double P, su2double rho );
    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("hs").
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );

    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

    void SetGamma_Trho ();

    void Set_Cv(su2double T, su2double v);

};

//#ifdef HAVE_FluidProp

/*!
 * \derived class CFluidProp
 * \brief Child class for defining the Peng-Robinson model.
 * \author: T.P. van der Stelt, M. Pini
 * \version 5.0.0 "Raven"
 */
class CFluidProp : public CFluidModel {
    
protected:
    string ThermoLib;     /*!< \brief Sub-library. */
    int nComp;            /*!< \brief Number of components. */
    string* Comp;         /*!< \brief Components. */
    double* Conc;         /*!< \brief Concentrations. */
    bool SinglePhaseOnly; /*!< \brief Single phase only: indicates that no phase equilibria are considered. */
    string TableName;     /*!< \brief Name of look-up table. */
    
    int ErrorLevel;       /*!< \brief Error level diagnostics flag. */

    su2double Gamma;             /*!< \brief Heat Capacity Ratio. */
    su2double Gas_Constant;
private:
    
    bool LuTSwitchedOn;   /*!< \brief LuT indicator. */

public:

    void SwitchLuTOff();
    void SwitchLuTOn();

    /*!
     * \brief Constructor of the class.
     */
    CFluidProp(void);
    
    /*!
     * \brief Constructor of the class.
     */
    CFluidProp(string thermolib, int ncomp, string* comp, double* conc, bool SinglePhaseOnly, 
               string TableName, double T_ref, double P_ref, double rho_ref, int ErrorLevel);
    
    /*!
     * \brief Destructor of the class.
     */
    virtual ~CFluidProp(void);
    
    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */
    void SetTDState_rhoe (double rho, double e );
    
    /*!
     * \brief Set the Dimensionless State using Pressure and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */
    void SetTDState_PT (double P, double T );
    
    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetTDState_Prho (double P, double rho );
    
    /*!
     * \brief Set the Dimensionless Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */
    void SetEnergy_Prho (double P, double rho );
    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("hs").
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (double h, double s );
    
    /*!
     * \brief virtual member that would be different for each gas model implemented
     * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (double rho, double T );

    /*!
     * \brief Set viscosity model.
     */
    void SetLaminarViscosityModel (CConfig *config);

    /*!
     * \brief Set thermal conductivity model.
     */
    void SetThermalConductivityModel (CConfig *config);

    /*!
     * \brief Return gamma value.
     */

    su2double GetGamma ();
};

//#endif


#include "fluid_model.inl"

