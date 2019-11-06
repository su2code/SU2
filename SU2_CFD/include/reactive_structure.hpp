/*!
 * \file transport_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#ifndef REACTIVE_STRUCTURE_HPP_
#define REACTIVE_STRUCTURE_HPP_
#endif /* REACTIVE_STRUCTURE_HPP_ */
#pragma once

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>


#include "stdio.h"
#include "math.h"

#include "../../Common/include/datatype_structure.hpp"
#include "../../Common/include/config_structure.hpp"

#include "cmutation.hpp"


#ifdef HAVE_MUTATIONPP
#include "mutation++.h"
#endif

using namespace std;


/*!
 * \class CViscosityModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CReactive {
protected:

  bool ionization;
  string MixtureFile;
  unsigned short nSpecies, iSpecies;
  CMutation *mutation;
  su2double *comp, *Cp_ks,  *Cp_trs, *Cp_ves, E, gamma, gammaFrozen, gammaEquilibrium, Tref, *hs;
  su2double  mu, *Ds;
  su2double  a, Density, *Xs, OmegaVT, OmegaCV;

  vector<su2double> Ms, Cv_trs, Cv_ks, Energies, Omega, hf, Energies_Species,  Cv_ves, Temp, lambda, Ws;

  
public:

    /*!
     * \brief Constructor of the class.
     */
    CReactive();

    /*!
     * \brief Destructor of the class.
     */
    virtual           ~CReactive();

    virtual void       InitializeMixture(CConfig *config);

    virtual bool       Get_Ionization();

    virtual vector<su2double> Get_MolarMass();

    virtual vector<su2double> Get_CvTraRotSpecies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_CvVibElSpecies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double* Get_CpTraRotSpecies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double* Get_CpVibElSpecies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

    //virtual su2double Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve);
  
    virtual vector<su2double> Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double* Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double* Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual vector<su2double> Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

    virtual su2double  Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

    virtual su2double  Get_Density(su2double T, su2double *Xs, su2double P);

    
};


class CReactiveMutation : public CReactive {
protected:
  
  Mutation::MixtureOptions opt;
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
             CReactiveMutation(string MixFile, string Transport);
  
  virtual    ~CReactiveMutation(void);
  
 
  void       InitializeMixture(CConfig *config);

  bool       Get_Ionization();

  vector<su2double> Get_MolarMass();

  vector<su2double> Get_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);
  
  vector<su2double> Get_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);
  
  su2double* Get_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);
  
  su2double  Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

  //su2double Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve);
  
  vector<su2double> Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve);
  
  su2double* Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

  su2double  Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_Density(su2double T, su2double *Xs, su2double P);

  
};

class CReactiveHardCode : public CReactive {
protected:

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CReactiveHardCode(void);
  
  /*!
   * \brief Constructor of the class.
   * \brief Destructor of the class.
   */
  virtual    ~CReactiveHardCode(void);
  
  void       InitializeMixture(CConfig *config);

  bool       Get_Ionization();

  vector<su2double> Get_MolarMass();

  vector<su2double> Get_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

  //su2double Get_CVEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);
  
  su2double  Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve);
  
  vector<su2double> Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double* Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

  vector<su2double> Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

  su2double  Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Get_Density(su2double T, su2double *Xs, su2double P);

 
};

#include "reactive_structure.inl"
