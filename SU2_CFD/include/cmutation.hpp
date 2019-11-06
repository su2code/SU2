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
 *  - Prof. Edwin van der Weide's group at the University of Twente.Fmixture
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

#ifndef CMUTATION_HPP_
#define CMUTATION_HPP_
#endif /* CMUTATION_HPP_ */
#pragma once

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

#include "../../Common/include/datatype_structure.hpp"
#include "../../Common/include/config_structure.hpp"

#ifdef HAVE_MUTATIONPP
#include "mutation++.h"
#endif

using namespace std;


/*!
 * \class CReactive
 * \brief Main class for defining the Transport-Physical ModelF
 * a child class for each particular Model (Power law, Puta, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CMutation {
protected:

   
   Mutation::Mixture mix; //!<Mixture object
   unsigned short ns, i;
   su2double  *Cp_ks, E,  gammaFrozen, gammaEquilibrium, Tref, *hs;
   su2double mu, *Ds;
   su2double  a, Density, *Xs;

   vector<su2double> MolarMass, Cv_ks, Energies, Omega, hf, Energies_species,  Temp,  lambda, rhos, temperatures, Ws;


public:

    /*!
     * \brief Constructor of the class.
     */
    CMutation(su2double *composition, su2double nSpecies, Mutation::MixtureOptions options);

    /*!
     * \brief Destructor of the class.
     */
    virtual    ~CMutation(void);

    vector<su2double> Mutation_MolarMass();

    vector<su2double> Mutation_SpeciesDensity(su2double *cs, su2double rho);

    void       Mutation_UpdateMixtureState(su2double* cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_CvModeSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double* Mutation_Get_CpModeSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_GammaFrozen(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_GammaEquilibrium(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_MixtureEnergy(su2double* cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve); 

    vector<su2double> Mutation_Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);
    
    vector<su2double> Mutation_Get_NetProductionRates(su2double* cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_EnergysourceTerm(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_ReferenceTemperature(su2double* cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_EnthalpiesFormation(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double* Mutation_Get_Enthalpies(su2double* cs, su2double rho, su2double T, su2double Tve);

    su2double* Mutation_Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

    vector<su2double> Mutation_Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

    su2double  Mutation_Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

    su2double  Mutation_Get_Density(su2double T, su2double *cs, su2double P);

    
    


};