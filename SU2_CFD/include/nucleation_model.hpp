/*!
 * \file nucleation_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

#ifndef NUCLEATION_MODEL_HPP_
#define NUCLEATION_MODEL_HPP_
#endif /* NUCLEATION_MODEL_HPP_ */
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

using namespace std;


/*!
 * \class CNucleationModel
 * \version 1.0
 */
class CNucleationModel {

protected: su2double J, G;

		   su2double Gamma, Gas_Constant;
		   double    Boltzmann, MolMass;



public:


    /*!
     * \brief Constructor of the class.
     */
    CNucleationModel();

    CNucleationModel(CConfig *config);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CNucleationModel(void);

    /*!
     * \brief return the nucleation rate.
     */
    su2double GetNucleationRate(void);

    /*!
     * \brief return the growth rate.
     */
    su2double GetGrowthRate(void);

    /*!
     * \brief compute the nucleation rate.
     */
    virtual   void SetNucleationRate(su2double P, su2double T, su2double rho,
            su2double h, su2double k, su2double mu, su2double *V_Liquid);

    /*!
     * \brief compute the growth rate.
     */
    virtual   void SetGrowthRate(su2double P, su2double T, su2double rho,
            su2double h, su2double k, su2double mu, su2double *V_Liquid);

};


/*!
 * \class CClassicalTheory
 * \version 1.0
 */
class CClassicalTheory : public CNucleationModel  {
protected:

  su2double Theta, J, G;
  su2double Lambda, Ni, Pr;

public:
  

  /*!
   * \brief Constructor of the class.
   */
  CClassicalTheory();
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CClassicalTheory(void);
  
  /*!
   * \brief compute the nucleation rate.
   */
  void SetNucleationRate(su2double P, su2double T, su2double rho,
          su2double h, su2double k, su2double mu, su2double *V_Liquid);

  /*!
   * \brief compute the growth rate.
   */
  void SetGrowthRate(su2double P, su2double T, su2double rho,
          su2double h, su2double k, su2double mu, su2double *V_Liquid);

};




#include "nucleation_model.inl"
