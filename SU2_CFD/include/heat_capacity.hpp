/*!
 * \file heat_capacity.hpp
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

#ifndef HEAT_CAPACITY_HPP_
#define HEAT_CAPACITY_HPP_
#endif /* HEAT_CAPACITY_HPP_ */
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
 * \class CHeatCapacity
 * \version 1.0
 */
class CHeatCapacity {

protected: su2double Gas_Constant, Gas_ConstantND, Cv0, Gamma, *coeff_Cp0, Tref;
           bool Constant_Gamma;
CConfig *config;


public:

     CHeatCapacity();
    /*!
     * \brief Constructor of the class.
     */
     CHeatCapacity(CConfig *config);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CHeatCapacity(void);

    /*!
     * \brief return liquid density value.
     */
    su2double Get_Cv0(void);

    /*!
     * \brief return liquid density value.
     */
    virtual void Set_Cv0 (su2double T);
};


class CHeatCapacity_Dimensionless : public CHeatCapacity {

public:

     CHeatCapacity_Dimensionless();

     CHeatCapacity_Dimensionless(CConfig *config);
    /*!
     * \brief Destructor of the class.
     */
     ~CHeatCapacity_Dimensionless(void);

    /*!
     * \brief return liquid density value.
     */
    void Set_Cv0 (su2double T);
};

class CHeatCapacity_Dimensional : public CHeatCapacity {

public:

     CHeatCapacity_Dimensional();

     CHeatCapacity_Dimensional(CConfig *config);
    /*!
     * \brief Destructor of the class.
     */
    ~CHeatCapacity_Dimensional(void);

    /*!
     * \brief return liquid density value.
     */
    void Set_Cv0 (su2double T);
};



#include "heat_capacity.inl"
