/*!
 * \file CBFMVariable.hpp
 * \brief Main class for defining the variables of the BFM solver.
 * \author E.C. Bunschoten
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "CVariable.hpp"

/*!
 * \class CBFMVariable
 * \brief Main class for defining the variables of the Body-Force Model solver.
 * \author E.C. Bunschoten
 */
class CBFMVariable : public CVariable {
    protected:
    unsigned short n_BFM_parameters{};
    vector<vector<su2double>> geometric_parameters;
    vector<vector<su2double>> Body_Force_Vector{};

    public:
    
    CBFMVariable();
    ~CBFMVariable();

    void SizeParameters(unsigned long nPoints, unsigned short nDim);

};