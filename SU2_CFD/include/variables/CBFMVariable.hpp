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
    using MatrixType = C2DContainer<unsigned long, su2double, StorageType::RowMajor,    64, DynamicSize, DynamicSize>;
    protected:
    unsigned short n_BFM_parameters{};
    MatrixType Geometric_Parameters;
    vector<vector<su2double>> Body_Force_Vector;
    vector<vector<su2double>> Relative_Velocity, proj_vector_axial, proj_vector_tangential, proj_vector_radial;
    vector<unsigned short> Body_Force_Factor{};
    MatrixType Blockage_Gradient;
    public:
    
    CBFMVariable(unsigned long nPoints, unsigned short nDim, unsigned short nBFMParams);
    ~CBFMVariable() override = default;

    void SetBodyForceFactor(unsigned long iPoint, unsigned short value){Body_Force_Factor.at(iPoint) = value;};
    inline unsigned short GetBodyForceFactor(unsigned long iPoint){return Body_Force_Factor.at(iPoint);};

    void SizeParameters(unsigned long nPoints, unsigned short nDim);

    su2double GetBodyForce(unsigned long iPoint, unsigned short iDim){return Body_Force_Vector.at(iPoint).at(iDim);};

    void SetBodyForce(unsigned long iPoint, unsigned short iDim, su2double value){Body_Force_Vector.at(iPoint).at(iDim) = value;}

    su2double GetRelativeVelocity(unsigned long iPoint, unsigned short iDim){return Relative_Velocity.at(iPoint).at(iDim);}

    void SetRelativeVelocity(unsigned long iPoint, unsigned short iDim, su2double value){Relative_Velocity.at(iPoint).at(iDim)= value;}

    su2double GetAxialProjection(unsigned long iPoint, unsigned short iDim){return proj_vector_axial.at(iPoint).at(iDim);};
    su2double GetTangentialProjection(unsigned long iPoint, unsigned short iDim){return proj_vector_tangential.at(iPoint).at(iDim);};
    su2double GetRadialProjection(unsigned long iPoint, unsigned short iDim){return proj_vector_radial.at(iPoint).at(iDim);};
    
    void SetAxialProjection(unsigned long iPoint, unsigned short iDim, su2double value){proj_vector_axial.at(iPoint).at(iDim)= value;}
    void SetTangentialProjection(unsigned long iPoint, unsigned short iDim, su2double value){proj_vector_tangential.at(iPoint).at(iDim) = value;}
    void SetRadialProjection(unsigned long iPoint, unsigned short iDim, su2double value){proj_vector_radial.at(iPoint).at(iDim)= value;}
    
    su2double GetBlockageGradient(unsigned long iPoint, unsigned short iDim){return GetAuxVarGradient(iPoint, I_BLOCKAGE_FACTOR, iDim);}
    void SetBlockageGradient(unsigned long iPoint, unsigned short iDim, su2double value){Blockage_Gradient(iPoint, iDim) = value;}

};