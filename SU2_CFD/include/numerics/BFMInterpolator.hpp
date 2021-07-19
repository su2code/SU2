/*!
 * \file ReadBFMInput.hpp
 * \brief Declarations of template (empty) numerics classes, these give
 *        an idea of the methods that need to be defined to implement
 *        new schemes in SU2, in practice you should look for a similar
 *        scheme and try to re-use functionality (not by copy-paste).
 * \author F. Palacios, T. Economon
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
#include "ReadBFMInput.hpp"
#include "../solvers/CSolver.hpp"

class BFMInterpolator
{
    private:
    vector<su2double> Rotation_Axis;
    vector<su2double> axial_direction;

    void Interp2D(su2double axis, su2double radius, unsigned long iPoint, ReadBFMInput *reader, CSolver*solver_container);
    bool RC_Inclusion(su2double axis, su2double radius, vector<su2double>ax_cell, vector<su2double>rad_cell);
    su2double DW_average(su2double axis, su2double radius, vector<su2double>ax_cell, vector<su2double>rad_cell, vector<su2double>val_cell);
    
    su2double Vector_Dot_Product(vector<su2double>a, vector<su2double>b){
        if(a.size() != b.size()){
            SU2_MPI::Error("Input vectors are of unequal size", CURRENT_FUNCTION);
        }
        su2double dot_product{0};
        for(size_t i=0; i<a.size(); ++i){
            dot_product += a.at(i) * b.at(i);
        }
        return dot_product;
    };

    su2double Vector_Dot_Product(vector<su2double>a, su2double *b){
        // if(a.size() != sizeof(b)/sizeof(b[0])){
        //     SU2_MPI::Error("Input vectors are of unequal size", CURRENT_FUNCTION);
        // }
        su2double dot_product{0};
        for(size_t i=0; i<a.size(); ++i){
            dot_product += a.at(i)*b[i];
        }
        return dot_product;
    };

    su2double Vector_Dot_Product(su2double *b, vector<su2double> a){
        // if(a.size() != sizeof(b)/sizeof(b[0])){
        //     SU2_MPI::Error("Input vectors are of unequal size", CURRENT_FUNCTION);
        // }
        return Vector_Dot_Product(a, b);
        // su2double dot_product{0};
        // for(size_t i=0; i<a.size(); ++i){
        //     dot_product += a.at(i)*b[i];
        // }
        // return dot_product;
    };
    
    vector<su2double> Vector_Cross_Product(vector<su2double>a, vector<su2double>b){
        vector<su2double> c;
        if(a.size() != b.size()){
            SU2_MPI::Error("Input vectors are of unequal size", CURRENT_FUNCTION);
        }
        if(a.size() > 2){
            c.resize(a.size());
            for(size_t i=0; i<a.size(); ++i){
                c.at(i) = a.at((i+1) % a.size())*b.at((i+2)%b.size()) - a.at((i+2) % a.size())*b.at((i+1) % b.size());
            }
        }else{
            c.resize(1);
            c.at(0) = a.at(0)*b.at(1) - b.at(0)*a.at(1);
        }
        
        return c;
    };

    public:
    BFMInterpolator();
    BFMInterpolator(ReadBFMInput *reader, CSolver *solver_container, CGeometry *geometry, CConfig *config);

    void Interpolate(ReadBFMInput *reader, CSolver *solver_container, CGeometry *geometry);

    ~BFMInterpolator();


};