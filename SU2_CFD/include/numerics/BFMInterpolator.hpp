/*!
 * \file ReadBFMInput.hpp
 * \brief Declarations and inlines of the classes used to interpolate
 *        the blade geometric parameters onto the mesh nodes.
 * \author E.C.Bunschoten
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
    su2double *Rotation_Axis;   // BFM rotation axis (Cartesian)
    su2double *axial_direction; // Axial direction unit vector (Cartesian)
    su2double *right_up_axis;

    su2double BFM_radius;
    /*! 
    * \brief 2D interpolation of blade geometric parameters onto mesh nodes.
    * \param [in] axis - Axial coordinate of mesh node.
    * \param [in] radius - Radial coordinate of mesh node.
    * \param [in] iPoint - Mesh node index.
    * \param [in] reader - Pointer to input file reader class.
    * \param [in] solver_container - Pointer to solver container class. 
    */
    void Interp2D(su2double axis, su2double radius, unsigned long iPoint, const ReadBFMInput * reader, CSolver * solver_container);

    void Interp3D(su2double axis, su2double radius, su2double theta, unsigned long iPoint, const ReadBFMInput * reader, CSolver * solver_container);

    bool Interp_ax_rad(su2double axis, su2double radius, unsigned long iPoint, const ReadBFMInput * reader, unsigned long iRow, unsigned long iTang, vector<su2double> *interp_solution);
    /*!
    * \brief Ray-cast inclusion algorithm which checks for query point inclusion whithin cell.
    * \param [in] axis - Axial coordinate of query point.
    * \param [in] radius - Radial coordinate of query point.
    * \param [in] ax_cell - Axial boundary coordinates of query cell.
    * \param [in] rad_cell - Radial boundary coordinates of query cell.
    * \returns Inclusion of query point into query cell.
    */
    bool RC_Inclusion(su2double axis, su2double radius, vector<su2double>ax_cell, vector<su2double>rad_cell);

    /*!
    * \brief Distance-weighted interpolation algorithm.
    * \param [in] axis - Axial coordinate of query point.
    * \param [in] radius - Radial coordinate of query point.
    * \param [in] ax_cell - Axial boundary coordinates of query cell.
    * \param [in] rad_cell - Radial boundary coordinates of query cell.
    * \param [in] val_cell - Data points at cell boundary.
    * \returns Interpolated data value at query point.
    */
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
    BFMInterpolator(const ReadBFMInput *reader, const CSolver * solver_container, const CGeometry * geometry, const CConfig * config);

    /*! 
    * \brief Interpolates the blade shape parameters onto the mesh nodes.
    * \param[in] reader - pointer to BFM input file reader class.
    * \param[in] solver_container - pointer to BFM solver class.
    * \param[in] geometry - pointer to geometry class.
    */
    void Interpolate(const ReadBFMInput * reader, CSolver * solver_container, const CGeometry * geometry);

    ~BFMInterpolator() {delete [] axial_direction; delete [] Rotation_Axis;};


};