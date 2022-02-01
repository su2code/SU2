/*!
 * \file ReadBFMInput.cpp
 * \brief 
 * \author E.C. Bunschoten
 * \version 7.2.1 "Blackbird"
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

#include "../../include/numerics/BFMInterpolator.hpp"
#include "../../include/solvers/CSolver.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

BFMInterpolator::BFMInterpolator(const ReadBFMInput * reader, const CSolver * solver_container, const CGeometry * geometry, const CConfig * config)
{
    unsigned long nPoint = geometry->GetnPoint();
    unsigned short nRows = reader->GetNBladeRows();
    Rotation_Axis = new su2double[3];
    axial_direction = new su2double[3];
    right_up_axis = new su2double[3];
    right_up_axis[0] = 0;
    right_up_axis[1] = 0;
    right_up_axis[2] = 1;

    if(geometry->GetnDim() == 2){
        BFM_radius = config->GetBFM_Radius();
    }
    for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
        Rotation_Axis[iDim] = config->GetBFM_Axis(iDim);
    }
    for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
        axial_direction[iDim] = Rotation_Axis[iDim];
    }

}

void BFMInterpolator::Interpolate(const ReadBFMInput * reader, CSolver * solver_container, const CGeometry *geometry){
    /*
    This function interpolates the BFM blade parameters to the mesh nodes. The BFM parameters listed in the reader
    are defined on a different grid than the SU2 mesh. Therefore, interpolation is required. Currently, ray-casting 
    is used to check for cell inclusion and distance-weighted averaging is used to interpolate the BFM parameters 
    from the reader mesh nodes. This works well for blade rows with repeating blades (all blades in the crown having
    the same geometry), but this might have to be changed to a more efficient lookup method in the future.
    */
    su2double *Coord_Cart, *Coord_Cyl, *Coord_rad;  // Cartesian, cylindrical cell coordinates
    su2double ax, rad, tang;                 // Axial, radial, and tangential coordinates
    su2double radial[geometry->GetnDim()];   // Radial projection vector      

    unsigned short barwidth = 70;              // Progress bar width
    su2double progress{0};          // Interpolation progress

    auto rank = SU2_MPI::GetRank();
    auto size = SU2_MPI::GetSize();

    // Constructing the perpendicular axis from the axial direction and right-up axis.
    su2double *perp_axis = new su2double[3];
    GeometryToolbox::CrossProduct(axial_direction, right_up_axis, perp_axis);

    // Looping over the solver mesh nodes and interpolating BFM parameters
    for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint){

        // Updating the interpolation progress bar
        if(rank == MASTER_NODE){
            progress = su2double(iPoint) / geometry->GetnPoint();
            unsigned short pos = floor(barwidth * progress);
            cout << "[";
            for(unsigned short iBar=0; iBar<barwidth; ++iBar){
                if(iBar < pos) cout << "=";
                else cout << " ";
            }
            cout << "] "<< floor(100*progress) << " %\r";
            cout.flush();
        }

        // Getting solver mesh Cartesian coordinates
        Coord_Cart = geometry->nodes->GetCoord(iPoint);

        // Computing axial projection of solver mesh coordinates
        ax = GeometryToolbox::DotProduct(geometry->GetnDim(), Coord_Cart, axial_direction);
        
        if(geometry->GetnDim() == 2){
            rad = BFM_radius;
        }else{
            // Computing radial projection of solver mesh coordinates
            for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
                radial[iDim] = Coord_Cart[iDim] - ax*axial_direction[iDim];
            }
            
            // Computing radius of solver mesh coordinates
            rad = GeometryToolbox::Norm(geometry->GetnDim(), radial);
        }
        
        // Computing the coordinates along the right-up axis and the perpendicular axis.
        su2double right_up_coord = GeometryToolbox::DotProduct(3, right_up_axis, radial);
        su2double perp_coord = GeometryToolbox::DotProduct(3, perp_axis, radial);

        // Computing the tangential angle according to the left-hand rule.
        tang = atan2(perp_coord, right_up_coord) + PI_NUMBER;

        // Interpolating the BFM parameters according to the axial, radial, and tangential coordinates of the current mesh node.
        Interp3D(ax, rad, tang, iPoint, reader, solver_container);
        
    }
    if(rank == MASTER_NODE) cout << endl;
}

void BFMInterpolator::Interp3D(su2double axis, su2double radius, su2double theta, unsigned long iPoint, const ReadBFMInput * reader, CSolver * solver_container){
    /*
    Perform 3D interpolation according to the meridional channel shape presented in the BFM input file. 
    */

    vector<su2double> interp_solution_lower, interp_solution_upper; // interpolated blade parameters on enclosing tangential sections
    interp_solution_lower.resize(N_BFM_PARAMS);
    interp_solution_upper.resize(N_BFM_PARAMS);

    su2double theta_sec, theta_next_sec, theta_upper, theta_lower, var_interp;
    unsigned long iTang_lower, iTang_upper;

    bool found, single_section;
    // Looping over blade rows
    for(unsigned long iRow = 0; iRow < reader->GetNBladeRows(); iRow++){

        // If the number of tangential sections for the current blade row is equal to 1, only axial-radial interpolation is performed
        unsigned long nTang{reader->GetNTangentialPoints(iRow)};
        single_section = (nTang == 1);

        // Going over the tangential sections to find the enclosing sections.
        for(unsigned long iTang=0; iTang<nTang; iTang++){
            theta_sec = iTang * 2*PI_NUMBER/nTang;
            theta_next_sec = (iTang+1) * 2*PI_NUMBER/nTang;
            if((theta >= theta_sec) && (theta <= theta_next_sec)){
                theta_lower = theta_sec;
                theta_upper = theta_next_sec;
                iTang_lower = iTang;
                iTang_upper = (iTang + 1) % nTang;
            }
        }

        // Perform axial-radial interpolation on the enclosing tangential section(s)
        found = Interp_ax_rad(axis, radius, iPoint, reader, iRow, iTang_lower, &interp_solution_lower);
        if(!single_section)
            found = Interp_ax_rad(axis, radius, iPoint, reader, iRow, iTang_upper, &interp_solution_upper);
        
    }

    // Looping over BFM variables to store interpolated values in the auxilary variable.
    for(size_t iVar=0; iVar<N_BFM_PARAMS; iVar++){
        // If the mesh node lies within the meridional channel contour of the blade, the interpolated values are assigned to the BFM parameters.
        if(found){
            // In case of multiple tangential sections, linear interpolation between the enclosing sections is done.
            if(!single_section){
                var_interp = interp_solution_lower[iVar] + (interp_solution_upper[iVar] - interp_solution_lower[iVar])*(theta - theta_lower)/(theta_upper - theta_lower);
            }else{
                var_interp = interp_solution_lower[iVar];
            }
            solver_container->GetNodes()->SetAuxVar(iPoint, iVar, var_interp);
        // In case the node lies outside the bladed region, default values are assigned.
        }else{
            if(iVar == I_BLOCKAGE_FACTOR){
                solver_container->GetNodes()->SetAuxVar(iPoint, iVar, 1);
            }else
                solver_container->GetNodes()->SetAuxVar(iPoint, iVar, 0);
        }
        
    }
}

bool BFMInterpolator::Interp_ax_rad(su2double axis, su2double radius, unsigned long iPoint, const ReadBFMInput * reader, unsigned long iRow, unsigned long iTang, vector<su2double>* interp_solution){
    /*
    Axial-radial interpolation on the bladed region of blade row iRow
    */
    unsigned long iRad{0};
    unsigned long iAx{0};
    bool found{false};
    vector<su2double> ax_cell{0, 0, 0, 0};
    vector<su2double> rad_cell{0, 0, 0, 0};
    vector<su2double> val_cell{0, 0, 0, 0};

    // Looping over the number of spanwise sections
    iRad = 0;
    while((iRad < reader->GetNRadialPoints(iRow)-1) && !found){
        // Looping over the number of axial sections
        iAx = 0;
        while((iAx < reader->GetNAxialPoints(iRow)-1) && !found){

            // Defining a rectangle between concurrent axial and spanwise points
            ax_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_AXIAL_COORDINATE);
            ax_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, I_AXIAL_COORDINATE);
            ax_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, I_AXIAL_COORDINATE);
            ax_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, I_AXIAL_COORDINATE);
            rad_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_RADIAL_COORDINATE);
            rad_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, I_RADIAL_COORDINATE);
            rad_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, I_RADIAL_COORDINATE);
            rad_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, I_RADIAL_COORDINATE);
            
            // Checking if the mesh node point lies within the rectangle
            found = RC_Inclusion(axis, radius, ax_cell, rad_cell);
            if(found){
                // In case the node lies within the rectangle, the BFM parameters are interpolated onto the node
                // through distance-weighted averaging.
                for(unsigned short iVar=0; iVar<N_BFM_PARAMS; ++iVar){
                    if((iVar != I_ROTATION_FACTOR) && (iVar != I_BODY_FORCE_FACTOR) && (iVar != I_BLADE_COUNT)){
                        val_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, iVar);
                        val_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, iVar);
                        val_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, iVar);
                        val_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, iVar);
                        
                        interp_solution->at(iVar) = DW_average(axis, radius, ax_cell, rad_cell, val_cell);
                        
                    }
                    interp_solution->at(I_ROTATION_FACTOR) = reader->GetRotationFactor(iRow);
                    interp_solution->at(I_BLADE_COUNT) = reader->GetBladeCount(iRow);
                    interp_solution->at(I_BODY_FORCE_FACTOR) = 1;
                }
            }
            ++iAx;
        }
        ++iRad;
        
    }
    return found;
}

bool BFMInterpolator::RC_Inclusion(su2double axis, su2double radius, vector<su2double> ax_cell, vector<su2double> rad_cell){
    unsigned short n_int {0};
	bool inside = false;
	unsigned short n_nodes = ax_cell.size();
	unsigned short i_next;
	su2double ax_min {0};
	su2double determinant, S, R;
	for(unsigned short i=0; i<n_nodes; i++){
        if(ax_cell.at(i) < ax_min) ax_min = ax_cell.at(i);
	}
	ax_min -= 1;
	for(unsigned short i=0; i<n_nodes; i++){
		i_next = (i+1) % n_nodes;
		determinant = (axis - ax_min)*(rad_cell.at(i_next) - rad_cell.at(i));
		if(determinant != 0){
			S = (1 / determinant) * (axis - ax_min)*(radius - rad_cell.at(i));
			R = (1 / determinant) * ((rad_cell.at(i) - rad_cell.at(i_next))*(ax_min - ax_cell.at(i)) + 
			(ax_cell.at(i_next) - ax_cell.at(i))*(radius - rad_cell.at(i)));
			if((S >= 0 && S <= 1) && (R >= 0 && R <= 1)){
				n_int ++;
			}
		}
	}
	if(n_int % 2 != 0){
		inside = true;
	}
	return inside;
}

su2double BFMInterpolator::DW_average(su2double axis, su2double radius, vector<su2double> ax_cell, vector<su2double> rad_cell, vector<su2double> val_cell){
    su2double enumerator{0};
    su2double denomintor{0};
    su2double distance{0};
    for(size_t i_node=0; i_node<ax_cell.size(); ++i_node){
        su2double cell_coords[2] = {axis, radius};
        su2double node_coords[2] = {ax_cell.at(i_node), rad_cell.at(i_node)};
        distance = GeometryToolbox::Distance(2, cell_coords, node_coords);
        if(distance == 0){
            return val_cell.at(i_node);
        }
        enumerator += (1/distance)*val_cell.at(i_node);
        denomintor += (1/distance);
    }
    return enumerator/denomintor;
}
