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
    Rotation_Axis = new su2double[geometry->GetnDim()];
    axial_direction = new su2double[geometry->GetnDim()];
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
    su2double *Coord_Cart, *Coord_Cyl;       // Cartesian, cylindrical cell coordinates
    su2double ax, rad, tang;                 // Axial, radial, and tangential coordinates
    su2double radial[geometry->GetnDim()];   // Radial projection vector      

    int barwidth = 70;              // Progress bar width
    su2double progress{0};          // Interpolation progress

    auto rank = SU2_MPI::GetRank();
    auto size = SU2_MPI::GetSize();

    // Looping over the solver mesh nodes and interpolating BFM parameters
    for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint){

        // Updating the interpolation progress bar
        if(rank == MASTER_NODE){
            progress = su2double(iPoint) / geometry->GetnPoint();
            int pos = barwidth * progress;
            cout << "[";
            for(int iBar=0; iBar<barwidth; ++iBar){
                if(iBar < pos) cout << "=";
                else cout << " ";
            }
            cout << "] "<< int(100*progress) << " %\r";
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

        // rad = 0;
        // for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
        //     rad += pow(radial[iDim], 2);
        // }
        // rad = sqrt(rad);
        tang = 0;
        
        Interp2D(ax, rad, iPoint, reader, solver_container);
        
        
    }
    if(rank == MASTER_NODE) cout << endl;
}

void BFMInterpolator::Interp2D(su2double axis, su2double radius, unsigned long iPoint, const ReadBFMInput * reader, CSolver * solver_container){
    unsigned short iRow{0};
    unsigned long iTang{0};
    unsigned long iRad{0};
    unsigned long iAx{0};
    bool found{false};
    vector<su2double> ax_cell{0, 0, 0, 0};
    vector<su2double> rad_cell{0, 0, 0, 0};
    vector<su2double> val_cell{0, 0, 0, 0};
    while((iRow < reader->GetNBladeRows()) && !found){
        iRad = 0;
        while((iRad < reader->GetNRadialPoints(iRow)-1) && !found){
            iAx = 0;
            while((iAx < reader->GetNAxialPoints(iRow)-1) && !found){
                ax_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_AXIAL_COORDINATE);
                ax_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, I_AXIAL_COORDINATE);
                ax_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, I_AXIAL_COORDINATE);
                ax_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, I_AXIAL_COORDINATE);
                rad_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_RADIAL_COORDINATE);
                rad_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, I_RADIAL_COORDINATE);
                rad_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, I_RADIAL_COORDINATE);
                rad_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, I_RADIAL_COORDINATE);
                
                found = RC_Inclusion(axis, radius, ax_cell, rad_cell);
                if(found){
                    for(unsigned short iVar=0; iVar<N_BFM_PARAMS; ++iVar){
                        if((iVar != I_ROTATION_FACTOR) && (iVar != I_BODY_FORCE_FACTOR) && (iVar != I_BLADE_COUNT)){
                            val_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, iVar);
                            val_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, iVar);
                            val_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, iVar);
                            val_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, iVar);
                            
                            solver_container->GetNodes()->SetAuxVar(iPoint, iVar, DW_average(axis, radius, ax_cell, rad_cell, val_cell));
                        }
                    }
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_BODY_FORCE_FACTOR, 1);
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_ROTATION_FACTOR, reader->GetRotationFactor(iRow));
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_BLADE_COUNT, reader->GetBladeCount(iRow));
                }
                ++iAx;
            }
            ++iRad;
            
        }
        ++iRow;
        
    }
    

    if(!found){
        for(unsigned short iVar=0; iVar<N_BFM_PARAMS; ++iVar){
            
            solver_container->GetNodes()->SetAuxVar(iPoint, iVar, 0);
            if(iVar == I_BLOCKAGE_FACTOR){
                solver_container->GetNodes()->SetAuxVar(iPoint, iVar, 1);
            }
        }
    }
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
        //distance = sqrt(pow((axis - ax_cell.at(i_node)), 2) + pow((radius - rad_cell.at(i_node)), 2));
        if(distance == 0){
            return val_cell.at(i_node);
        }
        enumerator += (1/distance)*val_cell.at(i_node);
        denomintor += (1/distance);
    }
    return enumerator/denomintor;
}
