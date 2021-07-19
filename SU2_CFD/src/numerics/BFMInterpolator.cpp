/*!
 * \file ReadBFMInput.cpp
 * \brief 
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

#include "../../include/numerics/BFMInterpolator.hpp"
#include "../../include/solvers/CSolver.hpp"

BFMInterpolator::BFMInterpolator()
{

}
BFMInterpolator::BFMInterpolator(ReadBFMInput *reader, CSolver*solver_container, CGeometry *geometry, CConfig *config)
{
    unsigned long nPoint = geometry->GetnPoint();
    unsigned short nRows = reader->GetNBladeRows();
    Rotation_Axis.resize(geometry->GetnDim());
    axial_direction.resize(geometry->GetnDim());
    su2double rotation_rate{0};
    // for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
    //     Rotation_Axis.at(iDim) = config->GetRotation_Rate(iDim);
    //     rotation_rate += pow(config->GetRotation_Rate(iDim), 2);
    // }
    // rotation_rate = sqrt(rotation_rate);
    Rotation_Axis.at(0) = 0;
    Rotation_Axis.at(1) = 0;
    Rotation_Axis.at(2) = 1;
    rotation_rate = 1;
    for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
        axial_direction.at(iDim) = Rotation_Axis.at(iDim)/rotation_rate;
    }

}

void BFMInterpolator::Interpolate(ReadBFMInput *reader, CSolver *solver_container, CGeometry *geometry){

    su2double *Coord_Cart, *Coord_Cyl;
    su2double ax, rad, tang;
    vector<su2double> axial_vector{}, Coord_Cart_v{};
    vector<su2double> radial_vector{};

    axial_vector.resize(geometry->GetnDim());
    radial_vector.resize(geometry->GetnDim());
    Coord_Cart_v.resize(geometry->GetnDim());
    for(unsigned long iPoint=0; iPoint<geometry->GetnPoint(); ++iPoint){
        Coord_Cart = geometry->nodes->GetCoord(iPoint);
        ax = Vector_Dot_Product(Coord_Cart, axial_direction);
        
        for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
            axial_vector.at(iDim) = ax*axial_direction.at(iDim);
            Coord_Cart_v.at(iDim) = Coord_Cart[iDim];
            radial_vector.at(iDim) = Coord_Cart_v.at(iDim) - axial_vector.at(iDim);
        }
        //radial_vector = Vector_Cross_Product(axial_vector, Coord_Cart_v);
        
        rad = 0;
        if(radial_vector.size() == 1){
            rad = radial_vector.at(0);
            tang = 0;
        }else{
            for(unsigned short iDim=0; iDim<geometry->GetnDim(); ++iDim){
                rad += pow(radial_vector.at(iDim), 2);
            }
            rad = sqrt(rad);
            tang = 0;
        }
        
        Interp2D(ax, rad, iPoint, reader, solver_container);
        
        
    }
}

void BFMInterpolator::Interp2D(su2double axis, su2double radius, unsigned long iPoint, ReadBFMInput *reader, CSolver*solver_container){
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
                //cout << axis << " " << radius << " " << ax_cell.at(0) << " " << ax_cell.at(1) << " " << ax_cell.at(2) << " " << ax_cell.at(3) << " " << rad_cell.at(0) << " " << rad_cell.at(1) << " " << rad_cell.at(2) << " " << rad_cell.at(3)<< endl;
                    
                if(found){
                    //cout << axis << " " << radius << " " << ax_cell.at(0) << " " << ax_cell.at(1) << " " << ax_cell.at(2) << " " << ax_cell.at(3) << " " << rad_cell.at(0) << " " << rad_cell.at(1) << " " << rad_cell.at(2) << " " << rad_cell.at(3) << endl;
                    for(unsigned short iVar=0; iVar<N_BFM_PARAMS; ++iVar){
                        if((iVar != I_ROTATION_FACTOR) && (iVar != I_BODY_FORCE_FACTOR) && (iVar != I_BLADE_COUNT)){
                            val_cell.at(0) = reader->GetBFMParameter(iRow, iTang, iRad, iAx, iVar);
                            val_cell.at(1) = reader->GetBFMParameter(iRow, iTang, iRad, iAx+1, iVar);
                            val_cell.at(2) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx+1, iVar);
                            val_cell.at(3) = reader->GetBFMParameter(iRow, iTang, iRad+1, iAx, iVar);
                            
                            solver_container->GetNodes()->SetAuxVar(iPoint, iVar, DW_average(axis, radius, ax_cell, rad_cell, val_cell));
                        }
                    }
                    
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_ROTATION_FACTOR, reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_ROTATION_FACTOR));
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_BODY_FORCE_FACTOR, 1);
                    solver_container->GetNodes()->SetAuxVar(iPoint, I_BLADE_COUNT, reader->GetBFMParameter(iRow, iTang, iRad, iAx, I_BLADE_COUNT));
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
        distance = sqrt(pow((axis - ax_cell.at(i_node)), 2) + pow((radius - rad_cell.at(i_node)), 2));
        if(distance == 0){
            return val_cell.at(i_node);
        }
        enumerator += (1/distance)*val_cell.at(i_node);
        denomintor += (1/distance);
    }
    return enumerator/denomintor;
}
BFMInterpolator::~BFMInterpolator()
{
    
}