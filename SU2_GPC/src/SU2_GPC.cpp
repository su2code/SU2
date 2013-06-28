/*!
 * \file SU2_GPC.cpp
 * \brief Main file of the Gradient Projection Code (SU2_GPC).
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009). 
 * \version 1.0.
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_GPC.hpp"

using namespace std;

int main(int argc, char *argv[])
{	
	unsigned short iChunk, iMarker, iDim;
	unsigned long iVertex, iPoint;
	double delta_eps, Gradient, *Normal, dS, *VarCoord, Sensitivity, dalpha_dx, dalpha_dy, 
	dalpha_dz, dx_deps, dy_deps, dz_deps, dalpha_deps;
	char *cstr;
	ofstream Gradient_file;
	bool *UpdatePoint;

	/*--- Definition of the Class for the dssefinition of the problem ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_GPC);
	else {
		char grid_file[200];
		strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_GPC);
	}

	/*--- Definition of the Class for the boundary of the geometry ---*/
	CGeometry *boundary; boundary = new CGeometry;
	boundary = new CBoundaryGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
	boundary->SetVertex();
	boundary->SetBoundControlVolume(config, ALLOCATE);
	boundary->SetBoundSensitivity(config->GetSurfAdjCoeff_FileName());

	UpdatePoint = new bool[boundary->GetnPoint()];

	/*--- Definition of the Class for surface deformation ---*/
	CSurfaceMovement *surf_def = new CSurfaceMovement();

	cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;

	/*--- Bump deformation for 2D problems ---*/
	if (boundary->GetnDim() == 2) {
		cout << "Perform 2D deformation of the surface." << endl;
		if (config->GetDesign_Variable(0) == HICKS_HENNE)        surf_def->SetHicksHenne(boundary, config, 0);
		if (config->GetDesign_Variable(0) == DISPLACEMENT)  surf_def->SetDisplacement(boundary, config, 0);
		if (config->GetDesign_Variable(0) == ROTATION) surf_def->SetRotation(boundary, config, 0);
		if (config->GetDesign_Variable(0) == NACA_4DIGITS) surf_def->SetNACA_4Digits(boundary, config);
		if (config->GetDesign_Variable(0) == PARABOLIC) surf_def->SetParabolic(boundary, config);
	}
	
	/*--- Free Form deformation for 3D problems ---*/
	if (boundary->GetnDim() == 3) {

		cout << "Perform 3D deformation of the surface." << endl;

		/*--- Definition of the FFD deformation class, only 1 chunk is defined ---*/
		unsigned short nChunk = 1;
		CFreeFormChunk** chunk; 
		chunk = new CFreeFormChunk*[nChunk];
		
		/*--- Read the FFD information fron the grid file ---*/
		surf_def->ReadFFDInfo(config, boundary, chunk, config->GetMesh_FileName());
		
		/*--- If the chunk was not defined in the input file ---*/
		if (!surf_def->GetChunkDefinition()) {
			
			/*--- Create a unitary chunk as baseline for other chunks shapes ---*/
			CFreeFormChunk chunk_unitary(1,1,1);
			chunk_unitary.SetUnitCornerPoints();
			
			/*--- Compute the control points of the unitary box, in this case the order is 2 and the degree is 1 ---*/
			chunk_unitary.SetControlPoints_Parallelepiped();
			
			for (iChunk = 0; iChunk < surf_def->GetnChunk(); iChunk++) {
				/*--- Compute the support control points for the final FFD using the unitary box ---*/
				chunk_unitary.SetSupportCP(chunk[iChunk]);
				
				/*--- Cambiar puntos de control de la caja unitaria y recalcular las coordenadas 
				 cartesianas de los puntos de apoyo ---*/
				chunk_unitary.SetSupportCPChange(chunk[iChunk]);
				
				/*--- Set the nodes that belong to a chunk ---*/
				chunk[iChunk]->SetChunkDomain(boundary, config, iChunk);
			}
			
			/*--- Compute the parametric coordinates ---*/
			surf_def->SetParametricCoord(boundary, config, chunk);
		}
		
		/*--- Apply the control point change ---*/
		if (config->GetDesign_Variable(0) == FFD_CONTROL_POINT) surf_def->SetFFDCPChange(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_DIHEDRAL_ANGLE) surf_def->SetFFDDihedralAngle(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_TWIST_ANGLE) surf_def->SetFFDTwistAngle(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_ROTATION) surf_def->SetFFDRotation(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_CAMBER) surf_def->SetFFDCamber(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_THICKNESS) surf_def->SetFFDThickness(boundary, config, chunk, 0);
		if (config->GetDesign_Variable(0) == FFD_VOLUME) surf_def->SetFFDVolume(boundary, config, chunk, 0);
		
		/*--- Recompute cartesian coordinates using the new control points position ---*/
		surf_def->SetCartesianCoord(boundary, config, chunk);

	}
		
	/*--- Continuos adjoint gradient computation ---*/
	cout << "Evaluate functional gradient using the continuous adjoint strategy." << endl;
	
	for (iPoint = 0; iPoint < boundary->GetnPoint(); iPoint++)
		UpdatePoint[iPoint] = true;
	
	delta_eps = config->GetDV_Value_New(0);
	
	Gradient = 0.0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			iPoint = boundary->vertex[iMarker][iVertex]->GetNode();
			if ((config->GetMarker_All_Moving(iMarker) == YES) && UpdatePoint[iPoint]) {
				
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
				dS = 0.0; for (iDim = 0; iDim < boundary->GetnDim(); iDim++) dS += Normal[iDim]*Normal[iDim];
				dS = sqrt(dS);
				VarCoord = boundary->vertex[iMarker][iVertex]->GetVarCoord();
				Sensitivity = boundary->vertex[iMarker][iVertex]->GetAuxVar();
				
				dalpha_dx = Normal[0] / dS; 
				dalpha_dy = Normal[1] / dS;

				dx_deps = VarCoord[0] / delta_eps;
				dy_deps = VarCoord[1] / delta_eps;

				dalpha_deps = -dalpha_dx*dx_deps - dalpha_dy*dy_deps;
				
				if (boundary->GetnDim() == 3) {
					dalpha_dz = Normal[2] / dS;
					dz_deps = VarCoord[2] / delta_eps;
					dalpha_deps -= dalpha_dz*dz_deps;
				}
				
				Gradient += Sensitivity * dalpha_deps;
				UpdatePoint[iPoint] = false;
			}
		}

	/*--- Write the gradient in a external file ---*/
	cstr = new char [config->GetObjFunc_Grad_FileName().size()+1];
	strcpy (cstr, config->GetObjFunc_Grad_FileName().c_str());
	Gradient_file.open(cstr, ios::out);
	if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT) {
		cout << "Lift coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Lift coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) {
		cout << "Drag coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Drag coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT) {
		cout << "Sideforce coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Sideforce coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == PRESSURE_COEFFICIENT) {
		cout << "Pressure coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Pressure coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) {
		cout << "Moment x coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Moment x coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) {
		cout << "Moment y coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Moment y coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) {
		cout << "Moment z coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Moment z coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == EFFICIENCY) {
		cout << "Efficiency coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Efficiency coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {
		cout << "Equivalent area coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Equivalent area coeff. grad. using cont. adj." << endl;
	}
	if (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE) {
		cout << "Near-field pressure coefficient gradient: "<< Gradient << "." << endl;
		Gradient_file << "Near-field pressure coeff. grad. using cont. adj." << endl;
	}
	
	Gradient_file << Gradient << endl;
	Gradient_file.close();

	delete [] UpdatePoint;

	return 1;
}
