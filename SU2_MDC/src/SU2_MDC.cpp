/*!
 * \file SU2_MDC.cpp
 * \brief Main file of Mesh Deformation Code (SU2_MDC).
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

#include "../include/SU2_MDC.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	unsigned short iChunk;
	char buffer[50], out_file[200]; 
	
	/*--- Definition of the configuration class ---*/
	CConfig *config;
	if (argc == 2) config = new CConfig(argv[1], SU2_MDC);
	else {
		char grid_file[200];
			strcpy (grid_file, "default.cfg");
		config = new CConfig(grid_file, SU2_MDC);
	}
	
	/*--- Definition of the geometry class ---*/
	CPhysicalGeometry *geometry = new CPhysicalGeometry(config, config->GetMesh_FileName(), config->GetMesh_FileFormat());
	
	cout << endl <<"----------------------- Preprocessing computations ----------------------" << endl;
	
	/*--- Compute elements surrounding points, points surrounding points, and elements surronding elements ---*/
	cout << "Setting local point and element connectivity." <<endl; 
	geometry->SetEsuP(); geometry->SetPsuP(); geometry->SetEsuE(); 
	
	/*--- Check the orientation before computing geometrical quantities ---*/
	cout << "Check numerical grid orientation." <<endl; 
	geometry->SetBoundVolume(); geometry->Check_Orientation(config); 
	
	/*--- Create the edge structure ---*/
	cout << "Identify edges and vertices." <<endl; 
	geometry->SetEdges(); geometry->SetVertex(config); geometry->SetCG();
	
	/*--- Create the control volume structures ---*/
	cout << "Set control volume structure." << endl; 
	geometry->SetControlVolume(config, ALLOCATE); geometry->SetBoundControlVolume(config, ALLOCATE);

	cout << endl <<"-------------------- Start numerical grid deformation -------------------" << endl;	

	/*--- Definition and initialization of the surface deformation class ---*/
	CSurfaceMovement *surface_mov = new CSurfaceMovement();
	surface_mov->CopyBoundary(geometry, config);

	/*--- Definition of the FFD deformation class, only 1 chunk is defined ---*/
	unsigned short nChunk = MAX_NUMBER_CHUNCK;
	CFreeFormChunk** chunk; 
	chunk = new CFreeFormChunk*[nChunk];

	/*--- Output original computational grids ---*/
	if (config->GetVisualize_Deformation()) {
		if (geometry->GetnDim() == 2) { strcpy (out_file, "original_volumetric_grid.vtk"); geometry->SetParaView(out_file); }
		if (geometry->GetnDim() == 3) { 
			strcpy (out_file, "original_volumetric_grid.vtk"); geometry->SetParaView(out_file);
			strcpy (out_file, "original_surface_grid.vtk"); geometry->SetBoundParaView(config, out_file);
		}
	}

	/*--- Bump deformation for 2D problems ---*/
	if (geometry->GetnDim() == 2) {
		cout << "Perform 2D deformation of the surface." << endl;
		
		/*--- Apply the design variables to the control point position ---*/
		for (unsigned short iDV = 0; iDV < config->GetnDV(); iDV++) {
			switch ( config->GetDesign_Variable(iDV) ) {
				case HICKS_HENNE :
					surface_mov->SetHicksHenne(geometry, config, iDV); break;
				case DISPLACEMENT :
					surface_mov->SetDisplacement(geometry, config, iDV); break;
				case ROTATION :
					surface_mov->SetRotation(geometry, config, iDV); break;
				case NACA_4DIGITS :
					surface_mov->SetNACA_4Digits(geometry, config); break;
				case PARABOLIC :
					surface_mov->SetParabolic(geometry, config); break;
			}
		}
	}

	/*--- Free Form deformation for 3D problems ---*/
	if (geometry->GetnDim() == 3) {
		if (config->GetDesign_Variable(0) != NO_DEFORMATION)
			cout << "Performing the deformation of the surface grid." << endl;
		else
			cout << "No deformation of the surface grid." << endl;

		/*--- Read the FFD information fron the grid file ---*/
		surface_mov->ReadFFDInfo(config, geometry, chunk, config->GetMesh_FileName());

		/*--- If the chunk was not defined in the input file ---*/
		if (!surface_mov->GetChunkDefinition()) {
			
			/*--- Create a unitary chunk as baseline for other chunks shapes ---*/
			CFreeFormChunk chunk_unitary(1,1,1);
			chunk_unitary.SetUnitCornerPoints();
			
			/*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/
			chunk_unitary.SetControlPoints_Parallelepiped();
			
			
			for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
				/*--- Compute the support control points for the final FFD using the unitary box ---*/
				chunk_unitary.SetSupportCP(chunk[iChunk]);
				
				/*--- Cambiar puntos de control de la caja unitaria y recalcular las coordenadas 
				 cartesianas de los puntos de apoyo ---*/
				chunk_unitary.SetSupportCPChange(chunk[iChunk]);
				
				/*--- Set the nodes that belong to a chunk, this is a first estimation based in boxes ---*/
				chunk[iChunk]->SetChunkDomain(geometry, config, iChunk);
			}
			
			/*--- Compute the parametric coordinates, it also find the points in 
			 the chunk using the parametrics cooridnates ---*/
			surface_mov->SetParametricCoord(geometry, config, chunk);
		}
		
		/*--- Plot original chunks ---*/
		for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
			if(config->GetOutput_FileFormat() == PARAVIEW) {
				if (surface_mov->GetnChunk() == 1) sprintf (buffer, "original_chunk.vtk"); 
				else sprintf (buffer, "original_%dchunk.vtk", int(iChunk));
				chunk[iChunk]->SetParaView(buffer);
			}
			if(config->GetOutput_FileFormat() == TECPLOT) {
				if (surface_mov->GetnChunk() == 1) sprintf (buffer, "original_chunk.plt"); 
				else sprintf (buffer, "original_%dchunk.plt", int(iChunk));
				chunk[iChunk]->SetTecplot(buffer);
			}
		}
		
		/*--- Apply the design variables to the control point position ---*/
		for (unsigned short iDV = 0; iDV < config->GetnDV(); iDV++) {
			switch ( config->GetDesign_Variable(iDV) ) {
				case FFD_CONTROL_POINT : 
					surface_mov->SetFFDCPChange(geometry, config, chunk, iDV); break;
				case FFD_DIHEDRAL_ANGLE :
					surface_mov->SetFFDDihedralAngle(geometry, config, chunk, iDV); break;
				case FFD_TWIST_ANGLE :
					surface_mov->SetFFDTwistAngle(geometry, config, chunk, iDV); break;
				case FFD_ROTATION :
					surface_mov->SetFFDRotation(geometry, config, chunk, iDV); break;
				case FFD_CAMBER :
					surface_mov->SetFFDCamber(geometry, config, chunk, iDV); break;
				case FFD_THICKNESS :
					surface_mov->SetFFDThickness(geometry, config, chunk, iDV); break;
				case FFD_VOLUME :
					surface_mov->SetFFDVolume(geometry, config, chunk, iDV); break;
			}
		}
		
		/*--- Recompute cartesian coordinates using the new control points position ---*/
		surface_mov->SetCartesianCoord(geometry, config, chunk);
		

		/*--- Plot deformed chunks ---*/
		if (config->GetDesign_Variable(0) != NO_DEFORMATION)
			for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++) {
				if(config->GetOutput_FileFormat() == PARAVIEW) {
					if (surface_mov->GetnChunk() == 1) sprintf (buffer, "deformed_chunk.vtk"); 
					else sprintf (buffer, "deformed_%dchunk.vtk", int(iChunk));
					chunk[iChunk]->SetParaView(buffer);
				}
				if(config->GetOutput_FileFormat() == TECPLOT) {
					if (surface_mov->GetnChunk() == 1) sprintf (buffer, "deformed_chunk.plt"); 
					else sprintf (buffer, "deformed_%dchunk.plt", int(iChunk));
					chunk[iChunk]->SetTecplot(buffer);
				}
			}
		
		/*--- Mark the volumetric zone that is going to be deformed (if a selective deformation is performed) ---*/
		for (iChunk = 0; iChunk < surface_mov->GetnChunk(); iChunk++)
			chunk[iChunk]->SetDeformationZone(geometry, config, iChunk);

	}

	/*--- Definition of the Class for grid movement ---*/
	CVolumetricMovement *grid_movement; grid_movement = new CVolumetricMovement(geometry);

	/*--- Grid deformation ---*/
	if (config->GetDesign_Variable(0) != NO_DEFORMATION) {
		cout << "Performing the deformation of the volumetric grid." << endl;
		if (config->GetKind_GridDef_Method() == ALGEBRAIC) grid_movement->AlgebraicMethod(geometry, config);
		if (config->GetKind_GridDef_Method() == SPRING) grid_movement->SpringMethod(geometry, config);
		if (config->GetKind_GridDef_Method() == TORSIONAL_SPRING) grid_movement->TorsionalSpringMethod(geometry, config);
	}
	else cout << "No deformation deformation of the volumetric grid." << endl;

	
	/*--- Output grid ---*/
	cout << "End and write output files." << endl;
	
	if (config->GetVisualize_Deformation()) {
		if (config->GetDesign_Variable(0) != NO_DEFORMATION) {
			strcpy (out_file, "deformed_volumetric_grid.vtk"); geometry->SetParaView(out_file);
			strcpy (out_file, "deformed_surface_grid.vtk"); geometry->SetBoundParaView(config, out_file);
		}
	}
	
	geometry->SetMeshFile(config,config->GetMesh_Out_FileName());
	if (geometry->GetnDim() == 3) surface_mov->WriteFFDInfo(geometry, config, chunk, config->GetMesh_Out_FileName());
	
	return 1;
	
}
