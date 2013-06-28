/*!
 * \file output_structure.cpp
 * \brief Main subroutines for output solver information.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.2
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2013 Aerospace Design Laboratory
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

#include "../include/output_structure.hpp"
string AssembleVariableNames(bool GridMovement, bool Incompressible, unsigned short KindSolver, unsigned short nVar_Consv, unsigned short dims, unsigned short *NVar);

void COutput::WriteTecplotBinary(CConfig *config,CGeometry *geometry, unsigned short iZone) { 

#ifndef NO_TECIO	

	double   t;
    INTEGER4 i, N, iVar, err, Debug, NPts, NElm, IsDouble, IMax, JMax, KMax;
    INTEGER4 ICellMax, JCellMax, KCellMax, ZoneType, StrandID, ParentZn, FileType;
    INTEGER4 *ShareFromZone, *connectivity, IsBlock, NumFaceConnections, FaceNeighborMode, ShareConnectivityFromZone;
	string buffer, variables;
	stringstream file;
	bool first_zone = true, unsteady = config->GetUnsteady_Simulation(), GridMovement = config->GetGrid_Movement();
	unsigned long iExtIter = config->GetExtIter();
	unsigned short NVar, dims = geometry->GetnDim(); 
	enum     FileType { FULL = 0, GRID = 1, SOLUTION = 2 };
	enum	 ZoneType { ORDERED=0, FELINESEG=1, FETRIANGLE=2, FEQUADRILATERAL=3, FETETRAHEDRON=4, FEBRICK=5, FEPOLYGON=6, FEPOLYHEDRON=7 };

	/*--- Consistent data for Tecplot zones ---*/
	Debug						= 0;
	IsDouble					= 1;
	NPts						= nGlobal_Poin;
	t							= iExtIter*config->GetDelta_UnstTimeND();
	KMax						= 0;
	ICellMax					= 0;
	JCellMax					= 0;
	KCellMax					= 0;
	StrandID					= 0;
	ParentZn					= 0;
	IsBlock						= 1;
	NumFaceConnections			= 0;
	FaceNeighborMode			= 0;
	ShareConnectivityFromZone	= 0; 

	/*--- Write Tecplot solution file ---*/
	if (!wrote_Tecplot_base) {

		file.str(string());
		buffer = config->GetFlow_FileName();
		
#ifndef NO_MPI
	/*--- Remove the domain number from the filename ---*/
	int nProcessor = MPI::COMM_WORLD.Get_size();
	if (nProcessor > 1) buffer.erase(buffer.end()-2, buffer.end());
#endif
		
		file << buffer << ".mesh.plt";
		FileType = GRID;

		if (dims == 2) variables = "x y";
		else if (dims == 3) variables = "x y z";
		else cout << "Error: wrong number of dimentsions: " << dims << endl;

		/*--- Open Tecplot file ---*/
		cout << "Opening Tecplot mesh file..." << endl;
		err = TECINI112((char *)config->GetFlow_FileName().c_str(),
                  (char *)variables.c_str(),
                  (char *)file.str().c_str(),
                  (char *)".",
                  &FileType,
                  &Debug,
                  &IsDouble);
		if (err) cout << "Error in opening Tecplot file" << endl;

		first_zone = true;
		ShareFromZone = new INTEGER4[dims]; 
		for (i = 0; i < dims; i++) ShareFromZone[i] = 0;

		if (nGlobal_Tria > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FETRIANGLE; NElm = nGlobal_Tria; N = NElm*N_POINTS_TRIANGLE;

			err = TECZNE112((char*)"Triangle Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  ShareFromZone,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}

			connectivity = new int[N];
			for (i = 0; i < N; i++) connectivity[i] = Conn_Tria[i] + 1;
			err = TECNOD112(connectivity);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
			delete [] connectivity; 

		}
		if (nGlobal_Quad > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FEQUADRILATERAL; NElm = nGlobal_Quad; N = NElm*N_POINTS_QUADRILATERAL;

			err = TECZNE112((char*)"Quadrilateral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  ShareFromZone,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}

			connectivity = new int[N];
			for (i = 0; i < N; i++) connectivity[i] = Conn_Quad[i] + 1;
			err = TECNOD112(connectivity);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
			delete [] connectivity; 

		}
		if (nGlobal_Tetr > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FETETRAHEDRON; NElm = nGlobal_Tetr; N = NElm*N_POINTS_TETRAHEDRON; 

			err = TECZNE112((char*)"Tetrahedral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}

			connectivity = new int[N];
			for (i = 0; i < N; i++) connectivity[i] = Conn_Tetr[i] + 1;
			err = TECNOD112(connectivity);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
			delete [] connectivity;
		}
		if (nGlobal_Hexa > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FEBRICK; NElm = nGlobal_Hexa; N = NElm*N_POINTS_HEXAHEDRON; 

			err = TECZNE112((char*)"Hexahedral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}

			connectivity = new int[N];
			for (i = 0; i < N; i++) connectivity[i] = Conn_Hexa[i] + 1;
			err = TECNOD112(connectivity);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
			delete [] connectivity;
		}
		if (nGlobal_Pyra > 0) {
			cout << "Pyramid element type not yet supported; no zone written." << endl;
		}
		if (nGlobal_Wedg > 0) {
			cout << "Wedge element type not yet supported; no zone written." << endl;
		}
		if (nGlobal_Line > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FELINESEG; NElm = nGlobal_Line; N = NElm*N_POINTS_LINE; 

			err = TECZNE112((char*)"Line Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[0] = 1;
				err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[1] = 1;
				if (geometry->GetnDim() == 3) {
					err = TECDAT112(&NPts, Coords[2], &IsDouble);
					ShareFromZone[2] = 1;
				}
				if (err) cout << "Error writing coordinates to Tecplot file" << endl;
				first_zone = false; 
			}

			connectivity = new int[N];
			for (i = 0; i < N; i++) connectivity[i] = Conn_Line[i] + 1;
			err = TECNOD112(connectivity);
			if (err) cout << "Error writing connectivity to Tecplot file" << endl;
			delete [] connectivity;
		}

		delete [] ShareFromZone;
		wrote_Tecplot_base = true; 

		err = TECEND112();
		if (err) cout << "Error in closing Tecplot file" << endl;

	}

	file.str(string());
	buffer = config->GetFlow_FileName();
		
#ifndef NO_MPI
	/*--- Remove the domain number from the filename ---*/
	int nProcessor = MPI::COMM_WORLD.Get_size();
	if (nProcessor > 1) buffer.erase(buffer.end()-2, buffer.end());
#endif
	
	file << buffer;

	if (unsteady) {
		if ((iExtIter >= 0) && (iExtIter < 10))			file << "_0000" << iExtIter;
		if ((iExtIter >= 10) && (iExtIter < 100))		file << "_000" << iExtIter;
		if ((iExtIter >= 100) && (iExtIter < 1000))		file << "_00" << iExtIter;
		if ((iExtIter >= 1000) && (iExtIter < 10000))	file << "_0" << iExtIter;
		if (iExtIter >= 10000)							file << iExtIter;
	}
	file << ".sol.plt";
	FileType = SOLUTION;
	variables = AssembleVariableNames(config->GetGrid_Movement(), config->GetIncompressible(), config->GetKind_Solver(), nVar_Consv, dims, &NVar);

	/*--- Open Tecplot file ---*/
	cout << "Opening Tecplot results file..." << endl;
	err = TECINI112((char *)config->GetFlow_FileName().c_str(),
					(char *)variables.c_str(),
					(char *)file.str().c_str(),
					(char *)".",
					&FileType,
					&Debug,
					&IsDouble);
	if (err) cout << "Error in opening Tecplot file" << endl;

	first_zone = true;
	ShareFromZone = new INTEGER4[NVar]; 
	for (i = 0; i < NVar; i++) ShareFromZone[i] = 0;

	if (nGlobal_Tria > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FETRIANGLE; NElm = nGlobal_Tria; N = NElm*N_POINTS_TRIANGLE;

			err = TECZNE112((char*)"Triangle Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  ShareFromZone,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				i = 0;
				if (GridMovement) {
					err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					if (dims == 3) {
						err = TECDAT112(&NPts, Coords[2], &IsDouble);
						if (err) cout << "Error writing coordinates to Tecplot file" << endl;
						ShareFromZone[i++] = 1;
					}
				}
				for (iVar = 0; iVar < nVar_Total; iVar++) {
					err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing data to Tecplot file" << endl;
				}

				first_zone = false; 
			}

		}
	if (nGlobal_Quad > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FEQUADRILATERAL; NElm = nGlobal_Quad; N = NElm*N_POINTS_QUADRILATERAL;

			err = TECZNE112((char*)"Quadrilateral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  ShareFromZone,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				i = 0;
				if (GridMovement) {
					err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					if (dims == 3) {
						err = TECDAT112(&NPts, Coords[2], &IsDouble);
						if (err) cout << "Error writing coordinates to Tecplot file" << endl;
						ShareFromZone[i++] = 1;
					}
				}
				for (iVar = 0; iVar < nVar_Total; iVar++) {
					err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing data to Tecplot file" << endl;
				}

				first_zone = false; 
			}

		}
	if (nGlobal_Tetr > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FETETRAHEDRON; NElm = nGlobal_Tetr; N = NElm*N_POINTS_TETRAHEDRON; 

			err = TECZNE112((char*)"Tetrahedral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				i = 0;
				if (GridMovement) {
					err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					if (dims == 3) {
						err = TECDAT112(&NPts, Coords[2], &IsDouble);
						if (err) cout << "Error writing coordinates to Tecplot file" << endl;
						ShareFromZone[i++] = 1;
					}
				}
				for (iVar = 0; iVar < nVar_Total; iVar++) {
					err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing data to Tecplot file" << endl;
				}

				first_zone = false; 
			}

		}
	if (nGlobal_Hexa > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FEBRICK; NElm = nGlobal_Hexa; N = NElm*N_POINTS_HEXAHEDRON; 

			err = TECZNE112((char*)"Hexahedral Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				i = 0;
				if (GridMovement) {
					err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					if (dims == 3) {
						err = TECDAT112(&NPts, Coords[2], &IsDouble);
						if (err) cout << "Error writing coordinates to Tecplot file" << endl;
						ShareFromZone[i++] = 1;
					}
				}
				for (iVar = 0; iVar < nVar_Total; iVar++) {
					err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing data to Tecplot file" << endl;
				}

				first_zone = false; 
			}

		}
	if (nGlobal_Pyra > 0) {
			cout << "Pyramid element type not yet supported; no zone written." << endl;
		}
	if (nGlobal_Wedg > 0) {
			cout << "Wedge element type not yet supported; no zone written." << endl;
		}
	if (nGlobal_Line > 0) {

			/*--- Write the zone header information ---*/
			ZoneType = FELINESEG; NElm = nGlobal_Line; N = NElm*N_POINTS_LINE; 

			err = TECZNE112((char*)"Line Elements",
                  &ZoneType,
                  &NPts,
                  &NElm,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &t,
                  &StrandID,
                  &ParentZn,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  0,         /* TotalNumFaceNodes */
                  0,         /* NumConnectedBoundaryFaces */
                  0,         /* TotalNumBoundaryConnections */
                  NULL,      /* PassiveVarList */
                  NULL,      /* ValueLocation */
                  NULL,      /* ShareVarFromZone */
                  &ShareConnectivityFromZone);
			if (err) cout << "Error writing Tecplot zone data" << endl;

			/*--- write node coordinates and data if not done already---*/
			if (first_zone) {
 
				i = 0;
				if (GridMovement) {
					err = TECDAT112(&NPts, Coords[0], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					err = TECDAT112(&NPts, Coords[1], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing coordinates to Tecplot file" << endl;
					if (dims == 3) {
						err = TECDAT112(&NPts, Coords[2], &IsDouble);
						if (err) cout << "Error writing coordinates to Tecplot file" << endl;
						ShareFromZone[i++] = 1;
					}
				}
				for (iVar = 0; iVar < nVar_Total; iVar++) {
					err = TECDAT112(&NPts, Data[iVar], &IsDouble); ShareFromZone[i++] = 1;
					if (err) cout << "Error writing data to Tecplot file" << endl;
				}

				first_zone = false; 
			}

		}

	delete [] ShareFromZone;

	err = TECEND112();
	if (err) cout << "Error in closing Tecplot file" << endl;

#else // Not built with Tecplot binary support

	cout << "Tecplot binary file requested but SU^2 was built without Tecio support. No file written" << "\n"; 

#endif

}

string AssembleVariableNames(bool GridMovement, bool Incompressible, unsigned short KindSolver, unsigned short nVar_Consv, unsigned short dims, unsigned short *NVar) {

	stringstream variables; 
	unsigned short iVar; 
	*NVar = 0;

	variables.str(string());
	if (GridMovement) {
		variables << "x y "; *NVar += 2; 
		if (dims == 3) {
			variables << "z "; *NVar += 1;
		}
	}
	for (iVar = 0; iVar < nVar_Consv; iVar++){ 
		variables << "Conservative_Variable_" << iVar+1 << " "; 
		*NVar += 1;
	}
	for (iVar = nVar_Consv; iVar < 2*nVar_Consv; iVar++) {
		variables << "Conservative_Residual_" << iVar-nVar_Consv+1 << " "; 
		*NVar += 1;
	}
	if (GridMovement) {
			variables << "Grid_Velocity_X " << "Grid_Velocity_Y "; *NVar += 2;
			if (dims == 3) {
				variables << "Grid_Velocity_Z "; *NVar += 1;
			}
	}

	variables << "Pressure " << "Mach "; *NVar += 2;
	if (!Incompressible) {
		switch (KindSolver) {

			/*--- Include temperature and laminar viscosity, if applicable ---*/
			case NAVIER_STOKES: variables << "Temperature " << "Viscosity "; *NVar += 2; break; 

			/*--- Include eddy viscosity, if applicable ---*/
			case RANS: variables << "Temperature " << "Viscosity " << "Eddy_Viscosity "; *NVar += 3; break; 
		}
	}

	return variables.str(); 

}