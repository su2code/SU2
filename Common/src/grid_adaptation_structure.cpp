/*!
 * \file grid_adaptation_structure.cpp
 * \brief Main subroutines for grid adaptation
 * \author F. Palacios
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../include/grid_adaptation_structure.hpp"
#include <math.h>

CGridAdaptation::CGridAdaptation(CGeometry *geometry, CConfig *config) {

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();
  
	unsigned long iPoint;
	
	nDim = geometry->GetnDim();
	
	switch (config->GetKind_Solver()) {			
			
		default:
			nVar = geometry->GetnDim()+2;
			break;			
	}

	ConsVar_Sol = new su2double* [geometry->GetnPoint()];
	AdjVar_Sol = new su2double* [geometry->GetnPoint()];
	LinVar_Sol = new su2double* [geometry->GetnPoint()];
	ConsVar_Res = new su2double* [geometry->GetnPoint()];
	AdjVar_Res = new su2double* [geometry->GetnPoint()];	
	LinVar_Res = new su2double* [geometry->GetnPoint()];
	Gradient = new su2double* [geometry->GetnPoint()];
	Gradient_Flow = new su2double* [geometry->GetnPoint()];
	Gradient_Adj = new su2double* [geometry->GetnPoint()];

	Index = new su2double [geometry->GetnPoint()];

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		ConsVar_Sol[iPoint] = new su2double [nVar];
		AdjVar_Sol[iPoint] = new su2double [nVar];
		LinVar_Sol[iPoint] = new su2double [nVar];
		ConsVar_Res[iPoint] = new su2double [nVar];
		LinVar_Res[iPoint] = new su2double [nVar];
		AdjVar_Res[iPoint] = new su2double [nVar];		
		Gradient[iPoint] = new su2double [nDim];
		Gradient_Flow[iPoint] = new su2double [nDim];		
		Gradient_Adj[iPoint] = new su2double [nDim];				
	}

}

CGridAdaptation::~CGridAdaptation(void) {
	
	unsigned short iVar, iDim;
	
	for (iVar = 0; iVar < nVar; iVar++) {
		delete [] ConsVar_Adapt[iVar]; 
		delete [] ConsVar_Sol[iVar]; 
		delete [] ConsVar_Res[iVar];
		delete [] AdjVar_Adapt[iVar]; 
		delete [] AdjVar_Sol[iVar]; 
		delete [] AdjVar_Res[iVar];
		delete [] LinVar_Adapt[iVar]; 
		delete [] LinVar_Sol[iVar]; 
		delete [] LinVar_Res[iVar];
	}
	
	for (iDim = 0; iDim < nDim; iDim++) {
		delete [] Gradient[iDim];
		delete [] Gradient_Flow[iDim];
		delete [] Gradient_Adj[iDim];		
	}
	delete [] ConsVar_Adapt; 
	delete [] ConsVar_Sol;
	delete [] ConsVar_Res;
	delete [] AdjVar_Adapt; 
	delete [] AdjVar_Sol;
	delete [] AdjVar_Res;
	delete [] LinVar_Adapt; 
	delete [] LinVar_Sol; 
	delete [] LinVar_Res;
	delete [] Gradient;
	delete [] Gradient_Flow;	
	delete [] Gradient_Adj;	
	delete [] Index;	
}

void CGridAdaptation::GetFlowSolution(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	unsigned short iVar;
  su2double dummy;

	string text_line;
		
	string mesh_filename = config->GetSolution_FlowFileName();
	ifstream restart_file;

	char *cstr = new char [mesh_filename.size()+1];

	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
    SU2_MPI::Error("There is no flow restart file!!", CURRENT_FUNCTION);
  }
	
  /*--- Read the header of the file ---*/
  getline(restart_file, text_line);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		getline(restart_file, text_line);
		istringstream point_line(text_line);
		
		point_line >> index;
    
    if (nDim == 2) point_line >> dummy >> dummy;
    else point_line >> dummy >> dummy >> dummy;
      
		for (iVar = 0; iVar < nVar; iVar ++)
			point_line >> ConsVar_Sol[iPoint][iVar];
	}
	restart_file.close();
}

void CGridAdaptation::GetFlowResidual(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	unsigned short iVar;
	
//	su2double dummy[5];
	su2double dummy;
	string text_line;
	
	string mesh_filename = config->GetSolution_FlowFileName();
	ifstream restart_file;
	
	char *cstr = new char [mesh_filename.size()+1];

	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
    SU2_MPI::Error(string("There is no flow restart file ") + mesh_filename, CURRENT_FUNCTION );
  }
	
  /*--- Read the header of the file ---*/
  getline(restart_file, text_line);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		getline(restart_file, text_line);
		istringstream point_line(text_line);

    point_line >> index;
    
    if (nDim == 2) point_line >> dummy >> dummy;
    else point_line >> dummy >> dummy >> dummy;
    
		for (iVar = 0; iVar < nVar; iVar++)
			point_line >> dummy;
		for (iVar = 0; iVar < nVar; iVar++)
			point_line >> ConsVar_Res[iPoint][iVar];
	}
	restart_file.close();
}


void CGridAdaptation::GetAdjSolution(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	unsigned short iVar;
  su2double dummy;
	string text_line;
	
	string copy, mesh_filename;
	ifstream restart_file;

  /*--- Get the adjoint solution file name ---*/
	mesh_filename = config->GetSolution_AdjFileName();
  mesh_filename = config->GetObjFunc_Extension(mesh_filename);
	
	restart_file.open(mesh_filename.c_str(), ios::in);
	if (restart_file.fail()) {
    SU2_MPI::Error(string("There is no adjoint restart file ") + mesh_filename, CURRENT_FUNCTION );
  }
	
  /*--- Read the header of the file ---*/
  getline(restart_file, text_line);
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		getline(restart_file, text_line);
		istringstream point_line(text_line);
		
    point_line >> index;
    
    if (nDim == 2) point_line >> dummy >> dummy;
    else point_line >> dummy >> dummy >> dummy;
    
		for (iVar = 0; iVar < nVar; iVar ++)
			point_line >> AdjVar_Sol[iPoint][iVar];
	}
	
	restart_file.close();
}


void CGridAdaptation::GetAdjResidual(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	string text_line;
	su2double dummy;

	string mesh_filename, copy;
	ifstream restart_file;

	char buffer[50], cstr[MAX_STRING_SIZE];
	mesh_filename = config->GetSolution_AdjFileName();
	copy.assign(mesh_filename);
  unsigned short lastindex = copy.find_last_of(".");
  copy = copy.substr(0, lastindex);
	strcpy (cstr, copy.c_str());
	if (config->GetnObj() > 1) {
	  SPRINTF (buffer, "_combo.dat");
	}
	else {
    if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT)        SPRINTF (buffer, "_cd.dat");
    if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT)        SPRINTF (buffer, "_cl.dat");
    if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT)   SPRINTF (buffer, "_csf.dat");
    if (config->GetKind_ObjFunc() == INVERSE_DESIGN_PRESSURE) SPRINTF (buffer, "_invpress.dat");
    if (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX) SPRINTF (buffer, "_invheat.dat");
    if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT)    SPRINTF (buffer, "_cmx.dat");
    if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT)    SPRINTF (buffer, "_cmy.dat");
    if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT)    SPRINTF (buffer, "_cmz.dat");
    if (config->GetKind_ObjFunc() == EFFICIENCY)              SPRINTF (buffer, "_eff.dat");
    if (config->GetKind_ObjFunc() == FORCE_X_COEFFICIENT)     SPRINTF (buffer, "_cfx.dat");
    if (config->GetKind_ObjFunc() == FORCE_Y_COEFFICIENT)     SPRINTF (buffer, "_cfy.dat");
    if (config->GetKind_ObjFunc() == FORCE_Z_COEFFICIENT)     SPRINTF (buffer, "_cfz.dat");
    if (config->GetKind_ObjFunc() == TOTAL_HEATFLUX)          SPRINTF (buffer, "_totheat.dat");
    if (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX)        SPRINTF (buffer, "_maxheat.dat");
    if (config->GetKind_ObjFunc() == SURFACE_TOTAL_PRESSURE)  SPRINTF (buffer, "_pt.dat");
    if (config->GetKind_ObjFunc() == SURFACE_STATIC_PRESSURE) SPRINTF (buffer, "_pe.dat");
    if (config->GetKind_ObjFunc() == SURFACE_MASSFLOW)        SPRINTF (buffer, "_mfr.dat");
    if (config->GetKind_ObjFunc() == SURFACE_MACH)            SPRINTF (buffer, "_mach.dat");
    if (config->GetKind_ObjFunc() == CUSTOM_OBJFUNC)          SPRINTF (buffer, "_custom.dat");
	}

	strcat(cstr, buffer);
	
	restart_file.open(cstr, ios::in);
	
	if (restart_file.fail()) {
	  if (rank == MASTER_NODE)
      SU2_MPI::Error(string("There is no flow restart file ") + mesh_filename, CURRENT_FUNCTION );
  }
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		getline(restart_file, text_line);
		istringstream point_line(text_line);
    
    point_line >> index;
    
    if (nDim == 2) point_line >> dummy >> dummy;
    else point_line >> dummy >> dummy >> dummy;
    
		if (nVar == 1) point_line >> dummy >>  AdjVar_Res[iPoint][0];
		if (nVar == 4) point_line >> dummy >> dummy >> dummy >> dummy >>
									AdjVar_Res[iPoint][0] >> AdjVar_Res[iPoint][1] >> AdjVar_Res[iPoint][2] >> 
									AdjVar_Res[iPoint][3];
		if (nVar == 5) point_line >> dummy >> dummy >> dummy >> dummy >> dummy >>
									AdjVar_Res[iPoint][0] >> AdjVar_Res[iPoint][1] >> AdjVar_Res[iPoint][2] >> 
									AdjVar_Res[iPoint][3] >> AdjVar_Res[iPoint][4];
	}
	restart_file.close();
}

void CGridAdaptation::SetComplete_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {	
		geometry->elem[iElem]->SetDivide (true);
	}
}

void CGridAdaptation::SetNo_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {	
		geometry->elem[iElem]->SetDivide (false);
	}
}

void CGridAdaptation::SetWake_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	su2double Coordx, Coordy, dist, wake = 0.5;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx = geometry->node[iPoint]->GetCoord(0);
			Coordy = geometry->node[iPoint]->GetCoord(1);
			dist = sqrt(Coordx*Coordx+Coordy*Coordy);
			if (dist < wake) {
				geometry->elem[iElem]->SetDivide (true);
			}
			if ((Coordx > 0) && ((Coordy > -wake) && (Coordy < wake))) {
				geometry->elem[iElem]->SetDivide (true);
			}
		}
}

void CGridAdaptation::SetSupShock_Refinement(CGeometry *geometry, CConfig *config) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	su2double Coordx, Coordy;
	su2double mu_1 = asin(1/config->GetMach()-0.1);
	su2double mu_2 = asin(1/(config->GetMach()-0.7));
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx = geometry->node[iPoint]->GetCoord(0);
			Coordy = geometry->node[iPoint]->GetCoord(1);
			if (Coordy < 0.0)
			if ((Coordx > fabs(Coordy/tan(mu_2))-0.25) && (Coordx < fabs(Coordy/tan(mu_1))+1.25)) {
				geometry->elem[iElem]->SetDivide (true);
			}
		}
}

long CGridAdaptation::CheckRectCode(bool *AdaptCode) {
	
	int Code = -1;
	
	// Default
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true)) { Code = 0; }
	
	// Combination 1:4
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 1; return Code;}
	
	// Combination 1:2
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 3; return Code;}
				
	return Code;
	
}

long CGridAdaptation::CheckRectExtCode(bool *AdaptCode) {
	
	int Code = -1;
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true)) {Code = 0;}
		
	// Combination 1R -> 3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 4; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 5; return Code;}
	
	// Combination 1R -> 4T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 7; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 8; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 9; return Code;}
		
	// Combination 1R -> 1R+3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 12; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 13; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 14; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 15; return Code;}
		
	return Code;
	
}

long CGridAdaptation::CheckTriangleCode(bool *AdaptCode) {
	
	int Code = -1;
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true)) {Code = 0;}
	
	// Combination 1T -> 3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true)) {Code = 1; return Code;}
		
	// Combination 1T -> 2T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true)) {Code = 4; return Code;}
		
	// Combination 1T -> 1R+1T (2D) or 3T (3D)
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false)) {Code = 5; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true)) {Code = 7; return Code;}
		
	return Code;
	
}

long CGridAdaptation::CheckTetraCode(bool *AdaptCode) {
	
	int Code = -1;
	unsigned short nDivEdges, iVar;
	
	// Default
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true) || (AdaptCode[4] == true) ||
			(AdaptCode[5] == true) ) { Code = 0; }
	
	// Combination 1:8
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) ) {Code = 1; return Code;}
	
	// Combination 1:4
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) ) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) ) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) ) {Code = 4; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) ) {Code = 5; return Code;}
	
	// Combinations with 1, 2, and 3 (no regular) divided edges.
	
	nDivEdges = 0;
	for (iVar = 0; iVar < 6; iVar ++)
		if (AdaptCode[iVar] == true) nDivEdges++;
	
	if ((nDivEdges == 1) || (nDivEdges == 2) || (nDivEdges == 3)) {Code = 6; return Code;}
	
	return Code;

}

long CGridAdaptation::CheckHexaCode(bool *AdaptCode) {

	int Code = -1;
	
	// Default
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true) || (AdaptCode[4] == true) ||
			(AdaptCode[5] == true) || (AdaptCode[6] == true) || (AdaptCode[7] == true) || (AdaptCode[8] == true) || (AdaptCode[9] == true) ||
			(AdaptCode[10] == true) || (AdaptCode[11] == true)) { Code = 0; }
	
	// Combination 1:8
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == true) && (AdaptCode[7] == true) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 1; return Code;}
	
	// Combination 1:4
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == true) && (AdaptCode[7] == false) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == false) && (AdaptCode[7] == true) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == true) && (AdaptCode[7] == true) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 4; return Code;}
	
	// Combination 1:2
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == false) && (AdaptCode[7] == false) && (AdaptCode[8] == true) && (AdaptCode[9] == true) &&
			(AdaptCode[10] == true) && (AdaptCode[11] == true)) {Code = 5; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true) && (AdaptCode[4] == false) &&
			(AdaptCode[5] == true) && (AdaptCode[6] == false) && (AdaptCode[7] == true) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false) && (AdaptCode[4] == true) &&
			(AdaptCode[5] == false) && (AdaptCode[6] == true) && (AdaptCode[7] == false) && (AdaptCode[8] == false) && (AdaptCode[9] == false) &&
			(AdaptCode[10] == false) && (AdaptCode[11] == false)) {Code = 7; return Code; }
			
		return Code;

}

long CGridAdaptation::CheckPyramCode(bool *AdaptCode) {
	
	int Code = -1;
	
	if ((AdaptCode[0] == true) || (AdaptCode[1] == true) || (AdaptCode[2] == true) || (AdaptCode[3] == true)) {Code = 0;}
	
	// Combination 1P -> 1P
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 1; return Code;}
	
	// Combination 1P -> 3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 2; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 3; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 4; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 5; return Code;}
	
	// Combination 1P -> 4T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == false)) {Code = 6; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 7; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 8; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 9; return Code;}
	
	// Combination 1P -> 2P
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 10; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 11; return Code;}
	
	// Combination 1P -> 1P+3T
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == false)) {Code = 12; return Code;}
	
	if ((AdaptCode[0] == false) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 13; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == false) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 14; return Code;}
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == false) && (AdaptCode[3] == true)) {Code = 15; return Code;}
	
	
	// Combination 1P -> 4P
	
	if ((AdaptCode[0] == true) && (AdaptCode[1] == true) && (AdaptCode[2] == true) && (AdaptCode[3] == true)) {Code = 16; return Code;}
	
	return Code;

}

void CGridAdaptation::RectDivision(long code , long *nodes, long **Division, long *nPart) {
	
	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[7];	Division[0][3] = nodes[0];	Division[0][4] = nodes[4]; 
		
		Division[1][1] = nodes[8];	Division[1][2] = nodes[4];	Division[1][3] = nodes[1]; Division[1][4] = nodes[5]; 
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[5];	Division[2][3] = nodes[2]; Division[2][4] = nodes[6]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3]; Division[3][4] = nodes[7]; 
				
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[6];	Division[0][4] = nodes[3]; 
		
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2]; Division[1][4] = nodes[6]; 
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; 
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[5];	Division[0][4] = nodes[7]; 
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3]; Division[1][4] = nodes[7]; 
				
	}	
	
}

void CGridAdaptation::RectExtDivision(long code , long *nodes, long **Division, long *nPart) {
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[4];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];
		
		Division[1][1] = nodes[4];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3];
		
		Division[2][1] = nodes[4];	Division[2][2] = nodes[3];	Division[2][3] = nodes[0];
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[5];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[3];
		
		Division[2][1] = nodes[3];	Division[2][2] = nodes[5];	Division[2][3] = nodes[2]; 
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[6];
		
		Division[2][1] = nodes[0];	Division[2][2] = nodes[6];	Division[2][3] = nodes[3];
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[7];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[7];
		
		Division[2][1] = nodes[7];	Division[2][2] = nodes[2];	Division[2][3] = nodes[3];
		
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[3];
		
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[5];
		
		Division[2][1] = nodes[4];	Division[2][2] = nodes[5];	Division[2][3] = nodes[3];
		
		Division[3][1] = nodes[5];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[5];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[6];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[2];	Division[2][3] = nodes[6];
		
		Division[3][1] = nodes[0];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 8) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[7];
		
		Division[1][1] = nodes[7];	Division[1][2] = nodes[1];	Division[1][3] = nodes[6];
		
		Division[2][1] = nodes[1];	Division[2][2] = nodes[2];	Division[2][3] = nodes[6]; 
		
		Division[3][1] = nodes[7];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 9) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[7];
		
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2];
		
		Division[2][1] = nodes[7];	Division[2][2] = nodes[4];	Division[2][3] = nodes[2];
		
		Division[3][1] = nodes[7];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 12) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[6];	Division[0][4] = nodes[3];
		
		Division[1][1] = nodes[4];	Division[1][2] = nodes[1];	Division[1][3] = nodes[5];
		
		Division[2][1] = nodes[4];	Division[2][2] = nodes[5];	Division[2][3] = nodes[6];
		
		Division[3][1] = nodes[5];	Division[3][2] = nodes[2];	Division[3][3] = nodes[6];
		
	}	
	
	if (code == 13) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[5]; Division[0][4] = nodes[7];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[2];	Division[1][3] = nodes[6];
		
		Division[2][1] = nodes[7];	Division[2][2] = nodes[5];	Division[2][3] = nodes[6];
		
		Division[3][1] = nodes[7];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 14) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[4];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[6];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[4];	Division[1][3] = nodes[7];
		
		Division[2][1] = nodes[7];	Division[2][2] = nodes[4];	Division[2][3] = nodes[6];
		
		Division[3][1] = nodes[7];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3];
		
	}	
	
	if (code == 15) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[7];	Division[0][2] = nodes[5];	Division[0][3] = nodes[2]; Division[0][4] = nodes[3];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[4];	Division[1][3] = nodes[7];
		
		Division[2][1] = nodes[4];	Division[2][2] = nodes[1];	Division[2][3] = nodes[5];
		
		Division[3][1] = nodes[4];	Division[3][2] = nodes[5];	Division[3][3] = nodes[7];
		
	}	

}

void CGridAdaptation::TriangleDivision(long code , long *nodes, long *edges, long **Division, long *nPart) {
	
	if (code == 1) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4; Division[2][0] = 4; Division[3][0] = 4;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[3];	Division[0][3] = nodes[5];
		
		Division[1][1] = nodes[3];	Division[1][2] = nodes[1];	Division[1][3] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[4];	Division[2][3] = nodes[2];

		Division[3][1] = nodes[3];	Division[3][2] = nodes[4];	Division[3][3] = nodes[5];
		return;
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[3];	Division[0][3] = nodes[2];
		
		Division[1][1] = nodes[3];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2];
		return;
	
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[4];	Division[1][3] = nodes[2];
		return;
	
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 4; Division[1][0] = 4;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[5];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2];
		return;
	
	}	
		
	if (edges == NULL) {
		
		if (code == 5) {
			// number of nodes at each new element
			Division[0][0] = 5; Division[1][0] = 4;
			*nPart = 2;
			
			// nodes that compose each element
			Division[0][1] = nodes[0];	Division[0][2] = nodes[3];	Division[0][3] = nodes[4]; Division[0][4] = nodes[2];
			
			Division[1][1] = nodes[3];	Division[1][2] = nodes[1];	Division[1][3] = nodes[4];
			return;

		}	
		
		if (code == 6) {
			// number of nodes at each new element
			Division[0][0] = 5; Division[1][0] = 4;
			*nPart = 2;
			
			// nodes that compose each element
			Division[0][1] = nodes[3];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[5];
			
			Division[1][1] = nodes[0];	Division[1][2] = nodes[3];	Division[1][3] = nodes[5];
			return;

		}	
		
		if (code == 7) {
			// number of nodes at each new element
			Division[0][0] = 5; Division[1][0] = 4;
			*nPart = 2;
			
			// nodes that compose each element
			Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[4]; Division[0][4] = nodes[5];
			
			Division[1][1] = nodes[5];	Division[1][2] = nodes[4];	Division[1][3] = nodes[2];
			return;

		}
		
	}
	else {
		
		unsigned short iDiv, nDiv, edge_div[6][3], max_iVar, iVar, nElem, iElem, iNode, iTriangle;
		bool set_0, set_1, new_div, new_triangle[8][10];
		long max_edge;
				
		nDiv = 0;
		do { new_div = false;
			
			// Compute the greatest edge index
			max_edge = 0; max_iVar = 0;
			for (iVar = 0; iVar < 3; iVar ++) {
				max_edge = max(max_edge, edges[iVar]);
				if ( max_edge == edges[iVar] ) max_iVar = iVar;
			}
			
			// If the edge is divided compose the vector with the information of the division
			if (edges[max_iVar] >= 0) {
				if (max_iVar == 0) { edge_div[nDiv][0] = 3; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 1; }
				if (max_iVar == 1) { edge_div[nDiv][0] = 4; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 2; }			
				if (max_iVar == 2) { edge_div[nDiv][0] = 5; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 2; }			
				nDiv++; new_div = true;
			}
			// In order to do not repeat the egde, restart the code
			edges[max_iVar] = -1;
		} while (new_div);
		
		
		//	Inicializa
		for (iVar = 0; iVar < 3; iVar ++) new_triangle[0][iVar] = true;
		for (iVar = 3; iVar < 6; iVar ++) new_triangle[0][iVar] = false;
		
		nElem = 1;
		for (iDiv = 0; iDiv < nDiv; iDiv++) {
			short target_elem = -1;
			
			for (iElem = 0; iElem < nElem; iElem++) {
				set_0 = false; set_1 = false;
				if (new_triangle[iElem][edge_div[iDiv][1]]) set_0 = true;
				if (new_triangle[iElem][edge_div[iDiv][2]]) set_1 = true;
				if (set_0 && set_1) target_elem = iElem;
			}
			
			if (target_elem != -1) {
				for (iNode = 0; iNode < 6; iNode++)
					new_triangle[nElem][iNode] = new_triangle[target_elem][iNode];
				new_triangle[target_elem][edge_div[iDiv][0]] = true;
				new_triangle[target_elem][edge_div[iDiv][2]] = false;
				new_triangle[nElem][edge_div[iDiv][0]] = true;
				new_triangle[nElem][edge_div[iDiv][1]] = false;
				nElem++;
			}
		}
		
		*nPart = nElem;
		
		for (iTriangle = 0; iTriangle < nElem; iTriangle++) {
			Division[iTriangle][0] = 3;
			iVar = 1;
			for (iNode = 0; iNode < 6; iNode++)
				if (new_triangle[iTriangle][iNode] == true) {
					if (iNode == 0) Division[iTriangle][iVar] = nodes[0];
					if (iNode == 1) Division[iTriangle][iVar] = nodes[1];
					if (iNode == 2) Division[iTriangle][iVar] = nodes[2];
					if (iNode == 3) Division[iTriangle][iVar] = nodes[3];
					if (iNode == 4) Division[iTriangle][iVar] = nodes[4];
					if (iNode == 5) Division[iTriangle][iVar] = nodes[5];
					iVar++;
				}
			
		}

	}
}

void CGridAdaptation::TetraDivision(long code , long *nodes, long *edges, long **Division, long *nPart) {

	
	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		Division[4][0] = 5; Division[5][0] = 5; Division[6][0] = 5; Division[7][0] = 5;
		*nPart = 8;
		
		// nodes that compose each element
		Division[0][1] = nodes[6];	Division[0][2] = nodes[5];	Division[0][3] = nodes[4]; Division[0][4] = nodes[0]; 
		Division[1][1] = nodes[8];	Division[1][2] = nodes[4];	Division[1][3] = nodes[7]; Division[1][4] = nodes[1]; 
		Division[2][1] = nodes[6];	Division[2][2] = nodes[8];	Division[2][3] = nodes[9]; Division[2][4] = nodes[3]; 
		Division[3][1] = nodes[9];	Division[3][2] = nodes[7];	Division[3][3] = nodes[5]; Division[3][4] = nodes[2]; 
		Division[4][1] = nodes[6];	Division[4][2] = nodes[8];	Division[4][3] = nodes[7]; Division[4][4] = nodes[9]; 
		Division[5][1] = nodes[6];	Division[5][2] = nodes[7];	Division[5][3] = nodes[5]; Division[5][4] = nodes[9]; 
		Division[6][1] = nodes[7];	Division[6][2] = nodes[8];	Division[6][3] = nodes[6]; Division[6][4] = nodes[4]; 
		Division[7][1] = nodes[5];	Division[7][2] = nodes[7];	Division[7][3] = nodes[6]; Division[7][4] = nodes[4]; 
		
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[6];	Division[0][4] = nodes[2]; 
		Division[1][1] = nodes[6];	Division[1][2] = nodes[4];	Division[1][3] = nodes[8]; Division[1][4] = nodes[2]; 
		Division[2][1] = nodes[6];	Division[2][2] = nodes[8];	Division[2][3] = nodes[3]; Division[2][4] = nodes[2]; 
		Division[3][1] = nodes[4];	Division[3][2] = nodes[1];	Division[3][3] = nodes[8]; Division[3][4] = nodes[2]; 
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[2];	Division[0][2] = nodes[7];	Division[0][3] = nodes[9]; Division[0][4] = nodes[0]; 
		Division[1][1] = nodes[7];	Division[1][2] = nodes[8];	Division[1][3] = nodes[9]; Division[1][4] = nodes[0]; 
		Division[2][1] = nodes[7];	Division[2][2] = nodes[1];	Division[2][3] = nodes[8]; Division[2][4] = nodes[0]; 
		Division[3][1] = nodes[9];	Division[3][2] = nodes[8];	Division[3][3] = nodes[3]; Division[3][4] = nodes[0]; 
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[2];	Division[0][2] = nodes[5];	Division[0][3] = nodes[9];	Division[0][4] = nodes[1]; 
		Division[1][1] = nodes[9];	Division[1][2] = nodes[5];	Division[1][3] = nodes[6]; Division[1][4] = nodes[1]; 
		Division[2][1] = nodes[5];	Division[2][2] = nodes[0];	Division[2][3] = nodes[6]; Division[2][4] = nodes[1]; 
		Division[3][1] = nodes[9];	Division[3][2] = nodes[6];	Division[3][3] = nodes[3]; Division[3][4] = nodes[1]; 
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[4];	Division[0][3] = nodes[5]; Division[0][4] = nodes[3]; 
		Division[1][1] = nodes[5];	Division[1][2] = nodes[4];	Division[1][3] = nodes[7]; Division[1][4] = nodes[3]; 
		Division[2][1] = nodes[2];	Division[2][2] = nodes[5];	Division[2][3] = nodes[7]; Division[2][4] = nodes[3]; 
		Division[3][1] = nodes[7];	Division[3][2] = nodes[4];	Division[3][3] = nodes[1]; Division[3][4] = nodes[3]; 
		
	}	
	
	if (code == 6) {
		
		unsigned short iDiv, nDiv, edge_div[6][3], max_iVar, iVar, nElem, iElem, iNode, iTetra;
		bool set_0, set_1, new_div, new_tetra[8][10];
		long max_edge;
		
		nDiv = 0;
		do { new_div = false;
			
			// Compute the greatest node at the divided edge
			max_edge = 0; max_iVar = 0;
			for (iVar = 0; iVar < 6; iVar ++) {
				max_edge = max(max_edge, edges[iVar]);
				if ( max_edge == edges[iVar] ) max_iVar = iVar;
			}
			
			// If the edge is divided compose the vector with the information of the division
			if (edges[max_iVar] >= 0) {
				if (max_iVar == 0) { edge_div[nDiv][0] = 4; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 1; }
				if (max_iVar == 1) { edge_div[nDiv][0] = 5; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 2; }			
				if (max_iVar == 2) { edge_div[nDiv][0] = 6; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 3; }			
				if (max_iVar == 3) { edge_div[nDiv][0] = 7; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 2; }			
				if (max_iVar == 4) { edge_div[nDiv][0] = 8; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 3; }			
				if (max_iVar == 5) { edge_div[nDiv][0] = 9; edge_div[nDiv][1] = 2; edge_div[nDiv][2] = 3; }			
				nDiv++; new_div = true;
			}
			// In order to do not repeat the egde, restart the code
			edges[max_iVar] = -1;
		} while (new_div);
		
		
		//	Inicializa
		for (iVar = 0; iVar < 4; iVar ++) new_tetra[0][iVar] = true;
		for (iVar = 4; iVar < 10; iVar ++) new_tetra[0][iVar] = false;
		
		nElem = 1;
		for (iDiv = 0; iDiv < nDiv; iDiv++) {
			for (iVar = 0; iVar < 3; iVar++) {
				short target_elem = -1;
				
				for (iElem = 0; iElem < nElem; iElem++) {
					set_0 = false; set_1 = false;
					if (new_tetra[iElem][edge_div[iDiv][1]]) set_0 = true;
					if (new_tetra[iElem][edge_div[iDiv][2]]) set_1 = true;
					if (set_0 && set_1) target_elem = iElem;
				}
				
				if (target_elem != -1) {
					for (iNode = 0; iNode < 10; iNode++)
						new_tetra[nElem][iNode] = new_tetra[target_elem][iNode];
					new_tetra[target_elem][edge_div[iDiv][0]] = true;
					new_tetra[target_elem][edge_div[iDiv][2]] = false;
					new_tetra[nElem][edge_div[iDiv][0]] = true;
					new_tetra[nElem][edge_div[iDiv][1]] = false;
					nElem++;
				}
			}
		}
		
		*nPart = nElem;
		
		for (iTetra = 0; iTetra < nElem; iTetra++) {
			Division[iTetra][0] = 4;
			iVar = 1;
			for (iNode = 0; iNode < 10; iNode++)
				if (new_tetra[iTetra][iNode] == true) {
					if (iNode == 0) Division[iTetra][iVar] = nodes[0];
					if (iNode == 1) Division[iTetra][iVar] = nodes[1];
					if (iNode == 2) Division[iTetra][iVar] = nodes[2];
					if (iNode == 3) Division[iTetra][iVar] = nodes[3];
					if (iNode == 4) Division[iTetra][iVar] = nodes[4];
					if (iNode == 5) Division[iTetra][iVar] = nodes[5];
					if (iNode == 6) Division[iTetra][iVar] = nodes[6];
					if (iNode == 7) Division[iTetra][iVar] = nodes[7];
					if (iNode == 8) Division[iTetra][iVar] = nodes[8];
					if (iNode == 9) Division[iTetra][iVar] = nodes[9];
					iVar++;
				}
			
//			cout <<"Boolean "<< new_tetra[iTetra][0] <<" "<< new_tetra[iTetra][1] <<" "<< new_tetra[iTetra][2] <<" "<< new_tetra[iTetra][3] <<" "<< new_tetra[iTetra][4] 
//			<<" "<< new_tetra[iTetra][5] <<" "<< new_tetra[iTetra][6] <<" "<< new_tetra[iTetra][7] <<" "<< new_tetra[iTetra][8] <<" "<< new_tetra[iTetra][9] << endl;

//			cout <<"Nodes "<< nodes[0] <<" "<< nodes[1] <<" "<< nodes[2] <<" "<< nodes[3] <<" "<< nodes[4] 
//			<<" "<< nodes[5] <<" "<< nodes[6] <<" "<< nodes[7] <<" "<< nodes[8] <<" "<< nodes[9] << endl;
			
//			cout <<"Tets "<< Division[iTetra][0] <<" "<< Division[iTetra][1] <<" "<< Division[iTetra][2] <<" "<< Division[iTetra][3] <<" "<< Division[iTetra][4] 
//			<<" "<< Division[iTetra][5] <<" "<< Division[iTetra][6] <<" "<< Division[iTetra][7] <<" "<< Division[iTetra][8] <<" "<< Division[iTetra][9] << endl;

//			cin.get();
		}
		
	}	

}

void CGridAdaptation::HexaDivision(long code , long *nodes, long **Division, long *nPart) {

	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		Division[4][0] = 9; Division[5][0] = 9; Division[6][0] = 9; Division[7][0] = 9;
		*nPart = 8;
		
		// nodes that compose each element
		Division[0][1] = nodes[20];	Division[0][2] = nodes[8];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[26]; Division[0][6] = nodes[25]; Division[0][7] = nodes[17]; Division[0][8] = nodes[22];
		
		Division[1][1] = nodes[20];	Division[1][2] = nodes[9];	Division[1][3] = nodes[2]; Division[1][4] = nodes[10]; 
		Division[1][5] = nodes[26]; Division[1][6] = nodes[22]; Division[1][7] = nodes[18]; Division[1][8] = nodes[23];
		
		Division[2][1] = nodes[20];	Division[2][2] = nodes[10];	Division[2][3] = nodes[3]; Division[2][4] = nodes[11]; 
		Division[2][5] = nodes[26]; Division[2][6] = nodes[23]; Division[2][7] = nodes[19]; Division[2][8] = nodes[24];
		
		Division[3][1] = nodes[20];	Division[3][2] = nodes[11];	Division[3][3] = nodes[0]; Division[3][4] = nodes[8]; 
		Division[3][5] = nodes[26]; Division[3][6] = nodes[24]; Division[3][7] = nodes[16]; Division[3][8] = nodes[25];
		
		Division[4][1] = nodes[26];	Division[4][2] = nodes[25];	Division[4][3] = nodes[17];	Division[4][4] = nodes[22]; 
		Division[4][5] = nodes[21]; Division[4][6] = nodes[12]; Division[4][7] = nodes[5]; Division[4][8] = nodes[13];
		
		Division[5][1] = nodes[26];	Division[5][2] = nodes[22];	Division[5][3] = nodes[18]; Division[5][4] = nodes[23]; 
		Division[5][5] = nodes[21]; Division[5][6] = nodes[13]; Division[5][7] = nodes[6]; Division[5][8] = nodes[14];
		
		Division[6][1] = nodes[26];	Division[6][2] = nodes[23];	Division[6][3] = nodes[19]; Division[6][4] = nodes[24]; 
		Division[6][5] = nodes[21]; Division[6][6] = nodes[14]; Division[6][7] = nodes[7]; Division[6][8] = nodes[15];
		
		Division[7][1] = nodes[26];	Division[7][2] = nodes[24];	Division[7][3] = nodes[16]; Division[7][4] = nodes[25]; 
		Division[7][5] = nodes[21]; Division[7][6] = nodes[15]; Division[7][7] = nodes[4]; Division[7][8] = nodes[12];
		
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[10]; 
		Division[0][5] = nodes[25]; Division[0][17] = nodes[25]; Division[0][7] = nodes[18]; Division[0][8] = nodes[23];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[8];	Division[1][3] = nodes[10]; Division[1][4] = nodes[3]; 
		Division[1][5] = nodes[16]; Division[1][6] = nodes[25]; Division[1][7] = nodes[23]; Division[1][8] = nodes[19];
		
		Division[2][1] = nodes[25];	Division[2][2] = nodes[17];	Division[2][3] = nodes[18]; Division[2][4] = nodes[23]; 
		Division[2][5] = nodes[12]; Division[2][6] = nodes[5]; Division[2][7] = nodes[6]; Division[2][8] = nodes[14];
		
		Division[3][1] = nodes[16];	Division[3][2] = nodes[25];	Division[3][3] = nodes[23]; Division[3][4] = nodes[19]; 
		Division[3][5] = nodes[4]; Division[3][6] = nodes[12]; Division[3][7] = nodes[14]; Division[3][8] = nodes[7];
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9; 
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[11];	Division[0][2] = nodes[0];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[24]; Division[0][6] = nodes[16]; Division[0][7] = nodes[17]; Division[0][8] = nodes[22];
		
		Division[1][1] = nodes[3];	Division[1][2] = nodes[11];	Division[1][9] = nodes[2]; Division[1][4] = nodes[2]; 
		Division[1][5] = nodes[19]; Division[1][6] = nodes[24]; Division[1][7] = nodes[22]; Division[1][8] = nodes[18];
		
		Division[2][1] = nodes[24];	Division[2][2] = nodes[16];	Division[2][3] = nodes[17]; Division[2][4] = nodes[22]; 
		Division[2][5] = nodes[15]; Division[2][6] = nodes[4]; Division[2][7] = nodes[5]; Division[2][8] = nodes[13];
		
		Division[3][1] = nodes[19];	Division[3][2] = nodes[24];	Division[3][3] = nodes[22]; Division[3][4] = nodes[18]; 
		Division[3][5] = nodes[7]; Division[3][6] = nodes[15]; Division[3][7] = nodes[13]; Division[3][8] = nodes[6];
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; Division[2][0] = 9; Division[3][0] = 9;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[20];	Division[0][2] = nodes[8];	Division[0][3] = nodes[1];	Division[0][4] = nodes[9]; 
		Division[0][5] = nodes[21]; Division[0][6] = nodes[12]; Division[0][7] = nodes[5]; Division[0][8] = nodes[13];
		
		Division[1][1] = nodes[20];	Division[1][2] = nodes[9];	Division[1][3] = nodes[2]; Division[1][4] = nodes[10]; 
		Division[1][5] = nodes[21]; Division[1][6] = nodes[13]; Division[1][7] = nodes[6]; Division[1][8] = nodes[14];
		
		Division[2][1] = nodes[20];	Division[2][2] = nodes[10];	Division[2][3] = nodes[0]; Division[2][4] = nodes[8]; 
		Division[2][5] = nodes[21]; Division[2][6] = nodes[15]; Division[2][7] = nodes[4]; Division[2][8] = nodes[12];
		
		Division[3][1] = nodes[20];	Division[3][2] = nodes[10];	Division[3][3] = nodes[3]; Division[3][4] = nodes[11]; 
		Division[3][5] = nodes[21]; Division[3][6] = nodes[14]; Division[3][7] = nodes[7]; Division[3][8] = nodes[15];
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[3]; 
		Division[0][5] = nodes[16]; Division[0][6] = nodes[17]; Division[0][7] = nodes[18]; Division[0][8] = nodes[19];
		
		Division[1][1] = nodes[16];	Division[1][2] = nodes[17];	Division[1][3] = nodes[18]; Division[1][4] = nodes[19]; 
		Division[1][5] = nodes[4]; Division[1][6] = nodes[5]; Division[1][7] = nodes[6]; Division[1][8] = nodes[7];
		
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[9];	Division[0][4] = nodes[11]; 
		Division[0][5] = nodes[4]; Division[0][6] = nodes[5]; Division[0][7] = nodes[13]; Division[0][8] = nodes[15];
		
		Division[1][1] = nodes[9];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3]; Division[1][4] = nodes[11]; 
		Division[1][5] = nodes[13]; Division[1][6] = nodes[6]; Division[1][7] = nodes[7]; Division[1][8] = nodes[15];
		
	}	
	
	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 9; Division[1][0] = 9; 
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[10]; 
		Division[0][5] = nodes[12]; Division[0][5] = nodes[5]; Division[0][7] = nodes[6]; Division[0][8] = nodes[14];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[8];	Division[1][3] = nodes[10]; Division[1][4] = nodes[3]; 
		Division[1][5] = nodes[4]; Division[1][6] = nodes[12]; Division[1][7] = nodes[14]; Division[1][8] = nodes[7];
		
	}	

}

void CGridAdaptation::PyramDivision(long code , long *nodes, long **Division, long *nPart) {

	if (code == 1) {
		
		// number of nodes at each new element
		Division[0][0] = 6;
		*nPart = 1;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2];	Division[0][4] = nodes[3]; Division[0][5] = nodes[4]; 
		
	}	
	
	if (code == 2) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[5];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[2];	Division[1][3] = nodes[3]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[3];	Division[2][3] = nodes[0]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 3) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[6];	Division[1][3] = nodes[3]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[3];	Division[2][2] = nodes[6];	Division[2][3] = nodes[2]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 4) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[7]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[0];	Division[2][2] = nodes[7];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 5) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5;
		*nPart = 3;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[8];	Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[1];	Division[1][2] = nodes[2];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[2];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
	}	
	
	if (code == 6) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[3]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[6]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[6];	Division[2][3] = nodes[3]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[6];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 7) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[6];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[6];	Division[2][2] = nodes[2];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[0];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 8) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[8]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[8];	Division[1][2] = nodes[1];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[1];	Division[2][2] = nodes[2];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 9) {
		// number of nodes at each new element
		Division[0][0] = 5; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[8]; Division[0][4] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[5];	Division[2][3] = nodes[2]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[2];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 10) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[7];	Division[0][4] = nodes[3]; Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[2];  Division[1][4] = nodes[7]; Division[1][5] = nodes[4];
		
	}	
	
	if (code == 11) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6;
		*nPart = 2;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6];	Division[0][4] = nodes[8]; Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[8];	Division[1][2] = nodes[6];	Division[1][3] = nodes[2];  Division[1][4] = nodes[3]; Division[1][5] = nodes[4];
		
	}	
	
	if (code == 12) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[5];	Division[0][3] = nodes[7];	Division[0][4] = nodes[3];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[5];	Division[1][2] = nodes[1];	Division[1][3] = nodes[6]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[6];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[6];	Division[3][2] = nodes[2];	Division[3][3] = nodes[7]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 13) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[0];	Division[0][2] = nodes[1];	Division[0][3] = nodes[6]; Division[0][4] = nodes[8];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[6];	Division[1][2] = nodes[2];	Division[1][3] = nodes[7]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[6];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 14) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[5];	Division[0][2] = nodes[1];	Division[0][3] = nodes[2]; Division[0][4] = nodes[7];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[8];	Division[2][2] = nodes[5];	Division[2][3] = nodes[7]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[8];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 15) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 5; Division[2][0] = 5; Division[3][0] = 5;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[8];	Division[0][2] = nodes[6];	Division[0][3] = nodes[2]; Division[0][4] = nodes[3];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[0];	Division[1][2] = nodes[5];	Division[1][3] = nodes[8]; Division[1][4] = nodes[4];
		
		Division[2][1] = nodes[5];	Division[2][2] = nodes[1];	Division[2][3] = nodes[6]; Division[2][4] = nodes[4]; 
		
		Division[3][1] = nodes[5];	Division[3][2] = nodes[6];	Division[3][3] = nodes[8]; Division[3][4] = nodes[4]; 
		
	}	
	
	if (code == 16) {
		// number of nodes at each new element
		Division[0][0] = 6; Division[1][0] = 6; Division[2][0] = 6; Division[3][0] = 6;
		*nPart = 4;
		
		// nodes that compose each element
		Division[0][1] = nodes[9];	Division[0][2] = nodes[8];	Division[0][3] = nodes[0]; Division[0][4] = nodes[5];	Division[0][5] = nodes[4];
		
		Division[1][1] = nodes[9];	Division[1][2] = nodes[5];	Division[1][3] = nodes[1]; Division[1][4] = nodes[6]; Division[1][5] = nodes[4];
		
		Division[2][1] = nodes[9];	Division[2][2] = nodes[6];	Division[2][3] = nodes[2]; Division[2][4] = nodes[7]; Division[2][5] = nodes[4]; 
		
		Division[3][1] = nodes[9];	Division[3][2] = nodes[7];	Division[3][3] = nodes[3]; Division[3][4] = nodes[8]; Division[3][5] = nodes[4];
		
	}	

}

void CGridAdaptation::SetHomothetic_Adaptation2D(CGeometry *geometry, CPhysicalGeometry *geo_adapt, 
																							 CConfig *config) {
	
	unsigned long iPoint, iElem, iEdge, ip_0, ip_1, ip_2, ip_3, iVertex;
	unsigned short iDim, iMarker, iVar;
	long no_0 = 0, no_1 = 0, no_2 = 0, no_3 = 0;
  unsigned short nMarker_Max = config->GetnMarker_Max();

	long *TriangleAdaptCode;
	long **TriangleEdgeIndex; bool **TriangleEdgeCode; long **TriangleEdgeNode;
	
	long *RectAdaptCode;
	long **RectEdgeIndex; bool **RectEdgeCode; long **RectEdgeNode;
	long **RectElemIndex; bool **RectElemCode; long **RectElemNode;
	
	bool Restart_Flow = false;
	bool Restart_Adjoint = false;
	
	if ((config->GetKind_Adaptation() == FULL_FLOW) ||
			(config->GetKind_Adaptation() == GRAD_FLOW) ||
			(config->GetKind_Adaptation() == FULL_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
			(config->GetKind_Adaptation() == REMAINING) ||
			(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Flow = true;
	
	if ((config->GetKind_Adaptation() == FULL_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
			(config->GetKind_Adaptation() == REMAINING) ||
			(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Adjoint = true;
		
	TriangleAdaptCode = new long[geometry->GetnElem()];
	TriangleEdgeIndex = new long*[geometry->GetnElem()];
	TriangleEdgeCode = new bool*[geometry->GetnElem()];
	TriangleEdgeNode = new long*[geometry->GetnElem()];
	
	RectAdaptCode = new long[geometry->GetnElem()];
	RectEdgeIndex = new long*[geometry->GetnElem()];
	RectEdgeCode = new bool*[geometry->GetnElem()];
	RectEdgeNode = new long*[geometry->GetnElem()];
	RectElemIndex = new long*[geometry->GetnElem()];
	RectElemCode = new bool*[geometry->GetnElem()];
	RectElemNode = new long*[geometry->GetnElem()];
	
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		TriangleEdgeIndex[iElem] = new long [3];
		TriangleEdgeCode[iElem] = new bool [3];
		TriangleEdgeNode[iElem] = new long [3];
		
		RectEdgeIndex[iElem] = new long [4];
		RectEdgeCode[iElem] = new bool [4];
		RectEdgeNode[iElem] = new long [4];
		RectElemIndex[iElem] = new long [1];
		RectElemCode[iElem] = new bool [1];
		RectElemNode[iElem] = new long [1];
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			
			TriangleEdgeIndex[iElem][0] = geometry->FindEdge(ip_0, ip_1); TriangleEdgeCode[iElem][0] = false; TriangleEdgeNode[iElem][0] = -1;
			TriangleEdgeIndex[iElem][1] = geometry->FindEdge(ip_1, ip_2); TriangleEdgeCode[iElem][1] = false; TriangleEdgeNode[iElem][1] = -1;
			TriangleEdgeIndex[iElem][2] = geometry->FindEdge(ip_2, ip_0); TriangleEdgeCode[iElem][2] = false; TriangleEdgeNode[iElem][2] = -1;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			
			RectEdgeIndex[iElem][0] = geometry->FindEdge(ip_0, ip_1); RectEdgeCode[iElem][0] = false; RectEdgeNode[iElem][0] = -1;
			RectEdgeIndex[iElem][1] = geometry->FindEdge(ip_1, ip_2); RectEdgeCode[iElem][1] = false; RectEdgeNode[iElem][1] = -1;
			RectEdgeIndex[iElem][2] = geometry->FindEdge(ip_2, ip_3); RectEdgeCode[iElem][2] = false; RectEdgeNode[iElem][2] = -1;
			RectEdgeIndex[iElem][3] = geometry->FindEdge(ip_3, ip_0); RectEdgeCode[iElem][3] = false; RectEdgeNode[iElem][3] = -1;
			
			RectElemIndex[iElem][0] = iElem; RectElemCode[iElem][0] = false;	RectElemNode[iElem][0] = -1;
		}
		
	}
	
	/*--- Initial edges that are going to be divided ---*/
  
	bool *DivEdge = new bool[geometry->GetnEdge()]; 
	bool *DivElem = new bool[geometry->GetnElem()];
	
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) DivEdge[iEdge] = false;

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++)
		DivElem[iElem] = geometry->elem[iElem]->GetDivide();
	
	/*--- Set the edge division in the reactangles and in the edge list. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (DivElem[iElem] == true) {
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				for (int iIndex = 0; iIndex < 3; iIndex++) {
					DivEdge[TriangleEdgeIndex[iElem][iIndex]] = true;
					TriangleEdgeCode[iElem][iIndex] = true;
				}
			}
			if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
				for (int iIndex = 0; iIndex < 4; iIndex++) {
					DivEdge[RectEdgeIndex[iElem][iIndex]] = true;
					RectEdgeCode[iElem][iIndex] = true;
				}
			}
		}
	}
	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) 
		if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
			for (unsigned long iBoundElem = 0; iBoundElem < geometry->GetnElem_Bound(iMarker); iBoundElem++) {			
				
				ip_0 = geometry->bound[iMarker][iBoundElem]->GetNode(0); 
				ip_1 = geometry->bound[iMarker][iBoundElem]->GetNode(1); 
				
				long edge = geometry->FindEdge(ip_0, ip_1);
				
				if (DivEdge[edge]) {
					
					unsigned long iv_1 = geometry->node[ip_1]->GetVertex(iMarker);
					unsigned long ip_1_Nearfield = geometry->vertex[iMarker][iv_1]->GetDonorPoint();
					unsigned long iv_0 = geometry->node[ip_0]->GetVertex(iMarker);
					unsigned long ip_0_Nearfield = geometry->vertex[iMarker][iv_0]->GetDonorPoint();
					
					long edge_Nearfield = geometry->FindEdge(ip_0_Nearfield, ip_1_Nearfield);
					
					DivEdge[edge_Nearfield] = true;
					
				}
				
			}	
	
	/*--- We must verify that all the elements have the right edges marked ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			for (int iIndex = 0; iIndex < 3; iIndex++) {
				if (DivEdge[TriangleEdgeIndex[iElem][iIndex]] == true) {
					TriangleEdgeCode[iElem][iIndex] = true;
				}
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			for (int iIndex = 0; iIndex < 4; iIndex++) {
				if (DivEdge[RectEdgeIndex[iElem][iIndex]] == true) {
					RectEdgeCode[iElem][iIndex] = true;
				}
			}
		}
	}
	
	/*--- Only those elements that verify certain rules will be marked for hexa adaptation...
	 the others will need a new point in the middle and RectExts ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			TriangleAdaptCode[iElem] = CheckTriangleCode(TriangleEdgeCode[iElem]);
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			RectAdaptCode[iElem] = CheckRectCode(RectEdgeCode[iElem]);
			
			/*--- Set the RectAdaptCode ---*/
      
			if (RectAdaptCode[iElem] == 1) {
				RectElemCode[iElem][0] = true;
			}
		}
	}
	
	/*--- Create the new nodes on the edges, on the faces, and in the element. ---*/
  
	long *NodeAtEdges = new long[geometry->GetnEdge()];
	long *NodeAtElem = new long[geometry->GetnElem()];
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) NodeAtEdges[iEdge] = -1;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) NodeAtElem[iElem] = -1;

	nPoint_new = geometry->GetnPoint();
	
	su2double **NewNodeCoord;
	NewNodeCoord = new su2double *[4*geometry->GetnPoint()];
	for (iPoint = 0; iPoint < 4*geometry->GetnPoint(); iPoint++)
		NewNodeCoord[iPoint] = new su2double[geometry->GetnDim()];
		
	if (Restart_Flow) {
		ConsVar_Adapt = new su2double *[4*geometry->GetnPoint()];
		for (iPoint = 0; iPoint < 4*geometry->GetnPoint(); iPoint++)
			ConsVar_Adapt[iPoint] = new su2double[nVar];
	}
	
	if (Restart_Adjoint) {
		AdjVar_Adapt = new su2double *[4*geometry->GetnPoint()];
		for (iPoint = 0; iPoint < 4*geometry->GetnPoint(); iPoint++)
			AdjVar_Adapt[iPoint] = new su2double[nVar];
	}
	
	/*--- Set the value of the variables ---*/
  
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		for (iVar = 0; iVar < nVar; iVar ++) {
			if (Restart_Flow) ConsVar_Adapt[iPoint][iVar] = ConsVar_Sol[iPoint][iVar];
			if (Restart_Adjoint) AdjVar_Adapt[iPoint][iVar] = AdjVar_Sol[iPoint][iVar];
		}
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			
			for (int iIndex = 0; iIndex < 3; iIndex++) {
				
				if (TriangleEdgeCode[iElem][iIndex] == true) {
					if (NodeAtEdges[TriangleEdgeIndex[iElem][iIndex]] != -1)
						TriangleEdgeNode[iElem][iIndex] = NodeAtEdges[TriangleEdgeIndex[iElem][iIndex]];
					
					if (NodeAtEdges[TriangleEdgeIndex[iElem][iIndex]] == -1) {
						
						NodeAtEdges[TriangleEdgeIndex[iElem][iIndex]] = nPoint_new;
						TriangleEdgeNode[iElem][iIndex] = nPoint_new;
						
						/*--- Compute the coordinates of the new node ---*/
            
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1;}
						if (iIndex == 1) {no_0 = ip_1; no_1 = ip_2;}
						if (iIndex == 2) {no_0 = ip_2; no_1 = ip_0;}
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.5*(geometry->node[no_0]->GetCoord(iDim)+geometry->node[no_1]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.5 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.5 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			
			for (int iIndex = 0; iIndex < 4; iIndex++) {
				
				if (RectEdgeCode[iElem][iIndex] == true) {
					if (NodeAtEdges[RectEdgeIndex[iElem][iIndex]] != -1)
						RectEdgeNode[iElem][iIndex] = NodeAtEdges[RectEdgeIndex[iElem][iIndex]];
					
					if (NodeAtEdges[RectEdgeIndex[iElem][iIndex]] == -1) {
						
						NodeAtEdges[RectEdgeIndex[iElem][iIndex]] = nPoint_new;
						RectEdgeNode[iElem][iIndex] = nPoint_new;
						
						/*--- Compute the coordinates of the new node ---*/
            
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1;}
						if (iIndex == 1) {no_0 = ip_1; no_1 = ip_2;}
						if (iIndex == 2) {no_0 = ip_2; no_1 = ip_3;}
						if (iIndex == 3) {no_0 = ip_3; no_1 = ip_0;}
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.5*(geometry->node[no_0]->GetCoord(iDim)+geometry->node[no_1]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.5 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.5 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
			
			for (int iIndex = 0; iIndex < 1; iIndex++) {
				
				if (RectElemCode[iElem][iIndex] == true) {
					
					if (NodeAtElem[RectElemIndex[iElem][iIndex]] != -1) 
						RectElemNode[iElem][iIndex] = NodeAtElem[RectElemIndex[iElem][iIndex]];
					
					if (NodeAtElem[RectElemIndex[iElem][iIndex]] == -1) {
						NodeAtElem[RectElemIndex[iElem][iIndex]] = nPoint_new;
						RectElemNode[iElem][iIndex] = nPoint_new;
						
						/*--- Compute the coordinates of the new node ---*/
            
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_2; no_3 = ip_3;}
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.25*(geometry->node[no_0]->GetCoord(iDim) +
																										geometry->node[no_1]->GetCoord(iDim) +
																										geometry->node[no_2]->GetCoord(iDim) +
																										geometry->node[no_3]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.25 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]+ConsVar_Adapt[no_2][iVar]+ConsVar_Adapt[no_3][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.25 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]+AdjVar_Adapt[no_2][iVar]+AdjVar_Adapt[no_3][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
		}
			
	}

	
	/*--- if Quadrilateral adapt code equals 0, then a semidivision is applied  ---*/
	long nSemiDivided = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			if (RectAdaptCode[iElem] == 0)
				nSemiDivided++;
		}
	}
	
	/*--- If semidivision, then divide add a new point, divide the quadrilateral into triangles,
   and find the right combination, it also create the new node (hexa).  ---*/
	long nRectExt = nSemiDivided;
	
	long *RectExtAdaptCode;
	long **RectExtNode;
	long **RectRectExtIndex;
	long *RectExtRectIndex;
	
	long **RectExtEdgeIndex;
	bool **RectExtEdgeCode;
	long **RectExtEdgeNode;
	
	RectExtAdaptCode = new long [nRectExt];
	RectExtNode = new long *[nRectExt];
	RectRectExtIndex = new long *[geometry->GetnElem()];
	RectExtRectIndex = new long [nRectExt];
	RectExtEdgeIndex = new long *[nRectExt];
	RectExtEdgeCode = new bool *[nRectExt];
	RectExtEdgeNode = new long *[nRectExt];
	
	for (long iRectExt = 0; iRectExt < nRectExt; iRectExt++) {
		RectExtNode[iRectExt] = new long [4];
		RectExtEdgeIndex[iRectExt] = new long [4];
		RectExtEdgeCode[iRectExt] = new bool [4];
		RectExtEdgeNode[iRectExt] = new long [4];
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			RectRectExtIndex[iElem] = new long [1];
		}
	}
	
	nRectExt = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			if (RectAdaptCode[iElem] == 0) {
				
				/*--- Write the edge combination on the base. ---*/
        
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				ip_3 = geometry->elem[iElem]->GetNode(3);
				
				/*--- Create the 1st RectExtid. ---*/
        
				RectRectExtIndex[iElem][0] = nRectExt; RectExtRectIndex[nRectExt] = iElem; 
				
				RectExtNode[nRectExt][0] = ip_0; RectExtNode[nRectExt][1] = ip_1;
				RectExtNode[nRectExt][2] = ip_2; RectExtNode[nRectExt][3] = ip_3;
				
				RectExtEdgeIndex[nRectExt][0] = RectEdgeIndex[iElem][0]; RectExtEdgeIndex[nRectExt][1] = RectEdgeIndex[iElem][1];
				RectExtEdgeIndex[nRectExt][2] = RectEdgeIndex[iElem][2]; RectExtEdgeIndex[nRectExt][3] = RectEdgeIndex[iElem][3];
				
				RectExtEdgeNode[nRectExt][0] = RectEdgeNode[iElem][0]; RectExtEdgeNode[nRectExt][1] = RectEdgeNode[iElem][1];
				RectExtEdgeNode[nRectExt][2] = RectEdgeNode[iElem][2]; RectExtEdgeNode[nRectExt][3] = RectEdgeNode[iElem][3];
				
				
				for (int iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[RectExtEdgeIndex[nRectExt][iIndex]] == true) RectExtEdgeCode[nRectExt][iIndex] = true; 
				nRectExt++;
				
			}
		}
	}
	
	/*--- Check the kind of RectExt partitioning that should be applied ---*/
  
	for (int iRectExt = 0; iRectExt < nRectExt; iRectExt ++) {
		RectExtAdaptCode[iRectExt] = CheckRectExtCode(RectExtEdgeCode[iRectExt]);
		if (RectExtAdaptCode[iRectExt] == 0) cout << "There is a problem with one RectExt" << endl;
	}
	
	/*--- Create new structure ---*/
  
	nElem_new = geometry->GetnElem();
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			if (TriangleAdaptCode[iElem] == 1) nElem_new = nElem_new + 3;
			if (TriangleAdaptCode[iElem] == 2) nElem_new = nElem_new + 1;
			if (TriangleAdaptCode[iElem] == 3) nElem_new = nElem_new + 1;
			if (TriangleAdaptCode[iElem] == 4) nElem_new = nElem_new + 1;
			if (TriangleAdaptCode[iElem] == 5) nElem_new = nElem_new + 1;
			if (TriangleAdaptCode[iElem] == 6) nElem_new = nElem_new + 1;
			if (TriangleAdaptCode[iElem] == 7) nElem_new = nElem_new + 1;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			if (RectAdaptCode[iElem] == 1) nElem_new = nElem_new + 3;
			if (RectAdaptCode[iElem] == 2) nElem_new = nElem_new + 1;
			if (RectAdaptCode[iElem] == 3) nElem_new = nElem_new + 1;
			if (RectAdaptCode[iElem] == 0) {
				long iRectExt = RectRectExtIndex[iElem][0];
				if (RectExtAdaptCode[iRectExt] == 2) nElem_new = nElem_new + 2;
				if (RectExtAdaptCode[iRectExt] == 3) nElem_new = nElem_new + 2;
				if (RectExtAdaptCode[iRectExt] == 4) nElem_new = nElem_new + 2;
				if (RectExtAdaptCode[iRectExt] == 5) nElem_new = nElem_new + 2;
				if (RectExtAdaptCode[iRectExt] == 6) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 7) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 8) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 9) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 12) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 13) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 14) nElem_new = nElem_new + 3;
				if (RectExtAdaptCode[iRectExt] == 15) nElem_new = nElem_new + 3;
			}
		}
	}
	
	/*--- New points ---*/
  
	geo_adapt->node = new CPoint*[nPoint_new];
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0), geometry->node[iPoint]->GetCoord(1), iPoint, config);

	for (iPoint = geometry->GetnPoint(); iPoint < nPoint_new; iPoint++)
		geo_adapt->node[iPoint] = new CPoint(NewNodeCoord[iPoint][0], NewNodeCoord[iPoint][1], iPoint, config);
	
	/*--- New elements ---*/
  
	geo_adapt->elem = new CPrimalGrid*[nElem_new];
	
	unsigned long iElemNew = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			if (TriangleAdaptCode[iElem] == -1) {
				geo_adapt->elem[iElemNew] = new CTriangle(geometry->elem[iElem]->GetNode(0), 
																							 geometry->elem[iElem]->GetNode(1), 
																							 geometry->elem[iElem]->GetNode(2), 2);
				iElemNew++;
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
			if (RectAdaptCode[iElem] == -1) {
				geo_adapt->elem[iElemNew] = new CQuadrilateral(geometry->elem[iElem]->GetNode(0),
																								geometry->elem[iElem]->GetNode(1),
																								geometry->elem[iElem]->GetNode(2),
																								geometry->elem[iElem]->GetNode(3), 2);
				iElemNew++;
			}
		}
	}
			
	long nodes[27]; long **Division; long nPart;
	
	Division = new long*[100];
	for (long iVar = 0; iVar < 100; iVar++)
		Division[iVar] = new long[100];
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			
			/*--- Triangle elements... ---*/
      
			if (TriangleAdaptCode[iElem] > 0) {
				
				/*--- First the corners ---*/
        
				nodes[0] = geometry->elem[iElem]->GetNode(0);		
				nodes[1] = geometry->elem[iElem]->GetNode(1); 
				nodes[2] = geometry->elem[iElem]->GetNode(2);
				
				/*--- Next the points that correspond to the broken edges. ---*/
        
				nodes[3] = TriangleEdgeNode[iElem][0]; 
				nodes[4] = TriangleEdgeNode[iElem][1];
				nodes[5] = TriangleEdgeNode[iElem][2];
				
				TriangleDivision(TriangleAdaptCode[iElem], nodes, NULL, Division, &nPart);
				for (long iPart = 0; iPart < nPart; iPart++) {
					
					/*--- Triangle case ---*/
          
					if (Division[iPart][0] == 4) {
						geo_adapt->elem[iElemNew] = new CTriangle(Division[iPart][1], 
																											Division[iPart][2], 
																											Division[iPart][3], 2);
						iElemNew++;
					}
					
					/*--- Quadrilateral case ---*/
          
					if (Division[iPart][0] == 5) {
						geo_adapt->elem[iElemNew] = new CQuadrilateral(Division[iPart][1],
																											 Division[iPart][2], 
																											 Division[iPart][3], 
																											 Division[iPart][4], 2);
						iElemNew++;
					}
				}
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
      
			/*--- Rect elements... ---*/
      
			if (RectAdaptCode[iElem] > 0) {
				
				/*--- First the corners ---*/
        
				nodes[0] = geometry->elem[iElem]->GetNode(0);		nodes[1] = geometry->elem[iElem]->GetNode(1);
				nodes[2] = geometry->elem[iElem]->GetNode(2);		nodes[3] = geometry->elem[iElem]->GetNode(3);
				
				/*--- Next the points that correspond to the broken edges. ---*/
        
				nodes[4] = RectEdgeNode[iElem][0]; nodes[5] = RectEdgeNode[iElem][1];
				nodes[6] = RectEdgeNode[iElem][2]; nodes[7] = RectEdgeNode[iElem][3];
				
				/*--- Next the points that correspond to the element. ---*/
        
				nodes[8] = RectElemNode[iElem][0];
				
				RectDivision(RectAdaptCode[iElem], nodes, Division, &nPart);
				for (long iPart = 0; iPart < nPart; iPart++) {
					geo_adapt->elem[iElemNew] = new CQuadrilateral(Division[iPart][1],
																										 Division[iPart][2], 
																										 Division[iPart][3], 
																										 Division[iPart][4], 2);
					iElemNew++;
				}
			}
			
			/*--- RectExt elements... ---*/
      
			if (RectAdaptCode[iElem] == 0) {
				long iRectExt = RectRectExtIndex[iElem][0];
				
				/*--- First the corners ---*/
        
				nodes[0] = RectExtNode[iRectExt][0];	
				nodes[1] = RectExtNode[iRectExt][1];	
				nodes[2] = RectExtNode[iRectExt][2];	
				nodes[3] = RectExtNode[iRectExt][3];	
				
				/*--- Next the points that correspond to the broken edges. ---*/
        
				nodes[4] = RectExtEdgeNode[iRectExt][0];
				nodes[5] = RectExtEdgeNode[iRectExt][1];
				nodes[6] = RectExtEdgeNode[iRectExt][2];
				nodes[7] = RectExtEdgeNode[iRectExt][3];
				
				RectExtDivision(RectExtAdaptCode[iRectExt], nodes, Division, &nPart);
				for (long iPart = 0; iPart < nPart; iPart++) {
					
					/*--- Triangle case ---*/
          
					if (Division[iPart][0] == 4) {
						geo_adapt->elem[iElemNew] = new CTriangle(Division[iPart][1], 
																											Division[iPart][2], 
																											Division[iPart][3], 2);
						iElemNew++;
					}
					
					/*--- Quadrilateral case ---*/
          
					if (Division[iPart][0] == 5) {
						geo_adapt->elem[iElemNew] = new CQuadrilateral(Division[iPart][1],
																											 Division[iPart][2], 
																											 Division[iPart][3], 
																											 Division[iPart][4], 2);
						iElemNew++;
					}
				}
			}
		}
	}
		
	geo_adapt->SetnElem(nElem_new);
	geo_adapt->SetnPoint(nPoint_new);
	geo_adapt->SetnPointDomain(nPoint_new);
	geo_adapt->SetnDim(nDim);
	
	/*--- Create boundary structure ---*/
  
	geo_adapt->SetnMarker(geometry->GetnMarker());
	geo_adapt->nElem_Bound = new unsigned long [geometry->GetnMarker()];
	geo_adapt->Tag_to_Marker = new string [nMarker_Max];
	geo_adapt->bound = new CPrimalGrid**[geometry->GetnMarker()];
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		long nNewBCcv = 0;
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {			
			ip_0 = geometry->bound[iMarker][iVertex]->GetNode(0);
			ip_1 = geometry->bound[iMarker][iVertex]->GetNode(1);
			if (DivEdge[geometry->FindEdge(ip_0, ip_1)]) nNewBCcv = nNewBCcv + 2;
			else nNewBCcv = nNewBCcv + 1;
		}
		geo_adapt->bound[iMarker] = new CPrimalGrid* [nNewBCcv];
		geo_adapt->SetnElem_Bound(iMarker, nNewBCcv);
		geo_adapt->SetMarker_Tag(iMarker, geometry->GetMarker_Tag(iMarker));
	}
	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		long nNewBCcv = 0;
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {			
			
			ip_0 = geometry->bound[iMarker][iVertex]->GetNode(0); geo_adapt->node[ip_0]->SetBoundary(geometry->GetnMarker());
			ip_1 = geometry->bound[iMarker][iVertex]->GetNode(1); geo_adapt->node[ip_1]->SetBoundary(geometry->GetnMarker());
			long ip_01 = NodeAtEdges[geometry->FindEdge(ip_0, ip_1)];
      
			if (ip_01 != -1) {
        
        
//        if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
//          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ||
//          (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)) {
//          
//          /*--- Recompute the coordinates using the NACA 4Digits analytical definition ---*/
//          
//          su2double Ya = 0.0 / 100.0; /*--- Maximum camber as a fraction of the chord
//                                    (100 m is the first of the four digits) ---*/
//          su2double Xa = 0.0 / 10.0; /*--- Location of maximum camber as a fraction of
//                                   the chord (10 p is the second digit in the NACA xxxx description) ---*/
//          su2double t = 12.0 / 100.0; /*--- Maximum thickness as a fraction of the
//                                    chord (so 100 t gives the last two digits in
//                                    the NACA 4-digit denomination) ---*/
//          
//          su2double *Coord = geo_adapt->node[ip_01]->GetCoord();
//          su2double *Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
//          
//          su2double Ycurv = 0.0;
//          if (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
//          else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
//          
//          su2double Yesp = 0.0;
//          Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
//                    1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
//          
//          if (Normal[1] > 0) Coord[1] = (Ycurv + Yesp);
//          if (Normal[1] < 0) Coord[1] = (Ycurv - Yesp);
//          
//        }
        
				geo_adapt->node[ip_01]->SetBoundary(geometry->GetnMarker());
				geo_adapt->bound[iMarker][nNewBCcv] = new CLine(ip_0, ip_01, 2);
				nNewBCcv++;
				geo_adapt->bound[iMarker][nNewBCcv] = new CLine(ip_01, ip_1, 2);
				nNewBCcv++;
			}
			else {
				geo_adapt->bound[iMarker][nNewBCcv] = new CLine(ip_0, ip_1, 2);
				nNewBCcv++;
			}
		}
	}
	
	delete [] DivEdge; 
	delete [] DivElem; 
	delete [] NodeAtEdges; 
	delete [] NodeAtElem; 
	
}

void CGridAdaptation::SetHomothetic_Adaptation3D(CGeometry *geometry, CPhysicalGeometry *geo_adapt, 
																								 CConfig *config) {

	unsigned long iPoint, iElem, iEdge, ip_0, ip_1, ip_2, ip_3, ip_4, ip_5, ip_6, ip_7, iVertex;
	unsigned short iDim, iMarker, iVar;
	long no_0 = 0, no_1 = 0, no_2 = 0, no_3 = 0, no_4 = 0, no_5 = 0, no_6 = 0, no_7 = 0;
	unsigned short counter;
  unsigned short nMarker_Max = config->GetnMarker_Max();

	long *TetraAdaptCode; 
	long **TetraEdgeIndex; bool **TetraEdgeCode; long **TetraEdgeNode;
	
	long *HexaAdaptCode; 
	long **HexaEdgeIndex; bool **HexaEdgeCode; long **HexaEdgeNode;
	long **HexaFaceIndex; bool **HexaFaceCode; long **HexaFaceNode;
	long **HexaElemIndex; bool **HexaElemCode; long **HexaElemNode;
	
	bool Restart_Flow = false;
	bool Restart_Adjoint = false;
	
	if ((config->GetKind_Adaptation() == FULL_FLOW) ||
			(config->GetKind_Adaptation() == GRAD_FLOW) ||
			(config->GetKind_Adaptation() == FULL_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
			(config->GetKind_Adaptation() == REMAINING) ||
			(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Flow = true;
	
	if ((config->GetKind_Adaptation() == FULL_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
			(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
			(config->GetKind_Adaptation() == REMAINING) ||
			(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Adjoint = true;
		
	TetraAdaptCode = new long[geometry->GetnElem()];
	TetraEdgeIndex = new long*[geometry->GetnElem()];
	TetraEdgeCode = new bool*[geometry->GetnElem()];
	TetraEdgeNode = new long*[geometry->GetnElem()];
	
	HexaAdaptCode = new long[geometry->GetnElem()];
	HexaEdgeIndex = new long*[geometry->GetnElem()];
	HexaEdgeCode = new bool*[geometry->GetnElem()];
	HexaEdgeNode = new long*[geometry->GetnElem()];
	HexaFaceIndex = new long*[geometry->GetnElem()];
	HexaFaceCode = new bool*[geometry->GetnElem()];
	HexaFaceNode = new long*[geometry->GetnElem()];
	HexaElemIndex = new long*[geometry->GetnElem()];
	HexaElemCode = new bool*[geometry->GetnElem()];
	HexaElemNode = new long*[geometry->GetnElem()];
		
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		TetraEdgeIndex[iElem] = new long [6];
		TetraEdgeCode[iElem] = new bool [6];
		TetraEdgeNode[iElem] = new long [6];
		
		HexaEdgeIndex[iElem] = new long [12];
		HexaEdgeCode[iElem] = new bool [12];
		HexaEdgeNode[iElem] = new long [12];
		HexaFaceIndex[iElem] = new long [6];
		HexaFaceCode[iElem] = new bool [6];
		HexaFaceNode[iElem] = new long [6];
		HexaElemIndex[iElem] = new long [1];
		HexaElemCode[iElem] = new bool [1];
		HexaElemNode[iElem] = new long [1];		
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			
			TetraEdgeIndex[iElem][0] = geometry->FindEdge(ip_0, ip_1); TetraEdgeCode[iElem][0] = false; TetraEdgeNode[iElem][0] = -1;
			TetraEdgeIndex[iElem][1] = geometry->FindEdge(ip_0, ip_2); TetraEdgeCode[iElem][1] = false; TetraEdgeNode[iElem][1] = -1;
			TetraEdgeIndex[iElem][2] = geometry->FindEdge(ip_0, ip_3); TetraEdgeCode[iElem][2] = false; TetraEdgeNode[iElem][2] = -1;
			TetraEdgeIndex[iElem][3] = geometry->FindEdge(ip_1, ip_2); TetraEdgeCode[iElem][3] = false; TetraEdgeNode[iElem][3] = -1;
			TetraEdgeIndex[iElem][4] = geometry->FindEdge(ip_1, ip_3); TetraEdgeCode[iElem][4] = false; TetraEdgeNode[iElem][4] = -1;
			TetraEdgeIndex[iElem][5] = geometry->FindEdge(ip_2, ip_3); TetraEdgeCode[iElem][5] = false; TetraEdgeNode[iElem][5] = -1;
			
		}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			ip_4 = geometry->elem[iElem]->GetNode(4);
			ip_5 = geometry->elem[iElem]->GetNode(5);
			ip_6 = geometry->elem[iElem]->GetNode(6);
			ip_7 = geometry->elem[iElem]->GetNode(7);
			
			HexaEdgeIndex[iElem][0] = geometry->FindEdge(ip_0, ip_1);		HexaEdgeCode[iElem][0] = false;		HexaEdgeNode[iElem][0] = -1;
			HexaEdgeIndex[iElem][1] = geometry->FindEdge(ip_1, ip_2);		HexaEdgeCode[iElem][1] = false;		HexaEdgeNode[iElem][1] = -1;
			HexaEdgeIndex[iElem][2] = geometry->FindEdge(ip_2, ip_3);		HexaEdgeCode[iElem][2] = false;		HexaEdgeNode[iElem][2] = -1;
			HexaEdgeIndex[iElem][3] = geometry->FindEdge(ip_3, ip_0);		HexaEdgeCode[iElem][3] = false;		HexaEdgeNode[iElem][3] = -1;
			HexaEdgeIndex[iElem][4] = geometry->FindEdge(ip_4, ip_5);		HexaEdgeCode[iElem][4] = false;		HexaEdgeNode[iElem][4] = -1;
			HexaEdgeIndex[iElem][5] = geometry->FindEdge(ip_5, ip_6);		HexaEdgeCode[iElem][5] = false;		HexaEdgeNode[iElem][5] = -1;
			HexaEdgeIndex[iElem][6] = geometry->FindEdge(ip_6, ip_7);		HexaEdgeCode[iElem][6] = false;		HexaEdgeNode[iElem][6] = -1;
			HexaEdgeIndex[iElem][7] = geometry->FindEdge(ip_7, ip_4);		HexaEdgeCode[iElem][7] = false;		HexaEdgeNode[iElem][7] = -1;
			HexaEdgeIndex[iElem][8] = geometry->FindEdge(ip_0, ip_4);		HexaEdgeCode[iElem][8] = false;		HexaEdgeNode[iElem][8] = -1;
			HexaEdgeIndex[iElem][9] = geometry->FindEdge(ip_1, ip_5);		HexaEdgeCode[iElem][9] = false;		HexaEdgeNode[iElem][9] = -1;
			HexaEdgeIndex[iElem][10] = geometry->FindEdge(ip_2, ip_6);		HexaEdgeCode[iElem][10] = false;	HexaEdgeNode[iElem][10] = -1;	
			HexaEdgeIndex[iElem][11] = geometry->FindEdge(ip_3, ip_7);		HexaEdgeCode[iElem][11] = false;	HexaEdgeNode[iElem][11] = -1;
			
//			HexaFaceIndex[iElem][0] = geometry->FindFace(iElem, ip_0, ip_1, ip_2, ip_3);		HexaFaceCode[iElem][0] = false;		HexaFaceNode[iElem][0] = -1;
//			HexaFaceIndex[iElem][1] = geometry->FindFace(iElem, ip_4, ip_5, ip_6, ip_7);		HexaFaceCode[iElem][1] = false;		HexaFaceNode[iElem][1] = -1;
//			HexaFaceIndex[iElem][2] = geometry->FindFace(iElem, ip_1, ip_2, ip_6, ip_5);		HexaFaceCode[iElem][2] = false;		HexaFaceNode[iElem][2] = -1;
//			HexaFaceIndex[iElem][3] = geometry->FindFace(iElem, ip_3, ip_2, ip_6, ip_7);		HexaFaceCode[iElem][3] = false;		HexaFaceNode[iElem][3] = -1;
//			HexaFaceIndex[iElem][4] = geometry->FindFace(iElem, ip_0, ip_3, ip_7, ip_4);		HexaFaceCode[iElem][4] = false;		HexaFaceNode[iElem][4] = -1;
//			HexaFaceIndex[iElem][5] = geometry->FindFace(iElem, ip_0, ip_1, ip_5, ip_4);		HexaFaceCode[iElem][5] = false;		HexaFaceNode[iElem][5] = -1;
			
			HexaElemIndex[iElem][0] = iElem; HexaElemCode[iElem][0] = false;	HexaElemNode[iElem][0] = -1;

		}
		
	}

	/*--- Remove pyramids and prisms in the adaptation process ---*/
	unsigned short iFace, iNode, ElemIndex;
	long jElem;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if ((geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) || 
				(geometry->elem[iElem]->GetVTK_Type() == PYRAMID) || 
				(geometry->elem[iElem]->GetVTK_Type() == PRISM)) {
			geometry->elem[iElem]->SetDivide(false);
			for (iFace = 0; iFace < geometry->elem[iElem]->GetnFaces(); iFace++)
				for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodesFace(iFace); iNode++) {
					iPoint = geometry->elem[iElem]->GetNode(geometry->elem[iElem]->GetFaces(iFace, iNode));
					for (ElemIndex = 0; ElemIndex < geometry->node[iPoint]->GetnElem(); ElemIndex++) {
						jElem = geometry->node[iPoint]->GetElem(ElemIndex);
						if (jElem != -1) geometry->elem[jElem]->SetDivide(false);
					}
				}
		}
	}
	
	/*--- Do not addapt the boundaries ---*/
	if (!config->GetAdaptBoundary()) {
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {			
				for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
					iPoint = geometry->bound[iMarker][iVertex]->GetNode(iNode);
					for (ElemIndex = 0; ElemIndex < geometry->node[iPoint]->GetnElem(); ElemIndex++) {
						jElem = geometry->node[iPoint]->GetElem(ElemIndex);
						if (jElem != -1) geometry->elem[jElem]->SetDivide(false);
					}
				}
			}
		}
	}
	
	/*--- Initial edges that are going to be divided ---*/
	bool *DivEdge = new bool[geometry->GetnEdge()]; 
	bool *DivElem = new bool[geometry->GetnElem()];
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) DivEdge[iEdge] = false;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++)
		DivElem[iElem] = geometry->elem[iElem]->GetDivide();
	
	/*--- Set the edge division in the reactangles and in the edge list ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (DivElem[iElem] == true) {
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				for (long iIndex = 0; iIndex < 6; iIndex++) {
					DivEdge[TetraEdgeIndex[iElem][iIndex]] = true;
					TetraEdgeCode[iElem][iIndex] = true;
				}
			}
			if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
				for (long iIndex = 0; iIndex < 12; iIndex++) {
					DivEdge[HexaEdgeIndex[iElem][iIndex]] = true;
					HexaEdgeCode[iElem][iIndex] = true;
				}
			}
		}
	}
	
	/*--- Only with tets, check if an element have more than 4 divided edges, the element should be completelly divided ---*/
	bool new_elem;	
	do { new_elem = false;
		for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				ip_3 = geometry->elem[iElem]->GetNode(3);
				
				counter = 0;
				if (DivEdge[TetraEdgeIndex[iElem][0]]) counter++;
				if (DivEdge[TetraEdgeIndex[iElem][1]]) counter++;
				if (DivEdge[TetraEdgeIndex[iElem][2]]) counter++;
				if (DivEdge[TetraEdgeIndex[iElem][3]]) counter++;
				if (DivEdge[TetraEdgeIndex[iElem][4]]) counter++;
				if (DivEdge[TetraEdgeIndex[iElem][5]]) counter++;
								
				if ((counter > 3) && (!DivElem[iElem])) { 
					DivEdge[geometry->FindEdge(ip_0, ip_1)] = true;
					DivEdge[geometry->FindEdge(ip_0, ip_2)] = true;
					DivEdge[geometry->FindEdge(ip_0, ip_3)] = true;
					DivEdge[geometry->FindEdge(ip_1, ip_2)] = true;
					DivEdge[geometry->FindEdge(ip_1, ip_3)] = true;
					DivEdge[geometry->FindEdge(ip_2, ip_3)] = true;
					DivElem[iElem] = true;
					new_elem = true;
				}				
			}
		}
	} while (new_elem);
		
	/*--- We must verify that all the elements have the right edges marked, 
	 with tets 4 and 5 edges are not allowed ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			for (long iIndex = 0; iIndex < 6; iIndex++) {
				if (DivEdge[TetraEdgeIndex[iElem][iIndex]] == true) {
					TetraEdgeCode[iElem][iIndex] = true;
				}
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			for (long iIndex = 0; iIndex < 12; iIndex++) {
				if (DivEdge[HexaEdgeIndex[iElem][iIndex]] == true) {
					HexaEdgeCode[iElem][iIndex] = true;
				}
			}
		}
	}
	
	// Only those elements that verify certain rules will be marked for hexa adaptation... 
	// the others will need a new point in the middle and Pyrams
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
		
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			TetraAdaptCode[iElem] = CheckTetraCode(TetraEdgeCode[iElem]);
//			if (TetraAdaptCode[iElem] != -1) cout << TetraAdaptCode[iElem] <<" "<< iElem << endl;
		}
		
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			
			HexaAdaptCode[iElem] = CheckHexaCode(HexaEdgeCode[iElem]);
//			if (HexaAdaptCode[iElem] != -1) cout << HexaAdaptCode[iElem] <<" "<< iElem << endl;
			
			// Set the HexaFaceCode, and HexaElemCode
			if (HexaAdaptCode[iElem] == 1) {
				HexaFaceCode[iElem][0] = true; HexaFaceCode[iElem][1] = true; HexaFaceCode[iElem][2] = true; 
				HexaFaceCode[iElem][3] = true; HexaFaceCode[iElem][4] = true; HexaFaceCode[iElem][5] = true;
				HexaElemCode[iElem][0] = true;
			}
			if (HexaAdaptCode[iElem] == 2) {
				HexaFaceCode[iElem][3] = true; HexaFaceCode[iElem][5] = true;
			}
			if (HexaAdaptCode[iElem] == 3) {
				HexaFaceCode[iElem][2] = true; HexaFaceCode[iElem][4] = true;
			}
			if (HexaAdaptCode[iElem] == 4) {
				HexaFaceCode[iElem][0] = true; HexaFaceCode[iElem][1] = true;
			}
		}
	}
		
	// Create the new nodes on the edges, on the faces, and in the element.
	long *NodeAtEdges = new long[geometry->GetnEdge()];
//	long *NodeAtElem = new long[geometry->GetnFace()];
	long *NodeAtElem = new long[geometry->GetnElem()];
	
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) NodeAtEdges[iEdge] = -1;
//	for (iFace = 0; iFace < geometry->GetnFace(); iFace++) NodeAtFaces[iEdge] = -1;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) NodeAtElem[iElem] = -1;
	
	nPoint_new = geometry->GetnPoint();
	
	su2double **NewNodeCoord;
	NewNodeCoord = new su2double *[10*geometry->GetnPoint()];
	for (iPoint = 0; iPoint < 10*geometry->GetnPoint(); iPoint++)
		NewNodeCoord[iPoint] = new su2double[geometry->GetnDim()];
	
	if (Restart_Flow) {
		ConsVar_Adapt = new su2double *[10*geometry->GetnPoint()];
		for (iPoint = 0; iPoint < 10*geometry->GetnPoint(); iPoint++)
			ConsVar_Adapt[iPoint] = new su2double[nVar];
	}
	
	if (Restart_Adjoint) {
		AdjVar_Adapt = new su2double *[10*geometry->GetnPoint()];
		for (iPoint = 0; iPoint < 10*geometry->GetnPoint(); iPoint++)
			AdjVar_Adapt[iPoint] = new su2double[nVar];
	}
	
	// Set the value of the variables
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		for (iVar = 0; iVar < nVar; iVar ++) {
			if (Restart_Flow) ConsVar_Adapt[iPoint][iVar] = ConsVar_Sol[iPoint][iVar];
			if (Restart_Adjoint) AdjVar_Adapt[iPoint][iVar] = AdjVar_Sol[iPoint][iVar];
		}
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			
			for (long iIndex = 0; iIndex < 6; iIndex++) {
				
				if (TetraEdgeCode[iElem][iIndex] == true) {
					if (NodeAtEdges[TetraEdgeIndex[iElem][iIndex]] != -1)
						TetraEdgeNode[iElem][iIndex] = NodeAtEdges[TetraEdgeIndex[iElem][iIndex]];
					
					if (NodeAtEdges[TetraEdgeIndex[iElem][iIndex]] == -1) {
						
						NodeAtEdges[TetraEdgeIndex[iElem][iIndex]] = nPoint_new;
						TetraEdgeNode[iElem][iIndex] = nPoint_new;
						
						// Compute the coordinates of the new node					
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1;}
						if (iIndex == 1) {no_0 = ip_0; no_1 = ip_2;}
						if (iIndex == 2) {no_0 = ip_0; no_1 = ip_3;}
						if (iIndex == 3) {no_0 = ip_1; no_1 = ip_2;}
						if (iIndex == 4) {no_0 = ip_1; no_1 = ip_3;}			
						if (iIndex == 5) {no_0 = ip_2; no_1 = ip_3;}
						
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.5*(geometry->node[no_0]->GetCoord(iDim)+geometry->node[no_1]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.5 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.5 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
		}
		
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			ip_3 = geometry->elem[iElem]->GetNode(3);
			ip_4 = geometry->elem[iElem]->GetNode(4);
			ip_5 = geometry->elem[iElem]->GetNode(5);
			ip_6 = geometry->elem[iElem]->GetNode(6);
			ip_7 = geometry->elem[iElem]->GetNode(7);
			
			for (long iIndex = 0; iIndex < 12; iIndex++) {
				
				if (HexaEdgeCode[iElem][iIndex] == true) {
					if (NodeAtEdges[HexaEdgeIndex[iElem][iIndex]] != -1)
						HexaEdgeNode[iElem][iIndex] = NodeAtEdges[HexaEdgeIndex[iElem][iIndex]];
					
					if (NodeAtEdges[HexaEdgeIndex[iElem][iIndex]] == -1) {
						
						NodeAtEdges[HexaEdgeIndex[iElem][iIndex]] = nPoint_new;
						HexaEdgeNode[iElem][iIndex] = nPoint_new;
						
						// Compute the coordinates of the new node					
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1;}
						if (iIndex == 1) {no_0 = ip_1; no_1 = ip_2;}
						if (iIndex == 2) {no_0 = ip_2; no_1 = ip_3;}
						if (iIndex == 3) {no_0 = ip_3; no_1 = ip_0;}
						if (iIndex == 4) {no_0 = ip_4; no_1 = ip_5;}
						if (iIndex == 5) {no_0 = ip_5; no_1 = ip_6;}
						if (iIndex == 6) {no_0 = ip_6; no_1 = ip_7;}
						if (iIndex == 7) {no_0 = ip_7; no_1 = ip_4;}
						if (iIndex == 8) {no_0 = ip_0; no_1 = ip_4;}
						if (iIndex == 9) {no_0 = ip_1; no_1 = ip_5;}
						if (iIndex == 10) {no_0 = ip_2; no_1 = ip_6;}
						if (iIndex == 11) {no_0 = ip_3; no_1 = ip_7;}
						
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.5*(geometry->node[no_0]->GetCoord(iDim)+geometry->node[no_1]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.5 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.5 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
			
/*			for (long iIndex = 0; iIndex < 6; iIndex++) {
				if (HexaFaceCode[iElem][iIndex] == true) {
					
					if (NodeAtFaces[HexaFaceIndex[iElem][iIndex]] != -1)
						HexaFaceNode[iElem][iIndex] = NodeAtFaces[HexaFaceIndex[iElem][iIndex]];
					
					if (NodeAtFaces[HexaFaceIndex[iElem][iIndex]] == -1) {
						NodeAtFaces[HexaFaceIndex[iElem][iIndex]] = nPoint_new;
						HexaFaceNode[iElem][iIndex] = nPoint_new;
						// Compute the coordinates of the new node
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_2; no_3 = ip_3;}
						if (iIndex == 1) {no_0 = ip_4; no_1 = ip_5; no_2 = ip_6; no_3 = ip_7;}
						if (iIndex == 2) {no_0 = ip_1; no_1 = ip_2; no_2 = ip_6; no_3 = ip_5;}
						if (iIndex == 3) {no_0 = ip_3; no_1 = ip_2; no_2 = ip_6; no_3 = ip_7;}
						if (iIndex == 4) {no_0 = ip_0; no_1 = ip_3; no_2 = ip_7; no_3 = ip_4;}
						if (iIndex == 5) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_5; no_3 = ip_4;}
						
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
								NewNodeCoord[nPoint_new][iDim] = 0.25*(x_no[no_0][iDim]+x_no[no_1][iDim]+x_no[no_2][iDim]+x_no[no_3][iDim]);
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.25 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]+ConsVar_Adapt[no_2][iVar]+ConsVar_Adapt[no_3][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.25 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]+AdjVar_Adapt[no_2][iVar]+AdjVar_Adapt[no_3][iVar]);
							if (Restart_Linear) LinVar_Adapt[nPoint_new][iVar] =  0.25 * (LinVar_Adapt[no_0][iVar]+LinVar_Adapt[no_1][iVar]+LinVar_Adapt[no_2][iVar]+LinVar_Adapt[no_3][iVar]);
						}
						
						nPoint_new++;
						
					}
				} 
 
			} */
			
			for (long iIndex = 0; iIndex < 1; iIndex++) {
				
				if (HexaElemCode[iElem][iIndex] == true) {
					
					if (NodeAtElem[HexaElemIndex[iElem][iIndex]] != -1) 
						HexaElemNode[iElem][iIndex] = NodeAtElem[HexaElemIndex[iElem][iIndex]];
					
					if (NodeAtElem[HexaElemIndex[iElem][iIndex]] == -1) {
						NodeAtElem[HexaElemIndex[iElem][iIndex]] = nPoint_new;
						HexaElemNode[iElem][iIndex] = nPoint_new;
						
						// Compute the coordinates of the new node
						if (iIndex == 0) {no_0 = ip_0; no_1 = ip_1; no_2 = ip_2; no_3 = ip_3; no_4 = ip_4; no_5 = ip_5; no_6 = ip_6; no_7 = ip_7;}
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++) 
							NewNodeCoord[nPoint_new][iDim] = 0.125*(geometry->node[no_0]->GetCoord(iDim) +
																											geometry->node[no_1]->GetCoord(iDim) +
																											geometry->node[no_2]->GetCoord(iDim) +
																											geometry->node[no_3]->GetCoord(iDim) +
																											geometry->node[no_4]->GetCoord(iDim) +
																											geometry->node[no_5]->GetCoord(iDim) +
																											geometry->node[no_6]->GetCoord(iDim) +
																											geometry->node[no_7]->GetCoord(iDim));
						
						for (iVar = 0; iVar < nVar; iVar ++) {
							if (Restart_Flow) ConsVar_Adapt[nPoint_new][iVar] =  0.125 * (ConsVar_Adapt[no_0][iVar]+ConsVar_Adapt[no_1][iVar]+ConsVar_Adapt[no_2][iVar]+ConsVar_Adapt[no_3][iVar]+
																																						ConsVar_Adapt[no_4][iVar]+ConsVar_Adapt[no_5][iVar]+ConsVar_Adapt[no_6][iVar]+ConsVar_Adapt[no_7][iVar]);
							if (Restart_Adjoint) AdjVar_Adapt[nPoint_new][iVar] =  0.125 * (AdjVar_Adapt[no_0][iVar]+AdjVar_Adapt[no_1][iVar]+AdjVar_Adapt[no_2][iVar]+AdjVar_Adapt[no_3][iVar]+
																																							AdjVar_Adapt[no_4][iVar]+AdjVar_Adapt[no_5][iVar]+AdjVar_Adapt[no_6][iVar]+AdjVar_Adapt[no_7][iVar]);
						}
						
						nPoint_new++;
					}
				}
			}
		}
		
	}
	
	
	// if Hexa adapt code equals 0, then a semidivision is applied
	long nSemiDivided = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			if (HexaAdaptCode[iElem] == 0)
				nSemiDivided++;
		}
	}
	
	// If semidivision, then divide add a new point (hexa), divide the hexahedron into Tetras, 
	// and find the right combination, it also create the new node (hexa).
	long nPyram = nSemiDivided*6;	
	
	long *PyramAdaptCode;
	long **PyramNode;	
	long **HexaPyramIndex;
	long *PyramHexaIndex;
	
	long **PyramEdgeIndex; 
	bool **PyramEdgeCode;
	long **PyramEdgeNode; 
	
	long **PyramFaceNode;
	
	long **PyramElemNode;
	
	PyramAdaptCode = new long [nPyram];
	PyramNode = new long *[nPyram];	
	HexaPyramIndex = new long *[geometry->GetnElem()];	
	PyramHexaIndex = new long [nPyram];	
	PyramEdgeIndex = new long *[nPyram]; 
	PyramEdgeCode = new bool *[nPyram];
	PyramEdgeNode = new long *[nPyram]; 
	PyramFaceNode = new long *[nPyram];
	PyramElemNode = new long *[nPyram];
	
	for (long iPyram = 0; iPyram < nPyram; iPyram++) {
		PyramNode[iPyram] = new long [4];
		PyramEdgeIndex[iPyram] = new long [4];
		PyramEdgeCode[iPyram] = new bool [4];
		PyramEdgeNode[iPyram] = new long [4];
		PyramFaceNode[iPyram] = new long [1];
		PyramElemNode[iPyram] = new long [1];
	}
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			HexaPyramIndex[iElem] = new long [6];
		}
	}
	
	nPyram = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			if (HexaAdaptCode[iElem] == 0) {
				
				// Write the edge combination on the base.			
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				ip_3 = geometry->elem[iElem]->GetNode(3);
				ip_4 = geometry->elem[iElem]->GetNode(4);
				ip_5 = geometry->elem[iElem]->GetNode(5);
				ip_6 = geometry->elem[iElem]->GetNode(6);
				ip_7 = geometry->elem[iElem]->GetNode(7);
				
				// Create the 1st pyramid.			
				HexaPyramIndex[iElem][0] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_1;
				PyramNode[nPyram][2] = ip_2; PyramNode[nPyram][3] = ip_3;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][0]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][1];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][2]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][3];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][0]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][1];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][2]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][3];
				
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				// Create the 2nd pyramid.			
				HexaPyramIndex[iElem][1] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_4; PyramNode[nPyram][1] = ip_7;
				PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_5;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][7]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][6];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][5]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][4];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][7]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][6];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][5]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][4];
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				// Create the 3th pyramid.			
				HexaPyramIndex[iElem][2] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_4;
				PyramNode[nPyram][2] = ip_5; PyramNode[nPyram][3] = ip_1;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][8]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][4];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][9]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][0];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][8]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][4];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][9]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][0];
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				// Create the 4th pyramid.		
				HexaPyramIndex[iElem][3] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_3; PyramNode[nPyram][1] = ip_2;
				PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_7;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][2]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][10];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][6]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][11];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][2]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][10];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][6]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][11];
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				// Create the 5th pyramid.			
				HexaPyramIndex[iElem][4] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_1; PyramNode[nPyram][1] = ip_5;
				PyramNode[nPyram][2] = ip_6; PyramNode[nPyram][3] = ip_2;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][9]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][5];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][10]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][1];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][9]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][5];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][10]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][1];
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				// Create the 6th pyramid.	
				HexaPyramIndex[iElem][5] = nPyram; PyramHexaIndex[nPyram] = iElem; PyramElemNode[nPyram][0] = nPoint_new; 
				
				PyramNode[nPyram][0] = ip_0; PyramNode[nPyram][1] = ip_3;
				PyramNode[nPyram][2] = ip_7; PyramNode[nPyram][3] = ip_4;
				
				PyramEdgeIndex[nPyram][0] = HexaEdgeIndex[iElem][3]; PyramEdgeIndex[nPyram][1] = HexaEdgeIndex[iElem][11];
				PyramEdgeIndex[nPyram][2] = HexaEdgeIndex[iElem][7]; PyramEdgeIndex[nPyram][3] = HexaEdgeIndex[iElem][8];
				
				PyramEdgeNode[nPyram][0] = HexaEdgeNode[iElem][3]; PyramEdgeNode[nPyram][1] = HexaEdgeNode[iElem][11];
				PyramEdgeNode[nPyram][2] = HexaEdgeNode[iElem][7]; PyramEdgeNode[nPyram][3] = HexaEdgeNode[iElem][8];
				
				for (long iIndex = 0; iIndex < 4; iIndex++)
					if (DivEdge[PyramEdgeIndex[nPyram][iIndex]] == true) PyramEdgeCode[nPyram][iIndex] = true; 
				nPyram++;
				
				nPoint_new++; 

				
			}
		}
	}
	
	//	Check the kind of Pyram partitioning that should be applied
	for (long iPyram = 0; iPyram < nPyram; iPyram ++) {
		PyramAdaptCode[iPyram] = CheckPyramCode(PyramEdgeCode[iPyram]);
		if (PyramAdaptCode[iPyram] == 0) cout << "There is a problem with one Pyram" << endl;
	}
	
	//  Create new structure
	nElem_new = geometry->GetnElem();
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			if (TetraAdaptCode[iElem] == 1) nElem_new = nElem_new + 7;
			if (TetraAdaptCode[iElem] == 2) nElem_new = nElem_new + 3;
			if (TetraAdaptCode[iElem] == 3) nElem_new = nElem_new + 3;
			if (TetraAdaptCode[iElem] == 4) nElem_new = nElem_new + 3;
			if (TetraAdaptCode[iElem] == 5) nElem_new = nElem_new + 3;
			if (TetraAdaptCode[iElem] == 6) nElem_new = nElem_new + 7;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			if (HexaAdaptCode[iElem] == 1) nElem_new = nElem_new + 7;
			if (HexaAdaptCode[iElem] == 2) nElem_new = nElem_new + 3;
			if (HexaAdaptCode[iElem] == 3) nElem_new = nElem_new + 3;
			if (HexaAdaptCode[iElem] == 4) nElem_new = nElem_new + 3;
			if (HexaAdaptCode[iElem] == 5) nElem_new = nElem_new + 1;
			if (HexaAdaptCode[iElem] == 6) nElem_new = nElem_new + 1;
			if (HexaAdaptCode[iElem] == 7) nElem_new = nElem_new + 1;
			if (HexaAdaptCode[iElem] == 0) {
				long iPyram;
				for (long iIndex = 0; iIndex < 6; iIndex++) {
					iPyram = HexaPyramIndex[iElem][0];
					if (PyramAdaptCode[iPyram] == 1) nElem_new = nElem_new + 0;
					if (PyramAdaptCode[iPyram] == 2) nElem_new = nElem_new + 2;
					if (PyramAdaptCode[iPyram] == 3) nElem_new = nElem_new + 2;
					if (PyramAdaptCode[iPyram] == 4) nElem_new = nElem_new + 2;
					if (PyramAdaptCode[iPyram] == 5) nElem_new = nElem_new + 2;
					if (PyramAdaptCode[iPyram] == 6) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 7) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 8) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 9) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 10) nElem_new = nElem_new + 1;
					if (PyramAdaptCode[iPyram] == 11) nElem_new = nElem_new + 1;
					if (PyramAdaptCode[iPyram] == 12) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 13) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 14) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 15) nElem_new = nElem_new + 3;
					if (PyramAdaptCode[iPyram] == 16) nElem_new = nElem_new + 3;
				}
			}
		}
	}
	
	// New points
	geo_adapt->node = new CPoint*[nPoint_new];
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
		geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0), geometry->node[iPoint]->GetCoord(1), geometry->node[iPoint]->GetCoord(2), iPoint, config);
	
	for (iPoint = geometry->GetnPoint(); iPoint < nPoint_new; iPoint++)
		geo_adapt->node[iPoint] = new CPoint(NewNodeCoord[iPoint][0], NewNodeCoord[iPoint][1], NewNodeCoord[iPoint][2], iPoint, config);
	
	// New elements
	geo_adapt->elem = new CPrimalGrid*[nElem_new];
	
	unsigned long iElemNew = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			if (TetraAdaptCode[iElem] == -1) {
				geo_adapt->elem[iElemNew] = new CTetrahedron(geometry->elem[iElem]->GetNode(0), 
																										 geometry->elem[iElem]->GetNode(1), 
																										 geometry->elem[iElem]->GetNode(2), 
																										 geometry->elem[iElem]->GetNode(3));
				iElemNew++;
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			if (HexaAdaptCode[iElem] == -1) {
				geo_adapt->elem[iElemNew] = new CHexahedron(geometry->elem[iElem]->GetNode(0),
																										geometry->elem[iElem]->GetNode(1),
																										geometry->elem[iElem]->GetNode(2),
																										geometry->elem[iElem]->GetNode(3),
																										geometry->elem[iElem]->GetNode(4),
																										geometry->elem[iElem]->GetNode(5),
																										geometry->elem[iElem]->GetNode(6),
																										geometry->elem[iElem]->GetNode(7));
				iElemNew++;
			}
		}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID) {
			geo_adapt->elem[iElemNew] = new CPyramid(geometry->elem[iElem]->GetNode(0),
																							 geometry->elem[iElem]->GetNode(1),
																							 geometry->elem[iElem]->GetNode(2),
																							 geometry->elem[iElem]->GetNode(3),
																							 geometry->elem[iElem]->GetNode(4));
			iElemNew++;
		}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM) {
			geo_adapt->elem[iElemNew] = new CPrism(geometry->elem[iElem]->GetNode(0),
																									geometry->elem[iElem]->GetNode(1),
																									geometry->elem[iElem]->GetNode(2),
																									geometry->elem[iElem]->GetNode(3),
																									geometry->elem[iElem]->GetNode(4),
																									geometry->elem[iElem]->GetNode(5));
			iElemNew++;
		}		
	}
		
	long nodes[27]; long **Division; long nPart;
	
	Division = new long*[100];
	for (long iVar = 0; iVar < 100; iVar++)
		Division[iVar] = new long[100];
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			
			// Tetra elements...
			if (TetraAdaptCode[iElem] > 0) {
				
				// First the corners
				nodes[0] = geometry->elem[iElem]->GetNode(0);		
				nodes[1] = geometry->elem[iElem]->GetNode(1); 
				nodes[2] = geometry->elem[iElem]->GetNode(2);
				nodes[3] = geometry->elem[iElem]->GetNode(3);
				
				// Next the points that correspond to the broken edges.
				nodes[4] = TetraEdgeNode[iElem][0]; 
				nodes[5] = TetraEdgeNode[iElem][1];
				nodes[6] = TetraEdgeNode[iElem][2];
				nodes[7] = TetraEdgeNode[iElem][3];
				nodes[8] = TetraEdgeNode[iElem][4];
				nodes[9] = TetraEdgeNode[iElem][5];
				
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				ip_3 = geometry->elem[iElem]->GetNode(3);
				
				long edges[6] = {-1, -1, -1, -1, -1 , -1};
				if (DivEdge[geometry->FindEdge(ip_0, ip_1)]) {edges[0] = geometry->FindEdge(ip_0, ip_1);}
				if (DivEdge[geometry->FindEdge(ip_0, ip_2)]) {edges[1] = geometry->FindEdge(ip_0, ip_2);}
				if (DivEdge[geometry->FindEdge(ip_0, ip_3)]) {edges[2] = geometry->FindEdge(ip_0, ip_3);}
				if (DivEdge[geometry->FindEdge(ip_1, ip_2)]) {edges[3] = geometry->FindEdge(ip_1, ip_2);}
				if (DivEdge[geometry->FindEdge(ip_1, ip_3)]) {edges[4] = geometry->FindEdge(ip_1, ip_3);}
				if (DivEdge[geometry->FindEdge(ip_2, ip_3)]) {edges[5] = geometry->FindEdge(ip_2, ip_3);}		

				
				TetraDivision(TetraAdaptCode[iElem], nodes, edges, Division, &nPart);
				
				for (long iPart = 0; iPart < nPart; iPart++) {
					
					// Tetra case
					geo_adapt->elem[iElemNew] = new CTetrahedron(Division[iPart][1], 
																											 Division[iPart][2], 
																											 Division[iPart][3], 
																											 Division[iPart][4]);
					iElemNew++;
				}
			}
		}
		
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
			// Hexa elements...
			if (HexaAdaptCode[iElem] > 0) {
				
				// First the corners
				nodes[0] = geometry->elem[iElem]->GetNode(0);		nodes[1] = geometry->elem[iElem]->GetNode(1);
				nodes[2] = geometry->elem[iElem]->GetNode(2);		nodes[3] = geometry->elem[iElem]->GetNode(3);
				nodes[4] = geometry->elem[iElem]->GetNode(4);		nodes[5] = geometry->elem[iElem]->GetNode(5);
				nodes[6] = geometry->elem[iElem]->GetNode(6);		nodes[7] = geometry->elem[iElem]->GetNode(7);
				
				// Next the points that correspond to the broken edges.
				nodes[8] = HexaEdgeNode[iElem][0]; nodes[9] = HexaEdgeNode[iElem][1];
				nodes[10] = HexaEdgeNode[iElem][2]; nodes[11] = HexaEdgeNode[iElem][3];
				nodes[12] = HexaEdgeNode[iElem][4]; nodes[13] = HexaEdgeNode[iElem][5];
				nodes[14] = HexaEdgeNode[iElem][6]; nodes[15] = HexaEdgeNode[iElem][7];
				nodes[16] = HexaEdgeNode[iElem][8]; nodes[17] = HexaEdgeNode[iElem][9];
				nodes[18] = HexaEdgeNode[iElem][10]; nodes[19] = HexaEdgeNode[iElem][11];
				
				// Next the points that correspond to the faces.
				nodes[20] = HexaFaceNode[iElem][0];
				nodes[21] = HexaFaceNode[iElem][1];
				nodes[22] = HexaFaceNode[iElem][2];
				nodes[23] = HexaFaceNode[iElem][3];
				nodes[24] = HexaFaceNode[iElem][4];
				nodes[25] = HexaFaceNode[iElem][5];
				
				// Next the points that correspond to the element.
				nodes[8] = HexaElemNode[iElem][0];
				
				HexaDivision(HexaAdaptCode[iElem], nodes, Division, &nPart);
				for (long iPart = 0; iPart < nPart; iPart++) {
					geo_adapt->elem[iElemNew] = new CHexahedron(Division[iPart][1], 
																											Division[iPart][2], 
																											Division[iPart][3], 
																											Division[iPart][4],
																											Division[iPart][5], 
																											Division[iPart][6], 
																											Division[iPart][7], 
																											Division[iPart][8]);
					iElemNew++;
				}
			}
			
			// Pyram elements...
			if (HexaAdaptCode[iElem] == 0) {
				long iPyram = HexaPyramIndex[iElem][0];
				
				// First the corners
				nodes[0] = PyramNode[iPyram][0];	
				nodes[1] = PyramNode[iPyram][1];	
				nodes[2] = PyramNode[iPyram][2];	
				nodes[3] = PyramNode[iPyram][3];	
				nodes[4] = PyramElemNode[iPyram][0];

				// Next the points that correspond to the broken edges.
				nodes[5] = PyramEdgeNode[iPyram][0];
				nodes[6] = PyramEdgeNode[iPyram][1];
				nodes[7] = PyramEdgeNode[iPyram][2];
				nodes[8] = PyramEdgeNode[iPyram][3];
				
				// Next the points that correspond to the base face.
				nodes[9] = PyramFaceNode[iPyram][0];
				
				PyramDivision(PyramAdaptCode[iPyram], nodes, Division, &nPart);
				for (long iPart = 0; iPart < nPart; iPart++) {
					
					// Tetra case
					if (Division[iPart][0] == 5) {
						geo_adapt->elem[iElemNew] = new CTetrahedron(Division[iPart][1], 
																												 Division[iPart][2], 
																												 Division[iPart][3], 
																												 Division[iPart][4]);
						iElemNew++;
					}
					
					// Pyram case
					if (Division[iPart][0] == 6) {
						geo_adapt->elem[iElemNew] = new CPyramid(Division[iPart][1], 
																												Division[iPart][2], 
																												Division[iPart][3], 
																												Division[iPart][4],
																												Division[iPart][5]);
						iElemNew++;
					}
				}
			}
		}
	}
	
	geo_adapt->SetnElem(iElemNew);
	geo_adapt->SetnPoint(nPoint_new);
	geo_adapt->SetnPointDomain(nPoint_new);
	geo_adapt->SetnDim(nDim);
	
	//  Create boundary structure
	geo_adapt->SetnMarker(geometry->GetnMarker());
	geo_adapt->nElem_Bound = new unsigned long [geometry->GetnMarker()];
	geo_adapt->Tag_to_Marker = new string [nMarker_Max];
	geo_adapt->bound = new CPrimalGrid**[geometry->GetnMarker()];

	// Conservative estimation of the number of boundary elements.
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		long nNewBCcv = 0;
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {			
			nNewBCcv = nNewBCcv + 4;
		}
		geo_adapt->bound[iMarker] = new CPrimalGrid* [nNewBCcv];
		geo_adapt->SetMarker_Tag(iMarker, geometry->GetMarker_Tag(iMarker));
	}

	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		long nNewBCcv = 0;
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {			
			
			bool TriangleEdgeCode[3] = { false, false, false };
			long TriangleAdaptCode, nodes[6];
			long nNodesBound = geometry->bound[iMarker][iVertex]->GetnNodes();
			
			ip_0 = geometry->bound[iMarker][iVertex]->GetNode(0); geo_adapt->node[ip_0]->SetBoundary(geometry->GetnMarker());
			ip_1 = geometry->bound[iMarker][iVertex]->GetNode(1); geo_adapt->node[ip_1]->SetBoundary(geometry->GetnMarker());
			ip_2 = geometry->bound[iMarker][iVertex]->GetNode(2); geo_adapt->node[ip_2]->SetBoundary(geometry->GetnMarker());
			if (nNodesBound == 4) {
				ip_3 = geometry->bound[iMarker][iVertex]->GetNode(3); geo_adapt->node[ip_3]->SetBoundary(geometry->GetnMarker());
				geo_adapt->bound[iMarker][nNewBCcv] = new CQuadrilateral(ip_0, ip_1, ip_2, ip_3, 3);
				nNewBCcv++;
			}
			else {
				// First the corners
				nodes[0] = geometry->bound[iMarker][iVertex]->GetNode(0);	
				nodes[1] = geometry->bound[iMarker][iVertex]->GetNode(1); 
				nodes[2] = geometry->bound[iMarker][iVertex]->GetNode(2);
				if (nNodesBound == 4) nodes[3] = geometry->bound[iMarker][iVertex]->GetNode(3);
				
				// Next the points that correspond to the broken edges.
				nodes[3] = NodeAtEdges[geometry->FindEdge(ip_0, ip_1)]; if (nodes[3] != -1) geo_adapt->node[nodes[3]]->SetBoundary(geometry->GetnMarker());
				nodes[4] = NodeAtEdges[geometry->FindEdge(ip_1, ip_2)]; if (nodes[4] != -1) geo_adapt->node[nodes[4]]->SetBoundary(geometry->GetnMarker());
				nodes[5] = NodeAtEdges[geometry->FindEdge(ip_0, ip_2)]; if (nodes[5] != -1) geo_adapt->node[nodes[5]]->SetBoundary(geometry->GetnMarker());
				
				long edges[3] = {-1, -1, -1};
				if (DivEdge[geometry->FindEdge(ip_0, ip_1)]) { edges[0] = geometry->FindEdge(ip_0, ip_1); TriangleEdgeCode[0] = true;}
				if (DivEdge[geometry->FindEdge(ip_1, ip_2)]) { edges[1] = geometry->FindEdge(ip_1, ip_2); TriangleEdgeCode[1] = true;}
				if (DivEdge[geometry->FindEdge(ip_2, ip_0)]) { edges[2] = geometry->FindEdge(ip_2, ip_0); TriangleEdgeCode[2] = true;}
				
				TriangleAdaptCode = CheckTriangleCode(TriangleEdgeCode);
				TriangleDivision(TriangleAdaptCode, nodes, edges, Division, &nPart);
				
				for (long iPart = 0; iPart < nPart; iPart++) {
					geo_adapt->bound[iMarker][nNewBCcv] = new CTriangle(Division[iPart][1], Division[iPart][2], Division[iPart][3], 3);
					nNewBCcv++;
				}
			}
		}
		geo_adapt->SetnElem_Bound(iMarker, nNewBCcv);
	}
	
	delete [] DivEdge; 
	delete [] DivElem;
	delete [] NodeAtEdges; 
	delete [] NodeAtElem;
}

void CGridAdaptation::SetIndicator_Flow(CGeometry *geometry, CConfig *config, unsigned short strength) {
  unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, max_elem_new;
  unsigned short iDim, iMarker;
  su2double Dual_Area, norm, Solution_Vertex, Solution_0, Solution_1, Solution_Average,
  DualArea, Partial_Res, Grad_Val, *Normal;
  su2double scale_area = config->GetDualVol_Power();
  
  /*--- Initialization ---*/
  nElem_new = 0;
  max_elem_new = SU2_TYPE::Int(0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
  for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
    geometry->elem[iElem]->SetDivide(false);
  }
  
  /*--- Compute the gradient of the first variable ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient[iPoint][iDim] = 0.0;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    Point_0 = geometry->edge[iEdge]->GetNode(0); Solution_0 = ConsVar_Sol[Point_0][0];
    Point_1 = geometry->edge[iEdge]->GetNode(1); Solution_1 = ConsVar_Sol[Point_1][0];
    Normal = geometry->edge[iEdge]->GetNormal();
    Solution_Average =  0.5 * ( Solution_0 + Solution_1);
    for (iDim = 0; iDim < nDim; iDim++) {
      Partial_Res = Solution_Average*Normal[iDim];
      Gradient[Point_0][iDim] = Gradient[Point_0][iDim] + Partial_Res;
      Gradient[Point_1][iDim] = Gradient[Point_1][iDim] - Partial_Res;
    }
  }
		
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        Point = geometry->vertex[iMarker][iVertex]->GetNode();
        Solution_Vertex = ConsVar_Sol[Point][0];
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++) {
          Partial_Res = Solution_Vertex*Normal[iDim];
          Gradient[Point][iDim] = Gradient[Point][iDim] - Partial_Res;
        }
      }
    }
  
  for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      DualArea = geometry->node[iPoint]->GetVolume();
      Grad_Val = Gradient[iPoint][iDim]/DualArea;
      Gradient[iPoint][iDim] = Grad_Val;
    }
  
  /*--- Compute the the adaptation index at each point ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    Dual_Area = geometry->node[iPoint]->GetVolume();
    norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm += Gradient[iPoint][iDim]*Gradient[iPoint][iDim];
    norm = sqrt(norm);
    Index[iPoint] = pow(Dual_Area, scale_area)*norm;
  }
  
  SetSensorElem(geometry, config, max_elem_new);
  
}


void CGridAdaptation::SetIndicator_Adj(CGeometry *geometry, CConfig *config, unsigned short strength) {
  su2double Dual_Area;
  unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, max_elem_new;
  unsigned short iDim, iMarker;
  su2double norm, Solution_Vertex, Solution_0, Solution_1, Solution_Average,
  DualArea, Partial_Res, Grad_Val, *Normal;
  su2double scale_area = config->GetDualVol_Power();
  
  // Initialization
  nElem_new = 0;
  max_elem_new = SU2_TYPE::Int(0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
  for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
    geometry->elem[iElem]->SetDivide(false);
  }
  
  
  // Compute the gradient of the density.
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient[iPoint][iDim] = 0.0;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    Point_0 = geometry->edge[iEdge]->GetNode(0); Solution_0 = AdjVar_Sol[Point_0][0];
    Point_1 = geometry->edge[iEdge]->GetNode(1); Solution_1 = AdjVar_Sol[Point_1][0];
    Normal = geometry->edge[iEdge]->GetNormal();
    Solution_Average =  0.5 * ( Solution_0 + Solution_1);
    for (iDim = 0; iDim < nDim; iDim++) {
      Partial_Res = Solution_Average*Normal[iDim];
      Gradient[Point_0][iDim] = Gradient[Point_0][iDim] + Partial_Res;
      Gradient[Point_1][iDim] = Gradient[Point_1][iDim] - Partial_Res;
    }
  }
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        Point = geometry->vertex[iMarker][iVertex]->GetNode();
        Solution_Vertex = AdjVar_Sol[Point][0];
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++) {
          Partial_Res = Solution_Vertex*Normal[iDim];
          Gradient[Point][iDim] = Gradient[Point][iDim] - Partial_Res;
        }
      }
    }
  
  for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      DualArea = geometry->node[iPoint]->GetVolume();
      Grad_Val = Gradient[iPoint][iDim]/DualArea;
      Gradient[iPoint][iDim] = Grad_Val;
    }
  
  
  // Compute the the adaptation index at each point.
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    Dual_Area = geometry->node[iPoint]->GetVolume();
    norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm += Gradient[iPoint][iDim]*Gradient[iPoint][iDim];
    norm = sqrt(norm);
    Index[iPoint] = pow(Dual_Area, scale_area)*norm;
  }
  
  SetSensorElem(geometry, config, max_elem_new);
  
}

void CGridAdaptation::SetIndicator_FlowAdj(CGeometry *geometry, CConfig *config) {
  su2double Dual_Area;
  unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, max_elem_new_flow, max_elem_new_adj;
  unsigned short iDim, iMarker;
  su2double norm, DualArea, Partial_Res, *Normal;
  su2double scale_area = config->GetDualVol_Power();
  
  // Initialization
  max_elem_new_flow = SU2_TYPE::Int(0.5*0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
  max_elem_new_adj =  SU2_TYPE::Int(0.5*0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
  for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
    geometry->elem[iElem]->SetDivide(false);
  }
  
  // Compute the gradient of the first variable.
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      Gradient_Flow[iPoint][iDim] = 0.0;
      Gradient_Adj[iPoint][iDim] = 0.0;
    }
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    Point_0 = geometry->edge[iEdge]->GetNode(0);
    Point_1 = geometry->edge[iEdge]->GetNode(1);
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iDim = 0; iDim < nDim; iDim++) {
      Partial_Res = 0.5 * ( ConsVar_Sol[Point_0][0] + ConsVar_Sol[Point_1][0] ) * Normal[iDim];
      Gradient_Flow[Point_0][iDim] = Gradient_Flow[Point_0][iDim] + Partial_Res;
      Gradient_Flow[Point_1][iDim] = Gradient_Flow[Point_1][iDim] - Partial_Res;
      
      Partial_Res = 0.5 * ( AdjVar_Sol[Point_0][0] + AdjVar_Sol[Point_1][0] ) * Normal[iDim];
      Gradient_Adj[Point_0][iDim] = Gradient_Adj[Point_0][iDim] + Partial_Res;
      Gradient_Adj[Point_1][iDim] = Gradient_Adj[Point_1][iDim] - Partial_Res;
    }
  }
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        Point = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++) {
          Gradient_Flow[Point][iDim] = Gradient_Flow[Point][iDim] - ConsVar_Sol[Point][0] * Normal[iDim];
          Gradient_Adj[Point][iDim] = Gradient_Adj[Point][iDim] - AdjVar_Sol[Point][0] * Normal[iDim];
        }
      }
    }
  
  for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      DualArea = geometry->node[iPoint]->GetVolume();
      Gradient_Flow[iPoint][iDim] = Gradient_Flow[iPoint][iDim]/DualArea;
      Gradient_Adj[iPoint][iDim] = Gradient_Adj[iPoint][iDim]/DualArea;
    }
  
  // Compute the the adaptation index at each point.
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    Dual_Area=geometry->node[iPoint]->GetVolume();
    norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm += Gradient_Flow[iPoint][iDim]*Gradient_Flow[iPoint][iDim];
    norm = sqrt(norm);
    Index[iPoint] = pow(Dual_Area, scale_area)*norm;
  }
  
  
  SetSensorElem(geometry, config, max_elem_new_flow);
  
  // Compute the the adaptation index at each point.
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    Dual_Area=geometry->node[iPoint]->GetVolume();
    norm = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm += Gradient_Adj[iPoint][iDim]*Gradient_Adj[iPoint][iDim];
    norm = sqrt(norm); 
    Index[iPoint] = pow(Dual_Area, scale_area)*norm;
  }
  
  SetSensorElem(geometry, config, max_elem_new_adj);
  
}

void CGridAdaptation::SetIndicator_Robust(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iElem, max_elem_new_flow, max_elem_new_adj;
	unsigned short iVar;
	su2double Dual_Area;
	su2double scale_area = config->GetDualVol_Power();
	
	// Inicializa la malla para la adaptacion
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
	}
	
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += ConsVar_Res[iPoint][iVar]*ConsVar_Res[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area, scale_area)*sqrt(Index[iPoint]);
	}
	
	max_elem_new_flow = SU2_TYPE::Int(0.5*0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
	SetSensorElem(geometry, config, max_elem_new_flow);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += AdjVar_Res[iPoint][iVar]*AdjVar_Res[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area, scale_area)*sqrt(Index[iPoint]);
	}
	
	max_elem_new_adj = SU2_TYPE::Int(0.5*0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));
	SetSensorElem(geometry, config, max_elem_new_adj);

}

void CGridAdaptation::SetIndicator_Computable(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iElem, max_elem_new;
	unsigned short iVar;
	su2double Dual_Area;
	su2double scale_area = config->GetDualVol_Power();
	
	max_elem_new = SU2_TYPE::Int(0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
	}
		
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += ConsVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar]*ConsVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area, scale_area)*sqrt(Index[iPoint]);
	}
	
	SetSensorElem(geometry, config, max_elem_new);

}

void CGridAdaptation::SetIndicator_Computable_Robust(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, iElem, max_elem_new;
	unsigned short iVar;
	su2double Dual_Area ;
	su2double scale_area = config->GetDualVol_Power();

	/*--- Initializate the numerical grid for the adaptation ---*/
	max_elem_new = SU2_TYPE::Int(0.01*config->GetNew_Elem_Adapt()*su2double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
	}
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += LinVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar]*LinVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area, scale_area)*sqrt(Index[iPoint]);
	}

	SetSensorElem(geometry, config, max_elem_new);
	
}

void CGridAdaptation::SetRestart_FlowSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_flowfilename) {
	
	unsigned long iPoint;
	unsigned short iVar, iDim;
		
	char *cstr = new char [mesh_flowfilename.size()+1];
	strcpy (cstr, mesh_flowfilename.c_str());
	
	ofstream restart_flowfile;
	restart_flowfile.open(cstr, ios::out);
	restart_flowfile.precision(15);
  
  restart_flowfile << "Restart file generated with SU2_MSH" << endl;

	for (iPoint = 0; iPoint < nPoint_new; iPoint++) {
		restart_flowfile << iPoint <<"\t";

    for (iDim = 0; iDim < nDim; iDim++)
			restart_flowfile << scientific << geo_adapt->node[iPoint]->GetCoord(iDim) <<"\t";
		for (iVar = 0; iVar < nVar; iVar++)
			restart_flowfile << scientific << ConsVar_Adapt[iPoint][iVar] <<"\t";

		restart_flowfile << endl;
	}
	restart_flowfile.close();
	
}

void CGridAdaptation::SetRestart_AdjSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_adjfilename) {
	
  char cstr[MAX_STRING_SIZE], buffer[50];
  unsigned short iDim, iVar;
  unsigned long iPoint;
	string copy;
	
	copy.assign(mesh_adjfilename);
  unsigned short lastindex = copy.find_last_of(".");
  copy = copy.substr(0, lastindex);
	strcpy (cstr, copy.c_str());
  if (config->GetnObj() > 1) {
    SPRINTF (buffer, "_combo.dat");
  }
  else {
    if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT)        SPRINTF (buffer, "_cd.dat");
    if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT)        SPRINTF (buffer, "_cl.dat");
    if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT)   SPRINTF (buffer, "_csf.dat");
    if (config->GetKind_ObjFunc() == INVERSE_DESIGN_PRESSURE) SPRINTF (buffer, "_invpress.dat");
    if (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX) SPRINTF (buffer, "_invheat.dat");
    if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT)    SPRINTF (buffer, "_cmx.dat");
    if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT)    SPRINTF (buffer, "_cmy.dat");
    if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT)    SPRINTF (buffer, "_cmz.dat");
    if (config->GetKind_ObjFunc() == EFFICIENCY)              SPRINTF (buffer, "_eff.dat");
    if (config->GetKind_ObjFunc() == FORCE_X_COEFFICIENT)     SPRINTF (buffer, "_cfx.dat");
    if (config->GetKind_ObjFunc() == FORCE_Y_COEFFICIENT)     SPRINTF (buffer, "_cfy.dat");
    if (config->GetKind_ObjFunc() == FORCE_Z_COEFFICIENT)     SPRINTF (buffer, "_cfz.dat");
    if (config->GetKind_ObjFunc() == TOTAL_HEATFLUX)          SPRINTF (buffer, "_totheat.dat");
    if (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX)        SPRINTF (buffer, "_maxheat.dat");
    if (config->GetKind_ObjFunc() == SURFACE_TOTAL_PRESSURE)  SPRINTF (buffer, "_pt.dat");
    if (config->GetKind_ObjFunc() == SURFACE_STATIC_PRESSURE) SPRINTF (buffer, "_pe.dat");
    if (config->GetKind_ObjFunc() == SURFACE_MASSFLOW)        SPRINTF (buffer, "_mfr.dat");
    if (config->GetKind_ObjFunc() == SURFACE_MACH)            SPRINTF (buffer, "_mach.dat");
    if (config->GetKind_ObjFunc() == CUSTOM_OBJFUNC)          SPRINTF (buffer, "_custom.dat");
  }
  
	strcat(cstr, buffer);
	
	ofstream restart_adjfile;
	restart_adjfile.open(cstr, ios::out);
	restart_adjfile.precision(15);
  
  restart_adjfile << "Restart file generated with SU2_MSH" << endl;
  
	for (iPoint = 0; iPoint < nPoint_new; iPoint++) {
		restart_adjfile << iPoint <<"\t";
    
    for (iDim = 0; iDim < nDim; iDim++)
			restart_adjfile << scientific << geo_adapt->node[iPoint]->GetCoord(iDim) <<"\t";
		for (iVar = 0; iVar < nVar; iVar++)
			restart_adjfile << scientific << AdjVar_Adapt[iPoint][iVar]<<"\t";
		restart_adjfile << endl;
    
	}
	restart_adjfile.close();
}

void CGridAdaptation::SetRestart_LinSolution(CConfig *config, CPhysicalGeometry *geo_adapt, string mesh_linfilename) {
	
	unsigned long iPoint;
	unsigned short iVar, iDim;
	
	char *cstr_ = new char [mesh_linfilename.size()+1];
	strcpy (cstr_, mesh_linfilename.c_str());
	
	ofstream restart_linfile;
	restart_linfile.open(cstr_, ios::out);
	restart_linfile.precision(15);

  restart_linfile << "Restart file generated with SU2_MSH" << endl;

	for (iPoint = 0; iPoint < nPoint_new; iPoint++) {
		restart_linfile << iPoint <<"\t";
    
    for (iDim = 0; iDim < nDim; iDim++)
			restart_linfile << scientific << geo_adapt->node[iPoint]->GetCoord(iDim) <<"\t";
		for (iVar = 0; iVar < nVar; iVar++)
			restart_linfile << scientific << LinVar_Adapt[iPoint][iVar]<<"\t";
		restart_linfile << endl;
    
	}
	restart_linfile.close();
}

void CGridAdaptation::SetSensorElem(CGeometry *geometry, CConfig *config, unsigned long max_elem) {
	su2double Max_Sensor, threshold;
	su2double *Sensor = new su2double[geometry->GetnElem()];
	unsigned long ip_0, ip_1, ip_2, ip_3, iElem, nElem_real;
  
	if (max_elem > geometry->GetnElem()) {
		cout << "WARNING: Attempted to adapt " << max_elem << " cells," << endl;
		cout << "  which is greater than the total number of cells, ";
		cout << geometry->GetnElem() << "." << endl;
		cout << "  Did you set the option NEW_ELEMS in the *.cfg file?" << endl;
  }
	
	/*--- Compute the the adaptation index at each element ---*/
	Max_Sensor = 0.0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		ip_0 = geometry->elem[iElem]->GetNode(0);
		ip_1 = geometry->elem[iElem]->GetNode(1);
		ip_2 = geometry->elem[iElem]->GetNode(2);
		Sensor[iElem] = (Index[ip_0]+Index[ip_1]+Index[ip_2])/3.0;
		if ((geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) ||
			(geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)) {
			ip_3 = geometry->elem[iElem]->GetNode(2);
			Sensor[iElem] = (Index[ip_0]+Index[ip_1]+Index[ip_2]+Index[ip_3])/4.0;
		}
		Max_Sensor = max(Max_Sensor, Sensor[iElem]);
	}
	
	/*--- Adimensionalization of the adaptation sensor ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		Sensor[iElem] = Sensor[iElem]/Max_Sensor;
	}
	
	/*--- Selection of the elements to be adapted ---*/
	threshold = 0.999;
	nElem_real = 0;
	while (nElem_real <= max_elem && threshold >= 0) {
		for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
			if ( Sensor[iElem] >= threshold && !geometry->elem[iElem]->GetDivide() ) {
				if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) nElem_real = nElem_real + 3;
				if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) nElem_real = nElem_real + 3;
				if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem_real = nElem_real + 7;
				geometry->elem[iElem]->SetDivide(true);
				if (nElem_real >= max_elem) break;
			}	
		threshold = threshold - 0.001;
	}

	if (threshold < 0) {
		cout << "WARNING: Tried to find " << max_elem;
		cout << " cells suitable for adaptation, but only found ";
		cout << nElem_real << endl;
		cout << "The following cell types are currently adaptable: " << endl;
		cout << "  + triangles" << endl;
		cout << "  + quadrilaterals" << endl;
		cout << "  + tetrahedrons" << endl;
		cout << "Your grid may have too high a percentage of other types." << endl;
	}
	
	cout << "Number of elements to adapt: " << nElem_real << endl;
	delete [] Sensor;
}
