/*!
 * \file grid_adaptation_structure.cpp
 * \brief Main subroutines for grid adaptation.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.1.
 */

#include "../include/grid_adaptation_structure.hpp"
#include <math.h>

CGridAdaptation::CGridAdaptation(CGeometry *geometry, CConfig *config) {

	unsigned long iPoint;
	
	nDim = geometry->GetnDim();
	if (config->GetKind_Solver() == ELECTRIC_POTENTIAL) nVar = 1;
	else nVar = geometry->GetnDim()+2;
	ConsVar_Sol = new double* [geometry->GetnPoint()];
	AdjVar_Sol = new double* [geometry->GetnPoint()];
	LinVar_Sol = new double* [geometry->GetnPoint()];
	ConsVar_Res = new double* [geometry->GetnPoint()];
	AdjVar_Res = new double* [geometry->GetnPoint()];	
	LinVar_Res = new double* [geometry->GetnPoint()];
	Gradient = new double* [geometry->GetnPoint()];
	Gradient_Flow = new double* [geometry->GetnPoint()];
	Gradient_Adj = new double* [geometry->GetnPoint()];

	Index = new double [geometry->GetnPoint()];

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		ConsVar_Sol[iPoint] = new double [nVar];
		AdjVar_Sol[iPoint] = new double [nVar];
		LinVar_Sol[iPoint] = new double [nVar];
		ConsVar_Res[iPoint] = new double [nVar];
		LinVar_Res[iPoint] = new double [nVar];
		AdjVar_Res[iPoint] = new double [nVar];		
		Gradient[iPoint] = new double [nDim];
		Gradient_Flow[iPoint] = new double [nDim];		
		Gradient_Adj[iPoint] = new double [nDim];				
	}

}

CGridAdaptation::~CGridAdaptation(void) {
	
	unsigned short iVar, iDim;
	
	for (iVar = 0; iVar < nVar; iVar++){
		delete [] ConsVar_Adapt[iVar], ConsVar_Sol[iVar], ConsVar_Res[iVar];
		delete [] AdjVar_Adapt[iVar], AdjVar_Sol[iVar], AdjVar_Res[iVar];
		delete [] LinVar_Adapt[iVar], LinVar_Sol[iVar], LinVar_Res[iVar];
	}
	
	for (iDim = 0; iDim < nDim; iDim++){
		delete [] Gradient[iDim];
		delete [] Gradient_Flow[iDim];
		delete [] Gradient_Adj[iDim];		
	}
	delete [] ConsVar_Adapt, ConsVar_Sol, ConsVar_Res;
	delete [] AdjVar_Adapt, AdjVar_Sol, AdjVar_Res;
	delete [] LinVar_Adapt, LinVar_Sol, LinVar_Res;
	delete [] Gradient;
	delete [] Gradient_Flow;	
	delete [] Gradient_Adj;	
	delete [] Index;	
}

void CGridAdaptation::GetFlowSolution(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	string text_line;
		
	string mesh_filename = config->GetSolution_FlowFileName();
	ifstream restart_file;

	char *cstr = new char [mesh_filename.size()+1];
	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
		cout << "There is no flow restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> ConsVar_Sol[iPoint][0];
		if (nVar == 4) point_line >> index >> ConsVar_Sol[iPoint][0] >> ConsVar_Sol[iPoint][1] >> ConsVar_Sol[iPoint][2] >> ConsVar_Sol[iPoint][3];
		if (nVar == 5) point_line >> index >> ConsVar_Sol[iPoint][0] >> ConsVar_Sol[iPoint][1] >> ConsVar_Sol[iPoint][2] >> ConsVar_Sol[iPoint][3] >> ConsVar_Sol[iPoint][4];
	}
	restart_file.close();
}

void CGridAdaptation::GetFlowResidual(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	double dummy[5];
	string text_line;
	
	string mesh_filename = config->GetSolution_FlowFileName();
	ifstream restart_file;
	
	char *cstr = new char [mesh_filename.size()+1];
	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
		cout << "There is no flow restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> dummy[0] >> ConsVar_Res[iPoint][0];
		if (nVar == 4) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> 
									ConsVar_Res[iPoint][0] >> ConsVar_Res[iPoint][1] >> ConsVar_Res[iPoint][2] >> 
									ConsVar_Res[iPoint][3];
		if (nVar == 5) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4] >> 
									ConsVar_Res[iPoint][0] >> ConsVar_Res[iPoint][1] >> ConsVar_Res[iPoint][2] >> 
									ConsVar_Res[iPoint][3] >> ConsVar_Res[iPoint][4];
	}
	restart_file.close();
}

void CGridAdaptation::GetLinResidual(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	double dummy[5];
	string text_line;
	
	string mesh_filename = config->GetSolution_LinFileName();
	ifstream restart_file;
	
	char *cstr = new char [mesh_filename.size()+1];
	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
		cout << "There is no linear restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> dummy[0] >> LinVar_Res[iPoint][0];
		if (nVar == 4) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> 
			LinVar_Res[iPoint][0] >> LinVar_Res[iPoint][1] >> LinVar_Res[iPoint][2] >> 
			LinVar_Res[iPoint][3];
		if (nVar == 5) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4] >> 
			LinVar_Res[iPoint][0] >> LinVar_Res[iPoint][1] >> LinVar_Res[iPoint][2] >> 
			LinVar_Res[iPoint][3] >> LinVar_Res[iPoint][4];
	}
	
	restart_file.close();
}

void CGridAdaptation::GetAdjSolution(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	string text_line;
	
	string copy, mesh_filename;
	ifstream restart_file;

	char buffer[50], cstr[200];
	mesh_filename = config->GetSolution_AdjFileName();
	copy.assign(mesh_filename);
	copy.erase (copy.end()-4, copy.end());
	strcpy (cstr, copy.c_str()); 
	if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) sprintf (buffer, "_cd.dat"); 
	if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT) sprintf (buffer, "_cl.dat");
	if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT) sprintf (buffer, "_csf.dat"); 
	if (config->GetKind_ObjFunc() == PRESSURE_COEFFICIENT) sprintf (buffer, "_cp.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) sprintf (buffer, "_cmx.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) sprintf (buffer, "_cmy.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) sprintf (buffer, "_cmz.dat"); 
	if (config->GetKind_ObjFunc() == EFFICIENCY) sprintf (buffer, "_eff.dat"); 
	if (config->GetKind_ObjFunc() == ELECTRIC_CHARGE) sprintf (buffer, "_cc.dat"); 
  if (config->GetKind_ObjFunc() == FORCE_X_COEFFICIENT) sprintf (buffer, "_cfx.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Y_COEFFICIENT) sprintf (buffer, "_cfy.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Z_COEFFICIENT) sprintf (buffer, "_cfz.dat"); 
  if (config->GetKind_ObjFunc() == THRUST_COEFFICIENT) sprintf (buffer, "_ct.dat"); 
	if (config->GetKind_ObjFunc() == TORQUE_COEFFICIENT) sprintf (buffer, "_cq.dat"); 
	if (config->GetKind_ObjFunc() == FIGURE_OF_MERIT) sprintf (buffer, "_merit.dat"); 
	strcat(cstr, buffer);
	
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
		cout << "There is no adjoint restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> AdjVar_Sol[iPoint][0];
		if (nVar == 4) point_line >> index >> AdjVar_Sol[iPoint][0] >> AdjVar_Sol[iPoint][1] >> AdjVar_Sol[iPoint][2] >> AdjVar_Sol[iPoint][3];
		if (nVar == 5) point_line >> index >> AdjVar_Sol[iPoint][0] >> AdjVar_Sol[iPoint][1] >> AdjVar_Sol[iPoint][2] >> AdjVar_Sol[iPoint][3] >> AdjVar_Sol[iPoint][4];
	}
	
	restart_file.close();
}

void CGridAdaptation::GetLinSolution(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, index;
	string text_line;
	
	string mesh_filename;
	mesh_filename = config->GetSolution_LinFileName();
	
	ifstream restart_file;
	
	char *cstr = new char [mesh_filename.size()+1];
	strcpy (cstr, mesh_filename.c_str());
	restart_file.open(cstr, ios::in);
	if (restart_file.fail()) {
		cout << "There is no linear restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> LinVar_Sol[iPoint][0];
		if (nVar == 4) point_line >> index >> LinVar_Sol[iPoint][0] >> LinVar_Sol[iPoint][1] >> LinVar_Sol[iPoint][2] >> LinVar_Sol[iPoint][3];
		if (nVar == 5) point_line >> index >> LinVar_Sol[iPoint][0] >> LinVar_Sol[iPoint][1] >> LinVar_Sol[iPoint][2] >> LinVar_Sol[iPoint][3] >> LinVar_Sol[iPoint][4];
	}
	
	restart_file.close();
}

void CGridAdaptation::GetAdjResidual(CGeometry *geometry, CConfig *config){
	unsigned long iPoint, index;
	string text_line;
	double dummy[5];

	string mesh_filename, copy;
	ifstream restart_file;
	
	char buffer[50], cstr[200];
	mesh_filename = config->GetSolution_AdjFileName();
	copy.assign(mesh_filename);
	copy.erase (copy.end()-4, copy.end());
	strcpy (cstr, copy.c_str()); 
	if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) sprintf (buffer, "_cd.dat"); 
	if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT) sprintf (buffer, "_cl.dat");
	if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT) sprintf (buffer, "_csf.dat"); 
	if (config->GetKind_ObjFunc() == PRESSURE_COEFFICIENT) sprintf (buffer, "_cp.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) sprintf (buffer, "_cmx.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) sprintf (buffer, "_cmy.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) sprintf (buffer, "_cmz.dat"); 
	if (config->GetKind_ObjFunc() == EFFICIENCY) sprintf (buffer, "_eff.dat"); 
	if (config->GetKind_ObjFunc() == ELECTRIC_CHARGE) sprintf (buffer, "_ec.dat"); 
  if (config->GetKind_ObjFunc() == FORCE_X_COEFFICIENT) sprintf (buffer, "_cfx.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Y_COEFFICIENT) sprintf (buffer, "_cfy.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Z_COEFFICIENT) sprintf (buffer, "_cfz.dat");
//  if (config->GetKind_ObjFunc() == ROTOR_EFFICIENCY) sprintf (buffer, "_rot_eff.dat");
	strcat(cstr, buffer);
	
	restart_file.open(cstr, ios::in);
	
	if (restart_file.fail()) {
		cout << "There is no flow restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get();
		exit(1); }
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
		getline(restart_file,text_line);
		istringstream point_line(text_line);
		if (nVar == 1) point_line >> index >> dummy[0] >>  AdjVar_Res[iPoint][0];		
		if (nVar == 4) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >>
									AdjVar_Res[iPoint][0] >> AdjVar_Res[iPoint][1] >> AdjVar_Res[iPoint][2] >> 
									AdjVar_Res[iPoint][3];
		if (nVar == 5) point_line >> index >> dummy[0] >> dummy[1] >> dummy[2] >> dummy[3] >> dummy[4] >> 
									AdjVar_Res[iPoint][0] >> AdjVar_Res[iPoint][1] >> AdjVar_Res[iPoint][2] >> 
									AdjVar_Res[iPoint][3] >> AdjVar_Res[iPoint][4];
	}
	restart_file.close();
}

void CGridAdaptation::SetComplete_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {	
		geometry->elem[iElem]->SetDivide (true);
		geometry->elem[iElem]->SetSemiDivide(false);
		geometry->elem[iElem]->SetDivStrength(1);
	}
}

void CGridAdaptation::SetNo_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {	
		geometry->elem[iElem]->SetDivide (false);
		geometry->elem[iElem]->SetSemiDivide(false);
		geometry->elem[iElem]->SetDivStrength(0);
	}
}

void CGridAdaptation::SetWake_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	double Coordx, Coordy, dist, wake = 0.5;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx = geometry->node[iPoint]->GetCoord(0);
			Coordy = geometry->node[iPoint]->GetCoord(1);
			dist = sqrt(Coordx*Coordx+Coordy*Coordy);
			if (dist < wake) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);
			}
			if ((Coordx > 0) && ((Coordy > -wake) && (Coordy < wake))) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);				
			}
		}
}

void CGridAdaptation::SetTwoPhase_Refinement(CGeometry *geometry, unsigned short strength) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	double Coordx, Coordy, wake = 0.1, offset = 1.0;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx = geometry->node[iPoint]->GetCoord(0);
			Coordy = geometry->node[iPoint]->GetCoord(1);
			if ((Coordy - offset > -wake) && (Coordy - offset < wake)) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);				
			}
		}
}

void CGridAdaptation::SetSupShock_Refinement(CGeometry *geometry, CConfig *config) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	double Coordx, Coordy;
	double mu_1 = asin(1/config->GetMach_FreeStreamND()-0.1);
	double mu_2 = asin(1/(config->GetMach_FreeStreamND()-0.7));
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx = geometry->node[iPoint]->GetCoord(0);
			Coordy = geometry->node[iPoint]->GetCoord(1);
			if (Coordy < 0.0)
			if ((Coordx > abs(Coordy/tan(mu_2))-0.25) && (Coordx < abs(Coordy/tan(mu_1))+1.25)) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);
			}
		}
}

void CGridAdaptation::SetNearField_Refinement(CGeometry *geometry, CConfig *config) {
	unsigned long iElem, iPoint;
	unsigned short iNode;
	double Coordx[4], Coordy[4], Coordz[4], NearField;
	
	NearField = config->GetPosition_Plane();
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetnNodes() > 4) { 
			cout << "The Near-Field adaptation is not prepared for elements with more than 4 points" << endl; 
			exit(1); 
		}
		
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			Coordx[iNode] = geometry->node[iPoint]->GetCoord(0);
			Coordy[iNode] = geometry->node[iPoint]->GetCoord(1);
			if (geometry->GetnDim() == 3) Coordz[iNode] = geometry->node[iPoint]->GetCoord(2);
		}
		
		if (geometry->GetnDim() == 2) {
			if ((Coordy[0] < NearField) && (Coordy[1] > NearField) || (Coordy[0] > NearField) && (Coordy[1] < NearField) ||
					(Coordy[0] < NearField) && (Coordy[2] > NearField) || (Coordy[0] > NearField) && (Coordy[2] < NearField) ||
					(Coordy[1] < NearField) && (Coordy[2] > NearField) || (Coordy[1] > NearField) && (Coordy[2] < NearField)) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);
			}
		}
		
		if (geometry->GetnDim() == 3) {
			if ((Coordz[0] < NearField) && (Coordz[1] > NearField) || (Coordz[0] > NearField) && (Coordz[1] < NearField) ||
					(Coordz[0] < NearField) && (Coordz[2] > NearField) || (Coordz[0] > NearField) && (Coordz[2] < NearField) ||
					(Coordz[0] < NearField) && (Coordz[3] > NearField) || (Coordz[0] > NearField) && (Coordz[3] < NearField) ||
					(Coordz[1] < NearField) && (Coordz[2] > NearField) || (Coordz[1] > NearField) && (Coordz[2] < NearField) ||
					(Coordz[1] < NearField) && (Coordz[3] > NearField) || (Coordz[1] > NearField) && (Coordz[3] < NearField)) {
				geometry->elem[iElem]->SetDivide (true);
				geometry->elem[iElem]->SetSemiDivide(false);
				geometry->elem[iElem]->SetDivStrength(1);
			}
		}
	}
}

void CGridAdaptation::SetTetraPattern (unsigned long edge_code[6], bool new_tetra[8][10], unsigned short &nTetra) {
	unsigned short iElem, nElem, iDiv, nDiv, iVar, iNode, edge_div[6][3], max_iVar;
	bool set_0, set_1, new_div, edge_code_bin[6];
	unsigned long max_node;
		
	
	for (iVar = 0; iVar < 6; iVar ++){
		if (edge_code[iVar] == 0) edge_code_bin[iVar] = 0;
		else edge_code_bin[iVar] = 1;
	}
	
	nDiv = 0;
	do { new_div = false;
		// Compute the greatest node at the divided edge
		max_node = 0; max_iVar = 0;
		for (iVar = 0; iVar < 6; iVar ++){
			max_node = max(max_node,edge_code[iVar]);
			if ( max_node == edge_code[iVar] ) max_iVar = iVar;
		}
		// If the edge is divided compse the vecto with the information of the division
		if (edge_code[max_iVar] > 0) {
			if (max_iVar == 0) { edge_div[nDiv][0] = 4; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 1; }
			if (max_iVar == 1) { edge_div[nDiv][0] = 5; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 2; }			
			if (max_iVar == 2) { edge_div[nDiv][0] = 6; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 3; }			
			if (max_iVar == 3) { edge_div[nDiv][0] = 7; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 2; }			
			if (max_iVar == 4) { edge_div[nDiv][0] = 8; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 3; }			
			if (max_iVar == 5) { edge_div[nDiv][0] = 9; edge_div[nDiv][1] = 2; edge_div[nDiv][2] = 3; }			
			nDiv++; new_div = true;
		}
		// In order to do not repeat the egde, restart the code
		edge_code[max_iVar] = 0;
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
		
			if (target_elem != -1){
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
	
// Special case (regular division of a face)
	if ((edge_code_bin[0] == 1) && (edge_code_bin[1] == 1) && (edge_code_bin[2] == 0) && (edge_code_bin[3] == 1) && (edge_code_bin[4] == 0) && (edge_code_bin[5] == 0)){
		new_tetra[0][0] = 1; new_tetra[0][1] = 0; new_tetra[0][2] = 0; new_tetra[0][3] = 1; 
		new_tetra[0][4] = 1; new_tetra[0][5] = 1; new_tetra[0][6] = 0; new_tetra[0][7] = 0; new_tetra[0][8] = 0; new_tetra[0][9] = 0;
		new_tetra[1][0] = 0; new_tetra[1][1] = 1; new_tetra[1][2] = 0; new_tetra[1][3] = 1; 
		new_tetra[1][4] = 1; new_tetra[1][5] = 0; new_tetra[1][6] = 0; new_tetra[1][7] = 1; new_tetra[1][8] = 0; new_tetra[1][9] = 0;
		new_tetra[2][0] = 0; new_tetra[2][1] = 0; new_tetra[2][2] = 0; new_tetra[2][3] = 1; 
		new_tetra[2][4] = 1; new_tetra[2][5] = 1; new_tetra[2][6] = 0; new_tetra[2][7] = 1; new_tetra[2][8] = 0; new_tetra[2][9] = 0;
		new_tetra[3][0] = 0; new_tetra[3][1] = 0; new_tetra[3][2] = 1; new_tetra[3][3] = 1; 
		new_tetra[3][4] = 0; new_tetra[3][5] = 1; new_tetra[3][6] = 0; new_tetra[3][7] = 1; new_tetra[3][8] = 0; new_tetra[3][9] = 0;		
		nElem = 4;
	}
	
	if ((edge_code_bin[0] == 0) && (edge_code_bin[1] == 0) && (edge_code_bin[2] == 0) && (edge_code_bin[3] == 1) && (edge_code_bin[4] == 1) && (edge_code_bin[5] == 1)){
		new_tetra[0][0] = 1; new_tetra[0][1] = 1; new_tetra[0][2] = 0; new_tetra[0][3] = 0; 
		new_tetra[0][4] = 0; new_tetra[0][5] = 0; new_tetra[0][6] = 0; new_tetra[0][7] = 1; new_tetra[0][8] = 1; new_tetra[0][9] = 0;
		new_tetra[1][0] = 1; new_tetra[1][1] = 0; new_tetra[1][2] = 1; new_tetra[1][3] = 0; 
		new_tetra[1][4] = 0; new_tetra[1][5] = 0; new_tetra[1][6] = 0; new_tetra[1][7] = 1; new_tetra[1][8] = 0; new_tetra[1][9] = 1;
		new_tetra[2][0] = 1; new_tetra[2][1] = 0; new_tetra[2][2] = 0; new_tetra[2][3] = 0; 
		new_tetra[2][4] = 0; new_tetra[2][5] = 0; new_tetra[2][6] = 0; new_tetra[2][7] = 1; new_tetra[2][8] = 1; new_tetra[2][9] = 1;
		new_tetra[3][0] = 1; new_tetra[3][1] = 0; new_tetra[3][2] = 0; new_tetra[3][3] = 1; 
		new_tetra[3][4] = 0; new_tetra[3][5] = 0; new_tetra[3][6] = 0; new_tetra[3][7] = 0; new_tetra[3][8] = 1; new_tetra[3][9] = 1;				
		nElem = 4;
	}
	
	if ((edge_code_bin[0] == 0) && (edge_code_bin[1] == 1) && (edge_code_bin[2] == 1) && (edge_code_bin[3] == 0) && (edge_code_bin[4] == 0) && (edge_code_bin[5] == 1)){
		new_tetra[0][0] = 0; new_tetra[0][1] = 1; new_tetra[0][2] = 1; new_tetra[0][3] = 0; 
		new_tetra[0][4] = 0; new_tetra[0][5] = 1; new_tetra[0][6] = 0; new_tetra[0][7] = 0; new_tetra[0][8] = 0; new_tetra[0][9] = 1;
		new_tetra[1][0] = 0; new_tetra[1][1] = 1; new_tetra[1][2] = 0; new_tetra[1][3] = 1; 
		new_tetra[1][4] = 0; new_tetra[1][5] = 0; new_tetra[1][6] = 1; new_tetra[1][7] = 0; new_tetra[1][8] = 0; new_tetra[1][9] = 1;
		new_tetra[2][0] = 0; new_tetra[2][1] = 1; new_tetra[2][2] = 0; new_tetra[2][3] = 0; 
		new_tetra[2][4] = 0; new_tetra[2][5] = 1; new_tetra[2][6] = 1; new_tetra[2][7] = 0; new_tetra[2][8] = 0; new_tetra[2][9] = 1;
		new_tetra[3][0] = 1; new_tetra[3][1] = 1; new_tetra[3][2] = 0; new_tetra[3][3] = 0; 
		new_tetra[3][4] = 0; new_tetra[3][5] = 1; new_tetra[3][6] = 1; new_tetra[3][7] = 0; new_tetra[3][8] = 0; new_tetra[3][9] = 0;				
		nElem = 4;
	}
	
	if ((edge_code_bin[0] == 1) && (edge_code_bin[1] == 0) && (edge_code_bin[2] == 1) && (edge_code_bin[3] == 0) && (edge_code_bin[4] == 1) && (edge_code_bin[5] == 0)){
		new_tetra[0][0] = 0; new_tetra[0][1] = 0; new_tetra[0][2] = 1; new_tetra[0][3] = 1; 
		new_tetra[0][4] = 0; new_tetra[0][5] = 0; new_tetra[0][6] = 1; new_tetra[0][7] = 0; new_tetra[0][8] = 1; new_tetra[0][9] = 0;
		new_tetra[1][0] = 1; new_tetra[1][1] = 0; new_tetra[1][2] = 1; new_tetra[1][3] = 0; 
		new_tetra[1][4] = 1; new_tetra[1][5] = 0; new_tetra[1][6] = 1; new_tetra[1][7] = 0; new_tetra[1][8] = 0; new_tetra[1][9] = 0;
		new_tetra[2][0] = 0; new_tetra[2][1] = 0; new_tetra[2][2] = 1; new_tetra[2][3] = 0; 
		new_tetra[2][4] = 1; new_tetra[2][5] = 0; new_tetra[2][6] = 1; new_tetra[2][7] = 0; new_tetra[2][8] = 1; new_tetra[2][9] = 0;
		new_tetra[3][0] = 0; new_tetra[3][1] = 1; new_tetra[3][2] = 1; new_tetra[3][3] = 0; 
		new_tetra[3][4] = 1; new_tetra[3][5] = 0; new_tetra[3][6] = 0; new_tetra[3][7] = 0; new_tetra[3][8] = 1; new_tetra[3][9] = 0;				
		nElem = 4;
	}
	
	nTetra = nElem;

}

void CGridAdaptation::SetTrianglePattern (unsigned long edge_code[3], bool new_triangle[4][6], unsigned short &nTriangle) {
	unsigned short iElem, nElem, iDiv, nDiv, iVar, iNode, edge_div[3][3], max_iVar;
	bool set_0, set_1, new_div;
	unsigned long max_node;
	
	nDiv = 0;
	do { new_div = false;
		// Compute the greatest node at the divided edge
		max_node = 0; max_iVar = 0;
		for (iVar = 0; iVar < 3; iVar ++){
			max_node = max(max_node,edge_code[iVar]);
			if ( max_node == edge_code[iVar] ) max_iVar = iVar;
		}
		// If the edge is divided compse the vecto with the information of the division
		if (edge_code[max_iVar] > 0) {
			if (max_iVar == 0) { edge_div[nDiv][0] = 3; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 1; }
			if (max_iVar == 1) { edge_div[nDiv][0] = 4; edge_div[nDiv][1] = 0; edge_div[nDiv][2] = 2; }			
			if (max_iVar == 2) { edge_div[nDiv][0] = 5; edge_div[nDiv][1] = 1; edge_div[nDiv][2] = 2; }			
			nDiv++; new_div = true;
		}
		// In order to do not repeat the egde, restart the code
		edge_code[max_iVar] = 0;
	} while (new_div);
	
	if (nDiv == 3) { // Regular division
		new_triangle[0][0] = 1; new_triangle[0][1] = 0; new_triangle[0][2] = 0; new_triangle[0][3] = 1; new_triangle[0][4] = 1; new_triangle[0][5] = 0; 
		new_triangle[1][0] = 0; new_triangle[1][1] = 1; new_triangle[1][2] = 0; new_triangle[1][3] = 1; new_triangle[1][4] = 0; new_triangle[1][5] = 1; 
		new_triangle[2][0] = 0; new_triangle[2][1] = 0; new_triangle[2][2] = 1; new_triangle[2][3] = 0; new_triangle[2][4] = 1; new_triangle[2][5] = 1; 
		new_triangle[3][0] = 0; new_triangle[3][1] = 0; new_triangle[3][2] = 0; new_triangle[3][3] = 1; new_triangle[3][4] = 1; new_triangle[3][5] = 1; 
		nElem = 4;
	}
	else {	// Irregular division
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
		
			if (target_elem != -1){
				for (iNode = 0; iNode < 6; iNode++)
					new_triangle[nElem][iNode] = new_triangle[target_elem][iNode];
				new_triangle[target_elem][edge_div[iDiv][0]] = true;
				new_triangle[target_elem][edge_div[iDiv][2]] = false;
				new_triangle[nElem][edge_div[iDiv][0]] = true;
				new_triangle[nElem][edge_div[iDiv][1]] = false;
				nElem++;
			}
		}
	}
	nTriangle = nElem;
}

void CGridAdaptation::SetHomothetic_Adaptation(CGeometry *geometry, CPhysicalGeometry *geo_adapt, 
										  CConfig *config){
	
	unsigned long iPoint, iElem, iEdge, iVertex, iPoint_new, iElem_new, ip_0, ip_1, ip_2, ip_3, 
					ip_01, ip_02, ip_03, ip_12, ip_13, ip_23, ip[4];
	unsigned short iDim, iVar, iMarker, counter, nTetra, iTetra, nTriangle, iTriangle, iNode;
	bool new_elem, new_tetra[8][10], new_triangle[4][6], PlaneAdaptation = false;
	double cog[3], NearField;
	long *EdgeNode = new long[geometry->GetnEdge()];
	bool *DivEdge = new bool[geometry->GetnEdge()];
	unsigned long *nElem_Bound_new = new unsigned long [geometry->GetnMarker()];
	unsigned long *iElem_Bound_new = new unsigned long [geometry->GetnMarker()];
	unsigned long (*Elem_NearField)[5] = NULL;
	unsigned long NearField_Trian = 0, NearField_Quad = 0;
	unsigned long nElem_NearField;
	
	NearField = config->GetPosition_Plane();
	
	bool Restart_Flow = false;
	bool Restart_Adjoint = false;
	bool Restart_Linear = false;
	
	if ((config->GetKind_Adaptation() == FULL_FLOW) ||
		(config->GetKind_Adaptation() == GRAD_FLOW) ||
		(config->GetKind_Adaptation() == FULL_ADJOINT) ||
		(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
		(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
		(config->GetKind_Adaptation() == FULL_LINEAR) ||
		(config->GetKind_Adaptation() == ROBUST) ||
		(config->GetKind_Adaptation() == REMAINING) ||
		(config->GetKind_Adaptation() == COMPUTABLE_ROBUST) ||
		(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Flow = true;
	
	if ((config->GetKind_Adaptation() == FULL_ADJOINT) ||
		(config->GetKind_Adaptation() == GRAD_ADJOINT) ||
		(config->GetKind_Adaptation() == GRAD_FLOW_ADJ) ||
		(config->GetKind_Adaptation() == ROBUST) ||
		(config->GetKind_Adaptation() == REMAINING) ||
		(config->GetKind_Adaptation() == COMPUTABLE_ROBUST) ||
		(config->GetKind_Adaptation() == COMPUTABLE)) Restart_Adjoint = true;

	if ((config->GetKind_Adaptation() == COMPUTABLE_ROBUST) ||
		(config->GetKind_Adaptation() == FULL_LINEAR)) Restart_Linear = true;

	
	//	Inicializa el vector que almacena los nuevos nodos de las aristas divididas
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge ++) EdgeNode[iEdge] = -1;
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge ++) DivEdge[iEdge] = false;
		
	//	Cuenta los elementos de cada tipo y marca las aristas divididas
	nElem_new = 0; nPoint_new = 0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
		if (geometry->elem[iElem]->GetDivStrength() == 1) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) 
				ip_3 = geometry->elem[iElem]->GetNode(3);
			
			DivEdge[geometry->FindEdge(ip_0,ip_1)] = true;
			DivEdge[geometry->FindEdge(ip_0,ip_2)] = true;
			DivEdge[geometry->FindEdge(ip_1,ip_2)] = true;
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) { 
				DivEdge[geometry->FindEdge(ip_0,ip_3)] = true;
				DivEdge[geometry->FindEdge(ip_1,ip_3)] = true;
				DivEdge[geometry->FindEdge(ip_2,ip_3)] = true;
			}
			geometry->elem[iElem]->SetDivStrength(1);
			geometry->elem[iElem]->SetDivide(true);
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				nElem_new = nElem_new + 3; nPoint_new = nPoint_new + 3; }
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) { 
				nElem_new = nElem_new + 7; nPoint_new = nPoint_new + 6; }
		}
	
	//	Para los tetrahedros. Divide todas las aristas de los elementos que tengan 4 o m‡s aristas divididas
	if (geometry->GetnDim() == 3) {
		do { new_elem = false;
			for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				ip_3 = geometry->elem[iElem]->GetNode(3);
			
				counter = 0;
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_0,ip_3)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_1,ip_3)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_2,ip_3)]) counter++;
						
				if ((counter > 3) && (!geometry->elem[iElem]->GetDivide())) { 
					DivEdge[geometry->FindEdge(ip_0,ip_1)] = true;
					DivEdge[geometry->FindEdge(ip_0,ip_2)] = true;
					DivEdge[geometry->FindEdge(ip_0,ip_3)] = true;
					DivEdge[geometry->FindEdge(ip_1,ip_2)] = true;
					DivEdge[geometry->FindEdge(ip_1,ip_3)] = true;
					DivEdge[geometry->FindEdge(ip_2,ip_3)] = true;
					geometry->elem[iElem]->SetDivStrength(1);
					geometry->elem[iElem]->SetDivide(true);
					geometry->elem[iElem]->SetSemiDivide(false);
					nElem_new = nElem_new + 7; nPoint_new = nPoint_new + 6;
					new_elem = true;
				}				
			}
		} while (new_elem);
	}
	
	if (geometry->GetnDim() == 2) {
		do { new_elem = false;
			for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
				ip_0 = geometry->elem[iElem]->GetNode(0);
				ip_1 = geometry->elem[iElem]->GetNode(1);
				ip_2 = geometry->elem[iElem]->GetNode(2);
				
				counter = 0;
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) counter++;
				if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) counter++;
				
				if ((counter == 3) && (!geometry->elem[iElem]->GetDivide())) { 
					DivEdge[geometry->FindEdge(ip_0,ip_1)] = true;
					DivEdge[geometry->FindEdge(ip_0,ip_2)] = true;
					DivEdge[geometry->FindEdge(ip_1,ip_2)] = true;
					geometry->elem[iElem]->SetDivStrength(1);
					geometry->elem[iElem]->SetDivide(true);
					geometry->elem[iElem]->SetSemiDivide(false);
					nElem_new = nElem_new + 3; nPoint_new = nPoint_new + 3;
					new_elem = true;
				}				
			}
		} while (new_elem);
	}
	
	//	Marca para adaptaci—n irregular todos los elementos con tres (s—lo tetrahedros), dos y una aristas divididas
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		ip_0 = geometry->elem[iElem]->GetNode(0);
		ip_1 = geometry->elem[iElem]->GetNode(1);
		ip_2 = geometry->elem[iElem]->GetNode(2);
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)
			ip_3 = geometry->elem[iElem]->GetNode(3);
		
		counter = 0;
		if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) counter++;
		if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) counter++;
		if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) counter++;
		
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
			if ((counter > 0) && (counter < 3)) {
				geometry->elem[iElem]->SetSemiDivide (true);
				geometry->elem[iElem]->SetDivStrength(0);
				geometry->elem[iElem]->SetDivide(false);
				nElem_new = nElem_new + 3; nPoint_new = nPoint_new + 3; // Peor de los casos
			}
		}
		
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
			if (DivEdge[geometry->FindEdge(ip_0,ip_3)]) counter++;
			if (DivEdge[geometry->FindEdge(ip_1,ip_3)]) counter++;
			if (DivEdge[geometry->FindEdge(ip_2,ip_3)]) counter++;		
			if ((counter > 0) && (counter < 4)) {
				geometry->elem[iElem]->SetSemiDivide (true);
				geometry->elem[iElem]->SetDivStrength(0);
				geometry->elem[iElem]->SetDivide(false);
				nElem_new = nElem_new + 7; nPoint_new = nPoint_new + 6; // Peor de los casos
			}
		}
	}				

	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		nElem_Bound_new[iMarker] = 4*geometry->GetnElem_Bound(iMarker); // Peor de los casos
	
	if (config->GetKind_Adaptation() == HORIZONTAL_PLANE) {
		PlaneAdaptation = true;
		Elem_NearField = new unsigned long[nPoint_new][5]; // Peor de los casos
		for (unsigned long iPoint_new = 0; iPoint_new < nPoint_new; iPoint_new++)
			Elem_NearField[iPoint_new][4] = 0;
		
		nElem_NearField = 0;
		counter = 0; 
	}

	//	Crea una copia de la malla inicial con espacio reservado
	geo_adapt->elem = new CPrimalGrid*[geometry->GetnElem() + nElem_new];
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)
			geo_adapt->elem[iElem] = new CTriangle(geometry->elem[iElem]->GetNode(0), 
													   geometry->elem[iElem]->GetNode(1), 
													   geometry->elem[iElem]->GetNode(2), 2);
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)
			geo_adapt->elem[iElem] = new CTetrahedron(geometry->elem[iElem]->GetNode(0),
												  geometry->elem[iElem]->GetNode(1),
												  geometry->elem[iElem]->GetNode(2),
												  geometry->elem[iElem]->GetNode(3));
	}
	
	geo_adapt->node = new CPoint*[geometry->GetnPoint() + nPoint_new];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		if (geometry->GetnDim() == 2)
			geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
												 geometry->node[iPoint]->GetCoord(1), iPoint, config);
		if (geometry->GetnDim() == 3)
			geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
												 geometry->node[iPoint]->GetCoord(1),
												 geometry->node[iPoint]->GetCoord(2), iPoint, config);
	}
	
	if (!PlaneAdaptation) geo_adapt->bound = new CPrimalGrid**[geometry->GetnMarker()];
	else geo_adapt->bound = new CPrimalGrid**[geometry->GetnMarker()+1];
	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++){
		geo_adapt->bound[iMarker] = new CPrimalGrid* [geometry->GetnElem_Bound(iMarker) + nElem_Bound_new[iMarker]];
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
				geo_adapt->bound[iMarker][iVertex] = new CLine(geometry->bound[iMarker][iVertex]->GetNode(0),
															   geometry->bound[iMarker][iVertex]->GetNode(1),2);
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
				geo_adapt->bound[iMarker][iVertex] = new CTriangle(geometry->bound[iMarker][iVertex]->GetNode(0),
															   geometry->bound[iMarker][iVertex]->GetNode(1), 
															   geometry->bound[iMarker][iVertex]->GetNode(2),3);
		}
	}

//	Crea e interpola las soluciones adjuntas en la malla adaptada
	if (Restart_Flow) {
		ConsVar_Adapt = new double* [geometry->GetnPoint() + nPoint_new];
		for (iPoint = 0; iPoint < geometry->GetnPoint()+ nPoint_new; iPoint ++)
			ConsVar_Adapt[iPoint] = new double [nVar];
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
			for (iVar = 0; iVar < nVar; iVar ++)
				ConsVar_Adapt[iPoint][iVar] = ConsVar_Sol[iPoint][iVar];
	}
	
	if (Restart_Adjoint) {
		AdjVar_Adapt = new double* [geometry->GetnPoint() + nPoint_new];
		for (iPoint = 0; iPoint < geometry->GetnPoint()+ nPoint_new; iPoint ++)
			AdjVar_Adapt[iPoint] = new double [nVar];
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
			for (iVar = 0; iVar < nVar; iVar ++)
				AdjVar_Adapt[iPoint][iVar] = AdjVar_Sol[iPoint][iVar];
	}
	
	if (Restart_Linear) {
		LinVar_Adapt = new double* [geometry->GetnPoint() + nPoint_new];
		for (iPoint = 0; iPoint < geometry->GetnPoint()+ nPoint_new; iPoint ++)
			LinVar_Adapt[iPoint] = new double [nVar];
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++)
			for (iVar = 0; iVar < nVar; iVar ++)
				LinVar_Adapt[iPoint][iVar] = LinVar_Sol[iPoint][iVar];
	}
		
	//	Indice inicial para empezar a a–adir elementos y nodos
	iPoint_new = geometry->GetnPoint()-1;
	iElem_new = geometry->GetnElem()-1;
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		iElem_Bound_new[iMarker] = geometry->GetnElem_Bound(iMarker)-1;

	//	Algoritmo de division de tetrahedros y triangulos
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetDivStrength() == 1) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)
				ip_3 = geometry->elem[iElem]->GetNode(3);
			
			// Antes de a~adir la arista chequea si ya ha sido divida y tenemos lozalizado el nodo de divisi'on
			if (EdgeNode[geometry->FindEdge(ip_0,ip_1)] == -1) {
				iPoint_new ++; ip_01 = iPoint_new; EdgeNode[geometry->FindEdge(ip_0,ip_1)] = ip_01;

				for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
					cog[iDim] = 0.5 * (geometry->node[ip_0]->GetCoord(iDim) + geometry->node[ip_1]->GetCoord(iDim));
				
				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_0]->GetCoord(1) < NearField) && (geometry->node[ip_1]->GetCoord(1) > NearField) || 
							(geometry->node[ip_0]->GetCoord(1) > NearField) && (geometry->node[ip_1]->GetCoord(1) < NearField)) {
							cog[1] = NearField;
							double x_0 = geometry->node[ip_0]->GetCoord(0);
							double x_1 = geometry->node[ip_1]->GetCoord(0);
							double y_0 = geometry->node[ip_0]->GetCoord(1);
							double y_1 = geometry->node[ip_1]->GetCoord(1);
							cog[0] = x_0 + (x_1-x_0)*(cog[1]-y_0)/(y_1-y_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_01;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_1]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_1]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_0]->GetCoord(0);
							double x_1 = geometry->node[ip_1]->GetCoord(0);
							double y_0 = geometry->node[ip_0]->GetCoord(1);
							double y_1 = geometry->node[ip_1]->GetCoord(1);
							double z_0 = geometry->node[ip_0]->GetCoord(2);
							double z_1 = geometry->node[ip_1]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_01;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
				if (geometry->GetnDim() == 2) geo_adapt->node[ip_01] = new CPoint(cog[0], cog[1], ip_01, config);
				if (geometry->GetnDim() == 3) geo_adapt->node[ip_01] = new CPoint(cog[0], cog[1], cog[2], ip_01, config);
				for (iVar = 0; iVar < nVar; iVar ++) {
					if (Restart_Flow) ConsVar_Adapt[ip_01][iVar] =  0.5 * (ConsVar_Adapt[ip_1][iVar]+ConsVar_Adapt[ip_0][iVar]);
					if (Restart_Adjoint) AdjVar_Adapt[ip_01][iVar] =  0.5 * (AdjVar_Adapt[ip_1][iVar]+AdjVar_Adapt[ip_0][iVar]);
					if (Restart_Linear) LinVar_Adapt[ip_01][iVar] =  0.5 * (LinVar_Adapt[ip_1][iVar]+LinVar_Adapt[ip_0][iVar]);
				}
			}
			else {
				ip_01 = EdgeNode[geometry->FindEdge(ip_0,ip_1)];

				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_0]->GetCoord(1) < NearField) && (geometry->node[ip_1]->GetCoord(1) > NearField) || 
								(geometry->node[ip_0]->GetCoord(1) > NearField) && (geometry->node[ip_1]->GetCoord(1) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_01;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_1]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_1]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_01;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
			}
			
			if (EdgeNode[geometry->FindEdge(ip_0,ip_2)] == -1) {
				iPoint_new ++; ip_02 = iPoint_new; EdgeNode[geometry->FindEdge(ip_0,ip_2)] = ip_02;

				for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
					cog[iDim] = 0.5 * (geometry->node[ip_0]->GetCoord(iDim) + geometry->node[ip_2]->GetCoord(iDim));
				
				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_0]->GetCoord(1) < NearField) && (geometry->node[ip_2]->GetCoord(1) > NearField) || 
								(geometry->node[ip_0]->GetCoord(1) > NearField) && (geometry->node[ip_2]->GetCoord(1) < NearField)) {
							cog[1] = NearField;
							double x_0 = geometry->node[ip_0]->GetCoord(0);
							double x_1 = geometry->node[ip_2]->GetCoord(0);
							double y_0 = geometry->node[ip_0]->GetCoord(1);
							double y_1 = geometry->node[ip_2]->GetCoord(1);
							cog[0] = x_0 + (x_1-x_0)*(cog[1]-y_0)/(y_1-y_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_02;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_2]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_2]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_0]->GetCoord(0);
							double x_1 = geometry->node[ip_2]->GetCoord(0);
							double y_0 = geometry->node[ip_0]->GetCoord(1);
							double y_1 = geometry->node[ip_2]->GetCoord(1);
							double z_0 = geometry->node[ip_0]->GetCoord(2);
							double z_1 = geometry->node[ip_2]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_02;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
				if (geometry->GetnDim() == 2) geo_adapt->node[ip_02] = new CPoint(cog[0], cog[1], ip_02, config);
				if (geometry->GetnDim() == 3) geo_adapt->node[ip_02] = new CPoint(cog[0], cog[1], cog[2], ip_02, config);
				for (iVar = 0; iVar < nVar; iVar ++){
					if (Restart_Flow) ConsVar_Adapt[ip_02][iVar] =  0.5 * (ConsVar_Adapt[ip_0][iVar]+ConsVar_Adapt[ip_2][iVar]);
					if (Restart_Adjoint) AdjVar_Adapt[ip_02][iVar] =  0.5 * (AdjVar_Adapt[ip_0][iVar]+AdjVar_Adapt[ip_2][iVar]);
					if (Restart_Linear) LinVar_Adapt[ip_02][iVar] =  0.5 * (LinVar_Adapt[ip_0][iVar]+LinVar_Adapt[ip_2][iVar]);
				}
			}
			else {
				ip_02 = EdgeNode[geometry->FindEdge(ip_0,ip_2)];

				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_0]->GetCoord(1) < NearField) && (geometry->node[ip_2]->GetCoord(1) > NearField) || 
								(geometry->node[ip_0]->GetCoord(1) > NearField) && (geometry->node[ip_2]->GetCoord(1) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_02;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_2]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_2]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_02;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
			}	
	
			if (EdgeNode[geometry->FindEdge(ip_1,ip_2)] == -1) {
				iPoint_new ++; ip_12 = iPoint_new; EdgeNode[geometry->FindEdge(ip_1,ip_2)] = ip_12;

				for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
					cog[iDim] = 0.5 * (geometry->node[ip_1]->GetCoord(iDim) + geometry->node[ip_2]->GetCoord(iDim));

				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_1]->GetCoord(1) < NearField) && (geometry->node[ip_2]->GetCoord(1) > NearField) || 
								(geometry->node[ip_1]->GetCoord(1) > NearField) && (geometry->node[ip_2]->GetCoord(1) < NearField)) {
							cog[1] = NearField;
							double x_0 = geometry->node[ip_1]->GetCoord(0);
							double x_1 = geometry->node[ip_2]->GetCoord(0);
							double y_0 = geometry->node[ip_1]->GetCoord(1);
							double y_1 = geometry->node[ip_2]->GetCoord(1);
							cog[0] = x_0 + (x_1-x_0)*(cog[1]-y_0)/(y_1-y_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_12;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_1]->GetCoord(2) < NearField) && (geometry->node[ip_2]->GetCoord(2) > NearField) || 
								(geometry->node[ip_1]->GetCoord(2) > NearField) && (geometry->node[ip_2]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_1]->GetCoord(0);
							double x_1 = geometry->node[ip_2]->GetCoord(0);
							double y_0 = geometry->node[ip_1]->GetCoord(1);
							double y_1 = geometry->node[ip_2]->GetCoord(1);
							double z_0 = geometry->node[ip_1]->GetCoord(2);
							double z_1 = geometry->node[ip_2]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_12;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
				if (geometry->GetnDim() == 2) geo_adapt->node[ip_12] = new CPoint(cog[0], cog[1], ip_12, config);
				if (geometry->GetnDim() == 3) geo_adapt->node[ip_12] = new CPoint(cog[0], cog[1], cog[2], ip_12, config);
				for (iVar = 0; iVar < nVar; iVar ++){
					if (Restart_Flow) ConsVar_Adapt[ip_12][iVar] =  0.5 * (ConsVar_Adapt[ip_1][iVar]+ConsVar_Adapt[ip_2][iVar]);
					if (Restart_Adjoint) AdjVar_Adapt[ip_12][iVar] =  0.5 * (AdjVar_Adapt[ip_1][iVar]+AdjVar_Adapt[ip_2][iVar]);
					if (Restart_Linear) LinVar_Adapt[ip_12][iVar] =  0.5 * (LinVar_Adapt[ip_1][iVar]+LinVar_Adapt[ip_2][iVar]);
				}
			}
			else {
				ip_12 = EdgeNode[geometry->FindEdge(ip_1,ip_2)];

				if (PlaneAdaptation) {
					if (geometry->GetnDim() == 2)
						if ((geometry->node[ip_1]->GetCoord(1) < NearField) && (geometry->node[ip_2]->GetCoord(1) > NearField) || 
								(geometry->node[ip_1]->GetCoord(1) > NearField) && (geometry->node[ip_2]->GetCoord(1) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_12;
							Elem_NearField[nElem_NearField][4]++;
						}
					if (geometry->GetnDim() == 3)
						if ((geometry->node[ip_1]->GetCoord(2) < NearField) && (geometry->node[ip_2]->GetCoord(2) > NearField) || 
								(geometry->node[ip_1]->GetCoord(2) > NearField) && (geometry->node[ip_2]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_12;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
			}

			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				if (EdgeNode[geometry->FindEdge(ip_0,ip_3)] == -1) {
					iPoint_new ++; ip_03 = iPoint_new; EdgeNode[geometry->FindEdge(ip_0,ip_3)] = ip_03;

					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						cog[iDim] = 0.5 * (geometry->node[ip_0]->GetCoord(iDim) + geometry->node[ip_3]->GetCoord(iDim));
					
					if (PlaneAdaptation)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_0]->GetCoord(0);
							double x_1 = geometry->node[ip_3]->GetCoord(0);
							double y_0 = geometry->node[ip_0]->GetCoord(1);
							double y_1 = geometry->node[ip_3]->GetCoord(1);
							double z_0 = geometry->node[ip_0]->GetCoord(2);
							double z_1 = geometry->node[ip_3]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_03;
							Elem_NearField[nElem_NearField][4]++;
						}

					geo_adapt->node[ip_03] = new CPoint(cog[0], cog[1], cog[2], ip_03, config);
					for (iVar = 0; iVar < nVar; iVar ++){
						if (Restart_Flow) ConsVar_Adapt[ip_03][iVar] =  0.5 * (ConsVar_Adapt[ip_0][iVar]+ConsVar_Adapt[ip_3][iVar]);
						if (Restart_Adjoint) AdjVar_Adapt[ip_03][iVar] =  0.5 * (AdjVar_Adapt[ip_0][iVar]+AdjVar_Adapt[ip_3][iVar]);
						if (Restart_Linear) LinVar_Adapt[ip_03][iVar] =  0.5 * (LinVar_Adapt[ip_0][iVar]+LinVar_Adapt[ip_3][iVar]);
					}
				} 
				else {
					ip_03 = EdgeNode[geometry->FindEdge(ip_0,ip_3)];
					if (PlaneAdaptation)
						if ((geometry->node[ip_0]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_0]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_03;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
			
				if (EdgeNode[geometry->FindEdge(ip_1,ip_3)] == -1) {
					iPoint_new ++; ip_13 = iPoint_new; EdgeNode[geometry->FindEdge(ip_1,ip_3)] = ip_13;

					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						cog[iDim] = 0.5 * (geometry->node[ip_1]->GetCoord(iDim) + geometry->node[ip_3]->GetCoord(iDim));
					
					if (PlaneAdaptation)
						if ((geometry->node[ip_1]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_1]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_1]->GetCoord(0);
							double x_1 = geometry->node[ip_3]->GetCoord(0);
							double y_0 = geometry->node[ip_1]->GetCoord(1);
							double y_1 = geometry->node[ip_3]->GetCoord(1);
							double z_0 = geometry->node[ip_1]->GetCoord(2);
							double z_1 = geometry->node[ip_3]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_13;
							Elem_NearField[nElem_NearField][4]++;
						}
					geo_adapt->node[ip_13] = new CPoint(cog[0], cog[1], cog[2], ip_13, config);
					for (iVar = 0; iVar < nVar; iVar ++){
						if (Restart_Flow) ConsVar_Adapt[ip_13][iVar] =  0.5 * (ConsVar_Adapt[ip_1][iVar]+ConsVar_Adapt[ip_3][iVar]);
						if (Restart_Adjoint) AdjVar_Adapt[ip_13][iVar] =  0.5 * (AdjVar_Adapt[ip_1][iVar]+AdjVar_Adapt[ip_3][iVar]);
						if (Restart_Linear) LinVar_Adapt[ip_13][iVar] =  0.5 * (LinVar_Adapt[ip_1][iVar]+LinVar_Adapt[ip_3][iVar]);
					}
				}
				else {
					ip_13 = EdgeNode[geometry->FindEdge(ip_1,ip_3)];
					if (PlaneAdaptation)
						if ((geometry->node[ip_1]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_1]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_13;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
						
				if (EdgeNode[geometry->FindEdge(ip_2,ip_3)] == -1) {
					iPoint_new ++; ip_23 = iPoint_new; EdgeNode[geometry->FindEdge(ip_2,ip_3)] = ip_23;

					for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
						cog[iDim] = 0.5 * (geometry->node[ip_2]->GetCoord(iDim) + geometry->node[ip_3]->GetCoord(iDim));
					
					if (PlaneAdaptation)
						if ((geometry->node[ip_2]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_2]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							cog[2] = NearField;
							double x_0 = geometry->node[ip_2]->GetCoord(0);
							double x_1 = geometry->node[ip_3]->GetCoord(0);
							double y_0 = geometry->node[ip_2]->GetCoord(1);
							double y_1 = geometry->node[ip_3]->GetCoord(1);
							double z_0 = geometry->node[ip_2]->GetCoord(2);
							double z_1 = geometry->node[ip_3]->GetCoord(2);
							cog[0] = x_0 + (x_1-x_0)*(cog[2]-z_0)/(z_1-z_0);
							cog[1] = y_0 + (y_1-y_0)*(cog[2]-z_0)/(z_1-z_0);
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_23;
							Elem_NearField[nElem_NearField][4]++;
						}
					geo_adapt->node[ip_23] = new CPoint(cog[0], cog[1], cog[2], ip_23, config);
					for (iVar = 0; iVar < nVar; iVar ++){
						if (Restart_Flow) ConsVar_Adapt[ip_23][iVar] =  0.5 * (ConsVar_Adapt[ip_2][iVar]+ConsVar_Adapt[ip_3][iVar]);
						if (Restart_Adjoint) AdjVar_Adapt[ip_23][iVar] =  0.5 * (AdjVar_Adapt[ip_2][iVar]+AdjVar_Adapt[ip_3][iVar]);
						if (Restart_Linear) LinVar_Adapt[ip_23][iVar] =  0.5 * (LinVar_Adapt[ip_2][iVar]+LinVar_Adapt[ip_3][iVar]);
					}
				}
				else {
					ip_23 = EdgeNode[geometry->FindEdge(ip_2,ip_3)];
					if (PlaneAdaptation)
						if ((geometry->node[ip_2]->GetCoord(2) < NearField) && (geometry->node[ip_3]->GetCoord(2) > NearField) || 
								(geometry->node[ip_2]->GetCoord(2) > NearField) && (geometry->node[ip_3]->GetCoord(2) < NearField)) {
							Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = ip_23;
							Elem_NearField[nElem_NearField][4]++;
						}
				}
			}
			
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				geo_adapt->elem[iElem] =				   new CTriangle(ip_01,ip_12,ip_02,2);
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTriangle(ip_0,ip_01,ip_02,2); 
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTriangle(ip_1,ip_12,ip_01,2);  
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTriangle(ip_2,ip_02,ip_12,2);
			}
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				geo_adapt->elem[iElem] =				   new CTetrahedron(ip_0,ip_01,ip_03,ip_02);		
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_1,ip_12,ip_01,ip_13);
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_2,ip_23,ip_12,ip_02);
				iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_3,ip_03,ip_13,ip_23);
				if ((!PlaneAdaptation) || ((PlaneAdaptation) && (Elem_NearField[nElem_NearField][4] != 4))) {
					iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_01,ip_02,ip_12,ip_03);
					iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_01,ip_13,ip_12,ip_03);
					iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_13,ip_23,ip_12,ip_03);
					iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_23,ip_02,ip_12,ip_03);
				}
			}
			
			if (PlaneAdaptation) {
				
				/*--- If we have a rectangle we must divide it into two 
				 triangles using the conectivity of the tetra ---*/
				if (Elem_NearField[nElem_NearField][4] == 4) {
					
					short combination = -1;
					unsigned long ip_0t = Elem_NearField[nElem_NearField][0];
					unsigned long ip_1t = Elem_NearField[nElem_NearField][1];
					unsigned long ip_2t = Elem_NearField[nElem_NearField][2];
					unsigned long ip_3t = Elem_NearField[nElem_NearField][3];
					
					bool triangle_0 = false; bool triangle_1 = false; bool triangle_2 = false; bool triangle_3 = false;
					
					if ((!triangle_0) && ((ip_0t == ip_01) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_01) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_01) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))) {
						combination = 0; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_01) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_01) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_01) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))
							&& ((ip_0t == ip_01) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))) {
						combination = 0; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_01) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_01) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_12) || (ip_2t == ip_03))) {
						combination = 0; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_12) || (ip_3t == ip_03))
							&& ((ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_12) || (ip_0t == ip_03))) {
						combination = 0; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_13) || (ip_0t == ip_23) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_13) || (ip_1t == ip_23) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_13) || (ip_2t == ip_23) || (ip_2t == ip_12) || (ip_2t == ip_03))) {
						combination = 0; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_13) || (ip_1t == ip_23) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_13) || (ip_2t == ip_23) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_13) || (ip_3t == ip_23) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_13) || (ip_2t == ip_23) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_13) || (ip_3t == ip_23) || (ip_3t == ip_12) || (ip_3t == ip_03))
							&& ((ip_0t == ip_13) || (ip_0t == ip_23) || (ip_0t == ip_12) || (ip_0t == ip_03))) {
						combination = 0; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_13) || (ip_0t == ip_23) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_13) || (ip_1t == ip_23) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_3t == ip_13) || (ip_3t == ip_23) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_23) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_23) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_23) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))) {
						combination = 0; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_23) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_2t == ip_23) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_23) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_23) || (ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_03))
							&& ((ip_3t == ip_23) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))
							&& ((ip_0t == ip_23) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))) {
						combination = 0; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_23) || (ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_03))
							&& ((ip_1t == ip_23) || (ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_03))
							&& ((ip_3t == ip_23) || (ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_03))) {
						combination = 0; triangle_3 = true; }
					
					if (combination == 0) goto stop_combination;
					
					if ((!triangle_0) && ((ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_01) || (ip_2t == ip_23))) {
						combination = 1; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_02) || (ip_2t == ip_12) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_01) || (ip_3t == ip_23))
							&& ((ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_01) || (ip_0t == ip_23))) {
						combination = 1; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_02) || (ip_0t == ip_12) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_02) || (ip_1t == ip_12) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_3t == ip_02) || (ip_3t == ip_12) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_01) || (ip_2t == ip_23))) {
						combination = 1; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_01) || (ip_3t == ip_23))
							&& ((ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_01) || (ip_0t == ip_23))) {
						combination = 1; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_13) || (ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_13) || (ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_13) || (ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_23))) {
						combination = 1; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_13) || (ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_13) || (ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_13) || (ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_13) || (ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_13) || (ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_23))
							&& ((ip_0t == ip_13) || (ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_23))) {
						combination = 1; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_13) || (ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_13) || (ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_3t == ip_13) || (ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_03) || (ip_0t == ip_02) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_03) || (ip_1t == ip_02) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_03) || (ip_2t == ip_02) || (ip_2t == ip_01) || (ip_2t == ip_23))) {
						combination = 1; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_03) || (ip_1t == ip_02) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_2t == ip_03) || (ip_2t == ip_02) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_03) || (ip_3t == ip_02) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_03) || (ip_2t == ip_02) || (ip_2t == ip_01) || (ip_2t == ip_23))
							&& ((ip_3t == ip_03) || (ip_3t == ip_02) || (ip_3t == ip_01) || (ip_3t == ip_23))
							&& ((ip_0t == ip_03) || (ip_0t == ip_02) || (ip_0t == ip_01) || (ip_0t == ip_23))) {
						combination = 1; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_03) || (ip_0t == ip_02) || (ip_0t == ip_01) || (ip_0t == ip_23))
							&& ((ip_1t == ip_03) || (ip_1t == ip_02) || (ip_1t == ip_01) || (ip_1t == ip_23))
							&& ((ip_3t == ip_03) || (ip_3t == ip_02) || (ip_3t == ip_01) || (ip_3t == ip_23))) {
						combination = 1; triangle_3 = true; }
					
					if (combination == 1) goto stop_combination;
					
					if ((!triangle_0) && ((ip_0t == ip_01) || (ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_01) || (ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_01) || (ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_02))) {
						combination = 2; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_01) || (ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_01) || (ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_01) || (ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_01) || (ip_2t == ip_12) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_01) || (ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_02))
							&& ((ip_0t == ip_01) || (ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_02))) {
						combination = 2; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_01) || (ip_0t == ip_12) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_01) || (ip_1t == ip_12) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_3t == ip_01) || (ip_3t == ip_12) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_3 = true; }
					
					if ((!triangle_0)
							&& ((ip_0t == ip_12) || (ip_0t == ip_23) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_12) || (ip_1t == ip_23) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_12) || (ip_2t == ip_23) || (ip_2t == ip_13) || (ip_2t == ip_02))) {
						combination = 2; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_12) || (ip_1t == ip_23) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_12) || (ip_2t == ip_23) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_12) || (ip_3t == ip_23) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_12) || (ip_2t == ip_23) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_12) || (ip_3t == ip_23) || (ip_3t == ip_13) || (ip_3t == ip_02))
							&& ((ip_0t == ip_12) || (ip_0t == ip_23) || (ip_0t == ip_13) || (ip_0t == ip_02))) {
						combination = 2; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_12) || (ip_0t == ip_23) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_12) || (ip_1t == ip_23) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_3t == ip_12) || (ip_3t == ip_23) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_23) || (ip_0t == ip_03) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_23) || (ip_1t == ip_03) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_23) || (ip_2t == ip_03) || (ip_2t == ip_13) || (ip_2t == ip_02))) {
						combination = 2; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_23) || (ip_1t == ip_03) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_23) || (ip_2t == ip_03) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_23) || (ip_3t == ip_03) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_23) || (ip_2t == ip_03) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_23) || (ip_3t == ip_03) || (ip_3t == ip_13) || (ip_3t == ip_02))
							&& ((ip_0t == ip_23) || (ip_0t == ip_03) || (ip_0t == ip_13) || (ip_0t == ip_02))) {
						combination = 2; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_23) || (ip_0t == ip_03) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_23) || (ip_1t == ip_03) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_3t == ip_23) || (ip_3t == ip_03) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_3 = true; }
					
					if ((!triangle_0) && ((ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_02))) {
						combination = 2; triangle_0 = true; }
					
					if ((!triangle_1) && ((ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_1 = true; }
					
					if ((!triangle_2) && ((ip_2t == ip_03) || (ip_2t == ip_01) || (ip_2t == ip_13) || (ip_2t == ip_02))
							&& ((ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_02))
							&& ((ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_02))) {
						combination = 2; triangle_2 = true; }
					
					if ((!triangle_3) && ((ip_0t == ip_03) || (ip_0t == ip_01) || (ip_0t == ip_13) || (ip_0t == ip_02))
							&& ((ip_1t == ip_03) || (ip_1t == ip_01) || (ip_1t == ip_13) || (ip_1t == ip_02))
							&& ((ip_3t == ip_03) || (ip_3t == ip_01) || (ip_3t == ip_13) || (ip_3t == ip_02))) {
						combination = 2; triangle_3 = true; }
					
					if (combination == 2) goto stop_combination;
					
				stop_combination:
					
					/*--- Change the combination in order to have a well defined plane ---*/
					
					if (combination == 0) {
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_01,ip_02,ip_12,ip_03);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_01,ip_13,ip_12,ip_03);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_13,ip_23,ip_12,ip_03);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_23,ip_02,ip_12,ip_03);
					}
					
					if (combination == 1) {
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_02,ip_12,ip_01,ip_23);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_12,ip_13,ip_01,ip_23);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_13,ip_03,ip_01,ip_23);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_03,ip_02,ip_01,ip_23);
					}
					
					if (combination == 2) {
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_01,ip_12,ip_13,ip_02);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_12,ip_23,ip_13,ip_02);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_23,ip_03,ip_13,ip_02);
						iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip_03,ip_01,ip_13,ip_02);
					}
					
					/*--- Create two triangles using the right combination ---*/
					bool JustOne = true;
					
					if (triangle_0) {
						Elem_NearField[nElem_NearField][0] = ip_0t;
						Elem_NearField[nElem_NearField][1] = ip_1t;
						Elem_NearField[nElem_NearField][2] = ip_2t;
						Elem_NearField[nElem_NearField][4] = 3;
						if (JustOne) { nElem_NearField++; JustOne = false; }
					}
					if (triangle_1) {
						Elem_NearField[nElem_NearField][0] = ip_1t;
						Elem_NearField[nElem_NearField][1] = ip_2t;
						Elem_NearField[nElem_NearField][2] = ip_3t;
						Elem_NearField[nElem_NearField][4] = 3;
						if (JustOne) { nElem_NearField++; JustOne = false; }					
					}
					if (triangle_2) {
						Elem_NearField[nElem_NearField][0] = ip_2t;
						Elem_NearField[nElem_NearField][1] = ip_3t;
						Elem_NearField[nElem_NearField][2] = ip_0t;
						Elem_NearField[nElem_NearField][4] = 3;
						if (JustOne) { nElem_NearField++; JustOne = false; }					
					}
					if (triangle_3) {
						Elem_NearField[nElem_NearField][0] = ip_0t;
						Elem_NearField[nElem_NearField][1] = ip_1t;
						Elem_NearField[nElem_NearField][2] = ip_3t;
						Elem_NearField[nElem_NearField][4] = 3;
						if (JustOne) { nElem_NearField++; JustOne = false; }
					}					
				}
				
				if (Elem_NearField[nElem_NearField][4] != 0)
					nElem_NearField++;
			}
		}
	}
	
	cout << "Regular adaptation finished. " <<iElem_new-geometry->GetnElem()+1<<" elements and "
		<<iPoint_new-geometry->GetnPoint()+1<<" points have been added." <<endl;
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		
		if (geometry->elem[iElem]->GetSemiDivide()) {
			ip_0 = geometry->elem[iElem]->GetNode(0);
			ip_1 = geometry->elem[iElem]->GetNode(1);
			ip_2 = geometry->elem[iElem]->GetNode(2);
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)
				ip_3 = geometry->elem[iElem]->GetNode(3);
						
			if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) {
				unsigned long edge_code[3] = {0, 0, 0};
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) edge_code[0] = geometry->FindEdge(ip_0,ip_1);
				if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) edge_code[1] = geometry->FindEdge(ip_0,ip_2);
				if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) edge_code[2] = geometry->FindEdge(ip_1,ip_2);
				SetTrianglePattern (edge_code, new_triangle, nTriangle);
				for (iTriangle = 0; iTriangle < nTriangle; iTriangle++) {
					iVar = 0;
					for (iNode = 0; iNode < 6; iNode++)
						if (new_triangle[iTriangle][iNode] == 1) {
							if (iNode == 0) ip[iVar] = ip_0;
							if (iNode == 1) ip[iVar] = ip_1;
							if (iNode == 2) ip[iVar] = ip_2;
							if (iNode == 3) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_1)];
							if (iNode == 4) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_2)];
							if (iNode == 5) ip[iVar] = EdgeNode[geometry->FindEdge(ip_1,ip_2)];
							iVar++;
						}	
					if (iTriangle == 0) geo_adapt->elem[iElem] = new CTriangle(ip[0],ip[1],ip[2],2);
					else { iElem_new ++; geo_adapt->elem[iElem_new] = new CTriangle(ip[0],ip[1],ip[2],2); }
				}
			}
			if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
				unsigned long edge_code[6] = {0, 0, 0, 0, 0 ,0};
				counter = 0;
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) {edge_code[0] = geometry->FindEdge(ip_0,ip_1); counter++;}
				if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) {edge_code[1] = geometry->FindEdge(ip_0,ip_2); counter++;}
				if (DivEdge[geometry->FindEdge(ip_0,ip_3)]) {edge_code[2] = geometry->FindEdge(ip_0,ip_3); counter++;}
				if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) {edge_code[3] = geometry->FindEdge(ip_1,ip_2); counter++;}
				if (DivEdge[geometry->FindEdge(ip_1,ip_3)]) {edge_code[4] = geometry->FindEdge(ip_1,ip_3); counter++;}
				if (DivEdge[geometry->FindEdge(ip_2,ip_3)]) {edge_code[5] = geometry->FindEdge(ip_2,ip_3); counter++;}
				
				SetTetraPattern (edge_code, new_tetra, nTetra);
				for (iTetra = 0; iTetra < nTetra; iTetra++) {
					iVar = 0;
					for (iNode = 0; iNode < 10; iNode++)
						if (new_tetra[iTetra][iNode] == 1) {
							if (iNode == 0) ip[iVar] = ip_0;
							if (iNode == 1) ip[iVar] = ip_1;
							if (iNode == 2) ip[iVar] = ip_2;
							if (iNode == 3) ip[iVar] = ip_3;
							if (iNode == 4) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_1)];
							if (iNode == 5) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_2)];
							if (iNode == 6) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_3)];
							if (iNode == 7) ip[iVar] = EdgeNode[geometry->FindEdge(ip_1,ip_2)];
							if (iNode == 8) ip[iVar] = EdgeNode[geometry->FindEdge(ip_1,ip_3)];
							if (iNode == 9) ip[iVar] = EdgeNode[geometry->FindEdge(ip_2,ip_3)];
							iVar++;
						}	
					if (iTetra == 0) geo_adapt->elem[iElem] = new CTetrahedron(ip[0],ip[1],ip[2],ip[3]);
					else { iElem_new ++; geo_adapt->elem[iElem_new] = new CTetrahedron(ip[0],ip[1],ip[2],ip[3]); }
				}
			}
		}
	}
	
	cout << "Irregular adaptation finished. " <<iElem_new-geometry->GetnElem()+1<<" elements have been added."<< endl;

	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++){
			ip_0 = geometry->bound[iMarker][iVertex]->GetNode(0);
			ip_1 = geometry->bound[iMarker][iVertex]->GetNode(1);

			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
				ip_2 = geometry->bound[iMarker][iVertex]->GetNode(2);
			
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE) {
				unsigned long edge_code[1] = {0}, ip_01;
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) {
					edge_code[0] = geometry->FindEdge(ip_0,ip_1);
					ip_01 = EdgeNode[edge_code[0]];
					geo_adapt->bound[iMarker][iVertex] = new CLine(ip_0,ip_01,2);
					iElem_Bound_new[iMarker]++; 
					geo_adapt->bound[iMarker][iElem_Bound_new[iMarker]] = new CLine(ip_01,ip_1,2);
				}
			}
				
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE) {
				unsigned long edge_code[3] = {0, 0, 0};
				counter = 0;
				if (DivEdge[geometry->FindEdge(ip_0,ip_1)]) { edge_code[0] = geometry->FindEdge(ip_0,ip_1); counter++;}
				if (DivEdge[geometry->FindEdge(ip_0,ip_2)]) { edge_code[1] = geometry->FindEdge(ip_0,ip_2); counter++;}
				if (DivEdge[geometry->FindEdge(ip_1,ip_2)]) { edge_code[2] = geometry->FindEdge(ip_1,ip_2); counter++;}
			
				if (counter != 0) {
					SetTrianglePattern (edge_code, new_triangle, nTriangle);
					for (iTriangle = 0; iTriangle < nTriangle; iTriangle++) {
						iVar = 0;
						for (iNode = 0; iNode < 6; iNode++)
							if (new_triangle[iTriangle][iNode] == 1) {
								if (iNode == 0) ip[iVar] = ip_0;
								if (iNode == 1) ip[iVar] = ip_1;
								if (iNode == 2) ip[iVar] = ip_2;
								if (iNode == 3) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_1)];
								if (iNode == 4) ip[iVar] = EdgeNode[geometry->FindEdge(ip_0,ip_2)];
								if (iNode == 5) ip[iVar] = EdgeNode[geometry->FindEdge(ip_1,ip_2)];
								iVar++;
							}										
						if (iTriangle == 0) geo_adapt->bound[iMarker][iVertex] = new CTriangle(ip[0],ip[1],ip[2],3);
						else { iElem_Bound_new[iMarker]++; geo_adapt->bound[iMarker][iElem_Bound_new[iMarker]] = new CTriangle(ip[0],ip[1],ip[2],3); }
					}
				}
			}
		}	
		cout << "Boundary adaptation finished. " <<iElem_Bound_new[iMarker]-geometry->GetnElem_Bound(iMarker)+1<< 
		" elems have been added to "<<iMarker<<"."<< endl;
	}

	/*--- There is an analytical definition of the surfaces ---*/
	if ((config->GetAnalytical_Surface() != NONE) && (config->GetKind_Adaptation() != HORIZONTAL_PLANE)) {
		double *Normal, Coord[3];
		
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
			for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
				
				ip_0 = geometry->bound[iMarker][iVertex]->GetNode(0);
				ip_1 = geometry->bound[iMarker][iVertex]->GetNode(1);
				Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

				if (EdgeNode[geometry->FindEdge(ip_0,ip_1)] != -1) {
					
					ip_01 = EdgeNode[geometry->FindEdge(ip_0,ip_1)];
					Coord[0] = geo_adapt->node[ip_01]->GetCoord(0);
					Coord[1] = geo_adapt->node[ip_01]->GetCoord(1);
					if (geometry->GetnDim() == 3) Coord[2] = geo_adapt->node[ip_01]->GetCoord(2);

					
					if ((config->GetMarker_All_Boundary(iMarker) == EULER_WALL) || (config->GetMarker_All_Boundary(iMarker) == NO_SLIP_WALL)) {
						if (config->GetAnalytical_Surface() == NACA0012_AIRFOIL) {
							double Ya = 0.0 / 100.0;
							double Xa = 0.0 / 10.0;
							double t = 12.0 / 100.0;
							
							double Ycurv, Yesp, y_NACA0012;
							if  (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
							else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
							
							Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
									  1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
							
							double y_1 = Ycurv + Yesp;
							double y_2 = Ycurv - Yesp;
							double d_1 = abs(y_1 - Coord[1]);
							double d_2 = abs(y_2 - Coord[1]);
							if (d_1 < d_2) y_NACA0012 = y_1;
							else y_NACA0012 = y_2;

							geo_adapt->node[ip_01]->SetCoord(1, y_NACA0012);
						}
						if (config->GetAnalytical_Surface() == NACA4412_AIRFOIL) {
							double Ya = 4.0 / 100.0;
							double Xa = 4.0 / 10.0;
							double t = 12.0 / 100.0;
							
							double Ycurv, Yesp, y_NACA4412;
							if  (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
							else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
							
							Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
									  1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
							
							double y_1 = Ycurv + Yesp;
							double y_2 = Ycurv - Yesp;
							double d_1 = abs(y_1 - Coord[1]);
							double d_2 = abs(y_2 - Coord[1]);
							if (d_1 < d_2) y_NACA4412 = y_1;
							else y_NACA4412 = y_2;

							geo_adapt->node[ip_01]->SetCoord(1, y_NACA4412);
						}
						if (config->GetAnalytical_Surface() == CYLINDER) {
							double x_Circle, y_Circle, radius = 0.5;
							x_Circle = abs(radius * cos(atan(Coord[1]/(Coord[0]-0.5))))*(Coord[0]-0.5)/abs((Coord[0]-0.5)+1E-16)+0.5;
							y_Circle = abs(radius * sin(atan(Coord[1]/(Coord[0]-0.5))))*Coord[1]/abs(Coord[1]+1E-16);
							
							geo_adapt->node[ip_01]->SetCoord(0, x_Circle);
							geo_adapt->node[ip_01]->SetCoord(1, y_Circle);
						}
						if (config->GetAnalytical_Surface() == BIPARABOLIC) {
							
							double c = 0.5;
							double t = 1.0 / 100.0;
							double y_BIPARABOLIC;
							
							if (Coord[1] > 0) {
								y_BIPARABOLIC =  t*(Coord[0]*Coord[0]-Coord[0])/(2.0*(c*c-c));
							}
							if (Coord[1] < 0) {
								y_BIPARABOLIC =  t*(Coord[0]-Coord[0]*Coord[0])/(2.0*(c*c-c));
							}
							
							geo_adapt->node[ip_01]->SetCoord(1, y_BIPARABOLIC);
							
						}
					}
					
					if (config->GetMarker_All_Boundary(iMarker) == FAR_FIELD) {
						if ((config->GetAnalytical_Surface() == NACA0012_AIRFOIL) || 
								(config->GetAnalytical_Surface() == NACA4412_AIRFOIL) ||
								(config->GetAnalytical_Surface() == BIPARABOLIC)) {
							double x_Circle, y_Circle, radius = 20.0;
							x_Circle = abs(radius * cos(atan(Coord[1]/Coord[0])))*Coord[0]/abs(Coord[0]+1E-16);
							y_Circle = abs(radius * sin(atan(Coord[1]/Coord[0])))*Coord[1]/abs(Coord[1]+1E-16);

							geo_adapt->node[ip_01]->SetCoord(0, x_Circle);
							geo_adapt->node[ip_01]->SetCoord(1, y_Circle);
						}
						if (config->GetAnalytical_Surface() == CYLINDER) {
							double x_Circle, y_Circle, radius = 15.0;
							x_Circle = abs(radius * cos(atan(Coord[1]/Coord[0])))*Coord[0]/abs(Coord[0]+1E-16);
							y_Circle = abs(radius * sin(atan(Coord[1]/Coord[0])))*Coord[1]/abs(Coord[1]+1E-16);
							
							geo_adapt->node[ip_01]->SetCoord(0, x_Circle);
							geo_adapt->node[ip_01]->SetCoord(1, y_Circle);
						}
					}
				}
			}
	}
				
/*	double x_original, y_original;
	long jVertex;
	int ibcbeg, ibcend;
	double t_inv[100], t[100], tval, y[100], y_inv[100], ybcbeg, ybcend, yppval, ypval, yval;
	
	for (iMarker=0; iMarker < geometry->GetnMarker(); iMarker++)
		for (iVertex=0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
			ip_0 = geo_adapt->bound[iMarker][iVertex]->GetNode(0);
			ip_1 = geo_adapt->bound[iMarker][iVertex]->GetNode(1);
			
			if (EdgeNode[geometry->FindEdge(ip_0,ip_1)] != -1) {
				double *ypp;
				ip_01 = EdgeNode[geometry->FindEdge(ip_0,ip_1)];
				x_original = geo_adapt->node[ip_01]->GetCoord(0);
				y_original = geo_adapt->node[ip_01]->GetCoord(1);
				
				int N = 0;
				for (jVertex = iVertex-7; jVertex < iVertex+7; jVertex++)				
					if ((jVertex >= 0) && (jVertex < geometry->GetnElem_Bound(iMarker))) {
						t[N] = geo_adapt->node[geometry->Vertex[iMarker][jVertex]->GetNode()]->GetCoord(0);	
						y[N] = geo_adapt->node[geometry->Vertex[iMarker][jVertex]->GetNode()]->GetCoord(1);
						N = N+1; }
				
				if (N == 0)
					for (jVertex = iVertex+7; jVertex < iVertex-7; jVertex--)
						if ((jVertex >= 0) && (jVertex < geometry->GetnElem_Bound(iMarker))) {
							t[N] = geo_adapt->node[geometry->Vertex[iMarker][jVertex]->GetNode()]->GetCoord(0);	
							y[N] = geo_adapt->node[geometry->Vertex[iMarker][jVertex]->GetNode()]->GetCoord(1);
							//							cout << iMarker <<" " <<  iVertex << " " << jVertex <<" "<< geometry->GetnElem_Bound(iMarker)  <<" "<< t[N] <<" "<< y[N]<< endl;
							N = N+1; }
				
				bool increasing_1 = true; bool increasing_2 = true; bool increasing_3 = true; bool increasing_4 = true;
				
				for (int i = 0; i<N-1; i++) if (t[i]>t[i+1]) increasing_1 = false;
				
				if (!increasing_1) {
					for (int i = 0; i<N; i++) t[i] = -t[i];
					for (int i = 0; i < N-1; i++) if (t[i]>t[i+1]) increasing_2 = false;
					if (!increasing_2) {
						for (int i = 0; i<N; i++) t[i] = -t[i];
						for (int i = 0; i<N; i++) { t_inv[i] = t[i]; y_inv[i] = y[i]; }
						for (int i = 0; i<N; i++) { t[i] = y_inv[i]; y[i] = t_inv[i]; }
						for (int i = 0; i<N-1; i++) if (t[i]>t[i+1]) increasing_3 = false;	
						if (!increasing_3) {
							for (int i = 0; i < N; i++) t[i] = -t[i];
							for (int i = 0; i < N-1; i++) if (t[i]>t[i+1]) increasing_4 = false;
							if (!increasing_4) 
								cout << "ERROR EN LA INTERPOLACION. REDUZC EN NUMERO DE SPLINES";
							
						}
					}
				}
				
				ibcbeg = 0; ibcend = 0; ybcbeg = 0.0; ybcend = 0.0;
				
				if (increasing_1 && increasing_2 && increasing_3 && increasing_4 && (N>10)) {
					ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );						
					yval = spline_cubic_val ( N, t, x_original, y, ypp, &ypval, &yppval );
					geo_adapt->node[ip_01]->SetCoord(0, x_original);
					geo_adapt->node[ip_01]->SetCoord(1, yval);
					delete [] ypp; }
				
				if (!increasing_1 && increasing_2 && increasing_3 && increasing_4 && (N>10)) {
					ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );											
					yval = spline_cubic_val ( N, t, -x_original, y, ypp, &ypval, &yppval );					
					geo_adapt->node[ip_01]->SetCoord(0, x_original);
					geo_adapt->node[ip_01]->SetCoord(1, yval);
					delete [] ypp; }
				
				if (!increasing_1 && !increasing_2 && increasing_3 && increasing_4 && (N>10)) {
					ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );											
					yval = spline_cubic_val ( N, t, y_original, y, ypp, &ypval, &yppval );
					geo_adapt->node[ip_01]->SetCoord(0, yval);
					geo_adapt->node[ip_01]->SetCoord(1, y_original);
					delete [] ypp; }
				
				if (!increasing_1 && !increasing_2 && !increasing_3 && increasing_4 && (N>10)) {
					ypp = spline_cubic_set ( N, t, y, ibcbeg, ybcbeg, ibcend, ybcend );											
					yval = spline_cubic_val ( N, t, -y_original, y, ypp, &ypval, &yppval );
					geo_adapt->node[ip_01]->SetCoord(0, yval);
					geo_adapt->node[ip_01]->SetCoord(1, y_original);
					delete [] ypp; }
				
				geo_adapt->bound[iMarker][iVertex] = new CLine(ip_0,ip_01,2);
				iElem_Bound_new[iMarker]++; 
				geo_adapt->bound[iMarker][iElem_Bound_new[iMarker]] = new CLine(ip_01,ip_1,2);
			}
		} */
	
	// Sumamos uno para pasar de indice a numero de elementos o nodos
	nElem_new = iElem_new + 1;
	nPoint_new = iPoint_new + 1;

	if (geometry->GetnDim() == 2) geo_adapt->SetnElem_Storage( nElem_new * 4 );
	if (geometry->GetnDim() == 3) geo_adapt->SetnElem_Storage( nElem_new * 5 );
	
	geo_adapt->SetnElem(nElem_new);
	geo_adapt->SetnPoint(nPoint_new);
	geo_adapt->SetnPointDomain(nPoint_new);
	geo_adapt->SetnDim(nDim);
		
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		nElem_Bound_new[iMarker] = iElem_Bound_new[iMarker]+1;
	
	/*--- Create the inner boundary ---*/

	if (!PlaneAdaptation) {
		geo_adapt->nElem_Bound = new unsigned long [geometry->GetnMarker()];
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			geo_adapt->SetnElem_Bound(iMarker, nElem_Bound_new[iMarker]);
		}
		
		geo_adapt->Tag_to_Marker = new string [MAX_INDEX_VALUE];	
		geo_adapt->SetnMarker(geometry->GetnMarker());
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
			geo_adapt->SetMarker_Tag(iMarker, geometry->GetMarker_Tag(iMarker));
	}
	else {
		
		geo_adapt->nElem_Bound = new unsigned long [geometry->GetnMarker()+1];
		geo_adapt->nElem_Bound_Storage = new unsigned long [geometry->GetnMarker()+1];
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			geo_adapt->SetnElem_Bound(iMarker, nElem_Bound_new[iMarker]);
			if (geometry->GetnDim() == 2) geo_adapt->SetnElem_Bound_Storage(iMarker, nElem_Bound_new[iMarker]*3);
			if (geometry->GetnDim() == 3) geo_adapt->SetnElem_Bound_Storage(iMarker, nElem_Bound_new[iMarker]*4);
		}
		
		geo_adapt->SetnElem_Bound(geometry->GetnMarker(), nElem_NearField);
		if (geometry->GetnDim() == 2) geo_adapt->SetnElem_Bound_Storage(geometry->GetnMarker(), nElem_NearField*3);
		
		geo_adapt->Tag_to_Marker = new string [MAX_INDEX_VALUE];	
		geo_adapt->SetnMarker(geometry->GetnMarker()+1);
		for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
			geo_adapt->SetMarker_Tag(iMarker, geometry->GetMarker_Tag(iMarker));
		
		geo_adapt->SetMarker_Tag(geometry->GetnMarker(), config->GetPlaneTag(0));
		config->SetMarker_All_Tag(geometry->GetnMarker(), config->GetPlaneTag(0));
		config->SetMarker_All_Plotting(geometry->GetnMarker(), YES);

		geo_adapt->bound[geometry->GetnMarker()] = new CPrimalGrid* [nElem_NearField];

		if (geometry->GetnDim() == 2)
			for (unsigned long iElem_NearField = 0; iElem_NearField < nElem_NearField; iElem_NearField++)
				geo_adapt->bound[geometry->GetnMarker()][iElem_NearField] = new CLine(Elem_NearField[iElem_NearField][0],Elem_NearField[iElem_NearField][1],2);
		
		
		if (geometry->GetnDim() == 3) {
			NearField_Trian = 0; NearField_Quad = 0;
			for (unsigned long iElem_NearField = 0; iElem_NearField < nElem_NearField; iElem_NearField++) {
				if (Elem_NearField[iElem_NearField][4] == 3) {
					NearField_Trian++;
					geo_adapt->bound[geometry->GetnMarker()][iElem_NearField] = 
					new CTriangle(Elem_NearField[iElem_NearField][0],Elem_NearField[iElem_NearField][1],
												Elem_NearField[iElem_NearField][2],3);
				}
				if (Elem_NearField[iElem_NearField][4] == 4) {
					NearField_Quad++;
					geo_adapt->bound[geometry->GetnMarker()][iElem_NearField] = 
					new CRectangle(Elem_NearField[iElem_NearField][0],Elem_NearField[iElem_NearField][1],
												 Elem_NearField[iElem_NearField][3],Elem_NearField[iElem_NearField][2],3);
				}
			}
		}
		
		if (geometry->GetnDim() == 3) geo_adapt->SetnElem_Bound_Storage(geometry->GetnMarker(), 
																																		NearField_Trian*4 + NearField_Quad*5);

		config->SetnMarker_All(geometry->GetnMarker()+1);
		
	}
	
	
	delete [] EdgeNode, DivEdge, nElem_Bound_new, iElem_Bound_new;

}

void CGridAdaptation::SetDomain_Interface(CGeometry *geometry, CPhysicalGeometry *geo_adapt, 
																					CConfig *config){
	
	unsigned long iPoint, iElem, iVertex;
	unsigned short iMarker, iNode;
	long *EdgeNode = new long[geometry->GetnEdge()];
	bool *DivEdge = new bool[geometry->GetnEdge()];
	unsigned long *nElem_Bound_new = new unsigned long [geometry->GetnMarker()];
	unsigned long *iElem_Bound_new = new unsigned long [geometry->GetnMarker()];
	unsigned long (*Elem_NearField)[5] = NULL;
	unsigned long NearField_Trian = 0, NearField_Quad = 0;
	unsigned long nElem_NearField, iElem_NearField;
	bool AddElem, AddPoint;
	double yCoord, zCoord, Dist;
	
	double eps = 10.0;
	double NearField = 200.0; //264.0; //	NearField = fabs(config->GetPosition_Plane());
	double yCut = 130.0; // 140.0;	
	double zCut = 0.0;

	nElem_NearField = geometry->GetnPoint();
	Elem_NearField = new unsigned long[nElem_NearField][5];
	for (iElem_NearField = 0; iElem_NearField < nElem_NearField; iElem_NearField++)
		Elem_NearField[iElem_NearField][4] = 0;
	nElem_NearField = 0;
	
	/*--- Create the new grid with the associated space ---*/
	geo_adapt->elem = new CPrimalGrid*[geometry->GetnElem()];
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)
			geo_adapt->elem[iElem] = new CTriangle(geometry->elem[iElem]->GetNode(0), 
																						 geometry->elem[iElem]->GetNode(1), 
																						 geometry->elem[iElem]->GetNode(2), 2);
		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)
			geo_adapt->elem[iElem] = new CTetrahedron(geometry->elem[iElem]->GetNode(0),
																								geometry->elem[iElem]->GetNode(1),
																								geometry->elem[iElem]->GetNode(2),
																								geometry->elem[iElem]->GetNode(3));
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)
			geo_adapt->elem[iElem] = new CHexahedron(geometry->elem[iElem]->GetNode(0),
																							 geometry->elem[iElem]->GetNode(1),
																							 geometry->elem[iElem]->GetNode(2),
																							 geometry->elem[iElem]->GetNode(3),
																							 geometry->elem[iElem]->GetNode(4),
																							 geometry->elem[iElem]->GetNode(5),
																							 geometry->elem[iElem]->GetNode(6),
																							 geometry->elem[iElem]->GetNode(7));
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)
			geo_adapt->elem[iElem] = new CPyramid(geometry->elem[iElem]->GetNode(0),
																								geometry->elem[iElem]->GetNode(1),
																								geometry->elem[iElem]->GetNode(2),
																								geometry->elem[iElem]->GetNode(3),
																								geometry->elem[iElem]->GetNode(4));
	}
	
	geo_adapt->node = new CPoint*[geometry->GetnPoint()];
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		if (geometry->GetnDim() == 2)
			geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
																					 geometry->node[iPoint]->GetCoord(1), iPoint, config);
		if (geometry->GetnDim() == 3)
			geo_adapt->node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
																					 geometry->node[iPoint]->GetCoord(1),
																					 geometry->node[iPoint]->GetCoord(2), iPoint, config);
	}
	
	geo_adapt->bound = new CPrimalGrid**[geometry->GetnMarker()+1];
	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++){
		geo_adapt->bound[iMarker] = new CPrimalGrid* [geometry->GetnElem_Bound(iMarker) + nElem_Bound_new[iMarker]];
		for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
				geo_adapt->bound[iMarker][iVertex] = new CLine(geometry->bound[iMarker][iVertex]->GetNode(0),
																											 geometry->bound[iMarker][iVertex]->GetNode(1),2);
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
				geo_adapt->bound[iMarker][iVertex] = new CTriangle(geometry->bound[iMarker][iVertex]->GetNode(0),
																													 geometry->bound[iMarker][iVertex]->GetNode(1), 
																													 geometry->bound[iMarker][iVertex]->GetNode(2),3);
			if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == RECTANGLE)
				geo_adapt->bound[iMarker][iVertex] = new CRectangle(geometry->bound[iMarker][iVertex]->GetNode(0),
																														geometry->bound[iMarker][iVertex]->GetNode(1),
																														geometry->bound[iMarker][iVertex]->GetNode(2), 
																														geometry->bound[iMarker][iVertex]->GetNode(3),3);
		}
	}
	
	geo_adapt->SetnElem_Storage(geometry->GetnElem_Storage());
	geo_adapt->SetnElem(geometry->GetnElem());
	geo_adapt->SetnPoint(geometry->GetnPoint());
	geo_adapt->SetnPointDomain(geometry->GetnPoint());
	geo_adapt->SetnDim(geometry->GetnDim());
	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		
		/*--- identification of the points that are in the inner domain ---*/
		AddElem = true;
		for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
			iPoint = geometry->elem[iElem]->GetNode(iNode);
			if (geometry->node[iPoint]->GetBoundary()) 
				geo_adapt->node[iPoint]->SetBoundary(true);
			
			yCoord = geometry->node[iPoint]->GetCoord(1);
			zCoord = geometry->node[iPoint]->GetCoord(2);
			Dist = sqrt (zCoord*zCoord + yCoord*yCoord);
			if (Dist > NearField + eps) AddElem = false;
		}
		
		/*--- Pick the point that compose the surface element ---*/
		if (AddElem) {
			Elem_NearField[nElem_NearField][4] = 0;
			for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
				iPoint = geometry->elem[iElem]->GetNode(iNode);				
				yCoord = geometry->node[iPoint]->GetCoord(1);
				zCoord = geometry->node[iPoint]->GetCoord(2);
				Dist = sqrt (zCoord*zCoord + yCoord*yCoord);
				AddPoint = false;
				
				if (((Dist + eps) > NearField) && ((Dist - eps) < NearField)) AddPoint = true;
				
				if ( (AddPoint) && (yCoord < yCut) && (zCoord < zCut) ) {
					geo_adapt->node[iPoint]->SetBoundary(true);
					Elem_NearField[nElem_NearField][Elem_NearField[nElem_NearField][4]] = iPoint;
					Elem_NearField[nElem_NearField][4]++;
				}
			}
		}
		if (Elem_NearField[nElem_NearField][4] == 4) nElem_NearField++;
		else Elem_NearField[nElem_NearField][4] = 0;
		
	}
	
	/*--- Create the inner boundary ---*/
	geo_adapt->SetnMarker(geometry->GetnMarker());
	geo_adapt->nElem_Bound = new unsigned long [geometry->GetnMarker()+1];
	geo_adapt->nElem_Bound_Storage = new unsigned long [geometry->GetnMarker()+1];
	geo_adapt->Tag_to_Marker = new string [MAX_INDEX_VALUE];	
	for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		geo_adapt->SetnElem_Bound(iMarker, geometry->GetnElem_Bound(iMarker));
		geo_adapt->SetnElem_Bound_Storage(iMarker, geometry->GetnElem_Bound_Storage(iMarker));
		geo_adapt->SetMarker_Tag(iMarker, geometry->GetMarker_Tag(iMarker));
	}
	
	/*--- New boundary (only for hexa) ---*/
	geo_adapt->SetnMarker(geometry->GetnMarker()+1);
	geo_adapt->SetnElem_Bound(geometry->GetnMarker(), nElem_NearField);
	geo_adapt->SetMarker_Tag(geometry->GetnMarker(), config->GetPlaneTag(0));
	geo_adapt->bound[geometry->GetnMarker()] = new CPrimalGrid* [nElem_NearField];
	
	NearField_Trian = 0; NearField_Quad = 0;
	for (iElem_NearField = 0; iElem_NearField < nElem_NearField; iElem_NearField++) {
		if (Elem_NearField[iElem_NearField][4] == 3) {
			NearField_Trian++;
			geo_adapt->bound[geometry->GetnMarker()][iElem_NearField] = 
			new CTriangle(Elem_NearField[iElem_NearField][0],Elem_NearField[iElem_NearField][1],
										Elem_NearField[iElem_NearField][2],3);
		}
		if (Elem_NearField[iElem_NearField][4] == 4) {
			NearField_Quad++;
			geo_adapt->bound[geometry->GetnMarker()][iElem_NearField] = 
			new CRectangle(Elem_NearField[iElem_NearField][0],Elem_NearField[iElem_NearField][1],
										 Elem_NearField[iElem_NearField][2],Elem_NearField[iElem_NearField][3],3);
		}
	}
	geo_adapt->SetnElem_Bound_Storage(geometry->GetnMarker(), NearField_Trian*4 + NearField_Quad*5);
	
	config->SetMarker_All_Tag(geometry->GetnMarker(), config->GetPlaneTag(0));
	config->SetMarker_All_Plotting(geometry->GetnMarker(), YES);
	config->SetnMarker_All(geometry->GetnMarker()+1);
	
	
	/*--- Move the surface points to the nearfield cylinder ---*/
	
	delete [] EdgeNode, DivEdge, nElem_Bound_new, iElem_Bound_new, (Elem_NearField)[5];
}

void CGridAdaptation::SetIndicator_Flow(CGeometry *geometry, CConfig *config, unsigned short strength){
	unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, nElem_real, max_elem_new;
	unsigned short iDim, iMarker;
	double Dual_Area, norm, Solution_Vertex, Solution_0, Solution_1, Solution_Average, 
			DualArea, Partial_Res, Grad_Val, *Normal;
	double scale_area = config->GetDualVol_Power();

	/*--- Initialization ---*/
	nElem_new = 0; nElem_real = 0;
	max_elem_new = int(0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivStrength(0);
		geometry->elem[iElem]->SetDivide(false);
		geometry->elem[iElem]->SetSemiDivide (false);
	}
	
	/*--- Compute the gradient of the first variable ---*/
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++)
			Gradient[iPoint][iDim] = 0.0;
	
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++){	
		Point_0 = geometry->edge[iEdge]->GetNode(0); Solution_0 = ConsVar_Sol[Point_0][0];
		Point_1 = geometry->edge[iEdge]->GetNode(1); Solution_1 = ConsVar_Sol[Point_1][0];
		Normal = geometry->edge[iEdge]->GetNormal();
		Solution_Average =  0.5 * ( Solution_0 + Solution_1);
		for(iDim = 0; iDim < nDim; iDim++) {
			Partial_Res = Solution_Average*Normal[iDim];
			Gradient[Point_0][iDim] = Gradient[Point_0][iDim] + Partial_Res;
			Gradient[Point_1][iDim] = Gradient[Point_1][iDim] - Partial_Res;
		}				
	}
		
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Solution_Vertex = ConsVar_Sol[Point][0];
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			for(iDim = 0; iDim < nDim; iDim++){
				Partial_Res = Solution_Vertex*Normal[iDim];
				Gradient[Point][iDim] = Gradient[Point][iDim] - Partial_Res;
			}
		}
		
	for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++){
			DualArea = geometry->node[iPoint]->GetVolume();
			Grad_Val = Gradient[iPoint][iDim]/DualArea;
			Gradient[iPoint][iDim] = Grad_Val;			
		}
	
	/*--- Compute the the adaptation index at each point ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		norm = 0.0;
		for(iDim = 0; iDim < nDim; iDim++) 
			norm += Gradient[iPoint][iDim]*Gradient[iPoint][iDim];
		norm = sqrt(norm); 
		Index[iPoint] = pow(Dual_Area,scale_area)*norm;
	}
	
	SetSensorElem(geometry, config, max_elem_new);
}


void CGridAdaptation::SetIndicator_Adj(CGeometry *geometry, CConfig *config, unsigned short strength){
	double Dual_Area;
	unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, nElem_real, max_elem_new;
	unsigned short iDim, iMarker;
	double norm, Solution_Vertex, Solution_0, Solution_1, Solution_Average, 
	DualArea, Partial_Res, Grad_Val, *Normal;
	double scale_area = config->GetDualVol_Power();
	
	// Initialization
	nElem_new = 0; nElem_real = 0; 
	max_elem_new = int(0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivStrength(0);
		geometry->elem[iElem]->SetDivide(false);
		geometry->elem[iElem]->SetSemiDivide (false);
	}
	
	
	// Compute the gradient of the density.
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++)
			Gradient[iPoint][iDim] = 0.0;
	
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++){	
		Point_0 = geometry->edge[iEdge]->GetNode(0); Solution_0 = AdjVar_Sol[Point_0][0];
		Point_1 = geometry->edge[iEdge]->GetNode(1); Solution_1 = AdjVar_Sol[Point_1][0];
		Normal = geometry->edge[iEdge]->GetNormal();
		Solution_Average =  0.5 * ( Solution_0 + Solution_1);
		for(iDim = 0; iDim < nDim; iDim++) {
			Partial_Res = Solution_Average*Normal[iDim];
			Gradient[Point_0][iDim] = Gradient[Point_0][iDim] + Partial_Res;
			Gradient[Point_1][iDim] = Gradient[Point_1][iDim] - Partial_Res;
		}				
	}
	
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Solution_Vertex = AdjVar_Sol[Point][0];
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			for(iDim = 0; iDim < nDim; iDim++){
				Partial_Res = Solution_Vertex*Normal[iDim];
				Gradient[Point][iDim] = Gradient[Point][iDim] - Partial_Res;
			}
		}
	
	for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++){
			DualArea = geometry->node[iPoint]->GetVolume();
			Grad_Val = Gradient[iPoint][iDim]/DualArea;
			Gradient[iPoint][iDim] = Grad_Val;			
		}
	
	
	// Compute the the adaptation index at each point.
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		norm = 0.0;
		for(iDim = 0; iDim < nDim; iDim++) 
			norm += Gradient[iPoint][iDim]*Gradient[iPoint][iDim];
		norm = sqrt(norm); 
		Index[iPoint] = pow(Dual_Area,scale_area)*norm;
	}
	
	SetSensorElem(geometry, config, max_elem_new);
	
}

void CGridAdaptation::SetIndicator_FlowAdj(CGeometry *geometry, CConfig *config){
	double Dual_Area;
	unsigned long Point = 0, Point_0 = 0, Point_1 = 0, iEdge, iVertex, iPoint, iElem, max_elem_new_flow, max_elem_new_adj;
	unsigned short iDim, iMarker;
	double norm, DualArea, Partial_Res, *Normal;
	double scale_area = config->GetDualVol_Power();
	
	// Initialization
	max_elem_new_flow = int(0.5*0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));
	max_elem_new_adj =  int(0.5*0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivStrength(0);
		geometry->elem[iElem]->SetDivide(false);
		geometry->elem[iElem]->SetSemiDivide (false);
	}
	
	// Compute the gradient of the first variable.
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++) {
			Gradient_Flow[iPoint][iDim] = 0.0;
			Gradient_Adj[iPoint][iDim] = 0.0;
		}

	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++){	
		Point_0 = geometry->edge[iEdge]->GetNode(0);
		Point_1 = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		for(iDim = 0; iDim < nDim; iDim++){
			Partial_Res = 0.5 * ( ConsVar_Sol[Point_0][0] + ConsVar_Sol[Point_1][0] ) * Normal[iDim];
			Gradient_Flow[Point_0][iDim] = Gradient_Flow[Point_0][iDim] + Partial_Res;
			Gradient_Flow[Point_1][iDim] = Gradient_Flow[Point_1][iDim] - Partial_Res;

			Partial_Res = 0.5 * ( AdjVar_Sol[Point_0][0] + AdjVar_Sol[Point_1][0] ) * Normal[iDim];
			Gradient_Adj[Point_0][iDim] = Gradient_Adj[Point_0][iDim] + Partial_Res;
			Gradient_Adj[Point_1][iDim] = Gradient_Adj[Point_1][iDim] - Partial_Res;			
		}				
	}
	
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			for(iDim = 0; iDim < nDim; iDim++){
				Gradient_Flow[Point][iDim] = Gradient_Flow[Point][iDim] - ConsVar_Sol[Point][0] * Normal[iDim];
				Gradient_Adj[Point][iDim] = Gradient_Adj[Point][iDim] - AdjVar_Sol[Point][0] * Normal[iDim];
			}
		}
	
	for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++){
			DualArea = geometry->node[iPoint]->GetVolume();
			Gradient_Flow[iPoint][iDim] = Gradient_Flow[iPoint][iDim]/DualArea;
			Gradient_Adj[iPoint][iDim] = Gradient_Adj[iPoint][iDim]/DualArea;
		}
	
	// Compute the the adaptation index at each point.
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area=geometry->node[iPoint]->GetVolume();
		norm = 0.0;
		for(iDim = 0; iDim < nDim; iDim++) 
			norm += Gradient_Flow[iPoint][iDim]*Gradient_Flow[iPoint][iDim];
		norm = sqrt(norm); 
		Index[iPoint] = pow(Dual_Area,scale_area)*norm;
	}
	
	
	SetSensorElem(geometry, config, max_elem_new_flow);

	// Compute the the adaptation index at each point.
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area=geometry->node[iPoint]->GetVolume();
		norm = 0.0;
		for(iDim = 0; iDim < nDim; iDim++) 
			norm += Gradient_Adj[iPoint][iDim]*Gradient_Adj[iPoint][iDim];
		norm = sqrt(norm); 
		Index[iPoint] = pow(Dual_Area,scale_area)*norm;
	}
	
	SetSensorElem(geometry, config, max_elem_new_adj);

}

void CGridAdaptation::SetIndicator_Robust(CGeometry *geometry, CConfig *config){
	unsigned long iPoint, iElem, max_elem_new_flow, max_elem_new_adj;
	unsigned short iVar;
	double Dual_Area;
	double scale_area = config->GetDualVol_Power();
	
	// Inicializa la malla para la adaptacion
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
		geometry->elem[iElem]->SetSemiDivide (false);
		geometry->elem[iElem]->SetDivStrength(0);
	}
	
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++){
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += ConsVar_Res[iPoint][iVar]*ConsVar_Res[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area,scale_area)*sqrt(Index[iPoint]);
	}
	
	max_elem_new_flow = int(0.5*0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));
	SetSensorElem(geometry, config, max_elem_new_flow);

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++){
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += AdjVar_Res[iPoint][iVar]*AdjVar_Res[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area,scale_area)*sqrt(Index[iPoint]);
	}
	
	max_elem_new_adj = int(0.5*0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));
	SetSensorElem(geometry, config, max_elem_new_adj);

}

void CGridAdaptation::SetIndicator_Computable(CGeometry *geometry, CConfig *config){
	unsigned long iPoint, iElem, max_elem_new;
	unsigned short iVar;
	double max_indicator = 0.0, Dual_Area; 
	double scale_area = config->GetDualVol_Power();
	
	max_elem_new = int(0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
		geometry->elem[iElem]->SetSemiDivide (false);
		geometry->elem[iElem]->SetDivStrength(0);
	}
		
	max_indicator = 0.0; 		
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += ConsVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar]*ConsVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area,scale_area)*sqrt(Index[iPoint]);
	}
	
	SetSensorElem(geometry, config, max_elem_new);

}

void CGridAdaptation::SetIndicator_Computable_Robust(CGeometry *geometry, CConfig *config){
	unsigned long iPoint, iElem, max_elem_new;
	unsigned short iVar;
	double max_indicator, Dual_Area ;
	double scale_area = config->GetDualVol_Power();

	/*--- Initializate the numerical grid for the adaptation ---*/
	max_elem_new = int(0.01*config->GetNew_Elem_Adapt()*double(geometry->GetnElem()));	
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		geometry->elem[iElem]->SetDivide (false);
		geometry->elem[iElem]->SetSemiDivide (false);
		geometry->elem[iElem]->SetDivStrength(0);
	}
	
	
	max_indicator = 0.0; 		
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		Dual_Area = geometry->node[iPoint]->GetVolume();
		Index[iPoint] = 0.0;
		for (iVar = 0; iVar < nVar; iVar++)
			Index[iPoint] += LinVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar]*LinVar_Res[iPoint][iVar]*AdjVar_Sol[iPoint][iVar];
		
		Index[iPoint] = pow(Dual_Area,scale_area)*sqrt(Index[iPoint]);
	}

	SetSensorElem(geometry, config, max_elem_new);
	
}

void CGridAdaptation::SetReStart_FlowSolution(CGeometry *geometry, CConfig *config, 
								 string mesh_flowfilename){
	
	unsigned long iPoint;
	unsigned short iVar;
		
	//	Rearranque del problema directo
	char *cstr = new char [mesh_flowfilename.size()+1];
	strcpy (cstr, mesh_flowfilename.c_str());
	
	ofstream restart_flowfile;
	restart_flowfile.open(cstr, ios::out);
	restart_flowfile.precision(15);
	for(iPoint = 0; iPoint < nPoint_new; iPoint++){
		restart_flowfile << iPoint <<"\t";
		for(iVar = 0; iVar < nVar; iVar++) {
			restart_flowfile << scientific << ConsVar_Adapt[iPoint][iVar]<<"\t";
		}
		restart_flowfile << endl;
	}
	restart_flowfile.close();
}

void CGridAdaptation::SetReStart_AdjSolution(CGeometry *geometry, CConfig *config, 
											 string mesh_adjfilename){
	char cstr[200], buffer[50];
	string copy;
	
	copy.assign(mesh_adjfilename);
	copy.erase (copy.end()-4, copy.end());
	strcpy (cstr, copy.c_str()); 
	if (config->GetKind_ObjFunc() == DRAG_COEFFICIENT) sprintf (buffer, "_cd.dat"); 
	if (config->GetKind_ObjFunc() == LIFT_COEFFICIENT) sprintf (buffer, "_cl.dat");
	if (config->GetKind_ObjFunc() == SIDEFORCE_COEFFICIENT) sprintf (buffer, "_csf.dat"); 
	if (config->GetKind_ObjFunc() == PRESSURE_COEFFICIENT) sprintf (buffer, "_cp.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_X_COEFFICIENT) sprintf (buffer, "_cmx.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Y_COEFFICIENT) sprintf (buffer, "_cmy.dat"); 
	if (config->GetKind_ObjFunc() == MOMENT_Z_COEFFICIENT) sprintf (buffer, "_cmz.dat"); 
	if (config->GetKind_ObjFunc() == EFFICIENCY) sprintf (buffer, "_eff.dat"); 
	if (config->GetKind_ObjFunc() == ELECTRIC_CHARGE) sprintf (buffer, "_cc.dat"); 
  if (config->GetKind_ObjFunc() == FORCE_X_COEFFICIENT) sprintf (buffer, "_cfx.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Y_COEFFICIENT) sprintf (buffer, "_cfy.dat"); 
	if (config->GetKind_ObjFunc() == FORCE_Z_COEFFICIENT) sprintf (buffer, "_cfz.dat");
//  if (config->GetKind_ObjFunc() == ROTOR_EFFICIENCY) sprintf (buffer, "_rot_eff.dat");
	strcat(cstr, buffer);
	
	ofstream restart_adjfile;
	restart_adjfile.open(cstr, ios::out);
	restart_adjfile.precision(15);
	for(unsigned long iPoint = 0; iPoint < nPoint_new; iPoint++){
		restart_adjfile << iPoint <<"\t";
		for(unsigned short iVar = 0; iVar < nVar; iVar++)
			restart_adjfile << scientific << AdjVar_Adapt[iPoint][iVar]<<"\t";
		restart_adjfile << endl;
	}
	restart_adjfile.close();
}

void CGridAdaptation::SetReStart_LinSolution(CGeometry *geometry, CConfig *config, 
											 string mesh_linfilename){
	
	unsigned long iPoint;
	unsigned short iVar;
	
	//	Rearranque del problema adjunto
	char *cstr_ = new char [mesh_linfilename.size()+1];
	strcpy (cstr_, mesh_linfilename.c_str());
	
	ofstream restart_linfile;
	restart_linfile.open(cstr_, ios::out);
	restart_linfile.precision(15);
	for(iPoint = 0; iPoint < nPoint_new; iPoint++){
		restart_linfile << iPoint <<"\t";
		for(iVar = 0; iVar < nVar; iVar++)
			restart_linfile << scientific << LinVar_Adapt[iPoint][iVar]<<"\t";
		restart_linfile << endl;
	}
	restart_linfile.close();
}

void CGridAdaptation::SetSensorElem(CGeometry *geometry, CConfig *config, unsigned long max_elem) {
	double Max_Sensor, threshold;
	double *Sensor = new double[geometry->GetnElem()];
	unsigned long ip_0, ip_1, ip_2, ip_3, iElem, nElem_real;
	
	/*--- Compute the the adaptation index at each element ---*/
	Max_Sensor = 0.0;
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		ip_0 = geometry->elem[iElem]->GetNode(0);
		ip_1 = geometry->elem[iElem]->GetNode(1);
		ip_2 = geometry->elem[iElem]->GetNode(2);
		Sensor[iElem] = (Index[ip_0]+Index[ip_1]+Index[ip_2])/3.0;
		if ((geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) ||
			(geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)) {
			ip_3 = geometry->elem[iElem]->GetNode(2);
			Sensor[iElem] = (Index[ip_0]+Index[ip_1]+Index[ip_2]+Index[ip_3])/4.0;
		}
		Max_Sensor = max(Max_Sensor,Sensor[iElem]);
	}
	
	/*--- Adimensionalization of the adaptation sensor ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
		Sensor[iElem] = Sensor[iElem]/Max_Sensor;
	}
	
	/*--- Selection of the elements to be adapted ---*/
	threshold = 0.999;
	nElem_real = 0;
	while (nElem_real <= max_elem) {
		for (iElem = 0; iElem < geometry->GetnElem(); iElem ++)
			if ( Sensor[iElem] >= threshold && !geometry->elem[iElem]->GetDivide() ) {
				if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE) nElem_real = nElem_real + 3;	
				if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE) nElem_real = nElem_real + 3;	
				if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem_real = nElem_real + 7;
				geometry->elem[iElem]->SetDivStrength(1);
				geometry->elem[iElem]->SetDivide(true);
				if (nElem_real >= max_elem) break;
			}	
		threshold = threshold - 0.001;
	}
	
	cout << "Number of elements to adapt: " << nElem_real << endl;
	delete [] Sensor;
}

void CGridAdaptation::WriteAdaptSensor(CGeometry *geometry, char mesh_filename[200])
{
	unsigned long iPoint;
	
	ofstream para_file;
	para_file.precision(15);
	para_file.open(mesh_filename, ios::out | ios::app);
		
	para_file << "POINT_DATA " << geometry->GetnPoint() << endl;
	para_file << "SCALARS Adapt_Sensor float 1" << endl;
	para_file << "LOOKUP_TABLE default" << endl;
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		para_file << scientific << Index[iPoint] << endl;
}
