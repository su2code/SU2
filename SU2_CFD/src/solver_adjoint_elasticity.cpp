/*!
 * \file solver_adjoint_elasticity.cpp
 * \brief Main subroutines for solving adjoint FEM elasticity problems.
 * \author R. Sanchez
 * \version 4.1.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/solver_structure.hpp"

/*------------------------------ RUBEN: TEMPORARY FUNCTIONS - ADJOINT STRUCTURE ------------------------------*/

void CFEM_ElasticitySolver::Run_Structural_Adjoint(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh,
		unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output){

	unsigned short iMarker = 1;

	cout << endl;
	cout << "--------------- I'm going to run the structural adjoint --------------------" << endl;

	Initialize_Structural_Adj(geometry, solver_container, numerics, config);

	Compute_RefGeom_Sensitivity(geometry, solver_container, numerics, config);

	BC_Clamped_Adj(geometry, solver_container, numerics[FEA_TERM], config, iMarker);

	Solve_System_Adj(geometry, solver_container, config);

	BC_Clamped_Post_Adj(geometry, solver_container, numerics[FEA_TERM], config, iMarker);

	Compute_DE_Sensitivity(geometry, solver_container, numerics, config);

	Compute_RefGeom_Gradient(geometry, solver_container, numerics, config);


}

void CFEM_ElasticitySolver::Initialize_Structural_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	cout << "STEP 1: Initialize the structural adjoint." << endl;

	unsigned long iPoint;
	unsigned short iVar;

	/*--- Initialize the gradient variable, and the residual and solution vectors ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++){
    	for (iVar = 0; iVar < nVar; iVar++){
        	node[iPoint]->SetGradient_Adj(iVar, 0.0);
    	}
		LinSysRes_Adj.SetBlock_Zero(iPoint);
		LinSysSol_Adj.SetBlock_Zero(iPoint);
    }


}

void CFEM_ElasticitySolver::Compute_RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	cout << "STEP 2: Compute the sensitivity of the objective function (I) respect to the state variables (x)." << endl;

	unsigned short iVar;
	unsigned long iPoint;
	su2double *reference_geometry, *current_solution;
	unsigned long realIndex;

    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	reference_geometry = node[iPoint]->GetReference_Geometry();
    	current_solution = node[iPoint]->GetSolution();

    	realIndex = geometry->node[iPoint]->GetGlobalIndex();

    	for (iVar = 0; iVar < nVar; iVar++){

    		Solution[iVar] = 2*(current_solution[iVar]-reference_geometry[iVar]);

    	}

//    	cout << "Node " << realIndex << ". dI/dx = (" << Solution[0] << "," << Solution[1] << ")." << endl;

    	LinSysRes_Adj.AddBlock(iPoint, Solution);

    }

}

void CFEM_ElasticitySolver::BC_Clamped_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker){

	cout << "STEP 3: Apply BC to the right hand side of the equations." << endl;

}

void CFEM_ElasticitySolver::Solve_System_Adj(CGeometry *geometry, CSolver **solver_container, CConfig *config){

	cout << "STEP 4: Solve the linear system for the adjoint variable." << endl;

	unsigned long IterLinSol = 0, iPoint, total_index;
	unsigned short iVar;

	/*--- Initialize residual and solution at the ghost points ---*/

	for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {

		for (iVar = 0; iVar < nVar; iVar++) {
		  total_index = iPoint*nVar + iVar;
		  LinSysRes[total_index] = 0.0;
		  LinSysSol[total_index] = 0.0;
		}

	 }


	cout << "Solving adjoint system..." << endl;

	CSysSolve femSystem;
	IterLinSol = femSystem.Solve(Jacobian, LinSysRes_Adj, LinSysSol_Adj, geometry, config);

//	cout << "Adjoint system solved with " << IterLinSol << " iterations." << endl;

	unsigned long realIndex;

    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	realIndex = geometry->node[iPoint]->GetGlobalIndex();

    	for (iVar = 0; iVar < nVar; iVar++){

    		Solution[iVar] = LinSysSol_Adj.GetBlock(iPoint, iVar);

        	node[iPoint]->SetSolution_Adj(iVar, Solution[iVar]);

    	}

//    	cout << "Node " << realIndex << ". PHI = (" << Solution[0] << "," << Solution[1] << ")." << endl;

    }

	/*--- The the number of iterations of the linear solver ---*/

	SetIterLinSolver(IterLinSol);

}

void CFEM_ElasticitySolver::BC_Clamped_Post_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker){

	cout << "STEP 5: Reinforce BC to the sensitivities." << endl;

}


void CFEM_ElasticitySolver::Compute_DE_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	cout << "STEP 6: compute the sensitivities of the structural equations respect to the design variables E." << endl;


	unsigned long iElem, iVar;
	unsigned short iNode, iDim, nNodes;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord, val_Sol;
	int EL_KIND;

	su2double *Ta = NULL;
	unsigned short NelNodes;

	unsigned long testID;

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[DE_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[DE_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		numerics[DE_ADJ]->Compute_NodalStress_Term(element_container[DE_TERM][EL_KIND], config);

		NelNodes = element_container[DE_TERM][EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			testID = geometry->node[indexNode[iNode]]->GetGlobalIndex();

			Ta = element_container[DE_TERM][EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++){
				node[indexNode[iNode]]->AddGradient_Adj(iVar,Ta[iVar]);
			}
//			cout << "For ID " << testID << " we add Ta = (" << Ta[0] << "," << Ta[1] << ")." << endl;

		}

	}


}


void CFEM_ElasticitySolver::Compute_RefGeom_Gradient(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	cout << "STEP 7: output the interest function I and the gradient dI/dv (v: design variables)." << endl;

	su2double *ref_geom, *curr_sol;
	su2double adj_var, sens_var;
	su2double objective_function, gradient_value;
	//string adj_name, sens_name, res_name, file_type;
	char adj_name[200], sens_name[200], res_name[200], file_type[200];

	strcpy (adj_name,"Adjoint_E");
	strcpy (sens_name,"Sensitivity_E");
	strcpy (res_name,"Results_E");
	strcpy (file_type,".txt");


	char buffer[50];

	int e_field_int = (int)config->Get_Electric_Field_Mod(0);

    if ((e_field_int >= 0)    && (e_field_int < 10))    SPRINTF (buffer, "_0000%d.dat", e_field_int);
    if ((e_field_int >= 10)   && (e_field_int < 100))   SPRINTF (buffer, "_000%d.dat",  e_field_int);
    if ((e_field_int >= 100)  && (e_field_int < 1000))  SPRINTF (buffer, "_00%d.dat",   e_field_int);
    if ((e_field_int >= 1000) && (e_field_int < 10000)) SPRINTF (buffer, "_0%d.dat",    e_field_int);

    string key_append = string(buffer);
    strcat(adj_name,buffer);
    strcat(sens_name,buffer);

    strcat(adj_name,file_type);
    strcat(sens_name,file_type);
    strcat(res_name,file_type);

	unsigned long iPoint, realIndex;
	unsigned short iVar;

	objective_function = 0.0;
	gradient_value = 0.0;

	ofstream myfile_adj, myfile_sens, myfile_res;
	myfile_adj.open (adj_name, ios::out);
	myfile_sens.open (sens_name, ios::out);
	myfile_res.open (res_name, ios::app);

	cout << "Point " << "adj_var   " << " " << "sens_var  " << endl;

    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	ref_geom = node[iPoint]->GetReference_Geometry();
    	curr_sol = node[iPoint]->GetSolution();

    	realIndex = geometry->node[iPoint]->GetGlobalIndex();

    	for (iVar = 0; iVar < nVar; iVar++){

    		objective_function += (curr_sol[iVar]-ref_geom[iVar])*(curr_sol[iVar]-ref_geom[iVar]);

    		adj_var = node[iPoint]->GetSolution_Adj(iVar);
    		sens_var = node[iPoint]->GetGradient_Adj(iVar);
    		myfile_adj << realIndex<< " " << adj_var << endl;
    		myfile_sens << realIndex<< " " << sens_var << endl;
    		gradient_value -= adj_var*sens_var;

    	}

    }

    myfile_res << e_field_int << " " << objective_function << " "<< gradient_value << endl;
	cout << "------------------------------- DONE ---------------------------------------" << endl;

	myfile_adj.close();
	myfile_sens.close();
	myfile_res.close();

}

/*------------------------------ RUBEN: END OF TEMPORARY FUNCTIONS - ADJOINT STRUCTURE ------------------------------*/





