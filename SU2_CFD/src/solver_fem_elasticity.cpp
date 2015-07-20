/*!
 * \file solution_template.cpp
 * \brief Main subroutines for solving direct FEM elasticity problems.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
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

CFEM_ElasticitySolver::CFEM_ElasticitySolver(void) : CSolver() {

	nElement = 0;
	nDim = 0;
	nMarker = 0;

	nPoint = 0;
	nPointDomain = 0;

	element_container = NULL;
	node = NULL;

}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long iPoint, iElem;
	unsigned short iVar, jVar, iDim, NodesElement = 0, nKindElements;

	nElement      = geometry->GetnElem();
	nDim          = geometry->GetnDim();
	nMarker       = geometry->GetnMarker();

	nPoint        = geometry->GetnPoint();
	nPointDomain  = geometry->GetnPointDomain();

	nKindElements = 2;

	element_container = new CElement*[nKindElements];
	node          	  = new CVariable*[nPoint];

	GradN_X = new double [nDim];
	GradN_x = new double [nDim];

	nVar = nDim;

	/*--- Define some auxiliary vectors related to the residual ---*/

	Residual = new double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max = new double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
	Point_Max_Coord = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Point_Max_Coord[iVar] = new double[nDim];
		for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
	}

	/*--- Define some auxiliary vectors related to the solution ---*/

	Solution   = new double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;

	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	/*--- Check for a restart, initialize from zero otherwise ---*/

	if (!restart) {
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
			node[iPoint] = new CFEM_LElasVariable(Solution, nDim, nVar, config);
		}
	}
	else {

		/* The restart from a file needs to be implemented */

	}

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_ij[iVar][jVar] = 0.0;
			}
	}

	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

	/*--- Here is where we assign the kind of each element ---*/

	for (iElem = 0; iElem < nElement; iElem++){

		/*--- As of now, only QUAD4 elements ---*/
		element_container[iElem] = new CQUAD4(nDim, iElem, config);

	}

	if (nDim == 2){
		element_container[EL_TRIA] = new CTRIA1(nDim, iElem, config);
		element_container[EL_QUAD] = new CQUAD4(nDim, iElem, config);
	}

}

CFEM_ElasticitySolver::~CFEM_ElasticitySolver(void) { }

void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;

	double *Kab;
	unsigned short NelNodes, jNode;

	cout << nElement << endl;

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        nNodes = 6;
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_QUAD]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_QUAD]->SetCurr_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics[VISC_TERM]->Compute_Tangent_Matrix(element_container[EL_QUAD]);

		NelNodes = element_container[EL_QUAD]->GetnNodes();
		for (iNode = 0; iNode < NelNodes; iNode++){
			for (jNode = 0; jNode < NelNodes; jNode++){
				Kab = element_container[EL_QUAD]->Get_Kab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
					}
				}

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

			}

		}

	}

	double checkJacobian;

	ofstream myfile;
	myfile.open ("newSolver.txt");

	for (iNode = 0; iNode < nPoint; iNode++){
		for (jNode = 0; jNode < nPoint; jNode++){
			myfile << "Node " << iNode << " " << jNode << endl;
			for (iVar = 0; iVar < nVar; iVar++){
				for (jVar = 0; jVar < nVar; jVar++){
					checkJacobian = Jacobian.GetBlock(iNode, jNode, iVar, jVar);
					myfile << checkJacobian << " " ;
				}
				myfile << endl;
			}
		}
	}
	myfile.close();

}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }
