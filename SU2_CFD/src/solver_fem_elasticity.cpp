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

	GradN_X = NULL;
	GradN_x = NULL;

	Jacobian_s_ij = NULL;
	Jacobian_k_ij = NULL;

	Res_Stress_i = NULL;

}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long iPoint, iElem = 0;
	unsigned short iVar, jVar, iDim, NodesElement = 0, nKindElements;

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

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



	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (stress contribution) ---*/

	Jacobian_s_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_s_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_s_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (incompressibility term) ---*/

	if (incompressible){
		Jacobian_k_ij = new double*[nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_k_ij[iVar] = new double [nVar];
				for (jVar = 0; jVar < nVar; jVar++) {
					Jacobian_k_ij[iVar][jVar] = 0.0;
				}
		}
	}
	else {
		Jacobian_k_ij = NULL;
	}

	/*--- Stress contribution to the node i ---*/
	Res_Stress_i = new double[nVar];


	/*--- Initialization of matrix structures ---*/
	if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;
	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
	MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

	/*--- Initialization of linear solver structures ---*/
	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

	/*--- Here is where we assign the kind of each element ---*/

	if (nDim == 2){
		if (incompressible){
			element_container[EL_TRIA] = new CTRIA1(nDim, config);
			element_container[EL_QUAD] = new CQUAD4P1(nDim, config);
		}
		else{
			element_container[EL_TRIA] = new CTRIA1(nDim, config);
			element_container[EL_QUAD] = new CQUAD4(nDim, config);
		}
	}
	else if (nDim == 3){
		if (incompressible){
			element_container[EL_TETRA] = new CTETRA1(nDim, config);
			element_container[EL_HEXA] = new CHEXA8P1(nDim, config);
		}
		else{
			element_container[EL_TETRA] = new CTETRA1(nDim, config);
			element_container[EL_HEXA] = new CHEXA8(nDim, config);
		}
	}


}

CFEM_ElasticitySolver::~CFEM_ElasticitySolver(void) {

	unsigned short iVar, nKindElements = 2;
	unsigned long iPoint;

	for (iPoint = 0; iPoint < nPoint; iPoint++){
		delete [] node[iPoint];
	}

	for (iVar = 0; iVar < nKindElements; iVar++){
		delete [] element_container[iVar];
	}

	for (iVar = 0; iVar < nVar; iVar++){
		delete [] Jacobian_s_ij[iVar];
		delete [] Jacobian_ij[iVar];
		delete [] Point_Max_Coord[iVar];
	}

	delete [] element_container;
	delete [] node;
	delete [] Jacobian_s_ij;
	delete [] Jacobian_ij;
	delete [] Res_Stress_i;
	delete [] Solution;
	delete [] GradN_X;
	delete [] GradN_x;

	delete [] Residual;
	delete [] Residual_RMS;
	delete [] Residual_Max;
	delete [] Point_Max;
	delete [] Point_Max_Coord;

}

void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {

//	unsigned long iPoint, iElem, iVar, jVar;
//	unsigned short iNode, iGauss, iDim;
//	unsigned short nNodes, nGauss;
//	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
//	double val_Coord, val_Sol;
//	int EL_KIND;
//
//	double Ks_ab;
//	double *Kab = NULL;
//	double *Kk_ab = NULL;
//	double *Ta = NULL;
//	unsigned short NelNodes, jNode;
//
//	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);
//
//	/*--- Loops over all the elements ---*/
//
//	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
//
//		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
//		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}
//
//		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
//		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
//		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
//		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}
//
//		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/
//
//		for (iNode = 0; iNode < nNodes; iNode++) {
//		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
////		  cout << "Elem: " << iElem << " iNode (renum):" << indexNode[iNode] << endl;
//		  for (iDim = 0; iDim < nDim; iDim++) {
//			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
////			  cout << "Coord[" << iDim << "]: " << val_Coord << endl;
//			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
//			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
//			  element_container[EL_KIND]->SetCurr_Coord(val_Coord, iNode, iDim);
//		  }
//		}
//
//		numerics[VISC_TERM]->Compute_Tangent_Matrix(element_container[EL_KIND]);
//
//		if (incompressible) numerics[VISC_TERM]->Compute_MeanDilatation_Term(element_container[EL_KIND]);
//
//		/*--- I can do it separately, or together within Compute_Tangent_Matrix (more efficient) ---*/
//		//numerics[VISC_TERM]->Compute_NodalStress_Term(element_container[EL_KIND]);
//
//		NelNodes = element_container[EL_KIND]->GetnNodes();
//
//		for (iNode = 0; iNode < NelNodes; iNode++){
//
//			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
//			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];
//
//			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);
//
//			for (jNode = 0; jNode < NelNodes; jNode++){
//
//				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);
//				Ks_ab = element_container[EL_KIND]->Get_Ks_ab(iNode,jNode);
//				if (incompressible) Kk_ab = element_container[EL_KIND]->Get_Kk_ab(iNode,jNode);
//
//				for (iVar = 0; iVar < nVar; iVar++){
//					Jacobian_s_ij[iVar][iVar] = Ks_ab;
//					for (jVar = 0; jVar < nVar; jVar++){
//						Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
//						if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
//					}
//				}
//
//				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
//
//				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
//
//				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
//
//			}
//
//		}
//
//	}
//
//
//	double checkJacobian;
//
//	ofstream myfile;
//	myfile.open ("newSolver.txt");
//
//	for (iNode = 0; iNode < nPoint; iNode++){
//		for (jNode = 0; jNode < nPoint; jNode++){
//			myfile << "Node " << iNode << " " << jNode << endl;
//			for (iVar = 0; iVar < nVar; iVar++){
//				for (jVar = 0; jVar < nVar; jVar++){
//					checkJacobian = Jacobian.GetBlock(iNode, jNode, iVar, jVar);
//					myfile << checkJacobian << " " ;
//				}
//				myfile << endl;
//			}
//		}
//	}
//	myfile.close();

}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

void CFEM_ElasticitySolver::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Ks_ab;
	double *Kab = NULL;
	double *Kk_ab = NULL;
	double *Ta = NULL;
	unsigned short NelNodes, jNode;

	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
//		  cout << "Elem: " << iElem << " iNode (renum):" << indexNode[iNode] << endl;
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
//			  cout << "Coord[" << iDim << "]: " << val_Coord << endl;
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		if (incompressible) numerics->Compute_MeanDilatation_Term(element_container[EL_KIND]);

		/*--- I can do it separately, or together within Compute_Tangent_Matrix (more efficient) ---*/
//		//numerics[VISC_TERM]->Compute_NodalStress_Term(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);
				Ks_ab = element_container[EL_KIND]->Get_Ks_ab(iNode,jNode);
				if (incompressible) Kk_ab = element_container[EL_KIND]->Get_Kk_ab(iNode,jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					Jacobian_s_ij[iVar][iVar] = Ks_ab;
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
						if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
					}
				}

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);

				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);

			}

		}

	}

	double checkJacobian;

	ofstream myfile;
	myfile.open ("newSolver_HOLA.txt");

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

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Ks_ab;
	double *Kab = NULL;
	double *Kk_ab = NULL;
	double *Ta = NULL;
	unsigned short NelNodes, jNode;

	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Loops over all the elements ---*/

	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

		if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     {nNodes = 3; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    {nNodes = 4; EL_KIND = EL_QUAD;}

		if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  {nNodes = 4; EL_KIND = EL_TETRA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      {nNodes = 5; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == PRISM)        {nNodes = 6; EL_KIND = EL_TRIA;}
		if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   {nNodes = 8; EL_KIND = EL_HEXA;}

		/*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

		for (iNode = 0; iNode < nNodes; iNode++) {
		  indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);
		  for (iDim = 0; iDim < nDim; iDim++) {
			  val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics->Compute_Mass_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);
				Ks_ab = element_container[EL_KIND]->Get_Ks_ab(iNode,jNode);
				if (incompressible) Kk_ab = element_container[EL_KIND]->Get_Kk_ab(iNode,jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					Jacobian_s_ij[iVar][iVar] = Ks_ab;
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
						if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
					}
				}

				MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);

				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);

			}

		}

	}

	double checkJacobian;

	ofstream myfile;
	myfile.open ("massMatrix.txt");

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

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) { }

void CFEM_ElasticitySolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) { }

void CFEM_ElasticitySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
		unsigned short iMesh) { }

