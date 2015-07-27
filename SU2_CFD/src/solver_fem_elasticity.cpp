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

	Jacobian_c_ij = NULL;
	Jacobian_s_ij = NULL;
	Jacobian_k_ij = NULL;

	MassMatrix_ij = NULL;

	mZeros_Aux = NULL;
	mId_Aux = NULL;

	Res_Stress_i = NULL;
	Res_Ext_Surf = NULL;
	Res_Time_Cont = NULL;

}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long iPoint, iElem = 0;
	unsigned short iVar, jVar, iDim, jDim, NodesElement = 0, nKindElements;

	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	double E = config->GetElasticyMod();

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
			node[iPoint] = new CFEM_ElasVariable(Solution, nDim, nVar, config);
		}
	}
	else {

		/* The restart from a file needs to be implemented */

	}



	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Term ij of the Mass Matrix ---*/

	MassMatrix_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		MassMatrix_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				MassMatrix_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (constitutive contribution) ---*/

	Jacobian_c_ij = new double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_c_ij[iVar] = new double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_c_ij[iVar][jVar] = 0.0;
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

	/*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
	Res_Ext_Surf = new double[nVar];

	/*--- Time integration contribution to the residual ---*/
	if (dynamic) {
		Res_Time_Cont = new double [nVar];
	}
	else {
		Res_Time_Cont = NULL;
	}

	/*--- Matrices to impose clamped boundary conditions (TODO: Initialize them conditionally). ---*/

	mZeros_Aux = new double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mZeros_Aux[iDim] = new double[nDim];

	mId_Aux = new double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mId_Aux[iDim] = new double[nDim];

	for(iDim = 0; iDim < nDim; iDim++){
		for (jDim = 0; jDim < nDim; jDim++){
			mZeros_Aux[iDim][jDim] = 0.0;
			mId_Aux[iDim][jDim] = 0.0;
		}
		mId_Aux[iDim][iDim] = E;		// TODO: This works for clamped boundary conditions...
	}


	/*--- Initialization of matrix structures ---*/
	if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;

	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

//	if (nonlinear_analysis) StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

	if (dynamic) {
		MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
		TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
		TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
	}


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
		delete [] Jacobian_c_ij[iVar];
		if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij[iVar];
		delete [] Point_Max_Coord[iVar];
		delete [] mZeros_Aux[iVar];
		delete [] mId_Aux[iVar];
	}

	delete [] element_container;
	delete [] node;
	delete [] Jacobian_s_ij;
	delete [] Jacobian_ij;
	delete [] Jacobian_c_ij;
	if (Jacobian_k_ij != NULL) delete[] Jacobian_k_ij;
	delete [] Res_Stress_i;
	delete [] Res_Ext_Surf;
	if (Res_Time_Cont != NULL) delete[] Res_Time_Cont;
	delete [] Solution;
	delete [] GradN_X;
	delete [] GradN_x;

	delete [] Residual;
	delete [] Residual_RMS;
	delete [] Residual_Max;
	delete [] Point_Max;
	delete [] Point_Max_Coord;

	delete [] mZeros_Aux;
	delete [] mId_Aux;

}

void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {


	unsigned long iPoint;
	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool first_iter = (Iteration == 0);													// Checks if it is the first iteration
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method


	/*--- Set vector entries to zero ---*/

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
		LinSysAux.SetBlock_Zero(iPoint);
		LinSysRes.SetBlock_Zero(iPoint);
		LinSysSol.SetBlock_Zero(iPoint);
	}

	/*--- Set matrix entries to zero ---*/

	/*
	 * If the problem is linear, we only need one Jacobian matrix in the problem, because
	 * it is going to be constant along the calculations. Therefore, we only initialize
	 * the Jacobian matrix once, at the beginning of the simulation.
	 *
	 * We don't need first_iter, because there is only one iteration per time step in linear analysis.
	 */
	if ((initial_calc) && (linear_analysis)){
		Jacobian.SetValZero();
	}

	/*
	 * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
	 * both for linear and nonlinear analysis. Only initialized once, at the first time step.
	 *
	 * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
	 */
	if ((dynamic) && (initial_calc) && (first_iter)) {
		MassMatrix.SetValZero();
	}

	/*
	 * If the problem is nonlinear, we need to initialize the Jacobian and the stiffness matrix at least at the beginning
	 * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
	 * iteration.
	 */

	if ((nonlinear_analysis) && ((newton_raphson) || (first_iter)))	{
		Jacobian.SetValZero();
//		StiffMatrix.SetValZero();
	}

	/*
	 * Some external forces may be considered constant over the time step.
	 */
	if (first_iter)	{
		for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Clear_SurfaceLoad_Res();
	}

}

void CFEM_ElasticitySolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned long Iteration) { }

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

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
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[EL_KIND]->Get_Kab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
					}
				}

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);

			}

		}

	}


}

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
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

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
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		if (incompressible) numerics->Compute_MeanDilatation_Term(element_container[EL_KIND]);

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
						Jacobian_c_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
						if (incompressible) Jacobian_k_ij[iVar][jVar] = Kk_ab[iVar*nVar+jVar];
					}
				}

				/*--- If the problem is dynamic, the jacobian of the problem will be a combination of the ---*/
				/*--- StiffMatrix and the Mass Matrix; otherwise, the Jacobian can be computed straight away ---*/
//				if (dynamic){
//					StiffMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
//					StiffMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
//					if (incompressible) StiffMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
//				}
//				else{
				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);
//				}

			}

		}

	}

}

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	double val_Coord, val_Sol;
	int EL_KIND;

	double Mab;
	unsigned short NelNodes, jNode;

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

			for (jNode = 0; jNode < NelNodes; jNode++){

				Mab = element_container[EL_KIND]->Get_Mab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					MassMatrix_ij[iVar][iVar] = Mab;
				}

				MassMatrix.AddBlock(indexNode[iNode], indexNode[jNode], MassMatrix_ij);

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
					checkJacobian = MassMatrix.GetBlock(iNode, jNode, iVar, jVar);
					myfile << checkJacobian << " " ;
				}
				myfile << endl;
			}
		}
	}
	myfile.close();

}

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {


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
			  val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
			  element_container[EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
			  element_container[EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
		  }
		}

		numerics->Compute_NodalStress_Term(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

		}

	}

	for (iDim = 0; iDim < nDim; iDim++) {
		val_Coord = geometry->node[0]->GetCoord(iDim);
		val_Sol = node[0]->GetSolution(iDim) + val_Coord;
		cout << val_Coord << " " << val_Sol << "    ";
	}
	cout << endl;

	double checkResidual;

	ofstream myfile;
	myfile.open ("residual_2.txt");

	for (iNode = 0; iNode < nPoint; iNode++){
			myfile << "Node " << iNode << endl;
			for (iVar = 0; iVar < nVar; iVar++){
					checkResidual = LinSysRes.GetBlock(iNode, iVar); myfile << checkResidual << endl;

			}
	}
	myfile.close();


}

void CFEM_ElasticitySolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CFEM_ElasticitySolver::Compute_IntegrationConstants(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	double Delta_t= config->GetDelta_DynTime();
	double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();

	/*--- Integration constants for Newmark scheme ---*/

	a_dt[0]= 1 / (alpha*pow(Delta_t,2.0));
	a_dt[1]= delta / (alpha*Delta_t);
	a_dt[2]= 1 / (alpha*Delta_t);
	a_dt[3]= 1 /(2*alpha) - 1;
	a_dt[4]= delta/alpha - 1;
	a_dt[5]= (Delta_t/2) * (delta/alpha - 2);
	a_dt[6]= Delta_t * (1-delta);
	a_dt[7]= delta * Delta_t;

}


void CFEM_ElasticitySolver::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	unsigned short iVar, jVar;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	unsigned short iNode, jNode;
	double checkJacobian;

	ofstream myfile;
	myfile.open ("jacobianPreClamped.txt");

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

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		/*--- Get node index ---*/

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (nDim == 2) {
			Solution[0] = 0.0;  Solution[1] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;
		}
		else {
			Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
			Residual[0] = 0.0;  Residual[1] = 0.0;  Residual[2] = 0.0;
		}

		node[iPoint]->SetSolution(Solution);

		if (dynamic){
			node[iPoint]->SetSolution_Vel(Solution);
			node[iPoint]->SetSolution_Accel(Solution);
		}

		LinSysRes.SetBlock(iPoint, Residual);

		/*--- STRONG ENFORCEMENT OF THE DISPLACEMENT BOUNDARY CONDITION ---*/

		/*--- Delete the columns for a particular node ---*/

		for (iVar = 0; iVar < nPoint; iVar++){
			if (iVar==iPoint) {
				Jacobian.SetBlock(iVar,iPoint,mId_Aux);
			}
			else {
				Jacobian.SetBlock(iVar,iPoint,mZeros_Aux);
			}
		}

		/*--- Delete the rows for a particular node ---*/
		for (jVar = 0; jVar < nPoint; jVar++){
			if (iPoint!=jVar) {
				Jacobian.SetBlock(jVar,iPoint,mZeros_Aux);
			}
		}

		/*--- If the problem is dynamic ---*/
		/*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/

		if(dynamic){

			node[iPoint]->SetSolution_time_n(Solution);
			node[iPoint]->SetSolution_Vel_time_n(Solution);
			node[iPoint]->SetSolution_Accel_time_n(Solution);

		}

	}

	myfile.open ("jacobianClamped.txt");

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

void CFEM_ElasticitySolver::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                        unsigned short val_marker) {

	unsigned long iPoint, iVertex;
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

		/*--- Get node index ---*/

		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

		if (nDim == 2) {
			Solution[0] = 0.0;  Solution[1] = 0.0;
		}
		else {
			Solution[0] = 0.0;  Solution[1] = 0.0;  Solution[2] = 0.0;
		}

		node[iPoint]->SetSolution(Solution);

		if (dynamic){
			node[iPoint]->SetSolution_Vel(Solution);
			node[iPoint]->SetSolution_Accel(Solution);
		}

	}

}

void CFEM_ElasticitySolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
		unsigned short iMesh) {

    unsigned short iVar;
	unsigned long iPoint, total_index;

	bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);		// Geometrically linear problems
	bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems



	/*--- Update solution and (if dynamic) advance velocity and acceleration vectors in time ---*/

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

		for (iVar = 0; iVar < nVar; iVar++) {

			/*--- Displacements component of the solution ---*/

			/*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/

			if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

			if (nonlinear)	node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

		}

	}

	/*--- MPI solution ---*/

	Set_MPI_Solution(geometry, config);

	/*---  Compute the residual Ax-f ---*/

	Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

	  /*--- Set maximum residual to zero ---*/

		for (iVar = 0; iVar < nVar; iVar++) {
			SetRes_RMS(iVar, 0.0);
			SetRes_Max(iVar, 0.0, 0);
		}

	  /*--- Compute the residual ---*/

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			for (iVar = 0; iVar < nVar; iVar++) {
				total_index = iPoint*nVar+iVar;
				AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
				AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
			}
		}

	  /*--- Compute the root mean square residual ---*/

	  SetResidual_RMS(geometry, config);

}

void CFEM_ElasticitySolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) {

	double a[3], b[3], AC[3], BD[3];
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
	double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	double Length_Elem = 0.0, Area_Elem = 0.0, Normal_Elem[3] = {0.0, 0.0, 0.0};
	unsigned short iDim;

	double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
	double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
	double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));

	double TotalLoad;

  bool Gradual_Load = config->GetGradual_Load();
	double CurrentTime=config->GetCurrent_DynTime();
	double ModAmpl, NonModAmpl;

  bool Ramp_Load = config->GetRamp_Load();
	double Ramp_Time = config->GetRamp_Time();

	if (Ramp_Load){
		ModAmpl=LoadDirVal*LoadDirMult*CurrentTime/Ramp_Time;
		NonModAmpl=LoadDirVal*LoadDirMult;
		TotalLoad=min(ModAmpl,NonModAmpl);
	}
	else if (Gradual_Load){
		ModAmpl=2*((1/(1+exp(-1*CurrentTime)))-0.5);
		TotalLoad=ModAmpl*LoadDirVal*LoadDirMult;
	}
	else{
		TotalLoad=LoadDirVal*LoadDirMult;
	}

	/*--- Compute the norm of the vector that was passed in the config file ---*/
	double Norm;
	if (nDim==2) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]);
	if (nDim==3) Norm=sqrt(Load_Dir_Local[0]*Load_Dir_Local[0]+Load_Dir_Local[1]*Load_Dir_Local[1]+Load_Dir_Local[2]*Load_Dir_Local[2]);

	for (iElem = 0; iElem < geometry->GetnElem_Bound(val_marker); iElem++) {

		Point_0 = geometry->bound[val_marker][iElem]->GetNode(0);     Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->bound[val_marker][iElem]->GetNode(1);     Coord_1 = geometry->node[Point_1]->GetCoord();
		if (nDim == 3) {

			Point_2 = geometry->bound[val_marker][iElem]->GetNode(2);	Coord_2 = geometry->node[Point_2]->GetCoord();
		    if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){
		    	Point_3 = geometry->bound[val_marker][iElem]->GetNode(3);	Coord_3 = geometry->node[Point_3]->GetCoord();
		    }

		}

		/*--- Compute area (3D), and length of the surfaces (2D) ---*/

		if (nDim == 2) {

			for (iDim = 0; iDim < nDim; iDim++) a[iDim] = Coord_0[iDim]-Coord_1[iDim];

			Length_Elem = sqrt(a[0]*a[0]+a[1]*a[1]);
			Normal_Elem[0] =   a[1];
			Normal_Elem[1] = -(a[0]);

		}

		if (nDim == 3) {

			if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){

				for (iDim = 0; iDim < nDim; iDim++) {
					a[iDim] = Coord_1[iDim]-Coord_0[iDim];
					b[iDim] = Coord_2[iDim]-Coord_0[iDim];
				}

				double Ni=0 , Nj=0, Nk=0;

				Ni=a[1]*b[2]-a[2]*b[1];
				Nj=-a[0]*b[2]+a[2]*b[0];
				Nk=a[0]*b[1]-a[1]*b[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);


				//Area_Elem = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);

			}

			else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

				for (iDim = 0; iDim < nDim; iDim++) {
					AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
					BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
				}

				double Ni=0 , Nj=0, Nk=0;

				Ni=AC[1]*BD[2]-AC[2]*BD[1];
				Nj=-AC[0]*BD[2]+AC[2]*BD[0];
				Nk=AC[0]*BD[1]-AC[1]*BD[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);

			}
		}

      if (nDim == 2) {

        Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;

//        LinSysRes.AddBlock(Point_0, Residual);
//        LinSysRes.AddBlock(Point_1, Residual);

        /*--- This will allow us to compute the boundary conditions only once ---*/
        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);

      }

      else {
    	  if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){

    		  Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

//    		  LinSysRes.AddBlock(Point_0, Residual);
//    		  LinSysRes.AddBlock(Point_1, Residual);
//    		  LinSysRes.AddBlock(Point_2, Residual);

    	      /*--- This will allow us to compute the boundary conditions only once ---*/
    	      node[Point_0]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_1]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_2]->Add_SurfaceLoad_Res(Residual);

    	  }
    	  else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

    		  Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

//    		  LinSysRes.AddBlock(Point_0, Residual);
//    		  LinSysRes.AddBlock(Point_1, Residual);
//    		  LinSysRes.AddBlock(Point_2, Residual);
//    		  LinSysRes.AddBlock(Point_3, Residual);

    	      /*--- This will allow us to compute the boundary conditions only once ---*/
    	      node[Point_0]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_1]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_2]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_3]->Add_SurfaceLoad_Res(Residual);


    	  }

      }

	}

}

void CFEM_ElasticitySolver::BC_Sine_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Pressure(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) { }

void CFEM_ElasticitySolver::ImplicitNewmark_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

	unsigned long iPoint, jPoint;
	unsigned short iVar, jVar;

	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool first_iter = (config->GetIntIter() == 0);
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method


	if (!dynamic){

		for (iPoint = 0; iPoint < nPoint; iPoint++){
			Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
		}

	}

	if (dynamic) {

		/*--- Add the mass matrix contribution to the Jacobian ---*/

		/*
		 * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
		 * of each time step. If the solution method is Newton Rapshon, we initialize it also at the beginning of each
		 * iteration.
		 */

		if ((nonlinear_analysis) && ((newton_raphson) || (first_iter))){
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
				for (jPoint = 0; jPoint < geometry->GetnPoint(); jPoint++){
					for(iVar = 0; iVar < nVar; iVar++){
						for (jVar = 0; jVar < nVar; jVar++){
							Jacobian_ij[iVar][jVar] = MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
						}
					}
				}

				Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
			}
		}

		/*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
		if (linear_analysis){
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				for (iVar = 0; iVar < nVar; iVar++){
					Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+		//a0*U(t)
								 	 a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+	//a2*U'(t)
								 	 a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
				}
				TimeRes_Aux.SetBlock(iPoint, Residual);
			}
		}
		else if (nonlinear_analysis){
			for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
				for (iVar = 0; iVar < nVar; iVar++){
					Residual[iVar] =   a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)  			//a0*U(t)
									 - a_dt[0]*node[iPoint]->GetSolution(iVar) 					//a0*U(t+dt)(k-1)
								 	 + a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)		//a2*U'(t)
								 	 + a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
				}
				TimeRes_Aux.SetBlock(iPoint, Residual);
			}
		}
		/*--- Once computed, compute M*TimeRes_Aux ---*/
		MassMatrix.MatrixVectorProduct(TimeRes_Aux,TimeRes,geometry,config);
		/*--- Add the components of M*TimeRes_Aux to the residual R(t+dt) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			Res_Time_Cont = TimeRes.GetBlock(iPoint);
			Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			LinSysRes.AddBlock(iPoint, Res_Time_Cont);
			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
		}
	}



}

void CFEM_ElasticitySolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config){

    unsigned short iVar;
	unsigned long iPoint, total_index, IterLinSol;

	unsigned short iNode, jNode, jVar;

	double checkJacobian;

	ofstream myfile;
	myfile.open ("Jacobian.txt");

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

//	myfile.open ("StiffMatrix.txt");
//
//	for (iNode = 0; iNode < nPoint; iNode++){
//		for (jNode = 0; jNode < nPoint; jNode++){
//			myfile << "Node " << iNode << " " << jNode << endl;
//			for (iVar = 0; iVar < nVar; iVar++){
//				for (jVar = 0; jVar < nVar; jVar++){
//					checkJacobian = StiffMatrix.GetBlock(iNode, jNode, iVar, jVar);
//					myfile << checkJacobian << " " ;
//				}
//				myfile << endl;
//			}
//		}
//	}
//	myfile.close();

	CSysSolve femSystem;
	IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

}
