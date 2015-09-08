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

	Total_CFEA = 0.0;
	WAitken_Dyn = 0.0;
	WAitken_Dyn_tn1 = 0.0;
	loadIncrement = 0.0;

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
	Res_FSI_Cont = NULL;

	nodeReactions = NULL;

	solutionPredictor = NULL;

	SolRest = NULL;

	normalVertex = NULL;
	stressTensor = NULL;

}

CFEM_ElasticitySolver::CFEM_ElasticitySolver(CGeometry *geometry, CConfig *config) : CSolver() {

	unsigned long iPoint, iElem = 0;
	unsigned short iVar, jVar, iDim, jDim, NodesElement = 0, nKindElements;

	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method
	bool fsi = config->GetFSI_Simulation();												// FSI simulation

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	su2double E = config->GetElasticyMod();

	nElement      = geometry->GetnElem();
	nDim          = geometry->GetnDim();
	nMarker       = geometry->GetnMarker();

	nPoint        = geometry->GetnPoint();
	nPointDomain  = geometry->GetnPointDomain();

	nKindElements = 2;

	element_container = new CElement*[nKindElements];
	node          	  = new CVariable*[nPoint];

	GradN_X = new su2double [nDim];
	GradN_x = new su2double [nDim];

	Total_CFEA			= 0.0;
	WAitken_Dyn      	= 0.0;
	WAitken_Dyn_tn1  	= 0.0;
	loadIncrement 		= 0.0;

  	SetFSI_ConvValue(0,0.0);
  	SetFSI_ConvValue(1,0.0);

	nVar = nDim;

	/*--- Define some auxiliary vectors related to the residual ---*/

	Residual = new su2double[nVar];          for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
	Residual_RMS = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
	Residual_Max = new su2double[nVar];      for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
	Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
	Point_Max_Coord = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Point_Max_Coord[iVar] = new su2double[nDim];
		for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
	}

	/*--- Define some auxiliary vectors related to the solution ---*/

	Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

	nodeReactions = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) nodeReactions[iVar]   = 0.0;

	/*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/


	unsigned short nSolVar;
	unsigned long index;
	string text_line, filename;
	ifstream restart_file;
	su2double dull_val;
	long Dyn_RestartIter;

	if (dynamic) nSolVar = 3 * nVar;
	else nSolVar = nVar;

	SolRest = new su2double[nSolVar];

	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	/*--- Check for a restart, initialize from zero otherwise ---*/

	if (!restart) {
		for (iVar = 0; iVar < nSolVar; iVar++) SolRest[iVar] = 0.0;
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
		}
	}
	else {

		/*--- Restart the solution from file information ---*/

		filename = config->GetSolution_FEMFileName();

	    if (dynamic) {

	      Dyn_RestartIter = SU2_TYPE::Int(config->GetDyn_RestartIter())-1;

	      filename = config->GetUnsteady_FileName(filename, Dyn_RestartIter);
	    }

		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no file ---*/

		if (restart_file.fail()) {
		    if (rank == MASTER_NODE)
			cout << "There is no FEM restart file!!" << endl;
			exit(EXIT_FAILURE);
		}

		/*--- In case this is a parallel simulation, we need to perform the
     	 Global2Local index transformation first. ---*/

		long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];

		/*--- First, set all indices to a negative value by default ---*/

	    for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
	      Global2Local[iPoint] = -1;

		/*--- Now fill array with the transform values only for local points ---*/

	    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
	      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;

	    /*--- Read all lines in the restart file ---*/

	    long iPoint_Local;
	    unsigned long iPoint_Global_Local = 0, iPoint_Global = 0; string text_line;
	    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

		/*--- The first line is the header ---*/

		getline (restart_file, text_line);

		while (getline (restart_file, text_line)) {
		   istringstream point_line(text_line);

		/*--- Retrieve local index. If this node from the restart file lives
       	   on a different processor, the value of iPoint_Local will be -1.
       	   Otherwise, the local index for this node on the current processor
       	   will be returned and used to instantiate the vars. ---*/

		   iPoint_Local = Global2Local[iPoint_Global];

			if (iPoint_Local >= 0) {
				if (dynamic){
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2] >> SolRest[3] >> SolRest[4] >> SolRest[5];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2] >> SolRest[3] >> SolRest[4] >> SolRest[5] >> SolRest[6] >> SolRest[7] >> SolRest[8];
				}
				else {
					if (nDim == 2) point_line >> index >> dull_val >> dull_val >> SolRest[0] >> SolRest[1];
					if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> SolRest[0] >> SolRest[1] >> SolRest[2];
				}

				node[iPoint_Local] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
				iPoint_Global_Local++;
			}
			iPoint_Global++;
		}

	    /*--- Detect a wrong solution file ---*/

	    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }

		#ifndef HAVE_MPI
	    	rbuf_NotMatching = sbuf_NotMatching;
		#else
	    	SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
		#endif

	    if (rbuf_NotMatching != 0) {
	      if (rank == MASTER_NODE) {
	        cout << endl << "The solution file " << filename.data() << " doesn't match with the mesh file!" << endl;
	        cout << "It could be empty lines at the end of the file." << endl << endl;
	      }
		#ifndef HAVE_MPI
	      	  exit(EXIT_FAILURE);
		#else
	      	  MPI_Barrier(MPI_COMM_WORLD);
	      	  MPI_Abort(MPI_COMM_WORLD,1);
	      	  MPI_Finalize();
		#endif
	    }

		/*--- Instantiate the variable class with an arbitrary solution
     	 at any halo/periodic nodes. The initial solution can be arbitrary,
     	 because a send/recv is performed immediately in the solver (Set_MPI_Solution()). ---*/

		for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CFEM_ElasVariable(SolRest, nDim, nVar, config);
		}


		/*--- Close the restart file ---*/

		restart_file.close();

		/*--- Free memory needed for the transformation ---*/

		delete [] Global2Local;

	}

	bool incompressible = (config->GetMaterialCompressibility() == INCOMPRESSIBLE_MAT);

	/*--- Term ij of the Mass Matrix ---*/

	MassMatrix_ij = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		MassMatrix_ij[iVar] = new su2double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				MassMatrix_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_ij[iVar] = new su2double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (constitutive contribution) ---*/

	Jacobian_c_ij = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_c_ij[iVar] = new su2double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_c_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (stress contribution) ---*/

	Jacobian_s_ij = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_s_ij[iVar] = new su2double [nVar];
			for (jVar = 0; jVar < nVar; jVar++) {
				Jacobian_s_ij[iVar][jVar] = 0.0;
			}
	}

	/*--- Term ij of the Jacobian (incompressibility term) ---*/

	if (incompressible){
		Jacobian_k_ij = new su2double*[nVar];
		for (iVar = 0; iVar < nVar; iVar++) {
			Jacobian_k_ij[iVar] = new su2double [nVar];
				for (jVar = 0; jVar < nVar; jVar++) {
					Jacobian_k_ij[iVar][jVar] = 0.0;
				}
		}
	}
	else {
		Jacobian_k_ij = NULL;
	}

	/*--- Stress contribution to the node i ---*/
	Res_Stress_i = new su2double[nVar];

	/*--- Contribution of the external surface forces to the residual (auxiliary vector) ---*/
	Res_Ext_Surf = new su2double[nVar];

	/*--- Contribution of the fluid tractions to the residual (auxiliary vector) ---*/
	if (fsi){
		Res_FSI_Cont = new su2double[nVar];
	}
	else {
		Res_FSI_Cont = NULL;
	}


	/*--- Time integration contribution to the residual ---*/
	if (dynamic) {
		Res_Time_Cont = new su2double [nVar];
	}
	else {
		Res_Time_Cont = NULL;
	}

	/*--- Matrices to impose clamped boundary conditions (TODO: Initialize them conditionally). ---*/

	mZeros_Aux = new su2double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mZeros_Aux[iDim] = new su2double[nDim];

	mId_Aux = new su2double *[nDim];
	for(iDim = 0; iDim < nDim; iDim++)
		mId_Aux[iDim] = new su2double[nDim];

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

	if (dynamic) {
		MassMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
		TimeRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);
		TimeRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
	}


	/*--- Initialization of linear solver structures ---*/
	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysReact.Initialize(nPoint, nPointDomain, nVar, 0.0);

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

	/*--- Initialize the auxiliary vector and matrix for the computation of the nodal Reactions ---*/

	normalVertex = new su2double [nDim];

	stressTensor = new su2double* [nDim];
	for (iVar = 0; iVar < nVar; iVar++){
		stressTensor[iVar] = new su2double [nDim];
	}

	/*---- Initialize the auxiliary vector for the solution predictor ---*/

	solutionPredictor = new su2double [nVar];

	/*--- Perform the MPI communication of the solution ---*/

	Set_MPI_Solution(geometry, config);

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
		delete [] stressTensor[iVar];
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
	delete [] SolRest;
	delete [] GradN_X;
	delete [] GradN_x;

	delete [] Residual;
	delete [] Residual_RMS;
	delete [] Residual_Max;
	delete [] Point_Max;
	delete [] Point_Max_Coord;

	delete [] mZeros_Aux;
	delete [] mId_Aux;

	delete [] nodeReactions;

	delete [] normalVertex;
	delete [] stressTensor;

}

void CFEM_ElasticitySolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {


	  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

	  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.

	  unsigned short nSolVar;

	  if (dynamic) nSolVar = 3 * nVar;
	  else nSolVar = nVar;

	#ifdef HAVE_MPI
	  int send_to, receive_from;
	  MPI_Status status;
	#endif

	  for (iMarker = 0; iMarker < nMarker; iMarker++) {

	    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
	        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

	      MarkerS = iMarker;  MarkerR = iMarker+1;

	#ifdef HAVE_MPI
	      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
	      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
	#endif

	      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
	      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;

	      /*--- Allocate Receive and send buffers  ---*/
	      Buffer_Receive_U = new su2double [nBufferR_Vector];
	      Buffer_Send_U = new su2double[nBufferS_Vector];

	      /*--- Copy the solution that should be sent ---*/
	      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
	        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
	        if (dynamic){
	            for (iVar = 0; iVar < nVar; iVar++){
	            	Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
	            	Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
	            }
	        }
	      }

	#ifdef HAVE_MPI

	      /*--- Send/Receive information using Sendrecv ---*/
	      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
	                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

	#else

	      /*--- Receive information without MPI ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
	        if (dynamic){
	            for (iVar = nVar; iVar < 3*nVar; iVar++)
	            	Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
	        }
	      }

	#endif

	      /*--- Deallocate send buffer ---*/
	      delete [] Buffer_Send_U;

	      /*--- Do the coordinate transformation ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

	        /*--- Find point and its type of transformation ---*/
	        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

	        /*--- Copy solution variables. ---*/
	        for (iVar = 0; iVar < nSolVar; iVar++)
	        	SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

	        /*--- Store received values back into the variable. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          node[iPoint]->SetSolution(iVar, SolRest[iVar]);

	        if (dynamic){

	            for (iVar = 0; iVar < nVar; iVar++){
	            	node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
	            	node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
	            }

	        }

	      }

	      /*--- Deallocate receive buffer ---*/
	      delete [] Buffer_Receive_U;

	    }

	  }

}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred(CGeometry *geometry, CConfig *config) {


	  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

	#ifdef HAVE_MPI
	  int send_to, receive_from;
	  MPI_Status status;
	#endif

	  for (iMarker = 0; iMarker < nMarker; iMarker++) {

	    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
	        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

	      MarkerS = iMarker;  MarkerR = iMarker+1;

	#ifdef HAVE_MPI
	      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
	      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
	#endif

	      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
	      nBufferS_Vector = nVertexS*nVar;     nBufferR_Vector = nVertexR*nVar;

	      /*--- Allocate Receive and send buffers  ---*/
	      Buffer_Receive_U = new su2double [nBufferR_Vector];
	      Buffer_Send_U = new su2double[nBufferS_Vector];

	      /*--- Copy the solution that should be sent ---*/
	      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
	        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
	      }

	#ifdef HAVE_MPI

	      /*--- Send/Receive information using Sendrecv ---*/
	      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
	                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

	#else

	      /*--- Receive information without MPI ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
	      }

	#endif

	      /*--- Deallocate send buffer ---*/
	      delete [] Buffer_Send_U;

	      /*--- Do the coordinate transformation ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

	        /*--- Find point and its type of transformation ---*/
	        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

	        /*--- Copy predicted solution variables back into the variables. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

	      }

	      /*--- Deallocate receive buffer ---*/
	      delete [] Buffer_Receive_U;

	    }

	  }

}

void CFEM_ElasticitySolver::Set_MPI_Solution_Pred_Old(CGeometry *geometry, CConfig *config) {

	  /*--- We are communicating the solution predicted, current and old, and the old solution ---*/
	  /*--- necessary for the Aitken relaxation ---*/

	  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
	  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
	  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

	  /*--- Analogous to the dynamic solution, in this case we need 3 * nVar variables per node ---*/
	  unsigned short nSolVar;
	  nSolVar = 3 * nVar;

	#ifdef HAVE_MPI
	  int send_to, receive_from;
	  MPI_Status status;
	#endif

	  for (iMarker = 0; iMarker < nMarker; iMarker++) {

	    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
	        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

	      MarkerS = iMarker;  MarkerR = iMarker+1;

	#ifdef HAVE_MPI
	      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
	      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
	#endif

	      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
	      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;

	      /*--- Allocate Receive and send buffers  ---*/
	      Buffer_Receive_U = new su2double [nBufferR_Vector];
	      Buffer_Send_U = new su2double[nBufferS_Vector];

	      /*--- Copy the solution that should be sent ---*/
	      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
	        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
	        for (iVar = 0; iVar < nVar; iVar++){
	        	Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
	        	Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred(iVar);
	        	Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Pred_Old(iVar);
	        }
	      }

	#ifdef HAVE_MPI

	      /*--- Send/Receive information using Sendrecv ---*/
	      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
	                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

	#else

	      /*--- Receive information without MPI ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
	        for (iVar = 0; iVar < nSolVar; iVar++)
	          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
	      }

	#endif

	      /*--- Deallocate send buffer ---*/
	      delete [] Buffer_Send_U;

	      /*--- Do the coordinate transformation ---*/
	      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

	        /*--- Find point and its type of transformation ---*/
	        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

	        /*--- Store received values back into the variable. ---*/
	        for (iVar = 0; iVar < nVar; iVar++){
	        	node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);
	        	node[iPoint]->SetSolution_Pred(iVar, Buffer_Receive_U[(iVar+nVar)*nVertexR+iVertex]);
	        	node[iPoint]->SetSolution_Pred_Old(iVar, Buffer_Receive_U[(iVar+2*nVar)*nVertexR+iVertex]);
	        }

	      }

	      /*--- Deallocate receive buffer ---*/
	      delete [] Buffer_Receive_U;

	    }

	  }

}


void CFEM_ElasticitySolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {


	unsigned long iPoint;
	bool initial_calc = (config->GetExtIter() == 0);									// Checks if it is the first calculation.
	bool first_iter = (config->GetIntIter() == 0);													// Checks if it is the first iteration
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.
	bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.
	bool newton_raphson = (config->GetKind_SpaceIteScheme_FEA() == NEWTON_RAPHSON);		// Newton-Raphson method
	bool restart = config->GetRestart();												// Restart analysis
	bool initial_calc_restart = (config->GetExtIter() == config->GetDyn_RestartIter()); // Initial calculation for restart

	/*--- Set vector entries to zero ---*/
    //TODO: nPoint or nPointDomain
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
	if ((initial_calc && linear_analysis)||
		(restart && initial_calc_restart && linear_analysis)){
		Jacobian.SetValZero();
	}

	/*
	 * If the problem is dynamic, we need a mass matrix, which will be constant along the calculation
	 * both for linear and nonlinear analysis. Only initialized once, at the first time step.
	 *
	 * The same with the integration constants, as for now we consider the time step to be constant.
	 *
	 * We need first_iter, because in nonlinear problems there are more than one subiterations in the first time step.
	 */
	if ((dynamic && initial_calc && first_iter) ||
		(dynamic && restart && initial_calc_restart && first_iter)) {
		MassMatrix.SetValZero();
		Compute_IntegrationConstants(config);
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

void CFEM_ElasticitySolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

	unsigned long iPoint, nPoint;
	bool incremental_load = config->GetIncrementalLoad();							// If an incremental load is applied

	nPoint = geometry[MESH_0]->GetnPoint();

	/*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

	if (incremental_load){
		for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
	}


}

void CFEM_ElasticitySolver::ResetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

	unsigned long iPoint, nPoint;
	bool incremental_load = config->GetIncrementalLoad();							// If an incremental load is applied

	nPoint = geometry[MESH_0]->GetnPoint();

	/*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

	if (incremental_load){
		for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_Solution();
	}

}

void CFEM_ElasticitySolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord, val_Sol;
	int EL_KIND;

	su2double Ks_ab;
	su2double *Kab = NULL;
	su2double *Kk_ab = NULL;
	su2double *Ta = NULL;
	unsigned short NelNodes, jNode;

	su2double checkJacobian, *checkCoord;

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
	su2double val_Coord, val_Sol;
	int EL_KIND;

	su2double Ks_ab;
	su2double *Kab = NULL;
	su2double *Kk_ab = NULL;
	su2double *Ta = NULL;
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

		/*--- If incompressible, we compute the Mean Dilatation term first so the volume is already computed ---*/

		if (incompressible) numerics->Compute_MeanDilatation_Term(element_container[EL_KIND]);

		numerics->Compute_Tangent_Matrix(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			/*--- Check if this is my node or not ---*/
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

				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_c_ij);
				Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_s_ij);
				if (incompressible) Jacobian.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_k_ij);

			}

		}

	}

}

void CFEM_ElasticitySolver::Compute_MassMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord, val_Sol;
	int EL_KIND;

	su2double Mab;
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

}

void CFEM_ElasticitySolver::Compute_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {


	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iGauss, iDim;
	unsigned short nNodes, nGauss;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord, val_Sol;
	int EL_KIND;

	su2double Ks_ab;
	su2double *Kab = NULL;
	su2double *Kk_ab = NULL;
	su2double *Ta = NULL;
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
	}

}

void CFEM_ElasticitySolver::Compute_NodalStress(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {

	unsigned long iPoint, iElem, iVar, jVar;
	unsigned short iNode, iDim, iStress;
	unsigned short nNodes, nStress;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord, val_Sol;
	int EL_KIND;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	if (nDim == 2) nStress = 3;
	else if (nDim == 3) nStress = 6;

	su2double *Ta = NULL;

	unsigned short NelNodes;

	/*--- Restart stress to avoid adding results from previous time steps ---*/

	 for (iPoint = 0; iPoint < nPointDomain; iPoint++){
		for (iStress = 0; iStress < nStress; iStress++){
				node[iPoint]->SetStress_FEM(iStress, 0.0);
		}
	}

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

		numerics->Compute_Averaged_NodalStress(element_container[EL_KIND]);

		NelNodes = element_container[EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			/*--- This only works if the problem is nonlinear ---*/
			Ta = element_container[EL_KIND]->Get_Kt_a(iNode);
			for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

			LinSysReact.AddBlock(indexNode[iNode], Res_Stress_i);

			for (iStress = 0; iStress < nStress; iStress++){
				node[indexNode[iNode]]->AddStress_FEM(iStress,
						(element_container[EL_KIND]->Get_NodalStress(iNode, iStress) /
								geometry->node[indexNode[iNode]]->GetnElem()) );
			}

		}

	}

	  su2double *Stress;
	  su2double VonMises_Stress, MaxVonMises_Stress = 0.0;
	  su2double Sxx,Syy,Szz,Sxy,Sxz,Syz,S1,S2;

	 /* --- For the number of nodes in the mesh ---*/
	  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		  /* --- Get the stresses, added up from all the elements that connect to the node ---*/

		  Stress  = node[iPoint]->GetStress_FEM();

		  /* --- Compute the stress averaged from all the elements connecting to the node and the Von Mises stress ---*/

		  if (geometry->GetnDim() == 2) {

			  Sxx=Stress[0];
			  Syy=Stress[1];
			  Sxy=Stress[2];

			  S1=(Sxx+Syy)/2+sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);
			  S2=(Sxx+Syy)/2-sqrt(((Sxx-Syy)/2)*((Sxx-Syy)/2)+Sxy*Sxy);

			  VonMises_Stress = sqrt(S1*S1+S2*S2-2*S1*S2);

		  }
		  else if (geometry->GetnDim() == 3) {

			  Sxx = Stress[0];
			  Syy = Stress[1];
			  Szz = Stress[3];

			  Sxy = Stress[2];
			  Sxz = Stress[4];
			  Syz = Stress[5];

			  VonMises_Stress = sqrt(0.5*(   pow(Sxx - Syy, 2.0)
											+ pow(Syy - Szz, 2.0)
											+ pow(Szz - Sxx, 2.0)
											+ 6.0*(Sxy*Sxy+Sxz*Sxz+Syz*Syz)
											));

		  }

		  node[iPoint]->SetVonMises_Stress(VonMises_Stress);

		  /*--- Compute the maximum value of the Von Mises Stress ---*/

		  MaxVonMises_Stress = max(MaxVonMises_Stress, VonMises_Stress);

	  }

  	su2double checkJacobian;
  	unsigned long jNode;

  	bool outputReactions = false;

	if (outputReactions) {

	  	ofstream myfile;
	  	myfile.open ("Reactions.txt");

		unsigned short iMarker;
		unsigned long iVertex;
		su2double val_Reaction;

		bool linear_analysis = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Linear analysis.
		bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.

		if (!dynamic){
			/*--- Loop over all the markers  ---*/
			for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
				switch (config->GetMarker_All_KindBC(iMarker)) {

					/*--- If it corresponds to a clamped boundary  ---*/

					case CLAMPED_BOUNDARY:

						myfile << "MARKER " << iMarker << ":" << endl;

						/*--- Loop over all the vertices  ---*/
						for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

						/*--- Get node index ---*/
						iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

						myfile << "Node " << iPoint << "." << " \t ";

						for (iDim = 0; iDim < nDim; iDim++){
							/*--- Retrieve coordinate ---*/
							val_Coord = geometry->node[iPoint]->GetCoord(iDim);
							myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
						}

						for (iVar = 0; iVar < nVar; iVar++){
							/*--- Retrieve reaction ---*/
							val_Reaction = LinSysReact.GetBlock(iPoint, iVar);
							myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
						}

						myfile << endl;
					}
				  myfile << endl;
				  break;
			}
		}
		else if (dynamic){

			switch (config->GetKind_TimeIntScheme_FEA()) {
				case (CD_EXPLICIT):
						  cout << "NOT IMPLEMENTED YET" << endl;
				  break;
				case (NEWMARK_IMPLICIT):

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

					/*--- Loop over all the markers  ---*/
					for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
						switch (config->GetMarker_All_KindBC(iMarker)) {

							/*--- If it corresponds to a clamped boundary  ---*/

							case CLAMPED_BOUNDARY:

							myfile << "MARKER " << iMarker << ":" << endl;

								/*--- Loop over all the vertices  ---*/
								for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

								/*--- Get node index ---*/
								iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

								myfile << "Node " << iPoint << "." << " \t ";

								for (iDim = 0; iDim < nDim; iDim++){
									/*--- Retrieve coordinate ---*/
									val_Coord = geometry->node[iPoint]->GetCoord(iDim);
									myfile << "X" << iDim + 1 << ": " << val_Coord << " \t " ;
								}

								/*--- Retrieve the time contribution ---*/
								Res_Time_Cont = TimeRes.GetBlock(iPoint);

								for (iVar = 0; iVar < nVar; iVar++){
									/*--- Retrieve reaction ---*/
									val_Reaction = LinSysReact.GetBlock(iPoint, iVar) + Res_Time_Cont[iVar];
									myfile << "F" << iVar + 1 << ": " << val_Reaction << " \t " ;
								}

								myfile << endl;
							}
						  myfile << endl;
						  break;
					}


				  break;
				case (GA_IMPLICIT):
						  cout << "NOT IMPLEMENTED YET" << endl;
				  break;
			  }

		}



		myfile.close();

	}

	#ifdef HAVE_MPI

	  /*--- Compute MaxVonMises_Stress using all the nodes ---*/

	  su2double MyMaxVonMises_Stress = MaxVonMises_Stress; MaxVonMises_Stress = 0.0;
	  SU2_MPI::Allreduce(&MyMaxVonMises_Stress, &MaxVonMises_Stress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	#endif

		  /*--- Set the value of the MaxVonMises_Stress as the CFEA coeffient ---*/

	  Total_CFEA = MaxVonMises_Stress;

}

void CFEM_ElasticitySolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CFEM_ElasticitySolver::Compute_IntegrationConstants(CConfig *config) {

	su2double Delta_t= config->GetDelta_DynTime();
	su2double delta = config->GetNewmark_delta(), alpha = config->GetNewmark_alpha();

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

	su2double tempCoord;

	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);

	unsigned short iNode, jNode;


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

//		for (iVar = 0; iVar < nVar; iVar++){
//			nodeReactions[iVar] = - 1.0 * LinSysRes.GetBlock(iPoint, iVar);
//		}
//
//		LinSysReact.SetBlock(iPoint,nodeReactions);

		/*--- Initialize the reaction vector ---*/
		LinSysReact.SetBlock(iPoint, Residual);


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
				Jacobian.SetBlock(iPoint,jVar,mZeros_Aux);
			}
		}

		/*--- If the problem is dynamic ---*/
		/*--- Enforce that in the previous time step all nodes had 0 U, U', U'' ---*/
		/*--- TODO: Do I really need to do this? ---*/

		if(dynamic){

			node[iPoint]->SetSolution_time_n(Solution);
			node[iPoint]->SetSolution_Vel_time_n(Solution);
			node[iPoint]->SetSolution_Accel_time_n(Solution);

		}

	}

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

	bool first_iter = (config->GetIntIter() == 0);
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);		// Nonlinear analysis.

	su2double solNorm = 0.0, tempCheck[3];

	if (nonlinear_analysis){

		/*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

		/*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
		/*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
		/*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

		if (first_iter){
			Conv_Ref[0] = 1.0;											// Position for the norm of the solution
			Conv_Ref[1] = max(LinSysRes.norm(), EPS);					// Position for the norm of the residual
			Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);		// Position for the energy tolerance

			/*--- Make sure the computation runs at least 2 iterations ---*/
			Conv_Check[0] = 1.0;
			Conv_Check[1] = 1.0;
			Conv_Check[2] = 1.0;
		}
		else {
			/*--- Compute the norm of the solution vector Uk ---*/
			for (iPoint = 0; iPoint < nPointDomain; iPoint++){
				for (iVar = 0; iVar < nVar; iVar++){
					solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
				}
			}
			Conv_Ref[0] = max(sqrt(solNorm), EPS);							// Norm of the solution vector

			Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];					// Norm of the delta-solution vector
			Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];					// Norm of the residual
			Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];	// Position for the energy tolerance
		}

	}
	else{

			/*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

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


			  /*--- MPI solution ---*/

			  Set_MPI_Solution(geometry, config);

			  /*--- Compute the root mean square residual ---*/

			  SetResidual_RMS(geometry, config);
	}

}

void CFEM_ElasticitySolver::BC_Normal_Displacement(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                         unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Normal_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) { }

void CFEM_ElasticitySolver::BC_Dir_Load(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
        unsigned short val_marker) {

	su2double a[3], b[3], AC[3], BD[3];
	unsigned long iElem, Point_0 = 0, Point_1 = 0, Point_2 = 0, Point_3=0;
	su2double *Coord_0 = NULL, *Coord_1= NULL, *Coord_2= NULL, *Coord_3= NULL;
	su2double Length_Elem = 0.0, Area_Elem = 0.0, Normal_Elem[3] = {0.0, 0.0, 0.0};
	unsigned short iDim;

	su2double LoadDirVal = config->GetLoad_Dir_Value(config->GetMarker_All_TagBound(val_marker));
	su2double LoadDirMult = config->GetLoad_Dir_Multiplier(config->GetMarker_All_TagBound(val_marker));
	su2double *Load_Dir_Local= config->GetLoad_Dir(config->GetMarker_All_TagBound(val_marker));

	su2double TotalLoad;

  bool Gradual_Load = config->GetGradual_Load();
	su2double CurrentTime=config->GetCurrent_DynTime();
	su2double ModAmpl, NonModAmpl;

  bool Ramp_Load = config->GetRamp_Load();
	su2double Ramp_Time = config->GetRamp_Time();

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
	su2double Norm;
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

				su2double Ni=0 , Nj=0, Nk=0;

				Ni=a[1]*b[2]-a[2]*b[1];
				Nj=-a[0]*b[2]+a[2]*b[0];
				Nk=a[0]*b[1]-a[1]*b[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);

			}

			else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

				for (iDim = 0; iDim < nDim; iDim++) {
					AC[iDim] = Coord_2[iDim]-Coord_0[iDim];
					BD[iDim] = Coord_3[iDim]-Coord_1[iDim];
				}

				su2double Ni=0 , Nj=0, Nk=0;

				Ni=AC[1]*BD[2]-AC[2]*BD[1];
				Nj=-AC[0]*BD[2]+AC[2]*BD[0];
				Nk=AC[0]*BD[1]-AC[1]*BD[0];

				Area_Elem = 0.5*sqrt(Ni*Ni+Nj*Nj+Nk*Nk);

			}
		}

      if (nDim == 2) {

        Residual[0] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
        Residual[1] = (1.0/2.0)*Length_Elem*TotalLoad*Load_Dir_Local[1]/Norm;

        node[Point_0]->Add_SurfaceLoad_Res(Residual);
        node[Point_1]->Add_SurfaceLoad_Res(Residual);

      }

      else {
    	  if (geometry->bound[val_marker][iElem]->GetVTK_Type() == TRIANGLE){

    		  Residual[0] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/3.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

    	      node[Point_0]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_1]->Add_SurfaceLoad_Res(Residual);
    	      node[Point_2]->Add_SurfaceLoad_Res(Residual);

    	  }
    	  else if (geometry->bound[val_marker][iElem]->GetVTK_Type() == RECTANGLE){

    		  Residual[0] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[0]/Norm;
    		  Residual[1] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[1]/Norm;
    		  Residual[2] = (1.0/4.0)*Area_Elem*TotalLoad*Load_Dir_Local[2]/Norm;

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
	bool fsi = config->GetFSI_Simulation();												// FSI simulation.

	bool restart = config->GetRestart();													// Restart solution
	bool initial_calc_restart = (config->GetExtIter() == config->GetDyn_RestartIter());	// Restart iteration

	bool incremental_load = config->GetIncrementalLoad();

	if (!dynamic){

		for (iPoint = 0; iPoint < nPointDomain; iPoint++){
			/*--- Add the external contribution to the residual    ---*/
			/*--- (the terms that are constant over the time step) ---*/
			if (incremental_load){
				for (iVar = 0; iVar < nVar; iVar++){
					Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
				}
			}
			else {
				Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			}

			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);
		}

	}

	if (dynamic) {

		/*--- Add the mass matrix contribution to the Jacobian ---*/

		/*
		 * If the problem is nonlinear, we need to add the Mass Matrix contribution to the Jacobian at the beginning
		 * of each time step. If the solution method is Newton Rapshon, we repeat this step at the beginning of each
		 * iteration, as the Jacobian is recomputed
		 *
		 * If the problem is linear, we add the Mass Matrix contribution to the Jacobian at the first calculation.
		 * From then on, the Jacobian is always the same matrix.
		 *
		 */

		if ((nonlinear_analysis && (newton_raphson || first_iter)) ||
			(linear_analysis && initial_calc) ||
			(linear_analysis && restart && initial_calc_restart)) {
			for (iPoint = 0; iPoint < nPointDomain; iPoint++){
				for (jPoint = 0; jPoint < nPointDomain; jPoint++){
					for(iVar = 0; iVar < nVar; iVar++){
						for (jVar = 0; jVar < nVar; jVar++){
							Jacobian_ij[iVar][jVar] = a_dt[0] * MassMatrix.GetBlock(iPoint, jPoint, iVar, jVar);
						}
					}
					Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
				}
			}
		}

		/*--- Loop over all points, and set aux vector TimeRes_Aux = a0*U+a2*U'+a3*U'' ---*/
		if (linear_analysis){
			for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
				for (iVar = 0; iVar < nVar; iVar++){
					Residual[iVar] = a_dt[0]*node[iPoint]->GetSolution_time_n(iVar)+		//a0*U(t)
								 	 a_dt[2]*node[iPoint]->GetSolution_Vel_time_n(iVar)+	//a2*U'(t)
								 	 a_dt[3]*node[iPoint]->GetSolution_Accel_time_n(iVar);	//a3*U''(t)
				}
				TimeRes_Aux.SetBlock(iPoint, Residual);
			}
		}
		else if (nonlinear_analysis){
			for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
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
		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
			/*--- Dynamic contribution ---*/
			/*--- TODO: Do I have to scale this one? I don't think so... ---*/
			Res_Time_Cont = TimeRes.GetBlock(iPoint);
			LinSysRes.AddBlock(iPoint, Res_Time_Cont);
			/*--- External surface load contribution ---*/
			if (incremental_load){
				for (iVar = 0; iVar < nVar; iVar++){
					Res_Ext_Surf[iVar] = loadIncrement * node[iPoint]->Get_SurfaceLoad_Res(iVar);
				}
			}
			else {
				Res_Ext_Surf = node[iPoint]->Get_SurfaceLoad_Res();
			}
			LinSysRes.AddBlock(iPoint, Res_Ext_Surf);

			/*--- Add FSI contribution ---*/
			if (fsi) {
				/*--- TODO: It may be worthy restricting the flow traction to the boundary elements... ---*/
				if (incremental_load){
					for (iVar = 0; iVar < nVar; iVar++){
						Res_FSI_Cont[iVar] = loadIncrement * node[iPoint]->Get_FlowTraction(iVar);
					}
				}
				else {
					Res_FSI_Cont = node[iPoint]->Get_FlowTraction();
				}
				LinSysRes.AddBlock(iPoint, Res_FSI_Cont);
			}
		}
	}


}

void CFEM_ElasticitySolver::ImplicitNewmark_Update(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

    unsigned short iVar;
	unsigned long iPoint, total_index;

	bool linear = (config->GetGeometricConditions() == SMALL_DEFORMATIONS);		// Geometrically linear problems
	bool nonlinear = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems
	bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);							// Dynamic simulations.

	unsigned short iNode, jNode, jVar;

	/*--- Update solution ---*/

	for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

		for (iVar = 0; iVar < nVar; iVar++) {

			/*--- Displacements component of the solution ---*/

			/*--- If it's a non-linear problem, the result is the DELTA_U, not U itself ---*/

			if (linear) node[iPoint]->SetSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

			if (nonlinear)	node[iPoint]->Add_DeltaSolution(iVar, LinSysSol[iPoint*nVar+iVar]);

		}

	}

	if (dynamic){

		for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

			for (iVar = 0; iVar < nVar; iVar++) {

				/*--- Acceleration component of the solution ---*/
				/*--- U''(t+dt) = a0*(U(t+dt)-U(t))+a2*(U'(t))+a3*(U''(t)) ---*/

				Solution[iVar]=a_dt[0]*(node[iPoint]->GetSolution(iVar) -
										node[iPoint]->GetSolution_time_n(iVar)) -
							   a_dt[2]* node[iPoint]->GetSolution_Vel_time_n(iVar) -
							   a_dt[3]* node[iPoint]->GetSolution_Accel_time_n(iVar);
			}

			/*--- Set the acceleration in the node structure ---*/

			node[iPoint]->SetSolution_Accel(Solution);

			for (iVar = 0; iVar < nVar; iVar++) {

				/*--- Velocity component of the solution ---*/
				/*--- U'(t+dt) = U'(t)+ a6*(U''(t)) + a7*(U''(t+dt)) ---*/

				Solution[iVar]=node[iPoint]->GetSolution_Vel_time_n(iVar)+
							   a_dt[6]* node[iPoint]->GetSolution_Accel_time_n(iVar) +
							   a_dt[7]* node[iPoint]->GetSolution_Accel(iVar);

			}

			/*--- Set the velocity in the node structure ---*/

			node[iPoint]->SetSolution_Vel(Solution);

		}

	}

	/*--- Perform the MPI communication of the solution ---*/

	Set_MPI_Solution(geometry, config);


}

void CFEM_ElasticitySolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config){


	unsigned long IterLinSol, iPoint, jPoint, total_index;
	unsigned short iVar;

	/*--- Initialize residual and solution at the ghost points ---*/

	for (iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++) {

		for (iVar = 0; iVar < nVar; iVar++) {
		  total_index = iPoint*nVar + iVar;
		  LinSysRes[total_index] = 0.0;
		  LinSysSol[total_index] = 0.0;
		}

	 }

	CSysSolve femSystem;
	IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

}



void CFEM_ElasticitySolver::SetFEA_Load(CSolver ***flow_solution, CGeometry **fea_geometry,
											CGeometry **flow_geometry, CConfig *fea_config,
											CConfig *flow_config, CNumerics *fea_numerics) {

	unsigned short nMarkerFSI, nMarkerStruct, nMarkerFlow;		// Number of markers on FSI problem, FEA and Flow side
	unsigned short iMarkerFSI, iMarkerStruct, iMarkerFlow;		// Variables for iteration over markers
	int Marker_Flow = -1, Marker_Struct = -1;

	unsigned long nVertexFEA, nVertexFlow;						// Number of vertices on FEA and Flow side
	unsigned long iVertex, iPoint;								// Variables for iteration over vertices and nodes
	unsigned long jVertex;										// Variables for iteration over vertices and nodes

	unsigned short iDim, jDim;

	// Check the kind of fluid problem
	bool compressible       = (flow_config->GetKind_Regime() == COMPRESSIBLE);
	bool incompressible     = (flow_config->GetKind_Regime() == INCOMPRESSIBLE);
	bool viscous_flow       = ((flow_config->GetKind_Solver() == NAVIER_STOKES) ||
								(flow_config->GetKind_Solver() == RANS) );

  	/*--- Redimensionalize the pressure ---*/

	su2double *Velocity_ND, *Velocity_Real;
	su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
	su2double factorForces;

    Velocity_Real = flow_config->GetVelocity_FreeStream();
    Density_Real = flow_config->GetDensity_FreeStream();

    Velocity_ND = flow_config->GetVelocity_FreeStreamND();
    Density_ND = flow_config->GetDensity_FreeStreamND();

	Velocity2_Real = 0.0;
	Velocity2_ND = 0.0;
    for (iDim = 0; iDim < nDim; iDim++){
    	Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    	Velocity2_ND += Velocity_ND[iDim]*Velocity_ND[iDim];
    }

    factorForces = Density_Real*Velocity2_Real/(Density_ND*Velocity2_ND);

	/*--- Apply a ramp to the transfer of the fluid loads ---*/

	su2double ModAmpl;
	su2double CurrentTime = fea_config->GetCurrent_DynTime();
	su2double Static_Time = fea_config->GetStatic_Time();

	bool Ramp_Load = fea_config->GetRamp_Load();
	su2double Ramp_Time = fea_config->GetRamp_Time();

	if (CurrentTime <= Static_Time){ ModAmpl=0.0; }
	else if((CurrentTime > Static_Time) &&
			(CurrentTime <= (Static_Time + Ramp_Time)) &&
			(Ramp_Load)){
		ModAmpl = (CurrentTime-Static_Time) / Ramp_Time;
		ModAmpl = max(ModAmpl,0.0);
		ModAmpl = min(ModAmpl,1.0);
	}
	else{ ModAmpl = 1.0; }

	/*--- Number of markers on the FSI interface ---*/

	nMarkerFSI = (fea_config->GetMarker_n_FSIinterface())/2;

	nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();		// Retrieve total number of markers on FEA side
	nMarkerFlow = flow_geometry[MESH_0]->GetnMarker();		// Retrieve total number of markers on Fluid side

	// Parameters for the calculations
	// Pn: Pressure
	// Pinf: Pressure_infinite
	// div_vel: Velocity divergence
	// Dij: Dirac delta
	su2double Pn = 0.0, Pinf = 0.0, div_vel = 0.0, Dij = 0.0;
	su2double Viscosity = 0.0, Density = 0.0;
	su2double **Grad_PrimVar;
	su2double Tau[3][3];

	unsigned long Point_Flow, Point_Struct;
	su2double *Normal_Flow, *Normal_Struct;

	su2double *tn_f;
	tn_f 				= new su2double [nVar];			// Fluid traction

#ifndef HAVE_MPI

	/*--- TODO: We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/
	/*--- This may be improved ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++){
		node[iPoint]->Clear_FlowTraction();
	}

	/*--- Loop over all the markers on the interface ---*/

	for (iMarkerFSI = 0; iMarkerFSI < nMarkerFSI; iMarkerFSI++){

		/*--- Identification of the markers ---*/

		/*--- Current structural marker ---*/
		for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++){
			if ( fea_config->GetMarker_All_FSIinterface(iMarkerStruct) == (iMarkerFSI+1)){
				Marker_Struct = iMarkerStruct;
			}
		}

		/*--- Current fluid marker ---*/
		for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++){
			if (flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == (iMarkerFSI+1)){
				Marker_Flow = iMarkerFlow;
			}
		}

		nVertexFEA = fea_geometry[MESH_0]->GetnVertex(Marker_Struct);		// Retrieve total number of vertices on FEA marker
		nVertexFlow = flow_geometry[MESH_0]->GetnVertex(Marker_Flow);  // Retrieve total number of vertices on Fluid marker

		/*--- Loop over the nodes in the fluid mesh, calculate the tf vector (unitary) ---*/
		/*--- Here, we are looping over the fluid, and we find the pointer to the structure (Point_Struct) ---*/
		for (iVertex = 0; iVertex < nVertexFlow; iVertex++){

			// Node from the flow mesh
			Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();

			// Normals at the vertex: these normals go inside the fluid domain.
			Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();

			// Corresponding node on the structural mesh
			Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();

			// Retrieve the values of pressure, viscosity and density
			if (incompressible){

				Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressureInc();
				Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosityInc();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetDensityInc();

				}
			}
			else if (compressible){

				Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
				Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

				if (viscous_flow){

					Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
					Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
					Density = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetDensity();

				}
			}

			// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
			for (iDim = 0; iDim < nDim; iDim++) {
				tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
			}

			// Calculate tn in the fluid nodes for the viscous term

			if (viscous_flow){

				// Divergence of the velocity
				div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
				if (incompressible) div_vel = 0.0;

				for (iDim = 0; iDim < nDim; iDim++) {

					for (jDim = 0 ; jDim < nDim; jDim++) {
						// Dirac delta
						Dij = 0.0; if (iDim == jDim) Dij = 1.0;

						// Viscous stress
						Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
								TWO3*Viscosity*div_vel*Dij;

						// Viscous component in the tn vector --> Units of force (non-dimensional).
						tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
					}
				}
			}

			// Rescale tn to SI units and apply time-dependent coefficient (static structure, ramp load, full load)

			for (iDim = 0; iDim < nDim; iDim++) {
				Residual[iDim] = tn_f[iDim]*factorForces*ModAmpl;
			}

			/*--- Set the Flow traction ---*/
			//node[Point_Struct]->Set_FlowTraction(Residual);
			/*--- Add to the Flow traction (to add values to corners...) ---*/
			node[Point_Struct]->Add_FlowTraction(Residual);
		}

	}

#else

    int rank = MASTER_NODE;
    int size = SINGLE_NODE;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    unsigned long nLocalVertexStruct = 0, nLocalVertexFlow = 0;
    unsigned long iVertexFlow = 0;
    unsigned long nPoint_Total = 0;

	unsigned long MaxLocalVertexStruct = 0, MaxLocalVertexFlow = 0;

	unsigned long nBuffer_FlowTraction = 0, nBuffer_StructTraction = 0;
	unsigned long nBuffer_DonorIndices = 0, nBuffer_SetIndex = 0;

	unsigned long Processor_Flow, Processor_Struct;

	int iProcessor, nProcessor = 0;

	/*--- Number of markers on the FSI interface ---*/

	nMarkerFSI     = (flow_config->GetMarker_n_FSIinterface())/2;
	nMarkerStruct  = fea_geometry[MESH_0]->GetnMarker();
	nMarkerFlow    = flow_geometry[MESH_0]->GetnMarker();

	nProcessor = size;

	/*--- TODO: We have to clear the traction before applying it, because we are "adding" to node and not "setting" ---*/
	/*--- This may be improved ---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++){
		node[iPoint]->Clear_FlowTraction();
	}

	/*--- Outer loop over the markers on the FSI interface: compute one by one ---*/
	/*--- The tags are always an integer greater than 1: loop from 1 to nMarkerFSI ---*/

	for (iMarkerFSI = 1; iMarkerFSI <= nMarkerFSI; iMarkerFSI++){

		Marker_Struct = -1;
		Marker_Flow = -1;

		/*--- Initialize pointer buffers inside the loop, so we can delete for each marker. ---*/
		unsigned long Buffer_Send_nVertexStruct[1], *Buffer_Recv_nVertexStruct = NULL;
		unsigned long Buffer_Send_nVertexFlow[1], *Buffer_Recv_nVertexFlow = NULL;

		/*--- The markers on the fluid and structural side are tagged with the same index.
		 *--- This is independent of the MPI domain decomposition.
		 *--- We need to loop over all markers on both sides and get the number of nodes
		 *--- that belong to each FSI marker for each processor ---*/

		/*--- On the structural side ---*/

		for (iMarkerStruct = 0; iMarkerStruct < nMarkerStruct; iMarkerStruct++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerStruct) equals the index we are looping at ---*/
			if ( fea_config->GetMarker_All_FSIinterface(iMarkerStruct) == iMarkerFSI ){
				/*--- We have identified the local index of the FEA marker ---*/
				/*--- Store the number of local points that belong to Marker_Struct on each processor ---*/
				/*--- This includes the halo nodes ---*/
				nLocalVertexStruct = fea_geometry[MESH_0]->GetnVertex(iMarkerStruct);
				/*--- Store the identifier for the structural marker ---*/
				Marker_Struct = iMarkerStruct;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the FEA markers ---*/
				nLocalVertexStruct = 0;
				Marker_Struct = -1;
			}
		}

		/*--- On the fluid side ---*/

		for (iMarkerFlow = 0; iMarkerFlow < nMarkerFlow; iMarkerFlow++){
			/*--- If the tag GetMarker_All_FSIinterface(iMarkerFlow) equals the index we are looping at ---*/
			if ( flow_config->GetMarker_All_FSIinterface(iMarkerFlow) == iMarkerFSI ){
				/*--- We have identified the local index of the Flow marker ---*/
				/*--- Store the number of local points that belong to Marker_Flow on each processor ---*/
				/*--- This includes the halo nodes ---*/
				nLocalVertexFlow = flow_geometry[MESH_0]->GetnVertex(iMarkerFlow);
				/*--- Store the identifier for the fluid marker ---*/
				Marker_Flow = iMarkerFlow;
				/*--- Exit the for loop: we have found the local index for iMarkerFSI on the FEA side ---*/
				break;
			}
			else {
				/*--- If the tag hasn't matched any tag within the Flow markers ---*/
				nLocalVertexFlow = 0;
				Marker_Flow = -1;
			}
		}

		Buffer_Send_nVertexStruct[0] = nLocalVertexStruct;							    // Retrieve total number of vertices on FEA marker
		Buffer_Send_nVertexFlow[0] = nLocalVertexFlow;								    // Retrieve total number of vertices on Flow marker
		if (rank == MASTER_NODE) Buffer_Recv_nVertexStruct = new unsigned long[size];   // Allocate memory to receive how many vertices are on each rank on the structural side
		if (rank == MASTER_NODE) Buffer_Recv_nVertexFlow = new unsigned long[size];     // Allocate memory to receive how many vertices are on each rank on the fluid side

		/*--- We receive MaxLocalVertexFEA as the maximum number of vertices in one single processor on the structural side---*/
		SU2_MPI::Allreduce(&nLocalVertexStruct, &MaxLocalVertexStruct, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
		/*--- We receive MaxLocalVertexFlow as the maximum number of vertices in one single processor on the fluid side ---*/
		SU2_MPI::Allreduce(&nLocalVertexFlow, &MaxLocalVertexFlow, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);

		/*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the structural side ---*/
		SU2_MPI::Gather(&Buffer_Send_nVertexStruct, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexStruct, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);
		/*--- We gather a vector in MASTER_NODE that determines how many elements are there on each processor on the fluid side ---*/
		SU2_MPI::Gather(&Buffer_Send_nVertexFlow, 1, MPI_UNSIGNED_LONG, Buffer_Recv_nVertexFlow, 1, MPI_UNSIGNED_LONG, MASTER_NODE, MPI_COMM_WORLD);

		/*--- We will be gathering the structural coordinates into the master node ---*/
		/*--- Then we will distribute them using a scatter operation into the appropriate fluid processor ---*/
		nBuffer_FlowTraction = MaxLocalVertexFlow * nDim;
		nBuffer_StructTraction = MaxLocalVertexStruct * nDim;

		/*--- We will be gathering donor index and donor processor (for flow -> donor = structure) ---*/
		/*--- Then we will pass on to the structural side the index (fea point) to the appropriate processor ---*/
		nBuffer_DonorIndices = 2 * MaxLocalVertexFlow;
		nBuffer_SetIndex = MaxLocalVertexStruct;

		/*--- Send and Recv buffers ---*/

		/*--- Buffers to send and receive the structural coordinates ---*/
		su2double *Buffer_Send_FlowTraction = new su2double[nBuffer_FlowTraction];
		su2double *Buffer_Recv_FlowTraction = NULL;

		/*--- Buffers to send and receive the donor index and processor ---*/
		long *Buffer_Send_DonorIndices = new long[nBuffer_DonorIndices];
		long *Buffer_Recv_DonorIndices = NULL;

		/*--- Buffers to send and receive the new fluid coordinates ---*/
		su2double *Buffer_Send_StructTraction = NULL;
		su2double *Buffer_Recv_StructTraction = new su2double[nBuffer_StructTraction];

		/*--- Buffers to send and receive the fluid index ---*/
		long *Buffer_Send_SetIndex = NULL;
		long *Buffer_Recv_SetIndex = new long[nBuffer_SetIndex];

		/*--- Prepare the receive buffers (1st step) and send buffers (2nd step) on the master node only. ---*/

		if (rank == MASTER_NODE) {
			Buffer_Recv_FlowTraction  = new su2double[size*nBuffer_FlowTraction];
			Buffer_Recv_DonorIndices = new long[size*nBuffer_DonorIndices];
			Buffer_Send_StructTraction = new su2double[size*nBuffer_StructTraction];
			Buffer_Send_SetIndex     = new long[size*nBuffer_SetIndex];
		}

		/*--- On the fluid side ---*/

		/*--- If this processor owns the marker we are looping at on the structural side ---*/

		/*--- First we initialize all of the indices and processors to -1 ---*/
		/*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
		for (iVertex = 0; iVertex < nBuffer_DonorIndices; iVertex++)
			Buffer_Send_DonorIndices[iVertex] = -1;

		if (Marker_Flow >= 0){

			/*--- We have identified the local index of the FEA marker ---*/
			/*--- We loop over all the vertices in that marker and in that particular processor ---*/

			for (iVertex = 0; iVertex < nLocalVertexFlow; iVertex++){

		        Point_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNode();

		        Point_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorPoint();

		        Processor_Struct = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetDonorProcessor();

				// Get the normal at the vertex: this normal goes inside the fluid domain.
				Normal_Flow = flow_geometry[MESH_0]->vertex[Marker_Flow][iVertex]->GetNormal();

				// Retrieve the values of pressure, viscosity and density
				if (incompressible){

					Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressureInc();
					// TODO: This can be taken out (one thing at a time!)
					Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

					if (viscous_flow){

						Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
						Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosityInc();
						Density = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetDensityInc();

					}
				}
				else if (compressible){

					Pn = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetPressure();
					// TODO: This can be taken out (one thing at a time!)
					Pinf = flow_solution[MESH_0][FLOW_SOL]->GetPressure_Inf();

					if (viscous_flow){

						Grad_PrimVar = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetGradient_Primitive();
						Viscosity = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetLaminarViscosity();
						Density = flow_solution[MESH_0][FLOW_SOL]->node[Point_Flow]->GetDensity();

					}
				}

				// Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
				for (iDim = 0; iDim < nDim; iDim++) {
					tn_f[iDim] = -(Pn-Pinf)*Normal_Flow[iDim];
				}

				// Calculate tn in the fluid nodes for the viscous term

				if (viscous_flow){

					// Divergence of the velocity
					div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_PrimVar[iDim+1][iDim];
					if (incompressible) div_vel = 0.0;

					for (iDim = 0; iDim < nDim; iDim++) {

						for (jDim = 0 ; jDim < nDim; jDim++) {
							// Dirac delta
							Dij = 0.0; if (iDim == jDim) Dij = 1.0;

							// Viscous stress
							Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[jDim+1][iDim] + Grad_PrimVar[iDim+1][jDim]) -
									TWO3*Viscosity*div_vel*Dij;

							// Viscous component in the tn vector --> Units of force (non-dimensional).
							tn_f[iDim] += Tau[iDim][jDim]*Normal_Flow[jDim];
						}
					}
				}

				for (iDim = 0; iDim < nDim; iDim++){
					Buffer_Send_FlowTraction[iVertex*nDim+iDim] = tn_f[iDim]*factorForces*ModAmpl;
				}
				/*--- If this processor owns the node ---*/
				if (flow_geometry[MESH_0]->node[Point_Flow]->GetDomain()){
					Buffer_Send_DonorIndices[2*iVertex]     = Point_Struct;
					Buffer_Send_DonorIndices[2*iVertex + 1] = Processor_Struct;
				}
				else{
					/*--- We set the values to be -1 to be able to identify them later as halo nodes ---*/
					Buffer_Send_DonorIndices[2*iVertex]     = -1;
					Buffer_Send_DonorIndices[2*iVertex + 1] = -1;
				}

			}
		}

		/*--- Once all the messages have been sent, we gather them all into the MASTER_NODE ---*/
		SU2_MPI::Gather(Buffer_Send_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, Buffer_Recv_FlowTraction, nBuffer_FlowTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Gather(Buffer_Send_DonorIndices, nBuffer_DonorIndices, MPI_LONG, Buffer_Recv_DonorIndices, nBuffer_DonorIndices, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

//		if (rank == MASTER_NODE){
//			cout << endl << "-----------------------------------------------------------" << endl;
//			cout << "For tag " << iMarkerFSI << ":" << endl;
//			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
//				cout << "The processor " << iProcessor << " has " << Buffer_Recv_nVertexStruct[iProcessor] << " nodes on the structural side and ";
//				cout << Buffer_Recv_nVertexFlow[iProcessor] << " nodes on the fluid side " << endl;
//			}
//			cout << "The max number of vertices is " << MaxLocalVertexStruct << " on the structural side and ";
//			cout << MaxLocalVertexFlow << " on the fluid side." << endl;
//
//			cout << "---------------- Check received buffers ---------------------" << endl;
//			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
//				long initialIndex, initialIndex2;
//				initialIndex = iProcessor*nBuffer_FlowTraction;
//				initialIndex2 = iProcessor*nBuffer_DonorIndices;
//				for (long iCheck = 0; iCheck < Buffer_Recv_nVertexStruct[iProcessor]; iCheck++){
//					cout << "From processor " << iProcessor << " we get coordinates (";
//						for (iDim = 0; iDim < nDim; iDim++)
//							cout << Buffer_Recv_FlowTraction[initialIndex+iCheck*nDim+iDim] << ",";
//					cout << "), the donor index for the flow " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2] ;
//					cout << " and the donor processor " << Buffer_Recv_DonorIndices[initialIndex2+iCheck*2+1] << endl;
//
//				}
//			}
//
//		}

		/*--- Counter to determine where in the array we have to set the information ---*/
		long *Counter_Processor_Struct = NULL;
		long iProcessor_Flow = 0, iIndex_Flow = 0;
		long iProcessor_Struct = 0, iPoint_Struct = 0, iIndex_Struct = 0;

		/*--- Now we pack the information to send it over to the different processors ---*/

		if (rank == MASTER_NODE){

			/*--- We set the counter to 0 ---*/
			Counter_Processor_Struct = new long[nProcessor];
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
				Counter_Processor_Struct[iProcessor] = 0;
			}

			/*--- First we initialize the index vector to -1 ---*/
			/*--- This helps on identifying halo nodes and avoids setting wrong values ---*/
			for (iVertex = 0; iVertex < nProcessor*nBuffer_SetIndex; iVertex++)
				Buffer_Send_SetIndex[iVertex] = -2;

			/*--- As of now we do the loop over the flow points ---*/
			/*--- The number of points for flow and structure does not necessarily have to match ---*/
			/*--- In fact, it's possible that a processor asks for nStruct nodes and there are only ---*/
			/*--- nFlow < nStruct available; this is due to halo nodes ---*/

			/*--- For every processor from which we have received information ---*/
			/*--- (This is, for every processor on the structural side) ---*/
			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){

				/*--- This is the initial index on the coordinates buffer for that particular processor on the structural side ---*/
				iProcessor_Flow = iProcessor*nBuffer_FlowTraction;
				/*--- This is the initial index on the donor index/processor buffer for that particular processor on the structural side ---*/
				iIndex_Flow = iProcessor*nBuffer_DonorIndices;

				/*--- For every vertex in the information retreived from iProcessor ---*/
				for (iVertex = 0; iVertex < Buffer_Recv_nVertexFlow[iProcessor]; iVertex++) {

					/*--- The processor and index for the flow are: ---*/
					Processor_Struct = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2+1];
					Point_Struct     = Buffer_Recv_DonorIndices[iIndex_Flow+iVertex*2];

					/*--- Load the buffer at the appropriate position ---*/
					/*--- This is determined on the fluid side by:
					 *--- Processor_Flow*nBuffer_StructTraction -> Initial position of the processor array (fluid side)
					 *--- +
					 *--- Counter_Processor_Struct*nDim -> Initial position of the nDim array for the particular point on the fluid side
					 *--- +
					 *--- iDim -> Position within the nDim array that corresponds to a point
					 *---
					 *--- While on the structural side is:
					 *--- iProcessor*nBuffer_FlowTraction -> Initial position on the processor array (structural side)
					 *--- +
					 *--- iVertex*nDim -> Initial position of the nDim array for the particular point on the structural side
					 */

					/*--- We check that we are not setting the value for a halo node ---*/
					if (Point_Struct != -1){
						iProcessor_Struct = Processor_Struct*nBuffer_StructTraction;
						iIndex_Struct = Processor_Struct*nBuffer_SetIndex;
						iPoint_Struct = Counter_Processor_Struct[Processor_Struct]*nDim;

						for (iDim = 0; iDim < nDim; iDim++)
							Buffer_Send_StructTraction[iProcessor_Struct + iPoint_Struct + iDim] = Buffer_Recv_FlowTraction[iProcessor_Flow + iVertex*nDim + iDim];

						/*--- We set the fluid index at an appropriate position matching the coordinates ---*/
						Buffer_Send_SetIndex[iIndex_Struct + Counter_Processor_Struct[Processor_Struct]] = Point_Struct;

						Counter_Processor_Struct[Processor_Struct]++;
					}

				}

			}

//			cout << "---------------- Check send buffers ---------------------" << endl;
//
//			for (iProcessor = 0; iProcessor < nProcessor; iProcessor++){
//				long initialIndex, initialIndex2;
//				initialIndex = iProcessor*nBuffer_StructTraction;
//				initialIndex2 = iProcessor*nBuffer_SetIndex;
//				for (long iCheck = 0; iCheck < Buffer_Recv_nVertexFlow[iProcessor]; iCheck++){
//					cout << "Processor " << iProcessor << " will receive the node " ;
//					cout << Buffer_Send_SetIndex[initialIndex2+iCheck] << " which corresponds to the coordinates ";
//					for (iDim = 0; iDim < nDim; iDim++)
//						cout << "x" << iDim << "=" << Buffer_Send_StructTraction[initialIndex + iCheck*nDim + iDim] << ", ";
//					cout << endl;
//				}
//
//			}

		}

		/*--- Once all the messages have been prepared, we scatter them all from the MASTER_NODE ---*/
		SU2_MPI::Scatter(Buffer_Send_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, Buffer_Recv_StructTraction, nBuffer_StructTraction, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Scatter(Buffer_Send_SetIndex, nBuffer_SetIndex, MPI_LONG, Buffer_Recv_SetIndex, nBuffer_SetIndex, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD);

		long indexPoint_iVertex, Point_Struct_Check;

		/*--- For the flow marker we are studying ---*/
		if (Marker_Struct >= 0){

			/*--- We have identified the local index of the Structural marker ---*/
			/*--- We loop over all the vertices in that marker and in that particular processor ---*/

			for (iVertex = 0; iVertex < nLocalVertexStruct; iVertex++){

				Point_Struct = fea_geometry[MESH_0]->vertex[Marker_Struct][iVertex]->GetNode();

				if (fea_geometry[MESH_0]->node[Point_Struct]->GetDomain()){
					/*--- Find the index of the point Point_Struct in the buffer Buffer_Recv_SetIndex ---*/
					indexPoint_iVertex = std::distance(Buffer_Recv_SetIndex, std::find(Buffer_Recv_SetIndex, Buffer_Recv_SetIndex + MaxLocalVertexStruct, Point_Struct));

					Point_Struct_Check = Buffer_Recv_SetIndex[indexPoint_iVertex];

					if (Point_Struct_Check < 0) {
						cout << "WARNING: A nonphysical point is being considered for traction transfer." << endl;
						exit(EXIT_FAILURE);
					}

					for (iDim = 0; iDim < nDim; iDim++)
						Residual[iDim] = Buffer_Recv_StructTraction[indexPoint_iVertex*nDim+iDim];

					/*--- Add to the Flow traction ---*/
					node[Point_Struct]->Add_FlowTraction(Residual);

				}

			}

		}

		delete [] Buffer_Send_FlowTraction;
		delete [] Buffer_Send_DonorIndices;
		delete [] Buffer_Recv_StructTraction;
		delete [] Buffer_Recv_SetIndex;

		if (rank == MASTER_NODE) {
			delete [] Buffer_Recv_nVertexStruct;
			delete [] Buffer_Recv_nVertexFlow;
			delete [] Buffer_Recv_FlowTraction;
			delete [] Buffer_Recv_DonorIndices;
			delete [] Buffer_Send_StructTraction;
			delete [] Buffer_Send_SetIndex;
			delete [] Counter_Processor_Struct;
		}

	}

#endif

	delete[] tn_f;


}

void CFEM_ElasticitySolver::SetFEA_Load_Int(CSolver ***flow_solution, CGeometry **fea_geometry,
											CGeometry **flow_geometry, CConfig *fea_config,
											CConfig *flow_config, CNumerics *fea_numerics){ }

void CFEM_ElasticitySolver::PredictStruct_Displacement(CGeometry **fea_geometry,
                            				CConfig *fea_config, CSolver ***fea_solution){

    unsigned short predOrder = fea_config->GetPredictorOrder();
	su2double Delta_t = fea_config->GetDelta_DynTime();
    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    su2double *solDisp, *solVel, *solVel_tn, *valPred;
    su2double *DisplacementDonor, *SolutionDonor;

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    //To nPointDomain: we need to communicate the predicted solution after setting it
    for (iPoint=0; iPoint < nPointDomain; iPoint++){
    	if (predOrder==0) fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    	else if (predOrder==1) {

    		solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
    		solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
    		valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

    		for (iDim=0; iDim<nDim; iDim++){
    			valPred[iDim] = solDisp[iDim] + Delta_t*solVel[iDim];
    		}

    	}
    	else if (predOrder==2) {

    		solDisp = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
    		solVel = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel();
    		solVel_tn = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Vel_time_n();
    		valPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

    		for (iDim=0; iDim<nDim; iDim++){
    			valPred[iDim] = solDisp[iDim] + 0.5*Delta_t*(3*solVel[iDim]-solVel_tn[iDim]);
    		}

    	}
    	else {
    		cout<< "Higher order predictor not implemented. Solving with order 0." << endl;
    		fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred();
    	}
    }


}

void CFEM_ElasticitySolver::ComputeAitken_Coefficient(CGeometry **fea_geometry, CConfig *fea_config,
        				  CSolver ***fea_solution, unsigned long iFSIIter){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    su2double *dispPred, *dispCalc, *dispPred_Old, *dispCalc_Old;
    su2double deltaU[3] = {0.0, 0.0, 0.0}, deltaU_p1[3] = {0.0, 0.0, 0.0};
    su2double delta_deltaU[3] = {0.0, 0.0, 0.0};
    su2double numAitk, denAitk, WAitken;
	su2double CurrentTime=fea_config->GetCurrent_DynTime();
	su2double Static_Time=fea_config->GetStatic_Time();
	su2double WAitkDyn_tn1, WAitkDyn_Max, WAitkDyn;

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    WAitken=fea_config->GetAitkenStatRelax();

	numAitk = 0.0;
	denAitk = 0.0;

	ofstream historyFile_FSI;
	bool writeHistFSI = fea_config->GetWrite_Conv_FSI();
	if (writeHistFSI){
		char cstrFSI[200];
		string filenameHistFSI = fea_config->GetConv_FileName_FSI();
		strcpy (cstrFSI, filenameHistFSI.data());
		historyFile_FSI.open (cstrFSI, std::ios_base::app);
	}


	/*--- Only when there is movement, and a dynamic coefficient is requested, it makes sense to compute the Aitken's coefficient ---*/

	if (CurrentTime > Static_Time) {

		if (iFSIIter == 0){

			WAitkDyn_tn1 = GetWAitken_Dyn_tn1();
			WAitkDyn_Max = fea_config->GetAitkenDynMaxInit();

			WAitkDyn = min(WAitkDyn_tn1, WAitkDyn_Max);

			/*--- Temporal fix, only for now ---*/
			WAitkDyn = max(WAitkDyn, 0.1);

			SetWAitken_Dyn(WAitkDyn);
			if (writeHistFSI){
				historyFile_FSI << " " << endl ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
				historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
				historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn ;
			}

		}
		else{
		    // To nPointDomain; we need to have communicated Solution, Solution_Old, Solution_Pred and Solution_Pred_Old beforehand
			for (iPoint = 0; iPoint < nPointDomain; iPoint++){

				dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
				dispPred_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred_Old();
				dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();
				dispCalc_Old = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Old();

				for (iDim = 0; iDim < nDim; iDim++){

					/*--- Compute the deltaU and deltaU_n+1 ---*/
					deltaU[iDim] = dispCalc_Old[iDim] - dispPred_Old[iDim];
					deltaU_p1[iDim] = dispCalc[iDim] - dispPred[iDim];

					/*--- Compute the difference ---*/
					delta_deltaU[iDim] = deltaU_p1[iDim] - deltaU[iDim];

					/*--- Add numerator and denominator ---*/
					numAitk += deltaU[iDim] * delta_deltaU[iDim];
					denAitk += delta_deltaU[iDim] * delta_deltaU[iDim];

				}

			}

				WAitkDyn = GetWAitken_Dyn();

			if (denAitk > 1E-8){
				WAitkDyn = - 1.0 * WAitkDyn * numAitk / denAitk ;
			}

				WAitkDyn = max(WAitkDyn, 0.1);
				WAitkDyn = min(WAitkDyn, 1.0);

				SetWAitken_Dyn(WAitkDyn);

				if (writeHistFSI){
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(4) << CurrentTime << "," ;
					historyFile_FSI << setiosflags(ios::fixed) << setprecision(1) << iFSIIter << "," ;
					historyFile_FSI << setiosflags(ios::scientific) << setprecision(4) << WAitkDyn << "," ;
				}

		}

	}

	if (writeHistFSI){historyFile_FSI.close();}

}

void CFEM_ElasticitySolver::SetAitken_Relaxation(CGeometry **fea_geometry,
        				  CConfig *fea_config, CSolver ***fea_solution){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    unsigned short RelaxMethod_FSI;
    su2double *dispPred, *dispCalc;
    su2double WAitken;
	su2double CurrentTime=fea_config->GetCurrent_DynTime();
	su2double Static_Time=fea_config->GetStatic_Time();

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();

    RelaxMethod_FSI = fea_config->GetRelaxation_Method_FSI();

	/*--- Only when there is movement it makes sense to update the solutions... ---*/

	if (CurrentTime > Static_Time) {

		if (RelaxMethod_FSI == NO_RELAXATION){
			WAitken = 1.0;
		}
		else if (RelaxMethod_FSI == FIXED_PARAMETER){
			WAitken = fea_config->GetAitkenStatRelax();
		}
		else if (RelaxMethod_FSI == AITKEN_DYNAMIC){
			WAitken = GetWAitken_Dyn();
		}
		else {
			WAitken = 1.0;
			cout << "No relaxation parameter used. " << endl;
		}

	    // To nPointDomain; we need to communicate the solutions (predicted, old and old predicted) after this routine
		for (iPoint=0; iPoint < nPointDomain; iPoint++){

			/*--- Retrieve pointers to the predicted and calculated solutions ---*/
			dispPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();
			dispCalc = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution();

			/*--- Set predicted solution as the old predicted solution ---*/
			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Pred_Old();

			/*--- Set calculated solution as the old solution (needed for dynamic Aitken relaxation) ---*/
			fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution_Old(dispCalc);

			/*--- Apply the Aitken relaxation ---*/
			for (iDim=0; iDim < nDim; iDim++){
				dispPred[iDim] = (1.0 - WAitken)*dispPred[iDim] + WAitken*dispCalc[iDim];
			}

		}

	}

}

void CFEM_ElasticitySolver::Update_StructSolution(CGeometry **fea_geometry,
        				  CConfig *fea_config, CSolver ***fea_solution){

    unsigned long iPoint, iDim;
    unsigned long nPoint, nDim;
    su2double *valSolutionPred, *valSolution;

    /*--- TODO: I don't think I'm using this routine right now, but I may need to use it after SetAitken_Relaxation... ---*/

    nPoint = fea_geometry[MESH_0]->GetnPoint();
    nDim = fea_geometry[MESH_0]->GetnDim();
    //TODO: nPoint or nPointDomain
    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	valSolutionPred = fea_solution[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_Pred();

		fea_solution[MESH_0][FEA_SOL]->node[iPoint]->SetSolution(valSolutionPred);

    }

}

