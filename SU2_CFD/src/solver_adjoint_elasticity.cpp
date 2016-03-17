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

CFEM_ElasticitySolver_Adj::CFEM_ElasticitySolver_Adj(void) : CSolver() {

	nElement   = 0;
	nMarker    = 0;
	nFEA_Terms = 0;

	direct_solver 		= NULL;
	GradN_X				= NULL;
	GradN_x				= NULL;

	Jacobian_c_ij		= NULL;
	Jacobian_s_ij		= NULL;

	mZeros_Aux			= NULL;
	mId_Aux				= NULL;

	element_container	= NULL;

}

CFEM_ElasticitySolver_Adj::CFEM_ElasticitySolver_Adj(CGeometry *geometry, CConfig *config, CSolver *direct_solver) : CSolver() {

	unsigned long iPoint;
	unsigned short iVar, jVar, iDim, jDim;
	unsigned short iTerm;

	unsigned short iZone = config->GetiZone();
	unsigned short nZone = geometry->GetnZone();

	su2double E = config->GetElasticyMod();

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif

	bool de_effects = config->GetDE_Effects();											// Test whether we consider dielectric elastomers
	bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Nonlinear analysis.

	/*--- Store a pointer to the direct solver: this is where the solution is stored ---*/
	this->direct_solver = direct_solver;

	nElement      = geometry->GetnElem();
	nDim          = geometry->GetnDim();
	nMarker       = geometry->GetnMarker();

	nPoint        = geometry->GetnPoint();
	nPointDomain  = geometry->GetnPointDomain();

	/*--- Number of different terms for FEA ---*/
	nFEA_Terms = 1;
	if (de_effects) nFEA_Terms++;
	/*--- First level: different possible terms of the equations ---*/
	element_container = new CElement** [MAX_TERMS];
	for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
		element_container[iTerm] = new CElement* [MAX_FE_KINDS];

	if (nDim == 2){
		element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
		element_container[FEA_TERM][EL_QUAD] = new CQUAD4(nDim, config);

		if (de_effects){
			element_container[DE_TERM][EL_TRIA] = new CTRIA1(nDim, config);
			element_container[DE_TERM][EL_QUAD] = new CQUAD1(nDim, config);
		}

	}
	else if (nDim == 3){
		element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
		element_container[FEA_TERM][EL_HEXA] = new CHEXA8(nDim, config);
	}
	GradN_X = new su2double [nDim];
	GradN_x = new su2double [nDim];

	nVar = nDim;

	/*--- Define some auxiliary vectors ---*/
	Residual = new su2double[nVar];		for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]	= 0.0;
	Solution = new su2double[nVar]; 	for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]	= 0.0;

	/*--- Term ij of the Jacobian ---*/

	Jacobian_ij = new su2double*[nVar];
	for (iVar = 0; iVar < nVar; iVar++) {
		Jacobian_ij[iVar] = new su2double [nVar];
		for (jVar = 0; jVar < nVar; jVar++) {
			Jacobian_ij[iVar][jVar] = 0.0;
		}
	}
	if (nonlinear_analysis){

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

	}
	else{
		Jacobian_c_ij = NULL;
		Jacobian_s_ij = NULL;
	}
	/*--- Matrices to impose clamped boundary conditions ---*/

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
		mId_Aux[iDim][iDim] = E;
	}

	unsigned short nSolVar;
	unsigned long index;
	string text_line, filename;
	ifstream restart_file;
	su2double dull_val;
	long Dyn_RestartIter;
	bool restart = (config->GetRestart() || config->GetRestart_Flow());

	/*--- Check for a restart, initialize from zero otherwise ---*/
	node = new CVariable*[nPoint];
	if (!restart) {
		for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
		for (iPoint = 0; iPoint < nPoint; iPoint++) {
			node[iPoint] = new CFEM_ElasVariable_Adj(Solution, nDim, nVar, config);
		}
	}
	else {

		/*--- Restart the solution from file information ---*/
		filename = config->GetSolution_FEMFileName();

		/*--- If multizone, append zone name ---*/
	    if (nZone > 1)
	      filename = config->GetMultizone_FileName(filename, iZone);

		restart_file.open(filename.data(), ios::in);

		/*--- In case there is no file ---*/

		if (restart_file.fail()) {
		    if (rank == MASTER_NODE)
			cout << "There is no FEM adjoint restart file!!" << endl;
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

				if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1];
				if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];

				node[iPoint_Local] = new CFEM_ElasVariable(Solution, nDim, nVar, config);
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
			node[iPoint] = new CFEM_ElasVariable_Adj(Solution, nDim, nVar, config);
		}


		/*--- Close the restart file ---*/

		restart_file.close();

		/*--- Free memory needed for the transformation ---*/

		delete [] Global2Local;

	}


	/*--- Initialize the different structures for the linear system ---*/
//	Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);
	LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);

	LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);				// Rv = dF/dv - dK/dV*x

	/*--- Initialize the different structures for the linear system ---*/
	LinSysSol_Direct.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes_dSdv.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes_FSens_DV.Initialize(nPoint, nPointDomain, nVar, 0.0);
	LinSysRes_Aux.Initialize(nPoint, nPointDomain, nVar, 0.0);

	Jacobian_ISens.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

	/*--- Check whether we need to read a reference geometry ---*/
	bool reference_geometry = config->GetRefGeom();
	/*--- Read the reference geometry file ---*/
	if (reference_geometry) Set_ReferenceGeometry(geometry, config);

	/*--- Perform the MPI communication of the solution ---*/
	Set_MPI_Solution(geometry, config);
}

CFEM_ElasticitySolver_Adj::~CFEM_ElasticitySolver_Adj(void) {

}

void CFEM_ElasticitySolver_Adj::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {

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
	      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

	      /*--- Allocate Receive and send buffers  ---*/
	      Buffer_Receive_U = new su2double [nBufferR_Vector];
	      Buffer_Send_U = new su2double[nBufferS_Vector];

	      /*--- Copy the solution that should be sended ---*/
	      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
	        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
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

	        /*--- Copy conserved variables before performing transformation. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

	        /*--- Copy solution variables back into buffer. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          node[iPoint]->SetSolution(iVar, Solution[iVar]);

	      }

	      /*--- Deallocate receive buffer ---*/
	      delete [] Buffer_Receive_U;

	    }

	  }

}

void CFEM_ElasticitySolver_Adj::Set_MPI_RefGeom(CGeometry *geometry, CConfig *config) {

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
	      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

	      /*--- Allocate Receive and send buffers  ---*/
	      Buffer_Receive_U = new su2double [nBufferR_Vector];
	      Buffer_Send_U = new su2double[nBufferS_Vector];

	      /*--- Copy the solution that should be sended ---*/
	      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
	        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
	        for (iVar = 0; iVar < nVar; iVar++)
	          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetReference_Geometry(iVar);
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

	        /*--- Copy conserved variables before performing transformation. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

	        /*--- Copy solution variables back into buffer. ---*/
	        for (iVar = 0; iVar < nVar; iVar++)
	          node[iPoint]->SetReference_Geometry(iVar, Solution[iVar]);

	      }

	      /*--- Deallocate receive buffer ---*/
	      delete [] Buffer_Receive_U;

	    }

	  }

}

void CFEM_ElasticitySolver_Adj::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics,
		unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {

	cout << "PREPROCESSING ADJOINT" << endl;

}

void CFEM_ElasticitySolver_Adj::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter){ }

void CFEM_ElasticitySolver_Adj::Set_ReferenceGeometry(CGeometry *geometry, CConfig *config){

	unsigned long iPoint;
	unsigned long index;

	unsigned short iVar;
	unsigned short iZone = config->GetiZone();
	unsigned short nZone = geometry->GetnZone();
	unsigned short file_format = config->GetRefGeom_FileFormat();

	string filename;
	su2double dull_val;
	ifstream reference_file;

	int rank = MASTER_NODE;
	#ifdef HAVE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#endif


	/*--- Restart the solution from file information ---*/

	filename = config->GetRefGeom_FEMFileName();

	/*--- If multizone, append zone name ---*/
    if (nZone > 1)
      filename = config->GetMultizone_FileName(filename, iZone);

	reference_file.open(filename.data(), ios::in);

	/*--- In case there is no file ---*/

	if (reference_file.fail()) {
	    if (rank == MASTER_NODE)
		cout << "There is no FEM reference geometry file!!" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Filename: " << filename << " and format " << file_format << "." << endl;

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

	getline (reference_file, text_line);

	while (getline (reference_file, text_line)) {
	   istringstream point_line(text_line);

	/*--- Retrieve local index. If this node from the restart file lives
   	   on a different processor, the value of iPoint_Local will be -1.
   	   Otherwise, the local index for this node on the current processor
   	   will be returned and used to instantiate the vars. ---*/

	   iPoint_Local = Global2Local[iPoint_Global];

		if (iPoint_Local >= 0) {

			if (nDim == 2) point_line >> index >> dull_val >> dull_val >> Solution[0] >> Solution[1];
			if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> Solution[0] >> Solution[1] >> Solution[2];

			for (iVar = 0; iVar < nVar; iVar++) node[iPoint_Local]->SetReference_Geometry(iVar, Solution[iVar]);

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

	/*--- Close the restart file ---*/

	reference_file.close();

	/*--- We need to communicate the reference geometry for the halo nodes. ---*/
    Set_MPI_RefGeom(geometry, config);

	/*--- Free memory needed for the transformation ---*/

	delete [] Global2Local;

}

void CFEM_ElasticitySolver_Adj::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	unsigned long iElem, iVar, jVar, iPoint;
	unsigned short iNode,  iDim, nNodes;
	unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
	su2double val_Coord;
	int EL_KIND;

	su2double *Kab = NULL;
	unsigned short NelNodes, jNode;

	/*--- For the solution of the adjoint problem, this routine computes dK/dv, being v the design variables ---*/

	/*--- TODO: This needs to be extended to n Design Variables. As of now, only for E. ---*/

	/*--- Set the value of the Jacobian and Auxiliary vectors to 0 ---*/
	Jacobian_ISens.SetValZero();
	LinSysRes_Aux.SetValZero();
	LinSysRes_dSdv.SetValZero();
	for (iVar = 0; iVar < nVar; iVar++){
		Residual[iVar] = 0.0;
	}

	/*--- Initialize the value of the vector -x, being x the solution vector of the direct solver---*/
	for (iPoint = 0; iPoint < nPoint; iPoint++){
		for (iVar = 0; iVar < nVar; iVar++){
			Residual[iVar] = direct_solver->node[iPoint]->GetSolution(iVar);
		}
		LinSysRes_Aux.SubtractBlock(iPoint, Residual);
	}

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
			  element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
		  }
		}

		numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);

		NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

		for (iNode = 0; iNode < NelNodes; iNode++){

			for (jNode = 0; jNode < NelNodes; jNode++){

				Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);

				for (iVar = 0; iVar < nVar; iVar++){
					for (jVar = 0; jVar < nVar; jVar++){
						Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
					}
				}

				Jacobian_ISens.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);

			}

		}

	}

	/*--- Once computed, compute Rv = dK/dE*(-x) =  ---*/
	Jacobian_ISens.MatrixVectorProduct(LinSysRes_Aux,LinSysRes_dSdv,geometry,config);

}

void CFEM_ElasticitySolver_Adj::Compute_StiffMatrix_NodalStressRes(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){ }

void CFEM_ElasticitySolver_Adj::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config){ }

void CFEM_ElasticitySolver_Adj::BC_Clamped(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker){ }

void CFEM_ElasticitySolver_Adj::BC_Clamped_Post(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                     unsigned short val_marker){ }

void CFEM_ElasticitySolver_Adj::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config){ }

void CFEM_ElasticitySolver_Adj::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
		unsigned short iMesh){ }

void CFEM_ElasticitySolver_Adj::RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CConfig *config){

	unsigned short iVar;
	unsigned long iPoint;
	su2double *reference_geometry, *current_solution;

	su2double objective_function = 0.0;

    for (iPoint = 0; iPoint < nPoint; iPoint++){

    	/*--- The reference geometry is stored in the adjoint variable ---*/
    	reference_geometry = node[iPoint]->GetReference_Geometry();

    	/*--- The current solution comes from the direct solver ---*/
    	current_solution = direct_solver->node[iPoint]->GetSolution();

    	for (iVar = 0; iVar < nVar; iVar++){

    		Solution[iVar] = 2*(current_solution[iVar]-reference_geometry[iVar]);
    		objective_function += (current_solution[iVar]-reference_geometry[iVar])*(current_solution[iVar]-reference_geometry[iVar]);

    	}

    	LinSysRes.AddBlock(iPoint, Solution);

    }

    val_I = objective_function;


}

void CFEM_ElasticitySolver_Adj::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config){

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

	/*--- Solve system of K*dx/dE=Rv ---*/
	/*--- LinSysSol_Direct is dx/dE ---*/
	/*--- LinSysRes_dSdv is Rv ---*/

	CSysSolve femSystem;
	IterLinSol = femSystem.Solve(direct_solver->Jacobian, LinSysRes_dSdv, LinSysSol_Direct, geometry, config);

	/*--- Solve system of K*PHI=dI/dx ---*/
	/*--- LinSysSol is PHI ---*/
	/*--- LinSysRes is dI/dx ---*/
	unsigned long 	IterLinSol1 = 0;
	CSysSolve femSystem2;
	IterLinSol1 = femSystem2.Solve(direct_solver->Jacobian, LinSysRes, LinSysSol, geometry, config);

//	ofstream myfile;
//	myfile.open ("Check_Adjoint.txt");
//
//	unsigned long iNode;
//	unsigned long realIndex;
//	su2double dxdE, Rv, PHI, dIdx;
//	myfile << "Node (" << "realIndex" << "," << "iVar" << "): " << "dxdE" << " " << "Rv" << " " << "PHI" << " " << "dIdx" << endl;
//
//	for (iNode = 0; iNode < nPoint; iNode++){
//			realIndex = geometry->node[iNode]->GetGlobalIndex();
//			for (iVar = 0; iVar < nVar; iVar++){
//				dxdE = LinSysSol_Direct.GetBlock(iNode, iVar);
//				Rv = LinSysRes.GetBlock(iNode, iVar);
//				PHI = LinSysSol.GetBlock(iNode, iVar);
//				dIdx = LinSysRes_ISens_X.GetBlock(iNode, iVar);
//				myfile << "Node (" << realIndex << "," << iVar << "): " << dxdE << " " << Rv << " " << PHI << " " << dIdx << endl;
//			}
//	}
//	myfile.close();

	su2double sensI_direct, sensI_adjoint;

	sensI_direct = dotProd(LinSysRes,LinSysSol_Direct);
	sensI_adjoint = dotProd(LinSysSol,LinSysRes_dSdv);

	su2double E_mod = config->GetElasticyMod();

	ofstream myfile_res;
	myfile_res.open ("Results_E.txt", ios::app);

    myfile_res << E_mod << " " << val_I << " "<< sensI_direct << " "  << sensI_adjoint << endl;

	myfile_res.close();

    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	for (iVar = 0; iVar < nVar; iVar++){

    		Solution[iVar] = LinSysSol.GetBlock(iPoint, iVar);

        	node[iPoint]->SetSolution(iVar, Solution[iVar]);

    	}

    }

	/*--- The the number of iterations of the linear solver ---*/

	SetIterLinSolver(IterLinSol);

}

void CFEM_ElasticitySolver_Adj::Run_Structural_Adjoint(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh,
		unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output){

//	unsigned short iMarker = 1;
//
//	cout << endl;
//	cout << "--------------- I'm going to run the structural adjoint --------------------" << endl;
//
//	Initialize_Structural_Adj(geometry, solver_container, numerics, config);
//
//	Compute_RefGeom_Sensitivity(geometry, solver_container, numerics, config);
//
//	BC_Clamped_Adj(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
//
//	Solve_System(geometry, solver_container, config);
//
//	BC_Clamped_Post_Adj(geometry, solver_container, numerics[FEA_TERM], config, iMarker);
//
//	Compute_DE_Sensitivity(geometry, solver_container, numerics, config);
//
//	Compute_RefGeom_Gradient(geometry, solver_container, numerics, config);

}


void CFEM_ElasticitySolver_Adj::Initialize_Structural_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

	cout << "STEP 1: Initialize the structural adjoint." << endl;

	unsigned long iPoint;
	unsigned short iVar;

	/*--- Initialize the gradient variable, and the residual and solution vectors ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++){
    	for (iVar = 0; iVar < nVar; iVar++){
        	node[iPoint]->SetGradient_Adj(iVar, 0.0);
    	}
		LinSysRes.SetBlock_Zero(iPoint);
		LinSysSol.SetBlock_Zero(iPoint);
    }


}

void CFEM_ElasticitySolver_Adj::Compute_RefGeom_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

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

    	LinSysRes.AddBlock(iPoint, Solution);

    }

}

void CFEM_ElasticitySolver_Adj::BC_Clamped_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker){

	cout << "STEP 3: Apply BC to the right hand side of the equations." << endl;

}

void CFEM_ElasticitySolver_Adj::Solve_System_Adj(CGeometry *geometry, CSolver **solver_container, CConfig *config){

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
	IterLinSol = femSystem.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

//	cout << "Adjoint system solved with " << IterLinSol << " iterations." << endl;

	unsigned long realIndex;

    for (iPoint=0; iPoint < nPointDomain; iPoint++){

    	realIndex = geometry->node[iPoint]->GetGlobalIndex();

    	for (iVar = 0; iVar < nVar; iVar++){

    		Solution[iVar] = LinSysSol.GetBlock(iPoint, iVar);

        	node[iPoint]->SetSolution_Adj(iVar, Solution[iVar]);

    	}

//    	cout << "Node " << realIndex << ". PHI = (" << Solution[0] << "," << Solution[1] << ")." << endl;

    }

	/*--- The the number of iterations of the linear solver ---*/

	SetIterLinSolver(IterLinSol);

}

void CFEM_ElasticitySolver_Adj::BC_Clamped_Post_Adj(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                             unsigned short val_marker){

	cout << "STEP 5: Reinforce BC to the sensitivities." << endl;

}


void CFEM_ElasticitySolver_Adj::Compute_DE_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

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


void CFEM_ElasticitySolver_Adj::Compute_RefGeom_Gradient(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config){

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





