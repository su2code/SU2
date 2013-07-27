/*!
 * \file solver_structure.cpp
 * \brief Main subrotuines for solving direct, adjoint and linearized problems.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.6
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

#include "../include/solver_structure.hpp"

CSolver::CSolver(void) {

  /*--- Array initialization ---*/
	Residual_RMS = NULL;
  Residual_Max = NULL;
	Residual = NULL;
	Residual_i = NULL;
	Residual_j = NULL;
  Point_Max = NULL;
	Solution = NULL;
	Solution_i = NULL;
	Solution_j = NULL;
	Vector = NULL;
	Vector_i = NULL;
	Vector_j = NULL;
	Res_Conv = NULL;
	Res_Visc = NULL;
	Res_Sour = NULL;
	Res_Conv_i = NULL;
	Res_Visc_i = NULL;
	Res_Conv_j = NULL;
	Res_Visc_j = NULL;
	Jacobian_i = NULL;
	Jacobian_j = NULL;
	Jacobian_ii = NULL;
	Jacobian_ij = NULL;
	Jacobian_ji = NULL;
	Jacobian_jj = NULL;
	Smatrix = NULL;
	cvector = NULL;
  node = NULL;

}

CSolver::~CSolver(void) {
/*  unsigned short iVar, iDim;
  unsigned long iPoint;
  
  if (Residual_RMS != NULL) delete [] Residual_RMS;
  if (Residual_Max != NULL) delete [] Residual_Max;
	if (Residual != NULL) delete [] Residual;
	if (Residual_i != NULL) delete [] Residual_i;
	if (Residual_j != NULL) delete [] Residual_j;
  if (Point_Max != NULL) delete [] Point_Max;
	if (Solution != NULL) delete [] Solution;
	if (Solution_i != NULL) delete [] Solution_i;
	if (Solution_j != NULL) delete [] Solution_j;
	if (Vector != NULL) delete [] Vector;
	if (Vector_i != NULL) delete [] Vector_i;
	if (Vector_j != NULL) delete [] Vector_j;
	if (Res_Conv != NULL) delete [] Res_Conv;
	if (Res_Visc != NULL) delete [] Res_Visc;
	if (Res_Sour != NULL) delete [] Res_Sour;
	if (Res_Conv_i != NULL) delete [] Res_Conv_i;
	if (Res_Visc_i != NULL) delete [] Res_Visc_i;
	if (Res_Visc_j != NULL) delete [] Res_Visc_j;
	if (Res_Sour_j != NULL) delete [] Res_Sour_j;
	if (rhs != NULL) delete [] rhs;
  
	if (Jacobian_i != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_i[iVar];
    delete [] Jacobian_i;
  }
  
	if (Jacobian_j != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_j[iVar];
    delete [] Jacobian_j;
  }
  
	if (Jacobian_MeanFlow_j != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_MeanFlow_j[iVar];
    delete [] Jacobian_MeanFlow_j;
  }
  
	if (Jacobian_ii != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_ii[iVar];
    delete [] Jacobian_ii;
  }
  
	if (Jacobian_ij != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_ij[iVar];
    delete [] Jacobian_ij;
  }
  
	if (Jacobian_ji != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_ji[iVar];
    delete [] Jacobian_ji;
  }
  
	if (Jacobian_jj != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_jj[iVar];
    delete [] Jacobian_jj;
  }
  
	if (Smatrix != NULL) {
    for (iDim = 0; iDim < nDim; iDim++)
      delete Smatrix[iDim];
    delete [] Smatrix;
  }
  
	if (cvector != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete cvector[iVar];
    delete [] cvector;
  }
  
  if (node != NULL) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      delete node[iPoint];
    }
    delete [] node;
  }
  
  //	delete [] **StiffMatrix_Elem;
  //	delete [] **StiffMatrix_Node;*/
  
}

void CSolver::SetResidual_RMS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;

#ifdef NO_MPI

	for (iVar = 0; iVar < nVar; iVar++)
    SetRes_RMS(iVar, max(EPS, sqrt(GetRes_RMS(iVar)/geometry->GetnPoint())));
  
#else
  
  int iProcessor;
  
	int nProcessor = MPI::COMM_WORLD.Get_size();
  double *sbuf_residual, *rbuf_residual;
  unsigned long *sbuf_point, *rbuf_point, Local_nPointDomain, Global_nPointDomain;
  
  /*--- Set the L2 Norm residual in all the processors ---*/
  sbuf_residual  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;
  
  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_RMS(iVar);
  Local_nPointDomain = geometry->GetnPointDomain();
  
  MPI::COMM_WORLD.Allreduce(sbuf_residual, rbuf_residual, nVar, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI::UNSIGNED_LONG, MPI::SUM);
  
  for (iVar = 0; iVar < nVar; iVar++)
    SetRes_RMS(iVar, max(EPS, sqrt(rbuf_residual[iVar]/Global_nPointDomain)));
  
  delete [] sbuf_residual;
  delete [] rbuf_residual;
  
  /*--- Set the Maximum residual in all the processors ---*/
  sbuf_residual = new double [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  sbuf_point = new unsigned long [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_point[iVar] = 0;
  
	rbuf_residual = new double [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_residual[iVar] = 0.0;
	rbuf_point = new unsigned long [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_point[iVar] = 0;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    sbuf_residual[iVar] = GetRes_Max(iVar);
    sbuf_point[iVar] = GetPoint_Max(iVar);
  }
  
	MPI::COMM_WORLD.Allgather(sbuf_residual, nVar, MPI::DOUBLE, rbuf_residual, nVar, MPI::DOUBLE);
  MPI::COMM_WORLD.Allgather(sbuf_point, nVar, MPI::UNSIGNED_LONG, rbuf_point, nVar, MPI::UNSIGNED_LONG);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      AddRes_Max(iVar, rbuf_residual[iProcessor*nVar+iVar], rbuf_point[iProcessor*nVar+iVar]);
    }
  }
  
  delete [] sbuf_residual;
  delete [] rbuf_residual;
  
  delete [] sbuf_point;
  delete [] rbuf_point;
  
#endif
  
}

void CSolver::SetGrid_Movement_Residual (CGeometry *geometry, CConfig *config) {

	unsigned short nDim = geometry->GetnDim();
	unsigned short nVar = GetnVar();
	double ProjGridVel, *Normal;

	//	Loop interior edges
	for(unsigned long iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	

		const unsigned long iPoint = geometry->edge[iEdge]->GetNode(0);
		const unsigned long jPoint = geometry->edge[iEdge]->GetNode(1);

		// Solution at each edge point
		double *Solution_i = node[iPoint]->GetSolution();
		double *Solution_j = node[jPoint]->GetSolution();

		for (unsigned short iVar = 0; iVar < nVar; iVar++)
			Solution[iVar] = 0.5* (Solution_i[iVar] + Solution_j[iVar]);

		// Grid Velocity at each edge point
		double *GridVel_i = geometry->node[iPoint]->GetGridVel();
		double *GridVel_j = geometry->node[jPoint]->GetGridVel();
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			Vector[iDim] = 0.5* (GridVel_i[iDim] + GridVel_j[iDim]);

		Normal = geometry->edge[iEdge]->GetNormal();
		//			dS = geometry->edge[iEdge]->GetArea_or_Length();

		ProjGridVel = 0.0;
		for (unsigned short iDim = 0; iDim < nDim; iDim++)
			ProjGridVel += Vector[iDim]*Normal[iDim];

		for (unsigned short iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = ProjGridVel*Solution[iVar];

		LinSysRes.SubtractBlock(iPoint, Residual);
		LinSysRes.AddBlock(jPoint, Residual);

	}

	//	Loop boundary edges
	for(unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for(unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			const unsigned long Point = geometry->vertex[iMarker][iVertex]->GetNode();

			// Solution at each edge point
			double *Solution = node[Point]->GetSolution();

			// Grid Velocity at each edge point
			double *GridVel = geometry->node[Point]->GetGridVel();

			// Summed normal components
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			//			dS = geometry->vertex[iMarker][iVertex]->GetArea_or_Length();

			ProjGridVel = 0.0;
			for (unsigned short iDim = 0; iDim < nDim; iDim++)
				ProjGridVel -= GridVel[iDim]*Normal[iDim];

			for (unsigned short iVar = 0; iVar < nVar; iVar++)
				Residual[iVar] = ProjGridVel*Solution[iVar];

			LinSysRes.AddBlock(Point, Residual);
		}
	}
}

void CSolver::SetAuxVar_Gradient_GG(CGeometry *geometry) {

	//	Internal variables
	unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
	unsigned short nDim = geometry->GetnDim(), iDim, iMarker;

	double AuxVar_Vertex, AuxVar_i, AuxVar_j, AuxVar_Average;
	double *Gradient, DualArea, Partial_Res, Grad_Val, *Normal;

	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		node[iPoint]->SetAuxVarGradientZero();		// Set Gradient to Zero

	//	Loop interior edges 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		AuxVar_i = node[iPoint]->GetAuxVar();
		AuxVar_j = node[jPoint]->GetAuxVar();

		Normal = geometry->edge[iEdge]->GetNormal();	
		AuxVar_Average =  0.5 * ( AuxVar_i + AuxVar_j);
		for(iDim = 0; iDim < nDim; iDim++) {
			Partial_Res = AuxVar_Average*Normal[iDim];
			node[iPoint]->AddAuxVarGradient(iDim, Partial_Res);
			node[jPoint]->SubtractAuxVarGradient(iDim, Partial_Res);					
		}				
	}

	//	Loop boundary edges
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			AuxVar_Vertex = node[Point]->GetAuxVar();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			for(iDim = 0; iDim < nDim; iDim++) {
				Partial_Res = AuxVar_Vertex*Normal[iDim];
				node[Point]->SubtractAuxVarGradient(iDim, Partial_Res);
			}
		}

	for (iPoint=0; iPoint<geometry->GetnPoint(); iPoint++)
		for(iDim = 0; iDim < nDim; iDim++) {
			Gradient = node[iPoint]->GetAuxVarGradient();
			DualArea = geometry->node[iPoint]->GetVolume();
			Grad_Val = Gradient[iDim]/(DualArea+EPS);
			node[iPoint]->SetAuxVarGradient(iDim,Grad_Val);				
		}	
}

void CSolver::SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config) {

	unsigned short iDim, jDim, iNeigh;
	unsigned short nDim = geometry->GetnDim();
	unsigned long iPoint, jPoint;
	double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j, weight;
	double *cvector;

	cvector = new double [nDim];

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

		Coord_i = geometry->node[iPoint]->GetCoord();
		AuxVar_i = node[iPoint]->GetAuxVar();

		/*--- Inizialization of variables ---*/
		for (iDim = 0; iDim < nDim; iDim++)
			cvector[iDim] = 0.0;
		double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();
			AuxVar_j = node[jPoint]->GetAuxVar();

			weight = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

			/*--- Sumations for entries of upper triangular matrix R ---*/
			r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
			r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
			r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
			if (nDim == 3) {
				r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
			}

			/*--- Entries of c:= transpose(A)*b ---*/
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/(weight);
		}

		/*--- Entries of upper triangular matrix R ---*/
		r11 = sqrt(r11);
		r12 = r12/(r11);
		r22 = sqrt(r22-r12*r12);
		if (nDim == 3) {
			r13 = r13/(r11);
			r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
			r33 = sqrt(r33-r23*r23-r13*r13);
		}

		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		if (nDim == 2) {
			double detR2 = (r11*r22)*(r11*r22);
			Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
			Smatrix[0][1] = -r11*r12/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = r11*r11/(detR2);
		}
		else {
			double detR2 = (r11*r22*r33)*(r11*r22*r33);
			double z11, z12, z13, z22, z23, z33;
			z11 = r22*r33;
			z12 = -r12*r33;
			z13 = r12*r23-r13*r22;
			z22 = r11*r33;
			z23 = -r11*r23;
			z33 = r11*r22;
			Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
			Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
			Smatrix[0][2] = (z13*z33)/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
			Smatrix[1][2] = (z23*z33)/(detR2);
			Smatrix[2][0] = Smatrix[0][2];
			Smatrix[2][1] = Smatrix[1][2];
			Smatrix[2][2] = (z33*z33)/(detR2);
		}

		/*--- Computation of the gradient: S*c ---*/
		double product;
		for (iDim = 0; iDim < nDim; iDim++) {
			product = 0.0;
			for (jDim = 0; jDim < nDim; jDim++)
				product += Smatrix[iDim][jDim]*cvector[jDim];
			if (geometry->node[iPoint]->GetDomain())
				node[iPoint]->SetAuxVarGradient(iDim,product);
		}
	}

	delete [] cvector;
}

void CSolver::SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config) {
	unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
	unsigned short iVar, iDim, iMarker;
	double *Solution_Vertex, *Solution_i, *Solution_j, Solution_Average, **Gradient, DualArea, 
	Partial_Res, Grad_Val, *Normal;
	
	/*--- Set Gradient to Zero ---*/
	for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		node[iPoint]->SetGradientZero();		

	/*--- Loop interior edges ---*/
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);

		Solution_i = node[iPoint]->GetSolution();
		Solution_j = node[jPoint]->GetSolution();
		Normal = geometry->edge[iEdge]->GetNormal();
		for(iVar = 0; iVar< nVar; iVar++) {
			Solution_Average =  0.5 * (Solution_i[iVar] + Solution_j[iVar]);
			for(iDim = 0; iDim < nDim; iDim++) {
				Partial_Res = Solution_Average*Normal[iDim];
				if (geometry->node[iPoint]->GetDomain())
					node[iPoint]->AddGradient(iVar, iDim, Partial_Res);
				if (geometry->node[jPoint]->GetDomain())
					node[jPoint]->SubtractGradient(iVar, iDim, Partial_Res);					
			}				
		}
	}

	/*--- Loop boundary edges ---*/
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Solution_Vertex = node[Point]->GetSolution();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			for(iVar = 0; iVar < nVar; iVar++)
				for(iDim = 0; iDim < nDim; iDim++) {
					Partial_Res = Solution_Vertex[iVar]*Normal[iDim];
					if (geometry->node[Point]->GetDomain())						
						node[Point]->SubtractGradient(iVar,iDim, Partial_Res);
				}
		}
	}

	/*--- Compute gradient ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
		for(iVar = 0; iVar < nVar; iVar++)
			for(iDim = 0; iDim < nDim; iDim++) {
				Gradient = node[iPoint]->GetGradient();
				DualArea = geometry->node[iPoint]->GetVolume();
				Grad_Val = Gradient[iVar][iDim] / (DualArea+EPS);
				node[iPoint]->SetGradient(iVar,iDim,Grad_Val);
			}
  
  /*--- Gradient MPI ---*/
  Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config) {
	unsigned short iDim, jDim, iVar, iNeigh;
	unsigned long iPoint, jPoint;
	double *Coord_i, *Coord_j, *Solution_i, *Solution_j,
	r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight, detR2, z11, z12, z13, 
	z22, z23, z33, product;
	double **cvector;

	cvector = new double* [nVar];
	for (iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];

	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

		Coord_i = geometry->node[iPoint]->GetCoord();
		Solution_i = node[iPoint]->GetSolution();

		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nVar; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] = 0.0;
		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();
			Solution_j = node[jPoint]->GetSolution();

			weight = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

			/*--- Sumations for entries of upper triangular matrix R ---*/
			r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
			r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
			r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
			if (nDim == 3) {
				r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
			}

			/*--- Entries of c:= transpose(A)*b ---*/
			for (iVar = 0; iVar < nVar; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/(weight);
		}

		/*--- Entries of upper triangular matrix R ---*/
		r11 = sqrt(r11);
		r12 = r12/(r11);
		r22 = sqrt(r22-r12*r12);
		if (nDim == 3) {
			r13 = r13/(r11);
			r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
			r33 = sqrt(r33-r23*r23-r13*r13);
		}
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		if (nDim == 2) {
			detR2 = (r11*r22)*(r11*r22);
			Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
			Smatrix[0][1] = -r11*r12/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = r11*r11/(detR2);
		}
		else {
			detR2 = (r11*r22*r33)*(r11*r22*r33);
			z11 = r22*r33;
			z12 = -r12*r33;
			z13 = r12*r23-r13*r22;
			z22 = r11*r33;
			z23 = -r11*r23;
			z33 = r11*r22;
			Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
			Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
			Smatrix[0][2] = (z13*z33)/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
			Smatrix[1][2] = (z23*z33)/(detR2);
			Smatrix[2][0] = Smatrix[0][2];
			Smatrix[2][1] = Smatrix[1][2];
			Smatrix[2][2] = (z33*z33)/(detR2);
		}
		/*--- Computation of the gradient: S*c ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				product = 0.0;
				for (jDim = 0; jDim < nDim; jDim++)
					product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
				node[iPoint]->SetGradient(iVar,iDim,product);
			}
		}
	}

	/*--- Deallocate memory ---*/
	for (iVar = 0; iVar < nVar; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
  
  /*--- Gradient MPI ---*/
  Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetRotVel_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iDim, jDim, iVar, iNeigh;
	unsigned long iPoint, jPoint;
	double *Coord_i, *Coord_j, *Solution_i, *Solution_j, Smatrix[3][3],
	r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight, detR2, z11, z12, z13,
	z22, z23, z33, product;
	double **cvector;
  
  /*--- Note that all nVar entries in this routine have been changed to nDim ---*/
	cvector = new double* [nDim];
	for (iVar = 0; iVar < nDim; iVar++)
		cvector[iVar] = new double [nDim];
  
	/*--- Loop over points of the grid ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
		Coord_i = geometry->node[iPoint]->GetCoord();
		Solution_i = geometry->node[iPoint]->GetRotVel();
    
		/*--- Inizialization of variables ---*/
		for (iVar = 0; iVar < nDim; iVar++)
			for (iDim = 0; iDim < nDim; iDim++)
				cvector[iVar][iDim] = 0.0;
		r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
		for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
			jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
			Coord_j = geometry->node[jPoint]->GetCoord();
			Solution_j = geometry->node[jPoint]->GetRotVel();
      
			weight = 0.0;
			for (iDim = 0; iDim < nDim; iDim++)
				weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
			/*--- Sumations for entries of upper triangular matrix R ---*/
			r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
			r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
			r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
			if (nDim == 3) {
				r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
				r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
				r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
			}
      
			/*--- Entries of c:= transpose(A)*b ---*/
			for (iVar = 0; iVar < nDim; iVar++)
				for (iDim = 0; iDim < nDim; iDim++)
					cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/(weight);
		}
    
		/*--- Entries of upper triangular matrix R ---*/
		r11 = sqrt(r11);
		r12 = r12/(r11);
		r22 = sqrt(r22-r12*r12);
		if (nDim == 3) {
			r13 = r13/(r11);
			r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
			r33 = sqrt(r33-r23*r23-r13*r13);
		}
		/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
		if (nDim == 2) {
			detR2 = (r11*r22)*(r11*r22);
			Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
			Smatrix[0][1] = -r11*r12/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = r11*r11/(detR2);
		}
		else {
			detR2 = (r11*r22*r33)*(r11*r22*r33);
			z11 = r22*r33;
			z12 = -r12*r33;
			z13 = r12*r23-r13*r22;
			z22 = r11*r33;
			z23 = -r11*r23;
			z33 = r11*r22;
			Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
			Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
			Smatrix[0][2] = (z13*z33)/(detR2);
			Smatrix[1][0] = Smatrix[0][1];
			Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
			Smatrix[1][2] = (z23*z33)/(detR2);
			Smatrix[2][0] = Smatrix[0][2];
			Smatrix[2][1] = Smatrix[1][2];
			Smatrix[2][2] = (z33*z33)/(detR2);
		}
		/*--- Computation of the gradient: S*c ---*/
		for (iVar = 0; iVar < nDim; iVar++) {
			for (iDim = 0; iDim < nDim; iDim++) {
				product = 0.0;
				for (jDim = 0; jDim < nDim; jDim++)
					product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
				geometry->node[iPoint]->SetRotVel_Grad(iVar,iDim,product);
			}
		}
	}
  
	/*--- Deallocate memory ---*/
	for (iVar = 0; iVar < nDim; iVar++)
		delete [] cvector[iVar];
	delete [] cvector;
  
  /*--- Gradient MPI ---*/
  // TO DO!!!
  //Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetSurface_Gradient(CGeometry *geometry, CConfig *config) {
	unsigned short iDim, jDim, iVar, iNeigh, iMarker, Boundary;
	unsigned short nDim = geometry->GetnDim();
	unsigned long iPoint, jPoint, iVertex;
	double *Coord_i, *Coord_j, *Solution_i, *Solution_j;
	double **Smatrix, **cvector;

	cvector = new double* [nVar];
	Smatrix = new double* [nDim];
	for (iVar = 0; iVar < nVar; iVar++)
		cvector[iVar] = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];

	/*--- Loop over boundary markers to select those for Euler or NS walls ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		switch (Boundary) {
		case EULER_WALL: case HEAT_FLUX: case ISOTHERMAL:

			/*--- Loop over points on the surface (Least-Squares approximation) ---*/
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Coord_i = geometry->node[iPoint]->GetCoord();
					Solution_i = node[iPoint]->GetSolution();

					/*--- Inizialization of variables ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						for (iDim = 0; iDim < nDim; iDim++)
							cvector[iVar][iDim] = 0.0;
					double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

					for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
						jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
						Coord_j = geometry->node[jPoint]->GetCoord();
						Solution_j = node[jPoint]->GetSolution();

						double weight = 0.0;
						for (iDim = 0; iDim < nDim; iDim++)
							weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

						/*--- Sumations for entries of upper triangular matrix R ---*/
						r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
						r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
						r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
						if (nDim == 3) {
							r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
							r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
							r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
							r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
						}

						/*--- Entries of c:= transpose(A)*b ---*/
						for (iVar = 0; iVar < nVar; iVar++)
							for (iDim = 0; iDim < nDim; iDim++)
								cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]
								                                                                 -Solution_i[iVar])/weight;
					}

					/*--- Entries of upper triangular matrix R ---*/
					r11 = sqrt(r11);
					r12 = r12/(r11+EPS);
					r22 = sqrt(r22-r12*r12);
					if (nDim == 3) {
						r13 = r13/(r11+EPS);
						r23 = r23_a/(r22+EPS) - r23_b*r12/(r11*r22+EPS);
						r33 = sqrt(r33-r23*r23-r13*r13);
					}
					/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
					if (nDim == 2) {
						double detR2 = (r11*r22)*(r11*r22);
						Smatrix[0][0] = (r12*r12+r22*r22)/(detR2+EPS);
						Smatrix[0][1] = -r11*r12/(detR2+EPS);
						Smatrix[1][0] = Smatrix[0][1];
						Smatrix[1][1] = r11*r11/(detR2+EPS);
					}
					else {
						double detR2 = (r11*r22*r33)*(r11*r22*r33);
						double z11, z12, z13, z22, z23, z33; // aux vars
						z11 = r22*r33;
						z12 = -r12*r33;
						z13 = r12*r23-r13*r22;
						z22 = r11*r33;
						z23 = -r11*r23;
						z33 = r11*r22;
						Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2+EPS);
						Smatrix[0][1] = (z12*z22+z13*z23)/(detR2+EPS);
						Smatrix[0][2] = (z13*z33)/(detR2+EPS);
						Smatrix[1][0] = Smatrix[0][1];
						Smatrix[1][1] = (z22*z22+z23*z23)/(detR2+EPS);
						Smatrix[1][2] = (z23*z33)/(detR2+EPS);
						Smatrix[2][0] = Smatrix[0][2];
						Smatrix[2][1] = Smatrix[1][2];
						Smatrix[2][2] = (z33*z33)/(detR2+EPS);
					}
					/*--- Computation of the gradient: S*c ---*/
					double product;
					for (iVar = 0; iVar < nVar; iVar++) {
						for (iDim = 0; iDim < nDim; iDim++) {
							product = 0.0;
							for (jDim = 0; jDim < nDim; jDim++)
								product += Smatrix[iDim][jDim]*cvector[iVar][jDim];
							node[iPoint]->SetGradient(iVar,iDim,product);
						}
					}
				}

			} /*--- End of loop over surface points ---*/
			break;
		default:
			break;
		}	
	}

	/*--- Memory deallocation ---*/	
	for (iVar = 0; iVar < nVar; iVar++)
		delete cvector[iVar];
	for (iDim = 0; iDim < nDim; iDim++)
		delete Smatrix[iDim];
	delete [] cvector;
	delete [] Smatrix;
}

void CSolver::SetAuxVar_Surface_Gradient(CGeometry *geometry, CConfig *config) {

	unsigned short iDim, jDim, iNeigh, iMarker, Boundary;
	unsigned short nDim = geometry->GetnDim();
	unsigned long iPoint, jPoint, iVertex;
	double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j;
	double **Smatrix, *cvector;

	Smatrix = new double* [nDim];
	cvector = new double [nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Smatrix[iDim] = new double [nDim];


	/*--- Loop over boundary markers to select those for Euler or NS walls ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		Boundary = config->GetMarker_All_Boundary(iMarker);
		switch (Boundary) {
      case EULER_WALL: case HEAT_FLUX: case ISOTHERMAL:
        
			/*--- Loop over points on the surface (Least-Squares approximation) ---*/
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				if (geometry->node[iPoint]->GetDomain()) {
					Coord_i = geometry->node[iPoint]->GetCoord();
					AuxVar_i = node[iPoint]->GetAuxVar();

					/*--- Inizialization of variables ---*/
					for (iDim = 0; iDim < nDim; iDim++)
						cvector[iDim] = 0.0;
					double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;

					for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
						jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
						Coord_j = geometry->node[jPoint]->GetCoord();
						AuxVar_j = node[jPoint]->GetAuxVar();

						double weight = 0;
						for (iDim = 0; iDim < nDim; iDim++)
							weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

						/*--- Sumations for entries of upper triangular matrix R ---*/
						r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
						r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
						r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
						if (nDim == 3) {
							r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
							r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
							r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
							r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
						}

						/*--- Entries of c:= transpose(A)*b ---*/
						for (iDim = 0; iDim < nDim; iDim++)
							cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/weight;
					}

					/*--- Entries of upper triangular matrix R ---*/
					r11 = sqrt(r11);
					r12 = r12/r11;
					r22 = sqrt(r22-r12*r12);
					if (nDim == 3) {
						r13 = r13/r11;
						r23 = r23_a/r22 - r23_b*r12/(r11*r22);
						r33 = sqrt(r33-r23*r23-r13*r13);
					}
					/*--- S matrix := inv(R)*traspose(inv(R)) ---*/
					if (nDim == 2) {
						double detR2 = (r11*r22)*(r11*r22);
						Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
						Smatrix[0][1] = -r11*r12/detR2;
						Smatrix[1][0] = Smatrix[0][1];
						Smatrix[1][1] = r11*r11/detR2;
					}
					else {
						double detR2 = (r11*r22*r33)*(r11*r22*r33);
						double z11, z12, z13, z22, z23, z33; // aux vars
						z11 = r22*r33;
						z12 = -r12*r33;
						z13 = r12*r23-r13*r22;
						z22 = r11*r33;
						z23 = -r11*r23;
						z33 = r11*r22;
						Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
						Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
						Smatrix[0][2] = (z13*z33)/detR2;
						Smatrix[1][0] = Smatrix[0][1];
						Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
						Smatrix[1][2] = (z23*z33)/detR2;
						Smatrix[2][0] = Smatrix[0][2];
						Smatrix[2][1] = Smatrix[1][2];
						Smatrix[2][2] = (z33*z33)/detR2;
					}
					/*--- Computation of the gradient: S*c ---*/
					double product;
					for (iDim = 0; iDim < nDim; iDim++) {
						product = 0.0;
						for (jDim = 0; jDim < nDim; jDim++)
							product += Smatrix[iDim][jDim]*cvector[jDim];
						node[iPoint]->SetAuxVarGradient(iDim, product);
					}
				}
			} /*--- End of loop over surface points ---*/
			break;
		default:
			break;
		}	
	}

	/*--- Memory deallocation ---*/	
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] cvector;
	delete [] Smatrix;
}

void CSolver::SetSolution_Limiter(CGeometry *geometry, CConfig *config) {
		
	unsigned long iEdge, iPoint, jPoint;
	unsigned short iVar, iDim;
	double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j, *Solution_i, *Solution_j, 
	dave, LimK, eps2, dm, dp, du, limiter;

	/*--- Initialize solution max and solution min in the entire domain --*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->SetSolution_Max(iVar, -EPS);
			node[iPoint]->SetSolution_Min(iVar, EPS);
		}
	}

	/*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
		
		/*--- Point identification, Normal vector and area ---*/
		iPoint = geometry->edge[iEdge]->GetNode(0); 
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		/*--- Get the conserved variables ---*/
		Solution_i = node[iPoint]->GetSolution();
		Solution_j = node[jPoint]->GetSolution();

		/*--- Compute the maximum, and minimum values for nodes i & j ---*/
		for (iVar = 0; iVar < nVar; iVar++) {
			du = (Solution_j[iVar] - Solution_i[iVar]);
			node[iPoint]->SetSolution_Min(iVar, min(node[iPoint]->GetSolution_Min(iVar), du));
			node[iPoint]->SetSolution_Max(iVar, max(node[iPoint]->GetSolution_Max(iVar), du));
			node[jPoint]->SetSolution_Min(iVar, min(node[jPoint]->GetSolution_Min(iVar), -du));
			node[jPoint]->SetSolution_Max(iVar, max(node[jPoint]->GetSolution_Max(iVar), -du));
		}
	}

	/*--- Initialize the limiter --*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		for (iVar = 0; iVar < nVar; iVar++) {
			node[iPoint]->SetLimiter(iVar, 2.0);
		}
	}
	
  switch (config->GetKind_SlopeLimit()) {
    
    /*--- Minmod (Roe 1984) limiter ---*/
    case MINMOD:

      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        
        iPoint     = geometry->edge[iEdge]->GetNode(0);
        jPoint     = geometry->edge[iEdge]->GetNode(1);
        Gradient_i = node[iPoint]->GetGradient();
        Gradient_j = node[jPoint]->GetGradient();
        Coord_i    = geometry->node[iPoint]->GetCoord();
        Coord_j    = geometry->node[jPoint]->GetCoord();
        
        for (iVar = 0; iVar < nVar; iVar++) {          
          
          /*--- Calculate the interface left gradient, delta- (dm) ---*/
          dm = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
          
          /*--- Calculate the interface right gradient, delta+ (dp) ---*/
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          
          limiter = max(0.0, min(1.0,dp/dm));
          
          if (limiter < node[iPoint]->GetLimiter(iVar))
            if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter(iVar, limiter);
          
          /*-- Repeat for point j on the edge ---*/
          dm = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
          
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          
          limiter = max(0.0, min(1.0,dp/dm));
          
          if (limiter < node[jPoint]->GetLimiter(iVar))
            if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter(iVar, limiter);   
        }
      }
      break;
    
    /*--- Venkatakrishnan (Venkatakrishnan 1994) limiter ---*/
    case VENKATAKRISHNAN:
      
      for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
        
        iPoint     = geometry->edge[iEdge]->GetNode(0);
        jPoint     = geometry->edge[iEdge]->GetNode(1);
        Solution_i = node[iPoint]->GetSolution();
        Solution_j = node[jPoint]->GetSolution();
        Gradient_i = node[iPoint]->GetGradient();
        Gradient_j = node[jPoint]->GetGradient();
        Coord_i    = geometry->node[iPoint]->GetCoord();
        Coord_j    = geometry->node[jPoint]->GetCoord();
				
        for (iVar = 0; iVar < nVar; iVar++) {
          
          /*-- Get limiter parameters from the configuration file ---*/
          dave = config->GetRefElemLength();
          LimK = config->GetLimiterCoeff();
          eps2 = pow((LimK*dave), 3.0);
          
          /*--- Calculate the interface left gradient, delta- (dm) ---*/
          dm = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
          
          /*--- Calculate the interface right gradient, delta+ (dp) ---*/
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          
          limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
          
          if (limiter < node[iPoint]->GetLimiter(iVar))
            if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter(iVar, limiter);
          
          /*-- Repeat for point j on the edge ---*/
          dave = config->GetRefElemLength();
          LimK = config->GetLimiterCoeff();
          eps2 = pow((LimK*dave), 3.0);
					
          dm = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
          
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          
          limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
          
          if (limiter < node[jPoint]->GetLimiter(iVar))
            if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter(iVar, limiter);
        }
      }
      break;
  }
  
  /*--- Limiter MPI ---*/
  Set_MPI_Solution_Limiter(geometry, config);
	
}

void CSolver::SetPressureLaplacian(CGeometry *geometry, double *PressureLaplacian) {

	unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
	unsigned short iMarker, iVar;
	double DualArea, Partial_Res, *Normal, Area;
	double **UxVar_Gradient, **UyVar_Gradient;

	UxVar_Gradient = new double* [geometry->GetnPoint()];
	UyVar_Gradient = new double* [geometry->GetnPoint()];
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		UxVar_Gradient[iPoint] = new double [2];
		UyVar_Gradient[iPoint] = new double [2];
	}

	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		for(iVar = 0; iVar < 2; iVar++) {
			UxVar_Gradient[iPoint][iVar] = 0.0;
			UyVar_Gradient[iPoint][iVar] = 0.0;
		}

	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
		
		Partial_Res =  0.5 * ( node[iPoint]->GetSolution(1) + node[jPoint]->GetSolution(1)) * Normal[0];
		UxVar_Gradient[iPoint][0] += Partial_Res;
		UxVar_Gradient[jPoint][0] -= Partial_Res;
		
		Partial_Res =  0.5 * ( node[iPoint]->GetSolution(1) + node[jPoint]->GetSolution(1)) * Normal[1];
		UxVar_Gradient[iPoint][1] += Partial_Res;
		UxVar_Gradient[jPoint][1] -= Partial_Res;
		
		Partial_Res =  0.5 * ( node[iPoint]->GetSolution(2) + node[jPoint]->GetSolution(2)) * Normal[0];
		UyVar_Gradient[iPoint][0] += Partial_Res;
		UyVar_Gradient[jPoint][0] -= Partial_Res;
		
		Partial_Res =  0.5 * ( node[iPoint]->GetSolution(2) + node[jPoint]->GetSolution(2)) * Normal[1];
		UyVar_Gradient[iPoint][1] += Partial_Res;
		UyVar_Gradient[jPoint][1] -= Partial_Res;
		
	}
	
	/*---	Loop boundary edges ---*/
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
			
			Partial_Res =  node[Point]->GetSolution(1) * Normal[0];
			UxVar_Gradient[Point][0] -= Partial_Res;
			
			Partial_Res =  node[Point]->GetSolution(1) * Normal[1];
			UxVar_Gradient[Point][1] -= Partial_Res;
			
			Partial_Res =  node[Point]->GetSolution(2) * Normal[0];
			UyVar_Gradient[Point][0] -= Partial_Res;
			
			Partial_Res =  node[Point]->GetSolution(2) * Normal[1];
			UyVar_Gradient[Point][1] -= Partial_Res;
		}

	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		DualArea = geometry->node[iPoint]->GetVolume();
		PressureLaplacian[iPoint] = (UxVar_Gradient[iPoint][0]*UxVar_Gradient[iPoint][0] + UyVar_Gradient[iPoint][1]*UyVar_Gradient[iPoint][1] +
		UxVar_Gradient[iPoint][1]*UyVar_Gradient[iPoint][0] + UxVar_Gradient[iPoint][0]*UyVar_Gradient[iPoint][1])/DualArea ;
	}
	

	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		delete[] UxVar_Gradient[iPoint];
		delete[] UyVar_Gradient[iPoint];
	}

	delete[] UxVar_Gradient;
	delete[] UyVar_Gradient;
	
}

void CSolver::Gauss_Elimination(double** A, double* rhs, unsigned long nVar) {
	unsigned long jVar, kVar, iVar;
	double weight, aux;
	
	if (nVar == 1)
		rhs[0] /= (A[0][0]+EPS*EPS);
	else {
		/*--- Transform system in Upper Matrix ---*/
		for (iVar = 1; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < iVar; jVar++) {
				weight = A[iVar][jVar]/(A[jVar][jVar]+EPS*EPS);
				for (kVar = jVar; kVar < nVar; kVar++)
					A[iVar][kVar] -= weight*A[jVar][kVar];
				rhs[iVar] -= weight*rhs[jVar];
			}
		}
		/*--- Backwards substitution ---*/
		rhs[nVar-1] = rhs[nVar-1]/(A[nVar-1][nVar-1]+EPS*EPS);
		for (short iVar = nVar-2; iVar >= 0; iVar--) {
			aux = 0;
			for (jVar = iVar+1; jVar < nVar; jVar++)
				aux += A[iVar][jVar]*rhs[jVar];
			rhs[iVar] = (rhs[iVar]-aux)/(A[iVar][iVar]+EPS*EPS);
			if (iVar == 0) break;
		}
	}
}

CBaselineSolver::CBaselineSolver(void) : CSolver() { }

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	unsigned long iPoint, index, flowIter, iPoint_Global;
  long iPoint_Local;
  char buffer[50];
  unsigned short iField, val_iZone, iVar;
  string Tag, text_line, AdjExt, UnstExt;
  
  unsigned short Kind_ObjFunc = config->GetKind_ObjFunc();
  unsigned short nZone = geometry->GetnZone();
  
	/*--- Define geometry constants in the solver structure ---*/
	nDim = geometry->GetnDim();
  
	/*--- Allocate the node variables ---*/
	node = new CVariable*[geometry->GetnPoint()];
  
  /*--- Restart the solution from file information ---*/
  ifstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
	if (config->GetAdjoint())
		filename = config->GetSolution_AdjFileName();
	else
		filename = config->GetSolution_FlowFileName();
  
	/*--- Remove restart filename extension ---*/
	filename.erase(filename.end()-4, filename.end());
  
	/*--- The adjoint problem requires a particular extension. ---*/
	if (config->GetAdjoint()) {
		switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT:      AdjExt = "_cd";   break;
      case LIFT_COEFFICIENT:      AdjExt = "_cl";   break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf";  break;
      case PRESSURE_COEFFICIENT:  AdjExt = "_cp";   break;
      case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx";  break;
      case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy";  break;
      case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz";  break;
      case EFFICIENCY:            AdjExt = "_eff";  break;
      case EQUIVALENT_AREA:       AdjExt = "_ea";   break;
      case NEARFIELD_PRESSURE:    AdjExt = "_nfp";  break;
      case FORCE_X_COEFFICIENT:   AdjExt = "_cfx";  break;
      case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy";  break;
      case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz";  break;
      case THRUST_COEFFICIENT:    AdjExt = "_ct";   break;
      case TORQUE_COEFFICIENT:    AdjExt = "_cq";   break;
      case FIGURE_OF_MERIT:       AdjExt = "_merit";break;
      case FREE_SURFACE:           AdjExt = "_fs";   break;
      case NOISE:                 AdjExt = "_fwh";  break;
      case HEAT_LOAD:                 AdjExt = "_Q";  break;
		}
		filename.append(AdjExt);
	}
  
	/*--- Multi-zone restart files. ---*/
	if (nZone > 1 && !(config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
    val_iZone = config->GetExtIter();
		filename.erase(filename.end()-4, filename.end());
		sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	}
  
	/*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		flowIter = val_iZone;
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		flowIter   = config->GetExtIter();
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	} else
    filename.append(".dat");
  
  restart_file.open(filename.data(), ios::in);

  /*--- In case there is no restart file ---*/
  if (restart_file.fail()) {
    cout << "There is no SU2 restart file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(1);
  }
  
  /*--- Output the file name to the console. ---*/
  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename << "." << endl;
  
  /*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
  long *Global2Local = new long[geometry->GetGlobal_nPointDomain()];
  
  /*--- First, set all indices to a negative value by default ---*/
  for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++)
    Global2Local[iPoint] = -1;
  
  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  
  
  /*--- Identify the number of fields (and names) in the restart file ---*/
  getline (restart_file, text_line);
  stringstream ss(text_line);
  while (ss >> Tag) {
    config->fields.push_back(Tag);
    if (ss.peek() == ',') ss.ignore();
  }
  
  /*--- Set the number of variables, one per field in the
   restart file (without including the PointID) ---*/
  nVar = config->fields.size() - 1;
  double Solution[nVar];
  
  /*--- Read all lines in the restart file ---*/
  iPoint_Global = 0;
  while (getline (restart_file, text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1.
     Otherwise, the local index for this node on the current processor
     will be returned and used to instantiate the vars. ---*/
    iPoint_Local = Global2Local[iPoint_Global];
    if (iPoint_Local >= 0) {
      
      /*--- The PointID is not stored --*/
      point_line >> index;
      
      /*--- Store the solution --*/
      for (iField = 0; iField < nVar; iField++)
        point_line >> Solution[iField];
      
      node[iPoint_Local] = new CBaselineVariable(Solution, nVar, config);
    }
    iPoint_Global++;
  }
  
  /*--- Instantiate the variable class with an arbitrary solution
   at any halo/periodic nodes. The initial solution can be arbitrary,
   because a send/recv is performed immediately in the solver. ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Solution[iVar] = 0.0;

  for(iPoint = geometry->GetnPointDomain(); iPoint < geometry->GetnPoint(); iPoint++)
    node[iPoint] = new CBaselineVariable(Solution, nVar, config);
  
  /*--- Close the restart file ---*/
  restart_file.close();
  
  /*--- Free memory needed for the transformation ---*/
  delete [] Global2Local;
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

void CBaselineSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
	unsigned short iVar, iMarker, iPeriodic_Index;
	unsigned long iVertex, iPoint, nVertex, nBuffer_Vector;
	double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi,
  psi, cosPsi, sinPsi, *newSolution = NULL, *Buffer_Receive_U = NULL;
	short SendRecv;
	int send_to, receive_from;
  
#ifndef NO_MPI
  
  MPI::COMM_WORLD.Barrier();
	double *Buffer_Send_U = NULL;
  
#endif
  
	newSolution = new double[nVar];

	/*--- Send-Receive boundary conditions ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == SEND_RECEIVE) {
			SendRecv = config->GetMarker_All_SendRecv(iMarker);
			nVertex = geometry->nVertex[iMarker];
			nBuffer_Vector = nVertex*nVar;
			send_to = SendRecv-1;
			receive_from = abs(SendRecv)-1;
      
#ifndef NO_MPI
      
			/*--- Send information using MPI  ---*/
			if (SendRecv > 0) {
        Buffer_Send_U = new double[nBuffer_Vector];
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Send_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
				}
        MPI::COMM_WORLD.Bsend(Buffer_Send_U, nBuffer_Vector, MPI::DOUBLE, send_to, 0); delete [] Buffer_Send_U;
			}
      
#endif

			/*--- Receive information  ---*/
			if (SendRecv < 0) {
        Buffer_Receive_U = new double [nBuffer_Vector];
        
#ifdef NO_MPI
        
				/*--- Receive information without MPI ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          for (iVar = 0; iVar < nVar; iVar++)
            Buffer_Receive_U[iVar*nVertex+iVertex] = node[iPoint]->GetSolution(iVar);
        }
        
#else
        
        MPI::COMM_WORLD.Recv(Buffer_Receive_U, nBuffer_Vector, MPI::DOUBLE, receive_from, 0);
        
#endif
        
				/*--- Do the coordinate transformation ---*/
				for (iVertex = 0; iVertex < nVertex; iVertex++) {
          
					/*--- Find point and its type of transformation ---*/
					iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
					iPeriodic_Index = geometry->vertex[iMarker][iVertex]->GetRotation_Type();
          
					/*--- Retrieve the supplied periodic information. ---*/
					angles = config->GetPeriodicRotation(iPeriodic_Index);
          
					/*--- Store angles separately for clarity. ---*/
					theta    = angles[0];   phi    = angles[1]; psi    = angles[2];
					cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
					sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
          
					/*--- Compute the rotation matrix. Note that the implicit
					 ordering is rotation about the x-axis, y-axis,
					 then z-axis. Note that this is the transpose of the matrix
					 used during the preprocessing stage. ---*/
					rotMatrix[0][0] = cosPhi*cosPsi; rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi; rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
					rotMatrix[0][1] = cosPhi*sinPsi; rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi; rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
					rotMatrix[0][2] = -sinPhi; rotMatrix[1][2] = sinTheta*cosPhi; rotMatrix[2][2] = cosTheta*cosPhi;
          
					/*--- Copy conserved variables before performing transformation. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						newSolution[iVar] = Buffer_Receive_U[iVar*nVertex+iVertex];
          
					/*--- Rotate the momentum components. ---*/
					if (nDim == 2) {
						newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[2*nVertex+iVertex];
						newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[2*nVertex+iVertex];
					}
					else {
						newSolution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[0][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[0][2]*Buffer_Receive_U[3*nVertex+iVertex];
						newSolution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[1][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[1][2]*Buffer_Receive_U[3*nVertex+iVertex];
						newSolution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertex+iVertex] + rotMatrix[2][1]*Buffer_Receive_U[2*nVertex+iVertex] + rotMatrix[2][2]*Buffer_Receive_U[3*nVertex+iVertex];
					}
          
					/*--- Copy transformed conserved variables back into buffer. ---*/
					for (iVar = 0; iVar < nVar; iVar++)
						Buffer_Receive_U[iVar*nVertex+iVertex] = newSolution[iVar];
          
          for (iVar = 0; iVar < nVar; iVar++)
            node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertex+iVertex]);
          
				}
        delete [] Buffer_Receive_U;
			}
		}
	}
	delete [] newSolution;
  
#ifndef NO_MPI
  
  MPI::COMM_WORLD.Barrier();
  
#endif
  
}

void CBaselineSolver::GetRestart(CGeometry *geometry, CConfig *config, unsigned short val_iZone) {
  
	int rank = MASTER_NODE;
#ifndef NO_MPI
	rank = MPI::COMM_WORLD.Get_rank();
#endif
  
	/*--- Restart the solution from file information ---*/
	string filename;
	unsigned long iPoint, index, flowIter;
	char buffer[50];
	string UnstExt, text_line, AdjExt;
	ifstream restart_file;
	unsigned short nZone = geometry->GetnZone();
  unsigned short Kind_ObjFunc = config->GetKind_ObjFunc();
  unsigned short iField;
  
  /*--- Retrieve filename from config ---*/
	if (config->GetAdjoint())
		filename = config->GetSolution_AdjFileName();
	else
		filename = config->GetSolution_FlowFileName();
  
	/*--- Remove restart filename extension ---*/
	filename.erase(filename.end()-4, filename.end());
  
	/*--- The adjoint problem requires a particular extension. ---*/
	if (config->GetAdjoint()) {
		switch (Kind_ObjFunc) {
      case DRAG_COEFFICIENT:      AdjExt = "_cd";   break;
      case LIFT_COEFFICIENT:      AdjExt = "_cl";   break;
      case SIDEFORCE_COEFFICIENT: AdjExt = "_csf";  break;
      case PRESSURE_COEFFICIENT:  AdjExt = "_cp";   break;
      case MOMENT_X_COEFFICIENT:  AdjExt = "_cmx";  break;
      case MOMENT_Y_COEFFICIENT:  AdjExt = "_cmy";  break;
      case MOMENT_Z_COEFFICIENT:  AdjExt = "_cmz";  break;
      case EFFICIENCY:            AdjExt = "_eff";  break;
      case EQUIVALENT_AREA:       AdjExt = "_ea";   break;
      case NEARFIELD_PRESSURE:    AdjExt = "_nfp";  break;
      case FORCE_X_COEFFICIENT:   AdjExt = "_cfx";  break;
      case FORCE_Y_COEFFICIENT:   AdjExt = "_cfy";  break;
      case FORCE_Z_COEFFICIENT:   AdjExt = "_cfz";  break;
      case THRUST_COEFFICIENT:    AdjExt = "_ct";   break;
      case TORQUE_COEFFICIENT:    AdjExt = "_cq";   break;
      case FIGURE_OF_MERIT:       AdjExt = "_merit";break;
      case FREE_SURFACE:          AdjExt = "_fs";   break;
      case NOISE:                 AdjExt = "_fwh";  break;
      case HEAT_LOAD:             AdjExt = "_Q";  break;
		}
		filename.append(AdjExt);
	}
  
	/*--- Multi-zone restart files. ---*/
	if (nZone > 1 && !(config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
		filename.erase(filename.end()-4, filename.end());
		sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	}
  
	/*--- For the unsteady adjoint, we integrate backwards through
   physical time, so load in the direct solution files in reverse. ---*/
	if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
		flowIter = val_iZone;
		if (int(val_iZone) < 10) sprintf (buffer, "_0000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 10) && (int(val_iZone) < 100)) sprintf (buffer, "_000%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 100) && (int(val_iZone) < 1000)) sprintf (buffer, "_00%d.dat", int(val_iZone));
		if ((int(val_iZone) >= 1000) && (int(val_iZone) < 10000)) sprintf (buffer, "_0%d.dat", int(val_iZone));
		if (int(val_iZone) >= 10000) sprintf (buffer, "_%d.dat", int(val_iZone));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	} else if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
		flowIter   = config->GetExtIter();
		if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
		if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
		if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
		if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
		if (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
		UnstExt = string(buffer);
		filename.append(UnstExt);
	} else
    filename.append(".dat");
  
  
  restart_file.open(filename.data(), ios::in);
  
	/*--- In case there is no file ---*/
	if (restart_file.fail()) {
		cout << "There is no SU2 restart file!!" << endl;
		cout << "Press any key to exit..." << endl;
		cin.get(); exit(1);
	}
  
  /*--- Output the file name to the console. ---*/
  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;
  
  /*--- Set the number of variables, one per field in the
   restart file (without including the PointID) ---*/
  nVar = config->fields.size() - 1;
  double Solution[nVar];
  
	/*--- In case this is a parallel simulation, we need to perform the
   Global2Local index transformation first. ---*/
	long *Global2Local = NULL;
	Global2Local = new long[geometry->GetGlobal_nPointDomain()];
	/*--- First, set all indices to a negative value by default ---*/
	for(iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
		Global2Local[iPoint] = -1;
	}
  
	/*--- Now fill array with the transform values only for local points ---*/
	for(iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
	}
  
	/*--- Read all lines in the restart file ---*/
	long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  
	/*--- The first line is the header ---*/
	getline (restart_file, text_line);
  
	while (getline (restart_file,text_line)) {
		istringstream point_line(text_line);
    
		/*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as
     initialized above. Otherwise, the local index for this node on the
     current processor will be returned and used to instantiate the vars. ---*/
		iPoint_Local = Global2Local[iPoint_Global];
		if (iPoint_Local >= 0) {
      
      /*--- The PointID is not stored --*/
      point_line >> index;
      
      /*--- Store the solution --*/
      for (iField = 0; iField < nVar; iField++)
        point_line >> Solution[iField];
      
      node[iPoint_Local]->SetSolution(Solution);
      
      
		}
		iPoint_Global++;
	}
  
	/*--- Close the restart file ---*/
	restart_file.close();
  
	/*--- Free memory needed for the transformation ---*/
	delete [] Global2Local;
  
}

CBaselineSolver::~CBaselineSolver(void) { }