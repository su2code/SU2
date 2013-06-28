/*!
 * \file solution_structure.cpp
 * \brief Main subrotuines for solving direct, adjoint and linearized problems.
 * \author Current Development: Stanford University.
 *         Original Structure: CADES 1.0 (2009).
 * \version 1.0.
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

#include "../include/solution_structure.hpp"

CSolution::CSolution(void) {}

CSolution::~CSolution(void) {}

void CSolution::Set_MultiSolution(CGeometry *geometry, short index) {
	unsigned long iPoint;
	
	if (index > 0) {
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint]->Set_OldSolution();
	}
	if (index < 0) {
		for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint]->Set_Solution();
	}
}

void CSolution::SetResidual_RKCoeff(CGeometry *geometry, CConfig *config, unsigned short iRKStep) {
	unsigned short iVar;
	unsigned long iPoint;
	double *local_Res_Conv, *local_Res_Visc, *local_Res_Total, *Res_Visc_km1, *Res_Visc_k, RK_BetaCoeff;
	
	Res_Visc_km1 = new double [nVar];
	Res_Visc_k = new double [nVar];
	local_Res_Total = new double [nVar];
	
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		local_Res_Conv = node[iPoint]->GetResConv();
		local_Res_Visc = node[iPoint]->GetResVisc();
		
		for (iVar = 0; iVar < nVar; iVar++) {
			if (iRKStep == 0) {
				Res_Visc_k[iVar] = local_Res_Visc[iVar];
				local_Res_Total[iVar] = local_Res_Conv[iVar] + Res_Visc_k[iVar]; }
			else {
				RK_BetaCoeff = config->Get_Beta_RKStep(iRKStep);
				Res_Visc_km1[iVar] = node[iPoint]->GetRes_Visc_RK(iVar, iRKStep-1);
				Res_Visc_k[iVar] = RK_BetaCoeff*local_Res_Visc[iVar] + (1.0 - RK_BetaCoeff) * Res_Visc_km1[iVar];
				local_Res_Total[iVar] = local_Res_Conv[iVar] + Res_Visc_k[iVar]; 
			}
		}
		node[iPoint]->SetRes_Visc_RK(Res_Visc_k, iRKStep);
		node[iPoint]->SetResidual(local_Res_Total);
	}
	
	delete [] Res_Visc_km1;
	delete [] Res_Visc_k;
	delete [] local_Res_Total;
}

void CSolution::SetResidual_Total(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh) {

	unsigned short iVar, nVar = GetnVar();
	unsigned long iPoint;
	double *Res_Conv, *Res_Visc;
	
	/*--- Add the convective and the viscous residual ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
		Res_Conv = node[iPoint]->GetResConv();
		Res_Visc = node[iPoint]->GetResVisc();
		for (iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = (Res_Conv[iVar]+Res_Visc[iVar]);
		node[iPoint]->SetResidual(Residual);		
	}
}

void CSolution::SetResidual_DualTime(CGeometry *geometry, CSolution **solution_container, CConfig *config, unsigned short iRKStep, unsigned short iMesh) {
	unsigned short iVar, jVar;
	unsigned long iPoint;
	double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep, Volume;
	
	bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
	bool incompressible = config->GetIncompressible();

	/*--- loop over points ---*/
	for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) { 
		
		/*--- Solution at time n-1, n and n+1 ---*/
		U_time_nM1 = node[iPoint]->GetSolution_time_n1();
		U_time_n = node[iPoint]->GetSolution_time_n();
		U_time_nP1 = node[iPoint]->GetSolution();
		
		/*--- Volume at time n-1 and n ---*/
		Volume_nM1 = geometry->node[iPoint]->GetVolume_n1();
		Volume_n = geometry->node[iPoint]->GetVolume_n();
		Volume_nP1 = geometry->node[iPoint]->GetVolume();
		
		Volume = geometry->node[iPoint]->GetVolume();
		
		/*--- Time Step ---*/
		TimeStep = config->GetUnst_TimeStep();
		
		/*--- Compute Residual ---*/
		for(iVar = 0; iVar < nVar; iVar++)
			Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
												+  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
		
		if (incompressible) Residual[0] = 0.0;
		
		/*--- Add Residual ---*/
		node[iPoint]->AddRes_Conv(Residual);
		
		if (implicit) {
			for (iVar = 0; iVar < nVar; iVar++) {
				for (jVar = 0; jVar < nVar; jVar++)
					Jacobian_i[iVar][jVar] = 0.0;
				Jacobian_i[iVar][iVar] = 3.0/(2.0*TimeStep);
			}
			if (incompressible) Jacobian_i[0][0] = 0.0;
			Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);	
		}
		
	}
}

void CSolution::SetGrid_Movement_Residual (CGeometry *geometry, CConfig *config) {

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
			
			node[iPoint]->SubtractRes_Conv(Residual);
			node[jPoint]->AddRes_Conv(Residual);

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

			node[Point]->AddResidual(Residual);
		}
	}
}

void CSolution::InitializeJacobianStructure(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, nPoint = geometry->GetnPoint();
	unsigned long *row_ptr, *col_ind, *vneighs, index, nnz;
	unsigned short iNeigh, nNeigh;

	row_ptr = new unsigned long [nPoint+1];

	row_ptr[0] = 0;
	for (iPoint = 0; iPoint < nPoint; iPoint++) 
		row_ptr[iPoint+1] = row_ptr[iPoint]+(geometry->node[iPoint]->GetnPoint()+1); // +1 -> to include diagonal element	
		
	nnz = row_ptr[nPoint]; 

	col_ind = new unsigned long [nnz];
	vneighs = new unsigned long [MAX_NEIGHBORS];

	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		nNeigh = geometry->node[iPoint]->GetnPoint();
		for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
			vneighs[iNeigh] = geometry->node[iPoint]->GetPoint(iNeigh);  // neighbors to point iPoint
		vneighs[nNeigh] = iPoint;  // to include relation of the point with itself (matrix diagonal)
		sort(vneighs,vneighs+nNeigh+1);
		index = row_ptr[iPoint];
		for (iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {
			col_ind[index] = vneighs[iNeigh];
			index++;
		}
	}
	bool preconditioner = false;
	Jacobian.SetIndexes(nPoint, nVar, row_ptr, col_ind, nnz, preconditioner);

	delete[] vneighs;
	// Don't delete *row_ptr, *col_ind because they are asigned to the Jacobian structure
}

void CSolution::InitializeStiffMatrixStructure(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, nPoint = geometry->GetnPoint();
	unsigned long *row_ptr, *col_ind, *vneighs, index, nnz;
	unsigned short iNeigh, nNeigh;  // Probably too much, even in 3D
	
	row_ptr = new unsigned long [nPoint+1];
	
	row_ptr[0] = 0;
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		row_ptr[iPoint+1] = row_ptr[iPoint]+(geometry->node[iPoint]->GetnPoint()+1); // +1 -> to include diagonal element	
	nnz = row_ptr[nPoint]; 
	
	col_ind = new unsigned long [nnz];
	vneighs = new unsigned long [MAX_NEIGHBORS];
	
	for (iPoint = 0; iPoint < nPoint; iPoint++) {
		nNeigh = geometry->node[iPoint]->GetnPoint();
		for (iNeigh = 0; iNeigh < nNeigh; iNeigh++)
			vneighs[iNeigh] = geometry->node[iPoint]->GetPoint(iNeigh);  // neighbors to point iPoint
		vneighs[nNeigh] = iPoint;  // to include relation of the point with itself (matrix diagonal)
		sort(vneighs,vneighs+nNeigh+1);
		index = row_ptr[iPoint];
		for (iNeigh = 0; iNeigh <= nNeigh; iNeigh++) {
			col_ind[index] = vneighs[iNeigh];
			index++;
		}
	}
	bool preconditioner = false;
	StiffMatrix.SetIndexes(nPoint, nVar, row_ptr, col_ind, nnz, preconditioner);
	
	delete[] vneighs; // Don't delete *row_ptr, *col_ind because they are asigned to the Jacobian structure
}

void CSolution::SetAuxVar_Gradient_GG(CGeometry *geometry) {
	
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

void CSolution::SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
	
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
//			if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
				for (iDim = 0; iDim < nDim; iDim++)
					weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
//			}
//			else weight = 1.0;
			
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

void CSolution::SetSolution_Gradient_GG(CGeometry *geometry) {
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
					Grad_Val = Gradient[iVar][iDim] / DualArea;
					node[iPoint]->SetGradient(iVar,iDim,Grad_Val);				
				}
}

void CSolution::SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config) {
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
//				if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
					for (iDim = 0; iDim < nDim; iDim++)
						weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
//				}
//				else weight = 1.0;
				
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
}


void CSolution::SetSurface_Gradient(CGeometry *geometry, CConfig *config) {
	
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
			case EULER_WALL: case NO_SLIP_WALL:
				
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
		delete [] cvector[iVar];
	for (iDim = 0; iDim < nDim; iDim++)
		delete [] Smatrix[iDim];
	delete [] cvector;
	delete [] Smatrix;
}

void CSolution::SetAuxVar_Surface_Gradient(CGeometry *geometry, CConfig *config) {
	
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
			case EULER_WALL: case NO_SLIP_WALL:
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

void CSolution::SetResidual_Smoothing (CGeometry *geometry, unsigned short val_nSmooth, double val_smooth_coeff) {
	
	double *Residual_Old, *Residual_Sum, Residual[5], *Residual_i, *Residual_j;
	unsigned short iVar, iSmooth, nneigh, iMarker;
	unsigned long iEdge, iPoint, jPoint, iVertex;
	
	/*--- Copy the original residual ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		node[iPoint]->SetResidual_Old(node[iPoint]->GetResidual());
	
	/*--- Jacobi iterations ---*/
	for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {

		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint]->SetResidualSumZero();
		
		/*--- Loop over Interior edges ---*/
		for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
			iPoint = geometry->edge[iEdge]->GetNode(0);
			Residual_i = node[iPoint]->GetResidual();
			
			jPoint = geometry->edge[iEdge]->GetNode(1);
			Residual_j = node[jPoint]->GetResidual();
			
			/*--- Accumulate nearest neighbor residual to Res_sum for each variable ---*/
			node[iPoint]->AddResidual_Sum(Residual_j);
			node[jPoint]->AddResidual_Sum(Residual_i);
		}
		
		/*--- Loop over all mesh points (Update Residuals with averaged sum) ---*/
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			nneigh = geometry->node[iPoint]->GetnPoint();
			Residual_Sum = node[iPoint]->GetResidual_Sum();
			Residual_Old = node[iPoint]->GetResidual_Old();
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] = (Residual_Old[iVar] + val_smooth_coeff*Residual_Sum[iVar]) 
								 / (1.0 + val_smooth_coeff*double(nneigh));
			}
			node[iPoint]->SetResidual(Residual);
		}
		
		/*--- Copy boundary values ---*/
		for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				Residual_Old = node[iPoint]->GetResidual_Old();
				node[iPoint]->SetResidual(Residual_Old);
			}
		}
	}
	
}

void CSolution::SetSolution_Smoothing (CGeometry *geometry, 
									   unsigned short val_nSmooth, double val_smooth_coeff) {
	//	Description: Perform a Jacobi approximation to an implicit residual smoothing    
	
	double *Residual_Old, *Residual_Sum, *Residual;
	//	const unsigned short nVar = geometry->GetnDim()+2;
	unsigned short iVar;
	unsigned long iPoint, iEdge, Point;
	Residual = new double [nVar];
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		double *Residual = node[iPoint]->GetSolution();
		node[iPoint]->SetResidual_Old(Residual);
	}
	
	//	Jacobi iterations
	for (unsigned short iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			node[iPoint]->SetResidualSumZero();
		
		//	 Loop over Interior edges
		for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
			const unsigned long iPoint = geometry->edge[iEdge]->GetNode(0);
			double *Residual_i = node[iPoint]->GetSolution();
			
			const unsigned long jPoint = geometry->edge[iEdge]->GetNode(1);
			double *Residual_j = node[jPoint]->GetSolution();
			
			//	Accumulate nearest neighbor residual to Res_sum for each variable
			node[iPoint]->AddResidual_Sum(Residual_j);
			node[jPoint]->AddResidual_Sum(Residual_i);
		}
		
		//	Loop over all mesh points (Update Residuals with averaged sum)
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			const unsigned short nneigh = geometry->node[iPoint]->GetnPoint();
			Residual_Sum = node[iPoint]->GetResidual_Sum();
			Residual_Old = node[iPoint]->GetResidual_Old();
			for (iVar = 0; iVar < nVar; iVar++) {
				Residual[iVar] =(Residual_Old[iVar] + 
								 val_smooth_coeff*Residual_Sum[iVar])
				/(1.0 + val_smooth_coeff*double(nneigh));
			}
			node[iPoint]->SetSolution(Residual);
		}
		
		//	Copy boundary values
		for(unsigned short iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for(unsigned long iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
				Point = geometry->vertex[iMarker][iVertex]->GetNode();
				Residual_Old = node[Point]->GetResidual_Old();
				node[Point]->SetSolution(Residual_Old);
			}
		}
	}
	
	delete [] Residual;
}

void CSolution::SetSolution_Limiter(CGeometry *geometry, CConfig *config) {
	unsigned long iPoint, jPoint;
	unsigned short iNeigh, nNeigh, iVar, iDim;
	double **Gradient_i, *Coord_i, *Coord_j, diff_coord, dist_ij, r_u, r_u_ij, 
	du_max, du_min, u_ij, *Solution_i, *Solution_j, eps2, dp, dm;
	
	double dx = 0.1;
	double LimK = config->GetLimiterCoeff();

	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		if (geometry->node[iPoint]->GetDomain()) {
			Solution_i = node[iPoint]->GetSolution();
			Gradient_i = node[iPoint]->GetGradient();
			Coord_i = geometry->node[iPoint]->GetCoord();
			nNeigh = geometry->node[iPoint]->GetnPoint();
			
			for (iVar = 0; iVar < nVar; iVar++) {
				
				/*--- Find max and min value of the variable in the control volume around the mesh point ---*/
				du_max = 1.0E-8; du_min = -1.0E-8;
				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
					jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
					Solution_j = node[jPoint]->GetSolution();
					du_max = max(du_max, Solution_j[iVar] - Solution_i[iVar]);
					du_min = min(du_min, Solution_j[iVar] - Solution_i[iVar]);
				}
				
				r_u = 1.0;
				for (iNeigh = 0; iNeigh < nNeigh; iNeigh++) {
					
					/*--- Unsconstrained reconstructed solution ---*/
					jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
					Solution_j = node[jPoint]->GetSolution();
					Coord_j = geometry->node[jPoint]->GetCoord();
					u_ij = Solution_i[iVar]; dist_ij = 0;
					for (iDim = 0; iDim < nDim; iDim++) {
						diff_coord = Coord_j[iDim]-Coord_i[iDim];
						u_ij += 0.5*diff_coord*Gradient_i[iVar][iDim];
					}
					
					/*--- Barth and Jespersen limiter ---*/
/*					r_u_ij = 1.0;
					if ((u_ij - Solution_i[iVar]) > 0.0) r_u_ij = min(1.0, du_max / (u_ij - Solution_i[iVar]));
					if ((u_ij - Solution_i[iVar]) < 0.0) r_u_ij = min(1.0, du_min / (u_ij - Solution_i[iVar]));*/
					
					/*--- Venkatakrishnan limiter ---*/
					if ((u_ij - Solution_i[iVar]) >= 0.0) dp = du_max;
					else	dp = du_min;
					dm = u_ij - Solution_i[iVar];
					
					eps2 = LimK * (dx * dx * dx);
					r_u_ij = (dp*dp+2.0*dm*dp + eps2)/(dp*dp+2*dm*dm+dm*dp + eps2);
					
					/*--- Take the smalest value of the limiter ---*/
					r_u = min(r_u, r_u_ij);
					
				}
				
				node[iPoint]->SetLimiter(iVar, r_u);
			} 
		}
	
/*	for (unsigned short iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_Boundary(iMarker) == INTERFACE_BOUNDARY) 
			for(unsigned long iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iVar = 0; iVar < nVar; iVar++)
					node[iPoint]->SetLimiter(iVar, 0.0);
			}*/
	
	
}

void CSolution::SetPressureLaplacian(CGeometry *geometry, double *PressureLaplacian) {
	
	/*---	Internal variables ---*/
	unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
	unsigned short iMarker, iVar;
	
	double NonLinVar_Vertex, NonLinVar_i, NonLinVar_j, PressVar_i, PressVar_j;
	double DualArea, Partial_Res, *Normal, Area;
	double **NonLinVar_Gradient, **PressVar_Gradient, **PressVar_Laplacian, *iCoord, *jCoord, Distance, Mean_Grad[3];
	
	NonLinVar_Gradient = new double* [geometry->GetnPoint()];
	PressVar_Gradient = new double* [geometry->GetnPoint()];
	PressVar_Laplacian = new double* [geometry->GetnPoint()];
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		NonLinVar_Gradient[iPoint] = new double [4];
		PressVar_Gradient[iPoint] = new double [2];
		PressVar_Laplacian[iPoint] = new double [2];
	}
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		for(iVar = 0; iVar < 4; iVar++)
			NonLinVar_Gradient[iPoint][iVar] = 0.0;
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		for(iVar = 0; iVar < 2; iVar++)
			PressVar_Gradient[iPoint][iVar] = 0.0;
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) 
		for(iVar = 0; iVar < 2; iVar++)
			PressVar_Laplacian[iPoint][iVar] = 0.0;
	
#ifndef CHECK
	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
		
		NonLinVar_i = 1.0*node[iPoint]->GetSolution(1)*node[iPoint]->GetSolution(1);
		NonLinVar_j = 1.0*node[jPoint]->GetSolution(1)*node[jPoint]->GetSolution(1);
		Partial_Res =  0.5 * ( NonLinVar_i + NonLinVar_j) * Normal[0];
		NonLinVar_Gradient[iPoint][0] += Partial_Res;
		NonLinVar_Gradient[jPoint][0] -= Partial_Res;
		
		NonLinVar_i = 1.0*node[iPoint]->GetSolution(1)*node[iPoint]->GetSolution(2);
		NonLinVar_j = 1.0*node[jPoint]->GetSolution(1)*node[jPoint]->GetSolution(2);
		Partial_Res =  0.5 * ( NonLinVar_i + NonLinVar_j) * Normal[1];
		NonLinVar_Gradient[iPoint][1] += Partial_Res;
		NonLinVar_Gradient[jPoint][1] -= Partial_Res;
		
		NonLinVar_i = 1.0*node[iPoint]->GetSolution(2)*node[iPoint]->GetSolution(1);
		NonLinVar_j = 1.0*node[jPoint]->GetSolution(2)*node[jPoint]->GetSolution(1);
		Partial_Res =  0.5 * ( NonLinVar_i + NonLinVar_j) * Normal[0];
		NonLinVar_Gradient[iPoint][2] += Partial_Res;
		NonLinVar_Gradient[jPoint][2] -= Partial_Res;
		
		NonLinVar_i = 1.0*node[iPoint]->GetSolution(2)*node[iPoint]->GetSolution(2);
		NonLinVar_j = 1.0*node[jPoint]->GetSolution(2)*node[jPoint]->GetSolution(2);
		Partial_Res =  0.5 * ( NonLinVar_i + NonLinVar_j) * Normal[1];
		NonLinVar_Gradient[iPoint][3] += Partial_Res;
		NonLinVar_Gradient[jPoint][3] -= Partial_Res;
		
	}
	
	/*---	Loop boundary edges ---*/
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);

			NonLinVar_Vertex = 1.0*node[Point]->GetSolution(1)*node[Point]->GetSolution(1);
			Partial_Res =  NonLinVar_Vertex * Normal[0];
			NonLinVar_Gradient[Point][0] -= Partial_Res;
			
			NonLinVar_Vertex = 1.0*node[Point]->GetSolution(1)*node[Point]->GetSolution(2);
			Partial_Res =  NonLinVar_Vertex * Normal[1];
			NonLinVar_Gradient[Point][1] -= Partial_Res;
			
			NonLinVar_Vertex = 1.0*node[Point]->GetSolution(2)*node[Point]->GetSolution(1);
			Partial_Res =  NonLinVar_Vertex * Normal[0];
			NonLinVar_Gradient[Point][2] -= Partial_Res;
			
			NonLinVar_Vertex = 1.0*node[Point]->GetSolution(2)*node[Point]->GetSolution(2);
			Partial_Res =  NonLinVar_Vertex * Normal[1];
			NonLinVar_Gradient[Point][3] -= Partial_Res;
		}
	
	/*---	Divide by the volume ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		DualArea = geometry->node[iPoint]->GetVolume();
		for(iVar = 0; iVar < 4; iVar++) {
			NonLinVar_Gradient[iPoint][iVar] /= DualArea;
		}
	}
		
	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		iCoord = geometry->node[iPoint]->GetCoord();
		jCoord = geometry->node[jPoint]->GetCoord();

		Normal = geometry->edge[iEdge]->GetNormal();	

		Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
		Distance = sqrt((iCoord[0]-jCoord[0])*(iCoord[0]-jCoord[0])+(iCoord[1]-jCoord[1])*(iCoord[1]-jCoord[1]));

		
		PressVar_i = NonLinVar_Gradient[iPoint][0] + NonLinVar_Gradient[iPoint][1];
		PressVar_j = NonLinVar_Gradient[jPoint][0] + NonLinVar_Gradient[jPoint][1];
		Mean_Grad[0] =  -0.5 * ( PressVar_i + PressVar_j);
		
	
		PressVar_i = NonLinVar_Gradient[iPoint][2] + NonLinVar_Gradient[iPoint][3];
		PressVar_j = NonLinVar_Gradient[jPoint][2] + NonLinVar_Gradient[jPoint][3];
		Mean_Grad[1] =  -0.5 * ( PressVar_i + PressVar_j);
		


		/*--- Compute vector going from iPoint to jPoint ---*/
		double Dist2 = 0; double ProjVector = 0; double Edge_Vector[3], Proj_Mean_Grad_Corrected;
		
		for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			Edge_Vector[iDim] = jCoord[iDim]-iCoord[iDim];
			Dist2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			ProjVector += Edge_Vector[iDim]*Normal[iDim];
		}
		ProjVector = ProjVector/Dist2;

		double NormalDerivative;
		
		NormalDerivative = (node[jPoint]->GetSolution(0)-node[iPoint]->GetSolution(0))/(Dist2*Dist2);
		
//		double TermA = (1.0*node[jPoint]->GetSolution(1)*node[jPoint]->GetSolution(1) - 1.0*node[iPoint]->GetSolution(1)*node[iPoint]->GetSolution(1))/(Edge_Vector[0] + EPS );
//		double TermB = (1.0*node[jPoint]->GetSolution(1)*node[jPoint]->GetSolution(2) - 1.0*node[iPoint]->GetSolution(1)*node[iPoint]->GetSolution(2))/(Edge_Vector[1] + EPS );
//		double TermC = (1.0*node[jPoint]->GetSolution(2)*node[jPoint]->GetSolution(1) - 1.0*node[iPoint]->GetSolution(2)*node[iPoint]->GetSolution(1))/(Edge_Vector[0] + EPS );
//		double TermD = (1.0*node[jPoint]->GetSolution(2)*node[jPoint]->GetSolution(2) - 1.0*node[iPoint]->GetSolution(2)*node[iPoint]->GetSolution(2))/(Edge_Vector[1] + EPS );

//		double NormalDerivative_New = ((TermA + TermB)*Normal[0] + (TermC + TermD)*Normal[1]) / Area ;
		double NormalDerivative_Press = (node[jPoint]->GetSolution(0)-node[iPoint]->GetSolution(0))/(Dist2*Dist2);
																																																																										
		NormalDerivative = NormalDerivative_Press;

		double Proj_Mean_Grad_Kappa = 0.0; double Proj_Mean_Grad_Edge = 0.0;
		for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			Proj_Mean_Grad_Kappa += Mean_Grad[iDim]*Normal[iDim];
			Proj_Mean_Grad_Edge += Mean_Grad[iDim]*Edge_Vector[iDim];
		}
		
		Proj_Mean_Grad_Corrected = Proj_Mean_Grad_Kappa;
    Proj_Mean_Grad_Corrected -= Proj_Mean_Grad_Edge*ProjVector - NormalDerivative*ProjVector*(Dist2*Dist2);
				
		Partial_Res =  Proj_Mean_Grad_Corrected ;
		
		PressVar_Laplacian[iPoint][0] += Partial_Res;
		PressVar_Laplacian[jPoint][0] -= Partial_Res;
		
	}
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		PressureLaplacian[iPoint] = PressVar_Laplacian[iPoint][0];
	}


#endif
	
#ifdef CHECK_1
	
	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();	
		
		PressVar_i = node[iPoint]->GetSolution(0);
		PressVar_j = node[jPoint]->GetSolution(0);
		
		Partial_Res =  0.5 * ( PressVar_i + PressVar_j) * Normal[0];
		PressVar_Gradient[iPoint][0] += Partial_Res;
		PressVar_Gradient[jPoint][0] -= Partial_Res;
		
		Partial_Res =  0.5 * ( PressVar_i + PressVar_j) * Normal[1];
		PressVar_Gradient[iPoint][1] += Partial_Res;
		PressVar_Gradient[jPoint][1] -= Partial_Res;
		
	}
	
	/*---	Loop boundary edges ---*/
	for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
		for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
			Point = geometry->vertex[iMarker][iVertex]->GetNode();
			Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
			
			PressVar_Vertex = node[Point]->GetSolution(0);
			
			Partial_Res =  PressVar_Vertex * Normal[0];
			PressVar_Gradient[Point][0] -= Partial_Res;
			
			Partial_Res =  PressVar_Vertex * Normal[1];
			PressVar_Gradient[Point][1] -= Partial_Res;
			
		}
	
	/*---	Divide by the volume ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		DualArea = geometry->node[iPoint]->GetVolume();
		for(iVar = 0; iVar < 2; iVar++) {
			PressVar_Gradient[iPoint][iVar] /= DualArea;
		}
	}
	
	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		Normal = geometry->edge[iEdge]->GetNormal();
		Area = sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
		
		iCoord = geometry->node[iPoint]->GetCoord();
		jCoord = geometry->node[jPoint]->GetCoord();

		/*--- Compute vector going from iPoint to jPoint ---*/
		double Dist2 = 0; double ProjVector = 0; double Edge_Vector[3], Mean_Grad[3], Proj_Mean_Grad_Corrected;
		
		for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			Edge_Vector[iDim] = jCoord[iDim]-iCoord[iDim];
			Dist2 += Edge_Vector[iDim]*Edge_Vector[iDim];
			ProjVector += Edge_Vector[iDim]*Normal[iDim];
		}
		ProjVector = ProjVector/Dist2;
		
		double Proj_Mean_Grad_Kappa = 0.0; double Proj_Mean_Grad_Edge = 0.0;
		for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
			Mean_Grad[iDim] = 0.5*(PressVar_Gradient[iPoint][iDim] + PressVar_Gradient[jPoint][iDim]);
			Proj_Mean_Grad_Kappa += Mean_Grad[iDim]*Normal[iDim];
			Proj_Mean_Grad_Edge += Mean_Grad[iDim]*Edge_Vector[iDim];
		}
		
		Proj_Mean_Grad_Corrected = Proj_Mean_Grad_Kappa;
    Proj_Mean_Grad_Corrected -= Proj_Mean_Grad_Edge*ProjVector - (node[jPoint]->GetSolution(0)-node[iPoint]->GetSolution(0))*ProjVector;

		Partial_Res =  Proj_Mean_Grad_Corrected ;
		
		PressVar_Laplacian[iPoint][0] += Partial_Res;
		PressVar_Laplacian[jPoint][0] -= Partial_Res;
				
	}
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		PressureLaplacian[iPoint] = PressVar_Laplacian[iPoint][0];
	}
	
	
#endif
	
	
#ifdef CHECK_2
	
	/*---	Loop interior edges ---*/ 
	for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
		iPoint = geometry->edge[iEdge]->GetNode(0);
		jPoint = geometry->edge[iEdge]->GetNode(1);
		
		iCoord = geometry->node[iPoint]->GetCoord();
		jCoord = geometry->node[jPoint]->GetCoord();
		
		Normal = geometry->edge[iEdge]->GetNormal();
		
		Area = 1.0; //sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]);
		Distance = 1.0; //sqrt((iCoord[0]-jCoord[0])*(iCoord[0]-jCoord[0])+(iCoord[1]-jCoord[1])*(iCoord[1]-jCoord[1]));
		
		PressVar_i = node[iPoint]->GetSolution(0);
		PressVar_j = node[jPoint]->GetSolution(0);
		
		Partial_Res =  Area * ( PressVar_j - PressVar_i ) / Distance ;
		
		PressVar_Laplacian[iPoint][0] += Partial_Res;
		PressVar_Laplacian[jPoint][0] -= Partial_Res;

	}
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		PressureLaplacian[iPoint] = PressVar_Laplacian[iPoint][0];
	}
	
#endif
	
#ifdef CHECK_2
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		PressureLaplacian[iPoint] = 0.0;
	}
#endif

	
	
	
/*	double *PressureLaplacian_Old, *PressureLaplacian_Sum, PressureLaplacian_i, PressureLaplacian_j;
	unsigned short iSmooth, nneigh;
	unsigned short nSmooth = 3;
	double smooth_coeff = 0.75;
	
	PressureLaplacian_Old = new double [geometry->GetnPoint()];
	PressureLaplacian_Sum = new double [geometry->GetnPoint()];
	
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
		PressureLaplacian_Old[iPoint]=PressureLaplacian[iPoint];
	
	for (iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
			PressureLaplacian_Old[iPoint] = 0.0;
		
		for(iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {	
			iPoint = geometry->edge[iEdge]->GetNode(0);
			PressureLaplacian_i = PressureLaplacian[iPoint];
			
			jPoint = geometry->edge[iEdge]->GetNode(1);
			PressureLaplacian_j = PressureLaplacian[jPoint];
			
			PressureLaplacian_Sum[iPoint] += PressureLaplacian_j;
			PressureLaplacian_Sum[jPoint] += PressureLaplacian_i;
		}
		
		for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
			nneigh = geometry->node[iPoint]->GetnPoint();
			
			PressureLaplacian[iPoint]= (PressureLaplacian_Old[iPoint] + smooth_coeff*PressureLaplacian_Sum[iPoint]) 
			/ (1.0 + smooth_coeff*double(nneigh));
			
		}
		
		for(iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
			for(iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++)
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				PressureLaplacian[iPoint]= PressureLaplacian_Old[iPoint];
		}
	}*/
	
	
	
	
	
	
	
	
	for(iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		delete[] NonLinVar_Gradient[iPoint];
		delete[] PressVar_Gradient[iPoint];
		delete[] PressVar_Laplacian[iPoint];
	}
	
	delete[] NonLinVar_Gradient;
	delete[] PressVar_Gradient;
	delete[] PressVar_Laplacian;

}
