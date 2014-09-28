/*!
 * \file solver_structure.cpp
 * \brief Main subrotuines for solving direct, adjoint and linearized problems.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.2 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

CSolver::CSolver(void) {
  
  /*--- Array initialization ---*/
  OutputHeadingNames = NULL;
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
  nOutputVariables = 0;
  
}

CSolver::~CSolver(void) {
  if( OutputHeadingNames != NULL){
    delete []OutputHeadingNames;
  }
  //  delete [] OutputHeadingNames;
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
  
#ifndef HAVE_MPI
  
  for (iVar = 0; iVar < nVar; iVar++)
    SetRes_RMS(iVar, max(EPS, sqrt(GetRes_RMS(iVar)/geometry->GetnPoint())));
  
#else
  
  int nProcessor, iProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  double *sbuf_residual, *rbuf_residual;
  unsigned long *sbuf_point, *rbuf_point, Local_nPointDomain, Global_nPointDomain;
  
  /*--- Set the L2 Norm residual in all the processors ---*/
  sbuf_residual  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;
  
  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_RMS(iVar);
  Local_nPointDomain = geometry->GetnPointDomain();
  
  
  MPI_Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  
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
  
  MPI_Allgather(sbuf_residual, nVar, MPI_DOUBLE, rbuf_residual, nVar, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(sbuf_point, nVar, MPI_UNSIGNED_LONG, rbuf_point, nVar, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
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
  double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j, weight, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, z11, z12, z13, z22, z23, z33, detR2, product;
  bool singular = false;
  
  double *cvector = new double [nDim];
  
  /*--- Loop over points of the grid ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    AuxVar_i = node[iPoint]->GetAuxVar();
    
    /*--- Inizialization of variables ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      cvector[iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      AuxVar_j = node[jPoint]->GetAuxVar();
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (fabs(weight) > EPS){
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
          cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/(weight);
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (fabs(r11) < EPS) r11 = EPS;
    r11 = sqrt(r11);
    r12 = r12/r11;
    r22 = sqrt(r22-r12*r12);
    if (fabs(r22) < EPS) r22 = EPS;
    if (nDim == 3) {
      r13 = r13/r11;
      r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
      r33 = sqrt(r33-r23*r23-r13*r13);
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (fabs(detR2) < EPS) singular = true;
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
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
    }
    
    /*--- Computation of the gradient: S*c ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      product = 0.0;
      for (jDim = 0; jDim < nDim; jDim++)
        product += Smatrix[iDim][jDim]*cvector[jDim];
      if (geometry->node[iPoint]->GetDomain())
        node[iPoint]->SetAuxVarGradient(iDim, product);
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
  bool singular = false;
  
  double **cvector = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    cvector[iVar] = new double [nDim];
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get consevative solution ---*/
    
    Solution_i = node[iPoint]->GetSolution();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      Solution_j = node[jPoint]->GetSolution();
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (fabs(weight) > EPS){
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
            cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/weight;
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (fabs(r11) < EPS) r11 = EPS;
    r11 = sqrt(r11);
    r12 = r12/(r11);
    r22 = sqrt(r22-r12*r12);
    if (fabs(r22) < EPS) r22 = EPS;
    if (nDim == 3) {
      r13 = r13/(r11);
      r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
      r33 = sqrt(r33-r23*r23-r13*r13);
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (fabs(detR2) < EPS) singular = true;
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
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

void CSolver::SetGridVel_Gradient(CGeometry *geometry, CConfig *config) {
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
    Solution_i = geometry->node[iPoint]->GetGridVel();
    
    /*--- Inizialization of variables ---*/
    for (iVar = 0; iVar < nDim; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        cvector[iVar][iDim] = 0.0;
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      Solution_j = geometry->node[jPoint]->GetGridVel();
      
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
        geometry->node[iPoint]->SetGridVel_Grad(iVar,iDim,product);
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
    Boundary = config->GetMarker_All_KindBC(iMarker);
    switch (Boundary) {
      case EULER_WALL:
      case HEAT_FLUX:
      case HEAT_FLUX_CATALYTIC:
      case HEAT_FLUX_NONCATALYTIC:
      case ISOTHERMAL:
      case ISOTHERMAL_CATALYTIC:
      case ISOTHERMAL_NONCATALYTIC:
        
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
  dave, LimK, eps1, eps2, dm, dp, du, ds, limiter, SharpEdge_Distance;
  
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
  
  /*--- Venkatakrishnan limiter ---*/
  
  if (config->GetKind_SlopeLimit() == VENKATAKRISHNAN) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetLimiterCoeff();
    eps1 = LimK*dave;
    eps2 = eps1*eps1*eps1;
    
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
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter(iVar))
          node[iPoint]->SetLimiter(iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter(iVar))
          node[jPoint]->SetLimiter(iVar, limiter);
      }
    }
  }
  
  /*--- Sharp edges limiter ---*/
  
  if (config->GetKind_SlopeLimit() == SHARP_EDGES) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    dave = config->GetRefElemLength();
    LimK = config->GetLimiterCoeff();
    eps1 = LimK*dave;
    eps2 = pow(eps1, 3.0);
    
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
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        SharpEdge_Distance = (geometry->node[iPoint]->GetSharpEdge_Distance() - config->GetSharpEdgesCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter(iVar))
          if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SetLimiter(iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        SharpEdge_Distance = (geometry->node[jPoint]->GetSharpEdge_Distance() - config->GetSharpEdgesCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter(iVar))
          if (geometry->node[jPoint]->GetDomain()) node[jPoint]->SetLimiter(iVar, limiter);
        
      }
    }
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

void CSolver::Aeroelastic(CSurfaceMovement *surface_movement, CGeometry *geometry, CConfig *config, unsigned long ExtIter) {
  
  /*--- Variables used for Aeroelastic case ---*/
  
  double Cl, Cd, Cn, Ct, Cm, Cn_rot;
  double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  double structural_solution[4]; //contains solution of typical section wing model.
  
  unsigned short iMarker, iMarker_Monitoring, Monitoring;
  string Marker_Tag, Monitoring_Tag;
  
  /*--- Loop over markers and find the ones being monitored. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if (Monitoring == YES) {
      
      /*--- Find the particular marker being monitored and get the forces acting on it. ---*/
      
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) {
          
          Cl = GetSurface_CLift(iMarker_Monitoring);
          Cd = GetSurface_CDrag(iMarker_Monitoring);
          
          /*--- For typical section wing model want the force normal to the airfoil (in the direction of the spring) ---*/
          Cn = Cl*cos(Alpha) + Cd*sin(Alpha);
          Ct = -Cl*sin(Alpha) + Cd*cos(Alpha);
          
          Cm = -1.0*GetSurface_CMz(iMarker_Monitoring);
          
          /*--- Calculate forces for the Typical Section Wing Model taking into account rotation ---*/
          
          /* --- Note that the calculation of the forces and the subsequent displacements ...
           is only correct for the airfoil that starts at the 0 degree position --- */
          
          if (config->GetKind_GridMovement(ZONE_0) == AEROELASTIC_RIGID_MOTION) {
            double Omega, dt, psi;
            dt = config->GetDelta_UnstTimeND();
            Omega  = (config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref());
            psi = Omega*(dt*ExtIter);
            
            /* --- Correct for the airfoil starting position (This is hardcoded in here) --- */
            if (Monitoring_Tag == "Airfoil1") {
              psi = psi + 0.0;
            }
            else if (Monitoring_Tag == "Airfoil2") {
              psi = psi + 2.0/3.0*PI_NUMBER;
            }
            else if (Monitoring_Tag == "Airfoil3") {
              psi = psi + 4.0/3.0*PI_NUMBER;
            }
            else
              cout << "WARNING: There is a marker that we are monitoring that doesn't match the values hardcoded above!" << endl;
            
            cout << Monitoring_Tag << " position " << psi*180.0/PI_NUMBER << " degrees. " << endl;
            
            Cn_rot = Cn*cos(psi) - Ct*sin(psi); //Note the signs are different for accounting for the AOA.
            Cn = Cn_rot;
          }
          
          //          Can compare the Cn values here, with the values in post-processing.
          //          std::cout << "Cn = " << Cn << endl;
          //          std::cout << "Cm = " << Cm << endl;
          
          /*--- Solve the aeroelastic equations for the particular marker(surface) ---*/
          SolveTypicalSectionWingModel(geometry, Cn, Cm, config, iMarker_Monitoring, structural_solution);
          
          break;
        }
      }
      
      /*--- Compute the new surface node locations ---*/
      surface_movement->AeroelasticDeform(geometry, config, ExtIter, iMarker, iMarker_Monitoring, structural_solution);
      
    }
    
  }
  
}

void CSolver::SetUpTypicalSectionWingModel(double (&PHI)[2][2],double (&lambda)[2], CConfig *config) {
  
  /*--- Retrieve values from the config file ---*/
  double w_h = config->GetAeroelastic_Frequency_Plunge();
  double w_a = config->GetAeroelastic_Frequency_Pitch();
  
  /*--- Geometrical Parameters */
  double x_a = 1.8;
  double r_a2 = 3.48;
  
  // Mass Matrix
  // double M[2][2] = {{1,x_a},{x_a,r_a2}};
  // Stiffness Matrix
  double K[2][2] = {{(w_h/w_a)*(w_h/w_a),0},{0,r_a2}};
  
  
  /* Eigenvector and Eigenvalue Matrices of the Generalized EigenValue Problem. */
  
  double LAMBDA[2][2];
  double y;
  y = sqrt(r_a2*pow(w_a,4) - 2*r_a2*pow(w_a,2)*pow(w_h,2) + r_a2*pow(w_h,4) + 4*pow(w_a,2)*pow(w_h,2)*pow(x_a,2));
  
  PHI[0][0] = (sqrt(r_a2)*y + r_a2*pow(w_a,2) - r_a2*pow(w_h,2))/(2*pow(w_h,2)*x_a);
  PHI[0][1] = -(sqrt(r_a2)*y - r_a2*pow(w_a,2) + r_a2*pow(w_h,2))/(2*pow(w_h,2)*x_a);
  PHI[1][0] = 1.0;
  PHI[1][1] = 1.0;
  
  LAMBDA[0][0] = (r_a2*pow(w_a,2) + r_a2*pow(w_h,2) - sqrt(r_a2)*y) / (2*pow(w_a,2)*(r_a2-pow(x_a,2)));
  LAMBDA[0][1] = 0;
  LAMBDA[1][0] = 0;
  LAMBDA[1][1] = (r_a2*pow(w_a,2) + r_a2*pow(w_h,2) + sqrt(r_a2)*y) / (2*pow(w_a,2)*(r_a2-pow(x_a,2)));
  
  /* Nondimesionalize the Eigenvectors such that PHI'*M*PHI = I and PHI'*K*PHI = LAMBDA */
  double temp1[2][2], temp2[2][2];
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      temp1[i][j] = 0;
      for (int k=0; k<2; k++) {
        temp1[i][j] += K[i][k]*PHI[k][j];
      }
    }
  }
  
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      temp2[i][j] = 0;
      for (int k=0; k<2; k++) {
        temp2[i][j] += PHI[k][i]*temp1[k][j]; //PHI transpose
      }
    }
  }
  
  //Modify the first column
  PHI[0][0] = 1/sqrt(temp2[0][0]/LAMBDA[0][0])*PHI[0][0];
  PHI[1][0] = 1/sqrt(temp2[0][0]/LAMBDA[0][0])*PHI[1][0];
  //Modify the second column
  PHI[0][1] = 1/sqrt(temp2[1][1]/LAMBDA[1][1])*PHI[0][1];
  PHI[1][1] = 1/sqrt(temp2[1][1]/LAMBDA[1][1])*PHI[1][1];
  
  //Eigenvalues
  lambda[0] = sqrt(LAMBDA[0][0]);
  lambda[1] = sqrt(LAMBDA[1][1]);
  
}

void CSolver::SolveTypicalSectionWingModel(CGeometry *geometry, double Cl, double Cm, CConfig *config, unsigned short iMarker, double (&displacements)[4]) {
  
  /*--- The aeroelastic model solved in this routine is the typical section wing model
   The details of the implementation can be found in J.J. Alonso "Fully-Implicit Time-Marching Aeroelastic Solutions" 1994.
   This routine is limited to 2 dimensional problems ---*/
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short nDim=geometry->GetnDim();
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      printf("\n\n   !!! Error !!!\n");
      printf("Grid movement kind Aeroelastic is only available in 2 dimensions.");
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }
  }
  
  /*--- Retrieve values from the config file ---*/
  double w_a = config->GetAeroelastic_Frequency_Pitch();
  double dt = config->GetDelta_UnstTime();
  dt = dt*w_a; //Non-dimensionalize the structural time.
  double Lref = config->GetLength_Ref();
  double b = Lref/2.0;  // airfoil semichord
  double Density_Inf  = config->GetDensity_FreeStreamND();
  double P_Inf = config->GetPressure_FreeStreamND();
  double Mach_Inf     = config->GetMach();
  double gamma = config->GetGamma();
  
  /*--- airfoil mass ratio ---*/
  double mu = 60;
  /*--- Structural Equation damping ---*/
  double xi[2] = {0.0,0.0};
  
  /*--- Flutter Speep Index ---*/
  double Vf = (Mach_Inf*sqrt(gamma*P_Inf/Density_Inf))/(b*w_a*sqrt(mu));
  
  /*--- Eigenvectors and Eigenvalues of the Generalized EigenValue Problem. ---*/
  double PHI[2][2];   // generalized eigenvectors.
  double w[2];        //generalized eigenvalues.
  SetUpTypicalSectionWingModel(PHI,w,config);
  
  /*--- Solving the Decoupled Aeroelastic Problem with second order time discretization Eq (9) ---*/
  
  /*--- Solution variables description. //x[j][i], j-entry, i-equation. // Time (n+1)->np1, n->n, (n-1)->n1 ---*/
  // This variable gets overwritten below. I'm just using this to initialize it.
  vector<vector<double> > x_np1 = config->GetAeroelastic_np1(iMarker);
  
  /*--- Values from previous movement of spring at true time step n+1
   We use this values because we are solving for delta changes not absolute changes ---*/
  vector<vector<double> > x_np1_old = config->GetAeroelastic_np1(iMarker);
  
  /*--- Values at previous timesteps. ---*/
  vector<vector<double> > x_n = config->GetAeroelastic_n(iMarker);
  vector<vector<double> > x_n1 = config->GetAeroelastic_n1(iMarker);
  
  /*--- Set up of variables used to solve the structural problem. ---*/
  double Q[2];
  double A_inv[2][2];
  double detA;
  double S1, S2;
  double RHS[2];
  double eta[2];
  double eta_dot[2];
  
  /*--- Forcing Term ---*/
  double cons = Vf*Vf/PI_NUMBER;
  double F[2] = {cons*(-Cl), cons*(2*Cm)};
  
  for (int i=0; i<2; i++) {
    Q[i] = 0;
    for (int k=0; k<2; k++) {
      Q[i] += PHI[k][i]*F[k]; //PHI transpose
    }
  }
  
  /*--- solve each decoupled equation (The inverse of the 2x2 matrix is provided) ---*/
  for (int i=0; i<2; i++) {
    /* Matrix Inverse */
    detA = 9.0/(4.0*dt*dt) + 3*w[i]*xi[i]/(dt) + w[i]*w[i];
    A_inv[0][0] = 1/detA * 3/(2.0*dt) + 2*xi[i]*w[i];
    A_inv[0][1] = 1/detA * 1;
    A_inv[1][0] = 1/detA * -w[i]*w[i];
    A_inv[1][1] = 1/detA * 3/(2.0*dt);
    
    /* Source Terms from previous iterations */
    S1 = (-4*x_n[0][i] + x_n1[0][i])/(2.0*dt);
    S2 = (-4*x_n[1][i] + x_n1[1][i])/(2.0*dt);
    
    /* Problem Right Hand Side */
    RHS[0] = -S1;
    RHS[1] = Q[i]-S2;
    
    /* Solve the equations */
    x_np1[0][i] = A_inv[0][0]*RHS[0] + A_inv[0][1]*RHS[1];
    x_np1[1][i] = A_inv[1][0]*RHS[0] + A_inv[1][1]*RHS[1];
    
    eta[i] = x_np1[0][i]-x_np1_old[0][i];  // For displacements, the change(deltas) is used.
    eta_dot[i] = x_np1[1][i]; // For velocities, absolute values are used.
  }
  
  /*--- Transform back from the generalized coordinates to get the actual displacements in plunge and pitch ---*/
  double q[2];
  double q_dot[2];
  for (int i=0; i<2; i++) {
    q[i] = 0;
    q_dot[i] = 0;
    for (int k=0; k<2; k++) {
      q[i] += PHI[i][k]*eta[k];
      q_dot[i] += PHI[i][k]*eta_dot[k];
    }
  }
  
  double dh = b*q[0];
  double dalpha = q[1];
  
  double h_dot = w_a*b*q_dot[0];
  double alpha_dot = w_a*q_dot[1];
  
  /*--- Set the solution of the structural equations ---*/
  displacements[0] = dh;
  displacements[1] = dalpha;
  displacements[2] = h_dot;
  displacements[3] = alpha_dot;
  
  /*--- Calculate the total plunge and total pitch displacements for the unsteady step by summing the displacement at each sudo time step ---*/
  double pitch, plunge;
  pitch = config->GetAeroelastic_pitch(iMarker);
  plunge = config->GetAeroelastic_plunge(iMarker);
  
  config->SetAeroelastic_pitch(iMarker ,pitch+dalpha);
  config->SetAeroelastic_plunge(iMarker ,plunge+dh/b);
  
  /*--- Set the Aeroelastic solution at time n+1. This gets update every sudo time step
   and after convering the sudo time step the solution at n+1 get moved to the solution at n
   in SetDualTime_Solver method ---*/
  
  config->SetAeroelastic_np1(iMarker, x_np1);
  
}

CBaselineSolver::CBaselineSolver(void) : CSolver() { }

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned long iPoint, index, iPoint_Global;
  long iPoint_Local;
  unsigned short iField, iVar;
  string Tag, text_line, AdjExt, UnstExt;
  unsigned long iExtIter = config->GetExtIter();
  
  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  
  /*--- Allocate the node variables ---*/
  node = new CVariable*[geometry->GetnPoint()];
  
  /*--- Restart the solution from file information ---*/
  ifstream restart_file;
  string filename;
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FlowFileName();
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetWrt_Unsteady() || config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, int(iExtIter));
  }
  
  /*--- Open the restart file ---*/
  restart_file.open(filename.data(), ios::in);
  
  /*--- In case there is no restart file ---*/
  if (restart_file.fail()) {
    cout << "SU2 flow file " << filename << " not found" << endl;
    exit(EXIT_FAILURE);
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
      
      /*--- Store the solution (starting with node coordinates) --*/
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
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  double rotMatrix[3][3], *transl, *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Solution = NULL;
  int send_to, receive_from;
  
  Solution = new double[nVar];
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new double [nBufferR_Vector];
      Buffer_Send_U = new double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
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
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        transl = config->GetPeriodicTranslate(iPeriodic_Index);
        angles = config->GetPeriodicRotation(iPeriodic_Index);
        
        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the spatial coordinates & momentum. ---*/
        if (nDim == 2) {
          /*--- Coords ---*/
          Solution[0] = (rotMatrix[0][0]*Buffer_Receive_U[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_U[1*nVertexR+iVertex] - transl[0]);
          Solution[1] = (rotMatrix[1][0]*Buffer_Receive_U[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_U[1*nVertexR+iVertex] - transl[1]);
          /*--- Momentum ---*/
          Solution[nDim+1] = (rotMatrix[0][0]*Buffer_Receive_U[(nDim+1)*nVertexR+iVertex] +
                              rotMatrix[0][1]*Buffer_Receive_U[(nDim+2)*nVertexR+iVertex]);
          Solution[nDim+2] = (rotMatrix[1][0]*Buffer_Receive_U[(nDim+1)*nVertexR+iVertex] +
                              rotMatrix[1][1]*Buffer_Receive_U[(nDim+2)*nVertexR+iVertex]);
        } else {
          /*--- Coords ---*/
          Solution[0] = (rotMatrix[0][0]*Buffer_Receive_U[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_U[1*nVertexR+iVertex] +
                         rotMatrix[0][2]*Buffer_Receive_U[2*nVertexR+iVertex] - transl[0]);
          Solution[1] = (rotMatrix[1][0]*Buffer_Receive_U[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_U[1*nVertexR+iVertex] +
                         rotMatrix[1][2]*Buffer_Receive_U[2*nVertexR+iVertex] - transl[1]);
          Solution[2] = (rotMatrix[2][0]*Buffer_Receive_U[0*nVertexR+iVertex] +
                         rotMatrix[2][1]*Buffer_Receive_U[1*nVertexR+iVertex] +
                         rotMatrix[2][2]*Buffer_Receive_U[2*nVertexR+iVertex] - transl[2]);
          /*--- Momentum ---*/
          Solution[nDim+1] = (rotMatrix[0][0]*Buffer_Receive_U[nDim+1*nVertexR+iVertex] +
                              rotMatrix[0][1]*Buffer_Receive_U[nDim+2*nVertexR+iVertex] +
                              rotMatrix[0][2]*Buffer_Receive_U[nDim+3*nVertexR+iVertex]);
          Solution[nDim+2] = (rotMatrix[1][0]*Buffer_Receive_U[nDim+1*nVertexR+iVertex] +
                              rotMatrix[1][1]*Buffer_Receive_U[nDim+2*nVertexR+iVertex] +
                              rotMatrix[1][2]*Buffer_Receive_U[nDim+3*nVertexR+iVertex]);
          Solution[nDim+3] = (rotMatrix[2][0]*Buffer_Receive_U[nDim+1*nVertexR+iVertex] +
                              rotMatrix[2][1]*Buffer_Receive_U[nDim+2*nVertexR+iVertex] +
                              rotMatrix[2][2]*Buffer_Receive_U[nDim+3*nVertexR+iVertex]);
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
  delete [] Solution;
  
}

void CBaselineSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter) {
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Restart the solution from file information ---*/
  string filename;
  unsigned long iPoint, index;
  string UnstExt, text_line, AdjExt;
  ifstream solution_file;
  unsigned short iField;
  unsigned long iExtIter = config->GetExtIter();
  
  /*--- Retrieve filename from config ---*/
  if (config->GetAdjoint()) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FlowFileName();
  }
  
  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetWrt_Unsteady() || config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    filename = config->GetUnsteady_FileName(filename, int(iExtIter));
  }
  
  /*--- Open the restart file ---*/
  solution_file.open(filename.data(), ios::in);
  
  /*--- In case there is no file ---*/
  if (solution_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no SU2 restart file!!" << endl;
    exit(EXIT_FAILURE);
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
  Global2Local = new long[geometry[ZONE_0]->GetGlobal_nPointDomain()];
  /*--- First, set all indices to a negative value by default ---*/
  for(iPoint = 0; iPoint < geometry[ZONE_0]->GetGlobal_nPointDomain(); iPoint++) {
    Global2Local[iPoint] = -1;
  }
  
  /*--- Now fill array with the transform values only for local points ---*/
  for(iPoint = 0; iPoint < geometry[ZONE_0]->GetnPointDomain(); iPoint++) {
    Global2Local[geometry[ZONE_0]->node[iPoint]->GetGlobalIndex()] = iPoint;
  }
  
  /*--- Read all lines in the restart file ---*/
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  
  /*--- The first line is the header ---*/
  getline (solution_file, text_line);
  
  while (getline (solution_file,text_line)) {
    istringstream point_line(text_line);
    
    /*--- Retrieve local index. If this node from the restart file lives
     on a different processor, the value of iPoint_Local will be -1, as
     initialized above. Otherwise, the local index for this node on the
     current processor will be returned and used to instantiate the vars. ---*/
    iPoint_Local = Global2Local[iPoint_Global];
    if (iPoint_Local >= 0) {
      
      /*--- The PointID is not stored --*/
      point_line >> index;
      
      /*--- Store the solution (starting with node coordinates) --*/
      for (iField = 0; iField < nVar; iField++)
        point_line >> Solution[iField];
      
      node[iPoint_Local]->SetSolution(Solution);
      
      
    }
    iPoint_Global++;
  }
  
  /*--- Close the restart file ---*/
  solution_file.close();
  
  /*--- Free memory needed for the transformation ---*/
  delete [] Global2Local;
  
}

CBaselineSolver::~CBaselineSolver(void) { }
