/*!
 * \file solver_structure.cpp
 * \brief Main subrotuines for solving direct, adjoint and linearized problems.
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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
  Residual_RMS       = NULL;
  Residual_Max       = NULL;
  Residual_BGS       = NULL;
  Residual_Max_BGS   = NULL;
  Residual           = NULL;
  Residual_i         = NULL;
  Residual_j         = NULL;
  Point_Max          = NULL;
  Point_Max_Coord    = NULL;
  Point_Max_BGS      = NULL;
  Point_Max_Coord_BGS = NULL;
  Solution           = NULL;
  Solution_i         = NULL;
  Solution_j         = NULL;
  Vector             = NULL;
  Vector_i           = NULL;
  Vector_j           = NULL;
  Res_Conv           = NULL;
  Res_Visc           = NULL;
  Res_Sour           = NULL;
  Res_Conv_i         = NULL;
  Res_Visc_i         = NULL;
  Res_Conv_j         = NULL;
  Res_Visc_j         = NULL;
  Jacobian_i         = NULL;
  Jacobian_j         = NULL;
  Jacobian_ii        = NULL;
  Jacobian_ij        = NULL;
  Jacobian_ji        = NULL;
  Jacobian_jj        = NULL;
  Smatrix            = NULL;
  Cvector            = NULL;
  Restart_Vars       = NULL;
  Restart_Data       = NULL;
  node               = NULL;
  nOutputVariables   = 0;
  
}

CSolver::~CSolver(void) {

  unsigned short iVar, iDim;
  unsigned long iPoint;
  
  /*--- Public variables, may be accessible outside ---*/

  if ( OutputHeadingNames != NULL) {
    delete [] OutputHeadingNames;
  }

  if (node != NULL) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      delete node[iPoint];
    }
    delete [] node;
  }

  /*--- Private ---*/

  if (Residual_RMS != NULL) delete [] Residual_RMS;
  if (Residual_Max != NULL) delete [] Residual_Max;
  if (Residual != NULL) delete [] Residual;
  if (Residual_i != NULL) delete [] Residual_i;
  if (Residual_j != NULL) delete [] Residual_j;
  if (Point_Max != NULL) delete [] Point_Max;

  if (Residual_BGS != NULL) delete [] Residual_BGS;
  if (Residual_Max_BGS != NULL) delete [] Residual_Max_BGS;
  if (Point_Max_BGS != NULL) delete [] Point_Max_BGS;

  if (Point_Max_Coord != NULL) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Point_Max_Coord[iVar];
    }
    delete [] Point_Max_Coord;
  }

  if (Point_Max_Coord_BGS != NULL) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Point_Max_Coord_BGS[iVar];
    }
    delete [] Point_Max_Coord_BGS;
  }

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


  if (Jacobian_i != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
    delete [] Jacobian_i;
  }

  if (Jacobian_j != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_j[iVar];
    delete [] Jacobian_j;
  }

  if (Jacobian_ii != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ii[iVar];
    delete [] Jacobian_ii;
  }

  if (Jacobian_ij != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ij[iVar];
    delete [] Jacobian_ij;
  }

  if (Jacobian_ji != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ji[iVar];
    delete [] Jacobian_ji;
  }

  if (Jacobian_jj != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_jj[iVar];
    delete [] Jacobian_jj;
  }

  if (Smatrix != NULL) {
    for (iDim = 0; iDim < nDim; iDim++)
      delete [] Smatrix[iDim];
    delete [] Smatrix;
  }

  if (Cvector != NULL) {
    for (iVar = 0; iVar < nVarGrad; iVar++)
      delete [] Cvector[iVar];
    delete [] Cvector;
  }

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;

}

void CSolver::SetResidual_RMS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;
  
#ifndef HAVE_MPI
  
  for (iVar = 0; iVar < nVar; iVar++) {
    
    if (GetRes_RMS(iVar) != GetRes_RMS(iVar)) {
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      exit(EXIT_FAILURE);
    }

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/geometry->GetnPoint())));
    
  }
  
#else
  
  int nProcessor, iProcessor, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  su2double *sbuf_residual, *rbuf_residual, *sbuf_coord, *rbuf_coord, *Coord;
  unsigned long *sbuf_point, *rbuf_point, Local_nPointDomain, Global_nPointDomain;
  unsigned short iDim;
  
  /*--- Set the L2 Norm residual in all the processors ---*/
  
  sbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;
  
  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_RMS(iVar);
  Local_nPointDomain = geometry->GetnPointDomain();
  
  
  SU2_MPI::Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  
  
  for (iVar = 0; iVar < nVar; iVar++) {
    
    if (rbuf_residual[iVar] != rbuf_residual[iVar]) {
      
      if (rank == MASTER_NODE)
        cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
      
    }
    
    SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf_residual[iVar]/Global_nPointDomain)));
    
  }
  
  delete [] sbuf_residual;
  delete [] rbuf_residual;
  
  /*--- Set the Maximum residual in all the processors ---*/
  sbuf_residual = new su2double [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  sbuf_point = new unsigned long [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_point[iVar] = 0;
  sbuf_coord = new su2double[nVar*nDim]; for (iVar = 0; iVar < nVar*nDim; iVar++) sbuf_coord[iVar] = 0.0;
  
  rbuf_residual = new su2double [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_residual[iVar] = 0.0;
  rbuf_point = new unsigned long [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_point[iVar] = 0;
  rbuf_coord = new su2double[nProcessor*nVar*nDim]; for (iVar = 0; iVar < nProcessor*nVar*nDim; iVar++) rbuf_coord[iVar] = 0.0;

  for (iVar = 0; iVar < nVar; iVar++) {
    sbuf_residual[iVar] = GetRes_Max(iVar);
    sbuf_point[iVar] = GetPoint_Max(iVar);
    Coord = GetPoint_Max_Coord(iVar);
    for (iDim = 0; iDim < nDim; iDim++)
      sbuf_coord[iVar*nDim+iDim] = Coord[iDim];
  }
  
  SU2_MPI::Allgather(sbuf_residual, nVar, MPI_DOUBLE, rbuf_residual, nVar, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_point, nVar, MPI_UNSIGNED_LONG, rbuf_point, nVar, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_coord, nVar*nDim, MPI_DOUBLE, rbuf_coord, nVar*nDim, MPI_DOUBLE, MPI_COMM_WORLD);

  for (iVar = 0; iVar < nVar; iVar++) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      AddRes_Max(iVar, rbuf_residual[iProcessor*nVar+iVar], rbuf_point[iProcessor*nVar+iVar], &rbuf_coord[iProcessor*nVar*nDim+iVar*nDim]);
    }
  }
  
  delete [] sbuf_residual;
  delete [] rbuf_residual;
  
  delete [] sbuf_point;
  delete [] rbuf_point;
  
  delete [] sbuf_coord;
  delete [] rbuf_coord;
  
#endif
  
}

void CSolver::SetResidual_BGS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;

#ifndef HAVE_MPI

  for (iVar = 0; iVar < nVar; iVar++) {

    if (GetRes_BGS(iVar) != GetRes_BGS(iVar)) {
      cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;
      exit(EXIT_FAILURE);
    }

    SetRes_BGS(iVar, max(EPS*EPS, sqrt(GetRes_BGS(iVar)/geometry->GetnPoint())));

  }

#else

  int nProcessor, iProcessor, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  su2double *sbuf_residual, *rbuf_residual, *sbuf_coord, *rbuf_coord, *Coord;
  unsigned long *sbuf_point, *rbuf_point, Local_nPointDomain, Global_nPointDomain;
  unsigned short iDim;

  /*--- Set the L2 Norm residual in all the processors ---*/

  sbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;

  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_BGS(iVar);
  Local_nPointDomain = geometry->GetnPointDomain();


  SU2_MPI::Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


  for (iVar = 0; iVar < nVar; iVar++) {

    if (rbuf_residual[iVar] != rbuf_residual[iVar]) {

      if (rank == MASTER_NODE)
        cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << endl;

      MPI_Abort(MPI_COMM_WORLD,1);

    }

    SetRes_BGS(iVar, max(EPS*EPS, sqrt(rbuf_residual[iVar]/Global_nPointDomain)));

  }

  delete [] sbuf_residual;
  delete [] rbuf_residual;

  /*--- Set the Maximum residual in all the processors ---*/
  sbuf_residual = new su2double [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  sbuf_point = new unsigned long [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_point[iVar] = 0;
  sbuf_coord = new su2double[nVar*nDim]; for (iVar = 0; iVar < nVar*nDim; iVar++) sbuf_coord[iVar] = 0.0;

  rbuf_residual = new su2double [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_residual[iVar] = 0.0;
  rbuf_point = new unsigned long [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_point[iVar] = 0;
  rbuf_coord = new su2double[nProcessor*nVar*nDim]; for (iVar = 0; iVar < nProcessor*nVar*nDim; iVar++) rbuf_coord[iVar] = 0.0;

  for (iVar = 0; iVar < nVar; iVar++) {
    sbuf_residual[iVar] = GetRes_Max_BGS(iVar);
    sbuf_point[iVar] = GetPoint_Max_BGS(iVar);
    Coord = GetPoint_Max_Coord_BGS(iVar);
    for (iDim = 0; iDim < nDim; iDim++)
      sbuf_coord[iVar*nDim+iDim] = Coord[iDim];
  }

  SU2_MPI::Allgather(sbuf_residual, nVar, MPI_DOUBLE, rbuf_residual, nVar, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_point, nVar, MPI_UNSIGNED_LONG, rbuf_point, nVar, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_coord, nVar*nDim, MPI_DOUBLE, rbuf_coord, nVar*nDim, MPI_DOUBLE, MPI_COMM_WORLD);

  for (iVar = 0; iVar < nVar; iVar++) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      AddRes_Max_BGS(iVar, rbuf_residual[iProcessor*nVar+iVar], rbuf_point[iProcessor*nVar+iVar], &rbuf_coord[iProcessor*nVar*nDim+iVar*nDim]);
    }
  }

  delete [] sbuf_residual;
  delete [] rbuf_residual;

  delete [] sbuf_point;
  delete [] rbuf_point;

  delete [] sbuf_coord;
  delete [] rbuf_coord;

#endif

}

void CSolver::SetGrid_Movement_Residual (CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, nDim = geometry->GetnDim(), iVar, nVar = GetnVar(), iMarker;
  unsigned long iVertex, iEdge;
  su2double ProjGridVel, *Normal;
  
  /*--- Loop interior edges ---*/
   
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    const unsigned long iPoint = geometry->edge[iEdge]->GetNode(0);
    const unsigned long jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Solution at each edge point ---*/
    
    su2double *Solution_i = node[iPoint]->GetSolution();
    su2double *Solution_j = node[jPoint]->GetSolution();
    
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = 0.5* (Solution_i[iVar] + Solution_j[iVar]);
    
    /*--- Grid Velocity at each edge point ---*/
    
    su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
    su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
    for (iDim = 0; iDim < nDim; iDim++)
      Vector[iDim] = 0.5* (GridVel_i[iDim] + GridVel_j[iDim]);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += Vector[iDim]*Normal[iDim];
    
    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = ProjGridVel*Solution[iVar];
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      const unsigned long Point = geometry->vertex[iMarker][iVertex]->GetNode();
      
      /*--- Solution at each edge point ---*/
      
      su2double *Solution = node[Point]->GetSolution();
      
      /*--- Grid Velocity at each edge point ---*/
      
      su2double *GridVel = geometry->node[Point]->GetGridVel();
      
      /*--- Summed normal components ---*/
      
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      
      ProjGridVel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjGridVel -= GridVel[iDim]*Normal[iDim];
      
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = ProjGridVel*Solution[iVar];
      
      LinSysRes.AddBlock(Point, Residual);
    }
  }
  
}

void CSolver::SetAuxVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
  
  unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
  unsigned short nDim = geometry->GetnDim(), iDim, iMarker;
  
  su2double AuxVar_Vertex, AuxVar_i, AuxVar_j, AuxVar_Average;
  su2double *Gradient, DualArea, Partial_Res, Grad_Val, *Normal;
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    node[iPoint]->SetAuxVarGradientZero();    // Set Gradient to Zero
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    AuxVar_i = node[iPoint]->GetAuxVar();
    AuxVar_j = node[jPoint]->GetAuxVar();
    
    Normal = geometry->edge[iEdge]->GetNormal();
    AuxVar_Average =  0.5 * ( AuxVar_i + AuxVar_j);
    for (iDim = 0; iDim < nDim; iDim++) {
      Partial_Res = AuxVar_Average*Normal[iDim];
      node[iPoint]->AddAuxVarGradient(iDim, Partial_Res);
      node[jPoint]->SubtractAuxVarGradient(iDim, Partial_Res);
    }
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      Point = geometry->vertex[iMarker][iVertex]->GetNode();
      AuxVar_Vertex = node[Point]->GetAuxVar();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = AuxVar_Vertex*Normal[iDim];
        node[Point]->SubtractAuxVarGradient(iDim, Partial_Res);
      }
    }
  
  for (iPoint=0; iPoint<geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      Gradient = node[iPoint]->GetAuxVarGradient();
      DualArea = geometry->node[iPoint]->GetVolume();
      Grad_Val = Gradient[iDim]/(DualArea+EPS);
      node[iPoint]->SetAuxVarGradient(iDim, Grad_Val);
    }
   
}

void CSolver::SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iNeigh;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j, weight, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, z11, z12, z13, z22, z23, z33, detR2, product;
  bool singular = false;
  
  su2double *Cvector = new su2double [nDim];
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    AuxVar_i = node[iPoint]->GetAuxVar();
    
    /*--- Inizialization of variables ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Cvector[iDim] = 0.0;
    
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
      
      if (fabs(weight) > EPS) {
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
          Cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/(weight);
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
        product += Smatrix[iDim][jDim]*Cvector[jDim];
      if (geometry->node[iPoint]->GetDomain())
        node[iPoint]->SetAuxVarGradient(iDim, product);
    }
  }
  
  delete [] Cvector;
  
}

void CSolver::SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
  unsigned short iVar, iDim, iMarker;
  su2double *Solution_Vertex, *Solution_i, *Solution_j, Solution_Average, **Gradient, DualArea,
  Partial_Res, Grad_Val, *Normal;
  
  /*--- Set Gradient to Zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    node[iPoint]->SetGradientZero();
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Solution_i = node[iPoint]->GetSolution();
    Solution_j = node[jPoint]->GetSolution();
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar< nVar; iVar++) {
      Solution_Average =  0.5 * (Solution_i[iVar] + Solution_j[iVar]);
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = Solution_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          node[iPoint]->AddGradient(iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          node[jPoint]->SubtractGradient(iVar, iDim, Partial_Res);
      }
    }
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      Point = geometry->vertex[iMarker][iVertex]->GetNode();
      Solution_Vertex = node[Point]->GetSolution();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      for (iVar = 0; iVar < nVar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++) {
          Partial_Res = Solution_Vertex[iVar]*Normal[iDim];
          if (geometry->node[Point]->GetDomain())
            node[Point]->SubtractGradient(iVar, iDim, Partial_Res);
        }
    }
  }
  
  /*--- Compute gradient ---*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++)
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++) {
        Gradient = node[iPoint]->GetGradient();
        DualArea = geometry->node[iPoint]->GetVolume();
        Grad_Val = Gradient[iVar][iDim] / (DualArea+EPS);
        node[iPoint]->SetGradient(iVar, iDim, Grad_Val);
      }
  
  /*--- Gradient MPI ---*/
  Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iVar, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, *Solution_i, *Solution_j,
  r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight, detR2, z11, z12, z13,
  z22, z23, z33, product;
  bool singular = false;
  
  su2double **Cvector = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get consevative solution ---*/
    
    Solution_i = node[iPoint]->GetSolution();
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;

    AD::StartPreacc();
    AD::SetPreaccIn(Solution_i, nVar);
    AD::SetPreaccIn(Coord_i, nDim);

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      Solution_j = node[jPoint]->GetSolution();

      AD::SetPreaccIn(Coord_j, nDim);
      AD::SetPreaccIn(Solution_j, nVar);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        if (nDim == 3) {
          r13   += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33   += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/weight;
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
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
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        node[iPoint]->SetGradient(iVar, iDim, product);
      }
    }

    AD::SetPreaccOut(node[iPoint]->GetGradient(), nVar, nDim);
    AD::EndPreacc();
  }
  
  /*--- Deallocate memory ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Cvector[iVar];
  delete [] Cvector;
  
  /*--- Gradient MPI ---*/
  
  Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetGridVel_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iDim, jDim, iVar, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, *Solution_i, *Solution_j, Smatrix[3][3],
  r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight, detR2, z11, z12, z13,
  z22, z23, z33, product;
  su2double **Cvector;
  
  /*--- Note that all nVar entries in this routine have been changed to nDim ---*/
  Cvector = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++)
    Cvector[iVar] = new su2double [nDim];
  
  /*--- Loop over points of the grid ---*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    Solution_i = geometry->node[iPoint]->GetGridVel();
    
    /*--- Inizialization of variables ---*/
    for (iVar = 0; iVar < nDim; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
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
          Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/(weight);
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
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        geometry->node[iPoint]->SetGridVel_Grad(iVar, iDim, product);
      }
    }
  }
  
  /*--- Deallocate memory ---*/
  for (iVar = 0; iVar < nDim; iVar++)
    delete [] Cvector[iVar];
  delete [] Cvector;
  
  /*--- Gradient MPI ---*/
  // TO DO!!!
  //Set_MPI_Solution_Gradient(geometry, config);
  
}

void CSolver::SetAuxVar_Surface_Gradient(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iNeigh, iMarker, Boundary;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint, iVertex;
  su2double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j;
  su2double **Smatrix, *Cvector;
  
  Smatrix = new su2double* [nDim];
  Cvector = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
  
  
  /*--- Loop over boundary markers to select those for Euler or NS walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    switch (Boundary) {
      case EULER_WALL:
      case HEAT_FLUX:
      case ISOTHERMAL:
        
        /*--- Loop over points on the surface (Least-Squares approximation) ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            Coord_i = geometry->node[iPoint]->GetCoord();
            AuxVar_i = node[iPoint]->GetAuxVar();
            
            /*--- Inizialization of variables ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Cvector[iDim] = 0.0;
            su2double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;
            
            for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
              Coord_j = geometry->node[jPoint]->GetCoord();
              AuxVar_j = node[jPoint]->GetAuxVar();
              
              su2double weight = 0;
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
                Cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/weight;
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
              su2double detR2 = (r11*r22)*(r11*r22);
              Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
              Smatrix[0][1] = -r11*r12/detR2;
              Smatrix[1][0] = Smatrix[0][1];
              Smatrix[1][1] = r11*r11/detR2;
            }
            else {
              su2double detR2 = (r11*r22*r33)*(r11*r22*r33);
              su2double z11, z12, z13, z22, z23, z33; // aux vars
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
            su2double product;
            for (iDim = 0; iDim < nDim; iDim++) {
              product = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                product += Smatrix[iDim][jDim]*Cvector[jDim];
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
  delete [] Cvector;
  delete [] Smatrix;
}

void CSolver::SetSolution_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Solution, *Solution_i, *Solution_j, *LocalMinSolution, *LocalMaxSolution,
  *GlobalMinSolution, *GlobalMaxSolution,
  dave, LimK, eps1, eps2, dm, dp, du, ds, y, limiter, SharpEdge_Distance;
  
  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();
  
  if (config->GetKind_SlopeLimit() == NO_LIMITER) {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->SetLimiter(iVar, 1.0);
      }
    }
    
  }
  
  else {
    
    /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        node[iPoint]->SetSolution_Max(iVar, -EPS);
        node[iPoint]->SetSolution_Min(iVar, EPS);
        node[iPoint]->SetLimiter(iVar, 2.0);
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
    
  }
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nVar, nDim);
      AD::SetPreaccIn(Gradient_j, nVar, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nVar; iVar++) {
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
          else dp = node[iPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[iPoint]->GetLimiter(iVar)) {
          node[iPoint]->SetLimiter(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter()[iVar]);
        }
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
          else dp = node[jPoint]->GetSolution_Min(iVar);
          limiter = dp/dm;
        }
        
        if (limiter < node[jPoint]->GetLimiter(iVar)) {
          node[jPoint]->SetLimiter(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter()[iVar]);
        }

      }
      
      AD::EndPreacc();
      
    }


    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        y =  node[iPoint]->GetLimiter(iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        node[iPoint]->SetLimiter(iVar, limiter);
      }
    }
    
  }

  /*--- Venkatakrishnan limiter ---*/
  
  if ((config->GetKind_SlopeLimit() == VENKATAKRISHNAN) || (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {
    
    /*--- Allocate memory for the max and min solution value --*/
    
    LocalMinSolution = new su2double [nVar]; GlobalMinSolution = new su2double [nVar];
    LocalMaxSolution = new su2double [nVar]; GlobalMaxSolution = new su2double [nVar];
    
    /*--- Compute the max value and min value of the solution ---*/
    
    Solution = node[iPoint]->GetSolution();
    for (iVar = 0; iVar < nVar; iVar++) {
      LocalMinSolution[iVar] = Solution[iVar];
      LocalMaxSolution[iVar] = Solution[iVar];
    }
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      
      /*--- Get the solution variables ---*/
      
      Solution = node[iPoint]->GetSolution();
      
      for (iVar = 0; iVar < nVar; iVar++) {
        LocalMinSolution[iVar] = min (LocalMinSolution[iVar], Solution[iVar]);
        LocalMaxSolution[iVar] = max (LocalMaxSolution[iVar], Solution[iVar]);
      }
      
    }
    
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(LocalMinSolution, GlobalMinSolution, nVar, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(LocalMaxSolution, GlobalMaxSolution, nVar, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    for (iVar = 0; iVar < nVar; iVar++) {
      GlobalMinSolution[iVar] = LocalMinSolution[iVar];
      GlobalMaxSolution[iVar] = LocalMaxSolution[iVar];
    }
#endif
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = node[iPoint]->GetGradient();
      Gradient_j = node[jPoint]->GetGradient();
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nVar, nDim);
      AD::SetPreaccIn(Gradient_j, nVar, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nVar; iVar++) {
        
        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          eps1 = LimK * (GlobalMaxSolution[iVar] - GlobalMinSolution[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }
        
        AD::SetPreaccIn(node[iPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[iPoint]->GetSolution_Min(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Max(iVar));
        AD::SetPreaccIn(node[jPoint]->GetSolution_Min(iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = node[iPoint]->GetSolution_Max(iVar);
        else dp = node[iPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter(iVar)) {
          node[iPoint]->SetLimiter(iVar, limiter);
          AD::SetPreaccOut(node[iPoint]->GetLimiter()[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter(iVar)) {
          node[jPoint]->SetLimiter(iVar, limiter);
          AD::SetPreaccOut(node[jPoint]->GetLimiter()[iVar]);
        }
      }
      
      AD::EndPreacc();

    }
    
    delete [] LocalMinSolution; delete [] GlobalMinSolution;
    delete [] LocalMaxSolution; delete [] GlobalMaxSolution;

  }
  
  /*--- Sharp edges limiter ---*/
  
  if (config->GetKind_SlopeLimit() == SHARP_EDGES) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetVenkat_LimiterCoeff();
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
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[iPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter(iVar))
          node[iPoint]->SetLimiter(iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[jPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter(iVar))
          node[jPoint]->SetLimiter(iVar, limiter);
        
      }
    }
  }
  
  /*--- Sharp edges limiter ---*/
  
  if (config->GetKind_SlopeLimit() == WALL_DISTANCE) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetVenkat_LimiterCoeff();
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
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[iPoint]->GetWall_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[iPoint]->GetLimiter(iVar))
          node[iPoint]->SetLimiter(iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = node[jPoint]->GetSolution_Max(iVar);
        else dp = node[jPoint]->GetSolution_Min(iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[jPoint]->GetWall_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < node[jPoint]->GetLimiter(iVar))
          node[jPoint]->SetLimiter(iVar, limiter);
        
      }
    }
  }

  
  /*--- Limiter MPI ---*/
  
  Set_MPI_Solution_Limiter(geometry, config);
  
}

void CSolver::SetPressureLaplacian(CGeometry *geometry, CConfig *config, su2double *PressureLaplacian) {
  
  unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
  unsigned short iMarker, iVar;
  su2double DualArea, Partial_Res, *Normal;
  su2double **UxVar_Gradient, **UyVar_Gradient;
  
  UxVar_Gradient = new su2double* [geometry->GetnPoint()];
  UyVar_Gradient = new su2double* [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    UxVar_Gradient[iPoint] = new su2double [2];
    UyVar_Gradient[iPoint] = new su2double [2];
  }
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
    for (iVar = 0; iVar < 2; iVar++) {
      UxVar_Gradient[iPoint][iVar] = 0.0;
      UyVar_Gradient[iPoint][iVar] = 0.0;
    }
  
  /*---  Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    Normal = geometry->edge[iEdge]->GetNormal();
    
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
  
  /*---  Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      Point = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      
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
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    delete[] UxVar_Gradient[iPoint];
    delete[] UyVar_Gradient[iPoint];
  }
  
  delete[] UxVar_Gradient;
  delete[] UyVar_Gradient;
  
}

void CSolver::Gauss_Elimination(su2double** A, su2double* rhs, unsigned short nVar) {
  
  short iVar, jVar, kVar;
  su2double weight, aux;
  
  if (nVar == 1)
    rhs[0] /= A[0][0];
  else {
    
    /*--- Transform system in Upper Matrix ---*/
    
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = A[iVar][jVar]/A[jVar][jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          A[iVar][kVar] -= weight*A[jVar][kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }
    
    /*--- Backwards substitution ---*/
    
    rhs[nVar-1] = rhs[nVar-1]/A[nVar-1][nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += A[iVar][jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux)/A[iVar][iVar];
      if (iVar == 0) break;
    }
  }
  
}

void CSolver::Aeroelastic(CSurfaceMovement *surface_movement, CGeometry *geometry, CConfig *config, unsigned long ExtIter) {
  
  /*--- Variables used for Aeroelastic case ---*/
  
  su2double Cl, Cd, Cn, Ct, Cm, Cn_rot;
  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  vector<su2double> structural_solution(4,0.0); //contains solution(displacements and rates) of typical section wing model.
  
  unsigned short iMarker, iMarker_Monitoring, Monitoring;
  string Marker_Tag, Monitoring_Tag;
  
  /*--- Loop over markers and find the ones being monitored. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if (Monitoring == YES) {
      
      /*--- Find the particular marker being monitored and get the forces acting on it. ---*/
      
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) {
          
          Cl = GetSurface_CL(iMarker_Monitoring);
          Cd = GetSurface_CD(iMarker_Monitoring);
          
          /*--- For typical section wing model want the force normal to the airfoil (in the direction of the spring) ---*/
          Cn = Cl*cos(Alpha) + Cd*sin(Alpha);
          Ct = -Cl*sin(Alpha) + Cd*cos(Alpha);
          
          Cm = GetSurface_CMz(iMarker_Monitoring);
          
          /*--- Calculate forces for the Typical Section Wing Model taking into account rotation ---*/
          
          /*--- Note that the calculation of the forces and the subsequent displacements ...
           is only correct for the airfoil that starts at the 0 degree position ---*/
          
          if (config->GetKind_GridMovement(ZONE_0) == AEROELASTIC_RIGID_MOTION) {
            su2double Omega, dt, psi;
            dt = config->GetDelta_UnstTimeND();
            Omega  = (config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref());
            psi = Omega*(dt*ExtIter);
            
            /*--- Correct for the airfoil starting position (This is hardcoded in here) ---*/
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

void CSolver::SetUpTypicalSectionWingModel(vector<vector<su2double> >& Phi, vector<su2double>& omega, CConfig *config) {
  
  /*--- Retrieve values from the config file ---*/
  su2double w_h = config->GetAeroelastic_Frequency_Plunge();
  su2double w_a = config->GetAeroelastic_Frequency_Pitch();
  su2double x_a = config->GetAeroelastic_CG_Location();
  su2double r_a = sqrt(config->GetAeroelastic_Radius_Gyration_Squared());
  su2double w = w_h/w_a;
  
  // Mass Matrix
  vector<vector<su2double> > M(2,vector<su2double>(2,0.0));
  M[0][0] = 1;
  M[0][1] = x_a;
  M[1][0] = x_a;
  M[1][1] = r_a*r_a;
  
  // Stiffness Matrix
  //  vector<vector<su2double> > K(2,vector<su2double>(2,0.0));
  //  K[0][0] = (w_h/w_a)*(w_h/w_a);
  //  K[0][1] = 0.0;
  //  K[1][0] = 0.0;
  //  K[1][1] = r_a*r_a;
  
  /* Eigenvector and Eigenvalue Matrices of the Generalized EigenValue Problem. */
  
  vector<vector<su2double> > Omega2(2,vector<su2double>(2,0.0));
  su2double aux; // auxiliary variable
  aux = sqrt(pow(r_a,2)*pow(w,4) - 2*pow(r_a,2)*pow(w,2) + pow(r_a,2) + 4*pow(x_a,2)*pow(w,2));
  Phi[0][0] = (r_a * (r_a - r_a*pow(w,2) + aux)) / (2*x_a*pow(w, 2));
  Phi[0][1] = (r_a * (r_a - r_a*pow(w,2) - aux)) / (2*x_a*pow(w, 2));
  Phi[1][0] = 1.0;
  Phi[1][1] = 1.0;
  
  Omega2[0][0] = (r_a * (r_a + r_a*pow(w,2) - aux)) / (2*(pow(r_a, 2) - pow(x_a, 2)));
  Omega2[0][1] = 0;
  Omega2[1][0] = 0;
  Omega2[1][1] = (r_a * (r_a + r_a*pow(w,2) + aux)) / (2*(pow(r_a, 2) - pow(x_a, 2)));
  
  /* Nondimesionalize the Eigenvectors such that Phi'*M*Phi = I and PHI'*K*PHI = Omega */
  // Phi'*M*Phi = D
  // D^(-1/2)*Phi'*M*Phi*D^(-1/2) = D^(-1/2)*D^(1/2)*D^(1/2)*D^(-1/2) = I,  D^(-1/2) = inv(sqrt(D))
  // Phi = Phi*D^(-1/2)
  
  vector<vector<su2double> > Aux(2,vector<su2double>(2,0.0));
  vector<vector<su2double> > D(2,vector<su2double>(2,0.0));
  // Aux = M*Phi
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      Aux[i][j] = 0;
      for (int k=0; k<2; k++) {
        Aux[i][j] += M[i][k]*Phi[k][j];
      }
    }
  }
  
  // D = Phi'*Aux
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      D[i][j] = 0;
      for (int k=0; k<2; k++) {
        D[i][j] += Phi[k][i]*Aux[k][j]; //PHI transpose
      }
    }
  }
  
  //Modify the first column
  Phi[0][0] = Phi[0][0] * 1/sqrt(D[0][0]);
  Phi[1][0] = Phi[1][0] * 1/sqrt(D[0][0]);
  //Modify the second column
  Phi[0][1] = Phi[0][1] * 1/sqrt(D[1][1]);
  Phi[1][1] = Phi[1][1] * 1/sqrt(D[1][1]);
  
  // Sqrt of the eigenvalues (frequency of vibration of the modes)
  omega[0] = sqrt(Omega2[0][0]);
  omega[1] = sqrt(Omega2[1][1]);
  
}

void CSolver::SolveTypicalSectionWingModel(CGeometry *geometry, su2double Cl, su2double Cm, CConfig *config, unsigned short iMarker, vector<su2double>& displacements) {
  
  /*--- The aeroelastic model solved in this routine is the typical section wing model
   The details of the implementation are similar to those found in J.J. Alonso 
   "Fully-Implicit Time-Marching Aeroelastic Solutions" 1994. ---*/
  
  /*--- Retrieve values from the config file ---*/
  su2double w_alpha = config->GetAeroelastic_Frequency_Pitch();
  su2double vf      = config->GetAeroelastic_Flutter_Speed_Index();
  su2double b       = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
  su2double dt      = config->GetDelta_UnstTimeND();
  dt = dt*w_alpha; //Non-dimensionalize the structural time.
  
  /*--- Structural Equation damping ---*/
  vector<su2double> xi(2,0.0);
  
  /*--- Eigenvectors and Eigenvalues of the Generalized EigenValue Problem. ---*/
  vector<vector<su2double> > Phi(2,vector<su2double>(2,0.0));   // generalized eigenvectors.
  vector<su2double> w(2,0.0);        // sqrt of the generalized eigenvalues (frequency of vibration of the modes).
  SetUpTypicalSectionWingModel(Phi, w, config);
  
  /*--- Solving the Decoupled Aeroelastic Problem with second order time discretization Eq (9) ---*/
  
  /*--- Solution variables description. //x[j][i], j-entry, i-equation. // Time (n+1)->np1, n->n, (n-1)->n1 ---*/
  vector<vector<su2double> > x_np1(2,vector<su2double>(2,0.0));
  
  /*--- Values from previous movement of spring at true time step n+1
   We use this values because we are solving for delta changes not absolute changes ---*/
  vector<vector<su2double> > x_np1_old = config->GetAeroelastic_np1(iMarker);
  
  /*--- Values at previous timesteps. ---*/
  vector<vector<su2double> > x_n = config->GetAeroelastic_n(iMarker);
  vector<vector<su2double> > x_n1 = config->GetAeroelastic_n1(iMarker);
  
  /*--- Set up of variables used to solve the structural problem. ---*/
  vector<su2double> f_tilde(2,0.0);
  vector<vector<su2double> > A_inv(2,vector<su2double>(2,0.0));
  su2double detA;
  su2double s1, s2;
  vector<su2double> rhs(2,0.0); //right hand side
  vector<su2double> eta(2,0.0);
  vector<su2double> eta_dot(2,0.0);
  
  /*--- Forcing Term ---*/
  su2double cons = vf*vf/PI_NUMBER;
  vector<su2double> f(2,0.0);
  f[0] = cons*(-Cl);
  f[1] = cons*(2*-Cm);
  
  //f_tilde = Phi'*f
  for (int i=0; i<2; i++) {
    f_tilde[i] = 0;
    for (int k=0; k<2; k++) {
      f_tilde[i] += Phi[k][i]*f[k]; //PHI transpose
    }
  }
  
  /*--- solve each decoupled equation (The inverse of the 2x2 matrix is provided) ---*/
  for (int i=0; i<2; i++) {
    /* Matrix Inverse */
    detA = 9.0/(4.0*dt*dt) + 3*w[i]*xi[i]/(dt) + w[i]*w[i];
    A_inv[0][0] = 1/detA * (3/(2.0*dt) + 2*xi[i]*w[i]);
    A_inv[0][1] = 1/detA * 1;
    A_inv[1][0] = 1/detA * -w[i]*w[i];
    A_inv[1][1] = 1/detA * 3/(2.0*dt);
    
    /* Source Terms from previous iterations */
    s1 = (-4*x_n[0][i] + x_n1[0][i])/(2.0*dt);
    s2 = (-4*x_n[1][i] + x_n1[1][i])/(2.0*dt);
    
    /* Problem Right Hand Side */
    rhs[0] = -s1;
    rhs[1] = f_tilde[i]-s2;
    
    /* Solve the equations */
    x_np1[0][i] = A_inv[0][0]*rhs[0] + A_inv[0][1]*rhs[1];
    x_np1[1][i] = A_inv[1][0]*rhs[0] + A_inv[1][1]*rhs[1];
    
    eta[i] = x_np1[0][i]-x_np1_old[0][i];  // For displacements, the change(deltas) is used.
    eta_dot[i] = x_np1[1][i]; // For velocities, absolute values are used.
  }
  
  /*--- Transform back from the generalized coordinates to get the actual displacements in plunge and pitch  q = Phi*eta ---*/
  vector<su2double> q(2,0.0);
  vector<su2double> q_dot(2,0.0);
  for (int i=0; i<2; i++) {
    q[i] = 0;
    q_dot[i] = 0;
    for (int k=0; k<2; k++) {
      q[i] += Phi[i][k]*eta[k];
      q_dot[i] += Phi[i][k]*eta_dot[k];
    }
  }
  
  su2double dh = b*q[0];
  su2double dalpha = q[1];
  
  su2double h_dot = w_alpha*b*q_dot[0];  //The w_a brings it back to actual time.
  su2double alpha_dot = w_alpha*q_dot[1];
  
  /*--- Set the solution of the structural equations ---*/
  displacements[0] = dh;
  displacements[1] = dalpha;
  displacements[2] = h_dot;
  displacements[3] = alpha_dot;
  
  /*--- Calculate the total plunge and total pitch displacements for the unsteady step by summing the displacement at each sudo time step ---*/
  su2double pitch, plunge;
  pitch = config->GetAeroelastic_pitch(iMarker);
  plunge = config->GetAeroelastic_plunge(iMarker);
  
  config->SetAeroelastic_pitch(iMarker , pitch+dalpha);
  config->SetAeroelastic_plunge(iMarker , plunge+dh/b);
  
  /*--- Set the Aeroelastic solution at time n+1. This gets update every sudo time step
   and after convering the sudo time step the solution at n+1 get moved to the solution at n
   in SetDualTime_Solver method ---*/
  
  config->SetAeroelastic_np1(iMarker, x_np1);
  
}

void CSolver::Restart_OldGeometry(CGeometry *geometry, CConfig *config) {

  /*--- This function is intended for dual time simulations ---*/

  unsigned long index;

  int Unst_RestartIter;
  ifstream restart_file_n;
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();
  string filename = config->GetSolution_FlowFileName();
  string filename_n;

  /*--- Auxiliary vector for storing the coordinates ---*/
  su2double *Coord;
  Coord = new su2double[nDim];

  /*--- Variables for reading the restart files ---*/
  string text_line;
  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0;
  unsigned short rbuf_NotMatching, sbuf_NotMatching;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- First, we load the restart file for time n ---*/

  /*-------------------------------------------------------------------------------------------*/

  /*--- Modify file name for an unsteady restart ---*/
  Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-1;
  filename_n = config->GetUnsteady_FileName(filename, Unst_RestartIter);

  /*--- Open the restart file, throw an error if this fails. ---*/

  restart_file_n.open(filename_n.data(), ios::in);
  if (restart_file_n.fail()) {
    if (rank == MASTER_NODE)
      cout << "There is no flow restart file!! " << filename_n.data() << "."<< endl;
    exit(EXIT_FAILURE);
  }

  /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
  iPoint_Global_Local = 0, iPoint_Global = 0;

  /*--- Read all lines in the restart file ---*/
  /*--- The first line is the header ---*/

  getline (restart_file_n, text_line);

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    getline (restart_file_n, text_line);
    
    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1];
      if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2];

      geometry->node[iPoint_Local]->SetCoord_n(Coord);

      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/

  rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  if (iPoint_Global_Local < geometry->GetnPointDomain()) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    if (rank == MASTER_NODE) {
      cout << endl << "The solution file " << filename_n.data() << " doesn't match with the mesh file!" << endl;
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

  restart_file_n.close();

  /*-------------------------------------------------------------------------------------------*/
  /*-------------------------------------------------------------------------------------------*/

  /*--- Now, we load the restart file for time n-1, if the simulation is 2nd Order ---*/

  if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND) {

    ifstream restart_file_n1;
    string filename_n1;

    /*--- Modify file name for an unsteady restart ---*/
    Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_RestartIter())-2;
    filename_n1 = config->GetUnsteady_FileName(filename, Unst_RestartIter);

    /*--- Open the restart file, throw an error if this fails. ---*/

    restart_file_n.open(filename_n1.data(), ios::in);
    if (restart_file_n.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no flow restart file!! " << filename_n1.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }

    /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
    iPoint_Global_Local = 0, iPoint_Global = 0;

    /*--- Read all lines in the restart file ---*/
    /*--- The first line is the header ---*/

    getline (restart_file_n, text_line);

    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file_n, text_line);
      
      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        if (nDim == 2) point_line >> index >> Coord[0] >> Coord[1];
        if (nDim == 3) point_line >> index >> Coord[0] >> Coord[1] >> Coord[2];

        geometry->node[iPoint_Local]->SetCoord_n1(Coord);

        iPoint_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    if (iPoint_Global_Local < geometry->GetnPointDomain()) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      if (rank == MASTER_NODE) {
        cout << endl << "The solution file " << filename_n1.data() << " doesn't match with the mesh file!" << endl;
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

    restart_file_n1.close();

  }

  /*--- It's necessary to communicate this information ---*/

  geometry->Set_MPI_OldCoord(config);
  
  delete [] Coord;

}

void CSolver::Read_SU2_Restart_ASCII(CGeometry *geometry, CConfig *config, string val_filename) {

  ifstream restart_file;
  string text_line, Tag;
  unsigned short iVar;
  long index, iPoint_Local = 0; unsigned long iPoint_Global = 0;
  int counter = 0;
  config->fields.clear();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  Restart_Vars = new int[5];

  /*--- First, check that this is not a binary restart file. ---*/

  char fname[100];
  strcpy(fname, val_filename.c_str());
  int magic_number;

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Attempt to read the first int, which should be our magic number. ---*/

  ret = fread(&magic_number, sizeof(int), 1, fhw);
  if (ret != 1) {
    cout << endl << "Error reading restart file." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
    cout << " SU2 reads/writes binary restart files by default." << endl;
    cout << " Note that backward compatibility for ASCII restart files is" << endl;
    cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
    exit(EXIT_FAILURE);
  }

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  int ierr;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    if (rank == MASTER_NODE)
      cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
  }

  /*--- Have the master attempt to read the magic number. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    if (rank == MASTER_NODE) {
      cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
  }

  MPI_File_close(&fhw);

#endif

  /*--- Open the restart file ---*/

  restart_file.open(val_filename.data(), ios::in);

  /*--- In case there is no restart file ---*/

  if (restart_file.fail()) {
    if (rank == MASTER_NODE)
      cout << "SU2 ASCII solution file " << val_filename << " not found." << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Identify the number of fields (and names) in the restart file ---*/

  getline (restart_file, text_line);
  stringstream ss(text_line);
  while (ss >> Tag) {
    config->fields.push_back(Tag);
    if (ss.peek() == ',') ss.ignore();
  }

  /*--- Set the number of variables, one per field in the
   restart file (without including the PointID) ---*/

  Restart_Vars[1] = (int)config->fields.size() - 1;

  /*--- Allocate memory for the restart data. ---*/

  Restart_Data = new passivedouble[Restart_Vars[1]*geometry->GetnPointDomain()];

  /*--- Read all lines in the restart file and extract data. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    getline (restart_file, text_line);

    istringstream point_line(text_line);

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- The PointID is not stored --*/

      point_line >> index;

      /*--- Store the solution (starting with node coordinates) --*/

      for (iVar = 0; iVar < Restart_Vars[1]; iVar++)
        point_line >> Restart_Data[counter*Restart_Vars[1] + iVar];

      /*--- Increment our local point counter. ---*/

      counter++;

    }
  }

}

void CSolver::Read_SU2_Restart_Binary(CGeometry *geometry, CConfig *config, string val_filename) {

  char str_buf[CGNS_STRING_SIZE], fname[100];
  unsigned short iVar;
  strcpy(fname, val_filename.c_str());
  int nRestart_Vars = 5, nFields;
  Restart_Vars = new int[5];
  config->fields.clear();

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- First, read the number of variables and points. ---*/

  ret = fread(Restart_Vars, sizeof(int), nRestart_Vars, fhw);
  if (ret != (unsigned long)nRestart_Vars) {
    cout << endl << "Error reading restart file." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
    cout << " SU2 reads/writes binary restart files by default." << endl;
    cout << " Note that backward compatibility for ASCII restart files is" << endl;
    cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Store the number of fields to be read for clarity. ---*/

  nFields = Restart_Vars[1];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. We pad the beginning of the
   variable string vector with the Point_ID tag that wasn't written. ---*/

  config->fields.push_back("Point_ID");
  for (iVar = 0; iVar < nFields; iVar++) {
    ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, fhw);
    if (ret != (unsigned long)CGNS_STRING_SIZE) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }
    config->fields.push_back(str_buf);
  }

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  Restart_Data = new passivedouble[nFields*geometry->GetnPointDomain()];

  /*--- Read in the data for the restart at all local points. ---*/

  ret = fread(Restart_Data, sizeof(passivedouble), nFields*geometry->GetnPointDomain(), fhw);
  if (ret != (unsigned long)nFields*geometry->GetnPointDomain()) {
    cout << endl << "Error reading restart file." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  MPI_Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  unsigned long iPoint_Global, index, iChar;
  string field_buf;

  int rank = MASTER_NODE, ierr;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    if (rank == MASTER_NODE)
      cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
  }

  /*--- First, read the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, read the
   variable string names here. Only the master rank reads the header. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, Restart_Vars, nRestart_Vars, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(Restart_Vars, nRestart_Vars, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    if (rank == MASTER_NODE) {
      cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
  }

  /*--- Store the number of fields to be read for clarity. ---*/

  nFields = Restart_Vars[1];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  char *mpi_str_buf = new char[nFields*CGNS_STRING_SIZE];
  if (rank == MASTER_NODE) {
    disp = nRestart_Vars*sizeof(int);
    MPI_File_read_at(fhw, disp, mpi_str_buf, nFields*CGNS_STRING_SIZE,
                     MPI_CHAR, MPI_STATUS_IGNORE);
  }

  /*--- Broadcast the string names of the variables. ---*/

  SU2_MPI::Bcast(mpi_str_buf, nFields*CGNS_STRING_SIZE, MPI_CHAR,
                 MASTER_NODE, MPI_COMM_WORLD);

  /*--- Now parse the string names and load into the config class in case
   we need them for writing visualization files (SU2_SOL). ---*/

  config->fields.push_back("Point_ID");
  for (iVar = 0; iVar < nFields; iVar++) {
    index = iVar*CGNS_STRING_SIZE;
    field_buf.append("\"");
    for (iChar = 0; iChar < CGNS_STRING_SIZE; iChar++) {
      str_buf[iChar] = mpi_str_buf[index + iChar];
    }
    field_buf.append(str_buf);
    field_buf.append("\"");
    config->fields.push_back(field_buf.c_str());
    field_buf.clear();
  }

  /*--- Free string buffer memory. ---*/

  delete [] mpi_str_buf;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- We need to ignore the 4 ints describing the nVar_Restart and nPoints,
   along with the string names of the variables. ---*/

  disp = nRestart_Vars*sizeof(int) + CGNS_STRING_SIZE*nFields*sizeof(char);

  /*--- Define a derived datatype for this rank's set of non-contiguous data
   that will be placed in the restart. Here, we are collecting each one of the
   points which are distributed throughout the file in blocks of nVar_Restart data. ---*/

  int *blocklen = new int[geometry->GetnPointDomain()];
  int *displace = new int[geometry->GetnPointDomain()];
  int counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    if (geometry->GetGlobal_to_Local_Point(iPoint_Global) > -1) {
      blocklen[counter] = nFields;
      displace[counter] = iPoint_Global*nFields;
      counter++;
    }
  }
  MPI_Type_indexed(geometry->GetnPointDomain(), blocklen, displace, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  Restart_Data = new passivedouble[nFields*geometry->GetnPointDomain()];

  /*--- Collective call for all ranks to read from their view simultaneously. ---*/

  MPI_File_read_all(fhw, Restart_Data, nFields*geometry->GetnPointDomain(), MPI_DOUBLE, &status);

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

  /*--- Free the derived datatype and release temp memory. ---*/

  MPI_Type_free(&filetype);

  delete [] blocklen;
  delete [] displace;
  
#endif
  
}

void CSolver::Read_SU2_Restart_Metadata(CGeometry *geometry, CConfig *config, bool adjoint_run, string val_filename) {

	su2double AoA_ = config->GetAoA();
	su2double AoS_ = config->GetAoS();
	su2double BCThrust_ = config->GetInitial_BCThrust();
	su2double dCD_dCL_ = config->GetdCD_dCL();
 su2double dCMx_dCL_ = config->GetdCMx_dCL();
 su2double dCMy_dCL_ = config->GetdCMy_dCL();
 su2double dCMz_dCL_ = config->GetdCMz_dCL();
	string::size_type position;
	unsigned long ExtIter_ = 0;
	ifstream restart_file;
	bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());

	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	if (config->GetRead_Binary_Restart()) {

		char fname[100];
		strcpy(fname, val_filename.c_str());
		int nVar_Buf = 5;
		int var_buf[5];
		int Restart_Iter = 0;
		passivedouble Restart_Meta_Passive[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
		su2double Restart_Meta[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#ifndef HAVE_MPI

		/*--- Serial binary input. ---*/

		FILE *fhw;
		fhw = fopen(fname,"rb");
    size_t ret;

		/*--- Error check for opening the file. ---*/

		if (!fhw) {
			cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
			exit(EXIT_FAILURE);
		}

		/*--- First, read the number of variables and points. ---*/

		ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Compute (negative) displacements and grab the metadata. ---*/

		fseek(fhw,-(sizeof(int) + 8*sizeof(passivedouble)), SEEK_END);

		/*--- Read the external iteration. ---*/

		ret = fread(&Restart_Iter, sizeof(int), 1, fhw);
    if (ret != 1) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

		/*--- Read the metadata. ---*/

		ret = fread(Restart_Meta_Passive, sizeof(passivedouble), 8, fhw);
    if (ret != 8) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

    for (unsigned short iVar = 0; iVar < 8; iVar++)
      Restart_Meta[iVar] = Restart_Meta_Passive[iVar];

		/*--- Close the file. ---*/

		fclose(fhw);

#else

		/*--- Parallel binary input using MPI I/O. ---*/

		MPI_File fhw;
		MPI_Offset disp;
		int rank = MASTER_NODE, ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		/*--- All ranks open the file using MPI. ---*/

		ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

		/*--- Error check opening the file. ---*/

		if (ierr) {
			if (rank == MASTER_NODE)
				cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Abort(MPI_COMM_WORLD,1);
			MPI_Finalize();
		}

		/*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

		if (rank == MASTER_NODE)
			MPI_File_read(fhw, var_buf, nVar_Buf, MPI_INT, MPI_STATUS_IGNORE);

		/*--- Broadcast the number of variables to all procs and store clearly. ---*/

		SU2_MPI::Bcast(var_buf, nVar_Buf, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      if (rank == MASTER_NODE) {
        cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
        cout << " SU2 reads/writes binary restart files by default." << endl;
        cout << " Note that backward compatibility for ASCII restart files is" << endl;
        cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }

    /*--- Access the metadata. ---*/

		if (rank == MASTER_NODE) {

      /*--- External iteration. ---*/

      disp = (nVar_Buf*sizeof(int) + var_buf[1]*CGNS_STRING_SIZE*sizeof(char) +
              var_buf[1]*var_buf[2]*sizeof(passivedouble));
      MPI_File_read_at(fhw, disp, &Restart_Iter, 1, MPI_INT, MPI_STATUS_IGNORE);

			/*--- Additional doubles for AoA, AoS, etc. ---*/

      disp = (nVar_Buf*sizeof(int) + var_buf[1]*CGNS_STRING_SIZE*sizeof(char) +
              var_buf[1]*var_buf[2]*sizeof(passivedouble) + 1*sizeof(int));
      MPI_File_read_at(fhw, disp, Restart_Meta_Passive, 8, MPI_DOUBLE, MPI_STATUS_IGNORE);

		}

		/*--- Communicate metadata. ---*/

		SU2_MPI::Bcast(&Restart_Iter, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

		/*--- Copy to a su2double structure (because of the SU2_MPI::Bcast
              doesn't work with passive data)---*/

		for (unsigned short iVar = 0; iVar < 8; iVar++)
			Restart_Meta[iVar] = Restart_Meta_Passive[iVar];

		SU2_MPI::Bcast(Restart_Meta, 8, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

		/*--- All ranks close the file after writing. ---*/

		MPI_File_close(&fhw);

#endif

		/*--- Store intermediate vals from file I/O in correct variables. ---*/

		ExtIter_  = Restart_Iter;
		AoA_      = Restart_Meta[0];
		AoS_      = Restart_Meta[1];
		BCThrust_ = Restart_Meta[2];
		dCD_dCL_  = Restart_Meta[3];
  dCMx_dCL_  = Restart_Meta[4];
  dCMy_dCL_  = Restart_Meta[5];
  dCMz_dCL_  = Restart_Meta[6];

	} else {

    /*--- First, check that this is not a binary restart file. ---*/

    char fname[100];
    strcpy(fname, val_filename.c_str());
    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      exit(EXIT_FAILURE);
    }

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      if (rank == MASTER_NODE)
        cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      if (rank == MASTER_NODE) {
        cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
        cout << " SU2 reads/writes binary restart files by default." << endl;
        cout << " Note that backward compatibility for ASCII restart files is" << endl;
        cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }
    
    MPI_File_close(&fhw);
    
#endif

    /*--- Carry on with ASCII metadata reading. ---*/

		restart_file.open(val_filename.data(), ios::in);
		if (restart_file.fail()) {
			if (rank == MASTER_NODE) {
				cout << " Warning: There is no restart file (" << val_filename.data() << ")."<< endl;
				cout << " Computation will continue without updating metadata parameters." << endl;
			}
		} else {

			unsigned long iPoint_Global = 0;
			string text_line;

			/*--- The first line is the header (General description) ---*/

			getline (restart_file, text_line);

			/*--- Space for the solution ---*/

			for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

				getline (restart_file, text_line);

			}

			/*--- Space for extra info (if any) ---*/

			while (getline (restart_file, text_line)) {

				/*--- External iteration ---*/

				position = text_line.find ("EXT_ITER=",0);
				if (position != string::npos) {
					text_line.erase (0,9); ExtIter_ = atoi(text_line.c_str());
				}

				/*--- Angle of attack ---*/

				position = text_line.find ("AOA=",0);
				if (position != string::npos) {
					text_line.erase (0,4); AoA_ = atof(text_line.c_str());
				}

				/*--- Sideslip angle ---*/

				position = text_line.find ("SIDESLIP_ANGLE=",0);
				if (position != string::npos) {
					text_line.erase (0,15); AoS_ = atof(text_line.c_str());
				}

				/*--- BCThrust angle ---*/

				position = text_line.find ("INITIAL_BCTHRUST=",0);
				if (position != string::npos) {
					text_line.erase (0,17); BCThrust_ = atof(text_line.c_str());
				}

				if (adjoint_run) {

					if (config->GetEval_dOF_dCX() == true) {

						/*--- dCD_dCL coefficient ---*/

       position = text_line.find ("DCD_DCL_VALUE=",0);
       if (position != string::npos) {
         text_line.erase (0,14); dCD_dCL_ = atof(text_line.c_str());
       }
       
       /*--- dCMx_dCL coefficient ---*/
       
       position = text_line.find ("DCMX_DCL_VALUE=",0);
       if (position != string::npos) {
         text_line.erase (0,15); dCMx_dCL_ = atof(text_line.c_str());
       }
       
       /*--- dCMy_dCL coefficient ---*/
       
       position = text_line.find ("DCMY_DCL_VALUE=",0);
       if (position != string::npos) {
         text_line.erase (0,15); dCMy_dCL_ = atof(text_line.c_str());
       }
       
       /*--- dCMz_dCL coefficient ---*/
       
       position = text_line.find ("DCMZ_DCL_VALUE=",0);
       if (position != string::npos) {
         text_line.erase (0,15); dCMz_dCL_ = atof(text_line.c_str());
       }
       
					}

				}

			}


			/*--- Close the restart meta file. ---*/

			restart_file.close();

		}
	}

	/*--- Load the metadata. ---*/

	/*--- Only from the direct problem ---*/

	if (!adjoint_run) {

		/*--- Angle of attack ---*/

		if (config->GetDiscard_InFiles() == false) {
			if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE)) {
				cout << fixed <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
				cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
			}
			config->SetAoA(AoA_ + config->GetAoA_Offset());
		}
		else {
			if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE))
				cout <<"WARNING: Discarding the AoA in the solution file." << endl;
		}

		/*--- Sideslip angle ---*/

		if (config->GetDiscard_InFiles() == false) {
			if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE)) {
				cout << fixed <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
				cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
			}
			config->SetAoS(AoS_ + config->GetAoS_Offset());
		}
		else {
			if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE))
				cout <<"WARNING: Discarding the AoS in the solution file." << endl;
		}

		/*--- BCThrust angle ---*/

		if (config->GetDiscard_InFiles() == false) {
			if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
				cout <<"WARNING: SU2 will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
			config->SetInitial_BCThrust(BCThrust_);
		}
		else {
			if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
				cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
		}


		/*--- The adjoint problem needs this information from the direct solution ---*/

		if (adjoint) {

			if (config->GetEval_dOF_dCX() == false) {

				if (config->GetDiscard_InFiles() == false) {

      if ((config->GetdCD_dCL() != dCD_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: SU2 will use the dCD/dCL provided in the direct solution file: " << dCD_dCL_ << "." << endl;
      config->SetdCD_dCL(dCD_dCL_);
      
      if ((config->GetdCMx_dCL() != dCMx_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: SU2 will use the dCMx/dCL provided in the direct solution file: " << dCMx_dCL_ << "." << endl;
      config->SetdCMx_dCL(dCMx_dCL_);
      
      if ((config->GetdCMy_dCL() != dCMy_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: SU2 will use the dCMy/dCL provided in the direct solution file: " << dCMy_dCL_ << "." << endl;
      config->SetdCMy_dCL(dCMy_dCL_);
      
      if ((config->GetdCMz_dCL() != dCMz_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: SU2 will use the dCMz/dCL provided in the direct solution file: " << dCMz_dCL_ << "." << endl;
      config->SetdCMz_dCL(dCMz_dCL_);

				}
				else {
      
      if ((config->GetdCD_dCL() != dCD_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: Discarding the dCD/dCL in the direct solution file." << endl;
      
      if ((config->GetdCMx_dCL() != dCMx_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: Discarding the dCMx/dCL in the direct solution file." << endl;
      
      if ((config->GetdCMy_dCL() != dCMy_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: Discarding the dCMy/dCL in the direct solution file." << endl;
      
      if ((config->GetdCMz_dCL() != dCMz_dCL_) &&  (rank == MASTER_NODE))
        cout <<"WARNING: Discarding the dCMz/dCL in the direct solution file." << endl;
      
    }

			}

		}

	}

	/*--- Only from the adjoint restart file ---*/

	else {

		/*--- The adjoint problem needs this information from the adjoint solution file ---*/

		if (config->GetEval_dOF_dCX() == true) {

			/*--- If it is a restart it will use the value that was stored in the adjoint solution file  ---*/

			if (config->GetRestart()) {

     /*--- dCD_dCL coefficient ---*/
     
     if ((config->GetdCD_dCL() != dCD_dCL_) &&  (rank == MASTER_NODE))
       cout <<"WARNING: SU2 will use the dCD/dCL provided in\nthe adjoint solution file: " << dCD_dCL_ << " ." << endl;
     config->SetdCD_dCL(dCD_dCL_);
     
     /*--- dCMx_dCL coefficient ---*/
     
     if ((config->GetdCMx_dCL() != dCMx_dCL_) &&  (rank == MASTER_NODE))
       cout <<"WARNING: SU2 will use the dCMx/dCL provided in\nthe adjoint solution file: " << dCMx_dCL_ << " ." << endl;
     config->SetdCMx_dCL(dCMx_dCL_);
     
     /*--- dCMy_dCL coefficient ---*/
     
     if ((config->GetdCMy_dCL() != dCMy_dCL_) &&  (rank == MASTER_NODE))
       cout <<"WARNING: SU2 will use the dCMy/dCL provided in\nthe adjoint solution file: " << dCMy_dCL_ << " ." << endl;
     config->SetdCMy_dCL(dCMy_dCL_);
     
     /*--- dCMz_dCL coefficient ---*/
     
     if ((config->GetdCMz_dCL() != dCMz_dCL_) &&  (rank == MASTER_NODE))
       cout <<"WARNING: SU2 will use the dCMz/dCL provided in\nthe adjoint solution file: " << dCMz_dCL_ << " ." << endl;
     config->SetdCMz_dCL(dCMz_dCL_);
     
			}


		}

	}


	/*--- External iteration ---*/

  if ((config->GetDiscard_InFiles() == false) &&
      (!config->GetContinuous_Adjoint() && !config->GetDiscrete_Adjoint()) || (adjoint && config->GetRestart()))
    config->SetExtIter_OffSet(ExtIter_);

}

CBaselineSolver::CBaselineSolver(void) : CSolver() { }

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Routines to access the number of variables and string names. ---*/

  SetOutputVariables(geometry, config);

  /*--- Initialize a zero solution and instantiate the CVariable class. ---*/

  Solution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = 0.0;
  }

  node = new CVariable*[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    node[iPoint] = new CBaselineVariable(Solution, nVar, config);
  }
  
}

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short nVar, vector<string> field_names) {

  unsigned long iPoint;
  unsigned short iVar;

  config->fields = field_names;

  Solution = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar] = 0.0;
  }

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Allocate the node variables ---*/

  node = new CVariable*[geometry->GetnPoint()];

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {

    node[iPoint] = new CBaselineVariable(Solution, nVar, config);

  }

}

void CBaselineSolver::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  /*--- Open the ASCII restart file and extract the nVar and field names. ---*/

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  string Tag, text_line, AdjExt, UnstExt;
  unsigned long iExtIter = config->GetExtIter();
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  ifstream restart_file;
  string filename;

  /*--- Retrieve filename from config ---*/

  if (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint()) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetSolution_FEMFileName();
  } else {
    filename = config->GetSolution_FlowFileName();
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1  || config->GetUnsteady_Simulation() == HARMONIC_BALANCE)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetWrt_Unsteady()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if (config->GetWrt_Dynamic()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  /*--- Read only the number of variables in the restart file. ---*/

  if (config->GetRead_Binary_Restart()) {

    char fname[100];
    strcpy(fname, filename.c_str());
    int nVar_Buf = 5;
    int var_buf[5];

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- First, read the number of variables and points. ---*/

    ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      exit(EXIT_FAILURE);
    }
    
    /*--- Close the file. ---*/

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int rank = MASTER_NODE, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      if (rank == MASTER_NODE)
        cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }

    /*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

    if (rank == MASTER_NODE) {
      MPI_File_read(fhw, var_buf, nVar_Buf, MPI_INT, MPI_STATUS_IGNORE);
    }

    /*--- Broadcast the number of variables to all procs and store more clearly. ---*/

    SU2_MPI::Bcast(var_buf, nVar_Buf, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      if (rank == MASTER_NODE) {
        cout << endl << endl << "Error: file " << fname << " is not a binary SU2 restart file." << endl;
        cout << " SU2 reads/writes binary restart files by default." << endl;
        cout << " Note that backward compatibility for ASCII restart files is" << endl;
        cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }

    /*--- All ranks close the file after writing. ---*/
    
    MPI_File_close(&fhw);

#endif

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];

  } else {

    /*--- First, check that this is not a binary restart file. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      cout << endl << "Error reading restart file." << endl;
      exit(EXIT_FAILURE);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
      cout << " SU2 reads/writes binary restart files by default." << endl;
      cout << " Note that backward compatibility for ASCII restart files is" << endl;
      cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      exit(EXIT_FAILURE);
    }

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      if (rank == MASTER_NODE)
        cout << endl << "Error: unable to open SU2 restart file " << fname << "." << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      if (rank == MASTER_NODE) {
        cout << endl << endl << "Error: file " << fname << " is a binary SU2 restart file, expected ASCII." << endl;
        cout << " SU2 reads/writes binary restart files by default." << endl;
        cout << " Note that backward compatibility for ASCII restart files is" << endl;
        cout << " possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options." << endl << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
    }
    
    MPI_File_close(&fhw);
    
#endif

    /*--- Open the restart file ---*/

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "SU2 solution file " << filename << " not found." << endl;

#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
    }
    
    /*--- Identify the number of fields (and names) in the restart file ---*/

    getline (restart_file, text_line);

    stringstream ss(text_line);
    while (ss >> Tag) {
      config->fields.push_back(Tag);
      if (ss.peek() == ',') ss.ignore();
    }

    /*--- Close the file (the solution date is read later). ---*/
    
    restart_file.close();

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = config->fields.size() - 1;

    /*--- Clear the fields vector since we'll read it again. ---*/

    config->fields.clear();

  }

}

void CBaselineSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR, GridVel_Index;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *transl, *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Solution = NULL;
  
  Solution = new su2double[nVar];

  GridVel_Index = 2*nDim;

  if (config->GetKind_Turb_Model() == SA) { GridVel_Index += 1; }
  else if (config->GetKind_Turb_Model() == SST) { GridVel_Index += 2; }
  if (config->GetKind_Regime() != INCOMPRESSIBLE) { GridVel_Index += 1; }
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
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

          if (config->GetGrid_Movement()) {
            Solution[GridVel_Index + 1] = (rotMatrix[0][0]*Buffer_Receive_U[(GridVel_Index+1)*nVertexR+iVertex] +
                                           rotMatrix[0][1]*Buffer_Receive_U[(GridVel_Index+2)*nVertexR+iVertex]);
            Solution[GridVel_Index + 2] = (rotMatrix[1][0]*Buffer_Receive_U[(GridVel_Index+1)*nVertexR+iVertex] +
                                           rotMatrix[1][1]*Buffer_Receive_U[(GridVel_Index+2)*nVertexR+iVertex]);
          }
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
          
          Solution[nDim+1] = (rotMatrix[0][0]*Buffer_Receive_U[(nDim+1)*nVertexR+iVertex] +
                              rotMatrix[0][1]*Buffer_Receive_U[(nDim+2)*nVertexR+iVertex] +
                              rotMatrix[0][2]*Buffer_Receive_U[(nDim+3)*nVertexR+iVertex]);
          Solution[nDim+2] = (rotMatrix[1][0]*Buffer_Receive_U[(nDim+1)*nVertexR+iVertex] +
                              rotMatrix[1][1]*Buffer_Receive_U[(nDim+2)*nVertexR+iVertex] +
                              rotMatrix[1][2]*Buffer_Receive_U[(nDim+3)*nVertexR+iVertex]);
          Solution[nDim+3] = (rotMatrix[2][0]*Buffer_Receive_U[(nDim+1)*nVertexR+iVertex] +
                              rotMatrix[2][1]*Buffer_Receive_U[(nDim+2)*nVertexR+iVertex] +
                              rotMatrix[2][2]*Buffer_Receive_U[(nDim+3)*nVertexR+iVertex]);

          if (config->GetGrid_Movement()) {
            Solution[GridVel_Index+1] = (rotMatrix[0][0]*Buffer_Receive_U[(GridVel_Index+1)*nVertexR+iVertex] +
                                         rotMatrix[0][1]*Buffer_Receive_U[(GridVel_Index+2)*nVertexR+iVertex] +
                                         rotMatrix[0][2]*Buffer_Receive_U[(GridVel_Index+3)*nVertexR+iVertex]);
            Solution[GridVel_Index+2] = (rotMatrix[1][0]*Buffer_Receive_U[(GridVel_Index+1)*nVertexR+iVertex] +
                                         rotMatrix[1][1]*Buffer_Receive_U[(GridVel_Index+2)*nVertexR+iVertex] +
                                         rotMatrix[1][2]*Buffer_Receive_U[(GridVel_Index+3)*nVertexR+iVertex]);
            Solution[GridVel_Index+3] = (rotMatrix[2][0]*Buffer_Receive_U[(GridVel_Index+1)*nVertexR+iVertex] +
                                         rotMatrix[2][1]*Buffer_Receive_U[(GridVel_Index+2)*nVertexR+iVertex] +
                                         rotMatrix[2][2]*Buffer_Receive_U[(GridVel_Index+3)*nVertexR+iVertex]);
          }
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

void CBaselineSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Restart the solution from file information ---*/

  string filename;
  unsigned long index;
  string UnstExt, text_line, AdjExt;
  ifstream solution_file;
  unsigned short iDim, iVar;
  unsigned long iExtIter = config->GetExtIter();
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool adjoint = ( config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint() ); 
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();
  bool grid_movement  = config->GetGrid_Movement();
  bool steady_restart = config->GetSteadyRestart();
  unsigned short turb_model = config->GetKind_Turb_Model();

  su2double *Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[iZone]->GetnDim();

  /*--- Retrieve filename from config ---*/

  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetSolution_FEMFileName();
  } else {
    filename = config->GetSolution_FlowFileName();
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1 )
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Unsteady problems require an iteration number to be appended. ---*/

  if (config->GetWrt_Unsteady() || config->GetUnsteady_Simulation() != HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if (config->GetWrt_Dynamic()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[iZone], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[iZone], config, filename);
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry[iZone]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[iZone]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      
      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);
     
      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      
      if (grid_movement && val_update_geo) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == SA || turb_model == SA_NEG) {
          index++;
        } else if (turb_model == SST) {
          index+=2;
        }
        
        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/
        
        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[iZone]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[iZone]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
    
  }

  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry[iZone], config);
  
  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (grid_movement && val_update_geo) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[iZone]->Set_MPI_Coord(config);
    geometry[iZone]->Set_MPI_GridVel(config);

  }
  
  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CBaselineSolver::LoadRestart_FSI(CGeometry *geometry, CSolver ***solver, CConfig *config, int val_iter) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Restart the solution from file information ---*/
  string filename;
  unsigned long index;
  string UnstExt, text_line, AdjExt;
  ifstream solution_file;
  unsigned short iVar;
  unsigned long iExtIter = config->GetExtIter();
  bool fem = (config->GetKind_Solver() == FEM_ELASTICITY);
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  /*--- Retrieve filename from config ---*/
  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else if (fem) {
    filename = config->GetSolution_FEMFileName();
  } else {
    filename = config->GetSolution_FlowFileName();
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  /*--- Unsteady problems require an iteration number to be appended. ---*/
  if (config->GetWrt_Unsteady() || config->GetUnsteady_Simulation() != HARMONIC_BALANCE) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  } else if (config->GetWrt_Dynamic()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(iExtIter));
  }

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry, config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry, config, filename);
  }

  unsigned short nVar_Local = Restart_Vars[1];
  su2double *Solution_Local = new su2double[nVar_Local];

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/
  
  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar_Local; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      node[iPoint_Local]->SetSolution(Solution);

      /*--- Increment the overall counter for how many points have been loaded. ---*/

      counter++;

    }

  }

  delete [] Solution_Local;

}

CBaselineSolver::~CBaselineSolver(void) { }
