/*!
 * \file grid_movement_structure.cpp
 * \brief Subroutines for doing the grid movement using different strategies.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.0.0 "eagle"
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

#include "../include/grid_movement_structure.hpp"
#include <list>

using namespace std;

CGridMovement::CGridMovement(void) { }

CGridMovement::~CGridMovement(void) { }

CVolumetricMovement::CVolumetricMovement(CGeometry *geometry) : CGridMovement() {
	
	nDim = geometry->GetnDim();
  
}

CVolumetricMovement::~CVolumetricMovement(void) {

}


void CVolumetricMovement::UpdateGridCoord(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim;
	unsigned long iPoint, total_index;
	double new_coord;
  
  /*--- Update the grid coordinates using the solution of the linear system
   after grid deformation (LinSysSol contains the x, y, z displacements). ---*/
  
	for (iPoint = 0; iPoint < nPoint; iPoint++)
		for (iDim = 0; iDim < nDim; iDim++) {
			total_index = iPoint*nDim + iDim;
			new_coord = geometry->node[iPoint]->GetCoord(iDim)+LinSysSol[total_index];
			if (fabs(new_coord) < EPS*EPS) new_coord = 0.0;
			geometry->node[iPoint]->SetCoord(iDim, new_coord);
		}
  
}

void CVolumetricMovement::UpdateDualGrid(CGeometry *geometry, CConfig *config) {
  
  /*--- After moving all nodes, update the dual mesh. Recompute the edges and
   dual mesh control volumes in the domain and on the boundaries. ---*/
  
	geometry->SetCG();
	geometry->SetControlVolume(config, UPDATE);
	geometry->SetBoundControlVolume(config, UPDATE);
  
}

void CVolumetricMovement::UpdateMultiGrid(CGeometry **geometry, CConfig *config) {
  
  unsigned short iMGfine, iMGlevel, nMGlevel = config->GetMGLevels();
  
  /*--- Update the multigrid structure after moving the finest grid,
   including computing the grid velocities on the coarser levels. ---*/
  
  for (iMGlevel = 1; iMGlevel <= nMGlevel; iMGlevel++) {
    iMGfine = iMGlevel-1;
    geometry[iMGlevel]->SetControlVolume(config,geometry[iMGfine], UPDATE);
    geometry[iMGlevel]->SetBoundControlVolume(config,geometry[iMGfine],UPDATE);
    geometry[iMGlevel]->SetCoord(geometry[iMGfine]);
    if (config->GetGrid_Movement())
      geometry[iMGlevel]->SetRestricted_GridVelocity(geometry[iMGfine],config);
  }
 
}

void CVolumetricMovement::SetVolume_Deformation(CGeometry *geometry, CConfig *config, bool UpdateGeo) {
	unsigned long IterLinSol, iGridDef_Iter;
  double MinVolume, NumError;
  bool Screen_Output = true;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Decide if there is going to be a screen output ---*/
  
  if (config->GetKind_SU2() == SU2_CFD) Screen_Output = false;

  /*--- Initialize the number of spatial dimensions, length of the state
   vector (same as spatial dimensions for grid deformation), and grid nodes. ---*/
  
  nDim   = geometry->GetnDim();
  nVar   = geometry->GetnDim();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize matrix, solution, and r.h.s. structures for the linear solver. ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  StiffMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry);
  
  /*--- Loop over the total number of grid deformation iterations. The surface
   deformation can be divided into increments to help with stability. In
   particular, the linear elasticity equations hold only for small deformations. ---*/
  
  for (iGridDef_Iter = 0; iGridDef_Iter < config->GetGridDef_Iter(); iGridDef_Iter++) {
    
    /*--- Initialize vector and sparse matrix ---*/
    
    LinSysSol.SetValZero();
    LinSysRes.SetValZero();
    StiffMatrix.SetValZero();
    
    /*--- Compute the stiffness matrix entries for all nodes/elements in the
     mesh. FEA uses a finite element method discretization of the linear
     elasticity equations (transfers element stiffnesses to point-to-point). ---*/
    
    MinVolume = SetFEAMethodContributions_Elem(geometry);
    
    /*--- Compute the tolerance of the linear solver using MinLength ---*/
    
    NumError = MinVolume * 0.001;
    
    /*--- Set the boundary displacements (as prescribed by the design variable
     perturbations controlling the surface shape) as a Dirichlet BC. ---*/
    
    SetBoundaryDisplacements(geometry, config);
    
    /*--- Fix the location of any points in the domain, if requested. ---*/
    
    if (config->GetHold_GridFixed())
      SetDomainDisplacements(geometry, config);
    
    /*--- Communicate any prescribed boundary displacements via MPI,
     so that all nodes have the same solution and r.h.s. entries
     across all paritions. ---*/
    
    StiffMatrix.SendReceive_Solution(LinSysSol, geometry, config);
    StiffMatrix.SendReceive_Solution(LinSysRes, geometry, config);
    
    /*--- Definition of the preconditioner matrix vector multiplication, and linear solver ---*/
    
    CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProduct(StiffMatrix, geometry, config);
    CPreconditioner* precond      = new CLU_SGSPreconditioner(StiffMatrix, geometry, config);
    CSysSolve *system             = new CSysSolve();
    
    /*--- Solve the linear system ---*/

    IterLinSol = system->FGMRES(LinSysRes, LinSysSol, *mat_vec, *precond, NumError, 500, Screen_Output);
    
    /*--- Deallocate memory needed by the Krylov linear solver ---*/
    
    delete system;
    delete mat_vec;
    delete precond;
    
    /*--- Update the grid coordinates and cell volumes using the solution
     of the linear system (usol contains the x, y, z displacements). ---*/
    
    UpdateGridCoord(geometry, config);
    if (UpdateGeo)
      UpdateDualGrid(geometry, config);
    
    /*--- Check for failed deformation (negative volumes). ---*/
    
    MinVolume = Check_Grid(geometry);
    
    if ((rank == MASTER_NODE) && Screen_Output) {
      cout << " Non-linear iter.: " << iGridDef_Iter+1 << "/" << config->GetGridDef_Iter()
      << ". Linear iter.: " << IterLinSol << ". Min vol.: " << MinVolume
      << ". Error: " << NumError << "." <<endl;
    }
    
  }
  
  /*--- Deallocate vectors for the linear system. ---*/
  
  LinSysSol.~CSysVector();
  LinSysRes.~CSysVector();
  StiffMatrix.~CSysMatrix();
  
}

double CVolumetricMovement::Check_Grid(CGeometry *geometry) {
  
	unsigned long iElem, ElemCounter = 0, PointCorners[8];
  double Area, Volume, MaxArea = -1E22, MaxVolume = -1E22, MinArea = 1E22, MinVolume = 1E22, CoordCorners[8][3];
  unsigned short nNodes, iNodes, iDim;
  bool RightVol;
  
  int rank = MASTER_NODE;
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
	/*--- Load up each triangle and tetrahedron to check for negative volumes. ---*/
  
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == WEDGE)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    /*--- Triangles ---*/
    
    if (nDim == 2) {
      
      if (nNodes == 3) Area = GetTriangle_Area(CoordCorners);
      if (nNodes == 4) Area = GetRectangle_Area(CoordCorners);
      
      if (Area >= -EPS) RightVol = true;
      else RightVol = false;;
      
      MaxArea = max(MaxArea, Area);
      MinArea = min(MinArea, Area);
      
    }
    
    /*--- Tetrahedra ---*/
    if (nDim == 3) {
      
      if (nNodes == 4) Volume = GetTetra_Volume(CoordCorners);
      if (nNodes == 5) Volume = GetPyram_Volume(CoordCorners);
      if (nNodes == 6) Volume = GetWedge_Volume(CoordCorners);
      if (nNodes == 8) Volume = GetHexa_Volume(CoordCorners);
      
      if (Volume >= -EPS) RightVol = true;
      else RightVol = false;;
      
      MaxVolume = max(MaxVolume, Volume);
      MinVolume = min(MinVolume, Volume);
      
    }
    
    if (!RightVol) ElemCounter++;
    
	}
  
#ifndef NO_MPI
  unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
  double MaxVolume_Local = MaxVolume; MaxVolume = 0.0;
  double MinVolume_Local = MinVolume; MinVolume = 0.0;
#ifdef WINDOWS
  MPI_Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&MinVolume_Local, &MinVolume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&MaxVolume_Local, &MaxVolume, 1, MPI::DOUBLE, MPI::MAX);
  MPI::COMM_WORLD.Allreduce(&MinVolume_Local, &MinVolume, 1, MPI::DOUBLE, MPI::MIN);
#endif
#endif
  
  if ((ElemCounter != 0) && (rank == MASTER_NODE))
    cout <<"There are " << ElemCounter << " elements with negative volume.\n" << endl;
  
  if (nDim == 2) return MinArea;
  else return MinVolume;
  
}

double CVolumetricMovement::SetFEAMethodContributions_Elem(CGeometry *geometry) {
  
	unsigned short iVar, iDim, nNodes, iNodes;
	unsigned long Point_0, Point_1, iElem, iEdge, ElemCounter = 0, PointCorners[8];
  double *Coord_0, *Coord_1, Length, MinLength = 1E10, **StiffMatrix_Elem, Scale, CoordCorners[8][3];
  double *Edge_Vector = new double [nDim];
  bool RightVol;
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Allocate maximum size (rectangle and hexahedron) ---*/
  
  if (nDim == 2) {
    StiffMatrix_Elem = new double* [8];
    for (iVar = 0; iVar < 8; iVar++)
      StiffMatrix_Elem[iVar] = new double [8];
  }
  if (nDim == 3) {
    StiffMatrix_Elem = new double* [24];
    for (iVar = 0; iVar < 24; iVar++)
      StiffMatrix_Elem[iVar] = new double [24];
  }
  
  /*--- Check the minimum edge length in the entire mesh. ---*/
  
	for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
		/*--- Points in edge and coordinates ---*/
		Point_0 = geometry->edge[iEdge]->GetNode(0);  Coord_0 = geometry->node[Point_0]->GetCoord();
		Point_1 = geometry->edge[iEdge]->GetNode(1);  Coord_1 = geometry->node[Point_1]->GetCoord();
    
		/*--- Compute Edge_Vector ---*/
		Length = 0;
		for (iDim = 0; iDim < nDim; iDim++) {
			Edge_Vector[iDim] = Coord_1[iDim] - Coord_0[iDim];
			Length += Edge_Vector[iDim]*Edge_Vector[iDim];
		}
		Length = sqrt(Length);
		MinLength = min(Length, MinLength);
    
	}
  
  /*--- Compute min volume in the entire mesh. ---*/
  Scale = Check_Grid(geometry);
  
	/*--- Compute contributions from each element by forming the stiffness matrix (FEA) ---*/
	for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    
    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)     nNodes = 3;
    if (geometry->elem[iElem]->GetVTK_Type() == RECTANGLE)    nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)  nNodes = 4;
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)      nNodes = 5;
    if (geometry->elem[iElem]->GetVTK_Type() == WEDGE)        nNodes = 6;
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)   nNodes = 8;
    
    for (iNodes = 0; iNodes < nNodes; iNodes++) {
      PointCorners[iNodes] = geometry->elem[iElem]->GetNode(iNodes);
      for (iDim = 0; iDim < nDim; iDim++) {
        CoordCorners[iNodes][iDim] = geometry->node[PointCorners[iNodes]]->GetCoord(iDim);
      }
    }
    
    if (nDim == 2) RightVol = SetFEA_StiffMatrix2D(geometry, StiffMatrix_Elem, CoordCorners, nNodes, Scale);
    if (nDim == 3) RightVol = SetFEA_StiffMatrix3D(geometry, StiffMatrix_Elem, CoordCorners, nNodes, Scale);

    AddFEA_StiffMatrix(geometry, StiffMatrix_Elem, PointCorners, nNodes);
    
    /*--- Create a list with the degenerate elements ---*/

    if (!RightVol) ElemCounter++;
      
	}
  
#ifndef NO_MPI
  unsigned long ElemCounter_Local = ElemCounter; ElemCounter = 0;
#ifdef WINDOWS
  MPI_Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&ElemCounter_Local, &ElemCounter, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
#endif
  
  if ((ElemCounter != 0) && (rank == MASTER_NODE))
    cout <<"There are " << ElemCounter << " degenerate elements in the original grid." << endl;
  
  /*--- Deallocate memory and exit ---*/
  
  if (nDim == 2) {
    for (iVar = 0; iVar < 8; iVar++)
      delete StiffMatrix_Elem[iVar];
    delete [] StiffMatrix_Elem;
  }
  if (nDim == 3) {
    for (iVar = 0; iVar < 24; iVar++)
      delete StiffMatrix_Elem[iVar];
    delete [] StiffMatrix_Elem;
  }
  
  delete [] Edge_Vector;
  
  /*--- If there are no degenerate cells, use the minimum volume instead ---*/
  if (ElemCounter == 0) MinLength = Scale;
  
#ifndef NO_MPI
  double MinLength_Local = MinLength; MinLength = 0.0;
#ifdef WINDOWS
  MPI_Allreduce(&MinLength_Local, &MinLength, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#else
  MPI::COMM_WORLD.Allreduce(&MinLength_Local, &MinLength, 1, MPI::DOUBLE, MPI::MIN);
#endif
#endif
      
	return MinLength;
}

double CVolumetricMovement::ShapeFunc_Hexa(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double a0, a1, a2, c0, c1, c2, xsj;
  double ss[3], xs[3][3], ad[3][3];
  double s0[8] = {-0.5, 0.5, 0.5,-0.5,-0.5, 0.5,0.5,-0.5};
  double s1[8] = {-0.5,-0.5, 0.5, 0.5,-0.5,-0.5,0.5, 0.5};
  double s2[8] = {-0.5,-0.5,-0.5,-0.5, 0.5, 0.5,0.5, 0.5};
  
  ss[0] = Xi;
  ss[1] = Eta;
  ss[2] = Mu;

  /*--- Shape functions ---*/
  
  for (i = 0; i < 8; i++) {
    a0 = 0.5+s0[i]*ss[0]; // shape function in xi-direction
    a1 = 0.5+s1[i]*ss[1]; // shape function in eta-direction
    a2 = 0.5+s2[i]*ss[2]; // shape function in mu-direction
    DShapeFunction[i][0] = s0[i]*a1*a2; // dN/d xi
    DShapeFunction[i][1] = s1[i]*a0*a2; // dN/d eta
    DShapeFunction[i][2] = s2[i]*a0*a1; // dN/d mu
    DShapeFunction[i][3] = a0*a1*a2; // actual shape function N
  }
  
  /*--- Jacobian transformation ---*/
   
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 8; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 8; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CVolumetricMovement::ShapeFunc_Tetra(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = Xi;
  DShapeFunction[1][3] = Eta;
  DShapeFunction[2][3] = Mu;
  DShapeFunction[3][3] = 1.0 - Xi - Eta - Mu;

  /*--- dN/d xi, dN/d eta, dN/d mu ---*/

  DShapeFunction[0][0] = 1.0;   DShapeFunction[0][1] = 0.0;   DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = 0.0;   DShapeFunction[1][1] = 1.0;   DShapeFunction[1][2] = 0.0;
  DShapeFunction[2][0] = 0.0;   DShapeFunction[2][1] = 0.0;   DShapeFunction[2][2] = 1.0;
  DShapeFunction[3][0] = -1.0;  DShapeFunction[3][1] = -1.0;  DShapeFunction[3][2] = -1.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CVolumetricMovement::ShapeFunc_Pyram(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  double Den = 4.0*(1.0 - Mu);
  
  DShapeFunction[0][3] = (-Xi+Eta+Mu-1.0)*(-Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[1][3] = (-Xi-Eta+Mu-1.0)*(Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[2][3] = (Xi+Eta+Mu-1.0)*(Xi-Eta+Mu-1.0)/Den;
  DShapeFunction[3][3] = (Xi+Eta+Mu-1.0)*(-Xi+Eta+Mu-1.0)/Den;
  DShapeFunction[4][3] = Mu;
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = 0.0; DShapeFunction[0][1] = 0.0; DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = 0.0; DShapeFunction[1][1] = 0.0; DShapeFunction[1][2] = 0.0;
  DShapeFunction[2][0] = 0.0; DShapeFunction[2][1] = 0.0; DShapeFunction[2][2] = 0.0;
  DShapeFunction[3][0] = 0.0; DShapeFunction[3][1] = 0.0; DShapeFunction[3][2] = 0.0;
  DShapeFunction[4][0] = 0.0; DShapeFunction[4][1] = 0.0; DShapeFunction[4][2] = 0.0;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 5; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 5; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CVolumetricMovement::ShapeFunc_Wedge(double Xi, double Eta, double Mu, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, c2, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.5*Eta*(1.0-Xi);
  DShapeFunction[1][3] = 0.5*Mu*(1.0-Xi);;
  DShapeFunction[2][3] = 0.5*(1.0-Eta-Mu)*(1.0-Xi);
  DShapeFunction[3][3] = 0.5*Eta*(Xi+1.0);
  DShapeFunction[4][3] = 0.5*Mu*(Xi+1.0);
  DShapeFunction[5][3] = 0.5*(1.0-Eta-Mu)*(Xi+1.0);
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = -0.5*Eta;            DShapeFunction[0][1] = 0.5*(1.0-Xi);      DShapeFunction[0][2] = 0.0;
  DShapeFunction[1][0] = -0.5*Mu;             DShapeFunction[1][1] = 0.0;               DShapeFunction[1][2] = 0.5*(1.0-Xi);
  DShapeFunction[2][0] = -0.5*(1.0-Eta-Mu);   DShapeFunction[2][1] = -0.5*(1.0-Xi);     DShapeFunction[2][2] = -0.5*(1.0-Xi);
  DShapeFunction[3][0] = 0.5*Eta;             DShapeFunction[3][1] = 0.5*(Xi+1.0);      DShapeFunction[3][2] = 0.0;
  DShapeFunction[4][0] = 0.5*Mu;              DShapeFunction[4][1] = 0.0;               DShapeFunction[4][2] = 0.5*(Xi+1.0);
  DShapeFunction[5][0] = 0.5*(1.0-Eta-Mu);    DShapeFunction[5][1] = -0.5*(Xi+1.0);     DShapeFunction[5][2] = -0.5*(Xi+1.0);
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 6; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1]*xs[2][2]-xs[1][2]*xs[2][1];
  ad[0][1] = xs[0][2]*xs[2][1]-xs[0][1]*xs[2][2];
  ad[0][2] = xs[0][1]*xs[1][2]-xs[0][2]*xs[1][1];
  ad[1][0] = xs[1][2]*xs[2][0]-xs[1][0]*xs[2][2];
  ad[1][1] = xs[0][0]*xs[2][2]-xs[0][2]*xs[2][0];
  ad[1][2] = xs[0][2]*xs[1][0]-xs[0][0]*xs[1][2];
  ad[2][0] = xs[1][0]*xs[2][1]-xs[1][1]*xs[2][0];
  ad[2][1] = xs[0][1]*xs[2][0]-xs[0][0]*xs[2][1];
  ad[2][2] = xs[0][0]*xs[1][1]-xs[0][1]*xs[1][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = xs[0][0]*ad[0][0]+xs[0][1]*ad[1][0]+xs[0][2]*ad[2][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 6; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]+xs[0][2]*DShapeFunction[k][2]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]+xs[1][2]*DShapeFunction[k][2]; // dN/dy
    c2 = xs[2][0]*DShapeFunction[k][0]+xs[2][1]*DShapeFunction[k][1]+xs[2][2]*DShapeFunction[k][2]; // dN/dz
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
    DShapeFunction[k][2] = c2; // store dN/dz instead of dN/d mu
  }
  
  return xsj;
  
}

double CVolumetricMovement::ShapeFunc_Triangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 1-Xi-Eta;
  DShapeFunction[1][3] = Xi;
  DShapeFunction[2][3] = Eta;
  
  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = -1.0;  DShapeFunction[0][1] = -1.0;
  DShapeFunction[1][0] = 1;     DShapeFunction[1][1] = 0.0;
  DShapeFunction[2][0] = 0;     DShapeFunction[2][1] = 1;
  
  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 3; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 3; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

double CVolumetricMovement::ShapeFunc_Rectangle(double Xi, double Eta, double CoordCorners[8][3], double DShapeFunction[8][4]) {
  
  int i, j, k;
  double c0, c1, xsj;
  double xs[3][3], ad[3][3];
  
  /*--- Shape functions ---*/
  
  DShapeFunction[0][3] = 0.25*(1.0-Xi)*(1.0-Eta);
  DShapeFunction[1][3] = 0.25*(1.0+Xi)*(1.0-Eta);
  DShapeFunction[2][3] = 0.25*(1.0+Xi)*(1.0+Eta);
  DShapeFunction[3][3] = 0.25*(1.0-Xi)*(1.0+Eta);

  /*--- dN/d xi, dN/d eta, dN/d mu ---*/
  
  DShapeFunction[0][0] = -0.25*(1.0-Eta); DShapeFunction[0][1] = -0.25*(1.0-Xi);
  DShapeFunction[1][0] =  0.25*(1.0-Eta); DShapeFunction[1][1] = -0.25*(1.0+Xi);
  DShapeFunction[2][0] =  0.25*(1.0+Eta); DShapeFunction[2][1] =  0.25*(1.0+Xi);
  DShapeFunction[3][0] = -0.25*(1.0+Eta); DShapeFunction[3][1] =  0.25*(1.0-Xi);

  /*--- Jacobian transformation ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = 0.0;
      for (k = 0; k < 4; k++) {
        xs[i][j] = xs[i][j]+CoordCorners[k][j]*DShapeFunction[k][i];
      }
    }
  }
  
  /*--- Adjoint to jacobian ---*/
  
  ad[0][0] = xs[1][1];
  ad[0][1] = -xs[0][1];
  ad[1][0] = -xs[1][0];
  ad[1][1] = xs[0][0];
  
  /*--- Determinant of jacobian ---*/
  
  xsj = ad[0][0]*ad[1][1]-ad[0][1]*ad[1][0];
  
  /*--- Jacobian inverse ---*/
  
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      xs[i][j] = ad[i][j]/xsj;
    }
  }
  
  /*--- Derivatives with repect to global coordinates ---*/
  
  for (k = 0; k < 4; k++) {
    c0 = xs[0][0]*DShapeFunction[k][0]+xs[0][1]*DShapeFunction[k][1]; // dN/dx
    c1 = xs[1][0]*DShapeFunction[k][0]+xs[1][1]*DShapeFunction[k][1]; // dN/dy
    DShapeFunction[k][0] = c0; // store dN/dx instead of dN/d xi
    DShapeFunction[k][1] = c1; // store dN/dy instead of dN/d eta
  }
  
  return xsj;
  
}

double CVolumetricMovement::GetTriangle_Area(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double a[3], b[3];
  double *Coord_0, *Coord_1, *Coord_2, Area;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

double CVolumetricMovement::GetRectangle_Area(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double a[3], b[3];
  double *Coord_0, *Coord_1, *Coord_2, Area;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area = 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    a[iDim] = Coord_0[iDim]-Coord_2[iDim];
    b[iDim] = Coord_1[iDim]-Coord_2[iDim];
  }
  
  Area += 0.5*fabs(a[0]*b[1]-a[1]*b[0]);
  
  return Area;
  
}

double CVolumetricMovement::GetTetra_Volume(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  double r1[3], r2[3], r3[3], CrossProduct[3], Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;
  
}

double CVolumetricMovement::GetPyram_Volume(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  double r1[3], r2[3], r3[3], CrossProduct[3], Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

double CVolumetricMovement::GetWedge_Volume(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  double r1[3], r2[3], r3[3], CrossProduct[3], Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[5];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[4];
  Coord_2 = CoordCorners[5];
  Coord_3 = CoordCorners[3];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

double CVolumetricMovement::GetHexa_Volume(double CoordCorners[8][3]) {
  
  unsigned short iDim;
  double *Coord_0, *Coord_1, *Coord_2, *Coord_3;
  double r1[3], r2[3], r3[3], CrossProduct[3], Volume;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[1];
  Coord_2 = CoordCorners[2];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume = (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[5];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[2];
  Coord_2 = CoordCorners[3];
  Coord_3 = CoordCorners[7];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[0];
  Coord_1 = CoordCorners[5];
  Coord_2 = CoordCorners[7];
  Coord_3 = CoordCorners[4];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  Coord_0 = CoordCorners[2];
  Coord_1 = CoordCorners[7];
  Coord_2 = CoordCorners[5];
  Coord_3 = CoordCorners[6];
  
  for (iDim = 0; iDim < nDim; iDim++) {
    r1[iDim] = Coord_1[iDim] - Coord_0[iDim];
    r2[iDim] = Coord_2[iDim] - Coord_0[iDim];
    r3[iDim] = Coord_3[iDim] - Coord_0[iDim];
  }
  
	CrossProduct[0] = (r1[1]*r2[2] - r1[2]*r2[1])*r3[0];
	CrossProduct[1] = (r1[2]*r2[0] - r1[0]*r2[2])*r3[1];
	CrossProduct[2] = (r1[0]*r2[1] - r1[1]*r2[0])*r3[2];
  
  Volume += (CrossProduct[0] + CrossProduct[1] + CrossProduct[2])/6.0;
  
  return Volume;

}

bool CVolumetricMovement::SetFEA_StiffMatrix3D(CGeometry *geometry, double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes, double scale) {
  
  double B_Matrix[6][24], D_Matrix[6][6], Aux_Matrix[24][6];
  double Xi = 0.0, Eta = 0.0, Mu = 0.0, Det;
  unsigned short iNode, iVar, jVar, kVar, iGauss, jGauss, kGauss;
  double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  double iWeight, jWeight, kWeight;
  unsigned short nVar = geometry->GetnDim();
  
  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Each element uses their own stiffness which is inversely
   proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
   is a modification to help allow rigid rotation of elements (see
   "Robust Mesh Deformation using the Linear Elasticity Equations" by
   R. P. Dwight. ---*/
  
  for (iGauss = 0; iGauss < 2; iGauss++) {
    for (jGauss = 0; jGauss < 2; jGauss++) {
      for (kGauss = 0; kGauss < 2; kGauss++) {
        
        if (iGauss == 0) { Xi = -0.577350269189626; iWeight = 1.0; }
        if (iGauss == 1) { Xi = 0.577350269189626; iWeight = 1.0; }
        if (jGauss == 0) { Eta = -0.577350269189626; jWeight = 1.0; }
        if (jGauss == 1) { Eta = 0.577350269189626; jWeight = 1.0; }
        if (kGauss == 0) { Mu = -0.577350269189626; kWeight = 1.0; }
        if (kGauss == 1) { Mu = 0.577350269189626; kWeight = 1.0; }
        
        if (nNodes == 4) Det = ShapeFunc_Tetra(Xi, Eta, Mu, CoordCorners, DShapeFunction);
        if (nNodes == 5) Det = ShapeFunc_Pyram(Xi, Eta, Mu, CoordCorners, DShapeFunction);
        if (nNodes == 6) Det = ShapeFunc_Wedge(Xi, Eta, Mu, CoordCorners, DShapeFunction);
        if (nNodes == 8) Det = ShapeFunc_Hexa(Xi, Eta, Mu, CoordCorners, DShapeFunction);
        
        /*--- Compute the B Matrix ---*/
        
        for (iVar = 0; iVar < 6; iVar++)
          for (jVar = 0; jVar < nNodes*nVar; jVar++)
            B_Matrix[iVar][jVar] = 0.0;
        
        for (iNode = 0; iNode < nNodes; iNode++) {
          B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
          B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
          B_Matrix[2][2+iNode*nVar] = DShapeFunction[iNode][2];
          
          B_Matrix[3][0+iNode*nVar] = DShapeFunction[iNode][1];
          B_Matrix[3][1+iNode*nVar] = DShapeFunction[iNode][0];
          
          B_Matrix[4][1+iNode*nVar] = DShapeFunction[iNode][2];
          B_Matrix[4][2+iNode*nVar] = DShapeFunction[iNode][1];
          
          B_Matrix[5][0+iNode*nVar] = DShapeFunction[iNode][2];
          B_Matrix[5][2+iNode*nVar] = DShapeFunction[iNode][0];
        }
        
      double E = scale / (iWeight * jWeight * kWeight * Det) ;
      double Mu = E;
      double Lambda = -E;
        
//        double E = 2E11;
//        double Nu = 0.30;
//        double Mu = E / (2.0*(1.0 + Nu));
//        double Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
        
        /*--- Compute the D Matrix (for plane strain and 3-D)---*/
        
        D_Matrix[0][0] = Lambda + 2.0*Mu;	D_Matrix[0][1] = Lambda;					D_Matrix[0][2] = Lambda;					D_Matrix[0][3] = 0.0;	D_Matrix[0][4] = 0.0;	D_Matrix[0][5] = 0.0;
        D_Matrix[1][0] = Lambda;					D_Matrix[1][1] = Lambda + 2.0*Mu;	D_Matrix[1][2] = Lambda;					D_Matrix[1][3] = 0.0;	D_Matrix[1][4] = 0.0;	D_Matrix[1][5] = 0.0;
        D_Matrix[2][0] = Lambda;					D_Matrix[2][1] = Lambda;					D_Matrix[2][2] = Lambda + 2.0*Mu;	D_Matrix[2][3] = 0.0;	D_Matrix[2][4] = 0.0;	D_Matrix[2][5] = 0.0;
        D_Matrix[3][0] = 0.0;							D_Matrix[3][1] = 0.0;							D_Matrix[3][2] = 0.0;							D_Matrix[3][3] = Mu;	D_Matrix[3][4] = 0.0;	D_Matrix[3][5] = 0.0;
        D_Matrix[4][0] = 0.0;							D_Matrix[4][1] = 0.0;							D_Matrix[4][2] = 0.0;							D_Matrix[4][3] = 0.0;	D_Matrix[4][4] = Mu;	D_Matrix[4][5] = 0.0;
        D_Matrix[5][0] = 0.0;							D_Matrix[5][1] = 0.0;							D_Matrix[5][2] = 0.0;							D_Matrix[5][3] = 0.0;	D_Matrix[5][4] = 0.0;	D_Matrix[5][5] = Mu;
        
        
        /*--- Compute the BT.D Matrix ---*/
        
        for (iVar = 0; iVar < nNodes*nVar; iVar++) {
          for (jVar = 0; jVar < 6; jVar++) {
            Aux_Matrix[iVar][jVar] = 0.0;
            for (kVar = 0; kVar < 6; kVar++)
              Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
          }
        }
        
        /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
         matrix using Gauss integration ---*/
        
        for (iVar = 0; iVar < nNodes*nVar; iVar++) {
          for (jVar = 0; jVar < nNodes*nVar; jVar++) {
            for (kVar = 0; kVar < 6; kVar++) {
              StiffMatrix_Elem[iVar][jVar] += iWeight * jWeight * kWeight * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * Det;
            }
          }
        }
        
      }
    }
    
  }
  
  return true;
  
}

bool CVolumetricMovement::SetFEA_StiffMatrix2D(CGeometry *geometry, double **StiffMatrix_Elem, double CoordCorners[8][3], unsigned short nNodes, double scale) {
  
  double B_Matrix[3][8], D_Matrix[3][3], Aux_Matrix[8][3];
  double Xi = 0.0, Eta = 0.0, Det;
  unsigned short iNode, iVar, jVar, kVar, iGauss, jGauss;
  double DShapeFunction[8][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
  double iWeight, jWeight;
  unsigned short nVar = geometry->GetnDim();
  
  for (iVar = 0; iVar < nNodes*nVar; iVar++) {
    for (jVar = 0; jVar < nNodes*nVar; jVar++) {
      StiffMatrix_Elem[iVar][jVar] = 0.0;
    }
  }
  
  /*--- Each element uses their own stiffness which is inversely
   proportional to the area/volume of the cell. Using Mu = E & Lambda = -E
   is a modification to help allow rigid rotation of elements (see
   "Robust Mesh Deformation using the Linear Elasticity Equations" by
   R. P. Dwight. ---*/
  
  for (iGauss = 0; iGauss < 2; iGauss++) {
    for (jGauss = 0; jGauss < 2; jGauss++) {
      
      if (iGauss == 0) { Xi = -0.577350269189626; iWeight = 1.0; }
      if (iGauss == 1) { Xi = 0.577350269189626; iWeight = 1.0; }
      if (jGauss == 0) { Eta = -0.577350269189626; jWeight = 1.0; }
      if (jGauss == 1) { Eta = 0.577350269189626; jWeight = 1.0; }
      
      if (nNodes == 3) Det = ShapeFunc_Triangle(Xi, Eta, CoordCorners, DShapeFunction);
      if (nNodes == 4) Det = ShapeFunc_Rectangle(Xi, Eta, CoordCorners, DShapeFunction);
      
      /*--- Compute the B Matrix ---*/
      
      for (iVar = 0; iVar < 3; iVar++)
        for (jVar = 0; jVar < nNodes*nVar; jVar++)
          B_Matrix[iVar][jVar] = 0.0;
      
      for (iNode = 0; iNode < nNodes; iNode++) {
        B_Matrix[0][0+iNode*nVar] = DShapeFunction[iNode][0];
        B_Matrix[1][1+iNode*nVar] = DShapeFunction[iNode][1];
        
        B_Matrix[2][0+iNode*nVar] = DShapeFunction[iNode][1];
        B_Matrix[2][1+iNode*nVar] = DShapeFunction[iNode][0];
      }
      
      double E = scale / (iWeight * jWeight * Det) ;
      double Mu = E;
      double Lambda = -E;
      
//      double E = 2E11;
//      double Nu = 0.30;
//      double Mu = E / (2.0*(1.0 + Nu));
//      double Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
      
      /*--- Compute the D Matrix (for plane strain and 3-D)---*/
      
      D_Matrix[0][0] = Lambda + 2.0*Mu;		D_Matrix[0][1] = Lambda;            D_Matrix[0][2] = 0.0;
      D_Matrix[1][0] = Lambda;            D_Matrix[1][1] = Lambda + 2.0*Mu;   D_Matrix[1][2] = 0.0;
      D_Matrix[2][0] = 0.0;               D_Matrix[2][1] = 0.0;               D_Matrix[2][2] = Mu;
      
      
      /*--- Compute the BT.D Matrix ---*/
      
      for (iVar = 0; iVar < nNodes*nVar; iVar++) {
        for (jVar = 0; jVar < 3; jVar++) {
          Aux_Matrix[iVar][jVar] = 0.0;
          for (kVar = 0; kVar < 3; kVar++)
            Aux_Matrix[iVar][jVar] += B_Matrix[kVar][iVar]*D_Matrix[kVar][jVar];
        }
      }
      
      /*--- Compute the BT.D.B Matrix (stiffness matrix), and add to the original
       matrix using Gauss integration ---*/
      
      for (iVar = 0; iVar < nNodes*nVar; iVar++) {
        for (jVar = 0; jVar < nNodes*nVar; jVar++) {
          for (kVar = 0; kVar < 3; kVar++) {
            StiffMatrix_Elem[iVar][jVar] += iWeight * jWeight * Aux_Matrix[iVar][kVar]*B_Matrix[kVar][jVar] * Det;
          }
        }
      }
      
    }
  }
  
  return true;
  
}

void CVolumetricMovement::AddFEA_StiffMatrix(CGeometry *geometry, double **StiffMatrix_Elem, unsigned long PointCorners[8], unsigned short nNodes) {
  unsigned short iVar, jVar, iDim, jDim;
  unsigned short nVar = geometry->GetnDim();

  double **StiffMatrix_Node;
  StiffMatrix_Node = new double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    StiffMatrix_Node[iVar] = new double [nVar];
  
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      StiffMatrix_Node[iVar][jVar] = 0.0;
  
  /*--- Transform the stiffness matrix for the hexahedral element into the
   contributions for the individual nodes relative to each other. ---*/
  
  for (iVar = 0; iVar < nNodes; iVar++) {
    for (jVar = 0; jVar < nNodes; jVar++) {
      
      for (iDim = 0; iDim < nVar; iDim++) {
        for (jDim = 0; jDim < nVar; jDim++) {
          StiffMatrix_Node[iDim][jDim] = StiffMatrix_Elem[(iVar*nVar)+iDim][(jVar*nVar)+jDim];
        }
      }

      StiffMatrix.AddBlock(PointCorners[iVar], PointCorners[jVar], StiffMatrix_Node);
      
    }
  }
  
  /*--- Deallocate memory and exit ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete StiffMatrix_Node[iVar];
  delete [] StiffMatrix_Node;
  
}

void CVolumetricMovement::SetBoundaryDisplacements(CGeometry *geometry, CConfig *config) {

	unsigned short iDim, nDim = geometry->GetnDim(), iMarker, axis = 0;
	unsigned long iPoint, total_index, iVertex;
	double *VarCoord, MeanCoord[3], VarIncrement = 1.0;
  
  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_MDC will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/
  
  unsigned short Kind_SU2 = config->GetKind_SU2();
  
  /*--- If requested (no by default) impose the surface deflections in
   increments and solve the grid deformation equations iteratively with
   successive small deformations. ---*/
  
  VarIncrement = 1.0/((double)config->GetGridDef_Iter());
	
	/*--- As initialization, set to zero displacements of all the surfaces except the symmetry
	 plane and the receive boundaries. ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if ((config->GetMarker_All_Boundary(iMarker) != SYMMETRY_PLANE)
        && (config->GetMarker_All_Boundary(iMarker) != SEND_RECEIVE)) {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					LinSysRes[total_index] = 0.0;
					LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
				}
			}
    }
  }
	
  /*--- Set to zero displacements of the normal component for the symmetry plane condition ---*/
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if ((config->GetMarker_All_Boundary(iMarker) == SYMMETRY_PLANE) && (nDim == 3)) {
      
			for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = 0.0;
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord = geometry->node[iPoint]->GetCoord();
				for (iDim = 0; iDim < nDim; iDim++)
					MeanCoord[iDim] += VarCoord[iDim]*VarCoord[iDim];
			}
			for (iDim = 0; iDim < nDim; iDim++) MeanCoord[iDim] = sqrt(MeanCoord[iDim]);
			
			if ((MeanCoord[0] <= MeanCoord[1]) && (MeanCoord[0] <= MeanCoord[2])) axis = 0;
			if ((MeanCoord[1] <= MeanCoord[0]) && (MeanCoord[1] <= MeanCoord[2])) axis = 1;
			if ((MeanCoord[2] <= MeanCoord[0]) && (MeanCoord[2] <= MeanCoord[1])) axis = 2;
						
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				total_index = iPoint*nDim + axis;
				LinSysRes[total_index] = 0.0;
				LinSysSol[total_index] = 0.0;
				StiffMatrix.DeleteValsRowi(total_index);
			}
		}
	}
  
	/*--- Set the known displacements, note that some points of the moving surfaces
   could be on on the symmetry plane, we should specify DeleteValsRowi again (just in case) ---*/
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_MDC))) {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					LinSysRes[total_index] = VarCoord[iDim] * VarIncrement;
					LinSysSol[total_index] = VarCoord[iDim] * VarIncrement;
          StiffMatrix.DeleteValsRowi(total_index);
				}
			}
    }
  }
  
  /*--- Don't move the nearfield plane ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		if (config->GetMarker_All_Boundary(iMarker) == NEARFIELD_BOUNDARY) {
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				for (iDim = 0; iDim < nDim; iDim++) {
					total_index = iPoint*nDim + iDim;
					LinSysRes[total_index] = 0.0;
					LinSysSol[total_index] = 0.0;
          StiffMatrix.DeleteValsRowi(total_index);
				}
			}
    }
  }

}

void CVolumetricMovement::SetDomainDisplacements(CGeometry *geometry, CConfig *config) {
	unsigned short iDim, nDim = geometry->GetnDim();
	unsigned long iPoint, total_index;
	double *Coord, MinCoordValues[3], MaxCoordValues[3], *Hold_GridFixed_Coord;
	
	Hold_GridFixed_Coord = config->GetHold_GridFixed_Coord();
	
	MinCoordValues[0] = Hold_GridFixed_Coord[0];
	MinCoordValues[1] = Hold_GridFixed_Coord[1];
	MinCoordValues[2] = Hold_GridFixed_Coord[2];
	MaxCoordValues[0] = Hold_GridFixed_Coord[3];
	MaxCoordValues[1] = Hold_GridFixed_Coord[4];
	MaxCoordValues[2] = Hold_GridFixed_Coord[5];

	/*--- Set to zero displacements of all the points that are not going to be moved
	 except the surfaces ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
		Coord = geometry->node[iPoint]->GetCoord();
		for (iDim = 0; iDim < nDim; iDim++) {
			if ((Coord[iDim] < MinCoordValues[iDim]) || (Coord[iDim] > MaxCoordValues[iDim])) {
				total_index = iPoint*nDim + iDim;
				LinSysRes[total_index] = 0.0;
				LinSysSol[total_index] = 0.0;
				StiffMatrix.DeleteValsRowi(total_index);
			}
		}
	}
}

void CVolumetricMovement::Rigid_Rotation(CGeometry *geometry, CConfig *config,
                                         unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
	/*--- Local variables ---*/
	unsigned short iDim, nDim; 
	unsigned long iPoint;
	double r[3], rotCoord[3], *Coord, Center[3], Omega[3], Lref, dt;
  double *GridVel, newGridVel[3];
	double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
	double dtheta, dphi, dpsi, cosTheta, sinTheta;
	double cosPhi, sinPhi, cosPsi, sinPsi;
	bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
	bool adjoint = config->GetAdjoint();

	/*--- Problem dimension and physical time step ---*/
	nDim = geometry->GetnDim();
	dt   = config->GetDelta_UnstTimeND();
	Lref = config->GetLength_Ref();

  /*--- For time-spectral, motion is the same in each zone (at each instance).
   *    This is used for calls to the config container ---*/
  if (time_spectral)
	  iZone = ZONE_0;
  
  /*--- For the unsteady adjoint, use reverse time ---*/
  if (adjoint) {
    /*--- Set the first adjoint mesh position to the final direct one ---*/
    if (iter == 0) dt = ((double)config->GetnExtIter()-1)*dt;
    /*--- Reverse the rotation direction for the adjoint ---*/
    else dt = -1.0*dt;
  } else {
    /*--- No rotation at all for the first direct solution ---*/
    if (iter == 0) dt = 0;
  }
  
  /*--- Center of rotation & angular velocity vector from config ---*/
  
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = (config->GetRotation_Rate_X(iZone)/config->GetOmega_Ref());
  Omega[1]  = (config->GetRotation_Rate_Y(iZone)/config->GetOmega_Ref());
  Omega[2]  = (config->GetRotation_Rate_Z(iZone)/config->GetOmega_Ref());

  /*-- Set dt for time-spectral cases ---*/
  if (time_spectral) {
	  /*--- period of oscillation & compute time interval using nTimeInstances ---*/
	  double period = config->GetTimeSpectral_Period();
	  dt = period * (double)iter/(double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta change in the angle about the x, y, & z axes. ---*/

  dtheta = Omega[0]*dt;
  dphi   = Omega[1]*dt;
  dpsi   = Omega[2]*dt;

  if (rank == MASTER_NODE && iter == 0) {
    cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
  }
  
	/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
	cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
	sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
	/*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
	rotMatrix[0][0] = cosPhi*cosPsi;
	rotMatrix[1][0] = cosPhi*sinPsi;
	rotMatrix[2][0] = -sinPhi;
  
	rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
	rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
	rotMatrix[2][1] = sinTheta*cosPhi;
  
	rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
	rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
	rotMatrix[2][2] = cosTheta*cosPhi;
  
	/*--- Loop over and rotate each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + Omega[1]*rotCoord[2] - Omega[2]*rotCoord[1];
    newGridVel[1] = GridVel[1] + Omega[2]*rotCoord[0] - Omega[0]*rotCoord[2];
    newGridVel[2] = GridVel[2] + Omega[0]*rotCoord[1] - Omega[1]*rotCoord[0];
    
    /*--- Store new node location & grid velocity. Add center. 
     Do not store the grid velocity if this is an adjoint calculation.---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,rotCoord[iDim] + Center[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
      
    }
  }
  
	/*--- After moving all nodes, update geometry class ---*/
  
	UpdateDualGrid(geometry, config);

}

void CVolumetricMovement::Rigid_Pitching(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Local variables ---*/
  double r[3], rotCoord[3],*Coord, Center[3], Omega[3], Ampl[3], Phase[3];
  double Lref, deltaT, alphaDot[3], *GridVel, newGridVel[3];
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
  double time_new, time_old;
  double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iDim;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool adjoint = config->GetAdjoint();
  
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND(); 
  Lref   = config->GetLength_Ref();

  /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
  if (time_spectral) {
	  iZone = ZONE_0;
  }

  /*--- Pitching origin, frequency, and amplitude from config. ---*/	
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = (config->GetPitching_Omega_X(iZone)/config->GetOmega_Ref());
  Omega[1]  = (config->GetPitching_Omega_Y(iZone)/config->GetOmega_Ref());
  Omega[2]  = (config->GetPitching_Omega_Z(iZone)/config->GetOmega_Ref());
  Ampl[0]   = config->GetPitching_Ampl_X(iZone)*DEG2RAD;
  Ampl[1]   = config->GetPitching_Ampl_Y(iZone)*DEG2RAD;
  Ampl[2]   = config->GetPitching_Ampl_Z(iZone)*DEG2RAD;
  Phase[0]   = config->GetPitching_Phase_X(iZone)*DEG2RAD;
  Phase[1]   = config->GetPitching_Phase_Y(iZone)*DEG2RAD;
  Phase[2]   = config->GetPitching_Phase_Z(iZone)*DEG2RAD;

  if (time_spectral) {    
	  /*--- period of oscillation & compute time interval using nTimeInstances ---*/
	  double period = config->GetTimeSpectral_Period();
	  deltaT = period/(double)(config->GetnTimeInstances());
  }

  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    if (time_spectral) {
    	/*--- For time-spectral, begin movement from the zero position ---*/
    	time_old = 0.0;
    } else {
    	time_old = time_new;
    	if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
    }
  }
  
	/*--- Compute delta change in the angle about the x, y, & z axes. ---*/
  
	dtheta = -Ampl[0]*(sin(Omega[0]*time_new + Phase[0]) - sin(Omega[0]*time_old + Phase[0]));
	dphi   = -Ampl[1]*(sin(Omega[1]*time_new + Phase[1]) - sin(Omega[1]*time_old + Phase[1]));
	dpsi   = -Ampl[2]*(sin(Omega[2]*time_new + Phase[2]) - sin(Omega[2]*time_old + Phase[2]));
  
  /*--- Angular velocity at the new time ---*/
  
  alphaDot[0] = -Omega[0]*Ampl[0]*cos(Omega[0]*time_new);
  alphaDot[1] = -Omega[1]*Ampl[1]*cos(Omega[1]*time_new);
  alphaDot[2] = -Omega[2]*Ampl[2]*cos(Omega[2]*time_new);

  if (rank == MASTER_NODE && iter == 0) {
      cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
      cout << ", " << Omega[2] << ") rad/s." << endl;
      cout << " Pitching amplitude: (" << Ampl[0]/DEG2RAD << ", ";
      cout << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
      cout << ") degrees."<< endl;
      cout << " Pitching phase lag: (" << Phase[0]/DEG2RAD << ", ";
      cout << Phase[1]/DEG2RAD <<", "<< Phase[2]/DEG2RAD;
      cout << ") degrees."<< endl;
  }
  
	/*--- Store angles separately for clarity. Compute sines/cosines. ---*/
  
	cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
	sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
  
	/*--- Compute the rotation matrix. Note that the implicit
   ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
  
	rotMatrix[0][0] = cosPhi*cosPsi;
	rotMatrix[1][0] = cosPhi*sinPsi;
	rotMatrix[2][0] = -sinPhi;
  
	rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
	rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
	rotMatrix[2][1] = sinTheta*cosPhi;
  
	rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
	rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
	rotMatrix[2][2] = cosTheta*cosPhi;
  
	/*--- Loop over and rotate each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    if (nDim == 2) r[nDim] = 0.0;
    
    /*--- Compute transformed point coordinates ---*/
    rotCoord[0] = rotMatrix[0][0]*r[0] 
                + rotMatrix[0][1]*r[1] 
                + rotMatrix[0][2]*r[2];
    
    rotCoord[1] = rotMatrix[1][0]*r[0] 
                + rotMatrix[1][1]*r[1] 
                + rotMatrix[1][2]*r[2];
    
    rotCoord[2] = rotMatrix[2][0]*r[0] 
                + rotMatrix[2][1]*r[1] 
                + rotMatrix[2][2]*r[2];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to 
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + alphaDot[1]*rotCoord[2] - alphaDot[2]*rotCoord[1];
    newGridVel[1] = GridVel[1] + alphaDot[2]*rotCoord[0] - alphaDot[0]*rotCoord[2];
    newGridVel[2] = GridVel[2] + alphaDot[0]*rotCoord[1] - alphaDot[1]*rotCoord[0];
    
    /*--- Store new node location & grid velocity. Add center location.
     Do not store the grid velocity if this is an adjoint calculation.---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,rotCoord[iDim]+Center[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim,newGridVel[iDim]);
    }
  }
  
	/*--- After moving all nodes, update geometry class ---*/
  
	UpdateDualGrid(geometry, config);
  
}

void CVolumetricMovement::Rigid_Plunging(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Local variables ---*/
  double deltaX[3], newCoord[3], Center[3], *Coord, Omega[3], Ampl[3], Lref;
  double *GridVel, newGridVel[3], xDot[3];
  double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool adjoint = config->GetAdjoint();
  
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
  if (time_spectral) {
	  iZone = ZONE_0;
  }
  
  /*--- Plunging frequency and amplitude from config. ---*/
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = (config->GetPlunging_Omega_X(iZone)/config->GetOmega_Ref());
  Omega[1]  = (config->GetPlunging_Omega_Y(iZone)/config->GetOmega_Ref());
  Omega[2]  = (config->GetPlunging_Omega_Z(iZone)/config->GetOmega_Ref());
  Ampl[0]   = config->GetPlunging_Ampl_X(iZone)/Lref;
  Ampl[1]   = config->GetPlunging_Ampl_Y(iZone)/Lref;
  Ampl[2]   = config->GetPlunging_Ampl_Z(iZone)/Lref;
  
  if (time_spectral) {
	  /*--- period of oscillation & time interval using nTimeInstances ---*/
	  double period = config->GetTimeSpectral_Period();
	  deltaT = period/(double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    if (time_spectral) {
    	/*--- For time-spectral, begin movement from the zero position ---*/
    	time_old = 0.0;
    } else {
    	time_old = time_new;
    	if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
    }
  }
  
	/*--- Compute delta change in the position in the x, y, & z directions. ---*/
	deltaX[0] = -Ampl[0]*(sin(Omega[0]*time_new) - sin(Omega[0]*time_old));
	deltaX[1] = -Ampl[1]*(sin(Omega[1]*time_new) - sin(Omega[1]*time_old));
	deltaX[2] = -Ampl[2]*(sin(Omega[2]*time_new) - sin(Omega[2]*time_old));
  
  /*--- Compute grid velocity due to plunge in the x, y, & z directions. ---*/
	xDot[0] = -Ampl[0]*Omega[0]*(cos(Omega[0]*time_new));
	xDot[1] = -Ampl[1]*Omega[1]*(cos(Omega[1]*time_new));
	xDot[2] = -Ampl[2]*Omega[2]*(cos(Omega[2]*time_new));
  
  if (rank == MASTER_NODE && iter == 0) {
    cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
    cout << ", " << Omega[2] << ") rad/s." << endl;
    cout << " Plunging amplitude: (" << Ampl[0] << ", ";
    cout << Ampl[1] << ", " << Ampl[2] <<  ") m."<< endl;
  }
  
	/*--- Loop over and move each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord   = geometry->node[iPoint]->GetCoord();
    GridVel = geometry->node[iPoint]->GetGridVel();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Cross Product of angular velocity and distance from center.
     Note that we have assumed the grid velocities have been set to
     an initial value in the plunging routine. ---*/
    
    newGridVel[0] = GridVel[0] + xDot[0];
    newGridVel[1] = GridVel[1] + xDot[1];
    newGridVel[2] = GridVel[2] + xDot[2];
    
    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,newCoord[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim,newGridVel[iDim]);
    }
  }
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/
  
  config->SetMotion_Origin_X(iZone,Center[0]+deltaX[0]);
  config->SetMotion_Origin_Y(iZone,Center[1]+deltaX[1]);
  config->SetMotion_Origin_Z(iZone,Center[2]+deltaX[2]);
  
  /*--- As the body origin may have moved, pring it to the console ---*/
  
  if (rank == MASTER_NODE) {
    cout << " Body origin: (" << Center[0]+deltaX[0];
    cout << ", " << Center[1]+deltaX[1] << ", " << Center[2]+deltaX[2];
    cout << ")." << endl;
  }
  
	/*--- After moving all nodes, update geometry class ---*/
	
  UpdateDualGrid(geometry, config);
  
}

void CVolumetricMovement::Rigid_Translation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Local variables ---*/
  double deltaX[3], newCoord[3], Center[3], *Coord, Lref;
  double xDot[3];
  double deltaT, time_new, time_old;
  unsigned short iDim, nDim = geometry->GetnDim();
  unsigned long iPoint;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool adjoint = config->GetAdjoint();
	
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- For time-spectral, motion is the same in each zone (at each instance). ---*/
  if (time_spectral) {
	  iZone = ZONE_0;
  }

  /*--- Get motion center and translation rates from config ---*/
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  xDot[0]   = config->GetTranslation_Rate_X(iZone);
  xDot[1]   = config->GetTranslation_Rate_Y(iZone);
  xDot[2]   = config->GetTranslation_Rate_Z(iZone);
  
  if (time_spectral) {
	  /*--- period of oscillation & time interval using nTimeInstances ---*/
	  double period = config->GetTimeSpectral_Period();
	  deltaT = period/(double)(config->GetnTimeInstances());
  }
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    if (time_spectral) {
    	/*--- For time-spectral, begin movement from the zero position ---*/
    	time_old = 0.0;
    } else {
    	time_old = time_new;
    	if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
    }
  }
  
	/*--- Compute delta change in the position in the x, y, & z directions. ---*/
	deltaX[0] = xDot[0]*deltaT;
	deltaX[1] = xDot[1]*deltaT;
	deltaX[2] = xDot[2]*deltaT;

  if (rank == MASTER_NODE) {
    cout << " New physical time: " << time_new << " seconds." << endl;
    if (iter == 0) {
    cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
    cout << ", " << xDot[2] << ") m/s." << endl;
    }
  }
  
	/*--- Loop over and move each node in the volume mesh ---*/
	for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Increment the node position using the delta values. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      newCoord[iDim] = Coord[iDim] + deltaX[iDim];
    
    /*--- Store new node location & grid velocity. Do not store the grid
     velocity if this is an adjoint calculation. ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      geometry->node[iPoint]->SetCoord(iDim,newCoord[iDim]);
      if (!adjoint) geometry->node[iPoint]->SetGridVel(iDim,xDot[iDim]);
    }
  }
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rigid translation. This
   new location will be used for subsequent pitching/rotation.---*/
  
  config->SetMotion_Origin_X(iZone,Center[0]+deltaX[0]);
  config->SetMotion_Origin_Y(iZone,Center[1]+deltaX[1]);
  config->SetMotion_Origin_Z(iZone,Center[2]+deltaX[2]);
  
	/*--- After moving all nodes, update geometry class ---*/
	
  UpdateDualGrid(geometry, config);
  
}

CSurfaceMovement::CSurfaceMovement(void) : CGridMovement() {
	nFFDBox = 0;
	FFDBoxDefinition = false;
}

CSurfaceMovement::~CSurfaceMovement(void) {}

void CSurfaceMovement::SetSurface_Deformation(CGeometry *geometry, CConfig *config) {
  unsigned short iFFDBox, iDV, iLevel, iChild, iParent, jFFDBox;
	char buffer_char[50];
	int rank = MASTER_NODE, iExtIter = 0;
	string FFDBoxTag;
  
  unsigned short nDim = geometry->GetnDim();
	
  /*--- Definition of the FFD deformation class ---*/
	FFDBox = new CFreeFormDefBox*[MAX_NUMBER_FFD];
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Arbitrary definition of surface coordinates from file. ---*/
  if (config->GetDesign_Variable(0) == SURFACE_FILE) {
    
    /*--- Check whether a surface file exists for input ---*/
    ofstream Surface_File;
    string filename = config->GetMotion_FileName();
    Surface_File.open(filename.c_str(), ios::in);
    
    /*--- A surface file does not exist, so write a new one for the
     markers that are specified as part of the motion. ---*/
    if (Surface_File.fail()) {
      
      if (rank == MASTER_NODE)
        cout << "No surface file found. Writing a new file: " << filename << "." << endl;
      
      Surface_File.open(filename.c_str(), ios::out);
      Surface_File.precision(15);
      unsigned long iMarker, jPoint, GlobalIndex, iVertex; double *Coords;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_DV(iMarker) == YES) {
          for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
            Coords = geometry->node[jPoint]->GetCoord();
            Surface_File << GlobalIndex << "\t" << Coords[0] << "\t" << Coords[1];
            if (nDim == 2) Surface_File << endl;
            else Surface_File << "\t" << Coords[2] << endl;
          }
        }
      }
      Surface_File.close();
      
      /*--- A surface file exists, so read in the coordinates ---*/
    } else {
      Surface_File.close();
      if (rank == MASTER_NODE)
        cout << "Updating the surface coordinates from the input file." << endl;
      SetExternal_Deformation(geometry, config, ZONE_0, iExtIter);
    }
    
    /*--- Spherical parameterization ---*/
  } else if (config->GetDesign_Variable(0) == SPHERICAL)  {
    if (rank == MASTER_NODE) cout << "Perform 3D deformation of the surface." << endl;
    SetSpherical(geometry, config, 0, false); // Note that the loop over the design variables is inside the subroutine
    
    /*--- Bump deformation for 2D problems ---*/
  } else if (nDim == 2) {
    
    /*--- Apply rotation, displacement and stretching design variables (this
     should be done before the bump function design variables) ---*/
    for (iDV = 0; iDV < config->GetnDV(); iDV++) {
			switch ( config->GetDesign_Variable(iDV) ) {
				case ROTATION :     SetRotation(geometry, config, iDV, false); break;
				case DISPLACEMENT : SetDisplacement(geometry, config, iDV, false); break;
			}
		}
    
		/*--- Apply the design variables to the control point position ---*/
		for (iDV = 0; iDV < config->GetnDV(); iDV++) {
			switch ( config->GetDesign_Variable(iDV) ) {
				case HICKS_HENNE :  SetHicksHenne(geometry, config, iDV, false); break;
				case COSINE_BUMP :  SetCosBump(geometry, config, iDV, false); break;
				case FOURIER :      SetFourier(geometry, config, iDV, false); break;
				case NACA_4DIGITS : SetNACA_4Digits(geometry, config); break;
				case PARABOLIC :    SetParabolic(geometry, config); break;
				case OBSTACLE :     SetObstacle(geometry, config); break;
				case AIRFOIL :      SetAirfoil(geometry, config); break;
        case SURFACE_FILE : SetExternal_Deformation(geometry, config, ZONE_0, iExtIter); break;
			}
		}
        
    /*--- Free Form Deformation for 3D problems ---*/
	} else if (nDim == 3) {
		    
    /*--- Read the FFD information fron the grid file ---*/
    ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName(), true);
    
    /*--- If the FFDBox was not defined in the input file ---*/
    if (!GetFFDBoxDefinition()) {
      
      if ((rank == MASTER_NODE) && (GetnFFDBox() != 0))
        cout << endl <<"----------------- FFD technique (cartesian -> parametric) ---------------" << endl;
      
      /*--- Create a unitary FFDBox as baseline for other FFDBoxs shapes ---*/
      CFreeFormDefBox FFDBox_unitary(1,1,1);
      FFDBox_unitary.SetUnitCornerPoints();
      
      /*--- Compute the control points of the unitary box, in this case the degree is 1 and the order is 2 ---*/
      FFDBox_unitary.SetControlPoints_Parallelepiped();
      
      for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
        /*--- Compute the support control points for the final FFD using the unitary box ---*/
        FFDBox_unitary.SetSupportCP(FFDBox[iFFDBox]);
        
        /*--- Compute control points in the support box ---*/
        FFDBox_unitary.SetSupportCPChange(FFDBox[iFFDBox]);
        
        /*--- Compute the parametric coordinates, it also find the points in
         the FFDBox using the parametrics coordinates ---*/
        SetParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
        
      }
      
    }
    
    /*--- Output original FFD FFDBox ---*/
    if (rank == MASTER_NODE) {
      for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
        sprintf (buffer_char, "original_FFDBox.plt");
        if (iFFDBox == 0) FFDBox[iFFDBox]->SetTecplot(buffer_char, true);
        else FFDBox[iFFDBox]->SetTecplot(buffer_char, false);
      }
    }
    
    if ((rank == MASTER_NODE) && (GetnFFDBox() != 0))
      cout << endl <<"----------------- FFD technique (parametric -> cartesian) ---------------" << endl;
    
    /*--- Loop over all the FFD boxes levels ---*/
    for (iLevel = 0; iLevel < GetnLevel(); iLevel++) {
      
      /*--- Loop over all FFD FFDBoxs ---*/
      for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
        
        /*--- Check the level of the FFD box ---*/
        if(FFDBox[iFFDBox]->GetLevel() == iLevel) {
          
          /*--- Compute the parametric coordinates of the child box
           control points (using the parent FFDBox)  ---*/
          for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
            FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
            for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
              if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
            SetParametricCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
          }
          
          /*--- Update the parametric coordinates if it is a child FFDBox ---*/
          if (iLevel > 0) UpdateParametricCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
          /*--- Apply the design variables to the control point position ---*/
          for (iDV = 0; iDV < config->GetnDV(); iDV++) {
            switch ( config->GetDesign_Variable(iDV) ) {
              case FFD_CONTROL_POINT : SetFFDCPChange(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_DIHEDRAL_ANGLE : SetFFDDihedralAngle(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_TWIST_ANGLE : SetFFDTwistAngle(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_ROTATION : SetFFDRotation(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_CAMBER : SetFFDCamber(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_THICKNESS : SetFFDThickness(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
              case FFD_VOLUME : SetFFDVolume(geometry, config, FFDBox[iFFDBox], iFFDBox, iDV, false); break;
            }
          }
          
          /*--- Recompute cartesian coordinates using the new control point location ---*/
          SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
          
          /*--- Reparametrization of the parent FFD box ---*/
          for (iParent = 0; iParent < FFDBox[iFFDBox]->GetnParentFFDBox(); iParent++) {
            FFDBoxTag = FFDBox[iFFDBox]->GetParentFFDBoxTag(iParent);
            for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
              if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
            UpdateParametricCoord(geometry, config, FFDBox[jFFDBox], jFFDBox);
          }
          
          /*--- Compute the new location of the control points of the child boxes
           (using the parent FFDBox) ---*/
          for (iChild = 0; iChild < FFDBox[iFFDBox]->GetnChildFFDBox(); iChild++) {
            FFDBoxTag = FFDBox[iFFDBox]->GetChildFFDBoxTag(iChild);
            for (jFFDBox = 0; jFFDBox < GetnFFDBox(); jFFDBox++)
              if (FFDBoxTag == FFDBox[jFFDBox]->GetTag()) break;
            GetCartesianCoordCP(geometry, config, FFDBox[iFFDBox], FFDBox[jFFDBox]);
          }
        }
      }
      
      /*--- Output the deformed FFDBoxs ---*/
      if (rank == MASTER_NODE) {
        for (iFFDBox = 0; iFFDBox < GetnFFDBox(); iFFDBox++) {
          sprintf (buffer_char, "deformed_FFDBox.plt");
          if (iFFDBox == 0) FFDBox[iFFDBox]->SetTecplot(buffer_char, true);
          else FFDBox[iFFDBox]->SetTecplot(buffer_char, false);
        }
      }
      
    }
	}
  
}

void CSurfaceMovement::CopyBoundary(CGeometry *geometry, CConfig *config) {
	unsigned short iMarker;
	unsigned long iVertex, iPoint;
	double *Coord;

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
			iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
			Coord = geometry->node[iPoint]->GetCoord();
			geometry->vertex[iMarker][iVertex]->SetCoord(Coord);
		}
}

void CSurfaceMovement::SetParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint;
	double *car_coord, *car_coord_new, *par_coord, guess[3], max_diff, 
	my_max_diff = 0.0, diff;
	int rank;
	
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
	guess[0] = 0.5; guess[1] = 0.5; guess[2] = 0.5;
		
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_DV(iMarker) == YES)
			for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
				car_coord = geometry->vertex[iMarker][iVertex]->GetCoord();
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				
				/*--- If the point is inside the FFD, compute the value of the parametric coordinate ---*/
				if (FFDBox->GetPointFFD(geometry, config, iPoint)) {
					
					/*--- Find the parametric coordinate ---*/
					par_coord = FFDBox->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
					
					/*--- If the parametric coordinates are in (0,1) the point belongs to the FFDBox ---*/
					if (((par_coord[0] >= - EPS) && (par_coord[0] <= 1.0 + EPS)) && 
							((par_coord[1] >= - EPS) && (par_coord[1] <= 1.0 + EPS)) && 
							((par_coord[2] >= - EPS) && (par_coord[2] <= 1.0 + EPS))) {
						
						/*--- Set the value of the parametric coordinate ---*/
						FFDBox->Set_MarkerIndex(iMarker);
						FFDBox->Set_VertexIndex(iVertex);
						FFDBox->Set_PointIndex(iPoint);
						FFDBox->Set_ParametricCoord(par_coord);
						FFDBox->Set_CartesianCoord(car_coord);						
						
						/*--- Compute the cartesian coordinates using the parametric coordinates 
						 to check that everithing is right ---*/
						car_coord_new = FFDBox->EvalCartesianCoord(par_coord);
						
						/*--- Compute max difference between original value and the recomputed value ---*/
						diff = 0.0; 
						for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
							diff += (car_coord_new[iDim]-car_coord[iDim])*(car_coord_new[iDim]-car_coord[iDim]);
						diff = sqrt(diff);
						my_max_diff = max(my_max_diff, diff);
						
						guess[0] = par_coord[0]; guess[1] = par_coord[1]; guess[2] = par_coord[2];
					}
				}
			}
		
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Allreduce(&my_max_diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#endif
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Compute parametric coord      | FFD box: " << FFDBox->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetParametricCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild) {
	unsigned short iOrder, jOrder, kOrder;
	double *car_coord, *par_coord, guess[3];
	int rank;

#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
	for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
				car_coord = FFDBoxChild->GetCoordControlPoints(iOrder, jOrder, kOrder);
				par_coord = FFDBoxParent->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
				FFDBoxChild->SetParCoordControlPoints(par_coord, iOrder, jOrder, kOrder);
			}

	if (rank == MASTER_NODE)
		cout << "Compute parametric coord (CP) | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() <<"."<< endl;


}

void CSurfaceMovement::GetCartesianCoordCP(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBoxParent, CFreeFormDefBox *FFDBoxChild) {
	unsigned short iOrder, jOrder, kOrder, iDim;
	double *car_coord, *par_coord;
	int rank;
	
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
		
	for (iOrder = 0; iOrder < FFDBoxChild->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < FFDBoxChild->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < FFDBoxChild->GetnOrder(); kOrder++) {
				par_coord = FFDBoxChild->GetParCoordControlPoints(iOrder, jOrder, kOrder);
				
				/*--- Clip the value of the parametric coordinates (just in case)  ---*/
				for (iDim = 0; iDim < 3; iDim++) {
					if (par_coord[iDim] >= 1.0) par_coord[iDim] = 1.0;
					if (par_coord[iDim] <= 0.0) par_coord[iDim] = 0.0;
				}

				car_coord = FFDBoxParent->EvalCartesianCoord(par_coord);
				FFDBoxChild->SetCoordControlPoints(car_coord, iOrder, jOrder, kOrder);
			}
	
	if (rank == MASTER_NODE)
		cout << "Update cartesian coord (CP)   | FFD parent box: " << FFDBoxParent->GetTag() << ". FFD child box: " << FFDBoxChild->GetTag() <<"."<< endl;

}


void CSurfaceMovement::UpdateParametricCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint, iSurfacePoints;
	double car_coord[3], *car_coord_new, *car_coord_old, *par_coord, *var_coord, guess[3], max_diff, 
	my_max_diff = 0.0, diff;
	int rank;
	
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
			
	/*--- Recompute the parametric coordinates ---*/
	for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
		
		/*--- Get the marker of the surface point ---*/
		iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);
		
		if (config->GetMarker_All_DV(iMarker) == YES) {
			
			/*--- Get the vertex of the surface point ---*/
			iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
			iPoint = FFDBox->Get_PointIndex(iSurfacePoints);
	
			/*--- Get the parametric and cartesians coordinates of the 
			 surface point (they don't mach) ---*/
			par_coord = FFDBox->Get_ParametricCoord(iSurfacePoints);
			
			/*--- Compute and set the cartesian coord using the variation computed 
			 with the previous deformation ---*/
			var_coord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
			car_coord_old = geometry->node[iPoint]->GetCoord();
			for (iDim = 0; iDim < 3; iDim++)
				car_coord[iDim] = car_coord_old[iDim] + var_coord[iDim];
			FFDBox->Set_CartesianCoord(car_coord, iSurfacePoints);

			/*--- Find the parametric coordinate using as guess the previous value ---*/	
			guess[0] = par_coord[0]; guess[1] = par_coord[1]; guess[2] = par_coord[2];
			par_coord = FFDBox->GetParametricCoord_Iterative(car_coord, guess, 1E-10, 99999);
					
			/*--- Set the new value of the parametric coordinates ---*/
			FFDBox->Set_ParametricCoord(par_coord, iSurfacePoints);
			
			/*--- Compute the cartesian coordinates using the parametric coordinates 
			 to check that everithing is right ---*/
			car_coord_new = FFDBox->EvalCartesianCoord(par_coord);
			
			/*--- Compute max difference between original value and the recomputed value ---*/
			diff = 0.0; 
			for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
				diff += (car_coord_new[iDim]-car_coord[iDim])*(car_coord_new[iDim]-car_coord[iDim]);
			diff = sqrt(diff);
			my_max_diff = max(my_max_diff, diff);
				
		}
	}
		
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Allreduce(&my_max_diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#endif
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Update parametric coord       | FFD box: " << FFDBox->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetCartesianCoord(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox) {
	double *car_coord_old, *car_coord_new, diff, my_max_diff = 0.0, max_diff,
	*par_coord, VarCoord[3];
	unsigned short iMarker, iDim;
	unsigned long iVertex, iPoint, iSurfacePoints;
	int rank;
	
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
	/*--- Recompute the cartesians coordinates ---*/
	for (iSurfacePoints = 0; iSurfacePoints < FFDBox->GetnSurfacePoint(); iSurfacePoints++) {
		
		/*--- Get the marker of the surface point ---*/
		iMarker = FFDBox->Get_MarkerIndex(iSurfacePoints);
		
		if (config->GetMarker_All_DV(iMarker) == YES) {
			
			/*--- Get the vertex of the surface point ---*/
			iVertex = FFDBox->Get_VertexIndex(iSurfacePoints);
			iPoint = FFDBox->Get_PointIndex(iSurfacePoints);

			/*--- Set to zero the variation of the coordinates ---*/
			for (iDim = 0; iDim < 3; iDim++) VarCoord[iDim] = 0.0;
			geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

			/*--- Get the parametric coordinate of the surface point ---*/
			par_coord = FFDBox->Get_ParametricCoord(iSurfacePoints);
			
			/*--- Compute the new cartesian coordinate, and set the value in 
			 the FFDBox structure ---*/
			car_coord_new = FFDBox->EvalCartesianCoord(par_coord);
			FFDBox->Set_CartesianCoord(car_coord_new, iSurfacePoints);
			
			/*--- Get the original cartesian coordinates of the surface point ---*/
			car_coord_old = geometry->node[iPoint]->GetCoord();

			/*--- Set the value of the variation of the coordinates ---*/
			for (iDim = 0; iDim < 3; iDim++) {
				VarCoord[iDim] = car_coord_new[iDim] - car_coord_old[iDim];
				if (fabs(VarCoord[iDim]) < EPS) VarCoord[iDim] = 0.0;
			}
			
			diff = sqrt((car_coord_new[0]-car_coord_old[0])*(car_coord_new[0]-car_coord_old[0]) +
									(car_coord_new[1]-car_coord_old[1])*(car_coord_new[1]-car_coord_old[1]) +
									(car_coord_new[2]-car_coord_old[2])*(car_coord_new[2]-car_coord_old[2]));
			
			my_max_diff = max(my_max_diff, diff);
			
			/*--- Set the variation of the coordinates ---*/
			geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

		}
	}
		
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Allreduce(&my_max_diff, &max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Allreduce(&my_max_diff, &max_diff, 1, MPI::DOUBLE, MPI::MAX); 	
#endif
#else
	max_diff = my_max_diff;
#endif
	
	if (rank == MASTER_NODE) 
		cout << "Update cartesian coord        | FFD box: " << FFDBox->GetTag() << ". Max diff: " << max_diff <<"."<< endl;
	
}

void CSurfaceMovement::SetFFDCPChange(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																			unsigned short iDV, bool ResetDef) {
	
	double movement[3], Ampl;
	unsigned short design_FFDBox, index[3];
		
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		Ampl = config->GetDV_Value(iDV);
		
		index[0] = int(config->GetParamDV(iDV, 1));
		index[1] = int(config->GetParamDV(iDV, 2)); 
		index[2] = int(config->GetParamDV(iDV, 3));
		
		movement[0] = config->GetParamDV(iDV, 4)*Ampl; 
		movement[1] = config->GetParamDV(iDV, 5)*Ampl; 
		movement[2] = config->GetParamDV(iDV, 6)*Ampl;
		
		if (ResetDef == true) FFDBox->SetOriginalControlPoints();
		FFDBox->SetControlPoints(index, movement);
		
	}
		
}

void CSurfaceMovement::SetFFDCamber(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																		unsigned short iDV, bool ResetDef) {
	double Ampl, movement[3];
	unsigned short design_FFDBox, index[3], kIndex;
	
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- Compute the variation of the design variable ---*/
		for (kIndex = 0; kIndex < 2; kIndex++) {
						
			Ampl = config->GetDV_Value(iDV);
			
			design_FFDBox = int(config->GetParamDV(iDV, 0));
			if (design_FFDBox > nFFDBox) { cout <<"The FFDBox ID is bigger than the number of FFDBoxs!!"<< endl; exit(1); }
			
			index[0] = int(config->GetParamDV(iDV, 1));
			index[1] = int(config->GetParamDV(iDV, 2)); 
			index[2] = kIndex;
			
			movement[0] = 0.0; movement[1] = 0.0; 
			if (kIndex == 0) movement[2] = Ampl;
			else movement[2] = Ampl;
			
			if (ResetDef == true) FFDBox->SetOriginalControlPoints();
			FFDBox->SetControlPoints(index, movement);
		}
		
	}
	
}

void CSurfaceMovement::SetFFDThickness(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																			 unsigned short iDV, bool ResetDef) {
	double Ampl, movement[3];
	unsigned short design_FFDBox, index[3], kIndex;
		
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- Compute the variation of the design variable ---*/
		for (kIndex = 0; kIndex < 2; kIndex++) {
			
			Ampl = config->GetDV_Value(iDV);
			
			design_FFDBox = int(config->GetParamDV(iDV, 0));
			
			index[0] = int(config->GetParamDV(iDV, 1));
			index[1] = int(config->GetParamDV(iDV, 2)); 
			index[2] = kIndex;
			
			movement[0] = 0.0; movement[1] = 0.0; 
			if (kIndex == 0) movement[2] = -Ampl;
			else movement[2] = Ampl;
			
			if (ResetDef == true) FFDBox->SetOriginalControlPoints();
			FFDBox->SetControlPoints(index, movement);
		}
		
	}
	
}

void CSurfaceMovement::SetFFDVolume(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																			 unsigned short iDV, bool ResetDef) {
	double Ampl, movement[3];
	unsigned short design_FFDBox, index[3];
			
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- Compute the variation of the design variable ---*/
		Ampl = config->GetDV_Value(iDV);
				
		index[0] = int(config->GetParamDV(iDV, 1));
		index[1] = int(config->GetParamDV(iDV, 2)); 
		index[2] = 0;
		
		movement[0] = 0.0; movement[1] = 0.0; 
		movement[2] = Ampl;
		
		if (ResetDef == true) FFDBox->SetOriginalControlPoints();
		FFDBox->SetControlPoints(index, movement);
		
	}
	
}


void CSurfaceMovement::SetFFDDihedralAngle(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																					 unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder, design_FFDBox, index[3];
	double movement[3];
			
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- The angle of rotation. ---*/
		double theta = config->GetDV_Value(iDV)*PI_NUMBER/180.0;
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
					movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1]*tan(theta);
					
					if (ResetDef == true) FFDBox->SetOriginalControlPoints();
					FFDBox->SetControlPoints(index, movement);
				}
		
	}

}

void CSurfaceMovement::SetFFDTwistAngle(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																				unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder;
	double  x, y, z, movement[3];
	unsigned short index[3], design_FFDBox;
	
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- xyz-coordinates of a point on the line of rotation. ---*/
		double a = config->GetParamDV(iDV, 1);
		double b = config->GetParamDV(iDV, 2);
		double c = config->GetParamDV(iDV, 3);
		
    /*--- xyz-coordinate of the line's direction vector. ---*/
		double u = config->GetParamDV(iDV, 4)-config->GetParamDV(iDV, 1);
		double v = config->GetParamDV(iDV, 5)-config->GetParamDV(iDV, 2);
		double w = config->GetParamDV(iDV, 6)-config->GetParamDV(iDV, 3);
		
		/*--- The angle of rotation. ---*/
		double theta = config->GetDV_Value(iDV)*PI_NUMBER/180.0;
		
		/*--- An intermediate value used in computations. ---*/
		double u2=u*u; double v2=v*v; double w2=w*w;     
		double l2 = u2 + v2 + w2; double l = sqrt(l2);
		double cosT; double sinT;  
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
					x = coord[0]; y = coord[1]; z = coord[2];
					
					double factor = 0.0; 
					if ( z < config->GetParamDV(iDV, 3) )
						factor = 0.0;
					if (( z >= config->GetParamDV(iDV, 3)) && ( z <= config->GetParamDV(iDV, 6)) )
						factor = (z-config->GetParamDV(iDV, 3)) / (config->GetParamDV(iDV, 6)-config->GetParamDV(iDV, 3));
					if ( z > config->GetParamDV(iDV, 6) )
						factor = 1.0;
					
					cosT = cos(theta*factor); 
					sinT = sin(theta*factor);  
					
					movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
					+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
					+ l*(-c*v + b*w - w*y + v*z)*sinT;
					movement[0] = movement[0]/l2 - x;
					
					movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z) 
					+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
					+ l*(c*u - a*w + w*x - u*z)*sinT;
					movement[1] = movement[1]/l2 - y;
					
					movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z) 
					+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
					+ l*(-b*u + a*v - v*x + u*y)*sinT;
					movement[2] = movement[2]/l2 - z;
					
					if (ResetDef == true) FFDBox->SetOriginalControlPoints();
					FFDBox->SetControlPoints(index, movement);		
				}
		
	}
	
}


void CSurfaceMovement::SetFFDRotation(CGeometry *geometry, CConfig *config, CFreeFormDefBox *FFDBox, unsigned short iFFDBox, 
																			unsigned short iDV, bool ResetDef) {
	unsigned short iOrder, jOrder, kOrder;
	double  movement[3], x, y, z;
	unsigned short index[3], design_FFDBox;
		
	design_FFDBox = int(config->GetParamDV(iDV, 0));
	
	if (design_FFDBox == iFFDBox) {
		
		/*--- xyz-coordinates of a point on the line of rotation. ---*/
		double a = config->GetParamDV(0,1);
		double b = config->GetParamDV(0,2);
		double c = config->GetParamDV(0,3);
		
		/*--- xyz-coordinate of the line's direction vector. ---*/
		double u = config->GetParamDV(0,4)-config->GetParamDV(0,1);
		double v = config->GetParamDV(0,5)-config->GetParamDV(0,2);
		double w = config->GetParamDV(0,6)-config->GetParamDV(0,3);
		
		/*--- The angle of rotation. ---*/
		double theta = config->GetDV_Value(0)*PI_NUMBER/180.0;
		
		/*--- An intermediate value used in computations. ---*/
		double u2=u*u; double v2=v*v; double w2=w*w;     
		double cosT = cos(theta); double sinT = sin(theta);  
		double l2 = u2 + v2 + w2; double l = sqrt(l2);
		
		/*--- Change the value of the control point if move is true ---*/
		for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
			for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
				for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
					index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
					double *coord = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
					x = coord[0]; y = coord[1]; z = coord[2];
					movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
					+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
					+ l*(-c*v + b*w - w*y + v*z)*sinT;
					movement[0] = movement[0]/l2 - x;
					
					movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z) 
					+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
					+ l*(c*u - a*w + w*x - u*z)*sinT;
					movement[1] = movement[1]/l2 - y;
					
					movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z) 
					+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
					+ l*(-b*u + a*v - v*x + u*y)*sinT;
					movement[2] = movement[2]/l2 - z;
					
					if (ResetDef == true) FFDBox->SetOriginalControlPoints();
					FFDBox->SetControlPoints(index, movement);		
				}
		
	}
	
}

void CSurfaceMovement::SetHicksHenne(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], VarCoord_[3], *Coord_, *Normal_, ek, fk, BumpSize = 1.0, BumpLoc = 0.0, Coord[3], Normal[3],
  xCoord, TPCoord[2] = {0.0, 0.0}, LPCoord[2] = {0.0, 0.0}, Distance, Chord, AoA, ValCos, ValSin;
  
	bool upper = true, double_surface = false;

	/*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
  
  /*--- Compute the angle of attack to apply the deformation ---*/
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      Coord_ = boundary->vertex[iMarker][0]->GetCoord();
      TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1];
      for (iVertex = 1; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        if (Coord_[0] > TPCoord[0]) { TPCoord[0] = Coord_[0]; TPCoord[1] = Coord_[1]; }
			}
		}
	}
  
#ifndef NO_MPI

  unsigned long *Buffer_Send_nVertex, *Buffer_Receive_nVertex;
	int iProcessor, nProcessor;
	double *Buffer_Send_Coord, *Buffer_Receive_Coord;

#ifdef WINDOWS
	MPI_Comm_size(MPI_COMM_WORLD,&nProcessor);
#else
	nProcessor = MPI::COMM_WORLD.Get_size();
#endif
  
	Buffer_Receive_Coord = new double [nProcessor*2];
  Buffer_Send_Coord = new double [2];
  
  Buffer_Send_Coord[0] = TPCoord[0]; Buffer_Send_Coord[1] = TPCoord[1];

#ifdef WINDOWS
	MPI_Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Allgather(Buffer_Send_Coord, 2, MPI::DOUBLE, Buffer_Receive_Coord, 2, MPI::DOUBLE);
#endif

  TPCoord[0] = Buffer_Receive_Coord[0]; TPCoord[1] = Buffer_Receive_Coord[1];
  for (iProcessor = 1; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    if (Coord[0] > TPCoord[0]) { TPCoord[0] = Coord[0]; TPCoord[1] = Coord[1]; }
  }
  
	delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif


  Chord = 0.0;
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
        Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
        Distance = sqrt(pow(Coord_[0] - TPCoord[0], 2.0) + pow(Coord_[1] - TPCoord[1], 2.0));
        if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord_[0]; LPCoord[1] = Coord_[1]; }
      }
    }
  }
  
#ifndef NO_MPI
 
#ifdef WINDOWS
	MPI_Comm_size(MPI_COMM_WORLD,&nProcessor);
#else
	nProcessor = MPI::COMM_WORLD.Get_size();
#endif
  
	Buffer_Receive_Coord = new double [nProcessor*2];
  Buffer_Send_Coord = new double [2];
  
  Buffer_Send_Coord[0] = LPCoord[0]; Buffer_Send_Coord[1] = LPCoord[1];

#ifdef WINDOWS
	MPI_Allgather(Buffer_Send_Coord, 2, MPI_DOUBLE, Buffer_Receive_Coord, 2, MPI_DOUBLE, MPI_COMM_WORLD);
#else
	MPI::COMM_WORLD.Allgather(Buffer_Send_Coord, 2, MPI::DOUBLE, Buffer_Receive_Coord, 2, MPI::DOUBLE);
#endif
  
  Chord = 0.0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    Coord[0] = Buffer_Receive_Coord[iProcessor*2 + 0];
    Coord[1] = Buffer_Receive_Coord[iProcessor*2 + 1];
    Distance = sqrt(pow(Coord[0] - TPCoord[0], 2.0) + pow(Coord[1] - TPCoord[1], 2.0));
    if (Chord < Distance) { Chord = Distance; LPCoord[0] = Coord[0]; LPCoord[1] = Coord[1]; }
  }
  
	delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  
#endif
  
  AoA = atan((LPCoord[1] - TPCoord[1]) / (TPCoord[0] - LPCoord[0]))*180/PI_NUMBER;
  
	/*--- Perform multiple airfoil deformation ---*/
  
	double Ampl = config->GetDV_Value(iDV);
	double xk = config->GetParamDV(iDV, 1);
	const double t2 = 3.0;
  
	if (config->GetParamDV(iDV, 0) == NO)  { upper = false; double_surface = true; }
	if (config->GetParamDV(iDV, 0) == YES) { upper = true; double_surface = true; }
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      
      
			if (config->GetMarker_All_DV(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
        
				Coord_ = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal_ = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- The Hicks Henne bump functions should be applied to a basic airfoil without AoA,
         and unitary chord, a tranformation is required ---*/
        ValCos = cos(AoA*PI_NUMBER/180.0);
        ValSin = sin(AoA*PI_NUMBER/180.0);
        
        Coord[0] = Coord_[0]*ValCos - Coord_[1]*ValSin;
        Coord[0] = max(0.0, Coord[0]); // Coord x should be always positive
        Coord[1] = Coord_[1]*ValCos + Coord_[0]*ValSin;
        
        Normal[0] = Normal_[0]*ValCos - Normal_[1]*ValSin;
        Normal[1] = Normal_[1]*ValCos + Normal_[0]*ValSin;
        
        /*--- Bump computation ---*/
				if (double_surface) {
					ek = log10(0.5)/log10(xk);
					fk = pow( sin( PI_NUMBER * pow(Coord[0],ek) ) , t2);
					/*--- Upper and lower surface ---*/
					if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }
					if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }
				}
				else {
					xCoord = Coord[0] - BumpLoc;
					ek = log10(0.5)/log10(xk/BumpSize);
					fk = pow( sin( PI_NUMBER * pow(xCoord/BumpSize,ek)),t2);
          
					/*--- Only one surface ---*/
					if ((xCoord <= 0.0) || (xCoord >= BumpSize)) VarCoord[1] =  0.0;
					else VarCoord[1] =  Ampl*fk;
				}
			}
      
      /*--- Apply the transformation to the coordinate variation ---*/
      ValCos = cos(-AoA*PI_NUMBER/180.0);
      ValSin = sin(-AoA*PI_NUMBER/180.0);
      
      VarCoord_[0] = VarCoord[0]*ValCos - VarCoord[1]*ValSin;
      VarCoord_[1] = VarCoord[1]*ValCos + VarCoord[0]*ValSin;

			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord_);
		}
	}
  
}

void CSurfaceMovement::SetSpherical(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
  
	unsigned long iVertex, iPoint, n;
	unsigned short iMarker, jDV;
	double VarCoord[3], *Coord, *Normal, Theta_Value, Radius_Value, Delta;
	double x, x2, y2, z2, r_yz, r_yz2, theta, r, cos_theta, sin_theta, cos_phi, sin_phi;
	vector<double> Theta_Spline, Radius_Spline, Radius2_Spline;
	int ControlPoint_Index;
	ofstream File;
  
  /*--- Read the baseline spline ---*/
  Theta_Spline.push_back(0.0);              Radius_Spline.push_back(0.1524);
  Theta_Spline.push_back(0.1963495408494);  Radius_Spline.push_back(0.1524);
  Theta_Spline.push_back(0.3926990816987);  Radius_Spline.push_back(0.1524);
  Theta_Spline.push_back(0.7853981633974);  Radius_Spline.push_back(0.1524);
  Theta_Spline.push_back(1.4137166941154);  Radius_Spline.push_back(0.1524);
  Theta_Spline.push_back(1.65766545);   Radius_Spline.push_back(0.15704997);
  
  
  /*--- Read the baseline spline ---*/
  //  Theta_Spline.push_back(0.0);              Radius_Spline.push_back(5.9);
  //  Theta_Spline.push_back(0.1963495408494);  Radius_Spline.push_back(5.0);
  //  Theta_Spline.push_back(0.3926990816987);  Radius_Spline.push_back(4.0);
  //  Theta_Spline.push_back(0.7853981633974);  Radius_Spline.push_back(3.0);
  //  Theta_Spline.push_back(1.570796326795);   Radius_Spline.push_back(2.6);
  
  if (ResetDef == true) { /*--- Only one deformation (iVD), reseting the spline at every design variable ---*/
    
    /*--- Modify the baseline control point positions using the the deformation ---*/
    ControlPoint_Index = int(config->GetParamDV(iDV, 0));
    
    Theta_Value = config->GetParamDV(iDV, 1);
    Radius_Value = config->GetParamDV(iDV, 2);
    Delta = config->GetDV_Value(iDV);
    
    Theta_Spline[ControlPoint_Index] += Delta*Theta_Value;
    Radius_Spline[ControlPoint_Index] += Delta*Radius_Value;
    
  }
  
  else { /*--- Aditive deformation, change all the points of the original spline,
          using all the design variables. ---*/
    
    for (jDV = 0; jDV < config->GetnDV(); jDV++) {
      
      /*--- Modify the baseline control point positions using the the deformation ---*/
      ControlPoint_Index = int(config->GetParamDV(jDV, 0));
      
      Theta_Value = config->GetParamDV(jDV, 1);
      Radius_Value = config->GetParamDV(jDV, 2);
      Delta = config->GetDV_Value(jDV);
      
      Theta_Spline[ControlPoint_Index] += Delta*Theta_Value;
      Radius_Spline[ControlPoint_Index] += Delta*Radius_Value;
      
    }
    
  }
  
	/*--- Create the spline ---*/
	n = Theta_Spline.size();
	Radius2_Spline.resize(n);
	boundary->SetSpline(Theta_Spline, Radius_Spline, n, 0, 0, Radius2_Spline);
  
	/*--- Loop over all the points in the surface to evaluate the coordinate change ---*/
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      
			if (config->GetMarker_All_DV(iMarker) == YES) {
        
        iPoint = boundary->vertex[iMarker][iVertex]->GetNode();
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Compute the coordinate variation due to a surface modification, given
         ControlPoint_Index, Theta, Radius, and Ampl ---*/
        
        //       x = Coord[0] - (6.5 + 2.9205);	// HARD-CODED VALUE; needs to be fixed (MC)
        //       if (x > 0) {
        
        x = 0.1524 - Coord[0];	// HARD-CODED VALUE; needs to be fixed (MC)
        //0.12855939
        if ((Coord[0] >= 0.0) && (Coord[0] <= 0.16602564))  {
          
          /*--- Compute the coordinate variation due to a surface modification, given
           ControlPoint_Index, Theta, Radius, and Ampl ---*/
          if (Coord[1] == 0.0 && Coord[2] == 0.0) {	// on axis
            r = boundary->GetSpline(Theta_Spline, Radius_Spline, Radius2_Spline, n, 0.0);
            VarCoord[0] = r - x;
            VarCoord[1] = 0.0;
            VarCoord[2] = 0.0;
          }
          else {
            x2 = x*x; y2 = Coord[1]*Coord[1]; z2 = Coord[2]*Coord[2];
            r_yz = sqrt(y2 + z2); r_yz2 = y2 + z2;
            theta = atan(r_yz/x);
            cos_theta = x/sqrt(x2 + r_yz2);
            sin_theta = r_yz/sqrt(x2 + r_yz2);
            cos_phi = Coord[1]/sqrt(z2 + y2);
            sin_phi = Coord[2]/sqrt(z2 + y2);
            r = boundary->GetSpline(Theta_Spline, Radius_Spline, Radius2_Spline, n, theta);
            VarCoord[0] = r*cos_theta - x;
            VarCoord[1] = r*sin_theta*cos_phi - Coord[1];
            VarCoord[2] = r*sin_theta*sin_phi - Coord[2];
          }
          
        }
        
				/*--- Set the coordinate variation due to a surface modification ---*/
        
				boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
        
			}
		}
	}
  
  
  Theta_Spline.clear();
  Radius_Spline.clear();
  Radius2_Spline.clear();
  
}

void CSurfaceMovement::SetRotation(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex;
	unsigned short iMarker;
	double VarCoord[3], *Coord;
	double  movement[3], x, y, z;
  
  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
  
	/*--- xyz-coordinates of a point on the line of rotation. */
  
	double a = config->GetParamDV(iDV, 0);
	double b = config->GetParamDV(iDV, 1);
	double c = 0.0;
	if (boundary->GetnDim() == 3) c = config->GetParamDV(0,2);
	
	/*--- xyz-coordinate of the line's direction vector. ---*/
  
	double u = config->GetParamDV(iDV, 3)-config->GetParamDV(iDV, 0);
	double v = config->GetParamDV(iDV, 4)-config->GetParamDV(iDV, 1);
	double w = 1.0;
	if (boundary->GetnDim() == 3) w = config->GetParamDV(iDV, 5)-config->GetParamDV(iDV, 2);
	
	/*--- The angle of rotation. ---*/
  
	double theta = config->GetDV_Value(iDV)*PI_NUMBER/180.0;
	
	/*--- An intermediate value used in computations. ---*/
  
	double u2=u*u; double v2=v*v; double w2=w*w;
	double cosT = cos(theta); double sinT = sin(theta);
	double l2 = u2 + v2 + w2; double l = sqrt(l2);
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				x = Coord[0]; y = Coord[1]; z = Coord[2];
				
				movement[0] = a*(v2 + w2) + u*(-b*v - c*w + u*x + v*y + w*z)
				+ (-a*(v2 + w2) + u*(b*v + c*w - v*y - w*z) + (v2 + w2)*x)*cosT
				+ l*(-c*v + b*w - w*y + v*z)*sinT;
				movement[0] = movement[0]/l2 - x;
				
				movement[1] = b*(u2 + w2) + v*(-a*u - c*w + u*x + v*y + w*z)
				+ (-b*(u2 + w2) + v*(a*u + c*w - u*x - w*z) + (u2 + w2)*y)*cosT
				+ l*(c*u - a*w + w*x - u*z)*sinT;
				movement[1] = movement[1]/l2 - y;
				
				movement[2] = c*(u2 + v2) + w*(-a*u - b*v + u*x + v*y + w*z)
				+ (-c*(u2 + v2) + w*(a*u + b*v - u*x - v*y) + (u2 + v2)*z)*cosT
				+ l*(-b*u + a*v - v*x + u*y)*sinT;
				if (boundary->GetnDim() == 3) movement[2] = movement[2]/l2 - z;
				else movement[2] = 0.0;
				
				VarCoord[0] = movement[0];
				VarCoord[1] = movement[1];
				if (boundary->GetnDim() == 3) VarCoord[2] = movement[2];
				
			}
			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetDisplacement(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex;
	unsigned short iMarker;
	double VarCoord[3];
	double Ampl = config->GetDV_Value(0);
	
  /*--- Reset airfoil deformation if first deformation or if it required by the solver ---*/
  
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
  
	double xDispl = config->GetParamDV(iDV, 0);
	double yDispl = config->GetParamDV(iDV, 1);
	double zDispl = 0;
	if (boundary->GetnDim() == 3) zDispl = config->GetParamDV(iDV, 2);
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
				VarCoord[0] = Ampl*xDispl;
				VarCoord[1] = Ampl*yDispl;
				if (boundary->GetnDim() == 3) VarCoord[2] = Ampl*zDispl;
			}
			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
		}
  
}

void CSurfaceMovement::SetCosBump(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal, fk, DesignSize = 2.0, DesignLoc = 1.0, xCoord, xCoord_Local;
  
	bool upper = true, double_surface = false;
  
	/*--- Reset airfoil deformation if first deformation ---*/
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
  
	/*--- Perform multiple airfoil deformation ---*/
	double Ampl = config->GetDV_Value(iDV);
	double BumpCenter = DesignLoc + config->GetParamDV(iDV, 1)*DesignSize;
	double BumpSize = config->GetParamDV(iDV, 2);
  
	if (config->GetParamDV(iDV, 0) == NO)  { upper = false; double_surface = true; }
	if (config->GetParamDV(iDV, 0) == YES) { upper = true; double_surface = true; }
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
        
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Bump computation ---*/
				if (double_surface) {
					xCoord = Coord[0];
          xCoord_Local = (xCoord - BumpCenter);
          
          if (fabs(xCoord_Local) < BumpSize) fk = 0.5*(1.0+cos(PI_NUMBER*xCoord_Local/BumpSize));
          else fk = 0.0;
					/*--- Upper and lower surface ---*/
					if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }
					if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }
				}
				else {
          
					xCoord = Coord[0];
          xCoord_Local = (xCoord - BumpCenter);
          
          //          if (fabs(xCoord_Local) < BumpSize*4.0) fk = 0.5*(1.0+cos(PI_NUMBER*xCoord_Local/BumpSize))*
          //            0.5*(1.0+cos(PI_NUMBER*xCoord_Local/(BumpSize*2.0)))*
          //            0.5*(1.0+cos(PI_NUMBER*xCoord_Local/(BumpSize*4.0)));
          //          else fk = 0.0;
          
          if (fabs(xCoord_Local) < BumpSize) fk = 0.5*(1.0+cos(PI_NUMBER*xCoord_Local/BumpSize));
          else fk = 0.0;
          
          //					xCoord_Local = (xCoord - BumpCenter)/(0.5*BumpSize);
          //					ek = -1.0 / (1.0 - xCoord_Local*xCoord_Local);
          //          if (fabs(xCoord_Local) < 1.0) fk = exp(ek)/exp(-1.0);
          //          else fk = 0.0;
          
					/*--- Only one surface ---*/
					VarCoord[1] =  Ampl*fk;
				}
			}
      
			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
		}
	}
}

void CSurfaceMovement::SetFourier(CGeometry *boundary, CConfig *config, unsigned short iDV, bool ResetDef) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal, fk, DesignSize = 2.0, DesignLoc = 1.0, xCoord, xCoord_Local;
  
	bool upper = true, double_surface = false;
  
	/*--- Reset airfoil deformation if first deformation ---*/
	if ((iDV == 0) || (ResetDef == true)) {
		for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
			for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
				VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
				boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
			}
	}
  
	/*--- Perform multiple airfoil deformation ---*/
	double Ampl = config->GetDV_Value(iDV);
  double T = DesignSize;
  double n = int(config->GetParamDV(iDV, 1));
  double omega = 2.0*PI_NUMBER/T;
  double omega_n = omega*n;
  
	if (config->GetParamDV(iDV, 0) == NO)  { upper = false; double_surface = true; }
	if (config->GetParamDV(iDV, 0) == YES) { upper = true; double_surface = true; }
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
        
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Bump computation ---*/
				if (double_surface) {
          
          xCoord = Coord[0];
          xCoord_Local = (xCoord - (DesignLoc+0.5*DesignSize));
          
          if ((xCoord_Local < -0.5*T) || (xCoord_Local > 0.5*T)) fk = 0.0;
          else {
            if (n == 0) fk = 0.5;
            else {
              if (int(config->GetParamDV(iDV, 2)) == 0) fk = cos(omega_n*xCoord_Local);
              else fk = sin(omega_n*xCoord_Local);
            }
          }
          
					/*--- Upper and lower surface ---*/
					if (( upper) && (Normal[1] > 0)) { VarCoord[1] =  Ampl*fk; }
					if ((!upper) && (Normal[1] < 0)) { VarCoord[1] = -Ampl*fk; }
				}
				else {
          
					xCoord = Coord[0];
          xCoord_Local = (xCoord - (DesignLoc+0.5*DesignSize));
          
          if ((xCoord_Local < -0.5*T) || (xCoord_Local > 0.5*T)) fk = 0.0;
          else {
            if (n == 0) fk = 0.5;
            else {
              if (int(config->GetParamDV(iDV, 2)) == 0) fk = cos(omega_n*xCoord_Local);
              else fk = sin(omega_n*xCoord_Local);
            }
          }
          
					/*--- Only one surface ---*/
					VarCoord[1] =  Ampl*fk;
				}
			}
      
			boundary->vertex[iMarker][iVertex]->AddVarCoord(VarCoord);
		}
	}
}

void CSurfaceMovement::Moving_Walls(CGeometry *geometry, CConfig *config,
                                    unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Local variables ---*/
  unsigned short iMarker, jMarker, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  double xDot[3] = {0.0,0.0,0.0}, *Coord, Center[3], Omega[3], r[3], GridVel[3];
	double L_Ref     = config->GetLength_Ref();
  double Omega_Ref = config->GetOmega_Ref();
  double Vel_Ref   = config->GetVelocity_Ref();
	string Marker_Tag;
  
  /*--- Store grid velocity for each node on the moving surface(s).
   Sum and store the x, y, & z velocities due to translation and rotation. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) == YES) {
      
      /*--- Identify iMarker from the list of those under MARKER_MOVING ---*/
      
      Marker_Tag = config->GetMarker_All_Tag(iMarker);
      jMarker    = config->GetMarker_Moving(Marker_Tag);
      
      /*--- Get prescribed wall speed from config for this marker ---*/
      
      Center[0] = config->GetMotion_Origin_X(jMarker);
      Center[1] = config->GetMotion_Origin_Y(jMarker);
      Center[2] = config->GetMotion_Origin_Z(jMarker);
      Omega[0]  = config->GetRotation_Rate_X(jMarker)/Omega_Ref;
      Omega[1]  = config->GetRotation_Rate_Y(jMarker)/Omega_Ref;
      Omega[2]  = config->GetRotation_Rate_Z(jMarker)/Omega_Ref;
      xDot[0]   = config->GetTranslation_Rate_X(jMarker)/Vel_Ref;
      xDot[1]   = config->GetTranslation_Rate_Y(jMarker)/Vel_Ref;
      xDot[2]   = config->GetTranslation_Rate_Z(jMarker)/Vel_Ref;
      
      if (rank == MASTER_NODE && iter == 0) {
        cout << " Storing grid velocity for marker: ";
        cout << Marker_Tag << "." << endl;
        cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
        cout << ", " << xDot[2] << ") m/s." << endl;
        cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
        cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
        cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
      }
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        /*--- Get the index and coordinates of the current point ---*/
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord  = geometry->node[iPoint]->GetCoord();
        
        /*--- Calculate non-dim. position from rotation center ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          r[iDim] = (Coord[iDim]-Center[iDim])/L_Ref;
        if (nDim == 2) r[nDim] = 0.0;
        
        /*--- Cross Product of angular velocity and distance from center to
         get the rotational velocity. Note that we are adding on the velocity
         due to pure translation as well. ---*/
        
        GridVel[0] = xDot[0] + Omega[1]*r[2] - Omega[2]*r[1];
        GridVel[1] = xDot[1] + Omega[2]*r[0] - Omega[0]*r[2];
        GridVel[2] = xDot[2] + Omega[0]*r[1] - Omega[1]*r[0];
        
        /*--- Store the moving wall velocity for this node ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          geometry->node[iPoint]->SetGridVel(iDim,GridVel[iDim]);
  
      }
		}
	}
}

void CSurfaceMovement::Surface_Translating(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
  
	double deltaT, time_new, time_old, Lref;
  double Center[3], VarCoord[3], xDot[3];
  unsigned short iMarker, jMarker, Moving;
  unsigned long iVertex;
  string Marker_Tag, Moving_Tag;
  int rank;
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<double>(iter-1)*deltaT;
  }
  
	/*--- Store displacement of each node on the translating surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to translate ---*/
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving(jMarker);
        Marker_Tag = config->GetMarker_All_Tag(iMarker);
        
        if (Marker_Tag == Moving_Tag) {

          /*--- Translation velocity from config. ---*/
          
          xDot[0]   = config->GetTranslation_Rate_X(jMarker);
          xDot[1]   = config->GetTranslation_Rate_Y(jMarker);
          xDot[2]   = config->GetTranslation_Rate_Z(jMarker);
          
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.
          
          if (rank == MASTER_NODE) {
            cout << " Storing translating displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Translational velocity: (" << xDot[0] << ", " << xDot[1];
              cout << ", " << xDot[2] << ") m/s." << endl;
            }
          }
          
          /*--- Compute delta change in the position in the x, y, & z directions. ---*/
          
          VarCoord[0] = xDot[0]*(time_new-time_old);
          VarCoord[1] = xDot[1]*(time_new-time_old);
          VarCoord[2] = xDot[2]*(time_new-time_old);
          
          for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
        same translational velocity, because we use the last VarCoord set ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the translation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
      Center[0] = config->GetMotion_Origin_X(jMarker) + VarCoord[0];
      Center[1] = config->GetMotion_Origin_Y(jMarker) + VarCoord[1];
      Center[2] = config->GetMotion_Origin_Z(jMarker) + VarCoord[2];
      config->SetMotion_Origin_X(jMarker, Center[0]);
      config->SetMotion_Origin_Y(jMarker, Center[1]);
      config->SetMotion_Origin_Z(jMarker, Center[2]);
    }
  }
  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the translation. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }
}

void CSurfaceMovement::Surface_Plunging(CGeometry *geometry, CConfig *config,
                                           unsigned long iter, unsigned short iZone) {
  
	double deltaT, time_new, time_old, Lref;
  double Center[3], VarCoord[3], Omega[3], Ampl[3];
  double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iMarker, jMarker, Moving;
  unsigned long iVertex;
  string Marker_Tag, Moving_Tag;
  int rank;
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<double>(iter-1)*deltaT;
  }
  
	/*--- Store displacement of each node on the plunging surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to plunge ---*/
  
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving(jMarker);
        Marker_Tag = config->GetMarker_All_Tag(iMarker);
        
        if (Marker_Tag == Moving_Tag) {
          
          /*--- Plunging frequency and amplitude from config. ---*/
          
          Omega[0]  = config->GetPlunging_Omega_X(jMarker)/config->GetOmega_Ref();
          Omega[1]  = config->GetPlunging_Omega_Y(jMarker)/config->GetOmega_Ref();
          Omega[2]  = config->GetPlunging_Omega_Z(jMarker)/config->GetOmega_Ref();
          Ampl[0]   = config->GetPlunging_Ampl_X(jMarker)/Lref;
          Ampl[1]   = config->GetPlunging_Ampl_Y(jMarker)/Lref;
          Ampl[2]   = config->GetPlunging_Ampl_Z(jMarker)/Lref;
          
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing plunging displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Plunging frequency: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s." << endl;
              cout << " Plunging amplitude: (" << Ampl[0]/DEG2RAD;
              cout << ", " << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
              cout << ") degrees."<< endl;
            }
          }
          
          /*--- Compute delta change in the position in the x, y, & z directions. ---*/
          
          VarCoord[0] = -Ampl[0]*(sin(Omega[0]*time_new) - sin(Omega[0]*time_old));
          VarCoord[1] = -Ampl[1]*(sin(Omega[1]*time_new) - sin(Omega[1]*time_old));
          VarCoord[2] = -Ampl[2]*(sin(Omega[2]*time_new) - sin(Omega[2]*time_old));
          
          for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
   same plunging movement, because we use the last VarCoord set ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the translation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
      Center[0] = config->GetMotion_Origin_X(jMarker) + VarCoord[0];
      Center[1] = config->GetMotion_Origin_Y(jMarker) + VarCoord[1];
      Center[2] = config->GetMotion_Origin_Z(jMarker) + VarCoord[2];
      config->SetMotion_Origin_X(jMarker, Center[0]);
      config->SetMotion_Origin_Y(jMarker, Center[1]);
      config->SetMotion_Origin_Z(jMarker, Center[2]);
    }
  }
  
  /*--- Set the moment computation center to the new location after
   incrementing the position with the translation. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    Center[0] = config->GetRefOriginMoment_X(jMarker) + VarCoord[0];
    Center[1] = config->GetRefOriginMoment_Y(jMarker) + VarCoord[1];
    Center[2] = config->GetRefOriginMoment_Z(jMarker) + VarCoord[2];
    config->SetRefOriginMoment_X(jMarker, Center[0]);
    config->SetRefOriginMoment_Y(jMarker, Center[1]);
    config->SetRefOriginMoment_Z(jMarker, Center[2]);
  }
}

void CSurfaceMovement::Surface_Pitching(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
	
	double deltaT, time_new, time_old, Lref, *Coord;
  double Center[3], VarCoord[3], Omega[3], Ampl[3], Phase[3], rotCoord[3], r[3];
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
  double DEG2RAD = PI_NUMBER/180.0;
  unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  string Marker_Tag, Moving_Tag;
  int rank;
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
  
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<double>(iter-1)*deltaT;
  }

	/*--- Store displacement of each node on the pitching surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to pitch ---*/

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving(jMarker);
        Marker_Tag = config->GetMarker_All_Tag(iMarker);
        
        if (Marker_Tag == Moving_Tag) {
          
          /*--- Pitching origin, frequency, and amplitude from config. ---*/
          
          Center[0] = config->GetMotion_Origin_X(jMarker);
          Center[1] = config->GetMotion_Origin_Y(jMarker);
          Center[2] = config->GetMotion_Origin_Z(jMarker);
          Omega[0]  = config->GetPitching_Omega_X(jMarker)/config->GetOmega_Ref();
          Omega[1]  = config->GetPitching_Omega_Y(jMarker)/config->GetOmega_Ref();
          Omega[2]  = config->GetPitching_Omega_Z(jMarker)/config->GetOmega_Ref();
          Ampl[0]   = config->GetPitching_Ampl_X(jMarker)*DEG2RAD;
          Ampl[1]   = config->GetPitching_Ampl_Y(jMarker)*DEG2RAD;
          Ampl[2]   = config->GetPitching_Ampl_Z(jMarker)*DEG2RAD;
          Phase[0]  = config->GetPitching_Phase_X(jMarker)*DEG2RAD;
          Phase[1]  = config->GetPitching_Phase_Y(jMarker)*DEG2RAD;
          Phase[2]  = config->GetPitching_Phase_Z(jMarker)*DEG2RAD;
          
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing pitching displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Pitching frequency: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
              cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
              cout << " Pitching amplitude about origin: (" << Ampl[0]/DEG2RAD;
              cout << ", " << Ampl[1]/DEG2RAD << ", " << Ampl[2]/DEG2RAD;
              cout << ") degrees."<< endl;
              cout << " Pitching phase lag about origin: (" << Phase[0]/DEG2RAD;
              cout << ", " << Phase[1]/DEG2RAD <<", "<< Phase[2]/DEG2RAD;
              cout << ") degrees."<< endl;
            }
          }
          
          /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
          
          dtheta = -Ampl[0]*(sin(Omega[0]*time_new + Phase[0])
                             - sin(Omega[0]*time_old + Phase[0]));
          dphi   = -Ampl[1]*(sin(Omega[1]*time_new + Phase[1])
                             - sin(Omega[1]*time_old + Phase[1]));
          dpsi   = -Ampl[2]*(sin(Omega[2]*time_new + Phase[2])
                             - sin(Omega[2]*time_old + Phase[2]));
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
          sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Index and coordinates of the current point ---*/
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Coord  = geometry->node[iPoint]->GetCoord();
            
            /*--- Calculate non-dim. position from rotation center ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
            if (nDim == 2) r[nDim] = 0.0;
            
            /*--- Compute transformed point coordinates ---*/
            
            rotCoord[0] = rotMatrix[0][0]*r[0]
            + rotMatrix[0][1]*r[1]
            + rotMatrix[0][2]*r[2] + Center[0];
            
            rotCoord[1] = rotMatrix[1][0]*r[0]
            + rotMatrix[1][1]*r[1]
            + rotMatrix[1][2]*r[2] + Center[1];
            
            rotCoord[2] = rotMatrix[2][0]*r[0]
            + rotMatrix[2][1]*r[1]
            + rotMatrix[2][2]*r[2] + Center[2];
            
            /*--- Calculate delta change in the x, y, & z directions ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              VarCoord[iDim] = (rotCoord[iDim]-Coord[iDim])/Lref;
            if (nDim == 2) VarCoord[nDim] = 0.0;
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  /*--- For pitching we don't update the motion origin and moment reference origin. ---*/
}

void CSurfaceMovement::Surface_Rotating(CGeometry *geometry, CConfig *config,
                                        unsigned long iter, unsigned short iZone) {
	
	double deltaT, time_new, time_old, Lref, *Coord;
  double Center[3], VarCoord[3], Omega[3], rotCoord[3], r[3], Center_Aux[3];
  double rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double dtheta, dphi, dpsi, cosTheta, sinTheta;
  double cosPhi, sinPhi, cosPsi, sinPsi;
  unsigned short iMarker, jMarker, Moving, iDim, nDim = geometry->GetnDim();
  unsigned long iPoint, iVertex;
  string Marker_Tag, Moving_Tag;
  int rank;
  
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
  /*--- Initialize the delta variation in coordinates ---*/
  VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
  
  /*--- Retrieve values from the config file ---*/
  
  deltaT = config->GetDelta_UnstTimeND();
  Lref   = config->GetLength_Ref();
  
  /*--- Compute delta time based on physical time step ---*/
  time_new = static_cast<double>(iter)*deltaT;
  if (iter == 0) {
    time_old = time_new;
  } else {
    time_old = static_cast<double>(iter-1)*deltaT;
  }
  
	/*--- Store displacement of each node on the rotating surface ---*/
    /*--- Loop over markers and find the particular marker(s) (surface) to rotate ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Moving = config->GetMarker_All_Moving(iMarker);
    if (Moving == YES) {
      for (jMarker = 0; jMarker<config->GetnMarker_Moving(); jMarker++) {
        
        Moving_Tag = config->GetMarker_Moving(jMarker);
        Marker_Tag = config->GetMarker_All_Tag(iMarker);
        
        if (Marker_Tag == Moving_Tag) {
          
          /*--- Rotation origin and angular velocity from config. ---*/
          
          Center[0] = config->GetMotion_Origin_X(jMarker);
          Center[1] = config->GetMotion_Origin_Y(jMarker);
          Center[2] = config->GetMotion_Origin_Z(jMarker);
          Omega[0]  = config->GetRotation_Rate_X(jMarker)/config->GetOmega_Ref();
          Omega[1]  = config->GetRotation_Rate_Y(jMarker)/config->GetOmega_Ref();
          Omega[2]  = config->GetRotation_Rate_Z(jMarker)/config->GetOmega_Ref();
          
          /*--- Print some information to the console. Be verbose at the first
           iteration only (mostly for debugging purposes). ---*/
          // Note that the MASTER_NODE might not contain all the markers being moved.

          if (rank == MASTER_NODE) {
            cout << " Storing rotating displacement for marker: ";
            cout << Marker_Tag << "." << endl;
            if (iter == 0) {
              cout << " Angular velocity: (" << Omega[0] << ", " << Omega[1];
              cout << ", " << Omega[2] << ") rad/s about origin: (" << Center[0];
              cout << ", " << Center[1] << ", " << Center[2] << ")." << endl;
            }
          }
          
          /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
          
          dtheta = Omega[0]*(time_new-time_old);
          dphi   = Omega[1]*(time_new-time_old);
          dpsi   = Omega[2]*(time_new-time_old);
          
          /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
          
          cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
          sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
          
          /*--- Compute the rotation matrix. Note that the implicit
           ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
          
          rotMatrix[0][0] = cosPhi*cosPsi;
          rotMatrix[1][0] = cosPhi*sinPsi;
          rotMatrix[2][0] = -sinPhi;
          
          rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
          rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
          rotMatrix[2][1] = sinTheta*cosPhi;
          
          rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
          rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
          rotMatrix[2][2] = cosTheta*cosPhi;
          
          for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            
            /*--- Index and coordinates of the current point ---*/
            
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            Coord  = geometry->node[iPoint]->GetCoord();
            
            /*--- Calculate non-dim. position from rotation center ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
            if (nDim == 2) r[nDim] = 0.0;
            
            /*--- Compute transformed point coordinates ---*/
            
            rotCoord[0] = rotMatrix[0][0]*r[0]
            + rotMatrix[0][1]*r[1]
            + rotMatrix[0][2]*r[2] + Center[0];
            
            rotCoord[1] = rotMatrix[1][0]*r[0]
            + rotMatrix[1][1]*r[1]
            + rotMatrix[1][2]*r[2] + Center[1];
            
            rotCoord[2] = rotMatrix[2][0]*r[0]
            + rotMatrix[2][1]*r[1]
            + rotMatrix[2][2]*r[2] + Center[2];
            
            /*--- Calculate delta change in the x, y, & z directions ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              VarCoord[iDim] = (rotCoord[iDim]-Coord[iDim])/Lref;
            if (nDim == 2) VarCoord[nDim] = 0.0;
            
            /*--- Set node displacement for volume deformation ---*/
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
            
          }
        }
      }
    }
  }
  
  /*--- When updating the origins it is assumed that all markers have the
   same rotation movement, because we use the last markers rotation matrix and center ---*/
  
  /*--- Set the mesh motion center to the new location after
   incrementing the position with the rotation. This new
   location will be used for subsequent mesh motion for the given marker.---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Moving(); jMarker++) {
    
    /*-- Check if we want to update the motion origin for the given marker ---*/
    
    if (config->GetMoveMotion_Origin(jMarker) == YES) {
      /*--- There is no movement in the first iteration ---*/
      if (iter !=0) {
        
        Center_Aux[0] = config->GetMotion_Origin_X(jMarker);
        Center_Aux[1] = config->GetMotion_Origin_Y(jMarker);
        Center_Aux[2] = config->GetMotion_Origin_Z(jMarker);
        
        /*--- Calculate non-dim. position from rotation center ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          r[iDim] = (Center_Aux[iDim]-Center[iDim])/Lref;
        if (nDim == 2) r[nDim] = 0.0;
        
        /*--- Compute transformed point coordinates ---*/
        
        rotCoord[0] = rotMatrix[0][0]*r[0]
        + rotMatrix[0][1]*r[1]
        + rotMatrix[0][2]*r[2] + Center[0];
        
        rotCoord[1] = rotMatrix[1][0]*r[0]
        + rotMatrix[1][1]*r[1]
        + rotMatrix[1][2]*r[2] + Center[1];
        
        rotCoord[2] = rotMatrix[2][0]*r[0]
        + rotMatrix[2][1]*r[1]
        + rotMatrix[2][2]*r[2] + Center[2];
        
        /*--- Calculate delta change in the x, y, & z directions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          VarCoord[iDim] = (rotCoord[iDim]-Center_Aux[iDim])/Lref;
        if (nDim == 2) VarCoord[nDim] = 0.0;
        config->SetMotion_Origin_X(jMarker, Center_Aux[0]+VarCoord[0]);
        config->SetMotion_Origin_Y(jMarker, Center_Aux[1]+VarCoord[1]);
        config->SetMotion_Origin_Z(jMarker, Center_Aux[2]+VarCoord[2]);
      }
    }
  }

  /*--- Set the moment computation center to the new location after
   incrementing the position with the rotation. ---*/
  
  for (jMarker=0; jMarker<config->GetnMarker_Monitoring(); jMarker++) {
    /*--- There is no movement in the first iteration ---*/
    if (iter !=0) {
      
      Center_Aux[0] = config->GetRefOriginMoment_X(jMarker);
      Center_Aux[1] = config->GetRefOriginMoment_Y(jMarker);
      Center_Aux[2] = config->GetRefOriginMoment_Z(jMarker);

      /*--- Calculate non-dim. position from rotation center ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        r[iDim] = (Center_Aux[iDim]-Center[iDim])/Lref;
      if (nDim == 2) r[nDim] = 0.0;
      
      /*--- Compute transformed point coordinates ---*/
      
      rotCoord[0] = rotMatrix[0][0]*r[0]
      + rotMatrix[0][1]*r[1]
      + rotMatrix[0][2]*r[2] + Center[0];
      
      rotCoord[1] = rotMatrix[1][0]*r[0]
      + rotMatrix[1][1]*r[1]
      + rotMatrix[1][2]*r[2] + Center[1];
      
      rotCoord[2] = rotMatrix[2][0]*r[0]
      + rotMatrix[2][1]*r[1]
      + rotMatrix[2][2]*r[2] + Center[2];
      
      /*--- Calculate delta change in the x, y, & z directions ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        VarCoord[iDim] = (rotCoord[iDim]-Center_Aux[iDim])/Lref;
      if (nDim == 2) VarCoord[nDim] = 0.0;
      
      config->SetRefOriginMoment_X(jMarker, Center_Aux[0]+VarCoord[0]);
      config->SetRefOriginMoment_Y(jMarker, Center_Aux[1]+VarCoord[1]);
      config->SetRefOriginMoment_Z(jMarker, Center_Aux[2]+VarCoord[2]);
    }
  }
}

void CSurfaceMovement::AeroelasticDeform(CGeometry *geometry, CConfig *config, unsigned short iMarker, unsigned short iMarker_Monitoring, double displacements[4]) {
  /* The sign conventions of these are those of the Typical Section Wing Model, below the signs are corrected */
  double dy = -displacements[0];           // relative plunge
  double dalpha = -displacements[1];       // relative pitch
  double Center[2];
  unsigned short iDim;
  double Lref = config->GetLength_Ref();
  double *Coord;
  unsigned long iPoint, iVertex;
  double x_new, y_new;
  double VarCoord[3];
  string Marker_Tag;

  /*--- Pitching origin from config. ---*/
        
  Center[0] = config->GetRefOriginMoment_X(iMarker_Monitoring);
  Center[1] = config->GetRefOriginMoment_Y(iMarker_Monitoring);
  
  for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    /*--- Coordinates of the current point ---*/
    Coord = geometry->node[iPoint]->GetCoord();
    
    /*--- Calculate non-dim. position from rotation center ---*/
    double r[2] = {0,0};
    for (iDim = 0; iDim < geometry->GetnDim(); iDim++)
        r[iDim] = (Coord[iDim]-Center[iDim])/Lref;
    
    /*--- Compute delta of transformed point coordinates ---*/
    // The deltas are needed for the FEA grid deformation Method.
    // rotation contribution + plunging contribution - previous position
    x_new = cos(dalpha)*r[0] - sin(dalpha)*r[1] -r[0];
    y_new = sin(dalpha)*r[0] + cos(dalpha)*r[1] -r[1] + dy;
    
    VarCoord[0] = x_new;
    VarCoord[1] = y_new;
    VarCoord[2] = 0.0;
    
    /*--- Store new delta node locations for the surface ---*/
    geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
  }
  /*--- Set the elastic axis to the new location after incrementing the position with the plunge ---*/
  config->SetRefOriginMoment_Y(iMarker_Monitoring,Center[1]+dy);
  
}

void CSurfaceMovement::SetBoundary_Flutter3D(CGeometry *geometry, CConfig *config,
                                             CFreeFormDefBox **FFDBox, unsigned long iter, unsigned short iZone) {
	
	double omega, deltaT, *vel;
  double alpha, alpha_new, alpha_old;
  double time_new, time_old;
  double Center[3], Omega[3], Ampl[3], Phase[3];
  double DEG2RAD = PI_NUMBER/180.0;
  int rank;
  bool adjoint = config->GetAdjoint();
    
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#else
	rank = MASTER_NODE;
#endif
	
  /*--- Retrieve values from the config file ---*/
  deltaT = config->GetDelta_UnstTimeND();
  vel    = config->GetVelocity_FreeStreamND();
  
  /*--- Pitching origin, frequency, and amplitude from config. ---*/
  Center[0] = config->GetMotion_Origin_X(iZone);
  Center[1] = config->GetMotion_Origin_Y(iZone);
  Center[2] = config->GetMotion_Origin_Z(iZone);
  Omega[0]  = (config->GetPitching_Omega_X(iZone)/config->GetOmega_Ref());
  Omega[1]  = (config->GetPitching_Omega_Y(iZone)/config->GetOmega_Ref());
  Omega[2]  = (config->GetPitching_Omega_Z(iZone)/config->GetOmega_Ref());
  Ampl[0]   = config->GetPitching_Ampl_X(iZone)*DEG2RAD;
  Ampl[1]   = config->GetPitching_Ampl_Y(iZone)*DEG2RAD;
  Ampl[2]   = config->GetPitching_Ampl_Z(iZone)*DEG2RAD;
  Phase[0]   = config->GetPitching_Phase_X(iZone)*DEG2RAD;
  Phase[1]   = config->GetPitching_Phase_Y(iZone)*DEG2RAD;
  Phase[2]   = config->GetPitching_Phase_Z(iZone)*DEG2RAD;
  
  /*--- Compute delta time based on physical time step ---*/
  if (adjoint) {
    /*--- For the unsteady adjoint, we integrate backwards through
     physical time, so perform mesh motion in reverse. ---*/ 
    unsigned long nFlowIter  = config->GetnExtIter();
    unsigned long directIter = nFlowIter - iter - 1;
    time_new = static_cast<double>(directIter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(directIter)+1.0)*deltaT;
  } else {
    /*--- Forward time for the direct problem ---*/
    time_new = static_cast<double>(iter)*deltaT;
    time_old = time_new;
    if (iter != 0) time_old = (static_cast<double>(iter)-1.0)*deltaT;
  }
	
  /*--- Update the pitching angle at this time step. Flip sign for
   nose-up positive convention. ---*/
  omega     = Omega[2];
  alpha_new = Ampl[2]*sin(omega*time_new);
  alpha_old = Ampl[2]*sin(omega*time_old);
  alpha     = (1E-10 + (alpha_new - alpha_old))*(-PI_NUMBER/180.0);
	
	if (rank == MASTER_NODE)
		cout << "New dihedral angle (alpha): " << alpha_new/DEG2RAD << " degrees." << endl;
	
	unsigned short iOrder, jOrder, kOrder;
	short iFFDBox;
	double movement[3];
	bool *move = new bool [nFFDBox];
	unsigned short *index = new unsigned short[3];
	
	move[0] = true; move[1] = true; move[2] = true;	
  
	/*--- Change the value of the control point if move is true ---*/
	for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
		if (move[iFFDBox])
			for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
						index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
						double *coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						movement[0] = 0.0; movement[1] = 0.0; movement[2] = coord[1]*tan(alpha);
						FFDBox[iFFDBox]->SetControlPoints(index, movement);
					}
	
	/*--- Recompute cartesian coordinates using the new control points position ---*/
	for (iFFDBox = 0; iFFDBox < nFFDBox; iFFDBox++)
		SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox);
	
}

void CSurfaceMovement::SetExternal_Deformation(CGeometry *geometry, CConfig *config, unsigned short iZone, unsigned long iter) {
  
  int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
  
  /*--- Local variables ---*/
  
	unsigned short iDim, nDim; 
	unsigned long iPoint, flowIter = 0;
  unsigned long jPoint, GlobalIndex;
	double VarCoord[3], *Coord_Old = NULL, *Coord_New = NULL, Center[3];
  double Lref   = config->GetLength_Ref();
  double NewCoord[3], rotMatrix[3][3] = {{0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
  double r[3], rotCoord[3];
  unsigned long iVertex;
  unsigned short iMarker;
  char buffer[50];
  string motion_filename, UnstExt, text_line;
  ifstream motion_file;
  bool unsteady = config->GetUnsteady_Simulation();
  bool adjoint = config->GetAdjoint();
  
	/*--- Load stuff from config ---*/
  
	nDim = geometry->GetnDim();
  motion_filename = config->GetMotion_FileName();
  
  /*--- Set the extension for the correct unsteady mesh motion file ---*/
  
  if (unsteady) {
    if (adjoint) {
      /*--- For the unsteady adjoint, we integrate backwards through
       physical time, so perform mesh motion in reverse. ---*/
      unsigned long nFlowIter = config->GetnExtIter() - 1;
      flowIter  = nFlowIter - iter;
      unsigned short lastindex = motion_filename.find_last_of(".");
      motion_filename = motion_filename.substr(0, lastindex);
      if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
      if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
      if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
      if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
      if  (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
      UnstExt = string(buffer);
      motion_filename.append(UnstExt);
    } else {
      /*--- Forward time for the direct problem ---*/
      flowIter = iter;
      unsigned short lastindex = motion_filename.find_last_of(".");
      motion_filename = motion_filename.substr(0, lastindex);
      if ((int(flowIter) >= 0) && (int(flowIter) < 10)) sprintf (buffer, "_0000%d.dat", int(flowIter));
      if ((int(flowIter) >= 10) && (int(flowIter) < 100)) sprintf (buffer, "_000%d.dat", int(flowIter));
      if ((int(flowIter) >= 100) && (int(flowIter) < 1000)) sprintf (buffer, "_00%d.dat", int(flowIter));
      if ((int(flowIter) >= 1000) && (int(flowIter) < 10000)) sprintf (buffer, "_0%d.dat", int(flowIter));
      if  (int(flowIter) >= 10000) sprintf (buffer, "_%d.dat", int(flowIter));
      UnstExt = string(buffer);
      motion_filename.append(UnstExt);
    }
    
    if (rank == MASTER_NODE)
      cout << "Reading in the arbitrary mesh motion from direct iteration " << flowIter << "." << endl;
  }
  
  /*--- Open the motion file ---*/

  motion_file.open(motion_filename.data(), ios::in);
  /*--- Throw error if there is no file ---*/
  if (motion_file.fail()) {
    cout << "There is no mesh motion file!" << endl;
    exit(1);
  }
  
  /*--- Read in and store the new mesh node locations ---*/ 
  
  while (getline(motion_file,text_line)) {
    istringstream point_line(text_line);
    if (nDim == 2) point_line >> iPoint >> NewCoord[0] >> NewCoord[1];
    if (nDim == 3) point_line >> iPoint >> NewCoord[0] >> NewCoord[1] >> NewCoord[2];
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_Moving(iMarker) == YES) {
        for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          GlobalIndex = geometry->node[jPoint]->GetGlobalIndex();
          if (GlobalIndex == iPoint) {
            geometry->vertex[iMarker][iVertex]->SetVarCoord(NewCoord);
            break;
          }
        }
      }
    }
  }
  /*--- Close the restart file ---*/
  motion_file.close();
  
  /*--- If rotating as well, prepare the rotation matrix ---*/
  
  if (config->GetGrid_Movement() &&
      config->GetKind_GridMovement(iZone) == EXTERNAL_ROTATION) {
    
    /*--- Variables needed only for rotation ---*/
    
    double Omega[3], dt;
    double dtheta, dphi, dpsi, cosTheta, sinTheta;
    double cosPhi, sinPhi, cosPsi, sinPsi;
    
    /*--- Center of rotation & angular velocity vector from config ---*/
    Center[0] = config->GetMotion_Origin_X(iZone);
    Center[1] = config->GetMotion_Origin_Y(iZone);
    Center[2] = config->GetMotion_Origin_Z(iZone);
    
    /*--- Angular velocity vector from config ---*/
    
    dt = static_cast<double>(iter)*config->GetDelta_UnstTimeND();
    Omega[0]  = config->GetRotation_Rate_X(iZone);
    Omega[1]  = config->GetRotation_Rate_Y(iZone);
    Omega[2]  = config->GetRotation_Rate_Z(iZone);
    
    /*--- For the unsteady adjoint, use reverse time ---*/
    if (adjoint) {
      /*--- Set the first adjoint mesh position to the final direct one ---*/
      if (iter == 0) dt = ((double)config->GetnExtIter()-1) * dt;
      /*--- Reverse the rotation direction for the adjoint ---*/
      else dt = -1.0*dt;
    } else {
      /*--- No rotation at all for the first direct solution ---*/
      if (iter == 0) dt = 0;
    }
    
    /*--- Compute delta change in the angle about the x, y, & z axes. ---*/
    
    dtheta = Omega[0]*dt;   
    dphi   = Omega[1]*dt; 
    dpsi   = Omega[2]*dt;
    
    /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
    
    cosTheta = cos(dtheta);  cosPhi = cos(dphi);  cosPsi = cos(dpsi);
    sinTheta = sin(dtheta);  sinPhi = sin(dphi);  sinPsi = sin(dpsi);
    
    /*--- Compute the rotation matrix. Note that the implicit
     ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
    
    rotMatrix[0][0] = cosPhi*cosPsi;
    rotMatrix[1][0] = cosPhi*sinPsi;
    rotMatrix[2][0] = -sinPhi;
    
    rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
    rotMatrix[2][1] = sinTheta*cosPhi;
    
    rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
    rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
    rotMatrix[2][2] = cosTheta*cosPhi;
    
  }
  
  /*--- Loop through to find only moving surface markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Moving(iMarker) == YES) {
      
      /*--- Loop over all surface points for this marker ---*/
      
      for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        /*--- Get current and new coordinates from file ---*/
        
        Coord_Old = geometry->node[iPoint]->GetCoord();
        Coord_New = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        
        /*--- If we're also rotating, multiply each point by the
         rotation matrix. It is assumed that the coordinates in
         Coord_Old have already been rotated using SetRigid_Rotation(). ---*/
        
        if (config->GetGrid_Movement() &&
            config->GetKind_GridMovement(iZone) == EXTERNAL_ROTATION) {
          
          /*--- Calculate non-dim. position from rotation center ---*/
          
          for (iDim = 0; iDim < nDim; iDim++)
            r[iDim] = (Coord_New[iDim]-Center[iDim])/Lref;
          if (nDim == 2) r[nDim] = 0.0;
          
          /*--- Compute transformed point coordinates ---*/
          
          rotCoord[0] = rotMatrix[0][0]*r[0] 
                      + rotMatrix[0][1]*r[1] 
                      + rotMatrix[0][2]*r[2] + Center[0];
          
          rotCoord[1] = rotMatrix[1][0]*r[0] 
                      + rotMatrix[1][1]*r[1] 
                      + rotMatrix[1][2]*r[2] + Center[1];
          
          rotCoord[2] = rotMatrix[2][0]*r[0] 
                      + rotMatrix[2][1]*r[1] 
                      + rotMatrix[2][2]*r[2] + Center[2];
          
          /*--- Copy rotated coords back to original array for consistency ---*/
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_New[iDim] = rotCoord[iDim];
        }
        
        /*--- Calculate delta change in the x, y, & z directions ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          VarCoord[iDim] = (Coord_New[iDim]-Coord_Old[iDim])/Lref;
        if (nDim == 2) VarCoord[nDim] = 0.0;

        /*--- Set position changes to be applied by the spring analogy ---*/
        geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
        
      }
    }	
  }
}

void CSurfaceMovement::SetNACA_4Digits(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal, Ycurv, Yesp;

	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();	}

	double Ya = config->GetParamDV(0,0) / 100.0; /*--- Maximum camber as a fraction of the chord 
					(100 m is the first of the four digits) ---*/
	double Xa = config->GetParamDV(0,1) / 10.0; /*--- Location of maximum camber as a fraction of 
					the chord (10 p is the second digit in the NACA xxxx description) ---*/
	double t = config->GetParamDV(0,2) / 100.0; /*--- Maximum thickness as a fraction of the
					  chord (so 100 t gives the last two digits in 
					  the NACA 4-digit denomination) ---*/
		
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
				
				if  (Coord[0] < Xa) Ycurv = (2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow(Xa,2.0));
				else Ycurv = ((1.0-2.0*Xa)+2.0*Xa*Coord[0]-pow(Coord[0],2.0))*(Ya/pow((1.0-Xa), 2.0));
				
				Yesp = t*(1.4845*sqrt(Coord[0])-0.6300*Coord[0]-1.7580*pow(Coord[0],2.0)+
						  1.4215*pow(Coord[0],3.0)-0.518*pow(Coord[0],4.0));
				
				if (Normal[1] > 0) VarCoord[1] =  (Ycurv + Yesp) - Coord[1];
				if (Normal[1] < 0) VarCoord[1] =  (Ycurv - Yesp) - Coord[1];

			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetParabolic(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, *Normal;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();	}
	
	double c = config->GetParamDV(0,0); /*--- Center of the parabola ---*/
	double t = config->GetParamDV(0,1) / 100.0; /*--- Thickness of the parabola ---*/
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				Normal = boundary->vertex[iMarker][iVertex]->GetNormal();
				
				if (Normal[1] > 0) {
					VarCoord[1] =  t*(Coord[0]*Coord[0]-Coord[0])/(2.0*(c*c-c)) - Coord[1];
				}
				if (Normal[1] < 0) {
					VarCoord[1] =  t*(Coord[0]-Coord[0]*Coord[0])/(2.0*(c*c-c)) - Coord[1];
				}
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetObstacle(CGeometry *boundary, CConfig *config) {
	unsigned long iVertex, Point;
	unsigned short iMarker;
	double VarCoord[3], *Coord, xCoord;
	
	if (config->GetnDV() != 1) { cout << "This kind of design variable is not prepared for multiple deformations."; cin.get();	}
	
	double H = config->GetParamDV(0,0); /*--- Non-dimensionalized height of the obstacle ---*/
	double L = config->GetParamDV(0,1); /*--- Non-dimensionalized length of the obstacle ---*/
	double xOffSet = 0.0; /*--- x offset ---*/

	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
			VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
			if (config->GetMarker_All_DV(iMarker) == YES) {
				Point = boundary->vertex[iMarker][iVertex]->GetNode();
				Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
				xCoord = Coord[0]-xOffSet;
				if ((xCoord > 0) && (xCoord < L))
					VarCoord[1] = (27.0/4.0)*(H/(L*L*L))*xCoord*(xCoord-L)*(xCoord-L);
				else 
					VarCoord[1] = 0.0;
			}
			boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
		}
}

void CSurfaceMovement::SetAirfoil(CGeometry *boundary, CConfig *config) {
  unsigned long iVertex, Point;
  unsigned short iMarker;
  double VarCoord[3], *Coord, NewYCoord, NewXCoord, *Coord_i, *Coord_ip1;
  double yp1, ypn;
  unsigned short iVar;
  unsigned long n_Airfoil = 0;
  double Airfoil_Coord[2];
  vector<double> Svalue, Xcoord, Ycoord, Xcoord2, Ycoord2, Xcoord_Aux, Ycoord_Aux;
  bool AddBegin = true, AddEnd = true;
  double x_i, x_ip1, y_i, y_ip1;
  char AirfoilFile[256];
  char AirfoilFormat[15];
  char MeshOrientation[15];
  char AirfoilClose[15];
  double AirfoilScale;
  double TrailingEdge = 0.95;
  unsigned short nUpper, nLower, iUpper, iLower;
  ifstream airfoil_file;
  string text_line;
  
  /*--- Get the SU2 module. SU2_CFD will use this routine for dynamically
   deforming meshes (MARKER_MOVING), while SU2_MDC will use it for deforming
   meshes after imposing design variable surface deformations (DV_MARKER). ---*/
  
  unsigned short Kind_SU2 = config->GetKind_SU2();
  
  /*--- Read the coordinates. Two main formats:
   - Selig are in an x,y format starting from trailing edge, along the upper surface to the leading
   edge and back around the lower surface to trailing edge.
   - Lednicer are upper surface points leading edge to trailing edge and then lower surface leading
   edge to trailing edge.
   ---*/
  
  /*--- Open the restart file, throw an error if this fails. ---*/
  
  cout << "Enter the name of file with the airfoil information: ";
  scanf ("%s", AirfoilFile);
  airfoil_file.open(AirfoilFile, ios::in);
  if (airfoil_file.fail()) {
    cout << "There is no airfoil file!! "<< endl;
    exit(1);
  }
  cout << "Enter the format of the airfoil (Selig or Lednicer): ";
  scanf ("%s", AirfoilFormat);

  cout << "Thickness scaling (1.0 means no scaling)?: ";
  scanf ("%lf", &AirfoilScale);
  
  cout << "Close the airfoil (Yes or No)?: ";
  scanf ("%s", AirfoilClose);
  
  cout << "Surface mesh orientation (clockwise, or anticlockwise): ";
  scanf ("%s", MeshOrientation);
  
  /*--- The first line is the header ---*/
  
  getline (airfoil_file, text_line);
  cout << "File info: " << text_line << endl;
  
  if (strcmp (AirfoilFormat,"Selig") == 0) {

    while (getline (airfoil_file, text_line)) {
      istringstream point_line(text_line);
      
      /*--- Read the x & y coordinates from this line of the file (anticlockwise) ---*/
      
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      
      /*--- Store the coordinates in vectors ---*/
      
      Xcoord.push_back(Airfoil_Coord[0]);
      Ycoord.push_back(Airfoil_Coord[1]*AirfoilScale);
    }
    
  }
  if (strcmp (AirfoilFormat,"Lednicer") == 0) {
    
    /*--- The second line is the number of points ---*/

    getline(airfoil_file, text_line);
    istringstream point_line(text_line);
    double Upper, Lower;
    point_line >> Upper >> Lower;
    
    nUpper = int(Upper);
    nLower = int(Lower);
  
    Xcoord.resize(nUpper+nLower-1);
    Ycoord.resize(nUpper+nLower-1);
    
    /*--- White line ---*/

    getline (airfoil_file, text_line);

    for (iUpper = 0; iUpper < nUpper; iUpper++) {
      getline (airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      Xcoord[nUpper-iUpper-1] = Airfoil_Coord[0];
      
      double factor;
      if (strcmp (AirfoilClose,"Yes") == 0) {
        double x = (Airfoil_Coord[0] - TrailingEdge) / (1.0-TrailingEdge);
        factor = (1.0-x)+sin(PI_NUMBER*(1.0-x))/PI_NUMBER;
        if (x < TrailingEdge) factor = 1.0;
      }
      else factor = 1.0;
      
      Ycoord[nUpper-iUpper-1] = Airfoil_Coord[1]*AirfoilScale*factor;
    }
    
    getline (airfoil_file, text_line);

    for (iLower = 0; iLower < nLower; iLower++) {
      getline (airfoil_file, text_line);
      istringstream point_line(text_line);
      point_line >> Airfoil_Coord[0] >> Airfoil_Coord[1];
      
      double factor;
      if (strcmp (AirfoilClose,"Yes") == 0) {
        double x = (Airfoil_Coord[0] - TrailingEdge) / (1.0-TrailingEdge);
        factor = (1.0-x)+sin(PI_NUMBER*(1.0-x))/PI_NUMBER;
        if (x < TrailingEdge) factor = 1.0;
      }
      else factor = 1.0;
      
      Xcoord[nUpper+iLower-1] = Airfoil_Coord[0];
      Ycoord[nUpper+iLower-1] = Airfoil_Coord[1]*AirfoilScale*factor;
    }
      
  }
  
  /*--- Check the coordinate (1,0) at the beginning and end of the file ---*/
  
  if (Xcoord[0] == 1.0) AddBegin = false;
  if (Xcoord[Xcoord.size()-1] == 1.0) AddEnd = false;
  
  if (AddBegin) { Xcoord.insert(Xcoord.begin(), 1.0);   Ycoord.insert(Ycoord.begin(), 0.0);}
  if (AddEnd) { Xcoord.push_back(1.0);                Ycoord.push_back(0.0);}
  
  /*--- Change the orientation (depend on the input file, and the mesh file) ---*/
  
  if (strcmp (MeshOrientation,"clockwise") == 0) {
    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord_Aux.push_back(Xcoord[iVar]);
      Ycoord_Aux.push_back(Ycoord[iVar]);
    }
    
    for (iVar = 0; iVar < Xcoord.size(); iVar++) {
      Xcoord[iVar] = Xcoord_Aux[Xcoord.size()-iVar-1];
      Ycoord[iVar] = Ycoord_Aux[Xcoord.size()-iVar-1];
    }
  }
  
  /*--- Compute the total arch length ---*/
  
  double Arch = 0.0;
  Svalue.push_back(Arch);

  for (iVar = 0; iVar < Xcoord.size()-1; iVar++) {
    x_i = Xcoord[iVar];  x_ip1 = Xcoord[iVar+1];
    y_i = Ycoord[iVar];  y_ip1 = Ycoord[iVar+1];
    Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
    Svalue.push_back(Arch);
  }
  x_i = Xcoord[Xcoord.size()-1];  x_ip1 = Xcoord[0];
  y_i = Ycoord[Xcoord.size()-1];  y_ip1 = Ycoord[0];
  Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
  
  /*--- Non dimensionalization ---*/
  
  for (iVar = 0; iVar < Svalue.size(); iVar++) { Svalue[iVar] /= Arch; }

  /*--- Close the restart file ---*/
  
  airfoil_file.close();
  
  /*--- Create a spline for X and Y coordiantes using the arch length ---*/
  
  n_Airfoil = Svalue.size();
  yp1 = (Xcoord[1]-Xcoord[0])/(Svalue[1]-Svalue[0]);
  ypn = (Xcoord[n_Airfoil-1]-Xcoord[n_Airfoil-2])/(Svalue[n_Airfoil-1]-Svalue[n_Airfoil-2]);
  
  Xcoord2.resize(n_Airfoil+1);
  boundary->SetSpline(Svalue, Xcoord, n_Airfoil, yp1, ypn, Xcoord2);
  
  n_Airfoil = Svalue.size();
  yp1 = (Ycoord[1]-Ycoord[0])/(Svalue[1]-Svalue[0]);
  ypn = (Ycoord[n_Airfoil-1]-Ycoord[n_Airfoil-2])/(Svalue[n_Airfoil-1]-Svalue[n_Airfoil-2]);
  
  Ycoord2.resize(n_Airfoil+1);
  boundary->SetSpline(Svalue, Ycoord, n_Airfoil, yp1, ypn, Ycoord2);
  
  NewXCoord = 0.0; NewYCoord = 0.0;
  
  double TotalArch = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
        ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_MDC))) {
      for (iVertex = 0; iVertex < boundary->nVertex[iMarker]-1; iVertex++) {
        Coord_i = boundary->vertex[iMarker][iVertex]->GetCoord();
        Coord_ip1 = boundary->vertex[iMarker][iVertex+1]->GetCoord();
        
        x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
        y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
        
        TotalArch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
      }
      Coord_i = boundary->vertex[iMarker][boundary->nVertex[iMarker]-1]->GetCoord();
      Coord_ip1 = boundary->vertex[iMarker][0]->GetCoord();
      x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
      y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
      TotalArch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i));
    }
  }
  
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Arch = 0.0;
    for (iVertex = 0; iVertex < boundary->nVertex[iMarker]; iVertex++) {
      VarCoord[0] = 0.0; VarCoord[1] = 0.0; VarCoord[2] = 0.0;
      if (((config->GetMarker_All_Moving(iMarker) == YES) && (Kind_SU2 == SU2_CFD)) ||
          ((config->GetMarker_All_DV(iMarker) == YES) && (Kind_SU2 == SU2_MDC))) {
        Point = boundary->vertex[iMarker][iVertex]->GetNode();
        Coord = boundary->vertex[iMarker][iVertex]->GetCoord();
        
        if (iVertex == 0) Arch = 0.0;
        else {
          Coord_i = boundary->vertex[iMarker][iVertex-1]->GetCoord();
          Coord_ip1 = boundary->vertex[iMarker][iVertex]->GetCoord();
          x_i = Coord_i[0]; x_ip1 = Coord_ip1[0];
          y_i = Coord_i[1]; y_ip1 = Coord_ip1[1];
          Arch += sqrt((x_ip1-x_i)*(x_ip1-x_i)+(y_ip1-y_i)*(y_ip1-y_i))/TotalArch;
        }
        
        NewXCoord = boundary->GetSpline(Svalue, Xcoord, Xcoord2, n_Airfoil, Arch);
        NewYCoord = boundary->GetSpline(Svalue, Ycoord, Ycoord2, n_Airfoil, Arch);
        
        /*--- Store the delta change in the x & y coordinates ---*/
        
        VarCoord[0] = NewXCoord - Coord[0];
        VarCoord[1] = NewYCoord - Coord[1];
      }

      boundary->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);
      
    }
  }

}

void CSurfaceMovement::ReadFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox, string val_mesh_filename, bool val_fullmesh) {
	string text_line, iTag;
	ifstream mesh_file;
	double coord[3];
	unsigned short degree[3], iFFDBox, iCornerPoints, iControlPoints, iMarker, iDegree, jDegree, kDegree, iChar, LevelFFDBox, nParentFFDBox, iParentFFDBox, nChildFFDBox, iChildFFDBox, nMarker;
	unsigned long iSurfacePoints, iPoint, jPoint, iVertex, nVertex, nPoint, iElem = 0, nElem, my_nSurfPoints, nSurfPoints;

  int rank = MASTER_NODE;

#ifndef NO_MPI
#ifdef WINDOWS
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
	rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
	
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.open(cstr, ios::in);
	if (mesh_file.fail()) {
		cout << "There is no geometry file (ReadFFDInfo)!!" << endl;
		exit(1);
	}
	
	while (getline (mesh_file, text_line)) {
		
		/*--- Read the inner elements ---*/
		string::size_type position = text_line.find ("NELEM=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nElem = atoi(text_line.c_str());
			for (iElem = 0; iElem < nElem; iElem++)  {
				getline(mesh_file, text_line);
			}
		}
		
		/*--- Read the inner points ---*/
		position = text_line.find ("NPOINT=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nPoint = atoi(text_line.c_str());
			for (iPoint = 0; iPoint < nPoint; iPoint++)  {
				getline(mesh_file, text_line);
			}
		}

    /*--- Read the boundaries  ---*/
		position = text_line.find ("NMARK=",0);
		if (position != string::npos) {
			text_line.erase (0,6); nMarker = atoi(text_line.c_str());
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        getline(mesh_file, text_line);
        getline(mesh_file, text_line);
        text_line.erase (0,13); nVertex = atoi(text_line.c_str());
        for (iVertex = 0; iVertex < nVertex; iVertex++)  {
          getline(mesh_file, text_line);
        }
      }
		}
    
    /*--- Read the FFDBox information  ---*/
		position = text_line.find ("FFD_NBOX=",0);
		if (position != string::npos) {
			text_line.erase (0,9);
			nFFDBox = atoi(text_line.c_str());
			if (rank == MASTER_NODE) cout << nFFDBox << " Free Form Deformation (FFD) FFDBoxs." << endl;
			unsigned short *nCornerPoints = new unsigned short[nFFDBox];
			unsigned short *nControlPoints = new unsigned short[nFFDBox];
			unsigned long *nSurfacePoints = new unsigned long[nFFDBox];
			
			getline (mesh_file,text_line);
			text_line.erase (0,11);
			nLevel = atoi(text_line.c_str());
			if (rank == MASTER_NODE) cout << nLevel << " Free Form Deformation (FFD) nested levels." << endl;

			for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
				
				/*--- Read the name of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,8);
				
				/*--- Remove extra data from the FFDBox name ---*/
				string::size_type position;
				for (iChar = 0; iChar < 20; iChar++) {
					position = text_line.find( " ", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\r", 0 );
					if(position != string::npos) text_line.erase (position,1);
					position = text_line.find( "\n", 0 );
					if(position != string::npos) text_line.erase (position,1);
				}
				
				string TagFFDBox = text_line.c_str();
				if (rank == MASTER_NODE) cout << "FFD box tag: " << TagFFDBox <<". ";

				/*--- Read the level of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,10);
				LevelFFDBox = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "FFD box level: " << LevelFFDBox <<". ";
				
				/*--- Read the degree of the FFD box ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,13); degree[0] = atoi(text_line.c_str());
				getline (mesh_file,text_line);
				text_line.erase (0,13); degree[1] = atoi(text_line.c_str());
				getline (mesh_file,text_line);
				text_line.erase (0,13); degree[2] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Degrees: " << degree[0] <<", " << degree[1] <<", "<< degree[2] <<". "<< endl;
				FFDBox[iFFDBox] = new CFreeFormDefBox(int(degree[0]), int(degree[1]), int(degree[2]));				
				FFDBox[iFFDBox]->SetTag(TagFFDBox); FFDBox[iFFDBox]->SetLevel(LevelFFDBox);

				/*--- Read the number of parents boxes ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,12);
				nParentFFDBox = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Number of parent boxes: " << nParentFFDBox <<". ";
				for (iParentFFDBox = 0; iParentFFDBox < nParentFFDBox; iParentFFDBox++) {
					getline(mesh_file, text_line);
					
					/*--- Remove extra data from the FFDBox name ---*/
					string::size_type position;
					for (iChar = 0; iChar < 20; iChar++) {
						position = text_line.find( " ", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\r", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\n", 0 );
						if(position != string::npos) text_line.erase (position,1);
					}
					
					string ParentFFDBox = text_line.c_str();
					FFDBox[iFFDBox]->SetParentFFDBox(ParentFFDBox);
				}
				
				/*--- Read the number of children boxes ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,13);
				nChildFFDBox = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Number of child boxes: " << nChildFFDBox <<"." << endl;
				for (iChildFFDBox = 0; iChildFFDBox < nChildFFDBox; iChildFFDBox++) {
					getline(mesh_file, text_line);
					
					/*--- Remove extra data from the FFDBox name ---*/
					string::size_type position;
					for (iChar = 0; iChar < 20; iChar++) {
						position = text_line.find( " ", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\r", 0 );
						if(position != string::npos) text_line.erase (position,1);
						position = text_line.find( "\n", 0 );
						if(position != string::npos) text_line.erase (position,1);
					}
					
					string ChildFFDBox = text_line.c_str();
					FFDBox[iFFDBox]->SetChildFFDBox(ChildFFDBox);
				}
								
				/*--- Read the number of the corner points ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,18); nCornerPoints[iFFDBox] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Corner points: " << nCornerPoints[iFFDBox] <<". ";
				
				/*--- Read the coordinates of the corner points ---*/
				for (iCornerPoints = 0; iCornerPoints < nCornerPoints[iFFDBox]; iCornerPoints++) {
					getline(mesh_file,text_line); istringstream FFDBox_line(text_line);
					FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2];
					FFDBox[iFFDBox]->SetCoordCornerPoints(coord, iCornerPoints);
				}
				
				/*--- Read the number of the control points ---*/
				getline (mesh_file,text_line);
				text_line.erase (0,19); nControlPoints[iFFDBox] = atoi(text_line.c_str());
				if (rank == MASTER_NODE) cout << "Control points: " << nControlPoints[iFFDBox] <<". ";
				
				/*--- Method to identify if there is a FFDBox definition ---*/
				if (nControlPoints[iFFDBox] != 0) FFDBoxDefinition = true;

				/*--- Read the coordinates of the control points ---*/
				for (iControlPoints = 0; iControlPoints < nControlPoints[iFFDBox]; iControlPoints++) {
					getline(mesh_file,text_line); istringstream FFDBox_line(text_line);
					FFDBox_line >> iDegree; FFDBox_line >> jDegree; FFDBox_line >> kDegree; 
					FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2]; 
					FFDBox[iFFDBox]->SetCoordControlPoints(coord, iDegree, jDegree, kDegree); 
				}
				
				getline (mesh_file,text_line);
				text_line.erase (0,19); nSurfacePoints[iFFDBox] = atoi(text_line.c_str());

				/*--- The the surface points parametric coordinates ---*/
        my_nSurfPoints = 0;
				for (iSurfacePoints = 0; iSurfacePoints < nSurfacePoints[iFFDBox]; iSurfacePoints++) {
					getline(mesh_file,text_line); istringstream FFDBox_line(text_line);
					FFDBox_line >> iTag; FFDBox_line >> iPoint;
					iMarker = config->GetTag_Marker_All(iTag);
					FFDBox_line >> coord[0]; FFDBox_line >> coord[1]; FFDBox_line >> coord[2];
          
          if (val_fullmesh) {  // With vertices information (mesh deformation).
            for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
              jPoint =  geometry->vertex[iMarker][iVertex]->GetNode();
              if (iPoint == jPoint) {
                FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
                FFDBox[iFFDBox]->Set_VertexIndex(iVertex);
                FFDBox[iFFDBox]->Set_PointIndex(iPoint);
                FFDBox[iFFDBox]->Set_ParametricCoord(coord);
                FFDBox[iFFDBox]->Set_CartesianCoord(geometry->node[iPoint]->GetCoord());
                my_nSurfPoints++;
              }
            }
            /*--- It is possible to remove some points in the FFD that are
             not associated with surface vertices, this is the case of send receive
             points that are on the surface, but the surface is not in the domain ---*/
					}
          else {  // Without vertices information (partitioning).
            FFDBox[iFFDBox]->Set_MarkerIndex(iMarker);
            FFDBox[iFFDBox]->Set_PointIndex(iPoint);
            FFDBox[iFFDBox]->Set_ParametricCoord(coord);
            my_nSurfPoints = nSurfacePoints[iFFDBox];
          }
				}
        
        nSurfacePoints[iFFDBox] = my_nSurfPoints;
        nSurfPoints = 0;
#ifndef NO_MPI
        if (config->GetKind_SU2() != SU2_DDC)
#ifdef WINDOWS
          MPI_Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
          MPI::COMM_WORLD.Allreduce(&my_nSurfPoints, &nSurfPoints, 1, MPI::UNSIGNED_LONG, MPI::SUM);
#endif
        else
          nSurfPoints = my_nSurfPoints;
#else
				nSurfPoints = my_nSurfPoints;
#endif
				
				if (rank == MASTER_NODE) cout << "Surface points: " << nSurfPoints <<"."<<endl;
        
			}
			
			delete [] nCornerPoints;
			delete [] nControlPoints;
			delete [] nSurfacePoints;		
		}
	}
	mesh_file.close();
  
	if (nFFDBox == 0) {
		if (rank == MASTER_NODE) cout <<"There is no FFD box definition. Just in case, review the .su2 file" << endl;
	}

}

void CSurfaceMovement::WriteFFDInfo(CGeometry *geometry, CConfig *config, string val_mesh_filename) {
	ofstream mesh_file;
	unsigned short iOrder, jOrder, kOrder, iFFDBox, iCornerPoints, iMarker, iParentFFDBox, iChildFFDBox;
	unsigned long iVertex, iPoint, iSurfacePoints;
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.precision(15);
	mesh_file.open(cstr, ios::out | ios::app);
	
	mesh_file << "FFD_NBOX= " << nFFDBox << endl;
	mesh_file << "FFD_NLEVEL= " << nLevel << endl;
	
	for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
		
		mesh_file << "FFD_TAG= " << FFDBox[iFFDBox]->GetTag() << endl;
		mesh_file << "FFD_LEVEL= " << FFDBox[iFFDBox]->GetLevel() << endl;

		mesh_file << "FFD_DEGREE_I= " << FFDBox[iFFDBox]->GetlOrder()-1 << endl;
		mesh_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetmOrder()-1 << endl;
		mesh_file << "FFD_DEGREE_K= " << FFDBox[iFFDBox]->GetnOrder()-1 << endl;
		
		mesh_file << "FFD_PARENTS= " << FFDBox[iFFDBox]->GetnParentFFDBox() << endl;
		for (iParentFFDBox = 0; iParentFFDBox < FFDBox[iFFDBox]->GetnParentFFDBox(); iParentFFDBox++)
			mesh_file << FFDBox[iFFDBox]->GetParentFFDBoxTag(iParentFFDBox) << endl;
		mesh_file << "FFD_CHILDREN= " << FFDBox[iFFDBox]->GetnChildFFDBox() << endl;
		for (iChildFFDBox = 0; iChildFFDBox < FFDBox[iFFDBox]->GetnChildFFDBox(); iChildFFDBox++)
			mesh_file << FFDBox[iFFDBox]->GetChildFFDBoxTag(iChildFFDBox) << endl;
		
		mesh_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() << endl;
		for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints(); iCornerPoints++) {
			double *coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
			mesh_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
		}

		/*--- No FFD definition ---*/
		if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
			mesh_file << "FFD_CONTROL_POINTS= 0" << endl;
			mesh_file << "FFD_SURFACE_POINTS= 0" << endl;				
		}
		else {
			mesh_file << "FFD_CONTROL_POINTS= " << FFDBox[iFFDBox]->GetnControlPoints() << endl;
			for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
						double *coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						mesh_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
					}
      mesh_file << "FFD_SURFACE_POINTS= " << FFDBox[iFFDBox]->GetnSurfacePoint() << endl;
      for (iSurfacePoints = 0; iSurfacePoints < FFDBox[iFFDBox]->GetnSurfacePoint(); iSurfacePoints++) {
        iMarker = FFDBox[iFFDBox]->Get_MarkerIndex(iSurfacePoints);
        iVertex = FFDBox[iFFDBox]->Get_VertexIndex(iSurfacePoints);
        iPoint = FFDBox[iFFDBox]->Get_PointIndex(iSurfacePoints);
        double *parcoord = FFDBox[iFFDBox]->Get_ParametricCoord(iSurfacePoints);
        mesh_file << scientific << config->GetMarker_All_Tag(iMarker) << "\t" << iPoint << "\t" << parcoord[0] << "\t" << parcoord[1] << "\t" << parcoord[2] << endl;
      }
		}
	}
	mesh_file.close();
}

void CSurfaceMovement::WriteFFDInfo(CGeometry *geometry, CConfig *config, CFreeFormDefBox **FFDBox, string val_mesh_filename) {
	ofstream mesh_file;
	unsigned short iOrder, jOrder, kOrder, iFFDBox, iCornerPoints, iMarker, iParentFFDBox, iChildFFDBox;
	unsigned long iPoint, iSurfacePoints;
	char *cstr = new char [val_mesh_filename.size()+1];
	strcpy (cstr, val_mesh_filename.c_str());
	
	mesh_file.precision(15);
	mesh_file.open(cstr, ios::out | ios::app);
	
	mesh_file << "FFD_NBOX= " << nFFDBox << endl;
	mesh_file << "FFD_NLEVEL= " << nLevel << endl;
	
	for (iFFDBox = 0 ; iFFDBox < nFFDBox; iFFDBox++) {
		
		mesh_file << "FFD_TAG= " << FFDBox[iFFDBox]->GetTag() << endl;
		mesh_file << "FFD_LEVEL= " << FFDBox[iFFDBox]->GetLevel() << endl;
    
		mesh_file << "FFD_DEGREE_I= " << FFDBox[iFFDBox]->GetlOrder()-1 << endl;
		mesh_file << "FFD_DEGREE_J= " << FFDBox[iFFDBox]->GetmOrder()-1 << endl;
		mesh_file << "FFD_DEGREE_K= " << FFDBox[iFFDBox]->GetnOrder()-1 << endl;
		
		mesh_file << "FFD_PARENTS= " << FFDBox[iFFDBox]->GetnParentFFDBox() << endl;
		for (iParentFFDBox = 0; iParentFFDBox < FFDBox[iFFDBox]->GetnParentFFDBox(); iParentFFDBox++)
			mesh_file << FFDBox[iFFDBox]->GetParentFFDBoxTag(iParentFFDBox) << endl;
		mesh_file << "FFD_CHILDREN= " << FFDBox[iFFDBox]->GetnChildFFDBox() << endl;
		for (iChildFFDBox = 0; iChildFFDBox < FFDBox[iFFDBox]->GetnChildFFDBox(); iChildFFDBox++)
			mesh_file << FFDBox[iFFDBox]->GetChildFFDBoxTag(iChildFFDBox) << endl;
		
		mesh_file << "FFD_CORNER_POINTS= " << FFDBox[iFFDBox]->GetnCornerPoints() << endl;
		for (iCornerPoints = 0; iCornerPoints < FFDBox[iFFDBox]->GetnCornerPoints(); iCornerPoints++) {
			double *coord = FFDBox[iFFDBox]->GetCoordCornerPoints(iCornerPoints);
			mesh_file << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
		}
    
		/*--- No FFD definition ---*/
		if (FFDBox[iFFDBox]->GetnControlPoints() == 0) {
			mesh_file << "FFD_CONTROL_POINTS= 0" << endl;
			mesh_file << "FFD_SURFACE_POINTS= 0" << endl;
		}
		else {
			mesh_file << "FFD_CONTROL_POINTS= " << FFDBox[iFFDBox]->GetnControlPoints() << endl;
			for (iOrder = 0; iOrder < FFDBox[iFFDBox]->GetlOrder(); iOrder++)
				for (jOrder = 0; jOrder < FFDBox[iFFDBox]->GetmOrder(); jOrder++)
					for (kOrder = 0; kOrder < FFDBox[iFFDBox]->GetnOrder(); kOrder++) {
						double *coord = FFDBox[iFFDBox]->GetCoordControlPoints(iOrder, jOrder, kOrder);
						mesh_file << iOrder << "\t" << jOrder << "\t" << kOrder << "\t" << coord[0] << "\t" << coord[1] << "\t" << coord[2] << endl;
					}

      /*--- Compute the number of points on the new surfaces, note that we are not
       adding the new ghost points (receive), which eventually are also inside the chunck ---*/
      unsigned long nSurfacePoint = 0;
      for (iSurfacePoints = 0; iSurfacePoints < FFDBox[iFFDBox]->GetnSurfacePoint(); iSurfacePoints++) {
        iPoint = FFDBox[iFFDBox]->Get_PointIndex(iSurfacePoints);
        
        /*--- Note that we cannot use && because of the range of 
         geometry->GetGlobal_to_Local_Point(iPoint) is smaller than iPoint ---*/
        
        if (iPoint <= geometry->GetMax_GlobalPoint()) {
          if (geometry->GetGlobal_to_Local_Point(iPoint) != -1) nSurfacePoint++;
        }
      }
      
      mesh_file << "FFD_SURFACE_POINTS= " << nSurfacePoint << endl;
      for (iSurfacePoints = 0; iSurfacePoints < FFDBox[iFFDBox]->GetnSurfacePoint(); iSurfacePoints++) {
        iMarker = FFDBox[iFFDBox]->Get_MarkerIndex(iSurfacePoints);
        iPoint = FFDBox[iFFDBox]->Get_PointIndex(iSurfacePoints);
        if (iPoint <= geometry->GetMax_GlobalPoint()) {
          if (geometry->GetGlobal_to_Local_Point(iPoint) != -1) {
            double *parCoord = FFDBox[iFFDBox]->Get_ParametricCoord(iSurfacePoints);
            mesh_file << scientific << config->GetMarker_All_Tag(iMarker) << "\t" << geometry->GetGlobal_to_Local_Point(iPoint) << "\t" << parCoord[0] << "\t" << parCoord[1] << "\t" << parCoord[2] << endl;
          }
        }
      }
			
		}
		
	}
	mesh_file.close();
}

CFreeFormDefBox::CFreeFormDefBox(void) : CGridMovement() { }

CFreeFormDefBox::CFreeFormDefBox(unsigned short val_lDegree, unsigned short val_mDegree, unsigned short val_nDegree) : CGridMovement() {
	unsigned short iCornerPoints, iOrder, jOrder, kOrder, iDim;
	
	/*--- Only for 3D problems and FFD with Hexahedron ---*/
	nDim = 3;
	nCornerPoints = 8;
	
	/*--- Allocate Corners points ---*/
	Coord_Corner_Points = new double* [nCornerPoints];
	for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
		Coord_Corner_Points[iCornerPoints] = new double [nDim];
	
	param_coord = new double[nDim]; param_coord_ = new double[nDim];
	cart_coord = new double[nDim]; cart_coord_ = new double[nDim];
	gradient = new double[nDim];

	lDegree = val_lDegree; lOrder = lDegree+1;
	mDegree = val_mDegree; mOrder = mDegree+1;
	nDegree = val_nDegree; nOrder = nDegree+1;
	nControlPoints = lOrder*mOrder*nOrder;
	
	Coord_Control_Points = new double*** [lOrder];
	ParCoord_Control_Points = new double*** [lOrder];
	Coord_Control_Points_Copy = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		Coord_Control_Points[iOrder] = new double** [mOrder];
		ParCoord_Control_Points[iOrder] = new double** [mOrder];
		Coord_Control_Points_Copy[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			Coord_Control_Points[iOrder][jOrder] = new double* [nOrder];
			ParCoord_Control_Points[iOrder][jOrder] = new double* [nOrder];
			Coord_Control_Points_Copy[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				Coord_Control_Points[iOrder][jOrder][kOrder] = new double [nDim];
				ParCoord_Control_Points[iOrder][jOrder][kOrder] = new double [nDim];
				Coord_Control_Points_Copy[iOrder][jOrder][kOrder] = new double [nDim];
			}
		}
	}
	
	/*--- Zero-initialization ---*/
	for (iOrder = 0; iOrder < lOrder; iOrder++) 
		for (jOrder = 0; jOrder < mOrder; jOrder++) 
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				for (iDim = 0; iDim < nDim; iDim++)
					Coord_Control_Points[iOrder][jOrder][kOrder][iDim] = 0.0;
}

CFreeFormDefBox::~CFreeFormDefBox(void) {
	unsigned short iOrder, jOrder, kOrder, iCornerPoints;
	
	for (iOrder = 0; iOrder < lOrder; iOrder++) 
		for (jOrder = 0; jOrder < mOrder; jOrder++) 
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				delete [] Coord_Control_Points[iOrder][jOrder][kOrder];
				delete [] ParCoord_Control_Points[iOrder][jOrder][kOrder];
				delete [] Coord_Control_Points_Copy[iOrder][jOrder][kOrder];
			}
	delete [] Coord_Control_Points;
	delete [] ParCoord_Control_Points;
	delete [] Coord_Control_Points_Copy;

	delete [] param_coord;
	delete [] cart_coord;
	delete [] gradient;
	
	for (iCornerPoints = 0; iCornerPoints < nCornerPoints; iCornerPoints++)
		delete [] Coord_Corner_Points[iCornerPoints];
	delete [] Coord_Corner_Points;
}

void  CFreeFormDefBox::SetUnitCornerPoints(void) {
	double coord [3];
	
	coord [0] = 0.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,0);
	coord [0] = 1.0; coord [1] = 0.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,1);
	coord [0] = 1.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,2);
	coord [0] = 0.0; coord [1] = 1.0; coord [2] = 0.0; this->SetCoordCornerPoints(coord,3);
	coord [0] = 0.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,4);
	coord [0] = 1.0; coord [1] = 0.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,5);
	coord [0] = 1.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,6);
	coord [0] = 0.0; coord [1] = 1.0; coord [2] = 1.0; this->SetCoordCornerPoints(coord,7);
}

void CFreeFormDefBox::SetControlPoints_Parallelepiped (void) {
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	/*--- Set base control points according to the notation of Vtk for hexahedrons ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_Control_Points	[0]			[0]			[0]			[iDim]	= Coord_Corner_Points[0][iDim];
		Coord_Control_Points	[lOrder-1]	[0]			[0]			[iDim]	= Coord_Corner_Points[1][iDim];
		Coord_Control_Points	[lOrder-1]	[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[2][iDim];
		Coord_Control_Points	[0]			[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[3][iDim];
		Coord_Control_Points	[0]			[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[4][iDim];
		Coord_Control_Points	[lOrder-1]	[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[5][iDim];
		Coord_Control_Points	[lOrder-1]	[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[6][iDim];
		Coord_Control_Points	[0]			[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[7][iDim];
	}
	
	/*--- Fill the rest of the cubic matrix of control points with uniform spacing (parallelepiped) ---*/
	for (iDegree = 0; iDegree <= lDegree; iDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (kDegree = 0; kDegree <= nDegree; kDegree++) {
				Coord_Control_Points[iDegree][jDegree][kDegree][0] = Coord_Corner_Points[0][0] 
				+ double(iDegree)/double(lDegree)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
				Coord_Control_Points[iDegree][jDegree][kDegree][1] = Coord_Corner_Points[0][1] 
				+ double(jDegree)/double(mDegree)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);				
				Coord_Control_Points[iDegree][jDegree][kDegree][2] = Coord_Corner_Points[0][2] 
				+ double(kDegree)/double(nDegree)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
			}
}

void CFreeFormDefBox::SetSupportCP(CFreeFormDefBox *FFDBox) {
	unsigned short iDim, iOrder, jOrder, kOrder;
	unsigned short lOrder = FFDBox->GetlOrder();
	unsigned short mOrder = FFDBox->GetmOrder();
	unsigned short nOrder = FFDBox->GetnOrder();
	
	Coord_SupportCP = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		Coord_SupportCP[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			Coord_SupportCP[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				Coord_SupportCP[iOrder][jOrder][kOrder] = new double [nDim];
		}
	}
	
	/*--- Set base support control points according to the notation of Vtk for hexahedrons ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_SupportCP	[0]			[0]			[0]			[iDim]	= Coord_Corner_Points[0][iDim];
		Coord_SupportCP	[lOrder-1]	[0]			[0]			[iDim]	= Coord_Corner_Points[1][iDim];
		Coord_SupportCP	[lOrder-1]	[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[2][iDim];
		Coord_SupportCP	[0]			[mOrder-1]	[0]			[iDim]	= Coord_Corner_Points[3][iDim];
		Coord_SupportCP	[0]			[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[4][iDim];
		Coord_SupportCP	[lOrder-1]	[0]			[nOrder-1]	[iDim]	= Coord_Corner_Points[5][iDim];
		Coord_SupportCP	[lOrder-1]	[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[6][iDim];
		Coord_SupportCP	[0]			[mOrder-1]	[nOrder-1]	[iDim]	= Coord_Corner_Points[7][iDim];
	}
	
	/*--- Fill the rest of the cubic matrix of support control points with uniform spacing  ---*/
	for (iOrder = 0; iOrder < lOrder; iOrder++)
		for (jOrder = 0; jOrder < mOrder; jOrder++)
			for (kOrder = 0; kOrder < nOrder; kOrder++) {
				Coord_SupportCP[iOrder][jOrder][kOrder][0] = Coord_Corner_Points[0][0] 
				+ double(iOrder)/double(lOrder-1)*(Coord_Corner_Points[1][0]-Coord_Corner_Points[0][0]);
				Coord_SupportCP[iOrder][jOrder][kOrder][1] = Coord_Corner_Points[0][1] 
				+ double(jOrder)/double(mOrder-1)*(Coord_Corner_Points[3][1]-Coord_Corner_Points[0][1]);				
				Coord_SupportCP[iOrder][jOrder][kOrder][2] = Coord_Corner_Points[0][2] 
				+ double(kOrder)/double(nOrder-1)*(Coord_Corner_Points[4][2]-Coord_Corner_Points[0][2]);
			}
}

void CFreeFormDefBox::SetSupportCPChange(CFreeFormDefBox *FFDBox) {
	unsigned short iDim, iOrder, jOrder, kOrder;
	double movement[3], *car_coord_old, *car_coord_new, *par_coord;
	unsigned short lOrder = FFDBox->GetlOrder();
	unsigned short mOrder = FFDBox->GetmOrder();
	unsigned short nOrder = FFDBox->GetnOrder();
	unsigned short *index = new unsigned short[nDim];

	double ****param_Coord_SupportCP = new double*** [lOrder];
	for (iOrder = 0; iOrder < lOrder; iOrder++) {
		param_Coord_SupportCP[iOrder] = new double** [mOrder];
		for (jOrder = 0; jOrder < mOrder; jOrder++) {
			param_Coord_SupportCP[iOrder][jOrder] = new double* [nOrder];
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				param_Coord_SupportCP[iOrder][jOrder][kOrder] = new double [nDim];
		}
	}
	
	for (iOrder = 0; iOrder < lOrder; iOrder++)
		for (jOrder = 0; jOrder < mOrder; jOrder++)
			for (kOrder = 0; kOrder < nOrder; kOrder++)
				for (iDim = 0; iDim < nDim; iDim++)
					param_Coord_SupportCP[iOrder][jOrder][kOrder][iDim] = 
					Coord_SupportCP[iOrder][jOrder][kOrder][iDim];
	
	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_Control_Points[0][0][0][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 0);
		Coord_Control_Points[1][0][0][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 1);
		Coord_Control_Points[1][1][0][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 2);
		Coord_Control_Points[0][1][0][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 3);
		Coord_Control_Points[0][0][1][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 4);
		Coord_Control_Points[1][0][1][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 5);
		Coord_Control_Points[1][1][1][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 6);
		Coord_Control_Points[0][1][1][iDim]	= FFDBox->GetCoordCornerPoints(iDim, 7);
	}
	
	for (iOrder = 0; iOrder < FFDBox->GetlOrder(); iOrder++)
		for (jOrder = 0; jOrder < FFDBox->GetmOrder(); jOrder++)
			for (kOrder = 0; kOrder < FFDBox->GetnOrder(); kOrder++) {
				par_coord = param_Coord_SupportCP[iOrder][jOrder][kOrder];
				car_coord_new = EvalCartesianCoord(par_coord);
				car_coord_old = FFDBox->GetCoordControlPoints(iOrder, jOrder, kOrder);
				index[0] = iOrder; index[1] = jOrder; index[2] = kOrder;
				movement[0] = car_coord_new[0] - car_coord_old[0]; 
				movement[1] = car_coord_new[1] - car_coord_old[1]; 
				movement[2] = car_coord_new[2] - car_coord_old[2]; 
				FFDBox->SetControlPoints(index, movement);
			}
}

void CFreeFormDefBox::SetTecplot(char FFDBox_filename[200], bool new_file) {
	ofstream FFDBox_file;
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	if (new_file) {
		FFDBox_file.open(FFDBox_filename, ios::out);
		FFDBox_file << "TITLE = \"Visualization of the FFD box\"" << endl;
		FFDBox_file << "VARIABLES = \"x\", \"y\", \"z\"" << endl;
	}
	else FFDBox_file.open(FFDBox_filename, ios::out | ios::app);

	FFDBox_file << "ZONE I="<<lDegree+1<<", J="<<mDegree+1<<", K="<<nDegree+1<<", DATAPACKING=POINT" << endl;
	
	FFDBox_file.precision(15);
	
	for (kDegree = 0; kDegree <= nDegree; kDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (iDegree = 0; iDegree <= lDegree; iDegree++) {
				for(iDim = 0; iDim < nDim; iDim++)
					FFDBox_file << scientific << Coord_Control_Points[iDegree][jDegree][kDegree][iDim] << "\t";
				FFDBox_file << "\n";
			}
		
	FFDBox_file.close();
}


double *CFreeFormDefBox::GetParametricCoord_Analytical(double *cart_coord) {
	unsigned short iDim;
	double *e1, *e2, *e3, *e12, *e23, *e13, *p;
	
	/*--- Auxiliary Basis Vectors of the deformed FFDBox ---*/
	e1 = new double[3]; e2 = new double[3]; e3 = new double[3];
	for (iDim = 0; iDim < nDim; iDim++) {
		e1[iDim] = Coord_Corner_Points[1][iDim]-Coord_Corner_Points[0][iDim];
		e2[iDim] = Coord_Corner_Points[3][iDim]-Coord_Corner_Points[0][iDim];
		e3[iDim] = Coord_Corner_Points[4][iDim]-Coord_Corner_Points[0][iDim];
	}
	
	/*--- Respective Cross-Products ---*/
	e12 = new double[3]; e23 = new double[3]; e13 = new double[3];
	CrossProduct(e1,e2,e12);
	CrossProduct(e1,e3,e13);
	CrossProduct(e2,e3,e23);
	
	/*--- p is Tranlated vector from the origin ---*/
	p = new double[3];
	for (iDim = 0; iDim < nDim; iDim++)
		p[iDim] = cart_coord[iDim] - Coord_Corner_Points[0][iDim];
	
	param_coord[0] = DotProduct(e23,p)/DotProduct(e23,e1);
	param_coord[1] = DotProduct(e13,p)/DotProduct(e13,e2);
	param_coord[2] = DotProduct(e12,p)/DotProduct(e12,e3);
	
	delete [] e1;
  delete [] e2;
  delete [] e3;
  delete [] e12;
  delete [] e23;
  delete [] e13;
  delete [] p;
	
	return param_coord;
}

double *CFreeFormDefBox::EvalCartesianCoord(double *param_coord) {
	unsigned short iDim, iDegree, jDegree, kDegree;
	
	for (iDim = 0; iDim < nDim; iDim++)
		cart_coord[iDim] = 0;
	
	for (iDegree = 0; iDegree <= lDegree; iDegree++)
		for (jDegree = 0; jDegree <= mDegree; jDegree++)
			for (kDegree = 0; kDegree <= nDegree; kDegree++)
				for (iDim = 0; iDim < nDim; iDim++) {
					cart_coord[iDim] += Coord_Control_Points[iDegree][jDegree][kDegree][iDim]
					* GetBernstein(lDegree, iDegree, param_coord[0])
					* GetBernstein(mDegree, jDegree, param_coord[1])
					* GetBernstein(nDegree, kDegree, param_coord[2]);
				}
	
	return cart_coord;
}

double CFreeFormDefBox::GetBernstein(short val_n, short val_i, double val_t) {
	double value;

	if (val_i > val_n) { value = 0; return value; }
	if (val_i == 0) {
		if (val_t == 0) value = 1;
		else if (val_t == 1) value = 0;
		else value = Binomial(val_n,val_i)*(pow(val_t, val_i)) * pow(1.0 - val_t, val_n - val_i);
	}
	else if (val_i == val_n) {
		if (val_t == 0) value = 0;
		else if (val_t == 1) value = 1;
		else value = pow(val_t,val_n);
	}
	else value = Binomial(val_n,val_i)*(pow(val_t,val_i)) * pow(1.0-val_t, val_n - val_i);
	
	return value;
}

double CFreeFormDefBox::GetBernsteinDerivative(short val_n, short val_i, 
											   double val_t, short val_order) {
	double value = 0.0;
	
	/*--- Verify this subroutine, it provides negative val_n, 
	 which is a wrong value for GetBernstein ---*/
	
	if (val_order == 0) { 
		value = GetBernstein(val_n, val_i, val_t); return value; 
	}
	
	if (val_i == 0) { 
		value = val_n*(-GetBernsteinDerivative(val_n-1, val_i, val_t, val_order-1)); return value; 
	}
	else {
		if (val_n == 0) { 
			value = val_t; return value; 
		}
		else {
			value = val_n*(GetBernsteinDerivative(val_n-1, val_i-1, val_t, val_order-1) - GetBernsteinDerivative(val_n-1, val_i, val_t, val_order-1));
			return value;
		}
	}

	return value;
}

double *CFreeFormDefBox::GetGradient_Analytical(double *val_coord, double *xyz) {
	unsigned short iDim, jDim, lmn[3];
	
	/*--- Set the Degree of the Berstein polynomials ---*/
	lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
	
	for (iDim = 0; iDim < nDim; iDim++) gradient[iDim] = 0;
	
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			gradient[jDim] += GetDerivative2(val_coord, iDim, xyz,  lmn) *  
			GetDerivative3(val_coord, iDim, jDim, lmn);
	
	return gradient;
}

double *CFreeFormDefBox::GetGradient_Numerical(double *uvw, double *xyz) {
	double delta = 1E-6, parametric[3], *coord_eval, functional_plus, functional_minus;
	
	parametric[0] = uvw[0] + delta;
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0] - delta;
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[0] = 0.5*(functional_plus-functional_minus)/delta;
	
	parametric[0] = uvw[0];
	parametric[1] = uvw[1] + delta;
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0];
	parametric[1] = uvw[1] - delta;
	parametric[2] = uvw[2]; 
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[1] = 0.5*(functional_plus-functional_minus)/delta;
	
	parametric[0] = uvw[0];
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2] + delta;
	coord_eval = EvalCartesianCoord(parametric);
	functional_plus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
					   (coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
					   (coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	parametric[0] = uvw[0];
	parametric[1] = uvw[1]; 
	parametric[2] = uvw[2] - delta;
	coord_eval = EvalCartesianCoord(parametric);
	functional_minus = ((coord_eval[0]-xyz[0])*(coord_eval[0]-xyz[0]) + 
						(coord_eval[1]-xyz[1])*(coord_eval[1]-xyz[1]) +
						(coord_eval[2]-xyz[2])*(coord_eval[2]-xyz[2]));
	gradient[2] = 0.5*(functional_plus-functional_minus)/delta;
	
	return gradient;
}

double *CFreeFormDefBox::GetParametricCoord_Iterative(double *xyz, double *guess, double tol, 
																										 unsigned long it_max) {
	double **Hessian, Indep_Term[3], under_relax = 1.0, MinNormError, NormError;
	unsigned short iDim, RandonCounter;
	unsigned long iter;
	
	/*--- Allocate the Hessian ---*/
	Hessian = new double* [nDim];
	for (iDim = 0; iDim < nDim; iDim++) {
		Hessian[iDim] = new double[nDim];
		param_coord[iDim] = guess[iDim];
		Indep_Term [iDim] = 0.0;
	}
	
	RandonCounter = 0; MinNormError = 1E6;
	
	for (iter = 0; iter < it_max; iter++) {
		
		/*--- The independent term of the solution of our system is -Gradient(sol_old) ---*/
		gradient = GetGradient_Analytical(param_coord, xyz);
		
		for (iDim = 0; iDim < nDim; iDim++) 
			Indep_Term[iDim] = -gradient[iDim];
						
		/*--- Relaxation of the Newton Method ---*/
		for (iDim = 0; iDim < nDim; iDim++) 
			Indep_Term[iDim] = under_relax * Indep_Term[iDim];
		
		/*--- Hessian = The Matrix of our system, getHessian(sol_old,xyz,...) ---*/
		GetHessian_Analytical(param_coord, xyz, Hessian);
		
		/*--- Gauss elimination algorithm. Solution will be stored on Indep_Term ---*/
		Gauss_Elimination(Hessian, Indep_Term, nDim);				
		
		/*--- Solution is in fact par_new-par_old; Must Update doing par_new=par_old + solution ---*/
		for (iDim = 0; iDim < nDim; iDim++) 
			param_coord[iDim] += Indep_Term[iDim];
		
		/*--- If the gradient is small, we have converged ---*/
		if ((fabs(Indep_Term[0]) < tol) && (fabs(Indep_Term[1]) < tol) && (fabs(Indep_Term[2]) < tol))	break;
		NormError = sqrt(Indep_Term[0]*Indep_Term[0] + Indep_Term[1]*Indep_Term[1] + Indep_Term[2]*Indep_Term[2]);
		MinNormError = min(NormError, MinNormError);
		
		/*--- If we have no convergence with 50 iterations probably we are out of the FFDBox, then 
		 we try with a ramdom choice ---*/
		if (((iter % 50) == 0) && (iter != 0)) {
			RandonCounter++;
			param_coord[0] = double(rand())/double(RAND_MAX);
			param_coord[1] = double(rand())/double(RAND_MAX);
			param_coord[2] = double(rand())/double(RAND_MAX);
		}
		
		if (RandonCounter == 100) {
			cout << "I can not localize this point: " << xyz[0] <<" "<< xyz[1] <<" "<< xyz[2] <<". Min Error: "<< MinNormError <<"."<< endl;
			param_coord[0] = 0.0; param_coord[1] = 0.0; param_coord[2] = 0.0;
			Indep_Term[0] = 0.0; Indep_Term[1] = 0.0; Indep_Term[2] = 0.0;
			break;
		}
	}
	
	/*--- There is no convergence of the point inversion algorithm ---*/
	if ((fabs(Indep_Term[0]) > tol) || (fabs(Indep_Term[1]) > tol) || (fabs(Indep_Term[2]) > tol))
		cout << "No Convergence Detected After " << iter << " Iterations" << endl;
	
	for (iDim = 0; iDim < nDim; iDim++) 
		delete [] Hessian[iDim];
	delete [] Hessian;
	
	/*--- Real Solution is now par_coord; Return it ---*/
	return param_coord;
}

unsigned short CFreeFormDefBox::Binomial (unsigned short n, unsigned short m) {
	unsigned short result;

	if ( (m == 0) || (m == n) ) result = 1;
	else result = Factorial(n) / (Factorial(n-m)*Factorial(m));
		
	return result;
}

unsigned long CFreeFormDefBox::BinomialOpt (unsigned long n, unsigned long m) {
	unsigned long b[100], i , j;
	if (n+1 > 100) cout << "ERROR!!! Increase the size of b in the BinomialOpt subroutine!" <<endl;
	
	b[0] = 1;
	for (i = 1; i <= n; ++i) {
		b[i] = 1;
		for(j = i-1U; j > 0; --j)
			b[j] += b[j-1U];
	}
	
	return b[m];
}

unsigned short CFreeFormDefBox::Factorial (unsigned short n) {
	
	if ( n > 1 ) n = n*Factorial(n-1);
	if ( n == 0 ) n = 1;
	
	return n;
}


bool CFreeFormDefBox::GetPointFFD(CGeometry *geometry, CConfig *config, unsigned long iPoint) {
	double *Coord;
	unsigned short iVar, jVar;
	bool Inside;
	
	unsigned short Index[5][7] = {
		{0, 1, 2, 5, 0, 1, 2},
		{0, 2, 7, 5, 0, 2, 7},
		{0, 2, 3, 7, 0, 2, 3},
		{0, 5, 7, 4, 0, 5, 7},
		{2, 7, 5, 6, 2, 7, 5}};
	
	Coord = geometry->node[iPoint]->GetCoord();
	
	/*--- 1st tetrahedron {V0, V1, V2, V5}
	 2nd tetrahedron {V0, V2, V7, V5}
	 3th tetrahedron {V0, V2, V3, V7}
	 4th tetrahedron {V0, V5, V7, V4}
	 5th tetrahedron {V2, V7, V5, V6} ---*/
	
	for (iVar = 0; iVar < 5; iVar++) {
		Inside = true;
		for (jVar = 0; jVar < 4; jVar++) {
			double Distance_Point = geometry->Point2Plane_Distance(Coord, 
																														 Coord_Corner_Points[Index[iVar][jVar+1]], 
																														 Coord_Corner_Points[Index[iVar][jVar+2]], 
																														 Coord_Corner_Points[Index[iVar][jVar+3]]);
			
			double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]], 
																															Coord_Corner_Points[Index[iVar][jVar+1]], 
																															Coord_Corner_Points[Index[iVar][jVar+2]], 
																															Coord_Corner_Points[Index[iVar][jVar+3]]);
			if (Distance_Point*Distance_Vertex < 0.0) Inside = false;					
		}
		if (Inside) break;
	}
	
	return Inside;

}

void CFreeFormDefBox::SetDeformationZone(CGeometry *geometry, CConfig *config, unsigned short iFFDBox) {
	double *Coord;
	unsigned short iMarker, iVar, jVar;
	unsigned long iVertex, iPoint;
	bool Inside;
	
	unsigned short Index[5][7] = {
		{0, 1, 2, 5, 0, 1, 2},
		{0, 2, 7, 5, 0, 2, 7},
		{0, 2, 3, 7, 0, 2, 3},
		{0, 5, 7, 4, 0, 5, 7},
		{2, 7, 5, 6, 2, 7, 5}};
	
	for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
		if (config->GetMarker_All_DV(iMarker) == YES)
			for(iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {	
				iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
				geometry->node[iPoint]->SetMove(false);
				
				Coord = geometry->node[iPoint]->GetCoord();
				
				/*--- 1st tetrahedron {V0, V1, V2, V5}
				 2nd tetrahedron {V0, V2, V7, V5}
				 3th tetrahedron {V0, V2, V3, V7}
				 4th tetrahedron {V0, V5, V7, V4}
				 5th tetrahedron {V2, V7, V5, V6} ---*/
				
				for (iVar = 0; iVar < 5; iVar++) {
					Inside = true;
					for (jVar = 0; jVar < 4; jVar++) {
						double Distance_Point = geometry->Point2Plane_Distance(Coord, 
																																	 Coord_Corner_Points[Index[iVar][jVar+1]], 
																																	 Coord_Corner_Points[Index[iVar][jVar+2]], 
																																	 Coord_Corner_Points[Index[iVar][jVar+3]]);				
						double Distance_Vertex = geometry->Point2Plane_Distance(Coord_Corner_Points[Index[iVar][jVar]], 
																																		Coord_Corner_Points[Index[iVar][jVar+1]], 
																																		Coord_Corner_Points[Index[iVar][jVar+2]], 
																																		Coord_Corner_Points[Index[iVar][jVar+3]]);
						if (Distance_Point*Distance_Vertex < 0.0) Inside = false;					
					}
					if (Inside) break;
				}
				
				if (Inside) {
					geometry->node[iPoint]->SetMove(true);
				}
				
			}
}

double CFreeFormDefBox::GetDerivative1 (double *uvw, unsigned short val_diff, unsigned short *ijk, unsigned short *lmn) {
	unsigned short iDim;
	double value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 1);
	
	for (iDim = 0; iDim < nDim; iDim++)
		if (iDim != val_diff)
			value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	
	return value;	
}

double CFreeFormDefBox::GetDerivative2 (double *uvw, unsigned short dim, double *xyz, unsigned short *lmn) {
	
	unsigned short iDegree, jDegree, kDegree;
	double value = 0.0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++)
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] 
				* GetBernstein(lmn[0], iDegree, uvw[0])
				* GetBernstein(lmn[1], jDegree, uvw[1])
				* GetBernstein(lmn[2], kDegree, uvw[2]);
	
	return 2.0*(value - xyz[dim]);	
}

double CFreeFormDefBox::GetDerivative3(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short *lmn) {
	unsigned short iDegree, jDegree, kDegree, ijk[3];
	double value = 0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
				ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] * 
				GetDerivative1(uvw, diff_this, ijk, lmn);
			}
	
	return value;
}

double CFreeFormDefBox::GetDerivative4 (double *uvw, unsigned short val_diff, unsigned short val_diff2,
																			 unsigned short *ijk, unsigned short *lmn) {
	unsigned short iDim;
	double value;
	
	if (val_diff == val_diff2) {
		value = GetBernsteinDerivative(lmn[val_diff], ijk[val_diff], uvw[val_diff], 2);
		for (iDim = 0; iDim < nDim; iDim++)
			if (iDim != val_diff)
				value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	}
	else {
		value = GetBernsteinDerivative(lmn[val_diff],  ijk[val_diff],  uvw[val_diff], 1) *
		GetBernsteinDerivative(lmn[val_diff2], ijk[val_diff2], uvw[val_diff2], 1);
		for (iDim = 0; iDim < nDim; iDim++)
			if ((iDim != val_diff) && (iDim != val_diff2))
				value *= GetBernstein(lmn[iDim], ijk[iDim], uvw[iDim]);
	}
	
	return value;
}

double CFreeFormDefBox::GetDerivative5(double *uvw, unsigned short dim, unsigned short diff_this, unsigned short diff_this_also, 
																			unsigned short *lmn) {
	
	unsigned short iDegree, jDegree, kDegree, ijk[3];
	double value = 0.0;
	
	for (iDegree = 0; iDegree <= lmn[0]; iDegree++)
		for (jDegree = 0; jDegree <= lmn[1]; jDegree++)
			for (kDegree = 0; kDegree <= lmn[2]; kDegree++) {
				ijk[0] = iDegree; ijk[1] = jDegree; ijk[2] = kDegree;
				value += Coord_Control_Points[iDegree][jDegree][kDegree][dim] *
				GetDerivative4(uvw, diff_this, diff_this_also, ijk, lmn);
			}
	
	return value;
}

void CFreeFormDefBox::GetHessian_Analytical(double *uvw, double *xyz, double **val_Hessian) {
	
	unsigned short iDim, jDim;
	unsigned short l, m, n, lmn[3];
	
	/*--- Set the Degree of the Berstein polynomials ---*/
	lmn[0] = lDegree; lmn[1] = mDegree; lmn[2] = nDegree;
	
	/*--- Berstein polynomials degrees ---*/
	l = lmn[0]; m = lmn[1]; n = lmn[2];
	
	for (iDim = 0; iDim < nDim; iDim++)
		for (jDim = 0; jDim < nDim; jDim++)
			val_Hessian[iDim][jDim] = 0.0;
	
	/*--- Note that being all the functions linear combinations of polynomials, they are C^\infty,
	 and the Hessian will be symmetric; no need to compute the under-diagonal part, for example ---*/
	for (iDim = 0; iDim < nDim; iDim++) {
		val_Hessian[0][0] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,0,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,0,lmn);
		
		val_Hessian[1][1] += 2.0 * GetDerivative3(uvw,iDim,1,lmn) * GetDerivative3(uvw,iDim,1,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,1,1,lmn);
		
		val_Hessian[2][2] += 2.0 * GetDerivative3(uvw,iDim,2,lmn) * GetDerivative3(uvw,iDim,2,lmn) + 
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,2,2,lmn);
		
		val_Hessian[0][1] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,1,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,1,lmn);
		
		val_Hessian[0][2] += 2.0 * GetDerivative3(uvw,iDim,0,lmn) * GetDerivative3(uvw,iDim,2,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,0,2,lmn);
		
		val_Hessian[1][2] += 2.0 * GetDerivative3(uvw,iDim,1,lmn) * GetDerivative3(uvw,iDim,2,lmn) +
		GetDerivative2(uvw,iDim,xyz,lmn) * GetDerivative5(uvw,iDim,1,2,lmn);
	}
	
	val_Hessian[1][0] = val_Hessian[0][1];
	val_Hessian[2][0] = val_Hessian[0][2];
	val_Hessian[2][1] = val_Hessian[1][2];
}

void CFreeFormDefBox::Gauss_Elimination(double** A, double* rhs, unsigned short nVar) {
	unsigned short jVar, kVar, iVar;
    double weight, aux;
	
	if (nVar == 1) {
    if (fabs(A[0][0]) < EPS) cout <<"Gauss' elimination error, value:" << abs(A[0][0]) << "." << endl;
		rhs[0] /= A[0][0];
  }
	else {
		/*--- Transform system in Upper Matrix ---*/
		for (iVar = 1; iVar < nVar; iVar++) {
			for (jVar = 0; jVar < iVar; jVar++) {
        if (fabs(A[jVar][jVar]) < EPS) cout <<"Gauss' elimination error, value:" << abs(A[jVar][jVar]) << "." << endl;
				weight = A[iVar][jVar]/A[jVar][jVar];
				for (kVar = jVar; kVar < nVar; kVar++)
					A[iVar][kVar] -= weight*A[jVar][kVar];
				rhs[iVar] -= weight*rhs[jVar];
			}
		}
		/*--- Backwards substitution ---*/
    if (fabs(A[nVar-1][nVar-1]) < EPS) cout <<"Gauss' elimination error, value:" << abs(A[nVar-1][nVar-1]) << "." << endl;
		rhs[nVar-1] = rhs[nVar-1]/A[nVar-1][nVar-1];
		for (short iVar = nVar-2; iVar >= 0; iVar--) {
			aux = 0;
			for (jVar = iVar+1; jVar < nVar; jVar++)
				aux += A[iVar][jVar]*rhs[jVar];
      if (fabs(A[iVar][iVar]) < EPS) cout <<"Gauss' elimination error, value:" << abs(A[iVar][iVar]) << "." << endl;
			rhs[iVar] = (rhs[iVar]-aux)/A[iVar][iVar];
			if (iVar == 0) break;
		}
	}
}
