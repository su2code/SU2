/*!
 * \file geometry_structure.cpp
 * \brief Main subroutines for creating the primal grid and multigrid structure.
 * \author F. Palacios
 * \version 3.2.5 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

#include "../include/geometry_structure.hpp"

CGeometry::CGeometry(void) {
  
  nEdge = 0;
  nPoint = 0;
  nElem = 0;
  
  nElem_Bound = NULL;
  Tag_to_Marker = NULL;
  elem = NULL;
  face = NULL;
  bound = NULL;
  node = NULL;
  edge = NULL;
  vertex = NULL;
  nVertex = NULL;
  newBound = NULL;
  nNewElem_Bound = NULL;
  Marker_All_SendRecv = NULL;
  
  //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].clear();
  //	PeriodicElem[MAX_NUMBER_PERIODIC].clear();
  //	OldBoundaryElems[MAX_NUMBER_MARKER].clear();
  //	XCoordList.clear();
  
  //	Xcoord_plane.clear();
  //	Ycoord_plane.clear();
  //	Zcoord_plane.clear();
  //	FaceArea_plane.clear();
  //	Plane_points.clear();
  
}

CGeometry::~CGeometry(void) {
  unsigned long iElem, iElem_Bound, iPoint, iFace, iVertex, iEdge;
  unsigned short iMarker;
  
  if (elem != NULL) {
    for (iElem = 0; iElem < nElem; iElem++)
      if (elem[iElem] != NULL) delete elem[iElem];
    delete[] elem;
  }
  
  if (bound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
      }
    }
    delete[] bound;
  }
  
  if (face != NULL) {
    for (iFace = 0; iFace < nFace; iFace ++)
      if (face[iFace] != NULL) delete face[iFace];
    delete[] face;
  }
  
  if (node != NULL) {
    for (iPoint = 0; iPoint < nPoint; iPoint ++)
      if (node[iPoint] != NULL) delete node[iPoint];
    delete[] node;
  }
  
  if (edge != NULL) {
    for (iEdge = 0; iEdge < nEdge; iEdge ++)
      if (edge[iEdge] != NULL) delete edge[iEdge];
    delete[] edge;
  }
  
  if (vertex != NULL)  {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        if (vertex[iMarker][iVertex] != NULL) delete vertex[iMarker][iVertex];
      }
    }
    delete[] vertex;
  }
  
  if (newBound != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (newBound[iMarker][iElem_Bound] != NULL) delete newBound[iMarker][iElem_Bound];
      }
    }
    delete[] newBound;
  }
  
  if (nElem_Bound != NULL) delete[] nElem_Bound;
  if (nVertex != NULL) delete[] nVertex;
  if (nNewElem_Bound != NULL) delete[] nNewElem_Bound;
  if (Marker_All_SendRecv != NULL) delete[] Marker_All_SendRecv;
  if (Tag_to_Marker != NULL) delete[] Tag_to_Marker;
  
  //	PeriodicPoint[MAX_NUMBER_PERIODIC][2].~vector();
  //	PeriodicElem[MAX_NUMBER_PERIODIC].~vector();
  //	OldBoundaryElems[MAX_NUMBER_MARKER].~vector();
  //	XCoordList.~vector();
  
  //	Xcoord_plane.~vector()
  //	Ycoord_plane.~vector()
  //	Zcoord_plane.~vector()
  //	FaceArea_plane.~vector()
  //	Plane_points.~vector()
  
}

double CGeometry::Point2Plane_Distance(double *Coord, double *iCoord, double *jCoord, double *kCoord) {
  double CrossProduct[3], iVector[3], jVector[3], distance, modulus;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim ++) {
    iVector[iDim] = jCoord[iDim] - iCoord[iDim];
    jVector[iDim] = kCoord[iDim] - iCoord[iDim];
  }
  
  CrossProduct[0] = iVector[1]*jVector[2] - iVector[2]*jVector[1];
  CrossProduct[1] = iVector[2]*jVector[0] - iVector[0]*jVector[2];
  CrossProduct[2] = iVector[0]*jVector[1] - iVector[1]*jVector[0];
  
  modulus = sqrt(CrossProduct[0]*CrossProduct[0]+CrossProduct[1]*CrossProduct[1]+CrossProduct[2]*CrossProduct[2]);
  
  distance = 0.0;
  for (iDim = 0; iDim < 3; iDim ++)
    distance += CrossProduct[iDim]*(Coord[iDim]-iCoord[iDim]);
  distance /= modulus;
  
  return distance;
  
}

long CGeometry::FindEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return node[first_point]->GetEdge(iNode);
  else {
    cout << "\n\n   !!! Error !!!\n" << endl;
    cout <<"Can't find the edge that connects "<< first_point <<" and "<< second_point <<"."<< endl;
    exit(EXIT_FAILURE);
    return -1;
  }
}

bool CGeometry::CheckEdge(unsigned long first_point, unsigned long second_point) {
  unsigned long iPoint = 0;
  unsigned short iNode;
  for (iNode = 0; iNode < node[first_point]->GetnPoint(); iNode++) {
    iPoint = node[first_point]->GetPoint(iNode);
    if (iPoint == second_point) break;
  }
  
  if (iPoint == second_point) return true;
  else return false;
  
}

void CGeometry::SetEdges(void) {
  unsigned long iPoint, jPoint;
  long iEdge;
  unsigned short jNode, iNode;
  long TestEdge = 0;
  
  nEdge = 0;
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
        if (node[jPoint]->GetPoint(jNode) == iPoint) {
          TestEdge = node[jPoint]->GetEdge(jNode);
          break;
        }
      if (TestEdge == -1) {
        node[iPoint]->SetEdge(nEdge, iNode);
        node[jPoint]->SetEdge(nEdge, jNode);
        nEdge++;
      }
    }
  
  edge = new CEdge*[nEdge];
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
      jPoint = node[iPoint]->GetPoint(iNode);
      iEdge = FindEdge(iPoint, jPoint);
      if (iPoint < jPoint) edge[iEdge] = new CEdge(iPoint, jPoint, nDim);
    }
}

void CGeometry::SetFaces(void) {
  //	unsigned long iPoint, jPoint, iFace;
  //	unsigned short jNode, iNode;
  //	long TestFace = 0;
  //
  //	nFace = 0;
  //	for(iPoint = 0; iPoint < nPoint; iPoint++)
  //		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			for(jNode = 0; jNode < node[jPoint]->GetnPoint(); jNode++)
  //				if (node[jPoint]->GetPoint(jNode) == iPoint) {
  //					TestFace = node[jPoint]->GetFace(jNode);
  //					break;
  //				}
  //			if (TestFace == -1) {
  //				node[iPoint]->SetFace(nFace, iNode);
  //				node[jPoint]->SetFace(nFace, jNode);
  //				nFace++;
  //			}
  //		}
  //
  //	face = new CFace*[nFace];
  //
  //	for(iPoint = 0; iPoint < nPoint; iPoint++)
  //		for(iNode = 0; iNode < node[iPoint]->GetnPoint(); iNode++) {
  //			jPoint = node[iPoint]->GetPoint(iNode);
  //			iFace = FindFace(iPoint, jPoint);
  //			if (iPoint < jPoint) face[iFace] = new CFace(iPoint,jPoint,nDim);
  //		}
}

void CGeometry::TestGeometry(void) {
  
  ofstream para_file;
  
  para_file.open("test_geometry.dat", ios::out);
  
  double *Normal = new double[nDim];
  
  for(unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {
    para_file << "Edge index: " << iEdge << endl;
    para_file << "   Point index: " << edge[iEdge]->GetNode(0) << "\t" << edge[iEdge]->GetNode(1) << endl;
    edge[iEdge]->GetNormal(Normal);
    para_file << "      Face normal : ";
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      para_file << Normal[iDim] << "\t";
    para_file << endl;
  }
  
  para_file << endl;
  para_file << endl;
  para_file << endl;
  para_file << endl;
  
  for(unsigned short iMarker =0; iMarker < nMarker; iMarker++) {
    para_file << "Marker index: " << iMarker << endl;
    for(unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      para_file << "   Vertex index: " << iVertex << endl;
      para_file << "      Point index: " << vertex[iMarker][iVertex]->GetNode() << endl;
      para_file << "      Point coordinates : ";
      for(unsigned short iDim = 0; iDim < nDim; iDim++) {
        para_file << node[vertex[iMarker][iVertex]->GetNode()]->GetCoord(iDim) << "\t";}
      para_file << endl;
      vertex[iMarker][iVertex]->GetNormal(Normal);
      para_file << "         Face normal : ";
      for(unsigned short iDim = 0; iDim < nDim; iDim++)
        para_file << Normal[iDim] << "\t";
      para_file << endl;
    }
  }
  
}

void CGeometry::SetSpline(vector<double> &x, vector<double> &y, unsigned long n, double yp1, double ypn, vector<double> &y2) {
  unsigned long i, k;
  double p, qn, sig, un, *u;
  
  u = new double [n];
  
  if (yp1 > 0.99e30)			// The lower boundary condition is set either to be "nat
    y2[0]=u[0]=0.0;			  // -ural"
  else {									// or else to have a specified first derivative.
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  
  for (i=2; i<=n-1; i++) {									//  This is the decomposition loop of the tridiagonal al-
    sig=(x[i-1]-x[i-2])/(x[i]-x[i-2]);		//	gorithm. y2 and u are used for tem-
    p=sig*y2[i-2]+2.0;										//	porary storage of the decomposed
    y2[i-1]=(sig-1.0)/p;										//	factors.
    
    double a1 = (y[i]-y[i-1])/(x[i]-x[i-1]); if (x[i] == x[i-1]) a1 = 1.0;
    double a2 = (y[i-1]-y[i-2])/(x[i-1]-x[i-2]); if (x[i-1] == x[i-2]) a2 = 1.0;
    u[i-1]= a1 - a2;
    u[i-1]=(6.0*u[i-1]/(x[i]-x[i-2])-sig*u[i-2])/p;
    
  }
  
  if (ypn > 0.99e30)						// The upper boundary condition is set either to be
    qn=un=0.0;									// "natural"
  else {												// or else to have a specified first derivative.
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-1; k>=1; k--)					// This is the backsubstitution loop of the tridiagonal
    y2[k-1]=y2[k-1]*y2[k]+u[k-1];	  // algorithm.
  
  delete[] u;
  
}

double CGeometry::GetSpline(vector<double>&xa, vector<double>&ya, vector<double>&y2a, unsigned long n, double x) {
  unsigned long klo, khi, k;
  double h, b, a, y;
  
  if (x < xa[0]) x = xa[0];       // Clip max and min values
  if (x > xa[n-1]) x = xa[n-1];
  
  klo = 1;										// We will find the right place in the table by means of
  khi = n;										// bisection. This is optimal if sequential calls to this
  while (khi-klo > 1) {			// routine are at random values of x. If sequential calls
    k = (khi+klo) >> 1;				// are in order, and closely spaced, one would do better
    if (xa[k-1] > x) khi = k;		// to store previous values of klo and khi and test if
    else klo=k;							// they remain appropriate on the next call.
  }								// klo and khi now bracket the input value of x
  h = xa[khi-1] - xa[klo-1];
  if (h == 0.0) h = EPS; // cout << "Bad xa input to routine splint" << endl;	// The xa’s must be distinct.
  a = (xa[khi-1]-x)/h;
  b = (x-xa[klo-1])/h;				// Cubic spline polynomial is now evaluated.
  y = a*ya[klo-1]+b*ya[khi-1]+((a*a*a-a)*y2a[klo-1]+(b*b*b-b)*y2a[khi-1])*(h*h)/6.0;
  
  return y;
}

unsigned short CGeometry::ComputeSegmentPlane_Intersection(double *Segment_P0, double *Segment_P1, double Variable_P0, double Variable_P1,
                                                           double *Plane_P0, double *Plane_Normal, double *Intersection, double &Variable_Interp) {
  double u[3], v[3], Denominator, Numerator, Aux, ModU;
  unsigned short iDim;
  
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P1[iDim] - Segment_P0[iDim];
    v[iDim] = Plane_P0[iDim] - Segment_P0[iDim];
  }
  
  ModU = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
  
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  
  if (fabs(Denominator) <= 0.0) return 0; // No intersection.
  
  Aux = Numerator / Denominator;
  
  if (Aux < 0.0 || Aux > 1.0) return 0; // No intersection.
  
  for (iDim = 0; iDim < 3; iDim++)
    Intersection[iDim] = Segment_P0[iDim] + Aux * u[iDim];
  
  
  /*--- Check that the intersection is in the segment ---*/
  for (iDim = 0; iDim < 3; iDim++) {
    u[iDim] = Segment_P0[iDim] - Intersection[iDim];
    v[iDim] = Segment_P1[iDim] - Intersection[iDim];
  }
  
  Variable_Interp = Variable_P0 + (Variable_P1 - Variable_P0)*sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])/ModU;
  
  Denominator = Plane_Normal[0]*u[0] + Plane_Normal[1]*u[1] + Plane_Normal[2]*u[2];
  Numerator = Plane_Normal[0]*v[0] + Plane_Normal[1]*v[1] + Plane_Normal[2]*v[2];
  
  Aux = Numerator * Denominator;
  
  if (Aux > 0.0) return 3; // Intersection outside the segment.
  else return 1;
  
}

void CGeometry::ComputeAirfoil_Section(double *Plane_P0, double *Plane_Normal, unsigned short iSection,
                                       double MinXCoord, double MaxXCoord, double *FlowVariable,
                                       vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil,
                                       vector<double> &Zcoord_Airfoil, vector<double> &Variable_Airfoil,
                                       bool original_surface, CConfig *config) {
  
  unsigned short iMarker, iNode, jNode, iDim, intersect;
  long MinDist_Point, MinDistAngle_Point;
  unsigned long iPoint, jPoint, iElem, Trailing_Point, Airfoil_Point, iVertex, jVertex;
  double Segment_P0[3] = {0.0, 0.0, 0.0}, Segment_P1[3] = {0.0, 0.0, 0.0}, Variable_P0 = 0.0, Variable_P1 = 0.0, Intersection[3] = {0.0, 0.0, 0.0}, Trailing_Coord, MinDist_Value, MinDistAngle_Value, Dist_Value,
  Airfoil_Tangent[3] = {0.0, 0.0, 0.0}, Segment[3] = {0.0, 0.0, 0.0}, Length, Angle_Value, MaxAngle = 30, *VarCoord = NULL, CosValue, Variable_Interp;
  vector<double> Xcoord, Ycoord, Zcoord, Variable;
  vector<unsigned long> Duplicate;
  vector<unsigned long>::iterator it;
  int rank = MASTER_NODE;
  double **Coord_Variation = NULL;
  
#ifdef HAVE_MPI
  unsigned long nLocalVertex, nGlobalVertex, MaxLocalVertex, *Buffer_Send_nVertex, *Buffer_Receive_nVertex, nBuffer;
  int nProcessor, iProcessor;
  double *Buffer_Send_Coord, *Buffer_Receive_Coord;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  Xcoord_Airfoil.clear();
  Ycoord_Airfoil.clear();
  Zcoord_Airfoil.clear();
  Variable_Airfoil.clear();
  
  /*--- Set the right plane in 2D (note the change in Y-Z plane) ---*/
  
  if (nDim == 2) {
    iSection = 0;
    Plane_P0[0] = 0.0;      Plane_P0[1] = 0.0;      Plane_P0[2] = 0.0;
    Plane_Normal[0] = 0.0;  Plane_Normal[1] = 1.0;  Plane_Normal[2] = 0.0;
  }
  
  /*--- the grid variation is stored using a vertices information,
   we should go from vertex to points ---*/
  
  if (original_surface == false) {
    
    Coord_Variation = new double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      Coord_Variation[iPoint] = new double [nDim];
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_GeoEval(iMarker) == YES) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          VarCoord = vertex[iMarker][iVertex]->GetVarCoord();
          iPoint = vertex[iMarker][iVertex]->GetNode();
          for (iDim = 0; iDim < nDim; iDim++)
            Coord_Variation[iPoint][iDim] = VarCoord[iDim];
        }
      }
    }
    
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_GeoEval(iMarker) == YES) {
      for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
          iPoint = bound[iMarker][iElem]->GetNode(iNode);
          for(jNode = 0; jNode < bound[iMarker][iElem]->GetnNodes(); jNode++) {
            jPoint = bound[iMarker][iElem]->GetNode(jNode);
            
            if ((jPoint > iPoint) && ((node[iPoint]->GetCoord(0) > MinXCoord) && (node[iPoint]->GetCoord(0) < MaxXCoord))) {
              
              Segment_P0[0] = 0.0;  Segment_P0[1] = 0.0;  Segment_P0[2] = 0.0;  Variable_P0 = 0.0;
              Segment_P1[0] = 0.0;  Segment_P1[1] = 0.0;  Segment_P1[2] = 0.0;  Variable_P1 = 0.0;
              
              for (iDim = 0; iDim < nDim; iDim++) {
                if (original_surface == true) {
                  Segment_P0[iDim] = node[iPoint]->GetCoord(iDim);
                  Segment_P1[iDim] = node[jPoint]->GetCoord(iDim);
                }
                else {
                  Segment_P0[iDim] = node[iPoint]->GetCoord(iDim) + Coord_Variation[iPoint][iDim];
                  Segment_P1[iDim] = node[jPoint]->GetCoord(iDim) + Coord_Variation[jPoint][iDim];
                }
              }
              
              if (FlowVariable != NULL) {
                Variable_P0 = FlowVariable[iPoint];
                Variable_P1 = FlowVariable[jPoint];
              }
              
              /*--- In 2D add the points directly (note the change between Y and Z coordinate) ---*/
              
              if (nDim == 2) {
                Xcoord.push_back(Segment_P0[0]);    Xcoord.push_back(Segment_P1[0]);
                Ycoord.push_back(Segment_P0[2]);    Ycoord.push_back(Segment_P1[2]);
                Zcoord.push_back(Segment_P0[1]);    Zcoord.push_back(Segment_P1[1]);
                Variable.push_back(Variable_P0);    Variable.push_back(Variable_P1);
              }
              /*--- In 3D compute the intersection ---*/
              
              else if (nDim == 3) {
                intersect = ComputeSegmentPlane_Intersection(Segment_P0, Segment_P1, Variable_P0, Variable_P1, Plane_P0, Plane_Normal, Intersection, Variable_Interp);
                if (intersect == 1) {
                  Xcoord.push_back(Intersection[0]);
                  Ycoord.push_back(Intersection[1]);
                  Zcoord.push_back(Intersection[2]);
                  Variable.push_back(Variable_Interp);
                }
              }
              
            }
            
          }
        }
      }
    }
  }
  
  if (original_surface == false) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete [] Coord_Variation[iPoint];
    delete [] Coord_Variation;
    
  }
  
  
#ifdef HAVE_MPI
  
  /*--- Copy the coordinates of all the points in the plane to the master node ---*/
  
  nLocalVertex = 0, nGlobalVertex = 0, MaxLocalVertex = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  Buffer_Send_nVertex = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  nLocalVertex = Xcoord.size();
  
  Buffer_Send_nVertex[0] = nLocalVertex;
  
  MPI_Allreduce(&nLocalVertex, &nGlobalVertex, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  Buffer_Send_Coord = new double [MaxLocalVertex*4];
  Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex*4];
  nBuffer = MaxLocalVertex*4;
  
  for (iVertex = 0; iVertex < nLocalVertex; iVertex++) {
    Buffer_Send_Coord[iVertex*4 + 0] = Xcoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 1] = Ycoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 2] = Zcoord[iVertex];
    Buffer_Send_Coord[iVertex*4 + 3] = Variable[iVertex];
  }
  
  MPI_Allgather(Buffer_Send_Coord,nBuffer,MPI_DOUBLE,Buffer_Receive_Coord,nBuffer,MPI_DOUBLE,MPI_COMM_WORLD);
  
  /*--- Clean the vectors before adding the new vertices only to the master node ---*/
  
  Xcoord.clear();
  Ycoord.clear();
  Zcoord.clear();
  Variable.clear();
  
  /*--- Copy the boundary to the master node vectors ---*/
  if (rank == MASTER_NODE) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Xcoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 0] );
        Ycoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 1] );
        Zcoord.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 2] );
        Variable.push_back( Buffer_Receive_Coord[ iProcessor*MaxLocalVertex*4 + iVertex*4 + 3] );
      }
    }
  }
  
  delete[] Buffer_Send_Coord;   delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex; delete[] Buffer_Receive_nVertex;
  
#endif
  
  if ((rank == MASTER_NODE) && (Xcoord.size() != 0)) {
    
    /*--- Create a list with the duplicated points ---*/
    
    for (iVertex = 0; iVertex < Xcoord.size()-1; iVertex++) {
      for (jVertex = iVertex+1; jVertex < Xcoord.size(); jVertex++) {
        Segment[0] = Xcoord[jVertex] - Xcoord[iVertex];
        Segment[1] = Ycoord[jVertex] - Ycoord[iVertex];
        Segment[2] = Zcoord[jVertex] - Zcoord[iVertex];
        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        if (Dist_Value < 1E-6) {
          Duplicate.push_back (jVertex);
        }
      }
    }
    
    sort(Duplicate.begin(), Duplicate.end());
    it = unique(Duplicate.begin(), Duplicate.end());
    Duplicate.resize(it - Duplicate.begin());
    
    /*--- Remove duplicated points (starting from the back) ---*/
    
    for (iVertex = Duplicate.size(); iVertex > 0; iVertex--) {
      Xcoord.erase (Xcoord.begin() + Duplicate[iVertex-1]);
      Ycoord.erase (Ycoord.begin() + Duplicate[iVertex-1]);
      Zcoord.erase (Zcoord.begin() + Duplicate[iVertex-1]);
      Variable.erase (Variable.begin() + Duplicate[iVertex-1]);
    }
    
    if (Xcoord.size() != 1) {
      
      /*--- Find the trailing edge ---*/
      
      Trailing_Point = 0; Trailing_Coord = Xcoord[0];
      for (iVertex = 1; iVertex < Xcoord.size(); iVertex++) {
        if (Xcoord[iVertex] > Trailing_Coord) {
          Trailing_Point = iVertex; Trailing_Coord = Xcoord[iVertex];
        }
      }
      
      /*--- Add the trailing edge to the list, and remove from the original list ---*/
      Xcoord_Airfoil.push_back(Xcoord[Trailing_Point]); Ycoord_Airfoil.push_back(Ycoord[Trailing_Point]); Zcoord_Airfoil.push_back(Zcoord[Trailing_Point]); Variable_Airfoil.push_back(Variable[Trailing_Point]);
      Xcoord.erase (Xcoord.begin() + Trailing_Point); Ycoord.erase (Ycoord.begin() + Trailing_Point); Zcoord.erase (Zcoord.begin() + Trailing_Point); Variable.erase (Variable.begin() + Trailing_Point);
      
      /*--- Find the next point using the right hand side rule ---*/
      MinDist_Value = 1E6; MinDist_Point = 0;
      for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
        Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[0];
        Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[0];
        Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[0];
        Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
        Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
        
        if ((Dist_Value < MinDist_Value) && (Segment[2] > 0.0)) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
      }
      
      Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);  Variable_Airfoil.push_back(Variable[MinDist_Point]);
      Xcoord.erase (Xcoord.begin() + MinDist_Point);    Ycoord.erase (Ycoord.begin() + MinDist_Point);    Zcoord.erase (Zcoord.begin() + MinDist_Point);    Variable.erase (Variable.begin() + MinDist_Point);
      
      /*--- Algorithm for the rest of the points ---*/
      do {
        
        /*--- Last added point in the list ---*/
        Airfoil_Point = Xcoord_Airfoil.size() - 1;
        
        /*--- Compute the slope of the curve ---*/
        Airfoil_Tangent[0] = Xcoord_Airfoil[Airfoil_Point] - Xcoord_Airfoil[Airfoil_Point-1];
        Airfoil_Tangent[1] = Ycoord_Airfoil[Airfoil_Point] - Ycoord_Airfoil[Airfoil_Point-1];
        Airfoil_Tangent[2] = Zcoord_Airfoil[Airfoil_Point] - Zcoord_Airfoil[Airfoil_Point-1];
        Length = sqrt(pow(Airfoil_Tangent[0], 2.0) + pow(Airfoil_Tangent[1], 2.0) + pow(Airfoil_Tangent[2], 2.0));
        Airfoil_Tangent[0] /= Length; Airfoil_Tangent[1] /= Length; Airfoil_Tangent[2] /= Length;
        
        /*--- Find the closest point with the right slope ---*/
        MinDist_Value = 1E6; MinDistAngle_Value = 180;
        MinDist_Point = -1; MinDistAngle_Point = -1;
        for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
          
          Segment[0] = Xcoord[iVertex] - Xcoord_Airfoil[Airfoil_Point];
          Segment[1] = Ycoord[iVertex] - Ycoord_Airfoil[Airfoil_Point];
          Segment[2] = Zcoord[iVertex] - Zcoord_Airfoil[Airfoil_Point];
          
          /*--- Compute the distance to each point ---*/
          Dist_Value = sqrt(pow(Segment[0], 2.0) + pow(Segment[1], 2.0) + pow(Segment[2], 2.0));
          
          /*--- Compute the angle of the point ---*/
          Segment[0] /= Dist_Value; Segment[1] /= Dist_Value; Segment[2] /= Dist_Value;
          
          /*--- Clip the value of the cosine, this is important due to the round errors ---*/
          CosValue = Airfoil_Tangent[0]*Segment[0] + Airfoil_Tangent[1]*Segment[1] + Airfoil_Tangent[2]*Segment[2];
          if (CosValue >= 1.0) CosValue = 1.0;
          if (CosValue <= -1.0) CosValue = -1.0;
          
          Angle_Value = acos(CosValue) * 180 / PI_NUMBER;
          
          if (Dist_Value < MinDist_Value) { MinDist_Point = iVertex; MinDist_Value = Dist_Value; }
          if ((Dist_Value < MinDistAngle_Value) && (Angle_Value < MaxAngle)) {MinDistAngle_Point = iVertex; MinDistAngle_Value = Dist_Value;}
          
        }
        
        if ( MinDistAngle_Point != -1) MinDist_Point = MinDistAngle_Point;
        
        /*--- Add and remove the min distance to the list ---*/
        Xcoord_Airfoil.push_back(Xcoord[MinDist_Point]);  Ycoord_Airfoil.push_back(Ycoord[MinDist_Point]);  Zcoord_Airfoil.push_back(Zcoord[MinDist_Point]);    Variable_Airfoil.push_back(Variable[MinDist_Point]);
        Xcoord.erase(Xcoord.begin() + MinDist_Point);     Ycoord.erase(Ycoord.begin() + MinDist_Point);     Zcoord.erase(Zcoord.begin() + MinDist_Point);       Variable.erase(Variable.begin() + MinDist_Point);
        
      } while (Xcoord.size() != 0);
      
      /*--- Clean the vector before using them again for storing the upper or the lower side ---*/
      
      Xcoord.clear(); Ycoord.clear(); Zcoord.clear(); Variable.clear();
      
    }
    
    
  }
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}

void CGeometry::ComputeSurf_Curvature(CConfig *config) {
  unsigned short iMarker, iNeigh_Point, iDim, iNode, iNeighbor_Nodes, Neighbor_Node;
  unsigned long Neighbor_Point, iVertex, iPoint, jPoint,iElem_Bound, iEdge, nLocalVertex, nGlobalVertex , MaxLocalVertex , *Buffer_Send_nVertex, *Buffer_Receive_nVertex, nBuffer, TotalnPointDomain;
  int iProcessor, nProcessor;
  vector<unsigned long> Point_NeighborList, Elem_NeighborList, Point_Triangle, Point_Edge, Point_Critical;
  vector<unsigned long>::iterator it;
  double U[3], V[3], W[3], Length_U, Length_V, Length_W, CosValue, Angle_Value, *K, *Angle_Defect, *Area_Vertex, *Angle_Alpha, *Angle_Beta, **NormalMeanK, MeanK, GaussK, MaxPrinK, MinPrinK, cot_alpha, cot_beta, delta, X1, X2, X3, Y1, Y2, Y3, radius, *Buffer_Send_Coord, *Buffer_Receive_Coord, *Coord, Dist, MinDist, MaxK, MinK, SigmaK;
  bool *Check_Edge;
  int rank;
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Allocate surface curvature ---*/
  K = new double [nPoint];
  
  if (nDim == 2) {
    
    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        
        /*--- Loop through all marker vertices again, this time also
         finding the neighbors of each node.---*/
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();
          
          if (node[iPoint]->GetDomain()) {
            /*--- Loop through neighbors. In 2-D, there should be 2 nodes on either
             side of this vertex that lie on the same surface. ---*/
            Point_Edge.clear();
            
            for (iNeigh_Point = 0; iNeigh_Point < node[iPoint]->GetnPoint(); iNeigh_Point++) {
              Neighbor_Point = node[iPoint]->GetPoint(iNeigh_Point);
              
              /*--- Check if this neighbor lies on the surface. If so,
               add to the list of neighbors. ---*/
              if (node[Neighbor_Point]->GetPhysicalBoundary()) {
                Point_Edge.push_back(Neighbor_Point);
              }
              
            }
            
            if (Point_Edge.size() == 2) {
              
              /*--- Compute the curvature using three points ---*/
              X1 = node[iPoint]->GetCoord(0);
              X2 = node[Point_Edge[0]]->GetCoord(0);
              X3 = node[Point_Edge[1]]->GetCoord(0);
              Y1 = node[iPoint]->GetCoord(1);
              Y2 = node[Point_Edge[0]]->GetCoord(1);
              Y3 = node[Point_Edge[1]]->GetCoord(1);
              
              radius = sqrt(((X2-X1)*(X2-X1) + (Y2-Y1)*(Y2-Y1))*
                            ((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))*
                            ((X3-X1)*(X3-X1) + (Y3-Y1)*(Y3-Y1)))/
              (2.0*fabs(X1*Y2+X2*Y3+X3*Y1-X1*Y3-X2*Y1-X3*Y2));
              
              K[iPoint] = 1.0/radius;
              node[iPoint]->SetCurvature(K[iPoint]);
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  else {
    
    Angle_Defect = new double [nPoint];
    Area_Vertex = new double [nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      Angle_Defect[iPoint] = 2*PI_NUMBER;
      Area_Vertex[iPoint] = 0.0;
    }
    
    Angle_Alpha = new double [nEdge];
    Angle_Beta = new double [nEdge];
    Check_Edge = new bool [nEdge];
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      Angle_Alpha[iEdge] = 0.0;
      Angle_Beta[iEdge] = 0.0;
      Check_Edge[iEdge] = true;
    }
    
    NormalMeanK = new double *[nPoint];
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      NormalMeanK[iPoint] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        NormalMeanK[iPoint][iDim] = 0.0;
      }
    }
    
    /*--- Loop over all the markers ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        
        /*--- Loop over all the boundary elements ---*/
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          
          /*--- Only triangles ---*/
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
            
            /*--- Loop over all the nodes of the boundary element ---*/
            for(iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
              
              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
              
              Point_Triangle.clear();
              
              for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                Neighbor_Point = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                Point_Triangle.push_back(Neighbor_Point);
              }
              
              iEdge = FindEdge(Point_Triangle[0], Point_Triangle[1]);
              
              for (iDim = 0; iDim < nDim; iDim++) {
                U[iDim] = node[Point_Triangle[0]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
                V[iDim] = node[Point_Triangle[1]]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim);
              }
              
              W[0] = 0.5*(U[1]*V[2]-U[2]*V[1]); W[1] = -0.5*(U[0]*V[2]-U[2]*V[0]); W[2] = 0.5*(U[0]*V[1]-U[1]*V[0]);
              
              Length_U = 0.0, Length_V = 0.0, Length_W = 0.0, CosValue = 0.0;
              for (iDim = 0; iDim < nDim; iDim++) { Length_U += U[iDim]*U[iDim]; Length_V += V[iDim]*V[iDim]; Length_W += W[iDim]*W[iDim]; }
              Length_U = sqrt(Length_U); Length_V = sqrt(Length_V); Length_W = sqrt(Length_W);
              for (iDim = 0; iDim < nDim; iDim++) { U[iDim] /= Length_U; V[iDim] /= Length_V; CosValue += U[iDim]*V[iDim]; }
              if (CosValue >= 1.0) CosValue = 1.0;
              if (CosValue <= -1.0) CosValue = -1.0;
              
              Angle_Value = acos(CosValue);
              Area_Vertex[iPoint] += Length_W;
              Angle_Defect[iPoint] -= Angle_Value;
              if (Angle_Alpha[iEdge] == 0.0) Angle_Alpha[iEdge] = Angle_Value;
              else Angle_Beta[iEdge] = Angle_Value;
              
            }
          }
        }
      }
    }
    
    /*--- Compute mean curvature ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE) {
            for(iNode = 0; iNode < bound[iMarker][iElem_Bound]->GetnNodes(); iNode++) {
              iPoint = bound[iMarker][iElem_Bound]->GetNode(iNode);
              
              for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem_Bound]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
                Neighbor_Node = bound[iMarker][iElem_Bound]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
                jPoint = bound[iMarker][iElem_Bound]->GetNode(Neighbor_Node);
                
                iEdge = FindEdge(iPoint, jPoint);
                
                if (Check_Edge[iEdge]) {
                  
                  Check_Edge[iEdge] = false;
                  
                  if (tan(Angle_Alpha[iEdge]) != 0.0) cot_alpha = 1.0/tan(Angle_Alpha[iEdge]); else cot_alpha = 0.0;
                  if (tan(Angle_Beta[iEdge]) != 0.0) cot_beta = 1.0/tan(Angle_Beta[iEdge]); else cot_beta = 0.0;
                  
                  /*--- iPoint, and jPoint ---*/
                  for (iDim = 0; iDim < nDim; iDim++) {
                    if (Area_Vertex[iPoint] != 0.0) NormalMeanK[iPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[iPoint]->GetCoord(iDim) - node[jPoint]->GetCoord(iDim)) / Area_Vertex[iPoint];
                    if (Area_Vertex[jPoint] != 0.0) NormalMeanK[jPoint][iDim] += 3.0 * (cot_alpha + cot_beta) * (node[jPoint]->GetCoord(iDim) - node[iPoint]->GetCoord(iDim)) / Area_Vertex[jPoint];
                  }
                }
                
              }
            }
          }
        }
      }
    }
    
    /*--- Compute Gauss, mean, max and min principal curvature,
     and set the list of critical points ---*/
    
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint  = vertex[iMarker][iVertex]->GetNode();
          
          if (node[iPoint]->GetDomain()) {
            
            if (Area_Vertex[iPoint] != 0.0) GaussK = 3.0*Angle_Defect[iPoint]/Area_Vertex[iPoint];
            else GaussK = 0.0;
            
            MeanK = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              MeanK += NormalMeanK[iPoint][iDim]*NormalMeanK[iPoint][iDim];
            MeanK = sqrt(MeanK);
            
            delta = max((MeanK*MeanK - GaussK), 0.0);
            
            MaxPrinK = MeanK + sqrt(delta);
            MinPrinK = MeanK - sqrt(delta);
            
            /*--- Store the curvature value ---*/
            K[iPoint] = MaxPrinK;
            node[iPoint]->SetCurvature(K[iPoint]);
          }
          
        }
      }
    }
    
    delete [] Angle_Defect;
    delete [] Area_Vertex;
    delete [] Angle_Alpha;
    delete [] Angle_Beta;
    delete [] Check_Edge;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      delete NormalMeanK[iPoint];
    delete [] NormalMeanK;
    
  }
  
  /*--- Sharp edge detection is based in the statistical
   distribution of the curvature ---*/
  
  MaxK = K[0]; MinK = K[0]; MeanK = 0.0; TotalnPointDomain = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          MaxK = max(MaxK, fabs(K[iPoint]));
          MinK = min(MinK, fabs(K[iPoint]));
          MeanK += fabs(K[iPoint]);
          TotalnPointDomain++;
        }
      }
    }
  }
  
#ifdef HAVE_MPI
  double MyMeanK = MeanK; MeanK = 0.0;
  double MyMaxK = MaxK; MaxK = 0.0;
  unsigned long MynPointDomain = TotalnPointDomain; TotalnPointDomain = 0;
  MPI_Allreduce(&MyMeanK, &MeanK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&MyMaxK, &MaxK, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&MynPointDomain, &TotalnPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  /*--- Compute the mean ---*/
  MeanK /= double(TotalnPointDomain);
  
  /*--- Compute the standard deviation ---*/
  SigmaK = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          SigmaK += (fabs(K[iPoint]) - MeanK) * (fabs(K[iPoint]) - MeanK);
        }
      }
    }
  }
  
#ifdef HAVE_MPI
  double MySigmaK = SigmaK; SigmaK = 0.0;
  MPI_Allreduce(&MySigmaK, &SigmaK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  
  SigmaK = sqrt(SigmaK/double(TotalnPointDomain));
  
  if (rank == MASTER_NODE)
    cout << "Max K: " << MaxK << ". Mean K: " << MeanK << ". Standard deviation K: " << SigmaK << "." << endl;
  
  Point_Critical.clear();
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          if (fabs(K[iPoint]) > MeanK + config->GetRefSharpEdges()*SigmaK) {
            Point_Critical.push_back(iPoint);
          }
        }
      }
    }
  }
  
  /*--- Variables and buffers needed for MPI ---*/
  
#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#else
  nProcessor = 1;
#endif
  
  nLocalVertex = 0, nGlobalVertex = 0, MaxLocalVertex = 0;
  Buffer_Send_nVertex    = new unsigned long [1];
  Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  /*--- Count the total number of critical edge nodes. ---*/
  nLocalVertex = Point_Critical.size();
  Buffer_Send_nVertex[0] = nLocalVertex;
  
  /*--- Communicate to all processors the total number of critical edge nodes. ---*/
#ifdef HAVE_MPI
  MPI_Allreduce(&nLocalVertex, &nGlobalVertex,  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex, &MaxLocalVertex, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
  MaxLocalVertex = nLocalVertex;
  nGlobalVertex = nLocalVertex;
  Buffer_Receive_nVertex[0] = nLocalVertex;
#endif
  
  
  /*--- Create and initialize to zero some buffers to hold the coordinates
   of the boundary nodes that are communicated from each partition (all-to-all). ---*/
  
  Buffer_Send_Coord     = new double [MaxLocalVertex*nDim];
  Buffer_Receive_Coord  = new double [nProcessor*MaxLocalVertex*nDim];
  nBuffer               = MaxLocalVertex*nDim;
  
  for (iVertex = 0; iVertex < MaxLocalVertex; iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
  }
  
  /*--- Retrieve and store the coordinates of the sharp edges boundary nodes on
   the local partition and broadcast them to all partitions. ---*/
  
  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    iPoint = Point_Critical[iVertex];
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = node[iPoint]->GetCoord(iDim);
  }
  
#ifdef HAVE_MPI
  MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
#else
  for (iVertex = 0; iVertex < Point_Critical.size(); iVertex++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Receive_Coord[iVertex*nDim+iDim] = Buffer_Send_Coord[iVertex*nDim+iDim];
    }
  }
#endif
  
  /*--- Loop over all interior mesh nodes on the local partition and compute
   the distances to each of the no-slip boundary nodes in the entire mesh.
   Store the minimum distance to the wall for each interior mesh node. ---*/
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    Coord = node[iPoint]->GetCoord();
    
    MinDist = 1E20;
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
        Dist = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Dist += (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim])*
          (Coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex+iVertex)*nDim+iDim]);
        }
        Dist = sqrt(Dist);
        if (Dist < MinDist) MinDist = Dist;
      }
    }
    node[iPoint]->SetSharpEdge_Distance(MinDist);
  }
  
  /*--- Deallocate Max curvature ---*/
  delete[] K;
  
  /*--- Deallocate the buffers needed for the MPI communication. ---*/
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
  
}

CPhysicalGeometry::CPhysicalGeometry() : CGeometry() {}

CPhysicalGeometry::CPhysicalGeometry(CConfig *config, unsigned short val_iZone, unsigned short val_nZone) : CGeometry() {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short iNode_Surface, iMarker, iDim;
  unsigned long Point_Surface, iElem_Surface, iPoint;
  double Mesh_Scale_Change = 1.0, NewCoord[3];
  int rank = MASTER_NODE;
  nZone = val_nZone;
  
  string val_mesh_filename = config->GetMesh_FileName();
  unsigned short val_format = config->GetMesh_FileFormat();
  
  /*--- Initialize counters for local/global points & elements ---*/
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl <<"---------------------- Read grid file information -----------------------" << endl;
  
  switch (val_format) {
    case SU2:
      Read_SU2_Format(config, val_mesh_filename, val_iZone, val_nZone);
      break;
    case CGNS:
      Read_CGNS_Format(config, val_mesh_filename, val_iZone, val_nZone);
      break;
    case NETCDF_ASCII:
      Read_NETCDF_Format(config, val_mesh_filename, val_iZone, val_nZone);
      break;
    default:
      cout << "Unrecognized mesh format specified!!" << endl;
#ifndef HAVE_MPI
      exit(EXIT_FAILURE);
#else
      MPI_Abort(MPI_COMM_WORLD,1);
      MPI_Finalize();
#endif
      break;
  }
  
  /*--- Loop over the surface element to set the boundaries ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL &&
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX &&
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[Point_Surface]->SetSolidBoundary(true);
      }
  
  /*--- Loop over the points element to re-scale the mesh, and plot it ---*/
  
  if (config->GetKind_SU2() == SU2_CFD) {
    
    Mesh_Scale_Change = config->GetMesh_Scale_Change();
    
    /*--- Change the scale of the mesh ---*/
    
    if (Mesh_Scale_Change != 1.0) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          NewCoord[iDim] = Mesh_Scale_Change*node[iPoint]->GetCoord(iDim);
        }
        node[iPoint]->SetCoord(NewCoord);
      }
    }
    
    /*--- Output the grid using SU2 format ---*/

    if (config->GetMesh_Output()) {
      SetMeshFile(config, config->GetMesh_Out_FileName());
      cout.precision(4);
      cout << "Scaled mesh by a factor of " << Mesh_Scale_Change << endl;
      cout << " and wrote to the output file: " << config->GetMesh_Out_FileName() << endl;
    }
    
    /*--- The US system uses feet, but SU2 assumes that the grid is in inches ---*/
    
    if (config->GetSystemMeasurements() == US) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          NewCoord[iDim] = node[iPoint]->GetCoord(iDim)/12.0;
        }
        node[iPoint]->SetCoord(NewCoord);
      }
    }
    
  }
  
}


CPhysicalGeometry::CPhysicalGeometry(CGeometry *geometry, CConfig *config) {
  
  unsigned long iter,  iPoint, jPoint, iElem, jElem, iVertex;
  unsigned long nElemTotal = 0, nPointTotal = 0, nPointDomainTotal = 0, nPointGhost = 0, nPointPeriodic = 0, nElemTriangle = 0, nElemRectangle = 0, nElemTetrahedron = 0, nElemHexahedron = 0, nElemWedge = 0, nElemPyramid = 0;
  unsigned long iElemTotal, iPointTotal, iPointGhost, iPointDomain, iPointPeriodic, iElemTriangle, iElemRectangle, iElemTetrahedron, iElemHexahedron, iElemWedge, iElemPyramid;
  unsigned long nBoundLineTotal = 0, iBoundLineTotal;
  unsigned long nBoundTriangleTotal = 0, iBoundTriangleTotal;
  unsigned long nBoundRectangleTotal = 0, iBoundRectangleTotal;
  unsigned long ReceptorColor = 0, DonorColor = 0, Transformation;
  unsigned long *nElem_Color = NULL, **Elem_Color = NULL, Max_nElem_Color = 0;
  unsigned long nTotalSendDomain_Periodic = 0, iTotalSendDomain_Periodic, nTotalReceivedDomain_Periodic = 0, iTotalReceivedDomain_Periodic, *nSendDomain_Periodic = NULL, *nReceivedDomain_Periodic = NULL;
  unsigned long Buffer_Send_nPointTotal = 0, Buffer_Send_nPointDomainTotal = 0, Buffer_Send_nPointGhost = 0, Buffer_Send_nPointPeriodic = 0;
  unsigned long Buffer_Send_nElemTotal, Buffer_Send_nElemTriangle = 0, Buffer_Send_nElemRectangle = 0, Buffer_Send_nElemTetrahedron = 0, Buffer_Send_nElemHexahedron = 0, Buffer_Send_nElemWedge = 0, Buffer_Send_nElemPyramid = 0;
  unsigned long Buffer_Send_nTotalSendDomain_Periodic = 0, Buffer_Send_nTotalReceivedDomain_Periodic = 0, *Buffer_Send_nSendDomain_Periodic = NULL, *Buffer_Send_nReceivedDomain_Periodic = NULL;
  unsigned long Buffer_Send_nBoundLineTotal = 0, Buffer_Send_nBoundTriangleTotal = 0, Buffer_Send_nBoundRectangleTotal = 0;
  unsigned long iVertexDomain, iBoundLine, iBoundTriangle, iBoundRectangle;
  
  /*--- Need to double-check these shorts in case we go to nprocs > ~32,000 ---*/
  unsigned short iNode, iDim, iMarker, jMarker, nMarkerDomain = 0, iMarkerDomain;
  unsigned short nDomain = 0, iDomain, jDomain, nPeriodic = 0, iPeriodic, overhead = 4, Buffer_Send_nMarkerDomain = 0, Buffer_Send_nDim = 0, Buffer_Send_nZone = 0, Buffer_Send_nPeriodic = 0;
  
  bool *MarkerIn = NULL, **VertexIn = NULL, CheckDomain;
  long vnodes_local[8], *Global_to_Local_Point = NULL;
  char Marker_All_TagBound[MAX_NUMBER_MARKER][MAX_STRING_SIZE], Buffer_Send_Marker_All_TagBound[MAX_NUMBER_MARKER][MAX_STRING_SIZE];
  vector<long> DomainList;
  short *Marker_All_SendRecv_Copy = NULL;
  string *Marker_All_TagBound_Copy = NULL;
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
  /*--- Some dynamic arrays so we're not allocating too much on the stack ---*/
  unsigned long *nVertexDomain = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *nBoundLine = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *nBoundTriangle = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *nBoundRectangle = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *Buffer_Send_nVertexDomain = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *Buffer_Send_nBoundLine = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *Buffer_Send_nBoundTriangle = new unsigned long[MAX_NUMBER_MARKER];
  unsigned long *Buffer_Send_nBoundRectangle = new unsigned long[MAX_NUMBER_MARKER];
  unsigned short *Buffer_Send_Marker_All_SendRecv = new unsigned short[MAX_NUMBER_MARKER];
  
#ifdef HAVE_MPI
  
  /*--- MPI initialization ---*/
  
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  /*--- MPI status and request arrays for non-blocking communications ---*/
  
  MPI_Status send_stat[31], recv_stat[31];
  MPI_Request send_req[31], recv_req[31];
  
#endif
  
  /*--- Define buffer vector interior domain ---*/
  
  double        *Buffer_Send_Coord = NULL,            *Buffer_Receive_Coord = NULL;
  unsigned long *Buffer_Send_Color = NULL,            *Buffer_Receive_Color = NULL;
  unsigned long *Buffer_Send_GlobalPointIndex = NULL, *Buffer_Receive_GlobalPointIndex = NULL;
  unsigned long *Buffer_Send_Triangle = NULL,         *Buffer_Receive_Triangle = NULL;
  unsigned long *Buffer_Send_Rectangle = NULL,        *Buffer_Receive_Rectangle = NULL;
  unsigned long *Buffer_Send_Tetrahedron = NULL,      *Buffer_Receive_Tetrahedron = NULL;
  unsigned long *Buffer_Send_Hexahedron = NULL,       *Buffer_Receive_Hexahedron = NULL;
  unsigned long *Buffer_Send_Wedge = NULL,            *Buffer_Receive_Wedge = NULL;
  unsigned long *Buffer_Send_Pyramid = NULL,          *Buffer_Receive_Pyramid = NULL;
  
  /*--- Define buffer vector boundary ---*/
  
  unsigned long *Buffer_Send_BoundLine = NULL,            *Buffer_Receive_BoundLine = NULL;
  unsigned long *Buffer_Send_BoundTriangle = NULL,        *Buffer_Receive_BoundTriangle = NULL;
  unsigned long *Buffer_Send_BoundRectangle = NULL,       *Buffer_Receive_BoundRectangle = NULL;
  unsigned long *Buffer_Send_Local2Global_Marker = NULL,  *Buffer_Receive_Local2Global_Marker = NULL;
  
  /*--- Define buffer vector periodic boundary conditions ---*/
  
  double *Buffer_Send_Center = NULL,    *Buffer_Receive_Center = NULL;
  double *Buffer_Send_Rotation = NULL,  *Buffer_Receive_Rotation = NULL;
  double *Buffer_Send_Translate = NULL, *Buffer_Receive_Translate = NULL;
  
  /*--- Define buffer vector periodic boundary conditions ---*/
  
  unsigned long *Buffer_Send_SendDomain_Periodic = NULL,          *Buffer_Receive_SendDomain_Periodic = NULL;
  unsigned long *Buffer_Send_SendDomain_PeriodicTrans = NULL,     *Buffer_Receive_SendDomain_PeriodicTrans = NULL;
  unsigned long *Buffer_Send_SendDomain_PeriodicReceptor = NULL,  *Buffer_Receive_SendDomain_PeriodicReceptor = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_Periodic = NULL,      *Buffer_Receive_ReceivedDomain_Periodic = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_PeriodicTrans = NULL, *Buffer_Receive_ReceivedDomain_PeriodicTrans = NULL;
  unsigned long *Buffer_Send_ReceivedDomain_PeriodicDonor = NULL, *Buffer_Receive_ReceivedDomain_PeriodicDonor = NULL;
  
  
  /*--- Basic dimensionalization ---*/
  
  nDomain = size;
  Marker_All_SendRecv = new short[MAX_NUMBER_MARKER];
  nSendDomain_Periodic = new unsigned long [nDomain];
  nReceivedDomain_Periodic = new unsigned long [nDomain];
  
  /*--- Auxiliar vector based on the original geometry ---*/
  
  if (rank == MASTER_NODE) {
    
    MarkerIn = new bool [geometry->GetnMarker()];
    
    VertexIn = new bool* [geometry->GetnMarker()];
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
      VertexIn[iMarker] = new bool [geometry->GetnElem_Bound(iMarker)];
    
    Global_to_Local_Point =  new long[geometry->GetnPoint()];
    
    Buffer_Send_nDim = geometry->GetnDim();
    Buffer_Send_nZone = geometry->GetnZone();
    
    Buffer_Send_nPeriodic = config->GetnPeriodicIndex();
    Buffer_Send_Center = new double[Buffer_Send_nPeriodic*3];
    Buffer_Send_Rotation  = new double[Buffer_Send_nPeriodic*3];
    Buffer_Send_Translate = new double[Buffer_Send_nPeriodic*3];
    
    Buffer_Send_nSendDomain_Periodic = new unsigned long [nDomain];
    Buffer_Send_nReceivedDomain_Periodic = new unsigned long [nDomain];
    
    /*--- Divide the elements in color list to speed up the grid partitioning ---*/
    
    nElem_Color = new unsigned long[nDomain];
    for (iDomain = 0; iDomain < nDomain; iDomain++) nElem_Color[iDomain] = 0;
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      DomainList.clear();
      for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        iDomain = geometry->node[iPoint]->GetColor();
        
        CheckDomain = true;
        for (jDomain = 0; jDomain < DomainList.size(); jDomain++) {
          if (DomainList[jDomain] == iDomain) { CheckDomain = false; break; }
        }
        
        /*--- If the element is not in the list, then add it ---*/
        if (CheckDomain) {
          DomainList.push_back(iDomain);
          nElem_Color[iDomain]++;
        }
        
      }
    }
    
    /*--- Find the maximum number of elements per color to allocate the list ---*/
    
    Max_nElem_Color = 0;
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
      if (nElem_Color[iDomain] > Max_nElem_Color) Max_nElem_Color = nElem_Color[iDomain];
    }
    
    /*--- Allocate the element color array ---*/
    
    Elem_Color = new unsigned long* [nDomain];
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
      Elem_Color[iDomain] =  new unsigned long[Max_nElem_Color];
      nElem_Color[iDomain] = 0;
    }
    
    /*--- Create the lement list based on the color ---*/
    
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      
      DomainList.clear();
      for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
        iPoint = geometry->elem[iElem]->GetNode(iNode);
        iDomain = geometry->node[iPoint]->GetColor();
        
        /*--- Check if the element has been already added to the color ---*/
        
        CheckDomain = true;
        for (jDomain = 0; jDomain < DomainList.size(); jDomain++) {
          if (DomainList[jDomain] == iDomain) { CheckDomain = false; break; }
        }
        
        if (CheckDomain) {
          DomainList.push_back(iDomain);
          Elem_Color[iDomain][nElem_Color[iDomain]] = iElem;
          nElem_Color[iDomain]++;
        }
        
      }
    }
    
    
    /*--- Create a local copy of config->GetMarker_All_SendRecv and
     config->GetMarker_All_TagBound in the master node ---*/
    
    Marker_All_SendRecv_Copy = new short [geometry->GetnMarker()];
    Marker_All_TagBound_Copy   = new string[geometry->GetnMarker()];
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      Marker_All_SendRecv_Copy[iMarker] = config->GetMarker_All_SendRecv(iMarker);
      Marker_All_TagBound_Copy[iMarker] = config->GetMarker_All_TagBound(iMarker);
    }
    
    
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    if (rank == MASTER_NODE) {

      /*--- Interior dimensionalization. Loop over the original grid to perform the
       dimensionalizaton of the domain variables ---*/
      
      Buffer_Send_nElemTotal = 0; Buffer_Send_nPointTotal = 0; Buffer_Send_nPointGhost = 0; Buffer_Send_nPointDomainTotal = 0; Buffer_Send_nPointPeriodic = 0;
      Buffer_Send_nElemTriangle = 0; Buffer_Send_nElemRectangle = 0; Buffer_Send_nElemTetrahedron = 0; Buffer_Send_nElemHexahedron = 0; Buffer_Send_nElemWedge = 0; Buffer_Send_nElemPyramid = 0;
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Global_to_Local_Point[iPoint] = -1;
      
      for (jElem = 0; jElem < nElem_Color[iDomain]; jElem++) {
        
        iElem = Elem_Color[iDomain][jElem];
        
        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = geometry->elem[iElem]->GetNode(iNode);
          if (Global_to_Local_Point[iPoint] == -1) {
            Global_to_Local_Point[iPoint] = 1;
            Buffer_Send_nPointTotal++;
            if ( geometry->node[iPoint]->GetColor() != iDomain ) Buffer_Send_nPointGhost++;
            else {
              if (iPoint > geometry->GetnPointDomain() - 1) { Buffer_Send_nPointGhost++;  Buffer_Send_nPointPeriodic++; }
              else Buffer_Send_nPointDomainTotal++;
            }
          }
        }
        
        switch(geometry->elem[iElem]->GetVTK_Type()) {
          case TRIANGLE: Buffer_Send_nElemTriangle++; break;
          case RECTANGLE: Buffer_Send_nElemRectangle++; break;
          case TETRAHEDRON: Buffer_Send_nElemTetrahedron++; break;
          case HEXAHEDRON: Buffer_Send_nElemHexahedron++; break;
          case WEDGE: Buffer_Send_nElemWedge++; break;
          case PYRAMID: Buffer_Send_nElemPyramid++; break;
        }
        Buffer_Send_nElemTotal++;
        
      }
      
      /*--- Boundary dimensionalization. Dimensionalization with physical boundaries, compute Buffer_Send_nMarkerDomain,
       Buffer_Send_nVertexDomain[nMarkerDomain] ---*/
      
      Buffer_Send_nMarkerDomain = 0; Buffer_Send_nBoundLineTotal = 0; Buffer_Send_nBoundTriangleTotal = 0; Buffer_Send_nBoundRectangleTotal = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        Buffer_Send_nVertexDomain[iMarker] = 0;
        Buffer_Send_nBoundLine[iMarker] = 0;
        Buffer_Send_nBoundTriangle[iMarker] = 0;
        Buffer_Send_nBoundRectangle[iMarker] = 0;
        
        Buffer_Send_Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
        sprintf(Buffer_Send_Marker_All_TagBound[iMarker], "%s", Marker_All_TagBound_Copy[iMarker].c_str());
        
      }
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) {
          MarkerIn[iMarker] = false; Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] = 0;
          
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            VertexIn[iMarker][iVertex] = false;
            for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
              iPoint = geometry->bound[iMarker][iVertex]->GetNode(iNode);
              if (geometry->node[iPoint]->GetColor() == iDomain) VertexIn[iMarker][iVertex] = true;
            }
            
            if (VertexIn[iMarker][iVertex]) {
              switch(geometry->bound[iMarker][iVertex]->GetVTK_Type()) {
                case LINE: Buffer_Send_nBoundLine[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundLineTotal++; break;
                case TRIANGLE: Buffer_Send_nBoundTriangle[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundTriangleTotal++; break;
                case RECTANGLE: Buffer_Send_nBoundRectangle[Buffer_Send_nMarkerDomain]++; Buffer_Send_nBoundRectangleTotal++; break;
              }
              
              Buffer_Send_nVertexDomain[Buffer_Send_nMarkerDomain] ++;
              MarkerIn[iMarker] = true;
              
            }
          }
          
          if (MarkerIn[iMarker]) { Buffer_Send_nMarkerDomain++; }
          
        }
      }
      
      /*--- Copy periodic information from the config file ---*/
      
      for (iPeriodic = 0; iPeriodic < Buffer_Send_nPeriodic; iPeriodic++) {
        for (iDim = 0; iDim < 3; iDim++) {
          Buffer_Send_Center[iDim+iPeriodic*3] = config->GetPeriodicCenter(iPeriodic)[iDim];
          Buffer_Send_Rotation[iDim+iPeriodic*3] = config->GetPeriodicRotation(iPeriodic)[iDim];
          Buffer_Send_Translate[iDim+iPeriodic*3] = config->GetPeriodicTranslate(iPeriodic)[iDim];
        }
      }
      
      /*--- Dimensionalization of the periodic auxiliar vectors ---*/
      
      for (jDomain = 0; jDomain < nDomain; jDomain++) {
        Buffer_Send_nSendDomain_Periodic[jDomain] = 0;
        Buffer_Send_nReceivedDomain_Periodic[jDomain] = 0;
      }
      Buffer_Send_nTotalSendDomain_Periodic = 0;
      Buffer_Send_nTotalReceivedDomain_Periodic = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
            if (iDomain == geometry->node[iPoint]->GetColor()) {
              
              if (config->GetMarker_All_SendRecv(iMarker) > 0) {
                
                /*--- Identify the color of the receptor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    ReceptorColor = geometry->node[jPoint]->GetColor();
                  }
                }
                
                Buffer_Send_nSendDomain_Periodic[ReceptorColor]++;
                Buffer_Send_nTotalSendDomain_Periodic++;
                
              }
              if (config->GetMarker_All_SendRecv(iMarker) < 0) {
                
                /*--- Identify the color of the donor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    DonorColor = geometry->node[jPoint]->GetColor();
                  }
                }
                
                Buffer_Send_nReceivedDomain_Periodic[DonorColor]++;
                Buffer_Send_nTotalReceivedDomain_Periodic++;
                
              }
            }
          }
        }
      }
      
      /*--- Allocate the buffer vectors in the appropiate domain (master, iDomain) ---*/
      
      Buffer_Send_Coord =             new double [Buffer_Send_nPointTotal*Buffer_Send_nDim];
      Buffer_Send_Color =             new unsigned long [Buffer_Send_nPointTotal];
      Buffer_Send_GlobalPointIndex =  new unsigned long [Buffer_Send_nPointTotal];
      Buffer_Send_Triangle =          new unsigned long [Buffer_Send_nElemTriangle*3];
      Buffer_Send_Rectangle =         new unsigned long [Buffer_Send_nElemRectangle*4];
      Buffer_Send_Tetrahedron =       new unsigned long [Buffer_Send_nElemTetrahedron*4];
      Buffer_Send_Hexahedron =        new unsigned long [Buffer_Send_nElemHexahedron*8];
      Buffer_Send_Wedge =             new unsigned long [Buffer_Send_nElemWedge*6];
      Buffer_Send_Pyramid =           new unsigned long [Buffer_Send_nElemPyramid*5];
      
      Buffer_Send_BoundLine =           new unsigned long [Buffer_Send_nBoundLineTotal*2];
      Buffer_Send_BoundTriangle =       new unsigned long [Buffer_Send_nBoundTriangleTotal*3];
      Buffer_Send_BoundRectangle =      new unsigned long [Buffer_Send_nBoundRectangleTotal*4];
      Buffer_Send_Local2Global_Marker = new unsigned long [Buffer_Send_nMarkerDomain];
      
      Buffer_Send_SendDomain_Periodic           = new unsigned long [Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_SendDomain_PeriodicTrans      = new unsigned long [Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_SendDomain_PeriodicReceptor   = new unsigned long [Buffer_Send_nTotalSendDomain_Periodic];
      Buffer_Send_ReceivedDomain_Periodic       = new unsigned long [Buffer_Send_nTotalReceivedDomain_Periodic];
      Buffer_Send_ReceivedDomain_PeriodicTrans  = new unsigned long [Buffer_Send_nTotalReceivedDomain_Periodic];
      Buffer_Send_ReceivedDomain_PeriodicDonor  = new unsigned long [Buffer_Send_nTotalReceivedDomain_Periodic];

      if (iDomain != MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        MPI_Isend(&Buffer_Send_nDim,               1,  MPI_UNSIGNED_SHORT,  iDomain, 0, MPI_COMM_WORLD, &send_req[0]);
        MPI_Isend(&Buffer_Send_nZone,              1,  MPI_UNSIGNED_SHORT,  iDomain, 1, MPI_COMM_WORLD, &send_req[1]);
        MPI_Isend(&Buffer_Send_nPointTotal,        1,  MPI_UNSIGNED_LONG,   iDomain, 2, MPI_COMM_WORLD, &send_req[2]);
        MPI_Isend(&Buffer_Send_nPointDomainTotal,  1,  MPI_UNSIGNED_LONG,   iDomain, 3, MPI_COMM_WORLD, &send_req[3]);
        MPI_Isend(&Buffer_Send_nPointGhost,        1,  MPI_UNSIGNED_LONG,   iDomain, 4, MPI_COMM_WORLD, &send_req[4]);
        MPI_Isend(&Buffer_Send_nPointPeriodic,     1,  MPI_UNSIGNED_LONG,   iDomain, 5, MPI_COMM_WORLD, &send_req[5]);
        MPI_Isend(&Buffer_Send_nElemTotal,         1,  MPI_UNSIGNED_LONG,   iDomain, 6, MPI_COMM_WORLD, &send_req[6]);
        MPI_Isend(&Buffer_Send_nElemTriangle,      1,  MPI_UNSIGNED_LONG,   iDomain, 7, MPI_COMM_WORLD, &send_req[7]);
        MPI_Isend(&Buffer_Send_nElemRectangle,     1,  MPI_UNSIGNED_LONG,   iDomain, 8, MPI_COMM_WORLD, &send_req[8]);
        MPI_Isend(&Buffer_Send_nElemTetrahedron,   1,  MPI_UNSIGNED_LONG,   iDomain, 9, MPI_COMM_WORLD, &send_req[9]);
        MPI_Isend(&Buffer_Send_nElemHexahedron,    1,  MPI_UNSIGNED_LONG,   iDomain, 10, MPI_COMM_WORLD, &send_req[10]);
        MPI_Isend(&Buffer_Send_nElemWedge,         1,  MPI_UNSIGNED_LONG,   iDomain, 11, MPI_COMM_WORLD, &send_req[11]);
        MPI_Isend(&Buffer_Send_nElemPyramid,       1,  MPI_UNSIGNED_LONG,   iDomain, 12, MPI_COMM_WORLD, &send_req[12]);
        
        MPI_Isend(&Buffer_Send_nBoundLineTotal,                    1, MPI_UNSIGNED_LONG,  iDomain, 13, MPI_COMM_WORLD, &send_req[13]);
        MPI_Isend(&Buffer_Send_nBoundTriangleTotal,                1, MPI_UNSIGNED_LONG,  iDomain, 14, MPI_COMM_WORLD, &send_req[14]);
        MPI_Isend(&Buffer_Send_nBoundRectangleTotal,               1, MPI_UNSIGNED_LONG,  iDomain, 15, MPI_COMM_WORLD, &send_req[15]);
        MPI_Isend(&Buffer_Send_nMarkerDomain,                      1, MPI_UNSIGNED_SHORT, iDomain, 16, MPI_COMM_WORLD, &send_req[16]);
        MPI_Isend(Buffer_Send_nVertexDomain,       MAX_NUMBER_MARKER, MPI_UNSIGNED_LONG,  iDomain, 17, MPI_COMM_WORLD, &send_req[17]);
        MPI_Isend(Buffer_Send_nBoundLine,          MAX_NUMBER_MARKER, MPI_UNSIGNED_LONG,  iDomain, 18, MPI_COMM_WORLD, &send_req[18]);
        MPI_Isend(Buffer_Send_nBoundTriangle,      MAX_NUMBER_MARKER, MPI_UNSIGNED_LONG,  iDomain, 19, MPI_COMM_WORLD, &send_req[19]);
        MPI_Isend(Buffer_Send_nBoundRectangle,     MAX_NUMBER_MARKER, MPI_UNSIGNED_LONG,  iDomain, 20, MPI_COMM_WORLD, &send_req[20]);
        MPI_Isend(Buffer_Send_Marker_All_SendRecv, MAX_NUMBER_MARKER, MPI_SHORT,          iDomain, 21, MPI_COMM_WORLD, &send_req[21]);
        MPI_Isend(Buffer_Send_Marker_All_TagBound,  MAX_NUMBER_MARKER*MAX_STRING_SIZE, MPI_CHAR,           iDomain, 22, MPI_COMM_WORLD, &send_req[22]);
        
        MPI_Isend(&Buffer_Send_nPeriodic,              1, MPI_UNSIGNED_SHORT, iDomain, 23, MPI_COMM_WORLD, &send_req[23]);
        MPI_Isend(Buffer_Send_Center,    nPeriodic*3, MPI_DOUBLE, iDomain, 24, MPI_COMM_WORLD, &send_req[24]);
        MPI_Isend(Buffer_Send_Rotation,  nPeriodic*3, MPI_DOUBLE, iDomain, 25, MPI_COMM_WORLD, &send_req[25]);
        MPI_Isend(Buffer_Send_Translate, nPeriodic*3, MPI_DOUBLE, iDomain, 26, MPI_COMM_WORLD, &send_req[26]);
        
        MPI_Isend(&Buffer_Send_nTotalSendDomain_Periodic,           1, MPI_UNSIGNED_LONG, iDomain, 27, MPI_COMM_WORLD, &send_req[27]);
        MPI_Isend(&Buffer_Send_nTotalReceivedDomain_Periodic,       1, MPI_UNSIGNED_LONG, iDomain, 28, MPI_COMM_WORLD, &send_req[28]);
        MPI_Isend(Buffer_Send_nSendDomain_Periodic,           nDomain, MPI_UNSIGNED_LONG, iDomain, 29, MPI_COMM_WORLD, &send_req[29]);
        MPI_Isend(Buffer_Send_nReceivedDomain_Periodic,       nDomain, MPI_UNSIGNED_LONG, iDomain, 30, MPI_COMM_WORLD, &send_req[30]);
        
        /*--- Wait for this set of non-blocking comm. to complete ---*/
        
        MPI_Waitall(31, send_req, send_stat);
        
#endif
        
        
      } else {
        
        /*--- We are the master node, so simply copy values into place ---*/
        
        nDim = Buffer_Send_nDim;
        nZone = Buffer_Send_nZone;
        
        nPeriodic = Buffer_Send_nPeriodic;
        
        nPointTotal       = Buffer_Send_nPointTotal;
        nPointDomainTotal = Buffer_Send_nPointDomainTotal;
        nPointGhost       = Buffer_Send_nPointGhost;
        nPointPeriodic    = Buffer_Send_nPointPeriodic;
        
        nElemTotal        = Buffer_Send_nElemTotal;
        nElemTriangle     = Buffer_Send_nElemTriangle;
        nElemRectangle    = Buffer_Send_nElemRectangle;
        nElemTetrahedron  = Buffer_Send_nElemTetrahedron;
        nElemHexahedron   = Buffer_Send_nElemHexahedron;
        nElemWedge        = Buffer_Send_nElemWedge;
        nElemPyramid      = Buffer_Send_nElemPyramid;
        
        nelem_triangle = nElemTriangle;
        nelem_quad = nElemRectangle;
        nelem_tetra = nElemTetrahedron;
        nelem_hexa = nElemHexahedron;
        nelem_wedge = nElemWedge;
        nelem_pyramid = nElemPyramid;

        nBoundLineTotal      = Buffer_Send_nBoundLineTotal;
        nBoundTriangleTotal  = Buffer_Send_nBoundTriangleTotal;
        nBoundRectangleTotal = Buffer_Send_nBoundRectangleTotal;
        nMarkerDomain        = Buffer_Send_nMarkerDomain;
        
        for (iMarker = 0; iMarker < MAX_NUMBER_MARKER; iMarker++) {
          nVertexDomain[iMarker] = Buffer_Send_nVertexDomain[iMarker];
          nBoundLine[iMarker] = Buffer_Send_nBoundLine[iMarker];
          nBoundTriangle[iMarker] = Buffer_Send_nBoundTriangle[iMarker];
          nBoundRectangle[iMarker] = Buffer_Send_nBoundRectangle[iMarker];
          Marker_All_SendRecv[iMarker] = Buffer_Send_Marker_All_SendRecv[iMarker];
          for (iter = 0; iter < MAX_STRING_SIZE; iter++)
            Marker_All_TagBound[iMarker][iter] = Buffer_Send_Marker_All_TagBound[iMarker][iter];
        }
        
        Buffer_Receive_Center    = new double[nPeriodic*3];
        Buffer_Receive_Rotation  = new double[nPeriodic*3];
        Buffer_Receive_Translate = new double[nPeriodic*3];
        
        for (iter = 0; iter < nPeriodic*3; iter++) {
          Buffer_Receive_Center[iter] =  Buffer_Send_Center[iter];
          Buffer_Receive_Rotation[iter] =  Buffer_Send_Rotation[iter];
          Buffer_Receive_Translate[iter] =  Buffer_Send_Translate[iter];
        }
        
        nTotalSendDomain_Periodic = Buffer_Send_nTotalSendDomain_Periodic;
        nTotalReceivedDomain_Periodic = Buffer_Send_nTotalReceivedDomain_Periodic;
        
        for (iter = 0; iter < nDomain; iter++) {
          nSendDomain_Periodic[iter] = Buffer_Send_nSendDomain_Periodic[iter];
          nReceivedDomain_Periodic[iter] = Buffer_Send_nReceivedDomain_Periodic[iter];
        }
        
      }
    }
    
    /*--- Receive the size of buffers---*/
    
    if (rank == iDomain) {
      
      /*--- Receive the size of buffers---*/
      
      if (rank != MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        MPI_Irecv(&nDim,              1, MPI_UNSIGNED_SHORT, MASTER_NODE, 0,  MPI_COMM_WORLD, &recv_req[0]);
        MPI_Irecv(&nZone,             1, MPI_UNSIGNED_SHORT, MASTER_NODE, 1,  MPI_COMM_WORLD, &recv_req[1]);
        MPI_Irecv(&nPointTotal,       1, MPI_UNSIGNED_LONG,  MASTER_NODE, 2,  MPI_COMM_WORLD, &recv_req[2]);
        MPI_Irecv(&nPointDomainTotal, 1, MPI_UNSIGNED_LONG,  MASTER_NODE, 3,  MPI_COMM_WORLD, &recv_req[3]);
        MPI_Irecv(&nPointGhost,       1, MPI_UNSIGNED_LONG,  MASTER_NODE, 4,  MPI_COMM_WORLD, &recv_req[4]);
        MPI_Irecv(&nPointPeriodic,    1, MPI_UNSIGNED_LONG,  MASTER_NODE, 5,  MPI_COMM_WORLD, &recv_req[5]);
        MPI_Irecv(&nElemTotal,        1, MPI_UNSIGNED_LONG,  MASTER_NODE, 6,  MPI_COMM_WORLD, &recv_req[6]);
        MPI_Irecv(&nElemTriangle,     1, MPI_UNSIGNED_LONG,  MASTER_NODE, 7,  MPI_COMM_WORLD, &recv_req[7]);
        MPI_Irecv(&nElemRectangle,    1, MPI_UNSIGNED_LONG,  MASTER_NODE, 8,  MPI_COMM_WORLD, &recv_req[8]);
        MPI_Irecv(&nElemTetrahedron,  1, MPI_UNSIGNED_LONG,  MASTER_NODE, 9,  MPI_COMM_WORLD, &recv_req[9]);
        MPI_Irecv(&nElemHexahedron,   1, MPI_UNSIGNED_LONG,  MASTER_NODE, 10,  MPI_COMM_WORLD, &recv_req[10]);
        MPI_Irecv(&nElemWedge,        1, MPI_UNSIGNED_LONG,  MASTER_NODE, 11, MPI_COMM_WORLD, &recv_req[11]);
        MPI_Irecv(&nElemPyramid,      1, MPI_UNSIGNED_LONG,  MASTER_NODE, 12, MPI_COMM_WORLD, &recv_req[12]);
        
        MPI_Irecv(&nBoundLineTotal,                    1,  MPI_UNSIGNED_LONG, MASTER_NODE, 13, MPI_COMM_WORLD, &recv_req[13]);
        MPI_Irecv(&nBoundTriangleTotal,                1,  MPI_UNSIGNED_LONG, MASTER_NODE, 14, MPI_COMM_WORLD, &recv_req[14]);
        MPI_Irecv(&nBoundRectangleTotal,               1,  MPI_UNSIGNED_LONG, MASTER_NODE, 15, MPI_COMM_WORLD, &recv_req[15]);
        MPI_Irecv(&nMarkerDomain,                      1, MPI_UNSIGNED_SHORT, MASTER_NODE, 16, MPI_COMM_WORLD, &recv_req[16]);
        MPI_Irecv(nVertexDomain,       MAX_NUMBER_MARKER,  MPI_UNSIGNED_LONG, MASTER_NODE, 17, MPI_COMM_WORLD, &recv_req[17]);
        MPI_Irecv(nBoundLine,          MAX_NUMBER_MARKER,  MPI_UNSIGNED_LONG, MASTER_NODE, 18, MPI_COMM_WORLD, &recv_req[18]);
        MPI_Irecv(nBoundTriangle,      MAX_NUMBER_MARKER,  MPI_UNSIGNED_LONG, MASTER_NODE, 19, MPI_COMM_WORLD, &recv_req[19]);
        MPI_Irecv(nBoundRectangle,     MAX_NUMBER_MARKER,  MPI_UNSIGNED_LONG, MASTER_NODE, 20, MPI_COMM_WORLD, &recv_req[20]);
        MPI_Irecv(Marker_All_SendRecv, MAX_NUMBER_MARKER,          MPI_SHORT, MASTER_NODE, 21, MPI_COMM_WORLD, &recv_req[21]);
        MPI_Irecv(Marker_All_TagBound,      MAX_NUMBER_MARKER*MAX_STRING_SIZE,       MPI_CHAR, MASTER_NODE, 22, MPI_COMM_WORLD, &recv_req[22]);
        MPI_Irecv(&nPeriodic, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, 23, MPI_COMM_WORLD, &recv_req[23]);
        
        /*--- Wait for the this set of non-blocking recv's to complete ---*/
        
        MPI_Waitall(24, recv_req, recv_stat);
        
#endif
        
        /*--- Update the number of elements (local) ---*/

        nelem_triangle = nElemTriangle;
        nelem_quad = nElemRectangle;
        nelem_tetra = nElemTetrahedron;
        nelem_hexa = nElemHexahedron;
        nelem_wedge = nElemWedge;
        nelem_pyramid = nElemPyramid;
        
        /*--- Marker_All_TagBound and Marker_All_SendRecv, set the same values in the config files of all the files ---*/
        
        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
          config->SetMarker_All_TagBound(iMarker, string(Marker_All_TagBound[iMarker]));
        }
        
        /*--- Periodic boundary conditions, set the values in the config files of all the files ---*/
        
        Buffer_Receive_Center    = new double[nPeriodic*3];
        Buffer_Receive_Rotation  = new double[nPeriodic*3];
        Buffer_Receive_Translate = new double[nPeriodic*3];
        
#ifdef HAVE_MPI
        
        MPI_Irecv(Buffer_Receive_Center,    nPeriodic*3, MPI_DOUBLE, MASTER_NODE, 24, MPI_COMM_WORLD, &recv_req[0]);
        MPI_Irecv(Buffer_Receive_Rotation,  nPeriodic*3, MPI_DOUBLE, MASTER_NODE, 25, MPI_COMM_WORLD, &recv_req[1]);
        MPI_Irecv(Buffer_Receive_Translate, nPeriodic*3, MPI_DOUBLE, MASTER_NODE, 26, MPI_COMM_WORLD, &recv_req[2]);
        
        MPI_Irecv(&nTotalSendDomain_Periodic,       1, MPI_UNSIGNED_LONG, MASTER_NODE, 27, MPI_COMM_WORLD, &recv_req[3]);
        MPI_Irecv(&nTotalReceivedDomain_Periodic,   1, MPI_UNSIGNED_LONG, MASTER_NODE, 28, MPI_COMM_WORLD, &recv_req[4]);
        MPI_Irecv(nSendDomain_Periodic,             nDomain, MPI_UNSIGNED_LONG, MASTER_NODE, 29, MPI_COMM_WORLD, &recv_req[5]);
        MPI_Irecv(nReceivedDomain_Periodic,         nDomain, MPI_UNSIGNED_LONG, MASTER_NODE, 30, MPI_COMM_WORLD, &recv_req[6]);
        
        /*--- Wait for this set of non-blocking comm. to complete ---*/
        
        MPI_Waitall(7, recv_req, recv_stat);
        
#endif
        
        config->SetnPeriodicIndex(nPeriodic);
        
        for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
          
          double* center = new double[3];       // Do not deallocate the memory
          double* rotation  = new double[3];    // Do not deallocate the memory
          double* translate = new double[3];    // Do not deallocate the memory
          
          for (iDim = 0; iDim < 3; iDim++) {
            center[iDim] = Buffer_Receive_Center[iDim+iPeriodic*3];
            rotation[iDim] = Buffer_Receive_Rotation[iDim+iPeriodic*3];
            translate[iDim] = Buffer_Receive_Translate[iDim+iPeriodic*3];
          }
          config->SetPeriodicCenter(iPeriodic, center);
          config->SetPeriodicRotation(iPeriodic, rotation);
          config->SetPeriodicTranslate(iPeriodic, translate);
        }
        
      }
      
      delete [] Buffer_Receive_Center;
      delete [] Buffer_Receive_Rotation;
      delete [] Buffer_Receive_Translate;
      
      /*--- Allocate the receive buffer vector ---*/
      
      Buffer_Receive_Coord =              new double [nPointTotal*nDim];
      Buffer_Receive_Color =              new unsigned long [nPointTotal];
      Buffer_Receive_GlobalPointIndex =   new unsigned long [nPointTotal];
      Buffer_Receive_Triangle =           new unsigned long [nElemTriangle*3];
      Buffer_Receive_Rectangle =          new unsigned long [nElemRectangle*4];
      Buffer_Receive_Tetrahedron =        new unsigned long [nElemTetrahedron*4];
      Buffer_Receive_Hexahedron =         new unsigned long [nElemHexahedron*8];
      Buffer_Receive_Wedge =              new unsigned long [nElemWedge*6];
      Buffer_Receive_Pyramid =            new unsigned long [nElemPyramid*5];
      Buffer_Receive_BoundLine =            new unsigned long [nBoundLineTotal*2];
      Buffer_Receive_BoundTriangle =        new unsigned long [nBoundTriangleTotal*3];
      Buffer_Receive_BoundRectangle =       new unsigned long [nBoundRectangleTotal*4];
      Buffer_Receive_Local2Global_Marker =  new unsigned long [nMarkerDomain];
      
      Buffer_Receive_SendDomain_Periodic           = new unsigned long [nTotalSendDomain_Periodic];
      Buffer_Receive_SendDomain_PeriodicTrans      = new unsigned long [nTotalSendDomain_Periodic];
      Buffer_Receive_SendDomain_PeriodicReceptor   = new unsigned long [nTotalSendDomain_Periodic];
      Buffer_Receive_ReceivedDomain_Periodic       = new unsigned long [nTotalReceivedDomain_Periodic];
      Buffer_Receive_ReceivedDomain_PeriodicTrans  = new unsigned long [nTotalReceivedDomain_Periodic];
      Buffer_Receive_ReceivedDomain_PeriodicDonor  = new unsigned long [nTotalReceivedDomain_Periodic];
      
    }
    
    /*--- Set the value of the Send buffers ---*/
    
    if (rank == MASTER_NODE) {
      
      /*--- Set the value of the interior geometry ---*/
      
      iElemTotal = 0; iPointTotal = 0; iPointDomain = 0; iPointPeriodic = Buffer_Send_nPointDomainTotal; iPointGhost = Buffer_Send_nPointDomainTotal + Buffer_Send_nPointPeriodic;
      iElemTriangle = 0; iElemRectangle = 0; iElemTetrahedron = 0; iElemHexahedron = 0; iElemWedge = 0; iElemPyramid = 0;
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) Global_to_Local_Point[iPoint] = -1;
      
      for (jElem = 0; jElem < nElem_Color[iDomain]; jElem++) {
        
        iElem = Elem_Color[iDomain][jElem];
        
        for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = geometry->elem[iElem]->GetNode(iNode);
          if (Global_to_Local_Point[iPoint] == -1) {
            
            if ( geometry->node[iPoint]->GetColor() == iDomain ) {
              if ( iPoint > geometry->GetnPointDomain() - 1) iPointTotal = iPointPeriodic;
              else iPointTotal = iPointDomain;
            }
            else iPointTotal = iPointGhost;
            
            Global_to_Local_Point[iPoint] = iPointTotal;
            Buffer_Send_Color[iPointTotal] = geometry->node[iPoint]->GetColor();
            Buffer_Send_GlobalPointIndex[iPointTotal] = iPoint;
            for (iDim = 0; iDim < Buffer_Send_nDim; iDim++)
              Buffer_Send_Coord[Buffer_Send_nDim*iPointTotal+iDim] = geometry->node[iPoint]->GetCoord(iDim);
            
            if ( geometry->node[iPoint]->GetColor() == iDomain ) {
              if ( iPoint > geometry->GetnPointDomain() - 1) iPointPeriodic++;
              else iPointDomain++;
            }
            else iPointGhost++;
            
          }
          
          vnodes_local[iNode] = Global_to_Local_Point[iPoint];
          
        }
        
        switch(geometry->elem[iElem]->GetVTK_Type()) {
          case TRIANGLE:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Triangle[3*iElemTriangle+iNode] = vnodes_local[iNode];
            iElemTriangle++; break;
          case RECTANGLE:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Rectangle[4*iElemRectangle+iNode] = vnodes_local[iNode];
            iElemRectangle++; break;
          case TETRAHEDRON:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Tetrahedron[4*iElemTetrahedron+iNode] = vnodes_local[iNode];
            iElemTetrahedron++; break;
          case HEXAHEDRON:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Hexahedron[8*iElemHexahedron+iNode] = vnodes_local[iNode];
            iElemHexahedron++; break;
          case WEDGE:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Wedge[6*iElemWedge+iNode] = vnodes_local[iNode];
            iElemWedge++; break;
          case PYRAMID:
            for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++)
              Buffer_Send_Pyramid[5*iElemPyramid+iNode] = vnodes_local[iNode];
            iElemPyramid++; break;
        }
        
        iElemTotal++;
        
      }
      
      /*--- Set the value of the boundary geometry ---*/
      
      iMarkerDomain = 0;
      iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundRectangleTotal = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) && (MarkerIn[iMarker])) {
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            
            if (VertexIn[iMarker][iVertex]) {
              
              for (iNode = 0; iNode < geometry->bound[iMarker][iVertex]->GetnNodes(); iNode++) {
                vnodes_local[iNode] = Global_to_Local_Point[geometry->bound[iMarker][iVertex]->GetNode(iNode)];
              }
              
              switch(geometry->bound[iMarker][iVertex]->GetVTK_Type()) {
                case LINE:
                  Buffer_Send_BoundLine[2*iBoundLineTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundLine[2*iBoundLineTotal+1] = vnodes_local[1];
                  iBoundLineTotal++;
                  break;
                case TRIANGLE:
                  Buffer_Send_BoundTriangle[3*iBoundTriangleTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundTriangle[3*iBoundTriangleTotal+1] = vnodes_local[1];
                  Buffer_Send_BoundTriangle[3*iBoundTriangleTotal+2] = vnodes_local[2];
                  iBoundTriangleTotal++;
                  break;
                case RECTANGLE:
                  Buffer_Send_BoundRectangle[4*iBoundRectangleTotal+0] = vnodes_local[0];
                  Buffer_Send_BoundRectangle[4*iBoundRectangleTotal+1] = vnodes_local[1];
                  Buffer_Send_BoundRectangle[4*iBoundRectangleTotal+2] = vnodes_local[2];
                  Buffer_Send_BoundRectangle[4*iBoundRectangleTotal+3] = vnodes_local[3];
                  iBoundRectangleTotal++;
                  break;
              }
            }
          }
          
          Buffer_Send_Local2Global_Marker[iMarkerDomain] = iMarker;
          iMarkerDomain++;
          
        }
      }
      
      /*--- Evaluate the number of already existing periodic boundary conditions ---*/
      
      iTotalSendDomain_Periodic = 0;
      iTotalReceivedDomain_Periodic = 0;
      
      for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
        
        if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
          
          for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
            
            iPoint = geometry->bound[iMarker][iVertex]->GetNode(0);
            Transformation = geometry->bound[iMarker][iVertex]->GetRotation_Type();
            
            if (iDomain == geometry->node[iPoint]->GetColor()) {
              
              /*--- If the information is going to be sended, find the
               domain of the receptor ---*/
              
              if (config->GetMarker_All_SendRecv(iMarker) > 0) {
                
                /*--- Identify the color of the receptor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker))) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    ReceptorColor = geometry->node[jPoint]->GetColor();
                  }
                }
                
                /*--- For each color of the receptor we will han an extra marker (+) ---*/
                
                Buffer_Send_SendDomain_Periodic[iTotalSendDomain_Periodic] = Global_to_Local_Point[iPoint];
                Buffer_Send_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic] = Transformation;
                Buffer_Send_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] = ReceptorColor;
                
                iTotalSendDomain_Periodic++;
                
              }
              
              /*--- If the information is goint to be received, find the domain if the donor ---*/
              
              if (config->GetMarker_All_SendRecv(iMarker) < 0) {
                
                /*--- Identify the color of the donor ---*/
                
                for (jMarker = 0; jMarker < geometry->GetnMarker(); jMarker++) {
                  if ((config->GetMarker_All_KindBC(jMarker) == SEND_RECEIVE) &&
                      (config->GetMarker_All_SendRecv(jMarker) == -config->GetMarker_All_SendRecv(iMarker) )) {
                    jPoint = geometry->bound[jMarker][iVertex]->GetNode(0);
                    DonorColor = geometry->node[jPoint]->GetColor();
                  }
                }
                
                /*--- For each color of the donor we will han an extra marker (-) ---*/
                
                Buffer_Send_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic] = Global_to_Local_Point[iPoint];
                Buffer_Send_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic] = Transformation;
                Buffer_Send_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] = DonorColor;
                
                iTotalReceivedDomain_Periodic++;
                
              }
            }
          }
        }
      }
      
      /*--- Send the buffers with the geometrical information ---*/
      
      if (iDomain != MASTER_NODE) {
        
#ifdef HAVE_MPI
        
        MPI_Isend(Buffer_Send_Coord,               Buffer_Send_nPointTotal*Buffer_Send_nDim,  MPI_DOUBLE, iDomain, 0,  MPI_COMM_WORLD, &send_req[0]);
        MPI_Isend(Buffer_Send_GlobalPointIndex,    Buffer_Send_nPointTotal,            MPI_UNSIGNED_LONG, iDomain, 1,  MPI_COMM_WORLD, &send_req[1]);
        MPI_Isend(Buffer_Send_Color,               Buffer_Send_nPointTotal,            MPI_UNSIGNED_LONG, iDomain, 2,  MPI_COMM_WORLD, &send_req[2]);
        MPI_Isend(Buffer_Send_Triangle,            Buffer_Send_nElemTriangle*3,        MPI_UNSIGNED_LONG, iDomain, 3,  MPI_COMM_WORLD, &send_req[3]);
        MPI_Isend(Buffer_Send_Rectangle,           Buffer_Send_nElemRectangle*4,       MPI_UNSIGNED_LONG, iDomain, 4,  MPI_COMM_WORLD, &send_req[4]);
        MPI_Isend(Buffer_Send_Tetrahedron,         Buffer_Send_nElemTetrahedron*4,     MPI_UNSIGNED_LONG, iDomain, 5,  MPI_COMM_WORLD, &send_req[5]);
        MPI_Isend(Buffer_Send_Hexahedron,          Buffer_Send_nElemHexahedron*8,      MPI_UNSIGNED_LONG, iDomain, 6,  MPI_COMM_WORLD, &send_req[6]);
        MPI_Isend(Buffer_Send_Wedge,               Buffer_Send_nElemWedge*6,           MPI_UNSIGNED_LONG, iDomain, 7,  MPI_COMM_WORLD, &send_req[7]);
        MPI_Isend(Buffer_Send_Pyramid,             Buffer_Send_nElemPyramid*5,         MPI_UNSIGNED_LONG, iDomain, 8,  MPI_COMM_WORLD, &send_req[8]);
        MPI_Isend(Buffer_Send_BoundLine,           Buffer_Send_nBoundLineTotal*2,      MPI_UNSIGNED_LONG, iDomain, 9,  MPI_COMM_WORLD, &send_req[9]);
        MPI_Isend(Buffer_Send_BoundTriangle,       Buffer_Send_nBoundTriangleTotal*3,  MPI_UNSIGNED_LONG, iDomain, 10, MPI_COMM_WORLD, &send_req[10]);
        MPI_Isend(Buffer_Send_BoundRectangle,      Buffer_Send_nBoundRectangleTotal*4, MPI_UNSIGNED_LONG, iDomain, 11, MPI_COMM_WORLD, &send_req[11]);
        MPI_Isend(Buffer_Send_Local2Global_Marker, Buffer_Send_nMarkerDomain,          MPI_UNSIGNED_LONG, iDomain, 12, MPI_COMM_WORLD, &send_req[12]);
        
        MPI_Isend(Buffer_Send_SendDomain_Periodic,          Buffer_Send_nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, iDomain, 13, MPI_COMM_WORLD, &send_req[13]);
        MPI_Isend(Buffer_Send_SendDomain_PeriodicTrans,     Buffer_Send_nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, iDomain, 14, MPI_COMM_WORLD, &send_req[14]);
        MPI_Isend(Buffer_Send_SendDomain_PeriodicReceptor,  Buffer_Send_nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, iDomain, 15, MPI_COMM_WORLD, &send_req[15]);
        MPI_Isend(Buffer_Send_ReceivedDomain_Periodic,      Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 16, MPI_COMM_WORLD, &send_req[16]);
        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicTrans, Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 17, MPI_COMM_WORLD, &send_req[17]);
        MPI_Isend(Buffer_Send_ReceivedDomain_PeriodicDonor, Buffer_Send_nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, iDomain, 18, MPI_COMM_WORLD, &send_req[18]);
        
        /*--- Wait for this set of non-blocking comm. to complete ---*/
        
        MPI_Waitall(19, send_req, send_stat);
        
#endif
        
      } else {
        
        for (iter = 0; iter < Buffer_Send_nPointTotal*Buffer_Send_nDim; iter++)
          Buffer_Receive_Coord[iter] = Buffer_Send_Coord[iter];
        
        for (iter = 0; iter < Buffer_Send_nPointTotal; iter++) {
          Buffer_Receive_GlobalPointIndex[iter] = Buffer_Send_GlobalPointIndex[iter];
          Buffer_Receive_Color[iter] = Buffer_Send_Color[iter];
        }
        
        for (iter = 0; iter < Buffer_Send_nElemTriangle*3; iter++)
          Buffer_Receive_Triangle[iter] =  Buffer_Send_Triangle[iter];
        
        for (iter = 0; iter < Buffer_Send_nElemRectangle*4; iter++)
          Buffer_Receive_Rectangle[iter] =  Buffer_Send_Rectangle[iter];
        
        for (iter = 0; iter < Buffer_Send_nElemTetrahedron*4; iter++)
          Buffer_Receive_Tetrahedron[iter] =  Buffer_Send_Tetrahedron[iter];
        
        for (iter = 0; iter < Buffer_Send_nElemHexahedron*8; iter++)
          Buffer_Receive_Hexahedron[iter] =  Buffer_Send_Hexahedron[iter];
        
        for (iter = 0; iter < Buffer_Send_nElemWedge*6; iter++)
          Buffer_Receive_Wedge[iter] =  Buffer_Send_Wedge[iter];
        
        for (iter = 0; iter < Buffer_Send_nElemPyramid*5; iter++)
          Buffer_Receive_Pyramid[iter] =  Buffer_Send_Pyramid[iter];
        
        for (iter = 0; iter < Buffer_Send_nBoundLineTotal*2; iter++)
          Buffer_Receive_BoundLine[iter] =  Buffer_Send_BoundLine[iter];
        
        for (iter = 0; iter < Buffer_Send_nBoundTriangleTotal*3; iter++)
          Buffer_Receive_BoundTriangle[iter] =  Buffer_Send_BoundTriangle[iter];
        
        for (iter = 0; iter < Buffer_Send_nBoundRectangleTotal*4; iter++)
          Buffer_Receive_BoundRectangle[iter] =  Buffer_Send_BoundRectangle[iter];
        
        for (iter = 0; iter < Buffer_Send_nMarkerDomain; iter++)
          Buffer_Receive_Local2Global_Marker[iter] =  Buffer_Send_Local2Global_Marker[iter];
        
        for (iter = 0; iter < Buffer_Send_nTotalSendDomain_Periodic; iter++) {
          Buffer_Receive_SendDomain_Periodic[iter] = Buffer_Send_SendDomain_Periodic[iter];
          Buffer_Receive_SendDomain_PeriodicTrans[iter] = Buffer_Send_SendDomain_PeriodicTrans[iter];
          Buffer_Receive_SendDomain_PeriodicReceptor[iter] = Buffer_Send_SendDomain_PeriodicReceptor[iter];
        }
        
        for (iter = 0; iter < Buffer_Send_nTotalReceivedDomain_Periodic; iter++) {
          Buffer_Receive_ReceivedDomain_Periodic[iter] = Buffer_Send_ReceivedDomain_Periodic[iter];
          Buffer_Receive_ReceivedDomain_PeriodicTrans[iter] = Buffer_Send_ReceivedDomain_PeriodicTrans[iter];
          Buffer_Receive_ReceivedDomain_PeriodicDonor[iter] = Buffer_Send_ReceivedDomain_PeriodicDonor[iter];
        }
        
      }
      
      delete[] Buffer_Send_Coord;
      delete[] Buffer_Send_GlobalPointIndex;
      delete[] Buffer_Send_Color;
      delete[] Buffer_Send_Triangle;
      delete[] Buffer_Send_Rectangle;
      delete[] Buffer_Send_Tetrahedron;
      delete[] Buffer_Send_Hexahedron;
      delete[] Buffer_Send_Wedge;
      delete[] Buffer_Send_Pyramid;
      delete[] Buffer_Send_BoundLine;
      delete[] Buffer_Send_BoundTriangle;
      delete[] Buffer_Send_BoundRectangle;
      delete[] Buffer_Send_Local2Global_Marker;
      
      delete[] Buffer_Send_SendDomain_Periodic;
      delete[] Buffer_Send_SendDomain_PeriodicTrans;
      delete[] Buffer_Send_SendDomain_PeriodicReceptor;
      delete[] Buffer_Send_ReceivedDomain_Periodic;
      delete[] Buffer_Send_ReceivedDomain_PeriodicTrans;
      delete[] Buffer_Send_ReceivedDomain_PeriodicDonor;
      
    }
    
    if (rank == iDomain) {
      
      if (rank != MASTER_NODE) {
        
        /*--- Receive the buffers with the geometrical information ---*/
        
#ifdef HAVE_MPI
        
        MPI_Irecv(Buffer_Receive_Coord,               nPointTotal*nDim,       MPI_DOUBLE,        MASTER_NODE, 0,  MPI_COMM_WORLD, &recv_req[0]);
        MPI_Irecv(Buffer_Receive_GlobalPointIndex,    nPointTotal,            MPI_UNSIGNED_LONG, MASTER_NODE, 1,  MPI_COMM_WORLD, &recv_req[1]);
        MPI_Irecv(Buffer_Receive_Color,               nPointTotal,            MPI_UNSIGNED_LONG, MASTER_NODE, 2,  MPI_COMM_WORLD, &recv_req[2]);
        MPI_Irecv(Buffer_Receive_Triangle,            nElemTriangle*3,        MPI_UNSIGNED_LONG, MASTER_NODE, 3,  MPI_COMM_WORLD, &recv_req[3]);
        MPI_Irecv(Buffer_Receive_Rectangle,           nElemRectangle*4,       MPI_UNSIGNED_LONG, MASTER_NODE, 4,  MPI_COMM_WORLD, &recv_req[4]);
        MPI_Irecv(Buffer_Receive_Tetrahedron,         nElemTetrahedron*4,     MPI_UNSIGNED_LONG, MASTER_NODE, 5,  MPI_COMM_WORLD, &recv_req[5]);
        MPI_Irecv(Buffer_Receive_Hexahedron,          nElemHexahedron*8,      MPI_UNSIGNED_LONG, MASTER_NODE, 6,  MPI_COMM_WORLD, &recv_req[6]);
        MPI_Irecv(Buffer_Receive_Wedge,               nElemWedge*6,           MPI_UNSIGNED_LONG, MASTER_NODE, 7,  MPI_COMM_WORLD, &recv_req[7]);
        MPI_Irecv(Buffer_Receive_Pyramid,             nElemPyramid*5,         MPI_UNSIGNED_LONG, MASTER_NODE, 8,  MPI_COMM_WORLD, &recv_req[8]);
        MPI_Irecv(Buffer_Receive_BoundLine,           nBoundLineTotal*2,      MPI_UNSIGNED_LONG, MASTER_NODE, 9,  MPI_COMM_WORLD, &recv_req[9]);
        MPI_Irecv(Buffer_Receive_BoundTriangle,       nBoundTriangleTotal*3,  MPI_UNSIGNED_LONG, MASTER_NODE, 10, MPI_COMM_WORLD, &recv_req[10]);
        MPI_Irecv(Buffer_Receive_BoundRectangle,      nBoundRectangleTotal*4, MPI_UNSIGNED_LONG, MASTER_NODE, 11, MPI_COMM_WORLD, &recv_req[11]);
        MPI_Irecv(Buffer_Receive_Local2Global_Marker, nMarkerDomain,          MPI_UNSIGNED_LONG, MASTER_NODE, 12, MPI_COMM_WORLD, &recv_req[12]);
        
        MPI_Irecv(Buffer_Receive_SendDomain_Periodic,          nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, MASTER_NODE, 13, MPI_COMM_WORLD, &recv_req[13]);
        MPI_Irecv(Buffer_Receive_SendDomain_PeriodicTrans,     nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, MASTER_NODE, 14, MPI_COMM_WORLD, &recv_req[14]);
        MPI_Irecv(Buffer_Receive_SendDomain_PeriodicReceptor,  nTotalSendDomain_Periodic,     MPI_UNSIGNED_LONG, MASTER_NODE, 15, MPI_COMM_WORLD, &recv_req[15]);
        MPI_Irecv(Buffer_Receive_ReceivedDomain_Periodic,      nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 16, MPI_COMM_WORLD, &recv_req[16]);
        MPI_Irecv(Buffer_Receive_ReceivedDomain_PeriodicTrans, nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 17, MPI_COMM_WORLD, &recv_req[17]);
        MPI_Irecv(Buffer_Receive_ReceivedDomain_PeriodicDonor, nTotalReceivedDomain_Periodic, MPI_UNSIGNED_LONG, MASTER_NODE, 18, MPI_COMM_WORLD, &recv_req[18]);
        
        /*--- Wait for this set of non-blocking recv's to complete ---*/
        
        MPI_Waitall(19, recv_req, recv_stat);
        
#endif
        
      }
      
      /*--- Create the domain structures for the points ---*/
      
      nPoint = nPointTotal; iPoint = 0;
      nPointDomain = nPointDomainTotal;
      FinestMGLevel = true;
      node = new CPoint*[nPoint];
      Local_to_Global_Point =  new unsigned long[nPoint];
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        Local_to_Global_Point[iPoint] = Buffer_Receive_GlobalPointIndex[iPoint];
        if ( nDim == 2 ) node[iPoint] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0], Buffer_Receive_Coord[iPoint*nDim+1], Local_to_Global_Point[iPoint], config);
        if ( nDim == 3 ) node[iPoint] = new CPoint(Buffer_Receive_Coord[iPoint*nDim+0], Buffer_Receive_Coord[iPoint*nDim+1], Buffer_Receive_Coord[iPoint*nDim+2], Local_to_Global_Point[iPoint], config);
        node[iPoint]->SetColor(Buffer_Receive_Color[iPoint]);
      }
      
      delete[] Buffer_Receive_Coord;
      delete[] Buffer_Receive_GlobalPointIndex;
      delete[] Buffer_Receive_Color;
      
      /*--- Create the domain structures for the elements ---*/
      
      nElem = nElemTotal; iElem = 0;
      elem = new CPrimalGrid*[nElem];
      
      for (iElemTriangle = 0; iElemTriangle < nElemTriangle; iElemTriangle++) {
        elem[iElem] = new CTriangle(Buffer_Receive_Triangle[iElemTriangle*3+0], Buffer_Receive_Triangle[iElemTriangle*3+1], Buffer_Receive_Triangle[iElemTriangle*3+2], 2);
        iElem++;
      }
      for (iElemRectangle = 0; iElemRectangle < nElemRectangle; iElemRectangle++) {
        elem[iElem] = new CRectangle(Buffer_Receive_Rectangle[iElemRectangle*4+0], Buffer_Receive_Rectangle[iElemRectangle*4+1], Buffer_Receive_Rectangle[iElemRectangle*4+2], Buffer_Receive_Rectangle[iElemRectangle*4+3], 2);
        iElem++;
      }
      for (iElemTetrahedron = 0; iElemTetrahedron < nElemTetrahedron; iElemTetrahedron++) {
        elem[iElem] = new CTetrahedron(Buffer_Receive_Tetrahedron[iElemTetrahedron*4+0], Buffer_Receive_Tetrahedron[iElemTetrahedron*4+1], Buffer_Receive_Tetrahedron[iElemTetrahedron*4+2], Buffer_Receive_Tetrahedron[iElemTetrahedron*4+3]);
        iElem++;
      }
      for (iElemHexahedron = 0; iElemHexahedron < nElemHexahedron; iElemHexahedron++) {
        elem[iElem] = new CHexahedron(Buffer_Receive_Hexahedron[iElemHexahedron*8+0], Buffer_Receive_Hexahedron[iElemHexahedron*8+1], Buffer_Receive_Hexahedron[iElemHexahedron*8+2], Buffer_Receive_Hexahedron[iElemHexahedron*8+3], Buffer_Receive_Hexahedron[iElemHexahedron*8+4], Buffer_Receive_Hexahedron[iElemHexahedron*8+5], Buffer_Receive_Hexahedron[iElemHexahedron*8+6], Buffer_Receive_Hexahedron[iElemHexahedron*8+7]);
        iElem++;
      }
      for (iElemWedge = 0; iElemWedge < nElemWedge; iElemWedge++) {
        elem[iElem] = new CWedge(Buffer_Receive_Wedge[iElemWedge*6+0], Buffer_Receive_Wedge[iElemWedge*6+1], Buffer_Receive_Wedge[iElemWedge*6+2], Buffer_Receive_Wedge[iElemWedge*6+3], Buffer_Receive_Wedge[iElemWedge*6+4], Buffer_Receive_Wedge[iElemWedge*6+5]);
        iElem++;
      }
      for (iElemPyramid = 0; iElemPyramid < nElemPyramid; iElemPyramid++) {
        elem[iElem] = new CPyramid(Buffer_Receive_Pyramid[iElemPyramid*5+0], Buffer_Receive_Pyramid[iElemPyramid*5+1], Buffer_Receive_Pyramid[iElemPyramid*5+2], Buffer_Receive_Pyramid[iElemPyramid*5+3], Buffer_Receive_Pyramid[iElemPyramid*5+4]);
        iElem++;
      }
      
      delete[] Buffer_Receive_Triangle;
      delete[] Buffer_Receive_Rectangle;
      delete[] Buffer_Receive_Tetrahedron;
      delete[] Buffer_Receive_Hexahedron;
      delete[] Buffer_Receive_Wedge;
      delete[] Buffer_Receive_Pyramid;
      
      /*--- Create the domain structures for the boundaries ---*/
      
      nMarker = nMarkerDomain;
      
      nElem_Bound = new unsigned long [MAX_NUMBER_MARKER];
      Local_to_Global_Marker = new unsigned short [MAX_NUMBER_MARKER];
      Tag_to_Marker = new string [MAX_NUMBER_MARKER];
      string *TagBound_Copy = new string [MAX_NUMBER_MARKER];
      short *SendRecv_Copy = new short [MAX_NUMBER_MARKER];
      
      for (iMarker = 0; iMarker < nMarker; iMarker++) nElem_Bound[iMarker] = nVertexDomain[iMarker];
      
      bound = new CPrimalGrid**[nMarker+(overhead*nDomain)];
      for (iMarker = 0; iMarker < nMarker; iMarker++) bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
      
      iBoundLineTotal = 0; iBoundTriangleTotal = 0; iBoundRectangleTotal = 0;
      
      
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        
        iVertexDomain = 0;
        
        for (iBoundLine = 0; iBoundLine < nBoundLine[iMarker]; iBoundLine++) {
          bound[iMarker][iVertexDomain] = new CLine(Buffer_Receive_BoundLine[iBoundLineTotal*2+0],
                                                    Buffer_Receive_BoundLine[iBoundLineTotal*2+1], 2);
          iVertexDomain++; iBoundLineTotal++;
        }
        for (iBoundTriangle = 0; iBoundTriangle < nBoundTriangle[iMarker]; iBoundTriangle++) {
          bound[iMarker][iVertexDomain] = new CTriangle(Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+0],
                                                        Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+1],
                                                        Buffer_Receive_BoundTriangle[iBoundTriangleTotal*3+2], 3);
          iVertexDomain++; iBoundTriangleTotal++;
        }
        for (iBoundRectangle = 0; iBoundRectangle < nBoundRectangle[iMarker]; iBoundRectangle++) {
          bound[iMarker][iVertexDomain] = new CRectangle(Buffer_Receive_BoundRectangle[iBoundRectangleTotal*4+0],
                                                         Buffer_Receive_BoundRectangle[iBoundRectangleTotal*4+1],
                                                         Buffer_Receive_BoundRectangle[iBoundRectangleTotal*4+2],
                                                         Buffer_Receive_BoundRectangle[iBoundRectangleTotal*4+3], 3);
          iVertexDomain++; iBoundRectangleTotal++;
        }
        
        Local_to_Global_Marker[iMarker] = Buffer_Receive_Local2Global_Marker[iMarker];
        
        /*--- Now each domain has the right information ---*/
        
        string Grid_Marker = config->GetMarker_All_TagBound(Local_to_Global_Marker[iMarker]);
        short SendRecv = config->GetMarker_All_SendRecv(Local_to_Global_Marker[iMarker]);
        TagBound_Copy[iMarker] = Grid_Marker;
        SendRecv_Copy[iMarker] = SendRecv;
        
      }
      
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        config->SetMarker_All_TagBound(iMarker, TagBound_Copy[iMarker]);
        config->SetMarker_All_SendRecv(iMarker, SendRecv_Copy[iMarker]);
      }
      
      /*--- Add the new periodic markers to the domain ---*/
      
      iTotalSendDomain_Periodic = 0;
      iTotalReceivedDomain_Periodic = 0;
      
      for (jDomain = 0; jDomain < nDomain; jDomain++) {
        
        if (nSendDomain_Periodic[jDomain] != 0) {
          nVertexDomain[nMarker] = 0;
          bound[nMarker] = new CPrimalGrid* [nSendDomain_Periodic[jDomain]];
          
          iVertex = 0;
          for (iTotalSendDomain_Periodic = 0; iTotalSendDomain_Periodic < nTotalSendDomain_Periodic; iTotalSendDomain_Periodic++) {
            if (Buffer_Receive_SendDomain_PeriodicReceptor[iTotalSendDomain_Periodic] == jDomain) {
              bound[nMarker][iVertex] = new CVertexMPI(Buffer_Receive_SendDomain_Periodic[iTotalSendDomain_Periodic], nDim);
              bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_SendDomain_PeriodicTrans[iTotalSendDomain_Periodic]);
              nVertexDomain[nMarker]++; iVertex++;
            }
          }
          
          Marker_All_SendRecv[nMarker] = jDomain+1;
          nElem_Bound[nMarker] = nVertexDomain[nMarker];
          nMarker++;
        }
        
        if (nReceivedDomain_Periodic[jDomain] != 0) {
          nVertexDomain[nMarker] = 0;
          bound[nMarker] = new CPrimalGrid* [nReceivedDomain_Periodic[jDomain]];
          
          iVertex = 0;
          for (iTotalReceivedDomain_Periodic = 0; iTotalReceivedDomain_Periodic < nTotalReceivedDomain_Periodic; iTotalReceivedDomain_Periodic++) {
            if (Buffer_Receive_ReceivedDomain_PeriodicDonor[iTotalReceivedDomain_Periodic] == jDomain) {
              bound[nMarker][iVertex] = new CVertexMPI(Buffer_Receive_ReceivedDomain_Periodic[iTotalReceivedDomain_Periodic], nDim);
              bound[nMarker][iVertex]->SetRotation_Type(Buffer_Receive_ReceivedDomain_PeriodicTrans[iTotalReceivedDomain_Periodic]);
              nVertexDomain[nMarker]++; iVertex++;
            }
          }
          
          Marker_All_SendRecv[nMarker] = -(jDomain+1);
          nElem_Bound[nMarker] = nVertexDomain[nMarker];
          nMarker++;
        }
        
      }
      
      delete[] TagBound_Copy;
      delete[] SendRecv_Copy;
      
      delete[] Buffer_Receive_BoundLine;
      delete[] Buffer_Receive_BoundTriangle;
      delete[] Buffer_Receive_BoundRectangle;
      delete[] Buffer_Receive_Local2Global_Marker;
      
      delete[] Buffer_Receive_SendDomain_Periodic;
      delete[] Buffer_Receive_SendDomain_PeriodicTrans;
      delete[] Buffer_Receive_SendDomain_PeriodicReceptor;
      delete[] Buffer_Receive_ReceivedDomain_Periodic;
      delete[] Buffer_Receive_ReceivedDomain_PeriodicTrans;
      delete[] Buffer_Receive_ReceivedDomain_PeriodicDonor;
      
    }
    
  }
  
  
  /*--- Set the value of Marker_All_SendRecv and Marker_All_TagBound in the config structure ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
  }
  
  /*--- Set the value of Global_nPoint and Global_nPointDomain ---*/
  
  unsigned long Local_nPoint = nPoint;
  unsigned long Local_nPointDomain = nPointDomain;
  
#ifdef HAVE_MPI
  
  MPI_Allreduce(&Local_nPoint, &Global_nPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  
  /*--- End of the MPI stuff, each node has the right piece of the grid ---*/
  
  MPI_Barrier(MPI_COMM_WORLD);
  
#else
  
  Global_nPoint = Local_nPoint;
  Global_nPointDomain = Local_nPointDomain;
  
#endif
  
  
  if (rank == MASTER_NODE) {
    
    delete [] MarkerIn;
    delete [] nElem_Color;
    
    delete [] Buffer_Send_Center;
    delete [] Buffer_Send_Rotation;
    delete [] Buffer_Send_Translate;
    
    delete [] Buffer_Send_nSendDomain_Periodic;
    delete [] Buffer_Send_nReceivedDomain_Periodic;
    
    delete [] Marker_All_SendRecv_Copy;
    delete [] Marker_All_TagBound_Copy;
    
    for (iDomain = 0; iDomain < nDomain; iDomain++) {
      delete[] Elem_Color[iDomain];
    }
    delete[] Elem_Color;
    delete[] Global_to_Local_Point;
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
      delete VertexIn[iMarker];
    delete[] VertexIn;
    
  }
  
  delete [] nSendDomain_Periodic;
  delete [] nReceivedDomain_Periodic;
  
  delete [] nVertexDomain;
  delete [] nBoundLine;
  delete [] nBoundTriangle;
  delete [] nBoundRectangle;
  delete [] Buffer_Send_nVertexDomain;
  delete [] Buffer_Send_nBoundLine;
  delete [] Buffer_Send_nBoundTriangle;
  delete [] Buffer_Send_nBoundRectangle;
  delete [] Buffer_Send_Marker_All_SendRecv;
  
}

CPhysicalGeometry::~CPhysicalGeometry(void) {
  
  if (Global_to_Local_Point != NULL) delete[] Global_to_Local_Point;
  if (Local_to_Global_Point != NULL) delete[] Local_to_Global_Point;
  if (Global_to_Local_Marker != NULL) delete[] Global_to_Local_Marker;
  if (Local_to_Global_Marker != NULL) delete[] Local_to_Global_Marker;
  
}


void CPhysicalGeometry::SetSendReceive(CConfig *config) {
  
  unsigned short Counter_Send, Counter_Receive, iMarkerSend, iMarkerReceive;
  unsigned long iVertex, LocalNode;
  
  unsigned long  nVertexDomain[MAX_NUMBER_MARKER], iPoint, jPoint, iElem;
  unsigned short nDomain, iNode, iDomain, jDomain, jNode;
  vector<unsigned long>::iterator it;
  
  vector<vector<unsigned long> > SendTransfLocal;	/*!< \brief Vector to store the type of transformation for this send point. */
  vector<vector<unsigned long> > ReceivedTransfLocal;	/*!< \brief Vector to store the type of transformation for this received point. */
	vector<vector<unsigned long> > SendDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
	vector<vector<unsigned long> > ReceivedDomainLocal; /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
  
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  /*--- MPI initialization ---*/
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  
  nDomain = size;
  
  SendTransfLocal.resize(nDomain);
  ReceivedTransfLocal.resize(nDomain);
  SendDomainLocal.resize(nDomain);
  ReceivedDomainLocal.resize(nDomain);
  
  /*--- Loop over the all the points of the element
   to find the points with different colours, and create the send/received list ---*/
  for (iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      iDomain = node[iPoint]->GetColor();
      
      if (iDomain == rank) {
        for(jNode = 0; jNode < elem[iElem]->GetnNodes(); jNode++) {
          jPoint = elem[iElem]->GetNode(jNode);
          jDomain = node[jPoint]->GetColor();
          
          /*--- If different color and connected by an edge, then we add them to the list ---*/
          if (iDomain != jDomain) {
            
            /*--- We send from iDomain to jDomain the value of iPoint, we save the
             global value becuase we need to sort the lists ---*/
            SendDomainLocal[jDomain].push_back(Local_to_Global_Point[iPoint]);
            /*--- We send from jDomain to iDomain the value of jPoint, we save the
             global value becuase we need to sort the lists ---*/
            ReceivedDomainLocal[jDomain].push_back(Local_to_Global_Point[jPoint]);
            
          }
        }
      }
    }
  }
  
  /*--- Sort the points that must be sended and delete repeated points, note
   that the sorting should be done with the global point (not the local) ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    it = unique( SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    SendDomainLocal[iDomain].resize( it - SendDomainLocal[iDomain].begin() );
  }
  
  /*--- Sort the points that must be received and delete repeated points, note
   that the sorting should be done with the global point (not the local) ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    it = unique( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    ReceivedDomainLocal[iDomain].resize( it - ReceivedDomainLocal[iDomain].begin() );
  }
  
  /*--- Create Global to Local Point array, note that the array is smaller (Max_GlobalPoint) than the total
   number of points in the simulation  ---*/
  Max_GlobalPoint = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (Local_to_Global_Point[iPoint] > Max_GlobalPoint)
      Max_GlobalPoint = Local_to_Global_Point[iPoint];
  }
  Global_to_Local_Point =  new long[Max_GlobalPoint+1]; // +1 to include the bigger point.
  
  /*--- Initialization of the array with -1 this is important for the FFD ---*/
  for (iPoint = 0; iPoint < Max_GlobalPoint+1; iPoint++)
    Global_to_Local_Point[iPoint] = -1;
  
  /*--- Set the value of some of the points ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    Global_to_Local_Point[Local_to_Global_Point[iPoint]] = iPoint;
  
  /*--- Add the new MPI send receive boundaries, reset the transformation, and save the local value ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = SendDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {
        SendDomainLocal[iDomain][iVertex] = Global_to_Local_Point[SendDomainLocal[iDomain][iVertex]];
        SendTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }
  
  /*--- Add the new MPI receive boundaries, reset the transformation, and save the local value ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = ReceivedDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {
        ReceivedDomainLocal[iDomain][iVertex] = Global_to_Local_Point[ReceivedDomainLocal[iDomain][iVertex]];
        ReceivedTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }
  
  /*--- First compute the Send/Receive boundaries ---*/
  Counter_Send = 0; 	Counter_Receive = 0;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (SendDomainLocal[iDomain].size() != 0) Counter_Send++;
  
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (ReceivedDomainLocal[iDomain].size() != 0) Counter_Receive++;
  
  iMarkerSend = nMarker - Counter_Send - Counter_Receive;
  iMarkerReceive = nMarker - Counter_Receive;
  
  /*--- First we do the send ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerSend); iVertex++) {
        LocalNode = SendDomainLocal[iDomain][iVertex];
        bound[iMarkerSend][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerSend][iVertex]->SetRotation_Type(SendTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerSend] = iDomain+1;
      iMarkerSend++;
    }
  }
  
  /*--- Second we do the receive ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerReceive); iVertex++) {
        LocalNode = ReceivedDomainLocal[iDomain][iVertex];
        bound[iMarkerReceive][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerReceive][iVertex]->SetRotation_Type(ReceivedTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerReceive] = -(iDomain+1);
      iMarkerReceive++;
    }
  }
  
}


void CPhysicalGeometry::SetBoundaries(CConfig *config) {
  
  unsigned long iElem_Bound, TotalElem, *nElem_Bound_Copy, iVertex_;
  string Grid_Marker;
  unsigned short iDomain, nDomain, iMarkersDomain, iLoop, *DomainCount, nMarker_Physical, Duplicate_SendReceive, *DomainSendCount, **DomainSendMarkers, *DomainReceiveCount, **DomainReceiveMarkers, nMarker_SendRecv, iMarker, iMarker_;
  CPrimalGrid*** bound_Copy;
  short *Marker_All_SendRecv_Copy;
  bool CheckStart;
  
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  /*--- MPI initialization ---*/
  MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif
  
  nDomain = size+1;
  
  /*--- Count the number of physical markers
   in the boundaries ---*/
  
  nMarker_Physical = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      nMarker_Physical++;
    }
  }
  
  /*--- Identify if there are markers that send/received with the same domain,
   they should be together---*/
  
  Duplicate_SendReceive = 0;
  for (iLoop = 0; iLoop < 2; iLoop++) {
    
    DomainCount = new unsigned short [nDomain];
    
    for (iDomain = 0; iDomain < nDomain; iDomain++)
      DomainCount[iDomain] = 0;
    
    if (iLoop == 0) {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == iDomain) DomainCount[iDomain]++;
    }
    else {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == -iDomain) DomainCount[iDomain]++;
    }
    
    for (iDomain = 0; iDomain < nDomain; iDomain++)
      if (DomainCount[iDomain] > 1) Duplicate_SendReceive++;
    
    delete [] DomainCount;
    
  }
  
  DomainSendCount = new unsigned short [nDomain];
  DomainSendMarkers = new unsigned short *[nDomain];
  DomainReceiveCount = new unsigned short [nDomain];
  DomainReceiveMarkers = new unsigned short *[nDomain];
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    DomainSendCount[iDomain] = 0;
    DomainSendMarkers[iDomain] = new unsigned short [nMarker];
    
    DomainReceiveCount[iDomain] = 0;
    DomainReceiveMarkers[iDomain] = new unsigned short [nMarker];
  }
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
        if (Marker_All_SendRecv[iMarker] == iDomain) {
          DomainSendMarkers[iDomain][DomainSendCount[iDomain]] = iMarker;
          DomainSendCount[iDomain]++;
        }
        if (Marker_All_SendRecv[iMarker] == -iDomain) {
          DomainReceiveMarkers[iDomain][DomainReceiveCount[iDomain]] = iMarker;
          DomainReceiveCount[iDomain]++;
        }
      }
    }
  }
  
  /*--- Create an structure to store the Send/Receive
   boundaries, because they require some reorganization ---*/
  
  nMarker_SendRecv = nMarker - nMarker_Physical - Duplicate_SendReceive;
  bound_Copy = new CPrimalGrid**[nMarker_Physical + nMarker_SendRecv];
  nElem_Bound_Copy = new unsigned long [nMarker_Physical + nMarker_SendRecv];
  Marker_All_SendRecv_Copy = new short [nMarker_Physical + nMarker_SendRecv];
  iMarker_ = nMarker_Physical;
  iVertex_ = 0;
  CheckStart = false;
  
  /*--- Copy and allocate the physical markers in the data structure ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      nElem_Bound_Copy[iMarker] = nElem_Bound[iMarker];
      bound_Copy[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == LINE)
          bound_Copy[iMarker][iElem_Bound] = new CLine(bound[iMarker][iElem_Bound]->GetNode(0),
                                                       bound[iMarker][iElem_Bound]->GetNode(1), 2);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE)
          
          bound_Copy[iMarker][iElem_Bound] = new CTriangle(bound[iMarker][iElem_Bound]->GetNode(0),
                                                           bound[iMarker][iElem_Bound]->GetNode(1),
                                                           bound[iMarker][iElem_Bound]->GetNode(2), 3);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == RECTANGLE)
          bound_Copy[iMarker][iElem_Bound] = new CRectangle(bound[iMarker][iElem_Bound]->GetNode(0),
                                                            bound[iMarker][iElem_Bound]->GetNode(1),
                                                            bound[iMarker][iElem_Bound]->GetNode(2),
                                                            bound[iMarker][iElem_Bound]->GetNode(3), 3);
      }
    }
  }
  
  
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    
    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/
    
    if (DomainSendCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];
    }
    
    for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }
      
    }
    
    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/
    
    if (DomainReceiveCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];
      
    }
    
    for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }
      
    }
    
  }
  
  delete [] DomainSendCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete DomainSendMarkers[iDomain];
  delete[] DomainSendMarkers;
  
  delete [] DomainReceiveCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete DomainReceiveMarkers[iDomain];
  delete[] DomainReceiveMarkers;
  
  /*--- Deallocate the bound variables ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    delete bound[iMarker];
  delete bound;
  
  /*--- Allocate the new bound variables, and set the number of markers ---*/
  
  bound = bound_Copy;
  nMarker = nMarker_Physical + nMarker_SendRecv;
  config->SetnMarker_All(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElem_Bound[iMarker] = nElem_Bound_Copy[iMarker];
  }
  for (iMarker = nMarker_Physical; iMarker < nMarker; iMarker++) {
    Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
    config->SetMarker_All_TagBound(iMarker, "SEND_RECEIVE");
  }
  
  /*--- Update config information storing the boundary information in the right place ---*/
  
  for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
    
    string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
    
    if (Marker_Tag != "SEND_RECEIVE") {
      
      /*--- Update config information storing the boundary information in the right place ---*/
      
      Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
      config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
      config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
      config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
      config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
      config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
      config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
      config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
      config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
      config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));

    }
    
    /*--- Send-Receive boundaries definition ---*/
    
    else {
      config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (config->GetMarker_All_SendRecv(iMarker) < 0)
          node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
      }
      
    }
    
    /*--- Loop over the surface element to set the boundaries ---*/
    
    unsigned long Point_Surface, iElem_Surface;
    unsigned short iNode_Surface;
    
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL &&
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX &&
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[Point_Surface]->SetSolidBoundary(true);
      }
    }
    
  }
  
}

void CPhysicalGeometry::Read_SU2_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
  
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short VTK_Type, iMarker, iChar, iCount = 0;
  unsigned long iElem_Bound = 0, iPoint = 0, ielem_div = 0, ielem = 0, *Local2Global = NULL, vnodes_edge[2], vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8],
  vnodes_wedge[6], vnodes_pyramid[5], dummyLong, GlobalIndex, iElem;
  char cstr[MAX_STRING_SIZE];
  double Coord_2D[2], Coord_3D[3], dummyDouble;
  string::size_type position;
  bool time_spectral = (config->GetUnsteady_Simulation() == TIME_SPECTRAL);
  int rank = MASTER_NODE, size = SINGLE_NODE;
  bool domain_flag = false;
  bool found_transform = false;
  nZone = val_nZone;
  
  /*--- Initialize counters for local/global points & elements ---*/
#ifdef HAVE_MPI
  unsigned long LocalIndex;
  unsigned long Local_nPoint, Local_nPointDomain;
  unsigned long Local_nElem;
  unsigned long Local_nElemTri, Local_nElemQuad, Local_nElemTet;
  unsigned long Local_nElemHex, Local_nElemWedge, Local_nElemPyramid;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  FinestMGLevel = true;
  Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_wedge    = 0; Global_nelem_wedge    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Open grid file ---*/
  
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  
  /*--- Check the grid ---*/
  
  if (mesh_file.fail()) {
    cout << "There is no geometry file (CPhysicalGeometry)!! " << cstr << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- If divided grid, we need the global index to
   perform the right element division, this is just a hack in the future we
   should read first the coordinates ---*/
  
  if (config->GetDivide_Element()) {
    
    unsigned long nDim_ = 0, nElem_ = 0, nPoint_ = 0;
    
    while (getline (mesh_file, text_line)) {
      
      position = text_line.find ("NDIME=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nDim_ = atoi(text_line.c_str());
      }
      
      position = text_line.find ("NELEM=",0);
      if (position != string::npos) {
        text_line.erase (0,6); nElem_ = atoi(text_line.c_str());
        for (iElem = 0; iElem < nElem_; iElem ++)
          getline(mesh_file, text_line);
      }
      
      position = text_line.find ("NPOIN=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        
        /*--- Now read the number of points and ghost points. ---*/
        
        stringstream  stream_line(text_line);
        stream_line >> nPoint_;
        
        Local2Global = new unsigned long [nPoint_];
        
        for (iPoint = 0; iPoint < nPoint_; iPoint ++) {
          getline(mesh_file, text_line);
          istringstream point_line(text_line);
          if (size == 1) { Local2Global[iPoint] = iPoint; }
          else {
            point_line >> dummyDouble; point_line >> dummyDouble; if (nDim_ == 3) point_line >> dummyDouble;
            point_line >> dummyLong; point_line >> Local2Global[iPoint];
          }
        }
        
      }
    }
    
    /*--- Close and open again the grid file ---*/
    
    mesh_file.close();
    strcpy (cstr, val_mesh_filename.c_str());
    mesh_file.open(cstr, ios::in);
    
  }
  
  /*--- If more than one, find the zone in the mesh file ---*/
  
  if (val_nZone > 1 || time_spectral) {
    if (time_spectral) {
      if (rank == MASTER_NODE) cout << "Reading time spectral instance " << val_iZone+1 << ":" << endl;
    } else {
      while (getline (mesh_file,text_line)) {
        /*--- Search for the current domain ---*/
        position = text_line.find ("IZONE=",0);
        if (position != string::npos) {
          text_line.erase (0,6);
          unsigned short jDomain = atoi(text_line.c_str());
          if (jDomain == val_iZone+1) {
            if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone+1 << ":" << endl;
            break;
          }
        }
      }
    }
  }
  
  /*--- Read grid file with format SU2 ---*/
  
  while (getline (mesh_file,text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      if (domain_flag == false) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
        if (rank == MASTER_NODE) {
          if (nDim == 2) cout << "Two dimensional problem." << endl;
          if (nDim == 3) cout << "Three dimensional problem." << endl;
        }
        domain_flag = true;
      } else { break; }
    }
    
    /*--- Read the information about inner elements ---*/
    
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem = atoi(text_line.c_str());
      
#ifdef HAVE_MPI
      if (config->GetKind_SU2() != SU2_PRT) {
        Local_nElem = nElem;
        MPI_Allreduce(&Local_nElem, &Global_nElem, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
        Local_nElem = nElem;
        Global_nElem = Local_nElem;
      }
#else
      Global_nElem = nElem;
#endif
      
      if ((rank == MASTER_NODE) && (size > SINGLE_NODE) && (config->GetKind_SU2() != SU2_PRT))
        cout << Global_nElem << " interior elements including halo cells. " << endl;
      else if (rank == MASTER_NODE)
        cout << Global_nElem << " interior elements. " << endl;
      
      /*--- Allocate space for elements ---*/
      
      if (!config->GetDivide_Element()) elem = new CPrimalGrid*[nElem];
      else {
        if (nDim == 2) elem = new CPrimalGrid*[2*nElem];
        if (nDim == 3) elem = new CPrimalGrid*[6*nElem];
      }
      
      unsigned short IndirectionPrism[6][6], CT_FromVTK_Prism[6], IndirectionHexa[9][9];
      unsigned short temp, zero, one, two, three, four, five, six, seven, eight, iNode, smallestNode = 0, lookupindex;
      unsigned long smallest;
      
      /*--- Indirection matrix for dividing prisms into tets, using vtk format, and conversion table for local numbering of prisms in vtk format ---*/
      
      CT_FromVTK_Prism[0] = 1;  CT_FromVTK_Prism[1] = 3;  CT_FromVTK_Prism[2] = 2;  CT_FromVTK_Prism[3] = 4;  CT_FromVTK_Prism[4] = 6;  CT_FromVTK_Prism[5] = 5;
      
      IndirectionPrism[0][0] = 0; IndirectionPrism[0][1] = 2; IndirectionPrism[0][2] = 1; IndirectionPrism[0][3] = 3; IndirectionPrism[0][4] = 5; IndirectionPrism[0][5] = 4;
      IndirectionPrism[1][0] = 2; IndirectionPrism[1][1] = 1; IndirectionPrism[1][2] = 0; IndirectionPrism[1][3] = 5; IndirectionPrism[1][4] = 4; IndirectionPrism[1][5] = 3;
      IndirectionPrism[2][0] = 1; IndirectionPrism[2][1] = 0; IndirectionPrism[2][2] = 2; IndirectionPrism[2][3] = 4; IndirectionPrism[2][4] = 3; IndirectionPrism[2][5] = 5;
      IndirectionPrism[3][0] = 3; IndirectionPrism[3][1] = 4; IndirectionPrism[3][2] = 5; IndirectionPrism[3][3] = 0; IndirectionPrism[3][4] = 1; IndirectionPrism[3][5] = 2;
      IndirectionPrism[4][0] = 5; IndirectionPrism[4][1] = 3; IndirectionPrism[4][2] = 4; IndirectionPrism[4][3] = 2; IndirectionPrism[4][4] = 0; IndirectionPrism[4][5] = 1;
      IndirectionPrism[5][0] = 4; IndirectionPrism[5][1] = 5; IndirectionPrism[5][2] = 3; IndirectionPrism[5][3] = 1; IndirectionPrism[5][4] = 2; IndirectionPrism[5][5] = 0;
      
      /*--- Indirection matrix for dividing hexahedron into tets ---*/
      
      IndirectionHexa[1][1] = 1; IndirectionHexa[1][2] = 2; IndirectionHexa[1][3] = 3; IndirectionHexa[1][4] = 4; IndirectionHexa[1][5] = 5; IndirectionHexa[1][6] = 6; IndirectionHexa[1][7] = 7; IndirectionHexa[1][8] = 8;
      IndirectionHexa[2][1] = 2; IndirectionHexa[2][2] = 1; IndirectionHexa[2][3] = 5; IndirectionHexa[2][4] = 6; IndirectionHexa[2][5] = 3; IndirectionHexa[2][6] = 4; IndirectionHexa[2][7] = 8; IndirectionHexa[2][8] = 7;
      IndirectionHexa[3][1] = 3; IndirectionHexa[3][2] = 2; IndirectionHexa[3][3] = 6; IndirectionHexa[3][4] = 7; IndirectionHexa[3][5] = 4; IndirectionHexa[3][6] = 1; IndirectionHexa[3][7] = 5; IndirectionHexa[3][8] = 8;
      IndirectionHexa[4][1] = 4; IndirectionHexa[4][2] = 1; IndirectionHexa[4][3] = 2; IndirectionHexa[4][4] = 3; IndirectionHexa[4][5] = 8; IndirectionHexa[4][6] = 5; IndirectionHexa[4][7] = 6; IndirectionHexa[4][8] = 7;
      IndirectionHexa[5][1] = 5; IndirectionHexa[5][2] = 1; IndirectionHexa[5][3] = 4; IndirectionHexa[5][4] = 8; IndirectionHexa[5][5] = 6; IndirectionHexa[5][6] = 2; IndirectionHexa[5][7] = 3; IndirectionHexa[5][8] = 7;
      IndirectionHexa[6][1] = 6; IndirectionHexa[6][2] = 2; IndirectionHexa[6][3] = 1; IndirectionHexa[6][4] = 5; IndirectionHexa[6][5] = 7; IndirectionHexa[6][6] = 3; IndirectionHexa[6][7] = 4; IndirectionHexa[6][8] = 8;
      IndirectionHexa[7][1] = 7; IndirectionHexa[7][2] = 3; IndirectionHexa[7][3] = 2; IndirectionHexa[7][4] = 6; IndirectionHexa[7][5] = 8; IndirectionHexa[7][6] = 4; IndirectionHexa[7][7] = 1; IndirectionHexa[7][8] = 5;
      IndirectionHexa[8][1] = 8; IndirectionHexa[8][2] = 4; IndirectionHexa[8][3] = 3; IndirectionHexa[8][4] = 7; IndirectionHexa[8][5] = 5; IndirectionHexa[8][6] = 1; IndirectionHexa[8][7] = 2; IndirectionHexa[8][8] = 6;
      
      
      /*--- Loop over all the volumetric elements ---*/
      
      while (ielem_div < nElem) {
        getline(mesh_file,text_line);
        istringstream elem_line(text_line);
        
        elem_line >> VTK_Type;
        
        switch(VTK_Type) {
          case TRIANGLE:
            
            elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
            elem[ielem] = new CTriangle(vnodes_triangle[0], vnodes_triangle[1], vnodes_triangle[2], 2);
            ielem_div++; ielem++; nelem_triangle++;
            break;
            
          case RECTANGLE:
            
            elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
            if (!config->GetDivide_Element()) {
              elem[ielem] = new CRectangle(vnodes_quad[0], vnodes_quad[1], vnodes_quad[2], vnodes_quad[3], 2);
              ielem++; nelem_quad++; }
            else {
              
              unsigned long index1, index2;
              bool division1 = false;
              index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
              index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
              if (index1 < index2) division1 = true;
              
              if (division1) {
                elem[ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 2);
                ielem++; nelem_triangle++;
                elem[ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 2);
                ielem++; nelem_triangle++;
              }
              else {
                elem[ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 2);
                ielem++; nelem_triangle++;
                elem[ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 2);
                ielem++; nelem_triangle++;
              }
              
            }
            
            ielem_div++;
            break;
            
          case TETRAHEDRON:
            
            elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
            elem[ielem] = new CTetrahedron(vnodes_tetra[0], vnodes_tetra[1], vnodes_tetra[2], vnodes_tetra[3]);
            ielem_div++; ielem++; nelem_tetra++;
            break;
            
          case HEXAHEDRON:
            
            elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2];
            elem_line >> vnodes_hexa[3]; elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5];
            elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
            
            if (!config->GetDivide_Element()) {
              elem[ielem] = new CHexahedron(vnodes_hexa[0], vnodes_hexa[1], vnodes_hexa[2], vnodes_hexa[3],
                                            vnodes_hexa[4], vnodes_hexa[5], vnodes_hexa[6], vnodes_hexa[7]);
              ielem++; nelem_hexa++;
            }
            else {
              
              smallest = Local2Global[vnodes_hexa[0]]; smallestNode = 0;
              for (iNode = 1; iNode < 8; iNode ++) {
                if ( smallest > Local2Global[vnodes_hexa[iNode]]) {
                  smallest = Local2Global[vnodes_hexa[iNode]];
                  smallestNode = iNode;
                }
              }
              
              one  = IndirectionHexa[smallestNode+1][1] - 1;
              two  = IndirectionHexa[smallestNode+1][2] - 1;
              three  = IndirectionHexa[smallestNode+1][3] - 1;
              four = IndirectionHexa[smallestNode+1][4] - 1;
              five = IndirectionHexa[smallestNode+1][5] - 1;
              six = IndirectionHexa[smallestNode+1][6] - 1;
              seven = IndirectionHexa[smallestNode+1][7] - 1;
              eight = IndirectionHexa[smallestNode+1][8] - 1;
              
              unsigned long index1, index2;
              unsigned short code1 = 0, code2 = 0, code3 = 0;
              
              index1 = min(Local2Global[vnodes_hexa[two]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[six]]);
              if (index1 < index2) code1 = 1;
              
              index1 = min(Local2Global[vnodes_hexa[four]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[three]], Local2Global[vnodes_hexa[eight]]);
              if (index1 < index2) code2 = 1;
              
              index1 = min(Local2Global[vnodes_hexa[five]], Local2Global[vnodes_hexa[seven]]);
              index2 = min(Local2Global[vnodes_hexa[six]], Local2Global[vnodes_hexa[eight]]);
              if (index1 < index2) code3 = 1;
              
              /*--- Rotation of 120 degrees ---*/
              
              if ((!code1 && !code2 && code3) || (code1 && code2 && !code3)) {
                temp = two; two = five; five = four; four = temp;
                temp = six; six = eight; eight = three; three = temp;
              }
              
              /*--- Rotation of 240 degrees ---*/
              
              if ((!code1 && code2 && !code3) || (code1 && !code2 && code3)) {
                temp = two; two = four; four = five; five = temp;
                temp = six; six = three; three = eight; eight = temp;
              }
              
              if ((code1 + code2 + code3) == 0) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[three], vnodes_hexa[eight], vnodes_hexa[six], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
              }
              if ((code1 + code2 + code3) == 1) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[eight], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[eight], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[three], vnodes_hexa[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[two], vnodes_hexa[three]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[eight], vnodes_hexa[seven], vnodes_hexa[three]);
                ielem++; nelem_tetra++;
                
              }
              if ((code1 + code2 + code3) == 2) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[five], vnodes_hexa[seven], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[four], vnodes_hexa[eight]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[two], vnodes_hexa[three], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[seven], vnodes_hexa[five], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[seven], vnodes_hexa[six]);
                ielem++; nelem_tetra++;
              }
              if ((code1 + code2 + code3) == 3) {
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[three], vnodes_hexa[four], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[four], vnodes_hexa[eight], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[eight], vnodes_hexa[five], vnodes_hexa[seven]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[one], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[six], vnodes_hexa[seven], vnodes_hexa[one]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_hexa[two], vnodes_hexa[seven], vnodes_hexa[three], vnodes_hexa[one]);
                ielem++; nelem_tetra++;
              }
              
            }
            ielem_div++;
            break;
            
          case WEDGE:
            
            elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2];
            elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
            
            if (!config->GetDivide_Element()) {
              
              elem[ielem] = new CWedge(vnodes_wedge[0],vnodes_wedge[1],vnodes_wedge[2],vnodes_wedge[3],vnodes_wedge[4],vnodes_wedge[5]);
              ielem++; nelem_wedge++;
              
            }
            else {
              
              smallest = Local2Global[vnodes_wedge[0]]; smallestNode = 0;
              for (iNode = 1; iNode < 6; iNode ++) {
                if ( smallest > Local2Global[vnodes_wedge[iNode]]) {
                  smallest = Local2Global[vnodes_wedge[iNode]];
                  smallestNode = iNode;
                }
              }
              
              lookupindex = (CT_FromVTK_Prism[smallestNode] - 1);
              zero  = IndirectionPrism[lookupindex][0];
              one  = IndirectionPrism[lookupindex][1];
              two  = IndirectionPrism[lookupindex][2];
              three = IndirectionPrism[lookupindex][3];
              four = IndirectionPrism[lookupindex][4];
              five = IndirectionPrism[lookupindex][5];
              
              unsigned long index1, index2;
              bool division = false;
              index1 = min(Local2Global[vnodes_wedge[one]], Local2Global[vnodes_wedge[five]]);
              index2 = min(Local2Global[vnodes_wedge[two]], Local2Global[vnodes_wedge[four]]);
              if (index1 < index2) division = true;
              
              if (division) {
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[five], vnodes_wedge[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                ielem++; nelem_tetra++;
              }
              else {
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[one], vnodes_wedge[two], vnodes_wedge[four]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[two], vnodes_wedge[five]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_wedge[zero], vnodes_wedge[four], vnodes_wedge[five], vnodes_wedge[three]);
                ielem++; nelem_tetra++;
              }
            }
            
            ielem_div++;
            break;
          case PYRAMID:
            
            elem_line >> vnodes_pyramid[0]; elem_line >> vnodes_pyramid[1]; elem_line >> vnodes_pyramid[2];
            elem_line >> vnodes_pyramid[3]; elem_line >> vnodes_pyramid[4];
            
            if (!config->GetDivide_Element()) {
              
              elem[ielem] = new CPyramid(vnodes_pyramid[0],vnodes_pyramid[1],vnodes_pyramid[2],vnodes_pyramid[3],vnodes_pyramid[4]);
              ielem++; nelem_pyramid++;
              
            }
            else {
              
              unsigned long index1, index2;
              bool division = false;
              index1 = min(Local2Global[vnodes_pyramid[0]], Local2Global[vnodes_pyramid[2]]);
              index2 = min(Local2Global[vnodes_pyramid[1]], Local2Global[vnodes_pyramid[3]]);
              if (index1 < index2) division = true;
              
              if (division) {
                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_pyramid[0], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
              }
              else {
                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[2], vnodes_pyramid[3], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
                elem[ielem] = new CTetrahedron(vnodes_pyramid[1], vnodes_pyramid[3], vnodes_pyramid[0], vnodes_pyramid[4]);
                ielem++; nelem_tetra++;
              }
              
            }
            
            ielem_div++;
            break;
        }
      }
      
      if (config->GetDivide_Element()) nElem = nelem_triangle + nelem_quad + nelem_tetra + nelem_hexa + nelem_wedge + nelem_pyramid;
      
      /*--- Communicate the number of each element type to all processors. ---*/
      
#ifdef HAVE_MPI
      if (config->GetKind_SU2() != SU2_PRT) {
        Local_nElemTri = nelem_triangle;
        Local_nElemQuad = nelem_quad;
        Local_nElemTet = nelem_tetra;
        Local_nElemHex = nelem_hexa;
        Local_nElemWedge = nelem_wedge;
        Local_nElemPyramid = nelem_pyramid;
        MPI_Allreduce(&Local_nElemTri, &Global_nelem_triangle, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nElemQuad, &Global_nelem_quad, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nElemTet, &Global_nelem_tetra, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nElemHex, &Global_nelem_hexa, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nElemWedge, &Global_nelem_wedge, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      }
      else {
        Local_nElemTri = nelem_triangle;
        Global_nelem_triangle = Local_nElemTri;
        Local_nElemQuad = nelem_quad;
        Global_nelem_quad = Local_nElemQuad;
        Local_nElemTet = nelem_tetra;
        Global_nelem_tetra = Local_nElemTet;
        Local_nElemHex = nelem_hexa;
        Global_nelem_hexa = Local_nElemHex;
        Local_nElemWedge = nelem_wedge;
        Global_nelem_wedge = Local_nElemWedge;
        Local_nElemPyramid = nelem_pyramid;
        Global_nelem_pyramid = Local_nElemPyramid;
      }
#else
      Global_nelem_triangle = nelem_triangle;
      Global_nelem_quad     = nelem_quad;
      Global_nelem_tetra    = nelem_tetra;
      Global_nelem_hexa     = nelem_hexa;
      Global_nelem_wedge    = nelem_wedge;
      Global_nelem_pyramid  = nelem_pyramid;
#endif
      
      /*--- Print information about the elements to the console ---*/
      
      if (rank == MASTER_NODE) {
        if (Global_nelem_triangle > 0)  cout << Global_nelem_triangle << " triangles." << endl;
        if (Global_nelem_quad > 0)      cout << Global_nelem_quad << " quadrilaterals." << endl;
        if (Global_nelem_tetra > 0)     cout << Global_nelem_tetra << " tetrahedra." << endl;
        if (Global_nelem_hexa > 0)      cout << Global_nelem_hexa << " hexahedra." << endl;
        if (Global_nelem_wedge > 0)     cout << Global_nelem_wedge << " prisms." << endl;
        if (Global_nelem_pyramid > 0)   cout << Global_nelem_pyramid << " pyramids." << endl;
      }
    }
    
    /*--- Read number of points ---*/
    
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      
      /*--- Check for ghost points. ---*/
      stringstream test_line(text_line);
      while (test_line >> dummyLong)
        iCount++;
      
      /*--- Now read and store the number of points and possible ghost points. ---*/
      
      stringstream  stream_line(text_line);
      if (iCount == 2) {
        stream_line >> nPoint;
        stream_line >> nPointDomain;
        
        /*--- Set some important point information for parallel simulations. ---*/
        
#ifdef HAVE_MPI
        if (config->GetKind_SU2() != SU2_PRT) {
          Local_nPoint = nPoint; Local_nPointDomain = nPointDomain;
          MPI_Allreduce(&Local_nPoint, &Global_nPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        }
        else {
          Local_nPoint = nPoint; Local_nPointDomain = nPointDomain;
          Global_nPoint = Local_nPoint;
          Global_nPointDomain = Local_nPointDomain;
        }
#else
        Global_nPoint = nPoint;
        Global_nPointDomain = nPointDomain;
#endif
        
        if (rank == MASTER_NODE)
          cout << Global_nPointDomain << " points, and " << Global_nPoint-Global_nPointDomain << " ghost points." << endl;
        
      }
      else if (iCount == 1) {
        stream_line >> nPoint;
        nPointDomain = nPoint;
        Global_nPointDomain = nPoint;
        Global_nPoint = nPoint;
        if (rank == MASTER_NODE) cout << nPoint << " points." << endl;
      }
      else {
        cout << "NPOIN improperly specified!!" << endl;
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }
      
      node = new CPoint*[nPoint];
      iPoint = 0;
      while (iPoint < nPoint) {
        getline(mesh_file,text_line);
        istringstream point_line(text_line);
        switch(nDim) {
          case 2:
            GlobalIndex = iPoint;
#ifndef HAVE_MPI
            point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
#else
            if (size > SINGLE_NODE) { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; point_line >> LocalIndex; point_line >> GlobalIndex; }
            else { point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; LocalIndex = iPoint; GlobalIndex = iPoint; }
#endif
            node[iPoint] = new CPoint(Coord_2D[0], Coord_2D[1], GlobalIndex, config);
            iPoint++; break;
          case 3:
            GlobalIndex = iPoint;
#ifndef HAVE_MPI
            point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
#else
            if (size > SINGLE_NODE) { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; point_line >> LocalIndex; point_line >> GlobalIndex; }
            else { point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; LocalIndex = iPoint; GlobalIndex = iPoint; }
#endif
            node[iPoint] = new CPoint(Coord_3D[0], Coord_3D[1], Coord_3D[2], GlobalIndex, config);
            iPoint++; break;
        }
      }
    }
    
    
    /*--- Read number of markers ---*/
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nMarker = atoi(text_line.c_str());
      if (size == 1) cout << nMarker << " surface markers." << endl;
      config->SetnMarker_All(nMarker);
      bound = new CPrimalGrid**[nMarker];
      nElem_Bound = new unsigned long [nMarker];
      Tag_to_Marker = new string [MAX_NUMBER_MARKER];
      
      for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
        getline (mesh_file,text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if(position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Physical boundaries definition ---*/
        if (Marker_Tag != "SEND_RECEIVE") {
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          if (size == 1)
            cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
          
          
          /*--- Allocate space for elements ---*/
          if (!config->GetDivide_Element()) bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          else {
            if (nDim == 2) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
            if (nDim == 3) bound[iMarker] = new CPrimalGrid* [2*nElem_Bound[iMarker]];;
          }
          
          nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type;
            switch(VTK_Type) {
              case LINE:
                
                if (nDim == 3) {
                  cout << "Please remove line boundary conditions from the mesh file!" << endl;
#ifndef HAVE_MPI
                  exit(EXIT_FAILURE);
#else
                  MPI_Abort(MPI_COMM_WORLD,1);
                  MPI_Finalize();
#endif
                }
                
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                bound[iMarker][ielem] = new CLine(vnodes_edge[0],vnodes_edge[1],2);
                ielem++; nelem_edge_bound++; break;
                
              case TRIANGLE:
                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                ielem++; nelem_triangle_bound++; break;
                
              case RECTANGLE:
                
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                
                if (!config->GetDivide_Element()) {
                  
                  bound[iMarker][ielem] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                  ielem++; nelem_quad_bound++;
                  
                }
                else {
                  
                  unsigned long index1, index2;
                  bool division1 = false;
                  index1 = min(Local2Global[vnodes_quad[0]], Local2Global[vnodes_quad[2]]);
                  index2 = min(Local2Global[vnodes_quad[1]], Local2Global[vnodes_quad[3]]);
                  if (index1 < index2) division1 = true;
                  
                  if (division1) {
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[0], vnodes_quad[2], vnodes_quad[1], 3);
                    ielem++; nelem_triangle_bound++;
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[3], vnodes_quad[2], vnodes_quad[0], 3);
                    ielem++; nelem_triangle_bound++;
                  }
                  else {
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[2], vnodes_quad[1], vnodes_quad[3], 3);
                    ielem++; nelem_triangle_bound++;
                    bound[iMarker][ielem] = new CTriangle(vnodes_quad[1], vnodes_quad[0], vnodes_quad[3], 3);
                    ielem++; nelem_triangle_bound++;
                  }
                  
                }
                
                break;
                
                
            }
          }
          if (config->GetDivide_Element()) nElem_Bound[iMarker] = nelem_edge_bound + nelem_triangle_bound + nelem_quad_bound;
          
          /*--- Update config information storing the boundary information in the right place ---*/
          Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
          config->SetMarker_All_TagBound(iMarker, Marker_Tag);
          config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
          config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
          config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
          config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
          config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
          config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
          config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
          config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
          config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
          config->SetMarker_All_SendRecv(iMarker, NONE);
          
        }
        
        /*--- Send-Receive boundaries definition ---*/
        else {
          unsigned long nelem_vertex = 0, vnodes_vertex;
          unsigned short transform;
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          
          nelem_vertex = 0; ielem = 0;
          getline (mesh_file,text_line); text_line.erase (0,8);
          config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
          config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
            
            bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
            bound[iMarker][ielem]->SetRotation_Type(transform);
            ielem++; nelem_vertex++;
            if (config->GetMarker_All_SendRecv(iMarker) < 0)
              node[vnodes_vertex]->SetDomain(false);
          }
          
        }
        
      }
    }
    /*--- Read periodic transformation info (center, rotation, translation) ---*/
    position = text_line.find ("NPERIODIC=",0);
    if (position != string::npos) {
      unsigned short nPeriodic, iPeriodic, iIndex;
      
      /*--- Set bool signifying that periodic transormations were found ---*/
      found_transform = true;
      
      /*--- Read and store the number of transformations. ---*/
      text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
      if (rank == MASTER_NODE) {
        if (nPeriodic - 1 != 0)
          cout << nPeriodic - 1 << " periodic transformations." << endl;
      }
      config->SetnPeriodicIndex(nPeriodic);
      
      /*--- Store center, rotation, & translation in that order for each. ---*/
      for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
        getline (mesh_file,text_line);
        position = text_line.find ("PERIODIC_INDEX=",0);
        if (position != string::npos) {
          text_line.erase (0,15); iIndex = atoi(text_line.c_str());
          if (iIndex != iPeriodic) {
            cout << "PERIODIC_INDEX out of order in SU2 file!!" << endl;
#ifndef HAVE_MPI
            exit(EXIT_FAILURE);
#else
            MPI_Abort(MPI_COMM_WORLD,1);
            MPI_Finalize();
#endif
          }
        }
        double* center    = new double[3];
        double* rotation  = new double[3];
        double* translate = new double[3];
        getline (mesh_file,text_line);
        istringstream cent(text_line);
        cent >> center[0]; cent >> center[1]; cent >> center[2];
        config->SetPeriodicCenter(iPeriodic, center);
        getline (mesh_file,text_line);
        istringstream rot(text_line);
        rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
        config->SetPeriodicRotation(iPeriodic, rotation);
        getline (mesh_file,text_line);
        istringstream tran(text_line);
        tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
        config->SetPeriodicTranslate(iPeriodic, translate);
        
      }
      
    }
    
  }
  
  /*--- If no periodic transormations were found, store default zeros ---*/
  if (!found_transform) {
    unsigned short nPeriodic = 1, iPeriodic = 0;
    config->SetnPeriodicIndex(nPeriodic);
    double* center    = new double[3];
    double* rotation  = new double[3];
    double* translate = new double[3];
    for (unsigned short iDim = 0; iDim < 3; iDim++) {
      center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
    }
    config->SetPeriodicCenter(iPeriodic, center);
    config->SetPeriodicRotation(iPeriodic, rotation);
    config->SetPeriodicTranslate(iPeriodic, translate);
  }
  
  /*--- Close the input file ---*/
  mesh_file.close();
  
  if (config->GetDivide_Element()) {
    if (Local2Global != NULL) delete [] Local2Global;
  }
  
}

void CPhysicalGeometry::Read_CGNS_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone){
  
  /*--- Original CGNS reader implementation by Thomas D. Economon,
   Francisco Palacios. Improvements for mixed-element meshes generated
   by ICEM added by Martin Spel (3D) & Shlomy Shitrit (2D), April 2014. ---*/
  
#ifdef HAVE_CGNS
  
  /*--- Local variables and initialization ---*/
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short VTK_Type, iMarker;
  unsigned long iPoint = 0, ielem_div = 0, ielem = 0, GlobalIndex;
  int rank = MASTER_NODE, size = SINGLE_NODE;
  nZone = val_nZone;
  
  /*--- Local variables which are needed when calling the CGNS mid-level API. ---*/
  unsigned long vnodes_cgns[8];
  double Coord_cgns[3];
  int fn, nbases = 0, nzones = 0, ngrids = 0, ncoords = 0, nsections = 0, file_type;
  int *vertices = NULL, *cells = NULL, nMarkers = 0, *boundVerts = NULL, npe;
  int interiorElems = 0, boundaryElems = 0, totalVerts = 0, prevVerts = 0;
  int cell_dim = 0, phys_dim = 0, nbndry, parent_flag;
  char basename[CGNS_STRING_SIZE], zonename[CGNS_STRING_SIZE];
  char coordname[CGNS_STRING_SIZE];
  cgsize_t* cgsize; cgsize = new cgsize_t[3];
  ZoneType_t zonetype;
  DataType_t datatype;
  double** coordArray = NULL;
  double*** gridCoords = NULL;
  ElementType_t elemType;
  cgsize_t range_min, range_max, startE, endE;
  range_min = 1;
  string currentElem;
  int** elemTypeVTK = NULL;
  int** elemIndex = NULL;
  int** nElems = NULL;
  int indexMax, elemMax; indexMax = elemMax = 0;
  cgsize_t**** connElems = NULL;
  cgsize_t* connElemTemp = NULL;
  cgsize_t ElementDataSize = 0;
  cgsize_t* parentData = NULL;
  int** dataSize = NULL;
  bool** isInternal = NULL;
  char*** sectionNames = NULL;
  
  /*--- Initialize counters for local/global points & elements ---*/
  FinestMGLevel = true;
  Global_nPoint  = 0; Global_nPointDomain = 0; Global_nElem = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_wedge    = 0; Global_nelem_wedge    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Check whether the supplied file is truly a CGNS file. ---*/
  if ( cg_is_cgns(val_mesh_filename.c_str(),&file_type) != CG_OK ) {
    printf( "\n\n   !!! Error !!!\n" );
    printf( " %s is not a CGNS file.\n", val_mesh_filename.c_str());
    printf( " Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }
  
  /*--- Open the CGNS file for reading. The value of fn returned
   is the specific index number for this file and will be
   repeatedly used in the function calls. ---*/
  if ( cg_open(val_mesh_filename.c_str(),CG_MODE_READ,&fn) ) cg_error_exit();
  cout << "Reading the CGNS file: " << val_mesh_filename.c_str() << endl;
  
  /*--- Get the number of databases. This is the highest node
   in the CGNS heirarchy. ---*/
  if ( cg_nbases(fn, &nbases) ) cg_error_exit();
  cout << "CGNS file contains " << nbases << " database(s)." << endl;
  
  /*--- Check if there is more than one database. Throw an
   error if there is because this reader can currently
   only handle one database. ---*/
  if ( nbases > 1 ) {
    printf("\n\n   !!! Error !!!\n" );
    printf("CGNS reader currently incapable of handling more than 1 database.");
    printf("Now exiting...\n\n");
    exit(EXIT_FAILURE);
  }
  
  /*--- Read the databases. Note that the indexing starts at 1. ---*/
  for ( int i = 1; i <= nbases; i++ ) {
    
    if ( cg_base_read(fn, i, basename, &cell_dim, &phys_dim) ) cg_error_exit();
    
    /*--- Get the number of zones for this base. ---*/
    if ( cg_nzones(fn, i, &nzones) ) cg_error_exit();
    cout << "Database " << i << ", " << basename << ": " << nzones;
    cout << " zone(s), cell dimension of " << cell_dim << ", physical ";
    cout << "dimension of " << phys_dim << "." << endl;
    
    /*--- Check if there is more than one zone. Throw an
     error if there is, because this reader can currently
     only handle one zone. This can/will be extended in the future. ---*/
    if ( nzones > 1 ) {
      printf("\n\n   !!! Error !!!\n" );
      printf("CGNS reader currently incapable of handling more than 1 zone.");
      printf("Now exiting...\n\n");
      exit(EXIT_FAILURE);
    }
    
    /*--- Initialize some data structures for  all zones. ---*/
    vertices   = new int[nzones];
    cells      = new int[nzones];
    boundVerts = new int[nzones];
    
    coordArray  = new double*[nzones];
    gridCoords  = new double**[nzones];
    elemTypeVTK = new int*[nzones];
    elemIndex   = new int*[nzones];
    nElems      = new int*[nzones];
    dataSize    = new int*[nzones];
    isInternal  = new bool*[nzones];
    nMarkers    = 0;
    
    sectionNames= new char**[nzones];
    connElems = new cgsize_t***[nzones];
    
    /*--- Loop over all zones in this base. Again, indexing starts at 1. ---*/
    for ( int j = 1; j <= nzones; j++ ) {
      
      /*--- Read the basic information for this zone, including
       the name and the number of vertices, cells, and
       boundary cells which are stored in the cgsize variable. ---*/
      if ( cg_zone_read(fn, i, j, zonename, cgsize) ) cg_error_exit();
      
      /*--- Rename the zone size information for clarity.
       NOTE: The number of cells here may be only the number of
       interior elements or it may be the total. This needs to
       be counted explicitly later. ---*/
      vertices[j-1]   = cgsize[0];
      cells[j-1]      = cgsize[1];
      boundVerts[j-1] = cgsize[2];
      
      /*--- Increment the total number of vertices from all zones. ---*/
      totalVerts += vertices[j-1];
      
      if ( cg_zone_type(fn, i, j, &zonetype) ) cg_error_exit();
      cout << "Zone " << j << ", " << zonename << ": " << vertices[j-1];
      cout << " vertices, " << cells[j-1] << " cells, " << boundVerts[j-1];
      cout << " boundary vertices." << endl;
      
      /*--- Retrieve the number of grids in this zone.
       For now, we know this is one, but to be more
       general, this will need to check and allow
       for a loop over all grids. ---*/
      
      if ( cg_ngrids(fn, i, j, &ngrids) ) cg_error_exit();
      
      /*--- Check the number of coordinate arrays stored
       in this zone. Should be 2 for 2-D grids and
       3 for 3-D grids. ---*/
      if ( ngrids > 1 ) {
        printf("\n\n   !!! Error !!!\n" );
        printf("CGNS reader currently handles only 1 grid per zone.");
        printf("Now exiting...\n\n");
        exit(EXIT_FAILURE);
      }
      
      if ( cg_ncoords( fn, i, j, &ncoords) ) cg_error_exit();
      cout << "Reading grid coordinates..." << endl;
      cout << "Number of coordinate dimensions is " << ncoords << "." << endl;
      
      /*--- Set the value of range_max to the total number
       of nodes in the unstructured mesh. Also allocate
       memory for the temporary array that will hold
       the grid coordinates as they are extracted. ---*/
      
      range_max       = cgsize[0];
      coordArray[j-1] = new double[range_max];
      
      /*--- Allocate memory for the 2-D array which will
       store the x, y, & z (if required) coordinates
       for writing into the SU2 mesh. ---*/
      
      gridCoords[j-1] = new double*[ncoords];
      for (int ii = 0; ii < ncoords; ii++) {
        *( gridCoords[j-1] + ii ) = new double[range_max];
      }
      
      /*--- Loop over each set of coordinates. Note again
       that the indexing starts at 1. ---*/
      
      for ( int k = 1; k <= ncoords; k++ ) {
        
        /*--- Read the coordinate info. This will retrieve the
         data type (either RealSingle or RealDouble) as
         well as the coordname which will specifiy the
         type of data that it is based in the SIDS convention.
         This might be "CoordinateX," for instance. ---*/
        
        if ( cg_coord_info(fn, i, j, k, &datatype, coordname) ) cg_error_exit();
        cout << "Reading " << coordname << " values from file." << endl;
        
        /*--- Always retrieve the grid coords in double precision. ---*/
        datatype = RealDouble;
        if ( cg_coord_read(fn, i, j, coordname, datatype, &range_min,
                           &range_max, coordArray[j-1]) ) cg_error_exit();
        
        /*--- Copy these coords into the 2-D array for storage until
         writing the SU2 mesh. ---*/
        
        for (int m = 0; m < range_max; m++ ) {
          gridCoords[j-1][k-1][m] = coordArray[j-1][m];
        }
        
      }
      
      /*--- Begin section for retrieving the connectivity info. ---*/
      
      cout << "Reading connectivity information..." << endl;
      
      /*--- First check the number of sections. ---*/
      
      if ( cg_nsections(fn, i, j, &nsections) ) cg_error_exit();
      cout << "Number of connectivity sections is " << nsections << "." << endl;
      
      /*--- Allocate several data structures to hold the various
       pieces of information describing each section. It is
       stored in this manner so that it can be written to
       SU2 memory later. ---*/
      
      elemTypeVTK[j-1] = new int[nsections];
      elemIndex[j-1]   = new int[nsections];
      nElems[j-1]      = new int[nsections];
      dataSize[j-1]    = new int[nsections];
      isInternal[j-1]  = new bool[nsections];
      
      sectionNames[j-1] = new char*[nsections];
      for (int ii = 0; ii < nsections; ii++) {
        sectionNames[j-1][ii]= new char[CGNS_STRING_SIZE];
      }
      
      /*--- Loop over each section. This will include the main
       connectivity information for the grid cells, as well
       as any boundaries which were labeled before export. ---*/
      
      for ( int s = 1; s <= nsections; s++ ) {
        
        /*--- Read the connectivity details for this section.
         Store the total number of elements in this section
         to be used later for memory allocation. ---*/
        
        if ( cg_section_read(fn, i, j, s, sectionNames[j-1][s-1], &elemType, &startE,
                             &endE, &nbndry, &parent_flag) ) cg_error_exit();
        nElems[j-1][s-1] = (int) (endE-startE+1);
        
        /*--- Read the total number of nodes that will be
         listed when reading this section. ---*/
        
        if ( cg_ElementDataSize(fn, i, j, s, &ElementDataSize) )
          cg_error_exit();
        dataSize[j-1][s-1] = ElementDataSize;
        
        /*--- Find the number of nodes required to represent
         this type of element. ---*/
        
        if ( cg_npe(elemType, &npe) ) cg_error_exit();
        elemIndex[j-1][s-1] = npe;
        
        /*--- Need to check the element type and correctly
         specify the VTK identifier for that element.
         SU2 recognizes elements by their VTK number. ---*/
        
        switch (elemType) {
          case NODE:
            currentElem      = "Vertex";
            elemTypeVTK[j-1][s-1] = 1;
            break;
          case BAR_2:
            currentElem      = "Line";
            elemTypeVTK[j-1][s-1] = 3;
            break;
          case BAR_3:
            currentElem      = "Line";
            elemTypeVTK[j-1][s-1] = 3;
            break;
          case TRI_3:
            currentElem      = "Triangle";
            elemTypeVTK[j-1][s-1] = 5;
            break;
          case QUAD_4:
            currentElem      = "Quadrilateral";
            elemTypeVTK[j-1][s-1] = 9;
            break;
          case TETRA_4:
            currentElem      = "Tetrahedron";
            elemTypeVTK[j-1][s-1] = 10;
            break;
          case HEXA_8:
            currentElem      = "Hexahedron";
            elemTypeVTK[j-1][s-1] = 12;
            break;
          case PENTA_6:
            currentElem      = "Prism";
            elemTypeVTK[j-1][s-1] = 13;
            break;
          case PYRA_5:
            currentElem      = "Pyramid";
            elemTypeVTK[j-1][s-1] = 14;
            break;
          case HEXA_20:
            printf( "\n\n   !!! Error !!!\n" );
            printf( " HEXA-20 element type not supported\n");
            printf(" Section %d, npe=%d\n", s, npe);
            printf(" startE %d, endE %d\n", startE, endE);
            printf( " Now exiting...\n\n");
            exit(EXIT_FAILURE);
            break;
          case MIXED:
            currentElem      = "Mixed";
            elemTypeVTK[j-1][s-1] = -1;
            break;
          default:
            printf( "\n\n   !!! Error !!!\n" );
            printf( " Unrecognized element type: (type %d, npe=%d)\n", elemType, npe);
            printf(" Section %d\n", s);
            printf(" startE %d, endE %d\n", startE, endE);
            printf( " Now exiting...\n\n");
            exit(EXIT_FAILURE);
            break;
        }
        
        /*--- Check if the elements in this section are part
         of the internal domain or are part of the boundary
         surfaces. This will be used later to separate the
         internal connectivity from the boundary connectivity.
         We will check for quad and tri elements for 3-D meshes
         because these will be the boundaries. Similarly, line
         elements will be boundaries to 2-D problems. ---*/
        
        if ( cell_dim == 2 ) {
          /*--- In 2-D check for line elements, VTK type 3. ---*/
          if (elemTypeVTK[j-1][s-1] == 3) {
            isInternal[j-1][s-1] = false;
            nMarkers++;
            boundaryElems += nElems[j-1][s-1];
          } else {
            isInternal[j-1][s-1] = true;
            interiorElems += nElems[j-1][s-1];
          }
        } else {
          /*--- In 3-D check for tri or quad elements, VTK types 5 or 9. ---*/
          switch (elemTypeVTK[j-1][s-1]) {
            case 5:
            case 9:
              isInternal[j-1][s-1] = false;
              nMarkers++;
              boundaryElems += nElems[j-1][s-1];
              break;
            case -1:
              /*--- MIXED element support (treated later). ---*/
              break;
            default:
              isInternal[j-1][s-1] = true;
              interiorElems += nElems[j-1][s-1];
              break;
          }
        }
        
        if (elemTypeVTK[j-1][s-1] == -1) {
          /*--- In case of mixed data type, allocate place for 8 nodes maximum
           (hex), plus element type. ---*/
          elemIndex[j-1][s-1] = 9;
        }
        
        /*--- Keep track of the sections with the largest
         number of elements and the number of nodes
         required to specify an element of a specific
         type. These max values will be used to allocate
         one large array. ---*/
        
        if ( elemIndex[j-1][s-1] > indexMax || s == 1 ) indexMax = elemIndex[j-1][s-1];
        if ( nElems[j-1][s-1] > elemMax || s == 1 )     elemMax  = nElems[j-1][s-1];
        
        /*--- Print some information to the console. ---*/
        
        cout << "Reading section " << sectionNames[j-1][s-1];
        cout << " of element type " << currentElem << "\n   starting at ";
        cout << startE << " and ending at " << endE << "." << endl;
        
      }
      
      /*--- Allocate memory to store all of the connectivity
       information in one large array. ---*/
      
      connElems[j-1] = new cgsize_t**[nsections];
      for (int ii = 0; ii < nsections; ii++) {
        connElems[j-1][ii] = new cgsize_t*[indexMax];
        for (int jj = 0; jj < indexMax; jj++) {
          connElems[j-1][ii][jj] = new cgsize_t[elemMax];
        }
      }
      
      for ( int s = 1; s <= nsections; s++ ) {
        
        connElemTemp = new cgsize_t[dataSize[j-1][s-1]];
        
        /*--- Retrieve the connectivity information and store. ---*/
        
        if ( cg_elements_read(fn, i, j, s, connElemTemp, parentData) )
          cg_error_exit();
        
        /*--- Copy these values into the larger array for
         storage until writing the SU2 file. ---*/
        
        if (elemTypeVTK[j-1][s-1] == -1) {
          
          /*--- In the zone, we assume we don't mix volumetric and surface
           elements, so if at least one surface element is found (TRIA, QUAD),
           all elements are assumed to be surface elements. Mixed-element
           support initially added here by Martin Spel. ---*/
          
          bool isBoundary = false;
          int counter = 0;
          
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
            ElementType_t elmt_type = ElementType_t (connElemTemp[counter]);
            cg_npe( elmt_type, &npe);
            
            /*--- Mixed element support for 2D added here by Shlomy Shitrit ---*/
            if (elmt_type == 5 || elmt_type == 9){
              if (cell_dim == 2) { isBoundary = false; }
              if (cell_dim == 3) { isBoundary = true;  }
            }
            
            counter++;
            connElems[j-1][s-1][0][ii] = elmt_type;
            for ( int jj = 0; jj < npe; jj++ ) {
              connElems[j-1][s-1][jj+1][ii] = connElemTemp[counter] + prevVerts;
              counter++;
            }
          }
          if (isBoundary) {
            isInternal[j-1][s-1] = false;
            nMarkers++;
            boundaryElems += nElems[j-1][s-1];
          } else if ( cell_dim == 3 ) {
            isInternal[j-1][s-1] = true;
            interiorElems += nElems[j-1][s-1];
          }
          
          
        } else {
          int counter = 0;
          for ( int ii = 0; ii < nElems[j-1][s-1]; ii++ ) {
            for ( int jj = 0; jj < elemIndex[j-1][s-1]; jj++ ) {
              connElems[j-1][s-1][jj][ii] = connElemTemp[counter] + prevVerts;
              counter++;
            }
          }
        }
        delete[] connElemTemp;
      }
      prevVerts += vertices[j-1];
      
    }
  }
  
  /*--- Close the CGNS file. ---*/
  
  if ( cg_close(fn) ) cg_error_exit();
  cout << "Successfully closed the CGNS file." << endl;
  
  /*--- Write a SU2 mesh if requested in the config file. ---*/
  if (config->GetCGNS_To_SU2()) {
    
    string fileNameSU2 = config->GetMesh_Out_FileName();
    cout << "Writing SU2 mesh file: " << fileNameSU2 << "." << endl;
    
    /*--- Open the solution file for writing. ---*/
    
    FILE *SU2File;
    if ( (SU2File = fopen(fileNameSU2.c_str(), "w")) != NULL ) {
      
      /*--- Write the dimension of the problem and the
       total number of elements first. Note that we
       need to use the interior elements here (not "cells"). ---*/
      
      fprintf( SU2File, "NDIME= %i\n", cell_dim);
      fprintf( SU2File, "NELEM= %i\n", interiorElems);
      
      /*--- Connectivity info for the internal domain. ---*/
      
      int counter = 0;
      for ( int k = 0; k < nzones; k++ ) {
        for ( int s = 0; s < nsections; s++ ) {
          if ( isInternal[k][s] ) {
            for ( int i = 0; i < nElems[k][s]; i++ ) {
              if (elemTypeVTK[k][s] == -1 ) {
                
                /*--- Mixed-element support. ---*/
                ElementType_t elmt_type = ElementType_t (connElems[k][s][0][i]);
                cg_npe( elmt_type, &npe);
                switch (elmt_type) {
                  case NODE:
                    fprintf( SU2File, "%2i\t", 1);
                    break;
                  case BAR_2:
                    fprintf( SU2File, "%2i\t", 3);
                    break;
                  case BAR_3:
                    fprintf( SU2File, "%2i\t", 3);
                    break;
                  case TRI_3:
                    fprintf( SU2File, "%2i\t", 5);
                    break;
                  case QUAD_4:
                    fprintf( SU2File, "%2i\t", 9);
                    break;
                  case TETRA_4:
                    fprintf( SU2File, "%2i\t", 10);
                    break;
                  case HEXA_8:
                    fprintf( SU2File, "%2i\t", 12);
                    break;
                  case PENTA_6:
                    fprintf( SU2File, "%2i\t", 13);
                    break;
                  case PYRA_5:
                    fprintf( SU2File, "%2i\t", 14);
                    break;
                  default:
                    fprintf( SU2File, "%2i\t", -1);
                    break;
                }
                for ( int j = 1; j < npe+1; j++ ) {
                  fprintf( SU2File, "%8i\t", connElems[k][s][j][i] - 1);
                }
                fprintf( SU2File, "%d\n", counter);
                counter++;
              } else {
                fprintf( SU2File, "%2i\t", elemTypeVTK[k][s]);
                for ( int j = 0; j < elemIndex[k][s]; j++ ) {
                  fprintf( SU2File, "%8i\t", connElems[k][s][j][i] - 1);
                }
                fprintf( SU2File, "%d\n", counter);
                counter++;
              }
            }
          }
        }
      }
      
      /*--- Now write the node coordinates. First write
       the total number of vertices. ---*/
      
      fprintf( SU2File, "NPOIN= %i\n", totalVerts);
      counter = 0;
      for ( int k = 0; k < nzones; k ++ ) {
        for ( int i = 0; i < vertices[k]; i++ ) {
          for ( int j = 0; j < cell_dim; j++ ) {
            fprintf( SU2File, "%.16le\t", gridCoords[k][j][i]);
          }
          fprintf( SU2File, "%d\n", counter);
          counter++;
        }
      }
      
      /*--- Lastly write the boundary information.
       These will write out in the same order
       that they are stored in the cgns file.
       Note that the connectivity
       values are decremented by 1 in order to
       match the indexing in SU2.              ---*/
      
      fprintf( SU2File, "NMARK= %i\n", nMarkers );
      counter = 0;
      for ( int k = 0; k < nzones; k ++ ) {
        for ( int s = 0; s < nsections; s++ ) {
          if ( !isInternal[k][s] ) {
            fprintf( SU2File, "MARKER_TAG= %s\n", sectionNames[k][s] );
            fprintf( SU2File, "MARKER_ELEMS= %i\n", nElems[k][s] );
            counter++;
            
            if (elemTypeVTK[k][s] == -1 ) {
              
              /*--- Mixed-element support. ---*/
              
              for ( int i = 0; i < nElems[k][s]; i++ ) {
                ElementType_t elmt_type = ElementType_t (connElems[k][s][0][i]);
                cg_npe( elmt_type, &npe);
                switch (elmt_type) {
                  case NODE:
                    fprintf( SU2File, "%2i\t", 1);
                    break;
                  case BAR_2:
                    fprintf( SU2File, "%2i\t", 3);
                    break;
                  case BAR_3:
                    fprintf( SU2File, "%2i\t", 3);
                    break;
                  case TRI_3:
                    fprintf( SU2File, "%2i\t", 5);
                    break;
                  case QUAD_4:
                    fprintf( SU2File, "%2i\t", 9);
                    break;
                  case TETRA_4:
                    fprintf( SU2File, "%2i\t", 10);
                    break;
                  case HEXA_8:
                    fprintf( SU2File, "%2i\t", 12);
                    break;
                  case PENTA_6:
                    fprintf( SU2File, "%2i\t", 13);
                    break;
                  case PYRA_5:
                    fprintf( SU2File, "%2i\t", 14);
                    break;
                  default: // error
                    fprintf( SU2File, "%2i\t", -1);
                    break;
                }
                for ( int j = 1; j < npe+1; j++ ) {
                  fprintf( SU2File, "%8i\t", connElems[k][s][j][i] - 1);
                }
                fprintf( SU2File, "\n");
              }
            } else {
              counter++;
              for ( int i = 0; i < nElems[k][s]; i++ ) {
                fprintf( SU2File, "%2i\t", elemTypeVTK[k][s]);
                for ( int j = 0; j < elemIndex[k][s]; j++ ) {
                  fprintf( SU2File, "%8i\t", connElems[k][s][j][i] - 1 );
                }
                fprintf( SU2File, "\n");
              }
            }
          }
        }
      }
    }
    
    /*--- Close the SU2 mesh file. ---*/
    
    fclose( SU2File );
    cout << "Successfully wrote the SU2 mesh file." << endl;
    
  }
  
  /*--- Load the data from the CGNS file into SU2 memory. ---*/
  
  if (rank == MASTER_NODE)
    cout << endl << "Loading CGNS data into SU2 data structures..." << endl;
  
  /*--- Read the dimension of the problem ---*/
  
  nDim = cell_dim;
  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Two dimensional problem." << endl;
    if (nDim == 3) cout << "Three dimensional problem." << endl;
  }
  
  /*--- Read the information about inner elements ---*/
  
  nElem = interiorElems;
  cout << nElem << " interior elements." << endl;
  Global_nElem = nElem;
  
  /*--- Allocate space for elements ---*/
  
  elem = new CPrimalGrid*[nElem];
  
  /*--- Loop over all the volumetric elements ---*/
  for ( int k = 0; k < nzones; k ++ ) {
    for ( int s = 0; s < nsections; s++ ) {
      if ( isInternal[k][s] ) {
        for ( int i = 0; i < nElems[k][s]; i++ ) {
          
          /*--- Get the VTK type for this element. Check for mixed
           elements. ---*/
          
          if (elemTypeVTK[k][s] == -1 ) {
            
            /*--- Mixed-element support. ---*/
            ElementType_t elmt_type = ElementType_t (connElems[k][s][0][i]);
            cg_npe( elmt_type, &npe);
            switch (elmt_type) {
              case NODE:
                VTK_Type = 1;
                break;
              case BAR_2:
                VTK_Type = 3;
                break;
              case BAR_3:
                VTK_Type = 3;
                break;
              case TRI_3:
                VTK_Type = 5;
                break;
              case QUAD_4:
                VTK_Type = 9;
                break;
              case TETRA_4:
                VTK_Type = 10;
                break;
              case HEXA_8:
                VTK_Type = 12;
                break;
              case PENTA_6:
                VTK_Type = 13;
                break;
              case PYRA_5:
                VTK_Type = 14;
                break;
              default: // error
                cout << "Kind of element not suppported!" << endl;
                break;
            }
            
            /*--- Transfer the nodes for this element. ---*/
            for ( int j = 1; j < npe+1; j++ ) {
              vnodes_cgns[j-1] = connElems[k][s][j][i] - 1;
            }
            
          } else {
            VTK_Type = elemTypeVTK[k][s];
            
            /*--- Transfer the nodes for this element. ---*/
            for ( int j = 0; j < elemIndex[k][s]; j++ ) {
              vnodes_cgns[j] = connElems[k][s][j][i] - 1;
            }
          }
          
          /* Instantiate this element. */
          switch(VTK_Type) {
            case TRIANGLE:
              elem[ielem] = new CTriangle(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],nDim);
              ielem_div++; ielem++; nelem_triangle++; break;
            case RECTANGLE:
              elem[ielem] = new CRectangle(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3],nDim);
              ielem++; nelem_quad++;
              ielem_div++;
              break;
            case TETRAHEDRON:
              elem[ielem] = new CTetrahedron(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3]);
              ielem_div++; ielem++; nelem_tetra++; break;
            case HEXAHEDRON:
              elem[ielem] = new CHexahedron(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3],vnodes_cgns[4],vnodes_cgns[5],vnodes_cgns[6],vnodes_cgns[7]);
              ielem++; nelem_hexa++; ielem_div++;
              break;
            case WEDGE:
              elem[ielem] = new CWedge(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3],vnodes_cgns[4],vnodes_cgns[5]);
              ielem++; nelem_wedge++; ielem_div++;
              break;
            case PYRAMID:
              elem[ielem] = new CPyramid(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3],vnodes_cgns[4]);
              ielem++; nelem_pyramid++; ielem_div++;
              break;
          }
        }
      }
    }
  }
  
  if (config->GetDivide_Element()) nElem = nelem_triangle + nelem_quad + nelem_tetra + nelem_hexa + nelem_wedge + nelem_pyramid;
  
  Global_nelem_triangle = nelem_triangle;
  Global_nelem_quad     = nelem_quad;
  Global_nelem_tetra    = nelem_tetra;
  Global_nelem_hexa     = nelem_hexa;
  Global_nelem_wedge    = nelem_wedge;
  Global_nelem_pyramid  = nelem_pyramid;
  
  /*--- Print information about the elements to the console ---*/
  if (size == 1) {
    if (Global_nelem_triangle > 0)
      cout << Global_nelem_triangle << " triangles." << endl;
    if (Global_nelem_quad > 0)
      cout << Global_nelem_quad << " quadrilaterals." << endl;
    if (Global_nelem_tetra > 0)
      cout << Global_nelem_tetra << " tetrahedra." << endl;
    if (Global_nelem_hexa > 0)
      cout << Global_nelem_hexa << " hexahedra." << endl;
    if (Global_nelem_wedge > 0)
      cout << Global_nelem_wedge << " prisms." << endl;
    if (Global_nelem_pyramid > 0)
      cout << Global_nelem_pyramid << " pyramids." << endl;
  }
  
  /*--- Read node coordinates. Note this assumes serial mode. ---*/
  nPoint = totalVerts;
  nPointDomain = nPoint;
  cout << nPoint << " points." << endl;
  node = new CPoint*[nPoint];
  for ( int k = 0; k < nzones; k++ ) {
    for ( int i = 0; i < vertices[k]; i++ ) {
      for ( int j = 0; j < cell_dim; j++ ) {
        Coord_cgns[j] = gridCoords[k][j][i];
      }
      switch(nDim) {
        case 2:
          GlobalIndex = i;
          node[iPoint] = new CPoint(Coord_cgns[0], Coord_cgns[1], GlobalIndex, config);
          iPoint++; break;
        case 3:
          GlobalIndex = i;
          node[iPoint] = new CPoint(Coord_cgns[0], Coord_cgns[1], Coord_cgns[2], GlobalIndex, config);
          iPoint++; break;
      }
    }
  }
  
  /*--- Set some important point information for parallel simulations. ---*/
  Global_nPoint = nPoint;
  Global_nPointDomain = nPointDomain;
  
  /*--- Read number of markers ---*/
  nMarker = nMarkers;
  cout << nMarker << " surface markers." << endl;
  config->SetnMarker_All(nMarker);
  bound = new CPrimalGrid**[nMarker];
  nElem_Bound = new unsigned long [nMarker];
  Tag_to_Marker = new string [MAX_NUMBER_MARKER];
  
  iMarker = 0;
  for ( int k = 0; k < nzones; k ++ ) {
    for ( int s = 0; s < nsections; s++ ) {
      if ( !isInternal[k][s] ) {
        nelem_edge_bound = 0; nelem_triangle_bound = 0; nelem_quad_bound = 0; ielem = 0;
        Marker_Tag = sectionNames[k][s];
        
        /*--- Remove whitespaces from the marker names ---*/
        Marker_Tag.erase(remove(Marker_Tag.begin(),Marker_Tag.end(),' '),Marker_Tag.end());
        
        if (Marker_Tag != "SEND_RECEIVE") {
          nElem_Bound[iMarker] = nElems[k][s];
          if (rank == MASTER_NODE)
            cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
          bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          
          for ( int i = 0; i < nElems[k][s]; i++ ) {
            
            /*--- Get the VTK type for this element. Check for mixed
             elements. ---*/
            
            if (elemTypeVTK[k][s] == -1 ) {
              
              /*--- Mixed-element support. ---*/
              ElementType_t elmt_type = ElementType_t (connElems[k][s][0][i]);
              cg_npe( elmt_type, &npe);
              switch (elmt_type) {
                case NODE: VTK_Type = 1; break;
                case BAR_2: VTK_Type = 3; break;
                case BAR_3: VTK_Type = 3; break;
                case TRI_3: VTK_Type = 5; break;
                case QUAD_4: VTK_Type = 9; break;
                case TETRA_4: VTK_Type = 10; break;
                case HEXA_8: VTK_Type = 12; break;
                case PENTA_6: VTK_Type = 13; break;
                case PYRA_5: VTK_Type = 14; break;
                default: cout << "Kind of element not suppported!" << endl; break;  // error
              }
              /*--- Transfer the nodes for this element. ---*/
              for ( int j = 1; j < npe+1; j++ ) {
                vnodes_cgns[j-1] = connElems[k][s][j][i] - 1;
              }
              
            }else {
              VTK_Type = elemTypeVTK[k][s];
              
              /* Transfer the nodes for this element. */
              for ( int j = 0; j < elemIndex[k][s]; j++ ) {
                vnodes_cgns[j] = connElems[k][s][j][i] - 1;
              }
            }
            
            
            switch(VTK_Type) {
              case LINE:
                
                if (nDim == 3) {
                  cout << "Please remove line boundary conditions from the mesh file!" << endl;
#ifndef HAVE_MPI
                  exit(EXIT_FAILURE);
#else
                  MPI_Abort(MPI_COMM_WORLD,1);
                  MPI_Finalize();
#endif
                }
                
                bound[iMarker][ielem] = new CLine(vnodes_cgns[0],vnodes_cgns[1],2);
                ielem++; nelem_edge_bound++; break;
              case TRIANGLE:
                bound[iMarker][ielem] = new CTriangle(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],3);
                ielem++; nelem_triangle_bound++; break;
              case RECTANGLE:
                bound[iMarker][ielem] = new CRectangle(vnodes_cgns[0],vnodes_cgns[1],vnodes_cgns[2],vnodes_cgns[3],3);
                ielem++; nelem_quad_bound++; break;
            }
          }
          
          /*--- Update config information storing the boundary information in the right place ---*/
          Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
          config->SetMarker_All_TagBound(iMarker, Marker_Tag);
          config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
          config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
          config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
          config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
          config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
          config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
          config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
          config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
          config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
          config->SetMarker_All_SendRecv(iMarker, NONE);

        }
        iMarker++;
      }
    }
  }
  
  /*--- Periodic transormations is not implement, store default zeros ---*/
  unsigned short nPeriodic = 1, iPeriodic = 0;
  config->SetnPeriodicIndex(nPeriodic);
  double* center    = new double[3];
  double* rotation  = new double[3];
  double* translate = new double[3];
  for (unsigned short iDim = 0; iDim < 3; iDim++) {
    center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
  }
  config->SetPeriodicCenter(iPeriodic, center);
  config->SetPeriodicRotation(iPeriodic, rotation);
  config->SetPeriodicTranslate(iPeriodic, translate);
  
  /*--- Deallocate temporary memory. ---*/
  delete[] vertices;
  delete[] cells;
  delete[] boundVerts;
  
  for ( int j = 0; j < nzones; j++) {
    delete[] coordArray[j];
    delete[] elemTypeVTK[j];
    delete[] elemIndex[j];
    delete[] nElems[j];
    delete[] dataSize[j];
    delete[] isInternal[j];
    delete[] sectionNames[j];
  }
  delete[] coordArray;
  delete[] elemTypeVTK;
  delete[] elemIndex;
  delete[] nElems;
  delete[] dataSize;
  delete[] isInternal;
  delete[] sectionNames;
  
  for ( int j = 0; j < nzones; j++) {
    for( int i = 0; i < ncoords; i++ ) {
      delete[] gridCoords[j][i];
    }
    delete[] gridCoords[j];
  }
  delete[] gridCoords;
  
  for ( int kk = 0; kk < nzones; kk++) {
    for (int ii = 0; ii < nsections; ii++) {
      for (int jj = 0; jj < indexMax; jj++) {
        delete[] connElems[kk][ii][jj];
      }
      delete connElems[kk][ii];
    }
    delete connElems[kk];
  }
  delete[] connElems;
  
#else
  cout << "SU2 built without CGNS support!!" << endl;
  cout << "To use CGNS, remove the -DNO_CGNS directive ";
  cout << "from the makefile and supply the correct path";
  cout << " to the CGNS library." << endl;
  exit(EXIT_FAILURE);
#endif
  
}

void CPhysicalGeometry::Read_NETCDF_Format(CConfig *config, string val_mesh_filename, unsigned short val_iZone, unsigned short val_nZone) {
  
  /*--- Local variables and initialization ---*/
  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short VTK_Type, iMarker;
  unsigned long ielem = 0,
  vnodes_triangle[3], vnodes_quad[4], vnodes_tetra[4], vnodes_hexa[8], vnodes_wedge[6];
  char cstr[MAX_STRING_SIZE];
  string::size_type position;
  nZone = val_nZone;
  
  /*--- Initialize counters for local/global points & elements ---*/
  
  FinestMGLevel = true;
  Global_nPoint = 0; Global_nPointDomain = 0; Global_nElem = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_wedge    = 0; Global_nelem_wedge    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;
  
  /*--- Throw error if not in serial mode. ---*/
#ifdef HAVE_MPI
  cout << "Parallel support with NETCDF format not yet implemented!!" << endl;
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif
  
  unsigned short Marker_Index, marker, icommas, iDim;
  unsigned long nmarker, ielem_triangle, ielem_hexa, ncoord, iSurfElem, ielem_wedge,
  ielem_quad, *marker_list, **surf_elem, ielem_surface, *surf_marker, nSurfElem;
  double coord;
  string::size_type position_;
  bool stop, add;
  
  ielem_surface = 0; nSurfElem = 0;
  surf_marker = NULL;
  surf_elem = NULL;
  
  
  nDim = 3; cout << "Three dimensional problem." << endl;
  
  /*--- Open grid file ---*/
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  if (mesh_file.fail()) {
    cout << "There is no geometry file (CPhysicalGeometry)!!" << endl;
    exit(EXIT_FAILURE);
  }
  
  while (getline (mesh_file, text_line)) {
    
    position = text_line.find ("no_of_elements = ",0);
    if (position != string::npos) {
      text_line.erase (0,17); nElem = atoi(text_line.c_str());
      cout << nElem << " inner elements to store." << endl;
      elem = new CPrimalGrid*[nElem]; }
    
    position = text_line.find ("no_of_surfaceelements = ",0);
    if (position != string::npos) {
      text_line.erase (0,24); nSurfElem = atoi(text_line.c_str());
      cout << nSurfElem << " surface elements to store." << endl;
      surf_elem = new unsigned long* [nSurfElem];
      for (ielem_surface = 0; ielem_surface < nSurfElem; ielem_surface++)
        surf_elem[ielem_surface] = new unsigned long [5];
      ielem_surface = 0;
    }
    
    position = text_line.find ("no_of_points = ",0);
    if (position != string::npos) {
      text_line.erase (0,15); nPoint = atoi(text_line.c_str()); nPointDomain = nPoint;
      cout << nPoint << " points to store." << endl;
      node = new CPoint*[nPoint]; }
    
    position = text_line.find ("no_of_tetraeders = ",0);
    if (position != string::npos) {
      text_line.erase (0,19); nelem_tetra = atoi(text_line.c_str());
      cout << nelem_tetra << " tetraeders elements to store." << endl; }
    
    position = text_line.find ("no_of_prisms = ",0);
    if (position != string::npos) {
      text_line.erase (0,15); nelem_wedge = atoi(text_line.c_str());
      cout << nelem_wedge << " prims elements to store." << endl; }
    
    position = text_line.find ("no_of_hexaeders = ",0);
    if (position != string::npos) {
      text_line.erase (0,18); nelem_hexa = atoi(text_line.c_str());
      cout << nelem_hexa << " hexaeders elements to store." << endl; }
    
    position = text_line.find ("no_of_surfacetriangles = ",0);
    if (position != string::npos) {
      text_line.erase (0,25); nelem_triangle = atoi(text_line.c_str());
      cout << nelem_triangle << " surface triangle elements to store." << endl; }
    
    position = text_line.find ("no_of_surfacequadrilaterals = ",0);
    if (position != string::npos) {
      text_line.erase (0,30); nelem_quad = atoi(text_line.c_str());
      cout << nelem_quad << " surface quadrilaterals elements to store." << endl; }
    
    position = text_line.find ("points_of_tetraeders =",0);
    if (position != string::npos) {
      for (unsigned long ielem_tetra = 0; ielem_tetra < nelem_tetra; ielem_tetra++) {
        getline(mesh_file,text_line);
        for (unsigned short icommas = 0; icommas < 15; icommas++) {
          position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
        }
        istringstream elem_line(text_line);
        VTK_Type = TETRAHEDRON;
        elem_line >> vnodes_tetra[0]; elem_line >> vnodes_tetra[1]; elem_line >> vnodes_tetra[2]; elem_line >> vnodes_tetra[3];
        elem[ielem] = new CTetrahedron(vnodes_tetra[0],vnodes_tetra[1],vnodes_tetra[2],vnodes_tetra[3]);
        ielem++;
      }
      cout << "finish tetrahedron element reading" << endl;
    }
    
    position = text_line.find ("points_of_prisms =",0);
    if (position != string::npos) {
      for (ielem_wedge = 0; ielem_wedge < nelem_wedge; ielem_wedge++) {
        getline(mesh_file,text_line);
        for (icommas = 0; icommas < 15; icommas++) {
          position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
        }
        istringstream elem_line(text_line);
        VTK_Type = WEDGE;
        elem_line >> vnodes_wedge[0]; elem_line >> vnodes_wedge[1]; elem_line >> vnodes_wedge[2]; elem_line >> vnodes_wedge[3]; elem_line >> vnodes_wedge[4]; elem_line >> vnodes_wedge[5];
        elem[ielem] = new CWedge(vnodes_wedge[0],vnodes_wedge[1],vnodes_wedge[2],vnodes_wedge[3],vnodes_wedge[4],vnodes_wedge[5]);
        ielem++;
      }
      cout << "finish prims element reading" << endl;
    }
    
    position = text_line.find ("points_of_hexaeders =",0);
    if (position != string::npos) {
      for (ielem_hexa = 0; ielem_hexa < nelem_hexa; ielem_hexa++) {
        getline(mesh_file,text_line);
        for (icommas = 0; icommas < 15; icommas++) {
          position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
        }
        istringstream elem_line(text_line);
        VTK_Type = HEXAHEDRON;
        elem_line >> vnodes_hexa[0]; elem_line >> vnodes_hexa[1]; elem_line >> vnodes_hexa[2]; elem_line >> vnodes_hexa[3];
        elem_line >> vnodes_hexa[4]; elem_line >> vnodes_hexa[5]; elem_line >> vnodes_hexa[6]; elem_line >> vnodes_hexa[7];
        elem[ielem] = new CHexahedron(vnodes_hexa[0],vnodes_hexa[1],vnodes_hexa[2],vnodes_hexa[3],vnodes_hexa[4],vnodes_hexa[5],vnodes_hexa[6],vnodes_hexa[7]);
        ielem++;
      }
      cout << "finish hexaeders element reading" << endl;
    }
    
    position = text_line.find ("points_of_surfacetriangles =",0);
    if (position != string::npos) {
      for (ielem_triangle = 0; ielem_triangle < nelem_triangle; ielem_triangle++) {
        getline(mesh_file,text_line);
        for (icommas = 0; icommas < 15; icommas++) {
          position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
        }
        istringstream elem_line(text_line);
        elem_line >> vnodes_triangle[0]; elem_line >> vnodes_triangle[1]; elem_line >> vnodes_triangle[2];
        surf_elem[ielem_surface][0]= 3;
        surf_elem[ielem_surface][1]= vnodes_triangle[0];
        surf_elem[ielem_surface][2]= vnodes_triangle[1];
        surf_elem[ielem_surface][3]= vnodes_triangle[2];
        ielem_surface++;
      }
      cout << "finish surface triangles element reading" << endl;
    }
    
    position = text_line.find ("points_of_surfacequadrilaterals =",0);
    if (position != string::npos) {
      for (ielem_quad = 0; ielem_quad < nelem_quad; ielem_quad++) {
        getline(mesh_file,text_line);
        for (icommas = 0; icommas < 15; icommas++) {
          position_ = text_line.find( ",", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 ); if(position_!=string::npos) text_line.erase (position_,1);
        }
        istringstream elem_line(text_line);
        elem_line >> vnodes_quad[0]; elem_line >> vnodes_quad[1]; elem_line >> vnodes_quad[2]; elem_line >> vnodes_quad[3];
        surf_elem[ielem_surface][0]= 4;
        surf_elem[ielem_surface][1]= vnodes_quad[0];
        surf_elem[ielem_surface][2]= vnodes_quad[1];
        surf_elem[ielem_surface][3]= vnodes_quad[2];
        surf_elem[ielem_surface][4]= vnodes_quad[3];
        ielem_surface++;
      }
      cout << "finish surface quadrilaterals element reading" << endl;
    }
    
    position = text_line.find ("boundarymarker_of_surfaces =",0);
    if (position != string::npos) {
      nmarker=0;
      stop = false;
      surf_marker = new unsigned long [nelem_triangle + nelem_quad];
      
      text_line.erase (0,29);
      for (icommas = 0; icommas < 50; icommas++) {
        position_ = text_line.find( ",", 0 );
        if(position_!=string::npos) text_line.erase (position_,1);
      }
      
      stringstream  point_line(text_line);
      while (point_line >> marker,!point_line.eof()) {
        surf_marker[nmarker] = marker;
        nmarker++; }
      
      while (!stop) {
        getline(mesh_file,text_line);
        for (icommas = 0; icommas < 50; icommas++) {
          position_ = text_line.find( ",", 0 );
          if(position_!=string::npos) text_line.erase (position_,1);
          position_ = text_line.find( ";", 0 );
          if(position_!=string::npos) text_line.erase (position_,1);
        }
        stringstream  point_line(text_line);
        while (point_line>> marker,!point_line.eof()) {
          surf_marker[nmarker] = marker;
          if (nmarker == nSurfElem-1) {stop = true; break;}
          nmarker++;
        }
      }
    }
    
    for (iDim = 0; iDim < nDim; iDim++) {
      ncoord = 0; stop = false;
      if (iDim == 0) position = text_line.find ("points_xc = ",0);
      if (iDim == 1) position = text_line.find ("points_yc = ",0);
      if (iDim == 2) position = text_line.find ("points_zc = ",0);
      
      if (position != string::npos) {
        text_line.erase (0,12);
        for (icommas = 0; icommas < 50; icommas++) {
          position_ = text_line.find( ",", 0 );
          if(position_!=string::npos) text_line.erase (position_,1);
        }
        stringstream  point_line(text_line);
        while (point_line>> coord,!point_line.eof()) {
          if (iDim==0) node[ncoord] = new CPoint(coord, 0.0, 0.0, ncoord, config);
          if (iDim==1) node[ncoord]->SetCoord(1, coord);
          if (iDim==2) node[ncoord]->SetCoord(2, coord);
          ncoord++; }
        while (!stop) {
          getline(mesh_file,text_line);
          for (icommas = 0; icommas < 50; icommas++) {
            position_ = text_line.find( ",", 0 );
            if(position_!=string::npos) text_line.erase (position_,1);
            position_ = text_line.find( ";", 0 );
            if(position_!=string::npos) text_line.erase (position_,1);
          }
          stringstream  point_line(text_line);
          while (point_line>> coord,!point_line.eof()) {
            if (iDim==0) node[ncoord] = new CPoint(coord, 0.0, 0.0, ncoord, config);
            if (iDim==1) node[ncoord]->SetCoord(1, coord);
            if (iDim==2) node[ncoord]->SetCoord(2, coord);
            if (ncoord == nPoint-1) {stop = true; break;}
            ncoord++;
          }
        }
        if (iDim==0) cout << "finish point xc reading" << endl;
        if (iDim==1) cout << "finish point yc reading" << endl;
        if (iDim==2) cout << "finish point zc reading" << endl;
      }
    }
  }
  
  
  /*--- Create a list with all the markers ---*/
  marker_list = new unsigned long [MAX_NUMBER_MARKER];
  marker_list[0] = surf_marker[0]; nMarker = 1;
  for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++) {
    add = true;
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (marker_list[iMarker] == surf_marker[iSurfElem]) {
        add = false; break; }
    if (add) {
      marker_list[nMarker] = surf_marker[iSurfElem];
      nMarker++;
    }
  }
  
  nElem_Bound = new unsigned long [nMarker];
  
  /*--- Compute the number of element per marker ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElem_Bound[iMarker] = 0;
    for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++)
      if (surf_marker[iSurfElem] == marker_list[iMarker])
        nElem_Bound[iMarker]++;
  }
  
  
  /*--- Realate the marker index with the position in the array of markers ---*/
  unsigned short *Index_to_Marker;
  Index_to_Marker = new unsigned short [MAX_NUMBER_MARKER];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Marker_Index = marker_list[iMarker];
    Index_to_Marker[Marker_Index] = iMarker;
  }
  
  bound = new CPrimalGrid**[nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Marker_Index = marker_list[iMarker];
    bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
    ielem_triangle = 0; ielem_quad = 0;
    for (iSurfElem = 0; iSurfElem < nSurfElem; iSurfElem++)
      if (surf_marker[iSurfElem] == Marker_Index) {
        if (surf_elem[iSurfElem][0] == 3) {
          vnodes_triangle[0] = surf_elem[iSurfElem][1]; vnodes_triangle[1] = surf_elem[iSurfElem][2]; vnodes_triangle[2] = surf_elem[iSurfElem][3];
          bound[iMarker][ielem_triangle+ielem_quad] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
          ielem_triangle ++;
        }
        if (surf_elem[iSurfElem][0] == 4) {
          vnodes_quad[0] = surf_elem[iSurfElem][1]; vnodes_quad[1] = surf_elem[iSurfElem][2]; vnodes_quad[2] = surf_elem[iSurfElem][3]; vnodes_quad[3] = surf_elem[iSurfElem][4];
          bound[iMarker][ielem_triangle+ielem_quad] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
          ielem_quad ++;
        }
      }
  }
  
  
  
  /*--- Update config information storing the boundary information in the right place ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    stringstream out;
    out << marker_list[iMarker];
    Marker_Tag = out.str();
    
    Tag_to_Marker = new string [MAX_NUMBER_MARKER];
    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
    config->SetMarker_All_TagBound(iMarker, Marker_Tag);
    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
    config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
    config->SetMarker_All_SendRecv(iMarker, NONE);

  }
  
}

void CPhysicalGeometry::Check_IntElem_Orientation(CConfig *config) {
  
  unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6,
  iElem;
  double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, *Coord_6, a[3], b[3], c[3], n[3], test;
  unsigned short iDim;
  
  /*--- Loop over all the elements ---*/
  
  for (iElem = 0; iElem < nElem; iElem++) {
    
    /*--- 2D grid, triangle case ---*/
    
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test = a[0]*b[1]-b[0]*a[1];
      
      if (test < 0.0) elem[iElem]->Change_Orientation();
    }
    
    /*--- 2D grid, rectangle case ---*/
    
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test_1 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
        b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]); }
      test_2 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
        b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]); }
      test_3 = a[0]*b[1]-b[0]*a[1];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]); }
      test_4 = a[0]*b[1]-b[0]*a[1];
      
      if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0))
        elem[iElem]->Change_Orientation();
    }
    
    /*--- 3D grid, tetrahedron case ---*/
    
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      if (test < 0.0) elem[iElem]->Change_Orientation();
      
    }
    
    /*--- 3D grid, wedge case ---*/
    
    if (elem[iElem]->GetVTK_Type() == WEDGE) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
      Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        c[iDim] = (Coord_4[iDim]-Coord_1[iDim])+
        (Coord_5[iDim]-Coord_2[iDim])+
        (Coord_6[iDim]-Coord_3[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_5[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_6[iDim]-Coord_4[iDim]);
        c[iDim] = (Coord_1[iDim]-Coord_4[iDim])+
        (Coord_2[iDim]-Coord_5[iDim])+
        (Coord_3[iDim]-Coord_6[iDim]); }
      
      /*--- The normal vector should point to the interior of the element ---*/
      
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0))
        elem[iElem]->Change_Orientation();
      
    }
    
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
          || (test_4 < 0.0)) elem[iElem]->Change_Orientation();
      
    }
    
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      
      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();
      
      for(iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];
      
      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      
      if ((test_1 < 0.0) || (test_2 < 0.0))
        elem[iElem]->Change_Orientation();
      
    }
    
  }
  
}

void CPhysicalGeometry::Check_BoundElem_Orientation(CConfig *config) {
  
  unsigned long Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
  iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface;
  double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, a[3], b[3], c[3], n[3], test;
  unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
  bool find;
  
  /*--- Surface elements ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      
      iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
      for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
        Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
        find = false;
        for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
          Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
          if (Point_Surface == Point_Domain) {find = true; break;}
        }
        if (!find) break;
      }
      
      /*--- 2D grid, line case ---*/
      
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Coord_3 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
        test = a[0]*b[1]-b[0]*a[1];
        
        if (test < 0.0) {
          bound[iMarker][iElem_Surface]->Change_Orientation();
          node[Point_1_Surface]->SetFlip_Orientation();
          node[Point_2_Surface]->SetFlip_Orientation();
        }
        
      }
      
      /*--- 3D grid, triangle case ---*/
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
        Coord_4 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        
        test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        if (test < 0.0) {
          bound[iMarker][iElem_Surface]->Change_Orientation();
          node[Point_1_Surface]->SetFlip_Orientation();
          node[Point_2_Surface]->SetFlip_Orientation();
          node[Point_3_Surface]->SetFlip_Orientation();
        }
        
      }
      
      if (bound[iMarker][iElem_Surface]->GetVTK_Type() == RECTANGLE) {
        
        Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
        Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
        Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
        Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
        Coord_5 = node[Point_Domain]->GetCoord();
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_1[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
          b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_2[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
          b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_3[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        for(iDim = 0; iDim < nDim; iDim++) {
          a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
          b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]);
          c[iDim] = Coord_5[iDim]-Coord_4[iDim]; }
        n[0] = a[1]*b[2]-b[1]*a[2];
        n[1] = -(a[0]*b[2]-b[0]*a[2]);
        n[2] = a[0]*b[1]-b[0]*a[1];
        test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
        
        if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0)) {
          bound[iMarker][iElem_Surface]->Change_Orientation();
          node[Point_1_Surface]->SetFlip_Orientation();
          node[Point_2_Surface]->SetFlip_Orientation();
          node[Point_3_Surface]->SetFlip_Orientation();
          node[Point_4_Surface]->SetFlip_Orientation();
        }
        
      }
    }
}

void CPhysicalGeometry::ComputeWall_Distance(CConfig *config) {
  
  double *coord, dist2, dist;
  unsigned short iDim, iMarker;
  unsigned long iPoint, iVertex, nVertex_SolidWall;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == MASTER_NODE)
    cout << "Computing wall distances." << endl;
  
#ifndef HAVE_MPI
  
  /*--- Compute the total number of nodes on no-slip boundaries ---*/
  
  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      nVertex_SolidWall += GetnVertex(iMarker);
  
  /*--- Allocate an array to hold boundary node coordinates ---*/
  
  double **Coord_bound;
  Coord_bound = new double* [nVertex_SolidWall];
  for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
    Coord_bound[iVertex] = new double [nDim];
  
  /*--- Retrieve and store the coordinates of the no-slip boundary nodes ---*/
  
  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
       (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
       (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++)
          Coord_bound[nVertex_SolidWall][iDim] = node[iPoint]->GetCoord(iDim);
        nVertex_SolidWall++;
      }
  }
  
  /*--- Loop over all interior mesh nodes and compute the distances to each
   of the no-slip boundary nodes. Store the minimum distance to the wall for
   each interior mesh node. ---*/
  
  if (nVertex_SolidWall != 0) {
    for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
      coord = node[iPoint]->GetCoord();
      dist = 1E20;
      for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++) {
        dist2 = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dist2 += (coord[iDim]-Coord_bound[iVertex][iDim])
          *(coord[iDim]-Coord_bound[iVertex][iDim]);
        if (dist2 < dist) dist = dist2;
      }
      node[iPoint]->SetWall_Distance(sqrt(dist));
    }
  }
  else {
    for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
      node[iPoint]->SetWall_Distance(0.0);
  }
  
  /*--- Deallocate the vector of boundary coordinates. ---*/
  
  for (iVertex = 0; iVertex < nVertex_SolidWall; iVertex++)
    delete[] Coord_bound[iVertex];
  delete[] Coord_bound;
  
  
#else
  
  /*--- Variables and buffers needed for MPI ---*/
  
  int iProcessor, nProcessor;
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long nLocalVertex_NS = 0, nGlobalVertex_NS = 0, MaxLocalVertex_NS = 0;
  unsigned long *Buffer_Send_nVertex    = new unsigned long [1];
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  /*--- Count the total number of nodes on no-slip boundaries within the
   local partition. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      nLocalVertex_NS += GetnVertex(iMarker);
  
  /*--- Communicate to all processors the total number of no-slip boundary
   nodes, the maximum number of no-slip boundary nodes on any single single
   partition, and the number of no-slip nodes on each partition. ---*/
  
  Buffer_Send_nVertex[0] = nLocalVertex_NS;
  MPI_Allreduce(&nLocalVertex_NS, &nGlobalVertex_NS,  1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_NS, &MaxLocalVertex_NS, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Create and initialize to zero some buffers to hold the coordinates
   of the boundary nodes that are communicated from each partition (all-to-all). ---*/
  
  double *Buffer_Send_Coord    = new double [MaxLocalVertex_NS*nDim];
  double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_NS*nDim];
  unsigned long nBuffer = MaxLocalVertex_NS*nDim;
  
  for (iVertex = 0; iVertex < MaxLocalVertex_NS; iVertex++)
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
  
  /*--- Retrieve and store the coordinates of the no-slip boundary nodes on
   the local partition and broadcast them to all partitions. ---*/
  
  nVertex_SolidWall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC)   )
      for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[nVertex_SolidWall*nDim+iDim] = node[iPoint]->GetCoord(iDim);
        nVertex_SolidWall++;
      }
  
  MPI_Allgather(Buffer_Send_Coord, nBuffer, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer, MPI_DOUBLE, MPI_COMM_WORLD);
  
  /*--- Loop over all interior mesh nodes on the local partition and compute
   the distances to each of the no-slip boundary nodes in the entire mesh.
   Store the minimum distance to the wall for each interior mesh node. ---*/
  
  nVertex_SolidWall = 0;
  for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
    nVertex_SolidWall += Buffer_Receive_nVertex[iProcessor];
  }
  
  if (nVertex_SolidWall != 0) {
    for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
      coord = node[iPoint]->GetCoord();
      dist = 1E20;
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
        for (iVertex = 0; iVertex < Buffer_Receive_nVertex[iProcessor]; iVertex++) {
          dist2 = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dist2 += (coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS+iVertex)*nDim+iDim])*
            (coord[iDim]-Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NS+iVertex)*nDim+iDim]);
          if (dist2 < dist) dist = dist2;
        }
      node[iPoint]->SetWall_Distance(sqrt(dist));
    }
  }
  else {
    for (iPoint = 0; iPoint < GetnPoint(); iPoint++)
      node[iPoint]->SetWall_Distance(0.0);
  }
  
  /*--- Deallocate the buffers needed for the MPI communication. ---*/
  
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
  
#endif
  
}

void CPhysicalGeometry::SetPositive_ZArea(CConfig *config) {
  unsigned short iMarker, Boundary, Monitoring;
  unsigned long iVertex, iPoint;
  double *Normal, PositiveZArea;
  int rank = MASTER_NODE;
  
#ifndef HAVE_MPI
  
  PositiveZArea = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if (((Boundary == EULER_WALL)              ||
         (Boundary == HEAT_FLUX)               ||
         (Boundary == HEAT_FLUX_CATALYTIC)     ||
         (Boundary == HEAT_FLUX_NONCATALYTIC)  ||
         (Boundary == ISOTHERMAL)              ||
         (Boundary == ISOTHERMAL_CATALYTIC)    ||
         (Boundary == ISOTHERMAL_NONCATALYTIC) ||
         (Boundary == LOAD_BOUNDARY)           ||
         (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES))
      for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          Normal = vertex[iMarker][iVertex]->GetNormal();
          if (Normal[nDim-1] < 0) PositiveZArea -= Normal[nDim-1];
        }
      }
  }
  
#else
  
  double TotalPositiveZArea;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  PositiveZArea = 0.0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    
    if (((Boundary == EULER_WALL)              ||
         (Boundary == HEAT_FLUX)               ||
         (Boundary == HEAT_FLUX_CATALYTIC)     ||
         (Boundary == HEAT_FLUX_NONCATALYTIC)  ||
         (Boundary == ISOTHERMAL)              ||
         (Boundary == ISOTHERMAL_CATALYTIC)    ||
         (Boundary == ISOTHERMAL_NONCATALYTIC) ||
         (Boundary == LOAD_BOUNDARY)           ||
         (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES))
      for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          Normal = vertex[iMarker][iVertex]->GetNormal();
          if (Normal[nDim-1] < 0) PositiveZArea -= Normal[nDim-1];
        }
      }
  }
  MPI_Reduce(&PositiveZArea, &TotalPositiveZArea, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == MASTER_NODE) PositiveZArea = TotalPositiveZArea;
  MPI_Bcast(&PositiveZArea, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
  
#endif
  
  if (config->GetRefAreaCoeff() == 0.0)
    config->SetRefAreaCoeff(PositiveZArea);
  
  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Area projection in the y-plane = "<< PositiveZArea << "." << endl;
    else cout << "Area projection in the z-plane = "<< PositiveZArea << "." << endl;
  }
  
}

void CPhysicalGeometry::SetPoint_Connectivity(void) {
  
  unsigned short Node_Neighbor, iNode, iNeighbor;
  unsigned long jElem, Point_Neighbor, iPoint, iElem;
  
  /*--- Loop over all the elements ---*/
  for(iElem = 0; iElem < nElem; iElem++)
  /*--- Loop over all the nodes of an element ---*/
    for(iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      /*--- Store the element into the point ---*/
      node[iPoint]->SetElem(iElem);
    }

  /*--- Loop over all the points ---*/
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    
  /*--- Loop over all elements shared by the point ---*/
    
    for(iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
      
      jElem = node[iPoint]->GetElem(iElem);
      
      /*--- If we find the point iPoint in the surronding element ---*/
      
      for(iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)
        
        if (elem[jElem]->GetNode(iNode) == iPoint)
          
        /*--- Localize the local index of the neighbor of iPoint in the element ---*/
          
          for(iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
            Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode,iNeighbor);
            Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);
            
            /*--- Store the point into the point ---*/
            
            node[iPoint]->SetPoint(Point_Neighbor);
          }
    }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
  
}

void CPhysicalGeometry::SetRCM_Ordering(CConfig *config) {
  unsigned long iPoint, AdjPoint, AuxPoint, AddPoint, iElem, iNode, jNode;
  vector<unsigned long> Queue, AuxQueue, Result;
  unsigned short Degree, MinDegree, iDim, iMarker;
  bool *inQueue;
  
  inQueue = new bool [nPoint];
  
  for(iPoint = 0; iPoint < nPoint; iPoint++)
    inQueue[iPoint] = false;
  
  /*--- Select the node with the lowest degree in the grid. ---*/
  
  MinDegree = node[0]->GetnNeighbor(); AddPoint = 0;
  for(iPoint = 1; iPoint < nPointDomain; iPoint++) {
    Degree = node[iPoint]->GetnPoint();
    if (Degree < MinDegree) { MinDegree = Degree; AddPoint = iPoint; }
  }
  
  /*--- Add the node in the first free position. ---*/
  
  Result.push_back(AddPoint); inQueue[AddPoint] = true;
  
  /*--- Loop until reorganize all the nodes ---*/
  
  do {
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
    AuxQueue.clear();
    for (iNode = 0; iNode < node[AddPoint]->GetnPoint(); iNode++) {
      AdjPoint = node[AddPoint]->GetPoint(iNode);
      if ((!inQueue[AdjPoint]) && (AdjPoint < nPointDomain)) {
        AuxQueue.push_back(AdjPoint);
      }
    }
    
    if (AuxQueue.size() != 0) {
      
      /*--- Sort the auxiliar queue based on the number of neighbors ---*/
      
      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        for (jNode = 0; jNode < AuxQueue.size() - 1 - iNode; jNode++) {
          if (node[AuxQueue[jNode]]->GetnPoint() > node[AuxQueue[jNode+1]]->GetnPoint()) {
            AuxPoint = AuxQueue[jNode];
            AuxQueue[jNode] = AuxQueue[jNode+1];
            AuxQueue[jNode+1] = AuxPoint;
          }
        }
      }
      
      Queue.insert(Queue.end(), AuxQueue.begin(), AuxQueue.end());
      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        inQueue[AuxQueue[iNode]] = true;
      }
      
    }
    
    /*--- Extract the first node from the queue and add it in the first free
     position. ---*/
    
    if (Queue.size() != 0) {
      AddPoint = Queue[0];
      Result.push_back(Queue[0]);
      Queue.erase (Queue.begin(), Queue.begin()+1);
    }
    
    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/
    
  } while (Queue.size() != 0);
  
  /*--- Check that all the points have been added ---*/
  
  for(iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (inQueue[iPoint] == false) Result.push_back(iPoint);
  }
  
  delete[] inQueue;
  
  reverse(Result.begin(), Result.end());
  
  /*--- Add the MPI points ---*/
  
  for(iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    Result.push_back(iPoint);
  }
  
  /*--- Reset old data structures ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->ResetElem();
    node[iPoint]->ResetPoint();
    node[iPoint]->ResetBoundary();
    node[iPoint]->SetPhysicalBoundary(false);
    node[iPoint]->SetSolidBoundary(false);
    node[iPoint]->SetDomain(true);
  }
  
  /*--- Set the new coordinates ---*/
  
  double **AuxCoord;
  unsigned long *AuxGlobalIndex;
  
  AuxGlobalIndex = new unsigned long [nPoint];
  AuxCoord = new double* [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    AuxCoord[iPoint] = new double [nDim];
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    AuxGlobalIndex[iPoint] = node[iPoint]->GetGlobalIndex();
    for (iDim = 0; iDim < nDim; iDim++) {
      AuxCoord[iPoint][iDim] = node[iPoint]->GetCoord(iDim);
    }
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetGlobalIndex(AuxGlobalIndex[Result[iPoint]]);
    for (iDim = 0; iDim < nDim; iDim++)
      node[iPoint]->SetCoord(iDim, AuxCoord[Result[iPoint]][iDim]);
  }
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    delete[] AuxCoord[iPoint];
  delete[] AuxCoord;
  delete[] AuxGlobalIndex;
  
  /*--- Set the new conectivities ---*/
  
  unsigned long *InvResult;
  InvResult = new unsigned long [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    InvResult[Result[iPoint]] = iPoint;
  
  for(iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      elem[iElem]->SetNode(iNode, InvResult[iPoint]);
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      
      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      if (Marker_Tag == "SEND_RECEIVE") {
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (config->GetMarker_All_SendRecv(iMarker) < 0)
            node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
        }
      }
      
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        bound[iMarker][iElem]->SetNode(iNode, InvResult[iPoint]);
        node[InvResult[iPoint]]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[InvResult[iPoint]]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL &&
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX &&
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[InvResult[iPoint]]->SetSolidBoundary(true);
      }
    }
  }
  
  
  delete[] InvResult;
  
  /*--- MPI Synchronization point ---*/
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}

void CPhysicalGeometry::SetElement_Connectivity(void) {
  unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
  unsigned long face_point, Test_Elem, iElem;
  
  /*--- Loop over all the elements, faces and nodes ---*/
  for(iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
      for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++) {
        face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iNode));
        /*--- Loop over all elements sharing the face point ---*/
        for(jElem = 0; jElem < node[face_point]->GetnElem(); jElem++) {
          Test_Elem = node[face_point]->GetElem(jElem);
          /*--- If it is a new element in this face ---*/
          if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
              (FindFace(iElem, Test_Elem, first_elem_face, second_elem_face))) {
            /*--- Localice which faces are sharing both elements ---*/
            elem[iElem]->SetNeighbor_Elements(Test_Elem,first_elem_face);
            /*--- Store the element for both elements ---*/
            elem[Test_Elem]->SetNeighbor_Elements(iElem,second_elem_face);
            
          }
        }
      }
}

void CPhysicalGeometry::SetBoundVolume(void) {
  unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
  unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
  bool CheckVol;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      
      /*--- Choose and arbitrary point from the surface --*/
      Point = bound[iMarker][iElem_Surface]->GetNode(0);
      CheckVol = false;
      
      for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++) {
        /*--- Look for elements surronding that point --*/
        cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) cont++;
            if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
          }
          if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
        }
        
        if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) {
          bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
          CheckVol = true;
          break;
        }
      }
      if (!CheckVol) {
        cout << "The surface element ("<< iMarker <<", "<< iElem_Surface << ") doesn't have an associated volume element." << endl;
        exit(EXIT_FAILURE);
      }
    }
}

void CPhysicalGeometry::SetVertex(CConfig *config) {
  unsigned long  iPoint, iVertex, iElem;
  unsigned short iMarker, iNode;
  
  /*--- Initialize the Vertex vector for each node of the grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create and compute the vector with the number of vertex per marker ---*/
  nVertex = new unsigned long [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    nVertex[iMarker] = 0;
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)) {
          iVertex = nVertex[iMarker];
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
  
  /*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create the bound vertex structure, note that the order
   is the same as in the input file, this is important for Send/Receive part ---*/
  vertex = new CVertex**[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [nVertex[iMarker]];
    nVertex[iMarker] = 0;
    
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)){
          iVertex = nVertex[iMarker];
          vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
          
          if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
            vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
          }
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
}

void CPhysicalGeometry::SetCG(void) {
  unsigned short nNode, iDim, iMarker, iNode;
  unsigned long elem_poin, edge_poin, iElem, iEdge;
  double **Coord;
  
  /*--- Compute the center of gravity for elements ---*/
  for(iElem = 0; iElem<nElem; iElem++) {
    nNode = elem[iElem]->GetnNodes();
    Coord = new double* [nNode];
    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      elem_poin = elem[iElem]->GetNode(iNode);
      Coord[iNode] = new double [nDim];
      for (iDim = 0; iDim < nDim; iDim++)
        Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
    }
    /*--- Compute the element CG coordinates ---*/
    elem[iElem]->SetCG(Coord);
    
    for (iNode = 0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
  
  /*--- Center of gravity for face elements ---*/
  for(iMarker = 0; iMarker < nMarker; iMarker++)
    for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      nNode = bound[iMarker][iElem]->GetnNodes();
      Coord = new double* [nNode];
      /*--- Store the coordinates for all the element nodes ---*/
      for (iNode = 0; iNode < nNode; iNode++) {
        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
        Coord[iNode] = new double [nDim];
        for (iDim = 0; iDim < nDim; iDim++)
          Coord[iNode][iDim]=node[elem_poin]->GetCoord(iDim);
      }
      /*--- Compute the element CG coordinates ---*/
      bound[iMarker][iElem]->SetCG(Coord);
      for (iNode=0; iNode < nNode; iNode++)
        if (Coord[iNode] != NULL) delete[] Coord[iNode];
      if (Coord != NULL) delete[] Coord;
    }
  
  /*--- Center of gravity for edges ---*/
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    nNode = edge[iEdge]->GetnNodes();
    Coord = new double* [nNode];
    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < nNode; iNode++) {
      edge_poin=edge[iEdge]->GetNode(iNode);
      Coord[iNode] = new double [nDim];
      for (iDim = 0; iDim<nDim; iDim++)
        Coord[iNode][iDim]=node[edge_poin]->GetCoord(iDim);
    }
    /*--- Compute the edge CG coordinates ---*/
    edge[iEdge]->SetCG(Coord);
    
    for (iNode=0; iNode < nNode; iNode++)
      if (Coord[iNode] != NULL) delete[] Coord[iNode];
    if (Coord != NULL) delete[] Coord;
  }
}

void CPhysicalGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
  unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
  unsigned long Neighbor_Point, iVertex, iPoint, iElem;
  long iEdge;
  double Area, *NormalFace = NULL;
  
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  
  double *Coord_Edge_CG = new double [nDim];
  double *Coord_Elem_CG = new double [nDim];
  double *Coord_Vertex = new double [nDim];
  
  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
  /*--- Loop over all the boundary elements ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
    /*--- Loop over all the nodes of the boundary ---*/
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        iVertex = node[iPoint]->GetVertex(iMarker);
        /*--- Loop over the neighbor nodes, there is a face for each one ---*/
        for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
          Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNeighbor_Nodes);
          Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
          /*--- Shared edge by the Neighbor Point and the point ---*/
          iEdge = FindEdge(iPoint, Neighbor_Point);
          for (iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
            Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
            Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
          }
          switch (nDim) {
            case 2:
              /*--- Store the 2D face ---*/
              if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex);
              if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG);
              break;
            case 3:
              /*--- Store the 3D face ---*/
              if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex);
              if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex);
              break;
          }
        }
      }
  
  delete[] Coord_Edge_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_Vertex;
  
  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CPhysicalGeometry::MatchInterface(CConfig *config) {
  double epsilon = 1.5e-1;
  
  unsigned short nMarker_InterfaceBound = config->GetnMarker_InterfaceBound();
  
  if (nMarker_InterfaceBound != 0) {
#ifndef HAVE_MPI
    
    unsigned short iMarker, jMarker;
    unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
    double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
    
    cout << "Set Interface boundary conditions." << endl;
    
    maxdist = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          Coord_i = node[iPoint]->GetCoord();
          
          mindist = 1E6;
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
            if ((config->GetMarker_All_KindBC(jMarker) == INTERFACE_BOUNDARY) && (iMarker != jMarker))
              for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                jPoint = vertex[jMarker][jVertex]->GetNode();
                Coord_j = node[jPoint]->GetCoord();
                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                if (dist < mindist) {mindist = dist; pPoint = jPoint;}
              }
          maxdist = max(maxdist, mindist);
          vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
          
          if (mindist > epsilon) {
            cout.precision(10);
            cout << endl;
            cout << "   Bad match for point " << iPoint << ".\tNearest";
            cout << " donor distance: " << scientific << mindist << ".";
            vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
            maxdist = min(maxdist, 0.0);
          }
          
        }
        cout <<"The max distance between points is: " << maxdist <<"."<< endl;
      }
    
#else
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    unsigned short iMarker, iDim;
    unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
    double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_Interface = 0, nGlobalVertex_Interface = 0, MaxLocalVertex_Interface = 0;
    int rank, nProcessor;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    
    unsigned long *Buffer_Send_nVertex = new unsigned long [1];
    unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
    
    if (rank == MASTER_NODE) cout << "Set Interface boundary conditions (if any)." << endl;
    
    /*--- Compute the number of vertex that have interfase boundary condition
     without including the ghost nodes ---*/
    
    nLocalVertex_Interface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_Interface ++;
        }
    
    Buffer_Send_nVertex[0] = nLocalVertex_Interface;
    
    /*--- Send Interface vertex information --*/
    
    MPI_Allreduce(&nLocalVertex_Interface, &nGlobalVertex_Interface, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nLocalVertex_Interface, &MaxLocalVertex_Interface, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    
    double *Buffer_Send_Coord = new double [MaxLocalVertex_Interface*nDim];
    unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_Interface];
    
    double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_Interface*nDim];
    unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_Interface];
    
    unsigned long nBuffer_Coord = MaxLocalVertex_Interface*nDim;
    unsigned long nBuffer_Point = MaxLocalVertex_Interface;
    
    for (iVertex = 0; iVertex < MaxLocalVertex_Interface; iVertex++) {
      Buffer_Send_Point[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
    
    /*--- Copy coordinates and point to the auxiliar vector --*/
    
    nLocalVertex_Interface = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_Interface] = iPoint;
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_Coord[nLocalVertex_Interface*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_Interface++;
          }
        }
    
    MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    
    /*--- Compute the closest point to a Near-Field boundary point ---*/
    
    maxdist_local = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) {
            
            /*--- Coordinates of the boundary point ---*/
            
            Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
            
            /*--- Loop over all the boundaries to find the pair ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Interface+jVertex];
                
                /*--- Compute the distance ---*/
                
                dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                  Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Interface+jVertex)*nDim+iDim];
                  dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
                } dist = sqrt(dist);
                
                if (((dist < mindist) && (iProcessor != rank)) ||
                    ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                  mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                }
              }
            
            /*--- Store the value of the pair ---*/
            
            maxdist_local = max(maxdist_local, mindist);
            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);
            
            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPoint << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
              maxdist_local = min(maxdist_local, 0.0);
            }
            
          }
        }
      }
    }
    
    MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
    
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Point;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Point;
    
    delete[] Buffer_Send_nVertex;
    delete[] Buffer_Receive_nVertex;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
#endif
    
  }
}

void CPhysicalGeometry::MatchNearField(CConfig *config) {
  double epsilon = 1e-1;
  
  unsigned short nMarker_NearfieldBound = config->GetnMarker_NearFieldBound();
  
  if (nMarker_NearfieldBound != 0) {
    
#ifndef HAVE_MPI
    
    unsigned short iMarker, jMarker;
    unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
    double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
    
    cout << "Set Near-Field boundary conditions. " << endl;
    
    maxdist = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          Coord_i = node[iPoint]->GetCoord();
          
          mindist = 1e10;
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++)
            if ((config->GetMarker_All_KindBC(jMarker) == NEARFIELD_BOUNDARY) && (iMarker != jMarker))
              for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
                jPoint = vertex[jMarker][jVertex]->GetNode();
                Coord_j = node[jPoint]->GetCoord();
                if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
                if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
                if (dist < mindist) { mindist = dist; pPoint = jPoint; }
              }
          maxdist = max(maxdist, mindist);
          vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
          
          if (mindist > epsilon) {
            cout.precision(10);
            cout << endl;
            cout << "   Bad match for point " << iPoint << ".\tNearest";
            cout << " donor distance: " << scientific << mindist << ".";
            vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
            maxdist = min(maxdist, 0.0);
          }
        }
        cout <<"The max distance between points is: " << maxdist <<"."<< endl;
      }
    
#else
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    unsigned short iMarker, iDim;
    unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
    double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_NearField = 0, nGlobalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
    int rank, nProcessor;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
    
    unsigned long *Buffer_Send_nVertex = new unsigned long [1];
    unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
    
    if (rank == MASTER_NODE) cout << "Set Near-Field boundary conditions." << endl;
    
    /*--- Compute the number of vertex that have nearfield boundary condition
     without including the ghost nodes ---*/
    
    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_NearField ++;
        }
    
    Buffer_Send_nVertex[0] = nLocalVertex_NearField;
    
    /*--- Send Near-Field vertex information --*/
    
    MPI_Allreduce(&nLocalVertex_NearField, &nGlobalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    
    double *Buffer_Send_Coord = new double [MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_NearField];
    
    double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    
    unsigned long nBuffer_Coord = MaxLocalVertex_NearField*nDim;
    unsigned long nBuffer_Point = MaxLocalVertex_NearField;
    
    for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
      Buffer_Send_Point[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }
    
    /*--- Copy coordinates and point to the auxiliar vector --*/
    
    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_NearField] = iPoint;
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_Coord[nLocalVertex_NearField*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_NearField++;
          }
        }
    
    MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    
    /*--- Compute the closest point to a Near-Field boundary point ---*/
    
    maxdist_local = 0.0;
    maxdist_global = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) {
            
            /*--- Coordinates of the boundary point ---*/
            
            Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
            
            /*--- Loop over all the boundaries to find the pair ---*/
            
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_NearField+jVertex];
                
                /*--- Compute the distance ---*/
                
                dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                  Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NearField+jVertex)*nDim+iDim];
                  dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
                } dist = sqrt(dist);
                
                if (((dist < mindist) && (iProcessor != rank)) ||
                    ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                  mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                }
              }
            
            /*--- Store the value of the pair ---*/
            
            maxdist_local = max(maxdist_local, mindist);
            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);
            
            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPoint << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
              maxdist_local = min(maxdist_local, 0.0);
            }
            
          }
        }
        
      }
    }
    
    MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
    
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Point;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Point;
    
    delete[] Buffer_Send_nVertex;
    delete[] Buffer_Receive_nVertex;
    
    MPI_Barrier(MPI_COMM_WORLD);
    
#endif
    
  }
  
}

void CPhysicalGeometry::MatchActuator_Disk(CConfig *config) {
  double epsilon = 1e-1;
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
  double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
  int iProcessor, pProcessor = 0;
  unsigned long nLocalVertex_ActDisk = 0, nGlobalVertex_ActDisk = 0, MaxLocalVertex_ActDisk = 0;
  int rank, nProcessor;
  unsigned short Beneficiary = 0, Donor = 0, iBC;
  
  unsigned short nMarker_ActDisk_Inlet = config->GetnMarker_ActDisk_Inlet();
  
  if (nMarker_ActDisk_Inlet != 0) {
    
    for (iBC = 0; iBC < 2; iBC++) {
      
      if (iBC == 0) { Beneficiary = ACTDISK_INLET; Donor = ACTDISK_OUTLET; }
      if (iBC == 1) { Beneficiary = ACTDISK_OUTLET; Donor = ACTDISK_INLET; }
      
#ifndef HAVE_MPI
      rank = MASTER_NODE;
      nProcessor = SINGLE_NODE;
#else
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif
      
      unsigned long *Buffer_Send_nVertex = new unsigned long [1];
      unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
      
      if ((iBC == 0) && (rank == MASTER_NODE)) cout << "Set Actuator Disk inlet boundary conditions." << endl;
      if ((iBC == 1) && (rank == MASTER_NODE)) cout << "Set Actuator Disk outlet boundary conditions." << endl;
      
      /*--- Compute the number of vertex that have an actuator disk outlet boundary condition
       without including the ghost nodes ---*/
      
      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            if (node[iPoint]->GetDomain()) nLocalVertex_ActDisk ++;
          }
        }
      }
      
      Buffer_Send_nVertex[0] = nLocalVertex_ActDisk;
      
      /*--- Send actuator disk vertex information --*/
      
#ifndef HAVE_MPI
      nGlobalVertex_ActDisk = nLocalVertex_ActDisk;
      MaxLocalVertex_ActDisk = nLocalVertex_ActDisk;
      Buffer_Receive_nVertex[0] = Buffer_Send_nVertex[0];
#else
      MPI_Allreduce(&nLocalVertex_ActDisk, &nGlobalVertex_ActDisk, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&nLocalVertex_ActDisk, &MaxLocalVertex_ActDisk, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
      
      /*--- Array dimensionalization --*/
      
      double *Buffer_Send_Coord = new double [MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Send_Point  = new unsigned long [MaxLocalVertex_ActDisk];
      
      double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      
      unsigned long nBuffer_Coord = MaxLocalVertex_ActDisk*nDim;
      unsigned long nBuffer_Point = MaxLocalVertex_ActDisk;
      
      for (iVertex = 0; iVertex < MaxLocalVertex_ActDisk; iVertex++) {
        Buffer_Send_Point[iVertex] = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      }
      
      /*--- Copy coordinates and point to the auxiliar vector --*/
      
      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            if (node[iPoint]->GetDomain()) {
              Buffer_Send_Point[nLocalVertex_ActDisk] = iPoint;
              for (iDim = 0; iDim < nDim; iDim++)
                Buffer_Send_Coord[nLocalVertex_ActDisk*nDim+iDim] = node[iPoint]->GetCoord(iDim);
              nLocalVertex_ActDisk++;
            }
          }
        }
      }
      
#ifndef HAVE_MPI
      for (unsigned long iBuffer_Coord = 0; iBuffer_Coord < nBuffer_Coord; iBuffer_Coord++)
        Buffer_Receive_Coord[iBuffer_Coord] = Buffer_Send_Coord[iBuffer_Coord];
      for (unsigned long iBuffer_Point = 0; iBuffer_Point < nBuffer_Point; iBuffer_Point++)
        Buffer_Receive_Point[iBuffer_Point] = Buffer_Send_Point[iBuffer_Point];
#else
      MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#endif
      
      /*--- Compute the closest point to an actuator disk inlet point ---*/
      
      maxdist_local = 0.0; maxdist_global = 0.0;
      
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Beneficiary) {
          
          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            if (node[iPoint]->GetDomain()) {
              
              /*--- Coordinates of the boundary point ---*/
              
              Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
              
              /*--- Loop over all the boundaries to find the pair ---*/
              
              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                  jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  
                  /*--- Compute the distance ---*/
                  
                  dist = 0.0;
                  for (iDim = 0; iDim < nDim; iDim++) {
                    Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_ActDisk+jVertex)*nDim+iDim];
                    dist += pow(Coord_j[iDim]-Coord_i[iDim], 2.0);
                  }
                  dist = sqrt(dist);
                  
                  if (dist < mindist) {
                    mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
                    if (dist == 0.0) break;
                  }
                  
                }
              }
              
              /*--- Store the value of the pair ---*/
              
              maxdist_local = max(maxdist_local, mindist);
              vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);
              
              if (mindist > epsilon) {
                cout.precision(10);
                cout << endl;
                cout << "   Bad match for point " << iPoint << ".\tNearest";
                cout << " donor distance: " << scientific << mindist << ".";
                vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
                maxdist_local = min(maxdist_local, 0.0);
              }
              
            }
          }
          
        }
      }
      
#ifndef HAVE_MPI
      maxdist_global = maxdist_local;
#else
      MPI_Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
#endif
      
      if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;
      
      delete[] Buffer_Send_Coord;
      delete[] Buffer_Send_Point;
      
      delete[] Buffer_Receive_Coord;
      delete[] Buffer_Receive_Point;
      
      delete[] Buffer_Send_nVertex;
      delete[] Buffer_Receive_nVertex;
      
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      
    }
  }
  
}

void CPhysicalGeometry::MatchZone(CConfig *config, CGeometry *geometry_donor, CConfig *config_donor,
                                  unsigned short val_iZone, unsigned short val_nZone) {
  
#ifndef HAVE_MPI
  
  unsigned short iMarker, jMarker;
  unsigned long iVertex, iPoint, jVertex, jPoint = 0, pPoint = 0;
  double *Coord_i, *Coord_j, dist = 0.0, mindist, maxdist;
  
  if (val_iZone == ZONE_0) cout << "Set zone boundary conditions (if any)." << endl;
  
  maxdist = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iPoint = vertex[iMarker][iVertex]->GetNode();
      Coord_i = node[iPoint]->GetCoord();
      
      mindist = 1E6;
      for (jMarker = 0; jMarker < config_donor->GetnMarker_All(); jMarker++)
        for (jVertex = 0; jVertex < geometry_donor->GetnVertex(jMarker); jVertex++) {
          jPoint = geometry_donor->vertex[jMarker][jVertex]->GetNode();
          Coord_j = geometry_donor->node[jPoint]->GetCoord();
          if (nDim == 2) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0));
          if (nDim == 3) dist = sqrt(pow(Coord_j[0]-Coord_i[0],2.0) + pow(Coord_j[1]-Coord_i[1],2.0) + pow(Coord_j[2]-Coord_i[2],2.0));
          if (dist < mindist) { mindist = dist; pPoint = jPoint; }
        }
      
      maxdist = max(maxdist, mindist);
      vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
      
    }
  }
  
#else
  MPI_Barrier(MPI_COMM_WORLD);
  
  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint, pPoint = 0, jVertex, jPoint;
  double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist;
  int iProcessor, pProcessor = 0;
  unsigned long nLocalVertex_Zone = 0, nGlobalVertex_Zone = 0, MaxLocalVertex_Zone = 0;
  int rank, nProcessor;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];
  
  if (val_iZone == ZONE_0) cout << "Set zone boundary conditions (if any)." << endl;
  
  nLocalVertex_Zone = 0;
  for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
      if (geometry_donor->node[iPoint]->GetDomain()) nLocalVertex_Zone ++;
    }
  
  Buffer_Send_nVertex[0] = nLocalVertex_Zone;
  
  /*--- Send Interface vertex information --*/
  
  MPI_Allreduce(&nLocalVertex_Zone, &nGlobalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nLocalVertex_Zone, &MaxLocalVertex_Zone, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  double *Buffer_Send_Coord = new double [MaxLocalVertex_Zone*nDim];
  unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_Zone];
  
  double *Buffer_Receive_Coord = new double [nProcessor*MaxLocalVertex_Zone*nDim];
  unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_Zone];
  
  unsigned long nBuffer_Coord = MaxLocalVertex_Zone*nDim;
  unsigned long nBuffer_Point = MaxLocalVertex_Zone;
  
  for (iVertex = 0; iVertex < MaxLocalVertex_Zone; iVertex++) {
    Buffer_Send_Point[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++)
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
  }
  
  /*--- Copy coordinates and point to the auxiliar vector --*/
  nLocalVertex_Zone = 0;
  for (iMarker = 0; iMarker < config_donor->GetnMarker_All(); iMarker++)
    for (iVertex = 0; iVertex < geometry_donor->GetnVertex(iMarker); iVertex++) {
      iPoint = geometry_donor->vertex[iMarker][iVertex]->GetNode();
      if (geometry_donor->node[iPoint]->GetDomain()) {
        Buffer_Send_Point[nLocalVertex_Zone] = iPoint;
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[nLocalVertex_Zone*nDim+iDim] = geometry_donor->node[iPoint]->GetCoord(iDim);
        nLocalVertex_Zone++;
      }
    }
  
  MPI_Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  
  /*--- Compute the closest point to a Near-Field boundary point ---*/
  maxdist = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iPoint = vertex[iMarker][iVertex]->GetNode();
      
      if (node[iPoint]->GetDomain()) {
        
        /*--- Coordinates of the boundary point ---*/
        Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;
        
        /*--- Loop over all the boundaries to find the pair ---*/
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
          for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
            jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_Zone+jVertex];
            
            /*--- Compute the distance ---*/
            dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
              Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Zone+jVertex)*nDim+iDim];
              dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
            } dist = sqrt(dist);
            
            if (((dist < mindist) && (iProcessor != rank)) ||
                ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
              mindist = dist; pProcessor = iProcessor; pPoint = jPoint;
            }
          }
        
        /*--- Store the value of the pair ---*/
        maxdist = max(maxdist, mindist);
        vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pProcessor);
        
        
      }
    }
  }
  
  delete[] Buffer_Send_Coord;
  delete[] Buffer_Send_Point;
  
  delete[] Buffer_Receive_Coord;
  delete[] Buffer_Receive_Point;
  
  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Receive_nVertex;
  
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}


void CPhysicalGeometry::SetControlVolume(CConfig *config, unsigned short action) {
  unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
  Volume, DomainVolume, my_DomainVolume, *NormalFace = NULL;
  bool change_face_orientation;
  int rank;
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE) {
    for(iEdge = 0; iEdge < nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
    for(iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetVolume (0.0);
  }
  
  Coord_Edge_CG = new double [nDim];
  Coord_FaceElem_CG = new double [nDim];
  Coord_Elem_CG = new double [nDim];
  Coord_FaceiPoint = new double [nDim];
  Coord_FacejPoint = new double [nDim];
  
  my_DomainVolume = 0.0;
  for(iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
      
      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
      
      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
        
        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }
        
        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }
        
        /*--- We define a direction (from the smalest index to the greatest) --*/
        change_face_orientation = false;
        if (face_iPoint > face_jPoint) change_face_orientation = true;
        iEdge = FindEdge(face_iPoint, face_jPoint);
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace,iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }
        
        switch (nDim) {
          case 2:
            /*--- Two dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_Elem_CG);
            Area = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Area); my_DomainVolume +=Area;
            Area = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Area); my_DomainVolume +=Area;
            break;
          case 3:
            /*--- Three dimensional problem ---*/
            if (change_face_orientation) edge[iEdge]->SetNodes_Coord(Coord_FaceElem_CG,Coord_Edge_CG,Coord_Elem_CG);
            else edge[iEdge]->SetNodes_Coord(Coord_Edge_CG,Coord_FaceElem_CG,Coord_Elem_CG);
            Volume = edge[iEdge]->GetVolume(Coord_FaceiPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            Volume = edge[iEdge]->GetVolume(Coord_FacejPoint,Coord_Edge_CG,Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            break;
        }
      }
    }
  
  /*--- Check if there is a normal with null area ---*/
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
  
#ifdef HAVE_MPI
  MPI_Allreduce(&my_DomainVolume, &DomainVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  DomainVolume = my_DomainVolume;
#endif
  
  if ((rank == MASTER_NODE) && (action == ALLOCATE)) {
    if (nDim == 2) cout <<"Area of the computational grid: "<< DomainVolume <<"."<< endl;
    if (nDim == 3) cout <<"Volume of the computational grid: "<< DomainVolume <<"."<< endl;
  }
  
  config->SetDomainVolume(DomainVolume);
  
  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
}

void CPhysicalGeometry::VisualizeControlVolume(CConfig *config, unsigned short action) {
  
  /*--- This routine is only meant for visualization in serial currently ---*/
#ifndef HAVE_MPI
  
  unsigned long face_iPoint = 0, face_jPoint = 0, iElem, iPoint_Viz;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint,
  *Coord_FacejPoint;
  bool change_face_orientation;
  int rank;
  int counter = 0;
  char cstr[MAX_STRING_SIZE], buffer[50];
  ofstream Tecplot_File;
  string mesh_filename;
  vector<double> X, Y, Z, X_n, Y_n, Z_n;
  double r1[3], r2[3], CrossProduct[3];
  double Area_Dual = 0.0;
  
  /*--- Only serial visualization ---*/
  
  rank = MASTER_NODE;
  
  /*--- Access the point number for control volume we want to vizualize ---*/
  
  iPoint_Viz = config->GetVisualize_CV();
  
  /*--- Allocate some structures for building the dual CVs ---*/
  
  Coord_Edge_CG     = new double [nDim];
  Coord_FaceElem_CG = new double [nDim];
  Coord_Elem_CG     = new double [nDim];
  Coord_FaceiPoint  = new double [nDim];
  Coord_FacejPoint  = new double [nDim];
  
  /*--- Loop over each face of each element ---*/
  
  CrossProduct[0] = 0.0; CrossProduct[1] = 0.0; CrossProduct[2] = 0.0;
  
  for(iElem = 0; iElem < nElem; iElem++) {
    
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {
      
      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);
      
      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {
        
        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }
        
        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }
        
        /*--- We define a direction (from the smallest index to the greatest) --*/
        change_face_orientation = false;
        if (face_iPoint > face_jPoint) change_face_orientation = true;
        iEdge = FindEdge(face_iPoint, face_jPoint);
        
        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edge[iEdge]->GetCG(iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace,iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }
        
        /*--- Print out the coordinates for a set of triangles making
         up a single dual control volume for visualization. ---*/

        if (face_iPoint == iPoint_Viz || face_jPoint == iPoint_Viz) {

          if (nDim == 2) {
            X.push_back(Coord_Elem_CG[0]); X.push_back(Coord_Edge_CG[0]);
            Y.push_back(Coord_Elem_CG[1]); Y.push_back(Coord_Edge_CG[1]);
          } else if (nDim == 3) {
            X.push_back(Coord_FaceElem_CG[0]); X.push_back(Coord_Edge_CG[0]); X.push_back(Coord_Elem_CG[0]);
            Y.push_back(Coord_FaceElem_CG[1]); Y.push_back(Coord_Edge_CG[1]); Y.push_back(Coord_Elem_CG[1]);
            Z.push_back(Coord_FaceElem_CG[2]); Z.push_back(Coord_Edge_CG[2]); Z.push_back(Coord_Elem_CG[2]);

            for (iDim = 0; iDim < nDim; iDim++) {
              r1[iDim] = Coord_FaceElem_CG[iDim]-Coord_Elem_CG[iDim];
              r2[iDim] = Coord_Edge_CG[iDim]-Coord_Elem_CG[iDim];
            }
            CrossProduct[0] += 0.5*(r1[1]*r2[2] - r1[2]*r2[1]);
            CrossProduct[1] += 0.5*(r1[2]*r2[0] - r1[0]*r2[2]);
            CrossProduct[2] += 0.5*(r1[0]*r2[1] - r1[1]*r2[0]);
          }
          counter++;
        }
      }
    }
  }
  
  Area_Dual = sqrt( CrossProduct[0]*CrossProduct[0]
                   +CrossProduct[1]*CrossProduct[1]
                   +CrossProduct[2]*CrossProduct[2]);

  /*--- Write a Tecplot file to visualize the CV ---*/
  
  strcpy(cstr,"dual_cv");
  sprintf (buffer, "_%d.dat", int(iPoint_Viz));
  strcat(cstr,buffer);
  
  Tecplot_File.open(cstr, ios::out);
  Tecplot_File << "TITLE= \"Visualization of the control volume\"" << endl;
  
  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*2 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  } if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*3 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }
  
  /*--- Write coordinates for the nodes in the order that they were found
   for each of the edges/triangles making up a dual control volume. ---*/
  
  for(vector<double>::size_type i = 0; i != X.size(); i++) {
    Tecplot_File << X[i] << "\t" << Y[i];
    if (nDim == 3) Tecplot_File << "\t" << Z[i];
    Tecplot_File << "\n";
  }
  
  /*--- Create a new connectivity table in the order the faces were found ---*/
  
  int j;
  for (int i= 0; i < counter; i++){
    if (nDim == 2) {
      j = i*2;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+2 <<"\t"<<j+2 << endl;
    } if (nDim == 3) {
      j = i*3;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+3 <<"\t"<<j+3 <<"\t";
      Tecplot_File << j+3<<"\t" <<j+3 <<"\t"<<j+3 <<"\t"<<j+3 << endl;
    }
  }
  
  Tecplot_File.close();
  X.clear();
  Y.clear();
  Z.clear();
  
  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
  
#endif

}

void CPhysicalGeometry::SetMeshFile (CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound;
  unsigned short iMarker, iNodes, iDim;
  unsigned short iPeriodic, nPeriodic = 0;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;
  double *center, *angles, *transl;
  
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  
  /*--- Open .su2 grid file ---*/
  
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Write dimension, number of elements and number of points ---*/
  
  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << elem[iElem]->GetNode(iNodes);
    output_file << "\t"<<iElem<< endl;
  }
  
  /*--- Write the node coordinates ---*/
  
  output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
  output_file.precision(15);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
#ifndef HAVE_MPI
    output_file << "\t" << iPoint << endl;
#else
    output_file << "\t" << iPoint << "\t" << node[iPoint]->GetGlobalIndex() << endl;
#endif
    
  }
  
  /*--- Loop through and write the boundary info ---*/
  
  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Ignore SEND_RECEIVE for the moment ---*/
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      
      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file	<< iElem_Bound << endl;
        }
      }
      
    } else if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type() << endl;
      }
      
    }
  }
  
  /*--- Get the total number of periodic transformations ---*/
  
  nPeriodic = config->GetnPeriodicIndex();
  output_file << "NPERIODIC= " << nPeriodic << endl;
  
  /*--- From iPeriodic obtain the iMarker ---*/
  
  for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
    
    /*--- Retrieve the supplied periodic information. ---*/
    
    center = config->GetPeriodicCenter(iPeriodic);
    angles = config->GetPeriodicRotation(iPeriodic);
    transl = config->GetPeriodicTranslate(iPeriodic);
    
    output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
    output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
    output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
    output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
    
  }
  
  
  output_file.close();
}

void CPhysicalGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, double val_smooth_coeff, CConfig *config) {
  unsigned short iSmooth, nneigh, iMarker;
  double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  double eps = 1E-6;
  bool NearField = false;
  
  Coord = new double [nDim];
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    double *Coord = node[iPoint]->GetCoord();
    node[iPoint]->SetCoord_Old(Coord);
  }
  
  /*--- Jacobi iterations ---*/
  for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetCoord_SumZero();
    
    
    /*--- Loop over Interior edges ---*/
    for(iEdge = 0; iEdge < nEdge; iEdge++) {
      iPoint = edge[iEdge]->GetNode(0);
      Coord_i = node[iPoint]->GetCoord();
      
      jPoint = edge[iEdge]->GetNode(1);
      Coord_j = node[jPoint]->GetCoord();
      
      /*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
      node[iPoint]->AddCoord_Sum(Coord_j);
      node[jPoint]->AddCoord_Sum(Coord_i);
      
    }
    
    /*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      nneigh = node[iPoint]->GetnPoint();
      Coord_Sum = node[iPoint]->GetCoord_Sum();
      Coord_Old = node[iPoint]->GetCoord_Old();
      
      if (nDim == 2) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
        if ((NearField) && ((Coord_Old[1] > Position_Plane-eps) && (Coord_Old[1] < Position_Plane+eps)))
          Coord[1] = Coord_Old[1];
      }
      
      if (nDim == 3) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*double(nneigh));
        Coord[2] =(Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) /(1.0 + val_smooth_coeff*double(nneigh));
        if ((NearField) && ((Coord_Old[2] > Position_Plane-eps) && (Coord_Old[2] < Position_Plane+eps)))
          Coord[2] = Coord_Old[2];
      }
      
      node[iPoint]->SetCoord(Coord);
    }
    
    /*--- Copy boundary values ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Coord_Old = node[iPoint]->GetCoord_Old();
        node[iPoint]->SetCoord(Coord_Old);
      }
  }
  
  delete[] Coord;
}

bool CPhysicalGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                                 unsigned short &face_second_elem) {
  
  /*--- Find repeated nodes between two elements to identify the common face ---*/
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short face_node, iFace, iNode, jNode, kNode, nNodesFace;
  vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
  vector<unsigned long>::iterator IterPoint;
  pair<vector <unsigned long>::iterator, vector <unsigned long>::iterator> mypair;
  bool face_first_found = false, face_second_found =false;
  
  if (first_elem == second_elem) return 0;
  
  kNode = 0;
  for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) {
    iPoint = elem[first_elem]->GetNode(iNode);
    for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) {
      jPoint = elem[second_elem]->GetNode(jNode);
      if (iPoint == jPoint) {
        CommonPoints.push_back(iPoint);
        break;
      }
    }
  }
  
  /*--- Sort point in face and check that the list is unique ---*/
  sort( CommonPoints.begin(), CommonPoints.end());
  IterPoint = unique( CommonPoints.begin(), CommonPoints.end());
  CommonPoints.resize( distance(CommonPoints.begin(), IterPoint) );
  
  /*--- Search the secuence in the first element ---*/
  for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[first_elem]->GetnNodesFace(iFace);
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[first_elem]->GetFaces(iFace, iNode);
      PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceFirst.begin(), PointFaceFirst.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
    if (mypair.first == PointFaceFirst.end()) {
      face_first_elem = iFace;
      face_first_found = true;
      break;
    }
    
    PointFaceFirst.erase (PointFaceFirst.begin(),PointFaceFirst.end());
  }
  
  /*--- Search the secuence in the second element ---*/
  for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[second_elem]->GetnNodesFace(iFace);
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[second_elem]->GetFaces(iFace,iNode);
      PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
    }
    
    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceSecond.begin(), PointFaceSecond.end());
    
    /*--- List comparison ---*/
    mypair = mismatch (PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
    if (mypair.first == PointFaceSecond.end()) {
      face_second_elem = iFace;
      face_second_found = true;
      break;
    }
    
    PointFaceSecond.erase (PointFaceSecond.begin(),PointFaceSecond.end());
  }
  
  if (face_first_found && face_second_found) return true;
  else return false;
  
}

void CPhysicalGeometry::SetTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file) {
  
  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;
  
  /*--- Open the tecplot file and write the header ---*/
  
  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);
  
  Tecplot_File << "ZONE T= ";
  if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
  else Tecplot_File << "\"Deformed grid\", C=RED, ";
  Tecplot_File << "NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING= POINT";
  if (nDim == 2) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
  if (nDim == 3) Tecplot_File << ", ZONETYPE= FEBRICK"<< endl;
  
  /*--- Adding coordinates ---*/
  
  for(iPoint = 0; iPoint < nPoint; iPoint++) {
    for(iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }
  
  /*--- Adding conectivity ---*/
  
  for(iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == WEDGE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }
  
  Tecplot_File.close();
}

void CPhysicalGeometry::SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {
  
  ofstream Tecplot_File;
  unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
  unsigned short Coord_i, iMarker;
  
  /*--- It is important to do a renumbering to don't add points
   that do not belong to the surfaces ---*/
  
  PointSurface = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (node[iPoint]->GetBoundary()) {
      PointSurface[iPoint] = nPointSurface;
      nPointSurface++;
    }
  
  /*--- Compute the total number of elements ---*/
  
  Total_nElem_Bound = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      Total_nElem_Bound += nElem_Bound[iMarker];
    }
  }
  
  /*--- Open the tecplot file and write the header ---*/
  
  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the surface grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);
  
  if (Total_nElem_Bound != 0) {
    
    /*--- Write the header of the file ---*/
    
    Tecplot_File << "ZONE T= ";
    if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
    else Tecplot_File << "\"Deformed grid\", C=RED, ";
    Tecplot_File << "NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING= POINT";
    if (nDim == 2) Tecplot_File << ", ZONETYPE= FELINESEG"<< endl;
    if (nDim == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
    
    /*--- Only write the coordiantes of the points that are on the surfaces ---*/
    
    if (nDim == 3) {
      for(iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for(Coord_i = 0; Coord_i < nDim-1; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
        }
    }
    else {
      for(iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()){
          for(Coord_i = 0; Coord_i < nDim; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << "\n";
        }
    }
    
    /*--- Write the cells using the new numbering ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
          if (nDim == 2) {
            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
          }
          if (nDim == 3) {
            if (bound[iMarker][iElem]->GetnNodes() == 3) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
            }
            if (bound[iMarker][iElem]->GetnNodes() == 4) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
            }
          }
        }
  }
  else {
    
    /*--- No elements in the surface ---*/
    
    if (nDim == 2) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      Tecplot_File << "0.0 0.0"<< endl;
      Tecplot_File << "1 1"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      Tecplot_File << "0.0 0.0 0.0"<< endl;
      Tecplot_File << "1 1 1 1"<< endl;
    }
  }
  
  /*--- Dealocate memory and close the file ---*/
  
  delete[] PointSurface;
  Tecplot_File.close();
  
}

void CPhysicalGeometry::SetBoundSTL(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {
  
  ofstream STL_File;
  unsigned long this_node, iNode, nNode, iElem;
  unsigned short iDim, iMarker;
  double p[3], u[3], v[3], n[3], a;
  
  /*---	STL format:
   solid NAME
   ...
   facet normal 0.00 0.00 1.00
   outer loop
   vertex  2.00  2.00  0.00
   vertex -1.00  1.00  0.00
   vertex  0.00 -1.00  0.00
   endloop
   endfacet
   ...
   end solid
   --- */
  
  /*--- Open the STL file ---*/
  
  if (new_file) STL_File.open(mesh_filename, ios::out);
  else STL_File.open(mesh_filename, ios::out | ios::app);
  
  /*--- Write the header of the file ---*/
  
  STL_File << "solid surface_mesh" << endl;
  
  /*--- Write facets of surface markers ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_Plotting(iMarker) == YES)
      for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
        
        /*--- number of nodes for this elemnt ---*/
        
        nNode = bound[iMarker][iElem]->GetnNodes();
        
        /*--- Calculate Normal Vector ---*/
        
        for (iDim=0; iDim<nDim; iDim++){
          p[0] = node[bound[iMarker][iElem]->GetNode(0)]      ->GetCoord(iDim);
          p[1] = node[bound[iMarker][iElem]->GetNode(1)]      ->GetCoord(iDim);
          p[2] = node[bound[iMarker][iElem]->GetNode(nNode-1)]->GetCoord(iDim);
          u[iDim] = p[1]-p[0];
          v[iDim] = p[2]-p[0];
        }
        
        n[0] = u[1]*v[2]-u[2]*v[1];
        n[1] = u[2]*v[0]-u[0]*v[2];
        n[2] = u[0]*v[1]-u[1]*v[0];
        a = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        
        /*--- Print normal vector ---*/
        
        STL_File << "  facet normal ";
        for (iDim=0; iDim<nDim; iDim++){
          STL_File << n[iDim]/a << " ";
        }
        STL_File << endl;
        
        /*--- STL Facet Loop --*/
        
        STL_File << "    outer loop" << endl;
        
        /*--- Print Nodes for Facet ---*/
        
        for(iNode=0; iNode<nNode; iNode++) {
          this_node = bound[iMarker][iElem]->GetNode(iNode);
          STL_File << "      vertex ";
          for (iDim = 0; iDim < nDim; iDim++)
            STL_File << node[this_node]->GetCoord(iDim) << " ";
          if (nDim==2)
            STL_File << 0.0 << " ";
          STL_File <<  endl;
        }
        STL_File << "    endloop" << endl;
        STL_File << "  endfacet" << endl;
      }
  
  /*--- Done with Surface Mesh ---*/
  
  STL_File << "endsolid" << endl;
  
  /*--- Close the file ---*/
  
  STL_File.close();
  
}

void CPhysicalGeometry::SetColorGrid(CConfig *config) {
  
#ifdef HAVE_MPI
  
#ifdef HAVE_METIS
  
  unsigned long iPoint, iElem, iElem_Triangle, iElem_Tetrahedron, nElem_Triangle,
  nElem_Tetrahedron, kPoint, jPoint, iVertex;
  unsigned short iMarker, iMaxColor = 0, iColor, MaxColor = 0, iNode, jNode;
  idx_t ne = 0, nn, *elmnts = NULL, etype, *epart = NULL, *npart = NULL, numflag, nparts, edgecut, *eptr;
  int rank, size;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  unsigned short nDomain = size;
  
  nElem_Triangle = 0;
  nElem_Tetrahedron = 0;
  for (iElem = 0; iElem < GetnElem(); iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) nElem_Triangle = nElem_Triangle + 1;
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) nElem_Triangle = nElem_Triangle + 2;
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) nElem_Tetrahedron = nElem_Tetrahedron + 1;
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) nElem_Tetrahedron = nElem_Tetrahedron + 5;
    if (elem[iElem]->GetVTK_Type() == PYRAMID) nElem_Tetrahedron = nElem_Tetrahedron + 2;
    if (elem[iElem]->GetVTK_Type() == WEDGE) nElem_Tetrahedron = nElem_Tetrahedron + 3;
  }
  
  if (GetnDim() == 2) {
    ne = nElem_Triangle;
    elmnts = new idx_t [ne*3];
    etype = 1;
  }
  if (GetnDim() == 3) {
    ne = nElem_Tetrahedron;
    elmnts = new idx_t [ne*4];
    etype = 2;
  }
  
  nn = nPoint;
  numflag = 0;
  nparts = nDomain;
  epart = new idx_t [ne];
  npart = new idx_t [nn];
  eptr  = new idx_t[ne+1];
  
  /*--- Initialize the color vector ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetColor(0);
  
  if (nparts > 1) {
    
    iElem_Triangle = 0; iElem_Tetrahedron = 0;
    for (iElem = 0; iElem < GetnElem(); iElem++) {
      if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(1);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(2);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
      }
      if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(1);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(2);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
        elmnts[3*iElem_Triangle+0]= elem[iElem]->GetNode(0);
        elmnts[3*iElem_Triangle+1]= elem[iElem]->GetNode(2);
        elmnts[3*iElem_Triangle+2]= elem[iElem]->GetNode(3);
        eptr[iElem_Triangle] = 3*iElem_Triangle;
        iElem_Triangle++;
      }
      if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(3);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(5);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(7);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(6);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(7);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(5);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == PYRAMID) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
      if (elem[iElem]->GetVTK_Type() == WEDGE) {
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(1);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(4);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(2);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(0);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(2);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(4);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
        elmnts[4*iElem_Tetrahedron+0]= elem[iElem]->GetNode(3);
        elmnts[4*iElem_Tetrahedron+1]= elem[iElem]->GetNode(4);
        elmnts[4*iElem_Tetrahedron+2]= elem[iElem]->GetNode(5);
        elmnts[4*iElem_Tetrahedron+3]= elem[iElem]->GetNode(2);
        eptr[iElem_Tetrahedron] = 4*iElem_Tetrahedron;
        iElem_Tetrahedron++;
      }
    }
    
    /*--- Add final value to element pointer array ---*/
    
    if (GetnDim() == 2) eptr[ne] = 3*ne;
    else eptr[ne] = 4*ne;
    
    METIS_PartMeshNodal(&ne, &nn, eptr, elmnts, NULL, NULL, &nparts, NULL, NULL, &edgecut, epart, npart);
    
    cout << "Finished partitioning using METIS. ("  << edgecut << " edge cuts)." << endl;
    
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetColor(npart[iPoint]);
  }
  
  delete[] epart;
  delete[] npart;
  delete[] elmnts;
  delete[] eptr;
  
#endif
#endif
  
}

void CPhysicalGeometry::GetQualityStatistics(double *statistics) {
  unsigned long jPoint, Point_2, Point_3, iElem;
  double *Coord_j, *Coord_2, *Coord_3;
  unsigned short iDim;
  
  statistics[0] = 1e06;
  statistics[1] = 0;
  
  /*--- Loop interior edges ---*/
  for (iElem = 0; iElem < this->GetnElem(); iElem++) {
    
    if ((this->GetnDim() == 2) && (elem[iElem]->GetVTK_Type() == TRIANGLE)) {
      
      jPoint = elem[iElem]->GetNode(0); Coord_j = node[jPoint]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      
      /*--- Compute sides of the triangle ---*/
      double a = 0, b = 0, c = 0;
      for (iDim = 0; iDim < nDim; iDim++) {
        a += (Coord_2[iDim]-Coord_j[iDim])*(Coord_2[iDim]-Coord_j[iDim]);
        b += (Coord_3[iDim]-Coord_j[iDim])*(Coord_3[iDim]-Coord_j[iDim]);
        c += (Coord_3[iDim]-Coord_2[iDim])*(Coord_3[iDim]-Coord_2[iDim]);
      }
      a = sqrt(a); b = sqrt(b); c = sqrt(c);
      
      /*--- Compute semiperimeter (s) and area ---*/
      double s = 0.5*(a + b + c);
      double Area = sqrt(s*(s-a)*(s-b)*(s-c));
      
      /*--- Compute radius of the circumcircle (R) and of the incircle (r) ---*/
      double R = (a*b*c) / (4.0*Area);
      double r = Area / s;
      double roR = r / R;
      
      /*--- Update statistics ---*/
      if (roR < statistics[0])
        statistics[0] = roR;
      statistics[1] += roR;
      
    }
  }
  statistics[1] /= this->GetnElem();
  
}

void CPhysicalGeometry::SetRotationalVelocity(CConfig *config) {
  
  unsigned long iPoint;
  double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Center of rotation & angular velocity vector from config ---*/
  
  Center[0] = config->GetMotion_Origin_X(ZONE_0);
  Center[1] = config->GetMotion_Origin_Y(ZONE_0);
  Center[2] = config->GetMotion_Origin_Z(ZONE_0);
  Omega[0]  = config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref();
  Omega[1]  = config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref();
  Omega[2]  = config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref();
  L_Ref     = config->GetLength_Ref();
  
  /*--- Print some information to the console ---*/
  
  if (rank == MASTER_NODE) {
    cout << " Rotational origin (x,y,z): ( " << Center[0] << ", " << Center[1];
    cout << ", " << Center[2] << " )" << endl;
    cout << " Angular velocity about x, y, z axes: ( " << Omega[0] << ", ";
    cout << Omega[1] << ", " << Omega[2] << " ) rad/s" << endl;
  }
  
  /*--- Loop over all nodes and set the rotational velocity ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Get the coordinates of the current node ---*/
    
    Coord = node[iPoint]->GetCoord();
    
    /*--- Calculate the non-dim. distance from the rotation center ---*/
    
    Distance[0] = (Coord[0]-Center[0])/L_Ref;
    Distance[1] = (Coord[1]-Center[1])/L_Ref;
    Distance[2] = (Coord[2]-Center[2])/L_Ref;
    
    /*--- Calculate the angular velocity as omega X r ---*/
    
    RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
    RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
    RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
    
    /*--- Store the grid velocity at this node ---*/
    
    node[iPoint]->SetGridVel(RotVel);
    
  }
  
}

void CPhysicalGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
  
  /*--- Local variables ---*/
  
  double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
  double TimeStep, GridVel = 0.0;
  unsigned long iPoint;
  unsigned short iDim;
  
  /*--- Compute the velocity of each node in the volume mesh ---*/
  
  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    
    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
    
    Coord_nM1 = node[iPoint]->GetCoord_n1();
    Coord_n   = node[iPoint]->GetCoord_n();
    Coord_nP1 = node[iPoint]->GetCoord();
    
    /*--- Unsteady time step ---*/
    
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
    
    for(iDim = 0; iDim < nDim; iDim++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                   + 1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
      
      /*--- Store grid velocity for this point ---*/
      
      node[iPoint]->SetGridVel(iDim, GridVel);
    }
  }
  
}

void CPhysicalGeometry::Set_MPI_Coord(CConfig *config)  {
  
  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_Coord = NULL, *Buffer_Send_Coord = NULL, *Coord = NULL, *newCoord = NULL;
  int send_to, receive_from;
  
  newCoord = new double[nDim];
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
      nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_Coord = new double [nBufferR_Vector];
      Buffer_Send_Coord = new double[nBufferS_Vector];
      
      /*--- Copy the coordinates that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = vertex[MarkerS][iVertex]->GetNode();
        Coord = node[iPoint]->GetCoord();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[iDim*nVertexS+iVertex] = Coord[iDim];
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Coord,nBufferS_Vector,MPI_DOUBLE,send_to,0,
                   Buffer_Receive_Coord,nBufferR_Vector,MPI_DOUBLE,receive_from,0,MPI_COMM_WORLD,&status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_Coord[iDim*nVertexR+iVertex] = Buffer_Send_Coord[iDim*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_Coord;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        
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
        
        /*--- Copy coordinates before performing transformation. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          newCoord[iDim] = Buffer_Receive_Coord[iDim*nVertexR+iVertex];
        
        /*--- Rotate the coordinates. ---*/
        
        if (nDim == 2) {
          newCoord[0] = (rotMatrix[0][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord[1*nVertexR+iVertex]);
          newCoord[1] = (rotMatrix[1][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord[1*nVertexR+iVertex]);
        }
        else {
          newCoord[0] = (rotMatrix[0][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[0][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[0][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
          newCoord[1] = (rotMatrix[1][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[1][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[1][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
          newCoord[2] = (rotMatrix[2][0]*Buffer_Receive_Coord[0*nVertexR+iVertex] +
                         rotMatrix[2][1]*Buffer_Receive_Coord[1*nVertexR+iVertex] +
                         rotMatrix[2][2]*Buffer_Receive_Coord[2*nVertexR+iVertex]);
        }
        
        /*--- Copy transformed coordinates back into buffer. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          node[iPoint]->SetCoord(iDim, newCoord[iDim]);
        
      }
      
      /*--- Deallocate receive buffer. ---*/
      
      delete [] Buffer_Receive_Coord;
      
    }
    
  }
  
  delete [] newCoord;
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}

void CPhysicalGeometry::Set_MPI_GridVel(CConfig *config)  {
  
  unsigned short iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_GridVel = NULL, *Buffer_Send_GridVel = NULL, *GridVel = NULL, *newGridVel = NULL;
  int send_to, receive_from;
  
  newGridVel = new double[nDim];
  
#ifdef HAVE_MPI
  MPI_Status status;
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
      nVertexS = nVertex[MarkerS];  nVertexR = nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nDim;        nBufferR_Vector = nVertexR*nDim;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_GridVel = new double [nBufferR_Vector];
      Buffer_Send_GridVel = new double[nBufferS_Vector];
      
      /*--- Copy the grid velocity that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = vertex[MarkerS][iVertex]->GetNode();
        GridVel = node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_GridVel[iDim*nVertexS+iVertex] = GridVel[iDim];
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_GridVel,nBufferS_Vector,MPI_DOUBLE,send_to,0,
                   Buffer_Receive_GridVel,nBufferR_Vector,MPI_DOUBLE,receive_from,0,MPI_COMM_WORLD,&status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Receive_GridVel[iDim*nVertexR+iVertex] = Buffer_Send_GridVel[iDim*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_GridVel;
      
      /*--- Do the coordinate transformation ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        
        iPoint = vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = vertex[MarkerR][iVertex]->GetRotation_Type();
        
        /*--- Retrieve the supplied periodic information. ---*/
        
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
        
        /*--- Copy grid velocity before performing transformation. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          newGridVel[iDim] = Buffer_Receive_GridVel[iDim*nVertexR+iVertex];
        
        if (nDim == 2) {
          newGridVel[0] = (rotMatrix[0][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[0][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex]);
          newGridVel[1] = (rotMatrix[1][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[1][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex]);
        }
        else {
          newGridVel[0] = (rotMatrix[0][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[0][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[0][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
          newGridVel[1] = (rotMatrix[1][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[1][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[1][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
          newGridVel[2] = (rotMatrix[2][0]*Buffer_Receive_GridVel[0*nVertexR+iVertex] +
                           rotMatrix[2][1]*Buffer_Receive_GridVel[1*nVertexR+iVertex] +
                           rotMatrix[2][2]*Buffer_Receive_GridVel[2*nVertexR+iVertex]);
        }
        
        /*--- Copy transformed grid velocity back into buffer. ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          node[iPoint]->SetGridVel(iDim, newGridVel[iDim]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive_GridVel;
      
    }
    
  }
  
  delete [] newGridVel;
  
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
}

void CPhysicalGeometry::SetPeriodicBoundary(CConfig *config) {
  
  unsigned short iMarker, jMarker, kMarker = 0, iPeriodic, iDim, nPeriodic = 0, VTK_Type;
  unsigned long iNode, iIndex, iVertex, iPoint, iElem, kElem;
  unsigned long jElem, kPoint = 0, jVertex = 0, jPoint = 0, pPoint = 0, nPointPeriodic, newNodes[4];
  vector<unsigned long>::iterator IterElem, IterPoint[MAX_NUMBER_PERIODIC][2], IterNewElem[MAX_NUMBER_MARKER];
  double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  dx, dy, dz, rotCoord[3], epsilon = 1e-10, mindist = 1e6, *Coord_i, *Coord_j, dist = 0.0;
  bool isBadMatch = false;
  
  /*--- It only create the mirror structure for the second boundary ---*/
  bool CreateMirror[10];
  CreateMirror[1] = false;
  CreateMirror[2] = true;
  
  /*--- Send an initial message to the console. ---*/
  cout << "Setting the periodic boundary conditions." << endl;
  
  /*--- Loop through each marker to find any periodic boundaries. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      
      /*--- Evaluate the number of periodic boundary conditions defined
       in the geometry file ---*/
      nPeriodic++;
      
      /*--- Get marker index of the periodic donor boundary. ---*/
      jMarker = config->GetMarker_Periodic_Donor(config->GetMarker_All_TagBound(iMarker));
      
      /*--- Write some info to the console. ---*/
      cout << "Checking " << config->GetMarker_All_TagBound(iMarker);
      cout << " boundary against periodic donor, " << config->GetMarker_All_TagBound(jMarker) << ". ";
      
      /*--- Retrieve the supplied periodic information. ---*/
      center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
      angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
      trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
      
      /*--- Store (center+trans) as it is constant and will be added on. ---*/
      translation[0] = center[0] + trans[0];
      translation[1] = center[1] + trans[1];
      translation[2] = center[2] + trans[2];
      
      /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
      theta = angles[0];
      phi   = angles[1];
      psi   = angles[2];
      
      cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
      sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
      
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
      
      /*--- Loop through all vertices and find/set the periodic point. ---*/
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- Retrieve node information for this boundary point. ---*/
        iPoint  = vertex[iMarker][iVertex]->GetNode();
        Coord_i = node[iPoint]->GetCoord();
        
        /*--- Get the position vector from rot center to point. ---*/
        dx = Coord_i[0] - center[0];
        dy = Coord_i[1] - center[1];
        if (nDim == 3) {
          dz = Coord_i[2] - center[2];
        } else {
          dz = 0.0;
        }
        
        /*--- Compute transformed point coordinates. ---*/
        rotCoord[0] = rotMatrix[0][0]*dx
        + rotMatrix[0][1]*dy
        + rotMatrix[0][2]*dz + translation[0];
        
        rotCoord[1] = rotMatrix[1][0]*dx
        + rotMatrix[1][1]*dy
        + rotMatrix[1][2]*dz + translation[1];
        
        rotCoord[2] = rotMatrix[2][0]*dx
        + rotMatrix[2][1]*dy
        + rotMatrix[2][2]*dz + translation[2];
        
        /*--- Perform a search to find the closest donor point. ---*/
        mindist = 1e10;
        for (jVertex = 0; jVertex < nVertex[jMarker]; jVertex++) {
          
          /*--- Retrieve information for this jPoint. ---*/
          jPoint = vertex[jMarker][jVertex]->GetNode();
          Coord_j = node[jPoint]->GetCoord();
          
          /*--- Check the distance between the computed periodic
           location and this jPoint. ---*/
          dist = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            dist += (Coord_j[iDim]-rotCoord[iDim])*(Coord_j[iDim]-rotCoord[iDim]);
          }
          dist = sqrt(dist);
          
          /*---  Store vertex information if this is the closest
           point found thus far. ---*/
          if (dist < mindist) { mindist = dist; pPoint = jPoint; }
        }
        
        /*--- Set the periodic point for this iPoint. ---*/
        vertex[iMarker][iVertex]->SetDonorPoint(pPoint);
        
        /*--- Print warning if the nearest point was not within
         the specified tolerance. Computation will continue. ---*/
        if (mindist > epsilon) {
          isBadMatch = true;
          cout.precision(10);
          cout << endl;
          cout << "   Bad match for point " << iPoint << ".\tNearest";
          cout << " donor distance: " << scientific << mindist << ".";
        }
      }
      
      /*--- Print final warning when finding bad matches. ---*/
      if (isBadMatch) {
        cout << endl;
        cout << "\n !!! Warning !!!" << endl;
        cout << "Bad matches found. Computation will continue, but be cautious.\n";
      }
      cout << endl;
      isBadMatch = false;
      
    }
  
  /*--- Create a vector to identify the points that belong to each periodic boundary condition ---*/
  bool *PeriodicBC = new bool [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) PeriodicBC[iPoint] = false;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        PeriodicBC[iPoint] = true;
      }
  
  /*--- Determine the new points that must be added to each periodic boundary,
   note that only one of the boundaries require the extra data ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      
      /*--- An integer identify the periodic boundary condition --*/
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- iPoint is the original point on the surface and jPoint is the
         equivalent point in the other periodic surface ---*/
        iPoint = vertex[iMarker][iVertex]->GetNode();
        jPoint = vertex[iMarker][iVertex]->GetDonorPoint();
        
        /*--- First the case in which it is necessary to create a mirror set of elements ---*/
        if (CreateMirror[iPeriodic]) {
          /*--- Now we must determine the neighbor points (including indirect ones) to the periodic points
           and store all the information (in this list we do not include the points
           that already belong to the periodic boundary), we also add the elements that
           share a point with the periodic boundary condition ---*/
          for (iIndex = 0; iIndex < node[jPoint]->GetnElem(); iIndex++) {
            iElem = node[jPoint]->GetElem(iIndex);
            PeriodicElem[iPeriodic].push_back(iElem);
            for (unsigned short iNode = 0; iNode <	elem[iElem]->GetnNodes(); iNode ++) {
              kPoint = elem[iElem]->GetNode(iNode);
              if (!PeriodicBC[kPoint]) PeriodicPoint[iPeriodic][0].push_back(kPoint);
            }
          }
        }
        /*--- Second the case where no new element is added, neither points ---*/
        else {
          PeriodicPoint[iPeriodic][0].push_back(jPoint);
          PeriodicPoint[iPeriodic][1].push_back(iPoint);
        }
      }
    }
  }
  
  /*--- Sort the points that must be sended and delete repeated points ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      sort( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
      IterPoint[iPeriodic][0] = unique( PeriodicPoint[iPeriodic][0].begin(), PeriodicPoint[iPeriodic][0].end());
      PeriodicPoint[iPeriodic][0].resize( IterPoint[iPeriodic][0] - PeriodicPoint[iPeriodic][0].begin() );
    }
  }
  
  /*--- Create a list of the points that receive the values (only the new points) ---*/
  nPointPeriodic = nPoint;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iPoint = 0; iPoint < PeriodicPoint[iPeriodic][0].size(); iPoint++) {
        PeriodicPoint[iPeriodic][1].push_back(nPointPeriodic);
        nPointPeriodic++;
      }
    }
  }
  
  /*--- Sort the elements that must be replicated in the periodic boundary
   and delete the repeated elements ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      sort( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
      IterElem = unique( PeriodicElem[iPeriodic].begin(), PeriodicElem[iPeriodic].end());
      PeriodicElem[iPeriodic].resize( IterElem - PeriodicElem[iPeriodic].begin() );
    }
  }
  
  /*--- Check all SEND points to see if they also lie on another boundary. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      
      /*--- iPoint is a node that lies on the current marker. ---*/
      iPoint = vertex[iMarker][iVertex]->GetNode();
      
      /*--- Search through SEND points to check for iPoint. ---*/
      for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
        if (CreateMirror[iPeriodic]) {
          
          /*--- jPoint is the SEND point. ---*/
          for (iElem = 0; iElem < PeriodicPoint[iPeriodic][0].size(); iElem++) {
            jPoint = PeriodicPoint[iPeriodic][0][iElem];
            
            /*--- If the two match, then jPoint lies on this boundary.
             However, we are concerned with the new points, so we
             will store kPoint instead. ---*/
            if (iPoint == jPoint) {
              kPoint = PeriodicPoint[iPeriodic][1][iElem];
              
              /*--- We also want the type of boundary element that this point
               was within, so that we know what type of element to add
               built from the new points. ---*/
              bool isJPoint, isPeriodic;
              for(jElem = 0; jElem < nElem_Bound[iMarker]; jElem++) {
                isJPoint = false; isPeriodic = false;
                for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
                  if (bound[iMarker][jElem]->GetNode(iNode) == jPoint) isJPoint = true;
                  if (PeriodicBC[bound[iMarker][jElem]->GetNode(iNode)]) isPeriodic = true;
                }
                
                /*--- If both points were found, store this element. ---*/
                if (isJPoint && isPeriodic) {
                  OldBoundaryElems[iMarker].push_back(jElem);
                }
                
              }
              
            }
          }
        }
      }
    }
  }
  
  /*--- Sort the elements that must be added and remove duplicates. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    sort( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
    IterNewElem[iMarker] = unique( OldBoundaryElems[iMarker].begin(), OldBoundaryElems[iMarker].end());
    OldBoundaryElems[iMarker].resize( IterNewElem[iMarker] - OldBoundaryElems[iMarker].begin() );
  }
  
  /*--- Create the new boundary elements. Points making up these new
   elements must either be SEND/RECEIVE or periodic points. ---*/
  nNewElem_Bound = new unsigned long[nMarker];
  newBound = new CPrimalGrid**[nMarker];
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    nNewElem_Bound[iMarker] = OldBoundaryElems[iMarker].size();
    newBound[iMarker]       = new CPrimalGrid*[nNewElem_Bound[iMarker]];
    
    /*--- Loop through all new elements to be added. ---*/
    for (iElem = 0; iElem < nNewElem_Bound[iMarker]; iElem++) {
      jElem = OldBoundaryElems[iMarker][iElem];
      
      /*--- Loop through all nodes of this element. ---*/
      for (iNode = 0; iNode < bound[iMarker][jElem]->GetnNodes(); iNode++) {
        pPoint = bound[iMarker][jElem]->GetNode(iNode);
        
        /*--- Check if this node is a send point. If so, the corresponding
         receive point will be used in making the new boundary element. ---*/
        for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
          for (kElem = 0; kElem < PeriodicPoint[iPeriodic][0].size(); kElem++) {
            if (pPoint == PeriodicPoint[iPeriodic][0][kElem]) newNodes[iNode] = PeriodicPoint[iPeriodic][1][kElem];
          }
        }
        
        /*--- Check if this node is a periodic point. If so, the corresponding
         periodic point will be used in making the new boundary element. ---*/
        if (PeriodicBC[pPoint]) {
          
          /*--- Find the corresponding periodic point. ---*/
          for (jMarker = 0; jMarker < config->GetnMarker_All(); jMarker++) {
            if (config->GetMarker_All_KindBC(jMarker) == PERIODIC_BOUNDARY) {
              for (iVertex = 0; iVertex < nVertex[jMarker]; iVertex++) {
                if (pPoint == vertex[jMarker][iVertex]->GetNode()) {kMarker = jMarker; jVertex = iVertex;}
              }
            }
          }
          newNodes[iNode] = vertex[kMarker][jVertex]->GetDonorPoint();
        }
      }
      
      /*--- Now instantiate the new element. ---*/
      VTK_Type = bound[iMarker][jElem]->GetVTK_Type();
      switch(VTK_Type) {
        case LINE:
          newBound[iMarker][iElem] = new CLine(newNodes[0],newNodes[1],2);
          break;
        case TRIANGLE:
          newBound[iMarker][iElem] = new CTriangle(newNodes[0],newNodes[1],newNodes[2],3);
          break;
        case RECTANGLE:
          newBound[iMarker][iElem] = new CRectangle(newNodes[0],newNodes[1],newNodes[2],newNodes[3],3);
          break;
      }
      
    }
  }
  
  delete[] PeriodicBC;
  
}

void CPhysicalGeometry::FindNormal_Neighbor(CConfig *config) {
  double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
  unsigned long Point_Normal, jPoint;
  unsigned short iNeigh, iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Normal = vertex[iMarker][iVertex]->GetNormal();
        
        /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
        Point_Normal = 0; cos_max = -1.0;
        for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
          jPoint = node[iPoint]->GetPoint(iNeigh);
          scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
          for(iDim = 0; iDim < nDim; iDim++) {
            diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
            scalar_prod += diff_coord*Normal[iDim];
            norm_vect += diff_coord*diff_coord;
            norm_Normal += Normal[iDim]*Normal[iDim];
          }
          norm_vect = sqrt(norm_vect);
          norm_Normal = sqrt(norm_Normal);
          cos_alpha = scalar_prod/(norm_vect*norm_Normal);
          
          /*--- Get maximum cosine ---*/
          if (cos_alpha >= cos_max) {
            Point_Normal = jPoint;
            cos_max = cos_alpha;
          }
        }
        vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
      }
    }
  }
}

void CPhysicalGeometry::SetGeometryPlanes(CConfig *config) {
  
  bool loop_on;
  unsigned short iMarker = 0;
  double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new double[nVertex_Wall];
  Ycoord = new double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new double[nVertex_Wall];
  FaceArea = new double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<double> XCoordList;
  vector<double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (nDim==3) delete[] Zcoord;
  delete[] FaceArea;
}

CMultiGridGeometry::CMultiGridGeometry(CGeometry ***geometry, CConfig **config_container, unsigned short iMesh, unsigned short iZone) : CGeometry() {
  
  /*--- CGeometry & CConfig pointers to the fine grid level for clarity. We may
   need access to the other zones in the mesh for zone boundaries. ---*/
  
  CGeometry *fine_grid = geometry[iZone][iMesh-1];
  CConfig *config = config_container[iZone];
  
  /*--- Local variables ---*/
  
  unsigned long iPoint, Index_CoarseCV, CVPoint, iElem, iVertex, jPoint, iteration, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, iParent, jVertex, *Buffer_Receive_Parent = NULL, *Buffer_Send_Parent = NULL, *Buffer_Receive_Children = NULL, *Buffer_Send_Children = NULL, *Parent_Remote = NULL, *Children_Remote = NULL, *Parent_Local = NULL, *Children_Local = NULL, Local_nPointCoarse, Local_nPointFine, Global_nPointCoarse, Global_nPointFine;;
  short marker_seed;
  int send_to, receive_from, rank;
  bool agglomerate_seed = true;
  unsigned short nChildren, iNode, counter, iMarker, jMarker, Marker_Boundary, priority, copy_marker[MAX_NUMBER_MARKER], MarkerS, MarkerR, *nChildren_MPI;
  vector<unsigned long> Suitable_Indirect_Neighbors, Aux_Parent;
  vector<unsigned long>::iterator it;
  
#ifndef HAVE_MPI
  rank = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;
#endif
  
  FinestMGLevel = false; // Set the boolean to indicate that this is a coarse multigrid level.
  nDim = fine_grid->GetnDim(); // Write the number of dimensions of the coarse grid.
  
  /*--- Create a queue system to deo the agglomeration
   1st) More than two markers ---> Vertices (never agglomerate)
   2nd) Two markers ---> Edges (agglomerate if same BC, never agglomerate if different BC)
   3rd) One marker ---> Surface (always agglomarate)
   4th) No marker ---> Internal Volume (always agglomarate) ---*/
  
  /*--- Set a marker to indicate indirect agglomeration ---*/
  
  if (iMesh == MESH_1) {
    
    for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++)
      fine_grid->node[iPoint]->SetAgglomerate_Indirect(false);
    
    for (iElem = 0; iElem < fine_grid->GetnElem(); iElem++) {
      if ((fine_grid->elem[iElem]->GetVTK_Type() == HEXAHEDRON) ||
          (fine_grid->elem[iElem]->GetVTK_Type() == RECTANGLE)) {
        for (iNode = 0; iNode < fine_grid->elem[iElem]->GetnNodes(); iNode++) {
          iPoint = fine_grid->elem[iElem]->GetNode(iNode);
          fine_grid->node[iPoint]->SetAgglomerate_Indirect(true);
        }
      }
    }
    
  }
  
  /*--- Create the coarse grid structure using as baseline the fine grid ---*/
  
  CMultiGridQueue MGQueue_InnerCV(fine_grid->GetnPoint());
  
  node = new CPoint*[fine_grid->GetnPoint()];
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- Create node structure ---*/
    
    node[iPoint] = new CPoint(nDim, iPoint, config);
    
    /*--- Set the indirect agglomeration to false ---*/
    
    node[iPoint]->SetAgglomerate_Indirect(false);
  }
  
  Index_CoarseCV = 0;
  
  /*--- The first step is the boundary agglomeration. ---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++) {
    Marker_Boundary = config->GetMarker_All_KindBC(iMarker);
    
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      
      /*--- If the element has not being previously agglomerated and it belongs
       to the physical domain, then the agglomeration is studied ---*/
      
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain()) &&
          (GeometricalCheck(iPoint, fine_grid, config))) {
        
        nChildren = 1;
        
        /*--- We set an index for the parent control volume ---*/
        
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        
        /*--- We add the seed point (child) to the parent control volume ---*/
        
        node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
        agglomerate_seed = true; counter = 0; marker_seed = iMarker;
        
        /*--- For a particular point in the fine grid we save all the markers
         that are in that point ---*/
        
        for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
          if (fine_grid->node[iPoint]->GetVertex(jMarker) != -1) {
            copy_marker[counter] = jMarker;
            counter++;
          }
        
        /*--- To aglomerate a vertex it must have only one physical bc!!
         This can be improved. If there is only a marker, it is a good
         candidate for agglomeration ---*/
        
        if (counter == 1) agglomerate_seed = true;
        
        /*--- If there are two markers, we will aglomerate if one of the
         marker is SEND_RECEIVE ---*/
        
        if (counter == 2) {
          if ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) ||
              (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) agglomerate_seed = true;
          else agglomerate_seed = false;
        }
        
        /*--- If there are more than 2 markers, the aglomeration will be discarted ---*/
        
        if (counter > 2) agglomerate_seed = false;
        
        /*--- If the seed can be agglomerated, we try to agglomerate more points ---*/
        
        if (agglomerate_seed) {
          
          /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
          
          for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
            
            CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
            
            /*--- The new point can be agglomerated ---*/
            
            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config))  {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
              nChildren++;
            }
            
          }
          
          Suitable_Indirect_Neighbors.clear();
          
          if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
            SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
          
          /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
          
          for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
            
            CVPoint = Suitable_Indirect_Neighbors[iNode];
            
            /*--- The new point can be agglomerated ---*/
            
            if (SetBoundAgglomeration(CVPoint, marker_seed, fine_grid, config))  {
              
              /*--- We set the value of the parent ---*/
              
              fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
              
              /*--- We set the indirect agglomeration information ---*/
              
              if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
                node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
              
              /*--- We set the value of the child ---*/
              
              node[Index_CoarseCV]->SetChildren_CV(nChildren,CVPoint);
              nChildren++;
            }
          }
          
          
        }
        
        /*--- Update the number of child of the control volume ---*/
        
        node[Index_CoarseCV]->SetnChildren_CV(nChildren);
        Index_CoarseCV++;
      }
    }
  }
  
  /*--- Agglomerate all the nodes that have more than one physical boundary condition,
   Maybe here we can add the posibility of merging the vertex that have the same number,
   and kind  of markers---*/
  
  for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker++)
    for (iVertex = 0; iVertex < fine_grid->GetnVertex(iMarker); iVertex++) {
      iPoint = fine_grid->vertex[iMarker][iVertex]->GetNode();
      if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
          (fine_grid->node[iPoint]->GetDomain())) {
        fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
        node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
        node[Index_CoarseCV]->SetnChildren_CV(1);
        Index_CoarseCV++;
      }
    }
  
  /*--- Update the queue with the results from the boundary agglomeration ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    
    /*--- The CV has been agglomerated, remove form the list ---*/
    
    if (fine_grid->node[iPoint]->GetAgglomerate() == true) {
      
      MGQueue_InnerCV.RemoveCV(iPoint);
      
    }
    
    else {
      
      /*--- Count the number of agglomerated neighbors, and modify the queue ---*/
      
      priority = 0;
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        if (fine_grid->node[jPoint]->GetAgglomerate() == true) priority++;
      }
      MGQueue_InnerCV.MoveCV(iPoint, priority);
    }
  }
  
  /*--- Agglomerate the domain nodes ---*/
  
  iteration = 0;
  while (!MGQueue_InnerCV.EmptyQueue() && (iteration < fine_grid->GetnPoint())) {
    
    iPoint = MGQueue_InnerCV.NextCV();
    iteration ++;
    
    /*--- If the element has not being previously agglomerated, belongs to the physical domain,
     and satisfies several geometrical criteria then the seed CV is acepted for agglomeration ---*/
    
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) &&
        (fine_grid->node[iPoint]->GetDomain()) &&
        (GeometricalCheck(iPoint, fine_grid, config))) {
      
      nChildren = 1;
      
      /*--- We set an index for the parent control volume ---*/
      
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      
      /*--- We add the seed point (child) to the parent control volume ---*/
      
      node[Index_CoarseCV]->SetChildren_CV(0, iPoint);
      
      /*--- Update the queue with the seed point (remove the seed and
       increase the priority of the neighbors) ---*/
      
      MGQueue_InnerCV.Update(iPoint, fine_grid);
      
      /*--- Now we do a sweep over all the nodes that surround the seed point ---*/
      
      for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
        
        CVPoint = fine_grid->node[iPoint]->GetPoint(iNode);
        
        /*--- Determine if the CVPoint can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[CVPoint]->GetDomain()) &&
            (GeometricalCheck(CVPoint, fine_grid, config))) {
          
          /*--- We set the value of the parent ---*/
          
          fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
          
          /*--- We set the value of the child ---*/
          
          node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
          nChildren++;
          
          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/
          
          MGQueue_InnerCV.Update(CVPoint, fine_grid);
          
        }
        
      }
      
      /*--- Subrotuine to identify the indirect neighbors ---*/
      
      Suitable_Indirect_Neighbors.clear();
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        SetSuitableNeighbors(&Suitable_Indirect_Neighbors, iPoint, Index_CoarseCV, fine_grid);
      
      /*--- Now we do a sweep over all the indirect nodes that can be added ---*/
      
      for (iNode = 0; iNode <	Suitable_Indirect_Neighbors.size(); iNode ++) {
        
        CVPoint = Suitable_Indirect_Neighbors[iNode];
        
        /*--- The new point can be agglomerated ---*/
        
        if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
            (fine_grid->node[CVPoint]->GetDomain()))  {
          
          /*--- We set the value of the parent ---*/
          
          fine_grid->node[CVPoint]->SetParent_CV(Index_CoarseCV);
          
          /*--- We set the indirect agglomeration information ---*/
          
          if (fine_grid->node[CVPoint]->GetAgglomerate_Indirect())
            node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
          
          /*--- We set the value of the child ---*/
          
          node[Index_CoarseCV]->SetChildren_CV(nChildren, CVPoint);
          nChildren++;
          
          /*--- Update the queue with the new control volume (remove the CV and
           increase the priority of the neighbors) ---*/
          
          MGQueue_InnerCV.Update(CVPoint, fine_grid);
          
        }
      }
      
      /*--- Update the number of control of childrens ---*/
      
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
    }
    else {
      
      /*--- The seed point can not be agglomerated because of size, domain, streching, etc.
       move the point to the lowest priority ---*/
      
      MGQueue_InnerCV.MoveCV(iPoint, -1);
    }
    
  }
  
  /*--- Add all the elements that have not being agglomerated, in the previous stage ---*/
  
  for (iPoint = 0; iPoint < fine_grid->GetnPoint(); iPoint ++) {
    if ((fine_grid->node[iPoint]->GetAgglomerate() == false) && (fine_grid->node[iPoint]->GetDomain())) {
      
      nChildren = 1;
      fine_grid->node[iPoint]->SetParent_CV(Index_CoarseCV);
      if (fine_grid->node[iPoint]->GetAgglomerate_Indirect())
        node[Index_CoarseCV]->SetAgglomerate_Indirect(true);
      node[Index_CoarseCV]->SetChildren_CV(0,iPoint);
      node[Index_CoarseCV]->SetnChildren_CV(nChildren);
      Index_CoarseCV++;
      
    }
  }
  
  nPointDomain = Index_CoarseCV;
  
  /*--- Check that there are no hanging nodes ---*/
  
  unsigned long iFinePoint, iFinePoint_Neighbor, iCoarsePoint, iCoarsePoint_Complete;
  unsigned short iChildren;
  
  /*--- Find the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
      }
    }
  }
  
  /*--- Detect isolated points and merge them with its correct neighbor ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    
    if (node[iCoarsePoint]->GetnPoint() == 1) {
      
      /*--- Find the neighbor of the isolated point. This neighbor is the right control volume ---*/
      
      iCoarsePoint_Complete = node[iCoarsePoint]->GetPoint(0);
      
      /*--- Add the children to the connected control volume (and modify it parent indexing).
       Identify the child CV from the finest grid and added to the correct control volume.
       Set the parent CV of iFinePoint. Instead of using the original
       (iCoarsePoint) one use the new one (iCoarsePoint_Complete) ---*/
      
      nChildren = node[iCoarsePoint_Complete]->GetnChildren_CV();
      
      for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        node[iCoarsePoint_Complete]->SetChildren_CV(nChildren, iFinePoint);
        nChildren++;
        fine_grid->node[iFinePoint]->SetParent_CV(iCoarsePoint_Complete);
      }
      
      /*--- Update the number of children control volumes ---*/
      
      node[iCoarsePoint_Complete]->SetnChildren_CV(nChildren);
      node[iCoarsePoint]->SetnChildren_CV(0);
      
    }
  }
  
  //  unsigned long iPointFree = nPointDomain-1;
  //  iCoarsePoint = 0;
  //
  //  do {
  //
  //    if (node[iCoarsePoint]->GetnChildren_CV() == 0) {
  //
  //      while (node[iPointFree]->GetnChildren_CV() == 0) {
  //        Index_CoarseCV--;
  //        iPointFree--;
  //      }
  //
  //      nChildren = node[iPointFree]->GetnChildren_CV();
  //      for (iChildren = 0; iChildren <  nChildren; iChildren ++) {
  //        iFinePoint = node[iPointFree]->GetChildren_CV(iChildren);
  //        node[iCoarsePoint]->SetChildren_CV(iChildren, iFinePoint);
  //        fine_grid->node[iFinePoint]->SetParent_CV(iCoarsePoint);
  //      }
  //      node[iCoarsePoint]->SetnChildren_CV(nChildren);
  //      node[iPointFree]->SetnChildren_CV(0);
  //
  //      Index_CoarseCV--;
  //      iPointFree--;
  //
  //    }
  //
  //    iCoarsePoint++;
  //
  //  } while ((iCoarsePoint-1) < Index_CoarseCV);
  //
  //  nPointDomain = Index_CoarseCV;
  
  /*--- Reset the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPointDomain; iCoarsePoint ++) {
    node[iCoarsePoint]->ResetPoint();
  }
  
  /*--- Dealing with MPI parallelization, the objective is that the received nodes must be agglomerated
   in the same way as the donor nodes. Send the node agglomeration information of the donor
   (parent and children), Sending only occurs with MPI ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
      
      MarkerS = iMarker;  MarkerR = iMarker+1;
      
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
      
      nVertexS = fine_grid->nVertex[MarkerS];   nVertexR = fine_grid->nVertex[MarkerR];
      nBufferS_Vector = nVertexS;               nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      
      Buffer_Receive_Children = new unsigned long [nBufferR_Vector];
      Buffer_Send_Children = new unsigned long [nBufferS_Vector];
      
      Buffer_Receive_Parent = new unsigned long [nBufferR_Vector];
      Buffer_Send_Parent = new unsigned long [nBufferS_Vector];
      
      /*--- Copy the information that should be sended ---*/
      
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = fine_grid->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Children[iVertex] = iPoint;
        Buffer_Send_Parent[iVertex] = fine_grid->node[iPoint]->GetParent_CV();
      }
      
#ifdef HAVE_MPI
      /*--- Send/Receive information using Sendrecv ---*/
      MPI_Sendrecv(Buffer_Send_Children,nBufferS_Vector,MPI_UNSIGNED_LONG,send_to,0,
                   Buffer_Receive_Children,nBufferR_Vector,MPI_UNSIGNED_LONG,receive_from,0,MPI_COMM_WORLD,&status);
      MPI_Sendrecv(Buffer_Send_Parent,nBufferS_Vector,MPI_UNSIGNED_LONG,send_to,1,
                   Buffer_Receive_Parent,nBufferR_Vector,MPI_UNSIGNED_LONG,receive_from,1,MPI_COMM_WORLD,&status);
#else
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Children[iVertex] = Buffer_Send_Children[iVertex];
        Buffer_Receive_Parent[iVertex] = Buffer_Send_Parent[iVertex];
      }
#endif
      
      /*--- Deallocate send buffer ---*/
      
      delete [] Buffer_Send_Children;
      delete [] Buffer_Send_Parent;
      
      /*--- Create a list of the parent nodes without repeated parents ---*/
      
      Aux_Parent.clear();
      for (iVertex = 0; iVertex < nVertexR; iVertex++)
        Aux_Parent.push_back (Buffer_Receive_Parent[iVertex]);
      
      sort(Aux_Parent.begin(), Aux_Parent.end());
      it = unique(Aux_Parent.begin(), Aux_Parent.end());
      Aux_Parent.resize(it - Aux_Parent.begin());
      
      /*--- Allocate some structures ---*/
      
      Parent_Remote = new unsigned long[nVertexR];
      Children_Remote = new unsigned long[nVertexR];
      Parent_Local = new unsigned long[nVertexR];
      Children_Local = new unsigned long[nVertexR];
      
      /*--- Create the local vector and remote for the parents and the children ---*/
      
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        Parent_Remote[iVertex] = Buffer_Receive_Parent[iVertex];
        
        /*--- We use the same sorting as in the donor domain ---*/
        
        for (jVertex = 0; jVertex < Aux_Parent.size(); jVertex++) {
          if (Parent_Remote[iVertex] == Aux_Parent[jVertex]) {
            Parent_Local[iVertex] = jVertex + Index_CoarseCV;
            break;
          }
        }
        
        Children_Remote[iVertex] = Buffer_Receive_Children[iVertex];
        Children_Local[iVertex] = fine_grid->vertex[MarkerR][iVertex]->GetNode();
        
      }
      
      Index_CoarseCV += Aux_Parent.size();
      
      nChildren_MPI = new unsigned short [Index_CoarseCV];
      for (iParent = 0; iParent < Index_CoarseCV; iParent++)
        nChildren_MPI[iParent] = 0;
      
      /*--- Create the final structure ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Be careful, it is possible that a node change the agglomeration configuration, the priority
         is always, when receive the information ---*/
        
        fine_grid->node[Children_Local[iVertex]]->SetParent_CV(Parent_Local[iVertex]);
        node[Parent_Local[iVertex]]->SetChildren_CV(nChildren_MPI[Parent_Local[iVertex]], Children_Local[iVertex]);
        nChildren_MPI[Parent_Local[iVertex]]++;
        node[Parent_Local[iVertex]]->SetnChildren_CV(nChildren_MPI[Parent_Local[iVertex]]);
        node[Parent_Local[iVertex]]->SetDomain(false);
        
      }
      
      /*--- Deallocate auxiliar structures ---*/
      
      delete[] nChildren_MPI;
      delete[] Parent_Remote;
      delete[] Children_Remote;
      delete[] Parent_Local;
      delete[] Children_Local;
      
      /*--- Deallocate receive buffer ---*/
      
      delete [] Buffer_Receive_Children;
      delete [] Buffer_Receive_Parent;
      
    }
    
  }
  
  /*--- Update the number of points after the MPI agglomeration ---*/
  
  nPoint = Index_CoarseCV;
  
  /*--- Console output with the summary of the agglomeration ---*/
  
  Local_nPointCoarse = nPoint;
  Local_nPointFine = fine_grid->GetnPoint();
  
#ifdef HAVE_MPI
  MPI_Allreduce(&Local_nPointCoarse, &Global_nPointCoarse, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Local_nPointFine, &Global_nPointFine, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nPointCoarse = Local_nPointCoarse;
  Global_nPointFine = Local_nPointFine;
#endif
  
  double Coeff = 1.0, CFL = 0.0, factor = 1.5;
  
  if (iMesh != MESH_0) {
    if (nDim == 2) Coeff = pow(double(Global_nPointFine)/double(Global_nPointCoarse), 1./2.);
    if (nDim == 3) Coeff = pow(double(Global_nPointFine)/double(Global_nPointCoarse), 1./3.);
    CFL = factor*config->GetCFL(iMesh-1)/Coeff;
    config->SetCFL(iMesh, CFL);
  }
  
  double ratio = double(Global_nPointFine)/double(Global_nPointCoarse);
  
  if (((nDim == 2) && (ratio < 2.5)) ||
      ((nDim == 3) && (ratio < 2.5))) {
    config->SetMGLevels(iMesh-1);
  }
  else {
    if (rank == MASTER_NODE) {
      if (iMesh == 1) cout <<"MG level: "<< iMesh-1 <<"-> CVs: " << Global_nPointFine << ". Agglomeration rate 1/1.00. CFL "<< config->GetCFL(iMesh-1) <<"." << endl;
      cout <<"MG level: "<< iMesh <<"-> CVs: " << Global_nPointCoarse << ". Agglomeration rate 1/" << ratio <<". CFL "<< CFL <<"." << endl;
    }
  }
  
}


CMultiGridGeometry::~CMultiGridGeometry(void) {
  
}

bool CMultiGridGeometry::SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config) {
  
  bool agglomerate_CV = false;
  unsigned short counter, jMarker;
  unsigned short copy_marker[MAX_NUMBER_MARKER];

  /*--- Basic condition, the element has not being previously agglomerated, it belongs to the domain, 
   and has passed some basic geometrical check ---*/
  
  if ((fine_grid->node[CVPoint]->GetAgglomerate() == false) &&
      (fine_grid->node[CVPoint]->GetDomain()) &&
      (GeometricalCheck(CVPoint, fine_grid, config))) {
    
    /*--- If the element belong to the boundary, we must be careful ---*/
    
    if (fine_grid->node[CVPoint]->GetBoundary()) {
      
      /*--- Identify the markers of the vertex that we want to agglomerate ---*/
      
      counter = 0;
      for (jMarker = 0; jMarker < fine_grid->GetnMarker(); jMarker ++)
        if (fine_grid->node[CVPoint]->GetVertex(jMarker) != -1) {
          copy_marker[counter] = jMarker;
          counter++;
        }
      
      /*--- The basic condition is that the aglomerated vertex must have the same physical marker,
       but eventually a send-receive condition ---*/
      
      /*--- Only one marker in the vertex that is going to be aglomerated ---*/
      
      if (counter == 1) {
        
        /*--- We agglomerate if there is only a marker and is the same marker as the seed marker ---*/
        
        if (copy_marker[0] == marker_seed)
          agglomerate_CV = true;
        
        /*--- If there is only a marker, but the marker is the SEND_RECEIVE ---*/
        
        if (config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE)
          agglomerate_CV = true;
        
      }
      
      /*--- If there are two markers in the vertex that is going to be aglomerated ---*/
      
      if (counter == 2) {
        
        /*--- First we verify that the seed is a physical boundary ---*/
        
        if (config->GetMarker_All_KindBC(marker_seed) != SEND_RECEIVE) {
          
          /*--- Then we check that one of the marker is equal to the seed marker, and the other is send/receive ---*/
          
          if (((copy_marker[0] == marker_seed) && (config->GetMarker_All_KindBC(copy_marker[1]) == SEND_RECEIVE)) ||
              ((config->GetMarker_All_KindBC(copy_marker[0]) == SEND_RECEIVE) && (copy_marker[1] == marker_seed)))
            agglomerate_CV = true;
        }
        
      }
      
    }
    
    /*--- If the element belong to the domain, it is allways aglomerated ---*/
    
    else { agglomerate_CV = true; }
    
  }
  
  return agglomerate_CV;
}


bool CMultiGridGeometry::GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config) {
  
  double max_dimension = 1.2;
  
  /*--- Evaluate the total size of the element ---*/
  
  bool Volume = true;
  double ratio = pow(fine_grid->node[iPoint]->GetVolume(), 1.0/double(nDim))*max_dimension;
  double limit = pow(config->GetDomainVolume(), 1.0/double(nDim));
  if ( ratio > limit ) Volume = false;
  
  /*--- Evaluate the stretching of the element ---*/
  
  bool Stretching = true;
  
  /*	unsigned short iNode, iDim;
   unsigned long jPoint;
   double *Coord_i = fine_grid->node[iPoint]->GetCoord();
   double max_dist = 0.0 ; double min_dist = 1E20;
   for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
   jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
   double *Coord_j = fine_grid->node[jPoint]->GetCoord();
   double distance = 0.0;
   for (iDim = 0; iDim < nDim; iDim++)
   distance += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
   distance = sqrt(distance);
   max_dist = max(distance, max_dist);
   min_dist = min(distance, min_dist);
   }
   if ( max_dist/min_dist > 100.0 ) Stretching = false;*/
  
  return (Stretching && Volume);
  
}

void CMultiGridGeometry::SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                                              unsigned long Index_CoarseCV, CGeometry *fine_grid) {
  
  unsigned long jPoint, kPoint, lPoint;
  unsigned short iNode, jNode, iNeighbor, jNeighbor, kNode;
  bool SecondNeighborSeed, ThirdNeighborSeed;
  vector<unsigned long>::iterator it;
  
  /*--- Create a list with the first neighbors, including the seed ---*/
  
  vector<unsigned long> First_Neighbor_Points;
  First_Neighbor_Points.push_back(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    First_Neighbor_Points.push_back(jPoint);
  }
  
  /*--- Create a list with the second neighbors, without first, and seed neighbors ---*/
  
  vector<unsigned long> Second_Neighbor_Points, Second_Origin_Points, Suitable_Second_Neighbors;
  
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    
    for (jNode = 0; jNode <	fine_grid->node[jPoint]->GetnPoint(); jNode ++) {
      kPoint = fine_grid->node[jPoint]->GetPoint(jNode);
      
      /*--- Check that the second neighbor do not belong to the first neighbor or the seed ---*/
      
      SecondNeighborSeed = true;
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (kPoint == First_Neighbor_Points[iNeighbor]) {
          SecondNeighborSeed = false; break;
        }
      
      if (SecondNeighborSeed) {
        Second_Neighbor_Points.push_back(kPoint);
        Second_Origin_Points.push_back(jPoint);
      }
      
    }
  }
  
  /*---  Identify those second neighbors that are repeated (candidate to be added) ---*/
  
  for (iNeighbor = 0; iNeighbor <	Second_Neighbor_Points.size(); iNeighbor ++)
    
    for (jNeighbor = 0; jNeighbor <	Second_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      
      if ((Second_Neighbor_Points[iNeighbor] == Second_Neighbor_Points[jNeighbor]) &&
          (Second_Origin_Points[iNeighbor] != Second_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        Suitable_Indirect_Neighbors->push_back(Second_Neighbor_Points[iNeighbor]);
        
        /*--- Create alist with the suitable second neighbor, that we will use
         to compute the third neighbors --*/
        
        Suitable_Second_Neighbors.push_back(Second_Neighbor_Points[iNeighbor]);
        
      }
  
  
  /*--- Remove repeated from the suitable second neighbors ---*/
  
  sort(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  it = unique(Suitable_Second_Neighbors.begin(), Suitable_Second_Neighbors.end());
  Suitable_Second_Neighbors.resize(it - Suitable_Second_Neighbors.begin());
  
  /*--- Remove repeated from first neighbors ---*/
  
  sort(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  it = unique(First_Neighbor_Points.begin(), First_Neighbor_Points.end());
  First_Neighbor_Points.resize(it - First_Neighbor_Points.begin());
  
  /*--- Create a list with the third neighbors, without first, second, and seed neighbors ---*/
  
  vector<unsigned long> Third_Neighbor_Points, Third_Origin_Points;
  
  for (jNode = 0; jNode <	Suitable_Second_Neighbors.size(); jNode ++) {
    kPoint = Suitable_Second_Neighbors[jNode];
    
    for (kNode = 0; kNode <	fine_grid->node[kPoint]->GetnPoint(); kNode ++) {
      lPoint = fine_grid->node[kPoint]->GetPoint(kNode);
      
      /*--- Check that the third neighbor do not belong to the first neighbors or the seed ---*/
      
      ThirdNeighborSeed = true;
      
      for (iNeighbor = 0; iNeighbor <	First_Neighbor_Points.size(); iNeighbor ++)
        if (lPoint == First_Neighbor_Points[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      /*--- Check that the third neighbor do not belong to the second neighbors ---*/
      
      for (iNeighbor = 0; iNeighbor <	Suitable_Second_Neighbors.size(); iNeighbor ++)
        if (lPoint == Suitable_Second_Neighbors[iNeighbor]) {
          ThirdNeighborSeed = false;
          break;
        }
      
      if (ThirdNeighborSeed) {
        Third_Neighbor_Points.push_back(lPoint);
        Third_Origin_Points.push_back(kPoint);
      }
      
    }
  }
  
  /*---  Identify those third neighbors that are repeated (candidate to be added) ---*/
  
  for (iNeighbor = 0; iNeighbor <	Third_Neighbor_Points.size(); iNeighbor ++)
    for (jNeighbor = 0; jNeighbor <	Third_Neighbor_Points.size(); jNeighbor ++)
      
    /*--- Repeated second neighbor with different origin ---*/
      
      if ((Third_Neighbor_Points[iNeighbor] == Third_Neighbor_Points[jNeighbor]) &&
          (Third_Origin_Points[iNeighbor] != Third_Origin_Points[jNeighbor]) &&
          (iNeighbor < jNeighbor)) {
        
        Suitable_Indirect_Neighbors->push_back(Third_Neighbor_Points[iNeighbor]);
        
      }
  
  /*--- Remove repeated from Suitable Indirect Neighbors List ---*/
  
  sort(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  it = unique(Suitable_Indirect_Neighbors->begin(), Suitable_Indirect_Neighbors->end());
  Suitable_Indirect_Neighbors->resize(it - Suitable_Indirect_Neighbors->begin());
  
}

void CMultiGridGeometry::SetPoint_Connectivity(CGeometry *fine_grid) {
  
  unsigned long iFinePoint, iFinePoint_Neighbor, iParent, iCoarsePoint;
  unsigned short iChildren, iNode;
  
  /*--- Set the point surrounding a point ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    for (iChildren = 0; iChildren <  node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if (iParent != iCoarsePoint) node[iCoarsePoint]->SetPoint(iParent);
      }
    }
  }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    node[iCoarsePoint]->SetnNeighbor(node[iCoarsePoint]->GetnPoint());
  
}

void CMultiGridGeometry::SetVertex(CGeometry *fine_grid, CConfig *config) {
  unsigned long  iVertex, iFinePoint, iCoarsePoint;
  unsigned short iMarker, iMarker_Tag, iChildren;
  
  nMarker = fine_grid->GetnMarker();
  
  /*--- If any children node belong to the boundary then the entire control
   volume will belong to the boundary ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      if (fine_grid->node[iFinePoint]->GetBoundary()) {
        node[iCoarsePoint]->SetBoundary(nMarker);
        break;
      }
    }
  
  vertex = new CVertex**[nMarker];
  nVertex = new unsigned long [nMarker];
  
  Tag_to_Marker = new string [MAX_NUMBER_MARKER];
  for (iMarker_Tag = 0; iMarker_Tag < MAX_NUMBER_MARKER; iMarker_Tag++)
    Tag_to_Marker[iMarker_Tag] = fine_grid->GetMarker_Tag(iMarker_Tag);
  
  /*--- Compute the number of vertices to do the dimensionalization ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    if (node[iCoarsePoint]->GetBoundary()) {
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < nMarker; iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            node[iCoarsePoint]->SetVertex(iVertex,iMarker);
            nVertex[iMarker]++;
          }
        }
      }
    }
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [fine_grid->GetnVertex(iMarker)+1];
    nVertex[iMarker] = 0;
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iMarker = 0; iMarker < nMarker; iMarker ++)
        node[iCoarsePoint]->SetVertex(-1,iMarker);
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) nVertex[iMarker] = 0;
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    if (node[iCoarsePoint]->GetBoundary())
      for (iChildren = 0; iChildren <	node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        for (iMarker = 0; iMarker < fine_grid->GetnMarker(); iMarker ++) {
          if ((fine_grid->node[iFinePoint]->GetVertex(iMarker) != -1) && (node[iCoarsePoint]->GetVertex(iMarker) == -1)) {
            iVertex = nVertex[iMarker];
            vertex[iMarker][iVertex] = new CVertex(iCoarsePoint, nDim);
            node[iCoarsePoint]->SetVertex(iVertex,iMarker);
            
            /*--- Set the transformation to apply ---*/
            unsigned long ChildVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
            unsigned short RotationKind = fine_grid->vertex[iMarker][ChildVertex]->GetRotation_Type();
            vertex[iMarker][iVertex]->SetRotation_Type(RotationKind);
            nVertex[iMarker]++;
          }
        }
      }
}

void CMultiGridGeometry::MatchNearField(CConfig *config) {
  
#ifndef HAVE_MPI
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
      }
  
#else
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, rank);
        }
      }
  
#endif
  
}

void CMultiGridGeometry::MatchActuator_Disk(CConfig *config) {
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int iProcessor;
  
#ifndef HAVE_MPI
  iProcessor = MASTER_NODE;
#else
  MPI_Comm_rank(MPI_COMM_WORLD, &iProcessor);
#endif
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) == ACTDISK_INLET) ||
        (config->GetMarker_All_KindBC(iMarker) == ACTDISK_OUTLET)) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iProcessor);
        }
      }
    }
  }
  
}

void CMultiGridGeometry::MatchInterface(CConfig *config) {
  
#ifndef HAVE_MPI
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        vertex[iMarker][iVertex]->SetDonorPoint(iPoint);
      }
  
#else
  
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        if (node[iPoint]->GetDomain()) {
          vertex[iMarker][iVertex]->SetDonorPoint(iPoint, rank);
        }
      }
  
#endif
  
}


void CMultiGridGeometry::SetControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  
  unsigned long iFinePoint,iFinePoint_Neighbor, iCoarsePoint, iEdge, iParent;
  long FineEdge, CoarseEdge;
  unsigned short iChildren, iNode, iDim;
  bool change_face_orientation;
  double *Normal, Coarse_Volume, Area, *NormalFace = NULL;
  Normal = new double [nDim];
  
  /*--- Compute the area of the coarse volume ---*/
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++) {
    node[iCoarsePoint]->SetVolume(0.0);
    Coarse_Volume = 0.0;
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      Coarse_Volume += fine_grid->node[iFinePoint]->GetVolume();
    }
    node[iCoarsePoint]->SetVolume(Coarse_Volume);
  }
  
  /*--- Update or not the values of faces at the edge ---*/
  if (action != ALLOCATE) {
    for(iEdge=0; iEdge < nEdge; iEdge++)
      edge[iEdge]->SetZeroValues();
  }
  
  for (iCoarsePoint = 0; iCoarsePoint < nPoint; iCoarsePoint ++)
    for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
      iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
      
      for (iNode = 0; iNode < fine_grid->node[iFinePoint]->GetnPoint(); iNode ++) {
        iFinePoint_Neighbor = fine_grid->node[iFinePoint]->GetPoint(iNode);
        iParent = fine_grid->node[iFinePoint_Neighbor]->GetParent_CV();
        if ((iParent != iCoarsePoint) && (iParent < iCoarsePoint)) {
          
          FineEdge = fine_grid->FindEdge(iFinePoint, iFinePoint_Neighbor);
          
          change_face_orientation = false;
          if (iFinePoint < iFinePoint_Neighbor) change_face_orientation = true;
          
          CoarseEdge = FindEdge(iParent, iCoarsePoint);
          
          fine_grid->edge[FineEdge]->GetNormal(Normal);
          
          if (change_face_orientation) {
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            edge[CoarseEdge]->AddNormal(Normal);
          }
          else {
            edge[CoarseEdge]->AddNormal(Normal);
          }
        }
      }
    }
  delete[] Normal;
  
  /*--- Check if there is a normal with null area ---*/
  
  for (iEdge = 0; iEdge < nEdge; iEdge++) {
    NormalFace = edge[iEdge]->GetNormal();
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
  }
  
}

void CMultiGridGeometry::SetBoundControlVolume(CConfig *config, CGeometry *fine_grid, unsigned short action) {
  unsigned long iCoarsePoint, iFinePoint, FineVertex, iVertex;
  unsigned short iMarker, iChildren, iDim;
  double *Normal, Area, *NormalFace = NULL;
  
  Normal = new double [nDim];
  
  if (action != ALLOCATE) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();
  }
  
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for(iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iCoarsePoint = vertex[iMarker][iVertex]->GetNode();
      for (iChildren = 0; iChildren < node[iCoarsePoint]->GetnChildren_CV(); iChildren ++) {
        iFinePoint = node[iCoarsePoint]->GetChildren_CV(iChildren);
        if (fine_grid->node[iFinePoint]->GetVertex(iMarker)!=-1) {
          FineVertex = fine_grid->node[iFinePoint]->GetVertex(iMarker);
          fine_grid->vertex[iMarker][FineVertex]->GetNormal(Normal);
          vertex[iMarker][iVertex]->AddNormal(Normal);
        }
      }
    }
  
  delete[] Normal;
  
  /*--- Check if there is a normal with null area ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }
  
}

void CMultiGridGeometry::SetCoord(CGeometry *geometry) {
  unsigned long Point_Fine, Point_Coarse;
  unsigned short iChildren, iDim;
  double Area_Parent, Area_Children;
  double *Coordinates_Fine, *Coordinates;
  Coordinates = new double[nDim];
  
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    Area_Parent = node[Point_Coarse]->GetVolume();
    for (iDim = 0; iDim < nDim; iDim++) Coordinates[iDim] = 0.0;
    for (iChildren = 0; iChildren < node[Point_Coarse]->GetnChildren_CV(); iChildren++) {
      Point_Fine = node[Point_Coarse]->GetChildren_CV(iChildren);
      Area_Children = geometry->node[Point_Fine]->GetVolume();
      Coordinates_Fine = geometry->node[Point_Fine]->GetCoord();
      for (iDim = 0; iDim < nDim; iDim++)
        Coordinates[iDim] += Coordinates_Fine[iDim]*Area_Children/Area_Parent;
    }
    for (iDim = 0; iDim < nDim; iDim++)
      node[Point_Coarse]->SetCoord(iDim,Coordinates[iDim]);
  }
  delete[] Coordinates;
}

void CMultiGridGeometry::SetRotationalVelocity(CConfig *config) {
  
  unsigned long iPoint_Coarse;
  double RotVel[3], Distance[3], *Coord, Center[3], Omega[3], L_Ref;
  
  /*--- Center of rotation & angular velocity vector from config. ---*/
  
  Center[0] = config->GetMotion_Origin_X(ZONE_0);
  Center[1] = config->GetMotion_Origin_Y(ZONE_0);
  Center[2] = config->GetMotion_Origin_Z(ZONE_0);
  Omega[0]  = config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref();
  Omega[1]  = config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref();
  Omega[2]  = config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref();
  L_Ref     = config->GetLength_Ref();
  
  /*--- Loop over all nodes and set the rotational velocity. ---*/
  
  for (iPoint_Coarse = 0; iPoint_Coarse < GetnPoint(); iPoint_Coarse++) {
    
    /*--- Get the coordinates of the current node ---*/
    
    Coord = node[iPoint_Coarse]->GetCoord();
    
    /*--- Calculate the non-dim. distance from the rotation center ---*/
    
    Distance[0] = (Coord[0]-Center[0])/L_Ref;
    Distance[1] = (Coord[1]-Center[1])/L_Ref;
    Distance[2] = (Coord[2]-Center[2])/L_Ref;
    
    /*--- Calculate the angular velocity as omega X r ---*/
    
    RotVel[0] = Omega[1]*(Distance[2]) - Omega[2]*(Distance[1]);
    RotVel[1] = Omega[2]*(Distance[0]) - Omega[0]*(Distance[2]);
    RotVel[2] = Omega[0]*(Distance[1]) - Omega[1]*(Distance[0]);
    
    /*--- Store the grid velocity at this node ---*/
    
    node[iPoint_Coarse]->SetGridVel(RotVel);
    
  }
  
}

void CMultiGridGeometry::SetGridVelocity(CConfig *config, unsigned long iter) {
  
  /*--- Local variables ---*/
  
  double *Coord_nP1 = NULL, *Coord_n = NULL, *Coord_nM1 = NULL;
  double TimeStep, GridVel = 0.0;
  unsigned long Point_Coarse;
  unsigned short iDim;
  
  /*--- Compute the velocity of each node in the volume mesh ---*/
  
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    
    /*--- Coordinates of the current point at n+1, n, & n-1 time levels ---*/
    
    Coord_nM1 = node[Point_Coarse]->GetCoord_n1();
    Coord_n   = node[Point_Coarse]->GetCoord_n();
    Coord_nP1 = node[Point_Coarse]->GetCoord();
    
    /*--- Unsteady time step ---*/
    
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute mesh velocity with 1st or 2nd-order approximation ---*/
    
    for(iDim = 0; iDim < nDim; iDim++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        GridVel = ( Coord_nP1[iDim] - Coord_n[iDim] ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        GridVel = ( 3.0*Coord_nP1[iDim] - 4.0*Coord_n[iDim]
                   +  1.0*Coord_nM1[iDim] ) / (2.0*TimeStep);
      
      /*--- Store grid velocity for this point ---*/
      
      node[Point_Coarse]->SetGridVel(iDim,GridVel);
      
    }
  }
}

void CMultiGridGeometry::SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) {
  
  /*--- Local variables ---*/
  unsigned short iDim, iChild;
  unsigned long Point_Coarse, Point_Fine;
  double Area_Parent, Area_Child, Grid_Vel[3], *Grid_Vel_Fine;
  
  /*--- Loop over all coarse mesh points ---*/
  for (Point_Coarse = 0; Point_Coarse < GetnPoint(); Point_Coarse++) {
    Area_Parent = node[Point_Coarse]->GetVolume();
    
    /*--- Zero out the grid velocity ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Grid_Vel[iDim] = 0.0;
    
    /*--- Loop over all of the children for this coarse CV and compute
     a grid velocity based on the values in the child CVs (fine mesh). ---*/
    for (iChild = 0; iChild < node[Point_Coarse]->GetnChildren_CV(); iChild++) {
      Point_Fine    = node[Point_Coarse]->GetChildren_CV(iChild);
      Area_Child    = fine_mesh->node[Point_Fine]->GetVolume();
      Grid_Vel_Fine = fine_mesh->node[Point_Fine]->GetGridVel();
      for (iDim = 0; iDim < nDim; iDim++)
        Grid_Vel[iDim] += Grid_Vel_Fine[iDim]*Area_Child/Area_Parent;
    }
    
    /*--- Set the grid velocity for this coarse node. ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      node[Point_Coarse]->SetGridVel(iDim, Grid_Vel[iDim]);
  }
}


void CMultiGridGeometry::FindNormal_Neighbor(CConfig *config) {
  
  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {
      
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        iPoint = vertex[iMarker][iVertex]->GetNode();
        
        /*--- If the node belong to the domain ---*/
        if (node[iPoint]->GetDomain()) {
          
          /*--- Compute closest normal neighbor ---*/
          double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord;
          unsigned long Point_Normal = 0, jPoint;
          unsigned short iNeigh;
          double *Normal = vertex[iMarker][iVertex]->GetNormal();
          cos_max = -1.0;
          for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
            jPoint = node[iPoint]->GetPoint(iNeigh);
            scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
            for(iDim = 0; iDim < nDim; iDim++) {
              diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
              scalar_prod += diff_coord*Normal[iDim];
              norm_vect += diff_coord*diff_coord;
              norm_Normal += Normal[iDim]*Normal[iDim];
            }
            norm_vect = sqrt(norm_vect);
            norm_Normal = sqrt(norm_Normal);
            cos_alpha = scalar_prod/(norm_vect*norm_Normal);
            
            /*--- Get maximum cosine (not minimum because normals are oriented inwards) ---*/
            if (cos_alpha >= cos_max) {
              Point_Normal = jPoint;
              cos_max = cos_alpha;
            }
          }
          vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
        }
      }
    }
  }
}


void CMultiGridGeometry::SetGeometryPlanes(CConfig *config) {
  bool loop_on;
  unsigned short iMarker = 0;
  double auxXCoord, auxYCoord, auxZCoord,	*Face_Normal = NULL, auxArea, *Xcoord = NULL, *Ycoord = NULL, *Zcoord = NULL, *FaceArea = NULL;
  unsigned long jVertex, iVertex,ixCoord, iPoint, iVertex_Wall, nVertex_Wall = 0;
  
  /*--- Compute the total number of points on the near-field ---*/
  nVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      nVertex_Wall += nVertex[iMarker];
  
  
  /*--- Create an array with all the coordinates, points, pressures, face area,
   equivalent area, and nearfield weight ---*/
  Xcoord = new double[nVertex_Wall];
  Ycoord = new double[nVertex_Wall];
  if (nDim == 3)	Zcoord = new double[nVertex_Wall];
  FaceArea = new double[nVertex_Wall];
  
  /*--- Copy the boundary information to an array ---*/
  iVertex_Wall = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX)               ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_CATALYTIC)     ||
        (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX_NONCATALYTIC)  ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)              ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_CATALYTIC)    ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL_NONCATALYTIC) ||
        (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)                )
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Xcoord[iVertex_Wall] = node[iPoint]->GetCoord(0);
        Ycoord[iVertex_Wall] = node[iPoint]->GetCoord(1);
        if (nDim==3) Zcoord[iVertex_Wall] = node[iPoint]->GetCoord(2);
        Face_Normal = vertex[iMarker][iVertex]->GetNormal();
        FaceArea[iVertex_Wall] = fabs(Face_Normal[nDim-1]);
        iVertex_Wall ++;
      }
  
  
  //vector<double> XCoordList;
  vector<double>::iterator IterXCoordList;
  
  for (iVertex = 0; iVertex < nVertex_Wall; iVertex++)
    XCoordList.push_back(Xcoord[iVertex]);
  
  sort( XCoordList.begin(), XCoordList.end());
  IterXCoordList = unique( XCoordList.begin(), XCoordList.end());
  XCoordList.resize( IterXCoordList - XCoordList.begin() );
  
  /*--- Create vectors and distribute the values among the different PhiAngle queues ---*/
  Xcoord_plane.resize(XCoordList.size());
  Ycoord_plane.resize(XCoordList.size());
  if (nDim==3) Zcoord_plane.resize(XCoordList.size());
  FaceArea_plane.resize(XCoordList.size());
  Plane_points.resize(XCoordList.size());
  
  
  double dist_ratio;
  unsigned long iCoord;
  
  /*--- Distribute the values among the different PhiAngles ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (node[iPoint]->GetDomain()) {
      loop_on = true;
      for (ixCoord = 0; ixCoord < XCoordList.size()-1 && loop_on; ixCoord++) {
        dist_ratio = (node[iPoint]->GetCoord(0) - XCoordList[ixCoord])/(XCoordList[ixCoord+1]- XCoordList[ixCoord]);
        if (dist_ratio >= 0 && dist_ratio <= 1.0) {
          if (dist_ratio <= 0.5) iCoord = ixCoord;
          else iCoord = ixCoord+1;
          Xcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(0) );
          Ycoord_plane[iCoord].push_back(node[iPoint]->GetCoord(1) );
          if (nDim==3) Zcoord_plane[iCoord].push_back(node[iPoint]->GetCoord(2) );
          FaceArea_plane[iCoord].push_back(node[iPoint]->GetVolume());   ///// CHECK AREA CALCULATION
          Plane_points[iCoord].push_back(iPoint );
          loop_on = false;
        }
      }
    }
  }
  
  unsigned long auxPoint;
  /*--- Order the arrays in ascending values of y ---*/
  for (ixCoord = 0; ixCoord < XCoordList.size(); ixCoord++)
    for (iVertex = 0; iVertex < Xcoord_plane[ixCoord].size(); iVertex++)
      for (jVertex = 0; jVertex < Xcoord_plane[ixCoord].size() - 1 - iVertex; jVertex++)
        if (Ycoord_plane[ixCoord][jVertex] > Ycoord_plane[ixCoord][jVertex+1]) {
          auxXCoord = Xcoord_plane[ixCoord][jVertex]; Xcoord_plane[ixCoord][jVertex] = Xcoord_plane[ixCoord][jVertex+1]; Xcoord_plane[ixCoord][jVertex+1] = auxXCoord;
          auxYCoord = Ycoord_plane[ixCoord][jVertex]; Ycoord_plane[ixCoord][jVertex] = Ycoord_plane[ixCoord][jVertex+1]; Ycoord_plane[ixCoord][jVertex+1] = auxYCoord;
          auxPoint = Plane_points[ixCoord][jVertex]; Plane_points[ixCoord][jVertex] = Plane_points[ixCoord][jVertex+1]; Plane_points[ixCoord][jVertex+1] = auxPoint;
          if (nDim==3) {
            auxZCoord = Zcoord_plane[ixCoord][jVertex]; Zcoord_plane[ixCoord][jVertex] = Zcoord_plane[ixCoord][jVertex+1]; Zcoord_plane[ixCoord][jVertex+1] = auxZCoord;
          }
          auxArea = FaceArea_plane[ixCoord][jVertex]; FaceArea_plane[ixCoord][jVertex] = FaceArea_plane[ixCoord][jVertex+1]; FaceArea_plane[ixCoord][jVertex+1] = auxArea;
        }
  
  /*--- Delete structures ---*/
  delete[] Xcoord; delete[] Ycoord;
  if (nDim==3) delete[] Zcoord;
  delete[] FaceArea;
}


CBoundaryGeometry::CBoundaryGeometry(CConfig *config, string val_mesh_filename, unsigned short val_format) : CGeometry() {
  
  string text_line;
  ifstream mesh_file;
  unsigned short iNode_Surface, VTK_Type, iMarker, iChar, iCount = 0, val_iZone = 1, val_nZone = 1;
  unsigned long Point_Surface, iElem_Surface, iElem_Bound = 0, iPoint = 0, iElem = 0, ielem = 0,
  nelem_edge = 0, nelem_triangle = 0, nelem_quad = 0,
  vnodes_edge[2], vnodes_triangle[3],
  vnodes_quad[4], dummy, GlobalIndex;
  string Marker_Tag;
  char cstr[MAX_STRING_SIZE];
  int rank = MASTER_NODE, size = SINGLE_NODE;
  bool domain_flag = false;
  bool found_transform = false;
  nZone = val_nZone;
  double Coord_2D[2], Coord_3D[3];
  string::size_type position;
  
#ifdef HAVE_MPI
  unsigned long LocalIndex;
  unsigned long Local_nPoint, Local_nPointDomain, Global_nPoint = 0;
  unsigned long Local_nElem, Global_nElem = 0;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  Global_nPointDomain = 0;
  FinestMGLevel = true;
  
  if(rank == MASTER_NODE)
    cout << endl <<"---------------------- Read grid file information -----------------------" << endl;
  
  
  /*--- Determine whether there are multiplze zones first ---*/
  strcpy (cstr, val_mesh_filename.c_str());
  mesh_file.open(cstr, ios::in);
  if (mesh_file.fail()) {
    cout << "There is no geometry file (CBoundaryGeometry)!" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- If more than one, find the domain in the mesh file ---*/
  if (val_nZone > 1) {
    while (getline (mesh_file,text_line)) {
      /*--- Search for the current domain ---*/
      position = text_line.find ("IZONE=",0);
      if (position != string::npos) {
        text_line.erase (0,6);
        unsigned short jDomain = atoi(text_line.c_str());
        if (jDomain == val_iZone) {
          if (rank == MASTER_NODE) cout << "Reading zone " << val_iZone << ":" << endl;
          break;
        }
      }
    }
  }
  
  /*--- Read grid file with format SU2 ---*/
  while (getline (mesh_file,text_line)) {
    
    /*--- Read the dimension of the problem ---*/
    position = text_line.find ("NDIME=",0);
    if (position != string::npos) {
      if (domain_flag == false) {
        text_line.erase (0,6); nDim = atoi(text_line.c_str());
        if (rank == MASTER_NODE) {
          if (nDim == 2) cout << "Two dimensional problem." << endl;
          if (nDim == 3) cout << "Three dimensional problem." << endl;
        }
        domain_flag = true;
      } else {
        break;
      }
    }
    
    /*--- Read the information about inner elements ---*/
    position = text_line.find ("NELEM=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nElem = atoi(text_line.c_str());
      if (size == 1)
        cout << nElem << " interior elements. ";
      
#ifdef HAVE_MPI
      Local_nElem = nElem;
      MPI_Allreduce(&Local_nElem, &Global_nElem, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
      while (iElem < nElem) {
        getline(mesh_file,text_line);
        iElem++;
      }
    }
    
    /*--- Read number of points ---*/
    position = text_line.find ("NPOIN=",0);
    if (position != string::npos) {
      text_line.erase (0,6);
      
      /*--- Check for ghost points. ---*/
      stringstream test_line(text_line);
      while (test_line >> dummy)
        iCount++;
      
      /*--- Now read and store the number of points and possible ghost points. ---*/
      stringstream  stream_line(text_line);
      if (iCount == 2) {
        stream_line >> nPoint;
        stream_line >> nPointDomain;
        if (size == 1)
          cout << nPoint << " points, and " << nPoint-nPointDomain << " ghost points." << endl;
        
        /*--- Set some important point information for parallel simulations. ---*/
#ifdef HAVE_MPI
        Local_nPoint = nPoint; Local_nPointDomain = nPointDomain;
        MPI_Allreduce(&Local_nPoint, &Global_nPoint, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
        Global_nPointDomain = nPointDomain;
#endif
      }
      else if (iCount == 1) {
        stream_line >> nPoint;
        nPointDomain = nPoint;
        Global_nPointDomain = nPoint;
        if (rank == MASTER_NODE) cout << nPoint << " points." << endl;
      }
      else {
        cout << "NPOIN improperly specified!!" << endl;
#ifndef HAVE_MPI
        exit(EXIT_FAILURE);
#else
        MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Finalize();
#endif
      }
      
      node = new CPoint*[nPoint];
      while (iPoint < nPoint) {
        getline(mesh_file,text_line);
        istringstream point_line(text_line);
        switch(nDim) {
          case 2:
            GlobalIndex = iPoint;
#ifndef HAVE_MPI
            point_line >> Coord_2D[0]; point_line >> Coord_2D[1];
#else
            point_line >> Coord_2D[0]; point_line >> Coord_2D[1]; point_line >> LocalIndex;
            if (size > SINGLE_NODE) point_line >> GlobalIndex;
#endif
            node[iPoint] = new CPoint(Coord_2D[0], Coord_2D[1], GlobalIndex, config);
            iPoint++; break;
          case 3:
            GlobalIndex = iPoint;
#ifndef HAVE_MPI
            point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2];
#else
            point_line >> Coord_3D[0]; point_line >> Coord_3D[1]; point_line >> Coord_3D[2]; point_line >> LocalIndex;
            if (size > SINGLE_NODE) point_line >> GlobalIndex;
#endif
            node[iPoint] = new CPoint(Coord_3D[0], Coord_3D[1], Coord_3D[2], GlobalIndex, config);
            iPoint++; break;
        }
      }
    }
    
    /*--- Read number of markers ---*/
    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nMarker = atoi(text_line.c_str());
      if (size == 1) cout << nMarker << " surface markers." << endl;
      config->SetnMarker_All(nMarker);
      bound = new CPrimalGrid**[nMarker];
      nElem_Bound = new unsigned long [nMarker];
      Tag_to_Marker = new string [MAX_NUMBER_MARKER];
      
      for (iMarker = 0 ; iMarker < nMarker; iMarker++) {
        getline (mesh_file,text_line);
        text_line.erase (0,11);
        string::size_type position;
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 );
          if(position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 );
          if(position != string::npos) text_line.erase (position,1);
        }
        Marker_Tag = text_line.c_str();
        
        /*--- Physical boundaries definition ---*/
        if (Marker_Tag != "SEND_RECEIVE") {
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          if (size == 1)
            cout << nElem_Bound[iMarker]  << " boundary elements in index "<< iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;
          bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          
          nelem_edge = 0; nelem_triangle = 0; nelem_quad = 0; ielem = 0;
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type;
            switch(VTK_Type) {
              case LINE:
                bound_line >> vnodes_edge[0]; bound_line >> vnodes_edge[1];
                bound[iMarker][ielem] = new CLine(vnodes_edge[0],vnodes_edge[1],2);
                ielem++; nelem_edge++; break;
              case TRIANGLE:
                bound_line >> vnodes_triangle[0]; bound_line >> vnodes_triangle[1]; bound_line >> vnodes_triangle[2];
                bound[iMarker][ielem] = new CTriangle(vnodes_triangle[0],vnodes_triangle[1],vnodes_triangle[2],3);
                ielem++; nelem_triangle++; break;
              case RECTANGLE:
                bound_line >> vnodes_quad[0]; bound_line >> vnodes_quad[1]; bound_line >> vnodes_quad[2]; bound_line >> vnodes_quad[3];
                bound[iMarker][ielem] = new CRectangle(vnodes_quad[0],vnodes_quad[1],vnodes_quad[2],vnodes_quad[3],3);
                ielem++; nelem_quad++; break;
            }
          }
          
          /*--- Update config information storing the boundary information in the right place ---*/
          Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
          config->SetMarker_All_TagBound(iMarker, Marker_Tag);
          config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
          config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
          config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
          config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
          config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
          config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
          config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
          config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
          config->SetMarker_All_Out_1D(iMarker, config->GetMarker_CfgFile_Out_1D(Marker_Tag));
          config->SetMarker_All_SendRecv(iMarker, NONE);
          
        }
        
        /*--- Send-Receive boundaries definition ---*/
        else {
          unsigned long nelem_vertex = 0, vnodes_vertex;
          unsigned short transform;
          getline (mesh_file,text_line);
          text_line.erase (0,13); nElem_Bound[iMarker] = atoi(text_line.c_str());
          bound[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];
          
          nelem_vertex = 0; ielem = 0;
          getline (mesh_file,text_line); text_line.erase (0,8);
          config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
          config->SetMarker_All_SendRecv(iMarker, atoi(text_line.c_str()));
          
          for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
            getline(mesh_file,text_line);
            istringstream bound_line(text_line);
            bound_line >> VTK_Type; bound_line >> vnodes_vertex; bound_line >> transform;
            
            bound[iMarker][ielem] = new CVertexMPI(vnodes_vertex, nDim);
            bound[iMarker][ielem]->SetRotation_Type(transform);
            ielem++; nelem_vertex++;
            if (config->GetMarker_All_SendRecv(iMarker) < 0)
              node[vnodes_vertex]->SetDomain(false);
            
          }
          
        }
        
      }
    }
    
    /*--- Read periodic transformation info (center, rotation, translation) ---*/
    position = text_line.find ("NPERIODIC=",0);
    if (position != string::npos) {
      unsigned short nPeriodic, iPeriodic, iIndex;
      
      /*--- Set bool signifying that periodic transormations were found ---*/
      found_transform = true;
      
      /*--- Read and store the number of transformations. ---*/
      text_line.erase (0,10); nPeriodic = atoi(text_line.c_str());
      if (rank == MASTER_NODE) {
        if (nPeriodic - 1 != 0)
          cout << nPeriodic - 1 << " periodic transformations." << endl;
      }
      config->SetnPeriodicIndex(nPeriodic);
      
      /*--- Store center, rotation, & translation in that order for each. ---*/
      for (iPeriodic = 0; iPeriodic < nPeriodic; iPeriodic++) {
        getline (mesh_file,text_line);
        position = text_line.find ("PERIODIC_INDEX=",0);
        if (position != string::npos) {
          text_line.erase (0,15); iIndex = atoi(text_line.c_str());
          if (iIndex != iPeriodic) {
            cout << "PERIODIC_INDEX out of order in SU2 file!!" << endl;
#ifndef HAVE_MPI
            exit(EXIT_FAILURE);
#else
            MPI_Abort(MPI_COMM_WORLD,1);
            MPI_Finalize();
#endif
          }
        }
        double* center    = new double[3];
        double* rotation  = new double[3];
        double* translate = new double[3];
        getline (mesh_file,text_line);
        istringstream cent(text_line);
        cent >> center[0]; cent >> center[1]; cent >> center[2];
        config->SetPeriodicCenter(iPeriodic, center);
        getline (mesh_file,text_line);
        istringstream rot(text_line);
        rot >> rotation[0]; rot >> rotation[1]; rot >> rotation[2];
        config->SetPeriodicRotation(iPeriodic, rotation);
        getline (mesh_file,text_line);
        istringstream tran(text_line);
        tran >> translate[0]; tran >> translate[1]; tran >> translate[2];
        config->SetPeriodicTranslate(iPeriodic, translate);
      }
      
    }
    
  }
  
  /*--- If no periodic transormations were found, store default zeros ---*/
  if (!found_transform) {
    unsigned short nPeriodic = 1, iPeriodic = 0;
    config->SetnPeriodicIndex(nPeriodic);
    double* center    = new double[3];
    double* rotation  = new double[3];
    double* translate = new double[3];
    for (unsigned short iDim = 0; iDim < 3; iDim++) {
      center[iDim] = 0.0; rotation[iDim] = 0.0; translate[iDim] = 0.0;
    }
    config->SetPeriodicCenter(iPeriodic, center);
    config->SetPeriodicRotation(iPeriodic, rotation);
    config->SetPeriodicTranslate(iPeriodic, translate);
  }
  
  /*--- Close the input file ---*/
  mesh_file.close();
  
#ifdef HAVE_MPI
  if ((size > SINGLE_NODE) && (rank == MASTER_NODE))
    cout << Global_nElem << " interior elements (incl. halo cells). " << Global_nPoint << " points (incl. ghost points) " << endl;
#endif
  
  /*--- Loop over the surface element to set the boundaries ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++)
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);
        
        if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL &&
            config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX &&
            config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)
          node[Point_Surface]->SetSolidBoundary(true);
      }
  
}


CBoundaryGeometry::~CBoundaryGeometry(void) {
  
}

void CBoundaryGeometry::SetVertex(void) {
  unsigned long  iPoint, iVertex, iElem;
  unsigned short iMarker, iNode;
  
  /*--- Initialize the Vertex vector for each node of the grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create and compute the vector with the number of vertex per marker ---*/
  nVertex = new unsigned long [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    nVertex[iMarker] = 0;
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if (node[iPoint]->GetVertex(iMarker) == -1) {
          iVertex = nVertex[iMarker];
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
  
  /*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1,iMarker);
  
  /*--- Create the bound vertex structure, note that the order
   is the same as in the input file, this is important for Send/Receive part ---*/
  vertex = new CVertex**[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [nVertex[iMarker]];
    nVertex[iMarker] = 0;
    /*--- Initialize the number of Bound Vertex for each Marker ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        /*--- Set the vertex in the node information ---*/
        if (node[iPoint]->GetVertex(iMarker) == -1) {
          iVertex = nVertex[iMarker];
          vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);
          node[iPoint]->SetVertex(nVertex[iMarker],iMarker);
          nVertex[iMarker]++;
        }
      }
  }
}

void CBoundaryGeometry::SetPoint_Connectivity(void) {
  
  unsigned short Node_Neighbor, iNode, iNeighbor, iMarker, jMarker;
  unsigned long jElem, Point_Neighbor, iPoint, iElem, iVertex, kElem;
  
  /*--- Loop over all of the markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Loop over all of the elements on this marker ---*/
    
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      
      /*--- Loop over all the nodes of an element ---*/
      
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        
        /*--- Store this element as a neighbor of this point ---*/
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        node[iPoint]->SetElem(iElem);
        
      }
    }
  }

  /*--- Loop over all of the markers ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    /*--- Loop over all the vertices on this boundary marker ---*/
    
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      
      iPoint = vertex[iMarker][iVertex]->GetNode();
      
      /*--- Loop over all elements shared by the point ---*/
      
      for(iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {
        
        jElem = node[iPoint]->GetElem(iElem);
        
        for (jMarker = 0; jMarker < nMarker; jMarker++) {
          for (kElem = 0; kElem < nElem_Bound[jMarker]; kElem++) {
            
            /*--- We've found this boundary element (could be in another marker) ---*/
            if (kElem == jElem) {
              
              /*--- If we find the point iPoint in the surrounding element ---*/
              
              for(iNode = 0; iNode < bound[jMarker][kElem]->GetnNodes(); iNode++)
                
                if (bound[jMarker][kElem]->GetNode(iNode) == iPoint)
                  
                /*--- Localize the local index of the neighbor of iPoint in the element ---*/
                  
                  for(iNeighbor = 0; iNeighbor < bound[jMarker][kElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
                    Node_Neighbor = bound[jMarker][kElem]->GetNeighbor_Nodes(iNode,iNeighbor);
                    Point_Neighbor = bound[jMarker][kElem]->GetNode(Node_Neighbor);
                    
                    /*--- Store the point into the point ---*/
                    
                    node[iPoint]->SetPoint(Point_Neighbor);
                  }
            }
          }
        }
      }
    }
  }
  
  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      iPoint = vertex[iMarker][iVertex]->GetNode();
      node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());
    }
  }
  
}

void CBoundaryGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
  unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
  unsigned long Neighbor_Point, iVertex, iPoint, iElem;
  double **Coord;
  unsigned short nNode;
  unsigned long elem_poin;
  
  /*--- Center of gravity for face elements ---*/
  for(iMarker = 0; iMarker < nMarker; iMarker++)
    for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
      nNode = bound[iMarker][iElem]->GetnNodes();
      Coord = new double* [nNode];
      /*--- Store the coordinates for all the element nodes ---*/
      for (iNode = 0; iNode < nNode; iNode++) {
        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
        Coord[iNode] = new double [nDim];
        for (iDim = 0; iDim < nDim; iDim++)
          Coord[iNode][iDim] = node[elem_poin]->GetCoord(iDim);
      }
      /*--- Compute the element CG coordinates ---*/
      bound[iMarker][iElem]->SetCG(Coord);
      for (iNode=0; iNode < nNode; iNode++)
        if (Coord[iNode] != NULL) delete[] Coord[iNode];
      if (Coord != NULL) delete[] Coord;
    }
  
  double *Coord_Edge_CG = new double [nDim];
  double *Coord_Elem_CG = new double [nDim];
  double *Coord_Vertex = new double [nDim];
  
  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++)
  /*--- Loop over all the boundary elements ---*/
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
    /*--- Loop over all the nodes of the boundary ---*/
      for(iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        iVertex = node[iPoint]->GetVertex(iMarker);
        /*--- Loop over the neighbor nodes, there is a face for each one ---*/
        for(iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
          Neighbor_Node  = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode,iNeighbor_Nodes);
          Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);
          /*--- Shared edge by the Neighbor Point and the point ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = 0.5*(node[iPoint]->GetCoord(iDim) + node[Neighbor_Point]->GetCoord(iDim));
            Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
            Coord_Vertex[iDim]  = node[iPoint]->GetCoord(iDim);
          }
          
          vertex[iMarker][iVertex]->SetCoord(Coord_Vertex);
          
          switch (nDim) {
            case 2:
              /*--- Store the 2D face (ojo hay cambio de sentido para ajustarse al sentido del contorno de nodo 0 al 1) ---*/
              if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex);
              if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG);
              break;
            case 3:
              /*--- Store the 3D face (ojo hay cambio de sentido para ajustarse al sentido del contorno de nodo 0 al 1) ---*/
              if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex);
              if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex);
          }
        }
      }
  
  delete[] Coord_Edge_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_Vertex;
}

void CBoundaryGeometry::SetBoundSensitivity(CConfig *config) {
  unsigned short iMarker, icommas;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  double Sensitivity;
  bool *PointInDomain;
  int rank = MASTER_NODE;
  int size = SINGLE_NODE;
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  
  nPointLocal = nPoint;
#ifdef HAVE_MPI
  MPI_Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nPointGlobal = nPointLocal;
#endif
  
  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];
  
  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    if (config->GetMarker_All_DV(iMarker) == YES)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        
        /*--- The sensitivity file uses the global numbering ---*/
#ifndef HAVE_MPI
        iPoint = vertex[iMarker][iVertex]->GetNode();
#else
        iPoint = node[vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();
#endif
        if (vertex[iMarker][iVertex]->GetNode() < GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          vertex[iMarker][iVertex]->SetAuxVar(0.0);
        }
      }
  
  /*--- Time-average any unsteady surface sensitivities ---*/
  unsigned long iExtIter, nExtIter;
  double delta_T, total_T;
  if (config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) {
    nExtIter = config->GetUnst_AdjointIter();
    delta_T  = config->GetDelta_UnstTimeND();
    total_T  = (double)nExtIter*delta_T;
  } else if (config->GetUnsteady_Simulation() == TIME_SPECTRAL) {
    
    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/
    double period = config->GetTimeSpectral_Period();
    nExtIter  = config->GetnTimeInstances();
    delta_T   = period/(double)nExtIter;
    total_T   = period;
    
  } else {
    nExtIter = 1;
    delta_T  = 1.0;
    total_T  = 1.0;
  }
  
  for (iExtIter = 0; iExtIter < nExtIter; iExtIter++) {
    
    /*--- Prepare to read surface sensitivity files (CSV) ---*/
    string text_line;
    ifstream Surface_file;
    char buffer[50];
    char cstr[MAX_STRING_SIZE];
    string surfadj_filename = config->GetSurfAdjCoeff_FileName();
    
    /*--- Remove the domain number from the surface csv filename ---*/
    if (size > SINGLE_NODE) {
      if ((rank+1 >= 0) && (rank+1 < 10)) surfadj_filename.erase (surfadj_filename.end()-2, surfadj_filename.end());
      if ((rank+1 >= 10) && (rank+1 < 100)) surfadj_filename.erase (surfadj_filename.end()-3, surfadj_filename.end());
      if ((rank+1 >= 100) && (rank+1 < 1000)) surfadj_filename.erase (surfadj_filename.end()-4, surfadj_filename.end());
      if ((rank+1 >= 1000) && (rank+1 < 10000)) surfadj_filename.erase (surfadj_filename.end()-5, surfadj_filename.end());
    }
    strcpy (cstr, surfadj_filename.c_str());
    
    /*--- Write file name with extension if unsteady or steady ---*/
    if ((config->GetUnsteady_Simulation() && config->GetWrt_Unsteady()) ||
        (config->GetUnsteady_Simulation() == TIME_SPECTRAL)) {
      if ((int(iExtIter) >= 0)    && (int(iExtIter) < 10))    sprintf (buffer, "_0000%d.csv", int(iExtIter));
      if ((int(iExtIter) >= 10)   && (int(iExtIter) < 100))   sprintf (buffer, "_000%d.csv",  int(iExtIter));
      if ((int(iExtIter) >= 100)  && (int(iExtIter) < 1000))  sprintf (buffer, "_00%d.csv",   int(iExtIter));
      if ((int(iExtIter) >= 1000) && (int(iExtIter) < 10000)) sprintf (buffer, "_0%d.csv",    int(iExtIter));
      if  (int(iExtIter) >= 10000) sprintf (buffer, "_%d.csv", int(iExtIter));
    }
    else
      sprintf (buffer, ".csv");
    
    strcat (cstr, buffer);
    
    /*--- Read the sensitivity file ---*/
    string::size_type position;
    
    Surface_file.open(cstr, ios::in);
    getline(Surface_file,text_line);
    
    while (getline(Surface_file,text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if(position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      point_line >> iPoint >> Sensitivity;
      
      if (PointInDomain[iPoint]) {
        
        /*--- Find the vertex for the Point and Marker ---*/
        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];
        
        /*--- Increment the auxiliary variable with the contribution of
         this unsteady timestep. For steady problems, this reduces to
         a single sensitivity value multiplied by 1.0. ---*/
        vertex[iMarker][iVertex]->AddAuxVar(Sensitivity*(delta_T/total_T));
      }
      
    }
    Surface_file.close();
  }
  
  delete[] Point2Vertex;
}

double CBoundaryGeometry::Compute_MaxThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, CConfig *config, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
  unsigned long iVertex, jVertex, n, Trailing_Point, Leading_Point;
  double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, zp1, zpn, MaxThickness_Value = 0, MaxThickness_Location, Thickness, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, MaxDistance, Distance, AoA;
  vector<double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  
  /*--- Translate to the origin ---*/
  
  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }
  
  /*--- Rotate the airfoil ---*/
  
  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
    
  }
  
  /*--- Identify upper and lower side, and store the value of the normal --*/
  
  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
    
    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
    
    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
    
    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
    
    unsigned short index = 2;
    if ((config->GetAxis_Orientation() == Z_AXIS) && (nDim == 3)) index = 0;
    
    if (Normal[index] >= 0.0) {
      Xcoord.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord.push_back(Zcoord_Airfoil_[iVertex]);
    }
    
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) {
      if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
        auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex+1]; Ycoord[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex+1]; Zcoord[jVertex+1] = auxZCoord;
      }
    }
  }
  
  n = Xcoord.size();
  zp1 = (Zcoord[1]-Zcoord[0])/(Xcoord[1]-Xcoord[0]);
  zpn = (Zcoord[n-1]-Zcoord[n-2])/(Xcoord[n-1]-Xcoord[n-2]);
  Z2coord.resize(n+1);
  SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);
  
  /*--- Compute the thickness (we add a fabs because we can not guarantee the
   right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
  
  MaxThickness_Value = 0.0; MaxThickness_Location = 0.0;
  for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    if (Zcoord_Normal[iVertex] < 0.0) {
      Thickness = fabs(Zcoord_Airfoil_[iVertex] - GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil_[iVertex]));
      if (Thickness > MaxThickness_Value) { MaxThickness_Value = Thickness; MaxThickness_Location = Xcoord_Airfoil_[iVertex]; }
    }
  }
  
  return MaxThickness_Value;
  
}

double CBoundaryGeometry::Compute_AoA(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
  unsigned long iVertex, Trailing_Point, Leading_Point;
  double MaxDistance, Distance, AoA = 0.0;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  
  return AoA;
  
}

double CBoundaryGeometry::Compute_Chord(double *Plane_P0, double *Plane_Normal, unsigned short iSection, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
  unsigned long iVertex, Trailing_Point, Leading_Point;
  double MaxDistance, Distance, Chord = 0.0;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  MaxDistance = 0.0; Trailing_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  Chord = MaxDistance;
  
  return Chord;
  
}

double CBoundaryGeometry::Compute_Thickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, double Location, CConfig *config, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
  unsigned long iVertex, jVertex, n_Upper, n_Lower, Trailing_Point, Leading_Point;
  double Thickness_Location, Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Thickness_Value = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, zp1, zpn, Chord, MaxDistance, Distance, AoA;
  vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Z2coord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  Chord = MaxDistance;
  
  /*--- Translate to the origin ---*/
  
  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }
  
  /*--- Rotate the airfoil ---*/
  
  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
  }
  
  /*--- Identify upper and lower side, and store the value of the normal --*/
  
  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
    
    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
    
    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
    
    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
    
    unsigned short index = 2;
    if ((config->GetAxis_Orientation() == Z_AXIS) && (nDim == 3)) index = 0;
    
    if (Normal[index] >= 0.0) {
      Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
    }
    else {
      Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
    }
    
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
        auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
      }
    }
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
        auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
      }
    }
  }
  
  n_Upper = Xcoord_Upper.size();
  zp1 = (Zcoord_Upper[1]-Zcoord_Upper[0])/(Xcoord_Upper[1]-Xcoord_Upper[0]);
  zpn = (Zcoord_Upper[n_Upper-1]-Zcoord_Upper[n_Upper-2])/(Xcoord_Upper[n_Upper-1]-Xcoord_Upper[n_Upper-2]);
  Z2coord_Upper.resize(n_Upper+1);
  SetSpline(Xcoord_Upper, Zcoord_Upper, n_Upper, zp1, zpn, Z2coord_Upper);
  
  n_Lower = Xcoord_Lower.size();
  zp1 = (Zcoord_Lower[1]-Zcoord_Lower[0])/(Xcoord_Lower[1]-Xcoord_Lower[0]);
  zpn = (Zcoord_Lower[n_Lower-1]-Zcoord_Lower[n_Lower-2])/(Xcoord_Lower[n_Lower-1]-Xcoord_Lower[n_Lower-2]);
  Z2coord_Lower.resize(n_Lower+1);
  SetSpline(Xcoord_Lower, Zcoord_Lower, n_Lower, zp1, zpn, Z2coord_Lower);
  
  /*--- Compute the thickness (we add a fabs because we can not guarantee the
   right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/
  
  Thickness_Location = - Chord*(1.0-Location);
  
  Thickness_Value = fabs(GetSpline(Xcoord_Upper, Zcoord_Upper, Z2coord_Upper, n_Upper, Thickness_Location) - GetSpline(Xcoord_Lower, Zcoord_Lower, Z2coord_Lower, n_Lower, Thickness_Location));
  
  return Thickness_Value;
  
}

double CBoundaryGeometry::Compute_Area(double *Plane_P0, double *Plane_Normal, unsigned short iSection, CConfig *config, vector<double> &Xcoord_Airfoil, vector<double> &Ycoord_Airfoil, vector<double> &Zcoord_Airfoil, bool original_surface) {
  unsigned long iVertex, jVertex;
  double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Area_Value = 0.0, Area_Value_Upper = 0.0, Area_Value_Lower = 0.0, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue;
  vector<double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  unsigned long Trailing_Point, Leading_Point;
  double MaxDistance, Distance, AoA;
  
  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  
  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }
  
  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  
  /*--- Translate to the origin ---*/
  
  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }
  
  /*--- Rotate the airfoil ---*/
  
  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);
  
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
    
  }
  
  /*--- Identify upper and lower side, and store the value of the normal --*/
  
  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;
    
    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;
    
    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];
    
    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);
    
    unsigned short index = 2;
    if ((config->GetAxis_Orientation() == Z_AXIS) && (nDim == 3)) index = 0;
    
    if (Normal[index] >= 0.0) {
      Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
    }
    else {
      Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
    }
    
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
        auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
      }
    }
  }
  
  /*--- Order the arrays using the X component ---*/
  
  for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
        auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
      }
    }
  }
  
  /*--- Compute total area ---*/
  
  Area_Value = 0.0;
  Area_Value_Upper = 0.0;
  Area_Value_Lower = 0.0;

  for (iVertex = 0; iVertex < Xcoord_Upper.size()-1; iVertex++)
    Area_Value_Upper += (Xcoord_Upper[iVertex+1] - Xcoord_Upper[iVertex]) * 0.5*(Zcoord_Upper[iVertex+1] + Zcoord_Upper[iVertex]);
  for (iVertex = 0; iVertex < Xcoord_Lower.size()-1; iVertex++)
    Area_Value_Lower += (Xcoord_Lower[iVertex+1] - Xcoord_Lower[iVertex]) * 0.5*(Zcoord_Lower[iVertex+1] + Zcoord_Lower[iVertex]);
  
  Area_Value = fabs(Area_Value_Upper - Area_Value_Lower);
  return Area_Value;
  
}


double CBoundaryGeometry::Compute_Volume(CConfig *config, bool original_surface) {
  
  int rank = MASTER_NODE;
  
  /*--- MPI initialization ---*/
  
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  unsigned short iPlane, nPlane;
  double Volume = 0.0, MinPlane, MaxPlane, MinXCoord, MaxXCoord, dPlane,
  **Plane_P0, **Plane_Normal, *Area;
  vector<double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  
  /*--- Make a large number of section cuts for approximating volume ---*/
  
  nPlane = config->GetnVolSections();
  
  /*--- Allocate memory for the section cutting ---*/
  Area         = new double [nPlane];
  Plane_P0     = new double*[nPlane];
  Plane_Normal = new double*[nPlane];
  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    Plane_P0[iPlane] = new double[3];
    Plane_Normal[iPlane] = new double[3];
  }
  
  MinPlane = config->GetSection_Location(0); MaxPlane = config->GetSection_Location(1);
  MinXCoord = -1E6; MaxXCoord = 1E6;
  dPlane = fabs((MaxPlane - MinPlane)/double(nPlane-1));
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
    Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;
    Plane_Normal[iPlane][config->GetAxis_Orientation()] = 1.0;
    Plane_P0[iPlane][config->GetAxis_Orientation()] = MinPlane + iPlane*dPlane;
  }
  
  /*--- Allocate some vectors for storing airfoil coordinates ---*/
  
  Xcoord_Airfoil   = new vector<double>[nPlane];
  Ycoord_Airfoil   = new vector<double>[nPlane];
  Zcoord_Airfoil   = new vector<double>[nPlane];
  Variable_Airfoil = new vector<double>[nPlane];
  
  /*--- Create the section slices through the geometry ---*/
  
  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane], iPlane,
                           MinXCoord, MaxXCoord, NULL, Xcoord_Airfoil[iPlane],
                           Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                           Variable_Airfoil[iPlane], original_surface, config);
  }
  
  /*--- Compute the area at each section ---*/
  
  if (rank == MASTER_NODE) {
    
    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      Area[iPlane] = 0.0;
      if (Xcoord_Airfoil[iPlane].size() != 0) {
        Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane],
                                    iPlane, config, Xcoord_Airfoil[iPlane],
                                    Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                    original_surface);
      }
    }
    
    /*--- Compute the volume using a composite Simpson's rule ---*/
    
    Volume = 0.0;
    for (iPlane = 0; iPlane < nPlane-2; iPlane+=2) {
      if (Xcoord_Airfoil[iPlane].size() != 0) {
        Volume += (1.0/3.0)*dPlane*(Area[iPlane] + 4.0*Area[iPlane+1] + Area[iPlane+2]);
      }
    }
    
  }
  
  /*--- Free memory for the section cuts ---*/
  
  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  delete [] Variable_Airfoil;
  
  for(iPlane = 0; iPlane < nPlane; iPlane++ ) {
    delete Plane_P0[iPlane];
    delete Plane_Normal[iPlane];
  }
  delete [] Plane_P0;
  delete [] Plane_Normal;
  delete [] Area;
  
  /*--- Return the volume and exit ---*/
  
  return Volume;
}

void CBoundaryGeometry::SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {
  
  ofstream Tecplot_File;
  unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
  unsigned short Coord_i, iMarker;
  
  /*--- It is important to do a renumbering to don't add points
   that do not belong to the surfaces ---*/
  
  PointSurface = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (node[iPoint]->GetBoundary()) {
      PointSurface[iPoint] = nPointSurface;
      nPointSurface++;
    }
  
  /*--- Compute the total number of elements ---*/
  
  Total_nElem_Bound = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      Total_nElem_Bound += nElem_Bound[iMarker];
    }
  }
  
  /*--- Open the tecplot file and write the header ---*/
  
  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the surface grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\",\"Curvature\"  " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\",\"Curvature\"  " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);
  
  if (Total_nElem_Bound != 0) {
    
    /*--- Write the header of the file ---*/
    
    Tecplot_File << "ZONE T= ";
    if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
    else Tecplot_File << "\"Deformed grid\", C=RED, ";
    Tecplot_File << "NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING= POINT";
    if (nDim == 2) Tecplot_File << ", ZONETYPE= FELINESEG"<< endl;
    if (nDim == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
    
    /*--- Only write the coordinates of the points that are on the surfaces ---*/
    
    //    if (nDim == 3) {
    //      for(iPoint = 0; iPoint < nPoint; iPoint++)
    //        if (node[iPoint]->GetBoundary()) {
    //          for(Coord_i = 0; Coord_i < nDim-1; Coord_i++)
    //            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
    //          Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
    //        }
    //    }
    //    else {
    //      for(iPoint = 0; iPoint < nPoint; iPoint++)
    //        if (node[iPoint]->GetBoundary()){
    //          for(Coord_i = 0; Coord_i < nDim; Coord_i++)
    //            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
    //          Tecplot_File << "\n";
    //        }
    //    }
    
    for(iPoint = 0; iPoint < nPoint; iPoint++)
      if (node[iPoint]->GetBoundary()) {
        for(Coord_i = 0; Coord_i < nDim; Coord_i++)
          Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
        Tecplot_File << node[iPoint]->GetCurvature() << "\n";
      }
    
    
    /*--- Write the cells using the new numbering ---*/
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for(iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
          if (nDim == 2) {
            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
          }
          if (nDim == 3) {
            if (bound[iMarker][iElem]->GetnNodes() == 3) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
            }
            if (bound[iMarker][iElem]->GetnNodes() == 4) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
            }
          }
        }
  }
  else {
    
    /*--- No elements in the surface ---*/
    
    if (nDim == 2) {
      Tecplot_File << "ZONE NODES= 1, C=BLACK, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      Tecplot_File << "0.0 0.0"<< endl;
      Tecplot_File << "1 1"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "ZONE NODES= 1, C=RED, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      Tecplot_File << "0.0 0.0 0.0"<< endl;
      Tecplot_File << "1 1 1 1"<< endl;
    }
  }
  
  /*--- Dealocate memory and close the file ---*/
  
  delete[] PointSurface;
  Tecplot_File.close();
  
}

CPeriodicGeometry::CPeriodicGeometry(CGeometry *geometry, CConfig *config) {
  unsigned long nElem_new, nPoint_new, jPoint, iPoint, iElem, jElem, iVertex,
  nelem_triangle = 0, nelem_quad = 0, nelem_tetra = 0, nelem_hexa = 0, nelem_wedge = 0,
  nelem_pyramid = 0, iIndex, newElementsBound = 0;
  unsigned short  iMarker, nPeriodic = 0, iPeriodic;
  double *center, *angles, rotMatrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}},
  translation[3], *trans, theta, phi, psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi,
  dx, dy, dz, rotCoord[3], *Coord_i;
  
  /*--- It only create the mirror structure for the second boundary ---*/
  bool CreateMirror[10];
  CreateMirror[1] = false;
  CreateMirror[2] = true;
  
  /*--- Compute the number of periodic bc on the geometry ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      nPeriodic++;
  
  /*--- Write the number of dimensions of the problem ---*/
  nDim = geometry->GetnDim();
  
  /*--- Copy the new boundary element information from the geometry class.
   Be careful, as these are pointers to vectors/objects. ---*/
  nNewElem_BoundPer = geometry->nNewElem_Bound;
  newBoundPer       = geometry->newBound;
  
  /*--- Count the number of new boundary elements. ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    newElementsBound += nNewElem_BoundPer[iMarker];
  
  /*--- Loop over the original grid to perform the dimensionalizaton of the new vectors ---*/
  nElem_new = 0; nPoint_new = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      nElem_new += geometry->PeriodicElem[iPeriodic].size();
      nPoint_new += geometry->PeriodicPoint[iPeriodic][0].size();
    }
  }
  
  cout << "Number of new points: " << nPoint_new << "." << endl;
  cout << "Number of new interior elements: " << nElem_new << "." << endl;
  cout << "Number of new boundary elements added to preexisting markers: " << newElementsBound << "." << endl;
  
  /*--- Create a copy of the original grid ---*/
  elem = new CPrimalGrid*[geometry->GetnElem() + nElem_new];
  for (iElem = 0; iElem < geometry->GetnElem(); iElem ++) {
    switch(geometry->elem[iElem]->GetVTK_Type()) {
      case TRIANGLE:
        elem[iElem] = new CTriangle(geometry->elem[iElem]->GetNode(0),
                                    geometry->elem[iElem]->GetNode(1),
                                    geometry->elem[iElem]->GetNode(2), 2);
        nelem_triangle++;
        break;
        
      case RECTANGLE:
        elem[iElem] = new CRectangle(geometry->elem[iElem]->GetNode(0),
                                     geometry->elem[iElem]->GetNode(1),
                                     geometry->elem[iElem]->GetNode(2),
                                     geometry->elem[iElem]->GetNode(3), 2);
        nelem_quad++;
        break;
        
      case TETRAHEDRON:
        elem[iElem] = new CTetrahedron(geometry->elem[iElem]->GetNode(0),
                                       geometry->elem[iElem]->GetNode(1),
                                       geometry->elem[iElem]->GetNode(2),
                                       geometry->elem[iElem]->GetNode(3));
        nelem_tetra++;
        break;
        
      case HEXAHEDRON:
        elem[iElem] = new CHexahedron(geometry->elem[iElem]->GetNode(0),
                                      geometry->elem[iElem]->GetNode(1),
                                      geometry->elem[iElem]->GetNode(2),
                                      geometry->elem[iElem]->GetNode(3),
                                      geometry->elem[iElem]->GetNode(4),
                                      geometry->elem[iElem]->GetNode(5),
                                      geometry->elem[iElem]->GetNode(6),
                                      geometry->elem[iElem]->GetNode(7));
        nelem_hexa++;
        break;
        
      case WEDGE:
        elem[iElem] = new CWedge(geometry->elem[iElem]->GetNode(0),
                                 geometry->elem[iElem]->GetNode(1),
                                 geometry->elem[iElem]->GetNode(2),
                                 geometry->elem[iElem]->GetNode(3),
                                 geometry->elem[iElem]->GetNode(4),
                                 geometry->elem[iElem]->GetNode(5));
        nelem_wedge++;
        break;
        
      case PYRAMID:
        elem[iElem] = new CPyramid(geometry->elem[iElem]->GetNode(0),
                                   geometry->elem[iElem]->GetNode(1),
                                   geometry->elem[iElem]->GetNode(2),
                                   geometry->elem[iElem]->GetNode(3),
                                   geometry->elem[iElem]->GetNode(4));
        nelem_pyramid++;
        break;
        
    }
  }
  
  /*--- Create a list with all the points and the new index ---*/
  unsigned long *Index = new unsigned long [geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) Index[iPoint] = 0;
  
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        iPoint =  geometry->PeriodicPoint[iPeriodic][0][iIndex];
        Index[iPoint] = geometry->PeriodicPoint[iPeriodic][1][iIndex];
      }
    }
  }
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
        Index[iPoint] = jPoint;
      }
  
  /*--- Add the new elements due to the periodic boundary condtion ---*/
  iElem = geometry->GetnElem();
  
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicElem[iPeriodic].size(); iIndex++) {
        jElem = geometry->PeriodicElem[iPeriodic][iIndex];
        
        switch(geometry->elem[jElem]->GetVTK_Type()) {
          case TRIANGLE:
            elem[iElem] = new CTriangle(Index[geometry->elem[jElem]->GetNode(0)],
                                        Index[geometry->elem[jElem]->GetNode(1)],
                                        Index[geometry->elem[jElem]->GetNode(2)], 2);
            iElem++; nelem_triangle++;
            break;
            
          case RECTANGLE:
            elem[iElem] = new CRectangle(Index[geometry->elem[jElem]->GetNode(0)],
                                         Index[geometry->elem[jElem]->GetNode(1)],
                                         Index[geometry->elem[jElem]->GetNode(2)],
                                         Index[geometry->elem[jElem]->GetNode(3)], 2);
            iElem++; nelem_quad++;
            break;
            
          case TETRAHEDRON:
            elem[iElem] = new CTetrahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                           Index[geometry->elem[jElem]->GetNode(1)],
                                           Index[geometry->elem[jElem]->GetNode(2)],
                                           Index[geometry->elem[jElem]->GetNode(3)]);
            iElem++; nelem_tetra++;
            break;
            
          case HEXAHEDRON:
            elem[iElem] = new CHexahedron(Index[geometry->elem[jElem]->GetNode(0)],
                                          Index[geometry->elem[jElem]->GetNode(1)],
                                          Index[geometry->elem[jElem]->GetNode(2)],
                                          Index[geometry->elem[jElem]->GetNode(3)],
                                          Index[geometry->elem[jElem]->GetNode(4)],
                                          Index[geometry->elem[jElem]->GetNode(5)],
                                          Index[geometry->elem[jElem]->GetNode(6)],
                                          Index[geometry->elem[jElem]->GetNode(7)]);
            iElem++; nelem_hexa++;
            break;
            
          case WEDGE:
            elem[iElem] = new CWedge(Index[geometry->elem[jElem]->GetNode(0)],
                                     Index[geometry->elem[jElem]->GetNode(1)],
                                     Index[geometry->elem[jElem]->GetNode(2)],
                                     Index[geometry->elem[jElem]->GetNode(3)],
                                     Index[geometry->elem[jElem]->GetNode(4)],
                                     Index[geometry->elem[jElem]->GetNode(5)]);
            iElem++; nelem_wedge++;
            break;
            
          case PYRAMID:
            elem[iElem] = new CPyramid(Index[geometry->elem[jElem]->GetNode(0)],
                                       Index[geometry->elem[jElem]->GetNode(1)],
                                       Index[geometry->elem[jElem]->GetNode(2)],
                                       Index[geometry->elem[jElem]->GetNode(3)],
                                       Index[geometry->elem[jElem]->GetNode(4)]);
            iElem++; nelem_pyramid++;
            break;
            
        }
      }
    }
  }
  
  nElem = geometry->GetnElem() + nElem_new;
  
  /*--- Add the old points ---*/
  node = new CPoint*[geometry->GetnPoint() + nPoint_new];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    if (geometry->GetnDim() == 2)
      node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                geometry->node[iPoint]->GetCoord(1), iPoint, config);
    if (geometry->GetnDim() == 3)
      node[iPoint] = new CPoint(geometry->node[iPoint]->GetCoord(0),
                                geometry->node[iPoint]->GetCoord(1),
                                geometry->node[iPoint]->GetCoord(2), iPoint, config);
  }
  
  /*--- Add the new points due to the periodic boundary condtion (only in the mirror part) ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (CreateMirror[iPeriodic]) {
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        
        /*--- From iPeriodic obtain the iMarker ---*/
        for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
          if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
        
        /*--- Retrieve the supplied periodic information. ---*/
        center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
        angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
        trans  = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
        
        /*--- Store center - trans as it is constant and will be added on.
         Note the subtraction, as this is the inverse translation. ---*/
        translation[0] = center[0] - trans[0];
        translation[1] = center[1] - trans[1];
        translation[2] = center[2] - trans[2];
        
        /*--- Store angles separately for clarity. Compute sines/cosines.
         Note the negative sign, as this is the inverse rotation. ---*/
        theta = -angles[0];
        phi   = -angles[1];
        psi   = -angles[2];
        
        cosTheta = cos(theta);  cosPhi = cos(phi);  cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);  sinPsi = sin(psi);
        
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
        
        /*--- Retrieve node information for this boundary point. ---*/
        iPoint = geometry->PeriodicPoint[iPeriodic][0][iIndex];
        jPoint = geometry->PeriodicPoint[iPeriodic][1][iIndex];
        Coord_i = geometry->node[iPoint]->GetCoord();
        
        /*--- Get the position vector from rot center to point. ---*/
        dx = Coord_i[0] - center[0];
        dy = Coord_i[1] - center[1];
        if (nDim == 3) {
          dz = Coord_i[2] - center[2];
        } else {
          dz = 0.0;
        }
        
        /*--- Compute transformed point coordinates. ---*/
        rotCoord[0] = rotMatrix[0][0]*dx + rotMatrix[0][1]*dy + rotMatrix[0][2]*dz + translation[0];
        rotCoord[1] = rotMatrix[1][0]*dx + rotMatrix[1][1]*dy + rotMatrix[1][2]*dz + translation[1];
        rotCoord[2] = rotMatrix[2][0]*dx + rotMatrix[2][1]*dy + rotMatrix[2][2]*dz + translation[2];
        
        /*--- Save the new points with the new coordinates. ---*/
        if (geometry->GetnDim() == 2)
          node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], jPoint, config);
        if (geometry->GetnDim() == 3)
          node[jPoint] = new CPoint(rotCoord[0], rotCoord[1], rotCoord[2], jPoint, config);
        
      }
    }
  }
  
  nPoint = geometry->GetnPoint() + nPoint_new;
  
  /*--- Add the old boundary, reserving space for two new bc (send/recive periodic bc) ---*/
  nMarker = geometry->GetnMarker() + 2;
  nElem_Bound = new unsigned long [nMarker];
  bound = new CPrimalGrid**[nMarker];
  Tag_to_Marker = new string [MAX_NUMBER_MARKER];
  config->SetnMarker_All(nMarker);
  
  /*--- Copy the olf boundary ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    
    bound[iMarker] = new CPrimalGrid* [geometry->GetnElem_Bound(iMarker)];
    
    for (iVertex = 0; iVertex < geometry->GetnElem_Bound(iMarker); iVertex++) {
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == LINE)
        bound[iMarker][iVertex] = new CLine(geometry->bound[iMarker][iVertex]->GetNode(0),
                                            geometry->bound[iMarker][iVertex]->GetNode(1), 2);
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == TRIANGLE)
        bound[iMarker][iVertex] = new CTriangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                geometry->bound[iMarker][iVertex]->GetNode(1),
                                                geometry->bound[iMarker][iVertex]->GetNode(2), 3);
      if (geometry->bound[iMarker][iVertex]->GetVTK_Type() == RECTANGLE)
        bound[iMarker][iVertex] = new CRectangle(geometry->bound[iMarker][iVertex]->GetNode(0),
                                                 geometry->bound[iMarker][iVertex]->GetNode(1),
                                                 geometry->bound[iMarker][iVertex]->GetNode(2),
                                                 geometry->bound[iMarker][iVertex]->GetNode(3), 3);
    }
    
    nElem_Bound[iMarker] = geometry->GetnElem_Bound(iMarker);
    Tag_to_Marker[iMarker] = geometry->GetMarker_Tag(iMarker);
    
  }
  
  delete[] Index;
  
}

CPeriodicGeometry::~CPeriodicGeometry(void) {
  unsigned long iElem_Bound;
  unsigned short iMarker;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
      if (newBoundPer[iMarker][iElem_Bound] != NULL) delete newBoundPer[iMarker][iElem_Bound];
    }
  }
  if (newBoundPer != NULL) delete[] newBoundPer;
  
  if (nNewElem_BoundPer != NULL) delete[] nNewElem_BoundPer;
  
}

void CPeriodicGeometry::SetPeriodicBoundary(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, iPeriodic, nPeriodic = 0, iMarkerSend, iMarkerReceive;
  unsigned long iVertex, Counter_Send = 0, Counter_Receive = 0, iIndex;
  
  /*--- Compute the number of periodic bc on the geometry ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      nPeriodic++;
  
  /*--- First compute the Send/Receive boundaries, count the number of points ---*/
  Counter_Send = 0; 	Counter_Receive = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
      Counter_Send += geometry->PeriodicPoint[iPeriodic][0].size();
    if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
      Counter_Receive += geometry->PeriodicPoint[iPeriodic][1].size();
  }
  
  /*--- Adimensionalization of the new boundaries ---*/
  iMarkerSend = nMarker - 2; iMarkerReceive = nMarker - 1;
  config->SetMarker_All_SendRecv(iMarkerSend,1);
  config->SetMarker_All_SendRecv(iMarkerReceive,-1);
  nElem_Bound[iMarkerSend] = Counter_Send;
  nElem_Bound[iMarkerReceive] = Counter_Receive;
  bound[iMarkerSend] = new CPrimalGrid* [Counter_Send];
  bound[iMarkerReceive] = new CPrimalGrid* [Counter_Receive];
  
  /*--- First we do the send ---*/
  iVertex = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
    if (geometry->PeriodicPoint[iPeriodic][0].size() != 0)
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][0].size(); iIndex++) {
        bound[iMarkerSend][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][0][iIndex], nDim);
        bound[iMarkerSend][iVertex]->SetRotation_Type(iPeriodic);
        iVertex++;
      }
  
  /*--- Second we do the receive ---*/
  iVertex = 0;
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++)
    if (geometry->PeriodicPoint[iPeriodic][1].size() != 0)
      for (iIndex = 0; iIndex < geometry->PeriodicPoint[iPeriodic][1].size(); iIndex++) {
        bound[iMarkerReceive][iVertex] = new CVertexMPI(geometry->PeriodicPoint[iPeriodic][1][iIndex], nDim);
        bound[iMarkerReceive][iVertex]->SetRotation_Type(iPeriodic);
        iVertex++;
      }
  
}

void CPeriodicGeometry::SetMeshFile(CGeometry *geometry, CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound, GhostPoints;
  unsigned short iMarker, iNodes, iDim;
  unsigned short iMarkerReceive, iPeriodic, nPeriodic = 0;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;
  double *center, *angles, *transl;
  
  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());
  
  /*--- Open .su2 grid file ---*/
  output_file.precision(15);
  output_file.open(cstr, ios::out);
  
  /*--- Ghost points, look at the nodes in the send receive ---*/
  iMarkerReceive = nMarker - 1;
  GhostPoints = nElem_Bound[iMarkerReceive];
  
  /*--- Change the numbering to guarantee that the all the receive
   points are at the end of the file ---*/
  unsigned long OldnPoint = geometry->GetnPoint();
  unsigned long NewSort[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    NewSort[iPoint] = iPoint;
  }
  
  unsigned long Index = OldnPoint-1;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      if (config->GetMarker_All_SendRecv(iMarker) < 0) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (bound[iMarker][iElem_Bound]->GetNode(0) < geometry->GetnPoint()) {
            NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] = Index;
            NewSort[Index] = bound[iMarker][iElem_Bound]->GetNode(0);
            Index--;
          }
        }
      }
    }
  }
  
  
  /*--- Write dimension, number of elements and number of points ---*/
  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << NewSort[elem[iElem]->GetNode(iNodes)];
    output_file << "\t"<<iElem<< endl;
  }
  
  output_file << "NPOIN= " << nPoint << "\t" << nPoint - GhostPoints << endl;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[NewSort[iPoint]]->GetCoord(iDim) ;
    output_file << "\t" << iPoint << endl;
  }
  
  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      
      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker] + nNewElem_BoundPer[iMarker] << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
        for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
          output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
        iNodes = bound[iMarker][iElem_Bound]->GetnNodes()-1;
        output_file << NewSort[bound[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
      }
      
      /*--- Write any new elements at the end of the list. ---*/
      if (nNewElem_BoundPer[iMarker] > 0) {
        for (iElem_Bound = 0; iElem_Bound < nNewElem_BoundPer[iMarker]; iElem_Bound++) {
          output_file << newBoundPer[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1; iNodes++)
            output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << "\t" ;
          iNodes = newBoundPer[iMarker][iElem_Bound]->GetnNodes()-1;
          output_file << NewSort[newBoundPer[iMarker][iElem_Bound]->GetNode(iNodes)] << endl;
        }
      }
      
    }
    
    if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      
      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        NewSort[bound[iMarker][iElem_Bound]->GetNode(0)] << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type()  << endl;
      }
    }
  }
  
  /*--- Compute the number of periodic bc on the geometry ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
      nPeriodic++;
  
  output_file << "NPERIODIC= " << nPeriodic + 1 << endl;
  
  /*--- Periodic 0 correspond with no movement of the surface ---*/
  output_file << "PERIODIC_INDEX= 0" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  output_file << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << "\t" << "0.000000000000000e+00" << endl;
  
  /*--- From iPeriodic obtain the iMarker ---*/
  for (iPeriodic = 1; iPeriodic <= nPeriodic; iPeriodic++) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (iPeriodic == config->GetMarker_All_PerBound(iMarker)) break;
    
    /*--- Retrieve the supplied periodic information. ---*/
    center = config->GetPeriodicRotCenter(config->GetMarker_All_TagBound(iMarker));
    angles = config->GetPeriodicRotAngles(config->GetMarker_All_TagBound(iMarker));
    transl = config->GetPeriodicTranslation(config->GetMarker_All_TagBound(iMarker));
    
    output_file << "PERIODIC_INDEX= " << iPeriodic << endl;
    output_file << center[0] << "\t" << center[1] << "\t" << center[2] << endl;
    output_file << angles[0] << "\t" << angles[1] << "\t" << angles[2] << endl;
    output_file << transl[0] << "\t" << transl[1] << "\t" << transl[2] << endl;
    
  }
  
  
  output_file.close();
  
}

void CPeriodicGeometry::SetTecPlot(char mesh_filename[MAX_STRING_SIZE]) {
  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;
  
  Tecplot_File.open(mesh_filename, ios::out);
  Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << endl;
  
  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  }
  if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }
  
  for(iPoint = 0; iPoint < nPoint; iPoint++) {
    for(iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }
  
  for(iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == RECTANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == WEDGE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }
  
  Tecplot_File.close();
}

CMultiGridQueue::CMultiGridQueue(unsigned long val_npoint) {
  unsigned long iPoint;
  
  nPoint = val_npoint;
  Priority = new short[nPoint];
  RightCV = new bool[nPoint];
  
  QueueCV.resize(1);
  
  /*--- Queue initialization with all the points in the finer grid ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    QueueCV[0].push_back(iPoint);
    Priority[iPoint] = 0;
    RightCV[iPoint] = true;
  }
  
}

CMultiGridQueue::~CMultiGridQueue(void) {
  
  delete[] Priority;
  delete[] RightCV;
  
}

void CMultiGridQueue::AddCV(unsigned long val_new_point, unsigned short val_number_neighbors) {
  
  unsigned short Max_Neighbors = QueueCV.size()-1;
  
  /*--- Basic check ---*/
  if (val_new_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Resize the list ---*/
  if (val_number_neighbors > Max_Neighbors)
    QueueCV.resize(val_number_neighbors+1);
  
  /*--- Find the point in the queue ---*/
  bool InQueue = false;
  if (Priority[val_new_point] == val_number_neighbors) InQueue = true;
  
  if (!InQueue) {
    /*--- Add the control volume, and update the priority list ---*/
    QueueCV[val_number_neighbors].push_back(val_new_point);
    Priority[val_new_point] = val_number_neighbors;
  }
  
}

void CMultiGridQueue::RemoveCV(unsigned long val_remove_point) {
  unsigned short iPoint;
  bool check;
  
  /*--- Basic check ---*/
  if (val_remove_point > nPoint) {
    cout << "The index of the CV is greater than the size of the priority list." << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Find priority of the Control Volume ---*/
  short Number_Neighbors = Priority[val_remove_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_remove_point <<" is not in the priority list. (RemoveCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Find the point in the queue ---*/
  vector<unsigned long>::iterator ItQueue = find(QueueCV[Number_Neighbors].begin(),
                                                 QueueCV[Number_Neighbors].end(),
                                                 val_remove_point);
  if( ItQueue != QueueCV[Number_Neighbors].end() ) QueueCV[Number_Neighbors].erase(ItQueue);
  
  Priority[val_remove_point] = -1;
  
  /*--- Check that the size of the queue is the right one ---*/
  unsigned short Size_QueueCV = 0;
  check = false;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { Size_QueueCV = iPoint; check = true;}
  
  /*--- Resize the queue, if check = false, the queue is empty, at least
   we need one element in the queue ---*/
  if (check) QueueCV.resize(Size_QueueCV+1);
  else QueueCV.resize(1);
  
}

void CMultiGridQueue::MoveCV(unsigned long val_move_point, short val_number_neighbors) {
  unsigned short Priority;
  
  if (val_number_neighbors < 0) {
    val_number_neighbors = 0;
    RightCV[val_move_point] = false;
  }
  else {
    Priority = val_number_neighbors;
    RightCV[val_move_point] = true;
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_move_point);
  
  /*--- Add a new control volume ---*/
  AddCV(val_move_point, val_number_neighbors);
  
}

void CMultiGridQueue::IncrPriorityCV(unsigned long val_incr_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_incr_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_incr_point <<" is not in the priority list. (IncrPriorityCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  /*--- Remove the control volume ---*/
  RemoveCV(val_incr_point);
  
  /*--- Increase the priority ---*/
  AddCV(val_incr_point, Number_Neighbors+1);
  
}

void CMultiGridQueue::RedPriorityCV(unsigned long val_red_point) {
  
  /*--- Find the priority list ---*/
  short Number_Neighbors = Priority[val_red_point];
  if (Number_Neighbors == -1) {
    cout << "The CV "<< val_red_point <<" is not in the priority list. (RedPriorityCV)" << endl;
    exit(EXIT_FAILURE);
  }
  
  if (Number_Neighbors != 0) {
    
    /*--- Remove the control volume ---*/
    RemoveCV(val_red_point);
    
    /*--- Increase the priority ---*/
    AddCV(val_red_point, Number_Neighbors-1);
    
  }
  
}

void CMultiGridQueue::VisualizeQueue(void) {
  unsigned short iPoint;
  unsigned long jPoint;
  
  cout << endl;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++) {
    cout << "Number of neighbors " << iPoint <<": ";
    for (jPoint = 0; jPoint < QueueCV[iPoint].size(); jPoint ++) {
      cout << QueueCV[iPoint][jPoint] << " ";
    }
    cout << endl;
  }
  
}

void CMultiGridQueue::VisualizePriority(void) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++)
    cout << "Control Volume: " << iPoint <<" Priority: " << Priority[iPoint] << endl;
  
}

long CMultiGridQueue::NextCV(void) {
  if (QueueCV.size() != 0) return QueueCV[QueueCV.size()-1][0];
  else return -1;
}

bool CMultiGridQueue::EmptyQueue(void) {
  unsigned short iPoint;
  
  /*--- In case there is only the no agglomerated elements,
   check if they can be agglomerated or we have already finished ---*/
  bool check = true;
  
  if ( QueueCV.size() == 1 ) {
    for (iPoint = 0; iPoint < QueueCV[0].size(); iPoint ++) {
      if (RightCV[QueueCV[0][iPoint]]) { check = false; break; }
    }
  }
  else {
    for (iPoint = 1; iPoint < QueueCV.size(); iPoint ++)
      if (QueueCV[iPoint].size() != 0) { check = false; break;}
  }
  
  return check;
}

unsigned long CMultiGridQueue::TotalCV(void) {
  unsigned short iPoint;
  unsigned long TotalCV;
  
  TotalCV = 0;
  for (iPoint = 0; iPoint < QueueCV.size(); iPoint ++)
    if (QueueCV[iPoint].size() != 0) { TotalCV += QueueCV[iPoint].size(); }
  
  return TotalCV;
}

void CMultiGridQueue::Update(unsigned long iPoint, CGeometry *fine_grid) {
  unsigned short iNode;
  unsigned long jPoint;
  
  RemoveCV(iPoint);
  for (iNode = 0; iNode <	fine_grid->node[iPoint]->GetnPoint(); iNode ++) {
    jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
    if (fine_grid->node[jPoint]->GetAgglomerate() == false)
      IncrPriorityCV(jPoint);
  }
  
}
