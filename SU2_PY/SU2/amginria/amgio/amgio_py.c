#include "amgio.h"
#include "Python.h"


int py_ConvertSU2toInria( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 0; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertSU2SolToGMF (mshopt);
	
	
	return 1;
}


int py_ConvertInriatoSU2( char *MshNam, char *SolNam, char *OutNam ) 
{
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->OutNam,OutNam);
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	mshopt->clean = 0; // remove unconnected vertices
	
	if ( !CheckOptions(mshopt) ) {
		return 0;
	}
	
	return ConvertGMFtoSU2Sol (mshopt);
	
	
	return 1;
}


int py_SplitSolution(char *SolNam, int dim, char *prefix, char *adap_sensor)
{
	
	int SizMsh[GmfMaxSizMsh+1];
	memset(SizMsh,0,sizeof(int)*(GmfMaxSizMsh+1));
	
	Mesh *Msh = AllocMesh(SizMsh);
	
	Msh->NbrVer = GetSU2SolSize(SolNam);
	
	LoadSU2Solution(SolNam, Msh);
	
	Msh->Dim = dim;
	SplitSolution(Msh, prefix, adap_sensor);
	
}





void py_ReadMesh (char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg, PyObject *pyHex, 
PyObject *pyQua, PyObject *pyPyr, PyObject *pyPri, 
 PyObject *pySol, PyObject *pySolHeader,  PyObject *pyMarkers)
{
	int i, j, d;
	
	Options *mshopt = AllocOptions();
	
	strcpy(mshopt->InpNam,MshNam);
	strcpy(mshopt->SolNam,SolNam);
	
	//--- Open mesh/solution file
	
	Mesh *Msh = NULL;
	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);
	
	for (i=1; i<=Msh->NbrVer; i++){
		for (d=0; d<3; d++)
			PyList_Append(pyVer, PyFloat_FromDouble(Msh->Ver[i][d]));
	}
	
	for (i=1; i<=Msh->NbrTri; i++){
		for (j=0; j<4; j++)
			PyList_Append(pyTri, PyInt_FromLong(Msh->Tri[i][j]));
	}
	
	for (i=1; i<=Msh->NbrTet; i++){
		for (j=0; j<5; j++)
			PyList_Append(pyTet, PyInt_FromLong(Msh->Tet[i][j]));
	}
	
	for (i=1; i<=Msh->NbrEfr; i++){
		for (j=0; j<3; j++)
			PyList_Append(pyEdg, PyInt_FromLong(Msh->Efr[i][j]));
	}
	
	for (i=1; i<=Msh->NbrHex; i++){
		for (j=0; j<9; j++)
			PyList_Append(pyHex, PyInt_FromLong(Msh->Hex[i][j]));
	}
	
	for (i=1; i<=Msh->NbrQua; i++){
		for (j=0; j<5; j++)
			PyList_Append(pyQua, PyInt_FromLong(Msh->Qua[i][j]));
	}
	
	for (i=1; i<=Msh->NbrPyr; i++){
		for (j=0; j<6; j++)
			PyList_Append(pyPyr, PyInt_FromLong(Msh->Pyr[i][j]));
	}
	
	for (i=1; i<=Msh->NbrPri; i++){
		for (j=0; j<7; j++)
			PyList_Append(pyPri, PyInt_FromLong(Msh->Pri[i][j]));
	}
	
	//--- First row of Markers contains dimension
	PyList_Append(pyMarkers, PyInt_FromLong(Msh->Dim));
	for (i=1; i<=Msh->NbrMarkers; i++){
		PyList_Append(pyMarkers, PyString_FromString(Msh->Markers[i]));
	}
	
	for (i=0; i<=Msh->SolSiz; i++){
		PyList_Append(pySolHeader, PyString_FromString(Msh->SolTag[i]));
	}
	
	if ( Msh->Sol ) {
		
		//--- Output solution
		int iVer;
		for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
			for (i=0; i<Msh->SolSiz; i++) {
				PyList_Append(pySol, PyFloat_FromDouble(Msh->Sol[iVer*Msh->SolSiz+i]));
			}
		}
		
	}
	
	if ( Msh )
 		FreeMesh(Msh);
	
}


void py_WriteMesh(char *MshNam, char *SolNam, PyObject *pyVer, PyObject *pyTri, PyObject *pyTet, PyObject *pyEdg,  PyObject *pyHex, 
PyObject *pyQua, PyObject *pyPyr, PyObject *pyPri, PyObject *pySol, PyObject *pyMarkers, int Dim)
{
	int i, j;
	Mesh *Msh= NULL;
	int SizMsh[GmfMaxKwd+1];
	
	int is[5], siz, ref, idx;
	double crd[3];
	
	int NbrMarkers = 0;
	
	for (i=0; i<GmfMaxKwd; i++)
		SizMsh[i] = 0;
	
	//--- Get mesh size

	if ( PyList_Check(pyVer) )
		SizMsh[GmfVertices] = PyList_Size(pyVer);
	
	if ( PyList_Check(pyTri) )
		SizMsh[GmfTriangles] = PyList_Size(pyTri);
	
	if ( PyList_Check(pyTet) )
		SizMsh[GmfTetrahedra] = PyList_Size(pyTet);
	
	if ( PyList_Check(pyEdg) )
		SizMsh[GmfEdges] = PyList_Size(pyEdg);

	if ( PyList_Check(pyHex) )
		SizMsh[GmfHexahedra] = PyList_Size(pyHex);
	
	if ( PyList_Check(pyQua) )
		SizMsh[GmfQuadrilaterals] = PyList_Size(pyQua);
	
	if ( PyList_Check(pyPyr) )
		SizMsh[GmfPyramids] = PyList_Size(pyPyr);
	
	if ( PyList_Check(pyPri) )
		SizMsh[GmfPrisms] = PyList_Size(pyPri);

	if ( PyList_Check(pyMarkers) )
		NbrMarkers = PyList_Size(pyMarkers);
	
	//--- Allocate mesh
	
	Msh = AllocMesh(SizMsh);
	
	Msh->Dim = Dim;
	
	//--- Fill mesh
	
	if ( PyList_Check(pyTri) )
  {
			siz = PyList_Size(pyTri);
			
			for (i=0; i<siz/4; i++)
      {
				idx = 4*i;
				
				for (j=0; j<3; j++) {
	       	PyObject *oo = PyList_GetItem(pyTri,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyTri,idx+3);
				ref = (int) PyInt_AS_LONG(oo);
				
				Msh->NbrTri++;
				AddTriangle(Msh,Msh->NbrTri,is,ref);
				
				//printf("-- Add tri %d : %d %d %d (ref %d)\n", Msh->NbrTri, is[0], is[1], is[2], ref);
				//exit(1);
      }
  }
	
	if ( PyList_Check(pyTet) )
  {
			siz = PyList_Size(pyTet);
			
			for (i=0; i<siz/5; i++)
      {
				idx = 5*i;
				
				for (j=0; j<5; j++) {
	       	PyObject *oo = PyList_GetItem(pyTet,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				Msh->NbrTet++;
				AddTetrahedron(Msh,Msh->NbrTet,is,is[4]);
				
      }
  }
		
	if ( PyList_Check(pyEdg) )
  {
			siz = PyList_Size(pyEdg);
			
			for (i=0; i<siz/3; i++)
      {
				idx = 3*i;
				
				for (j=0; j<2; j++) {
	       	PyObject *oo = PyList_GetItem(pyEdg,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyEdg,idx+2);
				ref = (int) PyInt_AS_LONG(oo);
				
				Msh->NbrEfr++;
				AddEdge(Msh,Msh->NbrEfr,is,ref);
      }
  }
	
	
	if ( PyList_Check(pyHex) )
  {
			siz = PyList_Size(pyHex);
			
			for (i=0; i<siz/9; i++)
      {
				idx = 9*i;
				
				for (j=0; j<8; j++) {
	       	PyObject *oo = PyList_GetItem(pyHex,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyHex,idx+8);
				ref = (int) PyInt_AS_LONG(oo);
				
				Msh->NbrHex++;
				AddHexahedron(Msh,Msh->NbrHex,is,ref);
      }
  }
	
	if ( PyList_Check(pyQua) )
  {
			siz = PyList_Size(pyQua);
			
			for (i=0; i<siz/5; i++)
      {
				idx = 5*i;
				
				for (j=0; j<4; j++) {
	       	PyObject *oo = PyList_GetItem(pyQua,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyQua,idx+4);
				ref = (int) PyInt_AS_LONG(oo);
								
				Msh->NbrQua++;
				AddQuadrilateral(Msh,Msh->NbrQua,is,ref);
      }
  }
		
		
	if ( PyList_Check(pyPyr) )
  {
			siz = PyList_Size(pyPyr);
			
			for (i=0; i<siz/6; i++)
      {
				idx = 6*i;
				
				for (j=0; j<5; j++) {
	       	PyObject *oo = PyList_GetItem(pyPyr,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyPyr,idx+5);
				ref = (int) PyInt_AS_LONG(oo);
				
				Msh->NbrPyr++;
				AddPyramid(Msh,Msh->NbrPyr,is,ref);
      }
  }
	
	if ( PyList_Check(pyPri) )
  {
			siz = PyList_Size(pyPri);
			
			for (i=0; i<siz/7; i++)
      {
				idx = 7*i;
				
				for (j=0; j<6; j++) {
	       	PyObject *oo = PyList_GetItem(pyPri,idx+j);
	       	if ( PyInt_Check(oo) )
	       	{
						is[j] = (int) PyInt_AS_LONG(oo);
	       	}
				}
				
				PyObject *oo = PyList_GetItem(pyPri,idx+6);
				ref = (int) PyInt_AS_LONG(oo);
				
				Msh->NbrPri++;
				AddPrism(Msh,Msh->NbrPri,is,ref);
      }
  }	
	
	if ( PyList_Check(pyVer) )
  {
			siz = PyList_Size(pyVer);
			
			for (i=0; i<siz/3; i++)
      {
				idx = 3*i;
				
				for (j=0; j<3; j++) {
	       	PyObject *oo = PyList_GetItem(pyVer,idx+j);
	       	if ( PyFloat_Check(oo) )
	       	{
						crd[j] = (double) PyFloat_AS_DOUBLE(oo);
	       	}
				}
				Msh->NbrVer++;
				AddVertex(Msh,Msh->NbrVer,crd);
				
				//printf("ADD VERTEX %d : %lf %lf %lf\n", Msh->NbrVer, crd[0], crd[1], crd[2]);
				//exit(1);
      }
  }
	
	if ( PyList_Check(pyMarkers) )
  {
		for (i=0; i<NbrMarkers; i++){
			PyObject *oo = PyList_GetItem(pyMarkers,i);
			strcpy(Msh->Markers[i], (char*) PyString_AS_STRING(oo));	
		}
		Msh->NbrMarkers = NbrMarkers;
	}
	
	
	//--- Get Solution size and check it matches the number of vertices
	
	if ( PyList_Check(pySol) )
		siz = PyList_Size(pySol);
	
	if ( siz > 0 ) {
			
		if ( siz%Msh->NbrVer == 0 ) {
			
			Msh->SolSiz = siz/Msh->NbrVer;
			Msh->NbrFld = Msh->SolSiz;
			Msh->FldTab = (int*) malloc(sizeof(int)*Msh->SolSiz);
			for (j=0; j<Msh->NbrFld; j++){
				Msh->FldTab[j] = GmfSca;
				sprintf(Msh->SolTag[j], "scalar_%d", j);
			}
			Msh->Sol = (double*) malloc(sizeof(double)*(Msh->NbrVer+1)*Msh->SolSiz);
			memset(Msh->Sol, 0, sizeof(double)*(Msh->NbrVer+1)*Msh->SolSiz);
			
			
			Msh->Sol[0] = 0.0;
			for (i=0; i<siz; i++)
      {
       	PyObject *oo = PyList_GetItem(pySol,i);
       	if ( PyFloat_Check(oo) )
       	{
					Msh->Sol[i+Msh->SolSiz] = (double) PyFloat_AS_DOUBLE(oo);
       	}
			}
		}
		else {
			printf("  ## ERROR py_WriteMesh: Inconsistent solution provided. Skip.\n");
		}
		
	}
	
	//--- Write Mesh
	
	int FilTyp = GetInputFileType(MshNam);
  char *ptr = NULL;
	char BasNam[1024], BasNamSol[1024], OutSol[1024];
	
	// --- Get BasNam
	
  strcpy(BasNam,MshNam);
	
  ptr = strstr(BasNam,".su2");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
  ptr = strstr(BasNam,".meshb");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
	
	strcpy(BasNamSol,SolNam);
	
  ptr = strstr(BasNamSol,".dat");	
  if ( ptr != NULL )
    BasNamSol[ptr-BasNamSol]='\0';
  ptr = strstr(BasNamSol,".solb");	
  if ( ptr != NULL )
    BasNamSol[ptr-BasNamSol]='\0';
	
	if ( FilTyp != FILE_SU2 ) {
		WriteGMFMesh(BasNam, Msh, 1);
		if ( Msh->Sol ) {
			sprintf(OutSol, "%s.solb", BasNamSol);
			if ( ! WriteGMFSolutionItf(OutSol, Msh) ) {
				printf("  ## ERROR : Output solution FAILED.\n");
			}
		}
	}
	else {
		WriteSU2Mesh(BasNam, Msh);
		if ( Msh->Sol ) {		
			sprintf(OutSol, "%s.dat", BasNamSol);
			WriteSU2Solution (OutSol, Msh, Msh->Sol, Msh->NbrVer,  Msh->SolSiz, Msh->SolTag);
		}
	}		
}

void py_WriteSolution(char *SolNam, PyObject *pyVer, PyObject *pySol, PyObject *pySolHeader, int NbrVer, int Dim)
{
	
	int siz, i, j, idx;
	int SolSiz=0, NbrFld=0, NbrTag=0;
	int *FldTab = NULL;
	
	double  *Sol = NULL;
	double3 *Ver = NULL;
	
	char SolTag[100][256];
	
	if ( PyList_Check(pySol) )
		siz = PyList_Size(pySol);
	
	if ( PyList_Check(pySolHeader) )
		NbrTag = PyList_Size(pySolHeader);
	
	int FilTyp = GetInputFileType(SolNam);
	
	if ( siz > 0 ) {
			
		if ( siz%NbrVer == 0 ) {
			
			SolSiz = siz/NbrVer;
			NbrFld = SolSiz;
			FldTab = (int*) malloc(sizeof(int)*SolSiz);
			for (j=0; j<NbrFld; j++){
				FldTab[j] = GmfSca;
				
				if ( NbrTag == NbrFld ) {
	       	PyObject *oo = PyList_GetItem(pySolHeader,j);
	       	if ( PyFloat_Check(oo) )
	       	{
						sprintf(SolTag[j], "%s", (char*) PyString_AS_STRING(oo));
	       	}
				}
				else 
					sprintf(SolTag[j], "scalar_%d", j);
			}
			
			Sol = (double*) malloc(sizeof(double)*(NbrVer+1)*SolSiz);
			memset(Sol, 0, sizeof(double)*(NbrVer+1)*SolSiz);
			
			Sol[0] = 0.0;
			for (i=0; i<siz; i++)
      {
       	PyObject *oo = PyList_GetItem(pySol,i);
       	if ( PyFloat_Check(oo) )
       	{
					Sol[i+SolSiz] = (double) PyFloat_AS_DOUBLE(oo);
       	}
			}
		}
		else {
			printf("  ## ERROR py_WriteSolution: Inconsistent solution provided. Skip.\n");
			printf("siz %d NbrVer %d -> %d\n", siz, NbrVer, siz%NbrVer);
			exit(1);
		}
		
		
		if ( PyList_Check(pyVer) )
	  {
				siz = PyList_Size(pyVer);
				
				if ( NbrVer != siz/3 ) {
					printf("  ## ERROR py_WriteSolution: Inconsistent number of vertices. Skip.\n");
					exit(1);
				}
				
				Ver = (double3*) malloc(sizeof(double3)*(NbrVer+1));
			
				for (i=0; i<siz/3; i++)
	      {
					idx = 3*i;
				
					for (j=0; j<3; j++) {
		       	PyObject *oo = PyList_GetItem(pyVer,idx+j);
		       	if ( PyFloat_Check(oo) )
		       	{
							Ver[i+1][j] = (double) PyFloat_AS_DOUBLE(oo);
		       	}
						
					}
					
					//printf("ADD VERTEX %d : %lf %lf %lf\n", Msh->NbrVer, crd[0], crd[1], crd[2]);
					//exit(1);
	      }
	  }
		
		if ( FilTyp == FILE_GMFSOL ) {
			WriteGMFSolution(SolNam, Sol, SolSiz, NbrVer, Dim, NbrFld, FldTab);
		}
		else {
			WriteSU2Solution_2 (SolNam, Dim, NbrVer, Ver, Sol, SolSiz, SolTag);
		}
		
	}	
	
	if ( Sol )
		free(Sol);
	
	if ( Ver )
		free(Ver);
	
}

