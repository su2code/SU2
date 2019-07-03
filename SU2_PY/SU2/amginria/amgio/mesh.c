#include "amgio.h"

/*
Victorien Menier Feb 2016
*/

void AddEdge(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrEfr ) {
    printf("  ## ERROR : Max number of boundary edges reached.\n");
		printf("MAX = %d, idx = %d \n", Msh->MaxNbrEfr, idx);
    exit(1);
  }
  Msh->Efr[idx][0] = is[0];
  Msh->Efr[idx][1] = is[1];
  Msh->Efr[idx][2] = ref;
}


void AddHexahedron(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrHex ) {
    printf("  ## ERROR : Max number of hexas reached.\n");
    exit(1);
  }
  Msh->Hex[idx][0] = is[0];
  Msh->Hex[idx][1] = is[1];
	Msh->Hex[idx][2] = is[2];
	Msh->Hex[idx][3] = is[3];
	Msh->Hex[idx][4] = is[4];
  Msh->Hex[idx][5] = is[5];
	Msh->Hex[idx][6] = is[6];
	Msh->Hex[idx][7] = is[7];
  Msh->Hex[idx][8] = ref;
}


void AddQuadrilateral(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrQua ) {
    printf("  ## ERROR : Max number of quads reached.\n");
    exit(1);
  }
  Msh->Qua[idx][0] = is[0];
  Msh->Qua[idx][1] = is[1];
	Msh->Qua[idx][2] = is[2];
	Msh->Qua[idx][3] = is[3];
  Msh->Qua[idx][4] = ref;
}

void AddTetrahedron(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrTet ) {
    printf("  ## ERROR : Max number of tetra reached.\n");
    exit(1);
  }
  Msh->Tet[idx][0] = is[0];
  Msh->Tet[idx][1] = is[1];
	Msh->Tet[idx][2] = is[2];
	Msh->Tet[idx][3] = is[3];
  Msh->Tet[idx][4] = ref;
}

void AddPyramid(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrPyr ) {
    printf("  ## ERROR : Max number of pyramids reached.\n");
    exit(1);
  }
  Msh->Pyr[idx][0] = is[0];
  Msh->Pyr[idx][1] = is[1];
	Msh->Pyr[idx][2] = is[2];
	Msh->Pyr[idx][3] = is[3];
	Msh->Pyr[idx][4] = is[4];
  Msh->Pyr[idx][5] = ref;
}


void AddPrism(Mesh *Msh, int idx, int *is, int ref)
{
  if ( idx > Msh->MaxNbrPri ) {
    printf("  ## ERROR : Max number of prisms reached.\n");
    exit(1);
  }
  Msh->Pri[idx][0] = is[0];
  Msh->Pri[idx][1] = is[1];
	Msh->Pri[idx][2] = is[2];
	Msh->Pri[idx][3] = is[3];
	Msh->Pri[idx][4] = is[4];
	Msh->Pri[idx][5] = is[5];
  Msh->Pri[idx][6] = ref;
}

void AddTriangle(Mesh *Msh, int idxTri, int *is, int ref)
{
  if ( idxTri > Msh->MaxNbrTri ) {
    printf("  ## ERROR : Max number of triangles reached (%d, max %d).\n", idxTri, Msh->MaxNbrTri);
    exit(1);
  }
  Msh->Tri[idxTri][0] = is[0];
  Msh->Tri[idxTri][1] = is[1];
  Msh->Tri[idxTri][2] = is[2];
  Msh->Tri[idxTri][3] = ref;
}


void AddVertex(Mesh *Msh, int idxVer, double *Crd)
{	
  Msh->Ver[idxVer][0] = Crd[0];
  Msh->Ver[idxVer][1] = Crd[1];
		
	if ( Msh->Dim == 3 )
		Msh->Ver[idxVer][2] = Crd[2];
}

Mesh* AllocMesh (int * SizMsh)
{
	Mesh *Msh=NULL;
	Msh = (Mesh*)malloc(sizeof(struct S_Mesh));
	
	Msh->MaxNbrVer = SizMsh[GmfVertices];
	Msh->MaxNbrEfr = SizMsh[GmfEdges];
	Msh->MaxNbrTri = SizMsh[GmfTriangles];
	Msh->MaxNbrTet = SizMsh[GmfTetrahedra];
	Msh->MaxNbrQua = SizMsh[GmfQuadrilaterals];
	Msh->MaxNbrHex = SizMsh[GmfHexahedra];
	Msh->MaxNbrPyr = SizMsh[GmfPyramids];
	Msh->MaxNbrPri = SizMsh[GmfPrisms];
	Msh->Dim       = SizMsh[GmfDimension];
	
	Msh->NbrVer = 0;
	Msh->NbrEfr = 0;
	Msh->NbrTri = 0;
	Msh->NbrTet = 0;
	Msh->NbrQua = 0;
	Msh->NbrHex = 0;
	Msh->NbrPri = 0;
	Msh->NbrPyr = 0;
	
	Msh->Ver = NULL;
	Msh->Efr = NULL;
	Msh->Tri = NULL;
	Msh->Tet = NULL;
	Msh->Qua = NULL;
	Msh->Hex = NULL;
	Msh->Pri = NULL;
	Msh->Pyr = NULL;
	
	Msh->FldTab = NULL;
	
	Msh->FilTyp = -1;
	
	Msh->Sol = NULL;
	
	Msh->NbrMarkers = 0;
	
	if ( Msh->MaxNbrVer > 0 ) {
		Msh->Ver = (double3*)malloc(sizeof(double3)*(Msh->MaxNbrVer+1));
	}
	
	if ( Msh->MaxNbrEfr > 0 ) {
		Msh->Efr = (int3*)malloc(sizeof(int3)*(Msh->MaxNbrEfr+1));
	}
	
	if ( Msh->MaxNbrTri > 0 ) {
		Msh->Tri = (int4*)malloc(sizeof(int4)*(Msh->MaxNbrTri+1));
	}
	
	if ( Msh->MaxNbrTet > 0 ) {
		Msh->Tet = (int5*)malloc(sizeof(int5)*(Msh->MaxNbrTet+1));
	}
	
	if ( Msh->MaxNbrHex > 0 ) {
		Msh->Hex = (int9*)malloc(sizeof(int9)*(Msh->MaxNbrHex+1));
	}
	
	if ( Msh->MaxNbrQua > 0 ) {
		Msh->Qua = (int5*)malloc(sizeof(int5)*(Msh->MaxNbrQua+1));
	}
	
	if ( Msh->MaxNbrPri > 0 ) {
		Msh->Pri = (int7*)malloc(sizeof(int7)*(Msh->MaxNbrPri+1));
	}
	
	if ( Msh->MaxNbrPyr > 0 ) {
		Msh->Pyr = (int6*)malloc(sizeof(int6)*(Msh->MaxNbrPyr+1));
	}
	
	
	
	return Msh;		
}


int cmp_int2(const void *a, const void *b)
{
  int *ia = (int*) a;
  int *ib = (int*) b;
  
  if ( ia[0] < ib[0] ) {
    return -1;
  }
  else if ( ia[0] > ib[0] ) {
    return 1;
  }
  else {
    
    if ( ia[1] < ib[1] ) {
      return -1;
    }
    else if ( ia[1] > ib[1 ] ) {
      return 1;
    }
    else {
      return 0;
    }
  }
  
}


double Dist (double *crd0, double *crd1)
{
  double len;
  len = (crd1[0]-crd0[0])*(crd1[0]-crd0[0]) + (crd1[1]-crd0[1])*(crd1[1]-crd0[1]);
  len = sqrt(len);
  return len;
}


int FreeMesh(Mesh *Msh)
{
	
	
	if ( !Msh ) return 0;
	
	if ( Msh->Ver ) {
		free(Msh->Ver);
	}
	
	if ( Msh->Efr ) {
		free(Msh->Efr);
	}
	
	if ( Msh->Pyr ) {
		free(Msh->Pyr);
	}
	
	if ( Msh->Pri ) {
		free(Msh->Pri);
	}
	
	if ( Msh->Tri ) {
		free(Msh->Tri);
	}
	
	if ( Msh->Tet ) {
		free(Msh->Tet);
	}
	
	if ( Msh->Sol ) {
		free(Msh->Sol);
	}
	
	if ( Msh->FldTab ) {
		free(Msh->FldTab);
	}
	
	free(Msh);
	
	return 1;	
}


int imin(int n, int *idx)
{
  int i;
  int im   = 0;
  for(i=1; i<n; i++) {
    if ( idx[i] < idx[im] ) {
      im = i;
    }
  }
  return im;  
}

void outwardNormalCrd(double2 crd0, double2 crd1, double *vno)
{
  vno[0]  = crd1[1]-crd0[1];  // define outward normal, how to verify this ? 
  vno[1]  = crd0[0]-crd1[0];
  vno[2]  = sqrt(vno[0]*vno[0]+vno[1]*vno[1]);
}


void PrintMeshInfo (Mesh *Msh)
{
	int i;
	printf("\n-------------------------\n");
	printf("Mesh: %s\n", Msh->MshNam);
	printf("  Vertices   : %d\n", Msh->NbrVer);
	if ( Msh->NbrTri > 0 ) printf("  Triangles  : %d\n", Msh->NbrTri);
	if ( Msh->NbrQua > 0 ) printf("  Quads      : %d\n", Msh->NbrQua);
	if ( Msh->NbrEfr > 0 ) printf("  Edges      : %d\n", Msh->NbrEfr);
	if ( Msh->NbrTet > 0 ) printf("  Tetrahedra : %d\n", Msh->NbrTet);
	if ( Msh->NbrPri > 0 ) printf("  Prisms     : %d\n", Msh->NbrPri);
	if ( Msh->NbrHex > 0 ) printf("  Hexahedra  : %d\n", Msh->NbrHex);
	if ( Msh->NbrPyr > 0 ) printf("  Pyramids   : %d\n", Msh->NbrPyr);
	
	if ( Msh->Sol ){
		printf("Solution: %s\n", Msh->SolNam);
		printf("  Solution fields: %d\n", Msh->NbrFld);
		for (i=0; i<Msh->NbrFld; i++)
			printf("    %s\n", Msh->SolTag[i]);
	}
	
	printf("-------------------------\n\n");
		
}


int RemoveUnconnectedVertices(Mesh *Msh)
{

	int iVer, iTri, NbrVer, iEfr, j;
	
	int *Tag = (int*)malloc(sizeof(int)*(Msh->NbrVer+1));
	memset(Tag,0,sizeof(int)*(Msh->NbrVer+1));
	
	//if ( Msh->Dim == 2 ) {
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
			for (j=0; j<3; j++) 
				Tag[Msh->Tri[iTri][j]] = 1;
		}
	//}
	//else {
	//	printf("  ## ERROR : RemoveUnconnectedVertices : Dim 3 not implemented\n");
	//	exit(1);
	//}
	
	NbrVer=0;
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {	
		if ( Tag[iVer] != 1 ) continue;
		
		NbrVer++;
		Tag[iVer] = NbrVer;
		
		Msh->Ver[NbrVer][0] = Msh->Ver[iVer][0];
		Msh->Ver[NbrVer][1] = Msh->Ver[iVer][1];
		Msh->Ver[NbrVer][2] = Msh->Ver[iVer][2];
		
	}

	for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
		for (j=0; j<3; j++) 
			Msh->Tri[iTri][j] = Tag[Msh->Tri[iTri][j]];
	}
	
	for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++) {
		for (j=0; j<2; j++) 
			Msh->Efr[iEfr][j] = Tag[Msh->Efr[iEfr][j]];
	}	
	
	printf(" -- Info : %d unconnected vertices removed.\n", Msh->NbrVer-NbrVer);
	
	Msh->NbrVer = NbrVer;
	
	if ( Tag )
		free(Tag);
	
	return 1;
}

Mesh *SetupMeshAndSolution (char *MshNam, char *SolNam)
{
	Mesh *Msh = NULL;
	int SizMsh[GmfMaxSizMsh+1];
	int FilTyp = GetInputFileType(MshNam);
	int val;
	
	if ( !FilTyp ) {
		printf("  ## ERROR SetupMeshAndSolution : Unknown mesh format.\n");
		return NULL;
	}
	
	if ( FilTyp == FILE_SU2 )
		AddSU2MeshSize(MshNam, SizMsh);
	else if ( FilTyp == FILE_GMF )
		AddGMFMeshSize(MshNam, SizMsh);
	
	Msh = AllocMesh(SizMsh);
	
	if ( FilTyp == FILE_SU2 ) {
		LoadSU2Mesh(MshNam, Msh);
		if ( strcmp(SolNam,"") )
			val = LoadSU2Solution(SolNam, Msh);
		
		if ( val != 1 ) {
			Msh->Sol = NULL;
		}
		
	}	
	else if ( FilTyp == FILE_GMF ){
		LoadGMFMesh(MshNam, Msh);
		if ( strcmp(SolNam,"") )
			LoadGMFSolution(SolNam, Msh);
	}
	return Msh;
}


void switchHexIdx(int *idx, int *swi)
{
  int im;
  im = imin(8,idx);  
    
  switch( im ) { 
  
  case 0:
    swi[0] = idx[0];
    swi[1] = idx[1];
    swi[2] = idx[2];
    swi[3] = idx[3];
    swi[4] = idx[4];
    swi[5] = idx[5];
    swi[6] = idx[6];
    swi[7] = idx[7];
    break;
    
  case 1:
    swi[0] = idx[1];
    swi[1] = idx[2];
    swi[2] = idx[3];
    swi[3] = idx[0];
    swi[4] = idx[5];
    swi[5] = idx[6];
    swi[6] = idx[7];
    swi[7] = idx[4];
    break;
  
  case 2:
    swi[0] = idx[2];
    swi[1] = idx[3];
    swi[2] = idx[0];
    swi[3] = idx[1];
    swi[4] = idx[6];
    swi[5] = idx[7];
    swi[6] = idx[4];
    swi[7] = idx[5];
    break;
  
  case 3:
    swi[0] = idx[3];
    swi[1] = idx[0];
    swi[2] = idx[1];
    swi[3] = idx[2];
    swi[4] = idx[7];
    swi[5] = idx[4];
    swi[6] = idx[5];
    swi[7] = idx[6];
    break; 
    
  case 4:
    swi[0] = idx[4];
    swi[1] = idx[7];
    swi[2] = idx[6];
    swi[3] = idx[5];
    swi[4] = idx[0];
    swi[5] = idx[3];
    swi[6] = idx[2];
    swi[7] = idx[1];
    break; 
    
  case 5:
    swi[0] = idx[5];
    swi[1] = idx[4];
    swi[2] = idx[7];
    swi[3] = idx[6];
    swi[4] = idx[1];
    swi[5] = idx[0];
    swi[6] = idx[3];
    swi[7] = idx[2];
    break;      
   
  case 6:
    swi[0] = idx[6];
    swi[1] = idx[5];
    swi[2] = idx[4];
    swi[3] = idx[7];
    swi[4] = idx[2];
    swi[5] = idx[1];
    swi[6] = idx[0];
    swi[7] = idx[3];
    break;
    
  case 7:
    swi[0] = idx[7];
    swi[1] = idx[6];
    swi[2] = idx[5];
    swi[3] = idx[4];
    swi[4] = idx[3];
    swi[5] = idx[2];
    swi[6] = idx[1];
    swi[7] = idx[0];
    break;         
  }   
}


void switchQuaIdx(int *idx, int *swi)
{
  int im = 0;
  
  im = imin(4,idx);  
  
  switch( im ) { 
  
  case 0:
    swi[0] = idx[0];
    swi[1] = idx[1];
    swi[2] = idx[2];
    swi[3] = idx[3];
    break;
    
  case 1:
    swi[0] = idx[1];
    swi[1] = idx[2];
    swi[2] = idx[3];
    swi[3] = idx[0];
    break;
  
  case 2:
    swi[0] = idx[2];
    swi[1] = idx[3];
    swi[2] = idx[0];
    swi[3] = idx[1];
    break;
  
  case 3:
    swi[0] = idx[3];
    swi[1] = idx[0];
    swi[2] = idx[1];
    swi[3] = idx[2];
    break;  
   
  }
  
}

void switchTetIdx(int *idx, int *swi)
{
  
  if ( idx[1] < idx[0] ) {  
    
    if ( idx[2] < idx[1] ) {   
      if ( idx[3] < idx[2] ) {  
        swi[0] = idx[3];   swi[1] = idx[2];
        swi[2] = idx[1];   swi[3] = idx[0];
      }
      else if ( idx[0] < idx[3] ){  
        swi[0] = idx[2];   swi[1] = idx[0];
        swi[2] = idx[1];   swi[3] = idx[3];
      }
      else { 
        swi[0] = idx[2];   swi[1] = idx[1];
        swi[2] = idx[3];   swi[3] = idx[0];
      }
    }   
    else if ( idx[0] < idx[2] ) {   
      if ( idx[3] < idx[1] ) {  
        swi[0] = idx[3];   swi[1] = idx[1];
        swi[2] = idx[0];   swi[3] = idx[2];
      }
      else if ( idx[2] < idx[3] ){ 
        swi[0] = idx[1];   swi[1] = idx[2];
        swi[2] = idx[0];   swi[3] = idx[3];
      }
      else {  
        swi[0] = idx[1];   swi[1] = idx[0];
        swi[2] = idx[3];   swi[3] = idx[2];
      }
    }
    
    else { 
      if ( idx[3] < idx[1] ) {  
        swi[0] = idx[3];   swi[1] = idx[2];
        swi[2] = idx[1];   swi[3] = idx[0];
      }
      else if ( idx[0] < idx[3] ){  
        swi[0] = idx[1];   swi[1] = idx[2];
        swi[2] = idx[0];   swi[3] = idx[3];
      }
      else {  
        swi[0] = idx[1];   swi[1] = idx[3];
        swi[2] = idx[2];   swi[3] = idx[0];
      }
    }   

  }
  
  else {   
    
    if ( idx[2] < idx[0] ) {  
      if ( idx[3] < idx[2] ) { 
        swi[0] = idx[3];   swi[1] = idx[0];
        swi[2] = idx[2];   swi[3] = idx[1];
      }
      else if ( idx[1] < idx[3] ){  
        swi[0] = idx[2];   swi[1] = idx[0];
        swi[2] = idx[1];   swi[3] = idx[3];
      }
      else { 
        swi[0] = idx[2];   swi[1] = idx[3];
        swi[2] = idx[0];   swi[3] = idx[1];
      }
    }   
    else if ( idx[1] < idx[2] ) {   
      if ( idx[3] < idx[0] ) { 
        swi[0] = idx[3];   swi[1] = idx[1];
        swi[2] = idx[0];   swi[3] = idx[2];
      }
      else if ( idx[2] < idx[3] ){  
        swi[0] = idx[0];   swi[1] = idx[1];
        swi[2] = idx[2];   swi[3] = idx[3];
      }
      else {  
        swi[0] = idx[0];   swi[1] = idx[3];
        swi[2] = idx[1];   swi[3] = idx[2];
      }
    }
    
    else {  
      if ( idx[3] < idx[0] ) { 
        swi[0] = idx[3];   swi[1] = idx[0];
        swi[2] = idx[2];   swi[3] = idx[1];
      }
      else if ( idx[1] < idx[3] ){ 
        swi[0] = idx[0];   swi[1] = idx[1];
        swi[2] = idx[2];   swi[3] = idx[3];
      }
      else { 
        swi[0] = idx[0];   swi[1] = idx[2];
        swi[2] = idx[3];   swi[3] = idx[1];
      }
    }   
    
  }

  //printf(" in  : %d %d %d %d \n",idx[0],idx[1],idx[2],idx[3]);
  //printf(" out : %d %d %d %d \n",swi[0],swi[1],swi[2],swi[3]);
  
}


void switchTriIdx(int *idx, int *swi)
{
  if ( idx[1] < idx[0] ) { 
    if ( idx[2] < idx[1] ) {  
      swi[0] = idx[2];
      swi[1] = idx[0];
      swi[2] = idx[1];
    }   
    else { 
      swi[0] = idx[1];
      swi[1] = idx[2];
      swi[2] = idx[0];
    }   
  }
  
  else {
    if ( idx[2] < idx[0] ) {  
        swi[0] = idx[2];
        swi[1] = idx[0];
        swi[2] = idx[1];
    }
    else {  
      swi[0] = idx[0];
      swi[1] = idx[1];
      swi[2] = idx[2];
    }
  }

}


int GetInputFileType (char *FilNam) 
{
  char *ext=NULL;  
	
  if ( !FilNam ) {
    printf("\n ## ERROR GetInputFileType : No file name provided.\n\n");
   	exit(1);
  }
  
  ext = strrchr(FilNam, '.');
  
  if (!ext) {
    return 0;    
  } else {      
      ext = ext+1;
      if ( strcmp(ext,"su2") == 0  ) {
        return FILE_SU2;
      }
      else if ( strcmp(ext,"dat") == 0  ) {
        return FILE_DAT;
      }
      else if ( strcmp(ext,"mesh") == 0 || strcmp(ext,"meshb") == 0 ) {
        return FILE_GMF;
      }
      else if ( strcmp(ext,"sol") == 0 || strcmp(ext,"solb") == 0 ) {
        return FILE_GMFSOL;
      }
      else if ( strcmp(ext,"geo") == 0  ) {
        return FILE_GEO;
      }
      else {
				return 0;
      }
  }
}


//--- Transforms all letters to lower case
int Str2Lower(char *buff)
{
  int iChr;
  
  for (iChr=0; iChr<strlen(buff); iChr++) 
    buff[iChr] = tolower( buff[iChr] );

  return 1;
}

//--- Removes all occurences of char c from str
void StrRemoveChars (char* s, char ch) {
	char *p = s;
	while (*s) {
	    if (*s != ch)
	        *p++ = *s;
	    s++;
	}
	*p = 0;
}
