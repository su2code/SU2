#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <setjmp.h>

#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)

#include "libmesh6.h"
#include "mesh.h"
#include "option.h"

#define PI_NUMBER 3.14159265359

// Prototypes

//--- SU2io.c
int  AddSU2MeshSize (char *FilNam, int *SizMsh) ;
int  GetSU2KeywordValue (FILE *FilHdl, char *Kwd);
int  GetSU2KeywordValueStr (FILE *FilHdl, char *Kwd, char *StrVal);
int  LoadSU2Elements(FILE *FilHdl, Mesh *Msh);
int  LoadSU2Mesh (char *FilNam, Mesh *Msh);
int  LoadSU2Solution(char *SolNam, Mesh *Msh);
int  LoadSU2Vertices (FILE *FilHdl, Mesh *Msh);
void WriteSU2Mesh(char *nam, Mesh *Msh);
int  WriteSU2Solution (char *SolNam, Mesh *Msh, double *Sol, int NbrVer, int SolSiz, char SolTag[100][256]);
int  GetSU2SolSize(char *SolNam);
int  WriteSU2Solution_2 (char *SolNam, int Dim, int NbrVer, double3 *Ver,  double *Sol, int SolSiz, char SolTag[100][256]);

//--- GMFio.c
int AddGMFMeshSize (char *MshNam, int *SizMsh);
int LoadGMFMesh (char *MshNam, Mesh *Msh);
int LoadGMFSolution(char *SolNam, Mesh *Msh);
int WriteGMFMesh(char *nam, Mesh *Msh, int OptBin);
int WriteGMFSolution(char *SolNam, double *Sol, int SolSiz, int NbrVer, int Dim, int NbrFld, int* FldTab);
int WriteGMFSolutionItf(char *SolNam, Mesh *Msh);

//--- option.c
Options* AllocOptions(void);
int      CheckOptions (Options *mshopt);
int      GetBasNam (char *InpNam, char *BasNam);
void     PrintOptions (Options *mshopt);

//---- mesh.c
Mesh* AllocMesh (int * SizMsh);
int   cmp_int2(const void *a, const void *b);
int   FreeMesh (Mesh *Msh);
int   GetMeshSize (char *MshNam, int *SizMsh);
Mesh *SetupMeshAndSolution (char *MshNam, char *SolNam);
int   RemoveUnconnectedVertices (Mesh *Msh);
void  AddEdge(Mesh *Msh, int idx, int *is, int ref);
void  AddHexahedron(Mesh *Msh, int idx, int *is, int ref);
void  AddQuadrilateral(Mesh *Msh, int idx, int *is, int ref);
void  AddTetrahedron(Mesh *Msh, int idx, int *is, int ref);
void  AddPyramid(Mesh *Msh, int idx, int *is, int ref);
void  AddPrism(Mesh *Msh, int idx, int *is, int ref);
void  AddTriangle(Mesh *Msh, int idxTri, int *is, int ref);
void  AddVertex(Mesh *Msh, int idxVer, double *Crd);
int   imin(int n, int *idx);
void  PrintMeshInfo (Mesh *Msh);
void  switchHexIdx(int *idx, int *swi);
void  switchQuaIdx(int *idx, int *swi);
void  switchTetIdx(int *idx, int *swi);
void  switchTriIdx(int *idx, int *swi);
int  GetInputFileType (char *FilNam);
int  Str2Lower(char *buff);
void StrRemoveChars (char* str, char c);

//--- convert.c
int  SplitSolution (Mesh *Msh, char *prefix, char *adap_sensor);
