#include "amgio.h"

/*
Victorien Menier Feb 2016
*/

int AddSU2MeshSize(char *FilNam, int *SizMsh) 
{
	int i, NbrElt, iElt, typ, CptElt;
  int NbrTri, NbrTet, NbrHex, NbrPyr, NbrRec, NbrLin, NbrWed, NbrP2Tri, NbrP2Lin;
  int NbrMark, iMark;
  FILE *FilHdl = NULL;
  char    str[1024];
  

  for (i=0; i<GmfMaxSizMsh; i++)
    SizMsh[i] = 0;
	
	
	if (  GetInputFileType (FilNam) != FILE_SU2 ) {
		printf("NOT SU2, return\n");
		return 0;
	}
		
	
	FilHdl = fopen (FilNam, "r");
	
	if ( !FilHdl ) {
		fprintf(stderr, "  -- Info : Tried to open %s. Failed.\n", FilNam );
		return 0;
	}
	
	//--- Dim ?

	rewind(FilHdl);
	SizMsh[GmfDimension] = GetSU2KeywordValue (FilHdl, "NDIME=");	
		
  //--- Elements?
  
  NbrTri = NbrTet = NbrHex = NbrPyr = NbrLin = NbrRec = NbrWed = NbrP2Tri = NbrP2Lin =  0;
	
	NbrElt = GetSU2KeywordValue (FilHdl, "NELEM=");
	
  for (iElt=0; iElt<NbrElt; iElt++) {
    fscanf(FilHdl, "%d", &typ);
    
    if ( typ == SU2_TRIANGLE ) {   
      NbrTri++; 
    }
    else if ( typ == SU2_TETRAHEDRAL ) {
      NbrTet++; 
    }
    else if ( typ == SU2_HEXAHEDRAL ) {
      NbrHex++;
    }
    else if ( typ == SU2_PYRAMID ) {
      NbrPyr++;
    }
    else if ( typ == SU2_RECTANGLE ) {
      NbrRec++; 
    }
    else if ( typ == SU2_WEDGE ) {
      NbrWed++;
    }
    else if ( typ == SU2_LINE ) {
      NbrLin++;
    }
		else if ( typ == SU2_TRIANGLEP2 ) {
			NbrP2Tri++;
		}
		else if ( typ == SU2_LINEP2 ) {
      NbrP2Lin++;
    }
		else {
			printf("  ## ERROR : AddSU2MeshSize: Unknown element type %d\n", typ);
			return 0;
		}
		fgets (str, sizeof str, FilHdl);
  }//for iElt
	
 	rewind(FilHdl);
	SizMsh[GmfVertices] = GetSU2KeywordValue (FilHdl, "NPOIN=");
  
  //--- Boundary Elements?
	NbrMark = 0;	
	
	rewind(FilHdl);
	NbrMark = GetSU2KeywordValue (FilHdl, "NMARK=");
	
	
	
	for (iMark=1; iMark<=NbrMark; iMark++) {
		
		GetSU2KeywordValueStr (FilHdl, "MARKER_TAG=", str);
		
		if ( !strcmp(str,"SEND_RECEIVE") ) {
			printf("      Tag %s was ignored.\n", str);
			continue;
		}
		
		CptElt = GetSU2KeywordValue (FilHdl, "MARKER_ELEMS=");
		
		
	  for (iElt=0; iElt<CptElt; iElt++) {
	    fscanf(FilHdl, "%d", &typ);
	    if ( typ == SU2_TRIANGLE ) {
	      NbrTri++;
	    }
	    else if ( typ == SU2_RECTANGLE ) {
	      NbrRec++;
	    }
	    else if ( typ == SU2_LINE ) {
	      NbrLin++;
	    }
	    else if ( typ == SU2_LINEP2 ) {
	      NbrP2Lin++;
	    }
			else {
				printf("  ## ERROR : AddSU2MeshSize : Unknown boundary element type %d\n", typ);
				return 0;
			}
	  	fgets (str, sizeof str, FilHdl);
		}
		
	}
	
  SizMsh[GmfTriangles]      = NbrTri+NbrP2Tri;
  SizMsh[GmfTetrahedra]     = NbrTet;
  SizMsh[GmfEdges]          = NbrLin+NbrP2Lin;
  SizMsh[GmfPrisms]         = NbrWed;
  SizMsh[GmfQuadrilaterals] = NbrRec;
  SizMsh[GmfPyramids]       = NbrPyr;
  SizMsh[GmfHexahedra]      = NbrHex;
	SizMsh[GmfTrianglesP2]    = NbrP2Tri;
	SizMsh[GmfFileType]       = FILE_SU2;
	SizMsh[GmfEdgesP2]        = NbrP2Lin;
	
	if ( FilHdl )
  	fclose(FilHdl);
	
	return 1;
}




int GetSU2KeywordValue (FILE *FilHdl, char *Kwd)
{
	
	size_t lenKwd=0, len=0;
	int buf=0, res=0;
	char str[1024], str2[1024], kwd[1024];
	
	if ( !FilHdl || !Kwd ) return 0;
	
	sprintf(kwd,"%s",Kwd);
	
	lenKwd = strlen(kwd);
		
	//rewind(FilHdl);
  do
	{
		res = fscanf(FilHdl, "%s", str);
	}while( (res != EOF) && strncmp(str, kwd, lenKwd) );
	
  if ((res == EOF)) {
		fprintf(stderr,"  ## ERROR: INVALID SU2 FILE (CHECK KEYWORD: %s).\n", Kwd);
		return 0;
	}
	
	len = strlen(str);
	
	if ( len > lenKwd ) {
		strncpy(str2, &str[lenKwd], len-lenKwd+1);
		//sprintf(str2, "%s", str[lenKwd]);
			
		buf = atoi(str2);
	}
	else {
		fscanf(FilHdl, "%d", &buf);
	}
	
	fgets (str, sizeof str, FilHdl);
	
	
	return buf;
}


int GetSU2KeywordValueStr (FILE *FilHdl, char *Kwd, char *StrVal)
{
	
	size_t lenKwd=0, len=0;
	int  res=0;
	char str[1024], kwd[1024], buf[1024];
	
	if ( !FilHdl || !Kwd ) return 0;
	
	strcpy(kwd,Kwd);
	lenKwd = strlen(kwd);
		
  do
	{
		res = fscanf(FilHdl, "%s", str);
	}while( (res != EOF) && strncmp(str, kwd, lenKwd) );
	
  if ((res == EOF)) {
		fprintf(stderr,"  ## ERROR: INVALID SU2 FILE (CHECK KEYWORD: %s).\n", Kwd);
		return 0;
	}
	
	len = strlen(str);
	
	if ( len > lenKwd ) {
		strncpy(StrVal, &str[lenKwd], len-lenKwd);
	}
	else {
		fscanf(FilHdl, "%s", buf);
		sprintf(StrVal, "%s", buf);
	}
	
	fgets (str, sizeof str, FilHdl);
	
	return 1;
}


int LoadSU2Elements(FILE *FilHdl, Mesh *Msh)
{
	int  ref=1;
	char   str[1024];

	int iMark, NbrMark=0, CptElt;
	int iElt, NbrElt=0, typ, is[8], swi[8], buf, s, idx, res;

  Msh->NbrTri = Msh->NbrTet = Msh->NbrHex =  Msh->NbrEfr = Msh->NbrQua = Msh->NbrPyr = Msh->NbrPri = 0;
  
  rewind(FilHdl);
  do
	{
		res = fscanf(FilHdl, "%s", str);
	}while( (res != EOF) && strcmp(str, "NELEM=") );
	
	fscanf(FilHdl, "%d", &NbrElt);
  fgets (str, sizeof str, FilHdl);
	
  idx=0;
	
  for (iElt=0; iElt<NbrElt; iElt++) {
    fscanf(FilHdl, "%d", &typ);

    if ( typ == SU2_TRIANGLE ) {
      
      for (s=0; s<3; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
				if ( swi[s] > Msh->NbrVer ) {
					printf("  ## ERROR LoadSU2Elements: vertex out of bound (vid=%d)\n", swi[s]);
					return 0;
				}
      }

			fscanf(FilHdl, "%d", &buf);
      
      Msh->NbrTri++;

			if ( Msh->NbrTri > Msh->MaxNbrTri ) {
				printf("  ## ERROR LoadSU2Elements: triangle out of bound (tid=%d, max=%d)\n", Msh->NbrTri, Msh->MaxNbrTri);
				return 0;
			}
			
			switchTriIdx(swi,is);
      AddTriangle(Msh,Msh->NbrTri,is,ref);

    }
    else if ( typ == SU2_TETRAHEDRAL ) {
      for (s=0; s<4; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
      }
			
      idx++;
			Msh->NbrTet++;
			
			if ( Msh->NbrTet > Msh->MaxNbrTet ) {
				printf("  ## ERROR LoadSU2Elements: tetra out of bound (tid=%d)\n", Msh->NbrTet);
				return 0;
			}
			
			switchTetIdx(swi,is);
			AddTetrahedron(Msh,Msh->NbrTet,is,ref);
    }
    else if ( typ == SU2_HEXAHEDRAL ) {
      for (s=0; s<8; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
      }
      fscanf(FilHdl, "%d", &idx);
      Msh->NbrHex++;
			
			if ( Msh->NbrHex > Msh->MaxNbrHex ) {
				printf("  ## ERROR LoadSU2Elements: hexahedron out of bound (hid=%d)\n", Msh->NbrHex);
				return 0;
			}
			
			switchHexIdx(swi,is);
			AddHexahedron(Msh,Msh->NbrHex,is,ref);
    }
    else if ( typ == SU2_PYRAMID ) {
      for (s=0; s<5; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
				is[s] = buf+1;
      }
      fscanf(FilHdl, "%d", &idx);
			Msh->NbrPyr++;
    
			if ( Msh->NbrPyr > Msh->MaxNbrPyr ) {
				printf("  ## ERROR LoadSU2Elements: pyramid out of bound (id=%d)\n", Msh->NbrPyr);
				return 0;
			}
    
			//switchPyrIdx(swi,is);
      AddPyramid(Msh,Msh->NbrPyr,is,ref);
			
			int i=Msh->NbrPyr;
			if ( i== 1 )
				printf("PYR %d : %d %d %d %d %d\n", i, Msh->Pyr[i][0], Msh->Pyr[i][1], Msh->Pyr[i][2], Msh->Pyr[i][3], Msh->Pyr[i][4]);
    
    }
    else if ( typ == SU2_RECTANGLE ) {
      for (s=0; s<4; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
      }
			
      Msh->NbrQua++;
    
			if ( Msh->NbrQua > Msh->MaxNbrQua ) {
				printf("  ## ERROR LoadSU2Elements: quad out of bound (id=%d)\n", Msh->NbrQua);
				return 0;
			}
    
			switchQuaIdx(swi,is);
      //DefaultQuadrilateral(Msh,Msh->NbrQua,&Msh->Qua[NbrQua],is,ref);
			AddQuadrilateral(Msh,Msh->NbrQua,is,ref);
    }
    else if ( typ == SU2_WEDGE ) {
			
      for (s=0; s<6; s++) {
        fscanf(FilHdl, "%d", &buf);
        swi[s] = buf+1;
				is[s] = buf+1;
      }
      fscanf(FilHdl, "%d", &idx);
     	Msh->NbrPri++;
    
			if ( Msh->NbrPri > Msh->MaxNbrPri ) {
				printf("  ## ERROR LoadSU2Elements: prism out of bound (id=%d)\n", Msh->NbrPri);
				return 0;
			}
    
			//switchPriIdx(swi,is);
      //DefaultPrism(Msh,Msh->NbrPri,&Msh->Pri[NbrPri],is,ref);
			AddPrism(Msh,Msh->NbrPri,is,ref);
    }
    //else if ( typ == SU2_LINE ) {
    //  for (s=0; s<2; s++) {
    //    fscanf(FilHdl, "%d", &buf);
    //    swi[s] = buf+1;
    //  }
		//	NbrEfr++;
		//				
		//	if ( NbrEfr > Msh->MaxNbrEfr ) {
		//		printf("  ## ERROR LoadSU2Elements: boundary edge out of bound (id=%d, max=%d)\n", NbrEfr,Msh->MaxNbrEfr );
		//		return 0;
		//	}
    //
    //   DefaultBdyEdge(Msh,Msh->NbrEfr,&Msh->Efr[NbrEfr],swi,ref);  
    //}
    //else if ( typ == SU2_LINEP2 ) {
    //  
    //  fscanf(FilHdl, "%d %d %d", &swi[0], &swi[1], &swi[2]);
		//	for (s=0; s<3; s++) 
    //  	swi[s]++;
    //
		//	NbrP2Efr++;
		//	
		//	if ( NbrP2Efr > Msh->MaxNbrP2Efr ) {
		//		printf("  ## ERROR LoadSU2Elements: P2 boundary edge out of bound (id=%d, max=%d)\n", NbrP2Efr,Msh->MaxNbrP2Efr );
		//		return 0;
		//	}
		//	
		//	switchP2EfrIdx(swi,is);
		//	NbrEfr++;
		//	if ( NbrEfr > Msh->MaxNbrEfr ) {
		//		printf("  ## ERROR LoadSU2Elements: P1 boundary edge out of bound (id=%d, max=%d)\n", NbrEfr,Msh->MaxNbrEfr );
		//		return 0;
		//	}
    //  DefaultBdyEdge(Msh,NbrEfr,&Msh->Efr[NbrEfr],is,ref);
    //  DefaultP2BdyEdge(Msh,NbrP2Efr,&Msh->Efr[NbrP2Efr],is[2],NbrEfr,&Msh->Efr[NbrEfr]);  
    //}
		//else if ( typ == SU2_TRIANGLEP2 ) {
    //
    //  //fscanf(FilHdl, "%d %d %d %d %d %d", &swi[0],&swi[3],&swi[1],&swi[4],&swi[2],&swi[5]);
		//	//fscanf(FilHdl, "%d %d %d %d %d %d", &swi[3],&swi[0],&swi[4],&swi[1],&swi[5],&swi[2]);
		//	//fscanf(FilHdl, "%d %d %d %d %d %d", &swi[0],&swi[1],&swi[2],&swi[3],&swi[4],&swi[5]);
		//	fscanf(FilHdl, "%d %d %d %d %d %d", &swi[0],&swi[3],&swi[1],&swi[5],&swi[4],&swi[2]);
		//	for (s=0; s<6; s++)	
    //  	swi[s]++;
    //  
		//	
		//	NbrP2Tri++;
    //
		//	if ( NbrP2Tri > Msh->MaxNbrP2Tri ) {
		//		printf("  ## ERROR LoadSU2Elements: P2 triangle out of bound (id=%d, max=%d)\n", NbrP2Tri,Msh->MaxNbrP2Tri );
		//		return 0;
		//	}
		//	
		//	switchP2TriIdx(swi,is);
		//	
		//	
		//	NbrTri++;
		//	
		//	if ( NbrTri > Msh->MaxNbrTri ) {
		//		printf("  ## ERROR LoadSU2Elements: P1 triangle out of bound (id=%d, max=%d)\n", NbrTri,Msh->MaxNbrTri );
		//		return 0;
		//	}
		//	
    //  DefaultTriangle(Msh,NbrTri,&Msh->Tri[NbrTri],is,ref);
    //  DefaultP2Triangle(Msh,NbrP2Tri,&Msh->P2Tri[NbrP2Tri],&is[3],NbrTri,&Msh->Tri[NbrTri]);
    //
    //}
	  //
		else {
			printf("  ## ERROR : Unknown element type %d\n", typ);
			return 0;
		}
    
    fgets (str, sizeof str, FilHdl); 

  }//for iElt

	//--- Read boundary elements
	
	//rewind(FilHdl);
  //do
	//{
	//	res = fscanf(FilHdl, "%s", str);
	//}while( (res != EOF) && strcmp(str, "NMARK=") );
	//fscanf(FilHdl, "%d", &NbrMark);
	//fgets (str, sizeof str, FilHdl); 
	
	rewind(FilHdl);
 	NbrMark = GetSU2KeywordValue (FilHdl, "NMARK=");
		
  for (iMark=1; iMark<=NbrMark; iMark++) {
		
		GetSU2KeywordValueStr (FilHdl, "MARKER_TAG=", str);
		//printf("      Tag %s becomes reference %d\n", str, iMark);
		
		
		if ( !strcmp(str,"SEND_RECEIVE") ) {
			printf("      Tag %s was ignored.\n", str);
			continue;
		}
		
		CptElt = GetSU2KeywordValue (FilHdl, "MARKER_ELEMS=");
		
    for (iElt=0; iElt<CptElt; iElt++) {
		
      fscanf(FilHdl, "%d", &typ);
      
      if ( typ == SU2_TRIANGLE ) {
        for (s=0; s<3; s++) {
          fscanf(FilHdl, "%d", &buf);
          swi[s] = buf+1;
        }

	      Msh->NbrTri++;

				if ( Msh->NbrTri > Msh->MaxNbrTri ) {
					printf("  ## ERROR LoadSU2Elements: triangle out of bound (tid=%d, max=%d)\n", Msh->NbrTri, Msh->MaxNbrTri);
					return 0;
				}

				switchTriIdx(swi,is);
	      AddTriangle(Msh,Msh->NbrTri,is,iMark);
      }
      else if ( typ == SU2_RECTANGLE ) {
								
	      for (s=0; s<4; s++) {
	        fscanf(FilHdl, "%d", &buf);
	        swi[s] = buf+1;
	      }
				
	      Msh->NbrQua++;

				if ( Msh->NbrQua > Msh->MaxNbrQua ) {
					printf("  ## ERROR LoadSU2Elements: quad out of bound (id=%d)\n", Msh->NbrQua);
					return 0;
				}

				switchQuaIdx(swi,is);
	      //DefaultQuadrilateral(Msh,Msh->NbrQua,&Msh->Qua[NbrQua],is,iMark);
				AddQuadrilateral(Msh,Msh->NbrQua,is,iMark);
      }
      else if ( typ == SU2_LINE ) {
	      for (s=0; s<2; s++) {
	        fscanf(FilHdl, "%d", &buf);
	        swi[s] = buf+1;
	      }
				Msh->NbrEfr++;
				
				if ( Msh->NbrEfr > Msh->MaxNbrEfr ) {
					printf("  ## ERROR LoadSU2Elements: boundary edge out of bound (id=%d, max=%d)\n", Msh->NbrEfr, Msh->MaxNbrEfr);
					return 0;
				}
				
	      //DefaultBdyEdge(Msh,Msh->NbrEfr,&Msh->Efr[NbrEfr],swi,iMark);
				AddEdge(Msh,Msh->NbrEfr,swi,iMark);
      }
	    //else if ( typ == SU2_LINEP2 ) {
      //
	    //  fscanf(FilHdl, "%d %d %d", &swi[0], &swi[2], &swi[1]);
			//	for (s=0; s<3; s++) 
	    //  	swi[s]++;
      //
			//	NbrP2Efr++;
			//	if ( NbrP2Efr > Msh->MaxNbrP2Efr ) {
			//		printf("  ## ERROR LoadSU2Elements: P2 boundary edge out of bound (id=%d, max=%d)\n", NbrP2Efr,Msh->MaxNbrP2Efr );
			//		return 0;
			//	}
			//	
			//	switchP2EfrIdx(swi,is);
			//	
			//	NbrEfr++;
			//	if ( NbrEfr > Msh->MaxNbrEfr ) {
			//		printf("  ## ERROR LoadSU2Elements: P1 boundary edge out of bound (id=%d, max=%d)\n", NbrEfr,Msh->MaxNbrEfr );
			//		return 0;
			//	}
			//	
	    //  DefaultBdyEdge(Msh,NbrEfr,&Msh->Efr[NbrEfr],is,iMark);
	    //  DefaultP2BdyEdge(Msh,NbrP2Efr,&Msh->Efr[NbrP2Efr],is[2],NbrEfr,&Msh->Efr[NbrEfr]);  
	    //}
			else {
				printf("  ## ERROR : Unknown element type %d\n", typ);
				return 0;
			}
			fgets (str, sizeof str, FilHdl);
			
    }
	}
  
	
	
	return 1;
}


int LoadSU2Mesh(char *FilNam, Mesh *Msh)
{
  
	FILE *FilHdl=NULL;
  
  if ( (Msh == NULL) || (FilNam == NULL) ) {
    printf("  ## ERROR: LOADMESH: MESH/FILE NAME NOT ALLOCATED \n");
    return 0; 
  }
	
	if (  GetInputFileType (FilNam) != FILE_SU2 )
		return 0;
	
	FilHdl = fopen (FilNam, "r");
	
	if ( !FilHdl ) {
		fprintf(stderr, "  -- Info : Tried to open %s. Failed.\n", FilNam );
		return 0;
	}
	
  //printf("  %%%% %s OPENED (READ)\n",FilNam);
	strcpy(Msh->MshNam, FilNam);
	Msh->FilTyp = FILE_SU2;
	
	rewind(FilHdl);
	Msh->Dim = GetSU2KeywordValue (FilHdl, "NDIME=");
	
	rewind(FilHdl);
	
	if ( Msh->Dim != 2 && Msh->Dim != 3 ) {
		fprintf(stderr, "  ## ERROR SU2: WRONG DIMENSION NUMBER FOR MESH FILE %s (DIM=%d).\n", FilNam, Msh->Dim);
		return 0;
	}
	
	//--- Read vertices
	LoadSU2Vertices (FilHdl, Msh);
	
	//--- Read Elements
	LoadSU2Elements(FilHdl, Msh);
	
	if ( FilHdl )
		fclose(FilHdl);

  return 1;
}


int GetSU2SolSize(char *SolNam)
{ 
	int NbrLin=0;
	char *tok=NULL, *lin=NULL;
	char str[1024];
	
	int i, iVer, iVerMax;
	size_t  len = 0;
	FILE *FilHdl=NULL;
		
	FilHdl = fopen(SolNam,"r");
	
	// Skip header
	getline(&lin, &len, FilHdl);

	//--- Count
	
	iVerMax = -1;
	
	NbrLin=0;
	while ( getline(&lin, &len, FilHdl) != -1 ) {
		NbrLin++;
		tok = strtok (lin, "	,");
		iVer = atoi(tok)+1;
		iVerMax = max(iVer, iVerMax);		
	}
	return iVerMax;
}


int LoadSU2Solution(char *SolNam, Mesh *Msh)
{
  int NbrLin=0;
	char *tok=NULL, *lin=NULL;
	char str[1024];
	
	int i, iVer, skip=0, SolSiz=0, idx, idxVer;
	size_t  len = 0;
	FILE *FilHdl=NULL;
	
	double *Sol = NULL;
	
	if ( Msh->Sol )
	{
		printf("  ## ERROR LoadSU2Solution : Msh->Sol already allocated.\n");
		return 0;
	}
	
  if ( (Msh == NULL) || (SolNam == NULL) ) {
    printf("  ## ERROR: LoadSU2Solution : MESH/FILE NAME NOT ALLOCATED \n");
    return 0; 
  }
	
  //sprintf(BasNam,"%s", FilNam);
  //ptr = strstr(BasNam,".su2");	
  //if ( ptr != NULL )
  //  BasNam[ptr-BasNam]='\0';
	//sprintf(SolNam, "%s.dat", BasNam);
	
	FilHdl = fopen(SolNam,"r");
	
	if ( !FilHdl ) {
		fprintf(stderr,"  ## ERROR: LOADSU2SOL: UNABLE TO OPEN %s. \n", SolNam);
    return 0;
	}
	
	//printf("  %%%% %s OPENED\n",SolNam);
	
	strcpy(Msh->SolNam, SolNam);
	
	if ( getline(&lin, &len, FilHdl) != -1 ) {
		tok = strtok (lin, "	,");
		skip = 0;
		SolSiz = 0;
		while ( tok ) {
			if ( !strcmp(tok,"\"PointID\"") || !strcmp(tok,"\"x\"") || !strcmp(tok,"\"y\"") || !strcmp(tok,"\"z\"")   ) {
				tok = strtok (NULL, "	,");
				skip++;
				continue;
			}
			
			strcpy(Msh->SolTag[SolSiz], tok);
			//Str2Lower(Msh->SolTag[SolSiz]);
			StrRemoveChars(Msh->SolTag[SolSiz], '\"');
			StrRemoveChars(Msh->SolTag[SolSiz], '\n');
			SolSiz++;
			tok = strtok (NULL, "	,");
		}
  }
	
	//--- Allocate Msh->Sol
	
	Msh->Sol = (double*) malloc(sizeof(double)*(Msh->NbrVer+1)*SolSiz);
	memset(Msh->Sol, 0, sizeof(double)*(Msh->NbrVer+1)*SolSiz);
	
	Sol = Msh->Sol;
	Msh->SolSiz = SolSiz;
	
	//--- Set Msh->FldTab
	// Note: Set all fields as scalars (later: vectors for velocity)
	Msh->FldTab = (int*) malloc(sizeof(int)*SolSiz);
	Msh->NbrFld = SolSiz;
	for (i=0; i<SolSiz; i++)
		Msh->FldTab[i] = GmfSca;
	
	//--- Load vertex solution
	
	NbrLin=0;
	while ( getline(&lin, &len, FilHdl) != -1 ) {
		
		NbrLin++;
		tok = strtok (lin, "	,");
		
		i=0, idx=0;
		while ( tok ) {
			
			if ( i == 0 ) {
				iVer = atoi(tok)+1;
				idxVer = iVer*SolSiz;
				if ( iVer > Msh->NbrVer ) {
					fprintf(stderr,"  ## ERROR: LOADSU2SOL: VERTEX OUT OF BOUND. \n");
			    return 0;
				}
			}
			else if ( i >= skip ) {
				Sol[idxVer+idx] = atof(tok);
				if ( Sol[idxVer+idx] != Sol[idxVer+idx] ) 
					Sol[idxVer+idx] = 0;
				idx++;
			}
			
			tok = strtok (NULL, "	,");
			i++;
			
			if ( idx == SolSiz )
				break;
		}
	}
	
	
	if ( FilHdl )
		fclose(FilHdl);
		
	if ( NbrLin != Msh->NbrVer ) {
		fprintf(stderr,"  ## ERROR: LOADSU2SOL (%s): INCONSISTENT NUMBER OF VERTICES. \n", SolNam);
		return 0;
	}

	return 1;
}



int LoadSU2Vertices(FILE *FilHdl, Mesh *Msh)
{
	int iVer, d, ref;
	double crd[3], bufDbl;
	char str[1024];
	
	rewind(FilHdl);
	
	Msh->NbrVer = GetSU2KeywordValue (FilHdl, "NPOIN=");
	
	if ( Msh->NbrVer > Msh->MaxNbrVer ) {
		printf("  ## ERROR: LoadSU2Vertices: INCONSISTENT NUMBER OF VERTICES.\n");
		return 0;
	}
	
  for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
		
    crd[2] = 0;
		
    for (d=0; d<Msh->Dim; d++) {
      fscanf(FilHdl, "%lf", &bufDbl);
      crd[d] = bufDbl;
    }
    
    fscanf(FilHdl, "%d", &ref);
    fgets (str, sizeof str, FilHdl);
    
		AddVertex(Msh, iVer, crd);
  }
	
	return 1;
}




void WriteSU2Mesh(char *nam, Mesh *Msh)
{
  int       i, j, s, idx;
  int       iVer,iTri,iEfr, iTet, iHex, iPri, iPyr, iQua, NbrElt=0;
  char      OutNam[512];
  
  int Dim = Msh->Dim;
	
	int3 *BdrTag=NULL;
	int NbrBdr, NbrTag, start, iTag, cpt;
	
	FILE *OutFil=NULL;
	
	GetBasNam (nam, OutNam);
	strcat(OutNam, ".su2");
 
	OutFil = fopen(OutNam, "wb");
	
	if ( !OutFil ) {
		printf("  ## ERROR Write SU2: Can't open %s\n", OutNam);
	}
 	
  //printf("  %%%% %s OPENED (WRITE)\n",OutNam);
  
  fprintf(OutFil, "NDIME= %d\n", Dim);

	if ( Msh->Dim == 2 ) {
		NbrElt = Msh->NbrTri+Msh->NbrQua;
	}
	else {
		NbrElt = Msh->NbrTet+Msh->NbrHex+Msh->NbrPri+Msh->NbrPyr;
	}

  fprintf(OutFil, "NELEM= %d\n", NbrElt);
	
	if ( Msh->Dim == 2 ){
		//--- Write triangles
		for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
  	  fprintf(OutFil, "%d ", SU2_TRIANGLE); 
			for (i=0; i<3; ++i) {
  	    fprintf(OutFil, "%d ",Msh->Tri[iTri][i]-1);
  	  }
  	  fprintf(OutFil, "%d\n", iTri-1); 
  	}
		
		//--- Write quads
		for (iQua=1; iQua<=Msh->NbrQua; iQua++) {
  	  fprintf(OutFil, "%d ", SU2_RECTANGLE); 
			for (i=0; i<4; ++i) {
  	    fprintf(OutFil, "%d ",Msh->Qua[iQua][i]-1);
  	  }
  	  fprintf(OutFil, "%d\n", iQua-1); 
  	}
		
	}
	
	
	for (iTet=1; iTet<=Msh->NbrTet; iTet++) {
    fprintf(OutFil, "%d ", SU2_TETRAHEDRAL); 
		for (i=0; i<4; ++i) {
      fprintf(OutFil, "%d ",Msh->Tet[iTet][i]-1);
    }
    fprintf(OutFil, "%d\n", iTet-1); 
  }
	
	for (i=1; i<=Msh->NbrHex; i++) {
    fprintf(OutFil, "%d ", SU2_HEXAHEDRAL); 
		for (j=0; j<8; ++j) {
      fprintf(OutFil, "%d ",Msh->Hex[i][j]-1);
    }
    fprintf(OutFil, "%d\n", i-1); 
  }
	
	for (i=1; i<=Msh->NbrPri; i++) {
    fprintf(OutFil, "%d ", SU2_WEDGE); 
		for (j=0; j<6; ++j) {
      fprintf(OutFil, "%d ",Msh->Pri[i][j]-1);
    }
    fprintf(OutFil, "%d\n", i-1); 
  }
	
	for (i=1; i<=Msh->NbrPyr; i++) {		
    fprintf(OutFil, "%d ", SU2_PYRAMID); 
		for (j=0; j<5; ++j) {
      fprintf(OutFil, "%d ",Msh->Pyr[i][j]-1);
    }
		
		if ( i== 1 )
			printf("PYR %d : %d %d %d %d %d\n", i, Msh->Pyr[i][0], Msh->Pyr[i][1], Msh->Pyr[i][2], Msh->Pyr[i][3], Msh->Pyr[i][4]);
    fprintf(OutFil, "%d\n", i-1); 
  }
	
  //--- Write vertices
  fprintf(OutFil, "NPOIN= %d\n", Msh->NbrVer);

	if ( Msh->Dim == 2 ){
  	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
    	fprintf(OutFil, "%.16le %.16le %d \n", Msh->Ver[iVer][0], Msh->Ver[iVer][1], iVer-1);
		}
  }
	else {
		for (iVer=1; iVer<=Msh->NbrVer; iVer++) {
    	fprintf(OutFil, "%.16le %.16le %.16le %d \n", Msh->Ver[iVer][0], Msh->Ver[iVer][1], Msh->Ver[iVer][2], iVer-1);
		}
	}

	//--- Write bdry elements
	
	if ( Msh->Dim == 2 ) {
		BdrTag = (int3*)malloc(sizeof(int3)*Msh->NbrEfr);
	  for (iEfr=1; iEfr<=Msh->NbrEfr; iEfr++) {
	    BdrTag[iEfr-1][0] = Msh->Efr[iEfr][2];
	    BdrTag[iEfr-1][1] = iEfr;
	  }  
	  NbrBdr = Msh->NbrEfr;
	}
	else {
		BdrTag = (int3*)malloc(sizeof(int3)*(Msh->NbrTri+Msh->NbrQua));
		
		idx = 0;
		
	  for (iTri=1; iTri<=Msh->NbrTri; iTri++) {
	    BdrTag[idx][0] = Msh->Tri[iTri][3];
	    BdrTag[idx][1] = iTri;
			BdrTag[idx][2] = SU2_TRIANGLE;
			idx++;
	  }  
	  for (iQua=1; iQua<=Msh->NbrQua; iQua++) {
	    BdrTag[idx][0] = Msh->Qua[iQua][4];
	    BdrTag[idx][1] = iQua;
			BdrTag[idx][2] = SU2_RECTANGLE;
			idx++;
	  }  
	  NbrBdr = Msh->NbrTri+Msh->NbrQua;
	}
	
	if ( NbrBdr > 0 ) {
		
		qsort(BdrTag, NbrBdr, sizeof(int3), cmp_int2 );
	  NbrTag = 1;
	
	  for (i=1; i<NbrBdr; i++) {
	    if ( BdrTag[i][0] != BdrTag[i-1][0] ) {
	      NbrTag++;
	    }
	  }
		
    fprintf(OutFil, "NMARK= %d\n", NbrTag);

    start = 0;
    for (iTag=0; iTag<NbrTag; iTag++) {
      fprintf(OutFil, "MARKER_TAG= %d\n", BdrTag[start][0]);
      
      cpt=1;
      for (i=start+1; i<NbrBdr; i++) {
        if ( BdrTag[i][0] != BdrTag[i-1][0] ) {
          break;
        }
        cpt++;
      }
      
      fprintf(OutFil, "MARKER_ELEMS= %d\n", cpt);
      
      for (i=start; i<start+cpt; i++) {
	
				if ( Msh->Dim == 2 ) {
					iEfr = BdrTag[i][1];
	      	fprintf(OutFil, "%d ", SU2_LINE);

	      	for (s=0; s<2; s++) {
	      	  iVer = Msh->Efr[iEfr][s]-1;
	      	  fprintf(OutFil, "%d ", iVer);
	      	}
	      	fprintf(OutFil, "\n");
				}
				else {
					
					if ( BdrTag[i][2] == SU2_TRIANGLE ) {
						iTri = BdrTag[i][1];
	      		fprintf(OutFil, "%d ", SU2_TRIANGLE);
	      		for (s=0; s<3; s++) {
	      		  iVer = Msh->Tri[iTri][s]-1;
	      		  fprintf(OutFil, "%d ", iVer);
	      		}
	      		fprintf(OutFil, "\n");
					}
					else if ( BdrTag[i][2] == SU2_RECTANGLE )  {
						iQua = BdrTag[i][1];
	      		fprintf(OutFil, "%d ", SU2_RECTANGLE);
	      		for (s=0; s<4; s++) {
	      		  iVer = Msh->Qua[iQua][s]-1;
	      		  fprintf(OutFil, "%d ", iVer);
	      		}
	      		fprintf(OutFil, "\n");
					}
			}
      		
      	  
     }
      
      start = start+cpt;
      
    }//for iTag 
    
	}
  else
		fprintf(OutFil, "NMARK= 0\n");
	
  //--- close mesh file
	if ( OutFil )
		fclose(OutFil);
		
	if ( BdrTag )
		free(BdrTag);
  
  return;
}



int WriteSU2Solution (char *SolNam, Mesh *Msh, double *Sol, int NbrVer, int SolSiz, char SolTag[100][256])
{
  int       i, s, d;
  int       iVer,idxVer;
	
	int2 *BdrTag=NULL;
	int NbrBdr, NbrTag, start, iTag, cpt;
	
	FILE *OutFil=NULL; 
	OutFil = fopen(SolNam, "wb");
	
	if ( !OutFil ) {
		printf("  ## ERROR WriteSU2Solution: Can't open %s\n", SolNam);
	}
 	
  //printf("  %%%% %s OPENED (WRITE)\n",SolNam);
	
	//--- Write header
	
	fprintf(OutFil,"\"PointID\"	\"x\"	\"y\"	");
	if ( Msh->Dim == 3 ) fprintf(OutFil,"\"z\"	");	
	for (i=0; i<SolSiz; i++) {
		fprintf(OutFil, "\"%s\"	", SolTag[i]);
	}
	fprintf(OutFil, "\n");
	
	//--- Write solution at vertices
  
	for (iVer=1; iVer<=NbrVer; iVer++) {
		fprintf(OutFil, "%d	", iVer-1);
		for (d=0; d<Msh->Dim; d++)
			fprintf(OutFil, "%.15le	", Msh->Ver[iVer][d]);
		
		idxVer = iVer*SolSiz;
		for (i=0; i<SolSiz; i++) {
			fprintf(OutFil, "%.15le	", Sol[idxVer+i]);
		}
		fprintf(OutFil, "\n");
	}
  
  //--- close mesh file
	if ( OutFil )
		fclose(OutFil);
	
  
  return 1;
}


