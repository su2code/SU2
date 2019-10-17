#include "amgio.h"


int ConvertGMFtoSU2Sol (Options *mshopt)
{
	Mesh *Msh = NULL;
	char OutSol[1024];

	FILE *FilHdl=NULL;
	char *tok=NULL, *lin=NULL;

	int skip=0, SolSiz=0;
	size_t  len = 0;

	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);

	if ( Msh->FilTyp != FILE_GMF ) {
		printf("  ## ERROR : Input mesh file must be a .mesh (GMF) (FilTyp=%d)\n", Msh->FilTyp);
		return 0;
	}

	//PrintMeshInfo (Msh);
	
	//--- Get header information from reference file

	if ( strcmp(mshopt->HeaderNam, "") ) {

		FilHdl = fopen(mshopt->HeaderNam,"r");

		if ( !FilHdl ) {
			printf("  ## ERROR: ConvertGMFtoSU2Sol: UNABLE TO OPEN %s. \n", mshopt->HeaderNam);
	    return 0;
		}

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
				StrRemoveChars(Msh->SolTag[SolSiz], '\"');
				StrRemoveChars(Msh->SolTag[SolSiz], '\n');
				SolSiz++;

				if ( SolSiz > Msh->SolSiz ) {
					printf("  ## ERROR: ConvertGMFtoSU2Sol: Provided header size does not match.\n");
			    return 0;
				}

				tok = strtok (NULL, "	,");
			}
	  }

	}

	if ( mshopt->clean == 1 )
		RemoveUnconnectedVertices(Msh);
	
	if ( mshopt->Dim == 2 )
		Msh->Dim = 2;
	
	WriteSU2Mesh(mshopt->OutNam, Msh);
	
	if ( Msh->Sol ) {
		sprintf(OutSol, "%s.csv", mshopt->OutNam);
		WriteSU2Solution (OutSol, Msh->Dim, Msh->NbrVer, Msh->Ver, Msh->Sol, Msh->SolSiz, Msh->SolTag);
	}

	if ( Msh )
 		FreeMesh(Msh);

	return 1;
}



int ConvertSU2SolToGMF (Options *mshopt)
{
	Mesh *Msh = NULL;
	char OutSol[1024];

	Msh = SetupMeshAndSolution (mshopt->InpNam, mshopt->SolNam);

	if ( Msh->FilTyp != FILE_SU2 ) {
		printf("  ## ERROR : Input mesh file must be a .su2.\n");
		return 0;
	}

	//PrintMeshInfo (Msh);

	WriteGMFMesh(mshopt->OutNam, Msh, 1);

	if ( Msh->Sol ) {
		sprintf(OutSol, "%s.solb", mshopt->OutNam);
		if ( ! WriteGMFSolutionItf(OutSol, Msh) ) {
			printf("  ## ERROR : outputmach FAILED.\n");
		}
	}

	if ( Msh )
 		FreeMesh(Msh);

	return 1;
}



int SplitSolution (Mesh *Msh, char *prefix, char *adap_sensor)
{
	int NbrFld = 1, i, iVer, idx;
	int FldTab[10]; 
	double *OutSol = NULL;
	char OutNam[256];
	
	int pres_flag=0, mach_flag=0, temp_flag=0, goal_flag=0;
		
	if (!strcmp(adap_sensor, "MACH_PRES")) {
		NbrFld = 2;
		pres_flag = mach_flag = 1;
	}
	else if (!strcmp(adap_sensor, "MACH")) {
		NbrFld = 1;
		mach_flag = 1;
	}
	else if (!strcmp(adap_sensor, "PRES")) {
		NbrFld = 1;
		pres_flag = 1;
	}
	else if (!strcmp(adap_sensor, "GOAL_ECC")) {
		NbrFld = 2;
		goal_flag = mach_flag = 1;
	}
	else {
		printf("## ERROR SplitSolution: Unknown adap_sensor.\n");
		exit(1);
	}
	
	for (i=0; i<NbrFld; i++) {
		FldTab[i] = 1;
	}
	
	OutSol = (double*)malloc(sizeof(double)*(Msh->NbrVer+1)*NbrFld);
	
	//--- Get fields indices

	int iMach = -1;
	int iPres = -1;
	int iTemp = -1;
	int iGoal = -1;
		
	for (i=0; i<Msh->NbrFld; i++) {
		if ( !strcmp(Msh->SolTag[i], "Mach") && mach_flag == 1 ) {
			iMach = i;
		}
		if ( !strcmp(Msh->SolTag[i], "Pressure") && pres_flag == 1 ) {
			iPres = i;
		}
		if ( !strcmp(Msh->SolTag[i], "Temperature") && temp_flag == 1 ) {
			iTemp = i;
		}
		if ( !strcmp(Msh->SolTag[i], "Adaptation_Parameter") && goal_flag == 1 ) {
			iGoal = i;
		}
	}
	
	if ( iMach < 0 ) {
		printf("  ## ERROR OutputMach : Mach index not found.\n");
		return 0;
	}
	
	if ( iPres < 0 ) {
		printf("  ## ERROR OutputMach : Pres index not found.\n");
		return 0;
	}
	
	for (iVer=1; iVer<=Msh->NbrVer; iVer++) {		
		
		idx = iVer*NbrFld-1;
		
		if ( mach_flag == 1 ){
			idx++;
			OutSol[idx] = Msh->Sol[iVer*Msh->SolSiz+iMach];
		}
		
		if ( pres_flag == 1 ){
			idx++;
			OutSol[idx] = Msh->Sol[iVer*Msh->SolSiz+iPres];
		}

		if ( goal_flag == 1 ){
			idx++;
			OutSol[idx] = Msh->Sol[iVer*Msh->SolSiz+iGoal];
		}
	}
	
	sprintf(OutNam, "%s_sensor.solb", prefix);
	if ( ! WriteGMFSolution(OutNam, OutSol, NbrFld, Msh->NbrVer, Msh->Dim, NbrFld, FldTab) ) {
		printf("  ## ERROR : Output of solution failed.\n");
	}
	
	if ( OutSol )
		free(OutSol);
	
	if ( Msh )
 		FreeMesh(Msh);
	
	return 1;
}
