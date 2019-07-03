#include "amgio.h"

/*
Victorien Menier Feb 2016
*/


Options* AllocOptions()
{
  Options *mshopt = (Options*)malloc(sizeof(struct T_Options));

  mshopt->Mod       = -1;
	mshopt->InpFilTyp = 0;
	
	mshopt->flagSol = 0;
	
	mshopt->clean = 0;
	
	strcpy(mshopt->InpNam, "");
	strcpy(mshopt->OutNam, "");
	strcpy(mshopt->BasNam, "");
	strcpy(mshopt->SolNam, "");
	
	strcpy(mshopt->HeaderNam, "");
	
  return mshopt;
}


int CheckOptions (Options *mshopt)
{
	int SolTyp=-1;
	
	if ( !strcmp(mshopt->InpNam, "") ) {
		printf("  ## ERROR CheckOptions : An input name must be provided.\n");
		return 0;
	}
	
	//--- Get mesh file extension
	mshopt->InpFilTyp = GetInputFileType(mshopt->InpNam);
	
	if ( !mshopt->InpFilTyp ) {
		printf("  ## ERROR CheckOptions : Unknown input file extension.\n");
		return 0;
	}
		
	if ( !strcmp(mshopt->OutNam, "") ) {
		GetBasNam (mshopt->InpNam, mshopt->BasNam);
		sprintf(mshopt->OutNam, "%s.o", mshopt->BasNam);
		printf("  -- Info: Output file name set to %s\n", mshopt->OutNam);
	}
	
	
	if ( strcmp(mshopt->SolNam, "") ) {
		SolTyp = GetInputFileType(mshopt->SolNam);
		
		if ( mshopt->InpFilTyp == FILE_SU2 && SolTyp != FILE_DAT ) {
			printf("  ## ERROR : Wrong format for solution file (.dat expected)\n");
			return 0;
		}
		
		if ( mshopt->InpFilTyp == FILE_GMF && SolTyp != FILE_GMFSOL ) {
			printf("  ## ERROR : Wrong format for solution file (.sol[b] expected)\n");
			return 0;
		}
		
	}
	
	
	return 1;
}

void PrintOptions (Options *mshopt)
{
	
	printf("\n--- Options Summary ---\n");
	printf("  Input mesh file : %s\n", mshopt->InpNam);
	printf("  Input sol  file : %s\n", (!strcmp(mshopt->SolNam,"")?"(None)":mshopt->SolNam));
	printf("  Output file     : %s\n", mshopt->OutNam);
	printf("-------------------------\n");
	
}


int GetBasNam (char *InpNam, char *BasNam)
{
	
  char *ptr = NULL;
  
  if ( !InpNam || !BasNam ) {
    printf("  ## ERROR GetBasNam.\n");
		return 0;
  }
	
  strcpy(BasNam,InpNam);
	
  ptr = strstr(BasNam,".mesh");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';
    
  ptr = strstr(BasNam,".sol");	
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';

  ptr = strstr(BasNam,".solb");
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';

  ptr = strstr(BasNam,".su2");
  if ( ptr != NULL )
    BasNam[ptr-BasNam]='\0';

	return 1;	
}