/*
Victorien Menier Feb 2016
*/

typedef struct T_Options
{
	
	int  Mod;              // module
	
	int  InpFilTyp; // input file's type: .mesh, .su2, etc.
	
	char InpNam[1024];   // Input mesh file name
	char BasNam[1024];   // Base name from InpNam
	char OutNam[1024];   // Output file name
	char SolNam[1024];   // Input solution name
	
	char HeaderNam[1024]; // SU2 solution file name used to copy the header information
	
	int flagSol;  // Output solution?
	
	int clean; // clean mesh?
	int Dim; // force dim to be 2?
	
	double Box[6];
	
} Options;

Options* AllocOptions(void);
