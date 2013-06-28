#ifndef _CALC_H_
#define _CALC_H_

#include <stdio.h>
#include "vecsym.h"

extern int cgnsFile;

/*--- base data ---*/

extern int NumBases, cgnsBase;
extern char BaseName[33];
extern int CellDim, PhyDim;
extern int BaseClass, BaseUnits[5];

/*--- zone data ---*/

extern int NumZones, cgnsZone;
extern char ZoneName[33];
extern int ZoneType, ZoneDims[6];
extern int ZoneClass, ZoneUnits[5];
extern int GridClass, GridUnits[5];

/*--- solution data ---*/

extern int NumSolns, cgnsSoln;
extern char SolnName[33];
extern int SolnLocation;
extern int SolnDims[3], SolnRind[6];
extern int SolnClass, SolnUnits[5];

/*--- solution field data ---*/

typedef struct {
    int units[5];
    float exps[5];
} Units;

typedef struct {
    char name[33];
    int type;
    int id;
    int valid;
    int len;
    int datatype;
    int dataclass;
    int hasunits;
    int units[5];
    int hasconv;
    double dataconv[2];
    int hasexp;
    double exponent[5];
    VECDATA *vd;
} Variable;

/*----- variables -----*/

extern Variable *variables;
extern int NumVariables;

/*----- reference conditions -----*/

extern Variable *reference;
extern int NumReference;

/*----- mesh -----*/

extern Variable *coordinates;
extern int NumCoordinates;

/*----- external functions -----*/

#ifdef EXTERN_FUNCS
extern void add_funcs (
    void
);
#endif

/*----- functions -----*/

void cgnsCalcFatal (     /* terminate with error message */
    char *errmsg         /* error message */
);

void cgnsCalcError (     /* print error message */
    char *errmsg         /* error message */
);

void cgnsCalcReset (     /* reset calculator (symbol table) */
    void
);

int cgnsCalcInit (       /* load CGNS file and initialize */
    char *cgnsfile,      /* CGNS file */
    int modify,          /* set for modify mode */
    void (*errhandler)(  /* calculator error callback */
        int errnum,      /* error number */
        char *errmsg,    /* error message */
        int pos,         /* location in string */
        char *str        /* string being parsed */
    )
);

void cgnsCalcDone (      /* close CGNS file */
    void
);

int cgnsCalcBase (       /* set base for calculations */
    int base             /* base number */
);

int cgnsCalcZone (       /* set zone for calculations */
    int zone             /* zone number */
);

int cgnsCalcSoln (       /* set solution for calculations */
    int soln             /* solution number */
);

int cgnsCalcCheck (      /* parse command and check for errors */
    char *expression     /* expression to be parsed */
);

VECSYM *cgnsCalcCommand (/* parse command string and return results */
    char *expression     /* expression to be parsed */
);

Variable *cgnsCalcVarGet (/* return a variable */
    char *varname
);

void cgnsCalcVarList (   /* print variables */
    FILE *fp             /* output file (NULL gives stdout) */
);

void cgnsCalcSymList (   /* print list of symbols */
    FILE *fp             /* output file (NULL gives stdout) */
);

#endif  /* _CALC_H_ */

