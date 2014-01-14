/*
*****************************************************************
*****************************************************************
*******                                                  ********
****** (C) Copyright 2004-2006  by Tecplot, Inc.         ********
****** (C) Copyright 1989-2003  by AMTEC ENGINEERING INC.********
*******       All Rights Reserved.                       ********
*******                                                  ********
*****************************************************************
*****************************************************************
*/


/*
****************************************************************
****************** BEGIN DEVELOPMENT NOTES *********************
****************************************************************

BEGIN CODELOG PLTVIEW
V 09/04/98
V ****************************************************************
V *                 Build 1.0 9-04-98                            *
V ****************************************************************
END CODELOG

*********************************************************************
* IMPORTANT NOTE: Only development notes for "pltview" stand-alone  *
*                 belong in this file. See "ADDONVER.h" for changes *
*                 related to the add-on.                            *
*********************************************************************

****************************************************************
*  V in column 1 marks date information.                       *
*  C in column 1 marks notes on new changes.                   *
*  B in column 1 marks notes on bug fixes.                     *
****************************************************************

****************************************************************
****************** END DEVELOPMENT NOTES ***********************
****************************************************************
*/

#if defined ADDON
#include "TECADDON.h"
#include "GUIDEFS.h"
#include "TECGUI.h"
#define READTEC                TecUtilReadBinaryData
#define SHOWINFO(S)            TecGUITextAppendString(Output_T_D1,S);
#define ERRMSG(S)              TecUtilDialogErrMsg(S)
#define ALLOC_ARRAY(N,Type,S)  (Type *)TecUtilStringAlloc((N)*sizeof(Type),"debug info")
#define FREE_ARRAY(N,S)           TecUtilStringDealloc((char **)&N)
#define STRINGLISTGETSTRING(S,N)  TecUtilStringListGetString(S,N)
#define STRINGLISTGETCOUNT(S)     TecUtilStringListGetCount(S)
#define STRINGLISTDEALLOC(S)      TecUtilStringListDealloc(S)

#define AUXDATADEALLOC(S)         TecUtilAuxDataDealloc(S)
#define AUXDATAGETNUMITEMS(S)     TecUtilAuxDataGetNumItems(S)
#define AUXDATAGETITEMBYINDEX(S,Index,Name,Value,Type,Retain) \
    TecUtilAuxDataGetItemByIndex(S,(Index),Name,Value,Type,Retain)

#else
#include "MASTER.h"
#include "GLOBAL.h"
#include "TASSERT.h"
#define ALLOC_ARRAY(N,Type,S)     (Type *)TecAlloc((N)*sizeof(Type))
#define FREE_ARRAY(N,S)           TecFree((void *)N)
#include "ARRLIST.h"
#include "STRLIST.h"
#include "AUXDATA.h"
#include "DATAUTIL.h"

#define READTEC                   ReadTec
#define SHOWINFO(S)               printf("%s",S);
#define ERRMSG(S)                 printf("Err: %s\n",S);
#define STRINGLISTGETSTRING(S,N)  StringListGetString(S,(N)-1)
#define STRINGLISTGETCOUNT(S)     StringListCount(S)
#define STRINGLISTDEALLOC(S)      StringListDealloc(S)

#define AUXDATADEALLOC(S)         AuxDataDealloc(S)
#define AUXDATAGETNUMITEMS(S)     AuxDataGetNumItems(S)
#define AUXDATAGETITEMBYINDEX(S,Index,Name,Value,Type,Retain) \
    AuxDataGetItemByIndex(S,(Index)-1,Name,Value,Type,Retain)

#endif


static int GetNumPtsPerCell(ZoneType_e ZoneType)
{
    int NumPts = 0;
    switch (ZoneType)
    {
        case ZoneType_FETriangle : NumPts = 3; break;
        case ZoneType_FEQuad     : NumPts = 4; break;
        case ZoneType_FETetra    : NumPts = 4; break;
        case ZoneType_FEBrick    : NumPts = 8; break;
        default : NumPts = 0;
    }
    return (NumPts);
}

static int GetNumPts(ZoneType_e ZoneType,
                     LgIndex_t  NumPtsI,
                     LgIndex_t  NumPtsJ,
                     LgIndex_t  NumPtsK)
{
    int NumPts = 0;
    switch (ZoneType)
    {
        case ZoneType_FETriangle  :
        case ZoneType_FEQuad      :
        case ZoneType_FETetra     :
        case ZoneType_FEBrick     :
        case ZoneType_FEPolygon   :
        case ZoneType_FEPolyhedron: NumPts = NumPtsI; break;
        default : NumPts = NumPtsI*NumPtsJ*NumPtsK;
    }
    return (NumPts);
}



static void DeallocHeaderInfo(char         **DataSetTitle,
                              StringList_pa *VarNames,
                              StringList_pa *ZoneNames,
                              LgIndex_t    **NumPtsI,
                              LgIndex_t    **NumPtsJ,
                              LgIndex_t    **NumPtsK,
                              ZoneType_e   **ZoneType,
                              StringList_pa *UserRec,
                              AuxData_pa    *DatasetAuxData)
{
    if (*DataSetTitle)
        FREE_ARRAY(*DataSetTitle, "data set title");
    if (*VarNames)
        STRINGLISTDEALLOC(VarNames);
    if (*ZoneNames)
        STRINGLISTDEALLOC(ZoneNames);
    if (*NumPtsI)
        FREE_ARRAY(*NumPtsI, "NumPtsI Array");
    if (*NumPtsJ)
        FREE_ARRAY(*NumPtsJ, "NumPtsJ Array");
    if (*NumPtsK)
        FREE_ARRAY(*NumPtsK, "NumPtsK Array");
    if (*ZoneType)
        FREE_ARRAY(*ZoneType, "ZoneType Array");
    if (*UserRec)
        STRINGLISTDEALLOC(UserRec);
    if (*DatasetAuxData)
        AUXDATADEALLOC(DatasetAuxData);

    *DataSetTitle  = NULL;
    *VarNames      = NULL;
    *ZoneNames     = NULL;
    *NumPtsI       = NULL;
    *NumPtsJ       = NULL;
    *NumPtsK       = NULL;
    *ZoneType      = NULL;
    *UserRec       = NULL;
}

#define MAXCHARSINFOLINE 5000

void ShowDatasetAuxData(AuxData_pa DatasetAuxData, char* InfoLine) 
{
    REQUIRE(VALID_REF(InfoLine));
    if (DatasetAuxData != NULL)
    {
        SHOWINFO("Dataset Auxiliary Data:\n");
        SHOWINFO("-----------------------\n");
        LgIndex_t NumItems = AUXDATAGETNUMITEMS(DatasetAuxData);
        for (LgIndex_t ii = 0; ii < NumItems; ++ii)
        {
#if defined ADDON
            char* Name;
#else
            const char* Name;
#endif
            ArbParam_t Value;
            AuxDataType_e Type;
            Boolean_t Retain;
            AUXDATAGETITEMBYINDEX(DatasetAuxData,
                ii+1,
                &Name,
                &Value,
                &Type,
                &Retain);
            sprintf(InfoLine, "  Name: %s\n", Name);
            SHOWINFO(InfoLine);
#if defined ADDON
            // The TecUtil layer returns copies which must be deallocated
            TecUtilStringDealloc(&Name);
#endif
            if (Type == AuxDataType_String)
            {
                char* ValueString = reinterpret_cast<char*>(Value);
                sprintf(InfoLine, "    Value : %s\n", ValueString);
                SHOWINFO(InfoLine);
                SHOWINFO("    Type  : String\n");
#if defined ADDON
                // The TecUtil layer returns copies which must be deallocated
                TecUtilStringDealloc(&ValueString);
#endif
            }
            sprintf(InfoLine, "    Retain: %s\n", Retain ? "True" : "False");
            SHOWINFO(InfoLine);
        }
        SHOWINFO("\n");
    }
}


void ReportFileInfo(char     *FName,
                    Boolean_t LoadRawData,
                    Boolean_t AllocateRawDataSpaceLocally)
{
    short          IVersion;
    EntIndex_t     NumZones;
    EntIndex_t     NumVars;
    char          *DataSetTitle   = NULL;
    StringList_pa  VarNames       = NULL;
    StringList_pa  ZoneNames      = NULL;
    LgIndex_t     *NumPtsI        = NULL;
    LgIndex_t     *NumPtsJ        = NULL;
    LgIndex_t     *NumPtsK        = NULL;
    ZoneType_e    *ZoneType       = NULL;
    StringList_pa  UserRec        = NULL;
    AuxData_pa     DatasetAuxData = NULL;
    int            CZ, CV;
    char           InfoLine[MAXCHARSINFOLINE+1];
    double       **VDataBase    = NULL;
    NodeMap_t    **NodeMap      = NULL;

    /*
     * Load in the header information only.
     */

    if (!READTEC(TRUE,
                 FName,
                 &IVersion,
                 &DataSetTitle,
                 &NumZones,
                 &NumVars,
                 &VarNames,
                 &ZoneNames,
                 &NumPtsI,
                 &NumPtsJ,
                 &NumPtsK,
                 &ZoneType,
                 &UserRec,
                 &DatasetAuxData,
                 FALSE,
                 (NodeMap_t ***)NULL,
                 (double ***)NULL))
    {
        sprintf(InfoLine, "Cannot read file \"%s\"\nor file is not a Tecplot binary data file.\n", FName);
        ERRMSG(InfoLine);
    }
    else
    {
        Boolean_t  IsOk = TRUE;
        if (LoadRawData)
        {
            if (AllocateRawDataSpaceLocally)
            {
                int NumPts;
                VDataBase = ALLOC_ARRAY(NumZones * NumVars, double *, "vdatabase array");
                for (CZ = 0; CZ < NumZones; CZ++)
                    for (CV = 0; CV < NumVars; CV++)
                    {
                        NumPts = GetNumPts(ZoneType[CZ], NumPtsI[CZ], NumPtsJ[CZ], NumPtsK[CZ]);
                        if (NumPts >= 1)
                            VDataBase[CZ*NumVars+CV] = ALLOC_ARRAY(NumPts, double, "vdatabase array");
                        else
                            VDataBase[CZ*NumVars+CV] = NULL;
                    }

                NodeMap = ALLOC_ARRAY(NumZones, NodeMap_t *, "nodemap array");
                for (CZ = 0; CZ < NumZones; CZ++)
                {
                    if (ZoneType[CZ] == ZoneType_Ordered)
                        NodeMap[CZ] = NULL;
                    else
                    {
                        int PtsPerCell = GetNumPtsPerCell(ZoneType[CZ]);
                        NodeMap[CZ] = ALLOC_ARRAY(PtsPerCell * NumPtsJ[CZ],
                                                  NodeMap_t, "zone nodemap");
                    }
                }
            }
            else
            {
                VDataBase = NULL;
                NodeMap   = NULL;
            }

            /*
             * NOTE: If any of the above alloc's failed then no big deal.  ReadTec
             *       itself "skips" vars if memory was not allocated for it.
             */

            DeallocHeaderInfo(&DataSetTitle,
                              &VarNames,
                              &ZoneNames,
                              &NumPtsI,
                              &NumPtsJ,
                              &NumPtsK,
                              &ZoneType,
                              &UserRec,
                              &DatasetAuxData);

            /*
             * Reread the datafile.  This time load in the header AND the raw data
             * Note that VDataBase may be preallocated or may be left up to ReadTec
             * to allocate (See AllocateRawDataSpaceLocally above).
             */
            if (!READTEC(FALSE,
                         FName,
                         &IVersion,
                         &DataSetTitle,
                         &NumZones,
                         &NumVars,
                         &VarNames,
                         &ZoneNames,
                         &NumPtsI,
                         &NumPtsJ,
                         &NumPtsK,
                         &ZoneType,
                         &UserRec,
                         &DatasetAuxData,
                         AllocateRawDataSpaceLocally,
                         &NodeMap,
                         &VDataBase))
            {
                if (IVersion > 99)
                {
                    sprintf(InfoLine,
                            "Error: ***\n"
                            "    This add-on can only display raw nodal data\n"
                            "    and it appears to contain cell centered data.\n");
                    SHOWINFO(InfoLine);
                }
                else
                {
                    sprintf(InfoLine,
                            "Cannot Read File, %s.\n"
                            "File may not be a tecplot binary data file.",
                            FName);
                    ERRMSG(InfoLine);
                }
                IsOk = FALSE;
            }
        }

        SHOWINFO("\n");
        sprintf(InfoLine, "FileName    : %s\n", FName);
        SHOWINFO(InfoLine);
        sprintf(InfoLine, "File Version: %3.1f\n", IVersion / 10.0);
        SHOWINFO(InfoLine);

        /* if the file contains filetype, then retrieve that separately since ReadTec should not be changed */
        if (IVersion >= 109)
        {
            DataFileType_e FileType = DataFileType_Full;
            char           FileTypeStr[32];
            FILE          *F = NULL;

            /* open the file and get the filetype */
            F = fopen(FName, "rb");
            if (F)
            {
                char    Buffer[8];
                Int32_t One;
                Int32_t FileTypeInt;

                /* 8 bytes for magic# and version and 4 bytes for Int32 */
                fread(Buffer, sizeof(Buffer[0]), 8, F);
                fread(&One, sizeof(One), 1, F);
                fread(&FileTypeInt, sizeof(FileTypeInt), 1, F);
                FileType = (DataFileType_e)FileTypeInt;
                fclose(F);
                F = NULL;
            }
            /* map the filetype */
            switch (FileType)
            {
                case DataFileType_Full:
                    strcpy(FileTypeStr, "Full");
                    break;
                case DataFileType_Grid:
                    strcpy(FileTypeStr, "Grid");
                    break;
                case DataFileType_Solution:
                    strcpy(FileTypeStr, "Solution");
                    break;
                default:
                    IsOk = FALSE;
                    CHECK(FALSE);
                    break;
            }
            sprintf(InfoLine, "File Type   : %s\n", FileTypeStr);
            SHOWINFO(InfoLine);
        }

        sprintf(InfoLine, "DataSetTitle: %s\n", DataSetTitle ? DataSetTitle : " ");
        SHOWINFO(InfoLine);

        ShowDatasetAuxData(DatasetAuxData, InfoLine);

        sprintf(InfoLine, "NumZones    : %d\n", (int)NumZones);
        SHOWINFO(InfoLine);
        sprintf(InfoLine, "NumVars     : %d\n", (int)NumVars);
        SHOWINFO(InfoLine);
        if (IsOk && (NumZones > 0))
        {
            SHOWINFO("Var Names   : ");
            for (CZ = 0; CZ < NumVars; CZ++)
            {
                char *VarName = STRINGLISTGETSTRING(VarNames, CZ + 1);
                sprintf(InfoLine, "%s", VarName ? VarName : "NULL");
                if (CZ < NumVars - 1)
                    strcat(InfoLine, ",");
                else
                    strcat(InfoLine, "\n\n");
                SHOWINFO(InfoLine);
                if (VarName)
                    FREE_ARRAY(VarName, "VarName array");
            }
            SHOWINFO("ZoneName            IMax    JMax    KMax    Node     Face     Elmt     EType\n");
            SHOWINFO("-------------------------------------------------------------------------------\n");

            for (CZ = 0; CZ < NumZones; CZ++)
            {
                char *ZoneName = STRINGLISTGETSTRING(ZoneNames, CZ + 1);
                if (ZoneType[CZ] != ZoneType_Ordered)
                {
                    if (ZoneType[CZ] == ZoneType_FEPolygon ||
                        ZoneType[CZ] == ZoneType_FEPolyhedron)
                        sprintf(InfoLine, "%-20s ---     ---     ---     %-8ld %-8ld %-8ld ",
                                (ZoneName ? ZoneName : "NULL"),
                                (long)NumPtsI[CZ],
                                (long)NumPtsK[CZ],
                                (long)NumPtsJ[CZ]);
                    else
                        sprintf(InfoLine, "%-20s ---     ---     ---     %-8ld ---      %-8ld ",
                                (ZoneName ? ZoneName : "NULL"),
                                (long)NumPtsI[CZ],
                                (long)NumPtsJ[CZ]);
                    SHOWINFO(InfoLine);
                    switch (ZoneType[CZ])
                    {
                        case ZoneType_FETriangle  : SHOWINFO("Tri\n");    break;
                        case ZoneType_FEQuad      : SHOWINFO("Quad\n");   break;
                        case ZoneType_FETetra     : SHOWINFO("Tetra\n");  break;
                        case ZoneType_FEBrick     : SHOWINFO("Brick\n");  break;
                        case ZoneType_FELineSeg   : SHOWINFO("LineSeg\n"); break;
                        case ZoneType_FEPolygon   : SHOWINFO("Polygon\n"); break;
                        case ZoneType_FEPolyhedron: SHOWINFO("Polyhed\n"); break;
                        default: CHECK(FALSE); break;
                    }
                }
                else
                {
                    sprintf(InfoLine, "%-20s %-7ld %-7ld %-7ld ---      ---      ---      ---\n",
                            (ZoneName ? ZoneName : "NULL"),
                            (long)NumPtsI[CZ],
                            (long)NumPtsJ[CZ],
                            (long)NumPtsK[CZ]);
                    SHOWINFO(InfoLine);
                }
                if (ZoneName)
                    FREE_ARRAY(ZoneName, "ZoneName Array");
            }
        }

        if (IsOk)
        {
            for (CZ = 1; CZ <= STRINGLISTGETCOUNT(UserRec); CZ++)
            {
                char *S = STRINGLISTGETSTRING(UserRec, CZ);
                if (S)
                {
                    int L;
                    strcpy(InfoLine, "UserRec: ");
                    L = (int)strlen(InfoLine);
                    strncat(&InfoLine[L], S, MAXCHARSINFOLINE - L - 2);
                    strcat(&InfoLine[strlen(InfoLine)], "\n");
                    SHOWINFO(InfoLine);
                    FREE_ARRAY(S, "temp string");
                }
            }


            if (LoadRawData)
            {
                for (CZ = 0; CZ < NumZones; CZ++)
                {
                    int CV;
                    for (CV = 0; CV < NumVars; CV++)
                        if (VDataBase[CZ*NumVars+CV])
                        {
                            int I;
                            int NumPts = GetNumPts(ZoneType[CZ], NumPtsI[CZ], NumPtsJ[CZ], NumPtsK[CZ]);
                            int  SLen = 0;
                            sprintf(InfoLine, "\n\nVariable data for zone %d, Var %d\n", CZ + 1, CV + 1);
                            SHOWINFO(InfoLine);
                            InfoLine[0] = '\0';
                            for (I = 0; I < NumPts; I++)
                            {
                                char PString[50];
                                if (SLen + 50 > MAXCHARSINFOLINE)
                                {
                                    SHOWINFO(InfoLine);
                                    InfoLine[0] = '\0';
                                    SLen = 0;
                                }

                                sprintf(PString, "%lG ", VDataBase[CZ*NumVars+CV][I]);
                                strcat(InfoLine, PString);
                                SLen += (int)strlen(PString);

                                if ((I % 5) == 4)
                                {
                                    strcat(InfoLine, "\n");
                                    SLen++;
                                }
                            }
                            if (*InfoLine)
                                SHOWINFO(InfoLine);
                            FREE_ARRAY(VDataBase[CZ*NumVars+CV], "vdatabase double");
                        }
                    if (NodeMap[CZ])
                    {
                        int I, J;
                        int PtsPerCell = GetNumPtsPerCell(ZoneType[CZ]);
                        int SLen = 0;
                        SHOWINFO("\nConnectivity list:\n");
                        InfoLine[0] = '\0';
                        for (J = 0; J < NumPtsJ[CZ]; J++)
                        {
                            if (SLen + 200 > MAXCHARSINFOLINE)
                            {
                                SHOWINFO(InfoLine);
                                InfoLine[0] = '\0';
                                SLen = 0;
                            }
                            for (I = 0; I < PtsPerCell; I++)
                            {
                                char NString[20];
                                sprintf(NString, "%u ", (unsigned int)NodeMap[CZ][J*PtsPerCell+I] + 1);
                                strcat(InfoLine, NString);
                                SLen += (int)strlen(NString);
                            }
                            strcat(InfoLine, "\n");
                            SLen++;
                        }
                        if (*InfoLine)
                            SHOWINFO(InfoLine);
                        FREE_ARRAY(NodeMap[CZ], "nodemap");
                    }
                }
                FREE_ARRAY(NodeMap, "Nodemap base array");
                FREE_ARRAY(VDataBase, "vdatabase base array");
            }
        }

        SHOWINFO("\n\n");

        DeallocHeaderInfo(&DataSetTitle,
                          &VarNames,
                          &ZoneNames,
                          &NumPtsI,
                          &NumPtsJ,
                          &NumPtsK,
                          &ZoneType,
                          &UserRec,
                          &DatasetAuxData);
    }
}



#if !defined ADDON
int main(int argc, char *(argv[]))
{
    short CurFile;

    if (argc == 1)
    {
        printf("Err: Need:  pltview file1 [file2] ...\n");
        exit(-1);
    }

    for (CurFile = 1; CurFile < argc; CurFile++)
        ReportFileInfo(argv[CurFile], FALSE, FALSE);

    return 0;
}
#endif
