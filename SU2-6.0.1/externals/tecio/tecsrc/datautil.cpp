#include "stdafx.h"
#include "MASTER.h"
#define TECPLOTENGINEMODULE
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

/*
 * datautil.c:
 *
 * version 1.00 : 12/10/91 (cm) changes made for manual
 * version 1.01 : 12/30/91 Get and ReturnHugeBlock are now ptr to function
 * version 6.00 : 04/21/92 updated to match version 6 of tecplot.
 * version 6.30 : 10/15/92 updated to match binary file version 6.3
 * version 6.30a: 05/04/93 (cm) minor changes to prototypes
 * version      : 11/01/93 (cm) put in D4GW stuff
 * version 6.30b: 12/27/93 (cm) fixed missing NumKPts in DumpZone
 * version 6.30c: 12/27/93 (cm) put back in D4GW stuff
BEGIN CODELOG TECXXX
C 03/06/96 (BDP)
C   Update to V7
C
C 03/14/97 (BDP)
C   Added code to main tecplot source.  Now can
C   be built stand alone or added so TecUtil_ functions
C   can access.
C 06/02/98 (bdp)
C   v75 coding.  Also removed Array of ZoneSpec_s
C   structs in favor of zonenames, i,j,and k dimensions
C   and zonetype array.
END CODELOG
 */



#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "SYSTEM.h"
#include "ALLOC.h"
#include "TECXXX.h"
#include "ARRLIST.h"
#include "SET.h"
#include "DATASET.h"
#include "FILESTREAM.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "DATAIO.h"
#include "DATAIO4.h"
#include "DATAUTIL.h"
#include "STRLIST.h"
#include "Q_MSG.h"
#if defined MAKEARCHIVE
#define INITMODULE
#endif
#include "INPUT.h"

using namespace tecplot::strutil;

#if defined MAKEARCHIVE
#define ANGLEEPSILON 1.0e-10

void InitInputSpecs(void)
{
    LineThicknessInputSpec.Type                                = Input_Double;
    LineThicknessInputSpec.Min                                 = 0.000001;
    LineThicknessInputSpec.Max                                 = 1.0;
    LineThicknessInputSpec.InterfaceAdjust.ScaleFact           = 100.0;
    LineThicknessInputSpec.SuffixModifier                      = NULL;

    PatternLengthInputSpec.Type                                = Input_Double;
    PatternLengthInputSpec.Min                                 = 0.0001;
    PatternLengthInputSpec.Max                                 = 1.0;
    PatternLengthInputSpec.InterfaceAdjust.ScaleFact           = 100.0;
    PatternLengthInputSpec.SuffixModifier                      = NULL;

    TextBoxMarginInputSpec.Type                                = Input_Double;
    TextBoxMarginInputSpec.Min                                 = 0.0;
    TextBoxMarginInputSpec.Max                                 = 20.0;
    TextBoxMarginInputSpec.InterfaceAdjust.ScaleFact           = 100.0;
    TextBoxMarginInputSpec.SuffixModifier                      = NULL;

    TextLineSpacingInputSpec.Type                              = Input_Double;
    TextLineSpacingInputSpec.Min                               = 0.0;
    TextLineSpacingInputSpec.Max                               = 5.0;
    TextLineSpacingInputSpec.InterfaceAdjust.ScaleFact         = 1.0;
    TextLineSpacingInputSpec.SuffixModifier                    = NULL;


    ArrowheadSizeInputSpec.Type                                = Input_Double;
    ArrowheadSizeInputSpec.Min                                 = 0.0;
    ArrowheadSizeInputSpec.Max                                 = 0.5;
    ArrowheadSizeInputSpec.InterfaceAdjust.ScaleFact           = 100.0;
    ArrowheadSizeInputSpec.SuffixModifier                      = NULL;

    TextAngleInputSpec.Type                                    = Input_Double;
    TextAngleInputSpec.Min                                     = -PI - ANGLEEPSILON;
    TextAngleInputSpec.Max                                     =  PI + ANGLEEPSILON;
    TextAngleInputSpec.InterfaceAdjust.ScaleFact               = DEGPERRADIANS;
    TextAngleInputSpec.SuffixModifier                          = NULL;

    ArrowheadAngleInputSpec.Type                               = Input_Double;
    ArrowheadAngleInputSpec.Min                                = 1.0 / DEGPERRADIANS - ANGLEEPSILON;
    ArrowheadAngleInputSpec.Max                                = PIOVER2 + ANGLEEPSILON;
    ArrowheadAngleInputSpec.InterfaceAdjust.ScaleFact          = DEGPERRADIANS;
    ArrowheadAngleInputSpec.SuffixModifier                     = NULL;
}
#endif




void LocalReadBlock(FileStream_s   *FileStream,
                    double         *CurVPtr,
                    FieldDataType_e FieldDataTypeInFile,
                    HgIndex_t       NumValues,
                    Boolean_t      *IsOk)
{
    REQUIRE(VALID_REF(IsOk) && VALID_BOOLEAN(*IsOk));
    REQUIRE(!(*IsOk) || VALID_REF(FileStream));
    REQUIRE(!(*IsOk) || VALID_FIELD_DATA_TYPE(FieldDataTypeInFile));

    if (*IsOk)
    {
        Boolean_t DoRead = (CurVPtr != NULL);
        Boolean_t ReadByBlock = (FieldDataType_Double == FieldDataTypeInFile) || !DoRead;
        if (ReadByBlock)
        {
            ReadPureBlock(FileStream,
                          DoRead,
                          (void *)CurVPtr,
                          FieldDataTypeInFile,
                          0,
                          NumValues,
                          IsOk);
        }
        else
        {
            LgIndex_t N;
            for (N = 0; *IsOk && (N < NumValues); N++)
            {
                double D = GetNextValue(FileStream, FieldDataTypeInFile, -LARGEDOUBLE, LARGEDOUBLE, IsOk);
                if (DoRead)
                    CurVPtr[N] = D;
            }
        }
    }
}




/*
 *
 * NOTE: ReadTec only allocates space for NodeMap and VDataBase
 *       if RawDataSpaceAllocated == FALSE and GetHeaderInfoOnly
 *       is FALSE.
 *
 *       Also note that all data read in by ReadTec is currently
 *       limited to be in double precision.
 *
 */


Boolean_t STDCALL ReadTec(Boolean_t       GetHeaderInfoOnly,
                          char           *FName,
                          short          *IVersion,
                          char          **DataSetTitle,
                          EntIndex_t     *NumZones,
                          EntIndex_t     *NumVars,
                          StringList_pa  *VarNames,
                          StringList_pa  *ZoneNames,
                          LgIndex_t     **NumPtsI,
                          LgIndex_t     **NumPtsJ,
                          LgIndex_t     **NumPtsK,
                          ZoneType_e    **ZoneType,
                          StringList_pa  *UserRec,
                          AuxData_pa     *DatasetAuxData,
                          Boolean_t       RawDataSpaceAllocated,
                          NodeMap_t    ***NodeMap,
                          double       ***VDataBase)
{
    Boolean_t     InputIsOk         = FALSE;
    ArrayList_pa  ZoneSpecList      = NULL;
    LgIndex_t    *FNNumBndryConns   = NULL; /* [NumZones] */
    FileStream_s *ReadTecFileStream = NULL;
    Set_pa       *IsVarCellCentered = NULL; /* [NumZones] */

    REQUIRE(VALID_BOOLEAN(GetHeaderInfoOnly));
    REQUIRE(VALID_NON_ZERO_LEN_STR(FName));
    REQUIRE(VALID_REF(IVersion));
    REQUIRE(VALID_REF(DataSetTitle) || DataSetTitle == NULL);
    REQUIRE(VALID_REF(NumZones));
    REQUIRE(VALID_REF(NumVars));
    REQUIRE(VarNames  == NULL || VALID_REF(VarNames));
    REQUIRE(ZoneNames == NULL || VALID_REF(ZoneNames));
    REQUIRE(NumPtsI   == NULL || VALID_REF(NumPtsI));
    REQUIRE(NumPtsJ   == NULL || VALID_REF(NumPtsJ));
    REQUIRE(NumPtsK   == NULL || VALID_REF(NumPtsK));
    REQUIRE(ZoneType  == NULL || VALID_REF(ZoneType));
    REQUIRE(UserRec   == NULL || VALID_REF(UserRec));
    REQUIRE(VALID_BOOLEAN(RawDataSpaceAllocated));
    REQUIRE(IMPLICATION(!GetHeaderInfoOnly && RawDataSpaceAllocated,
                        VALID_REF(NodeMap) && VALID_REF(VDataBase)));

#if defined MAKEARCHIVE
    InitInputSpecs();
#endif

    InputIsOk  = OpenBinaryFileAndCheckMagicNumber(&ReadTecFileStream,
                                                   FName,
                                                   0,
                                                   IVersion);

    if (InputIsOk)
        InputIsOk = ReadDataFileHeader(ReadTecFileStream,
                                       *IVersion,
                                       FALSE,
                                       NumZones,
                                       NumVars,
                                       (SmInteger_t *)NULL,
                                       DataSetTitle,
                                       (Text_s **)NULL,
                                       (Geom_s **)NULL,
                                       (StringList_pa  **)NULL,
                                       UserRec,
                                       DatasetAuxData,
                                       &IsVarCellCentered,
                                       (Boolean_t *)NULL,
                                       (Boolean_t *)NULL,
                                       &ZoneSpecList,
                                       VarNames,
                                       (ArrayList_pa *)NULL,
                                       (Set_pa *)NULL,
                                       &FNNumBndryConns,
                                       (DataFileType_e *)NULL);



    if (InputIsOk)
    {
        if (*NumZones == 0)
            *NumVars = 0;
        else if (*IVersion > 112)
        {
            /*
             * This may not be true but we put it hear to remind us to make
             * updates to this code when we change the version number.
             */
            ErrMsg(translate("ReadTec does not yet support version %d "
                             "Tecplot binary data files."), *IVersion);
            InputIsOk = FALSE;
        }
        else if (!GetHeaderInfoOnly)
        {
            EntIndex_t Z;
            for (Z = 0; Z < *NumZones && InputIsOk; Z++)
            {
                InputIsOk = (MemberCount(IsVarCellCentered[Z]) == 0);
                if (!InputIsOk)
                    ErrMsg(translate("Cell centered data not supported by ReadTec."));
            }
        }
    }

    if (IsVarCellCentered != NULL)
    {
        EntIndex_t Z;
        for (Z = 0; Z < *NumZones; Z++)
            DeallocSet(&IsVarCellCentered[Z]);
        FREE_ARRAY(IsVarCellCentered, "Array of IsVarCellCentered sets");
    }

    if (InputIsOk)
    {
        EntIndex_t Z;
        /*
         *  Allocate space for the zone info pieces.
         */
        if (ZoneNames)
            *ZoneNames = StringListAlloc();
        if (NumPtsI)
            *NumPtsI  = ALLOC_ARRAY(*NumZones, LgIndex_t, "numptsi");
        if (NumPtsJ)
            *NumPtsJ  = ALLOC_ARRAY(*NumZones, LgIndex_t, "numptsj");
        if (NumPtsK)
            *NumPtsK  = ALLOC_ARRAY(*NumZones, LgIndex_t, "numptsk");
        if (ZoneType)
            *ZoneType = ALLOC_ARRAY(*NumZones, ZoneType_e, "zonetype");
        for (Z = 0; Z < *NumZones; Z++)
        {
            ZoneSpec_s *ZoneSpec = GetZoneSpec(ZoneSpecList, Z);
            if (ZoneSpec != NULL)
            {
                if (ZoneNames && *ZoneNames)
                    StringListAppendString(*ZoneNames, ZoneSpec->Name);

                if (NumPtsI && *NumPtsI)
                    (*NumPtsI)[Z] = ZoneSpec->NumPtsI;

                if (NumPtsJ && *NumPtsJ)
                    (*NumPtsJ)[Z] = ZoneSpec->NumPtsJ;

                if (NumPtsK && *NumPtsK)
                    (*NumPtsK)[Z] = ZoneSpec->NumPtsK;

                if (ZoneType && *ZoneType)
                    (*ZoneType)[Z] = ZoneSpec->Type;
            }
            else
            {
                if (ZoneNames && *ZoneNames)
                    StringListAppendString(*ZoneNames, NULL);

                if (NumPtsI && *NumPtsI)
                    (*NumPtsI)[Z] = 0;

                if (NumPtsJ && *NumPtsJ)
                    (*NumPtsJ)[Z] = 0;

                if (NumPtsK && *NumPtsK)
                    (*NumPtsK)[Z] = 0;

                if (ZoneType && *ZoneType)
                    (*ZoneType)[Z] = ZoneType_Invalid;
            }
        }
    }
    if (!GetHeaderInfoOnly && InputIsOk && (*NumZones > 0))
    {
        EntIndex_t      *VarSharesFromZone          = NULL; /* [NumVars] */
        Boolean_t       *IsVarPassive               = NULL; /* [NumVars] */
        EntIndex_t      *ConnectivitySharesFromZone = NULL; /* [NumZones] */
        FieldDataType_e *VarType = NULL;
        int              CurZone;
        int              CurVar;
        LgIndex_t        NumIPts = 0;
        LgIndex_t        NumJPts = 0;
        LgIndex_t        NumKPts = 0;
        LgIndex_t        TotalNumPts;
        LgIndex_t        I, J;

        if ((*NumZones > 0) && !RawDataSpaceAllocated)
        {
            *VDataBase   = ALLOC_ARRAY(*NumZones * (*NumVars), double *, "vdatabase array");
            if (*VDataBase == NULL)
            {
                ErrMsg(translate("Cannot allocate space for field data"));
                InputIsOk = FALSE;
            }
            else
            {
                int I;
                for (I = 0; I < *NumZones*(*NumVars); I++)
                    (*VDataBase)[I] = NULL;
            }

            if (InputIsOk)
            {
                *NodeMap = ALLOC_ARRAY(*NumZones, NodeMap_t *, "nodemap array");
                if (*NodeMap == NULL)
                {
                    ErrMsg(translate("Cannot allocate space for nodemap"));
                    InputIsOk = FALSE;
                }
                else
                {
                    int I;
                    for (I = 0; I < *NumZones; I++)
                        (*NodeMap)[I] = NULL;
                }
            }
        }

        if (InputIsOk)
        {
            VarType           = ALLOC_ARRAY(*NumVars + 1, FieldDataType_e, "Var Type");
            VarSharesFromZone = ALLOC_ARRAY(*NumVars + 1, EntIndex_t, "VarSharesFromZone");
            IsVarPassive      = ALLOC_ARRAY(*NumVars + 1, Boolean_t, "IsVarPassive");

            ConnectivitySharesFromZone = ALLOC_ARRAY(*NumZones, EntIndex_t, "ConnectivitySharesFromZone");
            InputIsOk = (VarType                    != NULL &&
                         VarSharesFromZone          != NULL &&
                         IsVarPassive               != NULL &&
                         ConnectivitySharesFromZone != NULL);
        }

        /* for each zone */
        for (CurZone = 0; CurZone < *NumZones && InputIsOk; CurZone++)
        {
            double X1 = GetNextValue(ReadTecFileStream, FieldDataType_Float, 0.0, 1000.0, &InputIsOk);
            if (InputIsOk && (X1 == ZoneMarker))
            {
                ZoneSpec_s *CurZoneSpec   = GetZoneSpec(ZoneSpecList, CurZone);
                Boolean_t   ZoneIsFinite  = (CurZoneSpec->Type != ZoneType_Ordered);
                Boolean_t   ZoneIsFEPoly  = (CurZoneSpec->Type == ZoneType_FEPolygon ||
                                             CurZoneSpec->Type == ZoneType_FEPolyhedron);
                Boolean_t   InBlockFormat = CurZoneSpec->ZoneLoadInfo.IsInBlockFormat;
                for (J = 0; J < *NumVars; J++)
                {
                    VarSharesFromZone[J] = -1; /* eumulate DupVar: no DupVar */
                    VarType[J]           = FieldDataType_Float;
                    IsVarPassive[J]      = FALSE;
                }

                /* dupvars */
                if (*IVersion > 45 && *IVersion < 101 && InputIsOk)
                {
                    EntIndex_t NumDupVars, ZZ;

                    NumDupVars = (EntIndex_t)GetIoFileInt(ReadTecFileStream, *IVersion, 0, (LgIndex_t) * NumVars, &InputIsOk);
                    for (J = 0; J < NumDupVars; J++)
                    {
                        ZZ = (EntIndex_t)GetIoFileInt(ReadTecFileStream, *IVersion, 0, *NumVars, &InputIsOk) - 1;
                        VarSharesFromZone[ZZ] = CurZone - 1; /* emulate DupVar: share from previous zone */
                    }
                    /* Can't duplicate from the first zone */
                    if ((NumDupVars > 0) && (CurZone == 0))
                    {
                        ErrMsg(translate("Cannot duplicate variables from the first zone since there are "
                                         "no previous zones to duplicate from."));
                        InputIsOk = FALSE;
                    }
                }

                /* get the data type for each variable */
                if (*IVersion >= 70 && InputIsOk)
                {
                    for (J = 0; J < *NumVars; J++)
                    {
                        VarType[J] = (FieldDataType_e)GetIoFileInt(ReadTecFileStream, *IVersion,
                                                                   0,
                                                                   (LgIndex_t)FieldDataType_Bit,
                                                                   &InputIsOk);
                        if (!InputIsOk)
                        {
                            ErrMsg(translate("Invalid data type - binary input file corrupted"));
                            InputIsOk = FALSE;
                        }
                    }
                }

                if (InputIsOk)
                {
                    NumIPts = CurZoneSpec->NumPtsI;
                    NumJPts = CurZoneSpec->NumPtsJ;
                    NumKPts = CurZoneSpec->NumPtsK;
                }

                if (ZoneIsFinite)
                    TotalNumPts = NumIPts;
                else
                    TotalNumPts = (NumIPts * NumJPts * NumKPts);

                for (CurVar = 0; CurVar < *NumVars && InputIsOk; CurVar++)
                {
                    if (!RawDataSpaceAllocated && TotalNumPts >= 1)
                    {
                        /*
                         * The calling program did not allocate space for the
                         * data so do it here.
                         */
                        (*VDataBase)[CurVar+CurZone*(*NumVars)] =
                            ALLOC_ARRAY(TotalNumPts, double, "raw data");
                    }
                }

                if (*IVersion >= 105 && InputIsOk)
                {
                    /* passive variables */
                    if ((Boolean_t)GetIoFileInt(ReadTecFileStream, *IVersion, 0, 1, &InputIsOk) && InputIsOk)
                    {
                        for (CurVar = 0; CurVar < *NumVars && InputIsOk; CurVar++)
                        {
                            IsVarPassive[CurVar] = (Boolean_t)GetIoFileInt(ReadTecFileStream,
                                                                           *IVersion,
                                                                           0, 1, &InputIsOk);
                        }
                    }
                }

                if (*IVersion >= 101 && InputIsOk)
                {
                    /* variable sharing: equivalent to DupVar for ReadTec */
                    if ((Boolean_t)GetIoFileInt(ReadTecFileStream, *IVersion, 0, 1, &InputIsOk) && InputIsOk)
                    {
                        for (CurVar = 0; CurVar < *NumVars && InputIsOk; CurVar++)
                        {
                            EntIndex_t SharedZone = GetIoFileInt(ReadTecFileStream, *IVersion,
                                                                 -1, MaxNumZonesOrVars - 1,
                                                                 &InputIsOk);
                            if (SharedZone != -1 && InputIsOk)
                                VarSharesFromZone[CurVar] = SharedZone;
                        }
                    }

                    /* face neighbor or FE node connectivity sharing */
                    if (InputIsOk)
                    {
                        EntIndex_t SharedZone = GetIoFileInt(ReadTecFileStream, *IVersion,
                                                             -1, MaxNumZonesOrVars - 1,
                                                             &InputIsOk);
                        if (InputIsOk)
                            ConnectivitySharesFromZone[CurZone] = SharedZone;
                    }
                }

                /*
                 * Non-shared variable min/max (but not for Zombie zones).
                 */
                if (*IVersion >= 103 && InputIsOk)
                {
                    for (CurVar = 0; CurVar < *NumVars && InputIsOk; CurVar++)
                    {
                        if (VarSharesFromZone[CurVar] == -1 && !IsVarPassive[CurVar])
                        {
                            /*
                             * Currently ReadTec doesn't do anything with the
                             * min/max values.
                             */
                            GetNextValue(ReadTecFileStream, FieldDataType_Double,
                                         -LARGEDOUBLE, LARGEDOUBLE,
                                         &InputIsOk);
                            GetNextValue(ReadTecFileStream, FieldDataType_Double,
                                         -LARGEDOUBLE, LARGEDOUBLE,
                                         &InputIsOk);
                        }
                    }
                }

                if (InBlockFormat)
                {
                    CurVar = -1;
                    while (InputIsOk && ((CurVar + 1) < *NumVars))
                    {
                        CurVar++;
                        if ((CurVar < *NumVars) && (TotalNumPts > 0))
                        {
                            double *CurVPtr  = (*VDataBase)[CurVar+CurZone*(*NumVars)];
                            J = 0;
                            if (VarSharesFromZone[CurVar] != -1)
                            {
                                LgIndex_t M;
                                EntIndex_t SourceZone = VarSharesFromZone[CurVar];
                                double *SourceVPtr = (*VDataBase)[CurVar+SourceZone*(*NumVars)];
                                for (M = 0; M < TotalNumPts; M++)
                                    CurVPtr[M] = SourceVPtr[M];
                            }
                            else if (!IsVarPassive[CurVar])
                            {
                                LocalReadBlock(ReadTecFileStream,
                                               CurVPtr,
                                               VarType[CurVar],
                                               TotalNumPts,
                                               &InputIsOk);
                            }
                        }
                    }
                    if (!InputIsOk)
                        ErrMsg(translate("Invalid raw data section of binary file"));
                }
                else if (TotalNumPts > 0)
                {
                    /*
                     * Zone is not empty and is in POINT format
                     */
                    J = -1;
                    if (InputIsOk)
                    {
                        LgIndex_t N;
                        N = 0;
                        while (InputIsOk && (N < TotalNumPts))
                        {
                            EntIndex_t CurVar;
                            for (CurVar = 0; InputIsOk && (CurVar < *NumVars); CurVar++)
                            {
                                double *CurVPtr  = (*VDataBase)[CurVar+CurZone*(*NumVars)];
                                if (VarSharesFromZone[CurVar] != -1)
                                {
                                    EntIndex_t SourceZone = VarSharesFromZone[CurVar];
                                    double *SourceVPtr = (*VDataBase)[CurVar+SourceZone*(*NumVars)];
                                    CurVPtr[N] = SourceVPtr[N];
                                }
                                else if (!IsVarPassive[CurVar])
                                {
                                    double D = GetNextValue(ReadTecFileStream,
                                                            VarType[CurVar],
                                                            -LARGEDOUBLE,
                                                            LARGEDOUBLE,
                                                            &InputIsOk);

                                    if (InputIsOk && CurVPtr)
                                        CurVPtr[N] = D;
                                }
                            }

                            if (!InputIsOk)
                                ErrMsg(translate("Binary datafile corrupted!"));
                            N++;
                        }
                    }
                }

                if (InputIsOk && *IVersion < 101)
                {
                    if (ZoneIsFinite)
                    {
                        /*
                         * Pre-version 101 had FE connectivity sharing,
                         * FECONNECT, information here.
                         */
                        Boolean_t DupConnectivity;
                        if (*IVersion > 61)
                            DupConnectivity = GetIoFileInt(ReadTecFileStream, *IVersion, 0, 1, &InputIsOk);
                        else
                            DupConnectivity = FALSE;

                        if (DupConnectivity)
                            ConnectivitySharesFromZone[CurZone] = CurZone - 1; /* previous zone */
                        else
                            ConnectivitySharesFromZone[CurZone] = -1;
                    }
                    else
                        ConnectivitySharesFromZone[CurZone] = -1;
                }

                if (InputIsOk && ZoneIsFinite && !ZoneIsFEPoly)
                {
                    Boolean_t   SkipNodemap;
                    NodeMap_t  *NM = NULL;
                    NodeMap_t  *ONM = NULL;
                    /*
                     *  Allocate the nodemap ptr if necessary Note that if
                     *  RawDataSpaceAllocated is TRUE then (*NodeMap)[CurZone]
                     *  can either contain a valid address (read the connectivity
                     *  list) or be NULL (skip the list).
                     */
                    if (!RawDataSpaceAllocated && NumKPts*NumJPts >= 1)
                    {
                        (*NodeMap)[CurZone] = ALLOC_ARRAY(NumKPts * NumJPts, NodeMap_t, "node map");
                        if ((*NodeMap)[CurZone] == NULL)
                            ErrMsg(translate("Cannot allocate space for connectivity list",
                                             "See the Tecplot User's Manual for a definition of 'connectivity list'"));
                    }

                    if (InputIsOk)
                        NM = (*NodeMap)[CurZone];

                    SkipNodemap = (NM == NULL);

                    if (InputIsOk && ConnectivitySharesFromZone[CurZone] != -1)
                    {
                        EntIndex_t SourceZone = ConnectivitySharesFromZone[CurZone];
                        if (SourceZone >= CurZone)
                        {
                            ErrMsg(translate("Zone %d is attempting to share connectivity "
                                             "with a zone that has not yet been loaded."),
                                   CurZone + 1);
                            InputIsOk = FALSE;
                        }
                        else
                        {
                            ONM = (*NodeMap)[SourceZone];
                            if (ONM == NULL)
                            {
                                ErrMsg(translate("Zone %d is attempting to share connectivity "
                                                 "with a zone that is not finite element."),
                                       CurZone + 1);
                                InputIsOk = FALSE;
                            }
                        }
                    }

                    if (InputIsOk)
                    {
                        /* load the FE node connectivity */
                        for (J = 0; J < NumJPts; J++)
                            for (I = 0; I < NumKPts; I++)
                            {
                                LgIndex_t M;
                                LgIndex_t L = J * NumKPts + I;
                                if (ConnectivitySharesFromZone[CurZone] != -1)
                                    M = ONM[L];
                                else
                                    M = GetNextI(ReadTecFileStream, &InputIsOk) - 1;
                                if (!SkipNodemap)
                                    NM[L] = M;
                            }
                    }
                }

                /* skip over the face neighbor connectivity */
                if (*IVersion >= 101 && InputIsOk)
                {
                    EntIndex_t SharedZone = ConnectivitySharesFromZone[CurZone];
                    if (SharedZone == -1 && FNNumBndryConns[CurZone] != 0)
                    {
                        LgIndex_t Connection = 0;
                        while (Connection < FNNumBndryConns[CurZone] && InputIsOk)
                        {
                            /*
                             * Face neighbor connection have the following format for both
                             * ASCII and binary:
                             *
                             *   FaceNeighborMode_LocalOneToOne     3         cz,fz,cz
                             *   FaceNeighborMode_LocalOneToMany    nz+4      cz,fz,oz,nz,cz1,cz2,...,czn
                             *   FaceNeighborMode_GlobalOneToOne    4         cz,fz,ZZ,CZ
                             *   FaceNeighborMode_GlobalOneToMany   2*nz+4    cz,fz,oz,nz,ZZ1,CZ1,ZZ2,CZ2,...,ZZn,CZn
                             *
                             *   Where:
                             *       cz = cell in current zone
                             *       fz = face of cell in current zone
                             *       oz = face obsuration flag (only applies to one-to-many):
                             *              0 = face partially obscured
                             *              1 = face entirely obscured
                             *       nz = number of cell or zone/cell associations (only applies to one-to-many)
                             *       ZZ = remote Zone
                             *       CZ = cell in remote zone
                             */
                            (void)GetNextI(ReadTecFileStream, &InputIsOk); /* read cz */
                            if (!InputIsOk)
                                ErrMsg(translate("Unexpected end-of-file while reading face neighbor data."));

                            (void)GetNextI(ReadTecFileStream, &InputIsOk); /* read fz */

                            if (InputIsOk)
                            {
                                /*
                                 * read FaceNeighborMode_LocalOneToOne:   cz ||
                                 *      FaceNeighborMode_LocalOneToMany:  oz ||
                                 *      FaceNeighborMode_GlobalOneToOne:  ZZ ||
                                 *      FaceNeighborMode_GlobalOneToMany: oz
                                 */
                                if (CurZoneSpec->FNMode == FaceNeighborMode_LocalOneToOne)
                                    (void)GetNextI(ReadTecFileStream, &InputIsOk);
                                else if (CurZoneSpec->FNMode == FaceNeighborMode_LocalOneToMany)
                                    (void)GetNextI(ReadTecFileStream, &InputIsOk);
                                else if (CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToOne)
                                    (void)GetNextI(ReadTecFileStream, &InputIsOk);
                                else if (CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToMany)
                                    (void)GetNextI(ReadTecFileStream, &InputIsOk);
                                else
                                    CHECK(FALSE);

                                if (CurZoneSpec->FNMode != FaceNeighborMode_LocalOneToOne && InputIsOk)
                                {
                                    LgIndex_t NumAssociations = 0;
                                    /*
                                     * read FaceNeighborMode_LocalOneToMany:  nz ||
                                     *      FaceNeighborMode_GlobalOneToOne:  CZ ||
                                     *      FaceNeighborMode_GlobalOneToMany: nz
                                     */
                                    if (CurZoneSpec->FNMode == FaceNeighborMode_LocalOneToMany)
                                        NumAssociations = GetNextI(ReadTecFileStream, &InputIsOk);
                                    else if (CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToOne)
                                        (void)GetNextI(ReadTecFileStream, &InputIsOk);
                                    else if (CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToMany)
                                        NumAssociations = GetNextI(ReadTecFileStream, &InputIsOk);
                                    else
                                        CHECK(FALSE);

                                    if (CurZoneSpec->FNMode != FaceNeighborMode_GlobalOneToOne && InputIsOk)
                                    {
                                        LgIndex_t Assoc;
                                        if (CurZoneSpec->FNMode == FaceNeighborMode_LocalOneToMany)
                                            for (Assoc = 0; Assoc < NumAssociations && InputIsOk; Assoc++)
                                                (void)GetNextI(ReadTecFileStream, &InputIsOk); /* read czn */
                                        else if (CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToMany)
                                            for (Assoc = 0; Assoc < NumAssociations && InputIsOk; Assoc++)
                                            {
                                                (void)GetNextI(ReadTecFileStream, &InputIsOk); /* read ZZn */
                                                (void)GetNextI(ReadTecFileStream, &InputIsOk); /* read CZn */
                                            }
                                        else
                                            CHECK(FALSE);

                                        if (InputIsOk)
                                            Connection += NumAssociations;
                                    }
                                    else if (InputIsOk) /* CurZoneSpec->FNMode == FaceNeighborMode_GlobalOneToOne */
                                        Connection += 1;
                                }
                                else if (InputIsOk) /* CurZoneSpec->FNMode == FaceNeighborMode_LocalOneToOne */
                                    Connection += 1;

                                if (!InputIsOk)
                                    ErrMsg(translate("Corrupt input file: invalid face neighbors."));
                            }
                        }
                    }
                }/* face neighbor connectivity */
                /* skip over face map section */
                if (ZoneIsFEPoly                              &&
                    *IVersion >= 110                          &&
                    ConnectivitySharesFromZone[CurZone] != -1 &&
                    InputIsOk)
                {
                    if (!InBlockFormat)
                    {
                        ErrMsg(translate("Poly zones must be in block format"));
                        InputIsOk = FALSE;
                    }
                    if (InputIsOk)
                    {
                        HgIndex_t NumFaces = CurZoneSpec->NumPtsK;
                        if (*IVersion == 110) // ...version 111 moved these to the zone header
                        {
                            CurZoneSpec->NumFaceNodes  = GetNextI(ReadTecFileStream, &InputIsOk);
                            CurZoneSpec->NumFaceBndryFaces = GetNextI(ReadTecFileStream, &InputIsOk);
                            CurZoneSpec->NumFaceBndryItems = GetNextI(ReadTecFileStream, &InputIsOk);
                        }
                        HgIndex_t TotalNumFaceNodes  = CurZoneSpec->NumFaceNodes;
                        HgIndex_t TotalNumBndryFaces = CurZoneSpec->NumFaceBndryFaces;
                        HgIndex_t TotalNumBndryItems = CurZoneSpec->NumFaceBndryItems;
                        if (CurZoneSpec->Type == ZoneType_FEPolyhedron)
                            ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, NumFaces + 1, &InputIsOk);
                        if (InputIsOk)
                            ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, TotalNumFaceNodes, &InputIsOk);
                        if (InputIsOk)
                            ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, NumFaces, &InputIsOk);
                        if (InputIsOk)
                            ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, NumFaces, &InputIsOk);
                        if (TotalNumBndryFaces > 0)
                        {
                            if (InputIsOk)
                                ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, TotalNumBndryFaces + 1, &InputIsOk);
                            if (InputIsOk)
                                ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, TotalNumBndryItems, &InputIsOk);
                            if (InputIsOk)
                            {
                                if (*IVersion >= 112)
                                    ReadInt32Block(ReadTecFileStream, FALSE, NULL, 0, TotalNumBndryItems, &InputIsOk);
                                else
                                    ReadInt16Block(ReadTecFileStream, FALSE, NULL, 0, TotalNumBndryItems, &InputIsOk);
                            }
                        }
                    }
                }/* face map section */
            }
            else
            {
                ErrMsg(translate("Corrupt input file"));
                InputIsOk = FALSE;
            }
        }

        if (VarSharesFromZone)
            FREE_ARRAY(VarSharesFromZone, "VarSharesFromZone");
        if (IsVarPassive)
            FREE_ARRAY(IsVarPassive, "IsVarPassive");
        if (ConnectivitySharesFromZone)
            FREE_ARRAY(ConnectivitySharesFromZone, "ConnectivitySharesFromZone");
        if (VarType)
            FREE_ARRAY(VarType, "VarType");

        if (!InputIsOk && !RawDataSpaceAllocated)
        {
            int I;
            if (*VDataBase)
            {
                for (I = 0; I < *NumZones*(*NumVars); I++)
                {
                    if ((*VDataBase)[I])
                        FREE_ARRAY((*VDataBase)[I], "vdatabase array");
                }
                FREE_ARRAY(*VDataBase, "vdatabase pointer array");
            }


            if (*NodeMap)
            {
                for (I = 0; I < *NumZones; I++)
                {
                    if ((*NodeMap)[I])
                        FREE_ARRAY((*NodeMap)[I], "connectivity list");
                }
                FREE_ARRAY(*NodeMap, "connectivity pointer array");
            }
        }
    } /*Reading Raw Data*/

    if (FNNumBndryConns != NULL)
        FREE_ARRAY(FNNumBndryConns, "FNNumBndryConns");
    if (ZoneSpecList)
        ArrayListDealloc(&ZoneSpecList, ZoneSpecItemDestructor, 0);

    if (ReadTecFileStream)
    {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
        TP_FCLOSE(ReadTecFileStream->File);
        free(ReadTecFileStream);
#endif
    }
    return (InputIsOk);
}


void * STDCALL TecAlloc(size_t size)
{
    return (void *)ALLOC_ARRAY(size, char, "TecAlloc");
}

void STDCALL TecFree(void *ptr)
{
    /* Hack to remove delete warning... */
    char *Tmp = (char *)ptr;
    FREE_ARRAY(Tmp, "TecAlloc");
}


