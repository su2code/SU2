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

#define DATAIOMODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "TranslatedString.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
#ifdef MSWIN
/* Disable warning about conversion from long to short.
   Even if we have disabled it in stdafx.h,
   we still need to disable here also for
   tecio, which includes this cpp file. */

#pragma warning (disable : 4244)
#endif
/*
 * Temp text and geom buffers.
 */
static Geom_s TempGeom;
static Text_s TempText;
#endif

#include "DATASET0.h"
#include "SET.h"
#include "FILESTREAM.h"
#include "DATAIO.h"
#include "DATAIO4.h"
#include "STRUTIL.h"
#include "AUXDATA.h"
#include "ARRLIST.h"
#include "STRLIST.h"
#include "ALLOC.h"
#include "DATASET.h"
#include "SYSTEM.h"
#include "Q_MSG.h"

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
using namespace tecplot::strutil;

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#endif /* TECPLOTKERNEL */


/*
 * Multi-Purpose datafile header reader.  This is designed so that
 * not all parts of the header is loaded and so that the task of loading
 * the header information is separated out from installing of that
 * information into a dataset.
 *
 */

Boolean_t ReadDataFileHeader(FileStream_s    *FileStream,
                             short            IVersion,
                             Boolean_t        ShowDataIOStatus,
                             EntIndex_t      *NumZones,
                             EntIndex_t      *NumVars,
                             SmInteger_t     *NumCustomLabelSets,
                             char           **DataSetTitle,
                             Text_s         **BaseText,
                             Geom_s         **BaseGeom,
                             StringList_pa  **CustomLabelBase,
                             StringList_pa   *UserRec,
                             AuxData_pa      *DataSetAuxData,
                             Set_pa         **IsVarCellCentered,  /* Create an Array dim by zones */
                             Boolean_t       *HasText,
                             Boolean_t       *HasGeoms,
                             ArrayList_pa    *ZoneSpecList,
                             StringList_pa   *VarNames,
                             ArrayList_pa    *VarAuxDataList, /*<AuxData_pa>[NumVars]*/
                             Set_pa          *IsRawFNAvailable, /* classic data only */
                             LgIndex_t      **FNNumBndryConns,  /* classic data only */
                             DataFileType_e  *FileType)
{
    Boolean_t    IsOk = TRUE;
    Boolean_t    SentError = FALSE;
    double       X1;
    int          Pass;
    FileOffset_t InitialFilePosition;

    REQUIRE(VALID_REF(FileStream) && VALID_REF(FileStream->File));
    REQUIRE(IVersion > 0);
    REQUIRE(VALID_BOOLEAN(ShowDataIOStatus));
    REQUIRE(VALID_REF(NumZones));
    REQUIRE(VALID_REF(NumVars));
    REQUIRE(VALID_REF(DataSetTitle)       || (DataSetTitle == NULL));
    REQUIRE(VALID_REF(BaseText)           || (BaseText == NULL));
    REQUIRE(VALID_REF(BaseGeom)           || (BaseGeom == NULL));
    REQUIRE(VALID_REF(HasText)            || (HasText == NULL));
    REQUIRE(VALID_REF(HasGeoms)           || (HasGeoms == NULL));
    REQUIRE(VALID_REF(ZoneSpecList)       || (ZoneSpecList == NULL));
    REQUIRE(VALID_REF(VarNames)           || (VarNames == NULL));
    REQUIRE(VALID_REF(NumCustomLabelSets) || (NumCustomLabelSets == NULL));
    REQUIRE(VALID_REF(UserRec)            || (UserRec == NULL));
    REQUIRE((VALID_REF(DataSetAuxData) &&
             (VALID_REF(*DataSetAuxData) || *DataSetAuxData == NULL)) ||
            DataSetAuxData == NULL);
    REQUIRE((VALID_REF(VarAuxDataList) &&
             (VALID_REF(*VarAuxDataList) || *VarAuxDataList == NULL)) ||
            VarAuxDataList == NULL);
    REQUIRE(VALID_REF(IsVarCellCentered)  || (IsVarCellCentered == NULL));
    REQUIRE((VALID_REF(CustomLabelBase) && VALID_REF(NumCustomLabelSets)) || (CustomLabelBase == NULL));
    REQUIRE(VALID_REF(IsRawFNAvailable) || (IsRawFNAvailable == NULL));
    REQUIRE(VALID_REF(FNNumBndryConns)  || (FNNumBndryConns == NULL));
    REQUIRE(VALID_REF(FileType)         || (FileType == NULL));

    if (DataSetTitle)
        *DataSetTitle = NULL;
    if (BaseText)
        *BaseText = NULL;
    if (BaseGeom)
        *BaseGeom = NULL;
    if (HasText)
        *HasText = FALSE;
    if (HasGeoms)
        *HasGeoms = FALSE;
    if (ZoneSpecList)
        *ZoneSpecList = NULL;
    if (VarNames)
        *VarNames = NULL;
    if (NumCustomLabelSets)
        *NumCustomLabelSets = 0;
    if (CustomLabelBase)
        *CustomLabelBase = NULL;
    if (DataSetAuxData != NULL)
    {
        /*
         * Note unlike most of the other output only parameters that we nullify,
         * DataSetAuxData is both an input and output parameter therefore we do
         * not nullify it. The CHECK is here for clarity.
         */
        CHECK(VALID_REF(*DataSetAuxData) || *DataSetAuxData == NULL);
    }
    if (VarAuxDataList != NULL)
    {
        /*
         * Note unlike most of the other output only parameters that we nullify,
         * VarAuxDataList is both an input and output parameter therefore we do
         * not nullify it. The CHECK is here for clarity.
         */
        CHECK(VALID_REF(*VarAuxDataList) || *VarAuxDataList == NULL);
    }
    if (UserRec)
        *UserRec = NULL;
    if (IsVarCellCentered)
        *IsVarCellCentered = NULL;

    if (IsRawFNAvailable)
        *IsRawFNAvailable = NULL;

    if (FNNumBndryConns)
        *FNNumBndryConns = NULL;

    if (FileType)
        *FileType = DataFileType_Full;

    /*
     * Pass 1 is used only to count up the number of zones and custom label sets,
     * Also determine if there are any preset zone colors.
     */

    InitialFilePosition = TP_FTELL(FileStream->File);

    for (Pass = 1; IsOk && (Pass <= 2); Pass++)
    {
        if (Pass == 2)
        {
            if (TP_FSEEK(FileStream->File, InitialFilePosition, SEEK_SET) != 0)
                IsOk = FALSE;

            if (IsOk && (*NumZones > 0 && ZoneSpecList != NULL && *ZoneSpecList == NULL))
            {
                *ZoneSpecList = ArrayListAlloc(*NumZones, ArrayListType_VoidPtr,
                                               ZoneOrVarListAdjustCapacityRequest, 0);
                IsOk = (*ZoneSpecList != NULL);
            }
            if (IsOk && (CustomLabelBase != NULL    &&
                         *CustomLabelBase == NULL &&
                         *NumCustomLabelSets > 0))
            {
                *CustomLabelBase = ALLOC_ARRAY(*NumCustomLabelSets, StringList_pa, "CustomLabel Sets");
                IsOk = (*CustomLabelBase != NULL);
                if (IsOk)
                {
                    SmInteger_t N;
                    for (N = 0; N < *NumCustomLabelSets; N++)
                        (*CustomLabelBase)[N] = NULL;
                }
            }
            if (IsOk && (UserRec != NULL && *UserRec == NULL))
            {
                *UserRec = StringListAlloc();
                IsOk = (Boolean_t)(*UserRec != NULL);
            }
            if (IsOk && (DataSetAuxData != NULL && *DataSetAuxData == NULL))
            {
                *DataSetAuxData = AuxDataAlloc();
                IsOk = (Boolean_t)(*DataSetAuxData != NULL);
            }
            if (IsOk && (VarAuxDataList != NULL && *VarAuxDataList == NULL) && *NumVars > 0)
            {
                *VarAuxDataList = ArrayListAlloc(0, ArrayListType_VoidPtr,
                                                 ZoneOrVarListAdjustCapacityRequest, 0);
                IsOk = (*VarAuxDataList != NULL &&
                        ArrayListSetVoidPtr(*VarAuxDataList, *NumVars - 1, NULL));
            }
            if (IsOk            &&
                (*NumZones > 0) &&
                (IsVarCellCentered != NULL) &&
                (*IsVarCellCentered == NULL))
            {
                /*
                 * First construct the array of sets...
                 */
                *IsVarCellCentered = ALLOC_ARRAY(*NumZones, Set_pa, "Array of IsVarCellCentered sets");
                if (*IsVarCellCentered)
                {
                    EntIndex_t Z;
                    for (Z = 0; IsOk && (Z < *NumZones); Z++)
                    {
                        /*
                         * Now allocate a set for each zone
                         */
                        (*IsVarCellCentered)[Z] = AllocSet(FALSE);
                        IsOk = (Boolean_t)((*IsVarCellCentered)[Z] != NULL);
                    }
                }
                else
                    IsOk = FALSE;
            }
            if (IsOk && *NumZones > 0 && IsRawFNAvailable != NULL)
            {
                *IsRawFNAvailable = AllocSet(FALSE);
                IsOk = (*IsRawFNAvailable != NULL);
            }
            if (IsOk && *NumZones > 0 && FNNumBndryConns != NULL)
            {
                *FNNumBndryConns = ALLOC_ARRAY(*NumZones, LgIndex_t, "Array of FNNumBndryConns");
                IsOk = (*FNNumBndryConns != NULL);
                if (IsOk)
                    for (LgIndex_t i = 0; i < *NumZones; i++)
                        (*FNNumBndryConns)[i] = 0;
            }
        }

        if (NumCustomLabelSets != NULL)
            *NumCustomLabelSets = 0;

        EntIndex_t TotalNumZones = *NumZones; /* ...only meaningful for pass 2 */

        *NumZones = 0;
        *NumVars  = 0;

        if (IsOk)
        {
            char *S = NULL;
            int  INumVars;

            IsOk = ReadInDataFileTypeTitleAndVarNames(FileStream,
                                                      IVersion,
                                                      ((Pass == 2) ? &S : (char **)NULL),
                                                      ((Pass == 2) ? FileType : (DataFileType_e *)NULL),
                                                      &INumVars,
                                                      ((Pass == 2) ? VarNames : (StringList_pa *)NULL));

            if (IsOk)
                *NumVars = (EntIndex_t)INumVars;

            if ((Pass == 2) && S && IsOk && DataSetTitle)
                *DataSetTitle = S;
            else if (S != NULL)
                FREE_ARRAY(S, "data set title");
        }

        if (IsOk)
        {
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no dialog feedback */
          //  LgIndex_t NumGeoms = 0;
          //  LgIndex_t NumTexts = 0;
#endif

            if (IsOk)
                X1 = GetNextValue(FileStream, FieldDataType_Float, 0.0, 1000.0, &IsOk);

            while (IsOk && (X1 != EndHeaderMarker))
            {
                if (X1 == ZoneMarker)
                {
                    ZoneSpec_s *ZoneSpec = ZoneSpecAlloc();
                    Boolean_t OkToLoad = (Pass == 2 &&
                                          IsVarCellCentered != NULL);
                    IsOk = (ZoneSpec != NULL);
                    if (IsOk)
                    {
                        Boolean_t LocalIsRawFNAvailable;
                        LgIndex_t LocalFNNumBndryConns;
                        IsOk = ReadInZoneHeader(FileStream, IVersion, ZoneSpec,
                                                OkToLoad ? (*IsVarCellCentered)[*NumZones] : NULL,
                                                *NumVars, &LocalIsRawFNAvailable,
                                                &LocalFNNumBndryConns);
                        if (IsOk && OkToLoad && IsRawFNAvailable != NULL)
                        {
                            if (LocalIsRawFNAvailable)
                                IsOk = AddToSet(*IsRawFNAvailable, *NumZones, FALSE);
                        }
                        if (IsOk && OkToLoad && FNNumBndryConns != NULL)
                            (*FNNumBndryConns)[*NumZones] = LocalFNNumBndryConns;
                    }

                    if (IsOk                 &&
                        ZoneSpecList != NULL &&
                        Pass == 2)
                    {
                        IsOk = (ZoneSpec->ParentZone == BAD_SET_VALUE ||
                                (ZoneSpec->ParentZone != *NumZones &&
                                 (0 <= ZoneSpec->ParentZone && ZoneSpec->ParentZone < TotalNumZones)));
                        if (IsOk)
                        {
                            ArrayListItem_u CurZoneSpecItem;
                            CurZoneSpecItem.VoidPtr = (void *)ZoneSpec;
                            ArrayListSetItem(*ZoneSpecList, *NumZones,
                                             CurZoneSpecItem,
                                             ZoneSpecItemDestructor, 0);
                        }
                        else
                        {
                            if (ZoneSpec->ParentZone == *NumZones)
                                ErrMsg(translate("Parent zone assignment for zone %d "
                                                 "may not be self referencing."),
                                       *NumZones + 1);
                            else
                                ErrMsg(translate("Parent zone assignment for zone %d "
                                                 "must be to an existing zone within the datafile."),
                                       *NumZones + 1);
                            ZoneSpecDealloc(&ZoneSpec);
                            SentError = TRUE;
                        }
                    }
                    else
                        ZoneSpecDealloc(&ZoneSpec);

                    if (IsOk)
                        (*NumZones)++;
                    if (*NumZones > MaxNumZonesOrVars)
                    {
                        ErrMsg(translate("Exceeding Tecplot's current zone limit of %d. "
                                         "Reduce the number of zones being loaded."), MaxNumZonesOrVars);
                        IsOk = FALSE;
                        SentError = TRUE;
                    }
                }
                else if (X1 == GeomMarker)
                {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
                    IsOk = ReadInGeometry(FileStream,
                                          IVersion,
                                          FALSE,
                                          &TempGeom,
                                          5000);
#endif
                    if (IsOk)
                    {
                        if (Pass == 1)
                        {
                            if (HasGeoms)
                                *HasGeoms = TRUE;
                        }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
                    }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#endif
                }
                else if (X1 == TextMarker)
                {
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
                    IsOk = ReadInText(FileStream,
                                      IVersion,
                                      FALSE,
                                      &TempText,
                                      200);
#endif
                    if (IsOk)
                    {
                        if (Pass == 1)
                        {
                            if (HasText)
                                *HasText = TRUE;
                        }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
                    }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#endif
                }
                else if (X1 == CustomLabelMarker)
                {
                    Boolean_t OkToLoad;

                    OkToLoad = (Pass == 2)                                &&
                               NumCustomLabelSets                         &&
                               (*NumCustomLabelSets < MaxCustomLabelSets) &&
                               CustomLabelBase;

                    IsOk = ReadInCustomLabels(FileStream,
                                              IVersion,
                                              OkToLoad,
                                              (OkToLoad ? &(*CustomLabelBase)[*NumCustomLabelSets] : NULL));
                    if (IsOk && NumCustomLabelSets)
                        (*NumCustomLabelSets)++;
                }
                else if (X1 == UserRecMarker)
                {
                    Boolean_t OkToLoad;
                    char     *CurUserRec = NULL;

                    OkToLoad = (Boolean_t)((Pass == 2) && UserRec);

                    IsOk = ReadInUserRec(FileStream,
                                         IVersion,
                                         500,
                                         OkToLoad ? &CurUserRec : (char **)NULL);
                    if (IsOk && OkToLoad)
                        IsOk = StringListAppendString(*UserRec, CurUserRec);
                    if (CurUserRec)
                        FREE_ARRAY(CurUserRec, "temp user rec");
                    CurUserRec = NULL;
                }
                else if (X1 == DataSetAuxMarker)
                {
                    Boolean_t OkToLoad;
                    CHECK(IVersion >= 101);
                    OkToLoad = (Pass == 2 &&
                                DataSetAuxData != NULL);
                    IsOk = ReadInAuxData(FileStream, IVersion,
                                         OkToLoad ? *DataSetAuxData : NULL);
                    if (!IsOk)
                    {
                        ErrMsg(translate("Invalid DATASETAUXDATA record in binary datafile"));
                        SentError = TRUE;
                    }
                }
                else if (X1 == VarAuxMarker)
                {
                    Boolean_t OkToLoad;
                    LgIndex_t VarNum;
                    CHECK(IVersion >= 102);
                    OkToLoad = (Pass == 2 &&
                                VarAuxDataList != NULL);
                    VarNum = GetIoFileInt(FileStream, IVersion, 0, *NumVars - 1, &IsOk);
                    if (IsOk)
                    {
                        AuxData_pa VarAuxData;
                        if (OkToLoad)
                        {
                            VarAuxData = (AuxData_pa)ArrayListGetVoidPtr(*VarAuxDataList, VarNum);
                            if (VarAuxData == NULL)
                            {
                                VarAuxData = AuxDataAlloc();
                                IsOk = (VarAuxData != NULL &&
                                        ArrayListSetVoidPtr(*VarAuxDataList, VarNum, VarAuxData));
                            }
                        }
                        else
                            VarAuxData = NULL;

                        IsOk = IsOk && ReadInAuxData(FileStream, IVersion, VarAuxData);
                        if (!IsOk)
                        {
                            ErrMsg(translate("Invalid VARAUXDATA record in binary datafile"));
                            SentError = TRUE;
                        }
                    }
                    else
                    {
                        ErrMsg(translate("Invalid VARAUXDATA variable number association"));
                        SentError = TRUE;
                    }
                }
                else
                    IsOk = FALSE;
                if (IsOk)
                    X1 = GetNextValue(FileStream, FieldDataType_Float, 0.0, 1000.0, &IsOk);
            }
        }
    }

    /*
     * Old plt files that did not contain data still contained variable name
     * definitions in the header.  This is no longer necessary and in fact can
     * cause confusion in the data read options dialog.  If the number of zones
     * in the datafile is zero set the number of variables to 0 and dealloc the
     * variable name list.
     */

    if (IsOk && (*NumZones == 0) && (*NumVars > 0))
    {
        *NumVars = 0;
        if (VarNames && *VarNames)
        {
            StringListDealloc(VarNames);
        }
    }


    if (!IsOk)
    {
        if (ZoneSpecList && *ZoneSpecList)
            ArrayListDealloc(ZoneSpecList, ZoneSpecItemDestructor, 0);
        if (DataSetTitle && *DataSetTitle)
        {
            FREE_ARRAY(*DataSetTitle, "DataSetTitle");
            *DataSetTitle = NULL;
        }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        if (VarNames && *VarNames)
        {
            StringListDealloc(VarNames);
        }
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        if (UserRec && *UserRec)
            StringListDealloc(UserRec);
    }

    /* If there was an error, get rid of the auxiliary data list. */
    if ((DataSetAuxData != NULL && *DataSetAuxData != NULL) &&
        (!IsOk))
        AuxDataDealloc(DataSetAuxData);

    if (!IsOk && !SentError)
        ErrMsg(translate("Invalid header in binary datafile"));

    /*
     * NOTE: Do not close the file.  Some calling functions will continue
     *       to read from this point on.
     */

    ENSURE((VarNames == NULL) || (*VarNames == NULL) || StringListValid(*VarNames));
    ENSURE(IMPLICATION(UserRec != NULL,
                       (*UserRec == NULL ||
                        StringListValid(*UserRec))));
    ENSURE(IMPLICATION(DataSetAuxData != NULL,
                       (*DataSetAuxData == NULL ||
                        VALID_REF(*DataSetAuxData))));
    return (IsOk);
}



#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


Boolean_t OpenBinaryFileAndCheckMagicNumber(FileStream_s **FileStream,
                                            char          *FName,
                                            FileOffset_t   StartOffset,
                                            short         *IVersion)
{
    Boolean_t Result = TRUE;
    REQUIRE(VALID_REF(FileStream));
    REQUIRE(*FileStream == NULL);
    REQUIRE(VALID_REF(FName));
    REQUIRE(StartOffset >= 0);
    REQUIRE(VALID_REF(IVersion));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
    FILE *File = TP_FOPEN(FName, "rb");
    if (File == NULL)
        Result = FALSE;
#endif
    if (Result)
    {
        *FileStream = FileStreamAlloc(File, TRUE);
        Result = (*FileStream != NULL);
    }
    Result = Result && (TP_FSEEK((*FileStream)->File, StartOffset, SEEK_SET) == 0);
    if (Result)
    {
        *IVersion = GetInputVersion(*FileStream);
        /*
         * After version 71 we started listing individual valid
         * versions. Before that time we just listed ranges. Also,
         * note that versions 72, 73, and 74 were invalid.
         */
        Result = (/* past valid plt file version ranges follow: */
                     (40 <= *IVersion && *IVersion <=  71) ||
                     (*IVersion ==  75)                    ||
                     (100 <= *IVersion && *IVersion <= TecplotBinaryFileVersion));

        /*
         * This check is put here to make sure that the above code gets visited
         * when the TecplotBinaryFileVersion number changes. When the version
         * changes the "past valid plt file version ranges" above and the number
         * compared to the TecplotBinaryFileVersion below may need to be
         * adjusted such as when we skip a consecutive number as we did between
         * version 71 and 75, and between 75 and 100.
         */
        CHECK(TecplotBinaryFileVersion == 112);
    }

    ENSURE(VALID_BOOLEAN(Result));
    return (Result);
}




#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
            #if !defined NO_ASSERTS
            #endif
                #if !defined NO_ASSERTS
                #endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-M 12/13/2005: ENGINE-P2 - no status feedback */
#endif
#if !defined ENGINE /* TODO(RMS)-H 12/12/2005: ENGINE: refactor to use just the Interrupted flag as-is */
#else
#endif
#if 0 /* we changed this behavior... not sure when */
#endif
#endif /* TECPLOTKERNEL */
