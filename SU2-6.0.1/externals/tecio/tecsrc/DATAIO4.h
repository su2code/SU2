#ifndef DATAIO4_H
#define DATAIO4_H
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#include <set>

#if defined EXTERN
#undef EXTERN
#endif
#if defined DATAIO4MODULE
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN double GetNextValue(FileStream_s    *FileStream,
                           FieldDataType_e  FieldDataType,
                           double           Min,
                           double           Max,
                           Boolean_t       *IsOk);
EXTERN LgIndex_t GetNextI(FileStream_s *FileStream,
                          Boolean_t    *IsOk);
EXTERN LgIndex_t GetIoFileInt(FileStream_s *FileStream,
                              short         Version,
                              LgIndex_t     Min,
                              LgIndex_t     Max,
                              Boolean_t    *IsOk);
EXTERN Boolean_t ReadInString(FileStream_s  *FileStream,
                              short          IVersion,
                              int            MaxCharacters,
                              char         **S,
                              Boolean_t      ProcessData);
EXTERN void ReadByteBlock(FileStream_s *FileStream,
                          Boolean_t     DoRead,
                          Byte_t       *Buffer,
                          HgIndex_t     StartIndex,
                          HgIndex_t     NumValues,
                          Boolean_t    *IsOk);
EXTERN void ReadInt16Block(FileStream_s *FileStream,
                           Boolean_t     DoRead,
                           Int16_t      *Buffer,
                           HgIndex_t     StartIndex,
                           HgIndex_t     NumValues,
                           Boolean_t    *IsOk);
EXTERN void ReadInt16BlockToInt32(FileStream_s *FileStream,
                                  Boolean_t     DoRead,
                                  Int32_t      *Buffer,
                                  HgIndex_t     StartIndex,
                                  HgIndex_t     NumValues,
                                  Boolean_t    *IsOk);
EXTERN void ReadInt32Block(FileStream_s *FileStream,
                           Boolean_t     DoRead,
                           Int32_t      *Buffer,
                           HgIndex_t     StartIndex,
                           HgIndex_t     NumValues,
                           Boolean_t    *IsOk);
EXTERN void ReadPureBlock(FileStream_s   *FileStream,
                          Boolean_t       DoRead,
                          void           *Buffer,
                          FieldDataType_e FieldDataType,
                          HgIndex_t       StartIndex,
                          HgIndex_t       NumValues,
                          Boolean_t      *IsOk);
EXTERN void ReadBlock(FileStream_s   *FileStream,
                      FieldData_pa    FieldData,
                      Boolean_t       DoRead,
                      FieldDataType_e FieldDataTypeInFile,
                      HgIndex_t       StartIndex,
                      HgIndex_t       EndIndex,
                      Boolean_t      *IsOk);
EXTERN void ReadClassicOrderedCCBlock(FileStream_s    *DataFileStream,
                                      FieldData_pa     FieldData,
                                      FieldDataType_e  FieldDataTypeInFile,
                                      LgIndex_t        NumIPtsInFile,
                                      LgIndex_t        NumJPtsInFile,
                                      LgIndex_t        NumKPtsInFile,
                                      Boolean_t       *IsOk);
EXTERN Boolean_t ReadInDataFileTypeTitleAndVarNames(FileStream_s   *FileStream,
                                                    short           IVersion,
                                                    char          **DataSetTitle,
                                                    DataFileType_e *FileType,
                                                    int            *NumVars,
                                                    StringList_pa  *VarNames);
EXTERN Boolean_t ReadInZoneHeader(FileStream_s *FileStream,
                                  short         IVersion,
                                  ZoneSpec_s   *ZoneSpec,
                                  Set_pa        IsVarCellCentered,
                                  EntIndex_t    NumVars,
                                  Boolean_t    *IsRawFNAvailable,
                                  LgIndex_t    *FNNumBndryConns);
EXTERN Boolean_t ReadInCustomLabels(FileStream_s  *FileStream,
                                    short          IVersion,
                                    Boolean_t      OkToLoad,
                                    StringList_pa *CustomLabelBase);
EXTERN Boolean_t ReadInUserRec(FileStream_s  *FileStream,
                               short          IVersion,
                               int            MaxCharactersAllowed,
                               char         **UserRec);
EXTERN Boolean_t ReadInAuxData(FileStream_s *FileStream,
                               short         IVersion,
                               AuxData_pa    AuxData);
EXTERN Boolean_t ReadInGeometry(FileStream_s *FileStream,
                                short         IVersion,
                                Boolean_t     OkToLoad,
                                Geom_s       *G,
                                LgIndex_t     MaxDataPts);
EXTERN Boolean_t ReadInText(FileStream_s *FileStream,
                            short         IVersion,
                            Boolean_t     OkToLoad,
                            Text_s       *T,
                            LgIndex_t     MaxTextLen);
/*
 * STDCALL since PreplotAsciiDatafile is sent to RegisterDataSetReader
 * which can also be used by addons.
 */
EXTERN Boolean_t STDCALL PreplotAsciiDatafile(char  *CurFName,
                                              char  *BinaryFName,
                                              char **MessageString);
EXTERN short GetInputVersion(FileStream_s *FileStream);

EXTERN Boolean_t WriteBinaryInt16BlockUnaligned(FileStream_s *FileStream,
                                                Byte_t       *Int16Values,
                                                HgIndex_t     NumValues,
                                                Boolean_t     ValuesInNativeOrdering);
EXTERN Boolean_t WriteBinaryInt32BlockUnaligned(FileStream_s *FileStream,
                                                Byte_t       *Int32Values,
                                                HgIndex_t     NumValues,
                                                Boolean_t     ValuesInNativeOrdering);
EXTERN Boolean_t WriteBinaryByteBlock(FileStream_s    *FileStream,
                                      const Byte_t    *ByteValues,
                                      const HgIndex_t  NumValues);
EXTERN Boolean_t WriteBinaryInt16(FileStream_s *FileStream,
                                  Int16_t       Value);
EXTERN Boolean_t WriteBinaryInt32(FileStream_s *FileStream,
                                  Int32_t       Value);
EXTERN Boolean_t WriteBinaryReal(FileStream_s    *FileStream,
                                 double           RR,
                                 FieldDataType_e  FieldDataType);
EXTERN Boolean_t WriteFieldDataType(FileStream_s    *FileStream,
                                    FieldDataType_e  FDT,
                                    Boolean_t        WriteBinary);
EXTERN Boolean_t WriteBinaryFieldDataBlock(FileStream_s *FileStream,
                                           FieldData_pa  D,
                                           LgIndex_t     StartI,
                                           LgIndex_t     NumValues);
EXTERN Boolean_t WriteCCFieldDataBlock(FileStream_s *FileStream,
                                       FieldData_pa  FieldData,
                                       Boolean_t     IsOrderedData,
                                       LgIndex_t     NumIPts,
                                       LgIndex_t     NumJPts,
                                       LgIndex_t     NumKPts,
                                       Boolean_t     WriteBinary,
                                       SmInteger_t   AsciiPrecision);
EXTERN Boolean_t DumpDatafileString(FileStream_s *FileStream,
                                    const char   *S,
                                    Boolean_t     WriteBinary);
bool DumpGeometry(FileStream_s* FileStream,
                  Geom_s const* Geom,
                  Boolean_t     WriteBinary,
                  Boolean_t     WriteGridDataAsPolar);
bool DumpText(FileStream_s* FileStream,
              Text_s const* Text,
              Boolean_t     WriteBinary,
              Boolean_t     WriteGridDataAsPolar);
EXTERN Boolean_t DumpCustomAxisLabels(FileStream_s  *FileStream,
                                      Boolean_t      WriteBinary,
                                      StringList_pa  LabelBase);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

EXTERN Boolean_t WriteBinaryMagic(FileStream_s *FileStream);

bool writeBinaryVersionNumber(FileStream_s& fileStream,
                              int           versionNumber);

#endif //DATAIO4_H
