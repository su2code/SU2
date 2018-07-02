#if defined EXTERN
#undef EXTERN
#endif
#if defined DATAIOMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN Boolean_t OpenBinaryFileAndCheckMagicNumber(FileStream_s **FileStream,
                                                   char          *FName,
                                                   FileOffset_t   StartOffset,
                                                   short         *IVersion);

EXTERN Boolean_t ReadDataFileHeader(FileStream_s    *FileStream,
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
                                    Set_pa         **IsVarCellCentered,
                                    Boolean_t       *HasText,
                                    Boolean_t       *HasGeoms,
                                    ArrayList_pa    *ZoneSpecList,
                                    StringList_pa   *VarNames,
                                    ArrayList_pa    *VarAuxDataList, /*<AuxData_pa>[NumVars]*/
                                    Set_pa          *IsRawFNAvailable, /* classic data only */
                                    LgIndex_t      **FNNumBndryConns,  /* classic data only */
                                    DataFileType_e  *FileType);


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
