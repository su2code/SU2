#ifndef DATASET_h__
#define DATASET_h__

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
 * DataSet functions involving zones, vars and the
 * DataSet_s structure.  See dataset0.c for low level
 * dataset functions and dataset2 for higher level
 * functions.
 */


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined USE_MACROS_FOR_FUNCTIONS
#else
#endif
#endif /* TECPLOTKERNEL */

Boolean_t FieldDataItemDestructor(void       *ItemRef,
                                  ArbParam_t  ClientData);
Boolean_t ZoneSpecItemDestructor(void       *ItemRef,
                                 ArbParam_t  ClientData);
LgIndex_t ZoneOrVarListAdjustCapacityRequest(ArrayList_pa ZoneOrVarArrayList,
                                             LgIndex_t    CurrentCapacity,
                                             LgIndex_t    RequestedCapacity,
                                             ArbParam_t   ClientData);
void CleanoutZoneSpec(ZoneSpec_s *ZoneSpec);
void ZoneSpecExcludeBndryConnsFromMetrics(ZoneSpec_s* ZoneSpec);
ZoneSpec_s *ZoneSpecAlloc(void);
void ZoneSpecDealloc(ZoneSpec_s **ZoneSpec);
void SetZoneSpecDefaults(ZoneSpec_s *ZoneSpec);

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined USE_MACROS_FOR_FUNCTIONS
#else
#endif
#endif

#define GetZoneSpec(ZoneSpecList,Zone) ((ZoneSpec_s *)ArrayListGetVoidPtr(ZoneSpecList,Zone))
#define GetZoneAuxData(DataSet, Zone) (GetZoneSpec((DataSet)->ZoneSpecList, (Zone))->AuxData)
#define GetVarSpec(VarSpecList,Var) ((VarSpec_s *)ArrayListGetVoidPtr(VarSpecList,Var))
#define GetVarAuxData(DataSet, Var) (GetVarSpec((DataSet)->VarSpecList, (Var))->AuxData)
#define GetStrandInfo(StrandInfoList, StrandID) ((StrandInfo_s *)ArrayListGetVoidPtr(StrandInfoList,StrandID))


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* defined TECPLOTKERNEL */

#endif // DATASET_h__
