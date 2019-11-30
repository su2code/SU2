 #pragma once
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___884
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
EXTERN void ___2889(void);
 #if !defined NO_ASSERTS && !defined NO_DEBUG_FIELDVALUES && !defined DEBUG_FIELDVALUES
 #define DEBUG_FIELDVALUES
 #endif
struct ___1362 { void*                    ___816; ___1383 ___1781; ___1384 ___3506; ___372                isOrderedData; FieldDataType_e          ___4335; ValueLocation_e          ___4327; ___2227 iDim; ___2227 jDim; ___2227 kDim; }; inline void* ___1720(___1361 ___1352) { return ___1352->___816; } inline void ___3487(___1361 ___1352, void* ___3271) { ___1352->___816 = ___3271; } inline ___1383 ___1696(___1361 ___1352) { return ___1352->___1781; } inline ___1384 ___1723(___1361 ___1352) { return ___1352->___3506; } inline ___2227 GetFieldDataIDim(___1361 ___1352) { return ___1352->iDim; } inline ___2227 GetFieldDataJDim(___1361 ___1352) { return ___1352->jDim; } inline ___2227 GetFieldDataKDim(___1361 ___1352) { return ___1352->kDim; } inline ___2227 ___1717(___1361 ___1352) { if (___1352->isOrderedData) { if (___1352->___4327 == ___4330) { return ___1352->iDim * ___1352->jDim * ___1352->kDim; } else if (___1352->iDim == 0 && ___1352->jDim == 0 && ___1352->kDim == 0) { return 0; } else { ___2227 iDim = ___1352->iDim; ___2227 jDim = ___1352->jDim; ___2227 kDim = ___1352->kDim; if (___1352->kDim > 1) --kDim; else if (___1352->jDim > 1) --jDim; else if (___1352->iDim > 1) --iDim; return iDim * jDim * kDim; } } else if (___1352->___4327 == ___4330) { return ___1352->iDim; } else { return ___1352->jDim; } } inline ValueLocation_e ___1729(___1361 ___1352) { return ___1352->___4327; } inline FieldDataType_e ___1726(___1361 ___1352) { return ___1352->___4335; } EXTERN double STDCALL getUniformFieldValueAdapter(___1361 ___1351, ___2227 point); EXTERN double STDCALL ___1742(const ___1361 ___1308, ___2227 ___3249); EXTERN double STDCALL ___1741(const ___1361 ___1308, ___2227 ___3249); EXTERN double STDCALL ___1744(const ___1361 ___1308, ___2227 ___3249); EXTERN double STDCALL ___1743(const ___1361 ___1308, ___2227 ___3249); EXTERN double STDCALL ___1740(const ___1361 ___1308, ___2227 ___3249); EXTERN double STDCALL ___1739(const ___1361 ___1308, ___2227 ___3249); EXTERN ___1383 DetermineFieldDataGetFunction(___1361 ___1351); EXTERN ___1384 DetermineFieldDataSetFunction(___1361 ___1351); inline bool ___2019(___1361 ___1352) { return ___1720(___1352) != NULL; } typedef uint32_t ___1437; typedef uint64_t ___1111; typedef uint16_t ___1961; typedef uint32_t ___1966; typedef uint64_t ___1971; inline float*       ___1690(___1361 ___1308)     { return (float*)      ___1720(___1308); } inline ___1437*  ___1693(___1361 ___1308)  { return (___1437*) ___1720(___1308); } inline double*      ___1684(___1361 ___1308)    { return (double*)     ___1720(___1308); } inline ___1111* ___1687(___1361 ___1308) { return (___1111*)___1720(___1308); } inline int64_t*     ___1711(___1361 ___1308)     { return (int64_t*)    ___1720(___1308); } inline ___1971*  ___1714(___1361 ___1308)  { return (___1971*) ___1720(___1308); } inline int32_t*     ___1705(___1361 ___1308)     { return (int32_t*)    ___1720(___1308); } inline ___1966*  ___1708(___1361 ___1308)  { return (___1966*) ___1720(___1308); } inline int16_t*     ___1699(___1361 ___1308)     { return (int16_t*)    ___1720(___1308); } inline ___1961*  ___1702(___1361 ___1308)  { return (___1961*) ___1720(___1308); } inline uint8_t*     ___1681(___1361 ___1308)      { return (uint8_t*)    ___1720(___1308); } inline uint16_t*    ___1672(___1361 ___1308)     { return (uint16_t*)   ___1720(___1308); } inline uint32_t*    ___1675(___1361 ___1308)     { return (uint32_t*)   ___1720(___1308); } inline uint64_t*    ___1678(___1361 ___1308)     { return (uint64_t*)   ___1720(___1308); } inline void*        ___1732(___1361 ___1308)      { return               ___1720(___1308); } EXTERN ___1361 ___28(___2227       ___2842, FieldDataType_e ___4236, ___372       ___3571); EXTERN void ___938(___1361 *___3449); EXTERN void ___433(___1361  ___1353, double       *min_ptr, double       *max_ptr, ___2227     ___3684, ___1929 *___1928);
EXTERN void ___679( FieldDataType_e ___4335, void*           ___1122, ___2227       ___1127, void const*     ___3657, ___2227       ___3665, ___2227       ___3659); EXTERN void ___3911(FieldDataType_e  ___4335, void            *___3657, ___2227        ___3665, ___2227        ___3659, ___2227        ___3664); EXTERN void ___3912(FieldDataType_e  ___4335, void            *___3657, ___2227        ___3665, ___2227        ___3659, ___2227        ___3664); EXTERN void ___674(___1361 ___1121, ___2227    ___1128, ___1361 ___3656, ___2227    ___3666, ___2227    ___3660); EXTERN void ___673(___1361 ___1121, ___1361 ___3656); EXTERN void ___675(___1361 ___1121, ___2227    ___1124, ___1361 ___3656, ___2227    ___3661); EXTERN void SetFieldDataArrayBytesToZero(___1361 ___1308); EXTERN void ___3486(___1361 ___1353); inline double ___1735(___1361 ___1308, ___2227 ___3249) {
 #if !defined GET_FIELD_VALUE_BY_VIRTUAL_FUNCTION && \
 !defined GET_FIELD_VALUE_BY_FLOAT_ONLY_TEST && \
 !defined GET_FIELD_VALUE_BY_DOUBLE_ONLY_TEST && \
 !defined GET_FIELD_VALUE_BY_FLOAT_AND_DOUBLE_TEST
 #if !defined NO_ASSERTS || defined DEBUG_FIELDVALUES
 #define GET_FIELD_VALUE_BY_VIRTUAL_FUNCTION
 #else
 #define GET_FIELD_VALUE_BY_FLOAT_AND_DOUBLE_TEST
 #endif
 #endif
 #if defined GET_FIELD_VALUE_BY_VIRTUAL_FUNCTION
return ___1696(___1308)(___1308,___3249);
 #elif defined GET_FIELD_VALUE_BY_FLOAT_ONLY_TEST
return ___1696(___1308) == ___1742 ? ___1690(___1308)[___3249] : ___1696(___1308)(___1308,___3249);
 #elif defined GET_FIELD_VALUE_BY_DOUBLE_ONLY_TEST
return ___1696(___1308) == ___1741 ? ___1684(___1308)[___3249] : ___1696(___1308)(___1308,___3249);
 #elif defined GET_FIELD_VALUE_BY_FLOAT_AND_DOUBLE_TEST
return ___1696(___1308) == ___1742 ? ___1690(___1308)[___3249] : (___1696(___1308) == ___1741 ? ___1684(___1308)[___3249] : ___1696(___1308)(___1308,___3249));
 #else
 #error "Need to define one of GET_FIELD_VALUE_BY_XXX constants"
 #endif
} inline void ___3490(___1361 ___1308, ___2227 ___3249, double ___4298) { ___1723(___1308)(___1308,___3249,___4298); }
 #if defined _DEBUG
 #define USEFUNCTIONSFORNODEVALUES
 #endif
EXTERN ___372 ___1355( ___1361 ___1352, ___372    ___3571); EXTERN void ___1358(___1361 ___1352); EXTERN void ___1357(___1361 *___1352, ___372     ___1103);
