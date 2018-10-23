 #ifndef DATAIO4_H
 #define DATAIO4_H
#include <set>
 #if defined EXTERN
 #undef EXTERN
 #endif
 #if defined ___856
 #define EXTERN
 #else
 #define EXTERN extern
 #endif
EXTERN double ___1762(___1405    *___1401, FieldDataType_e  ___1364, double           ___2470, double           ___2346, ___372       *___2040); template <typename SRC_INT_TYPE> EXTERN ___2227 ___1750( ___1405* ___1400, short         pltFileVersion, ___2227     minValue, ___2227     maxValue, ___372*    ___2039); EXTERN void ___3277(___1405 *___1401, ___372     ___1099, uint8_t       *___417, ___2227     ___3685, ___2227     ___2842, ___372    *___2040); EXTERN void ___3290(___1405 *___1401, ___372     ___1099, int16_t      *___417, ___2227     ___3685, ___2227     ___2842, ___372    *___2040); EXTERN void ___3291(___1405 *___1401, ___372     ___1099, int32_t      *___417, ___2227     ___3685, ___2227     ___2842, ___372    *___2040); EXTERN void ___3292(___1405 *___1401, ___372     ___1099, int32_t      *___417, ___2227     ___3685, ___2227     ___2842, ___372    *___2040); EXTERN void ReadInt64Block(___1405 *___1401, ___372     ___1099, int64_t      *___417, ___2227     ___3685, ___2227     ___2842, ___372    *___2040); EXTERN void ___3296(___1405   *___1401, ___372       ___1099, void           *___417, FieldDataType_e ___1364, ___2227       ___3685, ___2227       ___2842, ___372      *___2040); EXTERN void ___3276(___1405   *___1401, ___1361    ___1352, ___372       ___1099, FieldDataType_e ___1367, ___2227       ___3685, ___2227       EndIndex, ___372      *___2040); EXTERN void ___3278(___1405    *___843, ___1361     ___1352, FieldDataType_e  ___1367, ___2227        ___2810, ___2227        ___2815, ___2227        ___2818, ___372       *___2040); EXTERN ___372 ___3287(___1405   *___1401, short           ___2104, char          **___903, DataFileType_e *___1408, int            *NumVars, ___3839  *___4366); EXTERN ___372 ___3295(___1405 *___1401, short         ___2104, ___4683   *___4677, ___3501        ___2075, ___1172    NumVars, ___90    ___263, ___372    *___2053, ___2227    *___1441); EXTERN ___372 ___3286(___1405  *___1401, short          ___2104, ___372      ___2868, ___3839 *___791); EXTERN ___372 ___3294(___1405  *___1401, short          ___2104, int            ___2352, char         **___4285); EXTERN ___372 ___3285(___1405 *___1401, short         ___2104, ___264    ___230); EXTERN ___372 ___3288(___1405 *___1401, short         ___2104, ___372     ___2868, ___1632       *G, ___2227     ___2366); EXTERN ___372 ___3293(___1405 *___1401, short         ___2104, ___372     ___2868, ___4118       *T, ___2227     ___2387); EXTERN ___372 STDCALL ___3172(char  *___692, char  *___359, char **___2456); EXTERN short ___1749(___1405 *___1401); EXTERN ___372 ___4489( ___1405*  ___1401, uint8_t const* ___1963, ___2227      ___2842, ___372      ___4333); EXTERN ___372 ___4491( ___1405*  ___1401, uint8_t const* ___1968, ___2227      ___2842, ___372      ___4333); EXTERN ___372 WriteBinaryInt64BlockUnaligned( ___1405*  ___1401, uint8_t const* Int64Values, ___2227      ___2842, ___372      ___4333); EXTERN ___372 ___4486( ___1405*  ___1401, uint8_t const* ___429, ___2227      ___2842); EXTERN ___372 ___4488(___1405 *___1401, int16_t       ___4315); EXTERN ___372 ___4490(___1405 *___1401, int32_t       ___4315); EXTERN ___372 WriteBinaryInt64(___1405 *___1401, int64_t       ___4315); EXTERN ___372 ___4493(___1405    *___1401, double           ___3425, FieldDataType_e  ___1364); EXTERN ___372 ___4513(___1405    *___1401, FieldDataType_e  ___1310, ___372        ___4485); EXTERN ___372 ___4487(___1405 *___1401, ___1361  D, ___2227     ___3683, ___2227     ___2842); EXTERN ___372 ___4495(___1405 *___1401, ___1361  ___1352, ___372     ___2043, ___2227     ___2809, ___2227     ___2814, ___2227     ___2817, ___372     ___4485, int32_t   ___200); EXTERN ___372 ___1131(___1405 *___1401, const char   *S, ___372     ___4485); bool ___1132(___1405* ___1401, ___1632 const* ___1556, ___372     ___4485, ___372     ___4525); bool ___1133(___1405* ___1401, ___4118 const* Text, ___372     ___4485,
___372     ___4525); EXTERN ___372 ___1130(___1405  *___1401, ___372      ___4485, ___3839  ___2170); EXTERN ___372 ___4492(___1405 *___1401); bool ___4494(___1405& ___1400, int           ___4409);
 #endif 
