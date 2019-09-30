 #ifndef _GLOBAL_H
 #define _GLOBAL_H
 #if !defined _MASTER_H_ && defined TECPLOTKERNEL
 #error "Must include MASTER.h before including GLOBAL.h"
 #endif
#include "StandardIntegralTypes.h"
#include <limits>
 #if defined EXTERN
 #undef EXTERN
 #endif
 #define EXTERN extern
 #define EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
 #if defined ___4226
 #undef ___4226
 #endif
 #if defined ___1305
 #undef ___1305
 #endif
 #if defined MIN
 #undef MIN
 #endif
 #if defined MAX
 #undef MAX
 #endif
 #if defined ___3419
 #undef ___3419
 #endif
 #if defined ___3420
 #undef ___3420
 #endif
 #if defined ___4227
 #undef ___4227
 #endif
 #define ___4226                  ((___372)1)
 #define ___1305                 ((___372)0)
 #if defined MSWIN
 #define STDABSINT64 _abs64
 #else
 #define STDABSINT64 std::abs
 #endif
 #define ___1(X)                ((X) >= 0 ? (X) : -(X) )
 #define MAX(X,Y)              ((X) > (Y) ? (X) : (Y) )
 #define MIN(X,Y)              ((X) < (Y) ? (X) : (Y) )
 #define ___352(X)      ((X) == ___4455 ? ___364 : ___4455)
 #define ___3423(X)      ((BYTE)((X)+0.499))
 #define ___3422(X)             ((short)((X)+0.499))
 #define ROUND32(X)            ((int32_t)((X)+0.499))
 #define ___3421(X)             ((___2227)((X)+0.499))
 #define ___3420(X)             ((X) >= 0 ? ((int)((X)+0.499)) : ((int)((X)-0.499)))
 #define ___4227(X)              ((short) (X))
 #define ___3267(___3266)       (180.*(___3266)/___3003)
 #define ___955(___953)       (___3003*(___953)/180.)
 # define ___444(C) ((char)(('a'<=(C)&&(C)<='z') ? ((C)+('A'-'a')) : (C))) 
 #define ___2017(S)      ( ((const char*)(S))[0] == '\0' )
 #define ___2082(C)       ((C == ' ') || (C == '\t') || (C == '\n'))
 #define ___2055(C)        ((C == ' ') || (C == '\t') || (C == ','))
 #define ___488(___4314,___2324,___1830) ((___4314)<(___2324) ? (___2324) : (___4314) > (___1830) ? (___1830) : (___4314))
 #define ___1974(n, d) (((int)(n)+(int)(d)-1)/(int)(d))
 #define ___1854(___801,___1832,___2106,___2135) ((___1832) + \
 ((___2106)*(___801)->___2809) + \
 ((___2135)*(___801)->___2807))
 #define ___1842(___801,N) ((N) % (___801)->___2809)
 #define ___2112(___801,N) (((N) % (___801)->___2807)/(___801)->___2809)
 #define ___2157(___801,N) ((N)/(___801)->___2807)
 #define ___1927(index,___3603) ((___3603) == 1 || ((index) % (___3603) == 0))
 #define ___3913(___4236,A,B)      do {___4236 T = (A); (A) = (B); (B) = T;} while (___1305)
 #define ___3914(A,B)   ___3913(double, (A), (B))
 #define ___1482(x)          (___372)((x) > 0)
 #define ___1807(F)      (((F)->___3112 == ___3115))
 #if defined ___4278
 #undef ___4278
 #endif
 #define ___4278(___2971) (void)___2971
 #define ___650(___4298) \
 ( (___4298) >= ___3632 \
 ? ( (___4298) < ___2180 \
 ? (float)(___4298) \
 : (float)___2180 \
 ) \
 : ( (___4298) <= -___3632  \
 ? ( (___4298) > -___2180 \
 ? (float)(___4298) \
 : (float)-___2180 \
 ) \
 : (float)0.0 \
 ) \
 )
 #define ___489(___4298) \
 ( (___4298) >= ___3628 \
 ? ( (___4298) < ___2179 \
 ? (double)(___4298) \
 : (double)___2179 \
 ) \
 : ( (___4298) <= -___3628  \
 ? ( (___4298) > -___2179 \
 ? (double)(___4298) \
 : (double)-___2179 \
 ) \
 : (double)0.0 \
 ) \
 )
 #define ___652(___4298) \
 ( (___4298) >= 1.0 \
 ? ( (___4298) < ___2182 \
 ? (int32_t)(___4298) \
 : (int32_t)___2182 \
 ) \
 : ( (___4298) <= -1.0  \
 ? ( (___4298) > (int32_t)-___2182 \
 ? (int32_t)(___4298) \
 : (int32_t)-___2182 \
 ) \
 : (int32_t)0.0 \
 ) \
 )
 #define ___651(___4298) \
 ( (___4298) >= 1.0 \
 ? ( (___4298) < ___2181 \
 ? (int16_t)(___4298) \
 : (int16_t)___2181 \
 ) \
 : ( (___4298) <= -1.0  \
 ? ( (___4298) > (int16_t)-___2181 \
 ? (int16_t)(___4298) \
 : (int16_t)-___2181 \
 ) \
 : (int16_t)0.0 \
 ) \
 )
 #define CONVERT_DOUBLE_TO_UINT8(___4298) \
 ( (___4298) >= 1.0 \
 ? ( (___4298) < 255.0 \
 ? (uint8_t)(___4298) \
 : 255 \
 ) \
 : 0 \
 )
 #define ___667(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[0]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[1]; \
 } while (___1305)
 #define ___670(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[1]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[0]; \
 } while (___1305)
 #define ___668(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[0]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[1]; \
 ((uint8_t *)(___1123))[2] = ((uint8_t *)(___3658))[2]; \
 ((uint8_t *)(___1123))[3] = ((uint8_t *)(___3658))[3]; \
 } while (___1305)
 #define ___671(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[3]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[2]; \
 ((uint8_t *)(___1123))[2] = ((uint8_t *)(___3658))[1]; \
 ((uint8_t *)(___1123))[3] = ((uint8_t *)(___3658))[0]; \
 } while (___1305)
 #define ___669(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[0]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[1]; \
 ((uint8_t *)(___1123))[2] = ((uint8_t *)(___3658))[2]; \
 ((uint8_t *)(___1123))[3] = ((uint8_t *)(___3658))[3]; \
 ((uint8_t *)(___1123))[4] = ((uint8_t *)(___3658))[4]; \
 ((uint8_t *)(___1123))[5] = ((uint8_t *)(___3658))[5]; \
 ((uint8_t *)(___1123))[6] = ((uint8_t *)(___3658))[6]; \
 ((uint8_t *)(___1123))[7] = ((uint8_t *)(___3658))[7]; \
 } while (___1305)
 #define ___672(___1123, ___3658) \
 do { \
   \
   \
 ((uint8_t *)(___1123))[0] = ((uint8_t *)(___3658))[7]; \
 ((uint8_t *)(___1123))[1] = ((uint8_t *)(___3658))[6]; \
 ((uint8_t *)(___1123))[2] = ((uint8_t *)(___3658))[5]; \
 ((uint8_t *)(___1123))[3] = ((uint8_t *)(___3658))[4]; \
 ((uint8_t *)(___1123))[4] = ((uint8_t *)(___3658))[3]; \
 ((uint8_t *)(___1123))[5] = ((uint8_t *)(___3658))[2]; \
 ((uint8_t *)(___1123))[6] = ((uint8_t *)(___3658))[1]; \
 ((uint8_t *)(___1123))[7] = ((uint8_t *)(___3658))[0]; \
 } while (___1305)
 #define ___3365(___417) \
 do { \
 uint8_t ___424 = ((uint8_t *)(___417))[0]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==2); \
 ((uint8_t *)(___417))[0] = ((uint8_t *)(___417))[1]; \
 ((uint8_t *)(___417))[1] = ___424; \
 } while (___1305)
 #define ___3366(___417) \
 do { \
 uint16_t ___828 = ((uint16_t *)(___417))[0]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==2); \
 ((uint16_t *)(___417))[0] = (((___828)<<8) | \
 ((___828&0xff))); \
 } while (___1305)
 #define ___3364 ___3365
 #define ___3368(___417) \
 do { \
 uint8_t ___424 = ((uint8_t *)(___417))[0]; \
 uint8_t ___425 = ((uint8_t *)(___417))[1]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==4); \
 ((uint8_t *)(___417))[0] = ((uint8_t *)(___417))[3]; \
 ((uint8_t *)(___417))[1] = ((uint8_t *)(___417))[2]; \
 ((uint8_t *)(___417))[2] = ___425; \
 ((uint8_t *)(___417))[3] = ___424; \
 } while (___1305)
 #define ___3369(___417) \
 do { \
 uint32_t ___828 = *((uint32_t *)(___417)); \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==4); \
 *((uint32_t *)(___417)) = (((___828)<<24)            | \
 ((___828&0x0000ff00)<<8)  | \
 ((___828&0x00ff0000)>>8)  | \
 ((___828)>>24)); \
 } while (___1305)
 #if defined MSWIN
 #define ___3367 ___3369
 #else
 #define ___3367 ___3368
 #endif
 #define ___3371(___417) \
 do { \
 uint8_t ___424 = ((uint8_t *)(___417))[0]; \
 uint8_t ___425 = ((uint8_t *)(___417))[1]; \
 uint8_t ___426 = ((uint8_t *)(___417))[2]; \
 uint8_t ___427 = ((uint8_t *)(___417))[3]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==8); \
 ((uint8_t *)(___417))[0] = ((uint8_t *)(___417))[7]; \
 ((uint8_t *)(___417))[1] = ((uint8_t *)(___417))[6]; \
 ((uint8_t *)(___417))[2] = ((uint8_t *)(___417))[5]; \
 ((uint8_t *)(___417))[3] = ((uint8_t *)(___417))[4]; \
 ((uint8_t *)(___417))[4] = ___427; \
 ((uint8_t *)(___417))[5] = ___426; \
 ((uint8_t *)(___417))[6] = ___425; \
 ((uint8_t *)(___417))[7] = ___424; \
 } while (___1305)
 #define ___3372(___417) \
 do { \
 uint16_t ___829 = ((uint16_t *)(___417))[0]; \
 uint16_t ___830 = ((uint16_t *)(___417))[1]; \
 uint16_t ___831 = ((uint16_t *)(___417))[2]; \
 uint16_t ___832 = ((uint16_t *)(___417))[3]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==8); \
 ((uint16_t *)(___417))[0] = (((___832)<<8) | \
 ((___832&0xff))); \
 ((uint16_t *)(___417))[1] = (((___831)<<8) | \
 ((___831&0xff))); \
 ((uint16_t *)(___417))[2] = (((___830)<<8) | \
 ((___830&0xff))); \
 ((uint16_t *)(___417))[3] = (((___829)<<8) | \
 ((___829&0xff))); \
 } while (___1305)
 #define ___3373(___417) \
 do { \
 uint32_t ___829 = ((uint32_t *)(___417))[0]; \
 uint32_t ___830 = ((uint32_t *)(___417))[1]; \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==8); \
 ((uint32_t *)(___417))[0] = (((___830)<<24)           | \
 ((___830&0x0000ff00)<<8) | \
 ((___830&0x00ff0000)>>8) | \
 ((___830)>>24)); \
 ((uint32_t *)(___417))[1] = (((___829)<<24)           | \
 ((___829&0x0000ff00)<<8) | \
 ((___829&0x00ff0000)>>8) | \
 ((___829)>>24)); \
 } while (___1305)
 #define ___3374(___417) \
 do { \
 uint64_t ___828 = *((uint64_t *)(___417)); \
 ___478(sizeof(*(___417))==1 || sizeof(*(___417))==8); \
 *((uint64_t *)(___417)) = (((___828)<<56) | \
 ((___828&0x000000000000ff00)<<40) | \
 ((___828&0x0000000000ff0000)<<24) | \
 ((___828&0x00000000ff000000)<<8)  | \
 ((___828&0x000000ff00000000)>>8)  | \
 ((___828&0x0000ff0000000000)>>24) | \
 ((___828&0x00ff000000000000)>>40) | \
 ((___828)>>56)); \
 } while (___1305)
 #if defined MSWIN
 #define ___3370 ___3373
 #else
 #define ___3370 ___3371
 #endif
 #if defined MSWIN
 #   define STDCALL __stdcall
 #else
 #   define STDCALL
 #endif 
 #if defined (__cplusplus)
 #   define EXTERNC extern "C"
 #   define TP_GLOBAL_NAMESPACE ::
 #else
 #   define EXTERNC
 #   define TP_GLOBAL_NAMESPACE
 #endif 
 #if defined MAKEARCHIVE
 #define tpsdkbase_API
 #else
#include "tpsdkbase_Exports.h"
 #endif
 #define ___2291 EXTERNC tpsdkbase_API
 #if defined MSWIN
 # define ___1235 EXTERNC _declspec ( dllexport )
 #else
 # define ___1235 EXTERNC
 #endif 
 #define ___1236 ___1235
 #define ___1939           ___1940
 #define ___3969 "InitTecAddOn113"
 #if defined (__cplusplus) && !defined _DEBUG
 # define ___1941 inline
 #else
 # define ___1941 static
 #endif 
 #if defined (MSWIN) ||\
 defined (___1985) ||\
 defined (LINUX) ||\
 defined (___3891) ||\
 defined (___534) ||\
 defined (DEC) ||\
 defined (__LITTLE_ENDIAN__)
 #define ___2326
 #endif
 #if defined( ___2326 )
 # define ___3910(___1977) (!(___1977))
 #else
 # define ___3910(___1977) (___1977)
 #endif
 #if defined DECALPHA   || \
 defined LINUXALPHA || \
 defined LINUX64    || \
 defined MAC64      || \
 defined ___1833  || \
 defined ___3886        || \
 defined ___1831         || \
 defined ___534
 #define ___2320
 #endif
 #define ___2185              ((size_t)-1)
 #define ___2183               (std::numeric_limits<int64_t>::max()-1)
 #define ___2182               (std::numeric_limits<int32_t>::max()-1)
 #define ___2181               (std::numeric_limits<int16_t>::max()-1)
 #define ___2184                (std::numeric_limits<int8_t>::max()-1)
 #define ___2193              (std::numeric_limits<uint64_t>::max()-1)
 #define ___2192              (std::numeric_limits<uint32_t>::max()-1)
 #define ___2191              (std::numeric_limits<uint16_t>::max()-1)
 #define ___2194               (std::numeric_limits<uint8_t>::max()-1)
 #ifdef INDEX_16_BIT
 #define ___2373                 ((___2227)___2181)
 #else
 #define ___2373                 ((___2227)(___2183>>1))
 #endif
 #define ___2391               ((___1172)___2182)
 #define ___2179              1.0e+150
 #define ___3628              1.0e-150
 #define ___2190          150
 #define ___3631         -150
 #define ___3629           ___3628
 #define ___2189    308
 #define ___3630   -307
 #define ___2188            1.0e+308
 #define ___2180               std::numeric_limits<float>::max() 
 #define ___3632               std::numeric_limits<float>::min() 
 #define ___3633            1.0e-307
 #define ___1189                      3
 #define ___2296                      0.69314718055994530942
 #define ___2295                     2.30258509299404568402
 #define ___3088                  1.57079632679489661923
 #define ___4235                    6.28318530717958647692
 #if defined ___3003
 #undef ___3003
 #endif
 #define ___3003                       3.14159265358979323846
 #define ___58             1.0e-10
 #define ___2187             (4*___3003+___58)
 #define ___954            57.295779513082323
 #define ___506                2.54
 #define ___3143            72.0
 #define ___1460             192
 #define ___1447         128
 #define ___1458             64
 #define ___333            ((___3493)-1)
 #define ___334           (static_cast<___1172>(-1))
 #define ___2420      (0)
 #define ___2421       (-1)
 #define ___1991        0
 #define NO_EMBEDDED_LPK_IMAGE_NUMBER 0
 #define ___332              ___333
 #define ___3637       0
 #define ___328             (-1.0)
 #define ___2473  4
 #define ___4312(___3786) (0 <= (___3786) && (___3786) < ___2391)
 #define ___3788          (-1)
 #define ___3787         (-2)
 #define ___2347 1
 #define ___2348 6
 #define ___1987 -1
 #define ___4301(clipPlane) (0 <= clipPlane && clipPlane < ___2348)
 #define    ___2382           ___2391
 #define    ___2389                    5
 #define    ___2390                    5
 #define    ___2371              50
 #define    ___2384       720
 #define    ___2370                   2048
 #define    ___2365          10
 #define    ___2369                20000
 #define    ___2362        16
 #define    ___2383           16
 #define    ___2472      12
 #define    ___2359      360
 #define    ___2388    8
 #define    ___2363            8
 #define    ___2374         8
 #define    ___2375 3
 #define    ___2386              8
 #define    ___952      15
 #define    ___948        ((int32_t)0)
 #define    ___947      ((int32_t)0)
 #define    ___951 ((int32_t)0)
 #define    ___331             ((int32_t)-1)
 #define    ___4279          ((int32_t)0)
 #define    ___327          ((int32_t)-1)
 #define ___4307(___1817) (((((int32_t)___1817) >= 0) && (((int32_t)___1817) < ___2374)))
 #define ___4311(___1817)      (((((int32_t)___1817) >= 0) && (((int32_t)___1817) < ___2386)))
 #define    ___2349  6
 #define    ___2355       256
 #define    ___2358            128
 #define    ___2356            128
 #define    ___2357        128
 #define    ___2353     32000
 #define    ___2354       1024
 #define    ___2381               16
 #define    ___2350             5
 #define    ___2360  50
 #define    ___2385     800
 #define    ___2367         100
 #define    ___2368      100
 #define    ___2361         20
 #define    ___2482     0.5
 #define    ___2476             0.25
 #define    ___2475            0.25
 #define    ___2471             0.1
 #define    ___329              255
 #define MAX_NODES_PER_CLASSIC_FACE    4
 #define MAX_NODES_PER_CLASSIC_ELEMENT 8
 #define ___3868 16
 #define AuxData_Common_Incompressible               "Common.Incompressible"
 #define AuxData_Common_Density                      "Common.Density"
 #define AuxData_Common_SpecificHeat                 "Common.SpecificHeat"
 #define AuxData_Common_SpecificHeatVar              "Common.SpecificHeatVar"
 #define AuxData_Common_GasConstant                  "Common.GasConstant"
 #define AuxData_Common_GasConstantVar               "Common.GasConstantVar"
 #define AuxData_Common_Gamma                        "Common.Gamma"
 #define AuxData_Common_GammaVar                     "Common.GammaVar"
 #define AuxData_Common_Viscosity                    "Common.Viscosity"
 #define AuxData_Common_ViscosityVar                 "Common.ViscosityVar"
 #define AuxData_Common_Conductivity                 "Common.Conductivity"
 #define AuxData_Common_ConductivityVar              "Common.ConductivityVar"
 #define AuxData_Common_AngleOfAttack                "Common.AngleOfAttack"
 #define AuxData_Common_SpeedOfSound                 "Common.SpeedOfSound"
 #define AuxData_Common_ReferenceU                   "Common.ReferenceU"
 #define AuxData_Common_ReferenceV                   "Common.ReferenceV"
 #define AuxData_Common_XVar                         "Common.XVar"
 #define AuxData_Common_YVar                         "Common.YVar"
 #define AuxData_Common_ZVar                         "Common.ZVar"
 #define AuxData_Common_CVar                         "Common.CVar"
 #define AuxData_Common_UVar                         "Common.UVar"
 #define AuxData_Common_VVar                         "Common.VVar"
 #define AuxData_Common_WVar                         "Common.WVar"
 #define AuxData_Common_VectorVarsAreVelocity        "Common.VectorVarsAreVelocity"
 #define AuxData_Common_PressureVar                  "Common.PressureVar"
 #define AuxData_Common_TemperatureVar               "Common.TemperatureVar"
 #define AuxData_Common_DensityVar                   "Common.DensityVar"
 #define AuxData_Common_StagnationEnergyVar          "Common.StagnationEnergyVar"
 #define AuxData_Common_MachNumberVar                "Common.MachNumberVar"
 #define AuxData_Common_ReferenceMachNumber          "Common.ReferenceMachNumber"
 #define AuxData_Common_ReferenceW                   "Common.ReferenceW"
 #define AuxData_Common_PrandtlNumber                "Common.PrandtlNumber"
 #define AuxData_Common_Axisymmetric                 "Common.Axisymmetric"
 #define AuxData_Common_AxisOfSymmetryVarAssignment  "Common.AxisOfSymmetryVarAssignment"
 #define AuxData_Common_AxisValue                    "Common.AxisValue"
 #define AuxData_Common_SteadyState                  "Common.SteadyState"
 #define AuxData_Common_TurbulentKineticEnergyVar    "Common.TurbulentKineticEnergyVar"
 #define AuxData_Common_TurbulentDissipationRateVar  "Common.TurbulentDissipationRateVar"
 #define AuxData_Common_TurbulentViscosityVar        "Common.TurbulentViscosityVar"
 #define AuxData_Common_TurbulentFrequencyVar        "Common.TurbulentFrequencyVar"
 #define AuxData_Common_Gravity                      "Common.Gravity"
 #define AuxData_Common_IsBoundaryZone               "Common.IsBoundaryZone"
 #define AuxData_Common_BoundaryCondition            "Common.BoundaryCondition"
 #define AuxData_Common_Time                         "Common.Time"
 #define AuxData_Common_Mean                         "Common.Mean"
 #define AuxData_Common_Median                       "Common.Median"
 #define AuxData_Common_Variance                     "Common.Variance"
 #define AuxData_Common_StdDev                       "Common.StdDev"
 #define AuxData_Common_AvgDev                       "Common.AvgDev"
 #define AuxData_Common_GeoMean                      "Common.GeoMean"
 #define AuxData_Common_ChiSqre                      "Common.ChiSqre"
 #define    ___364           ((___516)0)
 #define    ___3301             ((___516)1)
 #define    ___1810           ((___516)2)
 #define    ___366            ((___516)3)
 #define    ___799            ((___516)4)
 #define    ___4587          ((___516)5)
 #define    ___3256          ((___516)6)
 #define    ___4455           ((___516)7)
 #define    ___745         ((___516)8)
 #define    ___756         ((___516)9)
 #define    ___767         ((___516)10)
 #define    ___778         ((___516)11)
 #define    ___786         ((___516)12)
 #define    ___787         ((___516)13)
 #define    ___788         ((___516)14)
 #define    ___789         ((___516)15)
 #define    ___790         ((___516)16)
 #define    ___735         ((___516)17)
 #define    ___736         ((___516)18)
 #define    ___737         ((___516)19)
 #define    ___738         ((___516)20)
 #define    ___739         ((___516)21)
 #define    ___740         ((___516)22)
 #define    ___741         ((___516)23)
 #define    ___742         ((___516)24)
 #define    ___743         ((___516)25)
 #define    ___744         ((___516)26)
 #define    ___746         ((___516)27)
 #define    ___747         ((___516)28)
 #define    ___748         ((___516)29)
 #define    ___749         ((___516)30)
 #define    ___750         ((___516)31)
 #define    ___751         ((___516)32)
 #define    ___752         ((___516)33)
 #define    ___753         ((___516)34)
 #define    ___754         ((___516)35)
 #define    ___755         ((___516)36)
 #define    ___757         ((___516)37)
 #define    ___758         ((___516)38)
 #define    ___759         ((___516)39)
 #define    ___760         ((___516)40)
 #define    ___761         ((___516)41)
 #define    ___762         ((___516)42)
 #define    ___763         ((___516)43)
 #define    ___764         ((___516)44)
 #define    ___765         ((___516)45)
 #define    ___766         ((___516)46)
 #define    ___768         ((___516)47)
 #define    ___769         ((___516)48)
 #define    ___770         ((___516)49)
 #define    ___771         ((___516)50)
 #define    ___772         ((___516)51)
 #define    ___773         ((___516)52)
 #define    ___774         ((___516)53)
 #define    ___775         ((___516)54)
 #define    ___776         ((___516)55)
 #define    ___777         ((___516)56)
 #define    ___779         ((___516)57)
 #define    ___780         ((___516)58)
 #define    ___781         ((___516)59)
 #define    ___782         ((___516)60)
 #define    ___783         ((___516)61)
 #define    ___784         ((___516)62)
 #define    ___785         ((___516)63)
 #define    ___2662      ((___516)(-1))
 #define    ___2708         ((___516)(-2))
 #define    ___2655     ((___516)(-3))
 #define    ___2656     ((___516)(-4))
 #define    ___2657     ((___516)(-5))
 #define    ___3375        ((___516)(-6))
 #define    ___2658     ((___516)(-7))
 #define    ___2659     ((___516)(-8))
 #define    ___2660     ((___516)(-9))
 #define    ___2661     ((___516)(-10))
 #define    ___1988    ((___516)(-255))
 #define    ___1422  ___745
 #define    ___2196   ___785
 #define    ___2791   (___2196-___1422+1)
 #define    ___1420   ___364
 #define    ___2195    ___2196
 #define    ___2766    (___2195-___1420+1)
 #define    ___2872   ((___516)255)
 #define    ___1808                 (___2195+1)
 #define    ___815             (___2195+2) 
 #define    ___4571             (___2195+3)
 #define    ___1423    ___1808
 #define    ___2200     ___4571
 #define    ___2808    (___2200-___1423+1)
 #define    ___2788      (___1547.___2241.___2379+1)
 #define    ___2785 (___2766+___2808+___2788)
 #define    ___342      (0)
 #define    ___1981  (___2766)
 #define    ___614    (___2766+___2808)
 #define    ___2142           (short)31
 #define    ___2151           (short)13
 #define    ___2146              (short)27
 #define    ___2141        (short)8
 #define    ___2143        (short)127
 #define    ___2147        (short)29
 #define    ___2152       (short)30
 #define    ___2155          (short)11
 #define    ___2145        (short)10
 #define    ___2149             (short)43
 #define    ___2148            (short)45
 #define    ___2153      (short)128 
 #define    ___2154           (short)19  
 #define    ___2150         (short)18  
 #define    ___2144        (short)2   
 #define ___4649        299.0
 #define ___1617        399.0
 #define ___4112        499.0
 #define ___792 599.0
 #define ___4286     699.0
 #define ___887  799.0
 #define ___4339      899.0
 #define EndHeaderMarker   357.0
 #define    ___293          ___789+1
 #define    ___2334       ___789+2
 #define    ___2477       ___789+3
 #define    ___2339    ___789+4
 #define    ___3797    ___789+5
 #define    ___515   ___789+6
 #define    ___392      ___789+7
 #define    ___2172         ___789+8
 #define    ___2826   ___789+9
 #define    ___229    ___789+10
 #define    ___1990       ___789+99
 #define    ___1424   ___293
 #define    ___2202    ___2172
 #define    ___960           0.0001
 #define    ___326         NULL
 #ifdef MSWIN
 # define ___1090 "\\"
 #else
 # define ___1090 "/"
 #endif
 #define ___4198  fread
 #define ___4202 fwrite
 #if defined MSWIN
 #if !defined snprintf
 #define snprintf _snprintf
 #endif
 #define ___4196   fflush
 #define ___4195   fclose
 #define ___4206  TP_GLOBAL_NAMESPACE remove
 #define ___4204   TP_GLOBAL_NAMESPACE ___3399
 #define ___4205    TP_GLOBAL_NAMESPACE _stat
 #define ___4203  TP_GLOBAL_NAMESPACE getenv
 #if defined _WIN64
 #define ___4200(___3792,___2865,whence) TP_GLOBAL_NAMESPACE _fseeki64((___3792),(__int64)(___2865),(whence))
 #define ___4201                       TP_GLOBAL_NAMESPACE _ftelli64
 #else
 #define ___4200(___3792, ___2865, whence) TP_GLOBAL_NAMESPACE fseek((___3792), (long)(___2865), (whence))
 #define ___4201                         TP_GLOBAL_NAMESPACE ftell
 #endif
 #else
 #define FileStat_s struct _stat
 #define ___4204  TP_GLOBAL_NAMESPACE ___3398
 #define ___4206 TP_GLOBAL_NAMESPACE unlink
 #define ___4195 TP_GLOBAL_NAMESPACE fclose
 #define ___4196 TP_GLOBAL_NAMESPACE fflush
 #define ___4200  TP_GLOBAL_NAMESPACE fseeko
 #define ___4201  TP_GLOBAL_NAMESPACE ftello
 #define ___4205   TP_GLOBAL_NAMESPACE stat
 #define _stat     TP_GLOBAL_NAMESPACE stat 
 #define ___4203 TP_GLOBAL_NAMESPACE getenv
 #endif
typedef    ___2227 ___2732;
 #if defined CRAY
typedef char *___90;
 #elif defined ___2320
typedef long ___90;
 #elif defined MSWIN
typedef INT_PTR ___90;
 #else
typedef int ___90;
 #endif
typedef ___90 ___4264; typedef int64_t ___1398; typedef uint64_t ___2405; typedef    int32_t          ___516; typedef    int16_t          ___3881; typedef    char             ___372; typedef    char            *___4654; typedef    char            *___4367; typedef    char            *___2325; typedef    ___2227        ___1825; typedef    ___2227        ___3461[___2371]; typedef    double           BasicSize_t[___2350]; typedef    double          *___4355; typedef    ___2227        ___3493; typedef    unsigned long ___3480; typedef    ___3480 ___3483; typedef    ___3483* ___3481;
 #define ___1313                               (1L << 0) 
 #define ___1314                         (1L << 1)
 #define ___1312                               (1L << 2) 
 #define ___1338                               (1L << 3) 
 #define ___1316                        (1L << 4)
 #define ___1337                 (1L << 5)
 #define FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_1  (1L << 6)
 #define FEATURE_NUMBER_OF_ZONES_GREATER_THAN_1   (1L << 7)
 #define FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_5  (1L << 8)
 #define FEATURE_NUMBER_OF_ZONES_GREATER_THAN_5   (1L << 9)
 #define FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_10 (1L << 10)
 #define FEATURE_NUMBER_OF_ZONES_GREATER_THAN_10  (1L << 11)
 #define FEATURE_READ_NONOEM_TECPLOT_DATA         (1L << 12) 
 #define FEATURE_FOREIGN_DATALOADERS              (1L << 13) 
 #define ___1318            (1L << 14) 
 #define ___1324                     (1L << 15) 
 #define ___1329            (1L << 16) 
 #define ___1323                 (1L << 17) 
 #define ___1334                      (1L << 18) 
 #define ___1335                (1L << 19) 
 #define ___1321                     (1L << 20) 
 #define ___1320                        (1L << 21) 
 #define ___1336                          (1L << 22) 
 #define ___1319                    (1L << 23) 
 #define ___1315                          (1L << 24) 
 #define FEATURE_DATASETSIZE                      (1L << 25) 
 #define FEATURE_RPC                              (1L << 26) 
 #define FEATURE_DATALOADERS_EXCEPT_ALLOWED       (1L << 27) 
 #define FEATURE_NUMBER_OF_PAGES_GREATER_THAN_1   (1L << 28) 
 #define FEATURE_BATCH_MODE                       (1L << 29) 
 #define FEATURE_SIMPLEZONECREATION               (1L << 30) 
 #define NUM_POSSIBLE_INHIBITED_FEATURES 31
 #define ___2163 (___1313                               |\
 ___1314                         |\
 ___1312                               |\
 ___1338                               |\
 ___1316                        |\
 ___1337                 |\
 FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_1  |\
 FEATURE_NUMBER_OF_ZONES_GREATER_THAN_1   |\
 FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_5  |\
 FEATURE_NUMBER_OF_ZONES_GREATER_THAN_5   |\
 FEATURE_NUMBER_OF_FRAMES_GREATER_THAN_10 |\
 FEATURE_NUMBER_OF_ZONES_GREATER_THAN_10  |\
 FEATURE_READ_NONOEM_TECPLOT_DATA         |\
 FEATURE_FOREIGN_DATALOADERS              |\
 ___1318            |\
 ___1324                     |\
 ___1329            |\
 ___1323                 |\
 ___1334                      |\
 ___1335                |\
 ___1321                     |\
 ___1320                        |\
 ___1336                          |\
 ___1319                    |\
 ___1315                          |\
 FEATURE_DATASETSIZE                      |\
 FEATURE_RPC                              |\
 FEATURE_DATALOADERS_EXCEPT_ALLOWED       |\
 FEATURE_NUMBER_OF_PAGES_GREATER_THAN_1   |\
 FEATURE_BATCH_MODE                       |\
 FEATURE_SIMPLEZONECREATION)
 #define ___4302(___1311) (((___1311) & ___2163) != 0)
 #define ___4303(___2344) (((___2344) & ~___2163)==0)
typedef    uint64_t  ___1322; typedef    uint64_t  ___1325;
 #define MAX_UTF8_BYTES 4
 #define ___3916 1 + MAX_UTF8_BYTES + 1
typedef    char             ___3917[___3916]; typedef int32_t ___1295; typedef int32_t ___1146; typedef int32_t ___1257; typedef enum { Projection_Orthographic, Projection_Perspective, END_Projection_e, Projection_Invalid = ___329 } Projection_e; typedef enum { ___3091, ___3092, ___3093, END_PlacementPlaneOrientation_e, ___3090 = ___329 } PlacementPlaneOrientation_e; typedef enum { ___3849, ___3852, ___3850, END_StringMode_e, ___3851 = ___329 } StringMode_e; typedef enum { ___3595, ___3593, END_SidebarSizing_e, ___3594 = ___329 } SidebarSizing_e; typedef enum { ___3590, ___3591,    /**@internal TP_NOTAVAILABLE*/ ___3592,      /**@internal TP_NOTAVAILABLE*/ ___3588,   /**@internal TP_NOTAVAILABLE*/ END_SidebarLocation_e, ___3589 = ___329 } SidebarLocation_e; typedef enum { ___2415, ___2418, ___2416, ___2417, END_MenuItem_e, ___2414 = ___329 } MenuItem_e; typedef enum { ___3670, ___3669, ___3679, ___3676, ___3673, ___3668, ___3671, ___3680, ___3678, ___3672, ___3667, ___3675, ___3677, END_StandardMenu_e, ___3674 = ___329 } StandardMenu_e; typedef enum { ___1381, ___1378, ___1382, ___1379, END_FieldProbeDialogPage_e, ___1380 = ___329 } FieldProbeDialogPage_e; enum BooleanCache_e { ___367, ___369, ___370, END_BooleanCache_e, ___368 = ___329 }; enum LinePickLocation_e { ___2276, ___2277, ___2275, ___2274, ___2272, END_LinePickLocation_e, ___2273 = ___329 }; enum ViewDest_e { ___4422, ___4424, END_ViewDest_e, ___4423 = ___329 }; enum DataSetReaderOrigin_e { ___901, ___899, END_DataSetReaderOrigin_e, ___900 = ___329 }; enum ExtendedCurveFitOrigin_e { ___1249, ___1247, END_ExtendedCurveFitOrigin_e, ___1248 = ___329 }; enum CollapsedStatus_e { ___513, ___509, ___508, ___510, ___511, END_CollapsedStatus_e, ___512 = ___329 }; typedef enum { ___4244, ___4249, ___4253, ___4245, ___4254, ___4255, ___4251, ___4250, ___4242, ___4243, ___4252, ___4246, ___4248, END_UndoStateCategory_e, ___4247 = ___329 } UndoStateCategory_e; typedef enum { ___2294, ___2292, END_LinkType_e, ___2293 = ___329 } LinkType_e; typedef enum { ___1509, ___1511, END_FrameCollection_e, ___1510 = ___329 } FrameCollection_e; typedef enum { SurfaceGenerationMethod_AllowQuads, SurfaceGenerationMethod_AllTriangles, SurfaceGenerationMethod_AllPolygons, SurfaceGenerationMethod_Auto, END_SurfaceGenerationMethod_e, SurfaceGenerationMethod_Invalid = ___329 } SurfaceGenerationMethod_e; typedef enum { ___2215, ___2216, ___2217, END_LegendProcess_e, ___2218 = ___329 } LegendProcess_e; typedef enum { ___3382, ___3378, ___3377, ___3381, ___3379, ___3376, END_RGBLegendOrientation_e, ___3380 = ___329 } RGBLegendOrientation_e; struct ___3391 { uint8_t  ___3265; uint8_t  G; uint8_t  B; }; typedef ___3391 ___3388[256]; typedef enum { ___3758, ___3757, ___3766, ___3765, ___3739, ___3711, ___3736, ___3749, ___3706, ___3735, ___3699, ___3719, ___3700, ___3726, ___3725, ___3747, ___3764, ___3756, ___3720, ___3718,
___3760, ___3697, ___3701, ___3748, ___3734, ___3732, ___3741, ___3742, ___3743, ___3744, ___3716, ___3753, ___3750, ___3705, ___3704, StateChange_Text, ___3713, ___3707, ___3710, ___3746, ___3745, ___3690, ___3692, ___3691, ___3759, ___3751, ___3714, ___3755, ___3754, ___3740, ___3737, ___3733, ___3712, ___3738, ___3721, ___3689, ___3688, ___3724, ___3723, ___3722, ___3767, ___3717, ___3702, ___3698, StateChange_OpenLayout, StateChange_MacroLoaded, StateChange_PerformingUndoBegin, StateChange_PerformingUndoEnd, StateChange_SolutiontimeChangeBlockEnd, END_StateChange_e, ___3715 = ___329, ___3703         = ___3716, ___3709          = ___3753, ___3708         = ___3750, ___3762         = ___3719, ___3763               = ___3720, ___3761 = ___3718 } StateChange_e; typedef enum { AnimationType_LineMap, AnimationType_Time, AnimationType_Zone, AnimationType_IJKBlanking, AnimationType_IJKPlanes, AnimationType_IsoSurfaces, AnimationType_Slices, AnimationType_ContourLevels, AnimationType_Streamtraces, END_AnimationType_e, AnimationType_Invalid = ___329 } AnimationType_e; typedef enum { ___3730, ___3731, ___3728, ___3729, END_StateChangeMode_e, ___3727 = ___329 } StateChangeMode_e; typedef enum { ___3695, ___3693, ___3694, END_StateChangeCallbackAPI_e, ___3696 = ___329 } StateChangeCallbackAPI_e; typedef enum { ___86, ___84, ___87, END_AppMode_e, ___85 = ___329 } AppMode_e; typedef enum { ___2209, ___2211, ___2208, END_LayoutPackageObject_e, ___2210 = ___329 } LayoutPackageObject_e; typedef enum { ___4356, ___4357, END_VarLoadMode_e, ___4358 = ___329 } VarLoadMode_e; typedef enum { ___1903, ___1904, END_ImageSelection_e, ___1902 = ___329 } ImageSelection_e; typedef enum { ___2232, ___2235, ___2234, END_LibraryType_e, ___2233 = ___329 } LibraryType_e;   /**@internal TP_NOPYTHON*/ typedef enum { ___220, ___223, ___222, ___224, ___219, ___215, ___216, ___218, ___217, END_AssignOp_e, ___221 = ___329 } AssignOp_e; typedef enum { ___981, ___1000, ___1032, ___1088, ___1050, ___1051, ___1082, ___1047, ___1048, ___1036, ___1037, ___1057, ___1058, ___1029, ___1087, ___1045, ___1019, ___1001, ___1030, ___1031, ___979, ___1068, ___1052, ___1071, ___1073, ___1069, ___1022, ___1066, ___983, ___1084, ___1086, ___1083, ___1085, ___1063, ___1061, ___1062, ___1054, ___1053, ___1027, ___1018, ___996, ___1025, ___978, ___1081, ___1043, ___992, ___1070, ___1067, ___1078, ___1055, ___985, ___987, ___986, ___998, ___1034, ___988, ___989, ___994, ___995, ___1002, ___1004, ___1005,
___1009, ___1008, ___1010, ___1011, ___1003, ___1007, ___1006, ___1026, ___1021, ___1023, ___1080, ___991, ___990, ___993, ___1040, ___1039, ___1056, ___1075, ___1074, ___1079, ___1046, ___982, ___1035, ___1065, ___1060, ___1059, ___1028, ___1016, ___1015, ___1049, ___999, ___984, ___1072, ___1076, ___1041, ___1017, ___1013, ___980, Dialog_FourierTransform, Dialog_RotateData, Dialog_AxialDuplicate, END_Dialog_e, ___1020 = ___329, ___1042 = ___1088 } Dialog_e;   /**@internal TP_NOPYTHON*/ typedef enum { ___49, ___48, ___50, ___46, ___45, ___47, ___42, ___41, ___43, END_AnchorAlignment_e, ___44 = ___329 } AnchorAlignment_e; enum PositionAtAnchor_e { ___3165, ___3166, ___3163, END_PositionAtAnchor_e, ___3164 = ___329 }; struct ___1044 { AnchorAlignment_e  ___40; ___372          ___51; ___372          ___55; int32_t            ___2483; ___2227          ___1993; ___2227          ___2126; PositionAtAnchor_e ___3162; ___372          ___1819; };
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___3222, ___3223, ___3224, ___3225, ___3226, ___3227, ___3228, ___3229, ___3230, ___3231, ___3232, END_ProcessXYMode_e, ___3221 = ___329 } ProcessXYMode_e;
 #endif
typedef enum { ___719, ___722, ___721, END_CurveInfoMode_e, ___720 = ___329 } CurveInfoMode_e; enum ProcessLineMapMode_e { ___3206, ___3215, ___3207, ___3214, ___3204, ___3210, ___3209, ___3212, ___3208, ___3213, ___3205, ___3218, ___3220, ___3216, ___3211, ___3219, END_ProcessLineMapMode_e, ___3217 = ___329 }; typedef enum { ___3864, ___3863, END_StyleBase_e, ___3865 = ___329 } StyleBase_e; typedef enum { ___3282, ___3280, ___3284,
 #if defined ENABLE_ORPHANED_DATASETS
___3283,
 #endif
END_ReadDataOption_e, ___3281 = ___329 } ReadDataOption_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___2719, ___2721, ___2722, END_NodeLabel_e, ___2720 = ___329 } NodeLabel_e;
 #endif
typedef enum { ___2175, ___2177, ___2178, END_LabelType_e, ___2176 = ___329 } LabelType_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___3875, ___3873, ___3877, ___3874, END_SubBoundaryEditOption_e, ___3876 = ___329 } SubBoundaryEditOption_e;
 #endif
typedef enum { ___374, ___373, ___377, ___375, END_BorderAction_e, ___376 = ___329 } BorderAction_e; typedef enum { ___3130, ___3131, ___3132, ___3123, ___3133, ___3134, ___3135, ___3139, ___3140, ___3126, ___3128, ___3129, ___3136, ___3124, ___3137, ___3138, ___3127, END_PointerStyle_e, ___3125 = ___329 } PointerStyle_e; typedef enum { ___711, ___709, ___695, ___696, ___708, ___716, ___703, ___713, ___714, ___700, ___704, ___705, ___707, ___697, ___710, ___712, ___701, ___715, ___702, ___706, ___698, END_CursorStyle_e, ___699 = ___329 } CursorStyle_e; typedef enum { ___3076, ___3085, ___3077, ___3082, ___3084, ___3086, ___3087, ___3079, ___3080, ___3078, ___3083, END_PickSubPosition_e, ___3081 = ___329 } PickSubPosition_e; typedef enum { ___3964, ___3963, ___3962, ___3960, END_TecEngInitReturnCode_e, ___3961 = ___329 } TecEngInitReturnCode_e; typedef enum { ___1792, ___1793, ___1794, ___1789, ___1790, END_GetValueReturnCode_e, ___1791 = ___329, ___1787              = ___1792, ___1788 = ___1793, ___1795     = ___1794, ___1782    = ___1789, ___1784 = ___1790, ___1785         = ___1791 } GetValueReturnCode_e; typedef enum { ___3535, ___3524, ___3530, ___3532, ___3533, ___3534, ___3537, ___3538, ___3521, ___3531, ___3527, ___3522, ___3523, ___3536, ___3525, ___3528, END_SetValueReturnCode_e, ___3526 = ___3524, ___3529 = ___329, ___3519                       = ___3535, ___3509           = ___3524, ___3514     = ___3530, ___3516   = ___3532, ___3517     = ___3533, ___3518  = ___3534, ___3539          = ___3537, ___3540         = ___3538, ___3505            = ___3521, ___3515         = ___3531, ___3512      = ___3527, ___3507            = ___3522, ___3508            = ___3523, ___3520 = ___3536, ___3510      = ___3525, ___3511                  = ___3526, ___3513                  = ___3529 } SetValueReturnCode_e; typedef enum { ___817, ___818, END_DataAlterMode_e, ___819 = ___329 } DataAlterMode_e; typedef enum { ___827, ___825, ___826, ___823, ___824, ___822, ___820, END_DataAlterReturnCode_e,
___821 = ___329 } DataAlterReturnCode_e; typedef enum { ___2854, ___2855, ___2852, ___2856, ___2851, END_ObjectAlign_e, ___2853 = ___329 } ObjectAlign_e; typedef enum { ___2167, ___2166, ___2169, END_LabelAlignment_e, ___2168 = ___329 } LabelAlignment_e; typedef enum { ___4425, ___4421, ___4415, ___4433, ___4419, ___4437, ___4438, ___4428, ___4420, ___4431, ___4432, ___4434, ___4430, ___4417, ___4429, ___4416, ___4418, ___4426, END_View_e, ___4427 = ___329 } View_e; typedef enum { ___4470, ___4468, ___4469, ___4473, ___4472, ___4476, ___4474, ___4475, END_WorkspaceView_e, ___4471 = ___329 } WorkspaceView_e; typedef enum { ___192, ___189, ___190, END_ArrowheadStyle_e, ___191 = ___329, ___186   = ___192, ___183  = ___189, ___184  = ___190, ___185 = ___191 } ArrowheadStyle_e; typedef enum { ___181, ___177, ___179, ___178, END_ArrowheadAttachment_e, ___180 = ___329, ___182        = ___181, ___172 = ___177, ___174       = ___179, ___173  = ___178, ___175     = ___180 } ArrowheadAttachment_e; typedef enum { ___498, ___497, END_Clipping_e, ___499 = ___329 } Clipping_e; typedef enum { ___3771, ___3772, ___3773, ___3776, ___3775, END_StatusInfo_e, ___3774 = ___329 } StatusInfo_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___1515, ___1518, ___1519, ___1520, ___1517, END_FrameMode_e, ___1516 = ___329, ___1512   = ___1515, ___1526  = ___1518, ___1527    = ___1519, ___1528      = ___1520, ___1525  = ___1517, ___1514 = ___1516 } FrameMode_e;
 #endif
typedef enum { ___3113, ___3115, ___3114, ___3121, ___3118, ___3117, END_PlotType_e, ___3116 = ___329 } PlotType_e;
 #define ___4310(___3112) ( VALID_ENUM((___3112), PlotType_e) && \
 ((___3112) != ___3113) )
 #define ___4308(___3112) ( (___3112) == ___3121 || \
 (___3112) == ___3117 )
 #define ___4305(___3112) ( (___3112) == ___3114 || \
 (___3112) == ___3115 )
 #define ___3119(___3112) ___4305((___3112))
 #define ___3120(___3112) ___4308((___3112))
 #define ___4313(___3112) ( (___3112) == ___3118 || \
 (___3112) == ___3121 || \
 (___3112) == ___3114 || \
 (___3112) == ___3115 )
typedef enum { ___555, ___556, END_ContLineCreateMode_e, ___554 = ___329 } ContLineCreateMode_e; typedef enum { ___3055, ___3048, ___3043, ___3066, ___3049, ___3065, ___3047, ___3046, ___3060, ___3053, ___3058, ___3057, ___3063, ___3064, ___3056, ___3071, ___3070, ___3062, ___3061, ___3051, ___3059, ___3054, ___3044, END_PickObjects_e, ___3050 = ___329, ___3038                   = ___3055, ___3032                  = ___3048, ___3029                   = ___3043, ___3028      = ___3066, ___3033                   = ___3049, ___3072                   = ___3065, ___3031          = ___3047, ___3030           = ___3046, ___3045          = ___3060, ___3036             = ___3053, ___3073               = ___3053, ___3041        = ___3058, ___3040 = ___3057, ___3068    = ___3063, ___3069    = ___3064, ___3039                  = ___3056, ___3075                   = ___3071, ___3074              = ___3070, ___3067         = ___3062, ___3052               = ___3061, ___3035          = ___3051, ___3042              = ___3059, ___3037            = ___3054, ___3034                = ___3050 } PickObjects_e; enum SingleEditState_e { ___3598, ___3597, ___3599, END_SingleEditState_e, ___1143 = ___329 }; enum AxisSubPosition_e { ___307, ___306, ___309, ___305, ___310, ___311, END_AxisSubPosition_e, ___308 = ___329, ___304 = ___307, ___303 = ___309, ___313 = ___307, ___312 = ___311 }; typedef enum { ___300, ___299, ___302, END_AxisSubObject_e, ___301 = ___329 } AxisSubObject_e; typedef enum { ___2511, ___2512, ___2510, END_MouseButtonClick_e, ___2509 = ___329 } MouseButtonClick_e; typedef enum { ___2514, ___2535, ___2526, ___2534, ___2525, ___2515, ___2524, ___2530, ___2532, ___2537, ___2528, ___2536, ___2527, ___2516, ___2523, ___2531, ___2533, ___2538, ___2529, ___2520, ___2522, ___2519, ___2521, END_MouseButtonDrag_e, ___2513 = ___329, ___2518  = ___2524, ___2517 = ___2523 } MouseButtonDrag_e; struct ___2508 { MouseButtonClick_e ___421; MouseButtonDrag_e  ___3596; MouseButtonDrag_e  ___646; MouseButtonDrag_e  ___31; MouseButtonDrag_e  ___3557; MouseButtonDrag_e  ___644; MouseButtonDrag_e  ___649; MouseButtonDrag_e  ___35; MouseButtonDrag_e  ___645; }; struct ___2506 { ___2508 ___2468; ___2508 ___3394; };
 #define ___2751   0
 #define ___2214   1
 #define ___2469 2
 #define ___3396  3
typedef enum { ___33, ___34, END_AltMouseButtonMode_e, ___32 = ___329 } AltMouseButtonMode_e; typedef enum { ___2557, ___2568, ___2539, ___2578, ___2573, ___2558, ___2572, ___2551, ___2554, ___2549, ___2552, ___2550, ___2553, ___2545, ___2563, ___2559, ___2564, ___2565, ___2566, ___2567, ___2543, ___2541, ___2542, ___2571, ___2570, ___2548, ___2547, ___2546, ___2544, ___2569, ___2556, ___2574, ___2575, ___2576, ___2577, ___2561, ___2562, ___2540, END_MouseButtonMode_e, ___2555 = ___329, ___2560 = ___2563, ___2595                = ___2557, ___2603                = ___2568, ___2507                = ___2539, ___2613                  = ___2578, ___2608             = ___2573, ___2596                 = ___2558, ___2607                  = ___2572, ___2590          = ___2551, ___2593            = ___2554, ___2588            = ___2549, ___2591         = ___2552, ___2589           = ___2550, ___2592            = ___2553, ___2584           = ___2545, ___2598       = ___2563, ___2597      = ___2559, ___2599           = ___2564, ___2600           = ___2565, ___2601           = ___2566, ___2602           = ___2567, ___2582          = ___2543, ___2580            = ___2541, ___2581         = ___2542, ___2606          = ___2571, ___2605         = ___2570, ___2587         = ___2548, ___2586           = ___2547, ___2585 = ___2546, ___2583    = ___2544, ___2604                 = ___2569, ___2609                 = ___2574, ___2610                 = ___2575, ___2611                 = ___2576, ___2612                 = ___2577, ___2594               = ___2555 } MouseButtonMode_e; typedef enum { ___976, ___975, ___977, END_DetailsButtonState_e, ___974 = ___329 } DetailsButtonState_e; typedef enum { ___1193, ___1194, ___1192, ___1198, ___1195, ___1197, END_Event_e, ___1196 = ___329 } Event_e; typedef enum { ___2857, ___2859, ___2861, ___2860, END_ObjectDrawMode_e, ___2858 = ___329 } ObjectDrawMode_e; typedef enum { ___4163, ___4165, END_ThreeDViewChangeDrawLevel_e, ___4164 = ___329 } ThreeDViewChangeDrawLevel_e; typedef enum { ___2746, ___2748, END_NonCurrentFrameRedrawLevel_e, ___2747 = ___329 } NonCurrentFrameRedrawLevel_e; typedef enum { ___3312, ___3315, ___3313, ___3316, ___3305, ___3306, ___3307, ___3303, ___3304, ___3311, ___3310, ___3309, END_RedrawReason_e, ___3308 = ___329, ___3314 = ___3312, ___3317 = ___3315
} RedrawReason_e; typedef enum { ___3418, ___3417, ___3416, END_RotationMode_e, ___3415 = ___329 } RotationMode_e; typedef enum { ___3409, ___3410, ___3411, ___3405, ___3406, ___3401, ___3407, ___3408, ___3403, ___3400, ___3402,   /**@internal TP_NOTAVAILABLE*/ END_RotateAxis_e, ___3404 = ___329 } RotateAxis_e; typedef enum { ___3412, ___3414, END_RotateOriginLocation_e, ___3413 = ___329 } RotateOriginLocation_e; typedef enum { ___2885, ___2887, END_OriginResetLocation_e, ___2886 = ___329 } OriginResetLocation_e; typedef enum { ___3614, ___3615, ___3613, ___3612, END_SliceSource_e, ___3611 = ___329 } SliceSource_e; typedef enum { ___1950, ___1949, ___1946, ___1945, ___1942, ___1948, ___1953, ___1943, END_Input_e, ___1947 = ___329 } Input_e; typedef enum { ___3252, ___3254, ___3255, END_PtSelection_e, ___3253 = ___329 } PtSelection_e; typedef enum { ___1119, ___1118, ___1120, END_Drift_e, ___1117 = ___329 } Drift_e; typedef enum { ___966, ___967, ___971, ___970, ___969, END_DerivPos_e, ___968 = ___329 } DerivPos_e; typedef enum { ___2242, ___2244, END_LinearInterpMode_e, ___2243 = ___329 } LinearInterpMode_e; typedef enum { ___4442, ___4443, END_VolumeCellInterpolationMode_e, ___4441 = ___329 } VolumeCellInterpolationMode_e; typedef enum { ___3156, ___3154, END_PolyCellInterpolationMode_e, ___3155 = ___329 } PolyCellInterpolationMode_e; typedef enum { ___549, ___548, END_ConstraintOp2Mode_e, ___547 = ___329 } ConstraintOp2Mode_e; typedef enum { ___877, ___879, END_DataProbeVarLoadMode_e, ___878 = ___329 } DataProbeVarLoadMode_e; typedef enum { ___3936, ___3934, END_SZLSubzoneLoadModeForStreams_e, ___3935 = ___329 } SZLSubzoneLoadModeForStreams_e; typedef enum { ___4316, ___4317, ___4320, END_ValueBlankCellMode_e, ___4318 = ___329, ___4319 = ___4320 } ValueBlankCellMode_e; typedef enum { ___4321, ___4324, ___4322, END_ValueBlankMode_e, ___4323 = ___329 } ValueBlankMode_e; typedef enum { ___454, ___455, ___452, ___456, END_CellBlankedCond_e, ___453 = ___329 } CellBlankedCond_e; typedef enum { ___3331, ___3328, ___3330, ___3327, ___3326, ___3332, END_RelOp_e, ___3329 = ___329 } RelOp_e; typedef enum { ___1846, ___1845, END_IJKBlankMode_e, ___1847 = ___329 } IJKBlankMode_e; typedef enum { ___3109, ___3111, ___3108, END_PlotApproximationMode_e, ___3110 = ___329 } PlotApproximationMode_e; typedef enum { ___3652, ___3653, ___3650, END_SphereScatterRenderQuality_e, ___3651 = ___329 } SphereScatterRenderQuality_e; typedef enum { ExtractMode_SingleZone, ExtractMode_OneZonePerConnectedRegion, ExtractMode_OneZonePerSourceZone, END_ExtractMode_e, ExtractMode_Invalid = ___329 } ExtractMode_e; typedef enum { Resulting1DZoneType_IOrderedIfPossible, Resulting1DZoneType_FELineSegment, Resulting1DZoneType_Unused, END_Resulting1DZoneType_e, ResultingDZoneType_Invalid = ___329 } Resulting1DZoneType_e; typedef enum { ___2990, ___2988, ___2989, ___2984, END_FillPat_e, ___2985 = ___329 } FillPat_e; typedef enum { ___4224, ___4222, ___4223, ___4220, END_Translucency_e, ___4221 = ___329
} Translucency_e; typedef enum { ___3889, ___3890, ___3887, END_SunRaster_e, ___3888 = ___329 } SunRaster_e; typedef enum { ___384, ___387, ___386, END_BoundaryCondition_e, ___385 = ___329 } BoundaryCondition_e; typedef enum { ___289, ___292, ___291, END_AxisMode_e, ___290 = ___329 } AxisMode_e; typedef enum { ___3263, ___3261, ___3264, END_QuickColorMode_e, ___3262 = ___329 } QuickColorMode_e; typedef enum { ___1414, ___1417, ___1416, ___1415, END_FillMode_e, ___1413 = ___329 } FillMode_e; typedef enum { ___2271, ___2267, ___2265, ___2268, ___2270, ___2266, END_LinePattern_e, ___2269 = ___329 } LinePattern_e; typedef enum { ___2129, ___2130, ___2127, END_LineJoin_e, ___2128 = ___329 } LineJoin_e; typedef enum { ___442, ___445, ___446, END_LineCap_e, ___443 = ___329 } LineCap_e; typedef enum { ___1587, ___1589, ___1590, ___1583, ___1584, ___1588, ___1585, END_GeomForm_e, ___1586 = ___329, GeomType_LineSegs = ___1587, GeomType_Rectangle = ___1589, GeomType_Square = ___1590, GeomType_Circle = ___1583, GeomType_Ellipse = ___1584, GeomType_LineSegs3D = ___1588, GeomType_Image = ___1585, END_GeomType_e = END_GeomForm_e, GeomType_Invalid = ___1586 } GeomForm_e; typedef GeomForm_e GeomType_e; typedef enum { ___4348, ___4347, END_VariableDerivationMethod_e, ___4349 = ___329 } VariableDerivationMethod_e; typedef enum { ___270, END_AuxDataType_e, ___269 = ___329 } AuxDataType_e; typedef enum { ___259, ___253, ___254, ___258, ___256, ___257, AuxDataLocation_Layout, END_AuxDataLocation_e, ___255 = ___329 } AuxDataLocation_e; typedef enum { ___4704, ___4702, ___4700, ___4701, ___4695, ___4696, ___4698, ___4699, ___4697, END_ZoneType_e, ___4703 = ___329 } ZoneType_e; typedef enum { ___4659, ___4664, ___4666, ___4660, ___4662, ___4665, ___4661, END_ZoneOrder_e, ___4663 = ___329 } ZoneOrder_e; typedef enum { ___853, ___854, ___851, ___852, END_DataFormat_e, ___855 = ___329 } DataFormat_e; typedef enum { ___874, ___876, END_DataPacking_e, ___875 = ___329 } DataPacking_e; typedef enum { ProbeObject_None, ProbeObject_Streamtrace, ProbeObject_StreamtraceMarker, ProbeObject_Slice, ProbeObject_IsoSurface, ProbeObject_FieldZone, END_ProbeObject_e, ProbeObject_First = ProbeObject_Streamtrace, ProbeObject_Last  = ProbeObject_FieldZone, ProbeObject_Invalid = ___329 } ProbeObject_e; typedef enum { ProbeNearest_Position, ProbeNearest_Node, END_ProbeNearest_e, ProbeNearest_Invalid = ___329 } ProbeNearest_e; typedef enum { ___2992, ___2993, ___2996, ___2995, ___2991, ___2997, ___2998, ___2999, END_PrinterDriver_e, ___2994 = ___329 } PrinterDriver_e; typedef enum { ___1888, ___1905, ___1883, ___1885, END_EPSPreviewImage_e, ___1887 = ___329 } EPSPreviewImage_e; typedef enum { ___4172, ___4174, END_TIFFByteOrder_e, ___4173 = ___329 } TIFFByteOrder_e; typedef enum { ___2133, ___2132, END_JPEGEncoding_e, ___2131 = ___329 } JPEGEncoding_e; typedef enum { ___1433, ___1432, ___1430, END_FlashImageType_e, ___1431 = ___329, ___1429 = ___1430 } FlashImageType_e; typedef enum { ___1426, ___1428, END_FlashCompressionType_e, ___1427 = ___329 } FlashCompressionType_e; typedef enum { ___1225, ___1229, ___1226, ___1227, ___1233, ___1224, ___1216, ___1217, ___1223, ___1213, ___1220, ___1230, ___1211, ___1222, ___1210, ___1212,    /**@internal TP_NOTAVAILABLE*/ ___1219, ___1214, ___1232, ___1228, ___1215, ___1221, ___1231,
END_ExportFormat_e, ___1218 = ___329 } ExportFormat_e; typedef enum { ___275, ___277, ___278, END_AVICompression_e, ___276 = ___329 } AVICompression_e; typedef enum { ___65, ___59, ___64, ___60, ___61, ___63, ___66, END_AnimationDest_e, ___62 = ___329 } AnimationDest_e; typedef enum { ___69, ___67, ___71, ___68, END_AnimationOperation_e, ___70 = ___329 } AnimationOperation_e; typedef enum { ___73, ___78, ___72, ___79, ___75, ___77, ___76, END_AnimationStep_e, ___74 = ___329 } AnimationStep_e; typedef enum { ___4605, ___4603, ___4606, END_ZoneAnimationMode_e, ___4604 = ___329 } ZoneAnimationMode_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___361, ___360, ___363, END_BitDumpRegion_e, ___362 = ___329 } BitDumpRegion_e;
 #endif
typedef enum { ___1239, ___1238, ___1241, END_ExportRegion_e, ___1240 = ___329 } ExportRegion_e; typedef enum { ___2956, ___2954, ___2951, ___2950, ___2952, ___2953, END_PaperSize_e, ___2955 = ___329, ___2949  = ___2956, ___2947  = ___2954, ___2944      = ___2951, ___2943      = ___2950, ___2945 = ___2952, ___2946 = ___2953, ___2948 = ___2955 } PaperSize_e; typedef enum { ___2958, ___2961, ___2970, ___2965, ___2959, ___2962, ___2967, ___2969, ___2968, ___2957, ___2966, ___2964, ___2963, END_PaperUnitSpacing_e, ___2960 = ___329 } PaperUnitSpacing_e; typedef enum { ___2920, ___2921, ___2918, END_Palette_e, ___2919 = ___329 } Palette_e; typedef enum { ___3191, ___3189, END_PrintRenderType_e, ___3190 = ___329 } PrintRenderType_e; typedef enum { ___4269, ___4268, ___4271, ___4272, ___4267, END_Units_e, ___4270 = ___329 } Units_e; typedef enum { ___659, ___660, END_CoordScale_e, ___658 = ___329, ___3436 = ___659, ___3437 = ___660, ___3435 = ___658 } CoordScale_e;
 #define ___1753(___3265) ( ((___3265) < ___3628) ? ___3631 : ( ((___3265) > ___2179) ? ___2190 : ___2316((___3265)) ) )
typedef enum { CoordSys_Grid, CoordSys_Frame, ___661, ___664, ___665, ___662, CoordSys_Grid3D, ___666, END_CoordSys_e, ___663 = ___329 } CoordSys_e; typedef enum { ___3444, ___3446, END_Scope_e, ___3445 = ___329 } Scope_e; typedef enum { ___4050, ___4045, ___4055, ___4052, ___4051, ___4053, ___4047, ___4046, ___4048, ___4054, END_TextAnchor_e, ___4049 = ___329 } TextAnchor_e; typedef enum { TextType_Regular, TextType_LaTeX, END_TextType_e, TextType_Invalid = ___329 } TextType_e; typedef enum { ___4075, ___4063, ___4069, END_TextBox_e, ___4070 = ___329 } TextBox_e; typedef enum { ___1647, ___1637, ___1639, ___1645, ___1641, ___1638, ___1635, ___1636, ___1646, ___1642, ___1644, ___1643, GeomShape_LineArt, END_GeomShape_e, ___1640 = ___329 } GeomShape_e; typedef enum { ___348, ___347, ___346, ___345, ___343, END_BasicSize_e, ___344 = ___329 } BasicSize_e; typedef enum { ___2249, ___2246, ___2247, ___2251, ___2252, ___2250, END_LineForm_e, ___2248 = ___329 } LineForm_e; typedef enum { ___729, ___731, ___726, ___732, ___733, ___730, ___727, END_CurveType_e, ___728 = ___329, ___725 = ___731 } CurveType_e; typedef enum { ___3455, ___3457, ___3456, END_Script_e, ___3454 = ___329 } Script_e; typedef enum { ___1454, ___1455, ___1449, ___1459, ___1470, ___1466, ___1468, ___1467, ___1469, ___1445, ___1446, ___1448, Font_HelveticaItalic, Font_HelveticaItalicBold, Font_CourierItalic, Font_CourierItalicBold, END_Font_e, ___1456 = ___329 } Font_e; typedef enum { ___1465, ___1464, ___1461, ___1462, END_FontStyle_e, ___1463 = ___329 } FontStyle_e; typedef enum { ___4233, ___4232, END_TwoDDrawOrder_e, ___4234 = ___329 } TwoDDrawOrder_e; typedef enum { ___1114, ___1115, END_DrawOrder_e, ___1116 = ___329 } DrawOrder_e; typedef enum { ___3806, ___3807, ___3809, ___3810, ___3811, ___3808, END_Streamtrace_e, ___3804 = ___329 } Streamtrace_e; typedef enum { ___3794, ___3796, ___3793, END_StreamDir_e, ___3795 = ___329 } StreamDir_e; typedef enum { ___1092, ___1093, ___1094, ___1095, DistributionRegion_SurfacesOfSuppliedZones, END_DistributionRegion_e, ___1091 = ___329 } DistributionRegion_e; typedef enum { ___2045, ___2047, ___2049, ___2048, END_IsoSurfaceSelection_e, ___2046 = ___329 } IsoSurfaceSelection_e; typedef enum { ___4328, ___4330, END_ValueLocation_e, ___4329 = ___329 } ValueLocation_e; typedef enum { OffsetDataType_32Bit, OffsetDataType_64Bit, END_OffsetDataType_e, OffsetDataType_Invalid = ___329 } OffsetDataType_e;
 #define VALID_32OR64BIT_OFFSET_TYPE(t) ((t) == OffsetDataType_32Bit || (t) == OffsetDataType_64Bit)
typedef enum { ___1371,   /**@internal TP_NOTAVAILABLE*/ FieldDataType_Float, FieldDataType_Double, FieldDataType_Int32, FieldDataType_Int16, FieldDataType_Byte, ___1365, END_FieldDataType_e, ___1366,     /**@internal TP_NOTAVAILABLE*/ ___1368,   /**@internal TP_NOTAVAILABLE*/
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
___1370 = FieldDataType_Int32, ___1373 = FieldDataType_Int16,
 #endif
___1369 = ___329 } FieldDataType_e;
 #define VALID_FIELD_DATA_TYPE(___1364) (VALID_ENUM((___1364),FieldDataType_e) && \
 (___1364)!=___1371)
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___2431, ___2426, ___2424, END_MeshPlotType_e, ___2425 = ___329 } MeshPlotType_e;
 #endif
typedef enum { ___2430, ___2429, ___2427, END_MeshType_e, ___2428 = ___329 } MeshType_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___636, ___616, ___637, ___557, ___615, END_ContourPlotType_e, ___617 = ___329 } ContourPlotType_e;
 #endif
typedef enum { ___641, ___639, ___642, ___638, ___643, END_ContourType_e, ___640 = ___329 } ContourType_e; typedef enum { ___567, ___558, ___559, ___560, ___561, ___562, ___563, ___564, ___565, END_ContourColoring_e, ___566 = ___329 } ContourColoring_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___4400, ___4396, ___4399, ___4397, END_VectorPlotType_e, ___4398 = ___329 } VectorPlotType_e;
 #endif
typedef enum { ___4405, ___4401, ___4404, ___4402, END_VectorType_e, ___4403 = ___329 } VectorType_e; typedef enum { ___3547, ___3546, ___3544, ___3543, ___3542, END_ShadePlotType_e, ___3545 = ___329 } ShadePlotType_e; typedef enum { ___2239, ___2236, ___2238, END_LightingEffect_e, ___2237 = ___329 } LightingEffect_e; enum IJKLines_e { ___1857, ___1859, ___1860, END_IJKLines_e, ___1858 = ___329, ___2283       = ___1857, ___2285       = ___1859, ___2286       = ___1860, ___2284 = ___1858 }; typedef enum { ___1850, ___1848, ___1851, END_IJKCellType_e, ___1849 = ___329 } IJKCellType_e; typedef enum { ___1867, ___1872, ___1874, ___1866, ___1868, ___1873, ___1870, ___1869, ___1876, ___1875, END_IJKPlanes_e, ___1871 = ___329, ___3098       = ___1867, ___3103       = ___1872, ___3105       = ___1874, ___3099      = ___1868, ___3104      = ___1873, ___3101      = ___1870, ___3100     = ___1869, ___3097    = ___1866, ___3107  = ___1876, ___3106  = ___1875, ___3102 = ___1871 } IJKPlanes_e; typedef enum { ___3898, ___3899, ___3904, ___3906, ___3907, ___3901, ___3905, ___3902, ___3900, ___3897, ___3908, END_SurfacesToPlot_e, ___3903 = ___329 } SurfacesToPlot_e; typedef enum { ___3150, ___3147, ___3149, ___3145, ___3146, END_PointsToPlot_e,
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
___3151 = ___3150, ___3144          = ___3147,
 #endif
___3148 = ___329 } PointsToPlot_e; typedef enum { ___3622, ___3624, ___3626, ___3619, ___3620, ___3621, ___3617, ___3616, END_SliceSurface_e, ___3618 = ___329 } SliceSurface_e;   /* pytecplot defines this separately with CVar removed. */ typedef enum { ___503, ___501, ___500, END_ClipPlane_e, ___502 = ___329 } ClipPlane_e; typedef enum { ___3606, ___3605, END_SkipMode_e, ___3607 = ___329 } SkipMode_e; typedef enum { ___1139, ___1141, ___1140, END_EdgeType_e, ___1142 = ___329 } EdgeType_e;
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef enum { ___391, ___390, ___389, ___383, END_BoundPlotType_e, ___388 = ___329 } BoundPlotType_e;
 #endif
typedef enum { ___397, ___396, ___395, ___393, END_BoundaryType_e, ___394 = ___329 } BoundaryType_e; typedef enum { ___382, ___381, ___380, ___378, END_BorderLocation_e, ___379 = ___329 } BorderLocation_e; typedef enum { ___610, ___581, ___584, ___578, ___613, ___612, ___611, ___592, ___569, ___570, ___568, ___571, ___572, ___573, ___575, ___576, ___577, ___579, ___582, ___574, ___583, ___585, ___586, ___587, ___588, ___589, ___590, ___591, ___593, ___594, ___595, ___597, ___596, ___598, ___599, ___600, ___604, ___601, ___602, ___603, ___605, ___606, ___607, ___608, ___609, END_ContourColorMap_e, ___580 = ___329, ___530    = ___610, ___527    = ___581, ___528       = ___584, ___525    = ___578, ___533         = ___613, ___532      = ___612, ___531     = ___611, ___529   = ___592, ___526      = ___580 } ContourColorMap_e; typedef enum { ___1184, ___1179, ___1182, ___1183, ___1180, ___1185, ___1178, END_ErrorBar_e, ___1181 = ___329 } ErrorBar_e; typedef enum { ___635, ___634, ___632, END_ContourLineMode_e, ___633 = ___329 } ContourLineMode_e; enum Panel_e { ___2922, ___2923, ___2942, ___2933, ___2928, ___2924, ___2926, ___2932, ___2925, ___2939, ___2940, ___2938, ___2930, Panel_IsoSurfaceVector, ___2931, ___2934, ___2937, ___2936, ___2935, ___2927, ___2941, END_Panel_e, ___2929 = ___329 }; typedef enum { ___2443, ___2447, ___2444, ___2446,     /**@internal TP_NOPYTHON*/ ___2449, ___2450, ___2448, END_MessageBoxType_e, ___2445 = ___329, ___2433           = ___2443, ___2451         = ___2447, ___2434     = ___2444, ___2436        = ___2446, ___2453           = ___2449, ___2454     = ___2450, ___2452 = ___2448, ___2435         = ___2445 } MessageBoxType_e; typedef enum { ___2441, ___2439, ___2437, ___2440, END_MessageBoxReply_e, ___2438 = ___329 } MessageBoxReply_e; typedef enum { ___2772, ___2771, ___2770, ___2767, ___2776, ___2768, ___2774, ___2775, ___2769, ___2777, END_NumberFormat_e, ___2773 = ___329 } NumberFormat_e; typedef NumberFormat_e ValueFormat_e; typedef enum { ___324, ___325, ___323, END_BackingStoreMode_e, ___322 = ___329 } BackingStoreMode_e; typedef enum { ___4168, ___4170, ___4167, END_TickDirection_e, ___4169 = ___329 } TickDirection_e; typedef enum { ___320, ___318, ___321, END_AxisTitlePosition_e, ___319 = ___329 } AxisTitlePosition_e; typedef enum { ___315, ___317, ___316, END_AxisTitleMode_e, ___314 = ___329 } AxisTitleMode_e; typedef enum { ___288, ___286, ___285,
___284, AxisAlignment_WithSpecificAngle, ___283, ___280, ___281, ___282, END_AxisAlignment_e, ___279 = ___329 } AxisAlignment_e; typedef enum { ___1539, ___1540, END_FunctionDependency_e, ___1536 = ___329, ___1538 = ___1539, ___1537 = ___1540 } FunctionDependency_e; typedef enum { LegendShow_Always, LegendShow_Never, ___2219, END_LegendShow_e, ___2220 = ___329, ___2222 = LegendShow_Always, ___2221  = LegendShow_Never } LegendShow_e; typedef enum { ___2259, LineMapSort_ByIndependentVar, LineMapSort_ByDependentVar, LineMapSort_BySpecificVar, END_LineMapSort_e, ___2258 = ___329, ___2257 = LineMapSort_ByIndependentVar, ___2256   = LineMapSort_ByDependentVar, ___2260    = LineMapSort_BySpecificVar, } LineMapSort_e; typedef enum { ___551, ___552, ___550, END_ContLegendLabelLocation_e, ___553 = ___329 } ContLegendLabelLocation_e; typedef enum { ___4144, ___4146, ___4143, END_ThetaMode_e, ___4145 = ___329 } ThetaMode_e; typedef enum { ___4213, ___4216, ___4214, ___4215, END_Transform_e, ___4212 = ___329 } Transform_e; typedef enum { ___4464, ___4465, ___4462, ___4461, END_WindowFunction_e, ___4463 = ___329 } WindowFunction_e; typedef enum { ___2205, ___2206, ___2204, END_LaunchDialogMode_e, ___2203 = ___329 } LaunchDialogMode_e; typedef enum { ___3468, ___3467, ___3465, ___3471, ___3470, END_SelectFileOption_e, ___3466 = ___329 } SelectFileOption_e; typedef enum { ___356, ___357, ___358, ___354 = ___358, END_BinaryFileVersion_e, ___355 = ___329 } BinaryFileVersion_e; typedef enum { ___4414, ___4412, ___4411, END_ViewActionDrawMode_e, ___4413 = ___329 } ViewActionDrawMode_e; typedef enum { ___2898, ___2899, ___2897, ___2904, ___2905, ___2901, ___2903, ___2902, END_PageAction_e, ___2900 = ___329 } PageAction_e; typedef enum { ___1507, ___1503, ___1501, ___1490, ___1492, ___1505, ___1502, ___1506, ___1489, ___1487, ___1488, ___1484, ___1485, ___1486, ___1497, ___1498, ___1499, ___1494, ___1495, ___1496, ___1508, FrameAction_DeleteByNumber, FrameAction_Reset, END_FrameAction_e, ___1493 = ___329, ___1500 = ___1503, ___1504 = ___1506, ___1491 = ___1490 } FrameAction_e; typedef enum { ___1108, ___1107, ___1109, END_DoubleBufferAction_e, ___1106 = ___329 } DoubleBufferAction_e; typedef enum { ___3010, ___3004, ___3005, ___3015, ___3013, ___3012, ___3011, ___3018, ___3019, ___3023, ___3017, ___3021, ___3020, ___3022, ___3014, ___3009, ___3008, ___3007, ___3006, END_PickAction_e, ___3016 = ___329 } PickAction_e; typedef enum { ___629, ___628, ___631, END_ContourLevelsInitializationMode_e, ___630 = ___329 } ContourLevelsInitializationMode_e; typedef enum { ___621, ___625, ___623, ___626,
___627, ___622, END_ContourLevelAction_e, ___624 = ___329 } ContourLevelAction_e; typedef enum { ___618, ___619, END_ContourLabelAction_e, ___620 = ___329 } ContourLabelAction_e; typedef enum { ___3798, ___3799, ___3800, ___3803, ___3802, END_StreamtraceAction_e, ___3801 = ___329 } StreamtraceAction_e; typedef enum { ___519, ___517, ___520, END_ColorMapControlAction_e, ___518 = ___329 } ColorMapControlAction_e; typedef enum { ___522, ___521, END_ColorMapDistribution_e, ___523 = ___329 } ColorMapDistribution_e; typedef enum { ___3387, ___3386, ___3385, ___3384, END_RGBMode_e, ___3383 = ___329 } RGBMode_e; typedef enum { ___4037, ___4038, END_TecUtilErr_e, ___4036 = ___329 } TecUtilErr_e; enum AxisShape_e { ___298, ___296, ___297, ___294, END_AxisShape_e, ___295 = ___329 }; enum RunMode_e { ___3426, ___3427, ___3428, END_RunMode_e, ___3429 = ___329 }; typedef enum { ___1208, ___1204, ___1209, ___1203, ___1202, ___1207, ___1206, END_ExportCustReturnCode_e, ___1205 = ___329 } ExportCustReturnCode_e; typedef enum { ___808, ___807, ___810, ___811, ___813, ___814, ___812, END_CZType_e, ___809 = ___329 } CZType_e; typedef enum { ___1290, ___1289, ___1287, ___1286, END_FaceNeighborMode_e, ___1288 = ___329 } FaceNeighborMode_e; typedef enum { ___2914, ___2916, ___2915, ___2912, END_PageRenderDest_e, ___2913 = ___329 } PageRenderDest_e; enum RenderDest_e { ___3348, ___3341, ___3340, ___3339, ___3338, ___3346, ___3343, END_RenderDest_e, ___3344 = ___329, ___3342 = ___3348, ___3345 = ___3338 }; typedef enum { ___3781, ___3782, ___3784, END_Stipple_e, ___3783 = ___329 } Stipple_e; typedef enum { ___845, ___846, ___848, END_DataFileType_e, ___847 = ___329 } DataFileType_e; typedef enum { ___540, ___541, END_ConditionAwakeReason_e, ___539 = ___329 } ConditionAwakeReason_e; typedef enum { ___3200, ___3201, ___3198, END_ProbeStatus_e, ___3199 = ___329 } ProbeStatus_e; typedef enum { ___1523, ___1524, END_FrameSizePosUnits_e, ___1522 = ___329 } FrameSizePosUnits_e; typedef enum { ___1814, ___1816, ___1815, END_Gridline_e, ___1813 = ___329 } Gridline_e; typedef enum { ___3170, ___3168, END_PositionMarkerBy_e, ___3169 = ___329 } PositionMarkerBy_e; typedef enum { ___2298, ___2299, ___2300, END_LoaderCallbackVersion_e, ___2297 = ___329 } LoaderCallbackVersion_e; typedef enum { LoaderAdvancedOptions_NotAvailable, LoaderAdvancedOptions_Allow, LoaderAdvancedOptions_ForceLaunch, END_LoaderAdvancedOptions_e, LoaderAdvancedOptions_Invalid = ___329 } LoaderAdvancedOptions_e; typedef enum { ___3025, ___3027, ___3024, PickCollectMode_HomogeneousAdd, END_PickCollectMode_e, ___3026 = ___329 } PickCollectMode_e; typedef enum { ___398, ___400, END_BoundingBoxMode_e, ___399 = ___329 } BoundingBoxMode_e; typedef enum { ImageRenderingStrategy_Auto, ImageRenderingStrategy_OpenGL, ImageRenderingStrategy_Mesa, END_ImageRenderingStrategy_e, ImageRenderingStrategy_Invalid = ___329 } ImageRenderingStrategy_e; typedef enum { PreTranslateData_Auto,
PreTranslateData_On, PreTranslateData_Off, END_PreTranslateData_e, PreTranslateData_Invalid = ___329 } PreTranslateData_e; typedef struct ___2665* ___2664; typedef struct ___3655* ___3654; typedef void*(STDCALL *___4151)(___90 ___4150); typedef struct ___543* ___542; typedef struct ___2122* ___2120; typedef void (STDCALL *___4160)(___90 ___2124); typedef struct StringList_s* ___3839; typedef struct Menu_s*       ___2419; typedef struct ___135*  ___134; typedef struct LineSegmentProbeResult_s* ___2282; typedef enum { ___1900, ___1892, ___1897, ___1898, ___1901, ___1891, ___1893, ___1894, ___1899, ___1895, END_ImageResizeFilter_e, ___1896 = ___329 } ImageResizeFilter_e; typedef enum { ___4379, ___4374, ___4377, ___4375, ___4378, END_VarStatus_e, ___4376 = ___329 } VarStatus_e; typedef enum { ElementOrientation_Standard, ElementOrientation_Reversed, ElementOrientation_Arbitrary, END_ElementOrientation_e, ElementOrientation_Invalid = ___329 } ElementOrientation_e; typedef struct ___3502* ___3501; struct ___4579 { double X; double Y; }; typedef struct { double X; double Y; double Z; } ___4582; typedef ___4582 ___4581[8]; struct ___3250 { double ___3248; double ___4141; double ___30; }; namespace tecplot { class ___4237; } struct ___1550 { double ___4292; double ___4294; double ___4296; }; struct ___4149 { double ___4141; double ___3265; }; union ___54 { ___1550 ___1548; ___4582         ___4580; ___4149      ___4147; }; struct ___841 { char*                  ___3183; char*                  ___4040; DataFileType_e         ___1408; ___1398           ___842; ___3839          ___4363; ___1172             ___2847; ___1172             NumVars; double                 ___3639; struct ___841* ___2703; }; struct ___3872 { ___372 ___1919; ___372 ___1914; ___372 ___1922; ___372 ___1921; ___372 ___1916; ___372 ___1917; ___372 ___1911; ___372 ___1920; ___372 ___1912; ___372 UpdateInvalidContourLevels; ___372 ___1913; ___372 ___537; ___372 ___2423; ___372 ___1915; ___372 ___4284; }; struct ___2044 { ___372 ___3563; ___372 ___3577; ___372 ___3566; ___372 ___3586; ___372 ___3583; ___372 ___4283; ___372 ___4287; }; struct ___3610 { ___372 ___3563; ___372 ___3577; ___372 ___3566; ___372 ___3586; ___372 ___3583; ___372 ___3569; ___372 ___4283; ___372 ___4287; }; struct ___3805 { ___372 ___3563; ___372 ___3578; ___372 ___3567; ___372 ___3564; ___372 ___3577; ___372 ___3566; ___372 ___3583; ___372 ___3576; ___372 ___4283; ___372 ___4287; }; struct ___1374 {
 #if 0 
___372       ___3563;
 #endif
TwoDDrawOrder_e ___4231; ___372       ___3577; ___372       ___3566; ___372       ___3586; ___372       ___3582; ___372       ___3583; ___372       ___3569; ___372       ___4283; ___372       ___4287; }; struct ___803 { ___372 ___3577; ___372 ___3566; ___372 ___3586; ___372 ___3582; ___372 ___3583; ___372 ___3569; ___372 ___4283; ___372 ___4287; }; struct ___2279 {
 #if 0 
___372       ___3563;
 #endif
___372 ___3575; ___372 ___3585; ___372 ___3565; ___372 ___3573; }; union ___1980 { double    ___3434; ___2227 ___3555; }; typedef ___372 (*___3885)(TP_IN_OUT double* ___4315, const char*       ___3883); struct ___1951 { Input_e           ___4236; double            ___2470; double            ___2346; ___1980 ___1978; ___3885 ___3884; }; struct ___3390 { ___516 ___3265; ___516 G; ___516 B; bool operator==(___3390 const& ___3392) const { return ___3265 == ___3392.___3265 && G == ___3392.G && B == ___3392.B; } }; struct ___648 { double ___524; ___3390  ___2212; ___3390  ___4211; bool operator==(___648 const& ___3392) const { return ___524 == ___3392.___524 && ___2212          == ___3392.___2212          && ___4211         == ___3392.___4211; } }; struct ___1191 { int       ___1832; int       ___2106; int       ___2197; int       ___2201; int       ___337; int       ___338; int       ___422; Event_e   ___1190; ___372 ___2057; ___372 ___1999; ___372 ___2011; ___372 ___4450; ___372 ___4448; ___372 ___4449; }; struct ___2329 { ___2325          ___2332; struct ___2329* ___2702; }; struct ___1976 { ___2227 ___4567; ___2227 ___4584; ___2227 ___4568; ___2227 ___4585; }; struct ___3299 { double ___4567; double ___4584; double ___4568; double ___4585; }; struct ___1879 { ___2227  ___1832; ___2227  ___2106; ___2227  ___2135; }; struct ___1174 { ___1172 ___2470; ___1172 ___2346; ___1172 ___3604; }; struct ___1929 { ___2227 ___2470; ___2227 ___2346; ___2227 ___3604; }; struct ___4123 { Font_e                    ___1444; double                    ___1827; Units_e                   ___3601; };
 #define ___202(S)       (((S)->___4282 == ___1305) && (___1457::___1958().___1443((S)->___4238) == ___1449))
 #define ___203(S)        (((S)->___4282 == ___1305) && (___1457::___1958().___1443((S)->___4238) == ___1459))
 #define ___204(S) (((S)->___4282 == ___1305) && (___1457::___1958().___1443((S)->___4238) == ___1470))
struct ___205 { ___372                ___4282; tecplot::___4237 const* ___4238; ___3917             ___472; }; struct ___3919 { GeomShape_e  ___1634; ___372    ___2003; ___205 ___201; };
 #ifdef NOT_USED
struct AddOnList_s { int ___1129; };
 #endif
typedef struct AddOnList_s* ___11; typedef struct ___2730* ___2727; typedef struct ___3870*   ___3869; typedef struct ___834*  ___833; typedef struct ___3867* ___3866; typedef struct ___2755*       ___2754; typedef struct ___834* ElementOrientation_pa;
 #define ___1989 (-1)
 #define ___2749 (-1)
 #define ___2750    (-1)
 #define ___4273  (___2749-2) 
typedef struct ___1293* ___1292; typedef struct ___1272* ___1271; typedef struct ___1152* ___1151; typedef struct ___2743* ___2742; enum RecordingLangauge_e { RecordingLanguage_TecplotMacro, RecordingLanguage_Python, END_RecordingLangauge_e, RecordingLanguage_Invalid = ___329 }; enum FaceNeighborMemberArray_e { ___1282, ___1280, ___1278, ___1279, ___1277, ___1281, END_FaceNeighborMemberArray_e, ___1283 = ___329 }; int const ___1291 = (int)END_FaceNeighborMemberArray_e; enum FaceMapMemberArray_e { ___1266, ___1267, ___1265, ___1268, ___1264, ___1262, ___1263, END_FaceMapMemberArray_e, ___1269 = ___329 }; int const ___1270 = (int)END_FaceMapMemberArray_e; enum ElemToFaceMapMemberArray_e { ___1147, ___1148, END_ElemToFaceMapMemberArray_e, ___1149 = ___329 }; int const ___1150 = (int)END_ElemToFaceMapMemberArray_e; enum NodeToElemMapMemberArray_e { ___2739, ___2740, END_NodeToElemMapMemberArray_e, ___2738 = ___329 }; int const ___2741 = (int)END_NodeToElemMapMemberArray_e; typedef struct ___1362* ___1361; typedef struct AuxData_s* ___264; typedef enum { ___929, ___930, ___931, END_DataValueStructure_e, ___933 = (END_DataValueStructure_e - 1), ___932 = ___329 } DataValueStructure_e; typedef enum { ___871, ___872, END_DataNodeStructure_e, ___873 = ___329 } DataNodeStructure_e; typedef enum { ___4361, ___4359, END_VarLockMode_e, ___4360 = ___329 } VarLockMode_e; typedef enum { ___1376, ___1377, END_FieldMapMode_e, ___1375 = ___329 } FieldMapMode_e; typedef enum { ___4274, ___4277, ___4276, END_UnloadStrategy_e, ___4275 = ___329 } UnloadStrategy_e; struct ___4648 { ___516       ___3174; ___372          ___2027; }; struct ___4683 { ___4264         ___4263; ___4654         ___2686; ___1172         ___2975; ___1172         ___3786; double             ___3641; ___2227          ___2830; ___2227          ___2831; ___2227          ___2832; ___2227          ___1834; ___2227          ___2107; ___2227          ___2136; ZoneType_e         ___4236; ___372          ___2064; ___4648     ___4647; ___264         ___230; ___372          ___419; FaceNeighborMode_e ___1440; ___372          ___228; ___2227          ___2805; ___2227          ___2799; ___2227          ___2800; }; typedef struct GenericImage_s* ___1554; struct ___4077 { TextBox_e        ___411; double           ___2338; double           ___2290; ___516     ___351; ___516     ___1410; }; struct ___4118 { ___4264       ___4263; ___54      ___52; CoordSys_e       ___3167; ___1172       ___4600; ___372        ___227; ___516     ___351; ___4123      ___4121; ___4077        ___401; double           ___57; TextAnchor_e     ___39; double           ___2288; Scope_e          ___3443; char*            ___2331; Clipping_e       ___496; char*            Text; TextType_e       TextType; struct ___4118*   ___2705; struct ___4118*   ___3182; }; struct ___1552 { ___1361  ___4293; ___1361  ___4295; ___1361  ___4297; }; struct ___3153 { ___1361  ___4142; ___1361  ___3275; }; struct ___448 { ___1361  ___4569; ___1361  ___4586; ___1361  ___4595; }; union ___1575 { ___1552   ___1548; ___448 ___4580; ___3153     ___4147; }; struct ___1632 { ___1632(); struct GeomAnchor { public: GeomAnchor(___1632& ___1555); double XOrTheta() const; void setXOrTheta(double ___4314); double YOrR() const; void setYOrR(double); double Z() const; void setZ(double); ___54 anchorPosition() const; bool operator==(GeomAnchor const&) const; bool operator!=(GeomAnchor const&) const; GeomAnchor& operator=(GeomAnchor const&); private: ___1632& m_outer; double m_XOrTheta;
double m_YOrR; mutable double m_Z; mutable bool m_zPosHasBeenAssigned; }; ___4264              ___4263; GeomType_e              ___1652; CoordSys_e              ___3167; GeomAnchor              position; ___372               ___227; ___1172              ___4600; ___516            ___351; ___372               ___2023; ___516            ___1410; LinePattern_e           ___2264; double                  ___2987; double                  ___2290; Scope_e                 ___3443; DrawOrder_e             ___1113; Clipping_e              ___496; FieldDataType_e         ___907; char                   *___2331; ArrowheadStyle_e        ___188; ArrowheadAttachment_e   ___176; double                  ___187; double                  ___171; int32_t                 ___2794; char*                   ___1884; char*                   WorldFileName; ___2227               EmbeddedLpkImageNumber; ___372               ___2333; double                  ___3089; int32_t                 ___2836; ___3461           ___2838; ___1575              ___1573; ImageResizeFilter_e     ___1890; ___1554         ___1882; struct ___1632*          ___2704; struct ___1632*          ___3176; double _WorldFileAssignedWidth; double _WorldFileAssignedHeight; double _WorldFileAssignedXPos; double _WorldFileAssignedYPos; }; typedef struct ___4118* ___4114; typedef struct ___1632* ___1624; typedef enum { MarchingCubeAlgorithm_Classic, MarchingCubeAlgorithm_ClassicPlus, MarchingCubeAlgorithm_MC33, END_MarchingCubeAlgorithm_e, MarchingCubeAlgorithm_Invalid = ___329 } MarchingCubeAlgorithm_e; typedef enum { DataStoreStrategy_Auto, DataStoreStrategy_Heap, END_DataStoreStrategy_e, DataStoreStrategy_Invalid = ___329 } DataStoreStrategy_e; typedef enum { ArgListArgType_ArbParamPtr, ArgListArgType_DoublePtr, ArgListArgType_ArbParam, ArgListArgType_Array, ArgListArgType_Double, ArgListArgType_Function, ArgListArgType_Int, ArgListArgType_Set, ArgListArgType_String, ArgListArgType_StringList, ArgListArgType_StringPtr, END_ArgListArgType_e, ArgListArgType_Invalid = ___329 } ArgListArgType_e; typedef enum { StateModernizationLevel_Latest, StateModernizationLevel_LatestFocus, StateModernizationLevel_2006, StateModernizationLevel_2012, StateModernizationLevel_2012Focus, StateModernizationLevel_2013, StateModernizationLevel_2013Focus, StateModernizationLevel_2014, StateModernizationLevel_2014Focus, StateModernizationLevel_2016R1, StateModernizationLevel_2016R1Focus, StateModernizationLevel_2018R3, StateModernizationLevel_2018R3Focus, END_StateModernizationLevel_e, StateModernizationLevel_Invalid = ___329 } StateModernizationLevel_e; typedef ___372 (STDCALL *___2908)(___3839 ___2907, ___90    ___3325); typedef void (STDCALL *___2909)(___90 ___2906, ___90 ___3325); typedef void (STDCALL *___2910)(___90 ___2906, ___90 ___3325); typedef ___372 (STDCALL *___2862)(int32_t                  ___4459, int32_t                  ___1827, ExportRegion_e           ___1237, ImageRenderingStrategy_e imageRenderingStrategy, ___90               ___3325, TP_OUT ___90*       ___1886); typedef void (STDCALL *___2863)(___90 ___1886, ___90 ___3325); typedef ___372 (STDCALL *___2864)(___90           ___1886, int32_t              ___3424, ___90           ___3325, TP_ARRAY_OUT uint8_t* ___3300, TP_ARRAY_OUT uint8_t* ___1809, TP_ARRAY_OUT uint8_t* ___365); typedef void (STDCALL *OffscreenImageClearCacheCallback_pf)(___90 ___3325);
 #if defined MSWIN
typedef ___372 (STDCALL *___4467)(HDC        ___3188, ___90 ___1886, Palette_e  ___2917, ___90 ___3325);
 #endif 
 #if defined MSWIN
typedef HDC(STDCALL *___4466)(___90 ___3325);
 #endif 
typedef ___372 (STDCALL *___3336)(PageRenderDest_e ___2911, ___90       ___3337, ___90       ___3325); typedef ___372 (STDCALL *___3350)(___90 ___2906, ___90 ___3325); typedef void (STDCALL *___3347)(___90      ___2906, ___90      ___3325, TP_OUT int32_t* ___4459, TP_OUT int32_t* ___1827); typedef void (STDCALL *___3909)(___90 ___2906, ___90 ___3325); typedef void (STDCALL *___2156)(___90        ___3325, TP_OUT ___372* ___2058, TP_OUT ___372* ___2000, TP_OUT ___372* ___2010); typedef ___372 (STDCALL *___2579)(int        ___420, ___90 ___3325); typedef void (STDCALL *___4445)(___372  ___2, ___90 ___3325); typedef void (STDCALL *___335)(CursorStyle_e ___694, ___90    ___3349, ___90    ___3325); typedef void (STDCALL *___3472)(double ___4567, double ___4584, double ___4568, double ___4585, ___372     ___514, ___90    ___3325); typedef void (STDCALL *___3202)(___90 ___3325); typedef void (STDCALL *___1168)(___90 ___3325); typedef ___372 (STDCALL *___1024)(___90 ___3325); typedef void (STDCALL *___997)(___90 ___3325); typedef void (STDCALL *___1101)(___90     ___3325, TP_OUT double* ___1838, TP_OUT double* ___2110); typedef void (STDCALL *___3451)(___90  ___3325, TP_OUT int* ___4460, TP_OUT int* ___1828); typedef ___372(STDCALL *___1014)(const char* ___1959, char**      ___4331, ___372   ___3173, ___90  ___3324); typedef MessageBoxReply_e(STDCALL *___1033)(const char*      ___2456, MessageBoxType_e ___2442, ___90       ___3325); typedef ___372 (STDCALL *___1064)(SelectFileOption_e ___1038, const char*        ___1077, const char*        ___950, const char*        ___949, TP_GIVES char**    ___3363, ___90         ___3325); typedef void (STDCALL *___3777)(const char* ___3780, ___90  ___3325); typedef void (STDCALL *___3581)(const char* ___1421, const char* ___3458, ___90  ___3325); typedef void (STDCALL *___3243)(int        ___3246, ___90 ___3325); typedef void (STDCALL *___3245)(___372  ___3580, ___372  ___2029, ___90 ___3325); typedef void (STDCALL *___3244)(___90 ___3325); typedef ___372 (STDCALL *___12)(___90 ___494); typedef ___372 (STDCALL *___4176)(___12  ___4175, ___90             ___494, uint32_t               ___1984, ___90             ___3325); typedef void (STDCALL *___2618)(const char* ___2207, ___372 operationSucceeded, ___90  ___3325); typedef const char* (STDCALL *___2880)(___90  ___3325); typedef void (STDCALL *___2410)(___90 ___3325); typedef void (STDCALL *___2411)(___90 ___3325); typedef ___372 (STDCALL *___2412)(___90 ___3325); typedef ___372 (STDCALL *___2413)(___90 ___3325); typedef void (STDCALL *___3193)(___372 ___2034); typedef void (STDCALL *___3194)(___372  ___4451, ___372  ___2034, ___90 ___494); typedef void (STDCALL *___1137)(void); typedef ___372 (STDCALL *___1112)(RedrawReason_e ___3302, ___90     ___494); typedef ___372 (STDCALL *WriteLayoutPreWriteCallback_pf)(___3501 PageList, ___90 ___494); typedef int (STDCALL *___3844)(const char* ___3814, const char* ___3815, ___90  ___494); typedef double(STDCALL *___1383)(const ___1361 ___1309, ___2227          ___3249); typedef void (STDCALL *___1384)(___1361 ___1309, ___2227    ___3249, double       ___4298); typedef ___372 (STDCALL *___2311)(___1361 ___1352); typedef ___372 (STDCALL *___2312)(___1361 ___1352); typedef void (STDCALL *___2310)(___1361 ___1352);
typedef ___372 (STDCALL *___2308)(___2727 ___2724); typedef ___372 (STDCALL *___2309)(___2727 ___2724); typedef void (STDCALL *___2307)(___2727 ___2724); typedef ___372 (STDCALL *___2305)(___1292 ___1275); typedef ___372 (STDCALL *___2306)(___1292 ___1275); typedef void (STDCALL *___2304)(___1292 ___1275); typedef ___372 (STDCALL *___2302)(___1271 ___1260); typedef ___372 (STDCALL *___2303)(___1271 ___1260); typedef void (STDCALL *___2301)(___1271 ___1260); typedef void (STDCALL *___1250)(___2227 ___2829, double*   ___4577, double*   ___4594); typedef void (STDCALL *___3469)(void); typedef void (STDCALL *___653)(const char*  ___3177, const char*  ___3178, const ___3501 ___3179); typedef ___372 (STDCALL *___888)(char*           ___850, char*           ___4041, TP_GIVES char** ___2456); typedef ___372 (STDCALL *___898)(___3839 ___1960); typedef ___372 (STDCALL *___861)(___3839 ___1960, ___90    ___494); typedef ___1137 ___862; typedef void (STDCALL * ___863)( ___3839 ___3463, ___372     ___21, ___372     ___82, ___90    ___494); typedef void (STDCALL * ___864)( ___3839 ___3463, ___372     ___21, ___90    ___494); typedef ___372 (STDCALL *___897)(___3839  ___1960); typedef ___372 (STDCALL *___860)(___3839 ___1960, ___90    ___494); typedef ___372 (STDCALL * ___934)(const char*   ___3462, ___90    ___494); typedef void (STDCALL *___1664)(___3501        ___2255, ___3839 ___3464); typedef void (STDCALL *___1660)(___1172      ___2253, char*           ___724, TP_GIVES char** ___0); typedef ___372 (STDCALL *___1663)(___1361    ___3274, ___1361    ___3273, CoordScale_e    ___1930, CoordScale_e    ___963, ___2227       ___2833, ___1172      ___2253, char*           ___724, TP_GIVES char** ___723); typedef ___372 (STDCALL *___1752)(___1361   ___3274, ___1361   ___3273, CoordScale_e   ___1930, CoordScale_e   ___963, ___2227      ___2833, ___2227      ___2790, ___1172     ___2253, char*          ___724, TP_OUT double* ___1925, TP_OUT double* ___962);
 #if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
typedef ___1752 ___1801;
 #endif
typedef ___372 (STDCALL *___1770)(___1361   ___3274, ___1361   ___3273, CoordScale_e   ___1930, CoordScale_e   ___963, ___2227      ___2833, ___2227      ___2790, ___1172     ___2254, char*          ___724, double         ___3195, TP_OUT double* ___3192);
 #if defined MSWIN
typedef ___372 (STDCALL *___3175)(MSG *___3122);
 #endif
typedef ___372 (STDCALL *___1136)(double          ___4315, ___90      ___494, TP_GIVES char** ___2174); typedef void (STDCALL *___2874)(___90 ___494); typedef ___372 (STDCALL *___3452)(const char *___3453, ___90  ___494); typedef ___372 (STDCALL * ___2281)(___2227         ___4452, ___1172        ___4600, ___2227         ___450, ___2227         ___1252, TP_GIVES double * ___3161, ___90        ___494); typedef ___372 (STDCALL *___4240)(const uint32_t* ___2213, const uint32_t* ___3395, ___90      ___494); typedef ___372 (STDCALL *___4241)(const uint64_t* ___2213, const uint64_t* ___3395, ___90      ___494); typedef void (*TUAbort_pf)(const char* error_message); typedef ___3839 (STDCALL * MatchVariablesCallback_pf)(___3839         existingVariables, ___3839 const * incomingVariableLists, ___2227             numIncomingVariableLists, ___90            clientData);
 #define ___1478          0
 #define ___1477          1
 #define ___1480              2
 #define ___1479 3
typedef struct ViewState_s* ___3432; typedef struct ViewState_s* ___4436; typedef struct ProbeInfo_s* ___3197; static const char* const ___4031 = "support@tecplot.com";
 #define ___4035 0 
 #define TECUTIL_DEFAULT_TEXT_ID -1
typedef ___90     TextID_t; typedef ___90     GeomID_t;
 #endif 
