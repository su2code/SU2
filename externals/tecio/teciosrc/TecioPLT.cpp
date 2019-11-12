#include "stdafx.h"
#include "MASTER.h"
using std::vector;
#include "GLOBAL.h"
#include "TASSERT.h"
#include "SYSTEM.h"
#include "FILESTREAM.h"
#include "TecplotVersion.h"
#include "CHARTYPE.h"
#include "DATAIO4.h"
#include "DATASET0.h"
#include "TecioPLT.h"
#include "DATAUTIL.h"
#include "ALLOC.h"
#include <vector>
#include "FileSystem.h"
 #if !defined MAKEARCHIVE
#include "AUXDATA.h"
 #endif 
 #if defined MSWIN
#include <io.h>
 #endif
 #if defined UNIXX
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
 #endif
 #define TECIO_STRINGIZE(s) #s
 #if defined MAKEARCHIVE
 #if defined MSWIN && defined _DEBUG
 #define ___3184(s) do { OutputDebugStringA(s); } while (0)
 #define ___3185(s,a1) do { char ___416[512]; sprintf(___416,s,a1); OutputDebugStringA(___416); } while (0)
 #define ___3186(s,a1,a2) do { char ___416[512]; sprintf(___416,s,a1,a2); OutputDebugStringA(___416); } while (0)
 #define ___3187(s,a1,a2,a3) do { char ___416[512]; sprintf(___416,s,a1,a2,a3); OutputDebugStringA(___416); } while (0)
 #else
 #define ___3184(s) printf(s)
 #define ___3185(s,a1) printf(s,a1)
 #define ___3186(s,a1,a2) printf(s,a1,a2)
 #define ___3187(s,a1,a2,a3) printf(s,a1,a2,a3)
 #endif
 #else
 #if defined MSWIN
 #define ___3184(s) ((void)0)
 #define ___3185(s,a1) ((void)0)
 #define ___3186(s,a1,a2) ((void)0)
 #define ___3187(s,a1,a2,a3) ((void)0)
 #else
 #define ___3184(s) printf(s)
 #define ___3185(s,a1) printf(s,a1)
 #define ___3186(s,a1,a2) printf(s,a1,a2)
 #define ___3187(s,a1,a2,a3) printf(s,a1,a2,a3)
 #endif
 #endif
 #define ___2380    10
 #define BYTES_PER_CHUNK 4096
 #define TECIO_NO_NEIGHBORING_ELEM 0
 #define TECIO_NO_NEIGHBORING_ZONE 0
namespace {
 #if defined MAKEARCHIVE
___2227            DebugLevel[___2380] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 #endif
int32_t             ___2042[___2380]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; int32_t             NumErrs[___2380] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; int32_t             ___2847[___2380] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; int32_t             NumVars[___2380]; char*                DestFName[___2380] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }; char*                BlckFName[___2380] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }; ___1405*        BlckFile[___2380]; ___1405*        HeadFile[___2380]; vector<___1398> MinMaxOffset[___2380]; vector<double>       VarMinValue[___2380]; vector<double>       VarMaxValue[___2380]; ___372            DoWriteForeign = ___1305; ___372            IsWritingNative[___2380]; ___372            ___2006[___2380]; int32_t             ___4693[___2380]; ___2227            ___1906[___2380]; ___2227            ___2114[___2380]; ___2227            ___2159[___2380]; vector<___2227>    ___4193[___2380]; ___2227            ___4191[___2380]; ___2227            ___4190[___2380]; struct PolyZoneWriteInfo { PolyZoneWriteInfo() : numFacesWritten(0) , faceNodeSum(0) , numBoundaryFacesWritten(0) , boundaryConnectionSum(0) , numFaceNodesOffset(0) , faceNodesOffset(0) , leftElemsOffset(0) , rightElemsOffset(0) , connectionCountsOffset(0) , connectionElemsOffset(0) , connectionZonesOffset(0) {} ___2227 numFacesWritten; ___2227 faceNodeSum; ___2227 numBoundaryFacesWritten; ___2227 boundaryConnectionSum; ___1398 numFaceNodesOffset; ___1398 faceNodesOffset; ___1398 leftElemsOffset; ___1398 rightElemsOffset; ___1398 connectionCountsOffset; ___1398 connectionElemsOffset; ___1398 connectionZonesOffset; }; vector<PolyZoneWriteInfo> PolyZoneWriteInfos[___2380]; ___2227            ICellMax[___2380]; ___2227            JCellMax[___2380]; ___2227            KCellMax[___2380]; vector<int32_t>     ___2803[___2380]; int32_t             ___1285[___2380]; vector<int32_t>     FaceNeighborsOrMapWritten[___2380]; int32_t             NumIndices[___2380]; ___2227            NumDataValuesWritten[___2380]; ___2227            NumOrderedCCDataValuesWritten[___2380]; ___2227            NumDataValuesToWrite[___2380]; vector<___2227>    NumRunningVarValues[___2380]; vector<___372>    IsSharedVar[___2380]; vector<___372>    IsPassiveVar[___2380]; int32_t             ___734[___2380]; int32_t             ___718[___2380]; int32_t             ___1364[___2380]; vector<___372>    IsCellCentered[___2380]; ___372            HasFECONNECT[___2380]; int32_t             FileTypes[___2380]; vector<___2227>    NumConnectivityNodes[___2380]; vector<___2227>    NumConnectivityNodesWritten[___2380]; vector<___372>    ConnectivityWritten[___2380]; }
 #define ___2882 0
 #define ___1340 1
 #define ___1348 2
 #define ___1343 3
 #define ___1347 4
 #define FEBRICK 5
 #define ___1341 6
 #define ___1342 7
 #define ___1535 0
 #define ___1811 1
 #define ___3638 2
 #if defined MAKEARCHIVE
namespace { char const* ZoneTypes[] = { "ORDERED", "FELINESEG", "FETRIANGLE", "FEQUADRILATERAL", "FETETRAHEDRON", "FEBRICK", "FEPOLYGON", "FEPOLYHEDRON" }; }
 #endif 
namespace { void WriteErr( int32_t    ___1397, char const* routine_name) {
 #if defined MAKEARCHIVE
{ ___3186("Err: (%s) Write failure on file %d.\n", routine_name, ___1397 + 1); }
 #else
{ ___4278(routine_name); }
 #endif
NumErrs[___1397]++; } } namespace { ___1405* OpenFileStream( char const* FilePath, char const* AccessMode, ___372   ___2007) { REQUIRE(VALID_REF(FilePath)); REQUIRE(VALID_REF(AccessMode)); ___1405 *___3359 = NULL; FILE         *File   = tecplot::filesystem::fileOpen(FilePath, AccessMode); if (File != NULL) { ___3359 = ___1402(File, ___2007); if (___3359 == NULL) ___4195(File); } ENSURE((VALID_REF(___3359) && VALID_REF(___3359->File)) || ___3359 == NULL); return ___3359; } } namespace { void CloseFileStream(___1405** ___1401) { REQUIRE(VALID_REF(___1401)); REQUIRE(VALID_REF(*___1401) || *___1401 == NULL); if (*___1401 != NULL) { ___4195((*___1401)->File); ___1403(___1401); } ENSURE(*___1401 == NULL); } } namespace { char getBestTerminatorChar(char const* CompoundStr) { REQUIRE(VALID_REF(CompoundStr)); if (strchr(CompoundStr, '\n') != NULL) return '\n'; else if (strchr(CompoundStr, ',') != NULL) return ','; else return ' '; } } int32_t ___3982(int32_t ___1397) { return ___2847[___1397]; } int32_t ___3981(int32_t ___1397) { return NumVars[___1397]; } int32_t ___3977( int32_t        ___1397, char const*     ___4178, char const*     ___4351, char const*     ___1439, char const*     ___3448, int32_t const* ___1408, int32_t const* ___942, int32_t const* ___4439) {
 #if defined MAKEARCHIVE
___1935();
 #endif
if (___2042[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECINI142) file %d is already open.\n", ___1397+1);
 #endif
return (-1); } ___2847[___1397] = 0;
 #if defined MAKEARCHIVE
{ DebugLevel[___1397] = *___942; }
 #else
{ ___4278(___942); }
 #endif
___478(VarMinValue[___1397].empty()); ___478(VarMaxValue[___1397].empty()); ___478(NumRunningVarValues[___1397].empty()); ___478(IsSharedVar[___1397].empty()); ___478(IsPassiveVar[___1397].empty()); ___478(IsCellCentered[___1397].empty()); ___478(MinMaxOffset[___1397].empty()); ___478(___4193[___1397].empty()); ___478(PolyZoneWriteInfos[___1397].empty()); ___478(___2803[___1397].empty()); ___478(FaceNeighborsOrMapWritten[___1397].empty()); ___478(NumConnectivityNodes[___1397].empty()); ___478(NumConnectivityNodesWritten[___1397].empty()); ___478(ConnectivityWritten[___1397].empty()); ___734[___1397] = -1; size_t ___2165 = 0; if (___1439 != NULL) ___2165 = strlen(___1439); if (___2165 == 0) {
 #if defined MAKEARCHIVE
___3185("Err: (TECINI142) Bad file name for file %d.\n", ___1397+1);
 #endif
return (-1); } DestFName[___1397] = ___23(___2165 + 1, char, "data set fname"); strcpy(DestFName[___1397], ___1439); char RName[80];
 #if defined (___1100)
{ sprintf(RName, "BLCKFILE.%03d", (int)(___1397 + 1)); }
 #else
{ sprintf(RName, "tp%1dXXXXXX", ___1397 + 1); }
 #endif
___2165 = strlen(RName); if (___3448 != NULL) ___2165 += strlen(___3448) + 1; BlckFName[___1397] = ___23(___2165 + 1, char, "data set fname"); if (___3448 != NULL) { strcpy(BlckFName[___1397], ___3448);
 #if defined ___1100 || defined MSWIN
{ strcat(BlckFName[___1397], "\\"); }
 #else
{ strcat(BlckFName[___1397], "/"); }
 #endif
} else { BlckFName[___1397][0] = '\0'; } strcat(BlckFName[___1397], RName); ___478(strlen(BlckFName[___1397]) <= ___2165);
 #if defined MSWIN
{ if (!_mktemp(BlckFName[___1397])) { WriteErr(___1397, "Error creating new file name."); return (-1); } }
 #elif defined UNIXX
{ mode_t OrigUmask = umask(0022); int FileDesc = mkstemp(BlckFName[___1397]); if (FileDesc != -1) close(FileDesc); umask(OrigUmask); }
 #endif
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) { ___3186("Scratch File #%d: %s\n", ___1397 + 1, BlckFName[___1397]); ___3186("Dest    File #%d: %s\n", ___1397 + 1, DestFName[___1397]); }
 #endif
IsWritingNative[___1397] = !DoWriteForeign; HeadFile[___1397] = OpenFileStream(DestFName[___1397], "wb", IsWritingNative[___1397]); BlckFile[___1397] = OpenFileStream(BlckFName[___1397], "wb", IsWritingNative[___1397]); if (BlckFile[___1397] == NULL) {
 #if defined MAKEARCHIVE
___3184("Err: (TECINI142) Cannot open scratch file for output.\n"); ___3184("     Check permissions in scratch directory.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (HeadFile[___1397] == NULL) {
 #if defined MAKEARCHIVE
___3184("Err: (TECINI142) Cannot open plot file.  Check permissions.\n");
 #endif
NumErrs[___1397]++; return (-1); } ___4494(*HeadFile[___1397], TecplotSDKBinaryFileVersion); ___4492(HeadFile[___1397]); if (*___1408 >= ___1535 && *___1408 <= ___3638) FileTypes[___1397] = *___1408; else {
 #if defined MAKEARCHIVE
___3184("Err: (TECINI142) Bad filetype argument.  Check documentation.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (!___4490(HeadFile[___1397], (___2227)FileTypes[___1397])) { WriteErr(___1397, "TECINI142"); return (-1); } if (!___1131(HeadFile[___1397], ___4178, ___4226)) { WriteErr(___1397, "TECINI142"); return (-1); } NumVars[___1397] = 0; char const* ___685 = ___4351; { char terminator = getBestTerminatorChar(___685); while (*___685) { while (*___685 && *___685 == ' ') ___685++; if (*___685) { NumVars[___1397]++; while (*___685 && *___685 != terminator) ___685++; if (*___685) ___685++; } } }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3185("NumVars=%d\n", NumVars[___1397]);
 #endif
try { VarMinValue[___1397].resize(NumVars[___1397]); VarMaxValue[___1397].resize(NumVars[___1397]); NumRunningVarValues[___1397].resize(NumVars[___1397]); IsSharedVar[___1397].resize(NumVars[___1397]); IsPassiveVar[___1397].resize(NumVars[___1397]); IsCellCentered[___1397].resize(NumVars[___1397]); } catch (std::bad_alloc const&) {
 #if defined MAKEARCHIVE
___3184("Err: (TECINI142) Memory allocation error.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (!___4490(HeadFile[___1397], (___2227)NumVars[___1397])) { WriteErr(___1397, "TECINI142"); return (-1); } ___685 = ___4351; { char TString[___2356+1]; int ___1832; char terminator = getBestTerminatorChar(___685); while (*___685) { while (*___685 && *___685 == ' ') ___685++; if (*___685) { ___1832 = 0; while (*___685 && *___685 != terminator) { TString[___1832++] = *___685++; } if (*___685) ___685++; ___1832--; while (___1832 >= 0 && TString[___1832] == ' ') ___1832--; TString[___1832+1] = '\0'; if (!___1131(HeadFile[___1397], TString, ___4226)) { WriteErr(___1397, "TECINI142"); return (-1); } } } } ___2042[___1397] = 1; if (*___4439) ___1364[___1397] = FieldDataType_Double; else ___1364[___1397] = FieldDataType_Float; return (0); } namespace { int CheckData( int32_t    ___1397, char const* routine_name) { if (NumDataValuesToWrite[___1397] != NumDataValuesWritten[___1397]) {
 #if defined MAKEARCHIVE
{ ___3186("Err: (%s) Wrong number of data values in file %d:\n", routine_name, ___1397 + 1); ___3186("     %lld data values for Zone %d were processed,\n", (long long)NumDataValuesWritten[___1397], ___734[___1397] + 1); ___3185("     %lld data values were expected.\n", (long long)NumDataValuesToWrite[___1397]); }
 #else
{ ___4278(routine_name); }
 #endif
NumErrs[___1397]++; return (-1); } return (0); } } namespace { int CheckFile( int32_t    ___1397, char const* routine_name) { if ((___1397 < 0) || (___1397 >= ___2380)) { ___3187("Err: (%s) Attempt to use invalid file %d file must be between 1 and %d inclusive.\n", routine_name, ___1397+1, ___2380); return (-1); } if (!___2042[___1397]) { ___3186("Err: (%s) Attempt to use file %d that hasn't been initialized with TECINI142.\n", routine_name, ___1397+1); return (-1); } return (0); } } namespace { void AdvanceToNextVarWithValues(int32_t ___1397) { do { ___718[___1397]++; } while (___718[___1397] < NumVars[___1397] && (IsSharedVar[___1397][___718[___1397]] || IsPassiveVar[___1397][___718[___1397]])); } } int32_t ___3992( int32_t        ___1397, char const*     ___4598, int32_t const* ___4599, int32_t const* ___1910, int32_t const* ___2117, int32_t const* ___2162, int32_t const* ___1836, int32_t const* ___2109, int32_t const* ___2138, double const*   ___3641, int32_t const* ___3786, int32_t const* ___2975, int32_t const* ___2005, int32_t const* ___2801, int32_t const* ___1440, int32_t const* ___2805, int32_t const* ___2799, int32_t const* ___2798, int32_t const* ___2983, int32_t const* ___4327, int32_t const* ___3552, int32_t const* ___3550) { int        ___1832; if (CheckFile(___1397, "TECZNE142") < 0) return (-1); if (___734[___1397] > -1) { if (CheckData(___1397, "TECZNE142") < 0) return (-1); } if (NumVars[___1397] == 0) { WriteErr(___1397, "TECZNE142");
 #if defined MAKEARCHIVE
___3185("Err: (TECZNE142) Cannot write out zones if numvars is equal to zero (file %d).\n", ___1397 + 1);
 #endif
return (-1); } if (___734[___1397] > ___2382 - 2) { WriteErr(___1397, "TECZNE142");
 #if defined MAKEARCHIVE
___3186("Err: (TECZNE142) Exceeded max number of zones (%d) in file %d.\n", ___2382, ___1397 + 1);
 #endif
return (-1); } if (*___3786 < -1) {
 #if defined MAKEARCHIVE
___3186("Err: (TECZNE142) Invalid StrandID supplied for file %d, zone %d.\n", ___1397 + 1, ___734[___1397] + 1 + 1);
 #endif
return (-1); } if (*___2975 < 0) {
 #if defined MAKEARCHIVE
___3186("Err: (TECZNE142) Invalid ParentZone supplied for file %d, zone %d.\n", ___1397 + 1, ___734[___1397] + 1 + 1);
 #endif
return (-1); } if (*___2005 != 1) {
 #if defined MAKEARCHIVE
___3186("Err: (TECZNE142) Point data is not currently allowed. " " Please use block format for file %d, zone %d.\n", ___1397 + 1, ___734[___1397] + 1 + 1);
 #endif
return (-1); } NumDataValuesWritten[___1397]          = 0; NumOrderedCCDataValuesWritten[___1397] = 0; ___734[___1397]++; try { MinMaxOffset[___1397].resize(___734[___1397] + 1); ___4193[___1397].resize(___734[___1397] + 1); PolyZoneWriteInfos[___1397].resize(___734[___1397] + 1); ___2803[___1397].resize(___734[___1397] + 1); FaceNeighborsOrMapWritten[___1397].resize(___734[___1397] + 1); NumConnectivityNodes[___1397].resize(___734[___1397] + 1); NumConnectivityNodesWritten[___1397].resize(___734[___1397] + 1); ConnectivityWritten[___1397].resize(___734[___1397] + 1); } catch (std::bad_alloc const&) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZNE142) Memory allocation error.\n");
 #endif
NumErrs[___1397]++; return (-1); } ___4693[___1397] = *___4599; ___1906[___1397] = *___1910; ___2114[___1397] = *___2117; ___2159[___1397] = *___2162; ICellMax[___1397] = *___1836; JCellMax[___1397] = *___2109; KCellMax[___1397] = *___2138; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___1305; NumConnectivityNodesWritten[___1397][___734[___1397]] = 0; ConnectivityWritten[___1397][___734[___1397]] = ___1305; if (___4693[___1397] == ___4697) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZNE142) FE Mixed Volume Zone type not yet supported.\n");
 #endif
NumErrs[___1397]++; return (-1); } else if (___4693[___1397] == ___4698 || ___4693[___1397] == ___4699) { ___2803[___1397][___734[___1397]] = 0; ___1285[___1397]   = 0; NumConnectivityNodes[___1397][___734[___1397]] = 0; ___2006[___1397]                             = ___4226; ___4193[___1397][___734[___1397]] = *___2805; ___4191[___1397]              = *___2799; ___4190[___1397]              = *___2798; if (*___3550) { if (*___2799 > 0 || *___2798 > 0) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZNE142) Cannot share poly zone connectivity if there are connected boundary faces.\n");
 #endif
NumErrs[___1397]++; return (-1); } else { FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; } } } else { ___2006[___1397]                              = (*___2005 != 0); ___2803[___1397][___734[___1397]] = *___2801; ___1285[___1397]                     = *___1440; ___4193[___1397][___734[___1397]] = 0; ___4191[___1397]              = 0; ___4190[___1397]              = 0; } ___4493(HeadFile[___1397], (double)___4649, FieldDataType_Float); if (!___1131(HeadFile[___1397], ___4598, ___4226)) { WriteErr(___1397, "TECZNE142"); return (-1); } if (___734[___1397] == 0) { if (___3552) { bool aVarIsSharing = false; for (int32_t ___4291 = 0; ___4291 < NumVars[___1397]; ++___4291) { if (___3552[___4291]) { aVarIsSharing = true; break; } } if (aVarIsSharing) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZNE142) First zone cannot share variables because there is no zone to share from.\n");
 #endif
NumErrs[___1397]++; return (-1); } } if (*___3550) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZNE142) First zone cannot share connectivity because there is no zone to share from.\n");
 #endif
NumErrs[___1397]++; return (-1); } } switch (___4693[___1397]) { case ___2882: NumIndices[___1397] = 0; break; case ___1340: NumIndices[___1397] = 2; break; case ___1348: NumIndices[___1397] = 3; break; case ___1343: NumIndices[___1397] = 4; break; case ___1347: NumIndices[___1397] = 4; break; case FEBRICK: NumIndices[___1397] = 8; break; } if (___4693[___1397] != ___4698    && ___4693[___1397] != ___4699 && *___3550 == 0            && FileTypes[___1397] != ___3638) { NumConnectivityNodes[___1397][___734[___1397]] = NumIndices[___1397] * static_cast<___2227>(___2114[___1397]); } ___4490(HeadFile[___1397], (___2227)(*___2975) - 1); ___4490(HeadFile[___1397], (___2227)(*___3786) - 1); ___4493(HeadFile[___1397], *___3641, FieldDataType_Double); ___4490(HeadFile[___1397], (___2227) - 1); ___4490(HeadFile[___1397], ___4693[___1397]); NumDataValuesToWrite[___1397] = 0; for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { IsSharedVar[___1397][___1832]  = (___3552 != NULL && ___3552[___1832] != 0); IsPassiveVar[___1397][___1832] = (___2983   != NULL && ___2983[___1832]   == 1); } ___4490(HeadFile[___1397], (___2227)(___4327 != NULL ? 1 : 0)); if (___4327) { for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { ___2227  NumNodes; ___2227  NumCells; if (___4693[___1397] == ___2882) { NumNodes = ___1906[___1397] * ___2114[___1397] * ___2159[___1397]; NumCells = (MAX(___1906[___1397] - 1, 1) * MAX(___2114[___1397] - 1, 1) * MAX(___2159[___1397] - 1, 1)); } else { NumNodes = ___1906[___1397]; NumCells = ___2114[___1397]; } if (___1832 == 0) NumRunningVarValues[___1397][___1832] = 0; else NumRunningVarValues[___1397][___1832] = NumRunningVarValues[___1397][___1832-1]; IsCellCentered[___1397][___1832] = (___4327[___1832] == ___4328); if (___4327[___1832] == ___4328) { ___4490(HeadFile[___1397], (___2227)1); if (!IsSharedVar[___1397][___1832] && !IsPassiveVar[___1397][___1832]) { NumDataValuesToWrite[___1397]   += NumCells; NumRunningVarValues[___1397][___1832] += NumCells; } } else if (___4327[___1832] == ___4330) { ___4490(HeadFile[___1397], (___2227)0); if (!IsSharedVar[___1397][___1832] && !IsPassiveVar[___1397][___1832]) { NumDataValuesToWrite[___1397]   += NumNodes; NumRunningVarValues[___1397][___1832] += NumNodes; } } else {
 #if defined MAKEARCHIVE
___3186("Err: (TECZNE142) Bad zone value location for file %d, variable %d.\n", ___1397 + 1, ___1832 + 1);
 #endif
NumErrs[___1397]++; return(-1); } } } else { ___2227 NumNodes; if (___4693[___1397] == ___2882) { NumNodes = ___1906[___1397] * ___2114[___1397] * ___2159[___1397]; } else { NumNodes = ___1906[___1397]; } for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { if (___1832 == 0) NumRunningVarValues[___1397][___1832] = 0; else NumRunningVarValues[___1397][___1832] = NumRunningVarValues[___1397][___1832-1]; IsCellCentered[___1397][___1832] = ___1305; if (!IsSharedVar[___1397][___1832] && !IsPassiveVar[___1397][___1832]) { NumDataValuesToWrite[___1397]   += NumNodes; NumRunningVarValues[___1397][___1832] += NumNodes; } } } ___4490(HeadFile[___1397], (___2227)0); ___4490(HeadFile[___1397], (___2227)___2803[___1397][___734[___1397]]); if (___2803[___1397][___734[___1397]] > 0) { ___4490(HeadFile[___1397], (___2227)___1285[___1397]); if (___4693[___1397] != ___2882) ___4490(HeadFile[___1397], (___2227)0); } ___4490(HeadFile[___1397], (___2227)___1906[___1397]); if (___4693[___1397] == ___1341 || ___4693[___1397] == ___1342) { ___4490(HeadFile[___1397], (___2227)___2159[___1397]); ___4490(HeadFile[___1397], (___2227)___4193[___1397][___734[___1397]]); ___478(TecplotSDKBinaryFileVersion == 112); if (___4191[___1397] > 0) { if (___4190[___1397] < ___4191[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECZNE142) There must be at least 1 boundary connection for each boundary face in zone %d.\n", ___734[___1397] + 1); ___3186("     %lld boundary faces and %lld boundary connections were specified.\n", (long long)___4191[___1397], (long long)___4190[___1397]);
 #endif
NumErrs[___1397]++; return(-1); } ___4490(HeadFile[___1397], (___2227)___4191[___1397] + 1); } else ___4490(HeadFile[___1397], (___2227)___4191[___1397]); ___4490(HeadFile[___1397], (___2227)___4190[___1397]); } ___4490(HeadFile[___1397], (___2227)___2114[___1397]); if (___4693[___1397] == ___2882) { ___4490(HeadFile[___1397], (___2227)___2159[___1397]); } else { ___4490(HeadFile[___1397], (___2227)ICellMax[___1397]); ___4490(HeadFile[___1397], (___2227)JCellMax[___1397]); ___4490(HeadFile[___1397], (___2227)KCellMax[___1397]); } ___4490(HeadFile[___1397], (___2227)0); ___4493(BlckFile[___1397], (double)___4649, FieldDataType_Float); for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { if (!___4513(BlckFile[___1397], (FieldDataType_e)___1364[___1397], ___4226)) { WriteErr(___1397, "TECZNE142"); return (-1); } } if (___2983) { ___4490(BlckFile[___1397], 1); for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) ___4490(BlckFile[___1397], ___2983[___1832]); } else ___4490(BlckFile[___1397], 0); ___718[___1397] = -1; AdvanceToNextVarWithValues(___1397); if (___3552) { ___4490(BlckFile[___1397], 1); for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) ___4490(BlckFile[___1397], ___3552[___1832] - 1); } else ___4490(BlckFile[___1397], 0); ___4490(BlckFile[___1397], *___3550 - 1); MinMaxOffset[___1397][___734[___1397]] = (___1398)___4201(BlckFile[___1397]->File); for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { VarMinValue[___1397][___1832] = ___2179; VarMaxValue[___1397][___1832] = -___2179; if (!IsSharedVar[___1397][___1832] && !IsPassiveVar[___1397][___1832]) { ___4493(BlckFile[___1397], 0.0, FieldDataType_Double); ___4493(BlckFile[___1397], 0.0, FieldDataType_Double); } }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) { ___3185("Writing Zone %d:\n", ___734[___1397] + 1); ___3185("      Title = %s\n", ___4598); ___3185("      Type  = %s\n", ZoneTypes[___4693[___1397]]); ___3185("      IMax  = %lld\n", (long long)___1906[___1397]); ___3185("      JMax  = %lld\n", (long long)___2114[___1397]); ___3185("      KMax  = %lld\n", (long long)___2159[___1397]); if (___3552) { char ___1134[1024] = ""; for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { if (___1832 > 0) strcat(___1134, ","); sprintf(&___1134[strlen(___1134)], "%d", ___3552[___1832]); } ___3185("      DupList = %s\n", ___1134); } }
 #endif
return (0); } namespace { void RewritePendingMinMaxValues(int32_t ___1397) { ___1398 CurrentOffset = (___1398)___4201(BlckFile[___1397]->File); ___4200(BlckFile[___1397]->File, MinMaxOffset[___1397][___734[___1397]], SEEK_SET); int ___1832; for (___1832 = 0; ___1832 < NumVars[___1397]; ___1832++) { if (!IsSharedVar[___1397][___1832] && !IsPassiveVar[___1397][___1832]) { ___4493(BlckFile[___1397], VarMinValue[___1397][___1832], FieldDataType_Double); ___4493(BlckFile[___1397], VarMaxValue[___1397][___1832], FieldDataType_Double); } } ___4200(BlckFile[___1397]->File, CurrentOffset, SEEK_SET); } } int32_t ___3972( int32_t        ___1397, int32_t const* N, void const*     ___816, int32_t const* ___2014) { ___2227  ___1832; double    *dptr = (double *)___816; float     *fptr = (float *)___816; if (CheckFile(___1397, "TECDAT142") < 0) return (-1);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397] && (*N > 1)) ___3186("Writing %d values to file %d.\n", *N, ___1397 + 1);
 #endif
for (___1832 = 0; ___1832 < *N; ___1832++) { double ___4315 = (*___2014 == 1 ? dptr[___1832] : fptr[___1832]); if (___4315 < VarMinValue[___1397][___718[___1397]]) VarMinValue[___1397][___718[___1397]] = ___4315; if (___4315 > VarMaxValue[___1397][___718[___1397]]) VarMaxValue[___1397][___718[___1397]] = ___4315; if (!___4493(BlckFile[___1397], ___4315, (FieldDataType_e)___1364[___1397])) { WriteErr(___1397, "TECDAT142"); return (-1); } if (IsCellCentered[___1397][___718[___1397]] && ___4693[___1397] == ___2882) { ___478(___2006[___1397]); ___2227 PIndex = (NumOrderedCCDataValuesWritten[___1397]); ___2227 FinalIMax = MAX(___1906[___1397] - 1, 1); ___2227 FinalJMax = MAX(___2114[___1397] - 1, 1); ___2227 FinalKMax = MAX(___2159[___1397] - 1, 1); ___2227 IIndex = (PIndex % ___1906[___1397]); ___2227 JIndex = ((PIndex % (___1906[___1397] * ___2114[___1397])) / ___1906[___1397]); ___2227 KIndex = (PIndex / (___1906[___1397] * ___2114[___1397])); ___2227 IMaxAdjust = 0; ___2227 JMaxAdjust = 0; ___2227 KMaxAdjust = 0; if (___2159[___1397] > 1) KMaxAdjust = 1; else if (___2114[___1397] > 1) JMaxAdjust = 1; else if (___1906[___1397] > 1) IMaxAdjust = 1; if (IIndex + 1 == FinalIMax && FinalIMax < ___1906[___1397] - IMaxAdjust) { NumOrderedCCDataValuesWritten[___1397]++; if (!___4493(BlckFile[___1397], 0.0, (FieldDataType_e)___1364[___1397])) { WriteErr(___1397, "TECDAT142"); return (-1); } } if (IIndex + 1 == FinalIMax && (JIndex + 1 == FinalJMax && FinalJMax < ___2114[___1397] - JMaxAdjust)) { ___2227 II; for (II = 1; II <= ___1906[___1397] - IMaxAdjust; II++) { NumOrderedCCDataValuesWritten[___1397]++; if (!___4493(BlckFile[___1397], 0.0, (FieldDataType_e)___1364[___1397])) { WriteErr(___1397, "TECDAT142"); return (-1); } } } if (IIndex + 1 == FinalIMax && JIndex + 1 == FinalJMax && (KIndex + 1 == FinalKMax && FinalKMax < ___2159[___1397] - KMaxAdjust)) { ___2227 JJ, II; for (JJ = 1; JJ <= ___2114[___1397] - JMaxAdjust; JJ++) for (II = 1; II <= ___1906[___1397] - IMaxAdjust; II++) { NumOrderedCCDataValuesWritten[___1397]++; if (!___4493(BlckFile[___1397], 0.0, (FieldDataType_e)___1364[___1397])) { WriteErr(___1397, "TECDAT142"); return (-1); } } } NumOrderedCCDataValuesWritten[___1397]++; } NumDataValuesWritten[___1397]++; if (___2006[___1397]) { if (NumRunningVarValues[___1397][___718[___1397]] == NumDataValuesWritten[___1397]) { AdvanceToNextVarWithValues(___1397); if (___718[___1397] < NumVars[___1397]       && IsCellCentered[___1397][___718[___1397]] && ___4693[___1397] == ___2882) NumOrderedCCDataValuesWritten[___1397] = 0; } } else { AdvanceToNextVarWithValues(___1397); if (___718[___1397] >= NumVars[___1397]) { ___718[___1397] = -1; AdvanceToNextVarWithValues(___1397); } }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397] > 1) ___3186("%lld %G\n", (long long)NumDataValuesWritten[___1397], ___4315);
 #endif
} if (HasFECONNECT[___1397] && (NumDataValuesToWrite[___1397] == NumDataValuesWritten[___1397])) { if (!___4490(BlckFile[___1397], (___2227)1)) { WriteErr(___1397, "TECDAT142"); return (-1); } } if (NumDataValuesToWrite[___1397] == NumDataValuesWritten[___1397]) RewritePendingMinMaxValues(___1397); return (0); } int32_t ___3979( int32_t        ___1397, int32_t const* ___2689) { ___2227 ___2165 = NumConnectivityNodes[___1397][___734[___1397]]; ConnectivityWritten[___1397][___734[___1397]] = ___4226; if (CheckFile(___1397, "TECNOD142") < 0) return (-1); if (___4693[___1397] == ___1341 || ___4693[___1397] == ___1342) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNOD142) Cannot call TECNOD142 for polygonal or polyhedral zones.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (HasFECONNECT[___1397]) { return (-1); } if (FileTypes[___1397] == ___3638) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNOD142) Cannot call TECNOD142 if file type is SOLUTIONFILE.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (___4693[___1397] == ___2882) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNOD142) Cannot call TECNOD142 if zone type is ORDERED.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (CheckData(___1397, "TECNOD142") < 0) return (-1); for (___2227 ___1832 = 0; ___1832 < ___2165; ___1832++) { if ((___2689[___1832] > ___1906[___1397]) || (___2689[___1832] < 1)) {
 #if defined MAKEARCHIVE
___3185("Err: (TECNOD142) Invalid node map value at position %ld:\n", ___1832); ___3186("     node map value = %lld, max value = %lld.\n", (long long)___2689[___1832], (long long)___1906[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } if (!___4490(BlckFile[___1397], ___2689[___1832] - 1)) { WriteErr(___1397, "TECNOD142"); return (-1); } } return (0); } int32_t ___3980( int32_t        ___1397, int32_t const* N, int32_t const* ___2689) { if (CheckFile(___1397, "TECNODE142") < 0) return (-1); if (___4693[___1397] == ___1341 || ___4693[___1397] == ___1342) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNODE142) Cannot call TECNODE142 for polygonal or polyhedral zones.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (HasFECONNECT[___1397]) { return (-1); } if (FileTypes[___1397] == ___3638) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNODE142) Cannot call TECNODE142 if file type is SOLUTIONFILE.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (___4693[___1397] == ___2882) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNODE142) Cannot call TECNODE142 if zone type is ORDERED.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (CheckData(___1397, "TECNODE142") < 0) return (-1); if (NumConnectivityNodesWritten[___1397][___734[___1397]] + *N > NumConnectivityNodes[___1397][___734[___1397]]) {
 #if defined MAKEARCHIVE
___3184("Err: (TECNODE142) Connectivity Nodes chunk exceeds the total number of  Connectivity Nodes:\n"); ___3186("     Nodes written = %ld, Current chunk size = %d, ", NumConnectivityNodesWritten[___1397][___734[___1397]], *N); ___3185("total connectivity nodes = %ld.\n", NumConnectivityNodes[___1397][___734[___1397]]);
 #endif
NumErrs[___1397]++; return (-1); } for (___2227 ___1832 = 0; ___1832 < *N; ___1832++) { if ((___2689[___1832] > ___1906[___1397]) || (___2689[___1832] < 1)) {
 #if defined MAKEARCHIVE
___3185("Err: (TECNODE142) Invalid node map value at position %lld:\n", (long long)___1832); ___3186("     node map value = %lld, max value = %lld.\n", (long long)___2689[___1832], (long long)___1906[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } if (!___4490(BlckFile[___1397], ___2689[___1832] - 1)) { WriteErr(___1397, "TECNODE142"); return (-1); } NumConnectivityNodesWritten[___1397][___734[___1397]]++; } if (NumConnectivityNodesWritten[___1397][___734[___1397]] == NumConnectivityNodes[___1397][___734[___1397]]) { ConnectivityWritten[___1397][___734[___1397]] = ___4226; } return (0); } int32_t ___3973(int32_t ___1397) { int RetVal = 0; if (FileTypes[___1397] != ___3638) { for (int ZoneIndex = 0; (RetVal == 0) && (ZoneIndex <= ___734[___1397]); ZoneIndex++) { if (((NumConnectivityNodes[___1397][ZoneIndex] > 0) && (ConnectivityWritten[___1397][ZoneIndex] == ___1305))) {
 #if defined MAKEARCHIVE
___3185("Err: (TECEND142) File %d is being closed without writing connectivity data.\n", ___1397 + 1); if (NumConnectivityNodesWritten[___1397][ZoneIndex] == 0) ___3185("     Zone %d was defined with a Classic FE zone type but TECNOD142() was not called.\n", ZoneIndex + 1); else ___3185("     Zone %d was defined with a Classic FE zone type but TECNODE142() was not called for all node chunks.\n", ZoneIndex + 1);
 #endif
NumErrs[___1397]++; RetVal = -1; } if (((___2803[___1397][ZoneIndex] > 0) && (FaceNeighborsOrMapWritten[___1397][ZoneIndex] == ___1305))) {
 #if defined MAKEARCHIVE
___3185("Err: (TECEND142) File %d is being closed without writing face neighbor data.\n", ___1397 + 1); ___3186("     %d connections were specified for zone %d but TECFACE142() was not called.\n", ___2803[___1397][ZoneIndex], ZoneIndex + 1);
 #endif
NumErrs[___1397]++; RetVal = -1; } else if (((___4193[___1397][ZoneIndex] > 0) && (FaceNeighborsOrMapWritten[___1397][ZoneIndex] == ___1305))) {
 #if defined MAKEARCHIVE
___3185("Err: (TECEND142) File %d is being closed without writing face map data.\n", ___1397 + 1); ___3186("     %lld face nodes were specified for zone %d but TECPOLYFACE142() and/or\n" "     TECPOLYBCONN142() was not called for all face chunks.\n", (long long)___4193[___1397][ZoneIndex], ZoneIndex + 1);
 #endif
NumErrs[___1397]++; RetVal = -1; } } } if (RetVal == 0) { if (CheckFile(___1397, "TECEND142") < 0) RetVal = -1; } if (RetVal == 0) { if (CheckData(___1397, "TECEND142") < 0) RetVal = -1; } if (RetVal == 0) if (!___4493(HeadFile[___1397], EndHeaderMarker, FieldDataType_Float)) { WriteErr(___1397, "TECEND142"); RetVal = -1; } CloseFileStream(&BlckFile[___1397]); if (RetVal == 0) { BlckFile[___1397] = OpenFileStream(BlckFName[___1397], "rb", IsWritingNative[___1397]); char ___416[BYTES_PER_CHUNK]; size_t bytesRead = 0; while ((RetVal == 0) && (feof(BlckFile[___1397]->File) == 0)) { bytesRead = fread((void*)___416, 1, BYTES_PER_CHUNK, BlckFile[___1397]->File); if (ferror(BlckFile[___1397]->File) == 0) { if (bytesRead != fwrite((void*)___416, 1, bytesRead, HeadFile[___1397]->File)) {
 #if defined MAKEARCHIVE
___3185("Err: (TECEND142) Write failure during repack on file %d.\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; RetVal = -1; } } else {
 #if defined MAKEARCHIVE
___3185("Err: (TECEND142) Write failure during repack on file %d.\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; RetVal = -1; } } CloseFileStream(&BlckFile[___1397]); } ___4206(BlckFName[___1397]); CloseFileStream(&HeadFile[___1397]);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) { ___3185("File %d closed.\n", ___1397 + 1); if (NumErrs[___1397]) { ___3184("********************************************\n"); ___3185("      %d Errors occurred on this file\n", NumErrs[___1397]); ___3184("********************************************\n"); } }
 #endif
NumErrs[___1397] = 0; ___2042[___1397] = 0; if (DestFName[___1397]) ___1530(DestFName[___1397], "data set fname"); if (BlckFName[___1397]) ___1530(BlckFName[___1397], "data set fname"); BlckFName[___1397] = NULL; DestFName[___1397] = NULL; VarMinValue[___1397].clear(); VarMaxValue[___1397].clear(); NumRunningVarValues[___1397].clear(); IsSharedVar[___1397].clear(); IsPassiveVar[___1397].clear(); IsCellCentered[___1397].clear(); MinMaxOffset[___1397].clear(); ___4193[___1397].clear(); PolyZoneWriteInfos[___1397].clear(); ___2803[___1397].clear(); FaceNeighborsOrMapWritten[___1397].clear(); NumConnectivityNodes[___1397].clear(); NumConnectivityNodesWritten[___1397].clear(); ConnectivityWritten[___1397].clear(); return RetVal; } namespace { void GetNextLabel( char const** ___685, char*        NextLabel) { int N = 0; char *NPtr = NextLabel; *NPtr = '\0'; while ((**___685) && (**___685 !='"')) (*___685)++; if (**___685) (*___685)++; while ((N < 60) && (**___685) && (**___685 != '"')) { if (**___685 == '\\') { (*___685)++; } *NPtr++ = **___685; N++; (*___685)++; } if (**___685) (*___685)++; *NPtr = '\0'; } } int32_t ___3978( int32_t    ___1397, char const* S) { char const* ___685 = S; ___2227   N = 0; char        Label[60]; if (CheckFile(___1397, "TECLAB142") < 0) return (-1);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3184("\nInserting Custom Labels:\n");
 #endif
do { GetNextLabel(&___685, Label); if (*Label) N++; } while (*Label); if (N == 0) {
 #if defined MAKEARCHIVE
___3185("Err: (TECLAB142) Invalid custom label string: %s\n", (S ? S : " "));
 #endif
NumErrs[___1397]++; return (-1); } ___4493(HeadFile[___1397], ___792, FieldDataType_Float); if (!___4490(HeadFile[___1397], (___2227)N)) { WriteErr(___1397, "TECLAB142"); return (-1); } ___685 = S; do { GetNextLabel(&___685, Label); if (*Label) { if (!___1131(HeadFile[___1397], Label, ___4226)) { WriteErr(___1397, "TECLAB142"); return (-1); }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) printf("          %s\n", Label);
 #endif
} } while (*Label); return (0); } int32_t ___3989( int32_t    ___1397, char const* S) { if (CheckFile(___1397, "TECUSR142") < 0) return (-1);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3185("\nInserting UserRec: %s\n", S);
 #endif
if ((S == NULL) || (*S == '\0')) {
 #if defined MAKEARCHIVE
___3184("Err: (TECUSR142) Invalid TECUSR142 string\n");
 #endif
NumErrs[___1397]++; return (-1); } ___4493(HeadFile[___1397], ___4286, FieldDataType_Float); if (!___1131(HeadFile[___1397], S, ___4226)) {
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) printf("Err: (TECUSR142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } return (0); } int32_t ___3976( int32_t        ___1397, double const*   ___4575, double const*   ___4592, double const*   ___4716, int32_t const* ___3160, int32_t const* ___227, int32_t const* ___4600, int32_t const* Color, int32_t const* ___1412, int32_t const* ___2023, int32_t const* ___1652, int32_t const* ___2264, double const*   ___2987, double const*   ___2290, int32_t const* ___2794, int32_t const* ___188, int32_t const* ___176, double const*   ___187, double const*   ___171, int32_t const* ___3443, int32_t const* ___496, int32_t const* ___2836, int32_t const* ___2838, float const*    ___4573, float const*    ___4590, float const*    ___4597, char const*     mfc) { int    ___1832, RetVal; int    RawDataSize = 0; double Fract; ___1632 ___1556; if (CheckFile(___1397, "TECGEO142") < 0) return (-1); ___1556.___3167 = (CoordSys_e) * ___3160; if (___1556.___3167 == CoordSys_Frame) Fract = 0.01; else Fract = 1.0; ___1556.position.setXOrTheta((*___4575) * Fract); ___1556.position.setYOrR((*___4592) * Fract); ___1556.position.setZ((*___4716) * Fract); ___1556.___227           = *___227 != 0; ___1556.___4600                   = *___4600 - 1; ___1556.___351                 = (___516) * Color; ___1556.___1410             = (___516) * ___1412; ___1556.___2023               = *___2023 != 0; ___1556.___1652               = (GeomType_e) * ___1652; ___1556.___2264            = (LinePattern_e) * ___2264; ___1556.___2987          = *___2987 / 100.0; ___1556.___2290          = *___2290 / 100.0; ___1556.___2794          = (int32_t)*___2794; ___1556.___188         = (ArrowheadStyle_e) * ___188; ___1556.___176    = (ArrowheadAttachment_e) * ___176; ___1556.___187          = *___187 / 100.0; ___1556.___171         = *___171 / ___954; ___1556.___3443                  = (Scope_e) * ___3443; ___1556.___1113              = ___1114; ___1556.___496               = (Clipping_e) * ___496; ___1556.___2836            = (int32_t)*___2836; ___1556.___2331   = const_cast<char*>(mfc); ___1556.___1884          = NULL; ___1556.WorldFileName          = NULL; ___1556.EmbeddedLpkImageNumber = 0; ___1556.___2333    = ___4226; ___1556.___3089       = 1.0; ___1556.___1890      = ___1900; if (___1556.___1652 == GeomType_LineSegs3D) { ___1556.___1652         = GeomType_LineSegs; ___1556.___3167 = CoordSys_Grid3D; }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3184("\nInserting Geometry\n");
 #endif
switch (___1556.___1652) { case GeomType_LineSegs : { int ___1832; RawDataSize = 0; for (___1832 = 0; ___1832 < *___2836; ___1832++) { ___1556.___2838[___1832] = ___2838[___1832]; RawDataSize += ___2838[___1832]; } } break; case GeomType_Rectangle : case GeomType_Square : case GeomType_Circle : case GeomType_Ellipse : { RawDataSize = 1; } break; case GeomType_Image : { ___478(___1305); } break; default : { ___478(___1305); } break; } ___1556.___907                = FieldDataType_Float; ___1556.___1573.___1548.___4293 = ___28(RawDataSize, FieldDataType_Float, ___4226); ___1556.___1573.___1548.___4295 = ___28(RawDataSize, FieldDataType_Float, ___4226); ___1556.___1573.___1548.___4297 = ___28(RawDataSize, FieldDataType_Float, ___4226); for (___1832 = 0; ___1832 < RawDataSize; ___1832++) { ___3490(___1556.___1573.___1548.___4293, ___1832, (double)___4573[___1832]*Fract); ___3490(___1556.___1573.___1548.___4295, ___1832, (double)___4590[___1832]*Fract); ___3490(___1556.___1573.___1548.___4297, ___1832, (double)___4597[___1832]*Fract); } if (___1132(HeadFile[___1397], &___1556, ___4226, ___1305)) RetVal = 0; else RetVal = -1; ___938(&___1556.___1573.___1548.___4293); ___938(&___1556.___1573.___1548.___4295); ___938(&___1556.___1573.___1548.___4297); return RetVal; } int32_t ___3988( int32_t        ___1397, double const*   ___4575, double const*   ___4592, double const*   ___4714, int32_t const* ___3160, int32_t const* ___227, int32_t const* ___4600, int32_t const* ___353, int32_t const* ___1453, double const*   ___1451, int32_t const* ___411, double const*   ___409, double const*   ___407, int32_t const* ___403, int32_t const* ___405, double const*   ___57, int32_t const* ___39, double const*   ___2288, int32_t const* ___4081, int32_t const* ___3443, int32_t const* ___496, char const*     ___3813, char const*     mfc) { int    RetVal; ___4118 Text; double Fract; if (CheckFile(___1397, "TECTXT142") < 0) return (-1); Text.___3167    = (CoordSys_e) * ___3160; if (Text.___3167 == CoordSys_Frame) Fract = 0.01; else Fract = 1.0; Text.___52.___1548.___4292 = (*___4575) * Fract; Text.___52.___1548.___4294 = (*___4592) * Fract; Text.___52.___1548.___4296 = (*___4714) * Fract; Text.___227         = *___227 != 0; Text.___4600                 = *___4600 - 1; Text.___351               = static_cast<___516>(*___4081); Text.___4121.___1444       = static_cast<Font_e>(*___353); Text.___4121.___3601  = static_cast<Units_e>(*___1453); if (Text.___4121.___3601 == ___4268) Text.___4121.___1827   = (*___1451) / 100.0; else Text.___4121.___1827   = *___1451; Text.___401.___411          = static_cast<TextBox_e>(*___411); Text.___401.___2338           = *___409 / 100.0; Text.___401.___2290    = *___407 / 100.0; Text.___401.___351           = static_cast<___516>(*___403); Text.___401.___1410       = static_cast<___516>(*___405); Text.___39               = static_cast<TextAnchor_e>(*___39); Text.___2288          = *___2288; Text.___57                = *___57 / ___954; Text.___3443                = static_cast<Scope_e>(*___3443); Text.Text                 = const_cast<char*>(___3813); Text.___2331 = const_cast<char*>(mfc); Text.___496             = static_cast<Clipping_e>(*___496);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3185("\nInserting Text: %s\n", ___3813);
 #endif
if (___1133(HeadFile[___1397], &Text, ___4226, ___1305)) RetVal = 0; else RetVal = -1; return RetVal; } void ___3975( int32_t        ___1397, int32_t const* ___2891) { ___4278(___1397); REQUIRE(VALID_REF(___2891)); DoWriteForeign = (*___2891 != 0); }
 #if defined MAKEARCHIVE
namespace { ___372 ___250( char      ___472, ___372 ___2030) { REQUIRE(0 <= ___472 && "Char <= 127"); REQUIRE(VALID_BOOLEAN(___2030)); ___372 IsValidNameChar = (___472 == '_' || tecplot::___1998(___472)); if (!___2030) IsValidNameChar = (IsValidNameChar || ___472 == '.'     || tecplot::___2012(___472)); ENSURE(VALID_BOOLEAN(IsValidNameChar)); return IsValidNameChar; } }
 #endif 
 #if defined MAKEARCHIVE
namespace { ___372 ___249(char const* ___2686) { REQUIRE(VALID_REF(___2686)); ___372 IsValidName = ___250(*___2686, ___4226); for (char const* NPtr = ___2686; IsValidName && *NPtr != '\0'; ++NPtr) { IsValidName = ___250(*NPtr, ___1305); } ENSURE(VALID_BOOLEAN(IsValidName)); return IsValidName; } }
 #endif 
int32_t ___3971( int32_t    ___1397, char const* ___2686, char const* ___4315) { if (CheckFile(___1397, "TECAUXSTR142") < 0) return (-1);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3186("\nInserting data set aux data: '%s' = '%s'\n", ___2686, ___4315);
 #endif
if ((___2686 == NULL) || !___249(___2686)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECAUXSTR142) Invalid Name string\n");
 #endif
NumErrs[___1397]++; return (-1); } if ((___4315 == NULL) || (*___4315 == '\0')) {
 #if defined MAKEARCHIVE
___3184("Err: (TECAUXSTR142) Invalid Value string\n");
 #endif
NumErrs[___1397]++; return (-1); } if (!___4493(HeadFile[___1397], ___887, FieldDataType_Float)  || !___1131(HeadFile[___1397], ___2686, ___4226  ) || !___4490(HeadFile[___1397], (___2227)___270) || !___1131(HeadFile[___1397], ___4315, ___4226  )) {
 #if defined MAKEARCHIVE
printf("Err: (TECAUXSTR142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } return (0); } int32_t ___3991( int32_t    ___1397, char const* ___2686, char const* ___4315) { if (CheckFile(___1397, "TECZAUXSTR142") < 0) return (-1); if (___734[___1397] == -1) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZAUXSTR142) Must call TECZNE142 prior to TECZAUXSTR142\n");
 #endif
NumErrs[___1397]++; return (-1); }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3186("\nInserting zone aux data: '%s' = '%s'\n", ___2686, ___4315);
 #endif
if ((___2686 == NULL) || !___249(___2686)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZAUXSTR142) Invalid Name string\n");
 #endif
NumErrs[___1397]++; return (-1); } if ((___4315 == NULL) || (*___4315 == '\0')) {
 #if defined MAKEARCHIVE
___3184("Err: (TECZAUXSTR142) Invalid Value string\n");
 #endif
NumErrs[___1397]++; return (-1); } if (___4200(HeadFile[___1397]->File, -4, SEEK_CUR) || !___4490(HeadFile[___1397], 1)  || !___1131(HeadFile[___1397], ___2686, ___4226  ) || !___4490(HeadFile[___1397], (___2227)___270) || !___1131(HeadFile[___1397], ___4315, ___4226  ) || !___4490(HeadFile[___1397], 0)) {
 #if defined MAKEARCHIVE
printf("Err: (TECZAUXSTR142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } return (0); } int32_t ___3990( int32_t        ___1397, int32_t const* ___4337, char const*     ___2686, char const*     ___4315) { if (CheckFile(___1397, "TECVAUXSTR142") < 0) return (-1);
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3186("\nInserting variable aux data: '%s' = '%s'\n", ___2686, ___4315);
 #endif
if ((___2686 == NULL) || !___249(___2686)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECVAUXSTR142) Invalid Name string\n");
 #endif
NumErrs[___1397]++; return (-1); } if ((___4315 == NULL) || (*___4315 == '\0')) {
 #if defined MAKEARCHIVE
___3184("Err: (TECVAUXSTR142) Invalid Value string\n");
 #endif
NumErrs[___1397]++; return (-1); } if (!___4493(HeadFile[___1397], ___4339, FieldDataType_Float)  || !___4490(HeadFile[___1397], *___4337 - 1) || !___1131(HeadFile[___1397], ___2686, ___4226  ) || !___4490(HeadFile[___1397], (___2227)___270) || !___1131(HeadFile[___1397], ___4315, ___4226  )) {
 #if defined MAKEARCHIVE
printf("Err: (TECVAUXSTR142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } return (0); } int32_t ___3974( int32_t        ___1397, int32_t const* ___1258) { FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; if (CheckFile(___1397, "TECFACE142") < 0) return (-1); if (___4693[___1397] == ___1341 || ___4693[___1397] == ___1342) {
 #if defined MAKEARCHIVE
___3184("Err: (TECFACE142) Cannot call TECFACE142 for polygonal or polyhedral zones.\n");
 #endif
NumErrs[___1397]++; return (-1); } if (FileTypes[___1397] == ___3638) {
 #if defined MAKEARCHIVE
___3184("Err: (TECFACE142) Cannot call TECFACE142 if the file type is SOLUTIONFILE.\n");
 #endif
NumErrs[___1397]++; return (-1); }
 #if defined MAKEARCHIVE
if (DebugLevel[___1397]) ___3184("\nInserting face neighbor data\n");
 #endif
if (___1258 == NULL) {
 #if defined MAKEARCHIVE
___3184("Err: (TECFACE142) Invalid array\n");
 #endif
NumErrs[___1397]++; return (-1); } int32_t const* Ptr = ___1258; int32_t i = 0; while (i < ___2803[___1397][___734[___1397]]) { int32_t n; int32_t NumNum = 0; switch (___1285[___1397]) { case ___1290: NumNum = 3; i++; break; case ___1289: NumNum = 4 + Ptr[3]; i += Ptr[3]; break; case ___1287: NumNum = 4; i++; break; case ___1286: NumNum = 4 + 2 * Ptr[3]; i += Ptr[3]; break; default: ___478(___1305); break; } n = 0; if (___1285[___1397] == ___1289 || ___1285[___1397] == ___1286) { if (!___4490(BlckFile[___1397], Ptr[n++] - 1) || !___4490(BlckFile[___1397], Ptr[n++] - 1) || !___4490(BlckFile[___1397], Ptr[n++])   || !___4490(BlckFile[___1397], Ptr[n++])) {
 #if defined MAKEARCHIVE
printf("Err: (TECFACE142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } } for (; n < NumNum; n++) if (!___4490(BlckFile[___1397], Ptr[n] - 1)) {
 #if defined MAKEARCHIVE
printf("Err: (TECFACE142) Write failure for file %d\n", ___1397 + 1);
 #endif
NumErrs[___1397]++; return (-1); } Ptr += NumNum; } return (0); } int32_t ___3986(int32_t ___1397) { return ___4190[___1397]; } int32_t ___3987(int32_t ___1397) { return ___4191[___1397]; } int32_t ___3983( int32_t        ___1397, int32_t const* ___1294, int32_t const* ___1297, int32_t const* ___1259, int32_t const* ___1303, int32_t const* ___1253, int32_t const* ___1254, int32_t const* ___1256) { int32_t ___2806 = ___2159[___1397]; int32_t ___3359 = 0; ___2227 ___1926; ___2227 MinNeighborValue = TECIO_NO_NEIGHBORING_ELEM; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; if (___2806 == 0 || (___4693[___1397] != ___1341 && ___4693[___1397] != ___1342)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLY142) The zone type must be FEPOLYGON or FEPOLYHEDRON and have NumFaces (KMax) > 0.\n"); ___3185("     NumFaces = %d\n", ___2806);
 #endif
NumErrs[___1397]++; return (-1); } if (___4693[___1397] == ___1342) { if (___4193[___1397][___734[___1397]] <= 0) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLY142) TotalNumFaceNodes MUST be specified for polyhedral zones.\n"); ___3185("     TotalNumFaceNodes = %lld\n", (long long)___4193[___1397][___734[___1397]]);
 #endif
NumErrs[___1397]++; return (-1); } } else { if (___4193[___1397][___734[___1397]] != (2 * ___2806)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLY142) TotalNumFaceNodes is specified for the polygonal zone but is not equal to 2 * NumFaces.\n"); ___3186("     TotalNumFaceNodes = %lld.  If specified, it must be 2 * %d.", (long long)___4193[___1397][___734[___1397]], ___2806);
 #endif
NumErrs[___1397]++; return (-1); } } if ((___4191[___1397] > 0  && ___4190[___1397] > 0) || (___4191[___1397] == 0 && ___4190[___1397] == 0)) { if (___4191[___1397] > 0) MinNeighborValue = -___4191[___1397]; } else {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLY142) TotalNumFaceBndryFaces and TotalNumFaceBndryConns must both be 0 or both be > 0.\n"); ___3186("     TotalNumFaceBndryFaces = %lld, TotalNumFaceBndryConns = %lld\n", (long long)___4191[___1397], (long long)___4190[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } if (___3359 == 0) { if (___4693[___1397] == ___1342) { int32_t FaceNodeSum = 0; if (!___4490(BlckFile[___1397], 0)) ___3359 = -1; for (___1926 = 0; (___3359 == 0) && (___1926 < ___2806); ___1926++) { FaceNodeSum += ___1294[___1926]; if (___1294[___1926] < 3) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid face node count value at face %lld.  There must be at least 3 nodes in a face.\n", (long long)___1926 + 1); ___3185("     Face node count value = %lld.\n", (long long)___1294[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (FaceNodeSum > ___4193[___1397][___734[___1397]]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) The running face node count exceeds the TotalNumFaceNodes (%lld) specified.\n", (long long)___4193[___1397][___734[___1397]]); ___3185("     Face node count value = %lld.\n", (long long)___1294[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], FaceNodeSum)) ___3359 = -1; } } } for (___1926 = 0; (___3359 == 0) && (___1926 < ___4193[___1397][___734[___1397]]); ___1926++) { if (___1297[___1926] < 1 || ___1297[___1926] > ___1906[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid face node value at node %lld:\n", (long long)___1926 + 1); ___3186("     face node value = %lld, valid values are are 1 to %lld (inclusive).\n", (long long)___1297[___1926], (long long)___1906[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1297[___1926] - 1)) ___3359 = -1; } for (___1926 = 0; (___3359 == 0) && (___1926 < ___2806); ___1926++) { if (___1259[___1926] < MinNeighborValue || ___1259[___1926] > ___2114[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid left neighbor value at face %lld:\n", (long long)___1926); ___3186("     left neighbor value = %lld, min value = %lld,", (long long)___1259[___1926], (long long)MinNeighborValue); ___3185(" max value = %lld.\n", (long long)___2114[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1259[___1926] - 1)) ___3359 = -1; } for (___1926 = 0; (___3359 == 0) && (___1926 < ___2806); ___1926++) { if (___1303[___1926] < MinNeighborValue || ___1303[___1926] > ___2114[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid right neighbor value at face %lld:\n", (long long)___1926); ___3186("     right neighbor value = %lld, min value = %lld,", (long long)___1303[___1926], (long long)MinNeighborValue); ___3185(" max value = %lld.\n", (long long)___2114[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1303[___1926] - 1)) ___3359 = -1; if (___3359 == 0 && (___1259[___1926] == TECIO_NO_NEIGHBORING_ELEM && ___1303[___1926] == TECIO_NO_NEIGHBORING_ELEM)) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Both left and right neighbors are set to no neighboring element at face %lld.\n", (long long)___1926);
 #endif
NumErrs[___1397]++; return (-1); } } if (___3359 == 0 && ___4191[___1397] > 0) { ___478(TecplotSDKBinaryFileVersion == 112); if (!(___4490(BlckFile[___1397], 0) && ___4490(BlckFile[___1397], 0))) ___3359 = -1; int32_t BndryConnCount = 0; for (___1926 = 0; (___3359 == 0) && (___1926 < ___4191[___1397]); ___1926++) { BndryConnCount += ___1253[___1926]; if (___1253[___1926] < 0 || BndryConnCount > ___4190[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid boundary connection count at boundary face %lld:\n", (long long)___1926 + 1); ___3185("     boundary connection count = %lld.\n", (long long)___1253[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], BndryConnCount)) ___3359 = -1; } if (BndryConnCount != ___4190[___1397]) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLY142) Invalid number of boundary connections:\n"); ___3186("     number of boundary connections written = %lld, total number of boundary connections = %lld.", (long long)BndryConnCount, (long long)___4190[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } BndryConnCount = 0; for (___1926 = 0; (___3359 == 0) && (___1926 < ___4191[___1397]); ___1926++) { for (___2227 BIndex = 0; (___3359 == 0) && (BIndex < ___1253[___1926]); BIndex++) { if (BIndex > 0 && ___1254[BndryConnCount] == TECIO_NO_NEIGHBORING_ELEM) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Partially obscured faces must specify no neighboring element first. See boundary connections for face %lld.\n", (long long)___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } if (___1254[BndryConnCount] < TECIO_NO_NEIGHBORING_ELEM) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid boundary element value at boundary connections for face %lld:\n", (long long)___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } if (___1254[BndryConnCount] == TECIO_NO_NEIGHBORING_ELEM && ___1256[BndryConnCount] != TECIO_NO_NEIGHBORING_ZONE) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid boundary element/zone pair at boundary connections for face %lld:\n", (long long)___1926 + 1); ___3184("     Boundary elements specified as no neighboring element must also specify no neighboring zone.\n");
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1254[BndryConnCount] - 1)) ___3359 = -1; BndryConnCount++; } } BndryConnCount = 0; for (___1926 = 0; (___3359 == 0) && (___1926 < ___4191[___1397]); ___1926++) { for (___2227 BIndex = 0; (___3359 == 0) && (BIndex < ___1253[___1926]); BIndex++) { if (___1256[BndryConnCount] < TECIO_NO_NEIGHBORING_ZONE) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid boundary zone value at boundary connections for face %lld:\n", (long long)___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1256[BndryConnCount] - 1)) ___3359 = -1; BndryConnCount++; } } } if (___3359 != 0) { ___3359 = -1; WriteErr(___1397, "TECPOLY142"); } return ___3359; } namespace { void calculateInitialPolyZoneOffsets(int32_t ___1397) { PolyZoneWriteInfo& polyZoneWriteInfo = PolyZoneWriteInfos[___1397][___734[___1397]]; polyZoneWriteInfo.numFaceNodesOffset = (___1398)___4201(BlckFile[___1397]->File); if (___4693[___1397] == ___1342) { polyZoneWriteInfo.faceNodesOffset = polyZoneWriteInfo.numFaceNodesOffset + sizeof(int32_t) * (___2159[___1397] + 1); } else { polyZoneWriteInfo.faceNodesOffset = polyZoneWriteInfo.numFaceNodesOffset; } polyZoneWriteInfo.leftElemsOffset = polyZoneWriteInfo.faceNodesOffset + sizeof(int32_t) * ___4193[___1397][___734[___1397]]; polyZoneWriteInfo.rightElemsOffset = polyZoneWriteInfo.leftElemsOffset + sizeof(int32_t) * ___2159[___1397]; polyZoneWriteInfo.connectionCountsOffset = polyZoneWriteInfo.rightElemsOffset + sizeof(int32_t) * ___2159[___1397]; polyZoneWriteInfo.connectionElemsOffset = polyZoneWriteInfo.connectionCountsOffset + sizeof(int32_t) * (___4191[___1397] + 2); polyZoneWriteInfo.connectionZonesOffset = polyZoneWriteInfo.connectionElemsOffset + sizeof(int32_t) * ___4190[___1397]; } } namespace { int32_t checkForPolyFacePreconditions(int32_t ___1397) { if (CheckFile(___1397, "TECPOLYFACE142") < 0) return (-1); if (FaceNeighborsOrMapWritten[___1397][___734[___1397]]) { if (NumErrs[___1397] == 0) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYFACE142) TECPOLYFACE142 called after all specified faces were already written.\n");
 #endif
NumErrs[___1397]++; } return (-1); } if (___2159[___1397] == 0 || (___4693[___1397] != ___1341 && ___4693[___1397] != ___1342)) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYFACE142) The zone type must be FEPOLYGON or FEPOLYHEDRON and have KMax (NumFaces) > 0.\n"); ___3185("     KMax = %lld\n", (long long)___2159[___1397]);
 #endif
NumErrs[___1397]++; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; return (-1); } if (___4693[___1397] == ___1342) { if (___4193[___1397][___734[___1397]] <= 0) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYFACE142) TotalNumFaceNodes MUST be specified for polyhedral zones.\n"); ___3185("     TotalNumFaceNodes = %lld\n", (long long)___4193[___1397][___734[___1397]]);
 #endif
NumErrs[___1397]++; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; return (-1); } } else { if (___4193[___1397][___734[___1397]] != (2 * ___2159[___1397])) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYFACE142) TotalNumFaceNodes is specified for the polygonal zone but is not equal to 2 * KMax.\n"); ___3186("     TotalNumFaceNodes = %lld.  If specified, it must be 2 * %lld.", (long long)___4193[___1397][___734[___1397]], (long long)___2159[___1397]);
 #endif
NumErrs[___1397]++; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; return (-1); } } return 0; } } namespace { int32_t getMinNeighborValue( int32_t   ___1397, int32_t&   MinNeighborValue) { MinNeighborValue = TECIO_NO_NEIGHBORING_ELEM; if ((___4191[___1397] > 0  && ___4190[___1397] > 0) || (___4191[___1397] == 0 && ___4190[___1397] == 0)) { if (___4191[___1397] > 0) MinNeighborValue = -___4191[___1397]; } else {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYFACE142) TotalNumFaceBndryFaces and TotalNumFaceBndryConns must both be 0 or both be > 0.\n"); ___3186("     TotalNumFaceBndryFaces = %lld, TotalNumFaceBndryConns = %lld\n", (long long)___4191[___1397], (long long)___4190[___1397]);
 #endif
NumErrs[___1397]++; FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; return (-1); } return 0; } } namespace { int32_t writePolyhedralFaceNodeOffsets( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, int32_t           ___2806, int32_t const*    ___1294) { { ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.numFaceNodesOffset, SEEK_SET); if (polyZoneWriteInfo.faceNodeSum == 0) if (!___4490(BlckFile[___1397], 0)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } for (___2227 ___1926 = 0; ___1926 < ___2806; ___1926++) { polyZoneWriteInfo.faceNodeSum += ___1294[___1926]; if (___1294[___1926] < 3) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid face node count value at face %lld.  There must be at least 3 nodes in a face.\n", (long long)___1926 + 1); ___3185("     Face node count value = %lld.\n", (long long)___1294[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (polyZoneWriteInfo.faceNodeSum > ___4193[___1397][___734[___1397]]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) The running face node count exceeds the TotalNumFaceNodes (%lld) specified.\n", (long long)___4193[___1397][___734[___1397]]); ___3185("     Face node count value = %lld.\n", (long long)___1294[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], polyZoneWriteInfo.faceNodeSum)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } } polyZoneWriteInfo.numFaceNodesOffset = (___1398)___4201(BlckFile[___1397]->File); } return 0; } } namespace { int32_t writeFaceNodes( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, ___2227          BeginningFaceNodeSum, int32_t           ___2806, int32_t const*    ___1297) { ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.faceNodesOffset, SEEK_SET); ___2227 LocalFaceNodeSum; if (___4693[___1397] == ___1342) LocalFaceNodeSum = polyZoneWriteInfo.faceNodeSum - BeginningFaceNodeSum; else LocalFaceNodeSum = 2 * (___2806); for (___2227 ___1926 = 0; ___1926 < LocalFaceNodeSum; ___1926++) { if (___1297[___1926] < 1 || ___1297[___1926] > ___1906[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLY142) Invalid face node value at node %lld:\n", (long long)___1926 + 1); ___3186("     face node value = %lld, valid values are are 1 to %lld (inclusive).\n", (long long)___1297[___1926], (long long)___1906[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1297[___1926] - 1)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } } polyZoneWriteInfo.faceNodesOffset = (___1398)___4201(BlckFile[___1397]->File); return 0; } } namespace { int32_t checkElements( int32_t        ___1397, int32_t        ___2806, int32_t const* ___1259, int32_t const* ___1303, int32_t        MinNeighborValue) { for (___2227 ___1926 = 0; ___1926 < ___2806; ___1926++) { if (___1259[___1926] < MinNeighborValue || ___1259[___1926] > ___2114[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYFACE142) Invalid left neighbor value at face %lld:\n", (long long)___1926); ___3186("     left neighbor value = %lld, min value = %lld,", (long long)___1259[___1926], (long long)MinNeighborValue); ___3185(" max value = %lld.\n", (long long)___2114[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (___1303[___1926] < MinNeighborValue || ___1303[___1926] > ___2114[___1397]) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYFACE142) Invalid right neighbor value at face %lld:\n", (long long)___1926); ___3186("     right neighbor value = %lld, min value = %lld,", (long long)___1303[___1926], (long long)MinNeighborValue); ___3185(" max value = %lld.\n", (long long)___2114[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } else if (___1259[___1926] == TECIO_NO_NEIGHBORING_ELEM && ___1303[___1926] == TECIO_NO_NEIGHBORING_ELEM) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYFACE142) Both left and right neighbors are set to no neighboring element at face %lld.\n", (long long)___1926);
 #endif
NumErrs[___1397]++; return (-1); } } return 0; } } namespace { int32_t writeDecrementedIntegerArray( int32_t        ___1397, int32_t        ___2812, int32_t const* ___2099) { for (int32_t i = 0; i < ___2812; ++i) { if (!___4490(BlckFile[___1397], ___2099[i] - 1)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } } return 0; } } namespace { int32_t writeElements( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, int32_t           ___2806, int32_t const*    ___1259, int32_t const*    ___1303, int32_t           MinNeighborValue) { if (checkElements(___1397, ___2806, ___1259, ___1303, MinNeighborValue) != 0) return (-1); ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.leftElemsOffset, SEEK_SET); if (writeDecrementedIntegerArray(___1397, ___2806, ___1259) != 0) return -1; polyZoneWriteInfo.leftElemsOffset = (___1398)___4201(BlckFile[___1397]->File); ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.rightElemsOffset, SEEK_SET); if (writeDecrementedIntegerArray(___1397, ___2806, ___1303) != 0) return -1; polyZoneWriteInfo.rightElemsOffset = (___1398)___4201(BlckFile[___1397]->File); return 0; } } namespace { int32_t writeBoundaryConnectionCounts( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, int32_t           ___2778, int32_t const*    ___1253) { ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.connectionCountsOffset, SEEK_SET); if (polyZoneWriteInfo.boundaryConnectionSum == 0) { ___478(TecplotSDKBinaryFileVersion == 112); if (!(___4490(BlckFile[___1397], 0) && ___4490(BlckFile[___1397], 0))) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } } for (___2227 ___1926 = 0; ___1926 < ___2778; ___1926++) { polyZoneWriteInfo.boundaryConnectionSum += ___1253[___1926]; if (___1253[___1926] < 0) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYBCONN142) Invalid boundary connection count at position %lld:\n", (long long)___1926 + 1); ___3185("     boundary connection count = %lld.\n", (long long)___1253[___1926]);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], polyZoneWriteInfo.boundaryConnectionSum)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } } if (polyZoneWriteInfo.boundaryConnectionSum > ___4190[___1397]) {
 #if defined MAKEARCHIVE
___3184("Err: (TECPOLYBCONN142) Invalid number of boundary connections:\n"); ___3186("     number of boundary connections written = %lld, total number of boundary connections = %lld.", (long long)polyZoneWriteInfo.boundaryConnectionSum, (long long)___4190[___1397]);
 #endif
NumErrs[___1397]++; return (-1); } polyZoneWriteInfo.connectionCountsOffset = (___1398)___4201(BlckFile[___1397]->File); return 0; } } namespace { int32_t writeBoundaryConnectionElements( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, int32_t           ___2778, int32_t const*    ___1253, int32_t const*    ___1254, int32_t const*    ___1256) { ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.connectionElemsOffset, SEEK_SET); ___2227 LocalBndryConnCount = 0; for (___2227 ___1926 = 0; ___1926 < ___2778; ___1926++) { for (___2227 BIndex = 0; BIndex < ___1253[___1926]; BIndex++) { if (BIndex > 0 && ___1254[LocalBndryConnCount] == TECIO_NO_NEIGHBORING_ELEM) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYBCONN142) Partially obscured faces must specify no neighboring element first. See boundary connections for face %lld.\n", (long long)polyZoneWriteInfo.numBoundaryFacesWritten + ___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } if (___1254[LocalBndryConnCount] < TECIO_NO_NEIGHBORING_ELEM) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYBCONN142) Invalid boundary element value at boundary connections for face %lld:\n", (long long)polyZoneWriteInfo.numBoundaryFacesWritten + ___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } if (___1254[LocalBndryConnCount] == TECIO_NO_NEIGHBORING_ELEM && ___1256[LocalBndryConnCount] != TECIO_NO_NEIGHBORING_ZONE) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYBCONN142) Invalid boundary element/zone pair at boundary connections for face %lld:\n", (long long)___1926 + 1); ___3184("     Boundary elements specified as no neighboring element must also specify no neighboring zone.\n");
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1254[LocalBndryConnCount] - 1)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } LocalBndryConnCount++; } } polyZoneWriteInfo.connectionElemsOffset = (___1398)___4201(BlckFile[___1397]->File); return 0; } } namespace { int32_t writeBoundaryConnectionZones( int32_t           ___1397, PolyZoneWriteInfo& polyZoneWriteInfo, int32_t           ___2778, int32_t const*    ___1253, int32_t const*    ___1256) { ___4200(BlckFile[___1397]->File, polyZoneWriteInfo.connectionZonesOffset, SEEK_SET); int32_t LocalBndryConnCount = 0; for (___2227 ___1926 = 0; ___1926 < ___2778; ___1926++) { for (___2227 BIndex = 0; BIndex < ___1253[___1926]; BIndex++) { if (___1256[LocalBndryConnCount] < TECIO_NO_NEIGHBORING_ZONE) {
 #if defined MAKEARCHIVE
___3185("Err: (TECPOLYBCONN142) Invalid boundary zone value at boundary connections for face %lld:\n", (long long)___1926 + 1);
 #endif
NumErrs[___1397]++; return (-1); } else if (!___4490(BlckFile[___1397], ___1256[LocalBndryConnCount] - 1)) { WriteErr(___1397, "TECPOLYFACE142"); return -1; } LocalBndryConnCount++; } } polyZoneWriteInfo.connectionZonesOffset = (___1398)___4201(BlckFile[___1397]->File); return 0; } } int32_t ___3985( int32_t        ___1397, int32_t const* ___2806, int32_t const* ___1294, int32_t const* ___1297, int32_t const* ___1259, int32_t const* ___1303) { int32_t ___3359 = checkForPolyFacePreconditions(___1397); if (___3359 != 0) return ___3359; int32_t MinNeighborValue; ___3359 = getMinNeighborValue(___1397, MinNeighborValue); if (___3359 != 0) return ___3359; if (PolyZoneWriteInfos[___1397][___734[___1397]].numFacesWritten == 0 && PolyZoneWriteInfos[___1397][___734[___1397]].numBoundaryFacesWritten == 0) { calculateInitialPolyZoneOffsets(___1397); } PolyZoneWriteInfo& polyZoneWriteInfo = PolyZoneWriteInfos[___1397][___734[___1397]]; ___2227 BeginningFaceNodeSum = polyZoneWriteInfo.faceNodeSum; if (___4693[___1397] == ___1342) ___3359 = writePolyhedralFaceNodeOffsets(___1397, polyZoneWriteInfo, *___2806, ___1294); if (___3359 != 0) return ___3359; ___3359 = writeFaceNodes(___1397, polyZoneWriteInfo, BeginningFaceNodeSum, *___2806, ___1297); if (___3359 != 0) return ___3359; ___3359 = writeElements(___1397, polyZoneWriteInfo, *___2806, ___1259, ___1303, MinNeighborValue); if (___3359 != 0) return ___3359; polyZoneWriteInfo.numFacesWritten += *___2806; if (polyZoneWriteInfo.numFacesWritten == ___2159[___1397] && polyZoneWriteInfo.numBoundaryFacesWritten == ___4191[___1397]) { FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; } return 0; } int32_t ___3984( int32_t        ___1397, int32_t const* ___2778, int32_t const* ___1253, int32_t const* ___1254, int32_t const* ___1256) { if (CheckFile(___1397, "TECPOLYBCONN142") < 0) return (-1); if (PolyZoneWriteInfos[___1397][___734[___1397]].numFacesWritten == 0 && PolyZoneWriteInfos[___1397][___734[___1397]].numBoundaryFacesWritten == 0) { calculateInitialPolyZoneOffsets(___1397); } PolyZoneWriteInfo& polyZoneWriteInfo = PolyZoneWriteInfos[___1397][___734[___1397]]; if (___4191[___1397] > 0) { if (writeBoundaryConnectionCounts(___1397, polyZoneWriteInfo, *___2778, ___1253) != 0) return -1; if (writeBoundaryConnectionElements(___1397, polyZoneWriteInfo, *___2778, ___1253, ___1254, ___1256) != 0) return -1; if (writeBoundaryConnectionZones(___1397, polyZoneWriteInfo, *___2778, ___1253, ___1256) != 0) return -1; } polyZoneWriteInfo.numBoundaryFacesWritten += *___2778; if (polyZoneWriteInfo.numFacesWritten == ___2159[___1397] && polyZoneWriteInfo.numBoundaryFacesWritten == ___4191[___1397]) { FaceNeighborsOrMapWritten[___1397][___734[___1397]] = ___4226; } return 0; }
