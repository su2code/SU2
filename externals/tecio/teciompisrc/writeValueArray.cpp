#include "writeValueArray.h"
#include "ThirdPartyHeadersBegin.h"
#include <utility>
#include "ThirdPartyHeadersEnd.h"
#include "AsciiOutputInfo.h"
#include "SzlFileLoader.h"
namespace tecplot { namespace ___3933 {
 #if defined ASCII_ANNOTATE_TYPES
template <typename T> struct AsciiTypeString; template <> struct AsciiTypeString<char>          { static char const* typeString; }; template <> struct AsciiTypeString<uint8_t>       { static char const* typeString; }; template <> struct AsciiTypeString<int16_t>       { static char const* typeString; }; template <> struct AsciiTypeString<uint16_t>      { static char const* typeString; }; template <> struct AsciiTypeString<int32_t>       { static char const* typeString; }; template <> struct AsciiTypeString<uint32_t>      { static char const* typeString; }; template <> struct AsciiTypeString<uint64_t>      { static char const* typeString; }; template <> struct AsciiTypeString<float>         { static char const* typeString; }; template <> struct AsciiTypeString<double>        { static char const* typeString; }; template <> struct AsciiTypeString<StdPairUInt8>  { static char const* typeString; }; template <> struct AsciiTypeString<StdPairInt16>  { static char const* typeString; }; template <> struct AsciiTypeString<StdPairInt32>  { static char const* typeString; }; template <> struct AsciiTypeString<StdPairFloat>  { static char const* typeString; }; template <> struct AsciiTypeString<StdPairDouble> { static char const* typeString; }; char const* AsciiTypeString<char>::typeString          = "char"; char const* AsciiTypeString<uint8_t>::typeString       = "uint8_t"; char const* AsciiTypeString<int16_t>::typeString       = "int16_t"; char const* AsciiTypeString<uint16_t>::typeString      = "uint16_t"; char const* AsciiTypeString<int32_t>::typeString       = "int32_t"; char const* AsciiTypeString<uint32_t>::typeString      = "uint32_t"; char const* AsciiTypeString<uint64_t>::typeString      = "uint64_t"; char const* AsciiTypeString<float>::typeString         = "float"; char const* AsciiTypeString<double>::typeString        = "double"; char const* AsciiTypeString<StdPairUInt8>::typeString  = "2*uint8_t"; char const* AsciiTypeString<StdPairInt16>::typeString  = "2*int16_t"; char const* AsciiTypeString<StdPairInt32>::typeString  = "2*int32_t"; char const* AsciiTypeString<StdPairFloat>::typeString  = "2*float"; char const* AsciiTypeString<StdPairDouble>::typeString = "2*double";
 #endif
template ___372 encodeAsciiValue<uint32_t, false, 0>(char *str, int ___418, uint32_t const& ___4298);
 #define SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(T) \
 template<> \
 ___372 encodeAsciiValue<std::pair<T, T>, false, 0>(char* str, int ___418, std::pair<T, T> const& minMax) \
 { \
 REQUIRE(VALID_REF(str)); \
 REQUIRE(___418 > 0); \
 REQUIRE(minMax.first == minMax.first); \
 REQUIRE(minMax.second == minMax.second); \
 std::string ___1474 = std::string(___199<T, false>::___1474) + " " + ___199<T, false>::___1474; \
 int nChars = snprintf(str, ___418, ___1474.c_str(), \
 -___199<T, false>::size, minMax.first, -___199<T, false>::size, minMax.second); \
 ___372 ___2039 = (nChars > 0 && nChars < ___418); \
 ENSURE(VALID_BOOLEAN(___2039)); \
 return ___2039; \
 }
SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(uint8_t) SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(int16_t) SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(int32_t) SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(float) SPECIALIZE_MIN_MAX_ENCODING_FOR_TYPE(double) template <typename T, bool ___2025, int base> ___372 ___4563(FileWriterInterface& file, char const*          ___972, ___81           ___1251, size_t               ___2797, T const*             ___4299, size_t               ___4334  ) { ___372 ___2039 = ___4226; REQUIRE(file.___2041()); REQUIRE(VALID_DESCRIPTION(___972)); REQUIRE("extraID could have any value. NO_EXTRA_ID show only the description"); REQUIRE(___2797>0); REQUIRE(VALID_REF(___4299)); if ( file.___2002() ) { ASSERT_ONLY(___1393 beginningLocation = file.fileLoc()); std::string ___1418 = ___972; if ( ___1251 != ___2745 ) ___1418 += ___4187(___1251+1);
 #if defined ASCII_ANNOTATE_TYPES
___1418.append(1, ' ').append(AsciiTypeString<T>::typeString);
 #endif
if ( ___2797 == 1 ) file.fprintf("%*s  ", -___206, ___1418.c_str()); else file.fprintf("%*s\r\n", -___206, ___1418.c_str()); const int buffSize = 100; char buff[buffSize]; std::string outputBuffer; for ( size_t pos = 1; pos <= ___2797; ++pos ) { ___2039 = ___2039 && encodeAsciiValue<T, ___2025, base>(buff, buffSize, ___4299[pos-1]); outputBuffer.append(buff); if( (pos % ___4334) == 0 || (pos == ___2797 ) ) outputBuffer.append("\r\n"); else outputBuffer.append("  "); } file.fwrite(outputBuffer.c_str(), sizeof(char), outputBuffer.size()); ASSERT_ONLY(___1393 endingLocation = file.fileLoc()); ASSERT_ONLY(___1393 outputSize = (___1393)___206 + 2 + outputBuffer.size()); ___478(endingLocation - beginningLocation == outputSize); } else { file.fwrite(___4299, sizeof(T), ___2797); } ENSURE(VALID_BOOLEAN(___2039)); return ___2039; } template <typename OutType> ___372 ___4528( FileWriterInterface&        file, char const*                 ___972, ___81                  ___1251, size_t                      ___2797, ___2479 const*               ___4299, size_t                      ___4334  ) { ___2240<std::pair<OutType, OutType> > outputArray; ___372 ___2039 = outputArray.alloc(___2797); if (___2039) { for (size_t i = 0; i < ___2797; ++i) { outputArray[i].first = static_cast<OutType>(___4299[i].minValue()); outputArray[i].second = static_cast<OutType>(___4299[i].maxValue()); } ___2039 = ___4563<std::pair<OutType, OutType>, false, 0>(file, ___972, ___1251, ___2797, &outputArray[0], ___4334); } return ___2039; } template <typename T, bool ___2025> uint64_t arrayValueSizeInFile(bool ___2002) { if (___2002) return ___199<T, ___2025>::size + ASCII_SPACING_LEN; return sizeof(T); } template <typename T, bool ___2025> uint64_t arraySizeInFile(size_t ___2797, bool ___2002) { uint64_t charsPerNumber = arrayValueSizeInFile<T, ___2025>(___2002); ___478(charsPerNumber > 0); uint64_t ___3358 = static_cast<uint64_t>(___2797) * charsPerNumber; if (___2002) ___3358 += static_cast<uint64_t>(___206) + ASCII_SPACING_LEN; return ___3358; } template <typename T, bool ___2025> uint64_t valueSizeInFile(bool ___2002) { return arraySizeInFile<T, ___2025>(1, ___2002); }
 #define INSTANTIATE_FOR_TYPE(T, ___2025) \
 template ___372 ___4563< T, ___2025, 0 >( \
 FileWriterInterface& file, \
 char const*          ___972, \
 ___81           ___1251, \
 size_t               ___2797, \
 T const*             ___4299, \
 size_t               ___4334); \
 template uint64_t arrayValueSizeInFile< T, ___2025 > (bool ___2002); \
 template uint64_t valueSizeInFile< T, ___2025 > (bool ___2002); \
 template uint64_t arraySizeInFile< T, ___2025 > (size_t ___2797, bool ___2002);
INSTANTIATE_FOR_TYPE(char, false) INSTANTIATE_FOR_TYPE(uint8_t, false) INSTANTIATE_FOR_TYPE(uint8_t, true) INSTANTIATE_FOR_TYPE(int16_t, false) INSTANTIATE_FOR_TYPE(uint16_t, false) INSTANTIATE_FOR_TYPE(uint16_t, true) INSTANTIATE_FOR_TYPE(int32_t, false) INSTANTIATE_FOR_TYPE(uint32_t, false) INSTANTIATE_FOR_TYPE(uint64_t, false) INSTANTIATE_FOR_TYPE(uint64_t, true) INSTANTIATE_FOR_TYPE(float, false) INSTANTIATE_FOR_TYPE(double, false) INSTANTIATE_FOR_TYPE(StdPairFloat, false) INSTANTIATE_FOR_TYPE(StdPairDouble, false) INSTANTIATE_FOR_TYPE(StdPairInt32, false) INSTANTIATE_FOR_TYPE(StdPairInt16, false) INSTANTIATE_FOR_TYPE(StdPairUInt8, false) template ___372 ___4563< uint32_t, false, 1 >(FileWriterInterface& file, char const* ___972, ___81 ___1251, size_t ___2797, uint32_t const* ___4299, size_t ___4334);
 #define INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(T) \
 template ___372 ___4528 <T>( \
 FileWriterInterface&        file, \
 char const*                 ___972, \
 ___81                  ___1251, \
 size_t                      ___2797, \
 ___2479 const*               ___4299, \
 size_t                      ___4334);
INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(uint8_t) INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(int16_t) INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(int32_t) INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(float) INSTANTIATE_MIN_MAX_OUTPUT_FOR_TYPE(double) }}
