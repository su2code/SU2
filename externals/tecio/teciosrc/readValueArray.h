 #pragma once
#include "MASTER.h"
#include "GLOBAL.h"
#include "FileReaderInterface.h"
#include "LightweightVector.h"
#include "ioDescription.h"
 #define STRING_SIZE 100
 #define STRING_FORMAT "%99s"
namespace tecplot { namespace ___3933 { template<typename T, bool ___2025, int baseValue> ___372 readValues( ___1399& file, size_t               ___2796, T*                   ___4299, IODescription const& ___972); template<typename T, bool ___2025, int base> ___372 readValueArray( ___1399&  file, size_t                ___2865, size_t                ___2795, ___2240<T>& ___4299, IODescription const&  ___972); template<typename T> ___372 readMinMaxArray( ___1399& file, size_t               ___2865, size_t               ___2795, ___2481&         ___2480, IODescription const& ___972); template<typename T, bool ___2025> ___372 readValue( ___1399& file, T&                   ___4298, IODescription const& ___972); template<typename T, bool ___2025> ___372 readAndVerifyValue( ___1399& file, T const              expectedVal, IODescription const& ___972); ___372 readString( ___1399& file, size_t               length, ___473&           ___4298, IODescription const& ___972); ___372 readStringArray( ___1399& file, size_t               ___2865, size_t               ___2812, ___3816&         itemNameArray, IODescription const& ___972); }}
