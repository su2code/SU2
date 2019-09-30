 #pragma once
#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
inline double clampToDataTypeRange(double ___4298, FieldDataType_e ___1363) { switch (___1363) { case FieldDataType_Float: return ___650(___4298); case FieldDataType_Double: return ___489(___4298); case FieldDataType_Int32: return ___652(___4298); case FieldDataType_Int16: return ___651(___4298); case FieldDataType_Byte: return CONVERT_DOUBLE_TO_UINT8(___4298); case ___1365: return (___4298 < 1.0 ? 0.0 : 1.0); default: ___478(___1305); return ___2179; } }
