#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "stdafx.h"
#include "MASTER.h"

#define TECPLOTENGINEMODULE

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#define DATASET0MODULE
#include "GLOBAL.h"
#include "TASSERT.h"
#include "Q_UNICODE.h"
#include "ALLOC.h"
#include "ARRLIST.h"
#include "DATASET.h"
#include "SET.h"
#include "FILESTREAM.h"
#include "Q_MSG.h"
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#include "DATASET0.h"

using namespace tecplot::strutil;
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/*
 * Low level dataset functions.  No references to zones, vars or
 * the DataSet_s master structure here.
 */


/*
 */
void OutOfMemoryMsg(void)
{
    ErrMsg(translate("Cannot allocate enough memory for this operation."));
} /* OutOfMemoryMsg() */


/**
 */
FieldData_pa FieldDataAlloc(Boolean_t doTrackVarSharing)
{
    FieldData_pa Result;

    Result = (FieldData_pa)ALLOC_ITEM(struct _FieldData_a, "FieldDataPtr");
    if (Result != NULL)
    {
        Result->Data = NULL;

        #if defined TECPLOTKERNEL /* TecIO doesn't require these features yet */
/* CORE SOURCE CODE REMOVED */
        #else /* ...for TecIO only */
        Result->GetValueCallback[0] = NULL;
        Result->SetValueCallback[0] = NULL;
        #endif

        #if defined TECPLOTKERNEL /* TecIO doesn't require these features yet */
/* CORE SOURCE CODE REMOVED */
        #endif

        Result->Type             = FieldDataType_Invalid;
        Result->ValueLocation    = ValueLocation_Invalid;
        #if defined TECPLOTKERNEL /* TecIO doesn't require these features yet */
/* CORE SOURCE CODE REMOVED */
        #endif
        Result->NumValues        = 0;
        #if defined TECPLOTKERNEL /* TecIO doesn't require these features yet */
/* CORE SOURCE CODE REMOVED */
        #endif
    }

    ENSURE(VALID_REF(Result) || Result == NULL);
    return Result;
}

/**
 * Most clients should not call this function but FieldDataCleanup() instead.
 * An exception to this would be Tecplot's own storable load-on-demand
 * functions.
 */
void FieldDataDeallocData(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
        if (FieldData->Data != NULL)
        {
            /* Hack to remove 'deleting void* is undefined' warning... */
            char *Tmp = (char *)FieldData->Data;
            FREE_ARRAY(Tmp, "FieldData _Data");
            FieldData->Data = NULL;
        }

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

    ENSURE(FieldData->Data == NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 */
void FieldDataCleanup(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
    FieldDataDeallocData(FieldData);
# endif
}

/**
 */
void FieldDataDealloc(FieldData_pa *FieldData,
                      Boolean_t     DoTrackVarSharing)
{
    REQUIRE(VALID_REF(FieldData));
    REQUIRE(VALID_REF(*FieldData) || *FieldData == NULL);
    #if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
    #endif

    if (*FieldData != NULL)
    {
        #if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
        #endif
        {
            FieldDataCleanup(*FieldData);

            #if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
            #endif

            FREE_ITEM(*FieldData, "field data");
        }
        *FieldData = NULL;
    }

    ENSURE(*FieldData == NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

template <typename T>
static void copyTypedValueArray(void*     DstArray,
                                LgIndex_t DstStart,
                                void*     SrcArray,
                                LgIndex_t SrcStart,
                                LgIndex_t SrcEnd)
{
    REQUIRE(VALID_REF(DstArray));
    REQUIRE(DstStart >= 0);
    REQUIRE(VALID_REF(SrcArray));
    REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd+1);
    REQUIRE(DstArray != SrcArray);

    size_t const numBytes = sizeof(T) * (SrcEnd - SrcStart + 1);
    if (numBytes != 0)
    {
        T const* SrcPtr = ((T const*)SrcArray) + SrcStart;
        T* DstPtr       = ((T*)DstArray) + DstStart;
        memcpy(DstPtr, SrcPtr, numBytes);
    }
}

/**
 * DstArray and SrcArray are aligned on proper word boundaries.
 */
void CopyTypedValueArray(FieldDataType_e  ValueType,
                         void            *DstArray,
                         LgIndex_t        DstStart,
                         void            *SrcArray,
                         LgIndex_t        SrcStart,
                         LgIndex_t        SrcEnd)
{
    REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
            ValueType != FieldDataType_Bit);
    REQUIRE(VALID_REF(DstArray));
    REQUIRE(DstStart >= 0);
    REQUIRE(VALID_REF(SrcArray));
    REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd+1);
    REQUIRE(DstArray != SrcArray);

    switch (ValueType)
    {
        case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove check */
        case FieldDataType_Double :
        {
            CHECK(sizeof(UInt64_t) == 8 && sizeof(double) == 8);
            copyTypedValueArray<UInt64_t>(DstArray,
                                          DstStart,
                                          SrcArray,
                                          SrcStart,
                                          SrcEnd);
        } break;
        case FieldDataType_Float :
        case FieldDataType_Int32 :
        {
            CHECK(sizeof(UInt32_t) == 4 && sizeof(float) == 4);
            copyTypedValueArray<UInt32_t>(DstArray,
                                          DstStart,
                                          SrcArray,
                                          SrcStart,
                                          SrcEnd);
        } break;
        case FieldDataType_Int16 :
        {
            CHECK(sizeof(UInt16_t) == 2);
            copyTypedValueArray<UInt16_t>(DstArray,
                                          DstStart,
                                          SrcArray,
                                          SrcStart,
                                          SrcEnd);
        } break;
        case FieldDataType_Byte :
        {
            copyTypedValueArray<Byte_t>(DstArray,
                                        DstStart,
                                        SrcArray,
                                        SrcStart,
                                        SrcEnd);
        } break;
        default : CHECK(FALSE);
    }
}

/**
 * SrcArray is aligned on proper word boundaries.
 */
void SwapBytesInTypedValueArray(FieldDataType_e  ValueType,
                                void            *SrcArray,
                                LgIndex_t        SrcStart,
                                LgIndex_t        SrcEnd,
                                LgIndex_t        SrcSkip)
{
    REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
            ValueType != FieldDataType_Bit);
    REQUIRE(VALID_REF(SrcArray));
    REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
    REQUIRE(SrcSkip > 0);

    switch (ValueType)
    {
        case FieldDataType_Int64: CHECK(FALSE); /* Future work: remove CHECK */
        case FieldDataType_Double:
        {
            /* swap 8 bytes blocks */
            UInt64_t *SrcPtr = ((UInt64_t *)SrcArray) + SrcStart;
            UInt64_t *SrcPtrEnd = ((UInt64_t *)SrcArray) + SrcEnd;
            CHECK(sizeof(UInt64_t) == 8 && sizeof(double) == 8);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_8_BYTES(SrcPtr);
                SrcPtr += SrcSkip;
            }
        } break;
        case FieldDataType_Float:
        case FieldDataType_Int32:
        {
            /* swap 4 bytes blocks */
            UInt32_t *SrcPtr = ((UInt32_t *)SrcArray) + SrcStart;
            UInt32_t *SrcPtrEnd = ((UInt32_t *)SrcArray) + SrcEnd;
            CHECK(sizeof(UInt32_t) == 4 && sizeof(float) == 4);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_4_BYTES(SrcPtr);
                SrcPtr += SrcSkip;
            }
        } break;
        case FieldDataType_Int16:
        {
            /* swap 4 bytes blocks */
            UInt16_t *SrcPtr = ((UInt16_t *)SrcArray) + SrcStart;
            UInt16_t *SrcPtrEnd = ((UInt16_t *)SrcArray) + SrcEnd;
            CHECK(sizeof(UInt16_t) == 2);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_2_BYTES(SrcPtr);
                SrcPtr += SrcSkip;
            }
        } break;
        case FieldDataType_Byte:
        case FieldDataType_Bit:
        {
            /* nothing to do */
        } break;
        default: CHECK(FALSE);
    }
}

/**
 * Same as SwapBytesInTypedValueArray, but does extra work.  Doesn't assume
 * DstArray and SrcArray are aligned on proper word boundaries.
 */
void SwapBytesInUnalignedTypedValueArray(FieldDataType_e  ValueType,
                                         void            *SrcArray,
                                         LgIndex_t        SrcStart,
                                         LgIndex_t        SrcEnd,
                                         LgIndex_t        SrcSkip)
{
    REQUIRE(VALID_FIELD_DATA_TYPE(ValueType) &&
            ValueType != FieldDataType_Bit);
    REQUIRE(VALID_REF(SrcArray));
    REQUIRE(0 <= SrcStart && SrcStart <= SrcEnd);
    REQUIRE(SrcSkip > 0);

    switch (ValueType)
    {
        case FieldDataType_Int64: CHECK(FALSE); /* Future work: remove CHECK */
        case FieldDataType_Double:
        {
            /* swap 8-byte blocks */
            Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart * sizeof(UInt64_t);
            Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd * sizeof(UInt64_t);
            size_t byte_skip = SrcSkip * sizeof(UInt64_t);
            CHECK(sizeof(UInt64_t) == 8 && sizeof(double) == 8);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_8_BYTES_1_AT_A_TIME(SrcPtr);
                SrcPtr += byte_skip;
            }
        } break;
        case FieldDataType_Float:
        case FieldDataType_Int32:
        {
            /* swap 4-byte blocks */
            Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart * sizeof(UInt32_t);
            Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd * sizeof(UInt32_t);
            size_t byte_skip = SrcSkip * sizeof(UInt32_t);
            CHECK(sizeof(UInt32_t) == 4 && sizeof(float) == 4);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_4_BYTES_1_AT_A_TIME(SrcPtr);
                SrcPtr += byte_skip;
            }
        } break;
        case FieldDataType_Int16:
        {
            /* swap 2-byte blocks */
            Byte_t *SrcPtr = ((Byte_t *)SrcArray) + SrcStart * sizeof(UInt16_t);
            Byte_t *SrcPtrEnd = ((Byte_t *)SrcArray) + SrcEnd * sizeof(UInt16_t);
            size_t byte_skip = SrcSkip * sizeof(UInt16_t);
            CHECK(sizeof(UInt16_t) == 2);
            while (SrcPtr <= SrcPtrEnd)
            {
                REVERSE_2_BYTES_1_AT_A_TIME(SrcPtr);
                SrcPtr += byte_skip;
            }
        } break;
        case FieldDataType_Byte:
        case FieldDataType_Bit:
        {
            /* No swapping required. */
        } break;
        default: CHECK(FALSE);
    }
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined DEBUG_FIELDVALUES
# define DEBUG_FIELDVALUES_BAD_VALUE 0x11
static unsigned char BadValueStr[] =
{
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE,
    DEBUG_FIELDVALUES_BAD_VALUE
};
/**
 * If Tecplot is responsible for managing (allocating and deallocating) the
 * raw data then
 */
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else
#   define FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, type) \
                   ((sizeof(type) < 4) /* cannot make reliably test with less than four bytes */ || \
                    memcmp(BadValueStr,((char *)((fd)->Data))+sizeof(type)*(pt), sizeof(type)) != 0)
# endif
#else
# define FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, type) TRUE
#endif



/**
 * Used in macros, thus not static
 */
double STDCALL GetFieldValueForFloat(const FieldData_pa fd,
                                     LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));

    double Result = (double)GetFieldDataFloatPtr(fd)[pt];

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 * Used in macros, thus not static
 */
double STDCALL GetFieldValueForDouble(const FieldData_pa fd,
                                      LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));

    double Result = GetFieldDataDoublePtr(fd)[pt];

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
double STDCALL GetFieldValueForInt32(const FieldData_pa fd,
                                     LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));

    double Result = (double)GetFieldDataInt32Ptr(fd)[pt];

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
double STDCALL GetFieldValueForInt16(const FieldData_pa fd,
                                     LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));

    double Result = (double)GetFieldDataInt16Ptr(fd)[pt];

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
double STDCALL GetFieldValueForByte(const FieldData_pa fd,
                                    LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(fd->Type == FieldDataType_Byte);
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Byte_t));

    double Result = (double)GetFieldDataBytePtr(fd)[pt];

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
double STDCALL GetFieldValueForBit(const FieldData_pa fd,
                                   LgIndex_t          pt)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(fd->Type == FieldDataType_Bit);
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt / 8, Byte_t));

    LgIndex_t ByteOffset = pt / 8;
    Byte_t    BitMask    = (01 << (pt % 8));

    Byte_t *byte_array = GetFieldDataBytePtr(fd);
    double Result = (byte_array[ByteOffset] & BitMask) ? 1.0 : 0.0;

    return Result;
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
FieldValueGetFunction_pf DetermineFieldDataGetFunction(FieldDataType_e DataType,
                                                       Boolean_t       IsFragmented)
{
    FieldValueGetFunction_pf Result;

    REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
    REQUIRE(VALID_BOOLEAN(IsFragmented));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
        switch (DataType)
        {
            case FieldDataType_Float :
            {
                Result = GetFieldValueForFloat;
            } break;
            case FieldDataType_Double :
            {
                Result = GetFieldValueForDouble;
            } break;
            case FieldDataType_Int32 :
            {
                Result = GetFieldValueForInt32;
            } break;
            case FieldDataType_Int16 :
            {
                Result = GetFieldValueForInt16;
            } break;
            case FieldDataType_Byte :
            {
                Result = GetFieldValueForByte;
            } break;
            case FieldDataType_Bit :
            {
                Result = GetFieldValueForBit;
            } break;
            default :
            {
                CHECK(FALSE);
                Result = NULL; /* satisfy compiler */
            } break;
        }
    }
    return (Result);
}

/**
 */
static void STDCALL SetFieldValueForFloat(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    GetFieldDataFloatPtr(fd)[pt] = CONVERT_DOUBLE_TO_FLOAT(val);

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, float));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static void STDCALL SetFieldValueForDouble(FieldData_pa fd,
                                           LgIndex_t    pt,
                                           double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    GetFieldDataDoublePtr(fd)[pt] = CLAMP_DOUBLE(val);

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, double));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static void STDCALL SetFieldValueForInt32(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    GetFieldDataInt32Ptr(fd)[pt] = CONVERT_DOUBLE_TO_INT32(val);

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int32_t));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static void STDCALL SetFieldValueForInt16(FieldData_pa fd,
                                          LgIndex_t    pt,
                                          double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    GetFieldDataInt16Ptr(fd)[pt] = CONVERT_DOUBLE_TO_INT16(val);

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Int16_t));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static void STDCALL SetFieldValueForByte(FieldData_pa fd,
                                         LgIndex_t    pt,
                                         double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(fd->Type == FieldDataType_Byte);
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    if (val < 1.0)
        GetFieldDataBytePtr(fd)[pt] = 0;
    else if (val > 255.0)
        GetFieldDataBytePtr(fd)[pt] = 255;
    else
        GetFieldDataBytePtr(fd)[pt] = (Byte_t)val;

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt, Byte_t));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
static void STDCALL SetFieldValueForBit(FieldData_pa fd,
                                        LgIndex_t    pt,
                                        double       val)
{
    REQUIRE(VALID_REF(fd));
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    REQUIRE(fd->Type == FieldDataType_Bit);
    REQUIRE(0 <= pt && pt < GetFieldDataNumValues(fd));
    REQUIRE("val can have any value");

    LgIndex_t ByteOffset = pt / 8;
    Byte_t    BitMask    = (01 << (pt % 8));

    if (val < 1.0)
        GetFieldDataBytePtr(fd)[ByteOffset] &= ~BitMask;
    else
        GetFieldDataBytePtr(fd)[ByteOffset] |= BitMask;

    ENSURE(FIELD_DATA_VALUE_IS_INITIALIZED(fd, pt / 8, Byte_t));
}


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/**
 */
FieldValueSetFunction_pf DetermineFieldDataSetFunction(FieldDataType_e DataType,
                                                       Boolean_t       IsFragmented)
{
    FieldValueSetFunction_pf Result;

    REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
    REQUIRE(VALID_BOOLEAN(IsFragmented));

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
    {
        switch (DataType)
        {
            case FieldDataType_Float :
            {
                Result = SetFieldValueForFloat;
            } break;
            case FieldDataType_Double :
            {
                Result = SetFieldValueForDouble;
            } break;
            case FieldDataType_Int32 :
            {
                Result = SetFieldValueForInt32;
            } break;
            case FieldDataType_Int16 :
            {
                Result = SetFieldValueForInt16;
            } break;
            case FieldDataType_Byte :
            {
                Result = SetFieldValueForByte;
            } break;
            case FieldDataType_Bit :
            {
                Result = SetFieldValueForBit;
            } break;
            default :
            {
                CHECK(FALSE);
                Result = NULL; /* satisfy compiler */
            } break;
        }
    }
    return (Result);
}


/**
 */
Int64_t FieldDataGetBytesNeeded(LgIndex_t       NumValues,
                                FieldDataType_e DataType)
{
    Int64_t Result = 0; /* ...quite compiler */

    REQUIRE(NumValues >= 0);
    REQUIRE(VALID_FIELD_DATA_TYPE(DataType));

    switch (DataType)
    {
        case FieldDataType_Float:  Result = ((Int64_t)NumValues)*sizeof(float);        break;
        case FieldDataType_Double: Result = ((Int64_t)NumValues)*sizeof(double);       break;
        case FieldDataType_Int32:  Result = ((Int64_t)NumValues)*sizeof(LgIndex_t);    break;
        case FieldDataType_Int16:  Result = ((Int64_t)NumValues)*sizeof(SmInteger_t);  break;
        case FieldDataType_Byte:   Result = ((Int64_t)NumValues)*sizeof(Byte_t);       break;
        case FieldDataType_Bit:    Result = ((Int64_t)(NumValues+7)/8)*sizeof(Byte_t); break;
        default: CHECK(FALSE); break;
    }

    ENSURE(Result >= 0);
    return Result;
}

/**
 * On the SGI, HP, Sun and Itanium Linux machines 64 bit objects such as
 * doubles must be 8 byte aligned while on all other machines 32 bit alignment
 * suffices. Some allow 1 byte alignment but we won't bother with that.
 */
#if defined IRISX || defined HPUX || defined SUNX
# define SIZEOF_LARGEST_OBJECT_TO_ALIGN sizeof(Int64_t)
#else
# define SIZEOF_LARGEST_OBJECT_TO_ALIGN sizeof(Int32_t)
#endif

/**
 */
Boolean_t IsOffsetAlignedForFieldDataType(FieldDataType_e FieldDataType,
                                          Int64_t         Offset)
{
    REQUIRE(VALID_FIELD_DATA_TYPE(FieldDataType));
    REQUIRE(Offset >= 0);

    Int64_t SizeOfType = FieldDataGetBytesNeeded(1, FieldDataType);
    if (SizeOfType > (Int64_t)SIZEOF_LARGEST_OBJECT_TO_ALIGN)
        SizeOfType = SIZEOF_LARGEST_OBJECT_TO_ALIGN;

    Boolean_t HasValidAlignment = (Offset % SizeOfType == 0);

    ENSURE(VALID_BOOLEAN(HasValidAlignment));
    return HasValidAlignment;
}

/**
 */
Int64_t GetAlignedOffsetForFieldDataType(FieldDataType_e FieldDataType,
                                         Int64_t         Offset)
{
    REQUIRE(VALID_FIELD_DATA_TYPE(FieldDataType));
    REQUIRE(Offset >= 0);

    Int64_t SizeOfType = FieldDataGetBytesNeeded(1, FieldDataType);
    if (SizeOfType > (Int64_t)SIZEOF_LARGEST_OBJECT_TO_ALIGN)
        SizeOfType = SIZEOF_LARGEST_OBJECT_TO_ALIGN;

    Int64_t NumBytesPastAlignment = (Offset % SizeOfType);
    Int64_t Result = Offset - NumBytesPastAlignment;

    ENSURE(0 <= Result && Result <= Offset);
    ENSURE(IsOffsetAlignedForFieldDataType(FieldDataType, Result));
    return Result;
}

/**
 */
void FieldDataDefineData(FieldData_pa    FieldData,
                         LgIndex_t       NumValues,
                         FieldDataType_e DataType,
                         ValueLocation_e ValueLocation)
{
    REQUIRE(VALID_REF(FieldData));
    REQUIRE(NumValues >= 0);
    REQUIRE(VALID_FIELD_DATA_TYPE(DataType));
    REQUIRE(VALID_ENUM(ValueLocation, ValueLocation_e));

    /*
     * Remove any old data (transformed UVW is one example that calls this
     * function with a non-null data pointer when switching the value location
     * when style changes the value location and therefore the amount of data
     * allocated.)
     */
    FieldDataCleanup(FieldData);

    /*
     * The reference count is not modified here. This function only allocates the
     * structure and makes adjustments to the some of the members. The reference
     * count was initialized when the structure was initially created and the
     * structure may be shared before the data portion is even allocated.
     */
    FieldData->NumValues                = NumValues;
    FieldData->Type                     = DataType;
    FieldData->ValueLocation            = ValueLocation;
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# else /* ...for TecIO only */
    FieldData->GetValueCallback[0] = (void *)DetermineFieldDataGetFunction(DataType, FALSE);
    FieldData->SetValueCallback[0] = (void *)DetermineFieldDataSetFunction(DataType, FALSE);
#endif

    ENSURE(FieldData->Data == NULL);
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 */
Boolean_t FieldDataAllocData(FieldData_pa FieldData,
                             Boolean_t    ShowErrMsg)
{
    REQUIRE(VALID_REF(FieldData));
    REQUIRE(FieldData->Type != FieldDataType_Invalid); /* ...must call FieldDataDefineData first */
    REQUIRE(FieldData->Data == NULL);
    REQUIRE(VALID_BOOLEAN(ShowErrMsg));

    /*
     * The size of size_t may be smaller than our unsigned 64 bit integer value
     * so we might have to squeeze it down possibly loosing precision.
     */
    Int64_t ActualBytesNeeded = FieldDataGetBytesNeeded(FieldData->NumValues, FieldData->Type);
    size_t  BytesToAllocate   = (size_t)ActualBytesNeeded;

    /*
     * 64 bit architectures are effectively unlimited in their allocation size
     * while 32 architectures are limited to 4GB (some may limit further to 2GB
     * which will be borne out by the call to malloc).
     */
    CHECK(sizeof(size_t) == 4 || sizeof(size_t) == 8);
    Boolean_t IsOk = (FieldData->NumValues <= MAXINDEX &&
                      IMPLICATION(sizeof(size_t) == 4,
                                  ActualBytesNeeded <= (Int64_t)0xffffffff));
    if (IsOk)
    {
        if (FieldData->NumValues > 0)
        {
            FieldData->Data = (void *)ALLOC_ARRAY(BytesToAllocate, char, "FieldData's Data");
            #if defined DEBUG_FIELDVALUES
            {
                if (FieldData->Data != NULL)
                    memset(FieldData->Data, DEBUG_FIELDVALUES_BAD_VALUE, BytesToAllocate);
            }
            #endif
            /*
             * For bit type data zero the last byte in the data array. We do
             * this because NumValues is probably not a multiple of 8 bits and
             * thus the valid bit values will not occupy all bits of the last
             * byte. By zeroing the unused bits at the end of the array we
             * produce consistent data files when written to disk.
             */
            if (FieldData->Type == FieldDataType_Bit && FieldData->Data != NULL)
                ((char*)FieldData->Data)[BytesToAllocate-1] = '\0';
        }
        IsOk = (FieldData->NumValues == 0 ||
                FieldData->Data != NULL);
        if (!IsOk && ShowErrMsg)
            OutOfMemoryMsg();
    }
    else if (ShowErrMsg)
        ErrMsg(translate("Storage limit (%ld) exceeded for a single variable."), (long)MAXINDEX);

# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif

    ENSURE(VALID_REF(FieldData->Data) || FieldData->Data == NULL);
    ENSURE(VALID_BOOLEAN(IsOk));
    return IsOk;
}

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if !defined NO_ASSERTS
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

/**
 * Allocates a field data pointer with space for "num_pts" of field data type
 * "field_data_type" nodal values.
 *
 * IMPORTANT:
 *   This field data may NOT be used for zones but only for things like
 *   geometries or other temporary field data that will never be placed
 *   into a COB or zone.
 */
FieldData_pa AllocScratchNodalFieldDataPtr(LgIndex_t       NumValues,
                                           FieldDataType_e Type,
                                           Boolean_t       ShowErrMsg)
{
    REQUIRE(NumValues >= 0);
    REQUIRE(VALID_FIELD_DATA_TYPE(Type));
    REQUIRE(VALID_BOOLEAN(ShowErrMsg));

    FieldData_pa Result = FieldDataAlloc(FALSE);
    if (Result != NULL)
    {
        FieldDataDefineData(Result, NumValues, Type, ValueLocation_Nodal);
        if (!FieldDataAllocData(Result, ShowErrMsg))
            FieldDataDealloc(&Result, FALSE);
    }
    else if (ShowErrMsg)
        OutOfMemoryMsg();

    ENSURE(VALID_REF(Result) || Result == NULL);
# if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
# endif
    ENSURE(IMPLICATION(Result != NULL,
                       (Result->NumValues >= 0                &&
                        IMPLICATION(Result->NumValues != 0,
                                    VALID_REF(Result->Data))  &&
                        VALID_FIELD_DATA_TYPE(Result->Type))));

    return Result;
}


/**
 * Frees memory allocated with AllocScratchNodalFieldDataPtr().
 *
 * @param ScratchFieldData
 *   Scratch field data pointer to deallocate. This should NEVER be a field
 *   data from a zone or COB. See note in AllocScratchNodalFieldDataPtr().
 */
void DeallocScratchNodalFieldDataPtr(FieldData_pa *FieldDataRef)
{
    FieldDataDealloc(FieldDataRef,
                     FALSE); /* DoTrackVarSharing */
}


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldDataType_e GetFieldDataType_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    FieldDataType_e Result = GetFieldDataType_MACRO(FieldData);

    ENSURE(VALID_FIELD_DATA_TYPE(Result));
    return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldValueGetFunction_pf GetFieldDataGetFunction_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    FieldValueGetFunction_pf Result = GetFieldDataGetFunction_MACRO(FieldData);

    ENSURE(VALID_FN_REF(Result));
    return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
FieldValueSetFunction_pf GetFieldDataSetFunction_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    FieldValueSetFunction_pf Result = GetFieldDataSetFunction_MACRO(FieldData);

    ENSURE(VALID_FN_REF(Result));
    return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
LgIndex_t GetFieldDataNumValues_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    LgIndex_t Result = GetFieldDataNumValues_MACRO(FieldData);

    ENSURE(Result >= 0);
    return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
ValueLocation_e GetFieldDataValueLocation_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    ValueLocation_e Result = GetFieldDataValueLocation_MACRO(FieldData);

    ENSURE(Result == ValueLocation_Invalid || /* i.e. pending assignment */
           VALID_ENUM(Result, ValueLocation_e));
    return Result;
}
#endif /* !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS */


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
Boolean_t IsFieldDataDirectAccessAllowed_FUNC(FieldData_pa FieldData)
{
    REQUIRE(VALID_REF(FieldData));

    Boolean_t Result = IsFieldDataDirectAccessAllowed_MACRO(FieldData);

    ENSURE(VALID_BOOLEAN(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
float *GetFieldDataFloatPtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    float *Result = GetFieldDataFloatPtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
double *GetFieldDataDoublePtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    double *Result = GetFieldDataDoublePtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
Int64_t *GetFieldDataInt64Ptr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    Int64_t *Result = GetFieldDataInt64Ptr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
Int32_t *GetFieldDataInt32Ptr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    Int32_t *Result = GetFieldDataInt32Ptr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 */
Int16_t *GetFieldDataInt16Ptr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    Int16_t *Result = GetFieldDataInt16Ptr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * No byte ordering or alignment issues with single bytes (which are also used with the "Bit" type)
 */
Byte_t *GetFieldDataBytePtr_FUNC(FieldData_pa fd)
{
    /*
     * This function gets called for Byte and Bit types, but we cannot REQUIRE
     * those types because it is also used for non-aligned values.  We can't
     * check for non-aligned because we might be copying aligned bytes to a
     * non-aligned location.
     */
    REQUIRE(VALID_REF(fd));

    Byte_t *Result = GetFieldDataBytePtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/**
 * Gets a ptr to 2-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt16_t *GetFieldData2BytePtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));
    REQUIRE(fd->Type == FieldDataType_Int16);

    UInt16_t *Result = GetFieldData2BytePtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * Gets a ptr to 4-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt32_t *GetFieldData4BytePtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));
    REQUIRE(fd->Type == FieldDataType_Int32 || fd->Type == FieldDataType_Float);

    UInt32_t *Result = GetFieldData4BytePtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * Gets a ptr to 8-byte blocks regardless of byte ordering, but still has to
 * worry about byte alignment
 */
UInt64_t *GetFieldData8BytePtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));
    REQUIRE(fd->Type == FieldDataType_Int64 || fd->Type == FieldDataType_Double);

    UInt64_t *Result = GetFieldData8BytePtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif


#if !defined USE_MACROS_FOR_FIELD_DATA_FUNCTIONS
/*
 * WARNING: GetFieldDataVoidPtr checks nothing, and thus should only be
 * used with extreme caution (that is, checking the alignment
 * and byte order by hand).
 */
void *GetFieldDataVoidPtr_FUNC(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));

    void *Result = GetFieldDataVoidPtr_MACRO(fd);

    ENSURE(VALID_REF(Result));
    return Result;
}
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


/*
 */
void CopyFieldValue(FieldData_pa  dst,
                    LgIndex_t     dstindex,
                    FieldData_pa  src,
                    LgIndex_t     srcindex)
{
    REQUIRE(VALID_REF(dst));
    REQUIRE(VALID_REF(src));
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
    REQUIRE(dstindex >= 0 && dstindex < GetFieldDataNumValues(dst) &&
            srcindex >= 0 && srcindex < GetFieldDataNumValues(src));

    Boolean_t DoBruteForceCopy = TRUE;

    if (IsFieldDataDirectAccessAllowed(src) &&
        IsFieldDataDirectAccessAllowed(dst) &&
        GetFieldDataType(src) == GetFieldDataType(dst))
    {
        switch (GetFieldDataType(src))
        {
            case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove and let fall through */
            case FieldDataType_Double :
            {
                CHECK(sizeof(UInt64_t) == 8 && sizeof(double) == 8);
                UInt64_t *dst_ptr = GetFieldData8BytePtr(dst) + dstindex;
                UInt64_t *src_ptr = GetFieldData8BytePtr(src) + srcindex;
                *dst_ptr = *src_ptr;
                DoBruteForceCopy = FALSE;
            } break;
            case FieldDataType_Float :
            case FieldDataType_Int32 :
            {
                CHECK(sizeof(UInt32_t) == 4 && sizeof(float) == 4);
                UInt32_t *dst_ptr = GetFieldData4BytePtr(dst) + dstindex;
                UInt32_t *src_ptr = GetFieldData4BytePtr(src) + srcindex;
                *dst_ptr = *src_ptr;
                DoBruteForceCopy = FALSE;
            } break;
            case FieldDataType_Int16 :
            {
                CHECK(sizeof(UInt16_t) == 2);
                UInt16_t *dst_ptr = GetFieldData2BytePtr(dst) + dstindex;
                UInt16_t *src_ptr = GetFieldData2BytePtr(src) + srcindex;
                *dst_ptr = *src_ptr;
            } break;
            case FieldDataType_Byte :
            {
                GetFieldDataBytePtr(dst)[dstindex] = GetFieldDataBytePtr(src)[srcindex];
                DoBruteForceCopy = FALSE;
            } break;
            case FieldDataType_Bit : break; /* handle below */
            default : CHECK(FALSE); /* Future work: when more complex types are added, remove this CHECK */
        }
    }

    if (DoBruteForceCopy)
    {
        double val = GetFieldValue(src, srcindex);
        SetFieldValue(dst, dstindex, val);
    }
} /* CopyFieldValue() */


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined TECPLOTKERNEL
#endif /* TECPLOTKERNEL */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


/*
 */
void SetFieldDataPtrToAllZeros(FieldData_pa fd)
{
    REQUIRE(VALID_REF(fd));
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

    LgIndex_t NumValues = GetFieldDataNumValues(fd);

    /*
     * memset each byte to 0 works for floats and doubles and works regardless
     * of byte ordering or alignment.
     */
    size_t NumBytesToMemSet = 0;
    if (IsFieldDataDirectAccessAllowed(fd))
    {
        switch (GetFieldDataType(fd))
        {
            case FieldDataType_Int64 : CHECK(FALSE); /* Future work: remove CHECK */
            case FieldDataType_Double :
            {
                CHECK(sizeof(UInt64_t) == 8 && sizeof(double) == 8);
                NumBytesToMemSet = NumValues * sizeof(Int64_t);
            } break;
            case FieldDataType_Int32 :
            case FieldDataType_Float :
            {
                CHECK(sizeof(UInt32_t) == 4 && sizeof(float) == 4);
                NumBytesToMemSet = NumValues * sizeof(Int32_t);
            } break;
            case FieldDataType_Int16 :
            {
                CHECK(sizeof(UInt16_t) == 2);
                NumBytesToMemSet = NumValues * sizeof(Int16_t);
            } break;
            case FieldDataType_Byte :
            {
                NumBytesToMemSet = NumValues * sizeof(Byte_t);
            } break;
            case FieldDataType_Bit :
            {
                NumBytesToMemSet = ((NumValues + 7) / 8) * sizeof(Byte_t);
            } break;
            default :
            {
                CHECK(FALSE);
            } break;
        }
    }

    if (NumBytesToMemSet > 0)
    {
        char* fd_data = (char*)GetFieldDataVoidPtr(fd);
        #if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
        #else
            memset(fd_data, 0, NumBytesToMemSet);
        #endif
    }
    else
    {
        int ii;
        for (ii = 0; ii < NumValues; ii++)
            SetFieldValue(fd, ii, 0.0);
    }

} /* SetFieldDataPtrToAllZeros() */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */
