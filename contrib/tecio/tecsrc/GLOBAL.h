/* BEGINREMOVEFROMADDON */
/* NOTE: All code contained between comments that look like
 *             BEGINREMOVEFROMADDON
 *             ENDREMOVEFROMADDON
 * are pulled out to create the GLOBAL.h file used in addons.
 */
/* ENDREMOVEFROMADDON */

/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2010 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/

#ifndef _GLOBAL_H
#define _GLOBAL_H

#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_MAINMODULE && defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
#define EXTERN extern
#endif

#define EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/* BEGINREMOVEFROMADDON */
/*
 * The reason for wrapping this test with "begin and end remove from addon" key
 * words is so that the ADK users doesn't have to see this mess.
 */
#if !defined COREAPI && \
    !defined TECUTILMMODULE && \
    !defined TECUTILOMODULE && \
    !defined TECUTILQMODULE && \
    !defined TECUTILSMODULE
/* we don't want Tecplot internals using deprecated interfaces */
#  undef EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
#endif
/* ENDREMOVEFROMADDON */


/****************************************************************
 *                                                              *
 *                          MACROS                              *
 *                                                              *
 ****************************************************************/
#if defined TRUE
#undef TRUE
#endif
#if defined FALSE
#undef FALSE
#endif
#if defined MIN
#undef MIN
#endif
#if defined MAX
#undef MAX
#endif
#if defined ROUND
#undef ROUND
#endif
#if defined ROUND2
#undef ROUND2
#endif
#if defined TRUNC
#undef TRUNC
#endif

#define TRUE                  ((Boolean_t)1)
#define FALSE                 ((Boolean_t)0)

/****************************************************************
 *                                                              *
 *                           MACROS                             *
 *                                                              *
 ****************************************************************/
#define ABS(X)                ((X) >= 0 ? (X) : -(X) )
#define MAX(X,Y)              ((X) > (Y) ? (X) : (Y) )
#define MIN(X,Y)              ((X) < (Y) ? (X) : (Y) )
#define BESTSHOWCOLOR(X)      ((X) == White_C ? Black_C : White_C)
#define ROUND_TO_BYTE(X)      ((BYTE)((X)+0.499))
#define ROUNDS(X)             ((short)((X)+0.499))
#define ROUNDL(X)             ((LgIndex_t)((X)+0.499))
#define ROUND2(X)             ((X) >= 0 ? ((int)((X)+0.499)) : ((int)((X)-0.499)))
#define TRUNC(X)              ((short) (X))
#define RAD_TO_DEG(rad)       (180.*(rad)/PI)
#define DEG_TO_RAD(deg)       (PI*(deg)/180.)

# define CAPITAL(C) ( ('a'<=(C)&&(C)<='z') ? ((C)+('A'-'a')) : (C) ) /* okay for UNICODE */

#include "TASSERT.h"

#if defined TECPLOTKERNEL && defined MSWIN
/* CORE SOURCE CODE REMOVED */
#else
#define ISEMPTYSTRING(S)      ( ((const char*)(S))[0] == '\0' )
#endif

#define ISWHITESPACE(C)       ((C == ' ') || (C == '\t') || (C == '\n'))
#define ISSEPARATOR(C)        ((C == ' ') || (C == '\t') || (C == ','))
/* clamp the input to the specified range */
#define CLAMP(value,low,high) ((value)<(low) ? (low) : (value) > (high) ? (high) : (value))
/* integer division rounds any fraction up (for example n=16,d=3 results in 6) */
#define INTEGER_DIVIDE_AND_ROUND_UP(n, d) (((int)(n)+(int)(d)-1)/(int)(d))

/* BEGINREMOVEFROMADDON */
/**
 * Calcualtes the cell's primary corner or cell centered index from the I, J,
 * and K indices.
 *
 * Consider this IJ zone dimensioned 4 by 3:
 * @verbatim
    +-------+-------+-------+-------+
    |       |       |       |       |
    |  <8>  |  <9>  |  <10> |  <11> |  <--- ghost cells
    |       |       |       |       |
    |8      |9      |10     |11     |
    +-------+-------+-------+-------+
    |       |       |       |       |
    |  <4>  |  <5>  |  <6>  |  <7>  |
    |       |       |       |       |
    |4      |5      |6      |7      |
    +-------+-------+-------+-------+
    |       |       |       |       |
    |  <0>  |  <1>  |  <2>  |  <3>  |
    |       |       |       |       |
    |0      |1      |2      |3      |
    +-------+-------+-------+-------+
                                 .
                                /|\
                                 |
                                 |
                            ghost cells
@endverbatim
 */
#define IJKINDEX(CZData,I,J,K) ((I) + \
                                ((J)*(CZData)->NumIPts) + \
                                ((K)*(CZData)->NumIJPts))

/**
 * Calculates the I indice from the cell's primary corner or cell centered
 * index. See IJKINDEX() for a picture.
 */
#define IINDEX(CZData,N) ((N) % (CZData)->NumIPts)

/**
 * Calculates the J indice from the cell's primary corner or cell centered
 * index. See IJKINDEX() for a picture.
 */
#define JINDEX(CZData,N) (((N) % (CZData)->NumIJPts)/(CZData)->NumIPts)

/**
 * Calculates the K indice from the cell's primary corner or cell centered
 * index. See IJKINDEX() for a picture.
 */
#define KINDEX(CZData,N) ((N)/(CZData)->NumIJPts)
/* ENDREMOVEFROMADDON */

/* */
#define SWITCH(Type,A,B)      do {Type T = (A); (A) = (B); (B) = T;} while (FALSE)
#define SWITCH_DOUBLES(A,B)   SWITCH(double, (A), (B))
#define FPRINTFOK(x)          (Boolean_t)((x) > 0)
#define GRAPHICSARE3D(F)      ((F->PlotType == PlotType_Cartesian3D))

/* convenience macros for implication, P -> Q, and equivalence, P <-> Q. */
#define IMPLICATION(P,Q) (!(P) || (Q))
#define EQUIVALENCE(P,Q) ((P) == (Q))

/* suppress compiler warnings about unused parameters */
#if defined UNUSED
#undef UNUSED
#endif
#define UNUSED(param) (void)param

/**
 * Converts a double into a float value
 *
 * param val
 *     double value to be converted
 */
#define CONVERT_DOUBLE_TO_FLOAT(val) \
  ( (val) >= SMALLFLOAT \
    ? ( (val) < LARGEFLOAT \
        ? (float)(val) \
        : (float)LARGEFLOAT \
      ) \
    : ( (val) <= -SMALLFLOAT  \
        ? ( (val) > -LARGEFLOAT \
            ? (float)(val) \
            : (float)-LARGEFLOAT \
          ) \
        : (float)0.0 \
      ) \
  )


/**
 * Clamps a double at the limits of Tecplot's precision
 *
 * param val
 *     double value to be clamped
 */
#define CLAMP_DOUBLE(val) \
  ( (val) >= SMALLDOUBLE \
    ? ( (val) < LARGEDOUBLE \
        ? (double)(val) \
        : (double)LARGEDOUBLE \
      ) \
    : ( (val) <= -SMALLDOUBLE  \
        ? ( (val) > -LARGEDOUBLE \
            ? (double)(val) \
            : (double)-LARGEDOUBLE \
          ) \
        : (double)0.0 \
      ) \
  )


/**
 * Converts a double into a 4-byte (signed) integer value
 *
 * param val
 *     double value to be converted
 */
#define CONVERT_DOUBLE_TO_INT32(val) \
  ( (val) >= 1.0 \
    ? ( (val) < MAXINT32 \
        ? (Int32_t)(val) \
        : (Int32_t)MAXINT32 \
      ) \
    : ( (val) <= -1.0  \
        ? ( (val) > (Int32_t)-MAXINT32 \
            ? (Int32_t)(val) \
            : (Int32_t)-MAXINT32 \
          ) \
        : (Int32_t)0.0 \
      ) \
  )


/**
 * Converts a double into a 2-byte (signed) integer value
 *
 * param val
 *     double value to be converted
 */
#define CONVERT_DOUBLE_TO_INT16(val) \
  ( (val) >= 1.0 \
    ? ( (val) < MAXINT16 \
        ? (Int16_t)(val) \
        : (Int16_t)MAXINT16 \
      ) \
    : ( (val) <= -1.0  \
        ? ( (val) > (Int16_t)-MAXINT16 \
            ? (Int16_t)(val) \
            : (Int16_t)-MAXINT16 \
          ) \
        : (Int16_t)0.0 \
      ) \
  )

/**
 * Copies two bytes from SrcBuffer to DstBuffer without causing a page
 * fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the two bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the two bytes from
 */
#define COPY_2_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[0]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[1]; \
        } while (FALSE)

/**
 * Copies two bytes from SrcBuffer to DstBuffer swapping the bytes
 * as it copies.  Will not cause a page fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the two bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the two bytes from
 */
#define COPY_AND_REVERSE_2_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[1]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[0]; \
        } while (FALSE)

/**
 * Copies four bytes from SrcBuffer to DstBuffer without causing a page
 * fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the four bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the four bytes from
 */
#define COPY_4_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[0]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[1]; \
          ((Byte_t *)(DstBuffer))[2] = ((Byte_t *)(SrcBuffer))[2]; \
          ((Byte_t *)(DstBuffer))[3] = ((Byte_t *)(SrcBuffer))[3]; \
        } while (FALSE)

/**
 * Copies four bytes from SrcBuffer to DstBuffer swapping the bytes
 * as it copies.  Will not cause a page fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the four bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the four bytes from
 */
#define COPY_AND_REVERSE_4_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[3]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[2]; \
          ((Byte_t *)(DstBuffer))[2] = ((Byte_t *)(SrcBuffer))[1]; \
          ((Byte_t *)(DstBuffer))[3] = ((Byte_t *)(SrcBuffer))[0]; \
        } while (FALSE)

/**
 * Copies four bytes from SrcBuffer to DstBuffer without causing a page
 * fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the four bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the four bytes from
 */
#define COPY_8_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[0]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[1]; \
          ((Byte_t *)(DstBuffer))[2] = ((Byte_t *)(SrcBuffer))[2]; \
          ((Byte_t *)(DstBuffer))[3] = ((Byte_t *)(SrcBuffer))[3]; \
          ((Byte_t *)(DstBuffer))[4] = ((Byte_t *)(SrcBuffer))[4]; \
          ((Byte_t *)(DstBuffer))[5] = ((Byte_t *)(SrcBuffer))[5]; \
          ((Byte_t *)(DstBuffer))[6] = ((Byte_t *)(SrcBuffer))[6]; \
          ((Byte_t *)(DstBuffer))[7] = ((Byte_t *)(SrcBuffer))[7]; \
        } while (FALSE)

/**
 * Copies eight bytes from SrcBuffer to DstBuffer swapping the bytes
 * as it copies.  Will not cause a page fault due to misaligned words.
 *
 * param DstBuffer
 *     Pointer the buffer to send the four bytes to
 * param SrcBuffer
 *     Pointer the buffer to get the four bytes from
 */
#define COPY_AND_REVERSE_8_UNALIGNED_BYTES(DstBuffer, SrcBuffer) \
        do { \
          /* cannot check sizeof(SrcBuffer) or sizeof(DstBuffer) because they are */ \
          /* most likely single byte pointers into unaligned blocks of data */ \
          ((Byte_t *)(DstBuffer))[0] = ((Byte_t *)(SrcBuffer))[7]; \
          ((Byte_t *)(DstBuffer))[1] = ((Byte_t *)(SrcBuffer))[6]; \
          ((Byte_t *)(DstBuffer))[2] = ((Byte_t *)(SrcBuffer))[5]; \
          ((Byte_t *)(DstBuffer))[3] = ((Byte_t *)(SrcBuffer))[4]; \
          ((Byte_t *)(DstBuffer))[4] = ((Byte_t *)(SrcBuffer))[3]; \
          ((Byte_t *)(DstBuffer))[5] = ((Byte_t *)(SrcBuffer))[2]; \
          ((Byte_t *)(DstBuffer))[6] = ((Byte_t *)(SrcBuffer))[1]; \
          ((Byte_t *)(DstBuffer))[7] = ((Byte_t *)(SrcBuffer))[0]; \
        } while (FALSE)

/**
 * Reverses the byte order of the specified 2 byte buffer.
 *
 * param Buffer
 *     Pointer to the 2 bytes needing byte order reversal.
 */
#define REVERSE_2_BYTES_1_AT_A_TIME(Buffer) \
          do { \
            Byte_t Byte0 = ((Byte_t *)(Buffer))[0]; \
            CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==2); \
            ((Byte_t *)(Buffer))[0] = ((Byte_t *)(Buffer))[1]; \
            ((Byte_t *)(Buffer))[1] = Byte0; \
          } while (FALSE)

#define REVERSE_2_BYTES_2_AT_A_TIME(Buffer) \
          do { \
            UInt16_t data_bits = ((UInt16_t *)(Buffer))[0]; \
            CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==2); \
            ((UInt16_t *)(Buffer))[0] = (((data_bits)<<8) | \
                                         ((data_bits&0xff))); \
          } while (FALSE)

/* REVERSE_2_BYTES_2_AT_A_TIME may actually be slower, needs testing. */
#define REVERSE_2_BYTES REVERSE_2_BYTES_1_AT_A_TIME

/**
 * Reverses the byte order of the specified 4 byte buffer.
 *
 * param Buffer
 *     Pointer to the 4 bytes needing byte order reversal.
 *
 * How this works:
 *
 *   ABCD
 *   D--- <<24  (1)
 *
 *   ABCD
 *   --C- &0x0000ff00
 *   -C-- <<8   (2)
 *
 *   ABCD
 *   -B-- &0x00ff0000
 *   --B- >>8   (3)
 *
 *   ABCD
 *   ---A >>24  (4)
 *
 * (1) | (2) | (3) | (4) = DCBA.
 *
 */
#define REVERSE_4_BYTES_1_AT_A_TIME(Buffer) \
          do { \
            Byte_t Byte0 = ((Byte_t *)(Buffer))[0]; \
            Byte_t Byte1 = ((Byte_t *)(Buffer))[1]; \
            CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==4); \
            ((Byte_t *)(Buffer))[0] = ((Byte_t *)(Buffer))[3]; \
            ((Byte_t *)(Buffer))[1] = ((Byte_t *)(Buffer))[2]; \
            ((Byte_t *)(Buffer))[2] = Byte1; \
            ((Byte_t *)(Buffer))[3] = Byte0; \
          } while (FALSE)

#define REVERSE_4_BYTES_4_AT_A_TIME(Buffer) \
          do { \
            UInt32_t data_bits = *((UInt32_t *)(Buffer)); \
            CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==4); \
            *((UInt32_t *)(Buffer)) = (((data_bits)<<24)            | \
                                       ((data_bits&0x0000ff00)<<8)  | \
                                       ((data_bits&0x00ff0000)>>8)  | \
                                       ((data_bits)>>24)); \
          } while (FALSE)

#if defined MSWIN
/*
 * The DevStuido compiler seems to be the only one that can truly handle this
 * when optimization is turned on.
 */
#define REVERSE_4_BYTES REVERSE_4_BYTES_4_AT_A_TIME
#else
#define REVERSE_4_BYTES REVERSE_4_BYTES_1_AT_A_TIME
#endif

/**
 * Reverses the byte order of the specified 8 byte buffer.
 *
 * param Buffer
 *     Pointer to the 8 bytes needing byte order reversal.
 */
#define REVERSE_8_BYTES_1_AT_A_TIME(Buffer) \
        do { \
            Byte_t Byte0 = ((Byte_t *)(Buffer))[0]; \
            Byte_t Byte1 = ((Byte_t *)(Buffer))[1]; \
            Byte_t Byte2 = ((Byte_t *)(Buffer))[2]; \
            Byte_t Byte3 = ((Byte_t *)(Buffer))[3]; \
            CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==8); \
            ((Byte_t *)(Buffer))[0] = ((Byte_t *)(Buffer))[7]; \
            ((Byte_t *)(Buffer))[1] = ((Byte_t *)(Buffer))[6]; \
            ((Byte_t *)(Buffer))[2] = ((Byte_t *)(Buffer))[5]; \
            ((Byte_t *)(Buffer))[3] = ((Byte_t *)(Buffer))[4]; \
            ((Byte_t *)(Buffer))[4] = Byte3; \
            ((Byte_t *)(Buffer))[5] = Byte2; \
            ((Byte_t *)(Buffer))[6] = Byte1; \
            ((Byte_t *)(Buffer))[7] = Byte0; \
        } while (FALSE)

#define REVERSE_8_BYTES_2_AT_A_TIME(Buffer) \
        do { \
          UInt16_t data_bits_0 = ((UInt16_t *)(Buffer))[0]; \
          UInt16_t data_bits_1 = ((UInt16_t *)(Buffer))[1]; \
          UInt16_t data_bits_2 = ((UInt16_t *)(Buffer))[2]; \
          UInt16_t data_bits_3 = ((UInt16_t *)(Buffer))[3]; \
          CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==8); \
          ((UInt16_t *)(Buffer))[0] = (((data_bits_3)<<8) | \
                                       ((data_bits_3&0xff))); \
          ((UInt16_t *)(Buffer))[1] = (((data_bits_2)<<8) | \
                                       ((data_bits_2&0xff))); \
          ((UInt16_t *)(Buffer))[2] = (((data_bits_1)<<8) | \
                                       ((data_bits_1&0xff))); \
          ((UInt16_t *)(Buffer))[3] = (((data_bits_0)<<8) | \
                                       ((data_bits_0&0xff))); \
        } while (FALSE)

#define REVERSE_8_BYTES_4_AT_A_TIME(Buffer) \
        do { \
          UInt32_t data_bits_0 = ((UInt32_t *)(Buffer))[0]; \
          UInt32_t data_bits_1 = ((UInt32_t *)(Buffer))[1]; \
          CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==8); \
          ((UInt32_t *)(Buffer))[0] = (((data_bits_1)<<24)           | \
                                       ((data_bits_1&0x0000ff00)<<8) | \
                                       ((data_bits_1&0x00ff0000)>>8) | \
                                       ((data_bits_1)>>24)); \
          ((UInt32_t *)(Buffer))[1] = (((data_bits_0)<<24)           | \
                                       ((data_bits_0&0x0000ff00)<<8) | \
                                       ((data_bits_0&0x00ff0000)>>8) | \
                                       ((data_bits_0)>>24)); \
        } while (FALSE)

#define REVERSE_8_BYTES_8_AT_A_TIME(Buffer) \
        do { \
          UInt64_t data_bits = *((UInt64_t *)(Buffer)); \
          CHECK(sizeof(*(Buffer))==1 || sizeof(*(Buffer))==8); \
          *((UInt64_t *)(Buffer)) = (((data_bits)<<56) | \
                                     ((data_bits&0x000000000000ff00)<<40) | \
                                     ((data_bits&0x0000000000ff0000)<<24) | \
                                     ((data_bits&0x00000000ff000000)<<8)  | \
                                     ((data_bits&0x000000ff00000000)>>8)  | \
                                     ((data_bits&0x0000ff0000000000)>>24) | \
                                     ((data_bits&0x00ff000000000000)>>40) | \
                                     ((data_bits)>>56)); \
        } while (FALSE)


#if defined MSWIN
/*
 * The DevStuido compiler seems to be the only one that can truly handle this
 * when optimization is turned on.
 */
#define REVERSE_8_BYTES REVERSE_8_BYTES_4_AT_A_TIME
#else
#define REVERSE_8_BYTES REVERSE_8_BYTES_1_AT_A_TIME
#endif


/****************************************************************
 *                                                              *
 *             ADD-ON MSWIN IMPORT/EXPORT DEFINITIONS            *
 *                                                              *
 ****************************************************************/
#if defined MSWIN
#  define STDCALL __stdcall
#else
#  define STDCALL
#endif /* MSWIN */

#if defined (__cplusplus)
# define EXTERNC extern "C"
#else
# define EXTERNC
#endif /* __cplusplus */

#if defined MSWIN
#if defined AMTEC_INTERNAL_MAKELIBTEC || defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#  else
#    define TECPLOT_DLLAPI _declspec ( dllimport )
#  endif
#else
#  define TECPLOT_DLLAPI
#endif

#define LINKTOADDON EXTERNC TECPLOT_DLLAPI


/*
 *
 * Usage:
 * EXPORTFROMADDON void STDCALL InitTecAddOn(void) { ... }
 *
 */
#if defined MSWIN
# define EXPORTFROMADDON EXTERNC _declspec ( dllexport )
#else
# define EXPORTFROMADDON EXTERNC
#endif /* MSWIN */

#define EXPORTFROMDLL EXPORTFROMADDON

#define InitTecAddOn           InitTecAddOn113
#define TEC_INIT_FUNCTION_NAME "InitTecAddOn113"

/* BEGINREMOVEFROMADDON */
/* Use INLINE for static functions that could be optimized as inline. */
#if defined (__cplusplus) && !defined _DEBUG
# define INLINE inline
#else
# define INLINE static
#endif /* __cplusplus */
/* ENDREMOVEFROMADDON */


/* BEGINREMOVEFROMADDON */
#if defined (MSWIN) ||\
    defined (INTERX) ||\
    defined (LINUX) ||\
    defined (SUNSOLARIS86X) ||\
    defined (COMPAQALPHA) ||\
    defined (DEC) ||\
    defined (__LITTLE_ENDIAN__)
#define MACHINE_DOES_INTEL_ORDER
#endif

#if defined( MACHINE_DOES_INTEL_ORDER )
# define SwapBytes(IntelOrder) (!(IntelOrder))
#else
# define SwapBytes(IntelOrder) (IntelOrder)
#endif
/* ENDREMOVEFROMADDON */

#if defined DECALPHA   || \
    defined LINUXALPHA || \
    defined LINUX64    || \
    defined MAC64      || \
    defined IBMRS6000  || \
    defined SUN        || \
    defined HP         || \
    defined COMPAQALPHA
#define LONGIS64
#endif

/****************************************************************
 *                                                              *
 *                       HARD CONSTANTS                         *
 *                                                              *
 ****************************************************************/
#define LARGEMEMORY              ((size_t)-1)

/* BEGINREMOVEFROMADDON */
/* Tclinterp add-on barfs on these huge integer constants */
/* Note: Tecplot is conservative by one on LARGEINTs max */
#define LARGEINT64               9223372036854775806LL
/* ENDREMOVEFROMADDON */
#define LARGEINT32               2147483646
#define LARGEINT16               32766
#define LARGEINT8                126

/* BEGINREMOVEFROMADDON */
#define LARGEUINT64              18446744073709551614ULL
/* ENDREMOVEFROMADDON */
#define LARGEUINT32              4294967294U
#define LARGEUINT16              65534U
#define LARGEUINT8               254U

#ifdef INDEX_16_BIT
#define MAXINDEX               ((LgIndex_t)LARGEINT16)
#else
#define MAXINDEX               ((LgIndex_t)LARGEINT32)
#endif
#define MAXZONEMAP               MAXINDEX
#define LARGEDOUBLE              1.0e+150
#define SMALLDOUBLE              1.0e-150
#define LARGESTEXPONENT          150
#define SMALLESTEXPONENT         -150

#define SMALLESTDOUBLE           SMALLDOUBLE

#define LARGESTDOUBLEEXPONENT    308
#define SMALLESTDOUBLEEXPONENT   -307
#define LARGESTDOUBLE            1.0e+308
#define LARGEFLOAT               3.40282347E+38
#define SMALLFLOAT               1.17549435E-38
#define SMALLSTDOUBLE            1.0e-307

/* Visual Studio 2008 defines MAXINT32, MAXINT16 which collide with ours */
#if defined MAXINT32
#undef MAXINT32
#endif
#if defined MAXINT16
#undef MAXINT16
#endif

#define MAXINT32                 LARGEINT32
#define MAXINT16                 LARGEINT16
#define ETX                      3
#define LN2                      0.69314718055994530942
#define LN10                     2.30258509299404568402
#define PIOVER2                  1.57079632679489661923
#define TWOPI                    6.28318530717958647692
#if defined PI
#undef PI
#endif
#define PI                       3.14159265358979323846
#define ANGLEEPSILON             1.0e-10
#define LARGESTANGLE             (4*PI+ANGLEEPSILON)
#define DEGPERRADIANS            57.295779513082323
#define CMPERINCH                2.54
#define POINTSPERINCH            72.0
#define FONTMOVEMARK             192
#define FONTDECISIONMARK         128
#define FONTLINEMARK             64
#define BAD_SET_VALUE            ((SetIndex_t)-1)
#define MENU_POSITION_FIRST      (0)
#define MENU_POSITION_LAST       (-1)
#define INVALID_UNIQUE_ID        0

#define BADSETVALUE              BAD_SET_VALUE
#define SOLID_TRANSLUCENCY       0
#define BAD_DISTANCE             (-1.0)
/* MIN_CIRCUMFERENTIAL_INDEX is the min J dimension for circular zones */
#define MIN_CIRCUMFERENTIAL_INDEX  4

#define VALID_STRAND_ID(StrandID) (0 <= (StrandID) && (StrandID) < MAXZONEMAP)
#define STRAND_ID_STATIC          (-1)
#define STRAND_ID_PENDING         (-2)

/*
 * Need 3 passes for "Rest of pie" method but can only use 3 clip planes
 * Need only 1 pass for "Piece of pie" method and can use 6 clip planes
*/
#define MAX_ALLOWABLE_CLIPPASSES 1
#define MAX_ALLOWABLE_CLIPPLANES 6
#define INVALID_CLIP_PLANE -1
#define VALID_CLIP_PLANE(clipPlane) (0 <= clipPlane && clipPlane < MAX_ALLOWABLE_CLIPPLANES)

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined _DEBUG
#else
#endif
#if 0 /* NOTUSED */
#endif
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */


/*
 * NOTE: If you change TecplotBinaryFileVersion, you MUST also:
 *
 * 1. Update preplot:
 *    - Change this define symbol in preplot.cpp
 *    - Change version number in the data file format in the comments in preplot.cpp
 *    - Change the version number of Preplot itself in preplot.cpp
 * 2. Maintain the ability to write the old plt file format:
 *    - Add a new entry to BinaryFileVersion_e
 *    - Add a concrete class of the VersionWriterInterface, and update
 *      VersionWriterAbstractFactory to return the correct instance for the previous and
 *      new BinaryFileVersion_e
 *    - Abstract away the difference in the two versions behind an interface (if one does
 *      not yet exist) and create concrete implementations that can write the old and the
 *      new versions. For a trivial example of this, see FileTypeWriterInterface and its
 *      associated factory and concrete classes.
 */
#define TecplotBinaryFileVersion    112 /* NOTE: only change this when we change the binary file format */

#define    MaxNumZonesOrVars           MAXZONEMAP
#define    MaxXAxes                    5
#define    MaxYAxes                    5
#define    MaxGeoSegments              50
#define    MaxPtsCircleOrEllipse       720
#define    MaxFrames                   2048
#define    MaxCustomLabelSets          10
#define    MaxFontMoves                20000
#define    MaxColorMapOverrides        16
#define    MaxValueBlankConstraints    8
#define    MaxContourGroups            8
#define    MaxIsoSurfaceGroups         8
#define    MaxIsoSurfaceSpecificLevels 3
#define    MaxSliceGroups              8

#define    MaxColorMapGroups         8
#define    DefaultNumContLevels      15


#define    DefaultColorMapGroup      ((SmInteger_t)0)
#define    BADGROUPNUMBER            ((SmInteger_t)-1)
#define    UNUSEDGROUPNUMBER         ((SmInteger_t)0)

#define VALID_ISOSURFACE_GROUP(Group) (((((SmInteger_t)Group) >= 0) && (((SmInteger_t)Group) < MaxIsoSurfaceGroups)))
#define VALID_SLICE_GROUP(Group)      (((((SmInteger_t)Group) >= 0) && (((SmInteger_t)Group) < MaxSliceGroups)))
#define VALID_COLORMAP_GROUP(Group)   (((((SmInteger_t)Group) >= 0) && (((SmInteger_t)Group) < MaxColorMapGroups)))

#define    MAX_AUTO_COLOR_SEQUENCE_VALUES  6


/*
 * If any of these values changes its corresponding value in preplot.c must
 * change to match it so that files created by preplot and Tecplot are
 * consistent.
 */
#define    MaxChrsDatasetTitle       256
#define    MaxChrsZnTitle            128
#define    MaxChrsVarName            128
#define    MaxChrsZnOrVarName        128
/* currently limited to MaxLineIndex in preplot.c */
#define    MaxChrsAuxValueString     32000

#define    MaxNumViews               16
#define    MaxBasicSizes             5
#define    MaxColorMapControlPoints  50
#define    MaxRawColorMapEntries     800
#define    MaxDataSetReaders         100
#define    MaxExtendedCurveFits      100
#define    MaxColorMapCycles         20


/* Dimension Limits */

#define    MinPaperDimInWorkArea     0.5
#define    MinFrameWidth             0.25
#define    MinFrameHeight            0.25
#define    MinAxisLength             0.1


#define    BadEnumValue              255

/* BEGINREMOVEFROMADDON */
/* define class element limits */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif
/* ENDREMOVEFROMADDON */

/*
 * Floating point values are written to layouts with a certain precision.
 * A high precision is necessary in some cases (like streamtrace starting locations)
 * This used to be set to 12 which was not high enough.   It is now set to 16 which
 * appears to be sufficient.   This also seems to jive with the number of digits of
 * precision that are found in "ieee double precision" values which is 53 bits or
 * equivalent to approximately 16 digits. -bdp
 *
 */
#define STYLE_FLOAT_PRECISION 16


/*
 * Auxiliary data common names.
 *
 *      Define Name                                 Data Name                               Data Type    Data Location
 *      ------------------------------------------  ------------------------------------    ---------    -------------
 */
#define AuxData_Common_Incompressible               "Common.Incompressible"              /* Boolean_t    Dataset */
#define AuxData_Common_Density                      "Common.Density"                     /* double       Dataset */
#define AuxData_Common_SpecificHeat                 "Common.SpecificHeat"                /* double       Dataset */
#define AuxData_Common_SpecificHeatVar              "Common.SpecificHeatVar"             /* int          Dataset */
#define AuxData_Common_GasConstant                  "Common.GasConstant"                 /* double       Dataset */
#define AuxData_Common_GasConstantVar               "Common.GasConstantVar"              /* int          Dataset */
#define AuxData_Common_Gamma                        "Common.Gamma"                       /* double       Dataset */
#define AuxData_Common_GammaVar                     "Common.GammaVar"                    /* int          Dataset */
#define AuxData_Common_Viscosity                    "Common.Viscosity"                   /* double       Dataset */
#define AuxData_Common_ViscosityVar                 "Common.ViscosityVar"                /* int          Dataset */
#define AuxData_Common_Conductivity                 "Common.Conductivity"                /* double       Dataset */
#define AuxData_Common_ConductivityVar              "Common.ConductivityVar"             /* int          Dataset */
#define AuxData_Common_AngleOfAttack                "Common.AngleOfAttack"               /* double       Dataset */
#define AuxData_Common_SpeedOfSound                 "Common.SpeedOfSound"                /* double       Dataset */
#define AuxData_Common_ReferenceU                   "Common.ReferenceU"                  /* double       Dataset */
#define AuxData_Common_ReferenceV                   "Common.ReferenceV"                  /* double       Dataset */
#define AuxData_Common_XVar                         "Common.XVar"                        /* int          Dataset */
#define AuxData_Common_YVar                         "Common.YVar"                        /* int          Dataset */
#define AuxData_Common_ZVar                         "Common.ZVar"                        /* int          Dataset */
#define AuxData_Common_CVar                         "Common.CVar"                        /* int          Dataset */
#define AuxData_Common_UVar                         "Common.UVar"                        /* int          Dataset */
#define AuxData_Common_VVar                         "Common.VVar"                        /* int          Dataset */
#define AuxData_Common_WVar                         "Common.WVar"                        /* int          Dataset */
#define AuxData_Common_VectorVarsAreVelocity        "Common.VectorVarsAreVelocity"       /* Boolean_t    Dataset */
#define AuxData_Common_PressureVar                  "Common.PressureVar"                 /* int          Dataset */
#define AuxData_Common_TemperatureVar               "Common.TemperatureVar"              /* int          Dataset */
#define AuxData_Common_DensityVar                   "Common.DensityVar"                  /* int          Dataset */
#define AuxData_Common_StagnationEnergyVar          "Common.StagnationEnergyVar"         /* int          Dataset */
#define AuxData_Common_MachNumberVar                "Common.MachNumberVar"               /* int          Dataset */
#define AuxData_Common_ReferenceMachNumber          "Common.ReferenceMachNumber"         /* double       Dataset */
#define AuxData_Common_ReferenceW                   "Common.ReferenceW"                  /* double       Dataset */
#define AuxData_Common_PrandtlNumber                "Common.PrandtlNumber"               /* double       DataSet */
#define AuxData_Common_Axisymmetric                 "Common.Axisymmetric"                /* Boolean_t    Dataset */
#define AuxData_Common_AxisOfSymmetryVarAssignment  "Common.AxisOfSymmetryVarAssignment" /* int          Dataset */
#define AuxData_Common_AxisValue                    "Common.AxisValue"                   /* double       Dataset */
#define AuxData_Common_SteadyState                  "Common.SteadyState"                 /* Boolean_t    Dataset */
#define AuxData_Common_TurbulentKineticEnergyVar    "Common.TurbulentKineticEnergyVar"   /* int          Dataset */
#define AuxData_Common_TurbulentDissipationRateVar  "Common.TurbulentDissipationRateVar" /* int          Dataset */
#define AuxData_Common_TurbulentViscosityVar        "Common.TurbulentViscosityVar"       /* int          Dataset */
#define AuxData_Common_TurbulentFrequencyVar        "Common.TurbulentFrequencyVar"       /* int          Dataset */
#define AuxData_Common_Gravity                      "Common.Gravity"                     /* double       Dataset */
#define AuxData_Common_IsBoundaryZone               "Common.IsBoundaryZone"              /* Boolean_t    Zone */
#define AuxData_Common_BoundaryCondition            "Common.BoundaryCondition"           /* BCondition   Zone */
#define AuxData_Common_Time                         "Common.Time"                        /* double       Zone */
#define AuxData_Common_Mean                         "Common.Mean"                        /* double       Variable */
#define AuxData_Common_Median                       "Common.Median"                      /* double       Variable */
#define AuxData_Common_Variance                     "Common.Variance"                    /* double       Variable */
#define AuxData_Common_StdDev                       "Common.StdDev"                      /* double       Variable */
#define AuxData_Common_AvgDev                       "Common.AvgDev"                      /* double       Variable */
#define AuxData_Common_GeoMean                      "Common.GeoMean"                     /* double       Variable */
#define AuxData_Common_ChiSqre                      "Common.ChiSqre"                     /* double       Variable */







/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined THREED
#endif
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

/* Tecplot Add-on Custom Products */

/* BEGINREMOVEFROMADDON */
/* In activeX, the color constants are an enum type,
   so the activeX source code parser handles these as
   a special case, and the types do not need to be
   indicated as with the other hard #define constants */
/* ENDREMOVEFROMADDON */

#define    Black_C           ((ColorIndex_t)0)
#define    Red_C             ((ColorIndex_t)1)
#define    Green_C           ((ColorIndex_t)2)
#define    Blue_C            ((ColorIndex_t)3)
#define    Cyan_C            ((ColorIndex_t)4)
#define    Yellow_C          ((ColorIndex_t)5)
#define    Purple_C          ((ColorIndex_t)6)
#define    White_C           ((ColorIndex_t)7)

#define    Custom1_C         ((ColorIndex_t)8)
#define    Custom2_C         ((ColorIndex_t)9)
#define    Custom3_C         ((ColorIndex_t)10)
#define    Custom4_C         ((ColorIndex_t)11)
#define    Custom5_C         ((ColorIndex_t)12)
#define    Custom6_C         ((ColorIndex_t)13)
#define    Custom7_C         ((ColorIndex_t)14)
#define    Custom8_C         ((ColorIndex_t)15)
#define    Custom9_C         ((ColorIndex_t)16)

#define    Custom10_C         ((ColorIndex_t)17)
#define    Custom11_C         ((ColorIndex_t)18)
#define    Custom12_C         ((ColorIndex_t)19)
#define    Custom13_C         ((ColorIndex_t)20)
#define    Custom14_C         ((ColorIndex_t)21)
#define    Custom15_C         ((ColorIndex_t)22)
#define    Custom16_C         ((ColorIndex_t)23)
#define    Custom17_C         ((ColorIndex_t)24)
#define    Custom18_C         ((ColorIndex_t)25)
#define    Custom19_C         ((ColorIndex_t)26)

#define    Custom20_C         ((ColorIndex_t)27)
#define    Custom21_C         ((ColorIndex_t)28)
#define    Custom22_C         ((ColorIndex_t)29)
#define    Custom23_C         ((ColorIndex_t)30)
#define    Custom24_C         ((ColorIndex_t)31)
#define    Custom25_C         ((ColorIndex_t)32)
#define    Custom26_C         ((ColorIndex_t)33)
#define    Custom27_C         ((ColorIndex_t)34)
#define    Custom28_C         ((ColorIndex_t)35)
#define    Custom29_C         ((ColorIndex_t)36)

#define    Custom30_C         ((ColorIndex_t)37)
#define    Custom31_C         ((ColorIndex_t)38)
#define    Custom32_C         ((ColorIndex_t)39)
#define    Custom33_C         ((ColorIndex_t)40)
#define    Custom34_C         ((ColorIndex_t)41)
#define    Custom35_C         ((ColorIndex_t)42)
#define    Custom36_C         ((ColorIndex_t)43)
#define    Custom37_C         ((ColorIndex_t)44)
#define    Custom38_C         ((ColorIndex_t)45)
#define    Custom39_C         ((ColorIndex_t)46)

#define    Custom40_C         ((ColorIndex_t)47)
#define    Custom41_C         ((ColorIndex_t)48)
#define    Custom42_C         ((ColorIndex_t)49)
#define    Custom43_C         ((ColorIndex_t)50)
#define    Custom44_C         ((ColorIndex_t)51)
#define    Custom45_C         ((ColorIndex_t)52)
#define    Custom46_C         ((ColorIndex_t)53)
#define    Custom47_C         ((ColorIndex_t)54)
#define    Custom48_C         ((ColorIndex_t)55)
#define    Custom49_C         ((ColorIndex_t)56)

#define    Custom50_C         ((ColorIndex_t)57)
#define    Custom51_C         ((ColorIndex_t)58)
#define    Custom52_C         ((ColorIndex_t)59)
#define    Custom53_C         ((ColorIndex_t)60)
#define    Custom54_C         ((ColorIndex_t)61)
#define    Custom55_C         ((ColorIndex_t)62)
#define    Custom56_C         ((ColorIndex_t)63)

#define    MultiColor_C      ((ColorIndex_t)(-1))
#define    NoColor_C         ((ColorIndex_t)(-2))
#define    MultiColor2_C     ((ColorIndex_t)(-3))
#define    MultiColor3_C     ((ColorIndex_t)(-4))
#define    MultiColor4_C     ((ColorIndex_t)(-5))
#define    RGBColor_C        ((ColorIndex_t)(-6))
#define    MultiColor5_C     ((ColorIndex_t)(-7))
#define    MultiColor6_C     ((ColorIndex_t)(-8))
#define    MultiColor7_C     ((ColorIndex_t)(-9))
#define    MultiColor8_C     ((ColorIndex_t)(-10))
#define    InvalidColor_C    ((ColorIndex_t)(-255))

#define    FirstCustomColor  Custom1_C
#define    LastCustomColor   Custom56_C
#define    NumCustomColors   (LastCustomColor-FirstCustomColor+1)

#define    FirstBasicColor   Black_C
#define    LastBasicColor    LastCustomColor
#define    NumBasicColors    (LastBasicColor-FirstBasicColor+1)

/* BEGINREMOVEFROMADDON */

/*
 * V8 and earlier used this for MultiColor_C.  We adjust this
 * to the new value in the SetValue layer so old addons work.
 */
#define    OldMultiColor_C   ((ColorIndex_t)255)
/*
 * Gray is only used in the interface for workspace background and
 * for insensitive buttons in Motif.
 * True Black and True White are also interface only.  They draw
 * true black or true white - regardless of what the user has set
 * the RGB values for the black and white basic colors.
 * XOrColor_C is also for interface only.
 */
#define    Gray_C                 (LastBasicColor+1)
#define    DarkGray_C             (LastBasicColor+2) /* Used for inactive frame border color */
#define    XOrColor_C             (LastBasicColor+3)
#define    FirstInterfaceColor    Gray_C
#define    LastInterfaceColor     XOrColor_C

#define    NumInterfaceColors    (LastInterfaceColor-FirstInterfaceColor+1)
#define    NumContourShades      (GeneralBase.Limits.MaxNumContourLevels+1)
#define    NumColorsInColorTable (NumBasicColors+NumInterfaceColors+NumContourShades)
#define    BasicColorOffset      (0)
#define    InterfaceColorOffset  (NumBasicColors)
#define    ContourColorOffset    (NumBasicColors+NumInterfaceColors)

#define    BadKey           (short)31
#define    Plus             (short)43
#define    Minus            (short)45
#define    RetKey           (short)13
#define    DeleteKey        (short)127
#define    ShiftDelete      (short)128
#define    BackSpace        (short)8
#define    LeftArrow        (short)29
#define    RightArrow       (short)30
#define    UpArrow          (short)11
#define    DownArrow        (short)10
#define    Toggle           (short)19
#define    Esc              (short)27
#define    RegFrame         (short)18
#define    DoBitDump        (short)2


/* File Markers */
#define ZoneMarker        299.0
#define GeomMarker        399.0
#define TextMarker        499.0
#define CustomLabelMarker 599.0
#define UserRecMarker     699.0
#define DataSetAuxMarker  799.0
#define VarAuxMarker      899.0
#define EndHeaderMarker   357.0


/*
 * Additional objects that have plotter
 * pens assigned to them.
 */
#define    AxisPen          Custom8_C+1
#define    MajGridPen       Custom8_C+2
#define    MinGridPen       Custom8_C+3
#define    MarkerGridPen    Custom8_C+4
#define    StreamlinePen    Custom8_C+5
#define    ColoredLinePen   Custom8_C+6
#define    BoundaryPen      Custom8_C+7
#define    LabelPen         Custom8_C+8
#define    NumPlotterPens   Custom8_C+9
/* AutoSelectPen will select the correct pen from Black_C thru Custom8_C or ColoredLinePen */
#define    AutoSelectPen    Custom8_C+10
#define    InvalidPen       Custom8_C+99

#define    FirstObjectPen   AxisPen
#define    LastObjectPen    LabelPen

#define    DelZFactor           0.0001

#define    BadBaseValue         NULL


/*
 *   NOTES ON TYPEDEFS:
 *
 *   TYPEDEF TYPE          Suffix
 *   ------------          ------
 *   simple                _t
 *   enumerated            _e
 *   structure             _s
 *   union                 _u
 *   abstract              _a
 *   pointer to simple     _pt
 *   pointer to enumerated _pe
 *   pointer to structure  _ps
 *   pointer to union      _pu
 *   pointer to abstract   _pa
 *   pointer to function   _pf
 *
 *
 *   The only execption is char * typedef's these use _t
 *
 *   Abstract types are intentionally made to be
 *   obscure.  The programmer should not have to know
 *   what the underlying structure really is for abstract
 *   types.
 *
 */


#ifdef MSWIN
# define DIR_SEPARATOR "\\"
#else
# define DIR_SEPARATOR "/"
#endif

/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
/* ENDREMOVEFROMADDON */
    #define TP_FREAD  fread
    #define TP_FWRITE fwrite
/* BEGINREMOVEFROMADDON */
#endif
/* ENDREMOVEFROMADDON */

#if defined MSWIN
#define TP_FFLUSH   fflush
#define TP_FCLOSE   fclose

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
/* ENDREMOVEFROMADDON */
#define TP_UNLINK   remove
#define TP_RMDIR    _rmdir
#define TP_FOPEN    ::fopen
#define TP_FREOPEN  ::freopen
#define TP_STAT     ::_stat
#define TP_GETENV   ::getenv
/* BEGINREMOVEFROMADDON */
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

#if defined _WIN64
#define TP_FSEEK(stream,offset,whence) _fseeki64((stream),(__int64)(offset),(whence))
#define TP_FTELL                       _ftelli64
#else
#define TP_FSEEK(stream, offset, whence) fseek((stream), (long)(offset), (whence))
#define TP_FTELL                         ftell
#endif

#else
#define TP_RMDIR    rmdir
#define TP_UNLINK   unlink
#define TP_FOPEN    fopen
#define TP_FREOPEN  freopen
#define TP_FCLOSE   fclose
#define TP_FFLUSH   fflush
#define TP_FSEEK    fseeko
#define TP_FTELL    ftello
#define TP_STAT     stat
#define _stat       stat /* ...make the UNIXX and MSWIN platforms have the same syntax to use "struct _stat" */
#define TP_GETENV   getenv
#endif

/****************************************************************
 *                                                              *
 *                          SIMPLE TYPEDEFS                     *
 *                                                              *
 ****************************************************************/



/* How to define UInt64_t/Int64_t is platform specific, but they are always 8-bytes */
#if defined MSWIN
typedef    unsigned __int64     UInt64_t;
typedef    __int64              Int64_t;
#else
#if defined CRAY
typedef    unsigned int       UInt64_t;
typedef    int                Int64_t;
#else
#if defined LONGIS64
typedef unsigned long      UInt64_t;
typedef long               Int64_t;
#else
typedef unsigned long long UInt64_t;
typedef long long          Int64_t;
#endif
#endif
#endif

#if defined LONGIS64
typedef    unsigned int    UInt32_t;
typedef    int             Int32_t;
typedef    int             LgInteger_t;
#else
typedef    unsigned int    UInt32_t;
typedef    int             Int32_t;
typedef    int             LgInteger_t;
#endif

typedef    short           Int16_t;
typedef    unsigned short  UInt16_t;
typedef    signed char     Int8_t;
typedef    unsigned char   UInt8_t;

#ifdef INDEX_16_BIT
typedef  Int16_t         LgIndex_t;
#else
typedef  Int32_t         LgIndex_t;
#endif
typedef    LgIndex_t       NodeMap_t;
typedef    LgIndex_t       ScreenDim_t;

/**
 * ArbParam_t type is used for passing arbitrary integers or pointers in
 * parameters. HgIndex_t is used for counting node maps and other things that
 * may individually be LgIndex_t, but in total exceed 32-bit.
 * The general rule is that these are 4 bytes on "32-bit" machines
 * and 8 bytes on "64-bit" machines.
 */
#if defined CRAY
typedef char *ArbParam_t;
typedef long HgIndex_t;
#elif defined LONGIS64
typedef long ArbParam_t;
typedef long HgIndex_t;
#elif defined MSWIN
typedef INT_PTR ArbParam_t;
typedef INT_PTR HgIndex_t;
#else
typedef int ArbParam_t;
typedef int HgIndex_t;
#endif

typedef    ArbParam_t      UniqueID_t;

/* 64 bit offset used to hold file offset and size values. */
typedef Int64_t FileOffset_t;

/**
 * 64 bit offset for memory mapped I/O.
 */
typedef UInt64_t MemMapOffset_t;

/*
 *  SmInteger must be at least a short....
 */

typedef    unsigned char    Byte_t;
typedef    short            SmInteger_t;

/**
 * A number of color index constants are \#defined. These include Black_C, Red_C,
 * Green_C, Blue_C, Cyan_C, Yellow_C, Purple_C, White_C, Custom1_C through
 * Custom56_C, MultiColor_C, NoColor_C, MulitiColor2_C through MulitiColor8_C,
 * RGBColor_C, and InvalidColor_C.  
 */ 
typedef    SmInteger_t      ColorIndex_t;

#ifdef INDEX_16_BIT
typedef  Int16_t          EntIndex_t;
#else
typedef  Int32_t          EntIndex_t;
#endif
typedef    Int16_t          SubZoneIndex_t;

typedef    char             Boolean_t;
typedef    char            *ZoneName_t;
typedef    char            *VarName_t;
typedef    char            *LString_t;

typedef    LgIndex_t        Strand_t;
typedef    LgIndex_t        HeapLength_t;
typedef    LgIndex_t        SegPtsArray_t[MaxGeoSegments];
typedef    double           BasicSize_t[MaxBasicSizes];
typedef    double          *VarList_t;

typedef    HgIndex_t        SetIndex_t;

typedef    unsigned long    SetData_t;
typedef    SetData_t       *SetData_pt;

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

/* The following list identifies items that can be inhibited. */
#define FEATURE_3D                       (1L << 0)
#define FEATURE_3DVOLUME                 (1L << 1)
#define FEATURE_2D                       (1L << 2)
#define FEATURE_XY                       (1L << 3)
#define FEATURE_DATAALTER                (1L << 4)
#define FEATURE_UNSTRUCTUREDDATA         (1L << 5)
#define FEATURE_MULTIPLEFRAMES1          (1L << 6)
#define FEATURE_MULTIPLEZONES1           (1L << 7)
#define FEATURE_MULTIPLEFRAMES5          (1L << 8)
#define FEATURE_MULTIPLEZONES5           (1L << 9)
#define FEATURE_MULTIPLEFRAMES10         (1L << 10)
#define FEATURE_MULTIPLEZONES10          (1L << 11)
#define FEATURE_READNONOEMDATA           (1L << 12) /* Added 07/22/2000 */
#define FEATURE_DATALOADERS              (1L << 13) /* Added 07/22/2000 */
#define FEATURE_DATALOADERS_EXCEPTONE    (1L << 14) /* Added 11/26/2001 */
#define FEATURE_LOADONDEMAND             (1L << 15) /* Added 09/13/2007 */
#define FEATURE_MULTIPLEWORKERTHREADS    (1L << 16) /* Added 09/13/2007 */
#define FEATURE_ISOSURFACEGROUPS         (1L << 17) /* Added 09/21/2007 */
#define FEATURE_SLICEGROUPS              (1L << 18) /* Added 09/21/2007 */
#define FEATURE_STREAMTRACEGROUPS        (1L << 19) /* Added 09/25/2007 not used yet */
#define FEATURE_FEPOLYHEDRON             (1L << 20) /* Added 09/25/2007 */
#define FEATURE_FEPOLYGON                (1L << 21) /* Added 09/25/2007 */

/*
 * KnowFeaturesToInhibit must be updated whenever a new
 * item is added above.
 */
#define KnownFeaturesToInhibit (FEATURE_3D                       |\
                                  FEATURE_3DVOLUME                 |\
                                  FEATURE_2D                       |\
                                  FEATURE_XY                       |\
                                  FEATURE_DATAALTER                |\
                                  FEATURE_UNSTRUCTUREDDATA         |\
                                  FEATURE_MULTIPLEFRAMES1          |\
                                  FEATURE_MULTIPLEZONES1           |\
                                  FEATURE_MULTIPLEFRAMES5          |\
                                  FEATURE_MULTIPLEZONES5           |\
                                  FEATURE_MULTIPLEFRAMES10         |\
                                  FEATURE_MULTIPLEZONES10          |\
                                  FEATURE_READNONOEMDATA           |\
                                  FEATURE_DATALOADERS              |\
                                  FEATURE_DATALOADERS_EXCEPTONE    |\
                                  FEATURE_LOADONDEMAND             |\
                                  FEATURE_MULTIPLEWORKERTHREADS    |\
                                  FEATURE_ISOSURFACEGROUPS         |\
                                  FEATURE_SLICEGROUPS              |\
                                  FEATURE_STREAMTRACEGROUPS        |\
                                  FEATURE_FEPOLYHEDRON             |\
                                  FEATURE_FEPOLYGON)

#define VALID_FEATURE_INHIBIT_FLAG(feature) (((feature) & KnownFeaturesToInhibit) != 0)
#define VALID_FEATURE_INHIBIT_MASK(mask) (((mask) & ~KnownFeaturesToInhibit)==0)



/* The following are used by the OEM libs, so they need
   to be outside of TECPLOTKERNEL */
typedef    unsigned long  FeatureFlag_t;
typedef    unsigned long  FeatureMask_t;

/* ENDREMOVEFROMADDON */

typedef    char             SymbolChar_t[3];

/**
 * Face node offset used for identifying which node of a polytope face is
 * desired.
 */
typedef LgIndex_t FaceNodeOffset_t;

/**
 * Element face offset used for identifying which face of a polytope element is
 * desired.
 */
typedef LgIndex_t ElemFaceOffset_t;

/**
 * Face boundary item offset used for identifying which boundary item of a
 * polytope face is desired.
 */
typedef LgIndex_t FaceBndryItemOffset_t;

/****************************************************************
 *                                                              *
 *                     ENUMERATED TYPEDEFS                      *
 *                                                              *
 ****************************************************************/

typedef enum
{
    PlacementPlaneOrientation_X,
    PlacementPlaneOrientation_Y,
    PlacementPlaneOrientation_Z,
    END_PlacementPlaneOrientation_e,
    PlacementPlaneOrientation_Invalid = BadEnumValue
} PlacementPlaneOrientation_e;

typedef enum
{
    StringMode_ASCII,
    StringMode_UTF8,
    StringMode_Blend,
    END_StringMode_e,
    StringMode_Invalid = BadEnumValue

} StringMode_e;

typedef enum
{
    SidebarSizing_MaxOfAll,
    SidebarSizing_Dynamic,
    END_SidebarSizing_e,
    SidebarSizing_Invalid = BadEnumValue

} SidebarSizing_e;

typedef enum
{
    SidebarLocation_Left,
    SidebarLocation_Right,  /* Not allowed at this time */
    SidebarLocation_Top,    /* Not allowed at this time */
    SidebarLocation_Bottom, /* Not allowed at this time */
    END_SidebarLocation_e,
    SidebarLocation_Invalid = BadEnumValue

} SidebarLocation_e;

typedef enum
{
    MenuItem_Option,
    MenuItem_Toggle,
    MenuItem_Separator,
    MenuItem_SubMenu,
    END_MenuItem_e,
    MenuItem_Invalid = BadEnumValue
} MenuItem_e;

typedef enum
{
    StandardMenu_File,
    StandardMenu_Edit,
    StandardMenu_View,
    StandardMenu_Plot,
    StandardMenu_Insert,
    StandardMenu_Data,
    StandardMenu_Frame,
    StandardMenu_Workspace, /* deprecated: use Options instead */
    StandardMenu_Tools,
    StandardMenu_Help,
    StandardMenu_Animate,
    StandardMenu_Options,
    StandardMenu_Scripting,
    END_StandardMenu_e,
    StandardMenu_Invalid = BadEnumValue
} StandardMenu_e;

typedef enum
{
    FieldProbeDialogPage_NodalValues,
    FieldProbeDialogPage_CellCenteredValues,
    FieldProbeDialogPage_ZoneCellInfo,
    FieldProbeDialogPage_FaceNeighbors,
    END_FieldProbeDialogPage_e,
    FieldProbeDialogPage_Invalid = BadEnumValue
} FieldProbeDialogPage_e;

/* BEGINREMOVEFROMADDON */

/* used for caches of boolean type */
typedef enum
{
    BooleanCache_False, /* Value is cached and is FALSE */
    BooleanCache_True,  /* Value is cached and is TRUE */
    BooleanCache_Uncached, /* Value is not cached.  Value is unknown. */
    END_BooleanCache_e,
    BooleanCache_Invalid = BadEnumValue
} BooleanCache_e;

/*
 * For determining pick location along a line
 */
typedef enum
{
    LinePickLocation_None,
    LinePickLocation_StartHandle,
    LinePickLocation_MidLineOnHandle,
    LinePickLocation_MidLineOffHandles,
    LinePickLocation_EndHandle,
    END_LinePickLocation_e,
    LinePickLocation_Invalid = BadEnumValue
} LinePickLocation_e;


/*
 * Defines destination for setting up views: hardware (ie, OpenGL) or
 * software (ie, internal transformation matrices).
 */
typedef enum
{
    ViewDest_Hardware,
    ViewDest_Software,
    END_ViewDest_e,
    ViewDest_Invalid = BadEnumValue
} ViewDest_e;

/* used for identifying the origin of the dataset reader */
typedef enum
{
    DataSetReaderOrigin_Native,  /* created by Tecplot   */
    DataSetReaderOrigin_Foreign, /* created by an add-on */
    END_DataSetReaderOrigin_e,
    DataSetReaderOrigin_Invalid = BadEnumValue
} DataSetReaderOrigin_e;

/* used for identifying the origin of the extended curve fit */
typedef enum
{
    ExtendedCurveFitOrigin_Native,  /* created by Tecplot   */
    ExtendedCurveFitOrigin_Foreign, /* created by an add-on */
    END_ExtendedCurveFitOrigin_e,
    ExtendedCurveFitOrigin_Invalid = BadEnumValue
} ExtendedCurveFitOrigin_e;

typedef enum
{
    CollapsedStatus_NotCollapsed,
    CollapsedStatus_CollapsedToPoint,
    CollapsedStatus_CollapsedToLine,
    CollapsedStatus_CollapsedToSegmentedLine,
    CollapsedStatus_CollapsedToTriangle,
    END_CollapsedStatus_e,
    CollapsedStatus_Invalid = BadEnumValue
} CollapsedStatus_e;
/* ENDREMOVEFROMADDON */

/**
 */
typedef enum
{
    UndoStateCategory_Frame,
    UndoStateCategory_Picked,            /* picked changes, not the pick itself */
    UndoStateCategory_Text,
    UndoStateCategory_Geom,
    UndoStateCategory_View,
    UndoStateCategory_WorkspaceView,
    UndoStateCategory_Style,            /* style less text and geometries */
    UndoStateCategory_SpecificStyle,    /* meaning that specific undo style will be added by the caller */
    UndoStateCategory_Data,
    UndoStateCategory_DataAndStyle,
    UndoStateCategory_StyleIncTextGeom, /* style including text and geometires */
    UndoStateCategory_GlobalStyle,      /* style less field, map, text and geometries */
    UndoStateCategory_PageAction,
    END_UndoStateCategory_e,
    UndoStateCategory_Invalid = BadEnumValue
} UndoStateCategory_e;


/*
 * Used only for Action_PropagateLinking
 */
typedef enum
{
    LinkType_WithinFrame,
    LinkType_BetweenFrames,
    END_LinkType_e,
    LinkType_Invalid = BadEnumValue
} LinkType_e;

typedef enum
{
    FrameCollection_All,
    FrameCollection_Picked,
    END_FrameCollection_e,
    FrameCollection_Invalid = BadEnumValue
} FrameCollection_e;



typedef enum
{
    LegendProcess_DrawLegend,
    LegendProcess_EraseLegend,
    LegendProcess_GetExtents,
    END_LegendProcess_e,
    LegendProcess_Invalid = BadEnumValue
} LegendProcess_e;


typedef enum
{
    RGBLegendOrientation_RGB,
    RGBLegendOrientation_GBR,
    RGBLegendOrientation_BRG,
    RGBLegendOrientation_RBG,
    RGBLegendOrientation_GRB,
    RGBLegendOrientation_BGR,
    END_RGBLegendOrientation_e,
    RGBLegendOrientation_Invalid = BadEnumValue
} RGBLegendOrientation_e;



/* BEGINREMOVEFROMADDON */
/* Used by some of the image exporters/animators */
typedef struct
{
    Byte_t  R;
    Byte_t  G;
    Byte_t  B;
} RGBTriple_s;

typedef RGBTriple_s RGBPalette_t[256];

/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
/* The tag on the following line is so that the Windows
   build script can parse all of the current state changes
   out of this file, and compare them to the state changes
   found in the main.c template file.
   Do not change or delete the line below.*/
/*StateChange_e_BeginDef*/
/* ENDREMOVEFROMADDON */

typedef enum
{
    StateChange_VarsAltered,
    StateChange_VarsAdded,
    StateChange_ZonesDeleted,
    StateChange_ZonesAdded,
    StateChange_NodeMapsAltered,
    StateChange_FrameDeleted,
    StateChange_NewTopFrame,               /* deprecated: use NewActiveFrame and/or FrameOrderChange */
    StateChange_Style,
    StateChange_DataSetReset,
    StateChange_NewLayout,
    StateChange_CompleteReset,             /* deprecated: no longer broadcast */
    StateChange_LineMapAssignment,         /* was StateChange_XYMapAssignment */
    StateChange_ContourLevels,
    StateChange_ModalDialogLaunch,
    StateChange_ModalDialogDismiss,
    StateChange_QuitTecplot,
    StateChange_ZoneName,
    StateChange_VarName,
    StateChange_LineMapName,               /* was StateChange_XYMapName */
    StateChange_LineMapAddDeleteOrReorder, /* was StateChange_XYMapAddDeleteOrReorder */
    StateChange_View,
    StateChange_ColorMap,
    StateChange_ContourVar,
    StateChange_Streamtrace,
    StateChange_NewAxisVariables,
    StateChange_MouseModeUpdate,
    StateChange_PickListCleared,
    StateChange_PickListGroupSelect,
    StateChange_PickListSingleSelect,
    StateChange_PickListStyle,
    StateChange_DataSetFileName,
    StateChange_UnsuspendInterface,        /* was StateChange_DrawGraphicsOn */
    StateChange_SuspendInterface,          /* was StateChange_DrawGraphicsOff */
    StateChange_DataSetLockOn,
    StateChange_DataSetLockOff,
    StateChange_Text,
    StateChange_Geom,
    StateChange_DataSetTitle,
    StateChange_DrawingInterrupted,
    StateChange_PrintPreviewLaunch,
    StateChange_PrintPreviewDismiss,
    StateChange_AuxDataAdded,
    StateChange_AuxDataDeleted,
    StateChange_AuxDataAltered,
    StateChange_VarsDeleted,
    StateChange_TecplotIsInitialized,
    StateChange_ImageExported,
    StateChange_VariableLockOn,
    StateChange_VariableLockOff,
    StateChange_PageDeleted,
    StateChange_NewTopPage,
    StateChange_NewActiveFrame,
    StateChange_FrameOrderChanged,
    END_StateChange_e,
    StateChange_Invalid = BadEnumValue,
    /* deprecated values */
    StateChange_DrawGraphicsOn          = StateChange_UnsuspendInterface,
    StateChange_DrawGraphicsOff         = StateChange_SuspendInterface,
    StateChange_XYMapAssignment         = StateChange_LineMapAssignment,
    StateChange_XYMapName               = StateChange_LineMapName,
    StateChange_XYMapAddDeleteOrReorder = StateChange_LineMapAddDeleteOrReorder
} StateChange_e;

typedef enum
{
    StateChangeMode_v75,
    StateChangeMode_v80,
    StateChangeMode_v100,
    StateChangeMode_v113,
    END_StateChangeMode_e,
    StateChangeMode_Invalid = BadEnumValue
} StateChangeMode_e;

typedef enum
{
    StateChangeCallbackAPI_Classic,
    StateChangeCallbackAPI_ChangeOnly,
    StateChangeCallbackAPI_ChangePlusClient,
    END_StateChangeCallbackAPI_e,
    StateChangeCallbackAPI_Invalid = BadEnumValue
} StateChangeCallbackAPI_e;

typedef enum
{
    AppMode_Normal,
    AppMode_Demo,
    AppMode_OEM,
    END_AppMode_e,
    AppMode_Invalid = BadEnumValue
} AppMode_e;

typedef enum
{
    ProductFlavor_TecplotFocus,
    ProductFlavor_Tecplot360,
    ProductFlavor_TecplotRS,
    ProductFlavor_TecplotSDK,
    END_ProductFlavor_e,
    ProductFlavor_Invalid = BadEnumValue,
    ProductFlavor_Focus = ProductFlavor_TecplotFocus, /* deprecated */
    ProductFlavor_360   = ProductFlavor_Tecplot360,   /* deprecated */
    ProductFlavor_RS    = ProductFlavor_TecplotRS,    /* deprecated */
    ProductFlavor_SDK   = ProductFlavor_TecplotSDK    /* deprecated */
} ProductFlavor_e;

typedef enum
{
    LayoutPackageObject_Image,
    LayoutPackageObject_Layout,
    LayoutPackageObject_Data,
    END_LayoutPackageObject_e,
    LayoutPackageObject_Invalid = BadEnumValue
} LayoutPackageObject_e;

typedef enum
{
    VarLoadMode_ByName,
    VarLoadMode_ByPosition,
    END_VarLoadMode_e,
    VarLoadMode_Invalid = BadEnumValue
} VarLoadMode_e;

typedef enum
{
    ImageSelection_OnePerFrame,
    ImageSelection_WorkspaceOnly,
    END_ImageSelection_e,
    ImageSelection_Invalid = BadEnumValue
} ImageSelection_e;

typedef enum
{
    LibraryType_Foreign,
    LibraryType_V7Standard,
    LibraryType_V7ActiveX,
    END_LibraryType_e,
    LibraryType_Invalid = BadEnumValue
} LibraryType_e; /* <help> "Add-on types" */


typedef enum
{
    AssignOp_Equals,
    AssignOp_PlusEquals,
    AssignOp_MinusEquals,
    AssignOp_TimesEquals,
    AssignOp_DivideEquals,
    AssignOp_ConvertFromCm,
    AssignOp_ConvertFromIn,
    AssignOp_ConvertFromPt,
    AssignOp_ConvertFromPix,
    END_AssignOp_e,
    AssignOp_Invalid = BadEnumValue
} AssignOp_e;

typedef enum
{
    Dialog_ColorMap,
    Dialog_Equation,
    Dialog_MacroViewer,
    Dialog_ZoneMapStyle, /* was Dialog_PlotAttributes*/
    Dialog_QuickEdit,
    Dialog_QuickMacroPanel,
    Dialog_ValueBlanking,
    Dialog_Probe,          /* used for dialog positioning only */
    Dialog_ProbeAt,
    Dialog_NewLayout,
    Dialog_OpenLayout,
    Dialog_Save,
    Dialog_SaveAs,
    Dialog_LoadData,
    Dialog_WriteData,
    Dialog_Print,
    Dialog_Import,
    Dialog_Export,
    Dialog_MacroPlay,
    Dialog_MacroRecord,
    Dialog_AxisEdit,
    Dialog_SpatialVars,
    Dialog_Reset3DAxes,
    Dialog_ThreeDAxisLimits,
    Dialog_ThreeDOrientationAxis,
    Dialog_Streamtraces,
    Dialog_IsoSurfaces,
    Dialog_Slices,
    Dialog_Contour,
    Dialog_VectorLength,
    Dialog_VectorVars,
    Dialog_VectorArrowheads,
    Dialog_VectorReferenceVector,
    Dialog_ScatterSizeAndFont,
    Dialog_ScatterLegend,
    Dialog_ScatterReferenceSymbol,
    Dialog_RGBColorVarsAndRange,
    Dialog_RGBColorLegend,
    Dialog_LineMapLegend,
    Dialog_IJKBlanking,
    Dialog_DepthBlanking,
    Dialog_LightSource,
    Dialog_Advanced3DControl,
    Dialog_TwoDDrawOrder,
    Dialog_PolarDrawingOptions,
    Dialog_DataLabels,
    Dialog_StyleLinking,
    Dialog_Smooth,
    Dialog_TransformCoordinates,
    Dialog_Rotate2DData,
    Dialog_Create1DLine,
    Dialog_CreateRectangularZone,
    Dialog_CreateCircularZone,
    Dialog_DuplicateZone,
    Dialog_MirrorZone,
    Dialog_CreateZoneFromPolylines,
    Dialog_CreateZoneFromValues,
    Dialog_DeleteVariables,
    Dialog_DeleteZones,
    Dialog_ExtractContourLines,
    Dialog_ExtractFEBoundary,
    Dialog_ExtractIsoSurfaces,
    Dialog_ExtractSlices,
    Dialog_ExtractSliceFromPlane,
    Dialog_ExtractStreamtraces,
    Dialog_ExtractSubZone,
    Dialog_ExtractDiscretePoints,
    Dialog_ExtractPointsFromPolyline,
    Dialog_ExtractPointsFromGeometry,
    Dialog_LinearInterpolation,
    Dialog_InverseDistanceInterpolation,
    Dialog_KrigingInterpolation,
    Dialog_Triangulate,
    Dialog_DataInfo,
    Dialog_CurveInfo,
    Dialog_DataSpreadsheet,
    Dialog_PaperSetup,
    Dialog_OrderFrames,
    Dialog_RulerGrid,
    Dialog_ThreeDViewRotate,
    Dialog_ThreeDViewDetails,
    Dialog_TranslateMagnify,
    Dialog_PrintPreview,
    Dialog_ColorPreferences,
    Dialog_MiscPreferences,
    Dialog_SizePreferences,
    Dialog_SaveConfiguration,
    Dialog_SaveColorMap,
    Dialog_LoadColorMap,
    Dialog_HelpAboutTecplot,
    Dialog_HelpAboutAddOns,
    Dialog_Publish,
    Dialog_EditFrame,
    Dialog_CopyToClipboard,
    Dialog_ThreeDEdge,
    Dialog_TimeDetails,
    Dialog_Performance,
    Dialog_HelpTecplotLicensing,
    Dialog_GeomDetails,
    Dialog_BasicColorLegend,
    Dialog_FourierTransform,
    END_Dialog_e,
    Dialog_Invalid = BadEnumValue,
    /* deprecated values */
    Dialog_PlotAttributes = Dialog_ZoneMapStyle
} Dialog_e; /* <help> "Tecplot dialog types" */

typedef enum
{
    AnchorAlignment_TopLeft,
    AnchorAlignment_TopCenter,
    AnchorAlignment_TopRight,
    AnchorAlignment_MiddleLeft,
    AnchorAlignment_MiddleCenter,
    AnchorAlignment_MiddleRight,
    AnchorAlignment_BottomLeft,
    AnchorAlignment_BottomCenter,
    AnchorAlignment_BottomRight,
    END_AnchorAlignment_e,
    AnchorAlignment_Invalid = BadEnumValue
} AnchorAlignment_e;

/* BEGINREMOVEFROMADDON */
typedef enum
{
    PositionAtAnchor_Never,
    PositionAtAnchor_Once,
    PositionAtAnchor_Always,
    END_PositionAtAnchor_e,
    PositionAtAnchor_Invalid = BadEnumValue
} PositionAtAnchor_e;
/* ENDREMOVEFROMADDON */

/* BEGINREMOVEFROMADDON */
typedef struct
{
    AnchorAlignment_e  AnchorAlignment;
    Boolean_t          AnchorHorizontalInside;
    Boolean_t          AnchorVerticalInside;
    SmInteger_t        MinVisibilityPercentage;
    LgIndex_t          IOffset;
    LgIndex_t          JOffset;
    PositionAtAnchor_e PositionAtAnchor;
    Boolean_t          HasBeenPositioned; /* not persistent */
} DialogPosition_s;
/* ENDREMOVEFROMADDON */


#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref CurveInfoMode_e instead.
 */
typedef enum
{
    ProcessXYMode_NotUsed1,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed2,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed3,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed4,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed5,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed6,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed7,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed8,         /* deprecated: do not use                     */
    ProcessXYMode_NotUsed9,         /* deprecated: do not use                     */
    ProcessXYMode_WriteCurveCoef,   /* deprecated: use CurveInfoMode_Coefficients */
    ProcessXYMode_WriteCurvePoints, /* deprecated: use CurveInfoMode_RawData      */
    END_ProcessXYMode_e,
    ProcessXYMode_Invalid = BadEnumValue
} ProcessXYMode_e;
#endif

typedef enum
{
    CurveInfoMode_Coefficients, /* ProcessXYMode_WriteCurveCoef   */
    CurveInfoMode_RawData,      /* ProcessXYMode_WriteCurvePoints */
    CurveInfoMode_Macro,        /* ProcessXYMode_WriteCurveCoefMacro */
    END_CurveInfoMode_e,
    CurveInfoMode_Invalid = BadEnumValue
} CurveInfoMode_e;

/* BEGINREMOVEFROMADDON */
typedef enum
{
    ProcessLineMapMode_Draw,
    ProcessLineMapMode_GetXYMinMax,
    ProcessLineMapMode_GetDataMinMax,
    ProcessLineMapMode_GetSinglePick,
    ProcessLineMapMode_CheckOnlyForGroupPick,
    ProcessLineMapMode_GetGroupPick,
    ProcessLineMapMode_GetFirstValidDataPoint,
    ProcessLineMapMode_GetNearestPoint,
    ProcessLineMapMode_GetDependentValue,
    ProcessLineMapMode_GetRSquaredGoodness,
    ProcessLineMapMode_DisplayCurveCoef,
    ProcessLineMapMode_WriteCurveCoef,
    ProcessLineMapMode_WriteCurvePoints,
    ProcessLineMapMode_InsertLabels,
    ProcessLineMapMode_GetIndependentValue,
    ProcessLineMapMode_WriteCurveCoefMacro,
    END_ProcessLineMapMode_e,
    ProcessLineMapMode_Invalid = BadEnumValue
} ProcessLineMapMode_e;
/* ENDREMOVEFROMADDON */

typedef enum
{
    StyleBase_Factory,
    StyleBase_Config,
    END_StyleBase_e,
    StyleBase_Invalid = BadEnumValue
} StyleBase_e;


typedef enum
{
    ReadDataOption_NewData,
    ReadDataOption_AppendData,
    ReadDataOption_ReplaceData,
    END_ReadDataOption_e,
    ReadDataOption_Invalid = BadEnumValue
} ReadDataOption_e;

#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref LabelType_e instead.
 */
typedef enum
{
    NodeLabel_Index,         /* deprecated: use LabelType_Index         */
    NodeLabel_VarValue,      /* deprecated: use LabelType_VarValue      */
    NodeLabel_XAndYVarValue, /* deprecated: use LabelType_XAndYVarValue */
    END_NodeLabel_e,
    NodeLabel_Invalid = BadEnumValue
} NodeLabel_e;
#endif

typedef enum
{
    LabelType_Index,         /* NodeLabel_Index         */
    LabelType_VarValue,      /* NodeLabel_VarValue      */
    LabelType_XAndYVarValue, /* NodeLabel_XAndYVarValue */
    END_LabelType_e,
    LabelType_Invalid = BadEnumValue
} LabelType_e;


#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref BorderAction_e instead.
 */
typedef enum
{
    SubBoundaryEditOption_All,      /* deprecated: use BorderAction_AddAll  */
    SubBoundaryEditOption_Add,      /* deprecated: use BorderAction_Add     */
    SubBoundaryEditOption_Remove,   /* deprecated: use BorderAction_Remove  */
    SubBoundaryEditOption_AddOnly,  /* deprecated: use BorderAction_AddOnly */
    END_SubBoundaryEditOption_e,
    SubBoundaryEditOption_Invalid = BadEnumValue
} SubBoundaryEditOption_e;
#endif

typedef enum
{
    BorderAction_AddAll,   /* SubBoundaryEditOption_All     */
    BorderAction_Add,      /* SubBoundaryEditOption_Add     */
    BorderAction_Remove,   /* SubBoundaryEditOption_Remove  */
    BorderAction_AddOnly,  /* SubBoundaryEditOption_AddOnly */
    END_BorderAction_e,
    BorderAction_Invalid = BadEnumValue
} BorderAction_e;


typedef enum
{
    PointerStyle_NotUsed1,
    PointerStyle_NotUsed2,
    PointerStyle_NotUsed3,
    PointerStyle_AllDirections,
    PointerStyle_NotUsed4,
    PointerStyle_NotUsed5,
    PointerStyle_NotUsed6,
    PointerStyle_UpperLeftBracket,
    PointerStyle_UpperRightBracket,
    PointerStyle_LeftBracket,
    PointerStyle_LowerLeftBracket,
    PointerStyle_LowerRightBracket,
    PointerStyle_RightBracket,
    PointerStyle_BottomBracket,
    PointerStyle_TopBracket,
    PointerStyle_UpDown,
    PointerStyle_LeftRight,
    END_PointerStyle_e,
    PointerStyle_Invalid = BadEnumValue
} PointerStyle_e;

typedef enum
{
    CursorStyle_Undefined,
    CursorStyle_StandardArrow,
    CursorStyle_AdjusterArrow,
    CursorStyle_AllDirections,
    CursorStyle_Rotate,
    CursorStyle_Zoom,
    CursorStyle_Locate,
    CursorStyle_UpperLeftBracket,
    CursorStyle_UpperRightBracket,
    CursorStyle_LeftBracket,
    CursorStyle_LowerLeftBracket,
    CursorStyle_LowerRightBracket,
    CursorStyle_RightBracket,
    CursorStyle_BottomBracket,
    CursorStyle_TopBracket,
    CursorStyle_UpDown,
    CursorStyle_LeftRight,
    CursorStyle_Waiting,
    END_CursorStyle_e,
    CursorStyle_Invalid = BadEnumValue
} CursorStyle_e;


typedef enum
{
    PickSubPosition_All,
    PickSubPosition_Top,
    PickSubPosition_Bottom,
    PickSubPosition_Left,
    PickSubPosition_Right,
    PickSubPosition_TopLeft,
    PickSubPosition_TopRight,
    PickSubPosition_BottomLeft,
    PickSubPosition_BottomRight,
    PickSubPosition_BottomAndTop,
    PickSubPosition_LeftAndRight,
    END_PickSubPosition_e,
    PickSubPosition_Invalid = BadEnumValue
} PickSubPosition_e;

typedef enum
{
    TecEngInitReturnCode_Ok,
    TecEngInitReturnCode_LicenseIsInvalid,
    TecEngInitReturnCode_LicenseExpired,
    TecEngInitReturnCode_InternalInitializationError,
    END_TecEngInitReturnCode_e,
    TecEngInitReturnCode_Invalid = BadEnumValue
} TecEngInitReturnCode_e;

typedef enum
{
    GetValueReturnCode_Ok,
    GetValueReturnCode_ResultTypeError,
    GetValueReturnCode_SyntaxError,
    GetValueReturnCode_ContextError,
    GetValueReturnCode_DeprecatedError,
    END_GetValueReturnCode_e,
    GetValueReturnCode_Invalid = BadEnumValue,
    /* deprecated values */
    GetValue_Ok              = GetValueReturnCode_Ok,              /* deprecated */
    GetValue_ResultTypeError = GetValueReturnCode_ResultTypeError, /* deprecated */
    GetValue_SyntaxError     = GetValueReturnCode_SyntaxError,     /* deprecated */
    GetValue_ContextError    = GetValueReturnCode_ContextError,    /* deprecated */
    GetValue_DeprecatedError = GetValueReturnCode_DeprecatedError, /* deprecated */
    GetValue_Invalid         = GetValueReturnCode_Invalid          /* deprecated */
} GetValueReturnCode_e;

typedef enum
{
    SetValueReturnCode_Ok,
    SetValueReturnCode_DuplicateValue,
    SetValueReturnCode_InvalidCommandOption,
    SetValueReturnCode_NoAttachedDatasetError,
    SetValueReturnCode_NoAttachedFrameError,
    SetValueReturnCode_NotAllowedInConfigError,
    SetValueReturnCode_ValueRangeError,
    SetValueReturnCode_ValueSyntaxError,
    SetValueReturnCode_AssignOpError,
    SetValueReturnCode_InvalidVarOrZone,
    SetValueReturnCode_InternalMemoryError,
    SetValueReturnCode_ContextError1,
    SetValueReturnCode_ContextError2,
    SetValueReturnCode_OnlyAllowedInConfigError,
    SetValueReturnCode_FeatureNotAvailable,
    END_SetValueReturnCode_e,
    /* BEGINREMOVEFROMADDON */
    /* For now this value is only used in Tecplot code.
     * the value is here as an option for the future. */
    SetValueReturnCode_Ignored = SetValueReturnCode_DuplicateValue,
    /* ENDREMOVEFROMADDON */
    SetValueReturnCode_Invalid = BadEnumValue,
    /* deprecated values */
    SetValue_Ok                       = SetValueReturnCode_Ok,                       /* deprecated */
    SetValue_DuplicateValue           = SetValueReturnCode_DuplicateValue,           /* deprecated */
    SetValue_InvalidCommandOption     = SetValueReturnCode_InvalidCommandOption,     /* deprecated */
    SetValue_NoAttachedDatasetError   = SetValueReturnCode_NoAttachedDatasetError,   /* deprecated */
    SetValue_NoAttachedFrameError     = SetValueReturnCode_NoAttachedFrameError,     /* deprecated */
    SetValue_NotAllowedInConfigError  = SetValueReturnCode_NotAllowedInConfigError,  /* deprecated */
    SetValue_ValueRangeError          = SetValueReturnCode_ValueRangeError,          /* deprecated */
    SetValue_ValueSyntaxError         = SetValueReturnCode_ValueSyntaxError,         /* deprecated */
    SetValue_AssignOpError            = SetValueReturnCode_AssignOpError,            /* deprecated */
    SetValue_InvalidVarOrZone         = SetValueReturnCode_InvalidVarOrZone,         /* deprecated */
    SetValue_InternalMemoryError      = SetValueReturnCode_InternalMemoryError,      /* deprecated */
    SetValue_ContextError1            = SetValueReturnCode_ContextError1,            /* deprecated */
    SetValue_ContextError2            = SetValueReturnCode_ContextError2,            /* deprecated */
    SetValue_OnlyAllowedInConfigError = SetValueReturnCode_OnlyAllowedInConfigError, /* deprecated */
    SetValue_FeatureNotAvailable      = SetValueReturnCode_FeatureNotAvailable,      /* deprecated */
    /* BEGINREMOVEFROMADDON */
    SetValue_Ignored                  = SetValueReturnCode_Ignored,                  /* deprecated */
    /* ENDREMOVEFROMADDON */
    SetValue_Invalid                  = SetValueReturnCode_Invalid                   /* deprecated */
} SetValueReturnCode_e;


typedef enum
{
    ObjectAlign_LeftJustify,
    ObjectAlign_RightJustify,
    ObjectAlign_Center,
    ObjectAlign_Top,
    ObjectAlign_Bottom,
    END_ObjectAlign_e,
    ObjectAlign_Invalid = BadEnumValue
} ObjectAlign_e;


/*
 * For 3D axis labels only.
 */
typedef enum
{
    LabelAlignment_ByAngle,
    LabelAlignment_AlongAxis,
    LabelAlignment_PerpendicularToAxis,
    END_LabelAlignment_e,
    LabelAlignment_Invalid = BadEnumValue
} LabelAlignment_e; /* <help> Label alignment for 3D axis labels only" */

/*
 * View_SetMagnification added 02/24/03 so all plot types
 * can behave the same way "do a 'centered' magnification change".
 * Line plots will still accept View_Scale option and zoom towards
 * the corner so old macros/addons still work.
 */
typedef enum
{
    View_Fit,
    View_DataFit,
    View_AxisFit,
    View_Scale,   /* deprecated, Use SetMagnification */
    View_Center,
    View_Translate,
    View_Zoom,
    View_Last,
    View_Copy,
    View_Paste,
    View_Push,  /* End of V9 enums */
    View_SetMagnification,
    View_NiceFit,
    View_AxisNiceFit,
    View_MakeCurrentViewNice,
    View_AxisMakeCurrentValuesNice,
    View_AxisResetToEntireCircle,
    View_FitSurfaces,
    END_View_e,
    View_Invalid = BadEnumValue
} View_e;



typedef enum
{
    WorkspaceView_FitSelectedFrames,
    WorkspaceView_FitAllFrames,
    WorkspaceView_FitPaper,
    WorkspaceView_Maximize,
    WorkspaceView_LastView,
    WorkspaceView_Zoom,
    WorkspaceView_Translate,
    WorkspaceView_UnMaximize,
    END_WorkspaceView_e,
    WorkspaceView_Invalid = BadEnumValue
} WorkspaceView_e;


typedef enum
{
    ArrowheadStyle_Plain,
    ArrowheadStyle_Filled,
    ArrowheadStyle_Hollow,
    END_ArrowheadStyle_e,
    ArrowheadStyle_Invalid = BadEnumValue,
    /* deprecated values */
    Arrowhead_Plain   = ArrowheadStyle_Plain,  /* deprecated */
    Arrowhead_Filled  = ArrowheadStyle_Filled, /* deprecated */
    Arrowhead_Hollow  = ArrowheadStyle_Hollow, /* deprecated */
    Arrowhead_Invalid = ArrowheadStyle_Invalid /* deprecated */
} ArrowheadStyle_e;


typedef enum
{
    ArrowheadAttachment_None,
    ArrowheadAttachment_AtBeginning,
    ArrowheadAttachment_AtEnd,
    ArrowheadAttachment_AtBothEnds,
    END_ArrowheadAttachment_e,
    ArrowheadAttachment_Invalid = BadEnumValue,
    /* deprecated values */
    ArrowheadAttach_None        = ArrowheadAttachment_None,        /* deprecated */
    ArrowheadAttach_AtBeginning = ArrowheadAttachment_AtBeginning, /* deprecated */
    ArrowheadAttach_AtEnd       = ArrowheadAttachment_AtEnd,       /* deprecated */
    ArrowheadAttach_AtBothEnds  = ArrowheadAttachment_AtBothEnds,  /* deprecated */
    ArrowheadAttach_Invalid     = ArrowheadAttachment_Invalid      /* deprecated */
} ArrowheadAttachment_e;

typedef enum
{
    Clipping_ClipToViewport,
    Clipping_ClipToFrame,
    END_Clipping_e,
    Clipping_Invalid = BadEnumValue
} Clipping_e;

typedef enum
{
    StatusInfo_Hover,
    StatusInfo_Identify,
    StatusInfo_Instruction,
    StatusInfo_Working,
    StatusInfo_PercentDone,
    END_StatusInfo_e,
    StatusInfo_Invalid = BadEnumValue
} StatusInfo_e;


#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref PlotType_e instead.
 */
typedef enum
{
    FrameMode_Empty,    /* deprecated: use PlotType_Automatic   */
    FrameMode_ThreeD,   /* deprecated: use PlotType_Cartesian3D */
    FrameMode_TwoD,     /* deprecated: use PlotType_Cartesian2D */
    FrameMode_XY,       /* deprecated: use PlotType_XYLine      */
    FrameMode_Sketch,   /* deprecated: use PlotType_Sketch      */
    END_FrameMode_e,
    FrameMode_Invalid = BadEnumValue,
    /* deprecated values */
    Frame_Empty   = FrameMode_Empty,    /* deprecated */
    Frame_ThreeD  = FrameMode_ThreeD,   /* deprecated */
    Frame_TwoD    = FrameMode_TwoD,     /* deprecated */
    Frame_XY      = FrameMode_XY,       /* deprecated */
    Frame_Sketch  = FrameMode_Sketch,   /* deprecated */
    Frame_Invalid = FrameMode_Invalid   /* deprecated */
} FrameMode_e;
#endif


typedef enum
{
    PlotType_Automatic,   /* Frame_Empty  */
    PlotType_Cartesian3D, /* Frame_ThreeD */
    PlotType_Cartesian2D, /* Frame_TwoD   */
    PlotType_XYLine,      /* Frame_XY     */
    PlotType_Sketch,      /* Frame_Sketch */
    PlotType_PolarLine,
    END_PlotType_e,
    PlotType_Invalid = BadEnumValue
} PlotType_e;


#define VALID_PLOTTYPE(PlotType) ( VALID_ENUM((PlotType), PlotType_e) && \
                                   ((PlotType) != PlotType_Automatic) )
#define VALID_LINEPLOT_PLOTTYPE(PlotType) ( (PlotType) == PlotType_XYLine || \
                                            (PlotType) == PlotType_PolarLine )
#define VALID_FIELDPLOT_PLOTTYPE(PlotType) ( (PlotType) == PlotType_Cartesian2D || \
                                             (PlotType) == PlotType_Cartesian3D )
#define PLOTTYPE_USES_FIELDZONES(PlotType) VALID_FIELDPLOT_PLOTTYPE((PlotType))
#define PLOTTYPE_USES_LINEMAPS(PlotType) VALID_LINEPLOT_PLOTTYPE((PlotType))
#define VALID_V9_PLOTTYPE(PlotType) ( (PlotType) == PlotType_Sketch || \
                                      (PlotType) == PlotType_XYLine || \
                                      (PlotType) == PlotType_Cartesian2D || \
                                      (PlotType) == PlotType_Cartesian3D )


typedef enum
{
    ContLineCreateMode_OneZonePerContourLevel,
    ContLineCreateMode_OneZonePerIndependentPolyline,
    END_ContLineCreateMode_e,
    ContLineCreateMode_Invalid = BadEnumValue
} ContLineCreateMode_e;


typedef enum
{
    PickObjects_None,
    PickObjects_Frame,
    PickObjects_Axis,
    PickObjects_ThreeDOrientationAxis,
    PickObjects_Geom,
    PickObjects_Text,
    PickObjects_ContourLegend,
    PickObjects_ContourLabel,
    PickObjects_ScatterLegend,
    PickObjects_LineLegend,
    PickObjects_ReferenceVector,
    PickObjects_ReferenceScatterSymbol,
    PickObjects_StreamtracePosition,
    PickObjects_StreamtraceTermLine,
    PickObjects_Paper,
    PickObjects_Zone,
    PickObjects_XYMapping, /* deprecated: use PickObject_LineMapping */
    PickObjects_StreamtraceCOB,
    PickObjects_SliceCOB,
    PickObjects_IsoSurfaceCOB,
    PickObjects_RGBLegend,
    PickObjects_LineMapping,
    PickObjects_BasicColorLegend,
    END_PickObjects_e,
    PickObjects_Invalid = BadEnumValue,
    /* deprecated values */
    PickObject_None                   = PickObjects_None,                   /* deprecated */
    PickObject_Frame                  = PickObjects_Frame,                  /* deprecated */
    PickObject_Axis                   = PickObjects_Axis,                   /* deprecated */
    PickObject_3DOrientationAxis      = PickObjects_ThreeDOrientationAxis,  /* deprecated */
    PickObject_Geom                   = PickObjects_Geom,                   /* deprecated */
    PickObject_Text                   = PickObjects_Text,                   /* deprecated */
    PickObject_ContourLegend          = PickObjects_ContourLegend,          /* deprecated */
    PickObject_ContourLabel           = PickObjects_ContourLabel,           /* deprecated */
    PickObject_ScatterLegend          = PickObjects_ScatterLegend,          /* deprecated */
    PickObject_LineLegend             = PickObjects_LineLegend,             /* deprecated */
    PickObject_XYLegend               = PickObjects_LineLegend,             /* deprecated */
    PickObject_ReferenceVector        = PickObjects_ReferenceVector,        /* deprecated */
    PickObject_ReferenceScatterSymbol = PickObjects_ReferenceScatterSymbol, /* deprecated */
    PickObject_StreamtracePosition    = PickObjects_StreamtracePosition,    /* deprecated */
    PickObject_StreamtraceTermLine    = PickObjects_StreamtraceTermLine,    /* deprecated */
    PickObject_Paper                  = PickObjects_Paper,                  /* deprecated */
    PickObject_Zone                   = PickObjects_Zone,                   /* deprecated */
    PickObject_XYMapping              = PickObjects_XYMapping,              /* deprecated */
    PickObject_StreamtraceCOB         = PickObjects_StreamtraceCOB,         /* deprecated */
    PickObject_SliceCOB               = PickObjects_SliceCOB,               /* deprecated */
    PickObject_IsoSurfaceCOB          = PickObjects_IsoSurfaceCOB,          /* deprecated */
    PickObject_RGBLegend              = PickObjects_RGBLegend,              /* deprecated */
    PickObject_LineMapping            = PickObjects_LineMapping,            /* deprecated */
    PickObject_Invalid                = PickObjects_Invalid                 /* deprecated */
} PickObjects_e;


/* BEGINREMOVEFROMADDON */
typedef enum
{
    SingleEditState_NotEditing,
    SingleEditState_ActivelyEditing,
    SingleEditState_WasEditing,
    END_SingleEditState_e,
    EditingInvalid = BadEnumValue
} SingleEditState_e;


typedef enum
{
    AxisSubObject_GridArea,
    AxisSubObject_AxisLine,
    AxisSubObject_Title,
    END_AxisSubObject_e,
    AxisSubObject_Invalid = BadEnumValue
} AxisSubObject_e;

typedef enum
{
    AxisSubPosition_GridMinBorder,
    AxisSubPosition_GridMaxBorder,
    AxisSubPosition_MainAxisLine,
    AxisSubPosition_BackAxisLine,
    AxisSubPosition_PerpAxisLine,
    AxisSubPosition_PerpBackAxisLine,
    END_AxisSubPosition_e,
    AxisSubPosition_Invalid = BadEnumValue,
    AxisSubPosition_2DStart = AxisSubPosition_GridMinBorder,
    AxisSubPosition_2DEnd = AxisSubPosition_MainAxisLine,
    AxisSubPosition_PolarStart = AxisSubPosition_GridMinBorder,
    AxisSubPosition_PolarEnd = AxisSubPosition_PerpBackAxisLine
} AxisSubPosition_e;
/* ENDREMOVEFROMADDON */

/*
 * NOTE: The _NoOp value is not at the top so this
 *       enumeration aligns with the old AltMouseButtonMode_e
 *       enumeration.
 */
typedef enum
{
    MouseButtonClick_Redraw,
    MouseButtonClick_RevertToSelect,
    MouseButtonClick_NoOp,
    END_MouseButtonClick_e,
    MouseButtonClick_Invalid = BadEnumValue
} MouseButtonClick_e;


typedef enum
{
    MouseButtonDrag_NoOp,
    MouseButtonDrag_ZoomPaper,
    MouseButtonDrag_TranslatePaper,
    MouseButtonDrag_ZoomData,
    MouseButtonDrag_TranslateData,
    MouseButtonDrag_RlrBallRtatData,
    MouseButtonDrag_SpherZRtatData,     /* Was SpherRtatData*/
    MouseButtonDrag_XRotateData,
    MouseButtonDrag_YRotateData,
    MouseButtonDrag_ZRotateData,
    MouseButtonDrag_TwistRotateData,
    MouseButtonDrag_ZoomViewer,
    MouseButtonDrag_TranslateViewer,
    MouseButtonDrag_RlrBallRtatVwr,
    MouseButtonDrag_SpherZRotateVwr,    /* Was SpherRotateVwr*/
    MouseButtonDrag_XRotateViewer,
    MouseButtonDrag_YRotateViewer,
    MouseButtonDrag_ZRotateViewer,
    MouseButtonDrag_TwistRotateViewer,
    MouseButtonDrag_SpherXRtatData,
    MouseButtonDrag_SpherYRtatData,
    MouseButtonDrag_SpherXRotateVwr,
    MouseButtonDrag_SpherYRotateVwr,
    END_MouseButtonDrag_e,
    MouseButtonDrag_Invalid = BadEnumValue,
    /* deprecated values */
    MouseButtonDrag_SpherRtatData  = MouseButtonDrag_SpherZRtatData,
    MouseButtonDrag_SpherRotateVwr = MouseButtonDrag_SpherZRotateVwr
} MouseButtonDrag_e;


/* BEGINREMOVEFROMADDON */
typedef struct
{
    MouseButtonClick_e ButtonClick;
    MouseButtonDrag_e  SimpleDrag;
    MouseButtonDrag_e  ControlledDrag;
    MouseButtonDrag_e  AltedDrag;
    MouseButtonDrag_e  ShiftedDrag;
    MouseButtonDrag_e  ControlAltedDrag;
    MouseButtonDrag_e  ControlShiftedDrag;
    MouseButtonDrag_e  AltShiftedDrag;
    MouseButtonDrag_e  ControlAltShiftedDrag;
} MouseButtonAction_s;


typedef struct
{
    MouseButtonAction_s MiddleButton;
    MouseButtonAction_s RightButton;
} MouseActions_s;
/* ENDREMOVEFROMADDON */


typedef enum  /* deprecated */
{
    AltMouseButtonMode_Regen,
    AltMouseButtonMode_RevertToSelect,
    END_AltMouseButtonMode_e,
    AltMouseButtonMode_Invalid = BadEnumValue
} AltMouseButtonMode_e;


typedef enum
{
    MouseButtonMode_NoMode,
    MouseButtonMode_Select,
    MouseButtonMode_Adjust,
    MouseButtonMode_Zoom,
    MouseButtonMode_Translate,
    MouseButtonMode_Probe,
    MouseButtonMode_Text,
    MouseButtonMode_GeomPolyline,
    MouseButtonMode_GeomSquare,
    MouseButtonMode_GeomCircle,
    MouseButtonMode_GeomRectangle,
    MouseButtonMode_GeomEllipse,
    MouseButtonMode_GeomSpline,
    MouseButtonMode_CreateFrame,
    MouseButtonMode_RotateSphericalZ,   /* Was MouseButtonMode_RotateSpherical */
    MouseButtonMode_RotateRollerBall,
    MouseButtonMode_RotateTwist,
    MouseButtonMode_RotateXAxis,
    MouseButtonMode_RotateYAxis,
    MouseButtonMode_RotateZAxis,
    MouseButtonMode_ContourLabel,
    MouseButtonMode_ContourAdd,
    MouseButtonMode_ContourDelete,
    MouseButtonMode_StreamPoints,
    MouseButtonMode_StreamEndLine,
    MouseButtonMode_ExtractPoints,
    MouseButtonMode_ExtractLine,
    MouseButtonMode_CreateRectangularZone,
    MouseButtonMode_CreateCircularZone,
    MouseButtonMode_Slice,
    MouseButtonMode_LightSource,
    MouseButtonMode_User1,
    MouseButtonMode_User2,
    MouseButtonMode_User3,
    MouseButtonMode_User4,
    MouseButtonMode_RotateSphericalX,
    MouseButtonMode_RotateSphericalY,
    MouseButtonMode_AdvancedAdjust,
    END_MouseButtonMode_e,
    MouseButtonMode_Invalid = BadEnumValue,
    /* deprecated values */
    MouseButtonMode_RotateSpherical = MouseButtonMode_RotateSphericalZ,
    Mouse_NoMode                = MouseButtonMode_NoMode,                /* deprecated */
    Mouse_Select                = MouseButtonMode_Select,                /* deprecated */
    Mouse_Adjust                = MouseButtonMode_Adjust,                /* deprecated */
    Mouse_Zoom                  = MouseButtonMode_Zoom,                  /* deprecated */
    Mouse_Translate             = MouseButtonMode_Translate,             /* deprecated */
    Mouse_Probe                 = MouseButtonMode_Probe,                 /* deprecated */
    Mouse_Text                  = MouseButtonMode_Text,                  /* deprecated */
    Mouse_GeomPolyline          = MouseButtonMode_GeomPolyline,          /* deprecated */
    Mouse_GeomSquare            = MouseButtonMode_GeomSquare,            /* deprecated */
    Mouse_GeomCircle            = MouseButtonMode_GeomCircle,            /* deprecated */
    Mouse_GeomRectangle         = MouseButtonMode_GeomRectangle,         /* deprecated */
    Mouse_GeomEllipse           = MouseButtonMode_GeomEllipse,           /* deprecated */
    Mouse_GeomSpline            = MouseButtonMode_GeomSpline,            /* deprecated */
    Mouse_CreateFrame           = MouseButtonMode_CreateFrame,           /* deprecated */
    Mouse_RotateSpherical       = MouseButtonMode_RotateSphericalZ,      /* deprecated */
    Mouse_RotateRollerBall      = MouseButtonMode_RotateRollerBall,      /* deprecated */
    Mouse_RotateTwist           = MouseButtonMode_RotateTwist,           /* deprecated */
    Mouse_RotateXAxis           = MouseButtonMode_RotateXAxis,           /* deprecated */
    Mouse_RotateYAxis           = MouseButtonMode_RotateYAxis,           /* deprecated */
    Mouse_RotateZAxis           = MouseButtonMode_RotateZAxis,           /* deprecated */
    Mouse_ContourLabel          = MouseButtonMode_ContourLabel,          /* deprecated */
    Mouse_ContourAdd            = MouseButtonMode_ContourAdd,            /* deprecated */
    Mouse_ContourDelete         = MouseButtonMode_ContourDelete,         /* deprecated */
    Mouse_StreamPoints          = MouseButtonMode_StreamPoints,          /* deprecated */
    Mouse_StreamEndLine         = MouseButtonMode_StreamEndLine,         /* deprecated */
    Mouse_ExtractPoints         = MouseButtonMode_ExtractPoints,         /* deprecated */
    Mouse_ExtractLine           = MouseButtonMode_ExtractLine,           /* deprecated */
    Mouse_CreateRectangularZone = MouseButtonMode_CreateRectangularZone, /* deprecated */
    Mouse_CreateCircularZone    = MouseButtonMode_CreateCircularZone,    /* deprecated */
    Mouse_Slice                 = MouseButtonMode_Slice,                 /* deprecated */
    Mouse_User1                 = MouseButtonMode_User1,                 /* deprecated */
    Mouse_User2                 = MouseButtonMode_User2,                 /* deprecated */
    Mouse_User3                 = MouseButtonMode_User3,                 /* deprecated */
    Mouse_User4                 = MouseButtonMode_User4,                 /* deprecated */
    Mouse_Invalid               = MouseButtonMode_Invalid                /* deprecated */
} MouseButtonMode_e;


typedef enum
{
    DetailsButtonState_QuickEdit,
    DetailsButtonState_ObjectDetails,
    DetailsButtonState_ToolDetails,
    END_DetailsButtonState_e,
    DetailsButtonState_Invalid = BadEnumValue
} DetailsButtonState_e;


typedef enum
{
    Event_ButtonPress,
    Event_ButtonRelease,
    Event_ButtonDoublePress,
    Event_Motion,
    Event_Drag,
    Event_KeyPress,
    END_Event_e,
    Event_Invalid = BadEnumValue
} Event_e;


typedef enum
{
    ObjectDrawMode_DrawFirst,
    ObjectDrawMode_Move,
    ObjectDrawMode_Remove,
    ObjectDrawMode_Place,
    END_ObjectDrawMode_e,
    ObjectDrawMode_Invalid = BadEnumValue
} ObjectDrawMode_e;


typedef enum
{
    ThreeDViewChangeDrawLevel_Full,
    ThreeDViewChangeDrawLevel_Trace,
    END_ThreeDViewChangeDrawLevel_e,
    ThreeDViewChangeDrawLevel_Invalid = BadEnumValue
} ThreeDViewChangeDrawLevel_e; /* <help> "ThreeDViewChangeDrawLevel is deprecated. Use PlotApproximateMode.\n"*/

typedef enum
{
    NonCurrentFrameRedrawLevel_Full,
    NonCurrentFrameRedrawLevel_Trace,
    END_NonCurrentFrameRedrawLevel_e,
    NonCurrentFrameRedrawLevel_Invalid = BadEnumValue
} NonCurrentFrameRedrawLevel_e; /* <help> "NonCurrentFrameRedrawLevel is deprecated. Use PlotApproximateMode.\n"*/


/**
 * Enumerates the redraw reasons and is passed as an argument to registered
 * draw event callbacks.
 *
 *   - RedrawReason_UserReqRedrawActiveFrame:\n
 *     The full draw event is in response to the "redraw" action function.
 *
 *   - RedrawReason_UserReqTraceActiveFrame:\n
 *     The approximate draw event is in response to the "redraw" action function.
 *
 *   - RedrawReason_UserReqRedrawAllFrames:\n
 *     The full draw event is in response to the "redraw all" action function.
 *
 *   - RedrawReason_UserReqTraceAllFrames:\n
 *     The approximate draw event is in response to the "redraw all" action function.
 *
 *   - RedrawReason_InteractiveDataViewChange:\n
 *     The draw event is in response to an interactive data view change such as
 *     rotate, translate, zoom, etc.
 *
 *   - RedrawReason_InteractivePaperViewChange:\n
 *     The draw event is in response to an interactive paper translate view or
 *     paper zoom view change.
 *
 *   - RedrawReason_InteractiveStyleChange:\n
 *     The draw event is in response to an interactive style changes such as
 *     dragging a contour level or a slice.
 *
 *   - RedrawReason_Animation:\n
 *     The draw event is in response to an animation.
 *
 *   - RedrawReason_AutoRedraw:\n
 *     The draw event is in response to an auto redraw.
 *
 *   - RedrawReason_RedrawForcedViewUpdate:\n
 *     The draw event is in response to forced view update when auto redraw is
 *     off such as a view fit or movement of the frame.
 *
 *   - RedrawReason_RedrawForcedStyleUpdate:\n
 *     The draw event is in response to forced view update when auto redraw is
 *     off such as deleting a contour level.
 *
 *   - RedrawReason_PreFullRedrawTraceOfAllFrames:\n
 *     The draw event is an approximate redraw done prior to a full redraw.
 *
 * @sa TecUtilEventAddPreDrawCallback(), TecUtilEventAddPostDrawCallback()
 */
typedef enum
{
    RedrawReason_UserReqRedrawActiveFrame,
    RedrawReason_UserReqTraceActiveFrame,
    RedrawReason_UserReqRedrawAllFrames,
    RedrawReason_UserReqTraceAllFrames,
    RedrawReason_InteractiveDataViewChange,
    RedrawReason_InteractivePaperViewChange,
    RedrawReason_InteractiveStyleChange,
    RedrawReason_Animation,
    RedrawReason_AutoRedraw,
    RedrawReason_RedrawForcedViewUpdate,
    RedrawReason_RedrawForcedStyleUpdate,
    RedrawReason_PreFullRedrawTraceOfAllFrames,
    END_RedrawReason_e,
    RedrawReason_Invalid = BadEnumValue,
    RedrawReason_UserReqRedrawCurrentFrame = RedrawReason_UserReqRedrawActiveFrame,
    RedrawReason_UserReqTraceCurrentFrame = RedrawReason_UserReqTraceActiveFrame
} RedrawReason_e;

typedef enum
{
    RotationMode_XYZAxis,
    RotationMode_Spherical,
    RotationMode_RollerBall,
    END_RotationMode_e,
    RotationMode_Invalid = BadEnumValue
} RotationMode_e;

typedef enum
{
    RotateAxis_X,
    RotateAxis_Y,
    RotateAxis_Z,
    RotateAxis_Psi,
    RotateAxis_Theta,
    RotateAxis_Alpha,
    RotateAxis_Twist,
    RotateAxis_VertRollerBall,
    RotateAxis_HorzRollerBall,
    RotateAxis_AboutVector,
    /* BEGINREMOVEFROMADDON */
    RotateAxis_DontCare, /* internal use only */
    /* ENDREMOVEFROMADDON */
    END_RotateAxis_e,
    RotateAxis_Invalid = BadEnumValue
} RotateAxis_e;

typedef enum
{
    RotateOriginLocation_DefinedOrigin,
    RotateOriginLocation_Viewer,
    END_RotateOriginLocation_e,
    RotateOriginLocation_Invalid = BadEnumValue
} RotateOriginLocation_e;

/*
 * NOTE: This is only used with the $!Reset3DOrigin command.
 */
typedef enum
{
    OriginResetLocation_DataCenter,
    OriginResetLocation_ViewCenter,
    END_OriginResetLocation_e,
    OriginResetLocation_Invalid = BadEnumValue
} OriginResetLocation_e;

/*
 * NOTE: This is only used with the $!CreateSliceZoneFromPlane command.
 */
typedef enum
{
    SliceSource_SurfaceZones,
    SliceSource_VolumeZones,
    SliceSource_SurfacesOfVolumeZones,
    SliceSource_LinearZones,
    END_SliceSource_e,
    SliceSource_Invalid = BadEnumValue
} SliceSource_e;





typedef enum
{
    Input_SmInteger,
    Input_Short,
    Input_Integer,
    Input_Float,
    Input_Double,
    Input_Radians,
    Input_TimeDateDouble,
    Input_ElapsedTimeDouble,
    END_Input_e,
    Input_Invalid = BadEnumValue
} Input_e;



typedef enum
{
    PtSelection_All,
    PtSelection_NearestN,
    PtSelection_OctantN,
    END_PtSelection_e,
    PtSelection_Invalid = BadEnumValue
} PtSelection_e;



typedef enum
{
    Drift_None,
    Drift_Linear,
    Drift_Quad,
    END_Drift_e,
    Drift_Invalid = BadEnumValue
} Drift_e;



/* atpoint is simple boundary condition.
   atpointb2 is better boundary condition.
*/
typedef enum
{
    DerivPos_atpoint,
    DerivPos_atpointb2,
    DerivPos_kphalf,
    DerivPos_jphalf,
    DerivPos_iphalf,
    END_DerivPos_e,
    DerivPos_Invalid = BadEnumValue
} DerivPos_e; /*<help>"atpoint is the simple boundary condition\n"*/
/*<help>"atpointb2 is a better boundary condition"*/


typedef enum
{
    LinearInterpMode_DontChange,
    LinearInterpMode_SetToConst,
    END_LinearInterpMode_e,
    LinearInterpMode_Invalid = BadEnumValue
} LinearInterpMode_e;

typedef enum
{
    VolumeCellInterpolationMode_PiecewiseLinear,
    VolumeCellInterpolationMode_TriLinear,
    END_VolumeCellInterpolationMode_e,
    VolumeCellInterpolationMode_Invalid = BadEnumValue
} VolumeCellInterpolationMode_e;

typedef enum
{
    PolyCellInterpolationMode_UseCCValue,
    PolyCellInterpolationMode_AverageNodes,
    END_PolyCellInterpolationMode_e,
    PolyCellInterpolationMode_Invalid = BadEnumValue
} PolyCellInterpolationMode_e;

typedef enum
{
    ConstraintOp2Mode_UseVar,
    ConstraintOp2Mode_UseConstant,
    END_ConstraintOp2Mode_e,
    ConstraintOp2Mode_Invalid = BadEnumValue
} ConstraintOp2Mode_e;

/**
 * Controls how data is loaded for interactive probe events.
 * DataProbeVarLoadMode_IncrementallyLoadAll will load as much data as possible within
 * load-on-demand time/space thresholds. DataProbeVarLoadMode_LoadRequiredVarsOnly will
 * load only the data necessary to complete the probe request.
 * DataProbeVarLoadMode_IncrementallyLoadAll is the default.
 */
typedef enum
{
    DataProbeVarLoadMode_IncrementallyLoadAll,
    DataProbeVarLoadMode_LoadRequiredVarsOnly,
    END_DataProbeVarLoadMode_e,
    DataProbeVarLoadMode_Invalid = BadEnumValue
} DataProbeVarLoadMode_e;

typedef enum
{
    ValueBlankCellMode_AllCorners,
    ValueBlankCellMode_AnyCorner,
    ValueBlankCellMode_PrimaryValue,
    END_ValueBlankCellMode_e,
    ValueBlankCellMode_Invalid = BadEnumValue,
    /* deprecated values */
    ValueBlankCellMode_PrimaryCorner = ValueBlankCellMode_PrimaryValue
} ValueBlankCellMode_e;


/*
 * deprecated: ValueBlankMode_e enumeration will not be supported after
 *             version 8. This API was retained for add-on developers
 *             using the TecUtilStyleSetLowLevel API.
 */
typedef enum
{
    ValueBlankMode_AndRule,
    ValueBlankMode_OrRule,
    ValueBlankMode_CornerRule,
    END_ValueBlankMode_e,
    ValueBlankMode_Invalid = BadEnumValue
} ValueBlankMode_e; /*<help>"DEPRECATED: ValueBlankMode_e will not be supported after version 8"*/


typedef enum
{
    CellBlankedCond_NotBlanked,
    CellBlankedCond_PartiallyBlanked,
    CellBlankedCond_EntirelyBlanked,
    CellBlankedCond_Uncertain,
    END_CellBlankedCond_e,
    CellBlankedCond_Invalid = BadEnumValue
} CellBlankedCond_e;


typedef enum
{
    RelOp_LessThanOrEqual,
    RelOp_GreaterThanOrEqual,
    RelOp_LessThan,
    RelOp_GreaterThan,
    RelOp_EqualTo,
    RelOp_NotEqualTo,
    END_RelOp_e,
    RelOp_Invalid = BadEnumValue
} RelOp_e;



typedef enum
{
    IJKBlankMode_BlankInterior,
    IJKBlankMode_BlankExterior,
    END_IJKBlankMode_e,
    IJKBlankMode_Invalid = BadEnumValue
} IJKBlankMode_e;


typedef enum
{
    PlotApproximationMode_Automatic,
    PlotApproximationMode_NonCurrentAlwaysApproximated,
    PlotApproximationMode_AllFramesAlwaysApproximated,
    END_PlotApproximationMode_e,
    PlotApproximationMode_Invalid = BadEnumValue
} PlotApproximationMode_e;

typedef enum
{
    SphereScatterRenderQuality_Low,
    SphereScatterRenderQuality_Medium,
    SphereScatterRenderQuality_High,
    END_SphereScatterRenderQuality_e,
    SphereScatterRenderQuality_Invalid = BadEnumValue
} SphereScatterRenderQuality_e;

/*
 * NOTE: FillPat_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by Translucency_e.
 */
typedef enum
{
    Pattern_Solid,
    Pattern_LowTranslucent,
    Pattern_MedTranslucent,
    Pattern_HighTranslucent,
    END_FillPat_e,
    Pattern_Invalid = BadEnumValue
} FillPat_e; /*<help>"DEPRECATED: Replaced by Translucency_e"*/


typedef enum
{
    Translucency_Solid,
    Translucency_Low,
    Translucency_Medium,
    Translucency_High,
    END_Translucency_e,
    Translucency_Invalid = BadEnumValue
} Translucency_e;



typedef enum
{
    SunRaster_OldFormat,
    SunRaster_Standard,
    SunRaster_ByteEncoded,
    END_SunRaster_e,
    SunRaster_Invalid = BadEnumValue
} SunRaster_e;


typedef enum
{
    BoundaryCondition_Fixed,
    BoundaryCondition_ZeroGradient,
    BoundaryCondition_Zero2nd,
    END_BoundaryCondition_e,
    BoundaryCondition_Invalid = BadEnumValue
} BoundaryCondition_e;



/* Note:
 *   In 2D: AxisMode_Independent and AxisMode_XYDependent are used;
 *   in 3D: AxisMode_Independent, AxisMode_XYZDependent, and AxisMode_XYDependent are used.
 */
typedef enum
{
    AxisMode_Independent,
    AxisMode_XYZDependent,
    AxisMode_XYDependent,
    END_AxisMode_e,
    AxisMode_Invalid = BadEnumValue
} AxisMode_e;/*<help>"In 2D AxisMode_Independent and AxisMode_XYDependent are used\n"*/
/*<help>"In 3D AxisMode_Independent, "*/
/*<help>"AxisMode_XYZDependent, and AxisMode_XYDependent are used."*/

typedef enum
{
    Quick_LineColor,
    Quick_FillColor,
    Quick_TextColor,
    END_QuickColorMode_e,
    Quick_Invalid = BadEnumValue
} QuickColorMode_e;


typedef enum
{
    FillMode_None,
    FillMode_UseSpecificColor,
    FillMode_UseLineColor,
    FillMode_UseBackgroundColor,
    END_FillMode_e,
    FillMode_Invalid = BadEnumValue
} FillMode_e;


typedef enum
{
    LinePattern_Solid,
    LinePattern_Dashed,
    LinePattern_DashDot,
    LinePattern_Dotted,
    LinePattern_LongDash,
    LinePattern_DashDotDot,
    END_LinePattern_e,
    LinePattern_Invalid = BadEnumValue
} LinePattern_e;



typedef enum
{
    Join_Miter,
    Join_Round,
    Join_Bevel,
    END_LineJoin_e,
    Join_Invalid = BadEnumValue
} LineJoin_e;



typedef enum
{
    Cap_Flat,
    Cap_Round,
    Cap_Square,
    END_LineCap_e,
    Cap_Invalid = BadEnumValue
} LineCap_e;



typedef enum
{
    GeomForm_LineSegs,
    GeomForm_Rectangle,
    GeomForm_Square,
    GeomForm_Circle,
    GeomForm_Ellipse,
    GeomForm_LineSegs3D, /* deprecated: use GeomForm_LineSegs with CoordSys_Grid3D */
    GeomForm_Image,
    END_GeomForm_e,
    GeomForm_Invalid = BadEnumValue,
    /* new value names */
    GeomType_LineSegs = GeomForm_LineSegs,
    GeomType_Rectangle = GeomForm_Rectangle,
    GeomType_Square = GeomForm_Square,
    GeomType_Circle = GeomForm_Circle,
    GeomType_Ellipse = GeomForm_Ellipse,
    GeomType_LineSegs3D = GeomForm_LineSegs3D, /* deprecated: use GeomType_LineSegs with CoordSys_Grid3D */
    GeomType_Image = GeomForm_Image,
    END_GeomType_e = END_GeomForm_e,
    GeomType_Invalid = GeomForm_Invalid
} GeomForm_e;

typedef GeomForm_e GeomType_e;

typedef enum
{
    VariableDerivationMethod_Fast,
    VariableDerivationMethod_Accurate,
    END_VariableDerivationMethod_e,
    VariableDerivationMethod_Invalid = BadEnumValue
} VariableDerivationMethod_e;

/**
 */
typedef enum
{
    AuxDataType_String,
    END_AuxDataType_e,
    AuxDataType_Invalid = BadEnumValue
} AuxDataType_e;

/**
 */
typedef enum
{
    AuxDataLocation_Zone,
    AuxDataLocation_DataSet,
    AuxDataLocation_Frame,
    AuxDataLocation_Var,
    AuxDataLocation_LineMap,
    AuxDataLocation_Page,
    END_AuxDataLocation_e,
    AuxDataLocation_Invalid = BadEnumValue
} AuxDataLocation_e;


/* Note: This replaces Element_e */
typedef enum
{
    ZoneType_Ordered,
    ZoneType_FETriangle,
    ZoneType_FEQuad,
    ZoneType_FETetra,
    ZoneType_FEBrick,
    ZoneType_FELineSeg,
    ZoneType_FEPolygon,
    ZoneType_FEPolyhedron,
    END_ZoneType_e,
    ZoneType_Invalid = BadEnumValue
} ZoneType_e;

typedef enum
{
    ZoneOrder_I,
    ZoneOrder_J,
    ZoneOrder_K,
    ZoneOrder_IJ,
    ZoneOrder_IK,
    ZoneOrder_JK,
    ZoneOrder_IJK,
    END_ZoneOrder_e,
    ZoneOrder_Invalid = BadEnumValue
} ZoneOrder_e;

/* deprecated: replaced by ZoneType_e DataPacking_e */
typedef enum
{
    DataFormat_IJKBlock,
    DataFormat_IJKPoint,
    DataFormat_FEBlock,
    DataFormat_FEPoint,
    END_DataFormat_e,
    DataFormat_Invalid = BadEnumValue
} DataFormat_e;

typedef enum
{
    DataPacking_Block,
    DataPacking_Point,
    END_DataPacking_e,
    DataPacking_Invalid = BadEnumValue
} DataPacking_e;



typedef enum
{
    PD_HPGL,
    PD_HPGL2,
    PD_PS,
    PD_LASERG, /* deprecated */
    PD_EPS,
    PD_WINDOWS, /* Windows Print Driver */
    PD_WMF, /* Windows MetaFile (used from Export only) */
    PD_X3D,
    END_PrinterDriver_e,
    PD_Invalid = BadEnumValue
} PrinterDriver_e;



typedef enum
{
    Image_None,
    Image_TIFF,
    Image_EPSI2,
    Image_FRAME,
    END_EPSPreviewImage_e,
    Image_Invalid = BadEnumValue
} EPSPreviewImage_e;

typedef enum
{
    TIFFByteOrder_Intel,
    TIFFByteOrder_Motorola,
    END_TIFFByteOrder_e,
    TIFFByteOrder_Invalid = BadEnumValue
} TIFFByteOrder_e;

typedef enum
{
    JPEGEncoding_Standard,
    JPEGEncoding_Progressive,
    END_JPEGEncoding_e,
    JPEGEncoding_Invalid = BadEnumValue
} JPEGEncoding_e;


typedef enum
{
    FlashImageType_Lossless,
    FlashImageType_JPEG,
    FlashImageType_Color256,
    END_FlashImageType_e,
    FlashImageType_Invalid = BadEnumValue,
    /* deprecated values */
    FlashImageType_256Color = FlashImageType_Color256
} FlashImageType_e;

typedef enum
{
    FlashCompressionType_BestSpeed,
    FlashCompressionType_SmallestSize,
    END_FlashCompressionType_e,
    FlashCompressionType_Invalid = BadEnumValue
} FlashCompressionType_e;


typedef enum
{
    ExportFormat_RasterMetafile,
    ExportFormat_TIFF,
    ExportFormat_SGI,
    ExportFormat_SunRaster,
    ExportFormat_XWindows,
    ExportFormat_PSImage,       /* deprecated */
    ExportFormat_HPGL,
    ExportFormat_HPGL2,
    ExportFormat_PS,
    ExportFormat_EPS,
    ExportFormat_LaserGraphics, /* deprecated */
    ExportFormat_WindowsMetafile,
    ExportFormat_BMP,
    ExportFormat_PNG,
    ExportFormat_AVI,
    ExportFormat_Custom,  /* May be used in a future version */
    ExportFormat_JPEG,
    ExportFormat_Flash,
    ExportFormat_X3D,
    ExportFormat_TecplotViewer,
    ExportFormat_FLV,
    ExportFormat_MPEG4,
    ExportFormat_WMV,
    END_ExportFormat_e,
    ExportFormat_Invalid = BadEnumValue
} ExportFormat_e;

typedef enum
{
    AVICompression_ColorPreserving,
    AVICompression_LinePreserving,
    AVICompression_LosslessUncompressed,
    END_AVICompression_e,
    AVICompression_Invalid = BadEnumValue
} AVICompression_e;

typedef enum
{
    AnimationDest_Screen,
    AnimationDest_AVI,
    AnimationDest_RM,
    AnimationDest_Flash,
    AnimationDest_FLV,
    AnimationDest_MPEG4,
    AnimationDest_WMV,
    END_AnimationDest_e,
    AnimationDest_Invalid = BadEnumValue
} AnimationDest_e;



typedef enum
{
    AnimationOperation_Forward,
    AnimationOperation_Backward,
    AnimationOperation_Loop,
    AnimationOperation_Bounce,
    END_AnimationOperation_e,
    AnimationOperation_Invalid = BadEnumValue
} AnimationOperation_e;

typedef enum
{
    AnimationStep_First,
    AnimationStep_Second,
    AnimationStep_Current,
    AnimationStep_SecondToLast,
    AnimationStep_Last,
    AnimationStep_Previous,
    AnimationStep_Next,
    END_AnimationStep_e,
    AnimationStep_Invalid = BadEnumValue
} AnimationStep_e;

typedef enum
{
    ZoneAnimationMode_StepByNumber,
    ZoneAnimationMode_GroupStepByNumber,
    ZoneAnimationMode_StepByTime,
    END_ZoneAnimationMode_e,
    ZoneAnimationMode_Invalid = BadEnumValue
} ZoneAnimationMode_e;

#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref ExportRegion_e instead.
 */
typedef enum
{
    BitDumpRegion_CurrentFrame,
    BitDumpRegion_AllFrames,
    BitDumpRegion_WorkArea,
    END_BitDumpRegion_e,
    BitDumpRegion_Invalid = BadEnumValue
} BitDumpRegion_e;
#endif

typedef enum
{
    ExportRegion_CurrentFrame,
    ExportRegion_AllFrames,
    ExportRegion_WorkArea,
    END_ExportRegion_e,
    ExportRegion_Invalid = BadEnumValue
} ExportRegion_e;

typedef enum
{
    PaperSize_Letter,
    PaperSize_Double,
    PaperSize_A4,
    PaperSize_A3,
    PaperSize_Custom1,
    PaperSize_Custom2,
    END_PaperSize_e,
    PaperSize_Invalid = BadEnumValue,
    /* deprecated values */
    Paper_Letter  = PaperSize_Letter,  /* deprecated */
    Paper_Double  = PaperSize_Double,  /* deprecated */
    Paper_A4      = PaperSize_A4,      /* deprecated */
    Paper_A3      = PaperSize_A3,      /* deprecated */
    Paper_Custom1 = PaperSize_Custom1, /* deprecated */
    Paper_Custom2 = PaperSize_Custom2, /* deprecated */
    Paper_Invalid = PaperSize_Invalid  /* deprecated */
} PaperSize_e;



typedef enum
{
    PaperUnitSpacing_HalfCentimeter,
    PaperUnitSpacing_OneCentimeter,
    PaperUnitSpacing_TwoCentimeters,
    PaperUnitSpacing_QuarterInch,
    PaperUnitSpacing_HalfInch,
    PaperUnitSpacing_OneInch,
    PaperUnitSpacing_TenPoints,
    PaperUnitSpacing_TwentyFourPoints,
    PaperUnitSpacing_ThirtySixPoints,
    PaperUnitSpacing_FiftyPoints,
    PaperUnitSpacing_SeventyTwoPoints,
    PaperUnitSpacing_OneTenthInch,
    PaperUnitSpacing_OneTenthCentimeter,
    END_PaperUnitSpacing_e,
    PaperUnitSpacing_Invalid = BadEnumValue
} PaperUnitSpacing_e;


typedef enum
{
    Palette_Monochrome,
    Palette_PenPlotter,
    Palette_Color,
    END_Palette_e,
    Palette_Invalid = BadEnumValue
} Palette_e;


typedef enum
{
    PrintRenderType_Vector,
    PrintRenderType_Image,
    END_PrintRenderType_e,
    PrintRenderType_Invalid = BadEnumValue
} PrintRenderType_e;


typedef enum
{
    Units_Grid,
    Units_Frame,
    Units_Point,
    Units_Screen,
    Units_AxisPercentage,
    END_Units_e,
    Units_Invalid = BadEnumValue
} Units_e;


typedef enum
{
    CoordScale_Linear,
    CoordScale_Log,
    END_CoordScale_e,
    CoordScale_Invalid = BadEnumValue,
    /* old names for the same values */
    Scale_Linear = CoordScale_Linear,
    Scale_Log = CoordScale_Log,
    Scale_Invalid = CoordScale_Invalid
} CoordScale_e;

/* BEGINREMOVEFROMADDON */
#define GetLog10(R) ( ((R) < SMALLDOUBLE) ? SMALLESTEXPONENT : ( ((R) > LARGEDOUBLE) ? LARGESTEXPONENT : log10((R)) ) )
/* ENDREMOVEFROMADDON */

typedef enum
{
    CoordSys_Grid,
    CoordSys_Frame,
    CoordSys_FrameOffset,
    CoordSys_Paper,
    CoordSys_Screen,
    CoordSys_Hardcopy,
    CoordSys_Grid3D,
    CoordSys_Workspace,
    END_CoordSys_e,
    CoordSys_Invalid = BadEnumValue
} CoordSys_e;

/*
 *  NOTE:  CoordSys_FrameOffset always is stored in inches internally.
 *         in stylesheet this may be written in other units if
 *         appropriate suffix is added.
 *
 */



typedef enum
{
    Scope_Global,
    Scope_Local,
    END_Scope_e,
    Scope_Invalid = BadEnumValue
} Scope_e;


typedef enum
{
    TextAnchor_Left,
    TextAnchor_Center,
    TextAnchor_Right,
    TextAnchor_MidLeft,
    TextAnchor_MidCenter,
    TextAnchor_MidRight,
    TextAnchor_HeadLeft,
    TextAnchor_HeadCenter,
    TextAnchor_HeadRight,
    TextAnchor_OnSide,
    END_TextAnchor_e,
    TextAnchor_Invalid = BadEnumValue
} TextAnchor_e;



typedef enum
{
    TextBox_None,
    TextBox_Filled,
    TextBox_Hollow,
    END_TextBox_e,
    TextBox_Invalid = BadEnumValue
} TextBox_e;



typedef enum
{
    GeomShape_Square,
    GeomShape_Del,
    GeomShape_Grad,
    GeomShape_RTri,
    GeomShape_LTri,
    GeomShape_Diamond,
    GeomShape_Circle,
    GeomShape_Cube,
    GeomShape_Sphere,
    GeomShape_Octahedron,
    GeomShape_Point,
    END_GeomShape_e,
    GeomShape_Invalid = BadEnumValue
} GeomShape_e;


typedef enum
{
    BasicSize_Tiny,
    BasicSize_Small,
    BasicSize_Medium,
    BasicSize_Large,
    BasicSize_Huge,
    END_BasicSize_e,
    BasicSize_Invalid = BadEnumValue
} BasicSize_e;



/*
 * NOTE: LineForm_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by CurveType_e.
 */
typedef enum
{
    LineForm_LineSeg,
    LineForm_CurvFit,
    LineForm_EToRFit,
    LineForm_PowerFit,
    LineForm_Spline,
    LineForm_ParaSpline,
    END_LineForm_e,
    LineForm_Invalid = BadEnumValue
} LineForm_e;


typedef enum
{
    CurveType_LineSeg,
    CurveType_PolynomialFit,
    CurveType_EToRFit,
    CurveType_PowerFit,
    CurveType_Spline,
    CurveType_ParaSpline,
    CurveType_Extended,
    END_CurveType_e,
    CurveType_Invalid = BadEnumValue,
    CurveType_CurvFit = CurveType_PolynomialFit
} CurveType_e;

typedef enum
{
    Script_None,
    Script_Super,
    Script_Sub,
    END_Script_e,
    Script_Invalid = BadEnumValue
} Script_e;


typedef enum
{
    Font_Helvetica,
    Font_HelveticaBold,
    Font_Greek,
    Font_Math,
    Font_UserDefined,
    Font_Times,
    Font_TimesItalic,
    Font_TimesBold,
    Font_TimesItalicBold,
    Font_Courier,
    Font_CourierBold,
    Font_Extended,
    END_Font_e,
    Font_Invalid = BadEnumValue
} Font_e;

typedef enum
{
    FontStyle_Regular,
    FontStyle_Italic,
    FontStyle_Bold,
    FontStyle_BoldItalic,
    END_FontStyle_e,
    FontStyle_Invalid = BadEnumValue
} FontStyle_e;

typedef enum
{
    TwoDDrawOrder_ByZone,
    TwoDDrawOrder_ByLayer,
    END_TwoDDrawOrder_e,
    TwoDDrawOrder_Invalid = BadEnumValue
} TwoDDrawOrder_e;

typedef enum
{
    DrawOrder_AfterData,
    DrawOrder_BeforeData,
    END_DrawOrder_e,
    DrawOrder_Invalid = BadEnumValue
} DrawOrder_e;

/*
 *
 * NOTE: Streamtrace_TwoDLine is new.  All 2D
 *       streamtraces are assigned this value.
 */
typedef enum
{
    Streamtrace_SurfaceLine,
    Streamtrace_SurfaceRibbon,
    Streamtrace_VolumeLine,
    Streamtrace_VolumeRibbon,
    Streamtrace_VolumeRod,
    Streamtrace_TwoDLine,
    END_Streamtrace_e,
    Streamtrace_Invalid = BadEnumValue
} Streamtrace_e;



typedef enum
{
    StreamDir_Forward,
    StreamDir_Reverse,
    StreamDir_Both,
    END_StreamDir_e,
    StreamDir_Invalid = BadEnumValue
} StreamDir_e;

typedef enum
{
    DistributionRegion_Point,
    DistributionRegion_Rake,
    DistributionRegion_SurfacesOfActiveZones,
    DistributionRegion_SurfacesOfSelectedObjects,
    END_DistributionRegion_e,
    DistributionRegion_Invalid = BadEnumValue
} DistributionRegion_e;

typedef enum
{
    IsoSurfaceSelection_AllContourLevels,
    IsoSurfaceSelection_OneSpecificValue,
    IsoSurfaceSelection_TwoSpecificValues,
    IsoSurfaceSelection_ThreeSpecificValues,
    END_IsoSurfaceSelection_e,
    IsoSurfaceSelection_Invalid = BadEnumValue
} IsoSurfaceSelection_e;


typedef enum
{
    ValueLocation_CellCentered,
    ValueLocation_Nodal,
    END_ValueLocation_e,
    ValueLocation_Invalid = BadEnumValue
} ValueLocation_e;

typedef enum
{
    FieldDataType_Reserved, /* never use */
    FieldDataType_Float,
    FieldDataType_Double,
    FieldDataType_Int32,
    FieldDataType_Int16,
    FieldDataType_Byte,
    FieldDataType_Bit,
    END_FieldDataType_e,
    FieldDataType_IJKFunction,   /* Not used yet */
    FieldDataType_Int64, /* Not used yet */
#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
    FieldDataType_LongInt = FieldDataType_Int32,
    FieldDataType_ShortInt = FieldDataType_Int16,
#endif
    FieldDataType_Invalid = BadEnumValue
} FieldDataType_e;

#define VALID_FIELD_DATA_TYPE(FieldDataType) (VALID_ENUM((FieldDataType),FieldDataType_e) && \
                                              (FieldDataType)!=FieldDataType_Reserved)

#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref MeshType_e instead.
 */
typedef enum
{
    Mesh_Wireframe,  /* deprecated: use MeshType_Wireframe  */
    Mesh_Overlay,    /* deprecated: use MeshType_Overlay    */
    Mesh_HiddenLine, /* deprecated: use MeshType_HiddenLine */
    END_MeshPlotType_e,
    Mesh_Invalid = BadEnumValue
} MeshPlotType_e;
#endif

typedef enum
{
    MeshType_Wireframe,  /* Mesh_Wireframe  */
    MeshType_Overlay,    /* Mesh_Overlay    */
    MeshType_HiddenLine, /* Mesh_HiddenLine */
    END_MeshType_e,
    MeshType_Invalid = BadEnumValue
} MeshType_e;




#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref ContourType_e instead.
 */
typedef enum
{
    Contour_Lines,       /* deprecated: use ContourType_Lines        */
    Contour_Flood,       /* deprecated: use ContourType_Flood        */
    Contour_Overlay,     /* deprecated: use ContourType_Overlay      */
    Contour_AverageCell, /* deprecated: use ContourType_AverageCell  */
    Contour_CornerCell,  /* deprecated: use ContourType_PrimaryValue */
    END_ContourPlotType_e,
    Contour_Invalid = BadEnumValue
} ContourPlotType_e;
#endif


typedef enum
{
    ContourType_Lines,         /* Contour_Lines       */
    ContourType_Flood,         /* Contour_Flood       */
    ContourType_Overlay,       /* Contour_Overlay     */
    ContourType_AverageCell,   /* Contour_AverageCell */
    ContourType_PrimaryValue,  /* Contour_CornerCell  */
    END_ContourType_e,
    ContourType_Invalid = BadEnumValue
} ContourType_e;

typedef enum
{
    ContourColoring_RGB,
    ContourColoring_Group1,
    ContourColoring_Group2,
    ContourColoring_Group3,
    ContourColoring_Group4,
    ContourColoring_Group5,
    ContourColoring_Group6,
    ContourColoring_Group7,
    ContourColoring_Group8,
    END_ContourColoring_e,
    ContourColoring_Invalid = BadEnumValue
} ContourColoring_e;

#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref VectorType_e instead.
 */
typedef enum
{
    Vector_TailAtPoint, /* deprecated: use VectorType_TailAtPoint */
    Vector_HeadAtPoint, /* deprecated: use VectorType_HeadAtPoint */
    Vector_MidAtPoint,  /* deprecated: use VectorType_MidAtPoint  */
    Vector_HeadOnly,    /* deprecated: use VectorType_HeadOnly    */
    END_VectorPlotType_e,
    Vector_Invalid = BadEnumValue
} VectorPlotType_e;
#endif


typedef enum
{
    VectorType_TailAtPoint, /* Vector_TailAtPoint */
    VectorType_HeadAtPoint, /* Vector_HeadAtPoint */
    VectorType_MidAtPoint,  /* Vector_MidAtPoint  */
    VectorType_HeadOnly,    /* Vector_HeadOnly    */
    END_VectorType_e,
    VectorType_Invalid = BadEnumValue
} VectorType_e;


/*
 * NOTE: ShadePlotType_e is deprecated.  It must be retained to maintain
 *       backward compatibility with the TecUtil layer however.
 *       This has been replaced by LightingEffect_e.
 */
typedef enum
{
    Shade_SolidColor,
    Shade_Paneled,
    Shade_Gouraud,
    Shade_ColoredPaneled,
    Shade_ColoredGouraud,
    END_ShadePlotType_e,
    Shade_Invalid = BadEnumValue
} ShadePlotType_e;

/*
 * NOTE: LightingEffect_None is deprecated.  It must remain
 *       in the list to allow macro processing of older
 *       (i.e. early v9) macros.
 */
typedef enum
{
    LightingEffect_Paneled,
    LightingEffect_Gouraud,
    LightingEffect_None,
    END_LightingEffect_e,
    LightingEffect_Invalid = BadEnumValue
} LightingEffect_e;

typedef enum
{
    IJKLines_I,
    IJKLines_J,
    IJKLines_K,
    END_IJKLines_e,
    IJKLines_Invalid = BadEnumValue,
    /* deprecated values */
    Lines_I       = IJKLines_I,       /* deprecated */
    Lines_J       = IJKLines_J,       /* deprecated */
    Lines_K       = IJKLines_K,       /* deprecated */
    Lines_Invalid = IJKLines_Invalid  /* deprecated */
} IJKLines_e;

typedef enum
{
    IJKCellType_Planes,
    IJKCellType_FacePlanes,
    IJKCellType_Volume,
    END_IJKCellType_e,
    IJKCellType_Invalid = BadEnumValue
} IJKCellType_e;


/*
 *  Ver 6 used PlaneSet.  Ver 7 uses CellType and Planes variables.
 *
 *   "PlaneSet" in version 6    vs.  IJKPlanes in v7:
 *
 *   'A' = AllPlanes                 CellType = IJKCellType_Volume
 *   'd','e','f','C' = ComboPlanes   CellType = IJKCellType_Planes, IJKPlanes = depends on defC
 *   'F' = Faces Planes Only         CellType = IJKCellType_FacePlanes
 *   'I' = I-Planes                  CellType = IJKCellType_Planes, IJKPlanes = IJKPlanes_I
 *   'J' = J-Planes                  CellType = IJKCellType_Planes, IJKPlanes = IJKPlanes_J
 *   'K' = K-Planes                  CellType = IJKCellType_Planes, IJKPlanes = IJKPlanes_K
 *
 *
 * NOTE: IJKPlanes_e is still used internally in tecplot (and in the TecUtil layer).
 *       it has been relagated to communicating which planes of an IJK zone are in
 *       use.
 *
 */

typedef enum
{
    IJKPlanes_I,
    IJKPlanes_J,
    IJKPlanes_K,
    IJKPlanes_Face,  /* used on the panel heap */
    IJKPlanes_IJ,    /* deprecated */
    IJKPlanes_JK,    /* deprecated */
    IJKPlanes_IK,    /* deprecated */
    IJKPlanes_IJK,   /* deprecated */
    IJKPlanes_Volume,
    IJKPlanes_Unused,
    END_IJKPlanes_e,
    IJKPlanes_Invalid = BadEnumValue,
    /* deprecated values */
    Planes_I       = IJKPlanes_I,      /* deprecated */
    Planes_J       = IJKPlanes_J,      /* deprecated */
    Planes_K       = IJKPlanes_K,      /* deprecated */
    Planes_IJ      = IJKPlanes_IJ,     /* deprecated */
    Planes_JK      = IJKPlanes_JK,     /* deprecated */
    Planes_IK      = IJKPlanes_IK,     /* deprecated */
    Planes_IJK     = IJKPlanes_IJK,    /* deprecated */
    Planes_Face    = IJKPlanes_Face,   /* deprecated */
    Planes_Volume  = IJKPlanes_Volume, /* deprecated */
    Planes_Unused  = IJKPlanes_Unused, /* deprecated */
    Planes_Invalid = IJKPlanes_Invalid /* deprecated */
} IJKPlanes_e;



typedef enum
{
    SurfacesToPlot_BoundaryFaces,
    SurfacesToPlot_ExposedCellFaces,
    SurfacesToPlot_IPlanes,
    SurfacesToPlot_JPlanes,
    SurfacesToPlot_KPlanes,
    SurfacesToPlot_IJPlanes,
    SurfacesToPlot_JKPlanes,
    SurfacesToPlot_IKPlanes,
    SurfacesToPlot_IJKPlanes,
    SurfacesToPlot_All,
    SurfacesToPlot_None,
    END_SurfacesToPlot_e,
    SurfacesToPlot_Invalid = BadEnumValue
} SurfacesToPlot_e;

typedef enum
{
    PointsToPlot_SurfaceNodes,  /* was _SurfacesOnly */
    PointsToPlot_AllNodes,      /* was _All          */
    PointsToPlot_SurfaceCellCenters,
    PointsToPlot_AllCellCenters,
    PointsToPlot_AllConnected,
    END_PointsToPlot_e,
#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
    PointsToPlot_SurfacesOnly = PointsToPlot_SurfaceNodes, /* deprecated */
    PointsToPlot_All          = PointsToPlot_AllNodes,     /* deprecated */
#endif
    PointsToPlot_Invalid = BadEnumValue
} PointsToPlot_e;


typedef enum
{
    SliceSurface_XPlanes,
    SliceSurface_YPlanes,
    SliceSurface_ZPlanes,
    SliceSurface_IPlanes,
    SliceSurface_JPlanes,
    SliceSurface_KPlanes,
    END_SliceSurface_e,
    SliceSurface_Invalid = BadEnumValue
} SliceSurface_e;


typedef enum
{
    ClipPlane_None,
    ClipPlane_BelowPrimarySlice,
    ClipPlane_AbovePrimarySlice,
    END_ClipPlane_e,
    ClipPlane_Invalid = BadEnumValue
} ClipPlane_e;

typedef enum
{
    Skip_ByIndex,
    Skip_ByFrameUnits,
    END_SkipMode_e,
    Skip_Invalid = BadEnumValue
} SkipMode_e;


typedef enum
{
    EdgeType_Borders,
    EdgeType_Creases,
    EdgeType_BordersAndCreases,
    END_EdgeType_e,
    EdgeType_Invalid = BadEnumValue
} EdgeType_e;

#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref BorderLocation_e instead.
 */
typedef enum
{
    Boundary_None, /* deprecated: use BoundaryType_None */
    Boundary_Min,  /* deprecated: use BoundaryType_Min  */
    Boundary_Max,  /* deprecated: use BoundaryType_Max  */
    Boundary_Both, /* deprecated: use BoundaryType_Both */
    END_BoundPlotType_e,
    Boundary_Invalid = BadEnumValue
} BoundPlotType_e;
#endif

typedef enum
{
    BoundaryType_None, /* Boundary_None */
    BoundaryType_Min,  /* Boundary_Min  */
    BoundaryType_Max,  /* Boundary_Max  */
    BoundaryType_Both, /* Boundary_Both */
    END_BoundaryType_e,
    BoundaryType_Invalid = BadEnumValue
} BoundaryType_e;  /* deprecated */

typedef enum
{
    BorderLocation_None, /* Boundary_None */
    BorderLocation_Min,  /* Boundary_Min  */
    BorderLocation_Max,  /* Boundary_Max  */
    BorderLocation_Both, /* Boundary_Both */
    END_BorderLocation_e,
    BorderLocation_Invalid = BadEnumValue
} BorderLocation_e;

typedef enum
{
    ContourColorMap_SmRainbow,
    ContourColorMap_LgRainbow,
    ContourColorMap_Modern,
    ContourColorMap_GrayScale,
    ContourColorMap_Wild,
    ContourColorMap_UserDef,
    ContourColorMap_TwoColor,
    ContourColorMap_RawUserDef,
    ContourColorMap_DivBuRd,
    ContourColorMap_DivBuYlRd,
    ContourColorMap_DivBrBG,
    ContourColorMap_DivOrPu,
    ContourColorMap_DivPiYG,
    ContourColorMap_DivPRGn,
    ContourColorMap_Doppler,
    ContourColorMap_ElevAboveGrnd,
    ContourColorMap_ElevAbsolute,
    ContourColorMap_HotMetal,
    ContourColorMap_Magma,
    ContourColorMap_DkRainbow,
    ContourColorMap_MdRainbow,
    ContourColorMap_QualAccent,
    ContourColorMap_QualDark1,
    ContourColorMap_QualDark2,
    ContourColorMap_QualPaired,
    ContourColorMap_QualPastel1,
    ContourColorMap_QualPastel2,
    ContourColorMap_QualPastel3,
    ContourColorMap_SeqBlue,
    ContourColorMap_SeqBuGn,
    ContourColorMap_SeqBuPu,
    ContourColorMap_SeqGreen,
    ContourColorMap_SeqGnBu,
    ContourColorMap_SeqOrange,
    ContourColorMap_SeqOrRd,
    ContourColorMap_SeqPiPu,
    ContourColorMap_SeqPurple,
    ContourColorMap_SeqPuBu,
    ContourColorMap_SeqPuBuGn,
    ContourColorMap_SeqPuRd,
    ContourColorMap_SeqRed,
    ContourColorMap_SeqYlGn,
    ContourColorMap_SeqYlGnBu,
    ContourColorMap_SeqYlOrBr,
    ContourColorMap_SeqYlOrRd,
    END_ContourColorMap_e,
    ContourColorMap_Invalid = BadEnumValue,
    /* deprecated values */
    ColorMap_SmRainbow    = ContourColorMap_SmRainbow,    /* deprecated */
    ColorMap_LgRainbow    = ContourColorMap_LgRainbow,    /* deprecated */
    ColorMap_Modern       = ContourColorMap_Modern,       /* deprecated */
    ColorMap_GrayScale    = ContourColorMap_GrayScale,    /* deprecated */
    ColorMap_Wild         = ContourColorMap_Wild,         /* deprecated */
    ColorMap_UserDef      = ContourColorMap_UserDef,      /* deprecated */
    ColorMap_TwoColor     = ContourColorMap_TwoColor,     /* deprecated */
    ColorMap_RawUserDef   = ContourColorMap_RawUserDef,   /* deprecated */
    ColorMap_Invalid      = ContourColorMap_Invalid       /* deprecated */
} ContourColorMap_e;



typedef enum
{
    ErrorBar_Up,
    ErrorBar_Down,
    ErrorBar_Left,
    ErrorBar_Right,
    ErrorBar_Horz,
    ErrorBar_Vert,
    ErrorBar_Cross,
    END_ErrorBar_e,
    ErrorBar_Invalid = BadEnumValue
} ErrorBar_e;



typedef enum
{
    ContourLineMode_UseZoneLineType,
    ContourLineMode_SkipToSolid,
    ContourLineMode_DashNegative,
    END_ContourLineMode_e,
    ContourLineMode_Invalid = BadEnumValue
} ContourLineMode_e;


/* BEGINREMOVEFROMADDON */
typedef enum
{
    Panel_Bad,
    Panel_Cell,                 /* FieldZone */
    Panel_Vector,               /* FieldZone */
    Panel_Scatter,              /* FieldZone */
    Panel_IJKBorderLine,        /* FieldZone IJK border lines */
    Panel_CellEdge,             /* FieldZone border lines and creases */
    Panel_FEBoundaryCell,       /* FieldZone */
    Panel_NodeLabel,            /* FieldZone */
    Panel_CellLabel,            /* FieldZone */
    Panel_StreamtraceCell,      /* Streamtrace COB          */
    Panel_StreamtraceMarker,    /* StreamtraceMarker COB (Scatter Symbol) */
    Panel_StreamtraceArrowhead, /* StreamtraceArrowhead COB (Vector) */
    Panel_IsoSurfaceCell,       /* IsoSurface COB */
    Panel_IsoSurfaceCellEdge,   /* IsoSurface COB border lines and creases (border lines and creases not currently used) */
    Panel_SliceCell,            /* Slice COB */
    Panel_SliceVector,          /* Slice COB */
    Panel_SliceIJKBorderLine,   /* Slice COB IJK border lines */
    Panel_SliceCellEdge,        /* Slice COB border lines and creases (creases not currently used) */
    Panel_Geom,                 /* Misc */
    Panel_Text,                 /* Misc */
    END_Panel_e,
    Panel_Invalid = BadEnumValue
} Panel_e;
/* ENDREMOVEFROMADDON */


typedef enum
{
    MessageBoxType_Error,
    MessageBoxType_Warning,
    MessageBoxType_Information,
    MessageBoxType_Question,   /* Ok, Cancel buttons */
    MessageBoxType_YesNo,
    MessageBoxType_YesNoCancel,
    MessageBoxType_WarningOkCancel,
    END_MessageBoxType_e,
    MessageBoxType_Invalid = BadEnumValue,
    /* deprecated values */
    MessageBox_Error           = MessageBoxType_Error,           /* deprecated */
    MessageBox_Warning         = MessageBoxType_Warning,         /* deprecated */
    MessageBox_Information     = MessageBoxType_Information,     /* deprecated */
    MessageBox_Question        = MessageBoxType_Question,        /* deprecated */
    MessageBox_YesNo           = MessageBoxType_YesNo,           /* deprecated */
    MessageBox_YesNoCancel     = MessageBoxType_YesNoCancel,     /* deprecated */
    MessageBox_WarningOkCancel = MessageBoxType_WarningOkCancel, /* deprecated */
    MessageBox_Invalid         = MessageBoxType_Invalid          /* deprecated */
} MessageBoxType_e;


typedef enum
{
    MessageBoxReply_Yes,
    MessageBoxReply_No,
    MessageBoxReply_Cancel,
    MessageBoxReply_Ok,
    END_MessageBoxReply_e,
    MessageBoxReply_Invalid = BadEnumValue
} MessageBoxReply_e;

typedef enum
{
    NumberFormat_Integer,
    NumberFormat_FixedFloat,
    NumberFormat_Exponential,
    NumberFormat_BestFloat,
    NumberFormat_SuperScript,
    NumberFormat_CustomLabel,
    NumberFormat_LogSuperScript,
    NumberFormat_RangeBestFloat,
    NumberFormat_DynamicLabel,
    NumberFormat_TimeDate,
    END_NumberFormat_e,
    NumberFormat_Invalid = BadEnumValue
} NumberFormat_e;

/* For backward compatibility with v9- */
typedef NumberFormat_e ValueFormat_e;


typedef enum
{
    BackingStoreMode_QuickAndDirty,
    BackingStoreMode_RealTimeUpdate,
    BackingStoreMode_PeriodicUpdate,
    END_BackingStoreMode_e,
    BackingStoreMode_Invalid = BadEnumValue
} BackingStoreMode_e;


typedef enum
{
    TickDirection_In,
    TickDirection_Out,
    TickDirection_Centered,
    END_TickDirection_e,
    TickDirection_Invalid = BadEnumValue
} TickDirection_e;

/* This enumerated type is no longer used as of Tecplot V10. */
typedef enum
{
    AxisTitlePosition_Left,
    AxisTitlePosition_Center,
    AxisTitlePosition_Right,
    END_AxisTitlePosition_e,
    AxisTitlePosition_Invalid = BadEnumValue
} AxisTitlePosition_e;

typedef enum
{
    AxisTitleMode_NoTitle,
    AxisTitleMode_UseVarName,
    AxisTitleMode_UseText,
    END_AxisTitleMode_e,
    AxisTitleMode_Invalid = BadEnumValue
} AxisTitleMode_e;

typedef enum
{
    AxisAlignment_WithViewport,
    AxisAlignment_WithOpposingAxisValue,
    AxisAlignment_WithGridMin,
    AxisAlignment_WithGridMax,
    AxisAlignment_WithSpecificAngle,
    AxisAlignment_WithGridAreaTop,
    AxisAlignment_WithGridAreaBottom,
    AxisAlignment_WithGridAreaLeft,
    AxisAlignment_WithGridAreaRight,
    END_AxisAlignment_e,
    AxisAlignment_Invalid = BadEnumValue
} AxisAlignment_e;

typedef enum
{
    FunctionDependency_XIndependent,
    FunctionDependency_YIndependent,
    END_FunctionDependency_e,
    FunctionDependency_Invalid = BadEnumValue,
    FunctionDependency_ThetaIndependent = FunctionDependency_XIndependent,
    FunctionDependency_RIndependent = FunctionDependency_YIndependent
} FunctionDependency_e;

typedef enum
{
    LegendShow_Yes,
    LegendShow_No,
    LegendShow_Auto,
    END_LegendShow_e,
    LegendShow_Invalid = BadEnumValue
} LegendShow_e;

typedef enum
{
    LineMapSort_None,
    LineMapSort_IndependentVar,
    LineMapSort_DependentVar,
    LineMapSort_SpecificVar,
    END_LineMapSort_e,
    LineMapSort_Invalid = BadEnumValue
} LineMapSort_e;

typedef enum
{
    ContLegendLabelLocation_ContourLevels,
    ContLegendLabelLocation_Increment,
    ContLegendLabelLocation_ColorMapDivisions,
    END_ContLegendLabelLocation_e,
    ContLegendLabelLocation_Invalid = BadEnumValue
} ContLegendLabelLocation_e;

typedef enum
{
    ThetaMode_Degrees,
    ThetaMode_Radians,
    ThetaMode_Arbitrary,
    END_ThetaMode_e,
    ThetaMode_Invalid = BadEnumValue
} ThetaMode_e;

typedef enum
{
    Transform_PolarToRect,
    Transform_SphericalToRect,
    Transform_RectToPolar,
    Transform_RectToSpherical,
    END_Transform_e,
    Transform_Invalid = BadEnumValue
} Transform_e;

typedef enum
{
    WindowFunction_Rectangular,
    WindowFunction_Triangular,
    WindowFunction_Hann,
    WindowFunction_Hamming,
    END_WindowFunction_e,
    WindowFunction_Invalid = BadEnumValue
} WindowFunction_e;

typedef enum
{
    LaunchDialogMode_ModalSync,
    LaunchDialogMode_Modeless,
    LaunchDialogMode_ModalAsync,
    END_LaunchDialogMode_e,
    LaunchDialogMode_Invalid = BadEnumValue
} LaunchDialogMode_e;


typedef enum
{
    SelectFileOption_ReadSingleFile,
    SelectFileOption_ReadMultiFile,
    SelectFileOption_AllowMultiFileRead,
    SelectFileOption_WriteFile,
    SelectFileOption_SelectDirectory,
    END_SelectFileOption_e,
    SelectFileOption_Invalid = BadEnumValue
} SelectFileOption_e;

typedef enum
{
    BinaryFileVersion_Tecplot2006,
    BinaryFileVersion_Tecplot2008,
    BinaryFileVersion_Tecplot2009,
    BinaryFileVersion_Current,
    END_BinaryFileVersion_e,
    BinaryFileVersion_Invalid = BadEnumValue
} BinaryFileVersion_e;

/*   CURRENTLY NOT USED .... */
typedef enum
{
    ViewActionDrawMode_NoDraw,
    ViewActionDrawMode_DrawTrace,
    ViewActionDrawMode_DrawFull,
    END_ViewActionDrawMode_e,
    ViewActionDrawMode_Invalid = BadEnumValue
} ViewActionDrawMode_e;

typedef enum
{
    PageAction_Create,
    PageAction_Delete,
    PageAction_Clear,
    PageAction_SetCurrentToNext,
    PageAction_SetCurrentToPrev,
    PageAction_SetCurrentByName,
    PageAction_SetCurrentByUniqueID,
    END_PageAction_e,
    PageAction_Invalid = BadEnumValue
} PageAction_e;

typedef enum
{
    FrameAction_PushTop,
    FrameAction_PopByNumber,
    FrameAction_PopAtPosition,
    FrameAction_DeleteActive,
    FrameAction_FitAllToPaper,
    FrameAction_PushByName,
    FrameAction_PopByName,
    FrameAction_PushByNumber,
    FrameAction_ActivateTop,
    FrameAction_ActivateNext,
    FrameAction_ActivatePrevious,
    FrameAction_ActivateAtPosition,
    FrameAction_ActivateByName,
    FrameAction_ActivateByNumber,
    FrameAction_MoveToTopActive,
    FrameAction_MoveToTopByName,
    FrameAction_MoveToTopByNumber,
    FrameAction_MoveToBottomActive,
    FrameAction_MoveToBottomByName,
    FrameAction_MoveToBottomByNumber,
    END_FrameAction_e,
    FrameAction_Invalid = BadEnumValue,
    FrameAction_Pop = FrameAction_PopByNumber,
    FrameAction_Push = FrameAction_PushByNumber,
    FrameAction_DeleteTop = FrameAction_DeleteActive
} FrameAction_e;

typedef enum
{
    DoubleBufferAction_On,
    DoubleBufferAction_Off,
    DoubleBufferAction_Swap,
    END_DoubleBufferAction_e,
    DoubleBufferAction_Invalid = BadEnumValue
} DoubleBufferAction_e;

/*
 * PickAction_CheckToAdd had the side effects of popping a frame that was selected
 * only if not collecting.  Pick_AddAtPosition avoids this.
 */
typedef enum
{
    PickAction_CheckToAdd, /* deprecated: use Pick_AddAtPosition */
    PickAction_AddAll,
    PickAction_AddAllInRegion,
    PickAction_Edit,
    PickAction_Cut,
    PickAction_Copy,
    PickAction_Clear,
    PickAction_Paste,
    PickAction_PasteAtPosition,
    PickAction_Shift,
    PickAction_Magnify,
    PickAction_Push,
    PickAction_Pop,
    PickAction_SetMouseMode,
    PickAction_DeselectAll,
    PickAction_AddZones,
    PickAction_AddXYMaps, /* deprecated: use PickAction_AddLineMaps */
    PickAction_AddLineMaps,
    PickAction_AddAtPosition,
    END_PickAction_e,
    PickAction_Invalid = BadEnumValue
} PickAction_e;


typedef enum
{
    ContourLevelAction_Add,
    ContourLevelAction_New,
    ContourLevelAction_DeleteRange,
    ContourLevelAction_Reset,
    ContourLevelAction_ResetToNice,
    ContourLevelAction_DeleteNearest,
    END_ContourLevelAction_e,
    ContourLevelAction_Invalid = BadEnumValue
} ContourLevelAction_e;

typedef enum
{
    ContourLabelAction_Add,
    ContourLabelAction_DeleteAll,
    END_ContourLabelAction_e,
    ContourLabelAction_Invalid = BadEnumValue
} ContourLabelAction_e;

typedef enum
{
    StreamtraceAction_Add,
    StreamtraceAction_DeleteAll,
    StreamtraceAction_DeleteRange,
    StreamtraceAction_SetTerminationLine,
    StreamtraceAction_ResetDeltaTime,
    END_StreamtraceAction_e,
    StreamtraceAction_Invalid = BadEnumValue
} StreamtraceAction_e;

typedef enum
{
    ColorMapControlAction_RedistributeControlPoints,
    ColorMapControlAction_CopyCannedColorMap,
    ColorMapControlAction_ResetToFactoryDefaults,
    END_ColorMapControlAction_e,
    ColorMapControlAction_Invalid = BadEnumValue
} ColorMapControlAction_e;

typedef enum
{
    ColorMapDistribution_Continuous,
    ColorMapDistribution_Banded,
    END_ColorMapDistribution_e,
    ColorMapDistribution_Invalid = BadEnumValue
} ColorMapDistribution_e;

typedef enum
{
    RGBMode_SpecifyRGB,
    RGBMode_SpecifyRG,
    RGBMode_SpecifyRB,
    RGBMode_SpecifyGB,
    END_RGBMode_e,
    RGBMode_Invalid = BadEnumValue
} RGBMode_e;

typedef enum
{
    TecUtilErr_None,
    TecUtilErr_Undetermined,
    END_TecUtilErr_e,
    TecUtilErr_Invalid = BadEnumValue
} TecUtilErr_e;

/* BEGINREMOVEFROMADDON */
/* deprecated type from alpha/beta v10 */
typedef enum
{
    AxisShape_Ray,
    AxisShape_LineTwoDirections,
    AxisShape_LShape,
    AxisShape_CrossOrBox,
    END_AxisShape_e,
    AxisShape_Invalid = BadEnumValue
} AxisShape_e;

/* licensing enums : keep hidden */
typedef enum
{
    RunMode_Demo,
    RunMode_Eval,
    RunMode_Full,
    /**/
    END_RunMode_e,
    /**/
    RunMode_Invalid = BadEnumValue
} RunMode_e;

/* ENDREMOVEFROMADDON */

typedef enum /* Custom exporter error message */
{
    ExportCustReturnCode_Ok,
    ExportCustReturnCode_Failed,
    ExportCustReturnCode_TecplotLocked,
    ExportCustReturnCode_ExporterNotLoaded,
    ExportCustReturnCode_ExportCallbackFailed,
    ExportCustReturnCode_NotAnImageExporter,
    ExportCustReturnCode_NotAFieldDataExporter,
    END_ExportCustReturnCode_e,
    ExportCustReturnCode_Invalid = BadEnumValue
} ExportCustReturnCode_e;

/**
 * COB/Zone types.
 */
typedef enum
{
    CZType_FieldDataZone,
    CZType_FEBoundaryCOB,
    CZType_IsoSurfaceCOB,
    CZType_SliceCOB,
    CZType_StreamtraceCOB,
    CZType_StreamtraceMarkerCOB,
    CZType_StreamtraceArrowheadCOB,
    END_CZType_e,
    CZType_Invalid = BadEnumValue
} CZType_e;

/**
 */
typedef enum
{
    FaceNeighborMode_LocalOneToOne,
    FaceNeighborMode_LocalOneToMany,
    FaceNeighborMode_GlobalOneToOne,
    FaceNeighborMode_GlobalOneToMany,
    END_FaceNeighborMode_e,
    FaceNeighborMode_Invalid = BadEnumValue
} FaceNeighborMode_e;


/**
 * Page render destinations.
 */
typedef enum
{
    PageRenderDest_None,
    PageRenderDest_OnScreen,
    PageRenderDest_OffScreen,
    END_PageRenderDest_e,
    PageRenderDest_Invalid = BadEnumValue
} PageRenderDest_e;

/* BEGINREMOVEFROMADDON */
/*
 * Destination for all internal rendering (VDI/Gr) functions. For external
 * linkage we translate RenderDest_WorkArea to PageRenderDest_OnScreen,
 * RenderDest_OffscreenBitmap to PageRenderDest_OffScreen and
 * RenderDest_Invalid to PageRenderDest_None.
 */
typedef enum
{
    RenderDest_WorkArea, /* Do not move from start of screen entries */
    RenderDest_ExampleText,
    RenderDest_ExampleLightSourcePosition,
    RenderDest_ExampleColorMap,
    RenderDest_ExampleBasicColor, /* Do not move from end of screen entries */
    RenderDest_OffscreenBitmap,
    RenderDest_Hardcopy,
    END_RenderDest_e,
    RenderDest_Invalid = BadEnumValue,
    /*
     * These next two are optimizations to make the
     * RDT_IsScreen() macro as efficient as possible.
     */
    RenderDest_FirstScreenEntry = RenderDest_WorkArea,
    RenderDest_LastScreenEntry = RenderDest_ExampleBasicColor
} RenderDest_e;
/* ENDREMOVEFROMADDON */

typedef enum
{
    Stipple_All,
    Stipple_Critical,
    Stipple_None,
    END_Stipple_e,
    Stipple_Invalid = BadEnumValue
} Stipple_e;

typedef enum
{
    DataFileType_Full,
    DataFileType_Grid,
    DataFileType_Solution,
    END_DataFileType_e,
    DataFileType_Invalid = BadEnumValue
} DataFileType_e;

typedef enum
{
    ConditionAwakeReason_Signaled,
    ConditionAwakeReason_TimedOut,
    END_ConditionAwakeReason_e,
    ConditionAwakeReason_Invalid = BadEnumValue
} ConditionAwakeReason_e;

typedef enum
{
    ProbeStatus_Normal,
    ProbeStatus_Terminated,
    ProbeStatus_Exited,
    END_ProbeStatus_e,
    ProbeStatus_Invalid = BadEnumValue
} ProbeStatus_e;

typedef enum
{
    FrameSizePosUnits_Paper,
    FrameSizePosUnits_Workspace,
    END_FrameSizePosUnits_e,
    FrameSizePosUnits_Invalid = BadEnumValue
} FrameSizePosUnits_e;

typedef enum
{
    Gridline_Major,
    Gridline_Minor,
    Gridline_Marker,
    END_Gridline_e,
    Gridline_Invalid = BadEnumValue
} Gridline_e;

/* Used by MarkerGridlineDetail */
typedef enum
{
    PositionMarkerBy_SolutionTime,
    PositionMarkerBy_Constant,
    END_PositionMarkerBy_e,
    PositionMarkerBy_Invalid = BadEnumValue
} PositionMarkerBy_e;


/****************************************************************
 *                                                              *
 *                     STRUCTURE TYPEDEFS                       *
 *                                                              *
 ****************************************************************/

/*
 * These are defined to work with pthreads, more work for WINAPI needed
 */
typedef struct _Mutex_a* Mutex_pa;

typedef struct _SpinLock_a* SpinLock_pa; /* NOTE: some platforms may not have spin locks; see use of HAVE_SPINLOCKS define */

typedef void*(STDCALL *ThreadFunction_pf)(ArbParam_t ThreadData);

typedef struct _Condition_a* Condition_pa;

typedef struct _JobControl_s* JobControl_pa;

typedef void (STDCALL *ThreadPoolJob_pf)(ArbParam_t JobData);

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined USE_OOSTYLE
#endif
#endif /* TECPLOTKERNEL */
/* ENDREMOVEFROMADDON */

typedef struct _StringList_s *StringList_pa;
typedef struct _Menu_s       *Menu_pa;
/* BEGINREMOVEFROMADDON */
typedef struct _ArrayList_s  *ArrayList_pa;
/* ENDREMOVEFROMADDON */
typedef struct _LineSegmentProbeResult_s *LineSegProbeResult_pa;

typedef enum
{
    ImageResizeFilter_Texture,
    ImageResizeFilter_Box,
    ImageResizeFilter_Lanczos2,
    ImageResizeFilter_Lanczos3,
    ImageResizeFilter_Triangle,
    ImageResizeFilter_Bell,
    ImageResizeFilter_BSpline,
    ImageResizeFilter_Cubic,
    ImageResizeFilter_Mitchell,
    ImageResizeFilter_Gaussian,
    END_ImageResizeFilter_e,
    ImageResizeFilter_Invalid = BadEnumValue
} ImageResizeFilter_e;

typedef enum
{
    VarStatus_Passive,
    VarStatus_Custom,
    VarStatus_Map,
    VarStatus_Heap,
    VarStatus_NotLoaded,
    END_VarStatus_e,
    VarStatus_Invalid = BadEnumValue
} VarStatus_e;



/* BEGINREMOVEFROMADDON */

/* here until GR and GRHW layer can be rearranged. */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#   if !defined NO_ASSERTS
#   endif
#endif /* TECPLOTKERNEL */

/* ENDREMOVEFROMADDON */

typedef struct _Set_a *Set_pa;

typedef struct
{
    double X;
    double Y;
    double Z;
} XYZ_s;

/* BEGINREMOVEFROMADDON */

typedef struct  
{
    double Psi;
    double Theta;
    double Alpha;
} PTA_s;

namespace tecplot
{
    class Typeface;
}

typedef struct _Generic3Var_s
{
    double V1;
    double V2;
    double V3;
} Generic3Var_s;

typedef struct _ThetaR_s
{
    double Theta;
    double R;
} ThetaR_s;

/*
 * This union is designed to allow different plot types
 * to access the same values by different names.  In
 * C++ we could use member access functions, or we
 * could have used macros, but instead we use this
 * union.  NOTE: This only works if all the structures
 * have the same alignment.
 */
typedef union _AnchorPos_u
{
    Generic3Var_s Generic;
    XYZ_s         XYZ;
    ThetaR_s      ThetaR;
} AnchorPos_u;

typedef struct _DataFileInfo_s
{
    char          *PrimaryFName;
    char          *TempBinaryFName;
    DataFileType_e FileType;
    FileOffset_t   DataFileOffset;
    StringList_pa  VarName;
    EntIndex_t     NumZones;
    EntIndex_t     NumVars;
    double         SolutionFileTime;
    struct _DataFileInfo_s *NextFile;
} DataFileInfo_s;

typedef struct _StylesheetIOFlags_s
{
    Boolean_t IncludePlotStyle;
    Boolean_t IncludeFieldAndMapStyle;     /* Only used for undo */
    Boolean_t IncludeUniqueIDs;            /* Only used for undo */
    Boolean_t IncludeText;
    Boolean_t IncludeGeom;
    Boolean_t IncludeGeomImageData;
    Boolean_t IncludeAuxData;
    Boolean_t IncludeStreamPositions;
    Boolean_t IncludeContourLevels;
    Boolean_t IncludeFactoryDefaults;      /* Only used when writing */
    Boolean_t CompressStyleCommands;       /* Only used when writing */
    Boolean_t MergeStyle;                  /* Only used when reading */
    Boolean_t IncludeFrameSizeAndPosition; /* Only used when reading */
    Boolean_t UseRelativePaths;
} StylesheetIOFlags_s;


/**
 */
typedef struct
{
    Boolean_t Show; /* power switch */
    Boolean_t ShowMesh;
    Boolean_t ShowContour;
    Boolean_t ShowShade;
    Boolean_t UseLightingEffect;
    Boolean_t UseTranslucency;
} IsoSurfaceLayers_s;

/**
 */
typedef struct
{
    Boolean_t Show; /* power switch */
    Boolean_t ShowMesh;
    Boolean_t ShowContour;
    Boolean_t ShowVector;
    Boolean_t ShowShade;
    Boolean_t ShowEdge;
    Boolean_t UseLightingEffect;
    Boolean_t UseTranslucency;
} SliceLayers_s;

/**
 */
typedef struct
{
    Boolean_t Show; /* power switch */
    Boolean_t ShowPaths;
    Boolean_t ShowDashes;
    Boolean_t ShowArrowheads;
    Boolean_t ShowMesh;
    Boolean_t ShowContour;
    Boolean_t ShowShade;
    Boolean_t ShowMarkers;
    Boolean_t UseLightingEffect;
    Boolean_t UseTranslucency;
} StreamtraceLayers_s;

/**
 */
typedef struct
{
#if 0 /* in the future we may add a main power switch */
    Boolean_t       Show; /* power switch */
#endif
    TwoDDrawOrder_e TwoDDrawOrder;
    Boolean_t       ShowMesh;
    Boolean_t       ShowContour;
    Boolean_t       ShowVector;
    Boolean_t       ShowScatter;
    Boolean_t       ShowShade;
    Boolean_t       ShowEdge;
    Boolean_t       UseLightingEffect;
    Boolean_t       UseTranslucency;
} FieldLayers_s;

/**
 * General purpose field layers structure used for low level drawing code only.
 * SetupXxxx is responsible for populating this general field layers structure
 * from the specific layer structures above for CZInfo.
 */
typedef struct
{
    Boolean_t ShowMesh;
    Boolean_t ShowContour;
    Boolean_t ShowVector;
    Boolean_t ShowScatter;
    Boolean_t ShowShade;
    Boolean_t ShowEdge;
    Boolean_t UseLightingEffect;
    Boolean_t UseTranslucency;
} CZFieldLayers_s;

/**
 */
typedef struct _LinePlotLayers_s
{
#if 0 /* in the future we may add a main power switch */
    Boolean_t       Show; /* power switch */
#endif
    Boolean_t ShowLines;
    Boolean_t ShowSymbols;
    Boolean_t ShowBarCharts;
    Boolean_t ShowErrorBars;
} LinePlotLayers_s;


typedef union _InterfaceAdjust_u
{
    double    ScaleFact;
    LgIndex_t Shift;
} InterfaceAdjust_u;

typedef Boolean_t (*SuffixModifier_pf)(TP_IN_OUT double* Value,
                                       const char*       Suffix);

typedef struct _InputSpecs_s
{
    Input_e           Type;
    double            Min;
    double            Max;
    InterfaceAdjust_u InterfaceAdjust;
    SuffixModifier_pf SuffixModifier;
} InputSpec_s;


typedef struct _RGB_s
{
    ColorIndex_t R;
    ColorIndex_t G;
    ColorIndex_t B;
} RGB_s;


typedef struct _ControlPoint_s
{
    double ColorMapFraction;
    RGB_s  LeadRGB;
    RGB_s  TrailRGB;
} ControlPoint_s;


typedef struct _ColorMapBand_s
{
    short          NumControlPoints;
    ControlPoint_s ControlPoint[MaxColorMapControlPoints];
} ColorMapBand_s;


typedef struct _EventAction_s
{
    int       I;
    int       J;
    int       LastI;
    int       LastJ;
    int       BaseI;
    int       BaseJ;
    int       ButtonOrKey;
    Event_e   Event;
    Boolean_t IsShifted;
    Boolean_t IsAlted;
    Boolean_t IsControlled;
    Boolean_t WasShiftedOnButtonPress;
    Boolean_t WasAltedOnButtonPress;
    Boolean_t WasControlledOnButtonPress;
} EventAction_s;

typedef struct _MacroCmd_s
{
    LString_t           MacroLine;
    struct _MacroCmd_s *NextCmd;
} MacroCmd_s;


typedef struct _IntegerRect_s
{
    LgIndex_t X1;
    LgIndex_t Y1;
    LgIndex_t X2;
    LgIndex_t Y2;
} IntegerRect_s;


typedef struct _Rect_s
{
    double X1;
    double Y1;
    double X2;
    double Y2;
} Rect_s;

typedef struct _XY_s
{
    double X;
    double Y;
} XY_s;

typedef struct _IJKSkip_s
{
    LgIndex_t  I;
    LgIndex_t  J;
    LgIndex_t  K;
} IJKSkip_s;



/*
 *
 *  NOTE ON RANGES (Ent and Index)
 *
 *  Min, Max and Skip all use the following assignment logic:
 *
 *              0 = First element
 *             -1 = mxindex value,  (X[mxindex-1] in c)
 *             -n = mxindex-n+1  value (X[mxindex+n] in c)
 *              n = n+1 value (X[n] in c)
 *
 */

/*
 *  2/28/95:  NOTE:  EntRange_s is no longer used but may be
 *                   needed later.
 */

typedef struct _EntRange_s
{
    EntIndex_t Min;
    EntIndex_t Max;
    EntIndex_t Skip;
} EntRange_s;


typedef struct _IndexRange_s
{
    LgIndex_t Min;
    LgIndex_t Max;
    LgIndex_t Skip;
} IndexRange_s;


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if defined (THREED)
#endif
#endif /* TECPLOTKERNEL */

typedef struct _TextShape_s
{
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#else
    Font_e                    Font;
#endif
    double                    Height;
    Units_e                   SizeUnits;
} TextShape_s;

#define AsciiShapeFontIsGreek(S)       (((S)->useBaseFont == FALSE) && (FontLibrary::instance().font((S)->typefaceOverride) == Font_Greek))
#define AsciiShapeFontIsMath(S)        (((S)->useBaseFont == FALSE) && (FontLibrary::instance().font((S)->typefaceOverride) == Font_Math))
#define AsciiShapeFontIsUserDefined(S) (((S)->useBaseFont == FALSE) && (FontLibrary::instance().font((S)->typefaceOverride) == Font_UserDefined))


typedef struct
{
    Boolean_t                useBaseFont;
    tecplot::Typeface const* typefaceOverride;
    SymbolChar_t             Char;
} AsciiShape_s;

typedef struct _SymbolShape_s
{
    GeomShape_e  GeomShape;
    Boolean_t    IsAscii;
    AsciiShape_s AsciiShape;
} SymbolShape_s;

#ifdef NOT_USED
struct _AddOnList_a
{
    /* added temporarily so Windows makelibtec works */
    int dummy;
};
#endif

/* ENDREMOVEFROMADDON */

typedef struct _AddOnList_a *AddOn_pa;

typedef struct _NodeMap_a *NodeMap_pa;

/* BEGINREMOVEFROMADDON */
typedef struct _StylePointState_a   *StylePointState_pa;
typedef struct _DataElementState_a  *DataElementState_pa;
typedef struct _StyleElementState_a *StyleElementState_pa;
typedef struct _NormalCache_a       *NormalCache_pa;
/* ENDREMOVEFROMADDON */


#define INVALID_INDEX (-1)

/* used to indicate that no neighboring element or zone exists */
#define NO_NEIGHBORING_ELEMENT (-1)
#define NO_NEIGHBORING_ZONE    (-1)

typedef struct _FaceNeighbor_a *FaceNeighbor_pa;

/**
 */
typedef struct _FaceMap_a *FaceMap_pa;

/**
 */
typedef struct _ElemToFaceMap_a *ElemToFaceMap_pa;

/**
 */
typedef struct _NodeToElemMap_a *NodeToElemMap_pa;

/* BEGINREMOVEFROMADDON */

/*
 * Enumerates the face neighbor array members to make indexed members
 * identifiable.
 */
typedef enum
{
    FaceNeighborMemberArray_CellFaceNbrs,
    FaceNeighborMemberArray_BndryConnectNbrsCompObscure,
    FaceNeighborMemberArray_BndryConnectFaceToCellsMap,
    FaceNeighborMemberArray_BndryConnectIsPerfectNbr,
    FaceNeighborMemberArray_BndryConnectCellList,
    FaceNeighborMemberArray_BndryConnectZoneList,
    END_FaceNeighborMemberArray_e,
    FaceNeighborMemberArray_Invalid = BadEnumValue
} FaceNeighborMemberArray_e;

int const FaceNeighborNumMemberArrays = (int)END_FaceNeighborMemberArray_e;

/*
 * Enumerates the face map's array members to make indexed members
 * identifiable.
 */
typedef enum
{
    FaceMapMemberArray_FaceNodeOffsets,
    FaceMapMemberArray_FaceNodes,
    FaceMapMemberArray_FaceLeftElems,
    FaceMapMemberArray_FaceRightElems,
    FaceMapMemberArray_FaceBndryItemOffsets,
    FaceMapMemberArray_FaceBndryItemElems,
    FaceMapMemberArray_FaceBndryItemElemZones,
    END_FaceMapMemberArray_e,
    FaceMapMemberArray_Invalid = BadEnumValue
} FaceMapMemberArray_e;

const int FaceMapNumMemberArrays = (int)END_FaceMapMemberArray_e;

/*
 * Enumerates the element to face map's array members to make indexed members
 * identifiable.
 */
typedef enum
{
    ElemToFaceMapMemberArray_ElemFaceOffsets,
    ElemToFaceMapMemberArray_ElemFaces,
    END_ElemToFaceMapMemberArray_e,
    ElemToFaceMapMemberArray_Invalid = BadEnumValue
} ElemToFaceMapMemberArray_e;

const int ElemToFaceMapNumMemberArrays = (int)END_ElemToFaceMapMemberArray_e;

/*
 * Enumerates the element map's array members to make indexed members
 * identifiable.
 */
typedef enum
{
    NodeToElemMapMemberArray_NodeElemOffsets,
    NodeToElemMapMemberArray_NodeElems,
    END_NodeToElemMapMemberArray_e,
    NodeToElemMapMemberArray_Invalid = BadEnumValue
} NodeToElemMapMemberArray_e;

const int NodeToElemMapNumMemberArrays = (int)END_NodeToElemMapMemberArray_e;

/* ENDREMOVEFROMADDON */


typedef struct _FieldData_a *FieldData_pa;

/**
 */
typedef struct _AuxData_s  *AuxData_pa;


/**
 * Enumerates the data value structure of a variable in a data file.
 * For all but ordered cell centered data the classic, classic padded and
 * classic plus formats are identical. All values are laid out contiguously
 * in the file. The number of values written depends upon the value location:
 *
 *   - FE nodal:\n
 *     The number of values equals the number of data points.
 *   - FE cell centered:\n
 *     The number of values equals the number of elements.
 *   - Ordered nodal:\n
 *     The number of values equals the number of data points.
 *   - Ordered cell centered:\n
 *     There are three formats:
 *     -# Classic (binary version < 103):\n
 *          Classic is a compressed format of ordered cell centered data in
 *          that it does not include ghost cells. The cell index of each cell
 *          does not correspond to the lowest corner point index of each cell
 *          as it does internally in Tecplot.\n
 *          The number of values in the data file is calculated as follows:
 *          @code
 *            NumValues = MAX(IMax-1,1) * MAX(JMax-1,1) * MAX(KMax-1,1);
 *          @endcode
 *          Where IMax, JMax, and KMax are the maximum point dimensions of the
 *          zone.
 *     -# Classic padded (binary version < 104):\n
 *          Classic padded is an intermediary format that was available only
 *          within Tecplot, Inc. The cell centered data includes the ghost cells
 *          and each cell index corresponds to the lowest corner point index of
 *          each cell.\n
 *          The number of values in the data file (including ghost cells) is
 *          calculated as follows:
 *          @code
 *            NumValues = IMax * JMax * KMax;
 *          @endcode
 *          Where IMax, JMax, and KMax are the maximum point dimensions of the
 *          zone. The contents of the ghost cells is undefined and should not
 *          be used.
 *     -# Classic plus (binary version >= 104):\n
 *          Classic plus is similar to classic padded except that it does not
 *          include the ghost cells of the slowest moving index greater than
 *          one.\n
 *          The number of values in the data file (including ghost cells) is
 *          calculated as follows:
 *          @code
 *            FinalIMax = IMax;
 *            FinalJMax = JMax;
 *            FinalKMax = KMax;
 *
 *            // decrement the max index of the slowest moving index greater than 1
 *            if (KMax > 1)
 *              FinalKMax--;
 *            else if (JMax > 1)
 *              FinalJMax--;
 *            else if (IMax > 1)
 *              FinalIMax--;
 *
 *            NumValues = FinalIMax * FinalJMax * FinalKMax;
 *          @endcode
 *          Where IMax, JMax, and KMax are the maximum point dimensions of the
 *          zone. The contents of the ghost cells is undefined and should not
 *          be used.
 */
typedef enum
{
    DataValueStructure_Classic,
    DataValueStructure_ClassicPadded,
    DataValueStructure_ClassicPlus,
    END_DataValueStructure_e,
    /* BEGINREMOVEFROMADDON */
    DataValueStructure_Latest = (END_DataValueStructure_e - 1),
    /* ENDREMOVEFROMADDON */
    DataValueStructure_Invalid = BadEnumValue
} DataValueStructure_e;

/**
 * Enumerates the data node structure of a node map in a data file. The classic
 * format uses 1 based nodes while the classic plus format uses zero based
 * node.
 */
typedef enum
{
    DataNodeStructure_Classic,     /* ones based node maps */
    DataNodeStructure_ClassicPlus, /* zero based node maps */
    END_DataNodeStructure_e,
    DataNodeStructure_Invalid = BadEnumValue
} DataNodeStructure_e;

/**
 * Enumerates the variable locking modes. The \ref VarLockMode_ValueChange mode
 * prevents modification of the values in a variable but permits deletion, and
 * the \ref VarLockMode_Delete mode prevents deletion of a varaible but permits
 * modification.
 */
typedef enum
{
    VarLockMode_ValueChange,
    VarLockMode_Delete,
    END_VarLockMode_e,
    VarLockMode_Invalid = BadEnumValue
} VarLockMode_e;

typedef enum
{
    FieldMapMode_UseStrandID,
    FieldMapMode_UseZoneSet,
    END_FieldMapMode_e,
    FieldMapMode_Invalid = BadEnumValue
} FieldMapMode_e;

typedef enum
{
    UnloadStrategy_Auto,
    UnloadStrategy_NeverUnload,
    UnloadStrategy_MinimizeMemoryUse,
    END_UnloadStrategy_e,
    UnloadStrategy_Invalid = BadEnumValue
} UnloadStrategy_e;

/* BEGINREMOVEFROMADDON */



typedef struct
{
    ColorIndex_t       PresetZoneColor;
    Boolean_t          IsInBlockFormat;
} ZoneLoadInfo_s;

/*
 * Note: For FE Data, NumPtsI = Number of data points.
 *                    NumPtsJ = Number of elements.
 *                    NumPtsK = Number of points per element.
 */

typedef struct _ZoneSpec_s
{
    UniqueID_t         UniqueID;
    ZoneName_t         Name;
    EntIndex_t         ParentZone;
    Strand_t           StrandID;
    double             SolutionTime;
    LgIndex_t          NumPtsI;  /* ...NumDataPts */
    LgIndex_t          NumPtsJ;  /* ...NumElements */
    LgIndex_t          NumPtsK;  /* ...NumPtsPerElem or NumFaces */
    LgIndex_t          ICellDim; /* ...currently not used */
    LgIndex_t          JCellDim; /* ...currently not used */
    LgIndex_t          KCellDim; /* ...currently not used */
    ZoneType_e         Type;
    ZoneLoadInfo_s     ZoneLoadInfo;
    AuxData_pa         AuxData;
    Boolean_t          BuildZoneOptInfo;

    /* classic data only */
    FaceNeighborMode_e FNMode;
    Boolean_t          FNAreCellFaceNbrsSupplied; /* ...meaning we don't need to update them */

    /* polytope data only */
    LgIndex_t          NumFaceNodes;
    LgIndex_t          NumFaceBndryFaces;
    LgIndex_t          NumFaceBndryItems;
} ZoneSpec_s;



#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */

typedef struct _GenericImage_a *GenericImage_pa;

typedef struct _TextBox_s
{
    TextBox_e        BoxType;       /* Used to be textbox */
    double           Margin;        /* Used to be textboxmargin */
    double           LineThickness; /* Used to be textboxmargin */
    ColorIndex_t     BColor;        /* Used to be textboxcolor */
    ColorIndex_t     FillBColor;    /* Used to be textboxfillcolor */
} TextBox_s;


typedef struct _Text_s
{
    UniqueID_t       UniqueID; /* Not used yet */
    AnchorPos_u      AnchorPos;
    CoordSys_e       PositionCoordSys;
    EntIndex_t       Zone;
    Boolean_t        AttachToZone; /* New */
    ColorIndex_t     BColor;       /* Used to be TextColor */
    TextShape_s      TextShape;
    TextBox_s        Box;          /* Box items used to be here*/
    double           Angle;        /* NOTE: short in v6, now in rad */
    TextAnchor_e     Anchor;       /* New */
    double           LineSpacing;  /* New */
    Scope_e          Scope;
    char            *MacroFunctionCommand;
    Clipping_e       Clipping;
    char            *Text;
    struct _Text_s  *NextText;
    struct _Text_s  *PrevText;
} Text_s;


typedef struct _GenericGeomData_s
{
    FieldData_pa  V1Base;
    FieldData_pa  V2Base;
    FieldData_pa  V3Base;
} GenericGeomData_s;

typedef struct _PolarGeomData_s
{
    FieldData_pa  ThetaBase;
    FieldData_pa  RBase;
} PolarGeomData_s;

typedef struct _CartesianGeomData_s
{
    FieldData_pa  XBase;
    FieldData_pa  YBase;
    FieldData_pa  ZBase;
} CartesianGeomData_s;

/*
 * This union is designed to allow different plottypes
 * to access the same values by different names.  In
 * C++ we could use member access functions, or we
 * could have used macros, but instead we use this
 * union.  NOTE: This only works if all the structures
 * have the same alignment.
 */
typedef union _GeomData_u
{
    GenericGeomData_s   Generic;
    CartesianGeomData_s XYZ;
    PolarGeomData_s     ThetaR;
} GeomData_u;

typedef struct _Geom_s
{
    UniqueID_t              UniqueID;
    GeomType_e              GeomType;
    CoordSys_e              PositionCoordSys;
    AnchorPos_u             AnchorPos;
    Boolean_t               AttachToZone;
    EntIndex_t              Zone;
    ColorIndex_t            BColor;
    Boolean_t               IsFilled;
    ColorIndex_t            FillBColor;
    LinePattern_e           LinePattern;
    double                  PatternLength;
    double                  LineThickness;
    Scope_e                 Scope;
    DrawOrder_e             DrawOrder;
    Clipping_e              Clipping;
    FieldDataType_e         DataType;
    char                   *MacroFunctionCommand;
    ArrowheadStyle_e        ArrowheadStyle;
    ArrowheadAttachment_e   ArrowheadAttachment;
    double                  ArrowheadSize;
    double                  ArrowheadAngle;
    SmInteger_t             NumEllipsePts;
    char                   *ImageFileName;
    LgIndex_t               ImageNumber; /* used only to locate images within .lpk files */
    Boolean_t               MaintainAspectRatio;
    double                  PixelAspectRatio; /* VerticalPixelsPerHorizontalPixel */
    SmInteger_t             NumSegments;
    SegPtsArray_t           NumSegPts;
    GeomData_u              GeomData;
    ImageResizeFilter_e     ImageResizeFilter;
    /* Internal Scratch */
    GenericImage_pa         _ImageData;
    struct _Geom_s         *_NextGeom;
    struct _Geom_s         *_PrevGeom;
} Geom_s;


typedef struct _Text_s  *Text_pa;
typedef struct _Geom_s  *Geom_pa;


#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
    #if defined USE_VBOs
    #endif
#if defined USE_OOSTYLE
#endif
#if defined USE_OOSTYLE
#endif
#endif /* TECPLOTKERNEL */

/* ENDREMOVEFROMADDON */
/* - NO DOXYGEN COMMENT GENERATION -
 * Page creation callback is responsible for creating a RenderHandler for the page and
 * calling @ref TecEngPageCreateNew(ArbParam_t RenderHandle)
 *
 * The RenderHandler type can be anything, for example,a pointer to a class instance that will
 * be responsible for handling requests from the engine to perform operations on
 * a page.
 *
 * @param PageConstructionHints a string list of construction hints that can be used for deciding
 * how the page should be displayed in an application's UI. The construction hints could have been
 * restored from a saved layout file or passed to @ref TecUtilPageCreateNew function.
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @return TRUE if page create request was handled and TecEngPageCreateNew() returned TRUE.
 *
 * @sa TecEngPageCreateRegisterCallback, TecEngPageCreateNew
 *
 * @since
 *   11.0-5-014
 */
typedef Boolean_t (STDCALL *PageCreateCallback_pf)(StringList_pa PageConstructionHints,
                                                   ArbParam_t    RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Page destruction callback responsible for destroying a page.
 *
 * @param PageClientData
 *   Data associated with a page that was returned from the PageCreateCallback_pf
 *   callback function.   You will get a different value for each page.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @sa TecEngPageDestroyRegisterCallback, PageCreateCallback_pf
 *
 * @since
 *   11.0-5-014
 */
typedef void (STDCALL *PageDestroyCallback_pf)(ArbParam_t PageClientData,
                                               ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for informing the parent application of a new current page.
 * Note that this could be done via a state change monitor but a more secure method
 * is needed as state changes may be shut down from time to time.
 *
 * @param PageClientData
 *   Data associated with a page that was returned from the PageCreateCallback_pf
 *   callback function.   You will get a different value for each page.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @since
 *   11.0-5-017
 */
typedef void (STDCALL *PageNewCurrentCallback_pf)(ArbParam_t PageClientData,
                                                  ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for creation of an offscreen image.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @param ImageHandle handle to a newly created image. This is an output parameter.
 *
 * @return TRUE if an offscreen image was created successfully.
 *
 * @since
 *   11.2-0-054
 */
typedef Boolean_t (STDCALL *OffscreenImageCreateCallback_pf)(ScreenDim_t        Width,
                                                             ScreenDim_t        Height,
                                                             ArbParam_t         RegistrationClientData,
                                                             TP_OUT ArbParam_t* ImageHandle);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for destruction of an offscreen image.
 *
 * @param ImageHandle handle to an offscreen image to be destroyed.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @since
 *   11.2-0-054
 */
typedef void (STDCALL *OffscreenImageDestroyCallback_pf)(ArbParam_t ImageHandle,
                                                         ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for returning RGB values for a row.
 *
 * @param ImageHandle
 *     Handle to an off-screen image from which RGB values to be retrieved.
 *
 * @param Row
 *     Row for which RGB values to be retrieved.
 *
 * @param RedArray
 *     Array to receive the red byte values for the specified Row. The number of values in
 *     the array must equal the width of the image. The array address is maintained by the
 *     Tecplot Engine until the image is destroyed however it is reused for each invocation
 *     of this callback.
 *
 * @param GreenArray
 *     Array to receive the green byte values for the specified Row. The number of values in
 *     the array must equal the width of the image. The array address is maintained by the
 *     Tecplot Engine until the image is destroyed however it is reused for each invocation
 *     of this callback.
 *
 * @param BlueArray
 *     Array to receive the blue byte values for the specified Row. The number of values in
 *     the array must equal the width of the image. The array address is maintained by the
 *     Tecplot Engine until the image is destroyed however it is reused for each invocation
 *     of this callback.
 *
 * @param RegistrationClientData
 *     Data associated with the registration of this function.   This will always return
 *     the value supplied in the original registration of this function.
 *
 * @return TRUE if successful, FALSE otherwise.
 *
 * @since
 *     11.2-0-054
 */
typedef Boolean_t (STDCALL *OffscreenImageGetRGBRowCallback_pf)(ArbParam_t           ImageHandle,
                                                                ScreenDim_t          Row,
                                                                ArbParam_t           RegistrationClientData,
                                                                TP_ARRAY_OUT Byte_t* RedArray,
                                                                TP_ARRAY_OUT Byte_t* GreenArray,
                                                                TP_ARRAY_OUT Byte_t* BlueArray);

#if defined MSWIN
/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for printing an image on the specified printer DC
 *
 * @param PrintDC a device context of a printer on which the printing should be performed.
 *
 * @param ImageHandle handle to an image to print.
 *
 * @param Palette specifies if an image should be printed as a color or monochrome image.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @return TRUE if the printing operation was successfull.
 *
 * @since
 *   11.2-0-463
 */
typedef Boolean_t (STDCALL *WinPrintImageCallback_pf)(HDC        PrintDC,
                                                      ArbParam_t ImageHandle,
                                                      Palette_e  Palette,
                                                      ArbParam_t RegistrationClientData);

#endif /* MSWIN */

#if defined MSWIN
/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for providing a printer context.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function.   This will always return
 *   the value supplied in the original registration of this function.
 *
 * @return HDC context of the destination printer.
 *
 * @since
 *   11.2-0-468
 */
typedef HDC(STDCALL *WinPrinterGetContextCallback_pf)(ArbParam_t RegistrationClientData);

#endif /* MSWIN */

/* - NO DOXYGEN COMMENT GENERATION -
 * Render destination callback responsible for switching the render destination
 * of the OpenGL drawing state when requested by the Tecplot engine.
 *
 * @since
 *   11.0-0-397
 *
 * @param PageRenderDest
 *   Enumeration of page render destination of interest.
 *
 * @param RenderDestClientData
 *   Data associated with a render destination, such as returned from the PageCreateCallback_pf or
 *   OffscreenImageCreate_pf callback functions.
 *
 * @param RegistrationClientData
 *   Data associated with the registration of this function. This will always return
 *   the value supplied in the original registration of this function.
 *
 * @return
 *   TRUE if render destination was set successfully. FALSE, otherwise.
 *
 * @sa TecEngRenderDestRegisterCallback
 */
typedef Boolean_t (STDCALL *RenderDestCallback_pf)(PageRenderDest_e PageRenderDest,
                                                   ArbParam_t       RenderDestClientData,
                                                   ArbParam_t       RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Render query callback responsible for informing Tecplot if the page
 * associated with the PageClientData should be rendered into.
 *
 * @since
 *   11.0-5-018
 *
 * @param PageClientData
 *   Data associated with a page that was returned from the
 *   PageCreateCallback_pf callback function.
 * @param RegistrationClientData
 *   Data associated with the registration of this function. This will always
 *   return the value supplied in the original registration of this function.
 *
 *
 * @return
 *   TRUE if Tecplot should render to the page identified by the
 *   PageClientData, FALSE otherwise.
 *
 * @sa TecEngRenderQueryRegisterCallback
 */
typedef Boolean_t (STDCALL *RenderQueryCallback_pf)(ArbParam_t PageClientData,
                                                    ArbParam_t RegistrationClientData);
/* - NO DOXYGEN COMMENT GENERATION -
 * Render destination size callback responsible for returning the size of the
 * specified render destination when requested by the Tecplot engine.
 *
 * @since
 *   11.0-0-397
 *
 * @param PageClientData
 *   Data associated with a page that was returned from the
 *   PageCreateCallback_pf callback function.
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 * @param Width
 *   Pointer who's contents should receive the width of the current render
 *   destination.
 * @param Height
 *   Pointer who's contents should receive the height of the current render
 *   destination.
 *
 * @sa TecEngRenderDestSizeRegisterCallback
 */
typedef void (STDCALL *RenderDestSizeCallback_pf)(ArbParam_t        PageClientData,
                                                  ArbParam_t        RegistrationClientData,
                                                  TP_OUT LgIndex_t* Width,
                                                  TP_OUT LgIndex_t* Height);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for swapping the front and back buffers for the current
 * OpenGL drawing state's render destination when requested by the Tecplot
 * engine.
 *
 * @since
 *   11.0-0-397
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecUtilpBuffersRegisterCallback
 */
typedef void (STDCALL *SwapBuffersCallback_pf)(ArbParam_t RegistrationClientData);


/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for querying of key states.
 *
 * @since
 *   11.0-0-399
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 * @param IsShiftKeyDown
 *   Boolean pointer. If non-NULL, set the boolean to TRUE if the Shift key is
 *   down or FALSE if it is up.
 * @param IsAltKeyDown
 *   Boolean pointer. If non-NULL, set the boolean to TRUE if the Alt key is
 *   down or FALSE if it is up.
 * @param IsCntrlKeyDown
 *   Boolean pointer. If non-NULL, set the boolean to TRUE if the Cntrl key is
 *   down or FALSE if it is up.
 *
 * @sa TecEngKeyStateRegisterCallback
 */
typedef void (STDCALL *KeyStateCallback_pf)(ArbParam_t        RegistrationClientData,
                                            TP_OUT Boolean_t* IsShiftKeyDown,
                                            TP_OUT Boolean_t* IsAltKeyDown,
                                            TP_OUT Boolean_t* IsCntrlKeyDown);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for querying of a mouse button state.
 *
 * @since
 *   11.0-0-424
 *
 * @param Button
 *   Mouse button number to query. Button numbers start at one.
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @return
 *   TRUE if the specified mouse button is down, FALSE otherwise.
 *
 * @sa TecEngMouseButtonStateRegisterCallback
 */
typedef Boolean_t (STDCALL *MouseButtonStateCallback_pf)(int        Button,
                                                         ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for setting wait cursor when requested by the kernel
 *
 * @since
 *   11.2-0-302
 *
 * @param Activate
 *   TRUE if the kernel is requesting that the wait cursor be activated.
 *   FALSE if the kernel is requesting that the wait cursor be deactivated.
 * @param RegistractionClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngWaitCursorStateRegisterCallback
 */
typedef void (STDCALL *WaitCursorStateCallback_pf)(Boolean_t  Activate,
                                                   ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for setting cursor style when requested by the kernel
 *
 * @since
 *   11.2-0-302
 *
 * @param CursorStyle
 *   The cursor style which the kernel is requesting.
 * @param RenderHandle
 *   Handle to page where new cursor shape is being set.
 * @param RegistractionClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngBaseCursorStyleRegisterCallback
 */
typedef void (STDCALL *BaseCursorStyleCallback_pf)(CursorStyle_e CursorStyle,
                                                   ArbParam_t    RenderHandle,
                                                   ArbParam_t    RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for processing events when the Tecplot engine is busy
 * peforming a requested operation. This callback will be called at regular
 * intervals to repair the interface and if required check for interrupts. Very
 * little work should be done by this function.
 *
 * @since
 *   11.0-0-415
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngProcessBusyEventsRegisterCallback, TecUtilInterrupt
 */
typedef void (STDCALL *ProcessBusyEventsCallback_pf)(ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for launching a dialog.
 *
 * @since
 *   11.0-0-415
 *
 * @param RegistrationClientData
 *   Client data that was registered with this launch dialog callback.
 *
 * @return
 *   TRUE if the dialog was launched, FALSE if it could not be launched
 *   programmatically.
 *
 * @sa TecUtilDialogLaunch, TecUtilDialogDrop
 */
typedef Boolean_t (STDCALL *DialogLaunchCallback_pf)(ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for dropping a dialog.
 *
 * @since
 *   11.0-0-407
 *
 * @param RegistrationClientData
 *   Client data that was registered with this drop dialog callback.
 *
 * @sa TecUtilDialogLaunch, TecUtilDialogDrop
 */
typedef void (STDCALL *DialogDropCallback_pf)(ArbParam_t RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for querying of the physical display's horizontal and
 * vertical dot pitch.
 *
 * @since
 *   11.0-0-407
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 * @param IDotsPerCm
 *   Pointer who's contents should receive the physical display's horizontal
 *   dot pitch in terms of the number of dots per centimeter.
 * @param JDotsPerCm
 *   Pointer who's contents should receive the physical display's vertical
 *   dot pitch in terms of the number of dots per centimeter.
 *
 * @sa TecEngDotPitchRegisterCallback
 */
typedef void (STDCALL *DotPitchCallback_pf)(ArbParam_t     RegistrationClientData,
                                            TP_OUT double* IDotsPerCm,
                                            TP_OUT double* JDotsPerCm);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for querying of the physical display's width and
 * height in pixels.
 *
 * @since
 *   11.2-0-471
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 * @param WidthInPixels
 *   Pointer who's contents should receive the physical display's width
 *   in pixels. NULL may be passed.
 * @param HeightInPixels
 *   Pointer who's contents should receive the physical display's height
 *   in pixels. NULL may be passed.
 *
 * @sa TecEngScreenSizeRegisterCallback
 */
typedef void (STDCALL *ScreenSizeCallback_pf)(ArbParam_t  RegistrationClientData,
                                              TP_OUT int* WidthInPixels,
                                              TP_OUT int* HeightInPixels);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for displaying a message box dialog and returning the
 * user's response.
 *
 * @since
 *   11.0-0-415
 *
 * @param MessageString
 *   Message string to display in the dialog.
 * @param MessageBoxType
 *   Type of message box to display.
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @return
 *   Result of user's response to the dialog.
 *
 * @sa TecEngDialogMessageBoxRegisterCallback
 */
typedef MessageBoxReply_e(STDCALL *DialogMessageBoxCallback_pf)(const char*      MessageString,
                                                                MessageBoxType_e MessageBoxType,
                                                                ArbParam_t       RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback responsible for displaying a status line
 *
 * @since
 *   11.2-0-085
 *
 * @param StatusString
 *   Message string to display in the dialog.
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngStatusLineRegisterCallback
 */
typedef void (STDCALL *StatusLineCallback_pf)(const char* StatusString,
                                              ArbParam_t  RegistrationClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback that will be called with the updated progress status.
 *
 * @since 11.2-0-098
 *
 *
 * @param ProgressStatus
 *   Percentage of the progress.
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngProgressMonitorRegisterCallback
 */
typedef void (STDCALL *ProgressMonitorCallback_pf)(int        ProgressStatus,
                                                   ArbParam_t RegistrationClientData);
/* - NO DOXYGEN COMMENT GENERATION -
 * Callback that will be called with Tecplot Engine is about to perform a lengthy operation.
 * The client that registers such the callback may present a user with a progress bar,
 * if the ShowProgressBar argument is TRUE, and a stop button that would interrupt the operation by
 * calling TecUtilInterrupt().
 *
 * @since 11.2-0-098
 *
 * @param ShowProgressBar
 *   Boolean indicating if the progress steps can be monitored for an operation. If TRUE, Tecplot Engine will be calling
 *   the registered ProgressMonitorCallback_pf function with the updated progress status.
 *
 * @param IsInterruptible
 *   Boolean indicating if the operation can be interrupted before completion.
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngProgressMonitorRegisterCallback
 */
typedef void (STDCALL *ProgressMonitorStartCallback_pf)(Boolean_t  ShowProgressBar,
                                                        Boolean_t  IsInterruptible,
                                                        ArbParam_t RegistrationClientData);
/* - NO DOXYGEN COMMENT GENERATION -
 * Callback tht will be called with Tecplot Engine has finished performing a lengthy operation.
 * At this point, client may hide progress bar that was shown during handling of ProgressMonitorStartCallback callback and
 * disable or hide the stop button.
 *
 * @since 11.2-0-098
 *
 * @param RegistrationClientData
 *   Client data that was registered with the callback.
 *
 * @sa TecEngProgressMonitorRegisterCallback
 */
typedef void (STDCALL *ProgressMonitorFinishCallback_pf)(ArbParam_t RegistrationClientData);

/*********************************************************
 * Add-on Timers
 *********************************************************/
/**
 * This is called when a registered timer fires.
 *
 * @par Limitation:
 *   Unix and Linux versions of Tecplot currently do not fire timer events when
 *   Tecplot is running in batch mode (with the -b flag). This behavior
 *   limitation is subject to change.
 *
 * @param ClientData
 *   Arbitrary client data.
 *
 * @return
 *   Return TRUE if the timer should be reinstated.   Return FALSE
 *   to stop subsequent callbacks.
 *
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyAddOnTimerCallback(
 *   &                     ClientDataPtr)
 *    POINTER (ClientDataPtr,DummyClientData)
 * </FortranSyntax>
 */
typedef Boolean_t (STDCALL *AddOnTimerCallback_pf)(ArbParam_t ClientData);

/* - NO DOXYGEN COMMENT GENERATION -
 * Callback that will be called when Tecplot Engine has requested an event timer to be created.
 *
 * @since 12.0.1.5642
 *
 * @param ClientData
 *   ClientData that should be sent in the callback.
 *
 * @param TimerCallback
 *   Callback to fire when the timer interval has expired.
 *
 * @param Interval
 *   The time (in milliseconds) after which the timer callback should be called.
 *
 * @param RegistrationClientData
 *   Client data that was registered via TecEngTimerRegisterCallback.
 *
 * @return
 *   Return TRUE if the timer was successfully created, FALSE if not.
 */
typedef Boolean_t (STDCALL *TimerCallback_pf)(AddOnTimerCallback_pf  TimerCallback,
                                              ArbParam_t             ClientData,
                                              UInt32_t               Interval,
                                              ArbParam_t             RegistrationClientData);

/**
 * This function is called when the user activates a menu item
 * added via TecUtilMenuInsertOption or TecUtilMenuInsertToggle.
 *
 * @param RegistrationClientData
 *   Arbitrary client data.
 */
typedef void (STDCALL *MenuActivateCallback_pf)(ArbParam_t RegistrationClientData);

/**
 * This function is called when the a menu is deleted.
 *
 * @param RegistrationClientData
 *   Arbitrary client data.
 */
typedef void (STDCALL *MenuDeleteCallback_pf)(ArbParam_t RegistrationClientData);

/**
 * This function is called to determine the sensitivity for a menu item (option,
 * toggle or submenu).
 *
 * @param RegistrationClientData
 *   Arbitrary client data.
 *
 * @return
 *   Return TRUE if the menu item should be sensitive to user input,
 *   or FALSE if it should be insensitive to user input (gray).
 */
typedef Boolean_t (STDCALL *MenuGetSensitivityCallback_pf)(ArbParam_t RegistrationClientData);

/**
 * This function is called to determine the checked state for a toggle menu item.
 *
 * @param RegistrationClientData
 *   Arbitrary client data.
 *
 * @return
 *   Return TRUE if the toggle should be checked,
 *   or FALSE if it should be unchecked.
 */
typedef Boolean_t (STDCALL *MenuGetToggleStateCallback_pf)(ArbParam_t RegistrationClientData);


/**
 * This function is called when the user performs a probe event.
 *
 * @param IsNearestPoint
 *   This is TRUE if the previous probe event was a nearest point probe.
 *   This is FALSE if it was an interpolated probe.
 *
 * <FortranSyntax>
 *   SUBROUTINE MyProbeDestinationCallback(
 *              IsNearestPoint)
 *   INTEGER*4 IsNearestPoint
 * </FortranSyntax>
 */
typedef void (STDCALL *ProbeDestination_pf)(Boolean_t IsNearestPoint);


/**
 * This function type called when a probe callback is installed via
 * TecUtilProbeInstallCallbackX.
 *
 * @param WasSuccessful
 *   This is TRUE if the previous probe event was successful.
 *   This is FALSE if it was the probe failed. Probe events may fail if the
 *   user probes in a region of the plot that contains no data.
 *
 * @param IsNearestPoint
 *   This is TRUE if the previous probe event was a nearest point probe.
 *   This is FALSE if it was an interpolated probe.
 *
 * @param ClientData
 *   Arbitrary client data.
 *
 */
typedef void (STDCALL *ProbeDestinationX_pf)(Boolean_t  WasSuccessful,
                                             Boolean_t  IsNearestPoint,
                                             ArbParam_t ClientData);


/**
 * DynamicMenu Functions are called upon a user selecting
 * a menu item added via TecUtilMenuAddOption.
 *
 * <FortranSyntax>
 *   SUBROUTINE MyDynamicMenuCallback()
 * </FortranSyntax>
 */
typedef void (STDCALL *DynamicMenuCallback_pf)(void);

/**
 * This callback signature is used to perform redraw events.
 *
 * @since
 *   11.0-0-363
 *
 * @param RedrawReason
 *   An enumerated value describing the reason for the re-draw event.
 * @param ClientData
 *   Client data that was registered with the callback.
 *
 * @return
 *   TRUE if successfull, FALSE otherwise.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION DrawEventCallback(
 *   &                     RedrawReason,
 *   &                     ClientDataPtr)
 *    INTEGER*4 RedrawReason
 *    POINTER   (ClientDataPtr,ClientData)
 * </FortranSyntax>
 *
 * @sa TecUtilEventAddPreDrawCallback(), TecUtilEventAddPostDrawCallback()
 */
typedef Boolean_t (STDCALL *DrawEventCallback_pf)(RedrawReason_e RedrawReason,
                                                  ArbParam_t     ClientData);


/**
 * Compares two strings from a list string. Note that either string may be NULL
 * as StringLists allow for NULL elements.
 *
 * @param String1
 *   String to compare against String2.
 * @param String2
 *   String to compare against String1.
 * @param ClientData
 *   Contextual information that was passed to the 'StringListSort' function.
 *
 * @return
 *   - A value less than zero if String1 is less than String2.
 *   - A value of zero if String1 is equal to String2.
 *   - A value greater than zero if String1 is greater than String2.
 */
typedef int (STDCALL *StringListStringComparator_pf)(const char* String1,
                                                     const char* String2,
                                                     ArbParam_t  ClientData);

/**
 * Gets a value at the specified point index using, if necessary, the private
 * client data retrieved from the field data handle.
 *
 * @par Note:
 *   This callback is called asynchronously. This callback should NOT
 *   lock/unlock Tecplot.
 *
 * @since
 *   10.0-3-128
 *
 * @param FD
 *   Field data handle for which to set the value.  This
 *   FieldValueGetFunction_pf must have been retrieved from this field data
 *   handle via TecUtilDataValueRefGetGetFunc.
 *
 * @param pt
 *   Zero-based index into the field data.
 *
 * @return
 *   Value for that index, always passed as a double precision floating-point
 *   value regardless of the data type of the field data handle.
 *
 * @sa TecUtilDataValueCustomLOD(), TecUtilDataValueGetClientData()
 */
typedef double(STDCALL *FieldValueGetFunction_pf)(const FieldData_pa FD,
                                                  LgIndex_t          pt);

/**
 * Sets a value at the specified index using the private client data retrieved
 * from the field data handle.
 *
 * @par Note:
 *   This callback is called asynchronously. This callback should NOT
 *   lock/unlock Tecplot.
 *
 * @since
 *   10.0-3-128
 *
 * @param FD
 *   Field data handle for which to set the value.  This
 *   FieldValueSetFunction_pf must have been retrieved from this field data
 *   handle via TecUtilDataValueRefGetSetFunc.
 *
 * @param pt
 *   Zero-based index into the field data.
 *
 * @param val
 *   New value for that index, always passed as a double precision
 *   floating-point value regardless of the data type of the field data handle.
 *
 * @sa TecUtilDataValueCustomLOD(), TecUtilDataValueGetClientData()
 */
typedef void (STDCALL *FieldValueSetFunction_pf)(FieldData_pa FD,
                                                 LgIndex_t    pt,
                                                 double       val);

/**
 * Callback responsible for loading the specified variable for Tecplot using
 * the private client data retrieved from the field data handle.
 *
 * @par Note:
 *   This callback is called asynchronously. With the exception of calls to
 *   modify the field data all calls back to Tecplot through the TecUtil layer
 *   should be limited to queries.
 *
 * @since
 *   11.0-0-001
 *
 * @param FieldData
 *   Field data handle of the variable load.
 *
 * @return
 *   TRUE if the variable was loaded, FALSE if unable to do so.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       LgIndex_t  NumValues;
 *       ... other information needed to load variable data
 *     } MyVariableClientData_s;
 *
 *   Boolean_t STDCALL MyVariableLoader(FieldData_pa FieldData)
 *   {
 *     REQUIRE(VALID_REF(FieldData));
 *
 *     MyVariableClientData_s *MyClientData = (MyVariableClientData_s *)TecUtilDataValueGetClientData(FieldData);
 *
 *     // open the data file
 *     FILE *MyDataFile = fopen(MyClientData->DataFileName, "rb");
 *     Boolean_t IsOk = (MyDataFile != NULL);
 *
 *     // seek to the place in the file where the variable data is located
 *     IsOk = IsOk && (fseek(MyDataFile, MyClientData->SeekOffset, SEEK_SET) == 0);
 *     if (IsOk)
 *       {
 *         // load the data into the variable's field data
 *         IsOk = ReadMyDataInfoVariable(MyDataFile, MyClientData, FieldData);
 *       }
 *
 *     // cleanup
 *     if (MyDataFile != NULL)
 *       fclose(MyDataFile);
 *
 *     ENSURE(VALID_BOOLEAN(IsOk));
 *     return IsOk;
 *   }
 * @endcode
 *
 * @sa TecUtilDataValueCustomLOD(), TecUtilDataValueGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandVarLoad_pf)(FieldData_pa FieldData);

/**
 * Callback responsible for performing private actions associated with a
 * variable being unloaded using the private client data retrieved from the
 * field data handle. Whenever possible the callback should honor Tecplot's
 * request to unload the variable by returning TRUE. This callback is
 * responsible for performing private actions associated with a variable being
 * unloaded.
 *
 * Most add-ons should simply supply NULL for this callback thereby instructing
 * Tecplot to handle the unloading (and subsequent reloading) of the variable
 * without the intervention of the add-on.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.0-0-001
 *
 * @param FieldData
 *   Field data handle of the variable Tecplot wants to unload.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       LgIndex_t  NumValues;
 *       ... other information needed to load variable data
 *     } MyVariableClientData_s;
 *
 *   Boolean_t STDCALL MyVariableUnload(FieldData_pa FieldData)
 *   {
 *     REQUIRE(VALID_REF(FieldData));
 *
 *     // We don't have any private data to cleanup (i.e in addition to the
 *     // private client data which we don't cleanup here) so all we have to do
 *     // is return TRUE or FALSE letting Tecplot know that it can or can not
 *     // unload the variable.
 *     Boolean_t Result = TRUE; // ...tell Tecplot to go ahead and unload the variable
 *
 *     ENSURE(VALID_BOOLEAN(Result));
 *     return Result;
 *   }
 * @endcode
 *
 * @return
 *   TRUE if the variable can be unloaded, FALSE otherwise. The add-on should
 *   if at all possible honor the request to unload the variable. Most add-ons
 *   should return TRUE.
 *
 * @sa TecUtilDataValueCustomLOD(), TecUtilDataValueGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandVarUnload_pf)(FieldData_pa FieldData);

/**
 * Callback responsible for performing private actions associated with a
 * variable being cleaned up using the private client data retrieved from the
 * field data handle. Most add-ons will need to register this callback in order
 * to cleanup privately allocated client data.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.0-0-001
 *
 * @param FieldData
 *   Field data handle of the variable being cleaned up.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       LgIndex_t  NumValues;
 *       ... other information needed to load variable data
 *     } MyVariableClientData_s;
 *
 *   void STDCALL MyVariableCleanup(FieldData_pa FieldData)
 *   {
 *     REQUIRE(VALID_REF(FieldData));
 *
 *     MyVariableClientData_s *MyClientData = (MyVariableClientData_s *)TecUtilDataValueGetClientData(FieldData);
 *
 *     // cleanup privately allocated resources
 *     free(MyClientData->DataFileName);
 *     free(MyClientData);
 *   }
 * @endcode
 *
 * @sa TecUtilDataValueCustomLOD(), TecUtilDataValueGetClientData()
 */
typedef void (STDCALL *LoadOnDemandVarCleanup_pf)(FieldData_pa FieldData);

/**
 * Callback responsible for loading the specified node mapping for Tecplot
 * using the private client data retrieved from the node mapping handle.
 *
 * @par Note:
 *   This callback is called asynchronously. With the exception of calls to
 *   modify the node mapping, all calls back to Tecplot through the TecUtil
 *   layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param NodeMap
 *   Handle of the node mapping.
 *
 * @return
 *   TRUE if the node mapping was loaded, FALSE if unable to do so.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ... other information needed to load node map data
 *     } MyNodeMapClientData_s;
 *
 *   Boolean_t STDCALL MyNodeMapLoader(NodeMap_pa NodeMap)
 *   {
 *     REQUIRE(VALID_REF(NodeMap));
 *
 *     MyNodeMapClientData_s *MyClientData =
 *             (MyNodeMapClientData_s *)TecUtilDataNodeGetClientData(NodeMap);
 *
 *     // open the data file
 *     FILE *MyDataFile = fopen(MyClientData->DataFileName, "rb");
 *     Boolean_t IsOk = (MyDataFile != NULL);
 *
 *     // seek to the place in the file where the node map data is located
 *     IsOk = IsOk && (fseek(MyDataFile, MyClientData->SeekOffset, SEEK_SET) == 0);
 *     if (IsOk)
 *       {
 *         // load the data into the zone's node map
 *         IsOk = ReadMyNodeMapDataIntoZone(MyDataFile, MyClientData, NodeMap);
 *       }
 *
 *     // cleanup
 *     if (MyDataFile != NULL)
 *       fclose(MyDataFile);
 *
 *     ENSURE(VALID_BOOLEAN(IsOk));
 *     return IsOk;
 *   }
 * @endcode
 *
 * @sa TecUtilDataNodeCustomLOD(), TecUtilDataNodeGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandNodeMapLoad_pf)(NodeMap_pa NodeMap);

/**
 * Callback responsible for performing private actions associated with a
 * node mapping being unloaded using the private client data retrieved from the
 * node mapping handle. Whenever possible the callback should honor Tecplot's
 * request to unload the node mapping by returning TRUE.
 *
 * Most add-ons should simply supply NULL for this callback thereby instructing
 * Tecplot to handle the unloading (and subsequent reloading) of the node mapping
 * without the intervention of the add-on.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param NodeMap
 *   Node mapping handle of the node mapping Tecplot wants to unload.
 *
 * @code
 *   Boolean_t STDCALL MyNodeMapUnload(NodeMap_pa NodeMap)
 *   {
 *     REQUIRE(VALID_REF(NodeMap));
 *
 *     // We don't have any private data to cleanup (i.e in addition to the
 *     // private client data which we don't cleanup here) so all we have to do
 *     // is return TRUE or FALSE letting Tecplot know that it can or can not
 *     // unload the variable.
 *     Boolean_t Result = TRUE; // ...tell Tecplot to go ahead and unload the node mapping
 *
 *     ENSURE(VALID_BOOLEAN(Result));
 *     return Result;
 *   }
 * @endcode
 *
 * @return
 *   TRUE if the node mapping can be unloaded, FALSE otherwise. The add-on should
 *   if at all possible honor the request to unload the node mapping. Most add-ons
 *   should return TRUE.
 *
 * @sa TecUtilDataNodeCustomLOD(), TecUtilDataNodeGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandNodeMapUnload_pf)(NodeMap_pa NodeMap);

/**
 * Callback responsible for performing private actions associated with a
 * node mapping being cleaned up using the private client data retrieved from the
 * node mapping handle. Most add-ons will need to register this callback in order
 * to cleanup privately allocated client data.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param NodeMap
 *   Node Mapping data handle of the node mapping being cleaned up.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ... other information needed to load node map data
 *     } MyNodeMapClientData_s;
 *
 *   void STDCALL MyNodeMapCleanup(NodeMap_pa NodeMap)
 *   {
 *     REQUIRE(VALID_REF(NodeMap));
 *
 *     MyNodeMapClientData_s *MyClientData = (MyNodeMapClientData_s *)TecUtilDataNodeGetClientData(NodeMap);
 *
 *     // cleanup privately allocated resources
 *     free(MyClientData->DataFileName);
 *     free(MyClientData);
 *   }
 * @endcode
 *
 * @sa TecUtilDataNodeCustomLOD(), TecUtilDataNodeGetClientData()
 */
typedef void (STDCALL *LoadOnDemandNodeMapCleanup_pf)(NodeMap_pa NodeMap);

/**
 * Callback responsible for loading the specified face neighbor for Tecplot
 * using the private client data retrieved from the face neighbor handle.
 *
 * @par Note:
 *   This callback is called asynchronously. With the exception of calls to
 *   modify the face neighbors, all calls back to Tecplot through the TecUtil
 *   layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param FaceNeighbor
 *   Handle of the face neighbors.
 *
 * @return
 *   TRUE if the face neighbors was loaded, FALSE if unable to do so.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ...other information needed to load face neighbor data
 *     } MyFaceNeighborClientData_s;
 *
 *   Boolean_t STDCALL MyFaceNeighborLoader(FaceNeighbor_pa FaceNeighbor)
 *   {
 *     REQUIRE(VALID_REF(FaceNeighbor));
 *
 *     MyFaceNeighborClientData_s *MyClientData =
 *             (MyFaceNeighborClientData_s*)TecUtilDataFaceNbrGetClientData(FaceNeighbor);
 *
 *     // open the data file
 *     FILE *MyDataFile = fopen(MyClientData->DataFileName, "rb");
 *     Boolean_t IsOk = (MyDataFile != NULL);
 *
 *     // seek to the place in the file where the face neighbor data is located
 *     IsOk = IsOk && (fseek(MyDataFile, MyClientData->SeekOffset, SEEK_SET) == 0);
 *     if (IsOk)
 *       {
 *         // load the data into the zone's face neighbor
 *         IsOk = ReadMyFaceNeighborDataIntoZone(MyDataFile, MyClientData, FaceNeighbor);
 *       }
 *
 *     // cleanup
 *     if (MyDataFile != NULL)
 *       fclose(MyDataFile);
 *
 *     ENSURE(VALID_BOOLEAN(IsOk));
 *     return IsOk;
 *   }
 * @endcode
 *
 * @sa TecUtilDataFaceNbrCustomLOD(), TecUtilDataFaceNbrGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandFaceNeighborLoad_pf)(FaceNeighbor_pa FaceNeighbor);

/**
 * Callback responsible for performing private actions associated with a
 * face neighbors being unloaded using the private client data retrieved from
 * the face neighbor handle. Whenever possible the callback should honor
 * Tecplot's request to unload the face neighbors by returning TRUE.
 *
 * Most add-ons should simply supply NULL for this callback thereby instructing
 * Tecplot to handle the unloading (and subsequent reloading) of the face
 * neighbors without the intervention of the add-on.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param FaceNeighbor
 *   Face neighbor handle of the face neighbors Tecplot wants to unload.
 *
 * @code
 *   Boolean_t STDCALL MyFaceNeighborUnload(FaceNeighbor_pa FaceNeighbor)
 *   {
 *     REQUIRE(VALID_REF(FaceNeighbor));
 *
 *     // We don't have any private data to cleanup (i.e in addition to the
 *     // private client data which we don't cleanup here) so all we have to do
 *     // is return TRUE or FALSE letting Tecplot know that it can or can not
 *     // unload the variable.
 *     Boolean_t Result = TRUE; // ...tell Tecplot to go ahead and unload the face neighbors
 *
 *     ENSURE(VALID_BOOLEAN(Result));
 *     return Result;
 *   }
 * @endcode
 *
 * @return
 *   TRUE if the face neighbors can be unloaded, FALSE otherwise. The add-on
 *   should if at all possible honor the request to unload the face neighbors.
 *   Most add-ons should return TRUE.
 *
 * @sa TecUtilDataFaceNbrCustomLOD(), TecUtilDataFaceNbrGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandFaceNeighborUnload_pf)(FaceNeighbor_pa FaceNeighbor);

/**
 * Callback responsible for performing private actions associated with a face
 * neighbors being cleaned up using the private client data retrieved from the
 * face neighbor handle. Most add-ons will need to register this callback in
 * order to cleanup privately allocated client data.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.3-0-010
 *
 * @param FaceNeighbor
 *   Face neighbor data handle of the Face neighbors being cleaned up.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ... other information needed to load face neighbor data
 *     } MyFaceNeighborClientData_s;
 *
 *   void STDCALL MyFaceNeighborCleanup(FaceNeighbor_pa FaceNeighbor)
 *   {
 *     REQUIRE(VALID_REF(FaceNeighbor));
 *
 *     MyFaceNeighborClientData_s *MyClientData = (MyFaceNeighborClientData_s *)TecUtilDataFaceNbrGetClientData(FaceNeighbor);
 *
 *     // cleanup privately allocated resources
 *     free(MyClientData->DataFileName);
 *     free(MyClientData);
 *   }
 * @endcode
 *
 * @sa TecUtilDataFaceNbrCustomLOD(), TecUtilDataFaceNbrGetClientData()
 */
typedef void (STDCALL *LoadOnDemandFaceNeighborCleanup_pf)(FaceNeighbor_pa FaceNeighbor);

/**
 * Callback responsible for loading the specified face mapping for Tecplot
 * using the private client data retrieved from the face mapping handle.
 *
 * @par Note:
 *   This callback is called asynchronously. With the exception of calls to
 *   modify the face mapping, all calls back to Tecplot through the TecUtil
 *   layer should be limited to queries.
 *
 * @since
 *   11.2-1-0
 *
 * @param FaceMap
 *   Handle of the face mapping.
 *
 * @return
 *   TRUE if the face mapping was loaded, FALSE if unable to do so.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ... other information needed to load face map data
 *     } MyFaceMapClientData_s;
 *
 *   Boolean_t STDCALL MyFaceMapLoader(FaceMap_pa FaceMap)
 *   {
 *     REQUIRE(VALID_REF(FaceMap));
 *
 *     MyFaceMapClientData_s *MyClientData =
 *             (MyFaceMapClientData_s *)TecUtilDataFaceMapGetClientData(FaceMap);
 *
 *     // open the data file
 *     FILE *MyDataFile = fopen(MyClientData->DataFileName, "rb");
 *     Boolean_t IsOk = (MyDataFile != NULL);
 *
 *     // seek to the place in the file where the face map data is located
 *     IsOk = IsOk && (fseek(MyDataFile, MyClientData->SeekOffset, SEEK_SET) == 0);
 *     if (IsOk)
 *       {
 *         // load the data into the zone's face map
 *         IsOk = ReadMyFaceMapDataIntoZone(MyDataFile, MyClientData, FaceMap);
 *       }
 *
 *     // cleanup
 *     if (MyDataFile != NULL)
 *       fclose(MyDataFile);
 *
 *     ENSURE(VALID_BOOLEAN(IsOk));
 *     return IsOk;
 *   }
 * @endcode
 *
 * @sa TecUtilDataFaceMapCustomLOD(), TecUtilDataFaceMapGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandFaceMapLoad_pf)(FaceMap_pa FaceMap);

/**
 * Callback responsible for performing private actions associated with a
 * face mapping being unloaded using the private client data retrieved from the
 * face mapping handle. Whenever possible the callback should honor Tecplot's
 * request to unload the face mapping by returning TRUE.
 *
 * Most add-ons should simply supply NULL for this callback thereby instructing
 * Tecplot to handle the unloading (and subsequent reloading) of the face mapping
 * without the intervention of the add-on.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.2-1-0
 *
 * @param FaceMap
 *   Face mapping handle of the face mapping Tecplot wants to unload.
 *
 * @code
 *   Boolean_t STDCALL MyFaceMapUnload(FaceMap_pa FaceMap)
 *   {
 *     REQUIRE(VALID_REF(FaceMap));
 *
 *     // We don't have any private data to cleanup (i.e in addition to the
 *     // private client data which we don't cleanup here) so all we have to do
 *     // is return TRUE or FALSE letting Tecplot know that it can or can not
 *     // unload the variable.
 *     Boolean_t Result = TRUE; // ...tell Tecplot to go ahead and unload the face mapping
 *
 *     ENSURE(VALID_BOOLEAN(Result));
 *     return Result;
 *   }
 * @endcode
 *
 * @return
 *   TRUE if the face mapping can be unloaded, FALSE otherwise. The add-on should
 *   if at all possible honor the request to unload the face mapping. Most add-ons
 *   should return TRUE.
 *
 * @sa TecUtilDataFaceMapCustomLOD(), TecUtilDataFaceMapGetClientData()
 */
typedef Boolean_t (STDCALL *LoadOnDemandFaceMapUnload_pf)(FaceMap_pa FaceMap);

/**
 * Callback responsible for performing private actions associated with a
 * face mapping being cleaned up using the private client data retrieved from the
 * face mapping handle. Most add-ons will need to register this callback in order
 * to cleanup privately allocated client data.
 *
 * @par Note:
 *   This callback is called asynchronously. All calls back to Tecplot through
 *   the TecUtil layer should be limited to queries.
 *
 * @since
 *   11.2-1-0
 *
 * @param FaceMap
 *   Face Mapping data handle of the face mapping being cleaned up.
 *
 * @code
 *   typedef struct
 *     {
 *       char      *DataFileName;
 *       long       SeekOffset;
 *       ... other information needed to load face map data
 *     } MyFaceMapClientData_s;
 *
 *   void STDCALL MyFaceMapCleanup(FaceMap_pa FaceMap)
 *   {
 *     REQUIRE(VALID_REF(FaceMap));
 *
 *     MyFaceMapClientData_s *MyClientData = (MyFaceMapClientData_s *)TecUtilDataFaceMapGetClientData(FaceMap);
 *
 *     // cleanup privately allocated resources
 *     free(MyClientData->DataFileName);
 *     free(MyClientData);
 *   }
 * @endcode
 *
 * @sa TecUtilDataFaceMapCustomLOD(), TecUtilDataFaceMapGetClientData()
 */
typedef void (STDCALL *LoadOnDemandFaceMapCleanup_pf)(FaceMap_pa FaceMap);


/**
 * ExtractDestination functions are called upon successful completion of an
 * extract polyline or extract discrete points operation.
 *
 * @param NumPts
 *   Number of points extracted.
 *
 * @param XValues
 *   Double precision array of X-Coordinates of the extracted polyline.
 *
 * @param YValues
 *   Double precision array of Y-Coordinates of the extracted polyline.
 *
 * <FortranSyntax>
 *   INTEGER*4 FUNCTION MyExtractDestinationCallback(
 *  &                   NumPts,
 *  &                   XValues,
 *  &                   YValues)
 *   INTEGER*4 NumPts
 *   REAL*8    XValues
 *   REAL*8    YValues
 * </FortranSyntax>
 */
typedef void (STDCALL *ExtractDestination_pf)(LgIndex_t NumPts,
                                              double*   XValues,
                                              double*   YValues);



/**
 * SelectFileOptionsCallback Functions are called when the
 * "Options" button is pressed in the modal file selection
 * dialog.
 *
 * <FortranSyntax>
 *   SUBROUTINE MySelectFileOptionsCallback()
 * </FortranSyntax>
 */
typedef void (STDCALL *SelectFileOptionsCallback_pf)(void);




/**
 * Post data load instruction callback for "Converter-Plus" addons.
 *
 * @param PreviousInstructions
 *   The previous set of instructions used by the converter.
 *
 * @param PreviousRawData
 *   The previous raw data associated with the instructions.
 *
 * @param PreviousZones
 *   Set of zones loaded with the previous instructions.
 *
 * <FortranSyntax>
 *    SUBROUTINE MyConverterPostReadCallback(
 *   &                   PreviousInstructions,
 *   &                   PreviousRawData,
 *   &                   PreviousZones)
 *    CHARACTER*(*)   CommandString
 *    CHARACTER*(*)   ErrMsgString
 *    POINTER         (PreviousZones,DummyPreviousZonesData)
 * </FortranSyntax>
 *
 */
typedef void (STDCALL *ConverterPostReadCallback_pf)(const char*  PreviousInstructions,
                                                     const char*  PreviousRawData,
                                                     const Set_pa PreviousZones);


/**
 * Callback registered by your addon to convert a foreign datafile into a
 * Tecplot Binary datafile format.
 *
 * @return
 *   Return TRUE if the conversion is successful. Otherwise return FALSE.
 *   If FALSE is returned then *MessageString is assumed to contain an error
 *   message.
 *
 * @param DataFName
 *   Name of the original foreign data file to be converted.
 *
 * @param TempBinFName
 *   Name of the temporary binary datafile that is created (by your converter).
 *
 * @param MessageString
 *   Reference to a string. If an error occurs during conversion allocate space
 *   for an error message and copy the message string into that allocated
 *   space otherwise be sure to assign *MessageString to NULL. If
 *   *MessageString is non NULL Tecplot will release the allocated memory when
 *   finished.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyDataSetConverterCallback(
 *   &                   DataFName,
 *   &                   TempBinFName,
 *   &                   MessageString)
 *    CHARACTER*(*)   DataFName
 *    CHARACTER*(*)   TempBinFName
 *    CHARACTER*(*)   MessageString
 * </FortranSyntax>
 *
 */
typedef Boolean_t (STDCALL *DataSetConverter_pf)(char*           DataFName,
                                                 char*           TempBinFName,
                                                 TP_GIVES char** MessageString);







/**
 * Callback registered by your addon to process foreign loader instructions.
 * When called, it must parse the supplied instructions and load the data into Tecplot.
 *
 * @return
 *   Return TRUE if the data is loaded successfully. Otherwise, FALSE.
 *
 * @param Instructions
 *   This contains all of the instructions needed to load the data.
 *
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyDataSetLoaderCallback(
 *   &                   Instructions)
 *    POINTER        (Instructions,DummyInstructionsData)
 * </FortranSyntax>
 */
typedef Boolean_t (STDCALL *DataSetLoader_pf)(StringList_pa Instructions);





/**
 * Callback used to provide the ability to override data loader instructions
 * while processing a layout.
 *
 * @return
 *   Return TRUE if the instructions are successfully replaced or left alone.
 *   Return FALSE if the user cancels the operation.
 *
 * @param Instructions
 *   The original instructions needed to load the data.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyDataSetLoaderInstOverCallback(
 *   &                   Instructions)
 *    POINTER        (Instructions,DummyInstructionsData)
 * </FortranSyntax>
 *
 */
typedef Boolean_t (STDCALL *DataSetLoaderInstructionOverride_pf)(StringList_pa  Instructions);



/**
 *  Callback used to assign extended curve settings.
 *  This is called when the user presses the "Curve Settings"
 *  button in the mapping style dialog.
 *
 *  @param LineMapSet
 *    Set of line maps currently selected.
 *  @param SelectedLineMapSettings
 *    A string list of the curve settings for the Line-maps that are selected in the
 *    Line mappings dialog.
 *
 * <FortranSyntax>
 *   SUBROUTINE MyGetCurveSettingsCallback(
 *  &                LineMapSet,
 *  &                SelectedLineMapSettings)
 *    POINTER    (LineMapSet,DummyLineMapData)
 *    POINTER    (SelectedLineMapSettings,DummyLineMapSettings)
 * </FortranSyntax>
 */
typedef void (STDCALL *GetCurveSettingsCallback_pf)(Set_pa        LineMapSet,
                                                    StringList_pa SelectedLineMapSettings);




/**
 * Callback function that returns an abbreviated version of the curve settings
 * for a particular Line Map for display in the Line Mappings dialog.
 *
 * @param LineMap
 *   The map number that is currently being operated on.
 * @param CurveSettings
 *   The string that Tecplot maintains which contains the extended curve fit
 *   settings for the current Line-map. This argument may be NULL indicating
 *   that defaults should be used.
 * @param AbbreviatedSettings
 *   The short form of the CurveSettings that is allocated and returned from
 *   your function and used by Tecplot. This must be allocated by the addon
 *   using TecUtilStringAlloc().
 *
 * <FortranSyntax>
 *   SUBROUTINE MyGetAbrevSettingsStringCallback(
 *  &                LineMap,
 *  &                CurveSettings,
 *  &                AbbreviatedSettings),
 *    INTEGER*4  LineMap
 *    CHARACTER*(*) CurveSettings
 *    CHARACTER*(*) AbbreviatedSettings
 * </FortranSyntax>
 */
typedef void (STDCALL *GetAbbreviatedSettingsStringCallback_pf)(EntIndex_t      LineMap,
                                                                char*           CurveSettings,
                                                                TP_GIVES char** AbbreviatedSettings);




/**
 * This function returns a string (CurveInfoString) for Tecplot to display
 * information about a particular curve in the curve info dialog.
 *
 * @param RawIndV
 *   The handle to the raw field data of the independent variable.
 * @param RawDepV
 *   The handle to the raw field data of the dependent variable.
 * @param IndVCoordScale
 *   An enumerated variable whose values are Scale_linear when the independent variable
 *   axis has a linear scale and Scale_log when it has a log scale.
 * @param DepVCoordScale
 *   An enumerated variable whose values are Scale_linear when the dependent variable axis
 *   has a linear scale and Scale_log when it has a log scale.
 * @param NumRawPts
 *   number of raw field data values.
 * @param LineMap
 *   The map number that is currently being operated on.
 * @param CurveSettings
 *   The curve settings string for the current Line-map. This argument may be
 *   NULL indicating that defaults should be used.
 * @param CurveInfoString
 *   The string that is allocated and returned by your function and be
 *   presented in the Data/XY-Plot Curve Info dialog. The CurveInfoString must
 *   be allocated by the addon using TecUtilStringAlloc().
 *
 * @return
 *   Return TRUE if the curve info string can be generated, otherwise FALSE.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyGetCurveInfoStringCallback(
 *   &                   RawIndV,
 *   &                   RawDepV,
 *   &                   IndVCoordScale,
 *   &                   DepVCoordScale,
 *   &                   NumRawPts,
 *   &                   LineMap,
 *   &                   CurveSettings,
 *   &                   CurveInfoString)
 *    POINTER       (RawIndV,DummyRawIndVData)
 *    POINTER       (RawDepV,DummyRawDepVData)
 *    INTEGER*4     IndVCoordScale
 *    INTEGER*4     DepVCoordScale
 *    INTEGER*4     NumRawPts
 *    INTEGER*4     LineMap
 *    CHARACTER*(*) CurveSettings
 *    CHARACTER*(*) CurveInfoString
 * </FortranSyntax>
 */
typedef Boolean_t (STDCALL *GetCurveInfoStringCallback_pf)(FieldData_pa    RawIndV,
                                                           FieldData_pa    RawDepV,
                                                           CoordScale_e    IndVCoordScale,
                                                           CoordScale_e    DepVCoordScale,
                                                           LgIndex_t       NumRawPts,
                                                           EntIndex_t      LineMap,
                                                           char*           CurveSettings,
                                                           TP_GIVES char** CurveInfoString);

/**
 * Callback function used to calculate data points for an extended curve fit.
 *
 * @return
 *   Return TRUE if the curve can be calculated, otherwise FALSE.
 *
 * @param RawIndV
 *   The handle to the raw field data of the independent variable.
 * @param RawDepV
 *   The handle to the raw field data of the dependent variable.
 * @param IndVCoordScale
 *   An enumerated variable whose values are Scale_linear when the independent variable
 *   axis has a linear scale and Scale_log when it has a log scale.
 * @param DepVCoordScale
 *   An enumerated variable whose values are Scale_linear when the dependent variable axis
 *   has a linear scale and Scale_log when it has a log scale.
 * @param NumRawPts
 *   number of raw field data values.
 * @param NumCurvePts
 *   The number of points that will construct the curve fit.
 * @param LineMap
 *   The line map to operated on.
 * @param CurveSettings
 *   The curve settings string for the current Line-map. This argument may be
 *   NULL indicating that defaults should be used.
 * @param IndCurveValues
 *   A pre-allocated array of size NumCurvePts which the addon will populate with
 *   the independent values for the curve fit
 * @param DepCurveValues.
 *   A pre-allocated array of size NumCurvePts which the add-on will populate
 *   with the dependent values for the curve fit.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyGetLinePlotDataPointsCallback(
 *   &                   RawIndV,
 *   &                   RawDepV,
 *   &                   IndVCoordScale,
 *   &                   DepVCoordScale,
 *   &                   NumRawPts,
 *   &                   NumCurvePts,
 *   &                   LineMap,
 *   &                   CurveSettings,
 *   &                   IndCurveValues,
 *   &                   DepCurveValues)
 *    POINTER       (RawIndV,DummyRawIndVData)
 *    POINTER       (RawDepV,DummyRawDepVData)
 *    INTEGER*4     IndVCoordScale
 *    INTEGER*4     DepVCoordScale
 *    INTEGER*4     NumRawPts
 *    INTEGER*4     NumCurvePts
 *    INTEGER*4     LineMap
 *    CHARACTER*(*) CurveSettings
 *    REAL*8        IndCurveValues()
 *    REAL*8        DepCurveValues()
 * </FortranSyntax>
 */
typedef Boolean_t (STDCALL *GetLinePlotDataPointsCallback_pf)(FieldData_pa   RawIndV,
                                                              FieldData_pa   RawDepV,
                                                              CoordScale_e   IndVCoordScale,
                                                              CoordScale_e   DepVCoordScale,
                                                              LgIndex_t      NumRawPts,
                                                              LgIndex_t      NumCurvePts,
                                                              EntIndex_t     LineMap,
                                                              char*          CurveSettings,
                                                              TP_OUT double* IndCurveValues,
                                                              TP_OUT double* DepCurveValues);
#if defined EXPORT_DEPRECATED_INTERFACES_TO_ADK_ONLY
/**
 * @deprecated
 *     Please use \ref GetLinePlotDataPointsCallback_pf instead.
 */
typedef GetLinePlotDataPointsCallback_pf GetXYDataPointsCallback_pf;
#endif




/**
 * A Callback function used to obtain an interpolated dependent value for an
 * extended curve fit given an independent value.
 *
 * @return
 *   Return TRUE if it is possible to obtain the interpolated value, otherwise FALSE.
 *
 * @param RawIndV
 *   handle to the raw field data of the independent variable.
 * @param RawDepV
 *   The handle to the raw field data of the dependent variable.
 * @param IndVCoordScale
 *   An enumerated variable whose values are Scale_linear when the independent variable
 *   axis has a linear scale and Scale_log when it has a log scale.
 * @param DepVCoordScale
 *   An enumerated variable whose values are Scale_linear when the dependent variable axis
 *   has a linear scale and Scale_log when it has a log scale.
 * @param NumRawPts
 *   The number of field data values.
 * @param NumCurvePts
 *   The number of points used to construct the curve fit.
 * @param LineMapNum
 *   The line map number currently being operated on.
 * @param CurveSettings
 *   The curve settings string for the current Line-map. This argument may be
 *   NULL indicating that defaults should be used.
 * @param ProbeIndValue
 *   The independent value location of the probe (supplied).
 * @param ProbeDepValue
 *   Reference to the calculated dependent value location of the probe.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyGetProbeValueCallback(
 *   &                   RawIndV,
 *   &                   RawDepV,
 *   &                   IndVCoordScale,
 *   &                   DepVCoordScale,
 *   &                   NumRawPts,
 *   &                   NumCurvePts,
 *   &                   LineMapNum,
 *   &                   CurveSettings,
 *   &                   CurveInfoString,
 *   &                   ProbeIndValue,
 *   &                   ProbeDepValue)
 *    POINTER       (RawIndV,DummyRawIndVData)
 *    POINTER       (RawDepV,DummyRawDepVData)
 *    INTEGER*4     IndVCoordScale
 *    INTEGER*4     DepVCoordScale
 *    INTEGER*4     NumRawPts
 *    INTEGER*4     NumCurvePts
 *    INTEGER*4     LineMapNum
 *    CHARACTER*(*) CurveSettings
 *    REAL*8        ProbeIndValue
 *    REAL*8        ProbeDepValue
 * </FortranSyntax>
 *
 */
typedef Boolean_t (STDCALL *GetProbeValueCallback_pf)(FieldData_pa   RawIndV,
                                                      FieldData_pa   RawDepV,
                                                      CoordScale_e   IndVCoordScale,
                                                      CoordScale_e   DepVCoordScale,
                                                      LgIndex_t      NumRawPts,
                                                      LgIndex_t      NumCurvePts,
                                                      EntIndex_t     LineMapNum,
                                                      char*          CurveSettings,
                                                      double         ProbeIndValue,
                                                      TP_OUT double* ProbeDepValue);



#if defined MSWIN
typedef Boolean_t (STDCALL *PreTranslateMessage_pf)(MSG *pMsg);
#endif


/**
 * Callback function pointer for providing a Dynamic Axis labels.
 * @since
 *    10.0-6-015
 * @param Value
 *     Value that corresponds to a tick label that will be drwan.
 *
 * @param ClientData
 *     Convenience storage of user client data.
 *
 * @param LabelString
 *    Output label for the tick mark.
 *    This must be allocated by the addon using TecUtilStringAlloc().
 *
 * @return
 *    Returns TRUE if the LabelString has been successfully allocated.
 *    Otherwise, FALSE is returned.
 */
typedef Boolean_t (STDCALL *DynamicLabelCallback_pf)(double          Value,
                                                     ArbParam_t      ClientData,
                                                     TP_GIVES char** LabelString);

/**
 * This is called when Tecplot is idle.
 *
 * @par Note:
 *   Tecplot is never idle when running in batch mode (with the -b flag).
 *
 * @param ClientData
 *   Arbitrary client data.
 *
 * <FortranSyntax>
 *    INTEGER*4 FUNCTION MyOnIdleCallback(
 *   &                     ClientDataPtr)
 *    POINTER (ClientDataPtr,DummyClientData)
 * </FortranSyntax>
 *
 */
typedef void (STDCALL *OnIdleCallback_pf)(ArbParam_t ClientData);

/**
 * Callback responsible for executing the specified script file.
 *
 * @since
 *     11.0-2-005
 *
 * @param ScriptFileName
 *     Relative or absolute file name of the script to execute. If the path
 *     is relative it is relative to the current working directory.
 * @param ClientData
 *     Client data registered with the callback.
 *
 * @return
 *     TRUE if the script executed successfully, FALSE otherwise.
 *
 * @sa TecUtilScriptExecRegisterCallback
 */
typedef Boolean_t (STDCALL *ScriptExecCallback_pf)(const char *ScriptFileName,
                                                   ArbParam_t  ClientData);

/**
 * This callback is called by TecUtilLineSegProbe() each time a cell face
 * is about to be passed through in the course of a probe.
 *
 * @return
 *   Return FALSE if you want to stop the probe at the current face.
 *   Return TRUE if you want the probe to progress through the face toward
 *   the specified end point.
 *
 * @param WhichEndingPosition
 *   Which ending position of the call to TecUtilLineSegProbe() is currently
 *   being probed. It will be between 1 and NumEndingPositions, inclusive, where
 *   NumEndingPositions is the number of ending positions passed into
 *   TecUtilLineSegProbe().
 *
 * @param Zone
 *   The zone that the probe is currently passing through.
 *
 * @param Cell
 *   The number of the cell that the probe is currently passing through.
 *
 * @param Face
 *   The face number of the face that the probe is about to pass through.
 *   For ordered and classic finite-element zones, this will be between
 *   1 and the number of faces per cell, inclusive. For polygonal and
 *   polyhedral zones, this is the overall zone face number, and will be
 *   between 1 and the total number of faces in the zone, inclusive.
 *
 * @param Position
 *   An array contining the X, Y, Z location of the point at which the
 *   probe trajectory has intercepted the face. In 2-D, it contains X and Y.
 *
 * @param ClientData
 *   The client data that was passed into TecUtilLineSegProbe().
 *
 * <FortranSyntax>
 * INTEGER*4 FUNCTION LineSegProbeCallback()
 * </FortranSyntax>
 */
typedef Boolean_t (STDCALL * LineSegProbeCallback_pf)(LgIndex_t         WhichEndingPosition,
                                                      EntIndex_t        Zone,
                                                      LgIndex_t         Cell,
                                                      LgIndex_t         Face,
                                                      TP_GIVES double * Position,
                                                      ArbParam_t        ClientData);

/* BEGINREMOVEFROMADDON */
#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#if 0 /* NOTUSED */
#endif
#if !defined NO_ASSERTS
#endif
#if defined MSWIN
#endif /* MSWIN */
#if !defined (MSWIN)
#endif
#if defined Q_MAINMODULE
#else
#endif
#if 0 /* NOTUSED */
#endif
#endif /* TECPLOTKERNEL */

#if defined TECPLOTKERNEL
/* CORE SOURCE CODE REMOVED */
#endif /* TECPLOTKERNEL */


/* ENDREMOVEFROMADDON */
struct _ViewState_a;
typedef struct _ViewState_a *SavedView_pa, *ViewState_pa;

/* define Tecplot support email address in one place */
static char const* const tecplotSupportEmailAddress = "support@tecplot.com";

#endif /* _GLOBAL_H */
