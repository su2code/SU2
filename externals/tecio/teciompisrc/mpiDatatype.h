 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <mpi.h>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
 #if defined MAC64 
 #if !defined MPI_INT8_T
 #define MPI_INT8_T MPI_CHAR
 #endif
 #if !defined MPI_INT16_T
 #define MPI_INT16_T MPI_SHORT
 #endif
 #if !defined MPI_UINT16_T
 #define MPI_UINT16_T MPI_UNSIGNED_SHORT
 #endif
 #if !defined MPI_INT32_T
 #define MPI_INT32_T MPI_INT
 #endif
 #if !defined MPI_UINT32_T
 #define MPI_UINT32_T MPI_UNSIGNED
 #endif
 #if !defined MPI_INT64_T
 #define MPI_INT64_T MPI_LONG_LONG
 #endif
 #if !defined MPI_UINT64_T
 #define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
 #endif
 #endif
template <typename T> MPI_Datatype mpiDatatype();
