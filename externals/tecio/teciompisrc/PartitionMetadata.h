 #pragma once
#include "ClassMacros.h"
#include "SzlFileLoader.h"
namespace tecplot { namespace ___3933 { class PartitionMetadata { UNCOPYABLE_CLASS(PartitionMetadata) public: ___1392 m_cszConnectivityFileLocs; ___1392 m_nszConnectivityFileLocs; ___2090::___2980 m_numRefPartitions; PartitionArray           m_refPartitions; RefSubzoneOffsetArray m_cszNumRefNszs; RefSubzoneOffsetArray m_nszNumRefCszs; UInt8Array m_cszIncludesPtn; UInt8Array m_nszIncludesPtn; ___1392 m_cszMinMaxFileLocs; ___1392 m_nszMinMaxFileLocs; ___1392 m_szDataStartFileLocs; FileLoc2DArray m_varSzFileLoc; public: PartitionMetadata() {} }; }}
