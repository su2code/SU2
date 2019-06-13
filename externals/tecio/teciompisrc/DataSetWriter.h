 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_ptr.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
#include "ZoneInfoCache.h"
#include "ZoneVarMetadata.h"
namespace tecplot { namespace ___3933 { class ___37; class FileWriterInterface; class ItemSetIterator; class DataSetWriter { public: DataSetWriter( ___37*               ___36, ___3501                    vars, ___3501                    ___4671, ___1844 const&                maxIJKSubzoneSize, ___2090::ItemOffset_t maxFESubzoneSize, bool                      flushToDisk = false); virtual ~DataSetWriter(); virtual ___372 writeDataSet( FileWriterInterface& szpltFile, ___1392&        szpltZoneHeaderFileLocs); void replaceDataSource( ___37* ___36, ___3501      vars, ___3501      ___4671); void clear( int32_t        numZonesToRetain, int32_t const* zonesToRetain); ___4707 const& ___4706() { return *m_zoneVarMetadata; } protected: void getZoneSharing( std::vector<___372>& ___4564, ___372&              ___4499, ___4636             zone, ___4636             ___341, DataFileType_e          ___844) const; ___37*                        ___2337; boost::scoped_ptr<ItemSetIterator> m_varIter; boost::scoped_ptr<ItemSetIterator> m_zoneIter; ZoneInfoCache                      ___2680; boost::scoped_ptr<___4707> m_zoneVarMetadata; bool const                         m_flushingToDisk; }; }}
