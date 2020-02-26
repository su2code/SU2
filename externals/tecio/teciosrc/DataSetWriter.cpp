#include "DataSetWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/unordered_set.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "CodeContract.h"
#include "FileWriterInterface.h"
#include "ItemSetIterator.h"
#include "SzlFileLoader.h"
#include "ZoneWriterAbstract.h"
#include "ZoneWriterFactory.h"
#include "zoneUtil.h"
#include "ZoneVarMetadata.h"
namespace tecplot { namespace ___3933 { DataSetWriter::DataSetWriter( ___37*               ___36, ___3501                    vars, ___3501                    ___4671, ___1844 const&                maxIJKSubzoneSize, ___2090::ItemOffset_t maxFESubzoneSize, bool                      flushingToDisk  ) : ___2337(___36) , m_varIter(new ItemSetIterator(*___36, ___36->___896() ? ___36->___890() : 0, vars)) , m_zoneIter(new ItemSetIterator(*___36, ___36->___896() ? ___36->___891() : 0, ___4671)) , ___2680(___36, maxIJKSubzoneSize, maxFESubzoneSize) , m_zoneVarMetadata(new ___4707(*___36, *m_varIter, *m_zoneIter, flushingToDisk  )) , m_flushingToDisk(flushingToDisk) { REQUIRE(VALID_REF(___36)); } DataSetWriter::~DataSetWriter() {} ___372 DataSetWriter::writeDataSet( FileWriterInterface& szpltFile, ___1392&        szpltZoneHeaderFileLocs) { REQUIRE(szpltFile.___2041()); if (!___2337->___896()) return ___4226; ___372 ___2039 = ___4226;
 #if defined _WIN32 && _MSC_FULL_VER < 190022816
unsigned int oldOutputFormat = _set_output_format(_TWO_DIGIT_EXPONENT);
 #endif
m_zoneIter->reset(); ___4636 const ___341 = m_zoneIter->baseItem(); while (___2039 && m_zoneIter->hasNext()) { ___4636 const ___904 = m_zoneIter->next(); if (!___2337->___4638(___904+1)) continue; try { std::vector<___372> ___4564; ___372 ___4499; getZoneSharing(___4564, ___4499, ___904, ___341, szpltFile.___844());
 #if defined OUTPUT_TIMES
uint64_t ___3687 = ___4709::___717();
 #endif
___4711 ___4710(___2680, *___2337); boost::shared_ptr<___4709> ___4708 = ___4710.___4708(*m_varIter, ___904, m_flushingToDisk ? 0 : ___341, ___4564, ___4499);
 #if defined OUTPUT_TIMES
uint64_t ___1167 = ___4709::___717(); ___1931(NULL, "%g seconds partitioning zone.", (double)(___1167 - ___3687) / 1000.0);
 #endif
___2039 = ___4708->writeZone(szpltFile, szpltFile.fileLoc()); ___4636 const fileZone = ___904 - ___341; if (___2039) { szpltZoneHeaderFileLocs[fileZone] = ___4708->getZoneHeaderFilePosition();
 #if defined DO_SUBZONE_HISTOGRAM || defined DO_ITEMANDSUBZONE_HISTOGRAM
{ if (___4644(*___2337, ___904) && !zoneIsPartitioned(*___2337, ___904)) { boost::shared_ptr<___1350> zoneInfo = ___2680.getFEZoneInfo(___904, ___341);
 #ifdef DO_SUBZONE_HISTOGRAM
extern ___372 OutputSubzoneHistograms( char const*       szpltFileName, ___37&       ___36, ___4636       zone, boost::shared_ptr<___1350 const> ___1349); OutputSubzoneHistograms(szpltFile.___1394().c_str(), *___2337, ___904, zoneInfo);
 #endif
 #ifdef DO_ITEMANDSUBZONE_HISTOGRAM
extern ___372 OutputItemAndSubzoneHistograms( char const*       szpltFileName, ___37&       ___36, ___4636       zone, boost::shared_ptr<___1350 const> ___1349); OutputItemAndSubzoneHistograms(szpltFile.___1394().c_str(), *___2337, ___904, zoneInfo);
 #endif
} }
 #endif
} m_varIter->reset(); ___4352 const baseVar = m_varIter->baseItem(); while (m_varIter->hasNext()) { ___4352 const datasetVar = m_varIter->next(); ___4352 const fileVar = datasetVar - baseVar; m_zoneVarMetadata->m_vzMinMaxes[fileVar][fileZone] = ___4708->varMinMax(datasetVar); } if (szpltFile.___844() == ___845 && !m_flushingToDisk) ___2680.remove(___904); } catch(std::exception const& e) { ___2039 = ___1186(e.what()); } if (!___2039) ___2680.clear(); }
 #if defined _WIN32 && _MSC_FULL_VER < 190022816
_set_output_format(oldOutputFormat);
 #endif
return ___2039; } void DataSetWriter::replaceDataSource( ___37* ___36, ___3501      vars, ___3501      ___4671) { REQUIRE(VALID_REF_OR_NULL(___36)); ___2337 = ___36; ___2680.replaceDataSource(___36); if (___36 == NULL) { m_varIter.reset(); m_zoneIter.reset(); m_zoneVarMetadata.reset(); } else { m_varIter.reset(new ItemSetIterator(*___2337, ___36->___896() ? ___36->___890() : 0, vars)); m_zoneIter.reset(new ItemSetIterator(*___2337, ___36->___896() ? ___36->___891() : 0, ___4671)); m_zoneVarMetadata.reset(new ___4707(*___36, *m_varIter, *m_zoneIter, m_flushingToDisk)); } } void DataSetWriter::clear( int32_t        numZonesToRetain, int32_t const* zonesToRetain) { REQUIRE(IMPLICATION(numZonesToRetain > 0, VALID_REF(zonesToRetain))); if (numZonesToRetain == 0) { ___2680.clear(); } else { boost::unordered_set<int32_t> retainZoneSet; for (int32_t i = 0; i < numZonesToRetain; ++i) retainZoneSet.insert(zonesToRetain[i]); m_zoneIter->reset(); while (m_zoneIter->hasNext()) { ___4636 const ___904 = m_zoneIter->next(); if (___2337->___4638(___904 + 1) && retainZoneSet.find(___904 + 1) == retainZoneSet.end()) ___2680.remove(___904); } } } void DataSetWriter::getZoneSharing( std::vector<___372>& ___4564, ___372&              ___4499, ___4636             ___904, ___4636             ___341, DataFileType_e          ___844) const { REQUIRE(0 <= ___904 && ___2337->___4638(___904 + 1)); REQUIRE(0 <= ___341 && ___341 <= ___904); REQUIRE(___4564.empty()); m_varIter->reset(); ___4352 const baseVar = m_varIter->baseItem(); ___4636 const fileZone = ___904 - ___341; while (m_varIter->hasNext()) { ___4352 const fileVar = m_varIter->next() - baseVar; ___4564.push_back( !m_zoneVarMetadata->m_vzIsPassive[fileVar][fileZone] && m_zoneVarMetadata->m_vzShareVarWithZone[fileVar][fileZone] == -1); } ___4499 = (___844 != ___848 && m_zoneVarMetadata->m_zoneShareConnectivityWithZone[fileZone] == -1); } }}
