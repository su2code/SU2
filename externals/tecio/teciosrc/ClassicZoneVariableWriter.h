 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <string>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3933 { class ___37; class FileWriterInterface; class ___1352; class ItemSetIterator; class ClassicZoneVariableWriter { public: ClassicZoneVariableWriter( ItemSetIterator& varIter, ___4636      zone, ___4636      ___341, ___37&      ___36); static uint64_t varHeaderSizeInFile(bool ___2002); uint64_t varSizeInFile(___4352 ___4336, bool ___2002) const; ___372 writeVarHeader( FileWriterInterface& file, ValueLocation_e      ___4326, ___4352           ___4336); ___372 writeVar( FileWriterInterface& szpltFile, ___1352 const&     ___1351); private: ItemSetIterator&  m_varIter; ___4636 const ___2677; ___4636 const m_baseZone; ___37&       ___2337; std::string const m_zoneNumberLabel; }; }}
