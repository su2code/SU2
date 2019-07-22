#include "ClassicOrderedZoneWriter.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/assign.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "AltTecUtil.h"
#include "checkPercentDone.h"
#include "FieldData.h"
#include "fileStuff.h"
#include "ItemSetIterator.h"
#include "writeValueArray.h"
namespace tecplot { namespace ___3933 { ClassicOrderedZoneWriter::ClassicOrderedZoneWriter( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499, ___37&                   ___36) : ClassicZoneWriterAbstract(varIter, zone, ___341, ___4564, ___4499, ___36) , m_faceNeighborGenerator(___36) , m_faceNeighborWriter(m_faceNeighborGenerator, zone, ___341) {} ClassicOrderedZoneWriter::~ClassicOrderedZoneWriter() {} ___372 ClassicOrderedZoneWriter::writeZoneConnectivity(FileWriterInterface& szpltFile) { ___372 ___2039 = ___4226; m_zoneFileLocations.___2498 = ___330; if (m_writeConnectivity) { m_zoneFileLocations.___2663 = szpltFile.fileLoc(); ___2039 = m_faceNeighborWriter.write(szpltFile); } else { m_zoneFileLocations.___2663 = ___330; } return ___2039; } uint64_t ClassicOrderedZoneWriter::zoneConnectivityFileSize(bool ___2002) { if (m_writeConnectivity) return m_faceNeighborWriter.sizeInFile(___2002); else return 0; } }}
