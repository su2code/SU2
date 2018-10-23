 #pragma once
#include "MASTER.h"
#include "GLOBAL.h"
#include "ThirdPartyHeadersBegin.h"
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "fileio.h"
struct AuxData_s { typedef boost::shared_ptr<AuxData_s> Ptr; struct AuxDataItem { AuxDataItem(AuxDataType_e auxDataType, ___372 ___3361, std::string const& ___2685, std::string const& ___4314) : m_auxDataType(auxDataType) , m_retain(___3361) , ___2495(___2685) , ___2667(___4314) {} AuxDataType_e m_auxDataType; ___372 m_retain; std::string ___2495; std::string ___2667; AuxDataItem() {} void writeToFile(std::ofstream& outputFile, bool ___4480) const { tecplot::tecioszl::writeScalar(outputFile, m_auxDataType, ___4480); tecplot::tecioszl::writeScalar(outputFile, m_retain, ___4480); tecplot::tecioszl::___4544(outputFile, ___2495, ___4480); tecplot::tecioszl::___4544(outputFile, ___2667, ___4480); } AuxDataItem(std::ifstream& inputFile, bool readASCII) { tecplot::tecioszl::readScalar(inputFile, (unsigned int&)m_auxDataType, readASCII); tecplot::tecioszl::readScalar(inputFile, m_retain, readASCII); tecplot::tecioszl::readString(inputFile, ___2495, readASCII); tecplot::tecioszl::readString(inputFile, ___2667, readASCII); } }; std::vector<AuxDataItem> m_auxDataItems; AuxData_s() {} void writeToFile(std::ofstream& outputFile, bool ___4480) const { tecplot::tecioszl::writeVectorOfObjects(outputFile, m_auxDataItems, ___4480); } static Ptr makePtr(std::ifstream& inputFile, bool readASCII) { Ptr ___3358 = boost::make_shared<AuxData_s>(); tecplot::tecioszl::readVectorOfObjects(inputFile, ___3358->m_auxDataItems, readASCII); return ___3358; } };
