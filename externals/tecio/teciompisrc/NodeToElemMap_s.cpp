#include "NodeToElemMap_s.h"
#include "ThirdPartyHeadersBegin.h"
 #if !defined TECIOMPI
#include <boost/atomic.hpp>
 #endif
#include <boost/ref.hpp>
#include <boost/shared_ptr.hpp>
#include <algorithm>
#include "ThirdPartyHeadersEnd.h"
 #if !defined TECIOMPI
#include "JobControl_s.h"
 #endif
using namespace tecplot::___3933;
 #define MIN_CELLS_FOR_MULTITHREAD 100000 
 #if !defined TECIOMPI
namespace { struct NodeToElemData { ___2730 const& ___2723; boost::scoped_array<tecplot::___3933::___465> const& elemIndex; boost::scoped_array<tecplot::___3933::___465>& elem; boost::scoped_array<boost::atomic<int> >& count; ___465 begin; ___465 end; NodeToElemData( ___2730 const& ___2723, boost::scoped_array<tecplot::___3933::___465> const& elemIndex, boost::scoped_array<tecplot::___3933::___465>& elem, boost::scoped_array<boost::atomic<int> >& count, ___465 begin, ___465 end) : ___2723(___2723) , elemIndex(elemIndex) , elem(elem) , count(count) , begin(begin) , end(end) {} }; void STDCALL fillElemArrayForCellRange(___90 ___2124) { NodeToElemData* nodeToElemData = reinterpret_cast<NodeToElemData*>(___2124); for(___465 ___449 = nodeToElemData->begin; ___449 < nodeToElemData->end; ++___449) { for(int32_t ___681 = 0; ___681 < nodeToElemData->___2723.___2500; ++___681) { ___2718 ___2709 = nodeToElemData->___2723.___4314(___449 * nodeToElemData->___2723.___2500 + ___681); ___465 ind = nodeToElemData->elemIndex[___2709] + nodeToElemData->count[___2709]++; ___478(nodeToElemData->elem[ind] == 0); nodeToElemData->elem[ind] = ___449; } } } }
 #endif
___2743::___2743(___2730 const& ___2723, ___2718 nodeCount) : m_nodeCount(nodeCount) { size_t indexSize = static_cast<size_t>(nodeCount + 1); m_elemIndex.reset(new ___465[indexSize]); memset(&m_elemIndex[0], 0, indexSize * sizeof(m_elemIndex[0])); size_t arraySize = ___2723.___2500 * static_cast<size_t>(___2723.___2392); m_elem.reset(new ___465[arraySize]); memset(&m_elem[0], 0, arraySize * sizeof(m_elem[0])); boost::scoped_array<int> count(new int[m_nodeCount]); memset(&count[0], 0, m_nodeCount * sizeof(count[0])); for(size_t i = 0; i < arraySize; ++i) { int64_t ___2709 = ___2723.___4314(i); ___478(0 <= ___2709 && ___2709 < m_nodeCount); ++count[___2709]; } m_elemIndex[0] = 0; for(___2718 ___2709 = 0; ___2709 < m_nodeCount; ++___2709) { m_elemIndex[___2709 + 1] = m_elemIndex[___2709] + count[___2709]; count[___2709] = 0; }
 #if !defined TECIOMPI
int numThreads = 1; if (___2723.___2392 >= MIN_CELLS_FOR_MULTITHREAD) numThreads = ___2122::___2827(); if (numThreads == 1) {
 #endif
for(___465 ___449 = 0; ___449 < ___2723.___2392; ++___449) { for(int32_t ___681 = 0; ___681 < ___2723.___2500; ++___681) { ___2718 ___2709 = ___2723.___4314(___449 * ___2723.___2500 + ___681); ___465 ind = m_elemIndex[___2709] + count[___2709]++; ___478(m_elem[ind] == 0); m_elem[ind] = ___449; } }
 #if !defined TECIOMPI
} else { boost::scoped_array<boost::atomic<int> > atomiccount(new boost::atomic<int>[m_nodeCount]); for(___2718 i = 0; i < m_nodeCount; ++i) atomiccount[i] = 0; std::vector<boost::shared_ptr<NodeToElemData> > nodeToElemData; for(int i = 0; i < numThreads; ++i) { ___465 begin = static_cast<___465>((size_t)___2723.___2392 * i / numThreads); ___465 end = static_cast<___465>((size_t)___2723.___2392 * (i + 1) / numThreads); nodeToElemData.push_back(boost::make_shared<NodeToElemData>(boost::cref(___2723), boost::cref(m_elemIndex), boost::ref(m_elem), boost::ref(atomiccount), begin, end)); } ___2122 ___2119; for(int i = 0; i < numThreads; ++i) ___2119.addJob(fillElemArrayForCellRange, reinterpret_cast<___90>(nodeToElemData[i].get())); ___2119.wait(); }
 #endif
}
