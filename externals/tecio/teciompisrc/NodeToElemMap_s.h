 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <boost/scoped_array.hpp>
#include "ThirdPartyHeadersEnd.h"
#include "MASTER.h"
#include "GLOBAL.h"
#include "CodeContract.h"
#include "NodeMap_s.h"
struct ___2743 { tecplot::___3933::___2718 m_nodeCount; boost::scoped_array<tecplot::___3933::___465> m_elemIndex; boost::scoped_array<tecplot::___3933::___465> m_elem; ___2743(___2730 const& ___2723, tecplot::___3933::___2718 nodeCount); tecplot::___3933::___465 cellCountForNode(tecplot::___3933::___2718 ___2709) { REQUIRE(0 <= ___2709 && ___2709 < m_nodeCount); return m_elemIndex[___2709 + 1] - m_elemIndex[___2709]; } };
