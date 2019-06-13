 #pragma once
#include "ThirdPartyHeadersBegin.h"
#include <vector>
#include <boost/shared_ptr.hpp>
#include <vector>
#include "ThirdPartyHeadersEnd.h"
#include "basicTypes.h"
namespace tecplot { namespace ___3933 { class ___37; class ZoneInfoCache; class ___4709; class ItemSetIterator; class ___4711 { public: ___4711( ZoneInfoCache& zoneInfoCache, ___37& ___36); boost::shared_ptr<___4709> ___4708( ItemSetIterator&              varIter, ___4636                   zone, ___4636                   ___341, std::vector<___372> const& ___4564, ___372                     ___4499); protected: ZoneInfoCache& ___2680; ___37& ___2337; }; }}
