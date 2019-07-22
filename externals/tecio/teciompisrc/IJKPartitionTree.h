 #pragma once
 #if defined (___1989)
 #undef ___1989 
 #endif
#include "ThirdPartyHeadersBegin.h"
#include <iterator>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "ThirdPartyHeadersEnd.h"
 #if defined (___1989)
 #undef ___1989
 #endif
 #define ___1989 (-1)
#include "IJK.h"
#include "ItemAddress.h"
namespace tecplot { namespace ___3933 { int const RTREE_PAGE_SIZE = 8; typedef boost::geometry::model::point<int64_t, 3, boost::geometry::cs::cartesian> ___1853; typedef boost::geometry::model::box<___1853> ___1855; typedef std::pair<___1855, tecplot::___2090::___2980> ___1864; class ___1863 : public boost::geometry::index::rtree<___1864, boost::geometry::index::quadratic<RTREE_PAGE_SIZE> > { public: ___1863() {} explicit ___1863(std::vector<___1864> const& ___2981) : boost::geometry::index::rtree<___1864, boost::geometry::index::quadratic<RTREE_PAGE_SIZE> >(___2981) {} void ___13(___2090::___2980 ___2977, ___1844 const& ___2474, ___1844 const& ___2364) { ___1853 ___2478(___2474.i(), ___2474.___2105(), ___2474.___2134()); ___1853 ___2372(___2364.i(), ___2364.___2105(), ___2364.___2134()); insert(std::make_pair(___1855(___2478, ___2372), ___2977)); } bool find(___1844 const& ___2474, ___1844 const& ___2364, std::vector<___1864>& ___2099) { REQUIRE(___2099.empty()); ___1853 ___2478(___2474.i(), ___2474.___2105(), ___2474.___2134()); ___1853 ___2372(___2364.i(), ___2364.___2105(), ___2364.___2134()); query(boost::geometry::index::intersects(___1855(___2478, ___2372)), std::back_inserter(___2099)); return !___2099.empty(); } }; }}
