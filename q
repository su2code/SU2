[1mdiff --git a/Common/include/geometry/CMultiGridGeometry.hpp b/Common/include/geometry/CMultiGridGeometry.hpp[m
[1mindex e1ba5a80fd..b5529ed4fe 100644[m
[1m--- a/Common/include/geometry/CMultiGridGeometry.hpp[m
[1m+++ b/Common/include/geometry/CMultiGridGeometry.hpp[m
[36m@@ -92,31 +92,13 @@[m [mclass CMultiGridGeometry final : public CGeometry {[m
    * \param[in] fine_grid - Fine grid geometry.[m
    * \param[in] config - Configuration.[m
    * \param[in] iMesh - Multigrid level being built, used to label the summary.[m
[32m+[m[32m   * \param[in] mixedBC - Nodes that must stay on their own, from FindMixedBoundaryNodes.[m
[32m+[m[32m   * \param[in] onPhysBoundary - Nodes carrying a physical boundary condition, excluding SEND_RECEIVE.[m
    * \return Summary of the paving, empty except on the master rank.[m
    */[m
[31m-  string AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid, const CConfig* config,[m
[31m-                                  unsigned short iMesh);[m
[31m-[m
[31m-  /*![m
[31m-   * \brief Weakest and strongest dual-grid coupling at each node, and the neighbour across the[m
[31m-   *        strongest edge. Their ratio is the local cell aspect ratio, available on every MG level.[m
[31m-   */[m
[31m-  struct CNodeStiffness {[m
[31m-    vector<su2double> wMin, wMax;    /*!< \brief Weakest and strongest edge coupling at each node. */[m
[31m-    vector<unsigned long> jStiffest; /*!< \brief Neighbour across the strongest edge. */[m
[31m-[m
[31m-    /*!< \brief Local aspect ratio at a node, 1.0 where it could not be measured. */[m
[31m-    su2double AspectRatio(unsigned long iPoint) const {[m
[31m-      return (wMin[iPoint] > 0.0) ? wMax[iPoint] / wMin[iPoint] : su2double(1.0);[m
[31m-    }[m
[31m-  };[m
[31m-[m
[31m-  /*![m
[31m-   * \brief Measure the dual-grid coupling at every node of a grid.[m
[31m-   * \param[in] fine_grid - Grid to measure.[m
[31m-   * \return Weakest/strongest coupling per node.[m
[31m-   */[m
[31m-  CNodeStiffness ComputeNodeStiffness(const CGeometry* fine_grid) const;[m
[32m+[m[32m  string PaveAdvancingFronts(unsigned long& Index_CoarseCV, const CGeometry* fine_grid, const CConfig* config,[m
[32m+[m[32m                                  unsigned short iMesh, const vector<char>& mixedBC,[m
[32m+[m[32m                                  const vector<char>& onPhysBoundary);[m
 [m
   /*![m
    * \brief Boundary nodes that seed a front, with the direction each starts marching in.[m
[36m@@ -124,6 +106,7 @@[m [mclass CMultiGridGeometry final : public CGeometry {[m
   struct CFrontSeeds {[m
     vector<unsigned long> node;                    /*!< \brief Seed node on the boundary. */[m
     vector<std::array<su2double, MAXNDIM>> normal; /*!< \brief Unit normal there, pointing into the domain. */[m
[32m+[m[32m    unsigned long nRefusedCurvature = 0; /*!< \brief Euler wall nodes the curvature limit kept out. */[m
   };[m
 [m
   /*![m
[36m@@ -131,10 +114,9 @@[m [mclass CMultiGridGeometry final : public CGeometry {[m
    *        boundary carrying a stretched layer normal to itself.[m
    * \param[in] fine_grid - Fine grid geometry.[m
    * \param[in] config - Definition of the particular problem.[m
[31m-   * \param[in] stiff - Node coupling from ComputeNodeStiffness.[m
    * \return Seed nodes and their inward boundary normals.[m
    */[m
[31m-  CFrontSeeds SeedFrontNodes(const CGeometry* fine_grid, const CConfig* config, const CNodeStiffness& stiff) const;[m
[32m+[m[32m  CFrontSeeds SeedFrontNodes(const CGeometry* fine_grid, const CConfig* config) const;[m
 [m
   /*![m
    * \brief Partition the seed nodes into compact surface patches by repeated pairwise matching. Each[m
[36m@@ -143,14 +125,19 @@[m [mclass CMultiGridGeometry final : public CGeometry {[m
    * \param[in] fine_grid - Fine grid geometry.[m
    * \param[in] config - Definition of the particular problem.[m
    * \param[in] mixedBC - Nodes that must stay on their own, from FindMixedBoundaryNodes.[m
[31m-   * \return One vector of indices into seeds.node per patch.[m
[32m+[m[32m   * \return One vector of indices into seeds.node per patch, at most two entries in 2D, four in 3D.[m
    */[m
   vector<vector<unsigned long>> BuildFrontPatches(const CFrontSeeds& seeds, const CGeometry* fine_grid,[m
                                                   const CConfig* config, const vector<char>& mixedBC) const;[m
 [m
[32m+[m[32m  string pavingReport; /*!< \brief Paving summary for this level. */[m
[32m+[m
  public:[m
[31m-  /*!< \brief Paving summary for this level. */[m
[31m-  string pavingReport;[m
[32m+[m[32m  /*![m
[32m+[m[32m   * \brief Get the paving summary for this level, for console output.[m
[32m+[m[32m   * \return Summary text, empty except on the master rank.[m
[32m+[m[32m   */[m
[32m+[m[32m  const string& GetPavingReport() const { return pavingReport; }[m
 [m
   /*--- This is to suppress Woverloaded-virtual, omitting it has no negative impact. ---*/[m
   using CGeometry::SetBoundControlVolume;[m
[1mdiff --git a/Common/include/geometry/dual_grid/CPoint.hpp b/Common/include/geometry/dual_grid/CPoint.hpp[m
[1mindex 73e6f8fb04..bac83ccbfd 100644[m
[1m--- a/Common/include/geometry/dual_grid/CPoint.hpp[m
[1m+++ b/Common/include/geometry/dual_grid/CPoint.hpp[m
[36m@@ -656,6 +656,15 @@[m [mclass CPoint {[m
     Children_CV[iPoint][nchildren_CV] = children_CV;[m
   }[m
 [m
[32m+[m[32m  /*![m
[32m+[m[32m   * \brief Set all the children control volumes of an agglomerated control volume at once.[m
[32m+[m[32m   * \param[in] iPoint - Index of the point.[m
[32m+[m[32m   * \param[in] children_CV - Indices of the children control volumes.[m
[32m+[m[32m   */[m
[32m+[m[32m  inline void SetChildren_CV(unsigned long iPoint, const vector<unsigned long>& children_CV) {[m
[32m+[m[32m    Children_CV[iPoint] = children_CV;[m
[32m+[m[32m  }[m
[32m+[m
   /*![m
    * \brief Get the parent control volume of an agglomerated control volume.[m
    * \param[in] iPoint - Index of the point.[m
[36m@@ -673,6 +682,14 @@[m [mclass CPoint {[m
     return Children_CV[iPoint][nchildren_CV];[m
   }[m
 [m
[32m+[m[32m  /*![m
[32m+[m[32m   * \brief Get the children control volumes of an agglomerated control volume. Only the first[m
[32m+[m[32m   *        GetnChildren_CV entries are meaningful, the storage may be longer.[m
[32m+[m[32m   * \param[in] iPoint - Index of the point.[m
[32m+[m[32m   * \return Indices of the children control volumes.[m
[32m+[m[32m   */[m
[32m+[m[32m  inline const vector<unsigned long>& GetChildren_CV(unsigned long iPoint) const { return Children_CV[iPoint]; }[m
[32m+[m
   /*![m
    * \brief Get information about if a control volume has been agglomerated.[m
    * \param[in] iPoint - Index of the point.[m
[1mdiff --git a/Common/src/CConfig.cpp b/Common/src/CConfig.cpp[m
[1mindex 29161a7433..6a0469295e 100644[m
[1m--- a/Common/src/CConfig.cpp[m
[1m+++ b/Common/src/CConfig.cpp[m
[36m@@ -2069,7 +2069,8 @@[m [mvoid CConfig::SetConfig_Options() {[m
   /*!\brief MG_MIN_MESHSIZE\n DESCRIPTION: Minimum number of CVs on the coarsest multigrid level, checked per MPI rank[m
    * (i.e. on the smallest partition). Levels that would produce fewer CVs on any rank are not created. DEFAULT: 500 \ingroup Config*/[m
   addUnsignedLongOption("MG_MIN_MESHSIZE", MGOptions.MG_Min_MeshSize, 500);[m
[31m-  /*!\brief MG_IMPLICIT_LINES\n DESCRIPTION: Enable agglomeration along implicit lines from wall seeds. DEFAULT: NO \ingroup Config*/[m
[32m+[m[32m  /*!\brief MG_IMPLICIT_LINES\n DESCRIPTION: Pave the coarse grid with advancing fronts raised from boundaries[m
[32m+[m[32m   * that carry a stretched layer normal to themselves. DEFAULT: NO \ingroup Config*/[m
   addBoolOption("MG_IMPLICIT_LINES", MGOptions.MG_Implicit_Lines, false);[m
   /*!\brief MG_STARTUP_ITER\n DESCRIPTION: Max number of iterations spent on each mesh during the Full[m
    * Multigrid (FMG) startup phase. DEFAULT: 100 \ingroup Config*/[m
[1mdiff --git a/Common/src/geometry/CMultiGridGeometry.cpp b/Common/src/geometry/CMultiGridGeometry.cpp[m
[1mindex 17ef93a31d..177add7c96 100644[m
[1m--- a/Common/src/geometry/CMultiGridGeometry.cpp[m
[1m+++ b/Common/src/geometry/CMultiGridGeometry.cpp[m
[36m@@ -30,6 +30,13 @@[m
 #include "../../include/toolboxes/printing_toolbox.hpp"[m
 #include "../../../Common/include/toolboxes/geometry_toolbox.hpp"[m
 [m
[32m+[m[32m#include <algorithm>[m
[32m+[m[32m#include <array>[m
[32m+[m[32m#include <limits>[m
[32m+[m[32m#include <map>[m
[32m+[m[32m#include <sstream>[m
[32m+[m[32m#include <utility>[m
[32m+[m
 namespace {[m
 [m
 /*--- Euler wall nodes are not agglomerated where the surface turns by more than this, in[m
[36m@@ -118,14 +125,6 @@[m [mCMultiGridGeometry::CMultiGridGeometry(CGeometry* fine_grid, CConfig* config, un[m
     }[m
   }[m
 [m
[31m-  /*--- STEP 0: pave the domain with advancing fronts rising from the boundaries. The coarse CVs it[m
[31m-   *    creates occupy [firstLineCV, endLineCV). ---*/[m
[31m-  const auto firstLineCV = Index_CoarseCV;[m
[31m-  if (config->GetMGOptions().MG_Implicit_Lines) {[m
[31m-    pavingReport = AgglomerateImplicitLines(Index_CoarseCV, fine_grid, config, iMesh);[m
[31m-  }[m
[31m-  const auto endLineCV = Index_CoarseCV;[m
[31m-[m
   /*--- Points carrying a physical boundary condition. This does not include SEND_RECEIVE. ---*/[m
   vector<char> onPhysBoundary(fine_grid->GetnPoint(), 0);[m
   for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {[m
[36m@@ -136,6 +135,14 @@[m [mCMultiGridGeometry::CMultiGridGeometry(CGeometry* fine_grid, CConfig* config, un[m
   /*--- Nodes where two different boundary conditions meet.  ---*/[m
   const auto mixedBC = FindMixedBoundaryNodes(fine_grid, config);[m
 [m
[32m+[m[32m  /*--- STEP 0: pave the domain with advancing fronts rising from the boundaries. The coarse CVs it[m
[32m+[m[32m   *    creates occupy [firstLineCV, endLineCV). ---*/[m
[32m+[m[32m  const auto firstLineCV = Index_CoarseCV;[m
[32m+[m[32m  if (config->GetMGOptions().MG_Implicit_Lines) {[m
[32m+[m[32m    pavingReport = PaveAdvancingFronts(Index_CoarseCV, fine_grid, config, iMesh, mixedBC, onPhysBoundary);[m
[32m+[m[32m  }[m
[32m+[m[32m  const auto endLineCV = Index_CoarseCV;[m
[32m+[m
   /*--- STEP 1: The first step is the boundary agglomeration. ---*/[m
   for (auto iMarker = 0u; iMarker < fine_grid->GetnMarker(); iMarker++) {[m
     /*--- Skip periodic boundaries: do not agglomerate on periodic markers. ---*/[m
[36m@@ -721,10 +728,8 @@[m [mCMultiGridGeometry::CMultiGridGeometry(CGeometry* fine_grid, CConfig* config, un[m
       for (auto iCoarsePoint = 0ul; iCoarsePoint < nPointDomain; iCoarsePoint++) {[m
         const auto iNew = newIndex[iCoarsePoint];[m
         if ((iNew == NO_INDEX) || (iNew == iCoarsePoint)) continue;[m
[31m-        const auto nChildren = nodes->GetnChildren_CV(iCoarsePoint);[m
[31m-        for (auto iChildren = 0u; iChildren < nChildren; iChildren++)[m
[31m-          nodes->SetChildren_CV(iNew, iChildren, nodes->GetChildren_CV(iCoarsePoint, iChildren));[m
[31m-        nodes->SetnChildren_CV(iNew, nChildren);[m
[32m+[m[32m        nodes->SetChildren_CV(iNew, nodes->GetChildren_CV(iCoarsePoint));[m
[32m+[m[32m        nodes->SetnChildren_CV(iNew, nodes->GetnChildren_CV(iCoarsePoint));[m
         nodes->SetAgglomerate_Indirect(iNew, nodes->GetAgglomerate_Indirect(iCoarsePoint));[m
       }[m
       for (auto iCoarsePoint = nKept; iCoarsePoint < nPointDomain; iCoarsePoint++)[m
[36m@@ -1539,43 +1544,6 @@[m [msu2double CMultiGridGeometry::ComputeLocalCurvature(const CGeometry* fine_grid,[m
   return max_angle;[m
 }[m
 [m
[31m-CMultiGridGeometry::CNodeStiffness CMultiGridGeometry::ComputeNodeStiffness(const CGeometry* fine_grid) const {[m
[31m-  /*--- Coupling across the dual face between a node and a neighbour, so the ratio of largest to[m
[31m-   *    smallest weight at a node is the local aspect ratio. ---*/[m
[31m-  const auto nPointFine = fine_grid->GetnPoint();[m
[31m-[m
[31m-  CNodeStiffness stiff;[m
[31m-  stiff.wMin.assign(nPointFine, 0.0);[m
[31m-  stiff.wMax.assign(nPointFine, 0.0);[m
[31m-  stiff.jStiffest.assign(nPointFine, std::numeric_limits<unsigned long>::max());[m
[31m-[m
[31m-  for (auto iPoint = 0ul; iPoint < nPointFine; ++iPoint) {[m
[31m-    su2double wmin = std::numeric_limits<su2double>::max(), wmax = 0.0;[m
[31m-    auto jStiffest = std::numeric_limits<unsigned long>::max();[m
[31m-[m
[31m-    for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(iPoint); ++iNeigh) {[m
[31m-      const auto jPoint = fine_grid->nodes->GetPoint(iPoint, iNeigh);[m
[31m-      const auto iEdge = fine_grid->nodes->GetEdge(iPoint, iNeigh);[m
[31m-      const su2double area = GeometryToolbox::Norm(nDim, fine_grid->edges->GetNormal(iEdge));[m
[31m-      const su2double w =[m
[31m-          0.5 * area * (1.0 / fine_grid->nodes->GetVolume(iPoint) + 1.0 / fine_grid->nodes->GetVolume(jPoint));[m
[31m-      if (w > wmax) {[m
[31m-        wmax = w;[m
[31m-        jStiffest = jPoint;[m
[31m-      }[m
[31m-      wmin = std::min(wmin, w);[m
[31m-    }[m
[31m-[m
[31m-    /*--- A node with no neighbours keeps the zeroed defaults, so AspectRatio reads 1. ---*/[m
[31m-    if (jStiffest != std::numeric_limits<unsigned long>::max()) {[m
[31m-      stiff.wMin[iPoint] = wmin;[m
[31m-      stiff.wMax[iPoint] = wmax;[m
[31m-      stiff.jStiffest[iPoint] = jStiffest;[m
[31m-    }[m
[31m-  }[m
[31m-  return stiff;[m
[31m-}[m
[31m-[m
 namespace {[m
 [m
 /*--- Unit normal of a boundary at a vertex, false if the marker does not reach iPoint. Boundary[m
[36m@@ -1669,8 +1637,8 @@[m [munsigned long TagOfSet(const CGeometry* grid, const vector<unsigned long>& set)[m
 [m
 }  // namespace[m
 [m
[31m-CMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeometry* fine_grid, const CConfig* config,[m
[31m-                                                                   const CNodeStiffness& stiff) const {[m
[32m+[m[32mCMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeometry* fine_grid,[m
[32m+[m[32m                                                                   const CConfig* config) const {[m
   /*--- Fraction of a marker's nodes that must sit in a layer before the whole marker may seed. ---*/[m
   constexpr passivedouble QUALIFIED_FRACTION = 0.5;[m
   constexpr passivedouble ANGLE_THRESHOLD_DEG = 30.0;[m
[36m@@ -1691,11 +1659,28 @@[m [mCMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeomet[m
 [m
 [m
   /*--- True if the mesh at iPoint is stretched along the boundary normal, i.e. this boundary has a[m
[31m-   *    layer growing off it the way a viscous wall does. ---*/[m
[32m+[m[32m   *    layer growing off it the way a viscous wall does. The coupling across the dual face between[m
[32m+[m[32m   *    iPoint and a neighbour is measured here, so the ratio of largest to smallest weight is the[m
[32m+[m[32m   *    local aspect ratio. Only boundary nodes are ever asked, at most twice each. ---*/[m
   auto hasLayerNormalTo = [&](unsigned long iPoint, const su2double* unitNormal) {[m
[31m-    const auto jStiffest = stiff.jStiffest[iPoint];[m
[32m+[m[32m    su2double wMin = std::numeric_limits<su2double>::max(), wMax = 0.0;[m
[32m+[m[32m    auto jStiffest = NO_POINT;[m
[32m+[m[32m    for (auto iNeigh = 0u; iNeigh < fine_grid->nodes->GetnPoint(iPoint); ++iNeigh) {[m
[32m+[m[32m      const auto jPoint = fine_grid->nodes->GetPoint(iPoint, iNeigh);[m
[32m+[m[32m      const auto iEdge = fine_grid->nodes->GetEdge(iPoint, iNeigh);[m
[32m+[m[32m      const su2double area = GeometryToolbox::Norm(nDim, fine_grid->edges->GetNormal(iEdge));[m
[32m+[m[32m      const su2double w =[m
[32m+[m[32m          0.5 * area * (1.0 / fine_grid->nodes->GetVolume(iPoint) + 1.0 / fine_grid->nodes->GetVolume(jPoint));[m
[32m+[m[32m      if (w > wMax) {[m
[32m+[m[32m        wMax = w;[m
[32m+[m[32m        jStiffest = jPoint;[m
[32m+[m[32m      }[m
[32m+[m[32m      wMin = std::min(wMin, w);[m
[32m+[m[32m    }[m
[32m+[m
[32m+[m[32m    /*--- A node with no neighbours has no aspect ratio to measure. ---*/[m
     if (jStiffest == NO_POINT) return false;[m
[31m-    if (stiff.AspectRatio(iPoint) < MIN_AR) return false;[m
[32m+[m[32m    if (((wMin > 0.0) ? wMax / wMin : su2double(1.0)) < MIN_AR) return false;[m
 [m
     su2double vec[MAXNDIM] = {0.0};[m
     GeometryToolbox::Distance(nDim, fine_grid->nodes->GetCoord(jStiffest), fine_grid->nodes->GetCoord(iPoint), vec);[m
[36m@@ -1719,8 +1704,10 @@[m [mCMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeomet[m
       /*--- A front claims its seed before the boundary pass runs, so the Euler wall curvature[m
        *    limit is applied here too, on the same terms. ---*/[m
       if ((config->GetMarker_All_KindBC(iMarker) == EULER_WALL) && !fine_grid->boundIsStraight[iMarker] &&[m
[31m-          (ComputeLocalCurvature(fine_grid, iPoint, iMarker) >= EULER_WALL_MAX_CURVATURE))[m
[32m+[m[32m          (ComputeLocalCurvature(fine_grid, iPoint, iMarker) >= EULER_WALL_MAX_CURVATURE)) {[m
[32m+[m[32m        seeds.nRefusedCurvature++;[m
         continue;[m
[32m+[m[32m      }[m
 [m
       std::array<su2double, MAXNDIM> n0{};[m
       for (unsigned short d = 0; d < nDim; ++d) n0[d] = Normal[d];[m
[36m@@ -1788,8 +1775,10 @@[m [mCMultiGridGeometry::CFrontSeeds CMultiGridGeometry::SeedFrontNodes(const CGeomet[m
 vector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFrontSeeds& seeds,[m
                                                                     const CGeometry* fine_grid, const CConfig* config,[m
                                                                     const vector<char>& mixedBC) const {[m
[31m-  /*--- Repeated pairwise matching, one round per doubling, partitions the seeds into compact[m
[31m-   *    patches: a boundary edge in 2D, a boundary quadrilateral in 3D. ---*/[m
[32m+[m[32m  /*--- Repeated pairwise matching groups the seeds into connected patches of 1 to max_group seeds,[m
[32m+[m[32m   *    smaller wherever no partner was found. A patch is any connected shape the boundary gives:[m
[32m+[m[32m   *    a square, a strip, a triangle or a single node. What follows keys off the patch size and[m
[32m+[m[32m   *    the layer isomorphism check, never off an assumed shape. ---*/[m
   const auto nSeeds = seeds.node.size();[m
   const unsigned long max_group = (nDim == 2) ? 2 : 4;[m
 [m
[36m@@ -1828,13 +1817,14 @@[m [mvector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFront[m
     groups.push_back({si});[m
   }[m
 [m
[31m-  const unsigned nRounds = (max_group <= 1) ? 0 : ((max_group <= 2) ? 1 : 2);[m
[32m+[m[32m  const unsigned nRounds = (max_group <= 2) ? 1 : 2;[m
 [m
   /*--- One admissible merge of two groups, weighted by how many seed-to-seed adjacencies they[m
    *    share: 2 for a group lying alongside, 1 for one continuing in the same direction. ---*/[m
   struct CMerge {[m
     unsigned long g, h;       /*!< \brief The two groups, g < h. */[m
     unsigned long weight;     /*!< \brief Shared adjacencies: 2 makes a square, 1 makes a strip. */[m
[32m+[m[32m    unsigned long size;       /*!< \brief Seeds the merged group would hold. */[m
     unsigned long keyG, keyH; /*!< \brief Their global-index keys, the deterministic tie-break. */[m
   };[m
 [m
[36m@@ -1869,15 +1859,17 @@[m [mvector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFront[m
           if (nShared[h]++ == 0) touched.push_back(h);[m
         }[m
       for (auto h : touched) {[m
[31m-        merges.push_back({g, h, nShared[h], gkey[g], gkey[h]});[m
[32m+[m[32m        merges.push_back({g, h, nShared[h], static_cast<unsigned long>(groups[g].size() + groups[h].size()),[m
[32m+[m[32m                          gkey[g], gkey[h]});[m
         nShared[h] = 0;[m
       }[m
     }[m
 [m
     /*--- Best merges first over all groups at once, so every square is considered before the first[m
[31m-     *    strip. ---*/[m
[32m+[m[32m     *    strip. At equal weight the bigger patch wins, to fill max_group before settling for less. ---*/[m
     std::sort(merges.begin(), merges.end(), [](const CMerge& a, const CMerge& b) {[m
       if (a.weight != b.weight) return a.weight > b.weight;[m
[32m+[m[32m      if (a.size != b.size) return a.size > b.size;[m
       if (a.keyG != b.keyG) return a.keyG < b.keyG;[m
       return a.keyH < b.keyH;[m
     });[m
[36m@@ -1905,13 +1897,14 @@[m [mvector<vector<unsigned long>> CMultiGridGeometry::BuildFrontPatches(const CFront[m
   return groups;[m
 }[m
 [m
[31m-string CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,[m
[31m-                                                    const CConfig* config, unsigned short iMesh) {[m
[32m+[m[32mstring CMultiGridGeometry::PaveAdvancingFronts(unsigned long& Index_CoarseCV, const CGeometry* fine_grid,[m
[32m+[m[32m                                                    const CConfig* config, unsigned short iMesh,[m
[32m+[m[32m                                                    const vector<char>& mixedBC,[m
[32m+[m[32m                                                    const vector<char>& onPhysBoundary) {[m
   /*--- Paving by advancing fronts. Each boundary patch rises into the domain keeping its footprint,[m
    *    stopping at a boundary or where the next layer is not isomorphic to the current one. ---*/[m
   const auto starting_Index_CoarseCV = Index_CoarseCV;[m
   const auto nPointFine = fine_grid->GetnPoint();[m
[31m-  const auto nMarkerFine = fine_grid->GetnMarker();[m
   constexpr auto NO_POINT = std::numeric_limits<unsigned long>::max();[m
   const short int maxAgglomSize = (nDim == 2) ? 4 : 8;[m
 [m
[36m@@ -1922,22 +1915,10 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
   /*--- Weight of the new step direction when the front's marching direction is updated. ---*/[m
   constexpr passivedouble DIR_BLEND = 0.5;[m
 [m
[31m-  const auto stiff = ComputeNodeStiffness(fine_grid);[m
[31m-[m
   /*--- PHASE 1. SeedFrontNodes is collective and must be reached by every rank. ---*/[m
[31m-  const auto seeds = SeedFrontNodes(fine_grid, config, stiff);[m
[31m-  const auto mixedBC = FindMixedBoundaryNodes(fine_grid, config);[m
[32m+[m[32m  const auto seeds = SeedFrontNodes(fine_grid, config);[m
   const auto patches = BuildFrontPatches(seeds, fine_grid, config, mixedBC);[m
 [m
[31m-  /*--- Nodes on a boundary carrying a boundary condition, which a front must not grow into. CPoint's[m
[31m-   *    Boundary flag cannot be used, it is also set by SEND_RECEIVE. ---*/[m
[31m-  vector<char> onPhysicalBoundary(nPointFine, 0);[m
[31m-  for (auto iMarker = 0u; iMarker < nMarkerFine; iMarker++) {[m
[31m-    if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) continue;[m
[31m-    for (auto iVertex = 0ul; iVertex < fine_grid->GetnVertex(iMarker); iVertex++)[m
[31m-      onPhysicalBoundary[fine_grid->vertex[iMarker][iVertex]->GetNode()] = 1;[m
[31m-  }[m
[31m-[m
   /*==================================================================================================[m
    *  PHASE 2 - advance every front, one layer per round.[m
    *================================================================================================*/[m
[36m@@ -2007,7 +1988,7 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
 [m
       const su2double dot = GeometryToolbox::DotProduct(nDim, vec, marchDir);[m
       const bool admissible =[m
[31m-          !(onPhysicalBoundary[jPoint] && EntersBoundary(fine_grid, config, nDim, jPoint, vec, cos_boundary)) &&[m
[32m+[m[32m          !(onPhysBoundary[jPoint] && EntersBoundary(fine_grid, config, nDim, jPoint, vec, cos_boundary)) &&[m
           GeometricalCheck(jPoint, fine_grid, config);[m
 [m
       /*--- Halo parents are assigned by the owning rank through the MPI relay, so a halo node is[m
[36m@@ -2041,9 +2022,32 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
   vector<unsigned long> newLayer;[m
 [m
   /*--- Summed over all ranks for the one-line report at the end. ---*/[m
[31m-  enum { P_STACKS, P_LAYERS, P_COVERED, P_COUNT };[m
[32m+[m[32m  enum {[m
[32m+[m[32m    P_STACKS,       /*!< \brief Fronts started. */[m
[32m+[m[32m    P_LAYERS,       /*!< \brief Layers laid by all fronts. */[m
[32m+[m[32m    P_COVERED,      /*!< \brief Fine nodes given a coarse parent by the paving. */[m
[32m+[m[32m    P_RETIRED,      /*!< \brief Fronts that died on a failed layer rather than running out. */[m
[32m+[m[32m    P_SEED_CURV,    /*!< \brief Euler wall seeds refused by the curvature limit. */[m
[32m+[m[32m    P_HAND_SENT,    /*!< \brief Footprints offered across a partition interface. */[m
[32m+[m[32m    P_HAND_SPLIT,   /*!< \brief ...of those, ones whose nodes are owned by more than one rank. */[m
[32m+[m[32m    P_HAND_CONTEST, /*!< \brief Halo nodes two fronts reached for in the same round. */[m
[32m+[m[32m    P_HAND_TAKEN,   /*!< \brief Inherited footprints adopted. */[m
[32m+[m[32m    P_HAND_DROPPED, /*!< \brief Inherited footprints refused as not free or not connected. */[m
[32m+[m[32m    P_PATCH1,       /*!< \brief Seed patches by width, 1 to 4. ---*/[m
[32m+[m[32m    P_PATCH2,[m
[32m+[m[32m    P_PATCH3,[m
[32m+[m[32m    P_PATCH4,[m
[32m+[m[32m    P_COUNT[m
[32m+[m[32m  };[m
   unsigned long ct[P_COUNT] = {0};[m
 [m
[32m+[m[32m  /*--- Patch widths, the footprint every front starts from. ---*/[m
[32m+[m[32m  for (const auto& patch : patches) {[m
[32m+[m[32m    if (patch.empty()) continue;[m
[32m+[m[32m    ct[P_PATCH1 + std::min<size_t>(patch.size(), 4) - 1]++;[m
[32m+[m[32m  }[m
[32m+[m[32m  ct[P_SEED_CURV] = seeds.nRefusedCurvature;[m
[32m+[m
   /*--- One footprint node arriving from a neighbouring rank, to be regrouped by tag. ---*/[m
   struct CInherited {[m
     unsigned long tag;      /*!< \brief The front it belongs to, the same name on both ranks. */[m
[36m@@ -2072,10 +2076,9 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
   /*--- Turn everything buffered for this front into one coarse control volume. ---*/[m
   auto emit = [&](unsigned long f) {[m
     if (fronts[f].pending.empty()) return;[m
[31m-    for (unsigned long c = 0; c < fronts[f].pending.size(); ++c) {[m
[31m-      const auto p = fronts[f].pending[c];[m
[32m+[m[32m    nodes->SetChildren_CV(Index_CoarseCV, fronts[f].pending);[m
[32m+[m[32m    for (auto p : fronts[f].pending) {[m
       fine_grid->nodes->SetParent_CV(p, Index_CoarseCV);[m
[31m-      nodes->SetChildren_CV(Index_CoarseCV, c, p);[m
       if (fine_grid->nodes->GetAgglomerate_Indirect(p)) nodes->SetAgglomerate_Indirect(Index_CoarseCV, true);[m
     }[m
     nodes->SetnChildren_CV(Index_CoarseCV, static_cast<unsigned short>(fronts[f].pending.size()));[m
[36m@@ -2122,6 +2125,44 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
     emit(f);[m
   }[m
 [m
[32m+[m[32m  /*--- The handover carries a marching direction, but only to make discrete choices, so it is[m
[32m+[m[32m   *    exchanged as a passive type and never taped. ---*/[m
[32m+[m[32m  using CPassiveMPI = SelectMPIWrapper<passivedouble>::W;[m
[32m+[m
[32m+[m[32m  /*--- One SEND_RECEIVE marker pair per neighbour, with where its vertices sit in the flat[m
[32m+[m[32m   *    exchange buffers. ---*/[m
[32m+[m[32m  struct CHandoverPair {[m
[32m+[m[32m    unsigned short markerS, markerR;[m
[32m+[m[32m    int send_to, receive_from;[m
[32m+[m[32m    unsigned long nVertexS, nVertexR, offS, offR;[m
[32m+[m[32m  };[m
[32m+[m[32m  vector<CHandoverPair> handPairs;[m
[32m+[m[32m  unsigned long nSendTotal = 0, nRecvTotal = 0;[m
[32m+[m[32m  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {[m
[32m+[m[32m    if (!((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)))[m
[32m+[m[32m      continue;[m
[32m+[m[32m    CHandoverPair hp;[m
[32m+[m[32m    hp.markerS = iMarker;[m
[32m+[m[32m    hp.markerR = iMarker + 1;[m
[32m+[m[32m    hp.send_to = config->GetMarker_All_SendRecv(hp.markerS) - 1;[m
[32m+[m[32m    hp.receive_from = abs(config->GetMarker_All_SendRecv(hp.markerR)) - 1;[m
[32m+[m[32m    hp.nVertexS = fine_grid->nVertex[hp.markerS];[m
[32m+[m[32m    hp.nVertexR = fine_grid->nVertex[hp.markerR];[m
[32m+[m[32m    hp.offS = nSendTotal;[m
[32m+[m[32m    hp.offR = nRecvTotal;[m
[32m+[m[32m    nSendTotal += hp.nVertexS;[m
[32m+[m[32m    nRecvTotal += hp.nVertexR;[m
[32m+[m[32m    handPairs.push_back(hp);[m
[32m+[m[32m  }[m
[32m+[m
[32m+[m[32m  /*--- Reused every round: the handovers bucketed by the receive marker they cross, and the[m
[32m+[m[32m   *    exchange buffers packed against the halo vertices of every pair at once. ---*/[m
[32m+[m[32m  vector<vector<std::pair<unsigned long, unsigned long>>> handByMarker(config->GetnMarker_All());[m
[32m+[m[32m  vector<unsigned long> tagOut(nRecvTotal), tagIn(nSendTotal);[m
[32m+[m[32m  vector<passivedouble> dirOut(nRecvTotal * nDim), dirIn(nSendTotal * nDim);[m
[32m+[m[32m  vector<CPassiveMPI::Request> handReq(4 * handPairs.size());[m
[32m+[m[32m  vector<CPassiveMPI::Status> handStat(4 * handPairs.size());[m
[32m+[m
   for (unsigned long layer = 1;; ++layer) {[m
     /*--- Every rank runs the same number of rounds, each ends in a collective handover[m
      *    exchange. ---*/[m
[36m@@ -2195,8 +2236,7 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
       }[m
     }[m
 [m
[31m-    /*--- (b) Contention resolved from bids that were all collected before any was granted, so the[m
[31m-     *    outcome does not depend on the order the fronts are visited in. ---*/[m
[32m+[m[32m    /*--- (b) Place every bid; nothing is granted until (c) reads the settled table. ---*/[m
     auto better = [](const CStep& a, const CStep& b) {[m
       if (a.score != b.score) return a.score > b.score;[m
       if (a.dist != b.dist) return a.dist < b.dist;[m
[36m@@ -2255,6 +2295,7 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
 [m
       if (fronts[f].failed) {[m
         fronts[f].alive = 0;[m
[32m+[m[32m        ct[P_RETIRED]++;[m
         emit(f);[m
         continue;[m
       }[m
[36m@@ -2287,42 +2328,60 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
 [m
     /*--- (d) Hand stacks across partition interfaces. The footprint is sent to the owning rank,[m
      *    packed against the receive marker and sent to the rank that marker receives from. ---*/[m
[31m-    for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {[m
[31m-      if (!((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) && (config->GetMarker_All_SendRecv(iMarker) > 0)))[m
[31m-        continue;[m
[31m-[m
[31m-      const auto MarkerS = iMarker, MarkerR = iMarker + 1;[m
[31m-      const auto send_to = config->GetMarker_All_SendRecv(MarkerS) - 1;[m
[31m-      const auto receive_from = abs(config->GetMarker_All_SendRecv(MarkerR)) - 1;[m
[31m-      const auto nVertexS = fine_grid->nVertex[MarkerS];[m
[31m-      const auto nVertexR = fine_grid->nVertex[MarkerR];[m
[31m-[m
[31m-      /*--- Packed against the halo vertices. Tag and direction go separately, the AD MPI wrapper has[m
[31m-       *    no byte type to send a struct. ---*/[m
[31m-      vector<unsigned long> tagOut(nVertexR, 0), tagIn(nVertexS, 0);[m
[31m-      vector<su2double> dirOut(nVertexR * nDim, 0.0), dirIn(nVertexS * nDim, 0.0);[m
[31m-[m
[31m-      for (auto& F : fronts)[m
[31m-        for (auto p : F.handTo) {[m
[31m-          if (haloMarker[p] != static_cast<int>(MarkerR)) continue;[m
[31m-          const auto v = haloVertex[p];[m
[31m-          /*--- Two fronts reaching for one node: the lower tag takes it. ---*/[m
[31m-          if ((tagOut[v] != 0) && (tagOut[v] <= F.handTag)) continue;[m
[31m-          tagOut[v] = F.handTag;[m
[31m-          for (unsigned short d = 0; d < nDim; ++d) dirOut[v * nDim + d] = F.dir[d];[m
[31m-        }[m
[31m-[m
[31m-      SU2_MPI::Sendrecv(tagOut.data(), nVertexR, MPI_UNSIGNED_LONG, receive_from, 2, tagIn.data(), nVertexS,[m
[31m-                        MPI_UNSIGNED_LONG, send_to, 2, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);[m
[31m-      SU2_MPI::Sendrecv(dirOut.data(), nVertexR * nDim, MPI_DOUBLE, receive_from, 3, dirIn.data(), nVertexS * nDim,[m
[31m-                        MPI_DOUBLE, send_to, 3, SU2_MPI::GetComm(), MPI_STATUS_IGNORE);[m
[32m+[m[32m    for (auto& bucket : handByMarker) bucket.clear();[m
[32m+[m[32m    for (unsigned long f = 0; f < fronts.size(); ++f) {[m
[32m+[m[32m      int firstMarker = -1;[m
[32m+[m[32m      bool split = false;[m
[32m+[m[32m      for (auto p : fronts[f].handTo) {[m
[32m+[m[32m        if (haloMarker[p] < 0) continue;[m
[32m+[m[32m        if (firstMarker < 0) firstMarker = haloMarker[p];[m
[32m+[m[32m        else if (haloMarker[p] != firstMarker) split = true;[m
[32m+[m[32m        handByMarker[haloMarker[p]].push_back({f, p});[m
[32m+[m[32m      }[m
[32m+[m[32m      /*--- A footprint owned by more than one rank cannot travel whole. ---*/[m
[32m+[m[32m      if (firstMarker >= 0) {[m
[32m+[m[32m        ct[P_HAND_SENT]++;[m
[32m+[m[32m        ct[P_HAND_SPLIT] += split;[m
[32m+[m[32m      }[m
[32m+[m[32m    }[m
 [m
[31m-      for (auto iVertex = 0ul; iVertex < nVertexS; iVertex++) {[m
[31m-        if (tagIn[iVertex] == 0) continue;[m
[31m-        inherited.push_back({tagIn[iVertex], fine_grid->vertex[MarkerS][iVertex]->GetNode(), {}});[m
[31m-        for (unsigned short d = 0; d < nDim; ++d) inherited.back().dir[d] = dirIn[iVertex * nDim + d];[m
[32m+[m[32m    /*--- Pack every pair, then exchange them all at once so a round costs one wait, not one[m
[32m+[m[32m     *    round trip per neighbour. Tag and direction go separately, there is no byte type to[m
[32m+[m[32m     *    send a struct with. ---*/[m
[32m+[m[32m    std::fill(tagOut.begin(), tagOut.end(), 0ul);[m
[32m+[m[32m    std::fill(dirOut.begin(), dirOut.end(), passivedouble(0.0));[m
[32m+[m
[32m+[m[32m    for (const auto& hp : handPairs)[m
[32m+[m[32m      for (const auto& hand : handByMarker[hp.markerR]) {[m
[32m+[m[32m        const auto& F = fronts[hand.first];[m
[32m+[m[32m        const auto v = hp.offR + haloVertex[hand.second];[m
[32m+[m[32m        /*--- Two fronts reaching for one node: the lower tag takes it. ---*/[m
[32m+[m[32m        if (tagOut[v] != 0) ct[P_HAND_CONTEST]++;[m
[32m+[m[32m        if ((tagOut[v] != 0) && (tagOut[v] <= F.handTag)) continue;[m
[32m+[m[32m        tagOut[v] = F.handTag;[m
[32m+[m[32m        for (unsigned short d = 0; d < nDim; ++d) dirOut[v * nDim + d] = SU2_TYPE::GetValue(F.dir[d]);[m
       }[m
[32m+[m
[32m+[m[32m    unsigned long nReq = 0;[m
[32m+[m[32m    for (const auto& hp : handPairs) {[m
[32m+[m[32m      CPassiveMPI::Irecv(&tagIn[hp.offS], hp.nVertexS, MPI_UNSIGNED_LONG, hp.send_to, 2, CPassiveMPI::GetComm(),[m
[32m+[m[32m                         &handReq[nReq++]);[m
[32m+[m[32m      CPassiveMPI::Irecv(&dirIn[hp.offS * nDim], hp.nVertexS * nDim, MPI_DOUBLE, hp.send_to, 3,[m
[32m+[m[32m                         CPassiveMPI::GetComm(), &handReq[nReq++]);[m
[32m+[m[32m      CPassiveMPI::Isend(&tagOut[hp.offR], hp.nVertexR, MPI_UNSIGNED_LONG, hp.receive_from, 2,[m
[32m+[m[32m                         CPassiveMPI::GetComm(), &handReq[nReq++]);[m
[32m+[m[32m      CPassiveMPI::Isend(&dirOut[hp.offR * nDim], hp.nVertexR * nDim, MPI_DOUBLE, hp.receive_from, 3,[m
[32m+[m[32m                         CPassiveMPI::GetComm(), &handReq[nReq++]);[m
     }[m
[32m+[m[32m    if (nReq > 0) CPassiveMPI::Waitall(static_cast<int>(nReq), handReq.data(), handStat.data());[m
[32m+[m
[32m+[m[32m    for (const auto& hp : handPairs)[m
[32m+[m[32m      for (auto iVertex = 0ul; iVertex < hp.nVertexS; iVertex++) {[m
[32m+[m[32m        const auto v = hp.offS + iVertex;[m
[32m+[m[32m        if (tagIn[v] == 0) continue;[m
[32m+[m[32m        inherited.push_back({tagIn[v], fine_grid->vertex[hp.markerS][iVertex]->GetNode(), {}});[m
[32m+[m[32m        for (unsigned short d = 0; d < nDim; ++d) inherited.back().dir[d] = dirIn[v * nDim + d];[m
[32m+[m[32m      }[m
 [m
     /*--- A front that handed its whole footprint over is finished here, one that handed over only a[m
      *    piece keeps marching on what was left. ---*/[m
[36m@@ -2352,7 +2411,7 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
         if (claimed[p] || fine_grid->nodes->GetAgglomerate(p) || !GeometricalCheck(p, fine_grid, config)) ok = false;[m
         layer0.push_back(p);[m
       }[m
[31m-      /*--- The footprint has to arrive whole and connected. ---*/[m
[32m+[m[32m      /*--- Whatever arrived must be free and form one connected layer. ---*/[m
       if (ok && !IsConnectedLayer(fine_grid, layer0)) ok = false;[m
 [m
       if (ok) {[m
[36m@@ -2364,7 +2423,10 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
           claimed[p] = 1;[m
         }[m
         ct[P_LAYERS]++;[m
[32m+[m[32m        ct[P_HAND_TAKEN]++;[m
         if (fronts[nf].pendingLayers >= fronts[nf].nBlock) emit(nf);[m
[32m+[m[32m      } else {[m
[32m+[m[32m        ct[P_HAND_DROPPED]++;[m
       }[m
       i = j;[m
     }[m
[36m@@ -2397,5 +2459,12 @@[m [mstring CMultiGridGeometry::AgglomerateImplicitLines(unsigned long& Index_CoarseC[m
   out << "  MG level " << iMesh << " paving: " << tot[P_STACKS] << " fronts from " << pairTot[1] << " seeds, "[m
       << pairTot[0] << " CVs covering " << tot[P_COVERED] << " nodes in " << tot[P_LAYERS] << " layers, depth "[m
       << depthMin << " to " << depthMax << "\n";[m
[32m+[m[32m  out << "    patches 1/2/3/4 wide: " << tot[P_PATCH1] << "/" << tot[P_PATCH2] << "/" << tot[P_PATCH3] << "/"[m
[32m+[m[32m      << tot[P_PATCH4] << ", " << tot[P_SEED_CURV] << " seeds refused on curvature, " << tot[P_RETIRED][m
[32m+[m[32m      << " fronts retired early\n";[m
[32m+[m[32m  if (tot[P_HAND_SENT] > 0)[m
[32m+[m[32m    out << "    handovers: " << tot[P_HAND_SENT] << " offered (" << tot[P_HAND_SPLIT] << " split across ranks, "[m
[32m+[m[32m        << tot[P_HAND_CONTEST] << " nodes contested), " << tot[P_HAND_TAKEN] << " adopted, " << tot[P_HAND_DROPPED][m
[32m+[m[32m        << " dropped\n";[m
   return out.str();[m
 }[m
[1mdiff --git a/SU2_CFD/src/drivers/CDriver.cpp b/SU2_CFD/src/drivers/CDriver.cpp[m
[1mindex fa16c68543..c5b4eab97a 100644[m
[1m--- a/SU2_CFD/src/drivers/CDriver.cpp[m
[1m+++ b/SU2_CFD/src/drivers/CDriver.cpp[m
[36m@@ -827,7 +827,7 @@[m [mvoid CDriver::InitializeGeometryFVM(CConfig *config, CGeometry **&geometry) {[m
       geometry[iMGlevel] = nullptr;[m
       break;[m
     }[m
[31m-    pavingReports += coarse_grid->pavingReport;[m
[32m+[m[32m    pavingReports += coarse_grid->GetPavingReport();[m
 [m
     /*--- Compute points surrounding points. ---*/[m
 [m
[36m@@ -857,6 +857,16 @@[m [mvoid CDriver::InitializeGeometryFVM(CConfig *config, CGeometry **&geometry) {[m
   /*--- Held back so they do not interleave with the multigrid level table. ---*/[m
   if (rank == MASTER_NODE) cout << pavingReports;[m
 [m
[32m+[m[32m  /*--- MG_MIN_MESHSIZE is a per-rank floor, so the levels actually built fall with rank count. ---*/[m
[32m+[m[32m  if ((rank == MASTER_NODE) && (requestedMGlevels > 0)) {[m
[32m+[m[32m    if (config->GetnMGLevels() == 0)[m
[32m+[m[32m      cout << "\nWARNING: no multigrid levels used, reduce MG_MIN_MESHSIZE or the number of MPI ranks\n"[m
[32m+[m[32m              "         if you want multigrid.\n" << endl;[m
[32m+[m[32m    else[m
[32m+[m[32m      cout << config->GetnMGLevels() << " multigrid levels used, maximum allowed is " << requestedMGlevels[m
[32m+[m[32m           << ". Change MGLEVEL or MG_MIN_MESHSIZE to use a different number." << endl;[m
[32m+[m[32m  }[m
[32m+[m
   if (config->GetWrt_MultiGrid()) geometry[MESH_0]->ColorMGLevels(config->GetnMGLevels(), geometry);[m
 [m
   /*--- For unsteady simulations, initialize the grid volumes[m
[1mdiff --git a/SU2_CFD/src/integration/CMultiGridIntegration.cpp b/SU2_CFD/src/integration/CMultiGridIntegration.cpp[m
[1mindex 30c93177e1..bd208ed113 100644[m
[1m--- a/SU2_CFD/src/integration/CMultiGridIntegration.cpp[m
[1m+++ b/SU2_CFD/src/integration/CMultiGridIntegration.cpp[m
[36m@@ -63,24 +63,6 @@[m [minline passivedouble ComputeLinSysResRMS(const CSolver* solver) {[m
   return sqrt(result);[m
 }[m
 [m
[31m-/*!\cond PRIVATE[m
[31m- *  Prolongate a coarse-grid field onto the fine grid via constant injection: every fine[m
[31m- *  child gets its parent's value. \c getCoarse returns the coarse-grid block of a point[m
[31m- *  and \c setFine writes it to a fine-grid point.[m
[31m- \endcond */[m
[31m-template <class GetCoarse, class SetFine>[m
[31m-void ProlongateField(CGeometry* geo_coarse, GetCoarse getCoarse, SetFine setFine) {[m
[31m-[m
[31m-  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPoint(), omp_get_num_threads()))[m
[31m-  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {[m
[31m-    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {[m
[31m-      auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);[m
[31m-      setFine(Point_Fine, getCoarse(Point_Coarse));[m
[31m-    }[m
[31m-  }[m
[31m-  END_SU2_OMP_FOR[m
[31m-}[m
[31m-[m
 }  // anonymous namespace[m
 [m
 void CMultiGridIntegration::adaptDampingFactors(CConfig* config, passivedouble crossCycleRatio) {[m
[36m@@ -964,11 +946,16 @@[m [mvoid CMultiGridIntegration::GetProlongated_Correction(unsigned short RunTime_EqS[m
   /*--- Interpolate the coarse-grid correction onto the fine[m
    *    grid and store in LinSysRes. ---*/[m
 [m
[31m-  ProlongateField(geo_coarse,[m
[31m-                  [&](unsigned long iPoint) { return sol_coarse->GetNodes()->GetSolution_Old(iPoint); },[m
[31m-                  [&](unsigned long Point_Fine, const su2double* value) {[m
[31m-                    sol_fine->LinSysRes.SetBlock(Point_Fine, value);[m
[31m-                  });[m
[32m+[m[32m  /*--- Halos too: the correction smoother reads them before its first exchange. ---*/[m
[32m+[m[32m  SU2_OMP_FOR_STAT(roundUpDiv(geo_coarse->GetnPoint(), omp_get_num_threads()))[m
[32m+[m[32m  for (auto Point_Coarse = 0ul; Point_Coarse < geo_coarse->GetnPoint(); Point_Coarse++) {[m
[32m+[m[32m    const auto* Correction = sol_coarse->GetNodes()->GetSolution_Old(Point_Coarse);[m
[32m+[m[32m    for (auto iChildren = 0u; iChildren < geo_coarse->nodes->GetnChildren_CV(Point_Coarse); iChildren++) {[m
[32m+[m[32m      const auto Point_Fine = geo_coarse->nodes->GetChildren_CV(Point_Coarse, iChildren);[m
[32m+[m[32m      sol_fine->LinSysRes.SetBlock(Point_Fine, Correction);[m
[32m+[m[32m    }[m
[32m+[m[32m  }[m
[32m+[m[32m  END_SU2_OMP_FOR[m
 [m
 }[m
 [m
[1mdiff --git a/config_template.cfg b/config_template.cfg[m
[1mindex 582a3272a4..06349adba8 100644[m
[1m--- a/config_template.cfg[m
[1m+++ b/config_template.cfg[m
[36m@@ -1734,7 +1734,8 @@[m [mMG_SMOOTH_OUTPUT= NO[m
 % so the effective global limit scales with the number of ranks.[m
 MG_MIN_MESHSIZE= 500[m
 %[m
[31m-% Enable agglomeration along implicit lines seeded from viscous walls (NO, YES)[m
[32m+[m[32m% Pave the coarse grid with advancing fronts raised from boundaries that carry a[m
[32m+[m[32m% stretched layer normal to themselves (NO, YES)[m
 MG_IMPLICIT_LINES= NO[m
 %[m
 % Number of iterations spent on each mesh during the Full Multigrid (FMG) startup phase. After[m
