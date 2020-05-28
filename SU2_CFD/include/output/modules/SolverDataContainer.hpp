#pragma once

#include "../../../Common/include/CConfig.hpp"

class CSolver;
class CGeometry;
class CConfig;

struct IterationInfo{

  unsigned long Iter    ;
  unsigned long TimeIter;

};

struct SolverData{

  CConfig* config;
  CGeometry* geometry;
  CSolver** solver;

};


struct PointInfo{

  unsigned long iPoint;
  unsigned long iVertex;
  uint iMarker;

};

