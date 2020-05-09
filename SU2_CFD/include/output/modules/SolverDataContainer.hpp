#pragma once

#include "../../../Common/include/CConfig.hpp"

class CSolver;
class CVertex;
class CGeometry;
class CConfig;
class CSurfaceElementFEM;

class OutputData{

public:
  inline virtual ~OutputData() {};
};

class SolverDataContainer : public OutputData{
public:
  const CGeometry *geometry;
  const CConfig   *config;
  CSolver  **solver;
  CVertex*  vertex;
  unsigned long iPoint;
  unsigned long Iter, TimeIter;

  const CSurfaceElementFEM *surfElem;
  const su2double *sol;
  su2double weight;
  int i;
  inline ~SolverDataContainer(){}

};
