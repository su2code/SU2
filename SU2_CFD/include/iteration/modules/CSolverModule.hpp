#pragma once
#include "../../Common/include/toolboxes/CModule.hpp"

class CGeometry;
class CConfig;
class CSolver;
class COutput;

class CSolverModule : public CModule<CSolverModule>{

public:
  CSolverModule(bool enabled_ = true) : CModule(enabled_) {};

  virtual void PreIterationHook(COutput *output, CConfig* config, CGeometry* geometry, CSolver** solver) {};
  CREATE_ACTION(PreIterationHook, COutput *output, CConfig* config, CGeometry* geometry, CSolver** solver);
  virtual void PostIterationHook(COutput *output, CConfig* config, CGeometry* geometry, CSolver** solver) {};
  CREATE_ACTION(PostIterationHook, COutput *output, CConfig* config, CGeometry* geometry, CSolver** solver);
};

