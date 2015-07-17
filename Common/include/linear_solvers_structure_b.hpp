#pragma once

#include "config_structure.hpp"

#ifdef CODI_REVERSE_TYPE
class CSysSolve_b{

public:
  static void Solve_b(AD::CheckpointHandler *data);
  static void Delete_b(AD::CheckpointHandler *data);
};
#endif
