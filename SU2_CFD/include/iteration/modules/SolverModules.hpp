#pragma once

#include "CFixedCLModule.hpp"

using Modules    = CModuleList<CConfig, CFixedCLModule>;
using ModulesPtr = std::unique_ptr<Modules>;

extern ModulesPtr modules;