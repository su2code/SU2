/*!
 * \file metricUtils.hpp
 * \brief Utility functions for mesh adaptation metric computation.
 * \author B. Munguía
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include "../../../Common/include/CConfig.hpp"
#include "../solvers/CSolver.hpp"
#include "../variables/CPrimitiveIndices.hpp"
#include "../../../Common/include/parallelization/mpi_structure.hpp"

namespace MetricUtils {

/*!
 * \brief Build a map of supported computed sensor names to their point-wise evaluator lambdas.
 *
 * Computed sensors are officially-supported computed quantities that are not stored directly
 * as primitive variables. Supported sensors:
 *   - VELOCITY  velocity magnitude (all flow solvers)
 *   - MACH      local Mach number  (compressible flow solvers only)
 *
 * Each lambda has the signature \c su2double(const su2double* prim) where \p prim is a pointer
 * to the primitive variable array for a single point.
 *
 * \param[in] primitive_map  - Map of primitive variable names to their array indices.
 * \param[in] nDim           - Number of spatial dimensions.
 * \param[in] incompressible - True for incompressible solvers (MACH is excluded).
 * \return Map of computed sensor names to evaluator lambdas.
 */
inline std::map<std::string, std::function<su2double(const su2double*)>>
ComputedNameToFunctionMap(const std::map<std::string, unsigned short>& primitive_map,
                          unsigned long nDim, bool incompressible) {
  using SensorFn = std::function<su2double(const su2double*)>;
  std::map<std::string, SensorFn> computed_map;

  const unsigned short vel_x = primitive_map.at("VELOCITY_X");
  const unsigned short vel_y = primitive_map.at("VELOCITY_Y");
  const unsigned short vel_z = (nDim == 3) ? primitive_map.at("VELOCITY_Z")
                                            : std::numeric_limits<unsigned short>::max();

  /*--- Velocity magnitude: supported by all flow solvers ---*/
  computed_map["VELOCITY"] = [vel_x, vel_y, vel_z, nDim](const su2double* prim) -> su2double {
    su2double vel2 = prim[vel_x] * prim[vel_x] + prim[vel_y] * prim[vel_y];
    if (nDim == 3) vel2 += prim[vel_z] * prim[vel_z];
    return std::sqrt(vel2);
  };

  /*--- Mach number: compressible solvers only ---*/
  if (!incompressible) {
    const unsigned short a_idx = primitive_map.at("SOUND_SPEED");
    computed_map["MACH"] = [vel_x, vel_y, vel_z, a_idx, nDim](const su2double* prim) -> su2double {
      su2double vel2 = prim[vel_x] * prim[vel_x] + prim[vel_y] * prim[vel_y];
      if (nDim == 3) vel2 += prim[vel_z] * prim[vel_z];
      return std::sqrt(vel2) / std::max(std::abs(prim[a_idx]), su2double(1e-20));
    };
  }

  return computed_map;
}

/*!
 * \brief Resolve sensor names to solver and variable indices, storing them directly in each solver.
 *
 * Maps sensor names from config (e.g., "DENSITY", "TEMPERATURE") to their corresponding
 * solver index and variable index within that solver. Works for both flow solvers
 * (using primitive variables) and non-flow solvers (using solution fields).
 * Directly stores the resolved sensors in each solver via SetMetricSensors().
 *
 * \param[in] config - Configuration containing sensor names
 * \param[in] geometry - Geometry for dimension info
 * \param[in,out] solver_container - Array of solvers [iSol] - indices will be set
 * \return True if all sensors were successfully resolved
 */
inline bool ResolveSensorIndices(
    const CConfig* config,
    const CGeometry* geometry,
    CSolver** solver_container) {

  const int rank = SU2_MPI::GetRank();

  /*--- Get sensor names from config ---*/
  std::vector<std::string> sensor_names = config->GetMetric_SensorList();

  /*--- Group sensors by solver (in config order) ---*/
  std::map<unsigned short, std::vector<CSolver::MetricSensorInfo>> sensors_by_solver;

  /*--- Build maps of available variables across all solvers.
   *    var_map:     name -> (iSol, prim_idx)  for PRIMITIVE sensors
   *    computed_map: name -> evaluator lambda   for COMPUTED sensors ---*/
  std::map<std::string, std::pair<unsigned short, unsigned short>> var_map;
  std::map<std::string, std::function<su2double(const su2double*)>> computed_map;

  for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++) {
    if (solver_container[iSol] == nullptr) continue;

    std::string solver_name = solver_container[iSol]->GetSolverName();

    /*--- For flow solvers, get primitive variable names ---*/
    if (solver_name.find("FLOW") != std::string::npos ||
        solver_name.find("EULER") != std::string::npos ||
        solver_name.find("NAVIER_STOKES") != std::string::npos ||
        solver_name.find("RANS") != std::string::npos ||
        solver_name.find("INC") != std::string::npos ||
        solver_name.find("NEMO") != std::string::npos) {

      const auto nDim = geometry->GetnDim();
      const auto nSpecies = config->GetnSpecies();
      const bool incompressible = config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE;
      const bool nemo = config->GetKind_FluidModel() == ENUM_FLUIDMODEL::MUTATIONPP ||
                        config->GetKind_FluidModel() == ENUM_FLUIDMODEL::SU2_NONEQ;

      CPrimitiveIndices<unsigned short> indices(incompressible, nemo, nDim, nSpecies);
      std::map<std::string, unsigned short> primitive_map = PrimitiveNameToIndexMap(indices);

      for (const auto& [varname, varidx] : primitive_map) {
        var_map[varname] = std::make_pair(iSol, varidx);
      }

      /*--- Build computed sensor map. VELOCITY is supported by all flow solvers;
       *    MACH is excluded for incompressible (see ComputedNameToFunctionMap). ---*/
      computed_map = ComputedNameToFunctionMap(primitive_map, nDim, incompressible);
    } else {
      /*--- For non-flow solvers, use solution fields ---*/
      std::vector<std::string> solution_fields = solver_container[iSol]->GetSolutionFields();

      /*--- Remove PointID if present (first entry) ---*/
      if (!solution_fields.empty() && solution_fields[0].find("PointID") != std::string::npos) {
        solution_fields.erase(solution_fields.begin());
      }

      /*--- Remove quotation marks and build map ---*/
      unsigned short idx = 0;
      for (auto& field_name : solution_fields) {
        if (field_name.size() >= 2 && field_name.front() == '\"' && field_name.back() == '\"') {
          field_name = field_name.substr(1, field_name.size() - 2);
        }
        var_map[field_name] = std::make_pair(iSol, idx++);
      }
    }
  }

  /*--- Resolve each sensor in config order, grouping by solver.
   *    Three categories are distinguished:
   *      PRIMITIVE: index into primitive variable array (prim_idx is valid).
   *      COMPUTED:  supported computed quantity (e.g. Mach); evaluator lambda stored in fn.
   *      CUSTOM:    filled externally via Python wrapper (CDriverBase::AdaptSensors).
   *    Config order is preserved so iSensor=0 is always the first listed sensor
   *    (whose Hessian drives the metric and is normalised). ---*/
  for (const auto& sensor_name : sensor_names) {
    auto it = var_map.find(sensor_name);
    if (it != var_map.end()) {
      /*--- PRIMITIVE sensor ---*/
      const unsigned short iSol = it->second.first;
      const unsigned short var_index = it->second.second;
      sensors_by_solver[iSol].push_back(
          {var_index, sensor_name, SensorType::PRIMITIVE, {}});
      if (rank == MASTER_NODE) {
        std::cout << "  Primitive sensor '" << sensor_name << "' resolved to solver "
                  << iSol << ", variable index " << var_index << std::endl;
      }
    } else {
      auto dit = computed_map.find(sensor_name);
      if (dit != computed_map.end()) {
        /*--- COMPUTED sensor ---*/
        sensors_by_solver[FLOW_SOL].push_back(
            {std::numeric_limits<unsigned short>::max(), sensor_name,
             SensorType::COMPUTED, dit->second});
        if (rank == MASTER_NODE) {
          std::cout << "  Computed sensor '" << sensor_name << "' registered." << std::endl;
        }
      } else {
        /*--- CUSTOM sensor (must be populated via Python wrapper) ---*/
        if (rank == MASTER_NODE) {
          std::cout << "  Custom sensor '" << sensor_name << "' detected." << std::endl;
        }
        sensors_by_solver[FLOW_SOL].push_back(
            {std::numeric_limits<unsigned short>::max(), sensor_name,
             SensorType::CUSTOM, {}});
      }
    }
  }

  /*--- Set sensor list directly in each solver ---*/
  for (const auto& [iSol, sensors] : sensors_by_solver) {
    solver_container[iSol]->SetMetricSensors(sensors);
  }

  /*--- Return true if at least one sensor was resolved or reserved ---*/
  return !sensors_by_solver.empty();
}

/*!
 * \brief Allocate metric sensor arrays in all solvers.
 *
 * Allocates the necessary arrays for metric computation in each solver that has
 * metric sensors. Should be called after ResolveSensorIndices() has set the indices.
 *
 * \param[in,out] solver_container - Array of solvers [iSol]
 */
inline void InitializeMetrics(CSolver** solver_container) {

  /*--- Allocate arrays in each solver that has metric sensors ---*/
  for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++) {
    if (solver_container[iSol] == nullptr) continue;

    const auto nSensors = solver_container[iSol]->GetnMetricSensor();
    if (nSensors > 0) {
      solver_container[iSol]->AllocateMetricSensorArrays(nSensors);
    }
  }
}

/*!
 * \brief Get total number of sensors in all solvers.
 * \param[in] solver_container - Array of solvers [iSol]
 * \return Number of sensors in all solvers.
 */
inline unsigned short TotalNumSensors(CSolver** solver_container) {
  unsigned short num_sensor = 0;
  for (unsigned short iSol = 0; iSol < MAX_SOLS; iSol++) {
    if (solver_container[iSol] != nullptr) {
      num_sensor += solver_container[iSol]->GetnMetricSensor();
    }
  }
  return num_sensor;
}

} // namespace MetricUtils
