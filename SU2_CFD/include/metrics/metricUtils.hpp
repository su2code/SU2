/*!
 * \file metricUtils.hpp
 * \brief Utility functions for mesh adaptation metric computation.
 * \author B. Munguía
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#include <set>
#include <string>
#include <vector>
#include "../../../Common/include/CConfig.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
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
 * \brief Uppercase copy of a sensor name; sensor names match case-insensitively.
 */
inline std::string SensorKey(std::string name) {
  PrintingToolbox::StringToUpperCase(name);
  return name;
}

/*!
 * \brief Build a map of solution variables of the non-flow solvers that can be used as sensors.
 *
 * The names are the volume-output keys of those variables (e.g. NU_TILDE, TKE). The table must
 * follow the solution layout of each solver:
 *   - TURB_SOL:    SA family: NU_TILDE; k-w family: TKE, DISSIPATION
 *   - TRANS_SOL:   LM: INTERMITTENCY, RE_THETA_T
 *   - SPECIES_SOL: species transport: SPECIES_<i>
 *
 * \param[in] config - Definition of the particular problem.
 * \param[in] solver_container - Array of solvers [iSol]; absent solvers contribute no entries.
 * \return Map of variable names to (solver slot, solution index).
 */
inline std::map<std::string, std::pair<unsigned short, unsigned short>>
SolverFieldNameToIndexMap(const CConfig* config, CSolver** solver_container) {
  std::map<std::string, std::pair<unsigned short, unsigned short>> field_map;

  if (solver_container[TURB_SOL] != nullptr) {
    switch (TurbModelFamily(config->GetKind_Turb_Model())) {
      case TURB_FAMILY::SA:
        field_map["NU_TILDE"] = {TURB_SOL, 0};
        break;
      case TURB_FAMILY::KW:
        field_map["TKE"] = {TURB_SOL, 0};
        field_map["DISSIPATION"] = {TURB_SOL, 1};
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (solver_container[TRANS_SOL] != nullptr && config->GetKind_Trans_Model() == TURB_TRANS_MODEL::LM) {
    field_map["INTERMITTENCY"] = {TRANS_SOL, 0};
    field_map["RE_THETA_T"] = {TRANS_SOL, 1};
  }

  if (solver_container[SPECIES_SOL] != nullptr &&
      config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < solver_container[SPECIES_SOL]->GetnVar(); iVar++) {
      field_map["SPECIES_" + std::to_string(iVar)] = {SPECIES_SOL, iVar};
    }
  }

  return field_map;
}

/*!
 * \brief Resolve the METRIC_SENSOR names and store them, in config order, in the flow solver.
 *
 * Each name is matched case-insensitively, in this order:
 *   - CUSTOM:       listed in METRIC_CUSTOM_SENSOR; filled through the Python wrapper (CDriverBase::SetAdaptSensor).
 *   - PRIMITIVE:    a primitive variable of the flow solver.
 *   - COMPUTED:     a supported computed quantity (see ComputedNameToFunctionMap).
 *   - SOLVER_FIELD: a solution variable of a non-flow solver (see SolverFieldNameToIndexMap),
 *                   copied into the flow solver before the Hessians are computed.
 * Any other name is an error. All sensors live in the flow solver, so iSensor = 0 is the first
 * listed sensor, whose Hessian drives the metric.
 *
 * \param[in] config - Configuration containing sensor names
 * \param[in] geometry - Geometry for dimension info
 * \param[in,out] solver_container - Array of solvers [iSol]; the flow solver receives the sensors
 * \return True if at least one sensor was resolved
 */
inline bool ResolveSensorIndices(
    const CConfig* config,
    const CGeometry* geometry,
    CSolver** solver_container) {

  const int rank = SU2_MPI::GetRank();
  constexpr auto no_index = std::numeric_limits<unsigned short>::max();

  if (solver_container[FLOW_SOL] == nullptr || !config->GetFluidProblem()) {
    SU2_MPI::Error("COMPUTE_METRIC requires a flow solver.", CURRENT_FUNCTION);
  }

  /*--- Flow primitives and computed quantities ---*/
  const auto nDim = geometry->GetnDim();
  const bool incompressible = config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE;
  CPrimitiveIndices<unsigned short> indices(incompressible, config->GetNEMOProblem(), nDim, config->GetnSpecies());
  const auto primitive_map = PrimitiveNameToIndexMap(indices);

  std::map<std::string, unsigned short> primitives;
  for (const auto& [name, idx] : primitive_map) primitives[SensorKey(name)] = idx;

  std::map<std::string, std::function<su2double(const su2double*)>> computed;
  for (const auto& [name, fn] : ComputedNameToFunctionMap(primitive_map, nDim, incompressible))
    computed[SensorKey(name)] = fn;

  const auto solver_fields = SolverFieldNameToIndexMap(config, solver_container);

  const std::vector<std::string> custom_names = config->GetMetric_CustomSensorList();
  std::set<std::string> custom;
  for (const auto& name : custom_names) custom.insert(SensorKey(name));

  /*--- Resolve each sensor in config order ---*/
  const std::vector<std::string> sensor_names = config->GetMetric_SensorList();
  std::vector<CSolver::MetricSensorInfo> sensors;
  for (const auto& sensor_name : sensor_names) {
    const std::string key = SensorKey(sensor_name);
    CSolver::MetricSensorInfo info;
    info.name = sensor_name;
    info.var_idx = no_index;
    std::string description;

    if (custom.count(key)) {
      info.type = SensorType::CUSTOM;
      description = "custom (Python wrapper)";
    } else if (primitives.count(key)) {
      info.type = SensorType::PRIMITIVE;
      info.var_idx = primitives.at(key);
      description = "primitive variable " + std::to_string(info.var_idx);
    } else if (computed.count(key)) {
      info.type = SensorType::COMPUTED;
      info.fn = computed.at(key);
      description = "computed";
    } else if (solver_fields.count(key)) {
      info.type = SensorType::SOLVER_FIELD;
      info.iSol = solver_fields.at(key).first;
      info.var_idx = solver_fields.at(key).second;
      description = "solver " + std::to_string(info.iSol) + ", solution variable " + std::to_string(info.var_idx);
    } else {
      std::string available;
      for (const auto& entry : primitives) available += " " + entry.first;
      for (const auto& entry : computed) available += " " + entry.first;
      for (const auto& entry : solver_fields) available += " " + entry.first;
      SU2_MPI::Error("Unknown METRIC_SENSOR " + sensor_name + ". Available sensors:" + available +
                     ". List sensors filled through the Python wrapper in METRIC_CUSTOM_SENSOR.",
                     CURRENT_FUNCTION);
    }

    if (rank == MASTER_NODE) {
      std::cout << "  Sensor '" << sensor_name << "': " << description << "." << std::endl;
    }
    sensors.push_back(info);
  }

  solver_container[FLOW_SOL]->SetMetricSensors(sensors);
  return !sensors.empty();
}

/*!
 * \brief Allocate the metric sensor arrays of the flow solver.
 *
 * Should be called after ResolveSensorIndices() has set the sensors.
 *
 * \param[in,out] solver_container - Array of solvers [iSol]
 */
inline void InitializeMetrics(CSolver** solver_container) {
  const auto nSensors = solver_container[FLOW_SOL]->GetnMetricSensor();
  if (nSensors > 0) {
    solver_container[FLOW_SOL]->AllocateMetricSensorArrays(nSensors);
  }
}

/*!
 * \brief Get the number of metric sensors (all held by the flow solver).
 * \param[in] solver_container - Array of solvers [iSol]
 * \return Number of sensors.
 */
inline unsigned short TotalNumSensors(CSolver** solver_container) {
  return solver_container[FLOW_SOL]->GetnMetricSensor();
}

} // namespace MetricUtils
