/*!
 * \file UnitQuadTestCase.hpp
 * \brief Simple unit quad test to be used in unit tests.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include <string>

#include "../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../SU2_CFD/include/solvers/CNSSolver.hpp"

struct UnitQuadTestCase {
  std::string config_options =
      "SOLVER= NAVIER_STOKES\n"
      "KIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD\n"
      "MESH_FORMAT= BOX\n"
      "INIT_OPTION=TD_CONDITIONS\n"
      "MACH_NUMBER=0.5\n"
      "MARKER_HEATFLUX= (y_minus, 0.0, y_plus, 0.0)\n"
      "MARKER_CUSTOM= ( x_minus, x_plus, z_plus, z_minus)\n"
      "VISCOSITY_MODEL= CONSTANT_VISCOSITY\n"
      "MESH_BOX_SIZE=5,5,5\n"
      "MESH_BOX_LENGTH=1,1,1\n"
      "MESH_BOX_OFFSET=0,0,0\n"
      "REF_ORIGIN_MOMENT_X=0.0\n"
      "REF_ORIGIN_MOMENT_Y=0.0\n"
      "REF_ORIGIN_MOMENT_Z=0.0\n";
  std::unique_ptr<CConfig> config;
  std::unique_ptr<CGeometry> geometry;
  CSolver** solver{nullptr};
  streambuf* orig_buf{nullptr};
  UnitQuadTestCase() : orig_buf(cout.rdbuf()) {}

  /*!
   * \brief Add a line to the base config string stream
   * \param[in] optionLine - String containing the option(s)
   */
  void AddOption(const std::string& optionLine) { config_options += optionLine + "\n"; }

  /*!
   * \brief Initialize the config structure
   */
  void InitConfig() {
    cout.rdbuf(nullptr);
    stringstream ss(config_options);
    config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_COMPONENT::SU2_CFD, false));
    cout.rdbuf(orig_buf);
  }

  /*!
   * \brief Initialize the solver array
   */
  void InitSolver() {
    cout.rdbuf(nullptr);
    solver = CSolverFactory::CreateSolverContainer(config.get()->GetKind_Solver(), config.get(), geometry.get(), 0);
    cout.rdbuf(orig_buf);
  }

  /*!
   * \brief Initialize the geometry
   */
  void InitGeometry() {
    cout.rdbuf(nullptr);
    {
      auto aux_geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(config.get(), 0, 1));
      geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux_geometry.get(), config.get()));
    }
    geometry->SetSendReceive(config.get());
    geometry->SetBoundaries(config.get());
    geometry->SetPoint_Connectivity();
    geometry->SetElement_Connectivity();
    geometry->SetBoundVolume();
    geometry->Check_IntElem_Orientation(config.get());
    geometry->Check_BoundElem_Orientation(config.get());
    geometry->SetEdges();
    geometry->SetVertex(config.get());
    geometry->SetControlVolume(config.get(), ALLOCATE);
    geometry->SetBoundControlVolume(config.get(), ALLOCATE);
    geometry->FindNormal_Neighbor(config.get());
    geometry->SetGlobal_to_Local_Point();
    geometry->PreprocessP2PComms(geometry.get(), config.get());

    cout.rdbuf(orig_buf);
  }

  /*!
   * \brief Desctructor
   */
  ~UnitQuadTestCase() {
    if (solver != nullptr) delete solver[FLOW_SOL];
    delete[] solver;
  }
};
