#pragma once

#include <string>

#include "../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../SU2_CFD/include/solvers/CNSSolver.hpp"

struct UnitQuadTestCase {

  const std::string config_options = "SOLVER= NAVIER_STOKES\n"
                                     "KIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD\n"
                                     "MESH_FORMAT= BOX\n"
                                     "INIT_OPTION=TD_CONDITIONS\n"
                                     "MACH_NUMBER=0.5\n"
                                     "MARKER_HEATFLUX= (y_minus, 0.0, y_plus, 0.0)\n"
                                     "MARKER_CUSTOM= ( x_minus, x_plus, z_plus, z_minus)\n"
                                     "VISCOSITY_MODEL= CONSTANT_VISCOSITY\n"
                                     "CONV_FIELD=RMS_DENSITY\n"
                                     "MESH_BOX_SIZE=5,5,5\n"
                                     "MESH_BOX_LENGTH=1,1,1\n"
                                     "MESH_BOX_OFFSET=0,0,0\n"
                                     "REF_ORIGIN_MOMENT_X=0.0\n"
                                     "REF_ORIGIN_MOMENT_Y=0.0\n"
                                     "REF_ORIGIN_MOMENT_Z=0.0\n";
  std::stringstream ss;
  std::unique_ptr<CConfig>   config;
  std::unique_ptr<CGeometry> geometry;
  CSolver** solver{nullptr};
  streambuf* orig_buf{nullptr};
  UnitQuadTestCase() : ss(config_options), orig_buf(cout.rdbuf()) {}

  void InitConfig(){
    cout.rdbuf(nullptr);
    config = std::unique_ptr<CConfig>(new CConfig(ss, SU2_CFD, false));
    cout.rdbuf(orig_buf);
  }

  void InitSolver(){
    cout.rdbuf(nullptr);
    solver = new CSolver*[MAX_SOLS];
    solver[FLOW_SOL] = new CNSSolver(geometry.get(), config.get(), 0);
    cout.rdbuf(orig_buf);

  }

  void InitGeometry(){
    cout.rdbuf(nullptr);
    {
      auto aux_geometry = std::unique_ptr<CGeometry>(new CPhysicalGeometry(config.get(), 0, 1));
      geometry     = std::unique_ptr<CGeometry>(new CPhysicalGeometry(aux_geometry.get(), config.get()));
    }
    geometry->SetSendReceive(config.get());
    geometry->SetBoundaries(config.get());
    geometry->SetPoint_Connectivity();
    geometry->SetElement_Connectivity();
    geometry->SetBoundVolume();

    geometry->SetEdges();
    geometry->SetVertex(config.get());
    geometry->SetCoord_CG();
    geometry->SetControlVolume(config.get(), ALLOCATE);
    geometry->SetBoundControlVolume(config.get(), ALLOCATE);
    geometry->FindNormal_Neighbor(config.get());
    geometry->SetGlobal_to_Local_Point();
    geometry->PreprocessP2PComms(geometry.get(),config.get());

    cout.rdbuf(orig_buf);

  }

  ~UnitQuadTestCase(){
    if (solver != nullptr){
      delete solver[FLOW_SOL];
      delete [] solver;
    }
  }

};