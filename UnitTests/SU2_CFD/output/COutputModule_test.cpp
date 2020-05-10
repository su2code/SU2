#include "catch.hpp"
#include <sstream>

#include "../../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../SU2_CFD/include/solvers/CNSSolver.hpp"
#include "../../../SU2_CFD/include/output/modules/CModuleManager.hpp"


struct __TestCase__ {

  const std::string config_options = "SOLVER= NAVIER_STOKES\n"
                                     "KIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD\n"
                                     "MESH_FORMAT= BOX\n"
                                     "INIT_OPTION=TD_CONDITIONS\n"
                                     "MACH_NUMBER=0.5\n"
                                     "MARKER_CUSTOM= ( x_minus, x_plus, y_minus, y_plus,z_plus, z_minus)\n"
                                     "VISCOSITY_MODEL= CONSTANT_VISCOSITY\n"
                                     "MESH_BOX_SIZE=30,30,30\n"
//                                     "DIRECT_DIFF= MACH\n"
                                     "MESH_BOX_LENGTH=1,1,1\n"
                                     "MESH_BOX_OFFSET=0,0,0\n"
                                     "REF_ORIGIN_MOMENT_X=0.0\n"
                                     "REF_ORIGIN_MOMENT_Y=0.0\n"
                                     "REF_ORIGIN_MOMENT_Z=0.0\n";
  std::stringstream ss;
  CConfig config;
  CPhysicalGeometry aux_geometry, geometry;
  CSolver** solver;
  __TestCase__() : ss(config_options),
    config(ss, SU2_CFD, false),
    aux_geometry(&config, 0, 1),
    geometry(&aux_geometry, &config) {
    streambuf* orig_buf = cout.rdbuf();

    cout.rdbuf(nullptr);
    geometry.SetSendReceive(&config);
    geometry.SetBoundaries(&config);
    geometry.SetPoint_Connectivity();
    geometry.SetElement_Connectivity();
    geometry.SetBoundVolume();

    geometry.SetEdges();
    geometry.SetVertex(&config);
    geometry.SetCoord_CG();
    geometry.SetControlVolume(&config, ALLOCATE);
    geometry.SetBoundControlVolume(&config, ALLOCATE);
    geometry.FindNormal_Neighbor(&config);
    geometry.SetGlobal_to_Local_Point();
    geometry.PreprocessP2PComms(&geometry,&config);

    solver = new CSolver*[MAX_SOLS];
    solver[FLOW_SOL] = new CNSSolver(&geometry, &config, 0);
    cout.rdbuf(orig_buf);

  }
  ~__TestCase__(){
    delete solver[FLOW_SOL];
    delete [] solver;
  }

};

TEST_CASE("Aerodynamics Module", "[Output Module]"){


  __TestCase__ TestCase;
  CModuleManager<ModuleList<CCommonModule, CAerodynamicsModule, CFlowCoefficientModule, CConvergenceModule, CResidualModule, CDirectDiffModule, CUserFunctionModule>> modules(&TestCase.config);

  CGeometry* geometry = &TestCase.geometry;
  TestCase.solver[FLOW_SOL]->SetInitialCondition(&geometry, &TestCase.solver, &TestCase.config, 0);

  TestCase.solver[FLOW_SOL]->Preprocessing(geometry, TestCase.solver, &TestCase.config, 0, 0, RUNTIME_FLOW_SYS, false);

  SolverDataContainer data;

  data.solver = TestCase.solver;
  data.geometry = &TestCase.geometry;
  data.config = &TestCase.config;

  modules.LoadData(&data);

  for (const auto& fields : modules.GetHistoryFields().GetReferencesAll()){
    cout << fields->first << " " << fields->second.value << endl;
  }
//  for (const auto& fields : modules.GetVolumeFields().GetReferencesAll()){
//    cout << fields->first << " " << fields->second.value << endl;
//  }
//  cout << TestCase.config.GetVelocity_FreeStream()[0] << endl;
}