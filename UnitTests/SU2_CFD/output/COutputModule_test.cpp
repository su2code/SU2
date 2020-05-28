#include "catch.hpp"
#include <sstream>

#include "../../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../SU2_CFD/include/solvers/CNSSolver.hpp"
#include "../../../SU2_CFD/include/output/modules/CModuleManager.hpp"
#include "../../../SU2_CFD/include/output/modules/CCommonModule.hpp"
#include "../../../SU2_CFD/include/output/modules/CAerodynamicsModule.hpp"
#include "../../../SU2_CFD/include/output/modules/CFVMBaseModule.hpp"

struct __TestCase__ {

  const std::string config_options = "SOLVER= NAVIER_STOKES\n"
                                     "KIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD\n"
                                     "MESH_FORMAT= BOX\n"
                                     "INIT_OPTION=TD_CONDITIONS\n"
                                     "MACH_NUMBER=0.5\n"
                                     "MARKER_CUSTOM= ( x_minus, x_plus, y_minus, y_plus,z_plus, z_minus)\n"
                                     "VISCOSITY_MODEL= CONSTANT_VISCOSITY\n"
                                     "MESH_BOX_SIZE=30,30,30\n"
                                     "MESH_BOX_LENGTH=1,1,1\n"
                                     "MESH_BOX_OFFSET=0,0,0\n"
                                     "REF_ORIGIN_MOMENT_X=0.0\n"
                                     "REF_ORIGIN_MOMENT_Y=0.0\n"
                                     "REF_ORIGIN_MOMENT_Z=0.0\n";
  std::stringstream ss;
  CConfig *config;
  CGeometry *aux_geometry,
                    *geometry;
  CSolver** solver;
  __TestCase__() : ss(config_options) {
    streambuf* orig_buf = cout.rdbuf();

    cout.rdbuf(nullptr);

    config       = new CConfig(ss, SU2_CFD, false);
    aux_geometry = new CPhysicalGeometry(config, 0, 1);
    geometry     = new CPhysicalGeometry(aux_geometry, config);
    delete aux_geometry;

    geometry->SetSendReceive(config);
    geometry->SetBoundaries(config);
    geometry->SetPoint_Connectivity();
    geometry->SetElement_Connectivity();
    geometry->SetBoundVolume();

    geometry->SetEdges();
    geometry->SetVertex(config);
    geometry->SetCoord_CG();
    geometry->SetControlVolume(config, ALLOCATE);
    geometry->SetBoundControlVolume(config, ALLOCATE);
    geometry->FindNormal_Neighbor(config);
    geometry->SetGlobal_to_Local_Point();
    geometry->PreprocessP2PComms(geometry,config);

    solver = new CSolver*[MAX_SOLS];
    solver[FLOW_SOL] = new CNSSolver(geometry, config, 0);
    cout.rdbuf(orig_buf);

  }
  ~__TestCase__(){
    delete solver[FLOW_SOL];
    delete [] solver;
    delete geometry;
    delete config;
  }

};

std::unique_ptr<__TestCase__> TestCase;

TEST_CASE("Common module", "[Output Module]"){

  TestCase = std::unique_ptr<__TestCase__>(new __TestCase__());
  TestCase->solver[FLOW_SOL]->SetInitialCondition(&TestCase->geometry, &TestCase->solver, TestCase->config, 0);

  CModuleManager<ModuleList<CCommonModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{100,0});

  REQUIRE(modules.GetHistoryFields().GetValueByKey("INNER_ITER") == Approx(100));

}

TEST_CASE("FVM Base Module", "[Output Module]"){
  CModuleManager<ModuleList<CFVMBaseModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{0,0});

  auto Area =   modules.GetHistoryFields().GetValueByKey("AREA@x_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@x_plus")
              + modules.GetHistoryFields().GetValueByKey("AREA@y_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@y_plus")
              + modules.GetHistoryFields().GetValueByKey("AREA@z_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@z_plus");

  REQUIRE(modules.GetHistoryFields().GetValueByKey("AREA") == Approx(Area));

}

TEST_CASE("Aerodynamics Module", "[Output Module]"){


//  __TestCase__ TestCase;
  CModuleManager<ModuleList<CAerodynamicsModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  TestCase->solver[FLOW_SOL]->Preprocessing(TestCase->geometry, TestCase->solver, TestCase->config, 0, 0, RUNTIME_FLOW_SYS, false);

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{0,0});

  REQUIRE(modules.GetHistoryFields().GetValueByKey("LIFT") == Approx(0.245231));

}