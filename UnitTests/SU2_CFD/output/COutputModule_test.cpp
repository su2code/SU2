#include "catch.hpp"
#include <sstream>

#include "../../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../SU2_CFD/include/solvers/CNSSolver.hpp"
#include "../../../SU2_CFD/include/output/modules/CModuleManager.hpp"
#include "../../../SU2_CFD/include/output/modules/CCommonModule.hpp"
#include "../../../SU2_CFD/include/output/CFlowCompOutput.hpp"
#include "../../../SU2_CFD/include/output/modules/CConvergenceModule.hpp"
#include "../../../SU2_CFD/include/output/modules/CAerodynamicsModule.hpp"
#include "../../../SU2_CFD/include/output/modules/CFVMBaseModule.hpp"
#include "../../../SU2_CFD/include/output/modules/CResidualModule.hpp"

struct __TestCase__ {

  const std::string config_options = "SOLVER= NAVIER_STOKES\n"
                                     "KIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD\n"
                                     "MESH_FORMAT= BOX\n"
                                     "INIT_OPTION=TD_CONDITIONS\n"
                                     "MACH_NUMBER=0.5\n"
                                     "MARKER_CUSTOM= ( x_minus, x_plus, y_minus, y_plus,z_plus, z_minus)\n"
                                     "VISCOSITY_MODEL= CONSTANT_VISCOSITY\n"
                                     "CONV_FIELD=RMS_DENSITY, DRAG\n"
                                     "CONV_STARTITER=5\n"
                                     "CONV_CAUCHY_ELEMS=3\n"
                                     "CONV_CAUCHY_EPS=1e-3\n"
                                     "MESH_BOX_SIZE=5,5,5\n"
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

  CHECK(modules.GetHistoryFields().GetValueByKey("INNER_ITER") == Approx(100));

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

  CHECK(modules.GetHistoryFields().GetValueByKey("AREA") == Approx(Area));

}

TEST_CASE("Aerodynamics Module", "[Output Module]"){

  CModuleManager<ModuleList<CAerodynamicsModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  TestCase->solver[FLOW_SOL]->Preprocessing(TestCase->geometry, TestCase->solver, TestCase->config, 0, 0, RUNTIME_FLOW_SYS, false);

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{0,0});

  CHECK(modules.GetHistoryFields().GetValueByKey("LIFT") == Approx(0.5338212762));
  CHECK(modules.GetHistoryFields().GetValueByKey("DRAG") == Approx(-2.3919536834));

}

TEST_CASE("Convergence Module", "[Output Module]"){

  CModuleManager<ModuleList<CConvergenceModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  modules.Clear();

  modules.GetHistoryFields().AddItem("RMS_DENSITY", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("DRAG", COutputField("", ScreenOutputFormat::SCIENTIFIC, "DRAG", "", FieldType::AUTO_COEFFICIENT));

  modules.Init(TestCase->config);

  REQUIRE(modules.GetHistoryFields().CheckKey("CAUCHY_DRAG"));

  modules.GetHistoryFields().GetItemByKey("RMS_DENSITY").value = -15.0;
  auto val1 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.001;

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{0,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == 1.0);
  auto val2 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.002;

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{1,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == Approx((0 + 0 + abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetValueByKey("CONVERGENCE") == 0.0);

  auto val3 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.003;

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{50,0});

  auto val4 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.004;

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{51,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == Approx((abs(val4-val3)/abs(val4) +
                                                                   abs(val3-val2)/abs(val3) +
                                                                   abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetValueByKey("CONVERGENCE") == 1.0);

}

TEST_CASE("Residual Module", "[Residual Module]"){

  CModuleManager<ModuleList<CResidualModule>>  modules(TestCase->config, TestCase->geometry->GetnDim());

  modules.Clear();

  modules.GetHistoryFields().AddItem("RMS_DENSITY", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("RMS_MOMENTUM_X", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("RMS_MOMENTUM_Y", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("RMS_MOMENTUM_Z", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("RMS_ENERGY", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));

  auto rms_density = modules.GetHistoryFields().GetItemByKey("RMS_DENSITY").value    = -5.765;
  auto rms_mom_x   = modules.GetHistoryFields().GetItemByKey("RMS_MOMENTUM_X").value = -6.52;
  auto rms_mom_y   = modules.GetHistoryFields().GetItemByKey("RMS_MOMENTUM_Y").value = -8.52;
  auto rms_mom_z   = modules.GetHistoryFields().GetItemByKey("RMS_MOMENTUM_Z").value = -3.32;
  auto rms_energy  = modules.GetHistoryFields().GetItemByKey("RMS_ENERGY").value     = -4.64;

  modules.Init(TestCase->config);

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{0,0});

  REQUIRE(modules.GetHistoryFields().CheckKey("AVG_RMS_RES"));
  CHECK(modules.GetHistoryFields().GetValueByKey("AVG_RMS_RES") == (rms_density+rms_mom_x+rms_mom_y+rms_mom_z+rms_energy)/5);

  su2double reduction = -1.453;
  rms_density += reduction;
  modules.GetHistoryFields().GetItemByKey("RMS_DENSITY").value    = rms_density;

  modules.LoadData({TestCase->config, TestCase->geometry, TestCase->solver},{10,0});

  REQUIRE(modules.GetHistoryFields().CheckKey("REL_RMS_DENSITY"));
  CHECK(modules.GetHistoryFields().GetValueByKey("REL_RMS_DENSITY") == Approx(reduction));

}

