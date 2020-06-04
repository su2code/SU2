#include "catch.hpp"
#include <sstream>

#include "../../UnitQuadTestCase.hpp"

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

std::unique_ptr<UnitQuadTestCase> ModuleTestCase;

TEST_CASE("Common module", "[Output Module]"){

  ModuleTestCase = std::unique_ptr<UnitQuadTestCase>(new UnitQuadTestCase());

  ModuleTestCase->AddOption("CONV_FIELD=RMS_DENSITY, DRAG");
  ModuleTestCase->AddOption("CONV_STARTITER=5");
  ModuleTestCase->AddOption("CONV_CAUCHY_EPS=1e-3");
  ModuleTestCase->AddOption("CONV_CAUCHY_ELEMS=3");

  ModuleTestCase->InitConfig();
  ModuleTestCase->InitGeometry();
  ModuleTestCase->InitSolver();

  CGeometry* geometry = ModuleTestCase->geometry.get();
  ModuleTestCase->solver[FLOW_SOL]->SetInitialCondition(&geometry, &ModuleTestCase->solver, ModuleTestCase->config.get(), 0);

  CModuleManager<ModuleList<CCommonModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{100,0});

  CHECK(modules.GetHistoryFields().GetValueByKey("INNER_ITER") == Approx(100));

}

TEST_CASE("FVM Base Module", "[Output Module]"){
  CModuleManager<ModuleList<CFVMBaseModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  auto Area =   modules.GetHistoryFields().GetValueByKey("AREA@x_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@x_plus")
              + modules.GetHistoryFields().GetValueByKey("AREA@y_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@y_plus")
              + modules.GetHistoryFields().GetValueByKey("AREA@z_minus")
              + modules.GetHistoryFields().GetValueByKey("AREA@z_plus");

  CHECK(modules.GetHistoryFields().GetValueByKey("AREA") == Approx(Area));

}

TEST_CASE("Aerodynamics Module", "[Output Module]"){

  CModuleManager<ModuleList<CAerodynamicsModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  ModuleTestCase->solver[FLOW_SOL]->Preprocessing(ModuleTestCase->geometry.get(), ModuleTestCase->solver, ModuleTestCase->config.get(), 0, 0, RUNTIME_FLOW_SYS, false);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  CHECK(modules.GetHistoryFields().GetValueByKey("LIFT") == Approx(0.5338212762));
  CHECK(modules.GetHistoryFields().GetValueByKey("DRAG") == Approx(-2.3919536834));

}

TEST_CASE("Convergence Module", "[Output Module]"){

  CModuleManager<ModuleList<CConvergenceModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.Clear();

  modules.GetHistoryFields().AddItem("RMS_DENSITY", COutputField("", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL));
  modules.GetHistoryFields().AddItem("DRAG", COutputField("", ScreenOutputFormat::SCIENTIFIC, "DRAG", "", FieldType::AUTO_COEFFICIENT));

  modules.Init(ModuleTestCase->config.get());

  REQUIRE(modules.GetHistoryFields().CheckKey("CAUCHY_DRAG"));

  modules.GetHistoryFields().GetItemByKey("RMS_DENSITY").value = -15.0;
  auto val1 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.001;

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == 1.0);
  auto val2 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.002;

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{1,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == Approx((0 + 0 + abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetValueByKey("CONVERGENCE") == 0.0);

  auto val3 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.003;

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{50,0});

  auto val4 = modules.GetHistoryFields().GetItemByKey("DRAG").value = 1.004;

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{51,0});

  CHECK(modules.GetHistoryFields()["CAUCHY_DRAG"].value == Approx((abs(val4-val3)/abs(val4) +
                                                                   abs(val3-val2)/abs(val3) +
                                                                   abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetValueByKey("CONVERGENCE") == 1.0);

}

TEST_CASE("Residual Module", "[Residual Module]"){

  CModuleManager<ModuleList<CResidualModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

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

  modules.Init(ModuleTestCase->config.get());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  REQUIRE(modules.GetHistoryFields().CheckKey("AVG_RMS_RES"));
  CHECK(modules.GetHistoryFields().GetValueByKey("AVG_RMS_RES") == (rms_density+rms_mom_x+rms_mom_y+rms_mom_z+rms_energy)/5);

  su2double reduction = -1.453;
  rms_density += reduction;
  modules.GetHistoryFields().GetItemByKey("RMS_DENSITY").value    = rms_density;

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{10,0});

  REQUIRE(modules.GetHistoryFields().CheckKey("REL_RMS_DENSITY"));
  CHECK(modules.GetHistoryFields().GetValueByKey("REL_RMS_DENSITY") == Approx(reduction));

}

