#include "catch.hpp"
#include <sstream>

#include "../../UnitQuadTestCase.hpp"

#include "../../../SU2_CFD/include/solvers/CSolverFactory.hpp"
#include "../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../SU2_CFD/include/solvers/CNSSolver.hpp"
#include "../../../SU2_CFD/include/output/module_manager/CModuleManager.hpp"
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

  CModuleManager<ModuleList<CCommonModule>, ModuleList<>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{100,0});

  CHECK(modules.GetHistoryFields().GetFieldValue("INNER_ITER") == Approx(100));

}

TEST_CASE("FVM Base Module", "[Output Module]"){
  CModuleManager<ModuleList<CFVMBaseModule>, ModuleList<>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  auto Area =   modules.GetHistoryFields().GetFieldValue("AREA@x_minus")
              + modules.GetHistoryFields().GetFieldValue("AREA@x_plus")
              + modules.GetHistoryFields().GetFieldValue("AREA@y_minus")
              + modules.GetHistoryFields().GetFieldValue("AREA@y_plus")
              + modules.GetHistoryFields().GetFieldValue("AREA@z_minus")
              + modules.GetHistoryFields().GetFieldValue("AREA@z_plus");

  CHECK(modules.GetHistoryFields().GetFieldValue("AREA") == Approx(Area));

}

TEST_CASE("Aerodynamics Module", "[Output Module]"){

  CModuleManager<ModuleList<CAerodynamicsModule>, ModuleList<>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  ModuleTestCase->solver[FLOW_SOL]->Preprocessing(ModuleTestCase->geometry.get(), ModuleTestCase->solver, ModuleTestCase->config.get(), 0, 0, RUNTIME_FLOW_SYS, false);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});


  CHECK(modules.GetHistoryFields().GetFieldValue("FORCE_Y") == modules.GetHistoryFields().GetFieldValue("SIDEFORCE"));
  CHECK(modules.GetHistoryFields().GetFieldValue("SIDEFORCE") == Approx(-1.9387861639));
  CHECK(modules.GetHistoryFields().GetFieldValue("SIDEFORCE") == modules.GetHistoryFields().GetFieldValue("SIDEFORCE@y_minus") +
        modules.GetHistoryFields().GetFieldValue("SIDEFORCE@y_plus"));
}

TEST_CASE("Convergence Module", "[Output Module]"){

  CModuleManager<ModuleList<>, ModuleList<CConvergenceModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.Clear();

  modules.GetHistoryFields().AddField("RMS_DENSITY", "", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);
  modules.GetHistoryFields().AddField("DRAG", "", ScreenOutputFormat::SCIENTIFIC, "DRAG", "", FieldType::AUTO_COEFFICIENT);

  modules.Init(ModuleTestCase->config.get());

  REQUIRE(modules.GetHistoryFields().GetCollection().CheckKey("CAUCHY_DRAG"));

  modules.GetHistoryFields().SetFieldValue("RMS_DENSITY", -15.0);
  su2double val1 = 1.001;
  modules.GetHistoryFields().SetFieldValue("DRAG", val1);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  CHECK(modules.GetHistoryFields().GetFieldValue("CAUCHY_DRAG") == 1.0);
  su2double val2 = 1.002;
  modules.GetHistoryFields().SetFieldValue("DRAG", 1.002);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{1,0});

  CHECK(modules.GetHistoryFields().GetFieldValue("CAUCHY_DRAG") == Approx((0 + 0 + abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetFieldValue("CONVERGENCE") == 0.0);

  su2double val3 = 1.003;
  modules.GetHistoryFields().SetFieldValue("DRAG", val3);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{50,0});

  su2double val4 = 1.004;
  modules.GetHistoryFields().SetFieldValue("DRAG", val4);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{51,0});

  CHECK(modules.GetHistoryFields().GetFieldValue("CAUCHY_DRAG") == Approx((abs(val4-val3)/abs(val4) +
                                                                           abs(val3-val2)/abs(val3) +
                                                                           abs(val2-val1)/abs(val2))/3));

  CHECK(modules.GetHistoryFields().GetFieldValue("CONVERGENCE") == 1.0);

}

TEST_CASE("Residual Module", "[Residual Module]"){

  CModuleManager<ModuleList<>, ModuleList<CResidualModule>>  modules(ModuleTestCase->config.get(), ModuleTestCase->geometry->GetnDim());

  modules.Clear();

  modules.GetHistoryFields().AddField("RMS_DENSITY", "", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);
  modules.GetHistoryFields().AddField("RMS_MOMENTUM_X","", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);
  modules.GetHistoryFields().AddField("RMS_MOMENTUM_Y", "", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);
  modules.GetHistoryFields().AddField("RMS_MOMENTUM_Z", "", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);
  modules.GetHistoryFields().AddField("RMS_ENERGY", "", ScreenOutputFormat::SCIENTIFIC, "RMS_RES", "", FieldType::RESIDUAL);

  su2double rms_density = -5.765;
  modules.GetHistoryFields().SetFieldValue("RMS_DENSITY", rms_density);
  su2double rms_mom_x   = -6.52;
  modules.GetHistoryFields().SetFieldValue("RMS_MOMENTUM_X", rms_mom_x);
  su2double rms_mom_y   = -8.52;
  modules.GetHistoryFields().SetFieldValue("RMS_MOMENTUM_Y", rms_mom_y);
  su2double rms_mom_z   = -3.32;
  modules.GetHistoryFields().SetFieldValue("RMS_MOMENTUM_Z", rms_mom_z);
  su2double rms_energy  = -4.64;
  modules.GetHistoryFields().SetFieldValue("RMS_ENERGY", rms_energy);

  modules.Init(ModuleTestCase->config.get());

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{0,0});

  REQUIRE(modules.GetHistoryFields().GetCollection().CheckKey("AVG_RMS_RES"));
  CHECK(modules.GetHistoryFields().GetFieldValue("AVG_RMS_RES") == (rms_density+rms_mom_x+rms_mom_y+rms_mom_z+rms_energy)/5);

  su2double reduction = -1.453;
  rms_density += reduction;
  modules.GetHistoryFields().SetFieldValue("RMS_DENSITY", rms_density);

  modules.LoadData({ModuleTestCase->config.get(), ModuleTestCase->geometry.get(), ModuleTestCase->solver},{10,0});

  REQUIRE(modules.GetHistoryFields().GetCollection().CheckKey("REL_RMS_DENSITY"));
  CHECK(modules.GetHistoryFields().GetFieldValue("REL_RMS_DENSITY") == Approx(reduction));

}

