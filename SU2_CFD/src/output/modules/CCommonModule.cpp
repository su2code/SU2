#include "../../../include/output/modules/CCommonModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

void CCommonModule::DefineHistoryFields(COutFieldCollection &fieldCollection){


  fieldCollection.AddItem("TIME_ITER", COutputField("Time iter", ScreenOutputFormat::INTEGER,
                                                          "ITER", FieldType::DEFAULT, "Time iteration index"));
  fieldCollection.AddItem("OUTER_ITER", COutputField("Outer iter", ScreenOutputFormat::INTEGER,
                                                           "ITER", FieldType::DEFAULT, "Outer iteration index"));
  fieldCollection.AddItem("INNER_ITER", COutputField("Inner iter", ScreenOutputFormat::INTEGER,
                                                           "ITER", FieldType::DEFAULT, "Inner iteration index"));

  fieldCollection.AddItem("CUR_TIME", COutputField("Cur_time", ScreenOutputFormat::SCIENTIFIC,
                                                           "TIME_DOMAIN", FieldType::DEFAULT, "Current physical time (s)"));
  fieldCollection.AddItem("TIME_STEP", COutputField("Time_Step", ScreenOutputFormat::SCIENTIFIC,
                                                           "TIME_DOMAIN", FieldType::DEFAULT, "Current time step (s)"));

  fieldCollection.AddItem("WALL_TIME", COutputField("Time(sec)", ScreenOutputFormat::SCIENTIFIC,
                                                          "WALL_TIME", FieldType::DEFAULT, "Average wall-clock time"));

  fieldCollection.AddItem("NONPHYSICAL_POINTS", COutputField("Time(sec)", ScreenOutputFormat::INTEGER,
                                                                   "NONPHYSICAL_POINTS", FieldType::DEFAULT,
                                                                    "The number of non-physical points in the solution"));

}

void CCommonModule::LoadHistoryData(COutFieldCollection &fieldCollection){

  unsigned long curOuterIter = solverData.config->GetOuterIter();
  unsigned long curInnerIter = solverData.config->GetInnerIter();
  su2double TimeStep = solverData.config->GetDelta_UnstTimeND()*solverData.config->GetTime_Ref();
  su2double curTime = fieldCollection.GetValueByKey("CUR_TIME");

  fieldCollection.SetValueByKey("TIME_STEP", TimeStep);

   if (SU2_TYPE::Int(TimeStep) != static_cast<int>(solverData.config->GetTimeIter())) {
     fieldCollection.SetValueByKey("CUR_TIME", curTime*TimeStep);
   }

   fieldCollection.SetValueByKey("TIME_ITER", solverData.config->GetTimeIter());
   fieldCollection.SetValueByKey("OUTER_ITER", curOuterIter);
   fieldCollection.SetValueByKey("INNER_ITER", curInnerIter);

   su2double StopTime, UsedTime;

   StopTime = SU2_MPI::Wtime();

   UsedTime = (StopTime - solverData.config->Get_StartTime())/((curOuterIter + 1) * (curInnerIter+1));

   fieldCollection.SetValueByKey("WALL_TIME", UsedTime);
   fieldCollection.SetValueByKey("NONPHYSICAL_POINTS", solverData.config->GetNonphysical_Points());
}