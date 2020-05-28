#include "../../../include/output/modules/CCommonModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

void CCommonModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields){

  historyFields.AddField("TIME_ITER", "Time iter", ScreenOutputFormat::INTEGER,
                         "ITER", "Time iteration index", FieldType::DEFAULT);
  historyFields.AddField("OUTER_ITER", "Outer iter", ScreenOutputFormat::INTEGER,
                         "ITER", "Outer iteration index", FieldType::DEFAULT);
  historyFields.AddField("INNER_ITER", "Inner iter", ScreenOutputFormat::INTEGER,
                         "ITER", "Inner iteration index", FieldType::DEFAULT);

  historyFields.AddField("CUR_TIME", "Cur_time", ScreenOutputFormat::SCIENTIFIC,
                         "TIME_DOMAIN", "Current physical time (s)", FieldType::DEFAULT);
  historyFields.AddField("TIME_STEP", "Time_Step", ScreenOutputFormat::SCIENTIFIC,
                         "TIME_DOMAIN", "Current time step (s)", FieldType::DEFAULT);

  historyFields.AddField("WALL_TIME", "Time(sec)", ScreenOutputFormat::SCIENTIFIC,
                         "WALL_TIME", "Average wall-clock time", FieldType::DEFAULT);

  historyFields.AddField("NONPHYSICAL_POINTS", "Time(sec)", ScreenOutputFormat::INTEGER,
                         "NONPHYSICAL_POINTS", "The number of non-physical points in the solution", FieldType::DEFAULT);

}

void CCommonModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                    const IterationInfo& iterationInfo){

  const auto* config = solverData.config;
  const auto curInnerIter = iterationInfo.Iter;
  const auto curTimeIter  = iterationInfo.TimeIter;

  su2double TimeStep = config->GetDelta_UnstTimeND()*config->GetTime_Ref();
  su2double curTime = historyFields.GetFieldValue("CUR_TIME");

  historyFields.SetFieldValue("TIME_STEP", TimeStep);

   if (SU2_TYPE::Int(historyFields.GetFieldValue("TIME_ITER")) != static_cast<int>(curTimeIter)) {
     historyFields.SetFieldValue("CUR_TIME", curTime + TimeStep);
   }

   historyFields.SetFieldValue("TIME_ITER", curTimeIter);
   historyFields.SetFieldValue("INNER_ITER", curInnerIter);

   su2double StopTime, UsedTime;

   StopTime = SU2_MPI::Wtime();

   UsedTime = (StopTime - config->Get_StartTime())/(curInnerIter+1);

   historyFields.SetFieldValue("WALL_TIME", UsedTime);
   historyFields.SetFieldValue("NONPHYSICAL_POINTS", config->GetNonphysical_Points());
}