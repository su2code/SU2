#include "../../../include/output/modules/CTimeConvergenceModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

constexpr static char avgPrefix[] = "TAVG_";

CTimeConvergenceModule::CTimeConvergenceModule(CConfig* config) : CSolverOutputModule(config->GetTime_Domain()){
  /*--- Initialize time convergence monitoring structure ---*/

  nWndCauchy_Elems = config->GetWnd_Cauchy_Elems();
  wndCauchyEps     = config->GetWnd_Cauchy_Eps();

  wndConvFields.reserve(config->GetnWndConv_Field());
  for (unsigned short iField = 0; iField < config->GetnWndConv_Field(); iField++){
    wndConvFields.emplace_back(config->GetWndConv_Field(iField));
  }

  WndOld_Func = vector<su2double>(wndConvFields.size());
  WndNew_Func = vector<su2double>(wndConvFields.size());
  WndCauchy_Serie = vector<vector<su2double>>(wndConvFields.size(), vector<su2double>(nWndCauchy_Elems, 0.0));
  WndCauchy_Value = 0.0;
  TimeConvergence = false;

  /*--- Check that the number of cauchy elems is not too large ---*/

  if (nWndCauchy_Elems > 1000){
    SU2_MPI::Error("Number of Time Cauchy Elems must be smaller than 1000", CURRENT_FUNCTION);
  }

}

void CTimeConvergenceModule::DefineHistoryFields(COutFieldCollection &fieldCollection){

  fieldCollection.AddItem("TIME_CONVERGENCE", COutputField("Time Convergence", ScreenOutputFormat::INTEGER,
                                                            "CONVERGENCE", "Time Convergence indicator", FieldType::DEFAULT));
}

void CTimeConvergenceModule::DefineHistoryFieldModifier(COutFieldCollection &fieldCollection){

  const auto& coefficentFields = fieldCollection.GetFieldsByType({FieldType::COEFFICIENT});

  for (auto field : coefficentFields){
    fieldCollection.AddItem(avgPrefix + field->first, COutputField("tavg["  + field->second.fieldName + "]",
                            field->second.screenFormat, avgPrefix   + field->second.outputGroup, "Time averaged values.", FieldType::AUTO_COEFFICIENT));
  }

  const auto& wndConvergenceFields = fieldCollection.GetFieldsByKey(wndConvFields);

  for (const auto& field : wndConvergenceFields){
    fieldCollection.AddItem("CAUCHY_" + field->first, COutputField("Cauchy["  + field->second.fieldName + "]",
                                                                   ScreenOutputFormat::SCIENTIFIC, "CAUCHY",
                                                                   "Cauchy residual value of field set with WND_CONV_FIELD." ,
                                                                   FieldType::DEFAULT));
  }
}

void CTimeConvergenceModule::LoadHistoryData(COutFieldCollection &fieldCollection){

  bool Inner_IterConv = fieldCollection.GetValueByKey("CONVERGENCE") ||
     solverData.config->GetnInner_Iter()-1 <=  solverData.Iter; //Check, if Inner_Iter is converged

  if (Inner_IterConv){
    for (const auto& field : fieldCollection.GetFieldsByType({FieldType::COEFFICIENT})){
      windowedTimeAverages[field->first].addValue(field->second.value, solverData.TimeIter,  solverData.config->GetStartWindowIteration());
      fieldCollection.SetValueByKey(avgPrefix + field->first, windowedTimeAverages[field->first].WindowedUpdate(solverData.config->GetKindWindow()));
    }
  }

  if(solverData.TimeIter == 0){
    for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
      const string WndConv_Field= wndConvFields[iField_Conv];
      if (fieldCollection.CheckKey(WndConv_Field)){
        fieldCollection.SetValueByKey("CAUCHY_"+ WndConv_Field, 1.0);
      }
    }
  }

  if(Inner_IterConv && solverData.TimeIter >= solverData.config->GetStartWindowIteration()){
    TimeConvergence = true;
    unsigned short iCounter;

    for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
      const string WndConv_Field= wndConvFields[iField_Conv];

      if (fieldCollection.CheckKey(WndConv_Field)){
        bool fieldConverged = false;

        su2double monitor = fieldCollection.GetValueByKey(WndConv_Field);

        /*--- Cauchy based convergence criteria ---*/

        if (fieldCollection.GetItemByKey(WndConv_Field).fieldType == FieldType::AUTO_COEFFICIENT) { //TAVG values are AUTO_COEFF
          if (solverData.TimeIter == solverData.config->GetStartWindowIteration()){
            for (iCounter = 0; iCounter < nWndCauchy_Elems; iCounter++){
              WndCauchy_Serie[iField_Conv][iCounter] = 0.0;
            }
            WndNew_Func[iField_Conv] = monitor;
          }
          WndOld_Func[iField_Conv] = WndNew_Func[iField_Conv];
          WndNew_Func[iField_Conv] = monitor;
          WndCauchy_Func = fabs(WndNew_Func[iField_Conv] - WndOld_Func[iField_Conv]);
          WndCauchy_Serie[iField_Conv][solverData.TimeIter % nWndCauchy_Elems] = WndCauchy_Func;
          WndCauchy_Value = 1.0;

          if (solverData.TimeIter >= nWndCauchy_Elems+solverData.config->GetStartWindowIteration()){
            WndCauchy_Value = 0.0;
            for (iCounter = 0; iCounter < nWndCauchy_Elems; iCounter++){
              WndCauchy_Value += WndCauchy_Serie[iField_Conv][iCounter];
            }
            WndCauchy_Value /= nWndCauchy_Elems;
          }
          if (WndCauchy_Value >= wndCauchyEps){fieldConverged = false;}
          else{fieldConverged = true;}

          /*--- Start monitoring only if the current iteration is larger than the
           *  number of cauchy elements and the number of start-up iterations ---*/

          if (solverData.TimeIter <  solverData.config->GetStartWindowIteration() + max(solverData.config->GetWnd_StartConv_Iter(), nWndCauchy_Elems)){
            fieldConverged = false;
          }
          fieldCollection.SetValueByKey("CAUCHY_" + WndConv_Field, WndCauchy_Value);
        }
        TimeConvergence = fieldConverged && TimeConvergence;

        /*--- Stop the simulation in case a nan appears, do not save the solution ---*/

        if (monitor != monitor){
          SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);}
      }
    }

    /*--- Do not apply any convergence criterion if the option is disabled. */
    if(!solverData.config->GetWnd_Cauchy_Crit()){TimeConvergence = false;}
    if(wndConvFields.empty()){TimeConvergence = false;}
  }
  fieldCollection.SetValueByKey("TIME_CONVERGENCE", TimeConvergence);

}