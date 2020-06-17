#include "../../../include/output/modules/CTimeConvergenceModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

constexpr static char avgPrefix[] = "TAVG_";

CTimeConvergenceModule::CTimeConvergenceModule(CConfig* config, int nDim) : CModifierModule(nDim, config->GetTime_Domain()){
  /*--- Initialize time convergence monitoring structure ---*/

  nWndCauchy_Elems = config->GetWnd_Cauchy_Elems();
  wndCauchyEps     = config->GetWnd_Cauchy_Eps();
  startWindowIteration = config->GetStartWindowIteration();
  kindWindow = config->GetKindWindow();
  startWindowConvIteration = config->GetWnd_StartConv_Iter();
  windowCauchyCrit = config->GetWnd_Cauchy_Crit();
  wndConvFields.reserve(config->GetnWndConv_Field());

  nIter = config->GetnInner_Iter();

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


void CTimeConvergenceModule::DefineHistoryFieldModifier(CHistoryOutFieldManager &historyFields){

  WndOld_Func = vector<su2double>(wndConvFields.size());
  WndNew_Func = vector<su2double>(wndConvFields.size());
  WndCauchy_Serie = vector<vector<su2double>>(wndConvFields.size(), vector<su2double>(nWndCauchy_Elems, 0.0));
  WndCauchy_Value = 0.0;
  TimeConvergence = false;

  historyFields.AddField("TIME_CONVERGENCE", "Time Convergence", ScreenOutputFormat::INTEGER,
                         "CONVERGENCE", "Time Convergence indicator", FieldType::DEFAULT);

  const auto& coefficentFields = historyFields.GetCollection().GetFieldsByType({FieldType::COEFFICIENT,
                                                                                FieldType::AUTO_COEFFICIENT,
                                                                                FieldType::PER_SURFACE_COEFFICIENT});

  for (auto field : coefficentFields){
    historyFields.AddField(avgPrefix + field->first, "tavg["  + field->second.fieldName + "]",
                           field->second.screenFormat, avgPrefix   + field->second.outputGroup, "Time averaged values.", FieldType::AUTO_COEFFICIENT);
  }

  const auto& wndConvergenceFields = historyFields.GetCollection().GetFieldsByKey(wndConvFields);

  for (const auto& field : wndConvergenceFields){
    if (field->first.substr(0,4) != "TAVG"){
      SU2_MPI::Error("Option values for CONV_WINDOW_FIELD must be time-averaged quantities", CURRENT_FUNCTION);
    }
    historyFields.AddField("CAUCHY_" + field->first,"Cauchy["  + field->second.fieldName + "]",
                           ScreenOutputFormat::SCIENTIFIC, "CAUCHY",
                           "Cauchy residual value of field set with WND_CONV_FIELD." ,
                           FieldType::DEFAULT);
  }
}

void CTimeConvergenceModule::LoadHistoryDataModifier(CHistoryOutFieldManager &historyFields, const IterationInfo &iterationInfo){

  const auto Iter     = iterationInfo.Iter;
  const auto TimeIter = iterationInfo.TimeIter;

  bool Inner_IterConv = historyFields.GetFieldValue("CONVERGENCE") ||
                        nIter-1 <=  Iter; //Check, if Inner_Iter is converged

  if (Inner_IterConv){
    for (const auto& field : historyFields.GetCollection().GetFieldsByType({FieldType::COEFFICIENT, FieldType::PER_SURFACE_COEFFICIENT})){
      windowedTimeAverages[field->first].addValue(field->second.value, TimeIter,  startWindowIteration);
      historyFields.SetFieldValue(avgPrefix + field->first, windowedTimeAverages[field->first].WindowedUpdate(kindWindow));
    }
  }

  if(TimeIter == 0){
    for (unsigned short iField_Conv = 0; iField_Conv < wndConvFields.size(); iField_Conv++){
      const string WndConv_Field= wndConvFields[iField_Conv];
      if (historyFields.GetCollection().CheckKey(WndConv_Field)){
        historyFields.SetFieldValue("CAUCHY_"+ WndConv_Field, 1.0);
      }
    }
  }

  if(Inner_IterConv && TimeIter >= startWindowIteration){
    TimeConvergence = true;
    unsigned short iCounter;

    const auto& convFields = historyFields.GetCollection().GetFieldsByKey(wndConvFields);
    int iField_Conv = 0;
    for (const auto& field : convFields){

      const string WndConv_Field= wndConvFields[iField_Conv];

      bool fieldConverged = false;

      /*--- Cauchy based convergence criteria ---*/

      if (field->second.fieldType == FieldType::AUTO_COEFFICIENT) { //TAVG values are AUTO_COEFF
        if (TimeIter == startWindowIteration){
          for (iCounter = 0; iCounter < nWndCauchy_Elems; iCounter++){
            WndCauchy_Serie[iField_Conv][iCounter] = 0.0;
          }
          WndNew_Func[iField_Conv] = field->second.value;
        }
        WndOld_Func[iField_Conv] = WndNew_Func[iField_Conv];
        WndNew_Func[iField_Conv] = field->second.value;
        WndCauchy_Func = fabs(WndNew_Func[iField_Conv] - WndOld_Func[iField_Conv]);
        WndCauchy_Serie[iField_Conv][TimeIter % nWndCauchy_Elems] = WndCauchy_Func;
        WndCauchy_Value = 1.0;

        if (TimeIter >= nWndCauchy_Elems + startWindowIteration){
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

        if (TimeIter <  startWindowIteration + max(startWindowConvIteration, nWndCauchy_Elems)){
          fieldConverged = false;
        }
        historyFields.SetFieldValue("CAUCHY_" + WndConv_Field, WndCauchy_Value);
      }
      TimeConvergence = fieldConverged && TimeConvergence;

      /*--- Stop the simulation in case a nan appears, do not save the solution ---*/

      if (field->second.value != field->second.value){
        SU2_MPI::Error("SU2 has diverged (NaN detected).", CURRENT_FUNCTION);}

      iField_Conv++;
    }

    /*--- Do not apply any convergence criterion if the option is disabled. */
    if(windowCauchyCrit){TimeConvergence = false;}
    if(wndConvFields.empty()){TimeConvergence = false;}
  }
  historyFields.SetFieldValue("TIME_CONVERGENCE", TimeConvergence);

}