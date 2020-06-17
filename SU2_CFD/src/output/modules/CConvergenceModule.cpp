#include "../../../include/output/modules/CConvergenceModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/toolboxes/printing_toolbox.hpp"

CConvergenceModule::CConvergenceModule(CConfig* config, int nDim) : CModifierModule(nDim){

  /*--- Initialize convergence monitoring structure ---*/

  nCauchy_Elems = config->GetCauchy_Elems();
  cauchyEps = config->GetCauchy_Eps();
  minLogResidual = config->GetMinLogResidual();
  convStartIter = config->GetStartConv_Iter();

  for (unsigned short iField = 0; iField < config->GetnConv_Field(); iField++){
    convFields.emplace_back(config->GetConv_Field(iField));
  }

  newFunc = vector<su2double>(convFields.size());
  oldFunc = vector<su2double>(convFields.size());
  cauchySerie = vector<vector<su2double>>(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
  cauchyValue = 0.0;
  convergence = false;
  partConvNew.resize(convFields.size());
  partConvOld.resize(convFields.size());

  if (nCauchy_Elems > 1000){
    SU2_MPI::Error("Number of Cauchy Elems must be smaller than 1000", CURRENT_FUNCTION);
  }
}


void CConvergenceModule::DefineHistoryFieldModifier(CHistoryOutFieldManager &historyFields){

  newFunc = vector<su2double>(convFields.size());
  oldFunc = vector<su2double>(convFields.size());
  cauchySerie = vector<vector<su2double>>(convFields.size(), vector<su2double>(nCauchy_Elems, 0.0));
  cauchyValue = 0.0;
  convergence = false;

  historyFields.AddField("CONVERGENCE", "Convergence", ScreenOutputFormat::INTEGER,
                         "CONVERGENCE", "Convergence indicator", FieldType::DEFAULT);

  for (unsigned short iField_Conv = 0; iField_Conv < convFields.size(); iField_Conv++){
    const auto &convField = convFields[iField_Conv];
    if (!historyFields.GetCollection().CheckKey(convField)){
      convFields.erase(convFields.begin() + iField_Conv);
    }
  }

  /*--- Filter convergence fields which are coefficients ---*/

  const auto& convergenceFields = COutFieldCollection::GetFieldsByType({FieldType::COEFFICIENT, FieldType::AUTO_COEFFICIENT, FieldType::PER_SURFACE_COEFFICIENT},
                                                                       historyFields.GetCollection().GetFieldsByKey(convFields));
  for (const auto& field : convergenceFields){

    historyFields.AddField("CAUCHY_" + field->first, "Cauchy["  + field->second.fieldName + "]",
                           ScreenOutputFormat::SCIENTIFIC,
                           "CAUCHY", "Cauchy residual value of field set with CONV_FIELD.", FieldType::DEFAULT);
  }
}

void CConvergenceModule::LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields,
                                                 const IterationInfo& iterationInfo){

  const auto Iter = iterationInfo.Iter;
  unsigned short iCounter;

  convergence = true;
  const auto& outField = historyFields.GetCollection().GetFieldsByKey(convFields);

  int iField_Conv = 0;

  for (const auto& field : outField){

    bool fieldConverged = false;

    su2double monitor = field->second.value;

    /*--- Cauchy based convergence criteria ---*/

    if (field->second.fieldType == FieldType::COEFFICIENT ||
        field->second.fieldType == FieldType::AUTO_COEFFICIENT ||
        field->second.fieldType == FieldType::PER_SURFACE_COEFFICIENT) {

      if (Iter == 0){
        for (iCounter = 0; iCounter < nCauchy_Elems; iCounter++){
          cauchySerie[iField_Conv][iCounter] = 0.0;
        }
        newFunc[iField_Conv] = monitor;
      }

      oldFunc[iField_Conv] = newFunc[iField_Conv];
      newFunc[iField_Conv] = monitor;
      cauchyFunc = fabs(newFunc[iField_Conv] - oldFunc[iField_Conv])/fabs(monitor);

      cauchySerie[iField_Conv][Iter % nCauchy_Elems] = cauchyFunc;
      cauchyValue = 0.0;
      for (iCounter = 0; iCounter < nCauchy_Elems; iCounter++)
        cauchyValue += cauchySerie[iField_Conv][iCounter];

      cauchyValue /= nCauchy_Elems;

      if (cauchyValue >= cauchyEps) { fieldConverged = false;}
      else { fieldConverged = true;}

      /*--- Start monitoring only if the current iteration
         *  is larger than the number of cauchy elements and
         * the number of start-up iterations --- */

      if (Iter < max(convStartIter, nCauchy_Elems)){
        fieldConverged = false;
      }

      historyFields.SetFieldValue("CAUCHY_" + field->first, cauchyValue);

      if(Iter == 0){
        historyFields.SetFieldValue("CAUCHY_" + field->first, 1.0);
      }
    }


    /*--- Residual based convergence criteria ---*/

    if (field->second.fieldType == FieldType::RESIDUAL ||
        field->second.fieldType == FieldType::AUTO_RESIDUAL) {

      /*--- Check the convergence ---*/

      if (Iter != 0 && (monitor <= minLogResidual)) { fieldConverged = true;  }
      else { fieldConverged = false; }

    }

    /*--- Do not apply any convergence criteria of the number
     of iterations is less than a particular value ---*/

    if (Iter < convStartIter) {
      fieldConverged = false;
    }

    partConvNew[iField_Conv] = fieldConverged;

    convergence = fieldConverged && convergence;

    iField_Conv++;
  }

  if (convFields.empty()) convergence = false;

  /*--- Apply the same convergence criteria to all the processors ---*/

  unsigned int sbuf_conv = 0, rbuf_conv = 0;

  /*--- Convergence criteria ---*/

  sbuf_conv = convergence;
  SU2_MPI::Reduce(&sbuf_conv, &rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);

  /*-- Compute global convergence criteria in the master node --*/

  sbuf_conv = 0;
  if (SU2_MPI::GetRank() == MASTER_NODE) {
    if (rbuf_conv == static_cast<unsigned int>(SU2_MPI::GetSize())) sbuf_conv = 1;
    else sbuf_conv = 0;
  }

  SU2_MPI::Bcast(&sbuf_conv, 1, MPI_UNSIGNED_SHORT, MASTER_NODE, MPI_COMM_WORLD);

  if (sbuf_conv == 1) { convergence = true; }
  else { convergence = false;  }

  historyFields.SetFieldValue("CONVERGENCE", convergence);

}

void CConvergenceModule::PrintToStream(std::ostream *stream, const CHistoryOutFieldManager& historyOutFields){

  if ((partConvNew != partConvOld) || convergence){
    PrintingToolbox::CTablePrinter  ConvSummary(stream);

    ConvSummary.AddColumn("Convergence Field", 28);
    ConvSummary.AddColumn("Value", 14);
    ConvSummary.AddColumn("Criterion", 14);
    ConvSummary.AddColumn("Converged",12);
    ConvSummary.SetAlign(PrintingToolbox::CTablePrinter::CENTER);
    ConvSummary.PrintHeader();

    const auto& convFieldRef = historyOutFields.GetCollection().GetFieldsByKey(convFields);

    for (const auto& field : COutFieldCollection::GetFieldsByType({FieldType::COEFFICIENT}, convFieldRef)){
      ConvSummary << field->second.fieldName
                  << field->second.value
                  << " < " + PrintingToolbox::to_string(cauchyEps)
                  << (field->second.value < cauchyEps ? "Yes" : "No");
    }

    for (const auto& field : COutFieldCollection::GetFieldsByType({FieldType::RESIDUAL,
                                                                  FieldType::AUTO_RESIDUAL}, convFieldRef)){
      ConvSummary << field->second.fieldName
                  << field->second.value
                  << " < " + PrintingToolbox::to_string(minLogResidual)
                  << (field->second.value < minLogResidual ? "Yes" : "No");
    }

    partConvOld = partConvNew;

    ConvSummary.PrintFooter();
  }
}