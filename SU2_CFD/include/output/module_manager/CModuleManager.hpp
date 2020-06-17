#pragma once
#include "../modules/COutputModule.hpp"
#include "CModuleManagerBase.hpp"

#include "../filewriter/CParallelDataSorter.hpp"
#include "../../../include/solvers/CSolver.hpp"

template<typename ModuleList, typename ModifierModuleList>
class CModuleManager : public CModuleManagerBase<SolverData, IterationInfo> {

  ModuleList modules;
  ModifierModuleList modifierModules;


public:
  explicit CModuleManager(CConfig* config, int nDim);

  void SetHistoryFields(CConfig* config);

  void SetVolumeFields(CConfig* config);

  void LoadData(const SolverData& solverData, const IterationInfo& iterationInfo) override;

  void IntegrateCoefficients(const SolverData& solverData, const IterationInfo& iterationInfo, unsigned short iMarker);

  void CommunicateSurfaceIntegrals();

  void Init(CConfig* config);

  void LoadVolumeDataAtPoint(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter* sorter) override;

  void LoadSurfaceDataAtVertex(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter* sorter) override;

  std::string GetPerSurfaceName(const std::string& fieldName, const std::string& markerName){
    return fieldName + "@" + markerName;
  }

  void PrintToStream(std::ostream *stream) override {
    CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionPrintToStream, stream, historyFieldsAll);
    CModifierModule::for_each(modifierModules, CModifierModule::ActionPrintToStream, stream, historyFieldsAll);
  }
};



template<typename ModuleList, typename ModifierModuleList>
CModuleManager<ModuleList, ModifierModuleList>::CModuleManager(CConfig* config, int nDim) : modules(config, nDim), modifierModules(config, nDim) {

  Init(config);

}

template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::Init(CConfig* config) {

  /*--- Set the history output fields for all modules ---*/

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionDefineHistoryFields, historyFieldsAll);

  /*--- Set the volume output fields for all modules ---*/

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionDefineVolumeFields,volumeFieldsAll);

  CModifierModule::for_each(modifierModules, CModifierModule::ActionDefineVolumeFieldModifier, volumeFieldsAll);

  SetHistoryFields(config);

  CModifierModule::for_each(modifierModules, CModifierModule::ActionDefineHistoryFieldModifier, historyFieldsAll);

}

template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::SetHistoryFields(CConfig *config){

  /*--- Loop through all volume output fields which are surface integrals
   * and add them to the history field collection ---*/

  const auto &surfaceIntegralVolume = volumeFieldsAll.GetCollection().GetFieldsByType({FieldType::SURFACE_INTEGRATE});

  for (const auto& field : surfaceIntegralVolume){
    historyFieldsAll.AddField(field->first, field->second.fieldName,
                              ScreenOutputFormat::SCIENTIFIC,
                              field->second.outputGroup,
                              "Integrated volume field (total)", FieldType::COEFFICIENT);
  }

  std::vector<std::string> markerNames;
  for (int iMarker_CfgFile = 0; iMarker_CfgFile < config->GetnMarker_CfgFile(); iMarker_CfgFile++){

    const std::string markerName = config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);

    markerNames.push_back(markerName);
  }
  /*--- Loop through all coefficients and add a field for each marker ---*/

  const auto& coefficients = historyFieldsAll.GetCollection().GetFieldsByType({FieldType::COEFFICIENT});

  for (const auto& field : coefficients){
    for (const auto& marker : markerNames){
      historyFieldsAll.AddField(GetPerSurfaceName(field->first, marker), GetPerSurfaceName(field->second.fieldName, marker),
                                ScreenOutputFormat::SCIENTIFIC,
                                GetPerSurfaceName(field->second.outputGroup, marker),
                                "Integrated volume field",
                                FieldType::PER_SURFACE_COEFFICIENT);
    }
  }
}

template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::LoadData(const SolverData& solverData, const IterationInfo& iterationInfo){

  const auto& coefficients      = historyFieldsAll.GetCollection().GetFieldsByType({FieldType::COEFFICIENT});

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadHistoryData, historyFieldsAll, solverData, iterationInfo);

  for (const auto& field : historyFieldsAll.GetCollection().GetFieldsByType({FieldType::PER_SURFACE_COEFFICIENT})){ field->second.value = 0.0; }

  const auto* config = solverData.config;

  for (int iMarker_CfgFile = 0; iMarker_CfgFile < config->GetnMarker_CfgFile(); iMarker_CfgFile++){
    const string markerNameCfg =  config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);
    for (int iMarker = 0; iMarker <  config->GetnMarker_All(); iMarker++) {
      const string markerName =  config->GetMarker_All_TagBound(iMarker);
      if (markerName == markerNameCfg){

        IntegrateCoefficients(solverData, iterationInfo, iMarker);

        CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadHistoryDataPerSurface, historyFieldsAll, solverData, iterationInfo);

        for (const auto& field : coefficients){
          historyFieldsAll.GetCollection().GetItemByKey(GetPerSurfaceName(field->first, markerName)).value = field->second.value;
        }
      }
    }
  }

  CommunicateSurfaceIntegrals();

  for (const auto& field : coefficients){
    field->second.value = 0.0;
    for (int iMarker_CfgFile = 0; iMarker_CfgFile <  config->GetnMarker_CfgFile(); iMarker_CfgFile++){
      const string markerNameCfg = config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);
      field->second.value += historyFieldsAll.GetFieldValue(GetPerSurfaceName(field->first, markerNameCfg));
    }
  }

  CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadHistoryDataModifier, historyFieldsAll, iterationInfo);

}


template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::LoadVolumeDataAtPoint(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter *sorter){

  volumeFieldsAll.SetCaching(true);

  const auto* geometry = solverData.geometry;

  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++){

    if (geometry->nodes->GetDomain(iPoint)){

      volumeFieldsAll.StartCaching();

      CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                    solverData, iterationInfo, PointInfo{iPoint, 0, 0});

      CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadVolumeDataModifier, volumeFieldsAll,
                                iterationInfo);

      if (sorter != nullptr){
        for (const auto& field : volumeFieldsAll.GetCollection().GetReferencesAll()){
          if (field->second.offset != -1){
            sorter->SetUnsorted_Data(iPoint, field->second.offset, field->second.value);
          }
        }
      }
    }
  }
  volumeFieldsAll.SetCaching(false);

}

template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::LoadSurfaceDataAtVertex(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter *sorter){

  const auto* config = solverData.config;
  const auto* geometry = solverData.geometry;

  volumeFieldsAll.SetCaching(true);

  for (uint iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++){

      unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      if (geometry->nodes->GetDomain(iPoint)){

        volumeFieldsAll.StartCaching();

        CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                      solverData, iterationInfo, PointInfo{iPoint, iVertex, iMarker});

        CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadVolumeDataModifier, volumeFieldsAll,
                                  iterationInfo);

        CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadSurfaceData, volumeFieldsAll,
                                      solverData, iterationInfo, PointInfo{iPoint, iVertex, iMarker});

        CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadSurfaceDataModifier, volumeFieldsAll,
                                  iterationInfo);

        if (sorter != nullptr){
          for (const auto& field : volumeFieldsAll.GetCollection().GetReferencesAll()){
            if (field->second.offset != -1){
              sorter->SetUnsorted_Data(iPoint, field->second.offset, field->second.value);
            }
          }
        }
      }
    }
  }
  volumeFieldsAll.SetCaching(false);

}
template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::IntegrateCoefficients(const SolverData& solverData, const IterationInfo& iterationInfo,
                                                       unsigned short iMarker){

  const auto* geometry = solverData.geometry;

  const auto& surfaceIntegralFields = volumeFieldsAll.GetCollection().GetFieldsByType({FieldType::SURFACE_INTEGRATE});

  for (const auto& field : historyFieldsAll.GetCollection().GetFieldsByType({FieldType::COEFFICIENT})){ field->second.value = 0.0; }
  volumeFieldsAll.SetCaching(true);
  historyFieldsAll.SetCaching(true);
  for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++){

    unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
    if (geometry->nodes->GetDomain(iPoint)){

      volumeFieldsAll.StartCaching();
      historyFieldsAll.StartCaching();

      CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                    solverData, iterationInfo, PointInfo{iPoint, iVertex, iMarker});

      CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadVolumeDataModifier, volumeFieldsAll,
                                iterationInfo);

      /*--- Fill all volume fields with the values at the current surface point --- */

      CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadSurfaceData, volumeFieldsAll,
                                    solverData, iterationInfo, PointInfo{iPoint, iVertex, iMarker});

      CModifierModule::for_each(modifierModules, CModifierModule::ActionLoadSurfaceDataModifier, volumeFieldsAll,
                                iterationInfo);

      /*--- For all volume fields marked as integral, add the contribution to the corresponding history field ---*/

      for (const auto& field : surfaceIntegralFields){
        historyFieldsAll.SetFieldValue(field->first, historyFieldsAll.GetFieldValue(field->first) + field->second.value);
      }
    }

  }
  volumeFieldsAll.SetCaching(false);
  historyFieldsAll.SetCaching(false);

}

template<typename ModuleList, typename ModifierModuleList>
void CModuleManager<ModuleList, ModifierModuleList>::CommunicateSurfaceIntegrals(){

  const auto& perSurfaceCoefficients = historyFieldsAll.GetCollection().GetFieldsByType({FieldType::PER_SURFACE_COEFFICIENT});

  std::vector<su2double> sendBuffer(perSurfaceCoefficients.size());
  std::vector<su2double> recvBuffer(perSurfaceCoefficients.size());

  for (unsigned int iField = 0; iField < perSurfaceCoefficients.size(); iField++){
    sendBuffer[iField] = perSurfaceCoefficients[iField]->second.value;
  }

  SU2_MPI::Allreduce(sendBuffer.data(), recvBuffer.data(), perSurfaceCoefficients.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned int iField = 0; iField < perSurfaceCoefficients.size(); iField++){
    perSurfaceCoefficients[iField]->second.value = recvBuffer[iField];
  }

}
