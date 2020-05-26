#pragma once
#include "COutputModule.hpp"


#include "COutFieldManager.hpp"

#include "../filewriter/CParallelDataSorter.hpp"

#include "../../../include/solvers/CSolver.hpp"

class CGeometry;
class CConfig;
class CSolver;

template <typename ...T>
class CModuleManagerBase {

protected:
  CHistoryOutFieldManager historyFieldsAll;
  CVolumeOutFieldManager volumeFieldsAll;

public:
  virtual void LoadData(const T&... args) = 0;

  virtual void LoadVolumeDataAtPoint(const T&... args, CParallelDataSorter* sorter) = 0;
  virtual void LoadSurfaceDataAtVertex(const T&... args, CParallelDataSorter* sorter) = 0;

  COutFieldCollection& GetHistoryFields()  {return historyFieldsAll.GetCollection();}
  COutFieldCollection& GetVolumeFields()   {return volumeFieldsAll.GetCollection();}
};



template<typename... Modules>
class ModuleList : public std::tuple<Modules...>{
public:
  explicit ModuleList(CConfig* config, int nDim): std::tuple<Modules...>(Modules(config, nDim)...){}
};

template<typename ModuleList>
class CModuleManager : public CModuleManagerBase<SolverData, IterationInfo> {

  ModuleList modules;

  COutFieldCollection::InsertionVector surfaceIntegralHistory;
  COutFieldCollection::InsertionVector surfaceIntegralVolume;

public:
  explicit CModuleManager(CConfig* config, int nDim);

  void SetHistoryFields(CConfig* config);

  void SetVolumeFields(CConfig* config);

  void LoadData(const SolverData& solverData, const IterationInfo& iterationInfo) override;

  void IntegrateCoefficients(const SolverData& solverData, const IterationInfo& iterationInfo, unsigned short iMarker, const string& markerName);

  void CommunicateIntegrals(const string& markerName);

  void LoadVolumeDataAtPoint(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter* sorter) override;

  void LoadSurfaceDataAtVertex(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter* sorter) override;

  std::string GetPerSurfaceName(const std::string& fieldName, const std::string& markerName){
    return fieldName + "@" + markerName;
  }
};

template<typename ModuleList>
CModuleManager<ModuleList>::CModuleManager(CConfig* config, int nDim) : modules(config, nDim) {

  /*--- Set the history output fields for all modules ---*/

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionDefineHistoryFields, historyFieldsAll);

  /*--- Set the volume output fields for all modules ---*/

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionDefineVolumeFields,volumeFieldsAll);

  SetHistoryFields(config);

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionDefineHistoryFieldModifier, historyFieldsAll);

}

template<typename ModuleList>
void CModuleManager<ModuleList>::SetHistoryFields(CConfig *config){

  /*--- Loop through all volume output fields which are surface integrals
   * and add them to the history field collection ---*/

  surfaceIntegralVolume = volumeFieldsAll.GetCollection().GetFieldsByType({FieldType::SURFACE_INTEGRATE});
  surfaceIntegralHistory.clear();

  for (const auto& field : surfaceIntegralVolume){
    auto newField = historyFieldsAll.AddField(field->first, field->second.fieldName,
                                              ScreenOutputFormat::SCIENTIFIC,
                                              field->second.outputGroup,
                                              "Integrated volume field (total)", FieldType::COEFFICIENT);
    surfaceIntegralHistory.push_back(newField);
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

template<typename ModuleList>
void CModuleManager<ModuleList>::LoadData(const SolverData& solverData, const IterationInfo& iterationInfo){

//  SolverDataContainer* solverData = dynamic_cast<SolverDataContainer*>(data);

  const auto& coefficients      = historyFieldsAll.GetCollection().GetFieldsByType({FieldType::COEFFICIENT});

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadHistoryData, historyFieldsAll, solverData, iterationInfo);

  const auto* config = std::get<0>(solverData);

  for (int iMarker_CfgFile = 0; iMarker_CfgFile < config->GetnMarker_CfgFile(); iMarker_CfgFile++){
    const string markerNameCfg =  config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);
    for (int iMarker = 0; iMarker <  config->GetnMarker_All(); iMarker++) {
      const string markerName =  config->GetMarker_All_TagBound(iMarker);
      if (markerName == markerNameCfg){

        IntegrateCoefficients(solverData, iterationInfo, iMarker, markerName);

        CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadHistoryDataPerSurface, historyFieldsAll, solverData, iterationInfo);

        for (const auto& field : coefficients){
          historyFieldsAll.GetCollection().GetItemByKey(GetPerSurfaceName(field->first, markerName)).value = field->second.value;
        }
      }
    }
    CommunicateIntegrals(markerNameCfg);
  }

  for (const auto& field : surfaceIntegralHistory){
    field->second.value = 0.0;
    for (int iMarker_CfgFile = 0; iMarker_CfgFile <  config->GetnMarker_CfgFile(); iMarker_CfgFile++){
      const string markerNameCfg = config->GetMarker_CfgFile_TagBound(iMarker_CfgFile);
      field->second.value += historyFieldsAll.GetCollection().GetValueByKey(GetPerSurfaceName(field->first, markerNameCfg));
    }
  }

  CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadHistoryDataModifier, historyFieldsAll, solverData, iterationInfo);

}


template<typename ModuleList>
void CModuleManager<ModuleList>::LoadVolumeDataAtPoint(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter *sorter){

  volumeFieldsAll.GetCollection().SetCaching(true);

  const auto* geometry = get<1>(solverData);

  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++){

    volumeFieldsAll.GetCollection().StartCaching();

    CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                  solverData, iterationInfo, make_tuple(iPoint, 0, 0));

    if (sorter != nullptr){
      for (const auto& field : volumeFieldsAll.GetCollection().GetReferencesAll()){
        if (field->second.offset != -1){
          sorter->SetUnsorted_Data(iPoint, field->second.offset, field->second.value);
        }
      }
    }
  }
  volumeFieldsAll.GetCollection().SetCaching(false);

}

template<typename ModuleList>
void CModuleManager<ModuleList>::LoadSurfaceDataAtVertex(const SolverData& solverData, const IterationInfo& iterationInfo, CParallelDataSorter *sorter){

  const auto* config = std::get<0>(solverData);
  const auto* geometry = get<1>(solverData);

  volumeFieldsAll.GetCollection().SetCaching(true);

  for (int iMarker = 0; iMarker <  config->GetnMarker_All(); iMarker++) {
    for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++){

      unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

      volumeFieldsAll.GetCollection().StartCaching();

      CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                    solverData, iterationInfo, make_tuple(iPoint, iVertex, iMarker));

      CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadSurfaceData, volumeFieldsAll,
                                    solverData, iterationInfo, make_tuple(iPoint, iVertex, iMarker));

      if (sorter != nullptr){
        for (const auto& field : volumeFieldsAll.GetCollection().GetReferencesAll()){
          if (field->second.offset != -1){
            sorter->SetUnsorted_Data(iPoint, field->second.offset, field->second.value);
          }
        }
      }
    }
  }
  volumeFieldsAll.GetCollection().SetCaching(false);

}
template<typename ModuleList>
void CModuleManager<ModuleList>::IntegrateCoefficients(const SolverData& solverData, const IterationInfo& iterationInfo,
                                                       unsigned short iMarker, const string& markerName){

  const auto* geometry = get<1>(solverData);

  for (const auto& field : surfaceIntegralHistory){ field->second.value = 0.0; }
  volumeFieldsAll.GetCollection().SetCaching(true);
  historyFieldsAll.GetCollection().SetCaching(true);
  for (auto iVertex = 0ul; iVertex < geometry->GetnVertex(iMarker); iVertex++){
    volumeFieldsAll.GetCollection().StartCaching();
    historyFieldsAll.GetCollection().StartCaching();
    unsigned long iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

    CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadVolumeData, volumeFieldsAll,
                                  solverData, iterationInfo, make_tuple(iPoint, iVertex, iMarker));

    /*--- Fill all volume fields with the values at the current surface point --- */

    CSolverOutputModule::for_each(modules, CSolverOutputModule::ActionLoadSurfaceData, volumeFieldsAll,
                                  solverData, iterationInfo, make_tuple(iPoint, iVertex, iMarker));

    /*--- For all volume fields marked as integral, add the contribution to the corresponding history field ---*/

    for (const auto& field : surfaceIntegralVolume){
      historyFieldsAll.GetCollection().GetItemByKey(field->first).value += field->second.value;
    }

  }
  volumeFieldsAll.GetCollection().SetCaching(false);
  historyFieldsAll.GetCollection().SetCaching(false);

  for (const auto& field : surfaceIntegralHistory){
    historyFieldsAll.GetCollection().GetItemByKey(GetPerSurfaceName(field->first, markerName)).value = field->second.value;
  }

}

template <typename ModuleList>
void CModuleManager<ModuleList>::CommunicateIntegrals(const string& markerName){
  std::vector<su2double> sendBuffer(surfaceIntegralHistory.size());
  std::vector<su2double> recvBuffer(surfaceIntegralHistory.size());

  for (unsigned int iField = 0; iField < surfaceIntegralHistory.size(); iField++){
    sendBuffer[iField] = surfaceIntegralHistory[iField]->second.value;
  }

  SU2_MPI::Allreduce(sendBuffer.data(), recvBuffer.data(), surfaceIntegralHistory.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (unsigned int iField = 0; iField < surfaceIntegralHistory.size(); iField++){
    surfaceIntegralHistory[iField]->second.value = recvBuffer[iField];
  }

  for (const auto& field : surfaceIntegralHistory){
    historyFieldsAll.GetCollection().SetValueByKey(GetPerSurfaceName(field->first, markerName), field->second.value);
  }
}
