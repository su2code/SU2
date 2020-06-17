#pragma once
#include "../fields/COutFieldManager.hpp"

class CConfig;
class CParallelDataSorter;

template <typename ...T>
class CModuleManagerBase {

protected:
  CHistoryOutFieldManager historyFieldsAll;
  CVolumeOutFieldManager volumeFieldsAll;

public:
  virtual void LoadData(const T&... args) = 0;

  virtual void LoadVolumeDataAtPoint(const T&... args, CParallelDataSorter* sorter) = 0;
  virtual void LoadSurfaceDataAtVertex(const T&... args, CParallelDataSorter* sorter) = 0;

  virtual void PrintToStream(std::ostream* stream) = 0;

  CHistoryOutFieldManager& GetHistoryFields() {return historyFieldsAll;}
  CVolumeOutFieldManager& GetVolumeFields() {return volumeFieldsAll;}

  void Clear(){
    historyFieldsAll = CHistoryOutFieldManager();
    volumeFieldsAll  = CVolumeOutFieldManager();
  }
};


template<typename... Modules>
class ModuleList : public std::tuple<Modules...>{
public:
  ModuleList(CConfig* config, int nDim): std::tuple<Modules...>(Modules(config, nDim)...){}
};