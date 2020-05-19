#include "../fields/COutFieldCollection.hpp"
#include <tuple>
#include <utility>
#pragma once
#include "SolverDataContainer.hpp"
#include "COutFieldManager.hpp"

#define CREATE_ACTION(name, func_args...) \
  static struct Action##name {} Action##name; \
  template <typename type, typename... Targs> \
  static void Action(type& Module, struct Action##name, Targs& ...args){ Module.name(args...); }


template<class Module>
class COutputModule {

protected:

  bool enabled;

public:

  COutputModule(bool enabled_ = true) :
    enabled(enabled_){}


  bool IsEnabled() { return enabled; }

  virtual void UpdateData(const OutputData* data) = 0;

  template<std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
  inline static typename std::enable_if<I == sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...> &, TAction, Targs&...) { }

  template<std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
  inline static typename std::enable_if<I < sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...>& t, TAction action, Targs&... args)
  {
    if (std::get<I>(t).IsEnabled()){
      Module::Action(std::get<I>(t), action, args...);
    }
    for_each<I + 1, TAction, Tp...>(t, action, args...);
  }

};

class CSolverOutputModule : public COutputModule<CSolverOutputModule>{

protected:
  SolverDataContainer solverData;
  int nDim;

public:
  CSolverOutputModule(int nDim, bool enabled_ = true) : COutputModule(enabled_), nDim(nDim) {};

  void UpdateData(const OutputData *data) override {
    solverData = *dynamic_cast<const SolverDataContainer*>(data);
  }
  CREATE_ACTION(UpdateData, const OutputData *data);

  virtual void DefineHistoryFields(CHistoryOutFieldManager& historyFields) {};
  CREATE_ACTION(DefineHistoryFields, CHistoryOutFieldManager& historyFields);
  virtual void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) {};
  CREATE_ACTION(DefineHistoryFieldModifier, CHistoryOutFieldManager& historyFields);
  virtual void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {};
  CREATE_ACTION(DefineVolumeFields, CVolumeOutFieldManager& volumeFields);
  virtual void LoadHistoryData(CHistoryOutFieldManager& historyFields) {};
  CREATE_ACTION(LoadHistoryData, CHistoryOutFieldManager& historyFields);
  virtual void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields) {};
  CREATE_ACTION(LoadHistoryDataModifier, CHistoryOutFieldManager& historyFields);
  virtual void LoadHistoryDataPerSurface(CHistoryOutFieldManager& historyFields) {};
  CREATE_ACTION(LoadHistoryDataPerSurface, CHistoryOutFieldManager& historyFields);
  virtual void LoadVolumeData(CVolumeOutFieldManager& volumeFields) {};
  CREATE_ACTION(LoadVolumeData, CVolumeOutFieldManager& volumeFields);
  virtual void LoadSurfaceData(CVolumeOutFieldManager& volumeFields) {};
  CREATE_ACTION(LoadSurfaceData, CVolumeOutFieldManager& volumeFields);


};

#undef CREATE_ACTION

