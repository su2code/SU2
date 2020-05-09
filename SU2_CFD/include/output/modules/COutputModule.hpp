#include "../fields/COutFieldCollection.hpp"
#include <tuple>
#include <utility>
#pragma once
#include "SolverDataContainer.hpp"


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

public:
  CSolverOutputModule(bool enabled_ = true) : COutputModule(enabled_) {};

  void UpdateData(const OutputData *data) override {
    solverData = *dynamic_cast<const SolverDataContainer*>(data);
  }
  CREATE_ACTION(UpdateData, const OutputData *data);

  virtual void DefineHistoryFields(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(DefineHistoryFields, COutFieldCollection& fieldCollection);
  virtual void DefineHistoryFieldModifier(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(DefineHistoryFieldModifier, COutFieldCollection& fieldCollection);
  virtual void DefineVolumeFields(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(DefineVolumeFields, COutFieldCollection& fieldCollection);
  virtual void LoadHistoryData(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(LoadHistoryData, COutFieldCollection& fieldCollection);
  virtual void LoadHistoryDataModifier(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(LoadHistoryDataModifier, COutFieldCollection& fieldCollection);
  virtual void LoadHistoryDataPerSurface(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(LoadHistoryDataPerSurface, COutFieldCollection& fieldCollection);
  virtual void LoadVolumeData(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(LoadVolumeData);
  virtual void LoadSurfaceData(COutFieldCollection& fieldCollection) {};
  CREATE_ACTION(LoadSurfaceData, COutFieldCollection& fieldCollection);


};

#undef CREATE_ACTION

