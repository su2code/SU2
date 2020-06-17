#include <tuple>
#include <utility>

#pragma once

#include "../fields/COutFieldManager.hpp"
#include "SolverDataContainer.hpp"

#define CREATE_ACTION(name, func_args...)                                 \
  static struct Action##name {                                            \
  } Action##name;                                                         \
  template <typename type, typename... Targs>                             \
  static void Action(type& Module, struct Action##name, Targs&... args) { \
    Module.name(args...);                                                 \
  }

template <class Module>
class COutputModule {
 protected:
  bool enabled;

 public:
  COutputModule(bool enabled_ = true) : enabled(enabled_) {}

  bool IsEnabled() { return enabled; }

  template <std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
  inline static typename std::enable_if<I == sizeof...(Tp), void>::type for_each(std::tuple<Tp...>&, TAction,
                                                                                 Targs&&...) {}

  template <std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
      inline static typename std::enable_if <
      I<sizeof...(Tp), void>::type for_each(std::tuple<Tp...>& t, TAction action, Targs&&... args) {
    if (std::get<I>(t).IsEnabled()) {
      Module::Action(std::get<I>(t), action, args...);
    }
    for_each<I + 1, TAction, Tp...>(t, action, args...);
  }
};

class CSolverOutputModule : public COutputModule<CSolverOutputModule> {
 protected:
  int nDim;

 public:
  explicit CSolverOutputModule(int nDim, bool enabled_ = true) : COutputModule(enabled_), nDim(nDim){};

  virtual void DefineHistoryFields(CHistoryOutFieldManager& historyFields) {}
  CREATE_ACTION(DefineHistoryFields, CHistoryOutFieldManager& historyFields);
  virtual void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) {}
  CREATE_ACTION(DefineHistoryFieldModifier, CHistoryOutFieldManager& historyFields);
  virtual void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {}
  CREATE_ACTION(DefineVolumeFields, CVolumeOutFieldManager& volumeFields);
  virtual void LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                               const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadHistoryData, CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                const IterationInfo& iterationInfo);
  virtual void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                       const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadHistoryDataModifier, CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                const IterationInfo& iterationInfo);
  virtual void LoadHistoryDataPerSurface(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                         const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadHistoryDataPerSurface, HistoryOutFieldManager& historyFields, const SolverData& solverData,
                const IterationInfo& iterationInfo);
  virtual void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                              const IterationInfo& iterationInfo, const PointInfo& pointInfo) {}
  CREATE_ACTION(LoadVolumeData, CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                const IterationInfo& iterationInfo, const PointInfo& pointInfo);
  virtual void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                               const IterationInfo& iterationInfo, const PointInfo& pointInfo) {}
  CREATE_ACTION(LoadSurfaceData, CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                const IterationInfo& iterationInfo, const PointInfo& pointInfo);

  virtual void PrintToStream(std::ostream* stream, const CHistoryOutFieldManager& historyFields) {}
  CREATE_ACTION(PrintToStream, std::ostream* stream, const CHistoryOutFieldManager& historyFields);
};

class CModifierModule : public COutputModule<CModifierModule> {
 protected:
  int nDim;

 public:
  explicit CModifierModule(int nDim, bool enabled_ = true) : COutputModule(enabled_), nDim(nDim){};

  virtual void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) {}
  CREATE_ACTION(DefineHistoryFieldModifier, CHistoryOutFieldManager& historyFields);
  virtual void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields, const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadHistoryDataModifier, CHistoryOutFieldManager& historyFields, const IterationInfo& iterationInfo);
  virtual void DefineVolumeFieldModifier(CVolumeOutFieldManager& volumeFields) {}
  CREATE_ACTION(DefineVolumeFieldModifier, CVolumeOutFieldManager& volumeFields);
  virtual void LoadVolumeDataModifier(CVolumeOutFieldManager& volumeFields, const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadVolumeDataModifier, CVolumeOutFieldManager& volumeFields, const IterationInfo& iterationInfo);
  virtual void LoadSurfaceDataModifier(CVolumeOutFieldManager& volumeFields, const IterationInfo& iterationInfo) {}
  CREATE_ACTION(LoadSurfaceDataModifier, CVolumeOutFieldManager& volumeFields, const IterationInfo& iterationInfo);
  virtual void PrintToStream(std::ostream* stream, const CHistoryOutFieldManager& historyFields) {}
  CREATE_ACTION(PrintToStream, std::ostream* stream, const CHistoryOutFieldManager& historyFields);
};

#undef CREATE_ACTION
