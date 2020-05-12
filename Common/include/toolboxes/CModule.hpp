#include <tuple>
#include <utility>
#pragma once


#define CREATE_ACTION(name, func_args...) \
  static struct Action##name {} Action##name; \
  template <typename type, typename... Targs> \
  static void Action(type& Module, struct Action##name, Targs& ...args){ Module.name(args...); }


template<class DerivedModule>
class CModule {

protected:

  bool enabled;

public:

  CModule(bool enabled_ = true) :
    enabled(enabled_){}

  bool IsEnabled() { return enabled; }

  template<std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
  inline static typename std::enable_if<I == sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...> &, TAction, Targs&...) { }

  template<std::size_t I = 0, typename TAction, typename... Tp, typename... Targs>
  inline static typename std::enable_if<I < sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...>& t, TAction action, Targs&... args)
  {
    if (std::get<I>(t).IsEnabled()){
      DerivedModule::Action(std::get<I>(t), action, args...);
    }
    for_each<I + 1, TAction, Tp...>(t, action, args...);
  }

};


template<typename T, typename... Modules>
class CModuleList : public std::tuple<Modules...>{
public:
  explicit CModuleList(T* t): std::tuple<Modules...>(Modules(t)...){};
};

