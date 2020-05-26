#pragma once
#include "../../datatype_structure.hpp"
#include <fstream>
#include <list>
#include "../printing_toolbox.hpp"
namespace Parser{


  enum class CustomFunction {
    SURFACE_SUM
  };


  class Scope {

  private:
    struct ScopeImpl;
    ScopeImpl* Impl;
  public:

    Scope();

    ~Scope();

    // Asssignment Operator and Copy Constructor

     Scope(const Scope& other) = delete;
     Scope& operator=(Scope rhs) = delete;

    void addVariable(const std::string& name, su2double &var);

    void addStringVar(const std::string& name, std::string &var);

    bool addExpression(const std::string &expr);

    std::string GetError();

    std::vector<su2double> EvalExpressions();

    void addCustomFunction(const std::string& name, CustomFunction type);

  };
};
