#include "exprtk.hpp"
#include "../../../include/toolboxes/parser/ExpressionParser.hpp"
#include "../../../include/mpi_structure.hpp"

namespace Parser {

  typedef exprtk::symbol_table<su2double>        SymbolTable;
  typedef exprtk::expression<su2double>          Expression;
  typedef exprtk::parser<su2double>              Parser;
  typedef exprtk::function_compositor<su2double> Compositor;
  typedef typename Compositor::function          Function;

#include "./CustomFunctions.hpp"

  struct Scope::ScopeImpl{
    SymbolTable symbolTable;

    Parser parser;
    std::list<Expression> expressions;

    std::list<std::unique_ptr<exprtk::igeneric_function<su2double>>> functions;

    ~ScopeImpl() = default;
  };


  Scope::Scope() : Impl(new ScopeImpl()){}

  Scope::~Scope() { delete Impl; }

  void Scope::addVariable(const std::string &name, su2double &var){
    Impl->symbolTable.add_variable(name, var);
  }

  void Scope::addStringVar(const std::string &name, std::string &var){
    Impl->symbolTable.add_stringvar(name, var, false);
  }

  bool Scope::addExpression(const std::string &expr){
    Impl->expressions.emplace_back(Expression(Impl->symbolTable));
    return Impl->parser.compile(expr, Impl->expressions.back());
  }

  std::vector<su2double> Scope::EvalExpressions(){
    std::vector<su2double> vals;
    for (const auto& expr : Impl->expressions){
      vals.push_back(expr.value());
    }
    return vals;
  }

  std::string Scope::GetError(){
    return Impl->parser.error();
  }

  void Scope::addCustomFunction(const std::string& name, CustomFunction type){

    switch (type) {
      case CustomFunction::SURFACE_SUM:
        Impl->functions.emplace_back(new SurfaceSum{&Impl->symbolTable, name});
        Impl->symbolTable.add_function(name, *Impl->functions.back());
        break;
    }
  }

}  // namespace Parser