#pragma once

#include "shunting-yard.h"
#include "shunting-yard-exceptions.h"
#include "statement.h"
#include "block.h"
#include "../datatype_structure.hpp"

class CExpressionParser {

private:
  BlockStatement code;
  const TokenMap* scope;
  const char* rest;
  std::string name;
  const packToken* tokenRef;

public:
  CExpressionParser() : scope(nullptr) {}

  CExpressionParser(const TokenMap* scope_) : scope(scope_){}

  bool CompileAndExec(const std::string& codeAsString, const std::string& funcName){
    name = funcName;
    code.compile(codeAsString.c_str(), &rest, *scope);
    code.exec(*scope);
    tokenRef = scope->find(name);
    if (tokenRef) return true;
    else return false;
  }

  void Compile(const std::string& codeAsString){
    code.compile(codeAsString.c_str(), &rest, *scope);
  }

  void ExecCode(){
    code.exec(*scope);
  }

  su2double Eval(){
    return tokenRef->asFunc()->exec(*scope).asDouble();
  }
};
