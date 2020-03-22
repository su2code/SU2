#pragma once

#include "shunting-yard.h"
#include "statement.h"
#include "block.h"
#include "../datatype_structure.hpp"

class CExpressionParser {

private:
  BlockStatement code;
  const TokenMap* scope;
  const char* rest;

public:
  CExpressionParser() : scope(nullptr) {}

  CExpressionParser(const TokenMap* scope_) : scope(scope_){}

  void Compile(const std::string& codeAsString){
    code.compile(codeAsString.c_str(), &rest, *scope);
  }


  void ExecCode(){
    code.exec(*scope);
  }

  su2double Eval( const std::string& name){
    return (*scope).find(name)->asFunc()->exec(*scope).asDouble();
  }
};
