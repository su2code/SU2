#include <cstdlib>
#include <cstring> // For strchr
#include <sstream>
#include <iostream>
#include <algorithm>

#include "shunting-yard-exceptions.h"
#include "../../../include/toolboxes/interpreter/FunctionStatement.hpp"

namespace interpreter {

  constexpr char FuncDeclaration::literal[];

  const std::map<std::string, FunctionType> FuncTypeMap = {
    {"",            FunctionType::DEFAULT},
    {"historyfield", FunctionType::HISTFIELD},
    {"volumefield", FunctionType::VOLUMEFIELD},
    {"surfaceintegral", FunctionType::SURFACEINTEGRAL}
  };

  /* * * * * Utility functions: * * * * */

  std::string parseName(const char** source, char end_char = '\0') {
    std::stringstream ss;
    const char* code = *source;

    // Parse the function name:
    if (isalpha(*code) || *code == '_') {
      ss << *code;
      ++code;

      while (isalnum(*code) || *code == '_') {
        ss << *code;
        ++code;
      }

      // Find the beggining of the non space character:
      while (isspace(*code) && *code != end_char) ++code;
    } else {
      throw syntax_error("Expected variable name!");
    }

    *source = code;

    return ss.str();
  }
  /* * * * * FuncDeclaration Statement * * * * */

  void FuncDeclaration::_compile(const char* code, const char** rest,
                                 TokenMap parent_scope) {

    // Find the start of the name:
    while (isspace(*code)) ++code;

    // Parse the function name:
    try {
      name = parseName(&code);
    } catch(syntax_error e) {
      throw syntax_error("Missing name after `function` key-word!");
    }

    try {
      type = parseName(&code);
    }  catch (syntax_error) {
      /*--- Do nothing ---*/
    }

    if (!type.empty()){
      std::string tmp = name;
      name = type;
      type = tmp;
    }

    _compile(name, code, rest, parent_scope);

    // Ignore white spaces:
    while (isspace(*code)) ++code;
  }

  void FuncDeclaration::_compile(std::string name, const char* code,
                                 const char** rest, TokenMap parent_scope) {
    // Make sure its empty:
    args.clear();

    if (*code != '(') {
      throw syntax_error("Expected argument list after `function` reserved word!");
    }

    // Find the next non-blank character:
    ++code;
    while (isspace(*code)) ++code;

    // Parse each argument of the block:
    while (*code && *code != ')') {

      if (*code == ',') {
        throw syntax_error("Empty item on argument list!");
      } else {

        try {
          // Parse the argument name:
          args.push_back(parseName(&code));
        } catch (syntax_error e) {
          throw syntax_error("Invalid argument name!");
        }

        if(*code == ',') {
          ++code;
        } else if (*code != ')') {
          if (isalpha(*code) || *code == '_') {
            syntax_error("Missing ',' in argument list!");
          } else {
            syntax_error("Invalid character in argument list!");
          }
        }
      }

      // Discard blank spaces:
      while (isspace(*code)) ++code;
    }

    if (*code != ')') {
      throw syntax_error("Expected function body before end of file!");
    } else {
      ++code;
    }

    body.compile(code, &code, parent_scope);

    if (rest) *rest = code;
  }

  returnState FuncDeclaration::_exec(TokenMap scope) const {
    scope[name] = this->asFunc();

    return NORMAL;
  }

  packToken FuncDeclaration::asFunc() const {
    FunctionType funcType;
    try {
      funcType = FuncTypeMap.at(type);
    } catch (std::out_of_range& e){
      throw syntax_error("Undefined user function type \"" + type + "\".");
    }
    if ((funcType != FunctionType::DEFAULT) && (args.size() != 0)){
      throw syntax_error("No arguments allowed for specific user function type " + type);
    }
    return packToken(new UserFunction(args, body, name, funcType));
  }

  std::vector<interpreter::UserFunction*> GetUserFunctions(const TokenMap& scope,
                                                           std::list<interpreter::FunctionType> functionType) {
    std::vector<interpreter::UserFunction*> userFunctions;
    for (auto iter = scope.map().begin();
         iter != scope.map().end(); iter++){
      interpreter::UserFunction* function = dynamic_cast<interpreter::UserFunction*>(iter->second.token());
      if (function != NULL){
        if (std::any_of(functionType.begin(), functionType.end(),
                        [&](const interpreter::FunctionType& item){return item == function->getType();})){
          userFunctions.push_back(function);
        }
      }
    }
    return userFunctions;
  }

}