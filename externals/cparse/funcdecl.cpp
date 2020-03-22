#include <cstdlib>
#include <cstring> // For strchr
#include <sstream>
#include <iostream>

#include "shunting-yard.h"
#include "statement.h"
#include "block.h"
#include "funcdecl.h"
#include "shunting-yard-exceptions.h"
#include "range.h"

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
  return packToken(new UserFunction(args, body, name));
}

void ReturnStatement::_compile(const char* code, const char** rest,
                               TokenMap parent_scope) {
  while (isspace(*code)) ++code;

  if (strchr(";}\n", *code)) {
    expr.compile("None");
  } else {
    expr.compile(code, parent_scope, ";}\n", &code);
    value_omitted = false;
  }

  if (*code && *code != '}') ++code;

  if (rest) *rest = code;
}

returnState ReturnStatement::_exec(TokenMap scope) const {
  return returnState(RETURN, expr.eval(scope), value_omitted);
}