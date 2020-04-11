#include "shunting-yard.h"
#include "statement.h"
#include "block.h"
#include "expression.h"
#include "funcdecl.h"
#include <cstring>
#include "shunting-yard-exceptions.h"

// Decide what type of statement to build:
Statement* buildStatement(const char** source, TokenMap scope) {
  const char* code = *source;
  const char* _template;
  uint i;

  switch (*code) {
    case '{':
      return new BlockStatement(code, source, scope);
    case 'i':
      if (code[1] == 'f' && !(isalnum(code[2]) || code[2] == '_'))
        return new IfStatement(code+2, source, scope);
      break;
    case 'r':
      _template = "return";
      for (i = 1; i < 6; ++i)
        if (code[i] != _template[i]) break;

      if (i == 6 && !(isalnum(code[i]) || code[i] == '_'))
        return new ReturnStatement(code+6, source, scope);
      break;
    case 'f':

      _template = "function";
      for (i = 1; i < 8; ++i)
        if (code[i] != _template[i]) break;

      if (i == 8 && !(isalnum(code[i]) || code[i] == '_'))
        return new FuncDeclaration(code+8, source, scope);

  }

  return new ExpStatement(code, source, scope);
}

void BlockStatement::cleanList(codeBlock_t* list) {
  for(auto stmt : *list) {
    delete stmt;
  }

  list->clear();
}

BlockStatement::BlockStatement(const BlockStatement& other) {
  for(const Statement* stmt : other.list) {
    list.push_back(stmt->clone());
  }
}

BlockStatement& BlockStatement::operator=(const BlockStatement& other) {
  cleanList(&list);
  for(const Statement* stmt : other.list) {
    list.push_back(stmt->clone());
  }
  return *this;
}

BlockStatement::~BlockStatement() {
  cleanList(&list);
}

void BlockStatement::_compile(const char* code, const char** rest,
                              TokenMap parent_scope) {
  // Make sure the list is empty:
  cleanList(&list);

  while (isspace(*code)) ++code;

  if (*code == '{') {

    // Find the next non-blank character:
    ++code;
    while (isspace(*code)) ++code;

    // Parse each statement of the block:
    while (*code && *code != '}') {
      // Ignore empty statements:
      if (strchr(";\n", *code)) {
        ++code;
      } else {
        list.push_back(buildStatement(&code, parent_scope));
      }

      // Discard blank spaces:
      while (isspace(*code)) ++code;
    }

    if (*code == '}') {
      ++code;
    } else {
      throw syntax_error("Missing a '}' somewhere on the code!");
    }
  } else {
    list.push_back(buildStatement(&code, parent_scope));
  }

  if (rest) *rest = code;
}

returnState BlockStatement::_exec(TokenMap scope) const {
  returnState rs;
  for(const auto stmt : list) {
    rs = stmt->exec(scope);
    if (rs.type != NORMAL) return rs;
  }

  return rs;
}

void IfStatement::_compile(const char* code, const char** rest,
                           TokenMap parent_scope) {

  while (isspace(*code)) ++code;

  if (*code != '(') {
    throw syntax_error("Expected '(' after `if` statement!");
  }

  // Parse the condition:
  cond.compile(code+1, parent_scope, ")", &code);

  if (*code != ')') {
    throw syntax_error("Missing ')' after `if` statement!");
  }

  _then.compile(code+1, &code, parent_scope);

  while (isspace(*code)) ++code;

  // Check if the next word is else:
  static const char* str = "else";
  for (int i = 0; i < 4; ++i) {
    if (str[i] != code[i]) {
      if(rest) *rest = code;
      return;
    }
  }
  if (isalnum(code[4]) || code[4] == '_') {
    if(rest) *rest = code;
    return;
  }

  _else.compile(code+4, &code, parent_scope);

  if(rest) *rest = code;
}

returnState IfStatement::_exec(TokenMap scope) const {
  if (cond.eval(scope).asBool()) {
    return _then.exec(scope);
  } else {
    return _else.exec(scope);
  }
}