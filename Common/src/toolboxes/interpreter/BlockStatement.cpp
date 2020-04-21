#include <cstring> // For strchr

#include "shunting-yard-exceptions.h"
#include "../../../include/toolboxes/interpreter/BlockStatement.hpp"

namespace interpreter {

  constexpr char BlockStatement::literal[];

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
          list.push_back(Statement::buildStatement(&code, parent_scope));
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
      list.push_back(Statement::buildStatement(&code, parent_scope));
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

}
