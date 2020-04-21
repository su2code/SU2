#include "Statement.hpp"
#pragma once
namespace interpreter {

  class ExpStatement : public Statement {
    calculator expr;

  private:
    void _compile(const char* code, const char** rest, TokenMap parent_scope);
    returnState _exec(TokenMap scope) const;

  public:
    ExpStatement() {}
    ExpStatement(const char* code, const char** rest = 0,
                 TokenMap parent_scope = &TokenMap::empty) {
      _compile(code, rest, parent_scope);
    }
    virtual Statement* clone() const {
      return new ExpStatement(*this);
    }
  };

}
