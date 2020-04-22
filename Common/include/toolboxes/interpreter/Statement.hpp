#include "shunting-yard.h"
#include "shunting-yard-exceptions.h"
#pragma once

namespace interpreter {

  extern GlobalScope globalScope;

  enum returnType {
    NORMAL, FINISH, RETURN,
    YIELD, BREAK, CONTINUE, THROW
  };

  struct returnState {
    uint8_t type;
    packToken value;
    bool value_omitted = true;
    returnState() : type(NORMAL), value(packToken::None()) {}
    returnState(const returnType& type) : type(type), value(packToken::None()) {}
    returnState(packToken value) : type(NORMAL), value(value) {}
    returnState(const returnType& type, packToken value, bool value_omitted = true)
      : type(type), value(value), value_omitted(value_omitted) {}
  };

  class Statement {
  protected:
    virtual void _compile(const char* code, const char** rest,
                          TokenMap parent_scope) = 0;
    virtual returnState _exec(TokenMap scope) const = 0;

  public:
    virtual ~Statement() {}
    void compile(const char* code, const char** rest = 0,
                 TokenMap parent_scope = &TokenMap::empty) {
      return _compile(code, rest, parent_scope);
    }

    returnState exec(TokenMap scope) const { return _exec(scope); }

    virtual Statement* clone() const = 0;

    static Statement* buildStatement(const char** source, TokenMap scope);

  };

}
