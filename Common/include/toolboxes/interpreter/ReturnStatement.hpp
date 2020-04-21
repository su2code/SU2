#include "Statement.hpp"

namespace interpreter {

  class ReturnStatement : public Statement {
  protected:
    calculator expr;
    bool value_omitted = true;

  protected:
    void _compile(const char* code, const char** rest, TokenMap parent_scope);
    returnState _exec(TokenMap scope) const;

  public:
    ReturnStatement() {}
    ReturnStatement(const char* code, const char** rest = 0,
                    TokenMap parent_scope = &TokenMap::empty) {
      _compile(code, rest, parent_scope);
    }
    virtual Statement* clone() const {
      return new ReturnStatement(*this);
    }
    static constexpr char literal[] = "return";

  };

}
