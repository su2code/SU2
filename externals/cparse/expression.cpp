#include "shunting-yard.h"
#include "statement.h"
#include "expression.h"

void ExpStatement::_compile(const char* code, const char** rest,
                            TokenMap parent_scope) {
  expr.compile(code, parent_scope, ";}\n", &code);

  // Skip the delimiter character:
  if (*code && *code != '}') ++code;

  if (rest) *rest = code;
}

returnState ExpStatement::_exec(TokenMap scope) const {
  return expr.eval(scope);
}