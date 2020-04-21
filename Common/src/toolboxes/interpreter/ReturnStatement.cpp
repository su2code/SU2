#include <cstring>

#include "../../../include/toolboxes/interpreter/ReturnStatement.hpp"

namespace interpreter {

  constexpr char ReturnStatement::literal[];

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

}
