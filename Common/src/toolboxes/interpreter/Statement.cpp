#include<cstring>

#include "../../../include/toolboxes/interpreter/BlockStatement.hpp"
#include "../../../include/toolboxes/interpreter/ExpressionStatement.hpp"
#include "../../../include/toolboxes/interpreter/FunctionStatement.hpp"
#include "../../../include/toolboxes/interpreter/ReturnStatement.hpp"
#define STATIC_CPARSE_STARTUP
#include "./builtin-features.inc"

namespace interpreter {

  GlobalScope globalScope;

  bool checkStatement(const char* code, const char *literal, uint size){
    uint i=0;
    /*--- Check that the full name of the token matches the statement literal ---*/

    for (i = 0; i < size; i++){
      if (code[i] != literal[i]) return false;
    }

    /*--- Check that the name does not continue ---*/

    if ((i == size) && !(isalnum(code[i]) || code[i] == '_'))
      return true;

    return false;
  };


  Statement* Statement::buildStatement(const char** source, TokenMap scope) {
    const char* code = *source;

    // Decide what type of statement to build by looking at the first character

    switch (*code) {
      case BlockStatement::literal[0]:
        return new BlockStatement(code, source, scope);
        break;
      case ReturnStatement::literal[0]:
        {
          const uint size = strlen(ReturnStatement::literal);
          if (checkStatement(code, ReturnStatement::literal, size)){
            return new ReturnStatement(code+size, source, scope);
          }
        }
        break;
      case FuncDeclaration::literal[0]:
        {
          const uint size = strlen( FuncDeclaration::literal);
          if (checkStatement(code,  FuncDeclaration::literal, size)){
            return new FuncDeclaration(code+size, source, scope);
          }
        }
    }

    return new ExpStatement(code, source, scope);
  }

}