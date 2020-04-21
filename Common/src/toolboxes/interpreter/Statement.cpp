#include<cstring>

#include "../../../include/toolboxes/interpreter/BlockStatement.hpp"
#include "../../../include/toolboxes/interpreter/ExpressionStatement.hpp"
#include "../../../include/toolboxes/interpreter/FunctionStatement.hpp"
#include "../../../include/toolboxes/interpreter/ReturnStatement.hpp"

namespace interpreter {

  bool checkStatement(const char* code, const char *literal, uint size){
    uint i=0;
    for (i = 0; i < size; i++){
      if (code[i] != literal[i]) return false;
    }

    if ((i == size) && !(isalnum(code[i]) || code[i] == '_'))
      return true;

    return false;
  };


  // Decide what type of statement to build:
  Statement* Statement::buildStatement(const char** source, TokenMap scope) {
    const char* code = *source;

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