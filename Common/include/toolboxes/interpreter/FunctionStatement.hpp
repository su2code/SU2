#include "BlockStatement.hpp"

#pragma once

namespace interpreter {

  enum class FunctionType{
    DEFAULT,
    HISTFIELD,
    VOLUMEFIELD,
    SURFACEINTEGRAL,
    VOLUMEINTEGRAL
  };

  extern const std::map<std::string, FunctionType> FuncTypeMap;

  class FuncDeclaration : public Statement {
    args_t args;
    BlockStatement body;
    std::string name;
    std::string type;

  private:
    void _compile(const char* code, const char** rest, TokenMap parent_scope);
    returnState _exec(TokenMap scope) const;

    void _compile(std::string name, const char* code, const char** rest = 0,
                  TokenMap parent_scope = &TokenMap::empty);
  public:
    FuncDeclaration() {}
    FuncDeclaration(const char* code, const char** rest = 0,
                    TokenMap parent_scope = &TokenMap::empty) {
      _compile(code, rest, parent_scope);
    }
    FuncDeclaration(const std::string& name,
                    const char* code, const char** rest = 0,
                    TokenMap parent_scope = &TokenMap::empty) {
      _compile(name, code, rest, parent_scope);
    }


  public:
    packToken asFunc() const;

    virtual Statement* clone() const {
      return new FuncDeclaration(*this);
    }

    static constexpr char literal[] = "def";
  };

  class UserFunction : public Function {
    args_t _args;
    BlockStatement body;
    const std::string _name;
    FunctionType _type;

  public:
    const std::string name() const { return _name; }

  public:
    UserFunction(args_t args, BlockStatement body, std::string name, FunctionType type)
      : _args(args), body(body), _name(name), _type(type) {}

    virtual const args_t args() const { return _args; }
    virtual packToken exec(TokenMap scope) const {
      returnState st;
      st = body.exec(scope);
      if (st.type == RETURN || st.type == FINISH) {
        return st.value;
      } else {
        return packToken::None();
      }
    }

    FunctionType getType() const {
      return _type;
    }

    virtual TokenBase* clone() const {
      return new UserFunction(static_cast<const UserFunction&>(*this));
    }
  };

  std::vector<interpreter::UserFunction*> GetUserFunctions(const TokenMap& scope, std::list<interpreter::FunctionType> functionType);

}
