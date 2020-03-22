#ifndef FUNCTION_H
#define FUNCTION_H

class FuncDeclaration : public Statement {
  args_t args;
  BlockStatement body;
  std::string name;

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
};

/* * * * * Return statements: * * * * */

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
};


#endif // FUNCTION_H
