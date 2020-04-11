#ifndef BLOCK_H
#define BLOCK_H

class BlockStatement : public Statement {
  typedef std::list<Statement*> codeBlock_t;
  codeBlock_t list;

 private:
  void cleanList(codeBlock_t* list);

 private:
  void _compile(const char* code, const char** rest, TokenMap parent_scope);
  returnState _exec(TokenMap scope) const;

 public:
  BlockStatement() {}
  BlockStatement(const char* code, const char** rest = 0,
                 TokenMap parent_scope = &TokenMap::empty) {
    _compile(code, rest, parent_scope);
  }
  BlockStatement(const char* code,
                 TokenMap parent_scope = &TokenMap::empty) {
    _compile(code, 0, parent_scope);
  }
  BlockStatement(const BlockStatement& other);
  ~BlockStatement();
  BlockStatement& operator=(const BlockStatement& other);
  uint32_t size() { return list.size(); }

  virtual Statement* clone() const {
    return new BlockStatement(*this);
  }
};


class IfStatement : public Statement {
  calculator cond;
  BlockStatement _then;
  BlockStatement _else;

 private:
  void _compile(const char* code, const char** rest, TokenMap parent_scope);
  returnState _exec(TokenMap scope) const;

 public:
  IfStatement() {}
  IfStatement(const char* code, const char** rest = 0,
              TokenMap parent_scope = &TokenMap::empty) {
    _compile(code, rest, parent_scope);
  }
  virtual Statement* clone() const {
    return new IfStatement(*this);
  }
};
#endif // BLOCK_H
