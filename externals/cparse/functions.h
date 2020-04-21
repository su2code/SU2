#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <list>
#include <string>

typedef std::list<std::string> args_t;

class Function : public TokenBase {
 public:
  static packToken call(packToken _this, const Function* func,
                        TokenList* args, TokenMap scope);
 public:
  Function() : TokenBase(FUNC) {}
  virtual ~Function() {}

 public:
  virtual const std::string name() const = 0;
  virtual const args_t args() const = 0;
  virtual packToken exec(TokenMap scope) const = 0;
  virtual TokenBase* clone() const = 0;
};

class CppFunction : public Function {
 public:
  packToken (*func)(TokenMap, const std::string&);
  args_t _args;
  std::string _name;

  CppFunction(packToken (*func)(TokenMap, const std::string&), const args_t args,
              std::string name = "");
  CppFunction(packToken (*func)(TokenMap, const std::string&), unsigned int nargs,
              const char** args, std::string name = "");
  CppFunction(packToken (*func)(TokenMap, const std::string&), std::string name = "");

  virtual const std::string name() const { return _name; }
  virtual const args_t args() const { return _args; }
  virtual packToken exec(TokenMap scope) const { return func(scope, _name); }

  virtual TokenBase* clone() const {
    return new CppFunction(static_cast<const CppFunction&>(*this));
  }
};

#endif  // FUNCTIONS_H_
