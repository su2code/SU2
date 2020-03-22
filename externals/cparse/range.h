
#ifndef GENERATORS_H_

#include <iostream>

#include "./shunting-yard.h"
#include "./statement.h"
#include "./block.h"

class Range : public Iterator {
 private:
  struct Startup;
  struct RangeIterator;

 public:
  int64_t from, to, step, i;
  packToken last;

 public:
  Range(int64_t f, int64_t t, int64_t s) : from(f), to(t), step(s), i(f) {}

 public:
  packToken* next();
  void reset();

 public:
  virtual TokenBase* clone() const {
    return new Range(*this);
  }
};

class UserFunction : public Function {
  args_t _args;
  BlockStatement body;
  const std::string _name;

 public:
  const std::string name() const { return _name; }

 public:
  UserFunction(args_t args, BlockStatement body, std::string name = "")
               : _args(args), body(body), _name(name) {}

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

  virtual TokenBase* clone() const {
    return new UserFunction(static_cast<const UserFunction&>(*this));
  }
};

#endif  // GENERATORS_H_
