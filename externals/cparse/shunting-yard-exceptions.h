
#ifndef SHUNTING_YARD_EXCEPTIONS_H_
#define SHUNTING_YARD_EXCEPTIONS_H_

#include "./shunting-yard.h"

#include <string>
#include <stdexcept>

class msg_exception : public std::exception {
 protected:
  const std::string msg;
 public:
  msg_exception(const std::string& msg) : msg(msg) {}
  ~msg_exception() throw() {}
  const char* what() const throw() {
    return msg.c_str();
  }
};

struct bad_cast : public msg_exception {
  bad_cast(const std::string& msg) : msg_exception(msg) {}
};

struct syntax_error : public msg_exception {
  syntax_error(const std::string& msg) : msg_exception(msg) {}
};

struct type_error : public msg_exception {
  type_error(const std::string& msg) : msg_exception(msg) {}
};

struct undefined_operation : public msg_exception {
  undefined_operation(const std::string& op, const TokenBase* left, const TokenBase* right)
                      : undefined_operation(op, packToken(left->clone()), packToken(right->clone())) {}
  undefined_operation(const std::string& op, const packToken& left, const packToken& right)
    : msg_exception("Unexpected operation with operator '" + op + "' and operands: " + left.str() + " and " + right.str() + ".") {}
};

#endif  // SHUNTING_YARD_EXCEPTIONS_H_
