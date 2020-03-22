#include "./shunting-yard.h"
#include "./shunting-yard-exceptions.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <exception>
#include <string>
#include <stack>
#include <utility>  // For std::pair
#include <cstring>  // For strchr()

/* * * * * Operation class: * * * * */

// Convert a type into an unique mask for bit wise operations:
uint32_t Operation::mask(tokType_t type) {
  if (type == ANY_TYPE) {
    return 0xFFFF;
  } else {
    return ((type & 0xE0) << 24) | (1 << (type & 0x1F));
  }
}

// Build a mask for each pair of operands
opID_t Operation::build_mask(tokType_t left, tokType_t right) {
  opID_t result = mask(left);
  return (result << 32) | mask(right);
}

/* * * * * Utility functions: * * * * */

bool match_op_id(opID_t id, opID_t mask) {
  uint64_t result = id & mask;
  uint32_t* val = reinterpret_cast<uint32_t*>(&result);
  if (val[0] && val[1]) return true;
  return false;
}

TokenBase* exec_operation(const packToken& left, const packToken& right,
                          evaluationData* data, const std::string& OP_MASK) {
  auto it = data->opMap.find(OP_MASK);
  if (it == data->opMap.end()) return 0;
  for (const Operation& operation : it->second) {
    if (match_op_id(data->opID, operation.getMask())) {
      try {
        return operation.exec(left, right, data).release();
      } catch (const Operation::Reject& e) {
        continue;
      }
    }
  }

  return 0;
}

inline std::string normalize_op(std::string op) {
  if (op[0] == 'L' || op[0] == 'R') {
    op.erase(0, 1);
    return op;
  } else {
    return op;
  }
}

// Use this function to discard a reference to an object
// And obtain the original TokenBase*.
// Please note that it only deletes memory if the token
// is of type REF.
TokenBase* resolve_reference(TokenBase* b, TokenMap* scope = 0) {
  if (b->type & REF) {
    // Resolve the reference:
    RefToken* ref = static_cast<RefToken*>(b);
    TokenBase* value = ref->resolve(scope);

    delete ref;
    return value;
  }

  return b;
}

/* * * * * Static containers: * * * * */

// Build configurations once only:
Config_t& calculator::Default() {
  static Config_t conf;
  return conf;
}

typeMap_t& calculator::type_attribute_map() {
  static typeMap_t type_map;
  return type_map;
}

/* * * * * rpnBuilder Class: * * * * */

void rpnBuilder::cleanRPN(TokenQueue_t* rpn) {
  while (rpn->size()) {
    delete resolve_reference(rpn->front());
    rpn->pop();
  }
}

/**
 * Consume operators with precedence >= than op
 * and add them to the RPN
 *
 * The algorithm works as follows:
 *
 * Let p(o) denote the precedence of an operator o.
 *
 * If the token is an operator, o1, then
 *   While there is an operator token, o2, at the top
 *       and p(o1) >= p(o2) (`>` for Right to Left associativity)
 *     then:
 *
 *     pop o2 off the stack onto the output queue.
 *   Push o1 on the stack.
 */
void rpnBuilder::handle_opStack(const std::string& op) {
  std::string cur_op;

  // If it associates from left to right:
  if (opp.assoc(op) == 0) {
    while (!opStack.empty() &&
        opp.prec(op) >= opp.prec(opStack.top())) {
      cur_op = normalize_op(opStack.top());
      rpn.push(new Token<std::string>(cur_op, OP));
      opStack.pop();
    }
  } else {
    while (!opStack.empty() &&
        opp.prec(op) > opp.prec(opStack.top())) {
      cur_op = normalize_op(opStack.top());
      rpn.push(new Token<std::string>(cur_op, OP));
      opStack.pop();
    }
  }
}

void rpnBuilder::handle_binary(const std::string& op) {
  // Handle OP precedence
  handle_opStack(op);
  // Then push the current op into the stack:
  opStack.push(op);
}

// Convert left unary operators to binary and handle them:
void rpnBuilder::handle_left_unary(const std::string& unary_op) {
  this->rpn.push(new TokenUnary());
  // Only put it on the stack and wait to check op precedence:
  opStack.push(unary_op);
}

// Convert right unary operators to binary and handle them:
void rpnBuilder::handle_right_unary(const std::string& unary_op) {
  // Handle OP precedence:
  handle_opStack(unary_op);
  // Add the unary token:
  this->rpn.push(new TokenUnary());
  // Then add the current op directly into the rpn:
  rpn.push(new Token<std::string>(normalize_op(unary_op), OP));
}

// Find out if op is a binary or unary operator and handle it:
void rpnBuilder::handle_op(const std::string& op) {
  // If its a left unary operator:
  if (this->lastTokenWasOp) {
    if (opp.exists("L"+op)) {
      handle_left_unary("L"+op);
      this->lastTokenWasUnary = true;
      this->lastTokenWasOp = op[0];
    } else {
      cleanRPN(&(this->rpn));
      throw std::domain_error(
          "Unrecognized unary operator: '" + op + "'.");
    }

  // If its a right unary operator:
  } else if (opp.exists("R"+op)) {
    handle_right_unary("R"+op);

    // Set it to false, since we have already added
    // an unary token and operand to the stack:
    this->lastTokenWasUnary = false;
    this->lastTokenWasOp = false;

  // If it is a binary operator:
  } else {
    if (opp.exists(op)) {
      handle_binary(op);
    } else {
      cleanRPN(&(rpn));
      throw std::domain_error(
          "Undefined operator: `" + op + "`!");
    }

    this->lastTokenWasUnary = false;
    this->lastTokenWasOp = op[0];
  }
}

void rpnBuilder::handle_token(TokenBase* token) {
  rpn.push(token);
  lastTokenWasOp = false;
  lastTokenWasUnary = false;
}

void rpnBuilder::open_bracket(const std::string& bracket) {
  opStack.push(bracket);
  lastTokenWasOp = bracket[0];
  lastTokenWasUnary = false;
  ++bracketLevel;
}

void rpnBuilder::close_bracket(const std::string& bracket) {
  if (lastTokenWasOp == bracket[0]) {
    rpn.push(new Tuple());
  }

  std::string cur_op;
  while (opStack.size() && opStack.top() != bracket) {
    cur_op = normalize_op(opStack.top());
    rpn.push(new Token<std::string>(cur_op, OP));
    opStack.pop();
  }

  if (opStack.size() == 0) {
    rpnBuilder::cleanRPN(&rpn);
    throw syntax_error("Extra '" + bracket + "' on the expression!");
  }

  opStack.pop();
  lastTokenWasOp = false;
  lastTokenWasUnary = false;
  --bracketLevel;
}

/* * * * * RAII_TokenQueue_t struct  * * * * */

// Used to make sure an rpn is dealloc'd correctly
// even when an exception is thrown.
//
// Note: This is needed because C++ does not
// allow a try-finally block.
struct calculator::RAII_TokenQueue_t : TokenQueue_t {
  RAII_TokenQueue_t() {}
  RAII_TokenQueue_t(const TokenQueue_t& rpn) : TokenQueue_t(rpn) {}
  ~RAII_TokenQueue_t() { rpnBuilder::cleanRPN(this); }

  RAII_TokenQueue_t(const RAII_TokenQueue_t& rpn) {
    throw std::runtime_error("You should not copy this class!");
  }
  RAII_TokenQueue_t& operator=(const RAII_TokenQueue_t& rpn) {
    throw std::runtime_error("You should not copy this class!");
  }
};

/* * * * * calculator class * * * * */

TokenQueue_t calculator::toRPN(const char* expr,
                               TokenMap vars, const char* delim,
                               const char** rest, Config_t config) {
  rpnBuilder data(vars, config.opPrecedence);
  char* nextChar;

  static char c = '\0';
  if (!delim) delim = &c;

  while (*expr && isspace(*expr) && !strchr(delim, *expr)) ++expr;

  if (*expr == '\0' || strchr(delim, *expr)) {
    throw std::invalid_argument("Cannot build a calculator from an empty expression!");
  }

  // In one pass, ignore whitespace and parse the expression into RPN
  // using Dijkstra's Shunting-yard algorithm.
  while (*expr && (data.bracketLevel || !strchr(delim, *expr))) {
    if (isdigit(*expr)) {
      // If the token is a number, add it to the output queue.
      int64_t _int = strtoll(expr, &nextChar, 10);

      // If the number was not a float:
      if (!strchr(".eE", *nextChar)) {
        data.handle_token(new Token<int64_t>(_int, INT));
      } else {
        doubleType digit = strtod(expr, &nextChar);
        data.handle_token(new Token<doubleType>(digit, REAL));
      }

      expr = nextChar;
    } else if (rpnBuilder::isvarchar(*expr)) {
      rWordParser_t* parser;

      // If the token is a variable, resolve it and
      // add the parsed number to the output queue.
      std::string key = rpnBuilder::parseVar(expr, &expr);

      if ((parser=config.parserMap.find(key))) {
        // Parse reserved words:
        try {
          parser(expr, &expr, &data);
        } catch (...) {
          rpnBuilder::cleanRPN(&data.rpn);
          throw;
        }
      } else {
        packToken* value = vars.find(key);

        if (value) {
          // Save a reference token:
          TokenBase* copy = (*value)->clone();
          data.handle_token(new RefToken(key, copy));
        } else {
          // Save the variable name:
          data.handle_token(new Token<std::string>(key, VAR));
        }
      }
    } else if (*expr == '\'' || *expr == '"') {
      // If it is a string literal, parse it and
      // add to the output queue.
      char quote = *expr;

      ++expr;
      std::stringstream ss;
      while (*expr && *expr != quote && *expr != '\n') {
        if (*expr == '\\') {
          switch (expr[1]) {
          case 'n':
            expr+=2;
            ss << '\n';
            break;
          case 't':
            expr+=2;
            ss << '\t';
            break;
          default:
            if (strchr("\"'\n", expr[1])) ++expr;
            ss << *expr;
            ++expr;
          }
        } else {
          ss << *expr;
          ++expr;
        }
      }

      if (*expr != quote) {
        std::string squote = (quote == '"' ? "\"": "'");
        rpnBuilder::cleanRPN(&data.rpn);
        throw syntax_error("Expected quote (" + squote +
                           ") at end of string declaration: " + squote + ss.str() + ".");
      }
      ++expr;
      data.handle_token(new Token<std::string>(ss.str(), STR));
    } else {
      // Otherwise, the variable is an operator or paranthesis.
      switch (*expr) {
      case '(':
        // If it is a function call:
        if (data.lastTokenWasOp == false) {
          // This counts as a bracket and as an operator:
          data.handle_op("()");
          // Add it as a bracket to the op stack:
        }
        data.open_bracket("(");
        ++expr;
        break;
      case '[':
        if (data.lastTokenWasOp == false) {
          // If it is an operator:
          data.handle_op("[]");
        } else {
          // If it is the list constructor:
          // Add the list constructor to the rpn:
          data.handle_token(new CppFunction(&TokenList::default_constructor, "list"));

          // We make the program see it as a normal function call:
          data.handle_op("()");
        }
        // Add it as a bracket to the op stack:
        data.open_bracket("[");
        ++expr;
        break;
      case '{':
        // Add a map constructor call to the rpn:
        data.handle_token(new CppFunction(&TokenMap::default_constructor, "map"));

        // We make the program see it as a normal function call:
        data.handle_op("()");
        data.open_bracket("{");
        ++expr;
        break;
      case ')':
        data.close_bracket("(");
        ++expr;
        break;
      case ']':
        data.close_bracket("[");
        ++expr;
        break;
      case '}':
        data.close_bracket("{");
        ++expr;
        break;
      default:
        {
          // Then the token is an operator

          const char* start = expr;
          std::stringstream ss;
          ss << *expr;
          ++expr;
          while (*expr && ispunct(*expr) && !strchr("+-'\"()[]{}_", *expr)) {
            ss << *expr;
            ++expr;
          }
          std::string op = ss.str();

          // Check if the word parser applies:
          rWordParser_t* parser = config.parserMap.find(op);

          // Evaluate the meaning of this operator in the following order:
          // 1. Is there a word parser for it?
          // 2. Is it a valid operator?
          // 3. Is there a character parser for its first character?
          if (parser) {
            // Parse reserved operators:
            try {
              parser(expr, &expr, &data);
            } catch (...) {
              rpnBuilder::cleanRPN(&data.rpn);
              throw;
            }
          } else if (data.opp.exists(op)) {
            data.handle_op(op);
          } else if ((parser=config.parserMap.find(op[0]))) {
            expr = start+1;
            try {
              parser(expr, &expr, &data);
            } catch (...) {
              rpnBuilder::cleanRPN(&data.rpn);
              throw;
            }
          } else {
            rpnBuilder::cleanRPN(&data.rpn);
            throw syntax_error("Invalid operator: " + op);
          }
        }
      }
    }
    // Ignore spaces but stop on delimiter if not inside brackets.
    while (*expr && isspace(*expr)
           && (data.bracketLevel || !strchr(delim, *expr))) ++expr;
  }

  // Check for syntax errors (excess of operators i.e. 10 + + -1):
  if (data.lastTokenWasUnary) {
    rpnBuilder::cleanRPN(&data.rpn);
    throw syntax_error("Expected operand after unary operator `" + data.opStack.top() + "`");
  }

  std::string cur_op;
  while (!data.opStack.empty()) {
    cur_op = normalize_op(data.opStack.top());
    data.rpn.push(new Token<std::string>(cur_op, OP));
    data.opStack.pop();
  }

  // In case one of the custom parsers left an empty expression:
  if (data.rpn.size() == 0) data.rpn.push(new TokenNone());
  if (rest) *rest = expr;
  return data.rpn;
}

packToken calculator::calculate(const char* expr, TokenMap vars,
                                const char* delim, const char** rest) {
  // Convert to RPN with Dijkstra's Shunting-yard algorithm.
  RAII_TokenQueue_t rpn = calculator::toRPN(expr, vars, delim, rest);

  TokenBase* ret = calculator::calculate(rpn, vars);

  return packToken(resolve_reference(ret));
}

void cleanStack(std::stack<TokenBase*> st) {
  while (st.size() > 0) {
    delete resolve_reference(st.top());
    st.pop();
  }
}

TokenBase* calculator::calculate(const TokenQueue_t& rpn, TokenMap scope,
                                 const Config_t& config) {
  evaluationData data(rpn, scope, config.opMap);

  // Evaluate the expression in RPN form.
  std::stack<TokenBase*> evaluation;
  while (!data.rpn.empty()) {
    TokenBase* base = data.rpn.front()->clone();
    data.rpn.pop();

    // Operator:
    if (base->type == OP) {
      data.op = static_cast<Token<std::string>*>(base)->val;
      delete base;

      /* * * * * Resolve operands Values and References: * * * * */

      if (evaluation.size() < 2) {
        cleanStack(evaluation);
        throw std::domain_error("Invalid equation.");
      }
      TokenBase* r_token = evaluation.top(); evaluation.pop();
      TokenBase* l_token = evaluation.top(); evaluation.pop();

      if (r_token->type & REF) {
        data.right.reset(static_cast<RefToken*>(r_token));
        r_token = data.right->resolve(&data.scope);
      } else if (r_token->type == VAR) {
        packToken key = static_cast<Token<std::string>*>(r_token)->val;
        data.right.reset(new RefToken(key));
      } else {
        data.right.reset(new RefToken());
      }

      if (l_token->type & REF) {
        data.left.reset(static_cast<RefToken*>(l_token));
        l_token = data.left->resolve(&data.scope);
      } else if (l_token->type == VAR) {
        packToken key = static_cast<Token<std::string>*>(l_token)->val;
        data.left.reset(new RefToken(key));
      } else {
        data.left.reset(new RefToken());
      }

      if (l_token->type == FUNC && data.op == "()") {
        // * * * * * Resolve Function Calls: * * * * * //

        Function* l_func = static_cast<Function*>(l_token);

        // Collect the parameter tuple:
        Tuple right;
        if (r_token->type == TUPLE) {
          right = *static_cast<Tuple*>(r_token);
        } else {
          right = Tuple(r_token);
        }
        delete r_token;

        packToken _this;
        if (data.left->origin->type != NOTYPE) {
          _this = data.left->origin;
        } else {
          _this = data.scope;
        }

        // Execute the function:
        packToken ret;
        try {
          ret = Function::call(_this, l_func, &right, data.scope);
        } catch (...) {
          cleanStack(evaluation);
          delete l_func;
          throw;
        }

        delete l_func;
        evaluation.push(ret->clone());
      } else {
        // * * * * * Resolve All Other Operations: * * * * * //

        data.opID = Operation::build_mask(l_token->type, r_token->type);
        packToken l_pack(l_token);
        packToken r_pack(r_token);
        TokenBase* result = 0;

        try {
          // Resolve the operation:
          result = exec_operation(l_pack, r_pack, &data, data.op);
          if (!result) {
            result = exec_operation(l_pack, r_pack, &data, ANY_OP);
          }
        } catch (...) {
          cleanStack(evaluation);
          throw;
        }

        if (result) {
          evaluation.push(result);
        } else {
          cleanStack(evaluation);
          throw undefined_operation(data.op, l_pack, r_pack);
        }
      }
    } else if (base->type == VAR) {  // Variable
      packToken* value = NULL;
      std::string key = static_cast<Token<std::string>*>(base)->val;

      value = data.scope.find(key);

      if (value) {
        TokenBase* copy = (*value)->clone();
        evaluation.push(new RefToken(key, copy));
        delete base;
      } else {
        evaluation.push(base);
      }
    } else {
      evaluation.push(base);
    }
  }

  return evaluation.top();
}

/* * * * * Non Static Functions * * * * */

calculator::~calculator() {
  rpnBuilder::cleanRPN(&this->RPN);
}

calculator::calculator(const calculator& calc) {
  TokenQueue_t _rpn = calc.RPN;

  // Deep copy the token list, so everything can be
  // safely deallocated:
  while (!_rpn.empty()) {
    TokenBase* base = _rpn.front();
    _rpn.pop();
    this->RPN.push(base->clone());
  }
}

// Work as a sub-parser:
// - Stops at delim or '\0'
// - Returns the rest of the string as char* rest
calculator::calculator(const char* expr, TokenMap vars, const char* delim,
                       const char** rest, const Config_t& config) {
  this->RPN = calculator::toRPN(expr, vars, delim, rest, config);
}

void calculator::compile(const char* expr, TokenMap vars, const char* delim,
                         const char** rest) {
  // Make sure it is empty:
  rpnBuilder::cleanRPN(&this->RPN);

  this->RPN = calculator::toRPN(expr, vars, delim, rest, Config());
}

packToken calculator::eval(TokenMap vars, bool keep_refs) const {
  TokenBase* value = calculate(this->RPN, vars, Config());
  packToken p = packToken(value->clone());
  if (keep_refs) {
    return packToken(value);
  } else {
    return packToken(resolve_reference(value));
  }
}

calculator& calculator::operator=(const calculator& calc) {
  // Make sure the RPN is empty:
  rpnBuilder::cleanRPN(&this->RPN);

  // Deep copy the token list, so everything can be
  // safely deallocated:
  TokenQueue_t _rpn = calc.RPN;
  while (!_rpn.empty()) {
    TokenBase* base = _rpn.front();
    _rpn.pop();
    this->RPN.push(base->clone());
  }
  return *this;
}

/* * * * * For Debug Only * * * * */

std::string calculator::str() const {
  return str(this->RPN);
}

std::string calculator::str(TokenQueue_t rpn) {
  std::stringstream ss;

  ss << "calculator { RPN: [ ";
  while (rpn.size()) {
    ss << packToken(resolve_reference(rpn.front()->clone())).str();
    rpn.pop();

    ss << (rpn.size() ? ", ":"");
  }
  ss << " ] }";
  return ss.str();
}
