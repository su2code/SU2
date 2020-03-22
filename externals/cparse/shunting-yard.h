#ifndef SHUNTING_YARD_H_
#define SHUNTING_YARD_H_
#include <iostream>

#include <map>
#include <stack>
#include <string>
#include <queue>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <memory>
#include <utility>

#if defined CODI_REVERSE_TYPE
#include "codi.hpp"
typedef codi::RealReverse doubleType;
#elif defined CODI_FORWARD_TYPE
#include "codi.hpp"
typedef codi::RealForward doubleType;
#else
typedef double doubleType;
#endif

/*
 * About tokType enum:
 *
 * The 3 left most bits (0x80, 0x40 and 0x20) of the Token Type
 * are reserved for denoting Numerals, Iterators and References.
 * If you want to define your own type please mind this bits.
 */
typedef uint8_t tokType_t;
typedef uint64_t opID_t;
enum tokType {
  // Internal types:
  NOTYPE, OP, UNARY, VAR,

  // Base types:
  // Note: The mask system accepts at most 29 (32-3) different base types.
  STR, FUNC,

  // Numerals:
  NUM = 0x20,   // Everything with the bit 0x20 set is a number.
  REAL = 0x21,  // == 0x20 + 0x1 => Real numbers.
  INT = 0x22,   // == 0x20 + 0x2 => Integral numbers.
  BOOL = 0x23,  // == 0x20 + 0x3 => Boolean Type.

  // Complex types:
  IT = 0x40,      // Everything with the bit 0x40 set are iterators.
  LIST = 0x41,    // == 0x40 + 0x01 => Lists are iterators.
  TUPLE = 0x42,   // == 0x40 + 0x02 => Tuples are iterators.
  STUPLE = 0x43,  // == 0x40 + 0x03 => ArgTuples are iterators.
  MAP = 0x44,     // == 0x40 + 0x04 => Maps are Iterators

  // References are internal tokens used by the calculator:
  REF = 0x80,

  // Mask used when defining operations:
  ANY_TYPE = 0xFF
};

#define ANY_OP ""

struct TokenBase {
  tokType_t type;

  virtual ~TokenBase() {}
  TokenBase() {}
  TokenBase(tokType_t type) : type(type) {}

  virtual TokenBase* clone() const = 0;
};

template<class T> class Token : public TokenBase {
 public:
  T val;
  Token(T t, tokType_t type) : TokenBase(type), val(t) {}
  virtual TokenBase* clone() const {
    return new Token(*this);
  }
};

struct TokenNone : public TokenBase {
  TokenNone() : TokenBase(NOTYPE) {}
  virtual TokenBase* clone() const {
    return new TokenNone(*this);
  }
};

struct TokenUnary : public TokenBase {
  TokenUnary() : TokenBase(UNARY) {}
  virtual TokenBase* clone() const {
    return new TokenUnary(*this);
  }
};

class packToken;
typedef std::queue<TokenBase*> TokenQueue_t;
class OppMap_t {
  // Set of operators that should be evaluated from right to left:
  std::set<std::string> RtoL;
  // Map of operators precedence:
  std::map<std::string, int> pr_map;

 public:
  OppMap_t() {
    // These operations are hard-coded inside the calculator,
    // thus their precedence should always be defined:
    pr_map["[]"] = -1; pr_map["()"] = -1;
    pr_map["["] = 0x7FFFFFFF; pr_map["("] = 0x7FFFFFFF; pr_map["{"] = 0x7FFFFFFF;
    RtoL.insert("=");
  }

  void add(const std::string& op, int precedence) {
    if (precedence < 0) {
      RtoL.insert(op);
      precedence = -precedence;
    }

    pr_map[op] = precedence;
  }

  void addUnary(const std::string& op, int precedence) {
    add("L"+op, precedence);

    // Also add a binary operator with same precedence so
    // it is possible to verify if an op exists just by checking
    // the binary set of operators:
    if (!exists(op)) {
      add(op, precedence);
    }
  }

  void addRightUnary(const std::string& op, int precedence) {
    add("R"+op, precedence);

    // Also add a binary operator with same precedence so
    // it is possible to verify if an op exists just by checking
    // the binary set of operators:
    if (!exists(op)) {
      add(op, precedence);
    } else {
      // Note that using a unary and binary operators with
      // the same left operand is ambiguous and that the unary
      // operator will take precedence.
      //
      // So only do it if you know the expected left operands
      // have distinct types.
    }
  }

  int prec(const std::string& op) const { return pr_map.at(op); }
  bool assoc(const std::string& op) const { return RtoL.count(op); }
  bool exists(const std::string& op) const { return pr_map.count(op); }
};

class TokenMap;
class TokenList;
class Tuple;
class STuple;
class Function;
#include "./packToken.h"

// Define the Tuple, TokenMap and TokenList classes:
#include "./containers.h"

// Define the `Function` class
// as well as some built-in functions:
#include "./functions.h"

// This struct was created to expose internal toRPN() variables
// to custom parsers, in special to the rWordParser_t functions.
struct rpnBuilder {
  TokenQueue_t rpn;
  std::stack<std::string> opStack;
  uint8_t lastTokenWasOp = true;
  bool lastTokenWasUnary = false;
  TokenMap scope;
  const OppMap_t& opp;

  // Used to make sure the expression won't
  // end inside a bracket evaluation just because
  // found a delimiter like '\n' or ')'
  uint32_t bracketLevel = 0;

  rpnBuilder(TokenMap scope, const OppMap_t& opp) : scope(scope), opp(opp) {}

 public:
  static void cleanRPN(TokenQueue_t* rpn);

 public:
  void handle_op(const std::string& op);
  void handle_token(TokenBase* token);
  void open_bracket(const std::string& bracket);
  void close_bracket(const std::string& bracket);

  // * * * * * Static parsing helpers: * * * * * //

  // Check if a character is the first character of a variable:
  static inline bool isvarchar(const char c) {
    return isalpha(c) || c == '_' || c == '@';
  }

  static inline std::string parseVar(const char* expr, const char** rest = 0) {
    std::stringstream ss;
    ss << *expr;
    ++expr;
    while (rpnBuilder::isvarchar(*expr) || isdigit(*expr)) {
      ss << *expr;
      ++expr;
    }
    if (rest) *rest = expr;
    return ss.str();
  }

 private:
  void handle_opStack(const std::string& op);
  void handle_binary(const std::string& op);
  void handle_left_unary(const std::string& op);
  void handle_right_unary(const std::string& op);
};

class RefToken;
class opMap_t;
struct evaluationData {
  TokenQueue_t rpn;
  TokenMap scope;
  const opMap_t& opMap;

  std::unique_ptr<RefToken> left;
  std::unique_ptr<RefToken> right;

  std::string op;
  opID_t opID;

  evaluationData(TokenQueue_t rpn, TokenMap scope, const opMap_t& opMap)
                : rpn(rpn), scope(scope), opMap(opMap) {}
};

// The reservedWordParser_t is the function type called when
// a reserved word or character is found at parsing time.
typedef void rWordParser_t(const char* expr, const char** rest,
                           rpnBuilder* data);
typedef std::map<std::string, rWordParser_t*> rWordMap_t;
typedef std::map<char, rWordParser_t*> rCharMap_t;

struct parserMap_t {
  rWordMap_t wmap;
  rCharMap_t cmap;

  // Add reserved word:
  void add(const std::string& word, const rWordParser_t* parser) {
    wmap[word] = parser;
  }

  // Add reserved character:
  void add(char c, const rWordParser_t* parser) {
    cmap[c] = parser;
  }

  rWordParser_t* find(const std::string text) {
    rWordMap_t::iterator w_it;

    if ((w_it=wmap.find(text)) != wmap.end()) {
      return w_it->second;
    }

    return 0;
  }

  rWordParser_t* find(char c) {
    rCharMap_t::iterator c_it;

    if ((c_it=cmap.find(c)) != cmap.end()) {
      return c_it->second;
    }

    return 0;
  }
};

// The RefToken keeps information about the context
// in which a variable was originally evaluated
// and allow a final value to be correctly resolved
// afterwards.
class RefToken : public TokenBase {
  packToken original_value;

 public:
  packToken key;
  packToken origin;
  RefToken(packToken k, TokenBase* v, packToken m = packToken::None()) :
    TokenBase(v->type | REF), original_value(v), key(std::forward<packToken>(k)), origin(std::forward<packToken>(m)) {}
  RefToken(packToken k = packToken::None(), packToken v = packToken::None(), packToken m = packToken::None()) :
    TokenBase(v->type | REF), original_value(std::forward<packToken>(v)), key(std::forward<packToken>(k)), origin(std::forward<packToken>(m)) {}

  TokenBase* resolve(TokenMap* localScope = 0) const {
    TokenBase* result = 0;

    // Local variables have no origin == tokType::NONE,
    // thus, require a localScope to be resolved:
    if (origin->type == NOTYPE && localScope) {
      // Get the most recent value from the local scope:
      packToken* r_value = localScope->find(key.asString());
      if (r_value) {
        result = (*r_value)->clone();
      }
    }

    // In last case return the compilation-time value:
    return result ? result : original_value->clone();
  }

  virtual TokenBase* clone() const {
    return new RefToken(*this);
  }
};

struct opSignature_t {
  tokType_t left; std::string op; tokType_t right;
  opSignature_t(const tokType_t L, const std::string op, const tokType_t R)
               : left(L), op(op), right(R) {}
};

class Operation {
 public:
  typedef packToken (*opFunc_t)(const packToken& left, const packToken& right,
                                evaluationData* data);

 public:
  // Use this exception to reject an operation.
  // Without stoping the operation matching process.
  struct Reject : public std::exception {};

 public:
  static inline uint32_t mask(tokType_t type);
  static opID_t build_mask(tokType_t left, tokType_t right);

 private:
  opID_t _mask;
  opFunc_t _exec;

 public:
  Operation(opSignature_t sig, opFunc_t func)
           : _mask(build_mask(sig.left, sig.right)), _exec(func) {}

 public:
  opID_t getMask() const { return _mask; }
  packToken exec(const packToken& left, const packToken& right,
                 evaluationData* data) const {
    return _exec(left, right, data);
  }
};

typedef std::map<tokType_t, TokenMap> typeMap_t;
typedef std::vector<Operation> opList_t;
struct opMap_t : public std::map<std::string, opList_t> {
  void add(const opSignature_t sig, Operation::opFunc_t func) {
    (*this)[sig.op].push_back(Operation(sig, func));
  }

  std::string str() const {
    if (this->size() == 0) return "{}";

    std::string result = "{ ";
    for (const auto& pair : (*this)) {
      result += "\"" + pair.first + "\", ";
    }
    result.pop_back();
    result.pop_back();
    return result + " }";
  }
};

struct Config_t {
  parserMap_t parserMap;
  OppMap_t opPrecedence;
  opMap_t opMap;

  Config_t() {}
  Config_t(parserMap_t p, OppMap_t opp, opMap_t opMap)
          : parserMap(p), opPrecedence(opp), opMap(opMap) {}
};

class calculator {
 public:
  static Config_t& Default();

 public:
  static typeMap_t& type_attribute_map();

 public:
  static packToken calculate(const char* expr, TokenMap vars = &TokenMap::empty,
                             const char* delim = 0, const char** rest = 0);

 public:
  static TokenBase* calculate(const TokenQueue_t& RPN, TokenMap scope,
                              const Config_t& config = Default());
  static TokenQueue_t toRPN(const char* expr, TokenMap vars,
                            const char* delim = 0, const char** rest = 0,
                            Config_t config = Default());

 public:
  // Used to dealloc a TokenQueue_t safely.
  struct RAII_TokenQueue_t;

 protected:
  virtual const Config_t Config() const { return Default(); }

 private:
  TokenQueue_t RPN;

 public:
  virtual ~calculator();
  calculator() { this->RPN.push(new TokenNone()); }
  calculator(const calculator& calc);
  calculator(const char* expr, TokenMap vars = &TokenMap::empty,
             const char* delim = 0, const char** rest = 0,
             const Config_t& config = Default());
  void compile(const char* expr, TokenMap vars = &TokenMap::empty,
               const char* delim = 0, const char** rest = 0);
  packToken eval(TokenMap vars = &TokenMap::empty, bool keep_refs = false) const;

  // Serialization:
  std::string str() const;
  static std::string str(TokenQueue_t rpn);

  // Operators:
  calculator& operator=(const calculator& calc);
};

#endif  // SHUNTING_YARD_H_
