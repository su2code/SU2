#include <cmath>

namespace builtin_functions {

/* * * * * Built-in Functions: * * * * */

packToken default_print(TokenMap scope) {
  // Get the argument list:
  TokenList list = scope["args"].asList();

  bool first = true;
  for (packToken item : list.list()) {
    if (first) {
      first = false;
    } else {
      std::cout << " ";
    }

    if (item->type == STR) {
      std::cout << item.asString();
    } else {
      std::cout << item.str();
    }
  }

  std::cout << std::endl;

  return packToken::None();
}

packToken default_sum(TokenMap scope) {
  // Get the arguments:
  TokenList list = scope["args"].asList();

  if (list.list().size() == 1 && list.list().front()->type == LIST) {
    list = list.list().front().asList();
  }

  doubleType sum = 0;
  for (packToken num : list.list()) {
    sum += num.asDouble();
  }

  return sum;
}

packToken default_eval(TokenMap scope) {
  std::string code = scope["value"].asString();
  // Evaluate it as a calculator expression:
  return calculator::calculate(code.c_str(), scope);
}

packToken default_float(TokenMap scope) {
  packToken tok = scope["value"];
  if (tok->type & NUM) return tok.asDouble();

  // Convert it to doubleType:
  char* rest;
  const std::string& str = tok.asString();
  errno = 0;
  doubleType ret = strtod(str.c_str(), &rest);

  if (str == rest) {
    throw std::runtime_error("Could not convert \"" + str + "\" to float!");
  } else if (errno) {
    std::range_error("Value too big or too small to fit a Double!");
  }
  return ret;
}

packToken default_int(TokenMap scope) {
  packToken tok = scope["value"];
  if (tok->type & NUM) return tok.asInt();

  // Convert it to doubleType:
  char* rest;
  const std::string& str = tok.asString();
  errno = 0;
  int64_t ret = strtol(str.c_str(), &rest, 10);

  if (str == rest) {
    throw std::runtime_error("Could not convert \"" + str + "\" to integer!");
  } else if (errno) {
    std::range_error("Value too big or too small to fit an Integer!");
  }
  return ret;
}

packToken default_str(TokenMap scope) {
  // Return its string representation:
  packToken tok = scope["value"];
  if (tok->type == STR) return tok;
  return tok.str();
}

packToken default_type(TokenMap scope) {
  packToken tok = scope["value"];
  packToken* p_type;

  switch (tok->type) {
  case NOTYPE: return "none";
  case VAR: return "variable";
  case REAL: return "real";
  case INT: return "integer";
  case BOOL: return "boolean";
  case STR: return "string";
  case FUNC: return "function";
  case IT: return "iterable";
  case TUPLE: return "tuple";
  case STUPLE: return "argument tuple";
  case LIST: return "list";
  case MAP:
    p_type = tok.asMap().find("__type__");
    if (p_type && (*p_type)->type == STR) {
      return *p_type;
    } else  {
      return "map";
    }
  default: return "unknown_type";
  }
}

packToken default_sqrt(TokenMap scope) {
  // Get a single argument:
  doubleType number = sqrt(scope["num"].asDouble());

  return number;
}
packToken default_sin(TokenMap scope) {
  // Get a single argument:
  doubleType number = sin(scope["num"].asDouble());

  return number;
}
packToken default_cos(TokenMap scope) {
  // Get a single argument:
  doubleType number = cos(scope["num"].asDouble());

  return number;
}
packToken default_tan(TokenMap scope) {
  // Get a single argument:
  doubleType number = tan(scope["num"].asDouble());

  return number;
}
packToken default_abs(TokenMap scope) {
  // Get a single argument:
  doubleType number =  abs(scope["num"].asDouble());

  return number;
}

const args_t pow_args = {"number", "exp"};
packToken default_pow(TokenMap scope) {
  // Get two arguments:
  doubleType number = scope["number"].asDouble();
  doubleType exp = pow(number, scope["exp"].asDouble());

  return exp;
}

/* * * * * default constructor functions * * * * */

packToken default_list(TokenMap scope) {
  // Get the arguments:
  TokenList list = scope["args"].asList();

  // If the only argument is iterable:
  if (list.list().size() == 1 && list.list()[0]->type & IT) {
    TokenList new_list;
    Iterator* it = static_cast<Iterable*>(list.list()[0].token())->getIterator();

    packToken* next = it->next();
    while (next) {
      new_list.list().push_back(*next);
      next = it->next();
    }

    delete it;
    return new_list;
  } else {
    return list;
  }
}

packToken default_map(TokenMap scope) {
  return scope["kwargs"];
}

/* * * * * Object inheritance tools: * * * * */

packToken default_extend(TokenMap scope) {
  packToken tok = scope["value"];

  if (tok->type == MAP) {
    return tok.asMap().getChild();
  } else {
    throw std::runtime_error(tok.str() + " is not extensible!");
  }
}

// Example of replacement function for packToken::str():
std::string packToken_str(const TokenBase* base, uint32_t nest) {
  const Function* func;

  // Find the TokenMap with the type specific functions
  // for the type of the base token:
  const TokenMap* typeFuncs;
  if (base->type == MAP) {
    typeFuncs = static_cast<const TokenMap*>(base);
  } else {
    typeFuncs = &calculator::type_attribute_map()[base->type];
  }

  // Check if this type has a custom stringify function:
  const packToken* p_func = typeFuncs->find("__str__");
  if (p_func && (*p_func)->type == FUNC) {
    // Return the result of this function passing the
    // nesting level as first (and only) argument:
    func = p_func->asFunc();
    packToken _this = packToken(base->clone());
    TokenList args;
    args.push(static_cast<int64_t>(nest));
    return Function::call(_this, func, &args, TokenMap()).asString();
  }

  // Return "" to ask for the normal `packToken::str()`
  // function to complete the job.
  return "";
}

struct Startup {
  Startup() {
    TokenMap& global = TokenMap::default_global();

    global["print"] = CppFunction(&default_print, "print");
    global["sum"] = CppFunction(&default_sum, "sum");
    global["sqrt"] = CppFunction(&default_sqrt, {"num"}, "sqrt");
    global["sin"] = CppFunction(&default_sin, {"num"}, "sin");
    global["cos"] = CppFunction(&default_cos, {"num"}, "cos");
    global["tan"] = CppFunction(&default_tan, {"num"}, "tan");
    global["abs"] = CppFunction(&default_abs, {"num"}, "abs");
    global["pow"] = CppFunction(&default_pow, pow_args, "pow");
    global["float"] = CppFunction(&default_float, {"value"}, "float");
    global["real"] = CppFunction(&default_float, {"value"}, "real");
    global["int"] = CppFunction(&default_int, {"value"}, "int");
    global["str"] = CppFunction(&default_str, {"value"}, "str");
    global["eval"] = CppFunction(&default_eval, {"value"}, "eval");
    global["type"] = CppFunction(&default_type, {"value"}, "type");
    global["extend"] = CppFunction(&default_extend, {"value"}, "extend");

    // Default constructors:
    global["list"] = CppFunction(&default_list, "list");
    global["map"] = CppFunction(&default_map, "map");

    // Set the custom str function to `packToken_str()`
    packToken::str_custom() = packToken_str;
  }
} __CPARSE_STARTUP;

}  // namespace builtin_functions
