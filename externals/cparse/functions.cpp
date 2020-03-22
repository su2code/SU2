#include <string>

#include "./shunting-yard.h"
#include "./functions.h"
#include "./shunting-yard-exceptions.h"

/* * * * * class Function * * * * */
packToken Function::call(packToken _this, const Function* func,
                         TokenList* args, TokenMap scope) {
  // Build the local namespace:
  TokenMap kwargs;
  TokenMap local = scope.getChild();

  args_t arg_names = func->args();

  TokenList_t::iterator args_it = args->list().begin();
  args_t::const_iterator names_it = arg_names.begin();

  /* * * * * Parse positional arguments: * * * * */

  while (args_it != args->list().end() && names_it != arg_names.end()) {
    // If the positional argument list is over:
    if ((*args_it)->type == STUPLE) break;

    // Else add it to the local namespace:
    local[*names_it] = *args_it;

    ++args_it;
    ++names_it;
  }

  /* * * * * Parse extra positional arguments: * * * * */

  TokenList arglist;
  for (; args_it != args->list().end(); ++args_it) {
    // If there is a keyword argument:
    if ((*args_it)->type == STUPLE) break;
    // Else add it to arglist:
    arglist.list().push_back(*args_it);
  }

  /* * * * * Parse keyword arguments: * * * * */

  for (; args_it != args->list().end(); ++args_it) {
    packToken& arg = *args_it;

    if (arg->type != STUPLE) {
      throw syntax_error("Positional argument follows keyword argument");
    }

    STuple* st = static_cast<STuple*>(arg.token());

    if (st->list().size() != 2) {
      throw syntax_error("Keyword tuples must have exactly 2 items!");
    }

    if (st->list()[0]->type != STR) {
      throw syntax_error("Keyword first argument should be of type string!");
    }

    // Save it:
    std::string key = st->list()[0].asString();
    packToken& value = st->list()[1];
    kwargs[key] = value;
  }

  /* * * * * Set missing positional arguments: * * * * */

  for (; names_it != arg_names.end(); ++names_it) {
    // If not set by a keyword argument:
    auto kw_it = kwargs.map().find(*names_it);
    if (kw_it == kwargs.map().end()) {
      local[*names_it] = packToken::None();
    } else {
      local[*names_it] = kw_it->second;
      kwargs.map().erase(kw_it);
    }
  }

  /* * * * * Set built-in variables: * * * * */

  local["this"] = _this;
  local["args"] = arglist;
  local["kwargs"] = kwargs;

  return func->exec(local);
}

/* * * * * class CppFunction * * * * */

CppFunction::CppFunction(packToken (*func)(TokenMap), const args_t args,
                         std::string name)
                         : func(func), _args(args) {
  this->_name = name;
}

CppFunction::CppFunction(packToken (*func)(TokenMap), unsigned int nargs,
                         const char** args, std::string name)
                         : func(func) {
  this->_name = name;
  // Add all strings to args list:
  for (uint32_t i = 0; i < nargs; ++i) {
    this->_args.push_back(args[i]);
  }
}

// Build a function with no named args:
CppFunction::CppFunction(packToken (*func)(TokenMap), std::string name)
                         : func(func) {
  this->_name = name;
}
