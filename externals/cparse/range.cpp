#include <cstring>

#include "./range.h"

packToken* Range::next() {
    int64_t value = i;
    if ((step > 0 && value >= to) || (step < 0 && value <= to)) {
      i = from;
      return NULL;
    } else {
      i += step;
      last = static_cast<doubleType>(value);
      return &last;
    }
  }

void Range::reset() { i=from; }

/* * * * * Built-in range function: * * * * */

const char* range_args[] = {"from", "to", "step"};
packToken default_range(TokenMap scope) {
  int64_t to, step, from;

  packToken* p_from = scope.find("from");
  packToken* p_to = scope.find("to");
  packToken* p_step = scope.find("step");

  if ((*p_from)->type & NUM) {
    from = p_from->asInt();
  } else if ((*p_from)->type == NOTYPE) {
    throw std::invalid_argument("range() expects at least 1 argument!");
  } else {
    throw std::invalid_argument("range() expects only numbers!");
  }

  if ((*p_to)->type == NOTYPE) {
    to = from;
    from = 0;
    step = 1;
  } else if ((*p_to)->type & NUM) {
    to = p_to->asInt();

    if ((*p_step)->type & NUM) {
      step = p_step->asInt();
    } else if ((*p_step)->type == NOTYPE) {
      step = 1;
    } else {
      throw std::invalid_argument("range() expects only numbers!");
    }
  } else {
    throw std::invalid_argument("range() expects only numbers!");
  }

  return packToken(new Range(from, to, step));
}

/* * * * * Range Startup class * * * * */

struct Range::Startup {
  Startup() {
    TokenMap& global = TokenMap::default_global();
    global["range"] = CppFunction(default_range, 3, range_args, "range");
  }
} iterator_startup;
