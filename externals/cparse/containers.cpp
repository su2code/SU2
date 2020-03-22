#include <string>

#include "./shunting-yard.h"

#include "./containers.h"
#include "./functions.h"

/* * * * * Initialize TokenMap * * * * */

// Using the "Construct On First Use Idiom"
// to avoid the "static initialization order fiasco",
// for more information read:
//
// - https://isocpp.org/wiki/faq/ctors#static-init-order
//
TokenMap& TokenMap::base_map() {
  static TokenMap _base_map(0);
  return _base_map;
}

TokenMap& TokenMap::default_global() {
  static TokenMap global_map(base_map());
  return global_map;
}

packToken TokenMap::default_constructor(TokenMap scope) {
  return scope["kwargs"];
}

TokenMap TokenMap::empty = TokenMap(&default_global());


/* * * * * Iterator functions * * * * */

Iterator*  Iterator::getIterator() const {
  return static_cast<Iterator*>(this->clone());
}

/* * * * * TokenMap iterator implemented functions * * * * */

packToken* TokenMap::MapIterator::next() {
  if (it != map.end()) {
    last = packToken(it->first);
    ++it;
    return &last;
  } else {
    it = map.begin();
    return NULL;
  }
}

void TokenMap::MapIterator::reset() { it = map.begin(); }

/* * * * * TokenList functions: * * * * */

packToken TokenList::default_constructor(TokenMap scope) {
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

/* * * * * TokenList iterator implemented functions * * * * */

packToken* TokenList::ListIterator::next() {
  if (i < list->size()) {
    return &list->at(i++);
  } else {
    i = 0;
    return NULL;
  }
}

void TokenList::ListIterator::reset() { i = 0; }

/* * * * * MapData_t struct: * * * * */

MapData_t::MapData_t(TokenMap* p) : parent(p ? new TokenMap(*p) : 0) {}
MapData_t::MapData_t(const MapData_t& other) {
  map = other.map;
  if (other.parent) {
    parent = new TokenMap(*(other.parent));
  } else {
    parent = 0;
  }
}
MapData_t::~MapData_t() { if (parent) delete parent; }

MapData_t& MapData_t::operator=(const MapData_t& other) {
  if (this != &other) {
    if (parent) delete parent;
    map = other.map;
    parent = other.parent;
  }
  return *this;
}

/* * * * * TokenMap Class: * * * * */

packToken* TokenMap::find(const std::string& key) {
  TokenMap_t::iterator it = map().find(key);

  if (it != map().end()) {
    return &it->second;
  } else if (parent()) {
    return parent()->find(key);
  } else {
    return 0;
  }
}

const packToken* TokenMap::find(const std::string& key) const {
  TokenMap_t::const_iterator it = map().find(key);

  if (it != map().end()) {
    return &it->second;
  } else if (parent()) {
    return parent()->find(key);
  } else {
    return 0;
  }
}

TokenMap* TokenMap::findMap(const std::string& key) {
  TokenMap_t::iterator it = map().find(key);

  if (it != map().end()) {
    return this;
  } else if (parent()) {
    return parent()->findMap(key);
  } else {
    return 0;
  }
}

void TokenMap::assign(std::string key, TokenBase* value) {
  if (value) {
    value = value->clone();
  } else {
    throw std::invalid_argument("TokenMap assignment expected a non NULL argument as value!");
  }

  packToken* variable = find(key);

  if (variable) {
    (*variable) = packToken(value);
  } else {
    map()[key] = packToken(value);
  }
}

void TokenMap::insert(std::string key, TokenBase* value) {
  (*this)[key] = packToken(value->clone());
}

packToken& TokenMap::operator[](const std::string& key) {
  return map()[key];
}

TokenMap TokenMap::getChild() {
  return TokenMap(this);
}

void TokenMap::erase(std::string key) {
  map().erase(key);
}
