
#ifndef CONTAINERS_H_
#define CONTAINERS_H_

#include <map>
#include <list>
#include <vector>
#include <string>
#include <memory>

template <typename T>
class Container {
 protected:
  std::shared_ptr<T> ref;

 public:
  Container() : ref(std::make_shared<T>()) {}
  Container(const T& t) : ref(std::make_shared<T>(t)) {}

 public:
  operator T*() const { return ref.get(); }
  friend bool operator==(Container<T> first, Container<T> second) {
    return first.ref == second.ref;
  }
};

class Iterator;

struct Iterable : public TokenBase {
  virtual ~Iterable() {}
  Iterable() {}
  Iterable(tokType_t type) : TokenBase(type) {}

  virtual Iterator* getIterator() const = 0;
};

// Iterator super class.
struct Iterator : public Iterable {
  Iterator() : Iterable(IT) {}
  virtual ~Iterator() {}
  // Return the next position of the iterator.
  // When it reaches the end it should return NULL
  // and reset the iterator automatically.
  virtual packToken* next() = 0;
  virtual void reset() = 0;

  Iterator* getIterator() const;
};

class TokenMap;
typedef std::map<std::string, packToken> TokenMap_t;

struct MapData_t {
  TokenMap_t map;
  TokenMap* parent;
  MapData_t(TokenMap* p);
  MapData_t(const MapData_t& other);
  ~MapData_t();

  MapData_t& operator=(const MapData_t& other);
};

struct TokenMap : public Container<MapData_t>, public Iterable {
  // Static factories:
  static TokenMap empty;
  static TokenMap& base_map();
  static TokenMap& default_global();
  static packToken default_constructor(TokenMap scope);

 public:
  // Attribute getters for the `MapData_t` content:
  TokenMap_t& map() const { return ref->map; }
  TokenMap* parent() const { return ref->parent; }

 public:
  // Implement the Iterable Interface:
  struct MapIterator : public Iterator {
    const TokenMap_t& map;
    TokenMap_t::const_iterator it = map.begin();
    packToken last;

    MapIterator(const TokenMap_t& map) : map(map) {}

    packToken* next();
    void reset();

    TokenBase* clone() const {
      return new MapIterator(*this);
    }
  };

  Iterator* getIterator() const {
    return new MapIterator(map());
  }

 public:
  TokenMap(TokenMap* parent = &TokenMap::base_map())
          : Container(parent), Iterable(MAP) {
    // For the TokenBase super class
    this->type = MAP;
  }
  TokenMap(const TokenMap& other) : Container(other) {
    this->type = MAP;
  }

  virtual ~TokenMap() {}

 public:
  // Implement the TokenBase abstract class
  TokenBase* clone() const {
    return new TokenMap(*this);
  }

 public:
  packToken* find(const std::string& key);
  const packToken* find(const std::string& key) const;
  TokenMap* findMap(const std::string& key);
  void assign(std::string key, TokenBase* value);
  void insert(std::string key, TokenBase* value);

  TokenMap getChild();

  packToken& operator[](const std::string& str);

  void erase(std::string key);
};

// Build a TokenMap which is a child of default_global()
struct GlobalScope : public TokenMap {
  GlobalScope() : TokenMap(&TokenMap::default_global()) {}
};

typedef std::vector<packToken> TokenList_t;

struct TokenList : public Container<TokenList_t>, public Iterable {
  static packToken default_constructor(TokenMap scope);

 public:
  // Attribute getter for the `TokenList_t` content:
  TokenList_t& list() const { return *ref; }

 public:
  struct ListIterator : public Iterator {
    TokenList_t* list;
    uint64_t i = 0;

    ListIterator(TokenList_t* list) : list(list) {}

    packToken* next();
    void reset();

    TokenBase* clone() const {
      return new ListIterator(*this);
    }
  };

  Iterator* getIterator() const {
    return new ListIterator(&list());
  }

 public:
  TokenList() { this->type = LIST; }
  virtual ~TokenList() {}

  packToken& operator[](const uint64_t idx) const {
    if (list().size() <= idx) {
      throw std::out_of_range("List index out of range!");
    }
    return list()[idx];
  }

  void push(packToken val) const { list().push_back(val); }
  packToken pop() const {
    packToken back = list().back();
    list().pop_back();
    return back;
  }

 public:
  // Implement the TokenBase abstract class
  TokenBase* clone() const {
    return new TokenList(*this);
  }
};

class Tuple : public TokenList {
 public:
  Tuple() { this->type = TUPLE; }
  Tuple(const TokenBase* first) {
    this->type = TUPLE;
    list().push_back(packToken(first->clone()));
  }
  Tuple(const packToken first) : Tuple(first.token()) {}

  Tuple(const TokenBase* first, const TokenBase* second) {
    this->type = TUPLE;
    list().push_back(packToken(first->clone()));
    list().push_back(packToken(second->clone()));
  }
  Tuple(const packToken first, const packToken second)
       : Tuple(first.token(), second.token()) {}

 public:
  // Implement the TokenBase abstract class
  TokenBase* clone() const {
    return new Tuple(*this);
  }
};

// This Special Tuple is to be used only as syntactic sugar, and
// constructed only with the operator `:`, i.e.:
// - passing key-word arguments: func(1, 2, optional_arg:10)
// - slicing lists or strings: my_list[2:10:2] (not implemented)
//
// STuple means one of:
// - Special Tuple, Syntactic Tuple or System Tuple
//
// I haven't decided yet. Suggestions accepted.
class STuple : public Tuple {
 public:
  STuple() { this->type = STUPLE; }
  STuple(const TokenBase* first) {
    this->type = STUPLE;
    list().push_back(packToken(first->clone()));
  }
  STuple(const packToken first) : STuple(first.token()) {}

  STuple(const TokenBase* first, const TokenBase* second) {
    this->type = STUPLE;
    list().push_back(packToken(first->clone()));
    list().push_back(packToken(second->clone()));
  }
  STuple(const packToken first, const packToken second)
       : STuple(first.token(), second.token()) {}

 public:
  // Implement the TokenBase abstract class
  TokenBase* clone() const {
    return new STuple(*this);
  }
};

#endif  // CONTAINERS_H_
