
namespace builtin_typeSpecificFunctions {

/* * * * * MAP Type built-in functions * * * * */

const args_t map_pop_args = {"key", "default"};
packToken map_pop(TokenMap scope) {
  TokenMap map = scope["this"].asMap();
  std::string key = scope["key"].asString();

  // Check if the item is available and remove it:
  if (map.map().count(key)) {
    packToken value = map[key];
    map.erase(key);
    return value;
  }

  // If not available return the default value or None
  packToken* def = scope.find("default");
  if (def) {
    return *def;
  } else {
    return packToken::None();
  }
}

packToken map_len(TokenMap scope) {
  TokenMap map = scope.find("this")->asMap();
  return map.map().size();
}

packToken default_instanceof(TokenMap scope) {
  TokenMap _super = scope["value"].asMap();
  TokenMap* _this = scope["this"].asMap().parent();

  TokenMap* parent = _this;
  while (parent) {
    if ((*parent) == _super) {
      return true;
    }

    parent = parent->parent();
  }

  return false;
}

/* * * * * LIST Type built-in functions * * * * */

const args_t push_args = {"item"};
packToken list_push(TokenMap scope) {
  packToken* list = scope.find("this");
  packToken* token = scope.find("item");

  // If "this" is not a list it will throw here:
  list->asList().list().push_back(*token);

  return *list;
}

const args_t list_pop_args = {"pos"};
packToken list_pop(TokenMap scope) {
  TokenList list = scope.find("this")->asList();
  packToken* token = scope.find("pos");

  int64_t pos;

  if ((*token)->type & NUM) {
    pos = token->asInt();

    // So that pop(-1) is the same as pop(last_idx):
    if (pos < 0) pos = list.list().size()-pos;
  } else {
    pos = list.list().size()-1;
  }

  packToken result = list.list()[pos];

  // Erase the item from the list:
  // Note that this operation is optimal if pos == list.size()-1
  list.list().erase(list.list().begin() + pos);

  return result;
}

packToken list_len(TokenMap scope) {
  TokenList list = scope.find("this")->asList();
  return list.list().size();
}

packToken list_join(TokenMap scope) {
  TokenList list = scope["this"].asList();
  std::string chars = scope["chars"].asString();
  std::stringstream result;

  std::vector<packToken>::const_iterator it = list.list().begin();
  result << it->asString();
  for (++it; it != list.list().end(); ++it) {
    result << chars << it->asString();
  }

  return result.str();
}

/* * * * * STR Type built-in functions * * * * */

packToken string_len(TokenMap scope) {
  std::string str = scope["this"].asString();
  return static_cast<int64_t>(str.size());
}

packToken string_lower(TokenMap scope) {
  std::string str = scope["this"].asString();
  std::string out;
  for (char c : str) {
    out.push_back(tolower(c));
  }
  return out;
}

packToken string_upper(TokenMap scope) {
  std::string str = scope["this"].asString();
  std::string out;
  for (char c : str) {
    out.push_back(toupper(c));
  }
  return out;
}

packToken string_strip(TokenMap scope) {
  std::string str = scope["this"].asString();

  std::string::const_iterator it = str.begin();
  while (it != str.end() && isspace(*it)) ++it;

  std::string::const_reverse_iterator rit = str.rbegin();
  while (rit.base() != it && isspace(*rit)) ++rit;

  return std::string(it, rit.base());
}

packToken string_split(TokenMap scope) {
  TokenList list;
  std::string str = scope["this"].asString();
  std::string split_chars = scope["chars"].asString();

  // Split the string:
  size_t start = 0;
  size_t i = str.find(split_chars, 0);
  size_t size = split_chars.size();
  while (i < str.size()) {
    // Add a new item:
    list.push(std::string(str, start, i-start));
    // Resume search:
    start = i + size;
    i = str.find(split_chars, start);
  }

  // Add a new item:
  list.push(std::string(str, start, str.size()-start));

  return list;
}

/* * * * * Type-Specific Functions Startup: * * * * */

struct Startup {
  Startup() {
    TokenMap& base_list = calculator::type_attribute_map()[LIST];
    base_list["push"] = CppFunction(list_push, push_args, "push");
    base_list["pop"] = CppFunction(list_pop, list_pop_args, "pop");
    base_list["len"] = CppFunction(list_len, "len");
    base_list["join"] = CppFunction(list_join, {"chars"}, "join");

    TokenMap& base_str = calculator::type_attribute_map()[STR];
    base_str["len"] = CppFunction(&string_len, "len");
    base_str["lower"] = CppFunction(&string_lower, "lower");
    base_str["upper"] = CppFunction(&string_upper, "upper");
    base_str["strip"] = CppFunction(&string_strip, "strip");
    base_str["split"] = CppFunction(&string_split, {"chars"}, "split");

    TokenMap& base_map = TokenMap::base_map();
    base_map["pop"] = CppFunction(map_pop, map_pop_args, "pop");
    base_map["len"] = CppFunction(map_len, "len");
    base_map["instanceof"] = CppFunction(&default_instanceof,
                                         {"value"}, "instanceof");
  }
} __CPARSE_STARTUP;

}  // namespace builtin_typeSpecificFunctions
