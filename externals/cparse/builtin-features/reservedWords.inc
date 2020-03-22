
namespace builtin_reservedWords {

// Literal Tokens: True, False and None:
packToken trueToken = packToken(true);
packToken falseToken = packToken(false);
packToken noneToken = packToken::None();

void True(const char* expr, const char** rest, rpnBuilder* data) {
  data->handle_token(trueToken->clone());
}

void False(const char* expr, const char** rest, rpnBuilder* data) {
  data->handle_token(falseToken->clone());
}

void None(const char* expr, const char** rest, rpnBuilder* data) {
  data->handle_token(noneToken->clone());
}

void LineComment(const char* expr, const char** rest, rpnBuilder* data) {
  while (*expr && *expr != '\n') ++expr;
  *rest = expr;
}

void SlashStarComment(const char* expr, const char** rest, rpnBuilder* data) {
  while (*expr && !(expr[0] == '*' && expr[1] == '/')) ++expr;
  if (*expr == '\0') {
    throw syntax_error("Unexpected end of file after '/*' comment!");
  }
  // Drop the characters `*/`:
  expr += 2;
  *rest = expr;
}

void KeywordOperator(const char* expr, const char** rest, rpnBuilder* data) {
  // Convert any STuple like `a : 10` to `'a': 10`:
  if (data->rpn.back()->type == VAR) {
    data->rpn.back()->type = STR;
  }
  data->handle_op(":");
}

void DotOperator(const char* expr, const char** rest, rpnBuilder* data) {
  data->handle_op(".");

  while (*expr && isspace(*expr)) ++expr;

  // If it did not find a valid variable name after it:
  if (!rpnBuilder::isvarchar(*expr)) {
    throw syntax_error("Expected variable name after '.' operator");
  }

  // Parse the variable name and save it as a string:
  std::string key = rpnBuilder::parseVar(expr, rest);
  data->handle_token(new Token<std::string>(key, STR));
}

struct Startup {
  Startup() {
    parserMap_t& parser = calculator::Default().parserMap;
    parser.add("True", &True);
    parser.add("False", &False);
    parser.add("None", &None);
    parser.add("#", &LineComment);
    parser.add("//", &LineComment);
    parser.add("/*", &SlashStarComment);
    parser.add(":", &KeywordOperator);
    parser.add(':', &KeywordOperator);
    parser.add(".", &DotOperator);
    parser.add('.', &DotOperator);
  }
} __CPARSE_STARTUP;

}  // namespace builtin_reservedWords
