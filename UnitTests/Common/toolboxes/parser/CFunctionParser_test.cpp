#include "catch.hpp"
#include "../../../../Common/include/toolboxes/parser/CFunctionParser.hpp"

TEST_CASE("Function parser", "[Parser]"){

  const std::string functions = "def global test1( a ,b ):"
                                "  a + b;"
                                "end"
                                "def global test2( arg1 ,arg2 ):"
                                "  var x:= 0;"
                                "  if (arg1 > arg2) {"
                                "     x+= arg1;"
                                "  }"
                                "end";

  Parser::CFunctionParser funcParser;

  funcParser.ParseFunctions(functions);

  const auto rawFunctions = funcParser.GetFunctions({"global"});

  REQUIRE(rawFunctions[0].name == "test1");
  REQUIRE(rawFunctions[0].expr == "a + b;");
  REQUIRE(rawFunctions[0].args[0] == "a");
  REQUIRE(rawFunctions[0].args[1] == "b");
  REQUIRE(rawFunctions[0].type == "global");
  REQUIRE(rawFunctions[1].name == "test2");
  REQUIRE(rawFunctions[1].expr == "var x:= 0;  if (arg1 > arg2) {     x+= arg1;  }");
  REQUIRE(rawFunctions[1].args[0] == "arg1");
  REQUIRE(rawFunctions[1].args[1] == "arg2");
  REQUIRE(rawFunctions[1].type == "global");

}
