#include "catch.hpp"
#include "../../../../Common/include/toolboxes/parser/ExpressionParser.hpp"

TEST_CASE("Expression parser and scope", "[Parser]"){

  Parser::Scope scope;

  su2double x, y;

  scope.addVariable("x", x);
  scope.addVariable("y", y);
  REQUIRE(scope.addExpression("x + y;"));
  REQUIRE(scope.addExpression("x * y;"));
  REQUIRE(scope.addExpression("x/y"));
  REQUIRE(scope.addExpression("sqrt(x)/log(y)"));

  x = 5.5;
  y = 10.2;

  auto result = scope.EvalExpressions();


  REQUIRE(result[0] == x+y);
  REQUIRE(result[1] == x*y);
  REQUIRE(result[2] == x/y);
  REQUIRE(result[3] == sqrt(x)/log(y));

}


TEST_CASE("Custom function sum", "[Parser]"){

  Parser::Scope scope;

  std::string s1 = "surface1";
  std::string s2 = "surface2";

  su2double x_surf1, x_surf2, y_surf1, y_surf2;

  scope.addVariable("x@surface1", x_surf1);
  scope.addVariable("y@surface1", y_surf1);
  scope.addVariable("x@surface2", x_surf2);
  scope.addVariable("y@surface2", y_surf2);

  scope.addStringVar("surface1", s1);
  scope.addStringVar("surface2", s2);

  scope.addCustomFunction("x@", Parser::CustomFunction::SURFACE_SUM);
  scope.addCustomFunction("y@", Parser::CustomFunction::SURFACE_SUM);

  REQUIRE(scope.addExpression("x@(surface1, surface2);"));
  REQUIRE(scope.addExpression("y@(surface1, surface2);"));

  x_surf1 = 10.4;
  y_surf1 = 45.23;
  x_surf2 = 0.35;
  y_surf2 = 421.042;

  auto result = scope.EvalExpressions();

  REQUIRE(result[0] == x_surf1 + x_surf2);
  REQUIRE(result[1] == y_surf1 + y_surf2);
}
