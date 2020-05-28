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
  passivedouble x_dot = 53.2;
  passivedouble y_dot = 21.5;
  SU2_TYPE::SetDerivative(x, x_dot);
  SU2_TYPE::SetDerivative(y, y_dot);

  const auto result = scope.EvalExpressions();

  REQUIRE(SU2_TYPE::GetValue(result[0]) == x+y);
  REQUIRE(SU2_TYPE::GetDerivative(result[0]) == x_dot + y_dot );
  REQUIRE(SU2_TYPE::GetValue(result[1]) == x*y);
  REQUIRE(SU2_TYPE::GetDerivative(result[1]) == x_dot*y + x*y_dot);
  REQUIRE(SU2_TYPE::GetValue(result[2]) == x/y);
  REQUIRE(SU2_TYPE::GetDerivative(result[2]) == x_dot/y - x*y_dot/pow(y,2));
  REQUIRE(SU2_TYPE::GetValue(result[3]) == sqrt(x)/log(y));
  REQUIRE(SU2_TYPE::GetDerivative(result[3]) == 0.5/(sqrt(x)*log(y))*x_dot - sqrt(x)*y_dot/(y*pow(log(y),2)));

}
