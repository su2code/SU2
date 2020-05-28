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
  AD::Reset();
  AD::StartRecording();

  AD::RegisterInput(x);
  AD::RegisterInput(y);

  auto result = scope.EvalExpressions();

  AD::RegisterOutput(result[0]);
  AD::RegisterOutput(result[1]);
  AD::RegisterOutput(result[2]);
  AD::RegisterOutput(result[3]);

  AD::StopRecording();

  passivedouble r1_bar = 152;
  passivedouble r2_bar = 52.063;
  passivedouble r3_bar = 0.432;
  passivedouble r4_bar = 5.6788;

  SU2_TYPE::SetDerivative(result[0], r1_bar);
  SU2_TYPE::SetDerivative(result[1], r2_bar);
  SU2_TYPE::SetDerivative(result[2], r3_bar);
  SU2_TYPE::SetDerivative(result[3], r4_bar);

  AD::ComputeAdjoint();

  CHECK(SU2_TYPE::GetValue(result[0]) == x+y);
  CHECK(SU2_TYPE::GetValue(result[1]) == x*y);
  CHECK(SU2_TYPE::GetValue(result[2]) == x/y);
  CHECK(SU2_TYPE::GetValue(result[3]) == sqrt(x)/log(y));

  CHECK(SU2_TYPE::GetDerivative(x) == 0.5/(sqrt(x)*log(y))*r4_bar + r3_bar/y + r2_bar*y + r1_bar);
  CHECK(SU2_TYPE::GetDerivative(y) ==  - sqrt(x)*r4_bar/(y*pow(log(y),2)) - x*r3_bar/pow(y,2) + x*r2_bar + r1_bar);

}
