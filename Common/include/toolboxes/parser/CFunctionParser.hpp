#pragma once
#include <string>
#include <vector>

namespace Parser{

class CFunctionParser {

public:

  struct RawExpression{
    std::string type;
    std::string name;
    std::string expr;
    std::vector<std::string> args;
  };

  CFunctionParser() = default;

  void ParseFunctionsFromFile(std::string fileName);
  std::vector<RawExpression> GetExpressions(std::vector<std::string> types);

private:
  std::vector<CFunctionParser::RawExpression> expressionsFromFile;

};


}