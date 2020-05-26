#include "../../../include/toolboxes/parser/CFunctionParser.hpp"
#include "../../../include/toolboxes/printing_toolbox.hpp"
#include "../../../include/mpi_structure.hpp"

#include <fstream>

namespace Parser{

void CFunctionParser::ParseFunctionsFromFile(std::string fileName) {

  std::ifstream functionFile;
  std::string UserFunctionCode;
  functionFile.open(fileName);
  if (functionFile.is_open()) {
    std::string line;
    while (!functionFile.eof()) {
      getline(functionFile, line);
      line = PrintingToolbox::trim(line);
      if (line.front() != '%')
        UserFunctionCode += line;
    }
  }

  auto size = UserFunctionCode.size();
  for (decltype(size) iChar = 0; iChar < size; iChar++) {
    if ((size - iChar >= 4) && UserFunctionCode.substr(iChar, 4) == "def ") {
      auto beginDef = iChar;
      auto endDef = UserFunctionCode.find_first_of("}", iChar);
      if (endDef == std::string::npos) {
        SU2_MPI::Error("Error while parsing function file. Missing \"}\".", CURRENT_FUNCTION);
      }

      auto defString = UserFunctionCode.substr(beginDef, endDef - beginDef + 1);

      auto posOpenBracket = defString.find_first_of("(");
      auto posCloseBracket = defString.find_first_of(")");

      if (posOpenBracket == std::string::npos) {
        SU2_MPI::Error("Error while parsing function file. Missing \"(\".", CURRENT_FUNCTION);
      }
      if (posCloseBracket == std::string::npos) {
        SU2_MPI::Error("Error while parsing function file. Missing \")\".", CURRENT_FUNCTION);
      }

      auto keyWords = PrintingToolbox::split(defString.substr(0, posOpenBracket), ' ');

      auto name = std::string();
      auto type = std::string();
      auto args =
          PrintingToolbox::split(defString.substr(posOpenBracket + 1, posCloseBracket - posOpenBracket - 1), ',');
      if (keyWords.size() == 2) {
        name = PrintingToolbox::trim(keyWords[1]);
        type = "global";
        for (auto arg : args) std::cout << arg << std::endl;
      } else if (keyWords.size() == 3) {
        name = PrintingToolbox::trim(keyWords[2]);
        type = PrintingToolbox::trim(keyWords[1]);
        if (args.size() != 0) {
          SU2_MPI::Error("No arguments allowed for definitions with specific type " + type, CURRENT_FUNCTION);
        }
      }

      auto posOpenExpBracket = defString.find_first_of("{", posCloseBracket + 1);
      if (posOpenExpBracket == std::string::npos) {
        SU2_MPI::Error("Error while parsing function file. Missing \"{\".", CURRENT_FUNCTION);
      }
      auto posCloseExpBracket = defString.find_first_of("}", posOpenExpBracket + 1);
      if (posCloseExpBracket == std::string::npos) {
        SU2_MPI::Error("Error while parsing function file. Missing \"}\".", CURRENT_FUNCTION);
      }
      auto expression = defString.substr(posOpenExpBracket + 1, posCloseExpBracket - posOpenExpBracket - 1);

      expressionsFromFile.push_back({type, name, expression, args});

      iChar = endDef;

    } else {
      SU2_MPI::Error("Error while parsing function file. Undefined keyword " +
                     PrintingToolbox::split(UserFunctionCode.substr(iChar), ' ')[0],
          CURRENT_FUNCTION);
    }

  }

  functionFile.close();

}

std::vector<CFunctionParser::RawExpression> CFunctionParser::GetExpressions(std::vector<std::string> types){

  std::vector<RawExpression> expressions;

  for (const auto &type : types){
    for (const auto& expr : expressionsFromFile){
      if (expr.type == type) expressions.push_back(expr);
    }
  }

  return expressions;
}
}