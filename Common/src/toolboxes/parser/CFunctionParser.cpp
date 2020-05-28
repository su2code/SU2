#include "../../../include/toolboxes/parser/CFunctionParser.hpp"

#include <fstream>

#include "../../../include/mpi_structure.hpp"
#include "../../../include/toolboxes/printing_toolbox.hpp"

namespace Parser {

  void CFunctionParser::ParseFunctionsFromFile(const std::string &fileName){
    std::ifstream functionFile;
    std::string UserFunctionCode;
    functionFile.open(fileName);
    if (functionFile.is_open()) {
      std::string line;
      while (!functionFile.eof()) {
        getline(functionFile, line);
        line = PrintingToolbox::trim(line);
        if (!line.empty()) line = PrintingToolbox::split(line, '%')[0];
        if (line.front() != '%') UserFunctionCode += line;
      }
    }
    functionFile.close();

    errorPrefix = "While parsing functions from file " + fileName;
    ParseFunctions(UserFunctionCode);
    errorPrefix.clear();

  }

void CFunctionParser::ParseFunctions(const std::string& functions) {

  expressionsFromFile.clear();
  if (errorPrefix.empty()){
    errorPrefix = "While parsing functions ";
  }

  auto size = functions.size();
  for (decltype(size) iChar = 0; iChar < size; iChar++) {
    if (functions[iChar] == ' ')
      continue;
    if ((size - iChar >= 4) && functions.substr(iChar, 4) == "def ") {
      auto beginDef = iChar;
      auto endDef = functions.find("end", iChar);
      if (endDef == std::string::npos) {
        SU2_MPI::Error(errorPrefix + ": Missing \"end\" keyword.", CURRENT_FUNCTION);
      }
      endDef += 3;

      auto defString = functions.substr(beginDef, endDef - beginDef);

      auto posOpenBracket = defString.find_first_of("(");
      auto posCloseBracket = defString.find_first_of(")");

      if (posOpenBracket == std::string::npos) {
        SU2_MPI::Error(errorPrefix + ": Missing \"(\".", CURRENT_FUNCTION);
      }
      if (posCloseBracket == std::string::npos) {
        SU2_MPI::Error(errorPrefix + ": Missing \")\".", CURRENT_FUNCTION);
      }

      auto keyWords = PrintingToolbox::split(defString.substr(0, posOpenBracket), ' ');

      auto name = std::string();
      auto type = std::string();
      auto args =
          PrintingToolbox::split(defString.substr(posOpenBracket + 1, posCloseBracket - posOpenBracket - 1), ',');
      for (auto& arg : args) {
        arg = PrintingToolbox::trim(arg);
        if (arg.find_first_of(" ") != std::string::npos){
          SU2_MPI::Error(errorPrefix + ": \",\" missing between arguments",
                         CURRENT_FUNCTION);
        }
      }
      if (keyWords.size() == 2) {
        name = PrintingToolbox::trim(keyWords[1]);
        type = "global";
      } else if (keyWords.size() == 3) {
        name = PrintingToolbox::trim(keyWords[2]);
        type = PrintingToolbox::trim(keyWords[1]);
        if (type != "global" && args.size() != 0) {
          SU2_MPI::Error(errorPrefix + ": No arguments allowed for definitions with specific type " + type,
                         CURRENT_FUNCTION);
        }
      }

      auto posOpenExpBracket = defString.find_first_of(":", posCloseBracket + 1);
      if (posOpenExpBracket == std::string::npos) {
        SU2_MPI::Error(errorPrefix + ": Missing \":\".", CURRENT_FUNCTION);
      }
      auto posCloseExpBracket = defString.find("end", posOpenExpBracket + 1);
      if (posCloseExpBracket == std::string::npos) {
        SU2_MPI::Error(errorPrefix + ": Missing \"end\" keyword", CURRENT_FUNCTION);
      }
      auto expression = defString.substr(posOpenExpBracket + 1, posCloseExpBracket - posOpenExpBracket - 1);
      expression = PrintingToolbox::trim(expression);
      expressionsFromFile.push_back({type, name, expression, args});

      iChar = endDef - 1;

    } else {
      SU2_MPI::Error(errorPrefix + ": Undefined keyword \"" +
                         PrintingToolbox::split(functions.substr(iChar), ' ')[0] + "\"",
                     CURRENT_FUNCTION);
    }
  }
}

std::vector<CFunctionParser::RawFunction> CFunctionParser::GetFunctions(const std::vector<std::string>& types) {
  std::vector<RawFunction> expressions;

  for (const auto& type : types) {
    for (const auto& expr : expressionsFromFile) {
      if (expr.type == type) expressions.push_back(expr);
    }
  }

  return expressions;
}
}  // namespace Parser
