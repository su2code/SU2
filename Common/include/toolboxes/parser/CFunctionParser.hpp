#pragma once
#include <string>
#include <vector>

namespace Parser {

class CFunctionParser {
 public:
  /*!
   * \brief Struct that represents a raw function stored as a string
   */
  struct RawFunction {
    std::string type;               //<! The type of the function
    std::string name;               //<! The name of the function
    std::string expr;               //<! The expression of the function
    std::vector<std::string> args;  //<! Arguments of the function
  };

  CFunctionParser() = default;

  /*!
   * \brief Parse functions from a file and store them
   * \param[in] fileName - The name of the file containing the functions
   */
  void ParseFunctionsFromFile(const std::string& fileName);

  /*!
   * \brief Parse functions from a string and store them
   * \param[in] functions - String containing the functions
   */
  void ParseFunctions(const std::string& functions);

  /*!
   * \brief Get all functions parsed using the ::ParseFunctionsFromFile call
   * \param[in] types - The types of functions to return
   * \return vector containing all functions of specified type
   */
  std::vector<RawFunction> GetFunctions(const std::vector<std::string>& types);

 private:
  /*!
   * \brief Vector containing all functions parsed using the ::ParseFunctionsFromFile call
   */
  std::vector<CFunctionParser::RawFunction> expressionsFromFile;

  /*!
   * \brief Prefix string that is appended to error messages
   */
  std::string errorPrefix;
};
}  // namespace Parser
