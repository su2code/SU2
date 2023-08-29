/*!
 * \file printing_toolbox.hpp
 * \brief Header file for the printing toolbox.
 * \author T. Albring
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "../basic_types/datatype_structure.hpp"

namespace PrintingToolbox {

/*!
 * \class CTablePrinter
 * \brief Class for writing output in a table.
 *
 * For example the output
 *
 * +-------------------------------------------+
 * |  MG Level| Presmooth|PostSmooth|CorrectSmo|
 * +-------------------------------------------+
 * |         0|         1|         0|         0|
 * |         1|         1|         0|         0|
 * |         2|         1|         0|         0|
 * |         3|         1|         0|         0|
 * +-------------------------------------------+
 *
 *
 * can be generated with the code
 *
 * CTablePrinter MGTable(&std::cout);
 * MGTable.AddColumn("MG Level",      10);
 * MGTable.AddColumn("Presmooth",     10);
 * MGTable.AddColumn("PostSmooth",    10);
 * MGTable.AddColumn("CorrectSmooth", 10);
 * MGTable.PrintHeader();
 * for (unsigned short iLevel = 0; iLevel < nMGLevels+1; iLevel++) {
 *   MGTable << iLevel << MG_PreSmooth[iLevel] << MG_PostSmooth[iLevel] << MG_CorrecSmooth[iLevel];
 * }
 * MGTable.PrintFooter();
 *
 *
 * \author T. Albring
 */
class CTablePrinter {
 public:
  CTablePrinter(std::ostream* output, const std::string& separator = "|");

  enum alignment { CENTER, LEFT, RIGHT };

  /*!
   * \brief Get number of columns of the table
   * \return Number of columns.
   */
  int GetNumColumns() const;

  /*!
   * \brief Get total width of the table.
   * \return Total width of the table.
   */
  int GetTableWidth() const;

  /*!
   * \brief Set the separator between columns (outer decoration)
   * \param[in] separator - The separation character.
   */
  void SetSeparator(const std::string& separator);

  /*!
   * \brief Set the separator between columns (inner decoration)
   * \param[in] separator - The separation character.
   */
  void SetInnerSeparator(const std::string& inner_separator);

  /*!
   * \brief Set the alignment of the table entries (CENTER only works for the header at the moment).
   * \param[in] align_ - The alignment (CENTER, LEFT, RIGHT).
   */
  void SetAlign(int align_);

  /*!
   * \brief Set whether to print the line at the bottom of the table.
   * \param[in] print - If TRUE, the bottom line is printed.
   */
  void SetPrintHeaderBottomLine(bool print);

  /*!
   * \brief Set whether to print the line at the top of the table.
   * \param[in] print - If TRUE, the top line is printed.
   */
  void SetPrintHeaderTopLine(bool print);

  /*!
   * \brief Add a column to the table by specifiying the header name and the width.
   * \param[in] header_name - The name printed in the header.
   * \param[in] column_width - The width of the column.
   */
  void AddColumn(const std::string& header_name, int column_width);

  /*!
   * \brief Print the header.
   */
  void PrintHeader();

  /*!
   * \brief Print the footer.
   */
  void PrintFooter();

  /*!
   * \brief Set the floating point precision.
   */
  void SetPrecision(int precision);

  template <typename T>
  CTablePrinter& operator<<(T input) {
    int indent = 0;

    /* --- Set the left separator --- */
    if (j_ == 0) *out_stream_ << separator_;

    /* --- Determine and set the current alignment in the stream --- */
    if (align_ == LEFT)
      *out_stream_ << std::left;
    else if (align_ == RIGHT || align_ == CENTER)
      *out_stream_ << std::right;

    /*--- Print the current column value to the stream --- */
    *out_stream_ << std::setw(column_widths_.at(j_) - indent) << std::setprecision(precision_) << input;

    /*--- Reset the column counter and if it is the last column,
     * add also a line break ---*/
    if (j_ == GetNumColumns() - 1) {
      *out_stream_ << std::setw(indent + 1 + (int)separator_.size()) << separator_ + "\n";
      i_ = i_ + 1;
      j_ = 0;
    } else {
      *out_stream_ << std::setw(indent + (int)inner_separator_.size()) << inner_separator_;
      j_ = j_ + 1;
    }

    return *this;
  }

 private:
  /*!
   * \brief Print a horizontal line.
   */
  void PrintHorizontalLine();

  std::ostream* out_stream_;                /*< \brief The output stream. */
  std::vector<std::string> column_headers_; /*< \brief Vector of column header names. */
  std::vector<int> column_widths_;          /*< \brief Vector of column widths. */
  std::string separator_;                   /*< \brief Column separator char. */
  std::string inner_separator_;             /*< \brief Inner column separator char. */

  int precision_; /*< \brief Floating point precision */

  int i_; /*< \brief Index of the current row. */
  int j_; /*< \brief Index of the current column. */

  int table_width_;              /*< \brief The total width of the table. */
  int align_;                    /*< \brief The current alignment. */
  bool print_header_top_line_,   /*< \brief Printing the header top line. */
      print_header_bottom_line_; /*< \brief Printing the header bottom line. */
};

inline void PrintScreenFixed(std::ostream& stream, su2double val, unsigned short field_width) {
  stream.precision(6);
  stream.setf(std::ios::fixed, std::ios::floatfield);
  stream.width(field_width);
  stream << std::right << val;
  stream.unsetf(std::ios::fixed);
}

inline void PrintScreenScientific(std::ostream& stream, su2double val, unsigned short field_width) {
  stream.precision(4);
  stream.setf(std::ios::scientific, std::ios::floatfield);
  stream.width(field_width);
  stream << std::right << val;
  stream.unsetf(std::ios::scientific);
}

inline void PrintScreenInteger(std::ostream& stream, unsigned long val, unsigned short field_width) {
  stream.width(field_width);
  stream << std::right << val;
}

inline void PrintScreenPercent(std::ostream& stream, su2double val, unsigned short field_width) {
  stream.precision(2);
  stream.setf(std::ios::fixed, std::ios::floatfield);
  stream.width(field_width - 1);
  stream << std::right << val << "%";
  stream.unsetf(std::ios::fixed);
}

inline std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

inline int stoi(const std::string s) {
  std::istringstream ss(s);
  int number;
  ss >> number;
  return number;
}

inline su2double stod(const std::string s) {
  std::istringstream ss(s);
  su2double number;
  ss >> number;
  return number;
}

inline std::string to_string(const su2double number) {
  std::stringstream ss;

  ss << number;

  return ss.str();
}

const static char* ws = " \t\n\r\f\v";

// trim from end of string (right)
inline std::string& rtrim(std::string& s, const char* t = ws) {
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = ws) {
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = ws) { return ltrim(rtrim(s, t), t); }

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in,out] str - string we want to convert
 */
inline void StringToUpperCase(std::string& str) { std::transform(str.begin(), str.end(), str.begin(), ::toupper); }

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in] str - string we want a copy of converted to uppercase
 * \return a copy of str in uppercase
 */
inline std::string StringToUpperCase(const std::string& str) {
  std::string upp_str(str);
  std::transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::toupper);
  return upp_str;
}

}  // namespace PrintingToolbox
