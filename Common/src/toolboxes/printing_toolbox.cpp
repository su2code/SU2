/*!
 * \file printing_toolbox.cpp
 * \brief Printing tools
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

#include <stdexcept>
#include <iomanip>
#include "../../include/toolboxes/printing_toolbox.hpp"

PrintingToolbox::CTablePrinter::CTablePrinter(std::ostream* output, const std::string& separator) {
  out_stream_ = output;
  i_ = 0;
  j_ = 0;
  separator_ = separator;
  table_width_ = 0;
  print_header_bottom_line_ = true;
  print_header_top_line_ = true;
  align_ = RIGHT;
  inner_separator_ = separator;
  precision_ = 6;
}

int PrintingToolbox::CTablePrinter::GetNumColumns() const { return (int)column_headers_.size(); }

int PrintingToolbox::CTablePrinter::GetTableWidth() const { return table_width_; }

void PrintingToolbox::CTablePrinter::SetSeparator(const std::string& separator) { separator_ = separator; }

void PrintingToolbox::CTablePrinter::SetInnerSeparator(const std::string& inner_separator) {
  inner_separator_ = inner_separator;
}

void PrintingToolbox::CTablePrinter::SetPrintHeaderBottomLine(bool print) { print_header_bottom_line_ = print; }

void PrintingToolbox::CTablePrinter::SetPrintHeaderTopLine(bool print) { print_header_top_line_ = print; }

void PrintingToolbox::CTablePrinter::SetAlign(int align) { align_ = align; }

void PrintingToolbox::CTablePrinter::SetPrecision(int precision) {
  precision_ = precision;
  out_stream_->precision(precision_);
}

void PrintingToolbox::CTablePrinter::AddColumn(const std::string& header_name, int column_width) {
  if (column_width < 4) {
    throw std::invalid_argument("Column size has to be >= 4");
  }

  column_headers_.push_back(header_name);
  column_widths_.push_back(column_width);
  table_width_ += column_width + separator_.size();  // for the separator
}

void PrintingToolbox::CTablePrinter::PrintHorizontalLine() {
  *out_stream_ << "+";  // the left bar

  for (int i = 0; i < table_width_ - 1; ++i) *out_stream_ << "-";

  *out_stream_ << "+";  // the right bar
  *out_stream_ << "\n";
}

void PrintingToolbox::CTablePrinter::PrintHeader() {
  if (print_header_top_line_) PrintHorizontalLine();
  *out_stream_ << separator_;
  int indent = 0;
  for (int i = 0; i < GetNumColumns(); ++i) {
    std::stringstream ss;

    ss << column_headers_.at(i).substr(0, column_widths_.at(i));

    indent = 0;

    if (align_ == LEFT)
      *out_stream_ << std::left;
    else if (align_ == RIGHT)
      *out_stream_ << std::right;
    else if (align_ == CENTER) {
      *out_stream_ << std::right;
      indent = (int)(column_widths_.at(i) - ss.str().size()) / 2;
    }

    *out_stream_ << std::setw(column_widths_.at(i) - indent) << column_headers_.at(i).substr(0, column_widths_.at(i));
    if (i != GetNumColumns() - 1) {
      *out_stream_ << std::setw(indent + (int)inner_separator_.size()) << inner_separator_;
    }
  }

  *out_stream_ << std::setw((int)separator_.size() + indent) << separator_;
  *out_stream_ << "\n";
  if (print_header_bottom_line_) {
    PrintHorizontalLine();
  }
}

void PrintingToolbox::CTablePrinter::PrintFooter() { PrintHorizontalLine(); }
