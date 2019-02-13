#include <stdexcept>
#include <iomanip>
#include <stdexcept>
#include "../../include/toolboxes/printing_toolbox.hpp"

PrintingToolbox::CTablePrinter::CTablePrinter(std::ostream * output, const std::string & separator){
  out_stream_ = output;
  i_ = 0;
  j_ = 0;
  separator_ = separator;
  table_width_ = 0;
  print_header_bottom_line_ = true;
  print_header_top_line_    = true;
  align_ = RIGHT;
}

PrintingToolbox::CTablePrinter::~CTablePrinter(){

}

int PrintingToolbox::CTablePrinter::GetNumColumns() const {
  return (int)column_headers_.size();
}

int PrintingToolbox::CTablePrinter::GetTableWidth() const {
  return table_width_;
}

void PrintingToolbox::CTablePrinter::SetSeparator(const std::string &separator){
  separator_ = separator;
}

void PrintingToolbox::CTablePrinter::SetPrintHeaderBottomLine(bool print){
  print_header_bottom_line_ = print;
}

void PrintingToolbox::CTablePrinter::SetPrintHeaderTopLine(bool print){
  print_header_top_line_ = print;
}

void PrintingToolbox::CTablePrinter::SetAlign(int align){
  align_ = align;
}

void PrintingToolbox::CTablePrinter::AddColumn(const std::string & header_name, int column_width){
  if (column_width < 4){
    throw std::invalid_argument("Column size has to be >= 4");
  }

  column_headers_.push_back(header_name);
  column_widths_.push_back(column_width);
  table_width_ += column_width + separator_.size(); // for the separator  
}

void PrintingToolbox::CTablePrinter::PrintHorizontalLine() {
  *out_stream_ << "+"; // the left bar

  for (int i=0; i<table_width_-1; ++i)
    *out_stream_ << "-";

  *out_stream_ << "+"; // the right bar
  *out_stream_ << "\n";
}

void PrintingToolbox::CTablePrinter::PrintHeader(){

  
  
  if (print_header_top_line_) PrintHorizontalLine();
  *out_stream_ << "|";
  int indent = 0;
  for (int i=0; i<GetNumColumns(); ++i){
    
    std::stringstream ss;
    
    ss << column_headers_.at(i).substr(0, column_widths_.at(i));
    
    indent = 0;

    if(align_ == LEFT)
      *out_stream_ << std::left;
    else if (align_ == RIGHT)
      *out_stream_ << std::right; 
    else if (align_ == CENTER) {
      *out_stream_ << std::right; 
      indent = (int)(column_widths_.at(i) - ss.str().size()) / 2;
    }

    *out_stream_ << std::setw(column_widths_.at(i) - indent) << column_headers_.at(i).substr(0, column_widths_.at(i));
    if (i != GetNumColumns()-1){
      *out_stream_ << std::setw(1+indent) << separator_;
    }
  }

  *out_stream_ << std::setw(2+indent) << "|\n";
  if (print_header_bottom_line_) PrintHorizontalLine();
}

void PrintingToolbox::CTablePrinter::PrintFooter(){
  PrintHorizontalLine();
}

