#include <stdexcept>
#include <iomanip>
#include <stdexcept>
#include "../include/toolboxes/printing_toolbox.hpp"

PrintingToolbox::CTablePrinter::CTablePrinter(std::ostream * output, const std::string & separator){
  out_stream_ = output;
  i_ = 0;
  j_ = 0;
  separator_ = separator;
  table_width_ = 0;
  flush_left_ = false;
}

PrintingToolbox::CTablePrinter::~CTablePrinter(){

}

int PrintingToolbox::CTablePrinter::get_num_columns() const {
  return column_headers_.size();
}

int PrintingToolbox::CTablePrinter::get_table_width() const {
  return table_width_;
}

void PrintingToolbox::CTablePrinter::set_separator(const std::string &separator){
  separator_ = separator;
}

void PrintingToolbox::CTablePrinter::set_flush_left(){
  flush_left_ = true;
}

void PrintingToolbox::CTablePrinter::set_flush_right(){
  flush_left_ = false;
}

/** \brief Add a column to our table
 ** 
 ** \param header_name Name to be print for the header
 ** \param column_width the width of the column (has to be >=5)
 ** */
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
  PrintHorizontalLine();
  *out_stream_ << "|";

  for (int i=0; i<get_num_columns(); ++i){

    if(flush_left_)
      *out_stream_ << std::left;
    else
      *out_stream_ << std::right; 

    *out_stream_ << std::setw(column_widths_.at(i)) << column_headers_.at(i).substr(0, column_widths_.at(i));
    if (i != get_num_columns()-1){
      *out_stream_ << separator_;
    }
  }

  *out_stream_ << "|\n";
  PrintHorizontalLine();
}

void PrintingToolbox::CTablePrinter::PrintFooter(){
  PrintHorizontalLine();
}

//TablePrinter& TablePrinter::operator<<(float input){
//  OutputDecimalNumber<float>(input);
//  return *this;
//}

//TablePrinter& TablePrinter::operator<<(double input){
//  OutputDecimalNumber<double>(input);
//  return *this;
//}

