#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

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
class CTablePrinter{
public:
  CTablePrinter(std::ostream * output, const std::string & separator = "|");
  ~CTablePrinter();
  
  enum alignment {
    CENTER,
    LEFT,
    RIGHT
  };

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
  void SetSeparator(const std::string & separator);
  
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
  void AddColumn(const std::string & header_name, int column_width);
  
  /*!
   * \brief Print the header.
   */
  void PrintHeader();
  
  /*!
   * \brief Print the footer.
   */
  void PrintFooter();

  template<typename T> CTablePrinter& operator<<(T input){
    
    int indent = 0;
   
    /* --- Set the left separator --- */
    if (j_ == 0)
      *out_stream_ << "|";

    /* --- Determine and set the current alignment in the stream --- */
    if(align_ == LEFT)
      *out_stream_ << std::left;
    else if (align_ == RIGHT || align_ == CENTER)
      *out_stream_ << std::right; 

    /*--- Print the current column value to the stream --- */
    *out_stream_ << std::setw(column_widths_.at(j_) - indent)
                 << input;

    /*--- Reset the column counter and if it is the last column, 
     * add also a line break ---*/
    if (j_ == GetNumColumns()-1){
      *out_stream_ << std::setw(indent+2) << "|\n";
      i_ = i_ + 1;
      j_ = 0;
    } else {
      *out_stream_ << std::setw(indent+1) << separator_;
      j_ = j_ + 1;
    }

    return *this;
  }

private:
  
  /*!
   * \brief Print a horizontal line.
   */
  void PrintHorizontalLine();

  std::ostream * out_stream_;               /*< \brief The output stream. */
  std::vector<std::string> column_headers_; /*< \brief Vector of column header names. */
  std::vector<int> column_widths_;          /*< \brief Vector of column widths. */
  std::string separator_;                   /*< \brief Column separator char. */

  int i_; /*< \brief Index of the current row. */
  int j_; /*< \brief Index of the current column. */

  int table_width_;  /*< \brief The total width of the table. */
  int align_;        /*< \brief The current alignment. */
  bool print_header_top_line_,  /*< \brief Printing the header top line. */
  print_header_bottom_line_;   /*< \brief Printing the header bottom line. */
};

}
