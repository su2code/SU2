#pragma once
#include "../../../Common/include/mpi_structure.hpp"
#include "../../../Common/include/option_structure.hpp"
#include <sys/stat.h>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>

#include "../../output/filewriter/CParallelDataSorter.hpp"

using namespace std;

class CFileWriter{
protected:

  std::vector<std::string> fieldnames;
  unsigned short nDim;
  
  int rank, size;
  
  std::string file_ext;
  
  su2double StartTime, StopTime, UsedTime, Bandwidth, file_size;
  
public:
  /*!
   * \brief CFileWriter
   * \param fields
   * \param nDim
   */  
  CFileWriter(std::vector<std::string> fields, unsigned short nDim);
  
  /*!
   * \brief ~CFileWriter
   */
  virtual ~CFileWriter();
  
  /*!
   * \brief Write_Data
   * \param filename
   * \param data_sorter
   */
  virtual void Write_Data(std::string filename, CParallelDataSorter* data_sorter){}
  
  /*!
   * \brief Get_Bandwidth
   */
  su2double Get_Bandwidth(){return Bandwidth;}
  
  /*!
   * \brief Get_Filesize
   */
  su2double Get_Filesize(){return file_size;}
  
  /*!
   * \brief Get_UsedTime
   * \return 
   */
  su2double Get_UsedTime(){return UsedTime;}
 
protected:
  /*!
   * \brief Determine_Filesize
   * \param filename
   * \return 
   */
  unsigned long Determine_Filesize(std::string filename){
      struct stat stat_buf;
      int rc = stat(filename.c_str(), &stat_buf);
      return rc == 0 ? stat_buf.st_size : -1;
  }
  
};

