#pragma once
#include "../../../Common/include/geometry_structure.hpp"
#include "CParallelDataSorter.hpp"
#include <sys/stat.h>

class CFileWriter{
protected:

  vector<string> fieldnames;
  unsigned short nDim;
  
  int rank, size;
  
  string file_ext;
  
  su2double StartTime, StopTime, UsedTime, Bandwidth, file_size;
  
public:
  /*!
   * \brief CFileWriter
   * \param fields
   * \param nDim
   */  
  CFileWriter(vector<string> fields, unsigned short nDim);
  
  /*!
   * \brief ~CFileWriter
   */
  virtual ~CFileWriter();
  
  /*!
   * \brief Write_Data
   * \param filename
   * \param data_sorter
   */
  virtual void Write_Data(string filename, CParallelDataSorter* data_sorter){}
  
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

