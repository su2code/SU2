#pragma once

#include "../mpi_structure.hpp"
#include "../option_structure.hpp"

#define SU2_INFO_RANK(rank)\
  CSU2Logging(rank, loguru::Verbosity_INFO, __FILE__, __LINE__)

#define SU2_WARN_RANK(rank)\
  CSU2Logging(rank, loguru::Verbosity_WARNING, __FILE__, __LINE__)

#define SU2_INFO_ALL\
  CSU2Logging(-1, loguru::Verbosity_INFO, __FILE__, __LINE__)

#define SU2_WARN_ALL\
  CSU2Logging(-1, loguru::Verbosity_WARNING, __FILE__, __LINE__)

#define SU2_INFO \
  CSU2Logging(MASTER_NODE, loguru::Verbosity_INFO, __FILE__, __LINE__)
  
#define SU2_WARN \
  CSU2Logging(MASTER_NODE, loguru::Verbosity_WARNING, __FILE__, __LINE__)
  
class CSU2Logging {
public:
  
//  static void Info(string info, unsigned short rank_only = MASTER_NODE){
//    if (SU2_MPI::GetRank() == rank_only){
//     LOG_F(INFO, info.c_str());
//    }
//  }
  
  CSU2Logging(short rank, loguru::NamedVerbosity verb, const char* file, unsigned line) : 
      _file(file), _line(line), rank_only(rank), _verb(verb){}
  
  ~CSU2Logging(){
    if (SU2_MPI::GetRank() == rank_only || rank_only == -1)    
      loguru::StreamLogger(_verb, _file, _line) << _ss.str();
  }
  
  template<typename T>
  CSU2Logging& operator<<(const T& t){
    _ss << t;
    return *this;
  }
  
  CSU2Logging& operator<<(std::ostream&(*f)(std::ostream&))
  {
    f(_ss);
    return *this;
  }
  
private:
  const char* _file;
  unsigned _line;
  std::stringstream _ss;  
  short rank_only;
  loguru::NamedVerbosity _verb;
};

#define CHECK_WITH_INFO_S(cond, info)                                                              \
	LOGURU_PREDICT_TRUE((cond) == true)                                                            \
		? (void)0                                                                                  \
		: loguru::Voidify() & loguru::AbortLogger("CHECK FAILED:  " info "  ", __FILE__, __LINE__)

