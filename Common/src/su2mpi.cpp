#include "../include/su2mpi.hpp"
/*
namespace SU2MPI {
  const int MASTER_NODE = 0;
  // Safetly exits with MPI
  void FinalizeAndExit1(){
#ifdef NO_MPI
    exit(1);
#else
#ifdef WINDOWS
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#else
    MPI::COMM_WORLD.Abort(1);
    MPI::Finalize();
#endif
#endif
  }
  
  !<\brief Prints to the head node and exits (using MPI if applicable)

  void PrintAndFinalize(std::string str){
    int rank = Rank();
    if (rank == MASTER_NODE){
      std::cout << str << std::endl;
    }
    FinalizeAndExit1();
  }
  
  !<\brief Returns the rank of the processor (always SU2MPI::MASTER_NODE if no MPI)

  int Rank(){
    int rank = MASTER_NODE;
#ifndef NO_MPI
#ifdef WINDOWS
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank = MPI::COMM_WORLD.Get_rank();
#endif
#endif
    return rank;
  }
}
*/