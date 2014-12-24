/*!
 * \file su2mpi.cpp
 * \brief Header for caller functions of the turbulence models.
 * \author B. Tracey
 * \version 3.2.6 "eagle"
 *
 * Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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

#include "../include/su2mpi.hpp"
/*
namespace SU2MPI {
  const int MASTER_NODE = 0;
  // Safetly exits with MPI
  void FinalizeAndExit1(){
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
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
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
  }
}
*/