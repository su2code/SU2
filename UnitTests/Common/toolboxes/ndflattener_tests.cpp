/*!
 * \file ndflattener_tests.cpp
 * \brief Unit tests for NdFlattener template classes.
 * \author M. Aehle
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "catch.hpp"
#include "../../Common/include/toolboxes/ndflattener.hpp"

TEST_CASE("NdFlattener Test", "[NdFlattener]"){

  int rank; SU2_MPI::Comm_rank(SU2_MPI::GetComm(), &rank);
  int size; SU2_MPI::Comm_size(SU2_MPI::GetComm(), &size);

  /*-- Provide non-flat array --*/
  su2double** A = new su2double*[2];
  A[0] = new su2double[2]; A[0][0] = 0.0; A[0][1] = 1.0;
  A[1] = new su2double[3+rank];
  for(int i=0; i<3+rank; i++)
    A[1][i] = 2.0 + rank + i;

  /*-- Accessor --*/
  auto f = std::make_pair( (size_t)2,  [rank,A](int i) {
    return std::make_pair( (size_t)(i==0?2:(3+rank)), [rank,A,i](int j){
      return  A[i][j];
    });
  });

  /*-- Read into flattening structure --*/
  NdFlattener<2> nd2(f);

  /*-- Modify A -> this should not alter nd2 at this point --*/
  A[0][0] = 0.5;

  /*-- Check structure --*/
  REQUIRE( nd2.getNChildren() == 2 );
  REQUIRE( nd2.get( 0, 0) == 0.0 );
  REQUIRE( nd2[0][1] == 1.0 );
  REQUIRE( nd2.getNChildren(1) == 3 + rank );
  for(int i=0; i<3+rank; i++){
    REQUIRE( nd2[1][i] == 2.0 + rank + i );
  }

  /*-- gather flattening structures of all processes --*/
  NdFlattener<3> nd3(Get_Nd_MPI_Env(), &nd2);

  /*-- Check gathered structure --*/
  for(int r=0; r<size; r++){
    REQUIRE( nd3.getNChildren(r) == 2 );
    REQUIRE( nd3.get(r, 0, 0) == 0.0 );
    REQUIRE( nd3[r][0][1] == 1.0 );
    REQUIRE( nd3.getNChildren(r,1) == 3 + rank );
    for(int i=0; i<3+rank; i++){
      REQUIRE( nd3[r][1][i] == 2.0 + rank + i );
    }
  }
    
  /*-- reread modified A and check again --*/
  nd2.initialize_or_refresh(f);
  nd3.refresh(Get_Nd_MPI_Env(), &nd2);
  for(int r=0; r<size; r++){
    REQUIRE( nd3.getNChildren(r) == 2 );
    REQUIRE( nd3.get(r, 0, 0) == 0.5 );
    REQUIRE( nd3[r][0][1] == 1.0 );
    REQUIRE( nd3.getNChildren(r,1) == 3 + rank );
    for(int i=0; i<3+rank; i++){
      REQUIRE( nd3[r][1][i] == 2.0 + rank + i );
    }
  }

  /*-- free stuff --*/
  delete[] A[0]; delete[] A[1]; delete[] A;
}
