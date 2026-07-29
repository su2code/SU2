/*!
 * \file gpu_ast.hpp
 * \brief definition of the Abstract Syntax Tree for evaluating CVecExpr expression trees on GPU with a single generic interpreter kernel
 * \author D. Di giusto
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once

#include <iostream>
#include "../parallelization/mpi_structure.hpp"

/*!
 * \brief enumeration of the possible node-operations on GPU to build abstract syntax trees
 */
enum class GPUOpType { VEC, SCALAR, ADD, SUB, MUL, DIV, NEG };

/*!
 * \brief GPU Abstract Syntax Tree node structure
 */
template <typename ScalarType>
struct GPUASTNode {
  GPUOpType op = GPUOpType::SCALAR;  // operation for this node
  int left = -1;                     // index of left child node-operation, -1 for leaf
  int right = -1;                    // index of right child node-operation, -1 for leaf
  const ScalarType* d_ptr = nullptr; // device pointer for vector
  ScalarType val = ScalarType(0);    // constant value for scalar
};

/*!
 * \brief arbitrary max limit on the number of nodes in one Abstract Syntax tree
 */
constexpr int MAX_AST_NODES = 16;

/*!
 * \brief struct for the abstract syntax tree payload sent from CPU to GPU during one kernel launch
 */
template <class ScalarType>
struct GPUASTPayload {
  GPUASTNode<ScalarType> nodes[MAX_AST_NODES]; // array hosting all nodes for a single expression-operator
  int root_idx = 0;
  int node_count = 0;

  // allocation method 
  int NewNode() {
    if (node_count >= MAX_AST_NODES) {
      SU2_MPI::Error("GPU AST node budget exceeded, you should increase MAX_AST_NODES!", __FUNCTION__);
    }
    return node_count++;
  }
};