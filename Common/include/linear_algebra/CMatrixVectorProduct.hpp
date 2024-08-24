/*!
 * \file CMatrixVectorProduct.hpp
 * \brief Headers for the classes related to sparse matrix-vector product wrappers.
 *        The actual operations are currently implemented mostly by CSysMatrix.
 * \author F. Palacios, J. Hicken, T. Economon
 * \version 8.0.1 "Harrier"
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

#include "../CConfig.hpp"
#include "../geometry/CGeometry.hpp"
#include "CSysVector.hpp"
#include "CSysMatrix.hpp"

/*!
 * \class CMatrixVectorProduct
 * \ingroup SpLinSys
 * \brief Abstract base class for defining matrix-vector products
 * \author J. Hicken.
 *
 * The Krylov-subspace solvers require only matrix-vector products and
 * not the actual matrix/Jacobian.  We need some way to indicate which
 * function will perform the product.  However, sometimes the
 * functions that define the product will require different numbers
 * and types of inputs.  For example, the forward-difference
 * approximation to a Jacobian-vector product requires the vector that
 * defines the Jacobian and a perturbation parameter.  The
 * CMatrixVectorProduct class is used to derive child classes that can
 * handle the different types of matrix-vector products and still be
 * passed to a single implementation of the Krylov solvers.
 * This abstraction may also be used to define matrix-free products.
 * 
 * There is also the use of a dummy class being made to select the 
 * correct function as defined by the user while deciding between
 * CPU or GPU execution. This dummy class calls the correct member
 * functions from its derived classes to map the suitable path of 
 * execution - CPU or GPU.
 */
template <class ScalarType>
class CMatrixVectorProduct {
 public:
  virtual ~CMatrixVectorProduct() = 0;
  virtual void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const = 0;
};
template <class ScalarType>
CMatrixVectorProduct<ScalarType>::~CMatrixVectorProduct() {}

/*!
 * \class executionPath
 * \brief Dummy super class that holds the correct member functions in its child classes
 */

template <class ScalarType>
class executionPath{
public:
   virtual void mat_vec_prod(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v, CGeometry* geometry, const CConfig* config, const CSysMatrix<ScalarType>& matrix) = 0;
};

/*!
 * \class cpuExecution
 * \brief Derived class containing the CPU Matrix Vector Product Function
 */
template <class ScalarType>
class cpuExecution : public executionPath<ScalarType> {
public:
    void mat_vec_prod(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v, CGeometry* geometry, const CConfig* config, const CSysMatrix<ScalarType>& matrix) override
    {
      matrix.MatrixVectorProduct(u, v, geometry, config);
    }
};

/*!
 * \class gpuExecution
 * \brief Derived class containing the GPU Matrix Vector Product Function
 */
template <class ScalarType> 
class gpuExecution : public executionPath<ScalarType> {
public:
     void mat_vec_prod(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v, CGeometry* geometry, const CConfig* config, const CSysMatrix<ScalarType>& matrix) override
    {
      #ifdef HAVE_CUDA
      matrix.GPUMatrixVectorProduct(u, v, geometry, config);
      #else
      matrix.MatrixVectorProduct(u, v, geometry, config);
      #endif
    }
};

/*!
 * \class CSysMatrixVectorProduct
 * \ingroup SpLinSys
 * \brief Specialization of matrix-vector product that uses CSysMatrix class
 */
template <class ScalarType>
class CSysMatrixVectorProduct final : public CMatrixVectorProduct<ScalarType> {
 private:
  const CSysMatrix<ScalarType>& matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry;                  /*!< \brief geometry associated with the matrix. */
  const CConfig* config;                /*!< \brief config of the problem. */
  executionPath<ScalarType>* exec;      /*!< \brief interface that decides which path of execution to choose from. */

 public:
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CSysMatrixVectorProduct(const CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref,
                                 const CConfig* config_ref)
      : matrix(matrix_ref), geometry(geometry_ref), config(config_ref) 
      {
         if(config->GetCUDA())
         {
            exec = new gpuExecution<ScalarType>;
         }
         else
         {
            exec = new cpuExecution<ScalarType>;
         }
      }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CSysMatrixVectorProduct() = delete;

  /*!
   * \brief operator that defines the CSysMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
  
   exec->mat_vec_prod(u, v, geometry, config, matrix);

  }
};

