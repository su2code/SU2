/*!
 * \file CMatrixVectorProduct.hpp
 * \brief Headers for the classes related to sparse matrix-vector product wrappers.
 *        The actual operations are currently implemented mostly by CSysMatrix.
 * \author F. Palacios, J. Hicken, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../config_structure.hpp"
#include "../geometry_structure.hpp"
#include "CSysVector.hpp"
#include "CSysMatrix.hpp"

/*!
 * \class CMatrixVectorProduct
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
 */
template<class ScalarType>
class CMatrixVectorProduct {
public:
  virtual ~CMatrixVectorProduct() = 0; ///< class destructor
  virtual void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v)
  const = 0; ///< matrix-vector product operation
};
template<class ScalarType>
CMatrixVectorProduct<ScalarType>::~CMatrixVectorProduct() {}


/*!
 * \class CSysMatrixVectorProduct
 * \brief Specialization of matrix-vector product that uses CSysMatrix class
 */
template<class ScalarType>
class CSysMatrixVectorProduct : public CMatrixVectorProduct<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry;                   /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config;                       /*!< \brief pointer to matrix that defines the config. */

  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CSysMatrixVectorProduct();

public:

  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CSysMatrixVectorProduct(CSysMatrix<ScalarType> & matrix_ref,
                                 CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \brief destructor of the class
   */
  ~CSysMatrixVectorProduct() {}

  /*!
   * \brief operator that defines the CSysMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->MatrixVectorProduct(u, v, geometry, config);
  }
};


/*!
 * \class CSysMatrixVectorProductTransposed
 * \brief Specialization of matrix-vector product that uses CSysMatrix class for transposed products
 */
template<class ScalarType>
class CSysMatrixVectorProductTransposed : public CMatrixVectorProduct<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the product. */
  CGeometry* geometry;                   /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config;                       /*!< \brief pointer to matrix that defines the config. */

  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CSysMatrixVectorProductTransposed();

public:

  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the products
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CSysMatrixVectorProductTransposed(CSysMatrix<ScalarType> & matrix_ref,
                                           CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \brief destructor of the class
   */
  ~CSysMatrixVectorProductTransposed() {}

  /*!
   * \brief operator that defines the CSysMatrix-CSysVector product
   * \param[in] u - CSysVector that is being multiplied by the sparse matrix
   * \param[out] v - CSysVector that is the result of the product
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->MatrixVectorProductTransposed(u, v, geometry, config);
  }
};
