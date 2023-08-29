/*!
 * \file CPreconditioner.hpp
 * \brief Classes related to linear preconditioner wrappers.
 *        The actual operations are currently implemented mostly by CSysMatrix.
 * \author F. Palacios, J. Hicken, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

/// \addtogroup SpLinSys
/// @{

/*!
 * \class CPreconditioner
 * \brief Abstract base class for defining a preconditioning operation.
 * \author J. Hicken.
 *
 * See the remarks regarding the CMatrixVectorProduct class. The same
 * idea applies here to the preconditioning operation.
 */
template <class ScalarType>
class CPreconditioner {
 public:
  /*!
   * \brief Destructor of the class
   */
  virtual ~CPreconditioner() = 0;

  /*!
   * \brief Overload of operator (), applies the preconditioner to "u" storing the result in "v".
   */
  virtual void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const = 0;

  /*!
   * \brief Generic "preprocessing" hook derived classes may implement to build the preconditioner.
   */
  virtual void Build() {}

  /*!
   * \brief Return true to identify the identity preconditioner, may allow some solvers to be more efficient.
   */
  virtual bool IsIdentity() const { return false; }

  /*!
   * \brief Factory method.
   */
  static CPreconditioner* Create(ENUM_LINEAR_SOLVER_PREC kind, CSysMatrix<ScalarType>& jacobian, CGeometry* geometry,
                                 const CConfig* config);
};
template <class ScalarType>
CPreconditioner<ScalarType>::~CPreconditioner() {}

/*!
 * \class CJacobiPreconditioner
 * \brief Specialization of preconditioner that uses CSysMatrix class.
 */
template <class ScalarType>
class CJacobiPreconditioner final : public CPreconditioner<ScalarType> {
 private:
  CSysMatrix<ScalarType>& sparse_matrix; /*!< \brief Pointer to matrix that defines the preconditioner. */
  CGeometry* geometry;                   /*!< \brief Pointer to geometry associated with the matrix. */
  const CConfig* config;                 /*!< \brief Pointer to problem configuration. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] matrix_ref - Matrix reference that will be used to define the preconditioner.
   * \param[in] geometry_ref - Geometry associated with the problem.
   * \param[in] config_ref - Config of the problem.
   */
  inline CJacobiPreconditioner(CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref, const CConfig* config_ref)
      : sparse_matrix(matrix_ref) {
    if ((geometry_ref == nullptr) || (config_ref == nullptr))
      SU2_MPI::Error("Preconditioner needs to be built with valid references.", CURRENT_FUNCTION);
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid Pointers.
   */
  CJacobiPreconditioner() = delete;

  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
    sparse_matrix.ComputeJacobiPreconditioner(u, v, geometry, config);
  }

  /*!
   * \note Request the associated matrix to build the preconditioner.
   */
  inline void Build() override { sparse_matrix.BuildJacobiPreconditioner(); }
};

/*!
 * \class CILUPreconditioner
 * \brief Specialization of preconditioner that uses CSysMatrix class
 */
template <class ScalarType>
class CILUPreconditioner final : public CPreconditioner<ScalarType> {
 private:
  CSysMatrix<ScalarType>& sparse_matrix; /*!< \brief Pointer to matrix that defines the preconditioner. */
  CGeometry* geometry;                   /*!< \brief Pointer to geometry associated with the matrix. */
  const CConfig* config;                 /*!< \brief Pointer to problem configuration. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] matrix_ref - Matrix reference that will be used to define the preconditioner.
   * \param[in] geometry_ref - Geometry associated with the problem.
   * \param[in] config_ref - Config of the problem.
   */
  inline CILUPreconditioner(CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref, const CConfig* config_ref)
      : sparse_matrix(matrix_ref) {
    if ((geometry_ref == nullptr) || (config_ref == nullptr))
      SU2_MPI::Error("Preconditioner needs to be built with valid references.", CURRENT_FUNCTION);
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid Pointers.
   */
  CILUPreconditioner() = delete;

  /*!
   * \brief Operator that defines the preconditioner operation.
   * \param[in] u - CSysVector that is being preconditioned.
   * \param[out] v - CSysVector that is the result of the preconditioning.
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
    sparse_matrix.ComputeILUPreconditioner(u, v, geometry, config);
  }

  /*!
   * \note Request the associated matrix to build the preconditioner.
   */
  inline void Build() override { sparse_matrix.BuildILUPreconditioner(); }
};

/*!
 * \class CLU_SGSPreconditioner
 * \brief Specialization of preconditioner that uses CSysMatrix class.
 */
template <class ScalarType>
class CLU_SGSPreconditioner final : public CPreconditioner<ScalarType> {
 private:
  CSysMatrix<ScalarType>& sparse_matrix; /*!< \brief Pointer to matrix that defines the preconditioner. */
  CGeometry* geometry;                   /*!< \brief Pointer to geometry associated with the matrix. */
  const CConfig* config;                 /*!< \brief Pointer to problem configuration. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] matrix_ref - Matrix reference that will be used to define the preconditioner.
   * \param[in] geometry_ref - Geometry associated with the problem.
   * \param[in] config_ref - Config of the problem.
   */
  inline CLU_SGSPreconditioner(CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref, const CConfig* config_ref)
      : sparse_matrix(matrix_ref) {
    if ((geometry_ref == nullptr) || (config_ref == nullptr))
      SU2_MPI::Error("Preconditioner needs to be built with valid references.", CURRENT_FUNCTION);
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid Pointers.
   */
  CLU_SGSPreconditioner() = delete;

  /*!
   * \brief operator that defines the preconditioner operation.
   * \param[in] u - CSysVector that is being preconditioned.
   * \param[out] v - CSysVector that is the result of the preconditioning.
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
    sparse_matrix.ComputeLU_SGSPreconditioner(u, v, geometry, config);
  }
};

/*!
 * \class CLineletPreconditioner
 * \brief Specialization of preconditioner that uses CSysMatrix class.
 */
template <class ScalarType>
class CLineletPreconditioner final : public CPreconditioner<ScalarType> {
 private:
  CSysMatrix<ScalarType>& sparse_matrix; /*!< \brief Pointer to matrix that defines the preconditioner. */
  CGeometry* geometry;                   /*!< \brief Pointer to geometry associated with the matrix. */
  const CConfig* config;                 /*!< \brief Pointer to problem configuration. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] matrix_ref - Matrix reference that will be used to define the preconditioner.
   * \param[in] geometry_ref - Geometry associated with the problem.
   * \param[in] config_ref - Config of the problem.
   */
  inline CLineletPreconditioner(CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref, const CConfig* config_ref)
      : sparse_matrix(matrix_ref) {
    if ((geometry_ref == nullptr) || (config_ref == nullptr))
      SU2_MPI::Error("Preconditioner needs to be built with valid references.", CURRENT_FUNCTION);
    geometry = geometry_ref;
    config = config_ref;
  }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid Pointers.
   */
  CLineletPreconditioner() = delete;

  /*!
   * \brief Operator that defines the preconditioner operation.
   * \param[in] u - CSysVector that is being preconditioned.
   * \param[out] v - CSysVector that is the result of the preconditioning.
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
    sparse_matrix.ComputeLineletPreconditioner(u, v, geometry, config);
  }

  /*!
   * \note Request the associated matrix to build the preconditioner.
   */
  inline void Build() override { sparse_matrix.BuildLineletPreconditioner(geometry, config); }
};

/*!
 * \class CPastixPreconditioner
 * \brief Specialization of preconditioner that uses PaStiX to factorize a CSysMatrix.
 */
template <class ScalarType>
class CPastixPreconditioner final : public CPreconditioner<ScalarType> {
 private:
  CSysMatrix<ScalarType>& sparse_matrix; /*!< \brief Pointer to the matrix. */
  CGeometry* geometry;                   /*!< \brief Geometry associated with the problem. */
  const CConfig* config;                 /*!< \brief Configuration of the problem. */
  unsigned short kind_fact;              /*!< \brief The type of factorization desired. */

 public:
  /*!
   * \brief Constructor of the class
   * \param[in] matrix_ref - Matrix reference that will be used to define the preconditioner.
   * \param[in] geometry_ref - Associated geometry.
   * \param[in] config_ref - Problem configuration.
   * \param[in] kind_factorization - Type of factorization required.
   */
  inline CPastixPreconditioner(CSysMatrix<ScalarType>& matrix_ref, CGeometry* geometry_ref, const CConfig* config_ref,
                               unsigned short kind_factorization)
      : sparse_matrix(matrix_ref) {
    if ((geometry_ref == nullptr) || (config_ref == nullptr))
      SU2_MPI::Error("Preconditioner needs to be built with valid references.", CURRENT_FUNCTION);
    geometry = geometry_ref;
    config = config_ref;
    kind_fact = kind_factorization;
  }

  /*!
   * \note This class cannot be default constructed as that would leave us with invalid Pointers.
   */
  CPastixPreconditioner() = delete;

  /*!
   * \brief Operator that defines the preconditioner operation.
   * \param[in] u - CSysVector that is being preconditioned.
   * \param[out] v - CSysVector that is the result of the preconditioning.
   */
  inline void operator()(const CSysVector<ScalarType>& u, CSysVector<ScalarType>& v) const override {
    sparse_matrix.ComputePastixPreconditioner(u, v, geometry, config);
  }

  /*!
   * \note Request the associated matrix to build the preconditioner.
   */
  inline void Build() override { sparse_matrix.BuildPastixPreconditioner(geometry, config, kind_fact); }
};

template <class ScalarType>
CPreconditioner<ScalarType>* CPreconditioner<ScalarType>::Create(ENUM_LINEAR_SOLVER_PREC kind,
                                                                 CSysMatrix<ScalarType>& jacobian, CGeometry* geometry,
                                                                 const CConfig* config) {
  CPreconditioner<ScalarType>* prec = nullptr;

  switch (kind) {
    case JACOBI:
      prec = new CJacobiPreconditioner<ScalarType>(jacobian, geometry, config);
      break;
    case LINELET:
      prec = new CLineletPreconditioner<ScalarType>(jacobian, geometry, config);
      break;
    case LU_SGS:
      prec = new CLU_SGSPreconditioner<ScalarType>(jacobian, geometry, config);
      break;
    case ILU:
      prec = new CILUPreconditioner<ScalarType>(jacobian, geometry, config);
      break;
    case PASTIX_ILU:
    case PASTIX_LU_P:
    case PASTIX_LDLT_P:
      prec = new CPastixPreconditioner<ScalarType>(jacobian, geometry, config, kind);
      break;
  }

  return prec;
}

/// @}
