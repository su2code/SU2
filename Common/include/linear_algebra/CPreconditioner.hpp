/*!
 * \file CPreconditioner.hpp
 * \brief Headers for the classes related to linear preconditioner wrappers.
 *        The actual operations are currently implemented mostly by CSysMatrix.
 * \author F. Palacios, J. Hicken, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include "../vector_structure.hpp"
#include "../matrix_structure.hpp"
#include "../geometry_structure.hpp"

/*!
 * \class CPreconditioner
 * \brief abstract base class for defining preconditioning operation
 * \author J. Hicken.
 *
 * See the remarks regarding the CMatrixVectorProduct class. The same
 * idea applies here to the preconditioning operation.
 */
template<class ScalarType>
class CPreconditioner {
public:
  virtual ~CPreconditioner() = 0; ///< class destructor
  virtual void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v)
  const = 0; ///< preconditioning operation
};
template<class ScalarType>
CPreconditioner<ScalarType>::~CPreconditioner() {}


/*!
 * \class CJacobiPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CJacobiPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CJacobiPreconditioner();
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CJacobiPreconditioner(CSysMatrix<ScalarType> & matrix_ref,
                               CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }
  
  /*!
   * \brief destructor of the class
   */
  ~CJacobiPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->ComputeJacobiPreconditioner(u, v, geometry, config);
  }
};


/*!
 * \class CILUPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CILUPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CILUPreconditioner();
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CILUPreconditioner(CSysMatrix<ScalarType> & matrix_ref,
                            CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }
  
  /*!
   * \brief destructor of the class
   */
  ~CILUPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->ComputeILUPreconditioner(u, v, geometry, config);
  }
};


/*!
 * \class CLU_SGSPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CLU_SGSPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CLU_SGSPreconditioner();
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CLU_SGSPreconditioner(CSysMatrix<ScalarType> & matrix_ref,
                               CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }
  
  /*!
   * \brief destructor of the class
   */
  ~CLU_SGSPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->ComputeLU_SGSPreconditioner(u, v, geometry, config);
  }
};


/*!
 * \class CLineletPreconditioner
 * \brief specialization of preconditioner that uses CSysMatrix class
 */
template<class ScalarType>
class CLineletPreconditioner : public CPreconditioner<ScalarType> {
private:
  CSysMatrix<ScalarType>* sparse_matrix; /*!< \brief pointer to matrix that defines the preconditioner. */
  CGeometry* geometry; /*!< \brief pointer to matrix that defines the geometry. */
  CConfig* config; /*!< \brief pointer to matrix that defines the config. */
  
  /*!
   * \brief Default constructor of the class
   * \note This class cannot be default constructed as that would leave us with invalid pointers.
   */
  CLineletPreconditioner();
  
public:
  
  /*!
   * \brief constructor of the class
   * \param[in] matrix_ref - matrix reference that will be used to define the preconditioner
   * \param[in] geometry_ref - geometry associated with the problem
   * \param[in] config_ref - config of the problem
   */
  inline CLineletPreconditioner(CSysMatrix<ScalarType> & matrix_ref,
                                CGeometry *geometry_ref, CConfig *config_ref) {
    sparse_matrix = &matrix_ref;
    geometry = geometry_ref;
    config = config_ref;
  }
  
  /*!
   * \brief destructor of the class
   */
  ~CLineletPreconditioner() {}
  
  /*!
   * \brief operator that defines the preconditioner operation
   * \param[in] u - CSysVector that is being preconditioned
   * \param[out] v - CSysVector that is the result of the preconditioning
   */
  inline void operator()(const CSysVector<ScalarType> & u, CSysVector<ScalarType> & v) const {
    sparse_matrix->ComputeLineletPreconditioner(u, v, geometry, config);
  }
};
