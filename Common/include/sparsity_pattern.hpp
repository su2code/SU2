#pragma once

#include "geometry_structure.hpp"

/*!
 * \class CSparsityPattern
 * \brief Base class for defining a sparsity pattern for the CRS matrix structure. 
 * \ingroup crs_matrix
 * \author F. Palacios, T. Albring
 */
class CSparsityPattern{
  
private:
  unsigned long nPoint;
  unsigned long nnz;                 /*!< \brief Number of possible nonzero entries in the matrix. */
  unsigned long *row_ptr;            /*!< \brief Pointers to the first element in each row. */
  unsigned long *col_ind;            /*!< \brief Column index for each of the elements in val(). */
public:
  /*!
   * \brief Constructor of the class
   */
  CSparsityPattern();
  
  /*!
   * \brief Destructor of the class
   */
  ~CSparsityPattern();
  
  /*!
   * \brief Returns the index of the first element in row i.
   * \param[in] i - Index of the row.
   * \return Index of the first element.
   */
  inline unsigned long GetRowPointer(unsigned long i)  const {assert(i < nPoint+1); return row_ptr[i];}
  
  /*!
   * \brief Returns the index of the column of element j.
   * \param[in] j - Index of the element.
   * \return Column of the element with index j.
   */
  inline unsigned long GetColumnIndex(unsigned long j) const {assert(j < nnz); return col_ind[j];}
  
  /*!
   * \brief Returns the number of possible non-zero entries.
   * \return Number of non-zero entries.
   */
  inline unsigned long GetnNonZero()                   const {return nnz;}
  
  /*!
   * \brief Returns the index in the value array of a matrix entry.
   * \param[in] iPoint - Row index
   * \param[in] jPoint - Column index
   * \return The index where this matrix entry can be found in the value array.
   */
  unsigned long GetIndex(const unsigned long iPoint, const unsigned long jPoint);
  
};

/*!
 * \class CDualMeshSparsity
 * \brief Class for defining a sparsity pattern based on the dual mesh structure of the mesh. 
 * \ingroup crs_matrix
 * \author F. Palacios, T. Albring
 */
class CDualMeshSparsity : CSparsityPattern {
  
public:
  
  /*!
   * \brief Constructor of the class
   * \param[in] geometry - Geometrical definition
   * \param[in] fill_in - Fill-in size for ILU(p) preconditioner.
   */
  CDualMeshSparsity(CGeometry *geometry, unsigned short fill_in);
  
  /*!
   * \brief Destructor of the class
   */
  ~CDualMeshSparsity();
  
private:
  
  /*!
   * \brief Determines neighbouring points
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iPoint - Base point to compute neighbours.
   * \param[in] deep_level - Deep level for the recursive algorithm.
   * \param[in] fill_level - ILU fill in level.
   * \param[in] vneighs - Storage the neighbours points to iPoint.
   */
  void SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level,
                     vector<unsigned long> & vneighs);
  
};

/*!
 * \class CPrimalMeshSparsity
 * \brief Class for defining a sparsity pattern based on the primal mesh structure of the mesh. 
 * \ingroup crs_matrix
 * \author F. Palacios, T. Albring
 */
class CPrimalMeshSparsity : CSparsityPattern {
  
public:
  
  /*!
   * \brief Constructor of the class
   * \param[in] geometry - Geometrical definition
   * \param[in] fill_in - Fill-in size for ILU(p) preconditioner.
   */
  CPrimalMeshSparsity(CGeometry *geometry, unsigned short fill_in);
  
  /*!
   * \brief Destructor of the class
   */
  ~CPrimalMeshSparsity();
  
private:
  
  /*!
   * \brief Determines neighbouring points
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] iPoint - Base point to compute neighbours.
   * \param[in] deep_level - Deep level for the recursive algorithm.
   * \param[in] fill_level - ILU fill in level.
   * \param[in] vneighs - Storage the neighbours points to iPoint.
   */
  void SetNeighbours(CGeometry *geometry, unsigned long iPoint, unsigned short deep_level, unsigned short fill_level,
                     vector<unsigned long> & vneighs);
};

/*!
 * \class CILU0MeshSparsity
 * \brief Special class for defining a sparsity pattern for the ILU0 preconditioner. 
 * \ingroup crs_matrix
 * \author F. Palacios, T. Albring
 */
class CILU0Sparsity : CSparsityPattern {
  
public:
  
  /*!
   * \brief Constructor of the class
   * \param[in] geometry - Geometrical definition
   * \param[in] sparsity_base - Sparsity pattern of the base matrix.
   */
  CILU0Sparsity(CGeometry* geometry, CConfig* config, CSparsityPattern* sparsity_base);
  
  /*!
   * \brief Destructor of the class
   */
  ~CILU0Sparsity();
};

