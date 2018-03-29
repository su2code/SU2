/*!
 * \file matrix_interface_structure.hpp
 * \brief Headers of the main subroutines for creating the sparse matrices-by-blocks.
 *        The subroutines and functions are in the <i>matrix_structure.cpp</i> file.
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 6.0.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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
#include "matrix_structure.hpp"

template<class CalcType, class BaseType, template<class, class> class MatrixType >
class TCMatrixBase {
  
  typedef MatrixType<CalcType, BaseType> Matrix; 
  
private:
  unsigned short BlockSize;
  
  CalcType* BlockLin_CalcType;
  BaseType* BlockLin_BaseType;
  CalcType** Block_CalcType;
  BaseType** Block_BaseType;
  
  BaseType** ToBaseType(CalcType** block, unsigned short BlockSizeI, unsigned short BlockSizeJ);
  
  CalcType** ToCalcType(BaseType** block, unsigned short BlockSizeI, unsigned short BlockSizeJ);

  BaseType* ToBaseType(CalcType* block, unsigned short BlockSize);
  
  CalcType* ToCalcType(BaseType* block, unsigned short BlockSize);
  
  BaseType ToBaseType(CalcType val);
  
  CalcType ToCalcType(BaseType val);
  
  Matrix* Cast(){return static_cast<Matrix*>(this);}
  
public:
  
  TCMatrixBase();
  
  ~TCMatrixBase();
  
  virtual void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn,
                          CSparsityPattern *pattern);
  
  virtual void SetValZero() = 0;

  virtual void DeleteValsRowi(unsigned long i) = 0;
  
  virtual void SetVal2Diag(unsigned long iPoint, BaseType val) = 0;
  
  virtual void AddVal2Diag(unsigned long iPoint, BaseType val) = 0;
  
  virtual void AddBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block) = 0;
  
  virtual void SubtractBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block) = 0;
  
  virtual const BaseType* GetBlock(unsigned long iPoint, unsigned long jPoint) const = 0;
  
  virtual BaseType GetBlock(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar) const = 0;
  
  virtual void SetBlock(unsigned long iPoint, unsigned long jPoint, BaseType* block) = 0;
  
  virtual void SetBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block) = 0;
  
  virtual void MatrixVectorProduct() = 0;    
  
  virtual void Build_Preconditioner(unsigned short kind_prec, bool transpose) = 0;
};

template<class CalcType, class BaseType >
class TCMatrixInterface : TCMatrixBase<CalcType, BaseType, TCMatrixInterface>{
  
  typedef Matrix<CalcType, CalcType> Matrix;
  
public:
  
  TCMatrixInternal();
  
  ~TCMatrixInternal();
  
  void Initialize(unsigned long nPoint, unsigned long nPointDomain, unsigned short nVar, unsigned short nEqn,
                  CSparsityPattern *pattern);
  
  void SetValZero();

  void DeleteValsRowi(unsigned long i);
  
  void SetVal2Diag(unsigned long iPoint, BaseType val);
  
  void AddVal2Diag(unsigned long iPoint, BaseType val);
 
  void AddBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  
  void SubtractBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  
  const BaseType* GetBlock(unsigned long iPoint, unsigned long jPoint) const;
  
  BaseType GetBlock(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar) const;
  
  void SetBlock(unsigned long iPoint, unsigned long jPoint, BaseType* block);
  
  void SetBlock(unsigned long iPoint, unsigned long jPoint, BaseType** block);
  
  void MatrixVectorProduct();    
  
  void Build_Preconditioner(unsigned short kind_prec, bool transpose);
  
};

#include "matrix_interface_structure.inl"