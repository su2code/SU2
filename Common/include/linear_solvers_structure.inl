/*!
 * \file linear_solvers_structure.inl
 * \brief inline subroutines of the <i>linear_solvers_structure.hpp</i> file.
 * \author J. Hicken, F. Palacios, T. Economon
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
#include "linear_solvers_structure.hpp"


template<class CalcType>
inline CalcType TCLinSolver_FGMRES<CalcType>::Sign(const CalcType & x, const CalcType & y) const {
  if (y == 0.0)
    return 0.0;
  else {
//    return (y < 0 ? -fabs(x) : fabs(x));
    if (y < 0) return -fabs(x);
    else return fabs(x);
  }
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetValZero_Matrix(){
  Matrix.SetValZero();
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetValZero_Rhs(){
  LinSysRes.SetValZero();
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetValZero_Sol(){
  LinSysSol.SetValZero();
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::DeleteValsRowi(unsigned long i){
  Matrix.DeleteValsRowi(i);
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetVal2Diag_Matrix(unsigned long iPoint, BaseType val){
  Matrix.SetVal2Diag(iPoint, convert.ToCalcType(val));
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::AddVal2Diag_Matrix(unsigned long iPoint, BaseType val){
  Matrix.AddVal2Diag(iPoint,convert.ToCalcType(val));
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::AddBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block){
  Matrix.AddBlock(iPoint, jPoint, convert.ToCalcType(block, BlockSize, BlockSize));
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SubtractBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block){
  Matrix.SubtractBlock(iPoint, jPoint,  convert.ToCalcType(block, BlockSize, BlockSize));
}

template<class CalcType, class BaseType>
inline BaseType* TCSysSolve<CalcType, BaseType>::GetBlock_Matrix(unsigned long iPoint, unsigned long jPoint){
  return convert.ToBaseType(Matrix.GetBlock(iPoint, jPoint), BlockSize*BlockSize);
}

template<class CalcType, class BaseType>
inline BaseType TCSysSolve<CalcType, BaseType>::GetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, unsigned short iVar, unsigned short jVar){
 return convert.ToBaseType(Matrix.GetBlock(iPoint, jPoint, iVar, jVar));
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType* block){
  Matrix.SetBlock(iPoint, jPoint, convert.ToCalcType(block, BlockSize*BlockSize));
}

template<class CalcType, class BaseType>
inline void TCSysSolve<CalcType, BaseType>::SetBlock_Matrix(unsigned long iPoint, unsigned long jPoint, BaseType** block){
  Matrix.SetBlock(iPoint, jPoint, convert.ToCalcType(block, BlockSize, BlockSize));
}


template<>
inline su2double** Convert<su2double, su2double>::ToBaseType(su2double **block, unsigned short BlockSizeI, unsigned short BlockSizeJ){
  return block;
}

template<>
inline su2double** Convert<su2double, su2double>::ToCalcType(su2double **block, unsigned short BlockSizeI, unsigned short BlockSizeJ){
  return block;
}

template<>
inline su2double* Convert<su2double, su2double>::ToBaseType(su2double *block, unsigned short BlockSize){
  return block;
}

template<>
inline su2double* Convert<su2double, su2double>::ToCalcType(su2double *block, unsigned short BlockSize){
  return block;
}

template<>
inline su2double Convert<su2double, su2double>::ToBaseType(su2double val){
  return val;
}

template<>
inline su2double Convert<su2double, su2double>::ToCalcType(su2double val){
  return val;
}

#ifdef CODI_REVERSE_TYPE
template<>
inline su2double** Convert<passivedouble, su2double>::ToBaseType(passivedouble **block, unsigned short BlockSizeI, unsigned short BlockSizeJ){
  for (unsigned short i = 0; i < BlockSizeI; i++){
    for (unsigned short j = 0; j< BlockSizeJ; j++){    
    Block_BaseType[i][j] = su2double(block[i][j]);
    }
  }
  return Block_BaseType;
}

template<>
inline passivedouble** Convert<passivedouble, su2double>::ToCalcType(su2double **block, unsigned short BlockSizeI, unsigned short BlockSizeJ){
  for (unsigned short i = 0; i < BlockSizeI; i++){
    for (unsigned short j = 0; j< BlockSizeJ; j++){    
    Block_CalcType[i][j] = SU2_TYPE::GetValue(block[i][j]);
    }
  }
  return Block_CalcType;
}

template<>
inline su2double* Convert<passivedouble, su2double>::ToBaseType(passivedouble *block, unsigned short BlockSize){
  for (unsigned short i = 0; i < BlockSize; i++){
    BlockLin_BaseType[i] = su2double(block[i]);
  }
  return BlockLin_BaseType;
}

template<>
inline passivedouble* Convert<passivedouble, su2double>::ToCalcType(su2double *block, unsigned short BlockSize){
  for (unsigned short i = 0; i < BlockSize; i++){
    BlockLin_CalcType[i] = SU2_TYPE::GetValue(block[i]);
  }
  return BlockLin_CalcType;
}

template<>
inline su2double Convert<passivedouble, su2double>::ToBaseType(passivedouble val){
  return su2double(val);
}

template<>
inline passivedouble Convert<passivedouble, su2double>::ToCalcType(su2double val){
  return  SU2_TYPE::GetValue(val);
}
#endif

