/*!
 * \file vector_interface_structure.hpp
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

template<class CalcType>
class TCVectorInterface {
  
public:
  TCVectorInterface();
  
  ~TCVectorInterface();
  
  virtual void SetValZero() = 0;
  
  virtual void Initialize(const unsigned long & numBlk, const unsigned long & numBlkDomain, const unsigned short & numVar, const CalcType & val = 0.0) = 0;
  
  virtual void Equals_AX(const CalcType & a, TCVectorInterface<CalcType> & x) = 0;
  
  virtual void Plus_AX(const CalcType & a, TCVectorInterface<CalcType> & x) = 0;
  
  virtual void Equals_AX_Plus_BY(const CalcType & a, TCVectorInterface<CalcType> & x, const CalcType & b, TCVectorInterface<CalcType> & y) = 0;
  
  virtual TCVectorInterface<CalcType> & operator=(const TCVectorInterface<CalcType> & u) = 0;
  
  virtual TCVectorInterface<CalcType> & operator=(const CalcType & val) = 0;
  
  virtual TCVectorInterface<CalcType> & operator+(const TCVectorInterface<CalcType> & u) const = 0;
  
  virtual TCVectorInterface<CalcType> & operator+=(const TCVectorInterface<CalcType> & u) = 0;
  
  virtual TCVectorInterface<CalcType> & operator-(const TCVectorInterface<CalcType> & u) const = 0;

  virtual TCVectorInterface<CalcType> & operator-=(const TCVectorInterface<CalcType> & u) = 0;
  
  virtual TCVectorInterface<CalcType> & operator*(const CalcType & val) const = 0;
  
  virtual TCVectorInterface<CalcType> & operator*=(const CalcType & val) = 0;
  
  virtual TCVectorInterface<CalcType> & operator/(const CalcType & val) const = 0;
  
  virtual TCVectorInterface<CalcType> & operator/=(const CalcType & val) = 0;
  
  virtual CalcType & operator[](const unsigned long & i) = 0;
  
  virtual const CalcType & operator[](const unsigned long & i) const = 0;
  
  virtual CalcType norm() const = 0;
  
  virtual void CopyToArray(CalcType* u_array) = 0;
  
  virtual void SubtractBlock(unsigned long val_ipoint, CalcType *val_residual) = 0;
  
  virtual void AddBlock(unsigned long val_ipoint, CalcType *val_residual) = 0;

  virtual void SetBlock(unsigned long val_ipoint, unsigned short val_var, CalcType val_residual) = 0;
  
  virtual void SetBlock(unsigned long val_ipoint, CalcType *val_residual) = 0;
  
  virtual void SetBlock_Zero(unsigned long val_ipoint) = 0;
  
  virtual void SetBlock_Zero(unsigned long val_ipoint, unsigned short val_var) = 0;
	
  virtual const CalcType *GetBlock(unsigned long val_ipoint) = 0;
	
  virtual CalcType GetBlock(unsigned long val_ipoint, unsigned short val_var);
  
  template<class T>
  friend T dotProd(const TCVectorInterface<T> & u, const TCVectorInterface<T> & v);
  
  virtual void SendReceive(CGeometry *geometry, CConfig *config) = 0;
  
  virtual void SendReceive_Reverse(CGeometry *geometry, CConfig *config) = 0;
};