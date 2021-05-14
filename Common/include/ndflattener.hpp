/*!
 * \file ndflattener.hpp
 * \brief Flatten pointer-to-pointer-... arrays for MPI communication
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

#pragma once

#include <iostream>
#include <utility>
#include <cassert>
#include <sstream>



template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
struct MPI_Environment {
  MPI_Allgatherv_type MPI_Allgatherv;
  MPI_Datatype_type mpi_data;
  MPI_Datatype_type mpi_index;
  MPI_Communicator_type comm;
  int rank; int size;
};

/*!
 * \class NdFlattener
 * \brief Serialize pointer-to-pointer-... array into one 1D array, keeping track
 * of the offsets in few additional 1D arrays.
 *
 * The pointer-to-pointer-... array can be provided by a nested lambda function
 * ('recursive function') or by gathering such arrays from MPI processes ('collective
 * communication'). After initializing an NdFlattener with either of these data, 
 * it can be refreshed in the same way after the the pointer-to-pointer-... array's
 * values (but not its structure) have changed.
 *
 * \tparam Data - Type of stored array data
 * \tparam K - number of indices
 * \tparam Index - Type of index
 */

/* --- Implementation details---
 * If your array has K indices, instantiate this class with template parameter K,
 * which is derived recursively from this class with template parameter (K-1). 
 * In each layer, there is an array: Of type Data for K=1, and of type Index for K>1.
 *
 * The data array of K=1 contains the values of the array A in lexicographic ordering:
 * [0]...[0][0][0], [0]...[0][0][1], ..., [0]...[0][0][something], 
 * [0]...[0][1][0], [0]...[0][1][1], ..., [0]...[0][1][something],
 * ...,
 * [0]...[1][0][0], [0]...[1][0][1], ..., [0]...[1][0][something],
 * ... ...,
 * [1]...[0][0][0], [1]...[0][0][1], ..., [1]...[0][0][something],
 * ... ... ...
 * Let us call each row in this representation a "last-index list".
 * Note that they might have different lengths, so "something" might stand for different
 * values here. Last-index lists can also be empty.
 *
 * The indices array of K=2 contains the indices of the (K=1)-data array at which a new
 * last-index list starts. If a last-index list is empty, the are repetitive entries in the 
 * indices array.
 *
 * The indices array of K=3 contains the indices of the (K=2)-indices array at which the
 * last-but-two index increases by one or drops to zero. If the last-but-two index is increased
 * by more than one, there are repetitive entries.
 *
 * Etc. etc, up to the indices array at layer K.
 *
 * To form such a structure, we typically iterate twice through the pointer-to-pointer-... 
 * array: The first time we get to know how much space to reserve in each layer, then we 
 * allocate it and fill it with data during the second iteration.
 */

template<typename Data, size_t K, typename Index>
class NdFlattener: public NdFlattener<Data,K-1,Index>{
public:
  typedef NdFlattener<Data,K-1,Index> Base; // could also be named LowerLayer
  typedef NdFlattener<Data,K,Index> CurrentLayer;
  typedef typename Base::LowestLayer LowestLayer; // the K=1 class

private:
  /*! \brief Number of nodes in this layer.
   * 
   * For the layer K=1, nNodes will be the number of data points.
   * For a layer K>1, nNodes will be the number of sublists.
   */
  Index nNodes=0; 

  /*! \brief Iterator used at construction, runs from 0 to (nNodes-1). */
  Index iNode=0;

  /*! \brief Indices in the lower layer's indices or data array */
  Index* indices=nullptr;

  /*! \brief Return a const pointer to the base class (one layer down).
   */
  Base const* cast_to_Base() const {
    return (Base const*)this;
  }

  /*=== Outputting ===*/

public:
  /*! \brief Return a python-style list string.
   *
   * Like this: [[1, 2], [10, 20, 30]]
   */
  std::string toPythonString() const {
    return toPythonString_fromto(0, getNChildren());
  }
protected:
  /*! \brief Return a python-style list string, using the data of the
   * indices array between 'from' (inclusive) and 'to' (exclusive).
   *
   * Like this: [[1, 2], [10, 20, 30]]
   * \param[in] from - Beginning of the representation in the indices array.
   * \param[in] to - Ending of the representation in the indices array.
   */
  std::string toPythonString_fromto(Index from, Index to) const {
    std::string result = "[";
    for(Index i=from; i<to; ){
      result += Base::toPythonString_fromto(indices[i], indices[i+1]);
      if(++i<to) result += ", ";
    }
    result += "]";
    return result;
  }

public:

  /*! \brief Basic constructor.
   * 
   * Called recursively when a derived class (higher K) is constructed.
   * Allocation is done later with initialize_or_refresh.
   */
  NdFlattener() {};

  /*! \brief Initialize or refresh the NdFlattener.
   * \details Either a 'recursive function' or 'collective communication'
   * may be used. When the NdFlattener does not hold data yet, it is 
   * initialized, meaning that the data are collected and the indices arrays
   * are allocated and filled. Otherwise it is refreshed, meaning that the data are
   * recollected under the assumption that the indices arrays did not change.
   */
  template<class ...ARGS>
  void initialize_or_refresh(ARGS... args){
    if( initialized() ){
      refresh(args...);
    } else {
      initialize(args...);
    }
  }

  /*! \brief Initialization status of the NdFlattener.
   * \returns true if the NdFlattener has been initialized
   */
  bool initialized(){
    return (indices!=nullptr);
  }


  virtual ~NdFlattener(void) {
    if(indices!=nullptr)
      delete[] indices;
  }


protected:
  /*! \brief Allocate the indices array after \a nNodes has been determined.
   */
  void allocate() {
    indices = new Index[nNodes+1];
    indices[0] = 0;
    Base::allocate();
  }

  /*! \brief Set \a iNode to 0 in all layers.
   */
  void reset_iNode(){
    iNode = 0;
    Base::reset_iNode();
  }

  /*=== Construct from 'recursive function' ===*/
public:
  /*! \brief Initialize from a 'recursive function'.
   *
   * The function should return a pair. Its first entry is the number of children. 
   * Its second entry is a function with the same meaning, recursively 
   * one layer down.
   * \param f - the 'recursive function'
   */
  template<class f_type>
  void initialize(f_type f) {
    count_f(f);
    allocate();
    set_f(f, false);
  }

  /*! \brief Refresh the data according to the 'recursive function'
   *
   * The NdFlattener must have been constructed with a 'recursive function'.
   * Now refresh the values with another 'recursive function'. The subarray lengths 
   * resulting from both 'recursive functions' must coincide, as the indices arrays
   * are not changed.
   * 
   * \param f - the 'recursive function'
   * \tparam f_type - to allow for any type of the 'recursive function'
   */
  template<class f_type>
  void refresh(f_type f){
    reset_iNode();
    set_f(f, true);
  }

protected:
  /*! \brief Determine the space required for reading the 'recursive function'.
   *
   * \param f - the 'recursive function'
   */
  template<class f_type>
  void count_f(f_type f) {
    Index nChild = f.first;
    for(Index iChild=0; iChild<nChild; iChild++){
      nNodes++;
      Base::count_f(f.second(iChild));
    }
  }  

  /*! \brief Read the 'recursive function' into the allocated arrays.
   *
   * \param f - the 'recursive function'
   * \param refresh - if true, the object is already initialized and only the data
   *   in layer 1 have to be overwritten
   * \tparam f_type - to allow for any type of the 'recursive function'
   */
  template<class f_type>
  void set_f(f_type f, bool refresh) {
    Index nChild = f.first;
    for(Index iChild=0; iChild<nChild; iChild++){
      Base::set_f(f.second(iChild), refresh);
      if(!refresh){
        indices[iNode+1] = indices[iNode] + f.second(iChild).first;
      } else {
        assert( indices[iNode+1] == indices[iNode] + f.second(iChild).first );
      }
      iNode++;
    }
  }


  /*=== Construct with Allgatherv ===*/

public:
  /*! \brief Initialize a flattener with K indices by combining distributed flatteners with (K-1) indices each.
   *
   * The new first index will encode the rank of the process. Data is exchanged in MPI::Allgatherv-style
   * collective communication.
   * \param[in] mpi_env - The MPI environment used for communication.
   * \param[in] local_version - The local NdFlattener structure with (K-1) indices.
   */
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void initialize( 
    MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> const& mpi_env, 
    Base const* local_version 
  ) {
    Index** Nodes_all = new Index*[K]; // [k][r] is number of all nodes in layer (k+1), rank r in the new structure
    for(size_t k=0; k<K; k++)
      Nodes_all[k] = nullptr;
    Nodes_all[K-1] = new Index[mpi_env.size]; // {1, 1, ..., 1}
    int* displs = new int[mpi_env.size]; // {0, 1, ..., size-1}
    int* ones = new int[mpi_env.size]; // {1,1, ...}
    for(int r=0; r<mpi_env.size; r++){
      nNodes += Nodes_all[K-1][r] = 1;
      displs[r] = r;
      ones[r] = 1;
    }
    Base::count_g(mpi_env, Nodes_all, local_version, displs, ones); // set the lower layers' nNodes and Nodes_all[k]

    allocate();

    indices[0] = 0;
    for(int r=0; r<mpi_env.size; r++){
      indices[r+1] = indices[r] + Nodes_all[K-2][r];
    }
    Base::set_g(mpi_env, Nodes_all, local_version);
    
    for(size_t k=0; k<K; k++){
      delete[] Nodes_all[k];
    }
    delete[] Nodes_all;
    delete[] displs;
    delete[] ones;
  }

  /*! \brief Refresh the data by MPI collective communication.
   *
   * The NdFlattener must have been constructed by MPI collective communication.
   * Now refresh the values with another collective communication. The subarray lengths 
   * resulting from both collective communications must coincide, as the indices arrays
   * are not changed.
   * \param[in] mpi_env - The MPI environment used for communication.
   * \param[in] local_version - The local NdFlattener structure.
   */
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void refresh( 
    MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> const& mpi_env, 
    Base const* local_version 
  ) {
    Index* Nodes_all_0 = nullptr;
    int* displs = new int[mpi_env.size]; // {0, 1, ..., size-1}
    int* ones = new int[mpi_env.size]; // {1,1, ...}
    for(int r=0; r<mpi_env.size; r++){
      displs[r] = r;
      ones[r] = 1;
    }
    LowestLayer::count_g(mpi_env, &Nodes_all_0, (LowestLayer const*)local_version, displs, ones); 
    LowestLayer::set_g(mpi_env, &Nodes_all_0, (LowestLayer const*)local_version);

    delete[] Nodes_all_0; // allocated by count_g
    delete[] displs;
    delete[] ones;
  }

protected:
  /*! \brief Count the distributed flatteners' numbers of nodes, and set nNodes.
   *
   * \param[in] mpi_env - MPI environment for communication
   * \param[out] Nodes_all - [k][r] is set to number of nodes in layer (k+1), rank r.
   * \param[in] local_version - local instance to be send to the other processes
   * \param[in] displs - {0,1,...,size-1}
   * \param[in] ones - {1,1,...,1}
   */
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void count_g(MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> const& mpi_env, 
         Index** Nodes_all,
         CurrentLayer const* local_version,
         int const* displs, int const* ones )
  { 
    assert( Nodes_all[K-1]==nullptr);
    Nodes_all[K-1] = new Index[mpi_env.size];
    nNodes = 0;
    // gather numbers of nodes in the current layer from all processes
    mpi_env.MPI_Allgatherv( &(local_version->nNodes), 1, mpi_env.mpi_index, Nodes_all[K-1], ones, displs, mpi_env.mpi_index, mpi_env.comm );
    for(int r=0; r<mpi_env.size; r++){
      nNodes += Nodes_all[K-1][r];
    }
    Base::count_g(mpi_env, Nodes_all, local_version->cast_to_Base(), displs, ones);
  }

  /*! \brief Gather the distributed flatteners' data and index arrays into the allocated arrays.
   *
   * \param[in] mpi_env - MPI environment for communication
   * \param[in] Nodes_all - [k][r] is the number of nodes in layer (k+1), rank r.
   * \param[in] local_version - local instance to be send to the other processes
   */
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void set_g(MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> const& mpi_env, 
         Index** Nodes_all,
         CurrentLayer const* local_version )
  { 

    int* Nodes_all_K_as_int = new int[mpi_env.size];
    int* Nodes_all_k_cumulated = new int[mpi_env.size+1]; // [r] is number of nodes in the current layer, summed over all processes with rank below r 
    // plus one. Used as displacements in Allgatherv, but we do not want to transfer the initial zeros and rather the last element of indices, 
    // which is the local nNodes of the layer below. Note that MPI needs indices of type 'int'.
    Nodes_all_k_cumulated[0] = 1;
    for(int r=0; r<mpi_env.size; r++){
      Nodes_all_k_cumulated[r+1] = Nodes_all_k_cumulated[r] + Nodes_all[K-1][r];
      Nodes_all_K_as_int[r] = Nodes_all[K-1][r];
    }
    mpi_env.MPI_Allgatherv( local_version->indices+1, Nodes_all[K-1][mpi_env.rank], mpi_env.mpi_index, indices, Nodes_all_K_as_int, Nodes_all_k_cumulated, mpi_env.mpi_index, mpi_env.comm );
    // shift indices 
    for(int r=1; r<mpi_env.size; r++){
      Index first_entry_to_be_shifted = Nodes_all_k_cumulated[r];
      Index last_entry_to_be_shifted = Nodes_all_k_cumulated[r+1]-1;
      Index shift = indices[ first_entry_to_be_shifted - 1];
      for(Index i=first_entry_to_be_shifted; i<=last_entry_to_be_shifted; i++){
        indices[ i ] += shift;
      }
    }
    delete[] Nodes_all_K_as_int;
    delete[] Nodes_all_k_cumulated;

    Base::set_g(mpi_env, Nodes_all, local_version->cast_to_Base());
  }

    
    

  /*=== Access to data and numbers of children ===*/
protected:
  template<class ...ARGS>
  Data get_withoffset(Index offset, Index i1, ARGS... i2) const {
    return Base::get_withoffset(indices[offset+i1], i2...);
  }
  template<class ...ARGS>
  Index getNChildren_withoffset(Index offset, Index i1, ARGS... i2) const {
    return Base::getNChildren_withoffset(indices[offset+i1], i2...);
  }
  Index getNChildren_withoffset(Index offset, Index i) const {
    return indices[offset+i+1] - indices[offset+i];
  }
public:
  /*! \brief Index look-up.
   *
   * The number of parameters must be K.
   */
  template<class ...ARGS>
  Data get(ARGS... i) const {
    return get_withoffset(0, i...);
  }
  /*! \brief Look-up of length of the next-layer sublist.
   *
   * Specify less than K indices. When the function returns N,
   * the next index must lie inside {0, 1, ..., N-1}.
   */
  template<class ...ARGS>
  Index getNChildren(ARGS... i) const {
    return getNChildren_withoffset(0, i...);
  }
  Index getNChildren() const { // should not be called by recursion
    return nNodes;
  }


};

template<typename Data, typename Index>
class NdFlattener<Data, 1, Index> {
public:
  typedef NdFlattener<Data, 1, Index> CurrentLayer;
  typedef CurrentLayer LowestLayer;

private:
  Index nNodes=0;
  Index iNode=0;
  Data* data=nullptr;

  /*=== Outputting ===*/
protected:
  std::string toPythonString_fromto(Index from, Index to) const {
    std::stringstream result;
    result  << "[";
    for(Index i=from; i<to; ){
      result << data[i];
      if(++i<to) result << ", ";
    }
    result << "]";
    return result.str();
  }

public:
  NdFlattener(void) {}

  virtual ~NdFlattener(void) {
    if(data!=nullptr)
      delete[] data;
  }

protected:
  void allocate(){
    data = new Data[nNodes];
  }

  void reset_iNode(){
    iNode = 0;
  }

  /*=== Construct from 'recursive function' ===*/
protected:
  template<typename f_type>
  void count_f(f_type f){
    nNodes += f.first;
  }

  template<typename f_type>
  void set_f(f_type f, bool refresh){
    Index nChild = f.first;
    for(Index iChild=0; iChild<nChild; iChild++){
      data[iNode] = f.second(iChild);
      iNode++;
    }
  }
  
  /*=== Construct with Allgatherv ===*/
protected:
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void count_g(MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> const& mpi_env, 
         Index** Nodes_all,
         CurrentLayer const* local_version,
         int const* displs, int const* ones)
  { 
    assert( Nodes_all[0]==nullptr);
    Nodes_all[0] = new Index[mpi_env.size];
    nNodes = 0;
    // gather numbers of nodes in the current layer from all processes
    mpi_env.MPI_Allgatherv( &(local_version->nNodes), 1, mpi_env.mpi_index, Nodes_all[0], ones, displs, mpi_env.mpi_index, mpi_env.comm );
    for(int r=0; r<mpi_env.size; r++){
      nNodes += Nodes_all[0][r];
    }
  }

  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  void set_g(MPI_Environment<MPI_Allgatherv_type, MPI_Datatype_type, MPI_Communicator_type> mpi_env, 
         Index** Nodes_all,
         CurrentLayer const* local_version )
  { 


    int* Nodes_all_0_as_int = new int[mpi_env.size];
    int* Nodes_all_0_cumulated = new int[mpi_env.size+1];
    Nodes_all_0_cumulated[0] = 0;
    for(int r=0; r<mpi_env.size; r++){
      Nodes_all_0_cumulated[r+1] = Nodes_all_0_cumulated[r] + Nodes_all[0][r];
      Nodes_all_0_as_int[r] = Nodes_all[0][r];
    }

    mpi_env.MPI_Allgatherv( local_version->data, Nodes_all[0][mpi_env.rank], mpi_env.mpi_data, data, Nodes_all_0_as_int, Nodes_all_0_cumulated, mpi_env.mpi_data, mpi_env.comm );
    delete[] Nodes_all_0_as_int;
    delete[] Nodes_all_0_cumulated;
  }


  /*== Access to data and numbers of children ==*/
protected:
  Data get_withoffset(Index offset, Index i) const{
    return data[offset+i];
  }
public:
  Data get(Index i) const{
    return get_withoffset(0, i);
  }


};


