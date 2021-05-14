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
#include <vector>
#include "../parallelization/mpi_structure.hpp"

template<size_t K, typename Data=su2double, typename Index=unsigned long>
class NdFlattener;

namespace helpers {
  template<typename MPI_Allgatherv_type, typename MPI_Datatype_type, typename MPI_Communicator_type>
  struct NdFlattener_MPI_Environment {
    MPI_Allgatherv_type MPI_Allgatherv;
    MPI_Datatype_type mpi_data;
    MPI_Datatype_type mpi_index;
    MPI_Communicator_type comm;
    int rank; int size;
  };

  /*! \class IndexAccumulator
   * \brief Data structure holding an offset for the NdFlattener, to provide a []...[]-interface.
   * \details Derived from IndexAccumulator_Base, specifying the operator[] method:
   *  - For N==2, the structure has already read all indices but two. So after this method has read the last-but-one
   *    index, return a pointer to the corresponding section of the data array in layer K=1.
   *  - For N>2, more indices have to be read, return an IndexAccumulator<N-1>.
   * \tparam N - Number of missing parameters
   * \tparam K - number of indices of accessed NdFlattener
   * \tparam Data - Data type of accessed NdFlattener
   * \tparam Index - Index type of accessed NdFlattener
   */
  /*! \class IndexAccumulator_Base
   * \brief Parent class of IndexAccumulator.
   * \details IndexAccumulator provides the operator[] method.
   */
  template<size_t N, size_t K, typename Data, typename Index>
  class IndexAccumulator_Base {
  protected:
    typedef NdFlattener<N,Data,Index> Nd_type; /*!< \brief Access this level of the NdFlattener next. */
    const Nd_type* nd; /*!< \brief The accessed NdFlattener. */
    const Index offset; /*!< \brief Index in the currently accessed layer. */
    IndexAccumulator_Base(const Nd_type* nd, Index offset): nd(nd), offset(offset) {}
  };
  template<size_t N, size_t K, typename Data, typename Index>
  class IndexAccumulator : public IndexAccumulator_Base<N,K,Data,Index>{
  public:
    typedef IndexAccumulator_Base<N,K,Data,Index> Base;
    template<class ...ARGS> IndexAccumulator(ARGS... args): Base(args...) {};
    typedef IndexAccumulator<N-1,K,Data,Index> LookupType; /*!< Return type of operator[] */
    typedef const LookupType LookupType_const; /*!< Return type of operator[], const version */
    /* \brief Read one more index. */
    /*LookupType operator[] (Index i){
      return LookupType(static_cast<typename Base::Nd_type::Base*>(this->nd),this->nd->GetIndices()[this->offset+i]);
    }*/
    /*! \brief Read one more index, const version. */
    LookupType_const operator[] (Index i) const {
      return LookupType(static_cast<const typename Base::Nd_type::Base*>(this->nd),this->nd->GetIndices()[this->offset+i]);
    }
  };
  template<size_t K, typename Data, typename Index>
  class IndexAccumulator<2,K,Data,Index> : public IndexAccumulator_Base<2,K,Data,Index>{
  public:
    typedef IndexAccumulator_Base<2,K,Data,Index> Base;
    template<class ...ARGS> IndexAccumulator(ARGS... args): Base(args...) {};
    typedef Data* LookupType; /*!< Return type of operator[] */
    typedef const Data* LookupType_const; /*!< Return type of operator[], const version */
    /* \brief Read the last-but-one index.
     * \return Pointer to the corresponding section in the data array in layer K=1.
     */
    /*LookupType operator[] (Index i){
      return static_cast<const typename Base::Nd_type::Base*>(this->nd)->GetData() + this->nd->GetIndices()[this->offset+i];
    }*/
    /*! \brief Read the last-but-one index, const version.
     * \return Const pointer to the corresponding section in the data array in layer K=1.
     */
    LookupType_const operator[] (Index i) const {
      return static_cast<const typename Base::Nd_type::Base*>(this->nd)->GetData() + this->nd->GetIndices()[this->offset+i];
    }
  };

} // namespace helpers

static helpers::NdFlattener_MPI_Environment<decltype(&(SU2_MPI::Allgatherv)), decltype(MPI_INT), decltype(SU2_MPI::GetComm())>
Get_Nd_MPI_Env() {
  helpers::NdFlattener_MPI_Environment<decltype(&(SU2_MPI::Allgatherv)), decltype(MPI_INT), decltype(SU2_MPI::GetComm())> mpi_env;
  mpi_env.MPI_Allgatherv = &(SU2_MPI::Allgatherv);
  mpi_env.mpi_index = MPI_UNSIGNED_LONG;
  mpi_env.mpi_data = MPI_DOUBLE;
  mpi_env.comm=SU2_MPI::GetComm();
  SU2_MPI::Comm_rank(mpi_env.comm, &(mpi_env.rank));
  SU2_MPI::Comm_size(mpi_env.comm, &(mpi_env.size));
  return mpi_env;
}

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
 * \tparam K - number of indices
 * \tparam Data - Type of stored array data
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

template<size_t K, typename Data, typename Index>
class NdFlattener: public NdFlattener<K-1,Data,Index>{

public:
  typedef NdFlattener<K-1,Data,Index> Base; // could also be named LowerLayer
  typedef NdFlattener<K,Data,Index> CurrentLayer;
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
  std::vector<Index> indices;

  /*=== Getters ===*/
public:
  Index* GetIndices() {return indices.data();}
  const Index* GetIndices() const {return indices.data();}

  /*=== Outputting ===*/

public:
  /*! \brief Write in Python-list style.
   *
   * Like this: [[1, 2], [10, 20, 30]]
   */
  friend std::ostream& operator<<(std::ostream& output, NdFlattener const& nd) {
    nd.toPythonString_fromto(output, 0, nd.getNChildren());
    return output;
  }

protected:
  /*! \brief Write to stream in Python-list style, using the data of the
   * indices array between 'from' (inclusive) and 'to' (exclusive).
   *
   * Like this: [[1, 2], [10, 20, 30]]
   * \param[in] output - Stream
   * \param[in] from - Beginning of the representation in the indices array.
   * \param[in] to - Ending of the representation in the indices array.
   */
  void toPythonString_fromto(std::ostream& output, Index from, Index to) const {
    output << "[";
    for(Index i=from; i<to; ){
      Base::toPythonString_fromto(output, indices[i], indices[i+1]);
      if(++i<to) output << ", ";
    }
    output << "]";
  }

public:

  /*! \brief Basic constructor. Afterwards, initialization can be done with initialize_or_refresh.
   * 
   * Called recursively when a derived class (higher K) is constructed.
   */
  NdFlattener() {};

  /*! \brief Constructor which calls initialize_or_refresh.
   *
   */
  template<class... ARGS>
  NdFlattener(ARGS... args) {
    initialize_or_refresh(args...);
  };

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
    return indices.size()>0;
  }

protected:
  /*! \brief Allocate the indices array after \a nNodes has been determined.
   */
  void allocate() {
    indices.reserve(nNodes+1);
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
  template<typename MPI_Environment_type>
  void initialize( 
    MPI_Environment_type const& mpi_env,
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
  template<typename MPI_Environment_type>
  void refresh( 
    MPI_Environment_type const& mpi_env,
    Base const* local_version 
  ) {
    Index* Nodes_all_0 = nullptr;
    int* displs = new int[mpi_env.size]; // {0, 1, ..., size-1}
    int* ones = new int[mpi_env.size]; // {1,1, ...}
    for(int r=0; r<mpi_env.size; r++){
      displs[r] = r;
      ones[r] = 1;
    }
    LowestLayer::count_g(mpi_env, &Nodes_all_0, static_cast<LowestLayer const*>(local_version), displs, ones);
    LowestLayer::set_g(mpi_env, &Nodes_all_0, static_cast<LowestLayer const*>(local_version));

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
  template<typename MPI_Environment_type>
  void count_g(MPI_Environment_type const& mpi_env,
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
    Base::count_g(mpi_env, Nodes_all, static_cast<const Base*>(local_version), displs, ones);
  }

  /*! \brief Gather the distributed flatteners' data and index arrays into the allocated arrays.
   *
   * \param[in] mpi_env - MPI environment for communication
   * \param[in] Nodes_all - [k][r] is the number of nodes in layer (k+1), rank r.
   * \param[in] local_version - local instance to be send to the other processes
   */
  template<typename MPI_Environment_type>
  void set_g(MPI_Environment_type const& mpi_env,
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
    mpi_env.MPI_Allgatherv( local_version->indices.data()+1, Nodes_all[K-1][mpi_env.rank], mpi_env.mpi_index, indices.data(), Nodes_all_K_as_int, Nodes_all_k_cumulated, mpi_env.mpi_index, mpi_env.comm );
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

    Base::set_g(mpi_env, Nodes_all, static_cast<const Base*>(local_version));
  }

    
    

  /*=== Access to data and numbers of children ===*/
  // There are two ways to access data:
  // - nd.get(i_1, ..., i_K) for reading
  // - nd[i_1]...[i_{K-1}] returns a pointer to the data array, so you can
  //   read and write nd[i_1]...[i_K].
  // If you have indices i_1, ..., i_k and are interested in the bound for i_{k+1},
  // use getNChildren(i_1,...,i_k).
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
  /*! \brief Reading access.
   *
   * The number of parameters must be K.
   * If you also need write access, use the nd[i1]...[iK] interface.
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

  /* \brief Look-up with IndexAccumulator.
   */
  /*typename helpers::IndexAccumulator<K,K,Data,Index>::LookupType operator[](Index i0){
    return helpers::IndexAccumulator<K,K,Data,Index>(this,0)[i0];
  }*/
  /*! \brief Look-up with IndexAccumulator, const version.
   */
  typename helpers::IndexAccumulator<K,K,Data,Index>::LookupType_const operator[](Index i0) const {
    return helpers::IndexAccumulator<K,K,Data,Index>(this,0)[i0];
  }
};

template<typename Data, typename Index>
class NdFlattener<1, Data, Index> {
public:
  typedef NdFlattener<1, Data, Index> CurrentLayer;
  typedef CurrentLayer LowestLayer;

private:
  Index nNodes=0;
  Index iNode=0;
  std::vector<Data> data;


  /*=== Getters ===*/
public:
  Data* GetData() {return data.data();}
  const Data* GetData() const {return data.data();}

  /*=== Outputting ===*/
protected:
  void toPythonString_fromto(std::ostream& output, Index from, Index to) const {
    output  << "[";
    for(Index i=from; i<to; ){
      output << data[i];
      if(++i<to) output << ", ";
    }
    output << "]";
  }

public:
  NdFlattener(void) {}

protected:
  void allocate(){
    data.reserve(nNodes);
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
  template<typename MPI_Environment_type>
  void count_g(MPI_Environment_type const& mpi_env,
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

  template<typename MPI_Environment_type>
  void set_g(MPI_Environment_type const& mpi_env,
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

    mpi_env.MPI_Allgatherv( local_version->data.data(), Nodes_all[0][mpi_env.rank], mpi_env.mpi_data, data.data(), Nodes_all_0_as_int, Nodes_all_0_cumulated, mpi_env.mpi_data, mpi_env.comm );
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




