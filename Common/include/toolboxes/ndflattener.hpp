/*!
 * \file ndflattener.hpp
 * \brief Flatten pointer-to-pointer-... arrays for MPI communication
 * \author M. Aehle
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

#include <cassert>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "../containers/C2DContainer.hpp"
#include "../parallelization/mpi_structure.hpp"

/// \addtogroup Containers
/// @{

// --- Usage
/*! \page ndflattener_usage NdFlattener usage
 *
 * # What is an NdFlattener?
 *
 * NdFlattener is a container to store a deep copy of multidimensional data in a single object that provides a
 * generic []...[] interface. The data and their structure are stored in few contiguous arrays, allowing for an
 * efficient MPI gathering operation.
 *
 * # Which problem does it solve?
 * In SU2, there is a lot of multidimensional data: data depending on K>1 indices, like a wall property
 * depending on the rank, zone index, and marker index.
 *
 * The data is often stored non-contiguously as a "tree":
 * - an array of pointers, each one of which points to
 * - (repeat the previous line another (K-2) times)
 * - an array of data.
 *
 * The individual arrays' sizes are usually unrelated to each other. Therefore in addition to this tree of
 * pointers-to-...-to-data, usually there is a tree of pointers-to-...-to-size for each lower level.
 * Instead of arrays, there might be getter functions for the data and the sizes.
 *
 * In order to use multidimensional data behind such a non-generic interface as an argument in a function call,
 * we usually make the callee aware of the interface and expose all the trees: Either as separate pointers to arrays
 * (or functions), which is tedious, or by handing an object in which they naturally live, which exposes much more
 * information.
 *
 * If the callee is an MPI communication operation, this is not possible. Here, a frequent pattern is to copy the data
 * into a single contiguous array, which is then communicated alongside with the structure (lengths of the arrays).
 *
 * The NdFlattener container is a systematic way to store multidimensional data in a single contiguous array, and
 * the structure in few additional contiguous arrays. It exposes them through a generic []...[] interface and allows
 * for an efficient MPI_Allgatherv-like operation.
 *
 * To demonstrate the usage of NdFlattener, let us store a property of type su2double that depends on
 * the rank, zone index, and marker index.
 *
 * # Form the local NdFlattener
 * by a "recursive function" like this:
 *
 *     auto f_local =
 *     make_pair( nZone, [=](unsigned long iZone){
 *       return make_pair( (Geometry of iZone)->GetnMarker() , [=](unsigned long iMarker){
 *         return YOUR_PROPERTY(iZone, iMarker);
 *       });
 *     });
 *     NdFlattener<2> nd_local(f_local);
 *
 * It might be safer to explicitly capture what you need in the lambda expressions.
 *
 * You can also construct an NdFlattener without any arguments, and use NdFlattener::initialize(f_local).
 *
 * f_local is
 *  - a pair of an (exclusive) upper bound for the first index and a function that maps the first index to
 *  - a pair of an (exclusive) upper bound for the second index and a function that maps the second index to
 *  - (iterate this if there are more layers)
 *  - the value corresponding to these indices.
 * The template parameter "2" is the number of indices. The data type (default: su2double) and the index type
 * (default: unsigned long) are further optional template parameters.
 *
 * # Form the global NdFlattener
 * by "collective communication" like this:
 *
 *     NdFlattener<3> nd_global(Nd_MPI_Environment(), &nd_local);
 *
 * nd_global's first index is the rank, then the indices of nd_local follow.
 *
 * The struct Nd_MPI_Environment contains the relevant information about the parallel environment. The default
 * values are suitable to communicate su2double data. In order to communicate data of another MPI datatype,
 * use e.g. Nd_MPI_Environment(MPI_UNSIGNED_LONG).
 *
 * You can also construct an NdFlattener without any arguments, and use NdFlattener::initialize(mpi_env, &nd_local).
 *
 * # Refreshing
 * If only the data but not the indices (i.e. the structure of sublists' lengths) have changed, you can refresh
 * them with the call
 *
 *     nd_local.refresh(f_local); // , or
 *     nd_global(Nd_MPI_Environment(), &nd_local); // respectively.
 *
 * An NdFlattener constructed or initialized from a "recursive function" or "collective communication" must be
 * refreshed in the same way.
 *
 * We also provide a function NdFlattener::initialize_or_refresh, that initializes if not initialized, and refreshes
 * otherwise.
 *
 * # Look-up
 * You can access the NdFlattener via a "[]...[]" interface that resembles multidimensional arrays
 * or pointer-to-pointer-... arrays. If you provide all of the K indices, this returns a reference
 * to the data element. If you provide less than K indices, you obtain an IndexAccumulator object
 * to which you can later pass more indices. With respect to the above example:
 *
 *     std::cout << nd_global[rank][zone][marker]; // e.g. 4.0
 *     auto nd_g_rank = nd_global[rank];
 *     nd_g_rank[zone][marker] *= 2.0;
 *     std::cout << nd_g_rank[zone][marker]; // e.g. 8.0
 *
 * You may use the IndexAccumulator repeatedly as long as the accessed NdFlattener exists.
 * As NdFlatteners are always deep copies, the above block would not change nd_local or the non-contiguous
 * data from which it was copied.
 *
 * Each index is checked to lie in the correct range dictated by the NdFlattener and the
 * previous indices. You can obtain the minimal overflowing index with `size`:
 *
 *     auto nd_g_rank = nd_global[rank];
 *     std::cout << nd_g_rank[zone].size();
 *
 * When all indices except the last one are fixed, the corresponding data is stored contiguously and a pointer
 * to this 1D array can be retrieved by
 *
 *     nd_global[rank][zone].data()
 *
 * # Constness
 * The reference returned by the last `operator[]`, and the pointer returned by `data()`, are const if
 * - the accessed NdFlattener object was const, or
 * - one of the intermediate IndexAccumulator objects was qualified as const.
 *
 * This reflects the semantic meaning of const:
 *
 *     NdFlattener<2> nd2(...);
 *
 *     template<class T> set(T& t) { t[0] = 1.0; }
 *     set( nd2[0] ); // OK
 *
 *     template<class T> set_const(T const& t) { t[0] = 1.0; }
 *     set_const( nd2[0] ); // does not compile
 *
 * # Unit tests
 * The interface described here is tested in UnitTests/Common/toolboxes/ndflattener_tests.cpp, which may also
 * serve as an addition example of how to use the NdFlattener.
 */

//  --- Implementation details---
/*! \page ndflattener_implementationdetails Implementation details of NdFlattener
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
 * allocate it and during the second iteration we copy the data.
 */

template <size_t K, typename Data_t_ = su2double, typename Index_t_ = unsigned long>
class NdFlattener;

/*! Allgatherv wrapper defaulting to a copy operation if there is only a single MPI rank.
 *
 *  Introducing this was necessary because MPICH's Allgatherv behaved unexpectedly if there
 *  is only one MPI rank (seemingly ignoring displs[0] != 0).
 */
static inline void SU2_MPI_Allgatherv_safe(const void* sendbuf, int sendcount, SU2_MPI::Datatype sendtype,
                                           void* recvbuf, const int* recvcounts, const int* displs,
                                           SU2_MPI::Datatype recvtype, SU2_MPI::Comm comm) {
  if (SU2_MPI::GetSize() == 1) {
    SU2_MPI::CopyData(sendbuf, recvbuf, sendcount, sendtype, displs[0]);
  } else {
    SU2_MPI::Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
  }
}

/*! \struct Nd_MPI_Environment
 * \brief Contains information for the collective communication of NdFlatteners.
 *
 * The default arguments in the constructor are chosen in a sensible way:
 * - To communicate su2double data, just call Nd_MPI_Environment().
 * - To communicate unsigned long data, call Nd_MPI_Environment(MPI_UNSIGNED_LONG).
 */
struct Nd_MPI_Environment {
  using MPI_Allgather_t = decltype(&(SU2_MPI::Allgather));
  using MPI_Allgatherv_t = decltype(&(SU2_MPI::Allgatherv));
  using MPI_Datatype_t = typename SU2_MPI::Datatype;
  using MPI_Communicator_t = typename SU2_MPI::Comm;

  const MPI_Datatype_t mpi_data;
  const MPI_Datatype_t mpi_index;
  MPI_Communicator_t comm;
  MPI_Allgather_t MPI_Allgather_fun;
  MPI_Allgatherv_t MPI_Allgatherv_fun;
  const int rank;
  const int size;

  Nd_MPI_Environment(MPI_Datatype_t mpi_data = MPI_DOUBLE, MPI_Datatype_t mpi_index = MPI_UNSIGNED_LONG,
                     MPI_Communicator_t comm = SU2_MPI::GetComm(),
                     MPI_Allgather_t MPI_Allgather_fun = &(SU2_MPI::Allgather),
                     MPI_Allgatherv_t MPI_Allgatherv_fun = &(SU2_MPI_Allgatherv_safe))
      : mpi_data(mpi_data),
        mpi_index(mpi_index),
        comm(comm),
        MPI_Allgather_fun(MPI_Allgather_fun),
        MPI_Allgatherv_fun(MPI_Allgatherv_fun),
        rank(SU2_MPI::GetRank()),
        size(SU2_MPI::GetSize()) {}
};

namespace helpers {
/*! \class IndexAccumulator
 * \brief Data structure holding an offset for the NdFlattener, to provide a []...[]-interface.
 * \details Derived from IndexAccumulator_Base, specifying the operator[] method:
 *  - For N==1, the structure has already read all indices but one. So after this method has read the last
 *    index, return a (possibly const) reference to the data. An additional function data() should be provided.
 *  - For N>1, more indices have to be read, return an IndexAccumulator<N-1,Nd_t_lower>. Nd_t_lower is the
 *    NdFlattener type with one index less. Nd_t_lower is const-qualified if Nd_t is const-qualified, or if the
 *    index look-up is performed in a const IndexAccumulator.
 * \tparam N - Number of missing parameters
 * \tparam Nd_t - Type of the accessed NdFlattener, should be NdFlattener<N,...>
 * \tparam Check - if true, raise error if index to operator[] is not in range
 */
/*! \class IndexAccumulator_Base
 * \brief Parent class of IndexAccumulator.
 * \details IndexAccumulator provides the operator[] method.
 */
template <size_t N_, typename Nd_t_, bool Check>
class IndexAccumulator_Base {
 public:
  static constexpr size_t N = N_;
  using Nd_t = Nd_t_;
  using Index_t = typename Nd_t::Index_t;

  /*! \brief Return exclusive upper bound for next index.
   */
  Index_t size() const { return size_; }

 protected:
  Nd_t& nd;             /*!< \brief The accessed NdFlattener. */
  const Index_t offset; /*!< \brief Index in the currently accessed layer. */
  const Index_t size_;  /*!< \brief Exclusive upper bound for the next index. */

  IndexAccumulator_Base(Nd_t& nd, Index_t offset, Index_t size) : nd(nd), offset(offset), size_(size) {}

  FORCEINLINE void CheckBound(Index_t i) const {
    assert(i < size());
    if (Check) {
      if (i >= size()) SU2_MPI::Error("NdFlattener: Index out of range.", CURRENT_FUNCTION);
    }
  }
};

template <size_t N_, typename Nd_t_, bool Check = false>
class IndexAccumulator : public IndexAccumulator_Base<N_, Nd_t_, Check> {
 public:
  using Base = IndexAccumulator_Base<N_, Nd_t_, Check>;
  static constexpr size_t N = N_;
  using Nd_t = Nd_t_;
  using Index_t = typename Nd_t::Index_t;

  IndexAccumulator(Nd_t& nd, Index_t offset, Index_t size) : Base(nd, offset, size) {}

  /*! The Base of NdFlattener<K> is NdFlattener<K-1>, but do also preserve constness.
   */
  using Nd_Base_t = su2conditional_t<std::is_const<Nd_t>::value, const typename Nd_t::Base, typename Nd_t::Base>;
  /*! Return type of operator[]. */
  using LookupType = IndexAccumulator<N - 1, Nd_Base_t>;
  /*! Return type of operator[] const. */
  using LookupType_const = IndexAccumulator<N - 1, const typename Nd_t::Base>;
  using Base::CheckBound;
  using Base::nd;
  using Base::offset;
  using Base::size;

  /*! \brief Read one more index, checking whether it is in the range dictated by the NdFlattener and
   * previous indices. Non-const version.
   * \param[in] i - Index.
   */
  LookupType operator[](Index_t i) {
    CheckBound(i);
    const Index_t new_offset = nd.GetIndices()[offset + i];
    const Index_t new_size = nd.GetIndices()[offset + i + 1] - new_offset;
    return LookupType(nd, new_offset, new_size);
  }
  /*! \brief Read one more index, const version.
   * \param[in] i - Index.
   */
  LookupType_const operator[](Index_t i) const {
    CheckBound(i);
    const Index_t new_offset = nd.GetIndices()[offset + i];
    const Index_t new_size = nd.GetIndices()[offset + i + 1] - new_offset;
    return LookupType_const(nd, new_offset, new_size);
  }
};

template <typename Nd_t_, bool Check>
class IndexAccumulator<1, Nd_t_, Check> : public IndexAccumulator_Base<1, Nd_t_, Check> {
 public:
  using Base = IndexAccumulator_Base<1, Nd_t_, Check>;
  static constexpr size_t N = 1;
  using Nd_t = Nd_t_;
  using Index_t = typename Nd_t::Index_t;

  IndexAccumulator(Nd_t& nd, Index_t offset, Index_t size) : Base(nd, offset, size) {}

  /*! Return type of operator[].
   * \details Data type of NdFlattener, but do also preserve constness.
   */
  using LookupType = su2conditional_t<std::is_const<Nd_t>::value, const typename Nd_t::Data_t, typename Nd_t::Data_t>;
  using LookupType_const = const typename Nd_t::Data_t;
  using Base::CheckBound;
  using Base::nd;
  using Base::offset;
  using Base::size;

  /*! \brief Return (possibly const) reference to the corresponding data element, checking if the index is in its range.
   * Non-const version.
   * \details The returned reference is const if and only if the type of the accessed NdFlattener is const-qualified.
   * Even though the original NdFlattener might be non-const, marking any IndexAccumulator as const during index look-up
   * will invoke the const version of operator[], which turns the NdFlattener type to const.
   * \param[in] i - Last index.
   */
  LookupType& operator[](Index_t i) {
    CheckBound(i);
    return nd.data()[offset + i];
  }
  /*! \brief Return const reference to the corresponding data element, checking if the index is in its range.
   * \param[in] i - Last index.
   */
  LookupType_const& operator[](Index_t i) const {
    CheckBound(i);
    return nd.data()[offset + i];
  }

  /*! \brief Return (possibly const) pointer to data, non-const version.
   * \details If all indices except the last one are fixed, the corresponding data
   * is stored contiguously. Return a pointer to the beginning of the
   * block. If this IndexAccumulator was generated from a non-const NdFlattener, the
   * pointer is non-const, otherwise it is const.
   */
  LookupType* data() { return &(this->operator[](0)); }
  /*! \brief Return const pointer to data.
   */
  LookupType_const* data() const { return &(this->operator[](0)); }
};

}  // namespace helpers

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
 * \tparam Data_t - Type of stored array data
 * \tparam Index - Type of index
 */

template <size_t K_, typename Data_t_, typename Index_t_>
class NdFlattener : public NdFlattener<K_ - 1, Data_t_, Index_t_> {
 public:
  static constexpr size_t K = K_;
  using Data_t = Data_t_;
  using Index_t = Index_t_;

  using Base = NdFlattener<K - 1, Data_t, Index_t>;
  using CurrentLayer = NdFlattener<K, Data_t, Index_t>;
  using LowestLayer = typename Base::LowestLayer;  // the K=1 class

 private:
  /*! \brief Number of nodes in this layer.
   *
   * For the layer K=1, nNodes will be the number of data points.
   * For a layer K>1, nNodes will be the number of sublists.
   */
  Index_t nNodes = 0;

  /*! \brief Iterator used at construction, runs from 0 to (nNodes-1). */
  Index_t iNode = 0;

  /*! \brief Indices in the lower layer's indices or data array */
  std::vector<Index_t> indices;

  /*=== Getters ===*/
 public:
  Index_t* GetIndices() { return indices.data(); }
  const Index_t* GetIndices() const { return indices.data(); }

  /*=== Outputting ===*/

  /*! \brief Write in Python-list style.
   *
   * Like this: [[1, 2], [10, 20, 30]]
   */
  friend std::ostream& operator<<(std::ostream& output, NdFlattener const& nd) {
    nd.toPythonString_fromto(output, 0, nd.size());
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
  void toPythonString_fromto(std::ostream& output, Index_t from, Index_t to) const {
    output << "[";
    for (Index_t i = from; i < to;) {
      Base::toPythonString_fromto(output, indices[i], indices[i + 1]);
      if (++i < to) output << ", ";
    }
    output << "]";
  }

 public:
  /*! \brief Basic constructor. Afterwards, initialization can be done with initialize_or_refresh.
   *
   * Called recursively when a derived class (higher K) is constructed.
   */
  NdFlattener() {}

  /*! \brief Constructor which calls initialize_or_refresh.
   *
   */
  template <class... ARGS>
  NdFlattener(ARGS const&... args) {
    initialize_or_refresh(args...);
  }

  /*! \brief Initialize or refresh the NdFlattener.
   * \details Either a 'recursive function' or 'collective communication'
   * may be used. When the NdFlattener does not hold data yet, it is
   * initialized, meaning that the data are collected and the indices arrays
   * are allocated and filled. Otherwise it is refreshed, meaning that the data are
   * recollected under the assumption that the indices arrays did not change.
   */
  template <class... ARGS>
  void initialize_or_refresh(ARGS const&... args) {
    if (initialized()) {
      refresh(args...);
    } else {
      initialize(args...);
    }
  }

  /*! \brief Initialization status of the NdFlattener.
   * \returns true if the NdFlattener has been initialized
   */
  bool initialized() { return nNodes > 0; }

 protected:
  /*! \brief Allocate the indices array after \a nNodes has been determined.
   */
  void allocate() {
    indices.resize(nNodes + 1);
    indices[0] = 0;
    Base::allocate();
  }

  /*! \brief Set \a iNode to 0 in all layers.
   */
  void reset_iNode() {
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
  template <class f_type>
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
  template <class f_type>
  void refresh(f_type f) {
    reset_iNode();
    set_f(f, true);
  }

 protected:
  /*! \brief Determine the space required for reading the 'recursive function'.
   *
   * \param f - the 'recursive function'
   */
  template <class f_type>
  void count_f(f_type f) {
    Index_t nChild = f.first;
    for (Index_t iChild = 0; iChild < nChild; iChild++) {
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
  template <class f_type>
  void set_f(f_type f, bool refresh) {
    Index_t nChild = f.first;
    for (Index_t iChild = 0; iChild < nChild; iChild++) {
      Base::set_f(f.second(iChild), refresh);
      if (!refresh) {
        indices[iNode + 1] = indices[iNode] + f.second(iChild).first;
      } else {
        if (indices[iNode + 1] != indices[iNode] + f.second(iChild).first) {
          SU2_MPI::Error("NdFlattener: Structure has changed, cannot refresh.", CURRENT_FUNCTION);
        }
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
  template <typename MPI_Environment_type>
  void initialize(MPI_Environment_type const& mpi_env, Base const& local_version) {
    // [k][r] is number of all nodes in layer (k+1), rank r in the new structure
    su2matrix<Index_t> Nodes_all(K, mpi_env.size);
    // the first index decides on the rank, so there is exactly one node per rank
    for (int r = 0; r < mpi_env.size; r++) {
      nNodes += Nodes_all[K - 1][r] = 1;
    }
    // set the lower layers' nNodes and Nodes_all[k]
    Base::count_g(mpi_env, Nodes_all, local_version);

    allocate();

    indices[0] = 0;
    for (int r = 0; r < mpi_env.size; r++) {
      indices[r + 1] = indices[r] + Nodes_all[K - 2][r];
    }
    Base::set_g(mpi_env, Nodes_all, local_version);
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
  template <typename MPI_Environment_type>
  void refresh(MPI_Environment_type const& mpi_env, Base const& local_version) {
    su2matrix<Index_t> Nodes_all_0(1, mpi_env.size);
    LowestLayer::count_g(mpi_env, Nodes_all_0, local_version);
    LowestLayer::set_g(mpi_env, Nodes_all_0, local_version);
  }

 protected:
  /*! \brief Count the distributed flatteners' numbers of nodes, and set nNodes.
   *
   * \param[in] mpi_env - MPI environment for communication
   * \param[out] Nodes_all - [k][r] is set to number of nodes in layer (k+1), rank r.
   * \param[in] local_version - local instance to be send to the other processes
   */
  void count_g(Nd_MPI_Environment const& mpi_env, su2matrix<Index_t>& Nodes_all, CurrentLayer const& local_version) {
    nNodes = 0;
    // gather numbers of nodes in the current layer from all processes
    mpi_env.MPI_Allgather_fun(&(local_version.nNodes), 1, mpi_env.mpi_index, Nodes_all[K - 1], 1, mpi_env.mpi_index,
                              mpi_env.comm);
    for (int r = 0; r < mpi_env.size; r++) {
      nNodes += Nodes_all[K - 1][r];
    }
    Base::count_g(mpi_env, Nodes_all, local_version);
  }

  /*! \brief Gather the distributed flatteners' data and index arrays into the allocated arrays.
   *
   * \param[in] mpi_env - MPI environment for communication
   * \param[in] Nodes_all - [k][r] is the number of nodes in layer (k+1), rank r.
   * \param[in] local_version - local instance to be sent to the other processes
   */
  void set_g(Nd_MPI_Environment const& mpi_env, su2matrix<Index_t> const& Nodes_all,
             CurrentLayer const& local_version) {
    std::vector<int> Nodes_all_K_as_int(mpi_env.size);
    std::vector<int> Nodes_all_k_cumulated(mpi_env.size + 1);
    //< [r] is the number of nodes in the current layer, summed over all processes with rank below r, **plus one**.
    // Used as displacements in Allgatherv, as we do not want to transfer the initial zeros, but we want to transfer the
    // last element of indices, which is the local nNodes of the layer below. Note that MPI needs indices of type 'int'.
    Nodes_all_k_cumulated[0] = 1;
    for (int r = 0; r < mpi_env.size; r++) {
      Nodes_all_k_cumulated[r + 1] = Nodes_all_k_cumulated[r] + Nodes_all[K - 1][r];
      Nodes_all_K_as_int[r] = Nodes_all[K - 1][r];
    }
    mpi_env.MPI_Allgatherv_fun(local_version.indices.data() + 1, Nodes_all[K - 1][mpi_env.rank], mpi_env.mpi_index,
                               indices.data(), Nodes_all_K_as_int.data(), Nodes_all_k_cumulated.data(),
                               mpi_env.mpi_index, mpi_env.comm);
    // shift indices
    for (int r = 1; r < mpi_env.size; r++) {
      Index_t first_entry_to_be_shifted = Nodes_all_k_cumulated[r];
      Index_t last_entry_to_be_shifted = Nodes_all_k_cumulated[r + 1] - 1;
      Index_t shift = indices[first_entry_to_be_shifted - 1];
      for (Index_t i = first_entry_to_be_shifted; i <= last_entry_to_be_shifted; i++) {
        indices[i] += shift;
      }
    }

    Base::set_g(mpi_env, Nodes_all, local_version);
  }

  /*== Data access ==*/
 public:
  Index_t size() const {  // should not be called by recursion, is incorrect in lower layers!
    return nNodes;
  }

  /*! \brief Look-up with IndexAccumulator, non-const version.
   */
  helpers::IndexAccumulator<K - 1, NdFlattener<K - 1, Data_t, Index_t> > operator[](Index_t i0) {
    return helpers::IndexAccumulator<K, NdFlattener<K, Data_t, Index_t> >(*this, 0, size())[i0];
  }
  /*! \brief Look-up with IndexAccumulator, const version.
   */
  helpers::IndexAccumulator<K - 1, const NdFlattener<K - 1, Data_t, Index_t> > operator[](Index_t i0) const {
    return helpers::IndexAccumulator<K, const NdFlattener<K, Data_t, Index_t> >(*this, 0, size())[i0];
  }
};

template <typename Data_t_, typename Index_t_>
class NdFlattener<1, Data_t_, Index_t_> {
 public:
  static constexpr size_t K = 1;
  using Data_t = Data_t_;
  using Index_t = Index_t_;

  using CurrentLayer = NdFlattener<1, Data_t, Index_t>;
  using LowestLayer = CurrentLayer;

 private:
  Index_t nNodes = 0;
  Index_t iNode = 0;
  std::vector<Data_t> data_;

  /*=== Outputting ===*/
 protected:
  void toPythonString_fromto(std::ostream& output, Index_t from, Index_t to) const {
    output << "[";
    for (Index_t i = from; i < to;) {
      output << data_[i];
      if (++i < to) output << ", ";
    }
    output << "]";
  }

 public:
  NdFlattener(void) {}

  template <class... ARGS>
  NdFlattener(ARGS const&... args) {
    initialize_or_refresh(args...);
  }

  template <class... ARGS>
  void initialize_or_refresh(ARGS const&... args) {
    if (initialized()) {
      refresh(args...);
    } else {
      initialize(args...);
    }
  }

  bool initialized() { return nNodes > 0; }

  // Functionality to initialize/refresh from a recursive function
  // could be desirable also for N=1, in order to gather from such
  // NdFlatteners an NdFlattener with N=2.
  // Gathering an NdFlattener with N=1 is not meaningful however.

  template <class f_type>
  void initialize(f_type f) {
    count_f(f);
    allocate();
    set_f(f, false);
  }

  template <class f_type>
  void refresh(f_type f) {
    reset_iNode();
    set_f(f, true);
  }

 protected:
  void allocate() { data_.resize(nNodes); }

  void reset_iNode() { iNode = 0; }

  /*=== Construct from 'recursive function' ===*/
 protected:
  template <typename f_type>
  void count_f(f_type f) {
    nNodes += f.first;
  }

  template <typename f_type>
  void set_f(f_type f, bool refresh) {
    Index_t nChild = f.first;
    for (Index_t iChild = 0; iChild < nChild; iChild++) {
      data_[iNode] = f.second(iChild);
      iNode++;
    }
  }

  /*=== Construct with Allgatherv ===*/
 protected:
  void count_g(Nd_MPI_Environment const& mpi_env, su2matrix<Index_t>& Nodes_all, CurrentLayer const& local_version) {
    nNodes = 0;
    // gather numbers of nodes in the current layer from all processes
    mpi_env.MPI_Allgather_fun(&(local_version.nNodes), 1, mpi_env.mpi_index, Nodes_all[0], 1, mpi_env.mpi_index,
                              mpi_env.comm);
    for (int r = 0; r < mpi_env.size; r++) {
      nNodes += Nodes_all[0][r];
    }
  }

  void set_g(Nd_MPI_Environment const& mpi_env, su2matrix<Index_t> const& Nodes_all,
             CurrentLayer const& local_version) {
    std::vector<int> Nodes_all_0_as_int(mpi_env.size);
    std::vector<int> Nodes_all_0_cumulated(mpi_env.size + 1);
    Nodes_all_0_cumulated[0] = 0;
    for (int r = 0; r < mpi_env.size; r++) {
      Nodes_all_0_cumulated[r + 1] = Nodes_all_0_cumulated[r] + Nodes_all[0][r];
      Nodes_all_0_as_int[r] = Nodes_all[0][r];
    }

    mpi_env.MPI_Allgatherv_fun(local_version.data_.data(), Nodes_all[0][mpi_env.rank], mpi_env.mpi_data, data_.data(),
                               Nodes_all_0_as_int.data(), Nodes_all_0_cumulated.data(), mpi_env.mpi_data, mpi_env.comm);
  }

  /*== Access to data ==*/
  // Calling the following functions is a bit strange, because if you have only one level,
  // you do not really benefit from NdFlattener's functionality.
 public:
  Index_t size() const {  // should not be called by recursion, is incorrect in lower layers!
    return nNodes;
  }
  /*! \brief Data look-up, non-const version.
   */
  Data_t& operator[](Index_t i0) { return data_[i0]; }
  /*! \brief Data look-up, const version.
   */
  const Data_t& operator[](Index_t i0) const { return data_[i0]; }

  Data_t* data() { return data_.data(); }

  const Data_t* data() const { return data_.data(); }
};

/// @}
