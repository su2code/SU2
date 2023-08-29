/*!
 * \file CMultiGridQueue.hpp
 * \brief Header of the multigrid queue class for the FVM solver.
 *        The subroutines and functions are in the <i>CMultiGridQueue.cpp</i> file.
 * \author F. Palacios
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

#include <vector>
#include "../containers/CFastFindAndEraseQueue.hpp"
#include "CGeometry.hpp"

using namespace std;

/*!
 * \class CMultiGridQueue
 * \brief Class for a multigrid queue system for the finite volume solver.
 * \author F. Palacios
 */
class CMultiGridQueue {
 private:
  using QueueType = CFastFindAndEraseQueue<>;
  vector<QueueType>
      QueueCV;            /*!< \brief Queue structure to choose the next control volume in the agglomeration process. */
  vector<short> Priority; /*!< \brief The priority is based on the number of pre-agglomerated neighbors. */
  vector<char> RightCV;   /*!< \brief In the lowest priority there are some CV that can not be agglomerated, this is the
                             way to identify them. */
  const unsigned long nPoint = 0; /*!< \brief Total number of points. */

  /*!
   * \brief Throw error with error message that the point is not in the priority list.
   */
  void ThrowPointNotInListError(unsigned long iPoint) const;

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of control volumes.
   */
  CMultiGridQueue(unsigned long npoint);

  /*!
   * \brief No default construction.
   */
  CMultiGridQueue() = delete;

  /*!
   * \brief Add a new CV to the list.
   * \param[in] newPoint - Index of the new point.
   * \param[in] numberNeighbors - Number of neighbors of the new point.
   */
  void AddCV(unsigned long newPoint, short numberNeighbors);

  /*!
   * \brief Remove a CV from the list.
   * \param[in] removePoint - Index of the control volume to be removed.
   */
  void RemoveCV(unsigned long removePoint);

  /*!
   * \brief Change a CV from a list to a different list.
   * \param[in] movePoint - Index of the control volume to be moved.
   * \param[in] numberNeighbors - New number of neighbors of the control volume.
   */
  void MoveCV(unsigned long movePoint, short numberNeighbors);

  /*!
   * \brief Increase the priority of the CV.
   * \param[in] incrPoint - Index of the control volume.
   */
  void IncrPriorityCV(unsigned long incrPoint);

  /*!
   * \brief Increase the priority of the CV.
   * \param[in] redPoint - Index of the control volume.
   */
  void RedPriorityCV(unsigned long redPoint);

  /*!
   * \brief Visualize the control volume queue.
   */
  void VisualizeQueue(void) const;

  /*!
   * \brief Visualize the priority list.
   */
  void VisualizePriority(void) const;

  /*!
   * \brief Find a new seed control volume.
   * \return Index of the new control volume.
   */
  inline long NextCV(void) const {
    if (!QueueCV.empty())
      return QueueCV.back().front();
    else
      return -1;
  }

  /*!
   * \brief Check if the queue is empty.
   * \return <code>TRUE</code> or <code>FALSE</code> depending if the queue is empty.
   */
  bool EmptyQueue(void) const;

  /*!
   * \brief Total number of control volume in the queue.
   * \return Total number of control points.
   */
  unsigned long TotalCV(void) const;

  /*!
   * \brief Update the queue with the new control volume (remove the CV and
   *        increase the priority of the neighbors).
   * \param[in] updatePoint - Index of the new point.
   * \param[in] fineGrid - Fine grid geometry.
   */
  void Update(unsigned long updatePoint, CGeometry* fineGrid);
};
