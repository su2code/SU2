/*!
 * \file CMultiGridQueue.hpp
 * \brief Header of the multigrid queue class for the FVM solver.
 *        The subroutines and functions are in the <i>CMultiGridQueue.cpp</i> file.
 * \author F. Palacios
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "mpi_structure.hpp"
#include "geometry/CGeometry.hpp"

using namespace std;

/*!
 * \class CMultiGridQueue
 * \brief Class for a multigrid queue system for the finite volume solver.
 * \author F. Palacios
 */
class CMultiGridQueue {
  vector<vector<unsigned long> > QueueCV; /*!< \brief Queue structure to choose the next control volume in the agglomeration process. */
  short *Priority;                        /*!< \brief The priority is based on the number of pre-agglomerated neighbors. */
  bool *RightCV;                          /*!< \brief In the lowest priority there are some CV that can not be agglomerated, this is the way to identify them */
  unsigned long nPoint;                   /*!< \brief Total number of points. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_npoint - Number of control volumes.
   */
  CMultiGridQueue(unsigned long val_npoint);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CMultiGridQueue(void);
  
  /*!
   * \brief Add a new CV to the list.
   * \param[in] val_new_point - Index of the new point.
   * \param[in] val_number_neighbors - Number of neighbors of the new point.
   */
  void AddCV(unsigned long val_new_point, unsigned short val_number_neighbors);
  
  /*!
   * \brief Remove a CV from the list.
   * \param[in] val_remove_point - Index of the control volume to be removed.
   */
  void RemoveCV(unsigned long val_remove_point);
  
  /*!
   * \brief Change a CV from a list to a different list.
   * \param[in] val_move_point - Index of the control volume to be moved.
   * \param[in] val_number_neighbors - New number of neighbors of the control volume.
   */
  void MoveCV(unsigned long val_move_point, short val_number_neighbors);
  
  /*!
   * \brief Increase the priority of the CV.
   * \param[in] val_incr_point - Index of the control volume.
   */
  void IncrPriorityCV(unsigned long val_incr_point);
  
  /*!
   * \brief Increase the priority of the CV.
   * \param[in] val_red_point - Index of the control volume.
   */
  void RedPriorityCV(unsigned long val_red_point);
  
  /*!
   * \brief Visualize the control volume queue.
   */
  void VisualizeQueue(void);
  
  /*!
   * \brief Visualize the priority list.
   */
  void VisualizePriority(void);
  
  /*!
   * \brief Find a new seed control volume.
   * \return Index of the new control volume.
   */
  long NextCV(void);
  
  /*!
   * \brief Check if the queue is empty.
   * \return <code>TRUE</code> or <code>FALSE</code> depending if the queue is empty.
   */
  bool EmptyQueue(void);
  
  /*!
   * \brief Total number of control volume in the queue.
   * \return Total number of control points.
   */
  unsigned long TotalCV(void);
  
  /*!
   * \brief Update the queue with the new control volume (remove the CV and
   increase the priority of the neighbors).
   * \param[in] val_update_point - Index of the new point.
   * \param[in] fine_grid - Fine grid geometry.
   */
  void Update(unsigned long val_update_point, CGeometry *fine_grid);
  
};
