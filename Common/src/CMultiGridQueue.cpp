/*!
 * \file CMultiGridQueue.cpp
 * \brief Implementation of the multigrid queue class for the FVM solver.
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

#include "../include/CMultiGridQueue.hpp"
#include <numeric>


CMultiGridQueue::CMultiGridQueue(unsigned long npoint) :
  Priority(npoint,0),
  RightCV(npoint,true),
  nPoint(npoint) {

  /*--- Queue initialization with all the points in the finer grid ---*/
  QueueCV.emplace_back(nPoint);
  iota(QueueCV[0].begin(), QueueCV[0].end(), 0ul);

}

void CMultiGridQueue::ThrowPointNotInListError(unsigned long iPoint) const {
  char buf[200];
  SPRINTF(buf, "The CV %lu is not in the priority list.", iPoint);
  SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
}

void CMultiGridQueue::AddCV(unsigned long newPoint, short numberNeighbors) {

  const short maxNeighbors = QueueCV.size()-1;

  /*--- Basic check ---*/
  if (newPoint >= nPoint) {
    SU2_MPI::Error("The index of the CV is greater than the size of the priority list.", CURRENT_FUNCTION);
  }

  /*--- Resize the list ---*/
  if (numberNeighbors > maxNeighbors) {
    const auto newSize = numberNeighbors+1;
    if (QueueCV.capacity() < newSize)
      QueueCV.reserve(2*newSize);
    QueueCV.resize(newSize);
  }

  /*--- Find the point in the queue ---*/
  const bool inQueue = (Priority[newPoint] == numberNeighbors);

  if (!inQueue) {
    /*--- Add the control volume, and update the priority list ---*/
    QueueCV[numberNeighbors].push_back(newPoint);
    Priority[newPoint] = numberNeighbors;
  }

}

void CMultiGridQueue::RemoveCV(unsigned long removePoint) {

  /*--- Basic check ---*/
  if (removePoint >= nPoint) {
    SU2_MPI::Error("The index of the CV is greater than the size of the priority list." , CURRENT_FUNCTION);
  }

  /*--- Find priority of the Control Volume ---*/
  const auto numberNeighbors = Priority[removePoint];
  if (numberNeighbors == -1) ThrowPointNotInListError(removePoint);

  /*--- Find the point in the queue ---*/
  auto end = QueueCV[numberNeighbors].end();
  auto it = find(QueueCV[numberNeighbors].begin(), end, removePoint);
  if (it != end) QueueCV[numberNeighbors].erase(it);

  Priority[removePoint] = -1;

  /*--- Check that the size of the queue is the right one ---*/
  unsigned short sizeQueueCV = 0, iQ = 0;

  for (const auto& Q : QueueCV) {
    if (!Q.empty()) { sizeQueueCV = iQ; }
    iQ++;
  }

  /*--- Resize the queue, we need at least one element in it. ---*/
  QueueCV.resize(sizeQueueCV+1);

}

void CMultiGridQueue::MoveCV(unsigned long movePoint, short numberNeighbors) {

  if (numberNeighbors < 0) {
    numberNeighbors = 0;
    RightCV[movePoint] = false;
  }
  else {
    RightCV[movePoint] = true;
  }

  /*--- Remove the control volume ---*/
  RemoveCV(movePoint);

  /*--- Add a new control volume ---*/
  AddCV(movePoint, numberNeighbors);

}

void CMultiGridQueue::IncrPriorityCV(unsigned long incrPoint) {

  /*--- Find the priority list ---*/
  const short numberNeighbors = Priority[incrPoint];
  if (numberNeighbors == -1) ThrowPointNotInListError(incrPoint);

  /*--- Remove the control volume ---*/
  RemoveCV(incrPoint);

  /*--- Increase the priority ---*/
  AddCV(incrPoint, numberNeighbors+1);

}

void CMultiGridQueue::RedPriorityCV(unsigned long redPoint) {

  /*--- Find the priority list ---*/
  const short numberNeighbors = Priority[redPoint];
  if (numberNeighbors == -1) ThrowPointNotInListError(redPoint);

  if (numberNeighbors != 0) {
    /*--- Remove the control volume ---*/
    RemoveCV(redPoint);

    /*--- Decrease the priority ---*/
    AddCV(redPoint, numberNeighbors-1);
  }

}

void CMultiGridQueue::VisualizeQueue(void) const {

  cout << endl;
  unsigned short iQ = 0;
  for (const auto& Q : QueueCV) {
    cout << "Number of neighbors " << iQ << ": ";
    for (auto iPoint : Q) cout << iPoint << " ";
    cout << endl;
    iQ++;
  }
}

void CMultiGridQueue::VisualizePriority(void) const {

  for (auto iPoint = 0ul; iPoint < nPoint; ++iPoint) {
    cout << "Control Volume: " << iPoint << " Priority: " << Priority[iPoint] << endl;
  }
}

long CMultiGridQueue::NextCV(void) const {
  if (!QueueCV.empty()) return QueueCV.back()[0];
  else return -1;
}

bool CMultiGridQueue::EmptyQueue(void) const {

  /*--- In case there is only the no agglomerated elements (size 1), check
   *    if they can be agglomerated or if we have already finished. ---*/

  if (QueueCV.size() == 1) {
    for (auto iPoint : QueueCV[0])
      if (RightCV[iPoint]) return false;
  }
  else {
    for (size_t iQ = 1; iQ < QueueCV.size(); ++iQ)
      if (!QueueCV[iQ].empty()) return false;
  }
  return true;
}

unsigned long CMultiGridQueue::TotalCV(void) const {
  unsigned long TotalCV = 0;
  for (const auto& Q : QueueCV) TotalCV += Q.size();
  return TotalCV;
}

void CMultiGridQueue::Update(unsigned long updatePoint, CGeometry *fineGrid) {

  RemoveCV(updatePoint);

  for (auto iNode = 0u; iNode < fineGrid->node[updatePoint]->GetnPoint(); ++iNode) {
    const auto jPoint = fineGrid->node[updatePoint]->GetPoint(iNode);
    if (!fineGrid->node[jPoint]->GetAgglomerate())
      IncrPriorityCV(jPoint);
  }
}
