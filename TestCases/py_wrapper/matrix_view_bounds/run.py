#!/usr/bin/env python

## \file run.py
#  \brief Bounds checks of the Python wrapper matrix views.
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import subprocess
import sys
import pysu2

def Views(driver):
  """
  One view of each class: CPyWrapperMatrixView, CPyWrapperMarkerMatrixView, CPyWrapper3DMatrixView.
  """
  iSolver = driver.GetSolverIndices()['HEAT']
  iMarker = driver.GetMarkerIndices()['x_plus']
  return {'Solution': driver.Solution(iSolver),
          'MarkerSolution': driver.MarkerSolution(iSolver, iMarker),
          'Gradient': driver.Gradient(iSolver)}

def main():
  """
  Run the heat solver, then check that the views accept the last valid index and
  stop with an error at the first index past the end, in every dimension.
  Usage: run.py <cfg> [<view> <dimension>], the second form is the child process.
  """
  cfg = sys.argv[1]

  if len(sys.argv) == 4:
    # Child process: access the view with the first index past the end of one dimension.
    driver = pysu2.CSinglezoneDriver(cfg, 1, 0)
    view = Views(driver)[sys.argv[2]]
    iDim = int(sys.argv[3])
    index = [0] * len(view.Shape())
    index[iDim] = view.Shape()[iDim]
    print(f'{sys.argv[2]}{tuple(index)} returned {view(*index)}')
    return

  # SU2_MPI::Error ends the process, so each access past the end runs in a child process,
  # which must stop with the "out of bounds" error.
  failed = []
  for name, nDim in (('Solution', 2), ('MarkerSolution', 2), ('Gradient', 3)):
    for iDim in range(nDim):
      child = subprocess.run([sys.executable, __file__, cfg, name, str(iDim)],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
      if child.returncode == 0 or 'out of bounds' not in child.stdout:
        returned = [line for line in child.stdout.splitlines() if ' returned ' in line]
        failed.append(f'{name}, index past the end of dimension {iDim}: exit code {child.returncode}, '
                      f'no out of bounds error. {" ".join(returned)}')

  driver = pysu2.CSinglezoneDriver(cfg, 1, 0)
  driver.StartSolver()

  # The last valid index of each view returns the value set through it.
  for name, view in Views(driver).items():
    index = [n - 1 for n in view.Shape()]
    value = view(*index) + 1.0
    view.Set(*index, value)
    if view(*index) != value:
      failed.append(f'{name}{tuple(index)} returned {view(*index)} instead of {value}')

  driver.Finalize()

  for message in failed:
    print(f'FAILED: {message}')
  if failed:
    sys.exit(1)

if __name__ == '__main__':
  main()
