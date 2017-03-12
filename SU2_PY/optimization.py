#!/usr/bin/env python

## \file optimization.py
#  \brief Python script for performing generic optimization.
#  \author R. Sanchez
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

import os, sys, shutil, copy
from optparse import OptionParser

sys.path.append(os.environ['SU2_RUN'])
import SU2


# -------------------------------------------------------------------
#  Main
# -------------------------------------------------------------------

def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-r", "--name", dest="projectname", default='',
                      help="try to restart from project file NAME", metavar="NAME")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="CONTINUOUS_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)",
                      metavar="GRADIENT")
    parser.add_option("-o", "--optimization", dest="optimization", default="SLSQP",
                      help="OPTIMIZATION techique (SLSQP, CG, BFGS, POWELL)", metavar="OPTIMIZATION")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")

    (options, args) = parser.parse_args()

    # process inputs
    options.partitions = int(options.partitions)
    options.quiet = options.quiet.upper() == 'TRUE'
    options.gradient = options.gradient.upper()

    sys.stdout.write('\n-------------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )   Release 5.0.0 \"Raven\"                             |\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                      |\n')
    sys.stdout.write('|   |___/\\___//___|   Optimization Script                               |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Lead Dev.: Dr. Francisco Palacios, Francisco.D.Palacios@boeing.com|\n')
    sys.stdout.write('|                Dr. Thomas D. Economon, economon@stanford.edu          |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Developers:                                                       |\n')
    sys.stdout.write('| - Prof. Juan J. Alonso\'s group at Stanford University.                |\n')
    sys.stdout.write('| - Prof. Piero Colonna\'s group at Delft University of Technology.      |\n')
    sys.stdout.write('| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. |\n')
    sys.stdout.write('| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  |\n')
    sys.stdout.write('| - Prof. Rafael Palacios\' group at Imperial College London.            |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| Copyright (C) 2012-2017 SU2, the open-source CFD code.                |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is free software; you can redistribute it and/or                  |\n')
    sys.stdout.write('| modify it under the terms of the GNU Lesser General Public            |\n')
    sys.stdout.write('| License as published by the Free Software Foundation; either          |\n')
    sys.stdout.write('| version 2.1 of the License, or (at your option) any later version.    |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| SU2 is distributed in the hope that it will be useful,                |\n')
    sys.stdout.write('| but WITHOUT ANY WARRANTY; without even the implied warranty of        |\n')
    sys.stdout.write('| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |\n')
    sys.stdout.write('| Lesser General Public License for more details.                       |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('| You should have received a copy of the GNU Lesser General Public      |\n')
    sys.stdout.write('| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')

    optimization(options.filename,
                 options.projectname,
                 options.partitions,
                 options.gradient,
                 options.optimization,
                 options.quiet)


#: main()

def optimization(filename,
                 projectname='',
                 partitions=0,
                 gradient='CONTINUOUS_ADJOINT',
                 optimization='SLSQP',
                 quiet=False):


    # Object: opt. Contains config, optimization properties
    problem = SU2.io.Problem(filename)

    # Config is contained in problem.config
    problem.config.NUMBER_PART = partitions
    if quiet: problem.config.CONSOLE = 'CONCISE'

    its = int(problem.ITERATIONS)
    accu = float(problem.ACCURACY)
    x0 = copy.deepcopy(problem.x0)
    xb = copy.deepcopy(problem.xb)

    # State
    state = SU2.io.State()
    physics = SU2.io.read_physics(problem.config, problem.OBJECTIVE_FUNCTION)
    state.find_files(physics)

    # Project
    #if os.path.exists(projectname):
        #project = SU2.io.load_data(projectname)
        #project.config = config
    #else:
    project = SU2.opt.Project(problem, state)

    ## Optimize
    #if optimization == 'SLSQP':
        #SU2.opt.SLSQP(project, x0, xb, its, accu)
    #if optimization == 'CG':
    SU2.opt.CG(project, x0, xb, its, accu)
    #if optimization == 'BFGS':
        #SU2.opt.BFGS(project, x0, xb, its, accu)
    #if optimization == 'POWELL':
        #SU2.opt.POWELL(project, x0, xb, its, accu)

    ## rename project file
    #if projectname:
        #shutil.move('project.pkl', projectname)

    #return project


#: shape_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

