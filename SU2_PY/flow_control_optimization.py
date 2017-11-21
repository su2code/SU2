#!/usr/bin/env python 

## \file flow_control_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author B. Munguia
#  \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
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

    parser=OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-r", "--name", dest="projectname", default='',
                      help="try to restart from project file NAME", metavar="NAME")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="DISCRETE_ADJOINT",
                      help="Method for computing the GRADIENT (CONTINUOUS_ADJOINT, DISCRETE_ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-o", "--optimization", dest="optimization", default="SLSQP",
                      help="OPTIMIZATION techique (SLSQP, CG, BFGS, POWELL)", metavar="OPTIMIZATION")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")


    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    options.nzones      = int( options.nzones )
    
    sys.stdout.write('\n-------------------------------------------------------------------------\n')
    sys.stdout.write('|    ___ _   _ ___                                                      |\n')
    sys.stdout.write('|   / __| | | |_  )   Release 5.0.0 \"Raven\"                             |\n')
    sys.stdout.write('|   \\__ \\ |_| |/ /                                                      |\n')
    sys.stdout.write('|   |___/\\___//___|   Flow Control Optimization Script                |\n')
    sys.stdout.write('|                                                                       |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Original Developers: Dr. Francisco D. Palacios.                   |\n')
    sys.stdout.write('|                          Dr. Thomas D. Economon.                      |\n')
    sys.stdout.write('-------------------------------------------------------------------------\n')
    sys.stdout.write('| SU2 Developers:                                                       |\n')
    sys.stdout.write('| - Prof. Juan J. Alonso\'s group at Stanford University.                |\n')
    sys.stdout.write('| - Prof. Piero Colonna\'s group at Delft University of Technology.      |\n')
    sys.stdout.write('| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. |\n')
    sys.stdout.write('| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  |\n')
    sys.stdout.write('| - Prof. Rafael Palacios\' group at Imperial College London.            |\n')
    sys.stdout.write('| - Prof. Edwin van der Weide\' group at the University of Twente.       |\n')
    sys.stdout.write('| - Prof. Vincent Terrapon\' group at the University of Liege.           |\n')
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

    flow_control_optimization( options.filename    ,
                        	   options.projectname ,
                        	   options.partitions  ,
                        	   options.gradient    ,
                        	   options.optimization ,
                        	   options.quiet       ,
                        	   options.nzones      )
    
#: main()

def flow_control_optimization( filename                           ,
                        	   projectname = ''                   ,
                        	   partitions  = 0                    ,
                        	   gradient    = 'DISCRETE_ADJOINT' ,
                        	   optimization = 'SLSQP'             ,
                        	   quiet       = False                ,
                        	   nzones      = 1                    ):
  
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )
    if quiet: config.CONSOLE = 'CONCISE'
    config.GRADIENT_METHOD = gradient
    
    its              = int ( config.OPT_ITERATIONS )                      # number of opt iterations
    relax_factor     = float ( config.OPT_RELAX_FACTOR )                  # line search scale
    gradient_factor  = float ( config.OPT_GRADIENT_FACTOR )               # objective function and gradient scale
    def_dv           = config.DEFINITION_DV                               # complete definition of the desing variable
    n_dv             = sum(def_dv['SIZE'])                                # number of design variables
    accu             = float ( config.OPT_ACCURACY ) * gradient_factor    # optimizer accuracy
    x0          = [0.0]*n_dv # initial design

    # Bounds for design variables
    xb_low = [-1.0E10]*n_dv
    xb_up  = [1.0E10]*n_dv

    shape_upper = [float(x) for x in config.OPT_BOUND_UPPER.split(",")]
    shape_lower = [float(x) for x in config.OPT_BOUND_LOWER.split(",")]

    afc_upper = [float(x) for x in config.OPT_BOUND_UPPER_AFC.split(",")]
    afc_lower = [float(x) for x in config.OPT_BOUND_LOWER_AFC.split(",")]

    n_kind = len(def_dv['KIND'])
    i_off = 0
    i_afc = 0
    i_ffd = 0
    for i_kind in range(n_kind):
        if def_dv['KIND'][i_kind] == 'TRANSP_DV':
        	for i in range(def_dv['SIZE'][i_kind]):
        		xb_low[i_off+i] = afc_lower[i_afc]*def_dv['SCALE'][i_kind]
        		xb_up[i_off+i]  = afc_upper[i_afc]*def_dv['SCALE'][i_kind]
            i_afc = i_afc + 1
        else:
            for i in range(def_dv['SIZE'][i_kind]):
                xb_low[i_off+i] = shape_lower[i_ffd]*def_dv['SCALE'][i_kind]
                xb_up[i_off+i]  = shape_upper[i_ffd]*def_dv['SCALE'][i_kind]
            i_ffd = i_ffd + 1
        i_off = i_off + def_dv['SIZE'][i_kind]

    xb = zip(xb_low,xb_up) # design bounds
    
    # State
    state = SU2.io.State()
    state.find_files(config)
    
    # Project
    if os.path.exists(projectname):
        project = SU2.io.load_data(projectname)
        project.config = config
    else:
        project = SU2.opt.Project(config,state)
    
    # Optimize
    if optimization == 'SLSQP':
      SU2.opt.SLSQP(project,x0,xb,its,accu)
    if optimization == 'CG':
      SU2.opt.CG(project,x0,xb,its,accu)
    if optimization == 'BFGS':
      SU2.opt.BFGS(project,x0,xb,its,accu)
    if optimization == 'POWELL':
      SU2.opt.POWELL(project,x0,xb,its,accu)


    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)
    
    return project

#: flow_control_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
