#!/usr/bin/env python 

## \file shape_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author T. Economon, T. Lukaczyk, F. Palacios
#  \version 3.2.5 "eagle"
#
# Copyright (C) 2012-2014 SU2 <https://github.com/su2code>.
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
    parser.add_option("-p", "--oldpartitions", dest="oldpartitions", default="oldpartitions",
                      help="old number of PARTITIONS (use -n instead)", metavar="OLDPARTITIONS")
    parser.add_option("-g", "--gradient", dest="gradient", default="Adjoint",
                      help="Method for computing the GRADIENT (ADJOINT, FINDIFF, NONE)", metavar="GRADIENT")
    parser.add_option("-q", "--quiet", dest="quiet", default="False",
                      help="True/False Quiet all SU2 output (optimizer output only)", metavar="QUIET")
    parser.add_option("-c", "--cycle", dest="cycle", default=0,
                      help="number of mesh adaptation CYCLEs", metavar="CYCLE")
    parser.add_option("-i", "--its", dest="its", default=100,
                      help="number of ITERations", metavar="ITER")
    parser.add_option("-s", "--step", dest="step", default=1e-4,
                      help="finite difference STEP", metavar="STEP")
    
    (options, args)=parser.parse_args()
    
    # process inputs
    options.partitions  = int( options.partitions )
    options.cycle       = int( options.cycle )
    options.its         = int( options.its )
    options.step        = float( options.step )
    options.quiet       = options.quiet.upper() == 'TRUE'
    options.gradient    = options.gradient.upper()
    
    if options.oldpartitions != "oldpartitions":
        print ("\n IMPORTANT: -p is no longer available in SU2 v3.2.4, use -n flag instead \n")
        sys.exit()
    
    shape_optimization( options.filename    ,
                        options.projectname ,
                        options.partitions  ,
                        options.gradient    ,
                        options.quiet       ,
                        options.cycle       ,
                        options.its         ,
                        options.step         )
    
#: main()

def shape_optimization( filename                , 
                        projectname = ''        ,
                        partitions  = 0         , 
                        gradient    = 'ADJOINT' ,
                        quiet       = False     , 
                        cycle       = 0         ,
                        its         = 100       ,
                        step        = 1e-4       ):
    
    # TODO: findif step
    
    # Config
    config = SU2.io.Config(filename)
    config.NUMBER_PART = partitions
    if quiet: config.CONSOLE = 'CONCISE'
    config.GRADIENT_METHOD = gradient
    
    def_dv = config.DEFINITION_DV
    n_dv   = len(def_dv['KIND'])  
    x0     = [0.0]*n_dv # initial design
    xb     = []         # design bounds
    
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
    SU2.opt.SLSQP(project,x0,xb,its)
    
    # rename project file
    if projectname:
        shutil.move('project.pkl',projectname)
    
    return project

#: shape_optimization()


# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()

