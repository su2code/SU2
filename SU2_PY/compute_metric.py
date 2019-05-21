#!/usr/bin/env python

## \file compute_metric.py
#  \brief Python script for calculating a metric field
#  \author B. Munguia
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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

# imports
import numpy as np
from optparse import OptionParser
import os, sys, shutil, copy, os.path
sys.path.append(os.environ['SU2_RUN'])
import SU2
from SU2.io.redirect import output as redirect_output

def main():
# Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-c", "--compute",    dest="compute",    default="True",
                      help="COMPUTE direct and adjoint problem", metavar="COMPUTE")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet all SU2 output", metavar="QUIET")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")

    (options, args)=parser.parse_args()
    partitions = int( options.partitions )
    compute    = options.compute.upper() == 'TRUE'
    quiet      = options.quiet.upper() == 'TRUE'
    nzones     = int( options.nzones )

    # load config, start state
    config = SU2.io.Config(options.filename)
    config.NUMBER_PART = partitions
    config.NZONES      = int( nzones )

    state  = SU2.io.State()

    # find solution files if they exist
    state.find_files(config)

    # run flow and adjoint if requested
    if compute:
        print "\n\n =================== Computing Flow Solution =================== \n\n"

        # make copies
        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(state)

        konfig.MATH_PROBLEM = 'DIRECT'

        # send output to a folder
        folderName = 'PRIMAL/'
        if os.path.isdir(folderName):
           os.system('rm -R '+folderName)
        os.system('mkdir ' + folderName)
        sendOutputFiles(konfig, folderName)

        # run su2
        if quiet:
            log = 'PRIMAL/log_primal.out'
        else:
            log = None
        with redirect_output(log):   
            info = SU2.run.CFD(konfig)
            ztate.update(info)
    
            # Solution merging
            konfig.SOLUTION_FLOW_FILENAME = konfig.RESTART_FLOW_FILENAME
            info = SU2.run.merge(konfig)
            ztate.update(info)

        print "\n\n =================== Computing Adjoint Solution =================== \n\n"

        # make copies
        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(state)

        konfig.SOLUTION_FLOW_FILENAME = 'PRIMAL/' + konfig.RESTART_FLOW_FILENAME
        konfig.MATH_PROBLEM = 'DISCRETE_ADJOINT'

        # send output to a folder
        folderName = 'ADJOINT/'
        if os.path.isdir(folderName):
           os.system('rm -R '+folderName)
        os.system('mkdir ' + folderName)
        sendOutputFiles(konfig, folderName)

        # run su2
        if quiet:
            log = 'ADJOINT/log_adjoint.out'
        else:
            log = None
        with redirect_output(log):   
            info = SU2.run.CFD(konfig)
            ztate.update(info)
    
            # Solution merging
            konfig.SOLUTION_ADJ_FILENAME = konfig.RESTART_ADJ_FILENAME
            info = SU2.run.merge(konfig)
            ztate.update(info)

    else:
        # send output to a folder
        folderName = 'PRIMAL/'
        if os.path.isdir(folderName):
           os.system('rm -R '+folderName)
        os.system('mkdir ' + folderName)
        os.system('cp ' + konfig.SOLUTION_FLOW_FILENAME + ' ' + folderName)

        folderName = 'ADJOINT/'
        if os.path.isdir(folderName):
           os.system('rm -R '+folderName)
        os.system('mkdir ' + folderName)
        os.system('cp ' + konfig.SOLUTION_ADJ_FILENAME + ' ' + folderName)

    # compute metric field
    print "\n\n =================== Computing Metric Field =================== \n\n"

    # make copies
    konfig = copy.deepcopy(config)
    ztate  = copy.deepcopy(state)

    konfig.ERROR_ESTIMATE = 'YES'
    konfig.MATH_PROBLEM = 'DISCRETE_ADJOINT'

    konfig.SOLUTION_FLOW_FILENAME = 'PRIMAL/' + konfig.RESTART_FLOW_FILENAME
    konfig.SOLUTION_ADJ_FILENAME = 'ADJOINT/' + konfig.RESTART_ADJ_FILENAME

    folderName = 'METRIC/'
    if os.path.isdir(folderName):
       os.system('rm -R '+folderName)
    os.system('mkdir ' + folderName)
    sendOutputFiles(konfig, folderName)

    if quiet:
        log = 'METRIC/log_metric.out'
    else:
        log = None
    with redirect_output(log):   
        info = SU2.run.MET(konfig)

def sendOutputFiles( config, folderName = ''):
    
    config.CONV_FILENAME = folderName + config.CONV_FILENAME
    config.BREAKDOWN_FILENAME = folderName + config.BREAKDOWN_FILENAME

    config.RESTART_FLOW_FILENAME = folderName + config.RESTART_FLOW_FILENAME
    config.RESTART_ADJ_FILENAME = folderName + config.RESTART_ADJ_FILENAME

    config.INTERPOLATED_RESTART_FILENAME = folderName + config.INTERPOLATED_RESTART_FILENAME
    config.INTERPOLATED_RESTART_ADJ_FILENAME = folderName + config.INTERPOLATED_RESTART_ADJ_FILENAME

    config.ECC_RESTART_FILENAME = folderName + config.ECC_RESTART_FILENAME
    config.ECC_RESTART_ADJ_FILENAME = folderName + config.ECC_RESTART_ADJ_FILENAME

    config.VOLUME_FLOW_FILENAME = folderName + config.VOLUME_FLOW_FILENAME
    config.SURFACE_FLOW_FILENAME = folderName + config.SURFACE_FLOW_FILENAME
    config.VOLUME_ADJ_FILENAME = folderName + config.VOLUME_ADJ_FILENAME
    config.SURFACE_ADJ_FILENAME = folderName + config.SURFACE_ADJ_FILENAME


if __name__ == "__main__":
    main()
