#!/usr/bin/env python

## \file direct.py
#  \brief python package for running direct solutions
#  \author T. Lukaczyk, F. Palacios
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
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import copy

from .. import io  as su2io
from .merge     import merge     as su2merge
from .interface import CFD       as SU2_CFD

# ----------------------------------------------------------------------
#  Direct Simulation
# ----------------------------------------------------------------------

def direct ( config ): 
    """ info = SU2.run.direct(config)
        
        Runs an adjoint analysis with:
            SU2.run.decomp()
            SU2.run.CFD()
            SU2.run.merge()
            
        Assumptions:
            Does not rename restart filename to solution filename
            Adds 'direct' suffix to convergence filename
                        
        Outputs:
            info - SU2 State with keys:
                FUNCTIONS
                HISTORY.DIRECT
                FILES.DIRECT
                
        Updates:
            config.MATH_PROBLEM
            
        Executes in:
            ./
    """
    
    # local copy
    konfig = copy.deepcopy(config)

    # setup direct problem
    konfig['MATH_PROBLEM']  = 'DIRECT'
    konfig['CONV_FILENAME'] = konfig['CONV_FILENAME'] + '_direct'    
    
    direct_diff = konfig.get('DIRECT_DIFF','NO') == "YES"

    # Run Solution
    SU2_CFD(konfig)
    
    # multizone cases
    multizone_cases = su2io.get_multizone(konfig)

    # merge
    konfig['SOLUTION_FLOW_FILENAME'] = konfig['RESTART_FLOW_FILENAME']
    if 'FLUID_STRUCTURE_INTERACTION' in multizone_cases:
        konfig['SOLUTION_STRUCTURE_FILENAME'] = konfig['RESTART_STRUCTURE_FILENAME']
    su2merge(konfig)
    
    # filenames
    plot_format      = konfig['OUTPUT_FORMAT']
    plot_extension   = su2io.get_extension(plot_format)
    history_filename = konfig['CONV_FILENAME'] + plot_extension
    special_cases    = su2io.get_specialCases(konfig)
    
    # averaging final iterations
    final_avg = config.get('ITER_AVERAGE_OBJ',0)

    # get history and objectives
    history      = su2io.read_history( history_filename , config.NZONES)
    aerodynamics = su2io.read_aerodynamics( history_filename , config.NZONES, special_cases, final_avg )
    
    # update super config
    config.update({ 'MATH_PROBLEM' : konfig['MATH_PROBLEM']  })
                    
    # info out
    info = su2io.State()
    info.FUNCTIONS.update( aerodynamics )
    info.FILES.DIRECT = konfig['RESTART_FLOW_FILENAME']
    if 'EQUIV_AREA' in special_cases:
        info.FILES.WEIGHT_NF = 'WeightNF.dat'
    if 'INV_DESIGN_CP' in special_cases:
        info.FILES.TARGET_CP = 'TargetCp.dat'
    if 'INV_DESIGN_HEATFLUX' in special_cases:
        info.FILES.TARGET_HEATFLUX = 'TargetHeatFlux.dat'
    info.HISTORY.DIRECT = history
    
    return info
