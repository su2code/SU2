#!/usr/bin/env python

## \file find_cfl_number.py
#  \author Trent Lukaczyk
#
# SU2, Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
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
import numpy as np
import scipy as sp
from collections import OrderedDict
import SU2

# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

# Run the cfl experiment here (local or absolute path)
RUN_FOLDER = 'EVALUATIONS_EXPERIMENT_1/'

def main():

    # problem setup
    config_filename = 'config.cfg'
    max_iterations  = 10
    bisection_range = [ 0.001 , 20.0 ]
    
    # read config file and state
    config = SU2.io.Config(config_filename)
    state  = SU2.io.State(config=config)
    state.find_files(config)
    pull,link = state.pullnlink(config)
    
    # setup result data
    results = OrderedDict()
    results['config']  = copy.deepcopy(config)
    results['history'] = OrderedDict()
    
    # move to run folder
    with SU2.io.redirect_folder(RUN_FOLDER,pull,link):
        
        # run a bracketing algorithm
        sp.optimize.brent(
            func    = run_su2,
            args    = [config,state,results],
            brack   = ( bisection_range[0] , np.mean(bisection_range) , bisection_range[1] ),
            tol     = 0.001,
            maxiter = max_iterations,
        )
        
        # save the data
        SU2.io.save_data('experiment_results.pkl',results)    
    
    #: with run_folder
    
    return


# ----------------------------------------------------------------------
#   SU2 Run Wrapper
# ----------------------------------------------------------------------

def run_su2(x,config,state,results):
    
    # unpack inputs
    cfl    = x   
    konfig = copy.deepcopy(config)
    konfig.CFL_NUMBER = cfl
    ztate  = copy.deepcopy(state)
    pull,link = ztate.pullnlink(konfig)
    
    # work in a folder for this cfl number
    with SU2.io.redirect_folder('CFL_%.4f/',pull,link):
        
        # --------------------------------------------------------------
        #   RUN SU2
        # --------------------------------------------------------------
        try:
            info = SU2.run.direct(konfig)
        except: # SU2.DivergenceFailure:
            iterations = 900000
            info = ztate
        else:
            iterations = info.HISTORY.ITERATIONS[-1]

        # pack result
        result = copy.deepcopy(info)
        result.CONFIG     = konfig
        result.CFL_NUMBER = cfl   
        result.ITERATIONS = iterations
        
        # store result for this cfl number
        SU2.io.save_data('result.pkl',result)
        results['history'][x] = append(result)
    
    #: with cfl folder
    
    # write total experiment results
    SU2.io.save_data('experiment_results.pkl',results)
        
    return iterations
            
#: def run_su2()


# ----------------------------------------------------------------------
#   Call Main
# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()

