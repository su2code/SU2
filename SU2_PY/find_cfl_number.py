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
from scipy import optimize
import pylab as plt

import matplotlib
matplotlib.rc('font', size=18)

import SU2
from SU2.util import ordered_bunch, ordered_dict


# ----------------------------------------------------------------------
#   Inputs
# ----------------------------------------------------------------------

# Run the cfl experiment here (local or absolute path)
RUN_FOLDER      = 'EVALUATIONS_EXPERIMENT_1/'
CONFIG_FILENAME = 'config_NACA0012.cfg'
NUMBER_PART     = 4
MAX_ANALYSES    = 10
SEARCH_BRACKET  = [ 0.001, 10.0, 100.0 ]
AOA_VALUES      = [ 0.0, 10.0, 25.0 ]
LIMITER_VALUES  = [ 0.1, 1.0, 10.0 ]


# ----------------------------------------------------------------------
#   Main
# ----------------------------------------------------------------------

def main():

    # read config file and state
    print 'Read %s' % CONFIG_FILENAME
    config = SU2.io.Config(CONFIG_FILENAME)
    config.NUMBER_PART = NUMBER_PART
    
    # start state
    state  = SU2.io.State(config)
    pull,link = state.pullnlink(config)
    
    # move to run folder
    print 'Starting %s' % RUN_FOLDER
    with SU2.io.redirect_folder(RUN_FOLDER,pull,link):
        
        # start the grid experiemnt for AoA and Limiter
        grid_experiment(config,state)
        
    #: with run_folder
    
    return

#: def main()


# ----------------------------------------------------------------------
#   Grid Experiment
# ----------------------------------------------------------------------
        
def grid_experiment(config,state):
    
    # track the grid number
    i_exp = 1
    
    # run the grid
    for aoa in AOA_VALUES:
        for lim in LIMITER_VALUES:
            
            # local copy of config and state
            konfig = copy.deepcopy(config)
            ztate  = copy.deepcopy(state)
            
            # set the config
            konfig.AoA           = aoa
            konfig.LIMITER_COEFF = lim
            pull,link = state.pullnlink(konfig)
            
            # report
            print 'Grid Number %i Start' % i_exp
            print '  AoA = %.1f' % aoa
            print '  Lim = %.1f' % lim
            
            # move to grid point folder
            grid_folder = '%i_AOA=%i_LIM=%.1f' % (i_exp,aoa,lim)
            with SU2.io.redirect_folder(grid_folder,pull,link):
                
                # find the best cfl number
                find_cfl_number(konfig,state)
                
            #: with grid_folder
                
            # done
            print 'Grid Number %i Complete\n' % i_exp
            i_exp += 1
            
        #: for each lim
    #: for each aoa
    
    return

#: def grid_experiment()
    

# ----------------------------------------------------------------------
#   Find CFL Number
# ----------------------------------------------------------------------
    
def find_cfl_number(config,state):

    # setup result data
    results = ordered_bunch()
    results.CONFIG = copy.deepcopy(config)
    results.ANALYSES = ordered_dict()
        
    # run a bracketing algorithm
    sp.optimize.brent(
        func    = run_su2,
        args    = ( config, state, results ),
        brack   = tuple(SEARCH_BRACKET),
        tol     = 0.001,
        maxiter = MAX_ANALYSES-3, #analyses are used in the first bracket pass and aren't counted by this algorithm
    )
    
    ## run a search method
    #sp.optimize.fminbound(
        #func    = run_su2,
        #x1      = SEARCH_BRACKET[0],
        #x2      = SEARCH_BRACKET[2],
        #args    = ( config, state, results ),
        #xtol    = 0.001,
        #maxfun  = MAX_ANALYSES,
        #disp    = 0,
    #)
    
    
    # get final result
    cfl = np.array([ r.CFL_NUMBER for r in results.ANALYSES.values() ])
    its = np.array([ r.ITERATIONS for r in results.ANALYSES.values() ])
    ind = np.lexsort((-cfl,its))
    
    # done!
    print 'Final Result'
    print '  CFL        = %.4f' % cfl[ind[0]]
    print '  Iterations = %i'   % its[ind[0]]
    
    # save the data
    SU2.io.save_data('experiment_results.pkl',results)    
    
    return

#: def find_cfl_number()


# ----------------------------------------------------------------------
#   SU2 Run Wrapper
# ----------------------------------------------------------------------

def run_su2(x,config,state,results):
    
    # unpack inputs
    cfl    = x
    
    # check for previous result
    if results.ANALYSES.has_key(cfl):
        return results.ANALYSES[cfl].ITERATIONS
    
    # setup config
    konfig = copy.deepcopy(config)
    konfig.CFL_NUMBER = cfl
    ztate  = copy.deepcopy(state)
    pull,link = ztate.pullnlink(konfig)
    
    # some names
    i_analysis = len(results.ANALYSES.keys()) + 1
    run_folder = '%i_CFL_%.4f/'% (i_analysis,cfl)
    log_filename = 'log_DIRECT.out'
    
    # report
    print '  Analysis %i' % i_analysis
    print '    CFL = %.4f' % cfl
    
    # work in a folder for this cfl number
    with SU2.io.redirect_folder(run_folder,pull,link):
        
        # --------------------------------------------------------------
        #   Call SU2
        # --------------------------------------------------------------
        
        try:
            print '    Run SU2'
            with SU2.io.redirect_output(log_filename,sys.stderr):
                info = SU2.run.direct(konfig)
        except SU2.EvaluationFailure:
            iterations = 900000
            info = ztate
            print '    Failure'
        else:
            iterations = info.HISTORY.DIRECT.ITERATION[-1]
            print '    Success, iterations = %i' % iterations

        # pack result
        result = copy.deepcopy(info)
        result.CONFIG     = konfig
        result.CFL_NUMBER = cfl   
        result.ITERATIONS = iterations
        
        # store result for this cfl number
        SU2.io.save_data('result.pkl',result)
        results.ANALYSES[x] = result
    
    #: with cfl folder
    
    # write total experiment results
    SU2.io.save_data('experiment_results.pkl',results)
    plot_cfl_results(results)
        
    # make the third point always high, to get a propper bracket
    if i_analysis == 3:
        return 9000000    
    else:
        return iterations
            
#: def run_su2()


# ----------------------------------------------------------------------
#   Plot CFL Results
# ----------------------------------------------------------------------

def plot_cfl_results(results):
    
    cfl = np.array([ r.CFL_NUMBER for r in results.ANALYSES.values() ])
    its = np.array([ r.ITERATIONS for r in results.ANALYSES.values() ])
    evl = np.arange(len(cfl))
    
    i = np.argsort( cfl )
    cfl = cfl[i]
    its = its[i]
    evl = evl[i]
    
    its[its>1000] = np.nan
    
    plt.figure(1)
    plt.clf()
    
    plt.plot(cfl,its,'bo-',lw=2,ms=10)
    plt.xlabel('CFL Number')
    plt.ylabel('Iterations')
    for e,x,y in zip(evl,cfl,its):
        plt.text(x,y,'%i'%(e+1))
    
    plt.savefig('cfl_history.eps')
    
    plt.close()
    

# ----------------------------------------------------------------------
#   Call Main
# ----------------------------------------------------------------------

if __name__ == '__main__':
    main()

