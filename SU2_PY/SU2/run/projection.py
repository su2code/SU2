## \file projection.py
#  \brief python package for running gradient projection
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy

from .. import io   as su2io
from .. import util as su2util
from decompose import decompose as su2decomp
from interface import GPC as SU2_GPC


# ----------------------------------------------------------------------
#  Gradient Projection
# ----------------------------------------------------------------------

def projection( config, step = 1e-3 ):
    """ info = SU2.run.projection(config,step=1e-3)
        
        Runs an gradient projection with:
            SU2.run.decomp()
            SU2.run.GPC()
            
        Assumptions:
            Redundant decomposition if config.DECOMPOSED == True
            Writes tecplot file of gradients
            Adds objective suffix to gradient plot filename
            
        Inputs:
            config - an SU2 config
            step   - a float or list of floats for geometry sensitivity
                     finite difference step
            
        Outputs:
            info - SU2 State with keys:
                GRADIENTS.<config.ADJ_OBJFUNC>
                
        Updates:
            config.DECOMPOSED
            config.MATH_PROBLEM
            
        Executes in:
            ./
    """
    # local copy
    konfig = copy.deepcopy(config)
    
    # decompose
    su2decomp(konfig)
        
    # choose dv values 
    Definition_DV = konfig['DEFINITION_DV']
    n_DV          = len(Definition_DV['KIND'])
    if isinstance(step,list):
        assert len(step) == n_DV , 'unexpected step vector length'
    else:
        step = [step]*n_DV
    dv_old = [0.0]*n_DV # SU2_GPC input requirement, assumes linear superposition of design variables
    dv_new = step
    konfig.unpack_dvs(dv_new,dv_old)

    # filenames
    objective      = konfig['ADJ_OBJFUNC']    
    grad_filename  = konfig['GRAD_OBJFUNC_FILENAME']
    output_format  = konfig['OUTPUT_FORMAT']
    plot_extension = su2io.get_extension(output_format)
    adj_suffix     = su2io.get_adjointSuffix(objective)
    grad_plotname  = os.path.splitext(grad_filename)[0] + '_' + adj_suffix + plot_extension    

    # Run Projection
    SU2_GPC(konfig)
    
    # read raw gradients
    raw_gradients = su2io.read_gradients(grad_filename)
    os.remove(grad_filename)
    
    # Write Gradients
    data_plot = su2util.ordered_bunch()
    data_plot['VARIABLE']     = range(len(raw_gradients)) 
    data_plot['GRADIENT']     = raw_gradients             
    data_plot['FINDIFF_STEP'] = step                       
    su2util.write_plot(grad_plotname,output_format,data_plot)

    # gradient output dictionary
    gradients = { objective : raw_gradients }
    
    # info out
    info = su2io.State()
    info.GRADIENTS.update( gradients )
    
    return info
