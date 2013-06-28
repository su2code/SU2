#!/usr/bin/env python 

## \file project.py
#  \brief package for optimization projects
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.3
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

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy, glob, time
import numpy as np
from .. import io   as su2io
from .. import eval as su2eval
from .. import util as su2util
from ..io import redirect_folder

inf = 1.0e20


# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Project(object):
    """ project = SU2.opt.Project(self,config,state=None, 
                                  designs=[],folder='.')
        
        Starts a project class to manage multiple designs
        
        Runs multiple design classes, avoiding redundancy
        Looks for closest design on restart
        Currently only based on DV_VALUE_NEW
        Exposes all methods of SU2.eval.design
        
        Attributes:
             config  - base config
             state   - base state
             files   - base files
             designs - list of designs
             folder  - project working folder
             data    - project design data
             
        Methods:
            Optimizer Interface
            The following methods take a design vector for input
            as a list (shape n) or numpy array (shape n or nx1 or 1xn).
            Values are returned as floats or lists or lists of lists.
            See SU2.eval.obj_f, etc for more detail.
            
            obj_f(dvs)     - objective function              : float
            obj_df(dvs)    - objective function derivatives  : list
            con_ceq(dvs)   - equality constraints            : list
            con_dceq(dvs)  - equality constraint derivatives : list[list]
            con_cieq(dvs)  - inequality constraints          : list
            con_dcieq(dvs) - inequality constraint gradients : list[list]
            
            Fucntional Interface
            The following methods take an objective function name for input.
            func(func_name)                  - function of specified name
            grad(func_name,method='ADJOINT') - gradient of specified name
            
    """  
    
    _design_folder = 'DESIGNS/DSN_*'
    _design_number = '%3d'
    
    
    def __init__( self, config, state=None , 
                  designs=None, folder='.'  ):
        
        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = su2io.next_folder(folder)        
        if designs is None: designs = []
        
        print 'New Project: %s' % (folder)
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = su2io.state.default_state(state)
        state.find_files(config)
        
        self.config  = config      # base config
        self.state   = state       # base state
        self.files   = state.FILES # base files
        self.designs = designs     # design list
        self.folder  = folder      # project folder
        self.data    = su2util.ordered_bunch() # project design data
        
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
            folders = glob.glob(self._design_folder)
            if len(folders)>0:
                sys.stdout.write('Warning, removing old designs...')
                sys.stdout.flush()
                time.sleep(7)
                sys.stdout.write(' now\n')
                for f in folders: shutil.rmtree(f)
            #: if existing designs
        return
    
    def _eval(self,config,func,*args):
        """ evalautes a config, checking for existing designs
        """
        
        konfig = copy.deepcopy(config) # design config
        config = self.config           # project config
        state  = self.state            # project state
        folder = self.folder           # project folder
        
        # check folder
        assert os.path.exists(folder) , 'cannot find project folder %s' % folder        
        
        # list project files to pull and link
        pull,link = state.pullnlink(config)
        
        # project folder redirection, don't overwrite files
        with redirect_folder(folder,pull,link,force=False) as push:        
        
            # find closest design
            closest,delta = self.find_design(konfig)
            # found existing design
            if delta == 0.0 and closest:
                design = closest
            # start new design
            else:
                design = self.new_design(konfig,closest)
            #: if new design
            
            # run design
            vals = design._eval(func,*args)
                        
            # recompile design data
            self.compile_data()
            
            # plot data
            self.plot_data()
            
            # save project
            su2io.save_data('project.pkl',self)
            
        #: with redirect folder
        
        # done, return output
        return vals
    
    def _unpack_dvs(self,dvs):
        dvs = copy.deepcopy(dvs)
        konfig = copy.deepcopy( self.config )
        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        konfig.unpack_dvs(dvs)
        return konfig, dvs
    
    def obj_f(self,dvs):
        func = su2eval.obj_f
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
        
    def obj_df(self,dvs):
        func = su2eval.obj_df
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_ceq(self,dvs):
        func = su2eval.con_ceq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_dceq(self,dvs):
        func = su2eval.con_dceq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_cieq(self,dvs):
        func = su2eval.con_cieq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_dcieq(self,dvs):
        func = su2eval.con_dcieq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def func(self,func_name,config):
        func = su2eval.func
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name)
    
    def grad(self,func_name,method,config):
        func = su2eval.grad
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name,method)
    
    def user(self,user_func,config,*args):
        raise NotImplementedError
        #return self._eval(config, user_func,*args) 
        
        
    def find_design(self,config):
        """ looks for an existing or closest design
        """        
                
        designs = self.designs
        
        keys_check = ['DV_VALUE_NEW']
        
        if not designs: 
            return [] , inf
        
        diffs = []
        for this_design in designs:
            this_config = this_design.config
            distance = config.dist(this_config,keys_check)
            diffs.append(distance) 
                        
        #: for each design 
        
        # pick closest design
        i_min = np.argmin(diffs)
        delta  = diffs[i_min]
        closest = designs[i_min]
        
        return closest, delta 
    
    def new_design(self,config,closest=None):
        """ starts a new design
        """
        
        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(self.state)
        if closest is None: closest = []
        
        # use closest design as seed
        if closest:
            # copy useful state info
            seed_folder = closest.folder
            seed_files  = closest.files
            for key in seed_files.keys():
                # ignore mesh
                if key == 'MESH': continue 
                # build file path
                name = seed_files[key]
                name = os.path.join(seed_folder,name)
                # update pull files
                ztate.FILES[key] = name

        # start new design (pulls files to folder)
        design = su2eval.Design(konfig,ztate)
        
        # update local state filenames
        for key in design.files:
            name = design.files[key]
            name = os.path.split(name)[-1]
            design.files[key] = name
        
        # add design to project 
        self.designs.append(design)        
        
        return design
    
    def compile_data(self,default=np.nan):
        """ data = SU2.opt.Project.compile_data(default=np.nan)
            builds a Bunch() of design data
            
            Inputs:
                default - value for missing values
                
            Outputs:
                data - state with items filled with list of
                values ordered by each design iteration.
                Contents:
                    VARIABLES
                    FUNCTIONS
                    GRADIENTS
                    HISTORY
                    
                
        """
        
        data = su2io.State()
        data.VARIABLES = []
        del data.FILES
        
        n_dv = 0
        
        # populate fields
        for i,design in enumerate(self.designs):
            for key in design.state.FUNCTIONS.keys():
                data.FUNCTIONS[key] = []
            for key in design.state.GRADIENTS.keys():
                data.GRADIENTS[key] = []
            for key in design.state.HISTORY.keys():
                data.HISTORY[key] = []
            this_ndv = len( design.state.design_vector() )
            if i == 0:
                n_dv = this_ndv
            else:
                if n_dv != this_ndv:
                    print 'warning - different dv vector length during reduce_data()'
                
        # populate data
        for design in self.designs:
            this_designvector = design.state.design_vector()
            data.VARIABLES.append( this_designvector )
            for key in data.FUNCTIONS.keys():
                if design.state.FUNCTIONS.has_key(key):
                    new_func = design.state.FUNCTIONS[key]
                else:
                    new_func = default
                data.FUNCTIONS[key].append(new_func)
            for key in data.GRADIENTS.keys():
                if design.state.GRADIENTS.has_key(key):
                    new_grad = design.state.GRADIENTS[key]
                else:
                    new_grad = [default] * len( this_designvector )
                data.GRADIENTS[key].append(new_grad)
            for key in data.HISTORY.keys():
                if design.state.HISTORY.has_key(key):
                    new_func = design.state.HISTORY[key].ITERATION[-1]
                else:
                    new_func = 0.0
                data.HISTORY[key].append(new_func)                
        
        self.data = data
            
        return self.data
    
    def plot_data(self):
        """ writes a tecplot file for plotting design data
        """
        output_format = self.config.OUTPUT_FORMAT
        functions     = self.data.FUNCTIONS
        history       = self.data.HISTORY
        
        data_plot     = su2util.ordered_bunch()
        data_plot.EVALUATION = range(1,len(self.designs)+1)
        data_plot.update(functions)
        data_plot.update(history)
        
        su2util.write_plot('history_project.plt',output_format,data_plot)
        
    def __repr__(self):
        return '<Project> with %i <Design>' % len(self.designs)
    def __str__(self):
        output = self.__repr__()
        return output    